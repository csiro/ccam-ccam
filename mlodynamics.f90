! These subroutines handle dynamics for the Mixed-Layer-Ocean model

module mlodynamics

implicit none

private
public mlodiffusion,mlorouter,mlohadv,watbdy

real, dimension(:), allocatable, save :: watbdy

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine processes horizontal diffusion
subroutine mlodiffusion

use cc_mpi
use indices_m
use map_m
use mlo
use soil_m
use vecsuv_m

implicit none

include 'newmpar.h'
include 'const_phys.h'
include 'parm.h'

integer k,i
integer, parameter :: salfilt = 0 ! Additional salinity filter (0=off, 1=Katzfey)
real, parameter :: k_smag = 2. ! 0.4 in mom3, 2. in Griffies (2000)
real hdif
real, dimension(ifull+iextra) :: uc,vc,wc,ee,dz,u,v
real, dimension(ifull+iextra) :: t_kh,xfact,yfact
real, dimension(ifull) :: dudx,dvdx,dudy,dvdy,base
real, dimension(ifull) :: cc,ff,emi,ucc,vcc,wcc
logical, dimension(ifull+iextra) :: wtr

hdif=dt*(k_smag/pi)**2
ee=0.
where(land(1:ifull))
  ee(1:ifull)=1.
end where
call bounds(ee)
wtr=ee.lt.0.5

if (.not.any(wtr(1:ifull))) return

emi=1./em(1:ifull)**2
do k=1,wlev
  call mloexpdep(1,dz(1:ifull),k,0)
  call bounds(dz)
  dz=max(dz,1.E-3)

  ! neglect bathymetry following component for now
  u=0.
  v=0.
  call mloexport(2,u(1:ifull),k,0)
  call mloexport(3,v(1:ifull),k,0)
  call boundsuv(u,v)
  dudx=0.
  dvdx=0.
  dudy=0.
  dvdy=0.
  where (wtr(1:ifull).and.wtr(ie).and.wtr(iw))
    dudx=(u(ieu)-u(iwu))*0.5*em(1:ifull)/ds
    dvdx=(v(iev)-v(iwv))*0.5*em(1:ifull)/ds
  else where (wtr(1:ifull).and.wtr(ie))
    dudx=(u(ieu)-u(1:ifull))*em(1:ifull)/ds
    dvdx=(v(iev)-u(1:ifull))*em(1:ifull)/ds
  else where (wtr(1:ifull).and.wtr(iw))
    dudx=(u(1:ifull)-u(iwu))*em(1:ifull)/ds
    dvdx=(v(1:ifull)-v(iwv))*em(1:ifull)/ds
  end where
  where (wtr(1:ifull).and.wtr(in).and.wtr(is))
    dudy=(u(inu)-u(isu))*0.5*em(1:ifull)/ds
    dvdy=(v(inv)-v(isv))*0.5*em(1:ifull)/ds
  else where (wtr(1:ifull).and.wtr(ie))
    dudy=(u(inu)-u(1:ifull))*em(1:ifull)/ds
    dvdy=(v(inv)-u(1:ifull))*em(1:ifull)/ds
  else where (wtr(1:ifull).and.wtr(iw))
    dudy=(u(1:ifull)-u(isu))*em(1:ifull)/ds
    dvdy=(v(1:ifull)-v(isv))*em(1:ifull)/ds
  end where

  uc(1:ifull) = ax(1:ifull)*u(1:ifull) + bx(1:ifull)*v(1:ifull)
  vc(1:ifull) = ay(1:ifull)*u(1:ifull) + by(1:ifull)*v(1:ifull)
  wc(1:ifull) = az(1:ifull)*u(1:ifull) + bz(1:ifull)*v(1:ifull)
  call bounds(uc)
  call bounds(vc)
  call bounds(wc)

  ! Smagorinsky
  cc=(dudx-dvdy)**2+(dudy+dvdx)**2
  t_kh(1:ifull)=sqrt(cc)*hdif*emi  ! this one with em in D terms
  call bounds(t_kh)
  xfact=0.
  yfact=0.
  where (wtr(1:ifull).and.wtr(ie))
    xfact(1:ifull) = (t_kh(ie)+t_kh)*.5
  end where
  where (wtr(1:ifull).and.wtr(in))
    yfact(1:ifull) = (t_kh(in)+t_kh)*.5
  end where
  call boundsuv(xfact,yfact)
  
  base= ( emi +                         &
          xfact(1:ifull) +              & 
          xfact(iwu) +                  &
          yfact(1:ifull) +              &
          yfact(isv) )

  ucc = ( uc(1:ifull)*emi +             &
          xfact(1:ifull)*uc(ie) +       &
          xfact(iwu)*uc(iw) +           &
          yfact(1:ifull)*uc(in) +       &
          yfact(isv)*uc(is) ) / base
  vcc = ( vc(1:ifull)*emi +             &
          xfact(1:ifull)*vc(ie) +       &
          xfact(iwu)*vc(iw) +           &
          yfact(1:ifull)*vc(in) +       &
          yfact(isv)*vc(is) ) / base
  wcc = ( wc(1:ifull)*emi +             &
          xfact(1:ifull)*wc(ie) +       &
          xfact(iwu)*wc(iw) +           &
          yfact(1:ifull)*wc(in) +       &
          yfact(isv)*wc(is) ) / base
  u(1:ifull) = ax(1:ifull)*ucc + ay(1:ifull)*vcc + az(1:ifull)*wcc
  v(1:ifull) = bx(1:ifull)*ucc + by(1:ifull)*vcc + bz(1:ifull)*wcc
  
  call mloimport(2,u(1:ifull),k,0)
  call mloimport(3,v(1:ifull),k,0)
   
  do i=0,1
    ee=0.
    call mloexport(i,ee(1:ifull),k,0)
    call bounds(ee)
    ff = ( ee(1:ifull)*emi +             &
           xfact(1:ifull)*ee(ie) +       &
           xfact(iwu)*ee(iw) +           &
           yfact(1:ifull)*ee(in) +       &
           yfact(isv)*ee(is) ) / base
    call mloimport(i,ff,k,0)
  end do
  
  ! Jack Katzfey salinity filter
  if (salfilt.eq.1) then
    ee(1:ifull)=ff
    call bounds(ee)
    xfact=0.
    yfact=0.
    where (wtr(1:ifull).and.wtr(ie))
      xfact(1:ifull)=1.
    end where
    where (wtr(1:ifull).and.wtr(in))
      yfact(1:ifull)=1.
    end where
    call boundsuv(xfact,yfact)
    ff = ( ee(1:ifull)*emi +             &
           xfact(1:ifull)*ee(ie) +       &
           xfact(iwu)*ee(iw) +           &
           yfact(1:ifull)*ee(in) +       &
           yfact(isv)*ee(is) ) / base
    call mloimport(1,ff,k,0)
  end if
  
end do

return
end subroutine mlodiffusion

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates the river routing
! This version assumes that the orography is sufficently resolved
! to determine the river flow
subroutine mlorouter

use arrays_m
use cable_ccam 
use cc_mpi
use indices_m
use map_m
use mlo
use nsibd_m
use soil_m
use soilsnow_m

implicit none

include 'newmpar.h'
include 'const_phys.h'
include 'parm.h'
include 'soilv.h'

integer i,iq
integer, dimension(ifull,4) :: xp
real, dimension(ifull) :: newwat
real, dimension(ifull,4) :: dp,slope,mslope,vel,flow
real :: xx,yy

if (.not.allocated(watbdy)) then
  allocate(watbdy(ifull+iextra))
  watbdy=0.
end if

call bounds(watbdy)

newwat=watbdy(1:ifull)
xp(:,1)=in
xp(:,2)=ie
xp(:,3)=is
xp(:,4)=iw
do i=1,4
  dp(:,i)=0.5*(ds/em(1:ifull)+ds/em(xp(:,i)))
  slope(:,i)=(zs(1:ifull)-zs(xp(:,i)))/dp(:,i)
end do

! outflow
mslope=max(slope,0.)
vel=0.35*sqrt(mslope/0.00005) ! from Hal's Mk3.5 scheme
where (mslope.gt.1.E-10)
  vel=min(max(vel,0.15),5.)
elsewhere
  vel=1.E-10
end where
do i=1,4
  flow(:,i)=watbdy(1:ifull)/(dp(:,i)/(-vel(:,i)*dt)-1.) ! (kg/m^2)
end do
newwat=newwat+sum(flow,2)
  
! inflow
mslope=max(-slope,0.)
vel=0.35*sqrt(mslope/0.00005) ! from Hal's Mk3.5 scheme
where (mslope.gt.1.E-10)
  vel=min(max(vel,0.15),5.)
elsewhere
  vel=1.E-10
end where
do i=1,4
  flow(:,i)=watbdy(xp(:,i))/(dp(:,i)/(vel(:,i)*dt)-1.) ! (kg/m^2)
  flow(:,i)=flow(:,i)*em(1:ifull)*em(1:ifull)/(em(xp(:,i))*em(xp(:,i))) ! correct for changing grid size
end do
newwat=newwat+sum(flow,2)
  
! basin 
do iq=1,ifull
  if (all(slope(iq,:).lt.-1.E-10).and.land(iq)) then
  
    ! runoff is inserted into soil moisture since CCAM has discrete land and sea points
    xx=watbdy(iq)
    if (nsib.eq.4.or.nsib.eq.6.or.nsib.eq.7) then
      call cableinflow(iq,xx)
    else
      yy=min(xx,(ssat(isoilm(iq))-wb(iq,ms))*1000.*zse(ms))
      wb(iq,ms)=wb(iq,ms)+yy/(1000.*zse(ms))
      xx=max(xx-yy,0.)
    end if
    newwat(iq)=newwat(iq)-watbdy(iq)+xx

    ! runoff is removed (assume some sub-grid scale lake is present)
    !newwat(iq)=newwat(iq)-watbdy(iq)
  end if
end do

watbdy(1:ifull)=max(newwat,0.)

return
end subroutine mlorouter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine implements some basic hydrostatic dynamics
! (note that this subroutine is a work in progress)

subroutine mlohadv

use arrays_m
use cc_mpi
use indices_m
use latlong_m
use map_m
use mlo
use soil_m
use soilsnow_m
use vecsuv_m

implicit none

include 'newmpar.h'
include 'const_phys.h'
include 'dates.h'
include 'mpif.h'
include 'parm.h'

integer iq,l,ll,ii,intsch,ierr,totits
integer jyear,jmonth,jday,jhour,jmin,mstart,mins
integer tyear,jstart
integer, dimension(ifull,wlev) :: nface
integer, dimension(12) :: ndoy
real alpha,maxloclseta,maxglobseta
real delpos,delneg,alph_pm,alph_p
real, dimension(ifull+iextra) :: neta,dd,snu,snv,ntide,pice
real, dimension(ifull) :: rhobaru,rhobarv
real, dimension(ifull) :: div,seta,w_e,rhobarsav
real, dimension(ifull) :: tnu,tsu,tev,twv
real, dimension(ifull) :: rhou,rhov
real, dimension(ifull) :: sou,spu,squ,sru
real, dimension(ifull) :: sov,spv,sqv,srv
real, dimension(ifull) :: dpsdxu,dpsdyu,dpsdxv,dpsdyv
real, dimension(ifull) :: dttdxu,dttdyu,dttdxv,dttdyv
real, dimension(ifull) :: detadxu,detadyu,detadxv,detadyv
real, dimension(ifull) :: drhobardxu,drhobardyu,drhobardxv,drhobardyv
real, dimension(ifull) :: au,bu,cu,av,bv,cv
real, dimension(ifull+iextra,wlev) :: nu,nv,nt,ns
real, dimension(ifull+iextra,wlev) :: cou,cov,cow
real, dimension(ifull+iextra,wlev) :: dep,rhobar,rho,dz
real, dimension(ifull,wlev) :: w_u,w_v,w_t,w_s,dzbar,dum
real, dimension(ifull,wlev) :: nuh,nvh,xg,yg,uau,uav
real, dimension(ifull,wlev) :: kku,llu,mmu,nnu
real, dimension(ifull,wlev) :: kkv,llv,mmv,nnv
real, dimension(:,:), allocatable, save :: oldu1,oldu2,oldv1,oldv2
real*8, dimension(ifull,wlev) :: x3d,y3d,z3d
logical, dimension(ifull+iextra) :: wtr
logical lleap
integer, parameter :: lmax=2       ! 1=predictor-only, 2+=predictor-corrector
integer, parameter :: llmax=5000   ! iterations for surface height
real, parameter :: tol = 5.E-4     ! Tolerance for GS solver
real, parameter :: sal = 0.948     ! SAL parameter
real, parameter :: rhosn=330.      ! density snow
real, parameter :: rhoic=900.      ! density ice

data ndoy/0,31,59,90,120,151,181,212,243,273,304,334/

! new z levels for including free surface eta (effectively sigma-height levels)
! newz=-eta+oldz*(1+eta/maxdepth)
! where 0<=oldz<=maxdepth and -eta<=newz<=maxdepth

intsch=mod(ktau,2)
alpha = 0.9 ! Initial SOR weight

!Define land/sea mask
dd=0.
where(land(1:ifull))
  dd(1:ifull)=1.
end where
call bounds(dd,nrows=2)
wtr=dd.lt.0.5

cou=0.
cov=0.
cow=0.
rhobar=1030.
dep=0.
dz=0.
w_t=293.
w_s=0.
w_u=0.
w_v=0.
w_e=0.
nt=293.
ns=30.
nu=0.
nv=0.
neta=0.

! ADVECT WATER ------------------------------------------------------
do ii=1,wlev
  call mloexpdep(0,dep(1:ifull,ii),ii,0)
  call mloexpdep(1,dz(:,ii),ii,0)
  call mloexport(0,w_t(:,ii),ii,0)
  call mlofill(w_t(:,ii))
  call mloexport(1,w_s(:,ii),ii,0)
  call mlofill(w_s(:,ii))
  call mloexport(2,w_u(:,ii),ii,0)
  call mloexport(3,w_v(:,ii),ii,0)
end do
call mloexport(4,w_e,0,0)
call bounds(dep,nrows=2)
call bounds(dz)
ps(1:ifull)=1.e5*exp(psl(1:ifull))

dep=max(dep,1.E-3)
dz=max(dz,1.E-3/real(wlev))

! Calculate depth arrays
dzbar(:,1)=dz(1:ifull,1)
do ii=2,wlev
  dzbar(:,ii)=dzbar(:,ii-1)+dz(1:ifull,ii)
end do
dd(1:ifull)=dzbar(:,wlev)
call bounds(dd)

! ------ NEED TO ADD PRESSURE FROM SEA ICE ------
pice(1:ifull) = ps(1:ifull) + grav*fracice*(sicedep*rhoic+0.001*snowd*rhosn)
call bounds(pice,nrows=2)

! estimate tides (need to fix base date !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! )
jyear=kdate/10000
jstart=0
if (jyear.gt.1900) then
  do tyear=1900,jyear-1
    call mloleap(tyear,lleap)
    if (lleap) jstart=jstart+1
    jstart=jstart+365
  end do
else
  do tyear=1900,jyear,-1
    call mloleap(tyear,lleap)
    if (lleap) jstart=jstart-1
    jstart=jstart-1440
  end do
end if
mstart=mstart+720 ! base time is 12Z 31 Dec 1899
jmonth=(kdate-jyear*10000)/100
jday=kdate-jyear*10000-jmonth*100
jhour=ktime/100
jmin=ktime-jhour*100
mstart=mstart+1440*(ndoy(jmonth)+jday-1) + 60*jhour + jmin ! mins from start of year
if (jmonth.gt.2) then
  call mloleap(jyear,lleap)
  if (lleap) mstart=mstart+1440
end if
mins = mtimer + mstart
call mlotide(ntide(1:ifull),rlongg,rlatt,mins,jstart)
call bounds(ntide,nrows=2)

! should replace this with pre-computed weights !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
where (wtr(in).and.wtr(ine))
  tnu=0.5*(ntide(in)+ntide(ine)) 
else where (wtr(in))
  tnu=ntide(in)
else where (wtr(ine))
  tnu=ntide(ine)
else where
  tnu=0.5*(ntide(1:ifull)+ntide(ie))
end where
where (wtr(is).and.wtr(ise))
  tsu=0.5*(ntide(is)+ntide(ise)) 
else where (wtr(is))
  tsu=ntide(is)
else where (wtr(ise))
  tsu=ntide(ise)
else where
  tsu=0.5*(ntide(1:ifull)+ntide(ie))
end where
where (wtr(ie).and.wtr(ien))
  tev=0.5*(ntide(ie)+ntide(ien)) 
else where (wtr(ie))
  tev=ntide(ie)
else where (wtr(ien))
  tev=ntide(ien)
else where
  tev=0.5*(ntide(1:ifull)+ntide(in))
end where
where (wtr(iw).and.wtr(iwn))
  twv=0.5*(ntide(iw)+ntide(iwn)) 
else where (wtr(iw))
  twv=ntide(iw)
else where (wtr(iwn))
  twv=ntide(iwn)
else where
  twv=0.5*(ntide(1:ifull)+ntide(in))
end where
dttdxu=(ntide(ie)-ntide(1:ifull))*emu(1:ifull)/ds
dttdyu=0.5*(tnu-tsu)*emu(1:ifull)/ds
dttdxv=0.5*(tev-twv)*emv(1:ifull)/ds
dttdyv=(ntide(in)-ntide(1:ifull))*emv(1:ifull)/ds

! surface pressure gradients
where (wtr(in).and.wtr(ine))
  tnu=0.5*(pice(in)+pice(ine)) 
else where (wtr(in))
  tnu=pice(in)
else where (wtr(ine))
  tnu=pice(ine)
else where
  tnu=0.5*(pice(1:ifull)+pice(ie))
end where
where (wtr(is).and.wtr(ise))
  tsu=0.5*(pice(is)+pice(ise)) 
else where (wtr(is))
  tsu=pice(is)
else where (wtr(ise))
  tsu=pice(ise)
else where
  tsu=0.5*(pice(1:ifull)+pice(ie))
end where
where (wtr(ie).and.wtr(ien))
  tev=0.5*(pice(ie)+pice(ien)) 
else where (wtr(ie))
  tev=pice(ie)
else where (wtr(ien))
  tev=pice(ien)
else where
  tev=0.5*(pice(1:ifull)+pice(in))
end where
where (wtr(iw).and.wtr(iwn))
  twv=0.5*(pice(iw)+pice(iwn)) 
else where (wtr(iw))
  twv=pice(iw)
else where (wtr(iwn))
  twv=pice(iwn)
else where
  twv=0.5*(pice(1:ifull)+pice(in))
end where
dpsdxu=(pice(ie)-pice(1:ifull))*emu(1:ifull)/ds
dpsdyu=0.5*(tnu-tsu)*emu(1:ifull)/ds
dpsdxv=0.5*(tev-twv)*emv(1:ifull)/ds
dpsdyv=(pice(in)-pice(1:ifull))*emv(1:ifull)/ds

! save arrays used to estimate currents at t+1/2
if (.not.allocated(oldu1)) then
  allocate(oldu1(ifull,wlev),oldv1(ifull,wlev))
  allocate(oldu2(ifull,wlev),oldv2(ifull,wlev))
  oldu1=w_u
  oldv1=w_v
  oldu2=w_u
  oldv2=w_v
end if

! initialise t+1 variables for predictor-corrector loop
neta(1:ifull)=w_e
nt(1:ifull,:)=w_t
ns(1:ifull,:)=w_s

totits=0
do l=1,lmax ! predictor-corrector loop

  ! calculate rhobar (unstaggered)
  ! note that rho is estimated using the first guess for neta (not implicit)
  do ii=1,wlev
    dum(:,ii)=dep(1:ifull,ii)*max(1.+neta(1:ifull)/dd(1:ifull),0.01) ! depth relative to water surface, =newz+eta and ranges between 0 and dd+neta
  end do
  call mloexpdensity(rho(1:ifull,:),nt(1:ifull,:),ns(1:ifull,:),dum,pice(1:ifull),0)
  call bounds(rho)
  rhobar(1:ifull,1)=rho(1:ifull,1)*dz(1:ifull,1)
  do ii=2,wlev
    rhobar(1:ifull,ii)=rhobar(1:ifull,ii-1)+rho(1:ifull,ii)*dz(1:ifull,ii)
  end do
  rhobar(1:ifull,:)=rhobar(1:ifull,:)/dzbar
  call bounds(rhobar,nrows=2)
  if (l.eq.1) rhobarsav=rhobar(1:ifull,wlev)

  ! estimate currents at t+1/2 for semi-Lagrangian advection
  if (l.eq.1) then
    nuh=(15.*w_u-10.*oldu1+3.*oldu2)/8. ! U at t+1/2
    nvh=(15.*w_v-10.*oldv1+3.*oldv2)/8. ! V at t+1/2
  else
    nuh=0.5*nu(1:ifull,:)+0.5*w_u ! U at t+1/2
    nvh=0.5*nv(1:ifull,:)+0.5*w_v ! V at t+1/2
  end if

  ! Calculate depature points
  call mlodeps(nuh,nvh,nface,xg,yg,x3d,y3d,z3d)

  ! Vertical advection
  ! call mlovadv

  ! Convert (u,v) to cartesian coordinates (U,V,W)
  do ii=1,wlev
    cou(1:ifull,ii)=ax(1:ifull)*w_u(:,ii)+bx(1:ifull)*w_v(:,ii)
    cov(1:ifull,ii)=ay(1:ifull)*w_u(:,ii)+by(1:ifull)*w_v(:,ii)
    cow(1:ifull,ii)=az(1:ifull)*w_u(:,ii)+bz(1:ifull)*w_v(:,ii)
  end do

  ! Horizontal advection for U,V,W (need to vertically interpolate...)
  call mloints(cou,intsch,nface,xg,yg,2)
  call mloints(cov,intsch,nface,xg,yg,2)
  call mloints(cow,intsch,nface,xg,yg,2)

  ! Rotate vector to arrival point
  call mlorot(cou(1:ifull,:),cov(1:ifull,:),cow(1:ifull,:),x3d,y3d,z3d)

  ! Convert (U,V,W) back to conformal cubic coordinates
  do ii=1,wlev
    where (wtr(1:ifull))
      nu(1:ifull,ii)=ax(1:ifull)*cou(1:ifull,ii)+ay(1:ifull)*cov(1:ifull,ii)+az(1:ifull)*cow(1:ifull,ii)
      nv(1:ifull,ii)=bx(1:ifull)*cou(1:ifull,ii)+by(1:ifull)*cov(1:ifull,ii)+bz(1:ifull)*cow(1:ifull,ii)
    elsewhere
      nu(1:ifull,ii)=0.
      nv(1:ifull,ii)=0.
    end where
  end do
  
  ! Horizontal advection for T,S
  nt(1:ifull,:)=w_t
  ns(1:ifull,:)=w_s
  call mloints(nt,intsch,nface,xg,yg,2)
  call mloints(ns,intsch,nface,xg,yg,5)
  ns=max(ns,0.)

  ! stagger the following values for t=tstar
  ! This is the same as eps=1. in JLM's atmosphere semi-Lagrangian dynamics
  do ii=1,wlev
    uau(:,ii)=nu(1:ifull,ii)+dt*f(1:ifull)*nv(1:ifull,ii)
    uav(:,ii)=nv(1:ifull,ii)-dt*f(1:ifull)*nu(1:ifull,ii)
  end do
  call mlostaguv(uau,uav,cou(1:ifull,:),cov(1:ifull,:))
  
  sou=0.
  spu=0.
  squ=0.
  sru=0.
  sov=0.
  spv=0.
  sqv=0.
  srv=0.
  
  ! Precompute U,V current and integral terms
  ! staggered terms for currents
  do ii=1,wlev
    rhou=0.5*(rho(1:ifull,ii)+rho(ie,ii))
    rhov=0.5*(rho(1:ifull,ii)+rho(in,ii))
    rhobaru=(rhobar(1:ifull,ii)+rhobar(ie,ii))*0.5
    rhobarv=(rhobar(1:ifull,ii)+rhobar(in,ii))*0.5

    ! u^(t+1) = au^(t*) + bu*dpdx + cu*dpdy (staggered)
    ! v^(t+1) = av^(t*) + bv*dpdy + cv*dpdx (staggered)

    au=cou(1:ifull,ii)/(1.+dt*dt*fu(1:ifull)*fu(1:ifull))
    bu=-dt/(rhou*(1.+dt*dt*fu(1:ifull)*fu(1:ifull)))
    cu=dt*fu(1:ifull)*bu

    av=cov(1:ifull,ii)/(1.+dt*dt*fv(1:ifull)*fv(1:ifull))
    bv=-dt/(rhov*(1.+dt*dt*fv(1:ifull)*fv(1:ifull)))
    cv=-dt*fv(1:ifull)*bv

    where (wtr(in).and.wtr(ine))
      tnu=0.5*(rhobar(in,ii)+rhobar(ine,ii)) 
    else where (wtr(in))
      tnu=rhobar(in,ii)
    else where (wtr(ine))
      tnu=rhobar(ine,ii)
    else where
      tnu=0.5*(rhobar(1:ifull,ii)+rhobar(ie,ii))
    end where
    where (wtr(is).and.wtr(ise))
      tsu=0.5*(rhobar(is,ii)+rhobar(ise,ii)) 
    else where (wtr(is))
      tsu=rhobar(is,ii)
    else where (wtr(ise))
      tsu=rhobar(ise,ii)
    else where
      tsu=0.5*(rhobar(1:ifull,ii)+rhobar(ie,ii))
    end where
    where (wtr(ie).and.wtr(ien))
      tev=0.5*(rhobar(ie,ii)+rhobar(ien,ii)) 
    else where (wtr(ie))
      tev=rhobar(ie,ii)
    else where (wtr(ien))
      tev=rhobar(ien,ii)
    else where
      tev=0.5*(rhobar(1:ifull,ii)+rhobar(in,ii))
    end where
    where (wtr(iw).and.wtr(iwn))
      twv=0.5*(rhobar(iw,ii)+rhobar(iwn,ii)) 
    else where (wtr(iw))
      twv=rhobar(iw,ii)
    else where (wtr(iwn))
      twv=rhobar(iwn,ii)
    else where
      twv=0.5*(rhobar(1:ifull,ii)+rhobar(in,ii))
    end where
 
    drhobardxu=(rhobar(ie,ii)-rhobar(1:ifull,ii))*emu(1:ifull)/ds !-dzdx*drhobardz
    drhobardyu=0.5*(tnu-tsu)*emu(1:ifull)/ds                      !-dzdy*drhobardz
    drhobardxv=0.5*(tev-twv)*emv(1:ifull)/ds                      !-dzdx*drhobardz
    drhobardyv=(rhobar(in,ii)-rhobar(1:ifull,ii))*emv(1:ifull)/ds !-dzdy*drhobardz

    !nu=au+bu*dppdxu+cu*dppdyu
    !nv=av+bv*dppdyv+cv*dppdxv

    !dppdxu=dpsdxu+grav*depu*(1+etau/ddu)*drhobardxu+(grav*rhobaru+(1-sal)*grav)*detadxu+grav*dttdxu
    !dppdyu=dpsdyu+grav*depu*(1+etau/ddu)*drhobardyu+(grav*rhobaru+(1-sal)*grav)*detadyu+grav*dttdyu
    !dppdxv=dpsdxv+grav*depv*(1+etav/ddv)*drhobardxv+(grav*rhobarv+(1-sal)*grav)*detadxv+grav*dttdxv
    !dppdyv=dpsdyv+grav*depv*(1+etav/ddv)*drhobardyv+(grav*rhobarv+(1-sal)*grav)*detadyv+grav*dttdyv
    
    !nu=kku+llu*(1+etau/ddu)+mmu*detadxu+nnu*detadyu
    !nv=kkv+llv*(1+etav/ddv)+mmv*detadyv+nnv*detadxv
    
    !int nu dz = sou+spu*(1+etau/ddu)+squ*detadxu+sru*detadyu
    !int nv dz = sov+spv*(1+etav/ddv)+sqv*detadyv+srv*detadxv

    kku(:,ii)=au+bu*(dpsdxu+grav*rhou*dttdxu)+cu*(dpsdyu+grav*rhou*dttdyu)
    llu(:,ii)=(bu*grav*drhobardxu+cu*grav*drhobardyu)*0.5*(dep(1:ifull,ii)+dep(ie,ii))
    mmu(:,ii)=bu*grav*(rhobaru+rhou*(1.-sal))
    nnu(:,ii)=cu*grav*(rhobaru+rhou*(1.-sal))

    kkv(:,ii)=av+bv*(dpsdyv+grav*rhov*dttdyv)+cv*(dpsdxv+grav*rhov*dttdxv)
    llv(:,ii)=(bv*grav*drhobardyv+cv*grav*drhobardxv)*0.5*(dep(1:ifull,ii)+dep(in,ii))
    mmv(:,ii)=bv*grav*(rhobarv+rhov*(1.-sal))
    nnv(:,ii)=cv*grav*(rhobarv+rhov*(1.-sal))
    
    sou=sou+kku(:,ii)*0.5*(dz(1:ifull,ii)+dz(ie,ii))
    spu=spu+llu(:,ii)*0.5*(dz(1:ifull,ii)+dz(ie,ii))
    squ=squ+mmu(:,ii)*0.5*(dz(1:ifull,ii)+dz(ie,ii))
    sru=sru+nnu(:,ii)*0.5*(dz(1:ifull,ii)+dz(ie,ii))
    
    sov=sov+kkv(:,ii)*0.5*(dz(1:ifull,ii)+dz(in,ii))
    spv=spv+llv(:,ii)*0.5*(dz(1:ifull,ii)+dz(in,ii))
    sqv=sqv+mmv(:,ii)*0.5*(dz(1:ifull,ii)+dz(in,ii))
    srv=srv+nnv(:,ii)*0.5*(dz(1:ifull,ii)+dz(in,ii))
    
  end do

    ! update land boundaries
  where (.not.wtr(1:ifull).or..not.wtr(ie))
    sou=0.
    spu=0.
    squ=0.
    sru=0.
  end where
  where (.not.wtr(1:ifull).or..not.wtr(in))
    sov=0.
    spv=0.
    sqv=0.
    srv=0.
  end where

  ! Now iteratively solve for eta
  do ll=1,llmax

    ! calculate neta gradients
    call bounds(neta,nrows=2)
    where (wtr(in).and.wtr(ine))
      tnu=0.5*(neta(in)+neta(ine)) 
    else where (wtr(in))
      tnu=neta(in)
    else where (wtr(ine))
      tnu=neta(ine)
    else where
      tnu=0.5*(neta(1:ifull)+neta(ie))
    end where
    where (wtr(is).and.wtr(ise))
      tsu=0.5*(neta(is)+neta(ise)) 
    else where (wtr(is))
      tsu=neta(is)
    else where (wtr(ise))
      tsu=neta(ise)
    else where
      tsu=0.5*(neta(1:ifull)+neta(ie))
    end where
    where (wtr(ie).and.wtr(ien))
      tev=0.5*(neta(ie)+neta(ien)) 
    else where (wtr(ie))
      tev=neta(ie)
    else where (wtr(ien))
      tev=neta(ien)
    else where
      tev=0.5*(neta(1:ifull)+neta(in))
    end where
    where (wtr(iw).and.wtr(iwn))
      twv=0.5*(neta(iw)+neta(iwn)) 
    else where (wtr(iw))
      twv=neta(iw)
    else where (wtr(iwn))
      twv=neta(iwn)
    else where
      twv=0.5*(neta(1:ifull)+neta(in))
    end where
 
    detadxu=(neta(ie)-neta(1:ifull))*emu(1:ifull)/ds
    detadyu=0.5*(tnu-tsu)*emu(1:ifull)/ds
    detadxv=0.5*(tev-twv)*emv(1:ifull)/ds
    detadyv=(neta(in)-neta(1:ifull))*emv(1:ifull)/ds

    snu(1:ifull)=sou+spu*(1.+(neta(ie)+neta(1:ifull))/(dd(ie)+dd(1:ifull)))+squ*detadxu+sru*detadyu
    snv(1:ifull)=sov+spv*(1.+(neta(in)+neta(1:ifull))/(dd(in)+dd(1:ifull)))+sqv*detadyv+srv*detadxv
    snu(1:ifull)=snu(1:ifull)*(1.+(neta(ie)+neta(1:ifull))/(dd(ie)+dd(1:ifull)))
    snv(1:ifull)=snv(1:ifull)*(1.+(neta(in)+neta(1:ifull))/(dd(in)+dd(1:ifull)))

    call boundsuv(snu,snv)
    div=(snu(1:ifull)-snu(iwu)+snv(1:ifull)-snv(isv))*em(1:ifull)/ds !-dzdx*dsnudz-dzdy*dsnvdz

    ! Update neta
    seta=-neta(1:ifull)+w_e-dt*div

    if (alpha*maxval(abs(seta)).gt.1.) then ! attempt to bring model into balance
      alpha=min(0.1*alpha,1./maxval(abs(seta)))
      alpha=max(alpha,1.E-4)
    end if
    seta=max(seta,(5.-dd(1:ifull)-neta(1:ifull))/alpha) ! this should become a land point
    neta(1:ifull)=alpha*seta+neta(1:ifull)
    where (.not.wtr(1:ifull))
      neta(1:ifull)=0.
      seta=0.
    end where
    
    maxloclseta=maxval(abs(seta))
    call MPI_AllReduce(maxloclseta,maxglobseta,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
    if (maxglobseta.lt.tol.and.ll.gt.2) exit    
    !print *,"err  ",l,ll,maxglobseta,alpha
    !print *,"xnmag ",l,ll,maxval(neta(1:ifull)),minval(neta(1:ifull))
    !print *,"xndiv ",l,ll,maxval(div),minval(div)
    !print *,"xnloc ",l,ll,maxloc(neta(1:ifull)),minloc(neta(1:ifull))
    !print *,"netae",l,ll,neta(iw(6865)),neta(6865),neta(ie(6865))
    !print *,"netan",l,ll,neta(is(6865)),neta(6865),neta(in(6865))
    !print *,"dde  ",l,ll,dd(iw(6865)),dd(6865),dd(ie(6865))
    !print *,"ddn  ",l,ll,dd(is(6865)),dd(6865),dd(in(6865))
    !print *,"div,w_e ",l,ll,div(6865),w_e(6865),seta(6865)

    totits=totits+1
  end do


  ! Update currents once neta is calculated
  do ii=1,wlev
  
    ! update currents (staggered)
    nu(1:ifull,ii)=kku(:,ii)+llu(:,ii)*(1.+(neta(ie)+neta(1:ifull))/(dd(ie)+dd(1:ifull)))+mmu(:,ii)*detadxu+nnu(:,ii)*detadyu
    nv(1:ifull,ii)=kkv(:,ii)+llv(:,ii)*(1.+(neta(in)+neta(1:ifull))/(dd(in)+dd(1:ifull)))+mmv(:,ii)*detadyv+nnv(:,ii)*detadxv
    
    ! update land boundaries
    where (.not.wtr(1:ifull).or..not.wtr(ie))
      nu(1:ifull,ii)=0.
    end where
    where (.not.wtr(1:ifull).or..not.wtr(in))
      nv(1:ifull,ii)=0.
    end where
  end do

  ! unstagger nu and nv
  call mlounstaguv(nu(1:ifull,:),nv(1:ifull,:),nuh,nvh)
  nu(1:ifull,:)=nuh
  nv(1:ifull,:)=nvh
  
end do

! mass conservation
neta(1:ifull)=neta(1:ifull)*rhobar(1:ifull,wlev)
w_e=w_e*rhobarsav
dum(:,1)=neta(1:ifull)-w_e
call ccglobal_posneg(dum(:,1),delpos,delneg)
alph_p = sqrt( -delneg/max(1.e-20,delpos))
alph_pm=1./max(1.e-20,alph_p)
neta(1:ifull) = w_e + alph_p*max(0.,dum(:,1)) + alph_pm*min(0.,dum(:,1))
neta(1:ifull)=neta(1:ifull)/rhobar(1:ifull,wlev)
w_e=w_e/rhobarsav

call MPI_AllReduce(maxloclseta,maxglobseta,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
if (myid==0.and.(ktau.le.100.or.maxglobseta.gt.tol)) then
  write(6,*) "MLODYNAMICS ",totits,maxglobseta
end if

oldu2=oldu1
oldv2=oldv1
oldu1=w_u
oldv1=w_v

do ii=1,wlev
  call mloimport(0,nt(1:ifull,ii),ii,0)
  call mloimport(1,ns(1:ifull,ii),ii,0)
  call mloimport(2,nu(1:ifull,ii),ii,0)
  call mloimport(3,nv(1:ifull,ii),ii,0)
end do
call mloimport(4,neta(1:ifull),0,0)

! ADVECT ICE --------------------------------------------------------
!do ii=1,wlev
!  call mloexpice(w_u(:,1),8,0)
!  call mloexpice(w_v(:,1),9,0)
!end do

!if (.not.allocated(oldiceu1)) then
!  allocate(oldiceu1(ifull))
!  allocate(oldicev1(ifull))
!  oldiceu1=w_u(:,1)
!  oldicev1=w_v(:,1)
!end if

!do l=1,lmax
!
!  if (l.eq.1) then
!    nuh(1:ifull,1)=1.5*w_u(:,1)-0.5*oldiceu1 ! U at t+1/2
!    nvh(1:ifull,1)=1.5*w_v(:,1)-0.5*oldicev1 ! V at t+1/2
!  else
!    nuh(1:ifull,1)=0.5*nu(1:ifull,1)+0.5*w_u(:,1) ! U at t+1/2
!    nvh(1:ifull,1)=0.5*nv(1:ifull,1)+0.5*w_v(:,1) ! V at t+1/2
!  end if
!
!  ! Calculate depature points
!  call mlodeps(nuh(:,1),nvh(:,1),nface(:,1),xg(:,1),yg(:,1))
!
!! Account for ice deformation
!!...
!
!! Horizontal advector for ice fraction
!  nt(:,1)=i_fracice
!  call bounds(nt(:,1))
!  call mloints(nt(:,1),intsch,nface(:,1),xg(:,1),yg(:,1),2)
!  nt=max(nt,0.)
!  newfracice=nt(:,1)
!
!! Horizontal advector for ice volume
!  nt(:,1)=i_dic*i_fracice
!  call bounds(nt(:,1))
!  call mloints(nt(:,1),intsch,nface(:,1),xg(:,1),yg(:,1),2)
!  nt=max(nt,0.)
!  newdic=nt(:,1)/newfracice
!
!! Horizontal advector for snow volume
!  nt(:,1)=i_dsn*i_fracice
!  call bounds(nt(:,1))
!  call mloints(nt(:,1),intsch,nface(:,1),xg(:,1),yg(:,1),2)
!  nt=max(nt,0.)
!  newdsn=nt(:,1)/newfracice
!
!! Horizontal advector for energy store
!  nt(:,1)=i_sto*i_fracice
!  call bounds(nt(:,1))
!  call mloints(nt(:,1),intsch,nface(:,1),xg(:,1),yg(:,1),2)
!  nt=max(nt,0.)
!  newsto=nt(:,1)/newfracice
!
!! Horizontal advector for temperatures 
!  nt(:,1)=i_tsurf*i_fracice
!  call bounds(nt(:,1))
!  call mloints(nt(:,1),intsch,nface(:,1),xg(:,1),yg(:,1),2)
!  nt=max(nt,0.)
!  newtsurf=nt(:,1)/newfracice
!  nt(:,1)=i_t0*i_fracice
!  call bounds(nt(:,1))
!  call mloints(nt(:,1),intsch,nface(:,1),xg(:,1),yg(:,1),2)
!  nt=max(nt,0.)
!  newt0=nt(:,1)/newfracice
!  nt(:,1)=i_t1*i_fracice
!  call bounds(nt(:,1))
!  call mloints(nt(:,1),intsch,nface(:,1),xg(:,1),yg(:,1),2)
!  nt=max(nt,0.)
!  newt1=nt(:,1)/newfracice
!  nt(:,1)=i_t2*i_fracice
!  call bounds(nt(:,1))
!  call mloints(nt(:,1),intsch,nface(:,1),xg(:,1),yg(:,1),2)
!  nt=max(nt,0.)
!  newt2=nt(:,1)/newfracice
!  nt(:,1)=i_t3*i_fracice
!  call bounds(nt(:,1))
!  call mloints(nt(:,1),intsch,nface(:,1),xg(:,1),yg(:,1),2)
!  nt=max(nt,0.)
!  newt3=nt(:,1)/newfracice

!
!end do

return
end subroutine mlohadv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate depature points for MLO semi-Lagrangian advection
! (This subroutine is based on depts.f)

subroutine mlodeps(ubar,vbar,nface,xg,yg,x3d,y3d,z3d)

use cc_mpi
use mlo
use vecsuv_m
use xyzinfo_m

implicit none

include 'newmpar.h'
include 'const_phys.h'
include 'parm.h'

integer ii,intsch,n
integer, dimension(ifull,wlev), intent(out) :: nface
real, dimension(ifull,wlev), intent(in) :: ubar,vbar
real, dimension(ifull,wlev), intent(out) :: xg,yg
real*8, dimension(ifull,wlev), intent(out) :: x3d,y3d,z3d
real, dimension(ifull,wlev) :: uc,vc,wc
real, dimension(ifull+iextra,wlev) :: temp
integer, parameter :: nguess = 2

intsch=mod(ktau,2)

! departure point x, y, z is called x3d, y3d, z3d
! first find corresponding cartesian vels
do ii=1,wlev
  uc(:,ii)=(ax(1:ifull)*ubar(:,ii)+bx(1:ifull)*vbar(:,ii))*dt/rearth ! unit sphere 
  vc(:,ii)=(ay(1:ifull)*ubar(:,ii)+by(1:ifull)*vbar(:,ii))*dt/rearth ! unit sphere 
  wc(:,ii)=(az(1:ifull)*ubar(:,ii)+bz(1:ifull)*vbar(:,ii))*dt/rearth ! unit sphere 
  x3d(:,ii)=x-uc(:,ii) ! 1st guess
  y3d(:,ii)=y-vc(:,ii)
  z3d(:,ii)=z-wc(:,ii)
end do

! convert to grid point numbering
do ii=1,wlev
  call mlotoij5(x3d(:,ii),y3d(:,ii),z3d(:,ii),nface(:,ii),xg(:,ii),yg(:,ii))
end do
! Share off processor departure points.
call deptsync(nface,xg,yg)

do n=1,nguess
  temp(1:ifull,:) = uc
  call mloints(temp,intsch,nface,xg,yg,2)
  do ii=1,wlev
    x3d(:,ii) = x - 0.5*(uc(:,ii)+temp(1:ifull,ii)) ! n+1 guess
  end do
  temp(1:ifull,:) = vc
  call mloints(temp,intsch,nface,xg,yg,2)
  do ii=1,wlev
    y3d(:,ii) = y - 0.5*(vc(:,ii)+temp(1:ifull,ii)) ! n+1 guess
  end do
  temp(1:ifull,:) = wc
  call mloints(temp,intsch,nface,xg,yg,2)
  do ii=1,wlev
    z3d(:,ii) = z - 0.5*(wc(:,ii)+temp(1:ifull,ii)) ! n+1 guess
  end do

  do ii=1,wlev
    call mlotoij5(x3d(:,ii),y3d(:,ii),z3d(:,ii),nface(:,ii),xg(:,ii),yg(:,ii))
  end do
  !     Share off processor departure points.
  call deptsync(nface,xg,yg)
end do

return
end subroutine mlodeps

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate indices
! This code is from depts.f

subroutine mlotoij5(x3d,y3d,z3d,nface,xg,yg)

use bigxy4_m
use cc_mpi
use xyzinfo_m

implicit none

include 'newmpar.h'
include 'parm.h'
include 'parmgeom.h'

integer loop,iq,i,j,is,js
integer, dimension(ifull), intent(out) :: nface
real, dimension(ifull), intent(out) :: xg,yg
real, dimension(ifull) :: xstr,ystr,zstr
real, dimension(ifull) :: denxyz,xd,yd,zd
real, dimension(ifull) :: ri,rj
real dxx,dxy,dyx,dyy
real*8, dimension(ifull), intent(inout) :: x3d,y3d,z3d
real*8, dimension(ifull) :: den
real*8 alf,alfonsch
real*8, parameter :: one = 1.
integer, parameter :: nmaploop = 3

!     if necessary, transform (x3d, y3d, z3d) to equivalent
!     coordinates (xstr, ystr, zstr) on regular gnomonic panels
if(schmidt.eq.1.)then
   xstr=x3d
   ystr=y3d
   zstr=z3d
else      ! (schmidt.ne.1.)
   alf=(one-schmidt**2)/(one+schmidt**2)
   alfonsch=2.*schmidt/(one+schmidt**2)  ! same but bit more accurate
   den=one-alf*z3d ! to force real*8
   xstr=x3d*(alfonsch/den)
   ystr=y3d*(alfonsch/den)
   zstr=    (z3d-alf)/den
endif     ! (schmidt.ne.1.)

!      first deduce departure faces
!      instead calculate cubic coordinates
!      The faces are:
!      0: X=1   1: Z=1   2: Y=1   3: X=-1   4: Z=-1   5: Y=-1
denxyz=max( abs(xstr),abs(ystr),abs(zstr) )
xd=xstr/denxyz
yd=ystr/denxyz
zd=zstr/denxyz

where (abs(xstr-denxyz).lt.1.E-6)
  nface(:)    =0
  xg(:) =      yd
  yg(:) =      zd
elsewhere (abs(xstr+denxyz).lt.1.E-6)
  nface(:)    =3
  xg(:) =     -zd
  yg(:) =     -yd
elsewhere (abs(zstr-denxyz).lt.1.E-6)
  nface(:)    =1
  xg(:) =      yd
  yg(:) =     -xd
elsewhere (abs(zstr+denxyz).lt.1.E-6)
  nface(:)    =4
  xg(:) =      xd
  yg(:) =     -yd
elsewhere (abs(ystr-denxyz).lt.1.E-6)
  nface(:)    =2
  xg(:) =     -zd
  yg(:) =     -xd
elsewhere
  nface(:)    =5
  xg(:) =      xd
  yg(:) =      zd
end where

!     use 4* resolution grid il --> 4*il
xg=min(max(-.99999,xg),.99999)
yg=min(max(-.99999,yg),.99999)
!      first guess for ri, rj and nearest i,j
ri=1.+(1.+xg)*2.*real(il_g)
rj=1.+(1.+yg)*2.*real(il_g)
do loop=1,nmaploop
  do iq=1,ifull
    i=nint(ri(iq))
    j=nint(rj(iq))
    is=nint(sign(1.,ri(iq)-real(i)))
    js=nint(sign(1.,rj(iq)-real(j)))
!       predict new value for ri, rj
    dxx=xx4(i+is,j)-xx4(i,j)
    dyx=xx4(i,j+js)-xx4(i,j)
    dxy=yy4(i+is,j)-yy4(i,j)
    dyy=yy4(i,j+js)-yy4(i,j)       
    den(iq)=dxx*dyy-dyx*dxy
    ri(iq)=real(i)+real(is)*((xg(iq)-xx4(i,j))*dyy-(yg(iq)-yy4(i,j))*dyx)/den(iq)
    rj(iq)=real(j)+real(js)*((yg(iq)-yy4(i,j))*dxx-(xg(iq)-xx4(i,j))*dxy)/den(iq)
  end do
enddo  ! loop loop
!      expect xg, yg to range between .5 and il+.5
xg=.25*(ri+3.) -.5  ! -.5 for stag; back to normal ri, rj defn
yg=.25*(rj+3.) -.5  ! -.5 for stag

return
end subroutine mlotoij5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate depature points for semi-Lagrangian advection
! This code is from ints.f

subroutine mloints(s,intsch,nface,xg,yg,nfield)

use cc_mpi
use indices_m
use mlo

implicit none

include 'newmpar.h'
include 'parm.h'
include 'parmhor.h'
include 'mpif.h'

integer, intent(in) :: nfield,intsch
integer idel,iq,jdel,nn
integer i,j,k,n,ind,ip,jp,iproc,ierr
integer, dimension(ifull,wlev), intent(in) :: nface
real, dimension(ifull,wlev), intent(in) :: xg,yg
real, dimension(ifull+iextra,wlev), intent(inout) :: s
real, dimension(-1:ipan+2,-1:jpan+2,1:npan,wlev) :: sx
real, dimension(4) :: r
real a3,a4,c1,c2,c3,c4,cmax,cmin,sss,xxg,yyg
ind(i,j,n)=i+(j-1)*ipan+(n-1)*ipan*jpan  ! *** for n=1,npan

call bounds(s,nrows=2)

!======================== start of intsch=1 section ====================
if(intsch==1)then

  do k=1,wlev

    do n=1,npan         ! first simple copy into larger array
      do j=1,jpan
        do i=1,ipan
          sx(i,j,n,k) = s(ind(i,j,n),k)
        enddo         ! i loop
      enddo            ! j loop
    enddo               ! n loop
            
!   this is intsb           EW interps done first
!   first extend s arrays into sx - this one -1:il+2 & -1:il+2
    do n=1,npan
      do j=1,jpan
        sx(0,j,n,k) = s(iw(ind(1,j,n)),k)
        sx(-1,j,n,k) = s(iww(ind(1,j,n)),k)
        sx(ipan+1,j,n,k) = s(ie(ind(ipan,j,n)),k)
        sx(ipan+2,j,n,k) = s(iee(ind(ipan,j,n)),k)
      enddo            ! j loop
      do i=1,ipan
        sx(i,0,n,k) = s(is(ind(i,1,n)),k)
        sx(i,-1,n,k) = s(iss(ind(i,1,n)),k)
        sx(i,jpan+1,n,k) = s(in(ind(i,jpan,n)),k)
        sx(i,jpan+2,n,k) = s(inn(ind(i,jpan,n)),k)
      enddo            ! i loop
!        for ew interpolation, sometimes need (different from ns):
!            (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
!           (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)
      sx(-1,0,n,k) = s(lwws(n),k)
      sx(0,0,n,k) = s(lws(n),k)
      sx(0,-1,n,k) = s(lwss(n),k)
      sx(ipan+1,0,n,k) = s(les(n),k)
      sx(ipan+2,0,n,k) = s(lees(n),k)
      sx(ipan+1,-1,n,k) = s(less(n),k)
      sx(-1,jpan+1,n,k) = s(lwwn(n),k)
      sx(0,jpan+2,n,k) = s(lwnn(n),k)
      sx(ipan+2,jpan+1,n,k) = s(leen(n),k)
      sx(ipan+1,jpan+2,n,k) = s(lenn(n),k)
      sx(0,jpan+1,n,k)    = s(iwn(ind(1,jpan,n)),k)
      sx(ipan+1,jpan+1,n,k) = s(ien(ind(ipan,jpan,n)),k)
    enddo               ! n loop

    if(nfield<mh_bs)then
      do iq=1,ifull    ! non Berm-Stan option
!                 Convert face index from 0:npanels to array indices
        idel=int(xg(iq,k))
        xxg=xg(iq,k)-idel
        jdel=int(yg(iq,k))
        yyg=yg(iq,k)-jdel
        ! Now make them proper indices in this processor's region
        idel = idel - ioff(nface(iq,k))
        jdel = jdel - joff(nface(iq,k))
        n = nface(iq,k) + noff ! Make this a local index
        if ( idel < 0 .or. idel > ipan .or. jdel < 0 .or.jdel > jpan .or. n < 1 .or. n > npan ) then
          cycle      ! Will be calculated on another processor
        end if
        c1 = sx(idel-1,jdel,n,k)
        c2 = sx(idel  ,jdel,n,k)
        c3 = sx(idel+1,jdel,n,k)
        c4 = sx(idel+2,jdel,n,k)
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(2) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
        else
          r(2) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.) &
               -xxg*(1.+xxg)*c4/3.)                          &
               +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
        endif         !  (mhint==2)
        c1 = sx(idel-1,jdel+1,n,k)
        c2 = sx(idel  ,jdel+1,n,k)
        c3 = sx(idel+1,jdel+1,n,k)
        c4 = sx(idel+2,jdel+1,n,k)
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(3) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
        else
          r(3) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.) &
               -xxg*(1.+xxg)*c4/3.)                          &
               +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
        endif         !  (mhint==2)
        do nn=1,4,3   ! N.B.
          c2 = sx(idel  ,jdel+nn-2,n,k)
          c3 = sx(idel+1,jdel+nn-2,n,k)
          r(nn) = (1.-xxg)*c2 +xxg*c3
        enddo         ! nn loop
        if(mhint==2)then ! Bessel interp
          a4 = r(4)-r(1)+3.*(r(2)-r(3))
          a3 = r(1)-2.*r(2)+r(3)-a4
          s(iq,k) = r(2)+.5*yyg*(r(3)-r(1) +yyg*(a3+yyg*a4))
        else
          s(iq,k) = ((1.-yyg)*((2.-yyg)*        &
                    ((1.+yyg)*r(2)-yyg*r(1)/3.) &
                    -yyg*(1.+yyg)*r(4)/3.)      &
                    +yyg*(1.+yyg)*(2.-yyg)*r(3))/2.
        endif         !  (mhint==2)
      enddo            ! iq loop
    else                ! (nfield<mh_bs)
      do iq=1,ifull    ! Berm-Stan option here e.g. qg & gases
        idel=int(xg(iq,k))
        xxg=xg(iq,k)-idel
        jdel=int(yg(iq,k))
        yyg=yg(iq,k)-jdel
        ! Now make them proper indices in this processor's region
        idel = idel - ioff(nface(iq,k))
        jdel = jdel - joff(nface(iq,k))
        n = nface(iq,k) + noff ! Make this a local index
        if ( idel < 0 .or. idel > ipan .or. jdel < 0 .or. jdel > jpan .or. n < 1 .or. n > npan ) then
          cycle      ! Will be calculated on another processor
        end if
        c1 = sx(idel-1,jdel,n,k)
        c2 = sx(idel  ,jdel,n,k)
        c3 = sx(idel+1,jdel,n,k)
        c4 = sx(idel+2,jdel,n,k)
        cmin = min( 1.e20,c2,c3) ! Bermejo & Staniforth
        cmax = max(-1.e20,c2,c3) ! Bermejo & Staniforth
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(2) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
        else
          r(2) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.) &
               -xxg*(1.+xxg)*c4/3.)                          &
               +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
        endif         !  (mhint==2)
        c1 = sx(idel-1,jdel+1,n,k)
        c2 = sx(idel  ,jdel+1,n,k)
        c3 = sx(idel+1,jdel+1,n,k)
        c4 = sx(idel+2,jdel+1,n,k)
        cmin = min(cmin,c2,c3) ! Bermejo & Staniforth
        cmax = max(cmax,c2,c3) ! Bermejo & Staniforth
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(3) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
        else
          r(3) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.) &
               -xxg*(1.+xxg)*c4/3.)                          &
               +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
        endif         !  (mhint==2)
        do nn=1,4,3   ! N.B.
          c2 = sx(idel  ,jdel+nn-2,n,k)
          c3 = sx(idel+1,jdel+nn-2,n,k)
          r(nn) = (1.-xxg)*c2 +xxg*c3
        enddo         ! nn loop
        if(mhint==2)then ! Bessel interp
          a4 = r(4)-r(1)+3.*(r(2)-r(3))
          a3 = r(1)-2.*r(2)+r(3)-a4
          sss = r(2)+.5*yyg*(r(3)-r(1) +yyg*(a3+yyg*a4))
        else
          sss = ((1.-yyg)*((2.-yyg)*        &
                ((1.+yyg)*r(2)-yyg*r(1)/3.) &
                 -yyg*(1.+yyg)*r(4)/3.)     &
                 +yyg*(1.+yyg)*(2.-yyg)*r(3))/2.
        endif         !  (mhint==2)
        s(iq,k) = min(max(cmin,sss),cmax) ! Bermejo & Staniforth
      enddo            ! iq loop
    endif               ! (nfield<mh_bs)  .. else ..
            
  end do                 ! k

! Loop over points that need to be calculated for other processes

  if(nfield<mh_bs)then
    do iproc=0,nproc-1
      if ( iproc == myid ) then
        cycle
      end if
      do iq=1,drlen(iproc)
        !  Convert face index from 0:npanels to array indices
        ip = min(il_g,max(1,nint(dpoints(iproc)%a(2,iq))))    ! MJT memory
        jp = min(il_g,max(1,nint(dpoints(iproc)%a(3,iq))))    ! MJT memory
        n = nint(dpoints(iproc)%a(1,iq)) + noff ! Local index ! MJT memory
        !  Need global face index in fproc call
        idel = int(dpoints(iproc)%a(2,iq))                    ! MJT memory
        xxg = dpoints(iproc)%a(2,iq) - idel                   ! MJT memory
        jdel = int(dpoints(iproc)%a(3,iq))                    ! MJT memory
        yyg = dpoints(iproc)%a(3,iq) - jdel                   ! MJT memory
        k = nint(dpoints(iproc)%a(4,iq))                      ! MJT memory
        idel = idel - ioff(n-noff)
        jdel = jdel - joff(n-noff)
        c1 = sx(idel-1,jdel,n,k)
        c2 = sx(idel  ,jdel,n,k)
        c3 = sx(idel+1,jdel,n,k)
        c4 = sx(idel+2,jdel,n,k)
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(2) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
        else
          r(2) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.) &
               -xxg*(1.+xxg)*c4/3.)                          &
               +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
        endif         !  (mhint==2)
        c1 = sx(idel-1,jdel+1,n,k)
        c2 = sx(idel  ,jdel+1,n,k)
        c3 = sx(idel+1,jdel+1,n,k)
        c4 = sx(idel+2,jdel+1,n,k)
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(3) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
        else
          r(3) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.) &
               -xxg*(1.+xxg)*c4/3.)                          &
               +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
        endif         !  (mhint==2)
        do nn=1,4,3   ! N.B.
          c2 = sx(idel  ,jdel+nn-2,n,k)
          c3 = sx(idel+1,jdel+nn-2,n,k)
          r(nn) = (1.-xxg)*c2 +xxg*c3
        enddo         ! nn loop
        if(mhint==2)then ! Bessel interp
          a4 = r(4)-r(1)+3.*(r(2)-r(3))
          a3 = r(1)-2.*r(2)+r(3)-a4
          sextra(iproc)%a(iq) = r(2) + 0.5*yyg*(r(3)-r(1) +yyg*(a3+yyg*a4))
        else
          sextra(iproc)%a(iq) = ((1.-yyg)*((2.-yyg)* &        ! MJT memory
            ((1.+yyg)*r(2)-yyg*r(1)/3.)              &
            -yyg*(1.+yyg)*r(4)/3.)                   &
            +yyg*(1.+yyg)*(2.-yyg)*r(3))/2.
        endif         !  (mhint==2)
      enddo            ! iq loop
    end do              ! iproc loop
  else                   ! (nfield<mh_bs)
    do iproc=0,nproc-1
      if ( iproc == myid ) then
        cycle
      end if
      do iq=1,drlen(iproc)
        !  Convert face index from 0:npanels to array indices
        ip = min(il_g,max(1,nint(dpoints(iproc)%a(2,iq))))    ! MJT memory
        jp = min(il_g,max(1,nint(dpoints(iproc)%a(3,iq))))    ! MJT memory
        n = nint(dpoints(iproc)%a(1,iq)) + noff ! Local index ! MJT memory
        !  Need global face index in fproc call
        idel = int(dpoints(iproc)%a(2,iq))                    ! MJT memory
        xxg = dpoints(iproc)%a(2,iq) - idel                   ! MJT memory
        jdel = int(dpoints(iproc)%a(3,iq))                    ! MJT memory
        yyg = dpoints(iproc)%a(3,iq) - jdel                   ! MJT memory
        k = nint(dpoints(iproc)%a(4,iq))                      ! MJT memory
        idel = idel - ioff(n-noff)
        jdel = jdel - joff(n-noff)
        c1 = sx(idel-1,jdel,n,k)
        c2 = sx(idel  ,jdel,n,k)
        c3 = sx(idel+1,jdel,n,k)
        c4 = sx(idel+2,jdel,n,k)
        cmin = min( 1.e20,c2,c3) ! Bermejo & Staniforth
        cmax = max(-1.e20,c2,c3) ! Bermejo & Staniforth
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(2) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
        else
          r(2) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.) &
               -xxg*(1.+xxg)*c4/3.)                          &
               +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
        endif         !  (mhint==2)
        c1 = sx(idel-1,jdel+1,n,k)
        c2 = sx(idel  ,jdel+1,n,k)
        c3 = sx(idel+1,jdel+1,n,k)
        c4 = sx(idel+2,jdel+1,n,k)
        cmin = min(cmin,c2,c3) ! Bermejo & Staniforth
        cmax = max(cmax,c2,c3) ! Bermejo & Staniforth
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(3) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
        else
          r(3) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.) &
               -xxg*(1.+xxg)*c4/3.)                          &
               +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
        endif         !  (mhint==2)
        do nn=1,4,3   ! N.B.
          c2 = sx(idel  ,jdel+nn-2,n,k)
          c3 = sx(idel+1,jdel+nn-2,n,k)
          r(nn) = (1.-xxg)*c2 +xxg*c3
        enddo         ! nn loop
        if(mhint==2)then ! Bessel interp
          a4 = r(4)-r(1)+3.*(r(2)-r(3))
          a3 = r(1)-2.*r(2)+r(3)-a4
          sss = r(2)+.5*yyg*(r(3)-r(1) +yyg*(a3+yyg*a4))
        else
          sss = ((1.-yyg)*((2.-yyg)*        &
                ((1.+yyg)*r(2)-yyg*r(1)/3.) &
                -yyg*(1.+yyg)*r(4)/3.)      &
                +yyg*(1.+yyg)*(2.-yyg)*r(3))/2.
        endif         !  (mhint==2)
        sextra(iproc)%a(iq) = min(max(cmin,sss),cmax) ! Bermejo & Staniforth ! MJT memory
      enddo            ! iq loop
    end do              ! iproc loop
  endif                  ! (nfield<mh_bs)  .. else ..
!========================   end of intsch=1 section ====================
else     ! if(intsch==1)then
!======================== start of intsch=2 section ====================
  
!       this is intsc           NS interps done first
!       first extend s arrays into sx - this one -1:il+2 & -1:il+2
  do k=1,wlev              ! A single main k loop uses cache better
    do n=1,npan         ! first simple copy into larger array
      do j=1,jpan
        do i=1,ipan
          sx(i,j,n,k)=s(ind(i,j,n),k)
        enddo         ! i loop
      enddo            ! j loop
    enddo               ! n loop

    do n=1,npan
      do j=1,jpan
        sx(0,j,n,k) = s(iw(ind(1,j,n)),k)
        sx(-1,j,n,k) = s(iww(ind(1,j,n)),k)
        sx(ipan+1,j,n,k) = s(ie(ind(ipan,j,n)),k)
        sx(ipan+2,j,n,k) = s(iee(ind(ipan,j,n)),k)
      enddo            ! j loop
      do i=1,ipan
        sx(i,0,n,k) = s(is(ind(i,1,n)),k)
        sx(i,-1,n,k) = s(iss(ind(i,1,n)),k)
        sx(i,jpan+1,n,k) = s(in(ind(i,jpan,n)),k)
        sx(i,jpan+2,n,k) = s(inn(ind(i,jpan,n)),k)
      enddo            ! i loop
!        for ns interpolation, sometimes need (different from ew):
!            (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
!          (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)
      sx(-1,0,n,k)=s(lsww(n),k)
      sx(0,0,n,k) = s(lsw(n),k)
      sx(0,-1,n,k) = s(lssw(n),k)
      sx(ipan+2,0,n,k) = s(lsee(n),k)
      sx(ipan+1,-1,n,k) = s(lsse(n),k)
      sx(-1,jpan+1,n,k) = s(lnww(n),k)
      sx(0,jpan+1,n,k) = s(lnw(n),k)
      sx(0,jpan+2,n,k) = s(lnnw(n),k)
      sx(ipan+2,jpan+1,n,k) = s(lnee(n),k)
      sx(ipan+1,jpan+2,n,k) = s(lnne(n),k)
      sx(ipan+1,0,n,k)    = s(ise(ind(ipan,1,n)),k)
      sx(ipan+1,jpan+1,n,k) = s(ine(ind(ipan,jpan,n)),k)
    enddo               ! n loop

    if(nfield<mh_bs)then
      do iq=1,ifull    ! non Berm-Stan option
!       Convert face index from 0:npanels to array indices
        idel=int(xg(iq,k))
        xxg=xg(iq,k)-idel
        jdel=int(yg(iq,k))
        yyg=yg(iq,k)-jdel
        ! Now make them proper indices in this processor's region
        idel = idel - ioff(nface(iq,k))
        jdel = jdel - joff(nface(iq,k))
        n = nface(iq,k) + noff ! Make this a local index
        if ( idel < 0 .or. idel > ipan .or. jdel < 0 .or. jdel > jpan .or. n < 1 .or. n > npan ) then
          cycle      ! Will be calculated on another processor
        end if
        c1 = sx(idel,jdel-1,n,k)
        c2 = sx(idel,jdel  ,n,k)
        c3 = sx(idel,jdel+1,n,k)
        c4 = sx(idel,jdel+2,n,k)
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(2) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
        else
          r(2) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.) &
               -yyg*(1.+yyg)*c4/3.)                          &
               +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
        endif         !  (mhint==2)
        c1 = sx(idel+1,jdel-1,n,k)
        c2 = sx(idel+1,jdel  ,n,k)
        c3 = sx(idel+1,jdel+1,n,k)
        c4 = sx(idel+1,jdel+2,n,k)
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(3) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
        else
          r(3) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.) &
               -yyg*(1.+yyg)*c4/3.)                          &
               +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
        endif         !  (mhint==2)
        do nn=1,4,3   ! N.B.
          c2 = sx(idel+nn-2,jdel  ,n,k)
          c3 = sx(idel+nn-2,jdel+1,n,k)
          r(nn) = (1.-yyg)*c2 +yyg*c3
        enddo         ! nn loop
        if(mhint==2)then ! Bessel interp
          a4 = r(4)-r(1)+3.*(r(2)-r(3))
          a3 = r(1)-2.*r(2)+r(3)-a4
          s(iq,k) = r(2)+.5*xxg*(r(3)-r(1) +xxg*(a3+xxg*a4))
        else
          s(iq,k) = ((1.-xxg)*((2.-xxg)*        &
                    ((1.+xxg)*r(2)-xxg*r(1)/3.) &
                   -xxg*(1.+xxg)*r(4)/3.)       &
                   +xxg*(1.+xxg)*(2.-xxg)*r(3))/2.
        endif         !  (mhint==2)
      enddo            ! iq loop
    else                ! (nfield<mh_bs)
      do iq=1,ifull    ! Berm-Stan option here e.g. qg & gases
        idel=int(xg(iq,k))
        xxg=xg(iq,k)-idel
        jdel=int(yg(iq,k))
        yyg=yg(iq,k)-jdel
        ! Now make them proper indices in this processor's region
        idel = idel - ioff(nface(iq,k))
        jdel = jdel - joff(nface(iq,k))
        n = nface(iq,k) + noff ! Make this a local index
        if ( idel < 0 .or. idel > ipan .or. jdel < 0 .or. jdel > jpan .or. n < 1 .or. n > npan ) then
          cycle      ! Will be calculated on another processor
        end if
        c1 = sx(idel,jdel-1,n,k)
        c2 = sx(idel,jdel  ,n,k)
        c3 = sx(idel,jdel+1,n,k)
        c4 = sx(idel,jdel+2,n,k)
        cmin = min( 1.e20,c2,c3) ! Bermejo & Staniforth
        cmax = max(-1.e20,c2,c3) ! Bermejo & Staniforth
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(2) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
        else
          r(2) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.) &
               -yyg*(1.+yyg)*c4/3.)                          &
               +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
        endif         !  (mhint==2)
        c1 = sx(idel+1,jdel-1,n,k)
        c2 = sx(idel+1,jdel  ,n,k)
        c3 = sx(idel+1,jdel+1,n,k)
        c4 = sx(idel+1,jdel+2,n,k)
        cmin = min(cmin,c2,c3) ! Bermejo & Staniforth
        cmax = max(cmax,c2,c3) ! Bermejo & Staniforth
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(3) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
        else
          r(3) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.) &
               -yyg*(1.+yyg)*c4/3.)                          &
               +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
        endif         !  (mhint==2)
        do nn=1,4,3   ! N.B.
          c2 = sx(idel+nn-2,jdel  ,n,k)
          c3 = sx(idel+nn-2,jdel+1,n,k)
          r(nn) = (1.-yyg)*c2 +yyg*c3
        enddo         ! nn loop
        if(mhint==2)then ! Bessel interp
          a4 = r(4)-r(1)+3.*(r(2)-r(3))
          a3 = r(1)-2.*r(2)+r(3)-a4
          sss = r(2)+.5*xxg*(r(3)-r(1) +xxg*(a3+xxg*a4))
        else
          sss = ((1.-xxg)*((2.-xxg)*        &
                ((1.+xxg)*r(2)-xxg*r(1)/3.) &
                -xxg*(1.+xxg)*r(4)/3.)      &
                +xxg*(1.+xxg)*(2.-xxg)*r(3))/2.
        endif         !  (mhint==2)
        s(iq,k) = min(max(cmin,sss),cmax) ! Bermejo & Staniforth
      enddo            ! iq loop
    endif               ! (nfield<mh_bs)  .. else ..
  end do                 ! k

! For other processes
  if(nfield<mh_bs)then
    do iproc=0,nproc-1
      if ( iproc == myid ) then
        cycle
      end if
      do iq=1,drlen(iproc)
        !  Convert face index from 0:npanels to array indices
        ip = min(il_g,max(1,nint(dpoints(iproc)%a(2,iq))))    ! MJT memory
        jp = min(il_g,max(1,nint(dpoints(iproc)%a(3,iq))))    ! MJT memory
        n = nint(dpoints(iproc)%a(1,iq)) + noff ! Local index ! MJT memory
        !  Need global face index in fproc call
        idel = int(dpoints(iproc)%a(2,iq))                    ! MJT memory
        xxg = dpoints(iproc)%a(2,iq) - idel                   ! MJT memory
        jdel = int(dpoints(iproc)%a(3,iq))                    ! MJT memory
        yyg = dpoints(iproc)%a(3,iq) - jdel                   ! MJT memory
        k = nint(dpoints(iproc)%a(4,iq))                      ! MJT memory
        idel = idel - ioff(n-noff)
        jdel = jdel - joff(n-noff)
        c1 = sx(idel,jdel-1,n,k)
        c2 = sx(idel,jdel  ,n,k)
        c3 = sx(idel,jdel+1,n,k)
        c4 = sx(idel,jdel+2,n,k)
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(2) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
        else
          r(2) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.) &
               -yyg*(1.+yyg)*c4/3.)                          &
               +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
        endif         !  (mhint==2)
        c1 = sx(idel+1,jdel-1,n,k)
        c2 = sx(idel+1,jdel  ,n,k)
        c3 = sx(idel+1,jdel+1,n,k)
        c4 = sx(idel+1,jdel+2,n,k)
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(3) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
        else
          r(3) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.) &
               -yyg*(1.+yyg)*c4/3.)                          &
               +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
        endif         !  (mhint==2)
        do nn=1,4,3   ! N.B.
          c2 = sx(idel+nn-2,jdel  ,n,k)
          c3 = sx(idel+nn-2,jdel+1,n,k)
          r(nn) = (1.-yyg)*c2 +yyg*c3
        enddo         ! nn loop
        if(mhint==2)then ! Bessel interp
          a4 = r(4)-r(1)+3.*(r(2)-r(3))
          a3 = r(1)-2.*r(2)+r(3)-a4
          sextra(iproc)%a(iq) = r(2)+ 0.5*xxg*(r(3)-r(1) +xxg*(a3+xxg*a4))
        else
          sextra(iproc)%a(iq) = ((1.-xxg)*((2.-xxg)*  & ! MJT memory
            ((1.+xxg)*r(2)-xxg*r(1)/3.)               &
            -xxg*(1.+xxg)*r(4)/3.)                    &
            +xxg*(1.+xxg)*(2.-xxg)*r(3))/2.
        endif         !  (mhint==2)
      enddo            ! iq loop
    end do              ! iproc
  else                   ! (nfield<mh_bs)
    do iproc=0,nproc-1
      if ( iproc == myid ) then
        cycle
      end if
      do iq=1,drlen(iproc)
        !  Convert face index from 0:npanels to array indices
        ip = min(il_g,max(1,nint(dpoints(iproc)%a(2,iq))))    ! MJT memory
        jp = min(il_g,max(1,nint(dpoints(iproc)%a(3,iq))))    ! MJT memory
        n = nint(dpoints(iproc)%a(1,iq)) + noff ! Local index ! MJT memory
        !  Need global face index in fproc call
        idel = int(dpoints(iproc)%a(2,iq))                    ! MJT memory
        xxg = dpoints(iproc)%a(2,iq) - idel                   ! MJT memory
        jdel = int(dpoints(iproc)%a(3,iq))                    ! MJT memory
        yyg = dpoints(iproc)%a(3,iq) - jdel                   ! MJT memory
        k = nint(dpoints(iproc)%a(4,iq))                      ! MJT memory
        idel = idel - ioff(n-noff)
        jdel = jdel - joff(n-noff)
        c1 = sx(idel,jdel-1,n,k)
        c2 = sx(idel,jdel  ,n,k)
        c3 = sx(idel,jdel+1,n,k)
        c4 = sx(idel,jdel+2,n,k)
        cmin = min( 1.e20,c2,c3) ! Bermejo & Staniforth
        cmax = max(-1.e20,c2,c3) ! Bermejo & Staniforth
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(2) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
        else
          r(2) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.) &
               -yyg*(1.+yyg)*c4/3.)                          &
               +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
        endif         !  (mhint==2)
        c1 = sx(idel+1,jdel-1,n,k)
        c2 = sx(idel+1,jdel  ,n,k)
        c3 = sx(idel+1,jdel+1,n,k)
        c4 = sx(idel+1,jdel+2,n,k)
        cmin = min(cmin,c2,c3) ! Bermejo & Staniforth
        cmax = max(cmax,c2,c3) ! Bermejo & Staniforth
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(3) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
        else
          r(3) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.) &
               -yyg*(1.+yyg)*c4/3.)                          &
               +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
        endif         !  (mhint==2)
        do nn=1,4,3   ! N.B.
          c2 = sx(idel+nn-2,jdel  ,n,k)
          c3 = sx(idel+nn-2,jdel+1,n,k)
          r(nn) = (1.-yyg)*c2 +yyg*c3
        enddo         ! nn loop
        if(mhint==2)then ! Bessel interp
          a4 = r(4)-r(1)+3.*(r(2)-r(3))
          a3 = r(1)-2.*r(2)+r(3)-a4
          sss = r(2)+.5*xxg*(r(3)-r(1) +xxg*(a3+xxg*a4))
        else
          sss = ((1.-xxg)*((2.-xxg)*        &
                ((1.+xxg)*r(2)-xxg*r(1)/3.) &
                -xxg*(1.+xxg)*r(4)/3.)      &
                +xxg*(1.+xxg)*(2.-xxg)*r(3))/2.
        endif         !  (mhint==2)
        sextra(iproc)%a(iq) = min(max(cmin,sss),cmax) ! Bermejo & Staniforth ! MJT memory
      enddo            ! iq loop
    end do              ! iproc
  endif                  ! (nfield<mh_bs)  .. else ..

endif                     ! (intsch==1) .. else ..
!========================   end of intsch=1 section ====================

call intssync(s)

return
end subroutine mloints

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Rotate wind vector to arrival point

subroutine mlorot(cou,cov,cow,x3d,y3d,z3d)

use mlo
use xyzinfo_m

implicit none

include 'newmpar.h'

integer k
real, dimension(ifull,wlev), intent(inout) :: cou,cov,cow
real, dimension(ifull) :: vec1x,vec1y,vec1z,denb
real, dimension(ifull) :: vec2x,vec2y,vec2z,vecdot
real, dimension(ifull) :: vec3x,vec3y,vec3z,vdot1,vdot2
real*8, dimension(ifull,wlev), intent(in) :: x3d,y3d,z3d

do k=1,wlev
!         cross product n1xn2 into vec1
  vec1x = y3d(:,k)*z - y*z3d(:,k)
  vec1y = z3d(:,k)*x - z*x3d(:,k)
  vec1z = x3d(:,k)*y - x*y3d(:,k)
  denb = vec1x**2 + vec1y**2 + vec1z**2
!         N.B. rotation formula is singular for small denb,
!         but the rotation is unnecessary in this case
  where (denb>1.e-4)
    vecdot = x3d(:,k)*x + y3d(:,k)*y + z3d(:,k)*z
    vec2x = x3d(:,k)*vecdot - x
    vec2y = y3d(:,k)*vecdot - y
    vec2z = z3d(:,k)*vecdot - z
    vec3x = x3d(:,k) - vecdot*x
    vec3y = y3d(:,k) - vecdot*y
    vec3z = z3d(:,k) - vecdot*z
    vdot1 = (vec1x*cou(:,k) + vec1y*cov(:,k) + vec1z*cow(:,k))/denb
    vdot2 = (vec2x*cou(:,k) + vec2y*cov(:,k) + vec2z*cow(:,k))/denb
    cou(:,k) = vdot1*vec1x + vdot2*vec3x
    cov(:,k) = vdot1*vec1y + vdot2*vec3y
    cow(:,k) = vdot1*vec1z + vdot2*vec3z
  end where
end do ! k

return
end subroutine mlorot


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Stagger u and v

subroutine mlostaguv(u,v,uout,vout)

use cc_mpi
use indices_m
use mlo

implicit none

include 'newmpar.h'

integer k,itn
real, dimension(ifull,wlev), intent(in) :: u,v
real, dimension(ifull,wlev), intent(out) :: uout,vout
real, dimension(ifull+iextra,wlev) :: uin,vin
real, dimension(ifull+iextra,wlev) :: ua,va,ud,vd
integer, parameter :: itnmax=3

uin(1:ifull,:)=u
vin(1:ifull,:)=v
call boundsuv(uin,vin)

do k=1,wlev
  ud(1:ifull,k)= uin(iwu,k)/10.+uin(1:ifull,k)+uin(ieu,k)/2.
  vd(1:ifull,k)= vin(isv,k)/10.+vin(1:ifull,k)+vin(inv,k)/2.
enddo
call boundsuv(ud,vd)

do k=1,wlev
  ua(1:ifull,k)=ud(1:ifull,k)-ud(iwu,k)/2. ! 1st guess
  va(1:ifull,k)=vd(1:ifull,k)-vd(isv,k)/2. ! 1st guess
enddo

do itn=1,itnmax        ! each loop is a double iteration
  call boundsuv(ua,va,nrows=2)

  do k=1,wlev
    uin(1:ifull,k)=(ud(1:ifull,k)-.5*ud(iwu,k)-ua(ieu,k)/10. +ua(iwwu,k)/4.)/.95
    vin(1:ifull,k)=(vd(1:ifull,k)-.5*vd(isv,k)-va(inv,k)/10. +va(issv,k)/4.)/.95
  enddo

  call boundsuv(uin,vin,nrows=2)
  do k=1,wlev
    ua(1:ifull,k)=(ud(1:ifull,k)-.5*ud(iwu,k)-uin(ieu,k)/10. +uin(iwwu,k)/4.)/.95
    va(1:ifull,k)=(vd(1:ifull,k)-.5*vd(isv,k)-vin(inv,k)/10. +vin(issv,k)/4.)/.95
  end do
end do                 ! itn=1,itnmax

uout=ua(1:ifull,:)
vout=va(1:ifull,:)

return
end subroutine mlostaguv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Unstagger u and v

subroutine mlounstaguv(u,v,uout,vout)

use cc_mpi
use indices_m
use mlo

implicit none

include 'newmpar.h'

integer k,itn
real, dimension(ifull,wlev), intent(in) :: u,v
real, dimension(ifull,wlev), intent(out) :: uout,vout
real, dimension(ifull+iextra,wlev) :: uin,vin
real, dimension(ifull+iextra,wlev) :: ua,va,ud,vd
integer, parameter :: itnmax=3

uin(1:ifull,:)=u
vin(1:ifull,:)=v
call boundsuv(uin,vin)

do k=1,wlev
  ud(1:ifull,k)= uin(ieu,k)/10.+uin(1:ifull,k)+uin(iwu,k)/2.
  vd(1:ifull,k)= vin(inv,k)/10.+vin(1:ifull,k)+vin(isv,k)/2.
enddo
call boundsuv(ud,vd)

do k=1,wlev
  ua(1:ifull,k)=ud(1:ifull,k)-ud(ieu,k)/2. ! 1st guess
  va(1:ifull,k)=vd(1:ifull,k)-vd(inv,k)/2. ! 1st guess
enddo

do itn=1,itnmax        ! each loop is a double iteration
  call boundsuv(ua,va,nrows=2)

  do k=1,wlev
    uin(1:ifull,k)=(ud(1:ifull,k)-.5*ud(ieu,k)-ua(iwu,k)/10. +ua(ieeu,k)/4.)/.95
    vin(1:ifull,k)=(vd(1:ifull,k)-.5*vd(inv,k)-va(isv,k)/10. +va(innv,k)/4.)/.95
  enddo
  call boundsuv(uin,vin,nrows=2)
  do k=1,wlev
    ua(1:ifull,k)=(ud(1:ifull,k)-.5*ud(ieu,k)-uin(iwu,k)/10. +uin(ieeu,k)/4.)/.95
    va(1:ifull,k)=(vd(1:ifull,k)-.5*vd(inv,k)-vin(isv,k)/10. +vin(innv,k)/4.)/.95
  enddo
enddo                  ! itn=1,itnmax
      
uout=ua(1:ifull,:)
vout=va(1:ifull,:)

return
end subroutine mlounstaguv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates horizontal gradients for z-sigma
! coordinates

!subroutine truedelta(s,dep,wtr,dsdx,dsdy,d2sdx2,d2sdy2)
!
!use map_m
!use mlo
!use indices_m
!
!implicit none
!
!include 'newmpar.h'
!include 'parm.h'
!
!integer iq,ii,iqn,iqe,iqs,iqw
!integer idn,ide,ids,idw
!real sx,sn,se,ss,sw,ddx
!real, dimension(wlev) :: ddn,dde,dds,ddw
!real, dimension(ifull+iextra,wlev), intent(in) :: s,dep
!real, dimension(ifull,wlev), intent(out) :: dsdx,dsdy,d2sdx2,d2sdy2
!logical, dimension(ifull+iextra), intent(in) :: wtr
!
!dsdx=0.
!dsdy=0.
!d2sdx2=0.
!d2sdy2=0.
!
!do iq=1,ifull
!  if (wtr(iq)) then
!    iqn=in(iq)
!    iqe=ie(iq)
!    iqs=is(iq)
!    iqw=iw(iq)
!  
!    idn=1
!    ide=1
!    ids=1
!    idw=1
!    
!    ! now search for other levels
!    do ii=1,wlev
!      sx=s(iq,ii)
!      ddx=dep(iq,ii)
!      ddn=dep(iqn,:)
!      dde=dep(iqe,:)
!      dds=dep(iqs,:)
!      ddw=dep(iqw,:)
!    
!      call searchdelta(s(iqn,:),ddn,wtr(iqn),ddx,sx,idn,sn)
!      call searchdelta(s(iqe,:),dde,wtr(iqe),ddx,sx,ide,se)
!      call searchdelta(s(iqs,:),dds,wtr(iqs),ddx,sx,ids,ss)
!      call searchdelta(s(iqw,:),ddw,wtr(iqw),ddx,sx,idw,sw)
!    
!      dsdx(iq,ii)=0.5*(se-sw)*em(iq)/ds
!      d2sdx2(iq,ii)=(se-2.*sx+sw)*em(iq)*em(iq)/(ds*ds)
!      dsdy(iq,ii)=(sn-ss)*em(iq)*em(iq)/(ds*ds)
!      d2sdy2(iq,ii)=(sn-2.*sx+ss)*em(iq)*em(iq)/(ds*ds)
!    end do
!  end if
!end do
!
!return
!end subroutine truedelta
!
!subroutine stagtruedeltau(au,bu,cu,dep,wtr,daudx,dbudx,dcudx)
!
!use map_m
!use mlo
!use indices_m
!
!implicit none
!
!include 'newmpar.h'
!include 'parm.h'
!
!integer iq,iqe,iqw,iqwu,ide,idw,ii
!real sxau,sxbu,sxcu,ddx,seau,sebu,secu,swau,swbu,swcu
!real, dimension(wlev) :: dde,ddw
!real, dimension(ifull+iextra,wlev), intent(in) :: au,bu,cu,dep
!real, dimension(ifull,wlev), intent(out) :: daudx,dbudx,dcudx
!logical, dimension(ifull+iextra), intent(in) :: wtr
!
!daudx=0.
!dbudx=0.
!dcudx=0.
!
!do iq=1,ifull
!  if (wtr(iq)) then
!    iqe=ie(iq)
!    iqw=iw(iq)
!    iqwu=iwu(iq)
!
!    ide=1  
!    idw=1
!    ! now search for other levels
!    do ii=1,wlev
!      sxau=au(iq,ii)
!      sxbu=bu(iq,ii)
!      sxcu=cu(iq,ii)
!      
!      ddx=dep(iq,ii)
!      dde=0.5*(dep(iqe,:)+ddx)
!      ddw=0.5*(dep(iqw,:)+ddx)
!    
!      call searchdelta(au(iq,:),dde,wtr(iqe),ddx,0.,ide,seau)
!      call searchdelta(bu(iq,:),dde,wtr(iqe),ddx,0.,ide,sebu)
!      call searchdelta(cu(iq,:),dde,wtr(iqe),ddx,0.,ide,secu)
!      call searchdelta(au(iqwu,:),ddw,wtr(iqw),ddx,0.,idw,swau)
!      call searchdelta(bu(iqwu,:),ddw,wtr(iqw),ddx,0.,idw,swbu)
!      call searchdelta(cu(iqwu,:),ddw,wtr(iqw),ddx,0.,idw,swcu)
!    
!      daudx(iq,ii)=(seau-swau)*em(iq)/ds
!      dbudx(iq,ii)=(sebu-swbu)*em(iq)/ds
!      dcudx(iq,ii)=(secu-swcu)*em(iq)/ds
!    end do
!  end if
!end do
!
!return
!end subroutine stagtruedeltau
!
!subroutine stagtruedeltav(av,bv,cv,dep,wtr,davdy,dbvdy,dcvdy)
!
!use map_m
!use mlo
!use indices_m
!
!implicit none
!
!include 'newmpar.h'
!include 'parm.h'
!
!integer iq,iqn,iqs,iqsv,idn,ids,ii
!real sxav,sxbv,sxcv,ddx,snav,snbv,sncv,ssav,ssbv,sscv
!real, dimension(wlev) :: ddn,dds
!real, dimension(ifull+iextra,wlev), intent(in) :: av,bv,cv,dep
!real, dimension(ifull,wlev), intent(out) :: davdy,dbvdy,dcvdy
!logical, dimension(ifull+iextra), intent(in) :: wtr
!
!davdy=0.
!dbvdy=0.
!dcvdy=0.
!
!do iq=1,ifull
!  if (wtr(iq)) then
!    iqn=in(iq)
!    iqs=is(iq)
!    iqsv=isv(iq)
!
!    idn=1  
!    ids=1
!    ! now search for other levels
!    do ii=1,wlev
!      sxav=av(iq,ii)
!      sxbv=bv(iq,ii)
!      sxcv=cv(iq,ii)
!      
!      ddx=dep(iq,ii)
!      ddn=0.5*(dep(iqn,:)+ddx)
!      dds=0.5*(dep(iqs,:)+ddx)
!    
!      call searchdelta(av(iq,:),ddn,wtr(iqn),ddx,0.,idn,snav)
!      call searchdelta(bv(iq,:),ddn,wtr(iqn),ddx,0.,idn,snbv)
!      call searchdelta(cv(iq,:),ddn,wtr(iqn),ddx,0.,idn,sncv)
!      call searchdelta(av(iqsv,:),dds,wtr(iqs),ddx,0.,ids,ssav)
!      call searchdelta(bv(iqsv,:),dds,wtr(iqs),ddx,0.,ids,ssbv)
!      call searchdelta(cv(iqsv,:),dds,wtr(iqs),ddx,0.,ids,sscv)
!    
!      davdy(iq,ii)=(snav-ssav)*em(iq)/ds
!      dbvdy(iq,ii)=(snbv-ssbv)*em(iq)/ds
!      dcvdy(iq,ii)=(sncv-sscv)*em(iq)/ds
!    end do
!  end if
!end do
!
!return
!end subroutine stagtruedeltav
!
!subroutine stagtruedelta(s,dep,wtr,dsdxu,dsdyu,dsdxv,dsdyv)
!
!use map_m
!use mlo
!use indices_m
!
!implicit none
!
!include 'newmpar.h'
!include 'parm.h'
!
!integer iq,ii
!integer idn,ide,ids,idw,idne,idse,iden,idwn
!real sx,sn,se,ss,sw,sne,sse,sen,swn,ddx
!real, dimension(wlev) :: ddn,dde,dds,ddw,ddne,ddse,dden,ddwn
!real, dimension(ifull+iextra,wlev), intent(in) :: s,dep
!real, dimension(ifull,wlev), intent(out) :: dsdxu,dsdyu,dsdxv,dsdyv
!logical, dimension(ifull+iextra), intent(in) :: wtr
!
!dsdxu=0.
!dsdyu=0.
!dsdxv=0.
!dsdyv=0.
!
!do iq=1,ifull
!  if (wtr(iq)) then
!
!    idn=1
!    ide=1
!    ids=1
!    idw=1
!    idne=1
!    idse=1
!    iden=1
!    idwn=1
!
!    ! now search for other levels
!    do ii=1,wlev
!      sx=s(iq,ii)
!      ddx=dep(iq,ii)
!      ddn=dep(in(iq),:)
!      ddne=dep(ine(iq),:)
!      dds=dep(is(iq),:)
!      ddse=dep(ise(iq),:)
!      dde=dep(ie(iq),:)
!      dden=dep(ien(iq),:)
!      ddw=dep(iw(iq),:)
!      ddwn=dep(iwn(iq),:)
!
!    
!      call searchdelta(s(in(iq),:),ddn,wtr(in(iq)),ddx,sx,idn,sn)
!      call searchdelta(s(ine(iq),:),ddne,wtr(ine(iq)),ddx,sn,idne,sne)
!      call searchdelta(s(is(iq),:),dds,wtr(is(iq)),ddx,sx,ids,ss)
!      call searchdelta(s(ise(iq),:),ddse,wtr(ise(iq)),ddx,ss,idse,sse)
!      call searchdelta(s(ie(iq),:),dde,wtr(ie(iq)),ddx,sx,ide,se)
!      call searchdelta(s(ien(iq),:),dden,wtr(ien(iq)),ddx,se,iden,sen)
!      call searchdelta(s(iw(iq),:),ddw,wtr(iw(iq)),ddx,sx,idw,sw)
!      call searchdelta(s(iwn(iq),:),ddwn,wtr(iwn(iq)),ddx,sw,idwn,swn)
!    
!      call fourpoint(sx,sn,sne,sen,se,sse,ss,sw,swn,emu(iq),emv(iq),ds,dsdxu(iq,ii),dsdxv(iq,ii),dsdyu(iq,ii),dsdyv(iq,ii))
!    end do
!  end if
!end do
!
!return
!end subroutine stagtruedelta
!
!subroutine fourpoint(sx,sn,sne,sen,se,sse,ss,sw,swn,emu,emv,ds,dsdxu,dsdxv,dsdyu,dsdyv)
!
!implicit none
!
!real, intent(in) :: sx,sn,sne,sen,se,sse,ss,sw,swn,emu,emv,ds
!real, intent(out) :: dsdxu,dsdyu,dsdxv,dsdyv
!real tnu,tsu,tev,twv
!
!tnu=0.5*(sn+sne) 
!tsu=0.5*(ss+sse) 
!tev=0.5*(se+sen) 
!twv=0.5*(sw+swn) 
! 
!dsdxu=(se-sx)*emu/ds
!dsdyu=0.5*(tnu-tsu)*emu/ds
!dsdxv=0.5*(tev-twv)*emv/ds
!dsdyv=(sn-sx)*emv/ds
!
!return
!end subroutine fourpoint
!
!subroutine searchdelta(s,dep,wtr,dd,sx,id,ss)
!
!use mlo
!
!implicit none
!
!integer, intent(inout) :: id
!integer jj,fnd
!real, dimension(wlev), intent(in) :: s,dep
!real, intent(in) :: dd,sx
!real, intent(out) :: ss
!real xp
!logical, intent(in) :: wtr
!
!ss=sx
!
!if (.not.wtr) return
!
!if (dep(id).gt.dd+0.01) return
!
!fnd=-999
!do jj=id,wlev-1
!  if (dep(jj+1).gt.dd+0.01) then
!    fnd=jj
!    exit
!  end if
!end do
!if (fnd.lt.1) return
!
!! located a valid level to interpolate
!id=fnd
!xp=(dd-dep(id))/(dep(id+1)-dep(id))
!xp=max(min(xp,1.),0.)
!ss=(1.-xp)*s(id)+xp*s(id+1)
!
!return
!end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Tide

subroutine mlotide(eta,slon,slat,mtimer,jstart)

implicit none

include 'newmpar.h'

integer i
integer, intent(in) :: mtimer,jstart
real, dimension(ifull), intent(out) :: eta
real, dimension(ifull), intent(in) :: slon,slat
real :: stime
real, dimension(ifull) :: sinlon,coslon,sinlat,coslat
real, dimension(4) :: aa,ab,ba,bb,wa,wb
real, parameter :: pi = 3.1415927

aa(1)=0.141565
aa(2)=0.100661
aa(3)=0.046848
aa(4)=0.019273
ab(1)=0.242334
ab(2)=0.112743
ab(3)=0.046397
ab(4)=0.030684
ba(1)=0.736
ba(2)=0.695
ba(3)=0.706
ba(4)=0.695
bb(1)=0.693
bb(2)=0.693
bb(3)=0.693
bb(4)=0.693
wa(1)=0.7292117
wa(2)=0.6750774
wa(3)=0.7252295
wa(4)=0.6495854
wb(1)=1.405189
wb(2)=1.454441
wb(3)=1.378797
wb(4)=1.458423

stime=2.*pi*(real(mtimer)/1440.+real(jstart))
sinlon=sin(2.*slon)
coslon=cos(2.*slon)
coslat=cos(slat)**2
sinlat=sin(2.*slat)

eta=0.
do i=1,4
  eta=eta+ba(i)*aa(i)*coslat*(cos(wa(i)*stime)*coslon-sin(wa(i)*stime)*sinlon)
  eta=eta+bb(i)*ab(i)*sinlat*(cos(wb(i)*stime)*coslon-sin(wb(i)*stime)*sinlon)
end do

eta=-eta

return
end subroutine mlotide

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine fills ocean data over land points
subroutine mlofill(x)

use cc_mpi
use indices_m
use soil_m

implicit none

include 'newmpar.h'
include 'mpif.h'

integer globalc,localc,ierr,iq,lnum
real miss,lsum
real, dimension(ifull), intent(inout) :: x
real, dimension(ifull) :: yy
real, dimension(ifull+iextra) :: xx
logical, dimension(ifull+iextra) :: smap

miss=999999.

xx=miss
where (.not.land)
  xx(1:ifull)=x
end where
globalc=1

! technically we only need 1 pass of this fill to ensure all water points have a non-trival neighbour
do while (globalc.gt.0)
  call bounds(xx)
  smap=abs(xx-miss).lt.0.1
  yy=xx(1:ifull)
  do iq=1,ifull
    if (smap(iq)) then
      lsum=0.
      lnum=0
      if (.not.smap(in(iq))) then
        lsum=lsum+xx(in(iq))
        lnum=lnum+1
      end if
      if (.not.smap(ie(iq))) then
        lsum=lsum+xx(ie(iq))
        lnum=lnum+1
      end if
      if (.not.smap(is(iq))) then
        lsum=lsum+xx(is(iq))
        lnum=lnum+1
      end if
      if (.not.smap(iw(iq))) then
        lsum=lsum+xx(iw(iq))
        lnum=lnum+1
      end if
      if (lnum.gt.0) then
        yy(iq)=lsum/real(lnum)
      end if
    end if
  end do
  xx(1:ifull)=yy
  smap(1:ifull)=abs(xx(1:ifull)-miss).lt.0.1
  localc=count(smap(1:ifull))
  call MPI_AllReduce(localc,globalc,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
end do

x=xx(1:ifull)

return
end subroutine mlofill

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! check for leap year

subroutine mloleap(tyear,ttest)

implicit none

integer, intent(in) :: tyear
logical, intent(out) :: ttest

ttest=.false.
if (mod(tyear,4).eq.0) ttest=.true.
if (mod(tyear,100).eq.0) ttest=.false.
if (mod(tyear,400).eq.0) ttest=.true.

return
end subroutine mloleap

end module mlodynamics
