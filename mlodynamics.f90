! These subroutines handle dynamics for the Mixed-Layer-Ocean model

! - Horizontal diffusion
! - River routing
! - Ocean dynamics
! - Ice dynamics

! This module also links to mlo.f90 which solves for 1D physics of
! ocean and ice processes

module mlodynamics

implicit none

private
public mlodiffusion,mlorouter,mlohadv,watbdy

real, dimension(:), allocatable, save :: watbdy
integer, parameter :: salfilt=0    ! Additional salinity filter (0=off, 1=Katzfey)
integer, parameter :: usetide=1    ! tidal forcing (0=Off, 1=On)
integer, parameter :: icemode=2    ! Ice stress (0=free-drift, 1=incompressible, 2=cavitating)
integer, parameter :: lmax   =1    ! 1=predictor-only, 2+=predictor-corrector
real, parameter :: k_smag=0.4      ! horizontal diffusion (0.4 in mom3, 2. in Griffies (2000))
real, parameter :: rhosn =330.     ! density snow (kg m^-3)
real, parameter :: rhoic =900.     ! density ice  (kg m^-3)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine processes horizontal diffusion, based on Griffies (2000)
! and McGregor's hordifg.f routines for CCAM.
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

integer k,i,iq,ida,idb
real hdif,c1,c2,c3,c4,cxx,dxx
real, dimension(ifull+iextra,wlev) :: u,v,dep
real, dimension(ifull+iextra) :: uc,vc,wc,ee,gg
real, dimension(ifull+iextra) :: t_kh,xfact,yfact
real, dimension(ifull,wlev) :: dudx,dvdx,dudy,dvdy
real, dimension(ifull) :: base
real, dimension(ifull) :: cc,ff,emi,ucc,vcc,wcc
real, dimension(wlev) :: cc1,cc2,cc3,cc4,d1,d2
logical, dimension(ifull+iextra) :: wtr

hdif=dt*(k_smag/pi)**2
ee=1.
where(land(1:ifull))
  ee(1:ifull)=0.
end where
call bounds(ee)
wtr=ee.gt.0.5

emi=1./em(1:ifull)**2
u=0.
v=0.
dep=0.
do k=1,wlev
  call mloexpdep(0,dep(1:ifull,k),k,0)
  call mloexport(2,u(1:ifull,k),k,0)
  call mloexport(3,v(1:ifull,k),k,0)
end do
call bounds(dep)
call boundsuv(u,v)
dep=max(dep,1.E-3)

! Interpolate velocities to horizontal surfaces
dudx=0.
dvdx=0.
dudy=0.
dvdy=0.
do iq=1,ifull
  if (wtr(iq)) then
    ida=1
    idb=1
    cc1=u(ieu(iq),:)
    cc2=u(iwu(iq),:)
    cc3=v(iev(iq),:)
    cc4=v(iwv(iq),:)
    d1=dep(ie(iq),:)
    d2=dep(iw(iq),:)
    cxx=0.
    do k=1,wlev
      dxx=dep(iq,k)
      call searchdelta(cc1,d1,dxx,cxx,ida,c1)
      call searchdelta(cc2,d2,dxx,cxx,idb,c2)
      call searchdelta(cc3,d1,dxx,cxx,ida,c3)
      call searchdelta(cc4,d2,dxx,cxx,idb,c4)
      dudx(iq,k)=(c1/em(ie(iq))-c2/em(iw(iq)))*0.5*em(iq)*em(iq)/ds
      dvdx(iq,k)=(c3/em(ie(iq))-c4/em(iw(iq)))*0.5*em(iq)*em(iq)/ds
    end do
    ida=1
    idb=1
    cc1=u(inu(iq),:)
    cc2=u(isu(iq),:)
    cc3=v(inv(iq),:)
    cc4=v(isv(iq),:)
    d1=dep(in(iq),:)
    d2=dep(is(iq),:)
    cxx=0.
    do k=1,wlev
      dxx=dep(iq,k)
      call searchdelta(cc1,d1,dxx,cxx,ida,c1)
      call searchdelta(cc2,d2,dxx,cxx,idb,c2)
      call searchdelta(cc3,d1,dxx,cxx,ida,c3)
      call searchdelta(cc4,d2,dxx,cxx,idb,c4)
      dudy(iq,k)=(c1/em(in(iq))-c2/em(is(iq)))*0.5*em(iq)*em(iq)/ds
      dvdy(iq,k)=(c3/em(in(iq))-c4/em(is(iq)))*0.5*em(iq)*em(iq)/ds
    end do
  end if
end do

do k=1,wlev
  uc(1:ifull) = ax(1:ifull)*u(1:ifull,k) + bx(1:ifull)*v(1:ifull,k)
  vc(1:ifull) = ay(1:ifull)*u(1:ifull,k) + by(1:ifull)*v(1:ifull,k)
  wc(1:ifull) = az(1:ifull)*u(1:ifull,k) + bz(1:ifull)*v(1:ifull,k)
  call bounds(uc)
  call bounds(vc)
  call bounds(wc)

  ! Smagorinsky
  cc=(dudx(:,k)-dvdy(:,k))**2+(dudy(:,k)+dvdx(:,k))**2
  t_kh(1:ifull)=sqrt(cc)*hdif*emi  ! this one with em in D terms
  call bounds(t_kh)
  xfact(1:ifull) = (t_kh(ie)+t_kh)*.5*ee(1:ifull)*ee(ie)
  yfact(1:ifull) = (t_kh(in)+t_kh)*.5*ee(1:ifull)*ee(in)
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
  u(1:ifull,k) = ax(1:ifull)*ucc + ay(1:ifull)*vcc + az(1:ifull)*wcc
  v(1:ifull,k) = bx(1:ifull)*ucc + by(1:ifull)*vcc + bz(1:ifull)*wcc
  
  call mloimport(2,u(1:ifull,k),k,0)
  call mloimport(3,v(1:ifull,k),k,0)
   
  do i=0,1
    gg=0.
    call mloexport(i,gg(1:ifull),k,0)
    call bounds(gg)
    ff = ( gg(1:ifull)*emi +             &
           xfact(1:ifull)*gg(ie) +       &
           xfact(iwu)*gg(iw) +           &
           yfact(1:ifull)*gg(in) +       &
           yfact(isv)*gg(is) ) / base
    call mloimport(i,ff,k,0)
  end do
  
  ! Jack Katzfey salinity filter
  if (salfilt.eq.1) then
    gg(1:ifull)=ff
    call bounds(gg)
    xfact(1:ifull)=ee(1:ifull)*ee(ie)
    yfact(1:ifull)=ee(1:ifull)*ee(in)
    call boundsuv(xfact,yfact)
    ff = ( gg(1:ifull)*emi +             &
           xfact(1:ifull)*gg(ie) +       &
           xfact(iwu)*gg(iw) +           &
           yfact(1:ifull)*gg(in) +       &
           yfact(isv)*gg(is) ) / base
    call mloimport(1,ff,k,0)
  end if
  
end do

return
end subroutine mlodiffusion

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates the river routing.
! This version currently assumes that the orography is sufficently resolved
! to determine the river flow.  Plan to read in effective gradients between
! grid boxes from high resolution river flow datasets.
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

! To speed up the code, we use a (semi-)implicit solution rather than an iterative approach
! This avoids additional MPI calls.  Also, mass is conserved since abs(flow).le.watbdy/10.
! for sensible choices of dt (e.g., dt=1200 for 60 km resolution).

! Note that unlike Miller et al (1994), this scheme currently has multiple outflows from a given
! grid box (i.e., following Gordon in Mk3.5).  So if the slope of the grid box is divergent, then
! on average water leaves the grid box, whereas if the slop of the grid box is convergent, then on
! average water enters the grid box.

! The scheme currently does not support outflow from lakes.

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
vel=0.35*sqrt(mslope/0.00005) ! from Miller et al (1994)
where (mslope.gt.1.E-10)
  vel=min(max(vel,0.15),5.)
elsewhere
  vel=1.E-10
end where
do i=1,4
  flow(:,i)=-watbdy(1:ifull)/(dp(:,i)/(vel(:,i)*dt)+1.) ! (kg/m^2)
end do
newwat=newwat+sum(flow,2)
  
! inflow
mslope=max(-slope,0.)
vel=0.35*sqrt(mslope/0.00005) ! from Miller et al (1994)
where (mslope.gt.1.E-10)
  vel=min(max(vel,0.15),5.)
elsewhere
  vel=1.E-10
end where
do i=1,4
  flow(:,i)=watbdy(xp(:,i))/(dp(:,i)/(vel(:,i)*dt)+1.) ! (kg/m^2)
  flow(:,i)=flow(:,i)*em(1:ifull)*em(1:ifull)/(em(xp(:,i))*em(xp(:,i))) ! correct for changing grid box area
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
! This subroutine implements some basic hydrostatic dynamics for the
! ocean and ice. (note that this subroutine is a work in progress)

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

integer iq,l,ll,ii,intsch,ierr,totits,itotits
integer jyear,jmonth,jday,jhour,jmin,mstart,mins
integer tyear,jstart
integer, dimension(ifull,wlev) :: nface
integer, dimension(12) :: ndoy
real alpha,maxloclseta,maxglobseta,maxloclip,maxglobip
real delpos,delneg,alph_p,dumpp,dumpn
real, dimension(ifull+iextra) :: ee,neta,dd,snu,snv,pice,ip
real, dimension(ifull+iextra) :: nfracice,ndic,ndsn,nsto,niu,niv,ndum
real, dimension(ifull+iextra) :: imass,spu,squ,sru,spv,sqv,srv,sou,sov
real, dimension(ifull+iextra) :: ibu,ibv,icu,icv
real, dimension(ifull) :: i_u,i_v,i_sto,rhobaru,rhobarv
real, dimension(ifull) :: div,seta,w_e,imu,imv
real, dimension(ifull) :: tnu,tsu,tev,twv,rhou,rhov
real, dimension(ifull) :: dpsdxu,dpsdyu,dpsdxv,dpsdyv
real, dimension(ifull) :: dttdxu,dttdyu,dttdxv,dttdyv
real, dimension(ifull) :: detadxu,detadyu,detadxv,detadyv
real, dimension(ifull) :: dipdxu,dipdyu,dipdxv,dipdyv
real, dimension(ifull) :: au,bu,cu,av,bv,cv,odum,oeu,oev
real, dimension(ifull) :: nip,ipn,ipe,ips,ipw,ipmax
real, dimension(ifull) :: dsnudeta,dsnuwdeta,dsnvdeta,dsnvsdeta,ddivdeta
real, dimension(ifull) :: sssa,sssb,sssc,sssd
real, dimension(ifull) :: depu,depv,ddu,ddv,ddddxu,ddddxv,ddddyu,ddddyv
real, dimension(ifull+iextra,4) :: nit
real, dimension(ifull+iextra,wlev) :: nu,nv,nt,ns
real, dimension(ifull+iextra,wlev) :: cou,cov,cow
real, dimension(ifull+iextra,wlev) :: dep,rhobar,rho,dz
real, dimension(ifull,1) :: siu,siv
real, dimension(ifull,4) :: i_it
real, dimension(ifull,wlev) :: w_u,w_v,w_t,w_s,nw,dzbar,dum,dumdep,dumdz
real, dimension(ifull,wlev) :: nuh,nvh,xg,yg,uau,uav,stagu,stagv
real, dimension(ifull,wlev) :: kku,llu,mmu,nnu
real, dimension(ifull,wlev) :: kkv,llv,mmv,nnv
real, dimension(ifull,wlev) :: drhobardxu,drhobardyu,drhobardxv,drhobardyv
real, dimension(ifull,2) :: stwgt
real, dimension(:,:), allocatable, save :: oldu1,oldu2,oldv1,oldv2
real*8, dimension(ifull,wlev) :: x3d,y3d,z3d
logical, dimension(ifull+iextra) :: wtr
logical, dimension(ifull) :: stest
logical lleap
integer, parameter :: llmax=1000   ! iterations for calculating surface height
real, parameter :: tol  = 5.E-4    ! Tolerance for GS solver (water)
real, parameter :: itol = 5.E-1    ! Tolerance for GS solver (ice)
real, parameter :: sal = 0.948     ! SAL parameter for tidal forcing

data ndoy/0,31,59,90,120,151,181,212,243,273,304,334/

! new z levels for including free surface eta (effectively sigma-height levels)
! newz=-eta+oldz*(1+eta/maxdepth)
! newdz=olddz*(1+eta/maxdepth)
! where 0<=oldz<=maxdepth and -eta<=newz<=maxdepth
! depth below free suface is then newz+eta=oldz*(1+eta/maxdepth)
! horizontal derivatives are along constant oldz surfaces

! initial values for working arrays
cou=0.
cov=0.
cow=0.
rho=1030.
rhobar=1030.
pice=0.
imass=0.
dep=0.
dz=0.
w_t=293.
w_s=0.
w_u=0.
w_v=0.
w_e=0.
i_it=273.
i_sto=0.
i_u=0.
i_v=0.
nt=293.
ns=30.
nu=0.
nv=0.
neta=0.
nit=273.
nfracice=0.
ndic=0.
ndsn=0.
nsto=0.
niu=0.
niv=0.

! alternate semi-lagrangian interpolation
intsch=mod(ktau,2)

!Define land/sea mask
ee=1.
where(land(1:ifull))
  ee(1:ifull)=0.
end where
call bounds(ee,nrows=2)
wtr=ee.gt.0.5

! PRECOMPUTE WEIGHTS FOR CALCULATING STAGGERED GRADIENTS
stwgt=0.
where (wtr(in).and.wtr(ine).and.wtr(is).and.wtr(ise))
  stwgt(:,1)=1.
end where
where (wtr(ie).and.wtr(ien).and.wtr(iw).and.wtr(iwn))
  stwgt(:,2)=1.
end where

! IMPORT WATER AND ICE DATA -----------------------------------------
do ii=1,wlev
  call mloexpdep(0,dep(1:ifull,ii),ii,0)
  call mloexpdep(1,dz(:,ii),ii,0)
  call mloexport(0,w_t(:,ii),ii,0)
  call mloexport(1,w_s(:,ii),ii,0)
  call mloexport(2,w_u(:,ii),ii,0)
  call mloexport(3,w_v(:,ii),ii,0)
end do
call mloexport(4,w_e,0,0)
do ii=1,4
  call mloexpice(i_it(:,ii),ii,0)
end do
call mloexpice(fracice,5,0)
call mloexpice(sicedep,6,0)
call mloexpice(snowd,7,0)
call mloexpice(i_sto,8,0)
call mloexpice(i_u,9,0)
call mloexpice(i_v,10,0)
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
call bounds(dd,nrows=2)

! calculate bathymetry gradients
tnu=0.5*(dd(in)+dd(ine))
tsu=0.5*(dd(is)+dd(ise))
tev=0.5*(dd(ie)+dd(ien))
twv=0.5*(dd(iw)+dd(iwn))
ddddxu=(dd(ie)-dd(1:ifull))*emu(1:ifull)/ds
ddddyu=stwgt(:,1)*0.5*(tnu-tsu)*emu(1:ifull)/ds
ddddxv=stwgt(:,2)*0.5*(tev-twv)*emv(1:ifull)/ds
ddddyv=(dd(in)-dd(1:ifull))*emv(1:ifull)/ds

! estimate tidal forcing
if (usetide.eq.1) then
  jyear=kdate/10000
  jstart=0
  if (jyear.gt.1900) then
    do tyear=1900,jyear-1
      call mloleap(tyear,lleap)
      if (lleap) jstart=jstart+1
      jstart=jstart+365
    end do
  else if (jyear.lt.1900) then
    do tyear=1899,jyear,-1
      call mloleap(tyear,lleap)
      if (lleap) jstart=jstart-1
      jstart=jstart-365
    end do
  end if
  mstart=720 ! base time is 12Z 31 Dec 1899
  jmonth=(kdate-jyear*10000)/100
  jday=kdate-jyear*10000-jmonth*100
  jhour=ktime/100
  jmin=ktime-jhour*100
  mstart=mstart+1440*(ndoy(jmonth)+jday-1)+60*jhour+jmin ! mins from start of year
  if (jmonth.gt.2) then
    call mloleap(jyear,lleap)
    if (lleap) mstart=mstart+1440
  end if
  mins = mtimer + mstart
  call mlotide(ndum(1:ifull),rlongg,rlatt,mins,jstart)
  ndum(1:ifull)=ndum(1:ifull)*ee(1:ifull)
  call bounds(ndum,nrows=2)
  
  tnu=0.5*(ndum(in)+ndum(ine))
  tsu=0.5*(ndum(is)+ndum(ise))
  tev=0.5*(ndum(ie)+ndum(ien))
  twv=0.5*(ndum(iw)+ndum(iwn))
  dttdxu=(ndum(ie)-ndum(1:ifull))*emu(1:ifull)/ds
  dttdyu=stwgt(:,1)*0.5*(tnu-tsu)*emu(1:ifull)/ds
  dttdxv=stwgt(:,2)*0.5*(tev-twv)*emv(1:ifull)/ds
  dttdyv=(ndum(in)-ndum(1:ifull))*emv(1:ifull)/ds  
else
  dttdxu=0.
  dttdyu=0.
  dttdxv=0.
  dttdyv=0.
end if

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
nu(1:ifull,:)=w_u
nv(1:ifull,:)=w_v
nt(1:ifull,:)=w_t
ns(1:ifull,:)=w_s
where (wtr(1:ifull))
  nfracice(1:ifull)=fracice
  ndic(1:ifull)=sicedep
  ndsn(1:ifull)=snowd*0.001
end where
nit(1:ifull,:)=i_it
nsto(1:ifull)=i_sto
niu(1:ifull)=i_u
niv(1:ifull)=i_v

totits=0
itotits=0
do l=1,lmax ! predictor-corrector loop

  ! surface pressure gradients (including ice)
  imass(1:ifull)=ndic(1:ifull)*rhoic+ndsn(1:ifull)*rhosn ! ice mass per unit area (kg/m^2)
  pice(1:ifull)=ps(1:ifull)+grav*nfracice(1:ifull)*imass(1:ifull) ! pressure due to atmosphere and ice at top of water column
  call bounds(pice,nrows=2)
  tnu=0.5*(pice(in)+pice(ine))
  tsu=0.5*(pice(is)+pice(ise))
  tev=0.5*(pice(ie)+pice(ien))
  twv=0.5*(pice(iw)+pice(iwn))
  dpsdxu=(pice(ie)-pice(1:ifull))*emu(1:ifull)/ds
  dpsdyu=stwgt(:,1)*0.5*(tnu-tsu)*emu(1:ifull)/ds
  dpsdxv=stwgt(:,2)*0.5*(tev-twv)*emv(1:ifull)/ds
  dpsdyv=(pice(in)-pice(1:ifull))*emv(1:ifull)/ds


  ! calculate normalised density rhobar (unstaggered)
  odum=(1.+neta(1:ifull)/dd(1:ifull))
  do ii=1,wlev
    dumdep(:,ii)=dep(1:ifull,ii)*odum
    dumdz(:,ii)=dz(1:ifull,ii)*odum
  end do
  call mloexpdensity(rho(1:ifull,:),nt(1:ifull,:),ns(1:ifull,:),dumdep(1:ifull,:),dumdz(1:ifull,:),pice(1:ifull),0)
  call bounds(rho)
  rhobar(1:ifull,1)=rho(1:ifull,1)*dz(1:ifull,1)
  do ii=2,wlev
    rhobar(1:ifull,ii)=rhobar(1:ifull,ii-1)+rho(1:ifull,ii)*dz(1:ifull,ii)
  end do
  do ii=1,wlev
    rhobar(1:ifull,ii)=rhobar(1:ifull,ii)/dzbar(:,ii)
  end do
  call bounds(rhobar,nrows=2)

  ! ADVECT WATER ----------------------------------------------------
  ! Water currents are advected using sem-Lagrangian advection
  ! (i.e., taken from McGregor's CCAM advection routines).
  ! Velocity is set to zero at ocean boundaries, and changing
  ! bathymetry is accounted for in the integral of the divergence.

  ! nu,nv,nt,ns are now reset to begin the advection from the current time step
  nu(1:ifull,:)=w_u
  nv(1:ifull,:)=w_v
  ns(1:ifull,:)=w_s
  nt(1:ifull,:)=w_t

  ! estimate vertical velocity
  call mlostaguv(w_u,w_v,cou(1:ifull,:),cov(1:ifull,:))
  do ii=1,wlev
    cou(1:ifull,ii)=cou(1:ifull,ii)*ee(1:ifull)*ee(ie)
    cov(1:ifull,ii)=cov(1:ifull,ii)*ee(1:ifull)*ee(in)
  end do
  call boundsuv(cou,cov)
  call getww(cou,cov,dd(1:ifull),dz(1:ifull,:),ee(1:ifull),nw)

  ! Vertical advection (first call)
  call mlovadv(0.5*dt,nw,nu(1:ifull,:),nv(1:ifull,:),ns(1:ifull,:),nt(1:ifull,:),dep(1:ifull,:), &
               dz(1:ifull,:),wtr(1:ifull))

  ! estimate currents at t+1/2 for semi-Lagrangian advection
  if (l.eq.1) then
    nuh=(15.*w_u-10.*oldu1+3.*oldu2)/8. ! U at t+1/2
    nvh=(15.*w_v-10.*oldv1+3.*oldv2)/8. ! V at t+1/2
  else
    nuh=0.5*nu(1:ifull,:)+0.5*w_u ! U at t+1/2
    nvh=0.5*nv(1:ifull,:)+0.5*w_v ! V at t+1/2
  end if

  ! Calculate depature points along oldz levels
  call mlodeps(nuh,nvh,nface,xg,yg,x3d,y3d,z3d,dep,wtr)

  ! Convert (u,v) to cartesian coordinates (U,V,W)
  do ii=1,wlev
    cou(1:ifull,ii)=ax(1:ifull)*nu(1:ifull,ii)+bx(1:ifull)*nv(1:ifull,ii)
    cov(1:ifull,ii)=ay(1:ifull)*nu(1:ifull,ii)+by(1:ifull)*nv(1:ifull,ii)
    cow(1:ifull,ii)=az(1:ifull)*nu(1:ifull,ii)+bz(1:ifull)*nv(1:ifull,ii)
  end do

  ! Horizontal advection for U,V,W
  call mloints(cou,dep,intsch,nface,xg,yg,2,wtr)
  call mloints(cov,dep,intsch,nface,xg,yg,2,wtr)
  call mloints(cow,dep,intsch,nface,xg,yg,2,wtr)
 
  ! Rotate vector to arrival point
  call mlorot(cou(1:ifull,:),cov(1:ifull,:),cow(1:ifull,:),x3d,y3d,z3d)

  ! Convert (U,V,W) back to conformal cubic coordinates
  do ii=1,wlev
    nu(1:ifull,ii)=ax(1:ifull)*cou(1:ifull,ii)+ay(1:ifull)*cov(1:ifull,ii)+az(1:ifull)*cow(1:ifull,ii)
    nv(1:ifull,ii)=bx(1:ifull)*cou(1:ifull,ii)+by(1:ifull)*cov(1:ifull,ii)+bz(1:ifull)*cow(1:ifull,ii)
    nu(1:ifull,ii)=nu(1:ifull,ii)*ee(1:ifull)
    nv(1:ifull,ii)=nv(1:ifull,ii)*ee(1:ifull)
  end do

  ! Horizontal advection for T,S
  call mlointsb(nt,dep,intsch,nface,xg,yg,2,wtr)
  call mlointsb(ns,dep,intsch,nface,xg,yg,5,wtr)
  ns=max(ns,0.)

  ! FREE SURFACE CALCULATION ----------------------------------------

  ! stagger the following values for t=tstar
  ! This is the same as eps=0.5 in JLM's atmosphere semi-Lagrangian dynamics
  do ii=1,wlev
    uau(:,ii)=nu(1:ifull,ii)*(1.-0.25*dt*dt*f(1:ifull)*f(1:ifull))+dt*f(1:ifull)*nv(1:ifull,ii)
    uav(:,ii)=nv(1:ifull,ii)*(1.-0.25*dt*dt*f(1:ifull)*f(1:ifull))-dt*f(1:ifull)*nu(1:ifull,ii)
  end do
  call mlostaguv(uau,uav,cou(1:ifull,:),cov(1:ifull,:))
  stagu=cou(1:ifull,:)
  stagv=cov(1:ifull,:)

  ! interpolate density gradient to horizontal surfaces
  call stagtruedelta(rhobar,dep,wtr,drhobardxu,drhobardyu,drhobardxv,drhobardyv)
  do ii=1,wlev
    drhobardyu(:,ii)=drhobardyu(:,ii)*stwgt(:,1)
    drhobardxv(:,ii)=drhobardxv(:,ii)*stwgt(:,2)
  end do

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
  ndum(1:ifull)=w_e
  call bounds(ndum,nrows=2)
  tnu=0.5*(ndum(in)+ndum(ine))
  tsu=0.5*(ndum(is)+ndum(ise))
  tev=0.5*(ndum(ie)+ndum(ien))
  twv=0.5*(ndum(iw)+ndum(iwn))
  detadxu=(ndum(ie)-ndum(1:ifull))*emu(1:ifull)/ds
  detadyu=stwgt(:,1)*0.5*(tnu-tsu)*emu(1:ifull)/ds
  detadxv=stwgt(:,2)*0.5*(tev-twv)*emv(1:ifull)/ds
  detadyv=(ndum(in)-ndum(1:ifull))*emv(1:ifull)/ds
  oeu=0.5*(ndum(1:ifull)+ndum(ie))
  oev=0.5*(ndum(1:ifull)+ndum(in))
  ddu=0.5*(dd(1:ifull)+dd(ie))
  ddv=0.5*(dd(1:ifull)+dd(in))
  do ii=1,wlev
    depu=0.5*(dep(1:ifull,ii)+dep(ie,ii))
    depv=0.5*(dep(1:ifull,ii)+dep(in,ii))
    rhou=0.5*(rho(1:ifull,ii)+rho(ie,ii))
    rhov=0.5*(rho(1:ifull,ii)+rho(in,ii))
    rhobaru=(rhobar(1:ifull,ii)+rhobar(ie,ii))*0.5
    rhobarv=(rhobar(1:ifull,ii)+rhobar(in,ii))*0.5

    ! u^(t+1) = nu = au^(t*) + bu*dpdx + cu*dpdy (staggered)
    ! v^(t+1) = nv = av^(t*) + bv*dpdy + cv*dpdx (staggered)

    au=cou(1:ifull,ii)/(1.+0.25*dt*dt*fu(1:ifull)*fu(1:ifull))
    bu=-dt/(rhou*(1.+0.25*dt*dt*fu(1:ifull)*fu(1:ifull)))
    cu=0.5*dt*fu(1:ifull)*bu

    av=cov(1:ifull,ii)/(1.+0.25*dt*dt*fv(1:ifull)*fv(1:ifull))
    bv=-dt/(rhov*(1.+0.25*dt*dt*fv(1:ifull)*fv(1:ifull)))
    cv=-0.5*dt*fv(1:ifull)*bv

    !dppdxu=dpsdxu+grav*depu*(1+etau/ddu)*drhobardxu+grav*rhobaru*(depu/ddu)*(detadxu-(eta/ddu)*dddudxu)
    !dppdyu=dpsdyu+grav*depu*(1+etau/ddu)*drhobardyu+grav*rhobaru*(depu/ddu)*(detadyu-(eta/ddu)*dddudyu)
    !dppdxv=dpsdxv+grav*depv*(1+etav/ddv)*drhobardxv+grav*rhobarv*(depv/ddv)*(detadxv-(eta/ddv)*dddvdxv)
    !dppdyv=dpsdyv+grav*depv*(1+etav/ddu)*drhobardyv+grav*rhobarv*(depv/ddv)*(detadyv-(eta/ddv)*dddvdyv)
    
    !nu=kku+llu*etau+mmu*detadxu+nnu*detadyu (staggered)
    !nv=kkv+llv*etav+mmv*detadyv+nnv*detadxv (staggered)
    
    !int nu dz = sou+spu*etau+squ*detadxu+sru*detadyu
    !int nv dz = sov+spv*etav+sqv*detadyv+srv*detadxv

    if (usetide.eq.1) then
      llu(:,ii)=0.5*(bu*(drhobardxu(:,ii)-rhobaru*ddddxu/ddu)+cu*(drhobardyu(:,ii)-rhobaru*ddddyu/ddu)) &
                *grav*depu/ddu
      mmu(:,ii)=0.5*bu*grav*(rhobaru*depu/ddu+rhou*(1.-sal))
      nnu(:,ii)=0.5*cu*grav*(rhobaru*depu/ddu+rhou*(1.-sal))
      kku(:,ii)=au+bu*(dpsdxu+grav*rhou*dttdxu+grav*drhobardxu(:,ii)*depu) &
                  +cu*(dpsdyu+grav*rhou*dttdyu+grav*drhobardyu(:,ii)*depu) &
                  +llu(:,ii)*oeu+mmu(:,ii)*detadxu+nnu(:,ii)*detadyu

      llv(:,ii)=0.5*(bv*(drhobardyv(:,ii)-rhobarv*ddddyv/ddv)+cv*(drhobardxv(:,ii)-rhobarv*ddddxv/ddv)) &
                *grav*depv/ddv
      mmv(:,ii)=0.5*bv*grav*(rhobarv*depv/ddv+rhov*(1.-sal))
      nnv(:,ii)=0.5*cv*grav*(rhobarv*depv/ddv+rhov*(1.-sal))
      kkv(:,ii)=av+bv*(dpsdyv+grav*rhov*dttdyv+grav*drhobardyv(:,ii)*depv) &
                  +cv*(dpsdxv+grav*rhov*dttdxv+grav*drhobardxv(:,ii)*depv) &
                  +llv(:,ii)*oev+mmv(:,ii)*detadyv+nnv(:,ii)*detadxv
    else
      llu(:,ii)=0.5*(bu*(drhobardxu(:,ii)-rhobaru*ddddxu/ddu)+cu*(drhobardyu(:,ii)-rhobaru*ddddyu/ddu)) &
                *grav*depu/ddu
      mmu(:,ii)=0.5*bu*grav*rhobaru*depu/ddu
      nnu(:,ii)=0.5*cu*grav*rhobaru*depu/ddu
      kku(:,ii)=au+bu*(dpsdxu+grav*drhobardxu(:,ii)*depu)+cu*(dpsdyu+grav*drhobardyu(:,ii)*depu) &
                  +llu(:,ii)*oeu+mmu(:,ii)*detadxu+nnu(:,ii)*detadyu

      llv(:,ii)=0.5*(bv*(drhobardyv(:,ii)-rhobarv*ddddyv/ddv)+cv*(drhobardxv(:,ii)-rhobarv*ddddxv/ddv)) &
                *grav*depv/ddv
      mmv(:,ii)=0.5*bv*grav*rhobarv*depv/ddv
      nnv(:,ii)=0.5*cv*grav*rhobarv*depv/ddv
      kkv(:,ii)=av+bv*(dpsdyv+grav*drhobardyv(:,ii)*depv)+cv*(dpsdxv+grav*drhobardxv(:,ii)*depv) &
                  +llv(:,ii)*oev+mmv(:,ii)*detadyv+nnv(:,ii)*detadxv
    end if

    sou(1:ifull)=sou(1:ifull)+kku(:,ii)*0.5*(dz(1:ifull,ii)+dz(ie,ii))
    spu(1:ifull)=spu(1:ifull)+llu(:,ii)*0.5*(dz(1:ifull,ii)+dz(ie,ii))
    squ(1:ifull)=squ(1:ifull)+mmu(:,ii)*0.5*(dz(1:ifull,ii)+dz(ie,ii))
    sru(1:ifull)=sru(1:ifull)+nnu(:,ii)*0.5*(dz(1:ifull,ii)+dz(ie,ii))

    sov(1:ifull)=sov(1:ifull)+kkv(:,ii)*0.5*(dz(1:ifull,ii)+dz(in,ii))
    spv(1:ifull)=spv(1:ifull)+llv(:,ii)*0.5*(dz(1:ifull,ii)+dz(in,ii))
    sqv(1:ifull)=sqv(1:ifull)+mmv(:,ii)*0.5*(dz(1:ifull,ii)+dz(in,ii))
    srv(1:ifull)=srv(1:ifull)+nnv(:,ii)*0.5*(dz(1:ifull,ii)+dz(in,ii))

  end do

  ! update land boundaries
  sou(1:ifull)=sou(1:ifull)*ee(1:ifull)*ee(ie)
  spu(1:ifull)=spu(1:ifull)*ee(1:ifull)*ee(ie)
  squ(1:ifull)=squ(1:ifull)*ee(1:ifull)*ee(ie)
  sru(1:ifull)=sru(1:ifull)*ee(1:ifull)*ee(ie)
  sov(1:ifull)=sov(1:ifull)*ee(1:ifull)*ee(in)
  spv(1:ifull)=spv(1:ifull)*ee(1:ifull)*ee(in)
  sqv(1:ifull)=sqv(1:ifull)*ee(1:ifull)*ee(in)
  srv(1:ifull)=srv(1:ifull)*ee(1:ifull)*ee(in)

  call boundsuv(spu,spv)
  call boundsuv(squ,sqv)

  ! prep gradient terms to improve numerical stability
  sssa=spu(1:ifull)*0.5-squ(1:ifull)*emu(1:ifull)/ds
  sssb=spu(iwu)*0.5+squ(iwu)*emu(iwu)/ds
  sssc=spv(1:ifull)*0.5-sqv(1:ifull)*emv(1:ifull)/ds
  sssd=spv(isv)*0.5+sqv(isv)*emv(isv)/ds

  ! Iteratively solve for free surface height, eta
  alpha=0.1
  do ll=1,llmax

    ! 9-point version -----------------------------------------------
    ! would be nice to simplify this down to a 5-point stencil, rather than the effective 9-point used here
    ! However, when this is done (e.g., in unstaggered coordinates), then the non-linear terms arising from
    ! the integration of the column (i.e., the column height is a function of eta), results in a decoupling
    ! of the solution between adjacent grid points.  This 9-point version avoids this problem, but the two
    ! bounds calls and nrows=2 slows things down a little.

    ! calculate neta gradients
    call bounds(neta,nrows=2)
    tnu=0.5*(neta(in)+neta(ine))
    tsu=0.5*(neta(is)+neta(ise))
    tev=0.5*(neta(ie)+neta(ien))
    twv=0.5*(neta(iw)+neta(iwn))
    detadxu=(neta(ie)-neta(1:ifull))*emu(1:ifull)/ds
    detadyu=stwgt(:,1)*0.5*(tnu-tsu)*emu(1:ifull)/ds
    detadxv=stwgt(:,2)*0.5*(tev-twv)*emv(1:ifull)/ds
    detadyv=(neta(in)-neta(1:ifull))*emv(1:ifull)/ds
    ! column integrated velocity terms
    snu(1:ifull)=sou(1:ifull)+spu(1:ifull)*0.5*(neta(ie)+neta(1:ifull))+squ(1:ifull)*detadxu+sru(1:ifull)*detadyu
    snv(1:ifull)=sov(1:ifull)+spv(1:ifull)*0.5*(neta(in)+neta(1:ifull))+sqv(1:ifull)*detadyv+srv(1:ifull)*detadxv
    snu(1:ifull)=snu(1:ifull)*(1.+0.5*(neta(ie)+neta(1:ifull)+ndum(ie)+ndum(1:ifull))/(dd(ie)+dd(1:ifull))) ! conserve volume
    snv(1:ifull)=snv(1:ifull)*(1.+0.5*(neta(in)+neta(1:ifull)+ndum(in)+ndum(1:ifull))/(dd(in)+dd(1:ifull))) ! conserve volume
    call boundsuv(snu,snv)
    div=(snu(1:ifull)/emu(1:ifull)-snu(iwu)/emu(iwu)+snv(1:ifull)/emv(1:ifull)-snv(isv)/emv(isv)) &
        *em(1:ifull)*em(1:ifull)/ds
    ! Update neta using d(div)/d(neta) to improve stability
    dsnudeta=sssa*(1.+0.5*(neta(ie)+neta(1:ifull)+ndum(ie)+ndum(1:ifull))/(dd(ie)+dd(1:ifull))) &
             +0.5*snu(1:ifull)/(dd(ie)+dd(1:ifull))
    dsnuwdeta=sssb*(1.+0.5*(neta(iw)+neta(1:ifull)+ndum(iw)+ndum(1:ifull))/(dd(iw)+dd(1:ifull))) &
             +0.5*snu(iwu)/(dd(iw)+dd(1:ifull))
    dsnvdeta=sssc*(1.+0.5*(neta(in)+neta(1:ifull)+ndum(in)+ndum(1:ifull))/(dd(in)+dd(1:ifull))) &
             +0.5*snv(1:ifull)/(dd(in)+dd(1:ifull))
    dsnvsdeta=sssd*(1.+0.5*(neta(is)+neta(1:ifull)+ndum(is)+ndum(1:ifull))/(dd(is)+dd(1:ifull))) &
             +0.5*snv(isv)/(dd(is)+dd(1:ifull))
    ddivdeta=(dsnudeta/emu(1:ifull)-dsnuwdeta/emu(iwu)+dsnvdeta/emv(1:ifull)-dsnvsdeta/emv(isv))*em(1:ifull)*em(1:ifull)/ds
    seta=-neta(1:ifull)+w_e-div/(1./dt+ddivdeta)
    
    ! The following expression limits the minimum depth to 1m
    seta=max(seta,1.-dd(1:ifull)-neta(1:ifull))*ee(1:ifull) ! this should become a land point
    neta(1:ifull)=alpha*seta+neta(1:ifull)
    
    ! Break iterative loop when maximum error is below tol (expensive)
    maxloclseta=maxval(abs(seta))
    call MPI_AllReduce(maxloclseta,maxglobseta,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
    if (maxglobseta.lt.tol.and.ll.gt.2) exit

    totits=totits+1
  end do

  ! volume conservation for water ---------------------------------------
  ! (include mass conservation in mlo.f90 due to thermal change)
  odum=neta(1:ifull)-w_e
  call ccglobal_posneg(odum,delpos,delneg)
  alph_p = -delneg/max(delpos,1.e-30)
  alph_p = min(alph_p,sqrt(alph_p))  ! best option
  alph_p = max(alph_p,1.E-20)
  neta(1:ifull)=w_e+alph_p*max(0.,odum)+min(0.,odum)/alph_p

  call bounds(neta,nrows=2)
  tnu=0.5*(neta(in)+neta(ine))
  tsu=0.5*(neta(is)+neta(ise))
  tev=0.5*(neta(ie)+neta(ien))
  twv=0.5*(neta(iw)+neta(iwn))
  detadxu=(neta(ie)-neta(1:ifull))*emu(1:ifull)/ds
  detadyu=stwgt(:,1)*0.5*(tnu-tsu)*emu(1:ifull)/ds
  detadxv=stwgt(:,2)*0.5*(tev-twv)*emv(1:ifull)/ds
  detadyv=(neta(in)-neta(1:ifull))*emv(1:ifull)/ds

  ! Update currents once neta is calculated
  do ii=1,wlev
    ! update currents (staggered)
    nu(1:ifull,ii)=kku(:,ii)+llu(:,ii)*0.5*(neta(ie)+neta(1:ifull))+mmu(:,ii)*detadxu+nnu(:,ii)*detadyu
    nv(1:ifull,ii)=kkv(:,ii)+llv(:,ii)*0.5*(neta(in)+neta(1:ifull))+mmv(:,ii)*detadyv+nnv(:,ii)*detadxv
    ! update land boundaries
    nu(1:ifull,ii)=nu(1:ifull,ii)*ee(1:ifull)*ee(ie)
    nv(1:ifull,ii)=nv(1:ifull,ii)*ee(1:ifull)*ee(in)
  end do
  call boundsuv(nu,nv)
  ! update vertical velocity
  call getww(nu,nv,dd(1:ifull),dz(1:ifull,:),ee(1:ifull),nw)
  ! unstagger nu and nv\
  do ii=1,wlev
    where(.not.wtr(1:ifull))
      nu(1:ifull,ii)=stagu(:,ii) ! to help reversible staggering at land boundaries
      nv(1:ifull,ii)=stagv(:,ii)
    end where
  end do
  call mlounstaguv(nu(1:ifull,:),nv(1:ifull,:),nuh,nvh)
  nu(1:ifull,:)=nuh
  nv(1:ifull,:)=nvh

  ! Vertical advection (second call)
  call mlovadv(0.5*dt,nw,nu(1:ifull,:),nv(1:ifull,:),ns(1:ifull,:),nt(1:ifull,:),dep(1:ifull,:), &
               dz(1:ifull,:),wtr(1:ifull))

  ! UPDATE ICE DYNAMICS ---------------------------------------------
  ! Here we start by calculating the ice velocity and then advecting
  ! the various ice prognostic variables.

  ! convert to staggered coordinates
  uau(:,1)=i_u+dt*f(1:ifull)*i_v
  uav(:,1)=i_v-dt*f(1:ifull)*i_u
  call mlostaguv(uau(:,1:1),uav(:,1:1),siu,siv)

  ! Update ice velocities
  ! niu and niv hold the free drift solution (staggered).  Wind stress terms are updated in mlo.f90
  niu(1:ifull)=(siu(:,1)-dt*grav*(detadxu+dt*fu(1:ifull)*detadyu))/(1.+dt*dt*fu(1:ifull)*fu(1:ifull))
  niv(1:ifull)=(siv(:,1)-dt*grav*(detadyv-dt*fv(1:ifull)*detadxv))/(1.+dt*dt*fv(1:ifull)*fv(1:ifull))
  niu(1:ifull)=niu(1:ifull)*ee(1:ifull)*ee(ie)
  niv(1:ifull)=niv(1:ifull)*ee(1:ifull)*ee(in)
  call boundsuv(niu,niv)
  
  imass(1:ifull)=max(imass(1:ifull),10.)
  call bounds(imass)
  if (icemode.gt.0) then
    ! (staggered)
    ! niu(t+1) = niu + ibu*dipdx + icu*dipdy
    ! niv(t+1) = niv + ibv*dipdy + icv*dipdx
  
    imu=0.5*(imass(1:ifull)+imass(ie))
    imv=0.5*(imass(1:ifull)+imass(in))
 
    ibu(1:ifull)=-dt/(imu*(1.+dt*dt*fu(1:ifull)*fu(1:ifull)))*ee(1:ifull)*ee(ie)
    ibv(1:ifull)=-dt/(imv*(1.+dt*dt*fv(1:ifull)*fv(1:ifull)))*ee(1:ifull)*ee(in)
    icu(1:ifull)=dt*fu(1:ifull)*ibu(1:ifull)*ee(1:ifull)*ee(ie)
    icv(1:ifull)=-dt*fv(1:ifull)*ibv(1:ifull)*ee(1:ifull)*ee(in)
    call boundsuv(ibu,ibv)

    ! Iterative loop to estimate ice 'pressure'
    alpha=0.9
    ip=0. ! free drift solution
    do ll=1,llmax
  
      call bounds(ip,nrows=2)
      spu(1:ifull)=0.5*(ip(in)+ip(ine))
      squ(1:ifull)=0.5*(ip(is)+ip(ise))
      spv(1:ifull)=0.5*(ip(ie)+ip(ien))
      sqv(1:ifull)=0.5*(ip(iw)+ip(iwn))
      spu(1:ifull)=stwgt(:,1)*icu(1:ifull)*(spu(1:ifull)-squ(1:ifull))
      spv(1:ifull)=stwgt(:,2)*icv(1:ifull)*(spv(1:ifull)-sqv(1:ifull))
      call boundsuv(spu,spv)

      !dipdxu=(ip(ie)-ip(1:ifull))*emu(1:ifull)/ds
      !dipdyu=0.5*spu(1:ifull)*emu(1:ifull)/(ds*icu(1:ifull))
      !dipdxv=0.5*spv(1:ifull)*emv(1:ifull)/(ds*icv(1:ifull))
      !dipdyv=(ip(in)-ip(1:ifull))*emv(1:ifull)/ds

      !snu(1:ifull)=niu(1:ifull)+ibu*dipdxu+icu*dipdyu
      !snv(1:ifull)=niv(1:ifull)+ibv*dipdyv+icv*dipdxv
      !call boundsuv(snu,snv)
      !div=(snu(1:ifull)/emu(1:ifull)-snu(iwu)/emu(iwu)+snv(1:ifull)/emv(1:ifull)-snv(isv)/emv(isv)) &
      !    *em(1:ifull)*em(1:ifull)/ds
      div=(niu(1:ifull)/emu(1:ifull)-niu(iwu)/emu(iwu))*em(1:ifull)*em(1:ifull)/ds                          &
         +(ibu(1:ifull)*(ip(ie)-ip(1:ifull))-ibu(iwu)*(ip(1:ifull)-ip(iw)))*em(1:ifull)*em(1:ifull)/(ds*ds) &
         +(spu(1:ifull)-spu(iwu))*0.5*em(1:ifull)*em(1:ifull)/(ds*ds)                                       &
         +(niv(1:ifull)/emv(1:ifull)-niv(isv)/emv(isv))*em(1:ifull)*em(1:ifull)/ds                          &          
         +(ibv(1:ifull)*(ip(in)-ip(1:ifull))-ibv(isv)*(ip(1:ifull)-ip(is)))*em(1:ifull)*em(1:ifull)/(ds*ds) &
         +(spv(1:ifull)-spv(isv))*0.5*em(1:ifull)*em(1:ifull)/(ds*ds)

      nip=ip(1:ifull)
      where (div.lt.0..and.wtr(1:ifull))
        nip=((niu(1:ifull)/emu(1:ifull)-niu(iwu)/emu(iwu))*ds    &
            +(ibu(1:ifull)*ip(ie)+ibu(iwu)*ip(iw))               &
            +(spu(1:ifull)-spu(iwu))*0.5                         &
            +(niv(1:ifull)/emv(1:ifull)-niv(isv)/emv(isv))*ds    &          
            +(ibv(1:ifull)*ip(in)+ibv(isv)*ip(is))               &
            +(spv(1:ifull)-spv(isv))*0.5)                        &
           /(ibu(1:ifull)+ibu(iwu)+ibv(1:ifull)+ibv(isv))
      end where
      if (icemode.gt.1) then
        ! cavitating fluid
        ipmax=27500.*ndic(1:ifull)*exp(-20.*(1.-nfracice(1:ifull)))
        nip=max(min(nip,ipmax),0.)
      else
        ! incompressible fluid
        nip=max(nip,0.)
      end if
      maxloclip=maxval(abs(nip-ip(1:ifull)))
      ip(1:ifull)=alpha*nip+(1.-alpha)*ip(1:ifull)

      call MPI_AllReduce(maxloclip,maxglobip,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
      if (maxglobip.lt.itol.and.ll.gt.2) exit
 
      itotits=itotits+1
    end do
  else
    maxloclip=0.
    maxglobip=0.
  end if

  ! Update ice velocity with plastic terms
  call bounds(ip,nrows=2)
  spu(1:ifull)=0.5*(ip(in)+ip(ine))
  squ(1:ifull)=0.5*(ip(is)+ip(ise))
  spv(1:ifull)=0.5*(ip(ie)+ip(ien))
  sqv(1:ifull)=0.5*(ip(iw)+ip(iwn))
  dipdxu=(ip(ie)-ip(1:ifull))*emu(1:ifull)/ds
  dipdyu=stwgt(:,1)*0.5*(tnu-tsu)*emu(1:ifull)/ds
  dipdxv=stwgt(:,2)*0.5*(tev-twv)*emv(1:ifull)/ds
  dipdyv=(ip(in)-ip(1:ifull))*emv(1:ifull)/ds
  niu(1:ifull)=(niu(1:ifull)+ibu*dipdxu+icu*dipdyu)*ee(1:ifull)*ee(ie)
  niv(1:ifull)=(niv(1:ifull)+ibv*dipdyv+icv*dipdxv)*ee(1:ifull)*ee(in)
  call boundsuv(niu,niv)

  ! ADVECT ICE ------------------------------------------------------
  ! use simple upwind scheme to improve volume conservation

  ! Horizontal advection for ice area
  ndum(1:ifull)=fracice/(em(1:ifull)*em(1:ifull)) ! ndum is an area
  call bounds(ndum)
  odum=ndum(1:ifull)
  odum=odum+0.5*(niu(iwu)*(ndum(1:ifull)+ndum(iw))-abs(niu(iwu))*(ndum(1:ifull)-ndum(iw)))*emu(iwu)/(ds/dt+abs(niu(iwu))*emu(iwu))
  odum=odum-0.5*(niu(1:ifull)*(ndum(1:ifull)+ndum(ie))+abs(niu(1:ifull))*(ndum(1:ifull)-ndum(ie)))*emu(1:ifull)/(ds/dt+abs(niu(1:ifull))*emu(1:ifull))
  odum=odum+0.5*(niv(isv)*(ndum(1:ifull)+ndum(is))-abs(niv(isv))*(ndum(1:ifull)-ndum(is)))*emv(isv)/(ds/dt+abs(niv(isv))*emu(isv))
  odum=odum-0.5*(niv(1:ifull)*(ndum(1:ifull)+ndum(in))+abs(niv(1:ifull))*(ndum(1:ifull)-ndum(in)))*emv(1:ifull)/(ds/dt+abs(niv(1:ifull))*emv(1:ifull))
  ndum(1:ifull)=odum
  ndum=max(ndum,0.)
  nfracice(1:ifull)=ndum(1:ifull)*em(1:ifull)*em(1:ifull)
  nfracice(1:ifull)=min(max(nfracice(1:ifull),0.),1.)

! Horizontal advection for ice volume
  ndum(1:ifull)=sicedep*fracice/(em(1:ifull)*em(1:ifull)) ! now ndum is a volume
  call bounds(ndum)
  odum=ndum(1:ifull)
  odum=odum+0.5*(niu(iwu)*(ndum(1:ifull)+ndum(iw))-abs(niu(iwu))*(ndum(1:ifull)-ndum(iw)))*emu(iwu)/(ds/dt+abs(niu(iwu))*emu(iwu))
  odum=odum-0.5*(niu(1:ifull)*(ndum(1:ifull)+ndum(ie))+abs(niu(1:ifull))*(ndum(1:ifull)-ndum(ie)))*emu(1:ifull)/(ds/dt+abs(niu(1:ifull))*emu(1:ifull))
  odum=odum+0.5*(niv(isv)*(ndum(1:ifull)+ndum(is))-abs(niv(isv))*(ndum(1:ifull)-ndum(is)))*emv(isv)/(ds/dt+abs(niv(isv))*emu(isv))
  odum=odum-0.5*(niv(1:ifull)*(ndum(1:ifull)+ndum(in))+abs(niv(1:ifull))*(ndum(1:ifull)-ndum(in)))*emv(1:ifull)/(ds/dt+abs(niv(1:ifull))*emv(1:ifull))
  ndum(1:ifull)=odum
  ndum=max(ndum,0.)
  ndic(1:ifull)=ndum(1:ifull)*em(1:ifull)*em(1:ifull)/max(nfracice(1:ifull),1.E-20)

! Horizontal advection for snow volume
  ndum(1:ifull)=snowd*0.001*fracice/(em(1:ifull)*em(1:ifull)) ! now ndum is a volume
  call bounds(ndum)
  odum=ndum(1:ifull)
  odum=odum+0.5*(niu(iwu)*(ndum(1:ifull)+ndum(iw))-abs(niu(iwu))*(ndum(1:ifull)-ndum(iw)))*emu(iwu)/(ds/dt+abs(niu(iwu))*emu(iwu))
  odum=odum-0.5*(niu(1:ifull)*(ndum(1:ifull)+ndum(ie))+abs(niu(1:ifull))*(ndum(1:ifull)-ndum(ie)))*emu(1:ifull)/(ds/dt+abs(niu(1:ifull))*emu(1:ifull))
  odum=odum+0.5*(niv(isv)*(ndum(1:ifull)+ndum(is))-abs(niv(isv))*(ndum(1:ifull)-ndum(is)))*emv(isv)/(ds/dt+abs(niv(isv))*emu(isv))
  odum=odum-0.5*(niv(1:ifull)*(ndum(1:ifull)+ndum(in))+abs(niv(1:ifull))*(ndum(1:ifull)-ndum(in)))*emv(1:ifull)/(ds/dt+abs(niv(1:ifull))*emv(1:ifull))
  ndum(1:ifull)=odum
  ndum=max(ndum,0.)
  ndsn(1:ifull)=ndum(1:ifull)*em(1:ifull)*em(1:ifull)/max(nfracice(1:ifull),1.E-20)

! Horizontal advection for ice energy store
  ndum(1:ifull)=i_sto
  call bounds(ndum)
  odum=ndum(1:ifull)
  odum=odum+0.5*(niu(iwu)*(ndum(1:ifull)+ndum(iw))-abs(niu(iwu))*(ndum(1:ifull)-ndum(iw)))*emu(iwu)/(ds/dt+abs(niu(iwu))*emu(iwu))
  odum=odum-0.5*(niu(1:ifull)*(ndum(1:ifull)+ndum(ie))+abs(niu(1:ifull))*(ndum(1:ifull)-ndum(ie)))*emu(1:ifull)/(ds/dt+abs(niu(1:ifull))*emu(1:ifull))
  odum=odum+0.5*(niv(isv)*(ndum(1:ifull)+ndum(is))-abs(niv(isv))*(ndum(1:ifull)-ndum(is)))*emv(isv)/(ds/dt+abs(niv(isv))*emu(isv))
  odum=odum-0.5*(niv(1:ifull)*(ndum(1:ifull)+ndum(in))+abs(niv(1:ifull))*(ndum(1:ifull)-ndum(in)))*emv(1:ifull)/(ds/dt+abs(niv(1:ifull))*emv(1:ifull))
  ndum(1:ifull)=odum
  ndum=max(ndum,0.)
  nsto(1:ifull)=ndum(1:ifull)

! Horizontal advection for surface temperature
  ndum(1:ifull)=i_it(1:ifull,1)*fracice/(em(1:ifull)*em(1:ifull))
  call bounds(ndum)
  odum=ndum(1:ifull)
  odum=odum+0.5*(niu(iwu)*(ndum(1:ifull)+ndum(iw))-abs(niu(iwu))*(ndum(1:ifull)-ndum(iw)))*emu(iwu)/(ds/dt+abs(niu(iwu))*emu(iwu))
  odum=odum-0.5*(niu(1:ifull)*(ndum(1:ifull)+ndum(ie))+abs(niu(1:ifull))*(ndum(1:ifull)-ndum(ie)))*emu(1:ifull)/(ds/dt+abs(niu(1:ifull))*emu(1:ifull))
  odum=odum+0.5*(niv(isv)*(ndum(1:ifull)+ndum(is))-abs(niv(isv))*(ndum(1:ifull)-ndum(is)))*emv(isv)/(ds/dt+abs(niv(isv))*emu(isv))
  odum=odum-0.5*(niv(1:ifull)*(ndum(1:ifull)+ndum(in))+abs(niv(1:ifull))*(ndum(1:ifull)-ndum(in)))*emv(1:ifull)/(ds/dt+abs(niv(1:ifull))*emv(1:ifull))
  ndum(1:ifull)=odum
  ndum=max(ndum,0.)
  nit(1:ifull,1)=ndum(1:ifull)*em(1:ifull)*em(1:ifull)/max(nfracice(1:ifull),1.E-20)
  
! Horizontal advection of snow temperature
  ndum(1:ifull)=i_it(1:ifull,2)*fracice*snowd*0.001/(em(1:ifull)*em(1:ifull))
  call bounds(ndum)
  odum=ndum(1:ifull)
  odum=odum+0.5*(niu(iwu)*(ndum(1:ifull)+ndum(iw))-abs(niu(iwu))*(ndum(1:ifull)-ndum(iw)))*emu(iwu)/(ds/dt+abs(niu(iwu))*emu(iwu))
  odum=odum-0.5*(niu(1:ifull)*(ndum(1:ifull)+ndum(ie))+abs(niu(1:ifull))*(ndum(1:ifull)-ndum(ie)))*emu(1:ifull)/(ds/dt+abs(niu(1:ifull))*emu(1:ifull))
  odum=odum+0.5*(niv(isv)*(ndum(1:ifull)+ndum(is))-abs(niv(isv))*(ndum(1:ifull)-ndum(is)))*emv(isv)/(ds/dt+abs(niv(isv))*emu(isv))
  odum=odum-0.5*(niv(1:ifull)*(ndum(1:ifull)+ndum(in))+abs(niv(1:ifull))*(ndum(1:ifull)-ndum(in)))*emv(1:ifull)/(ds/dt+abs(niv(1:ifull))*emv(1:ifull))
  ndum(1:ifull)=odum
  ndum=max(ndum,0.)
  nit(1:ifull,2)=ndum(1:ifull)*em(1:ifull)*em(1:ifull)/max(ndsn(1:ifull)*nfracice(1:ifull),1.E-20)
  
! Horizontal advection of ice temperature
  do ii=3,4
    ndum(1:ifull)=i_it(1:ifull,ii)*fracice*sicedep/(em(1:ifull)*em(1:ifull))
    call bounds(ndum)
    odum=ndum(1:ifull)
    odum=odum+0.5*(niu(iwu)*(ndum(1:ifull)+ndum(iw))-abs(niu(iwu))*(ndum(1:ifull)-ndum(iw)))*emu(iwu)/(ds/dt+abs(niu(iwu))*emu(iwu))
    odum=odum-0.5*(niu(1:ifull)*(ndum(1:ifull)+ndum(ie))+abs(niu(1:ifull))*(ndum(1:ifull)-ndum(ie)))*emu(1:ifull)/(ds/dt+abs(niu(1:ifull))*emu(1:ifull))
    odum=odum+0.5*(niv(isv)*(ndum(1:ifull)+ndum(is))-abs(niv(isv))*(ndum(1:ifull)-ndum(is)))*emv(isv)/(ds/dt+abs(niv(isv))*emu(isv))
    odum=odum-0.5*(niv(1:ifull)*(ndum(1:ifull)+ndum(in))+abs(niv(1:ifull))*(ndum(1:ifull)-ndum(in)))*emv(1:ifull)/(ds/dt+abs(niv(1:ifull))*emv(1:ifull))
    ndum(1:ifull)=odum
    ndum=max(ndum,0.)
    nit(1:ifull,ii)=ndum(1:ifull)*em(1:ifull)*em(1:ifull)/max(ndic(1:ifull)*nfracice(1:ifull),1.E-20)
  end do

  where (nfracice(1:ifull).lt.1.E-8)
    nfracice(1:ifull)=0.
    ndic(1:ifull)=0.
    ndsn(1:ifull)=0.
    nsto(1:ifull)=0.
    nit(1:ifull,1)=i_it(:,1)
    nit(1:ifull,2)=i_it(:,2)
    nit(1:ifull,3)=i_it(:,3)
    nit(1:ifull,4)=i_it(:,4)
  end where
  
  where (nit(1:ifull,1).lt.100.)
    nit(1:ifull,1)=i_it(:,1)
    nit(1:ifull,2)=i_it(:,2)
    nit(1:ifull,3)=i_it(:,3)
    nit(1:ifull,4)=i_it(:,4)
  end where

  ! unstagger ice velocities
  siu(:,1)=niu(1:ifull)
  siv(:,1)=niv(1:ifull)
  call mlounstaguv(siu,siv,nuh(:,1:1),nvh(:,1:1))
  niu(1:ifull)=nuh(:,1)
  niv(1:ifull)=nvh(:,1)
  
end do

! temperature conservation
dum=0.
delpos=0.
delneg=0.
nt=min(max(nt,250.),400.)
do ii=1,wlev
  where(wtr(1:ifull))
    dum(:,ii)=nt(1:ifull,ii)-w_t(1:ifull,ii) ! increments
  end where
  odum=dz(1:ifull,ii)*dum(:,ii)/dd(1:ifull)
  ! cannot use 3d version, since it is hardwired to dsig
  call ccglobal_posneg(odum,dumpp,dumpn)
  delpos=delpos+dumpp
  delneg=delneg+dumpn
end do
alph_p = -delneg/max(delpos,1.e-30)
alph_p = min(alph_p,sqrt(alph_p))  ! best option
alph_p = max(alph_p,1.E-20)
nt(1:ifull,:)=w_t(1:ifull,:)+alph_p*max(0.,dum)+min(0.,dum)/alph_p

! salinity conservation
dum=0.
delpos=0.
delneg=0.
ns=max(ns,0.)
do ii=1,wlev
  where(wtr(1:ifull).and.w_s(1:ifull,ii).gt.1.)
    dum(:,ii)=ns(1:ifull,ii)-w_s(1:ifull,ii) ! increments
  end where
  odum=dz(1:ifull,ii)*dum(:,ii)/dd(1:ifull)
  ! cannot use 3d version, since it is hardwired to dsig
  call ccglobal_posneg(odum,dumpp,dumpn)
  delpos=delpos+dumpp
  delneg=delneg+dumpn
end do
alph_p = -delneg/max(delpos,1.e-30)
alph_p = min(alph_p,sqrt(alph_p))  ! best option
alph_p = max(alph_p,1.E-20)
ns(1:ifull,:)=w_s(1:ifull,:)+alph_p*max(0.,dum)+min(0.,dum)/alph_p

if (myid==0.and.(ktau.le.100.or.maxglobseta.gt.tol)) then
  write(6,*) "MLODYNAMICS ",totits,maxglobseta,itotits,maxglobip
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
do ii=1,4
  call mloimpice(nit(1:ifull,ii),ii,0)
end do
call mloimpice(nfracice(1:ifull),5,0)
call mloimpice(ndic(1:ifull),6,0)
call mloimpice(ndsn(1:ifull),7,0)
call mloimpice(nsto(1:ifull),8,0)
call mloimpice(niu(1:ifull),9,0)
call mloimpice(niv(1:ifull),10,0)
where (wtr(1:ifull))
  fracice=nfracice(1:ifull)
  sicedep=ndic(1:ifull)
  snowd=ndsn(1:ifull)*1000.
end where

return
end subroutine mlohadv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate depature points for MLO semi-Lagrangian advection
! (This subroutine is based on depts.f)

subroutine mlodeps(ubar,vbar,nface,xg,yg,x3d,y3d,z3d,dep,wtr)

use cc_mpi
use mlo
use vecsuv_m
use xyzinfo_m

implicit none

include 'newmpar.h'
include 'const_phys.h'
include 'parm.h'

integer ii,intsch,n,kx
integer, dimension(:,:), intent(out) :: nface
real, dimension(ifull,size(nface,2)), intent(in) :: ubar,vbar
real, dimension(ifull,size(nface,2)), intent(out) :: xg,yg
real*8, dimension(ifull,size(nface,2)), intent(out) :: x3d,y3d,z3d
real, dimension(ifull+iextra,size(nface,2)), intent(in) :: dep
real, dimension(ifull,size(nface,2)) :: uc,vc,wc
real, dimension(ifull+iextra,size(nface,2)) :: temp
logical, dimension(ifull+iextra), intent(in) :: wtr
integer, parameter :: nguess = 2

kx=size(nface,2)
intsch=mod(ktau,2)

! departure point x, y, z is called x3d, y3d, z3d
! first find corresponding cartesian vels
do ii=1,kx
  uc(:,ii)=(ax(1:ifull)*ubar(:,ii)+bx(1:ifull)*vbar(:,ii))*dt/rearth ! unit sphere 
  vc(:,ii)=(ay(1:ifull)*ubar(:,ii)+by(1:ifull)*vbar(:,ii))*dt/rearth ! unit sphere 
  wc(:,ii)=(az(1:ifull)*ubar(:,ii)+bz(1:ifull)*vbar(:,ii))*dt/rearth ! unit sphere 
  x3d(:,ii)=x-uc(:,ii) ! 1st guess
  y3d(:,ii)=y-vc(:,ii)
  z3d(:,ii)=z-wc(:,ii)
end do

! convert to grid point numbering
do ii=1,kx
  call mlotoij5(x3d(:,ii),y3d(:,ii),z3d(:,ii),nface(:,ii),xg(:,ii),yg(:,ii))
end do
! Share off processor departure points.
call deptsync(nface,xg,yg)

do n=1,nguess
  temp(1:ifull,:) = uc
  call mloints(temp,dep,intsch,nface,xg,yg,2,wtr)
  do ii=1,kx
    x3d(:,ii) = x - 0.5*(uc(:,ii)+temp(1:ifull,ii)) ! n+1 guess
  end do
  temp(1:ifull,:) = vc
  call mloints(temp,dep,intsch,nface,xg,yg,2,wtr)
  do ii=1,kx
    y3d(:,ii) = y - 0.5*(vc(:,ii)+temp(1:ifull,ii)) ! n+1 guess
  end do
  temp(1:ifull,:) = wc
  call mloints(temp,dep,intsch,nface,xg,yg,2,wtr)
  do ii=1,kx
    z3d(:,ii) = z - 0.5*(wc(:,ii)+temp(1:ifull,ii)) ! n+1 guess
  end do

  do ii=1,kx
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

subroutine mloints(s,d,intsch,nface,xg,yg,nfield,wtr)

use cc_mpi
use indices_m
use mlo

implicit none

include 'newmpar.h'
include 'parm.h'
include 'parmhor.h'
include 'mpif.h'

integer, intent(in) :: nfield,intsch
integer idel,iq,jdel,nn,kx
integer i,j,k,n,ind,ip,jp,iproc,ierr
integer, dimension(:,:), intent(in) :: nface
real, dimension(ifull,size(nface,2)), intent(in) :: xg,yg
real, dimension(ifull+iextra,size(nface,2)), intent(inout) :: s
real, dimension(ifull+iextra,size(nface,2)), intent(in) :: d
real, dimension(-1:ipan+2,-1:jpan+2,1:npan,size(nface,2)) :: sx,dx
real, dimension(0:nproc-1,maxbuflen) :: sdin,ddin
real, dimension(4) :: r
real a3,a4,c1,c2,c3,c4,cmax,cmin,sss,xxg,yyg
real dxx,cxx
real, dimension(size(nface,2)) :: cc1,cc2,cc3,cc4,d1,d2,d3,d4
logical, dimension(ifull+iextra), intent(in) :: wtr
ind(i,j,n)=i+(j-1)*ipan+(n-1)*ipan*jpan  ! *** for n=1,npan

kx=size(nface,2)

call bounds(s,nrows=2)

!======================== start of intsch=1 section ====================
if(intsch==1)then

    do n=1,npan         ! first simple copy into larger array
      do j=1,jpan
        do i=1,ipan
          sx(i,j,n,:) = s(ind(i,j,n),:)
          dx(i,j,n,:) = d(ind(i,j,n),:)
        enddo         ! i loop
      enddo           ! j loop
    enddo             ! n loop
            
!   this is intsb           EW interps done first
!   first extend s arrays into sx - this one -1:il+2 & -1:il+2
    do n=1,npan
      do j=1,jpan
        sx(0,j,n,:)      = s(iw(ind(1,j,n)),:)
        sx(-1,j,n,:)     = s(iww(ind(1,j,n)),:)
        sx(ipan+1,j,n,:) = s(ie(ind(ipan,j,n)),:)
        sx(ipan+2,j,n,:) = s(iee(ind(ipan,j,n)),:)
        dx(0,j,n,:)      = d(iw(ind(1,j,n)),:)
        dx(-1,j,n,:)     = d(iww(ind(1,j,n)),:)
        dx(ipan+1,j,n,:) = d(ie(ind(ipan,j,n)),:)
        dx(ipan+2,j,n,:) = d(iee(ind(ipan,j,n)),:)
      enddo            ! j loop
      do i=1,ipan
        sx(i,0,n,:)      = s(is(ind(i,1,n)),:)
        sx(i,-1,n,:)     = s(iss(ind(i,1,n)),:)
        sx(i,jpan+1,n,:) = s(in(ind(i,jpan,n)),:)
        sx(i,jpan+2,n,:) = s(inn(ind(i,jpan,n)),:)
        dx(i,0,n,:)      = d(is(ind(i,1,n)),:)
        dx(i,-1,n,:)     = d(iss(ind(i,1,n)),:)
        dx(i,jpan+1,n,:) = d(in(ind(i,jpan,n)),:)
        dx(i,jpan+2,n,:) = d(inn(ind(i,jpan,n)),:)
      enddo            ! i loop
!        for ew interpolation, sometimes need (different from ns):
!            (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
!           (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)
      sx(-1,0,n,:)          = s(lwws(n),:)
      sx(0,0,n,:)           = s(lws(n),:)
      sx(0,-1,n,:)          = s(lwss(n),:)
      sx(ipan+1,0,n,:)      = s(les(n),:)
      sx(ipan+2,0,n,:)      = s(lees(n),:)
      sx(ipan+1,-1,n,:)     = s(less(n),:)
      sx(-1,jpan+1,n,:)     = s(lwwn(n),:)
      sx(0,jpan+2,n,:)      = s(lwnn(n),:)
      sx(ipan+2,jpan+1,n,:) = s(leen(n),:)
      sx(ipan+1,jpan+2,n,:) = s(lenn(n),:)
      sx(0,jpan+1,n,:)      = s(iwn(ind(1,jpan,n)),:)
      sx(ipan+1,jpan+1,n,:) = s(ien(ind(ipan,jpan,n)),:)
      dx(-1,0,n,:)          = d(lwws(n),:)
      dx(0,0,n,:)           = d(lws(n),:)
      dx(0,-1,n,:)          = d(lwss(n),:)
      dx(ipan+1,0,n,:)      = d(les(n),:)
      dx(ipan+2,0,n,:)      = d(lees(n),:)
      dx(ipan+1,-1,n,:)     = d(less(n),:)
      dx(-1,jpan+1,n,:)     = d(lwwn(n),:)
      dx(0,jpan+2,n,:)      = d(lwnn(n),:)
      dx(ipan+2,jpan+1,n,:) = d(leen(n),:)
      dx(ipan+1,jpan+2,n,:) = d(lenn(n),:)
      dx(0,jpan+1,n,:)      = d(iwn(ind(1,jpan,n)),:)
      dx(ipan+1,jpan+1,n,:) = d(ien(ind(ipan,jpan,n)),:)
    enddo               ! n loop
  
    if(nfield<mh_bs)then
     do iq=1,ifull    ! non Berm-Stan option
      if (wtr(iq)) then
       do k=1,kx           
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
        !if (intmode.eq.1) then
          cc1 = sx(idel-1,jdel,n,:)
          cc2 = sx(idel  ,jdel,n,:)
          cc3 = sx(idel+1,jdel,n,:)
          cc4 = sx(idel+2,jdel,n,:)
          d1 = dx(idel-1,jdel,n,:)
          d2 = dx(idel  ,jdel,n,:)
          d3 = dx(idel+1,jdel,n,:)
          d4 = dx(idel+2,jdel,n,:)
          dxx = d(iq,k)
            cxx = 0.
            call searchdeltab(cc1,d1,dxx,cxx,c1)
            call searchdeltab(cc2,d2,dxx,cxx,c2)
            call searchdeltab(cc3,d3,dxx,cxx,c3)
            call searchdeltab(cc4,d4,dxx,cxx,c4)
        !else
        !  c1 = sx(idel-1,jdel,n,k)
        !  c2 = sx(idel  ,jdel,n,k)
        !  c3 = sx(idel+1,jdel,n,k)
        !  c4 = sx(idel+2,jdel,n,k)        
        !end if
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(2) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
        else
          r(2) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.) &
               -xxg*(1.+xxg)*c4/3.)                          &
               +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
        endif         !  (mhint==2)
        !if (intmode.eq.1) then
          cc1 = sx(idel-1,jdel+1,n,:)
          cc2 = sx(idel  ,jdel+1,n,:)
          cc3 = sx(idel+1,jdel+1,n,:)
          cc4 = sx(idel+2,jdel+1,n,:)
          d1 = dx(idel-1,jdel+1,n,:)
          d2 = dx(idel  ,jdel+1,n,:)
          d3 = dx(idel+1,jdel+1,n,:)
          d4 = dx(idel+2,jdel+1,n,:)
            call searchdeltab(cc1,d1,dxx,cxx,c1)
            call searchdeltab(cc2,d2,dxx,cxx,c2)
            call searchdeltab(cc3,d3,dxx,cxx,c3)
            call searchdeltab(cc4,d4,dxx,cxx,c4)
        !else
        !  c1 = sx(idel-1,jdel+1,n,k)
        !  c2 = sx(idel  ,jdel+1,n,k)
        !  c3 = sx(idel+1,jdel+1,n,k)
        !  c4 = sx(idel+2,jdel+1,n,k)
        !end if
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(3) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
        else
          r(3) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.) &
               -xxg*(1.+xxg)*c4/3.)                          &
               +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
        endif         !  (mhint==2)
        !if (intmode.eq.1) then
            cc1 = sx(idel  ,jdel-1,n,:)
            cc2 = sx(idel+1,jdel-1,n,:)
            cc3 = sx(idel  ,jdel+2,n,:)
            cc4 = sx(idel+1,jdel+2,n,:)
            d1 = dx(idel  ,jdel-1,n,:)
            d2 = dx(idel+1,jdel-1,n,:)
            d3 = dx(idel  ,jdel+2,n,:)
            d4 = dx(idel+1,jdel+2,n,:)
              call searchdeltab(cc1,d1,dxx,cxx,c1)
              call searchdeltab(cc2,d2,dxx,cxx,c2)
              call searchdeltab(cc3,d3,dxx,cxx,c3)
              call searchdeltab(cc4,d4,dxx,cxx,c4)
            r(1) = (1.-xxg)*c1 +xxg*c2
            r(4) = (1.-xxg)*c3 +xxg*c4
        !else
        !  do nn=1,4,3   ! N.B.
        !    c2 = sx(idel  ,jdel+nn-2,n,k)
        !    c3 = sx(idel+1,jdel+nn-2,n,k)
        !    r(nn) = (1.-xxg)*c2 +xxg*c3
        !  enddo         ! nn loop
        !end if
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
       enddo
      endif           ! if (wtr(iq)) then
     enddo            ! iq loop
    else                ! (nfield<mh_bs)
     do iq=1,ifull    ! Berm-Stan option here e.g. qg & gases
      if (wtr(iq)) then
       do k=1,kx      
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
        !if (intmode.eq.1) then
          cc1 = sx(idel-1,jdel,n,:)
          cc2 = sx(idel  ,jdel,n,:)
          cc3 = sx(idel+1,jdel,n,:)
          cc4 = sx(idel+2,jdel,n,:)
          d1 = dx(idel-1,jdel,n,:)
          d2 = dx(idel  ,jdel,n,:)
          d3 = dx(idel+1,jdel,n,:)
          d4 = dx(idel+2,jdel,n,:)
          dxx = d(iq,k)
            cxx = 0.
            call searchdeltab(cc1,d1,dxx,cxx,c1)
            call searchdeltab(cc2,d2,dxx,cxx,c2)
            call searchdeltab(cc3,d3,dxx,cxx,c3)
            call searchdeltab(cc4,d4,dxx,cxx,c4)
        !else
        !  c1 = sx(idel-1,jdel,n,k)
        !  c2 = sx(idel  ,jdel,n,k)
        !  c3 = sx(idel+1,jdel,n,k)
        !  c4 = sx(idel+2,jdel,n,k)
        !end if
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
        !if (intmode.eq.1) then
          cc1 = sx(idel-1,jdel+1,n,:)
          cc2 = sx(idel  ,jdel+1,n,:)
          cc3 = sx(idel+1,jdel+1,n,:)
          cc4 = sx(idel+2,jdel+1,n,:)
          d1 = dx(idel-1,jdel+1,n,:)
          d2 = dx(idel  ,jdel+1,n,:)
          d3 = dx(idel+1,jdel+1,n,:)
          d4 = dx(idel+2,jdel+1,n,:)
            call searchdeltab(cc1,d1,dxx,cxx,c1)
            call searchdeltab(cc2,d2,dxx,cxx,c2)
            call searchdeltab(cc3,d3,dxx,cxx,c3)
            call searchdeltab(cc4,d4,dxx,cxx,c4)
        !else
        !  c1 = sx(idel-1,jdel+1,n,k)
        !  c2 = sx(idel  ,jdel+1,n,k)
        !  c3 = sx(idel+1,jdel+1,n,k)
        !  c4 = sx(idel+2,jdel+1,n,k)
        !end if
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
        !if (intmode.eq.1) then
            cc1 = sx(idel  ,jdel-1,n,:)
            cc2 = sx(idel+1,jdel-1,n,:)
            cc3 = sx(idel  ,jdel+2,n,:)
            cc4 = sx(idel+1,jdel+2,n,:)
            d1 = dx(idel  ,jdel-1,n,:)
            d2 = dx(idel+1,jdel-1,n,:)
            d3 = dx(idel  ,jdel+2,n,:)
            d4 = dx(idel+1,jdel+2,n,:)
              call searchdeltab(cc1,d1,dxx,cxx,c1)
              call searchdeltab(cc2,d2,dxx,cxx,c2)
              call searchdeltab(cc3,d3,dxx,cxx,c3)
              call searchdeltab(cc4,d4,dxx,cxx,c4)
            r(1) = (1.-xxg)*c1 +xxg*c2
            r(4) = (1.-xxg)*c2 +xxg*c4
        !else
        !  do nn=1,4,3   ! N.B.
        !    c2 = sx(idel  ,jdel+nn-2,n,k)
        !    c3 = sx(idel+1,jdel+nn-2,n,k)
        !    r(nn) = (1.-xxg)*c2 +xxg*c3
        !  enddo         ! nn loop
        !end if
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
       enddo
      endif           ! if (wtr) then
     enddo            ! iq loop
    endif               ! (nfield<mh_bs)  .. else ..
            
! Loop over points that need to be calculated for other processes

  !if (intmode.eq.1) then
    call intssend(d)
    do iproc=0,nproc-1
      ddin(iproc,1:drlen(iproc))=sextra(iproc)%a(1:drlen(iproc))
    end do
    call intssend(s)
    do iproc=0,nproc-1
      sdin(iproc,1:drlen(iproc))=sextra(iproc)%a(1:drlen(iproc))
    end do
  !end if

  if(nfield<mh_bs)then
    do iproc=0,nproc-1
      if ( iproc == myid ) then
        cycle
      end if
      do iq=1,drlen(iproc)
       if (ddin(iproc,iq).gt.0.1) then
        !  Convert face index from 0:npanels to array indices
        ip = min(il_g,max(1,nint(dpoints(iproc)%a(2,iq))))
        jp = min(il_g,max(1,nint(dpoints(iproc)%a(3,iq))))
        n = nint(dpoints(iproc)%a(1,iq)) + noff ! Local index
        !  Need global face index in fproc call
        idel = int(dpoints(iproc)%a(2,iq))
        xxg = dpoints(iproc)%a(2,iq) - idel
        jdel = int(dpoints(iproc)%a(3,iq))
        yyg = dpoints(iproc)%a(3,iq) - jdel
        k = nint(dpoints(iproc)%a(4,iq))
        idel = idel - ioff(n-noff)
        jdel = jdel - joff(n-noff)
        !if (intmode.eq.1) then
          cc1 = sx(idel-1,jdel,n,:)
          cc2 = sx(idel  ,jdel,n,:)
          cc3 = sx(idel+1,jdel,n,:)
          cc4 = sx(idel+2,jdel,n,:)
          d1 = dx(idel-1,jdel,n,:)
          d2 = dx(idel  ,jdel,n,:)
          d3 = dx(idel+1,jdel,n,:)
          d4 = dx(idel+2,jdel,n,:)
          dxx = d(iq,k)
            cxx = 0.
            call searchdeltab(cc1,d1,dxx,cxx,c1)
            call searchdeltab(cc2,d2,dxx,cxx,c2)
            call searchdeltab(cc3,d3,dxx,cxx,c3)
            call searchdeltab(cc4,d4,dxx,cxx,c4)
        !else
        !  c1 = sx(idel-1,jdel,n,k)
        !  c2 = sx(idel  ,jdel,n,k)
        !  c3 = sx(idel+1,jdel,n,k)
        !  c4 = sx(idel+2,jdel,n,k)
        !end if
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(2) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
        else
          r(2) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.) &
               -xxg*(1.+xxg)*c4/3.)                          &
               +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
        endif         !  (mhint==2)
        !if(intmode.eq.1) then
          cc1 = sx(idel-1,jdel+1,n,:)
          cc2 = sx(idel  ,jdel+1,n,:)
          cc3 = sx(idel+1,jdel+1,n,:)
          cc4 = sx(idel+2,jdel+1,n,:)
          d1 = dx(idel-1,jdel+1,n,:)
          d2 = dx(idel  ,jdel+1,n,:)
          d3 = dx(idel+1,jdel+1,n,:)
          d4 = dx(idel+2,jdel+1,n,:)
            call searchdeltab(cc1,d1,dxx,cxx,c1)
            call searchdeltab(cc2,d2,dxx,cxx,c2)
            call searchdeltab(cc3,d3,dxx,cxx,c3)
            call searchdeltab(cc4,d4,dxx,cxx,c4)
        !else
        !  c1 = sx(idel-1,jdel+1,n,k)
        !  c2 = sx(idel  ,jdel+1,n,k)
        !  c3 = sx(idel+1,jdel+1,n,k)
        !  c4 = sx(idel+2,jdel+1,n,k)
        !end if
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(3) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
        else
          r(3) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.) &
               -xxg*(1.+xxg)*c4/3.)                          &
               +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
        endif         !  (mhint==2)
        !if (intmode.eq.1) then
            cc1 = sx(idel  ,jdel-1,n,:)
            cc2 = sx(idel+1,jdel-1,n,:)
            cc3 = sx(idel  ,jdel+2,n,:)
            cc4 = sx(idel+1,jdel+2,n,:)
            d1 = dx(idel  ,jdel-1,n,:)
            d2 = dx(idel+1,jdel-1,n,:)
            d3 = dx(idel  ,jdel+2,n,:)
            d4 = dx(idel+1,jdel+2,n,:)
              call searchdeltab(cc1,d1,dxx,cxx,c1)
              call searchdeltab(cc2,d2,dxx,cxx,c2)
              call searchdeltab(cc3,d3,dxx,cxx,c3)
              call searchdeltab(cc4,d4,dxx,cxx,c4)
            r(1) = (1.-xxg)*c1 +xxg*c2
            r(4) = (1.-xxg)*c3 +xxg*c4
        !else
        !  do nn=1,4,3   ! N.B.
        !    c2 = sx(idel  ,jdel+nn-2,n,k)
        !    c3 = sx(idel+1,jdel+nn-2,n,k)
        !    r(nn) = (1.-xxg)*c2 +xxg*c3
        !  enddo         ! nn loop
        !end if
        if(mhint==2)then ! Bessel interp
          a4 = r(4)-r(1)+3.*(r(2)-r(3))
          a3 = r(1)-2.*r(2)+r(3)-a4
          sextra(iproc)%a(iq) = r(2) + 0.5*yyg*(r(3)-r(1) +yyg*(a3+yyg*a4))
        else
          sextra(iproc)%a(iq) = ((1.-yyg)*((2.-yyg)* &
            ((1.+yyg)*r(2)-yyg*r(1)/3.)              &
            -yyg*(1.+yyg)*r(4)/3.)                   &
            +yyg*(1.+yyg)*(2.-yyg)*r(3))/2.
        endif         !  (mhint==2)
       endif           ! if (ddin(iproc,iq).gt.0.1) then
      enddo            ! iq loop
    end do              ! iproc loop
  else                   ! (nfield<mh_bs)
    do iproc=0,nproc-1
      if ( iproc == myid ) then
        cycle
      end if
      do iq=1,drlen(iproc)
       if (ddin(iproc,iq).gt.0.1) then
        !  Convert face index from 0:npanels to array indices
        ip = min(il_g,max(1,nint(dpoints(iproc)%a(2,iq))))
        jp = min(il_g,max(1,nint(dpoints(iproc)%a(3,iq))))
        n = nint(dpoints(iproc)%a(1,iq)) + noff ! Local index
        !  Need global face index in fproc call
        idel = int(dpoints(iproc)%a(2,iq))
        xxg = dpoints(iproc)%a(2,iq) - idel
        jdel = int(dpoints(iproc)%a(3,iq))
        yyg = dpoints(iproc)%a(3,iq) - jdel
        k = nint(dpoints(iproc)%a(4,iq))
        idel = idel - ioff(n-noff)
        jdel = jdel - joff(n-noff)
        !if (intmode.eq.1) then
          cc1 = sx(idel-1,jdel,n,:)
          cc2 = sx(idel  ,jdel,n,:)
          cc3 = sx(idel+1,jdel,n,:)
          cc4 = sx(idel+2,jdel,n,:)
          d1 = dx(idel-1,jdel,n,:)
          d2 = dx(idel  ,jdel,n,:)
          d3 = dx(idel+1,jdel,n,:)
          d4 = dx(idel+2,jdel,n,:)
          dxx = d(iq,k)
            cxx = 0.
            call searchdeltab(cc1,d1,dxx,cxx,c1)
            call searchdeltab(cc2,d2,dxx,cxx,c2)
            call searchdeltab(cc3,d3,dxx,cxx,c3)
            call searchdeltab(cc4,d4,dxx,cxx,c4)
        !else
        !  c1 = sx(idel-1,jdel,n,k)
        !  c2 = sx(idel  ,jdel,n,k)
        !  c3 = sx(idel+1,jdel,n,k)
        !  c4 = sx(idel+2,jdel,n,k)
        !end if
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
        !if (intmode.eq.1) then
          cc1 = sx(idel-1,jdel+1,n,:)
          cc2 = sx(idel  ,jdel+1,n,:)
          cc3 = sx(idel+1,jdel+1,n,:)
          cc4 = sx(idel+2,jdel+1,n,:)
          d1 = dx(idel-1,jdel+1,n,:)
          d2 = dx(idel  ,jdel+1,n,:)
          d3 = dx(idel+1,jdel+1,n,:)
          d4 = dx(idel+2,jdel+1,n,:)
            call searchdeltab(cc1,d1,dxx,cxx,c1)
            call searchdeltab(cc2,d2,dxx,cxx,c2)
            call searchdeltab(cc3,d3,dxx,cxx,c3)
            call searchdeltab(cc4,d4,dxx,cxx,c4)
        !else
        !  c1 = sx(idel-1,jdel+1,n,k)
        !  c2 = sx(idel  ,jdel+1,n,k)
        !  c3 = sx(idel+1,jdel+1,n,k)
        !  c4 = sx(idel+2,jdel+1,n,k)
        !end if
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
        !if (intmode.eq.1) then
            cc1 = sx(idel  ,jdel-1,n,:)
            cc2 = sx(idel+1,jdel-1,n,:)
            cc3 = sx(idel  ,jdel+2,n,:)
            cc4 = sx(idel+1,jdel+2,n,:)
            d1 = dx(idel  ,jdel-1,n,:)
            d2 = dx(idel+1,jdel-1,n,:)
            d3 = dx(idel  ,jdel+2,n,:)
            d4 = dx(idel+1,jdel+2,n,:)
              call searchdeltab(cc1,d1,dxx,cxx,c1)
              call searchdeltab(cc2,d2,dxx,cxx,c2)
              call searchdeltab(cc3,d3,dxx,cxx,c3)
              call searchdeltab(cc4,d4,dxx,cxx,c4)
            r(1) = (1.-xxg)*c1 +xxg*c2
            r(4) = (1.-xxg)*c3 +xxg*c4
        !else
        !  do nn=1,4,3   ! N.B.
        !    c2 = sx(idel  ,jdel+nn-2,n,k)
        !    c3 = sx(idel+1,jdel+nn-2,n,k)
        !    r(nn) = (1.-xxg)*c2 +xxg*c3
        !  enddo         ! nn loop
        !end if
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
       endif            ! if (ddin(iproc,iq).gt.0.1) then
      enddo            ! iq loop
    end do              ! iproc loop
  endif                  ! (nfield<mh_bs)  .. else ..
!========================   end of intsch=1 section ====================
else     ! if(intsch==1)then
!======================== start of intsch=2 section ====================
  
!       this is intsc           NS interps done first
!       first extend s arrays into sx - this one -1:il+2 & -1:il+2

    do n=1,npan         ! first simple copy into larger array
      do j=1,jpan
        do i=1,ipan
          sx(i,j,n,:)=s(ind(i,j,n),:)
          dx(i,j,n,:)=d(ind(i,j,n),:)
        enddo         ! i loop
      enddo            ! j loop
    enddo               ! n loop

    do n=1,npan
      do j=1,jpan
        sx(0,j,n,:) = s(iw(ind(1,j,n)),:)
        sx(-1,j,n,:) = s(iww(ind(1,j,n)),:)
        sx(ipan+1,j,n,:) = s(ie(ind(ipan,j,n)),:)
        sx(ipan+2,j,n,:) = s(iee(ind(ipan,j,n)),:)
        dx(0,j,n,:) = d(iw(ind(1,j,n)),:)
        dx(-1,j,n,:) = d(iww(ind(1,j,n)),:)
        dx(ipan+1,j,n,:) = d(ie(ind(ipan,j,n)),:)
        dx(ipan+2,j,n,:) = d(iee(ind(ipan,j,n)),:)
      enddo            ! j loop
      do i=1,ipan
        sx(i,0,n,:) = s(is(ind(i,1,n)),:)
        sx(i,-1,n,:) = s(iss(ind(i,1,n)),:)
        sx(i,jpan+1,n,:) = s(in(ind(i,jpan,n)),:)
        sx(i,jpan+2,n,:) = s(inn(ind(i,jpan,n)),:)
        dx(i,0,n,:) = d(is(ind(i,1,n)),:)
        dx(i,-1,n,:) = d(iss(ind(i,1,n)),:)
        dx(i,jpan+1,n,:) = d(in(ind(i,jpan,n)),:)
        dx(i,jpan+2,n,:) = d(inn(ind(i,jpan,n)),:)
      enddo            ! i loop
!        for ns interpolation, sometimes need (different from ew):
!            (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
!          (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)
      sx(-1,0,n,:)=s(lsww(n),:)
      sx(0,0,n,:) = s(lsw(n),:)
      sx(0,-1,n,:) = s(lssw(n),:)
      sx(ipan+2,0,n,:) = s(lsee(n),:)
      sx(ipan+1,-1,n,:) = s(lsse(n),:)
      sx(-1,jpan+1,n,:) = s(lnww(n),:)
      sx(0,jpan+1,n,:) = s(lnw(n),:)
      sx(0,jpan+2,n,:) = s(lnnw(n),:)
      sx(ipan+2,jpan+1,n,:) = s(lnee(n),:)
      sx(ipan+1,jpan+2,n,:) = s(lnne(n),:)
      sx(ipan+1,0,n,:)    = s(ise(ind(ipan,1,n)),:)
      sx(ipan+1,jpan+1,n,:) = s(ine(ind(ipan,jpan,n)),:)
      dx(-1,0,n,:)=d(lsww(n),:)
      dx(0,0,n,:) = d(lsw(n),:)
      dx(0,-1,n,:) = d(lssw(n),:)
      dx(ipan+2,0,n,:) = d(lsee(n),:)
      dx(ipan+1,-1,n,:) = d(lsse(n),:)
      dx(-1,jpan+1,n,:) = d(lnww(n),:)
      dx(0,jpan+1,n,:) = d(lnw(n),:)
      dx(0,jpan+2,n,:) = d(lnnw(n),:)
      dx(ipan+2,jpan+1,n,:) = d(lnee(n),:)
      dx(ipan+1,jpan+2,n,:) = d(lnne(n),:)
      dx(ipan+1,0,n,:)    = d(ise(ind(ipan,1,n)),:)
      dx(ipan+1,jpan+1,n,:) = d(ine(ind(ipan,jpan,n)),:)
    enddo               ! n loop
  
    if(nfield<mh_bs)then
     do iq=1,ifull    ! non Berm-Stan option
      if (wtr(iq)) then
       do k=1,kx
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
        !if (intmode.eq.1) then
          cc1 = sx(idel,jdel-1,n,:)
          cc2 = sx(idel,jdel  ,n,:)
          cc3 = sx(idel,jdel+1,n,:)
          cc4 = sx(idel,jdel+2,n,:)
          d1 = dx(idel,jdel-1,n,:)
          d2 = dx(idel,jdel  ,n,:)
          d3 = dx(idel,jdel+1,n,:)
          d4 = dx(idel,jdel+2,n,:)
          dxx = d(iq,k)
            cxx = 0.
            call searchdeltab(cc1,d1,dxx,cxx,c1)
            call searchdeltab(cc2,d2,dxx,cxx,c2)
            call searchdeltab(cc3,d3,dxx,cxx,c3)
            call searchdeltab(cc4,d4,dxx,cxx,c4)
        !else
        !  c1 = sx(idel,jdel-1,n,k)
        !  c2 = sx(idel,jdel  ,n,k)
        !  c3 = sx(idel,jdel+1,n,k)
        !  c4 = sx(idel,jdel+2,n,k)
        !end if
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(2) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
        else
          r(2) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.) &
               -yyg*(1.+yyg)*c4/3.)                          &
               +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
        endif         !  (mhint==2)
        !if (intmode.eq.1) then
          cc1 = sx(idel+1,jdel-1,n,:)
          cc2 = sx(idel+1,jdel  ,n,:)
          cc3 = sx(idel+1,jdel+1,n,:)
          cc4 = sx(idel+1,jdel+2,n,:)
          d1 = dx(idel+1,jdel-1,n,:)
          d2 = dx(idel+1,jdel  ,n,:)
          d3 = dx(idel+1,jdel+1,n,:)
          d4 = dx(idel+1,jdel+2,n,:)
            call searchdeltab(cc1,d1,dxx,cxx,c1)
            call searchdeltab(cc2,d2,dxx,cxx,c2)
            call searchdeltab(cc3,d3,dxx,cxx,c3)
            call searchdeltab(cc4,d4,dxx,cxx,c4)
        !else
        !  c1 = sx(idel+1,jdel-1,n,k)
        !  c2 = sx(idel+1,jdel  ,n,k)
        !  c3 = sx(idel+1,jdel+1,n,k)
        !  c4 = sx(idel+1,jdel+2,n,k)
        !end if
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(3) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
        else
          r(3) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.) &
               -yyg*(1.+yyg)*c4/3.)                          &
               +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
        endif         !  (mhint==2)
        !if (intmode.eq.1) then
            cc1 = sx(idel-1,jdel  ,n,:)
            cc2 = sx(idel-1,jdel+1,n,:)
            cc3 = sx(idel+2,jdel  ,n,:)
            cc4 = sx(idel+2,jdel+1,n,:)
            d1 = dx(idel-1,jdel  ,n,:)
            d2 = dx(idel-1,jdel+1,n,:)
            d3 = dx(idel+2,jdel  ,n,:)
            d4 = dx(idel+2,jdel+1,n,:)
              call searchdeltab(cc1,d1,dxx,cxx,c1)
              call searchdeltab(cc2,d2,dxx,cxx,c2)
              call searchdeltab(cc3,d3,dxx,cxx,c3)
              call searchdeltab(cc4,d4,dxx,cxx,c4)
            r(1) = (1.-yyg)*c1 +yyg*c2
            r(4) = (1.-yyg)*c3 +yyg*c4
        !else
        !  do nn=1,4,3   ! N.B.
        !    c2 = sx(idel+nn-2,jdel  ,n,k)
        !    c3 = sx(idel+nn-2,jdel+1,n,k)
        !    r(nn) = (1.-yyg)*c2 +yyg*c3
        !  enddo         ! nn loop
        !end if
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
       enddo
      endif           ! if (wtr) then
     enddo            ! iq loop
    else                ! (nfield<mh_bs)
     do iq=1,ifull    ! Berm-Stan option here e.g. qg & gases
      if (wtr(iq)) then
       do k=1,kx
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
        !if (intmode.eq.1) then
          cc1 = sx(idel,jdel-1,n,:)
          cc2 = sx(idel,jdel  ,n,:)
          cc3 = sx(idel,jdel+1,n,:)
          cc4 = sx(idel,jdel+2,n,:)
          d1 = dx(idel,jdel-1,n,:)
          d2 = dx(idel,jdel  ,n,:)
          d3 = dx(idel,jdel+1,n,:)
          d4 = dx(idel,jdel+2,n,:)
          dxx = d(iq,k)
            cxx = 0.
            call searchdeltab(cc1,d1,dxx,cxx,c1)
            call searchdeltab(cc2,d2,dxx,cxx,c2)
            call searchdeltab(cc3,d3,dxx,cxx,c3)
            call searchdeltab(cc4,d4,dxx,cxx,c4)
        !else
        !  c1 = sx(idel,jdel-1,n,k)
        !  c2 = sx(idel,jdel  ,n,k)
        !  c3 = sx(idel,jdel+1,n,k)
        !  c4 = sx(idel,jdel+2,n,k)
        !end if
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
        !if (intmode.eq.1) then
          cc1 = sx(idel+1,jdel-1,n,:)
          cc2 = sx(idel+1,jdel  ,n,:)
          cc3 = sx(idel+1,jdel+1,n,:)
          cc4 = sx(idel+1,jdel+2,n,:)
          d1 = dx(idel+1,jdel-1,n,:)
          d2 = dx(idel+1,jdel  ,n,:)
          d3 = dx(idel+1,jdel+1,n,:)
          d4 = dx(idel+1,jdel+2,n,:)
            call searchdeltab(cc1,d1,dxx,cxx,c1)
            call searchdeltab(cc2,d2,dxx,cxx,c2)
            call searchdeltab(cc3,d3,dxx,cxx,c3)
            call searchdeltab(cc4,d4,dxx,cxx,c4)
        !else
        !  c1 = sx(idel+1,jdel-1,n,k)
        !  c2 = sx(idel+1,jdel  ,n,k)
        !  c3 = sx(idel+1,jdel+1,n,k)
        !  c4 = sx(idel+1,jdel+2,n,k)
        !end if
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
        !if (intmode.eq.1) then
            cc1 = sx(idel-1,jdel  ,n,:)
            cc2 = sx(idel-1,jdel+1,n,:)
            cc3 = sx(idel+2,jdel  ,n,:)
            cc4 = sx(idel+2,jdel+1,n,:)            
            d1 = dx(idel-1,jdel  ,n,:)
            d2 = dx(idel-1,jdel+1,n,:)
            d3 = dx(idel+2,jdel  ,n,:)
            d4 = dx(idel+2,jdel+1,n,:)
              call searchdeltab(cc1,d1,dxx,cxx,c1)
              call searchdeltab(cc2,d2,dxx,cxx,c2)
              call searchdeltab(cc3,d3,dxx,cxx,c3)
              call searchdeltab(cc4,d4,dxx,cxx,c4)
            r(1) = (1.-yyg)*c1 +yyg*c2
            r(4) = (1.-yyg)*c3 +yyg*c4
        !else
        !  do nn=1,4,3   ! N.B.
        !    c2 = sx(idel+nn-2,jdel  ,n,k)
        !    c3 = sx(idel+nn-2,jdel+1,n,k)
        !    r(nn) = (1.-yyg)*c2 +yyg*c3
        !  enddo         ! nn loop
        !end if
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
       end do
      endif           ! if (wtr) then
     enddo            ! iq loop
    endif               ! (nfield<mh_bs)  .. else ..

! For other processes

  !if (intmode.eq.1) then
    call intssend(d)
    do iproc=0,nproc-1
      ddin(iproc,1:drlen(iproc))=sextra(iproc)%a(1:drlen(iproc))
    end do
    call intssend(s)
    do iproc=0,nproc-1
      sdin(iproc,1:drlen(iproc))=sextra(iproc)%a(1:drlen(iproc))
    end do
  !end if

  if(nfield<mh_bs)then
    do iproc=0,nproc-1
      if ( iproc == myid ) then
        cycle
      end if
      do iq=1,drlen(iproc)
       if (ddin(iproc,iq).gt.0.1) then
        !  Convert face index from 0:npanels to array indices
        ip = min(il_g,max(1,nint(dpoints(iproc)%a(2,iq))))
        jp = min(il_g,max(1,nint(dpoints(iproc)%a(3,iq))))
        n = nint(dpoints(iproc)%a(1,iq)) + noff ! Local index
        !  Need global face index in fproc call
        idel = int(dpoints(iproc)%a(2,iq))
        xxg = dpoints(iproc)%a(2,iq) - idel
        jdel = int(dpoints(iproc)%a(3,iq))
        yyg = dpoints(iproc)%a(3,iq) - jdel
        k = nint(dpoints(iproc)%a(4,iq))
        idel = idel - ioff(n-noff)
        jdel = jdel - joff(n-noff)
        !if (intmode.eq.1) then
          cc1 = sx(idel,jdel-1,n,:)
          cc2 = sx(idel,jdel  ,n,:)
          cc3 = sx(idel,jdel+1,n,:)
          cc4 = sx(idel,jdel+2,n,:)
          d1 = dx(idel,jdel-1,n,:)
          d2 = dx(idel,jdel  ,n,:)
          d3 = dx(idel,jdel+1,n,:)
          d4 = dx(idel,jdel+2,n,:)
          dxx = d(iq,k)
            cxx = 0.
            call searchdeltab(cc1,d1,dxx,cxx,c1)
            call searchdeltab(cc2,d2,dxx,cxx,c2)
            call searchdeltab(cc3,d3,dxx,cxx,c3)
            call searchdeltab(cc4,d4,dxx,cxx,c4)
        !else
        !  c1 = sx(idel,jdel-1,n,k)
        !  c2 = sx(idel,jdel  ,n,k)
        !  c3 = sx(idel,jdel+1,n,k)
        !  c4 = sx(idel,jdel+2,n,k)
        !end if
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(2) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
        else
          r(2) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.) &
               -yyg*(1.+yyg)*c4/3.)                          &
               +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
        endif         !  (mhint==2)
        !if (intmode.eq.1) then
          cc1 = sx(idel+1,jdel-1,n,:)
          cc2 = sx(idel+1,jdel  ,n,:)
          cc3 = sx(idel+1,jdel+1,n,:)
          cc4 = sx(idel+1,jdel+2,n,:)
          d1 = dx(idel+1,jdel-1,n,:)
          d2 = dx(idel+1,jdel  ,n,:)
          d3 = dx(idel+1,jdel+1,n,:)
          d4 = dx(idel+1,jdel+2,n,:)
            call searchdeltab(cc1,d1,dxx,cxx,c1)
            call searchdeltab(cc2,d2,dxx,cxx,c2)
            call searchdeltab(cc3,d3,dxx,cxx,c3)
            call searchdeltab(cc4,d4,dxx,cxx,c4)
        !else
        !  c1 = sx(idel+1,jdel-1,n,k)
        !  c2 = sx(idel+1,jdel  ,n,k)
        !  c3 = sx(idel+1,jdel+1,n,k)
        !  c4 = sx(idel+1,jdel+2,n,k)
        !end if
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(3) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
        else
          r(3) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.) &
               -yyg*(1.+yyg)*c4/3.)                          &
               +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
        endif         !  (mhint==2)
        !if (intmode.eq.1) then
            cc1 = sx(idel-1,jdel  ,n,:)
            cc2 = sx(idel-1,jdel+1,n,:)
            cc3 = sx(idel+2,jdel  ,n,:)
            cc4 = sx(idel+2,jdel+1,n,:)
            d1 = dx(idel-1,jdel  ,n,:)
            d2 = dx(idel-1,jdel+1,n,:)
            d3 = dx(idel+2,jdel  ,n,:)
            d4 = dx(idel+2,jdel+1,n,:)
              call searchdeltab(cc1,d1,dxx,cxx,c1)
              call searchdeltab(cc2,d2,dxx,cxx,c2)
              call searchdeltab(cc3,d3,dxx,cxx,c3)
              call searchdeltab(cc4,d4,dxx,cxx,c4)
            r(1) = (1.-yyg)*c1 +yyg*c2
            r(4) = (1.-yyg)*c3 +yyg*c4
        !else
        !  do nn=1,4,3   ! N.B.
        !    c2 = sx(idel+nn-2,jdel  ,n,k)
        !    c3 = sx(idel+nn-2,jdel+1,n,k)
        !    r(nn) = (1.-yyg)*c2 +yyg*c3
        !  enddo         ! nn loop
        !end if
        if(mhint==2)then ! Bessel interp
          a4 = r(4)-r(1)+3.*(r(2)-r(3))
          a3 = r(1)-2.*r(2)+r(3)-a4
          sextra(iproc)%a(iq) = r(2)+ 0.5*xxg*(r(3)-r(1) +xxg*(a3+xxg*a4))
        else
          sextra(iproc)%a(iq) = ((1.-xxg)*((2.-xxg)*  &
            ((1.+xxg)*r(2)-xxg*r(1)/3.)               &
            -xxg*(1.+xxg)*r(4)/3.)                    &
            +xxg*(1.+xxg)*(2.-xxg)*r(3))/2.
        endif         !  (mhint==2)
       endif           ! if (ddin(iproc,iq).gt.0.1) then
      enddo            ! iq loop
    end do              ! iproc
  else                   ! (nfield<mh_bs)
    do iproc=0,nproc-1
      if ( iproc == myid ) then
        cycle
      end if
      do iq=1,drlen(iproc)
       if (ddin(iproc,iq).gt.0.1) then
        !  Convert face index from 0:npanels to array indices
        ip = min(il_g,max(1,nint(dpoints(iproc)%a(2,iq))))
        jp = min(il_g,max(1,nint(dpoints(iproc)%a(3,iq))))
        n = nint(dpoints(iproc)%a(1,iq)) + noff ! Local index
        !  Need global face index in fproc call
        idel = int(dpoints(iproc)%a(2,iq))
        xxg = dpoints(iproc)%a(2,iq) - idel
        jdel = int(dpoints(iproc)%a(3,iq))
        yyg = dpoints(iproc)%a(3,iq) - jdel
        k = nint(dpoints(iproc)%a(4,iq))
        idel = idel - ioff(n-noff)
        jdel = jdel - joff(n-noff)
        !if (intmode.eq.1) then
          cc1 = sx(idel,jdel-1,n,:)
          cc2 = sx(idel,jdel  ,n,:)
          cc3 = sx(idel,jdel+1,n,:)
          cc4 = sx(idel,jdel+2,n,:)
          d1 = dx(idel,jdel-1,n,:)
          d2 = dx(idel,jdel  ,n,:)
          d3 = dx(idel,jdel+1,n,:)
          d4 = dx(idel,jdel+2,n,:)
          dxx = d(iq,k)
            cxx = 0.
            call searchdeltab(cc1,d1,dxx,cxx,c1)
            call searchdeltab(cc2,d2,dxx,cxx,c2)
            call searchdeltab(cc3,d3,dxx,cxx,c3)
            call searchdeltab(cc4,d4,dxx,cxx,c4)
        !else
        !  c1 = sx(idel,jdel-1,n,k)
        !  c2 = sx(idel,jdel  ,n,k)
        !  c3 = sx(idel,jdel+1,n,k)
        !  c4 = sx(idel,jdel+2,n,k)
        !end if
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
        !if (intmode.eq.1) then
          cc1 = sx(idel+1,jdel-1,n,:)
          cc2 = sx(idel+1,jdel  ,n,:)
          cc3 = sx(idel+1,jdel+1,n,:)
          cc4 = sx(idel+1,jdel+2,n,:)
          d1 = dx(idel+1,jdel-1,n,:)
          d2 = dx(idel+1,jdel  ,n,:)
          d3 = dx(idel+1,jdel+1,n,:)
          d4 = dx(idel+1,jdel+2,n,:)
            call searchdeltab(cc1,d1,dxx,cxx,c1)
            call searchdeltab(cc2,d2,dxx,cxx,c2)
            call searchdeltab(cc3,d3,dxx,cxx,c3)
            call searchdeltab(cc4,d4,dxx,cxx,c4)
        !else
        !  c1 = sx(idel+1,jdel-1,n,k)
        !  c2 = sx(idel+1,jdel  ,n,k)
        !  c3 = sx(idel+1,jdel+1,n,k)
        !  c4 = sx(idel+1,jdel+2,n,k)
        !end if
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
        !if (intmode.eq.1) then
            cc1 = sx(idel-1,jdel  ,n,:)
            cc2 = sx(idel-1,jdel+1,n,:)
            cc3 = sx(idel+2,jdel  ,n,:)
            cc4 = sx(idel+2,jdel+1,n,:)
            d1 = dx(idel-1,jdel  ,n,:)
            d2 = dx(idel-1,jdel+1,n,:)
            d3 = dx(idel+2,jdel  ,n,:)
            d4 = dx(idel+2,jdel+1,n,:)
              call searchdeltab(cc1,d1,dxx,cxx,c1)
              call searchdeltab(cc2,d2,dxx,cxx,c2)
              call searchdeltab(cc3,d3,dxx,cxx,c3)
              call searchdeltab(cc4,d4,dxx,cxx,c4)
            r(1) = (1.-yyg)*c1 +yyg*c2
            r(4) = (1.-yyg)*c3 +yyg*c4
        !else
        !  do nn=1,4,3   ! N.B.
        !    c2 = sx(idel+nn-2,jdel  ,n,k)
        !    c3 = sx(idel+nn-2,jdel+1,n,k)
        !    r(nn) = (1.-yyg)*c2 +yyg*c3
        !  enddo         ! nn loop
        !end if
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
       endif            ! if (ddin(iproc,iq).gt.0.1) then
      enddo            ! iq loop
    end do              ! iproc
  endif                  ! (nfield<mh_bs)  .. else ..

endif                     ! (intsch==1) .. else ..
!========================   end of intsch=1 section ====================

call intssync(s)

return
end subroutine mloints

subroutine mlointsb(s,d,intsch,nface,xg,yg,nfield,wtr)

use cc_mpi
use indices_m
use mlo

implicit none

include 'newmpar.h'
include 'parm.h'
include 'parmhor.h'
include 'mpif.h'

integer, intent(in) :: nfield,intsch
integer idel,iq,jdel,nn,kx
integer i,j,k,n,ind,ip,jp,iproc,ierr
integer, dimension(:,:), intent(in) :: nface
real, dimension(ifull,size(nface,2)), intent(in) :: xg,yg
real, dimension(ifull+iextra,size(nface,2)), intent(inout) :: s
real, dimension(ifull+iextra,size(nface,2)), intent(in) :: d
real, dimension(-1:ipan+2,-1:jpan+2,1:npan,size(nface,2)) :: sx,dx
real, dimension(0:nproc-1,maxbuflen) :: sdin,ddin
real, dimension(4) :: r
real a3,a4,c1,c2,c3,c4,cmax,cmin,sss,xxg,yyg
real dxx,cxx
real, dimension(size(nface,2)) :: cc1,cc2,cc3,cc4,d1,d2,d3,d4
logical, dimension(ifull+iextra), intent(in) :: wtr
ind(i,j,n)=i+(j-1)*ipan+(n-1)*ipan*jpan  ! *** for n=1,npan

kx=size(nface,2)

call bounds(s,nrows=2)

!======================== start of intsch=1 section ====================
if(intsch==1)then

    do n=1,npan         ! first simple copy into larger array
      do j=1,jpan
        do i=1,ipan
          sx(i,j,n,:) = s(ind(i,j,n),:)
          dx(i,j,n,:) = d(ind(i,j,n),:)
        enddo         ! i loop
      enddo           ! j loop
    enddo             ! n loop
            
!   this is intsb           EW interps done first
!   first extend s arrays into sx - this one -1:il+2 & -1:il+2
    do n=1,npan
      do j=1,jpan
        sx(0,j,n,:)      = s(iw(ind(1,j,n)),:)
        sx(-1,j,n,:)     = s(iww(ind(1,j,n)),:)
        sx(ipan+1,j,n,:) = s(ie(ind(ipan,j,n)),:)
        sx(ipan+2,j,n,:) = s(iee(ind(ipan,j,n)),:)
        dx(0,j,n,:)      = d(iw(ind(1,j,n)),:)
        dx(-1,j,n,:)     = d(iww(ind(1,j,n)),:)
        dx(ipan+1,j,n,:) = d(ie(ind(ipan,j,n)),:)
        dx(ipan+2,j,n,:) = d(iee(ind(ipan,j,n)),:)
      enddo            ! j loop
      do i=1,ipan
        sx(i,0,n,:)      = s(is(ind(i,1,n)),:)
        sx(i,-1,n,:)     = s(iss(ind(i,1,n)),:)
        sx(i,jpan+1,n,:) = s(in(ind(i,jpan,n)),:)
        sx(i,jpan+2,n,:) = s(inn(ind(i,jpan,n)),:)
        dx(i,0,n,:)      = d(is(ind(i,1,n)),:)
        dx(i,-1,n,:)     = d(iss(ind(i,1,n)),:)
        dx(i,jpan+1,n,:) = d(in(ind(i,jpan,n)),:)
        dx(i,jpan+2,n,:) = d(inn(ind(i,jpan,n)),:)
      enddo            ! i loop
!        for ew interpolation, sometimes need (different from ns):
!            (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
!           (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)
      sx(-1,0,n,:)          = s(lwws(n),:)
      sx(0,0,n,:)           = s(lws(n),:)
      sx(0,-1,n,:)          = s(lwss(n),:)
      sx(ipan+1,0,n,:)      = s(les(n),:)
      sx(ipan+2,0,n,:)      = s(lees(n),:)
      sx(ipan+1,-1,n,:)     = s(less(n),:)
      sx(-1,jpan+1,n,:)     = s(lwwn(n),:)
      sx(0,jpan+2,n,:)      = s(lwnn(n),:)
      sx(ipan+2,jpan+1,n,:) = s(leen(n),:)
      sx(ipan+1,jpan+2,n,:) = s(lenn(n),:)
      sx(0,jpan+1,n,:)      = s(iwn(ind(1,jpan,n)),:)
      sx(ipan+1,jpan+1,n,:) = s(ien(ind(ipan,jpan,n)),:)
      dx(-1,0,n,:)          = d(lwws(n),:)
      dx(0,0,n,:)           = d(lws(n),:)
      dx(0,-1,n,:)          = d(lwss(n),:)
      dx(ipan+1,0,n,:)      = d(les(n),:)
      dx(ipan+2,0,n,:)      = d(lees(n),:)
      dx(ipan+1,-1,n,:)     = d(less(n),:)
      dx(-1,jpan+1,n,:)     = d(lwwn(n),:)
      dx(0,jpan+2,n,:)      = d(lwnn(n),:)
      dx(ipan+2,jpan+1,n,:) = d(leen(n),:)
      dx(ipan+1,jpan+2,n,:) = d(lenn(n),:)
      dx(0,jpan+1,n,:)      = d(iwn(ind(1,jpan,n)),:)
      dx(ipan+1,jpan+1,n,:) = d(ien(ind(ipan,jpan,n)),:)
    enddo               ! n loop
  
    if(nfield<mh_bs)then
     do iq=1,ifull    ! non Berm-Stan option
      if (wtr(iq)) then
       do k=1,kx           
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
        !if (intmode.eq.1) then
          cc1 = sx(idel-1,jdel,n,:)
          cc2 = sx(idel  ,jdel,n,:)
          cc3 = sx(idel+1,jdel,n,:)
          cc4 = sx(idel+2,jdel,n,:)
          d1 = dx(idel-1,jdel,n,:)
          d2 = dx(idel  ,jdel,n,:)
          d3 = dx(idel+1,jdel,n,:)
          d4 = dx(idel+2,jdel,n,:)
          dxx = d(iq,k)
            cxx = s(iq,k)
            call searchdeltab(cc1,d1,dxx,cxx,c1)
            call searchdeltab(cc4,d4,dxx,cxx,c4)
            call searchdeltab(cc2,d2,dxx,c1,c2)
            call searchdeltab(cc3,d3,dxx,c4,c3)
        !else
        !  c1 = sx(idel-1,jdel,n,k)
        !  c2 = sx(idel  ,jdel,n,k)
        !  c3 = sx(idel+1,jdel,n,k)
        !  c4 = sx(idel+2,jdel,n,k)        
        !end if
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(2) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
        else
          r(2) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.) &
               -xxg*(1.+xxg)*c4/3.)                          &
               +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
        endif         !  (mhint==2)
        !if (intmode.eq.1) then
          cc1 = sx(idel-1,jdel+1,n,:)
          cc2 = sx(idel  ,jdel+1,n,:)
          cc3 = sx(idel+1,jdel+1,n,:)
          cc4 = sx(idel+2,jdel+1,n,:)
          d1 = dx(idel-1,jdel+1,n,:)
          d2 = dx(idel  ,jdel+1,n,:)
          d3 = dx(idel+1,jdel+1,n,:)
          d4 = dx(idel+2,jdel+1,n,:)
            call searchdeltab(cc1,d1,dxx,cxx,c1)
            call searchdeltab(cc4,d4,dxx,cxx,c4)
            call searchdeltab(cc2,d2,dxx,c1,c2)
            call searchdeltab(cc3,d3,dxx,c4,c3)
        !else
        !  c1 = sx(idel-1,jdel+1,n,k)
        !  c2 = sx(idel  ,jdel+1,n,k)
        !  c3 = sx(idel+1,jdel+1,n,k)
        !  c4 = sx(idel+2,jdel+1,n,k)
        !end if
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(3) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
        else
          r(3) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.) &
               -xxg*(1.+xxg)*c4/3.)                          &
               +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
        endif         !  (mhint==2)
        !if (intmode.eq.1) then
            cc1 = sx(idel  ,jdel-1,n,:)
            cc2 = sx(idel+1,jdel-1,n,:)
            cc3 = sx(idel  ,jdel+2,n,:)
            cc4 = sx(idel+1,jdel+2,n,:)
            d1 = dx(idel  ,jdel-1,n,:)
            d2 = dx(idel+1,jdel-1,n,:)
            d3 = dx(idel  ,jdel+2,n,:)
            d4 = dx(idel+1,jdel+2,n,:)
              call searchdeltab(cc1,d1,dxx,cxx,c1)
              call searchdeltab(cc4,d4,dxx,cxx,c4)
              call searchdeltab(cc2,d2,dxx,c1,c2)
              call searchdeltab(cc3,d3,dxx,c4,c3)
            r(1) = (1.-xxg)*c1 +xxg*c2
            r(4) = (1.-xxg)*c3 +xxg*c4
        !else
        !  do nn=1,4,3   ! N.B.
        !    c2 = sx(idel  ,jdel+nn-2,n,k)
        !    c3 = sx(idel+1,jdel+nn-2,n,k)
        !    r(nn) = (1.-xxg)*c2 +xxg*c3
        !  enddo         ! nn loop
        !end if
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
       enddo
      endif           ! if (wtr(iq)) then
     enddo            ! iq loop
    else                ! (nfield<mh_bs)
     do iq=1,ifull    ! Berm-Stan option here e.g. qg & gases
      if (wtr(iq)) then
       do k=1,kx      
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
        !if (intmode.eq.1) then
          cc1 = sx(idel-1,jdel,n,:)
          cc2 = sx(idel  ,jdel,n,:)
          cc3 = sx(idel+1,jdel,n,:)
          cc4 = sx(idel+2,jdel,n,:)
          d1 = dx(idel-1,jdel,n,:)
          d2 = dx(idel  ,jdel,n,:)
          d3 = dx(idel+1,jdel,n,:)
          d4 = dx(idel+2,jdel,n,:)
          dxx = d(iq,k)
            cxx = s(iq,k)
            call searchdeltab(cc1,d1,dxx,cxx,c1)
            call searchdeltab(cc4,d4,dxx,cxx,c4)
            call searchdeltab(cc2,d2,dxx,c1,c2)
            call searchdeltab(cc3,d3,dxx,c4,c3)
        !else
        !  c1 = sx(idel-1,jdel,n,k)
        !  c2 = sx(idel  ,jdel,n,k)
        !  c3 = sx(idel+1,jdel,n,k)
        !  c4 = sx(idel+2,jdel,n,k)
        !end if
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
        !if (intmode.eq.1) then
          cc1 = sx(idel-1,jdel+1,n,:)
          cc2 = sx(idel  ,jdel+1,n,:)
          cc3 = sx(idel+1,jdel+1,n,:)
          cc4 = sx(idel+2,jdel+1,n,:)
          d1 = dx(idel-1,jdel+1,n,:)
          d2 = dx(idel  ,jdel+1,n,:)
          d3 = dx(idel+1,jdel+1,n,:)
          d4 = dx(idel+2,jdel+1,n,:)
            call searchdeltab(cc1,d1,dxx,cxx,c1)
            call searchdeltab(cc4,d4,dxx,cxx,c4)
            call searchdeltab(cc2,d2,dxx,c1,c2)
            call searchdeltab(cc3,d3,dxx,c4,c3)
        !else
        !  c1 = sx(idel-1,jdel+1,n,k)
        !  c2 = sx(idel  ,jdel+1,n,k)
        !  c3 = sx(idel+1,jdel+1,n,k)
        !  c4 = sx(idel+2,jdel+1,n,k)
        !end if
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
        !if (intmode.eq.1) then
            cc1 = sx(idel  ,jdel-1,n,:)
            cc2 = sx(idel+1,jdel-1,n,:)
            cc3 = sx(idel  ,jdel+2,n,:)
            cc4 = sx(idel+1,jdel+2,n,:)
            d1 = dx(idel  ,jdel-1,n,:)
            d2 = dx(idel+1,jdel-1,n,:)
            d3 = dx(idel  ,jdel+2,n,:)
            d4 = dx(idel+1,jdel+2,n,:)
              call searchdeltab(cc1,d1,dxx,cxx,c1)
              call searchdeltab(cc4,d4,dxx,cxx,c4)
              call searchdeltab(cc2,d2,dxx,c1,c2)
              call searchdeltab(cc3,d3,dxx,c4,c3)
            r(1) = (1.-xxg)*c1 +xxg*c2
            r(4) = (1.-xxg)*c2 +xxg*c4
        !else
        !  do nn=1,4,3   ! N.B.
        !    c2 = sx(idel  ,jdel+nn-2,n,k)
        !    c3 = sx(idel+1,jdel+nn-2,n,k)
        !    r(nn) = (1.-xxg)*c2 +xxg*c3
        !  enddo         ! nn loop
        !end if
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
       enddo
      endif           ! if (wtr) then
     enddo            ! iq loop
    endif               ! (nfield<mh_bs)  .. else ..
            
! Loop over points that need to be calculated for other processes

  !if (intmode.eq.1) then
    call intssend(d)
    do iproc=0,nproc-1
      ddin(iproc,1:drlen(iproc))=sextra(iproc)%a(1:drlen(iproc))
    end do
    call intssend(s)
    do iproc=0,nproc-1
      sdin(iproc,1:drlen(iproc))=sextra(iproc)%a(1:drlen(iproc))
    end do
  !end if

  if(nfield<mh_bs)then
    do iproc=0,nproc-1
      if ( iproc == myid ) then
        cycle
      end if
      do iq=1,drlen(iproc)
       if (ddin(iproc,iq).gt.0.1) then
        !  Convert face index from 0:npanels to array indices
        ip = min(il_g,max(1,nint(dpoints(iproc)%a(2,iq))))
        jp = min(il_g,max(1,nint(dpoints(iproc)%a(3,iq))))
        n = nint(dpoints(iproc)%a(1,iq)) + noff ! Local index
        !  Need global face index in fproc call
        idel = int(dpoints(iproc)%a(2,iq))
        xxg = dpoints(iproc)%a(2,iq) - idel
        jdel = int(dpoints(iproc)%a(3,iq))
        yyg = dpoints(iproc)%a(3,iq) - jdel
        k = nint(dpoints(iproc)%a(4,iq))
        idel = idel - ioff(n-noff)
        jdel = jdel - joff(n-noff)
        !if (intmode.eq.1) then
          cc1 = sx(idel-1,jdel,n,:)
          cc2 = sx(idel  ,jdel,n,:)
          cc3 = sx(idel+1,jdel,n,:)
          cc4 = sx(idel+2,jdel,n,:)
          d1 = dx(idel-1,jdel,n,:)
          d2 = dx(idel  ,jdel,n,:)
          d3 = dx(idel+1,jdel,n,:)
          d4 = dx(idel+2,jdel,n,:)
          dxx = d(iq,k)
            cxx = s(iq,k)
            call searchdeltab(cc1,d1,dxx,cxx,c1)
            call searchdeltab(cc4,d4,dxx,cxx,c4)
            call searchdeltab(cc2,d2,dxx,c1,c2)
            call searchdeltab(cc3,d3,dxx,c4,c3)
        !else
        !  c1 = sx(idel-1,jdel,n,k)
        !  c2 = sx(idel  ,jdel,n,k)
        !  c3 = sx(idel+1,jdel,n,k)
        !  c4 = sx(idel+2,jdel,n,k)
        !end if
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(2) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
        else
          r(2) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.) &
               -xxg*(1.+xxg)*c4/3.)                          &
               +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
        endif         !  (mhint==2)
        !if(intmode.eq.1) then
          cc1 = sx(idel-1,jdel+1,n,:)
          cc2 = sx(idel  ,jdel+1,n,:)
          cc3 = sx(idel+1,jdel+1,n,:)
          cc4 = sx(idel+2,jdel+1,n,:)
          d1 = dx(idel-1,jdel+1,n,:)
          d2 = dx(idel  ,jdel+1,n,:)
          d3 = dx(idel+1,jdel+1,n,:)
          d4 = dx(idel+2,jdel+1,n,:)
            call searchdeltab(cc1,d1,dxx,cxx,c1)
            call searchdeltab(cc4,d4,dxx,cxx,c4)
            call searchdeltab(cc2,d2,dxx,c1,c2)
            call searchdeltab(cc3,d3,dxx,c4,c3)
        !else
        !  c1 = sx(idel-1,jdel+1,n,k)
        !  c2 = sx(idel  ,jdel+1,n,k)
        !  c3 = sx(idel+1,jdel+1,n,k)
        !  c4 = sx(idel+2,jdel+1,n,k)
        !end if
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(3) = c2+.5*xxg*(c3-c1 +xxg*(a3+xxg*a4))
        else
          r(3) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.) &
               -xxg*(1.+xxg)*c4/3.)                          &
               +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
        endif         !  (mhint==2)
        !if (intmode.eq.1) then
            cc1 = sx(idel  ,jdel-1,n,:)
            cc2 = sx(idel+1,jdel-1,n,:)
            cc3 = sx(idel  ,jdel+2,n,:)
            cc4 = sx(idel+1,jdel+2,n,:)
            d1 = dx(idel  ,jdel-1,n,:)
            d2 = dx(idel+1,jdel-1,n,:)
            d3 = dx(idel  ,jdel+2,n,:)
            d4 = dx(idel+1,jdel+2,n,:)
              call searchdeltab(cc1,d1,dxx,cxx,c1)
              call searchdeltab(cc4,d4,dxx,cxx,c4)
              call searchdeltab(cc2,d2,dxx,c1,c2)
              call searchdeltab(cc3,d3,dxx,c4,c3)
            r(1) = (1.-xxg)*c1 +xxg*c2
            r(4) = (1.-xxg)*c3 +xxg*c4
        !else
        !  do nn=1,4,3   ! N.B.
        !    c2 = sx(idel  ,jdel+nn-2,n,k)
        !    c3 = sx(idel+1,jdel+nn-2,n,k)
        !    r(nn) = (1.-xxg)*c2 +xxg*c3
        !  enddo         ! nn loop
        !end if
        if(mhint==2)then ! Bessel interp
          a4 = r(4)-r(1)+3.*(r(2)-r(3))
          a3 = r(1)-2.*r(2)+r(3)-a4
          sextra(iproc)%a(iq) = r(2) + 0.5*yyg*(r(3)-r(1) +yyg*(a3+yyg*a4))
        else
          sextra(iproc)%a(iq) = ((1.-yyg)*((2.-yyg)* &
            ((1.+yyg)*r(2)-yyg*r(1)/3.)              &
            -yyg*(1.+yyg)*r(4)/3.)                   &
            +yyg*(1.+yyg)*(2.-yyg)*r(3))/2.
        endif         !  (mhint==2)
       endif           ! if (ddin(iproc,iq).gt.0.1) then
      enddo            ! iq loop
    end do              ! iproc loop
  else                   ! (nfield<mh_bs)
    do iproc=0,nproc-1
      if ( iproc == myid ) then
        cycle
      end if
      do iq=1,drlen(iproc)
       if (ddin(iproc,iq).gt.0.1) then
        !  Convert face index from 0:npanels to array indices
        ip = min(il_g,max(1,nint(dpoints(iproc)%a(2,iq))))
        jp = min(il_g,max(1,nint(dpoints(iproc)%a(3,iq))))
        n = nint(dpoints(iproc)%a(1,iq)) + noff ! Local index
        !  Need global face index in fproc call
        idel = int(dpoints(iproc)%a(2,iq))
        xxg = dpoints(iproc)%a(2,iq) - idel
        jdel = int(dpoints(iproc)%a(3,iq))
        yyg = dpoints(iproc)%a(3,iq) - jdel
        k = nint(dpoints(iproc)%a(4,iq))
        idel = idel - ioff(n-noff)
        jdel = jdel - joff(n-noff)
        !if (intmode.eq.1) then
          cc1 = sx(idel-1,jdel,n,:)
          cc2 = sx(idel  ,jdel,n,:)
          cc3 = sx(idel+1,jdel,n,:)
          cc4 = sx(idel+2,jdel,n,:)
          d1 = dx(idel-1,jdel,n,:)
          d2 = dx(idel  ,jdel,n,:)
          d3 = dx(idel+1,jdel,n,:)
          d4 = dx(idel+2,jdel,n,:)
          dxx = d(iq,k)
            cxx = s(iq,k)
            call searchdeltab(cc1,d1,dxx,cxx,c1)
            call searchdeltab(cc4,d4,dxx,cxx,c4)
            call searchdeltab(cc2,d2,dxx,c1,c2)
            call searchdeltab(cc3,d3,dxx,c4,c3)
        !else
        !  c1 = sx(idel-1,jdel,n,k)
        !  c2 = sx(idel  ,jdel,n,k)
        !  c3 = sx(idel+1,jdel,n,k)
        !  c4 = sx(idel+2,jdel,n,k)
        !end if
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
        !if (intmode.eq.1) then
          cc1 = sx(idel-1,jdel+1,n,:)
          cc2 = sx(idel  ,jdel+1,n,:)
          cc3 = sx(idel+1,jdel+1,n,:)
          cc4 = sx(idel+2,jdel+1,n,:)
          d1 = dx(idel-1,jdel+1,n,:)
          d2 = dx(idel  ,jdel+1,n,:)
          d3 = dx(idel+1,jdel+1,n,:)
          d4 = dx(idel+2,jdel+1,n,:)
            call searchdeltab(cc1,d1,dxx,cxx,c1)
            call searchdeltab(cc4,d4,dxx,cxx,c4)
            call searchdeltab(cc2,d2,dxx,c1,c2)
            call searchdeltab(cc3,d3,dxx,c4,c3)
        !else
        !  c1 = sx(idel-1,jdel+1,n,k)
        !  c2 = sx(idel  ,jdel+1,n,k)
        !  c3 = sx(idel+1,jdel+1,n,k)
        !  c4 = sx(idel+2,jdel+1,n,k)
        !end if
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
        !if (intmode.eq.1) then
            cc1 = sx(idel  ,jdel-1,n,:)
            cc2 = sx(idel+1,jdel-1,n,:)
            cc3 = sx(idel  ,jdel+2,n,:)
            cc4 = sx(idel+1,jdel+2,n,:)
            d1 = dx(idel  ,jdel-1,n,:)
            d2 = dx(idel+1,jdel-1,n,:)
            d3 = dx(idel  ,jdel+2,n,:)
            d4 = dx(idel+1,jdel+2,n,:)
              call searchdeltab(cc1,d1,dxx,cxx,c1)
              call searchdeltab(cc4,d4,dxx,cxx,c4)
              call searchdeltab(cc2,d2,dxx,c1,c2)
              call searchdeltab(cc3,d3,dxx,c4,c3)
            r(1) = (1.-xxg)*c1 +xxg*c2
            r(4) = (1.-xxg)*c3 +xxg*c4
        !else
        !  do nn=1,4,3   ! N.B.
        !    c2 = sx(idel  ,jdel+nn-2,n,k)
        !    c3 = sx(idel+1,jdel+nn-2,n,k)
        !    r(nn) = (1.-xxg)*c2 +xxg*c3
        !  enddo         ! nn loop
        !end if
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
       endif            ! if (ddin(iproc,iq).gt.0.1) then
      enddo            ! iq loop
    end do              ! iproc loop
  endif                  ! (nfield<mh_bs)  .. else ..
!========================   end of intsch=1 section ====================
else     ! if(intsch==1)then
!======================== start of intsch=2 section ====================
  
!       this is intsc           NS interps done first
!       first extend s arrays into sx - this one -1:il+2 & -1:il+2

    do n=1,npan         ! first simple copy into larger array
      do j=1,jpan
        do i=1,ipan
          sx(i,j,n,:)=s(ind(i,j,n),:)
          dx(i,j,n,:)=d(ind(i,j,n),:)
        enddo         ! i loop
      enddo            ! j loop
    enddo               ! n loop

    do n=1,npan
      do j=1,jpan
        sx(0,j,n,:) = s(iw(ind(1,j,n)),:)
        sx(-1,j,n,:) = s(iww(ind(1,j,n)),:)
        sx(ipan+1,j,n,:) = s(ie(ind(ipan,j,n)),:)
        sx(ipan+2,j,n,:) = s(iee(ind(ipan,j,n)),:)
        dx(0,j,n,:) = d(iw(ind(1,j,n)),:)
        dx(-1,j,n,:) = d(iww(ind(1,j,n)),:)
        dx(ipan+1,j,n,:) = d(ie(ind(ipan,j,n)),:)
        dx(ipan+2,j,n,:) = d(iee(ind(ipan,j,n)),:)
      enddo            ! j loop
      do i=1,ipan
        sx(i,0,n,:) = s(is(ind(i,1,n)),:)
        sx(i,-1,n,:) = s(iss(ind(i,1,n)),:)
        sx(i,jpan+1,n,:) = s(in(ind(i,jpan,n)),:)
        sx(i,jpan+2,n,:) = s(inn(ind(i,jpan,n)),:)
        dx(i,0,n,:) = d(is(ind(i,1,n)),:)
        dx(i,-1,n,:) = d(iss(ind(i,1,n)),:)
        dx(i,jpan+1,n,:) = d(in(ind(i,jpan,n)),:)
        dx(i,jpan+2,n,:) = d(inn(ind(i,jpan,n)),:)
      enddo            ! i loop
!        for ns interpolation, sometimes need (different from ew):
!            (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
!          (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)
      sx(-1,0,n,:)=s(lsww(n),:)
      sx(0,0,n,:) = s(lsw(n),:)
      sx(0,-1,n,:) = s(lssw(n),:)
      sx(ipan+2,0,n,:) = s(lsee(n),:)
      sx(ipan+1,-1,n,:) = s(lsse(n),:)
      sx(-1,jpan+1,n,:) = s(lnww(n),:)
      sx(0,jpan+1,n,:) = s(lnw(n),:)
      sx(0,jpan+2,n,:) = s(lnnw(n),:)
      sx(ipan+2,jpan+1,n,:) = s(lnee(n),:)
      sx(ipan+1,jpan+2,n,:) = s(lnne(n),:)
      sx(ipan+1,0,n,:)    = s(ise(ind(ipan,1,n)),:)
      sx(ipan+1,jpan+1,n,:) = s(ine(ind(ipan,jpan,n)),:)
      dx(-1,0,n,:)=d(lsww(n),:)
      dx(0,0,n,:) = d(lsw(n),:)
      dx(0,-1,n,:) = d(lssw(n),:)
      dx(ipan+2,0,n,:) = d(lsee(n),:)
      dx(ipan+1,-1,n,:) = d(lsse(n),:)
      dx(-1,jpan+1,n,:) = d(lnww(n),:)
      dx(0,jpan+1,n,:) = d(lnw(n),:)
      dx(0,jpan+2,n,:) = d(lnnw(n),:)
      dx(ipan+2,jpan+1,n,:) = d(lnee(n),:)
      dx(ipan+1,jpan+2,n,:) = d(lnne(n),:)
      dx(ipan+1,0,n,:)    = d(ise(ind(ipan,1,n)),:)
      dx(ipan+1,jpan+1,n,:) = d(ine(ind(ipan,jpan,n)),:)
    enddo               ! n loop
  
    if(nfield<mh_bs)then
     do iq=1,ifull    ! non Berm-Stan option
      if (wtr(iq)) then
       do k=1,kx
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
        !if (intmode.eq.1) then
          cc1 = sx(idel,jdel-1,n,:)
          cc2 = sx(idel,jdel  ,n,:)
          cc3 = sx(idel,jdel+1,n,:)
          cc4 = sx(idel,jdel+2,n,:)
          d1 = dx(idel,jdel-1,n,:)
          d2 = dx(idel,jdel  ,n,:)
          d3 = dx(idel,jdel+1,n,:)
          d4 = dx(idel,jdel+2,n,:)
          dxx = d(iq,k)
            cxx = s(iq,k)
            call searchdeltab(cc1,d1,dxx,cxx,c1)
            call searchdeltab(cc4,d4,dxx,cxx,c4)
            call searchdeltab(cc2,d2,dxx,c1,c2)
            call searchdeltab(cc3,d3,dxx,c4,c3)
        !else
        !  c1 = sx(idel,jdel-1,n,k)
        !  c2 = sx(idel,jdel  ,n,k)
        !  c3 = sx(idel,jdel+1,n,k)
        !  c4 = sx(idel,jdel+2,n,k)
        !end if
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(2) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
        else
          r(2) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.) &
               -yyg*(1.+yyg)*c4/3.)                          &
               +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
        endif         !  (mhint==2)
        !if (intmode.eq.1) then
          cc1 = sx(idel+1,jdel-1,n,:)
          cc2 = sx(idel+1,jdel  ,n,:)
          cc3 = sx(idel+1,jdel+1,n,:)
          cc4 = sx(idel+1,jdel+2,n,:)
          d1 = dx(idel+1,jdel-1,n,:)
          d2 = dx(idel+1,jdel  ,n,:)
          d3 = dx(idel+1,jdel+1,n,:)
          d4 = dx(idel+1,jdel+2,n,:)
            call searchdeltab(cc1,d1,dxx,cxx,c1)
            call searchdeltab(cc4,d4,dxx,cxx,c4)
            call searchdeltab(cc2,d2,dxx,c1,c2)
            call searchdeltab(cc3,d3,dxx,c4,c3)
        !else
        !  c1 = sx(idel+1,jdel-1,n,k)
        !  c2 = sx(idel+1,jdel  ,n,k)
        !  c3 = sx(idel+1,jdel+1,n,k)
        !  c4 = sx(idel+1,jdel+2,n,k)
        !end if
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(3) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
        else
          r(3) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.) &
               -yyg*(1.+yyg)*c4/3.)                          &
               +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
        endif         !  (mhint==2)
        !if (intmode.eq.1) then
            cc1 = sx(idel-1,jdel  ,n,:)
            cc2 = sx(idel-1,jdel+1,n,:)
            cc3 = sx(idel+2,jdel  ,n,:)
            cc4 = sx(idel+2,jdel+1,n,:)
            d1 = dx(idel-1,jdel  ,n,:)
            d2 = dx(idel-1,jdel+1,n,:)
            d3 = dx(idel+2,jdel  ,n,:)
            d4 = dx(idel+2,jdel+1,n,:)
              call searchdeltab(cc1,d1,dxx,cxx,c1)
              call searchdeltab(cc4,d4,dxx,cxx,c4)
              call searchdeltab(cc2,d2,dxx,c1,c2)
              call searchdeltab(cc3,d3,dxx,c4,c3)
            r(1) = (1.-yyg)*c1 +yyg*c2
            r(4) = (1.-yyg)*c3 +yyg*c4
        !else
        !  do nn=1,4,3   ! N.B.
        !    c2 = sx(idel+nn-2,jdel  ,n,k)
        !    c3 = sx(idel+nn-2,jdel+1,n,k)
        !    r(nn) = (1.-yyg)*c2 +yyg*c3
        !  enddo         ! nn loop
        !end if
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
       enddo
      endif           ! if (wtr) then
     enddo            ! iq loop
    else                ! (nfield<mh_bs)
     do iq=1,ifull    ! Berm-Stan option here e.g. qg & gases
      if (wtr(iq)) then
       do k=1,kx
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
        !if (intmode.eq.1) then
          cc1 = sx(idel,jdel-1,n,:)
          cc2 = sx(idel,jdel  ,n,:)
          cc3 = sx(idel,jdel+1,n,:)
          cc4 = sx(idel,jdel+2,n,:)
          d1 = dx(idel,jdel-1,n,:)
          d2 = dx(idel,jdel  ,n,:)
          d3 = dx(idel,jdel+1,n,:)
          d4 = dx(idel,jdel+2,n,:)
          dxx = d(iq,k)
            cxx = s(iq,k)
            call searchdeltab(cc1,d1,dxx,cxx,c1)
            call searchdeltab(cc4,d4,dxx,cxx,c4)
            call searchdeltab(cc2,d2,dxx,c1,c2)
            call searchdeltab(cc3,d3,dxx,c4,c3)
        !else
        !  c1 = sx(idel,jdel-1,n,k)
        !  c2 = sx(idel,jdel  ,n,k)
        !  c3 = sx(idel,jdel+1,n,k)
        !  c4 = sx(idel,jdel+2,n,k)
        !end if
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
        !if (intmode.eq.1) then
          cc1 = sx(idel+1,jdel-1,n,:)
          cc2 = sx(idel+1,jdel  ,n,:)
          cc3 = sx(idel+1,jdel+1,n,:)
          cc4 = sx(idel+1,jdel+2,n,:)
          d1 = dx(idel+1,jdel-1,n,:)
          d2 = dx(idel+1,jdel  ,n,:)
          d3 = dx(idel+1,jdel+1,n,:)
          d4 = dx(idel+1,jdel+2,n,:)
            call searchdeltab(cc1,d1,dxx,cxx,c1)
            call searchdeltab(cc4,d4,dxx,cxx,c4)
            call searchdeltab(cc2,d2,dxx,c1,c2)
            call searchdeltab(cc3,d3,dxx,c4,c3)
        !else
        !  c1 = sx(idel+1,jdel-1,n,k)
        !  c2 = sx(idel+1,jdel  ,n,k)
        !  c3 = sx(idel+1,jdel+1,n,k)
        !  c4 = sx(idel+1,jdel+2,n,k)
        !end if
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
        !if (intmode.eq.1) then
            cc1 = sx(idel-1,jdel  ,n,:)
            cc2 = sx(idel-1,jdel+1,n,:)
            cc3 = sx(idel+2,jdel  ,n,:)
            cc4 = sx(idel+2,jdel+1,n,:)            
            d1 = dx(idel-1,jdel  ,n,:)
            d2 = dx(idel-1,jdel+1,n,:)
            d3 = dx(idel+2,jdel  ,n,:)
            d4 = dx(idel+2,jdel+1,n,:)
              call searchdeltab(cc1,d1,dxx,cxx,c1)
              call searchdeltab(cc4,d4,dxx,cxx,c4)
              call searchdeltab(cc2,d2,dxx,c1,c2)
              call searchdeltab(cc3,d3,dxx,c4,c3)
            r(1) = (1.-yyg)*c1 +yyg*c2
            r(4) = (1.-yyg)*c3 +yyg*c4
        !else
        !  do nn=1,4,3   ! N.B.
        !    c2 = sx(idel+nn-2,jdel  ,n,k)
        !    c3 = sx(idel+nn-2,jdel+1,n,k)
        !    r(nn) = (1.-yyg)*c2 +yyg*c3
        !  enddo         ! nn loop
        !end if
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
       end do
      endif           ! if (wtr) then
     enddo            ! iq loop
    endif               ! (nfield<mh_bs)  .. else ..

! For other processes

  !if (intmode.eq.1) then
    call intssend(d)
    do iproc=0,nproc-1
      ddin(iproc,1:drlen(iproc))=sextra(iproc)%a(1:drlen(iproc))
    end do
    call intssend(s)
    do iproc=0,nproc-1
      sdin(iproc,1:drlen(iproc))=sextra(iproc)%a(1:drlen(iproc))
    end do
  !end if

  if(nfield<mh_bs)then
    do iproc=0,nproc-1
      if ( iproc == myid ) then
        cycle
      end if
      do iq=1,drlen(iproc)
       if (ddin(iproc,iq).gt.0.1) then
        !  Convert face index from 0:npanels to array indices
        ip = min(il_g,max(1,nint(dpoints(iproc)%a(2,iq))))
        jp = min(il_g,max(1,nint(dpoints(iproc)%a(3,iq))))
        n = nint(dpoints(iproc)%a(1,iq)) + noff ! Local index
        !  Need global face index in fproc call
        idel = int(dpoints(iproc)%a(2,iq))
        xxg = dpoints(iproc)%a(2,iq) - idel
        jdel = int(dpoints(iproc)%a(3,iq))
        yyg = dpoints(iproc)%a(3,iq) - jdel
        k = nint(dpoints(iproc)%a(4,iq))
        idel = idel - ioff(n-noff)
        jdel = jdel - joff(n-noff)
        !if (intmode.eq.1) then
          cc1 = sx(idel,jdel-1,n,:)
          cc2 = sx(idel,jdel  ,n,:)
          cc3 = sx(idel,jdel+1,n,:)
          cc4 = sx(idel,jdel+2,n,:)
          d1 = dx(idel,jdel-1,n,:)
          d2 = dx(idel,jdel  ,n,:)
          d3 = dx(idel,jdel+1,n,:)
          d4 = dx(idel,jdel+2,n,:)
          dxx = d(iq,k)
            cxx = s(iq,k)
            call searchdeltab(cc1,d1,dxx,cxx,c1)
            call searchdeltab(cc4,d4,dxx,cxx,c4)
            call searchdeltab(cc2,d2,dxx,c1,c2)
            call searchdeltab(cc3,d3,dxx,c4,c3)
        !else
        !  c1 = sx(idel,jdel-1,n,k)
        !  c2 = sx(idel,jdel  ,n,k)
        !  c3 = sx(idel,jdel+1,n,k)
        !  c4 = sx(idel,jdel+2,n,k)
        !end if
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(2) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
        else
          r(2) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.) &
               -yyg*(1.+yyg)*c4/3.)                          &
               +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
        endif         !  (mhint==2)
        !if (intmode.eq.1) then
          cc1 = sx(idel+1,jdel-1,n,:)
          cc2 = sx(idel+1,jdel  ,n,:)
          cc3 = sx(idel+1,jdel+1,n,:)
          cc4 = sx(idel+1,jdel+2,n,:)
          d1 = dx(idel+1,jdel-1,n,:)
          d2 = dx(idel+1,jdel  ,n,:)
          d3 = dx(idel+1,jdel+1,n,:)
          d4 = dx(idel+1,jdel+2,n,:)
            call searchdeltab(cc1,d1,dxx,cxx,c1)
            call searchdeltab(cc4,d4,dxx,cxx,c4)
            call searchdeltab(cc2,d2,dxx,c1,c2)
            call searchdeltab(cc3,d3,dxx,c4,c3)
        !else
        !  c1 = sx(idel+1,jdel-1,n,k)
        !  c2 = sx(idel+1,jdel  ,n,k)
        !  c3 = sx(idel+1,jdel+1,n,k)
        !  c4 = sx(idel+1,jdel+2,n,k)
        !end if
        if(mhint==2)then ! Bessel interp
          a4 = c4-c1+3.*(c2-c3)
          a3 = c1-2.*c2+c3-a4
          r(3) = c2+.5*yyg*(c3-c1 +yyg*(a3+yyg*a4))
        else
          r(3) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*c2-yyg*c1/3.) &
               -yyg*(1.+yyg)*c4/3.)                          &
               +yyg*(1.+yyg)*(2.-yyg)*c3)/2.
        endif         !  (mhint==2)
        !if (intmode.eq.1) then
            cc1 = sx(idel-1,jdel  ,n,:)
            cc2 = sx(idel-1,jdel+1,n,:)
            cc3 = sx(idel+2,jdel  ,n,:)
            cc4 = sx(idel+2,jdel+1,n,:)
            d1 = dx(idel-1,jdel  ,n,:)
            d2 = dx(idel-1,jdel+1,n,:)
            d3 = dx(idel+2,jdel  ,n,:)
            d4 = dx(idel+2,jdel+1,n,:)
              call searchdeltab(cc1,d1,dxx,cxx,c1)
              call searchdeltab(cc4,d4,dxx,cxx,c4)
              call searchdeltab(cc2,d2,dxx,c1,c2)
              call searchdeltab(cc3,d3,dxx,c4,c3)
            r(1) = (1.-yyg)*c1 +yyg*c2
            r(4) = (1.-yyg)*c3 +yyg*c4
        !else
        !  do nn=1,4,3   ! N.B.
        !    c2 = sx(idel+nn-2,jdel  ,n,k)
        !    c3 = sx(idel+nn-2,jdel+1,n,k)
        !    r(nn) = (1.-yyg)*c2 +yyg*c3
        !  enddo         ! nn loop
        !end if
        if(mhint==2)then ! Bessel interp
          a4 = r(4)-r(1)+3.*(r(2)-r(3))
          a3 = r(1)-2.*r(2)+r(3)-a4
          sextra(iproc)%a(iq) = r(2)+ 0.5*xxg*(r(3)-r(1) +xxg*(a3+xxg*a4))
        else
          sextra(iproc)%a(iq) = ((1.-xxg)*((2.-xxg)*  &
            ((1.+xxg)*r(2)-xxg*r(1)/3.)               &
            -xxg*(1.+xxg)*r(4)/3.)                    &
            +xxg*(1.+xxg)*(2.-xxg)*r(3))/2.
        endif         !  (mhint==2)
       endif           ! if (ddin(iproc,iq).gt.0.1) then
      enddo            ! iq loop
    end do              ! iproc
  else                   ! (nfield<mh_bs)
    do iproc=0,nproc-1
      if ( iproc == myid ) then
        cycle
      end if
      do iq=1,drlen(iproc)
       if (ddin(iproc,iq).gt.0.1) then
        !  Convert face index from 0:npanels to array indices
        ip = min(il_g,max(1,nint(dpoints(iproc)%a(2,iq))))
        jp = min(il_g,max(1,nint(dpoints(iproc)%a(3,iq))))
        n = nint(dpoints(iproc)%a(1,iq)) + noff ! Local index
        !  Need global face index in fproc call
        idel = int(dpoints(iproc)%a(2,iq))
        xxg = dpoints(iproc)%a(2,iq) - idel
        jdel = int(dpoints(iproc)%a(3,iq))
        yyg = dpoints(iproc)%a(3,iq) - jdel
        k = nint(dpoints(iproc)%a(4,iq))
        idel = idel - ioff(n-noff)
        jdel = jdel - joff(n-noff)
        !if (intmode.eq.1) then
          cc1 = sx(idel,jdel-1,n,:)
          cc2 = sx(idel,jdel  ,n,:)
          cc3 = sx(idel,jdel+1,n,:)
          cc4 = sx(idel,jdel+2,n,:)
          d1 = dx(idel,jdel-1,n,:)
          d2 = dx(idel,jdel  ,n,:)
          d3 = dx(idel,jdel+1,n,:)
          d4 = dx(idel,jdel+2,n,:)
          dxx = d(iq,k)
            cxx = s(iq,k)
            call searchdeltab(cc1,d1,dxx,cxx,c1)
            call searchdeltab(cc4,d4,dxx,cxx,c4)
            call searchdeltab(cc2,d2,dxx,c1,c2)
            call searchdeltab(cc3,d3,dxx,c4,c3)
        !else
        !  c1 = sx(idel,jdel-1,n,k)
        !  c2 = sx(idel,jdel  ,n,k)
        !  c3 = sx(idel,jdel+1,n,k)
        !  c4 = sx(idel,jdel+2,n,k)
        !end if
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
        !if (intmode.eq.1) then
          cc1 = sx(idel+1,jdel-1,n,:)
          cc2 = sx(idel+1,jdel  ,n,:)
          cc3 = sx(idel+1,jdel+1,n,:)
          cc4 = sx(idel+1,jdel+2,n,:)
          d1 = dx(idel+1,jdel-1,n,:)
          d2 = dx(idel+1,jdel  ,n,:)
          d3 = dx(idel+1,jdel+1,n,:)
          d4 = dx(idel+1,jdel+2,n,:)
            call searchdeltab(cc1,d1,dxx,cxx,c1)
            call searchdeltab(cc4,d4,dxx,cxx,c4)
            call searchdeltab(cc2,d2,dxx,c1,c2)
            call searchdeltab(cc3,d3,dxx,c4,c3)
        !else
        !  c1 = sx(idel+1,jdel-1,n,k)
        !  c2 = sx(idel+1,jdel  ,n,k)
        !  c3 = sx(idel+1,jdel+1,n,k)
        !  c4 = sx(idel+1,jdel+2,n,k)
        !end if
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
        !if (intmode.eq.1) then
            cc1 = sx(idel-1,jdel  ,n,:)
            cc2 = sx(idel-1,jdel+1,n,:)
            cc3 = sx(idel+2,jdel  ,n,:)
            cc4 = sx(idel+2,jdel+1,n,:)
            d1 = dx(idel-1,jdel  ,n,:)
            d2 = dx(idel-1,jdel+1,n,:)
            d3 = dx(idel+2,jdel  ,n,:)
            d4 = dx(idel+2,jdel+1,n,:)
              call searchdeltab(cc1,d1,dxx,cxx,c1)
              call searchdeltab(cc4,d4,dxx,cxx,c4)
              call searchdeltab(cc2,d2,dxx,c1,c2)
              call searchdeltab(cc3,d3,dxx,c4,c3)
            r(1) = (1.-yyg)*c1 +yyg*c2
            r(4) = (1.-yyg)*c3 +yyg*c4
        !else
        !  do nn=1,4,3   ! N.B.
        !    c2 = sx(idel+nn-2,jdel  ,n,k)
        !    c3 = sx(idel+nn-2,jdel+1,n,k)
        !    r(nn) = (1.-yyg)*c2 +yyg*c3
        !  enddo         ! nn loop
        !end if
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
       endif            ! if (ddin(iproc,iq).gt.0.1) then
      enddo            ! iq loop
    end do              ! iproc
  endif                  ! (nfield<mh_bs)  .. else ..

endif                     ! (intsch==1) .. else ..
!========================   end of intsch=1 section ====================

call intssync(s)

return
end subroutine mlointsb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Rotate wind vector to arrival point

subroutine mlorot(cou,cov,cow,x3d,y3d,z3d)

use mlo
use xyzinfo_m

implicit none

include 'newmpar.h'

integer k,kx
real, dimension(:,:), intent(inout) :: cou,cov,cow
real, dimension(ifull) :: vec1x,vec1y,vec1z,denb
real, dimension(ifull) :: vec2x,vec2y,vec2z,vecdot
real, dimension(ifull) :: vec3x,vec3y,vec3z,vdot1,vdot2
real*8, dimension(ifull,size(cou,2)), intent(in) :: x3d,y3d,z3d

kx=size(cou,2)

do k=1,kx
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

integer k,itn,kx
real, dimension(:,:), intent(in) :: u
real, dimension(ifull,size(u,2)), intent(in) :: v
real, dimension(ifull,size(u,2)), intent(out) :: uout,vout
real, dimension(ifull+iextra,size(u,2)) :: uin,vin
real, dimension(ifull+iextra,size(u,2)) :: ua,va,ud,vd
integer, parameter :: itnmax=3

kx=size(u,2)

uin(1:ifull,:)=u
vin(1:ifull,:)=v
call boundsuv(uin,vin)

do k=1,kx
  ud(1:ifull,k)= uin(iwu,k)*0.1+uin(1:ifull,k)+uin(ieu,k)*0.5
  vd(1:ifull,k)= vin(isv,k)*0.1+vin(1:ifull,k)+vin(inv,k)*0.5
enddo
call boundsuv(ud,vd)

do k=1,kx
  ua(1:ifull,k)=ud(1:ifull,k)-ud(iwu,k)*0.5 ! 1st guess
  va(1:ifull,k)=vd(1:ifull,k)-vd(isv,k)*0.5 ! 1st guess
enddo

do itn=1,itnmax        ! each loop is a double iteration
  call boundsuv(ua,va,nrows=2)
  do k=1,kx
    uin(1:ifull,k)=(ud(1:ifull,k)-0.5*ud(iwu,k)-ua(ieu,k)*0.1 +ua(iwwu,k)*0.25)/.95
    vin(1:ifull,k)=(vd(1:ifull,k)-0.5*vd(isv,k)-va(inv,k)*0.1 +va(issv,k)*0.25)/.95
  enddo
  call boundsuv(uin,vin,nrows=2)
  do k=1,kx
    ua(1:ifull,k)=(ud(1:ifull,k)-0.5*ud(iwu,k)-uin(ieu,k)*0.1 +uin(iwwu,k)*0.25)/.95
    va(1:ifull,k)=(vd(1:ifull,k)-0.5*vd(isv,k)-vin(inv,k)*0.1 +vin(issv,k)*0.25)/.95
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

integer k,itn,kx
real, dimension(:,:), intent(in) :: u
real, dimension(ifull,size(u,2)), intent(in) :: v
real, dimension(ifull,size(u,2)), intent(out) :: uout,vout
real, dimension(ifull+iextra,size(u,2)) :: uin,vin
real, dimension(ifull+iextra,size(u,2)) :: ua,va,ud,vd
integer, parameter :: itnmax=3

kx=size(u,2)

uin(1:ifull,:)=u
vin(1:ifull,:)=v
call boundsuv(uin,vin)

do k=1,kx
  ud(1:ifull,k)= uin(ieu,k)*0.1+uin(1:ifull,k)+uin(iwu,k)*0.5
  vd(1:ifull,k)= vin(inv,k)*0.1+vin(1:ifull,k)+vin(isv,k)*0.5
enddo
call boundsuv(ud,vd)

do k=1,kx
  ua(1:ifull,k)=ud(1:ifull,k)-ud(ieu,k)*0.5 ! 1st guess
  va(1:ifull,k)=vd(1:ifull,k)-vd(inv,k)*0.5 ! 1st guess
enddo

do itn=1,itnmax        ! each loop is a double iteration
  call boundsuv(ua,va,nrows=2)
  do k=1,kx
    uin(1:ifull,k)=(ud(1:ifull,k)-0.5*ud(ieu,k)-ua(iwu,k)*0.1 +ua(ieeu,k)*0.25)/.95
    vin(1:ifull,k)=(vd(1:ifull,k)-0.5*vd(inv,k)-va(isv,k)*0.1 +va(innv,k)*0.25)/.95
  enddo
  call boundsuv(uin,vin,nrows=2)
  do k=1,kx
    ua(1:ifull,k)=(ud(1:ifull,k)-0.5*ud(ieu,k)-uin(iwu,k)*0.1 +uin(ieeu,k)*0.25)/.95
    va(1:ifull,k)=(vd(1:ifull,k)-0.5*vd(inv,k)-vin(isv,k)*0.1 +vin(innv,k)*0.25)/.95
  enddo
enddo                  ! itn=1,itnmax
      
uout=ua(1:ifull,:)
vout=va(1:ifull,:)

return
end subroutine mlounstaguv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine estimates vertical velocity

subroutine getww(cou,cov,dd,dz,ee,nw)

use indices_m
use map_m
use mlo

implicit none

include 'newmpar.h'
include 'parm.h'

integer ii
real, dimension(ifull,wlev), intent(out) :: nw
real, dimension(ifull,wlev), intent(in) :: dz
real, dimension(ifull+iextra,wlev), intent(in) :: cou,cov
real, dimension(ifull,wlev) :: kku,kkv
real, dimension(ifull), intent(in) :: dd,ee
real, dimension(ifull) :: div

div=0.
do ii=1,wlev
  kku(:,ii)=(cou(1:ifull,ii)/emu(1:ifull)-cou(iwu,ii)/emu(iwu))*em(1:ifull)*em(1:ifull)/ds
  kkv(:,ii)=(cov(1:ifull,ii)/emv(1:ifull)-cov(isv,ii)/emv(isv))*em(1:ifull)*em(1:ifull)/ds
  div=div+(kku(:,ii)+kkv(:,ii))*dz(:,ii)
end do
div=div/dd(1:ifull)
! nw is at half levels with nw(:,1) at the bottom of level 1
! positive nw is moving downwards to the ocean bottom
nw(:,1)=(div-kku(:,1)-kkv(:,1))*dz(:,1)
do ii=2,wlev
  nw(:,ii)=nw(:,ii-1)+(div-kku(:,ii)-kkv(:,ii))*dz(:,ii)
end do
do ii=1,wlev
  nw(:,ii)=nw(:,ii)*ee
end do

return
end subroutine getww

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine performs vertical advection

subroutine mlovadv(dtin,ww,uu,vv,ss,tt,depdum,dzdum,wtr)

use mlo

implicit none

include 'mpif.h'
include 'newmpar.h'

integer its,its_g,ii,l,iq,ierr
real, intent(in) :: dtin
real dtnew
real, dimension(ifull,wlev), intent(in) :: ww,depdum,dzdum
real, dimension(ifull,wlev), intent(inout) :: uu,vv,ss,tt
logical, dimension(ifull), intent(in) :: wtr

! reduce time step to ensure stability
dtnew=dtin
do iq=1,ifull
  if (wtr(iq)) then
    dtnew=min(dtnew,0.1*dzdum(iq,1)/max(abs(ww(iq,1)),1.E-8))
    do ii=2,wlev-1
      dtnew=min(dtnew,0.1*dzdum(iq,ii)/max(abs(ww(iq,ii)),abs(ww(iq,ii-1)),1.E-8))
    end do
    dtnew=min(dtnew,0.1*dzdum(iq,wlev)/max(abs(ww(iq,wlev-1)),1.E-8))
  end if
end do
its=int(dtin/dtnew)+1
call MPI_AllReduce(its, its_g, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr )
dtnew=dtin/real(its_g)

!if (its_g.gt.10) write(6,*) "MLOVERT ",its_g

do l=1,its_g
  call mlotvd(dtnew,ww,uu,depdum,dzdum)
  call mlotvd(dtnew,ww,vv,depdum,dzdum)
  call mlotvd(dtnew,ww,ss,depdum,dzdum)
  call mlotvd(dtnew,ww,tt,depdum,dzdum)
end do

return
end subroutine mlovadv

subroutine mlotvd(dtnew,ww,uu,dep,dz)

use mlo

implicit none

include 'newmpar.h'

integer ii
real, intent(in) :: dtnew
real, dimension(ifull,wlev), intent(in) :: ww,dep,dz
real, dimension(ifull,wlev), intent(inout) :: uu
real, dimension(ifull,wlev-1) :: ff
real, dimension(ifull,0:wlev) :: delu
real, dimension(ifull) :: fl,fh,cc,rr,xx,jj

delu=0.
do ii=1,wlev-1
  delu(:,ii)=uu(:,ii+1)-uu(:,ii)
end do

! TVD part
do ii=1,wlev-1
  fl=0.5*ww(:,ii)*(uu(:,ii)+uu(:,ii+1))+0.5*abs(ww(:,ii))*(uu(:,ii)-uu(:,ii+1))
  fh=ww(:,ii)*0.5*(uu(:,ii)+uu(:,ii+1)) &
     -0.5*(uu(:,ii+1)-uu(:,ii))*ww(:,ii)**2*dtnew/max(dep(:,ii+1)-dep(:,ii),1.E-10)
  xx=delu(:,ii)+sign(1.E-20,delu(:,ii))
  jj=sign(1.,ww(:,ii))
  rr=0.5*(jj*(delu(:,ii+1)-delu(:,ii-1))+abs(jj)*(delu(:,ii+1)+delu(:,ii-1)))/xx
  cc=max(0.,min(1.,2.*rr),min(2.,rr)) ! superbee
  ff(:,ii)=fl+cc*(fh-fl)
  !ff(:,ii)=ww(:,ii)*0.5*(uu(:,ii)+uu(:,ii+1)) ! explicit
end do
uu(:,1)=uu(:,1)+dtnew*(uu(:,1)*ww(:,1)-ff(:,1))/dz(:,1)
do ii=2,wlev-1
  uu(:,ii)=uu(:,ii)+dtnew*(uu(:,ii)*(ww(:,ii)-ww(:,ii-1))-(ff(:,ii)-ff(:,ii-1)))/dz(:,ii)
end do
uu(:,wlev)=uu(:,wlev)+dtnew*(-uu(:,wlev)*ww(:,wlev-1)+ff(:,wlev-1))/dz(:,wlev)

return
end subroutine mlotvd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine interpolates to a common horizontal surface

subroutine stagtruedelta(s,dep,wtr,dsdxu,dsdyu,dsdxv,dsdyv)

use map_m
use mlo
use indices_m

implicit none

include 'newmpar.h'
include 'parm.h'

integer iq,ii
integer idn,ide,ids,idw,idne,idse,iden,idwn
real sx,sn,se,ss,sw,sne,sse,sen,swn,ddx
real tnu,tsu,tev,twv
real, dimension(wlev) :: ddn,dde,dds,ddw,ddne,ddse,dden,ddwn
real, dimension(wlev) :: isin,isie,isis,isiw,isine,isise,isien,isiwn
real, dimension(ifull+iextra,wlev), intent(in) :: s,dep
real, dimension(ifull,wlev), intent(out) :: dsdxu,dsdyu,dsdxv,dsdyv
logical, dimension(ifull+iextra), intent(in) :: wtr

dsdxu=0.
dsdyu=0.
dsdxv=0.
dsdyv=0.

do iq=1,ifull
  if (wtr(iq)) then

    idn=1
    ide=1
    ids=1
    idw=1
    idne=1
    idse=1
    iden=1
    idwn=1

    ddn  =dep(in(iq),:)
    isin =s(in(iq),:)
    ddne =dep(ine(iq),:)
    isine=s(ine(iq),:)
    dds  =dep(is(iq),:)
    isis =s(is(iq),:)
    ddse =dep(ise(iq),:)
    isise=s(ise(iq),:)
    dde  =dep(ie(iq),:)
    isie =s(ie(iq),:)
    dden =dep(ien(iq),:)
    isien=s(ien(iq),:)
    ddw  =dep(iw(iq),:)
    isiw =s(iw(iq),:)
    ddwn =dep(iwn(iq),:)
    isiwn=s(iwn(iq),:)
    
    ! now search for other levels
    do ii=1,wlev
      sx=s(iq,ii)
      ddx=dep(iq,ii)
    
      call searchdelta(isin ,ddn ,ddx,sx,idn ,sn)
      call searchdelta(isine,ddne,ddx,sn,idne,sne)
      call searchdelta(isis ,dds ,ddx,sx,ids ,ss)
      call searchdelta(isise,ddse,ddx,ss,idse,sse)
      call searchdelta(isie ,dde ,ddx,sx,ide ,se)
      call searchdelta(isien,dden,ddx,se,iden,sen)
      call searchdelta(isiw ,ddw ,ddx,sx,idw ,sw)
      call searchdelta(isiwn,ddwn,ddx,sw,idwn,swn)
    
      tnu=0.5*(sn+sne)
      tsu=0.5*(ss+sse) 
      tev=0.5*(se+sen) 
      twv=0.5*(sw+swn) 
      dsdxu(iq,ii)=(se-sx)*emu(iq)/ds
      dsdyu(iq,ii)=0.5*(tnu-tsu)*emu(iq)/ds
      dsdxv(iq,ii)=0.5*(tev-twv)*emv(iq)/ds
      dsdyv(iq,ii)=(sn-sx)*emv(iq)/ds
    end do
  end if
end do

return
end subroutine stagtruedelta

subroutine searchdelta(s,dep,dd,sx,id,ss)

use mlo

implicit none

integer, intent(inout) :: id
integer jj,fnd
real, dimension(wlev), intent(in) :: s,dep
real, intent(in) :: dd,sx
real, intent(out) :: ss
real xp
!integer, parameter :: vertintp=1 ! vertical interpolation (0=linear, 1=quadratic)

!if (dep(1)-0.01.gt.dd) return
if (dep(wlev).lt.dd) then
  ss=sx
  return
end if

fnd=wlev-1
do jj=id,wlev-2
  if (dep(jj+1).gt.dd) then
    fnd=jj
    exit
  end if
end do

! located a valid level to interpolate
id=fnd

!select case(vertintp)
!  case(0)
    xp=(dd-dep(id))/(dep(id+1)-dep(id))
    ss=(1.-xp)*s(id)+xp*s(id+1)
!  case(1)
!     !id=max(id,2)
!     ss=s(id)+(dd-dep(id))/(dep(id+1)-dep(id-1))*((dep(id)-dep(id-1))*(s(id+1)-s(id))/(dep(id+1)-dep(id)) &
!             +(dep(id+1)-dep(id))*(s(id)-s(id-1))/(dep(id)-dep(id-1)))                                    &
!             +(dd-dep(id))**2/(dep(id+1)-dep(id-1))*((s(id+1)-s(id))/(dep(id+1)-dep(id))                  &
!             -(s(id)-s(id-1))/(dep(id)-dep(id-1)))
!  case default
!    write(6,*) "ERROR: Invalid vertintp ",vertintp
!    stop
!end select

return
end subroutine

subroutine searchdeltab(s,dep,dd,sx,ss)

use mlo

implicit none

integer id,afnd,bfnd
real, dimension(wlev), intent(in) :: s,dep
real, intent(in) :: dd,sx
real, intent(out) :: ss
real xp
!integer, parameter :: vertintp=1 ! vertical interpolation (0=linear, 1=quadratic)

!if (dep(1)-0.01.gt.dd) return
if (dep(wlev).lt.dd) then
  ss=sx
  return
end if

!afnd=1
!bfnd=wlev
!do while (afnd.lt.bfnd-1)
!  id=0.5*(afnd+bfnd)
!  if (dep(id).le.dd) then
!    afnd=id
!  else
!    bfnd=id
!  end if
!end do
!id=afnd

afnd=wlev-1
do id=1,wlev-2
  if (dep(id+1).gt.dd) then
    afnd=id
    exit
  end if
end do
id=afnd

!select case(vertintp)
!  case(0)
    xp=(dd-dep(id))/(dep(id+1)-dep(id))
    ss=(1.-xp)*s(id)+xp*s(id+1)
!  case(1)
!     !id=max(id,2)
!     ss=s(id)+(dd-dep(id))/(dep(id+1)-dep(id-1))*((dep(id)-dep(id-1))*(s(id+1)-s(id))/(dep(id+1)-dep(id)) &
!             +(dep(id+1)-dep(id))*(s(id)-s(id-1))/(dep(id)-dep(id-1)))                                    &
!             +(dd-dep(id))**2/(dep(id+1)-dep(id-1))*((s(id+1)-s(id))/(dep(id+1)-dep(id))                  &
!             -(s(id)-s(id-1))/(dep(id)-dep(id-1)))
!  case default
!    write(6,*) "ERROR: Invalid vertintp ",vertintp
!    stop
!end select

return
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Tide

subroutine mlotide(eta,slon,slat,mtimer,jstart)

implicit none

include 'newmpar.h'

integer i
integer, intent(in) :: mtimer,jstart
real, dimension(ifull), intent(out) :: eta
real, dimension(ifull), intent(in) :: slon,slat
real :: ltime,stime,ctime,mn,sn,pn
real, dimension(ifull) :: sinlon,coslon,sinlat,coslat
real, dimension(4) :: aa,ab,ba,bb
real, parameter :: pi = 3.1415927

! amplitude
aa(1)=0.141565 ! K1
aa(2)=0.100661 ! O1
aa(3)=0.046848 ! P1
aa(4)=0.019273 ! Q1
ab(1)=0.242334 ! M2
ab(2)=0.112743 ! S2
ab(3)=0.046397 ! N2
ab(4)=0.030684 ! K2
! Love number
ba(1)=0.736
ba(2)=0.695
ba(3)=0.706
ba(4)=0.695
bb(1)=0.693
bb(2)=0.693
bb(3)=0.693
bb(4)=0.693
! Frequency
!wa(1)=0.7292117
!wa(2)=0.6750774
!wa(3)=0.7252295
!wa(4)=0.6495854
!wb(1)=1.405189
!wb(2)=1.454441 ! exactly twice per day solar day for S2
!wb(3)=1.378797
!wb(4)=1.458423

stime=real(mtimer)/1440.+real(jstart) ! solar time
ctime=stime/36525. ! century

! relative to 12Z 31 Dec 1899
mn=270.43659+481276.89057*ctime+0.00198*ctime*ctime+0.000002*ctime*ctime*ctime ! moon
sn=279.69668+36000.76892*ctime+0.00030*ctime*ctime                             ! sun
pn=334.32956+4069.03404*ctime-0.01032*ctime*ctime-0.00001*ctime*ctime*ctime    ! lunar perigee
mn=mn*pi/180.
sn=sn*pi/180.
pn=pn*pi/180.
stime=(stime-real(int(stime)))*2.*pi
ltime=stime+sn-mn
! 360*stime=360*ltime-9.3+12.2*day

coslat=cos(slat)**2
sinlat=sin(2.*slat)

eta=0.
eta=eta+ab(1)*aa(1)*sinlat*cos(ltime+mn-slon)             ! K1 (note sign change)
eta=eta-ab(2)*aa(2)*sinlat*cos(ltime-mn-slon)             ! O1
eta=eta-ab(3)*aa(3)*sinlat*cos(stime-sn-slon)             ! P1
eta=eta-ab(4)*aa(4)*sinlat*cos(ltime+pn-slon)             ! Q1
eta=eta-bb(1)*ab(1)*coslat*cos(2.*ltime-2.*slon)          ! M2
eta=eta-bb(2)*ab(2)*coslat*cos(2.*stime-2.*slon)          ! S2
eta=eta-bb(3)*ab(3)*coslat*cos(2.*ltime-mn+pn-2.*slon)    ! N2
eta=eta-bb(4)*ab(4)*coslat*cos(2.*stime+2.*sn-2.*slon)    ! K2

!sinlon=sin(2.*slon)
!coslon=cos(2.*slon)
!eta=0.
!do i=1,4
!  eta=eta-ba(i)*aa(i)*coslat*(cos(wa(i)*stime)*coslon-sin(wa(i)*stime)*sinlon)
!  eta=eta-bb(i)*ab(i)*sinlat*(cos(wb(i)*stime)*coslon-sin(wb(i)*stime)*sinlon)
!end do

return
end subroutine mlotide

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine fills ocean data over land points
subroutine mlofill(x,tst)

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
logical, dimension(ifull), intent(in) :: tst

miss=999999.

xx=miss
where (tst)
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine adjusts salinity to maintain a specified average
!
!subroutine mlosalfix
!
!use map_m
!use mlo
!use soil_m
!
!implicit none
!
!include 'mpif.h'
!include 'newmpar.h'
!include 'parm.h'
!
!integer iq,ii,ierr
!real salsuml,salsumg,voll,volg,adj
!real, dimension(ifull,wlev) :: sal,dz
!
!if (salfix.eq.0) return
!if (nud_sss.ne.0) return
!
!do ii=1,wlev
!  call mloexpdep(1,dz(:,ii),ii,0)
!  call mloexport(1,sal(:,ii),ii,0)
!end do
!dz=max(dz,1.E-3/real(wlev))
!
!voll=0.
!salsuml=0.
!do iq=1,ifull
!  if (.not.land(iq)) then
!    voll=voll+sum(dz(iq,:))/(em(iq)*em(iq))
!    salsuml=salsuml+sum(sal(iq,:)*dz(iq,:))/(em(iq)*em(iq))
!  end if
!end do
!
!call MPI_AllReduce(voll,volg,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
!call MPI_AllReduce(salsuml,salsumg,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
!
!adj=34.72-salsumg/volg
!
!do ii=1,wlev
!  sal(:,ii)=max(sal(:,ii)+adj,0.)
!  call mloimport(1,sal(:,ii),ii,0)
!end do
!
!return
!end subroutine

end module mlodynamics
