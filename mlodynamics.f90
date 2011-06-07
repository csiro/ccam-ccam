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
public mlodiffusion,mlorouter,mlohadv,watbdy,mlosalfix

real, dimension(:), allocatable, save :: watbdy
integer, parameter :: salfix = 0 ! fix salinity to 35 psu (0=off, 1=on)

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

integer k,i
integer, parameter :: salfilt = 0 ! Additional salinity filter (0=off, 1=Katzfey)
real, parameter :: k_smag = 0.4 ! 0.4 in mom3, 2. in Griffies (2000)
real hdif
real, dimension(ifull+iextra) :: uc,vc,wc,ee,dz,u,v,gg
real, dimension(ifull+iextra) :: t_kh,xfact,yfact
real, dimension(ifull) :: dudx,dvdx,dudy,dvdy,base
real, dimension(ifull) :: cc,ff,emi,ucc,vcc,wcc
logical, dimension(ifull+iextra) :: wtr

hdif=dt*(k_smag/pi)**2
ee=1.
where(land(1:ifull))
  ee(1:ifull)=0.
end where
call bounds(ee)
wtr=ee.gt.0.5

if (.not.any(wtr(1:ifull))) return

emi=1./em(1:ifull)**2
do k=1,wlev
  call mloexpdep(1,dz(1:ifull),k,0)
  call bounds(dz)
  dz=max(dz,1.E-3)

  ! Diffusion is currently applied in zstar=const with neta=0 plane.
  ! This avoids problems with the free surface, but creates issues with
  ! steep bathymetry.  At 60km resolution, we encounter 1/60 gradients
  ! in bathymetry which seems okay.
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
    dudx=(u(ieu)/em(ie)-u(iwu)/em(iw))*0.5*em(1:ifull)*em(1:ifull)/ds
    dvdx=(v(iev)/em(ie)-v(iwv)/em(iw))*0.5*em(1:ifull)*em(1:ifull)/ds
  else where (wtr(1:ifull).and.wtr(ie))
    dudx=(u(ieu)*em(1:ifull)/em(ie)-u(1:ifull))*em(1:ifull)/ds
    dvdx=(v(iev)*em(1:ifull)/em(ie)-u(1:ifull))*em(1:ifull)/ds
  else where (wtr(1:ifull).and.wtr(iw))
    dudx=(u(1:ifull)-u(iwu)*em(1:ifull)/em(iw))*em(1:ifull)/ds
    dvdx=(v(1:ifull)-v(iwv)*em(1:ifull)/em(iw))*em(1:ifull)/ds
  end where
  where (wtr(1:ifull).and.wtr(in).and.wtr(is))
    dudy=(u(inu)/em(in)-u(isu)/em(is))*0.5*em(1:ifull)*em(1:ifull)/ds
    dvdy=(v(inv)/em(in)-v(isv)/em(is))*0.5*em(1:ifull)*em(1:ifull)/ds
  else where (wtr(1:ifull).and.wtr(ie))
    dudy=(u(inu)*em(1:ifull)/em(in)-u(1:ifull))*em(1:ifull)/ds
    dvdy=(v(inv)*em(1:ifull)/em(in)-u(1:ifull))*em(1:ifull)/ds
  else where (wtr(1:ifull).and.wtr(iw))
    dudy=(u(1:ifull)-u(isu)*em(1:ifull)/em(is))*em(1:ifull)/ds
    dvdy=(v(1:ifull)-v(isv)*em(1:ifull)/em(is))*em(1:ifull)/ds
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
  u(1:ifull) = ax(1:ifull)*ucc + ay(1:ifull)*vcc + az(1:ifull)*wcc
  v(1:ifull) = bx(1:ifull)*ucc + by(1:ifull)*vcc + bz(1:ifull)*wcc
  
  call mloimport(2,u(1:ifull),k,0)
  call mloimport(3,v(1:ifull),k,0)
   
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
  flow(:,i)=watbdy(1:ifull)/(dp(:,i)/(-vel(:,i)*dt)-1.) ! (kg/m^2)
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
integer, dimension(ifull) :: wwn,wwe,wws,www
real alpha,maxloclseta,maxglobseta,maxloclip,maxglobip
real delpos,delneg,alph_p,dumpp,dumpn
real, dimension(ifull+iextra) :: ee,neta,dd,snu,snv,ntide,pice
real, dimension(ifull+iextra) :: nfracice,ndic,ndsn,nsto,niu,niv,ndum
real, dimension(ifull+iextra) :: imass,spu,squ,sru,spv,sqv,srv
real, dimension(:), allocatable, save :: ip
real, dimension(ifull) :: i_u,i_v,i_sto,rhobaru,rhobarv
real, dimension(ifull) :: div,seta,w_e
real, dimension(ifull) :: tnu,tsu,tev,twv,rhou,rhov,sou,sov
real, dimension(ifull) :: dpsdxu,dpsdyu,dpsdxv,dpsdyv
real, dimension(ifull) :: dttdxu,dttdyu,dttdxv,dttdyv
real, dimension(ifull) :: detadxu,detadyu,detadxv,detadyv
real, dimension(ifull) :: dipdxu,dipdyu,dipdxv,dipdyv
real, dimension(ifull) :: au,bu,cu,av,bv,cv,odum
real, dimension(ifull) :: nip,iia,dibdx,dibdy,dicdx,dicdy,ipn,ipe,ips,ipw,ipmax
real, dimension(ifull) :: dsnudeta,dsnuwdeta,dsnvdeta,dsnvsdeta,ddivdeta
real, dimension(ifull) :: sssa,sssb,sssc,sssd
real, dimension(ifull+iextra,4) :: nit
real, dimension(ifull+iextra,wlev) :: nu,nv,nt,ns
real, dimension(ifull+iextra,wlev) :: cou,cov,cow
real, dimension(ifull+iextra,wlev) :: dep,rhobar,rho,dz
real, dimension(ifull,1) :: siu,siv
real, dimension(ifull,4) :: i_it
real, dimension(ifull,wlev) :: w_u,w_v,w_t,w_s,nw,dzbar,dum,dumb
real, dimension(ifull,wlev) :: nuh,nvh,xg,yg,uau,uav
real, dimension(ifull,wlev) :: kku,llu,mmu,nnu
real, dimension(ifull,wlev) :: kkv,llv,mmv,nnv
real, dimension(ifull,wlev) :: drhobardxu,drhobardyu,drhobardxv,drhobardyv
real, dimension(ifull,4,4) :: stwgt
real, dimension(:,:), allocatable, save :: oldu1,oldu2,oldv1,oldv2
real*8, dimension(ifull,wlev) :: x3d,y3d,z3d
logical, dimension(ifull+iextra) :: wtr
logical, dimension(ifull) :: stest
logical lleap
integer, parameter :: usetide=1    ! tidal forcing (0=Off, 1=On)
integer, parameter :: icemode=2    ! Ice stress (0=free-drift, 1=incompressible, 2=cavitating)
integer, parameter :: drmeth=1     ! rho gradient calculation (0=basic, 1=interpolated)
integer, parameter :: lmax=1       ! 1=predictor-only, 2+=predictor-corrector
integer, parameter :: llmax=500    ! iterations for calculating surface height
real, parameter :: tol = 5.E-4     ! Tolerance for GS solver (water)
real, parameter :: itol = 5.E-7    ! Tolerance for GS solver (ice)
real, parameter :: sal = 0.948     ! SAL parameter
real, parameter :: rhosn=330.      ! density snow
real, parameter :: rhoic=900.      ! density ice

data ndoy/0,31,59,90,120,151,181,212,243,273,304,334/

! new z levels for including free surface eta (effectively sigma-height levels)
! newz=-eta+oldz*(1+eta/maxdepth)
! where 0<=oldz<=maxdepth and -eta<=newz<=maxdepth
! hence newz is the true z position and oldz measures the depth below the free surface

intsch=mod(ktau,2)

!Define land/sea mask
ee=1.
where(land(1:ifull))
  ee(1:ifull)=0.
end where
call bounds(ee,nrows=2)
wtr=ee.gt.0.5

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

! PRECOMPUTE WEIGHTS FOR CALCULATING GRADIENTS
stwgt=0.
where (wtr(in).and.wtr(ine))
  stwgt(:,1,2)=0.5 !in
  stwgt(:,1,3)=0.5 !ine
else where (wtr(in))
  stwgt(:,1,2)=1.  !in
else where (wtr(ine))
  stwgt(:,1,3)=1.  !ine
else where
  stwgt(:,1,1)=0.5 !ifull
  stwgt(:,1,4)=0.5 !ie
end where
where (wtr(is).and.wtr(ise))
  stwgt(:,2,2)=0.5 !is
  stwgt(:,2,3)=0.5 !ise
else where (wtr(is))
  stwgt(:,2,2)=1.  !is
else where (wtr(ise))
  stwgt(:,2,3)=1.  !ise
else where
  stwgt(:,2,1)=0.5 !ifull
  stwgt(:,2,4)=0.5 !ie
end where
where (wtr(ie).and.wtr(ien))
  stwgt(:,3,2)=0.5 !ie
  stwgt(:,3,3)=0.5 !ien
else where (wtr(ie))
  stwgt(:,3,2)=1.  !ie
else where (wtr(ien))
  stwgt(:,3,3)=1.  !ien
else where
  stwgt(:,3,1)=0.5 !ifull
  stwgt(:,3,4)=0.5 !in
end where
where (wtr(iw).and.wtr(iwn))
  stwgt(:,4,2)=0.5 !iw
  stwgt(:,4,3)=0.5 !iwn
else where (wtr(iw))
  stwgt(:,4,2)=1.  !iw
else where (wtr(iwn))
  stwgt(:,4,3)=1.  !iwn
else where
  stwgt(:,4,1)=0.5 !ifull
  stwgt(:,4,4)=0.5 !in
end where

! precompute ice indices
wwn=in
wwe=ie
wws=is
www=iw
do iq=1,ifull
  if (.not.wtr(wwn(iq))) wwn(iq)=iq
  if (.not.wtr(wwe(iq))) wwe(iq)=iq
  if (.not.wtr(wws(iq))) wws(iq)=iq
  if (.not.wtr(www(iq))) www(iq)=iq
end do

! ADVECT WATER AND ICE ----------------------------------------------
do ii=1,wlev
  call mloexpdep(0,dep(1:ifull,ii),ii,0)
  call mloexpdep(1,dz(:,ii),ii,0)
  call mloexport(0,w_t(:,ii),ii,0)
  call mlofill(w_t(:,ii),wtr(1:ifull))
  call mloexport(1,w_s(:,ii),ii,0)
  call mlofill(w_s(:,ii),wtr(1:ifull))
  call mloexport(2,w_u(:,ii),ii,0)
  call mloexport(3,w_v(:,ii),ii,0)
end do
call mloexport(4,w_e,0,0)
do ii=1,4
  call mloexpice(i_it(:,ii),ii,0)
  call mlofill(i_it(:,ii),fracice.gt.1.E-10)
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
call bounds(dd)

if (usetide.eq.1) then
  ! estimate tides
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
  call mlotide(ntide(1:ifull),rlongg,rlatt,mins,jstart)
  ntide(1:ifull)=ntide(1:ifull)*ee(1:ifull)
  call bounds(ntide,nrows=2)
  
  tnu=stwgt(:,1,1)*ntide(1:ifull)+stwgt(:,1,2)*ntide(in)+stwgt(:,1,3)*ntide(ine)+stwgt(:,1,4)*ntide(ie)
  tsu=stwgt(:,2,1)*ntide(1:ifull)+stwgt(:,2,2)*ntide(is)+stwgt(:,2,3)*ntide(ise)+stwgt(:,2,4)*ntide(ie)
  tev=stwgt(:,3,1)*ntide(1:ifull)+stwgt(:,3,2)*ntide(ie)+stwgt(:,3,3)*ntide(ien)+stwgt(:,3,4)*ntide(in)
  twv=stwgt(:,4,1)*ntide(1:ifull)+stwgt(:,4,2)*ntide(iw)+stwgt(:,4,3)*ntide(iwn)+stwgt(:,4,4)*ntide(in)
  dttdxu=(ntide(ie)-ntide(1:ifull))*emu(1:ifull)/ds
  dttdyu=0.5*(tnu-tsu)*emu(1:ifull)/ds
  dttdxv=0.5*(tev-twv)*emv(1:ifull)/ds
  dttdyv=(ntide(in)-ntide(1:ifull))*emv(1:ifull)/ds  
else
  ntide=0.
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
if (.not.allocated(ip)) then
  allocate(ip(ifull+iextra))
  ip=0. ! free drift solution for ice
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

  ! surface pressure gradients
  imass(1:ifull)=ndic(1:ifull)*rhoic+ndsn(1:ifull)*rhosn ! ice mass per unit area (kg/m^2)
  pice(1:ifull) = ps(1:ifull) + grav*nfracice(1:ifull)*imass(1:ifull) ! pressure due to atmosphere and ice at top of water column
  call bounds(pice,nrows=2)
  tnu=stwgt(:,1,1)*pice(1:ifull)+stwgt(:,1,2)*pice(in)+stwgt(:,1,3)*pice(ine)+stwgt(:,1,4)*pice(ie)
  tsu=stwgt(:,2,1)*pice(1:ifull)+stwgt(:,2,2)*pice(is)+stwgt(:,2,3)*pice(ise)+stwgt(:,2,4)*pice(ie)
  tev=stwgt(:,3,1)*pice(1:ifull)+stwgt(:,3,2)*pice(ie)+stwgt(:,3,3)*pice(ien)+stwgt(:,3,4)*pice(in)
  twv=stwgt(:,4,1)*pice(1:ifull)+stwgt(:,4,2)*pice(iw)+stwgt(:,4,3)*pice(iwn)+stwgt(:,4,4)*pice(in)
  dpsdxu=(pice(ie)-pice(1:ifull))*emu(1:ifull)/ds
  dpsdyu=0.5*(tnu-tsu)*emu(1:ifull)/ds
  dpsdxv=0.5*(tev-twv)*emv(1:ifull)/ds
  dpsdyv=(pice(in)-pice(1:ifull))*emv(1:ifull)/ds


  ! calculate rhobar (unstaggered)
  ! note that rho is estimated using the current guess for neta
  do ii=1,wlev
    dum(:,ii)=dep(1:ifull,ii)*max(1.+neta(1:ifull)/dd(1:ifull),0.01) ! depth below free surface =newz+eta, and ranges between 0 and dd+neta
    dumb(:,ii)=dz(1:ifull,ii)*max(1.+neta(1:ifull)/dd(1:ifull),0.01)
  end do
  call mloexpdensity(rho(1:ifull,:),nt(1:ifull,:),ns(1:ifull,:),dum,dumb,pice(1:ifull),0)
  call bounds(rho)
  rhobar(1:ifull,1)=rho(1:ifull,1)*dumb(:,1)
  do ii=2,wlev
    rhobar(1:ifull,ii)=rhobar(1:ifull,ii-1)+rho(1:ifull,ii)*dumb(:,ii)
  end do
  do ii=1,wlev
    rhobar(1:ifull,ii)=rhobar(1:ifull,ii)/(dzbar(:,ii)*max(1.+neta(1:ifull)/dd(1:ifull),0.01))
  end do
  call bounds(rhobar,nrows=2)

  ! ADVECT WATER ----------------------------------------------------
  ! Water currents are advected using sem-Lagrangian advection
  ! (i.e., taken from McGregor's CCAM advection routines).
  ! Velocity is set to zero at ocean boundaries, and changing
  ! bathymetry is accounted for in the integral of the divergence.

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
  call getww(cou,cov,w_e,dd(1:ifull),dz(1:ifull,:),ee(1:ifull),nw)

  ! Vertical advection (first call)
  call mlovadv(0.5*dt,nw,nu(1:ifull,:),nv(1:ifull,:),ns(1:ifull,:),nt(1:ifull,:),dep(1:ifull,:),dz(1:ifull,:), &
               w_e,dd(1:ifull),wtr(1:ifull))

  ! Convert (u,v) to cartesian coordinates (U,V,W)
  do ii=1,wlev
    cou(1:ifull,ii)=ax(1:ifull)*nu(1:ifull,ii)+bx(1:ifull)*nv(1:ifull,ii)
    cov(1:ifull,ii)=ay(1:ifull)*nu(1:ifull,ii)+by(1:ifull)*nv(1:ifull,ii)
    cow(1:ifull,ii)=az(1:ifull)*nu(1:ifull,ii)+bz(1:ifull)*nv(1:ifull,ii)
  end do

  ! Horizontal advection for U,V,W
  call mloints(cou,intsch,nface,xg,yg,2)
  call mloints(cov,intsch,nface,xg,yg,2)
  call mloints(cow,intsch,nface,xg,yg,2)

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
  call mloints(nt,intsch,nface,xg,yg,2)
  call mloints(ns,intsch,nface,xg,yg,5)
  ns=max(ns,0.)
  
    
  ! PREP WATER ARRAYS -----------------------------------------------

  ! stagger the following values for t=tstar
  ! This is the same as eps=1. in JLM's atmosphere semi-Lagrangian dynamics
  do ii=1,wlev
    uau(:,ii)=nu(1:ifull,ii)+dt*f(1:ifull)*nv(1:ifull,ii)
    uav(:,ii)=nv(1:ifull,ii)-dt*f(1:ifull)*nu(1:ifull,ii)
  end do
  call mlostaguv(uau,uav,cou(1:ifull,:),cov(1:ifull,:))
  uau(:,1)=i_u+dt*f(1:ifull)*i_v
  uav(:,1)=i_v-dt*f(1:ifull)*i_u
  call mlostaguv(uau(:,1:1),uav(:,1:1),siu,siv)

  select case(drmeth)
    case(0)
      do ii=1,wlev
        tnu=stwgt(:,1,1)*rhobar(1:ifull,ii)+stwgt(:,1,2)*rhobar(in,ii)+stwgt(:,1,3)*rhobar(ine,ii)+stwgt(:,1,4)*rhobar(ie,ii)
        tsu=stwgt(:,2,1)*rhobar(1:ifull,ii)+stwgt(:,2,2)*rhobar(is,ii)+stwgt(:,2,3)*rhobar(ise,ii)+stwgt(:,2,4)*rhobar(ie,ii)
        tev=stwgt(:,3,1)*rhobar(1:ifull,ii)+stwgt(:,3,2)*rhobar(ie,ii)+stwgt(:,3,3)*rhobar(ien,ii)+stwgt(:,3,4)*rhobar(in,ii)
        twv=stwgt(:,4,1)*rhobar(1:ifull,ii)+stwgt(:,4,2)*rhobar(iw,ii)+stwgt(:,4,3)*rhobar(iwn,ii)+stwgt(:,4,4)*rhobar(in,ii)
        drhobardxu(:,ii)=(rhobar(ie,ii)-rhobar(1:ifull,ii))*emu(1:ifull)/ds
        drhobardyu(:,ii)=0.5*(tnu-tsu)*emu(1:ifull)/ds
        drhobardxv(:,ii)=0.5*(tev-twv)*emv(1:ifull)/ds
        drhobardyv(:,ii)=(rhobar(in,ii)-rhobar(1:ifull,ii))*emv(1:ifull)/ds
      end do
    case(1)
      ! interpolate density gradient to horizontal surfaces
      ! (note that dep should be -eta+dep*(1+eta/maxdep), so we assume eta is small compared to dep)
      call stagtruedelta(rhobar,dep,wtr,drhobardxu,drhobardyu,drhobardxv,drhobardyv)
    case default
      write(6,*) "ERROR: Invalid drmeth ",drmeth
      stop
  end select
  
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

    !nu=au+bu*dppdxu+cu*dppdyu
    !nv=av+bv*dppdyv+cv*dppdxv

    !dppdxu=dpsdxu+grav*depu*(1+etau/ddu)*drhobardxu+grav*(rhobaru+(1-sal))*detadxu+grav*dttdxu
    !dppdyu=dpsdyu+grav*depu*(1+etau/ddu)*drhobardyu+grav*(rhobaru+(1-sal))*detadyu+grav*dttdyu
    !dppdxv=dpsdxv+grav*depv*(1+etav/ddv)*drhobardxv+grav*(rhobarv+(1-sal))*detadxv+grav*dttdxv
    !dppdyv=dpsdyv+grav*depv*(1+etav/ddv)*drhobardyv+grav*(rhobarv+(1-sal))*detadyv+grav*dttdyv
    
    !nu=kku+llu*(1+etau/ddu)+mmu*detadxu+nnu*detadyu
    !nv=kkv+llv*(1+etav/ddv)+mmv*detadyv+nnv*detadxv
    
    !int nu dz = sou+spu*(1+etau/ddu)+squ*detadxu+sru*detadyu
    !int nv dz = sov+spv*(1+etav/ddv)+sqv*detadyv+srv*detadxv

    if (usetide.eq.1) then
      kku(:,ii)=au+bu*(dpsdxu+grav*rhou*dttdxu)+cu*(dpsdyu+grav*rhou*dttdyu)
      llu(:,ii)=(bu*grav*drhobardxu(:,ii)+cu*grav*drhobardyu(:,ii))*0.5*(dep(1:ifull,ii)+dep(ie,ii))
      mmu(:,ii)=bu*grav*(rhobaru+rhou*(1.-sal))
      nnu(:,ii)=cu*grav*(rhobaru+rhou*(1.-sal))

      kkv(:,ii)=av+bv*(dpsdyv+grav*rhov*dttdyv)+cv*(dpsdxv+grav*rhov*dttdxv)
      llv(:,ii)=(bv*grav*drhobardyv(:,ii)+cv*grav*drhobardxv(:,ii))*0.5*(dep(1:ifull,ii)+dep(in,ii))
      mmv(:,ii)=bv*grav*(rhobarv+rhov*(1.-sal))
      nnv(:,ii)=cv*grav*(rhobarv+rhov*(1.-sal))
    else
      kku(:,ii)=au+bu*dpsdxu+cu*dpsdyu
      llu(:,ii)=(bu*grav*drhobardxu(:,ii)+cu*grav*drhobardyu(:,ii))*0.5*(dep(1:ifull,ii)+dep(ie,ii))
      mmu(:,ii)=bu*grav*rhobaru
      nnu(:,ii)=cu*grav*rhobaru

      kkv(:,ii)=av+bv*dpsdyv+cv*dpsdxv
      llv(:,ii)=(bv*grav*drhobardyv(:,ii)+cv*grav*drhobardxv(:,ii))*0.5*(dep(1:ifull,ii)+dep(in,ii))
      mmv(:,ii)=bv*grav*rhobarv
      nnv(:,ii)=cv*grav*rhobarv
    end if
    
    sou=sou+kku(:,ii)*0.5*(dz(1:ifull,ii)+dz(ie,ii))
    spu(1:ifull)=spu(1:ifull)+llu(:,ii)*0.5*(dz(1:ifull,ii)+dz(ie,ii))
    squ(1:ifull)=squ(1:ifull)+mmu(:,ii)*0.5*(dz(1:ifull,ii)+dz(ie,ii))
    sru(1:ifull)=sru(1:ifull)+nnu(:,ii)*0.5*(dz(1:ifull,ii)+dz(ie,ii))
    
    sov=sov+kkv(:,ii)*0.5*(dz(1:ifull,ii)+dz(in,ii))
    spv(1:ifull)=spv(1:ifull)+llv(:,ii)*0.5*(dz(1:ifull,ii)+dz(in,ii))
    sqv(1:ifull)=sqv(1:ifull)+mmv(:,ii)*0.5*(dz(1:ifull,ii)+dz(in,ii))
    srv(1:ifull)=srv(1:ifull)+nnv(:,ii)*0.5*(dz(1:ifull,ii)+dz(in,ii))
    
  end do
  
  ! update land boundaries
  sou=sou*ee(1:ifull)*ee(ie)
  spu(1:ifull)=spu(1:ifull)*ee(1:ifull)*ee(ie)
  squ(1:ifull)=squ(1:ifull)*ee(1:ifull)*ee(ie)
  sru(1:ifull)=sru(1:ifull)*ee(1:ifull)*ee(ie)
  sov=sov*ee(1:ifull)*ee(in)
  spv(1:ifull)=spv(1:ifull)*ee(1:ifull)*ee(in)
  sqv(1:ifull)=sqv(1:ifull)*ee(1:ifull)*ee(in)
  srv(1:ifull)=srv(1:ifull)*ee(1:ifull)*ee(in)

  call boundsuv(spu,spv)
  call boundsuv(squ,sqv)
  call boundsuv(sru,srv)

  sssa=spu(1:ifull)/(dd(ie)+dd(1:ifull))-squ(1:ifull)*emu(1:ifull)/ds+sru(1:ifull)*0.5*emu(1:ifull)/ds*(stwgt(:,1,1)-stwgt(:,2,1))
  sssb=spu(iwu)/(dd(iw)+dd(1:ifull))+squ(iwu)*emu(iwu)/ds+sru(iwu)*0.5*emu(iwu)/ds*(stwgt(:,1,4)-stwgt(:,2,4))
  sssc=spv(1:ifull)/(dd(in)+dd(1:ifull))-sqv(1:ifull)*emv(1:ifull)/ds+srv(1:ifull)*0.5*emv(1:ifull)/ds*(stwgt(:,3,1)-stwgt(:,4,1))
  sssd=spv(isv)/(dd(is)+dd(1:ifull))+sqv(isv)*emv(isv)/ds+srv(isv)*0.5*emv(isv)/ds*(stwgt(:,3,4)-stwgt(:,4,4))

  ! Now iteratively solve for eta and ice stress terms
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
    tnu=stwgt(:,1,1)*neta(1:ifull)+stwgt(:,1,2)*neta(in)+stwgt(:,1,3)*neta(ine)+stwgt(:,1,4)*neta(ie)
    tsu=stwgt(:,2,1)*neta(1:ifull)+stwgt(:,2,2)*neta(is)+stwgt(:,2,3)*neta(ise)+stwgt(:,2,4)*neta(ie)
    tev=stwgt(:,3,1)*neta(1:ifull)+stwgt(:,3,2)*neta(ie)+stwgt(:,3,3)*neta(ien)+stwgt(:,3,4)*neta(in)
    twv=stwgt(:,4,1)*neta(1:ifull)+stwgt(:,4,2)*neta(iw)+stwgt(:,4,3)*neta(iwn)+stwgt(:,4,4)*neta(in)
 
    detadxu=(neta(ie)-neta(1:ifull))*emu(1:ifull)/ds
    detadyu=0.5*(tnu-tsu)*emu(1:ifull)/ds
    detadxv=0.5*(tev-twv)*emv(1:ifull)/ds
    detadyv=(neta(in)-neta(1:ifull))*emv(1:ifull)/ds

    snu(1:ifull)=sou+spu(1:ifull)*(1.+(neta(ie)+neta(1:ifull))/(dd(ie)+dd(1:ifull)))+squ(1:ifull)*detadxu+sru(1:ifull)*detadyu
    snv(1:ifull)=sov+spv(1:ifull)*(1.+(neta(in)+neta(1:ifull))/(dd(in)+dd(1:ifull)))+sqv(1:ifull)*detadyv+srv(1:ifull)*detadxv
    snu(1:ifull)=snu(1:ifull)*(1.+(neta(ie)+neta(1:ifull))/(dd(ie)+dd(1:ifull))) ! conserve volume
    snv(1:ifull)=snv(1:ifull)*(1.+(neta(in)+neta(1:ifull))/(dd(in)+dd(1:ifull))) ! conserve volume

    call boundsuv(snu,snv)
    div=(snu(1:ifull)/emu(1:ifull)-snu(iwu)/emu(iwu)+snv(1:ifull)/emv(1:ifull)-snv(isv)/emv(isv)) &
        *em(1:ifull)*em(1:ifull)/ds

    ! Update neta
    dsnudeta=sssa*(1.+(neta(ie)+neta(1:ifull))/(dd(ie)+dd(1:ifull)))+snu(1:ifull)/(dd(ie)+dd(1:ifull))
    dsnuwdeta=sssb*(1.+(neta(iw)+neta(1:ifull))/(dd(iw)+dd(1:ifull)))+snu(iwu)/(dd(iw)+dd(1:ifull))
    dsnvdeta=sssc*(1.+(neta(in)+neta(1:ifull))/(dd(in)+dd(1:ifull)))+snv(1:ifull)/(dd(in)+dd(1:ifull))
    dsnvsdeta=sssd*(1.+(neta(is)+neta(1:ifull))/(dd(is)+dd(1:ifull)))+snv(isv)/(dd(is)+dd(1:ifull))
    ddivdeta=(dsnudeta/emu(1:ifull)-dsnuwdeta/emu(iwu)+dsnvdeta/emv(1:ifull)-dsnvsdeta/emv(isv))*em(1:ifull)*em(1:ifull)/ds
    seta=-neta(1:ifull)+w_e-div/(1./dt+ddivdeta)

    ! The following expression limits the minimum depth to 1m
    seta=max(seta,(1.-dd(1:ifull)-neta(1:ifull))/alpha) ! this should become a land point
    seta=min(max(seta,-1.),1.)
    neta(1:ifull)=alpha*seta+neta(1:ifull)
    neta(1:ifull)=neta(1:ifull)*ee(1:ifull)
    seta=seta*ee(1:ifull)
    
    ! Break iterative loop when maximum error is below tol (expensive)
    maxloclseta=maxval(abs(seta))
    call MPI_AllReduce(maxloclseta,maxglobseta,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
    if (maxglobseta.lt.tol.and.ll.gt.2) exit    

    totits=totits+1
  end do

  ! Update currents once neta is calculated
  do ii=1,wlev
    ! update currents (staggered)
    nu(1:ifull,ii)=kku(:,ii)+llu(:,ii)*(1.+(neta(ie)+neta(1:ifull))/(dd(ie)+dd(1:ifull)))+mmu(:,ii)*detadxu+nnu(:,ii)*detadyu
    nv(1:ifull,ii)=kkv(:,ii)+llv(:,ii)*(1.+(neta(in)+neta(1:ifull))/(dd(in)+dd(1:ifull)))+mmv(:,ii)*detadyv+nnv(:,ii)*detadxv
    ! update land boundaries
    nu(1:ifull,ii)=nu(1:ifull,ii)*ee(1:ifull)*ee(ie)
    nv(1:ifull,ii)=nv(1:ifull,ii)*ee(1:ifull)*ee(in)
  end do
  call boundsuv(nu,nv)
  call getww(nu,nv,neta(1:ifull),dd(1:ifull),dz(1:ifull,:),ee(1:ifull),nw)
  ! unstagger nu and nv
  call mlounstaguv(nu(1:ifull,:),nv(1:ifull,:),nuh,nvh)
  nu(1:ifull,:)=nuh
  nv(1:ifull,:)=nvh

  ! Vertical advection (second call)
  call mlovadv(0.5*dt,nw,nu(1:ifull,:),nv(1:ifull,:),ns(1:ifull,:),nt(1:ifull,:),dep(1:ifull,:),dz(1:ifull,:), &
               neta(1:ifull),dd(1:ifull),wtr(1:ifull))

  ! UPDATE ICE DYNAMICS ---------------------------------------------
  ! Here we start by calculating the ice velocity and then advecting
  ! the various ice prognostic variables.

  ! Update ice velocities
  ! niu and niv hold the free drift solution
  niu(1:ifull)=(siu(:,1)-dt*grav*(detadxu+dt*fu(1:ifull)*detadyu))/(1.+dt*dt*fu(1:ifull)*fu(1:ifull))
  niv(1:ifull)=(siv(:,1)-dt*grav*(detadyv-dt*fv(1:ifull)*detadxv))/(1.+dt*dt*fv(1:ifull)*fv(1:ifull))
  niu(1:ifull)=niu(1:ifull)*ee(1:ifull)*ee(ie)
  niv(1:ifull)=niv(1:ifull)*ee(1:ifull)*ee(in)
  call boundsuv(niu,niv)
  
  imass(1:ifull)=max(imass(1:ifull),10.)
  if (icemode.gt.0) then
    call bounds(imass)
    iia=imass(1:ifull)*(niu(1:ifull)-niu(iwu)+niv(1:ifull)-niv(isv))*0.5*em(1:ifull)/ds
    !iib=dt               !ibu=ibv=iib
    !iic=dt*dt*f(1:ifull) !icu=-icv=iic
    dibdx=dt*imass(1:ifull)*(1./(imass(wwe)*em(wwe))-1./(imass(www)*em(www)))*0.5*em(1:ifull)**2/ds
    dicdx=dt*dt*imass(1:ifull)*(f(wwe)/(imass(wwe)*em(wwe))-f(www)/(imass(www)*em(www)))*0.5*em(1:ifull)**2/ds
    dibdy=dt*imass(1:ifull)*(1./(imass(wwn)*em(wwn))-1./(imass(wws)*em(wws)))*0.5*em(1:ifull)**2/ds
    dicdy=dt*dt*imass(1:ifull)*(f(wwn)/(imass(wwn)*em(wwn))-f(wws)/(imass(wws)*em(wws)))*0.5*em(1:ifull)**2/ds
  
    alpha=0.9
    do ll=1,llmax
  
      call bounds(ip)
      ipn=ip(wwn)/em(wwn)
      ipe=ip(wwe)/em(wwe)
      ips=ip(wws)/em(wws)
      ipw=ip(www)/em(www)
   
      div=iia-dt*(ipe-2.*ip(1:ifull)+ipw)*em(1:ifull)**3/(ds*ds)                          &
             -dibdx*(ipe-ipw)*0.5*em(1:ifull)**2/ds-dicdx*(ipn-ips)*0.5*em(1:ifull)**2/ds &
             -dt*(ipn-2.*ip(1:ifull)+ips)*em(1:ifull)**3/(ds*ds)                          &
             -dibdy*(ipn-ips)*0.5*em(1:ifull)**2/ds+dicdy*(ipe-ipw)*0.5*em(1:ifull)**2/ds
      nip=ip(1:ifull)
      where (div.lt.0..and.wtr(1:ifull))
        nip=(-iia+dibdx*(ipe-ipw)*0.5*em(1:ifull)**2/ds+dt*(ipe+ipw)*em(1:ifull)**3/(ds*ds) &
                 +dicdx*(ipn-ips)*0.5*em(1:ifull)**2/ds                                     &
                 +dibdy*(ipn-ips)*0.5*em(1:ifull)**2/ds+dt*(ipn+ips)*em(1:ifull)**3/(ds*ds) &
                 -dicdy*(ipe-ipw)*0.5*em(1:ifull)**2/ds)                                    &
            /(4.*dt*em(1:ifull)**3/(ds*ds))
      end where
      if (icemode.gt.1) then
        ipmax=27500.*ndic(1:ifull)*exp(-20.*(1.-nfracice(1:ifull)))
        stest=ip(1:ifull).lt.ipmax
        if (any(stest)) then
          maxloclip=maxval(-div,stest)
        else
          maxloclip=0.
        end if
        ip(1:ifull)=max(min(alpha*nip+(1.-alpha)*ip(1:ifull),ipmax),0.)
      else
        ip(1:ifull)=max(alpha*nip+(1.-alpha)*ip(1:ifull),0.)
        maxloclip=maxval(-div)
      end if
  
      
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
  tnu=stwgt(:,1,1)*ip(1:ifull)+stwgt(:,1,2)*ip(in)+stwgt(:,1,3)*ip(ine)+stwgt(:,1,4)*ip(ie)
  tsu=stwgt(:,2,1)*ip(1:ifull)+stwgt(:,2,2)*ip(is)+stwgt(:,2,3)*ip(ise)+stwgt(:,2,4)*ip(ie)
  tev=stwgt(:,3,1)*ip(1:ifull)+stwgt(:,3,2)*ip(ie)+stwgt(:,3,3)*ip(ien)+stwgt(:,3,4)*ip(in)
  twv=stwgt(:,4,1)*ip(1:ifull)+stwgt(:,4,2)*ip(iw)+stwgt(:,4,3)*ip(iwn)+stwgt(:,4,4)*ip(in)
  dipdxu=(ip(ie)-ip(1:ifull))*emu(1:ifull)/ds
  dipdyu=0.5*(tnu-tsu)*emu(1:ifull)/ds
  dipdxv=0.5*(tev-twv)*emv(1:ifull)/ds
  dipdyv=(ip(in)-ip(1:ifull))*emv(1:ifull)/ds
  
  ! niu and niv hold the corrected solution after the plastic terms are accounted for
  niu(1:ifull)=niu(1:ifull)-(dt*dipdxu/imass(1:ifull)+dt*dt*fu(1:ifull)*dipdyu/imass(1:ifull)) &
                            /(1.+dt*dt*fu(1:ifull)*fu(1:ifull))
  niv(1:ifull)=niv(1:ifull)-(dt*dipdyv/imass(1:ifull)-dt*dt*fv(1:ifull)*dipdxv/imass(1:ifull)) &
                            /(1.+dt*dt*fv(1:ifull)*fv(1:ifull))
  niu(1:ifull)=niu(1:ifull)*ee(1:ifull)*ee(ie)
  niv(1:ifull)=niv(1:ifull)*ee(1:ifull)*ee(in)
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

! salinity conservation
if (nud_sss.eq.0) then
  dum=0.
  delpos=0.
  delneg=0.
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
  ns(1:ifull,:)=w_s(1:ifull,:)+alph_p*max(0.,dum)+min(0.,dum)/max(1.,alph_p)
end if

! volume conservation for water ---------------------------------------
! (include mass conservation in mlo.f90 due to thermal change)
if (nud_sfh.eq.0) then
  odum=neta(1:ifull)-w_e
  call ccglobal_posneg(odum,delpos,delneg)
  alph_p = sqrt( -delneg/max(1.e-20,delpos) )
  neta(1:ifull) = w_e + alph_p*max(0.,odum) + min(0.,odum)/max(1.e-20,alph_p)
end if

call MPI_AllReduce(maxloclseta,maxglobseta,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(maxloclip,maxglobip,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
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

subroutine mlodeps(ubar,vbar,nface,xg,yg,x3d,y3d,z3d)

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
real, dimension(:,:), intent(in) :: ubar,vbar
real, dimension(:,:), intent(out) :: xg,yg
real*8, dimension(:,:), intent(out) :: x3d,y3d,z3d
real, dimension(ifull,size(ubar,2)) :: uc,vc,wc
real, dimension(ifull+iextra,size(ubar,2)) :: temp
integer, parameter :: nguess = 2

kx=size(ubar,2)
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
  call mloints(temp,intsch,nface,xg,yg,2)
  do ii=1,kx
    x3d(:,ii) = x - 0.5*(uc(:,ii)+temp(1:ifull,ii)) ! n+1 guess
  end do
  temp(1:ifull,:) = vc
  call mloints(temp,intsch,nface,xg,yg,2)
  do ii=1,kx
    y3d(:,ii) = y - 0.5*(vc(:,ii)+temp(1:ifull,ii)) ! n+1 guess
  end do
  temp(1:ifull,:) = wc
  call mloints(temp,intsch,nface,xg,yg,2)
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
integer idel,iq,jdel,nn,kx
integer i,j,k,n,ind,ip,jp,iproc,ierr
integer, dimension(:,:), intent(in) :: nface
real, dimension(:,:), intent(in) :: xg,yg
real, dimension(:,:), intent(inout) :: s
real, dimension(-1:ipan+2,-1:jpan+2,1:npan,size(s,2)) :: sx
real, dimension(4) :: r
real a3,a4,c1,c2,c3,c4,cmax,cmin,sss,xxg,yyg
ind(i,j,n)=i+(j-1)*ipan+(n-1)*ipan*jpan  ! *** for n=1,npan

kx=size(s,2)

call bounds(s,nrows=2)

!======================== start of intsch=1 section ====================
if(intsch==1)then

  do k=1,kx

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
          sextra(iproc)%a(iq) = ((1.-yyg)*((2.-yyg)* &
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
  do k=1,kx             ! A single main k loop uses cache better
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
          sextra(iproc)%a(iq) = ((1.-xxg)*((2.-xxg)*  &
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

integer k,kx
real, dimension(:,:), intent(inout) :: cou,cov,cow
real, dimension(ifull) :: vec1x,vec1y,vec1z,denb
real, dimension(ifull) :: vec2x,vec2y,vec2z,vecdot
real, dimension(ifull) :: vec3x,vec3y,vec3z,vdot1,vdot2
real*8, dimension(:,:), intent(in) :: x3d,y3d,z3d

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
real, dimension(:,:), intent(in) :: u,v
real, dimension(:,:), intent(out) :: uout,vout
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
real, dimension(:,:), intent(in) :: u,v
real, dimension(:,:), intent(out) :: uout,vout
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

subroutine getww(cou,cov,neta,dd,dz,ee,nw)

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
real, dimension(ifull,wlev) :: kku,kkv,dzdum
real, dimension(ifull), intent(in) :: neta,dd,ee
real, dimension(ifull) :: div

div=0.
do ii=1,wlev
  dzdum(:,ii)=dz(:,ii)*max(1.+neta/dd,0.01)
  kku(:,ii)=(cou(1:ifull,ii)/emu(1:ifull)-cou(iwu,ii)/emu(iwu))*em(1:ifull)*em(1:ifull)/ds
  kkv(:,ii)=(cov(1:ifull,ii)/emv(1:ifull)-cov(isv,ii)/emv(isv))*em(1:ifull)*em(1:ifull)/ds
  div=div+(kku(:,ii)+kkv(:,ii))*dzdum(:,ii)
end do
div=div/(neta(1:ifull)+dd(1:ifull))
! nw is at half levels with nw(:,1) at the bottom of level 1
! positive nw is moving downwards to the ocean bottom
nw(:,1)=(div-kku(:,1)-kkv(:,1))*dzdum(:,1)
do ii=2,wlev
  nw(:,ii)=nw(:,ii-1)+(div-kku(:,ii)-kkv(:,ii))*dzdum(:,ii)
end do
do ii=1,wlev
  nw(:,ii)=nw(:,ii)*ee
end do

return
end subroutine getww

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine performs vertical advection

subroutine mlovadv(dtin,ww,uu,vv,ss,tt,dep,dz,neta,dd,wtr)

use mlo

implicit none

include 'mpif.h'
include 'newmpar.h'

integer its,its_g,ii,l,iq,ierr
real, intent(in) :: dtin
real dtnew
real, dimension(ifull), intent(in) :: neta,dd
real, dimension(ifull,wlev), intent(in) :: ww,dep,dz
real, dimension(ifull,wlev), intent(inout) :: uu,vv,ss,tt
real, dimension(ifull,wlev) :: depdum,dzdum
logical, dimension(ifull), intent(in) :: wtr

do ii=1,wlev
  depdum(:,ii)=dep(:,ii)*max(1.+neta/dd,0.01)
  dzdum(:,ii)=dz(:,ii)*max(1.+neta/dd,0.01)
end do

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
  xx=delu(:,ii)+sign(1.E-10,delu(:,ii))
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

    ddn=dep(in(iq),:)
    isin=s(in(iq),:)
    ddne=dep(ine(iq),:)
    isine=s(ine(iq),:)
    dds=dep(is(iq),:)
    isis=s(is(iq),:)
    ddse=dep(ise(iq),:)
    isise=s(ise(iq),:)
    dde=dep(ie(iq),:)
    isie=s(ie(iq),:)
    dden=dep(ien(iq),:)
    isien=s(ien(iq),:)
    ddw=dep(iw(iq),:)
    isiw=s(iw(iq),:)
    ddwn=dep(iwn(iq),:)
    isiwn=s(iwn(iq),:)
    
    ! now search for other levels
    do ii=1,wlev
      sx=s(iq,ii)
      ddx=dep(iq,ii)
    
      call searchdelta(isin,ddn,wtr(in(iq)),ddx,sx,idn,sn)
      call searchdelta(isine,ddne,wtr(ine(iq)),ddx,sn,idne,sne)
      call searchdelta(isis,dds,wtr(is(iq)),ddx,sx,ids,ss)
      call searchdelta(isise,ddse,wtr(ise(iq)),ddx,ss,idse,sse)
      call searchdelta(isie,dde,wtr(ie(iq)),ddx,sx,ide,se)
      call searchdelta(isien,dden,wtr(ien(iq)),ddx,se,iden,sen)
      call searchdelta(isiw,ddw,wtr(iw(iq)),ddx,sx,idw,sw)
      call searchdelta(isiwn,ddwn,wtr(iwn(iq)),ddx,sw,idwn,swn)
    
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

subroutine searchdelta(s,dep,wtr,dd,sx,id,ss)

use mlo

implicit none

integer, intent(inout) :: id
integer jj,fnd
real, dimension(wlev), intent(in) :: s,dep
real, intent(in) :: dd,sx
real, intent(out) :: ss
real xp
logical, intent(in) :: wtr

ss=sx
if (.not.wtr) return

fnd=wlev-1
do jj=id,wlev-2
  if (dep(jj+1).gt.dd-0.01) then
    fnd=jj
    exit
  end if
end do

! located a valid level to interpolate
id=fnd
xp=(dd-dep(id))/(dep(id+1)-dep(id))
xp=max(min(xp,1.),0.)
ss=(1.-xp)*s(id)+xp*s(id+1)

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

subroutine mlosalfix

use map_m
use mlo
use soil_m

implicit none

include 'mpif.h'
include 'newmpar.h'
include 'parm.h'

integer iq,ii,ierr
real salsuml,salsumg,voll,volg,adj
real, dimension(ifull,wlev) :: sal,dz

if (salfix.eq.0) return
if (nud_sss.ne.0) return

do ii=1,wlev
  call mloexpdep(1,dz(:,ii),ii,0)
  call mloexport(1,sal(:,ii),ii,0)
end do
dz=max(dz,1.E-3/real(wlev))

voll=0.
salsuml=0.
do iq=1,ifull
  if (.not.land(iq)) then
    voll=voll+sum(dz(iq,:))/(em(iq)*em(iq))
    salsuml=salsuml+sum(sal(iq,:)*dz(iq,:))/(em(iq)*em(iq))
  end if
end do

call MPI_AllReduce(voll,volg,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
call MPI_AllReduce(salsuml,salsumg,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)

adj=34.72-salsumg/volg

do ii=1,wlev
  sal(:,ii)=max(sal(:,ii)+adj,0.)
  call mloimport(1,sal(:,ii),ii,0)
end do

return
end subroutine

end module mlodynamics
