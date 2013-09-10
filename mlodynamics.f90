! These subroutines handle dynamics for the Mixed-Layer-Ocean model

! - Horizontal diffusion
! - River routing
! - Ocean dynamics
! - Ice dynamics

! This module also links to mlo.f90 which solves for 1D physics of
! ocean and ice processes.  Currently the code assumes the hydrostatic
! approximation which is reasonably valid to 1km resolution.

! Ocean and sea-ice dynamics are based on the R-grid used by CCAM

module mlodynamics

implicit none

private
public mlodiffusion,mlorouter,mlohadv,mlodyninit,watbdy,salbdy,ipice
public oldu1,oldv1,oldu2,oldv2,gosig,gosigh,godsig,ocnsmag,ocneps,fixsal,fixheight
public nstagoffmlo,mstagf

real, dimension(:), allocatable, save :: watbdy,salbdy
real, dimension(:), allocatable, save :: ipice,ee,eeu,eev,dd,ddu,ddv,dfdyu,dfdxv
real, dimension(:), allocatable, save :: gosig,gosigh,godsig
real, dimension(:,:), allocatable, save :: oldu1,oldu2,oldv1,oldv2
real, dimension(:,:), allocatable, save :: stwgt
integer, save :: nstagoffmlo
integer, parameter :: usetide  =1    ! tidal forcing (0=off, 1=on)
integer, parameter :: icemode  =2    ! ice stress (0=free-drift, 1=incompressible, 2=cavitating)
integer, parameter :: basinmd  =3    ! basin mode (0=soil, 1=redistribute, 2=pile-up, 3=leak)
integer, parameter :: mstagf   =10   ! alternating staggering (0=off left, -1=off right, >0 alternating)
integer, parameter :: koff     =1    ! time split stagger relative to A-grid (koff=0) or C-grid (koff=1)
integer, parameter :: nf       =2    ! power for horizontal diffusion reduction factor
integer, parameter :: itnmax   =6    ! number of interations for staggering
integer, save      :: fixsal   =1    ! Conserve salinity (0=Usual, 1=Fixed average salinity at 34.72)
integer, save      :: fixheight=1    ! Conserve free surface height (0=Usual, 1=Fixed average Height at 0)
real, parameter :: rhosn      = 330.      ! density snow (kg m^-3)
real, parameter :: rhoic      = 900.      ! density ice  (kg m^-3)
real, parameter :: grav       = 9.80616   ! gravitational constant (m s^-2)
real, parameter :: delphi     = 150.      ! horizontal diffusion reduction factor gradient
real, save      :: ocnsmag    = 1.        ! horizontal diffusion (2. in Griffies (2000), 1.-1.4 in POM (Mellor 2004))
real, save      :: ocneps     = 0.1       ! semi-implicit off-centring term
real, parameter :: inflowbias = 10.       ! Bias height for inflow into ocean.  Levels above this will flow back onto land.
real, parameter :: maxicefrac = 0.99      ! Maximum ice fraction

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine initialises mlo dynamical arrays
!
subroutine mlodyninit

use cc_mpi
use indices_m
use map_m
use mlo
use soil_m

implicit none

include 'newmpar.h'
include 'parm.h'
include 'parmdyn.h'

integer ii,iq,ierr
real, dimension(ifull,0:wlev) :: dephl
real, dimension(ifull,wlev) :: dep,dz
real, dimension(ifull+iextra,3) :: dumx,dumy
real, dimension(ifull+iextra) :: ff
real, dimension(ifull) :: tnu,tsu,tev,twv,tee,tnn
real, dimension(wlev) :: sig,sigh,dsig
real, dimension(3*wlev) :: dumz,gdumz
logical, dimension(ifull+iextra) :: wtr

! staggering
nstagoffmlo=0

! river water height
allocate(watbdy(ifull+iextra),salbdy(ifull+iextra))
watbdy=0.
salbdy=0.

! prep land-sea mask
allocate(ee(ifull+iextra))
allocate(eeu(ifull+iextra),eev(ifull+iextra))
ee=0.
ff=0.
eeu=0.
eev=0. 
where(.not.land)
  ee(1:ifull)=1.
end where
where(land)
  ff(1:ifull)=1.
end where

! Calculate depth arrays (free suface term is included later)
allocate(dd(ifull+iextra))
allocate(ddu(ifull+iextra),ddv(ifull+iextra))
dd=0.
ddu=0.
ddv=0.
dep=0.
dz=0.
do ii=1,wlev
  call mloexpdep(0,dep(:,ii),ii,0)
  call mloexpdep(1,dz(:,ii),ii,0)
end do
dephl(:,0)=0.
dephl(:,1)=dz(:,1)
do ii=2,wlev
  dephl(:,ii)=dephl(:,ii-1)+dz(:,ii)
end do
dd(1:ifull)=dephl(:,wlev)
dd(1:ifull)=max(dd(1:ifull),1.E-8)

! update bounds values
dumx(1:ifull,1)=ee(1:ifull)
dumx(1:ifull,2)=dd(1:ifull)
dumx(1:ifull,3)=ff(1:ifull)
call bounds(dumx(:,1:3),nrows=2)
ee(ifull+1:ifull+iextra)=dumx(ifull+1:ifull+iextra,1)
dd(ifull+1:ifull+iextra)=dumx(ifull+1:ifull+iextra,2)
ff(ifull+1:ifull+iextra)=dumx(ifull+1:ifull+iextra,3)
wtr=abs(ee-1.)<0.5
where (ee>1.5.or.ee<=0.5)
  ee=0.
end where
where (ff>1.5.or.ff<=0.5)
  ff=0.
end where

! Precompute weights for calculating staggered gradients
allocate(stwgt(ifull,4))
stwgt=0.
where (wtr(in).and.wtr(ine).and.wtr(ien).and.wtr(ie).and.wtr(1:ifull))
  stwgt(:,1)=1.
end where
where (wtr(is).and.wtr(ise).and.wtr(ies).and.wtr(ie).and.wtr(1:ifull))
  stwgt(:,2)=1.
end where
stwgt(:,3)=stwgt(:,1)
where (wtr(iw).and.wtr(iwn).and.wtr(inw).and.wtr(in).and.wtr(1:ifull))
  stwgt(:,4)=1.
end where

if (abs(nmlo)>=3) then
  onedice=0 ! Turn off 1D ice model

  ! dynamics save arrays
  allocate(oldu1(ifull,wlev),oldv1(ifull,wlev))
  allocate(oldu2(ifull,wlev),oldv2(ifull,wlev))
  allocate(ipice(ifull+iextra))
  oldu1=0.
  oldv1=0.
  oldu2=oldu1
  oldv2=oldv1
  ipice=0.

  ! special gradient arrays
  allocate(dfdyu(ifull),dfdxv(ifull))
  tnu=0.5*(f(in)+f(ine))
  tee=0.5*(f(1:ifull)+f(ie))
  tsu=0.5*(f(is)+f(ise))
  tev=0.5*(f(ie)+f(ien))
  tnn=0.5*(f(1:ifull)+f(in))
  twv=0.5*(f(iw)+f(iwn))
  dfdyu=0.5*(stwgt(:,1)*(tnu-tee)+stwgt(:,2)*(tee-tsu))*emu(1:ifull)/ds
  dfdxv=0.5*(stwgt(:,3)*(tev-tnn)+stwgt(:,4)*(tnn-twv))*emv(1:ifull)/ds
end if

! update staggered bounds values
eeu(1:ifull)=ee(1:ifull)*ee(ie)
eev(1:ifull)=ee(1:ifull)*ee(in)
ddu(1:ifull)=0.5*(dd(1:ifull)+dd(ie))
ddv(1:ifull)=0.5*(dd(1:ifull)+dd(in))
where (eeu(1:ifull)==0.)
  ddu(1:ifull)=1.E-8
end where
where (eev(1:ifull)==0.)
  ddv(1:ifull)=1.E-8
end where
dumx(1:ifull,1)=eeu(1:ifull)
dumy(1:ifull,1)=eev(1:ifull)
dumx(1:ifull,2)=ddu(1:ifull)
dumy(1:ifull,2)=ddv(1:ifull)
call boundsuv(dumx(:,1:2),dumy(:,1:2),nrows=2)
eeu(ifull+1:ifull+iextra)=dumx(ifull+1:ifull+iextra,1)
eev(ifull+1:ifull+iextra)=dumy(ifull+1:ifull+iextra,1)
ddu(ifull+1:ifull+iextra)=dumx(ifull+1:ifull+iextra,2)
ddv(ifull+1:ifull+iextra)=dumy(ifull+1:ifull+iextra,2)

! sigma coordinates should be the same for all iq
allocate(gosig(wlev),gosigh(wlev),godsig(wlev))
sig=0.
sigh=0.
dsig=0.
do iq=1,ifull
  if (.not.land(iq)) then
    do ii=1,wlev
      sig(ii)=max(sig(ii),dep(iq,ii)/dd(iq))
      sigh(ii)=max(sigh(ii),dephl(iq,ii)/dd(iq))
      dsig(ii)=max(dsig(ii),dz(iq,ii)/dd(iq))
    end do
  end if
end do
do iq=1,ifull
  if (.not.land(iq)) then
    do ii=1,wlev
      if (abs(sig(ii)*dd(iq)-dep(iq,ii))>0.01.or.    &
          abs(sigh(ii)*dd(iq)-dephl(iq,ii))>0.01.or. &
          abs(dsig(ii)*dd(iq)-dz(iq,ii))>0.01) then
        write(6,*) "ERROR: MLO not configured for sigma levels"
        call ccmpi_abort(-1)
      end if
    end do
  end if
end do
dumz(1:wlev)=sig
dumz(wlev+1:2*wlev)=sigh
dumz(2*wlev+1:3*wlev)=dsig
call ccmpi_allreduce(dumz(1:3*wlev),gdumz(1:3*wlev),"max",comm_world)
gosig=gdumz(1:wlev)
gosigh=gdumz(wlev+1:2*wlev)
godsig=gdumz(2*wlev+1:3*wlev)

return
end subroutine mlodyninit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine processes horizontal diffusion, if mlohadv is
! not called.  Calling diffusion from mlodadv avoids additional
! unpacking and packing
subroutine mlodiffusion

use mlo

implicit none

include 'newmpar.h'

integer k
real, dimension(ifull,wlev) :: u,v,tt,ss
real, dimension(ifull) :: eta

! extract data from MLO
u=0.
v=0.
eta=0.
tt=0.
ss=0.
do k=1,wlev
  call mloexport(0,tt(:,k),k,0)
  call mloexport(1,ss(:,k),k,0)
  call mloexport(2,u(:,k),k,0)
  call mloexport(3,v(:,k),k,0)
end do
call mloexport(4,eta(1:ifull),0,0)

call mlodiffusion_main(u,v,u,v,eta,tt,ss)

return
end subroutine mlodiffusion

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine processes horizontal diffusion, based on Griffies (2000)
! and McGregor's hordifg.f routines for CCAM.
subroutine mlodiffusion_main(uauin,uavin,u,v,etain,tt,ss)

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

integer k,i,iq
real hdif,xp
real, dimension(ifull,wlev), intent(in) :: uauin,uavin,u,v,tt,ss
real, dimension(ifull), intent(in) :: etain
real, dimension(ifull+iextra,wlev) :: uau,uav
real, dimension(ifull+iextra,wlev) :: xfact,yfact
real, dimension(ifull+iextra,wlev+1) :: t_kh
real, dimension(ifull+iextra,wlev,3) :: duma
real, dimension(ifull) :: dudx,dvdx,dudy,dvdy
real, dimension(ifull,wlev) :: ft,fs,base,outu,outv
real, dimension(ifull+iextra) :: depadj,eta
real, dimension(ifull) :: tx_fact,ty_fact
real, dimension(ifull) :: cc,emi,nu,nv,nw

call start_log(waterdiff_begin)

! Define diffusion scale and grid spacing
hdif=dt*(ocnsmag/pi)**2
emi=1./(em(1:ifull)*em(1:ifull))

! extract data from MLO
do k=1,wlev
  uau(1:ifull,k)=uauin(:,k)*ee(1:ifull)
  uav(1:ifull,k)=uavin(:,k)*ee(1:ifull)
end do
call boundsuv_allvec(uau,uav)

! calculate diffusion following Smagorinsky
do k=1,wlev
  dudx=0.5*((uau(ieu,k)-uau(1:ifull,k))*emu(1:ifull)*eeu(1:ifull) &
           +(uau(1:ifull,k)-uau(iwu,k))*emu(iwu)*eeu(iwu))/ds
  dudy=0.5*((uau(inu,k)-uau(1:ifull,k))*emv(1:ifull)*eev(1:ifull) &
           +(uau(1:ifull,k)-uau(isu,k))*emv(isv)*eev(isv))/ds
  dvdx=0.5*((uav(iev,k)-uav(1:ifull,k))*emu(1:ifull)*eeu(1:ifull) &
           +(uav(1:ifull,k)-uav(iwv,k))*emu(iwu)*eeu(iwu))/ds
  dvdy=0.5*((uav(inv,k)-uav(1:ifull,k))*emv(1:ifull)*eev(1:ifull) &
           +(uav(1:ifull,k)-uav(isv,k))*emv(isv)*eev(isv))/ds

  cc=(dudx-dvdy)**2+(dudy+dvdx)**2
  t_kh(1:ifull,k)=sqrt(cc)*hdif*emi
end do
t_kh(1:ifull,wlev+1)=etain
call bounds(t_kh,nehalf=.true.)
eta(:)=t_kh(:,wlev+1)

! reduce diffusion errors where bathymetry gradients are strong
do k=1,wlev
  depadj=gosig(k)*max(dd+eta,minwater)
  tx_fact=1./(1.+(abs(depadj(ie)-depadj(1:ifull))/delphi)**nf)
  ty_fact=1./(1.+(abs(depadj(in)-depadj(1:ifull))/delphi)**nf)

  xfact(1:ifull,k)=0.5*(t_kh(1:ifull,k)+t_kh(ie,k))*tx_fact*eeu(1:ifull) ! reduction factor
  yfact(1:ifull,k)=0.5*(t_kh(1:ifull,k)+t_kh(in,k))*ty_fact*eev(1:ifull) ! reduction factor
end do
call boundsuv(xfact,yfact,stag=-9)

! Laplacian diffusion terms (closure #1)
do k=1,wlev
  duma(1:ifull,k,1)=ax(1:ifull)*u(1:ifull,k)+bx(1:ifull)*v(1:ifull,k)
  duma(1:ifull,k,2)=ay(1:ifull)*u(1:ifull,k)+by(1:ifull)*v(1:ifull,k)
  duma(1:ifull,k,3)=az(1:ifull)*u(1:ifull,k)+bz(1:ifull)*v(1:ifull,k)
end do
call bounds(duma(:,:,1:3))

! allow drag on momentum along coastlines (but not for scalars, see below)
do k=1,wlev

  base(:,k)=emi+xfact(1:ifull,k)+xfact(iwu,k)+yfact(1:ifull,k)+yfact(isv,k)

  nu = ( duma(1:ifull,k,1)*emi +                      &
         xfact(1:ifull,k)*duma(ie,k,1) +              &
         xfact(iwu,k)*duma(iw,k,1) +                  &
         yfact(1:ifull,k)*duma(in,k,1) +              &
         yfact(isv,k)*duma(is,k,1) ) / base(:,k)

  nv = ( duma(1:ifull,k,2)*emi +                      &
         xfact(1:ifull,k)*duma(ie,k,2) +              &
         xfact(iwu,k)*duma(iw,k,2) +                  &
         yfact(1:ifull,k)*duma(in,k,2) +              &
         yfact(isv,k)*duma(is,k,2) ) / base(:,k)

  nw = ( duma(1:ifull,k,3)*emi +                      &
         xfact(1:ifull,k)*duma(ie,k,3) +              &
         xfact(iwu,k)*duma(iw,k,3) +                  &
         yfact(1:ifull,k)*duma(in,k,3) +              &
         yfact(isv,k)*duma(is,k,3) ) / base(:,k)

  outu(:,k)=ax(1:ifull)*nu+ay(1:ifull)*nv+az(1:ifull)*nw
  outv(:,k)=bx(1:ifull)*nu+by(1:ifull)*nv+bz(1:ifull)*nw

end do

! Laplacian diffusion and viscosity terms (closure #2)
! call boundsuv(u,v,allvec=.true.)
!do k=1,wlev
!
!  base(:,k)=emi+xfact(1:ifull,k)+xfact(iwu,k)+yfact(1:ifull,k)+yfact(isv,k)
!
!  outu(:,k)=(u(1:ifull,k)*emi+2.*xfact(1:ifull,k)*u(ieu,k)+2.*xfact(iwu,k)*u(iwu,k) &
!    +yfact(1:ifull,k)*u(inu,k)+yfact(isv,k)*u(isu,k)                         &
!    +(yfact(1:ifull,k)-yfact(isv,k))*0.5*(v(iev,k)-v(iwv,k))                 &
!    +t_kh(1:ifull,k)*0.5*(v(inv,k)+v(iev,k)-v(isv,k)-v(iwv,k)))              &
!    /(emi+2.*xfact(1:ifull,k)+2.*xfact(iwu,k)+yfact(1:ifull,k)+yfact(isv,k))
!  outv(:,k)=(v(1:ifull,k)*emi+2.*yfact(1:ifull,k)*v(inv,k)+2.*yfact(isv,k)*v(isv,k) &
!    +xfact(1:ifull,k)*v(iev,k)+xfact(iwu,k)*v(iwv,k)                         &
!    +(xfact(1:ifull,k)-xfact(iwu,k))*0.5*(u(inu,k)-u(isu,k))                 &
!    +t_kh(1:ifull,k)*0.5*(u(inu,k)+u(ieu,k)-u(isu,k)-u(iwu,k)))              &
!    /(emi+2.*yfact(1:ifull,k)+2.*yfact(isv,k)+xfact(1:ifull,k)+xfact(iwu,k))
!
!end do

! Potential temperature
duma(1:ifull,:,1)=tt-290.
duma(1:ifull,:,2)=ss-34.72
call bounds(duma(:,:,1:2))
do k=1,wlev
  ft(:,k) = ( duma(1:ifull,k,1)*emi +                      &
              xfact(1:ifull,k)*duma(ie,k,1) +              &
              xfact(iwu,k)*duma(iw,k,1) +                  &
              yfact(1:ifull,k)*duma(in,k,1) +              &
              yfact(isv,k)*duma(is,k,1) ) / base(:,k)
  fs(:,k) = ( duma(1:ifull,k,2)*emi +                      &
              xfact(1:ifull,k)*duma(ie,k,2) +              &
              xfact(iwu,k)*duma(iw,k,2) +                  &
              yfact(1:ifull,k)*duma(in,k,2) +              &
              yfact(isv,k)*duma(is,k,2) ) / base(:,k)
end do
ft=ft+290.
fs=max(fs+34.72,0.)

do k=1,wlev
  call mloimport(0,ft(:,k),k,0)
  call mloimport(1,fs(:,k),k,0)
  call mloimport(2,outu(:,k),k,0)
  call mloimport(3,outv(:,k),k,0)
end do

call end_log(waterdiff_end)

return
end subroutine mlodiffusion_main

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

integer i,ii,iq,ierr,k
integer nit
integer, dimension(ifull,4) :: xp
real, dimension(ifull,wlev) :: sallvl
real, dimension(ifull+iextra) :: neta,netflx,cc
real, dimension(ifull+iextra) :: newwat
real, dimension(ifull+iextra,2) :: dum
real, dimension(ifull) :: newsal,cover
real, dimension(ifull) :: deta,sal,salin,depdum
real, dimension(ifull,4) :: idp,slope,mslope,vel,flow
real, dimension(ifull,4) :: fta,ftb,ftx,fty
real, dimension(2) :: dumb,gdumb
real xx,yy,ll,lssum,gssum,lwsum,gwsum,netf,netg,rate

! To speed up the code, we use a (semi-)implicit solution rather than an iterative approach
! This avoids additional MPI calls.

! Note that unlike Miller et al (1994), this scheme currently has multiple outflows from a given
! grid box (i.e., following Gordon in Mk3.5).  So if the slope of the grid box is divergent, then
! on average water leaves the grid box, whereas if the slope of the grid box is convergent, then
! on average water enters the grid box.

! This version supports salinity in rivers and overflow from lakes

! MJT notes - Usually for MPI we would update the water level on the boundaries (MPI), estimate
! the fluxes between grid boxes, then update the fluxes on the boundaries (MPI), then update the
! new water levels at t+1.  However, in this case we already know the fluxes due to the slope
! (avoiding a MPI) so the second MPI is simply to constrain water avaliability so that the net amount
! of water is conserved (i.e., netvel in the code below).

! Currently fluxes are explicit.  May need iterative loop for implicit solver.

! setup indices and grid spacing
xp(:,1)=in
xp(:,2)=ie
xp(:,3)=is
xp(:,4)=iw
idp(:,1)=emv(1:ifull)/ds
idp(:,2)=emu(1:ifull)/ds
idp(:,3)=emv(isv)/ds
idp(:,4)=emu(iwu)/ds

! update bounds and initial conditions
neta=0.
salin=0.
sallvl=0.
call mloexport(4,neta(1:ifull),0,0)
neta(1:ifull)=neta(1:ifull)-inflowbias*ee(1:ifull)
do ii=1,wlev
  call mloexport(1,sallvl(:,ii),ii,0)
  ! could modify the following to only operate above neta=0
  salin=salin+godsig(ii)*sallvl(:,ii)
end do
where (ee(1:ifull)>0.5)
  depdum=max(neta(1:ifull)+dd(1:ifull),minwater)
  ! collect any left over water over ocean
  salin=(salin*depdum+salbdy(1:ifull)*watbdy(1:ifull)*0.001)/(depdum+watbdy(1:ifull)*0.001)
  neta(1:ifull)=neta(1:ifull)+0.001*watbdy(1:ifull)
  ! move water from neta to watbdy for transport
  deta=1000.*max(neta(1:ifull),0.)
  salbdy(1:ifull)=salin
  watbdy(1:ifull)=deta
  neta(1:ifull)=min(neta(1:ifull),0.)
elsewhere
  deta=0.
  depdum=0.
end where
salbdy(1:ifull)=salbdy(1:ifull)*watbdy(1:ifull) ! rescale salinity to PSU*mm for advection

! update boundaries
dum(1:ifull,1)=watbdy(1:ifull)
dum(1:ifull,2)=salbdy(1:ifull)
call bounds(dum(:,1:2))
watbdy(ifull+1:ifull+iextra)=dum(ifull+1:ifull+iextra,1)
salbdy(ifull+1:ifull+iextra)=dum(ifull+1:ifull+iextra,2)
newwat(1:ifull+iextra)=watbdy(1:ifull+iextra)
newsal=salbdy(1:ifull)

! predictor-corrector
do nit=1,2
  
  if (nit==2) call bounds(newwat)

  ! calculate slopes
  ! Currently this is has an explicit dependence on watbdy
  do i=1,4
    where ((ee(1:ifull)*ee(xp(:,i)))>0.5)
      slope(:,i)=0. ! no orographic slope within ocean bounds
    elsewhere
      !slope(:,i)=(zs(1:ifull)-zs(xp(:,i)))/(grav*dp(:,i))         ! basic
      slope(:,i)=(zs(1:ifull)/grav+0.001*newwat(1:ifull) &
                 -zs(xp(:,i))/grav-0.001*newwat(xp(:,i)))*idp(:,i) ! flood
    end where
  end do

  newwat(1:ifull)=watbdy(1:ifull)

  ! Basic expression

  ! m = mass/area
  ! flow = m * vel / dx
  ! m(t+1)-m(t) = dt*sum(inflow)-dt*sum(outflow)

  ! outflow
  mslope=max(slope,0.)
  vel=min(0.35*sqrt(mslope/0.00005),5.) ! from Miller et al (1994)
  ! compute net outgoing flux for a grid box so that total water is conserved
  do i=1,4
    fta(:,i)=-dt*vel(:,i)*idp(:,i)     ! outgoing flux
  end do
  netflx(1:ifull)=sum(abs(fta),2)
  call bounds(netflx)
  
  ! water outflow
  do i=1,4
    where (netflx(1:ifull)>1.E-10)
      ftx(:,i)=-fta(:,i)/netflx(1:ifull) ! max fraction of total outgoing flux
      flow(:,i)=watbdy(1:ifull)*min(fta(:,i),ftx(:,i)) ! (kg/m^2)
    elsewhere
      flow(:,i)=0.
    end where
  end do
  newwat(1:ifull)=newwat(1:ifull)+sum(flow,2)

  ! inflow
  mslope=max(-slope,0.)
  vel=min(0.35*sqrt(mslope/0.00005),5.) ! from Miller et al (1994)

  ! water inflow
  do i=1,4
    ftb(:,i)=dt*vel(:,i)*idp(:,i)     ! incomming flux
    where (netflx(xp(:,i))>1.E-10)
      fty(:,i)=ftb(:,i)/netflx(xp(:,i)) ! max fraction of flux from outgoing cel
      flow(:,i)=watbdy(xp(:,i))*min(ftb(:,i),fty(:,i)) ! (kg/m^2)
      flow(:,i)=flow(:,i)*em(1:ifull)*em(1:ifull)/(em(xp(:,i))*em(xp(:,i))) ! change in gridbox area
    elsewhere
      flow(:,i)=0.
    end where
  end do
  newwat(1:ifull)=newwat(1:ifull)+sum(flow,2)

end do

! salinity outflow
do i=1,4
  where (netflx(1:ifull)>1.E-10)
    flow(:,i)=salbdy(1:ifull)*min(fta(:,i),ftx(:,i))
  elsewhere
    flow(:,i)=0.
  end where
end do
newsal=newsal+sum(flow,2)
! salinity inflow
do i=1,4
  where (netflx(xp(:,i))>1.E-10)
    flow(:,i)=salbdy(xp(:,i))*min(ftb(:,i),fty(:,i))
    flow(:,i)=flow(:,i)*em(1:ifull)*em(1:ifull)/(em(xp(:,i))*em(xp(:,i))) ! change in gridbox area
  elsewhere
    flow(:,i)=0.
  end where
end do
newsal=newsal+sum(flow,2)

watbdy(1:ifull)=max(newwat,0.)
salbdy(1:ifull)=max(newsal,0.)

! estimate grid box area covered by water
cover=min(0.001*watbdy(1:ifull)/minwater,1.)
rate=min(dt/(8.*3600.),1.) ! MJT suggestion
  
! basin
select case(basinmd)
  case(0)
    ! add water to soil moisture 
    if (nsib==6.or.nsib==7) then
      ! CABLE
      do iq=1,ifull
        if (all(slope(iq,:)<-1.E-10).and.land(iq)) then
          ! runoff is inserted into soil moisture since CCAM has discrete land and sea points
          xx=watbdy(iq)
          ll=cover(iq)
          call cableinflow(iq,xx,ll,rate)
          newwat(iq)=newwat(iq)+(xx-watbdy(iq))*(1.-sigmu(iq))
        end if
      end do
    else
      ! Standard land surface model
      do iq=1,ifull
        if (all(slope(iq,:)<-1.E-10).and.land(iq)) then
          ! runoff is inserted into soil moisture since CCAM has discrete land and sea points
          xx=watbdy(iq)
          do k=1,ms
            ll=max(sfc(isoilm(iq))-wb(iq,k),0.)*1000.*zse(k)
            ll=ll*rate*cover(iq)
            yy=min(xx,ll)
            wb(iq,k)=wb(iq,k)+yy/(1000.*zse(k))
            xx=max(xx-yy,0.)
          end do
          newwat(iq)=newwat(iq)+(xx-watbdy(iq))*(1.-sigmu(iq))
        end if
      end do
    end if
  case(1)
    ! distribute water over ocean (used when soil column is too shallow)
    lssum=0.
    lwsum=0.
    do iq=1,ifull
      if (all(slope(iq,:)<-1.E-10).and.land(iq)) then
        lssum=lssum+newwat(iq)/(em(iq)*em(iq)) ! kg of water / (ds*ds)
        newwat(iq)=0.
      elseif (.not.land(iq)) then
        lwsum=lwsum+1./(em(iq)*em(iq))   
      end if
    end do
    dumb(1)=lssum
    dumb(2)=lwsum
    call ccmpi_allreduce(dumb(1:2),gdumb(1:2),"sum",comm_world)
    gssum=gdumb(1)
    gwsum=gdumb(2)
    xx=gssum/max(gwsum,0.1)
    where (.not.land)
      neta(1:ifull)=neta(1:ifull)+xx*0.001
    end where
    call mloimport(4,neta(1:ifull),0,0)
  case(2)
    ! pile-up water
  case(3)
    ! leak
    if (nsib==6.or.nsib==7) then
      ! CABLE
      do iq=1,ifull
        if (land(iq).and.watbdy(iq)>0.) then
          xx=watbdy(iq)
          ll=cover(iq)
          call cableinflow(iq,xx,ll,rate)
          newwat(iq)=newwat(iq)+(xx-watbdy(iq))*(1.-sigmu(iq))
        end if
      end do
    else
      ! Standard land surface model
      do iq=1,ifull
        if (land(iq).and.watbdy(iq)>0.) then
          xx=watbdy(iq)
          do k=1,ms
            ll=max(sfc(isoilm(iq))-wb(iq,k),0.)*1000.*zse(k)
            ll=ll*rate*cover(iq)
            yy=min(xx,ll)
            wb(iq,k)=wb(iq,k)+yy/(1000.*zse(k))
            xx=max(xx-yy,0.)
          end do
          newwat(iq)=newwat(iq)+(xx-watbdy(iq))*(1.-sigmu(iq))
        end if
      end do
    end if
  case default
    write(6,*) "ERROR: Unsupported basinmd ",basinmd
    stop
end select

watbdy(1:ifull)=max(newwat,0.)
salbdy(1:ifull)=max(newsal,0.)

! case when water is absorbed by land-surface, but there is
! non-zero salinity left behind.
if (any(salbdy(1:ifull)>0..and.watbdy(1:ifull)==0.)) then
  write(6,*) "WARN: Patch river salinity"
end if

! update ocean
where (ee(1:ifull)>0.5)
  depdum=max(dd(1:ifull)+neta(1:ifull),minwater)
  ! salbdy is already multiplied by watbdy
  sal=(salin*depdum+0.001*salbdy(1:ifull))/(depdum+0.001*watbdy(1:ifull))
  neta(1:ifull)=neta(1:ifull)+0.001*watbdy(1:ifull)
  watbdy(1:ifull)=0.
  salbdy(1:ifull)=sal   ! ocean default
elsewhere (watbdy(1:ifull)>1.E-10)
  salbdy(1:ifull)=salbdy(1:ifull)/watbdy(1:ifull) ! rescale salinity back to PSU
  sal=0.
elsewhere
  salbdy(1:ifull)=0.    ! land default
  sal=0.
end where

! import ocean data
neta(1:ifull)=neta(1:ifull)+inflowbias*ee(1:ifull)
call mloimport(4,neta(1:ifull),0,0)
do ii=1,wlev
  where (ee(1:ifull)>0.5.and.salin>0.)
    ! rescale salinity profile
    sallvl(:,ii)=sallvl(:,ii)*sal/salin
  elsewhere (ee(1:ifull)>0.5)
    ! set new uniform salinity profile
    sallvl(:,ii)=sal
  end where
  call mloimport(1,sallvl(:,ii),ii,0)
end do

return
end subroutine mlorouter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine implements some basic hydrostatic dynamics for the
! ocean and ice.  The ocean component employs the R-grid design used
! in CCAM semi-Lagragain dynamics, but uses sigma z vertical
! coordinates.  The staggered/unstaggered pivoting has been modified
! for the land-sea boundary.  Sea-ice is advected using an upwind
! scheme.  Internal sea-ice pressure follows a cavitating fluid
! apporximation.
subroutine mlohadv

use arrays_m
use cc_mpi
use infile
use indices_m
use latlong_m
use map_m
use mlo
use nharrs_m, only : lrestart
use soil_m
use soilsnow_m
use vecsuv_m

implicit none

include 'newmpar.h'
include 'const_phys.h'
include 'parm.h'
include 'parmdyn.h'

integer leap
common/leap_yr/leap  ! 1 to allow leap years

integer iq,ll,ii,ierr,totits
integer jyear,jmonth,jday,jhour,jmin,mins
integer tyear,jstart,iip,iim
integer itc
integer, dimension(ifull,wlev) :: nface
real maxglobseta,maxglobip
real delpos,delneg,alph_p,fjd
real, dimension(ifull+iextra) :: neta,pice,imass,xodum
real, dimension(ifull+iextra) :: nfracice,ndic,ndsn,nsto,niu,niv,nis,tiu,tiv
real, dimension(ifull+iextra) :: snu,sou,spu,squ,ssu,snv,sov,spv,sqv,ssv
real, dimension(ifull+iextra) :: ibu,ibv,icu,icv,idu,idv,spnet,oeu,oev,tide
real, dimension(ifull+iextra) :: ipmax
real, dimension(ifull) :: i_u,i_v,i_sto,i_sal,rhobaru,rhobarv,ndum
real, dimension(ifull) :: pdiv,qdiv,sdiv,div,odiv,w_e
real, dimension(ifull) :: pdivb,qdivb,sdivb,odivb,xps
real, dimension(ifull) :: tnu,tsu,tev,twv,tee,tnn,rhou,rhov
real, dimension(ifull) :: dpsdxu,dpsdyu,dpsdxv,dpsdyv
real, dimension(ifull) :: dttdxu,dttdyu,dttdxv,dttdyv
real, dimension(ifull) :: detadxu,detadyu,detadxv,detadyv
real, dimension(ifull) :: dipdxu,dipdyu,dipdxv,dipdyv
real, dimension(ifull) :: au,bu,cu,av,bv,cv,odum
real, dimension(ifull) :: imu,imv
real, dimension(ifull) :: sue,suw,svn,svs,snuw,snvs
real, dimension(ifull) :: pue,puw,pvn,pvs
real, dimension(ifull) :: que,quw,qvn,qvs
real, dimension(ifull) :: gamm,piceu,picev,tideu,tidev,ipiceu,ipicev
real, dimension(ifull) :: dumf,dumg
real, dimension(ifull+iextra,wlev,3) :: cou
real, dimension(ifull+iextra,wlev+1) :: eou,eov
real, dimension(ifull+iextra,wlev) :: nu,nv,nt,ns,mps,xdzdum
real, dimension(ifull+iextra,wlev) :: rhobar,rho,dalpha,dbeta
real, dimension(ifull+iextra,wlev) :: ccu,ccv
real, dimension(ifull+iextra,3*wlev) :: dume
real, dimension(ifull+iextra,10) :: dumc,dumd
real, dimension(ifull+iextra,4) :: nit
real, dimension(ifull,4) :: i_it
real, dimension(ifull,wlev+1) :: tau,tav,ttau,ttav
real, dimension(ifull,wlev) :: w_u,w_v,w_t,w_s,dum
real, dimension(ifull,wlev) :: nuh,nvh,xg,yg,uau,uav
real, dimension(ifull,wlev) :: kku,llu,mmu,nnu,oou,ppu
real, dimension(ifull,wlev) :: kkv,llv,mmv,nnv,oov,ppv
real, dimension(ifull,wlev) :: drhobardxu,drhobardyu,drhobardxv,drhobardyv
real, dimension(ifull,wlev) :: depdum,dzdum
real, dimension(ifull,wlev) :: dumq,dumt,dums,dumr,duma,dumb
real, dimension(ifull,0:wlev) :: nw
real*8, dimension(ifull,wlev) :: x3d,y3d,z3d
logical, dimension(ifull+iextra) :: wtr
logical lleap

integer, parameter :: nxtrrho    = 1   ! Estimate rho at t+1 (0=off, 1=on)
real, parameter :: tol    = 5.E-5      ! Tolerance for SOR solver (water)
real, parameter :: itol   = 2.E1       ! Tolerance for SOR solver (ice)

! new z levels for including free surface eta (effectively sigma-depth levels)
! newz=-eta+oldz*(1+eta/maxdepth)
! newdz=olddz*(1+eta/maxdepth)
! where 0<=oldz<=maxdepth and -eta<=newz<=maxdepth
! depth below free suface is newz+eta=oldz*(1+eta/maxdepth)

! We use a modified form of JLM's trick where:
!   dP/dx+f dP/dy = dP/dx + d(fP)/dy - P df/dy
!   dP/dy-f dP/dx = dP/dy - d(fP)/dy + P df/dx
! which produces a 5-point stencil when solving for
! neta and ip.

! Instead of packing ocean points, we keep the original
! grid and define a land/sea mask.  This allows us to
! use the same index and map factor arrays from the
! atmospheric dynamical core.

call start_log(watermisc_begin)

if (myid==0.and.nmaxpr==1) then
  write(6,*) "mlohadv: Start"
end if

! Define land/sea mask
wtr=ee>0.5

! Default values
w_t=293.16
w_s=34.72
w_u=0.
w_v=0.
w_e=0.
i_it=273.16
i_sto=0.
i_u=0.
i_v=0.
i_sal=0.
rho=0.
nw=0.
pice=0.
imass=0.
tide=0.

! IMPORT WATER AND ICE DATA -----------------------------------------
if (myid==0.and.nmaxpr==1) then
  write(6,*) "mlohadv: Import"
end if
do ii=1,wlev
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
call mloexpice(i_sal,11,0)

! estimate tidal forcing (assumes leap days)
if (myid==0.and.nmaxpr==1) then
  write(6,*) "mlohadv: Tides"
end if
if (usetide==1) then
  call getzinp(fjd,jyear,jmonth,jday,jhour,jmin,mins,allleap=.true.)
  jstart=0
  if (jyear>1900) then
    do tyear=1900,jyear-1
      call mloleap(tyear,lleap)
      if (lleap) jstart=jstart+1
      jstart=jstart+365
    end do
  else if (jyear<1900) then
    do tyear=1899,jyear,-1
      call mloleap(tyear,lleap)
      if (lleap) jstart=jstart-1
      jstart=jstart-365
    end do
  end if
  mins=mins+720 ! base time is 12Z 31 Dec 1899
  call mlotide(tide(1:ifull),rlongg,rlatt,mins,jstart)
end if

! initialise t+1 variables with t data
neta(1:ifull)=w_e
nu(1:ifull,:)=w_u
nv(1:ifull,:)=w_v
nt(1:ifull,:)=w_t
ns(1:ifull,:)=w_s
nfracice(ifull+1:ifull+iextra)=0.
ndic(ifull+1:ifull+iextra)=0.
ndsn(ifull+1:ifull+iextra)=0.
where (wtr(1:ifull))
  nfracice(1:ifull)=fracice
  ndic(1:ifull)=sicedep
  ndsn(1:ifull)=snowd*0.001
elsewhere
  nfracice(1:ifull)=0.
  ndic(1:ifull)=0.
  ndsn(1:ifull)=0.
end where
nit(1:ifull,:)=i_it
nsto(1:ifull)=i_sto
niu(1:ifull)=i_u
niv(1:ifull)=i_v
nis(1:ifull)=i_sal

totits=0

! surface pressure gradients (including ice)
! (assume ice velocity is 'slow' compared to 'fast' change in neta)
imass(1:ifull)=ndic(1:ifull)*rhoic+ndsn(1:ifull)*rhosn            ! ice mass per unit area (kg/m^2), unstaggered at time t
pice(1:ifull)=ps(1:ifull)+grav*nfracice(1:ifull)*imass(1:ifull)   ! pressure due to atmosphere and ice at top of water column (unstaggered at t)

! Limit minimum ice mass for ice velocity calculation.  Hence we can estimate the ice velocity at
! grid points where the ice is not yet present.  
imass(1:ifull)=max(imass(1:ifull),10.)

! maximum pressure for cavitating fluid
select case(icemode)
  case(2)
    ! cavitating fluid
    ipmax(1:ifull)=27500.*ndic(1:ifull)*exp(-20.*(1.-nfracice(1:ifull)))*ee(1:ifull)
  case(1)
    ! incompressible fluid
    ipmax(1:ifull)=9.E9*ee(1:ifull)
  case DEFAULT
    ! free drift
    ipmax(1:ifull)=0.
end select

! update scalar bounds and various gradients
dumc(1:ifull,1)=neta(1:ifull)
dumc(1:ifull,2)=pice(1:ifull)
dumc(1:ifull,3)=tide(1:ifull)
dumc(1:ifull,4)=imass(1:ifull)
dumc(1:ifull,5)=ipmax(1:ifull)
call bounds(dumc(:,1:5),corner=.true.)
neta(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,1)
pice(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,2)
tide(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,3)
imass(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,4)
ipmax(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,5)
! surface pressure
tnu=0.5*(pice(in)*f(in)+pice(ine)*f(ine))
tee=0.5*(pice(1:ifull)*f(1:ifull)+pice(ie)*f(ie))
tsu=0.5*(pice(is)*f(is)+pice(ise)*f(ise))
tev=0.5*(pice(ie)*f(ie)+pice(ien)*f(ien))
tnn=0.5*(pice(1:ifull)*f(1:ifull)+pice(in)*f(in))
twv=0.5*(pice(iw)*f(iw)+pice(iwn)*f(iwn))
dpsdxu=(pice(ie)-pice(1:ifull))*emu(1:ifull)/ds ! staggered at time t
dpsdyu=0.5*(stwgt(:,1)*(tnu-tee)+stwgt(:,2)*(tee-tsu))*emu(1:ifull)/ds
dpsdxv=0.5*(stwgt(:,3)*(tev-tnn)+stwgt(:,4)*(tnn-twv))*emv(1:ifull)/ds
dpsdyv=(pice(in)-pice(1:ifull))*emv(1:ifull)/ds
piceu=0.5*(pice(1:ifull)+pice(ie))
picev=0.5*(pice(1:ifull)+pice(in))
! tides
tnu=0.5*(tide(in)*f(in)+tide(ine)*f(ine))
tee=0.5*(tide(1:ifull)*f(1:ifull)+tide(ie)*f(ie))
tsu=0.5*(tide(is)*f(is)+tide(ise)*f(ise))
tev=0.5*(tide(ie)*f(ie)+tide(ien)*f(ien))
tnn=0.5*(tide(1:ifull)*f(1:ifull)+tide(in)*f(in))
twv=0.5*(tide(iw)*f(iw)+tide(iwn)*f(iwn))
dttdxu=(tide(ie)-tide(1:ifull))*emu(1:ifull)/ds ! staggered
dttdyu=0.5*(stwgt(:,1)*(tnu-tee)+stwgt(:,2)*(tee-tsu))*emu(1:ifull)/ds
dttdxv=0.5*(stwgt(:,3)*(tev-tnn)+stwgt(:,4)*(tnn-twv))*emv(1:ifull)/ds
dttdyv=(tide(in)-tide(1:ifull))*emv(1:ifull)/ds
tideu=0.5*(tide(1:ifull)+tide(ie))
tidev=0.5*(tide(1:ifull)+tide(in))

call end_log(watermisc_end)
call start_log(waterdeps_begin)

! ADVECT WATER AND ICE ----------------------------------------------
! Water currents are advected using semi-Lagrangian advection
! based on McGregor's CCAM advection routines.
! Velocity is set to zero at ocean boundaries.

if (myid==0.and.nmaxpr==1) then
  write(6,*) "mlohadv: Departure points"
end if

! save arrays
if (ktau==1.and..not.lrestart) then
  oldu1=nu(1:ifull,:)
  oldv1=nv(1:ifull,:)
  oldu2=nu(1:ifull,:)
  oldv2=nv(1:ifull,:)
  ipice=0. ! ice free drift solution
end if

! Estimate currents at t+1/2 for semi-Lagrangian advection
do ii=1,wlev
  nuh(:,ii)=(15.*nu(1:ifull,ii)-10.*oldu1(:,ii)+3.*oldu2(:,ii))*ee(1:ifull)/8. ! U at t+1/2
  nvh(:,ii)=(15.*nv(1:ifull,ii)-10.*oldv1(:,ii)+3.*oldv2(:,ii))*ee(1:ifull)/8. ! V at t+1/2
end do

! Calculate depature points
call mlodeps(dt,nuh,nvh,nface,xg,yg,x3d,y3d,z3d,wtr)

oldu2=oldu1
oldv2=oldv1
oldu1=nu(1:ifull,:)
oldv1=nv(1:ifull,:)

! Calculate adjusted depths and thicknesses
xodum=max(dd(:)+neta(:),minwater)
do ii=1,wlev
  depdum(:,ii)=gosig(ii)*xodum(1:ifull)
  xdzdum(:,ii)=godsig(ii)*xodum(:)
  dzdum(:,ii)=xdzdum(1:ifull,ii)
end do

call end_log(waterdeps_end)
call start_log(watereos_begin)

if (myid==0.and.nmaxpr==1) then
  write(6,*) "mlohadv: density EOS 1"
end if

! Calculate normalised density rhobar (unstaggered at time t)
! (Assume free surface correction is small so that changes in the compression 
! effect due to neta can be neglected.  Consequently, the neta dependence is 
! separable in the iterative loop)
cou(1:ifull,:,1)=nt(1:ifull,:)
cou(1:ifull,:,2)=ns(1:ifull,:)
call bounds(cou(:,:,1:2),corner=.true.)
nt(ifull+1:ifull+iextra,:)=cou(ifull+1:ifull+iextra,:,1)
ns(ifull+1:ifull+iextra,:)=cou(ifull+1:ifull+iextra,:,2)
call mloexpdensity(rho,dalpha,dbeta,nt,ns,xdzdum,pice,0)
rhobar(:,1)=(rho(:,1)-1030.)*godsig(1)
do ii=2,wlev
  rhobar(:,ii)=rhobar(:,ii-1)+(rho(:,ii)-1030.)*godsig(ii)
end do
do ii=1,wlev
  rhobar(:,ii)=rhobar(:,ii)/gosigh(ii)+1030.
end do

! Calculate normalised density gradients
! method 2: Use potential temperature and salinity Jacobians (see Shchepetkin and McWilliams 2003)
call tsjacobi(nt,ns,dalpha,dbeta,neta,drhobardxu,drhobardyu,drhobardxv,drhobardyv)
drhobardxu(:,1)=drhobardxu(:,1)*godsig(1)
drhobardxv(:,1)=drhobardxv(:,1)*godsig(1)
drhobardyu(:,1)=drhobardyu(:,1)*godsig(1)
drhobardyv(:,1)=drhobardyv(:,1)*godsig(1)
do ii=2,wlev
  drhobardxu(:,ii)=drhobardxu(:,ii-1)+drhobardxu(:,ii)*godsig(ii)
  drhobardxv(:,ii)=drhobardxv(:,ii-1)+drhobardxv(:,ii)*godsig(ii)
  drhobardyu(:,ii)=drhobardyu(:,ii-1)+drhobardyu(:,ii)*godsig(ii)
  drhobardyv(:,ii)=drhobardyv(:,ii-1)+drhobardyv(:,ii)*godsig(ii)
end do
do ii=1,wlev
  drhobardxu(:,ii)=drhobardxu(:,ii)/gosigh(ii)
  drhobardxv(:,ii)=drhobardxv(:,ii)/gosigh(ii)
  drhobardyu(:,ii)=drhobardyu(:,ii)/gosigh(ii)
  drhobardyv(:,ii)=drhobardyv(:,ii)/gosigh(ii)
end do

call end_log(watereos_end)
call start_log(waterhadv_begin)

if (myid==0.and.nmaxpr==1) then
  write(6,*) "mlohadv: continuity equation"
end if

! Advect continuity equation to tstar
! Calculate velocity on C-grid for consistancy with iterative free surface calculation
! nw is at half levels with nw(:,1) at the bottom of level 1
! positive nw is moving downwards to the ocean bottom
! For now, assume Boussinesq fluid and treat density in continuity equation as constant
! true vertical velocity = nw-u*((1-sig)*deta/dx-sig*d(dd)/dx)-v*((1-sig)*deta/dy-sig*d(dd)/dy)-(1-sig)*deta/dt
dumt=nu(1:ifull,:)
dums=nv(1:ifull,:)
call mlostaguv(dumt,dums,duma,dumb)
eou(1:ifull,1:wlev)=duma
eov(1:ifull,1:wlev)=dumb
! surface height at staggered coordinate
oeu(1:ifull)=0.5*(neta(1:ifull)+neta(ie))*eeu(1:ifull) ! height at staggered coordinate
oev(1:ifull)=0.5*(neta(1:ifull)+neta(in))*eev(1:ifull) ! height at staggered coordinate
eou(1:ifull,wlev+1)=oeu(1:ifull)
eov(1:ifull,wlev+1)=oev(1:ifull)
call boundsuv(eou,eov,stag=-9)
oeu(ifull+1:ifull+iextra)=eou(ifull+1:ifull+iextra,wlev+1)
oev(ifull+1:ifull+iextra)=eov(ifull+1:ifull+iextra,wlev+1)
ccu(:,1)=eou(:,1)*godsig(1)
ccv(:,1)=eov(:,1)*godsig(1)
do ii=2,wlev
  ccu(:,ii)=(ccu(:,ii-1)+eou(:,ii)*godsig(ii))
  ccv(:,ii)=(ccv(:,ii-1)+eov(:,ii)*godsig(ii))
end do
sdiv=(ccu(1:ifull,wlev)*max(ddu(1:ifull)+oeu(1:ifull),0.)/emu(1:ifull)-ccu(iwu,wlev)*max(ddu(iwu)+oeu(iwu),0.)/emu(iwu)  &
     +ccv(1:ifull,wlev)*max(ddv(1:ifull)+oev(1:ifull),0.)/emv(1:ifull)-ccv(isv,wlev)*max(ddv(isv)+oev(isv),0.)/emv(isv)) &
     *em(1:ifull)*em(1:ifull)/ds
do ii=1,wlev-1
  div=(ccu(1:ifull,ii)*max(ddu(1:ifull)+oeu(1:ifull),0.)/emu(1:ifull)-ccu(iwu,ii)*max(ddu(iwu)+oeu(iwu),0.)/emu(iwu)  &
      +ccv(1:ifull,ii)*max(ddv(1:ifull)+oev(1:ifull),0.)/emv(1:ifull)-ccv(isv,ii)*max(ddv(isv)+oev(isv),0.)/emv(isv)) &
      *em(1:ifull)*em(1:ifull)/ds
  nw(:,ii)=(sdiv*gosigh(ii)-div)*ee(1:ifull)
end do
! compute contunity equation horizontal transport terms
do ii=1,wlev
  div=(eou(1:ifull,ii)*(ddu(1:ifull)+neta(1:ifull))/emu(1:ifull)-eou(iwu,ii)*(ddu(iwu)+neta(1:ifull))/emu(iwu)  &
      +eov(1:ifull,ii)*(ddv(1:ifull)+neta(1:ifull))/emv(1:ifull)-eov(isv,ii)*(ddv(isv)+neta(1:ifull))/emv(isv)) &
      *em(1:ifull)*em(1:ifull)/ds
  div=div+(nw(:,ii)-nw(:,ii-1))/godsig(ii) ! nw is at half levels
  mps(1:ifull,ii)=neta(1:ifull)-(1.-ocneps)*0.5*dt*div
  mps(1:ifull,ii)=mps(1:ifull,ii)*ee(1:ifull)
end do

if (myid==0.and.nmaxpr==1) then
  write(6,*) "mlohadv: pressure gradient terms"
end if

! ocean
! Prepare pressure gradient terms at t=t and incorporate into velocity field
do ii=1,wlev
  rhou=0.5*(rho(1:ifull,ii)+rho(ie,ii))
  rhov=0.5*(rho(1:ifull,ii)+rho(in,ii))
  rhobaru=0.5*(rhobar(1:ifull,ii)+rhobar(ie,ii))
  rhobarv=0.5*(rhobar(1:ifull,ii)+rhobar(in,ii))
  tau(:,ii)=grav*(gosig(ii)*max(oeu(1:ifull)+ddu(1:ifull),0.)*drhobardxu(:,ii) &
           +rhobaru*(neta(ie)-neta(1:ifull))*emu(1:ifull)/ds)/rhou ! staggered
  tav(:,ii)=grav*(gosig(ii)*max(oev(1:ifull)+ddv(1:ifull),0.)*drhobardyv(:,ii) &
           +rhobarv*(neta(in)-neta(1:ifull))*emv(1:ifull)/ds)/rhov
end do
! ice
tau(:,wlev+1)=grav*(neta(ie)-neta(1:ifull))*emu(1:ifull)/ds ! staggered
tav(:,wlev+1)=grav*(neta(in)-neta(1:ifull))*emv(1:ifull)/ds
call mlounstaguv(tau,tav,ttau,ttav,toff=1) ! trick from JLM
! ocean
do ii=1,wlev
  uau(:,ii)=nu(1:ifull,ii)+(1.-ocneps)*0.5*dt*( f(1:ifull)*nv(1:ifull,ii)-ttau(:,ii)) ! unstaggered
  uav(:,ii)=nv(1:ifull,ii)+(1.-ocneps)*0.5*dt*(-f(1:ifull)*nu(1:ifull,ii)-ttav(:,ii))
  uau(:,ii)=uau(:,ii)*ee(1:ifull)
  uav(:,ii)=uav(:,ii)*ee(1:ifull)
end do
! ice
snu(1:ifull)=i_u-dt*ttau(:,wlev+1)
snv(1:ifull)=i_v-dt*ttav(:,wlev+1)

call end_log(waterhadv_end)
call start_log(watervadv_begin)

if (myid==0.and.nmaxpr==1) then
  write(6,*) "mlohadv: water vertical advection 1"
end if

! Vertical advection (first call for 0.5*dt)
dums=ns(1:ifull,:)
dumt=nt(1:ifull,:)
duma=mps(1:ifull,:)
call mlovadv(0.5*dt,nw,uau,uav,dums,dumt,duma,depdum,dzdum,wtr(1:ifull),1)
ns(1:ifull,:)=dums
nt(1:ifull,:)=dumt
mps(1:ifull,:)=duma

call end_log(watervadv_end)
call start_log(waterhadv_begin)

if (myid==0.and.nmaxpr==1) then
  write(6,*) "mlohadv: horizontal advection"
end if

! Convert (u,v) to cartesian coordinates (U,V,W)
do ii=1,wlev
  cou(1:ifull,ii,1)=ax(1:ifull)*uau(:,ii)+bx(1:ifull)*uav(:,ii)
  cou(1:ifull,ii,2)=ay(1:ifull)*uau(:,ii)+by(1:ifull)*uav(:,ii)
  cou(1:ifull,ii,3)=az(1:ifull)*uau(:,ii)+bz(1:ifull)*uav(:,ii)
end do

! Horizontal advection for U,V,W
call mlob2ints(cou(:,:,1:3),nface,xg,yg,wtr)

! Rotate vector to arrival point
call mlorot(cou(:,:,1),cou(:,:,2),cou(:,:,3),x3d,y3d,z3d)

! Convert (U,V,W) back to conformal cubic coordinates
do ii=1,wlev
  uau(:,ii)=ax(1:ifull)*cou(1:ifull,ii,1)+ay(1:ifull)*cou(1:ifull,ii,2)+az(1:ifull)*cou(1:ifull,ii,3)
  uav(:,ii)=bx(1:ifull)*cou(1:ifull,ii,1)+by(1:ifull)*cou(1:ifull,ii,2)+bz(1:ifull)*cou(1:ifull,ii,3)
  uau(:,ii)=uau(:,ii)*ee(1:ifull)
  uav(:,ii)=uav(:,ii)*ee(1:ifull)
end do

! Horizontal advection for continuity terms, T and S
cou(1:ifull,:,1)=mps(1:ifull,:)
cou(1:ifull,:,2)=nt(1:ifull,:)-290.
cou(1:ifull,:,3)=ns(1:ifull,:)-34.72
call mlob2intsb(cou(:,:,1:3),nface,xg,yg,wtr)
mps(1:ifull,:)=cou(1:ifull,:,1)
nt(1:ifull,:) =cou(1:ifull,:,2)+290.
ns(1:ifull,:) =cou(1:ifull,:,3)+34.72

call end_log(waterhadv_end)
call start_log(watervadv_begin)

if (myid==0.and.nmaxpr==1) then
  write(6,*) "mlohadv: water vertical advection 2"
end if

! Vertical advection (second call for 0.5*dt)
! use explicit nw and depdum,dzdum from t=tau step (i.e., following JLM in CCAM atmospheric dynamics)
! Could use nuh and nvh to estimate nw at t+1/2, but also require an estimate of neta at t+1/2
dums=ns(1:ifull,:)
dumt=nt(1:ifull,:)
duma=mps(1:ifull,:)
call mlovadv(0.5*dt,nw,uau,uav,dums,dumt,duma,depdum,dzdum,wtr(1:ifull),2)
ns(1:ifull,:)=dums
nt(1:ifull,:)=dumt
mps(1:ifull,:)=duma

! Integrate advected mps
xps=mps(1:ifull,1)*godsig(1)
do ii=2,wlev
  xps=xps+mps(1:ifull,ii)*godsig(ii)
end do
xps=xps*ee(1:ifull)

call end_log(watervadv_end)
call start_log(watereos_begin)

if (myid==0.and.nmaxpr==1) then
  write(6,*) "mlohadv: density EOS 2"
end if

! Approximate normalised density rhobar at t+1 (unstaggered, using T and S at t+1)
if (nxtrrho==1) then
  cou(1:ifull,:,1)=nt(1:ifull,:)
  cou(1:ifull,:,2)=ns(1:ifull,:)
  call bounds(cou(:,:,1:2),corner=.true.)
  nt(ifull+1:ifull+iextra,:)=cou(ifull+1:ifull+iextra,:,1)
  ns(ifull+1:ifull+iextra,:)=cou(ifull+1:ifull+iextra,:,2)
  call mloexpdensity(rho,dalpha,dbeta,nt,ns,xdzdum,pice,0)
  rhobar(:,1)=(rho(:,1)-1030.)*godsig(1)
  do ii=2,wlev
    rhobar(:,ii)=rhobar(:,ii-1)+(rho(:,ii)-1030.)*godsig(ii)
  end do
  do ii=1,wlev
    rhobar(:,ii)=rhobar(:,ii)/gosigh(ii)+1030.
  end do

  ! update normalised density gradients
  ! method 2
  call tsjacobi(nt,ns,dalpha,dbeta,neta,drhobardxu,drhobardyu,drhobardxv,drhobardyv)
  drhobardxu(:,1)=drhobardxu(:,1)*godsig(1)
  drhobardxv(:,1)=drhobardxv(:,1)*godsig(1)
  drhobardyu(:,1)=drhobardyu(:,1)*godsig(1)
  drhobardyv(:,1)=drhobardyv(:,1)*godsig(1)
  do ii=2,wlev
    drhobardxu(:,ii)=drhobardxu(:,ii-1)+drhobardxu(:,ii)*godsig(ii)
    drhobardxv(:,ii)=drhobardxv(:,ii-1)+drhobardxv(:,ii)*godsig(ii)
    drhobardyu(:,ii)=drhobardyu(:,ii-1)+drhobardyu(:,ii)*godsig(ii)
    drhobardyv(:,ii)=drhobardyv(:,ii-1)+drhobardyv(:,ii)*godsig(ii)
  end do
  do ii=1,wlev
    drhobardxu(:,ii)=drhobardxu(:,ii)/gosigh(ii)
    drhobardxv(:,ii)=drhobardxv(:,ii)/gosigh(ii)
    drhobardyu(:,ii)=drhobardyu(:,ii)/gosigh(ii)
    drhobardyv(:,ii)=drhobardyv(:,ii)/gosigh(ii)
  end do
end if

! FREE SURFACE CALCULATION ----------------------------------------

call end_log(watereos_end)
call start_log(waterhelm_begin)

if (myid==0.and.nmaxpr==1) then
  write(6,*) "mlohadv: free surface"
end if

! Prepare integral terms
sou=0.
spu=0.
squ=0.
ssu=0.
sov=0.
spv=0.
sqv=0.
ssv=0.

! ocean
! Precompute U,V current and integral terms at t+1
do ii=1,wlev
  au=uau(:,ii)
  av=uav(:,ii)
  tau(:,ii)=uau(:,ii)+0.5*dt*(1.+ocneps)*f(1:ifull)*av
  tav(:,ii)=uav(:,ii)-0.5*dt*(1.+ocneps)*f(1:ifull)*au
end do
! ice
tau(:,wlev+1)=snu(1:ifull)+dt*f(1:ifull)*snv(1:ifull) ! unstaggered
tav(:,wlev+1)=snv(1:ifull)-dt*f(1:ifull)*snu(1:ifull)
call mlostaguv(tau,tav,ttau,ttav)
! ocean
odum=1./(1.+(1.+ocneps)*(1.+ocneps)*0.25*dt*dt*fu(1:ifull)*fu(1:ifull))
odum=odum*eeu(1:ifull)
dumf=-(1.+ocneps)*0.5*dt*odum
do ii=1,wlev
  ccu(1:ifull,ii)=ttau(:,ii)*odum ! staggered
end do
odum=1./(1.+(1.+ocneps)*(1.+ocneps)*0.25*dt*dt*fv(1:ifull)*fv(1:ifull))
odum=odum*eev(1:ifull)
dumg=-(1.+ocneps)*0.5*dt*odum
do ii=1,wlev
  ccv(1:ifull,ii)=ttav(:,ii)*odum ! staggered
end do
! ice
! niu and niv hold the free drift solution (staggered).  Wind stress terms are updated in mlo.f90
niu(1:ifull)=ttau(:,wlev+1)/(1.+dt*dt*fu(1:ifull)*fu(1:ifull)) ! staggered
niv(1:ifull)=ttav(:,wlev+1)/(1.+dt*dt*fv(1:ifull)*fv(1:ifull))

do ii=1,wlev
  rhou=0.5*(rho(1:ifull,ii)+rho(ie,ii))
  rhov=0.5*(rho(1:ifull,ii)+rho(in,ii))
  rhobaru=0.5*(rhobar(1:ifull,ii)+rhobar(ie,ii))
  rhobarv=0.5*(rhobar(1:ifull,ii)+rhobar(in,ii))

  ! Create arrays to calcuate u and v at t+1, based on pressure gradient at t+1

  ! u^(t+1) = nu = au^(t*) + bu*dpdx^(t+1) + cu*dpdy^(t+1) (staggered)
  ! v^(t+1) = nv = av^(t*) + bv*dpdy^(t+1) + cv*dpdx^(t+1) (staggered)

  au=ccu(1:ifull,ii)
  bu=dumf/rhou
  cu= (1.+ocneps)*0.5*dt*bu ! fu now included in dpdy

  av=ccv(1:ifull,ii)
  bv=dumg/rhov
  cv=-(1.+ocneps)*0.5*dt*bv ! fv now included in dpdx

  ! Note pressure gradients are along constant newz surfaces
  !dppdxu=dpsdxu+grav*sig*(etau+ddu)*drhobardxu+grav*rhobaru*detadxu
  !dppdyu=dpsdyu+grav*sig*(etau+ddu)*drhobardyu+grav*rhobaru*detadyu
  !dppdxv=dpsdxv+grav*sig*(etav+ddv)*drhobardxv+grav*rhobarv*detadxv
  !dppdyv=dpsdyv+grav*sig*(etav+ddv)*drhobardyv+grav*rhobarv*detadyv

  ! Create arrays for u and v at t+1 in terms of neta gradients
    
  !nu=kku+ppu+oou*etau+llu*(etau+ddu)+mmu*detadxu+nnu*detadyu (staggered)
  !nv=kkv+ppv+oov*etav+llv*(etav+ddv)+mmv*detadyv+nnv*detadxv (staggered)

  llu(:,ii)=(bu*drhobardxu(:,ii)+cu*drhobardyu(:,ii))*grav*gosig(ii) ! drhobardyu contains dfdyu term
  mmu(:,ii)=bu*grav*rhobaru
  nnu(:,ii)=cu*grav*rhobaru
  oou(:,ii)=-nnu(:,ii)*dfdyu
  kku(:,ii)=au+bu*(dpsdxu+grav*rhou*dttdxu)*2./(1.+ocneps) &
              -cu*dfdyu*(piceu+grav*rhou*tideu)*2./(1.+ocneps)
  ppu(:,ii)=cu*(dpsdyu+grav*rhou*dttdyu)*2./(1.+ocneps)

  llv(:,ii)=(bv*drhobardyv(:,ii)+cv*drhobardxv(:,ii))*grav*gosig(ii) ! drhobardxv contains dfdxv term
  mmv(:,ii)=bv*grav*rhobarv
  nnv(:,ii)=cv*grav*rhobarv
  oov(:,ii)=-nnv(:,ii)*dfdxv
  kkv(:,ii)=av+bv*(dpsdyv+grav*rhov*dttdyv)*2./(1.+ocneps) &
              -cv*dfdxv*(picev+grav*rhov*tidev)*2./(1.+ocneps)
  ppv(:,ii)=cv*(dpsdxv+grav*rhov*dttdxv)*2./(1.+ocneps)

  llu(:,ii)=llu(:,ii)*eeu(1:ifull)
  mmu(:,ii)=mmu(:,ii)*eeu(1:ifull)
  nnu(:,ii)=nnu(:,ii)*eeu(1:ifull)
  oou(:,ii)=oou(:,ii)*eeu(1:ifull)
  kku(:,ii)=kku(:,ii)*eeu(1:ifull)
  ppu(:,ii)=ppu(:,ii)*eeu(1:ifull)

  llv(:,ii)=llv(:,ii)*eev(1:ifull)
  mmv(:,ii)=mmv(:,ii)*eev(1:ifull)
  nnv(:,ii)=nnv(:,ii)*eev(1:ifull)
  oov(:,ii)=oov(:,ii)*eev(1:ifull)
  kkv(:,ii)=kkv(:,ii)*eev(1:ifull)
  ppv(:,ii)=ppv(:,ii)*eev(1:ifull)

  ! Pre-integrate arrays for u and v at t+1 (i.e., for calculating net divergence at t+1)

  !int nu dz = sou+ssu*etau+spu*(etau+ddu)+squ*detadxu
  !int nv dz = sov+ssv*etav+spv*(etav+ddv)+sqv*detadyv

  sou(1:ifull)=sou(1:ifull)+kku(:,ii)*godsig(ii)
  spu(1:ifull)=spu(1:ifull)+llu(:,ii)*godsig(ii)
  squ(1:ifull)=squ(1:ifull)+mmu(:,ii)*godsig(ii)
  ssu(1:ifull)=ssu(1:ifull)+oou(:,ii)*godsig(ii)

  sov(1:ifull)=sov(1:ifull)+kkv(:,ii)*godsig(ii)
  spv(1:ifull)=spv(1:ifull)+llv(:,ii)*godsig(ii)
  sqv(1:ifull)=sqv(1:ifull)+mmv(:,ii)*godsig(ii)
  ssv(1:ifull)=ssv(1:ifull)+oov(:,ii)*godsig(ii)

end do

! calculate terms for ice velocity at t+1

! (staggered)
! niu(t+1) = niu + idu*ipu + ibu*dipdx + icu*dipdy
! niv(t+1) = niv + idv*ipv + ibv*dipdy + icv*dipdx

imu=0.5*(imass(1:ifull)+imass(ie))
imv=0.5*(imass(1:ifull)+imass(in))
! note missing 0.5 as ip is the average of t and t+1
ibu(1:ifull)=-dt/(imu*(1.+dt*dt*fu(1:ifull)*fu(1:ifull)))
ibv(1:ifull)=-dt/(imv*(1.+dt*dt*fv(1:ifull)*fv(1:ifull)))
ibu(1:ifull)=ibu(1:ifull)*eeu(1:ifull)
ibv(1:ifull)=ibv(1:ifull)*eev(1:ifull)
icu(1:ifull)= dt*ibu(1:ifull) ! fu now included in dipdy
icv(1:ifull)=-dt*ibv(1:ifull) ! fv now includex in dipdx
idu(1:ifull)=-icu(1:ifull)*dfdyu
idv(1:ifull)=-icv(1:ifull)*dfdxv

dumc(1:ifull,1)=sou(1:ifull)
dumd(1:ifull,1)=sov(1:ifull)
dumc(1:ifull,2)=spu(1:ifull)
dumd(1:ifull,2)=spv(1:ifull)
dumc(1:ifull,3)=squ(1:ifull)
dumd(1:ifull,3)=sqv(1:ifull)
dumc(1:ifull,4)=ssu(1:ifull)
dumd(1:ifull,4)=ssv(1:ifull)
dumc(1:ifull,5)=ibu(1:ifull)
dumd(1:ifull,5)=ibv(1:ifull)
dumc(1:ifull,6)=idu(1:ifull)
dumd(1:ifull,6)=idv(1:ifull)
dumc(1:ifull,7)=niu(1:ifull)
dumd(1:ifull,7)=niv(1:ifull)
call boundsuv(dumc(:,1:7),dumd(:,1:7),stag=-9)
sou(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,1)
sov(ifull+1:ifull+iextra)=dumd(ifull+1:ifull+iextra,1)
spu(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,2)
spv(ifull+1:ifull+iextra)=dumd(ifull+1:ifull+iextra,2)
squ(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,3)
sqv(ifull+1:ifull+iextra)=dumd(ifull+1:ifull+iextra,3)
ssu(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,4)
ssv(ifull+1:ifull+iextra)=dumd(ifull+1:ifull+iextra,4)
ibu(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,5)
ibv(ifull+1:ifull+iextra)=dumd(ifull+1:ifull+iextra,5)
idu(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,6)
idv(ifull+1:ifull+iextra)=dumd(ifull+1:ifull+iextra,6)
niu(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,7)
niv(ifull+1:ifull+iextra)=dumd(ifull+1:ifull+iextra,7)

! prep ocean gradient terms
odiv=(sou(1:ifull)/emu(1:ifull)-sou(iwu)/emu(iwu)  &
     +sov(1:ifull)/emv(1:ifull)-sov(isv)/emv(isv)) &
     *em(1:ifull)*em(1:ifull)/ds
odivb=(sou(1:ifull)*ddu(1:ifull)/emu(1:ifull)-sou(iwu)*ddu(iwu)/emu(iwu)  &
      +sov(1:ifull)*ddv(1:ifull)/emv(1:ifull)-sov(isv)*ddv(isv)/emv(isv)) &
      *em(1:ifull)*em(1:ifull)/ds

pue=spu(1:ifull)*0.5/emu(1:ifull)
puw=spu(iwu)*0.5/emu(iwu)
pvn=spv(1:ifull)*0.5/emv(1:ifull)
pvs=spv(isv)*0.5/emv(isv)
pdiv=(pue-puw+pvn-pvs)*em(1:ifull)*em(1:ifull)/ds
pdivb=(pue*ddu(1:ifull)-puw*ddu(iwu)+pvn*ddv(1:ifull)-pvs*ddv(isv))*em(1:ifull)*em(1:ifull)/ds

sue=-squ(1:ifull)/ds
suw=squ(iwu)/ds
svn=-sqv(1:ifull)/ds
svs=sqv(isv)/ds
sdiv=(sue-suw+svn-svs)*em(1:ifull)*em(1:ifull)/ds
sdivb=(sue*ddu(1:ifull)-suw*ddu(iwu)+svn*ddv(1:ifull)-svs*ddv(isv))*em(1:ifull)*em(1:ifull)/ds

sue=squ(1:ifull)/ds
suw=-squ(iwu)/ds
svn=sqv(1:ifull)/ds
svs=-sqv(isv)/ds

que=ssu(1:ifull)*0.5/emu(1:ifull)
quw=ssu(iwu)*0.5/emu(iwu)
qvn=ssv(1:ifull)*0.5/emv(1:ifull)
qvs=ssv(isv)*0.5/emv(isv)
qdiv=(que-quw+qvn-qvs)*em(1:ifull)*em(1:ifull)/ds
qdivb=(que*ddu(1:ifull)-quw*ddu(iwu)+qvn*ddv(1:ifull)-qvs*ddv(isv))*em(1:ifull)*em(1:ifull)/ds

! Iteratively solve for free surface height, eta
! Iterative loop to estimate ice 'pressure'
if (precon<-9999) then
  ! Multi-grid
  call mlomg(tol,itol,neta,sue,svn,suw,svs,pue,pvn,puw,pvs,que,qvn,quw,qvs,                      &
             pdiv,pdivb,sdiv,sdivb,odiv,odivb,qdiv,qdivb,xps,                                    &
             ipice,ibu,ibv,idu,idv,niu,niv,sicedep,ipmax,totits,itc,maxglobseta,maxglobip)
else
  ! Usual SOR
  call mlosor(tol,itol,neta,sue,svn,suw,svs,pue,pvn,puw,pvs,que,qvn,quw,qvs,                     &
              pdiv,pdivb,sdiv,sdivb,odiv,odivb,qdiv,qdivb,xps,                                   &
              ipice,ibu,ibv,idu,idv,niu,niv,sicedep,ipmax,totits,itc,maxglobseta,maxglobip)
end if

if (myid==0.and.nmaxpr==1) then
  write(6,*) "mlohadv: free surface conservation"
end if

! volume conservation for water ---------------------------------------
if (nud_sfh==0) then
  if (fixheight==0) then
    odum=(neta(1:ifull)-w_e)*ee(1:ifull)
    call ccglobal_posneg(odum,delpos,delneg)
    alph_p = -delneg/max(delpos,1.E-20)
    alph_p = min(max(sqrt(alph_p),1.E-20),1.E20)
    neta(1:ifull)=w_e+max(0.,odum)*alph_p+min(0.,odum)/alph_p
  else
    odum=neta(1:ifull)*ee(1:ifull)
    call ccglobal_posneg(odum,delpos,delneg)
    alph_p = -delneg/max(delpos,1.E-20)
    alph_p = min(max(sqrt(alph_p),1.E-20),1.E20)
    neta(1:ifull)=max(0.,odum)*alph_p+min(0.,odum)/alph_p
  end if
end if

if (myid==0.and.nmaxpr==1) then
  write(6,*) "mlohadv: update currents"
end if

dumc(1:ifull,1)=neta(1:ifull)
dumc(1:ifull,2)=ipice(1:ifull)
call bounds(dumc(:,1:2),corner=.true.)
neta(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,1)
ipice(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,2)

tnu=0.5*(neta(in)*f(in)+neta(ine)*f(ine))
tee=0.5*(neta(1:ifull)*f(1:ifull)+neta(ie)*f(ie))
tsu=0.5*(neta(is)*f(is)+neta(ise)*f(ise))
tev=0.5*(neta(ie)*f(ie)+neta(ien)*f(ien))
tnn=0.5*(neta(1:ifull)*f(1:ifull)+neta(in)*f(in))
twv=0.5*(neta(iw)*f(iw)+neta(iwn)*f(iwn))
detadxu=(neta(ie)-neta(1:ifull))*emu(1:ifull)/ds
detadyu=0.5*(stwgt(:,1)*(tnu-tee)+stwgt(:,2)*(tee-tsu))*emu(1:ifull)/ds
detadxv=0.5*(stwgt(:,3)*(tev-tnn)+stwgt(:,4)*(tnn-twv))*emv(1:ifull)/ds
detadyv=(neta(in)-neta(1:ifull))*emv(1:ifull)/ds
oeu(1:ifull)=0.5*(neta(1:ifull)+neta(ie))
oev(1:ifull)=0.5*(neta(1:ifull)+neta(in))

! Update currents once neta is calculated
do ii=1,wlev
  ! update currents (staggered)
  nu(1:ifull,ii)=kku(:,ii)+ppu(:,ii)+oou(:,ii)*oeu+llu(:,ii)*max(oeu(1:ifull)+ddu(1:ifull),0.)+mmu(:,ii)*detadxu+nnu(:,ii)*detadyu
  nv(1:ifull,ii)=kkv(:,ii)+ppv(:,ii)+oov(:,ii)*oev+llv(:,ii)*max(oev(1:ifull)+ddv(1:ifull),0.)+mmv(:,ii)*detadyv+nnv(:,ii)*detadxv
end do

call end_log(waterhelm_end)
call start_log(wateriadv_begin)

! Update ice velocity with internal pressure terms
tnu=0.5*(ipice(in)*f(in)+ipice(ine)*f(ine))
tee=0.5*(ipice(1:ifull)*f(1:ifull)+ipice(ie)*f(ie))
tsu=0.5*(ipice(is)*f(is)+ipice(ise)*f(ise))
tev=0.5*(ipice(ie)*f(ie)+ipice(ien)*f(ien))
tnn=0.5*(ipice(1:ifull)*f(1:ifull)+ipice(in)*f(in))
twv=0.5*(ipice(iw)*f(iw)+ipice(iwn)*f(iwn))
dipdxu=(ipice(ie)-ipice(1:ifull))*emu(1:ifull)/ds
dipdyu=0.5*(stwgt(:,1)*(tnu-tee)+stwgt(:,2)*(tee-tsu))*emu(1:ifull)/ds
dipdxv=0.5*(stwgt(:,3)*(tev-tnn)+stwgt(:,4)*(tnn-twv))*emv(1:ifull)/ds
dipdyv=(ipice(in)-ipice(1:ifull))*emv(1:ifull)/ds
ipiceu=0.5*(ipice(1:ifull)+ipice(ie))
ipicev=0.5*(ipice(1:ifull)+ipice(in))
! tiu and tiv average the velocity over the timestep
tiu(1:ifull)=niu(1:ifull)+0.5*(idu(1:ifull)*ipiceu+ibu(1:ifull)*dipdxu+icu(1:ifull)*dipdyu)
tiv(1:ifull)=niv(1:ifull)+0.5*(idv(1:ifull)*ipicev+ibv(1:ifull)*dipdyv+icv(1:ifull)*dipdxv)
tiu(1:ifull)=tiu(1:ifull)*eeu(1:ifull)
tiv(1:ifull)=tiv(1:ifull)*eev(1:ifull)
niu(1:ifull)=niu(1:ifull)+idu(1:ifull)*ipiceu+ibu(1:ifull)*dipdxu+icu(1:ifull)*dipdyu
niv(1:ifull)=niv(1:ifull)+idv(1:ifull)*ipicev+ibv(1:ifull)*dipdyv+icv(1:ifull)*dipdxv
niu(1:ifull)=niu(1:ifull)*eeu(1:ifull)
niv(1:ifull)=niv(1:ifull)*eev(1:ifull)
call boundsuv(tiu,tiv,stag=-9)

! Normalisation factor for conserving ice flow in and out of gridbox
spnet(1:ifull)=(-min(tiu(iwu)*emu(iwu),0.)+max(tiu(1:ifull)*emu(1:ifull),0.) &
                -min(tiv(isv)*emv(isv),0.)+max(tiv(1:ifull)*emv(1:ifull),0.))/ds

! ADVECT ICE ------------------------------------------------------
! use simple upwind scheme

if (myid==0.and.nmaxpr==1) then
  write(6,*) "mlohadv: Ice advection"
end if

! Horizontal advection for ice area
dumc(1:ifull,1)=fracice/(em(1:ifull)*em(1:ifull)) ! dumc is an area
! Horizontal advection for ice volume
dumc(1:ifull,2)=sicedep*fracice/(em(1:ifull)*em(1:ifull)) ! dumc is a volume
! Horizontal advection for snow volume
dumc(1:ifull,3)=snowd*0.001*fracice/(em(1:ifull)*em(1:ifull)) ! dumc is a volume
! Horizontal advection for ice energy store
dumc(1:ifull,4)=i_sto*fracice/(em(1:ifull)*em(1:ifull))
! Horizontal advection for ice salinity
dumc(1:ifull,5)=i_sal*fracice*sicedep/(em(1:ifull)*em(1:ifull))
! Horizontal advection for surface temperature
call mloexpscalar(0,gamm,0)
dumc(1:ifull,6)=i_it(1:ifull,1)*fracice*gamm/(em(1:ifull)*em(1:ifull))
! Horizontal advection of snow temperature
dumc(1:ifull,7)=i_it(1:ifull,2)*fracice*snowd*0.001/(em(1:ifull)*em(1:ifull))
! Horizontal advection of ice temperature
do ii=3,4
  dumc(1:ifull,5+ii)=i_it(1:ifull,ii)*fracice*sicedep/(em(1:ifull)*em(1:ifull))
end do
! Conservation
dumc(1:ifull,10)=spnet(1:ifull)
call bounds(dumc(:,1:10))
spnet(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,10)
do ii=1,9
  call upwindadv(dumc(:,ii),tiu,tiv,spnet)
end do  
nfracice(1:ifull)=dumc(1:ifull,1)*em(1:ifull)*em(1:ifull)
nfracice(1:ifull)=min(max(nfracice(1:ifull),0.),maxicefrac)
ndic(1:ifull)=dumc(1:ifull,2)*em(1:ifull)*em(1:ifull)/max(nfracice(1:ifull),1.E-10)
ndsn(1:ifull)=dumc(1:ifull,3)*em(1:ifull)*em(1:ifull)/max(nfracice(1:ifull),1.E-10)
nsto(1:ifull)=dumc(1:ifull,4)*em(1:ifull)*em(1:ifull)/max(nfracice(1:ifull),1.E-10)
nis(1:ifull)=dumc(1:ifull,5)*em(1:ifull)*em(1:ifull)/max(ndic(1:ifull)*nfracice(1:ifull),1.E-10)
nit(1:ifull,1)=dumc(1:ifull,6)*em(1:ifull)*em(1:ifull)/max(gamm*nfracice(1:ifull),1.E-10)
nit(1:ifull,2)=dumc(1:ifull,7)*em(1:ifull)*em(1:ifull)/max(ndsn(1:ifull)*nfracice(1:ifull),1.E-10)
do ii=3,4
  nit(1:ifull,ii)=dumc(1:ifull,5+ii)*em(1:ifull)*em(1:ifull)/max(ndic(1:ifull)*nfracice(1:ifull),1.E-10)
end do

! populate grid points that have no sea ice
where (nfracice(1:ifull)<1.E-4)
  ndic(1:ifull)=sicedep
  ndsn(1:ifull)=snowd*0.001
  nsto(1:ifull)=i_sto
  nis(1:ifull)=i_sal
  nit(1:ifull,1)=i_it(:,1)
  nit(1:ifull,2)=i_it(:,2)
  nit(1:ifull,3)=i_it(:,3)
  nit(1:ifull,4)=i_it(:,4)
elsewhere (ndsn(1:ifull)<1.E-4)
  nit(1:ifull,2)=i_it(:,2)
end where

! ocean
! use nu and nv for unstagged currents
tau(:,1:wlev)=nu(1:ifull,:)
tav(:,1:wlev)=nv(1:ifull,:)
! ice
! unstagger ice velocities
tau(:,wlev+1)=niu(1:ifull)
tav(:,wlev+1)=niv(1:ifull)
call mlounstaguv(tau,tav,ttau,ttav)
! ocean
nu(1:ifull,:)=ttau(:,1:wlev)
nv(1:ifull,:)=ttav(:,1:wlev)
! ice
niu(1:ifull)=ttau(:,wlev+1)
niv(1:ifull)=ttav(:,wlev+1)

call end_log(wateriadv_end)
call start_log(watermisc_begin)

if (myid==0.and.nmaxpr==1) then
  write(6,*) "mlohadv: conserve salinity"
end if
  
  
! salinity conservation
if (nud_sss==0) then
  delpos=0.
  delneg=0.
  ndum=0.
  do ii=1,wlev
    ndum=ndum+w_s(1:ifull,ii)*godsig(ii)
  end do
  if (fixsal==0) then
    do ii=1,wlev
      where(wtr(1:ifull).and.ndum>0.)
        dum(:,ii)=ns(1:ifull,ii)-w_s(:,ii)
      elsewhere
        dum(:,ii)=0.
      end where
    end do
    call ccglobal_posneg(dum,delpos,delneg,dsigin=godsig)
    alph_p = -delneg/max(delpos,1.E-20)
    alph_p = min(sqrt(alph_p),alph_p)
    do ii=1,wlev
      where(wtr(1:ifull).and.ndum>0.)
        ns(1:ifull,ii)=w_s(:,ii)+max(0.,dum(:,ii))*alph_p+min(0.,dum(:,ii))/max(1.,alph_p)
      end where
    end do
  else
    do ii=1,wlev
      where(wtr(1:ifull).and.ndum>0.)
        dum(:,ii)=ns(1:ifull,ii)-34.72
      elsewhere
        dum(:,ii)=0.
      end where
    end do
    call ccglobal_posneg(dum,delpos,delneg,dsigin=godsig)
    alph_p = -delneg/max(delpos,1.E-20)
    alph_p = min(sqrt(alph_p),alph_p)
    do ii=1,wlev
      where(wtr(1:ifull).and.ndum>0.)
        ns(1:ifull,ii)=34.72+max(0.,dum(:,ii))*alph_p+min(0.,dum(:,ii))/max(1.,alph_p)
      end where
    end do
  end if
end if

if (myid==0.and.(ktau<=5.or.maxglobseta>tol.or.maxglobip>itol)) then
  write(6,*) "MLODYNAMICS ",totits,itc,maxglobseta,maxglobip
end if

call end_log(watermisc_end)

! DIFFUSION -------------------------------------------------------------------
if (myid==0.and.nmaxpr==1) then
  write(6,*) "mlohadv: diffusion"
end if

duma=nu(1:ifull,:)
dumb=nv(1:ifull,:)
uau=av_vmod*nu(1:ifull,:)+(1.-av_vmod)*oldu1
uav=av_vmod*nv(1:ifull,:)+(1.-av_vmod)*oldv1
dumt=nt(1:ifull,:)
dums=ns(1:ifull,:)
odum=neta(1:ifull)
call mlodiffusion_main(uau,uav,duma,dumb,odum,dumt,dums)
call mloimport(4,odum,0,0)


! EXPORT ----------------------------------------------------------------------
! Water data is exported in mlodiffusion

call start_log(watermisc_begin)

if (myid==0.and.nmaxpr==1) then
  write(6,*) "mlohadv: Export"
end if

do ii=1,4
  call mloimpice(nit(1:ifull,ii),ii,0)
end do
call mloimpice(nfracice(1:ifull),5,0)
call mloimpice(ndic(1:ifull),6,0)
call mloimpice(ndsn(1:ifull),7,0)
call mloimpice(nsto(1:ifull),8,0)
call mloimpice(niu(1:ifull),9,0)
call mloimpice(niv(1:ifull),10,0)
call mloimpice(nis(1:ifull),11,0)
where (wtr(1:ifull))
  fracice=nfracice(1:ifull)
  sicedep=ndic(1:ifull)
  snowd=ndsn(1:ifull)*1000.
end where

if (myid==0.and.nmaxpr==1) then
  write(6,*) "mlohadv: Finish"
end if

call end_log(watermisc_end)

return
end subroutine mlohadv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate depature points for MLO semi-Lagrangian advection
! (This subroutine is based on depts.f)

subroutine mlodeps(dtin,ubar,vbar,nface,xg,yg,x3d,y3d,z3d,wtr)

use cc_mpi
use mlo
use vecsuv_m
use xyzinfo_m

implicit none

include 'newmpar.h'
include 'const_phys.h'

integer ii,n,kx
integer, dimension(:,:), intent(out) :: nface
real, intent(in) :: dtin
real, dimension(ifull,size(nface,2)), intent(in) :: ubar,vbar
real, dimension(ifull,size(nface,2)), intent(out) :: xg,yg
real*8, dimension(ifull,size(nface,2)), intent(out) :: x3d,y3d,z3d
real, dimension(ifull,size(nface,2)) :: uc,vc,wc
real, dimension(ifull+iextra,size(nface,2),3) :: temp
logical, dimension(ifull+iextra), intent(in) :: wtr
integer, parameter :: nguess = 2

kx=size(nface,2)

! departure point x, y, z is called x3d, y3d, z3d
! first find corresponding cartesian vels
do ii=1,kx
  uc(:,ii)=(ax(1:ifull)*ubar(:,ii)+bx(1:ifull)*vbar(:,ii))*dtin/rearth ! unit sphere 
  vc(:,ii)=(ay(1:ifull)*ubar(:,ii)+by(1:ifull)*vbar(:,ii))*dtin/rearth ! unit sphere 
  wc(:,ii)=(az(1:ifull)*ubar(:,ii)+bz(1:ifull)*vbar(:,ii))*dtin/rearth ! unit sphere 
  x3d(:,ii)=x-uc(:,ii) ! 1st guess
  y3d(:,ii)=y-vc(:,ii)
  z3d(:,ii)=z-wc(:,ii)
end do

! convert to grid point numbering
call mlotoij5(x3d,y3d,z3d,nface,xg,yg)
! Share off processor departure points.
call deptsync(nface,xg,yg)

do n=1,nguess
  temp(1:ifull,:,1) = uc
  temp(1:ifull,:,2) = vc
  temp(1:ifull,:,3) = wc
  call mlob2ints(temp(:,:,1:3),nface,xg,yg,wtr)
  do ii=1,kx
    x3d(:,ii) = x - 0.5*(uc(:,ii)+temp(1:ifull,ii,1)) ! n+1 guess
    y3d(:,ii) = y - 0.5*(vc(:,ii)+temp(1:ifull,ii,2)) ! n+1 guess
    z3d(:,ii) = z - 0.5*(wc(:,ii)+temp(1:ifull,ii,3)) ! n+1 guess
  end do
  call mlotoij5(x3d,y3d,z3d,nface,xg,yg)
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
integer ii,kx
integer, dimension(:,:), intent(out) :: nface
real, dimension(ifull,size(nface,2)), intent(out) :: xg,yg
real, dimension(ifull) :: xstr,ystr,zstr
real, dimension(ifull) :: denxyz,xd,yd,zd
real, dimension(ifull) :: ri,rj
real dxx,dxy,dyx,dyy
real(kind=8), dimension(ifull,size(nface,2)), intent(inout) :: x3d,y3d,z3d
real(kind=8), dimension(ifull) :: den
real(kind=8) alf,alfonsch
real(kind=8), parameter :: one = 1.
integer, parameter :: nmaploop = 3

kx=size(nface,2)

alf=(one-schmidt**2)/(one+schmidt**2)
alfonsch=2.*schmidt/(one+schmidt**2)  ! same but bit more accurate

do ii=1,kx

  !     if necessary, transform (x3d, y3d, z3d) to equivalent
  !     coordinates (xstr, ystr, zstr) on regular gnomonic panels
  den=one-alf*z3d(:,ii) ! to force real*8
  xstr=x3d(:,ii)*(alfonsch/den)
  ystr=y3d(:,ii)*(alfonsch/den)
  zstr=    (z3d(:,ii)-alf)/den

!      first deduce departure faces
!      instead calculate cubic coordinates
!      The faces are:
!      0: X=1   1: Z=1   2: Y=1   3: X=-1   4: Z=-1   5: Y=-1
  denxyz=max( abs(xstr),abs(ystr),abs(zstr) )
  xd=xstr/denxyz
  yd=ystr/denxyz
  zd=zstr/denxyz

  where (abs(xstr-denxyz)<1.E-8)
    nface(:,ii)    =0
    xg(:,ii) =      yd
    yg(:,ii) =      zd
  elsewhere (abs(xstr+denxyz)<1.E-8)
    nface(:,ii)    =3
    xg(:,ii) =     -zd
    yg(:,ii) =     -yd
  elsewhere (abs(zstr-denxyz)<1.E-8)
    nface(:,ii)    =1
    xg(:,ii) =      yd
    yg(:,ii) =     -xd
  elsewhere (abs(zstr+denxyz)<1.E-8)
    nface(:,ii)    =4
    xg(:,ii) =      xd
    yg(:,ii) =     -yd
  elsewhere (abs(ystr-denxyz)<1.E-8)
    nface(:,ii)    =2
    xg(:,ii) =     -zd
    yg(:,ii) =     -xd
  elsewhere
    nface(:,ii)    =5
    xg(:,ii) =      xd
    yg(:,ii) =      zd
  end where

  !     use 4* resolution grid il --> 4*il
  xg(:,ii)=min(max(-.999999,xg(:,ii)),.999999)
  yg(:,ii)=min(max(-.999999,yg(:,ii)),.999999)
  !      first guess for ri, rj and nearest i,j
  ri=1.+(1.+xg(:,ii))*real(2*il_g)
  rj=1.+(1.+yg(:,ii))*real(2*il_g)
  do loop=1,nmaploop
    do iq=1,ifull
      i=nint(ri(iq))
      j=nint(rj(iq))
      is=nint(sign(1.,ri(iq)-real(i)))
      js=nint(sign(1.,rj(iq)-real(j)))
      ! predict new value for ri, rj
      dxx=xx4(i+is,j)-xx4(i,j)
      dyx=xx4(i,j+js)-xx4(i,j)
      dxy=yy4(i+is,j)-yy4(i,j)
      dyy=yy4(i,j+js)-yy4(i,j)       
      den(iq)=dxx*dyy-dyx*dxy
      ri(iq)=real(i)+real(is)*((xg(iq,ii)-xx4(i,j))*dyy-(yg(iq,ii)-yy4(i,j))*dyx)/den(iq)
      rj(iq)=real(j)+real(js)*((yg(iq,ii)-yy4(i,j))*dxx-(xg(iq,ii)-xx4(i,j))*dxy)/den(iq)
      
      ri(iq) = min(ri(iq),1.0+1.999999*real(2*il_g))
      ri(iq) = max(ri(iq),1.0+0.000001*real(2*il_g))
      rj(iq) = min(rj(iq),1.0+1.999999*real(2*il_g))
      rj(iq) = max(rj(iq),1.0+0.000001*real(2*il_g))
      
    end do
  end do  ! loop loop
  !      expect xg, yg to range between .5 and il+.5
  xg(:,ii)=0.25*(ri+3.)-0.5  ! -.5 for stag; back to normal ri, rj defn
  yg(:,ii)=0.25*(rj+3.)-0.5  ! -.5 for stag

end do

return
end subroutine mlotoij5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate depature points for semi-Lagrangian advection
! This code is from ints.f

! This version is for velocity as the interpolated value is zero over land
! We use bi-cubic based on JLM's ints.f routines
subroutine mlob2ints(s,nface,xg,yg,wtr,bilinear)

use cc_mpi
use indices_m
use mlo

implicit none

include 'newmpar.h'
include 'parm.h'
include 'parmhor.h'

integer idel,iq,jdel,kx
integer i,j,k,n,ind,ip,jp,iproc,ierr,intsch,ncount
integer ii,ntr,nn
integer, dimension(:,:), intent(in) :: nface
real, dimension(ifull,size(nface,2)), intent(in) :: xg,yg
real, dimension(:,:,:), intent(inout) :: s
real, dimension(-1:ipan+2,-1:jpan+2,1:npan,size(s,2),size(s,3)) :: sx
real, dimension(-1:2,-1:2,size(s,3)) :: sc
real, dimension(4) :: r
real xxg,yyg,cxx
real aab,aac,aad
logical, intent(in), optional :: bilinear
logical, dimension(ifull+iextra), intent(in) :: wtr
logical :: lmode
ind(i,j,n)=i+(j-1)*ipan+(n-1)*ipan*jpan  ! *** for n=1,npan

lmode=.true.
if (present(bilinear)) lmode=.not.bilinear

kx=size(s,2)
ntr=size(s,3)
intsch=mod(ktau,2)
cxx=-9999.
sx=cxx-1.
sc=cxx-1.

do iq=1,ifull
  if (.not.wtr(iq)) then
    s(iq,:,:)=cxx-1.
  end if
end do
call bounds(s,nrows=2)

!======================== start of intsch=1 section ====================
if(intsch==1)then

  do n=1,npan         ! first simple copy into larger array
    do j=1,jpan
      do i=1,ipan
        sx(i,j,n,:,:) = s(ind(i,j,n),:,:)
      end do         ! i loop
      sx(0,j,n,:,:)      = s(iw(ind(1,j,n)),:,:)
      sx(-1,j,n,:,:)     = s(iww(ind(1,j,n)),:,:)
      sx(ipan+1,j,n,:,:) = s(ie(ind(ipan,j,n)),:,:)
      sx(ipan+2,j,n,:,:) = s(iee(ind(ipan,j,n)),:,:)
    end do            ! j loop
    do i=1,ipan
      sx(i,0,n,:,:)      = s(is(ind(i,1,n)),:,:)
      sx(i,-1,n,:,:)     = s(iss(ind(i,1,n)),:,:)
      sx(i,jpan+1,n,:,:) = s(in(ind(i,jpan,n)),:,:)
      sx(i,jpan+2,n,:,:) = s(inn(ind(i,jpan,n)),:,:)
    end do            ! i loop
!   for ew interpolation, sometimes need (different from ns):
!       (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
!     (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)

    sx(-1,0,n,:,:)          = s(lwws(n),:,:)
    sx(0,0,n,:,:)           = s(iws(ind(1,1,n)),:,:)
    sx(0,-1,n,:,:)          = s(lwss(n),:,:)
    sx(ipan+1,0,n,:,:)      = s(ies(ind(ipan,1,n)),:,:)
    sx(ipan+2,0,n,:,:)      = s(lees(n),:,:)
    sx(ipan+1,-1,n,:,:)     = s(less(n),:,:)
    sx(-1,jpan+1,n,:,:)     = s(lwwn(n),:,:)
    sx(0,jpan+2,n,:,:)      = s(lwnn(n),:,:)
    sx(ipan+2,jpan+1,n,:,:) = s(leen(n),:,:)
    sx(ipan+1,jpan+2,n,:,:) = s(lenn(n),:,:)
    sx(0,jpan+1,n,:,:)      = s(iwn(ind(1,jpan,n)),:,:)
    sx(ipan+1,jpan+1,n,:,:) = s(ien(ind(ipan,jpan,n)),:,:)
    !sx(-1,0,n,:,:)          = s(iww(ind(1,1,n)),:,:)
    !sx(0,0,n,:,:)           = 0.5*(s(iw(ind(1,1,n)),:,:)+s(is(ind(1,1,n)),:,:))
    !sx(0,-1,n,:,:)          = s(iss(ind(1,1,n)),:,:)
    !sx(ipan+1,0,n,:,:)      = 0.5*(s(ie(ind(ipan,1,n)),:,:)+s(is(ind(ipan,1,n)),:,:))
    !sx(ipan+2,0,n,:,:)      = s(iee(ind(ipan,1,n)),:,:)
    !sx(ipan+1,-1,n,:,:)     = s(iss(ind(ipan,1,n)),:,:)
    !sx(-1,jpan+1,n,:,:)     = s(iww(ind(1,jpan,n)),:,:)
    !sx(0,jpan+2,n,:,:)      = s(inn(ind(1,jpan,n)),:,:)
    !sx(ipan+2,jpan+1,n,:,:) = s(iee(ind(ipan,jpan,n)),:,:)
    !sx(ipan+1,jpan+2,n,:,:) = s(inn(ind(ipan,jpan,n)),:,:)
    !sx(0,jpan+1,n,:,:)      = 0.5*(s(iw(ind(1,jpan,n)),:,:)+s(in(ind(1,jpan,n)),:,:))
    !sx(ipan+1,jpan+1,n,:,:) = 0.5*(s(ie(ind(ipan,jpan,n)),:,:)+s(in(ind(ipan,jpan,n)),:,:))
  end do               ! n loop

! Loop over points that need to be calculated for other processes
  do ii=1,neighnum
    iproc=neighlistsend(ii)
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

      sc(-1,0,:) = sx(idel-1,jdel,n,k,:)
      sc(0,0,:)  = sx(idel  ,jdel,n,k,:)
      sc(1,0,:)  = sx(idel+1,jdel,n,k,:)
      sc(2,0,:)  = sx(idel+2,jdel,n,k,:)

      sc(-1,1,:) = sx(idel-1,jdel+1,n,k,:)
      sc(0,1,:)  = sx(idel  ,jdel+1,n,k,:)
      sc(1,1,:)  = sx(idel+1,jdel+1,n,k,:)
      sc(2,1,:)  = sx(idel+2,jdel+1,n,k,:)

      sc(0,-1,:) = sx(idel  ,jdel-1,n,k,:)
      sc(1,-1,:) = sx(idel+1,jdel-1,n,k,:)

      sc(0,2,:) = sx(idel  ,jdel+2,n,k,:)
      sc(1,2,:) = sx(idel+1,jdel+2,n,k,:)

      ncount=count(sc(:,:,1)>cxx)
      if (ncount>=12.and.lmode) then
        ! bi-cubic interpolation
        do nn=1,ntr
          r(2) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*sc(0,0,nn)-xxg*sc(-1,0,nn)/3.)     &
               -xxg*(1.+xxg)*sc(2,0,nn)/3.)+xxg*(1.+xxg)*(2.-xxg)*sc(1,0,nn))/2.
          r(3) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*sc(0,1,nn)-xxg*sc(-1,1,nn)/3.)     &
               -xxg*(1.+xxg)*sc(2,1,nn)/3.)+xxg*(1.+xxg)*(2.-xxg)*sc(1,1,nn))/2.
          r(1) = (1.-xxg)*sc(0,-1,nn) +xxg*sc(1,-1,nn)
          r(4) = (1.-xxg)*sc(0,2,nn) +xxg*sc(1,2,nn)

          sextra(iproc)%a(nn+(iq-1)*ntr) = ((1.-yyg)*((2.-yyg)* &
               ((1.+yyg)*r(2)-yyg*r(1)/3.)                      & 
               -yyg*(1.+yyg)*r(4)/3.)                           &
               +yyg*(1.+yyg)*(2.-yyg)*r(3))/2.
        end do
      else
        ! bi-linear interpolation
        do nn=1,ntr
          where (sc(0:1,0:1,nn)<=cxx)
            sc(0:1,0:1,nn)=0.
          end where       
          aad=sc(1,1,nn)-sc(0,1,nn)-sc(1,0,nn)+sc(0,0,nn)
          aab=sc(1,0,nn)-sc(0,0,nn)
          aac=sc(0,1,nn)-sc(0,0,nn)
        
          sextra(iproc)%a(nn+(iq-1)*ntr)=aab*xxg+aac*yyg+aad*xxg*yyg+sc(0,0,nn)
        end do
      end if
    end do            ! iq loop
  end do              ! iproc loop
  
  call intssync_send(ntr)

  do iq=1,ifull
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
        if ( idel < 0 .or. idel > ipan .or. jdel < 0 .or. &
             jdel > jpan .or. n < 1 .or. n > npan ) then
          cycle      ! Will be calculated on another processor
        end if

        sc(0,-1,:) = sx(idel  ,jdel-1,n,k,:)
        sc(1,-1,:) = sx(idel+1,jdel-1,n,k,:)

        sc(-1,0,:) = sx(idel-1,jdel,n,k,:)
        sc(0,0,:)  = sx(idel  ,jdel,n,k,:)
        sc(1,0,:)  = sx(idel+1,jdel,n,k,:)
        sc(2,0,:)  = sx(idel+2,jdel,n,k,:)

        sc(-1,1,:) = sx(idel-1,jdel+1,n,k,:)
        sc(0,1,:)  = sx(idel  ,jdel+1,n,k,:)
        sc(1,1,:)  = sx(idel+1,jdel+1,n,k,:)
        sc(2,1,:)  = sx(idel+2,jdel+1,n,k,:)

        sc(0,2,:) = sx(idel  ,jdel+2,n,k,:)
        sc(1,2,:) = sx(idel+1,jdel+2,n,k,:)

        ncount=count(sc(:,:,1)>cxx)
        if (ncount>=12.and.lmode) then
          ! bi-cubic interpolation
          do nn=1,ntr
            r(2) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*sc(0,0,nn)-xxg*sc(-1,0,nn)/3.) &
                 -xxg*(1.+xxg)*sc(2,0,nn)/3.)+xxg*(1.+xxg)*(2.-xxg)*sc(1,0,nn))/2.
            r(3) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*sc(0,1,nn)-xxg*sc(-1,1,nn)/3.) &
                 -xxg*(1.+xxg)*sc(2,1,nn)/3.)+xxg*(1.+xxg)*(2.-xxg)*sc(1,1,nn))/2.
            r(1) = (1.-xxg)*sc(0,-1,nn)+xxg*sc(1,-1,nn)
            r(4) = (1.-xxg)*sc(0,2,nn) +xxg*sc(1,2,nn)

            s(iq,k,nn) = ((1.-yyg)*((2.-yyg)*         &
                 ((1.+yyg)*r(2)-yyg*r(1)/3.)          &
                 -yyg*(1.+yyg)*r(4)/3.)               &
                 +yyg*(1.+yyg)*(2.-yyg)*r(3))/2.
          end do
        else
          ! bi-linear interpolation along coastline
          do nn=1,ntr
            where (sc(0:1,0:1,nn)<=cxx)
              sc(0:1,0:1,nn)=0.
            end where
            aad=sc(1,1,nn)-sc(0,1,nn)-sc(1,0,nn)+sc(0,0,nn)
            aab=sc(1,0,nn)-sc(0,0,nn)
            aac=sc(0,1,nn)-sc(0,0,nn)
            s(iq,k,nn)=aab*xxg+aac*yyg+aad*xxg*yyg+sc(0,0,nn)
          end do
        end if
      end do ! k loop
    end if   ! wtr
  end do      ! iq loop
       
!========================   end of intsch=1 section ====================
else     ! if(intsch==1)then
!======================== start of intsch=2 section ====================
!       this is intsc           NS interps done first
!       first extend s arrays into sx - this one -1:il+2 & -1:il+2

  do n=1,npan         ! first simple copy into larger array
    do j=1,jpan
      do i=1,ipan
        sx(i,j,n,:,:) = s(ind(i,j,n),:,:)
      end do         ! i loop
      sx(0,j,n,:,:)      = s(iw(ind(1,j,n)),:,:)
      sx(-1,j,n,:,:)     = s(iww(ind(1,j,n)),:,:)
      sx(ipan+1,j,n,:,:) = s(ie(ind(ipan,j,n)),:,:)
      sx(ipan+2,j,n,:,:) = s(iee(ind(ipan,j,n)),:,:)
    end do            ! j loop
    do i=1,ipan
      sx(i,0,n,:,:)      = s(is(ind(i,1,n)),:,:)
      sx(i,-1,n,:,:)     = s(iss(ind(i,1,n)),:,:)
      sx(i,jpan+1,n,:,:) = s(in(ind(i,jpan,n)),:,:)
      sx(i,jpan+2,n,:,:) = s(inn(ind(i,jpan,n)),:,:)
    end do            ! i loop
!        for ns interpolation, sometimes need (different from ew):
!            (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
!          (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)

    sx(-1,0,n,:,:)          = s(lsww(n),:,:)
    sx(0,0,n,:,:)           = s(isw(ind(1,1,n)),:,:)
    sx(0,-1,n,:,:)          = s(lssw(n),:,:)
    sx(ipan+2,0,n,:,:)      = s(lsee(n),:,:)
    sx(ipan+1,-1,n,:,:)     = s(lsse(n),:,:)
    sx(-1,jpan+1,n,:,:)     = s(lnww(n),:,:)
    sx(0,jpan+1,n,:,:)      = s(inw(ind(1,jpan,n)),:,:)
    sx(0,jpan+2,n,:,:)      = s(lnnw(n),:,:)
    sx(ipan+2,jpan+1,n,:,:) = s(lnee(n),:,:)
    sx(ipan+1,jpan+2,n,:,:) = s(lnne(n),:,:)
    sx(ipan+1,0,n,:,:)      = s(ise(ind(ipan,1,n)),:,:)
    sx(ipan+1,jpan+1,n,:,:) = s(ine(ind(ipan,jpan,n)),:,:)
  end do               ! n loop

! For other processes
  do ii=1,neighnum
    iproc=neighlistsend(ii)
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

      sc(0,-1,:) = sx(idel,jdel-1,n,k,:)
      sc(0,0,:)  = sx(idel,jdel  ,n,k,:)
      sc(0,1,:)  = sx(idel,jdel+1,n,k,:)
      sc(0,2,:)  = sx(idel,jdel+2,n,k,:)

      sc(1,-1,:) = sx(idel+1,jdel-1,n,k,:)
      sc(1,0,:)  = sx(idel+1,jdel  ,n,k,:)
      sc(1,1,:)  = sx(idel+1,jdel+1,n,k,:)
      sc(1,2,:)  = sx(idel+1,jdel+2,n,k,:)

      sc(-1,0,:) = sx(idel-1,jdel  ,n,k,:)
      sc(-1,1,:) = sx(idel-1,jdel+1,n,k,:)
        
      sc(2,0,:) = sx(idel+2,jdel  ,n,k,:)
      sc(2,1,:) = sx(idel+2,jdel+1,n,k,:)

      ncount=count(sc(:,:,1)>cxx)
      if (ncount>=12.and.lmode) then
        ! bi-cubic interpolation
        do nn=1,ntr
          r(2) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*sc(0,0,nn)-yyg*sc(0,-1,nn)/3.) &
               -yyg*(1.+yyg)*sc(0,2,nn)/3.)+yyg*(1.+yyg)*(2.-yyg)*sc(0,1,nn))/2.
          r(3) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*sc(1,0,nn)-yyg*sc(1,-1,nn)/3.) &
               -yyg*(1.+yyg)*sc(1,2,nn)/3.)+yyg*(1.+yyg)*(2.-yyg)*sc(1,1,nn))/2.
          r(1) = (1.-yyg)*sc(-1,0,nn) +yyg*sc(-1,1,nn)
          r(4) = (1.-yyg)*sc(2,0,nn) +yyg*sc(2,1,nn)
          sextra(iproc)%a(nn+(iq-1)*ntr) = ((1.-xxg)*((2.-xxg)* &
                    ((1.+xxg)*r(2)-xxg*r(1)/3.)                 &
                       -xxg*(1.+xxg)*r(4)/3.)                   &
                       +xxg*(1.+xxg)*(2.-xxg)*r(3))/2.
        end do
      else
        ! bi-linear interpolation
        do nn=1,ntr
          where (sc(0:1,0:1,nn)<=cxx)
            sc(0:1,0:1,nn)=0.
          end where       
          aad=sc(1,1,nn)-sc(0,1,nn)-sc(1,0,nn)+sc(0,0,nn)
          aab=sc(1,0,nn)-sc(0,0,nn)
          aac=sc(0,1,nn)-sc(0,0,nn)
          sextra(iproc)%a(nn+(iq-1)*ntr)=aab*xxg+aac*yyg+aad*xxg*yyg+sc(0,0,nn)
        end do
      end if
    end do            ! iq loop
  end do              ! iproc

  call intssync_send(ntr)

  do iq=1,ifull
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
        if ( idel < 0 .or. idel > ipan .or. jdel < 0 .or. &
             jdel > jpan .or. n < 1 .or. n > npan ) then
          cycle      ! Will be calculated on another processor
        end if

        sc(0,-1,:) = sx(idel,jdel-1,n,k,:)
        sc(0,0,:)  = sx(idel,jdel  ,n,k,:)
        sc(0,1,:)  = sx(idel,jdel+1,n,k,:)
        sc(0,2,:)  = sx(idel,jdel+2,n,k,:)

        sc(1,-1,:) = sx(idel+1,jdel-1,n,k,:)
        sc(1,0,:)  = sx(idel+1,jdel  ,n,k,:)
        sc(1,1,:)  = sx(idel+1,jdel+1,n,k,:)
        sc(1,2,:)  = sx(idel+1,jdel+2,n,k,:)

        sc(-1,0,:) = sx(idel-1,jdel  ,n,k,:)
        sc(-1,1,:) = sx(idel-1,jdel+1,n,k,:)
        
        sc(2,0,:) = sx(idel+2,jdel  ,n,k,:)
        sc(2,1,:) = sx(idel+2,jdel+1,n,k,:)
        
        ncount=count(sc(:,:,1)>cxx)
        if (ncount>=12.and.lmode) then
          ! bi-cubic interpolation
          do nn=1,ntr
            r(2) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*sc(0,0,nn)-yyg*sc(0,-1,nn)/3.) &
                   -yyg*(1.+yyg)*sc(0,2,nn)/3.)+yyg*(1.+yyg)*(2.-yyg)*sc(0,1,nn))/2.
            r(3) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*sc(1,0,nn)-yyg*sc(1,-1,nn)/3.) &
                   -yyg*(1.+yyg)*sc(1,2,nn)/3.)+yyg*(1.+yyg)*(2.-yyg)*sc(1,1,nn))/2.
            r(1) = (1.-yyg)*sc(-1,0,nn)+yyg*sc(-1,1,nn)
            r(4) = (1.-yyg)*sc(2,0,nn) +yyg*sc(2,1,nn)

            s(iq,k,nn) = ((1.-xxg)*((2.-xxg)*    &
                 ((1.+xxg)*r(2)-xxg*r(1)/3.)     &
                 -xxg*(1.+xxg)*r(4)/3.)          &
                 +xxg*(1.+xxg)*(2.-xxg)*r(3))/2.
          end do
        else
          ! bi-linear interpolation
          do nn=1,ntr
            where (sc(0:1,0:1,nn)<=cxx)
              sc(0:1,0:1,nn)=0.
            end where      
            aad=sc(1,1,nn)-sc(0,1,nn)-sc(1,0,nn)+sc(0,0,nn)
            aab=sc(1,0,nn)-sc(0,0,nn)
            aac=sc(0,1,nn)-sc(0,0,nn)
            s(iq,k,nn)=aab*xxg+aac*yyg+aad*xxg*yyg+sc(0,0,nn)
          end do
        end if
      end do
    end if
  end do

endif                     ! (intsch==1) .. else ..
!========================   end of intsch=1 section ====================

call intssync_recv(s)

do iq=1,ifull
  if (.not.wtr(iq)) then
    s(iq,:,:)=0.
  end if
end do

return
end subroutine mlob2ints

! This version is for scalars which is same as above, but
! missing land values are filled from non-trivial values
! instead of being set to zero
subroutine mlob2intsb(s,nface,xg,yg,wtr,bilinear)

use cc_mpi
use indices_m
use mlo

implicit none

include 'newmpar.h'
include 'parm.h'
include 'parmhor.h'

integer idel,iq,jdel,kx
integer i,j,k,n,ind,ip,jp,iproc,ierr,intsch,ncount
integer ii,ntr,nn
integer, dimension(:,:), intent(in) :: nface
real, dimension(ifull,size(nface,2)), intent(in) :: xg,yg
real, dimension(:,:,:), intent(inout) :: s
real, dimension(ifull,size(s,2),size(s,3)) :: ssav
real, dimension(-1:ipan+2,-1:jpan+2,1:npan,size(s,2),size(s,3)) :: sx
real, dimension(-1:2,-1:2,size(s,3)) :: sc
real, dimension(0:1,0:1) :: scb
real, dimension(4) :: r
real xxg,yyg,aab,aac,aad
real cmax,cmin,cxx
logical, intent(in), optional :: bilinear
logical, dimension(ifull+iextra), intent(in) :: wtr
logical :: lmode
ind(i,j,n)=i+(j-1)*ipan+(n-1)*ipan*jpan  ! *** for n=1,npan

lmode=.true.
if (present(bilinear)) lmode=.not.bilinear

kx=size(s,2)
ntr=size(s,3)
intsch=mod(ktau,2)
cxx=-9999. ! missing value flag
sx=cxx-1.
sc=cxx-1.
ssav(1:ifull,:,:)=s(1:ifull,:,:)

do iq=1,ifull
  if (.not.wtr(iq)) then
    s(iq,:,:)=cxx-1. ! missing value flag
  end if
end do
call bounds(s,nrows=2)

!======================== start of intsch=1 section ====================
if(intsch==1)then

  do nn=1,ntr
    do k=1,kx
      do n=1,npan         ! first simple copy into larger array
        do j=1,jpan
          do i=1,ipan
            sx(i,j,n,k,nn) = s(ind(i,j,n),k,nn)
          end do         ! i loop
          sx(0,j,n,k,nn)      = s(iw(ind(1,j,n)),k,nn)
          sx(-1,j,n,k,nn)     = s(iww(ind(1,j,n)),k,nn)
          sx(ipan+1,j,n,k,nn) = s(ie(ind(ipan,j,n)),k,nn)
          sx(ipan+2,j,n,k,nn) = s(iee(ind(ipan,j,n)),k,nn)
        end do            ! j loop
        do i=1,ipan
          sx(i,0,n,k,nn)      = s(is(ind(i,1,n)),k,nn)
          sx(i,-1,n,k,nn)     = s(iss(ind(i,1,n)),k,nn)
          sx(i,jpan+1,n,k,nn) = s(in(ind(i,jpan,n)),k,nn)
          sx(i,jpan+2,n,k,nn) = s(inn(ind(i,jpan,n)),k,nn)
        end do            ! i loop
!   for ew interpolation, sometimes need (different from ns):
!       (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
!     (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)

        sx(-1,0,n,k,nn)          = s(lwws(n),k,nn)
        sx(0,0,n,k,nn)           = s(iws(ind(1,1,n)),k,nn)
        sx(0,-1,n,k,nn)          = s(lwss(n),k,nn)
        sx(ipan+1,0,n,k,nn)      = s(ies(ind(ipan,1,n)),k,nn)
        sx(ipan+2,0,n,k,nn)      = s(lees(n),k,nn)
        sx(ipan+1,-1,n,k,nn)     = s(less(n),k,nn)
        sx(-1,jpan+1,n,k,nn)     = s(lwwn(n),k,nn)
        sx(0,jpan+2,n,k,nn)      = s(lwnn(n),k,nn)
        sx(ipan+2,jpan+1,n,k,nn) = s(leen(n),k,nn)
        sx(ipan+1,jpan+2,n,k,nn) = s(lenn(n),k,nn)
        sx(0,jpan+1,n,k,nn)      = s(iwn(ind(1,jpan,n)),k,nn)
        sx(ipan+1,jpan+1,n,k,nn) = s(ien(ind(ipan,jpan,n)),k,nn)
      end do               ! n loop
    end do                 ! k loop
  end do                   ! nn loop

! Loop over points that need to be calculated for other processes
  do ii=1,neighnum
    iproc=neighlistsend(ii)
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

      sc(-1,0,:) = sx(idel-1,jdel,n,k,:)
      sc(0,0,:)  = sx(idel  ,jdel,n,k,:)
      sc(1,0,:)  = sx(idel+1,jdel,n,k,:)
      sc(2,0,:)  = sx(idel+2,jdel,n,k,:)

      sc(-1,1,:) = sx(idel-1,jdel+1,n,k,:)
      sc(0,1,:)  = sx(idel  ,jdel+1,n,k,:)
      sc(1,1,:)  = sx(idel+1,jdel+1,n,k,:)
      sc(2,1,:)  = sx(idel+2,jdel+1,n,k,:)

      sc(0,-1,:) = sx(idel  ,jdel-1,n,k,:)
      sc(1,-1,:) = sx(idel+1,jdel-1,n,k,:)

      sc(0,2,:) = sx(idel  ,jdel+2,n,k,:)
      sc(1,2,:) = sx(idel+1,jdel+2,n,k,:)

      ncount=count(sc(:,:,1)>cxx)
      if (ncount>=12.and.lmode) then
        ! bi-cubic interpolation
        do nn=1,ntr
          r(2) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*sc(0,0,nn)-xxg*sc(-1,0,nn)/3.) &
               -xxg*(1.+xxg)*sc(2,0,nn)/3.)+xxg*(1.+xxg)*(2.-xxg)*sc(1,0,nn))/2.
          r(3) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*sc(0,1,nn)-xxg*sc(-1,1,nn)/3.) &
               -xxg*(1.+xxg)*sc(2,1,nn)/3.)+xxg*(1.+xxg)*(2.-xxg)*sc(1,1,nn))/2.
          r(1) = (1.-xxg)*sc(0,-1,nn) +xxg*sc(1,-1,nn)
          r(4) = (1.-xxg)*sc(0,2,nn) +xxg*sc(1,2,nn)

          sextra(iproc)%a(nn+(iq-1)*ntr) = ((1.-yyg)*((2.-yyg)* &
               ((1.+yyg)*r(2)-yyg*r(1)/3.)                      &
               -yyg*(1.+yyg)*r(4)/3.)                           &
               +yyg*(1.+yyg)*(2.-yyg)*r(3))/2.
        end do
      else
        ! bi-linear interpolation
        do nn=1,ntr
          scb=sc(0:1,0:1,nn)
          call lfill(scb,cxx)
          sc(0:1,0:1,nn)=scb
          aad=scb(1,1)-scb(0,1)-scb(1,0)+scb(0,0)
          aab=scb(1,0)-scb(0,0)
          aac=scb(0,1)-scb(0,0)
          sextra(iproc)%a(nn+(iq-1)*ntr)=aab*xxg+aac*yyg+aad*xxg*yyg+scb(0,0)
        end do
      end if
      do nn=1,ntr
        cmax=maxval(sc(0:1,0:1,nn))
        cmin=minval(sc(0:1,0:1,nn))
        sextra(iproc)%a(nn+(iq-1)*ntr)=min(max(sextra(iproc)%a(nn+(iq-1)*ntr),cmin),cmax)
      end do
    end do            ! iq loop
  end do              ! iproc loop

  call intssync_send(ntr)

  do iq=1,ifull
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
        if ( idel < 0 .or. idel > ipan .or. jdel < 0 .or. &
             jdel > jpan .or. n < 1 .or. n > npan ) then
          cycle      ! Will be calculated on another processor
        end if

        sc(0,-1,:) = sx(idel  ,jdel-1,n,k,:)
        sc(1,-1,:) = sx(idel+1,jdel-1,n,k,:)

        sc(-1,0,:) = sx(idel-1,jdel,n,k,:)
        sc(0,0,:)  = sx(idel  ,jdel,n,k,:)
        sc(1,0,:)  = sx(idel+1,jdel,n,k,:)
        sc(2,0,:)  = sx(idel+2,jdel,n,k,:)

        sc(-1,1,:) = sx(idel-1,jdel+1,n,k,:)
        sc(0,1,:)  = sx(idel  ,jdel+1,n,k,:)
        sc(1,1,:)  = sx(idel+1,jdel+1,n,k,:)
        sc(2,1,:)  = sx(idel+2,jdel+1,n,k,:)

        sc(0,2,:) = sx(idel  ,jdel+2,n,k,:)
        sc(1,2,:) = sx(idel+1,jdel+2,n,k,:)

        ncount=count(sc(:,:,1)>cxx)
        if (ncount>=12.and.lmode) then
          ! bi-cubic interpolation
          do nn=1,ntr
            r(2) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*sc(0,0,nn)-xxg*sc(-1,0,nn)/3.) &
                 -xxg*(1.+xxg)*sc(2,0,nn)/3.)+xxg*(1.+xxg)*(2.-xxg)*sc(1,0,nn))/2.
            r(3) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*sc(0,1,nn)-xxg*sc(-1,1,nn)/3.) &
                 -xxg*(1.+xxg)*sc(2,1,nn)/3.)+xxg*(1.+xxg)*(2.-xxg)*sc(1,1,nn))/2.
            r(1) = (1.-xxg)*sc(0,-1,nn)+xxg*sc(1,-1,nn)
            r(4) = (1.-xxg)*sc(0,2,nn) +xxg*sc(1,2,nn)

            s(iq,k,nn) = ((1.-yyg)*((2.-yyg)*          &
                 ((1.+yyg)*r(2)-yyg*r(1)/3.)           &
                 -yyg*(1.+yyg)*r(4)/3.)                &
                 +yyg*(1.+yyg)*(2.-yyg)*r(3))/2.
          end do
        else
          ! bi-linear interpolation
          do nn=1,ntr
            scb=sc(0:1,0:1,nn)
            call lfill(scb,cxx)
            sc(0:1,0:1,nn)=scb
            aad=scb(1,1)-scb(0,1)-scb(1,0)+scb(0,0)
            aab=scb(1,0)-scb(0,0)
            aac=scb(0,1)-scb(0,0)
            s(iq,k,nn)=aab*xxg+aac*yyg+aad*xxg*yyg+scb(0,0)
          end do
        end if
        do nn=1,ntr
          cmax=maxval(sc(0:1,0:1,nn))
          cmin=minval(sc(0:1,0:1,nn))
          s(iq,k,nn)=min(max(s(iq,k,nn),cmin),cmax)
        end do
      end do ! k loop
    end if   ! wtr
  end do     ! iq loop
       
!========================   end of intsch=1 section ====================
else     ! if(intsch==1)then
!======================== start of intsch=2 section ====================
!       this is intsc           NS interps done first
!       first extend s arrays into sx - this one -1:il+2 & -1:il+2

  do nn=1,ntr
    do k=1,kx
      do n=1,npan         ! first simple copy into larger array
        do j=1,jpan
          do i=1,ipan
            sx(i,j,n,k,nn) = s(ind(i,j,n),k,nn)
          end do         ! i loop
          sx(0,j,n,k,nn)      = s(iw(ind(1,j,n)),k,nn)
          sx(-1,j,n,k,nn)     = s(iww(ind(1,j,n)),k,nn)
          sx(ipan+1,j,n,k,nn) = s(ie(ind(ipan,j,n)),k,nn)
          sx(ipan+2,j,n,k,nn) = s(iee(ind(ipan,j,n)),k,nn)
        end do            ! j loop
        do i=1,ipan
          sx(i,0,n,k,nn)      = s(is(ind(i,1,n)),k,nn)
          sx(i,-1,n,k,nn)     = s(iss(ind(i,1,n)),k,nn)
          sx(i,jpan+1,n,k,nn) = s(in(ind(i,jpan,n)),k,nn)
          sx(i,jpan+2,n,k,nn) = s(inn(ind(i,jpan,n)),k,nn)
        end do            ! i loop
!        for ns interpolation, sometimes need (different from ew):
!            (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
!          (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)

        sx(-1,0,n,k,nn)          = s(lsww(n),k,nn)
        sx(0,0,n,k,nn)           = s(isw(ind(1,1,n)),k,nn)
        sx(0,-1,n,k,nn)          = s(lssw(n),k,nn)
        sx(ipan+2,0,n,k,nn)      = s(lsee(n),k,nn)
        sx(ipan+1,-1,n,k,nn)     = s(lsse(n),k,nn)
        sx(-1,jpan+1,n,k,nn)     = s(lnww(n),k,nn)
        sx(0,jpan+1,n,k,nn)      = s(inw(ind(1,jpan,n)),k,nn)
        sx(0,jpan+2,n,k,nn)      = s(lnnw(n),k,nn)
        sx(ipan+2,jpan+1,n,k,nn) = s(lnee(n),k,nn)
        sx(ipan+1,jpan+2,n,k,nn) = s(lnne(n),k,nn)
        sx(ipan+1,0,n,k,nn)      = s(ise(ind(ipan,1,n)),k,nn)
        sx(ipan+1,jpan+1,n,k,nn) = s(ine(ind(ipan,jpan,n)),k,nn)
      end do               ! n loop
    end do                 ! k loop
  end do                   ! nn loop

! For other processes
  do ii=1,neighnum
    iproc=neighlistsend(ii)
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

      sc(0,-1,:) = sx(idel,jdel-1,n,k,:)
      sc(0,0,:)  = sx(idel,jdel  ,n,k,:)
      sc(0,1,:)  = sx(idel,jdel+1,n,k,:)
      sc(0,2,:)  = sx(idel,jdel+2,n,k,:)

      sc(1,-1,:) = sx(idel+1,jdel-1,n,k,:)
      sc(1,0,:)  = sx(idel+1,jdel  ,n,k,:)
      sc(1,1,:)  = sx(idel+1,jdel+1,n,k,:)
      sc(1,2,:)  = sx(idel+1,jdel+2,n,k,:)

      sc(-1,0,:) = sx(idel-1,jdel  ,n,k,:)
      sc(-1,1,:) = sx(idel-1,jdel+1,n,k,:)
        
      sc(2,0,:) = sx(idel+2,jdel  ,n,k,:)
      sc(2,1,:) = sx(idel+2,jdel+1,n,k,:)

      ncount=count(sc(:,:,1)>cxx)
      if (ncount>=12.and.lmode) then
        ! bi-cubic interpolation
        do nn=1,ntr
          r(2) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*sc(0,0,nn)-yyg*sc(0,-1,nn)/3.) &
               -yyg*(1.+yyg)*sc(0,2,nn)/3.)+yyg*(1.+yyg)*(2.-yyg)*sc(0,1,nn))/2.
          r(3) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*sc(1,0,nn)-yyg*sc(1,-1,nn)/3.) &
               -yyg*(1.+yyg)*sc(1,2,nn)/3.)+yyg*(1.+yyg)*(2.-yyg)*sc(1,1,nn))/2.
          r(1) = (1.-yyg)*sc(-1,0,nn) +yyg*sc(-1,1,nn)
          r(4) = (1.-yyg)*sc(2,0,nn) +yyg*sc(2,1,nn)
          sextra(iproc)%a(nn+(iq-1)*ntr) = ((1.-xxg)*((2.-xxg)* &
               ((1.+xxg)*r(2)-xxg*r(1)/3.)           &
               -xxg*(1.+xxg)*r(4)/3.)                &
               +xxg*(1.+xxg)*(2.-xxg)*r(3))/2.
        end do
      else
        ! bi-linear interpolation
        do nn=1,ntr
          scb=sc(0:1,0:1,nn)
          call lfill(scb,cxx)
          sc(0:1,0:1,nn)=scb        
          aad=scb(1,1)-scb(0,1)-scb(1,0)+scb(0,0)
          aab=scb(1,0)-scb(0,0)
          aac=scb(0,1)-scb(0,0)
          sextra(iproc)%a(nn+(iq-1)*ntr)=aab*xxg+aac*yyg+aad*xxg*yyg+scb(0,0)
        end do
      end if
      do nn=1,ntr
        cmax=maxval(sc(0:1,0:1,nn))
        cmin=minval(sc(0:1,0:1,nn))
        sextra(iproc)%a(nn+(iq-1)*ntr)=min(max(sextra(iproc)%a(nn+(iq-1)*ntr),cmin),cmax)
      end do
    end do            ! iq loop
  end do              ! iproc

  call intssync_send(ntr)

  do iq=1,ifull
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
        if ( idel < 0 .or. idel > ipan .or. jdel < 0 .or. &
             jdel > jpan .or. n < 1 .or. n > npan ) then
          cycle      ! Will be calculated on another processor
        end if

        sc(0,-1,:) = sx(idel,jdel-1,n,k,:)
        sc(0,0,:)  = sx(idel,jdel  ,n,k,:)
        sc(0,1,:)  = sx(idel,jdel+1,n,k,:)
        sc(0,2,:)  = sx(idel,jdel+2,n,k,:)

        sc(1,-1,:) = sx(idel+1,jdel-1,n,k,:)
        sc(1,0,:)  = sx(idel+1,jdel  ,n,k,:)
        sc(1,1,:)  = sx(idel+1,jdel+1,n,k,:)
        sc(1,2,:)  = sx(idel+1,jdel+2,n,k,:)

        sc(-1,0,:) = sx(idel-1,jdel  ,n,k,:)
        sc(-1,1,:) = sx(idel-1,jdel+1,n,k,:)
        
        sc(2,0,:) = sx(idel+2,jdel  ,n,k,:)
        sc(2,1,:) = sx(idel+2,jdel+1,n,k,:)
        
        ncount=count(sc(:,:,1)>cxx)
        if (ncount>=12.and.lmode) then
          ! bi-cubic interpolation
          do nn=1,ntr
            r(2) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*sc(0,0,nn)-yyg*sc(0,-1,nn)/3.) &
                 -yyg*(1.+yyg)*sc(0,2,nn)/3.)+yyg*(1.+yyg)*(2.-yyg)*sc(0,1,nn))/2.
            r(3) = ((1.-yyg)*((2.-yyg)*((1.+yyg)*sc(1,0,nn)-yyg*sc(1,-1,nn)/3.) &
                 -yyg*(1.+yyg)*sc(1,2,nn)/3.)+yyg*(1.+yyg)*(2.-yyg)*sc(1,1,nn))/2.
            r(1) = (1.-yyg)*sc(-1,0,nn)+yyg*sc(-1,1,nn)
            r(4) = (1.-yyg)*sc(2,0,nn) +yyg*sc(2,1,nn)

            s(iq,k,nn) = ((1.-xxg)*((2.-xxg)*       &
                 ((1.+xxg)*r(2)-xxg*r(1)/3.)        &
                 -xxg*(1.+xxg)*r(4)/3.)             &
                 +xxg*(1.+xxg)*(2.-xxg)*r(3))/2.
          end do
        else
          ! bi-linear interpolation
          do nn=1,ntr
            scb=sc(0:1,0:1,nn)
            call lfill(scb,cxx)
            sc(0:1,0:1,nn)=scb        
            aad=scb(1,1)-scb(0,1)-scb(1,0)+scb(0,0)
            aab=scb(1,0)-scb(0,0)
            aac=scb(0,1)-scb(0,0)
            s(iq,k,nn)=aab*xxg+aac*yyg+aad*xxg*yyg+scb(0,0)
          end do
        end if
        do nn=1,ntr
          cmax=maxval(sc(0:1,0:1,nn))
          cmin=minval(sc(0:1,0:1,nn))
          s(iq,k,nn)=min(max(s(iq,k,nn),cmin),cmax)
        end do
      end do
    end if
  end do

endif                     ! (intsch==1) .. else ..
!========================   end of intsch=1 section ====================

call intssync_recv(s)

do iq=1,ifull
  if (.not.wtr(iq)) then
    s(iq,:,:)=ssav(iq,:,:)
  end if
end do

where (s(1:ifull,:,:)<cxx+10.)
  s(1:ifull,:,:)=ssav
end where

return
end subroutine mlob2intsb

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
    vdot1 = (vec1x*cou(1:ifull,k) + vec1y*cov(1:ifull,k) + vec1z*cow(1:ifull,k))/denb
    vdot2 = (vec2x*cou(1:ifull,k) + vec2y*cov(1:ifull,k) + vec2z*cow(1:ifull,k))/denb
    cou(1:ifull,k) = vdot1*vec1x + vdot2*vec3x
    cov(1:ifull,k) = vdot1*vec1y + vdot2*vec3y
    cow(1:ifull,k) = vdot1*vec1z + vdot2*vec3z
  end where
end do ! k

return
end subroutine mlorot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Stagger u and v
! Modified to include zero velocity next to land points
subroutine mlostaguv(u,v,uout,vout)

use cc_mpi
use indices_m
use mlo

implicit none

include 'newmpar.h'
include 'parm.h'

integer k,itn,kx
real, dimension(:,:), intent(in) :: u
real, dimension(ifull,size(u,2)), intent(in) :: v
real, dimension(ifull,size(u,2)), intent(out) :: uout,vout
real, dimension(ifull+iextra,size(u,2)) :: uin,vin
real, dimension(ifull+iextra,size(u,2)) :: ua,va
real, dimension(ifull+iextra,size(u,2)) :: ud,vd
real, dimension(:,:), allocatable, save :: wtul,wtvl
real, dimension(:,:), allocatable, save :: wtur,wtvr
real, dimension(:,:), allocatable, save :: dtul,dtvl
real, dimension(:,:), allocatable, save :: dtur,dtvr
real, dimension(ifull) :: nud,nvd
real, dimension(ifull,0:3) :: nwtu,nwtv
real, dimension(:), allocatable, save :: stul,stvl
real, dimension(:), allocatable, save :: stur,stvr
logical, dimension(ifull) :: euetest,euwtest,evntest,evstest
logical, dimension(ifull) :: euewtest,evnstest,eutest,evtest
logical ltest

call start_log(ocnstag_begin)

if (.not.allocated(wtul)) then
  allocate(wtul(ifull+iextra,0:3),wtvl(ifull+iextra,0:3))
  allocate(wtur(ifull+iextra,0:3),wtvr(ifull+iextra,0:3))
  allocate(dtul(ifull,3),dtvl(ifull,3))
  allocate(dtur(ifull,3),dtvr(ifull,3))
  allocate(stul(ifull),stvl(ifull))
  allocate(stur(ifull),stvr(ifull))

  ! assign land arrays
  eutest=eeu(1:ifull)>0.5
  evtest=eev(1:ifull)>0.5
  euetest=eutest      .and.eeu(ieu)>0.5
  euwtest=eeu(iwu)>0.5.and.eutest
  evntest=evtest      .and.eev(inv)>0.5
  evstest=eev(isv)>0.5.and.evtest
  euewtest=euetest.and.euwtest
  evnstest=evntest.and.evstest

  ! assign weights (left)

! |   *   | X E   |  EE  |     unstaggered
! W       * X     E            staggered
  where (euewtest)
    wtul(1:ifull,0)=1.
    wtul(1:ifull,1)=-0.5
    wtul(1:ifull,2)=-0.1
    wtul(1:ifull,3)=0.
    dtul(:,1)=0.1
    dtul(:,2)=1.
    dtul(:,3)=0.5
    !ud(1:ifull,k)=uin(ieeu,k)*0.1+uin(ieu,k)+uin(1:ifull,k)*0.5
    !uin(1:ifull,k)=ud(:,k)-ua(ieu,k)*0.5-ua(iwu,k)*0.1

! #   *   | X E   |  EE  |     unstaggered
! 0       * X     E            staggered
  elsewhere (euetest)
    wtul(1:ifull,0)=1.
    wtul(1:ifull,1)=-0.5
    wtul(1:ifull,2)=0.
    wtul(1:ifull,3)=0.
    dtul(:,1)=0.1
    dtul(:,2)=1.
    dtul(:,3)=0.5

! |   *   | X E   #  ##  #     unstaggered
!         * X     0  ##  #     staggered
  elsewhere (euwtest)
    wtul(1:ifull,0)=1.
    wtul(1:ifull,1)=0.
    wtul(1:ifull,2)=0.
    wtul(1:ifull,3)=0.
    dtul(:,1)=0.
    dtul(:,2)=1.
    dtul(:,3)=1./3.

! #   *   |   E   #  ##  #     unstaggered
! #       *       #  ##  #     staggered
  elsewhere (eutest)
    wtul(1:ifull,0)=1.
    wtul(1:ifull,1)=0.
    wtul(1:ifull,2)=0.
    wtul(1:ifull,3)=0.
    dtul(:,1)=0.
    dtul(:,2)=0.5
    dtul(:,3)=0.5

  elsewhere
    wtul(1:ifull,0)=0.
    wtul(1:ifull,1)=0.
    wtul(1:ifull,2)=0.
    wtul(1:ifull,3)=0.
    dtul(:,1)=0.
    dtul(:,2)=0.
    dtul(:,3)=0.
            
  end where
  where (evnstest)
    wtvl(1:ifull,0)=1.
    wtvl(1:ifull,1)=-0.5
    wtvl(1:ifull,2)=-0.1
    wtvl(1:ifull,3)=0.
    dtvl(:,1)=0.1
    dtvl(:,2)=1.
    dtvl(:,3)=0.5
    !vin(1:ifull,k)=vd(:,k)-va(inv,k)*0.5-va(isv,k)*0.1
  elsewhere (evntest)
    wtvl(1:ifull,0)=1.
    wtvl(1:ifull,1)=-0.5
    wtvl(1:ifull,2)=0.
    wtvl(1:ifull,3)=0.
    dtvl(:,1)=0.1
    dtvl(:,2)=1.
    dtvl(:,3)=0.5
  elsewhere (evstest)
    wtvl(1:ifull,0)=1.
    wtvl(1:ifull,1)=0.
    wtvl(1:ifull,2)=0.
    wtvl(1:ifull,3)=0.
    dtvl(:,1)=0. 
    dtvl(:,2)=1.
    dtvl(:,3)=1./3.
  elsewhere (evtest)
    wtvl(1:ifull,0)=1.
    wtvl(1:ifull,1)=0.
    wtvl(1:ifull,2)=0.
    wtvl(1:ifull,3)=0.
    dtvl(:,1)=0.
    dtvl(:,2)=0.5
    dtvl(:,3)=0.5
  elsewhere
    wtvl(1:ifull,0)=0.
    wtvl(1:ifull,1)=0.
    wtvl(1:ifull,2)=0.
    wtvl(1:ifull,3)=0.
    dtvl(:,1)=0.
    dtvl(:,2)=0.
    dtvl(:,3)=0.
  end where

  ! Apply JLM's preconditioner
  call boundsuv(wtul,wtvl,stag=-10)
  where (abs(wtul(ieu,0))>1.E-4.and.abs(wtul(1:ifull,1))>1.E-4)
    stul=-wtul(1:ifull,1)/wtul(ieu,0)
    nwtu(1:ifull,0)=wtul(1:ifull,0)+stul*wtul(ieu,2)
    nwtu(1:ifull,1)=wtul(1:ifull,1)+stul*wtul(ieu,0)
    nwtu(1:ifull,2)=wtul(1:ifull,2)
    nwtu(1:ifull,3)=wtul(1:ifull,3)-stul*wtul(ieu,1)
  elsewhere
    stul=0.
    nwtu(1:ifull,0)=wtul(1:ifull,0)
    nwtu(1:ifull,1)=wtul(1:ifull,1)
    nwtu(1:ifull,2)=wtul(1:ifull,2)
    nwtu(1:ifull,3)=wtul(1:ifull,3)
  end where
  where (abs(wtvl(inv,0))>1.E-4.and.abs(wtvl(1:ifull,1))>1.E-4)
    stvl=-wtvl(1:ifull,1)/wtvl(inv,0)
    nwtv(1:ifull,0)=wtvl(1:ifull,0)+stvl*wtvl(inv,2)
    nwtv(1:ifull,1)=wtvl(1:ifull,1)+stvl*wtvl(inv,0)
    nwtv(1:ifull,2)=wtvl(1:ifull,2)
    nwtv(1:ifull,3)=wtvl(1:ifull,3)-stvl*wtvl(inv,1)
  elsewhere
    stvl=0.
    nwtv(1:ifull,0)=wtvl(1:ifull,0)
    nwtv(1:ifull,1)=wtvl(1:ifull,1)
    nwtv(1:ifull,2)=wtvl(1:ifull,2)
    nwtv(1:ifull,3)=wtvl(1:ifull,3)
  end where
  where (abs(wtul(1:ifull,0))<1.E-4)
    wtul(1:ifull,0)=1.
    wtul(1:ifull,1)=0.
    wtul(1:ifull,2)=0.
    wtul(1:ifull,3)=0.
  elsewhere
    wtul(1:ifull,0)=nwtu(1:ifull,0)
    wtul(1:ifull,1)=nwtu(1:ifull,1)
    wtul(1:ifull,2)=nwtu(1:ifull,2)
    wtul(1:ifull,3)=nwtu(1:ifull,3)
  end where
  where (abs(wtvl(1:ifull,0))<1.E-4)
    wtvl(1:ifull,0)=1.
    wtvl(1:ifull,1)=0.
    wtvl(1:ifull,2)=0.
    wtvl(1:ifull,3)=0.
  elsewhere
    wtvl(1:ifull,0)=nwtv(1:ifull,0)
    wtvl(1:ifull,1)=nwtv(1:ifull,1)
    wtvl(1:ifull,2)=nwtv(1:ifull,2)
    wtvl(1:ifull,3)=nwtv(1:ifull,3)
  end where

  ! normalise
  call boundsuv(wtul,wtvl,stag=-10)
  do k=1,3
    dtul(:,k)=dtul(:,k)/wtul(1:ifull,0)
    dtvl(:,k)=dtvl(:,k)/wtvl(1:ifull,0)
    wtul(1:ifull,k)=wtul(1:ifull,k)/wtul(1:ifull,0)
    wtvl(1:ifull,k)=wtvl(1:ifull,k)/wtvl(1:ifull,0)
  end do
  stul=stul*wtul(ieu,0)/wtul(1:ifull,0)
  stvl=stvl*wtvl(inv,0)/wtvl(1:ifull,0)
  wtul(1:ifull,0)=1.
  wtvl(1:ifull,0)=1.

  ! assign weights (right)

! |   W   |   * X |  E   |     unstaggered
!         W     X *      E     staggered
  where (euewtest)
    wtur(1:ifull,0)=1.
    wtur(1:ifull,1)=-0.1
    wtur(1:ifull,2)=-0.5
    wtur(1:ifull,3)=0.
    dtur(:,1)=0.1
    dtur(:,2)=1.
    dtur(:,3)=0.5
    !ud(1:ifull,k)=uin(iwu,k)*0.1+uin(1:ifull,k)+uin(ieu,k)*0.5
    !uin(1:ifull,k)=ud(:,k)-ua(iwu,k)*0.5-ua(ieu,k)*0.1

! |   W   |   * X |  E   #     unstaggered
!         W     X *      0     staggered
  elsewhere (euwtest)
    wtur(1:ifull,0)=1.
    wtur(1:ifull,1)=0.
    wtur(1:ifull,2)=-0.5
    wtur(1:ifull,3)=0.
    dtur(:,1)=0.1
    dtur(:,2)=1.
    dtur(:,3)=0.5

! #  ##   #   * X |   E   |     unstaggered
! #  ##   0     X *             staggered
  elsewhere (euetest)
    wtur(1:ifull,0)=1.
    wtur(1:ifull,1)=0.
    wtur(1:ifull,2)=0.
    wtur(1:ifull,3)=0.
    dtur(:,1)=0.
    dtur(:,2)=1.
    dtur(:,3)=1./3.

! #  ##   #   *   |  E   #     unstaggered
! #  ##   #       *      #     staggered
  elsewhere (eutest)
    wtur(1:ifull,0)=1.
    wtur(1:ifull,1)=0.
    wtur(1:ifull,2)=0.
    wtur(1:ifull,3)=0.
    dtur(:,1)=0.
    dtur(:,2)=0.5
    dtur(:,3)=0.5
  elsewhere
    wtur(1:ifull,0)=0.
    wtur(1:ifull,1)=0.
    wtur(1:ifull,2)=0.
    wtur(1:ifull,3)=0.
    dtur(:,1)=0.
    dtur(:,2)=0.
    dtur(:,3)=0.
      
  end where
  where (evnstest)
    wtvr(1:ifull,0)=1.
    wtvr(1:ifull,1)=-0.1
    wtvr(1:ifull,2)=-0.5
    wtvr(1:ifull,3)=0.
    dtvr(:,1)=0.1
    dtvr(:,2)=1.
    dtvr(:,3)=0.5
  elsewhere (evstest)
    wtvr(1:ifull,0)=1.
    wtvr(1:ifull,1)=0.
    wtvr(1:ifull,2)=-0.5
    wtvr(1:ifull,3)=0.
    dtvr(:,1)=0.1
    dtvr(:,2)=1.
    dtvr(:,3)=0.5
  elsewhere (evntest)
    wtvr(1:ifull,0)=1.
    wtvr(1:ifull,1)=0.
    wtvr(1:ifull,2)=0.
    wtvr(1:ifull,3)=0.
    dtvr(:,1)=0.
    dtvr(:,2)=1.
    dtvr(:,3)=1./3.
  elsewhere (evtest)
    wtvr(1:ifull,0)=1.
    wtvr(1:ifull,1)=0.
    wtvr(1:ifull,2)=0.
    wtvr(1:ifull,3)=0.
    dtvr(:,1)=0.
    dtvr(:,2)=0.5
    dtvr(:,3)=0.5
  elsewhere
    wtvr(1:ifull,0)=0.
    wtvr(1:ifull,1)=0.
    wtvr(1:ifull,2)=0.
    wtvr(1:ifull,3)=0.
    dtvr(:,1)=0.
    dtvr(:,2)=0.
    dtvr(:,3)=0.
  end where

  ! Apply JLM's preconditioner
  call boundsuv(wtur,wtvr,stag=-9)
  where (abs(wtur(iwu,0))>1.E-4.and.abs(wtur(1:ifull,2))>1.E-4)
    stur=-wtur(1:ifull,2)/wtur(iwu,0)
    nwtu(1:ifull,0)=wtur(1:ifull,0)+stur*wtur(iwu,1)
    nwtu(1:ifull,1)=wtur(1:ifull,1)
    nwtu(1:ifull,2)=wtur(1:ifull,2)+stur*wtur(iwu,0)
    nwtu(1:ifull,3)=wtur(1:ifull,3)-stur*wtur(iwu,2)
  elsewhere
    stur=0.
    nwtu(1:ifull,0)=wtur(1:ifull,0)
    nwtu(1:ifull,1)=wtur(1:ifull,1)
    nwtu(1:ifull,2)=wtur(1:ifull,2)
    nwtu(1:ifull,3)=wtur(1:ifull,3)
  end where
  where (abs(wtvr(isv,0))>1.E-4.and.abs(wtvr(1:ifull,2))>1.E-4)
    stvr=-wtvr(1:ifull,2)/wtvr(isv,0)
    nwtv(1:ifull,0)=wtvr(1:ifull,0)+stvr*wtvr(isv,1)
    nwtv(1:ifull,1)=wtvr(1:ifull,1)
    nwtv(1:ifull,2)=wtvr(1:ifull,2)+stvr*wtvr(isv,0)
    nwtv(1:ifull,3)=wtvr(1:ifull,3)-stvr*wtvr(isv,2)
  elsewhere
    stvr=0.
    nwtv(1:ifull,0)=wtvr(1:ifull,0)
    nwtv(1:ifull,1)=wtvr(1:ifull,1)
    nwtv(1:ifull,2)=wtvr(1:ifull,2)
    nwtv(1:ifull,3)=wtvr(1:ifull,3)
  end where
  where (abs(wtur(1:ifull,0))<1.E-4)
    wtur(1:ifull,0)=1.
    wtur(1:ifull,1)=0.
    wtur(1:ifull,2)=0.
    wtur(1:ifull,3)=0.
  elsewhere
    wtur(1:ifull,0)=nwtu(1:ifull,0)
    wtur(1:ifull,1)=nwtu(1:ifull,1)
    wtur(1:ifull,2)=nwtu(1:ifull,2)
    wtur(1:ifull,3)=nwtu(1:ifull,3)
  end where
  where (abs(wtvr(1:ifull,0))<1.E-4)
    wtvr(1:ifull,0)=1.
    wtvr(1:ifull,1)=0.
    wtvr(1:ifull,2)=0.
    wtvr(1:ifull,3)=0.
  elsewhere
    wtvr(1:ifull,0)=nwtv(1:ifull,0)
    wtvr(1:ifull,1)=nwtv(1:ifull,1)
    wtvr(1:ifull,2)=nwtv(1:ifull,2)
    wtvr(1:ifull,3)=nwtv(1:ifull,3)
  end where

  ! normalise
  call boundsuv(wtur,wtvr,stag=-9)
  do k=1,3
    dtur(:,k)=dtur(:,k)/wtur(1:ifull,0)
    dtvr(:,k)=dtvr(:,k)/wtvr(1:ifull,0)
    wtur(1:ifull,k)=wtur(1:ifull,k)/wtur(1:ifull,0)
    wtvr(1:ifull,k)=wtvr(1:ifull,k)/wtvr(1:ifull,0)
  end do
  stur=stur*wtur(iwu,0)/wtur(1:ifull,0)
  stvr=stvr*wtvr(isv,0)/wtvr(1:ifull,0)
  wtur(1:ifull,0)=1.
  wtvr(1:ifull,0)=1.

end if

kx=size(u,2)

do k=1,kx
  uin(1:ifull,k)=u(:,k)*ee(1:ifull)
  vin(1:ifull,k)=v(:,k)*ee(1:ifull)
end do

if (mstagf==0) then
  ltest=.true.
else if (mstagf<0) then
  ltest=.false.
else
  ! using ktau-1 ensures that the staggering is relative to the C grid
  ltest=mod(ktau-koff-nstagoffmlo,2*mstagf)<mstagf
end if

if (ltest) then

  call boundsuv(uin,vin,stag=1)
  do k=1,kx
    ud(1:ifull,k)=dtul(:,1)*uin(ieeu,k)+dtul(:,2)*uin(ieu,k)+dtul(:,3)*uin(1:ifull,k)
    vd(1:ifull,k)=dtvl(:,1)*vin(innv,k)+dtvl(:,2)*vin(inv,k)+dtvl(:,3)*vin(1:ifull,k)
  end do
  
  call boundsuv(ud,vd,stag=-10)
  do k=1,kx
    ! Apply JLM's preconditioner
    nud=ud(1:ifull,k)-stul*ud(ieu,k)
    nvd=vd(1:ifull,k)-stvl*vd(inv,k)
    ud(1:ifull,k)=nud
    vd(1:ifull,k)=nvd

    ! 1st guess
    !ua(1:ifull,k)=(uin(1:ifull,k)+uin(ieu,k))*0.5*eeu(1:ifull)
    !va(1:ifull,k)=(vin(1:ifull,k)+vin(inv,k))*0.5*eev(1:ifull)
    ua(1:ifull,k)=nud
    va(1:ifull,k)=nvd
  end do

  ! There are many ways to handle staggering near coastlines.
  ! This version supports the following properties:
  ! - Staggering is exactly reversible for the staggered (C grid) frame
  ! - Zero flux at land boundaries
  ! - The wave amplitude should be preserved

  do itn=1,itnmax        ! each loop is a double iteration
    call boundsuv(ua,va,stag=2)
    do k=1,kx
      uin(1:ifull,k)=ud(1:ifull,k)+wtul(:,1)*ua(ieu,k)+wtul(:,2)*ua(iwu,k)+wtul(:,3)*ua(ieeu,k)
      vin(1:ifull,k)=vd(1:ifull,k)+wtvl(:,1)*va(inv,k)+wtvl(:,2)*va(isv,k)+wtvl(:,3)*va(innv,k)
    end do
    call boundsuv(uin,vin,stag=2)
    do k=1,kx
      ua(1:ifull,k)=ud(1:ifull,k)+wtul(:,1)*uin(ieu,k)+wtul(:,2)*uin(iwu,k)+wtul(:,3)*uin(ieeu,k)
      va(1:ifull,k)=vd(1:ifull,k)+wtvl(:,1)*vin(inv,k)+wtvl(:,2)*vin(isv,k)+wtvl(:,3)*vin(innv,k)
    end do
  end do                 ! itn=1,itnmax

else

  call boundsuv(uin,vin)
  do k=1,kx
    ud(1:ifull,k)=dtur(:,1)*uin(iwu,k)+dtur(:,2)*uin(1:ifull,k)+dtur(:,3)*uin(ieu,k)
    vd(1:ifull,k)=dtvr(:,1)*vin(isv,k)+dtvr(:,2)*vin(1:ifull,k)+dtvr(:,3)*vin(inv,k)
  end do

  call boundsuv(ud,vd,stag=-9)
  do k=1,kx
    ! Apply JLM's preconditioner
    nud=ud(1:ifull,k)-stur*ud(iwu,k)
    nvd=vd(1:ifull,k)-stvr*vd(isv,k)
    ud(1:ifull,k)=nud
    vd(1:ifull,k)=nvd

    ! 1st guess
    !ua(1:ifull,k)=(uin(1:ifull,k)+uin(ieu,k))*0.5*eeu(1:ifull)
    !va(1:ifull,k)=(vin(1:ifull,k)+vin(inv,k))*0.5*eev(1:ifull)
    ua(1:ifull,k)=nud
    va(1:ifull,k)=nvd
  end do

  do itn=1,itnmax        ! each loop is a double iteration
    call boundsuv(ua,va,stag=3)
    do k=1,kx
      uin(1:ifull,k)=ud(1:ifull,k)+wtur(:,1)*ua(ieu,k)+wtur(:,2)*ua(iwu,k)+wtur(:,3)*ua(iwwu,k)
      vin(1:ifull,k)=vd(1:ifull,k)+wtvr(:,1)*va(inv,k)+wtvr(:,2)*va(isv,k)+wtvr(:,3)*va(issv,k)
    end do
    call boundsuv(uin,vin,stag=3)
    do k=1,kx
      ua(1:ifull,k)=ud(1:ifull,k)+wtur(:,1)*uin(ieu,k)+wtur(:,2)*uin(iwu,k)+wtur(:,3)*uin(iwwu,k)
      va(1:ifull,k)=vd(1:ifull,k)+wtvr(:,1)*vin(inv,k)+wtvr(:,2)*vin(isv,k)+wtvr(:,3)*vin(issv,k)
    end do
  end do                 ! itn=1,itnmax

end if

uout=ua(1:ifull,:)
vout=va(1:ifull,:)

call end_log(ocnstag_end)

return
end subroutine mlostaguv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Unstagger u and v
! Modified to include zero velocity over land points
subroutine mlounstaguv(u,v,uout,vout,toff)

use cc_mpi
use indices_m
use mlo

implicit none

include 'newmpar.h'
include 'parm.h'

integer, intent(in), optional :: toff
integer k,itn,kx,zoff
real, dimension(:,:), intent(in) :: u
real, dimension(ifull,size(u,2)), intent(in) :: v
real, dimension(ifull,size(u,2)), intent(out) :: uout,vout
real, dimension(ifull+iextra,size(u,2)) :: uin,vin
real, dimension(ifull+iextra,size(u,2)) :: ua,va
real, dimension(ifull+iextra,size(u,2)) :: ud,vd
real, dimension(ifull) :: nud,nvd
real, dimension(ifull,0:3) :: nwtu,nwtv
real, dimension(:,:), allocatable, save :: wtul,wtvl
real, dimension(:,:), allocatable, save :: wtur,wtvr
real, dimension(:,:), allocatable, save :: dtul,dtvl
real, dimension(:,:), allocatable, save :: dtur,dtvr
real, dimension(:), allocatable, save :: stul,stvl
real, dimension(:), allocatable, save :: stur,stvr
logical, dimension(ifull) :: eutest,evtest
logical, dimension(ifull) :: eetest,ewtest,entest,estest
logical, dimension(ifull) :: euetest,euwtest,evntest,evstest
logical, dimension(ifull) :: euewtest,evnstest
logical, dimension(ifull) :: eeetest,ewwtest,enntest,esstest
logical ltest

call start_log(ocnstag_begin)

if (.not.allocated(wtul)) then
  allocate(wtul(ifull+iextra,0:3),wtvl(ifull+iextra,0:3))
  allocate(wtur(ifull+iextra,0:3),wtvr(ifull+iextra,0:3))
  allocate(dtul(ifull,3),dtvl(ifull,3))
  allocate(dtur(ifull,3),dtvr(ifull,3))
  allocate(stul(ifull),stvl(ifull))
  allocate(stur(ifull),stvr(ifull))
  
  ! assign land arrays
  eetest=ee(1:ifull)*ee(ie)>0.5
  ewtest=ee(iw)*ee(1:ifull)>0.5
  entest=ee(1:ifull)*ee(in)>0.5
  estest=ee(is)*ee(1:ifull)>0.5

  eutest=eeu(iwu)>0.5
  evtest=eev(isv)>0.5
  euetest=eutest.and.eeu(1:ifull)>0.5
  evntest=evtest.and.eev(1:ifull)>0.5
  euwtest=eutest.and.eeu(iwwu)>0.5
  evstest=evtest.and.eev(issv)>0.5
  euewtest=euetest.and.euwtest
  evnstest=evntest.and.evstest
  eeetest=eetest.and.ee(iee)>0.5
  enntest=entest.and.ee(inn)>0.5

  ! assign weights (left)

!  |   W   | X *   |  E   |     unstaggered
! WW       W X     *            staggered
  where (euewtest)
    wtul(1:ifull,0)=1.
    wtul(1:ifull,1)=-0.1
    wtul(1:ifull,2)=-0.5
    wtul(1:ifull,3)=0.
    dtul(:,1)=0.1
    dtul(:,2)=1.
    dtul(:,3)=0.5
    !ud(1:ifull,k)=uin(iwwu,k)*0.1+uin(iwu,k)+uin(1:ifull,k)*0.5
    !uin(1:ifull,k)=ud(:,k)-ua(ieu,k)*0.1-ua(iwu,k)*0.5

! ##   W   | X *   |  E   |     unstaggered
! #0       W X     *            staggered
  elsewhere (euetest)
    wtul(1:ifull,0)=1.
    wtul(1:ifull,1)=-0.1
    wtul(1:ifull,2)=-0.5
    wtul(1:ifull,3)=0.
    dtul(:,1)=0.
    dtul(:,2)=1.
    dtul(:,3)=0.5
      
!  |   W   | X *   #  ##  #     unstaggered
!          W X     0  ##  #     staggered
  elsewhere (euwtest)
    wtul(1:ifull,0)=1.
    wtul(1:ifull,1)=0.
    wtul(1:ifull,2)=-1./3.
    wtul(1:ifull,3)=0.
    dtul(:,1)=0.
    dtul(:,2)=1.
    dtul(:,3)=0.

! ##   ##  #   * X |   E   |     unstaggered
! ##   ##  0     X *             staggered
  elsewhere (eeetest)
    wtul(1:ifull,0)=1.
    wtul(1:ifull,1)=-1./3.
    wtul(1:ifull,2)=0.
    wtul(1:ifull,3)=0.
    dtul(:,1)=0.
    dtul(:,2)=0.
    dtul(:,3)=1.

! ##   ##  #   *   |   E    #     unstaggered
! ##   ##  #       *       #     staggered
  elsewhere (eetest)
    wtul(1:ifull,0)=1.
    wtul(1:ifull,1)=0.
    wtul(1:ifull,2)=0.
    wtul(1:ifull,3)=0.
    dtul(:,1)=0.
    dtul(:,2)=0.
    dtul(:,3)=1.

! ##   W   |   *   #  ##  #     unstaggered
! ##       W       #  ##  #     staggered
  elsewhere (eutest)
    wtul(1:ifull,0)=1.
    wtul(1:ifull,1)=0.
    wtul(1:ifull,2)=0.
    wtul(1:ifull,3)=0.
    dtul(:,1)=0.
    dtul(:,2)=1.
    dtul(:,3)=0.

  elsewhere
    wtul(1:ifull,0)=0.
    wtul(1:ifull,1)=0.
    wtul(1:ifull,2)=0.
    wtul(1:ifull,3)=0.
    dtul(:,1)=0.
    dtul(:,2)=0.
    dtul(:,3)=0.
      
  end where
  where (evnstest)
    wtvl(1:ifull,0)=1.
    wtvl(1:ifull,1)=-0.1
    wtvl(1:ifull,2)=-0.5
    wtvl(1:ifull,3)=0.
    dtvl(:,1)=0.1
    dtvl(:,2)=1.
    dtvl(:,3)=0.5
    !vd(1:ifull,k)=vin(issv,k)*0.1+vin(isv,k)+vin(1:ifull,k)*0.5
    !vin(1:ifull,k)=vd(:,k)-va(inv,k)*0.1-va(isv,k)*0.5
  elsewhere (evntest)
    wtvl(1:ifull,0)=1.
    wtvl(1:ifull,1)=-0.1
    wtvl(1:ifull,2)=-0.5
    wtvl(1:ifull,3)=0.
    dtvl(:,1)=0.
    dtvl(:,2)=1.
    dtvl(:,3)=0.5
  elsewhere (evstest)
    wtvl(1:ifull,0)=1.
    wtvl(1:ifull,1)=0.
    wtvl(1:ifull,2)=-1./3.
    wtvl(1:ifull,3)=0.
    dtvl(:,1)=0.
    dtvl(:,2)=1.
    dtvl(:,3)=0.
  elsewhere (enntest)
    wtvl(1:ifull,0)=1.
    wtvl(1:ifull,1)=-1./3.
    wtvl(1:ifull,2)=0.
    wtvl(1:ifull,3)=0.
    dtvl(:,1)=0.
    dtvl(:,2)=0.
    dtvl(:,3)=1.
  elsewhere (entest)
    wtvl(1:ifull,0)=1.
    wtvl(1:ifull,1)=0.
    wtvl(1:ifull,2)=0.
    wtvl(1:ifull,3)=0.
    dtvl(:,1)=0.
    dtvl(:,2)=0.
    dtvl(:,3)=1.
  elsewhere (evtest)
    wtvl(1:ifull,0)=1.
    wtvl(1:ifull,1)=0.
    wtvl(1:ifull,2)=0.
    wtvl(1:ifull,3)=0.
    dtvl(:,1)=0.
    dtvl(:,2)=1.
    dtvl(:,3)=0.
  elsewhere
    wtvl(1:ifull,0)=0.
    wtvl(1:ifull,1)=0.
    wtvl(1:ifull,2)=0.
    wtvl(1:ifull,3)=0.
    dtvl(:,1)=0.
    dtvl(:,2)=0.
    dtvl(:,3)=0.
  end where
    
  ! Apply JLM's preconditioner
  call boundsuv(wtul,wtvl,stag=-9)
  where (abs(wtul(iwu,0))>1.E-4.and.abs(wtul(1:ifull,2))>1.E-4)
    stul=-wtul(1:ifull,2)/wtul(iwu,0)
    nwtu(:,0)=wtul(1:ifull,0)+stul*wtul(iwu,1)
    nwtu(:,1)=wtul(1:ifull,1)
    nwtu(:,2)=wtul(1:ifull,2)+stul*wtul(iwu,0)
    nwtu(:,3)=wtul(1:ifull,3)-stul*wtul(iwu,2)
  elsewhere
    stul=0.
    nwtu(:,0)=wtul(1:ifull,0)
    nwtu(:,1)=wtul(1:ifull,1)
    nwtu(:,2)=wtul(1:ifull,2)
    nwtu(:,3)=wtul(1:ifull,3)
  end where
  where (abs(wtvl(isv,0))>1.E-4.and.abs(wtvl(1:ifull,2))>1.E-4)
    stvl=-wtvl(1:ifull,2)/wtvl(isv,0)
    nwtv(:,0)=wtvl(1:ifull,0)+stvl*wtvl(isv,1)
    nwtv(:,1)=wtvl(1:ifull,1)
    nwtv(:,2)=wtvl(1:ifull,2)+stvl*wtvl(isv,0)
    nwtv(:,3)=wtvl(1:ifull,3)-stvl*wtvl(isv,2)
  elsewhere
    stvl=0.
    nwtv(:,0)=wtvl(1:ifull,0)
    nwtv(:,1)=wtvl(1:ifull,1)
    nwtv(:,2)=wtvl(1:ifull,2)
    nwtv(:,3)=wtvl(1:ifull,3)
  end where
  where (abs(wtul(1:ifull,0))<1.E-4)
    wtul(1:ifull,0)=1.
    wtul(1:ifull,1)=0.
    wtul(1:ifull,2)=0.
    wtul(1:ifull,3)=0.
  elsewhere
    wtul(1:ifull,0)=nwtu(1:ifull,0)
    wtul(1:ifull,1)=nwtu(1:ifull,1)
    wtul(1:ifull,2)=nwtu(1:ifull,2)
    wtul(1:ifull,3)=nwtu(1:ifull,3)
  end where
  where (abs(wtvl(1:ifull,0))<1.E-4)
    wtvl(1:ifull,0)=1.
    wtvl(1:ifull,1)=0.
    wtvl(1:ifull,2)=0.
    wtvl(1:ifull,3)=0.
  elsewhere
    wtvl(1:ifull,0)=nwtv(1:ifull,0)
    wtvl(1:ifull,1)=nwtv(1:ifull,1)
    wtvl(1:ifull,2)=nwtv(1:ifull,2)
    wtvl(1:ifull,3)=nwtv(1:ifull,3)
  end where

  ! normalise
  call boundsuv(wtul,wtvl,stag=-9)
  do k=1,3
    wtul(1:ifull,k)=wtul(1:ifull,k)/wtul(1:ifull,0)
    wtvl(1:ifull,k)=wtvl(1:ifull,k)/wtvl(1:ifull,0)
    dtul(:,k)=dtul(:,k)/wtul(1:ifull,0)
    dtvl(:,k)=dtvl(:,k)/wtvl(1:ifull,0)
  end do
  stul=stul*wtul(iwu,0)/wtul(1:ifull,0)
  stvl=stvl*wtvl(isv,0)/wtvl(1:ifull,0)
  wtul(1:ifull,0)=1.
  wtvl(1:ifull,0)=1.

  ! assign land arrays
  eutest=eeu(1:ifull)>0.5
  evtest=eev(1:ifull)>0.5
  euetest=eutest.and.eeu(ieu)>0.5
  evntest=evtest.and.eev(inv)>0.5
  euwtest=eutest.and.eeu(iwu)>0.5
  evstest=evtest.and.eev(isv)>0.5
  euewtest=euetest.and.euwtest
  evnstest=evntest.and.evstest
  ewwtest=ewtest.and.ee(iww)>0.5
  esstest=estest.and.ee(iss)>0.5

  ! assign weights (right) 
  
!  |   W   |   * X |  E   |     unstaggered
!          W     X *      E     staggered
  where (euewtest)
    wtur(1:ifull,0)=1.
    wtur(1:ifull,1)=-0.5
    wtur(1:ifull,2)=-0.1
    wtur(1:ifull,3)=0.
    dtur(:,1)=0.1
    dtur(:,2)=1.
    dtur(:,3)=0.5
    !ud(1:ifull,k)=uin(ieu,k)*0.1+uin(1:ifull,k)+uin(iwu,k)*0.5
    !uin(1:ifull,k)=ud(:,k)-ua(iwu,k)*0.1-ua(ieu,k)*0.5
        
!  |   W   |   * X |  E   #     unstaggered
!          W     X *      0     staggered
  elsewhere (euwtest)
    wtur(1:ifull,0)=1.
    wtur(1:ifull,1)=-0.5
    wtur(1:ifull,2)=-0.1
    wtur(1:ifull,3)=0.
    dtur(:,1)=0.
    dtur(:,2)=1.
    dtur(:,3)=0.5
      
!  #   ##  #   * X |   E   |     unstaggered
!  #   ##  0     X *             staggered
  elsewhere (euetest)
    wtur(1:ifull,0)=1.
    wtur(1:ifull,1)=-1./3.
    wtur(1:ifull,2)=0.
    wtur(1:ifull,3)=0.
    dtur(:,1)=0.
    dtur(:,2)=1.
    dtur(:,3)=0.

!  |   W   | X *   #  ##  #     unstaggered
!          W X     0  ##  #     staggered
  elsewhere (ewwtest)
    wtur(1:ifull,0)=1.
    wtur(1:ifull,1)=0.
    wtur(1:ifull,2)=-1./3.
    wtur(1:ifull,3)=0.
    dtur(:,1)=0.
    dtur(:,2)=0.
    dtur(:,3)=1.

!  #   W   |   *   #  ##  #     unstaggered
!  #       W       #  ##  #     staggered
  elsewhere (ewtest)
    wtur(1:ifull,0)=1.
    wtur(1:ifull,1)=0.
    wtur(1:ifull,2)=0.
    wtur(1:ifull,3)=0.
    dtur(:,1)=0.
    dtur(:,2)=0.
    dtur(:,3)=1.

!  #   ##  #   *   |  E   #     unstaggered
!  #   ##  #       *      #     staggered
  elsewhere (eutest)
    wtur(1:ifull,0)=1.
    wtur(1:ifull,1)=0.
    wtur(1:ifull,2)=0.
    wtur(1:ifull,3)=0.
    dtur(:,1)=0.
    dtur(:,2)=1.
    dtur(:,3)=0.
  elsewhere
    wtur(1:ifull,0)=0.
    wtur(1:ifull,1)=0.
    wtur(1:ifull,2)=0.
    wtur(1:ifull,3)=0.
    dtur(:,1)=0.
    dtur(:,2)=0.
    dtur(:,3)=0.

  end where
  where (evnstest)
    wtvr(1:ifull,0)=1.
    wtvr(1:ifull,1)=-0.5
    wtvr(1:ifull,2)=-0.1
    wtvr(1:ifull,3)=0.
    dtvr(:,1)=0.1
    dtvr(:,2)=1.
    dtvr(:,3)=0.5
    !vd(1:ifull,k)=vin(inv,k)*0.1+vin(1:ifull,k)+vin(isv,k)*0.5
    !vin(1:ifull,k)=vd(:,k)-va(isv,k)*0.1-va(inv,k)*0.5
  elsewhere (evstest)
    wtvr(1:ifull,0)=1.
    wtvr(1:ifull,1)=-0.5
    wtvr(1:ifull,2)=-0.1
    wtvr(1:ifull,3)=0.
    dtvr(:,1)=0.
    dtvr(:,2)=1.
    dtvr(:,3)=0.5
  elsewhere (evntest)
    wtvr(1:ifull,0)=1.
    wtvr(1:ifull,1)=-1./3.
    wtvr(1:ifull,2)=0.
    wtvr(1:ifull,3)=0.
    dtvr(:,1)=0.
    dtvr(:,2)=1.
    dtvr(:,3)=0.
  elsewhere (esstest)
    wtvr(1:ifull,0)=1.
    wtvr(1:ifull,1)=0.
    wtvr(1:ifull,2)=-1./3.
    wtvr(1:ifull,3)=0.
    dtvr(:,1)=0.
    dtvr(:,2)=0.
    dtvr(:,3)=1.
  elsewhere (estest)
    wtvr(1:ifull,0)=1.
    wtvr(1:ifull,1)=0.
    wtvr(1:ifull,2)=0.
    wtvr(1:ifull,3)=0.
    dtvr(:,1)=0.
    dtvr(:,2)=0.
    dtvr(:,3)=1.
  elsewhere (evtest)
    wtvr(1:ifull,0)=1.
    wtvr(1:ifull,1)=0.
    wtvr(1:ifull,2)=0.
    wtvr(1:ifull,3)=0.
    dtvr(:,1)=0.
    dtvr(:,2)=1.
    dtvr(:,3)=0.
  elsewhere
    wtvr(1:ifull,0)=0.
    wtvr(1:ifull,1)=0.
    wtvr(1:ifull,2)=0.
    wtvr(1:ifull,3)=0.
    dtvr(:,1)=0.
    dtvr(:,2)=0.
    dtvr(:,3)=0.
  end where

  ! Apply JLM's preconditioner
  call boundsuv(wtur,wtvr,stag=-10)
  where (abs(wtur(ieu,0))>1.E-4.and.abs(wtur(1:ifull,1))>1.E-4)
    stur=-wtur(1:ifull,1)/wtur(ieu,0)
    nwtu(1:ifull,0)=wtur(1:ifull,0)+stur*wtur(ieu,2)
    nwtu(1:ifull,1)=wtur(1:ifull,1)+stur*wtur(ieu,0)
    nwtu(1:ifull,2)=wtur(1:ifull,2)
    nwtu(1:ifull,3)=wtur(1:ifull,3)-stur*wtur(ieu,1)
  elsewhere
    stur=0.
    nwtu(1:ifull,0)=wtur(1:ifull,0)
    nwtu(1:ifull,1)=wtur(1:ifull,1)
    nwtu(1:ifull,2)=wtur(1:ifull,2)
    nwtu(1:ifull,3)=wtur(1:ifull,3)
  end where
  where (abs(wtvr(inv,0))>1.E-4.and.abs(wtvr(1:ifull,1))>1.E-4)
    stvr=-wtvr(1:ifull,1)/wtvr(inv,0)
    nwtv(1:ifull,0)=wtvr(1:ifull,0)+stvr*wtvr(inv,2)
    nwtv(1:ifull,1)=wtvr(1:ifull,1)+stvr*wtvr(inv,0)
    nwtv(1:ifull,2)=wtvr(1:ifull,2)
    nwtv(1:ifull,3)=wtvr(1:ifull,3)-stvr*wtvr(inv,1)
  elsewhere
    stvr=0.
    nwtv(1:ifull,0)=wtvr(1:ifull,0)
    nwtv(1:ifull,1)=wtvr(1:ifull,1)
    nwtv(1:ifull,2)=wtvr(1:ifull,2)
    nwtv(1:ifull,3)=wtvr(1:ifull,3)
  end where
  where (abs(wtur(1:ifull,0))<1.E-4)
    wtur(1:ifull,0)=1.
    wtur(1:ifull,1)=0.
    wtur(1:ifull,2)=0.
    wtur(1:ifull,3)=0.
  elsewhere
    wtur(1:ifull,0)=nwtu(1:ifull,0)
    wtur(1:ifull,1)=nwtu(1:ifull,1)
    wtur(1:ifull,2)=nwtu(1:ifull,2)
    wtur(1:ifull,3)=nwtu(1:ifull,3)
  end where
  where (abs(wtvr(1:ifull,0))<1.E-4)
    wtvr(1:ifull,0)=1.
    wtvr(1:ifull,1)=0.
    wtvr(1:ifull,2)=0.
    wtvr(1:ifull,3)=0.
  elsewhere
    wtvr(1:ifull,0)=nwtv(1:ifull,0)
    wtvr(1:ifull,1)=nwtv(1:ifull,1)
    wtvr(1:ifull,2)=nwtv(1:ifull,2)
    wtvr(1:ifull,3)=nwtv(1:ifull,3)
  end where

  ! normalise
  call boundsuv(wtur,wtvr,stag=-10)
  do k=1,3
    wtur(1:ifull,k)=wtur(1:ifull,k)/wtur(1:ifull,0)
    wtvr(1:ifull,k)=wtvr(1:ifull,k)/wtvr(1:ifull,0)
    dtur(:,k)=dtur(:,k)/wtur(1:ifull,0)
    dtvr(:,k)=dtvr(:,k)/wtvr(1:ifull,0)
  end do
  stur=stur*wtur(ieu,0)/wtur(1:ifull,0)
  stvr=stvr*wtvr(inv,0)/wtvr(1:ifull,0)
  wtur(1:ifull,0)=1.
  wtvr(1:ifull,0)=1.
  
end if

kx=size(u,2)

zoff=0
if (present(toff)) then
  if (toff==1) zoff=koff
end if

do k=1,kx
  uin(1:ifull,k)=u(:,k)*eeu(1:ifull)
  vin(1:ifull,k)=v(:,k)*eev(1:ifull)
end do

if (mstagf==0) then
  ltest=.true.
else if (mstagf<0) then
  ltest=.false.
else
  ltest=mod(ktau-zoff-nstagoffmlo,2*mstagf)<mstagf
end if

if (ltest) then
  
  call boundsuv(uin,vin,stag=5)
  do k=1,kx
    ud(1:ifull,k)=dtul(:,1)*uin(iwwu,k)+dtul(:,2)*uin(iwu,k)+dtul(:,3)*uin(1:ifull,k)
    vd(1:ifull,k)=dtvl(:,1)*vin(issv,k)+dtvl(:,2)*vin(isv,k)+dtvl(:,3)*vin(1:ifull,k)
  end do
  
  call boundsuv(ud,vd,stag=-9)
  do k=1,kx
    ! Apply JLM's preconditioner
    nud=ud(1:ifull,k)-stul*ud(iwu,k)
    nvd=vd(1:ifull,k)-stvl*vd(isv,k)
    ud(1:ifull,k)=nud
    vd(1:ifull,k)=nvd

    ! 1st guess
    !ua(1:ifull,k)=(uin(1:ifull,k)+uin(iwu,k))*0.5*ee(1:ifull)
    !va(1:ifull,k)=(vin(1:ifull,k)+vin(isv,k))*0.5*ee(1:ifull)
    ua(1:ifull,k)=nud
    va(1:ifull,k)=nvd
  end do

  do itn=1,itnmax        ! each loop is a double iteration
    call boundsuv(ua,va,stag=3)
    do k=1,kx
      uin(1:ifull,k)=ud(1:ifull,k)+wtul(:,1)*ua(ieu,k)+wtul(:,2)*ua(iwu,k)+wtul(:,3)*ua(iwwu,k)
      vin(1:ifull,k)=vd(1:ifull,k)+wtvl(:,1)*va(inv,k)+wtvl(:,2)*va(isv,k)+wtvl(:,3)*va(issv,k)
    end do
    call boundsuv(uin,vin,stag=3)
    do k=1,kx
      ua(1:ifull,k)=ud(1:ifull,k)+wtul(:,1)*uin(ieu,k)+wtul(:,2)*uin(iwu,k)+wtul(:,3)*uin(iwwu,k)
      va(1:ifull,k)=vd(1:ifull,k)+wtvl(:,1)*vin(inv,k)+wtvl(:,2)*vin(isv,k)+wtvl(:,3)*vin(issv,k)
    end do

  end do                  ! itn=1,itnmax
  
else

  call boundsuv(uin,vin)
  do k=1,kx
    ud(1:ifull,k)=dtur(:,1)*uin(ieu,k)+dtur(:,2)*uin(1:ifull,k)+dtur(:,3)*uin(iwu,k)
    vd(1:ifull,k)=dtvr(:,1)*vin(inv,k)+dtvr(:,2)*vin(1:ifull,k)+dtvr(:,3)*vin(isv,k)
  end do

  call boundsuv(ud,vd,stag=-10)
  do k=1,kx
    ! Apply JLM's preconditioner
    nud=ud(1:ifull,k)-stur*ud(ieu,k)
    nvd=vd(1:ifull,k)-stvr*vd(inv,k)
    ud(1:ifull,k)=nud
    vd(1:ifull,k)=nvd

    ! 1st guess
    !ua(1:ifull,k)=(uin(1:ifull,k)+uin(iwu,k))*0.5*ee(1:ifull)
    !va(1:ifull,k)=(vin(1:ifull,k)+vin(isv,k))*0.5*ee(1:ifull)
    ua(1:ifull,k)=nud
    va(1:ifull,k)=nvd
  end do

  do itn=1,itnmax        ! each loop is a double iteration
    call boundsuv(ua,va,stag=2)
    do k=1,kx
      uin(1:ifull,k)=ud(1:ifull,k)+wtur(:,1)*ua(ieu,k)+wtur(:,2)*ua(iwu,k)+wtur(:,3)*ua(ieeu,k)
      vin(1:ifull,k)=vd(1:ifull,k)+wtvr(:,1)*va(inv,k)+wtvr(:,2)*va(isv,k)+wtvr(:,3)*va(innv,k)
    end do
    call boundsuv(uin,vin,stag=2)
    do k=1,kx
      ua(1:ifull,k)=ud(1:ifull,k)+wtur(:,1)*uin(ieu,k)+wtur(:,2)*uin(iwu,k)+wtur(:,3)*uin(ieeu,k)
      va(1:ifull,k)=vd(1:ifull,k)+wtvr(:,1)*vin(inv,k)+wtvr(:,2)*vin(isv,k)+wtvr(:,3)*vin(innv,k)
    end do
  end do                  ! itn=1,itnmax

end if

uout=ua(1:ifull,:)
vout=va(1:ifull,:)

call end_log(ocnstag_end)

return
end subroutine mlounstaguv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine performs vertical advection based on JLMs TVD scheme

subroutine mlovadv(dtin,ww,uu,vv,ss,tt,mm,depdum,idzdum,wtr,cnum)

use cc_mpi
use mlo

implicit none

include 'newmpar.h'

integer, intent(in) :: cnum
integer ii,l,iq,ierr,its_g
integer, dimension(ifull) :: its
integer, dimension(ifull,wlev-1) :: kp,kx
real, intent(in) :: dtin
real, dimension(ifull) :: dtnew
real, dimension(ifull,0:wlev), intent(in) :: ww
real, dimension(ifull,wlev), intent(in) :: depdum,idzdum
real, dimension(ifull,wlev), intent(inout) :: uu,vv,ss,tt,mm
real, dimension(ifull,wlev) :: dzdum
logical, dimension(ifull), intent(in) :: wtr

dzdum=max(idzdum,1.E-10)

! reduce time step to ensure stability
dtnew(:)=dtin
do iq=1,ifull
  if (wtr(iq)) then
    do ii=1,wlev-1
      ! this trick works if dzdum(iq,ii)<dzdum(iq,ii+1)
      dtnew(iq)=min(dtnew(iq),0.3*dzdum(iq,ii)/max(abs(ww(iq,ii)),1.E-12))
    end do
  end if
end do
its(:)=int(dtin/(dtnew(:)+0.01))+1
its_g=maxval(its(:))
if (its_g>500) then
  write(6,*) "MLOVERT myid,cnum,its_g",myid,cnum,its_g
end if
dtnew(:)=dtin/real(its(:))

do ii=1,wlev-1
  kp(:,ii)=sign(1.,ww(:,ii))
  kx(:,ii)=ii+(1-kp(:,ii))/2 !  k for ww +ve,  k+1 for ww -ve
end do
  
tt=tt-290.
ss=ss-34.72
call mlotvd(its,dtnew,ww,uu,depdum,dzdum,kp,kx)
call mlotvd(its,dtnew,ww,vv,depdum,dzdum,kp,kx)
call mlotvd(its,dtnew,ww,ss,depdum,dzdum,kp,kx)
call mlotvd(its,dtnew,ww,tt,depdum,dzdum,kp,kx)
call mlotvd(its,dtnew,ww,mm,depdum,dzdum,kp,kx)
tt=tt+290.
ss=ss+34.72
ss=max(ss,0.)

return
end subroutine mlovadv

subroutine mlotvd(its,dtnew,ww,uu,depadj,dzadj,kp,kx)

use mlo

implicit none

include 'newmpar.h'

integer ii,i,iq
integer, dimension(ifull), intent(in) :: its
integer, dimension(ifull,wlev-1), intent(in) :: kp,kx
real, dimension(ifull), intent(in) :: dtnew
real, dimension(ifull,0:wlev), intent(in) :: ww
real, dimension(ifull,wlev), intent(in) :: depadj,dzadj
real, dimension(ifull,wlev), intent(inout) :: uu
real, dimension(ifull,wlev-1) :: ff
real, dimension(ifull,0:wlev) :: delu
real, dimension(ifull) :: fl,fh,cc,rr,xx

! f=(w*u) at half levels
! du/dt = u*dw/dz-df/dz = -w*du/dz

delu(:,0)=0.
do ii=1,wlev-1
  delu(:,ii)=uu(:,ii+1)-uu(:,ii)
end do
delu(:,wlev)=0.

! TVD part
do ii=1,wlev-1
  xx=delu(:,ii)+sign(1.E-20,delu(:,ii))
  do iq=1,ifull
    rr(iq)=delu(iq,ii-kp(iq,ii))/xx(iq)
    fl(iq)=ww(iq,ii)*uu(iq,kx(iq,ii))
  end do
  fh=ww(:,ii)*0.5*(uu(:,ii)+uu(:,ii+1)) &
     -0.5*(uu(:,ii+1)-uu(:,ii))*ww(:,ii)**2*dtnew/max(depadj(:,ii+1)-depadj(:,ii),1.E-10)
  cc=max(0.,min(1.,2.*rr),min(2.,rr)) ! superbee
  ff(:,ii)=fl+cc*(fh-fl)
  !ff(:,ii)=ww(:,ii)*0.5*(uu(:,ii)+uu(:,ii+1)) ! explicit
end do
uu(:,1)=uu(:,1)+dtnew*(uu(:,1)*ww(:,1)-ff(:,1))/dzadj(:,1)
do ii=2,wlev-1
  uu(:,ii)=uu(:,ii)+dtnew*(uu(:,ii)*(ww(:,ii)-ww(:,ii-1))-ff(:,ii)+ff(:,ii-1))/dzadj(:,ii)
end do
uu(:,wlev)=uu(:,wlev)+dtnew*(-uu(:,wlev)*ww(:,wlev-1)+ff(:,wlev-1))/dzadj(:,wlev)


do iq=1,ifull
  do i=2,its(iq)

    do ii=1,wlev-1
      delu(iq,ii)=uu(iq,ii+1)-uu(iq,ii)
    end do

    ! TVD part
    do ii=1,wlev-1
      rr(iq)=delu(iq,ii-kp(iq,ii))
      fl(iq)=ww(iq,ii)*uu(iq,kx(iq,ii))
      fh(iq)=ww(iq,ii)*0.5*(uu(iq,ii)+uu(iq,ii+1)) &
        -0.5*(uu(iq,ii+1)-uu(iq,ii))*ww(iq,ii)**2*dtnew(iq)/max(depadj(iq,ii+1)-depadj(iq,ii),1.E-10)
      xx(iq)=delu(iq,ii)+sign(1.E-20,delu(iq,ii))
      cc(iq)=max(0.,min(1.,2.*rr(iq)),min(2.,rr(iq))) ! superbee
      ff(iq,ii)=fl(iq)+cc(iq)*(fh(iq)-fl(iq))
    end do
  
    uu(iq,1)=uu(iq,1)+dtnew(iq)*(uu(iq,1)*ww(iq,1)-ff(iq,1))/dzadj(iq,1)
    do ii=2,wlev-1
      uu(iq,ii)=uu(iq,ii)+dtnew(iq)*(uu(iq,ii)*(ww(iq,ii)-ww(iq,ii-1))-ff(iq,ii)+ff(iq,ii-1))/dzadj(iq,ii)
    end do
    uu(iq,wlev)=uu(iq,wlev)+dtnew(iq)*(-uu(iq,wlev)*ww(iq,wlev-1)+ff(iq,wlev-1))/dzadj(iq,wlev)

  end do
end do

return
end subroutine mlotvd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Use potential temperature and salinity Jacobians to calculate
! density Jacobian

subroutine tsjacobi(nti,nsi,alphabar,betabar,neta,drhobardxu,drhobardyu,drhobardxv,drhobardyv)

use indices_m
use mlo, only : wlev

implicit none

include 'newmpar.h'

integer ii
real, dimension(ifull+iextra,wlev), intent(in) :: nti,nsi,alphabar,betabar
real, dimension(ifull+iextra), intent(in) :: neta
real, dimension(ifull,wlev), intent(out) :: drhobardxu,drhobardyu,drhobardxv,drhobardyv
real, dimension(ifull+iextra,wlev) :: nt,ns
real, dimension(ifull,wlev) :: absu,bbsu,absv,bbsv
real, dimension(ifull,wlev) :: dntdxu,dntdxv,dntdyu,dntdyv
real, dimension(ifull,wlev) :: dnsdxu,dnsdxv,dnsdyu,dnsdyv
integer, parameter :: rhogradmeth = 1 ! Density gradient method (0 = finite difference, 1 = local interpolation)

nt=min(max(271.,nti),373.)-290.
ns=min(max(0.,nsi),50.)-34.72

do ii=1,wlev
  absu(:,ii)=0.5*(alphabar(1:ifull,ii)+alphabar(ie,ii))*eeu(1:ifull)
  bbsu(:,ii)=0.5*(betabar(1:ifull,ii)+betabar(ie,ii))*eeu(1:ifull)
  absv(:,ii)=0.5*(alphabar(1:ifull,ii)+alphabar(in,ii))*eev(1:ifull)
  bbsv(:,ii)=0.5*(betabar(1:ifull,ii)+betabar(in,ii))*eev(1:ifull)
end do

select case(rhogradmeth)
  case(0) ! finite difference
    call finitedelta(nt,neta,dntdxu,dntdyu,dntdxv,dntdyv)
    call finitedelta(ns,neta,dnsdxu,dnsdyu,dnsdxv,dnsdyv)
  case(1) ! local interpolation
    call seekdelta(nt,neta,dntdxu,dntdyu,dntdxv,dntdyv)
    call seekdelta(ns,neta,dnsdxu,dnsdyu,dnsdxv,dnsdyv)
  case default
    write(6,*) "ERROR: Invalid choice of rhogradmeth ",rhogradmeth
    stop
end select

drhobardxu=-absu*dntdxu+bbsu*dnsdxu
drhobardxv=-absv*dntdxv+bbsv*dnsdxv
drhobardyu=-absu*dntdyu+bbsu*dnsdyu
drhobardyv=-absv*dntdyv+bbsv*dnsdyv

return
end subroutine tsjacobi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate gradients using a finite difference method

subroutine finitedelta(rhobar,neta,drhobardxu,drhobardyu,drhobardxv,drhobardyv)

use indices_m
use map_m
use mlo, only : wlev,minwater

implicit none

include 'newmpar.h'
include 'parm.h'

integer ii
real, dimension(ifull+iextra,wlev), intent (in) :: rhobar
real, dimension(ifull+iextra), intent(in) :: neta
real, dimension(ifull,wlev), intent(out) :: drhobardxu,drhobardyu,drhobardxv,drhobardyv
real, dimension(ifull) :: drhobardzu,drhobardzv,tnu,tsu,tev,twv,tee,tnn
real, dimension(ifull+iextra,0:wlev) :: rhobarhl,dephladj
real xp

! set-up level depths
dephladj(:,0)=0.
do ii=1,wlev
  dephladj(:,ii)=gosigh(ii)*max(dd(:)+neta(:),minwater)
end do

! estimate rhobar at half levels
rhobarhl(:,0)=rhobar(:,1)
do ii=1,wlev-1
  xp=gosigh(ii)-gosig(ii)
  rhobarhl(:,ii)=rhobar(:,ii)+xp*(rhobar(:,ii+1)-rhobar(:,ii))/(gosig(ii+1)-gosig(ii))
end do
rhobarhl(:,wlev)=rhobar(:,wlev)

do ii=1,wlev
  drhobardzu=(rhobarhl(1:ifull,ii)+rhobarhl(ie,ii)-rhobarhl(1:ifull,ii-1)-rhobarhl(ie,ii-1))          &
            /(dephladj(1:ifull,ii)+dephladj(ie,ii)-dephladj(1:ifull,ii-1)-dephladj(ie,ii-1))
  drhobardzv=(rhobarhl(1:ifull,ii)+rhobarhl(in,ii)-rhobarhl(1:ifull,ii-1)-rhobarhl(in,ii-1))          &
            /(dephladj(1:ifull,ii)+dephladj(in,ii)-dephladj(1:ifull,ii-1)-dephladj(in,ii-1))
  tnu=0.5*(( rhobar(in,ii)+drhobardzu*((1.-gosig(ii))*neta(in) -gosig(ii)*dd(in) ))*f(in) &
          +(rhobar(ine,ii)+drhobardzu*((1.-gosig(ii))*neta(ine)-gosig(ii)*dd(ine)))*f(ine))
  tee=0.5*(( rhobar(1:ifull,ii)+drhobardzu*((1.-gosig(ii))*neta(1:ifull) -gosig(ii)*dd(1:ifull) ))*f(1:ifull) &
          +(rhobar(ie,ii)+drhobardzu*((1.-gosig(ii))*neta(ie)-gosig(ii)*dd(ie)))*f(ie))
  tsu=0.5*(( rhobar(is,ii)+drhobardzu*((1.-gosig(ii))*neta(is) -gosig(ii)*dd(is) ))*f(is) &
          +(rhobar(ise,ii)+drhobardzu*((1.-gosig(ii))*neta(ise)-gosig(ii)*dd(ise)))*f(ise))
  tev=0.5*(( rhobar(ie,ii)+drhobardzv*((1.-gosig(ii))*neta(ie) -gosig(ii)*dd(ie) ))*f(ie) &
          +(rhobar(ien,ii)+drhobardzv*((1.-gosig(ii))*neta(ien)-gosig(ii)*dd(ien)))*f(ien))
  tnn=0.5*(( rhobar(1:ifull,ii)+drhobardzv*((1.-gosig(ii))*neta(1:ifull) -gosig(ii)*dd(1:ifull) ))*f(1:ifull) &
          +(rhobar(in,ii)+drhobardzv*((1.-gosig(ii))*neta(in)-gosig(ii)*dd(in)))*f(in))
  twv=0.5*(( rhobar(iw,ii)+drhobardzv*((1.-gosig(ii))*neta(iw) -gosig(ii)*dd(iw) ))*f(iw) &
          +(rhobar(iwn,ii)+drhobardzv*((1.-gosig(ii))*neta(iwn)-gosig(ii)*dd(iwn)))*f(iwn))
  drhobardxu(:,ii)=(rhobar(ie,ii)-rhobar(1:ifull,ii)+drhobardzu*((1.-gosig(ii))*(neta(ie)-neta(1:ifull)) &
                                                     -gosig(ii)*(dd(ie)-dd(1:ifull))))*emu(1:ifull)/ds
  drhobardxu(:,ii)=drhobardxu(:,ii)*eeu(1:ifull)                                                   
  drhobardyu(:,ii)=0.5*(stwgt(:,1)*(tnu-tee)+stwgt(:,2)*(tee-tsu))*emu(1:ifull)/ds-0.5*(rhobar(1:ifull,ii)+rhobar(ie,ii))*dfdyu
  drhobardxv(:,ii)=0.5*(stwgt(:,3)*(tev-tnn)+stwgt(:,4)*(tnn-twv))*emv(1:ifull)/ds-0.5*(rhobar(1:ifull,ii)+rhobar(in,ii))*dfdxv
  drhobardyv(:,ii)=(rhobar(in,ii)-rhobar(1:ifull,ii)+drhobardzv*((1.-gosig(ii))*(neta(in)-neta(1:ifull)) &
                                                     -gosig(ii)*(dd(in)-dd(1:ifull))))*emv(1:ifull)/ds
  drhobardyv(:,ii)=drhobardyv(:,ii)*eev(1:ifull)
end do

return
end subroutine finitedelta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate gradients using an interpolation method

subroutine seekdelta(rhobar,neta,drhobardxu,drhobardyu,drhobardxv,drhobardyv)

use indices_m
use map_m
use mlo, only : wlev,minwater

implicit none

include 'newmpar.h'
include 'parm.h'

integer iqq,ii
integer sdi,sde,sdw,sdn,sds,sden,sdse,sdne,sdwn
real ddux,ddvy,ri,re,rw,rn,rs,ren,rse,rne,rwn
real mxi,mxe,mxw,mxn,mxs,mxen,mxse,mxne,mxwn
real drhobardzu,drhobardzv
real, dimension(wlev) :: ddi,dde,ddw,ddn,dds,dden,ddse,ddne,ddwn
real, dimension(wlev) :: ssi,sse,ssw,ssn,sss,ssen,ssse,ssne,sswn
real, dimension(ifull+iextra,wlev), intent (in) :: rhobar
real, dimension(ifull+iextra), intent(in) :: neta
real, dimension(ifull,wlev), intent(out) :: drhobardxu,drhobardyu,drhobardxv,drhobardyv
real, dimension(ifull+iextra,0:wlev) :: rhobarhl,dephladj
real xp

! set-up level depths
dephladj(:,0)=0.
do ii=1,wlev
  dephladj(:,ii)=gosigh(ii)*max(dd(:)+neta(:),minwater)
end do

! estimate rhobar at half levels
rhobarhl(:,0)=rhobar(:,1)
do ii=1,wlev-1
  xp=gosigh(ii)-gosig(ii)
  rhobarhl(:,ii)=rhobar(:,ii)+xp*(rhobar(:,ii+1)-rhobar(:,ii))/(gosig(ii+1)-gosig(ii))
end do
rhobarhl(:,wlev)=rhobar(:,wlev)

drhobardxu=0.
drhobardyu=0.
drhobardyv=0.
drhobardxv=0.

do iqq=1,ifull
  mxi=dd(iqq)
  mxe=dd(ie(iqq))
  mxn=dd(in(iqq))
  ssi=rhobar(iqq,:)
  sse=rhobar(ie(iqq),:)
  ssn=rhobar(in(iqq),:)
  ddi=gosig(:)*mxi
  dde=gosig(:)*mxe
  ddn=gosig(:)*mxn

  if (eeu(iqq)>0.5) then ! water
    sdi=2
    sde=2
    sdn=2
    sds=2
    sdne=2
    sdse=2
    mxs=dd(is(iqq))
    mxne=dd(ine(iqq))
    mxse=dd(ise(iqq))
    sss=rhobar(is(iqq),:)
    ssne=rhobar(ine(iqq),:)
    ssse=rhobar(ise(iqq),:)
    dds=gosig(:)*mxs
    ddne=gosig(:)*mxne
    ddse=gosig(:)*mxse
    do ii=1,wlev
      ! process staggered u locations
      ! use scaled depths (i.e., assume neta is small compared to dd)
      ddux=gosig(ii)*ddu(iqq) ! seek depth
      if (mxi>=ddux.and.mxe>=ddux) then
        ! estimate vertical gradient
        drhobardzu=(rhobarhl(iqq,ii)+rhobarhl(ie(iqq),ii)-rhobarhl(iqq,ii-1)-rhobarhl(ie(iqq),ii-1))          &
                  /(dephladj(iqq,ii)+dephladj(ie(iqq),ii)-dephladj(iqq,ii-1)-dephladj(ie(iqq),ii-1))
        call seekval(ri,ssi(:),ddi(:),ddux,sdi)
        call seekval(re,sse(:),dde(:),ddux,sde)
        ! the following terms correct for neglecting neta in the above intepolation of depths
        ri=ri+drhobardzu*(1.-gosig(ii))*neta(iqq)
        re=re+drhobardzu*(1.-gosig(ii))*neta(ie(iqq))
        drhobardxu(iqq,ii)=eeu(iqq)*(re-ri)*emu(iqq)/ds
        if (mxn>=ddux.and.mxne>=ddux) then
          call seekval(rn, ssn(:), ddn(:), ddux,sdn)
          call seekval(rne,ssne(:),ddne(:),ddux,sdne)
          ! the following terms correct for neglecting neta in the above interpolation of depths
          rn=rn+drhobardzu*(1.-gosig(ii))*neta(in(iqq))
          rne=rne+drhobardzu*(1.-gosig(ii))*neta(ine(iqq))
          drhobardyu(iqq,ii)=drhobardyu(iqq,ii)+0.25*stwgt(iqq,1)*(rn*f(in(iqq))+rne*f(ine(iqq))-ri*f(iqq)-re*f(ie(iqq)))*emu(iqq)/ds
          drhobardyu(iqq,ii)=drhobardyu(iqq,ii)-0.125*stwgt(iqq,1)*(ri+re)*(f(in(iqq))+f(ine(iqq))-f(iqq)-f(ie(iqq)))*emu(iqq)/ds
        end if
        if (mxs>=ddux.and.mxse>=ddux) then
          call seekval(rs, sss(:), dds(:), ddux,sds)
          call seekval(rse,ssse(:),ddse(:),ddux,sdse)
          ! the following terms correct for neglecting neta in the above interpolation of depths
          rs=rs+drhobardzu*(1.-gosig(ii))*neta(is(iqq))
          rse=rse+drhobardzu*(1.-gosig(ii))*neta(ise(iqq))
          drhobardyu(iqq,ii)=drhobardyu(iqq,ii)+0.25*stwgt(iqq,2)*(ri*f(iqq)+re*f(ie(iqq))-rs*f(is(iqq))-rse*f(ise(iqq)))*emu(iqq)/ds
          drhobardyu(iqq,ii)=drhobardyu(iqq,ii)-0.125*stwgt(iqq,2)*(ri+re)*(f(iqq)+f(ie(iqq))-f(is(iqq))-f(ise(iqq)))*emu(iqq)/ds
        end if
      end if
    end do
  end if

  if (eev(iqq)>0.5) then ! water
    sdi=2
    sdn=2
    sde=2
    sdw=2
    sden=2
    sdwn=2
    mxw=dd(iw(iqq))
    mxen=dd(ien(iqq))
    mxwn=dd(iwn(iqq))
    ssw=rhobar(iw(iqq),:)
    ssen=rhobar(ien(iqq),:)
    sswn=rhobar(iwn(iqq),:)
    ddw=gosig(:)*mxw
    dden=gosig(:)*mxen
    ddwn=gosig(:)*mxwn
    do ii=1,wlev
      ! now process staggered v locations
      ddvy=gosig(ii)*ddv(iqq) ! seek depth

      if (mxi>=ddvy.and.mxn>=ddvy) then
        ! estimate vertical gradient
        drhobardzv=(rhobarhl(iqq,ii)+rhobarhl(in(iqq),ii)-rhobarhl(iqq,ii-1)-rhobarhl(in(iqq),ii-1))          &
                  /(dephladj(iqq,ii)+dephladj(in(iqq),ii)-dephladj(iqq,ii-1)-dephladj(in(iqq),ii-1))
        call seekval(ri,ssi(:),ddi(:),ddvy,sdi)
        call seekval(rn,ssn(:),ddn(:),ddvy,sdn)
        ri=ri+drhobardzv*(1.-gosig(ii))*neta(iqq)
        rn=rn+drhobardzv*(1.-gosig(ii))*neta(in(iqq))
        drhobardyv(iqq,ii)=eev(iqq)*(rn-ri)*emv(iqq)/ds
        if (mxe>=ddvy.and.mxen>=ddvy) then
          call seekval(re, sse(:), dde(:), ddvy,sde)
          call seekval(ren,ssen(:),dden(:),ddvy,sden)
          re=re+drhobardzv*(1.-gosig(ii))*neta(ie(iqq))
          ren=ren+drhobardzv*(1.-gosig(ii))*neta(ien(iqq))
          drhobardxv(iqq,ii)=drhobardxv(iqq,ii)+0.25*stwgt(iqq,3)*(re*f(ie(iqq))+ren*f(ien(iqq))-ri*f(iqq)-rn*f(in(iqq)))*emv(iqq)/ds
          drhobardxv(iqq,ii)=drhobardxv(iqq,ii)-0.125*stwgt(iqq,3)*(ri+rn)*(f(ie(iqq))+f(ien(iqq))-f(iqq)-f(in(iqq)))*emv(iqq)/ds
        end if
        if (mxw>=ddvy.and.mxwn>=ddvy) then
          call seekval(rw, ssw(:), ddw(:), ddvy,sdw)
          call seekval(rwn,sswn(:),ddwn(:),ddvy,sdwn)
          rw=rw+drhobardzv*(1.-gosig(ii))*neta(iw(iqq))
          rwn=rwn+drhobardzv*(1.-gosig(ii))*neta(iwn(iqq))
          drhobardxv(iqq,ii)=drhobardxv(iqq,ii)+0.25*stwgt(iqq,4)*(ri*f(iqq)+rn*f(in(iqq))-rw*f(iw(iqq))-rwn*f(iwn(iqq)))*emv(iqq)/ds
          drhobardxv(iqq,ii)=drhobardxv(iqq,ii)-0.125*stwgt(iqq,4)*(ri+rn)*(f(iqq)+f(in(iqq))-f(iw(iqq))-f(iwn(iqq)))*emv(iqq)/ds
        end if
      end if
    end do
  end if

end do
 
return
end subroutine seekdelta

subroutine seekval(rout,ssin,ddin,ddseek,sindx)

use mlo, only : wlev

implicit none

integer, intent(inout) :: sindx
integer nindx
real, intent(in) :: ddseek
real, intent(out) :: rout
real, dimension(wlev), intent(in) :: ssin,ddin
real xp

do nindx=sindx,wlev-1
  if (ddin(nindx)>=ddseek) exit
end do
sindx=nindx

xp=(ddseek-ddin(sindx-1))/(ddin(sindx)-ddin(sindx-1))
xp=max(min(xp,1.),0.)

rout=ssin(sindx-1)+xp*(ssin(sindx)-ssin(sindx-1))

return
end subroutine seekval

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate tidal potential

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

! amplitude (mom3)
aa(1)=0.141565 ! K1
aa(2)=0.100661 ! O1
aa(3)=0.046848 ! P1
aa(4)=0.019273 ! Q1
ab(1)=0.242334 ! M2
ab(2)=0.112743 ! S2
ab(3)=0.046397 ! N2
ab(4)=0.030684 ! K2
! Love number (mom3)
ba(1)=0.736
ba(2)=0.695
ba(3)=0.706
ba(4)=0.695
bb(1)=0.693
bb(2)=0.693
bb(3)=0.693
bb(4)=0.693
! Frequency (mom3)
!wa(1)=0.7292117
!wa(2)=0.6750774
!wa(3)=0.7252295
!wa(4)=0.6495854
!wb(1)=1.405189
!wb(2)=1.454441 ! exactly twice per day solar day for S2
!wb(3)=1.378797
!wb(4)=1.458423

stime=real(mtimer)/1440.+real(jstart) ! solar time (days)
ctime=stime/36525. ! century (relative to 12Z 31 Dec 1899)

mn=270.43659+481276.89057*ctime+0.00198*ctime*ctime+0.000002*ctime*ctime*ctime ! moon
sn=279.69668+36000.76892*ctime+0.00030*ctime*ctime                             ! sun
pn=334.32956+4069.03404*ctime-0.01032*ctime*ctime-0.00001*ctime*ctime*ctime    ! lunar perigee
mn=mn*pi/180.
sn=sn*pi/180.
pn=pn*pi/180.
stime=(stime-real(int(stime)))*2.*pi ! solar time (rad)
ltime=stime+sn-mn                    ! lunar time (rad)
! 360*stime=360*ltime-9.3+12.2*day

coslat=cos(slat)**2
sinlat=sin(2.*slat)

eta=0.
eta=eta+ba(1)*aa(1)*sinlat*cos(ltime+mn-slon)             ! K1 (note sign change)
eta=eta-ba(2)*aa(2)*sinlat*cos(ltime-mn-slon)             ! O1
eta=eta-ba(3)*aa(3)*sinlat*cos(stime-sn-slon)             ! P1
eta=eta-ba(4)*aa(4)*sinlat*cos(ltime+pn-slon)             ! Q1
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
! this fill version is for a simply connected 2x2 grid
subroutine lfill(sc,sx)

implicit none

integer ncount
real, intent(in) :: sx
real, dimension(0:1,0:1), intent(inout) :: sc
real, dimension(0:1,0:1) :: nc
logical globc

ncount=count(sc>-990.)

if (ncount>=4) return

if (ncount<=0) then
  sc=sx
  return
end if

nc=sc
globc=.false.
call trial(nc,sc,0,0,.false.,.true.,.true.,.false.,globc)
call trial(nc,sc,1,0,.false.,.false.,.true.,.true.,globc)
call trial(nc,sc,0,1,.true.,.true.,.false.,.false.,globc)
call trial(nc,sc,1,1,.true.,.false.,.false.,.true.,globc)
sc=nc
if (globc) then
  call trial(nc,sc,0,0,.false.,.true.,.true.,.false.,globc)
  call trial(nc,sc,1,0,.false.,.false.,.true.,.true.,globc)
  call trial(nc,sc,0,1,.true.,.true.,.false.,.false.,globc)
  call trial(nc,sc,1,1,.true.,.false.,.false.,.true.,globc)
  sc=nc
end if

return
end subroutine lfill

! This routine averages adjacent points for an average fill
subroutine trial(nc,sc,i,j,ln,le,ls,lw,globc)

implicit none

integer, intent(in) :: i,j
integer nec
real, dimension(0:1,0:1), intent(inout) :: nc
real, dimension(0:1,0:1), intent(in) :: sc
real new
logical, intent(inout) :: globc
logical, intent(in) :: ln,le,ls,lw

if (nc(i,j)>-990.) return

new=0.
nec=0
if (ln) then
  if (sc(i,j-1)>-990.) then
    new=new+sc(i,j-1)
    nec=nec+1
  end if
end if
if (le) then
  if (sc(i+1,j)>-990.) then
    new=new+sc(i+1,j)
    nec=nec+1
  end if
end if
if (ls) then
  if (sc(i,j+1)>-990.) then
    new=new+sc(i,j+1)
    nec=nec+1
  end if
end if
if (lw) then
  if (sc(i-1,j)>-990.) then
    new=new+sc(i-1,j)
    nec=nec+1
  end if
end if
if (nec>0) then
  nc(i,j)=new/real(nec)
else
  globc=.true.
end if
  
return
end subroutine trial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! test for leap year

subroutine mloleap(tyear,ttest)

implicit none

integer, intent(in) :: tyear
logical, intent(out) :: ttest

ttest=.false.
if (mod(tyear,4)==0)   ttest=.true.
if (mod(tyear,100)==0) ttest=.false.
if (mod(tyear,400)==0) ttest=.true.

return
end subroutine mloleap

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Advect ice using upwind scheme

subroutine upwindadv(dumc,niu,niv,spnet)

use indices_m
use map_m

implicit none

include 'newmpar.h'
include 'parm.h'

real, dimension(ifull+iextra), intent(inout) :: dumc
real, dimension(ifull+iextra), intent(in) :: niu,niv,spnet
real, dimension(ifull) :: odum,dumd 

dumd=dumc(1:ifull)
odum=0.5*dt*(niu(iwu)*(dumc(1:ifull)+dumc(iw))    -abs(niu(iwu))*(dumc(1:ifull)-dumc(iw))    )*emu(iwu)/ds
where (spnet(iw)>0.)
  odum=min(odum,dumc(iw)*max(niu(iwu)*emu(iwu),0.)/(ds*spnet(iw)))
end where
where (spnet(1:ifull)>0.)
  odum=max(odum,dumc(1:ifull)*min(niu(iwu)*emu(iwu),0.)/(ds*spnet(1:ifull)))
end where
dumd=dumd+odum
odum=-0.5*dt*(niu(1:ifull)*(dumc(1:ifull)+dumc(ie))+abs(niu(1:ifull))*(dumc(1:ifull)-dumc(ie)))*emu(1:ifull)/ds
where (spnet(ie)>0.)
  odum=min(odum,-dumc(ie)*min(niu(1:ifull)*emu(1:ifull),0.)/(ds*spnet(ie)))
end where
where (spnet(1:ifull)>0.)
  odum=max(odum,-dumc(1:ifull)*max(niu(1:ifull)*emu(1:ifull),0.)/(ds*spnet(1:ifull)))
end where
dumd=dumd+odum  
odum=0.5*dt*(niv(isv)*(dumc(1:ifull)+dumc(is))    -abs(niv(isv))*(dumc(1:ifull)-dumc(is))    )*emv(isv)/ds
where (spnet(is)>0.)
  odum=min(odum,dumc(is)*max(niv(isv)*emv(isv),0.)/(ds*spnet(is)))
end where
where (spnet(1:ifull)>0.)
  odum=max(odum,dumc(1:ifull)*min(niv(isv)*emv(isv),0.)/(ds*spnet(1:ifull)))
end where
dumd=dumd+odum
odum=-0.5*dt*(niv(1:ifull)*(dumc(1:ifull)+dumc(in))+abs(niv(1:ifull))*(dumc(1:ifull)-dumc(in)))*emv(1:ifull)/ds
where (spnet(in)>0.)
  odum=min(odum,-dumc(in)*min(niv(1:ifull)*emv(1:ifull),0.)/(ds*spnet(in)))
end where
where (spnet(1:ifull)>0.)
  odum=max(odum,-dumc(1:ifull)*max(niv(1:ifull)*emv(1:ifull),0.)/(ds*spnet(1:ifull)))
end where
dumd=dumd+odum
dumc(1:ifull)=dumd
dumc(1:ifull)=max(dumc(1:ifull),0.)
  
return
end subroutine upwindadv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Use SOR to solve for free surface and ice pressure

subroutine mlosor(tol,itol,neta,sue,svn,suw,svs,pue,pvn,puw,pvs,que,qvn,quw,qvs,                    &
                  pdiv,pdivb,sdiv,sdivb,odiv,odivb,qdiv,qdivb,xps,                                  &
                  ipice,ibu,ibv,idu,idv,niu,niv,sicedep,ipmax,totits,itc,maxglobseta,maxglobip)

use cc_mpi
use indices_m
use map_m

implicit none

include 'newmpar.h'
include 'parm.h'

integer, intent(out) :: totits,itc
!integer itstest,itsave1,itsave2
integer nx,ifx,ll,ierr
integer iq
real, intent(in) :: tol,itol
real, intent(out) :: maxglobseta,maxglobip
real maxloclseta,maxloclip
!real gd,ci,itserr1,itserr2
real, dimension(ifull+iextra), intent(inout) :: neta
real, dimension(ifull), intent(in) :: sue,svn,suw,svs
real, dimension(ifull), intent(in) :: pue,pvn,puw,pvs
real, dimension(ifull), intent(in) :: que,qvn,quw,qvs
real, dimension(ifull), intent(in) :: pdiv,pdivb,sdiv,sdivb,odiv,odivb,qdiv,qdivb,xps
real, dimension(ifull+iextra), intent(inout) :: ipice
real, dimension(ifull+iextra), intent(in) :: ibu,ibv,idu,idv
real, dimension(ifull+iextra), intent(in) :: niu,niv
real, dimension(ifull), intent(in) :: sicedep
real, dimension(ifull+iextra), intent(in) :: ipmax
real, dimension(ifull,2) :: zz,zzn,zzs,zze,zzw,rhs
real, dimension(ifull) :: yy,yyn,yys,yye,yyw
real, dimension(ifull) :: hh
real, dimension(ifull) :: au,bu,cu,seta,setab,setac,nip
real, dimension(ifull,maxcolour,2) :: zzc,zznc,zzsc,zzec,zzwc,rhsc
real, dimension(ifull,maxcolour) :: yyc,yync,yysc,yyec,yywc
real, dimension(ifull,maxcolour) :: hhc
real, dimension(ifull+iextra,2) :: dumc
real, dimension(2) :: dume,dumf

integer, parameter :: llmax      = 400 ! Iterations for calculating surface height

! Use SOR as it converge faster that Conjugate Gradient (problems with coastlines?)
! and it is easy to include the non-linear terms and constraints.  Possibly combine
! with a multi-grid method to improve convergence at large scales, although it is
! more likely that the dynamical core will be replaced with JLM's mass flux scheme
!itsave2=0
!itserr2=9.E9
!itstest=1
!itc=0

dumc(1:ifull,1)=neta(1:ifull)
dumc(1:ifull,2)=ipice(1:ifull)
call bounds(dumc(:,1:2))
neta(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,1)
ipice(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,2)

! zz*(DIV^2 neta) + yy*neta*(DIV^2 neta) + hh*neta = rhs

! ocean
zzn(:,1)= (1.+ocneps)*0.5*dt*(qvn+svn+pvn)*ddv(1:ifull)*em(1:ifull)*em(1:ifull)/ds
zzs(:,1)=-(1.+ocneps)*0.5*dt*(qvs+svs+pvs)*ddv(isv)    *em(1:ifull)*em(1:ifull)/ds
zze(:,1)= (1.+ocneps)*0.5*dt*(que+sue+pue)*ddu(1:ifull)*em(1:ifull)*em(1:ifull)/ds
zzw(:,1)=-(1.+ocneps)*0.5*dt*(quw+suw+puw)*ddu(iwu)    *em(1:ifull)*em(1:ifull)/ds
zz(:,1) = (1.+ocneps)*0.5*dt*(qdivb+sdivb+pdivb)

yyn= (1.+ocneps)*0.5*dt*(qvn+svn+pvn)*em(1:ifull)*em(1:ifull)/ds
yys=-(1.+ocneps)*0.5*dt*(qvs+svs+pvs)*em(1:ifull)*em(1:ifull)/ds
yye= (1.+ocneps)*0.5*dt*(que+sue+pue)*em(1:ifull)*em(1:ifull)/ds
yyw=-(1.+ocneps)*0.5*dt*(quw+suw+puw)*em(1:ifull)*em(1:ifull)/ds
yy = (1.+ocneps)*0.5*dt*(qdiv+sdiv+pdiv)

hh     =1.+(1.+ocneps)*0.5*dt*(odiv                                                           &
        +(pvn*dd(in)-pvs*dd(is)+pue*dd(ie)-puw*dd(iw))*em(1:ifull)**2/ds+pdiv*dd(1:ifull))
rhs(:,1)=xps-(1.+ocneps)*0.5*dt*(odivb                                                        &
        +(pvn*ddv(1:ifull)*dd(in)-pvs*ddv(isv)*dd(is)+pue*ddu(1:ifull)*dd(ie)-puw*ddu(iwu)*dd(iw))*em(1:ifull)**2/ds+pdivb*dd(1:ifull))

! ice
zzn(:,2)=(-idv(1:ifull)*0.5-ibv(1:ifull))/ds
zzs(:,2)=(-idv(isv)*0.5    -ibv(isv)    )/ds
zze(:,2)=(-idu(1:ifull)*0.5-ibu(1:ifull))/ds
zzw(:,2)=(-idu(iwu)*0.5    -ibu(iwu)    )/ds
zz(:,2) =(ibu(1:ifull)+ibu(iwu)+ibv(1:ifull)+ibv(isv)-0.5*(idu(1:ifull)+idu(iwu)+idv(1:ifull)+idv(isv)))/ds

rhs(:,2)=min(niu(1:ifull)/emu(1:ifull)-niu(iwu)/emu(iwu)+niv(1:ifull)/emv(1:ifull)-niv(isv)/emv(isv),0.)

do nx=1,maxcolour
  ifx=ifullx(nx)
  zznc(1:ifx,nx,:)=zzn(iqx(1:ifx,nx),:)
  zzsc(1:ifx,nx,:)=zzs(iqx(1:ifx,nx),:)
  zzec(1:ifx,nx,:)=zze(iqx(1:ifx,nx),:)
  zzwc(1:ifx,nx,:)=zzw(iqx(1:ifx,nx),:)
  zzc(1:ifx,nx,:) = zz(iqx(1:ifx,nx),:)
  yync(1:ifx,nx)  =yyn(iqx(1:ifx,nx))
  yysc(1:ifx,nx)  =yys(iqx(1:ifx,nx))
  yyec(1:ifx,nx)  =yye(iqx(1:ifx,nx))
  yywc(1:ifx,nx)  =yyw(iqx(1:ifx,nx))
  yyc(1:ifx,nx)   = yy(iqx(1:ifx,nx))
  hhc(1:ifx,nx)   = hh(iqx(1:ifx,nx))
  rhsc(1:ifx,nx,:)=rhs(iqx(1:ifx,nx),:)
end do

do ll=1,llmax

  do nx=1,maxcolour
  
    ifx=ifullx(nx)

    ! ocean
    ! 5-point version -----------------------------------------------

    ! For now, assume Boussinesq fluid and treat density in continuity equation as constant
    au(1:ifx)=yyc(1:ifx,nx)
    bu(1:ifx)=zzc(1:ifx,nx,1)+hhc(1:ifx,nx)+yync(1:ifx,nx)*neta(iqn(1:ifx,nx))+yysc(1:ifx,nx)*neta(iqs(1:ifx,nx)) &
                                           +yyec(1:ifx,nx)*neta(iqe(1:ifx,nx))+yywc(1:ifx,nx)*neta(iqw(1:ifx,nx))
    cu(1:ifx)=-rhsc(1:ifx,nx,1)+zznc(1:ifx,nx,1)*neta(iqn(1:ifx,nx))+zzsc(1:ifx,nx,1)*neta(iqs(1:ifx,nx)) &
                               +zzec(1:ifx,nx,1)*neta(iqe(1:ifx,nx))+zzwc(1:ifx,nx,1)*neta(iqw(1:ifx,nx))
    
    ! alternative form
    setac(1:ifx)=-neta(iqx(1:ifx,nx))-2.*cu(1:ifx)/(bu(1:ifx)+sqrt(bu(1:ifx)*bu(1:ifx)-4.*au(1:ifx)*cu(1:ifx)))
    
    ! The following expression limits the minimum depth
    ! (should not occur for typical eta values)
    seta(iqx(1:ifx,nx))=max(setac(1:ifx),-(neta(iqx(1:ifx,nx))+dd(iqx(1:ifx,nx)))) ! this should become a land point
    seta(iqx(1:ifx,nx))=seta(iqx(1:ifx,nx))*ee(iqx(1:ifx,nx))
    neta(iqx(1:ifx,nx))=neta(iqx(1:ifx,nx))+seta(iqx(1:ifx,nx))

    ! ice
    ! 5-point version -------------------------------------------------
 
    bu(1:ifx)=zzc(1:ifx,nx,2)
    cu(1:ifx)=-rhsc(1:ifx,nx,2)+zznc(1:ifx,nx,2)*ipice(iqn(1:ifx,nx))+zzsc(1:ifx,nx,2)*ipice(iqs(1:ifx,nx)) &
                               +zzec(1:ifx,nx,2)*ipice(iqe(1:ifx,nx))+zzwc(1:ifx,nx,2)*ipice(iqw(1:ifx,nx))
    nip(1:ifx)=0.
    where (bu(1:ifx)/=0.)
      nip(1:ifx)=-cu(1:ifx)/bu(1:ifx)
    end where
 
    ! cavitating fluid
    nip(1:ifx)=max(min(nip(1:ifx),ipmax(iqx(1:ifx,nx))),0.)

    setab(iqx(1:ifx,nx))=nip(1:ifx)-ipice(iqx(1:ifx,nx))
    ipice(iqx(1:ifx,nx))=ipice(iqx(1:ifx,nx))+setab(iqx(1:ifx,nx))

    dumc(iqx(1:ifx,nx),1)=neta(iqx(1:ifx,nx))
    dumc(iqx(1:ifx,nx),2)=ipice(iqx(1:ifx,nx))
    call bounds_colour(dumc(:,1:2),nx)
    neta(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,1)
    ipice(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,2)
  
  end do
  
  maxloclseta=maxval(abs(seta))
  maxloclip=maxval(abs(setab))

  !if (ll>=itstest) then
    dume(1)=maxloclseta
    dume(2)=maxloclip
    call ccmpi_allreduce(dume(1:2),dumf(1:2),"max",comm_world)
    maxglobseta=dumf(1)
    maxglobip=dumf(2)

    if (maxglobseta<tol.and.maxglobip<itol) exit
    
  !  itsave1=itsave2
  !  itsave2=ll
  !  itserr1=itserr2
  !  itserr2=log10(maxglobseta)
  !  
  !  gd=(itserr2-itserr1)/real(itsave2-itsave1)
  !  ci=itserr2-gd*real(itsave2)
  !  if (gd/=0.) then
  !    itstest=nint((log10(tol)-ci)/gd)
  !    itstest=max(itstest,ll+1)
  !  else
  !    itstest=ll+1
  !  end if
  !  itc=itc+1    
  !end if

end do
totits=ll
itc=totits

return
end subroutine mlosor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Use multi-grid to solve for free surface and ice pressure

subroutine mlomg(tol,itol,neta,sue,svn,suw,svs,pue,pvn,puw,pvs,que,qvn,quw,qvs,                     &
                  pdiv,pdivb,sdiv,sdivb,odiv,odivb,qdiv,qdivb,xps,                                  &
                  ipice,ibu,ibv,idu,idv,niu,niv,sicedep,ipmax,totits,itc,maxglobseta,maxglobip)

use indices_m
use mgsolve
use map_m

implicit none

include 'newmpar.h'
include 'parm.h'

integer, intent(out) :: totits,itc
real, intent(in) :: tol,itol
real, intent(out) :: maxglobseta,maxglobip
real, dimension(ifull+iextra), intent(inout) :: neta,ipice
real, dimension(ifull), intent(in) :: sue,svn,suw,svs
real, dimension(ifull), intent(in) :: pue,pvn,puw,pvs
real, dimension(ifull), intent(in) :: que,qvn,quw,qvs
real, dimension(ifull), intent(in) :: pdiv,pdivb,sdiv,sdivb,odiv,odivb,qdiv,qdivb,xps
real, dimension(ifull+iextra), intent(in) :: ibu,ibv,idu,idv
real, dimension(ifull+iextra), intent(in) :: niu,niv
real, dimension(ifull), intent(in) :: sicedep
real, dimension(ifull+iextra), intent(in) :: ipmax
real, dimension(ifull,2) :: zz,zzn,zzs,zze,zzw,rhs
real, dimension(ifull) :: yy,yyn,yys,yye,yyw
real, dimension(ifull) :: hh

call mgsor_init

! zz*(DIV^2 neta) + yy*neta*(DIV^2 neta) + hh*neta = rhs

! ocean
zzn(:,1)= (1.+ocneps)*0.5*dt*(qvn+svn+pvn)*ddv(1:ifull)*em(1:ifull)*em(1:ifull)/ds
zzs(:,1)=-(1.+ocneps)*0.5*dt*(qvs+svs+pvs)*ddv(isv)    *em(1:ifull)*em(1:ifull)/ds
zze(:,1)= (1.+ocneps)*0.5*dt*(que+sue+pue)*ddu(1:ifull)*em(1:ifull)*em(1:ifull)/ds
zzw(:,1)=-(1.+ocneps)*0.5*dt*(quw+suw+puw)*ddu(iwu)    *em(1:ifull)*em(1:ifull)/ds
zz(:,1) = (1.+ocneps)*0.5*dt*(qdivb+sdivb+pdivb)

yyn= (1.+ocneps)*0.5*dt*(qvn+svn+pvn)*em(1:ifull)*em(1:ifull)/ds
yys=-(1.+ocneps)*0.5*dt*(qvs+svs+pvs)*em(1:ifull)*em(1:ifull)/ds
yye= (1.+ocneps)*0.5*dt*(que+sue+pue)*em(1:ifull)*em(1:ifull)/ds
yyw=-(1.+ocneps)*0.5*dt*(quw+suw+puw)*em(1:ifull)*em(1:ifull)/ds
yy = (1.+ocneps)*0.5*dt*(qdiv+sdiv+pdiv)

hh     =1.+(1.+ocneps)*0.5*dt*(odiv                                                           &
        +(pvn*dd(in)-pvs*dd(is)+pue*dd(ie)-puw*dd(iw))*em(1:ifull)*em(1:ifull)/ds+pdiv*dd(1:ifull))
rhs(:,1)=xps(1:ifull)-(1.+ocneps)*0.5*dt*(odivb                                               &
        +(pvn*ddv(1:ifull)*dd(in)-pvs*ddv(isv)*dd(is)+pue*ddu(1:ifull)*dd(ie)-puw*ddu(iwu)*dd(iw))*em(1:ifull)*em(1:ifull)/ds+pdivb*dd(1:ifull))

! zz*(DIV^2 ipice) + yy*ipice = rhs

! ice
zzn(:,2)=(-idv(1:ifull)*0.5-ibv(1:ifull))/ds
zzs(:,2)=(-idv(isv)*0.5    -ibv(isv)    )/ds
zze(:,2)=(-idu(1:ifull)*0.5-ibu(1:ifull))/ds
zzw(:,2)=(-idu(iwu)        -ibu(iwu)    )/ds
zz(:,2) =(-0.5*(idu(1:ifull)+idu(iwu)+idv(1:ifull)+idv(isv))+ibu(1:ifull)+ibu(iwu)+ibv(1:ifull)+ibv(isv))/ds

rhs(:,2)=min(niu(1:ifull)/emu(1:ifull)-niu(iwu)/emu(iwu)+niv(1:ifull)/emv(1:ifull)-niv(isv)/emv(isv),0.)

where (zz(:,2)>=0.)
  zz(:,2) =-dt/(ds*10.) ! 10 is the minimum imass
  zzn(:,2)=0.
  zzs(:,2)=0.
  zze(:,2)=0.
  zzw(:,2)=0.
  rhs(:,2)=0.
end where

call mgmlo(neta,ipice,yy,yyn,yys,yye,yyw,zz,zzn,zzs,zze,zzw,hh,rhs,tol,itol,totits,maxglobseta,maxglobip, &
           ipmax,ee,dd)
itc=totits

return
end subroutine mlomg


end module mlodynamics
