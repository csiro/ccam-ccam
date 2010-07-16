
! This is a 1D, mixed layer ocean model for ensemble regioal climate simulations based on Large, et al (1994)
! (i.e., adapted from the GFDL code).  This code is also used for modelling lakes in CCAM.

! This version has a relatively thin 1st layer (0.5m) so as to reproduce a diurnal cycle in SST.  It also
! supports sea ice based on O'Farrell's sea ice model from Mk3.5.

! This version can assimilate SSTs from GCMs, using a convolution based digital filter (see nestin.f),
! which avoids problems with complex land-sea boundary conditions

! Order of calls:
!    call mloinit
!    call mloload
!    ...
!    main loop
!       ...
!       call mloalb4
!       ...
!       call mloeval
!       ...
!       call mloscrnout
!       ...
!    end loop
!    ...
!    call mlosave
!    call mloend


! To do list:
!   - Include ice for frozen lakes

module mlo

implicit none

private
public mloinit,mloend,mloeval,mloimport,mloexport,mloload,mlosave,mloregrid,mlodiag,mloalb2,mloalb4, &
       mloscrnout,wlev,depth,mlodwn,mlootf,ocndwn,ocnotf

type tdata
  real temp,sal,u,v
end type tdata
type tice
  real tn(0:3)
  real dic,fracice
end type tice
type tprog2
  real mixdepth,bf
  real watervisdiralb,watervisdifalb,waternirdiralb,waternirdifalb
  real icevisdiralb,icevisdifalb,icenirdiralb,icenirdifalb
  real zo,cd,fg,eg,wetfac
  real zoice,cdice,fgice,egice,wetfacice
  real tscrn,uscrn,qgscrn,u10
  integer mixind
end type tprog2
type tatm
  real sg,rg,rnd,f,vnratio,fbvis,fbnir
  real ps,p1,u,v,temp,qg,zmin
end type tatm
type tdiag2
  real b0,ustar,wu0,wv0,wt0,ws0,taux,tauy
end type tdiag2
type tdiag3
  real rho,nsq,rad,alpha,beta
end type tdiag3

! parameters
integer, parameter :: wlev = 20
integer, parameter :: iqx = 500
! model arrays
integer, save :: wfull
integer, dimension(:), allocatable, save :: wgrid
real, dimension(:,:), allocatable, save :: depth,dz
real, dimension(:,:), allocatable, save :: depth_hl
real, dimension(:,:), allocatable, save :: dz_hl
real, dimension(:,:,:), allocatable, save :: mlodwn,mlootf ! These variables are for CCAM onthefly.f
real, dimension(:), allocatable, save :: ocndwn,ocnotf     ! These variables are for CCAM onthefly.f
type(tdata), dimension(:,:), allocatable, save :: water
type(tice), dimension(:), allocatable, save :: ice
type(tprog2), dimension(:), allocatable, save :: pg

! mode
integer, parameter :: incradbf  = 1 ! include shortwave in buoyancy forcing
integer, parameter :: incradgam = 0 ! include shortwave in non-local term
integer, parameter :: salrelax  = 1 ! relax salinity to 34.72 PSU (used for single column mode)

! max depth
real, parameter :: mxd    = 977.6   ! Max depth (m)
real, parameter :: mindep = 1.      ! Thickness of first layer (m)

! model parameters
real, parameter :: ric     = 0.3    ! Critical Ri for diagnosing mixed layer depth
real, parameter :: epsilon = 0.1

! radiation parameters
!real, parameter :: r_1 = 0.58       ! VIS/(VIS+NIR) ratio - now taken from host radiation
!real, parameter :: r_2 = 1.-r_1
real, parameter :: mu_1 = 23.       ! VIS depth (m)
real, parameter :: mu_2 = 0.35      ! NIR depth (m)

! physical parameters
real, parameter :: vkar=0.4               ! von Karman constant
real, parameter :: lv=2.501e6             ! Latent heat of vaporisation (J kg^-1)
real, parameter :: lf=3.337e5             ! Latent heat of fusion (J kg^-1)
real, parameter :: ls=lv+lf               ! Latent heat of sublimation (J kg^-1)
real, parameter :: cp0=3990.              ! heat capacity of mixed layer (J kg^-1 K^-1)
real, parameter :: grav=9.80              ! graviational constant (m/s^2)
real, parameter :: sbconst=5.67e-8        ! Stefan-Boltzmann constant
real, parameter :: cdbot=2.4E-3           ! bottom drag coefficent
real, parameter :: tfi=271.2              ! Freezing temperature
real, parameter :: cp=1004.64             ! Specific heat of dry air at const P
real, parameter :: rdry=287.04            ! Specific gas const for dry air
real, parameter :: rvap=461.5             ! Gas constant for water vapor
! stability function parameters
real, parameter :: bprm=5.                ! 4.7 in rams
real, parameter :: chs=2.6                ! 5.3 in rams
real, parameter :: cms=5.                 ! 7.4 in rams
real, parameter :: fmroot=0.57735
real, parameter :: rimax=(1./fmroot-1.)/bprm
! momentum flux parameters
real, parameter :: charnck=0.018
real, parameter :: chn10=0.00125

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
allocate(ice(wfull))
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

water%temp=288.          ! K
water%sal=35.            ! PSU
water%u=0.               ! m/s
water%v=0.               ! m/s

ice%dic=0.              ! m
ice%fracice=0.          ! %
do ii=0,3
  ice(:)%tn(ii)=271.2   ! K
end do

pg%mixdepth=100. ! m
pg%bf=0.
pg%watervisdiralb=0.06
pg%watervisdifalb=0.06
pg%waternirdiralb=0.06
pg%waternirdifalb=0.06
pg%icevisdiralb=0.65
pg%icevisdifalb=0.65
pg%icenirdiralb=0.65
pg%icenirdifalb=0.65
pg%zo=0.001
pg%cd=0.
pg%fg=0.
pg%eg=0.
pg%wetfac=1.
pg%zoice=0.001
pg%cdice=0.
pg%fgice=0.
pg%egice=0.
pg%wetfacice=1.
pg%tscrn=273.2
pg%qgscrn=0.
pg%uscrn=0.
pg%u10=0.
pg%mixind=wlev-1

! MLO
!depth = (/   0.5,   1.9,   4.3,   8.5,  15.3,  25.3,  39.3,  57.9,  81.9, 112.1 &
!           149.1, 193.7, 246.5, 308.3, 379.9, 461.9, 555.1, 660.1, 777.7, 908.7 /)
! Mk3.5
!depth = (/   5.0,  15.0,  28.2,  42.0,  59.7,  78.5, 102.1, 127.9, 159.5, 194.6, &
!           237.0, 284.7, 341.7, 406.4, 483.2, 570.9, 674.9, 793.8, 934.1 /)

smxd=maxval(depin(wgrid))
smnd=minval(depin(wgrid))

do iqw=1,wfull
  call vgrid(depin(wgrid(iqw)),depth(iqw,:),depth_hl(iqw,:))
  if (smxd.eq.depin(wgrid(iqw))) then
    write(6,*) "MLO max depth ",depth(iqw,:)
    smxd=smxd+10.
  end if
  if (smnd.eq.depin(wgrid(iqw))) then
    write (6,*) "MLO min depth ",depth(iqw,:)
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

subroutine vgrid(depin,depthout,depth_hlout)

implicit none

integer ii
real, intent(in) :: depin
real, dimension(wlev), intent(out) :: depthout
real, dimension(wlev+1), intent(out) :: depth_hlout
real dd,x,al,bt

dd=min(mxd,max(mindep,depin))
x=real(wlev)
if (dd.gt.x) then
  al=(mindep*x-dd)/(x-x*x*x)
  bt=(mindep*x*x*x-dd)/(x*x*x-x)
  do ii=1,wlev+1
    x=real(ii-1)
    depth_hlout(ii)=al*x*x*x+bt*x ! ii is for half level ii-0.5
  end do
else
  depth_hlout(1)=0.
  depth_hlout(2)=mindep
  do ii=3,wlev+1
    x=(dd-mindep)*real(ii-2)/real(wlev-1)+mindep
    depth_hlout(ii)=x
  end do
end if
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
deallocate(ice)
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
!ice(:)%dic=datainb(wgrid,1)
!ice(:)%fracice=datainb(wgrid,2)
!do ii=0,3
!  ice(:)%tn(ii)=datainc(wgrid,ii,1)
!end do

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
!dataoutb(wgrid,1)=ice(:)%dic
!dataoutb(wgrid,2)=ice(:)%fracice
!do ii=0,3
!  dataoutc(wgrid,ii,1)=ice(:)%tn(ii)
!end do
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
! Return mixed layer depth

subroutine mlodiag(ifull,mld,diag)

implicit none

integer, intent(in) :: ifull,diag
real, dimension(ifull), intent(out) :: mld

mld=0.
if (wfull.eq.0) return
mld(wgrid)=pg(:)%mixdepth

return
end subroutine mlodiag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate screen diagnostics

subroutine mloscrnout(ifull,tscrn,qgscrn,uscrn,u10,diag)

implicit none

integer, intent(in) :: ifull,diag
real, dimension(ifull), intent(inout) :: tscrn,qgscrn,uscrn,u10

if (wfull.eq.0) return
tscrn(wgrid)=pg%tscrn
qgscrn(wgrid)=pg%qgscrn
uscrn(wgrid)=pg%uscrn
u10(wgrid)=pg%u10

return
end subroutine mloscrnout

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate albedo data (VIS + NIR)

subroutine mloalb2(istart,ifull,land,coszro,fracice,ovisalb,oniralb,diag)

implicit none

integer, intent(in) :: istart,ifull,diag
integer ifinish,ib,ie,iqw
real, dimension(ifull), intent(in) :: coszro,fracice
real, dimension(ifull), intent(inout) :: ovisalb,oniralb
real, dimension(ifull) :: watervis,waternir,icevis,icenir
logical, dimension(ifull), intent(in) :: land

watervis=.05/(coszro+0.15)
waternir=.05/(coszro+0.15)
icevis=0.85
icenir=0.45
where (.not.land)
  ovisalb=icevis*fracice+(1.-fracice)*watervis
  oniralb=icenir*fracice+(1.-fracice)*waternir
endwhere

if (wfull.eq.0) return

ifinish=istart+ifull-1
ib=wfull+1
do iqw=wfull,1,-1
  if (wgrid(iqw).lt.istart) exit
  ib=iqw
end do
if (ib.gt.wfull) return
ie=-1
do iqw=ib,wfull
  if (wgrid(iqw).gt.ifinish) exit
  ie=iqw
end do
if (ie.lt.ib) return

pg(ib:ie)%watervisdiralb=watervis(wgrid(ib:ie)-istart+1)
pg(ib:ie)%watervisdifalb=watervis(wgrid(ib:ie)-istart+1)
pg(ib:ie)%waternirdiralb=waternir(wgrid(ib:ie)-istart+1)
pg(ib:ie)%waternirdifalb=waternir(wgrid(ib:ie)-istart+1)
pg(ib:ie)%icevisdiralb=icevis(wgrid(ib:ie)-istart+1)
pg(ib:ie)%icevisdifalb=icevis(wgrid(ib:ie)-istart+1)
pg(ib:ie)%icenirdiralb=icenir(wgrid(ib:ie)-istart+1)
pg(ib:ie)%icenirdifalb=icenir(wgrid(ib:ie)-istart+1)

return
end subroutine mloalb2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate albedo data (VIS-DIR, VIS-DIF, NIR-DIR & NIR-DIF)

subroutine mloalb4(istart,ifull,land,coszro,fracice,ovisdir,ovisdif,onirdir,onirdif,diag)

implicit none

integer, intent(in) :: istart,ifull,diag
integer ifinish,ib,ie,iqw
real, dimension(ifull), intent(in) :: coszro,fracice
real, dimension(ifull), intent(inout) :: ovisdir,ovisdif,onirdir,onirdif
real, dimension(ifull) :: watervisdiralb,watervisdifalb,waternirdiralb,waternirdifalb
real, dimension(ifull) :: icevisdiralb,icevisdifalb,icenirdiralb,icenirdifalb
logical, dimension(ifull), intent(in) :: land

where (coszro.gt.0.)
  watervisdiralb=0.026/(coszro**1.7+0.065)+0.15*(coszro-0.1)*(coszro-0.5)*(coszro-1.)
elsewhere
  watervisdiralb=0.3925 
endwhere
watervisdifalb=0.06
waternirdiralb=watervisdiralb
waternirdifalb=0.06
icevisdiralb=0.85
icevisdifalb=0.85
icenirdiralb=0.45
icenirdifalb=0.45
where (.not.land)
  ovisdir=fracice*icevisdiralb+(1.-fracice)*watervisdiralb
  ovisdif=fracice*icevisdifalb+(1.-fracice)*watervisdifalb
  onirdir=fracice*icenirdiralb+(1.-fracice)*waternirdiralb
  onirdif=fracice*icenirdifalb+(1.-fracice)*waternirdifalb
endwhere

if (wfull.eq.0) return

ifinish=istart+ifull-1
ib=wfull+1
do iqw=wfull,1,-1
  if (wgrid(iqw).lt.istart) exit
  ib=iqw
end do
if (ib.gt.wfull) return
ie=-1
do iqw=ib,wfull
  if (wgrid(iqw).gt.ifinish) exit
  ie=iqw
end do
if (ie.lt.ib) return

pg(ib:ie)%watervisdiralb=watervisdiralb(wgrid(ib:ie)-istart+1)
pg(ib:ie)%watervisdifalb=watervisdifalb(wgrid(ib:ie)-istart+1)
pg(ib:ie)%waternirdiralb=waternirdiralb(wgrid(ib:ie)-istart+1)
pg(ib:ie)%waternirdifalb=waternirdifalb(wgrid(ib:ie)-istart+1)
pg(ib:ie)%icevisdiralb=icevisdiralb(wgrid(ib:ie)-istart+1)
pg(ib:ie)%icevisdifalb=icevisdifalb(wgrid(ib:ie)-istart+1)
pg(ib:ie)%icenirdiralb=icenirdiralb(wgrid(ib:ie)-istart+1)
pg(ib:ie)%icenirdifalb=icenirdifalb(wgrid(ib:ie)-istart+1)

return
end subroutine mloalb4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Regrid MLO data

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

subroutine mloeval(ifull,sst,zo,cd,fg,eg,wetfac, &
                   dt,zmin,sg,rg,precp,uatm,vatm,temp,qg,ps,p1, &
                   f,visnirratio,fbvis,fbnir,diag)

implicit none

integer, intent(in) :: ifull,diag
real, intent(in) :: dt,zmin
real, dimension(ifull), intent(in) :: sg,rg,precp,f,uatm,vatm,temp,qg,ps,p1,visnirratio,fbvis,fbnir
real, dimension(ifull), intent(inout) :: sst,zo,cd,fg,eg,wetfac
real, dimension(wfull) :: workb
type(tatm), dimension(wfull) :: atm

if (wfull.eq.0) return

atm%sg=sg(wgrid)
atm%rg=rg(wgrid)
atm%rnd=precp(wgrid)
!atm%f=f(wgrid)
atm%f=0. ! turn off coriolis terms when no geostrophic term
atm%vnratio=visnirratio(wgrid)
atm%fbvis=fbvis(wgrid)
atm%fbnir=fbnir(wgrid)
atm%u=uatm(wgrid)-water(:,1)%u
atm%v=vatm(wgrid)-water(:,1)%v
atm%temp=temp(wgrid)
atm%qg=qg(wgrid)
atm%ps=ps(wgrid)
atm%p1=p1(wgrid)
atm%zmin=zmin

call mlocalc(dt,atm,diag)
call mloice(dt,atm,diag)
call scrncalc(atm,diag)

sst(wgrid)=(1.-ice%fracice)*water(:,1)%temp+ice%fracice*ice%tn(0)
workb=(1.-ice%fracice)/log(atm%zmin/pg%zo)**2+ice%fracice/log(atm%zmin/pg%zoice)**2
zo(wgrid)=atm%zmin*exp(-1./sqrt(workb))
cd(wgrid)=(1.-ice%fracice)*pg%cd+ice%fracice*pg%cdice
fg(wgrid)=(1.-ice%fracice)*pg%fg+ice%fracice*pg%fgice
eg(wgrid)=(1.-ice%fracice)*pg%eg+ice%fracice*pg%egice
wetfac(wgrid)=(1.-ice%fracice)*pg%wetfac+ice%fracice*pg%wetfacice

return
end subroutine mloeval

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MLO calcs for water (no ice)

subroutine mlocalc(dt,atm,diag)

implicit none

integer, intent(in) :: diag
integer ii,iqw,pos(1)
real, intent(in) :: dt
real, dimension(wfull,wlev) :: km,ks,gammas
real, dimension(wfull,wlev) :: rhs
real, dimension(wfull,2:wlev) :: aa
real, dimension(wfull,wlev) :: bb,dd
real, dimension(wfull,1:wlev-1) :: cc
real, dimension(wfull) :: xp,xm,dumt0,umag
type(tatm), dimension(wfull), intent(in) :: atm
type(tdiag2), dimension(wfull) :: dg2
type(tdiag3), dimension(wfull,wlev) :: dg3
type(tdata), dimension(wfull,wlev) :: new

call fluxcalc(dg2,atm,diag)        ! prepare fluxes
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
rhs(:,1)=(1./dz_hl(:,2))*(ks(:,2)*gammas(:,2)-ks(:,1)*gammas(:,1))
do ii=2,wlev-1
  where (ii.lt.pg%mixind)
    rhs(:,ii)=(0.5/dz(:,ii))*(ks(:,ii+1)*gammas(:,ii+1)-ks(:,ii-1)*gammas(:,ii-1)) ! non-local
  elsewhere (ii.eq.pg%mixind)
    rhs(:,ii)=(1./dz_hl(:,ii))*(-ks(:,ii-1)*gammas(:,ii-1))
  elsewhere
    rhs(:,ii)=0.
  end where
end do
where (wlev.eq.pg%mixind)
  rhs(:,wlev)=(1./dz_hl(:,wlev))*(-ks(:,wlev-1)*gammas(:,wlev-1))
elsewhere
  rhs(:,wlev)=0.
end where
cc(:,1)=-0.5*(ks(:,1)+ks(:,2))/(dz_hl(:,2)*dz(:,1))
bb(:,1)=1./dt-cc(:,1)
! use -ve for BC as depth is down
dd(:,1)=water(:,1)%temp/dt+rhs(:,1)*dumt0+dg3(:,1)%rad/dz(:,1)-dg2%wt0/dz(:,1)
do ii=2,wlev-1
  aa(:,ii)=-0.5*(ks(:,ii)+ks(:,ii-1))/(dz_hl(:,ii)*dz(:,ii))
  cc(:,ii)=-0.5*(ks(:,ii+1)+ks(:,ii))/(dz_hl(:,ii+1)*dz(:,ii))
  bb(:,ii)=1./dt-aa(:,ii)-cc(:,ii)
  dd(:,ii)=water(:,ii)%temp/dt+rhs(:,ii)*dumt0+dg3(:,ii)%rad/dz(:,ii)
end do
aa(:,wlev)=-0.5*(ks(:,wlev)+ks(:,wlev-1))/(dz_hl(:,wlev)*dz(:,wlev))
bb(:,wlev)=1./dt-aa(:,wlev)
dd(:,wlev)=water(:,wlev)%temp/dt+rhs(:,wlev)*dumt0+dg3(:,wlev)%rad/dz(:,wlev)
call thomas(new%temp,aa,bb,cc,dd)

! SALINITY
do ii=1,wlev
  dd(:,ii)=water(:,ii)%sal/dt+rhs(:,ii)*dg2%ws0
end do
dd(:,1)=dd(:,1)-dg2%ws0/dz(:,1)
if (salrelax.eq.1) then ! relax salinity
  dd=dd+(34.72-water%sal)/(3600.*24.*365.25*10.)
end if
call thomas(new%sal,aa,bb,cc,dd)

new%sal=max(0.,new%sal) ! MJT suggestion

! split U diffusion term
cc(:,1)=-0.5*(km(:,1)+km(:,2))/(dz_hl(:,2)*dz(:,1))
bb(:,1)=1./dt-cc(:,1)
dd(:,1)=water(:,1)%u/dt-dg2%wu0/dz(:,1)
do ii=2,wlev-1
  aa(:,ii)=-0.5*(km(:,ii)+km(:,ii-1))/(dz_hl(:,ii)*dz(:,ii))
  cc(:,ii)=-0.5*(km(:,ii+1)+km(:,ii))/(dz_hl(:,ii+1)*dz(:,ii))
  bb(:,ii)=1./dt-aa(:,ii)-cc(:,ii)
  dd(:,ii)=water(:,ii)%u/dt
end do
aa(:,wlev)=-0.5*(km(:,wlev)+km(:,wlev-1))/(dz_hl(:,wlev)*dz(:,wlev))
bb(:,wlev)=1./dt-aa(:,wlev)
dd(:,wlev)=water(:,wlev)%u/dt
umag=sqrt(water(:,wlev)%u*water(:,wlev)%u+water(:,wlev)%v*water(:,wlev)%v)
where (depth_hl(:,wlev+1).lt.(mxd-0.1)) ! bottom drag
  bb(:,wlev)=bb(:,wlev)+cdbot*umag/dz(:,wlev)
end where
call thomas(new%u,aa,bb,cc,dd)

! split V diffusion term
do ii=1,wlev
  dd(:,ii)=water(:,ii)%v/dt
end do
dd(:,1)=dd(:,1)-dg2%wv0/dz(:,1)
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

new%u=max(-100.,min(new%u,100.)) ! MJT suggestion
new%v=max(-100.,min(new%v,100.)) ! MJT suggestion

water=new

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
real, dimension(wfull) :: numh,wm1,dnumhdz,dwm1ds,g1m,dg1mds
real, dimension(wfull) :: nush,ws1,dnushdz,dws1ds,g1s,dg1sds
real xp,cg
type(tdiag2), dimension(wfull), intent(in) :: dg2
type(tdiag3), dimension(wfull,wlev), intent(in) :: dg3
real, parameter :: ri0 = 0.7
real, parameter :: nu0 = 50.E-4
real, parameter :: numw = 1.E-4
real, parameter :: nusw = 0.1E-4

ri(:,1)=dg3(:,1)%nsq*(2.*dz_hl(:,2))**2 &
           /max((water(:,1)%u-water(:,2)%u)**2 &
               +(water(:,1)%v-water(:,2)%v)**2,1.E-10)
do ii=2,wlev-1
  ri(:,ii)=dg3(:,ii)%nsq*(2.*dz(:,ii))**2 &
           /max((water(:,ii-1)%u-water(:,ii+1)%u)**2 &
               +(water(:,ii-1)%v-water(:,ii+1)%v)**2,1.E-10)
end do
ri(:,wlev)=dg3(:,wlev)%nsq*(2.*dz_hl(:,wlev))**2 &
           /max((water(:,wlev-1)%u-water(:,wlev)%u)**2 &
               +(water(:,wlev-1)%v-water(:,wlev)%v)**2,1.E-10)

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
! w - internal ave
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
  numh(iqw)=(1.-xp)*num(iqw,pg(iqw)%mixind)+xp*num(iqw,pg(iqw)%mixind+1)
  wm1(iqw)=(1.-xp)*wm(iqw,pg(iqw)%mixind)+xp*wm(iqw,pg(iqw)%mixind+1)
  dnumhdz(iqw)=-(num(iqw,pg(iqw)%mixind+1)-num(iqw,pg(iqw)%mixind))/dz_hl(iqw,pg(iqw)%mixind+1)
  dwm1ds(iqw)=-pg(iqw)%mixdepth*(wm(iqw,pg(iqw)%mixind+1)-wm(iqw,pg(iqw)%mixind))/dz_hl(iqw,pg(iqw)%mixind+1)
  nush(iqw)=(1.-xp)*nus(iqw,pg(iqw)%mixind)+xp*nus(iqw,pg(iqw)%mixind+1)
  ws1(iqw)=(1.-xp)*ws(iqw,pg(iqw)%mixind)+xp*ws(iqw,pg(iqw)%mixind+1)
  dnushdz(iqw)=-(nus(iqw,pg(iqw)%mixind+1)-nus(iqw,pg(iqw)%mixind))/dz_hl(iqw,pg(iqw)%mixind+1)
  dws1ds(iqw)=-pg(iqw)%mixdepth*(ws(iqw,pg(iqw)%mixind+1)-ws(iqw,pg(iqw)%mixind))/dz_hl(iqw,pg(iqw)%mixind+1)
end do

wm1=max(wm1,1.E-10)
ws1=max(ws1,1.E-10)
  
g1m=numh/(pg%mixdepth*wm1)
dg1mds=-dnumhdz/wm1-numh*dwm1ds/(pg%mixdepth*wm1*wm1)
g1s=nush/(pg%mixdepth*ws1)
dg1sds=-dnushdz/ws1-nush*dws1ds/(pg%mixdepth*ws1*ws1)
  
a2m=-2.+3.*g1m-dg1mds
a3m=1.-2.*g1m+dg1mds
a2s=-2.+3.*g1s-dg1sds
a3s=1.-2.*g1s+dg1sds

!--------------------------------------------------------------------
! combine
do ii=1,wlev
  where (ii.le.pg%mixind)
    sigma=depth(:,ii)/pg%mixdepth
    km(:,ii)=pg%mixdepth*wm(:,ii)*sigma*(1.+sigma*(a2m+a3m*sigma))
    ks(:,ii)=pg%mixdepth*ws(:,ii)*sigma*(1.+sigma*(a2m+a3m*sigma))
  elsewhere
    km(:,ii)=num(:,ii)
    ks(:,ii)=nus(:,ii)
  end where
end do

km=max(km,numw) ! MJT suggestion
ks=max(ks,nusw) ! MJT suggestion

! non-local term
! gammas is the same for temp and sal when double-diffusion is not employed
cg=10.*vkar*(98.96*vkar*epsilon)**(1./3.)
do ii=1,wlev
  where (pg%bf.lt.0.) ! unstable
    gammas(:,ii)=cg/max(ws(:,ii)*pg%mixdepth,1.E-10)
  elsewhere           ! stable
    gammas(:,ii)=0.
  end where
end do

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
pg%mixdepth=min(pg%mixdepth,depth(:,wlev))

! recalculate index for mixdepth
pg%mixind=wlev-1
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
uuu=max(ustar**3,1.E-10)
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
real, dimension(wfull) :: visalb,niralb
real, dimension(wfull,wlev) :: rs
real, parameter :: density = 1035.
type(tdiag2), dimension(wfull), intent(inout) :: dg2
type(tdiag3), dimension(wfull,wlev), intent(inout) :: dg3
type(tatm), dimension(wfull), intent(in) :: atm

do ii=1,wlev
  dg3(:,ii)%rho=density
end do

t = max(water(:,1)%temp-273.16,-2.)
s = max(water(:,1)%sal,0.)
t2 = t*t
t3 = t2*t
t4 = t3*t
t5 = t4*t
s2 = s*s
s3 = s2*s
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

do i=1,nits
  do ii=1,wlev
    t = max(water(:,ii)%temp-273.16,-2.)
    s = max(water(:,ii)%sal,0.01)
    p1 = grav*depth(:,ii)*dg3(:,ii)%rho/100000.
    t2 = t*t
    t3 = t2*t
    t4 = t3*t
    t5 = t4*t
    s2 = s*s
    s3 = s2*s
    p2 = p1*p1
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
dg3(:,1)%nsq=-(grav/dg3(:,1)%rho)*(dg3(:,1)%rho-dg3(:,2)%rho)/dz_hl(:,2)
do ii=2,wlev-1
  dg3(:,ii)%nsq=-(grav/dg3(:,ii)%rho)*0.5*(dg3(:,ii-1)%rho-dg3(:,ii+1)%rho)/dz(:,ii)
end do
dg3(:,wlev)%nsq=-(grav/dg3(:,wlev)%rho)*(dg3(:,wlev-1)%rho-dg3(:,wlev)%rho)/dz_hl(:,wlev)

! shortwave
! use -ve as depth is down
visalb=pg%watervisdiralb*atm%fbvis+pg%watervisdifalb*(1.-atm%fbvis)
niralb=pg%waternirdiralb*atm%fbnir+pg%waternirdifalb*(1.-atm%fbnir)
hlrb=0.5*(dg3(:,1)%rho+dg3(:,2)%rho)
dg3(:,1)%rad=-atm%sg/cp0*(((1.-visalb)*atm%vnratio*exp(-depth_hl(:,2)/mu_1)+ &
                           (1.-niralb)*(1.-atm%vnratio)*exp(-depth_hl(:,2)/mu_2))/hlrb- &
                          ((1.-visalb)*atm%vnratio*exp(-depth_hl(:,1)/mu_1)+ &
                           (1.-niralb)*(1.-atm%vnratio)*exp(-depth_hl(:,1)/mu_2))/rho0)
do ii=2,wlev-1 
  hlra=0.5*(dg3(:,ii-1)%rho+dg3(:,ii)%rho)
  hlrb=0.5*(dg3(:,ii)%rho+dg3(:,ii+1)%rho)
  dg3(:,ii)%rad=-atm%sg/cp0*(((1.-visalb)*atm%vnratio*exp(-depth_hl(:,ii+1)/mu_1)+ &
                              (1.-niralb)*(1.-atm%vnratio)*exp(-depth_hl(:,ii+1)/mu_2))/hlrb- &
                             ((1.-visalb)*atm%vnratio*exp(-depth_hl(:,ii)/mu_1)+ &
                              (1.-niralb)*(1.-atm%vnratio)*exp(-depth_hl(:,ii)/mu_2))/hlra)
end do
hlra=0.5*(dg3(:,wlev-1)%rho+dg3(:,wlev)%rho)
dg3(:,wlev)%rad=atm%sg/cp0*((1.-visalb)*atm%vnratio*exp(-depth_hl(:,wlev)/mu_1)+ &
                            (1.-niralb)*(1.-atm%vnratio)*exp(-depth_hl(:,wlev)/mu_2))/hlra ! remainder

! Boundary conditions (use rho at level 1 for consistancy with shortwave radiation)
dg2%wu0=-dg2%taux/dg3(:,1)%rho                                                ! BC
dg2%wv0=-dg2%tauy/dg3(:,1)%rho                                                ! BC
dg2%wt0=-(-pg%fg-pg%eg+atm%rg-sbconst*water(:,1)%temp**4)/(dg3(:,1)%rho*cp0)  ! BC
dg2%ws0=(atm%rnd-pg%eg/lv)*water(:,1)%sal/rs(:,1)                             ! BC ! negect ice for now (included in seaice using split form)

dg2%ustar=max(sqrt(sqrt(dg2%wu0*dg2%wu0+dg2%wv0*dg2%wv0)),1.E-10)
dg2%b0=-grav*(dg3(:,1)%alpha*dg2%wt0-dg3(:,1)%beta*dg2%ws0)
!dg2%b0=-grav*(-drho0dt*dg2%wt0-drho0ds*dg2%ws0)

return
end subroutine getrho

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate fluxes between MLO and atmosphere (from CCAM sflux.f)

subroutine fluxcalc(dg2,atm,diag)

implicit none

integer, intent(in) :: diag
integer it
type(tatm), dimension(wfull), intent(in) :: atm
type(tdiag2), dimension(wfull), intent(inout) :: dg2
real, dimension(wfull) :: qsat,ri,vmag,rho,srcp
real, dimension(wfull) :: fm,fh,con,consea,afroot,af,daf
real, dimension(wfull) :: den,dfm,sig,factch,root
real, dimension(wfull) :: z1onzt,aft
real ztv

vmag=sqrt(atm%u*atm%u+atm%v*atm%v)
sig=atm%p1/atm%ps
srcp=sig**(rdry/cp)
rho=atm%ps/(rdry*water(:,1)%temp)
ztv=exp(vkar/sqrt(chn10))/10.  ! proper inverse of ztsea
z1onzt=300.*rdry*(1.-sig)*ztv/grav
aft=(vkar/log(z1onzt))**2    ! should give .00085 for csiro9

call getqsat(qsat,water(:,1)%temp,atm%ps)
qsat=0.98*qsat ! with Zeng 1998 for sea water
ri=min(grav*atm%zmin*(1.-water(:,1)%temp*srcp/atm%temp)/max(vmag,0.1)**2,rimax)
consea=vmag*charnck/grav
pg%zo=.001    ! first guess
where (ri>0.) ! stable water points 
  fm=1./(1.+bprm*ri)**2
  con=consea*fm*vmag
end where  
do it=1,3
  afroot=vkar/log(atm%zmin/pg%zo)
  af=afroot*afroot
  daf=2.*af*afroot/(vkar*pg%zo)
  where (ri>0.) ! stable water points                                                 
    pg%zo=pg%zo-(pg%zo-con*af)/(1.-con*daf)
  elsewhere     ! unstable water points
    con=cms*2.*bprm*sqrt(-ri*atm%zmin/pg%zo)
    den=1.+af*con
    fm=vmag-vmag*2.*bprm*ri/den
    dfm=2.*bprm*ri*con*(daf-.5*af/pg%zo)/den**2
    pg%zo=pg%zo-(pg%zo-consea*af*fm)/(1.-consea*(daf*fm+af*dfm))
  end where
  pg%zo=min(max(pg%zo,1.5e-5),13.)
enddo    ! it=1,3
afroot=vkar/log(atm%zmin/pg%zo)
af=afroot**2
factch=sqrt(pg%zo*ztv)

where (ri>0.)
  fm=1./(1.+bprm*ri)**2  ! no zo contrib for stable
  fh=fm
elsewhere        ! ri is -ve
  root=sqrt(-ri*atm%zmin/pg%zo)
  den=1.+cms*2.*bprm*af*root
  fm=1.-2.*bprm*ri/den
  den=1.+chs*2.*bprm*factch*aft*root
  fh=1.-2.*bprm*ri/den
end where

pg%wetfac=1.
pg%cd=af*fm
pg%fg=rho*aft*cp*fh*vmag*(water(:,1)%temp-atm%temp/srcp)
pg%eg=rho*aft*lv*fh*vmag*(qsat-atm%qg)
dg2%taux=rho*pg%cd*vmag*atm%u
dg2%tauy=rho*pg%cd*vmag*atm%v

return
end subroutine fluxcalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Estimate saturation mixing ratio (from TAPM)

subroutine getqsat(qsat,temp,ps)

implicit none

real, dimension(wfull), intent(in) :: temp,ps
real, dimension(wfull), intent(out) :: qsat
real, dimension(wfull) :: esatf

where (temp.ge.273.15)
  esatf = 610.*exp(lv/rvap*(1./273.15-1./min(temp,343.)))
elsewhere
  esatf = 610.*exp(ls/rvap*(1./273.15-1./max(temp,123.)))
endwhere
qsat = 0.622*esatf/(ps-0.378*esatf)

return
end subroutine getqsat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Pack sea ice for calcuation

! Neglect leads for now until MLO supports horizontal advection.
! Sea ice model is then only a multi-layer thermodynamical model.

subroutine mloice(dt,atm,diag)

implicit none

integer, intent(in) :: diag
integer nice,iqc,iqw,ii
integer, dimension(wfull) :: icegrid
real, intent(in) :: dt
logical, dimension(wfull) :: cice
type(tatm), dimension(wfull), intent(in) :: atm
type(tatm), dimension(wfull) :: atm_pack
type(tice), dimension(wfull) :: ice_pack
type(tdiag2), dimension(wfull) :: dg2

call mlonewice(diag)
call iceflux(dg2,atm,diag)

! identify sea ice for packing
cice=ice%dic.gt.0.
nice=count(cice)
iqc=0
do iqw=1,wfull
  if (cice(iqw)) then
    iqc=iqc+1
    icegrid(iqc)=iqw
  end if
end do

! pack ice data
ice_pack(1:nice)%dic=ice(icegrid(1:nice))%dic
ice_pack(1:nice)%fracice=ice(icegrid(1:nice))%fracice
do ii=0,3
  ice_pack(1:nice)%tn(ii)=ice(icegrid(1:nice))%tn(ii)
end do

! pack atmosphere data
atm_pack(1:nice)%vnratio=atm(icegrid(1:nice))%vnratio
atm_pack(1:nice)%fbvis=atm(icegrid(1:nice))%fbvis
atm_pack(1:nice)%fbnir=atm(icegrid(1:nice))%fbnir
atm_pack(1:nice)%sg=atm(icegrid(1:nice))%sg
atm_pack(1:nice)%rg=atm(icegrid(1:nice))%rg
atm_pack(1:nice)%rnd=atm(icegrid(1:nice))%rnd

call seaicecalc(nice,dt,ice_pack(1:nice),atm_pack(1:nice))

! unpack ice data
ice(icegrid(1:nice))%dic=ice_pack(1:nice)%dic
ice(icegrid(1:nice))%fracice=ice_pack(1:nice)%fracice
do ii=0,3
  ice(icegrid(1:nice))%tn(ii)=ice_pack(1:nice)%tn(ii)
end do

return
end subroutine mloice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Form seaice before flux calculations

subroutine mlonewice(diag)

implicit none

integer, intent(in) :: diag
logical, dimension(wfull) :: cice
!real, parameter :: hcap=
!real, parameter :: qice=

! formation
cice=water(:,1)%temp.lt.tfi !.and.ice%fracice.le.0.
!newdic=0.
!newtn=tfi
!where(cice)
!  newdic=(tfi-water(:,1)%temp)*hcap*dz(:,1)/qice
!endwhere
!where(newdic.gt.0.15)
!  newdic=0.15
!  newtn=water(:,1)%temp+0.15*qice/(hcap*dz(:,1))
!endwhere
where (cice) ! form new sea-ice
!  ice%dic=newdic
!  ice%fracice=1.
!  ice%tn(0)=273.05
!  ice%tn(1)=newtn
!  ice%tn(2)=newtn
!  ice%tn(3)=newtn
!  ice%sto=0.
!  ice%pl=0. ! set leads to zero until dynamics is working
  water(:,1)%temp=tfi
endwhere

!fracice(wgrid)=ice%fracice
!sicedep(wgrid)=ice%dic

return
end subroutine mlonewice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update sea ice prognostic variables

subroutine seaicecalc(nice,dt,ice_pack,atm_pack)

implicit none

integer, intent(in) :: nice
real, intent(in) :: dt
type(tice), dimension(nice), intent(inout) :: ice_pack
type(tatm), dimension(nice), intent(in) :: atm_pack

!alb=     atm_pack%vnratio*(pg_pack%icevisdiralb*atm_pack%fbvis &
!                          +pg_pack%icevisdifalb*(1.-atm_pack%fbvis))+ &
!    (1.-atm_pack%vnratio)*(pg_pack%icevisdifalb*atm_pack%fbvis &
!                          +pg_pack%icevisdifalb*(1.-atm_pack%fbvis))

!rgn=ice_pack%rg-sbconst*ice_pack%tn(0)**4


! calculate water temperature at ice bottom


! update seaice temperature


! update snow


! update salinity due to melting ice (split form)


return
end subroutine seaicecalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Determine ice fluxes

subroutine iceflux(dg2,atm,diag)

implicit none

integer, intent(in) :: diag
integer it
type(tatm), dimension(wfull), intent(in) :: atm
type(tdiag2), dimension(wfull), intent(inout) :: dg2
real, dimension(wfull) :: qsat,ri,vmag,rho,srcp
real, dimension(wfull) :: fm,fh,af,aft,logzo
real, dimension(wfull) :: den,sig,factch,root

vmag=sqrt(atm%u*atm%u+atm%v*atm%v)
sig=atm%p1/atm%ps
srcp=sig**(rdry/cp)
rho=atm%ps/(rdry*ice%tn(0))

pg%zoice=0.001
logzo=log(atm%zmin/pg%zoice)
af=(vkar/logzo)**2
aft=vkar**2/(logzo*(2.+logzo))
factch=sqrt(7.4)

call getqsat(qsat,ice%tn(0),atm%ps)
ri=min(grav*atm%zmin*(1.-ice%tn(0)*srcp/atm%temp)/max(vmag,0.1)**2,rimax)

where (ri>0.)
  fm=1./(1.+bprm*ri)**2  ! no zo contrib for stable
  fh=fm
elsewhere        ! ri is -ve
  root=sqrt(-ri*atm%zmin/pg%zoice)
  den=1.+cms*2.*bprm*af*root
  fm=1.-2.*bprm*ri/den
  den=1.+chs*2.*bprm*factch*aft*root
  fh=1.-2.*bprm*ri/den
end where

pg%wetfacice=1.+.008*min(ice%tn(0)-273.16,0.)
pg%cdice=af*fm
pg%fgice=rho*aft*cp*fh*vmag*(ice%tn(0)-atm%temp/srcp)
pg%egice=pg%wetfacice*rho*aft*lv*fh*vmag*(qsat-atm%qg)

return
end subroutine iceflux

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Determine screen diagnostics

subroutine scrncalc(atm,diag)

implicit none

integer, intent(in) :: diag
type(tatm), dimension(wfull), intent(in) :: atm
real, dimension(wfull) :: tscrn,qgscrn,uscrn,u10
real, dimension(wfull) :: smixr,qsat

! water
call getqsat(qsat,water(:,1)%temp,atm%ps)
smixr=pg%wetfac*0.98*qsat+(1.-pg%wetfac)*min(0.98*qsat,atm%qg)
call scrntile(tscrn,qgscrn,uscrn,u10,pg%zo,water(:,1)%temp,smixr,atm,diag)
pg%tscrn=tscrn
pg%qgscrn=qgscrn
pg%uscrn=uscrn
pg%u10=u10

! ice
call getqsat(qsat,ice%tn(0),atm%ps)
smixr=pg%wetfac*qsat+(1.-pg%wetfac)*min(qsat,atm%qg)
call scrntile(tscrn,qgscrn,uscrn,u10,pg%zoice,ice%tn(0),smixr,atm,diag)
pg%tscrn=(1.-ice%fracice)*pg%tscrn+ice%fracice*tscrn
pg%qgscrn=(1.-ice%fracice)*pg%qgscrn+ice%fracice*qgscrn
pg%uscrn=(1.-ice%fracice)*pg%uscrn+ice%fracice*uscrn
pg%u10=(1.-ice%fracice)*pg%u10+ice%fracice*u10

return
end subroutine scrncalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! screen diagnostic for individual tile

subroutine scrntile(tscrn,qgscrn,uscrn,u10,zo,stemp,smixr,atm,diag)
      
implicit none

integer, intent(in) :: diag
integer ic
type(tatm), dimension(wfull), intent(in) :: atm
real, dimension(wfull), intent(in) :: zo,stemp,smixr
real, dimension(wfull), intent(out) :: tscrn,qgscrn,uscrn,u10
real, dimension(wfull) :: umag,sig
real, dimension(wfull) :: lzom,lzoh,af,aft,ri,fm,fh,root
real, dimension(wfull) :: denma,denha,cd,thetav,sthetav
real, dimension(wfull) :: thetavstar,z_on_l,z0_on_l,zt_on_l
real, dimension(wfull) :: pm0,ph0,pm1,ph1,integralm,integralh
real, dimension(wfull) :: ustar,qstar,z10_on_l
real, dimension(wfull) :: neutral,neutral10,pm10
real, dimension(wfull) :: integralm10,tstar,scrp
integer, parameter ::  nc     = 5
real, parameter    ::  a_1    = 1.
real, parameter    ::  b_1    = 2./3.
real, parameter    ::  c_1    = 5.
real, parameter    ::  d_1    = 0.35
!real, parameter    ::  aa1    = 3.8
!real, parameter    ::  bb1    = 0.5
!real, parameter    ::  cc1    = 0.3
real, parameter    ::  lna    = 2.3
real, parameter    ::  z0     = 1.5
real, parameter    ::  z10    = 10.

umag=sqrt(atm%u*atm%u+atm%v*atm%v)
sig=atm%p1/atm%ps
scrp=sig**(rdry/cp)
thetav=atm%temp*(1.+0.61*atm%qg)/scrp
sthetav=stemp*(1.+0.61*smixr)

! Roughness length for heat
lzom=log(atm%zmin/zo)
lzoh=lna+lzom

! use Louis as first guess for Dyer and Hicks scheme
af=vkar*vkar/(lzom*lzom)
aft=vkar*vkar/(lzom*lzoh)
! umag is now constrained to be above umin
ri=min(grav*atm%zmin*(1.-sthetav/thetav)/umag**2,rimax)
where (ri>0.)
  fm=1./(1.+bprm*ri)**2
  fh=fm
elsewhere
  root=sqrt(-ri*exp(lzom))
  denma=1.+cms*2.*bprm*af*root
  denha=1.+chs*2.*bprm*aft*exp(0.5*lna)*root
  fm=1.-2.*bprm*ri/denma
  fh=1.-2.*bprm*ri/denha
end where
cd=af*fm

! Dyer and Hicks approach 
thetavstar=aft*fh*(thetav-sthetav)/sqrt(cd)
do ic=1,nc
  z_on_l=vkar*atm%zmin*grav*thetavstar/(thetav*cd*umag**2)
  z_on_l=min(z_on_l,10.)
  z0_on_l  = z_on_l*exp(-lzom)
  zt_on_l  = z0_on_l*exp(-lna)
  where (z_on_l.lt.0.)
    pm0     = (1.-16.*z0_on_l)**(-0.25)
    ph0     = (1.-16.*zt_on_l)**(-0.5)
    pm1     = (1.-16.*z_on_l)**(-0.25)
    ph1     = (1.-16.*z_on_l)**(-0.5)
    integralm = lzom-2.*log((1.+1./pm1)/(1.+1./pm0)) &
                -log((1.+1./pm1**2)/(1.+1./pm0**2)) &
                +2.*(atan(1./pm1)-atan(1./pm0))
    integralh = lzoh-2.*log((1.+1./ph1)/(1.+1./ph0))
  elsewhere
    !--------------Beljaars and Holtslag (1991) momentum & heat            
    pm0 = -(a_1*z0_on_l+b_1*(z0_on_l-(c_1/d_1))*exp(-d_1*z0_on_l) &
           +b_1*c_1/d_1)
    pm1 = -(a_1*z_on_l+b_1*(z_on_l-(c_1/d_1))*exp(-d_1*z_on_l) &
           +b_1*c_1/d_1)
    ph0 = -((1.+(2./3.)*a_1*zt_on_l)**1.5 &
           +b_1*(zt_on_l-(c_1/d_1))*exp(-d_1*zt_on_l) &
           +b_1*c_1/d_1-1.)
    ph1 = -((1.+(2./3.)*a_1*z_on_l)**1.5 &
           +b_1*(z_on_l-(c_1/d_1))*exp(-d_1*z_on_l) &
           +b_1*c_1/d_1-1.)
    integralm = lzom-(pm1-pm0)
    integralh = lzoh-(ph1-ph0)
  endwhere
! where (z_on_l.le.0.4)
  cd = (max(0.01,min(vkar*umag/integralm,2.))/umag)**2
! elsewhere
!   cd = (max(0.01,min(vkar*umag/(aa1*( ( z_on_l**bb1)*
!  &         (1.0+cc1* z_on_l**(1.-bb1))
!  &         -(z0_on_l**bb1)*(1.+cc1*z0_on_l**(1.-bb1)) )),2.))
!  &         /umag)**2
! endwhere
  thetavstar= vkar*(thetav-sthetav)/integralh
end do
ustar=sqrt(cd)*umag
tstar=vkar*(atm%temp/scrp-stemp)/integralh
qstar=vkar*(atm%qg-smixr)/integralh
      
! estimate screen diagnostics
z0_on_l=z0*z_on_l/atm%zmin
z10_on_l=z10*z_on_l/atm%zmin
z0_on_l=min(z0_on_l,10.)
z10_on_l=min(z10_on_l,10.)
neutral=log(atm%zmin/z0)
neutral10=log(atm%zmin/z10)
where (z_on_l.lt.0.)
  ph0     = (1.-16.*z0_on_l)**(-0.50)
  ph1     = (1.-16.*z_on_l)**(-0.50)
  pm0     = (1.-16.*z0_on_l)**(-0.25)
  pm10    = (1.-16.*z10_on_l)**(-0.25)
  pm1     = (1.-16.*z_on_l)**(-0.25)
  integralh = neutral-2.*log((1.+1./ph1)/(1.+1./ph0))
  integralm = neutral-2.*log((1.+1./pm1)/(1.+1./pm0)) &
                 -log((1.+1./pm1**2)/(1.+1./pm0**2)) &
                 +2.*(atan(1./pm1)-atan(1./pm0))
  integralm10 = neutral10-2.*log((1.+1./pm1)/(1.+1./pm10)) &
                -log((1.+1./pm1**2)/(1.+1./pm10**2)) &
                +2.*(atan(1./pm1)-atan(1./pm10))     
elsewhere
  !-------Beljaars and Holtslag (1991) heat function
  ph0  = -((1.+(2./3.)*a_1*z0_on_l)**1.5 &
         +b_1*(z0_on_l-(c_1/d_1)) &
         *exp(-d_1*z0_on_l)+b_1*c_1/d_1-1.)
  ph1  = -((1.+(2./3.)*a_1*z_on_l)**1.5 &
         +b_1*(z_on_l-(c_1/d_1)) &
         *exp(-d_1*z_on_l)+b_1*c_1/d_1-1.)
  pm0 = -(a_1*z0_on_l+b_1*(z0_on_l-(c_1/d_1))*exp(-d_1*z0_on_l) &
         +b_1*c_1/d_1)
  pm10 = -(a_1*z10_on_l+b_1*(z10_on_l-(c_1/d_1)) &
         *exp(-d_1*z10_on_l)+b_1*c_1/d_1)
  pm1  = -(a_1*z_on_l+b_1*(z_on_l-(c_1/d_1))*exp(-d_1*z_on_l) &
         +b_1*c_1/d_1)
  integralh = neutral-(ph1-ph0)
  integralm = neutral-(pm1-pm0)
  integralm10 = neutral10-(pm1-pm10)
endwhere
tscrn       = atm%temp-tstar/vkar*integralh*scrp
qgscrn      = atm%qg-qstar/vkar*integralh
qgscrn      = max(qgscrn,1.E-4)

uscrn=max(umag-ustar/vkar*integralm,0.)
u10=max(umag-ustar/vkar*integralm10,0.)

return
end subroutine scrntile

end module mlo
