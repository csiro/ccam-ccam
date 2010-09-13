
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

module mlo

implicit none

private
public mloinit,mloend,mloeval,mloimport,mloexport,mloload,mlosave,mloregrid,mlodiag,mloalb2,mloalb4, &
       mloscrnout,mloexpice,wlev,depth,mlodwn,mlootf,micdwn,micotf,ocndwn,ocnotf,sssb

! parameters
integer, parameter :: wlev = 20
integer, parameter :: iqx = 500
! model arrays
integer, save :: wfull,ifull,iqwt
logical, dimension(:), allocatable, save :: wpack
real, dimension(:,:), allocatable, save :: depth,dz
real, dimension(:,:), allocatable, save :: depth_hl
real, dimension(:,:), allocatable, save :: dz_hl
real, dimension(:,:,:), allocatable, save :: mlodwn,mlootf ! These variables are for CCAM onthefly.f
real, dimension(:,:), allocatable, save :: micdwn,micotf   ! These variables are for CCAM onthefly.f
real, dimension(:), allocatable, save :: ocndwn,ocnotf     ! These variables are for CCAM onthefly.f
real, dimension(:), allocatable, save :: sssb              ! To hold sea-surface salinity for nudging
real, dimension(:,:), allocatable, save :: w_temp,w_sal,w_u,w_v
real, dimension(:,:), allocatable, save :: i_tn
real, dimension(:), allocatable, save :: i_dic,i_dsn,i_fracice,i_tsurf,i_sto
real, dimension(:), allocatable, save :: p_mixdepth,p_bf
real, dimension(:), allocatable, save :: p_watervisdiralb,p_watervisdifalb,p_waternirdiralb,p_waternirdifalb
real, dimension(:), allocatable, save :: p_icevisdiralb,p_icevisdifalb,p_icenirdiralb,p_icenirdifalb
real, dimension(:), allocatable, save :: p_zo,p_cd,p_fg,p_eg
real, dimension(:), allocatable, save :: p_zoice,p_cdice,p_fgice,p_egice,p_wetfacice
real, dimension(:), allocatable, save :: p_tscrn,p_uscrn,p_qgscrn,p_u10
integer, dimension(:), allocatable, save :: p_mixind
  
! mode
integer, parameter :: incradbf  = 1 ! include shortwave in buoyancy forcing
integer, parameter :: incradgam = 0 ! include shortwave in non-local term
integer, parameter :: salrelax  = 0 ! relax salinity to 34.72 PSU (used for single column mode)
integer, parameter :: zomode    = 2 ! roughness calculation (0=Charnock (CSIRO9), 1=Charnock (zot=zom), 2=Beljaars)
! max depth
real, parameter :: mxd    = 977.6   ! Max depth (m)
real, parameter :: mindep = 1.      ! Thickness of first layer (m)
! model parameters
real, parameter :: ric     = 0.3    ! Critical Ri for diagnosing mixed layer depth
real, parameter :: epsilon = 0.1
! radiation parameters
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
real, parameter :: cp=1004.64             ! Specific heat of dry air at const P
real, parameter :: rdry=287.04            ! Specific gas const for dry air
real, parameter :: rvap=461.5             ! Gas constant for water vapor
! ice parameters
real, parameter :: himin=0.1              ! minimum ice thickness (m)
real, parameter :: rhosn=330.             ! density snow
real, parameter :: rhowt=1025.            ! density water (replace with dg3%rho ?)
real, parameter :: rhoic=900.             ! density ice
real, parameter :: qice=lf*rhoic          ! latent heat of fusion (J m^-3)
real, parameter :: qsnow=qice
real, parameter :: cpi=1.8837e6           ! Sp heat ice  (J/m**3/K)
real, parameter :: cps=6.9069e5           ! Sp heat snow (J/m**3/K)
real, parameter :: condsnw=0.30976        ! conductivity snow
real, parameter :: condice=2.03439        ! conductivity ice
real, parameter :: antarck=0.15e-4
! stability function parameters
real, parameter :: bprm=5.                ! 4.7 in rams
real, parameter :: chs=2.6                ! 5.3 in rams
real, parameter :: cms=5.                 ! 7.4 in rams
real, parameter :: fmroot=0.57735
real, parameter :: rimax=(1./fmroot-1.)/bprm

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Initialise MLO

subroutine mloinit(ifin,depin,diag)

implicit none

integer, intent(in) :: ifin,diag
integer iq,iqw,ii
real, dimension(ifin), intent(in) :: depin
real, dimension(ifin) :: deptmp
real smxd,smnd

if (diag.ge.1) write(6,*) "Initialising MLO"

ifull=ifin
allocate(wpack(ifull))
wpack=depin.gt.mindep
wfull=count(wpack)
if (wfull.eq.0) then
  deallocate(wpack)
  return
end if

allocate(w_temp(wfull,wlev),w_sal(wfull,wlev))
allocate(w_u(wfull,wlev),w_v(wfull,wlev))
allocate(i_tn(wfull,0:2),i_dic(wfull),i_dsn(wfull))
allocate(i_fracice(wfull),i_tsurf(wfull),i_sto(wfull))
allocate(p_mixdepth(wfull),p_bf(wfull))
allocate(p_watervisdiralb(wfull),p_watervisdifalb(wfull))
allocate(p_waternirdiralb(wfull),p_waternirdifalb(wfull))
allocate(p_icevisdiralb(wfull),p_icevisdifalb(wfull))
allocate(p_icenirdiralb(wfull),p_icenirdifalb(wfull))
allocate(p_zo(wfull),p_cd(wfull),p_fg(wfull),p_eg(wfull))
allocate(p_zoice(wfull),p_cdice(wfull),p_fgice(wfull),p_egice(wfull))
allocate(p_wetfacice(wfull),p_mixind(wfull))
allocate(p_tscrn(wfull),p_uscrn(wfull),p_qgscrn(wfull),p_u10(wfull))
allocate(depth(wfull,wlev),dz(wfull,wlev))
allocate(depth_hl(wfull,wlev+1))
allocate(dz_hl(wfull,2:wlev))

iqw=0
iqwt=0
do iq=1,ifull
  if (depin(iq).gt.mindep) then
    iqw=iqw+1
    if (iq.ge.iqx) then
      iqwt=iqw
      exit
    end if
  end if
end do

w_temp=288.             ! K
w_sal=35.               ! PSU
w_u=0.                  ! m/s
w_v=0.                  ! m/s

i_dic=0.                ! m
i_dsn=0.                ! m
i_fracice=0.            ! %
i_tsurf=271.2           ! K
i_tn=271.2              ! K
i_sto=0.

p_mixdepth=100. ! m
p_bf=0.
p_watervisdiralb=0.06
p_watervisdifalb=0.06
p_waternirdiralb=0.06
p_waternirdifalb=0.06
p_icevisdiralb=0.65
p_icevisdifalb=0.65
p_icenirdiralb=0.65
p_icenirdifalb=0.65
p_zo=0.001
p_cd=0.
p_fg=0.
p_eg=0.
p_zoice=0.001
p_cdice=0.
p_fgice=0.
p_egice=0.
p_wetfacice=1.
p_tscrn=273.2
p_qgscrn=0.
p_uscrn=0.
p_u10=0.
p_mixind=wlev-1

! MLO
!depth = (/   0.5,   1.9,   4.3,   8.5,  15.3,  25.3,  39.3,  57.9,  81.9, 112.1 &
!           149.1, 193.7, 246.5, 308.3, 379.9, 461.9, 555.1, 660.1, 777.7, 908.7 /)
! Mk3.5
!depth = (/   5.0,  15.0,  28.2,  42.0,  59.7,  78.5, 102.1, 127.9, 159.5, 194.6, &
!           237.0, 284.7, 341.7, 406.4, 483.2, 570.9, 674.9, 793.8, 934.1 /)

deptmp(1:wfull)=pack(depin,wpack)

smxd=maxval(deptmp(1:wfull))
smnd=minval(deptmp(1:wfull))

do iqw=1,wfull
  call vgrid(deptmp(iqw),depth(iqw,:),depth_hl(iqw,:))
  if (smxd.eq.deptmp(iqw)) then
    write(6,*) "MLO max depth ",depth(iqw,:)
    smxd=smxd+10.
  end if
  if (smnd.eq.deptmp(iqw)) then
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

deallocate(wpack)
deallocate(w_temp,w_sal,w_u,w_v)
deallocate(i_tn,i_dic,i_dsn)
deallocate(i_fracice,i_tsurf,i_sto)
deallocate(p_mixdepth,p_bf)
deallocate(p_watervisdiralb,p_watervisdifalb)
deallocate(p_waternirdiralb,p_waternirdifalb)
deallocate(p_icevisdiralb,p_icevisdifalb)
deallocate(p_icenirdiralb,p_icenirdifalb)
deallocate(p_zo,p_cd,p_fg,p_eg)
deallocate(p_zoice,p_cdice,p_fgice,p_egice)
deallocate(p_wetfacice,p_mixind)
deallocate(p_tscrn,p_uscrn,p_qgscrn,p_u10)
deallocate(depth,dz,depth_hl,dz_hl)

return
end subroutine mloend

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Load MLO data

subroutine mloload(datain,icein,diag)

implicit none

integer, intent(in) :: diag
integer ii
real, dimension(ifull,wlev,4), intent(in) :: datain
real, dimension(ifull,8), intent(in) :: icein

if (wfull.eq.0) return

do ii=1,wlev
  w_temp(:,ii)=pack(datain(:,ii,1),wpack)
  w_sal(:,ii)=pack(datain(:,ii,2),wpack)
  w_u(:,ii)=pack(datain(:,ii,3),wpack)
  w_v(:,ii)=pack(datain(:,ii,4),wpack)
end do
i_tsurf=pack(icein(:,1),wpack)
i_tn(:,0)=pack(icein(:,2),wpack)
i_tn(:,1)=pack(icein(:,3),wpack)
i_tn(:,2)=pack(icein(:,4),wpack)
i_fracice=pack(icein(:,5),wpack)
i_dic=pack(icein(:,6),wpack)
i_dsn=pack(icein(:,7),wpack)
i_sto=pack(icein(:,8),wpack)

return
end subroutine mloload

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Save MLO data

subroutine mlosave(dataout,depout,iceout,diag)

implicit none

integer, intent(in) :: diag
integer ii
real, dimension(ifull,wlev,4), intent(out) :: dataout
real, dimension(ifull,8), intent(out) :: iceout
real, dimension(ifull), intent(out) :: depout

if (wfull.eq.0) return

do ii=1,wlev
  dataout(:,ii,1)=unpack(w_temp(:,ii),wpack,dataout(:,ii,1))
  dataout(:,ii,2)=unpack(w_sal(:,ii),wpack,dataout(:,ii,2))
  dataout(:,ii,3)=unpack(w_u(:,ii),wpack,dataout(:,ii,3))
  dataout(:,ii,4)=unpack(w_v(:,ii),wpack,dataout(:,ii,4))
end do
iceout(:,1)=unpack(i_tsurf,wpack,iceout(:,1))
iceout(:,2)=unpack(i_tn(:,0),wpack,iceout(:,2))
iceout(:,3)=unpack(i_tn(:,1),wpack,iceout(:,3))
iceout(:,4)=unpack(i_tn(:,2),wpack,iceout(:,4))
iceout(:,5)=unpack(i_fracice,wpack,iceout(:,5))
iceout(:,6)=unpack(i_dic,wpack,iceout(:,6))
iceout(:,7)=unpack(i_dsn,wpack,iceout(:,7))
iceout(:,8)=0.
iceout(:,8)=unpack(i_sto,wpack,iceout(:,8))
depout=0.
depout=unpack(depth_hl(:,wlev+1),wpack,depout)

return
end subroutine mlosave

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Import sst for nudging

subroutine mloimport(mode,sst,ilev,diag)

implicit none

integer, intent(in) :: mode,ilev,diag
real, dimension(ifull), intent(in) :: sst

if (wfull.eq.0) return
select case(mode)
  case(0)
    w_temp(:,ilev)=pack(sst,wpack)
  case(1)
    w_sal(:,ilev)=pack(sst,wpack)
end select

return
end subroutine mloimport

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Export sst for nudging

subroutine mloexport(mode,sst,ilev,diag)

implicit none

integer, intent(in) :: mode,ilev,diag
real, dimension(ifull), intent(inout) :: sst

if (wfull.eq.0) return
select case(mode)
  case(0)
    sst=unpack(w_temp(:,ilev),wpack,sst)
  case(1)
    sst=unpack(w_sal(:,ilev),wpack,sst)
end select

return
end subroutine mloexport

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Export ice temp

subroutine mloexpice(tsn,ilev,diag)

implicit none

integer, intent(in) :: ilev,diag
real, dimension(ifull), intent(inout) :: tsn

if (wfull.eq.0) return
select case(ilev)
  case(1)
    tsn=unpack(i_tsurf,wpack,tsn)
  case(2)
    tsn=unpack(i_tn(:,0),wpack,tsn)
  case(3)
    tsn=unpack(i_tn(:,1),wpack,tsn)
  case(4)
    tsn=unpack(i_tn(:,2),wpack,tsn)
  case DEFAULT
    write(6,*) "ERROR: Invalid level ",ilev
    stop
end select

return
end subroutine mloexpice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Return mixed layer depth

subroutine mlodiag(mld,diag)

implicit none

integer, intent(in) :: diag
real, dimension(ifull), intent(out) :: mld

mld=0.
if (wfull.eq.0) return
mld=unpack(p_mixdepth,wpack,mld)

return
end subroutine mlodiag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate screen diagnostics

subroutine mloscrnout(tscrn,qgscrn,uscrn,u10,diag)

implicit none

integer, intent(in) :: diag
real, dimension(ifull), intent(inout) :: tscrn,qgscrn,uscrn,u10

if (wfull.eq.0) return
tscrn=unpack(p_tscrn,wpack,tscrn)
qgscrn=unpack(p_qgscrn,wpack,qgscrn)
uscrn=unpack(p_uscrn,wpack,uscrn)
u10=unpack(p_u10,wpack,u10)

return
end subroutine mloscrnout

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate albedo data (VIS + NIR)

subroutine mloalb2(istart,ifin,coszro,ovisalb,oniralb,diag)

implicit none

integer, intent(in) :: istart,ifin,diag
integer ifinish,ib,ie
real, dimension(ifin), intent(in) :: coszro
real, dimension(ifin), intent(inout) :: ovisalb,oniralb
real, dimension(wfull) :: watervis,waternir,icevis,icenir
real, dimension(wfull) :: costmp

if (wfull.eq.0) return

ifinish=istart+ifin-1
ib=count(wpack(1:istart-1))+1
ie=count(wpack(istart:ifinish))+ib-1
if (ie.lt.ib) return

costmp(ib:ie)=pack(coszro,wpack(istart:ifinish))

watervis(ib:ie)=.05/(costmp(ib:ie)+0.15)
waternir(ib:ie)=.05/(costmp(ib:ie)+0.15)
! need to factor snow age into albedo
icevis(ib:ie)=0.85
icenir(ib:ie)=0.45
ovisalb=unpack(icevis(ib:ie)*i_fracice(ib:ie)+(1.-i_fracice(ib:ie))*watervis(ib:ie),wpack(istart:ifinish),ovisalb)
oniralb=unpack(icenir(ib:ie)*i_fracice(ib:ie)+(1.-i_fracice(ib:ie))*waternir(ib:ie),wpack(istart:ifinish),oniralb)

p_watervisdiralb(ib:ie)=watervis(ib:ie)
p_watervisdifalb(ib:ie)=watervis(ib:ie)
p_waternirdiralb(ib:ie)=waternir(ib:ie)
p_waternirdifalb(ib:ie)=waternir(ib:ie)
p_icevisdiralb(ib:ie)=icevis(ib:ie)
p_icevisdifalb(ib:ie)=icevis(ib:ie)
p_icenirdiralb(ib:ie)=icenir(ib:ie)
p_icenirdifalb(ib:ie)=icenir(ib:ie)

return
end subroutine mloalb2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate albedo data (VIS-DIR, VIS-DIF, NIR-DIR & NIR-DIF)

subroutine mloalb4(istart,ifin,coszro,ovisdir,ovisdif,onirdir,onirdif,diag)

implicit none

integer, intent(in) :: istart,ifin,diag
integer ifinish,ib,ie
real, dimension(ifin), intent(in) :: coszro
real, dimension(ifin), intent(inout) :: ovisdir,ovisdif,onirdir,onirdif
real, dimension(wfull) :: watervisdiralb,watervisdifalb,waternirdiralb,waternirdifalb
real, dimension(wfull) :: icevisdiralb,icevisdifalb,icenirdiralb,icenirdifalb
real, dimension(wfull) :: costmp

if (wfull.eq.0) return

ifinish=istart+ifin-1
ib=count(wpack(1:istart-1))+1
ie=count(wpack(istart:ifinish))+ib-1
if (ie.lt.ib) return

costmp(ib:ie)=pack(coszro(istart:ifinish),wpack(istart:ifinish))

where (costmp(ib:ie).gt.0.)
  watervisdiralb(ib:ie)=0.026/(costmp(ib:ie)**1.7+0.065)+0.15*(costmp(ib:ie)-0.1)* &
                        (costmp(ib:ie)-0.5)*(costmp(ib:ie)-1.)
elsewhere
  watervisdiralb(ib:ie)=0.3925
end where
watervisdifalb(ib:ie)=0.06
waternirdiralb(ib:ie)=watervisdiralb(ib:ie)
waternirdifalb(ib:ie)=0.06
! need to factor snow age into albedo
icevisdiralb(ib:ie)=0.85
icevisdifalb(ib:ie)=0.85
icenirdiralb(ib:ie)=0.45
icenirdifalb(ib:ie)=0.45
ovisdir=unpack(i_fracice(ib:ie)*icevisdiralb(ib:ie)+(1.-i_fracice(ib:ie))*watervisdiralb(ib:ie),wpack(istart:ifinish),ovisdir)
ovisdif=unpack(i_fracice(ib:ie)*icevisdifalb(ib:ie)+(1.-i_fracice(ib:ie))*watervisdifalb(ib:ie),wpack(istart:ifinish),ovisdif)
onirdir=unpack(i_fracice(ib:ie)*icenirdiralb(ib:ie)+(1.-i_fracice(ib:ie))*waternirdiralb(ib:ie),wpack(istart:ifinish),onirdir)
onirdif=unpack(i_fracice(ib:ie)*icenirdifalb(ib:ie)+(1.-i_fracice(ib:ie))*waternirdifalb(ib:ie),wpack(istart:ifinish),onirdif)

p_watervisdiralb(ib:ie)=watervisdiralb(ib:ie)
p_watervisdifalb(ib:ie)=watervisdifalb(ib:ie)
p_waternirdiralb(ib:ie)=waternirdiralb(ib:ie)
p_waternirdifalb(ib:ie)=waternirdifalb(ib:ie)
p_icevisdiralb(ib:ie)=icevisdiralb(ib:ie)
p_icevisdifalb(ib:ie)=icevisdifalb(ib:ie)
p_icenirdiralb(ib:ie)=icenirdiralb(ib:ie)
p_icenirdifalb(ib:ie)=icenirdifalb(ib:ie)

return
end subroutine mloalb4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Regrid MLO data

subroutine mloregrid(depin,mlodat)

implicit none

integer iq,iqw,ii,pos(1),nn
real, dimension(ifull), intent(in) :: depin
real, dimension(ifull,wlev,4), intent(inout) :: mlodat
real, dimension(wfull,wlev,4) :: newdata,newdatb
real, dimension(wfull) :: deptmp
real, dimension(wlev) :: dpin
real, dimension(wlev+1) :: dp_hlin
real x

if (wfull.eq.0) return
deptmp=pack(depin,wpack)
if (all(abs(deptmp-depth(:,wlev)).lt.0.01)) return

do nn=1,4
  do ii=1,wlev
    newdata(:,ii,nn)=pack(mlodat(:,ii,nn),wpack)
  end do
end do

do iqw=1,wfull
  call vgrid(deptmp(iqw),dpin,dp_hlin)
  do ii=1,wlev
    if (depth(iqw,ii).ge.dpin(wlev)) then
      newdatb(iqw,ii,:)=newdata(iqw,wlev,:)
    else if (depth(iqw,ii).le.dpin(1)) then
      newdatb(iqw,ii,:)=newdata(iqw,1,:)
    else
      pos=maxloc(dpin,dpin.lt.depth(iqw,ii))
      pos(1)=max(1,min(wlev-1,pos(1)))
      x=(depth(iqw,ii)-dpin(pos(1)))/(dpin(pos(1)+1)-dpin(pos(1)))
      x=max(0.,min(1.,x))
      newdatb(iqw,ii,:)=newdata(iqw,pos(1)+1,:)*x+newdata(iqw,pos(1),:)*(1.-x)
    end if
  end do  
end do

do nn=1,4
  do ii=1,wlev
    mlodat(:,ii,nn)=unpack(newdatb(:,ii,nn),wpack,mlodat(:,ii,nn))
  end do
end do

return
end subroutine mloregrid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Pack atmospheric data for MLO eval

subroutine mloeval(sst,zo,cd,fg,eg,wetfac,epot,epan,fracice,siced,snowd, &
                   dt,zmin,zmins,sg,rg,precp,uatm,vatm,temp,qg,ps, &
                   f,visnirratio,fbvis,fbnir,diag)

implicit none

integer, intent(in) :: diag
real, intent(in) :: dt
real, dimension(ifull), intent(in) :: sg,rg,precp,f,uatm,vatm,temp,qg,ps,visnirratio,fbvis,fbnir,zmin,zmins
real, dimension(ifull), intent(inout) :: sst,zo,cd,fg,eg,wetfac,fracice,siced,epot,epan,snowd
real, dimension(wfull) :: workb
real, dimension(wfull) :: a_sg,a_rg,a_rnd,a_snd,a_f,a_vnratio,a_fbvis,a_fbnir,a_u,a_v,a_temp,a_qg,a_ps,a_zmin,a_zmins
real, dimension(wfull,wlev) :: d_rho,d_nsq,d_rad,d_alpha,d_beta
real, dimension(wfull) :: d_b0,d_ustar,d_wu0,d_wv0,d_wt0,d_ws0,d_taux,d_tauy
real, dimension(wfull) :: d_ftop,d_tb,d_fb,d_timelt,d_salflx,d_rs
integer, dimension(wfull) :: d_nk,d_did

if (wfull.eq.0) return

a_sg=pack(sg,wpack)
a_rg=pack(rg,wpack)
where (a_temp.ge.273.16)
  a_rnd=pack(precp,wpack)
  a_snd=0.
elsewhere
  a_rnd=0.
  a_snd=pack(precp,wpack)
end where

!a_f=pack(f,wpack)
a_f=0. ! turn off coriolis terms when no geostrophic term
a_vnratio=pack(visnirratio,wpack)
a_fbvis=pack(fbvis,wpack)
a_fbnir=pack(fbnir,wpack)
a_u=pack(uatm,wpack)
a_v=pack(vatm,wpack)
a_temp=pack(temp,wpack)
a_qg=pack(qg,wpack)
a_ps=pack(ps,wpack)
a_zmin=pack(zmin,wpack)
a_zmins=pack(zmins,wpack)

call fluxcalc(a_u,a_v,a_temp,a_qg,a_ps,a_zmin,a_zmins,d_taux,d_tauy,diag)                          ! ocean fluxes
call getrho(a_sg,a_rg,a_rnd,a_snd,a_vnratio,a_fbvis,a_fbnir,d_rho,d_nsq,d_rad,d_alpha,d_beta,            &
            d_b0,d_ustar,d_wu0,d_wv0,d_wt0,d_ws0,d_taux,d_tauy,d_rs)                               ! boundary conditions
call iceflux(dt,a_sg,a_rg,a_rnd,a_snd,a_vnratio,a_fbvis,a_fbnir,a_u,a_v,a_temp,a_qg,a_ps,a_zmin,a_zmins, &
             d_rho,d_ftop,d_tb,d_fb,d_timelt,d_salflx,d_nk,d_did,diag)                             ! ice fluxes
call mlonewice(d_rho,d_timelt,diag)                                                                ! create new ice
call mloice(dt,a_rnd,a_snd,d_alpha,d_beta,d_b0,d_wt0,d_ws0,d_ftop,d_tb,d_fb,d_timelt,d_salflx,d_rs,d_nk, &
            diag)                                                                                  ! update ice
call mlocalc(dt,a_f,d_rho,d_nsq,d_rad,d_alpha,d_b0,d_ustar,d_wu0,d_wv0,d_wt0,d_ws0,d_fb,d_salflx,d_did,  &
             diag)                                                                                 ! update water
call mlonewice(d_rho,d_timelt,diag)                                                                ! create new ice
call scrncalc(a_u,a_v,a_temp,a_qg,a_ps,a_zmin,a_zmins,diag)                                        ! screen diagnostics

sst=unpack((1.-i_fracice)*w_temp(:,1)+i_fracice*i_tsurf,wpack,sst)
workb=(1.-i_fracice)/log(a_zmin/p_zo)**2+i_fracice/log(a_zmin/p_zoice)**2
zo=unpack(a_zmin*exp(-1./sqrt(workb)),wpack,zo)
cd=unpack((1.-i_fracice)*p_cd+i_fracice*p_cdice,wpack,cd)
fg=unpack((1.-i_fracice)*p_fg+i_fracice*p_fgice,wpack,fg)
eg=unpack((1.-i_fracice)*p_eg+i_fracice*p_egice,wpack,eg)
wetfac=unpack((1.-i_fracice)+i_fracice*p_wetfacice,wpack,wetfac)
epan=unpack(p_eg,wpack,epan)
epot=unpack((1.-i_fracice)*p_eg+i_fracice*p_egice/p_wetfacice,wpack,epot)
fracice=unpack(i_fracice,wpack,fracice)
siced=unpack(i_dic,wpack,siced)
snowd=unpack(i_dsn,wpack,snowd)

return
end subroutine mloeval

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MLO calcs for water (no ice)

subroutine mlocalc(dt,a_f,d_rho,d_nsq,d_rad,d_alpha,d_b0,d_ustar,d_wu0,d_wv0,d_wt0,d_ws0,d_fb,d_salflx,d_did,diag)

implicit none

integer, intent(in) :: diag
integer ii,iqw
real, intent(in) :: dt
real, dimension(wfull,wlev) :: km,ks,gammas
real, dimension(wfull,wlev) :: rhs
double precision, dimension(wfull,2:wlev) :: aa
double precision, dimension(wfull,wlev) :: bb,dd
double precision, dimension(wfull,1:wlev-1) :: cc
real, dimension(wfull) :: xp,xm,dumt0,umag
real, dimension(wfull), intent(in) :: a_f
real, dimension(wfull,wlev), intent(in) :: d_rho,d_nsq,d_rad,d_alpha
real, dimension(wfull), intent(in) :: d_b0,d_ustar,d_wu0,d_wv0,d_wt0,d_ws0
real, dimension(wfull), intent(in) :: d_fb,d_salflx
integer, dimension(wfull), intent(in) :: d_did
real, dimension(wfull) :: newa,newb

call getmixdepth(d_rho,d_nsq,d_rad,d_alpha,d_b0,d_ustar,a_f) ! solve for mixed layer depth
call getstab(km,ks,gammas,d_nsq,d_ustar)                     ! solve for stability functions and non-local term

! TEMPERATURE
! use +ve for gamma terms as depth is down
if (incradgam.gt.0) then
  do iqw=1,wfull
    dumt0(iqw)=d_wt0(iqw)-sum(d_rad(iqw,1:p_mixind(iqw)))
  end do
else
  dumt0=d_wt0
end if
rhs(:,1)=(1./dz_hl(:,2))*(ks(:,2)*gammas(:,2)-ks(:,1)*gammas(:,1))
rhs(:,2:wlev)=0.
do ii=2,wlev-1
  where (ii.lt.p_mixind)
    rhs(:,ii)=(0.5/dz(:,ii))*(ks(:,ii+1)*gammas(:,ii+1)-ks(:,ii-1)*gammas(:,ii-1)) ! non-local
  elsewhere (ii.eq.p_mixind)
    rhs(:,ii)=(1./dz_hl(:,ii))*(-ks(:,ii-1)*gammas(:,ii-1))
  end where
end do
where (wlev.eq.p_mixind)
  rhs(:,wlev)=(1./dz_hl(:,wlev))*(-ks(:,wlev-1)*gammas(:,wlev-1))
end where
cc(:,1)=-dt*0.5*(ks(:,1)+ks(:,2))/(dz_hl(:,2)*dz(:,1))
bb(:,1)=1.-cc(:,1)
! use -ve for BC as depth is down
dd(:,1)=w_temp(:,1)+dt*rhs(:,1)*dumt0+dt*d_rad(:,1)/dz(:,1)-dt*d_wt0/dz(:,1)
do ii=2,wlev-1
  aa(:,ii)=-dt*0.5*(ks(:,ii)+ks(:,ii-1))/(dz_hl(:,ii)*dz(:,ii))
  cc(:,ii)=-dt*0.5*(ks(:,ii+1)+ks(:,ii))/(dz_hl(:,ii+1)*dz(:,ii))
  bb(:,ii)=1.-aa(:,ii)-cc(:,ii)
  dd(:,ii)=w_temp(:,ii)+dt*rhs(:,ii)*dumt0+dt*d_rad(:,ii)/dz(:,ii)
end do
aa(:,wlev)=-dt*0.5*(ks(:,wlev)+ks(:,wlev-1))/(dz_hl(:,wlev)*dz(:,wlev))
bb(:,wlev)=1.-aa(:,wlev)
dd(:,wlev)=w_temp(:,wlev)+dt*rhs(:,wlev)*dumt0+dt*d_rad(:,wlev)/dz(:,wlev)
do iqw=1,wfull
  dd(iqw,d_did(iqw))=dd(iqw,d_did(iqw))-dt*i_fracice(iqw)*d_fb(iqw) &
                     /(d_rho(iqw,d_did(iqw))*cp0*dz(iqw,d_did(iqw)))
end do
call thomas(w_temp,aa,bb,cc,dd)

! SALINITY
do ii=1,wlev
  dd(:,ii)=w_sal(:,ii)+dt*rhs(:,ii)*d_ws0
end do
dd(:,1)=dd(:,1)-dt*d_ws0/dz(:,1)
dd(:,1)=dd(:,1)-dt*i_fracice*d_salflx*w_sal(:,1)/(d_rho(:,1)*dz(:,1))
if (salrelax.eq.1) then ! relax salinity
  where (w_sal.gt.1.E-6)
    dd=dd+dt*(34.72-w_sal)/(3600.*24.*365.25*10.)
  end where
end if
call thomas(w_sal,aa,bb,cc,dd)
w_sal=max(0.,w_sal)

! split U diffusion term
cc(:,1)=-dt*0.5*(km(:,1)+km(:,2))/(dz_hl(:,2)*dz(:,1))
bb(:,1)=1.-cc(:,1)
dd(:,1)=w_u(:,1)-dt*d_wu0/dz(:,1)
do ii=2,wlev-1
  aa(:,ii)=-dt*0.5*(km(:,ii)+km(:,ii-1))/(dz_hl(:,ii)*dz(:,ii))
  cc(:,ii)=-dt*0.5*(km(:,ii+1)+km(:,ii))/(dz_hl(:,ii+1)*dz(:,ii))
  bb(:,ii)=1.-aa(:,ii)-cc(:,ii)
  dd(:,ii)=w_u(:,ii)
end do
aa(:,wlev)=-dt*0.5*(km(:,wlev)+km(:,wlev-1))/(dz_hl(:,wlev)*dz(:,wlev))
bb(:,wlev)=1.-aa(:,wlev)
dd(:,wlev)=w_u(:,wlev)
umag=sqrt(w_u(:,wlev)*w_u(:,wlev)+w_v(:,wlev)*w_v(:,wlev))
where (depth_hl(:,wlev+1).lt.(mxd-0.1)) ! bottom drag
  bb(:,wlev)=bb(:,wlev)+dt*cdbot*umag/dz(:,wlev)
end where
call thomas(w_u,aa,bb,cc,dd)

! split V diffusion term
do ii=1,wlev
  dd(:,ii)=w_v(:,ii)
end do
dd(:,1)=dd(:,1)-dt*d_wv0/dz(:,1)
call thomas(w_v,aa,bb,cc,dd)

! Split U and V coriolis terms
xp=1.+(0.5*dt*a_f)**2
xm=1.-(0.5*dt*a_f)**2
do ii=1,wlev
  newa=(w_u(:,ii)*xm+w_v(:,ii)*dt*a_f)/xp
  newb=(w_v(:,ii)*xm-w_u(:,ii)*dt*a_f)/xp
  w_u(:,ii)=newa
  w_v(:,ii)=newb
end do
w_u=max(-100.,min(w_u,100.)) ! MJT suggestion
w_v=max(-100.,min(w_v,100.)) ! MJT suggestion

return
end subroutine mlocalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! tri-diagonal solver

subroutine thomas(outo,aa,bbi,cc,ddi)

implicit none

integer ii
double precision, dimension(wfull,2:wlev), intent(in) :: aa
double precision, dimension(wfull,wlev), intent(in) :: bbi,ddi
double precision, dimension(wfull,wlev-1), intent(in) :: cc
real, dimension(wfull,wlev), intent(out) :: outo
double precision, dimension(wfull,wlev) :: bb,dd,out
double precision, dimension(wfull) :: n

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
outo=out

return
end subroutine thomas

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve for stability functions
! GFDL use a look-up table to speed-up the code...

subroutine getstab(km,ks,gammas,d_nsq,d_ustar)

implicit none

integer ii,iqw
real, dimension(wfull,wlev), intent(out) :: km,ks,gammas
real, dimension(wfull,wlev) :: num,nus,wm,ws,ri
real, dimension(wfull) :: sigma
real, dimension(wfull) :: a2m,a3m,a2s,a3s
real, dimension(wfull) :: numh,wm1,dnumhdz,dwm1ds,g1m,dg1mds
real, dimension(wfull) :: nush,ws1,dnushdz,dws1ds,g1s,dg1sds
real xp,cg
real, dimension(wfull,wlev), intent(in) :: d_nsq
real, dimension(wfull), intent(in) :: d_ustar
real, parameter :: ri0 = 0.7
real, parameter :: nu0 = 50.E-4
real, parameter :: numw = 1.E-4
real, parameter :: nusw = 0.1E-4

ri(:,1)=d_nsq(:,1)*(2.*dz_hl(:,2))**2 &
           /max((w_u(:,1)-w_u(:,2))**2 &
               +(w_v(:,1)-w_v(:,2))**2,1.E-10)
do ii=2,wlev-1
  ri(:,ii)=d_nsq(:,ii)*(2.*dz(:,ii))**2 &
           /max((w_u(:,ii-1)-w_u(:,ii+1))**2 &
               +(w_v(:,ii-1)-w_v(:,ii+1))**2,1.E-10)
end do
ri(:,wlev)=d_nsq(:,wlev)*(2.*dz_hl(:,wlev))**2 &
           /max((w_u(:,wlev-1)-w_u(:,wlev))**2 &
               +(w_v(:,wlev-1)-w_v(:,wlev))**2,1.E-10)

! diffusion ---------------------------------------------------------
num=0.
do ii=1,wlev
  ! s - shear
  where (ri(:,ii).lt.0.)
    num(:,ii)=nu0
  elsewhere (ri(:,ii).lt.ri0)
    num(:,ii)=nu0*(1.-(ri(:,ii)/ri0)**2)**3
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
  call getwx(wm(:,ii),ws(:,ii),depth(:,ii),p_bf,d_ustar,p_mixdepth)
end do
!--------------------------------------------------------------------

! calculate G profile -----------------------------------------------
! -ve as z is down
do iqw=1,wfull
  xp=(p_mixdepth(iqw)-depth(iqw,p_mixind(iqw)))/(depth(iqw,p_mixind(iqw)+1)-depth(iqw,p_mixind(iqw)))
  xp=max(0.,min(1.,xp))
  numh(iqw)=(1.-xp)*num(iqw,p_mixind(iqw))+xp*num(iqw,p_mixind(iqw)+1)
  wm1(iqw)=(1.-xp)*wm(iqw,p_mixind(iqw))+xp*wm(iqw,p_mixind(iqw)+1)
  dnumhdz(iqw)=-(num(iqw,p_mixind(iqw)+1)-num(iqw,p_mixind(iqw)))/dz_hl(iqw,p_mixind(iqw)+1)
  dwm1ds(iqw)=-p_mixdepth(iqw)*(wm(iqw,p_mixind(iqw)+1)-wm(iqw,p_mixind(iqw)))/dz_hl(iqw,p_mixind(iqw)+1)
  nush(iqw)=(1.-xp)*nus(iqw,p_mixind(iqw))+xp*nus(iqw,p_mixind(iqw)+1)
  ws1(iqw)=(1.-xp)*ws(iqw,p_mixind(iqw))+xp*ws(iqw,p_mixind(iqw)+1)
  dnushdz(iqw)=-(nus(iqw,p_mixind(iqw)+1)-nus(iqw,p_mixind(iqw)))/dz_hl(iqw,p_mixind(iqw)+1)
  dws1ds(iqw)=-p_mixdepth(iqw)*(ws(iqw,p_mixind(iqw)+1)-ws(iqw,p_mixind(iqw)))/dz_hl(iqw,p_mixind(iqw)+1)
end do

wm1=max(wm1,1.E-10)
ws1=max(ws1,1.E-10)
  
g1m=numh/(p_mixdepth*wm1)
dg1mds=-dnumhdz/wm1-numh*dwm1ds/(p_mixdepth*wm1*wm1)
g1s=nush/(p_mixdepth*ws1)
dg1sds=-dnushdz/ws1-nush*dws1ds/(p_mixdepth*ws1*ws1)
  
a2m=-2.+3.*g1m-dg1mds
a3m=1.-2.*g1m+dg1mds
a2s=-2.+3.*g1s-dg1sds
a3s=1.-2.*g1s+dg1sds

!--------------------------------------------------------------------
! combine
km=num
ks=nus
do ii=1,wlev
  where (ii.le.p_mixind)
    sigma=depth(:,ii)/p_mixdepth
    km(:,ii)=p_mixdepth*wm(:,ii)*sigma*(1.+sigma*(a2m+a3m*sigma))
    ks(:,ii)=p_mixdepth*ws(:,ii)*sigma*(1.+sigma*(a2m+a3m*sigma))
  end where
end do

km=max(km,numw) ! MJT suggestion
ks=max(ks,nusw) ! MJT suggestion

! non-local term
! gammas is the same for temp and sal when double-diffusion is not employed
cg=10.*vkar*(98.96*vkar*epsilon)**(1./3.)
gammas=0.
do ii=1,wlev
  where (p_bf.lt.0.) ! unstable
    gammas(:,ii)=cg/max(ws(:,ii)*p_mixdepth,1.E-10)
  end where
end do

return
end subroutine getstab

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate mixing layer depth

subroutine getmixdepth(d_rho,d_nsq,d_rad,d_alpha,d_b0,d_ustar,a_f)

implicit none

integer ii,iqw
real vtc,dvsq,vtsq,xp
real, dimension(wfull,wlev) :: ws,wm
real, dimension(wfull) :: dumbf,l,he
real, dimension(wlev) :: rib
real, dimension(wfull,wlev), intent(in) :: d_rho,d_nsq,d_rad,d_alpha
real, dimension(wfull), intent(in) :: d_b0,d_ustar
real, dimension(wfull), intent(in) :: a_f

vtc=1.8*sqrt(0.2/(98.96*epsilon))/(vkar**2*ric)

if (incradbf.gt.0) then
  do ii=1,wlev
    dumbf=d_b0+grav*sum(d_alpha(:,1:ii)*d_rad(:,1:ii),2)
    call getwx(wm(:,ii),ws(:,ii),depth(:,ii),dumbf,d_ustar,depth(:,ii))
  end do
else
  dumbf=d_b0
  do ii=1,wlev
    call getwx(wm(:,ii),ws(:,ii),depth(:,ii),dumbf,d_ustar,depth(:,ii))
  end do
end if

p_mixind=wlev-1
p_mixdepth=depth(:,wlev)
do iqw=1,wfull
  ! should be averaged over 0 < sigma < epsilon instead of 1st level
  rib(1)=0.
  do ii=2,wlev
    vtsq=depth(iqw,ii)*ws(iqw,ii)*sqrt(abs(d_nsq(iqw,ii)))*vtc
    dvsq=(w_u(iqw,1)-w_u(iqw,ii))**2+(w_v(iqw,1)-w_v(iqw,ii))**2
    rib(ii)=(depth(iqw,ii)-depth(iqw,1))*(1.-d_rho(iqw,1)/d_rho(iqw,ii))/max(dvsq+vtsq,1.E-10)
    if (rib(ii).gt.ric) then
      p_mixind(iqw)=ii-1
      xp=min(max((ric-rib(ii-1))/max(rib(ii)-rib(ii-1),1.E-10),0.),1.)
      p_mixdepth(iqw)=(1.-xp)*depth(iqw,ii-1)+xp*depth(iqw,ii)
      exit
    end if
  end do 
end do

! calculate buoyancy forcing
call getbf(d_rad,d_alpha,d_b0)

! impose limits for stable conditions
where(p_bf.gt.0.)
  l=d_ustar*d_ustar*d_ustar/(vkar*p_bf)
  p_mixdepth=min(p_mixdepth,l)
end where
where(a_f.gt.0.)
  he=0.7*d_ustar/abs(a_f)
  p_mixdepth=min(p_mixdepth,he)
end where
p_mixdepth=max(p_mixdepth,depth(:,1))
p_mixdepth=min(p_mixdepth,depth(:,wlev))

! recalculate index for mixdepth
p_mixind=wlev-1
do iqw=1,wfull
  do ii=2,wlev
    if (depth(iqw,ii).gt.p_mixdepth(iqw)) then
      p_mixind(iqw)=ii-1
      exit
    end if
  end do
end do

! recalculate buoyancy forcing
call getbf(d_rad,d_alpha,d_b0)

return
end subroutine getmixdepth

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Estimate bouyancy forcing

subroutine getbf(d_rad,d_alpha,d_b0)

implicit none

integer iqw
real, dimension(wfull,wlev), intent(in) :: d_rad,d_alpha
real, dimension(wfull), intent(in) :: d_b0

if (incradbf.gt.0) then
  do iqw=1,wfull
    p_bf(iqw)=d_b0(iqw)+grav*sum(d_alpha(iqw,1:p_mixind(iqw))*d_rad(iqw,1:p_mixind(iqw)))
  end do
else
  p_bf=d_b0
end if

return
end subroutine getbf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This calculates the stability functions

subroutine getwx(wm,ws,dep,bf,d_ustar,mixdp)

implicit none

real, dimension(wfull), intent(out) :: wm,ws
real, dimension(wfull), intent(in) :: bf,mixdp,dep
real, dimension(wfull) :: zeta,sig,invl,uuu
real, dimension(wfull), intent(in) :: d_ustar
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
uuu=max(d_ustar**3,1.E-10)
invl=vkar*bf ! invl = ustar**3/L or L=ustar**3/(vkar*bf)
invl=max(invl,-10.*uuu/mixdp) ! MJT suggestion
invl=min(invl,1.*uuu/mixdp)   ! MJT suggestion
zeta=sig*mixdp*invl

where (zeta.gt.0.)
  wm=vkar*d_ustar*uuu/(uuu+5.*zeta)
elsewhere (zeta.gt.zetam*uuu)
  wm=vkar*d_ustar*(1.-16.*zeta/uuu)**(1./4.)
elsewhere
  wm=vkar*(am*uuu-cm*zeta)**(1./3.)
end where

where (zeta.gt.0.)
  ws=vkar*d_ustar*uuu/(uuu+5.*zeta)
elsewhere (zeta.gt.zetas*uuu)
  ws=vkar*d_ustar*(1.-16.*zeta/uuu)**(1./2.)
elsewhere
  ws=vkar*(as*uuu-cs*zeta)**(1./3.)
end where

return
end subroutine getwx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate rho from equation of state
! From GFDL (MOM3)

subroutine getrho(a_sg,a_rg,a_rnd,a_snd,a_vnratio,a_fbvis,a_fbnir,d_rho,d_nsq,d_rad,d_alpha,d_beta, &
                  d_b0,d_ustar,d_wu0,d_wv0,d_wt0,d_ws0,d_taux,d_tauy,d_rs)

implicit none

integer ii,i
integer, parameter :: nits=1 ! iterate for density (nits=1 recommended)
real, dimension(wfull) :: t,s,p1,p2,t2,t3,t4,t5,s2,s3,s32
real, dimension(wfull) :: drho0dt,drho0ds,dskdt,dskds,sk,sks
real, dimension(wfull) :: drhodt,drhods,rs0,rho0,hlra,hlrb
real, dimension(wfull) :: visalb,niralb
real, dimension(wfull,wlev) :: rs
real, parameter :: density = 1035.
real, dimension(wfull,wlev), intent(inout) :: d_rho,d_nsq,d_rad,d_alpha,d_beta
real, dimension(wfull), intent(inout) :: d_b0,d_ustar,d_wu0,d_wv0,d_wt0,d_ws0,d_taux,d_tauy,d_rs
real, dimension(wfull), intent(in) :: a_sg,a_rg,a_rnd,a_snd,a_vnratio,a_fbvis,a_fbnir

d_rho=density

t = max(w_temp(:,1)-273.16,-2.)
s = max(w_sal(:,1),0.)
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
    t = max(w_temp(:,ii)-273.16,-2.)
    s = max(w_sal(:,ii),0.01)
    p1 = grav*depth(:,ii)*d_rho(:,ii)/100000.
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
       
    d_rho(:,ii)=rho0/(1.-p1/sk)
    rs(:,ii)=rs0/(1.-p1/sks) ! sal=0.
  
    drhodt=drho0dt/(1.-p1/sk)-rho0*p1*dskdt/((sk-p1)**2)
    drhods=drho0ds/(1.-p1/sk)-rho0*p1*dskds/((sk-p1)**2)
    d_alpha(:,ii)=-drhodt/d_rho(:,ii) ! MOM convention
    d_beta(:,ii)=drhods/d_rho(:,ii)   ! MOM convention
  end do

end do

! buoyancy frequency
d_nsq(:,1)=-(grav/d_rho(:,1))*(d_rho(:,1)-d_rho(:,2))/dz_hl(:,2)
do ii=2,wlev-1
  d_nsq(:,ii)=-(grav/d_rho(:,ii))*0.5*(d_rho(:,ii-1)-d_rho(:,ii+1))/dz(:,ii)
end do
d_nsq(:,wlev)=-(grav/d_rho(:,wlev))*(d_rho(:,wlev-1)-d_rho(:,wlev))/dz_hl(:,wlev)

! shortwave
! use -ve as depth is down
visalb=p_watervisdiralb*a_fbvis+p_watervisdifalb*(1.-a_fbvis)
niralb=p_waternirdiralb*a_fbnir+p_waternirdifalb*(1.-a_fbnir)
hlrb=0.5*(d_rho(:,1)+d_rho(:,2))
d_rad(:,1)=-a_sg/cp0*(((1.-visalb)*a_vnratio*exp(-depth_hl(:,2)/mu_1)+ &
                       (1.-niralb)*(1.-a_vnratio)*exp(-depth_hl(:,2)/mu_2))/hlrb- &
                      ((1.-visalb)*a_vnratio*exp(-depth_hl(:,1)/mu_1)+ &
                       (1.-niralb)*(1.-a_vnratio)*exp(-depth_hl(:,1)/mu_2))/rho0)
do ii=2,wlev-1 
  hlra=0.5*(d_rho(:,ii-1)+d_rho(:,ii))
  hlrb=0.5*(d_rho(:,ii)+d_rho(:,ii+1))
  d_rad(:,ii)=-a_sg/cp0*(((1.-visalb)*a_vnratio*exp(-depth_hl(:,ii+1)/mu_1)+ &
                          (1.-niralb)*(1.-a_vnratio)*exp(-depth_hl(:,ii+1)/mu_2))/hlrb- &
                         ((1.-visalb)*a_vnratio*exp(-depth_hl(:,ii)/mu_1)+ &
                          (1.-niralb)*(1.-a_vnratio)*exp(-depth_hl(:,ii)/mu_2))/hlra)
end do
hlra=0.5*(d_rho(:,wlev-1)+d_rho(:,wlev))
d_rad(:,wlev)=a_sg/cp0*((1.-visalb)*a_vnratio*exp(-depth_hl(:,wlev)/mu_1)+ &
                        (1.-niralb)*(1.-a_vnratio)*exp(-depth_hl(:,wlev)/mu_2))/hlra ! remainder

! Boundary conditions (use rho at level 1 for consistancy with shortwave radiation)
d_wu0=-d_taux/d_rho(:,1)                                            ! BC
d_wv0=-d_tauy/d_rho(:,1)                                            ! BC
d_wt0=-(-p_fg-p_eg+a_rg-sbconst*w_temp(:,1)**4)/(d_rho(:,1)*cp0)    ! BC
d_wt0=d_wt0-lf*a_snd/(d_rho(:,1)*cp0) ! melting snow
d_ws0=(a_rnd+a_snd-p_eg/lv)*w_sal(:,1)/rs(:,1)                      ! BC

d_rs=rs(:,1) ! save for sea-ice flux

d_ustar=max(sqrt(sqrt(d_wu0*d_wu0+d_wv0*d_wv0)),1.E-10)
d_b0=-grav*(d_alpha(:,1)*d_wt0-d_beta(:,1)*d_ws0)
!d_b0=-grav*(-drho0dt*d_wt0-drho0ds*d_ws0)

return
end subroutine getrho

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate fluxes between MLO and atmosphere (from CCAM sflux.f)

subroutine fluxcalc(a_u,a_v,a_temp,a_qg,a_ps,a_zmin,a_zmins,d_taux,d_tauy,diag)

implicit none

integer, intent(in) :: diag
integer it
real, dimension(wfull), intent(in) :: a_u,a_v,a_temp,a_qg,a_ps,a_zmin,a_zmins
real, dimension(wfull), intent(inout) :: d_taux,d_tauy
real, dimension(wfull) :: qsat,ri,vmag,rho,srcp
real, dimension(wfull) :: fm,fh,con,consea,afroot,af,daf
real, dimension(wfull) :: den,dfm,sig,factch,root
real, dimension(wfull) :: aft,atu,atv,dcs,afq,facqch,fq
real ztv
! momentum flux parameters
real, parameter :: charnck=0.018
real, parameter :: chn10=0.00125
! ..... high wind speed - rough sea
real, parameter :: zcom1 = 1.8e-2    ! Charnock's constant
real, parameter :: zcoh1 = 0.0       ! Beljaars 1994 values
real, parameter :: zcoq1 = 0.0
! ..... low wind speed - smooth sea
real, parameter :: gnu   = 1.5e-5
real, parameter :: zcom2 = 0.11
real, parameter :: zcoh2 = 0.40
real, parameter :: zcoq2 = 0.62

atu=a_u-w_u(:,1)
atv=a_v-w_v(:,1)
vmag=sqrt(atu*atu+atv*atv)
sig=exp(-grav*a_zmins/(rdry*a_temp))
srcp=sig**(rdry/cp)
rho=a_ps/(rdry*w_temp(:,1))

call getqsat(qsat,w_temp(:,1),a_ps)
if (zomode.eq.0) then ! CSIRO9
  qsat=0.98*qsat ! with Zeng 1998 for sea water
end if
ri=min(grav*(a_zmin**2/a_zmins)*(1.-w_temp(:,1)*srcp/a_temp)/max(vmag,0.1)**2,rimax)
select case(zomode)
  case(0,1) ! Charnock
    consea=vmag*charnck/grav
    p_zo=0.001    ! first guess
    do it=1,4
      afroot=vkar/log(a_zmin/p_zo)
      af=afroot*afroot
      daf=2.*af*afroot/(vkar*p_zo)
      where (ri>0.) ! stable water points                                                 
        fm=1./(1.+bprm*ri)**2
        con=consea*fm*vmag
        p_zo=p_zo-(p_zo-con*af)/(1.-con*daf)
      elsewhere     ! unstable water points
        con=cms*2.*bprm*sqrt(-ri*a_zmin/p_zo)
        den=1.+af*con
        fm=vmag-vmag*2.*bprm*ri/den
        dfm=vmag*2.*bprm*ri*(con*daf+af*cms*bprm*ri*a_zmin/(sqrt(-ri*a_zmin/p_zo)*p_zo*p_zo))/(den*den) ! MJT suggestion
        p_zo=p_zo-(p_zo-consea*af*fm)/(1.-consea*(daf*fm+af*dfm))
      end where
      p_zo=min(max(p_zo,1.5e-10),0.1)
    enddo    ! it=1,4
  case(2) ! Beljaars
    p_zo=0.001    ! first guess
    do it=1,4
      afroot=vkar/log(a_zmin/p_zo)
      af=afroot*afroot
      daf=2.*af*afroot/(vkar*p_zo)
      where (ri>0.) ! stable water points                                                 
        fm=1./(1.+bprm*ri)**2
        consea=zcom1*vmag**2*fm*af/grav+zcom2*gnu/max(vmag*sqrt(fm*af),gnu)
        dcs=(zcom1*vmag**2/grav-0.5*zcom2*gnu/(max(vmag*sqrt(fm*af),gnu)*fm*af))*(fm*daf)
      elsewhere     ! unstable water points
        con=cms*2.*bprm*sqrt(-ri*a_zmin/p_zo)
        den=1.+af*con
        fm=1.-2.*bprm*ri/den
        dfm=2.*bprm*ri*(con*daf+af*cms*bprm*ri*a_zmin/(sqrt(-ri*a_zmin/p_zo)*p_zo*p_zo))/(den*den) ! MJT suggestion
        consea=zcom1*vmag**2*af*fm/grav+zcom2*gnu/max(vmag*sqrt(fm*af),gnu)
        dcs=(zcom1*vmag**2/grav-0.5*zcom2*gnu/(max(vmag*sqrt(fm*af),gnu)*fm*af))*(fm*daf+dfm*af)
      end where
      p_zo=p_zo-(p_zo-consea)/(1.-dcs)      
      p_zo=min(max(p_zo,1.5e-6),0.1)
    enddo    ! it=1,4  
end select
afroot=vkar/log(a_zmin/p_zo)
af=afroot**2

select case(zomode)
  case(0) ! Charnock CSIRO9
    ztv=exp(vkar/sqrt(chn10))/10.      ! proper inverse of ztsea
    aft=vkar**2/(log(a_zmins*ztv)*log(a_zmins*ztv))
    afq=aft    
    factch=sqrt(p_zo*ztv)
    facqch=factch
  case(1) ! Charnock zot=zom
    aft=vkar**2/(log(a_zmins/p_zo)*log(a_zmins/p_zo))
    afq=aft
    factch=1.
    facqch=1.
  case(2) ! Beljaars
    aft=max(zcoh1+zcoh2*gnu/max(vmag*sqrt(fm*af),gnu),1.5E-6) ! dummy for ztv
    afq=max(zcoq1+zcoq2*gnu/max(vmag*sqrt(fm*af),gnu),1.5E-6) ! dummy for zqv
    factch=sqrt(p_zo/aft)
    facqch=sqrt(p_zo/afq)
    aft=vkar**2/(log(a_zmins/p_zo)*log(a_zmins/aft)) ! now true aft
    afq=vkar**2/(log(a_zmins/p_zo)*log(a_zmins/afq)) ! now true afq
end select

where (ri>0.)
  fm=1./(1.+bprm*ri)**2  ! no zo contrib for stable
  fh=fm
  fq=fm
elsewhere        ! ri is -ve
  root=sqrt(-ri*a_zmin/p_zo)
  den=1.+cms*2.*bprm*af*root
  fm=1.-2.*bprm*ri/den
  root=sqrt(-ri*a_zmins/p_zo)
  den=1.+chs*2.*bprm*factch*aft*root
  fh=1.-2.*bprm*ri/den
  den=1.+chs*2.*bprm*facqch*afq*root
  fq=1.-2.*bprm*ri/den  
end where

p_cd=af*fm
p_fg=rho*aft*cp*fh*vmag*(w_temp(:,1)-a_temp/srcp)
p_eg=rho*afq*lv*fq*vmag*(qsat-a_qg)
d_taux=rho*p_cd*vmag*atu
d_tauy=rho*p_cd*vmag*atv

return
end subroutine fluxcalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Estimate saturation mixing ratio (from TAPM)

subroutine getqsat(qsat,temp,a_ps)

implicit none

real, dimension(wfull), intent(in) :: temp
real, dimension(wfull), intent(out) :: qsat
real, dimension(wfull) :: esatf
real, dimension(wfull), intent(in) :: a_ps

where (temp.ge.273.15)
  esatf = 610.*exp(lv/rvap*(1./273.15-1./min(temp,343.)))
elsewhere
  esatf = 610.*exp(ls/rvap*(1./273.15-1./max(temp,123.)))
endwhere
qsat = 0.622*esatf/(a_ps-0.378*esatf)

return
end subroutine getqsat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Pack sea ice for calcuation

! Neglect leads for now until MLO supports horizontal advection.
! Sea ice model is then only a multi-layer thermodynamical model.

subroutine mloice(dt,a_rnd,a_snd,d_alpha,d_beta,d_b0,d_wt0,d_ws0,d_ftop,d_tb,d_fb,d_timelt,d_salflx,d_rs,d_nk,diag)

implicit none

integer, intent(in) :: diag
integer nice,ii
real, intent(in) :: dt
logical, dimension(wfull) :: cice
real, dimension(wfull), intent(in) :: a_rnd,a_snd
real, dimension(wfull) :: ap_rnd,ap_snd
real, dimension(wfull,0:2) :: ip_tn
real, dimension(wfull) :: ip_dic,ip_dsn,ip_fracice,ip_tsurf,ip_sto
real, dimension(wfull) :: pp_egice
real, dimension(wfull,wlev), intent(in) :: d_alpha,d_beta
real, dimension(wfull), intent(inout) :: d_b0,d_wt0,d_ws0,d_ftop,d_tb,d_fb,d_timelt,d_salflx,d_rs
integer, dimension(wfull), intent(inout) :: d_nk
real, dimension(wfull) :: dp_ftop,dp_tb,dp_fb,dp_timelt,dp_salflx
integer, dimension(wfull) :: dp_nk
! identify sea ice for packing
cice=i_dic.gt.1.E-6
nice=count(cice)

if (nice.eq.0) return

! pack data
do ii=0,2
  ip_tn(1:nice,ii)=pack(i_tn(:,ii),cice)
end do
ip_dic(1:nice)=pack(i_dic,cice)
ip_dsn(1:nice)=pack(i_dsn,cice)
ip_fracice(1:nice)=pack(i_fracice,cice)
ip_tsurf(1:nice)=pack(i_tsurf,cice)
ip_sto(1:nice)=pack(i_sto,cice)
ap_rnd(1:nice)=pack(a_rnd,cice)
ap_snd(1:nice)=pack(a_snd,cice)
pp_egice(1:nice)=pack(p_egice,cice)
dp_ftop(1:nice)=pack(d_ftop,cice)
dp_tb(1:nice)=pack(d_tb,cice)
dp_fb(1:nice)=pack(d_fb,cice)
dp_timelt(1:nice)=pack(d_timelt,cice)
dp_salflx(1:nice)=pack(d_salflx,cice)
dp_nk(1:nice)=pack(d_nk,cice)
! update ice prognostic variables
call seaicecalc(nice,dt,ip_tn(1:nice,:),ip_dic(1:nice),ip_dsn(1:nice),ip_fracice(1:nice),        &
                ip_tsurf(1:nice),ip_sto(1:nice),ap_rnd(1:nice),ap_rnd(1:nice),pp_egice(1:nice),  &
                dp_ftop(1:nice),dp_tb(1:nice),dp_fb(1:nice),dp_timelt(1:nice),dp_salflx(1:nice), &
                dp_nk(1:nice),diag)
! unpack ice data
do ii=0,2
  i_tn(:,ii)=unpack(ip_tn(1:nice,ii),cice,i_tn(:,ii))
end do
i_dic=unpack(ip_dic(1:nice),cice,i_dic)
i_dsn=unpack(ip_dsn(1:nice),cice,i_dsn)
i_fracice=unpack(ip_fracice(1:nice),cice,i_fracice)
i_tsurf=unpack(ip_tsurf(1:nice),cice,i_tsurf)
i_sto=unpack(ip_sto(1:nice),cice,i_sto)
d_ftop=unpack(dp_ftop(1:nice),cice,d_ftop)
d_tb=unpack(dp_tb(1:nice),cice,d_tb)
d_fb=unpack(dp_fb(1:nice),cice,d_fb)
d_timelt=unpack(dp_timelt(1:nice),cice,d_timelt)
d_salflx=unpack(dp_salflx(1:nice),cice,d_salflx)
d_nk=unpack(dp_nk(1:nice),cice,d_nk)

! update ice fraction
where (i_dic.lt.1.E-6)
  i_fracice=0.
  i_dic=0.
  i_dsn=0.
end where
d_ws0=d_ws0-i_fracice*d_salflx*w_sal(:,1)/d_rs ! actually salflx*(watersal-icesal)/density(icesal)
d_b0=-grav*(d_alpha(:,1)*d_wt0-d_beta(:,1)*d_ws0)

return
end subroutine mloice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Form seaice before flux calculations

subroutine mlonewice(d_rho,d_timelt,diag)

implicit none

integer, intent(in) :: diag
real, dimension(wfull) :: newdic,newtn,oldfrac
real, dimension(wfull,wlev), intent(in) :: d_rho
real, dimension(wfull), intent(in) :: d_timelt

! formation
newdic=max((d_timelt-w_temp(:,1))*cp0*d_rho(:,1)*dz(:,1)/qice,0.)
newtn=d_timelt
newdic=min(0.15,newdic)
newtn=w_temp(:,1)+newdic*qice/(cp0*d_rho(:,1)*dz(:,1))
oldfrac=i_fracice
where (newdic.gt.1.E-6) ! form new sea-ice
  i_dic=newdic*(1.-oldfrac)+i_dic*oldfrac
  i_dsn=i_dsn*oldfrac
  i_fracice=(1.-oldfrac)+i_fracice*oldfrac
  i_tsurf=273.05*(1.-oldfrac)+i_tsurf*oldfrac
  i_tn(:,0)=273.05*(1.-oldfrac)+i_tn(:,0)*oldfrac
  i_tn(:,1)=newtn*(1.-oldfrac)+i_tn(:,1)*oldfrac
  i_tn(:,2)=newtn*(1.-oldfrac)+i_tn(:,2)*oldfrac
  i_sto=i_sto*oldfrac
  w_temp(:,1)=d_timelt
endwhere

return
end subroutine mlonewice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update sea ice prognostic variables

subroutine seaicecalc(nice,dt,ip_tn,ip_dic,ip_dsn,ip_fracice,ip_tsurf,ip_sto,ap_rnd,ap_snd, &
                      pp_egice,dp_ftop,dp_tb,dp_fb,dp_timelt,dp_salflx,dp_nk,diag)

implicit none

integer, intent(in) :: nice,diag
integer iqi,nc0,nc1,nc2,nc3,ii
logical, dimension(nice) :: pq0,pq1,pq2,pq3
real, intent(in) :: dt
real, dimension(nice) :: xxx,excess
real, dimension(nice,0:2), intent(inout) :: ip_tn
real, dimension(nice), intent(inout) :: ip_dic,ip_dsn,ip_fracice,ip_tsurf,ip_sto
real, dimension(nice,0:2) :: it_tn
real, dimension(nice) :: it_dic,it_dsn,it_fracice,it_tsurf,it_sto
real, dimension(nice), intent(inout) :: dp_ftop,dp_tb,dp_fb,dp_timelt,dp_salflx
integer, dimension(nice), intent(inout) :: dp_nk
real, dimension(nice) :: dt_ftop,dt_tb,dt_fb,dt_timelt,dt_salflx
integer, dimension(nice) :: dt_nk
real, dimension(nice), intent(in) :: ap_rnd,ap_snd
real, dimension(nice), intent(in) :: pp_egice
real, dimension(nice) :: pt_egice

! snow fall
where ((ap_rnd+ap_snd).gt.0..and.ip_dsn.le.1.E-6)
  ip_tn(:,0)=ip_tn(:,1)
end where
ip_dsn=ip_dsn+dt*(ap_rnd+ap_snd)/1000.

! Pack different ice configurations
nc0=0 ! snow + 1-2 ice layers
pq0=.false.
nc1=0 ! no snow + 1-2 ice layers
pq1=.false.
nc2=0 ! snow + 1 ice layer
pq2=.false.
nc3=0 ! no snow + 1 ice layer
pq3=.false.
do iqi=1,nice
  if (dp_nk(iqi).gt.0.and.ip_dsn(iqi).gt.1.E-6) then
    nc0=nc0+1
    pq0(nc0)=.true.
  elseif (dp_nk(iqi).gt.0) then
    nc1=nc1+1
    pq1(nc1)=.true.
  elseif (ip_dsn(iqi).gt.1.E-6) then
    nc2=nc2+1
    pq2(nc2)=.true.
  else
    nc3=nc3+1
    pq3(nc3)=.true.
  end if
end do

! Update ice and snow temperatures
if (nc0.gt.0) then
  do ii=0,2
    it_tn(1:nc0,ii)=pack(ip_tn(:,ii),pq0)
  end do
  it_dic(1:nc0)=pack(ip_dic,pq0)
  it_dsn(1:nc0)=pack(ip_dsn,pq0)
  it_fracice(1:nc0)=pack(ip_fracice,pq0)
  it_tsurf(1:nc0)=pack(ip_tsurf,pq0)
  it_sto(1:nc0)=pack(ip_sto,pq0)
  pt_egice(1:nc0)=pack(pp_egice,pq0)
  dt_ftop(1:nc0)=pack(dp_ftop,pq0)
  dt_tb(1:nc0)=pack(dp_tb,pq0)
  dt_fb(1:nc0)=pack(dp_fb,pq0)
  dt_timelt(1:nc0)=pack(dp_timelt,pq0)
  dt_salflx(1:nc0)=pack(dp_salflx,pq0)
  dt_nk(1:nc0)=pack(dp_nk,pq0)
  call icetemps2(nc0,dt,it_tn(1:nc0,:),it_dic(1:nc0),it_dsn(1:nc0),it_fracice(1:nc0),it_tsurf(1:nc0),it_sto(1:nc0), &
                 dt_ftop(1:nc0),dt_tb(1:nc0),dt_fb(1:nc0),dt_timelt(1:nc0),dt_salflx(1:nc0),dt_nk(1:nc0),pt_egice(1:nc0),diag)
  do ii=0,2
    ip_tn(:,ii)=unpack(it_tn(1:nc0,ii),pq0,ip_tn(:,ii))
  end do
  ip_dic=unpack(it_dic(1:nc0),pq0,ip_dic)
  ip_dsn=unpack(it_dsn(1:nc0),pq0,ip_dsn)
  ip_fracice=unpack(it_fracice(1:nc0),pq0,ip_fracice)
  ip_tsurf=unpack(it_tsurf(1:nc0),pq0,ip_tsurf)
  ip_sto=unpack(it_sto(1:nc0),pq0,ip_sto)
  dp_ftop=unpack(dt_ftop(1:nc0),pq0,dp_ftop)
  dp_tb=unpack(dt_tb(1:nc0),pq0,dp_tb)
  dp_fb=unpack(dt_fb(1:nc0),pq0,dp_fb)
  dp_timelt=unpack(dt_timelt(1:nc0),pq0,dp_timelt)
  dp_salflx=unpack(dt_salflx(1:nc0),pq0,dp_salflx)
  dp_nk=unpack(dt_nk(1:nc0),pq0,dp_nk)
end if
if (nc1.gt.0) then
  do ii=0,2
    it_tn(1:nc1,ii)=pack(ip_tn(:,ii),pq1)
  end do
  it_dic(1:nc1)=pack(ip_dic,pq1)
  it_dsn(1:nc1)=pack(ip_dsn,pq1)
  it_fracice(1:nc1)=pack(ip_fracice,pq1)
  it_tsurf(1:nc1)=pack(ip_tsurf,pq1)
  it_sto(1:nc1)=pack(ip_sto,pq1)
  pt_egice(1:nc1)=pack(pp_egice,pq1)
  dt_ftop(1:nc1)=pack(dp_ftop,pq1)
  dt_tb(1:nc1)=pack(dp_tb,pq1)
  dt_fb(1:nc1)=pack(dp_fb,pq1)
  dt_timelt(1:nc1)=pack(dp_timelt,pq1)
  dt_salflx(1:nc1)=pack(dp_salflx,pq1)
  dt_nk(1:nc1)=pack(dp_nk,pq1)
  call icetempn2(nc1,dt,it_tn(1:nc1,:),it_dic(1:nc1),it_dsn(1:nc1),it_fracice(1:nc1),it_tsurf(1:nc1),it_sto(1:nc1), &
                 dt_ftop(1:nc1),dt_tb(1:nc1),dt_fb(1:nc1),dt_timelt(1:nc1),dt_salflx(1:nc1),dt_nk(1:nc1),pt_egice(1:nc1),diag)
  do ii=0,2
    ip_tn(:,ii)=unpack(it_tn(1:nc1,ii),pq1,ip_tn(:,ii))
  end do
  ip_dic=unpack(it_dic(1:nc1),pq1,ip_dic)
  ip_dsn=unpack(it_dsn(1:nc1),pq1,ip_dsn)
  ip_fracice=unpack(it_fracice(1:nc1),pq1,ip_fracice)
  ip_tsurf=unpack(it_tsurf(1:nc1),pq1,ip_tsurf)
  ip_sto=unpack(it_sto(1:nc1),pq1,ip_sto)
  dp_ftop=unpack(dt_ftop(1:nc1),pq1,dp_ftop)
  dp_tb=unpack(dt_tb(1:nc1),pq1,dp_tb)
  dp_fb=unpack(dt_fb(1:nc1),pq1,dp_fb)
  dp_timelt=unpack(dt_timelt(1:nc1),pq1,dp_timelt)
  dp_salflx=unpack(dt_salflx(1:nc1),pq1,dp_salflx)
  dp_nk=unpack(dt_nk(1:nc1),pq1,dp_nk)
end if
if (nc2.gt.0) then
  do ii=0,2
    it_tn(1:nc2,ii)=pack(ip_tn(:,ii),pq2)
  end do
  it_dic(1:nc2)=pack(ip_dic,pq2)
  it_dsn(1:nc2)=pack(ip_dsn,pq2)
  it_fracice(1:nc2)=pack(ip_fracice,pq2)
  it_tsurf(1:nc2)=pack(ip_tsurf,pq2)
  it_sto(1:nc2)=pack(ip_sto,pq2)
  pt_egice(1:nc2)=pack(pp_egice,pq2)
  dt_ftop(1:nc2)=pack(dp_ftop,pq2)
  dt_tb(1:nc2)=pack(dp_tb,pq2)
  dt_fb(1:nc2)=pack(dp_fb,pq2)
  dt_timelt(1:nc2)=pack(dp_timelt,pq2)
  dt_salflx(1:nc2)=pack(dp_salflx,pq2)
  dt_nk(1:nc2)=pack(dp_nk,pq2)
  call icetemps1(nc2,dt,it_tn(1:nc2,:),it_dic(1:nc2),it_dsn(1:nc2),it_fracice(1:nc2),it_tsurf(1:nc2),it_sto(1:nc2), &
                 dt_ftop(1:nc2),dt_tb(1:nc2),dt_fb(1:nc2),dt_timelt(1:nc2),dt_salflx(1:nc2),dt_nk(1:nc2),pt_egice(1:nc2),diag)
  do ii=0,2
    ip_tn(:,ii)=unpack(it_tn(1:nc2,ii),pq2,ip_tn(:,ii))
  end do
  ip_dic=unpack(it_dic(1:nc2),pq2,ip_dic)
  ip_dsn=unpack(it_dsn(1:nc2),pq2,ip_dsn)
  ip_fracice=unpack(it_fracice(1:nc2),pq2,ip_fracice)
  ip_tsurf=unpack(it_tsurf(1:nc2),pq2,ip_tsurf)
  ip_sto=unpack(it_sto(1:nc2),pq2,ip_sto)
  dp_ftop=unpack(dt_ftop(1:nc2),pq2,dp_ftop)
  dp_tb=unpack(dt_tb(1:nc2),pq2,dp_tb)
  dp_fb=unpack(dt_fb(1:nc2),pq2,dp_fb)
  dp_timelt=unpack(dt_timelt(1:nc2),pq2,dp_timelt)
  dp_salflx=unpack(dt_salflx(1:nc2),pq2,dp_salflx)
  dp_nk=unpack(dt_nk(1:nc2),pq2,dp_nk)
end if
if (nc3.gt.0) then
  do ii=0,2
    it_tn(1:nc3,ii)=pack(ip_tn(:,ii),pq3)
  end do
  it_dic(1:nc3)=pack(ip_dic,pq3)
  it_dsn(1:nc3)=pack(ip_dsn,pq3)
  it_fracice(1:nc3)=pack(ip_fracice,pq3)
  it_tsurf(1:nc3)=pack(ip_tsurf,pq3)
  it_sto(1:nc3)=pack(ip_sto,pq3)
  pt_egice(1:nc3)=pack(pp_egice,pq3)
  dt_ftop(1:nc3)=pack(dp_ftop,pq3)
  dt_tb(1:nc3)=pack(dp_tb,pq3)
  dt_fb(1:nc3)=pack(dp_fb,pq3)
  dt_timelt(1:nc3)=pack(dp_timelt,pq3)
  dt_salflx(1:nc3)=pack(dp_salflx,pq3)
  dt_nk(1:nc3)=pack(dp_nk,pq3)
  call icetempn1(nc3,dt,it_tn(1:nc3,:),it_dic(1:nc3),it_dsn(1:nc3),it_fracice(1:nc3),it_tsurf(1:nc3),it_sto(1:nc3), &
                 dt_ftop(1:nc3),dt_tb(1:nc3),dt_fb(1:nc3),dt_timelt(1:nc3),dt_salflx(1:nc3),dt_nk(1:nc3),pt_egice(1:nc3),diag)
  do ii=0,2
    ip_tn(:,ii)=unpack(it_tn(1:nc3,ii),pq3,ip_tn(:,ii))
  end do
  ip_dic=unpack(it_dic(1:nc3),pq3,ip_dic)
  ip_dsn=unpack(it_dsn(1:nc3),pq3,ip_dsn)
  ip_fracice=unpack(it_fracice(1:nc3),pq3,ip_fracice)
  ip_tsurf=unpack(it_tsurf(1:nc3),pq3,ip_tsurf)
  ip_sto=unpack(it_sto(1:nc3),pq3,ip_sto)
  dp_ftop=unpack(dt_ftop(1:nc3),pq3,dp_ftop)
  dp_tb=unpack(dt_tb(1:nc3),pq3,dp_tb)
  dp_fb=unpack(dt_fb(1:nc3),pq3,dp_fb)
  dp_timelt=unpack(dt_timelt(1:nc3),pq3,dp_timelt)
  dp_salflx=unpack(dt_salflx(1:nc3),pq3,dp_salflx)
  dp_nk=unpack(dt_nk(1:nc3),pq3,dp_nk)
end if

! white ice formation
xxx=ip_dic+ip_dsn-(rhosn*ip_dsn+rhoic*ip_dic)/rhowt
excess=max(ip_dsn-xxx,0.)*0.1                   ! m of water
dp_salflx=dp_salflx-1000.*excess/dt ! m of water
ip_tn(:,1)=(0.5*ip_dic*ip_tn(:,1)+excess*ip_tn(:,0)) &
              /(0.5*ip_dic+excess) ! Assume 2 levels of ice, hence 0.5*hi
ip_dsn=min(xxx,ip_dsn)
ip_dic=ip_dic+excess

! Snow depth limitation and conversion to ice
excess=max(ip_dsn-0.2,0.)*0.1                   ! m water
dp_salflx=dp_salflx-1000.*excess/dt ! m of water
ip_dsn=min(0.2,ip_dsn)
ip_dic=ip_dic+excess

! Ice depth limitation
!ip_dic=min(5.,ip_dic)

return
end subroutine seaicecalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate ice temperature for two layers and snow (from Mk3 seaice.f)

! Solve for snow surface temperature (implicit time scheme)
!
! Solar Sg(in), LWR Rg(out), Heat Flux Fg(out), Evap Eg(out)
!
! ####################### Surface energy flux f ##### and fa(up)##
! ...ts = Snow temp for thin layer near surface................. S
! -------------------------------------------------------------- N
!                                                         fs(up) O
! ...t0 = Mid point temp of snow layer.......................... W
!
! ## ti = Temp at snow/ice interface #############################
!                                                         f0(up) I
! ...t(1) = Mid point temp of first ice layer................... C
!                                                                E
! ########################################################fl(1)###
!                                                                I
! ...t(2) = Mid point temp of second ice layer.................. C
!                                                                E
! ################################################################

subroutine icetemps2(nc,dt,it_tn,it_dic,it_dsn,it_fracice,it_tsurf,it_sto,dt_ftop,dt_tb,dt_fb,dt_timelt,dt_salflx, &
                     dt_nk,pt_egice,diag)

implicit none

integer, intent(in) :: nc,diag
integer ii,iqi
real, intent(in) :: dt
real gamms,qneed,fact,hinew,con0,con1,con2
real htup,htdown
real, dimension(nc) :: ti,con,bot,rhin,fs
real, dimension(nc) :: subl,snmelt,dhs,dhb,hsold
real, dimension(nc) :: fsnmelt
real, dimension(nc,0:2) :: fl
real, dimension(0:3) :: tp
real, dimension(0:2) :: zp,z
real, dimension(nc,0:2), intent(inout) :: it_tn
real, dimension(nc), intent(inout) :: it_dic,it_dsn,it_fracice,it_tsurf,it_sto
real, dimension(nc), intent(inout) :: dt_ftop,dt_tb,dt_fb,dt_timelt,dt_salflx
integer, dimension(nc), intent(inout) :: dt_nk
real, dimension(nc), intent(in) :: pt_egice

! adjust ice thickness by the remaining reservoir in sto
where (dt_nk.eq.1.and.it_sto.gt.0.)
  it_dic=max(0.,it_dic-it_sto/qice)
end where
where (dt_nk.eq.1)
  it_sto=0.
  it_tn(:,2)=dt_tb
end where

! update surface temperature
rhin=real(dt_nk)/it_dic
con=2.*condsnw/it_dsn
! Ammendments if surface layer < 5cms
where (it_dsn.lt.0.05)
  con=1./(it_dsn/condsnw+it_dic/(condice*real(dt_nk)*2.))
  it_tn(:,0)=it_tn(:,1)
end where
fs=con*(it_tn(:,0)-it_tsurf)
gamms=4.857e4 ! density*specific heat*depth (for snow)
bot=gamms/dt+4.*sbconst*it_tsurf**3+con
it_tsurf=it_tsurf+(dt_ftop+fs)/bot
! Snow melt (snmelt >= 0)
snmelt=10.*max(0.,it_tsurf-273.16)*gamms/qsnow ! m of "snow"
it_tsurf=min(it_tsurf,273.16)                  ! melting condition ts=tsmelt
! Ammendments if surface layer < 5cms : Use adjusted snow depth
where (it_dsn.lt.0.05)
  fl(:,0)=fs ! Middle of first ice layer to snow interface                 ! MJT suggestion
elsewhere
  ti=(it_dic*condsnw*it_tn(:,0)+real(dt_nk)*it_dsn*condice*it_tn(:,1)) &
    /(it_dic*condsnw+real(dt_nk)*it_dsn*condice)
  fl(:,0)=2.*rhin*condice*(it_tn(:,1)-ti) ! Middle of first ice layer to snow interface    
end where

! Upward fluxes between various levels below surface
fl(:,1)=rhin*condice*(it_tn(:,2)-it_tn(:,1)) ! Between ice layers 2 and 1

! Surface evap/sublimation (can be >0 or <0)
! Note : dt*eg/hl in Kgm/m**2 => mms of water
subl=0.01*dt*pt_egice/lv   ! m of "snow"
dt_salflx=dt_salflx+pt_egice/lv
! Change the snow thickness
dhs=snmelt+subl
hsold=it_dsn
where(dhs.gt.it_dsn)
  snmelt=max(it_dsn-subl,0.)
  it_dic=it_dic-0.1*(dhs-it_dsn)
  it_dsn=0.
elsewhere
  it_dsn=it_dsn-dhs
end where

! Update the mid-level snow temperature
it_tn(:,0)=it_tn(:,0)+dt*(fl(:,0)-fs)/(max(it_dsn,1.E-6)*cps)
it_tn(:,0)=min(it_tn(:,0),273.16)

! Remove very thin snow
where (it_dsn.gt.0..and.it_dsn.lt.1.E-6)
  fsnmelt=it_dsn*lf/(0.01*dt) ! J/m2/sec ! heat flux to melt snow
  it_tn(:,1)=it_tn(:,1)-dt*fsnmelt*rhin/cpi
  it_dsn=0.
end where

do iqi=1,nc
  zp=0.
  z=0.

  if (dt_nk(iqi).gt.1) then
    ! update temperature in top layer of ice
    tp(1)=it_tn(iqi,1)+dt*(fl(iqi,1)-fl(iqi,0))*rhin(iqi)/cpi
    tp(1)=min(tp(1),dt_timelt(iqi))
    ! use stored heat in brine pockets to keep temperature at -0.1 until heat is used up
    qneed=(dt_timelt(iqi)-tp(1))*cpi/rhin(iqi) ! J/m**2
    if (it_sto(iqi).le.qneed) then
      fact=it_sto(iqi)/qneed
      tp(1)=dt_timelt(iqi)*fact+tp(1)*(1.-fact)
      tp(1)=min(tp(1),dt_timelt(iqi))
      it_sto(iqi)=0.
    else
      tp(1)=dt_timelt(iqi)
      it_sto(iqi)=it_sto(iqi)-qneed
    end if
  end if

  ! determine amount of bottom ablation or accretion
  fl(iqi,dt_nk(iqi))=condice*(dt_tb(iqi)-it_tn(iqi,dt_nk(iqi)))*2.*rhin(iqi)
  dhb(iqi)=dt*(fl(iqi,dt_nk(iqi))-dt_fb(iqi))/qice
  dhb(iqi)=min(0.01,dhb(iqi)) ! limit the ice growth for cases of poor initial conditions
  hinew=it_dic(iqi)+dhb(iqi)
  ! determine depths below surface of interim and final interfaces
  con0=hinew/real(dt_nk(iqi))
  do ii=1,dt_nk(iqi)
    z(ii)=z(ii-1)-con0
  enddo
  con0=it_dic(iqi)/real(dt_nk(iqi))
  zp(1)=-con0

  if (dhb(iqi).le.0.) then
    ! construct final values of temperature for case of bottom ablation
    zp(dt_nk(iqi))=z(dt_nk(iqi))
    tp(dt_nk(iqi))=it_tn(iqi,dt_nk(iqi))+dt/cpi*(fl(iqi,dt_nk(iqi))-fl(iqi,dt_nk(iqi)-1)) &
                   /(1./rhin(iqi)+dhb(iqi))
    tp(dt_nk(iqi))=min(tp(dt_nk(iqi)),dt_timelt(iqi))
    tp(dt_nk(iqi)+1)=dt_tb(iqi)
    do ii=1,dt_nk(iqi)
      if(z(ii).lt.zp(ii).and.zp(ii).lt.z(ii-1))then
        it_tn(iqi,ii)=((z(ii-1)-zp(ii))*tp(ii)+(zp(ii)-z(ii))*tp(ii+1))*real(dt_nk(iqi))/hinew
      elseif(z(ii).lt.zp(ii-1).and.zp(ii-1).lt.z(ii-1))then
        it_tn(iqi,ii)=((z(ii-1)-zp(ii-1))*tp(ii-1)+(zp(ii-1)-z(ii))*tp(ii))*real(dt_nk(iqi))/hinew
      else
        it_tn(iqi,ii)=tp(ii)
      endif
    enddo
  else
    ! construct final values of temperature for case of bottom accretion
    zp(dt_nk(iqi))=zp(dt_nk(iqi)-1)-con0
    tp(dt_nk(iqi))=it_tn(iqi,dt_nk(iqi))+dt*(fl(iqi,dt_nk(iqi))-fl(iqi,dt_nk(iqi)-1))*rhin(iqi)/cpi
    tp(dt_nk(iqi))=min(tp(dt_nk(iqi)),dt_timelt(iqi))
    tp(dt_nk(iqi)+1)=dt_tb(iqi)
    do ii=1,dt_nk(iqi)
      it_tn(iqi,ii)=((z(ii-1)-zp(ii))*tp(ii)+(zp(ii)-z(ii))*tp(ii+1))*real(dt_nk(iqi))/hinew
    enddo
  end if

  it_dic(iqi)=hinew
  ! test whether to change number of layers
  htup=real(dt_nk(iqi)+1)*himin
  htdown=real(dt_nk(iqi))*himin

  if (it_dic(iqi).ge.htup.and.dt_nk(iqi).lt.2) then
    ! code to increase number of layers
    dt_nk(iqi)=dt_nk(iqi)+1
    tp(1)=it_tn(iqi,1)
    tp(dt_nk(iqi))=it_tn(iqi,dt_nk(iqi)-1)
    it_tn(iqi,1:dt_nk(iqi))=tp(1:dt_nk(iqi))
  elseif (it_dic(iqi).le.htdown) then
    ! code to decrease number of layers
    dt_nk(iqi)=dt_nk(iqi)-1
    if (dt_nk(iqi).gt.0) then
      con1=it_dic(iqi)/real(dt_nk(iqi)+1)
      con2=it_dic(iqi)/real(dt_nk(iqi))
      do ii=1,dt_nk(iqi)
        z(ii) =z(ii-1) -con1
        zp(ii)=zp(ii-1)-con2
      enddo
      do ii=1,dt_nk(iqi)
        tp(ii)=((zp(ii-1)-z(ii))*it_tn(iqi,ii)+(z(ii)-zp(ii))*it_tn(iqi,ii+1))/con2
      enddo
      it_tn(iqi,1:dt_nk(iqi))=tp(1:dt_nk(iqi))
    end if
  end if
end do

return
end subroutine icetemps2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate ice temperature for two layers and no snow

! Solve for surface ice temperature (implicit time scheme)
!
! ####################### Surface energy flux f ##### and fa(up)##
! ...ti = Ice temp for thin layer near surface..................
! -------------------------------------------------------------- I
!                                                         f0(up) C
! ...t(1) = Mid point temp of first ice layer................... E
!                                                                 
! ########################################################fl(1)###
!                                                                I
! ...t(2) = Mid point temp of second ice layer.................. C
!                                                                E
! ################################################################

subroutine icetempn2(nc,dt,it_tn,it_dic,it_dsn,it_fracice,it_tsurf,it_sto,dt_ftop,dt_tb,dt_fb,dt_timelt,dt_salflx, &
                     dt_nk,pt_egice,diag)

implicit none

integer, intent(in) :: nc,diag
integer ii,iqi
real, intent(in) :: dt
real gamms,qneed,fact,hinew,con0,con1,con2
real htup,htdown
real, dimension(nc) :: con,rhin,qmax,ftopadj
real, dimension(nc) :: bot,subl,simelt,dhi,dhb
real, dimension(nc,0:2) :: fl
real, dimension(0:3) :: tp
real, dimension(0:2) :: zp,z
real, dimension(nc,0:2), intent(inout) :: it_tn
real, dimension(nc), intent(inout) :: it_dic,it_dsn,it_fracice,it_tsurf,it_sto
real, dimension(nc), intent(inout) :: dt_ftop,dt_tb,dt_fb,dt_timelt,dt_salflx
integer, dimension(nc), intent(inout) :: dt_nk
real, dimension(nc), intent(in) :: pt_egice

! adjust ice thickness by the remaining reservoir in sto
where (dt_nk.eq.1.and.it_sto.gt.0.)
  it_dic=max(0.,it_dic-it_sto/qice)
end where
where (dt_nk.eq.1)
  it_sto=0.
  it_tn(:,2)=dt_tb
end where

! Get total flux into thin surface ice layer by adding flux from below
rhin=real(dt_nk)/it_dic
con=condice*rhin*2.
ftopadj=con*(it_tn(:,1)-it_tsurf)
gamms=3.471e5 ! density*specific heat*depth (for ice)
bot=gamms/dt+4.*sbconst*it_tsurf**3+con
it_tsurf=it_tsurf+(dt_ftop+ftopadj)/bot
! Ice melt (simelt >= 0)
simelt=max(0.,it_tsurf-dt_timelt)*gamms/qice
it_tsurf=min(it_tsurf,dt_timelt) ! melting condition ti=timelt

! Upward fluxes between various levels below surface
fl(:,0)=ftopadj                                      ! Middle of first ice layer to thin surface layer = ftopadj
fl(:,1)=condice*(it_tn(:,2)-it_tn(:,1))*rhin ! Between ice layers 2 and 1

! Surface evap/sublimation (can be >0 or <0)
subl=0.001*dt*pt_egice/lv
dt_salflx=dt_salflx+pt_egice/lv
dhi=-subl-simelt
qmax=qice*0.5*(it_dic-himin)
where (it_sto.gt.qmax)
  dhi=dhi-(it_sto-qmax)/qice
  it_sto=qmax
end where
where (dhi.lt.0..and.dt_nk.lt.2)
  ! if ice is melting, update thickness when single level
  it_dic=max(0.,it_dic+dhi)
  dhi=0.
  rhin=1./it_dic
end where

do iqi=1,nc
  zp=0.
  z=0.
  
  if (dt_nk(iqi).eq.2) then
    tp(1)=it_tn(iqi,1)+dt/cpi*(fl(iqi,1)-fl(iqi,0))/(it_dic(iqi)/real(dt_nk(iqi))+dhi(iqi))
    tp(1)=min(tp(1),dt_timelt(iqi))
    qneed=(dt_timelt(iqi)-tp(1))*cpi/rhin(iqi) ! J/m**2
    if (it_sto(iqi).le.qneed) then
      fact=it_dic(iqi)/qneed
      tp(1)=dt_timelt(iqi)*fact+tp(1)*(1.-fact)
      it_sto(iqi)=0.
    else
      tp(1)=dt_timelt(iqi)
      it_sto(iqi)=it_sto(iqi)-qneed
    end if
  end if

  ! determine amount of bottom ablation or accretion
  fl(iqi,dt_nk(iqi))=condice*(dt_tb(iqi)-it_tn(iqi,dt_nk(iqi)))*2.*rhin(iqi)
  dhb(iqi)=dt*(fl(iqi,dt_nk(iqi))-dt_fb(iqi))/qice
  dhb(iqi)=min(0.01,dhb(iqi)) ! limit for poor initial conditions
  hinew=it_dic(iqi)+dhi(iqi)+dhb(iqi)
  
  ! determine depths below surface of interim and final interfaces
  con0=hinew/real(dt_nk(iqi))
  do ii=1,dt_nk(iqi)
    z(ii)=z(ii-1)-con0
  enddo
  con0=it_dic(iqi)/real(dt_nk(iqi))
  zp(1)=-con0-dhi(iqi)

  if (dhb(iqi).le.0.) then
    ! construct final values of temperature for case of bottom ablation
    zp(dt_nk(iqi))=z(dt_nk(iqi))
    tp(dt_nk(iqi))=it_tn(iqi,dt_nk(iqi))+dt*(fl(iqi,dt_nk(iqi))-fl(iqi,dt_nk(iqi)-1))/cpi &
                   /(it_dic(iqi)/real(dt_nk(iqi))+dhb(iqi))
    tp(dt_nk(iqi))=min(tp(dt_nk(iqi)),dt_timelt(iqi))
    !tp(0)=ti(iqi) ! MJT suggestion
    tp(dt_nk(iqi)+1)=dt_tb(iqi)
    do ii=1,dt_nk(iqi)
      if(z(ii).lt.zp(ii).and.zp(ii).lt.z(ii-1))then
        it_tn(iqi,ii)=((z(ii-1)-zp(ii))*tp(ii)+(zp(ii)-z(ii))*tp(ii+1))*real(dt_nk(iqi))/hinew
      elseif(z(ii).lt.zp(ii-1).and.zp(ii-1).lt.z(ii-1))then
        it_tn(iqi,ii)=((z(ii-1)-zp(ii-1))*tp(ii-1)+(zp(ii-1)-z(ii))*tp(ii))*real(dt_nk(iqi))/hinew
      else
        it_tn(iqi,ii)=tp(ii)
      endif
    enddo
  else
    ! construct final values of temperature for case of bottom accretion
    zp(dt_nk(iqi))=zp(dt_nk(iqi)-1)-con0
    tp(dt_nk(iqi))=it_tn(iqi,dt_nk(iqi))+dt*(fl(iqi,dt_nk(iqi))-fl(iqi,dt_nk(iqi)-1))*rhin(iqi)/cpi
    tp(dt_nk(iqi))=min(tp(dt_nk(iqi)),dt_timelt(iqi))
    tp(dt_nk(iqi)+1)=dt_tb(iqi)
    do ii=1,dt_nk(iqi)
      it_tn(iqi,ii)=((z(ii-1)-zp(ii))*tp(ii)+(zp(ii)-z(ii))*tp(ii+1))*real(dt_nk(iqi))/hinew
    enddo
  end if

  it_dic(iqi)=hinew
  ! test whether to change number of layers
  htup=real(dt_nk(iqi)+1)*himin
  htdown=real(dt_nk(iqi))*himin

  if (it_dic(iqi).ge.htup.and.dt_nk(iqi).lt.2) then
    ! code to increase number of layers
    dt_nk(iqi)=dt_nk(iqi)+1
    tp(1)=it_tn(iqi,1)
    tp(dt_nk(iqi))=it_tn(iqi,dt_nk(iqi)-1)
    do ii=1,dt_nk(iqi)
      it_tn(iqi,ii)=tp(ii)
    enddo
  else if (it_dic(iqi).le.htdown) then
    ! code to decrease number of layers
    dt_nk(iqi)=dt_nk(iqi)-1
    if (dt_nk(iqi).gt.0) then
      con1=it_dic(iqi)/real(dt_nk(iqi)+1)
      con2=it_dic(iqi)/real(dt_nk(iqi))
      do ii=1,dt_nk(iqi)
        z(ii) =z(ii-1) -con1
        zp(ii)=zp(ii-1)-con2
      enddo
      do ii=1,dt_nk(iqi)
        tp(ii)=((zp(ii-1)-z(ii))*it_tn(iqi,ii)+(z(ii)-zp(ii))*it_tn(iqi,ii+1))/con2
      enddo
      do ii=1,dt_nk(iqi)
        it_tn(iqi,ii)=tp(ii)
      enddo
    end if
  end if
end do

return
end subroutine icetempn2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate ice temperature for one layer and snow

! Solve for snow surface temperature (implicit time scheme)
!
! ####################### Surface energy flux f ##### and fa(up)#S
!                                                                N
! ...ts = Snow temp for thin layer near surface................. O
!                                                                W
! --------------------------------------------------------fs(up)
!                                                                I
! ...ti = Temp at snow/ice interface............................ C
!                                                                E
! ## tb ##########################################################

subroutine icetemps1(nc,dt,it_tn,it_dic,it_dsn,it_fracice,it_tsurf,it_sto,dt_ftop,dt_tb,dt_fb,dt_timelt,dt_salflx, &
                     dt_nk,pt_egice,diag)

implicit none

integer, intent(in) :: nc,diag
real, intent(in) :: dt
real, dimension(nc) :: ti,con,bot,f0
real, dimension(nc) :: subl,snmelt,dhs,dhb
real gamms
real, dimension(nc,0:2), intent(inout) :: it_tn
real, dimension(nc), intent(inout) :: it_dic,it_dsn,it_fracice,it_tsurf,it_sto
real, dimension(nc), intent(inout) :: dt_ftop,dt_tb,dt_fb,dt_timelt,dt_salflx
integer, dimension(nc), intent(inout) :: dt_nk
real, dimension(nc), intent(in) :: pt_egice

! Update tsurf and ti based on fluxes from above and below
con=1./(it_dsn/condsnw+it_dic/condice)
f0=con*(dt_tb-it_tsurf) ! flux from below
gamms=4.857e4                  ! density*specific heat*depth (for snow)
bot=gamms/dt+4.*sbconst*it_tsurf**3+con
it_tsurf=it_tsurf+(dt_ftop+f0)/bot
! Snow melt (snmelt >= 0)
snmelt=10.*max(0.,it_tsurf-273.16)*gamms/qsnow
it_tsurf=min(it_tsurf,273.16) ! melting condition ts=tsmelt
ti=(it_dic*condsnw*it_tsurf+it_dsn*condice*dt_tb)/(it_dic*condsnw+it_dsn*condice)

! Surface evap/sublimation (can be >0 or <0)
subl=0.01*dt*pt_egice/lv
dt_salflx=dt_salflx+pt_egice/lv
dhs=snmelt+subl ! Change the snow thickness
where (dhs.gt.it_dsn)
  it_dic=it_dic-0.1*(dhs-it_dsn)
  snmelt=max(it_dsn-subl,0.)
  it_dsn=0.
elsewhere
  it_dsn=it_dsn-dhs
end where
dhb=dt*(f0-dt_fb)/qice       ! Ice melt (simelt >= 0)
dhb=min(0.01,dhb)             ! limit for poor initial conditions
it_dic=it_dic+dhb

where (it_dic.lt.1.E-3)
  it_dic=0.
end where

where (it_dsn.lt.1.E-4) ! Remove very thin snow
  it_dsn=0.
end where

where (it_dic.ge.himin) ! increase ice thickness
  dt_nk=1
  it_sto=0.
end where

it_tn(:,0)=0.5*(it_tsurf+ti)
it_tn(:,1)=0.5*(ti+dt_tb)
it_tn(:,2)=dt_tb

return
end subroutine icetemps1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate ice temperature for one layer and no snow

! Solve for surface ice temperature (implicit time scheme)
!
! ####################### Surface energy flux f ##### and fa(up)#
!                                                                I
! ...ti = Ice temp for thin layer near surface.................. C
!                                                                E
! ## tb ##########################################################

subroutine icetempn1(nc,dt,it_tn,it_dic,it_dsn,it_fracice,it_tsurf,it_sto,dt_ftop,dt_tb,dt_fb,dt_timelt,dt_salflx, &
                     dt_nk,pt_egice,diag)

implicit none

integer, intent(in) :: nc,diag
real, intent(in) :: dt
real, dimension(nc) :: con,bot,subl,simelt
real, dimension(nc) :: dhi,dhb,f0
real gamms
real, dimension(nc,0:2), intent(inout) :: it_tn
real, dimension(nc), intent(inout) :: it_dic,it_dsn,it_fracice,it_tsurf,it_sto
real, dimension(nc), intent(inout) :: dt_ftop,dt_tb,dt_fb,dt_timelt,dt_salflx
integer, dimension(nc), intent(inout) :: dt_nk
real, dimension(nc), intent(in) :: pt_egice

! Update tsurf and ti based on fluxes from above and below
con=condice/it_dic
f0=con*(dt_tb-it_tsurf) ! flux from below
gamms=3.471e5                  ! density*specific heat*depth (for ice)
bot=gamms/dt+4.*sbconst*it_tsurf**3+con
it_tsurf=it_tsurf+(dt_ftop+f0)/bot
simelt=max(0.,it_tsurf-dt_timelt)*gamms/qice ! Ice melt (simelt >= 0)
it_tsurf=min(it_tsurf,dt_timelt)       ! melting condition ti=timelt

! Surface evap/sublimation (can be >0 or <0)
subl=0.001*dt*pt_egice/lv
dt_salflx=dt_salflx+pt_egice/lv
dhi=-subl-simelt
dhb=dt*(f0-dt_fb)/qice           ! determine amount of bottom ablation or accretion
dhb=min(0.01,dhb)                 ! limit for poor initial conditions
it_dic=it_dic+dhi+dhb ! update ice thickness

where (it_dic.lt.1.E-3)
  it_dic=0.
end where

where (it_dic.ge.himin) ! increase ice thickness
  dt_nk=1
  it_sto=0.
end where

it_tn(:,0)=it_tsurf
it_tn(:,1)=0.5*(it_tsurf+dt_tb)
it_tn(:,2)=dt_tb

return
end subroutine icetempn1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Determine ice fluxes

subroutine iceflux(dt,a_sg,a_rg,a_rnd,a_snd,a_vnratio,a_fbvis,a_fbnir,a_u,a_v,a_temp,a_qg,a_ps,a_zmin,a_zmins, &
                   d_rho,d_ftop,d_tb,d_fb,d_timelt,d_salflx,d_nk,d_did,diag)

implicit none

integer, intent(in) :: diag
integer iqw,ii
real, dimension(wfull), intent(in) :: a_sg,a_rg,a_rnd,a_snd,a_vnratio,a_fbvis,a_fbnir,a_u,a_v,a_temp,a_qg,a_ps,a_zmin,a_zmins
real, dimension(wfull,wlev), intent(in) :: d_rho
real, dimension(wfull), intent(inout) :: d_ftop,d_tb,d_fb,d_timelt,d_salflx
integer, dimension(wfull), intent(inout) :: d_nk,d_did
real, intent(in) :: dt
real, dimension(wfull) :: qsat,ri,vmag,rho,srcp
real, dimension(wfull) :: fm,fh,af,aft
real, dimension(wfull) :: den,sig,root
real, dimension(wfull) :: alb,qmax,eye
real x,factch

vmag=sqrt(a_u*a_u+a_v*a_v)
sig=exp(-grav*a_zmins/(rdry*a_temp))
srcp=sig**(rdry/cp)
rho=a_ps/(rdry*i_tsurf)

p_zoice=0.001
af=vkar**2/(log(a_zmin/p_zoice)*log(a_zmin/p_zoice))
aft=vkar**2/(log(a_zmins/p_zoice)*(2.3+log(a_zmins/p_zoice)))
factch=sqrt(7.4)

call getqsat(qsat,i_tsurf,a_ps)
ri=min(grav*(a_zmin**2/a_zmins)*(1.-i_tsurf*srcp/a_temp)/max(vmag,0.1)**2,rimax)

where (ri>0.)
  fm=1./(1.+bprm*ri)**2  ! no zo contrib for stable
  fh=fm
elsewhere        ! ri is -ve
  root=sqrt(-ri*a_zmin/p_zoice)
  den=1.+cms*2.*bprm*af*root
  fm=1.-2.*bprm*ri/den
  root=sqrt(-ri*a_zmins/p_zoice)
  den=1.+chs*2.*bprm*factch*aft*root
  fh=1.-2.*bprm*ri/den
end where

p_wetfacice=max(1.+.008*min(i_tsurf-273.16,0.),0.)
p_cdice=af*fm
p_fgice=rho*aft*cp*fh*vmag*(i_tsurf-a_temp/srcp)
p_egice=p_wetfacice*rho*aft*lv*fh*vmag*(qsat-a_qg)

! ice melting temperature (should depend on salinity)
d_timelt=273.16-0.054*w_sal(:,1) ! from CICE

! determine water temperature at bottom of ice
d_did=0
do iqw=1,wfull
  if (i_dic(iqw).gt.depth(iqw,1)) then
    d_tb(iqw)=w_temp(iqw,wlev)
    d_did(iqw)=wlev
    do ii=2,wlev
      if (i_dic(iqw).lt.depth(iqw,ii)) then
        x=(i_dic(iqw)-depth(iqw,ii-1))/dz_hl(iqw,ii)
        x=max(min(x,1.),0.)
        d_did(iqw)=ii-1+nint(x)
        exit
      end if
    end do
  else
    d_did(iqw)=1
  end if
  ! water temperature at bottom of ice
  d_tb(iqw)=w_temp(iqw,d_did(iqw))
  ! flux from water to ice
  x=sum(dz(iqw,1:d_did(iqw))) ! cumulative depth
  d_fb(iqw)=cp0*d_rho(iqw,d_did(iqw))*antarck*(d_tb(iqw)-d_timelt(iqw))/x
  !   If the actual ocean model top level temperature (sstocean) is less
  !   than freezing, then add an additional ice to ocean heat flux
  !   proportional to the deficit of this temperature below freezing.
  !   Use a heat transfer coefficient hfrz, typically 20-200 times
  !   larger than antarck above. This simulates formation of sea ice from
  !   the ocean as the water temperature drops below freezing, and thereby
  !   prevents the ocean temperature from falling much below freezing
  !   while conserving heat flux.  
!  if (d_tb(iqw).lt.271.3) then
!    d_fb(iqw)=d_fb(iqw)+cp0*d_rho(iqw,d_did(iqw))*6.E-4*(d_tb(iqw)-d_timelt(iqw))/(0.25*dz(iqw,d_did(iqw)))
!  end if
end do

! index of different ice configurations
d_nk=min(int(i_dic/himin),2)

! radiation
alb=     a_vnratio*(p_icevisdiralb*a_fbvis+p_icevisdifalb*(1.-a_fbvis))+ &
    (1.-a_vnratio)*(p_icevisdifalb*a_fbvis+p_icevisdifalb*(1.-a_fbvis))
qmax=qice*0.5*(i_dic-himin)
eye=0.
where (i_dsn.lt.1.E-6.and.i_sto.le.qmax.and.d_nk.gt.0)
  eye=0.35
end where
i_sto=i_sto+dt*eye*a_sg
d_ftop=-p_fgice-p_egice+a_rg-sbconst*i_tsurf**4+a_sg*(1.-alb)*(1.-eye)

! Add flux of heat due to converting any rain to snowfall over ice
d_ftop=d_ftop+lf*a_rnd ! rain (mm) to W/m**2

! Salinity flux due to melting ice (zero to zero for non-ice points)
d_salflx=0.

return
end subroutine iceflux

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Determine screen diagnostics

subroutine scrncalc(a_u,a_v,a_temp,a_qg,a_ps,a_zmin,a_zmins,diag)

implicit none

integer, intent(in) :: diag
real, dimension(wfull), intent(in) :: a_u,a_v,a_temp,a_qg,a_ps,a_zmin,a_zmins
real, dimension(wfull) :: tscrn,qgscrn,uscrn,u10
real, dimension(wfull) :: smixr,qsat

! water
call getqsat(qsat,w_temp(:,1),a_ps)
if (zomode.eq.0) then
  smixr=0.98*qsat
else
  smixr=qsat
end if
call scrntile(tscrn,qgscrn,uscrn,u10,p_zo,w_temp(:,1),smixr,a_u,a_v,a_temp,a_qg,a_zmin,a_zmins,diag)
p_tscrn=tscrn
p_qgscrn=qgscrn
p_uscrn=uscrn
p_u10=u10

! ice
call getqsat(qsat,i_tsurf,a_ps)
smixr=p_wetfacice*qsat+(1.-p_wetfacice)*min(qsat,a_qg)
call scrntile(tscrn,qgscrn,uscrn,u10,p_zoice,i_tsurf,smixr,a_u,a_v,a_temp,a_qg,a_zmin,a_zmins,diag)
p_tscrn=(1.-i_fracice)*p_tscrn+i_fracice*tscrn
p_qgscrn=(1.-i_fracice)*p_qgscrn+i_fracice*qgscrn
p_uscrn=(1.-i_fracice)*p_uscrn+i_fracice*uscrn
p_u10=(1.-i_fracice)*p_u10+i_fracice*u10

return
end subroutine scrncalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! screen diagnostic for individual tile

subroutine scrntile(tscrn,qgscrn,uscrn,u10,zo,stemp,smixr,a_u,a_v,a_temp,a_qg,a_zmin,a_zmins,diag)
      
implicit none

integer, intent(in) :: diag
integer ic
real, dimension(wfull), intent(in) :: a_u,a_v,a_temp,a_qg,a_zmin,a_zmins
real, dimension(wfull), intent(in) :: zo,stemp,smixr
real, dimension(wfull), intent(out) :: tscrn,qgscrn,uscrn,u10
real, dimension(wfull) :: umag,sig
real, dimension(wfull) :: lzom,lzoh,thetav,sthetav
real, dimension(wfull) :: thetavstar,z_on_l,z0_on_l,z0s_on_l,zt_on_l
real, dimension(wfull) :: pm0,ph0,pm1,ph1,integralm,integralh
real, dimension(wfull) :: ustar,qstar,z10_on_l
real, dimension(wfull) :: neutrals,neutral,neutral10,pm10
real, dimension(wfull) :: integralm10,tstar,scrp
integer, parameter ::  nc     = 5
real, parameter    ::  a_1    = 1.
real, parameter    ::  b_1    = 2./3.
real, parameter    ::  c_1    = 5.
real, parameter    ::  d_1    = 0.35
real, parameter    ::  lna    = 2.3
real, parameter    ::  z0     = 1.5
real, parameter    ::  z10    = 10.

umag=sqrt(a_u*a_u+a_v*a_v)
sig=exp(-grav*a_zmins/(rdry*a_temp))
scrp=sig**(rdry/cp)
thetav=a_temp*(1.+0.61*a_qg)/scrp
sthetav=stemp*(1.+0.61*smixr)

! Roughness length for heat
lzom=log(a_zmin/zo)
lzoh=lna+log(a_zmins/zo)

! Dyer and Hicks approach 
thetavstar=vkar*(thetav-sthetav)/lzoh
ustar=vkar*umag/lzom
do ic=1,nc
  z_on_l=vkar*a_zmins*grav*thetavstar/(thetav*ustar**2)
  z_on_l=min(z_on_l,10.)
  z0_on_l  = z_on_l*exp(-lzom)
  zt_on_l  = z_on_l*exp(-lzoh)
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
  thetavstar=vkar*(thetav-sthetav)/integralh
  ustar=vkar*umag/integralm
end do
tstar=vkar*(a_temp-stemp)/integralh
qstar=vkar*(a_qg-smixr)/integralh
      
! estimate screen diagnostics
z0s_on_l=z0*z_on_l/a_zmins
z0_on_l=z0*z_on_l/a_zmin
z10_on_l=z10*z_on_l/a_zmin
z0s_on_l=min(z0s_on_l,10.)
z0_on_l=min(z0_on_l,10.)
z10_on_l=min(z10_on_l,10.)
neutrals=log(a_zmins/z0)
neutral=log(a_zmin/z0)
neutral10=log(a_zmin/z10)
where (z_on_l.lt.0.)
  ph0     = (1.-16.*z0s_on_l)**(-0.50)
  ph1     = (1.-16.*z_on_l)**(-0.50)
  pm0     = (1.-16.*z0_on_l)**(-0.25)
  pm10    = (1.-16.*z10_on_l)**(-0.25)
  pm1     = (1.-16.*z_on_l)**(-0.25)
  integralh = neutrals-2.*log((1.+1./ph1)/(1.+1./ph0))
  integralm = neutral-2.*log((1.+1./pm1)/(1.+1./pm0)) &
                 -log((1.+1./pm1**2)/(1.+1./pm0**2)) &
                 +2.*(atan(1./pm1)-atan(1./pm0))
  integralm10 = neutral10-2.*log((1.+1./pm1)/(1.+1./pm10)) &
                -log((1.+1./pm1**2)/(1.+1./pm10**2)) &
                +2.*(atan(1./pm1)-atan(1./pm10))     
elsewhere
  !-------Beljaars and Holtslag (1991) heat function
  ph0  = -((1.+(2./3.)*a_1*z0s_on_l)**1.5 &
         +b_1*(z0s_on_l-(c_1/d_1)) &
         *exp(-d_1*z0s_on_l)+b_1*c_1/d_1-1.)
  ph1  = -((1.+(2./3.)*a_1*z0s_on_l)**1.5 &
         +b_1*(z0s_on_l-(c_1/d_1)) &
         *exp(-d_1*z0s_on_l)+b_1*c_1/d_1-1.)
  pm0 = -(a_1*z0_on_l+b_1*(z0_on_l-(c_1/d_1))*exp(-d_1*z0_on_l) &
         +b_1*c_1/d_1)
  pm10 = -(a_1*z10_on_l+b_1*(z10_on_l-(c_1/d_1)) &
         *exp(-d_1*z10_on_l)+b_1*c_1/d_1)
  pm1  = -(a_1*z_on_l+b_1*(z_on_l-(c_1/d_1))*exp(-d_1*z_on_l) &
         +b_1*c_1/d_1)
  integralh = neutrals-(ph1-ph0)
  integralm = neutral-(pm1-pm0)
  integralm10 = neutral10-(pm1-pm10)
endwhere
tscrn  = a_temp-tstar*integralh/vkar
qgscrn = a_qg-qstar*integralh/vkar
qgscrn = max(qgscrn,1.E-4)

uscrn=max(umag-ustar*integralm/vkar,0.)
u10  =max(umag-ustar*integralm10/vkar,0.)

return
end subroutine scrntile

end module mlo