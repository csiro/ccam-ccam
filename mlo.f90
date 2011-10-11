
! This is a 1D, mixed layer ocean model based on Large, et al (1994), for ensemble regional climate simulations.
! This code is also used for modelling lakes in ACCESS and CCAM.  In CCAM this module interfaces with
! mlodynamics.f90 for river routing, diffusion and advection routines.

! This version has a relatively thin 1st layer (e.g, 0.5m) so as to better reproduce a diurnal cycle in SST.  It also
! supports a thermodynamic model of sea ice based on O'Farrell's from Mk3.5.  We have included a free surface so that
! lakes can change depth, etc.

! This version can assimilate SSTs from GCMs, using a convolution based digital filter (see nestin.f),
! which avoids problems with complex land-sea boundary conditions.

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
       mloscrnout,mloimpice,mloexpice,mloexpdep,mloexpdensity,wlev,micdwn

! parameters
integer, parameter :: wlev = 20
integer, parameter :: iqx = 4148
! model arrays
integer, save :: wfull,ifull,iqwt
logical, dimension(:), allocatable, save :: wpack
real, dimension(:,:), allocatable, save :: depth,dz,depth_hl,dz_hl
real, dimension(:,:), allocatable, save :: micdwn          ! This variable is for CCAM onthefly.f
real, dimension(:,:), allocatable, save :: w_temp,w_sal,w_u,w_v
real, dimension(:,:), allocatable, save :: i_tn
real, dimension(:), allocatable, save :: w_eta
real, dimension(:), allocatable, save :: i_dic,i_dsn,i_fracice,i_tsurf,i_sto,i_u,i_v
real, dimension(:), allocatable, save :: p_mixdepth,p_bf
real, dimension(:), allocatable, save :: p_watervisdiralb,p_watervisdifalb,p_waternirdiralb,p_waternirdifalb
real, dimension(:), allocatable, save :: p_icevisdiralb,p_icevisdifalb,p_icenirdiralb,p_icenirdifalb
real, dimension(:), allocatable, save :: p_zo,p_zoh,p_zoq,p_cd,p_cds,p_fg,p_eg
real, dimension(:), allocatable, save :: p_zoice,p_zohice,p_zoqice,p_cdice,p_cdsice,p_fgice,p_egice,p_wetfacice
real, dimension(:), allocatable, save :: p_tscrn,p_uscrn,p_qgscrn,p_u10
integer, dimension(:), allocatable, save :: p_mixind
  
! mode
integer, parameter :: incradbf  = 1 ! include shortwave in buoyancy forcing
integer, parameter :: incradgam = 1 ! include shortwave in non-local term
integer, parameter :: zomode    = 2 ! roughness calculation (0=Charnock (CSIRO9), 1=Charnock (zot=zom), 2=Beljaars)
integer, parameter :: mixmeth   = 1 ! Refine mixed layer depth calculation (0=None, 1=Iterative)
integer, parameter :: salrelax  = 0 ! relax salinity to 34.72 PSU (used for single column mode)
integer, parameter :: deprelax  = 0 ! surface height (0=vary, 1=relax, 2=set to zero)
integer, parameter :: buoymeth  = 0 ! buoyancy calculation (0=density, 1=alpha/beta)
! max depth
real, parameter :: mxd    = 914.    ! Max depth (m)
real, parameter :: mindep = 0.5     ! Thickness of first layer (m)
real, parameter :: minwater = 1.    ! Minimum water height above bottom (m)
! model parameters
real, parameter :: ric     = 0.3    ! Critical Ri for diagnosing mixed layer depth
real, parameter :: epsilon = 0.1
! radiation parameters
real, parameter :: mu_1 = 23.       ! VIS depth (m) - Type I
real, parameter :: mu_2 = 0.35      ! NIR depth (m) - Type I
! physical parameters
real, parameter :: vkar=0.4               ! von Karman constant
real, parameter :: lv=2.501e6             ! Latent heat of vaporisation (J kg^-1)
real, parameter :: lf=3.337e5             ! Latent heat of fusion (J kg^-1)
real, parameter :: ls=lv+lf               ! Latent heat of sublimation (J kg^-1)
real, parameter :: grav=9.80              ! graviational constant (m/s^2)
real, parameter :: sbconst=5.67e-8        ! Stefan-Boltzmann constant
real, parameter :: cdbot=2.4E-3           ! bottom drag coefficent
real, parameter :: cp0=3990.              ! heat capacity of mixed layer (J kg^-1 K^-1)
real, parameter :: cp=1004.64             ! Specific heat of dry air at const P
real, parameter :: rdry=287.04            ! Specific gas const for dry air
real, parameter :: rvap=461.5             ! Gas constant for water vapor
! ice parameters
real, parameter :: himin=0.1              ! minimum ice thickness for multiple layers (m)
real, parameter :: icebreak=0.05          ! minimum ice thickness before breakup (1D model)
real, parameter :: fracbreak=0.05         ! minimum ice fraction (1D model)
real, parameter :: icemin=0.01            ! minimum ice thickness (m)
real, parameter :: icemax=4.              ! maximum ice thickness (m)
real, parameter :: rhoic=900.             ! density ice
real, parameter :: rhosn=330.             ! density snow
real, parameter :: rhowt=1025.            ! density water (replace with d_rho ?)
real, parameter :: qice=lf*rhoic          ! latent heat of fusion (J m^-3)
real, parameter :: qsnow=lf*rhosn
real, parameter :: cpi=1.8837e6           ! Sp heat ice  (J/m**3/K)
real, parameter :: cps=6.9069e5           ! Sp heat snow (J/m**3/K)
real, parameter :: condice=2.03439        ! conductivity ice
real, parameter :: condsnw=0.30976        ! conductivity snow
real, parameter :: gammi=3.471e5          ! specific heat*depth (for ice)  (J m^-2 K^-1)
real, parameter :: gamms=4.857e4          ! specific heat*depth (for snow) (J m^-2 K^-1)
!real, parameter :: emisice=0.95          ! emissivity of ice
real, parameter :: emisice=1.             ! emissivity of ice
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
real, dimension(wlev) :: dumdf
real, dimension(wlev+1) :: dumdh
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
allocate(w_eta(wfull))
allocate(i_tn(wfull,0:2),i_dic(wfull),i_dsn(wfull))
allocate(i_fracice(wfull),i_tsurf(wfull),i_sto(wfull))
allocate(i_u(wfull),i_v(wfull))
allocate(p_mixdepth(wfull),p_bf(wfull))
allocate(p_watervisdiralb(wfull),p_watervisdifalb(wfull))
allocate(p_waternirdiralb(wfull),p_waternirdifalb(wfull))
allocate(p_icevisdiralb(wfull),p_icevisdifalb(wfull))
allocate(p_icenirdiralb(wfull),p_icenirdifalb(wfull))
allocate(p_zo(wfull),p_zoh(wfull),p_zoq(wfull),p_cd(wfull),p_cds(wfull),p_fg(wfull),p_eg(wfull))
allocate(p_zoice(wfull),p_zohice(wfull),p_zoqice(wfull),p_cdice(wfull),p_cdsice(wfull),p_fgice(wfull),p_egice(wfull))
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
w_eta=0.                ! m

i_dic=0.                ! m
i_dsn=0.                ! m
i_fracice=0.            ! %
i_tsurf=271.2           ! K
i_tn=271.2              ! K
i_sto=0.
i_u=0.                  ! m/s
i_v=0.                  ! m/s

p_mixdepth=-1. ! m
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
p_zoh=0.001
p_zoq=0.001
p_cd=0.
p_cds=0.
p_fg=0.
p_eg=0.
p_zoice=0.001
p_zohice=0.001
p_zoqice=0.001
p_cdice=0.
p_cdsice=0.
p_fgice=0.
p_egice=0.
p_wetfacice=1.
p_tscrn=273.2
p_qgscrn=0.
p_uscrn=0.
p_u10=0.
p_mixind=wlev-1

! MLO - 20 level (shallow)
!depth = (/  0.25,   1.1,   3.0,   6.7,  12.8,  22.0,  35.1,  52.8,  76.7, 104.5, &
!           140.0, 182.9, 233.8, 293.4, 362.5, 441.8, 531.9, 633.5, 747.4, 874.3 /)
! MLO - 40 level (deep)
!depth = (/   0.5,   1.7,   3.7,   6.8,  11.5,  18.3,  27.7,  40.1,  56.0,  75.9, &
!           100.2, 129.4, 164.0, 204.3, 251.0, 304.4, 365.1, 433.4, 509.9, 595.0, &
!           689.2, 793.0, 906.7,1031.0,1166.2,1312.8,1471.3,1642.2,1825.9,2022.8, &
!          2233.5,2458.4,2698.0,2952.8,3233.1,3509.5,3812.5,4132.4,4469.9,4825.3 /)
! Mk3.5 - 31 level
!depth = (/   5.0,  15.0,  28.2,  42.0,  59.7,  78.5, 102.1, 127.9, 159.5, 194.6, &
!           237.0, 284.7, 341.7, 406.4, 483.2, 570.9, 674.9, 793.8, 934.1,1095.2, &
!          1284.7,1502.9,1758.8,2054.3,2400.0,2800.0,3200.0,3600.0,4000.0,4400.0, &
!          4800.0 /)

deptmp(1:wfull)=pack(depin,wpack)

smxd=maxval(deptmp(1:wfull))
smnd=minval(deptmp(1:wfull))

do iqw=1,wfull
  call vgrid(deptmp(iqw),dumdf,dumdh)
  depth(iqw,:)=dumdf
  depth_hl(iqw,:)=dumdh
  !if (smxd.eq.deptmp(iqw)) then
  !  write(6,*) "MLO max depth ",depth(iqw,:)
  !  smxd=smxd+10.
  !end if
  !if (smnd.eq.deptmp(iqw)) then
  !  write (6,*) "MLO min depth ",depth(iqw,:)
  !  smnd=smnd-10.
  !end if
end do
!depth_hl(:,wlev+1)=deptmp
!depth(:,wlev)=0.5*(depth_hl(:,wlev)+depth_hl(:,wlev+1))
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
!if (dd.gt.x) then
!  al=(mindep*x-dd)/(x-x*x*x)           ! hybrid levels
!  bt=(mindep*x*x*x-dd)/(x*x*x-x)       ! hybrid levels
  al=dd*(mindep*x/mxd-1.)/(x-x*x*x)     ! sigma levels
  bt=dd*(mindep*x*x*x/mxd-1.)/(x*x*x-x) ! sigma levels 
  do ii=1,wlev+1
    x=real(ii-1)
    depth_hlout(ii)=al*x*x*x+bt*x ! ii is for half level ii-0.5
  end do
!else
!  depth_hlout(1)=0.
!  depth_hlout(2)=mindep       ! hybrid levels
!  do ii=3,wlev+1
!    x=(dd-mindep)*real(ii-2)/real(wlev-1)+mindep             ! hybrid levels
!    depth_hlout(ii)=x
!  end do
!end if
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
deallocate(w_eta)
deallocate(i_tn,i_dic,i_dsn)
deallocate(i_fracice,i_tsurf,i_sto)
deallocate(i_u,i_v)
deallocate(p_mixdepth,p_bf)
deallocate(p_watervisdiralb,p_watervisdifalb)
deallocate(p_waternirdiralb,p_waternirdifalb)
deallocate(p_icevisdiralb,p_icevisdifalb)
deallocate(p_icenirdiralb,p_icenirdifalb)
deallocate(p_zo,p_zoh,p_zoq,p_cd,p_cds,p_fg,p_eg)
deallocate(p_zoice,p_zohice,p_zoqice,p_cdice,p_cdsice,p_fgice,p_egice)
deallocate(p_wetfacice,p_mixind)
deallocate(p_tscrn,p_uscrn,p_qgscrn,p_u10)
deallocate(depth,dz,depth_hl,dz_hl)

return
end subroutine mloend

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Load MLO data

subroutine mloload(datain,shin,icein,diag)

implicit none

integer, intent(in) :: diag
integer ii
real, dimension(ifull,wlev,4), intent(in) :: datain
real, dimension(ifull,10), intent(in) :: icein
real, dimension(ifull), intent(in) :: shin

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
i_u=pack(icein(:,9),wpack)
i_v=pack(icein(:,10),wpack)
w_eta=pack(shin,wpack)

return
end subroutine mloload

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Save MLO data

subroutine mlosave(dataout,depout,shout,iceout,diag)

implicit none

integer, intent(in) :: diag
integer ii
real, dimension(ifull,wlev,4), intent(inout) :: dataout
real, dimension(ifull,10), intent(inout) :: iceout
real, dimension(ifull), intent(inout) :: depout,shout

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
iceout(:,9)=unpack(i_u,wpack,iceout(:,9))
iceout(:,10)=unpack(i_v,wpack,iceout(:,10))
depout=0.
depout=unpack(depth_hl(:,wlev+1),wpack,depout)
shout=0.
shout=unpack(w_eta,wpack,shout)

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
  case(2)
    w_u(:,ilev)=pack(sst,wpack)
  case(3)
    w_v(:,ilev)=pack(sst,wpack)
  case(4)
    w_eta=pack(sst,wpack)
end select

return
end subroutine mloimport

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Export ice temp

subroutine mloimpice(tsn,ilev,diag)

implicit none

integer, intent(in) :: ilev,diag
real, dimension(ifull), intent(inout) :: tsn

if (wfull.eq.0) return
select case(ilev)
  case(1)
    i_tsurf=pack(tsn,wpack)
  case(2)
    i_tn(:,0)=pack(tsn,wpack)
  case(3)
    i_tn(:,1)=pack(tsn,wpack)
  case(4)
    i_tn(:,2)=pack(tsn,wpack)
  case(5)
    i_fracice=pack(tsn,wpack)
  case(6)
    i_dic=pack(tsn,wpack)
  case(7)
    i_dsn=pack(tsn,wpack)
  case(8)
    i_sto=pack(tsn,wpack)
  case(9)
    i_u=pack(tsn,wpack)
  case(10)
    i_v=pack(tsn,wpack)
  case DEFAULT
    write(6,*) "ERROR: Invalid mode ",ilev
    stop
end select

return
end subroutine mloimpice

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
  case(2)
    sst=unpack(w_u(:,ilev),wpack,sst)
  case(3)
    sst=unpack(w_v(:,ilev),wpack,sst)
  case(4)
    sst=unpack(w_eta,wpack,sst)
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
  case(5)
    tsn=unpack(i_fracice,wpack,tsn)
  case(6)
    tsn=unpack(i_dic,wpack,tsn)
  case(7)
    tsn=unpack(i_dsn,wpack,tsn)
  case(8)
    tsn=unpack(i_sto,wpack,tsn)
  case(9)
    tsn=unpack(i_u,wpack,tsn)
  case(10)
    tsn=unpack(i_v,wpack,tsn)
  case DEFAULT
    write(6,*) "ERROR: Invalid mode ",ilev
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
real, dimension(wfull) :: costmp,pond,snow

if (wfull.eq.0) return

ifinish=istart+ifin-1
ib=count(wpack(1:istart-1))+1
ie=count(wpack(istart:ifinish))+ib-1
if (ie.lt.ib) return

costmp(ib:ie)=pack(coszro,wpack(istart:ifinish))

!pond(ib:ie)=max(1.+.008*min(i_tsurf(ib:ie)-273.16,0.),0.)
pond(ib:ie)=0.
snow(ib:ie)=min(max(i_dsn(ib:ie)/0.05,0.),1.)

watervis(ib:ie)=.05/(costmp(ib:ie)+0.15)
waternir(ib:ie)=.05/(costmp(ib:ie)+0.15)
! need to factor snow age into albedo
icevis(ib:ie)=(0.85*(1.-pond(ib:ie))+0.75*pond(ib:ie))*(1.-snow(ib:ie))+(0.95*(1.-pond(ib:ie))+0.85*pond(ib:ie))*snow(ib:ie)
icenir(ib:ie)=(0.45*(1.-pond(ib:ie))+0.35*pond(ib:ie))*(1.-snow(ib:ie))+(0.65*(1.-pond(ib:ie))+0.55*pond(ib:ie))*snow(ib:ie)
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
real, dimension(wfull) :: costmp,pond,snow

if (wfull.eq.0) return

ifinish=istart+ifin-1
ib=count(wpack(1:istart-1))+1
ie=count(wpack(istart:ifinish))+ib-1
if (ie.lt.ib) return

costmp(ib:ie)=pack(coszro,wpack(istart:ifinish))

!pond(ib:ie)=max(1.+.008*min(i_tsurf(ib:ie)-273.16,0.),0.)
pond(ib:ie)=0.
snow(ib:ie)=min(max(i_dsn(ib:ie)/0.05,0.),1.)

where (costmp(ib:ie).gt.0.)
  p_watervisdiralb(ib:ie)=0.026/(costmp(ib:ie)**1.7+0.065)+0.15*(costmp(ib:ie)-0.1)* &
                          (costmp(ib:ie)-0.5)*(costmp(ib:ie)-1.)
elsewhere
  p_watervisdiralb(ib:ie)=0.3925
end where
p_watervisdifalb(ib:ie)=0.06
p_waternirdiralb(ib:ie)=p_watervisdiralb(ib:ie)
p_waternirdifalb(ib:ie)=0.06
! need to factor snow age into albedo
p_icevisdiralb(ib:ie)=(0.85*(1.-pond(ib:ie))+0.75*pond(ib:ie))*(1.-snow(ib:ie))+(0.95*(1.-pond(ib:ie))+0.85*pond(ib:ie))*snow(ib:ie)
p_icevisdifalb(ib:ie)=p_icevisdiralb(ib:ie)
p_icenirdiralb(ib:ie)=(0.45*(1.-pond(ib:ie))+0.35*pond(ib:ie))*(1.-snow(ib:ie))+(0.65*(1.-pond(ib:ie))+0.55*pond(ib:ie))*snow(ib:ie)
p_icenirdifalb(ib:ie)=p_icenirdiralb(ib:ie)
ovisdir=unpack(i_fracice(ib:ie)*p_icevisdiralb(ib:ie)+(1.-i_fracice(ib:ie))*p_watervisdiralb(ib:ie),wpack(istart:ifinish),ovisdir)
ovisdif=unpack(i_fracice(ib:ie)*p_icevisdifalb(ib:ie)+(1.-i_fracice(ib:ie))*p_watervisdifalb(ib:ie),wpack(istart:ifinish),ovisdif)
onirdir=unpack(i_fracice(ib:ie)*p_icenirdiralb(ib:ie)+(1.-i_fracice(ib:ie))*p_waternirdiralb(ib:ie),wpack(istart:ifinish),onirdir)
onirdif=unpack(i_fracice(ib:ie)*p_icenirdifalb(ib:ie)+(1.-i_fracice(ib:ie))*p_waternirdifalb(ib:ie),wpack(istart:ifinish),onirdif)

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
! Export ocean depth data

subroutine mloexpdep(mode,odep,ii,diag)

implicit none

integer, intent(in) :: mode,ii,diag
real, dimension(ifull), intent(out) :: odep

odep=0.
select case(mode)
  case(0)
    odep=unpack(depth(:,ii),wpack,odep)
  case(1)
    odep=unpack(dz(:,ii),wpack,odep)
  case default
    write(6,*) "ERROR: Unknown mloexpdep mode ",mode
    stop
end select
  
return
end subroutine mloexpdep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Extract density

subroutine mloexpdensity(odensity,tt,ss,dep,ddz,pxtr,diag)

implicit none

integer, intent(in) :: diag
integer ii
real, dimension(ifull), intent(in) :: pxtr
real, dimension(ifull,wlev), intent(in) :: tt,ss,dep,ddz
real, dimension(ifull,wlev), intent(out) :: odensity
real, dimension(ifull,wlev) :: alpha,beta,rs,rho0

call calcdensity(ifull,odensity,alpha,beta,rs,rho0,tt,ss,dep,ddz,pxtr)

return
end subroutine mloexpdensity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Pack atmospheric data for MLO eval

subroutine mloeval(sst,zo,cd,cds,fg,eg,wetfac,epot,epan,fracice,siced,snowd, &
                   dt,zmin,zmins,sg,rg,precp,uatm,vatm,temp,qg,ps,f,     &
                   visnirratio,fbvis,fbnir,inflow,diag,calcprog)

implicit none

integer, intent(in) :: diag
real, intent(in) :: dt
real, dimension(ifull), intent(in) :: sg,rg,precp,f,uatm,vatm,temp,qg,ps,visnirratio,fbvis,fbnir,inflow,zmin,zmins
real, dimension(ifull), intent(inout) :: sst,zo,cd,cds,fg,eg,wetfac,fracice,siced,epot,epan,snowd
real, dimension(wfull) :: workb
real, dimension(wfull) :: a_sg,a_rg,a_rnd,a_snd,a_f,a_vnratio,a_fbvis,a_fbnir,a_u,a_v,a_temp,a_qg,a_ps,a_zmin,a_zmins
real, dimension(wfull) :: a_inflow
real, dimension(wfull,wlev) :: d_rho,d_rs,d_nsq,d_rad,d_alpha,d_beta
real, dimension(wfull) :: d_b0,d_ustar,d_wu0,d_wv0,d_wt0,d_ws0,d_wm0,d_taux,d_tauy,d_ftop,d_bot,d_tb,d_zcr
real, dimension(wfull) :: d_fb,d_timelt,d_tauxicw,d_tauyicw,d_tauxica,d_tauyica
integer, dimension(wfull) :: d_nk
logical, intent(in) :: calcprog

if (wfull.eq.0) return

a_sg=pack(sg,wpack)
a_rg=pack(rg,wpack)
a_f=pack(f,wpack)
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
a_inflow=pack(inflow,wpack)
where (a_temp.ge.273.16)
  a_rnd=pack(precp,wpack)
  a_snd=0.
elsewhere
  a_rnd=0.
  a_snd=pack(precp,wpack)
end where

d_zcr=max((1.+w_eta/depth_hl(:,wlev+1)),0.01) ! adjust levels for free surface
call getrho(a_ps,d_rho,d_rs,d_alpha,d_beta,d_zcr)                                                    ! equation of state
call fluxcalc(dt,a_u,a_v,a_temp,a_qg,a_ps,a_zmin,a_zmins,d_rho,d_taux,d_tauy,d_zcr,diag)             ! water fluxes
call getwflux(a_sg,a_rg,a_rnd,a_snd,a_vnratio,a_fbvis,a_fbnir,a_inflow,d_rho,d_rs,d_nsq,d_rad,     &
              d_alpha,d_beta,d_b0,d_ustar,d_wu0,d_wv0,d_wt0,d_ws0,d_wm0,d_taux,d_tauy,d_zcr)         ! boundary conditions
d_timelt=273.16-0.054*w_sal(:,1) ! ice melting temperature from CICE
! ice formation is called for both calcprog=.true. and .false.
call mlonewice(dt,d_rho,d_timelt,d_wm0,d_zcr,diag)                                                   ! create new ice
call iceflux(dt,a_sg,a_rg,a_rnd,a_snd,a_vnratio,a_fbvis,a_fbnir,a_u,a_v,a_temp,a_qg,a_ps,a_zmin,   &
             a_zmins,d_rho,d_ftop,d_bot,d_tb,d_fb,d_timelt,d_nk,d_tauxica,d_tauyica,d_tauxicw,     &
             d_tauyicw,d_zcr,diag)                                                                   ! ice fluxes
if (calcprog) then
  call mloice(dt,a_rnd,a_snd,d_alpha,d_beta,d_b0,d_wu0,d_wv0,d_wt0,d_ws0,d_wm0,d_ftop,d_bot,d_tb,  &
              d_fb,d_timelt,d_tauxica,d_tauyica,d_tauxicw,d_tauyicw,d_ustar,d_rho,d_nk,diag)         ! update ice
  call mlocalc(dt,a_f,d_rho,d_nsq,d_rad,d_alpha,d_beta,d_b0,d_ustar,d_wu0,d_wv0,d_wt0,d_ws0,d_wm0, &
               d_zcr,diag)                                                                           ! update water
end if 
call scrncalc(a_u,a_v,a_temp,a_qg,a_ps,a_zmin,a_zmins,diag)                                          ! screen diagnostics

workb=emisice**0.25*i_tsurf
sst=unpack((1.-i_fracice)*w_temp(:,1)+i_fracice*workb,wpack,sst)
workb=(1.-i_fracice)/log(a_zmin/p_zo)**2+i_fracice/log(a_zmin/p_zoice)**2
zo=unpack(a_zmin*exp(-1./sqrt(workb)),wpack,zo)
cd=unpack((1.-i_fracice)*p_cd+i_fracice*p_cdice,wpack,cd)
cds=unpack((1.-i_fracice)*p_cds+i_fracice*p_cdsice,wpack,cds)
fg=unpack((1.-i_fracice)*p_fg+i_fracice*p_fgice,wpack,fg)
eg=unpack((1.-i_fracice)*p_eg+i_fracice*p_egice,wpack,eg)
wetfac=unpack((1.-i_fracice)+i_fracice*p_wetfacice,wpack,wetfac)
epan=unpack(p_eg,wpack,epan)
epot=unpack((1.-i_fracice)*p_eg+i_fracice*p_egice/p_wetfacice,wpack,epot)
fracice=0.
fracice=unpack(i_fracice,wpack,fracice)
siced=0.
siced=unpack(i_dic,wpack,siced)
snowd=unpack(i_dsn,wpack,snowd)

return
end subroutine mloeval

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MLO calcs for water (no ice)

subroutine mlocalc(dt,a_f,d_rho,d_nsq,d_rad,d_alpha,d_beta,d_b0,d_ustar,d_wu0,d_wv0,d_wt0,d_ws0,d_wm0, &
                   d_zcr,diag)

implicit none

integer, intent(in) :: diag
integer ii,iqw
real, intent(in) :: dt
real, dimension(wfull,wlev) :: km,ks,gammas
real, dimension(wfull,wlev) :: rhs
double precision, dimension(wfull,2:wlev) :: aa
double precision, dimension(wfull,wlev) :: bb,dd
double precision, dimension(wfull,1:wlev-1) :: cc
real, dimension(wfull,wlev), intent(in) :: d_rho,d_nsq,d_rad,d_alpha,d_beta
real, dimension(wfull) :: xp,xm,dumt0,umag
real, dimension(wfull), intent(in) :: a_f
real, dimension(wfull), intent(inout) :: d_b0,d_ustar,d_wu0,d_wv0,d_wt0,d_ws0,d_wm0,d_zcr

call getmixdepth(d_rho,d_nsq,d_rad,d_alpha,d_beta,d_b0,d_ustar,a_f,d_zcr) ! solve for mixed layer depth (calculated at full levels)
call getstab(km,ks,gammas,d_nsq,d_ustar,d_zcr)                     ! solve for stability functions and non-local term
                                                                   ! (calculated at half levels)

! POTENTIAL TEMPERATURE
if (incradgam.gt.0) then
  do iqw=1,wfull
    dumt0(iqw)=d_wt0(iqw)+sum(d_rad(iqw,1:p_mixind(iqw)))
  end do
else
  dumt0=d_wt0
end if
! +ve sign for rhs terms since z +ve is down
rhs(:,1)=ks(:,2)*gammas(:,2)/(dz(:,1)*d_zcr)
do ii=2,wlev-1
  rhs(:,ii)=(ks(:,ii+1)*gammas(:,ii+1)-ks(:,ii)*gammas(:,ii))/(dz(:,ii)*d_zcr)
end do
rhs(:,wlev)=-ks(:,wlev)*gammas(:,wlev)/(dz(:,wlev)*d_zcr)
cc(:,1)=-dt*ks(:,2)/(dz_hl(:,2)*dz(:,1)*d_zcr*d_zcr)
bb(:,1)=1.-cc(:,1)
dd(:,1)=w_temp(:,1)+dt*rhs(:,1)*dumt0-dt*d_rad(:,1)/(dz(:,1)*d_zcr)-dt*d_wt0/(dz(:,1)*d_zcr)
do ii=2,wlev-1
  aa(:,ii)=-dt*ks(:,ii)/(dz_hl(:,ii)*dz(:,ii)*d_zcr*d_zcr)
  cc(:,ii)=-dt*ks(:,ii+1)/(dz_hl(:,ii+1)*dz(:,ii)*d_zcr*d_zcr)
  bb(:,ii)=1.-aa(:,ii)-cc(:,ii)
  dd(:,ii)=w_temp(:,ii)+dt*rhs(:,ii)*dumt0-dt*d_rad(:,ii)/(dz(:,ii)*d_zcr)
end do
aa(:,wlev)=-dt*ks(:,wlev)/(dz_hl(:,wlev)*dz(:,wlev)*d_zcr*d_zcr)
bb(:,wlev)=1.-aa(:,wlev)
dd(:,wlev)=w_temp(:,wlev)+dt*rhs(:,wlev)*dumt0-dt*d_rad(:,wlev)/(dz(:,wlev)*d_zcr)
call thomas(w_temp,aa,bb,cc,dd)


! SALINITY
do ii=1,wlev
  dd(:,ii)=w_sal(:,ii)+dt*rhs(:,ii)*d_ws0
end do
dd(:,1)=dd(:,1)-dt*d_ws0/(dz(:,1)*d_zcr)
if (salrelax.eq.1) then ! relax salinity
  where (w_sal.gt.1.E-6)
    dd=dd+dt*(34.72-w_sal)/(3600.*24.*365.25)
  end where
end if
call thomas(w_sal,aa,bb,cc,dd)
w_sal=max(0.,w_sal)


! split U diffusion term
cc(:,1)=-dt*km(:,2)/(dz_hl(:,2)*dz(:,1)*d_zcr*d_zcr)
bb(:,1)=1.-cc(:,1)
dd(:,1)=w_u(:,1)-dt*d_wu0/(dz(:,1)*d_zcr)
do ii=2,wlev-1
  aa(:,ii)=-dt*km(:,ii)/(dz_hl(:,ii)*dz(:,ii)*d_zcr*d_zcr)
  cc(:,ii)=-dt*km(:,ii+1)/(dz_hl(:,ii+1)*dz(:,ii)*d_zcr*d_zcr)
  bb(:,ii)=1.-aa(:,ii)-cc(:,ii)
  dd(:,ii)=w_u(:,ii)
end do
aa(:,wlev)=-dt*km(:,wlev)/(dz_hl(:,wlev)*dz(:,wlev)*d_zcr*d_zcr)
bb(:,wlev)=1.-aa(:,wlev)
dd(:,wlev)=w_u(:,wlev)
umag=sqrt(w_u(:,wlev)*w_u(:,wlev)+w_v(:,wlev)*w_v(:,wlev))
! bottom drag
where (depth_hl(:,wlev+1).lt.mxd)
  bb(:,wlev)=bb(:,wlev)+dt*cdbot*umag/(dz(:,wlev)*d_zcr)
end where
call thomas(w_u,aa,bb,cc,dd)


! split V diffusion term
do ii=1,wlev
  dd(:,ii)=w_v(:,ii)
end do
dd(:,1)=dd(:,1)-dt*d_wv0/(dz(:,1)*d_zcr)
call thomas(w_v,aa,bb,cc,dd)


! --- Turn off coriolis terms as this is processed in mlodynamics.f90 ---
!! Split U and V coriolis terms
!xp=1.+(0.5*dt*a_f)**2
!xm=1.-(0.5*dt*a_f)**2
!do ii=1,wlev
!  newa=(w_u(:,ii)*xm+w_v(:,ii)*dt*a_f)/xp
!  newb=(w_v(:,ii)*xm-w_u(:,ii)*dt*a_f)/xp
!  w_u(:,ii)=newa
!  w_v(:,ii)=newb
!end do
w_u=max(-100.,min(w_u,100.)) ! MJT suggestion
w_v=max(-100.,min(w_v,100.)) ! MJT suggestion

! adjust surface height
select case(deprelax)
  case(0) ! free surface height
    w_eta=w_eta+dt*d_wm0
  case(1) ! relax surface height
    d_wm0=d_wm0-w_eta/(3600.*24.*365.25)
    w_eta=w_eta+dt*d_wm0
  case(2) ! fix surface height
    w_eta=0.
  case DEFAULT
    write(6,*) "ERROR: Invalid deprelax ",deprelax
    stop
end select

! min limit of 1m
w_eta=max(w_eta,-depth_hl(:,wlev+1)+minwater)

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

subroutine getstab(km,ks,gammas,d_nsq,d_ustar,d_zcr)

implicit none

integer ii,jj,iqw
integer, dimension(wfull) :: mixind_hl
real, dimension(wfull,wlev), intent(out) :: km,ks,gammas
real, dimension(wfull,wlev) :: num,nus,wm,ws,ri
real, dimension(wfull) :: sigma,d_depth_hl,d_depth_hlp1
real, dimension(wfull) :: a2m,a3m,a2s,a3s
real, dimension(wfull) :: numh,wm1,dnumhdz,dwm1ds,g1m,dg1mds
real, dimension(wfull) :: nush,ws1,dnushdz,dws1ds,g1s,dg1sds
real xp,cg
real, dimension(wfull,wlev), intent(in) :: d_nsq
real, dimension(wfull), intent(in) :: d_ustar
real, dimension(wfull), intent(in) :: d_zcr
real, parameter :: ri0 = 0.7
real, parameter :: nu0 = 50.E-4
real, parameter :: numw = 1.0E-4
real, parameter :: nusw = 0.1E-4

ri=0. ! ri(:,1) is not used
do ii=2,wlev
  ri(:,ii)=d_nsq(:,ii)*(dz_hl(:,ii)*d_zcr)**2    & 
           /max((w_u(:,ii-1)-w_u(:,ii))**2       &
               +(w_v(:,ii-1)-w_v(:,ii))**2,1.E-20)
end do

! diffusion ---------------------------------------------------------
num=0.
nus=0.
do ii=2,wlev
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
num(:,1)=num(:,2) ! to avoid problems calculating shallow mixed layer
nus(:,1)=nus(:,2)
!--------------------------------------------------------------------

! stability ---------------------------------------------------------
wm=0.
ws=0.
do ii=2,wlev
  d_depth_hl=depth_hl(:,ii)*d_zcr
  call getwx(wm(:,ii),ws(:,ii),d_depth_hl,p_bf,d_ustar,p_mixdepth)
end do
wm(:,1)=wm(:,2) ! to avoid problems calculating shallow mixed layer
ws(:,1)=ws(:,2)
!--------------------------------------------------------------------

! calculate G profile -----------------------------------------------
! -ve as z is down
do iqw=1,wfull
  mixind_hl(iqw)=p_mixind(iqw)
  jj=min(mixind_hl(iqw)+1,wlev-1)
  d_depth_hl(iqw)=depth_hl(iqw,jj)*d_zcr(iqw)
  if (p_mixdepth(iqw).gt.d_depth_hl(iqw)) mixind_hl(iqw)=jj
  d_depth_hl(iqw)=depth_hl(iqw,mixind_hl(iqw))*d_zcr(iqw)
  d_depth_hlp1(iqw)=depth_hl(iqw,mixind_hl(iqw)+1)*d_zcr(iqw)
  xp=(p_mixdepth(iqw)-d_depth_hl(iqw))/(d_depth_hlp1(iqw)-d_depth_hl(iqw))
  xp=max(0.,min(1.,xp))
  numh(iqw)=(1.-xp)*num(iqw,mixind_hl(iqw))+xp*num(iqw,mixind_hl(iqw)+1)
  nush(iqw)=(1.-xp)*nus(iqw,mixind_hl(iqw))+xp*nus(iqw,mixind_hl(iqw)+1)
  wm1(iqw)=(1.-xp)*wm(iqw,mixind_hl(iqw))+xp*wm(iqw,mixind_hl(iqw)+1)
  ws1(iqw)=(1.-xp)*ws(iqw,mixind_hl(iqw))+xp*ws(iqw,mixind_hl(iqw)+1)
  dnumhdz(iqw)=min(num(iqw,mixind_hl(iqw)+1)-num(iqw,mixind_hl(iqw)),0.)/(dz_hl(iqw,mixind_hl(iqw)+1)*d_zcr(iqw))
  dnushdz(iqw)=min(nus(iqw,mixind_hl(iqw)+1)-nus(iqw,mixind_hl(iqw)),0.)/(dz_hl(iqw,mixind_hl(iqw)+1)*d_zcr(iqw))
  !dwm1ds and dws1ds are now multipled by 1/(mixdepth*wx1*wx1)
  dwm1ds(iqw)=-5.*max(p_bf(iqw),0.)/max(d_ustar(iqw)**4,1.E-20)
  dws1ds(iqw)=dwm1ds(iqw)
end do

g1m=numh/(p_mixdepth*wm1)
g1s=nush/(p_mixdepth*ws1)
dg1mds=dnumhdz/wm1-numh*dwm1ds
dg1sds=dnushdz/ws1-nush*dws1ds
  
a2m=-2.+3.*g1m-dg1mds
a2s=-2.+3.*g1s-dg1sds
a3m=1.-2.*g1m+dg1mds
a3s=1.-2.*g1s+dg1sds

!--------------------------------------------------------------------
! calculate diffusion terms
km=num
ks=nus
do ii=2,wlev
  where (ii.le.mixind_hl)
    sigma=depth_hl(:,ii)*d_zcr/p_mixdepth
    km(:,ii)=max(p_mixdepth*wm(:,ii)*sigma*(1.+sigma*(a2m+a3m*sigma)),num(:,ii))
    ks(:,ii)=max(p_mixdepth*ws(:,ii)*sigma*(1.+sigma*(a2s+a3s*sigma)),nus(:,ii))
  end where
end do

!--------------------------------------------------------------------
! non-local term
! gammas is the same for temp and sal when double-diffusion is not employed
cg=10.*vkar*(98.96*vkar*epsilon)**(1./3.) ! Large (1994)
!cg=5.*vkar*(98.96*vkar*epsilon)**(1./3.) ! Bernie (2004)
gammas=0.
do ii=2,wlev
  where (p_bf.lt.0..and.ii.le.mixind_hl) ! unstable
    gammas(:,ii)=cg/max(ws(:,ii)*p_mixdepth,1.E-20)
  end where
end do

return
end subroutine getstab

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate mixed layer depth

subroutine getmixdepth(d_rho,d_nsq,d_rad,d_alpha,d_beta,d_b0,d_ustar,a_f,d_zcr)

implicit none

integer ii,jj,kk,iqw
real vtc,dvsq,vtsq,xp
real tnsq,tws,twu,twv,tdepth,tbuoy,trho
real oldxp,oldtrib,trib,newxp
real, dimension(wfull,wlev) :: ws,wm,dumbuoy,rib
real, dimension(wfull) :: dumbf,l,d_depth
real, dimension(wfull,wlev), intent(in) :: d_rho,d_nsq,d_rad,d_alpha,d_beta
real, dimension(wfull), intent(in) :: d_b0,d_ustar,d_zcr
real, dimension(wfull), intent(in) :: a_f
integer, parameter :: maxits = 3
real, parameter :: alpha = 0.9

vtc=1.8*sqrt(0.2/(98.96*epsilon))/(vkar*vkar*ric)

if (incradbf.gt.0) then
  do ii=1,wlev
    dumbf=d_b0-grav*sum(d_alpha(:,1:ii)*d_rad(:,1:ii),2) ! -ve sign is to account for sign of d_rad
    d_depth=depth(:,ii)*d_zcr
    call getwx(wm(:,ii),ws(:,ii),d_depth,dumbf,d_ustar,d_depth)
  end do
else
  dumbf=d_b0
  do ii=1,wlev
    d_depth=depth(:,ii)*d_zcr
    call getwx(wm(:,ii),ws(:,ii),d_depth,dumbf,d_ustar,d_depth)
  end do
end if

if (buoymeth.eq.0) then
  do ii=1,wlev
    dumbuoy(:,ii)=grav*(d_rho(:,ii)-d_rho(:,1))
  end do
else
  do ii=1,wlev
    dumbuoy(:,ii)=grav*(d_alpha(:,ii)*w_temp(:,ii)-d_beta(:,ii)*w_sal(:,ii))
  end do
  do ii=2,wlev
    dumbuoy(:,ii)=(dumbuoy(:,ii)-dumbuoy(:,1))
  end do
  dumbuoy(:,1)=0.
end if

p_mixind=wlev-1
p_mixdepth=depth(:,wlev)*d_zcr
rib=0.
do iqw=1,wfull
  do ii=2,wlev
    jj=min(ii+1,wlev)
    vtsq=depth(iqw,ii)*d_zcr(iqw)*ws(iqw,ii)*sqrt(0.5*max(d_nsq(iqw,ii)+d_nsq(iqw,jj),0.))*vtc
    dvsq=(w_u(iqw,1)-w_u(iqw,ii))**2+(w_v(iqw,1)-w_v(iqw,ii))**2
    rib(iqw,ii)=(depth(iqw,ii)-depth(iqw,1))*d_zcr(iqw)*dumbuoy(iqw,ii)/(max(dvsq+vtsq,1.E-20)*d_rho(iqw,ii))
    if (rib(iqw,ii).gt.ric.or.abs(rib(iqw,ii)).lt.rib(iqw,ii-1)) then
      p_mixind(iqw)=ii-1
      xp=min(max((ric-rib(iqw,ii-1))/max(rib(iqw,ii)-rib(iqw,ii-1),1.E-20),0.),1.)
      p_mixdepth(iqw)=((1.-xp)*depth(iqw,ii-1)+xp*depth(iqw,ii))*d_zcr(iqw)
      exit
    end if
  end do 
end do

if (mixmeth.eq.1) then
  ! Refine mixed-layer-depth calculation
  do iqw=1,wfull
    ii=p_mixind(iqw)+1
    jj=min(ii+1,wlev)
    oldxp=0.
    oldtrib=rib(iqw,ii-1)
    xp=(p_mixdepth(iqw)-depth(iqw,ii-1)*d_zcr(iqw))/((depth(iqw,ii)-depth(iqw,ii-1))*d_zcr(iqw))
    xp=max(xp,oldxp+0.1)
    do kk=1,maxits
      if (xp.lt.0.5) then
        tnsq=(1.-2.*xp)*0.5*max(d_nsq(iqw,ii-1)+d_nsq(iqw,ii),0.)+(2.*xp)*max(d_nsq(iqw,ii),0.)
      else
        tnsq=(2.-2.*xp)*max(d_nsq(iqw,ii),0.)+(2.*xp-1.)*0.5*max(d_nsq(iqw,ii)+d_nsq(iqw,jj),0.)
      end if
      tws=(1.-xp)*ws(iqw,ii-1)+xp*ws(iqw,ii)
      twu=(1.-xp)*w_u(iqw,ii-1)+xp*w_u(iqw,ii)
      twv=(1.-xp)*w_v(iqw,ii-1)+xp*w_v(iqw,ii)
      tdepth=((1.-xp)*depth(iqw,ii-1)+xp*depth(iqw,ii))*d_zcr(iqw)
      tbuoy=(1.-xp)*dumbuoy(iqw,ii-1)+xp*dumbuoy(iqw,ii)
      trho=(1.-xp)*d_rho(iqw,ii-1)+xp*d_rho(iqw,ii)
      vtsq=tdepth*tws*sqrt(tnsq)*vtc
      dvsq=(w_u(iqw,1)-twu)**2+(w_v(iqw,1)-twv)**2
      trib=(tdepth-depth(iqw,1)*d_zcr(iqw))*tbuoy/(max(dvsq+vtsq,1.E-20)*trho)
      if (abs(trib-oldtrib).gt.1.E-5) then
        newxp=xp-alpha*(trib-ric)*(xp-oldxp)/(trib-oldtrib) ! i.e., (trib-ric-oldtrib+ric)
        oldtrib=trib
        oldxp=xp
        xp=newxp
        xp=min(max(xp,0.),1.)
      else
        exit
      end if
    end do
    p_mixdepth(iqw)=((1.-xp)*depth(iqw,ii-1)+xp*depth(iqw,ii))*d_zcr(iqw)
  end do
end if

! calculate buoyancy forcing
call getbf(d_rad,d_alpha,d_b0)

! impose limits for stable conditions
where(p_bf.gt.1.E-10)
  l=d_ustar*d_ustar*d_ustar/(vkar*p_bf)
  p_mixdepth=min(p_mixdepth,l)
end where
where(p_bf.gt.1.E-10.and.abs(a_f).gt.1.E-10)
  l=0.7*d_ustar/abs(a_f)
  p_mixdepth=min(p_mixdepth,l)
end where
p_mixdepth=max(p_mixdepth,depth(:,1)*d_zcr)
p_mixdepth=min(p_mixdepth,depth(:,wlev)*d_zcr)

! recalculate index for mixdepth
p_mixind=wlev-1
do iqw=1,wfull
  do ii=2,wlev
    if (depth(iqw,ii)*d_zcr(iqw).gt.p_mixdepth(iqw)) then
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
    p_bf(iqw)=d_b0(iqw)-grav*sum(d_alpha(iqw,1:p_mixind(iqw))*d_rad(iqw,1:p_mixind(iqw))) ! -ve sign is to account for sign of d_rad
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
uuu=max(d_ustar**3,1.E-20)
invl=vkar*bf      ! invl = ustar**3/L or L=ustar**3/(vkar*bf)
zeta=sig*mixdp*invl
zeta=min(zeta,uuu) ! MJT suggestion

where (zeta.gt.0.)
  wm=vkar*d_ustar*uuu/(uuu+5.*zeta)
elsewhere (zeta.gt.zetam*uuu)
  wm=vkar*(d_ustar*uuu-16.*d_ustar*zeta)**(1./4.)
elsewhere
  wm=vkar*(am*uuu-cm*zeta)**(1./3.)
end where

where (zeta.gt.0.)
  ws=vkar*d_ustar*uuu/(uuu+5.*zeta)
elsewhere (zeta.gt.zetas*uuu)
  ws=vkar*(d_ustar*d_ustar-16.*zeta/d_ustar)**(1./2.)
elsewhere
  ws=vkar*(as*uuu-cs*zeta)**(1./3.)
end where

wm=max(wm,1.E-20)
ws=max(ws,1.E-20)

return
end subroutine getwx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate rho from equation of state
! From GFDL (MOM3)

subroutine getrho(a_ps,d_rho,d_rs,d_alpha,d_beta,d_zcr)

implicit none

integer ii
real, dimension(wfull) :: rho0,pxtr
real, dimension(wfull,wlev) :: d_depth,d_dz
real, dimension(wfull,wlev), intent(inout) :: d_rho,d_rs,d_alpha,d_beta
real, dimension(wfull), intent(in) :: a_ps
real, dimension(wfull), intent(inout) :: d_zcr

pxtr=a_ps+grav*i_fracice*(i_dic*rhoic+i_dsn*rhosn)
do ii=1,wlev
  d_depth(:,ii)=depth(:,ii)*d_zcr
  d_dz(:,ii)=dz(:,ii)*d_zcr
end do
call calcdensity(wfull,d_rho,d_alpha,d_beta,d_rs,rho0,w_temp,w_sal,d_depth,d_dz,pxtr)

return
end subroutine getrho

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate fluxes between water and atmosphere

subroutine getwflux(a_sg,a_rg,a_rnd,a_snd,a_vnratio,a_fbvis,a_fbnir,a_inflow,d_rho,d_rs,d_nsq,d_rad,d_alpha, &
                    d_beta,d_b0,d_ustar,d_wu0,d_wv0,d_wt0,d_ws0,d_wm0,d_taux,d_tauy,d_zcr)

implicit none

integer ii
real, dimension(wfull) :: visalb,niralb
real, dimension(wfull,wlev), intent(inout) :: d_rho,d_rs,d_nsq,d_rad,d_alpha,d_beta
real, dimension(wfull), intent(inout) :: d_b0,d_ustar,d_wu0,d_wv0,d_wt0,d_ws0,d_wm0,d_taux,d_tauy,d_zcr
real, dimension(wfull), intent(in) :: a_sg,a_rg,a_rnd,a_snd,a_vnratio,a_fbvis,a_fbnir,a_inflow

! buoyancy frequency (calculated at half levels)
do ii=2,wlev
  d_nsq(:,ii)=-2.*grav/(d_rho(:,ii-1)+d_rho(:,ii))*(d_rho(:,ii-1)-d_rho(:,ii))/(dz_hl(:,ii)*d_zcr)
end do
d_nsq(:,1)=2.*d_nsq(:,2)-d_nsq(:,3) ! not used

! shortwave
! use -ve as depth is down
visalb=p_watervisdiralb*a_fbvis+p_watervisdifalb*(1.-a_fbvis)
niralb=p_waternirdiralb*a_fbnir+p_waternirdifalb*(1.-a_fbnir)
d_rad(:,1)=a_sg/(cp0*d_rho(:,1))*(((1.-visalb)*a_vnratio*exp(-depth_hl(:,2)*d_zcr/mu_1)+       &
                                   (1.-niralb)*(1.-a_vnratio)*exp(-depth_hl(:,2)*d_zcr/mu_2))- &
                                   ((1.-visalb)*a_vnratio+(1.-niralb)*(1.-a_vnratio)))
do ii=2,wlev-1 
  d_rad(:,ii)=a_sg/(cp0*d_rho(:,ii))*(((1.-visalb)*a_vnratio*exp(-depth_hl(:,ii+1)*d_zcr/mu_1)+       &
                                       (1.-niralb)*(1.-a_vnratio)*exp(-depth_hl(:,ii+1)*d_zcr/mu_2))- &
                                      ((1.-visalb)*a_vnratio*exp(-depth_hl(:,ii)*d_zcr/mu_1)+         &
                                       (1.-niralb)*(1.-a_vnratio)*exp(-depth_hl(:,ii)*d_zcr/mu_2)))
end do
d_rad(:,wlev)=-a_sg/(cp0*d_rho(:,wlev))*((1.-visalb)*a_vnratio*exp(-depth_hl(:,wlev)*d_zcr/mu_1)+ &
                                         (1.-niralb)*(1.-a_vnratio)*exp(-depth_hl(:,wlev)*d_zcr/mu_2)) ! remainder
do ii=1,wlev
  d_rad(:,ii)=(1.-i_fracice)*d_rad(:,ii)
end do

! Boundary conditions (use rho at level 1 for consistancy with shortwave radiation)
d_wu0=-(1.-i_fracice)*d_taux/d_rho(:,1)                                            ! BC
d_wv0=-(1.-i_fracice)*d_tauy/d_rho(:,1)                                            ! BC
d_wt0=-(1.-i_fracice)*(-p_fg-p_eg+a_rg-sbconst*w_temp(:,1)**4)/(d_rho(:,1)*cp0)    ! BC
d_wt0=d_wt0+(1.-i_fracice)*lf*a_snd/(d_rho(:,1)*cp0) ! melting snow                ! BC
d_ws0=(1.-i_fracice)*(a_rnd+a_snd-p_eg/lv)*w_sal(:,1)/d_rs(:,1)                    ! BC
d_ws0=d_ws0+(a_inflow)*w_sal(:,1)/d_rs(:,1) ! currently salinity transfer between ice and water is neglected
d_wm0=((1.-i_fracice)*(a_rnd+a_snd-p_eg/lv)+a_inflow)*0.001 ! convert mm to m      ! BC

d_ustar=max(sqrt(sqrt(d_wu0*d_wu0+d_wv0*d_wv0)),1.E-6)
d_b0=-grav*(d_alpha(:,1)*d_wt0-d_beta(:,1)*d_ws0) ! -ve sign to compensate for sign of wt0 and ws0
                                                  ! This is the opposite sign used by Large for wb0, but
                                                  ! is same sign as Large used for Bf (+ve stable, -ve unstable)

return
end subroutine getwflux

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate density

subroutine calcdensity(full,d_rho,d_alpha,d_beta,rs,rho0,tt,ss,dep,ddz,pxtr)

implicit none

integer, intent(in) :: full
integer i,ii
!integer, parameter :: nits=1 ! iterate for density (nits=1 recommended)
real, dimension(full,wlev), intent(in) :: tt,ss,dep,ddz
real, dimension(full,wlev), intent(out) :: d_rho,d_alpha,d_beta,rs
real, dimension(full), intent(in) :: pxtr
real, dimension(full), intent(out) :: rho0
real, dimension(full) :: t,s,p1,p2,t2,t3,t4,t5,s2,s3,s32,ptot
real, dimension(full) :: drho0dt,drho0ds,dskdt,dskds,sk,sks
real, dimension(full) :: drhodt,drhods,rs0
real, parameter :: density = 1035.

d_rho=density

t = max(tt(:,1)-273.16,-2.)
s = max(ss(:,1),0.)
t2 = t*t
t3 = t2*t
t4 = t3*t
t5 = t4*t
s2 = s*s
s3 = s2*s
s32 = sqrt(s3)

rs0 = 999.842594 + 6.793952e-2*t(:)            &
       - 9.095290e-3*t2(:) + 1.001685e-4*t3(:) &
       - 1.120083e-6*t4(:) + 6.536332e-9*t5(:) ! density for sal=0.
rho0 = rs0+ s(:)*(0.824493 - 4.0899e-3*t(:)    &
       + 7.6438e-5*t2(:)                       &
       - 8.2467e-7*t3(:) + 5.3875e-9*t4(:))    &
       + s32(:)*(-5.72466e-3 + 1.0227e-4*t(:)  &
       - 1.6546e-6*t2(:)) + 4.8314e-4*s2(:)     ! + sal terms    
drho0dt=6.793952e-2                                  &
       - 2.*9.095290e-3*t(:) + 3.*1.001685e-4*t2(:)  &
       - 4.*1.120083e-6*t3(:) + 5.*6.536332e-9*t4(:) &
       + s(:)*( - 4.0899e-3 + 2.*7.6438e-5*t(:)      &
       - 3.*8.2467e-7*t2(:) + 4.*5.3875e-9*t3(:))    &
       + s32(:)*(1.0227e-4 - 2.*1.6546e-6*t(:))
drho0ds= (0.824493 - 4.0899e-3*t(:) + 7.6438e-5*t2(:) &
       - 8.2467e-7*t3(:) + 5.3875e-9*t4(:))           &
       + 1.5*sqrt(s(:))*(-5.72466e-3 + 1.0227e-4*t(:) &
       - 1.6546e-6*t2(:)) + 2.*4.8314e-4*s(:)

!do i=1,nits
  ptot=pxtr*1.E-5
  do ii=1,wlev
    t = max(tt(:,ii)-273.16,-2.)
    s = max(ss(:,ii),0.)
    p1   = ptot+grav*d_rho(:,ii)*0.5*ddz(:,ii)*1.E-5 ! hydrostatic approximation
    ptot = ptot+grav*d_rho(:,ii)*ddz(:,ii)*1.E-5
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
!    dskdp=(3.186519 + 2.212276e-2*t(:) &
!                - 2.984642e-4*t2(:) + 1.956415e-6*t3(:)) &
!                + 2.*p1*(2.102898e-4 - 1.202016e-5*t(:) &
!                + 1.394680e-7*t2(:)) &
!                + s(:)*(6.704388e-3  -1.847318e-4*t(:) &
!                + 2.059331e-7*t2(:)) &
!                + 1.480266e-4*s32(:) &
!                + 2.*p1(:)*s(:)*(-2.040237e-6 &
!                + 6.128773e-8*t(:) + 6.207323e-10*t2(:))
       
    d_rho(:,ii)=rho0/(1.-p1/sk)
    rs(:,ii)=rs0/(1.-p1/sks) ! sal=0.
  
    drhodt=drho0dt/(1.-p1/sk)-rho0*p1*dskdt/((sk-p1)**2) ! neglected dp1drho*drhodt terms
    drhods=drho0ds/(1.-p1/sk)-rho0*p1*dskds/((sk-p1)**2) ! neglected dp1drho*drhods terms
    
    d_alpha(:,ii)=-drhodt              ! Large et al (1993) convention
    d_beta(:,ii)=drhods                ! Large et al (1993) convention
  end do

!end do

return
end subroutine calcdensity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Convert temp to theta based on MOM3 routine
subroutine gettheta(theta,tt,ss,pxtr,ddz)

implicit none

integer ii
real, dimension(wfull,wlev), intent(out) :: theta
real, dimension(wfull,wlev), intent(in) :: tt,ss,ddz
real, dimension(wfull), intent(in) :: pxtr
real, dimension(wfull) :: b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11
real, dimension(wfull) :: t,t2,t3,s,s2,p,p2,potmp,ptot
real, dimension(wfull,wlev) :: d_rho
real, parameter :: density = 1035.

d_rho=density

ptot=pxtr*1.E-5
do ii=1,wlev
  t = max(tt(:,ii)-273.16,-2.)
  s = max(ss(:,ii),0.)
  p   = ptot+grav*d_rho(:,ii)*0.5*ddz(:,ii)*1.E-5 ! hydrostatic approximation
  ptot = ptot+grav*d_rho(:,ii)*ddz(:,ii)*1.E-5
    
  b1    = -1.60e-5*p
  b2    = 1.014e-5*p*t
  t2    = t*t
  t3    = t2*t
  b3    = -1.27e-7*p*t2
  b4    = 2.7e-9*p*t3
  b5    = 1.322e-6*p*s
  b6    = -2.62e-8*p*s*t
  s2    = s*s
  p2    = p*p
  b7    = 4.1e-9*p*s2
  b8    = 9.14e-9*p2
  b9    = -2.77e-10*p2*t
  b10   = 9.5e-13*p2*t2
  b11   = -1.557e-13*p2*p
  potmp = b1+b2+b3+b4+b5+b6+b7+b8+b9+b10+b11
  theta(:,ii) = t-potmp+273.16
end do

return
end subroutine gettheta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate fluxes between MLO and atmosphere (from CCAM sflux.f)

subroutine fluxcalc(dt,a_u,a_v,a_temp,a_qg,a_ps,a_zmin,a_zmins,d_rho,d_taux,d_tauy,d_zcr,diag)

implicit none

integer, intent(in) :: diag
integer it
real, intent(in) :: dt
real, dimension(wfull,wlev), intent(inout) :: d_rho
real, dimension(wfull), intent(in) :: a_u,a_v,a_temp,a_qg,a_ps,a_zmin,a_zmins
real, dimension(wfull), intent(inout) :: d_taux,d_tauy,d_zcr
real, dimension(wfull) :: qsat,dqdt,ri,vmag,rho,srcp
real, dimension(wfull) :: fm,fh,con,consea,afroot,af,daf
real, dimension(wfull) :: den,dfm,dden,dcon,sig,factch,root
real, dimension(wfull) :: aft,atu,atv,dcs,afq,facqch,fq
real, dimension(wfull) :: newwu,newwv,a,vmagn
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
vmagn=max(vmag,0.2)
sig=exp(-grav*a_zmins/(rdry*a_temp))
srcp=sig**(rdry/cp)
rho=a_ps/(rdry*w_temp(:,1))

call getqsat(qsat,dqdt,w_temp(:,1),a_ps)
if (zomode.eq.0) then ! CSIRO9
  qsat=0.98*qsat ! with Zeng 1998 for sea water
end if
ri=min(grav*(a_zmin**2/a_zmins)*(1.-w_temp(:,1)*srcp/a_temp)/vmagn**2,rimax)
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
        den=1.+af*cms*2.*bprm*sqrt(-ri*a_zmin/p_zo)
        fm=1.-2.*bprm*ri/den
        con=consea*fm*vmag
        dden=daf*cms*2.*bprm*sqrt(-ri*a_zmin/p_zo)+af*cms*bprm*ri*a_zmin/(sqrt(-ri*a_zmin/p_zo)*p_zo*p_zo)
        dfm=2.*bprm*ri*dden/(den**2)
        dcon=consea*dfm*vmag
        p_zo=p_zo-(p_zo-con*af)/(1.-dcon*af-con*daf)
      end where
      p_zo=min(max(p_zo,1.5e-5),13.)
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
      p_zo=min(max(p_zo,1.5e-5),13.)
    enddo    ! it=1,4
end select
afroot=vkar/log(a_zmin/p_zo)
af=afroot**2

select case(zomode)
  case(0) ! Charnock CSIRO9
    ztv=exp(vkar/sqrt(chn10))/10.
    aft=(vkar/log(a_zmins*ztv))**2
    afq=aft
    p_zoh=1./ztv
    p_zoq=p_zoh
    factch=sqrt(p_zo/p_zoh)
    facqch=factch
  case(1) ! Charnock zot=zom
    p_zoh=p_zo
    p_zoq=p_zo
    aft=vkar**2/(log(a_zmins/p_zo)*log(a_zmins/p_zo))
    afq=aft
    factch=1.
    facqch=1.
  case(2) ! Beljaars
    p_zoh=max(zcoh1+zcoh2*gnu/max(vmag*sqrt(fm*af),gnu),1.5E-7)
    p_zoq=max(zcoq1+zcoq2*gnu/max(vmag*sqrt(fm*af),gnu),1.5E-7)
    factch=sqrt(p_zo/p_zoh)
    facqch=sqrt(p_zo/p_zoq)
    aft=vkar**2/(log(a_zmins/p_zo)*log(a_zmins/p_zoh))
    afq=vkar**2/(log(a_zmins/p_zo)*log(a_zmins/p_zoq))
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

! update fluxes with atmosphere
p_cd=af*fm
p_cds=aft*fh
p_fg=rho*aft*cp*fh*vmag*(w_temp(:,1)-a_temp/srcp)
p_eg=rho*afq*lv*fq*vmag*(qsat-a_qg)
d_taux=rho*p_cd*vmag*atu
d_tauy=rho*p_cd*vmag*atv

! turn off lake evaporation when minimum depth is reached
! fg should be replaced with bare ground value
where (depth_hl(:,wlev+1)+w_eta.le.minwater)
  p_eg=0.
end where

return
end subroutine fluxcalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Estimate saturation mixing ratio (from CCAM)
! following version of getqsat copes better with T < 0C over ice
subroutine getqsat(qsat,dqdt,temp,ps)

implicit none

real, dimension(wfull), intent(in) :: temp,ps
real, dimension(wfull), intent(out) :: qsat,dqdt
real, dimension(0:220) :: table
real, dimension(wfull) :: esatf,tdiff,dedt,rx
integer, dimension(wfull) :: ix

table(0:4)=    (/ 1.e-9, 1.e-9, 2.e-9, 3.e-9, 4.e-9 /)                                !-146C
table(5:9)=    (/ 6.e-9, 9.e-9, 13.e-9, 18.e-9, 26.e-9 /)                             !-141C
table(10:14)=  (/ 36.e-9, 51.e-9, 71.e-9, 99.e-9, 136.e-9 /)                          !-136C
table(15:19)=  (/ 0.000000188, 0.000000258, 0.000000352, 0.000000479, 0.000000648 /)  !-131C
table(20:24)=  (/ 0.000000874, 0.000001173, 0.000001569, 0.000002090, 0.000002774 /)  !-126C
table(25:29)=  (/ 0.000003667, 0.000004831, 0.000006340, 0.000008292, 0.00001081 /)   !-121C
table(30:34)=  (/ 0.00001404, 0.00001817, 0.00002345, 0.00003016, 0.00003866 /)       !-116C
table(35:39)=  (/ 0.00004942, 0.00006297, 0.00008001, 0.0001014, 0.0001280 /)         !-111C
table(40:44)=  (/ 0.0001613, 0.0002026, 0.0002538, 0.0003170, 0.0003951 /)            !-106C
table(45:49)=  (/ 0.0004910, 0.0006087, 0.0007528, 0.0009287, 0.001143 /)             !-101C
table(50:55)=  (/ .001403, .001719, .002101, .002561, .003117, .003784 /)             !-95C
table(56:63)=  (/ .004584, .005542, .006685, .008049, .009672,.01160,.01388,.01658 /) !-87C
table(64:72)=  (/ .01977, .02353, .02796,.03316,.03925,.04638,.05472,.06444,.07577 /) !-78C
table(73:81)=  (/ .08894, .1042, .1220, .1425, .1662, .1936, .2252, .2615, .3032 /)   !-69C
table(82:90)=  (/ .3511, .4060, .4688, .5406, .6225, .7159, .8223, .9432, 1.080 /)    !-60C
table(91:99)=  (/ 1.236, 1.413, 1.612, 1.838, 2.092, 2.380, 2.703, 3.067, 3.476 /)    !-51C
table(100:107)=(/ 3.935,4.449, 5.026, 5.671, 6.393, 7.198, 8.097, 9.098 /)            !-43C
table(108:116)=(/ 10.21, 11.45, 12.83, 14.36, 16.06, 17.94, 20.02, 22.33, 24.88 /)    !-34C
table(117:126)=(/ 27.69, 30.79, 34.21, 37.98, 42.13, 46.69,51.70,57.20,63.23,69.85 /) !-24C 
table(127:134)=(/ 77.09, 85.02, 93.70, 103.20, 114.66, 127.20, 140.81, 155.67 /)      !-16C
table(135:142)=(/ 171.69, 189.03, 207.76, 227.96 , 249.67, 272.98, 298.00, 324.78 /)  !-8C
table(143:150)=(/ 353.41, 383.98, 416.48, 451.05, 487.69, 526.51, 567.52, 610.78 /)   !0C
table(151:158)=(/ 656.62, 705.47, 757.53, 812.94, 871.92, 934.65, 1001.3, 1072.2 /)   !8C
table(159:166)=(/ 1147.4, 1227.2, 1311.9, 1401.7, 1496.9, 1597.7, 1704.4, 1817.3 /)   !16C
table(167:174)=(/ 1936.7, 2063.0, 2196.4, 2337.3, 2486.1, 2643.0, 2808.6, 2983.1 /)   !24C
table(175:182)=(/ 3167.1, 3360.8, 3564.9, 3779.6, 4005.5, 4243.0, 4492.7, 4755.1 /)   !32C
table(183:190)=(/ 5030.7, 5320.0, 5623.6, 5942.2, 6276.2, 6626.4, 6993.4, 7377.7 /)   !40C
table(191:197)=(/ 7780.2, 8201.5, 8642.3, 9103.4, 9585.5, 10089.0, 10616.0 /)         !47C
table(198:204)=(/ 11166.0, 11740.0, 12340.0, 12965.0, 13617.0, 14298.0, 15007.0 /)    !54C
table(205:211)=(/ 15746.0, 16516.0, 17318.0, 18153.0, 19022.0, 19926.0, 20867.0 /)    !61C
table(212:218)=(/ 21845.0, 22861.0, 23918.0, 25016.0, 26156.0, 27340.0, 28570.0 /)    !68C
table(219:220)=(/ 29845.0, 31169.0 /)

tdiff=min(max( temp-123.16, 0.), 219.)
rx=tdiff-aint(tdiff)
ix=int(tdiff)
esatf=(1.-rx)*table(ix)+ rx*table(ix+1)
qsat=0.622*esatf/(ps-esatf)

! method #1
dedt=table(ix+1)-table(ix) ! divide by 1C
dqdt=qsat*dedt*ps/(esatf*(ps-esatf))
! method #2
!dqdt=qsat*lv/rvap*qsat*ps/(temp*temp*(ps-0.378*esatf))

return
end subroutine getqsat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Pack sea ice for calcuation

subroutine mloice(dt,a_rnd,a_snd,d_alpha,d_beta,d_b0,d_wu0,d_wv0,d_wt0,d_ws0,d_wm0,d_ftop,d_bot,d_tb,d_fb,d_timelt, &
                  d_tauxica,d_tauyica,d_tauxicw,d_tauyicw,d_ustar,d_rho,d_nk,diag)

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
real, dimension(wfull,wlev), intent(in) :: d_alpha,d_beta,d_rho
real, dimension(wfull), intent(inout) :: d_b0,d_wu0,d_wv0,d_wt0,d_ws0,d_wm0,d_ftop,d_bot,d_tb,d_fb,d_timelt
real, dimension(wfull), intent(inout) :: d_tauxicw,d_tauyicw,d_tauxica,d_tauyica,d_ustar
integer, dimension(wfull), intent(inout) :: d_nk
real, dimension(wfull) :: dp_ftop,dp_bot,dp_tb,dp_fb,dp_timelt,dp_salflx,dp_tauxica,dp_tauyica,dp_tauxicw,dp_tauyicw
real, dimension(wfull) :: dp_wtrflx
real, dimension(wfull) :: d_salflx,d_wtrflx
integer, dimension(wfull) :: dp_nk

d_salflx=0.
d_wtrflx=0.

! identify sea ice for packing
cice=i_dic.gt.icemin
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
dp_bot(1:nice)=pack(d_bot,cice)
dp_tb(1:nice)=pack(d_tb,cice)
dp_fb(1:nice)=pack(d_fb,cice)
dp_timelt(1:nice)=pack(d_timelt,cice)
dp_salflx(1:nice)=pack(d_salflx,cice)
dp_wtrflx(1:nice)=pack(d_wtrflx,cice)
dp_nk(1:nice)=pack(d_nk,cice)
dp_tauxica(1:nice)=pack(d_tauxica,cice)
dp_tauyica(1:nice)=pack(d_tauyica,cice)
dp_tauxicw(1:nice)=pack(d_tauxicw,cice)
dp_tauyicw(1:nice)=pack(d_tauyicw,cice)
! update ice prognostic variables
call seaicecalc(nice,dt,ip_tn(1:nice,:),ip_dic(1:nice),ip_dsn(1:nice),ip_fracice(1:nice),        &
                ip_tsurf(1:nice),ip_sto(1:nice),ap_rnd(1:nice),ap_snd(1:nice),pp_egice(1:nice),  &
                dp_ftop(1:nice),dp_bot(1:nice),dp_tb(1:nice),dp_fb(1:nice),dp_timelt(1:nice),    &
                dp_salflx(1:nice),dp_wtrflx(1:nice),dp_nk(1:nice),dp_tauxica(1:nice),            &
                dp_tauyica(1:nice),dp_tauxicw(1:nice),dp_tauyicw(1:nice),diag)
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
d_bot=unpack(dp_bot(1:nice),cice,d_bot)
d_tb=unpack(dp_tb(1:nice),cice,d_tb)
d_fb=unpack(dp_fb(1:nice),cice,d_fb)
d_timelt=unpack(dp_timelt(1:nice),cice,d_timelt)
d_salflx=unpack(dp_salflx(1:nice),cice,d_salflx)
d_wtrflx=unpack(dp_wtrflx(1:nice),cice,d_wtrflx)
d_nk=unpack(dp_nk(1:nice),cice,d_nk)

! update here because d_salflx is updated here
d_wu0=d_wu0-i_fracice*d_tauxicw/d_rho(:,1)
d_wv0=d_wv0-i_fracice*d_tauyicw/d_rho(:,1)
d_wt0=d_wt0+i_fracice*d_fb/(d_rho(:,1)*cp0)
d_ws0=d_ws0-i_fracice*d_salflx*w_sal(:,1)/d_rho(:,1) ! actually salflx*(watersal-icesal)/density(icesal)
d_wm0=d_wm0+i_fracice*d_wtrflx
d_ustar=max(sqrt(sqrt(d_wu0*d_wu0+d_wv0*d_wv0)),1.E-6)
d_b0=-grav*(d_alpha(:,1)*d_wt0-d_beta(:,1)*d_ws0) ! -ve sign is to account for sign of wt0 and ws0

return
end subroutine mloice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Form seaice before flux calculations

subroutine mlonewice(dt,d_rho,d_timelt,d_wm0,d_zcr,diag)

implicit none

integer, intent(in) :: diag
integer ii
real, intent(in) :: dt
real, dimension(wfull) :: newdic,newtn,gamm
real, dimension(wfull,wlev), intent(in) :: d_rho
real, dimension(wfull), intent(inout) :: d_timelt,d_wm0
real, dimension(wfull), intent(in) :: d_zcr
real, dimension(wfull) :: worka

! formation
newdic=max(d_timelt-w_temp(:,1),0.)*cp0*d_rho(:,1)*dz(:,1)*d_zcr/qice
newdic=min(0.15,newdic)
newdic=min(depth_hl(:,wlev+1)-minwater+w_eta,newdic)
newtn=w_temp(:,1)+newdic*qice/(cp0*d_rho(:,1)*dz(:,1)*d_zcr)
where (newdic.gt.1.1*icemin*fracbreak) ! form new sea-ice
  i_dic=i_dic*i_fracice+newdic*(1.-i_fracice)
  i_dsn=i_dsn*i_fracice
  i_tsurf=i_tsurf*i_fracice+newtn*(1.-i_fracice)
  i_tn(:,0)=i_tn(:,0)*i_fracice+newtn*(1.-i_fracice)
  i_tn(:,1)=i_tn(:,1)*i_fracice+newtn*(1.-i_fracice)
  i_tn(:,2)=i_tn(:,2)*i_fracice+newtn*(1.-i_fracice)
  i_sto=i_sto*i_fracice
  w_temp(:,1)=w_temp(:,1)*i_fracice+newtn*(1.-i_fracice)
  d_wm0=d_wm0-newdic*(1.-i_fracice)*rhoic/rhowt/dt
  i_fracice=1. ! this will be adjusted for fracbreak below  
endwhere

! removal
gamm=(gammi*i_dic+gamms*i_dsn)/max(i_dic+i_dsn,1.E-6) ! for energy conservation
where (i_dic.le.icemin)
  d_wm0=d_wm0+i_dic*i_fracice*rhoic/rhowt/dt
  w_temp(:,1)=(w_temp(:,1)+i_fracice*gamm*i_tsurf/(cp0*d_rho(:,1)*dz(:,1)*d_zcr))/(1.+i_fracice*gamm/(cp0*d_rho(:,1)*dz(:,1)*d_zcr))
  w_temp(:,1)=w_temp(:,1)-i_fracice*i_dic*qice/(cp0*d_rho(:,1)*dz(:,1)*d_zcr)
  w_temp(:,1)=w_temp(:,1)-i_fracice*i_dsn*qsnow/(cp0*d_rho(:,1)*dz(:,1)*d_zcr)  
  w_temp(:,1)=w_temp(:,1)+i_fracice*i_sto/(cp0*d_rho(:,1)*dz(:,1)*d_zcr)
  i_fracice=0.
  i_dic=0.
  i_dsn=0.
  i_sto=0.
end where

! 1D model of ice break-up
where (i_dic.lt.icebreak.and.i_fracice.gt.fracbreak)
  i_fracice=i_fracice*i_dic/icebreak
  i_dic=icebreak
end where
where (i_fracice.lt.1..and.i_dic.gt.icebreak)
   worka=min(i_dic/icebreak,1./max(i_fracice,fracbreak))
   worka=max(worka,1.)
   i_fracice=i_fracice*worka
   i_dic=i_dic/worka
end where
i_fracice=min(max(i_fracice,0.),1.)

return
end subroutine mlonewice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update sea ice prognostic variables

subroutine seaicecalc(nice,dt,ip_tn,ip_dic,ip_dsn,ip_fracice,ip_tsurf,ip_sto,ap_rnd,        &
                      ap_snd,pp_egice,dp_ftop,dp_bot,dp_tb,dp_fb,dp_timelt,                 &
                      dp_salflx,dp_wtrflx,dp_nk,dp_tauxica,dp_tauyica,dp_tauxicw,           &
                      dp_tauyicw,diag)

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
real, dimension(nice), intent(inout) :: dp_ftop,dp_bot,dp_tb,dp_fb,dp_timelt,dp_salflx,dp_wtrflx
real, dimension(nice), intent(inout) :: dp_tauxica,dp_tauyica,dp_tauxicw,dp_tauyicw
integer, dimension(nice), intent(inout) :: dp_nk
real, dimension(nice) :: dt_ftop,dt_bot,dt_tb,dt_fb,dt_timelt,dt_salflx,dt_wtrflx
integer, dimension(nice) :: dt_nk
real, dimension(nice), intent(in) :: ap_rnd,ap_snd
real, dimension(nice), intent(in) :: pp_egice
real, dimension(nice) :: pt_egice

! snow fall
ip_dsn=ip_dsn+dt*(ap_rnd+ap_snd)*0.001 ! convert mm to m

! Pack different ice configurations
pq0=dp_nk.gt.0.and.ip_dsn.gt.icemin ! snow + 1-2 ice layers
nc0=count(pq0)
pq1=dp_nk.gt.0.and.ip_dsn.le.icemin ! no snow + 1-2 ice layers
nc1=count(pq1)
pq2=dp_nk.le.0.and.ip_dsn.gt.icemin ! snow + 1 ice layer
nc2=count(pq2)
pq3=dp_nk.le.0.and.ip_dsn.le.icemin ! no snow + 1 ice layer
nc3=count(pq3)

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
  dt_bot(1:nc0)=pack(dp_bot,pq0)
  dt_tb(1:nc0)=pack(dp_tb,pq0)
  dt_fb(1:nc0)=pack(dp_fb,pq0)
  dt_timelt(1:nc0)=pack(dp_timelt,pq0)
  dt_salflx(1:nc0)=pack(dp_salflx,pq0)
  dt_wtrflx(1:nc0)=pack(dp_wtrflx,pq0)
  dt_nk(1:nc0)=pack(dp_nk,pq0)
  call icetemps2(nc0,dt,it_tn(1:nc0,:),it_dic(1:nc0),it_dsn(1:nc0),it_fracice(1:nc0),it_tsurf(1:nc0),it_sto(1:nc0), &
                 dt_ftop(1:nc0),dt_bot(1:nc0),dt_tb(1:nc0),dt_fb(1:nc0),dt_timelt(1:nc0),dt_salflx(1:nc0),          &
                 dt_wtrflx(1:nc0),dt_nk(1:nc0),pt_egice(1:nc0),diag)
  do ii=0,2
    ip_tn(:,ii)=unpack(it_tn(1:nc0,ii),pq0,ip_tn(:,ii))
  end do
  ip_dic=unpack(it_dic(1:nc0),pq0,ip_dic)
  ip_dsn=unpack(it_dsn(1:nc0),pq0,ip_dsn)
  ip_fracice=unpack(it_fracice(1:nc0),pq0,ip_fracice)
  ip_tsurf=unpack(it_tsurf(1:nc0),pq0,ip_tsurf)
  ip_sto=unpack(it_sto(1:nc0),pq0,ip_sto)
  dp_ftop=unpack(dt_ftop(1:nc0),pq0,dp_ftop)
  dp_bot=unpack(dt_bot(1:nc0),pq0,dp_bot)
  dp_tb=unpack(dt_tb(1:nc0),pq0,dp_tb)
  dp_fb=unpack(dt_fb(1:nc0),pq0,dp_fb)
  dp_timelt=unpack(dt_timelt(1:nc0),pq0,dp_timelt)
  dp_salflx=unpack(dt_salflx(1:nc0),pq0,dp_salflx)
  dp_wtrflx=unpack(dt_wtrflx(1:nc0),pq0,dp_wtrflx)
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
  dt_bot(1:nc1)=pack(dp_bot,pq1)
  dt_tb(1:nc1)=pack(dp_tb,pq1)
  dt_fb(1:nc1)=pack(dp_fb,pq1)
  dt_timelt(1:nc1)=pack(dp_timelt,pq1)
  dt_salflx(1:nc1)=pack(dp_salflx,pq1)
  dt_wtrflx(1:nc1)=pack(dp_wtrflx,pq1)
  dt_nk(1:nc1)=pack(dp_nk,pq1)
  call icetempn2(nc1,dt,it_tn(1:nc1,:),it_dic(1:nc1),it_dsn(1:nc1),it_fracice(1:nc1),it_tsurf(1:nc1),it_sto(1:nc1), &
                 dt_ftop(1:nc1),dt_bot(1:nc1),dt_tb(1:nc1),dt_fb(1:nc1),dt_timelt(1:nc1),dt_salflx(1:nc1),          &
                 dt_wtrflx(1:nc1),dt_nk(1:nc1),pt_egice(1:nc1),diag)
  do ii=0,2
    ip_tn(:,ii)=unpack(it_tn(1:nc1,ii),pq1,ip_tn(:,ii))
  end do
  ip_dic=unpack(it_dic(1:nc1),pq1,ip_dic)
  ip_dsn=unpack(it_dsn(1:nc1),pq1,ip_dsn)
  ip_fracice=unpack(it_fracice(1:nc1),pq1,ip_fracice)
  ip_tsurf=unpack(it_tsurf(1:nc1),pq1,ip_tsurf)
  ip_sto=unpack(it_sto(1:nc1),pq1,ip_sto)
  dp_ftop=unpack(dt_ftop(1:nc1),pq1,dp_ftop)
  dp_bot=unpack(dt_bot(1:nc1),pq1,dp_bot)
  dp_tb=unpack(dt_tb(1:nc1),pq1,dp_tb)
  dp_fb=unpack(dt_fb(1:nc1),pq1,dp_fb)
  dp_timelt=unpack(dt_timelt(1:nc1),pq1,dp_timelt)
  dp_salflx=unpack(dt_salflx(1:nc1),pq1,dp_salflx)
  dp_wtrflx=unpack(dt_wtrflx(1:nc1),pq1,dp_wtrflx)
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
  dt_bot(1:nc2)=pack(dp_bot,pq2)
  dt_tb(1:nc2)=pack(dp_tb,pq2)
  dt_fb(1:nc2)=pack(dp_fb,pq2)
  dt_timelt(1:nc2)=pack(dp_timelt,pq2)
  dt_salflx(1:nc2)=pack(dp_salflx,pq2)
  dt_wtrflx(1:nc2)=pack(dp_wtrflx,pq2)
  dt_nk(1:nc2)=pack(dp_nk,pq2)
  call icetemps1(nc2,dt,it_tn(1:nc2,:),it_dic(1:nc2),it_dsn(1:nc2),it_fracice(1:nc2),it_tsurf(1:nc2),it_sto(1:nc2), &
                 dt_ftop(1:nc2),dt_bot(1:nc2),dt_tb(1:nc2),dt_fb(1:nc2),dt_timelt(1:nc2),dt_salflx(1:nc2),          &
                 dt_wtrflx(1:nc2),dt_nk(1:nc2),pt_egice(1:nc2),diag)
  do ii=0,2
    ip_tn(:,ii)=unpack(it_tn(1:nc2,ii),pq2,ip_tn(:,ii))
  end do
  ip_dic=unpack(it_dic(1:nc2),pq2,ip_dic)
  ip_dsn=unpack(it_dsn(1:nc2),pq2,ip_dsn)
  ip_fracice=unpack(it_fracice(1:nc2),pq2,ip_fracice)
  ip_tsurf=unpack(it_tsurf(1:nc2),pq2,ip_tsurf)
  ip_sto=unpack(it_sto(1:nc2),pq2,ip_sto)
  dp_ftop=unpack(dt_ftop(1:nc2),pq2,dp_ftop)
  dp_bot=unpack(dt_bot(1:nc2),pq2,dp_bot)
  dp_tb=unpack(dt_tb(1:nc2),pq2,dp_tb)
  dp_fb=unpack(dt_fb(1:nc2),pq2,dp_fb)
  dp_timelt=unpack(dt_timelt(1:nc2),pq2,dp_timelt)
  dp_salflx=unpack(dt_salflx(1:nc2),pq2,dp_salflx)
  dp_wtrflx=unpack(dt_wtrflx(1:nc2),pq2,dp_wtrflx)
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
  dt_bot(1:nc3)=pack(dp_bot,pq3)
  dt_tb(1:nc3)=pack(dp_tb,pq3)
  dt_fb(1:nc3)=pack(dp_fb,pq3)
  dt_timelt(1:nc3)=pack(dp_timelt,pq3)
  dt_salflx(1:nc3)=pack(dp_salflx,pq3)
  dt_wtrflx(1:nc3)=pack(dp_wtrflx,pq3)
  dt_nk(1:nc3)=pack(dp_nk,pq3)
  call icetempn1(nc3,dt,it_tn(1:nc3,:),it_dic(1:nc3),it_dsn(1:nc3),it_fracice(1:nc3),it_tsurf(1:nc3),it_sto(1:nc3), &
                 dt_ftop(1:nc3),dt_bot(1:nc3),dt_tb(1:nc3),dt_fb(1:nc3),dt_timelt(1:nc3),dt_salflx(1:nc3),          &
                 dt_wtrflx(1:nc3),dt_nk(1:nc3),pt_egice(1:nc3),diag)
  do ii=0,2
    ip_tn(:,ii)=unpack(it_tn(1:nc3,ii),pq3,ip_tn(:,ii))
  end do
  ip_dic=unpack(it_dic(1:nc3),pq3,ip_dic)
  ip_dsn=unpack(it_dsn(1:nc3),pq3,ip_dsn)
  ip_fracice=unpack(it_fracice(1:nc3),pq3,ip_fracice)
  ip_tsurf=unpack(it_tsurf(1:nc3),pq3,ip_tsurf)
  ip_sto=unpack(it_sto(1:nc3),pq3,ip_sto)
  dp_ftop=unpack(dt_ftop(1:nc3),pq3,dp_ftop)
  dp_bot=unpack(dt_bot(1:nc3),pq3,dp_bot)
  dp_tb=unpack(dt_tb(1:nc3),pq3,dp_tb)
  dp_fb=unpack(dt_fb(1:nc3),pq3,dp_fb)
  dp_timelt=unpack(dt_timelt(1:nc3),pq3,dp_timelt)
  dp_salflx=unpack(dt_salflx(1:nc3),pq3,dp_salflx)
  dp_wtrflx=unpack(dt_wtrflx(1:nc3),pq3,dp_wtrflx)
  dp_nk=unpack(dt_nk(1:nc3),pq3,dp_nk)
end if

! the following are snow to ice processes
where (ip_dic.gt.icemin)
  xxx=ip_dic+ip_dsn-(rhosn*ip_dsn+rhoic*ip_dic)/rhowt ! white ice formation
  excess=max(ip_dsn-xxx,0.)*rhosn/rhowt               ! white ice formation
elsewhere
  xxx=0.
  excess=0.
end where
where (excess.gt.0.and.dp_nk.eq.2)
  ip_tn(:,1)=(0.5*ip_dic*ip_tn(:,1)+rhowt/rhoic*excess*ip_tn(:,0)*cps/cpi)/(0.5*ip_dic+rhowt/rhoic*excess) ! Assume 2 levels of ice, hence 0.5*hi
elsewhere (excess.gt.0..and.dp_nk.eq.1)
  ip_tn(:,1)=(ip_dic*ip_tn(:,1)+rhowt/rhoic*excess*ip_tn(:,0)*cps/cpi)/(ip_dic+rhowt/rhoic*excess)         ! Assume 1 level of ice
end where
excess=excess+max(ip_dsn-0.2,0.)*rhosn/rhowt        ! Snow depth limitation and conversion to ice
dp_salflx=dp_salflx-rhowt*excess/dt
ip_dsn=min(min(xxx,ip_dsn),0.2)
ip_dic=ip_dic+rhowt/rhoic*excess

! Ice depth limitation for poor initial conditions
xxx=max(ip_dic-icemax,0.)
ip_dic=min(icemax,ip_dic)
dp_wtrflx=dp_wtrflx+xxx*rhoic/rhowt/dt

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

subroutine icetemps2(nc,dt,it_tn,it_dic,it_dsn,it_fracice,it_tsurf,it_sto,dt_ftop,dt_bot,dt_tb,dt_fb,dt_timelt,dt_salflx, &
                     dt_wtrflx,dt_nk,pt_egice,diag)

implicit none

integer, intent(in) :: nc,diag
integer ii,iqi
real, intent(in) :: dt
real qneed,fact
real htup,htdown
real, dimension(nc) :: con,rhin,fs,conb
real, dimension(nc) :: subl,snmelt,dhs,dhb,ssubl,ssnmelt
real, dimension(nc,0:2) :: fl
real, dimension(nc,0:2), intent(inout) :: it_tn
real, dimension(nc), intent(inout) :: it_dic,it_dsn,it_fracice,it_tsurf,it_sto
real, dimension(nc), intent(inout) :: dt_ftop,dt_bot,dt_tb,dt_fb,dt_timelt,dt_salflx,dt_wtrflx
integer, dimension(nc), intent(inout) :: dt_nk
real, dimension(nc), intent(in) :: pt_egice

! update surface temperature
rhin=real(dt_nk)/it_dic
con=2.*condsnw/max(it_dsn,icemin)
! Ammendments if surface layer < 5cms
where (it_dsn.lt.0.05)
  con=1./(it_dsn/condsnw+0.5/(condice*rhin))
  it_tn(:,0)=it_tn(:,1)
end where
fs=con*(it_tn(:,0)-it_tsurf)
it_tsurf=it_tsurf+(dt_ftop+fs)/(dt_bot+gamms/dt+con)
! Snow melt (snmelt >= 0)
snmelt=max(0.,it_tsurf-273.16)*gamms/qsnow ! m of "snow"
it_tsurf=min(it_tsurf,273.16)              ! melting condition ts=tsmelt

!ti=(it_dic*condsnw*it_tn(:,0)+real(dt_nk)*it_dsn*condice*it_tn(:,1)) &
!  /(it_dic*condsnw+real(dt_nk)*it_dsn*condice)
!fl(:,0)=2.*rhin*condice*(it_tn(:,1)-ti)      ! Middle of first ice layer to snow interface
conb=2./(it_dsn/condsnw+1./(condice*rhin))
fl(:,0)=conb*(it_tn(:,1)-it_tn(:,0))
! Ammendments if surface layer < 5cms
where (it_dsn.lt.0.05)
  fl(:,0)=fs
end where
fl(:,1)=rhin*condice*(it_tn(:,2)-it_tn(:,1)) ! Between ice layers 2 and 1

! Update the mid-level snow temperature
where (it_dsn.ge.0.05)
  it_tn(:,0)=it_tn(:,0)+dt*(fl(:,0)-fs)/(it_dsn*cps)
end where

! Surface evap/sublimation (can be >0 or <0)
! Note : dt*eg/hl in Kgm/m**2 => mms of water
subl=dt*pt_egice*lf/(lv*qsnow)      ! m of "snow"
ssubl=min(subl,it_dsn)           ! snow component of sublimation
ssnmelt=min(snmelt,it_dsn-ssubl) ! snow component of melt
dt_salflx=dt_salflx-ssnmelt*rhosn/dt ! melt fresh water snow (no salt from egice when melting snow)
dt_wtrflx=dt_wtrflx+snmelt*rhosn/rhowt/dt
! Change the snow thickness
dhs=snmelt+subl
dhb=(rhosn/rhoic)*(subl-ssubl+snmelt-ssnmelt)
it_dic=max(it_dic-dhb,0.)
it_dsn=max(it_dsn-dhs,0.)

do iqi=1,nc
  if (dt_nk(iqi).eq.2) then
    ! update temperature in top layer of ice
    it_tn(iqi,1)=it_tn(iqi,1)+dt*(fl(iqi,1)-fl(iqi,0))*rhin(iqi)/cpi
    ! use stored heat in brine pockets to keep temperature at -0.1 until heat is used up
    if (it_sto(iqi).gt.0..and.it_tn(iqi,1).lt.dt_timelt(iqi)) then
      qneed=(dt_timelt(iqi)-it_tn(iqi,1))*cpi/rhin(iqi) ! J/m**2
      if (it_sto(iqi).lt.qneed) then
        fact=it_sto(iqi)/qneed
        fact=min(max(fact,0.),1.)
        it_tn(iqi,1)=dt_timelt(iqi)*fact+it_tn(iqi,1)*(1.-fact)
        it_sto(iqi)=0.
      else
        it_tn(iqi,1)=dt_timelt(iqi)
        it_sto(iqi)=it_sto(iqi)-qneed
      end if
    end if
  else
    dt_wtrflx(iqi)=dt_wtrflx(iqi)+min(it_sto(iqi)/qice,it_dic(iqi))*rhoic/rhowt/dt
    it_dic(iqi)=max(0.,it_dic(iqi)-it_sto(iqi)/qice)
    it_sto(iqi)=0.
  end if

  ! determine amount of bottom ablation or accretion
  fl(iqi,dt_nk(iqi))=condice*(dt_tb(iqi)-it_tn(iqi,dt_nk(iqi)))*2.*rhin(iqi)
  fl(iqi,dt_nk(iqi))=min(max(fl(iqi,dt_nk(iqi)),-1000.),1000.)
  dhb(iqi)=dt*(fl(iqi,dt_nk(iqi))-dt_fb(iqi))/qice
  dt_wtrflx(iqi)=dt_wtrflx(iqi)+min(-dhb(iqi),it_dic(iqi))*rhoic/rhowt/dt
  it_dic(iqi)=max(0.,it_dic(iqi)+dhb(iqi))
  it_tn(iqi,dt_nk(iqi))=it_tn(iqi,dt_nk(iqi))+dt/cpi*(fl(iqi,dt_nk(iqi))-fl(iqi,dt_nk(iqi)-1))*rhin(iqi)

  ! test whether to change number of layers
  htup=real(dt_nk(iqi)+1)*himin
  htdown=real(dt_nk(iqi))*himin
  if (it_dic(iqi).ge.htup.and.dt_nk(iqi).lt.2) then
    ! code to increase number of layers
    dt_nk(iqi)=2
    it_tn(iqi,2)=it_tn(iqi,1)
    ! snow depth has not changed, so split ice layer into two
  elseif (it_dic(iqi).le.htdown) then
    ! code to decrease number of layers
    dt_nk(iqi)=dt_nk(iqi)-1
    if (dt_nk(iqi).eq.1) then
      ! merge ice layers into one
      it_tn(iqi,1)=0.5*sum(it_tn(iqi,1:2))
      it_tn(iqi,2)=it_tn(iqi,1)
    else
      it_tsurf(iqi)=(it_tsurf(iqi)+cpi*it_dic(iqi)*it_tn(iqi,1)/gammi)/(1.+cpi*it_dic(iqi)/gammi)
      it_tsurf(iqi)=(it_tsurf(iqi)+cps*it_dsn(iqi)*it_tn(iqi,0)/gammi)/(1.+cps*it_dsn(iqi)/gammi)
      it_tn(iqi,1)=it_tsurf(iqi)
      it_tn(iqi,2)=it_tsurf(iqi)
    end if
  end if
end do

where (dt_nk.eq.1)
  it_tn(:,2)=it_tn(:,1)
end where

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

subroutine icetempn2(nc,dt,it_tn,it_dic,it_dsn,it_fracice,it_tsurf,it_sto,dt_ftop,dt_bot,dt_tb,dt_fb,dt_timelt,dt_salflx, &
                     dt_wtrflx,dt_nk,pt_egice,diag)

implicit none

integer, intent(in) :: nc,diag
integer ii,iqi
real, intent(in) :: dt
real qneed,fact
real htup,htdown
real, dimension(nc) :: con,rhin,qmax,ftopadj
real, dimension(nc) :: subl,simelt,dhi,dhb
real, dimension(nc,0:2) :: fl
real, dimension(nc,0:2), intent(inout) :: it_tn
real, dimension(nc), intent(inout) :: it_dic,it_dsn,it_fracice,it_tsurf,it_sto
real, dimension(nc), intent(inout) :: dt_ftop,dt_bot,dt_tb,dt_fb,dt_timelt,dt_salflx,dt_wtrflx
integer, dimension(nc), intent(inout) :: dt_nk
real, dimension(nc), intent(in) :: pt_egice

! Get total flux into thin surface ice layer by adding flux from below
rhin=real(dt_nk)/it_dic
con=condice*rhin*2.
ftopadj=con*(it_tn(:,1)-it_tsurf)
it_tsurf=it_tsurf+(dt_ftop+ftopadj)/(dt_bot+gammi/dt+con)
! Ice melt (simelt >= 0)
simelt=max(0.,it_tsurf-dt_timelt)*gammi/qice
it_tsurf=min(it_tsurf,dt_timelt) ! melting condition

! Upward fluxes between various levels below surface
fl(:,0)=ftopadj                              ! Middle of first ice layer to thin surface layer
fl(:,1)=condice*(it_tn(:,2)-it_tn(:,1))*rhin ! Between ice layers 2 and 1

! Surface evap/sublimation (can be >0 or <0)
subl=dt*pt_egice*lf/(lv*qice)
dt_salflx=dt_salflx+pt_egice/lv
dhi=-subl-simelt
qmax=qice*0.5*max(it_dic-himin,0.)
dhi=dhi-max(it_sto-qmax,0.)/qice
it_sto=min(qmax,it_sto)
where (dhi.lt.0..and.dt_nk.lt.2)
  ! if ice is melting, update thickness when single level
  dt_wtrflx=dt_wtrflx+min(-dhi,it_dic)*rhoic/rhowt/dt
  it_dic=max(0.,it_dic+dhi)
  dhi=0.
  rhin=1./max(it_dic,icemin)
end where

do iqi=1,nc
  if (dt_nk(iqi).eq.2) then
    it_tn(iqi,1)=it_tn(iqi,1)+dt/cpi*(fl(iqi,1)-fl(iqi,0))*rhin(iqi)
    if (it_sto(iqi).gt.0..and.it_tn(iqi,1).lt.dt_timelt(iqi)) then
      qneed=(dt_timelt(iqi)-it_tn(iqi,1))*cpi/rhin(iqi) ! J/m**2
      if (it_sto(iqi).lt.qneed) then
        fact=it_sto(iqi)/qneed
        fact=min(max(fact,0.),1.)
        it_tn(iqi,1)=dt_timelt(iqi)*fact+it_tn(iqi,1)*(1.-fact)
        it_sto(iqi)=0.
      else
        it_tn(iqi,1)=dt_timelt(iqi)
        it_sto(iqi)=it_sto(iqi)-qneed
      end if
    end if
  else
    dt_wtrflx(iqi)=dt_wtrflx(iqi)+min(it_sto(iqi)/qice,it_dic(iqi))*rhoic/rhowt/dt  
    it_dic(iqi)=max(0.,it_dic(iqi)-it_sto(iqi)/qice)
    it_sto(iqi)=0.
  end if

  ! determine amount of bottom ablation or accretion
  fl(iqi,dt_nk(iqi))=condice*(dt_tb(iqi)-it_tn(iqi,dt_nk(iqi)))*2.*rhin(iqi)
  fl(iqi,dt_nk(iqi))=min(max(fl(iqi,dt_nk(iqi)),-1000.),1000.)
  dhb(iqi)=dt*(fl(iqi,dt_nk(iqi))-dt_fb(iqi))/qice
  dt_wtrflx(iqi)=dt_wtrflx(iqi)+min(-dhi(iqi)-dhb(iqi),it_dic(iqi))*rhoic/rhowt/dt
  it_dic(iqi)=max(0.,it_dic(iqi)+dhi(iqi)+dhb(iqi))
  it_tn(iqi,dt_nk(iqi))=it_tn(iqi,dt_nk(iqi))+dt/cpi*(fl(iqi,dt_nk(iqi))-fl(iqi,dt_nk(iqi)-1))*rhin(iqi)

  ! test whether to change number of layers
  htup=real(dt_nk(iqi)+1)*himin
  htdown=real(dt_nk(iqi))*himin
  if (it_dic(iqi).ge.htup.and.dt_nk(iqi).lt.2) then
    ! code to increase number of layers
    dt_nk(iqi)=2
    it_tn(iqi,2)=it_tn(iqi,1)
    ! split ice layer into two
  elseif (it_dic(iqi).le.htdown) then
    ! code to decrease number of layers
    dt_nk(iqi)=dt_nk(iqi)-1
    if (dt_nk(iqi).eq.1) then
      ! merge ice layers into one
      it_tn(iqi,1)=0.5*sum(it_tn(iqi,1:2))
      it_tn(iqi,2)=it_tn(iqi,1)
    else
      it_tsurf(iqi)=(it_tsurf(iqi)+cpi*it_dic(iqi)*it_tn(iqi,1)/gammi)/(1.+cpi*it_dic(iqi)/gammi)
      it_tn(iqi,1)=it_tsurf(iqi)
      it_tn(iqi,2)=it_tsurf(iqi)
    end if
  end if

end do
it_tn(:,0)=it_tn(:,1)
where (dt_nk.eq.1)
  it_tn(:,2)=it_tn(:,1)
end where

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

subroutine icetemps1(nc,dt,it_tn,it_dic,it_dsn,it_fracice,it_tsurf,it_sto,dt_ftop,dt_bot,dt_tb,dt_fb,dt_timelt,dt_salflx, &
                     dt_wtrflx,dt_nk,pt_egice,diag)

implicit none

integer, intent(in) :: nc,diag
real, intent(in) :: dt
real, dimension(nc) :: con,f0
real, dimension(nc) :: subl,snmelt,dhs,dhb,ssubl,ssnmelt,gamm
real, dimension(nc,0:2), intent(inout) :: it_tn
real, dimension(nc), intent(inout) :: it_dic,it_dsn,it_fracice,it_tsurf,it_sto
real, dimension(nc), intent(inout) :: dt_ftop,dt_bot,dt_tb,dt_fb,dt_timelt,dt_salflx,dt_wtrflx
integer, dimension(nc), intent(inout) :: dt_nk
real, dimension(nc), intent(in) :: pt_egice

! Update tsurf and ti based on fluxes from above and below
con=1./(it_dsn/condsnw+max(it_dic,icemin)/condice)
f0=con*(dt_tb-it_tsurf) ! flux from below
f0=min(max(f0,-1000.),1000.)
gamm=(gammi*it_dic+gamms*it_dsn)/(it_dic+it_dsn) ! for energy conservation
it_tsurf=it_tsurf+(dt_ftop+f0)/(dt_bot+gamm/dt+con)
! Snow melt (snmelt >= 0)
snmelt=max(0.,it_tsurf-273.16)*gamm/qsnow ! note change from gamms to gamm
it_tsurf=min(it_tsurf,273.16)             ! melting condition ts=tsmelt

! Surface evap/sublimation (can be >0 or <0)
subl=dt*pt_egice*lf/(lv*qsnow)
ssubl=min(subl,it_dsn)           ! snow component of sublimation
ssnmelt=min(snmelt,it_dsn-ssubl) ! snow component of melt
dt_salflx=dt_salflx-ssnmelt*rhosn/dt
dt_wtrflx=dt_wtrflx+snmelt*rhosn/rhowt/dt
dhs=ssnmelt+ssubl ! Change the snow thickness
it_dic=it_dic-(rhosn/rhoic)*(snmelt-ssnmelt+subl-ssubl)
dhb=dt*(f0-dt_fb)/qice       ! Ice melt
dt_wtrflx=dt_wtrflx+min(-dhb,it_dic)*rhoic/rhowt/dt
it_dsn=max(it_dsn-dhs,0.)
it_dic=max(it_dic+dhb,0.)

where (it_dic.ge.himin) ! increase ice thickness
  dt_nk=1
end where

it_tn(:,0)=it_tsurf
it_tn(:,1)=it_tsurf
it_tn(:,2)=it_tsurf

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

subroutine icetempn1(nc,dt,it_tn,it_dic,it_dsn,it_fracice,it_tsurf,it_sto,dt_ftop,dt_bot,dt_tb,dt_fb,dt_timelt,dt_salflx, &
                     dt_wtrflx,dt_nk,pt_egice,diag)

implicit none

integer, intent(in) :: nc,diag
real, intent(in) :: dt
real, dimension(nc) :: con,subl,simelt
real, dimension(nc) :: dhi,dhb,f0
real, dimension(nc,0:2), intent(inout) :: it_tn
real, dimension(nc), intent(inout) :: it_dic,it_dsn,it_fracice,it_tsurf,it_sto
real, dimension(nc), intent(inout) :: dt_ftop,dt_bot,dt_tb,dt_fb,dt_timelt,dt_salflx,dt_wtrflx
integer, dimension(nc), intent(inout) :: dt_nk
real, dimension(nc), intent(in) :: pt_egice

! Update tsurf and ti based on fluxes from above and below
con=condice/max(it_dic,icemin)
f0=con*(dt_tb-it_tsurf) ! flux from below
f0=min(max(f0,-1000.),1000.)
it_tsurf=it_tsurf+(dt_ftop+f0)/(dt_bot+gammi/dt+con)
simelt=max(0.,it_tsurf-dt_timelt)*gammi/qice ! Ice melt (simelt >= 0)
it_tsurf=min(it_tsurf,dt_timelt)             ! melting condition ti=timelt

! Surface evap/sublimation (can be >0 or <0)
subl=dt*pt_egice*lf/(lv*qice)
dt_salflx=dt_salflx+pt_egice/lv
dhi=-subl-simelt
dhb=dt*(f0-dt_fb)/qice ! determine amount of bottom ablation or accretion
dt_wtrflx=dt_wtrflx+min(-dhi-dhb,it_dic)*rhoic/rhowt/dt
it_dic=max(it_dic+dhi+dhb,0.)  ! update ice thickness

where (it_dic.ge.himin) ! increase ice thickness
  dt_nk=1
end where

it_tn(:,0)=it_tsurf
it_tn(:,1)=it_tsurf
it_tn(:,2)=it_tsurf

return
end subroutine icetempn1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Determine ice fluxes

subroutine iceflux(dt,a_sg,a_rg,a_rnd,a_snd,a_vnratio,a_fbvis,a_fbnir,a_u,a_v,a_temp,a_qg,a_ps,a_zmin,a_zmins, &
                   d_rho,d_ftop,d_bot,d_tb,d_fb,d_timelt,d_nk,d_tauxica,d_tauyica,d_tauxicw,d_tauyicw,d_zcr,   &
                   diag)

implicit none

integer, intent(in) :: diag
integer ii,ll,iqw
real, dimension(wfull), intent(in) :: a_sg,a_rg,a_rnd,a_snd,a_vnratio,a_fbvis,a_fbnir,a_u,a_v,a_temp,a_qg,a_ps,a_zmin,a_zmins
real, dimension(wfull,wlev), intent(in) :: d_rho
real, dimension(wfull), intent(inout) :: d_ftop,d_bot,d_tb,d_fb,d_timelt,d_tauxicw,d_tauyicw,d_tauxica,d_tauyica,d_zcr
integer, dimension(wfull), intent(inout) :: d_nk
real, intent(in) :: dt
real, dimension(wfull) :: qsat,dqdt,ri,vmag,rho,srcp
real, dimension(wfull) :: fm,fh,af,aft
real, dimension(wfull) :: den,sig,root
real, dimension(wfull) :: alb,qmax,eye
real, dimension(wfull) :: uu,vv,du,dv,vmagn,icemagn
real, dimension(wfull) :: x,ustar,icemag,g,h,dgu,dgv,dhu,dhv,det
real, dimension(wfull) :: newiu,newiv,imass
real factch

uu=a_u-i_u
vv=a_v-i_v
vmag=sqrt(uu*uu+vv*vv)
vmagn=max(vmag,0.2)
sig=exp(-grav*a_zmins/(rdry*a_temp))
srcp=sig**(rdry/cp)
rho=a_ps/(rdry*i_tsurf)

p_zoice=0.001
p_zohice=0.001
p_zoqice=0.001
af=vkar**2/(log(a_zmin/p_zoice)*log(a_zmin/p_zoice))
aft=vkar**2/(log(a_zmins/p_zoice)*log(a_zmins/p_zohice))
factch=1.

call getqsat(qsat,dqdt,i_tsurf,a_ps)
ri=min(grav*(a_zmin**2/a_zmins)*(1.-i_tsurf*srcp/a_temp)/vmagn**2,rimax)

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

! egice is for evaporating (lv).  Melting is included above with lf.

p_wetfacice=max(1.+.008*min(i_tsurf-273.16,0.),0.)
p_cdice=af*fm
p_cdsice=aft*fh
p_fgice=rho*aft*cp*fh*vmag*(i_tsurf-a_temp/srcp)
p_egice=p_wetfacice*rho*aft*lv*fh*vmag*(qsat-a_qg)

! water temperature at bottom of ice
d_tb=w_temp(:,1)

! iterative method to estimate ice velocity after stress tensors are applied
imass=max(rhoic*i_dic+rhosn*i_dsn,10.) ! ice mass per unit area
newiu=i_u
newiv=i_v
do ll=1,5
  uu=a_u-newiu
  vv=a_v-newiv
  du=w_u(:,1)-newiu
  dv=w_v(:,1)-newiv
  vmagn=max(sqrt(uu*uu+vv*vv),0.2)
  icemagn=max(sqrt(du*du+dv*dv),0.0002)

  g=i_u-newiu+dt*(rho*p_cdice*vmagn*uu+d_rho(:,1)*0.00536*icemagn*du)/imass
  h=i_v-newiv+dt*(rho*p_cdice*vmagn*vv+d_rho(:,1)*0.00536*icemagn*dv)/imass
  
  dgu=-1.-dt*(rho*p_cdice*vmagn*(1.+(uu/vmagn)**2)+d_rho(:,1)*0.00536*icemagn*(1.+(du/icemagn)**2))/imass
  dhu=-dt*(rho*p_cdice*uu*vv/vmagn+d_rho(:,1)*0.00536*du*dv/icemagn)/imass
  dgv=dhu
  dhv=-1.-dt*(rho*p_cdice*vmagn*(1.+(vv/vmagn)**2)+d_rho(:,1)*0.00536*icemagn*(1.+(dv/icemagn)**2))/imass

  det=dgu*dhv-dgv*dhu
  where (abs(det).gt.1.E-20)
    newiu=newiu-0.9*( g*dhv-h*dgv)/det
    newiv=newiv-0.9*(-g*dhu+h*dgu)/det
  end where

  newiu=max(newiu,min(a_u,w_u(:,1),i_u))
  newiu=min(newiu,max(a_u,w_u(:,1),i_u))
  newiv=max(newiv,min(a_v,w_v(:,1),i_v))
  newiv=min(newiv,max(a_v,w_v(:,1),i_v))
end do

! momentum transfer
d_tauxica=rho*p_cdice*vmagn*uu
d_tauyica=rho*p_cdice*vmagn*vv
d_tauxicw=-d_rho(:,1)*0.00536*icemagn*du
d_tauyicw=-d_rho(:,1)*0.00536*icemagn*dv
ustar=sqrt(sqrt(d_tauxicw*d_tauxicw+d_tauyicw*d_tauyicw)/d_rho(:,1))
ustar=max(ustar,5.E-4)
d_fb=cp0*d_rho(:,1)*0.006*ustar*(d_tb-d_timelt)
d_fb=min(max(d_fb,-1000.),1000.)  

! update ice velocity
i_u=newiu !=i_u+dt*(d_tauxica-d_tauxicw)/imass
i_v=newiv !=i_v+dt*(d_tauyica-d_tauyicw)/imass

! index of different ice thickness configurations
d_nk=min(int(i_dic/himin),2)

! radiation
alb=     a_vnratio*(p_icevisdiralb*a_fbvis+p_icevisdifalb*(1.-a_fbvis))+ &
    (1.-a_vnratio)*(p_icevisdifalb*a_fbvis+p_icevisdifalb*(1.-a_fbvis))
qmax=qice*0.5*max(i_dic-himin,0.)
eye=0.
where (i_dsn.lt.icemin.and.i_sto.lt.qmax.and.d_nk.gt.0)
  eye=0.35
end where
i_sto=i_sto+dt*a_sg*(1.-alb)*eye
d_ftop=-p_fgice-p_egice*ls/lv+a_rg-emisice*sbconst*i_tsurf**4+a_sg*(1.-alb)*(1.-eye)
d_bot=4.*emisice*sbconst*i_tsurf**3+rho*aft*fh*vmag*(cp+p_wetfacice*ls*dqdt)

! Add flux of heat due to converting any rain to snowfall over ice
d_ftop=d_ftop+lf*a_rnd ! rain (mm) to W/m**2

return
end subroutine iceflux

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Determine screen diagnostics

subroutine scrncalc(a_u,a_v,a_temp,a_qg,a_ps,a_zmin,a_zmins,diag)

implicit none

integer, intent(in) :: diag
real, dimension(wfull), intent(in) :: a_u,a_v,a_temp,a_qg,a_ps,a_zmin,a_zmins
real, dimension(wfull) :: tscrn,qgscrn,uscrn,u10
real, dimension(wfull) :: smixr,qsat,dqdt,atu,atv

! water
call getqsat(qsat,dqdt,w_temp(:,1),a_ps)
if (zomode.eq.0) then
  smixr=0.98*qsat
else
  smixr=qsat
end if
atu=a_u-w_u(:,1)
atv=a_v-w_v(:,1)
call scrntile(tscrn,qgscrn,uscrn,u10,p_zo,p_zoh,p_zoq,w_temp(:,1),smixr,atu,atv,a_temp,a_qg,a_zmin,a_zmins,diag)
p_tscrn=tscrn
p_qgscrn=qgscrn
p_uscrn=uscrn+sqrt(w_u(:,1)*w_u(:,1)+w_v(:,1)*w_v(:,1))
p_u10=u10+sqrt(w_u(:,1)*w_u(:,1)+w_v(:,1)*w_v(:,1))

! ice
call getqsat(qsat,dqdt,i_tsurf,a_ps)
smixr=p_wetfacice*qsat+(1.-p_wetfacice)*min(qsat,a_qg)
atu=a_u-i_u
atv=a_v-i_v
call scrntile(tscrn,qgscrn,uscrn,u10,p_zoice,p_zohice,p_zoqice,i_tsurf,smixr,atu,atv,a_temp,a_qg,a_zmin,a_zmins,diag)
p_tscrn=(1.-i_fracice)*p_tscrn+i_fracice*tscrn
p_qgscrn=(1.-i_fracice)*p_qgscrn+i_fracice*qgscrn
p_uscrn=(1.-i_fracice)*p_uscrn+i_fracice*(uscrn+sqrt(i_u*i_u+i_v*i_v))
p_u10=(1.-i_fracice)*p_u10+i_fracice*(u10+sqrt(i_u*i_u+i_v*i_v))

return
end subroutine scrncalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! screen diagnostic for individual tile

subroutine scrntile(tscrn,qgscrn,uscrn,u10,zo,zoh,zoq,stemp,smixr,a_u,a_v,a_temp,a_qg,a_zmin,a_zmins,diag)
      
implicit none

integer, intent(in) :: diag
integer ic
real, dimension(wfull), intent(in) :: a_u,a_v,a_temp,a_qg,a_zmin,a_zmins
real, dimension(wfull), intent(in) :: zo,zoh,zoq,stemp,smixr
real, dimension(wfull), intent(out) :: tscrn,qgscrn,uscrn,u10
real, dimension(wfull) :: umag,sig
real, dimension(wfull) :: lzom,lzoh,lzoq,thetav,sthetav
real, dimension(wfull) :: thetavstar,z_on_l,zs_on_l,z0_on_l,z0s_on_l,zt_on_l,zq_on_l
real, dimension(wfull) :: pm0,ph0,pq0,pm1,ph1,integralm,integralh,integralq
real, dimension(wfull) :: ustar,qstar,z10_on_l
real, dimension(wfull) :: neutrals,neutral,neutral10,pm10
real, dimension(wfull) :: integralm10,tstar,scrp
integer, parameter ::  nc     = 5
real, parameter    ::  a_1    = 1.
real, parameter    ::  b_1    = 2./3.
real, parameter    ::  c_1    = 5.
real, parameter    ::  d_1    = 0.35
real, parameter    ::  z0     = 1.5
real, parameter    ::  z10    = 10.

umag=max(sqrt(a_u*a_u+a_v*a_v),0.2)
sig=exp(-grav*a_zmins/(rdry*a_temp))
scrp=sig**(rdry/cp)
thetav=a_temp*(1.+0.61*a_qg)/scrp
sthetav=stemp*(1.+0.61*smixr)

! Roughness length for heat
lzom=log(a_zmin/zo)
lzoh=log(a_zmins/zoh)
lzoq=log(a_zmins/zoq)

! Dyer and Hicks approach 
thetavstar=vkar*(thetav-sthetav)/lzoh
ustar=vkar*umag/lzom
do ic=1,nc
  z_on_l=a_zmin*vkar*grav*thetavstar/(thetav*ustar**2)
  z_on_l=min(z_on_l,10.)
  zs_on_l  = z_on_l*a_zmins/a_zmin  
  z0_on_l  = z_on_l*zo/a_zmin
  zt_on_l  = z_on_l*zoh/a_zmin
  zq_on_l  = z_on_l*zoq/a_zmin
  where (z_on_l.lt.0.)
    pm0     = (1.-16.*z0_on_l)**(-0.25)
    ph0     = (1.-16.*zt_on_l)**(-0.5)
    pq0     = (1.-16.*zq_on_l)**(-0.5)
    pm1     = (1.-16.*z_on_l)**(-0.25)
    ph1     = (1.-16.*zs_on_l)**(-0.5) !=pq1
    integralm = lzom-2.*log((1.+1./pm1)/(1.+1./pm0)) &
                -log((1.+1./pm1**2)/(1.+1./pm0**2)) &
                +2.*(atan(1./pm1)-atan(1./pm0))
    integralh = lzoh-2.*log((1.+1./ph1)/(1.+1./ph0))
    integralq = lzoq-2.*log((1.+1./ph1)/(1.+1./pq0))
  elsewhere
    !--------------Beljaars and Holtslag (1991) momentum & heat            
    pm0 = -(a_1*z0_on_l+b_1*(z0_on_l-(c_1/d_1))*exp(-d_1*z0_on_l) &
           +b_1*c_1/d_1)
    pm1 = -(a_1*z_on_l+b_1*(z_on_l-(c_1/d_1))*exp(-d_1*z_on_l) &
           +b_1*c_1/d_1)
    ph0 = -((1.+(2./3.)*a_1*zt_on_l)**1.5 &
           +b_1*(zt_on_l-(c_1/d_1))*exp(-d_1*zt_on_l) &
           +b_1*c_1/d_1-1.)
    ph1 = -((1.+(2./3.)*a_1*zs_on_l)**1.5 &
           +b_1*(zs_on_l-(c_1/d_1))*exp(-d_1*zs_on_l) &
           +b_1*c_1/d_1-1.) !=pq1
    pq0 = -((1.+(2./3.)*a_1*zq_on_l)**1.5 &
           +b_1*(zq_on_l-(c_1/d_1))*exp(-d_1*zq_on_l) &
           +b_1*c_1/d_1-1.)
    integralm = lzom-(pm1-pm0)
    integralh = lzoh-(ph1-ph0)
    integralq = lzoq-(ph1-pq0)
  endwhere
  thetavstar=vkar*(thetav-sthetav)/integralh
  ustar=vkar*umag/integralm
end do
tstar=vkar*(a_temp-stemp)/integralh
qstar=vkar*(a_qg-smixr)/integralq
      
! estimate screen diagnostics
z0s_on_l=z0*zs_on_l/a_zmins
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
  ph1     = (1.-16.*zs_on_l)**(-0.50)
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
  ph1  = -((1.+(2./3.)*a_1*zs_on_l)**1.5 &
         +b_1*(zs_on_l-(c_1/d_1)) &
         *exp(-d_1*zs_on_l)+b_1*c_1/d_1-1.)
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
