
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
       mloscrnout,mloextra,mloimpice,mloexpice,mloexpdep,mloexpdensity,mloexpmelt,mloexpscalar,wlev, &
       micdwn,mxd,mindep,minwater,onedice,mloimport3d,mloexport3d

! parameters
integer, save :: wlev = 20
integer, parameter :: iqx = 4148
! model arrays
integer, save :: wfull,ifull,iqwt
integer, dimension(:), allocatable, save :: wmap
integer, dimension(:), allocatable, save :: p_mixind
logical, dimension(:), allocatable, save :: wpack
real, dimension(:,:), allocatable, save :: depth,dz,depth_hl,dz_hl
real, dimension(:,:), allocatable, save :: micdwn          ! This variable is for CCAM onthefly.f
real, dimension(:,:), allocatable, save :: w_temp,w_sal,w_u,w_v
real, dimension(:,:), allocatable, save :: i_tn
real, dimension(:), allocatable, save :: w_eta
real, dimension(:), allocatable, save :: i_dic,i_dsn,i_fracice,i_tsurf,i_sto,i_u,i_v,i_sal
real, dimension(:), allocatable, save :: p_mixdepth,p_bf
real, dimension(:), allocatable, save :: p_watervisdiralb,p_watervisdifalb,p_waternirdiralb,p_waternirdifalb
real, dimension(:), allocatable, save :: p_icevisdiralb,p_icevisdifalb,p_icenirdiralb,p_icenirdifalb
real, dimension(:), allocatable, save :: p_zo,p_zoh,p_zoq,p_cd,p_cds,p_cdq,p_fg,p_eg
real, dimension(:), allocatable, save :: p_zoice,p_zohice,p_zoqice,p_cdice,p_cdsice,p_cdqice,p_fgice,p_egice,p_wetfacice
real, dimension(:), allocatable, save :: p_tscrn,p_uscrn,p_qgscrn,p_u10,p_taux,p_tauy,p_tauxica,p_tauyica
  
! mode
integer, parameter :: incradbf  = 1 ! include shortwave in buoyancy forcing
integer, parameter :: incradgam = 1 ! include shortwave in non-local term
integer, parameter :: zomode    = 2 ! roughness calculation (0=Charnock (CSIRO9), 1=Charnock (zot=zom), 2=Beljaars)
integer, parameter :: mixmeth   = 1 ! Refine mixed layer depth calculation (0=None, 1=Iterative)
integer, parameter :: deprelax  = 0 ! surface height (0=vary, 1=relax, 2=set to zero)
integer, save      :: onedice   = 1 ! use 1D ice model (0=Off, 1=On)
! model parameters
real, save :: mxd      = 5002.18 ! Max depth (m)
real, save :: mindep   = 1.      ! Thickness of first layer (m)
real, save :: minwater = 5.      ! Minimum water depth (m)
real, parameter :: ric     = 0.3    ! Critical Ri for diagnosing mixed layer depth
real, parameter :: epsilon = 0.1    ! Ratio of surface layer and mixed layer thickness
real, parameter :: minsfc  = 0.5    ! Minimum thickness to average surface layer properties (m)
real, parameter :: maxsal  = 50.    ! Maximum salinity used in density and melting point calculations (PSU)
real, parameter :: mu_1    = 23.    ! VIS depth (m) - Type I
real, parameter :: mu_2    = 0.35   ! NIR depth (m) - Type I
real, parameter :: fluxwgt = 0.7    ! Time filter for flux calculations
! physical parameters
real, parameter :: vkar=0.4               ! von Karman constant
real, parameter :: lv=2.501e6             ! Latent heat of vaporisation (J kg^-1)
real, parameter :: lf=3.337e5             ! Latent heat of fusion (J kg^-1)
real, parameter :: ls=lv+lf               ! Latent heat of sublimation (J kg^-1)
real, parameter :: grav=9.80              ! graviational constant (m s^-2)
real, parameter :: sbconst=5.67e-8        ! Stefan-Boltzmann constant
real, parameter :: cdbot=2.4e-3           ! bottom drag coefficent
real, parameter :: cp0=3990.              ! heat capacity of mixed layer (J kg^-1 K^-1)
real, parameter :: cp=1004.64             ! Specific heat of dry air at const P
real, parameter :: rdry=287.04            ! Specific gas const for dry air
real, parameter :: rvap=461.5             ! Gas constant for water vapor
! ice parameters
real, parameter :: himin=0.1              ! minimum ice thickness for multiple layers (m)
real, parameter :: icebreak=0.05          ! minimum ice thickness before breakup (1D model)
real, parameter :: fracbreak=0.05         ! minimum ice fraction (1D model)
real, parameter :: icemin=0.01            ! minimum ice thickness (m)
real, parameter :: icemax=10.             ! maximum ice thickness (m)
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
real, parameter :: icesal=10.             ! Maximum salinity for sea-ice
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

if (diag>=1) write(6,*) "Initialising MLO"

ifull=ifin
allocate(wpack(ifull))
wpack=depin>minwater
wfull=count(wpack)
if (wfull==0) then
  deallocate(wpack)
  return
end if

allocate(wmap(wfull))
allocate(w_temp(wfull,wlev),w_sal(wfull,wlev))
allocate(w_u(wfull,wlev),w_v(wfull,wlev))
allocate(w_eta(wfull))
allocate(i_tn(wfull,0:2),i_dic(wfull),i_dsn(wfull))
allocate(i_fracice(wfull),i_tsurf(wfull),i_sto(wfull))
allocate(i_u(wfull),i_v(wfull),i_sal(wfull))
allocate(p_mixdepth(wfull),p_bf(wfull))
allocate(p_watervisdiralb(wfull),p_watervisdifalb(wfull))
allocate(p_waternirdiralb(wfull),p_waternirdifalb(wfull))
allocate(p_icevisdiralb(wfull),p_icevisdifalb(wfull))
allocate(p_icenirdiralb(wfull),p_icenirdifalb(wfull))
allocate(p_zo(wfull),p_zoh(wfull),p_zoq(wfull),p_cd(wfull),p_cds(wfull),p_fg(wfull),p_eg(wfull))
allocate(p_zoice(wfull),p_zohice(wfull),p_zoqice(wfull),p_cdice(wfull),p_cdsice(wfull),p_fgice(wfull),p_egice(wfull))
allocate(p_cdq(wfull),p_cdqice(wfull))
allocate(p_wetfacice(wfull),p_mixind(wfull))
allocate(p_tscrn(wfull),p_uscrn(wfull),p_qgscrn(wfull),p_u10(wfull))
allocate(p_taux(wfull),p_tauy(wfull))
allocate(p_tauxica(wfull),p_tauyica(wfull))
allocate(depth(wfull,wlev),dz(wfull,wlev))
allocate(depth_hl(wfull,wlev+1))
allocate(dz_hl(wfull,2:wlev))

iqw=0
do iq=1,ifull
  if (wpack(iq)) then
    iqw=iqw+1
    wmap(iqw)=iq
  end if
end do

iqw=0
iqwt=0
do iq=1,ifull
  if (wpack(iq)) then
    iqw=iqw+1
    if (iq>=iqx) then
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
i_sal=0.                ! PSU

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
p_cdq=0.
p_fg=0.
p_eg=0.
p_zoice=0.001
p_zohice=0.001
p_zoqice=0.001
p_cdice=0.
p_cdsice=0.
p_cdqice=0.
p_fgice=0.
p_egice=0.
p_wetfacice=1.
p_tscrn=273.2
p_qgscrn=0.
p_uscrn=0.
p_u10=0.
p_mixind=wlev-1
p_taux=0.
p_tauy=0.
p_tauxica=0.
p_tauyica=0.

! MLO - 20 level
!depth = (/   0.5,   3.4,  11.9,  29.7,  60.7, 108.5, 176.9, 269.7, 390.6, 543.3, &
!           731.6, 959.2,1230.0,1547.5,1915.6,2338.0,2818.5,3360.8,3968.7,4645.8 /)
! MLO - 30 level
!depth = (/   0.5,   2.1,   5.3,  11.3,  21.1,  36.0,  56.9,  85.0, 121.5, 167.3, &
!           223.7, 291.7, 372.4, 467.0, 576.5, 702.1, 844.8,1005.8,1186.2,1387.1, &
!          1609.6,1854.7,2123.7,2417.6,2737.5,3084.6,3456.9,3864.5,4299.6,4766.2 /)
! MLO - 40 level
!depth = (/   0.5,   1.7,   3.7,   6.8,  11.5,  18.3,  27.7,  40.1,  56.0,  75.9, &
!           100.2, 129.4, 164.0, 204.3, 251.0, 304.4, 365.1, 433.4, 509.9, 595.0, &
!           689.2, 793.0, 906.7,1031.0,1166.2,1312.8,1471.3,1642.2,1825.9,2022.8, &
!          2233.5,2458.4,2698.0,2952.8,3233.1,3509.5,3812.5,4132.4,4469.9,4825.3 /)
! MLO - 50 level
!depth = (/   0.5,   1.7,   3.1,   5.2,   8.1,  12.1,  17.3,  24.2,  32.8,  43.4, &
!            56.3,  71.7,  89.9, 111.0, 135.3, 163.1, 194.5, 229.9, 269.5, 313.4, &
!           362.0, 415.5, 474.1, 538.1, 607.6, 683.0, 764.4, 852.2, 946.5,1047.6, &
!          1155.7,1271.0,1393.9,1524.5,1663.0,1809.8,1965.0,2128.9,2301.8,2483.8, &
!          2675.2,2876.2,3087.1,3308.1,3539.5,3781.5,4034.3,4298.2,4573.4,4860.1 /)
! Mk3.5 - 31 level
!depth = (/   5.0,  15.0,  28.2,  42.0,  59.7,  78.5, 102.1, 127.9, 159.5, 194.6, &
!           237.0, 284.7, 341.7, 406.4, 483.2, 570.9, 674.9, 793.8, 934.1,1095.2, &
!          1284.7,1502.9,1758.8,2054.3,2400.0,2800.0,3200.0,3600.0,4000.0,4400.0, &
!          4800.0 /)

deptmp(1:wfull)=depin(wmap)

!smxd=maxval(deptmp(1:wfull))
!smnd=minval(deptmp(1:wfull))

do iqw=1,wfull
  call vgrid(wlev,deptmp(iqw),dumdf,dumdh)
  depth(iqw,:)=dumdf
  depth_hl(iqw,:)=dumdh
  !if (smxd==deptmp(iqw)) then
  !  write(6,*) "MLO max depth ",depth(iqw,:)
  !  smxd=smxd+10.
  !end if
  !if (smnd==deptmp(iqw)) then
  !  write (6,*) "MLO min depth ",depth(iqw,:)
  !  smnd=smnd-10.
  !end if
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

subroutine vgrid(wlin,depin,depthout,depth_hlout)

implicit none

integer, intent(in) :: wlin
integer ii
real, intent(in) :: depin
real, dimension(wlin), intent(out) :: depthout
real, dimension(wlin+1), intent(out) :: depth_hlout
real dd,x,al,bt

dd=min(mxd,max(mindep,depin))
x=real(wlin)
!if (dd>x) then
!  al=(mindep*x-dd)/(x-x*x*x)           ! hybrid levels
!  bt=(mindep*x*x*x-dd)/(x*x*x-x)       ! hybrid levels
  al=dd*(mindep*x/mxd-1.)/(x-x*x*x)     ! sigma levels
  bt=dd*(mindep*x*x*x/mxd-1.)/(x*x*x-x) ! sigma levels 
  do ii=1,wlin+1
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
do ii=1,wlin
  depthout(ii)=0.5*(depth_hlout(ii)+depth_hlout(ii+1))
end do

return
end subroutine vgrid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Deallocate MLO arrays

subroutine mloend

implicit none

if (wfull==0) return

deallocate(wmap,wpack)
deallocate(w_temp,w_sal,w_u,w_v)
deallocate(w_eta)
deallocate(i_tn,i_dic,i_dsn)
deallocate(i_fracice,i_tsurf,i_sto)
deallocate(i_u,i_v,i_sal)
deallocate(p_mixdepth,p_bf)
deallocate(p_watervisdiralb,p_watervisdifalb)
deallocate(p_waternirdiralb,p_waternirdifalb)
deallocate(p_icevisdiralb,p_icevisdifalb)
deallocate(p_icenirdiralb,p_icenirdifalb)
deallocate(p_zo,p_zoh,p_zoq,p_cd,p_cds,p_fg,p_eg)
deallocate(p_zoice,p_zohice,p_zoqice,p_cdice,p_cdsice,p_fgice,p_egice)
deallocate(p_cdq,p_cdqice)
deallocate(p_wetfacice,p_mixind)
deallocate(p_tscrn,p_uscrn,p_qgscrn,p_u10)
deallocate(p_taux,p_tauy,p_tauxica,p_tauyica)
deallocate(depth,dz,depth_hl,dz_hl)

return
end subroutine mloend

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Load MLO data

subroutine mloload(datain,shin,icein,diag)

implicit none

integer, intent(in) :: diag
real, dimension(ifull,wlev,4), intent(in) :: datain
real, dimension(ifull,11), intent(in) :: icein
real, dimension(ifull), intent(in) :: shin

if (wfull==0) then
  return
else if (wfull==ifull) then
  w_temp(:,:)=datain(:,:,1)
  w_sal(:,:) =datain(:,:,2)
  w_u(:,:)   =datain(:,:,3)
  w_v(:,:)   =datain(:,:,4)
  i_tsurf(:)  =icein(:,1)
  i_tn(:,0)   =icein(:,2)
  i_tn(:,1)   =icein(:,3)
  i_tn(:,2)   =icein(:,4)
  i_fracice(:)=icein(:,5)
  i_dic(:)    =icein(:,6)
  i_dsn(:)    =icein(:,7)
  i_sto(:)    =icein(:,8)
  i_u(:)      =icein(:,9)
  i_v(:)      =icein(:,10)
  i_sal(:)    =icein(:,11)
  w_eta(:)    =shin(:)
else
  w_temp(:,:)=datain(wmap,:,1)
  w_sal(:,:) =datain(wmap,:,2)
  w_u(:,:)   =datain(wmap,:,3)
  w_v(:,:)   =datain(wmap,:,4)
  i_tsurf(:)  =icein(wmap,1)
  i_tn(:,0)   =icein(wmap,2)
  i_tn(:,1)   =icein(wmap,3)
  i_tn(:,2)   =icein(wmap,4)
  i_fracice(:)=icein(wmap,5)
  i_dic(:)    =icein(wmap,6)
  i_dsn(:)    =icein(wmap,7)
  i_sto(:)    =icein(wmap,8)
  i_u(:)      =icein(wmap,9)
  i_v(:)      =icein(wmap,10)
  i_sal(:)    =icein(wmap,11)
  w_eta(:)    =shin(wmap)
end if

return
end subroutine mloload

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Save MLO data

subroutine mlosave(dataout,depout,shout,iceout,diag)

implicit none

integer, intent(in) :: diag
real, dimension(ifull,wlev,4), intent(inout) :: dataout
real, dimension(ifull,11), intent(inout) :: iceout
real, dimension(ifull), intent(inout) :: depout,shout

if (wfull==0) then
  return
else if (wfull==ifull) then
  dataout(:,:,1)=w_temp(:,:)
  dataout(:,:,2)=w_sal(:,:)
  dataout(:,:,3)=w_u(:,:)
  dataout(:,:,4)=w_v(:,:)
  iceout(:,1)  =i_tsurf(:)
  iceout(:,2)  =i_tn(:,0)
  iceout(:,3)  =i_tn(:,1)
  iceout(:,4)  =i_tn(:,2)
  iceout(:,5)  =i_fracice(:)
  iceout(:,6)  =i_dic(:)
  iceout(:,7)  =i_dsn(:)
  iceout(:,8)  =i_sto(:)
  iceout(:,9)  =i_u(:)
  iceout(:,10) =i_v(:)
  iceout(:,11) =i_sal(:)
  depout(:)=depth_hl(:,wlev+1)
  shout(:)=w_eta(:)
else
  iceout(:,8)=0.
  depout=0.
  shout=0.

  dataout(wmap,:,1)=w_temp(:,:)
  dataout(wmap,:,2)=w_sal(:,:)
  dataout(wmap,:,3)=w_u(:,:)
  dataout(wmap,:,4)=w_v(:,:)
  iceout(wmap,1)  =i_tsurf(:)
  iceout(wmap,2)  =i_tn(:,0)
  iceout(wmap,3)  =i_tn(:,1)
  iceout(wmap,4)  =i_tn(:,2)
  iceout(wmap,5)  =i_fracice(:)
  iceout(wmap,6)  =i_dic(:)
  iceout(wmap,7)  =i_dsn(:)
  iceout(wmap,8)  =i_sto(:)
  iceout(wmap,9)  =i_u(:)
  iceout(wmap,10) =i_v(:)
  iceout(wmap,11) =i_sal(:)
  depout(wmap)=depth_hl(:,wlev+1)
  shout(wmap)=w_eta(:)
end if

return
end subroutine mlosave

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Import sst for nudging

subroutine mloimport(mode,sst,ilev,diag)

implicit none

integer, intent(in) :: mode,ilev,diag
real, dimension(ifull), intent(in) :: sst

if (wfull==0) then
  return
else if (wfull==ifull) then
  select case(mode)
    case(0)
      w_temp(:,ilev)=sst(:)
    case(1)
      w_sal(:,ilev)=sst(:)
    case(2)
      w_u(:,ilev)=sst(:)
    case(3)
      w_v(:,ilev)=sst(:)
    case(4)
      w_eta=sst(:)
  end select
else
  select case(mode)
    case(0)
      w_temp(:,ilev)=sst(wmap)
    case(1)
      w_sal(:,ilev)=sst(wmap)
    case(2)
      w_u(:,ilev)=sst(wmap)
    case(3)
      w_v(:,ilev)=sst(wmap)
    case(4)
      w_eta=sst(wmap)
  end select
end if

return
end subroutine mloimport

subroutine mloimport3d(mode,sst,diag)

implicit none

integer, intent(in) :: mode,diag
real, dimension(ifull,wlev), intent(in) :: sst

if (wfull==0) then
  return
else if (wfull==ifull) then
  select case(mode)
    case(0)
      w_temp(:,:)=sst(:,:)
    case(1)
      w_sal(:,:) =sst(:,:)
    case(2)
      w_u(:,:)   =sst(:,:)
    case(3)
      w_v(:,:)   =sst(:,:)
  end select
else
  select case(mode)
    case(0)
      w_temp(:,:)=sst(wmap,:)
    case(1)
      w_sal(:,:) =sst(wmap,:)
    case(2)
      w_u(:,:)   =sst(wmap,:)
    case(3)
      w_v(:,:)   =sst(wmap,:)
  end select
end if

return
end subroutine mloimport3d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Export ice temp

subroutine mloimpice(tsn,ilev,diag)

implicit none

integer, intent(in) :: ilev,diag
real, dimension(ifull), intent(inout) :: tsn

if (wfull==0) then
  return
else if (wfull==ifull) then
  select case(ilev)
    case(1)
      i_tsurf=tsn(:)
    case(2)
      i_tn(:,0)=tsn(:)
    case(3)
      i_tn(:,1)=tsn(:)
    case(4)
      i_tn(:,2)=tsn(:)
    case(5)
      i_fracice=tsn(:)
    case(6)
      i_dic=tsn(:)
    case(7)
      i_dsn=tsn(:)
    case(8)
      i_sto=tsn(:)
    case(9)
      i_u=tsn(:)
    case(10)
      i_v=tsn(:)
    case(11)
      i_sal=tsn(:)
    case DEFAULT
      write(6,*) "ERROR: Invalid mode ",ilev
      stop
  end select
else
  select case(ilev)
    case(1)
      i_tsurf=tsn(wmap)
    case(2)
      i_tn(:,0)=tsn(wmap)
    case(3)
      i_tn(:,1)=tsn(wmap)
    case(4)
      i_tn(:,2)=tsn(wmap)
    case(5)
      i_fracice=tsn(wmap)
    case(6)
      i_dic=tsn(wmap)
    case(7)
      i_dsn=tsn(wmap)
    case(8)
      i_sto=tsn(wmap)
    case(9)
      i_u=tsn(wmap)
    case(10)
      i_v=tsn(wmap)
    case(11)
      i_sal=tsn(wmap)
    case DEFAULT
      write(6,*) "ERROR: Invalid mode ",ilev
      stop
  end select
end if

return
end subroutine mloimpice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Export sst for nudging

subroutine mloexport(mode,sst,ilev,diag)

implicit none

integer, intent(in) :: mode,ilev,diag
real, dimension(ifull), intent(inout) :: sst

if (wfull==0) then
  return
else if (wfull==ifull) then
  select case(mode)
    case(0)
      sst(:)=w_temp(:,ilev)
    case(1)
      sst(:)=w_sal(:,ilev)
    case(2)
      sst(:)=w_u(:,ilev)
    case(3)
      sst(:)=w_v(:,ilev)
    case(4)
      sst(:)=w_eta
  end select
else
  select case(mode)
    case(0)
      sst(wmap)=w_temp(:,ilev)
    case(1)
      sst(wmap)=w_sal(:,ilev)
    case(2)
      sst(wmap)=w_u(:,ilev)
    case(3)
      sst(wmap)=w_v(:,ilev)
    case(4)
      sst(wmap)=w_eta
  end select
end if

return
end subroutine mloexport

subroutine mloexport3d(mode,sst,diag)

implicit none

integer, intent(in) :: mode,diag
real, dimension(ifull,wlev), intent(inout) :: sst

if (wfull==0) then
  return
else if (wfull==ifull) then
  select case(mode)
    case(0)
      sst(:,:)=w_temp(:,:)
    case(1)
      sst(:,:)=w_sal(:,:)
    case(2)
      sst(:,:)=w_u(:,:)
    case(3)
      sst(:,:)=w_v(:,:)
  end select
else
  select case(mode)
    case(0)
      sst(wmap,:)=w_temp(:,:)
    case(1)
      sst(wmap,:)=w_sal(:,:)
    case(2)
      sst(wmap,:)=w_u(:,:)
    case(3)
      sst(wmap,:)=w_v(:,:)
  end select
end if

return
end subroutine mloexport3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Export ice temp

subroutine mloexpice(tsn,ilev,diag)

implicit none

integer, intent(in) :: ilev,diag
real, dimension(ifull), intent(inout) :: tsn

if (wfull==0) then
  return
else if (wfull==ifull) then
  select case(ilev)
    case(1)
      tsn(:)=i_tsurf
    case(2)
      tsn(:)=i_tn(:,0)
    case(3)
      tsn(:)=i_tn(:,1)
    case(4)
      tsn(:)=i_tn(:,2)
    case(5)
      tsn(:)=i_fracice
    case(6)
      tsn(:)=i_dic
    case(7)
      tsn(:)=i_dsn
    case(8)
      tsn(:)=i_sto
    case(9)
      tsn(:)=i_u
    case(10)
      tsn(:)=i_v
    case(11)
      tsn(:)=i_sal
    case DEFAULT
      write(6,*) "ERROR: Invalid mode ",ilev
      stop
  end select
else
  select case(ilev)
    case(1)
      tsn(wmap)=i_tsurf
    case(2)
      tsn(wmap)=i_tn(:,0)
    case(3)
      tsn(wmap)=i_tn(:,1)
    case(4)
      tsn(wmap)=i_tn(:,2)
    case(5)
      tsn(wmap)=i_fracice
    case(6)
      tsn(wmap)=i_dic
    case(7)
      tsn(wmap)=i_dsn
    case(8)
      tsn(wmap)=i_sto
    case(9)
      tsn(wmap)=i_u
    case(10)
      tsn(wmap)=i_v
    case(11)
      tsn(wmap)=i_sal
    case DEFAULT
      write(6,*) "ERROR: Invalid mode ",ilev
      stop
  end select
end if

return
end subroutine mloexpice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Return mixed layer depth

subroutine mlodiag(mld,diag)

implicit none

integer, intent(in) :: diag
real, dimension(ifull), intent(out) :: mld

mld=0.
if (wfull==0) return
mld(wmap)=p_mixdepth

return
end subroutine mlodiag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Return roughness length for heat

subroutine mloextra(mode,zoh,zmin,diag)

implicit none

integer, intent(in) :: mode,diag
real, dimension(ifull), intent(out) :: zoh
real, dimension(ifull), intent(in) :: zmin
real, dimension(wfull) :: a_zmin
real, dimension(wfull) :: workb,workc

zoh=0.
if (wfull==0) return

select case(mode)
  case(0) ! zoh
    a_zmin=zmin(wmap)
    workb=(1.-i_fracice)/(log(a_zmin/p_zo)*log(a_zmin/p_zoh))+i_fracice/(log(a_zmin/p_zoice)*log(a_zmin/p_zohice))
    workc=(1.-i_fracice)/log(a_zmin/p_zo)**2+i_fracice/log(a_zmin/p_zoice)**2
    workc=sqrt(workc)
    zoh(wmap)=a_zmin*exp(-workc/workb)
  case(1) ! taux
    workb=(1.-i_fracice)*p_taux+i_fracice*p_tauxica
    zoh(wmap)=workb
  case(2) ! tauy
    workb=(1.-i_fracice)*p_tauy+i_fracice*p_tauyica
    zoh(wmap)=workb
  case(3) ! zoq
    a_zmin=zmin(wmap)
    workb=(1.-i_fracice)/(log(a_zmin/p_zo)*log(a_zmin/p_zoq))+i_fracice/(log(a_zmin/p_zoice)*log(a_zmin/p_zoqice))
    workc=(1.-i_fracice)/log(a_zmin/p_zo)**2+i_fracice/log(a_zmin/p_zoice)**2
    workc=sqrt(workc)
    zoh(wmap)=a_zmin*exp(-workc/workb)
  case default
    write(6,*) "ERROR: Invalid mode ",mode
    stop
end select

return
end subroutine mloextra

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate screen diagnostics

subroutine mloscrnout(tscrn,qgscrn,uscrn,u10,diag)

implicit none

integer, intent(in) :: diag
real, dimension(ifull), intent(inout) :: tscrn,qgscrn,uscrn,u10

if (wfull==0) return

tscrn(wmap) =p_tscrn
qgscrn(wmap)=p_qgscrn
uscrn(wmap) =p_uscrn
u10(wmap)   =p_u10

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

if (wfull==0) return

ifinish=istart+ifin-1
ib=count(wpack(1:istart-1))+1
ie=count(wpack(istart:ifinish))+ib-1
if (ie<ib) return

costmp(ib:ie)=coszro(wmap(ib:ie)-istart+1)

!pond(ib:ie)=max(1.+.008*min(i_tsurf(ib:ie)-273.16,0.),0.)
pond(ib:ie)=0.
!snow(ib:ie)=min(max(i_dsn(ib:ie)/0.05,0.),1.)
snow(ib:ie)=0.

watervis(ib:ie)=.05/(costmp(ib:ie)+0.15)
waternir(ib:ie)=.05/(costmp(ib:ie)+0.15)
! need to factor snow age into albedo
icevis(ib:ie)=(0.85*(1.-pond(ib:ie))+0.75*pond(ib:ie))*(1.-snow(ib:ie))+(0.95*(1.-pond(ib:ie))+0.85*pond(ib:ie))*snow(ib:ie)
icenir(ib:ie)=(0.45*(1.-pond(ib:ie))+0.35*pond(ib:ie))*(1.-snow(ib:ie))+(0.65*(1.-pond(ib:ie))+0.55*pond(ib:ie))*snow(ib:ie)
ovisalb(wmap(ib:ie)-istart+1)=icevis(ib:ie)*i_fracice(ib:ie)+(1.-i_fracice(ib:ie))*watervis(ib:ie)
oniralb(wmap(ib:ie)-istart+1)=icenir(ib:ie)*i_fracice(ib:ie)+(1.-i_fracice(ib:ie))*waternir(ib:ie)

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

if (wfull==0) return

ifinish=istart+ifin-1
ib=count(wpack(1:istart-1))+1
ie=count(wpack(istart:ifinish))+ib-1
if (ie<ib) return

costmp(ib:ie)=coszro(wmap(ib:ie)-istart+1)

!pond(ib:ie)=max(1.+.008*min(i_tsurf(ib:ie)-273.16,0.),0.)
pond(ib:ie)=0.
!snow(ib:ie)=min(max(i_dsn(ib:ie)/0.05,0.),1.)
snow(ib:ie)=0.

where (costmp(ib:ie)>0.)
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
ovisdir(wmap(ib:ie)-istart+1)=i_fracice(ib:ie)*p_icevisdiralb(ib:ie)+(1.-i_fracice(ib:ie))*p_watervisdiralb(ib:ie)
ovisdif(wmap(ib:ie)-istart+1)=i_fracice(ib:ie)*p_icevisdifalb(ib:ie)+(1.-i_fracice(ib:ie))*p_watervisdifalb(ib:ie)
onirdir(wmap(ib:ie)-istart+1)=i_fracice(ib:ie)*p_icenirdiralb(ib:ie)+(1.-i_fracice(ib:ie))*p_waternirdiralb(ib:ie)
onirdif(wmap(ib:ie)-istart+1)=i_fracice(ib:ie)*p_icenirdifalb(ib:ie)+(1.-i_fracice(ib:ie))*p_waternirdifalb(ib:ie)

return
end subroutine mloalb4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Regrid MLO data

subroutine mloregrid(wlin,depin,mloin,mlodat,mode)

implicit none

integer, intent(in) :: wlin,mode
integer iqw,ii,pos(1),nn
real, dimension(ifull), intent(in) :: depin
real, dimension(ifull,wlin), intent(in) :: mloin
real, dimension(ifull,wlev), intent(inout) :: mlodat
real, dimension(wfull,wlin) :: newdata
real, dimension(wfull,wlev) :: newdatb
real, dimension(wfull) :: deptmp
real, dimension(wlin) :: dpin,sgin
real, dimension(wlin+1) :: dp_hlin
real, dimension(wlev) :: sig
real x

if (wfull==0) return

deptmp=depin(wmap)
do ii=1,wlin
  newdata(:,ii)=mloin(wmap,ii)
end do

select case(mode)
  case(0,1)
    do iqw=1,wfull
      call vgrid(wlin,deptmp(iqw),dpin,dp_hlin)
      do ii=1,wlev
        if (depth(iqw,ii)>=dpin(wlin)) then
          newdatb(iqw,ii)=newdata(iqw,wlin)
        else if (depth(iqw,ii)<=dpin(1)) then
          newdatb(iqw,ii)=newdata(iqw,1)
        else
          pos=maxloc(dpin,dpin<depth(iqw,ii))
          pos(1)=max(1,min(wlin-1,pos(1)))
          x=(depth(iqw,ii)-dpin(pos(1)))/(dpin(pos(1)+1)-dpin(pos(1)))
          x=max(0.,min(1.,x))
          newdatb(iqw,ii)=newdata(iqw,pos(1)+1)*x+newdata(iqw,pos(1))*(1.-x)
        end if
      end do
    end do
  case(2,3)
    do iqw=1,wfull
      sig=depth(iqw,:)/depth_hl(iqw,wlev)
      call vgrid(wlin,deptmp(iqw),dpin,dp_hlin)
      sgin=dpin/dp_hlin(wlin)
      do ii=1,wlev
        if (sig(ii)>=sgin(wlin)) then
          newdatb(iqw,ii)=newdata(iqw,wlin)
        else if (sig(ii)<=sgin(1)) then
          newdatb(iqw,ii)=newdata(iqw,1)
        else
          pos=maxloc(sgin,sgin<sig(ii))
          pos(1)=max(1,min(wlin-1,pos(1)))
          x=(sig(ii)-sgin(pos(1)))/(sgin(pos(1)+1)-sgin(pos(1)))
          x=max(0.,min(1.,x))
          newdatb(iqw,ii)=newdata(iqw,pos(1)+1)*x+newdata(iqw,pos(1))*(1.-x)
        end if
      end do
    end do
  case default
    write(6,*) "ERROR: Unknown mode for mloregrid ",mode
    stop
end select

do ii=1,wlev
  mlodat(wmap,ii)=newdatb(:,ii)
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
if (wfull==0) return

select case(mode)
  case(0)
    odep(wmap)=depth(:,ii)
  case(1)
    odep(wmap)=dz(:,ii)
  case default
    write(6,*) "ERROR: Unknown mloexpdep mode ",mode
    stop
end select
  
return
end subroutine mloexpdep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Extract density

subroutine mloexpdensity(odensity,alpha,beta,tt,ss,ddz,pxtr,diag)

implicit none

integer, intent(in) :: diag
real, dimension(:,:), intent(in) :: tt
real, dimension(size(tt,1)), intent(in) :: pxtr
real, dimension(size(tt,1),size(tt,2)), intent(in) :: ss,ddz
real, dimension(size(tt,1),size(tt,2)), intent(out) :: odensity,alpha,beta
real, dimension(size(tt,1),size(tt,2)) :: rs
real, dimension(size(tt,1)) :: rho0

call calcdensity(odensity,alpha,beta,rs,rho0,tt,ss,ddz,pxtr)

return
end subroutine mloexpdensity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Extract melting temperature

subroutine mloexpmelt(omelt)

implicit none

real, dimension(ifull), intent(out) :: omelt
real, dimension(wfull) :: tmelt,d_zcr

omelt=273.16
if (wfull==0) return
d_zcr=max(1.+w_eta/depth_hl(:,wlev+1),minwater/depth_hl(:,wlev+1))
call calcmelt(tmelt,d_zcr)
omelt(wmap)=tmelt

return
end subroutine mloexpmelt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Extract scalars

subroutine mloexpscalar(mode,ogamm,diag)

implicit none

integer, intent(in) :: mode,diag
integer, dimension(wfull) :: nk
real, dimension(ifull), intent(out) :: ogamm
real, dimension(wfull) :: gamm

ogamm=1.
if (wfull==0) return

select case(mode)
  case(0)
    nk=min(int(i_dic/himin),2)
    where (nk>0.and.i_dsn>icemin) ! snow + 1-2 ice layers
      gamm=gamms
    elsewhere (nk>0.and.i_dsn<=icemin) ! no snow + 1-2 ice layers
      gamm=gammi
    elsewhere (nk<=0.and.i_dsn>icemin) ! snow + 1 ice layer
      gamm=(gammi*max(i_dic,icemin)+gamms*i_dsn)/(max(i_dic,icemin)+i_dsn)
    elsewhere  !(nk<=0.and.i_dsn<=icemin) ! no snow + 1 ice layer
      gamm=gammi
    end where
  case default
    write(6,*) "ERROR: Unknown mloexpscalar mode ",mode
    stop
end select

ogamm(wmap)=gamm

return
end subroutine mloexpscalar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Pack atmospheric data for MLO eval

subroutine mloeval(sst,zo,cd,cds,fg,eg,wetfac,epot,epan,fracice,siced,snowd, &
                   dt,zmin,zmins,sg,rg,precp,precs,uatm,vatm,temp,qg,ps,f,   &
                   visnirratio,fbvis,fbnir,inflow,diag,calcprog,oldu,oldv)

implicit none

integer, intent(in) :: diag
integer iq,iqw,ii
real, intent(in) :: dt
real, dimension(ifull), intent(in) :: sg,rg,precp,precs,f,uatm,vatm,temp,qg,ps,visnirratio,fbvis,fbnir,inflow,zmin,zmins
real, dimension(ifull), intent(inout) :: sst,zo,cd,cds,fg,eg,wetfac,fracice,siced,epot,epan,snowd
real, dimension(ifull), intent(in), optional :: oldu,oldv
real, dimension(wfull) :: workb,workc
real, dimension(wfull) :: a_sg,a_rg,a_rnd,a_snd,a_f,a_vnratio,a_fbvis,a_fbnir,a_u,a_v,a_temp,a_qg,a_ps,a_zmin,a_zmins
real, dimension(wfull) :: a_inflow
real, dimension(wfull,wlev) :: d_rho,d_rs,d_nsq,d_rad,d_alpha,d_beta
real, dimension(wfull) :: d_b0,d_ustar,d_wu0,d_wv0,d_wt0,d_ws0,d_ftop,d_tb,d_zcr
real, dimension(wfull) :: d_fb,d_timelt,d_tauxicw,d_tauyicw,d_neta,d_ndsn
real, dimension(wfull) :: d_ndic,d_nsto,d_oldu,d_oldv
integer, dimension(wfull) :: d_nk
logical, intent(in) :: calcprog

if (wfull==0) return

a_sg     =sg(wmap)
a_rg     =rg(wmap)
a_f      =f(wmap)
a_vnratio=visnirratio(wmap)
a_fbvis  =fbvis(wmap)
a_fbnir  =fbnir(wmap)
a_u      =uatm(wmap)
a_v      =vatm(wmap)
a_temp   =temp(wmap)
a_qg     =qg(wmap)
a_ps     =ps(wmap)
a_zmin   =zmin(wmap)
a_zmins  =zmins(wmap)
a_inflow =inflow(wmap)
a_rnd    =precp(wmap)-precs(wmap)
a_snd    =precs(wmap)
if (present(oldu).and.present(oldv)) then
  d_oldu=oldu(wmap)
  d_oldv=oldv(wmap)
else
  d_oldu=w_u(:,1)
  d_oldv=w_v(:,1)
end if

! adjust levels for free surface
d_zcr=max(1.+w_eta/depth_hl(:,wlev+1),minwater/depth_hl(:,wlev+1))
! calculate melting temperature
call calcmelt(d_timelt,d_zcr)

call getrho(a_ps,d_rho,d_rs,d_alpha,d_beta,d_zcr)                                                    ! equation of state
! ice formation is called for both calcprog=.true. and .false.
call mlonewice(dt,a_ps,d_rho,d_timelt,d_zcr,diag)                                                    ! create new ice

! split adjustment of free surface and ice thickness to ensure conservation
d_neta=w_eta+0.001*dt*(a_inflow+(1.-i_fracice)*(a_rnd+a_snd))
d_ndsn=i_dsn
d_ndic=i_dic
d_nsto=i_sto
where (i_fracice>0.)
  d_ndsn=d_ndsn+0.001*dt*(a_rnd+a_snd)
end where

call fluxcalc(dt,a_u,a_v,a_temp,a_qg,a_ps,a_zmin,a_zmins,d_rho,d_zcr,d_neta,d_oldu,d_oldv,diag)      ! water fluxes
call getwflux(a_sg,a_rg,a_rnd,a_snd,a_vnratio,a_fbvis,a_fbnir,a_inflow,d_rho,d_rs,d_nsq,d_rad,     &
              d_alpha,d_beta,d_b0,d_ustar,d_wu0,d_wv0,d_wt0,d_ws0,d_zcr)                             ! boundary conditions
call iceflux(dt,a_sg,a_rg,a_rnd,a_snd,a_vnratio,a_fbvis,a_fbnir,a_u,a_v,a_temp,a_qg,a_ps,a_zmin,   &
             a_zmins,d_rho,d_ftop,d_tb,d_fb,d_timelt,d_nk,d_tauxicw,d_tauyicw,d_zcr,d_ndsn,d_ndic, &
             d_nsto,d_oldu,d_oldv,diag)                                                              ! ice fluxes
if (calcprog) then
  call mloice(dt,a_rnd,a_snd,a_ps,d_alpha,d_beta,d_b0,d_wu0,d_wv0,d_wt0,d_ws0,d_ftop,d_tb,d_fb,    &
              d_timelt,d_tauxicw,d_tauyicw,d_ustar,d_rho,d_rs,d_nk,d_neta,d_ndsn,d_ndic,d_nsto,    &
              d_zcr,diag)                                                                            ! update ice
  call mlocalc(dt,a_f,d_rho,d_nsq,d_rad,d_alpha,d_beta,d_b0,d_ustar,d_wu0,d_wv0,d_wt0,d_ws0,d_zcr, &
               d_neta,diag)                                                                          ! update water
end if 
call scrncalc(a_u,a_v,a_temp,a_qg,a_ps,a_zmin,a_zmins,diag)                                          ! screen diagnostics

workb=emisice**0.25*i_tsurf
workc=(1.-i_fracice)/log(a_zmin/p_zo)**2+i_fracice/log(a_zmin/p_zoice)**2
fracice=0.
siced=0.
sst(wmap)    =(1.-i_fracice)*w_temp(:,1)+i_fracice*workb
zo(wmap)     =a_zmin*exp(-1./sqrt(workc))
cd(wmap)     =(1.-i_fracice)*p_cd +i_fracice*p_cdice
cds(wmap)    =(1.-i_fracice)*p_cds+i_fracice*p_cdsice
fg(wmap)     =(1.-i_fracice)*p_fg +i_fracice*p_fgice
eg(wmap)     =(1.-i_fracice)*p_eg +i_fracice*p_egice
wetfac(wmap) =(1.-i_fracice)      +i_fracice*p_wetfacice
epan(wmap)   =p_eg
epot(wmap)   =(1.-i_fracice)*p_eg +i_fracice*p_egice/p_wetfacice
fracice(wmap)=i_fracice
siced(wmap)  =i_dic
snowd(wmap)  =i_dsn

return
end subroutine mloeval

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MLO calcs for water (no ice)

subroutine mlocalc(dt,a_f,d_rho,d_nsq,d_rad,d_alpha,d_beta,d_b0,d_ustar,d_wu0,d_wv0,d_wt0,d_ws0,d_zcr, &
                   d_neta,diag)

implicit none

integer, intent(in) :: diag
integer ii,iqw
real, intent(in) :: dt
real, dimension(wfull,wlev) :: km,ks,gammas
real, dimension(wfull,wlev) :: rhs
real(kind=8), dimension(wfull,2:wlev) :: aa
real(kind=8), dimension(wfull,wlev) :: bb,dd
real(kind=8), dimension(wfull,1:wlev-1) :: cc
real, dimension(wfull,wlev), intent(in) :: d_rho,d_nsq,d_rad,d_alpha,d_beta
real, dimension(wfull) :: xp,xm,dumt0,umag
real, dimension(wfull), intent(in) :: a_f
real, dimension(wfull), intent(inout) :: d_b0,d_ustar,d_wu0,d_wv0,d_wt0,d_ws0,d_zcr,d_neta

call getmixdepth(d_rho,d_nsq,d_rad,d_alpha,d_beta,d_b0,d_ustar,a_f,d_zcr) ! solve for mixed layer depth (calculated at full levels)
call getstab(km,ks,gammas,d_nsq,d_ustar,d_zcr) ! solve for stability functions and non-local term (calculated at half levels)

! POTENTIAL TEMPERATURE
if (incradgam>0) then
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
where (depth_hl(:,wlev+1)<mxd)
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

! adjust surface height
select case(deprelax)
  case(0) ! free surface height
    w_eta=d_neta
  case(1) ! relax surface height
    w_eta=d_neta-dt*d_neta/(3600.*24.)
  case(2) ! fix surface height
    w_eta=0.
  case DEFAULT
    write(6,*) "ERROR: Invalid deprelax ",deprelax
    stop
end select

return
end subroutine mlocalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! tri-diagonal solver

subroutine thomas(outo,aai,bbi,cci,ddi)

implicit none

integer ii
real(kind=8), dimension(wfull,2:wlev), intent(in) :: aai
real(kind=8), dimension(wfull,wlev), intent(in) :: bbi,ddi
real(kind=8), dimension(wfull,wlev-1), intent(in) :: cci
real, dimension(wfull,wlev), intent(out) :: outo
real(kind=8), dimension(wfull,wlev) :: cc,dd
real(kind=8), dimension(wfull) :: n

cc(:,1)=cci(:,1)/bbi(:,1)
dd(:,1)=ddi(:,1)/bbi(:,1)

do ii=2,wlev-1
  n=bbi(:,ii)-cc(:,ii-1)*aai(:,ii)
  cc(:,ii)=cci(:,ii)/n
  dd(:,ii)=(ddi(:,ii)-dd(:,ii-1)*aai(:,ii))/n
end do
n=bbi(:,wlev)-cc(:,wlev-1)*aai(:,wlev)
dd(:,wlev)=(ddi(:,wlev)-dd(:,wlev-1)*aai(:,wlev))/n
outo(:,wlev)=dd(:,wlev)
do ii=wlev-1,1,-1
  outo(:,ii)=dd(:,ii)-cc(:,ii)*outo(:,ii+1)
end do

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
do ii=2,wlev
  ! s - shear
  where (ri(:,ii)<0.)
    num(:,ii)=nu0
  elsewhere (ri(:,ii)<ri0)
    num(:,ii)=nu0*(1.-(ri(:,ii)/ri0)**2)**3
  elsewhere
    num(:,ii)=0.
  endwhere
  ! d - double-diffusive
  ! double-diffusive mixing is neglected for now
  ! w - internal ave
  num(:,ii)=num(:,ii)+numw
  nus(:,ii)=num(:,ii)+nusw ! note that nus is a copy of num
end do
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
  if (p_mixdepth(iqw)>d_depth_hl(iqw)) mixind_hl(iqw)=jj
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
  dwm1ds(iqw)=-5.*max(p_bf(iqw),0.)/d_ustar(iqw)**4
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
do ii=2,wlev
  where (ii<=mixind_hl)
    sigma=depth_hl(:,ii)*d_zcr/p_mixdepth
    km(:,ii)=max(p_mixdepth*wm(:,ii)*sigma*(1.+sigma*(a2m+a3m*sigma)),num(:,ii))
    ks(:,ii)=max(p_mixdepth*ws(:,ii)*sigma*(1.+sigma*(a2s+a3s*sigma)),nus(:,ii))
  elsewhere
    km(:,ii)=num(:,ii)
    ks(:,ii)=nus(:,ii)
  end where
end do

!--------------------------------------------------------------------
! non-local term
! gammas is the same for temp and sal when double-diffusion is not employed
cg=10.*vkar*(98.96*vkar*epsilon)**(1./3.) ! Large (1994)
!cg=5.*vkar*(98.96*vkar*epsilon)**(1./3.) ! Bernie (2004)
do ii=2,wlev
  where (p_bf<0..and.ii<=mixind_hl) ! unstable
    gammas(:,ii)=cg/max(ws(:,ii)*p_mixdepth,1.E-20)
  elsewhere
    gammas(:,ii)=0.
  end where
end do

return
end subroutine getstab

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate mixed layer depth

subroutine getmixdepth(d_rho,d_nsq,d_rad,d_alpha,d_beta,d_b0,d_ustar,a_f,d_zcr)

implicit none

integer ii,jj,kk,iqw
integer, dimension(wfull) :: isf
real vtc,dvsq,vtsq,xp
real tnsq,tws,twu,twv,tdepth,tbuoy,trho
real oldxp,oldtrib,trib,newxp,deldz,aa,bb,dsf
real, dimension(wfull,wlev) :: ws,wm,dumbuoy,rib
real, dimension(wfull) :: dumbf,l,d_depth,usf,vsf,rsf
real, dimension(wfull,wlev), intent(in) :: d_rho,d_nsq,d_rad,d_alpha,d_beta
real, dimension(wfull), intent(in) :: d_b0,d_ustar,d_zcr
real, dimension(wfull), intent(in) :: a_f
integer, parameter :: maxits = 3

vtc=1.8*sqrt(0.2/(98.96*epsilon))/(vkar*vkar*ric)

! Modify buoyancy forcing with solar radiation
if (incradbf>0) then
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

! Estimate surface layer values
usf=0.
vsf=0.
rsf=0.
dsf=0.
isf=wlev-1
do iqw=1,wfull
  dsf=0.
  do ii=1,wlev-1
    aa=dz(iqw,ii)*d_zcr(iqw)
    bb=max(minsfc-dsf,0.)
    deldz=min(aa,bb)
    usf(iqw)=usf(iqw)+w_u(iqw,ii)*deldz
    vsf(iqw)=vsf(iqw)+w_v(iqw,ii)*deldz
    rsf(iqw)=rsf(iqw)+d_rho(iqw,ii)*deldz
    dsf=dsf+deldz
    if (bb<=0.) exit
  end do
  if (depth(iqw,ii)*d_zcr(iqw)>minsfc) then
    isf(iqw)=ii
  else
    isf(iqw)=ii+1
  end if
  usf(iqw)=usf(iqw)/dsf
  vsf(iqw)=vsf(iqw)/dsf
  rsf(iqw)=rsf(iqw)/dsf
end do

! Calculate local buoyancy
dumbuoy=0.
do iqw=1,wfull
  do ii=isf(iqw),wlev
    dumbuoy(iqw,ii)=grav*(d_rho(iqw,ii)-rsf(iqw))
  end do
end do

! Calculate mixed layer depth from critical Ri
p_mixind=wlev-1
p_mixdepth=depth(:,wlev)*d_zcr
rib=0.
do iqw=1,wfull
  do ii=isf(iqw),wlev
    jj=min(ii+1,wlev)
    vtsq=depth(iqw,ii)*d_zcr(iqw)*ws(iqw,ii)*sqrt(0.5*max(d_nsq(iqw,ii)+d_nsq(iqw,jj),0.))*vtc
    dvsq=(usf(iqw)-w_u(iqw,ii))**2+(vsf(iqw)-w_v(iqw,ii))**2
    rib(iqw,ii)=(depth(iqw,ii)*d_zcr(iqw)-minsfc)*dumbuoy(iqw,ii)/(max(dvsq+vtsq,1.E-20)*d_rho(iqw,ii))
    if (rib(iqw,ii)>ric) then
      jj=max(ii-1,1)
      p_mixind(iqw)=jj
      xp=min(max((ric-rib(iqw,jj))/max(rib(iqw,ii)-rib(iqw,jj),1.E-20),0.),1.)
      p_mixdepth(iqw)=((1.-xp)*depth(iqw,jj)+xp*depth(iqw,ii))*d_zcr(iqw)
      exit
    end if
  end do 
end do

! Refine mixed-layer-depth calculation by improving vertical profile of buoyancy
if (mixmeth==1) then
  do iqw=1,wfull
    ii=p_mixind(iqw)+1
    jj=min(ii+1,wlev)
    oldxp=0.
    oldtrib=rib(iqw,ii-1)
    xp=(p_mixdepth(iqw)-depth(iqw,ii-1)*d_zcr(iqw))/((depth(iqw,ii)-depth(iqw,ii-1))*d_zcr(iqw))
    do kk=1,maxits
      if (xp<0.5) then
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
      dvsq=(usf(iqw)-twu)**2+(vsf(iqw)-twv)**2
      trib=(tdepth-minsfc)*tbuoy/(max(dvsq+vtsq,1.E-20)*trho)
      if (abs(trib-oldtrib)<1.E-5) then
        exit
      end if
      newxp=xp-(trib-ric)*(xp-oldxp)/(trib-oldtrib) ! i.e., (trib-ric-oldtrib+ric)
      oldtrib=trib
      oldxp=xp
      xp=newxp
      xp=min(max(xp,0.),1.)
    end do
    p_mixdepth(iqw)=((1.-xp)*depth(iqw,ii-1)+xp*depth(iqw,ii))*d_zcr(iqw)
  end do
end if

! calculate buoyancy forcing
call getbf(d_rad,d_alpha,d_b0)

! impose limits for stable conditions
where(p_bf>1.E-10.and.abs(a_f)>1.E-10)
  l=0.7*d_ustar/abs(a_f)
elsewhere (p_bf>1.E-10)
  l=d_ustar*d_ustar*d_ustar/(vkar*p_bf)
elsewhere
  l=depth(:,wlev)*d_zcr
end where
p_mixdepth=min(p_mixdepth,l)
p_mixdepth=max(p_mixdepth,depth(:,1)*d_zcr)
p_mixdepth=min(p_mixdepth,depth(:,wlev)*d_zcr)

! recalculate index for mixdepth
p_mixind=wlev-1
do iqw=1,wfull
  do ii=2,wlev
    if (depth(iqw,ii)*d_zcr(iqw)>p_mixdepth(iqw)) then
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

if (incradbf>0) then
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
where (bf<=0..and.sig>epsilon) ! unstable
  sig=epsilon
end where
uuu=d_ustar**3
invl=vkar*bf      ! invl = ustar**3/L or L=ustar**3/(vkar*bf)
zeta=sig*mixdp*invl
zeta=min(zeta,uuu) ! MJT suggestion

where (zeta>0.)
  wm=vkar*d_ustar*uuu/(uuu+5.*zeta)
elsewhere (zeta>zetam*uuu)
  wm=vkar*(d_ustar*uuu-16.*d_ustar*zeta)**(1./4.)
elsewhere
  wm=vkar*(am*uuu-cm*zeta)**(1./3.)
end where

where (zeta>0.)
  ws=vkar*d_ustar*uuu/(uuu+5.*zeta)
elsewhere (zeta>zetas*uuu)
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
real, dimension(wfull,wlev) :: d_dz
real, dimension(wfull,wlev), intent(inout) :: d_rho,d_rs,d_alpha,d_beta
real, dimension(wfull), intent(in) :: a_ps
real, dimension(wfull), intent(inout) :: d_zcr

pxtr=a_ps+grav*i_fracice*(i_dic*rhoic+i_dsn*rhosn)
do ii=1,wlev
  d_dz(:,ii)=dz(:,ii)*d_zcr
end do
call calcdensity(d_rho,d_alpha,d_beta,d_rs,rho0,w_temp,w_sal,d_dz,pxtr)

return
end subroutine getrho

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate surface density (used for sea-ice melting)

subroutine getrho1(salin,a_ps,d_rho,d_zcr)

implicit none

integer ii
real, dimension(wfull) :: rho0,pxtr
real, dimension(wfull,1) :: d_dz,d_rs,d_alpha,d_beta,d_isal
real, dimension(wfull,1), intent(inout) :: d_rho
real, dimension(wfull), intent(in) :: a_ps,salin
real, dimension(wfull), intent(inout) :: d_zcr

pxtr=a_ps+grav*i_fracice*(i_dic*rhoic+i_dsn*rhosn)
d_dz(:,1)=dz(:,1)*d_zcr
d_isal(:,1)=salin
call calcdensity(d_rho(:,1:1),d_alpha(:,1:1),d_beta(:,1:1),d_rs(:,1:1),rho0,w_temp(:,1:1),d_isal(:,1:1),d_dz(:,1:1),pxtr)

return
end subroutine getrho1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate water boundary conditions

subroutine getwflux(a_sg,a_rg,a_rnd,a_snd,a_vnratio,a_fbvis,a_fbnir,a_inflow,d_rho,d_rs,d_nsq,d_rad,d_alpha, &
                    d_beta,d_b0,d_ustar,d_wu0,d_wv0,d_wt0,d_ws0,d_zcr)

implicit none

integer ii
real, dimension(wfull) :: visalb,niralb
real, dimension(wfull,wlev), intent(inout) :: d_rho,d_rs,d_nsq,d_rad,d_alpha,d_beta
real, dimension(wfull), intent(inout) :: d_b0,d_ustar,d_wu0,d_wv0,d_wt0,d_ws0,d_zcr
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

! Boundary conditions
d_wu0=-(1.-i_fracice)*p_taux/d_rho(:,1)
d_wv0=-(1.-i_fracice)*p_tauy/d_rho(:,1)
d_wt0=-(1.-i_fracice)*(-p_fg-p_eg+a_rg-sbconst*w_temp(:,1)**4)/(d_rho(:,1)*cp0)
d_wt0=d_wt0+(1.-i_fracice)*lf*a_snd/(d_rho(:,1)*cp0) ! melting snow
d_ws0=(1.-i_fracice)*(a_rnd+a_snd-p_eg/lv)*w_sal(:,1)/d_rs(:,1)
d_ws0=d_ws0+(a_inflow)*w_sal(:,1)/d_rs(:,1)

d_ustar=max(sqrt(sqrt(d_wu0*d_wu0+d_wv0*d_wv0)),1.E-6)
d_b0=-grav*(d_alpha(:,1)*d_wt0-d_beta(:,1)*d_ws0) ! -ve sign to compensate for sign of wt0 and ws0
                                                  ! This is the opposite sign used by Large for wb0, but
                                                  ! is same sign as Large used for Bf (+ve stable, -ve unstable)

return
end subroutine getwflux

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate density

subroutine calcdensity(d_rho,d_alpha,d_beta,rs,rho0,tt,ss,ddz,pxtr)

implicit none

integer full,wlx,i,ii
!integer, parameter :: nits=1 ! iterate for density (nits=1 recommended)
real, dimension(:,:), intent(in) :: tt
real, dimension(size(tt,1),size(tt,2)), intent(in) :: ss,ddz
real, dimension(size(tt,1),size(tt,2)), intent(out) :: d_rho,d_alpha,d_beta,rs
real, dimension(size(tt,1)), intent(in) :: pxtr
real, dimension(size(tt,1)), intent(out) :: rho0
real, dimension(size(tt,1)) :: t,s,p1,p2,t2,t3,t4,t5,s2,s3,s32,ptot
real, dimension(size(tt,1)) :: drho0dt,drho0ds,dskdt,dskds,sk,sks
real, dimension(size(tt,1)) :: drhodt,drhods,rs0
real, parameter :: density = 1035.

full=size(tt,1)
wlx=size(tt,2)
d_rho=density

t = min(max(tt(:,1)-273.16,-2.2),100.)
s = min(max(ss(:,1),0.),maxsal) ! limit max salinity for equation of state
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
  do ii=1,wlx
    t = min(max(tt(:,ii)-273.16,-2.2),100.)
    s = min(max(ss(:,ii),0.),maxsal)
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
    
    sks = 1.965933e4 + 1.444304e2*t(:) - 1.706103*t2(:)  &
                + 9.648704e-3*t3(:)  - 4.190253e-5*t4(:) &
                + p1(:)*(3.186519 + 2.212276e-2*t(:)     &
                - 2.984642e-4*t2(:) + 1.956415e-6*t3(:)) &
                + p2(:)*(2.102898e-4 - 1.202016e-5*t(:)  &
                + 1.394680e-7*t2(:)) ! sal=0.
    sk=sks+ s(:)*(52.84855 - 3.101089e-1*t(:)                   &
                + 6.283263e-3*t2(:) -5.084188e-5*t3(:))         &
                + s32(:)*(3.886640e-1 + 9.085835e-3*t(:)        &
                - 4.619924e-4*t2(:))                            &
                + p1(:)*s(:)*(6.704388e-3  -1.847318e-4*t(:)    &
                + 2.059331e-7*t2(:)) + 1.480266e-4*p1(:)*s32(:) &
                +p2(:)*s(:)*(-2.040237e-6 &
                + 6.128773e-8*t(:) + 6.207323e-10*t2(:)) ! + sal terms             
    dskdt= 1.444304e2 - 2.*1.706103*t(:)                       &
                + 3.*9.648704e-3*t2(:)  - 4.*4.190253e-5*t3(:) &
                + s(:)*( - 3.101089e-1                         &
                + 2.*6.283263e-3*t(:) -3.*5.084188e-5*t2(:))   &
                + s32(:)*(9.085835e-3                          &
                - 2.*4.619924e-4*t(:))                         &
                + p1(:)*(2.212276e-2                           &
                - 2.*2.984642e-4*t(:) + 3.*1.956415e-6*t2(:))  &
                + p1(:)*s(:)*(-1.847318e-4                     &
                + 2.*2.059331e-7*t(:))                         &
                + p2(:)*(- 1.202016e-5                         &
                + 2.*1.394680e-7*t(:)) +p2(:)*s(:)*(           &
                + 6.128773e-8 + 2.*6.207323e-10*t(:))
    dskds=(52.84855 - 3.101089e-1*t(:)                                  &
                + 6.283263e-3*t2(:) -5.084188e-5*t3(:))                 &
                + 1.5*sqrt(s(:))*(3.886640e-1 + 9.085835e-3*t(:)        &
                - 4.619924e-4*t2(:))                                    &
                + p1(:)*(6.704388e-3  -1.847318e-4*t(:)                 &
                + 2.059331e-7*t2(:)) + 1.5*1.480266e-4*p1(:)*sqrt(s(:)) &
                +p2(:)*(-2.040237e-6                                    &
                + 6.128773e-8*t(:) + 6.207323e-10*t2(:))
!    dskdp=(3.186519 + 2.212276e-2*t(:)                   &
!                - 2.984642e-4*t2(:) + 1.956415e-6*t3(:)) &
!                + 2.*p1*(2.102898e-4 - 1.202016e-5*t(:)  &
!                + 1.394680e-7*t2(:))                     &
!                + s(:)*(6.704388e-3  -1.847318e-4*t(:)   &
!                + 2.059331e-7*t2(:))                     &
!                + 1.480266e-4*s32(:)                     &
!                + 2.*p1(:)*s(:)*(-2.040237e-6            &
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
! Estimate melting point
subroutine calcmelt(d_timelt,d_zcr)

implicit none

integer iqw,ii
real, dimension(wfull), intent(out) :: d_timelt
real, dimension(wfull), intent(in) :: d_zcr
real, dimension(wfull) :: ssf
real dsf,aa,bb,deldz
logical lflag

ssf=0.
do iqw=1,wfull
  dsf=0.
  do ii=1,wlev-1
    aa=dz(iqw,ii)*d_zcr(iqw)
    bb=minsfc-dsf
    lflag=aa>=bb
    deldz=min(aa,bb)
    ssf(iqw)=ssf(iqw)+w_sal(iqw,ii)*deldz
    dsf=dsf+deldz
    if (lflag) exit
  end do
  ssf(iqw)=ssf(iqw)/dsf
end do
d_timelt=273.16-0.054*min(max(ssf,0.),maxsal) ! ice melting temperature from CICE

return
end subroutine calcmelt

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
! Calculate fluxes between MLO and atmosphere

subroutine fluxcalc(dt,a_u,a_v,a_temp,a_qg,a_ps,a_zmin,a_zmins,d_rho, &
                    d_zcr,d_neta,d_oldu,d_oldv,diag)

implicit none

integer, intent(in) :: diag
integer it
real, intent(in) :: dt
real, dimension(wfull,wlev), intent(inout) :: d_rho
real, dimension(wfull), intent(in) :: a_u,a_v,a_temp,a_qg,a_ps,a_zmin,a_zmins
real, dimension(wfull), intent(inout) :: d_zcr,d_neta,d_oldu,d_oldv
real, dimension(wfull) :: qsat,dqdt,ri,vmag,rho,srcp
real, dimension(wfull) :: fm,fh,fq,con,consea,afroot,af,daf
real, dimension(wfull) :: den,dfm,dden,dcon,sig,factch,root
real, dimension(wfull) :: aft,afq,atu,atv,dcs,facqch
real, dimension(wfull) :: vmagn,egmax,d_wavail
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

sig=exp(-grav*a_zmins/(rdry*a_temp))
srcp=sig**(rdry/cp)
atu=a_u-fluxwgt*w_u(:,1)-(1.-fluxwgt)*d_oldu
atv=a_v-fluxwgt*w_v(:,1)-(1.-fluxwgt)*d_oldv
vmag=sqrt(max(atu*atu+atv*atv,0.))
vmagn=max(vmag,0.01)
rho=a_ps/(rdry*w_temp(:,1))
ri=min(grav*(a_zmin*a_zmin/a_zmins)*(1.-w_temp(:,1)*srcp/a_temp)/vmagn**2,rimax)

call getqsat(qsat,dqdt,w_temp(:,1),a_ps)
if (zomode==0) then ! CSIRO9
  qsat=0.98*qsat ! with Zeng 1998 for sea water
end if

select case(zomode)
  case(0,1) ! Charnock
    consea=vmag*charnck/grav
    p_zo=0.001    ! first guess
    do it=1,4
      afroot=vkar/log(a_zmin/p_zo)
      af=afroot*afroot
      daf=2.*af*afroot/(vkar*p_zo)
      where (ri>=0.) ! stable water points                                     
        fm=1./(1.+bprm*ri)**2
        con=consea*fm*vmag
        p_zo=p_zo-(p_zo-con*af)/(1.-con*daf)
      elsewhere     ! unstable water points
        den=1.+af*cms*2.*bprm*sqrt(-ri*a_zmin/p_zo)
        fm=1.-2.*bprm*ri/den
        con=consea*fm*vmag
        dden=daf*cms*2.*bprm*sqrt(-ri*a_zmin/p_zo)+af*cms*bprm*sqrt(-ri)*a_zmin/(sqrt(a_zmin/p_zo)*p_zo*p_zo)
        dfm=2.*bprm*ri*dden/(den*den)
        dcon=consea*dfm*vmag
        p_zo=p_zo-(p_zo-con*af)/(1.-dcon*af-con*daf)
      end where
      p_zo=min(max(p_zo,1.5e-5),6.)
    enddo    ! it=1,4
  case(2) ! Beljaars
    p_zo=0.001    ! first guess
    do it=1,4
      afroot=vkar/log(a_zmin/p_zo)
      af=afroot*afroot
      daf=2.*af*afroot/(vkar*p_zo)
      where (ri>=0.) ! stable water points
        fm=1./(1.+bprm*ri)**2
        consea=zcom1*vmag*vmag*fm*af/grav+zcom2*gnu/max(vmag*sqrt(fm*af),gnu)
        dcs=(zcom1*vmag*vmag/grav-0.5*zcom2*gnu/(max(vmag*sqrt(fm*af),gnu)*fm*af))*(fm*daf)
      elsewhere     ! unstable water points
        con=cms*2.*bprm*sqrt(-ri*a_zmin/p_zo)
        den=1.+af*con
        fm=1.-2.*bprm*ri/den
        dfm=2.*bprm*ri*(con*daf+af*cms*bprm*sqrt(-ri)*a_zmin/(sqrt(a_zmin/p_zo)*p_zo*p_zo))/(den*den) ! MJT suggestion
        consea=zcom1*vmag*vmag*af*fm/grav+zcom2*gnu/max(vmag*sqrt(fm*af),gnu)
        dcs=(zcom1*vmag*vmag/grav-0.5*zcom2*gnu/(max(vmag*sqrt(fm*af),gnu)*fm*af))*(fm*daf+dfm*af)
      end where
      p_zo=p_zo-(p_zo-consea)/(1.-dcs)      
      p_zo=min(max(p_zo,1.5e-5),6.)
    enddo    ! it=1,4
end select
afroot=vkar/log(a_zmin/p_zo)
af=afroot*afroot

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
    aft=vkar*vkar/(log(a_zmins/p_zo)*log(a_zmins/p_zoh))
    afq=vkar*vkar/(log(a_zmins/p_zo)*log(a_zmins/p_zoq))
end select

where (ri>=0.)
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

! update drag
p_cd=af*fm
p_cds=aft*fh
p_cdq=afq*fq

! turn off lake evaporation when minimum depth is reached
! fg should be replaced with bare ground value
d_wavail=max(depth_hl(:,wlev+1)+d_neta-minwater,0.)
egmax=1000.*lv*d_wavail/(dt*max(1.-i_fracice,0.001))

! explicit estimate of fluxes
! (replace with implicit scheme if water becomes too shallow)
p_fg=rho*p_cds*cp*vmagn*(w_temp(:,1)-a_temp/srcp)
p_eg=min(rho*p_cdq*lv*vmagn*(qsat-a_qg),egmax)
p_taux=rho*p_cd*vmagn*atu
p_tauy=rho*p_cd*vmagn*atv

! update free surface after evaporation
d_neta=d_neta-0.001*dt*(1.-i_fracice)*p_eg/lv

return
end subroutine fluxcalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Estimate saturation mixing ratio (from CCAM)
! following version of getqsat copes better with T < 0C over ice
subroutine getqsat(qsat,dqdt,temp,ps)

implicit none

real, dimension(wfull), intent(in) :: temp,ps
real, dimension(wfull), intent(out) :: qsat,dqdt
real, dimension(0:220), save :: table
real, dimension(wfull) :: esatf,tdiff,dedt,rx
integer, dimension(wfull) :: ix
logical, save :: first=.true.

if (first) then
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
  first=.false.
end if

tdiff=min(max( temp-123.16, 0.), 219.)
rx=tdiff-aint(tdiff)
ix=int(tdiff)
esatf=(1.-rx)*table(ix)+ rx*table(ix+1)
qsat=0.622*esatf/max(ps-esatf,0.1)

! method #1
dedt=table(ix+1)-table(ix) ! divide by 1C
dqdt=qsat*dedt*ps/(esatf*max(ps-esatf,0.1))
! method #2
!dqdt=qsat*lv/rvap*qsat*ps/(temp*temp*(ps-0.378*esatf))

return
end subroutine getqsat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Pack sea ice for calcuation

subroutine mloice(dt,a_rnd,a_snd,a_ps,d_alpha,d_beta,d_b0,d_wu0,d_wv0,d_wt0,d_ws0,d_ftop,d_tb,d_fb,d_timelt, &
                  d_tauxicw,d_tauyicw,d_ustar,d_rho,d_rs,d_nk,d_neta,d_ndsn,d_ndic,d_nsto,d_zcr,diag)

implicit none

integer, intent(in) :: diag
integer nice,ii,iqw
real, intent(in) :: dt
real, dimension(wfull), intent(in) :: a_rnd,a_snd,a_ps
real, dimension(wfull) :: ap_rnd,ap_snd
real, dimension(wfull,0:2) :: ip_tn
real, dimension(wfull) :: ip_dic,ip_dsn,ip_fracice,ip_tsurf,ip_sto
real, dimension(wfull) :: pp_egice
real, dimension(wfull,wlev), intent(in) :: d_alpha,d_beta,d_rho,d_rs
real, dimension(wfull), intent(inout) :: d_b0,d_wu0,d_wv0,d_wt0,d_ws0,d_ftop,d_tb,d_fb,d_timelt
real, dimension(wfull), intent(inout) :: d_tauxicw,d_tauyicw,d_ustar,d_neta
real, dimension(wfull), intent(inout) :: d_ndsn,d_ndic,d_nsto,d_zcr
integer, dimension(wfull), intent(inout) :: d_nk
real, dimension(wfull) :: dp_ftop,dp_tb,dp_fb,dp_timelt,dp_salflxf,dp_salflxs
real, dimension(wfull) :: dp_wavail
real, dimension(wfull) :: d_salflxf,d_salflxs,d_wavail
real, dimension(wfull,1) :: d_ri
real, dimension(wfull) :: imass,ofracice,deld
integer, dimension(wfull) :: dp_nk,cmap

d_salflxf=0.
d_salflxs=0.
d_wavail=max(depth_hl(:,wlev+1)+d_neta-minwater,0.)
imass=max(rhoic*i_dic+rhosn*i_dsn,10.) ! ice mass per unit area
ofracice=i_fracice

! identify sea ice for packing
nice=0
do iqw=1,wfull
  if (d_ndic(iqw)>icemin) then
    nice=nice+1
    cmap(nice)=iqw
  end if
end do

if (nice>0) then

  ! pack data
  do ii=0,2
    ip_tn(1:nice,ii)=i_tn(cmap(1:nice),ii)
  end do
  ip_dic(1:nice)=d_ndic(cmap(1:nice))
  ip_dsn(1:nice)=d_ndsn(cmap(1:nice))
  ip_fracice(1:nice)=i_fracice(cmap(1:nice))
  ip_tsurf(1:nice)=i_tsurf(cmap(1:nice))
  ip_sto(1:nice)=d_nsto(cmap(1:nice))
  ap_rnd(1:nice)=a_rnd(cmap(1:nice))
  ap_snd(1:nice)=a_snd(cmap(1:nice))
  pp_egice(1:nice)=p_egice(cmap(1:nice))
  dp_ftop(1:nice)=d_ftop(cmap(1:nice))
  dp_tb(1:nice)=d_tb(cmap(1:nice))
  dp_fb(1:nice)=d_fb(cmap(1:nice))
  dp_timelt(1:nice)=d_timelt(cmap(1:nice))
  dp_salflxf(1:nice)=d_salflxf(cmap(1:nice))
  dp_salflxs(1:nice)=d_salflxs(cmap(1:nice))
  dp_nk(1:nice)=d_nk(cmap(1:nice))
  dp_wavail(1:nice)=d_wavail(cmap(1:nice))
  ! update ice prognostic variables
  call seaicecalc(nice,dt,ip_tn(1:nice,0),ip_tn(1:nice,1),ip_tn(1:nice,2),ip_dic(1:nice),ip_dsn(1:nice),ip_fracice(1:nice), &
                  ip_tsurf(1:nice),ip_sto(1:nice),ap_rnd(1:nice),ap_snd(1:nice),pp_egice(1:nice),dp_ftop(1:nice),           &
                  dp_tb(1:nice),dp_fb(1:nice),dp_timelt(1:nice),dp_salflxf(1:nice),dp_salflxs(1:nice),dp_nk(1:nice),        &
                  dp_wavail(1:nice),diag)
  ! unpack ice data
  do ii=0,2
    i_tn(cmap(1:nice),ii)=ip_tn(1:nice,ii)
  end do
  i_fracice=0.
  i_dic=0.
  i_dsn=0.
  i_dic(cmap(1:nice))=ip_dic(1:nice)
  i_dsn(cmap(1:nice))=ip_dsn(1:nice)
  i_fracice(cmap(1:nice))=ip_fracice(1:nice)
  i_tsurf(cmap(1:nice))=ip_tsurf(1:nice)
  i_sto(cmap(1:nice))=ip_sto(1:nice)
  d_salflxf(cmap(1:nice))=dp_salflxf(1:nice)
  d_salflxs(cmap(1:nice))=dp_salflxs(1:nice)
  d_nk(cmap(1:nice))=dp_nk(1:nice)
end if

! update ice velocities due to stress terms
i_u=i_u+dt*(p_tauxica-d_tauxicw)/imass
i_v=i_v+dt*(p_tauyica-d_tauyicw)/imass

! Remove excessive salinity from ice
where (i_dic>=icemin)
  deld=i_dic*(1.-icesal/max(i_sal,icesal))
  i_sal=i_sal*(1.-deld/i_dic)
elsewhere
  deld=0.
end where
d_salflxf=d_salflxf+deld*rhoic/dt
d_salflxs=d_salflxs-deld*rhoic/dt

! estimate density of water from ice melting
call getrho1(i_sal,a_ps,d_ri,d_zcr)

! update water boundary conditions
d_wu0=d_wu0-ofracice*d_tauxicw/d_rho(:,1)
d_wv0=d_wv0-ofracice*d_tauyicw/d_rho(:,1)
d_wt0=d_wt0+ofracice*d_fb/(d_rho(:,1)*cp0)
d_ws0=d_ws0-ofracice*(d_salflxf*w_sal(:,1)/d_rs(:,1)+d_salflxs*(w_sal(:,1)-i_sal)/d_ri(:,1))
d_ustar=max(sqrt(sqrt(d_wu0*d_wu0+d_wv0*d_wv0)),1.E-6)
d_b0=-grav*(d_alpha(:,1)*d_wt0-d_beta(:,1)*d_ws0) ! -ve sign is to account for sign of wt0 and ws0

! update free surface height with water flux from ice
d_neta=d_neta-dt*ofracice*(d_salflxf+d_salflxs)/rhowt

return
end subroutine mloice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Form seaice before flux calculations

subroutine mlonewice(dt,a_ps,d_rho,d_timelt,d_zcr,diag)

implicit none

integer, intent(in) :: diag
integer iqw,ii
real, intent(in) :: dt
real aa,bb,dsf,deldz,avet,aves,aveto,aveso,delt,dels
real, dimension(wfull,wlev) :: sdic
real, dimension(wfull) :: newdic,newsal,newtn,newsl,gamm,cdic
real, dimension(wfull,wlev), intent(in) :: d_rho
real, dimension(wfull), intent(inout) :: d_timelt,d_zcr
real, dimension(wfull), intent(in) :: a_ps
real, dimension(wfull) :: worka,maxnewice,d_wavail
real, dimension(wfull,1) :: d_ri
logical, dimension(wfull) :: lnewice
logical lflag

! limits on ice formation
d_wavail=max(depth_hl(:,wlev+1)+w_eta-minwater,0.)
where (i_fracice<1.)
  maxnewice=d_wavail*rhowt/rhoic/(1.-i_fracice)
elsewhere
  maxnewice=0.
end where
maxnewice=max(min(maxnewice,0.15),0.)

! formation
! Search over all levels to avoid problems with thin surface layers
sdic=0.
do iqw=1,wfull
  dsf=0.
  do ii=1,wlev-1
    aa=dz(iqw,ii)*d_zcr(iqw)
    bb=max(minsfc-dsf,0.)
    deldz=min(aa,bb)
    sdic(iqw,ii)=max(d_timelt(iqw)-w_temp(iqw,ii),0.)*cp0*d_rho(iqw,ii)*deldz/qice
    dsf=dsf+deldz
    if (bb<=0.) exit
  end do
end do
newdic=sum(sdic,2)
newdic=min(maxnewice,newdic)
newsal=min(w_sal(:,1),icesal)
lnewice=newdic>1.1*icemin*fracbreak
call getrho1(newsal,a_ps,d_ri,d_zcr)
cdic=0.
do ii=1,wlev
  sdic(:,ii)=max(min(sdic(:,ii),newdic-cdic),0.)
  cdic=cdic+sdic(:,ii)  
  newtn=w_temp(:,ii)+sdic(:,ii)*qice/(cp0*d_rho(:,ii)*dz(:,ii)*d_zcr)
  newsl=w_sal(:,ii)+sdic(:,ii)*rhoic*(w_sal(:,ii)-newsal)/(d_ri(:,1)*dz(:,ii)*d_zcr)
  where (lnewice)
    w_temp(:,ii)=w_temp(:,ii)*i_fracice+newtn*(1.-i_fracice)
    w_sal(:,ii)=w_sal(:,ii)*i_fracice+newsl*(1.-i_fracice)
  end where
end do
where (lnewice) ! form new sea-ice
  i_dic=i_dic*i_fracice+newdic*(1.-i_fracice)
  i_dsn=i_dsn*i_fracice
  i_tsurf=i_tsurf*i_fracice+w_temp(:,1)*(1.-i_fracice)
  i_tn(:,0)=i_tn(:,0)*i_fracice+i_tsurf*(1.-i_fracice)
  i_tn(:,1)=i_tn(:,1)*i_fracice+i_tsurf*(1.-i_fracice)
  i_tn(:,2)=i_tn(:,2)*i_fracice+i_tsurf*(1.-i_fracice)
  i_sto=i_sto*i_fracice
  i_sal=i_sal*i_fracice+newsal*(1.-i_fracice)
  w_eta=w_eta-newdic*(1.-i_fracice)*rhoic/rhowt
  i_fracice=1. ! this will be adjusted for fracbreak below  
endwhere

! 1D model of ice break-up
where (i_dic<icebreak.and.i_fracice>fracbreak)
  i_fracice=i_fracice*i_dic/icebreak
  i_dic=icebreak
end where
if (onedice==1) then
  where (i_fracice<1..and.i_dic>icebreak)
     worka=min(i_dic/icebreak,1./max(i_fracice,fracbreak))
     worka=max(worka,1.)
     i_fracice=i_fracice*worka
     i_dic=i_dic/worka
  end where
  i_fracice=min(max(i_fracice,0.),1.)
end if

! removal
call getrho1(i_sal,a_ps,d_ri,d_zcr)
gamm=(gammi*i_dic+gamms*i_dsn)/max(i_dic+i_dsn,1.E-6) ! for energy conservation
do iqw=1,wfull
  if (i_dic(iqw)<=icemin.and.i_fracice(iqw)>0.) then

    avet=0.
    aves=0.
    dsf=0.
    do ii=1,wlev
      aa=dz(iqw,ii)*d_zcr(iqw)
      bb=max(minsfc-dsf,0.)
      deldz=min(aa,bb)
      avet=avet+w_temp(iqw,ii)*deldz
      aves=aves+w_sal(iqw,ii)*deldz
      dsf=dsf+deldz
      if (bb<=0.) exit
    end do
    avet=avet/dsf
    aves=aves/dsf
    aveto=avet
    aveso=aves
  
    w_eta(iqw)=w_eta(iqw)+i_dic(iqw)*i_fracice(iqw)*rhoic/rhowt &
              +i_dsn(iqw)*i_fracice(iqw)*rhosn/rhowt
    aves=aves/(1.+i_dsn(iqw)*rhowt/(d_rho(iqw,1)*dsf))
    aves=(aves+i_dic(iqw)*i_sal(iqw)*rhoic/(d_ri(iqw,1)*dsf)) &
        /(1.+i_dic(iqw)*rhoic/(d_ri(iqw,1)*dsf))
    avet=(avet+i_fracice(iqw)*gamm(iqw)*i_tsurf(iqw)/(cp0*d_rho(iqw,1)*dsf)) &
        /(1.+i_fracice(iqw)*gamm(iqw)/(cp0*d_rho(iqw,1)*dsf))
    avet=avet-i_fracice(iqw)*i_dic(iqw)*qice/(cp0*d_rho(iqw,1)*dsf)
    avet=avet-i_fracice(iqw)*i_dsn(iqw)*qsnow/(cp0*d_rho(iqw,1)*dsf)
    avet=avet+i_fracice(iqw)*i_sto(iqw)/(cp0*d_rho(iqw,1)*dsf)
    
    delt=avet-aveto
    dels=aves-aveso
    dsf=0.
    do ii=1,wlev
      aa=dz(iqw,ii)*d_zcr(iqw)
      bb=max(minsfc-dsf,0.)
      deldz=min(aa,bb)
      w_temp(iqw,ii)=w_temp(iqw,ii)+delt*deldz/aa
      w_sal(iqw,ii)=w_sal(iqw,ii)+dels*deldz/aa      
      dsf=dsf+deldz
      if (bb<=0.) exit
    end do

    i_fracice(iqw)=0.
    i_dic(iqw)=0.
    i_dsn(iqw)=0.
    i_sto(iqw)=0.
    i_sal(iqw)=0.
  end if
end do

return
end subroutine mlonewice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update sea ice prognostic variables

subroutine seaicecalc(nice,dt,ip_tn0,ip_tn1,ip_tn2,ip_dic,ip_dsn,ip_fracice,ip_tsurf,ip_sto,ap_rnd, &
                      ap_snd,pp_egice,dp_ftop,dp_tb,dp_fb,dp_timelt,dp_salflxf,dp_salflxs,dp_nk,    &
                      dp_wavail,diag)

implicit none

integer, intent(in) :: nice,diag
integer iqi,nc0,nc1,nc2,nc3,ii
real, intent(in) :: dt
real, dimension(nice) :: xxx,excess
real, dimension(nice), intent(inout) :: ip_tn0,ip_tn1,ip_tn2
real, dimension(nice), intent(inout) :: ip_dic,ip_dsn,ip_fracice,ip_tsurf,ip_sto
real, dimension(nice) :: it_tn0,it_tn1,it_tn2
real, dimension(nice) :: it_dic,it_dsn,it_fracice,it_tsurf,it_sto
real, dimension(nice), intent(inout) :: dp_ftop,dp_tb,dp_fb,dp_timelt,dp_salflxf,dp_salflxs
real, dimension(nice), intent(inout) :: dp_wavail
integer, dimension(nice), intent(inout) :: dp_nk
real, dimension(nice) :: dt_ftop,dt_tb,dt_fb,dt_timelt,dt_salflxf,dt_salflxs,dt_wavail
real, dimension(nice), intent(in) :: ap_rnd,ap_snd
real, dimension(nice), intent(in) :: pp_egice
real, dimension(nice) :: pt_egice
integer, dimension(nice) :: dt_nk,pq0map,pq1map,pq2map,pq3map

! Pack different ice configurations
nc0=0
nc1=0
nc2=0
nc3=0
do iqi=1,nice
  if (dp_nk(iqi)>0.and.ip_dsn(iqi)>icemin) then        ! snow + 1-2 ice layers
    nc0=nc0+1
    pq0map(nc0)=iqi
  else if (dp_nk(iqi)>0.and.ip_dsn(iqi)<=icemin) then  ! no snow + 1-2 ice layers
    nc1=nc1+1
    pq1map(nc1)=iqi
  else if (dp_nk(iqi)<=0.and.ip_dsn(iqi)>icemin) then  ! snow + 1 ice layer
    nc2=nc2+1
    pq2map(nc2)=iqi
  else if (dp_nk(iqi)<=0.and.ip_dsn(iqi)<=icemin) then ! no snow + 1 ice layer
    nc3=nc3+1
    pq3map(nc3)=iqi
  end if
end do

! Update ice and snow temperatures
if (nc0>0) then
  it_tn0(1:nc0)=ip_tn0(pq0map(1:nc0))
  it_tn1(1:nc0)=ip_tn1(pq0map(1:nc0))
  it_tn2(1:nc0)=ip_tn2(pq0map(1:nc0))
  it_dic(1:nc0)=ip_dic(pq0map(1:nc0))
  it_dsn(1:nc0)=ip_dsn(pq0map(1:nc0))
  it_fracice(1:nc0)=ip_fracice(pq0map(1:nc0))
  it_tsurf(1:nc0)=ip_tsurf(pq0map(1:nc0))
  it_sto(1:nc0)=ip_sto(pq0map(1:nc0))
  pt_egice(1:nc0)=pp_egice(pq0map(1:nc0))
  dt_ftop(1:nc0)=dp_ftop(pq0map(1:nc0))
  dt_tb(1:nc0)=dp_tb(pq0map(1:nc0))
  dt_fb(1:nc0)=dp_fb(pq0map(1:nc0))
  dt_timelt(1:nc0)=dp_timelt(pq0map(1:nc0))
  dt_salflxf(1:nc0)=dp_salflxf(pq0map(1:nc0))
  dt_salflxs(1:nc0)=dp_salflxs(pq0map(1:nc0))
  dt_nk(1:nc0)=dp_nk(pq0map(1:nc0))
  dt_wavail(1:nc0)=dp_wavail(pq0map(1:nc0))
  call icetemps2(nc0,dt,it_tn0(1:nc0),it_tn1(1:nc0),it_tn2(1:nc0),it_dic(1:nc0),it_dsn(1:nc0),it_fracice(1:nc0), &
                 it_tsurf(1:nc0),it_sto(1:nc0),dt_ftop(1:nc0),dt_tb(1:nc0),dt_fb(1:nc0),dt_timelt(1:nc0),        &
                 dt_salflxf(1:nc0),dt_salflxs(1:nc0),dt_nk(1:nc0),dt_wavail(1:nc0),pt_egice(1:nc0),diag)
  ip_tn0(pq0map(1:nc0))=it_tn0(1:nc0)
  ip_tn1(pq0map(1:nc0))=it_tn1(1:nc0)
  ip_tn2(pq0map(1:nc0))=it_tn2(1:nc0)
  ip_dic(pq0map(1:nc0))=it_dic(1:nc0)
  ip_dsn(pq0map(1:nc0))=it_dsn(1:nc0)
  ip_fracice(pq0map(1:nc0))=it_fracice(1:nc0)
  ip_tsurf(pq0map(1:nc0))=it_tsurf(1:nc0)
  ip_sto(pq0map(1:nc0))=it_sto(1:nc0)
  dp_salflxf(pq0map(1:nc0))=dt_salflxf(1:nc0)
  dp_salflxs(pq0map(1:nc0))=dt_salflxs(1:nc0)
  dp_nk(pq0map(1:nc0))=dt_nk(1:nc0)
end if
if (nc1>0) then
  it_tn0(1:nc1)=ip_tn0(pq1map(1:nc1))
  it_tn1(1:nc1)=ip_tn1(pq1map(1:nc1))
  it_tn2(1:nc1)=ip_tn2(pq1map(1:nc1))
  it_dic(1:nc1)=ip_dic(pq1map(1:nc1))
  it_dsn(1:nc1)=ip_dsn(pq1map(1:nc1))
  it_fracice(1:nc1)=ip_fracice(pq1map(1:nc1))
  it_tsurf(1:nc1)=ip_tsurf(pq1map(1:nc1))
  it_sto(1:nc1)=ip_sto(pq1map(1:nc1))
  pt_egice(1:nc1)=pp_egice(pq1map(1:nc1))
  dt_ftop(1:nc1)=dp_ftop(pq1map(1:nc1))
  dt_tb(1:nc1)=dp_tb(pq1map(1:nc1))
  dt_fb(1:nc1)=dp_fb(pq1map(1:nc1))
  dt_timelt(1:nc1)=dp_timelt(pq1map(1:nc1))
  dt_salflxf(1:nc1)=dp_salflxf(pq1map(1:nc1))
  dt_salflxs(1:nc1)=dp_salflxs(pq1map(1:nc1))
  dt_nk(1:nc1)=dp_nk(pq1map(1:nc1))
  dt_wavail(1:nc1)=dp_wavail(pq1map(1:nc1))
  call icetempn2(nc1,dt,it_tn0(1:nc1),it_tn1(1:nc1),it_tn2(1:nc1),it_dic(1:nc1),it_dsn(1:nc1),it_fracice(1:nc1), &
                 it_tsurf(1:nc1),it_sto(1:nc1),dt_ftop(1:nc1),dt_tb(1:nc1),dt_fb(1:nc1),dt_timelt(1:nc1),        &
                 dt_salflxf(1:nc1),dt_salflxs(1:nc1),dt_nk(1:nc1),dt_wavail(1:nc1),pt_egice(1:nc1),diag)
  ip_tn0(pq1map(1:nc1))=it_tn0(1:nc1)
  ip_tn1(pq1map(1:nc1))=it_tn1(1:nc1)
  ip_tn2(pq1map(1:nc1))=it_tn2(1:nc1)
  ip_dic(pq1map(1:nc1))=it_dic(1:nc1)
  ip_dsn(pq1map(1:nc1))=it_dsn(1:nc1)
  ip_fracice(pq1map(1:nc1))=it_fracice(1:nc1)
  ip_tsurf(pq1map(1:nc1))=it_tsurf(1:nc1)
  ip_sto(pq1map(1:nc1))=it_sto(1:nc1)
  dp_salflxf(pq1map(1:nc1))=dt_salflxf(1:nc1)
  dp_salflxs(pq1map(1:nc1))=dt_salflxs(1:nc1)
  dp_nk(pq1map(1:nc1))=dt_nk(1:nc1)
end if
if (nc2>0) then
  it_tn1(1:nc2)=ip_tn1(pq2map(1:nc2))
  it_tn2(1:nc2)=ip_tn2(pq2map(1:nc2))
  it_dic(1:nc2)=ip_dic(pq2map(1:nc2))
  it_dsn(1:nc2)=ip_dsn(pq2map(1:nc2))
  it_fracice(1:nc2)=ip_fracice(pq2map(1:nc2))
  it_tsurf(1:nc2)=ip_tsurf(pq2map(1:nc2))
  it_sto(1:nc2)=ip_sto(pq2map(1:nc2))
  pt_egice(1:nc2)=pp_egice(pq2map(1:nc2))
  dt_ftop(1:nc2)=dp_ftop(pq2map(1:nc2))
  dt_tb(1:nc2)=dp_tb(pq2map(1:nc2))
  dt_fb(1:nc2)=dp_fb(pq2map(1:nc2))
  dt_timelt(1:nc2)=dp_timelt(pq2map(1:nc2))
  dt_salflxf(1:nc2)=dp_salflxf(pq2map(1:nc2))
  dt_salflxs(1:nc2)=dp_salflxs(pq2map(1:nc2))
  dt_nk(1:nc2)=dp_nk(pq2map(1:nc2))
  dt_wavail(1:nc2)=dp_wavail(pq2map(1:nc2))
  call icetemps1(nc2,dt,it_tn0(1:nc2),it_tn1(1:nc2),it_tn2(1:nc2),it_dic(1:nc2),it_dsn(1:nc2),it_fracice(1:nc2), &
                 it_tsurf(1:nc2),it_sto(1:nc2),dt_ftop(1:nc2),dt_tb(1:nc2),dt_fb(1:nc2),dt_timelt(1:nc2),        &
                 dt_salflxf(1:nc2),dt_salflxs(1:nc2),dt_nk(1:nc2),dt_wavail(1:nc2),pt_egice(1:nc2),diag)
  ip_tn0(pq2map(1:nc2))=it_tn0(1:nc2)
  ip_tn1(pq2map(1:nc2))=it_tn1(1:nc2)
  ip_tn2(pq2map(1:nc2))=it_tn2(1:nc2)
  ip_dic(pq2map(1:nc2))=it_dic(1:nc2)
  ip_dsn(pq2map(1:nc2))=it_dsn(1:nc2)
  ip_fracice(pq2map(1:nc2))=it_fracice(1:nc2)
  ip_tsurf(pq2map(1:nc2))=it_tsurf(1:nc2)
  ip_sto(pq2map(1:nc2))=it_sto(1:nc2)
  dp_salflxf(pq2map(1:nc2))=dt_salflxf(1:nc2)
  dp_salflxs(pq2map(1:nc2))=dt_salflxs(1:nc2)
  dp_nk(pq2map(1:nc2))=dt_nk(1:nc2)
end if
if (nc3>0) then
  it_tn1(1:nc3)=ip_tn1(pq3map(1:nc3))
  it_tn2(1:nc3)=ip_tn2(pq3map(1:nc3))
  it_dic(1:nc3)=ip_dic(pq3map(1:nc3))
  it_dsn(1:nc3)=ip_dsn(pq3map(1:nc3))
  it_fracice(1:nc3)=ip_fracice(pq3map(1:nc3))
  it_tsurf(1:nc3)=ip_tsurf(pq3map(1:nc3))
  it_sto(1:nc3)=ip_sto(pq3map(1:nc3))
  pt_egice(1:nc3)=pp_egice(pq3map(1:nc3))
  dt_ftop(1:nc3)=dp_ftop(pq3map(1:nc3))
  dt_tb(1:nc3)=dp_tb(pq3map(1:nc3))
  dt_fb(1:nc3)=dp_fb(pq3map(1:nc3))
  dt_timelt(1:nc3)=dp_timelt(pq3map(1:nc3))
  dt_salflxf(1:nc3)=dp_salflxf(pq3map(1:nc3))
  dt_salflxs(1:nc3)=dp_salflxs(pq3map(1:nc3))
  dt_nk(1:nc3)=dp_nk(pq3map(1:nc3))
  dt_wavail(1:nc3)=dp_wavail(pq3map(1:nc3))
  call icetempn1(nc3,dt,it_tn0(1:nc3),it_tn1(1:nc3),it_tn2(1:nc3),it_dic(1:nc3),it_dsn(1:nc3),it_fracice(1:nc3), &
                 it_tsurf(1:nc3),it_sto(1:nc3),dt_ftop(1:nc3),dt_tb(1:nc3),dt_fb(1:nc3),dt_timelt(1:nc3),        &
                 dt_salflxf(1:nc3),dt_salflxs(1:nc3),dt_nk(1:nc3),dt_wavail(1:nc3),pt_egice(1:nc3),diag)
  ip_tn0(pq3map(1:nc3))=it_tn0(1:nc3)
  ip_tn1(pq3map(1:nc3))=it_tn1(1:nc3)
  ip_tn2(pq3map(1:nc3))=it_tn2(1:nc3)
  ip_dic(pq3map(1:nc3))=it_dic(1:nc3)
  ip_dsn(pq3map(1:nc3))=it_dsn(1:nc3)
  ip_fracice(pq3map(1:nc3))=it_fracice(1:nc3)
  ip_tsurf(pq3map(1:nc3))=it_tsurf(1:nc3)
  ip_sto(pq3map(1:nc3))=it_sto(1:nc3)
  dp_salflxf(pq3map(1:nc3))=dt_salflxf(1:nc3)
  dp_salflxs(pq3map(1:nc3))=dt_salflxs(1:nc3)
  dp_nk(pq3map(1:nc3))=dt_nk(1:nc3)
end if

! the following are snow to ice processes
where (ip_dic>icemin)
  xxx=ip_dic+ip_dsn-(rhosn*ip_dsn+rhoic*ip_dic)/rhowt ! white ice formation
  excess=max(ip_dsn-xxx,0.)*rhosn/rhowt               ! white ice formation
elsewhere
  xxx=0.
  excess=0.
end where
where (excess>0..and.dp_nk==2)
  ip_tn1=(0.5*ip_dic*ip_tn1+rhowt/rhoic*excess*ip_tn0*cps/cpi)/(0.5*ip_dic+rhowt/rhoic*excess) ! Assume 2 levels of ice
elsewhere (excess>0..and.dp_nk==1)
  ip_tn1=(ip_dic*ip_tn1+rhowt/rhoic*excess*ip_tn0*cps/cpi)/(ip_dic+rhowt/rhoic*excess)         ! Assume 1 level of ice
elsewhere (excess>0.)
  ip_tsurf=ip_tsurf*(max(ip_dic,icemin)+rhowt/rhoic*excess*gamms/gammi)/(max(ip_dic,icemin)+rhowt/rhoic*excess)
end where
excess=excess+max(ip_dsn-0.2,0.)*rhosn/rhowt  ! Snow depth limitation and conversion to ice
ip_dsn=min(min(xxx,ip_dsn),0.2)
ip_dic=ip_dic+excess*rhowt/rhoic
dp_salflxf=dp_salflxf-excess*rhowt/dt
dp_salflxs=dp_salflxs+excess*rhowt/dt

! Ice depth limitation for poor initial conditions
xxx=max(ip_dic-icemax,0.)
ip_dic=min(icemax,ip_dic)
dp_salflxs=dp_salflxs-rhoic*xxx/dt

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

subroutine icetemps2(nc,dt,it_tn0,it_tn1,it_tn2,it_dic,it_dsn,it_fracice,it_tsurf,it_sto,dt_ftop,dt_tb,dt_fb, &
                     dt_timelt,dt_salflxf,dt_salflxs,dt_nk,dt_wavail,pt_egice,diag)

implicit none

integer, intent(in) :: nc,diag
integer ii,iqi
real, intent(in) :: dt
real qneed,fact
real htup,htdown
real, dimension(nc) :: con,rhin,fs,conb,tnew,qmax
real, dimension(nc) :: subl,snmelt,dhb,isubl,ssubl,smax,sbrine
real, dimension(nc,0:2) :: fl
real, dimension(nc), intent(inout) :: it_tn0,it_tn1,it_tn2
real, dimension(nc), intent(inout) :: it_dic,it_dsn,it_fracice,it_tsurf,it_sto
real, dimension(nc), intent(inout) :: dt_ftop,dt_tb,dt_fb,dt_timelt,dt_salflxf,dt_salflxs,dt_wavail
integer, dimension(nc), intent(inout) :: dt_nk
real, dimension(nc), intent(in) :: pt_egice

! update surface temperature
rhin=real(dt_nk)/it_dic
conb=2./(it_dsn/condsnw+1./(condice*rhin))
! Ammendments if surface layer < 5cms
where (it_dsn<0.05)
  con=1./(it_dsn/condsnw+0.5/(condice*rhin))
  it_tn1=it_tn1+cps*it_dsn*(it_tn0-it_tn1)/(cpi*rhin)
  it_tn0=it_tn1
elsewhere
  con=2.*condsnw/max(it_dsn,icemin)
end where
tnew=it_tsurf+con*(it_tn0-it_tsurf)/(gamms/dt+con)
fs=con*(it_tn0-tnew)
it_tsurf=it_tsurf+dt*(dt_ftop+fs)/gamms

! fluxes between snow and ice layers
where (it_dsn<0.05)
  fl(:,0)=fs
elsewhere
  fl(:,0)=conb*(it_tn1-it_tn0)
  it_tn0=it_tn0+dt*(fl(:,0)-fs)/(it_dsn*cps)
end where

! Surface evap/sublimation (can be >0 or <0)
! Note : dt*eg/hl in Kgm/m**2 => mms of water
subl=dt*pt_egice*lf/(lv*qsnow)     ! net snow+ice sublimation
ssubl=min(subl,it_dsn)             ! snow only component of sublimation
isubl=(subl-ssubl)*rhosn/rhoic     ! ice only component of sublimation
dt_salflxf=dt_salflxf+isubl*rhosn/dt
dt_salflxs=dt_salflxs-isubl*rhosn/dt
it_dsn=max(it_dsn-ssubl,0.)
it_dic=max(it_dic-isubl,0.)


! Limit maximum energy stored in brine pockets
qmax=qice*0.5*max(it_dic-himin,0.)
smax=max(it_sto-qmax,0.)/qice
smax=min(smax,it_dic)
it_sto=max(it_sto-smax*qice,0.)
dt_salflxs=dt_salflxs-smax*rhoic/dt
it_dic=max(it_dic-smax,0.)

! Snow melt (snmelt >= 0)
snmelt=max(0.,it_tsurf-273.16)*gamms/qsnow
snmelt=min(snmelt,it_dsn)
it_tsurf=it_tsurf-snmelt*qsnow/gamms
dt_salflxf=dt_salflxf-snmelt*rhosn/dt        ! melt fresh water snow (no salt when melting snow)
it_dsn=max(it_dsn-snmelt,0.)

! Ice melt
where (it_dsn<icemin)
  snmelt=max(0.,it_tsurf-dt_timelt)*gammi/qice
  snmelt=min(snmelt,it_dic)
  it_tsurf=it_tsurf-snmelt*qice/gammi
  dt_salflxs=dt_salflxs-snmelt*rhoic/dt
  it_dic=max(it_dic-snmelt,0.)
end where

! fluxes between ice layers
where (dt_nk>=2)
  fl(:,1)=condice*(it_tn2-it_tn1)*rhin     ! Between ice layers 2 and 1
  fl(:,2)=condice*(dt_tb-it_tn2)*2.*rhin   ! Between water and ice layer 2
  fl(:,2)=max(min(fl(:,2),1000.),-1000.)
  dhb=dt*(fl(:,2)-dt_fb)/qice              ! Excess flux between water and ice layer 2
  dhb=max(dhb,-it_dic)
  dhb=min(dhb,dt_wavail*rhowt/rhoic)
  fl(:,2)=dhb*qice/dt+dt_fb
elsewhere
  fl(:,1)=condice*(dt_tb-it_tn1)*2.*rhin   ! Between water and ice layer 1
  fl(:,1)=max(min(fl(:,1),1000.),-1000.)
  fl(:,2)=0.
  dhb=dt*(fl(:,1)-dt_fb)/qice              ! Excess flux between water and ice layer 1
  dhb=max(dhb,-it_dic)
  dhb=min(dhb,dt_wavail*rhowt/rhoic)
  fl(:,1)=dhb*qice/dt+dt_fb
end where
dt_salflxs=dt_salflxs+dhb*rhoic/dt
it_dic=max(it_dic+dhb,0.)

do iqi=1,nc

  if (dt_nk(iqi)==2) then
    ! update temperature in top layer of ice
    it_tn1(iqi)=it_tn1(iqi)+dt*(fl(iqi,1)-fl(iqi,0))*rhin(iqi)/cpi
    it_tn2(iqi)=it_tn2(iqi)+dt/cpi*(fl(iqi,2)-fl(iqi,1))*rhin(iqi)
    ! use stored heat in brine pockets to keep temperature at -0.1 until heat is used up
    if (it_sto(iqi)>0..and.it_tn1(iqi)<dt_timelt(iqi)) then
      qneed=(dt_timelt(iqi)-it_tn1(iqi))*cpi/rhin(iqi) ! J/m**2
      if (it_sto(iqi)<qneed) then
        fact=it_sto(iqi)/qneed
        fact=min(max(fact,0.),1.)
        it_tn1(iqi)=dt_timelt(iqi)*fact+it_tn1(iqi)*(1.-fact)
        it_sto(iqi)=0.
      else
        it_tn1(iqi)=dt_timelt(iqi)
        it_sto(iqi)=it_sto(iqi)-qneed
      end if
    end if
  else
    ! remove brine store for single layer
    sbrine=min(it_sto/qice,it_dic)
    it_sto=max(it_sto-sbrine*qice,0.)
    dt_salflxs=dt_salflxs-sbrine*rhoic/dt
    it_dic=max(it_dic-sbrine,0.)
    it_tn1(iqi)=it_tn1(iqi)+dt/cpi*(fl(iqi,1)-fl(iqi,0))*rhin(iqi)
  end if

  ! test whether to change number of layers
  htup=real(dt_nk(iqi)+1)*himin
  htdown=real(dt_nk(iqi))*himin
  if (it_dic(iqi)>=htup.and.dt_nk(iqi)<2) then
    ! code to increase number of layers
    dt_nk(iqi)=2
    it_tn2(iqi)=it_tn1(iqi)
    ! snow depth has not changed, so split ice layer into two
  elseif (it_dic(iqi)<=htdown) then
    ! code to decrease number of layers
    dt_nk(iqi)=dt_nk(iqi)-1
    if (dt_nk(iqi)==1) then
      ! merge ice layers into one
      it_tn1(iqi)=0.5*(it_tn1(iqi)+it_tn2(iqi))
      it_tn2(iqi)=it_tn1(iqi)
    else
      it_tsurf(iqi)=(it_tsurf(iqi)+cpi*it_dic(iqi)*it_tn1(iqi)/gammi)/(1.+cpi*it_dic(iqi)/gammi)
      it_tsurf(iqi)=(it_tsurf(iqi)+cps*it_dsn(iqi)*it_tn0(iqi)/gammi)/(1.+cps*it_dsn(iqi)/gammi)
      it_tn1(iqi)=it_tsurf(iqi)
      it_tn2(iqi)=it_tsurf(iqi)
    end if
  end if
  
end do

! update unassigned levels
where (dt_nk==1)
  it_tn2=it_tn1
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

subroutine icetempn2(nc,dt,it_tn0,it_tn1,it_tn2,it_dic,it_dsn,it_fracice,it_tsurf,it_sto,dt_ftop,dt_tb,dt_fb, &
                     dt_timelt,dt_salflxf,dt_salflxs,dt_nk,dt_wavail,pt_egice,diag)

implicit none

integer, intent(in) :: nc,diag
integer ii,iqi
real, intent(in) :: dt
real qneed,fact
real htup,htdown
real, dimension(nc) :: con,rhin,qmax,ftopadj
real, dimension(nc) :: subl,simelt,sbrine,smax,dhb,tnew
real, dimension(nc,0:2) :: fl
real, dimension(nc), intent(inout) :: it_tn0,it_tn1,it_tn2
real, dimension(nc), intent(inout) :: it_dic,it_dsn,it_fracice,it_tsurf,it_sto
real, dimension(nc), intent(inout) :: dt_ftop,dt_tb,dt_fb,dt_timelt,dt_salflxf,dt_salflxs,dt_wavail
integer, dimension(nc), intent(inout) :: dt_nk
real, dimension(nc), intent(in) :: pt_egice

! Get total flux into thin surface ice layer by adding flux from below
rhin=real(dt_nk)/it_dic
con=2.*condice*rhin
tnew=it_tsurf+con*(it_tn1-it_tsurf)/(gammi/dt+con)
ftopadj=con*(it_tn1-tnew)
it_tsurf=it_tsurf+dt*(dt_ftop+ftopadj)/gammi

! Surface evap/sublimation (can be >0 or <0)
subl=dt*pt_egice*lf/(lv*qice)
dt_salflxf=dt_salflxf+subl*rhoic/dt
dt_salflxs=dt_salflxs-subl*rhoic/dt
it_dic=max(it_dic-subl,0.)

! Limit maximum energy stored in brine pockets
qmax=qice*0.5*max(it_dic-himin,0.)
smax=max(it_sto-qmax,0.)/qice
smax=min(smax,it_dic)
it_sto=max(it_sto-smax*qice,0.)
dt_salflxs=dt_salflxs-smax*rhoic/dt
it_dic=max(it_dic-smax,0.)

! Ice melt (simelt >= 0)
simelt=max(0.,it_tsurf-dt_timelt)*gammi/qice
simelt=min(simelt,it_dic)
it_tsurf=it_tsurf-simelt*qice/gammi
dt_salflxs=dt_salflxs-simelt*rhoic/dt
it_dic=max(it_dic-simelt,0.)

! Upward fluxes between various levels below surface
fl(:,0)=ftopadj                            ! Middle of first ice layer to thin surface layer
where (dt_nk>=2)
  fl(:,1)=condice*(it_tn2-it_tn1)*rhin     ! Between ice layers 2 and 1
  fl(:,2)=condice*(dt_tb-it_tn2)*2.*rhin   ! Between water and ice layer 2
  fl(:,2)=max(min(fl(:,2),1000.),-1000.)
  dhb=dt*(fl(:,2)-dt_fb)/qice              ! Excess flux between water and ice layer 2
  dhb=max(dhb,-it_dic)
  dhb=min(dhb,dt_wavail*rhowt/rhoic)
  fl(:,2)=dhb*qice/dt+dt_fb
elsewhere
  fl(:,1)=condice*(dt_tb-it_tn1)*2.*rhin   ! Between water and ice layer 1
  fl(:,1)=max(min(fl(:,1),1000.),-1000.)
  fl(:,2)=0.
  dhb=dt*(fl(:,1)-dt_fb)/qice              ! Excess flux between water and ice layer 1
  dhb=max(dhb,-it_dic)
  dhb=min(dhb,dt_wavail*rhowt/rhoic)
  fl(:,1)=dhb*qice/dt+dt_fb
end where
dt_salflxs=dt_salflxs+dhb*rhoic/dt
it_dic=max(it_dic+dhb,0.)

do iqi=1,nc

  if (dt_nk(iqi)==2) then
    ! update temperature and brine store for middle layer
    it_tn1(iqi)=it_tn1(iqi)+dt/cpi*(fl(iqi,1)-fl(iqi,0))*rhin(iqi)
    it_tn2(iqi)=it_tn2(iqi)+dt/cpi*(fl(iqi,2)-fl(iqi,1))*rhin(iqi)
    if (it_sto(iqi)>0..and.it_tn1(iqi)<dt_timelt(iqi)) then
      qneed=(dt_timelt(iqi)-it_tn1(iqi))*cpi/rhin(iqi) ! J/m**2
      if (it_sto(iqi)<qneed) then
        fact=it_sto(iqi)/qneed
        fact=min(max(fact,0.),1.)
        it_tn1(iqi)=dt_timelt(iqi)*fact+it_tn1(iqi)*(1.-fact)
        it_sto(iqi)=0.
      else
        it_tn1(iqi)=dt_timelt(iqi)
        it_sto(iqi)=it_sto(iqi)-qneed
      end if
    end if
  else
    ! remove brine store for single layer
    sbrine=min(it_sto/qice,it_dic)
    it_sto=max(it_sto-sbrine*qice,0.)
    dt_salflxs=dt_salflxs-sbrine*rhoic/dt
    it_dic=max(it_dic-sbrine,0.)
    it_tn1(iqi)=it_tn1(iqi)+dt/cpi*(fl(iqi,1)-fl(iqi,0))*rhin(iqi)
  end if

  ! test whether to change number of layers
  htup=real(dt_nk(iqi)+1)*himin
  htdown=real(dt_nk(iqi))*himin
  if (it_dic(iqi)>=htup.and.dt_nk(iqi)<2) then
    ! code to increase number of layers
    dt_nk(iqi)=2
    it_tn2(iqi)=it_tn1(iqi)
    ! split ice layer into two
  elseif (it_dic(iqi)<=htdown) then
    ! code to decrease number of layers
    dt_nk(iqi)=dt_nk(iqi)-1
    if (dt_nk(iqi)==1) then
      ! merge ice layers into one
      it_tn1(iqi)=0.5*(it_tn1(iqi)+it_tn2(iqi))
      it_tn2(iqi)=it_tn1(iqi)
    else
      it_tsurf(iqi)=(it_tsurf(iqi)+cpi*it_dic(iqi)*it_tn1(iqi)/gammi)/(1.+cpi*it_dic(iqi)/gammi)
      it_tn1(iqi)=it_tsurf(iqi)
      it_tn2(iqi)=it_tsurf(iqi)
    end if
  end if

end do

! update unassigned levels
it_tn0=it_tn1
where (dt_nk==1)
  it_tn2=it_tn1
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

subroutine icetemps1(nc,dt,it_tn0,it_tn1,it_tn2,it_dic,it_dsn,it_fracice,it_tsurf,it_sto,dt_ftop,dt_tb,dt_fb, &
                     dt_timelt,dt_salflxf,dt_salflxs,dt_nk,dt_wavail,pt_egice,diag)

implicit none

integer, intent(in) :: nc,diag
real, intent(in) :: dt
real, dimension(nc) :: con,f0,tnew
real, dimension(nc) :: subl,snmelt,smax,dhb,isubl,ssubl,ssnmelt,gamm
real, dimension(nc), intent(inout) :: it_tn0,it_tn1,it_tn2
real, dimension(nc), intent(inout) :: it_dic,it_dsn,it_fracice,it_tsurf,it_sto
real, dimension(nc), intent(inout) :: dt_ftop,dt_tb,dt_fb,dt_timelt,dt_salflxf,dt_salflxs,dt_wavail
integer, dimension(nc), intent(inout) :: dt_nk
real, dimension(nc), intent(in) :: pt_egice

! Update tsurf and ti based on fluxes from above and below
con=1./(it_dsn/condsnw+max(it_dic,icemin)/condice)
gamm=(gammi*max(it_dic,icemin)+gamms*it_dsn)/(max(it_dic,icemin)+it_dsn) ! for energy conservation
tnew=it_tsurf+con*(dt_tb-it_tsurf)/(gamm/dt+con)
f0=con*(dt_tb-tnew) ! flux from below
f0=min(max(f0,-1000.),1000.)
dhb=dt*(f0-dt_fb)/qice
dhb=max(dhb,-it_dic)
dhb=min(dhb,dt_wavail*rhowt/rhoic)
f0=dhb*qice/dt+dt_fb
it_tsurf=it_tsurf+dt*(dt_ftop+f0)/gamm

! Surface evap/sublimation (can be >0 or <0)
subl=dt*pt_egice*lf/(lv*qsnow)
ssubl=min(subl,it_dsn)             ! snow only component of sublimation
isubl=(subl-ssubl)*rhosn/rhoic     ! ice only component of sublimation
dt_salflxf=dt_salflxf+isubl*rhosn/dt
dt_salflxs=dt_salflxs-isubl*rhosn/dt
it_dsn=max(it_dsn-ssubl,0.)
it_dic=max(it_dic-isubl,0.)

! Limit maximum energy stored in brine pockets
smax=it_sto/qice
smax=min(smax,it_dic)
it_sto=max(it_sto-smax*qice,0.)
dt_salflxs=dt_salflxs-smax*rhoic/dt
it_dic=max(it_dic-smax,0.)

! Snow melt (snmelt >= 0)
snmelt=max(0.,it_tsurf-273.16)*gamms/qsnow
snmelt=min(snmelt,it_dsn)
it_tsurf=it_tsurf-snmelt*qsnow/gamms
dt_salflxf=dt_salflxf-snmelt*rhosn/dt        ! melt fresh water snow (no salt when melting snow)
it_dsn=max(it_dsn-snmelt,0.)

! Ice melt
where (it_dsn<icemin)
  snmelt=max(0.,it_tsurf-dt_timelt)*gammi/qice
  snmelt=min(snmelt,it_dic)
  it_tsurf=it_tsurf-snmelt*qice/gammi
  dt_salflxs=dt_salflxs-snmelt*rhoic/dt
  it_dic=max(it_dic-snmelt,0.)
end where

! Bottom ablation or accretion
dt_salflxs=dt_salflxs+dhb*rhoic/dt
it_dic=max(it_dic+dhb,0.)

! Recalculate thickness index
where (it_dic>=himin)
  dt_nk=1
end where

! Update unassigned levels
it_tn0=it_tsurf
it_tn1=it_tsurf
it_tn2=it_tsurf

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

subroutine icetempn1(nc,dt,it_tn0,it_tn1,it_tn2,it_dic,it_dsn,it_fracice,it_tsurf,it_sto,dt_ftop,dt_tb,dt_fb, &
                     dt_timelt,dt_salflxf,dt_salflxs,dt_nk,dt_wavail,pt_egice,diag)

implicit none

integer, intent(in) :: nc,diag
real, intent(in) :: dt
real, dimension(nc) :: con,subl,simelt,smax,tnew
real, dimension(nc) :: dhb,f0,f0a
real, dimension(nc), intent(inout) :: it_tn0,it_tn1,it_tn2
real, dimension(nc), intent(inout) :: it_dic,it_dsn,it_fracice,it_tsurf,it_sto
real, dimension(nc), intent(inout) :: dt_ftop,dt_tb,dt_fb,dt_timelt,dt_salflxf,dt_salflxs,dt_wavail
integer, dimension(nc), intent(inout) :: dt_nk
real, dimension(nc), intent(in) :: pt_egice

! Update tsurf and ti based on fluxes from above and below
con=condice/max(it_dic,icemin)
tnew=it_tsurf+con*(dt_tb-it_tsurf)/(gammi/dt+con)
f0=con*(dt_tb-tnew) ! flux from below
f0a=min(max(f0,-1000.),1000.)
dhb=dt*(f0a-dt_fb)/qice
dhb=max(dhb,-it_dic)
dhb=min(dhb,dt_wavail*rhowt/rhoic)
f0=dhb*qice/dt+dt_fb
it_tsurf=it_tsurf+dt*(dt_ftop+f0)/gammi

! Surface evap/sublimation (can be >0 or <0)
subl=dt*pt_egice*lf/(lv*qice)
dt_salflxf=dt_salflxf+subl*rhoic/dt
dt_salflxs=dt_salflxs-subl*rhoic/dt
it_dic=max(it_dic-subl,0.)

! Limit maximum energy stored in brine pockets
smax=it_sto/qice
smax=min(smax,it_dic)
it_sto=max(it_sto-smax*qice,0.)
dt_salflxs=dt_salflxs-smax*rhoic/dt
it_dic=max(it_dic-smax,0.)

! Ice melt (simelt >= 0)
simelt=max(0.,it_tsurf-dt_timelt)*gammi/qice
simelt=min(simelt,it_dic)
it_tsurf=it_tsurf-simelt*qice/gammi
dt_salflxs=dt_salflxs-simelt*rhoic/dt
it_dic=max(it_dic-simelt,0.)

! Bottom ablation or accretion
dt_salflxs=dt_salflxs+dhb*rhoic/dt
it_dic=max(it_dic+dhb,0.)

! Recalculate thickness index
where (it_dic>=himin)
  dt_nk=1
end where

! Update unassigned levels
it_tn0=it_tsurf
it_tn1=it_tsurf
it_tn2=it_tsurf

return
end subroutine icetempn1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Determine ice fluxes

subroutine iceflux(dt,a_sg,a_rg,a_rnd,a_snd,a_vnratio,a_fbvis,a_fbnir,a_u,a_v,a_temp,a_qg,a_ps,a_zmin,a_zmins, &
                   d_rho,d_ftop,d_tb,d_fb,d_timelt,d_nk,d_tauxicw,d_tauyicw,d_zcr,d_ndsn,d_ndic,d_nsto,d_oldu, &
                   d_oldv,diag)

implicit none

integer, intent(in) :: diag
integer ll,iqw
real, dimension(wfull), intent(in) :: a_sg,a_rg,a_rnd,a_snd,a_vnratio,a_fbvis,a_fbnir,a_u,a_v,a_temp,a_qg,a_ps,a_zmin,a_zmins
real, dimension(wfull,wlev), intent(in) :: d_rho
real, dimension(wfull), intent(inout) :: d_ftop,d_tb,d_fb,d_timelt,d_tauxicw,d_tauyicw,d_zcr,d_ndsn,d_ndic
real, dimension(wfull), intent(inout) :: d_nsto,d_oldu,d_oldv
integer, dimension(wfull), intent(inout) :: d_nk
real, intent(in) :: dt
real, dimension(wfull) :: qsat,dqdt,ri,vmag,rho,srcp,tnew,qsatnew,gamm,bot
real, dimension(wfull) :: fm,fh,fq,af,aft,afq
real, dimension(wfull) :: den,sig,root
real, dimension(wfull) :: alb,qmax,eye
real, dimension(wfull) :: uu,vv,du,dv,vmagn,icemagn
real, dimension(wfull) :: x,ustar,icemag,g,h,dgu,dgv,dhu,dhv,det
real, dimension(wfull) :: newiu,newiv,imass,dtsurf
real factch

! Prevent unrealistic fluxes due to poor input surface temperature
dtsurf=min(i_tsurf,273.16)
uu=a_u-i_u
vv=a_v-i_v
vmag=sqrt(uu*uu+vv*vv)
vmagn=max(vmag,0.01)
sig=exp(-grav*a_zmins/(rdry*a_temp))
srcp=sig**(rdry/cp)
rho=a_ps/(rdry*dtsurf)

p_zoice=0.001
p_zohice=0.001
p_zoqice=0.001
af=vkar**2/(log(a_zmin/p_zoice)*log(a_zmin/p_zoice))
!factch=1.       ! following CSIRO9
factch=sqrt(7.4) ! following CCAM sflux
aft=af/(factch*factch)
afq=aft


call getqsat(qsat,dqdt,dtsurf,a_ps)
ri=min(grav*(a_zmin**2/a_zmins)*(1.-dtsurf*srcp/a_temp)/vmagn**2,rimax)

where (ri>=0.)
  fm=1./(1.+bprm*ri)**2  ! no zo contrib for stable
  fh=fm
  fq=fh
elsewhere        ! ri is -ve
  root=sqrt(-ri*a_zmin/p_zoice)
  den=1.+cms*2.*bprm*af*root
  fm=1.-2.*bprm*ri/den
  root=sqrt(-ri*a_zmins/p_zoice)
  den=1.+chs*2.*bprm*factch*aft*root
  fh=1.-2.*bprm*ri/den
  den=1.+chs*2.*bprm*factch*afq*root
  fq=1.-2.*bprm*ri/den
end where

! egice is for evaporating (lv).  Melting is included above with lf.

p_wetfacice=max(1.+.008*min(dtsurf-273.16,0.),0.)
p_cdice=af*fm
p_cdsice=aft*fh
p_cdqice=afq*fq
p_fgice=rho*p_cdsice*cp*vmag*(dtsurf-a_temp/srcp)
p_egice=p_wetfacice*rho*p_cdqice*lv*vmag*(qsat-a_qg)
p_egice=min(p_egice,d_ndic*qice*lv/(lf*dt))

! index of different ice thickness configurations
d_nk=min(int(d_ndic/himin),2)

! radiation
alb=     a_vnratio*(p_icevisdiralb*a_fbvis+p_icevisdifalb*(1.-a_fbvis))+ &
    (1.-a_vnratio)*(p_icevisdifalb*a_fbvis+p_icevisdifalb*(1.-a_fbvis))
qmax=qice*0.5*max(d_ndic-himin,0.)
where (d_ndsn<icemin.and.d_nsto<qmax.and.d_nk>0)
  eye=0.35
elsewhere
  eye=0.
end where
d_nsto=d_nsto+dt*a_sg*(1.-alb)*eye

! Explicit estimate of fluxes
d_ftop=-p_fgice-p_egice*ls/lv+a_rg-emisice*sbconst*dtsurf**4+a_sg*(1.-alb)*(1.-eye) ! first guess
bot=4.*emisice*sbconst*dtsurf**3+rho*vmag*(p_cdsice*cp+p_cdqice*p_wetfacice*lv*dqdt)
where (d_ndsn<=icemin)
  gamm=gammi
elsewhere (d_nk>0)
  gamm=gamms
elsewhere
  gamm=(gammi*max(d_ndic,icemin)+gamms*d_ndsn)/(max(d_ndic,icemin)+d_ndsn)
end where

! water temperature at bottom of ice
d_tb=w_temp(:,1)

! iterative method to estimate ice velocity after stress tensors are applied
imass=max(rhoic*i_dic+rhosn*i_dsn,10.) ! ice mass per unit area
newiu=i_u
newiv=i_v
do iqw=1,wfull
  do ll=1,20 ! max iterations
    uu(iqw)=a_u(iqw)-newiu(iqw)
    vv(iqw)=a_v(iqw)-newiv(iqw)
    du(iqw)=fluxwgt*w_u(iqw,1)+(1.-fluxwgt)*d_oldu(iqw)-newiu(iqw)
    dv(iqw)=fluxwgt*w_v(iqw,1)+(1.-fluxwgt)*d_oldv(iqw)-newiv(iqw)
  
    vmagn(iqw)=max(sqrt(uu(iqw)*uu(iqw)+vv(iqw)*vv(iqw)),0.01)
    icemagn(iqw)=max(sqrt(du(iqw)*du(iqw)+dv(iqw)*dv(iqw)),0.0001)

    g(iqw)=i_u(iqw)-newiu(iqw)+dt*(rho(iqw)*p_cdice(iqw)*vmagn(iqw)*uu(iqw)        &
                                +d_rho(iqw,1)*0.00536*icemagn(iqw)*du(iqw))/imass(iqw)
    h(iqw)=i_v(iqw)-newiv(iqw)+dt*(rho(iqw)*p_cdice(iqw)*vmagn(iqw)*vv(iqw)        &
                                +d_rho(iqw,1)*0.00536*icemagn(iqw)*dv(iqw))/imass(iqw)
  
    dgu(iqw)=-1.-dt*(rho(iqw)*p_cdice(iqw)*vmagn(iqw)*(1.+(uu(iqw)/vmagn(iqw))**2) &
               +d_rho(iqw,1)*0.00536*icemagn(iqw)*(1.+(du(iqw)/icemagn(iqw))**2))/imass(iqw)
    dhu(iqw)=-dt*(rho(iqw)*p_cdice(iqw)*uu(iqw)*vv(iqw)/vmagn(iqw)                 &
               +d_rho(iqw,1)*0.00536*du(iqw)*dv(iqw)/icemagn(iqw))/imass(iqw)
    dgv(iqw)=dhu(iqw)
    dhv(iqw)=-1.-dt*(rho(iqw)*p_cdice(iqw)*vmagn(iqw)*(1.+(vv(iqw)/vmagn(iqw))**2) &
               +d_rho(iqw,1)*0.00536*icemagn(iqw)*(1.+(dv(iqw)/icemagn(iqw))**2))/imass(iqw)

    det(iqw)=dgu(iqw)*dhv(iqw)-dgv(iqw)*dhu(iqw)

    newiu(iqw)=newiu(iqw)-0.9*( g(iqw)*dhv(iqw)-h(iqw)*dgv(iqw))/det(iqw)
    newiv(iqw)=newiv(iqw)-0.9*(-g(iqw)*dhu(iqw)+h(iqw)*dgu(iqw))/det(iqw)

    !newiu(iqw)=max(newiu(iqw),min(a_u(iqw),w_u(iqw,1),i_u(iqw)))
    !newiu(iqw)=min(newiu(iqw),max(a_u(iqw),w_u(iqw,1),i_u(iqw)))
    !newiv(iqw)=max(newiv(iqw),min(a_v(iqw),w_v(iqw,1),i_v(iqw)))
    !newiv(iqw)=min(newiv(iqw),max(a_v(iqw),w_v(iqw,1),i_v(iqw)))
  
    if (abs(g(iqw))<1.E-3.and.abs(h(iqw))<1.E-3) exit
  end do
end do

! momentum transfer
uu=a_u-newiu
vv=a_v-newiv
du=fluxwgt*w_u(:,1)+(1.-fluxwgt)*d_oldu-newiu
dv=fluxwgt*w_v(:,1)+(1.-fluxwgt)*d_oldv-newiv
vmagn=max(sqrt(uu*uu+vv*vv),0.01)
icemagn=max(sqrt(du*du+dv*dv),0.0001)
p_tauxica=rho*p_cdice*vmagn*uu
p_tauyica=rho*p_cdice*vmagn*vv
d_tauxicw=-d_rho(:,1)*0.00536*icemagn*du
d_tauyicw=-d_rho(:,1)*0.00536*icemagn*dv
!d_tauxicw=(i_u-newiu)*imass/dt+p_tauxica
!d_tauyicw=(i_v-newiv)*imass/dt+p_tauyica
ustar=sqrt(sqrt(d_tauxicw*d_tauxicw+d_tauyicw*d_tauyicw)/d_rho(:,1))
ustar=max(ustar,5.E-4)
d_fb=cp0*d_rho(:,1)*0.006*ustar*(d_tb-d_timelt)
d_fb=min(max(d_fb,-1000.),1000.)  

! Estimate fluxes to prevent overshoot
tnew=min(i_tsurf+d_ftop/(gamm/dt+bot),273.16)
call getqsat(qsatnew,dqdt,tnew,a_ps)
p_fgice=rho*p_cdsice*cp*vmag*(tnew-a_temp/srcp)
p_egice=p_wetfacice*rho*p_cdqice*lv*vmag*(qsatnew-a_qg)
p_egice=min(p_egice,d_ndic*qice*lv/(lf*dt))
d_ftop=-p_fgice-p_egice*ls/lv+a_rg-emisice*sbconst*tnew**4+a_sg*(1.-alb)*(1.-eye)

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
real, dimension(wfull) :: smixr,qsat,dqdt,atu,atv,dmag

! water
call getqsat(qsat,dqdt,w_temp(:,1),a_ps)
if (zomode==0) then
  smixr=0.98*qsat
else
  smixr=qsat
end if
atu=a_u-w_u(:,1)
atv=a_v-w_v(:,1)
call scrntile(tscrn,qgscrn,uscrn,u10,p_zo,p_zoh,p_zoq,w_temp(:,1),smixr,atu,atv,a_temp,a_qg,a_zmin,a_zmins,diag)
p_tscrn=tscrn
p_qgscrn=qgscrn
dmag=max(sqrt(atu*atu+atv*atv),0.01)
atu=(a_u-w_u(:,1))*uscrn/dmag+w_u(:,1)
atv=(a_v-w_v(:,1))*uscrn/dmag+w_v(:,1)
p_uscrn=sqrt(atu*atu+atv*atv)
atu=(a_u-w_u(:,1))*u10/dmag+w_u(:,1)
atv=(a_v-w_v(:,1))*u10/dmag+w_v(:,1)
p_u10=sqrt(atu*atu+atv*atv)

! ice
call getqsat(qsat,dqdt,i_tsurf,a_ps)
smixr=p_wetfacice*qsat+(1.-p_wetfacice)*min(qsat,a_qg)
atu=a_u-i_u
atv=a_v-i_v
call scrntile(tscrn,qgscrn,uscrn,u10,p_zoice,p_zohice,p_zoqice,i_tsurf,smixr,atu,atv,a_temp,a_qg,a_zmin,a_zmins,diag)
p_tscrn=(1.-i_fracice)*p_tscrn+i_fracice*tscrn
p_qgscrn=(1.-i_fracice)*p_qgscrn+i_fracice*qgscrn
dmag=max(sqrt(atu*atu+atv*atv),0.01)
atu=(a_u-i_u)*uscrn/dmag+i_u
atv=(a_v-i_v)*uscrn/dmag+i_v
p_uscrn=(1.-i_fracice)*p_uscrn+i_fracice*sqrt(atu*atu+atv*atv)
atu=(a_u-i_u)*u10/dmag+i_u
atv=(a_v-i_v)*u10/dmag+i_v
p_u10=(1.-i_fracice)*p_u10+i_fracice*sqrt(atu*atu+atv*atv)

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

umag=max(sqrt(a_u*a_u+a_v*a_v),0.01)
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
  z_on_l  = a_zmin*vkar*grav*thetavstar/(thetav*ustar*ustar)
  z_on_l  = min(z_on_l,10.)
  zs_on_l = z_on_l*a_zmins/a_zmin  
  z0_on_l = z_on_l*zo/a_zmin
  zt_on_l = z_on_l*zoh/a_zmin
  zq_on_l = z_on_l*zoq/a_zmin
  where (z_on_l<0.)
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
  thetavstar = vkar*(thetav-sthetav)/integralh
  ustar      = vkar*umag/integralm
end do
tstar = vkar*(a_temp-stemp)/integralh
qstar = vkar*(a_qg-smixr)/integralq
      
! estimate screen diagnostics
z0s_on_l  = z0*zs_on_l/a_zmins
z0_on_l   = z0*z_on_l/a_zmin
z10_on_l  = z10*z_on_l/a_zmin
z0s_on_l  = min(z0s_on_l,10.)
z0_on_l   = min(z0_on_l,10.)
z10_on_l  = min(z10_on_l,10.)
neutrals  = log(a_zmins/z0)
neutral   = log(a_zmin/z0)
neutral10 = log(a_zmin/z10)
where (z_on_l<0.)
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
