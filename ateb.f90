
! This code was originally based on the TEB scheme of Masson, Boundary-Layer Meteorology, 94, p357 (2000)
! The snow scheme is based on Douville, Royer and Mahfouf, Climate Dynamics, 12, p21 (1995)
! The in-canyon vegetation is based on Kowalczyk et al, DAR Tech Paper 32 (1994), but simplified by assiming sigmaf=1.

! The main changes include an alternative formulation for in-canyon aerodynamical resistances based on Harman, et al (2004)
! and Kanada et al (2007), combined with a second canyon wall for completeness.  The scheme includes nrefl order reflections
! in the canyon for both longwave and shortwave radiation (in TEB infinite reflections are used for shortwave and 1st order
! reflections in longwave). A big-leaf vegetation tile is included in the canyon using the Kowalczyk et al (1994) scheme but
! with a simplified soil moisture budget and no modelling of the soil temperature since sigmaf=1.  Snow is also included in
! the canyon and on roofs using a single-layer scheme based on Douville, et al (1995).  Time dependent traffic heat fluxes
! are based on Coutts, et al (2007).

! Usual practice is:
!   call atebinit     ! to initalise state arrays, etc (use tebdisable to disable calls to ateb subroutines)
!   call atebloadm    ! to load previous state arrays (from tebsavem)
!   call atebtype     ! to define urban type (or use tebfndef to define urban properties at each grid point)
!   ...
!   do t=1,tmax
!     ...
!     call atebnewangle ! store current solar zenith and azimuthal angle (use tebccangle for CCAM or
!                       ! use tebnewangle1 for a single grid point)
!     call atebcalc     ! calculates urban temperatures, fluxes, etc and blends with input
!     call atebalb      ! blends input and urban albedo (use tebalb1 for a single grid point)
!     call atebcd       ! blends input and urban drag coefficent (use tebzo for blending momentum and heat roughness lengths)
!     ...
!   end do
!   ...
!   call atebsavem    ! to save current state arrays (for use by tebloadm)
!   call atebend      ! to deallocate memory before quiting

! only atebinit and atebcalc are manditory.  All other subroutine calls are optional.

! URBAN TYPES:
 
! 1 = Urban               (TAPM 31)
! 2 = Urban (low)         (TAPM 32)
! 3 = Urban (medium)      (TAPM 33)
! 4 = Urban (high)        (TAPM 34)
! 5 = Urban (cbd)         (TAPM 35)
! 6 = Industrial (low)    (TAPM 36)
! 7 = Industrial (medium) (TAPM 37)
! 8 = Industrial (high)   (TAPM 38)

module ateb

implicit none

private
public atebinit,atebcalc,atebend,atebzo,atebload,atebsave,atebtype,atebfndef,atebalb1, &
       atebnewangle1,atebccangle,atebdisable,atebloadm,atebsavem,atebcd,vegmode, &
       atebdwn,atebscrnout

! type definitions
type tatm
  real :: sg,rg,rho,temp,mixr,ps,pa,umag,udir,rnd,snd
end type tatm
type tout
  real :: fg,eg,ts,wf
end type tout
type trad
  real :: roof,road,walle,wallw,veg,rfsn,rdsn
end type trad
type tdiag
  real :: roofdelta,roaddelta,vegdelta,rfsndelta,rdsndelta
  real :: tempc,mixrc,tempr,mixrr,sigd,sigr,rfdzmin
  real :: accool,canyonrgout,roofrgout,tran,evap,c1
  real :: totdepth,netemiss,netrad,topu
  real :: cwa,cwe,cww,cwr,cra,crr,crw
  real :: canyontemp,canyonmix
end type tdiag
type twall
  real, dimension(3) :: temp
end type twall
type tsurf
  real, dimension(3) :: temp
  real :: water,snow,den,alpha
end type tsurf
type tvege
  real :: water,moist
end type tvege
type tdata
  real, dimension(3) :: roofdepth,walldepth,roaddepth,vegdepth
  real, dimension(3) :: roofcp,wallcp,roadcp,vegcp
  real, dimension(3) :: rooflambda,walllambda,roadlambda,veglambda
  real :: hwratio,sigmabld,industryfg,trafficfg,bldheight,vangle,hangle
  real :: roofalpha,wallalpha,roadalpha,ctime,roofemiss,wallemiss,roademiss
  real :: bldtemp,sigmaveg,vegalpha,vegemiss
end type tdata
type tprog
  real :: lzom,lzoh,cndzmin,cduv,vegtemp
  real :: tscrn,qscrn,uscrn,u10
end type tprog

! state arrays
integer, save :: ufull,iqut
integer, dimension(:), allocatable, save :: ugrid,mgrid
real, dimension(:), allocatable, save :: sigmau
real, dimension(:,:), allocatable, save :: atebdwn ! These variables are for CCAM onthefly.f
type(tsurf), dimension(:), allocatable, save :: roof,road
type(twall), dimension(:), allocatable, save :: walle,wallw
type(tvege), dimension(:), allocatable, save :: veg
type(tdata), dimension(:), allocatable, save :: fn
type(tprog), dimension(:), allocatable, save :: pg
! model parameters
integer, parameter :: resmeth=1           ! Canyon sensible heat transfer (0=Masson, 1=Harman, 2=Kusaka)
integer, parameter :: zohmeth=2           ! Urban roughness length for heat (0=0.1*zom, 1=Kanda, 2=0.003*zom)
integer, parameter :: acmeth=1            ! AC heat pump into canyon (0=Off, 1=On)
integer, parameter :: nrefl=3             ! Number of canyon reflections (default=3)
integer, parameter :: vegmode=2           ! In-canyon vegetation mode (0=50%/50%, 1=100%/0%, 2=0%/100%, where out-canyon/in-canyon = X%/Y%)
integer, parameter :: scrnmeth=1          ! Screen diagnostic method (0=Slab, 1=Canopy)
integer, parameter :: iqt = 3432          ! Diagnostic point (in terms of host grid)
! sectant solver parameters
integer, parameter :: nfgits=4            ! Maximum number of iterations for calculating sensible heat flux (default=6)
integer, parameter :: negits=2            ! Maximum number of iterations for calculating latent heat flux (default=2)
real, parameter    :: tol=0.001           ! Sectant method tolarance for sensible heat flux
real, parameter    :: tolqg=1.E-8         ! Sectant method tolarance for latent heat flux
real, parameter    :: alpha=0.7           ! Weighting for determining the rate of convergence when calculating canyon temperatures
! physical parameters
real, parameter :: waterden=1000.         ! water density (kg m^-3)
real, parameter :: icelambda=2.22         ! conductance of ice (W m^-1 K^-1)
real, parameter :: aircp=1004.64          ! Heat capapcity of dry air (J kg^-1 K^-1)
real, parameter :: icecp=2100.            ! Heat capacity of ice (J kg^-1 K^-1)
real, parameter :: grav=9.80616           ! gravity (m s^-2)
real, parameter :: vkar=0.4               ! von Karman constant
real, parameter :: lv=2.501e6             ! Latent heat of vaporisation (J kg^-1)
real, parameter :: lf=3.337e5             ! Latent heat of fusion (J kg^-1)
real, parameter :: ls=lv+lf               ! Latent heat of sublimation (J kg^-1)
real, parameter :: pi=3.1415927           ! pi
real, parameter :: rd=287.04              ! Gas constant for dry air
real, parameter :: rv=461.5               ! Gas constant for water vapor
real, parameter :: sbconst=5.67e-8        ! Stefan-Boltzmann constant
! snow parameters
real, parameter :: zosnow=0.001           ! Roughness length for snow (m)
real, parameter :: snowemiss=1.           ! snow emissitivity
real, parameter :: maxsnowalpha=0.85      ! max snow albedo
real, parameter :: minsnowalpha=0.5       ! min snow albedo
real, parameter :: maxsnowden=300.        ! max snow density (kg m^-3)
real, parameter :: minsnowden=100.        ! min snow density (kg m^-3)
! generic urban parameters
real, parameter :: refheight=0.6          ! Displacement height as a fraction of building height (Kanda et al 2007)
real, parameter :: zomratio=0.1           ! Roughness length to building height ratio (default=10%)
real, parameter :: zocanyon=0.1           ! Roughness length of in-canyon surfaces (m)
real, parameter :: maxrfwater=1.          ! Maximum roof water (kg m^-2)
real, parameter :: maxrdwater=1.          ! Maximum road water (kg m^-2)
real, parameter :: maxrfsn=1.             ! Maximum roof snow (kg m^-2)
real, parameter :: maxrdsn=1.             ! Maximum road snow (kg m^-2)
! in-canyon vegetation parameters
real, parameter :: zoveg=0.3              ! Roughness length of in-canyon vegetation (m)
real, parameter :: vegrlai=1.2            ! In-canyon vegetation LAI
real, parameter :: vegrsmin=80./vegrlai   ! Unconstrained canopy stomatal resistance
real, parameter :: maxvgwater=0.1*vegrlai ! Maximum leaf water (kg m^-2)
real, parameter :: swilt=0.18             ! In-canyon soil wilting point
real, parameter :: sfc=0.26               ! In-canyon soil field capacity
real, parameter :: ssat=0.42              ! In-canyon soil saturation point
! stability function parameters
integer, parameter :: icmax=1             ! number of iterations for stability functions (default=5)
real, parameter    :: bprm=5.             ! 4.7 in rams
real, parameter    :: chs=2.6             ! 5.3 in rams
real, parameter    :: cms=5.              ! 7.4 in rams
real, parameter    :: fmroot=0.57735
real, parameter    :: rimax=(1./fmroot-1.)/bprm
real, parameter    :: a_1=1.
real, parameter    :: b_1=2./3.
real, parameter    :: c_1=5.
real, parameter    :: d_1=0.35


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine prepare the arrays used by the aTEB scheme
! This is a compulsory subroutine that must be called during
! model initalisation

subroutine atebinit(ifull,sigu,zmin,diag)

implicit none

integer, intent(in) :: ifull,diag
integer, dimension(ifull) :: utype
integer iqu,iq,ii
real, intent(in) :: zmin
real, dimension(ifull), intent(in) :: sigu

if (diag.ge.1) write(6,*) "Initialising aTEB"

ufull=count(sigu.gt.0.)
if (ufull.eq.0) return

allocate(ugrid(ufull),mgrid(ifull),fn(ufull),pg(ufull))
allocate(roof(ufull),road(ufull),walle(ufull),wallw(ufull))
allocate(sigmau(ufull),veg(ufull))

! define grid arrays
ugrid=0
mgrid=0
sigmau=0.
iqu=0
iqut=0
do iq=1,ifull
  if (sigu(iq).gt.0.) then
    iqu=iqu+1
    ugrid(iqu)=iq
    mgrid(iq)=iqu
    sigmau(iqu)=sigu(iq)
    if (iq.ge.iqt.and.iqut.eq.0) iqut=iqu    
  end if
end do

if (iqut.eq.0) then
  !write(6,*) "WARN: Cannot located aTEB diagnostic point.  iqut=1"
  iqut=1
end if

! Initialise state variables
do ii=1,3
  roof%temp(ii)=291.
  road%temp(ii)=291.
  walle%temp(ii)=291.
  wallw%temp(ii)=291.
end do
roof%water=0.
roof%snow=0.
roof%den=minsnowden
roof%alpha=maxsnowalpha
road%water=0.
road%snow=0.
road%den=minsnowden
road%alpha=maxsnowalpha
veg%moist=0.5*(ssat+swilt)
veg%water=0.

do ii=1,3
  fn%roofdepth(ii)=0.1
  fn%walldepth(ii)=0.1
  fn%roaddepth(ii)=0.1
  fn%roofcp(ii)=2.E6
  fn%wallcp(ii)=2.E6
  fn%roadcp(ii)=2.E6
  fn%rooflambda(ii)=2.
  fn%walllambda(ii)=2.
  fn%roadlambda(ii)=2.
end do
fn%hwratio=1.
fn%sigmabld=0.5
fn%sigmaveg=0.5
fn%industryfg=0.
fn%trafficfg=0.
fn%bldheight=10.
fn%roofalpha=0.2
fn%wallalpha=0.2
fn%roadalpha=0.2
fn%vegalpha=0.2
fn%roofemiss=0.97
fn%wallemiss=0.97
fn%roademiss=0.97
fn%vegemiss=0.97
fn%bldtemp=291.
fn%vangle=0.
fn%hangle=0.
fn%ctime=0.

utype=1 ! default urban
call atebtype(ifull,utype,diag)

pg%cndzmin=max(zmin,0.1*fn%bldheight+1.)   ! updated in atebcalc
pg%lzom=log(pg%cndzmin/(0.1*fn%bldheight)) ! updated in atebcalc
pg%lzoh=6.+pg%lzom ! (Kanda et al 2005)    ! updated in atebcalc
pg%cduv=(vkar/pg%lzom)**2                  ! updated in atebcalc
pg%vegtemp=291.                            ! updated in atebcalc
pg%tscrn=291.                              ! updated in atebcalc
pg%qscrn=0.                                ! updated in atebcalc
pg%uscrn=0.                                ! updated in atebcalc
pg%u10=0.                                  ! updated in atebcalc

return
end subroutine atebinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine deallocates arrays used by the TEB scheme

subroutine atebend(diag)

implicit none

integer, intent(in) :: diag

if (ufull.eq.0) return
if (diag.ge.1) write(6,*) "Deallocating aTEB arrays"

deallocate(ugrid,mgrid,fn,pg)
deallocate(roof,road,walle,wallw)
deallocate(sigmau,veg)

return
end subroutine atebend

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine loads aTEB state arrays (not compulsory)

subroutine atebload(ifull,urban,diag)

implicit none

integer, intent(in) :: ifull,diag
integer ii
real, dimension(ifull,22), intent(in) :: urban

if (ufull.eq.0) return
if (diag.ge.1) write(6,*) "Load aTEB state arrays"

do ii=1,3
  roof%temp(ii)=urban(ugrid,ii)
  walle%temp(ii)=urban(ugrid,ii+3)
  wallw%temp(ii)=urban(ugrid,ii+6)
  road%temp(ii)=urban(ugrid,ii+9)
end do
veg%moist=urban(ugrid,13)
roof%water=urban(ugrid,14)
road%water=urban(ugrid,15)
veg%water=urban(ugrid,16)
roof%snow=urban(ugrid,17)
road%snow=urban(ugrid,18)
roof%den=urban(ugrid,19)
road%den=urban(ugrid,20)
roof%alpha=urban(ugrid,21)
road%alpha=urban(ugrid,22)

return
end subroutine atebload

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! temperature only version of tebsave

subroutine atebloadm(ifull,urban,moist,diag)

implicit none

integer, intent(in) :: ifull,diag
integer ii
real, dimension(ifull,12), intent(in) :: urban
real, dimension(ifull), intent(in) :: moist

if (ufull.eq.0) return
if (diag.ge.1) write(6,*) "Load aTEB state arrays"

do ii=1,3
  roof%temp(ii)=urban(ugrid,ii)
  walle%temp(ii)=urban(ugrid,ii+3)
  wallw%temp(ii)=urban(ugrid,ii+6)
  road%temp(ii)=urban(ugrid,ii+9)
end do
veg%moist=moist(ugrid)

return
end subroutine atebloadm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine loads aTEB type arrays (not compulsory)

subroutine atebtype(ifull,itype,diag)

implicit none

integer, intent(in) :: ifull,diag
integer ii
integer, dimension(ifull), intent(in) :: itype
integer, parameter :: maxtype = 8
real, dimension(ufull) :: tsigveg,tsigmabld
! In-canyon vegetation fraction
real, dimension(maxtype), parameter ::     csigveg=(/ 0.30, 0.45, 0.38, 0.34, 0.05, 0.50, 0.40, 0.30 /)
! Area fraction occupied by buildings
real, dimension(maxtype), parameter ::   csigmabld=(/ 0.40, 0.40, 0.44, 0.46, 0.65, 0.35, 0.40, 0.50 /)
! Building height (m)
real, dimension(maxtype), parameter ::  cbldheight=(/   6.,   4.,   6.,   8.,  20.,   5.,  10.,  15. /)
! Building height to width ratio
real, dimension(maxtype), parameter ::    chwratio=(/  0.6,  0.4,  0.6,  0.8,   2.,  0.5,   1.,  1.5 /)
! Industral sensible heat flux (W m^-2)
real, dimension(maxtype), parameter :: cindustryfg=(/  12.,   8.,  12.,  16.,  28.,  20.,  40.,  60. /)
! Daily averaged traffic sensible heat flux (W m^-2)
real, dimension(maxtype), parameter ::  ctrafficfg=(/   3.,   2.,   3.,   4.,   7.,   5.,  10.,  15. /)
! Comfort temperature (K)
real, dimension(maxtype), parameter :: cbldtemp=(/ 291.16, 291.16, 291.16, 291.16, 291.16, 291.16, 291.16, 291.16 /)
! Roof albedo
real, dimension(maxtype), parameter ::  croofalpha=(/ 0.20, 0.23, 0.20, 0.17, 0.13, 0.20, 0.20, 0.20 /)
! Wall albedo
real, dimension(maxtype), parameter ::  cwallalpha=(/ 0.33, 0.37, 0.33, 0.29, 0.22, 0.33, 0.33, 0.33 /)
! Road albedo
real, dimension(maxtype), parameter ::  croadalpha=(/ 0.10, 0.11, 0.10, 0.09, 0.07, 0.10, 0.10, 0.10 /)
! Veg albedo
real, dimension(maxtype), parameter ::   cvegalpha=(/ 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20 /)
! Roof emissitivity
real, dimension(maxtype), parameter ::  croofemiss=(/ 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90 /)
! Wall emissitivity
real, dimension(maxtype), parameter ::  cwallemiss=(/ 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85 /) 
! Road emissitivity
real, dimension(maxtype), parameter ::  croademiss=(/ 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94 /)
! Veg emissitivity
real, dimension(maxtype), parameter ::   cvegemiss=(/ 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96 /)
! Roof depths (m)
real, dimension(maxtype,3), parameter :: croofdepth=reshape((/ 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, &
                                                               0.05, 0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 0.40, &
                                                               0.40, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10 /), &
                                                               (/maxtype,3/))
! Wall depths (m)
real, dimension(maxtype,3), parameter :: cwalldepth=reshape((/ 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, &
                                                               0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, &
                                                               0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050 /), &
                                                               (/maxtype,3/))
! Road depths (m)
real, dimension(maxtype,3), parameter :: croaddepth=reshape((/ 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, &
                                                               0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, &
                                                               1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00 /), &
                                                               (/maxtype,3/))
! Roof heat capacity (J m^-3 K^-1)
real, dimension(maxtype,3), parameter :: croofcp=reshape((/ 2.11E6, 2.11E6, 2.11E6, 2.11E6, 2.11E6, 2.11E6, 2.11E6, 2.11E6, &
                                                            2.11E6, 0.28E6, 0.28E6, 0.28E6, 0.28E6, 0.28E6, 0.28E6, 0.28E6, &
                                                            0.28E6, 0.29E6, 0.29E6, 0.29E6, 0.29E6, 0.29E6, 0.29E6, 0.29E6 /), &
                                                            (/maxtype,3/))
! Wall heat capacity (J m^-3 K^-1)
real, dimension(maxtype,3), parameter :: cwallcp=reshape((/ 1.55E6, 1.55E6, 1.55E6, 1.55E6, 1.55E6, 1.55E6, 1.55E6, 1.55E6, &
                                                            1.55E6, 1.55E6, 1.55E6, 1.55E6, 1.55E6, 1.55E6, 1.55E6, 1.55E6, &
                                                            0.29E6, 0.29E6, 0.29E6, 0.29E6, 0.29E6, 0.29E6, 0.29E6, 0.29E6 /), &
                                                            (/maxtype,3/))
! Road heat capacity (J m^-3 K^-1)
real, dimension(maxtype,3), parameter :: croadcp=reshape((/ 1.94E6, 1.94E6, 1.94E6, 1.94E6, 1.94E6, 1.94E6, 1.94E6, 1.94E6, &
                                                            1.28E6, 1.28E6, 1.28E6, 1.28E6, 1.28E6, 1.28E6, 1.28E6, 1.28E6, &
                                                            1.28E6, 1.28E6, 1.28E6, 1.28E6, 1.28E6, 1.28E6, 1.28E6, 1.28E6 /), &
                                                            (/maxtype,3/))
! Roof conductance (W m^-1 K^-1)
real, dimension(maxtype,3), parameter :: crooflambda=reshape((/ 1.51, 1.51, 1.51, 1.51, 1.51, 1.51, 1.51, 1.51, &
                                                                1.51, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, &
                                                                0.08, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05 /), &
                                                                (/maxtype,3/))
! Wall conductance (W m^-1 K^-1)
real, dimension(maxtype,3), parameter :: cwalllambda=reshape((/ 0.9338, 0.9338, 0.9338, 0.9338, 0.9338, 0.9338, 0.9338, 0.9338, &
                                                                0.9338, 0.9338, 0.9338, 0.9338, 0.9338, 0.9338, 0.9338, 0.9338, &
                                                                0.0500, 0.0500, 0.0500, 0.0500, 0.0500, 0.0500, 0.0500, 0.0500  &
                                                                /), (/maxtype,3/))
! Road conductance (W m^-1 K^-1)
real, dimension(maxtype,3), parameter :: croadlambda=reshape((/ 0.7454, 0.7454, 0.7454, 0.7454, 0.7454, 0.7454, 0.7454, 0.7454, &
                                                                0.2513, 0.2513, 0.2513, 0.2513, 0.2513, 0.2513, 0.2513, 0.2513, &
                                                                0.2513, 0.2513, 0.2513, 0.2513, 0.2513, 0.2513, 0.2513, 0.2513  &
                                                                /), (/maxtype,3/))
                                                                
if (ufull.eq.0) return
if (diag.ge.1) write(6,*) "Load aTEB building properties"

if (any(itype(ugrid).lt.1).or.any(itype(ugrid).gt.maxtype)) then
  write(6,*) "ERROR: Urban type is out of range"
  stop
end if

select case(vegmode)
  case(0)
    tsigveg=0.5*csigveg(itype(ugrid))
    tsigmabld=csigmabld(itype(ugrid))*(1.-0.5*csigveg(itype(ugrid)))/(1.-csigveg(itype(ugrid)))
  case(1)
    tsigveg=0.
    tsigmabld=csigmabld(itype(ugrid))/(1.-csigveg(itype(ugrid)))
  case(2)
    tsigveg=csigveg(itype(ugrid))
    tsigmabld=csigmabld(itype(ugrid))
  case DEFAULT
    write(6,*) "ERROR: Unsupported vegmode ",vegmode
    stop
end select
fn%sigmaveg=tsigveg/(1.-tsigmabld)
fn%sigmabld=tsigmabld
fn%hwratio=chwratio(itype(ugrid))
fn%industryfg=cindustryfg(itype(ugrid))
fn%trafficfg=ctrafficfg(itype(ugrid))
fn%bldheight=cbldheight(itype(ugrid))
fn%roofalpha=croofalpha(itype(ugrid))
fn%wallalpha=cwallalpha(itype(ugrid))
fn%roadalpha=croadalpha(itype(ugrid))
fn%vegalpha=cvegalpha(itype(ugrid))
fn%roofemiss=croofemiss(itype(ugrid))
fn%wallemiss=cwallemiss(itype(ugrid))
fn%roademiss=croademiss(itype(ugrid))
fn%vegemiss=cvegemiss(itype(ugrid))
fn%bldtemp=cbldtemp(itype(ugrid))
do ii=1,3
  fn%roofdepth(ii)=croofdepth(itype(ugrid),ii)
  fn%walldepth(ii)=cwalldepth(itype(ugrid),ii)
  fn%roaddepth(ii)=croaddepth(itype(ugrid),ii)
  fn%roofcp(ii)=croofcp(itype(ugrid),ii)
  fn%wallcp(ii)=cwallcp(itype(ugrid),ii)
  fn%roadcp(ii)=croadcp(itype(ugrid),ii)
  fn%rooflambda(ii)=crooflambda(itype(ugrid),ii)
  fn%walllambda(ii)=cwalllambda(itype(ugrid),ii)
  fn%roadlambda(ii)=croadlambda(itype(ugrid),ii)
end do

return
end subroutine atebtype

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine specifies the urban properties for each grid point
!

subroutine atebfndef(ifull,ifn,diag)

implicit none

integer, intent(in) :: ifull,diag
integer ii
real, dimension(ifull,45), intent(in) :: ifn

if (ufull.eq.0) return
if (diag.ge.1) write(6,*) "Load aTEB building properties"

fn%hwratio=ifn(ugrid,1)
fn%sigmabld=ifn(ugrid,2)
fn%sigmaveg=ifn(ugrid,3)/(1.-ifn(ugrid,2))
fn%industryfg=ifn(ugrid,4)
fn%trafficfg=ifn(ugrid,5)
fn%bldheight=ifn(ugrid,6)
fn%roofalpha=ifn(ugrid,7)
fn%wallalpha=ifn(ugrid,8)
fn%roadalpha=ifn(ugrid,9)
fn%vegalpha=ifn(ugrid,10)
fn%roofemiss=ifn(ugrid,11)
fn%wallemiss=ifn(ugrid,12)
fn%roademiss=ifn(ugrid,13)
fn%vegemiss=ifn(ugrid,14)
fn%bldtemp=ifn(ugrid,15)
do ii=1,3
  fn%roofdepth(ii)=ifn(ugrid,15+ii)
  fn%walldepth(ii)=ifn(ugrid,18+ii)
  fn%roaddepth(ii)=ifn(ugrid,21+ii)
  fn%roofcp(ii)=ifn(ugrid,24+ii)
  fn%wallcp(ii)=ifn(ugrid,30+ii)
  fn%roadcp(ii)=ifn(ugrid,33+ii)
  fn%rooflambda(ii)=ifn(ugrid,36+ii)
  fn%walllambda(ii)=ifn(ugrid,39+ii)
  fn%roadlambda(ii)=ifn(ugrid,42+ii)
end do

return
end subroutine atebfndef

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine saves aTEB state arrays (not compulsory)

subroutine atebsave(ifull,urban,diag)

implicit none

integer, intent(in) :: ifull,diag
integer ii
real, dimension(ifull,22), intent(inout) :: urban

if (ufull.eq.0) return
if (diag.ge.1) write(6,*) "Save aTEB state arrays"

do ii=1,3
  urban(ugrid,ii)=roof%temp(ii)
  urban(ugrid,ii+3)=walle%temp(ii)
  urban(ugrid,ii+6)=wallw%temp(ii)
  urban(ugrid,ii+9)=road%temp(ii)
end do
urban(ugrid,13)=veg%moist
urban(ugrid,14)=roof%water
urban(ugrid,15)=road%water
urban(ugrid,16)=veg%water
urban(ugrid,17)=roof%snow
urban(ugrid,18)=road%snow
urban(ugrid,19)=roof%den
urban(ugrid,20)=road%den
urban(ugrid,21)=roof%alpha
urban(ugrid,22)=road%alpha

return
end subroutine atebsave

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! temperature only version of tebsave

subroutine atebsavem(ifull,urban,moist,diag)

implicit none

integer, intent(in) :: ifull,diag
integer ii
real, dimension(ifull,12), intent(inout) :: urban
real, dimension(ifull), intent(inout) :: moist

if (ufull.eq.0) return
if (diag.ge.1) write(6,*) "Save aTEB state arrays"

do ii=1,3
  urban(ugrid,ii)=roof%temp(ii)
  urban(ugrid,ii+3)=walle%temp(ii)
  urban(ugrid,ii+6)=wallw%temp(ii)
  urban(ugrid,ii+9)=road%temp(ii)
end do
moist(ugrid)=veg%moist

return
end subroutine atebsavem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine blends urban momentum and heat roughness lengths
! (This version neglects the displacement height (e.g., for CCAM))
!

subroutine atebzo(ifull,zom,zoh,diag,raw)

implicit none

integer, intent(in) :: ifull,diag
real, dimension(ifull), intent(inout) :: zom,zoh
real, dimension(ufull) :: workb,workc
real, parameter :: zr=1.e-15 ! limits minimum roughness length for heat
logical, intent(in), optional :: raw
logical mode

if (ufull.eq.0) return
if (diag.ge.1) write(6,*) "Calculate urban roughness lengths"

mode=.false.
if (present(raw)) mode=raw

if (mode) then
  zom(ugrid)=pg%cndzmin*exp(-pg%lzom)
  zoh(ugrid)=pg%cndzmin*exp(-pg%lzoh)
else 
  ! evaluate at canyon displacement height (really the atmospheric model should provide a displacement height)
  workb=sqrt((1.-sigmau)/log(pg%cndzmin/zom(ugrid))**2+sigmau/pg%lzom**2)
  workc=(1.-sigmau)/(log(pg%cndzmin/zom(ugrid))*log(pg%cndzmin/zoh(ugrid)))+sigmau/(pg%lzom*pg%lzoh)
  workc=workc/workb
  workb=pg%cndzmin*exp(-1./workb)
  workc=max(pg%cndzmin*exp(-1./workc),zr)
  zom(ugrid)=workb
  zoh(ugrid)=workc
  if (any(zoh(ugrid).le.zr)) write(6,*) "WARN: minimum zoh reached"
end if

return
end subroutine atebzo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine blends the urban drag coeff
!

subroutine atebcd(ifull,cduv,diag)

implicit none

integer, intent(in) :: ifull,diag
real, dimension(ifull), intent(inout) :: cduv

if (ufull.eq.0) return

cduv(ugrid)=(1.-sigmau)*cduv(ugrid)+sigmau*pg%cduv

return
end subroutine atebcd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates the urban contrabution to albedo.
! (selected grid points only)

subroutine atebalb1(is,ifull,alb,diag,raw)

implicit none

integer, intent(in) :: is,ifull,diag
integer, dimension(ifull) :: lgrid
integer i,ucount,ib,ie
real, dimension(ifull), intent(inout) :: alb
real, dimension(ifull) :: ualb
logical, intent(in), optional :: raw
logical mode

if (ufull.eq.0) return

mode=.false.
if (present(raw)) mode=raw

ucount=count(mgrid(is:is+ifull-1).ge.1)
if (ucount.eq.0) return

ucount=0
lgrid=0
do i=1,ifull
  if (mgrid(is+i-1).ge.1) then
    ucount=ucount+1
    lgrid(ucount)=i
  end if
end do

ib=mgrid(is+lgrid(1)-1)
ie=ib+ucount-1

call atebalbcalc(ib,ucount,ualb(1:ucount),diag)

if (mode) then
  alb(lgrid(1:ucount))=ualb(1:ucount)
else
  alb(lgrid(1:ucount))=(1.-sigmau(ib:ie))*alb(lgrid(1:ucount))+sigmau(ib:ie)*ualb(1:ucount)
end if

return
end subroutine atebalb1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Albedo calculations

subroutine atebalbcalc(is,ifull,alb,diag)

implicit none

integer, intent(in) :: is,ifull,diag
integer ie
real, dimension(ifull), intent(out) :: alb
real, dimension(ifull) :: albu,albr,snowdelta
real, dimension(ifull) :: wallpsi,roadpsi
type(trad), dimension(ifull) :: sg

ie=ifull+is-1

snowdelta=road(is:ie)%snow/(road(is:ie)%snow+maxrdsn)
call getswcoeff(ifull,sg,wallpsi,roadpsi,fn(is:ie),road(is:ie)%alpha,snowdelta)
albu=1.-(fn(is:ie)%hwratio*(sg%walle+sg%wallw)*(1.-fn(is:ie)%wallalpha)+snowdelta*sg%rdsn*(1.-road(is:ie)%alpha) &
    +(1.-snowdelta)*((1.-fn(is:ie)%sigmaveg)*sg%road*(1.-fn(is:ie)%roadalpha)+fn(is:ie)%sigmaveg*sg%veg*(1.-fn(is:ie)%vegalpha)))
snowdelta=roof(is:ie)%snow/(roof(is:ie)%snow+maxrfsn)
albr=sg%roof*(1.-snowdelta)*fn(is:ie)%roofalpha+sg%rfsn*snowdelta*roof(is:ie)%alpha
alb=fn(is:ie)%sigmabld*albr+(1.-fn(is:ie)%sigmabld)*albu

return
end subroutine atebalbcalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine stores the zenith angle and the solar azimuth angle
! (single grid point)

subroutine atebnewangle1(is,ifull,cosin,azimuthin,ctimein)

implicit none

integer, intent(in) :: is,ifull
real, dimension(ifull), intent(in) :: cosin     ! cosine of zenith angle
real, dimension(ifull), intent(in) :: azimuthin ! azimuthal angle
real, dimension(ifull), intent(in) :: ctimein   ! local hour (0<=ctime<=1)

if (ufull.eq.0) return

where (mgrid(is:ifull+is-1).ge.1)
  fn(mgrid(is:ifull+is-1))%hangle=0.5*pi-azimuthin
  fn(mgrid(is:ifull+is-1))%vangle=acos(cosin)
  fn(mgrid(is:ifull+is-1))%ctime=ctimein
end where

return
end subroutine atebnewangle1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This version of tebnewangle is for CCAM and TAPM
!

subroutine atebccangle(is,ifull,cosin,rlon,rlat,fjd,slag,dt,sdlt)

implicit none

integer, intent(in) :: is,ifull
real, intent(in) :: fjd,slag,dt,sdlt
real cdlt
real, dimension(ifull), intent(in) :: cosin,rlon,rlat
real, dimension(ifull) :: hloc,x,y

! cosin = cosine of zenith angle
! rlon = longitude
! rlat = latitude
! fjd = day of year
! slag = sun lag angle
! sdlt = sin declination of sun

if (ufull.eq.0) return

cdlt=sqrt(min(max(1.-sdlt*sdlt,0.),1.))

where (mgrid(is:ifull+is-1).ge.1)
  ! from CCAM zenith.f
  hloc=2.*pi*fjd+slag+pi+rlon+dt*pi/86400.
  ! estimate azimuth angle
  x=sin(-hloc)*cdlt
  y=-cos(-hloc)*cdlt*sin(rlat)+cos(rlat)*sdlt
  !azimuth=atan2(x,y)
  fn(mgrid(is:ifull+is-1))%hangle=0.5*pi-atan2(x,y)
  fn(mgrid(is:ifull+is-1))%vangle=acos(cosin)
  fn(mgrid(is:ifull+is-1))%ctime=min(max(mod(0.5*hloc/pi-0.5,1.),0.),1.)
end where

return
end subroutine atebccangle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calcuates screen level diagnostics
!

subroutine atebscrnout(ifull,tscrn,qscrn,uscrn,u10,diag,raw)

implicit none

integer, intent(in) :: ifull,diag
real, dimension(ifull), intent(inout) :: tscrn,qscrn,uscrn,u10
logical, intent(in), optional :: raw
logical mode

if (ufull.eq.0) return

mode=.false.
if (present(raw)) mode=raw

if (mode) then
  tscrn(ugrid)=pg%tscrn
  qscrn(ugrid)=pg%qscrn
  uscrn(ugrid)=pg%uscrn
  u10(ugrid)=pg%u10
else
  tscrn(ugrid)=sigmau*pg%tscrn+(1.-sigmau)*tscrn(ugrid)
  qscrn(ugrid)=sigmau*pg%qscrn+(1.-sigmau)*qscrn(ugrid)
  uscrn(ugrid)=sigmau*pg%uscrn+(1.-sigmau)*uscrn(ugrid)
  u10(ugrid)=sigmau*pg%u10+(1.-sigmau)*u10(ugrid)
end if

return
end subroutine atebscrnout

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main routine for calculating urban flux contrabution

! ifull = number of horizontal grid points
! dt = model time step (sec)
! zmin = first model level height (m)
! sg = incomming short wave radiation (W/m^2)
! rg = incomming long wave radiation (W/m^2)
! rnd = incomming rainfall/snowfall rate (kg/(m^2 s))
! rho = atmospheric density at first model level (kg/m^3)
! temp = atmospheric temperature at first model level (K)
! mixr = atmospheric mioxing ratio at first model level (kg/kg)
! ps = surface pressure (Pa)
! pa = pressure at first model level (Pa)
! uu = U component of wind speed at first model level (m/s)
! vv = V component of wind speed at first model level (m/s)
! umin = minimum wind speed (m/s)
! ofg = Input/Output sensible heat flux (W/m^2)
! oeg = Input/Output latient heat flux (W/m^2)
! ots = Input/Output radiative/skin temperature (K)
! owf = Input/Output wetness fraction/surface water (%)
! diag = diagnostic message mode (0=off, 1=basic messages, 2=more detailed messages, etc)

subroutine atebcalc(ifull,ofg,oeg,ots,owf,dt,zmin,sg,rg,rnd,rho,temp,mixr,ps,pa,uu,vv,umin,diag,raw)

implicit none

integer, intent(in) :: ifull,diag
real, intent(in) :: dt,zmin,umin
real, dimension(ifull), intent(in) :: sg,rg,rnd,rho,temp,mixr,ps,pa,uu,vv
real, dimension(ifull), intent(inout) :: ofg,oeg,ots,owf
type (tatm), dimension(ufull) :: atm
type (tout), dimension(ufull) :: uo
logical, intent(in), optional :: raw
logical mode

if (ufull.eq.0) return

mode=.false.
if (present(raw)) mode=raw

atm%sg=sg(ugrid)
atm%rg=rg(ugrid)
atm%rho=rho(ugrid)
atm%temp=temp(ugrid)
atm%mixr=mixr(ugrid)
atm%ps=ps(ugrid)
atm%pa=pa(ugrid)
atm%umag=max(sqrt(uu(ugrid)**2+vv(ugrid)**2),umin)
atm%udir=atan2(vv(ugrid),uu(ugrid))
where (atm%temp.le.273.16) ! diagnose snow
  atm%snd=rnd(ugrid)
  atm%rnd=0.
elsewhere
  atm%rnd=rnd(ugrid)
  atm%snd=0.
end where

call atebeval(uo,dt,atm,zmin,diag)

if (mode) then
  ofg(ugrid)=uo%fg
  oeg(ugrid)=uo%eg
  ots(ugrid)=uo%ts
  owf(ugrid)=uo%wf
else
  ofg(ugrid)=(1.-sigmau)*ofg(ugrid)+sigmau*uo%fg
  oeg(ugrid)=(1.-sigmau)*oeg(ugrid)+sigmau*uo%eg
  ots(ugrid)=((1.-sigmau)*ots(ugrid)**4+sigmau*uo%ts**4)**0.25
  owf(ugrid)=(1.-sigmau)*owf(ugrid)+sigmau*uo%wf
end if

return
end subroutine atebcalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! urban flux calculations

! Basic loop is:
!  Short wave flux (nrefl reflections)
!  Long wave flux (nrefl reflections precomputed)
!  Estimate building roughness length for momentum
!  Canyon aerodynamic resistances
!  Solve canyon snow energy budget
!    Canyon snow temperature
!    Solve vegetation energy budget
!      Vegetation canopy temperature
!      Solve canyon sensible heat budget
!        Canyon temperature
!        Solve canyon latent heat budget
!          Canyon mixing ratio
!        End latent heat budget loop
!      End sensible heat budget loop
!    End vegetation energy budget loop
!  End canyon snow energy budget loop
!  Solve roof snow energy budget
!    Roof snow temperature
!  End roof snow energy budget loop
!  Roof longwave, sensible and latent heat fluxes
!  Update water on canyon surfaces
!  Update snow albedo and density
!  Update urban roof, road and wall temperatures
!  Estimate bulk roughness length for heat
!  Estimate bulk long wave flux and surface temperature
!  Estimate bulk sensible and latent heat fluxes

subroutine atebeval(uo,ddt,atm,zmin,diag)

implicit none

integer, intent(in) :: diag
integer iqu,k,ii,cns,cnr
integer, dimension(ufull) :: igs,igr
real, intent(in) :: ddt,zmin
real, dimension(ufull,2:3) :: ggaroof,ggawall,ggaroad
real, dimension(ufull,1:3) :: ggbroof,ggbwall,ggbroad
real, dimension(ufull,1:2) :: ggcroof,ggcwall,ggcroad
real, dimension(ufull,3) :: garoof,gawalle,gawallw,garoad
real, dimension(ufull) :: garfsn,gardsn
real, dimension(ufull) :: rdsntemp,rfsntemp,rdsnmelt,rfsnmelt
real, dimension(ufull) :: wallpsi,roadpsi,fgtop,egtop,qsatr
real, dimension(ufull) :: oldval,newval,cu,ctmax,ctmin,evctx,evct
real, dimension(ufull) :: ln,rn,we,ww,wr,zolog,a,xe,xw,cuven,n,zom,zonet,dis
real, dimension(ufull) :: ncwa,ncwe,ncww,ncwr,ncra,ncrr,ncrw
real, dimension(ufull) :: rcwa,rcwe,rcww,rcwr,rcra,rcrr,rcrw
real, dimension(ufull) :: p_fgtop,p_wallpsi,p_roadpsi
real, dimension(ufull) :: p_sntemp,p_gasn,p_snmelt
real, dimension(ufull) :: z_on_l
logical, parameter :: caleffzo=.false. ! estimate effective roughness length
type(tatm), dimension(ufull), intent(in) :: atm
type(tatm), dimension(ufull) :: p_atm
type(tout), dimension(ufull), intent(out) :: uo
type(trad), dimension(ufull) :: sg,rg,fg,eg
type(trad), dimension(ufull) :: p_sg,p_rg,p_fg,p_eg
type(trad), dimension(ufull) :: acond,p_acond
type(tdiag), dimension(ufull) :: dg,p_dg
type(tsurf), dimension(ufull) :: p_rodum
type(tvege), dimension(ufull) :: p_vegdum
type(twall), dimension(ufull) :: p_walledum,p_wallwdum
type(tdata), dimension(ufull) :: p_fn
type(tprog), dimension(ufull) :: p_pg

if (diag.ge.1) write(6,*) "Evaluating aTEB"

! null values (for SX compiler)
igs=0
igr=0
acond%roof=0.
acond%road=0.
acond%walle=0.
acond%wallw=0.
acond%veg=0.
acond%rfsn=0.
acond%rdsn=0.
rg=acond
fg=acond
eg=acond
dg%roofdelta=0.
dg%roaddelta=0.
dg%vegdelta=0.
dg%rfsndelta=0.
dg%rdsndelta=0.
dg%tempc=0.
dg%mixrc=0.
dg%tempr=0.
dg%mixrr=0.
dg%sigd=0.
dg%sigr=0.
dg%rfdzmin=0.
dg%accool=0.
dg%canyonrgout=0.
dg%roofrgout=0.
dg%tran=0.
dg%evap=0.
dg%c1=0.
dg%topu=0.
dg%canyontemp=0.
dg%canyonmix=0.

! limit prognostic state variables
do ii=1,3
  roof%temp(ii)=min(max(roof%temp(ii),200.),400.)
  walle%temp(ii)=min(max(walle%temp(ii),200.),400.)
  wallw%temp(ii)=min(max(wallw%temp(ii),200.),400.)
  road%temp(ii)=min(max(road%temp(ii),200.),400.)
end do
veg%moist=min(max(veg%moist,swilt),ssat)
roof%water=min(max(roof%water,0.),maxrfwater)
road%water=min(max(road%water,0.),maxrdwater)
veg%water=min(max(veg%water,0.),maxvgwater)
roof%snow=min(max(roof%snow,0.),maxrfsn)
road%snow=min(max(road%snow,0.),maxrdsn)
roof%den=min(max(roof%den,minsnowden),maxsnowden)
road%den=min(max(road%den,minsnowden),maxsnowden)
roof%alpha=min(max(roof%alpha,minsnowalpha),maxsnowalpha)
road%alpha=min(max(road%alpha,minsnowalpha),maxsnowalpha)

! new snowfall
where (atm%snd.gt.0.)
  roof%den=(roof%snow*roof%den+atm%snd*ddt*minsnowden)/(roof%snow+ddt*atm%snd)
  road%den=(road%snow*road%den+atm%snd*ddt*minsnowden)/(road%snow+ddt*atm%snd)
  roof%alpha=maxsnowalpha
  road%alpha=maxsnowalpha
end where
    
! water and snow cover fractions
dg%roofdelta=(roof%water/maxrfwater)**(2./3.)
dg%roaddelta=(road%water/maxrdwater)**(2./3.)
dg%vegdelta=(veg%water/maxvgwater)**(2./3.)
dg%rfsndelta=roof%snow/(roof%snow+maxrfsn)
dg%rdsndelta=road%snow/(road%snow+maxrdsn)

! canyon (displacement height at refheight*building height)
dg%sigd=atm%ps*exp(-grav*fn%bldheight*refheight/(rd*atm%temp))
dg%tempc=atm%temp*(dg%sigd/atm%pa)**(rd/aircp)
call getqsat(ufull,qsatr,dg%tempc,dg%sigd)
call getqsat(ufull,a,atm%temp,atm%pa)
dg%mixrc=atm%mixr*qsatr/a ! a=qsata

! roof (displacement height at building height)
dg%sigr=atm%ps*exp(-grav*fn%bldheight/(rd*atm%temp))
dg%tempr=atm%temp*(dg%sigr/atm%pa)**(rd/aircp)
call getqsat(ufull,qsatr,dg%tempr,dg%sigr)
dg%mixrr=atm%mixr*qsatr/a ! a=qsata

! total soil depth (for in-canyon vegetation)
dg%totdepth=0.
do ii=1,3
  dg%totdepth=dg%totdepth+fn%roaddepth(ii)
end do

! calculate c1 for soil (for in-canyon vegetation)
n=veg%moist/ssat
where (n.le.0.226)
  dg%c1=10.
elsewhere
  dg%c1=(1.78*n+0.253)/(2.96*n-0.581)
end where

! calculate heat pumped into canyon by air conditioning
if (acmeth.eq.1) then
  dg%accool=max(0.,2.*fn%rooflambda(3)*(roof%temp(3)-fn%bldtemp)/fn%roofdepth(3) &
                  +2.*fn%walllambda(3)*(walle%temp(3)-fn%bldtemp)/fn%walldepth(3) &
                  +2.*fn%walllambda(3)*(wallw%temp(3)-fn%bldtemp)/fn%walldepth(3))
else
  dg%accool=0.
end if

! calculate shortwave radiation
call getswcoeff(ufull,sg,wallpsi,roadpsi,fn,road%alpha,dg%rdsndelta)
sg%roof=(1.-fn%roofalpha)*sg%roof*atm%sg
sg%walle=(1.-fn%wallalpha)*sg%walle*atm%sg
sg%wallw=(1.-fn%wallalpha)*sg%wallw*atm%sg
sg%road=(1.-fn%roadalpha)*sg%road*atm%sg
sg%veg=(1.-fn%vegalpha)*sg%veg*atm%sg
sg%rfsn=(1.-roof%alpha)*sg%rfsn*atm%sg
sg%rdsn=(1.-road%alpha)*sg%rdsn*atm%sg

! Calculate long wave reflections to nrefl order (pregenerated before solvecanyon subroutine)
dg%netemiss=dg%rdsndelta*snowemiss+(1.-dg%rdsndelta)*((1.-fn%sigmaveg)*fn%roademiss+fn%sigmaveg*fn%vegemiss)
dg%cwa=wallpsi
dg%cra=roadpsi
dg%cwe=0.
dg%cww=1.-2.*wallpsi
dg%crw=0.5*(1.-roadpsi)
dg%crr=0.
dg%cwr=wallpsi
rcwa=dg%cwa
rcra=dg%cra
rcwe=dg%cwe
rcww=dg%cww
rcrw=dg%crw
rcrr=dg%crr
rcwr=dg%cwr
do k=1,nrefl
  ncwa=(1.-dg%netemiss)*wallpsi*rcra+(1.-fn%wallemiss)*(1.-2.*wallpsi)*rcwa
  ncra=(1.-fn%wallemiss)*(1.-roadpsi)*rcwa
  ncwe=(1.-dg%netemiss)*wallpsi*rcrw+(1.-fn%wallemiss)*(1.-2.*wallpsi)*rcww
  ncww=(1.-dg%netemiss)*wallpsi*rcrw+(1.-fn%wallemiss)*(1.-2.*wallpsi)*rcwe
  ncrw=(1.-fn%wallemiss)*(1.-roadpsi)*0.5*(rcww+rcwe)  
  ncwr=(1.-dg%netemiss)*wallpsi*rcrr+(1.-fn%wallemiss)*(1.-2.*wallpsi)*rcwr
  ncrr=(1.-fn%wallemiss)*(1.-roadpsi)*rcwr
  rcwa=ncwa
  rcra=ncra
  rcwe=ncwe
  rcww=ncww
  rcrw=ncrw
  rcrr=ncrr
  rcwr=ncwr
  dg%cwa=dg%cwa+rcwa
  dg%cra=dg%cra+rcra
  dg%cwe=dg%cwe+rcwe
  dg%cww=dg%cww+rcww
  dg%crw=dg%crw+rcrw
  dg%cwr=dg%cwr+rcwr
  dg%crr=dg%crr+rcrr
end do

! calculate in-canyon roughness length
zolog=1./sqrt(dg%rdsndelta/log(0.1*fn%bldheight/zosnow)**2 &
     +(1.-dg%rdsndelta)*(fn%sigmaveg/log(0.1*fn%bldheight/zoveg)**2 &
     +(1.-fn%sigmaveg)/log(0.1*fn%bldheight/zocanyon)**2))
zonet=0.1*fn%bldheight*exp(-1./zolog)

! Estimate urban roughness length
zom=zomratio*fn%bldheight
n=road%snow/(road%snow+maxrdsn+0.408*grav*zom)                 ! snow cover for urban roughness calc (Douville, et al 1995)
zom=(1.-n)*zom+n*zosnow                                        ! blend urban and snow roughness lengths (i.e., snow fills canyon)
dg%rfdzmin=max(abs(zmin-fn%bldheight),zocanyon+1.)             ! distance to roof displacement height
where (zmin.ge.fn%bldheight)
  pg%lzom=log(zmin/zom)
elsewhere ! lowest atmospheric model level is within the canopy.  Need to interact with boundary layer scheme.
  pg%lzom=log(fn%bldheight/zom)*exp(0.5*fn%hwratio*(1.-zmin/fn%bldheight))
end where
pg%cndzmin=zom*exp(pg%lzom)                                     ! distance to canyon displacement height

! calculate canyon wind speed and bulk transfer coefficents
select case(resmeth)
  case(0) ! Masson (2000)
    cu=exp(-0.25*fn%hwratio)
    acond%road=cu ! bulk transfer coefficents are updated in solvecanyon
    acond%walle=cu
    acond%wallw=cu
    acond%rdsn=cu
    acond%veg=cu
  case(1) ! Harman et al (2004) - modified to include 2nd wall
    call getincanwind(we,ww,wr,atm%udir,fn%bldheight,fn%bldheight/fn%hwratio,zonet)
    dis=max(max(max(0.1*fn%bldheight,zocanyon+0.1),zoveg+0.1),zosnow+0.1)
    zolog=log(dis/zocanyon)
    a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
    acond%road=a*wr                  ! road bulk transfer
    acond%walle=a*we                 ! east wall bulk transfer
    acond%wallw=a*ww                 ! west wall bulk transfer
    zolog=log(dis/zoveg)
    a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
    acond%veg=a*wr
    zolog=log(dis/zosnow)
    a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
    acond%rdsn=a*wr                  ! road snow bulk transfer
  case(2) ! Kusaka et al (2001)
    cu=exp(-0.25*fn%hwratio)
    where (cu.le.5.)
      acond%road=6.15+4.18*cu
    elsewhere
      acond%road=7.51*cu**0.78
    end where
    acond%walle=acond%road
    acond%wallw=acond%road
    acond%rdsn=acond%road
    acond%veg=acond%road
end select

! solve for road snow temperature -------------------------------
! includes solution to vegetation canopy temperature, canopy temperature, sensible heat flux and longwave radiation
cns=0
cnr=0
do iqu=1,ufull ! prepare index arrays
  if (dg(iqu)%rdsndelta.gt.0.) then
    cns=cns+1
    igs(cns)=iqu
  else
    cnr=cnr+1
    igr(cnr)=iqu
  end if
end do
if (cns.gt.0) then ! road snow
  p_rodum(1:cns)=road(igs(1:cns)) ! pack
  p_walledum(1:cns)=walle(igs(1:cns))
  p_wallwdum(1:cns)=wallw(igs(1:cns))
  p_vegdum(1:cns)=veg(igs(1:cns))
  p_dg(1:cns)=dg(igs(1:cns))
  p_sg(1:cns)=sg(igs(1:cns))
  p_rg(1:cns)=rg(igs(1:cns))
  p_fg(1:cns)=fg(igs(1:cns))
  p_eg(1:cns)=eg(igs(1:cns))
  p_atm(1:cns)=atm(igs(1:cns))
  p_acond(1:cns)=acond(igs(1:cns))
  p_wallpsi(1:cns)=wallpsi(igs(1:cns))
  p_roadpsi(1:cns)=roadpsi(igs(1:cns))
  p_fn(1:cns)=fn(igs(1:cns))
  p_pg(1:cns)=pg(igs(1:cns))
  ctmax(1:cns)=max(p_dg(1:cns)%tempc,p_rodum(1:cns)%temp(1),p_walledum(1:cns)%temp(1), &
                   p_wallwdum(1:cns)%temp(1),p_pg(1:cns)%vegtemp)+10. ! max road snow temp
  ctmin(1:cns)=min(p_dg(1:cns)%tempc,p_rodum(1:cns)%temp(1),p_walledum(1:cns)%temp(1), &
                   p_wallwdum(1:cns)%temp(1),p_pg(1:cns)%vegtemp)-10. ! min road snow temp
  p_sntemp(1:cns)=0.75*ctmax(1:cns)+0.25*ctmin(1:cns)
  call solverdsn(cns,evct(1:cns),p_rg(1:cns),p_fg(1:cns),p_fgtop(1:cns),p_eg(1:cns), &
         p_gasn(1:cns),p_snmelt(1:cns),p_sntemp(1:cns),p_rodum(1:cns),p_walledum(1:cns), &
         p_wallwdum(1:cns),p_vegdum(1:cns),p_dg(1:cns),p_sg(1:cns),p_atm(1:cns),ddt, &
         p_acond(1:cns),p_wallpsi(1:cns),p_roadpsi(1:cns),p_fn(1:cns),p_pg(1:cns))
  oldval(1:cns)=p_sntemp(1:cns)
  p_sntemp(1:cns)=0.25*ctmax(1:cns)+0.75*ctmin(1:cns)
  do k=1,nfgits ! sectant
    evctx(1:cns)=evct(1:cns)
    call solverdsn(cns,evct(1:cns),p_rg(1:cns),p_fg(1:cns),p_fgtop(1:cns),p_eg(1:cns), &
           p_gasn(1:cns),p_snmelt(1:cns),p_sntemp(1:cns),p_rodum(1:cns),p_walledum(1:cns), &
           p_wallwdum(1:cns),p_vegdum(1:cns),p_dg(1:cns),p_sg(1:cns),p_atm(1:cns),ddt, &
           p_acond(1:cns),p_wallpsi(1:cns),p_roadpsi(1:cns),p_fn(1:cns),p_pg(1:cns))
    evctx(1:cns)=evct(1:cns)-evctx(1:cns)
    if (all(abs(evctx(1:cns)).le.tol)) exit
    where (abs(evctx(1:cns)).gt.tol)
      newval(1:cns)=p_sntemp(1:cns)-alpha*evct(1:cns)*(p_sntemp(1:cns)-oldval(1:cns))/evctx(1:cns)
      oldval(1:cns)=p_sntemp(1:cns)
      p_sntemp(1:cns)=newval(1:cns)
    end where
    p_sntemp(1:cns)=min(max(p_sntemp(1:cns),ctmin(1:cns)),ctmax(1:cns))      
  end do    
  rg(igs(1:cns))=p_rg(1:cns) ! unpack
  fg(igs(1:cns))=p_fg(1:cns)
  fgtop(igs(1:cns))=p_fgtop(1:cns)
  eg(igs(1:cns))=p_eg(1:cns)
  pg(igs(1:cns))=p_pg(1:cns)
  gardsn(igs(1:cns))=p_gasn(1:cns)
  rdsnmelt(igs(1:cns))=p_snmelt(1:cns)
  rdsntemp(igs(1:cns))=p_sntemp(1:cns)
  acond(igs(1:cns))=p_acond(1:cns)
  dg(igs(1:cns))=p_dg(1:cns)
end if
if (cnr.gt.0) then ! no road snow
  p_rodum(1:cnr)=road(igr(1:cnr)) ! pack
  p_walledum(1:cnr)=walle(igr(1:cnr))
  p_wallwdum(1:cnr)=wallw(igr(1:cnr))
  p_vegdum(1:cnr)=veg(igr(1:cnr))
  p_dg(1:cnr)=dg(igr(1:cnr))
  p_sg(1:cnr)=sg(igr(1:cnr))
  p_rg(1:cnr)=rg(igr(1:cnr))
  p_fg(1:cnr)=fg(igr(1:cnr))
  p_eg(1:cnr)=eg(igr(1:cnr))
  p_atm(1:cnr)=atm(igr(1:cnr))
  p_acond(1:cnr)=acond(igr(1:cnr))
  p_wallpsi(1:cnr)=wallpsi(igr(1:cnr))
  p_roadpsi(1:cnr)=roadpsi(igr(1:cnr))
  p_sntemp(1:cnr)=p_rodum(1:cnr)%temp(1)
  p_fn(1:cnr)=fn(igr(1:cnr))
  p_pg(1:cnr)=pg(igr(1:cnr))
  call solverdsn(cnr,evct(1:cnr),p_rg(1:cnr),p_fg(1:cnr),p_fgtop(1:cnr),p_eg(1:cnr), &
         p_gasn(1:cnr),p_snmelt(1:cnr),p_sntemp(1:cnr),p_rodum(1:cnr),p_walledum(1:cnr), &
         p_wallwdum(1:cnr),p_vegdum(1:cnr),p_dg(1:cnr),p_sg(1:cnr),p_atm(1:cnr),ddt, &
         p_acond(1:cnr),p_wallpsi(1:cnr),p_roadpsi(1:cnr),p_fn(1:cnr),p_pg(1:cnr))
  rg(igr(1:cnr))=p_rg(1:cnr) ! unpack
  fg(igr(1:cnr))=p_fg(1:cnr)
  fgtop(igr(1:cnr))=p_fgtop(1:cnr)
  eg(igr(1:cnr))=p_eg(1:cnr)
  pg(igr(1:cnr))=p_pg(1:cnr)
  gardsn(igr(1:cnr))=0.
  rdsnmelt(igr(1:cnr))=0.
  rdsntemp(igr(1:cnr))=road(igr(1:cnr))%temp(1)
  acond(igr(1:cnr))=p_acond(1:cnr)
  dg(igr(1:cnr))=p_dg(1:cnr)
end if
! ---------------------------------------------------------------    

! solve for roof snow temperature -------------------------------
cns=0
cnr=0
rg%rfsn=0.
fg%rfsn=0.
eg%rfsn=0.
garfsn=0.
rfsnmelt=0.
rfsntemp=roof%temp(1)
acond%rfsn=0.
do iqu=1,ufull ! prepare index arrays
  if (dg(iqu)%rfsndelta.gt.0.) then
    cns=cns+1
    igs(cns)=iqu
  else
    cnr=cnr+1
    igr(cnr)=iqu
  end if
end do
if (cns.gt.0) then ! roof snow
  p_rodum(1:cns)=roof(igs(1:cns)) ! pack
  p_dg(1:cns)=dg(igs(1:cns))
  p_sg(1:cns)=sg(igs(1:cns))
  p_rg(1:cns)=rg(igs(1:cns))
  p_fg(1:cns)=fg(igs(1:cns))
  p_eg(1:cns)=eg(igs(1:cns))
  p_atm(1:cns)=atm(igs(1:cns))
  p_acond(1:cns)=acond(igs(1:cns))
  p_fn(1:cns)=fn(igs(1:cns))
  ctmax(1:cns)=max(p_dg(1:cns)%tempr,p_rodum(1:cns)%temp(1))+10. ! max roof snow temp
  ctmin(1:cns)=min(p_dg(1:cns)%tempr,p_rodum(1:cns)%temp(1))-10. ! min roof snow temp
  p_sntemp(1:cns)=0.75*ctmax(1:cns)+0.25*ctmin(1:cns)
  call solverfsn(cns,evct(1:cns),p_rg(1:cns),p_fg(1:cns),p_eg(1:cns),p_gasn(1:cns),p_snmelt(1:cns),p_sntemp(1:cns), &
                 p_rodum(1:cns),p_dg(1:cns),p_sg(1:cns),p_atm(1:cns),p_acond(1:cns),p_fn(1:cns),ddt)
  oldval(1:cns)=p_sntemp(1:cns)
  p_sntemp(1:cns)=0.25*ctmax(1:cns)+0.75*ctmin(1:cns)
  do k=1,nfgits ! sectant
    evctx(1:cns)=evct(1:cns)
    call solverfsn(cns,evct(1:cns),p_rg(1:cns),p_fg(1:cns),p_eg(1:cns),p_gasn(1:cns),p_snmelt(1:cns),p_sntemp(1:cns), &
                   p_rodum(1:cns),p_dg(1:cns),p_sg(1:cns),p_atm(1:cns),p_acond(1:cns),p_fn(1:cns),ddt)
    evctx(1:cns)=evct(1:cns)-evctx(1:cns)
    if (all(abs(evctx(1:cns)).le.tol)) exit
    where (abs(evctx(1:cns)).gt.tol)
      newval(1:cns)=p_sntemp(1:cns)-alpha*evct(1:cns)*(p_sntemp(1:cns)-oldval(1:cns))/evctx(1:cns)
      oldval(1:cns)=p_sntemp(1:cns)
      p_sntemp(1:cns)=newval(1:cns)
    end where
    p_sntemp(1:cns)=min(max(p_sntemp(1:cns),ctmin(1:cns)),ctmax(1:cns))      
  end do
  rg(igs(1:cns))%rfsn=p_rg(1:cns)%rfsn ! unpack
  fg(igs(1:cns))%rfsn=p_fg(1:cns)%rfsn
  eg(igs(1:cns))%rfsn=p_eg(1:cns)%rfsn
  garfsn(igs(1:cns))=p_gasn(1:cns)
  rfsnmelt(igs(1:cns))=p_snmelt(1:cns)
  acond(igs(1:cns))%rfsn=p_acond(1:cns)%rfsn
end if
!---------------------------------------------------------------- 

! calculate roof sensible and latent heat fluxes (without snow)
rg%roof=fn%roofemiss*(atm%rg-sbconst*roof%temp(1)**4)
dg%roofrgout=sbconst*(dg%rfsndelta*snowemiss*rfsntemp**4+(1.-dg%rfsndelta)*fn%roofemiss*roof%temp(1)**4) &
                    +(1.-dg%rfsndelta*snowemiss-(1.-dg%rfsndelta)*fn%roofemiss)*atm%rg
a=log(dg%rfdzmin/zocanyon)
xw=(dg%roofdelta*(1.-dg%rfsndelta)+dg%rfsndelta)*qsatr ! roof surface mixing ratio
we=roof%temp(1)*(1.+0.61*xw)
ww=dg%tempr*(1.+0.61*dg%mixrr)
! n and xe are dummy variables for cd and lzohroof
call getinvres(ufull,acond%roof,n,z_on_l,xe,a,dg%rfdzmin,we,ww,atm%umag,1) ! Assume zot=0.1*zom (i.e., Kanda et al 2007, small experiment)
acond%roof=acond%roof/atm%umag
fg%roof=aircp*atm%rho*(roof%temp(1)-dg%tempr)*acond%roof*atm%umag
call getqsat(ufull,qsatr,roof%temp(1),dg%sigr)
where (qsatr.lt.dg%mixrr)
  eg%roof=lv*atm%rho*(qsatr-dg%mixrr)*acond%roof*atm%umag
elsewhere
  ! MJT suggestion for max eg   
  eg%roof=lv*min(atm%rho*dg%roofdelta*(qsatr-dg%mixrr)*acond%roof*atm%umag,roof%water/ddt+atm%rnd+rfsnmelt)
end where

! join two walls into a single wall (testing only)
if (resmeth.eq.0.or.resmeth.eq.2) then
  do k=1,3
    n=0.5*(walle%temp(k)+wallw%temp(k))
    walle%temp(k)=n
    wallw%temp(k)=n
  end do
  n=0.5*(sg%walle+sg%wallw)
  sg%walle=n
  sg%wallw=n
  n=0.5*(rg%walle+rg%wallw)
  rg%walle=n
  rg%wallw=n
  n=0.5*(fg%walle+fg%wallw)
  fg%walle=n
  fg%wallw=n
  n=0.5*(eg%walle+eg%wallw)
  eg%walle=n
  eg%wallw=n
end if

! tridiagonal solver coefficents for calculating roof, road and wall temperatures
do k=2,3
  ggaroof(:,k)=-2./(fn%roofdepth(k-1)/fn%rooflambda(k-1)+fn%roofdepth(k)/fn%rooflambda(k))
  ggawall(:,k)=-2./(fn%walldepth(k-1)/fn%walllambda(k-1)+fn%walldepth(k)/fn%walllambda(k))
  ggaroad(:,k)=-2./(fn%roaddepth(k-1)/fn%roadlambda(k-1)+fn%roaddepth(k)/fn%roadlambda(k))
end do
ggbroof(:,1)=2./(fn%roofdepth(1)/fn%rooflambda(1)+fn%roofdepth(2)/fn%rooflambda(2))+fn%roofcp(1)*fn%roofdepth(1)/ddt
ggbwall(:,1)=2./(fn%walldepth(1)/fn%walllambda(1)+fn%walldepth(2)/fn%walllambda(2))+fn%wallcp(1)*fn%walldepth(1)/ddt
ggbroad(:,1)=2./(fn%roaddepth(1)/fn%roadlambda(1)+fn%roaddepth(2)/fn%roadlambda(2))+fn%roadcp(1)*fn%roaddepth(1)/ddt
ggbroof(:,2)=2./(fn%roofdepth(1)/fn%rooflambda(1)+fn%roofdepth(2)/fn%rooflambda(2)) &
             +2./(fn%roofdepth(2)/fn%rooflambda(2)+fn%roofdepth(3)/fn%rooflambda(3)) &
             +fn%roofcp(2)*fn%roofdepth(2)/ddt
ggbwall(:,2)=2./(fn%walldepth(1)/fn%walllambda(1)+fn%walldepth(2)/fn%walllambda(2)) &
             +2./(fn%walldepth(2)/fn%walllambda(2)+fn%walldepth(3)/fn%walllambda(3)) &
             +fn%wallcp(2)*fn%walldepth(2)/ddt
ggbroad(:,2)=2./(fn%roaddepth(1)/fn%roadlambda(1)+fn%roaddepth(2)/fn%roadlambda(2)) &
             +2./(fn%roaddepth(2)/fn%roadlambda(2)+fn%roaddepth(3)/fn%roadlambda(3)) &
             +fn%roadcp(2)*fn%roaddepth(2)/ddt
ggbroof(:,3)=2./(fn%roofdepth(2)/fn%rooflambda(2)+fn%roofdepth(3)/fn%rooflambda(3)) &
             +2.*fn%rooflambda(3)/fn%roofdepth(3)+fn%roofcp(3)*fn%roofdepth(3)/ddt
ggbwall(:,3)=2./(fn%walldepth(2)/fn%walllambda(2)+fn%walldepth(3)/fn%walllambda(3)) &
             +2.*fn%walllambda(3)/fn%walldepth(3)+fn%wallcp(3)*fn%walldepth(3)/ddt
ggbroad(:,3)=2./(fn%roaddepth(2)/fn%roadlambda(2)+fn%roaddepth(3)/fn%roadlambda(3))+fn%roadcp(3)*fn%roaddepth(3)/ddt
do k=1,2
  ggcroof(:,k)=-2./(fn%roofdepth(k)/fn%rooflambda(k)+fn%roofdepth(k+1)/fn%rooflambda(k+1))
  ggcwall(:,k)=-2./(fn%walldepth(k)/fn%walllambda(k)+fn%walldepth(k+1)/fn%walllambda(k+1))
  ggcroad(:,k)=-2./(fn%roaddepth(k)/fn%roadlambda(k)+fn%roaddepth(k+1)/fn%roadlambda(k+1))
end do
garoof(:,1)=(1.-dg%rfsndelta)*(sg%roof+rg%roof-fg%roof-eg%roof) &
            +dg%rfsndelta*garfsn+roof%temp(1)*fn%roofcp(1)*fn%roofdepth(1)/ddt
gawalle(:,1)=sg%walle+rg%walle-fg%walle-eg%walle+walle%temp(1)*fn%wallcp(1)*fn%walldepth(1)/ddt
gawallw(:,1)=sg%wallw+rg%wallw-fg%wallw-eg%wallw+wallw%temp(1)*fn%wallcp(1)*fn%walldepth(1)/ddt
garoad(:,1)=(1.-dg%rdsndelta)*(sg%road+rg%road-fg%road-eg%road) &
            +dg%rdsndelta*gardsn+road%temp(1)*fn%roadcp(1)*fn%roaddepth(1)/ddt
garoof(:,2)=roof%temp(2)*fn%roofcp(2)*fn%roofdepth(2)/ddt
gawalle(:,2)=walle%temp(2)*fn%wallcp(2)*fn%walldepth(2)/ddt
gawallw(:,2)=wallw%temp(2)*fn%wallcp(2)*fn%walldepth(2)/ddt
garoad(:,2)=road%temp(2)*fn%roadcp(2)*fn%roaddepth(2)/ddt
garoof(:,3)=2.*fn%rooflambda(3)*fn%bldtemp/fn%roofdepth(3)+roof%temp(3)*fn%roofcp(3)*fn%roofdepth(3)/ddt
gawalle(:,3)=2.*fn%walllambda(3)*fn%bldtemp/fn%walldepth(3)+walle%temp(3)*fn%wallcp(3)*fn%walldepth(3)/ddt
gawallw(:,3)=2.*fn%walllambda(3)*fn%bldtemp/fn%walldepth(3)+wallw%temp(3)*fn%wallcp(3)*fn%walldepth(3)/ddt
garoad(:,3)=road%temp(3)*fn%roadcp(3)*fn%roaddepth(3)/ddt
    
! tridiagonal solver (Thomas algorithm) to solve for roof, road and wall temperatures
do ii=2,3
  n=ggaroof(:,ii)/ggbroof(:,ii-1)
  ggbroof(:,ii)=ggbroof(:,ii)-n*ggcroof(:,ii-1)
  garoof(:,ii)=garoof(:,ii)-n*garoof(:,ii-1)
  n=ggawall(:,ii)/ggbwall(:,ii-1)
  ggbwall(:,ii)=ggbwall(:,ii)-n*ggcwall(:,ii-1)
  gawalle(:,ii)=gawalle(:,ii)-n*gawalle(:,ii-1)
  gawallw(:,ii)=gawallw(:,ii)-n*gawallw(:,ii-1)
  n=ggaroad(:,ii)/ggbroad(:,ii-1)
  ggbroad(:,ii)=ggbroad(:,ii)-n*ggcroad(:,ii-1)
  garoad(:,ii)=garoad(:,ii)-n*garoad(:,ii-1)
end do
roof%temp(3)=garoof(:,3)/ggbroof(:,3)
walle%temp(3)=gawalle(:,3)/ggbwall(:,3)
wallw%temp(3)=gawallw(:,3)/ggbwall(:,3)
road%temp(3)=garoad(:,3)/ggbroad(:,3)
do ii=2,1,-1
  roof%temp(ii)=(garoof(:,ii)-ggcroof(:,ii)*roof%temp(ii+1))/ggbroof(:,ii)
  walle%temp(ii)=(gawalle(:,ii)-ggcwall(:,ii)*walle%temp(ii+1))/ggbwall(:,ii)
  wallw%temp(ii)=(gawallw(:,ii)-ggcwall(:,ii)*wallw%temp(ii+1))/ggbwall(:,ii)
  road%temp(ii)=(garoad(:,ii)-ggcroad(:,ii)*road%temp(ii+1))/ggbroad(:,ii)
end do

! calculate water for canyon surfaces
n=max(atm%rnd-max(dg%evap/lv,0.)-max(maxvgwater-veg%water,0.)/ddt,0.) ! rainfall reaching the soil under vegetation
! note that since sigmaf=1, then there is no soil evaporation, only transpiration.  Evaporation only occurs from water on leaves.
veg%moist=veg%moist+ddt*dg%c1*(max(n+rdsnmelt*road%den/waterden,0.)-dg%tran/lv)/(waterden*dg%totdepth) ! fix units to kg/m^2
roof%water=roof%water+ddt*(atm%rnd-eg%roof/lv+rfsnmelt)
road%water=road%water+ddt*(atm%rnd-eg%road/lv+rdsnmelt)
veg%water=veg%water+ddt*(atm%rnd-dg%evap/lv)

! calculate snow
roof%snow=roof%snow+ddt*(atm%snd-eg%rfsn/ls-rfsnmelt)
road%snow=road%snow+ddt*(atm%snd-eg%rdsn/ls-rdsnmelt)
roof%den=roof%den+ddt*(0.24/86400.)*(maxsnowden-roof%den)
road%den=road%den+ddt*(0.24/86400.)*(maxsnowden-road%den)
where (rfsntemp.lt.273.16)
  roof%alpha=roof%alpha+ddt*(-0.008/86400.)
elsewhere
  roof%alpha=roof%alpha+ddt*(0.24/86400.)*(minsnowalpha-roof%alpha)
end where
where (rdsntemp.lt.273.16)
  road%alpha=road%alpha+ddt*(-0.008/86400.)
elsewhere
  road%alpha=road%alpha+ddt*(0.24/86400.)*(minsnowalpha-road%alpha)
end where

! limit temperatures to sensible values
do ii=1,3
  roof%temp(ii)=min(max(roof%temp(ii),200.),400.)
  walle%temp(ii)=min(max(walle%temp(ii),200.),400.)
  wallw%temp(ii)=min(max(wallw%temp(ii),200.),400.)
  road%temp(ii)=min(max(road%temp(ii),200.),400.)
end do
veg%moist=min(max(veg%moist,swilt),ssat)
roof%water=min(max(roof%water,0.),maxrfwater)
road%water=min(max(road%water,0.),maxrdwater)
veg%water=min(max(veg%water,0.),maxvgwater)
roof%snow=min(max(roof%snow,0.),maxrfsn)
road%snow=min(max(road%snow,0.),maxrdsn)
roof%den=min(max(roof%den,minsnowden),maxsnowden)
road%den=min(max(road%den,minsnowden),maxsnowden)
roof%alpha=min(max(roof%alpha,minsnowalpha),maxsnowalpha)
road%alpha=min(max(road%alpha,minsnowalpha),maxsnowalpha)

! combine snow and snow-free tiles
fg%roof=dg%rfsndelta*fg%rfsn+(1.-dg%rfsndelta)*fg%roof ! redefine as net roof fg
eg%roof=dg%rfsndelta*eg%rfsn+(1.-dg%rfsndelta)*eg%roof ! redefine as net roof eg
egtop=dg%rdsndelta*eg%rdsn+(1.-dg%rdsndelta)*((1.-fn%sigmaveg)*eg%road+fn%sigmaveg*eg%veg)

! calculate longwave, sensible heat latent heat and surface water (i.e., bc for RH) outputs
! estimate surface temp from outgoing longwave radiation
uo%ts=((fn%sigmabld*dg%roofrgout+(1.-fn%sigmabld)*dg%canyonrgout)/sbconst)**0.25
uo%fg=fn%sigmabld*fg%roof+(1.-fn%sigmabld)*fgtop+fn%industryfg
uo%eg=fn%sigmabld*eg%roof+(1.-fn%sigmabld)*egtop
n=max(min((veg%moist-swilt)/(sfc-swilt),1.),0.) ! veg wetfac (see sflux.f or cable_canopy.f90)
uo%wf=fn%sigmabld*(1.-dg%rfsndelta)*dg%roofdelta+(1.-fn%sigmabld)*(1.-dg%rdsndelta)* &
      ((1.-fn%sigmaveg)*dg%roaddelta+fn%sigmaveg*((1.-dg%vegdelta)*n+dg%vegdelta))

! (re)calculate heat roughness length for MOST (diagnostic only)
call getqsat(ufull,a,uo%ts,dg%sigd)
a=a*uo%wf
we=uo%ts*(1.+0.61*a)
ww=dg%tempc*(1.+0.61*dg%mixrc)
select case(zohmeth)
  case(0) ! Use veg formulation
    pg%lzoh=2.3+pg%lzom
    call getinvres(ufull,n,pg%cduv,z_on_l,pg%lzoh,pg%lzom,pg%cndzmin,we,ww,atm%umag,4)
  case(1) ! Use Kanda parameterisation
    call getinvres(ufull,n,pg%cduv,z_on_l,pg%lzoh,pg%lzom,pg%cndzmin,we,ww,atm%umag,2)
    ! Adjust roughness length for heat to account for in-canyon vegetation
    pg%lzoh=1./((1.-fn%sigmaveg*(1.-fn%sigmabld))/pg%lzoh+fn%sigmaveg*(1.-fn%sigmabld)/(2.3+pg%lzom))
    ! (re)calculate drag coeff
    call getinvres(ufull,n,pg%cduv,z_on_l,pg%lzoh,pg%lzom,pg%cndzmin,we,ww,atm%umag,4)
  case(2) ! Use Kanda parameterisation
    pg%lzoh=1./((1.-fn%sigmaveg*(1.-fn%sigmabld))/(6.+pg%lzom)+fn%sigmaveg*(1.-fn%sigmabld)/(2.3+pg%lzom)) ! Kanda (2005)
    call getinvres(ufull,n,pg%cduv,z_on_l,pg%lzoh,pg%lzom,pg%cndzmin,we,ww,atm%umag,4)
end select

! calculate screen level diagnostics
call scrncalc(uo%ts,atm%temp,a,atm%mixr,atm%umag,dg,rdsntemp,veg%moist,zonet)

!if (caleffzo) then ! verification
!  if ((uo(iqut)%ts-dg(iqut)%tempc)*uo(iqut)%fg.le.0.) then
!    write(6,*) "NO SOLUTION iqut ",iqut
!    write(62,*) xw(1),50.
!  else
!    oldval(1)=2.3+pg(iqut)%lzom
!    call getinvres(1,a(1),xe(1),z_on_l,oldval(1),pg(iqut)%lzom,pg(iqut)%cndzmin,uo(iqut)%ts,dg(iqut)%tempc,atm(iqut)%umag,4)
!    evct(1)=aircp*atm(iqut)%rho*(uo(iqut)%ts-dg(iqut)%tempc)*a(1)-uo(iqut)%fg
!    write(6,*) "iqut,its,adj,bal ",iqut,0,oldval(1)-pg(iqut)%lzom,evct(1)
!    n(1)=6.+pg(iqut)%lzom
!    do k=1,20 ! sectant
!      evctx(1)=evct(1)
!      call getinvres(1,a(1),xe(1),z_on_l,n(1),pg(iqut)%lzom,pg(iqut)%cndzmin,uo(iqut)%ts,dg(iqut)%tempc,atm(iqut)%umag,4)  
!      evct(1)=aircp*atm(iqut)%rho*(uo(iqut)%ts-dg(iqut)%tempc)*a(1)-uo(iqut)%fg
!      write(6,*) "iqut,its,adj,bal ",iqut,k,n(1)-pg(iqut)%lzom,evct(1)  
!      evctx(1)=evct(1)-evctx(1)
!      if (abs(evctx(1)).le.tol) exit
!      newval(1)=n(1)-evct(1)*(n(1)-oldval(1))/evctx(1)
!      newval(1)=min(max(newval(1),0.1),50.+pg(iqut)%lzom)
!      oldval(1)=n(1)
!      n(1)=newval(1)
!    end do
!    xw(1)=max(sqrt(pg(iqut)%cduv)*atm(iqut)%umag*zmin*exp(-pg(iqut)%lzom)/1.461E-5,10.)
!    write(62,*) xw(1),n(1)-pg(iqut)%lzom
!  end if
!end if
  
return
end subroutine atebeval

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Estimate saturation mixing ratio (from TAPM)

subroutine getqsat(ifull,qsat,temp,ps)

implicit none

integer, intent(in) :: ifull
real, dimension(ifull), intent(in) :: temp,ps
real, dimension(ifull), intent(out) :: qsat
real, dimension(ifull) :: esatf

where (temp.ge.273.15)
  esatf = 610.*exp(lv/rv*(1./273.15-1./min(temp,343.)))
elsewhere
  esatf = 610.*exp(ls/rv*(1./273.15-1./max(temp,123.)))
endwhere
qsat = 0.622*esatf/(ps-0.378*esatf)

return
end subroutine getqsat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Louis parameterisation based on CCAM sflux.f version (i.e., from CSIRO9),
! but modified for increased ration between momentum and heat roughness lengths
! Also, included Dyer and Hicks scheme based on TAPM Surf.f

subroutine getinvres(cn,invres,cd,z_on_l,olzoh,ilzom,zmin,sthetav,thetav,umag,mode)

implicit none

integer, intent(in) :: cn,mode

real, dimension(cn), intent(in) :: ilzom,zmin,sthetav,thetav,umag
real, dimension(cn), intent(out) :: invres,cd,z_on_l
real, dimension(cn), intent(inout) :: olzoh
real, dimension(cn) :: af,aft,ri,fm,fh,root,denma,denha,re,lna
real, dimension(cn) :: thetavstar,integralh
real, parameter :: nu = 1.461E-5
!real, parameter :: eta0 = 1.827E-5
!real, parameter :: t0 = 291.15
!real, parameter :: c = 120.
!eta=eta0*((t0+c)/(theta+c))*(theta/t0)**(2./3.)
!nu=eta/rho

! use Louis as first guess for Dyer and Hicks scheme (if used)
af=(vkar/ilzom)**2
! umag is now constrained to be above umin in tebcalc
ri=min(grav*zmin*(1.-sthetav/thetav)/umag**2,rimax)
where (ri>0.)
  fm=1./(1.+bprm*ri)**2
elsewhere
  root=sqrt(-ri*exp(ilzom))
  denma=1.+cms*2.*bprm*af*root
  fm=1.-2.*bprm *ri/denma
end where
cd=af*fm

select case(mode) ! roughness length for heat
  case(1) ! zot=zom/10.
    lna=2.3
  case(2) ! Kanda et al 2007
    re=max(sqrt(cd)*umag*zmin*exp(-ilzom)/nu,10.)
    !lna=2.46*re**0.25-2. !(Brutsaet, 1982)
    lna=1.29*re**0.25-2. !(Kanda et al, 2007)
  case(3) ! zot=zom (neglect molecular diffusion)
    lna=0.
  case(4) ! user defined
    lna=olzoh-ilzom
end select
olzoh=lna+ilzom

aft=vkar*vkar/(ilzom*olzoh)
where (ri>0.)
  fh=fm
elsewhere
  denha=1.+chs*2.*bprm*aft*exp(0.5*lna)*root
  fh=1.-2.*bprm*ri/denha
end where
invres=aft*fh*umag
 
thetavstar=aft*fh*(thetav-sthetav)/sqrt(cd)
call dyerhicks(cn,integralh,z_on_l,cd,thetavstar,thetav,sthetav,umag,zmin,ilzom,lna)
invres=(vkar/integralh)*sqrt(cd)*umag

return
end subroutine getinvres

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate stability functions using Dyerhicks

subroutine dyerhicks(cn,integralh,z_on_l,cd,thetavstar,thetav,sthetav,umag,zmin,ilzom,lna)

implicit none

integer, intent(in) :: cn
integer ic
real, dimension(cn), intent(in) :: thetav,sthetav,umag,zmin,ilzom,lna
real, dimension(cn), intent(inout) :: cd,thetavstar
real, dimension(cn), intent(out) :: integralh,z_on_l
real, dimension(cn) :: z0_on_l,zt_on_l,olzoh
real, dimension(cn) :: pm0,ph0,pm1,ph1,integralm
real, parameter :: aa1 = 3.8
real, parameter :: bb1 = 0.5
real, parameter :: cc1 = 0.3

olzoh=ilzom+lna

do ic=1,icmax
  z_on_l=vkar*zmin*grav*thetavstar/(thetav*cd*umag**2)
  z_on_l=min(z_on_l,10.)
  z0_on_l  = z_on_l*exp(-ilzom)
  zt_on_l  = z0_on_l*exp(-lna)
  where (z_on_l.lt.0.)
    pm0     = (1.-16.*z0_on_l)**(-0.25)
    ph0     = (1.-16.*zt_on_l)**(-0.5)
    pm1     = (1.-16.*z_on_l)**(-0.25)
    ph1     = (1.-16.*z_on_l)**(-0.5)
    integralm = ilzom-2.*log((1.+1./pm1)/(1.+1./pm0))-log((1.+1./pm1**2)/(1.+1./pm0**2)) &
               +2.*(atan(1./pm1)-atan(1./pm0))
    integralh = olzoh-2.*log((1.+1./ph1)/(1.+1./ph0))
  elsewhere
    !--------------Beljaars and Holtslag (1991) momentum & heat            
    pm0 = -(a_1*z0_on_l+b_1*(z0_on_l-(c_1/d_1))*exp(-d_1*z0_on_l)+b_1*c_1/d_1)
    pm1 = -(a_1*z_on_l+b_1*(z_on_l-(c_1/d_1))*exp(-d_1*z_on_l)+b_1*c_1/d_1)
    ph0 = -((1.+(2./3.)*a_1*zt_on_l)**1.5+b_1*(zt_on_l-(c_1/d_1))*exp(-d_1*zt_on_l)+b_1*c_1/d_1-1.)
    ph1 = -((1.+(2./3.)*a_1*z_on_l)**1.5+b_1*(z_on_l-(c_1/d_1))*exp(-d_1*z_on_l)+b_1*c_1/d_1-1.)
    integralm = ilzom-(pm1-pm0)    
    integralh = olzoh-(ph1-ph0)         
  endwhere
  !where (z_on_l.le.0.4)
    cd = (max(0.01,min(vkar*umag/integralm,2.))/umag)**2
  !elsewhere
  !  cd = (max(0.01,min(vkar*umag/(aa1*( ( z_on_l**bb1)*(1.0+cc1* z_on_l**(1.-bb1)) &
  !      -(z0_on_l**bb1)*(1.+cc1*z0_on_l**(1.-bb1)) )),2.))/umag)**2
  !endwhere
  thetavstar= vkar*(thetav-sthetav)/integralh
end do

return
end subroutine dyerhicks

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate shortwave radiation coefficents (modified to include 2nd wall)

subroutine getswcoeff(ifull,sg,wallpsi,roadpsi,ifn,rdsnalpha,rdsndelta)

implicit none

integer, intent(in) :: ifull
integer k
real, dimension(ifull), intent(in) :: rdsnalpha,rdsndelta
real, dimension(ifull), intent(out) :: wallpsi,roadpsi
real, dimension(ifull) :: thetazero,walles,wallws,roads,ta,tc,xa,ya,roadnetalpha
real, dimension(ifull) :: nwalles,nwallws,nroads
type(tdata), dimension(ifull), intent(in) :: ifn
type(trad), dimension(ifull), intent(out) :: sg

wallpsi=0.5*(ifn%hwratio+1.-sqrt(ifn%hwratio*ifn%hwratio+1.))/ifn%hwratio
roadpsi=sqrt(ifn%hwratio*ifn%hwratio+1.)-ifn%hwratio

! integrate through 180 instead of 360
where (ifn%vangle.ge.0.5*pi)
  walles=0.
  wallws=1./ifn%hwratio
  roads=0.
elsewhere
  ta=tan(ifn%vangle)
  thetazero=asin(1./max(ifn%hwratio*ta,1.))
  tc=2.*(1.-cos(thetazero))
  xa=min(max(ifn%hangle-thetazero,0.),pi)-max(ifn%hangle-pi+thetazero,0.)-min(ifn%hangle+thetazero,0.)
  ya=cos(max(min(0.,ifn%hangle),ifn%hangle-pi))-cos(max(min(thetazero,ifn%hangle),ifn%hangle-pi)) &
    +cos(min(0.,-ifn%hangle))-cos(min(thetazero,-ifn%hangle)) &
    +cos(max(0.,pi-ifn%hangle))-cos(max(thetazero,pi-ifn%hangle))
  ! note that these terms now include the azimuth angle
  walles=(xa/ifn%hwratio+ta*ya)/pi
  wallws=((pi-2.*thetazero-xa)/ifn%hwratio+ta*(tc-ya))/pi
  roads=(2.*thetazero-ifn%hwratio*ta*tc)/pi
end where

! Calculate short wave reflections to nrefl order
roadnetalpha=rdsndelta*rdsnalpha+(1.-rdsndelta)*((1.-ifn%sigmaveg)*ifn%roadalpha+ifn%sigmaveg*ifn%vegalpha)
sg%walle=walles
sg%wallw=wallws
sg%road=roads
do k=1,nrefl
  nwalles=roadnetalpha*wallpsi*roads+ifn%wallalpha*(1.-2.*wallpsi)*wallws
  nwallws=roadnetalpha*wallpsi*roads+ifn%wallalpha*(1.-2.*wallpsi)*walles
  nroads=ifn%wallalpha*(1.-roadpsi)*0.5*(walles+wallws)
  walles=nwalles
  wallws=nwallws
  roads=nroads
  sg%walle=sg%walle+walles
  sg%wallw=sg%wallw+wallws
  sg%road=sg%road+roads
end do
sg%roof=1.
sg%rfsn=1.
sg%rdsn=sg%road
sg%veg=sg%road

return
end subroutine getswcoeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solve for roof snow temperature

subroutine solverfsn(cn,evct,rg,fg,eg,garfsn,rfsnmelt,rfsntemp,iroof,dg,sg,atm,acond,ifn,ddt)

implicit none

integer, intent(in) :: cn
real, intent(in) :: ddt
real, dimension(cn), intent(out) :: evct,rfsnmelt,garfsn
real, dimension(cn), intent(in) :: rfsntemp
real, dimension(cn) :: cd,rfsnqsat,lzotdum,lzosnow,sndepth,snlambda,ldratio,we,ww,z_on_l
type(trad), dimension(cn), intent(inout) :: rg,fg,eg,acond
type(trad), dimension(cn), intent(in) :: sg
type(tatm), dimension(cn), intent(in) :: atm
type(tdiag), dimension(cn), intent(in) :: dg
type(tsurf), dimension(cn), intent(in) :: iroof
type(tdata), dimension(cn), intent(in) :: ifn

! snow conductance
sndepth=iroof%snow*waterden/iroof%den
snlambda=icelambda*(iroof%den/waterden)**1.88
ldratio=0.5*(sndepth/snlambda+ifn%roofdepth(1)/ifn%rooflambda(1))

! Update roof snow energy budget
lzosnow=log(dg%rfdzmin/zosnow)
call getqsat(cn,rfsnqsat,rfsntemp,dg%sigr)
we=rfsntemp*(1.+0.61*rfsnqsat)
ww=dg%tempr*(1.+0.61*dg%mixrr)
call getinvres(cn,acond%rfsn,cd,z_on_l,lzotdum,lzosnow,dg%rfdzmin,we,ww,atm%umag,1)
acond%rfsn=acond%rfsn/atm%umag
rfsnmelt=dg%rfsndelta*max(0.,rfsntemp-273.16)/(icecp*iroof%den*lf*ddt) 
rg%rfsn=snowemiss*(atm%rg-sbconst*rfsntemp**4)
fg%rfsn=aircp*atm%rho*(rfsntemp-dg%tempr)*acond%rfsn*atm%umag
eg%rfsn=ls*min(atm%rho*dg%rfsndelta*max(0.,rfsnqsat-dg%mixrr)*acond%rfsn*atm%umag,iroof%snow/ddt+atm%snd-rfsnmelt) ! MJT suggestion for max eg
garfsn=(rfsntemp-iroof%temp(1))/ldratio
evct=sg%rfsn+rg%rfsn-fg%rfsn-eg%rfsn-garfsn

return
end subroutine solverfsn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solve for road snow temperature (includes vegetation canopy temperature and canyon temperature)

subroutine solverdsn(cn,evct,rg,fg,fgtop,eg,gardsn,rdsnmelt,rdsntemp,iroad,iwalle,iwallw,iveg,dg,sg,atm &
                    ,ddt,acond,wallpsi,roadpsi,ifn,ipg)

implicit none

integer, intent(in) :: cn
integer k,ic
real, intent(in) :: ddt
real, dimension(cn), intent(out) :: evct,fgtop,gardsn,rdsnmelt
real, dimension(cn), intent(in) :: rdsntemp,wallpsi,roadpsi
real, dimension(cn) :: ctmax,ctmin,topinvres
real, dimension(cn) :: trafficout,ldratio,sndepth,snlambda
real, dimension(cn) :: oldval,newval,evctveg,evctx
type(trad), dimension(cn), intent(inout) :: rg,fg,eg,acond
type(trad), dimension(cn), intent(in) :: sg
type(tatm), dimension(cn), intent(in) :: atm
type(tdiag), dimension(cn), intent(inout) :: dg
type(tsurf), dimension(cn), intent(in) :: iroad
type(twall), dimension(cn), intent(in) :: iwalle,iwallw
type(tvege), dimension(cn), intent(in) :: iveg
type(tdata), dimension(cn), intent(in) :: ifn
type(tprog), dimension(cn), intent(inout) :: ipg

! traffic sensible heat flux
call gettraffic(cn,trafficout,ifn%trafficfg,ifn%ctime)
trafficout=trafficout/(1.-ifn%sigmabld)

! snow conductance
sndepth=iroad%snow*waterden/iroad%den
snlambda=icelambda*(iroad%den/waterden)**1.88
ldratio=0.5*(sndepth/snlambda+ifn%roaddepth(1)/ifn%roadlambda(1))

! solve for vegetation canopy and canyon temperature
! (we used to employ a multi-variable Newton-Raphson method here.  However, it is difficult to handle the
!  dependence of the wind speed at the canyon top (i.e., affecting aerodynamical resistances)
!  when the wind speed at the canyon top depends on the stability functions.  This multiply
!  nested version is slower, but more robust).
ctmax=max(dg%tempc,iroad%temp(1),iwalle%temp(1),iwallw%temp(1),rdsntemp)+10. ! max temp
ctmin=min(dg%tempc,iroad%temp(1),iwalle%temp(1),iwallw%temp(1),rdsntemp)-10. ! min temp
if (any(ifn%sigmaveg.gt.0.)) then ! in-canyon vegetation
  ipg%vegtemp=0.6*ctmax+0.4*ctmin
  call solverdvg(cn,evctveg,rg,fg,eg,gardsn,rdsnmelt,topinvres,rdsntemp,iroad,iwalle,iwallw,iveg,dg,sg,atm &
                    ,ddt,acond,wallpsi,roadpsi,ifn,ipg,trafficout,ldratio)
  oldval=ipg%vegtemp
  ipg%vegtemp=0.4*ctmax+0.6*ctmin
  do k=1,nfgits ! Sectant
    evctx=evctveg
    call solverdvg(cn,evctveg,rg,fg,eg,gardsn,rdsnmelt,topinvres,rdsntemp,iroad,iwalle,iwallw,iveg,dg,sg,atm &
                    ,ddt,acond,wallpsi,roadpsi,ifn,ipg,trafficout,ldratio)
    evctx=evctveg-evctx
    if (all(abs(evctx).le.tol)) exit
    where (abs(evctx).gt.tol)  
      newval=ipg%vegtemp-alpha*evctveg*(ipg%vegtemp-oldval)/evctx
      oldval=ipg%vegtemp
      ipg%vegtemp=newval
    end where
    ipg%vegtemp=min(max(ipg%vegtemp,ctmin),ctmax)
  end do
else ! no in-canyon vegetation
  call solverdvg(cn,evctveg,rg,fg,eg,gardsn,rdsnmelt,topinvres,rdsntemp,iroad,iwalle,iwallw,iveg,dg,sg,atm &
                    ,ddt,acond,wallpsi,roadpsi,ifn,ipg,trafficout,ldratio)
end if

! ---------------------------------------------------------------    

! balance sensible heat flux
fgtop=(dg%rdsndelta*fg%rdsn+(1.-dg%rdsndelta)*((1.-ifn%sigmaveg)*fg%road &
      +ifn%sigmaveg*fg%veg)+ifn%hwratio*(fg%walle+fg%wallw)+trafficout+dg%accool)
dg%canyontemp=fgtop/(aircp*atm%rho*topinvres)+dg%tempc

! road snow energy balance
evct=sg%rdsn+rg%rdsn-fg%rdsn-eg%rdsn-gardsn

return
end subroutine solverdsn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solve for vegetation canopy temperature

subroutine solverdvg(cn,evct,rg,fg,eg,gardsn,rdsnmelt,topinvres,rdsntemp,iroad,iwalle,iwallw,iveg,dg,sg,atm &
                    ,ddt,acond,wallpsi,roadpsi,ifn,ipg,trafficout,ldratio)

implicit none

integer, intent(in) :: cn
integer k
real, intent(in) :: ddt
real, dimension(cn), intent(out) :: evct,gardsn,rdsnmelt,topinvres
real, dimension(cn), intent(in) :: rdsntemp,wallpsi,roadpsi
real, dimension(cn), intent(in) :: ldratio,trafficout
real, dimension(cn) :: oldval,evctveg,evctx,newval,ctmax,ctmin
type(trad), dimension(cn), intent(inout) :: rg,fg,eg,acond
type(trad), dimension(cn), intent(in) :: sg
type(tatm), dimension(cn), intent(in) :: atm
type(tdiag), dimension(cn), intent(inout) :: dg
type(tsurf), dimension(cn), intent(in) :: iroad
type(twall), dimension(cn), intent(in) :: iwalle,iwallw
type(tvege), dimension(cn), intent(in) :: iveg
type(tdata), dimension(cn), intent(in) :: ifn
type(tprog), dimension(cn), intent(inout) :: ipg

! update canyon temperature
ctmax=max(dg%tempc,iroad%temp(1),iwalle%temp(1),iwallw%temp(1),rdsntemp)+10. ! max temp
ctmin=min(dg%tempc,iroad%temp(1),iwalle%temp(1),iwallw%temp(1),rdsntemp)-10. ! min temp
dg%canyontemp=0.4*ctmax+0.6*ctmin
call solvecanyon(cn,evctveg,rg,fg,eg,gardsn,rdsnmelt,topinvres,rdsntemp,iroad,iwalle, &
     iwallw,iveg,dg,sg,atm,ddt,acond,wallpsi,roadpsi,ifn,ipg,trafficout,ldratio)
oldval=dg%canyontemp
dg%canyontemp=0.6*ctmax+0.4*ctmin
do k=1,nfgits
  evctx=evctveg
  call solvecanyon(cn,evctveg,rg,fg,eg,gardsn,rdsnmelt,topinvres,rdsntemp,iroad,iwalle, &
       iwallw,iveg,dg,sg,atm,ddt,acond,wallpsi,roadpsi,ifn,ipg,trafficout,ldratio)
  evctx=evctveg-evctx
  if (all(abs(evctx).le.tol)) exit
  where (abs(evctx).gt.tol) 
    newval=dg%canyontemp-alpha*evctveg*(dg%canyontemp-oldval)/evctx
    oldval=dg%canyontemp
    dg%canyontemp=newval
  end where
  dg%canyontemp=min(max(dg%canyontemp,ctmin),ctmax)
end do

! vegetation energy budget
evct=sg%veg+rg%veg-fg%veg-eg%veg

return
end subroutine solverdvg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solve for canyon temperature

subroutine solvecanyon(cn,evct,rg,fg,eg,gardsn,rdsnmelt,topinvres,rdsntemp,iroad,iwalle,iwallw,iveg &
                    ,dg,sg,atm,ddt,acond,wallpsi,roadpsi,ifn,ipg,trafficout,ldratio)

implicit none

integer, intent(in) :: cn
integer ic
real, intent(in) :: ddt
real, dimension(cn), intent(out) :: gardsn,rdsnmelt,topinvres
real, dimension(cn), intent(out) :: evct
real, dimension(cn), intent(in) :: rdsntemp,wallpsi,roadpsi,trafficout,ldratio
real, dimension(cn) :: roadqsat,rdsnqsat,fgtop
real, dimension(cn) :: vegqsat,res,f1,f2,f3,f4,ff,qsat
real, dimension(cn) :: dumroaddelta,dumvegdelta,we,ww
real, dimension(cn) :: newcanyonmix,oldcanyonmix,evctmx,evctmxdif
real, dimension(cn) :: z_on_l
type(trad), dimension(cn), intent(inout) :: rg,fg,eg,acond
type(trad), dimension(cn), intent(in) :: sg
type(tatm), dimension(cn), intent(in) :: atm
type(tdiag), dimension(cn), intent(inout) :: dg
type(tsurf), dimension(cn), intent(in) :: iroad
type(twall), dimension(cn), intent(in) :: iwalle,iwallw
type(tvege), dimension(cn), intent(in) :: iveg
type(tdata), dimension(cn), intent(in) :: ifn
type(tprog), dimension(cn), intent(inout) :: ipg

! transpiration terms (developed by Eva in CCAM sflux.f and CSIRO9)
ff=1.1*sg%veg/(vegrlai*150.)
f1=(1.+ff)/(ff+vegrsmin*vegrlai/5000.)
f2=max(0.5*(sfc-swilt)/max(iveg%moist-swilt,0.01*(sfc-swilt)),1.)
f4=max(1.-.0016*(298.-dg%canyontemp)**2,0.05)
call getqsat(cn,qsat,dg%canyontemp,dg%sigd)

! estimate mixing ratio at canyon surfaces
call getqsat(cn,roadqsat,iroad%temp(1),dg%sigd) ! evaluate using pressure at displacement height
call getqsat(cn,vegqsat,ipg%vegtemp,dg%sigd)    ! evaluate using pressure at displacement height
call getqsat(cn,rdsnqsat,rdsntemp,dg%sigd)

! set-up initial guess for canyon mixing ratio
dg%canyonmix=0.
oldcanyonmix=dg%canyonmix
we=dg%canyontemp*(1.+0.61*dg%canyonmix)
ww=dg%tempc*(1.+0.61*dg%mixrc)
call getinvres(cn,topinvres,ipg%cduv,z_on_l,ipg%lzoh,ipg%lzom,ipg%cndzmin,we,ww,atm%umag,3)
call gettopu(cn,dg%topu,atm%umag,ipg%cndzmin,ifn%bldheight*(1.-refheight),ifn%hwratio,z_on_l,ipg%cduv)
where (roadqsat.lt.dg%canyonmix)
  dumroaddelta=1.
elsewhere
  dumroaddelta=dg%roaddelta
endwhere
where (vegqsat.lt.dg%canyonmix)
  dumvegdelta=1.
elsewhere
  dumvegdelta=dg%vegdelta
endwhere
! remaining transpiration terms
f3=max(1.-.00025*(qsat-dg%canyonmix*dg%sigd/0.622),0.05)
res=max(30.,vegrsmin*f1*f2/(f3*f4))
evctmx=(dg%rdsndelta*rdsnqsat*ls/lv*acond%rdsn*dg%topu &
       +(1.-dg%rdsndelta)*((1.-fn%sigmaveg)*dumroaddelta*roadqsat*acond%road*dg%topu &
       +ifn%sigmaveg*vegqsat*(dumvegdelta*acond%veg*dg%topu &
       +(1.-dumvegdelta)/(1./(acond%veg*dg%topu)+res)))+dg%mixrc*topinvres) &
       -dg%canyonmix*(dg%rdsndelta*ls/lv*acond%rdsn*dg%topu+(1.-dg%rdsndelta)*((1.-fn%sigmaveg)*dumroaddelta*acond%road*dg%topu &
       +ifn%sigmaveg*(dumvegdelta*acond%veg*dg%topu+(1.-dumvegdelta)/(1./(acond%veg*dg%topu)+res)))+topinvres)
dg%canyonmix=dg%mixrc

! Sectant solver for canyonmix
do ic=1,negits
  evctmxdif=evctmx ! store to calculate derivative
  we=dg%canyontemp*(1.+0.61*dg%canyonmix)
  ! solve for aerodynamical resistance between canyon and atmosphere (neglect molecular diffusion)
  call getinvres(cn,topinvres,ipg%cduv,z_on_l,ipg%lzoh,ipg%lzom,ipg%cndzmin,we,ww,atm%umag,3)
  call gettopu(cn,dg%topu,atm%umag,ipg%cndzmin,ifn%bldheight*(1.-refheight),ifn%hwratio,z_on_l,ipg%cduv)
  ! correction for dew
  where (roadqsat.lt.dg%canyonmix)
    dumroaddelta=1.
  elsewhere
    dumroaddelta=dg%roaddelta
  endwhere
  where (vegqsat.lt.dg%canyonmix)
    dumvegdelta=1.
  elsewhere
    dumvegdelta=dg%vegdelta
  endwhere
  ! remaining transpiration terms
  f3=max(1.-.00025*(qsat-dg%canyonmix*dg%sigd/0.622),0.05)
  res=max(30.,vegrsmin*f1*f2/(f3*f4))
  ! balance canyon latent heat budget
  evctmx=(dg%rdsndelta*rdsnqsat*ls/lv*acond%rdsn*dg%topu &
         +(1.-dg%rdsndelta)*((1.-fn%sigmaveg)*dumroaddelta*roadqsat*acond%road*dg%topu &
         +ifn%sigmaveg*vegqsat*(dumvegdelta*acond%veg*dg%topu &
         +(1.-dumvegdelta)/(1./(acond%veg*dg%topu)+res)))+dg%mixrc*topinvres) &
         -dg%canyonmix*(dg%rdsndelta*ls/lv*acond%rdsn*dg%topu+(1.-dg%rdsndelta)*((1.-fn%sigmaveg)*dumroaddelta*acond%road*dg%topu &
         +ifn%sigmaveg*(dumvegdelta*acond%veg*dg%topu+(1.-dumvegdelta)/(1./(acond%veg*dg%topu)+res)))+topinvres)
  evctmxdif=evctmx-evctmxdif
  if (all(abs(evctmxdif).eq.tolqg)) exit
  where (abs(evctmxdif).gt.tolqg)
    newcanyonmix=dg%canyonmix-alpha*evctmx*(dg%canyonmix-oldcanyonmix)/evctmxdif
    oldcanyonmix=dg%canyonmix
    dg%canyonmix=newcanyonmix
  end where
  dg%canyonmix=min(max(dg%canyonmix,0.),1.)
end do

! solve for canyon sensible heat flux ----------------------------------
if (resmeth.eq.0) then
  !cw=pg%cduv*atm%umag ! diagnose canyonw (from Masson 2000)
  acond%road=(11.8+4.2*sqrt(acond%road**2+ipg%cduv*atm%umag**2))/(aircp*atm%rho) ! From Rowley, et al (1930)
  acond%walle=acond%road
  acond%wallw=acond%road
  acond%rdsn=acond%road
  acond%veg=acond%road
end if
fg%walle=aircp*atm%rho*(iwalle%temp(1)-dg%canyontemp)*acond%walle*dg%topu
fg%wallw=aircp*atm%rho*(iwallw%temp(1)-dg%canyontemp)*acond%wallw*dg%topu
fg%road=aircp*atm%rho*(iroad%temp(1)-dg%canyontemp)*acond%road*dg%topu
fg%veg=aircp*atm%rho*(ipg%vegtemp-dg%canyontemp)*acond%veg*dg%topu
where (dg%rdsndelta.gt.0.)
  fg%rdsn=aircp*atm%rho*(rdsntemp-dg%canyontemp)*acond%rdsn*dg%topu
elsewhere
  fg%rdsn=0.
end where
fgtop=aircp*atm%rho*(dg%canyontemp-dg%tempc)*topinvres
! ---------------------------------------------------------------   

! calculate longwave radiation
dg%netrad=dg%rdsndelta*snowemiss*rdsntemp**4+(1.-dg%rdsndelta)*((1.-ifn%sigmaveg)*ifn%roademiss*iroad%temp(1)**4 &
                  +ifn%sigmaveg*ifn%vegemiss*ipg%vegtemp**4)    
rg%walle=ifn%wallemiss*(atm%rg*dg%cwa+sbconst*iwalle%temp(1)**4*(-1.+ifn%wallemiss*dg%cwe) & 
                  +sbconst*iwallw%temp(1)**4*ifn%wallemiss*dg%cww+sbconst*dg%netrad*dg%cwr)
rg%wallw=ifn%wallemiss*(atm%rg*dg%cwa+sbconst*iwallw%temp(1)**4*(-1.+ifn%wallemiss*dg%cwe) &
                  +sbconst*iwalle%temp(1)**4*ifn%wallemiss*dg%cww+sbconst*dg%netrad*dg%cwr)
rg%road=ifn%roademiss*(atm%rg*dg%cra+sbconst*(-iroad%temp(1)**4+dg%netrad*dg%crr) &
                  +sbconst*ifn%wallemiss*(iwalle%temp(1)**4+iwallw%temp(1)**4)*dg%crw)
rg%veg=ifn%vegemiss*(atm%rg*dg%cra+sbconst*(-ipg%vegtemp**4+dg%netrad*dg%crr) &
                  +sbconst*ifn%wallemiss*(iwalle%temp(1)**4+iwallw%temp(1)**4)*dg%crw)                  
where (dg%rdsndelta.gt.0.)
  rg%rdsn=snowemiss*(atm%rg*dg%cra+sbconst*(-rdsntemp**4+dg%netrad*dg%crr) &
                    +sbconst*ifn%wallemiss*(iwalle%temp(1)**4+iwallw%temp(1)**4)*dg%crw)
elsewhere
  rg%rdsn=0.
end where
dg%canyonrgout=atm%rg*(2.*ifn%hwratio*wallpsi*(1.-ifn%wallemiss)*dg%cwa+roadpsi*(1.-dg%netemiss)*dg%cwr) &
               +sbconst*ifn%wallemiss*iwalle%temp(1)**4*(ifn%hwratio*wallpsi*(1.+(1.-ifn%wallemiss)*(dg%cwe+dg%cww)) &
               +roadpsi*(1.-dg%netemiss)*dg%crw) &
               +sbconst*ifn%wallemiss*iwallw%temp(1)**4*(ifn%hwratio*wallpsi*(1.+(1.-ifn%wallemiss)*(dg%cwe+dg%cww)) &
               +roadpsi*(1.-dg%netemiss)*dg%crw) &
               +sbconst*dg%netrad*(2.*ifn%hwratio*wallpsi*(1.-ifn%wallemiss)*dg%cwr+roadpsi*(1.+(1.-dg%netemiss)*dg%crr))

! estimate snow melt
rdsnmelt=dg%rdsndelta*max(0.,rdsntemp-273.16)/(icecp*iroad%den*lf*ddt)

! calculate transpiration and evaporation of in-canyon vegetation
dg%tran=lv*min(max((1.-dumvegdelta)*atm%rho*(vegqsat-dg%canyonmix)/(1./(acond%veg*dg%topu)+res),0.), &
               max((iveg%moist-swilt)*dg%totdepth*waterden/(dg%c1*ddt),0.))
dg%evap=lv*min(dumvegdelta*atm%rho*(vegqsat-dg%canyonmix)*acond%veg*dg%topu,iveg%water/ddt+atm%rnd)
eg%veg=dg%evap+dg%tran

! calculate canyon latent heat fluxes
eg%road=lv*min(atm%rho*dumroaddelta*(roadqsat-dg%canyonmix)*acond%road*dg%topu &
             ,iroad%water/ddt+atm%rnd+(1.-fn%sigmaveg)*rdsnmelt)
where (dg%rdsndelta.gt.0.)
  eg%rdsn=ls*min(atm%rho*dg%rdsndelta*max(0.,rdsnqsat-dg%canyonmix)*acond%rdsn*dg%topu &
                ,iroad%snow/ddt+atm%snd-rdsnmelt)
  gardsn=(rdsntemp-iroad%temp(1))/ldratio ! use road temperature to represent canyon bottom surface temperature
                                          ! (i.e., we have ommited soil under vegetation temperature)
elsewhere
  eg%rdsn=0.
  gardsn=0.
end where
eg%walle=0.
eg%wallw=0.

! balance sensible heat flux
evct=fgtop-(dg%rdsndelta*fg%rdsn+(1.-dg%rdsndelta)*((1.-ifn%sigmaveg)*fg%road &
     +ifn%sigmaveg*fg%veg)+ifn%hwratio*(fg%walle+fg%wallw)+trafficout+dg%accool)

return
end subroutine solvecanyon


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define traffic flux weights during the diurnal cycle

subroutine gettraffic(cn,trafficout,trafficfg,ctime)

implicit none

integer, intent(in) :: cn
integer, dimension(cn) :: ip
real, dimension(cn), intent(out) :: trafficout
real, dimension(cn), intent(in) :: trafficfg,ctime
real, dimension(cn) :: rp
! traffic diurnal cycle weights approximated from Coutts et al (2007)
real, dimension(25), parameter :: trafficcycle = (/ 0.1, 0.1, 0.1, 0.1, 0.5, 1.0, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, &
                                                    1.5, 1.5, 1.5, 1.5, 1.5, 1.4, 1.2,  1., 0.8, 0.6, 0.4, 0.2, &
                                                    0.1 /) 

ip=int(24.*ctime)
rp=24.*ctime-real(ip)
where (ip.lt.1) ip=ip+24
trafficout=trafficfg*((1.-rp)*trafficcycle(ip)+rp*trafficcycle(ip+1))

return
end subroutine gettraffic

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate in-canyon wind speed for walls and road

subroutine getincanwind(ueast,uwest,ufloor,wdirin,h,w,z0)

implicit none

real, dimension(ufull), intent(out) :: ueast,uwest,ufloor
real, dimension(ufull), intent(in) :: wdirin,h,w,z0
real, dimension(ufull) :: a,b,wsuma,wsumb,fsum
real, dimension(ufull) :: theta1,wdir

where (wdirin.ge.0.)
  wdir=wdirin
elsewhere
  wdir=wdirin+pi
endwhere

theta1=acos(min(w/(3.*h),1.))
wsuma=0.
wsumb=0.
fsum=0.

! integrate jet on road, venting side
a=0.
b=max(0.,wdir-pi+theta1)
call integratewind(wsuma,wsumb,fsum,a,b,h,w,wdir,z0,0)

! integrate jet on wall, venting side
a=max(0.,wdir-pi+theta1)
b=max(0.,wdir-theta1)
call integratewind(wsuma,wsumb,fsum,a,b,h,w,wdir,z0,1)

! integrate jet on road, venting side
a=max(0.,wdir-theta1)
b=wdir
call integratewind(wsuma,wsumb,fsum,a,b,h,w,wdir,z0,0)

! integrate jet on road, recirculation side (A)
a=wdir
b=min(pi,wdir+theta1)
call integratewind(wsumb,wsuma,fsum,a,b,h,w,wdir,z0,0)

! integrate jet on wall, recirculation side
a=min(pi,wdir+theta1)
b=min(pi,wdir+pi-theta1)
call integratewind(wsumb,wsuma,fsum,a,b,h,w,wdir,z0,1)

! integrate jet on road, recirculation side (B)
a=min(pi,wdir+pi-theta1)
b=pi
call integratewind(wsumb,wsuma,fsum,a,b,h,w,wdir,z0,0)

! Determine orientation
where (wdirin.ge.0.)
  ueast=0.5*wsuma
  uwest=0.5*wsumb
elsewhere
  ueast=0.5*wsumb
  uwest=0.5*wsuma
end where
ufloor=0.5*fsum

return
end subroutine getincanwind

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! integrate winds

subroutine integratewind(wsuma,wsumb,fsum,a,b,h,w,wdir,z0,mode)

implicit none

integer, intent(in) :: mode
integer n
integer, parameter :: ntot=45
real, dimension(ufull), intent(in) :: a,b,h,w,wdir,z0
real, dimension(ufull), intent(inout) :: wsuma,wsumb,fsum
real, dimension(ufull) :: theta,dtheta,st,nw
real, dimension(ufull) :: duf,dur,duv

dtheta=(b-a)/real(ntot)
select case(mode)
  case(0) ! jet on road
    do n=1,ntot
      theta=dtheta*(real(n)-0.5)+a
      st=abs(sin(theta-wdir))
      nw=w/max(st,1.E-9)
      call winda(duf,dur,duv,h,nw,z0)
      wsuma=wsuma+dur*st*dtheta
      wsumb=wsumb+duv*st*dtheta
      fsum=fsum+duf*st*dtheta
    end do
  case(1) ! jet on wall
    do n=1,ntot
      theta=dtheta*(real(n)-0.5)+a
      st=abs(sin(theta-wdir))
      nw=w/max(st,1.E-9)
      call windb(duf,dur,duv,h,nw,z0)
      wsuma=wsuma+dur*st*dtheta
      wsumb=wsumb+duv*st*dtheta
      fsum=fsum+duf*st*dtheta
    end do
  case DEFAULT
    write(6,*) "ERROR: Unknown mode ",mode
    stop
end select

return
end subroutine integratewind

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate wind speed, jet on road

subroutine winda(uf,ur,uv,h,w,z0)

implicit none

real, dimension(ufull), intent(out) :: uf,ur,uv
real, dimension(ufull), intent(in) :: h,w,z0
real, dimension(ufull) :: a,u0,cuven,zolog

a=0.15*max(1.,3.*h/(2.*w))
u0=exp(-0.9*sqrt(13./4.))

zolog=log(max(h,z0+0.1)/z0)
cuven=log(max(0.1*h,z0+0.1)/z0)/log(max(h,z0+0.1)/z0)
cuven=max(cuven*(1.-3.*h/w),(u0/a)*(h/w)*(1.-exp(-a*(w/h-3.))))
uf=(u0/a)*(h/w)*(1.-exp(-3.*a))+cuven
!uf=(u0/a)*(h/w)*(2.-exp(-a*3.)-exp(-a*(w/h-3.)))
ur=(u0/a)*exp(-a*3.)*(1.-exp(-a))
! MJT suggestion
cuven=1.-1./zolog
uv=(u0/a)*exp(-a*(w/h-3.))*(1.-exp(-a))
uv=max(cuven,uv)
!uv=(u0/a)*exp(-a*(w/h-3.))*(1.-exp(-a))

return
end subroutine winda

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate wind speed, jet on wall

subroutine windb(uf,ur,uv,h,win,z0)

implicit none

real, dimension(ufull), intent(out) :: uf,ur,uv
real, dimension(ufull), intent(in) :: h,win,z0
real, dimension(ufull) :: a,dh,u0,w
real, dimension(ufull) :: zolog,cuven

w=min(win,1.5*h)

a=0.15*max(1.,3.*h/(2.*w))
dh=max(2.*w/3.-h,0.)
u0=exp(-0.9*sqrt(13./4.)*dh/h)

zolog=log(max(h,z0+0.1)/z0)
! MJT suggestion (cuven is multipled by dh to avoid divide by zero)
cuven=h-(h-dh)*log(max(h-dh,z0+0.1)/z0)/zolog-dh/zolog
! MJT cuven is back to the correct units of m/s
cuven=max(cuven/h,(u0/a)*(1.-exp(-a*dh/h)))

uf=(u0/a)*(h/w)*exp(-a*(1.-dh/h))*(1.-exp(-a*w/h))
ur=(u0/a)*exp(-a*(1.-dh/h+w/h))*(1.-exp(-a))
uv=(u0/a)*(1.-exp(-a*(1.-dh/h)))+cuven
!uv=(u0/a)*(2.-exp(-a*(1.-dh/h))-exp(-a*dh/h))

return
end subroutine windb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate wind speed at canyon top
subroutine gettopu(cn,topu,umag,zmin,bldheight,hwratio,z_on_l,cd)
      
implicit none

integer, intent(in) :: cn
real, dimension(cn), intent(in) :: umag,zmin
real, dimension(cn), intent(in) :: bldheight,hwratio
real, dimension(cn), intent(in) :: z_on_l,cd
real, dimension(cn), intent(out) :: topu
real, dimension(cn) :: z0_on_l
real, dimension(cn) :: pm0,pm1,integralm
real, dimension(cn) :: ustar,neutral

ustar=sqrt(cd)*umag

z0_on_l=min(bldheight,zmin)*z_on_l/zmin ! calculate at canyon top
z0_on_l=min(z0_on_l,10.)
neutral = log(zmin/min(bldheight,zmin))
where (z_on_l.lt.0.)
  pm0     = (1.-16.*z0_on_l)**(-0.25)
  pm1     = (1.-16.*z_on_l)**(-0.25)
  integralm = neutral-2.*log((1.+1./pm1)/(1.+1./pm0)) &
                -log((1.+1./pm1**2)/(1.+1./pm0**2)) &
                +2.*(atan(1./pm1)-atan(1./pm0))
elsewhere
  !-------Beljaars and Holtslag (1991) heat function
  pm0 = -(a_1*z0_on_l+b_1*(z0_on_l-(c_1/d_1))*exp(-d_1*z0_on_l)+b_1*c_1/d_1)
  pm1  = -(a_1*z_on_l+b_1*(z_on_l-(c_1/d_1))*exp(-d_1*z_on_l)+b_1*c_1/d_1)
  integralm = neutral-(pm1-pm0)
endwhere
where (bldheight.lt.zmin)
  topu=(2./pi)*(umag-ustar*integralm/vkar)
elsewhere ! within canyon
  topu=(2./pi)*umag*exp(0.5*hwratio*(1.-zmin/bldheight))
end where

topu=max(topu,0.1)

return
end subroutine gettopu

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate screen diagnostics
subroutine scrncalc(stemp,temp,smixr,mixr,umag,dg,rdsntemp,moist,zonet)
      
implicit none

integer ic
real, dimension(ufull), intent(in) :: stemp,temp,umag
real, dimension(ufull), intent(in) :: smixr,mixr,rdsntemp,moist,zonet
real, dimension(ufull) :: af,aft,ri,fm,fh,root
real, dimension(ufull) :: denma,denha,cd,thetav,sthetav
real, dimension(ufull) :: thetavstar,z_on_l,z0_on_l
real, dimension(ufull) :: pm0,ph0,pm1,ph1,integralm,integralh
real, dimension(ufull) :: ustar,qstar,z10_on_l
real, dimension(ufull) :: neutral,neutral10,pm10
real, dimension(ufull) :: integralm10
real, dimension(ufull) :: tstar,lna
real, dimension(ufull) :: utop,ttop,qtop,wf,tsurf,qsurf,n
type(tdiag), dimension(ufull), intent(in) :: dg
real, parameter    ::  z0     = 1.5
real, parameter    ::  z10    = 10.

thetav=temp*(1.+0.61*mixr)
sthetav=stemp*(1.+0.61*smixr)

! Roughness length for heat
lna=pg%lzoh-pg%lzom

! use Louis as first guess for Dyer and Hicks scheme
af=vkar*vkar/(pg%lzom*pg%lzom)
aft=vkar*vkar/(pg%lzom*pg%lzoh)
! umag is now constrained to be above umin
ri=min(grav*pg%cndzmin*(1.-sthetav/thetav)/umag**2,rimax)
where (ri>0.)
  fm=1./(1.+bprm*ri)**2
  fh=fm
elsewhere
  root=sqrt(-ri*exp(pg%lzom))
  denma=1.+cms*2.*bprm*af*root
  denha=1.+chs*2.*bprm*aft*exp(0.5*lna)*root
  fm=1.-2.*bprm *ri/denma
  fh=1.-2.*bprm*ri/denha
end where
cd=af*fm

! Dyer and Hicks approach 
thetavstar=aft*fh*(thetav-sthetav)/sqrt(cd)
call dyerhicks(ufull,integralh,z_on_l,cd,thetavstar,thetav,sthetav,umag,pg%cndzmin,pg%lzom,lna)
ustar=sqrt(cd)*umag
tstar=vkar*(temp-stemp)/integralh
qstar=vkar*(mixr-smixr)/integralh

select case(scrnmeth)
  case(0) ! estimate screen diagnostics (slab at displacement height approach)
    z0_on_l   = z0*z_on_l/pg%cndzmin
    z10_on_l  = z10*z_on_l/pg%cndzmin
    z0_on_l=min(z0_on_l,10.)
    z10_on_l=min(z10_on_l,10.)
    neutral = log(pg%cndzmin/z0)
    neutral10 = log(pg%cndzmin/z10)
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
      pm0 = -(a_1*z0_on_l+b_1*(z0_on_l-(c_1/d_1))*exp(-d_1*z0_on_l)+b_1*c_1/d_1)
      pm10 = -(a_1*z10_on_l+b_1*(z10_on_l-(c_1/d_1))*exp(-d_1*z10_on_l)+b_1*c_1/d_1)
      pm1  = -(a_1*z_on_l+b_1*(z_on_l-(c_1/d_1))*exp(-d_1*z_on_l)+b_1*c_1/d_1)
      integralh = neutral-(ph1-ph0)
      integralm = neutral-(pm1-pm0)
      integralm10 = neutral10-(pm1-pm10)
    endwhere
    pg%tscrn       = temp-tstar/vkar*integralh
    pg%qscrn       = mixr-qstar/vkar*integralh
    pg%uscrn=max(umag-ustar/vkar*integralm,0.)
    pg%u10=max(umag-ustar/vkar*integralm10,0.)
    
  case(1) ! estimate screen diagnostics (two step canopy approach)
    z0_on_l   = fn%bldheight*(1.-refheight)*z_on_l/pg%cndzmin ! calculate at canyon top
    z10_on_l  = max(z10-fn%bldheight*refheight,1.)*z_on_l/pg%cndzmin
    z0_on_l=min(z0_on_l,10.)
    z10_on_l=min(z10_on_l,10.)
    neutral = log(pg%cndzmin/(fn%bldheight*(1.-refheight)))
    neutral10 = log(pg%cndzmin/max(z10-fn%bldheight*refheight,1.))
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
      pm0 = -(a_1*z0_on_l+b_1*(z0_on_l-(c_1/d_1))*exp(-d_1*z0_on_l)+b_1*c_1/d_1)
      pm10 = -(a_1*z10_on_l+b_1*(z10_on_l-(c_1/d_1))*exp(-d_1*z10_on_l)+b_1*c_1/d_1)
      pm1  = -(a_1*z_on_l+b_1*(z_on_l-(c_1/d_1))*exp(-d_1*z_on_l)+b_1*c_1/d_1)
      integralh = neutral-(ph1-ph0)
      integralm = neutral-(pm1-pm0)
      integralm10 = neutral10-(pm1-pm10)
    endwhere
    ttop       = temp-tstar/vkar*integralh
    qtop       = mixr-qstar/vkar*integralh
    utop       = umag-ustar/vkar*integralm

    where (fn%bldheight.le.z10) ! above canyon
      pg%u10=max(umag-ustar/vkar*integralm10,0.)
    end where

    ! assume standard stability functions hold for urban canyon (needs more work)
    tsurf=dg%rdsndelta*rdsntemp+(1.-dg%rdsndelta)*((1.-fn%sigmaveg)*road%temp(1)+fn%sigmaveg*pg%vegtemp)
    n=max(min((veg%moist-swilt)/(sfc-swilt),1.),0.)
    wf=(1.-dg%rdsndelta)*((1.-fn%sigmaveg)*dg%roaddelta+fn%sigmaveg*((1.-dg%vegdelta)*n+dg%vegdelta))
    call getqsat(ufull,qsurf,tsurf,dg%sigd)
    qsurf=qsurf*wf
    thetav=ttop*(1.+0.61*qtop)
    sthetav=tsurf*(1.+0.61*qsurf)
    lna=2.3
    n=log(fn%bldheight/zonet)
    af=vkar*vkar/(n*n)
    aft=vkar*vkar/(n*(lna+n))
    ri=min(grav*fn%bldheight*(1.-sthetav/thetav)/umag**2,rimax)
    where (ri>0.)
      fm=1./(1.+bprm*ri)**2
      fh=fm
    elsewhere
      root=sqrt(-ri*fn%bldheight/zonet)
      denma=1.+cms*2.*bprm*af*root
      denha=1.+chs*2.*bprm*aft*exp(0.5*lna)*root
      fm=1.-2.*bprm *ri/denma
      fh=1.-2.*bprm*ri/denha
    end where
    cd=af*fm
    thetavstar=aft*fh*(thetav-sthetav)/sqrt(cd)
    call dyerhicks(ufull,integralh,z_on_l,cd,thetavstar,thetav,sthetav,utop,fn%bldheight,n,lna)
    ustar=sqrt(cd)*utop
    tstar=vkar*(ttop-tsurf)/integralh
    qstar=vkar*(qtop-qsurf)/integralh
    z0_on_l   = z0*z_on_l/fn%bldheight
    z10_on_l  = max(z10,fn%bldheight)*z_on_l/fn%bldheight
    z0_on_l=min(z0_on_l,10.)
    z10_on_l=min(z10_on_l,10.)
    neutral = log(fn%bldheight/z0)
    neutral10 = log(fn%bldheight/max(z10,fn%bldheight))
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
      pm0 = -(a_1*z0_on_l+b_1*(z0_on_l-(c_1/d_1))*exp(-d_1*z0_on_l)+b_1*c_1/d_1)
      pm10 = -(a_1*z10_on_l+b_1*(z10_on_l-(c_1/d_1))*exp(-d_1*z10_on_l)+b_1*c_1/d_1)
      pm1  = -(a_1*z_on_l+b_1*(z_on_l-(c_1/d_1))*exp(-d_1*z_on_l)+b_1*c_1/d_1)
      integralh = neutral-(ph1-ph0)
      integralm = neutral-(pm1-pm0)
      integralm10 = neutral10-(pm1-pm10)
    endwhere

    pg%tscrn       = temp-tstar/vkar*integralh
    pg%qscrn       = mixr-qstar/vkar*integralh
    pg%uscrn=max(umag-ustar/vkar*integralm,0.)
    where (fn%bldheight.gt.z10) ! within canyon
      pg%u10=max(umag-ustar/vkar*integralm10,0.)
    end where

end select
pg%qscrn       = max(pg%qscrn,1.E-4)
      
return
end subroutine scrncalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Disables aTEB so subroutine calls have no effect

subroutine atebdisable(diag)

implicit none

integer, intent(in) :: diag

if (diag.ge.1) write(6,*) "Disable aTEB"
ufull=0

return
end subroutine atebdisable

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module ateb
