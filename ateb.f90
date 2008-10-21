
! This code was originally based on the TEB scheme of Masson, Boundary-Layer Meteorology, 94, p357 (2000)
! The snow scheme is based on Douville, Royer and Mahfouf, Climate Dynamics, 12, p21 (1995)

! Usual practice is:
!   call tebinit     ! to initalise state arrays, etc (use tebdisable to disable calls to ateb subroutines)
!   call tebloadm    ! to load previous state arrays (from tebsavem)
!   call tebtype     ! to define urban type (or use tebfndef to define urban properties at each grid point)
!   ...
!   do t=1,tmax
!     ...
!     call tebnewangle ! store current solar zenith and azimuthal angle (use tebccangle for CCAM or
!                        use tebnewangle1 for a single grid point)
!     call tebcalc     ! calculates urban temperatures, fluxes, etc and blends with input
!     call tebalb      ! blends input and urban albedo (use tebalb1 for a single grid point)
!     call tebzo       ! blends input and urban roughness lengths for momentum and heat
!     ...
!   end do
!   ...
!   call tebsavem    ! to save current state arrays (for use by tebloadm)
!   call tebend      ! to deallocate memory before quiting

! only tebinit and tebcalc are manditory.  All other subroutine calls are optional.

! URBAN TYPES:
 
! 1 = Urban (TAPM 31)
! 2 = Urban (low) (TAPM 32)
! 3 = Urban (medium) (TAPM 33)
! 4 = Urban (high) (TAPM 34)
! 5 = Urban (cbd) (TAPM 35)
! 6 = Industrial (low) (TAPM 36)
! 7 = Industrial (medium) (TAPM 37)
! 8 = Industrial (high) (TAPM 38)

! NOTES: 
!  Below are some differences between the TEB (Masson 2000) scheme and aTEB:
!
! - aTEB uses two walls instead of the TEB single wall.  Also, only up to nrefl order reflections are used in aTEB for
!   longwave and short wave radation (nrefl=3 by default).  In TEB, infinite reflections are used for shortwave, but
!   only 1st order for long wave.  See Harman, et. al (2004) for a complete treatment of the reflections.
!
! - aTEB uses an iterative scheme to solve for the sensible heat flux budget in the canyon.
!
! - aTEB returns ln(zom/zot)=1.29*(ustar*zom/nu)**0.25-2. (Kanda et al, 2007) for the urban cover.  In the model, we use
!   ln(zom/zot)=6 for canyon and roof surfaces (Kanda et al 2005), but ln(zom/zot)=2. between the canyon and the atmosphere.
!
! - aTEB calculates aerodynamic resistances for the recirculation and ventilation regions of the canyon (see Harman, et. al 2004).
!   This approach takes advantage of the second wall temperature (i.e., the fluxes depend on the wind direction).
!
! Minor differences include the use of CSIRO9 stability coefficents (McGregor et al 1993), instead of Mascart et al (1995).
! Also there are changes to the displacement height (not at roof height),  urban properties are based on TAPM classes,
! the time dependence of the traffic sensible heat flux follows Coutts et al (2007), etc
!

module ateb

implicit none

private
public tebinit,tebcalc,tebend,tebzo,tebload,tebsave,tebtype,tebfndef,tebalb,tebalb1, &
       tebnewangle,tebnewangle1,tebccangle,tebdisable,tebloadm,tebsavem,tebcd

! type definitions
type tatm
  real :: sg,rg,rho,temp,mixr,ps,pa,umag,udir,rnd,snd
end type tatm
type tout
  real :: fg,eg,ts,wf
end type tout
type trad
  real :: roof,road,walle,wallw,rfsn,rdsn
end type trad
type tdiag
  real :: roofdelta,roaddelta,rfsndelta,rdsndelta
  real :: tempc,mixrc,tempr,mixrr,sigd,sigr,rfdzmin
  real :: accool,canyonrgout,roofrgout
end type tdiag
type tsurf
  real, dimension(3) :: temp
  real :: water,snow,den,alpha
end type tsurf
type twall
  real, dimension(3) :: temp
end type twall
type tdata
  real :: hwratio,sigmabld,industryfg,trafficfg,bldheight,zo,vangle,hangle
  real :: roofalpha,wallalpha,roadalpha,ctime
end type tdata
type tprog
  real :: lzom,lzoh,cndzmin,cduv
end type tprog

! state arrays
integer, save :: ufull
integer, dimension(:), allocatable, save :: ugrid,mgrid
real, dimension(:), allocatable, save :: sigmau
type(tsurf), dimension(:), allocatable, save :: roof,road,roofadj,roadadj
type(twall), dimension(:), allocatable, save :: walle,wallw,walleadj,wallwadj
type(tdata), dimension(:), allocatable, save :: fn
type(tprog), dimension(:), allocatable, save :: pg
! Building parameters
!real, dimension(3), parameter :: roofdepth =(/ 0.05,0.4,0.1 /)          ! depth (m) (Ecoclim)
real, dimension(3), parameter :: roofdepth =(/ 0.05,0.05,0.4 /)          ! depth (m)
real, dimension(3), parameter :: walldepth =(/ 0.02,0.125,0.05 /)
real, dimension(3), parameter :: roaddepth =(/ 0.05,0.1,1. /)
!real, dimension(3), parameter :: roofcp =(/ 2.11E6,0.28E6,0.29E6 /)     ! heat capacity (J m^-3 K^-1) (Ecoclim)
real, dimension(3), parameter :: roofcp =(/ 2.11E6,2.11E6,0.29E6 /)     ! heat capacity (J m^-3 K^-1)
real, dimension(3), parameter :: wallcp =(/ 1.55E6,1.55E6,0.29E6 /)
real, dimension(3), parameter :: roadcp =(/ 1.94E6,1.28E6,1.28E6 /)
!real, dimension(3), parameter :: rooflambda =(/ 1.51,0.08,0.05 /)       ! conductance (W m^-1 K^-1) (Ecoclim)
real, dimension(3), parameter :: rooflambda =(/ 1.51,1.51,0.05 /)       ! conductance (W m^-1 K^-1)
real, dimension(3), parameter :: walllambda =(/ 0.9338,0.9338,0.05 /)
real, dimension(3), parameter :: roadlambda =(/ 0.7454,0.2513,0.2513 /)
real, parameter :: roofemiss=0.90    ! emissitivity
real, parameter :: wallemiss=0.85 
real, parameter :: roademiss=0.94
real, parameter :: bldtemp=291.16    ! Comfort temperature (K) = 18deg C
real, parameter :: zocanyon=3.E-4    ! Roughness lengths for canyon surfaces (m)
real, parameter :: zoroof=zocanyon   ! Roughness length for rooftops (m)
real, parameter :: maxroofwater=1.   ! max water on roof (kg m^-2)
real, parameter :: maxroadwater=1.   ! max water on road (kg m^-2)
real, parameter :: maxrfsn=1.        ! max snow on roof (kg m^-2)
real, parameter :: maxrdsn=1.        ! max snow on road (kg m^-2)
real, parameter :: refheight=0.4     ! displacement height as fraction of building height (Kanda et al 2007)
integer, parameter :: resmeth=1      ! Canyon sensible heat transfer (0=Masson, 1=Harman, 2=Kusaka)
integer, parameter :: zohmeth=1      ! Urban roughness length for heat (0=Veg, 1=Kanda 2007, 2=Kanda 2005)
integer, parameter :: acmeth=1       ! AC heat pump into canyon (0=Off, 1=On)
integer, parameter :: nrefl=3        ! Number of canyon reflections (default=3)
! Other parameters (e.g., snow, water, physics, etc)
real, parameter :: waterden=1000.    ! water density (kg m^-3)
real, parameter :: icelambda=2.22    ! conductance of ice (W m^-1 K^-1)
real, parameter :: aircp=1004.64     ! Heat capapcity of dry air (J kg^-1 K^-1)
real, parameter :: icecp=2100.       ! Heat capacity of ice (J kg^-1 K^-1)
real, parameter :: grav=9.80616      ! gravity (m s^-2)
real, parameter :: vkar=0.4          ! von Karman constant
real, parameter :: lv=2.501e6        ! Latent heat of vaporisation
real, parameter :: lf=3.337e5        ! Latent heat of fusion
real, parameter :: ls=lv+lf          ! Latent heat of sublimation
real, parameter :: pi=3.1415927      ! pi
real, parameter :: rd=287.04         ! Gas constant for dry air
real, parameter :: sbconst=5.67e-8   ! Stefan-Boltzmann constant
real, parameter :: zosnow=0.001      ! Roughness length for snow (m)
real, parameter :: snowemiss=1.      ! snow emissitivity
real, parameter :: maxsnowalpha=0.85 ! max snow albedo
real, parameter :: minsnowalpha=0.5  ! min snow albedo
real, parameter :: maxsnowden=300.   ! max snow density (kg m^-3)
real, parameter :: minsnowden=100.   ! min snow density (kg m^-3)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine prepare the arrays used by the aTEB scheme
! This is a compulsory subroutine that must be called during
! model initalisation

subroutine tebinit(ifull,sigu,zmin,diag)

implicit none

integer, intent(in) :: ifull,diag
integer iqu,iq,ii
real, intent(in) :: zmin
real, dimension(ifull), intent(in) :: sigu
integer, dimension(ifull) :: utype

if (diag.ge.1) write(6,*) "Initialising aTEB"

ufull=count(sigu.gt.0.)
if (ufull.eq.0) return

allocate(ugrid(ufull),mgrid(ifull),fn(ufull),pg(ufull))
allocate(roof(ufull),road(ufull),walle(ufull),wallw(ufull))
allocate(roofadj(ufull),roadadj(ufull),walleadj(ufull),wallwadj(ufull))
allocate(sigmau(ufull))

! define grid arrays
ugrid=0
mgrid=0
sigmau=0.
iqu=0
do iq=1,ifull
  if (sigu(iq).gt.0.) then
    iqu=iqu+1
    ugrid(iqu)=iq
    mgrid(iq)=iqu
    sigmau(iqu)=sigu(iq)
  end if
end do

! Initialise state variables
do ii=1,3
  roof%temp(ii)=bldtemp
  road%temp(ii)=bldtemp
  walle%temp(ii)=bldtemp
  wallw%temp(ii)=bldtemp
  roofadj%temp(ii)=0.
  roadadj%temp(ii)=0.
  walleadj%temp(ii)=0.
  wallwadj%temp(ii)=0.
end do
roof%water=0.
road%water=0.
roof%snow=0.
road%snow=0.
roof%den=minsnowden
road%den=minsnowden
roof%alpha=maxsnowalpha
road%alpha=maxsnowalpha
roofadj%water=0.
roadadj%water=0.
roofadj%snow=0.
roadadj%snow=0.
roofadj%den=0.
roadadj%den=0.
roofadj%alpha=0.
roadadj%alpha=0.

fn%hwratio=0.
fn%sigmabld=0.
fn%industryfg=0.
fn%trafficfg=0.
fn%bldheight=0.
fn%zo=0.
fn%roofalpha=0.
fn%wallalpha=0.
fn%roadalpha=0.
fn%vangle=0.
fn%hangle=0.
fn%ctime=0.

utype=1 ! default urban
call tebtype(ifull,utype,diag)

pg%cndzmin=max(zmin,fn%zo+1.)           ! updated in tebcalc
pg%lzom=log(pg%cndzmin/fn%zo)           ! updated in tebcalc
pg%lzoh=6.+pg%lzom ! (Kanda et al 2005) ! updated in tebcalc
pg%cduv=(vkar/pg%lzom)**2               ! updated in tebcalc

return
end subroutine tebinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine deallocates arrays used by the TEB scheme

subroutine tebend(diag)

implicit none

integer, intent(in) :: diag

if (ufull.eq.0) return
if (diag.ge.1) write(6,*) "Deallocating aTEB arrays"

deallocate(ugrid,mgrid,fn)
deallocate(roof,road,walle,wallw)
deallocate(roofadj,roadadj,walleadj,wallwadj)
deallocate(sigmau)

return
end subroutine tebend

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine loads aTEB state arrays (not compulsory)

subroutine tebload(ifull,urban,diag)

implicit none

integer, intent(in) :: ifull,diag
integer ii
real, dimension(ifull,20), intent(in) :: urban

if (ufull.eq.0) return
if (diag.ge.1) write(6,*) "Load aTEB state arrays"

do ii=1,3
  roof%temp(ii)=urban(ugrid,ii)
  walle%temp(ii)=urban(ugrid,ii+3)
  wallw%temp(ii)=urban(ugrid,ii+6)
  road%temp(ii)=urban(ugrid,ii+9)
end do
roof%water=urban(ugrid,13)
road%water=urban(ugrid,14)
roof%snow=urban(ugrid,15)
road%snow=urban(ugrid,16)
roof%den=urban(ugrid,17)
road%den=urban(ugrid,18)
roof%alpha=urban(ugrid,19)
road%alpha=urban(ugrid,20)

return
end subroutine tebload

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! temperature only version of tebsave

subroutine tebloadm(ifull,urban,diag)

implicit none

integer, intent(in) :: ifull,diag
integer ii
real, dimension(ifull,12), intent(in) :: urban

if (ufull.eq.0) return
if (diag.ge.1) write(6,*) "Load aTEB state arrays"

do ii=1,3
  roof%temp(ii)=urban(ugrid,ii)
  walle%temp(ii)=urban(ugrid,ii+3)
  wallw%temp(ii)=urban(ugrid,ii+6)
  road%temp(ii)=urban(ugrid,ii+9)
end do

return
end subroutine tebloadm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine loads aTEB type arrays (not compulsory)

subroutine tebtype(ifull,itype,diag)

implicit none

integer, intent(in) :: ifull,diag
integer, dimension(ifull), intent(in) :: itype
integer, parameter :: maxtype = 8
! Urban fraction (defined in host model)
! real, dimension(maxtype), parameter :: sigu=(/ 0.6, 0.5, 0.6, 0.7, 0.95, 0.5, 0.6, 0.7 /)
! Building height to width ratio
real, dimension(maxtype), parameter :: chwratio=(/ 1., 0.4, 0.5, 0.6, 2., 0.5, 1., 1.5 /)
! Area fraction occupied by buildings
real, dimension(maxtype), parameter :: csigmabld=(/ 0.6, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7 /)
! Industral sensible heat flux (W m^-2)
real, dimension(maxtype), parameter :: cindustryfg=(/ 20., 8., 12., 16., 28., 20., 40., 60. /)
! Daily averaged traffic sensible heat flux (W m^-2)
real, dimension(maxtype), parameter :: ctrafficfg=(/ 5., 2., 3., 4., 7., 5., 10., 15. /)
! Building height (m)
real, dimension(maxtype), parameter :: cbldheight=(/ 10., 4., 6., 8., 20., 5., 10., 15. /)
! Effective roughness length (m)
real, dimension(maxtype), parameter :: czo=0.1*cbldheight
! Roof albedo
real, dimension(maxtype), parameter :: croofalpha=(/ 0.20, 0.23, 0.20, 0.17, 0.13, 0.20, 0.20, 0.20 /)
! Wall albedo
real, dimension(maxtype), parameter :: cwallalpha=(/ 0.33, 0.37, 0.33, 0.29, 0.22, 0.33, 0.33, 0.33 /)
! Road albedo
real, dimension(maxtype), parameter :: croadalpha=(/ 0.10, 0.11, 0.10, 0.09, 0.07, 0.10, 0.10, 0.10 /)

if (ufull.eq.0) return
if (diag.ge.1) write(6,*) "Load aTEB building properties"

if (any(itype.lt.1).or.any(itype.gt.maxtype)) then
  write(6,*) "ERROR: Urban type is out of range"
  stop
end if

fn%hwratio=chwratio(itype(ugrid))
fn%sigmabld=csigmabld(itype(ugrid))
fn%industryfg=cindustryfg(itype(ugrid))
fn%trafficfg=ctrafficfg(itype(ugrid))
fn%bldheight=cbldheight(itype(ugrid))
fn%zo=czo(itype(ugrid))
fn%roofalpha=croofalpha(itype(ugrid))
fn%wallalpha=cwallalpha(itype(ugrid))
fn%roadalpha=croadalpha(itype(ugrid))

return
end subroutine tebtype

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine specifies the urban properties for each grid point
!

subroutine tebfndef(ifull,ifn,diag)

implicit none

integer, intent(in) :: ifull,diag
real, dimension(ifull,9), intent(in) :: ifn

if (ufull.eq.0) return
if (diag.ge.1) write(6,*) "Load aTEB building properties"

fn%hwratio=ifn(ugrid,1)
fn%sigmabld=ifn(ugrid,2)
fn%industryfg=ifn(ugrid,3)
fn%trafficfg=ifn(ugrid,4)
fn%bldheight=ifn(ugrid,5)
fn%zo=ifn(ugrid,6)
fn%roofalpha=ifn(ugrid,7)
fn%wallalpha=ifn(ugrid,8)
fn%roadalpha=ifn(ugrid,9)

return
end subroutine tebfndef

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine saves aTEB state arrays (not compulsory)

subroutine tebsave(ifull,urban,diag)

implicit none

integer, intent(in) :: ifull,diag
integer ii
real, dimension(ifull,20), intent(inout) :: urban

if (ufull.eq.0) return
if (diag.ge.1) write(6,*) "Save aTEB state arrays"

do ii=1,3
  urban(ugrid,ii)=roof%temp(ii)
  urban(ugrid,ii+3)=walle%temp(ii)
  urban(ugrid,ii+6)=wallw%temp(ii)
  urban(ugrid,ii+9)=road%temp(ii)
end do
urban(ugrid,13)=roof%water
urban(ugrid,14)=road%water
urban(ugrid,15)=roof%snow
urban(ugrid,16)=road%snow
urban(ugrid,17)=roof%den
urban(ugrid,18)=road%den
urban(ugrid,19)=roof%alpha
urban(ugrid,20)=road%alpha

return
end subroutine tebsave

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! temperature only version of tebsave

subroutine tebsavem(ifull,urban,diag)

implicit none

integer, intent(in) :: ifull,diag
integer ii
real, dimension(ifull,12), intent(inout) :: urban

if (ufull.eq.0) return
if (diag.ge.1) write(6,*) "Save aTEB state arrays"

do ii=1,3
  urban(ugrid,ii)=roof%temp(ii)
  urban(ugrid,ii+3)=walle%temp(ii)
  urban(ugrid,ii+6)=wallw%temp(ii)
  urban(ugrid,ii+9)=road%temp(ii)
end do

return
end subroutine tebsavem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine blends urban momentum and heat roughness lengths
! (This version neglects the displacement height (e.g., for CCAM))
!

subroutine tebzo(ifull,zom,zoh,diag)

implicit none

integer, intent(in) :: ifull,diag
real, dimension(ifull), intent(inout) :: zom,zoh
real, dimension(ufull) :: workb,workc
real, parameter :: zr=1.e-15 ! limits minimum roughness length for heat

if (ufull.eq.0) return
if (diag.ge.1) write(6,*) "Blend urban roughness lengths"

! evaluate at canyon displacement height (really the atmospheric model should provide a displacement height)
workb=sqrt((1.-sigmau)/log(pg%cndzmin/zom(ugrid))**2+sigmau/pg%lzom**2)
workc=(1.-sigmau)/(log(pg%cndzmin/zom(ugrid))*log(pg%cndzmin/zoh(ugrid)))+sigmau/(pg%lzom*pg%lzoh)
workc=workc/workb
workb=pg%cndzmin*exp(-1./workb)
workc=max(pg%cndzmin*exp(-1./workc),zr)
zom(ugrid)=workb
zoh(ugrid)=workc

if (any(zoh(ugrid).le.zr)) write(6,*) "WARN: minimum zoh reached"

return
end subroutine tebzo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine blends the urban drag coeff
!

subroutine tebcd(ifull,cduv,diag)

implicit none

integer, intent(in) :: ifull,diag
real, dimension(ifull), intent(inout) :: cduv

if (ufull.eq.0) return

cduv(ugrid)=(1.-sigmau)*cduv(ugrid)+sigmau*pg%cduv

return
end subroutine tebcd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates the urban contrabution to albedo.
! (all grid points)

subroutine tebalb(ifull,alb,diag)

implicit none

integer, intent(in) :: ifull,diag
real, dimension(ifull), intent(inout) :: alb
real, dimension(ufull) :: ualb

if (ufull.eq.0) return
if (diag.ge.1) write(6,*) "Calculate urban albedo"

call tebalbcalc(1,ufull,ualb,diag)
alb(ugrid)=(1.-sigmau)*alb(ugrid)+sigmau*ualb

return
end subroutine tebalb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates the urban contrabution to albedo.
! (selected grid points only)

subroutine tebalb1(is,ifull,alb,diag)

implicit none

integer, intent(in) :: is,ifull,diag
integer, dimension(ifull) :: lgrid
integer i,ucount,ib,ie
real, dimension(ifull), intent(inout) :: alb
real, dimension(ifull) :: ualb

if (ufull.eq.0) return

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

call tebalbcalc(ib,ucount,ualb(1:ucount),diag)
alb(lgrid(1:ucount))=(1.-sigmau(ib:ie))*alb(lgrid(1:ucount)) &
                     +sigmau(ib:ie)*ualb(1:ucount)

return
end subroutine tebalb1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Albedo calculations

subroutine tebalbcalc(is,ifull,alb,diag)

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
albu=1.-(fn(is:ie)%hwratio*(sg%walle+sg%wallw)*(1.-fn(is:ie)%wallalpha) &
    +sg%road*((1.-snowdelta)*(1.-fn(is:ie)%roadalpha)+snowdelta*(1.-road(is:ie)%alpha)))
snowdelta=roof(is:ie)%snow/(roof(is:ie)%snow+maxrfsn)
albr=(1.-snowdelta)*fn(is:ie)%roofalpha+snowdelta*roof(is:ie)%alpha
alb=fn(is:ie)%sigmabld*albr+(1.-fn(is:ie)%sigmabld)*albu

return
end subroutine tebalbcalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine stores the zenith angle and the solar azimuth angle
! (all grid points)

! ifull = size of array
! cosin = array of cosine of zenith angles
! azimuthin = array of azimuthal angles
! diag = dialogue flag (0=off)

subroutine tebnewangle(ifull,cosin,azimuthin,ctimein,diag)

implicit none

integer, intent(in) :: ifull,diag
real, dimension(ifull), intent(in) :: cosin     ! cosine of zenith angle
real, dimension(ifull), intent(in) :: azimuthin ! azimuthal angle
real, dimension(ifull), intent(in) :: ctimein   ! local hour (0<=ctime<=1)

if (ufull.eq.0) return
if (diag.ge.1) write(6,*) "Update solar zenith angle and azimuth angle"

fn%hangle=0.5*pi-azimuthin(ugrid)
fn%vangle=acos(cosin(ugrid))
fn%ctime=ctimein(ugrid)  

return
end subroutine tebnewangle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine stores the zenith angle and the solar azimuth angle
! (single grid point)

! iq = grid point
! cosin = cosine of zenith angle
! azmiuthin = azimuthal angle (rad)

subroutine tebnewangle1(is,ifull,cosin,azimuthin,ctimein)

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
end subroutine tebnewangle1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This version of tebnewangle is for CCAM
!

subroutine tebccangle(is,ifull,cosin,rlon,rlat,fjd,slag,dhr,dlt)

implicit none

integer, intent(in) :: is,ifull
real, intent(in) :: fjd,slag,dhr,dlt
real, dimension(ifull), intent(in) :: cosin,rlon,rlat
real, dimension(ifull) :: hloc,x,y

if (ufull.eq.0) return

where (mgrid(is:ifull+is-1).ge.1)
  ! from CCAM zenith.f
  hloc=2.*pi*fjd+slag+pi+rlon+dhr*pi/24.
  ! estimate azimuth angle
  x=sin(-hloc)*cos(dlt)
  y=-cos(-hloc)*cos(dlt)*sin(rlat)+cos(rlat)*sin(dlt)
  !azimuth=atan2(x,y)
  fn(mgrid(is:ifull+is-1))%hangle=0.5*pi-atan2(x,y)
  fn(mgrid(is:ifull+is-1))%vangle=acos(cosin)
  fn(mgrid(is:ifull+is-1))%ctime=min(max(mod(hloc/(2.*pi),1.),0.),1.)
end where

return
end subroutine tebccangle


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
! ps = surface pressure (Pa)
! pa = pressure at first model level (Pa)
! umag = horizontal wind speed at first model level (m/s)
! ofg = Input/Output sensible heat flux (W/m^2)
! oeg = Input/Output latient heat flux (W/m^2)
! ots = Input/Output radiative/skin temperature (K)
! owf = Input/Output wetness fraction/surface water (%)

subroutine tebcalc(ifull,ofg,oeg,ots,owf,ddt,zmin,sg,rg,rnd,rho,temp,mixr,ps,pa,uu,vv,umin,diag)

implicit none

integer, intent(in) :: ifull,diag
real, intent(in) :: ddt,zmin,umin
real, dimension(ifull), intent(in) :: sg,rg,rnd,rho,temp,mixr,ps,pa,uu,vv
real, dimension(ifull), intent(inout) :: ofg,oeg,ots,owf
type (tatm), dimension(ufull) :: atm
type (tout), dimension(ufull) :: uo

if (ufull.eq.0) return

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

call tebeval(uo,ddt,atm,zmin,diag)

ofg(ugrid)=(1.-sigmau)*ofg(ugrid)+sigmau*uo%fg
oeg(ugrid)=(1.-sigmau)*oeg(ugrid)+sigmau*uo%eg
ots(ugrid)=(1.-sigmau)*ots(ugrid)+sigmau*uo%ts
owf(ugrid)=(1.-sigmau)*owf(ugrid)+sigmau*uo%wf

return
end subroutine tebcalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! urban flux calculations

subroutine tebeval(uo,ddt,atm,zmin,diag)

implicit none

integer, intent(in) :: diag
integer iqu,j,k,ii,cns,cnr
integer, dimension(ufull) :: igs,igr
real, intent(in) :: ddt,zmin
real, dimension(ufull,3) :: garoof,gawalle,gawallw,garoad
real, dimension(ufull) :: garfsn,gardsn
real, dimension(ufull) :: rdsntemp,rfsntemp,rdsnmelt,rfsnmelt
real, dimension(ufull) :: wallpsi,roadpsi,fgtop,egtop,qsatr
real, dimension(ufull) :: oldval,newval,topu,cu,ctmax,ctmin,evctx,evct
real, dimension(ufull) :: ln,rn,we,ww,wr,zolog,a,xe,xw,cuven,n,zom
real, dimension(ufull) :: p_sntemp,p_netldratio,p_wallpsi,p_roadpsi
real, dimension(ufull) :: p_fgtop,p_gasn,p_snmelt
logical, save :: firstcall=.true.
type(tatm), dimension(ufull), intent(in) :: atm
type(tatm), dimension(ufull) :: p_atm
type(tout), dimension(ufull), intent(out) :: uo
type(trad), dimension(ufull) :: sg,rg,fg,eg
type(trad), dimension(ufull) :: p_sg,p_rg,p_fg,p_eg
type(trad), dimension(ufull) :: acond,p_acond
type(tdiag), dimension(ufull) :: dg,p_dg
type(tsurf), dimension(ufull) :: roofdum,roaddum,croof,croad,rooforg,roadorg
type(tsurf), dimension(ufull) :: p_rodum
type(twall), dimension(ufull) :: walledum,wallwdum,cwalle,cwallw,walleorg,wallworg
type(twall), dimension(ufull) :: p_walledum,p_wallwdum
type(tdata), dimension(ufull) :: p_fn
type(tprog), dimension(ufull) :: p_pg

if (diag.ge.1) write(6,*) "Evaluating aTEB"

! null values
igs=0
igr=0
acond%roof=0.
acond%road=0.
acond%walle=0.
acond%wallw=0.
acond%rfsn=0.
acond%rdsn=0.
p_rg=acond
p_fg=acond
p_eg=acond
dg%roofdelta=0.
dg%roaddelta=0.
dg%rfsndelta=0.
dg%rfsndelta=0.
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

! canyon (displacement height at refheight*building height)
dg%sigd=atm%ps
dg%tempc=atm%temp*(dg%sigd/atm%pa)**(rd/aircp)
call getqsat(ufull,qsatr,dg%tempc,dg%sigd)
call getqsat(ufull,a,atm%temp,atm%pa)
dg%mixrc=atm%mixr*qsatr/a ! a=qsata

! roof (displacement height at building height)
dg%sigr=atm%ps*exp(grav*fn%bldheight*(refheight-1.)/(rd*atm%temp))
dg%tempr=atm%temp*(dg%sigr/atm%pa)**(rd/aircp)
call getqsat(ufull,qsatr,dg%tempr,dg%sigr)
dg%mixrr=atm%mixr*qsatr/a ! a=qsata

! new snowfall
where (atm%snd.gt.0.)
  roof%den=(roof%snow*roof%den+atm%snd*ddt*minsnowden)/(roof%snow+ddt*atm%snd)
  road%den=(road%snow*road%den+atm%snd*ddt*minsnowden)/(road%snow+ddt*atm%snd)
  roof%alpha=maxsnowalpha
  road%alpha=maxsnowalpha
end where

! prep predictor-corrector arrays
roofdum=roof
walledum=walle
wallwdum=wallw
roaddum=road

do j=1,2 ! predictor-corrector loop -------------------------------
  ! limit state variables
  do ii=1,3
    roofdum%temp(ii)=min(max(roofdum%temp(ii),200.),400.)
    walledum%temp(ii)=min(max(walledum%temp(ii),200.),400.)
    wallwdum%temp(ii)=min(max(wallwdum%temp(ii),200.),400.)
    roaddum%temp(ii)=min(max(roaddum%temp(ii),200.),400.)
  end do
  roofdum%water=min(max(roofdum%water,0.),maxroofwater)
  roaddum%water=min(max(roaddum%water,0.),maxroadwater)
  roofdum%snow=min(max(roofdum%snow,0.),maxrfsn)
  roaddum%snow=min(max(roaddum%snow,0.),maxrdsn)
  roofdum%den=min(max(roofdum%den,minsnowden),maxsnowden)
  roaddum%den=min(max(roaddum%den,minsnowden),maxsnowden)
  roofdum%alpha=min(max(roofdum%alpha,minsnowalpha),maxsnowalpha)
  roaddum%alpha=min(max(roaddum%alpha,minsnowalpha),maxsnowalpha)
    
  ! water and snow cover fractions
  dg%roofdelta=(roofdum%water/maxroofwater)**(2./3.)
  dg%roaddelta=(roaddum%water/maxroadwater)**(2./3.)
  dg%rfsndelta=roofdum%snow/(roofdum%snow+maxrfsn)
  dg%rdsndelta=roaddum%snow/(roaddum%snow+maxrdsn)

  ! Adjust canyon roughness to include snow
  n=roaddum%snow/(roaddum%snow+maxrdsn+0.408*grav*fn%zo)               ! snow cover for urban roughness calc (Douville, et al 1995)
  zom=(1.-n)*fn%zo+n*zosnow                                            ! blend urban and snow roughness length
  dg%rfdzmin=max(abs(zmin-fn%bldheight*(1.-refheight)),zoroof+1.)      ! distance to roof displacement height
  where (zmin.ge.fn%bldheight*(1.-refheight))
    pg%lzom=log(zmin/zom)
    topu=(2./pi)*atm%umag*log(fn%bldheight*(1.-refheight)/zom)/pg%lzom ! wind speed at canyon top
  elsewhere
    pg%lzom=log(fn%bldheight*(1.-refheight)/zom)*exp(0.5*fn%hwratio*(1.-zmin/(fn%bldheight*(1.-refheight))))
    topu=(2./pi)*atm%umag*exp(0.5*fn%hwratio*(1.-zmin/(fn%bldheight*(1.-refheight))))
  end where
  pg%cndzmin=zom*exp(pg%lzom)                                          ! distance to canyon displacement height

  ! calculate heat conduction (e.g., through walls)
  ! done here to estimate cooling heat flux (e.g., AC) into canyon
  do k=1,2
    garoof(:,k)=2.*(roofdum%temp(k)-roofdum%temp(k+1))/(roofdepth(k)/rooflambda(k)+roofdepth(k+1)/rooflambda(k+1))
    gawalle(:,k)=2.*(walledum%temp(k)-walledum%temp(k+1))/(walldepth(k)/walllambda(k)+walldepth(k+1)/walllambda(k+1))
    gawallw(:,k)=2.*(wallwdum%temp(k)-wallwdum%temp(k+1))/(walldepth(k)/walllambda(k)+walldepth(k+1)/walllambda(k+1))
    garoad(:,k)=2.*(roaddum%temp(k)-roaddum%temp(k+1))/(roaddepth(k)/roadlambda(k)+roaddepth(k+1)/roadlambda(k+1))
  end do
  garoof(:,3)=2.*rooflambda(3)*(roofdum%temp(3)-bldtemp)/(roofdepth(3))
  gawalle(:,3)=2.*walllambda(3)*(walledum%temp(3)-bldtemp)/(walldepth(3))
  gawallw(:,3)=2.*walllambda(3)*(wallwdum%temp(3)-bldtemp)/(walldepth(3))
  garoad(:,3)=0.
  if (acmeth.eq.1) then
    dg%accool=max(0.,garoof(:,3)+gawalle(:,3)+gawallw(:,3)) ! should be divided by efficency factor
  else
    dg%accool=0.
  end if

  ! calculate shortwave radiation (up to 2nd order reflections)
  call getswcoeff(ufull,sg,wallpsi,roadpsi,fn,roaddum%alpha,dg%rdsndelta)
  sg%roof=(1.-fn%roofalpha)*sg%roof*atm%sg
  sg%walle=(1.-fn%wallalpha)*sg%walle*atm%sg
  sg%wallw=(1.-fn%wallalpha)*sg%wallw*atm%sg
  sg%road=(1.-fn%roadalpha)*sg%road*atm%sg
  sg%rfsn=(1.-roofdum%alpha)*sg%rfsn*atm%sg
  sg%rdsn=(1.-roaddum%alpha)*sg%rdsn*atm%sg

  ! calculate canyon wind speed and bulk transfer coefficents
  select case(resmeth)
    case(0) ! Masson (2000)
      cu=topu*exp(-0.25*fn%hwratio)
      acond%road=cu ! bulk transfer coefficents are updated in solvecanyon
      acond%walle=cu
      acond%wallw=cu
      acond%rdsn=cu
    case(1) ! Harman et al (2004) - modified to include 2nd wall
      ln=fn%bldheight*max(0.,1./fn%hwratio-1.5)
      rn=max(0.,fn%bldheight*(1.-2.*ln/(3.*fn%bldheight)))
      ln=min(1.5*fn%bldheight,ln)
      cu=topu*exp(-0.9*sqrt(ln**2+(fn%bldheight-rn)**2)/fn%bldheight)
      a=0.15*max(1.,1.5*fn%hwratio)
      where (rn.ge.fn%bldheight) ! recirculation only
        xe=exp(-a)
        we=cu*(1.-xe)/a
        xw=exp(-a/fn%hwratio)
        wr=cu*fn%hwratio*xe*(1.-xw)/a
        ww=cu*xe*xw*(1.-exp(-a))/a
      elsewhere (rn.gt.0.) ! recirculation starts on east wall
        n=max(zocanyon,rn)
        cuven=topu*(1.-1./log(fn%bldheight/zocanyon)+n/((fn%bldheight-n)*(log(zocanyon/n) &
              /log(n/fn%bldheight)+1.)))
        cuven=max(cuven*(1.-rn/fn%bldheight),cu*(1.-exp(-a*(1.-rn/fn%bldheight)))/a)
        xe=exp(-a*rn/fn%bldheight)
        we=(cuven+cu*(1.-xe)/a)
        xw=exp(-a/fn%hwratio)
        wr=cu*fn%hwratio*xe*(1.-xw)/a
        ww=cu*xe*xw*(1.-exp(-a))/a
      elsewhere ! recirculation starts on road
        zolog=log(0.1*fn%bldheight/zocanyon)
        cuven=topu*zolog/(2.3+zolog)
        xe=exp(-a*(1./fn%hwratio-3.))
        cuven=fn%bldheight*max(cuven*(1./fn%hwratio-3.),cu*(1.-xe)/a)
        xw=exp(-a*3.)
        wr=fn%hwratio*(cuven/fn%bldheight+cu*(1.-xw)/a)
        cuven=topu*(1.-1./log(fn%bldheight/zocanyon)+zocanyon/(fn%bldheight-zocanyon))
        xe=cu*xe*(1.-exp(-a))/a
        we=max(cuven,xe)
        ww=cu*xw*(1.-exp(-a))/a
      end where
      zolog=log(0.1*fn%bldheight/zocanyon)
      a=vkar**2/(zolog*(6.+zolog)) ! Kanda et al (2005)
      n=abs(atm%udir)/pi
      acond%road=a*wr                ! road bulk transfer
      acond%walle=a*(n*ww+(1.-n)*we) ! east wall bulk transfer
      acond%wallw=a*(n*we+(1.-n)*ww) ! west wall bulk transfer
      zolog=log(0.1*fn%bldheight/zosnow)
      a=vkar**2/(zolog*(2.3+zolog))
      acond%rdsn=a*wr                ! road snow bulk transfer
    case(2) ! Kusaka et al (2001)
      cu=topu*exp(-0.25*fn%hwratio)
      where (cu.le.5.)
        acond%road=6.15+4.18*cu
      elsewhere
        acond%road=7.51*cu**0.78
      end where
      acond%walle=acond%road
      acond%wallw=acond%road
      acond%rdsn=acond%road  
  end select

  ! solve for road snow temperature -------------------------------
  ! includes solution to canyon temperature and longwave radiation
  cns=0
  cnr=0
  do iqu=1,ufull ! prepare index arrays
    if (roaddum(iqu)%snow.gt.0.) then
      cns=cns+1
      igs(cns)=iqu
    else
      cnr=cnr+1
      igr(cnr)=iqu
    end if
  end do
  if (cns.gt.0) then ! road snow
    p_rodum(1:cns)=roaddum(igs(1:cns)) ! pack
    p_walledum(1:cns)=walledum(igs(1:cns))
    p_wallwdum(1:cns)=wallwdum(igs(1:cns))
    p_dg(1:cns)=dg(igs(1:cns))
    p_sg(1:cns)=sg(igs(1:cns))
    p_atm(1:cns)=atm(igs(1:cns))
    p_acond(1:cns)=acond(igs(1:cns))
    p_wallpsi(1:cns)=wallpsi(igs(1:cns))
    p_roadpsi(1:cns)=roadpsi(igs(1:cns))
    p_fn(1:cns)=fn(igs(1:cns))
    p_pg(1:cns)=pg(igs(1:cns))
    ctmax(1:cns)=max(p_dg(1:cns)%tempc,p_rodum(1:cns)%temp(1),p_walledum(1:cns)%temp(1), &
                     p_wallwdum(1:cns)%temp(1))+5. ! max road snow temp
    ctmin(1:cns)=min(p_dg(1:cns)%tempc,p_rodum(1:cns)%temp(1),p_walledum(1:cns)%temp(1), &
                     p_wallwdum(1:cns)%temp(1))-5. ! min road snow temp
    p_sntemp(1:cns)=ctmax(1:cns)
    n(1:cns)=p_rodum(1:cns)%snow*waterden/p_rodum(1:cns)%den ! snow depth
    a(1:cns)=icelambda*(p_rodum(1:cns)%den/waterden)**1.88   ! snow lambda
    p_netldratio(1:cns)=0.5*(n(1:cns)/a(1:cns)+roaddepth(1)/roadlambda(1))
    call solverdsn(cns,evctx(1:cns),p_rg(1:cns),p_fg(1:cns),p_fgtop(1:cns),p_eg(1:cns), &
           p_gasn(1:cns),p_snmelt(1:cns),p_sntemp(1:cns),p_rodum(1:cns),p_walledum(1:cns), &
           p_wallwdum(1:cns),p_dg(1:cns),p_sg(1:cns),p_atm(1:cns),p_netldratio(1:cns),ddt, &
           p_acond(1:cns),p_wallpsi(1:cns),p_roadpsi(1:cns),p_fn(1:cns),p_pg(1:cns))
    p_sntemp(1:cns)=0.5*(ctmax(1:cns)+ctmin(1:cns))
    call solverdsn(cns,evct(1:cns),p_rg(1:cns),p_fg(1:cns),p_fgtop(1:cns),p_eg(1:cns), &
           p_gasn(1:cns),p_snmelt(1:cns),p_sntemp(1:cns),p_rodum(1:cns),p_walledum(1:cns), &
           p_wallwdum(1:cns),p_dg(1:cns),p_sg(1:cns),p_atm(1:cns),p_netldratio(1:cns),ddt, &
           p_acond(1:cns),p_wallpsi(1:cns),p_roadpsi(1:cns),p_fn(1:cns),p_pg(1:cns))
    where (evct(1:cns)*evctx(1:cns).lt.0.)
      ctmin(1:cns)=p_sntemp(1:cns)
    elsewhere
      ctmax(1:cns)=p_sntemp(1:cns)
    end where
    oldval(1:cns)=p_sntemp(1:cns)
    p_sntemp(1:cns)=0.5*(ctmax(1:cns)+ctmin(1:cns))
    do k=1,5 ! sectant
      evctx(1:cns)=evct(1:cns)
      call solverdsn(cns,evct(1:cns),p_rg(1:cns),p_fg(1:cns),p_fgtop(1:cns),p_eg(1:cns), &
             p_gasn(1:cns),p_snmelt(1:cns),p_sntemp(1:cns),p_rodum(1:cns),p_walledum(1:cns), &
             p_wallwdum(1:cns),p_dg(1:cns),p_sg(1:cns),p_atm(1:cns),p_netldratio(1:cns),ddt, &
             p_acond(1:cns),p_wallpsi(1:cns),p_roadpsi(1:cns),p_fn(1:cns),p_pg(1:cns))
      evctx(1:cns)=evct(1:cns)-evctx(1:cns)
      if (all(evctx(1:cns).eq.0.)) exit
      where (evctx(1:cns).ne.0.)
        newval(1:cns)=p_sntemp(1:cns)-evct(1:cns)*(p_sntemp(1:cns)-oldval(1:cns))/evctx(1:cns)
        oldval(1:cns)=p_sntemp(1:cns)
        p_sntemp(1:cns)=newval(1:cns)
      end where
    end do    
    p_sntemp(1:cns)=min(max(p_sntemp(1:cns),ctmin(1:cns)),ctmax(1:cns))
    rg(igs(1:cns))=p_rg(1:cns) ! unpack
    fg(igs(1:cns))=p_fg(1:cns)
    fgtop(igs(1:cns))=p_fgtop(1:cns)
    eg(igs(1:cns))=p_eg(1:cns)
    gardsn(igs(1:cns))=p_gasn(1:cns)
    rdsnmelt(igs(1:cns))=p_snmelt(1:cns)
    rdsntemp(igs(1:cns))=p_sntemp(1:cns)
    acond(igs(1:cns))=p_acond(1:cns)
    dg(igs(1:cns))%canyonrgout=p_dg(1:cns)%canyonrgout
  end if
  if (cnr.gt.0) then ! no road snow
    p_rodum(1:cnr)=roaddum(igr(1:cnr)) ! pack
    p_walledum(1:cnr)=walledum(igr(1:cnr))
    p_wallwdum(1:cnr)=wallwdum(igr(1:cnr))
    p_dg(1:cnr)=dg(igr(1:cnr))
    p_sg(1:cnr)=sg(igr(1:cnr))
    p_atm(1:cnr)=atm(igr(1:cnr))
    p_acond(1:cnr)=acond(igr(1:cnr))
    p_wallpsi(1:cnr)=wallpsi(igr(1:cnr))
    p_roadpsi(1:cnr)=roadpsi(igr(1:cnr))
    p_sntemp(1:cnr)=p_rodum(1:cnr)%temp(1)
    p_netldratio(1:cnr)=0.5*roaddepth(1)/roadlambda(1)
    p_fn(1:cnr)=fn(igr(1:cnr))
    p_pg(1:cnr)=pg(igr(1:cnr))
    call solverdsn(cnr,evct(1:cnr),p_rg(1:cnr),p_fg(1:cnr),p_fgtop(1:cnr),p_eg(1:cnr), &
           p_gasn(1:cnr),p_snmelt(1:cnr),p_sntemp(1:cnr),p_rodum(1:cnr),p_walledum(1:cnr), &
           p_wallwdum(1:cnr),p_dg(1:cnr),p_sg(1:cnr),p_atm(1:cnr),p_netldratio(1:cnr),ddt, &
           p_acond(1:cnr),p_wallpsi(1:cnr),p_roadpsi(1:cnr),p_fn(1:cnr),p_pg(1:cnr))
    rg(igr(1:cnr))=p_rg(1:cnr) ! unpack
    fg(igr(1:cnr))=p_fg(1:cnr)
    fgtop(igr(1:cnr))=p_fgtop(1:cnr)
    eg(igr(1:cnr))=p_eg(1:cnr)
    gardsn(igs(1:cns))=0.
    rdsnmelt(igr(1:cnr))=0.
    rdsntemp(igr(1:cnr))=roaddum(igr(1:cnr))%temp(1)
    acond(igr(1:cnr))=p_acond(1:cnr)
    dg(igr(1:cnr))%canyonrgout=p_dg(1:cnr)%canyonrgout
  end if
  ! ---------------------------------------------------------------    

  ! solve for roof snow temperature -------------------------------
  ! includes solution to longwave radiation
  cns=0
  cnr=0
  rg(igr(1:cnr))%rfsn=0.
  fg(igr(1:cnr))%rfsn=0.
  eg(igr(1:cnr))%rfsn=0.
  garfsn(igr(1:cnr))=0.
  rfsnmelt(igr(1:cnr))=0.
  rfsntemp(igr(1:cnr))=roaddum(igr(1:cnr))%temp(1)
  acond(igs(1:cns))%rfsn=0.
  do iqu=1,ufull ! prepare index arrays
    if (roofdum(iqu)%snow.gt.0.) then
      cns=cns+1
      igs(cns)=iqu
    else
      cnr=cnr+1
      igr(cnr)=iqu
    end if
  end do
  if (cns.gt.0) then ! roof snow
    p_rodum(1:cns)=roofdum(igs(1:cns)) ! pack
    p_dg(1:cns)=dg(igs(1:cns))
    p_sg(1:cns)=sg(igs(1:cns))
    p_atm(1:cns)=atm(igs(1:cns))
    p_acond(1:cns)=acond(igr(1:cns))
    ctmax(1:cns)=max(p_dg(1:cns)%tempr,p_rodum(1:cns)%temp(1))+5. ! max roof snow temp
    ctmin(1:cns)=min(p_dg(1:cns)%tempr,p_rodum(1:cns)%temp(1))-5. ! min roof snow temp
    p_sntemp(1:cns)=ctmax(1:cns)
    n(1:cns)=p_rodum(1:cns)%snow*waterden/p_rodum(1:cns)%den ! snow depth
    a(1:cns)=icelambda*(p_rodum(1:cns)%den/waterden)**1.88   ! snow lambda
    p_netldratio(1:cns)=0.5*(n(1:cns)/a(1:cns)+roofdepth(1)/rooflambda(1))
    call solverfsn(cns,evctx(1:cns),p_rg(1:cns),p_fg(1:cns),p_eg(1:cns),p_gasn(1:cns),p_snmelt(1:cns),p_sntemp(1:cns), &
                   p_rodum(1:cns),p_dg(1:cns),p_sg(1:cns),p_atm(1:cns),p_netldratio(1:cns),p_acond(1:cns),ddt)
    p_sntemp(1:cns)=0.5*(ctmax(1:cns)+ctmin(1:cns))
    call solverfsn(cns,evct(1:cns),p_rg(1:cns),p_fg(1:cns),p_eg(1:cns),p_gasn(1:cns),p_snmelt(1:cns),p_sntemp(1:cns), &
                   p_rodum(1:cns),p_dg(1:cns),p_sg(1:cns),p_atm(1:cns),p_netldratio(1:cns),p_acond(1:cns),ddt)
    where (evct(1:cns)*evctx(1:cns).lt.0.)
      ctmin(1:cns)=p_sntemp(1:cns)
    elsewhere
      ctmax(1:cns)=p_sntemp(1:cns)
    end where
    oldval(1:cns)=p_sntemp(1:cns)
    p_sntemp(1:cns)=0.5*(ctmax(1:cns)+ctmin(1:cns))
    do k=1,5 ! sectant
      evctx(1:cns)=evct(1:cns)
      call solverfsn(cns,evct(1:cns),p_rg(1:cns),p_fg(1:cns),p_eg(1:cns),p_gasn(1:cns),p_snmelt(1:cns),p_sntemp(1:cns), &
                     p_rodum(1:cns),p_dg(1:cns),p_sg(1:cns),p_atm(1:cns),p_netldratio(1:cns),p_acond(1:cns),ddt)
      evctx(1:cns)=evct(1:cns)-evctx(1:cns)
      if (all(evctx(1:cns).eq.0.)) exit
      where (evctx(1:cns).ne.0.)
        newval(1:cns)=p_sntemp(1:cns)-evct(1:cns)*(p_sntemp(1:cns)-oldval(1:cns))/evctx(1:cns)
        oldval(1:cns)=p_sntemp(1:cns)
        p_sntemp(1:cns)=newval(1:cns)
      end where
    end do    
    p_sntemp(1:cns)=min(max(p_sntemp(1:cns),ctmin(1:cns)),ctmax(1:cns))
    rg(igs(1:cns))%rfsn=p_rg(1:cns)%rfsn ! unpack
    fg(igs(1:cns))%rfsn=p_fg(1:cns)%rfsn
    eg(igs(1:cns))%rfsn=p_eg(1:cns)%rfsn
    garfsn(igs(1:cns))=p_gasn(1:cns)
    rfsnmelt(igs(1:cns))=p_snmelt(1:cns)
    rfsntemp(igs(1:cns))=p_sntemp(1:cns)
    acond(igs(1:cns))%rfsn=p_acond(1:cns)%rfsn
  end if
  !---------------------------------------------------------------- 

  ! calculate roof sensible and latent heat fluxes (without snow)
  rg%roof=roofemiss*(atm%rg-sbconst*roofdum%temp(1)**4)
  dg%roofrgout=sbconst*(dg%rfsndelta*snowemiss*rfsntemp**4+(1.-dg%rfsndelta)*roofemiss*roof%temp(1)**4) &
               +(1.-dg%rfsndelta*snowemiss-(1.-dg%rfsndelta)*roofemiss)*atm%rg
  a=log(dg%rfdzmin/zoroof)
  ! n and xe are dummy variables for cd and lzohroof
  call getinvres(ufull,acond%roof,n,xe,a,dg%rfdzmin,roofdum%temp(1),dg%tempr,atm%umag,3)
  fg%roof=aircp*atm%rho*(roofdum%temp(1)-dg%tempr)*acond%roof
  call getqsat(ufull,qsatr,roofdum%temp(1),dg%sigr)
  where (qsatr.lt.dg%mixrr)
    eg%roof=lv*atm%rho*(qsatr-dg%mixrr)*acond%roof
  elsewhere
    ! MJT suggestion for max eg   
    eg%roof=lv*min(atm%rho*dg%roofdelta*(qsatr-dg%mixrr)*acond%roof,roofdum%water+atm%rnd+rfsnmelt)
  end where

  ! calculate change in urban temperatures
  croof%temp(1)=((1.-dg%rfsndelta)*(sg%roof+rg%roof-fg%roof-eg%roof-garoof(:,1)) &
              +dg%rfsndelta*(garfsn-garoof(:,1)))/(roofcp(1)*roofdepth(1))
  cwalle%temp(1)=(sg%walle+rg%walle-fg%walle-eg%walle-gawalle(:,1))/(wallcp(1)*walldepth(1))
  cwallw%temp(1)=(sg%wallw+rg%wallw-fg%wallw-eg%wallw-gawallw(:,1))/(wallcp(1)*walldepth(1))
  croad%temp(1)=((1.-dg%rdsndelta)*(sg%road+rg%road-fg%road-eg%road-garoad(:,1)) &
              +dg%rdsndelta*(gardsn-garoad(:,1)))/(roadcp(1)*roaddepth(1))
  do k=2,3
    croof%temp(k)=(garoof(:,k-1)-garoof(:,k))/(roofcp(k)*roofdepth(k))
    cwalle%temp(k)=(gawalle(:,k-1)-gawalle(:,k))/(wallcp(k)*walldepth(k))
    cwallw%temp(k)=(gawallw(:,k-1)-gawallw(:,k))/(wallcp(k)*walldepth(k))
    croad%temp(k)=(garoad(:,k-1)-garoad(:,k))/(roadcp(k)*roaddepth(k))
  end do
      
  ! calculate change in urban water
  croof%water=atm%rnd-eg%roof/lv+rfsnmelt
  croad%water=atm%rnd-eg%road/lv+rdsnmelt
    
  ! calculate change in urban snow
  croof%snow=atm%snd-eg%rfsn/ls-rfsnmelt
  croad%snow=atm%snd-eg%rdsn/ls-rdsnmelt
  croof%den=(0.24/86400.)*(maxsnowden-roofdum%den)
  croad%den=(0.24/86400.)*(maxsnowden-roaddum%den)
  where (rfsntemp.lt.273.16)
    croof%alpha=-0.008/86400.
  elsewhere
    croof%alpha=(0.24/86400.)*(minsnowalpha-roofdum%alpha)
  end where
  where (rdsntemp.lt.273.16)
    croad%alpha=-0.008/86400.
  elsewhere
    croad%alpha=(0.24/86400.)*(minsnowalpha-roaddum%alpha)
  end where

  ! predictor-corrector scheme to update urban state arrays (e.g., temperature, water and snow)
  if (j.eq.1) then ! predictor
    if (firstcall) then
      roofadj=croof
      walleadj=cwalle
      wallwadj=cwallw
      roadadj=croad
      firstcall=.false.
    end if
    do ii=1,3
      roofdum%temp(ii)=roof%temp(ii)+ddt*(1.5*croof%temp(ii)-0.5*roofadj%temp(ii))
      walledum%temp(ii)=walle%temp(ii)+ddt*(1.5*cwalle%temp(ii)-0.5*walleadj%temp(ii))
      wallwdum%temp(ii)=wallw%temp(ii)+ddt*(1.5*cwallw%temp(ii)-0.5*wallwadj%temp(ii))
      roaddum%temp(ii)=road%temp(ii)+ddt*(1.5*croad%temp(ii)-0.5*roadadj%temp(ii))
    end do
    roofdum%water=roof%water+ddt*(1.5*croof%water-0.5*roofadj%water)
    roaddum%water=road%water+ddt*(1.5*croad%water-0.5*roadadj%water)
    roofdum%snow=roof%snow+ddt*(1.5*croof%snow-0.5*roofadj%snow)
    roaddum%snow=road%snow+ddt*(1.5*croad%snow-0.5*roadadj%snow)
    roofdum%den=roof%den+ddt*(1.5*croof%den-0.5*roofadj%den)
    roaddum%den=road%den+ddt*(1.5*croad%den-0.5*roadadj%den)
    roofdum%alpha=roof%alpha+ddt*(1.5*croof%alpha-0.5*roofadj%alpha)
    roaddum%alpha=road%alpha+ddt*(1.5*croad%alpha-0.5*roadadj%alpha)
    do ii=1,3
      rooforg%temp(ii)=croof%temp(ii)
      walleorg%temp(ii)=cwalle%temp(ii)
      wallworg%temp(ii)=cwallw%temp(ii)
      roadorg%temp(ii)=croad%temp(ii)
    end do
    rooforg%water=croof%water
    roadorg%water=croad%water
    rooforg%snow=croof%snow
    roadorg%snow=croad%snow
    rooforg%den=croof%den
    roadorg%den=croad%den
    rooforg%alpha=croof%alpha
    roadorg%alpha=croad%alpha
  else ! corrector
    do ii=1,3
      roofdum%temp(ii)=roof%temp(ii)+(ddt/12.)*(5.*croof%temp(ii)+8.*rooforg%temp(ii)-roofadj%temp(ii))
      walledum%temp(ii)=walle%temp(ii)+(ddt/12.)*(5.*cwalle%temp(ii)+8.*walleorg%temp(ii)-walleadj%temp(ii))
      wallwdum%temp(ii)=wallw%temp(ii)+(ddt/12.)*(5.*cwallw%temp(ii)+8.*wallworg%temp(ii)-wallwadj%temp(ii))
      roaddum%temp(ii)=road%temp(ii)+(ddt/12.)*(5.*croad%temp(ii)+8.*roadorg%temp(ii)-roadadj%temp(ii))
    end do
    roofdum%water=roof%water+(ddt/12.)*(5.*croof%water+8.*rooforg%water-roofadj%water)
    roaddum%water=road%water+(ddt/12.)*(5.*croad%water+8.*roadorg%water-roadadj%water)
    roofdum%snow=roof%snow+(ddt/12.)*(5.*croof%snow+8.*rooforg%snow-roofadj%snow)
    roaddum%snow=road%snow+(ddt/12.)*(5.*croad%snow+8.*roadorg%snow-roadadj%snow)
    roofdum%den=roof%den+(ddt/12.)*(5.*croof%den+8.*rooforg%den-roofadj%den)
    roaddum%den=road%den+(ddt/12.)*(5.*croad%den+8.*roadorg%den-roadadj%den)
    roofdum%alpha=roof%alpha+(ddt/12.)*(5.*croof%alpha+8.*rooforg%alpha-roofadj%alpha)
    roaddum%alpha=road%alpha+(ddt/12.)*(5.*croad%alpha+8.*roadorg%alpha-roadadj%alpha)
    do ii=1,3
      roofadj%temp(ii)=rooforg%temp(ii)
      walleadj%temp(ii)=walleorg%temp(ii)
      wallwadj%temp(ii)=wallworg%temp(ii)
      roadadj%temp(ii)=roadorg%temp(ii)
    end do
    roofadj%water=rooforg%water
    roadadj%water=roadorg%water
    roofadj%snow=rooforg%snow
    roadadj%snow=roadorg%snow
    roofadj%den=rooforg%den
    roadadj%den=roadorg%den
    roofadj%alpha=rooforg%alpha
    roadadj%alpha=roadorg%alpha
  end if

end do

! limit temperatures to sensible values
do ii=1,3
  roof%temp(ii)=min(max(roofdum%temp(ii),200.),400.)
  walle%temp(ii)=min(max(walledum%temp(ii),200.),400.)
  wallw%temp(ii)=min(max(wallwdum%temp(ii),200.),400.)
  road%temp(ii)=min(max(roaddum%temp(ii),200.),400.)
end do
roof%water=min(max(roofdum%water,0.),maxroofwater)
road%water=min(max(roaddum%water,0.),maxroadwater)
roof%snow=min(max(roofdum%snow,0.),maxrfsn)
road%snow=min(max(roaddum%snow,0.),maxrdsn)
roof%den=min(max(roofdum%den,minsnowden),maxsnowden)
road%den=min(max(roaddum%den,minsnowden),maxsnowden)
roof%alpha=min(max(roofdum%alpha,minsnowalpha),maxsnowalpha)
road%alpha=min(max(roaddum%alpha,minsnowalpha),maxsnowalpha)

! combine snow and snow free tiles
fg%roof=dg%rfsndelta*fg%rfsn+(1.-dg%rfsndelta)*fg%roof ! redefine as net roof fg
eg%roof=dg%rfsndelta*eg%rfsn+(1.-dg%rfsndelta)*eg%roof ! redefine as net roof eg
egtop=dg%rdsndelta*eg%rdsn+(1.-dg%rdsndelta)*eg%road

! calculate outputs
uo%fg=fn%sigmabld*fg%roof+(1.-fn%sigmabld)*fgtop+fn%industryfg
uo%eg=fn%sigmabld*eg%roof+(1.-fn%sigmabld)*egtop
uo%wf=fn%sigmabld*dg%roofdelta*(1.-dg%rfsndelta)+(1.-fn%sigmabld)*dg%roaddelta*(1.-dg%rdsndelta)
uo%ts=((fn%sigmabld*dg%roofrgout+(1.-fn%sigmabld)*dg%canyonrgout)/sbconst)**0.25

! calculate roughness length
call getinvres(ufull,a,pg%cduv,pg%lzoh,pg%lzom,pg%cndzmin,uo%ts,dg%tempc,atm%umag,zohmeth+1)
  
return
end subroutine tebeval

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Estimate saturation mixing ratio (from CCAM)

subroutine getqsat(ifull,qsat,temp,ps)

implicit none

integer, intent(in) :: ifull
real, dimension(ifull), intent(in) :: temp,ps
real, dimension(ifull), intent(out) :: qsat
real, dimension(ifull) :: esatf,tdiff
real, dimension(0:220), parameter :: table = (/ &
1.e-9, 1.e-9, 2.e-9, 3.e-9, 4.e-9, &
6.e-9, 9.e-9, 13.e-9, 18.e-9, 26.e-9, &
36.e-9, 51.e-9, 71.e-9, 99.e-9, 136.e-9, &
0.000000188, 0.000000258, 0.000000352, 0.000000479, 0.000000648, &
0.000000874, 0.000001173, 0.000001569, 0.000002090, 0.000002774, &
0.000003667, 0.000004831, 0.000006340, 0.000008292, 0.00001081, &
0.00001404, 0.00001817, 0.00002345, 0.00003016, 0.00003866, &
0.00004942, 0.00006297, 0.00008001, 0.0001014, 0.0001280, &
0.0001613, 0.0002026, 0.0002538, 0.0003170, 0.0003951, &
0.0004910, 0.0006087, 0.0007528, 0.0009287, 0.001143, &
.001403, .001719, .002101, .002561, .003117, .003784, &
.004584, .005542, .006685, .008049, .009672,.01160,.01388,.01658, &
.01977, .02353, .02796,.03316,.03925,.04638,.05472,.06444,.07577, &
.08894, .1042, .1220, .1425, .1662, .1936, .2252, .2615, .3032, &
.3511, .4060, .4688, .5406, .6225, .7159, .8223, .9432, 1.080, &
1.236, 1.413, 1.612, 1.838, 2.092, 2.380, 2.703, 3.067, 3.476, &
3.935,4.449, 5.026, 5.671, 6.393, 7.198, 8.097, 9.098, &
10.21, 11.45, 12.83, 14.36, 16.06, 17.94, 20.02, 22.33, 24.88, &
27.69, 30.79, 34.21, 37.98, 42.13, 46.69,51.70,57.20,63.23,69.85, &
77.09, 85.02, 93.70, 103.20, 114.66, 127.20, 140.81, 155.67, &
171.69, 189.03, 207.76, 227.96 , 249.67, 272.98, 298.00, 324.78, &
353.41, 383.98, 416.48, 451.05, 487.69, 526.51, 567.52, 610.78, &
656.62, 705.47, 757.53, 812.94, 871.92, 934.65, 1001.3, 1072.2, &
1147.4, 1227.2, 1311.9, 1401.7, 1496.9, 1597.7, 1704.4, 1817.3, &
1936.7, 2063.0, 2196.4, 2337.3, 2486.1, 2643.0, 2808.6, 2983.1, &
3167.1, 3360.8, 3564.9, 3779.6, 4005.5, 4243.0, 4492.7, 4755.1, &
5030.7, 5320.0, 5623.6, 5942.2, 6276.2, 6626.4, 6993.4, 7377.7, &
7780.2, 8201.5, 8642.3, 9103.4, 9585.5, 10089.0, 10616.0, &
11166.0, 11740.0, 12340.0, 12965.0, 13617.0, 14298.0, 15007.0, &
15746.0, 16516.0, 17318.0, 18153.0, 19022.0, 19926.0, 20867.0, &
21845.0, 22861.0, 23918.0, 25016.0, 26156.0, 27340.0, 28570.0, &
29845.0, 31169.0 /)

tdiff=min(max( temp-123.16, 0.), 219.)
esatf=(1.-(tdiff-aint(tdiff)))*table(int(tdiff))+ (tdiff-aint(tdiff))*table(int(tdiff)+1)
qsat=.622*esatf/(ps-esatf)

return
end subroutine getqsat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Originally based on CCAM sflux.f version (i.e., vegetation)
! Modified for increased ratio between momentum and heat roughness
! lengths over urban areas using Brutsaet (1982) parameterisation.

subroutine getinvres(cn,invres,cd,olzoh,ilzom,zmin,stemp,theta,umag,mode)

implicit none

integer, intent(in) :: cn,mode
real, dimension(cn), intent(in) :: ilzom,zmin,stemp,theta,umag
real, dimension(cn), intent(out) :: invres,cd
real, dimension(cn), intent(inout) :: olzoh
real, dimension(cn) :: af,aft,ri,fm,fh,root,denma,denha,re,lna
real, parameter :: bprm=5. ! 4.7 in rams
real, parameter :: chs=2.6 ! 5.3 in rams
real, parameter :: cms=5.  ! 7.4 in rams
real, parameter :: fmroot=0.57735
real, parameter :: rimax=(1./fmroot-1.)/bprm
real, parameter :: nu = 1.461E-5
!real, parameter :: eta0 = 1.827E-5
!real, parameter :: t0 = 291.15
!real, parameter :: c = 120.
!eta=eta0*((t0+c)/(theta+c))*(theta/t0)**(2./3.)
!nu=eta/rho

af=(vkar/ilzom)**2
! umag is now constrained to be above umin in tebcalc
ri=min(grav*zmin*(1.-stemp/theta)/umag**2,rimax)

where (ri>0.)
  fm=1./(1.+bprm*ri)**2
elsewhere
  root=sqrt(-ri*exp(ilzom))
  denma=1.+cms*2.*bprm*af*root
  fm=1.-2.*bprm *ri/denma
end where

cd=af*fm
select case(mode)
  case(1) ! standard veg
    lna=2.
  case(2) ! Kanda et al 2007
    re=max(sqrt(cd)*umag*zmin*exp(-ilzom)/nu,10.)
    !lna=2.46*re**0.25-2. !(Brutsaet, 1982)
    lna=1.29*re**0.25-2. !(Kanda et al, 2007)
  case(3) ! Kanda et al 2005
    lna=6.
  case(4) ! User defined
    lna=olzoh-ilzom
end select
olzoh=lna+ilzom    
aft=vkar**2/(ilzom*olzoh)

where (ri>0.)
  fh=fm
elsewhere
  denha=1.+chs*2.*bprm*aft*exp(0.5*lna)*root
  fh=1.-2.*bprm*ri/denha
end where

invres=aft*fh*umag

return
end subroutine getinvres

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

wallpsi=0.5*(ifn%hwratio+1.-sqrt(ifn%hwratio**2+1.))/ifn%hwratio
roadpsi=sqrt(ifn%hwratio**2+1.)-ifn%hwratio

! integrate through 180 instead of 360
where (ifn%vangle.ge.0.5*pi)
  walles=0.
  wallws=1./ifn%hwratio
  roads=0.
elsewhere
  ta=tan(ifn%vangle)
  thetazero=asin(1./max(ifn%hwratio*ta,1.))
  tc=2.*(1.-cos(thetazero))
  xa=max(ifn%hangle-thetazero,0.)-max(ifn%hangle-pi+thetazero,0.)-min(ifn%hangle+thetazero,0.)
  ya=cos(min(thetazero,max(ifn%hangle-pi,0.)))-cos(min(thetazero,abs(ifn%hangle))) &
     +cos(pi-thetazero)-cos(min(pi,max(ifn%hangle,pi-thetazero)))
  ! note that these terms now include the azimuth angle
  walles=(xa/ifn%hwratio+ta*ya)/pi
  wallws=((pi-2.*thetazero-xa)/ifn%hwratio+ta*(tc-ya))/pi
  roads=(2.*thetazero-ifn%hwratio*ta*tc)/pi
end where

! note that these terms are truncated to nrefl order reflections, compared to TEB which uses infinite reflections.
roadnetalpha=(1.-rdsndelta)*ifn%roadalpha+rdsndelta*rdsnalpha
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

return
end subroutine getswcoeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solve for roof snow temperature

subroutine solverfsn(cn,evct,rg,fg,eg,garfsn,rfsnmelt,rfsntemp,iroof,dg,sg,atm,ldratio,acond,ddt)

implicit none

integer, intent(in) :: cn
real, intent(in) :: ddt
real, dimension(cn), intent(out) :: evct,rfsnmelt,garfsn
real, dimension(cn), intent(in) :: rfsntemp,ldratio
real, dimension(cn) :: cd,rfsnqsat,lzotdum,lzosnow
type(trad), dimension(cn), intent(inout) :: rg,fg,eg,acond
type(trad), dimension(cn), intent(in) :: sg
type(tatm), dimension(cn), intent(in) :: atm
type(tdiag), dimension(cn), intent(in) :: dg
type(tsurf), dimension(cn), intent(in) :: iroof

lzosnow=log(dg%rfdzmin/zosnow)
call getinvres(cn,acond%rfsn,cd,lzotdum,lzosnow,dg%rfdzmin,rfsntemp,dg%tempr,atm%umag,1)
call getqsat(cn,rfsnqsat,rfsntemp,dg%sigr)
rfsnmelt=dg%rfsndelta*max(0.,rfsntemp-273.16)/(icecp*iroof%den*lf*ddt) 
rg%rfsn=snowemiss*(atm%rg-sbconst*rfsntemp**4)
fg%rfsn=aircp*atm%rho*(rfsntemp-dg%tempr)*acond%rfsn
eg%rfsn=ls*min(atm%rho*dg%rfsndelta*max(0.,rfsnqsat-dg%mixrr)*acond%rfsn,iroof%snow+atm%snd-rfsnmelt) ! MJT suggestion for max eg
garfsn=(rfsntemp-iroof%temp(1))/ldratio
evct=sg%rfsn+rg%rfsn-fg%rfsn-eg%rfsn-garfsn

return
end subroutine solverfsn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solve for road snow temperature (includes canyon temperature)

subroutine solverdsn(cn,evct,rg,fg,fgtop,eg,gardsn,rdsnmelt,rdsntemp,iroad,iwalle,iwallw,dg,sg,atm &
                    ,ldratio,ddt,acond,wallpsi,roadpsi,ifn,ipg)

implicit none

integer, intent(in) :: cn
integer k
real, intent(in) :: ddt
real, dimension(cn), intent(out) :: evct,fgtop,gardsn,rdsnmelt
real, dimension(cn), intent(in) :: rdsntemp,ldratio,wallpsi,roadpsi
real, dimension(cn) :: ctmax,ctmin,cevctx,cevct,oldval,newval,roadqsat,rdsnqsat,canyonmix,netrad,netemiss
real, dimension(cn) :: cwa,cwe,cww,cwr,cra,crr,crw,ncwa,ncwe,ncww,ncwr,ncra,ncrr,ncrw
real, dimension(cn) :: rcwa,rcwe,rcww,rcwr,rcra,rcrr,rcrw
real, dimension(cn) :: canyontemp,topinvres
type(trad), dimension(cn), intent(inout) :: rg,fg,eg,acond
type(trad), dimension(cn), intent(in) :: sg
type(tatm), dimension(cn), intent(in) :: atm
type(tdiag), dimension(cn), intent(inout) :: dg
type(tsurf), dimension(cn), intent(in) :: iroad
type(twall), dimension(cn), intent(in) :: iwalle,iwallw
type(tdata), dimension(cn), intent(in) :: ifn
type(tprog), dimension(cn), intent(in) :: ipg

! solve for canyon temperature ----------------------------------
ctmax=max(dg%tempc,iwalle%temp(1),iwallw%temp(1),iroad%temp(1),rdsntemp)+5. ! max canyon temp
ctmin=min(dg%tempc,iwalle%temp(1),iwallw%temp(1),iroad%temp(1),rdsntemp)-5. ! min canyon temp
call solvecanyon(cn,cevctx,fg,fgtop,topinvres,ctmax,dg,atm,iwalle%temp(1),iwallw%temp(1),iroad%temp(1) &
                ,rdsntemp,acond,ifn,ipg)
canyontemp=0.5*(ctmax+ctmin)
call solvecanyon(cn,cevct,fg,fgtop,topinvres,canyontemp,dg,atm,iwalle%temp(1),iwallw%temp(1),iroad%temp(1) &
                ,rdsntemp,acond,ifn,ipg)
where ((cevct*cevctx).lt.0.)
  ctmin=canyontemp
elsewhere
  ctmax=canyontemp
end where
oldval=canyontemp
canyontemp=0.5*(ctmax+ctmin)
do k=1,5 ! sectant
  cevctx=cevct
  call solvecanyon(cn,cevct,fg,fgtop,topinvres,canyontemp,dg,atm,iwalle%temp(1),iwallw%temp(1),iroad%temp(1) &
                  ,rdsntemp,acond,ifn,ipg)
  cevctx=cevct-cevctx
  if (all(cevctx.eq.0.)) exit
  where (cevctx.ne.0.)
    newval=canyontemp-cevct*(canyontemp-oldval)/cevctx
    oldval=canyontemp
    canyontemp=newval
  end where
end do
canyontemp=min(max(canyontemp,ctmin),ctmax)
! ---------------------------------------------------------------    

! note additional nrefl order reflections
netrad=dg%rdsndelta*snowemiss*rdsntemp**4+(1.-dg%rdsndelta)*roademiss*iroad%temp(1)**4
netemiss=dg%rdsndelta*snowemiss+(1.-dg%rdsndelta)*roademiss
cwa=wallpsi
cra=roadpsi
cwe=0.
cww=1.-2.*wallpsi
crw=0.5*(1.-roadpsi)
crr=0.
cwr=wallpsi
rcwa=cwa
rcra=cra
rcwe=cwe
rcww=cww
rcrw=crw
rcrr=crr
rcwr=cwr
do k=1,nrefl
  ncwa=(1.-netemiss)*wallpsi*rcra+(1.-wallemiss)*(1.-2.*wallpsi)*rcwa
  ncra=(1.-wallemiss)*(1.-roadpsi)*rcwa
  ncwe=(1.-netemiss)*wallpsi*rcrw+(1.-wallemiss)*(1.-2.*wallpsi)*rcww
  ncww=(1.-netemiss)*wallpsi*rcrw+(1.-wallemiss)*(1.-2.*wallpsi)*rcwe
  ncrw=(1.-wallemiss)*(1.-roadpsi)*0.5*(rcww+rcwe)  
  ncwr=(1.-netemiss)*wallpsi*rcrr+(1.-wallemiss)*(1.-2.*wallpsi)*rcwr
  ncrr=(1.-wallemiss)*(1.-roadpsi)*rcwr
  rcwa=ncwa
  rcra=ncra
  rcwe=ncwe
  rcww=ncww
  rcrw=ncrw
  rcrr=ncrr
  rcwr=ncwr
  cwa=cwa+rcwa
  cra=cra+rcra
  cwe=cwe+rcwe
  cww=cww+rcww
  crw=crw+rcrw
  cwr=cwr+rcwr
  crr=crr+rcrr
end do
    
rg%walle=wallemiss*(atm%rg*cwa+sbconst*iwalle%temp(1)**4*(-1.+wallemiss*cwe) & 
                  +sbconst*iwallw%temp(1)**4*wallemiss*cww+sbconst*netrad*cwr)
rg%wallw=wallemiss*(atm%rg*cwa+sbconst*iwallw%temp(1)**4*(-1.+wallemiss*cwe) &
                  +sbconst*iwalle%temp(1)**4*wallemiss*cww+sbconst*netrad*cwr)
rg%road=roademiss*(atm%rg*cra+sbconst*(-iroad%temp(1)**4+netrad*crr) &
                  +sbconst*wallemiss*(iwalle%temp(1)**4+iwallw%temp(1)**4)*crw)
where (dg%rdsndelta.gt.0.)
  rg%rdsn=snowemiss*(atm%rg*cra+sbconst*(-rdsntemp**4+netrad*crr) &
                    +sbconst*wallemiss*(iwalle%temp(1)**4+iwallw%temp(1)**4)*crw)
elsewhere
  rg%rdsn=0.
end where
dg%canyonrgout=atm%rg*(2.*ifn%hwratio*wallpsi*(1.-wallemiss)*cwa+roadpsi*(1.-netemiss)*cwr) &
               +sbconst*wallemiss*iwalle%temp(1)**4*(ifn%hwratio*wallpsi*(1.+(1.-wallemiss)*(cwe+cww)) &
               +roadpsi*(1.-netemiss)*crw) &
               +sbconst*wallemiss*iwallw%temp(1)**4*(ifn%hwratio*wallpsi*(1.+(1.-wallemiss)*(cwe+cww)) &
               +roadpsi*(1.-netemiss)*crw) &
               +sbconst*netrad*(2.*ifn%hwratio*wallpsi*(1.-wallemiss)*cwr+roadpsi*(1.+(1.-netemiss)*crr))

rdsnmelt=dg%rdsndelta*max(0.,rdsntemp-273.16)/(icecp*iroad%den*lf*ddt)
call getqsat(cn,roadqsat,iroad%temp(1),dg%sigd) ! evaluate using pressure at displacement height
call getqsat(cn,rdsnqsat,rdsntemp,dg%sigd)
canyonmix=((dg%rdsndelta*rdsnqsat*ls/lv*acond%rdsn+(1.-dg%rdsndelta)*dg%roaddelta*roadqsat*acond%road) &
           +dg%mixrc*topinvres)/((dg%rdsndelta*ls/lv*acond%rdsn+(1.-dg%rdsndelta)*dg%roaddelta*acond%road) &
           +topinvres)
where (roadqsat.lt.canyonmix)
  eg%road=lv*atm%rho*(roadqsat-canyonmix)*acond%road
elsewhere
  eg%road=lv*min(atm%rho*dg%roaddelta*(roadqsat-canyonmix)*acond%road,iroad%water+atm%rnd+rdsnmelt) ! MJT suggestion for max eg
end where
where (dg%rdsndelta.gt.0.)
  eg%rdsn=ls*min(atm%rho*dg%rdsndelta*max(0.,rdsnqsat-canyonmix)*acond%rdsn,iroad%snow+atm%snd-rdsnmelt) ! MJT suggestion for max eg
  gardsn=(rdsntemp-iroad%temp(1))/ldratio
elsewhere
  eg%rdsn=0.
  gardsn=0.
end where
eg%walle=0.
eg%wallw=0.

evct=sg%rdsn+rg%rdsn-fg%rdsn-eg%rdsn-gardsn

return
end subroutine solverdsn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solve for canyon temperature

subroutine solvecanyon(cn,evct,fg,fgtop,topinvres,ctemp,dg,atm,walletemp,wallwtemp,roadtemp,rdsntemp,acond,ifn,ipg)

implicit none

integer, intent(in) :: cn
real, dimension(cn), intent(out) :: evct,fgtop,topinvres
real, dimension(cn), intent(in) :: ctemp,walletemp,wallwtemp,roadtemp,rdsntemp
real, dimension(cn) :: trafficout,cduv,lzoh
type(trad), dimension(cn), intent(inout) :: fg,acond
type(tdiag), dimension(cn), intent(in) :: dg
type(tatm), dimension(cn), intent(in) :: atm
type(tdata), dimension(cn), intent(in) :: ifn
type(tprog), dimension(cn), intent(in) :: ipg

call getinvres(cn,topinvres,cduv,lzoh,ipg%lzom,ipg%cndzmin,ctemp,dg%tempc,atm%umag,1)
if (resmeth.eq.0) then
  !cw=pg%cduv*atm%umag ! diagnose canyonw (from Masson 2000)
  acond%road=(11.8+4.2*sqrt(acond%road**2+ipg%cduv*atm%umag**2))/(aircp*atm%rho) ! From Rowley, et al (1930)
  acond%walle=acond%road
  acond%wallw=acond%road
  acond%rdsn=acond%road
end if
fg%walle=aircp*atm%rho*(walletemp-ctemp)*acond%walle
fg%wallw=aircp*atm%rho*(wallwtemp-ctemp)*acond%wallw
fg%road=aircp*atm%rho*(roadtemp-ctemp)*acond%road
where (dg%rdsndelta.gt.0.)
  fg%rdsn=aircp*atm%rho*(rdsntemp-ctemp)*acond%rdsn
elsewhere
  fg%rdsn=0.
end where
fgtop=aircp*atm%rho*(ctemp-dg%tempc)*topinvres
call gettraffic(cn,trafficout,ifn%trafficfg,ifn%ctime)
evct=fgtop-(dg%rdsndelta*fg%rdsn+(1.-dg%rdsndelta)*fg%road+ifn%hwratio*(fg%walle+fg%wallw) &
           +trafficout/(1.-ifn%sigmabld)+dg%accool)

return
end subroutine solvecanyon

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gettraffic(cn,trafficout,trafficfg,ctime)

implicit none

integer, intent(in) :: cn
integer, dimension(cn) :: ip
real, dimension(cn), intent(out) :: trafficout
real, dimension(cn), intent(in) :: trafficfg,ctime
real, dimension(cn) :: rp
! traffic diurnal cycle weights approximated from Coutts et al (2007)
real, dimension(25), parameter :: trafficcycle = (/ 0.1, 0.1, 0.1, 0.1, 0.5, 1.0, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, &
                                                    1.5, 1.5, 1.5, 1.5, 1.5, 1.4, 1.2, 1., 0.8, 0.6, 0.4, 0.2, 0.1 /) 

ip=int(24.*ctime)
rp=ctime-real(ip)
where (ip.lt.1) ip=ip+24
trafficout=trafficfg*((1.-rp)*trafficcycle(ip)+rp*trafficcycle(ip+1))

return
end subroutine gettraffic

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine tebdisable(diag)

implicit none

integer, intent(in) :: diag

if (diag.ge.1) write(6,*) "Disable aTEB"
ufull=0

return
end subroutine tebdisable

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module ateb
