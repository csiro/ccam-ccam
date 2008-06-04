
! This code was originally based on the TEB scheme of Masson, Boundary-Layer Meteorology, 94, p357 (2000)
! The snow scheme is based on Douville, Royer and Mahfouf, Climate Dynamics, 12, p21 (1995)

! Usual practice is:
!   call tebinit     ! to initalise state arrays, etc (use tebdisable to disable calls to ateb subroutines)
!   call tebloadm    ! to load previous state arrays (from tebsavem)
!   call tebtype     ! to define urban type (use tebfndef to define urban properties at each grid point)
!   ...
!   call tebnewangle ! store solar zenith and azimuthal angle (use tebccangle for CCAM or
!                      use tebnewangle1 for a single grid point)
!   call tebalb      ! blends input and urban albedo (use tebalb1 for a single grid point)
!   call tebcalc     ! calculates urban temperatures, fluxes, etc and blends with input
!   call tebzo       ! blends input and urban roughness lengths for momentum and heat
!   ...
!   call tebsavem    ! to save current state arrays (for use by tebloadm)
!   call tebend      ! to deallocate memory before quiting

! only tebinit and tebcalc are manditory.  All other subroutine calls are optional.

! URBAN TYPES:
 
! 1 = Default urban (ecosystems 007,152,153,154 with differences in sigmau only)
! 2 = Dense urban (ecosystems 151)
! 3 = Industries (ecosystems 155)
! 4 = Road and rail (ecosystems 156)
! 5 = Ports (ecosystems 157)
! 6 = Airport (ecosystems 158)
! 7 = Construction (ecosystems 159)
! 8 = Urban parks (ecosystems 160)
! 9 = Sport facilities (ecosystems 161)

! NOTES: 
!  Below are differences with TEB (Masson 2000) scheme and aTEB:
!
! - Usually the zero model height in aTEB is set to the canyon displacement height, not the building height as in TEB.
!
! - aTEB uses two walls instead of the TEB single wall.  Also, only up to 2nd order reflections are used for
!   longwave and short wave radation.  In TEB, infinite reflections are used for shortwave, but only 2nd order for
!   long wave.  See Harman, et. al (2004) for a complete treatment of the reflections.
!
! - aTEB employs zot=zom*7.4*exp(-2.5*(ustar*zom/nu)**0.25) (Brutsaet, 1982)
!
! - aTEB calculates resistances for the recirculation and ventilation regions of the canyon (see Harman, et. al 2004).
!   This approach takes advantage of the second wall temperature (i.e., the fluxes depend on the wind direction).
!
! - aTEB plans to improve the representation of traffic and industry heat fluxes.
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
! parameters
real, dimension(3), parameter :: roofdepth =(/ 0.05,0.4,0.1 /)          ! depth (m)
real, dimension(3), parameter :: walldepth =(/ 0.02,0.125,0.05 /)
real, dimension(3), parameter :: roaddepth =(/ 0.05,0.1,1. /)
real, dimension(3), parameter :: roofcp =(/ 2.11E6,0.28E6,0.29E6 /)     ! heat capacity (J m^-3 K^-1)
real, dimension(3), parameter :: wallcp =(/ 1.55E6,1.55E6,0.29E6 /)
real, dimension(3), parameter :: roadcp =(/ 1.94E6,1.28E6,1.28E6 /)
real, dimension(3), parameter :: rooflambda =(/ 1.51,0.08,0.05 /)       ! conductance (W m^-1 K^-1)
real, dimension(3), parameter :: walllambda =(/ 0.9338,0.9338,0.05 /)
real, dimension(3), parameter :: roadlambda =(/ 0.7454,0.2513,0.2513 /)
real, parameter :: maxsnowden=300.   ! max snow density (kg m^-3)
real, parameter :: minsnowden=100.   ! min snow density (kg m^-3)
real, parameter :: waterden=1000.    ! water density (kg m^-3)
real, parameter :: icelambda=2.22    ! conductance of ice (W m^-1 K^-1)
real, parameter :: aircp=1004.64     ! Heat capapcity of dry air (J kg^-1 K^-1)
real, parameter :: icecp=2100.       ! Heat capacity of ice (J kg^-1 K^-1)
real, parameter :: bldtemp=291.16    ! Comfort temperature (K) = 18deg C
real, parameter :: grav=9.80616      ! gravity (m s^-2)
real, parameter :: vkar=0.4          ! von Karman constant
real, parameter :: lv=2.501e6        ! Latent heat of vaporisation
real, parameter :: lf=3.337e5        ! Latent heat of fusion
real, parameter :: ls=lv+lf          ! Latent heat of sublimation (2.834e6)
real, parameter :: pi=3.1415927      ! pi
real, parameter :: rd=287.04         ! Gas constant for dry air
real, parameter :: sbconst=5.67e-8   ! Stefan-Boltzmann constant
real, parameter :: zoroof=0.15       ! Roughness length for rooftops (m) (see Masson 2000)
real, parameter :: zosnow=0.001      ! Roughness length for snow (m)
real, parameter :: zocanyon=5.e-5    ! Roughness length for canyon surfaces (m) (see Harman et al 2004)
real, parameter :: roofemiss=0.90    ! emissitivity
real, parameter :: wallemiss=0.85 
real, parameter :: roademiss=0.94
real, parameter :: snowemiss=1.      ! snow emissitivity
real, parameter :: roofalpha=0.15    ! Albedo
real, parameter :: wallalpha=0.25
real, parameter :: roadalpha=0.08
real, parameter :: maxsnowalpha=0.85 ! max snow albedo
real, parameter :: minsnowalpha=0.5  ! min snow albedo
real, parameter :: maxroofwater=1.   ! max water on roof (kg m^-2)
real, parameter :: maxroadwater=1.   ! max water on road (kg m^-2)
real, parameter :: maxrfsn=1.        ! max snow on roof (kg m^-2)
real, parameter :: maxrdsn=1.        ! max snow on road (kg m^-2)
real, parameter :: refheight=2./3.   ! zero model height as fraction of building height
integer, parameter :: resmeth=1      ! Canyon resistances (0=Masson, 1=Harman)

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

ufull=count(sigmau.gt.0.)
if (ufull.eq.0) return

allocate(ugrid(ufull),mgrid(ifull),fn(ufull),pg(ufull))
allocate(roof(ufull),road(ufull),walle(ufull),wallw(ufull))
allocate(roofadj(ufull),roadadj(ufull),walleadj(ufull),wallwadj(ufull))
allocate(sigmau(ufull))

! define grid arrays
mgrid=0
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
  roof(:)%temp(ii)=bldtemp
  road(:)%temp(ii)=bldtemp
  walle(:)%temp(ii)=bldtemp
  wallw(:)%temp(ii)=bldtemp
end do
roof(:)%water=0.
road(:)%water=0.
roof(:)%snow=0.
road(:)%snow=0.
roof(:)%den=minsnowden
road(:)%den=minsnowden
roof(:)%alpha=maxsnowalpha
road(:)%alpha=maxsnowalpha

utype(:)=1                    ! default urban
call tebtype(ifull,utype,diag)

fn(:)%vangle=0.
fn(:)%hangle=0.

pg(:)%cndzmin=max(abs(zmin-fn(:)%bldheight*(2./3.-refheight)),fn(:)%zo+1.) ! updated in tebcalc
pg(:)%lzom=log(pg(:)%cndzmin/fn(:)%zo)                                     ! updated in tebcalc
pg(:)%lzoh=2.+pg(:)%lzom                                                   ! updated in tebcalc

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
  roof(:)%temp(ii)=urban(ugrid(:),ii)
  walle(:)%temp(ii)=urban(ugrid(:),ii+3)
  wallw(:)%temp(ii)=urban(ugrid(:),ii+6)
  road(:)%temp(ii)=urban(ugrid(:),ii+9)
end do
roof(:)%water=urban(ugrid(:),13)
road(:)%water=urban(ugrid(:),14)
roof(:)%snow=urban(ugrid(:),15)
road(:)%snow=urban(ugrid(:),16)
roof(:)%den=urban(ugrid(:),17)
road(:)%den=urban(ugrid(:),18)
roof(:)%alpha=urban(ugrid(:),19)
road(:)%alpha=urban(ugrid(:),20)

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
  roof(:)%temp(ii)=urban(ugrid(:),ii)
  walle(:)%temp(ii)=urban(ugrid(:),ii+3)
  wallw(:)%temp(ii)=urban(ugrid(:),ii+6)
  road(:)%temp(ii)=urban(ugrid(:),ii+9)
end do

return
end subroutine tebloadm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine loads aTEB type arrays (not compulsory)

subroutine tebtype(ifull,itype,diag)

implicit none

integer, intent(in) :: ifull,diag
integer, dimension(ifull), intent(in) :: itype
integer, parameter :: maxtype = 9
real, dimension(maxtype), parameter :: chwratio=(/ 0.21, 0.82, 0.4, 0.05, 0.82, 0.02, 0.005, 0.005, 0.11 /) ! Building height to width ratio
real, dimension(maxtype), parameter :: csigmabld=(/ 0.5, 0.5, 0.5, 0.1, 0.5, 0.1, 0.1, 0.1, 0.5 /)          ! Area fraction occupied by buildings
real, dimension(maxtype), parameter :: cindustryfg=(/ 5., 10., 20., 0., 20., 0., 0., 0., 0. /)              ! Industral sensible heat flux (W m^-2)
real, dimension(maxtype), parameter :: ctrafficfg=(/ 10., 20., 10., 30., 10., 10., 0., 0., 0. /)            ! Traffic sensible heat flux (W m^-2)
real, dimension(maxtype), parameter :: cbldheight=(/ 10., 30., 20., 5., 20., 10., 5., 5., 10. /)            ! Building height (m)
real, dimension(maxtype), parameter :: czo=(/ 1., 3., 2., 0.5, 2., 0.01, 0.1, 0.5, 1. /)                    ! Effective roughness length (m)

if (ufull.eq.0) return
if (diag.ge.1) write(6,*) "Load aTEB building properties"

if (any(itype(:).lt.1).or.any(itype(:).gt.maxtype)) then
  write(6,*) "ERROR: Urban type is out of range"
  stop
end if

fn(:)%hwratio=chwratio(itype(ugrid(:)))
fn(:)%sigmabld=csigmabld(itype(ugrid(:)))
fn(:)%industryfg=cindustryfg(itype(ugrid(:)))
fn(:)%trafficfg=ctrafficfg(itype(ugrid(:)))
fn(:)%bldheight=cbldheight(itype(ugrid(:)))
fn(:)%zo=czo(itype(ugrid(:)))

return
end subroutine tebtype

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine specifies the urban properties for each grid point
!

subroutine tebfndef(ifull,ifn,diag)

implicit none

integer, intent(in) :: ifull,diag
real, dimension(ifull,6), intent(in) :: ifn

if (diag.ge.1) write(6,*) "Load aTEB building properties"

fn(:)%hwratio=ifn(ugrid(:),1)
fn(:)%sigmabld=ifn(ugrid(:),2)
fn(:)%industryfg=ifn(ugrid(:),3)
fn(:)%trafficfg=ifn(ugrid(:),4)
fn(:)%bldheight=ifn(ugrid(:),5)
fn(:)%zo=ifn(ugrid(:),6)

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
  urban(ugrid(:),ii)=roof(:)%temp(ii)
  urban(ugrid(:),ii+3)=walle(:)%temp(ii)
  urban(ugrid(:),ii+6)=wallw(:)%temp(ii)
  urban(ugrid(:),ii+9)=road(:)%temp(ii)
end do
urban(ugrid(:),13)=roof(:)%water
urban(ugrid(:),14)=road(:)%water
urban(ugrid(:),15)=roof(:)%snow
urban(ugrid(:),16)=road(:)%snow
urban(ugrid(:),17)=roof(:)%den
urban(ugrid(:),18)=road(:)%den
urban(ugrid(:),19)=roof(:)%alpha
urban(ugrid(:),20)=road(:)%alpha

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
  urban(ugrid(:),ii)=roof(:)%temp(ii)
  urban(ugrid(:),ii+3)=walle(:)%temp(ii)
  urban(ugrid(:),ii+6)=wallw(:)%temp(ii)
  urban(ugrid(:),ii+9)=road(:)%temp(ii)
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

if (diag.ge.1) write(6,*) "Blend urban roughness lengths"

! evaluate at canyon displacement height (really the atmospheric model should provide a displacement height)
workb(:)=sqrt((1.-sigmau(:))/log(pg(:)%cndzmin/zom(ugrid(:)))**2+sigmau(:)/pg(:)%lzom**2)
workc(:)=(1.-sigmau(:))/(log(pg(:)%cndzmin/zom(ugrid(:)))*log(pg(:)%cndzmin/zoh(ugrid(:)))) &
        +sigmau(:)/(pg(:)%lzom*pg(:)%lzoh)
workc(:)=workc(:)/workb(:)
workb(:)=pg(:)%cndzmin*exp(-1./workb(:))
workc(:)=max(pg(:)%cndzmin*exp(-1./workc(:)),zr)
zom(ugrid(:))=workb(:)
zoh(ugrid(:))=workc(:)

if (any(zoh(ugrid(:)).le.zr)) write(6,*) "WARN: minimum zoh reached"

return
end subroutine tebzo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine blends the urban drag coeff
!

subroutine tebcd(ifull,cduv,diag)

implicit none

integer, intent(in) :: ifull,diag
real, dimension(ifull), intent(inout) :: cduv

cduv(ugrid(:))=(1.-sigmau(:))*cduv(ugrid(:))+sigmau(:)*pg(:)%cduv

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
alb(ugrid(:))=(1.-sigmau(:))*alb(ugrid(:))+sigmau(:)*ualb(:)

return
end subroutine tebalb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates the urban contrabution to albedo.
! (selected grid points only)

subroutine tebalb1(is,ifull,alb,diag)

implicit none

integer, intent(in) :: is,ifull,diag
integer, dimension(:), allocatable :: lgrid
integer i,ucount,ib,ie
real, dimension(ifull), intent(inout) :: alb
real, dimension(:), allocatable :: ualb

if (ufull.eq.0) return

ucount=count(mgrid(is:is+ifull-1).ge.1)
if (ucount.eq.0) return

allocate(lgrid(ucount),ualb(ucount))

ucount=0
do i=1,ifull
  if (mgrid(is+i-1).ge.1) then
    ucount=ucount+1
    lgrid(ucount)=i
  end if
end do

ib=mgrid(is+lgrid(1)-1)
ie=ib+ucount-1

call tebalbcalc(ib,ucount,ualb,diag)
alb(lgrid(:))=(1.-sigmau(ib:ie))*alb(lgrid(:))+sigmau(ib:ie)*ualb(:)

deallocate(lgrid,ualb)

return
end subroutine tebalb1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Albedo calculations

subroutine tebalbcalc(is,ifull,alb,diag)

implicit none

integer, intent(in) :: is,ifull,diag
integer i,iqu
real, dimension(ifull), intent(out) :: alb
real albu,albr,snowdelta
real wallpsi,roadpsi
type(trad) :: sg

do i=1,ifull
  iqu=is+i-1
  snowdelta=road(iqu)%snow/(road(iqu)%snow+maxrdsn)
  call getswcoeff(sg,wallpsi,roadpsi,fn(iqu),road(iqu)%alpha,snowdelta)
  albu=fn(iqu)%hwratio*(sg%walle+sg%wallw)*wallalpha+sg%road*((1.-snowdelta)*roadalpha+snowdelta*road(iqu)%alpha)
  snowdelta=roof(iqu)%snow/(roof(iqu)%snow+maxrfsn)
  albr=(1.-snowdelta)*roofalpha+snowdelta*roof(iqu)%alpha
  alb(i)=fn(iqu)%sigmabld*albr+(1.-fn(iqu)%sigmabld)*albu
end do

return
end subroutine tebalbcalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine stores the zenith angle and the solar azimuth angle
! (all grid points)

! ifull = size of array
! cosin = array of cosine of zenith angles
! azimuthin = array of azimuthal angles
! diag = dialogue flag (0=off)

subroutine tebnewangle(ifull,cosin,azimuthin,diag)

implicit none

integer, intent(in) :: ifull,diag
real, dimension(ifull), intent(in) :: cosin     ! cosine of zenith angle
real, dimension(ifull), intent(in) :: azimuthin ! azimuthal angle

if (ufull.eq.0) return
if (diag.ge.1) write(6,*) "Update solar zenith angle and azimuth angle"

fn(:)%hangle=0.5*pi-azimuthin(ugrid(:))
fn(:)%vangle=acos(cosin(ugrid(:)))  

return
end subroutine tebnewangle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine stores the zenith angle and the solar azimuth angle
! (single grid point)

! iq = grid point
! cosin = cosine of zenith angle
! azmiuthin = azimuthal angle (rad)

subroutine tebnewangle1(is,ifull,cosin,azimuthin)

implicit none

integer, intent(in) :: is,ifull
integer i,iqu
real, dimension(ifull), intent(in) :: cosin     ! cosine of zenith angle
real, dimension(ifull), intent(in) :: azimuthin ! azimuthal angle

if (ufull.eq.0) return

do i=1,ifull
  iqu=mgrid(is+i-1)
  if (iqu.ge.1) then
    fn(iqu)%hangle=0.5*pi-azimuthin(i)
    fn(iqu)%vangle=acos(cosin(i))  
  end if
end do

return
end subroutine tebnewangle1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This version of tebnewangle is for CCAM
!

subroutine tebccangle(is,ifull,cosin,rlon,rlat,fjd,slag,dhr,dlt)

implicit none

integer, intent(in) :: is,ifull
integer i,iq,iqu
real, intent(in) :: fjd,slag,dhr,dlt
real, dimension(ifull), intent(in) :: cosin,rlon,rlat
real hloc,x,y,azimuth

if (ufull.eq.0) return

do i=1,ifull
  iq=i+is-1
  iqu=mgrid(iq)
  if (iqu.ge.1) then
    ! from zenith.f
    hloc = 2.*pi*fjd+slag+pi+rlon(i)+dhr*pi/24.

    ! estimate azimuth angle
    x=sin(-hloc)*cos(dlt)
    y=-cos(-hloc)*cos(dlt)*sin(rlat(i))+cos(rlat(i))*sin(dlt)
    azimuth=atan2(x,y)

    fn(iqu)%hangle=0.5*pi-azimuth
    fn(iqu)%vangle=acos(cosin(i))
  end if
end do

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

subroutine tebcalc(ifull,ofg,oeg,ots,owf,ddt,zmin,sg,rg,rnd,rho,temp,mixr,ps,pa,uu,vv,diag)

implicit none

integer, intent(in) :: ifull,diag
real, intent(in) :: ddt,zmin
real, dimension(ifull), intent(in) :: sg,rg,rnd,rho,temp,mixr,ps,pa,uu,vv
real, dimension(ifull), intent(inout) :: ofg,oeg,ots,owf
type (tatm), dimension(ufull) :: atm
type (tout), dimension(ufull) :: uo

if (ufull.eq.0) return

atm(:)%sg=sg(ugrid(:))
atm(:)%rg=rg(ugrid(:))
atm(:)%rho=rho(ugrid(:))
atm(:)%temp=temp(ugrid(:))
atm(:)%mixr=mixr(ugrid(:))
atm(:)%ps=ps(ugrid(:))
atm(:)%pa=pa(ugrid(:))
atm(:)%umag=sqrt(uu(ugrid(:))**2+vv(ugrid(:))**2)
atm(:)%udir=atan2(vv(ugrid(:)),uu(ugrid(:)))
where (atm(:)%temp.le.273.16)
  atm(:)%snd=rnd(ugrid(:))
  atm(:)%rnd=0.
else where
  atm(:)%rnd=rnd(ugrid(:))
  atm(:)%snd=0.
end where

call tebeval(uo,ddt,atm,zmin,diag)

ofg(ugrid(:))=(1.-sigmau(:))*ofg(ugrid(:))+sigmau(:)*uo(:)%fg
oeg(ugrid(:))=(1.-sigmau(:))*oeg(ugrid(:))+sigmau(:)*uo(:)%eg
ots(ugrid(:))=(1.-sigmau(:))*ots(ugrid(:))+sigmau(:)*uo(:)%ts
owf(ugrid(:))=(1.-sigmau(:))*owf(ugrid(:))+sigmau(:)*uo(:)%wf

end subroutine tebcalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! urban flux calculations

subroutine tebeval(uo,ddt,atm,zmin,diag)

implicit none

integer, intent(in) :: diag
integer iqu,j,k
real, intent(in) :: ddt,zmin
real, dimension(4) :: acond
real, dimension(3) :: garoof,gawalle,gawallw,garoad
real rdsntemp,rfsntemp,rdsnmelt,rfsnmelt
real wallpsi,roadpsi,fgtop,egtop,garfsn,gardsn
real oldtemp,newtemp,canyontemp,topu,cu,ctmax,ctmin,evctx,evct
real cd,roofqsat,roadqsat,qsata,lzotdum,roofinvres
real nettemp,netldratio
real lr,wdt,ln,rn,we,ww,wr,zolog,a,xe,xw,cuven,n
real sndepth,snlambda,zom,lzoroof
logical, save :: firstcall=.true.
type(tatm), dimension(ufull), intent(in) :: atm
type(tout), dimension(ufull), intent(out) :: uo
type(trad) :: sg,rg,fg,eg
type(tdiag) :: dg
type(tsurf) :: roofdum,roaddum,croof,croad,rooforg,roadorg
type(twall) :: walledum,wallwdum,cwalle,cwallw,walleorg,wallworg

if (diag.ge.1) write(6,*) "Evaluating aTEB"

! new snowfall
where (atm(:)%snd.gt.0.)
  roof(:)%den=(roof(:)%snow*roof(:)%den+atm(:)%snd*ddt*minsnowden)/(roof(:)%snow+ddt*atm(:)%snd)
  road(:)%den=(road(:)%snow*road(:)%den+atm(:)%snd*ddt*minsnowden)/(road(:)%snow+ddt*atm(:)%snd)
  roof(:)%alpha=maxsnowalpha
  road(:)%alpha=maxsnowalpha
end where

do iqu=1,ufull

  ! canyon (displacement height at 2./3. building height)
  dg%sigd=atm(iqu)%ps*exp(grav*fn(iqu)%bldheight*(refheight-2./3.)/(rd*atm(iqu)%temp))
  dg%tempc=atm(iqu)%temp*(dg%sigd/atm(iqu)%pa)**(rd/aircp)
  call getqsat(roadqsat,dg%tempc,dg%sigd)
  call getqsat(qsata,atm(iqu)%temp,atm(iqu)%pa)
  dg%mixrc=atm(iqu)%mixr*roadqsat/qsata

  ! roof (displacement height at building height)
  dg%sigr=atm(iqu)%ps*exp(grav*fn(iqu)%bldheight*(refheight-1.)/(rd*atm(iqu)%temp))
  dg%tempr=atm(iqu)%temp*(dg%sigr/atm(iqu)%pa)**(rd/aircp)
  call getqsat(roofqsat,dg%tempr,dg%sigr)
  dg%mixrr=atm(iqu)%mixr*roofqsat/qsata

  ! prep predictor-corrector arrays
  roofdum=roof(iqu)
  walledum=walle(iqu)
  wallwdum=wallw(iqu)
  roaddum=road(iqu)

  do j=1,2 ! predictor-corrector loop -------------------------------
    ! limit state variables
    roofdum%temp(:)=min(max(roofdum%temp(:),200.),400.)
    walledum%temp(:)=min(max(walledum%temp(:),200.),400.)
    wallwdum%temp(:)=min(max(wallwdum%temp(:),200.),400.)
    roaddum%temp(:)=min(max(roaddum%temp(:),200.),400.)
    roofdum%water=min(max(roofdum%water,0.),maxroofwater)
    roaddum%water=min(max(roaddum%water,0.),maxroadwater)
    roofdum%snow=min(max(roofdum%snow,0.),maxrfsn)
    roaddum%snow=min(max(roaddum%snow,0.),maxrdsn)
    roofdum%den=min(max(roofdum%den,minsnowden),maxsnowden)
    roaddum%den=min(max(roaddum%den,minsnowden),maxsnowden)
    roofdum%alpha=min(max(roofdum%alpha,minsnowalpha),maxsnowalpha)
    roaddum%alpha=min(max(roaddum%alpha,minsnowalpha),maxsnowalpha)
    
    ! water and snow diagnostics
    dg%roofdelta=(roofdum%water/maxroofwater)**(2./3.)
    dg%roaddelta=(roaddum%water/maxroadwater)**(2./3.)
    dg%rfsndelta=roofdum%snow/(roofdum%snow+maxrfsn)
    dg%rdsndelta=roaddum%snow/(roaddum%snow+maxrdsn)

    ! Adjust canyon roughness to include snow
    n=roaddum%snow/(roaddum%snow+maxrdsn+0.408*grav*fn(iqu)%zo) ! snow cover for urban roughness calc (Douville, et al 1995)
    zom=(1.-n)*fn(iqu)%zo+n*zosnow                              ! blended urban and snow roughness lenght
    dg%rfdzmin=max(abs(zmin-fn(iqu)%bldheight*(1.-refheight)),zoroof+1.)      ! roof displacement height 
    pg(iqu)%cndzmin=max(abs(zmin-fn(iqu)%bldheight*(2./3.-refheight)),zom+1.) ! canyon displacement height  
    pg(iqu)%lzom=log(pg(iqu)%cndzmin/zom) ! log of urban roughness length

    ! calculate shortwave radiation (up to 2nd order reflections)
    call getswcoeff(sg,wallpsi,roadpsi,fn(iqu),roaddum%alpha,dg%rdsndelta)
    sg%roof=(1.-roofalpha)*sg%roof*atm(iqu)%sg
    sg%walle=(1.-wallalpha)*sg%walle*atm(iqu)%sg
    sg%wallw=(1.-wallalpha)*sg%wallw*atm(iqu)%sg
    sg%road=(1.-roadalpha)*sg%road*atm(iqu)%sg
    sg%rfsn=(1.-roofdum%alpha)*sg%rfsn*atm(iqu)%sg
    sg%rdsn=(1.-roaddum%alpha)*sg%rdsn*atm(iqu)%sg

    ! calculate canyon wind speed and aerodynamical conductance
    if (zmin.ge.fn(iqu)%bldheight*(1.-refheight)) then
      topu=(2./pi)*atm(iqu)%umag*log((fn(iqu)%bldheight/3.)/zom)/pg(iqu)%lzom                   ! wind speed at canyon top
    else
      topu=(2./pi)*atm(iqu)%umag*exp(0.5*fn(iqu)%hwratio*(1.-refheight-zmin/fn(iqu)%bldheight)) ! wind speed at canyon top
    end if
    select case(resmeth)
      case(0) ! Masson (2000)
        cu=topu*exp(-0.25*fn(iqu)%hwratio)
        acond(:)=cu ! conductances are updated in solvecanyon
      case(1) ! Harman et al (2004) - modified to include 2nd wall
        lr=3.*fn(iqu)%bldheight
        wdt=fn(iqu)%bldheight/fn(iqu)%hwratio
        ln=max(0.,wdt-0.5*lr)
        rn=max(0.,fn(iqu)%bldheight*(1.-2.*ln/lr))
        ln=min(0.5*lr,ln)
        cu=topu*exp(-0.9*sqrt(ln**2+(fn(iqu)%bldheight-rn)**2)/fn(iqu)%bldheight)
        a=0.15*max(1.,1.5*fn(iqu)%hwratio)
        if (rn.ge.fn(iqu)%bldheight) then ! recirculation only
          xe=exp(-a)
          we=cu*(1.-xe)/a
          xw=exp(-a*wdt/fn(iqu)%bldheight)
          wr=cu*fn(iqu)%bldheight*xe*(1.-xw)/(wdt*a)
          ww=cu*xe*xw*(1.-exp(-a))/a
        else if (rn.gt.0.) then ! recirculation starts on east wall
          n=max(zocanyon,rn)
          cuven=topu*(1.-1./log(fn(iqu)%bldheight/zocanyon)+n/((fn(iqu)%bldheight-n)*(log(zocanyon/n) &
                /log(n/fn(iqu)%bldheight)+1.)))
          cuven=max(cuven*(1.-rn/fn(iqu)%bldheight),cu*(1.-exp(-a*(1.-rn/fn(iqu)%bldheight)))/a)
          xe=exp(-a*rn/fn(iqu)%bldheight)
          we=(cuven+cu*(1.-xe)/a)
          xw=exp(-a*wdt/fn(iqu)%bldheight)
          wr=cu*fn(iqu)%bldheight*xe*(1.-xw)/(wdt*a)
          ww=cu*xe*xw*(1.-exp(-a))/a
        else ! recirculation starts on road
          zolog=log(0.1*fn(iqu)%bldheight/zocanyon)
          cuven=topu*zolog/(2.3+zolog)
          xe=exp(-a*(wdt-lr)/fn(iqu)%bldheight)
          cuven=max(cuven*(wdt-lr),cu*(1.-xe)*fn(iqu)%bldheight/a)
          xw=exp(-a*lr/fn(iqu)%bldheight)
          wr=(cuven+cu*(1.-xw)*fn(iqu)%bldheight/a)/wdt
          cuven=topu*(1.-1./log(fn(iqu)%bldheight/zocanyon)+zocanyon/(fn(iqu)%bldheight-zocanyon))
          xe=cu*xe*(1.-exp(-a))/a
          we=max(cuven,xe)
          ww=cu*xw*(1.-exp(-a))/a
        end if
        zolog=log(0.1*fn(iqu)%bldheight/zocanyon)
        xe=vkar**2/(zolog*(2.3+zolog))
        xw=abs(atm(iqu)%udir)/pi
        acond(1)=xe*wr                 ! road conductance
        acond(2)=xe*(xw*ww+(1.-xw)*we) ! east wall conductance
        acond(3)=xe*(xw*we+(1.-xw)*ww) ! west wall conductance
        zolog=log(0.1*fn(iqu)%bldheight/zosnow)
        xe=vkar**2/(zolog*(2.3+zolog))
        acond(4)=xe*wr                 ! road snow conductance
    end select

    ! solve for road snow temperature -------------------------------
    ! includes solution to canyon temperature and longwave radiation
    ctmax=max(dg%tempc,roaddum%temp(1),walledum%temp(1),wallwdum%temp(1))+5. ! max road snow temp
    ctmin=min(dg%tempc,roaddum%temp(1),walledum%temp(1),wallwdum%temp(1))-5. ! min road snow temp
    sndepth=roaddum%snow*waterden/roaddum%den
    snlambda=icelambda*(roaddum%den/waterden)**1.88
    netldratio=0.5*(sndepth/snlambda+roaddepth(1)/roadlambda(1))    
    rdsntemp=0.5*(ctmax+ctmin)    
    if (roaddum%snow.gt.0.) then ! road snow
      call solverdsn(evctx,canyontemp,rg,fg,fgtop,eg,gardsn,rdsnmelt,ctmax,roaddum,walledum,wallwdum,dg,sg,atm(iqu) &
                    ,netldratio,ddt,iqu,acond,wallpsi,roadpsi)
      call solverdsn(evct,canyontemp,rg,fg,fgtop,eg,gardsn,rdsnmelt,rdsntemp,roaddum,walledum,wallwdum,dg,sg,atm(iqu) &
                    ,netldratio,ddt,iqu,acond,wallpsi,roadpsi)
      if ((evct*evctx).lt.0.) then
        ctmin=rdsntemp
      else
        ctmax=rdsntemp
      end if
      oldtemp=rdsntemp
      rdsntemp=0.5*(ctmax+ctmin)
      do k=1,5 ! sectant
        evctx=evct
        call solverdsn(evct,canyontemp,rg,fg,fgtop,eg,gardsn,rdsnmelt,rdsntemp,roaddum,walledum,wallwdum,dg,sg,atm(iqu) &
                      ,netldratio,ddt,iqu,acond,wallpsi,roadpsi)
        evctx=evct-evctx
        if (evctx.eq.0.) exit    
        newtemp=rdsntemp-evct*(rdsntemp-oldtemp)/evctx
        oldtemp=rdsntemp
        rdsntemp=newtemp
      end do
      rdsntemp=min(max(rdsntemp,ctmin),ctmax)
    else ! no road snow
      call solverdsn(evct,canyontemp,rg,fg,fgtop,eg,gardsn,rdsnmelt,rdsntemp,roaddum,walledum,wallwdum,dg,sg,atm(iqu) &
                    ,netldratio,ddt,iqu,acond,wallpsi,roadpsi)
      rg%rdsn=0.
      fg%rdsn=0.
      eg%rdsn=0.
      gardsn=0.
      rdsnmelt=0.   
    end if
    ! ---------------------------------------------------------------    

    ! solve for roof snow temperature -------------------------------
    ! includes solution to longwave radiation
    ctmax=max(dg%tempr,roofdum%temp(1))+5. ! max roof snow temp
    ctmin=min(dg%tempr,roofdum%temp(1))-5. ! min roof snow temp
    sndepth=roofdum%snow*waterden/roofdum%den
    snlambda=icelambda*(roofdum%den/waterden)**1.88
    netldratio=0.5*(sndepth/snlambda+roofdepth(1)/rooflambda(1))
    rfsntemp=0.5*(ctmax+ctmin)    
    if (roofdum%snow.gt.0.) then ! roof snow
      call solverfsn(evctx,rg,fg,eg,garfsn,rfsnmelt,ctmax,roofdum,dg,sg,atm(iqu),netldratio,ddt)
      call solverfsn(evct,rg,fg,eg,garfsn,rfsnmelt,rfsntemp,roofdum,dg,sg,atm(iqu),netldratio,ddt)
      if ((evct*evctx).lt.0.) then
        ctmin=rfsntemp
      else
        ctmax=rfsntemp
      end if
      oldtemp=rfsntemp
      rfsntemp=0.5*(ctmax+ctmin)
      do k=1,5 ! sectant
        evctx=evct
      call solverfsn(evct,rg,fg,eg,garfsn,rfsnmelt,rfsntemp,roofdum,dg,sg,atm(iqu),netldratio,ddt)
        evctx=evct-evctx
        if (evctx.eq.0.) exit    
        newtemp=rfsntemp-evct*(rfsntemp-oldtemp)/evctx
        oldtemp=rfsntemp
        rfsntemp=newtemp
      end do
      rfsntemp=min(max(rfsntemp,ctmin),ctmax)
    else ! no roof snow
      rg%rfsn=0.
      fg%rfsn=0.
      eg%rfsn=0.
      garfsn=0.
      rfsnmelt=0.
    end if
    !---------------------------------------------------------------- 

    ! calculate roof sensible and latent heat fluxes (without snow)
    rg%roof=roofemiss*(atm(iqu)%rg-sbconst*roofdum%temp(1)**4)
    lzoroof=log(dg%rfdzmin/zoroof)
    call getinvres(roofinvres,cd,lzotdum,lzoroof,dg%rfdzmin,roofdum%temp(1),dg%tempr,atm(iqu)%umag)
    fg%roof=aircp*atm(iqu)%rho*(roofdum%temp(1)-dg%tempr)*roofinvres 
    call getqsat(roofqsat,roofdum%temp(1),dg%sigr)
    if (roofqsat.lt.dg%mixrr) then
      eg%roof=lv*atm(iqu)%rho*(roofqsat-dg%mixrr)*roofinvres
    else
      eg%roof=lv*min(atm(iqu)%rho*dg%roofdelta*(roofqsat-dg%mixrr)*roofinvres,roofdum%water+atm(iqu)%rnd+rfsnmelt) ! MJT suggestion for max eg   
    end if

    ! calculate heat conduction (e.g., through walls)
    do k=1,2
      garoof(k)=2.*(roofdum%temp(k)-roofdum%temp(k+1))/(roofdepth(k)/rooflambda(k)+roofdepth(k+1)/rooflambda(k+1))
      gawalle(k)=2.*(walledum%temp(k)-walledum%temp(k+1))/(walldepth(k)/walllambda(k)+walldepth(k+1)/walllambda(k+1))
      gawallw(k)=2.*(wallwdum%temp(k)-wallwdum%temp(k+1))/(walldepth(k)/walllambda(k)+walldepth(k+1)/walllambda(k+1))
      garoad(k)=2.*(roaddum%temp(k)-roaddum%temp(k+1))/(roaddepth(k)/roadlambda(k)+roaddepth(k+1)/roadlambda(k+1))
    end do
    garoof(3)=2.*rooflambda(3)*(roofdum%temp(3)-bldtemp)/(roofdepth(3))
    gawalle(3)=2.*walllambda(3)*(walledum%temp(3)-bldtemp)/(walldepth(3))
    gawallw(3)=2.*walllambda(3)*(wallwdum%temp(3)-bldtemp)/(walldepth(3))
    garoad(3)=0.
  
    ! calculate change in urban temperatures
    croof%temp(1)=((1.-dg%rfsndelta)*(sg%roof+rg%roof-fg%roof-eg%roof-garoof(1)) &
                +dg%rfsndelta*(garfsn-garoof(1)))/(roofcp(1)*roofdepth(1))
    cwalle%temp(1)=(sg%walle+rg%walle-fg%walle-eg%walle-gawalle(1))/(wallcp(1)*walldepth(1))
    cwallw%temp(1)=(sg%wallw+rg%wallw-fg%wallw-eg%wallw-gawallw(1))/(wallcp(1)*walldepth(1))
    croad%temp(1)=((1.-dg%rdsndelta)*(sg%road+rg%road-fg%road-eg%road-garoad(1)) &
                +dg%rdsndelta*(gardsn-garoad(1)))/(roadcp(1)*roaddepth(1))
    do k=2,3
      croof%temp(k)=(garoof(k-1)-garoof(k))/(roofcp(k)*roofdepth(k))
      cwalle%temp(k)=(gawalle(k-1)-gawalle(k))/(wallcp(k)*walldepth(k))
      cwallw%temp(k)=(gawallw(k-1)-gawallw(k))/(wallcp(k)*walldepth(k))
      croad%temp(k)=(garoad(k-1)-garoad(k))/(roadcp(k)*roaddepth(k))
    end do
    
    ! calculate change in urban water
    croof%water=atm(iqu)%rnd-eg%roof/lv+rfsnmelt
    croad%water=atm(iqu)%rnd-eg%road/lv+rdsnmelt
    
    ! calculate change in urban snow
    croof%snow=atm(iqu)%snd-eg%rfsn/ls-rfsnmelt
    croad%snow=atm(iqu)%snd-eg%rdsn/ls-rdsnmelt
    croof%den=(0.24/86400.)*(maxsnowden-roofdum%den)
    croad%den=(0.24/86400.)*(maxsnowden-roaddum%den)
    if (rfsntemp.lt.273.16) then
      croof%alpha=-0.008/86400.
    else
      croof%alpha=(0.24/86400.)*(minsnowalpha-roofdum%alpha)
    end if
    if (rdsntemp.lt.273.16) then
      croad%alpha=-0.008/86400.
    else
      croad%alpha=(0.24/86400.)*(minsnowalpha-roaddum%alpha)
    end if

    ! predictor-corrector scheme to update urban state arrays (e.g., temperature, water and snow)
    if (j.eq.1) then ! predictor
      if (firstcall) then
        roofadj(iqu)=roofdum
        walleadj(iqu)=walledum
        wallwadj(iqu)=wallwdum
        roadadj(iqu)=roaddum
        firstcall=.false.
      end if
      roofdum%temp(:)=roof(iqu)%temp(:)+ddt*(1.5*croof%temp(:)-0.5*roofadj(iqu)%temp(:))
      walledum%temp(:)=walle(iqu)%temp(:)+ddt*(1.5*cwalle%temp(:)-0.5*walleadj(iqu)%temp(:))
      wallwdum%temp(:)=wallw(iqu)%temp(:)+ddt*(1.5*cwallw%temp(:)-0.5*wallwadj(iqu)%temp(:))
      roaddum%temp(:)=road(iqu)%temp(:)+ddt*(1.5*croad%temp(:)-0.5*roadadj(iqu)%temp(:))
      roofdum%water=roof(iqu)%water+ddt*(1.5*croof%water-0.5*roofadj(iqu)%water)
      roaddum%water=road(iqu)%water+ddt*(1.5*croad%water-0.5*roadadj(iqu)%water)
      roofdum%snow=roof(iqu)%snow+ddt*(1.5*croof%snow-0.5*roofadj(iqu)%snow)
      roaddum%snow=road(iqu)%snow+ddt*(1.5*croad%snow-0.5*roadadj(iqu)%snow)
      roofdum%den=roof(iqu)%den+ddt*(1.5*croof%den-0.5*roofadj(iqu)%den)
      roaddum%den=road(iqu)%den+ddt*(1.5*croad%den-0.5*roadadj(iqu)%den)
      roofdum%alpha=roof(iqu)%alpha+ddt*(1.5*croof%alpha-0.5*roofadj(iqu)%alpha)
      roaddum%alpha=road(iqu)%alpha+ddt*(1.5*croad%alpha-0.5*roadadj(iqu)%alpha)
      rooforg%temp(:)=croof%temp(:)
      walleorg%temp(:)=cwalle%temp(:)
      wallworg%temp(:)=cwallw%temp(:)
      roadorg%temp(:)=croad%temp(:)
      rooforg%water=croof%water
      roadorg%water=croad%water
      rooforg%snow=croof%snow
      roadorg%snow=croad%snow
      rooforg%den=croof%den
      roadorg%den=croad%den
      rooforg%alpha=croof%alpha
      roadorg%alpha=croad%alpha
    else ! corrector
      roofdum%temp(:)=roof(iqu)%temp(:)+(ddt/12.)*(5.*croof%temp(:)+8.*rooforg%temp(:)-roofadj(iqu)%temp(:))
      walledum%temp(:)=walle(iqu)%temp(:)+(ddt/12.)*(5.*cwalle%temp(:)+8.*walleorg%temp(:)-walleadj(iqu)%temp(:))
      wallwdum%temp(:)=wallw(iqu)%temp(:)+(ddt/12.)*(5.*cwallw%temp(:)+8.*wallworg%temp(:)-wallwadj(iqu)%temp(:))
      roaddum%temp(:)=road(iqu)%temp(:)+(ddt/12.)*(5.*croad%temp(:)+8.*roadorg%temp(:)-roadadj(iqu)%temp(:))
      roofdum%water=roof(iqu)%water+(ddt/12.)*(5.*croof%water+8.*rooforg%water-roofadj(iqu)%water)
      roaddum%water=road(iqu)%water+(ddt/12.)*(5.*croad%water+8.*roadorg%water-roadadj(iqu)%water)
      roofdum%snow=roof(iqu)%snow+(ddt/12.)*(5.*croof%snow+8.*rooforg%snow-roofadj(iqu)%snow)
      roaddum%snow=road(iqu)%snow+(ddt/12.)*(5.*croad%snow+8.*roadorg%snow-roadadj(iqu)%snow)
      roofdum%den=roof(iqu)%den+(ddt/12.)*(5.*croof%den+8.*rooforg%den-roofadj(iqu)%den)
      roaddum%den=road(iqu)%den+(ddt/12.)*(5.*croad%den+8.*roadorg%den-roadadj(iqu)%den)
      roofdum%alpha=roof(iqu)%alpha+(ddt/12.)*(5.*croof%alpha+8.*rooforg%alpha-roofadj(iqu)%alpha)
      roaddum%alpha=road(iqu)%alpha+(ddt/12.)*(5.*croad%alpha+8.*roadorg%alpha-roadadj(iqu)%alpha)
      roofadj(iqu)%temp(:)=rooforg%temp(:)
      walleadj(iqu)%temp(:)=walleorg%temp(:)
      wallwadj(iqu)%temp(:)=wallworg%temp(:)
      roadadj(iqu)%temp(:)=roadorg%temp(:)
      roofadj(iqu)%water=rooforg%water
      roadadj(iqu)%water=roadorg%water
      roofadj(iqu)%snow=rooforg%snow
      roadadj(iqu)%snow=roadorg%snow
      roofadj(iqu)%den=rooforg%den
      roadadj(iqu)%den=roadorg%den
      roofadj(iqu)%alpha=rooforg%alpha
      roadadj(iqu)%alpha=roadorg%alpha
    end if

  end do

  ! limit temperatures to sensible values
  roof(iqu)%temp(:)=min(max(roofdum%temp(:),200.),400.)
  walle(iqu)%temp(:)=min(max(walledum%temp(:),200.),400.)
  wallw(iqu)%temp(:)=min(max(wallwdum%temp(:),200.),400.)
  road(iqu)%temp(:)=min(max(roaddum%temp(:),200.),400.)
  roof(iqu)%water=min(max(roofdum%water,0.),maxroofwater)
  road(iqu)%water=min(max(roaddum%water,0.),maxroadwater)
  roof(iqu)%snow=min(max(roofdum%snow,0.),maxrfsn)
  road(iqu)%snow=min(max(roaddum%snow,0.),maxrdsn)
  roof(iqu)%den=min(max(roofdum%den,minsnowden),maxsnowden)
  road(iqu)%den=min(max(roaddum%den,minsnowden),maxsnowden)
  roof(iqu)%alpha=min(max(roofdum%alpha,minsnowalpha),maxsnowalpha)
  road(iqu)%alpha=min(max(roaddum%alpha,minsnowalpha),maxsnowalpha)

  ! calculate outputs
  fg%roof=dg%rfsndelta*fg%rfsn+(1.-dg%rfsndelta)*fg%roof ! redefine as net fg
  eg%roof=dg%rfsndelta*eg%rfsn+(1.-dg%rfsndelta)*eg%roof ! redefine as net eg
  egtop=dg%rdsndelta*eg%rdsn+(1.-dg%rdsndelta)*eg%road
  nettemp=dg%rfsndelta*rfsntemp+(1.-dg%rfsndelta)*roof(iqu)%temp(1)
  uo(iqu)%fg=fn(iqu)%sigmabld*fg%roof+(1.-fn(iqu)%sigmabld)*fgtop+fn(iqu)%industryfg
  uo(iqu)%eg=fn(iqu)%sigmabld*eg%roof+(1.-fn(iqu)%sigmabld)*egtop
  uo(iqu)%ts=fn(iqu)%sigmabld*nettemp+(1.-fn(iqu)%sigmabld)*canyontemp !MJT - since this is what the atmosphere can 'see'
  uo(iqu)%wf=fn(iqu)%sigmabld*dg%roofdelta*(1.-dg%rfsndelta)+(1.-fn(iqu)%sigmabld)*dg%roaddelta*(1.-dg%rdsndelta)

end do
  
return
end subroutine tebeval

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Estimate saturation mixing ratio

subroutine getqsat(qsat,temp,ps)

implicit none

real, intent(in) :: temp,ps
real, intent(out) :: qsat
real esatf,tdiff
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

subroutine getinvres(invres,cd,olzot,ilzo,zmin,stemp,theta,umag)

implicit none

real, intent(in) :: ilzo,zmin,stemp,theta,umag
real, intent(out) :: invres,cd,olzot
real af,aft,ri,fm,fh,root,denma,denha,re,lna,zodum
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

af=(vkar/ilzo)**2
ri=min(grav*zmin*(1.-stemp/theta)/max(umag,0.2)**2,rimax)

if (ri>0.) then
  fm=1./(1.+bprm*ri)**2
else
  root=sqrt(-ri*exp(ilzo))
  denma=1.+cms*2.*bprm*af*root
  fm=1.-2.*bprm *ri/denma
endif

cd=af*fm
zodum=zmin*exp(-ilzo)
re=sqrt(cd)*umag*zodum/nu
lna=2.46*re**0.25-log(7.4) !(Brutsaet, 1982)
olzot=lna+ilzo
aft=vkar**2/(ilzo*olzot)

if (ri>0.) then
  fh=fm
else
  denha=1.+chs*2.*bprm*aft*exp(0.5*lna)*root
  fh=1.-2.*bprm *ri/denha
endif

invres=aft*fh*umag

return
end subroutine getinvres

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate shortwave radiation coefficents (modified to include 2nd wall)

subroutine getswcoeff(sg,wallpsi,roadpsi,ifn,rdsnalpha,rdsndelta)

implicit none

real, intent(in) :: rdsnalpha,rdsndelta
real, intent(out) :: wallpsi,roadpsi
type(tdata), intent(in) :: ifn
type(trad), intent(out) :: sg
real thetazero,walles,wallws,roads,ta,tc,tz,xa,xb,ya,yb,roadnetalpha

wallpsi=0.5*(ifn%hwratio+1.-sqrt(ifn%hwratio**2+1.))/ifn%hwratio
roadpsi=sqrt(ifn%hwratio**2+1.)-ifn%hwratio

! integrate through 180 instead of 360
if (ifn%vangle.ge.(0.5*pi)) then
  walles=0.
  wallws=1./ifn%hwratio
  roads=0.
else
  ta=tan(ifn%vangle)
  thetazero=asin(1./max(ifn%hwratio*ta,1.))
  tc=2.*(1.-cos(thetazero))
  tz=2.*thetazero
  xa=max(ifn%hangle-thetazero,0.)-max(ifn%hangle-pi+thetazero,0.)-min(ifn%hangle+thetazero,0.)
  xb=pi-tz-xa
  ya=cos(min(thetazero,max(ifn%hangle-pi,0.)))-cos(min(thetazero,abs(ifn%hangle))) &
     +cos(pi-thetazero)-cos(min(pi,max(ifn%hangle,pi-thetazero)))
  yb=tc-ya
  ! note that these terms now include the azimuth angle
  walles=(xa/ifn%hwratio+ta*ya)/pi
  wallws=(xb/ifn%hwratio+ta*yb)/pi
  roads=(tz-ifn%hwratio*ta*tc)/pi
end if

! note that these terms are truncated to 2nd order reflections, compared to TEB which uses infinite reflections.
roadnetalpha=(1.-rdsndelta)*roadalpha+rdsndelta*rdsnalpha
sg%walle=walles+roadnetalpha*wallpsi*roads+wallalpha*(1.-2.*wallpsi)*wallws+(wallalpha*(1.-2.*wallpsi))**2*walles &
        +roadnetalpha*wallalpha*wallpsi*(1.-roadpsi)*wallws+roadnetalpha*wallalpha*wallpsi*(1.-2.*wallpsi)*roads
sg%wallw=wallws+roadnetalpha*wallpsi*roads+wallalpha*(1.-2.*wallpsi)*walles+(wallalpha*(1.-2.*wallpsi))**2*wallws &
        +roadnetalpha*wallalpha*wallpsi*(1.-roadpsi)*walles+roadnetalpha*wallalpha*wallpsi*(1.-2.*wallpsi)*roads
sg%road=roads+wallalpha*(1.-roadpsi)*0.5*(walles+wallws)+wallalpha*roadnetalpha*wallpsi*(1.-roadpsi)*roads &
        +wallalpha**2*(1.-roadpsi)*(1.-2.*wallpsi)*0.5*(walles+wallws)
sg%roof=1.
sg%rfsn=sg%roof
sg%rdsn=sg%road

return
end subroutine getswcoeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solve for roof snow temperature

subroutine solverfsn(evct,rg,fg,eg,garfsn,rfsnmelt,rfsntemp,roof,dg,sg,atm,ldratio,ddt)

implicit none

real, intent(out) :: evct,rfsnmelt,garfsn
real, intent(in) :: rfsntemp,ldratio,ddt
real rfsninvres,cd,rfsnqsat,lzotdum,lzosnow
type(trad), intent(inout) :: rg,fg,eg
type(trad), intent(in) :: sg
type(tatm), intent(in) :: atm
type(tdiag), intent(in) :: dg
type(tsurf), intent(in) :: roof

lzosnow=log(dg%rfdzmin/zosnow)
call getinvres(rfsninvres,cd,lzotdum,lzosnow,dg%rfdzmin,rfsntemp,dg%tempr,atm%umag)
call getqsat(rfsnqsat,rfsntemp,dg%sigr)
rfsnmelt=dg%rfsndelta*max(0.,rfsntemp-273.16)/(icecp*roof%den*lf*ddt) 
rg%rfsn=snowemiss*(atm%rg-sbconst*rfsntemp**4)
fg%rfsn=aircp*atm%rho*(rfsntemp-dg%tempr)*rfsninvres
eg%rfsn=ls*min(atm%rho*dg%rfsndelta*max(0.,rfsnqsat-dg%mixrr)*rfsninvres,roof%snow+atm%snd-rfsnmelt) ! MJT suggestion for max eg
garfsn=(rfsntemp-roof%temp(1))/ldratio
evct=sg%rfsn+rg%rfsn-fg%rfsn-eg%rfsn-garfsn

return
end subroutine solverfsn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solve for road snow temperature (includes canyon temperature)

subroutine solverdsn(evct,canyontemp,rg,fg,fgtop,eg,gardsn,rdsnmelt,rdsntemp,road,walle,wallw,dg,sg,atm &
                    ,ldratio,ddt,iqu,acond,wallpsi,roadpsi)

implicit none

integer, intent(in) :: iqu
integer k
real, dimension(4), intent(inout) :: acond
real, intent(out) :: evct,canyontemp,fgtop,gardsn,rdsnmelt
real, intent(in) :: rdsntemp,ldratio,ddt,wallpsi,roadpsi
real ctmax,ctmin,cevctx,cevct,oldtemp,newtemp,topinvres,roadqsat,rdsnqsat,canyonmix,nettemp,netemiss
type(trad), intent(inout) :: rg,fg,eg
type(trad), intent(in) :: sg
type(tatm), intent(in) :: atm
type(tdiag), intent(in) :: dg
type(tsurf), intent(in) :: road
type(twall), intent(in) :: walle,wallw

! solve for canyon temperature ----------------------------------
ctmax=max(dg%tempc,walle%temp(1),wallw%temp(1),road%temp(1),rdsntemp)+5. ! max canyon temp
ctmin=min(dg%tempc,walle%temp(1),wallw%temp(1),road%temp(1),rdsntemp)-5. ! min canyon temp
call solvecanyon(cevctx,fg,fgtop,topinvres,ctmax,dg,atm,walle%temp(1),wallw%temp(1),road%temp(1) &
                ,rdsntemp,acond,iqu)
canyontemp=0.5*(ctmax+ctmin)
call solvecanyon(cevct,fg,fgtop,topinvres,canyontemp,dg,atm,walle%temp(1),wallw%temp(1),road%temp(1) &
                ,rdsntemp,acond,iqu)
if ((cevct*cevctx).lt.0.) then
  ctmin=canyontemp
else
  ctmax=canyontemp
end if
oldtemp=canyontemp
canyontemp=0.5*(ctmax+ctmin)
do k=1,5 ! sectant
  cevctx=cevct
  call solvecanyon(cevct,fg,fgtop,topinvres,canyontemp,dg,atm,walle%temp(1),wallw%temp(1),road%temp(1) &
                  ,rdsntemp,acond,iqu)
  cevctx=cevct-cevctx
  if (cevctx.eq.0.) exit    
  newtemp=canyontemp-cevct*(canyontemp-oldtemp)/cevctx
  oldtemp=canyontemp
  canyontemp=newtemp
end do
canyontemp=min(max(canyontemp,ctmin),ctmax)
! ---------------------------------------------------------------    

nettemp=dg%rdsndelta*snowemiss*rdsntemp**4+(1.-dg%rdsndelta)*roademiss*road%temp(1)**4
netemiss=dg%rdsndelta*snowemiss+(1.-dg%rdsndelta)*roademiss
rg%walle=wallemiss*(atm%rg*(wallpsi+(1.-netemiss)*wallpsi*roadpsi+(1.-wallemiss)*wallpsi*(1.-2.*wallpsi)) &
                  +sbconst*walle%temp(1)**4*(-1.+wallemiss*(1.-wallemiss)*(1.-2.*wallpsi)**2) &
                  +sbconst*wallw%temp(1)**4*(wallemiss*(1.-2.*wallpsi)+wallemiss*(1.-netemiss)*wallpsi*(1.-roadpsi)) &
                  +sbconst*nettemp*(wallpsi+(1.-wallemiss)*wallpsi*(1.-2.*wallpsi)))
rg%wallw=wallemiss*(atm%rg*(wallpsi+(1.-netemiss)*wallpsi*roadpsi+(1.-wallemiss)*wallpsi*(1.-2.*wallpsi)) &
                  +sbconst*wallw%temp(1)**4*(-1.+wallemiss*(1.-wallemiss)*(1.-2.*wallpsi)**2) &
                  +sbconst*walle%temp(1)**4*(wallemiss*(1.-2.*wallpsi)+wallemiss*(1.-netemiss)*wallpsi*(1.-roadpsi)) &
                  +sbconst*nettemp*(wallpsi+(1.-wallemiss)*wallpsi*(1.-2.*wallpsi)))
rg%road=roademiss*(atm%rg*(roadpsi+(1.-wallemiss)*(1.-roadpsi)*wallpsi) &
                  +sbconst*(-road%temp(1)**4+nettemp*(1.-wallemiss)*(1.-roadpsi)*wallpsi) &
                  +sbconst*0.5*(walle%temp(1)**4+wallw%temp(1)**4) &
                  *(wallemiss*(1.-roadpsi)+wallemiss*(1.-wallemiss)*(1.-roadpsi)*(1.-2.*wallpsi)))
rg%rdsn=snowemiss*(atm%rg*(roadpsi+(1.-wallemiss)*(1.-roadpsi)*wallpsi) &
                  +sbconst*(-rdsntemp**4+nettemp*(1.-wallemiss)*(1.-roadpsi)*wallpsi) &
                  +sbconst*0.5*(walle%temp(1)**4+wallw%temp(1)**4) &
                  *(wallemiss*(1.-roadpsi)+wallemiss*(1.-wallemiss)*(1.-roadpsi)*(1.-2.*wallpsi)))

rdsnmelt=dg%rdsndelta*max(0.,rdsntemp-273.16)/(icecp*road%den*lf*ddt)
call getqsat(roadqsat,road%temp(1),dg%sigd) ! evaluate using pressure at displacement height
call getqsat(rdsnqsat,rdsntemp,dg%sigd)
canyonmix=((dg%rdsndelta*rdsnqsat*ls/lv*acond(4)+(1.-dg%rdsndelta)*dg%roaddelta*roadqsat*acond(1))+dg%mixrc*topinvres) &
          /((dg%rdsndelta*ls/lv*acond(4)+(1.-dg%rdsndelta)*dg%roaddelta*acond(1))+topinvres)
if (roadqsat.lt.canyonmix) then
  eg%road=lv*atm%rho*(roadqsat-canyonmix)*acond(1)
else
  eg%road=lv*min(atm%rho*dg%roaddelta*(roadqsat-canyonmix)*acond(1),road%water+atm%rnd+rdsnmelt) ! MJT suggestion for max eg
end if
eg%rdsn=ls*min(atm%rho*dg%rdsndelta*max(0.,rdsnqsat-canyonmix)*acond(4),road%snow+atm%snd-rdsnmelt) ! MJT suggestion for max eg
eg%walle=0.
eg%wallw=0.
gardsn=(rdsntemp-road%temp(1))/ldratio
evct=sg%rdsn+rg%rdsn-fg%rdsn-eg%rdsn-gardsn

return
end subroutine solverdsn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solve for canyon temperature

subroutine solvecanyon(evct,fg,fgtop,topinvres,ctemp,dg,atm,walletemp,wallwtemp,roadtemp,rdsntemp,acond,iqu)

implicit none

integer, intent(in) :: iqu
real, dimension(4), intent(inout) :: acond
real, intent(out) :: evct,fgtop,topinvres
real, intent(in) :: ctemp,walletemp,wallwtemp,roadtemp,rdsntemp
real cw,rwi
type(trad), intent(inout) :: fg
type(tdiag), intent(in) :: dg
type(tatm), intent(in) :: atm

call getinvres(topinvres,pg(iqu)%cduv,pg(iqu)%lzoh,pg(iqu)%lzom,pg(iqu)%cndzmin,ctemp,dg%tempc,atm%umag)
if (resmeth.eq.0) then
  cw=sqrt(pg(iqu)%cduv)*atm%umag ! diagnose canyonw (from Masson 2000)
  acond(:)=(11.8+4.2*sqrt(acond(:)**2+cw**2))/(aircp*atm%rho) ! From Rowley, et al (1930)
end if
fg%walle=aircp*atm%rho*(walletemp-ctemp)*acond(2)
fg%wallw=aircp*atm%rho*(wallwtemp-ctemp)*acond(3)
fg%road=aircp*atm%rho*(roadtemp-ctemp)*acond(1)
fg%rdsn=aircp*atm%rho*(rdsntemp-ctemp)*acond(4)
fgtop=aircp*atm%rho*(ctemp-dg%tempc)*topinvres
evct=fgtop-(dg%rdsndelta*fg%rdsn+(1.-dg%rdsndelta)*fg%road+fn(iqu)%hwratio*(fg%walle+fg%wallw) &
           +fn(iqu)%trafficfg/(1.-fn(iqu)%sigmabld))

return
end subroutine solvecanyon

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
