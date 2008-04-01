
! This code was originally based on the TEB scheme of Masson, Boundary-Layer Meteorology, 94, p357 (2000)
! The snow scheme is based on Douville, Royer and Mahfouf, Climate Dynamics, 12, p21 (1995)

! Usual pratice is:
!   call tebinit     ! to initalise state arrays, etc (use tebdisable to disable calls to ateb subroutines)
!   call tebloadm    ! to load previous state arrays (from tebsavem)
!   call tebtype     ! to define urban type (or use tebfndef instead)
!   ...
!   call tebnewangle ! define new solar zenith and azimuthal angle (use tebccangle for CCAM or
!                      use tebnewangle1 for a single grid point)
!   call tebalb      ! blends input and urban albedo (use tebalb1 for a single grid point)
!   call tebcalc     ! modifies input fluxes, surface temperature, etc to include urban
!   call tebzo       ! blends input and urban roughness lengths for momentum and heat
!   ...
!   call tebsavem    ! to save current state arrays (for use by tebloadm)
!   call tebend      ! to deallocate memory before quiting

! only tebinit and tebcalc are manditory.  All other subroutine calls are optional.

! NOTES: 
!  Below are differences with TEB (Masson 2000) scheme and aTEB:
!
! - aTEB assumes that the zero model height is at canyon displacement height (not roof level as in the TEB scheme).
!
! - aTEB uses two walls instead of the TEB single wall.  Also, only up to 2nd order reflections are used for
!   longwave and short wave radation.  In TEB, infinite reflections are used for shortwave, but only 2nd order
!   for long wave.
!
! - aTEB employs a much larger ratio between heat and momentum roughness lengths compared to canopies for vegetation.
!   (e.g., zot=zom/74.).  This is typical for urban schemes.
!
! - aTEB plans to improve traffic and industry heat fluxes as well as improve aerodynamical resistances within the
!   canyon.
!

module ateb

implicit none

private
public tebinit,tebcalc,tebend,tebzo,tebload,tebsave,tebtype,tebalb,tebalb1,tebfndef, &
       tebnewangle,tebnewangle1,tebccangle,tebdisable,tebloadm,tebsavem

! state arrays
integer ufull,maxtype
integer, dimension(:), allocatable :: ugrid,mgrid
real, dimension(:,:), allocatable :: rooftemp,walletemp,wallwtemp,roadtemp
real, dimension(:,:), allocatable :: roofadjt,walleadjt,wallwadjt,roadadjt
real, dimension(:), allocatable :: roofwater,roadwater
real, dimension(:), allocatable :: roofadjw,roadadjw
real, dimension(:), allocatable :: rfsn,rdsn,rfsnden,rdsnden,rfsnalpha,rdsnalpha
real, dimension(:), allocatable :: rfsnadjw,rdsnadjw,rfsnadjd,rdsnadjd,rfsnadja,rdsnadja
real, dimension(:), allocatable :: fnsigmabld,fnhwratio,fnindustryfg,fntrafficfg
real, dimension(:), allocatable :: fnzo,fnbldheight
real, dimension(:), allocatable :: vangle,hangle
! parameters
real, dimension(3), parameter :: roofdepth =(/ 0.05,0.4,0.1 /)          ! depth
real, dimension(3), parameter :: walldepth =(/ 0.02,0.125,0.05 /)
real, dimension(3), parameter :: roaddepth =(/ 0.05,0.1,1. /)
real, dimension(3), parameter :: roofcp =(/ 2.11E6,0.28E6,0.29E6 /)     ! heat capacity
real, dimension(3), parameter :: wallcp =(/ 1.55E6,1.55E6,0.29E6 /)
real, dimension(3), parameter :: roadcp =(/ 1.94E6,1.28E6,1.28E6 /)
real, dimension(3), parameter :: rooflambda =(/ 1.51,0.08,0.05 /)       ! conductance
real, dimension(3), parameter :: walllambda =(/ 0.9338,0.9338,0.05 /)
real, dimension(3), parameter :: roadlambda =(/ 0.7454,0.2513,0.2513 /)
real, parameter :: maxsnowden=300.   ! max snow density
real, parameter :: minsnowden=100.   ! min snow density
real, parameter :: waterden=1000.    ! water density
real, parameter :: icelambda=2.22    ! conductance of ice
real, parameter :: aircp=1004.64     ! Heat capapcity of dry air
real, parameter :: icecp=2100.       ! Heat capacity of ice
real, parameter :: bldtemp=291.16    ! Comfort temperature = 18deg C
real, parameter :: grav=9.80616      ! gravity
real, parameter :: lv=2.501e6        ! Latent heat of vaporisation
real, parameter :: lf=3.337e5        ! Latent heat of fusion
real, parameter :: ls=2.834e6        ! Latent heat of sublimation
real, parameter :: pi=3.1415927      ! pi
real, parameter :: rd=287.04         ! Gas constant for dry air
real, parameter :: sbconst=5.67e-8   ! Stefan-Boltzmann constant
real, parameter :: zoroof=0.15       ! Roughness length for rooftops (see Masson 2000)
real, parameter :: zosnow=0.001      ! Roughness length for snow
real, parameter :: roofemiss=0.90    ! emissitivity
real, parameter :: wallemiss=0.85 
real, parameter :: roademiss=0.94
real, parameter :: snowemiss=1.      ! snow emissitivity
real, parameter :: roofalpha=0.15    ! Albedo
real, parameter :: wallalpha=0.25
real, parameter :: roadalpha=0.08
real, parameter :: maxsnowalpha=0.85 ! max snow albedo
real, parameter :: minsnowalpha=0.5  ! min snow albedo
real, parameter :: maxroofwater=1.   ! max water on roof (mm)
real, parameter :: maxroadwater=1.   ! max water on road (mm)
real, parameter :: maxrfsn=1.        ! max snow on roof (mm)
real, parameter :: maxrdsn=1.        ! max snow on road (mm)
real, parameter :: heatzoinc = 4.3   ! Eva's veg (=2.) or MJT suggestion for urban (=4.3) which determines the
                                     ! ratio between momentum and heat roughness lengths (zot=zo/exp(heatzoinc))

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine prepainvres the arrays used by the TEB scheme
! This is a compulsory subroutine that must be called during
! model initalisation

subroutine tebinit(ifull,sigmau,diag)

implicit none

integer, intent(in) :: ifull,diag
integer iqu,iq
integer, dimension(ifull) :: utype
real, dimension(ifull), intent(in) :: sigmau

if (diag.ne.0) write(6,*) "Initialising aTEB"

ufull=count(sigmau.gt.0.)
if (ufull.eq.0) return

allocate(ugrid(ufull),mgrid(ifull))
allocate(rooftemp(ufull,3),walletemp(ufull,3),wallwtemp(ufull,3),roadtemp(ufull,3))
allocate(roofadjt(ufull,3),walleadjt(ufull,3),wallwadjt(ufull,3),roadadjt(ufull,3))
allocate(roofadjw(ufull),roadadjw(ufull))
allocate(roofwater(ufull),roadwater(ufull))
allocate(rfsn(ufull),rdsn(ufull),rfsnden(ufull),rdsnden(ufull),rfsnalpha(ufull),rdsnalpha(ufull))
allocate(rfsnadjw(ufull),rdsnadjw(ufull),rfsnadjd(ufull),rdsnadjd(ufull),rfsnadja(ufull),rdsnadja(ufull))
allocate(fnhwratio(ufull),fnsigmabld(ufull),fnindustryfg(ufull))
allocate(fntrafficfg(ufull),fnbldheight(ufull),fnzo(ufull))
allocate(vangle(ufull),hangle(ufull))

! define grid arrays
mgrid=0
iqu=0
do iq=1,ifull
  if (sigmau(iq).gt.0.) then
    iqu=iqu+1
    ugrid(iqu)=iq
    mgrid(iq)=iqu
  end if
end do

! Initialise state variables
rooftemp=bldtemp
walletemp=bldtemp
wallwtemp=bldtemp
roadtemp=bldtemp
roofwater=0.
roadwater=0.
vangle=0.
hangle=0.
rfsn=0.
rdsn=0.
rfsnden=minsnowden
rdsnden=minsnowden
rfsnalpha=maxsnowalpha
rdsnalpha=minsnowalpha

utype=1
call tebtype(ifull,utype,0)

return
end subroutine tebinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine deallocates arrays used by the TEB scheme

subroutine tebend(diag)

implicit none

integer, intent(in) :: diag

if (ufull.eq.0) return
if (diag.ne.0) write(6,*) "Deallocating aTEB arrays"

deallocate(ugrid,mgrid)
deallocate(rooftemp,walletemp,wallwtemp,roadtemp)
deallocate(roofadjt,walleadjt,wallwadjt,roadadjt)
deallocate(roofadjw,roadadjw)
deallocate(roofwater,roadwater)
deallocate(rfsn,rdsn,rfsnden,rdsnden,rfsnalpha,rdsnalpha)
deallocate(rfsnadjw,rdsnadjw,rfsnadjd,rdsnadjd,rfsnadja,rdsnadja)
deallocate(fnhwratio,fnsigmabld,fnindustryfg)
deallocate(fntrafficfg,fnbldheight,fnzo)
deallocate(vangle,hangle)

return
end subroutine tebend

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine loads aTEB state arrays (not compulsory)

subroutine tebload(ifull,urban,diag)

implicit none

integer, intent(in) :: ifull,diag
real, dimension(ifull,20), intent(in) :: urban

if (ufull.eq.0) return
if (diag.ne.0) write(6,*) "Load aTEB state arrays"

rooftemp(:,:)=urban(ugrid(:),1:3)
walletemp(:,:)=urban(ugrid(:),4:6)
wallwtemp(:,:)=urban(ugrid(:),7:9)
roadtemp(:,:)=urban(ugrid(:),10:12)
roofwater(:)=urban(ugrid(:),13)
roadwater(:)=urban(ugrid(:),14)
rfsn(:)=urban(ugrid(:),15)
rdsn(:)=urban(ugrid(:),16)
rfsnden(:)=urban(ugrid(:),17)
rdsnden(:)=urban(ugrid(:),18)
rfsnalpha(:)=urban(ugrid(:),19)
rdsnalpha(:)=urban(ugrid(:),20)

return
end subroutine tebload

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! temperature only version of tebsave

subroutine tebloadm(ifull,urban,diag)

implicit none

integer, intent(in) :: ifull,diag
real, dimension(ifull,12), intent(in) :: urban

if (ufull.eq.0) return
if (diag.ne.0) write(6,*) "Load aTEB state arrays"

rooftemp(:,:)=urban(ugrid(:),1:3)
walletemp(:,:)=urban(ugrid(:),4:6)
wallwtemp(:,:)=urban(ugrid(:),7:9)
roadtemp(:,:)=urban(ugrid(:),10:12)

return
end subroutine tebloadm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine loads aTEB type arrays (not compulsory)

subroutine tebtype(ifull,itype,diag)

implicit none

integer, intent(in) :: ifull,diag
integer, dimension(ifull), intent(in) :: itype
integer iqu

if (ufull.eq.0) return
if (diag.ne.0) write(6,*) "Load aTEB type arrays"

do iqu=1,ufull
  select case(itype(ugrid(iqu)))
    case DEFAULT
      ! default urban (ecosystems 007,152,153,154 with differences in sigmau only)
      fnhwratio(iqu)=0.21
      fnsigmabld(iqu)=0.5
      fnindustryfg(iqu)=5.
      fntrafficfg(iqu)=10.
      fnbldheight(iqu)=10.
      fnzo(iqu)=1.
    case(2)
      ! Dense urban (ecosystems 151)
      fnhwratio(iqu)=0.82
      fnsigmabld(iqu)=0.5
      fnindustryfg(iqu)=10.
      fntrafficfg(iqu)=20.
      fnbldheight(iqu)=30.
      fnzo(iqu)=3.
    case(3)
      ! Industries (ecosystems 155)
      fnhwratio(iqu)=0.4
      fnsigmabld(iqu)=0.5
      fnindustryfg(iqu)=20.
      fntrafficfg(iqu)=10.
      fnbldheight(iqu)=20.
      fnzo(iqu)=2.
    case(4)
      ! Road and rail (ecosystems 156)
      fnhwratio(iqu)=0.05
      fnsigmabld(iqu)=0.1
      fnindustryfg(iqu)=0.
      fntrafficfg(iqu)=30.
      fnbldheight(iqu)=5.
      fnzo(iqu)=0.5
    case(5)
      ! Ports (ecosystems 157)
      fnhwratio(iqu)=0.82
      fnsigmabld(iqu)=0.5
      fnindustryfg(iqu)=20.
      fntrafficfg(iqu)=10.
      fnbldheight(iqu)=20.
      fnzo(iqu)=2.
    case(6)
      ! Airport (ecosystems 158)
      fnhwratio(iqu)=0.02
      fnsigmabld(iqu)=0.1
      fnindustryfg(iqu)=0.
      fntrafficfg(iqu)=10.
      fnbldheight(iqu)=10.
      fnzo(iqu)=0.01
    case(7)
      ! Construction (ecosystems 159)
      fnhwratio(iqu)=0.005
      fnsigmabld(iqu)=0.1
      fnindustryfg(iqu)=0.
      fntrafficfg(iqu)=0.
      fnbldheight(iqu)=5.
      fnzo(iqu)=0.1
    case(8)
      ! Urban parks (ecosystems 160)
      fnhwratio(iqu)=0.005
      fnsigmabld(iqu)=0.1
      fnindustryfg(iqu)=0.
      fntrafficfg(iqu)=0.
      fnbldheight(iqu)=5.
      fnzo(iqu)=0.5
    case(9)
      ! Sport facilities (ecosystems 161)
      fnhwratio(iqu)=0.11
      fnsigmabld(iqu)=0.5
      fnindustryfg(iqu)=0.
      fntrafficfg(iqu)=0.
      fnbldheight(iqu)=10.
      fnzo(iqu)=1.
  end select
end do

return
end subroutine tebtype

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine defines the various urban parameters at each
! grid point, instead of by type
!

subroutine tebfndef(ifull,hwratioi,sigmabldi,industryfgi,trafficfgi,bldheighti,zoi,diag)

implicit none

integer, intent(in) :: ifull,diag
real, dimension(ifull), intent(in) :: hwratioi,sigmabldi,industryfgi,trafficfgi,bldheighti,zoi

if (ufull.eq.0) return

fnhwratio(:)=hwratioi(ugrid(:))
fnsigmabld(:)=sigmabldi(ugrid(:))
fnindustryfg(:)=industryfgi(ugrid(:))
fntrafficfg(:)=trafficfgi(ugrid(:))
fnbldheight(:)=bldheighti(ugrid(:))
fnzo(:)=zoi(ugrid(:))

return
end subroutine tebfndef

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine saves aTEB state arrays (not compulsory)

subroutine tebsave(ifull,urban,diag)

implicit none

integer, intent(in) :: ifull,diag
real, dimension(ifull,20), intent(inout) :: urban

if (ufull.eq.0) return
if (diag.ne.0) write(6,*) "Save aTEB state arrays"

urban(ugrid(:),1:3)=rooftemp(:,:)
urban(ugrid(:),4:6)=walletemp(:,:)
urban(ugrid(:),7:9)=wallwtemp(:,:)
urban(ugrid(:),10:12)=roadtemp(:,:)
urban(ugrid(:),13)=roofwater(:)
urban(ugrid(:),14)=roadwater(:)
urban(ugrid(:),15)=rfsn(:)
urban(ugrid(:),16)=rdsn(:)
urban(ugrid(:),17)=rfsnden(:)
urban(ugrid(:),18)=rdsnden(:)
urban(ugrid(:),19)=rfsnalpha(:)
urban(ugrid(:),20)=rdsnalpha(:)

return
end subroutine tebsave

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! temperature only version of tebsave

subroutine tebsavem(ifull,urban,diag)

implicit none

integer, intent(in) :: ifull,diag
real, dimension(ifull,12), intent(inout) :: urban

if (ufull.eq.0) return
if (diag.ne.0) write(6,*) "Save aTEB state arrays"

urban(ugrid(:),1:3)=rooftemp(:,:)
urban(ugrid(:),4:6)=walletemp(:,:)
urban(ugrid(:),7:9)=wallwtemp(:,:)
urban(ugrid(:),10:12)=roadtemp(:,:)

return
end subroutine tebsavem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine blends urban momentum and heat roughness lengths
! (This version neglects the displacement height (e.g., for CCAM))
!

subroutine tebzo(ifull,zom,zoh,zmin,sigmau,diag)

implicit none

integer, intent(in) :: ifull,diag
real, intent(in) :: zmin
real, dimension(ifull), intent(inout) :: zom,zoh
real, dimension(ifull), intent(in) :: sigmau
real, dimension(ufull) :: workb,workc,workp

workp(:)=rdsn(:)/(rdsn(:)+maxrdsn+0.408*grav*fnzo(:))
workb(:)=(1.-workp(:))*fnzo(:)+workp(:)*zosnow
workb(:)=sqrt((1.-sigmau(ugrid(:)))/log(zmin/zom(ugrid(:)))**2+sigmau(ugrid(:))/log(zmin/workb(:))**2)
workc(:)=(1.-sigmau(ugrid(:)))/(log(zmin/zom(ugrid(:)))*log(zmin/zoh(ugrid(:)))) &
        +sigmau(ugrid(:))/(log(zmin/fnzo(:))*(heatzoinc+log(zmin/fnzo(:))))
workc(:)=workc(:)/workb(:)
zom(ugrid(:))=zmin*exp(-1./workb(:))
zoh(ugrid(:))=zmin*exp(-1./workc(:))

return
end subroutine tebzo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates the urban contrabution to albedo.
! (all grid points)

subroutine tebalb(ifull,alb,sigmau,diag)

implicit none

integer, intent(in) :: ifull,diag
integer iq,iqu
real, dimension(ifull), intent(inout) :: alb
real, dimension(ifull), intent(in) :: sigmau

if (ufull.eq.0) return
if (diag.ne.0) write(6,*) "Calculate urban albedo"

do iqu=1,ufull
  iq=ugrid(iqu)
  call tebalbcalc(alb(iq),iqu,sigmau(iq))
end do

return
end subroutine tebalb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates the urban contrabution to albedo.
! (one grid point only)

subroutine tebalb1(alb,iq,sigmau,diag)

implicit none

integer, intent(in) :: iq,diag
integer iqu
real, intent(inout) :: alb
real, intent(in) :: sigmau

if (ufull.eq.0) return

iqu=mgrid(iq)
if (iqu.ge.1) call tebalbcalc(alb,iqu,sigmau)

return
end subroutine tebalb1

subroutine tebalbcalc(alb,iqu,sigmau)

implicit none

integer, intent(in) :: iqu
real, intent(inout) :: alb
real, intent(in) :: sigmau
real albu,wallesg,wallwsg,roadsg,wallpsi,roadpsi
real albr,snowdelta

snowdelta=rdsn(iqu)/(rdsn(iqu)+maxrdsn)
call getswcoeff(wallesg,wallwsg,roadsg,wallpsi,roadpsi,fnhwratio(iqu),vangle(iqu),hangle(iqu) &
                ,snowdelta,rdsnalpha(iqu))
albr=(1.-snowdelta)*roadalpha+snowdelta*rdsnalpha(iqu)
albu=fnhwratio(iqu)*(wallesg+wallwsg)*wallalpha+roadsg*albr
snowdelta=rfsn(iqu)/(rfsn(iqu)+maxrfsn)
albr=(1.-snowdelta)*roofalpha+snowdelta*rfsnalpha(iqu)
albu=fnsigmabld(iqu)*albr+(1.-fnsigmabld(iqu))*albu
alb=(1.-sigmau)*alb+sigmau*albu

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
if (diag.ne.0) write(6,*) "Update solar zenith angle and azimuth angle"

hangle(:)=0.5*pi-azimuthin(ugrid(:))
vangle(:)=acos(cosin(ugrid(:)))  

return
end subroutine tebnewangle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine stores the zenith angle and the solar azimuth angle
! (single grid point)

! iq = grid point
! cosin = cosine of zenith angle
! azmiuthin = azimuthal angle (rad)

subroutine tebnewangle1(iq,cosin,azimuthin)

implicit none

integer, intent(in) :: iq
integer iqu
real, intent(in) :: cosin     ! cosine of zenith angle
real, intent(in) :: azimuthin ! azimuthal angle

if (ufull.eq.0) return

iqu=mgrid(iq)
if (iqu.ge.1) then
  hangle(iqu)=0.5*pi-azimuthin
  vangle(iqu)=acos(cosin)  
end if

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

    hangle(iqu)=0.5*pi-azimuth
    vangle(iqu)=acos(cosin(i))
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
! rnd = incomming rainfall rate (kg/(m^2 s))
! snd = incomming snowfall rate (kg/(m^2 s))
! rho = atmospheric density at first model level (kg/m^3)
! temp = atmospheric temperature at first model level (K)
! ps = surface pressure (Pa)
! pa = pressure at first model level (Pa)
! umag = horizontal wind speed at first model level (m/s)
! ofg = Input/Output sensible heat flux (W/m^2)
! oeg = Input/Output latient heat flux (W/m^2)
! ots = Input/Output surface temperature (K)
! owf = Input/Output wetness fraction/surface water (%)

subroutine tebcalc(ifull,ofg,oeg,ots,owf,ddt,zmin,sg,rg,rnd,snd,rho,temp,mixr,ps,pa,umag,sigmau,diag)

implicit none

integer, intent(in) :: ifull,diag
real, intent(in) :: ddt,zmin
real, dimension(ifull), intent(in) :: sg,rg,rnd,snd,rho,temp,mixr,ps,pa,umag,sigmau
real, dimension(ifull), intent(inout) :: ofg,oeg,ots,owf
real, dimension(ufull) :: usg,urg,urnd,usnd,urho,utemp,umixr,ups,upa,uumag
real, dimension(ufull) :: uofg,uoeg,uots,uowf

if (ufull.eq.0) return

usg(:)=sg(ugrid(:))
urg(:)=rg(ugrid(:))
urnd(:)=rnd(ugrid(:))
usnd(:)=snd(ugrid(:))
urho(:)=rho(ugrid(:))
utemp(:)=temp(ugrid(:))
umixr(:)=mixr(ugrid(:))
ups(:)=ps(ugrid(:))
upa(:)=pa(ugrid(:))
uumag(:)=umag(ugrid(:))

call tebeval(uofg,uoeg,uots,uowf,ddt,zmin,usg,urg,urnd,usnd,urho,utemp,umixr,ups,upa,uumag,diag)

ofg(ugrid(:))=(1.-sigmau(ugrid(:)))*ofg(ugrid(:))+sigmau(ugrid(:))*uofg(:)
oeg(ugrid(:))=(1.-sigmau(ugrid(:)))*oeg(ugrid(:))+sigmau(ugrid(:))*uoeg(:)
ots(ugrid(:))=(1.-sigmau(ugrid(:)))*ots(ugrid(:))+sigmau(ugrid(:))*uots(:)
owf(ugrid(:))=(1.-sigmau(ugrid(:)))*owf(ugrid(:))+sigmau(ugrid(:))*uowf(:)

end subroutine tebcalc

subroutine tebeval(ofg,oeg,ots,owf,ddt,zmin,sg,rg,rnd,snd,rho,temp,mixr,ps,pa,umag,diag)

implicit none

integer, intent(in) :: diag
integer iqu,j,k
real, intent(in) :: ddt,zmin
real, dimension(ufull), intent(in) :: sg,rg,rnd,snd,rho,temp,mixr,ps,pa,umag
real, dimension(ufull), intent(out) :: ofg,oeg,ots,owf
real, dimension(3) :: roofga,wallega,wallwga,roadga
real, dimension(3) :: roofdumtemp,walledumtemp,wallwdumtemp,roaddumtemp,roofdumt,walledumt,wallwdumt,roaddumt
real, dimension(3) :: rooforgt,walleorgt,wallworgt,roadorgt
real rooforgw,roadorgw,roofdumwat,roaddumwat,roofdumw,roaddumw
real roofsg,roofrg,rooffg,roofeg
real wallesg,wallwsg,wallerg,wallwrg,wallefg,wallwfg
real roadsg,roadrg,roadfg,roadeg,topfg,topeg
real wallpsi,roadpsi,roofdelta,roaddelta
real roofinvres
real oldtemp,newtemp,canyontemp,canyonu,topu
real cd,roofqsat,roadqsat,qsata
real ctmax,ctmin,evctx,evct
real sigc,tempc,mixrc,sigr,tempr,mixrr,rfdzmin,cndzmin,efftrafffg
real rfsnsg,rdsnsg,rfsnrg,rdsnrg,rfsnfg,rdsnfg,rfsneg,rdsneg
real rfsntemp,rdsntemp,rfsndum,rdsndum
real rfsndumden,rdsndumden,rfsndumalpha,rdsndumalpha
real rfsnorgw,rdsnorgw,rfsndumw,rdsndumw
real rfsnorgd,rdsnorgd,rfsndumd,rdsndumd,rfsnorga,rdsnorga,rfsnduma,rdsnduma
real rfsndepth,rdsndepth,rfsnlambda,rdsnlambda,rfsncp,rdsncp
real rfsndelta,rdsndelta
real rfsnga,rdsnga,rfsnmelt,rdsnmelt
real netemiss,nettemp,netldratio
logical firstcall
data firstcall/.true./
save firstcall

if (diag.ne.0) write(6,*) "Evaluating aTEB"

do iqu=1,ufull
  ! canyon (displacement height at zero model height)
  tempc=temp(iqu)*(ps(iqu)/pa(iqu))**(rd/aircp)
  call getqsat(roadqsat,tempc,ps(iqu))
  call getqsat(qsata,temp(iqu),pa(iqu))
  mixrc=mixr(iqu)*roadqsat/qsata

  ! roof (displacement height at 1/3 building height)
  sigr=exp(-grav*fnbldheight(iqu)/(3.*rd*temp(iqu))) 
  tempr=temp(iqu)*(ps(iqu)*sigr/pa(iqu))**(rd/aircp)
  call getqsat(roofqsat,tempr,ps(iqu)*sigr)
  mixrr=mixr(iqu)*roofqsat/qsata

  ! road (displacement height at -2/3 building height)
  sigc=exp(2.*grav*fnbldheight(iqu)/(3.*rd*temp(iqu))) ! >1 for road

  ! diagnose canyon wind speed (rotating through 2pi)
  if (zmin.ge.fnbldheight(iqu)/3.) then
    topu=umag(iqu)*log(fnbldheight(iqu)/(3.*fnzo(iqu)))/log(zmin/fnzo(iqu))
  else
    topu=umag(iqu)/exp(-0.5*fnhwratio(iqu)*(1./3.-zmin/fnbldheight(iqu)))
  end if
  canyonu=(2./pi)*topu*exp(-0.25*fnhwratio(iqu)) ! Masson (2000) diagnosed wind speed

  ! scale traffic sensible heat flux for canyon          
  efftrafffg=fntrafficfg(iqu)/(1.-fnsigmabld(iqu))

  ! new snowfall
  if (snd(iqu).gt.0.) then
    rfsnden(iqu)=(rfsn(iqu)*rfsnden(iqu)+snd(iqu)*ddt*minsnowden)/(rfsn(iqu)+ddt*snd(iqu))
    rdsnden(iqu)=(rdsn(iqu)*rdsnden(iqu)+snd(iqu)*ddt*minsnowden)/(rdsn(iqu)+ddt*snd(iqu))
    rfsnalpha(iqu)=maxsnowalpha
    rdsnalpha(iqu)=maxsnowalpha
  end if

  ! prep predictor-corrector arrays
  roofdumtemp(:)=rooftemp(iqu,:)
  walledumtemp(:)=walletemp(iqu,:)
  wallwdumtemp(:)=wallwtemp(iqu,:)
  roaddumtemp(:)=roadtemp(iqu,:)
  roofdumwat=roofwater(iqu)
  roaddumwat=roadwater(iqu)
  rfsndum=rfsn(iqu)
  rdsndum=rdsn(iqu)
  rfsndumden=rfsnden(iqu)
  rdsndumden=rdsnden(iqu)
  rfsndumalpha=rfsnalpha(iqu)
  rdsndumalpha=rdsnalpha(iqu)

  do j=1,2 ! corrector-predictor loop -------------------------------
    ! limit state variables
    roofdumtemp(:)=min(max(roofdumtemp(:),225.),375.)
    walledumtemp(:)=min(max(walledumtemp(:),225.),375.)
    wallwdumtemp(:)=min(max(wallwdumtemp(:),225.),375.)
    roaddumtemp(:)=min(max(roaddumtemp(:),225.),375.)
    roofdumwat=min(max(roofdumwat,0.),maxroofwater)
    roaddumwat=min(max(roaddumwat,0.),maxroadwater)
    rfsndum=min(max(rfsndum,0.),maxrfsn)
    rdsndum=min(max(rdsndum,0.),maxrdsn)
    rfsndumden=min(max(rfsndumden,minsnowden),maxsnowden)
    rdsndumden=min(max(rdsndumden,minsnowden),maxsnowden)
    rfsndumalpha=min(max(rfsndumalpha,minsnowalpha),maxsnowalpha)
    rdsndumalpha=min(max(rdsndumalpha,minsnowalpha),maxsnowalpha)
    
    ! water and snow diagnostics
    roofdelta=(roofdumwat/maxroofwater)**(2./3.)
    roaddelta=(roaddumwat/maxroadwater)**(2./3.)
    rfsndelta=rfsndum/(rfsndum+maxrfsn)
    rdsndelta=rdsndum/(rdsndum+maxrdsn)
    rfsndepth=rfsndum*waterden/rfsndumden
    rdsndepth=rdsndum*waterden/rdsndumden
    rfsnlambda=icelambda*(rfsndumden/waterden)**1.88
    rdsnlambda=icelambda*(rdsndumden/waterden)**1.88
    rfsncp=icecp*rfsndumden
    rdsncp=icecp*rdsndumden

    ! calculate shortwave radiation
    call getswcoeff(wallesg,wallwsg,roadsg,wallpsi,roadpsi,fnhwratio(iqu),vangle(iqu),hangle(iqu) &
                   ,rdsndelta,rdsndumalpha)
    roofsg=(1.-roofalpha)*sg(iqu)
    wallesg=(1.-wallalpha)*wallesg*sg(iqu)
    wallwsg=(1.-wallalpha)*wallwsg*sg(iqu)
    roadsg=(1.-roadalpha)*roadsg*sg(iqu)
    rfsnsg=(1.-rfsndumalpha)*sg(iqu)
    rdsnsg=(1.-rdsndumalpha)*roadsg*sg(iqu)
    
    ! calculate long wave radiation
    netemiss=rdsndelta*snowemiss+(1.-rdsndelta)*roademiss
    nettemp=rdsndelta*snowemiss*rdsntemp**4+(1.-rdsndelta)*roademiss*roaddumtemp(1)**4
    roofrg=roofemiss*(rg(iqu)-sbconst*roofdumtemp(1)**4)
    wallerg=wallemiss*(rg(iqu)*(wallpsi+(1.-netemiss)*wallpsi*roadpsi+(1.-wallemiss)*wallpsi*(1.-2.*wallpsi)) &
                      +sbconst*walledumtemp(1)**4*(-1.+wallemiss*(1.-wallemiss)*(1.-2.*wallpsi)**2) &
                      +sbconst*wallwdumtemp(1)**4*(wallemiss*(1.-2.*wallpsi)+wallemiss*(1.-netemiss)*wallpsi*(1.-roadpsi)) &
                      +sbconst*nettemp*(wallpsi+(1.-wallemiss)*wallpsi*(1.-2.*wallpsi)))
    wallwrg=wallemiss*(rg(iqu)*(wallpsi+(1.-netemiss)*wallpsi*roadpsi+(1.-wallemiss)*wallpsi*(1.-2.*wallpsi)) &
                      +sbconst*wallwdumtemp(1)**4*(-1.+wallemiss*(1.-wallemiss)*(1.-2.*wallpsi)**2) &
                      +sbconst*walledumtemp(1)**4*(wallemiss*(1.-2.*wallpsi)+wallemiss*(1.-netemiss)*wallpsi*(1.-roadpsi)) &
                      +sbconst*nettemp*(wallpsi+(1.-wallemiss)*wallpsi*(1.-2.*wallpsi)))
    roadrg=roademiss*(rg(iqu)*(roadpsi+(1.-wallemiss)*(1.-roadpsi)*wallpsi) &
                      +sbconst*(-roaddumtemp(1)**4+nettemp*(1.-wallemiss)*(1.-roadpsi)*wallpsi) &
                      +sbconst*0.5*(walledumtemp(1)**4+wallwdumtemp(1)**4) &
                      *(wallemiss*(1.-roadpsi)+wallemiss*(1.-wallemiss)*(1.-roadpsi)*(1.-2.*wallpsi)))
    rfsnrg=snowemiss*(rg(iqu)-sbconst*rfsntemp**4)
    rdsnrg=snowemiss*(rg(iqu)*(roadpsi+(1.-wallemiss)*(1.-roadpsi)*wallpsi) &
                      +sbconst*(-rdsntemp**4+nettemp*(1.-wallemiss)*(1.-roadpsi)*wallpsi) &
                      +sbconst*0.5*(walledumtemp(1)**4+wallwdumtemp(1)**4) &
                      *(wallemiss*(1.-roadpsi)+wallemiss*(1.-wallemiss)*(1.-roadpsi)*(1.-2.*wallpsi)))

    ! calculate distances between roof and canyon displacement heights
    rfdzmin=max(abs(zmin-fnbldheight(iqu)/3.),zoroof+1.) ! MJT suggestion
    cndzmin=max(abs(zmin),fnzo(iqu)+1.) ! MJT suggestion    

    ! solve for roof snow temperature -------------------------------
    netldratio=0.5*(rfsndepth/rfsnlambda+roofdepth(1)/rooflambda(1))
    ctmax=max(tempr,roofdumtemp(1))+5. ! max roof snow temp
    ctmin=min(tempr,roofdumtemp(1))-5. ! min roof snow temp
    rfsntemp=0.5*(ctmax+ctmin)    
    if (rfsndum.gt.0.) then
      call solverfsn(evctx,rfsnmelt,rfsnfg,rfsneg,rfsnga,ctmax,rfsndum,rfsnsg,rfsnrg,rfsndelta,rfsncp,netldratio &
                     ,ddt,rfdzmin,tempr,umag(iqu),rho(iqu),ps(iqu)*sigr,mixrr,snd(iqu),roofdumtemp(1))
      call solverfsn(evct,rfsnmelt,rfsnfg,rfsneg,rfsnga,rfsntemp,rfsndum,rfsnsg,rfsnrg,rfsndelta,rfsncp,netldratio &
                     ,ddt,rfdzmin,tempr,umag(iqu),rho(iqu),ps(iqu)*sigr,mixrr,snd(iqu),roofdumtemp(1))
      if ((evct*evctx).lt.0.) then
        ctmin=rfsntemp
      else
        ctmax=rfsntemp
      end if
      oldtemp=rfsntemp
      rfsntemp=0.5*(ctmax+ctmin)
      do k=1,5 ! sectant
        evctx=evct
        call solverfsn(evct,rfsnmelt,rfsnfg,rfsneg,rfsnga,rfsntemp,rfsndum,rfsnsg,rfsnrg,rfsndelta,rfsncp,netldratio &
                       ,ddt,rfdzmin,tempr,umag(iqu),rho(iqu),ps(iqu)*sigr,mixrr,snd(iqu),roofdumtemp(1))
        evctx=evct-evctx
        if (evctx.eq.0.) exit    
        newtemp=rfsntemp-evct*(rfsntemp-oldtemp)/evctx
        oldtemp=rfsntemp
        rfsntemp=newtemp
      end do
      rfsntemp=min(max(rfsntemp,ctmin),ctmax)
    else
      rfsnfg=0.
      rfsneg=0.
      rfsnga=0.
      rfsnmelt=0.
    end if
    !---------------------------------------------------------------- 

    ! solve for road snow temperature -------------------------------
    ! includes solution to canyon temperature
    netldratio=0.5*(rdsndepth/rdsnlambda+roaddepth(1)/roadlambda(1))
    ctmax=max(tempc,roaddumtemp(1))+5. ! max roof snow temp
    ctmin=min(tempc,roaddumtemp(1))-5. ! min roof snow temp
    rdsntemp=0.5*(ctmax+ctmin)    
    if (rdsndum.gt.0.) then
      call solverdsn(evctx,canyontemp,wallefg,wallwfg,roadfg,rdsnfg,topfg,roadeg,rdsneg,rdsnga,rdsnmelt,ctmax &
                    ,rdsndelta,rdsndum,rdsncp,rdsnsg,rdsnrg,walledumtemp(1),wallwdumtemp(1),roaddumtemp(1),cndzmin &
                    ,tempc,mixrc,umag(iqu),canyonu,rho(iqu),efftrafffg,ps(iqu),sigc,roaddelta,roaddumwat,rnd(iqu) &
                    ,snd(iqu),netldratio,ddt,iqu)
      call solverdsn(evct,canyontemp,wallefg,wallwfg,roadfg,rdsnfg,topfg,roadeg,rdsneg,rdsnga,rdsnmelt,rdsntemp &
                    ,rdsndelta,rdsndum,rdsncp,rdsnsg,rdsnrg,walledumtemp(1),wallwdumtemp(1),roaddumtemp(1),cndzmin &
                    ,tempc,mixrc,umag(iqu),canyonu,rho(iqu),efftrafffg,ps(iqu),sigc,roaddelta,roaddumwat,rnd(iqu) &
                    ,snd(iqu),netldratio,ddt,iqu)
      if ((evct*evctx).lt.0.) then
        ctmin=rdsntemp
      else
        ctmax=rdsntemp
      end if
      oldtemp=rdsntemp
      rdsntemp=0.5*(ctmax+ctmin)
      do k=1,5 ! sectant
        evctx=evct
        call solverdsn(evct,canyontemp,wallefg,wallwfg,roadfg,rdsnfg,topfg,roadeg,rdsneg,rdsnga,rdsnmelt,rdsntemp &
                      ,rdsndelta,rdsndum,rdsncp,rdsnsg,rdsnrg,walledumtemp(1),wallwdumtemp(1),roaddumtemp(1),cndzmin &
                      ,tempc,mixrc,umag(iqu),canyonu,rho(iqu),efftrafffg,ps(iqu),sigc,roaddelta,roaddumwat,rnd(iqu) &
                      ,snd(iqu),netldratio,ddt,iqu)
        evctx=evct-evctx
        if (evctx.eq.0.) exit    
        newtemp=rdsntemp-evct*(rdsntemp-oldtemp)/evctx
        oldtemp=rdsntemp
        rdsntemp=newtemp
      end do
      rdsntemp=min(max(rdsntemp,ctmin),ctmax)
    else
      call solverdsn(evct,canyontemp,wallefg,wallwfg,roadfg,rdsnfg,topfg,roadeg,rdsneg,rdsnga,rdsnmelt,rdsntemp &
                    ,rdsndelta,rdsndum,rdsncp,rdsnsg,rdsnrg,walledumtemp(1),wallwdumtemp(1),roaddumtemp(1),cndzmin &
                    ,tempc,mixrc,umag(iqu),canyonu,rho(iqu),efftrafffg,ps(iqu),sigc,roaddelta,roaddumwat,rnd(iqu) &
                    ,snd(iqu),netldratio,ddt,iqu)
      rdsnfg=0.
      rdsneg=0.
      rdsnga=0.
      rdsnmelt=0.   
    end if
    ! ---------------------------------------------------------------    

    ! calculate roof sensible and latent heat fluxes
    call getinvres(roofinvres,cd,zoroof,rfdzmin,roofdumtemp(1),tempr,umag(iqu))
    rooffg=aircp*rho(iqu)*(roofdumtemp(1)-tempr)*roofinvres 
    call getqsat(roofqsat,roofdumtemp(1),ps(iqu)*sigr)
    if (roofqsat.lt.mixrr) then
      roofeg=lv*rho(iqu)*(roofqsat-mixrr)*roofinvres
    else
      roofeg=lv*min(rho(iqu)*roofdelta*(roofqsat-mixrr)*roofinvres,roofdumwat+rnd(iqu)+rfsnmelt) ! MJT suggestion for max eg   
    end if

    ! calculate heat conduction (e.g., through walls)
    do k=1,2
      roofga(k)=2.*(roofdumtemp(k)-roofdumtemp(k+1))/(roofdepth(k)/rooflambda(k)+roofdepth(k+1)/rooflambda(k+1))
      wallega(k)=2.*(walledumtemp(k)-walledumtemp(k+1))/(walldepth(k)/walllambda(k)+walldepth(k+1)/walllambda(k+1))
      wallwga(k)=2.*(wallwdumtemp(k)-wallwdumtemp(k+1))/(walldepth(k)/walllambda(k)+walldepth(k+1)/walllambda(k+1))
      roadga(k)=2.*(roaddumtemp(k)-roaddumtemp(k+1))/(roaddepth(k)/roadlambda(k)+roaddepth(k+1)/roadlambda(k+1))
    end do
    roofga(3)=2.*rooflambda(3)*(roofdumtemp(3)-bldtemp)/(roofdepth(3))
    wallega(3)=2.*walllambda(3)*(walledumtemp(3)-bldtemp)/(walldepth(3))
    wallwga(3)=2.*walllambda(3)*(wallwdumtemp(3)-bldtemp)/(walldepth(3))
    roadga(3)=0.
  
    ! calculate change in urban temperatures
    roofdumt(1)=((1.-rfsndelta)*(roofsg+roofrg-rooffg-roofeg-roofga(1)) &
                +rfsndelta*(rfsnga-roofga(1)))/(roofcp(1)*roofdepth(1))
    walledumt(1)=(wallesg+wallerg-wallefg-wallega(1))/(wallcp(1)*walldepth(1))
    wallwdumt(1)=(wallwsg+wallwrg-wallwfg-wallwga(1))/(wallcp(1)*walldepth(1))
    roaddumt(1)=((1.-rdsndelta)*(roadsg+roadrg-roadfg-roadeg-roadga(1)) &
                +rdsndelta*(rdsnga-roadga(1)))/(roadcp(1)*roaddepth(1))
    do k=2,3
      roofdumt(k)=(roofga(k-1)-roofga(k))/(roofcp(k)*roofdepth(k))
      walledumt(k)=(wallega(k-1)-wallega(k))/(wallcp(k)*walldepth(k))
      wallwdumt(k)=(wallwga(k-1)-wallwga(k))/(wallcp(k)*walldepth(k))
      roaddumt(k)=(roadga(k-1)-roadga(k))/(roadcp(k)*roaddepth(k))
    end do
    
    ! calculate change in urban water
    roofdumw=rnd(iqu)-roofeg/lv+rfsnmelt
    roaddumw=rnd(iqu)-roadeg/lv+rfsnmelt
    
    ! calculate change in urban snow
    rfsndumw=snd(iqu)-rfsneg/ls-rfsnmelt
    rdsndumw=snd(iqu)-rdsneg/ls-rdsnmelt
    rfsndumd=(0.24/86400.)*(maxsnowden-rfsndumden)
    rdsndumd=(0.24/86400.)*(maxsnowden-rdsndumden)
    if (rfsntemp.lt.273.16) then
      rfsnduma=-0.008/86400.
    else
      rfsnduma=(0.24/86400.)*(minsnowalpha-rfsndumalpha)
    end if
    if (rdsntemp.lt.273.16) then
      rdsnduma=-0.008/86400.
    else
      rdsnduma=(0.24/86400.)*(minsnowalpha-rdsndumalpha)
    end if

    ! predictor-corrector scheme to update urban state arrays (e.g., temperature, water and snow)
    if (j.eq.1) then ! predictor
      if (firstcall) then
        roofadjt(iqu,:)=roofdumt(:)
        walleadjt(iqu,:)=walledumt(:)
        wallwadjt(iqu,:)=wallwdumt(:)
        roadadjt(iqu,:)=roaddumt(:)
        roofadjw(iqu)=roofdumw
        roadadjw(iqu)=roaddumw
        rfsnadjw(iqu)=rfsndumw
        rdsnadjw(iqu)=rdsndumw
        rfsnadjd(iqu)=rfsndumd
        rdsnadjd(iqu)=rdsndumd
        rfsnadja(iqu)=rfsnduma
        rdsnadja(iqu)=rdsnduma
        firstcall=.false.
      end if
      roofdumtemp(:)=rooftemp(iqu,:)+ddt*(1.5*roofdumt(:)-0.5*roofadjt(iqu,:))
      walledumtemp(:)=walletemp(iqu,:)+ddt*(1.5*walledumt(:)-0.5*walleadjt(iqu,:))
      wallwdumtemp(:)=wallwtemp(iqu,:)+ddt*(1.5*wallwdumt(:)-0.5*wallwadjt(iqu,:))
      roaddumtemp(:)=roadtemp(iqu,:)+ddt*(1.5*roaddumt(:)-0.5*roadadjt(iqu,:))
      rooforgt(:)=roofdumt(:)
      walleorgt(:)=walledumt(:)
      wallworgt(:)=wallwdumt(:)
      roadorgt(:)=roaddumt(:)
      roofdumwat=roofwater(iqu)+ddt*(1.5*roofdumw-0.5*roofadjw(iqu))
      roaddumwat=roadwater(iqu)+ddt*(1.5*roaddumw-0.5*roadadjw(iqu))
      rooforgw=roofdumw
      roadorgw=roaddumw
      rfsndum=rfsn(iqu)+ddt*(1.5*rfsndumw-0.5*rfsnadjw(iqu))
      rdsndum=rdsn(iqu)+ddt*(1.5*rdsndumw-0.5*rdsnadjw(iqu))
      rfsnorgw=rfsndumw
      rdsnorgw=rdsndumw
      rfsndumden=rfsnden(iqu)+ddt*(1.5*rfsndumd-0.5*rfsnadjd(iqu))
      rdsndumden=rdsnden(iqu)+ddt*(1.5*rdsndumd-0.5*rdsnadjd(iqu))
      rfsnorgd=rfsndumd
      rdsnorgd=rdsndumd
      rfsndumalpha=rfsnalpha(iqu)+ddt*(1.5*rfsnduma-0.5*rfsnadja(iqu))
      rdsndumalpha=rdsnalpha(iqu)+ddt*(1.5*rdsnduma-0.5*rdsnadja(iqu))
      rfsnorga=rfsnduma
      rdsnorga=rdsnduma
    else ! corrector
      roofdumtemp(:)=rooftemp(iqu,:)+(ddt/12.)*(5.*roofdumt(:)+8.*rooforgt(:)-roofadjt(iqu,:))
      walledumtemp(:)=walletemp(iqu,:)+(ddt/12.)*(5.*walledumt(:)+8.*walleorgt(:)-walleadjt(iqu,:))
      wallwdumtemp(:)=wallwtemp(iqu,:)+(ddt/12.)*(5.*wallwdumt(:)+8.*wallworgt(:)-wallwadjt(iqu,:))
      roaddumtemp(:)=roadtemp(iqu,:)+(ddt/12.)*(5.*roaddumt(:)+8.*roadorgt(:)-roadadjt(iqu,:))
      roofadjt(iqu,:)=rooforgt(:)
      walleadjt(iqu,:)=walleorgt(:)
      wallwadjt(iqu,:)=wallworgt(:)
      roadadjt(iqu,:)=roadorgt(:)
      roofdumwat=roofwater(iqu)+(ddt/12.)*(5.*roofdumw+8.*rooforgw-roofadjw(iqu))
      roaddumwat=roadwater(iqu)+(ddt/12.)*(5.*roaddumw+8.*roadorgw-roadadjw(iqu))
      roofadjw(iqu)=rooforgw
      roadadjw(iqu)=roadorgw
      rfsndum=rfsn(iqu)+(ddt/12.)*(5.*rfsndumw+8.*rfsnorgw-rfsnadjw(iqu))
      rdsndum=rdsn(iqu)+(ddt/12.)*(5.*rdsndumw+8.*rdsnorgw-rdsnadjw(iqu))
      rfsnadjw(iqu)=rfsnorgw
      rdsnadjw(iqu)=rdsnorgw
      rfsndumden=rfsnden(iqu)+(ddt/12.)*(5.*rfsndumd+8.*rfsnorgd-rfsnadjd(iqu))
      rdsndumden=rdsnden(iqu)+(ddt/12.)*(5.*rdsndumd+8.*rdsnorgd-rdsnadjd(iqu))
      rfsnadjd(iqu)=rfsnorgd
      rdsnadjd(iqu)=rdsnorgd
      rfsndumalpha=rfsnalpha(iqu)+(ddt/12.)*(5.*rfsnduma+8.*rfsnorga-rfsnadja(iqu))
      rdsndumalpha=rdsnalpha(iqu)+(ddt/12.)*(5.*rdsnduma+8.*rdsnorga-rdsnadja(iqu))
      rfsnadja(iqu)=rfsnorga
      rdsnadja(iqu)=rdsnorga
    end if

  end do

  ! limit temperatures to sensible values
  rooftemp(iqu,:)=min(max(roofdumtemp(:),225.),375.)
  walletemp(iqu,:)=min(max(walledumtemp(:),225.),375.)
  wallwtemp(iqu,:)=min(max(wallwdumtemp(:),225.),375.)
  roadtemp(iqu,:)=min(max(roaddumtemp(:),225.),375.)
  roofwater(iqu)=min(max(roofdumwat,0.),maxroofwater)
  roadwater(iqu)=min(max(roaddumwat,0.),maxroadwater)
  rfsn(iqu)=min(max(rfsndum,0.),maxrfsn)
  rdsn(iqu)=min(max(rdsndum,0.),maxrdsn)
  rfsnden(iqu)=min(max(rfsndumden,minsnowden),maxsnowden)
  rdsnden(iqu)=min(max(rdsndumden,minsnowden),maxsnowden)
  rfsnalpha(iqu)=min(max(rfsndumalpha,minsnowalpha),maxsnowalpha)
  rdsnalpha(iqu)=min(max(rdsndumalpha,minsnowalpha),maxsnowalpha)

  ! calculate outputs
  rooffg=rfsndelta*rfsnfg+(1.-rfsndelta)*rooffg
  roofeg=rfsndelta*rfsneg+(1.-rfsndelta)*roofeg
  topeg=rdsndelta*rdsneg+(1.-rdsndelta)*roadeg
  nettemp=rfsndelta*rfsntemp+(1.-rfsndelta)*rooftemp(iqu,1)
  ofg(iqu)=fnsigmabld(iqu)*rooffg+(1.-fnsigmabld(iqu))*topfg+fnindustryfg(iqu)
  oeg(iqu)=fnsigmabld(iqu)*roofeg+(1.-fnsigmabld(iqu))*topeg
  ots(iqu)=fnsigmabld(iqu)*nettemp+(1.-fnsigmabld(iqu))*canyontemp !MJT - since this is what the atmosphere can 'see'
  owf(iqu)=fnsigmabld(iqu)*roofdelta*(1.-rfsndelta)+(1.-fnsigmabld(iqu))*roaddelta*(1.-rdsndelta)

end do
  
return
end subroutine tebeval

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
! Based on Eva's sflux.f version (i.e., vegetation)
! Modified for increased ratio between momentum and heat roughness
! lengths over urban areas.

subroutine getinvres(invres,cd,zo,zmin,stemp,theta,umag)

implicit none

real, intent(in) :: zo,zmin,stemp,theta,umag
real, intent(out) :: invres,cd
real zolog,af,aft,ri,fm,fh,root,denma,denha
real, parameter :: bprm=5. ! 4.7 in rams
real, parameter :: chs=2.6 ! 5.3 in rams
real, parameter :: cms=5.  ! 7.4 in rams
real, parameter :: vkar=0.4
real, parameter :: fmroot=0.57735
real, parameter :: rimax=(1./fmroot-1.)/bprm

zolog=log(zmin/zo)
af=(vkar/zolog)**2
aft=vkar**2/(zolog*(heatzoinc+zolog))
ri=min(grav*zmin*(1.-stemp/theta)/max(umag,0.2)**2,rimax)

if (ri>0.) then
  fm=1./(1.+bprm*ri)**2
  fh=fm
else
  root=sqrt(-ri*zmin/zo)
  denma=1.+cms*2.*bprm*af*root
  fm=1.-2.*bprm *ri/denma
  denha=1.+chs*2.*bprm*aft*sqrt(exp(heatzoinc))*root
  fh=1.-2.*bprm *ri/denha
endif

cd=af*fm
invres=aft*fh*umag

return
end subroutine getinvres

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate shortwave radiation coefficents (modified to include 2nd wall)

subroutine getswcoeff(wallesg,wallwsg,roadsg,wallpsi,roadpsi,hwratio,vangle,hangle &
                     ,rdsndelta,rdsnalpha)

implicit none

real, intent(in) :: hwratio,vangle,hangle,rdsndelta,rdsnalpha
real, intent(out) :: wallesg,wallwsg,roadsg,wallpsi,roadpsi
real thetazero,walles,wallws,roads,ta,tc,tz,xa,xb,ya,yb,roadnetalpha

wallpsi=0.5*(hwratio+1.-sqrt(hwratio**2+1.))/hwratio
roadpsi=sqrt(hwratio**2+1.)-hwratio

! integrate through 180 instead of 360
if (vangle.ge.(0.5*pi)) then
  walles=0.
  wallws=1./hwratio
  roads=0.
else
  ta=tan(vangle)
  thetazero=asin(1./max(hwratio*ta,1.))
  tc=2.*(1.-cos(thetazero))
  tz=2.*thetazero
  xa=max(hangle-thetazero,0.)-max(hangle-pi+thetazero,0.)+pi-min(hangle+pi+thetazero,pi)
  xb=pi-tz-xa
  ya=2.-cos(min(thetazero,max(abs(hangle),0.)))-cos(min(thetazero,max(hangle-pi+thetazero,0.)))
  yb=tc-ya
  ! note that these terms now include the azimuth angle
  walles=(xa/hwratio+ta*ya)/pi
  wallws=(xb/hwratio+ta*yb)/pi
  roads=(tz-hwratio*ta*tc)/pi
end if

! note that these terms are truncated to 2nd order reflections, compared to TEB which uses infinte reflections.
roadnetalpha=(1.-rdsndelta)*roadalpha+rdsndelta*rdsnalpha
wallesg=walles+roadnetalpha*wallpsi*roads+wallalpha*(1.-2.*wallpsi)*wallws+(wallalpha*(1.-2.*wallpsi))**2*walles &
        +roadnetalpha*wallalpha*wallpsi*(1.-roadpsi)*wallws+roadnetalpha*wallalpha*wallpsi*(1.-2.*wallpsi)*roads
wallwsg=wallws+roadnetalpha*wallpsi*roads+wallalpha*(1.-2.*wallpsi)*walles+(wallalpha*(1.-2.*wallpsi))**2*wallws &
        +roadnetalpha*wallalpha*wallpsi*(1.-roadpsi)*walles+roadnetalpha*wallalpha*wallpsi*(1.-2.*wallpsi)*roads
roadsg=roads+wallalpha*(1.-roadpsi)*0.5*(walles+wallws)+wallalpha*roadnetalpha*wallpsi*(1.-roadpsi)*roads &
        +wallalpha**2*(1.-roadpsi)*(1.-2.*wallpsi)*0.5*(walles+wallws)

return
end subroutine getswcoeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solve for canyon temperature

subroutine solvecanyon(evct,wallefg,wallwfg,roadfg,rdsnfg,topfg,rwinvres,topinvres,zo,zmin,ctemp &
  ,theta,umag,cu,walletemp,wallwtemp,roadtemp,rdsntemp,rho,efftrafficfg,hwratio,rdsndelta)

implicit none

real, intent(in) :: zo,zmin,ctemp,theta,umag,cu,walletemp,wallwtemp,roadtemp,rdsntemp,rho,efftrafficfg
real, intent(in) :: hwratio,rdsndelta
real, intent(out) :: evct,wallefg,wallwfg,roadfg,rdsnfg,topfg,rwinvres,topinvres
real cd,cw,rwi

call getinvres(topinvres,cd,zo,zmin,ctemp,theta,umag)
cw=sqrt(cd)*umag ! diagnose canyonw (from Masson 2000)
rwi=11.8+4.2*sqrt(cu**2+cw**2) ! From Rowley, et al (1930)
rwinvres=rwi/(aircp*rho)
wallefg=(walletemp-ctemp)*rwi
wallwfg=(wallwtemp-ctemp)*rwi
roadfg=(roadtemp-ctemp)*rwi
rdsnfg=(rdsntemp-ctemp)*rwi
topfg=aircp*rho*(ctemp-theta)*topinvres
evct=topfg-(rdsndelta*rdsnfg+(1.-rdsndelta)*roadfg+hwratio*(wallefg+wallwfg)+efftrafficfg)

return
end subroutine solvecanyon


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solve for roof snow temperature

subroutine solverfsn(evct,rfsnmelt,rfsnfg,rfsneg,rfsnga,rfsntemp,rfsn,rfsnsg,rfsnrg,rfsndelta,rfsncp,ldratio &
                     ,ddt,rfdzmin,tempr,umag,rho,psr,mixrr,snd,rooftemp)

implicit none

real, intent(out) :: evct,rfsnmelt,rfsnfg,rfsneg,rfsnga
real, intent(in) :: rfsntemp,rfsn,rfdzmin,tempr,umag,rfsndelta,rho,rfsncp,ddt,psr,mixrr,snd,ldratio,rooftemp
real, intent(in) :: rfsnsg,rfsnrg
real rfsninvres,cd,rfsnqsat

call getinvres(rfsninvres,cd,zosnow,rfdzmin,rfsntemp,tempr,umag)
call getqsat(rfsnqsat,rfsntemp,psr)
rfsnmelt=rfsndelta*max(0.,rfsntemp-273.16)/(rfsncp*lf*ddt) 
rfsnfg=aircp*rho*(rfsntemp-tempr)*rfsninvres
rfsneg=ls*min(rho*rfsndelta*max(0.,rfsnqsat-mixrr)*rfsninvres,rfsn+snd-rfsnmelt) ! MJT suggestion for max eg
rfsnga=(rfsntemp-rooftemp)/ldratio
evct=rfsnsg+rfsnrg-rfsnfg-rfsneg-rfsnga

return
end subroutine solverfsn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solve for road snow temperature (includes canyon temperature)

subroutine solverdsn(evct,canyontemp,wallefg,wallwfg,roadfg,rdsnfg,topfg,roadeg,rdsneg,rdsnga,rdsnmelt,rdsntemp &
                    ,rdsndelta,rdsn,rdsncp,rdsnsg,rdsnrg,walletemp,wallwtemp,roadtemp,cndzmin,tempc,mixrc,umag &
                    ,canyonu,rho,efftrafffg,ps,sigc,roaddelta,roadwat,rnd,snd,ldratio,ddt,iqu)

implicit none

integer, intent(in) :: iqu
real, intent(out) :: evct,canyontemp,wallefg,wallwfg,roadfg,rdsnfg,topfg,roadeg,rdsneg,rdsnga,rdsnmelt
real, intent(in) :: tempc,mixrc,walletemp,wallwtemp,roadtemp,cndzmin,umag,canyonu,rho,efftrafffg
real, intent(in) :: rdsntemp,rdsndelta,rdsn,rdsncp,rdsnsg,rdsnrg,ldratio,ps,roaddelta,roadwat,rnd,snd,ddt,sigc
real ctmax,ctmin,cevctx,cevct,oldtemp,newtemp,rwinvres,topinvres,roadqsat,rdsnqsat,canyonmix,tc,qc,qa,tcm
integer k

! solve for canyon temperature ----------------------------------
ctmax=max(tempc,walletemp,wallwtemp,roadtemp,rdsntemp)+5. ! max canyon temp
ctmin=min(tempc,walletemp,wallwtemp,roadtemp,rdsntemp)-5. ! min canyon temp
call solvecanyon(cevctx,wallefg,wallwfg,roadfg,rdsnfg,topfg,rwinvres,topinvres,fnzo(iqu) &
  ,cndzmin,ctmax,tempc,umag,canyonu,walletemp,wallwtemp,roadtemp &
  ,rdsntemp,rho,efftrafffg,fnhwratio(iqu),rdsndelta)
canyontemp=0.5*(ctmax+ctmin)
call solvecanyon(cevct,wallefg,wallwfg,roadfg,rdsnfg,topfg,rwinvres,topinvres,fnzo(iqu) &
  ,cndzmin,canyontemp,tempc,umag,canyonu,walletemp,wallwtemp,roadtemp &
  ,rdsntemp,rho,efftrafffg,fnhwratio(iqu),rdsndelta)
if ((cevct*cevctx).lt.0.) then
  ctmin=canyontemp
else
  ctmax=canyontemp
end if
oldtemp=canyontemp
canyontemp=0.5*(ctmax+ctmin)
do k=1,5 ! sectant
  cevctx=cevct
  call solvecanyon(cevct,wallefg,wallwfg,roadfg,rdsnfg,topfg,rwinvres,topinvres,fnzo(iqu) &
    ,cndzmin,canyontemp,tempc,umag,canyonu,walletemp,wallwtemp,roadtemp &
    ,rdsntemp,rho,efftrafffg,fnhwratio(iqu),rdsndelta)
  cevctx=cevct-cevctx
  if (cevctx.eq.0.) exit    
  newtemp=canyontemp-cevct*(canyontemp-oldtemp)/cevctx
  oldtemp=canyontemp
  canyontemp=newtemp
end do
canyontemp=min(max(canyontemp,ctmin),ctmax)
! ---------------------------------------------------------------    

rdsnmelt=rdsndelta*max(0.,rdsntemp-273.16)/(rdsncp*lf*ddt)
tc=canyontemp*sigc**(rd/aircp)
call getqsat(qc,tc,ps*sigc)
call getqsat(qa,canyontemp,ps)
call getqsat(roadqsat,roadtemp,ps*sigc)
call getqsat(rdsnqsat,rdsntemp,ps*sigc)
canyonmix=((rdsndelta*rdsnqsat*ls/lv+(1.-rdsndelta)*roaddelta*roadqsat)*rwinvres+mixrc*topinvres) &
          /((rdsndelta*ls/lv+(1.-rdsndelta)*roaddelta)*(qc/qa)*rwinvres+topinvres)
tcm=canyonmix*qc/qa
if (roadqsat.lt.tcm) then
  roadeg=lv*rho*(roadqsat-tcm)*rwinvres
else
  roadeg=lv*min(rho*roaddelta*(roadqsat-tcm)*rwinvres,roadwat+rnd+rdsnmelt) ! MJT suggestion for max eg
end if
rdsneg=ls*min(rho*rdsndelta*max(0.,rdsnqsat-tcm)*rwinvres,rdsn+snd-rdsnmelt) ! MJT suggestion for max eg
rdsnga=(rdsntemp-roadtemp)/ldratio
evct=rdsnsg+rdsnrg-rdsnfg-rdsneg-rdsnga

return
end subroutine solverdsn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine tebdisable(diag)

implicit none

integer, intent(in) :: diag

if (diag.ne.0) write(6,*) "Disable aTEB"
ufull=0

return
end subroutine tebdisable

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module ateb
