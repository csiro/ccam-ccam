
! This code was originally based on the TEB scheme of Masson, Boundary-Layer Meteorology, 94, p357 (2000)

! Usual pratice is:
!   call tebinit     ! to initalise state arrays, etc (use tebdisable to disable calls to ateb subroutines)
!   call tebload     ! to load previous state arrays (from tebsave)
!   call tebtype     ! to define urban type (or use tebfndef instead)
!   ...
!   call tebnewangle ! define new solar zenith and azimuthal angle (use tebccangle for CCAM or
!                      use tebnewangle1 for a single grid point)
!   call tebalb      ! modifies input albedo to include urban (use tebalb1 for a single grid point)
!   call tebcalc     ! modifies input fluxes, surface temperature, etc to include urban
!   call tebzo       ! blends input and urban momentum and heat roughness lengths
!   ...
!   call tebsave     ! to save current state arrays (for use by tebload)
!   call tebend      ! to deallocate memory before quiting

! only tebinit and tebcalc are manditory.  All other subroutine calls are optional.

! NOTES: 
!  Below are differences with TEB (Masson 2000) scheme and aTEB:
!
! - Currently snow is neglected for urban cover.
!
! - aTEB assumes that the surface is at street level (not roof level as in the TEB scheme).
!
! - aTEB uses two walls instead of the TEB single wall.  Also, only up to 2nd order reflections are used for
!   longwave and short wave radation.  In TEB, infinite reflections are used for shortwave, but only 2nd order
!   for long wave.
!

module ateb

implicit none

private
public tebinit,tebcalc,tebend,tebzo,tebzod,tebload,tebsave,tebtype,tebalb,tebalb1,tebfndef, &
       tebnewangle,tebnewangle1,tebccangle,tebdisable

! state arrays
integer ufull,maxtype
integer, dimension(:), allocatable :: ugrid,mgrid
real, dimension(:,:), allocatable :: rooftemp,walletemp,wallwtemp,roadtemp
real, dimension(:,:), allocatable :: roofadjt,walleadjt,wallwadjt,roadadjt
real, dimension(:), allocatable :: roofwater,roadwater
real, dimension(:), allocatable :: roofadjw,roadadjw
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
real, parameter :: aircp=1004.64   ! Specific heat of dry air
real, parameter :: bldtemp=291.16  ! Comfort temperature = 18deg C
real, parameter :: grav=9.80616    ! gravity
real, parameter :: lv=2.5104e6     ! Latent heat of vaporisation
real, parameter :: pi=3.1415927    ! pi
real, parameter :: rd=287.04       ! Gas constant for dry air
real, parameter :: sbconst=5.67e-8 ! Stefan-Boltzmann constant
real, parameter :: zoroof=0.15     ! Roughness length for rooftops (see Masson 2000)
real, parameter :: roofemiss=0.90  ! emissitivity
real, parameter :: wallemiss=0.85 
real, parameter :: roademiss=0.94
real, parameter :: roofalpha=0.15  ! Albedo
real, parameter :: wallalpha=0.25
real, parameter :: roadalpha=0.08
real, parameter :: maxroofwater=1. ! max water on roof (mm)
real, parameter :: maxroadwater=1. ! max water on road (mm)
real, parameter :: heatzoinc = 4.3 ! Eva's veg (=2.) or MJT suggestion for urban (=4.3) which determines the
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
real, dimension(ifull,14), intent(in) :: urban

if (ufull.eq.0) return
if (diag.ne.0) write(6,*) "Load aTEB state arrays"

rooftemp(:,:)=urban(ugrid(:),1:3)
walletemp(:,:)=urban(ugrid(:),4:6)
wallwtemp(:,:)=urban(ugrid(:),7:9)
roadtemp(:,:)=urban(ugrid(:),10:12)
roofwater(:)=urban(ugrid(:),13)
roadwater(:)=urban(ugrid(:),14)

return
end subroutine tebload

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
real, dimension(ifull,14), intent(inout) :: urban

if (ufull.eq.0) return
if (diag.ne.0) write(6,*) "Save aTEB state arrays"

urban(ugrid(:),1:3)=rooftemp(:,:)
urban(ugrid(:),4:6)=walletemp(:,:)
urban(ugrid(:),7:9)=wallwtemp(:,:)
urban(ugrid(:),10:12)=roadtemp(:,:)
urban(ugrid(:),13)=roofwater(:)
urban(ugrid(:),14)=roadwater(:)

return
end subroutine tebsave

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
real, dimension(ufull) :: workb,workc

workb(:)=sqrt((1.-sigmau(ugrid(:)))/log(zmin/zom(ugrid(:)))**2+sigmau(ugrid(:))/log(zmin/fnzo(:))**2)
workc(:)=(1.-sigmau(ugrid(:)))/(log(zmin/zom(ugrid(:)))*log(zmin/zoh(ugrid(:)))) &
        +sigmau(ugrid(:))/(log(zmin/fnzo(:))*(heatzoinc+log(zmin/fnzo(:))))
workc(:)=workc(:)/workb(:)
zom(ugrid(:))=zmin*exp(-1./workb(:))
zoh(ugrid(:))=zmin*exp(-1./workc(:))

return
end subroutine tebzo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine blends urban momentum and heat roughness lengths
! (This version includes the displacement height)
!

subroutine tebzod(ifull,zom,zoh,d,zmin,sigmau,rtg,diag)

implicit none

integer, intent(in) :: ifull,diag
real, intent(in) :: zmin
real, dimension(ifull), intent(inout) :: zom,zoh,d
real, dimension(ifull), intent(in) :: sigmau
real, dimension(ufull) :: worka,workb,workc,workd
logical, intent(in) :: rtg ! reduce displacement height to zero

where (zmin.ge.fnbldheight(:))
  worka(:)=log((zmin-2.*fnbldheight(:)/3.)/fnzo(:))
elsewhere ! MJT suggestion for inside canopy (using Masson (2000) diagnostic wind speed)
  worka(:)=exp(-0.5*fnhwratio(:)*(1.-zmin/fnbldheight(:)))*log(fnbldheight(:)/(3.*fnzo(:)))
end where

if (rtg) then ! reduce to ground level
  workd(:)=0. 
else ! MJT suggestion
  workd(:)=(1.-sigmau(ugrid(:)))*d(ugrid(:))+sigmau(ugrid(:))*fnbldheight(:)*2./3.
end if

workb(:)=sqrt((1.-sigmau(ugrid(:)))/log((zmin-d(ugrid(:)))/zom(ugrid(:)))**2 &
        +sigmau(ugrid(:))/worka(:)**2)
workc(:)=(1.-sigmau(ugrid(:)))/(log((zmin-d(ugrid(:)))/zom(ugrid(:)))*log((zmin-d(ugrid(:)))/zoh(ugrid(:)))) &
        +sigmau(ugrid(:))/(worka(:)*(heatzoinc+worka(:)))
workc(:)=workc(:)/workb(:)
zom(ugrid(:))=(zmin-workd(:))*exp(-1./workb(:))
zoh(ugrid(:))=(zmin-workd(:))*exp(-1./workc(:))
d(ugrid(:))=workd(:)

return
end subroutine tebzod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates the urban contrabution to albedo.
! (all grid points)

subroutine tebalb(ifull,alb,sigmau,diag)

implicit none

integer, intent(in) :: ifull,diag
integer iq,iqu
real, dimension(ifull), intent(inout) :: alb
real, dimension(ifull), intent(in) :: sigmau
real albu,wallesg,wallwsg,roadsg,wallpsi,roadpsi

if (ufull.eq.0) return
if (diag.ne.0) write(6,*) "Calculate urban albedo"

do iqu=1,ufull
  iq=ugrid(iqu)
  call getswcoeff(wallesg,wallwsg,roadsg,wallpsi,roadpsi,fnhwratio(iqu),vangle(iqu),hangle(iqu))
  albu=fnhwratio(iqu)*(wallesg+wallwsg)*(1.-wallalpha)+roadsg*(1.-roadalpha)
  albu=1.-(fnsigmabld(iqu)*(1.-roofalpha)+(1.-fnsigmabld(iqu))*albu)
  alb(iq)=(1.-sigmau(iqu))*alb(iq)+sigmau(iqu)*albu
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
real albu,wallesg,wallwsg,roadsg,wallpsi,roadpsi

if (ufull.eq.0) return

iqu=mgrid(iq)
if (iqu.ge.1) then
  call getswcoeff(wallesg,wallwsg,roadsg,wallpsi,roadpsi,fnhwratio(iqu),vangle(iqu),hangle(iqu))
  albu=fnhwratio(iqu)*(wallesg+wallwsg)*(1.-wallalpha)+roadsg*(1.-roadalpha)
  albu=1.-(fnsigmabld(iqu)*(1.-roofalpha)+(1.-fnsigmabld(iqu))*albu)
  alb=(1.-sigmau)*alb+sigmau*albu
end if

return
end subroutine tebalb1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine stores the zenith angle and the solar azimuth angle
! (all grid points)
!

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
!

subroutine tebnewangle1(iq,cosin,azimuthin,diag)

implicit none

integer, intent(in) :: iq,diag
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

subroutine tebccangle(is,ifull,cosin,rlon,rlat,fjd,slag,dhr,dlt,diag)

implicit none

integer, intent(in) :: is,ifull,diag
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
! dt = model time step
! zmin = first model level height (m)
! sg = incomming short wave radiation
! rg = incomming long wave radiation
! rnd = incomming rainfall rate (kg/(m^2 s))
! rho = atmospheric density at first model level
! temp = atmospheric temperature at first model level
! ps = surface pressure
! pa = pressure at first model level
! umag = horizontal wind speed at first model level
! ofg = Input/Output sensible heat flux
! oeg = Input/Output latient heat flux
! ots = Input/Output surface temperature
! owf = Input/Output wetness fraction (surface water)

subroutine tebcalc(ifull,ofg,oeg,ots,owf,ddt,zmin,sg,rg,rnd,rho,temp,mixr,ps,pa,umag,sigmau,diag)

implicit none

integer, intent(in) :: ifull,diag
real, intent(in) :: ddt,zmin
real, dimension(ifull), intent(in) :: sg,rg,rnd,rho,temp,mixr,ps,pa,umag,sigmau
real, dimension(ifull), intent(inout) :: ofg,oeg,ots,owf
real, dimension(ufull) :: usg,urg,urnd,urho,utemp,umixr,ups,upa,uumag
real, dimension(ufull) :: uofg,uoeg,uots,uowf

if (ufull.eq.0) return

usg(:)=sg(ugrid(:))
urg(:)=rg(ugrid(:))
urnd(:)=rnd(ugrid(:))
urho(:)=rho(ugrid(:))
utemp(:)=temp(ugrid(:))
umixr(:)=mixr(ugrid(:))
ups(:)=ps(ugrid(:))
upa(:)=pa(ugrid(:))
uumag(:)=umag(ugrid(:))

call tebeval(uofg,uoeg,uots,uowf,ddt,zmin,usg,urg,urnd,urho,utemp,umixr,ups,upa,uumag,diag)

ofg(ugrid(:))=(1.-sigmau(ugrid(:)))*ofg(ugrid(:))+sigmau(ugrid(:))*uofg(:)
oeg(ugrid(:))=(1.-sigmau(ugrid(:)))*oeg(ugrid(:))+sigmau(ugrid(:))*uoeg(:)
ots(ugrid(:))=(1.-sigmau(ugrid(:)))*ots(ugrid(:))+sigmau(ugrid(:))*uots(:)
owf(ugrid(:))=(1.-sigmau(ugrid(:)))*owf(ugrid(:))+sigmau(ugrid(:))*uowf(:)

end subroutine tebcalc

subroutine tebeval(ofg,oeg,ots,owf,ddt,zmin,sg,rg,rnd,rho,temp,mixr,ps,pa,umag,diag)

implicit none

integer, intent(in) :: diag
integer iqu,j,k
real, intent(in) :: ddt,zmin
real, dimension(ufull), intent(in) :: sg,rg,rnd,rho,temp,mixr,ps,pa,umag
real, dimension(ufull), intent(out) :: ofg,oeg,ots,owf
real, dimension(3) :: roofga,wallega,wallwga,roadga
real, dimension(3) :: roofdumtemp,walledumtemp,wallwdumtemp,roaddumtemp,roofdumt,walledumt,wallwdumt,roaddumt
real, dimension(3) :: rooforgt,walleorgt,wallworgt,roadorgt
real rooforgw,roadorgw,roofdumwat,roaddumwat,roofdumw,roaddumw
real roofsg,roofrg,rooffg,roofeg
real wallesg,wallwsg,wallerg,wallwrg,wallefg,wallwfg
real roadsg,roadrg,roadfg,roadeg,topfg,topeg
real wallpsi,roadpsi,roofdelta,roaddelta
real roofinvres,rwinvres,topinvres
real oldcanyontemp,newcanyontemp,canyontemp,canyonmix,canyonu,topu
real cd,roofqsat,roadqsat,qsata
real ctmax,ctmin,evctx,evct
real sigc,tempc,mixrc,sigr,tempr,mixrr,dzmin,efftrafffg
logical firstcall
data firstcall/.true./
save firstcall

if (diag.ne.0) write(6,*) "Evaluating aTEB"

do iqu=1,ufull
  ! calculate shortwave radiation
  call getswcoeff(wallesg,wallwsg,roadsg,wallpsi,roadpsi,fnhwratio(iqu),vangle(iqu),hangle(iqu))
  roofsg=(1.-roofalpha)*sg(iqu)
  wallesg=(1.-wallalpha)*wallesg*sg(iqu)
  wallwsg=(1.-wallalpha)*wallwsg*sg(iqu)
  roadsg=(1.-roadalpha)*roadsg*sg(iqu)

  ! canyon (displacement height at 2/3 building height)
  sigc=exp(-2.*grav*fnbldheight(iqu)/(3.*rd*temp(iqu)))
  tempc=temp(iqu)*(ps(iqu)*sigc/pa(iqu))**(rd/aircp)
  call getqsat(roadqsat,tempc,ps(iqu)*sigc)
  call getqsat(qsata,temp(iqu),pa(iqu))
  mixrc=mixr(iqu)*roadqsat/qsata

  ! roof (displacement height at building height)
  sigr=exp(-grav*fnbldheight(iqu)/(rd*temp(iqu))) 
  tempr=temp(iqu)*(ps(iqu)*sigr/pa(iqu))**(rd/aircp)
  call getqsat(roofqsat,tempr,ps(iqu)*sigr)
  mixrr=mixr(iqu)*roofqsat/qsata

  ! diagnose canyon wind speed (rotating through 2pi)
  if (zmin.ge.fnbldheight(iqu)) then
    topu=umag(iqu)*log(fnbldheight(iqu)/(3.*fnzo(iqu)))/log((zmin-fnbldheight(iqu)*2./3.)/fnzo(iqu))
  else
    topu=umag(iqu)/exp(-0.5*fnhwratio(iqu)*(1.-zmin/fnbldheight(iqu)))
  end if
  canyonu=(2./pi)*topu*exp(-0.25*fnhwratio(iqu)) ! Masson (2000) diagnosed wind speed
                                                 ! Raupach (1992) exponential profile may be better

  ! scale traffic sensible heat flux for canyon          
  efftrafffg=fntrafficfg(iqu)/(1.-fnsigmabld(iqu))

  ! prep predictor-corrector arrays
  roofdumtemp(:)=rooftemp(iqu,:)
  walledumtemp(:)=walletemp(iqu,:)
  wallwdumtemp(:)=wallwtemp(iqu,:)
  roaddumtemp(:)=roadtemp(iqu,:)
  roofdumwat=roofwater(iqu)
  roaddumwat=roadwater(iqu)

  do j=1,2 ! corrector-predictor loop
    roofdumtemp(:)=min(max(roofdumtemp(:),225.),375.)
    walledumtemp(:)=min(max(walledumtemp(:),225.),375.)
    wallwdumtemp(:)=min(max(wallwdumtemp(:),225.),375.)
    roaddumtemp(:)=min(max(roaddumtemp(:),225.),375.)
    roofdumwat=min(max(roofdumwat,0.),maxroofwater)
    roaddumwat=min(max(roaddumwat,0.),maxroadwater)
    
    ! calculate long wave radiation
    roofrg=roofemiss*(rg(iqu)-sbconst*roofdumtemp(1)**4)
    wallerg=wallemiss*(rg(iqu)*(wallpsi+(1.-roademiss)*wallpsi*roadpsi+(1.-wallemiss)*wallpsi*(1.-2.*wallpsi)) &
                      +sbconst*walledumtemp(1)**4*(-1.+wallemiss*(1.-wallemiss)*(1.-2.*wallpsi)**2) &
                      +sbconst*wallwdumtemp(1)**4*(wallemiss*(1.-2.*wallpsi)+wallemiss*(1.-roademiss)*wallpsi*(1.-roadpsi)) &
                      +sbconst*roaddumtemp(1)**4*(wallpsi*roademiss+(1.-wallemiss)*wallpsi*(1.-2.*wallpsi)))
    wallwrg=wallemiss*(rg(iqu)*(wallpsi+(1.-roademiss)*wallpsi*roadpsi+(1.-wallemiss)*wallpsi*(1.-2.*wallpsi)) &
                      +sbconst*wallwdumtemp(1)**4*(-1.+wallemiss*(1.-wallemiss)*(1.-2.*wallpsi)**2) &
                      +sbconst*walledumtemp(1)**4*(wallemiss*(1.-2.*wallpsi)+wallemiss*(1.-roademiss)*wallpsi*(1.-roadpsi)) &
                      +sbconst*roaddumtemp(1)**4*(wallpsi*roademiss+(1.-wallemiss)*wallpsi*(1.-2.*wallpsi)))
    roadrg=roademiss*(rg(iqu)*(roadpsi+(1.-wallemiss)*(1.-roadpsi)*wallpsi) &
                      +sbconst*roaddumtemp(1)**4*(-1.+(1.-wallemiss)*(1.-roadpsi)*wallpsi) &
                      +sbconst*0.5*(walledumtemp(1)**4+wallwdumtemp(1)**4) &
                      *(wallemiss*(1.-roadpsi)+wallemiss*(1.-wallemiss)*(1.-roadpsi)*(1.-2.*wallpsi)))

    ! calculate roof sensible heat flux
    dzmin=max(abs(zmin-fnbldheight(iqu)),zoroof+1.) ! MJT suggestion for symmetry
    call getinvres(roofinvres,cd,zoroof,dzmin,roofdumtemp(1),tempr,umag(iqu)) 
    rooffg=aircp*rho(iqu)*(roofdumtemp(1)-tempr)*roofinvres
    
    ! calculate canyon sensible heat fluxes 
    dzmin=max(abs(zmin-2.*fnbldheight(iqu)/3.),fnzo(iqu)+1.) ! MJT suggestion for symmetry
    ctmax=max(tempc,walledumtemp(1),wallwdumtemp(1),roaddumtemp(1))+1. ! max canyon temp
    ctmin=min(tempc,walledumtemp(1),wallwdumtemp(1),roaddumtemp(1))    ! min canyon temp
    call solvecanyon(evctx,wallefg,wallwfg,roadfg,topfg,rwinvres,topinvres,fnzo(iqu) &
      ,dzmin,ctmax,tempc,umag(iqu),canyonu,walledumtemp(1),wallwdumtemp(1),roaddumtemp(1) &
      ,rho(iqu),efftrafffg,fnhwratio(iqu))
    canyontemp=0.5*(ctmax+ctmin)
    do k=1,2 ! bisect
      call solvecanyon(evct,wallefg,wallwfg,roadfg,topfg,rwinvres,topinvres,fnzo(iqu) &
        ,dzmin,canyontemp,tempc,umag(iqu),canyonu,walledumtemp(1),wallwdumtemp(1),roaddumtemp(1) &
        ,rho(iqu),efftrafffg,fnhwratio(iqu))
      if ((evct*evctx).lt.0.) then
        ctmin=canyontemp
      else
        ctmax=canyontemp
        evctx=evct
      end if
      oldcanyontemp=canyontemp
      canyontemp=0.5*(ctmax+ctmin)
    end do
    do k=1,5 ! sectant
      evctx=evct
      call solvecanyon(evct,wallefg,wallwfg,roadfg,topfg,rwinvres,topinvres,fnzo(iqu) &
        ,dzmin,canyontemp,tempc,umag(iqu),canyonu,walledumtemp(1),wallwdumtemp(1),roaddumtemp(1) &
        ,rho(iqu),efftrafffg,fnhwratio(iqu))
      evctx=evct-evctx
      if (evctx.eq.0.) exit    
      newcanyontemp=canyontemp-evct*(canyontemp-oldcanyontemp)/evctx
      oldcanyontemp=canyontemp
      canyontemp=newcanyontemp
    end do
    canyontemp=min(max(canyontemp,ctmin),ctmax)    

    ! calculate latent heat fluxes
    call getqsat(roofqsat,roofdumtemp(1),ps(iqu)*sigr)
    call getqsat(roadqsat,roaddumtemp(1),ps(iqu))
    roofdelta=(roofdumwat/maxroofwater)**(2./3.)
    roaddelta=(roaddumwat/maxroadwater)**(2./3.)
    canyonmix=(roaddelta*roadqsat*rwinvres+mixrc*topinvres)/(roaddelta*rwinvres+topinvres)
    if (roofqsat.lt.mixrr) then
      roofeg=lv*rho(iqu)*(roofqsat-mixrr)*roofinvres ! MJT suggestion for dew
    else
      roofeg=lv*min(rho(iqu)*roofdelta*(roofqsat-mixrr)*roofinvres,roofdumwat+rnd(iqu)) ! MJT suggestion for max eg   
    end if
    if (roadqsat.lt.mixrc) then ! equivilent to roadqsat.lt.canyonmix    
      roadeg=lv*rho(iqu)*(roadqsat-canyonmix)*rwinvres ! MJT suggestion for dew
    else
      roadeg=lv*min(rho(iqu)*roaddelta*(roadqsat-canyonmix)*rwinvres,roaddumwat+rnd(iqu)) ! MJT suggestion for max eg
    end if
    topeg=roadeg
  
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
    roofdumt(1)=(roofsg+roofrg-rooffg-roofeg-roofga(1))/(roofcp(1)*roofdepth(1))
    walledumt(1)=(wallesg+wallerg-wallefg-wallega(1))/(wallcp(1)*walldepth(1))
    wallwdumt(1)=(wallwsg+wallwrg-wallwfg-wallwga(1))/(wallcp(1)*walldepth(1))
    roaddumt(1)=(roadsg+roadrg-roadfg-roadeg-roadga(1))/(roadcp(1)*roaddepth(1))
    do k=2,3
      roofdumt(k)=(roofga(k-1)-roofga(k))/(roofcp(k)*roofdepth(k))
      walledumt(k)=(wallega(k-1)-wallega(k))/(wallcp(k)*walldepth(k))
      wallwdumt(k)=(wallwga(k-1)-wallwga(k))/(wallcp(k)*walldepth(k))
      roaddumt(k)=(roadga(k-1)-roadga(k))/(roadcp(k)*roaddepth(k))
    end do
    roofdumw=(rnd(iqu)-roofeg/lv)
    roaddumw=(rnd(iqu)-roadeg/lv)

    ! predictor-corrector scheme to update urban temperatures
    if (j.eq.1) then ! predictor
      if (firstcall) then
        roofadjt(iqu,:)=roofdumt(:)
        walleadjt(iqu,:)=walledumt(:)
        wallwadjt(iqu,:)=wallwdumt(:)
        roadadjt(iqu,:)=roaddumt(:)
        roofadjw(iqu)=roofdumw
        roadadjw(iqu)=roaddumw
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
      roadorgw=roadorgw
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
    end if

  end do

  ! limit temperatures to sensible values     
  rooftemp(iqu,:)=min(max(roofdumtemp(:),225.),375.)
  walletemp(iqu,:)=min(max(walledumtemp(:),225.),375.)
  wallwtemp(iqu,:)=min(max(wallwdumtemp(:),225.),375.)
  roadtemp(iqu,:)=min(max(roaddumtemp(:),225.),375.)
  roofwater(iqu)=min(max(roofdumwat,0.),maxroofwater)
  roadwater(iqu)=min(max(roaddumwat,0.),maxroadwater)

  ! calculate outputs
  ofg(iqu)=fnsigmabld(iqu)*rooffg+(1.-fnsigmabld(iqu))*topfg+fnindustryfg(iqu)
  oeg(iqu)=fnsigmabld(iqu)*roofeg+(1.-fnsigmabld(iqu))*topeg
  ots(iqu)=fnsigmabld(iqu)*rooftemp(iqu,1)+(1.-fnsigmabld(iqu))*canyontemp !MJT - since this is what the atmosphere can 'see'
  owf(iqu)=fnsigmabld(iqu)*roofdelta+(1.-fnsigmabld(iqu))*roaddelta

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
! Based on Eva's sflux.f (i.e., vegetation)
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

subroutine getswcoeff(wallesg,wallwsg,roadsg,wallpsi,roadpsi,hwratio,vangle,hangle)

implicit none

real, intent(in) :: hwratio,vangle,hangle
real, intent(out) :: wallesg,wallwsg,roadsg,wallpsi,roadpsi
real thetazero,walles,wallws,roads,ta,tc,tz,xa,xb,ya,yb

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
  walles=((1./hwratio)*xa+ta*ya)/pi
  wallws=((1./hwratio)*xb+ta*yb)/pi
  roads=(tz-hwratio*ta*tc)/pi
end if

! note that these terms are truncated to 2nd order reflections, compared to TEB which uses infinte reflections.
wallesg=walles+roadalpha*wallpsi*roads+wallalpha*(1.-2.*wallpsi)*wallws+(wallalpha*(1.-2.*wallpsi))**2*walles &
        +roadalpha*wallalpha*wallpsi*(1.-roadpsi)*wallws+roadalpha*wallalpha*wallpsi*(1.-2.*wallpsi)*roads
wallwsg=wallws+roadalpha*wallpsi*roads+wallalpha*(1.-2.*wallpsi)*walles+(wallalpha*(1.-2.*wallpsi))**2*wallws &
        +roadalpha*wallalpha*wallpsi*(1.-roadpsi)*walles+roadalpha*wallalpha*wallpsi*(1.-2.*wallpsi)*roads
roadsg=roads+wallalpha*(1.-roadpsi)*0.5*(walles+wallws)+wallalpha*roadalpha*wallpsi*(1.-roadpsi)*roads &
        +wallalpha**2*(1.-roadpsi)*(1.-2.*wallpsi)*0.5*(walles+wallws)

return
end subroutine getswcoeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine solvecanyon(evct,wallefg,wallwfg,roadfg,topfg,rwinvres,topinvres,zo,zmin,ctemp &
  ,theta,umag,cu,walletemp,wallwtemp,roadtemp,rho,efftrafficfg,hwratio)

implicit none

real, intent(in) :: zo,zmin,ctemp,theta,umag,cu,walletemp,wallwtemp,roadtemp,rho,efftrafficfg,hwratio
real, intent(out) :: evct,wallefg,wallwfg,roadfg,topfg,rwinvres,topinvres
real cd,cw,rwi

call getinvres(topinvres,cd,zo,zmin,ctemp,theta,umag) 
cw=sqrt(cd)*umag ! diagnose canyonw (from Masson 2000)
rwi=11.8+4.2*sqrt(cu**2+cw**2) ! From Rowley, et al (1930)
rwinvres=rwi/(aircp*rho)
wallefg=(walletemp-ctemp)*rwi
wallwfg=(wallwtemp-ctemp)*rwi
roadfg=(roadtemp-ctemp)*rwi
topfg=aircp*rho*(ctemp-theta)*topinvres
evct=topfg-(roadfg+hwratio*(wallefg+wallwfg)+efftrafficfg)

return
end subroutine solvecanyon

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
