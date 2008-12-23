
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
!     call tebzo       ! blends input and urban roughness lengths for momentum and heat (use tebcd for
!                        blending the drag coefficent)
!     ...
!   end do
!   ...
!   call tebsavem    ! to save current state arrays (for use by tebloadm)
!   call tebend      ! to deallocate memory before quiting

! only tebinit and tebcalc are manditory.  All other subroutine calls are optional.

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
public atebinit,atebcalc,atebend,atebzo,atebload,atebsave,atebtype,atebfndef,atebalb,atebalb1, &
       atebnewangle,atebnewangle1,atebccangle,atebdisable,atebloadm,atebsavem,atebcd

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
  real, dimension(3) :: roofdepth,walldepth,roaddepth
  real, dimension(3) :: roofcp,wallcp,roadcp
  real, dimension(3) :: rooflambda,walllambda,roadlambda
  real :: hwratio,sigmabld,industryfg,trafficfg,bldheight,vangle,hangle
  real :: roofalpha,wallalpha,roadalpha,ctime,roofemiss,wallemiss,roademiss
  real :: bldtemp
end type tdata
type tprog
  real :: lzom,lzoh,cndzmin,cduv
end type tprog

! state arrays
integer, save :: ufull
integer, dimension(:), allocatable, save :: ugrid,mgrid
real, dimension(:), allocatable, save :: sigmau
type(tsurf), dimension(:), allocatable, save :: roof,road,roofadj,roadadj
type(tsurf), dimension(:), allocatable, save :: roofadj2,roadadj2,roofadj3,roadadj3
type(twall), dimension(:), allocatable, save :: walle,wallw,walleadj,wallwadj
type(twall), dimension(:), allocatable, save :: walleadj2,wallwadj2,walleadj3,wallwadj3
type(tdata), dimension(:), allocatable, save :: fn
type(tprog), dimension(:), allocatable, save :: pg
! model parameters
integer, parameter :: resmeth=1      ! Canyon sensible heat transfer (0=Masson, 1=Harman, 2=Kusaka)
integer, parameter :: zohmeth=1      ! Urban roughness length for heat (0=0.1*zom, 1=Kanda)
integer, parameter :: acmeth=1       ! AC heat pump into canyon (0=Off, 1=On)
integer, parameter :: nrefl=3        ! Number of canyon reflections (default=3)
integer, parameter :: nfgits=6       ! Maximum number of iterations for calculating sensible heat flux (default=6)
integer, parameter :: npgits=10      ! Maximum number of iterations for calculating prognostic variables (default=10 for dt=1200s)
integer, parameter :: stabfn=1       ! Stability function scheme (0=Louis, 1=Dyer and Hicks)
real, parameter :: tol=0.001         ! Minimum change for sectant method
real, parameter :: maxtol=0.1        ! Maximum change in temperature to terminate predictor-corrector loop (K)
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
! generic urban parameters
real, parameter :: refheight=0.4     ! Displacement height as a fraction of building height (Kanda et al 2007)
real, parameter :: zomratio=0.1      ! Roughness length to building height ratio (default=10%)
real, parameter :: zocanyon=0.1      ! Roughness length of canyon surfaces (m)
real, parameter :: maxrfwater=1.     ! Maximum roof water (kg m^-2)
real, parameter :: maxrdwater=1.     ! Maximum road water (kg m^-2)
real, parameter :: maxrfsn=1.        ! Maximum roof snow (kg m^-2)
real, parameter :: maxrdsn=1.        ! Maximum road snow (kg m^-2)

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
allocate(roofadj(ufull),roadadj(ufull),walleadj(ufull),wallwadj(ufull))
allocate(roofadj2(ufull),roadadj2(ufull),walleadj2(ufull),wallwadj2(ufull))
allocate(roofadj3(ufull),roadadj3(ufull),walleadj3(ufull),wallwadj3(ufull))
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
  roof%temp(ii)=291.
  road%temp(ii)=291.
  walle%temp(ii)=291.
  wallw%temp(ii)=291.
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
fn%industryfg=0.
fn%trafficfg=0.
fn%bldheight=10.
fn%roofalpha=0.2
fn%wallalpha=0.2
fn%roadalpha=0.2
fn%roofemiss=0.97
fn%wallemiss=0.97
fn%roademiss=0.97
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

return
end subroutine atebinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine deallocates arrays used by the TEB scheme

subroutine atebend(diag)

implicit none

integer, intent(in) :: diag

if (ufull.eq.0) return
if (diag.ge.1) write(6,*) "Deallocating aTEB arrays"

deallocate(ugrid,mgrid,fn)
deallocate(roof,road,walle,wallw)
deallocate(roofadj,roadadj,walleadj,wallwadj)
deallocate(roofadj2,roadadj2,walleadj2,wallwadj2)
deallocate(roofadj3,roadadj3,walleadj3,wallwadj3)
deallocate(sigmau)

return
end subroutine atebend

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine loads aTEB state arrays (not compulsory)

subroutine atebload(ifull,urban,diag)

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
end subroutine atebload

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! temperature only version of tebsave

subroutine atebloadm(ifull,urban,diag)

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
end subroutine atebloadm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine loads aTEB type arrays (not compulsory)

subroutine atebtype(ifull,itype,diag)

implicit none

integer, intent(in) :: ifull,diag
integer ii
integer, dimension(ifull), intent(in) :: itype
integer, parameter :: maxtype = 8
! Urban fraction (defined in host model)
!real, dimension(maxtype), parameter :: csigu=(/  0.6,  0.5,  0.6,  0.7, 0.95,  0.5,  0.6,  0.7 /)
! Building height (m)
real, dimension(maxtype), parameter :: cbldheight=(/ 10.,  4.,  6.,  8., 20.,  5., 10., 15. /)
! Building height to width ratio
real, dimension(maxtype), parameter ::   chwratio=(/  1., 0.4, 0.6, 0.8,  2., 0.5,  1., 1.5 /)
! Area fraction occupied by buildings
real, dimension(maxtype), parameter ::  csigmabld=(/ 0.5, 0.7, 0.7, 0.7, 0.4, 0.7, 0.5, 0.4 /)
! Industral sensible heat flux (W m^-2)
real, dimension(maxtype), parameter :: cindustryfg=(/ 20.,  8., 12., 16., 28., 20., 40., 60. /)
! Daily averaged traffic sensible heat flux (W m^-2)
real, dimension(maxtype), parameter ::  ctrafficfg=(/  5.,  2.,  3.,  4.,  7.,  5., 10., 15. /)
! Comfort temperature (K)
real, dimension(maxtype), parameter :: cbldtemp=(/ 291.16, 291.16, 291.16, 291.16, 291.16, 291.16, 291.16, 291.16 /)
! Roof albedo
real, dimension(maxtype), parameter :: croofalpha=(/ 0.20, 0.23, 0.20, 0.17, 0.13, 0.20, 0.20, 0.20 /)
! Wall albedo
real, dimension(maxtype), parameter :: cwallalpha=(/ 0.33, 0.37, 0.33, 0.29, 0.22, 0.33, 0.33, 0.33 /)
! Road albedo
real, dimension(maxtype), parameter :: croadalpha=(/ 0.10, 0.11, 0.10, 0.09, 0.07, 0.10, 0.10, 0.10 /)
! Roof emissitivity
real, dimension(maxtype), parameter :: croofemiss=(/ 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90 /)
! Wall emissitivity
real, dimension(maxtype), parameter :: cwallemiss=(/ 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85 /) 
! Road emissitivity
real, dimension(maxtype), parameter :: croademiss=(/ 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94 /)
! Roof depths (m)
real, dimension(maxtype,3), parameter :: croofdepth=reshape((/ 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, &
                                                               0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 0.40, &
                                                               0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10 /), &
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
                                                            0.28E6, 0.28E6, 0.28E6, 0.28E6, 0.28E6, 0.28E6, 0.28E6, 0.28E6, &
                                                            0.29E6, 0.29E6, 0.29E6, 0.29E6, 0.29E6, 0.29E6, 0.29E6, 0.29E6 /), &
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
                                                                0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, &
                                                                0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05 /), &
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

if (any(itype.lt.1).or.any(itype.gt.maxtype)) then
  write(6,*) "ERROR: Urban type is out of range"
  stop
end if

fn%hwratio=chwratio(itype(ugrid))
fn%sigmabld=csigmabld(itype(ugrid))
fn%industryfg=cindustryfg(itype(ugrid))
fn%trafficfg=ctrafficfg(itype(ugrid))
fn%bldheight=cbldheight(itype(ugrid))
fn%roofalpha=croofalpha(itype(ugrid))
fn%wallalpha=cwallalpha(itype(ugrid))
fn%roadalpha=croadalpha(itype(ugrid))
fn%roofemiss=croofemiss(itype(ugrid))
fn%wallemiss=cwallemiss(itype(ugrid))
fn%roademiss=croademiss(itype(ugrid))
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
real, dimension(ifull,39), intent(in) :: ifn

if (ufull.eq.0) return
if (diag.ge.1) write(6,*) "Load aTEB building properties"

fn%hwratio=ifn(ugrid,1)
fn%sigmabld=ifn(ugrid,2)
fn%industryfg=ifn(ugrid,3)
fn%trafficfg=ifn(ugrid,4)
fn%bldheight=ifn(ugrid,5)
fn%roofalpha=ifn(ugrid,6)
fn%wallalpha=ifn(ugrid,7)
fn%roadalpha=ifn(ugrid,8)
fn%roofemiss=ifn(ugrid,9)
fn%wallemiss=ifn(ugrid,10)
fn%roademiss=ifn(ugrid,11)
fn%bldtemp=ifn(ugrid,12)
do ii=1,3
  fn%roofdepth(ii)=ifn(ugrid,12+ii)
  fn%walldepth(ii)=ifn(ugrid,15+ii)
  fn%roaddepth(ii)=ifn(ugrid,18+ii)
  fn%roofcp(ii)=ifn(ugrid,21+ii)
  fn%wallcp(ii)=ifn(ugrid,24+ii)
  fn%roadcp(ii)=ifn(ugrid,27+ii)
  fn%rooflambda(ii)=ifn(ugrid,30+ii)
  fn%walllambda(ii)=ifn(ugrid,33+ii)
  fn%roadlambda(ii)=ifn(ugrid,36+ii)
end do

return
end subroutine atebfndef

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine saves aTEB state arrays (not compulsory)

subroutine atebsave(ifull,urban,diag)

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
end subroutine atebsave

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! temperature only version of tebsave

subroutine atebsavem(ifull,urban,diag)

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
end subroutine atebsavem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine blends urban momentum and heat roughness lengths
! (This version neglects the displacement height (e.g., for CCAM))
!

subroutine atebzo(ifull,zom,zoh,diag)

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
! (all grid points)

subroutine atebalb(ifull,alb,diag)

implicit none

integer, intent(in) :: ifull,diag
real, dimension(ifull), intent(inout) :: alb
real, dimension(ufull) :: ualb

if (ufull.eq.0) return
if (diag.ge.1) write(6,*) "Calculate urban albedo"

call atebalbcalc(1,ufull,ualb,diag)
alb(ugrid)=(1.-sigmau)*alb(ugrid)+sigmau*ualb

return
end subroutine atebalb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates the urban contrabution to albedo.
! (selected grid points only)

subroutine atebalb1(is,ifull,alb,diag)

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

call atebalbcalc(ib,ucount,ualb(1:ucount),diag)
alb(lgrid(1:ucount))=(1.-sigmau(ib:ie))*alb(lgrid(1:ucount)) &
                     +sigmau(ib:ie)*ualb(1:ucount)

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
albu=1.-(fn(is:ie)%hwratio*(sg%walle+sg%wallw)*(1.-fn(is:ie)%wallalpha) &
    +sg%road*((1.-snowdelta)*(1.-fn(is:ie)%roadalpha)+snowdelta*(1.-road(is:ie)%alpha)))
snowdelta=roof(is:ie)%snow/(roof(is:ie)%snow+maxrfsn)
albr=(1.-snowdelta)*fn(is:ie)%roofalpha+snowdelta*roof(is:ie)%alpha
alb=fn(is:ie)%sigmabld*albr+(1.-fn(is:ie)%sigmabld)*albu

return
end subroutine atebalbcalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine stores the zenith angle and the solar azimuth angle
! (all grid points)

! ifull = size of array
! cosin = array of cosine of zenith angles
! azimuthin = array of azimuthal angles
! diag = dialogue flag (0=off)

subroutine atebnewangle(ifull,cosin,azimuthin,ctimein,diag)

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
end subroutine atebnewangle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine stores the zenith angle and the solar azimuth angle
! (single grid point)

! iq = grid point
! cosin = cosine of zenith angle
! azmiuthin = azimuthal angle (rad)

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
! This version of tebnewangle is for CCAM
!

subroutine atebccangle(is,ifull,cosin,rlon,rlat,fjd,slag,dhr,dlt)

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
  fn(mgrid(is:ifull+is-1))%ctime=min(max(mod(hloc/(2.*pi)-0.5,1.),0.),1.)
end where

return
end subroutine atebccangle


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

subroutine atebcalc(ifull,ofg,oeg,ots,owf,dt,zmin,sg,rg,rnd,rho,temp,mixr,ps,pa,uu,vv,umin,diag)

implicit none

integer, intent(in) :: ifull,diag
real, intent(in) :: dt,zmin,umin
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

call atebeval(uo,dt,atm,zmin,diag)

ofg(ugrid)=(1.-sigmau)*ofg(ugrid)+sigmau*uo%fg
oeg(ugrid)=(1.-sigmau)*oeg(ugrid)+sigmau*uo%eg
ots(ugrid)=(1.-sigmau)*ots(ugrid)+sigmau*uo%ts
owf(ugrid)=(1.-sigmau)*owf(ugrid)+sigmau*uo%wf

return
end subroutine atebcalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! urban flux calculations

! Basic loop is:
!  Predictor/corrector loop
!    Heat flux through walls/roofs/roads
!    Short wave flux (nrefl reflections)
!    Estimate buld roughness length for momentum
!    Canyon aerodynamic resistances
!    Balance canyon snow temperature
!      Balance canyon sensible heat flux
!        Canyon sensible and latent heat fluxes
!      End sensible heat flux loop
!      Canyon long wave flux (nrefl reflections)
!    End canyon snow temprature loop
!    Roof sensible and latent heat fluxes
!    Roof long wave flux
!    Update water on canyon surfaces
!    Update snow albedo and density
!    Update urban temperatures
!  End predictor/corrector loop
!  Estimate bulk roughness length for heat
!  Estimate bulk long wave flux and surface temperature
!  Estimate bulk sensible and latent heat fluxes

subroutine atebeval(uo,ddt,atm,zmin,diag)

implicit none

integer, intent(in) :: diag
integer iqu,j,k,ii,cns,cnr
integer, dimension(ufull) :: igs,igr
integer, save :: callnumber=0
real, intent(in) :: ddt,zmin
real maxchange
real, dimension(ufull,3) :: garoof,gawalle,gawallw,garoad
real, dimension(ufull) :: garfsn,gardsn
real, dimension(ufull) :: rdsntemp,rfsntemp,rdsnmelt,rfsnmelt
real, dimension(ufull) :: wallpsi,roadpsi,fgtop,egtop,qsatr
real, dimension(ufull) :: oldval,newval,topu,cu,ctmax,ctmin,evctx,evct
real, dimension(ufull) :: ln,rn,we,ww,wr,zolog,a,xe,xw,cuven,n,zom
real, dimension(ufull) :: p_sntemp,p_netldratio,p_wallpsi,p_roadpsi
real, dimension(ufull) :: p_fgtop,p_gasn,p_snmelt
real, parameter :: temptol = 0.1
type(tatm), dimension(ufull), intent(in) :: atm
type(tatm), dimension(ufull) :: p_atm
type(tout), dimension(ufull), intent(out) :: uo
type(trad), dimension(ufull) :: sg,rg,fg,eg
type(trad), dimension(ufull) :: p_sg,p_rg,p_fg,p_eg
type(trad), dimension(ufull) :: acond,p_acond
type(tdiag), dimension(ufull) :: dg,p_dg
type(tsurf), dimension(ufull) :: roofdum,roaddum,croof,croad,rooforg,roadorg
type(tsurf), dimension(ufull) :: nroof,nroad,nroofold,nroadold,roofold,roadold
type(tsurf), dimension(ufull) :: p_rodum
type(twall), dimension(ufull) :: walledum,wallwdum,cwalle,cwallw,walleorg,wallworg
type(twall), dimension(ufull) :: nwalle,nwallw,nwalleold,nwallwold,walleold,wallwold
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
rg=acond
fg=acond
eg=acond
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

j=1
maxchange=tol+1.
do while ((j.le.npgits).and.(maxchange.gt.temptol)) ! predictor-corrector loop -------------------------------
  ! limit state variables
  do ii=1,3
    roofdum%temp(ii)=min(max(roofdum%temp(ii),200.),400.)
    walledum%temp(ii)=min(max(walledum%temp(ii),200.),400.)
    wallwdum%temp(ii)=min(max(wallwdum%temp(ii),200.),400.)
    roaddum%temp(ii)=min(max(roaddum%temp(ii),200.),400.)
  end do
  roofdum%water=min(max(roofdum%water,0.),maxrfwater)
  roaddum%water=min(max(roaddum%water,0.),maxrdwater)
  roofdum%snow=min(max(roofdum%snow,0.),maxrfsn)
  roaddum%snow=min(max(roaddum%snow,0.),maxrdsn)
  roofdum%den=min(max(roofdum%den,minsnowden),maxsnowden)
  roaddum%den=min(max(roaddum%den,minsnowden),maxsnowden)
  roofdum%alpha=min(max(roofdum%alpha,minsnowalpha),maxsnowalpha)
  roaddum%alpha=min(max(roaddum%alpha,minsnowalpha),maxsnowalpha)
    
  ! water and snow cover fractions
  dg%roofdelta=(roofdum%water/maxrfwater)**(2./3.)
  dg%roaddelta=(roaddum%water/maxrdwater)**(2./3.)
  dg%rfsndelta=roofdum%snow/(roofdum%snow+maxrfsn)
  dg%rdsndelta=roaddum%snow/(roaddum%snow+maxrdsn)

  ! Estimate urban roughness length
  zom=zomratio*fn%bldheight
  ! Adjust canyon roughness to include snow
  n=roaddum%snow/(roaddum%snow+maxrdsn+0.408*grav*zom)                 ! snow cover for urban roughness calc (Douville, et al 1995)
  zom=(1.-n)*zom+n*zosnow                                              ! blend urban and snow roughness length
  dg%rfdzmin=max(abs(zmin-fn%bldheight*(1.-refheight)),zocanyon+1.)    ! distance to roof displacement height
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
    garoof(:,k)=2.*(roofdum%temp(k)-roofdum%temp(k+1))/(fn%roofdepth(k)/fn%rooflambda(k)+fn%roofdepth(k+1)/fn%rooflambda(k+1))
    gawalle(:,k)=2.*(walledum%temp(k)-walledum%temp(k+1))/(fn%walldepth(k)/fn%walllambda(k)+fn%walldepth(k+1)/fn%walllambda(k+1))
    gawallw(:,k)=2.*(wallwdum%temp(k)-wallwdum%temp(k+1))/(fn%walldepth(k)/fn%walllambda(k)+fn%walldepth(k+1)/fn%walllambda(k+1))
    garoad(:,k)=2.*(roaddum%temp(k)-roaddum%temp(k+1))/(fn%roaddepth(k)/fn%roadlambda(k)+fn%roaddepth(k+1)/fn%roadlambda(k+1))
  end do
  garoof(:,3)=2.*fn%rooflambda(3)*(roofdum%temp(3)-fn%bldtemp)/fn%roofdepth(3)
  gawalle(:,3)=2.*fn%walllambda(3)*(walledum%temp(3)-fn%bldtemp)/fn%walldepth(3)
  gawallw(:,3)=2.*fn%walllambda(3)*(wallwdum%temp(3)-fn%bldtemp)/fn%walldepth(3)
  garoad(:,3)=0.
  if (acmeth.eq.1) then
    dg%accool=max(0.,garoof(:,3)+gawalle(:,3)+gawallw(:,3))
  else
    dg%accool=0.
  end if

  ! calculate shortwave radiation
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
      rn=max(0.,fn%bldheight-ln/1.5)
      ln=min(1.5*fn%bldheight,ln)
      cu=topu*exp(-0.9*sqrt(ln**2+(fn%bldheight-rn)**2)/fn%bldheight)
      a=0.15*max(1.,1.5*fn%hwratio)
      zolog=log(0.1*fn%bldheight/zocanyon)
      where (rn.gt.0.) ! recirculation starts on east wall
        n=min(max(0.1*fn%bldheight,rn),fn%bldheight)
        ! MJT suggestion (cuven is multipled by (fn%bldheight-rn) to avoid divide by zero)
        cuven=topu*((fn%bldheight-n*log(n/zocanyon)/(2.3+zolog))-(fn%bldheight-n)/(2.3+zolog))
        ! cuven is back to the correct units of m/s)
        cuven=max(cuven/fn%bldheight,cu*(1.-exp(-a*(1.-rn/fn%bldheight)))/a)
        xe=exp(-a*rn/fn%bldheight)
        we=(cuven+cu*(1.-xe)/a)
        xw=exp(-a/fn%hwratio)
        wr=cu*fn%hwratio*xe*(1.-xw)/a
        ww=cu*xe*xw*(1.-exp(-a))/a
      elsewhere ! recirculation starts on road
        cuven=topu*zolog/(2.3+zolog)
        n=max(1./fn%hwratio-3.,0.)
        xe=exp(-a*n)
        cuven=fn%bldheight*max(cuven*n,cu*(1.-xe)/a)
        xw=exp(-a*3.)
        wr=fn%hwratio*(cuven/fn%bldheight+cu*(1.-xw)/a)
        ! MJT suggestion
        cuven=topu*(1./0.9-(1.+0.1/0.9*zolog)/(2.3+zolog))
        xe=cu*xe*(1.-exp(-a))/a
        we=max(cuven,xe)
        ww=cu*xw*(1.-exp(-a))/a
      end where
      a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
      n=abs(atm%udir)/pi
      acond%road=a*wr                ! road bulk transfer
      acond%walle=a*(n*ww+(1.-n)*we) ! east wall bulk transfer
      acond%wallw=a*(n*we+(1.-n)*ww) ! west wall bulk transfer
      zolog=log(0.1*fn%bldheight/zosnow)
      a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
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
  ! includes solution to sensible heat flux and longwave radiation
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
                     p_wallwdum(1:cns)%temp(1))+5. ! max road snow temp
    ctmin(1:cns)=min(p_dg(1:cns)%tempc,p_rodum(1:cns)%temp(1),p_walledum(1:cns)%temp(1), &
                     p_wallwdum(1:cns)%temp(1))-5. ! min road snow temp
    p_sntemp(1:cns)=ctmax(1:cns)
    n(1:cns)=p_rodum(1:cns)%snow*waterden/p_rodum(1:cns)%den ! snow depth
    a(1:cns)=icelambda*(p_rodum(1:cns)%den/waterden)**1.88   ! snow lambda
    p_netldratio(1:cns)=0.5*(n(1:cns)/a(1:cns)+p_fn(1:cns)%roaddepth(1)/p_fn(1:cns)%roadlambda(1))
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
    do k=1,nfgits ! sectant
      evctx(1:cns)=evct(1:cns)
      call solverdsn(cns,evct(1:cns),p_rg(1:cns),p_fg(1:cns),p_fgtop(1:cns),p_eg(1:cns), &
             p_gasn(1:cns),p_snmelt(1:cns),p_sntemp(1:cns),p_rodum(1:cns),p_walledum(1:cns), &
             p_wallwdum(1:cns),p_dg(1:cns),p_sg(1:cns),p_atm(1:cns),p_netldratio(1:cns),ddt, &
             p_acond(1:cns),p_wallpsi(1:cns),p_roadpsi(1:cns),p_fn(1:cns),p_pg(1:cns))
      evctx(1:cns)=evct(1:cns)-evctx(1:cns)
      if (all(abs(evctx(1:cns)).le.tol)) exit
      where (abs(evctx(1:cns)).gt.tol)
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
    p_netldratio(1:cnr)=0.5*p_fn(1:cnr)%roaddepth(1)/p_fn(1:cnr)%roadlambda(1)
    call solverdsn(cnr,evct(1:cnr),p_rg(1:cnr),p_fg(1:cnr),p_fgtop(1:cnr),p_eg(1:cnr), &
           p_gasn(1:cnr),p_snmelt(1:cnr),p_sntemp(1:cnr),p_rodum(1:cnr),p_walledum(1:cnr), &
           p_wallwdum(1:cnr),p_dg(1:cnr),p_sg(1:cnr),p_atm(1:cnr),p_netldratio(1:cnr),ddt, &
           p_acond(1:cnr),p_wallpsi(1:cnr),p_roadpsi(1:cnr),p_fn(1:cnr),p_pg(1:cnr))
    rg(igr(1:cnr))=p_rg(1:cnr) ! unpack
    fg(igr(1:cnr))=p_fg(1:cnr)
    fgtop(igr(1:cnr))=p_fgtop(1:cnr)
    eg(igr(1:cnr))=p_eg(1:cnr)
    gardsn(igr(1:cnr))=0.
    rdsnmelt(igr(1:cnr))=0.
    rdsntemp(igr(1:cnr))=roaddum(igr(1:cnr))%temp(1)
    acond(igr(1:cnr))=p_acond(1:cnr)
    dg(igr(1:cnr))%canyonrgout=p_dg(1:cnr)%canyonrgout
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
  rfsntemp=roaddum%temp(1)
  acond%rfsn=0.
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
    p_rg(1:cns)=rg(igs(1:cns))
    p_fg(1:cns)=fg(igs(1:cns))
    p_eg(1:cns)=eg(igs(1:cns))
    p_atm(1:cns)=atm(igs(1:cns))
    p_acond(1:cns)=acond(igs(1:cns))
    p_fn(1:cns)=fn(igs(1:cns))
    ctmax(1:cns)=max(p_dg(1:cns)%tempr,p_rodum(1:cns)%temp(1))+5. ! max roof snow temp
    ctmin(1:cns)=min(p_dg(1:cns)%tempr,p_rodum(1:cns)%temp(1))-5. ! min roof snow temp
    p_sntemp(1:cns)=ctmax(1:cns)
    n(1:cns)=p_rodum(1:cns)%snow*waterden/p_rodum(1:cns)%den ! snow depth
    a(1:cns)=icelambda*(p_rodum(1:cns)%den/waterden)**1.88   ! snow lambda
    p_netldratio(1:cns)=0.5*(n(1:cns)/a(1:cns)+p_fn(1:cns)%roofdepth(1)/p_fn(1:cns)%rooflambda(1))
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
    do k=1,nfgits ! sectant
      evctx(1:cns)=evct(1:cns)
      call solverfsn(cns,evct(1:cns),p_rg(1:cns),p_fg(1:cns),p_eg(1:cns),p_gasn(1:cns),p_snmelt(1:cns),p_sntemp(1:cns), &
                     p_rodum(1:cns),p_dg(1:cns),p_sg(1:cns),p_atm(1:cns),p_netldratio(1:cns),p_acond(1:cns),ddt)
      evctx(1:cns)=evct(1:cns)-evctx(1:cns)
      if (all(abs(evctx(1:cns)).le.tol)) exit
      where (abs(evctx(1:cns)).gt.tol)
        newval(1:cns)=p_sntemp(1:cns)-evct(1:cns)*(p_sntemp(1:cns)-oldval(1:cns))/evctx(1:cns)
        oldval(1:cns)=p_sntemp(1:cns)
        p_sntemp(1:cns)=newval(1:cns)
      end where
    end do
    rfsntemp(igs(1:cns))=min(max(p_sntemp(1:cns),ctmin(1:cns)),ctmax(1:cns))
    rg(igs(1:cns))%rfsn=p_rg(1:cns)%rfsn ! unpack
    fg(igs(1:cns))%rfsn=p_fg(1:cns)%rfsn
    eg(igs(1:cns))%rfsn=p_eg(1:cns)%rfsn
    garfsn(igs(1:cns))=p_gasn(1:cns)
    rfsnmelt(igs(1:cns))=p_snmelt(1:cns)
    acond(igs(1:cns))%rfsn=p_acond(1:cns)%rfsn
  end if
  !---------------------------------------------------------------- 

  ! calculate roof sensible and latent heat fluxes (without snow)
  rg%roof=fn%roofemiss*(atm%rg-sbconst*roofdum%temp(1)**4)
  dg%roofrgout=sbconst*(dg%rfsndelta*snowemiss*rfsntemp**4+(1.-dg%rfsndelta)*fn%roofemiss*roof%temp(1)**4) &
               +(1.-dg%rfsndelta*snowemiss-(1.-dg%rfsndelta)*fn%roofemiss)*atm%rg
  a=log(dg%rfdzmin/zocanyon)
  ! n and xe are dummy variables for cd and lzohroof
  call getinvres(ufull,acond%roof,n,xe,a,dg%rfdzmin,roofdum%temp(1),dg%tempr,atm%umag,1) ! Assume zot=0.1*zom (i.e., Kanda et al 2007, small experiment)
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
              +dg%rfsndelta*(garfsn-garoof(:,1)))/(fn%roofcp(1)*fn%roofdepth(1))
  cwalle%temp(1)=(sg%walle+rg%walle-fg%walle-eg%walle-gawalle(:,1))/(fn%wallcp(1)*fn%walldepth(1))
  cwallw%temp(1)=(sg%wallw+rg%wallw-fg%wallw-eg%wallw-gawallw(:,1))/(fn%wallcp(1)*fn%walldepth(1))
  croad%temp(1)=((1.-dg%rdsndelta)*(sg%road+rg%road-fg%road-eg%road-garoad(:,1)) &
              +dg%rdsndelta*(gardsn-garoad(:,1)))/(fn%roadcp(1)*fn%roaddepth(1))
  do k=2,3
    croof%temp(k)=(garoof(:,k-1)-garoof(:,k))/(fn%roofcp(k)*fn%roofdepth(k))
    cwalle%temp(k)=(gawalle(:,k-1)-gawalle(:,k))/(fn%wallcp(k)*fn%walldepth(k))
    cwallw%temp(k)=(gawallw(:,k-1)-gawallw(:,k))/(fn%wallcp(k)*fn%walldepth(k))
    croad%temp(k)=(garoad(:,k-1)-garoad(:,k))/(fn%roadcp(k)*fn%roaddepth(k))
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
    if (callnumber.le.2) then
      if (callnumber.eq.0) then
        roofadj=croof
        walleadj=cwalle
        wallwadj=cwallw
        roadadj=croad
        roofadj2=croof
        walleadj2=cwalle
        wallwadj2=cwallw
        roadadj2=croad
        roofadj3=croof
        walleadj3=cwalle
        wallwadj3=cwallw
        roadadj3=croad
      end if
      do ii=1,3
        nroof%temp(ii)=-roofdum%temp(ii)+roof%temp(ii)+ddt*(1.5*croof%temp(ii)-0.5*roofadj%temp(ii))
        nwalle%temp(ii)=-walledum%temp(ii)+walle%temp(ii)+ddt*(1.5*cwalle%temp(ii)-0.5*walleadj%temp(ii))
        nwallw%temp(ii)=-wallwdum%temp(ii)+wallw%temp(ii)+ddt*(1.5*cwallw%temp(ii)-0.5*wallwadj%temp(ii))
        nroad%temp(ii)=-roaddum%temp(ii)+road%temp(ii)+ddt*(1.5*croad%temp(ii)-0.5*roadadj%temp(ii))
      end do
      nroof%water=-roofdum%water+roof%water+ddt*(1.5*croof%water-0.5*roofadj%water)
      nroad%water=-roaddum%water+road%water+ddt*(1.5*croad%water-0.5*roadadj%water)
      nroof%snow=-roofdum%snow+roof%snow+ddt*(1.5*croof%snow-0.5*roofadj%snow)
      nroad%snow=-roaddum%snow+road%snow+ddt*(1.5*croad%snow-0.5*roadadj%snow)
      nroof%den=-roofdum%den+roof%den+ddt*(1.5*croof%den-0.5*roofadj%den)
      nroad%den=-roaddum%den+road%den+ddt*(1.5*croad%den-0.5*roadadj%den)
      nroof%alpha=-roofdum%alpha+roof%alpha+ddt*(1.5*croof%alpha-0.5*roofadj%alpha)
      nroad%alpha=-roaddum%alpha+road%alpha+ddt*(1.5*croad%alpha-0.5*roadadj%alpha)
    else
      do ii=1,3
        nroof%temp(ii)=-roofdum%temp(ii)+roof%temp(ii)+ &
          (ddt/24.)*(55.*croof%temp(ii)-59.*roofadj%temp(ii)+37.*roofadj2%temp(ii)-9.*roofadj3%temp(ii))
        nwalle%temp(ii)=-walledum%temp(ii)+walle%temp(ii)+ &
          (ddt/24.)*(55.*cwalle%temp(ii)-59.*walleadj%temp(ii)+37.*walleadj2%temp(ii)-9.*walleadj3%temp(ii))
        nwallw%temp(ii)=-wallwdum%temp(ii)+wallw%temp(ii)+ &
          (ddt/24.)*(55.*cwallw%temp(ii)-59.*wallwadj%temp(ii)+37.*wallwadj2%temp(ii)-9.*wallwadj3%temp(ii))
        nroad%temp(ii)=-roaddum%temp(ii)+road%temp(ii)+ &
          (ddt/24.)*(55.*croad%temp(ii)-59.*roadadj%temp(ii)+37.*roadadj2%temp(ii)-9.*roadadj3%temp(ii))
      end do
      nroof%water=-roofdum%water+roof%water+(ddt/24.)*(55.*croof%water-59.*roofadj%water+37.*roofadj2%water-9.*roofadj3%water)
      nroad%water=-roaddum%water+road%water+(ddt/24.)*(55.*croad%water-59.*roadadj%water+37.*roadadj2%water-9.*roadadj3%water)
      nroof%snow=-roofdum%snow+roof%snow+(ddt/24.)*(55.*croof%snow-59.*roofadj%snow+37.*roofadj2%snow-9.*roofadj3%snow)
      nroad%snow=-roaddum%snow+road%snow+(ddt/24.)*(55.*croad%snow-59.*roadadj%snow+37.*roadadj2%snow-9.*roadadj3%snow)
      nroof%den=-roofdum%den+roof%den+(ddt/24.)*(55.*croof%den-59.*roofadj%den+37.*roofadj2%den-9.*roofadj3%den)
      nroad%den=-roaddum%den+road%den+(ddt/24.)*(55.*croad%den-59.*roadadj%den+37.*roadadj2%den-9.*roadadj3%den)
      nroof%alpha=-roofdum%alpha+roof%alpha+(ddt/24.)*(55.*croof%alpha-59.*roofadj%alpha+37.*roofadj2%alpha-9.*roofadj3%alpha)
      nroad%alpha=-roaddum%alpha+road%alpha+(ddt/24.)*(55.*croad%alpha-59.*roadadj%alpha+37.*roadadj2%alpha-9.*roadadj3%alpha)
    end if
    do ii=1,3
      n=nroof%temp(ii)+roofdum%temp(ii)
      roofold%temp(ii)=roofdum%temp(ii)
      roofdum%temp(ii)=n
      n=nwalle%temp(ii)+walledum%temp(ii)
      walleold%temp(ii)=walledum%temp(ii)
      walledum%temp(ii)=n
      n=nwallw%temp(ii)+wallwdum%temp(ii)
      wallwold%temp(ii)=wallwdum%temp(ii)
      wallwdum%temp(ii)=n
      n=nroad%temp(ii)+roaddum%temp(ii)    
      roadold%temp(ii)=roaddum%temp(ii)
      roaddum%temp(ii)=n
    end do
    n=nroof%water+roofdum%water      
    roofold%water=roofdum%water
    roofdum%water=n
    n=nroad%water+roaddum%water      
    roadold%water=roaddum%water
    roaddum%water=n
    n=nroof%snow+roofdum%snow
    roofold%snow=roofdum%snow
    roofdum%snow=n
    n=nroad%snow+roaddum%snow      
    roadold%snow=roaddum%snow
    roaddum%snow=n
    n=nroof%den+roofdum%den      
    roofold%den=roofdum%den
    roofdum%den=n
    n=nroad%den+roaddum%den      
    roadold%den=roaddum%den
    roaddum%den=n
    n=nroof%alpha+roofdum%alpha      
    roofold%alpha=roofdum%alpha
    roofdum%alpha=n
    n=nroad%alpha+roaddum%alpha      
    roadold%alpha=roaddum%alpha
    roaddum%alpha=n
    rooforg=croof
    walleorg=cwalle
    wallworg=cwallw
    roadorg=croad
  else ! corrector
    if (callnumber.le.2) then
      do ii=1,3
        nroof%temp(ii)=-roofdum%temp(ii)+roof%temp(ii)+(ddt/12.)*(5.*croof%temp(ii)+8.*rooforg%temp(ii)-roofadj%temp(ii))
        nwalle%temp(ii)=-walledum%temp(ii)+walle%temp(ii)+(ddt/12.)*(5.*cwalle%temp(ii)+8.*walleorg%temp(ii)-walleadj%temp(ii))
        nwallw%temp(ii)=-wallwdum%temp(ii)+wallw%temp(ii)+(ddt/12.)*(5.*cwallw%temp(ii)+8.*wallworg%temp(ii)-wallwadj%temp(ii))
        nroad%temp(ii)=-roaddum%temp(ii)+road%temp(ii)+(ddt/12.)*(5.*croad%temp(ii)+8.*roadorg%temp(ii)-roadadj%temp(ii))
      end do
      nroof%water=-roofdum%water+roof%water+(ddt/12.)*(5.*croof%water+8.*rooforg%water-roofadj%water)
      nroad%water=-roaddum%water+road%water+(ddt/12.)*(5.*croad%water+8.*roadorg%water-roadadj%water)
      nroof%snow=-roofdum%snow+roof%snow+(ddt/12.)*(5.*croof%snow+8.*rooforg%snow-roofadj%snow)
      nroad%snow=-roaddum%snow+road%snow+(ddt/12.)*(5.*croad%snow+8.*roadorg%snow-roadadj%snow)
      nroof%den=-roofdum%den+roof%den+(ddt/12.)*(5.*croof%den+8.*rooforg%den-roofadj%den)
      nroad%den=-roaddum%den+road%den+(ddt/12.)*(5.*croad%den+8.*roadorg%den-roadadj%den)
      nroof%alpha=-roofdum%alpha+roof%alpha+(ddt/12.)*(5.*croof%alpha+8.*rooforg%alpha-roofadj%alpha)
      nroad%alpha=-roaddum%alpha+road%alpha+(ddt/12.)*(5.*croad%alpha+8.*roadorg%alpha-roadadj%alpha)
    else
      do ii=1,3
        nroof%temp(ii)=-roofdum%temp(ii)+roof%temp(ii)+ &
          (ddt/24.)*(9.*croof%temp(ii)+19.*rooforg%temp(ii)-5.*roofadj%temp(ii)+roofadj2%temp(ii))
        nwalle%temp(ii)=-walledum%temp(ii)+walle%temp(ii)+ &
          (ddt/24.)*(9.*cwalle%temp(ii)+19.*walleorg%temp(ii)-5.*walleadj%temp(ii)+walleadj2%temp(ii))
        nwallw%temp(ii)=-wallwdum%temp(ii)+wallw%temp(ii)+ &
          (ddt/24.)*(9.*cwallw%temp(ii)+19.*wallworg%temp(ii)-5.*wallwadj%temp(ii)+wallwadj2%temp(ii))
        nroad%temp(ii)=-roaddum%temp(ii)+road%temp(ii)+ &
          (ddt/24.)*(9.*croad%temp(ii)+19.*roadorg%temp(ii)-5.*roadadj%temp(ii)+roadadj2%temp(ii))
      end do
      nroof%water=-roofdum%water+roof%water+(ddt/24.)*(9.*croof%water+19.*rooforg%water-5.*roofadj%water+roofadj2%water)
      nroad%water=-roaddum%water+road%water+(ddt/24.)*(9.*croad%water+19.*roadorg%water-5.*roadadj%water+roadadj2%water)
      nroof%snow=-roofdum%snow+roof%snow+(ddt/24.)*(9.*croof%snow+19.*rooforg%snow-5.*roofadj%snow+roofadj2%snow)
      nroad%snow=-roaddum%snow+road%snow+(ddt/24.)*(9.*croad%snow+19.*roadorg%snow-5.*roadadj%snow+roadadj2%snow)
      nroof%den=-roofdum%den+roof%den+(ddt/24.)*(9.*croof%den+19.*rooforg%den-5.*roofadj%den+roofadj2%den)
      nroad%den=-roaddum%den+road%den+(ddt/24.)*(9.*croad%den+19.*roadorg%den-5.*roadadj%den+roadadj2%den)
      nroof%alpha=-roofdum%alpha+roof%alpha+(ddt/24.)*(9.*croof%alpha+19.*rooforg%alpha-5.*roofadj%alpha+roofadj2%alpha)
      nroad%alpha=-roaddum%alpha+road%alpha+(ddt/24.)*(9.*croad%alpha+19.*roadorg%alpha-5.*roadadj%alpha+roadadj2%alpha)    
    end if
    maxchange=0.
    do ii=1,3
      where (abs(nroof%temp(ii)-nroofold%temp(ii)).gt.tol)
        n=roofdum%temp(ii)-nroof%temp(ii)*(roofdum%temp(ii)-roofold%temp(ii))/(nroof%temp(ii)-nroofold%temp(ii))
        roofold%temp(ii)=roofdum%temp(ii)
        roofdum%temp(ii)=n
      elsewhere
        roofold%temp(ii)=roofdum%temp(ii)
      end where
      maxchange=max(maxchange,maxval(abs(roofdum%temp(ii)-roofold%temp(ii))))
      where (abs(nwalle%temp(ii)-nwalleold%temp(ii)).gt.tol)
        n=walledum%temp(ii)-nwalle%temp(ii)*(walledum%temp(ii)-walleold%temp(ii))/(nwalle%temp(ii)-nwalleold%temp(ii))
        walleold%temp(ii)=walledum%temp(ii)
        walledum%temp(ii)=n
      elsewhere
        walleold%temp(ii)=walledum%temp(ii)
      end where
      maxchange=max(maxchange,maxval(abs(walledum%temp(ii)-walleold%temp(ii))))
      where (abs(nwallw%temp(ii)-nwallwold%temp(ii)).gt.tol)
        n=wallwdum%temp(ii)-nwallw%temp(ii)*(wallwdum%temp(ii)-wallwold%temp(ii))/(nwallw%temp(ii)-nwallwold%temp(ii))
        wallwold%temp(ii)=wallwdum%temp(ii)
        wallwdum%temp(ii)=n
      elsewhere
        wallwold%temp(ii)=wallwdum%temp(ii)
      end where
      maxchange=max(maxchange,maxval(abs(wallwdum%temp(ii)-wallwold%temp(ii))))
      where (abs(nroad%temp(ii)-nroadold%temp(ii)).gt.tol)
        n=roaddum%temp(ii)-nroad%temp(ii)*(roaddum%temp(ii)-roadold%temp(ii))/(nroad%temp(ii)-nroadold%temp(ii))
        roadold%temp(ii)=roaddum%temp(ii)
        roaddum%temp(ii)=n
      elsewhere
        roadold%temp(ii)=roaddum%temp(ii)
      end where
      maxchange=max(maxchange,maxval(abs(roaddum%temp(ii)-roadold%temp(ii))))
    end do
    where (abs(nroof%water-nroofold%water).gt.tol)
      n=roofdum%water-nroof%water*(roofdum%water-roofold%water)/(nroof%water-nroofold%water)
      roofold%water=roofdum%water
      roofdum%water=n
    elsewhere
      roofold%water=roofdum%water
    end where
    where (abs(nroad%water-nroadold%water).gt.tol)
      n=roaddum%water-nroad%water*(roaddum%water-roadold%water)/(nroad%water-nroadold%water)
      roadold%water=roaddum%water
      roaddum%water=n
    elsewhere
      roadold%water=roaddum%water
    end where
    where (abs(nroof%snow-nroofold%snow).gt.tol)
      n=roofdum%snow-nroof%snow*(roofdum%snow-roofold%snow)/(nroof%snow-nroofold%snow)
      roofold%snow=roofdum%snow
      roofdum%snow=n
    elsewhere
      roofold%snow=roofdum%snow
    end where
    where (abs(nroad%snow-nroadold%snow).gt.tol)
      n=roaddum%snow-nroad%snow*(roaddum%snow-roadold%snow)/(nroad%snow-nroadold%snow)
      roadold%snow=roaddum%snow
      roaddum%snow=n
    elsewhere
      roadold%snow=roaddum%snow
    end where
    where (abs(nroof%den-nroofold%den).gt.tol)
      n=roofdum%den-nroof%den*(roofdum%den-roofold%den)/(nroof%den-nroofold%den)
      roofold%den=roofdum%den
      roofdum%den=n
    elsewhere
      roofold%den=roofdum%den
    end where
    where (abs(nroad%den-nroadold%den).gt.tol)
      n=roaddum%den-nroad%den*(roaddum%den-roadold%den)/(nroad%den-nroadold%den)
      roadold%den=roaddum%den
      roaddum%den=n
    elsewhere
      roadold%den=roaddum%den
    end where
    where (abs(nroof%alpha-nroofold%alpha).gt.tol)
      n=roofdum%alpha-nroof%alpha*(roofdum%alpha-roofold%alpha)/(nroof%alpha-nroofold%alpha)
      roofold%alpha=roofdum%alpha
      roofdum%alpha=n
    elsewhere
      roofold%alpha=roofdum%alpha
    end where
    where (abs(nroad%alpha-nroadold%alpha).gt.tol)
      n=roaddum%alpha-nroad%alpha*(roaddum%alpha-roadold%alpha)/(nroad%alpha-nroadold%alpha)
      roadold%alpha=roaddum%alpha
      roaddum%alpha=n
    elsewhere
      roadold%alpha=roaddum%alpha
    end where
  end if

  nroofold=nroof
  nwalleold=nwalle
  nwallwold=nwallw
  nroadold=nroad

  j=j+1
end do

roofadj3=roofadj2
walleadj3=walleadj2
wallwadj3=wallwadj2
roadadj3=roadadj2
roofadj2=roofadj
walleadj2=walleadj
wallwadj2=wallwadj
roadadj2=roadadj
roofadj=rooforg
walleadj=walleorg
wallwadj=wallworg
roadadj=roadorg
callnumber=min(callnumber+1,3)

! limit temperatures to sensible values
do ii=1,3
  roof%temp(ii)=min(max(roofdum%temp(ii),200.),400.)
  walle%temp(ii)=min(max(walledum%temp(ii),200.),400.)
  wallw%temp(ii)=min(max(wallwdum%temp(ii),200.),400.)
  road%temp(ii)=min(max(roaddum%temp(ii),200.),400.)
end do
roof%water=min(max(roofdum%water,0.),maxrfwater)
road%water=min(max(roaddum%water,0.),maxrdwater)
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
! estimate surface temp from outgoing longwave radiation
uo%ts=((fn%sigmabld*dg%roofrgout+(1.-fn%sigmabld)*dg%canyonrgout)/sbconst)**0.25
uo%fg=fn%sigmabld*fg%roof+(1.-fn%sigmabld)*fgtop+fn%industryfg
uo%eg=fn%sigmabld*eg%roof+(1.-fn%sigmabld)*egtop
uo%wf=fn%sigmabld*dg%roofdelta*(1.-dg%rfsndelta)+(1.-fn%sigmabld)*dg%roaddelta*(1.-dg%rdsndelta)

! calculate roughness length for MOST
call getinvres(ufull,a,pg%cduv,pg%lzoh,pg%lzom,pg%cndzmin,uo%ts,dg%tempc,atm%umag,zohmeth+1)
  
return
end subroutine atebeval

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Louis parameterisation based on CCAM sflux.f version (i.e., from CSIRO9),
! but modified for increased ration between momentum and heat roughness lengths
! Also, included Dyer and Hicks scheme based on TAPM Surf.f

subroutine getinvres(cn,invres,cd,olzoh,ilzom,zmin,stemp,theta,umag,mode)

implicit none

integer, intent(in) :: cn,mode
integer ic
integer, parameter :: nc=5
real, dimension(cn), intent(in) :: ilzom,zmin,stemp,theta,umag
real, dimension(cn), intent(out) :: invres,cd,olzoh
real, dimension(cn) :: af,aft,ri,fm,fh,root,denma,denha,re,lna
real, dimension(cn) :: z_on_l,z0_on_l,zt_on_l,pm0,pm1,ph0,ph1
real, dimension(cn) :: integralm,integralh,pvstar
real, parameter :: bprm=5. ! 4.7 in rams
real, parameter :: chs=2.6 ! 5.3 in rams
real, parameter :: cms=5.  ! 7.4 in rams
real, parameter :: fmroot=0.57735
real, parameter :: rimax=(1./fmroot-1.)/bprm
real, parameter :: nu = 1.461E-5
real, parameter :: a_1 = 1.
real, parameter :: b_1 = 2./3.
real, parameter :: c_1 = 5.
real, parameter :: d_1 = 0.35
real, parameter :: aa1 = 3.8
real, parameter :: bb1 = 0.5
real, parameter :: cc1 = 0.3
!real, parameter :: eta0 = 1.827E-5
!real, parameter :: t0 = 291.15
!real, parameter :: c = 120.
!eta=eta0*((t0+c)/(theta+c))*(theta/t0)**(2./3.)
!nu=eta/rho

! use Louis as first guess for Dyer and Hicks scheme (if used)
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

select case(mode) ! roughness length for heat
  case(1) ! zot=zom/10.
    lna=2.3
  case(2) ! Kanda et al 2007
    re=max(sqrt(cd)*umag*zmin*exp(-ilzom)/nu,10.)
    !lna=2.46*re**0.25-2. !(Brutsaet, 1982)
    lna=1.29*re**0.25-2. !(Kanda et al, 2007)
  case(3) ! zot=zom (neglect molecular diffusion)
    lna=0.
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
 
 if (stabfn.eq.1) then ! from TAPM
  pvstar=aft*fh*(theta-stemp)/sqrt(cd)
  do ic=1,nc
    z_on_l=vkar*zmin*grav*pvstar/(theta*cd*umag**2)
    z_on_l=min(z_on_l,10.)
    z0_on_l  = z_on_l*exp(-ilzom)
    zt_on_l  = z0_on_l/exp(lna)
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
    where (z_on_l.le.0.4)
      cd = (max(0.01,min(vkar*umag/integralm,2.))/umag)**2
    elsewhere
      cd = (max(0.01,min(vkar*umag/(aa1*( ( z_on_l**bb1)*(1.0+cc1* z_on_l**(1.-bb1)) &
          -(z0_on_l**bb1)*(1.+cc1*z0_on_l**(1.-bb1)) )),2.))/umag)**2
    endwhere
    pvstar= vkar*(theta-stemp)/integralh
  end do
  invres=(vkar/integralh)*sqrt(cd)*umag
end if

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
! solve for road snow temperature (includes sensible heat flux)

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
real, dimension(cn) :: canyontemp,topinvres,trafficout
type(trad), dimension(cn), intent(inout) :: rg,fg,eg,acond
type(trad), dimension(cn), intent(in) :: sg
type(tatm), dimension(cn), intent(in) :: atm
type(tdiag), dimension(cn), intent(inout) :: dg
type(tsurf), dimension(cn), intent(in) :: iroad
type(twall), dimension(cn), intent(in) :: iwalle,iwallw
type(tdata), dimension(cn), intent(in) :: ifn
type(tprog), dimension(cn), intent(in) :: ipg

! traffic sensible heat flux
call gettraffic(cn,trafficout,ifn%trafficfg,ifn%ctime)
trafficout=trafficout/(1.-ifn%sigmabld)

! solve for canyon temperature ----------------------------------
ctmax=max(dg%tempc,iwalle%temp(1),iwallw%temp(1),iroad%temp(1),rdsntemp)+5. ! max canyon temp
ctmin=min(dg%tempc,iwalle%temp(1),iwallw%temp(1),iroad%temp(1),rdsntemp)-5. ! min canyon temp
call solvecanyon(cn,cevctx,fg,fgtop,topinvres,ctmax,dg,atm,iwalle%temp(1),iwallw%temp(1),iroad%temp(1) &
                ,rdsntemp,acond,ifn,ipg,trafficout)
canyontemp=0.5*(ctmax+ctmin)
call solvecanyon(cn,cevct,fg,fgtop,topinvres,canyontemp,dg,atm,iwalle%temp(1),iwallw%temp(1),iroad%temp(1) &
                ,rdsntemp,acond,ifn,ipg,trafficout)
where ((cevct*cevctx).lt.0.)
  ctmin=canyontemp
elsewhere
  ctmax=canyontemp
end where
oldval=canyontemp
canyontemp=0.5*(ctmax+ctmin)
do k=1,nfgits ! sectant
  cevctx=cevct
  call solvecanyon(cn,cevct,fg,fgtop,topinvres,canyontemp,dg,atm,iwalle%temp(1),iwallw%temp(1),iroad%temp(1) &
                  ,rdsntemp,acond,ifn,ipg,trafficout)
  cevctx=cevct-cevctx
  if (all(abs(cevctx).le.tol)) exit
  where (abs(cevctx).gt.tol)
    newval=canyontemp-cevct*(canyontemp-oldval)/cevctx
    oldval=canyontemp
    canyontemp=newval
  end where
end do
canyontemp=min(max(canyontemp,ctmin),ctmax)
! ---------------------------------------------------------------    

! note additional nrefl order reflections
netrad=dg%rdsndelta*snowemiss*rdsntemp**4+(1.-dg%rdsndelta)*ifn%roademiss*iroad%temp(1)**4
netemiss=dg%rdsndelta*snowemiss+(1.-dg%rdsndelta)*ifn%roademiss
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
  ncwa=(1.-netemiss)*wallpsi*rcra+(1.-ifn%wallemiss)*(1.-2.*wallpsi)*rcwa
  ncra=(1.-ifn%wallemiss)*(1.-roadpsi)*rcwa
  ncwe=(1.-netemiss)*wallpsi*rcrw+(1.-ifn%wallemiss)*(1.-2.*wallpsi)*rcww
  ncww=(1.-netemiss)*wallpsi*rcrw+(1.-ifn%wallemiss)*(1.-2.*wallpsi)*rcwe
  ncrw=(1.-ifn%wallemiss)*(1.-roadpsi)*0.5*(rcww+rcwe)  
  ncwr=(1.-netemiss)*wallpsi*rcrr+(1.-ifn%wallemiss)*(1.-2.*wallpsi)*rcwr
  ncrr=(1.-ifn%wallemiss)*(1.-roadpsi)*rcwr
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
    
rg%walle=ifn%wallemiss*(atm%rg*cwa+sbconst*iwalle%temp(1)**4*(-1.+ifn%wallemiss*cwe) & 
                  +sbconst*iwallw%temp(1)**4*ifn%wallemiss*cww+sbconst*netrad*cwr)
rg%wallw=ifn%wallemiss*(atm%rg*cwa+sbconst*iwallw%temp(1)**4*(-1.+ifn%wallemiss*cwe) &
                  +sbconst*iwalle%temp(1)**4*ifn%wallemiss*cww+sbconst*netrad*cwr)
rg%road=ifn%roademiss*(atm%rg*cra+sbconst*(-iroad%temp(1)**4+netrad*crr) &
                  +sbconst*ifn%wallemiss*(iwalle%temp(1)**4+iwallw%temp(1)**4)*crw)
where (dg%rdsndelta.gt.0.)
  rg%rdsn=snowemiss*(atm%rg*cra+sbconst*(-rdsntemp**4+netrad*crr) &
                    +sbconst*ifn%wallemiss*(iwalle%temp(1)**4+iwallw%temp(1)**4)*crw)
elsewhere
  rg%rdsn=0.
end where
dg%canyonrgout=atm%rg*(2.*ifn%hwratio*wallpsi*(1.-ifn%wallemiss)*cwa+roadpsi*(1.-netemiss)*cwr) &
               +sbconst*ifn%wallemiss*iwalle%temp(1)**4*(ifn%hwratio*wallpsi*(1.+(1.-ifn%wallemiss)*(cwe+cww)) &
               +roadpsi*(1.-netemiss)*crw) &
               +sbconst*ifn%wallemiss*iwallw%temp(1)**4*(ifn%hwratio*wallpsi*(1.+(1.-ifn%wallemiss)*(cwe+cww)) &
               +roadpsi*(1.-netemiss)*crw) &
               +sbconst*netrad*(2.*ifn%hwratio*wallpsi*(1.-ifn%wallemiss)*cwr+roadpsi*(1.+(1.-netemiss)*crr))

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

subroutine solvecanyon(cn,evct,fg,fgtop,topinvres,ctemp,dg,atm,walletemp,wallwtemp,roadtemp,rdsntemp,acond,ifn,ipg,trafficout)

implicit none

integer, intent(in) :: cn
real, dimension(cn), intent(out) :: evct,fgtop,topinvres
real, dimension(cn), intent(in) :: ctemp,walletemp,wallwtemp,roadtemp,rdsntemp,trafficout
real, dimension(cn) :: cduv,lzoh
type(trad), dimension(cn), intent(inout) :: fg,acond
type(tdiag), dimension(cn), intent(in) :: dg
type(tatm), dimension(cn), intent(in) :: atm
type(tdata), dimension(cn), intent(in) :: ifn
type(tprog), dimension(cn), intent(in) :: ipg

! Here we neglect the molecular diffusion and assume zot=zom.  This is similar to Kanada et al 2007, when the screen temp Te is used.
call getinvres(cn,topinvres,cduv,lzoh,ipg%lzom,ipg%cndzmin,ctemp,dg%tempc,atm%umag,3)
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
evct=fgtop-(dg%rdsndelta*fg%rdsn+(1.-dg%rdsndelta)*fg%road+ifn%hwratio*(fg%walle+fg%wallw) &
           +trafficout+dg%accool)

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
