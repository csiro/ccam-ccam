! Conformal Cubic Atmospheric Model
    
! Copyright 2015 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
! This file is part of the Conformal Cubic Atmospheric Model (CCAM)
!
! CCAM is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! CCAM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with CCAM.  If not, see <http://www.gnu.org/licenses/>.

!------------------------------------------------------------------------------   
    
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
!   call atebinit            ! to initalise state arrays, etc (use tebdisable to disable calls to ateb subroutines)
!   call atebloadm           ! to load previous state arrays (from tebsavem)
!   call atebtype            ! to define urban type (or use tebfndef to define urban properties at each grid point)
!   ...
!   do t=1,tmax
!     ...
!     call atebnewangle1     ! store current solar zenith and azimuthal angle (use atebccangle for CCAM)
!     call atebalb1(split=1) ! returns urban albedo for direct component of shortwave radiation
!     call atebalb1(split=2) ! returns urban albedo for diffuse component of shortwave radiation
!     ...
!     call atebfbeam         ! store fraction of direct shortwave radiation (or use atebspitter to estimate fraction)
!     call atebalb1          ! returns net urban albedo (i.e., default split=0)
!     ...
!     call atebcalc          ! calculates urban temperatures, fluxes, etc and blends with input
!     call atebcd            ! returns urban drag coefficent (or use atebzo for roughness length)
!     call atebscrnout       ! returns screen level diagnostics
!     ...
!   end do
!   ...
!   call atebsavem           ! to save current state arrays (for use by tebloadm)
!   call atebend             ! to deallocate memory before quiting

! only atebinit and atebcalc are manditory.  All other subroutine calls are optional.


! DEPENDICES:

! atebalb1(split=1)                   depends on     atebnewangle1 (or atebccangle)
! atebalb1(split=2)                   depends on     atebnewangle1 (or atebccangle)
! atebalb1 (i.e., default split=0)    depends on     atebnewangle1 (or atebccangle) and atebfbeam (or atebspitter)
! atebcalc                            depends on     atebnewangle1 (or atebccangle) and atebfbeam (or atebspitter)  
! atebcd (or atebzo)                  depends on     atebcalc
! atebscrnout                         depends on     atebcalc


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
       atebnewangle1,atebccangle,atebdisable,atebloadm,atebsavem,atebcd,vegmode,       &
       atebdwn,atebscrnout,atebfbeam,atebspitter,atebsigmau

! state arrays
integer, save :: ufull,ifull,iqut
logical, dimension(:), allocatable, save :: upack
real, dimension(:), allocatable, save :: sigmau
real, dimension(:,:), allocatable, save :: atebdwn ! These variables are for CCAM onthefly.f
real, dimension(:,:), allocatable, save :: f_roofdepth,f_walldepth,f_roaddepth
real, dimension(:,:), allocatable, save :: f_roofcp,f_wallcp,f_roadcp,f_vegcp
real, dimension(:,:), allocatable, save :: f_rooflambda,f_walllambda,f_roadlambda,f_veglambda
real, dimension(:), allocatable, save :: f_hwratio,f_sigmabld,f_industryfg,f_trafficfg,f_bldheight,f_vangle,f_hangle,f_fbeam
real, dimension(:), allocatable, save :: f_roofalpha,f_wallalpha,f_roadalpha,f_ctime,f_roofemiss,f_wallemiss,f_roademiss
real, dimension(:), allocatable, save :: f_bldtemp,f_sigmavegc,f_vegalphac,f_vegemissc
real, dimension(:), allocatable, save :: f_vegalphar,f_vegemissr,f_sigmavegr,f_vegdepthr
real, dimension(:), allocatable, save :: f_zovegc,f_vegrlaic,f_vegrsminc,f_zovegr,f_vegrlair,f_vegrsminr
real, dimension(:), allocatable, save :: f_swilt,f_sfc,f_ssat
real, dimension(:), allocatable, save :: p_lzom,p_lzoh,p_cndzmin,p_cduv,p_cdtq,p_vegtempc,p_vegtempr
real, dimension(:), allocatable, save :: p_tscrn,p_qscrn,p_uscrn,p_u10,p_emiss
real, dimension(:), allocatable, save :: p_roofskintemp,p_walleskintemp,p_wallwskintemp,p_roadskintemp

type roofroaddata
  real, dimension(:,:), allocatable :: temp
  real, dimension(:), allocatable   :: surfwater
  real, dimension(:), allocatable   :: leafwater
  real, dimension(:), allocatable   :: soilwater
  real, dimension(:), allocatable   :: snow
  real, dimension(:), allocatable   :: den
  real, dimension(:), allocatable   :: alpha
end type roofroaddata

type walldata
  real, dimension(:,:), allocatable :: temp
end type walldata

type(roofroaddata), save :: roof, road
type(walldata), save     :: walle, wallw


! model parameters
integer, parameter :: nmlfile=11      ! Read configuration from nml file (0=off, >0 unit number (default=11))
integer, save :: resmeth=1            ! Canyon sensible heat transfer (0=Masson, 1=Harman (varying width), 2=Kusaka,
                                      ! 3=Harman (fixed width))
integer, save :: useonewall=0         ! Combine both wall energy budgets into a single wall (0=two walls, 1=single wall) 
integer, save :: zohmeth=1            ! Urban roughness length for heat (0=0.1*zom, 1=Kanda, 2=0.003*zom)
integer, save :: acmeth=1             ! AC heat pump into canyon (0=Off, 1=On)
integer, save :: nrefl=3              ! Number of canyon reflections for radiation (default=3)
integer, save :: vegmode=2            ! In-canyon vegetation mode (0=50%/50%, 1=100%/0%, 2=0%/100%, where out/in=X/Y.
                                      ! Negative values are X=abs(vegmode))
integer, save :: scrnmeth=1           ! Screen diagnostic method (0=Slab, 1=Hybrid, 2=Canyon)
integer, save :: wbrelaxc=0           ! Relax canyon soil moisture for irrigation (0=Off, 1=On)
integer, save :: wbrelaxr=0           ! Relax roof soil moisture for irrigation (0=Off, 1=On)
integer, save :: iqt=314              ! Diagnostic point (in terms of host grid)
! sectant solver parameters
integer, save :: ncyits=6             ! Number of iterations for balancing canyon sensible and latent heat fluxes
integer, save :: nfgits=3             ! Number of iterations for balancing veg and snow energy budgets (default=3)
real, save    :: tol=0.001            ! Sectant method tolarance for sensible heat flux
real, save    :: alpha=1.             ! Weighting for determining the rate of convergence when calculating canyon temperatures
! physical parameters
real, parameter :: waterden=1000.     ! water density (kg m^-3)
real, parameter :: icelambda=2.22     ! conductance of ice (W m^-1 K^-1)
real, parameter :: aircp=1004.64      ! Heat capapcity of dry air (J kg^-1 K^-1)
real, parameter :: icecp=2100.        ! Heat capacity of ice (J kg^-1 K^-1)
real, parameter :: grav=9.80616       ! gravity (m s^-2)
real, parameter :: vkar=0.4           ! von Karman constant
real, parameter :: lv=2.501e6         ! Latent heat of vaporisation (J kg^-1)
real, parameter :: lf=3.337e5         ! Latent heat of fusion (J kg^-1)
real, parameter :: ls=lv+lf           ! Latent heat of sublimation (J kg^-1)
real, parameter :: pi=3.14159265      ! pi (must be rounded down for shortwave)
real, parameter :: rd=287.04          ! Gas constant for dry air
real, parameter :: rv=461.5           ! Gas constant for water vapor
real, parameter :: sbconst=5.67e-8    ! Stefan-Boltzmann constant
! snow parameters
real, save :: zosnow=0.001            ! Roughness length for snow (m)
real, save :: snowemiss=1.            ! snow emissitivity
real, save :: maxsnowalpha=0.85       ! max snow albedo
real, save :: minsnowalpha=0.5        ! min snow albedo
real, save :: maxsnowden=300.         ! max snow density (kg m^-3)
real, save :: minsnowden=100.         ! min snow density (kg m^-3)
! generic urban parameters
real, save :: refheight=0.6           ! Displacement height as a fraction of building height (Kanda et al 2007)
real, save :: zomratio=0.1            ! Ratio of roughness length to building height (default=0.1 or 10%)
real, save :: zocanyon=0.01           ! Roughness length of in-canyon surfaces (m)
real, save :: zoroof=0.01             ! Roughness length of roof surfaces (m)
real, save :: maxrfwater=1.           ! Maximum roof water (kg m^-2)
real, save :: maxrdwater=1.           ! Maximum road water (kg m^-2)
real, save :: maxrfsn=1.              ! Maximum roof snow (kg m^-2)
real, save :: maxrdsn=1.              ! Maximum road snow (kg m^-2)
real, save :: maxvwatf=0.1            ! Factor multiplied to LAI to predict maximum leaf water (kg m^-2)
real, save :: r_si=0.                 ! Building interior surface resistance (W^-1 m^2 K)
! atmosphere stability parameters
integer, save :: icmax=5              ! number of iterations for stability functions (default=5)
real, save    :: a_1=1.
real, save    :: b_1=2./3.
real, save    :: c_1=5.
real, save    :: d_1=0.35

interface getqsat
  module procedure getqsat_s,getqsat_v
end interface getqsat
interface getinvres
  module procedure getinvres_s,getinvres_v
end interface getinvres

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine prepare the arrays used by the aTEB scheme
! This is a compulsory subroutine that must be called during
! model initalisation

subroutine atebinit(ifin,sigu,diag)

implicit none

integer, intent(in) :: ifin,diag
integer, dimension(ifin) :: utype
integer iqu,iq
real, dimension(ifin), intent(in) :: sigu

if (diag>=1) write(6,*) "Initialising aTEB"

ifull=ifin
allocate(upack(ifull))
upack=sigu>0.
ufull=count(upack)
if (ufull==0) then
  deallocate(upack)
  return
end if

allocate(f_roofdepth(ufull,4),f_walldepth(ufull,4),f_roaddepth(ufull,4))
allocate(f_roofcp(ufull,4),f_wallcp(ufull,4),f_roadcp(ufull,4),f_vegcp(ufull,4))
allocate(f_rooflambda(ufull,4),f_walllambda(ufull,4),f_roadlambda(ufull,4),f_veglambda(ufull,4))
allocate(f_hwratio(ufull),f_sigmabld(ufull),f_industryfg(ufull),f_trafficfg(ufull),f_bldheight(ufull),f_vangle(ufull))
allocate(f_hangle(ufull),f_fbeam(ufull),f_roofalpha(ufull),f_wallalpha(ufull),f_roadalpha(ufull),f_ctime(ufull))
allocate(f_roofemiss(ufull),f_wallemiss(ufull),f_roademiss(ufull),f_bldtemp(ufull),f_sigmavegc(ufull),f_vegalphac(ufull))
allocate(f_vegemissc(ufull),f_sigmavegr(ufull),f_vegdepthr(ufull),f_vegalphar(ufull),f_vegemissr(ufull))
allocate(f_zovegc(ufull),f_vegrlaic(ufull),f_vegrsminc(ufull),f_zovegr(ufull),f_vegrlair(ufull),f_vegrsminr(ufull))
allocate(f_swilt(ufull),f_sfc(ufull),f_ssat(ufull))
allocate(p_lzom(ufull),p_lzoh(ufull),p_cndzmin(ufull),p_cduv(ufull),p_cdtq(ufull),p_vegtempc(ufull),p_vegtempr(ufull))
allocate(p_tscrn(ufull),p_qscrn(ufull),p_uscrn(ufull),p_u10(ufull),p_emiss(ufull))
allocate(p_roofskintemp(ufull),p_walleskintemp(ufull),p_wallwskintemp(ufull),p_roadskintemp(ufull))
allocate(roof%temp(ufull,4),roof%surfwater(ufull),roof%snow(ufull))
allocate(roof%den(ufull),roof%alpha(ufull))
allocate(road%temp(ufull,4),road%surfwater(ufull),road%snow(ufull))
allocate(road%den(ufull),road%alpha(ufull))
allocate(walle%temp(ufull,4),wallw%temp(ufull,4))
allocate(road%leafwater(ufull),road%soilwater(ufull),roof%leafwater(ufull),roof%soilwater(ufull))
allocate(sigmau(ufull))

! define grid arrays
iqu=0
do iq=1,ifull
  if (upack(iq)) then
    iqu=iqu+1
    sigmau(iqu)=sigu(iq)
  end if
end do

iqu=0
iqut=0
do iq=1,ifull
  if (upack(iq)) then
    iqu=iqu+1
    if (iq>=iqt) then
      iqut=iqu
      exit
    end if
  end if
end do

if (iqut==0) then
  !write(6,*) "WARN: Cannot located aTEB diagnostic point.  iqut=1"
  iqut=1
end if

! Initialise state variables
roof%temp=291.
road%temp=291.
walle%temp=291.
wallw%temp=291.
roof%surfwater=0.
roof%snow=0.
roof%den=minsnowden
roof%alpha=maxsnowalpha
road%surfwater=0.
road%snow=0.
road%den=minsnowden
road%alpha=maxsnowalpha

f_roofdepth=0.1
f_walldepth=0.1
f_roaddepth=0.1
f_vegdepthr=0.1
f_roofcp=2.E6
f_wallcp=2.E6
f_roadcp=2.E6
f_rooflambda=2.
f_walllambda=2.
f_roadlambda=2.
f_hwratio=1.
f_sigmabld=0.5
f_sigmavegc=0.5
f_sigmavegr=0.
f_industryfg=0.
f_trafficfg=0.
f_bldheight=10.
f_roofalpha=0.2
f_wallalpha=0.2
f_roadalpha=0.2
f_vegalphac=0.2
f_vegalphar=0.2
f_roofemiss=0.97
f_wallemiss=0.97
f_roademiss=0.97
f_vegemissc=0.97
f_vegemissr=0.97
f_bldtemp=291.
f_vangle=0.
f_hangle=0.
f_ctime=0.
f_fbeam=1.
f_zovegc=0.1
f_vegrlaic=1.
f_vegrsminc=200.
f_zovegr=0.1
f_vegrlair=1.
f_vegrsminr=200.
f_swilt=0.
f_sfc=0.5
f_ssat=1.

utype=1 ! default urban
call atebtype(utype,diag)

p_cndzmin=max(10.,0.1*f_bldheight+2.)   ! updated in atebcalc
p_lzom=log(p_cndzmin/(0.1*f_bldheight)) ! updated in atebcalc
p_lzoh=6.+p_lzom ! (Kanda et al 2005)   ! updated in atebcalc
p_cduv=(vkar/p_lzom)**2                 ! updated in atebcalc
p_cdtq=vkar**2/(p_lzom*p_lzoh)          ! updated in atebcalc
p_vegtempc=291.                         ! updated in atebcalc
p_vegtempr=291.                         ! updated in atebcalc
p_tscrn=291.                            ! updated in atebcalc
p_qscrn=0.                              ! updated in atebcalc
p_uscrn=0.                              ! updated in atebcalc
p_u10=0.                                ! updated in atebcalc
p_emiss=0.97                            ! updated in atebcalc
p_roofskintemp=291.                     ! updated in atebcalc
p_walleskintemp=291.                    ! updated in atebcalc
p_wallwskintemp=291.                    ! updated in atebcalc
p_roadskintemp=291.                     ! updated in atebcalc

road%soilwater=0.5*(f_ssat+f_swilt)
road%leafwater=0.
roof%soilwater=f_swilt
roof%leafwater=0.

return
end subroutine atebinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine deallocates arrays used by the TEB scheme

subroutine atebend(diag)

implicit none

integer, intent(in) :: diag

if (ufull==0) return
if (diag>=1) write(6,*) "Deallocating aTEB arrays"

deallocate(upack)
deallocate(f_roofdepth,f_walldepth,f_roaddepth)
deallocate(f_roofcp,f_wallcp,f_roadcp,f_vegcp)
deallocate(f_rooflambda,f_walllambda,f_roadlambda,f_veglambda)
deallocate(f_hwratio,f_sigmabld,f_industryfg,f_trafficfg,f_bldheight,f_vangle)
deallocate(f_hangle,f_fbeam,f_roofalpha,f_wallalpha,f_roadalpha,f_ctime)
deallocate(f_roofemiss,f_wallemiss,f_roademiss,f_bldtemp,f_sigmavegc,f_vegalphac)
deallocate(f_vegemissc,f_sigmavegr,f_vegdepthr,f_vegalphar,f_vegemissr)
deallocate(f_zovegc,f_vegrlaic,f_vegrsminc,f_zovegr,f_vegrlair,f_vegrsminr)
deallocate(f_swilt,f_sfc,f_ssat)
deallocate(p_lzom,p_lzoh,p_cndzmin,p_cduv,p_cdtq,p_vegtempc,p_vegtempr)
deallocate(p_tscrn,p_qscrn,p_uscrn,p_u10,p_emiss)
deallocate(p_roofskintemp,p_walleskintemp,p_wallwskintemp,p_roadskintemp)
deallocate(roof%temp,roof%surfwater,roof%snow)
deallocate(roof%den,roof%alpha)
deallocate(road%temp,road%surfwater,road%snow)
deallocate(road%den,road%alpha)
deallocate(walle%temp,wallw%temp)
deallocate(sigmau,road%leafwater,road%soilwater,roof%leafwater,roof%soilwater)

return
end subroutine atebend

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine loads aTEB state arrays (not compulsory)

subroutine atebload(urban,diag)

implicit none

integer, intent(in) :: diag
integer ii
real, dimension(ifull,28), intent(in) :: urban

if (ufull==0) return
if (diag>=1) write(6,*) "Load aTEB state arrays"

do ii=1,4
  roof%temp(:,ii) =pack(urban(:,ii),   upack)
  walle%temp(:,ii)=pack(urban(:,ii+4), upack)
  wallw%temp(:,ii)=pack(urban(:,ii+8), upack)
  road%temp(:,ii) =pack(urban(:,ii+12),upack)
end do
road%soilwater=pack(urban(:,17),upack)
roof%soilwater=pack(urban(:,18),upack)
roof%surfwater=pack(urban(:,19),upack)
road%surfwater=pack(urban(:,20),upack)
road%leafwater=pack(urban(:,21),upack)
roof%leafwater=pack(urban(:,22),upack)
roof%snow     =pack(urban(:,23),upack)
road%snow     =pack(urban(:,24),upack)
roof%den      =pack(urban(:,25),upack)
road%den      =pack(urban(:,26),upack)
roof%alpha    =pack(urban(:,27),upack)
road%alpha    =pack(urban(:,28),upack)

return
end subroutine atebload

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! temperature only version of tebsave

subroutine atebloadm(urban,moist,diag)

implicit none

integer, intent(in) :: diag
integer ii
real, dimension(ifull,16), intent(in) :: urban
real, dimension(ifull,2), intent(in) :: moist

if (ufull==0) return
if (diag>=1) write(6,*) "Load aTEB state arrays"

do ii=1,4
  roof%temp(:,ii) =pack(urban(:,ii),   upack)
  walle%temp(:,ii)=pack(urban(:,ii+4), upack)
  wallw%temp(:,ii)=pack(urban(:,ii+8), upack)
  road%temp(:,ii) =pack(urban(:,ii+12),upack)
end do
road%soilwater=pack(moist(:,1),upack)
roof%soilwater=pack(moist(:,2),upack)

return
end subroutine atebloadm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine loads aTEB type arrays (not compulsory)

subroutine atebtype(itype,diag)

implicit none

integer, intent(in) :: diag
integer ii,ierr
integer, dimension(ifull), intent(in) :: itype
integer, dimension(ufull) :: itmp
integer, parameter :: maxtype = 8
real x
real, dimension(ufull) :: tsigveg,tsigmabld
! In-canyon vegetation fraction
real, dimension(maxtype) ::    csigvegc=(/ 0.38, 0.45, 0.38, 0.34, 0.05, 0.40, 0.30, 0.20 /)
! Green roof vegetation fraction
real, dimension(maxtype) ::    csigvegr=(/ 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /)
! Area fraction occupied by buildings
real, dimension(maxtype) ::   csigmabld=(/ 0.45, 0.40, 0.45, 0.46, 0.65, 0.40, 0.45, 0.50 /)
! Building height (m)
real, dimension(maxtype) ::  cbldheight=(/   6.,   4.,   6.,   8.,  18.,   4.,   8.,  12. /)
! Building height to width ratio
real, dimension(maxtype) ::    chwratio=(/  0.4,  0.2,  0.4,  0.6,   2.,  0.5,   1.,  1.5 /)
! Industral sensible heat flux (W m^-2)
real, dimension(maxtype) :: cindustryfg=(/   0.,   0.,   0.,   0.,   0.,  10.,  20.,  30. /)
! Daily averaged traffic sensible heat flux (W m^-2)
real, dimension(maxtype) ::  ctrafficfg=(/  1.5,  1.5,  1.5,  1.5,  1.5,  1.5,  1.5,  1.5 /)
! Comfort temperature (K)
real, dimension(maxtype) :: cbldtemp=(/ 291.16, 291.16, 291.16, 291.16, 291.16, 291.16, 291.16, 291.16 /)
! Roof albedo
real, dimension(maxtype) ::  croofalpha=(/ 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20 /)
! Wall albedo
real, dimension(maxtype) ::  cwallalpha=(/ 0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30 /)
! Road albedo
real, dimension(maxtype) ::  croadalpha=(/ 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10 /)
! Canyon veg albedo
real, dimension(maxtype) ::  cvegalphac=(/ 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20 /)
! Roof veg albedo
real, dimension(maxtype) ::  cvegalphar=(/ 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20 /)
! Roof emissitivity
real, dimension(maxtype) ::  croofemiss=(/ 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90 /)
! Wall emissitivity
real, dimension(maxtype) ::  cwallemiss=(/ 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85 /) 
! Road emissitivity
real, dimension(maxtype) ::  croademiss=(/ 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94 /)
! Canyon veg emissitivity
real, dimension(maxtype) ::  cvegemissc=(/ 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96 /)
! Roof veg emissitivity
real, dimension(maxtype) ::  cvegemissr=(/ 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96 /)
! Green roof soil depth
real, dimension(maxtype) ::   cvegdeptr=(/ 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10 /)
! Roof depths (m)
real, dimension(maxtype,4) :: croofdepth=reshape((/ 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, &
                                                    0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, &
                                                    0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 0.40, &
                                                    0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10 /), &
                                                    (/maxtype,4/))
! Wall depths (m)
real, dimension(maxtype,4) :: cwalldepth=reshape((/ 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, &
                                                    0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, &
                                                    0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, &
                                                    0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05 /), &
                                                    (/maxtype,4/))
! Road depths (m)
real, dimension(maxtype,4) :: croaddepth=reshape((/ 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, &
                                                    0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, &
                                                    0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45, &
                                                    3.50, 3.50, 3.50, 3.50, 3.50, 3.50, 3.50, 3.50 /), &
                                                    (/maxtype,4/))
! Roof heat capacity (J m^-3 K^-1)
real, dimension(maxtype,4) :: croofcp=reshape((/ 2.11E6, 2.11E6, 2.11E6, 2.11E6, 2.11E6, 2.11E6, 2.11E6, 2.11E6, &
                                                 2.11E6, 2.11E6, 2.11E6, 2.11E6, 2.11E6, 2.11E6, 2.11E6, 2.11E6, &
                                                 0.28E6, 0.28E6, 0.28E6, 0.28E6, 0.28E6, 0.28E6, 0.28E6, 0.28E6, &
                                                 0.29E6, 0.29E6, 0.29E6, 0.29E6, 0.29E6, 0.29E6, 0.29E6, 0.29E6 /), &
                                                 (/maxtype,4/))
! Wall heat capacity (J m^-3 K^-1)
real, dimension(maxtype,4) :: cwallcp=reshape((/ 1.55E6, 1.55E6, 1.55E6, 1.55E6, 1.55E6, 1.55E6, 1.55E6, 1.55E6, &
                                                 1.55E6, 1.55E6, 1.55E6, 1.55E6, 1.55E6, 1.55E6, 1.55E6, 1.55E6, &
                                                 1.55E6, 1.55E6, 1.55E6, 1.55E6, 1.55E6, 1.55E6, 1.55E6, 1.55E6, &
                                                 0.29E6, 0.29E6, 0.29E6, 0.29E6, 0.29E6, 0.29E6, 0.29E6, 0.29E6 /), &
                                                 (/maxtype,4/))
! Road heat capacity (J m^-3 K^-1)
real, dimension(maxtype,4) :: croadcp=reshape((/ 1.94E6, 1.94E6, 1.94E6, 1.94E6, 1.94E6, 1.94E6, 1.94E6, 1.94E6, &
                                                 1.94E6, 1.94E6, 1.94E6, 1.94E6, 1.94E6, 1.94E6, 1.94E6, 1.94E6, &
                                                 1.28E6, 1.28E6, 1.28E6, 1.28E6, 1.28E6, 1.28E6, 1.28E6, 1.28E6, &
                                                 1.28E6, 1.28E6, 1.28E6, 1.28E6, 1.28E6, 1.28E6, 1.28E6, 1.28E6 /), &
                                                 (/maxtype,4/))
! Roof conductance (W m^-1 K^-1)
real, dimension(maxtype,4) :: crooflambda=reshape((/ 1.5100, 1.5100, 1.5100, 1.5100, 1.5100, 1.5100, 1.5100, 1.5100, &
                                                     1.5100, 1.5100, 1.5100, 1.5100, 1.5100, 1.5100, 1.5100, 1.5100, &
                                                     0.0800, 0.0800, 0.0800, 0.0800, 0.0800, 0.0800, 0.0800, 0.0800, &
                                                     0.0500, 0.0500, 0.0500, 0.0500, 0.0500, 0.0500, 0.0500, 0.0500 /), &
                                                     (/maxtype,4/))
! Wall conductance (W m^-1 K^-1)
real, dimension(maxtype,4) :: cwalllambda=reshape((/ 0.9338, 0.9338, 0.9338, 0.9338, 0.9338, 0.9338, 0.9338, 0.9338, &
                                                     0.9338, 0.9338, 0.9338, 0.9338, 0.9338, 0.9338, 0.9338, 0.9338, &
                                                     0.9338, 0.9338, 0.9338, 0.9338, 0.9338, 0.9338, 0.9338, 0.9338, &
                                                     0.0500, 0.0500, 0.0500, 0.0500, 0.0500, 0.0500, 0.0500, 0.0500 /), &
                                                     (/maxtype,4/))
! Road conductance (W m^-1 K^-1)
real, dimension(maxtype,4) :: croadlambda=reshape((/ 0.7454, 0.7454, 0.7454, 0.7454, 0.7454, 0.7454, 0.7454, 0.7454, &
                                                     0.7454, 0.7454, 0.7454, 0.7454, 0.7454, 0.7454, 0.7454, 0.7454, &
                                                     0.2513, 0.2513, 0.2513, 0.2513, 0.2513, 0.2513, 0.2513, 0.2513, &
                                                     0.2513, 0.2513, 0.2513, 0.2513, 0.2513, 0.2513, 0.2513, 0.2513 /), &
                                                     (/maxtype,4/))
! Roughness length of in-canyon vegetation (m)
real, dimension(maxtype) ::    czovegc=(/   0.1,   0.1,   0.1,   0.1,   0.1,   0.1,   0.1,   0.1 /)
! In-canyon vegetation LAI
real, dimension(maxtype) ::  cvegrlaic=(/   2.0,   2.0,   2.0,   2.0,   2.0,   2.0,   2.0,   2.0 /)
! Unconstrained canopy stomatal resistance
real, dimension(maxtype) :: cvegrsminc=(/ 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0 /)
! Roughness length of green roof vegetation (m)
real, dimension(maxtype) ::    czovegr=(/   0.1,   0.1,   0.1,   0.1,   0.1,   0.1,   0.1,   0.1 /)
! Green roof vegetation LAI
real, dimension(maxtype) ::  cvegrlair=(/   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0 /)
! Unconstrained canopy stomatal resistance
real, dimension(maxtype) :: cvegrsminr=(/ 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0 /)
! Soil wilting point (m^3 m^-3)
real, dimension(maxtype) ::     cswilt=(/  0.18,  0.18,  0.18,  0.18,  0.18,  0.18,  0.18,  0.18 /)
! Soil field capacity (m^3 m^-3)
real, dimension(maxtype) ::       csfc=(/  0.26,  0.26,  0.26,  0.26,  0.26,  0.26,  0.26,  0.26 /)
! Soil saturation point (m^3 m^-3)
real, dimension(maxtype) ::      cssat=(/  0.42,  0.42,  0.42,  0.42,  0.42,  0.42,  0.42,  0.42 /)
                                                     
namelist /atebnml/  resmeth,useonewall,zohmeth,acmeth,nrefl,vegmode,scrnmeth,wbrelaxc,wbrelaxr,iqt
namelist /atebsnow/ zosnow,snowemiss,maxsnowalpha,minsnowalpha,maxsnowden,minsnowden
namelist /atebgen/  refheight,zomratio,zocanyon,zoroof,maxrfwater,maxrdwater,maxrfsn,maxrdsn,maxvwatf
namelist /atebtile/ czovegc,cvegrlaic,cvegrsminc,czovegr,cvegrlair,cvegrsminr,cswilt,csfc,cssat,       &
                    cvegemissc,cvegemissr,cvegdeptr,cvegalphac,cvegalphar,csigvegc,csigvegr,           &
                    csigmabld,cbldheight,chwratio,cindustryfg,ctrafficfg,cbldtemp,croofalpha,          &
                    cwallalpha,croadalpha,croofemiss,cwallemiss,croademiss,croofdepth,cwalldepth,      &
                    croaddepth,croofcp,cwallcp,croadcp,crooflambda,cwalllambda,croadlambda
                                                                
if (ufull==0) return
if (diag>=1) write(6,*) "Load aTEB building properties"

itmp=pack(itype,upack)
if ((minval(itmp)<1).or.(maxval(itmp)>maxtype)) then
  write(6,*) "ERROR: Urban type is out of range"
  stop
end if

if (nmlfile/=0) then
  open(unit=nmlfile,file='ateb.nml',action="read",iostat=ierr)
  if (ierr==0) then
    write(6,*) "Reading ateb.nml"
    read(nmlfile,nml=atebnml,iostat=ierr)
    read(nmlfile,nml=atebsnow,iostat=ierr)
    read(nmlfile,nml=atebgen,iostat=ierr)
    read(nmlfile,nml=atebtile,iostat=ierr)
    close(nmlfile)  
  end if
end if

select case(vegmode)
  case(0)
    tsigveg=0.5*csigvegc(itmp)/(1.-0.5*csigvegc(itmp))
    tsigmabld=csigmabld(itmp)/(1.-0.5*csigvegc(itmp))
    sigmau=sigmau*(1.-0.5*csigvegc(itmp))
  case(1)
    tsigveg=0.
    tsigmabld=csigmabld(itmp)/(1.-csigvegc(itmp))
    sigmau=sigmau*(1.-csigvegc(itmp))
  case(2)
    tsigveg=csigvegc(itmp)
    tsigmabld=csigmabld(itmp)
  case DEFAULT
    if (vegmode<0) then
      x=real(abs(vegmode))/100.
      x=max(min(x,1.),0.)
      tsigveg=x*csigvegc(itmp)/(1.-(1.-x)*csigvegc(itmp))
      tsigmabld=csigmabld(itmp)/(1.-(1.-x)*csigvegc(itmp))
      sigmau=sigmau*(1.-(1.-x)*csigvegc(itmp))
    else
      write(6,*) "ERROR: Unsupported vegmode ",vegmode
      stop
    end if
end select
f_sigmavegc=max(min(tsigveg/(1.-tsigmabld),1.),0.)
f_sigmavegr=max(min(csigvegr(itmp),1.),0.)
f_sigmabld=max(min(tsigmabld,1.),0.)
f_hwratio=chwratio(itmp)*f_sigmabld/(1.-f_sigmabld) ! MJT suggested new definition
f_industryfg=cindustryfg(itmp)
f_trafficfg=ctrafficfg(itmp)
f_bldheight=cbldheight(itmp)
f_roofalpha=croofalpha(itmp)
f_wallalpha=cwallalpha(itmp)
f_roadalpha=croadalpha(itmp)
f_vegalphac=cvegalphac(itmp)
f_vegalphar=cvegalphar(itmp)
f_roofemiss=croofemiss(itmp)
f_wallemiss=cwallemiss(itmp)
f_roademiss=croademiss(itmp)
f_vegemissc=cvegemissc(itmp)
f_vegemissr=cvegemissr(itmp)
f_bldtemp=cbldtemp(itmp)
f_vegdepthr=cvegdeptr(itmp)
do ii=1,4
  f_roofdepth(:,ii)=croofdepth(itmp,ii)
  f_walldepth(:,ii)=cwalldepth(itmp,ii)
  f_roaddepth(:,ii)=croaddepth(itmp,ii)
  f_roofcp(:,ii)=croofcp(itmp,ii)
  f_wallcp(:,ii)=cwallcp(itmp,ii)
  f_roadcp(:,ii)=croadcp(itmp,ii)
  f_rooflambda(:,ii)=crooflambda(itmp,ii)
  f_walllambda(:,ii)=cwalllambda(itmp,ii)
  f_roadlambda(:,ii)=croadlambda(itmp,ii)
end do
f_zovegc=czovegc(itmp)
f_vegrlaic=cvegrlaic(itmp)
f_vegrsminc=cvegrsminc(itmp)/max(f_vegrlaic,1.E-8)
f_zovegr=czovegr(itmp)
f_vegrlair=cvegrlair(itmp)
f_vegrsminr=cvegrsminr(itmp)/max(f_vegrlair,1.E-8)
f_swilt=cswilt(itmp)
f_sfc=csfc(itmp)
f_ssat=cssat(itmp)

return
end subroutine atebtype

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine specifies the urban properties for each grid point
!

subroutine atebfndef(ifn,diag)

implicit none

integer, intent(in) :: diag
integer ii
real, dimension(ifull,63), intent(in) :: ifn

if (ufull==0) return
if (diag>=1) write(6,*) "Load aTEB building properties"

f_hwratio   =pack(ifn(:,1),upack)
f_sigmabld  =pack(ifn(:,2),upack)
f_sigmavegc =pack(ifn(:,3)/(1.-ifn(:,2)),upack)
f_sigmavegr =pack(ifn(:,4),upack)
f_industryfg=pack(ifn(:,5),upack)
f_trafficfg =pack(ifn(:,6),upack)
f_bldheight =pack(ifn(:,7),upack)
f_roofalpha =pack(ifn(:,8),upack)
f_wallalpha =pack(ifn(:,9),upack)
f_roadalpha =pack(ifn(:,10),upack)
f_vegalphac =pack(ifn(:,11),upack)
f_vegalphac =pack(ifn(:,12),upack)
f_roofemiss =pack(ifn(:,13),upack)
f_wallemiss =pack(ifn(:,14),upack)
f_roademiss =pack(ifn(:,15),upack)
f_vegemissc =pack(ifn(:,16),upack)
f_vegemissr =pack(ifn(:,17),upack)
f_bldtemp   =pack(ifn(:,18),upack)
do ii=1,4
  f_roofdepth(:,ii) =pack(ifn(:,18+ii),upack)
  f_walldepth(:,ii) =pack(ifn(:,22+ii),upack)
  f_roaddepth(:,ii) =pack(ifn(:,26+ii),upack)
  f_roofcp(:,ii)    =pack(ifn(:,30+ii),upack)
  f_wallcp(:,ii)    =pack(ifn(:,34+ii),upack)
  f_roadcp(:,ii)    =pack(ifn(:,38+ii),upack)
  f_rooflambda(:,ii)=pack(ifn(:,42+ii),upack)
  f_walllambda(:,ii)=pack(ifn(:,46+ii),upack)
  f_roadlambda(:,ii)=pack(ifn(:,50+ii),upack)
end do
f_zovegc   =pack(ifn(:,55),upack)
f_vegrlaic =pack(ifn(:,56),upack)
f_vegrsminc=pack(ifn(:,57),upack)
f_zovegr   =pack(ifn(:,58),upack)
f_vegrlair =pack(ifn(:,59),upack)
f_vegrsminr=pack(ifn(:,60),upack)
f_swilt    =pack(ifn(:,61),upack)
f_sfc      =pack(ifn(:,62),upack)
f_ssat     =pack(ifn(:,63),upack)

return
end subroutine atebfndef

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine saves aTEB state arrays (not compulsory)

subroutine atebsave(urban,diag)

implicit none

integer, intent(in) :: diag
integer ii
real, dimension(ifull,28), intent(inout) :: urban

if (ufull==0) return
if (diag>=1) write(6,*) "Save aTEB state arrays"

do ii=1,4
  urban(:,ii)   =unpack(roof%temp(:,ii),upack,urban(:,ii))
  urban(:,ii+4) =unpack(walle%temp(:,ii),upack,urban(:,ii+4))
  urban(:,ii+8) =unpack(wallw%temp(:,ii),upack,urban(:,ii+8))
  urban(:,ii+12)=unpack(road%temp(:,ii),upack,urban(:,ii+12))
end do
urban(:,17)=unpack(road%soilwater(:),upack,urban(:,17))
urban(:,18)=unpack(roof%soilwater(:),upack,urban(:,18))
urban(:,19)=unpack(roof%surfwater(:),upack,urban(:,19))
urban(:,20)=unpack(road%surfwater(:),upack,urban(:,20))
urban(:,21)=unpack(road%leafwater(:),upack,urban(:,21))
urban(:,22)=unpack(roof%leafwater(:),upack,urban(:,22))
urban(:,23)=unpack(roof%snow(:), upack,urban(:,23))
urban(:,24)=unpack(road%snow(:), upack,urban(:,24))
urban(:,25)=unpack(roof%den(:),  upack,urban(:,25))
urban(:,26)=unpack(road%den(:),  upack,urban(:,26))
urban(:,27)=unpack(roof%alpha(:),upack,urban(:,27))
urban(:,28)=unpack(road%alpha(:),upack,urban(:,28))

return
end subroutine atebsave

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! temperature only version of tebsave

subroutine atebsavem(urban,moist,diag)

implicit none

integer, intent(in) :: diag
integer ii
real, dimension(ifull,16), intent(inout) :: urban
real, dimension(ifull,2), intent(inout) :: moist

if (ufull==0) return
if (diag>=1) write(6,*) "Save aTEB state arrays"

do ii=1,4
  urban(:,ii)   =unpack(roof%temp(:,ii), upack,urban(:,ii))
  urban(:,ii+4) =unpack(walle%temp(:,ii),upack,urban(:,ii+4))
  urban(:,ii+8) =unpack(wallw%temp(:,ii),upack,urban(:,ii+8))
  urban(:,ii+12)=unpack(road%temp(:,ii), upack,urban(:,ii+12))
end do
moist(:,1)=unpack(road%soilwater(:),upack,moist(:,1))
moist(:,2)=unpack(roof%soilwater(:),upack,moist(:,2))

return
end subroutine atebsavem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine blends urban momentum and heat roughness lengths
! (This version neglects the displacement height (e.g., for CCAM))
!

subroutine atebzo(zom,zoh,zoq,diag,raw)

implicit none

integer, intent(in) :: diag
real, dimension(ifull), intent(inout) :: zom,zoh,zoq
real, dimension(ufull) :: workb,workc,workd,zmtmp,zhtmp,zqtmp
real, parameter :: zr=1.e-15 ! limits minimum roughness length for heat
logical, intent(in), optional :: raw
logical mode

if (ufull==0) return
if (diag>=1) write(6,*) "Calculate urban roughness lengths"

mode=.false.
if (present(raw)) mode=raw

if (mode) then
  zom=unpack(p_cndzmin*exp(-p_lzom),upack,zom)
  zoh=unpack(p_cndzmin*exp(-p_lzoh),upack,zoh)
  zoq=unpack(p_cndzmin*exp(-p_lzoh),upack,zoq)
else 
  ! evaluate at canyon displacement height (really the atmospheric model should provide a displacement height)
  zmtmp=pack(zom,upack)
  zhtmp=pack(zoh,upack)
  zqtmp=pack(zoq,upack)
  workb=sqrt((1.-sigmau)/log(p_cndzmin/zmtmp)**2+sigmau/p_lzom**2)
  workc=(1.-sigmau)/(log(p_cndzmin/zmtmp)*log(p_cndzmin/zhtmp))+sigmau/(p_lzom*p_lzoh)
  workc=workc/workb
  workd=(1.-sigmau)/(log(p_cndzmin/zmtmp)*log(p_cndzmin/zqtmp))+sigmau/(p_lzom*p_lzoh)
  workd=workd/workb
  workb=p_cndzmin*exp(-1./workb)
  workc=max(p_cndzmin*exp(-1./workc),zr)
  workd=max(p_cndzmin*exp(-1./workd),zr)
  zom=unpack(workb,upack,zom)
  zoh=unpack(workc,upack,zoh)
  zoq=unpack(workd,upack,zoq)
  if (minval(workc)<=zr) write(6,*) "WARN: minimum zoh reached"
end if

return
end subroutine atebzo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine blends the urban drag coeff
!

subroutine atebcd(cduv,cdtq,diag)

implicit none

integer, intent(in) :: diag
real, dimension(ifull), intent(inout) :: cduv,cdtq
real, dimension(ufull) :: ctmp

if (ufull==0) return

ctmp=pack(cduv,upack)
ctmp=(1.-sigmau)*ctmp+sigmau*p_cduv
cduv=unpack(ctmp,upack,cduv)

ctmp=pack(cdtq,upack)
ctmp=(1.-sigmau)*ctmp+sigmau*p_cdtq
cdtq=unpack(ctmp,upack,cdtq)

return
end subroutine atebcd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Store fraction of direct radiation
!

subroutine atebfbeam(is,ifin,fbeam,diag)

implicit none

integer, intent(in) :: is,ifin,diag
integer ifinish,ib,ie,ucount
real, dimension(ifin), intent(in) :: fbeam

if (ufull==0) return

ifinish=is+ifin-1
ucount=count(upack(is:ifinish))
if (ucount==0) return

ib=count(upack(1:is-1))+1
ie=ucount+ib-1
f_fbeam(ib:ie)=pack(fbeam,upack(is:ifinish))

return
end subroutine atebfbeam

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Use Spitter et al (1986) method to estimate fraction of direct
! shortwave radiation (from CABLE v1.4)
!

subroutine atebspitter(is,ifin,fjd,sg,cosin,diag)

implicit none

integer, intent(in) :: is,ifin,diag
integer ib,ie,ucount,ifinish
real, dimension(ifin), intent(in) :: sg,cosin
real, dimension(ufull) :: tmpr,tmpk,tmprat
real, dimension(ufull) :: lsg,lcosin
real, intent(in) :: fjd
real, parameter :: solcon = 1370.

if (ufull==0) return

ifinish=is+ifin-1
ucount=count(upack(is:ifinish))
if (ucount==0) return

ib=count(upack(1:is-1))+1
ie=ucount+ib-1

lsg(ib:ie)   =pack(sg,upack(is:ifinish))
lcosin(ib:ie)=pack(cosin,upack(is:ifinish))

tmpr(ib:ie)=0.847+lcosin(ib:ie)*(1.04*lcosin(ib:ie)-1.61)
tmpk(ib:ie)=(1.47-tmpr(ib:ie))/1.66
where (lcosin(ib:ie)>1.0e-10 .and. lsg(ib:ie)>10.)
  tmprat(ib:ie)=lsg(ib:ie)/(solcon*(1.+0.033*cos(2.*pi*(fjd-10.)/365.))*lcosin(ib:ie))
elsewhere
  tmprat(ib:ie)=0.
end where
where (tmprat(ib:ie)>tmpk(ib:ie))
  f_fbeam(ib:ie)=max(1.-tmpr(ib:ie),0.)
elsewhere (tmprat(ib:ie)>0.35)
  f_fbeam(ib:ie)=min(1.66*tmprat(ib:ie)-0.4728,1.)
elsewhere (tmprat(ib:ie)>0.22)
  f_fbeam(ib:ie)=6.4*(tmprat(ib:ie)-0.22)**2
elsewhere
  f_fbeam(ib:ie)=0.
end where

return
end subroutine atebspitter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates the urban contrabution to albedo.
! (selected grid points only)

! raw   (.false.=blend, .true.=output only)
! split (0=net albedo, 1=direct albedo, 2=diffuse albedo)

subroutine atebalb1(is,ifin,alb,diag,raw,split)

implicit none

integer, intent(in) :: is,ifin,diag
integer ucount,ib,ie,ifinish,albmode
integer, intent(in), optional :: split
real, dimension(ifin), intent(inout) :: alb
real, dimension(ufull) :: ualb,utmp
logical, intent(in), optional :: raw
logical outmode

if (ufull==0) return

outmode=.false.
if (present(raw)) outmode=raw

albmode=0
if (present(split)) albmode=split

ifinish=is+ifin-1
ucount=count(upack(is:ifinish))
if (ucount==0) return

ib=count(upack(1:is-1))+1
ie=ucount+ib-1
call atebalbcalc(ib,ucount,ualb(ib:ie),albmode,diag)

if (outmode) then
  alb(:)=unpack(ualb(ib:ie),upack(is:ifinish),alb)
else
  utmp(ib:ie)=pack(alb,upack(is:ifinish))
  utmp(ib:ie)=(1.-sigmau(ib:ie))*utmp(ib:ie)+sigmau(ib:ie)*ualb(ib:ie)
  alb(:)=unpack(utmp(ib:ie),upack(is:ifinish),alb)
end if

return
end subroutine atebalb1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Albedo calculations

subroutine atebalbcalc(is,ifin,alb,albmode,diag)

implicit none

integer, intent(in) :: is,ifin,diag,albmode
integer ie
real, dimension(ifin), intent(out) :: alb
real, dimension(ifin) :: albu,albr,snowdelta
real, dimension(ifin) :: wallpsi,roadpsi
real, dimension(ifin) :: sg_roof,sg_vegr,sg_road,sg_walle,sg_wallw,sg_vegc,sg_rfsn,sg_rdsn
real, dimension(ifin) :: dumfbeam,effhwratio,effbldheight

ie=ifin+is-1

select case(albmode)
  case default ! net albedo
    dumfbeam=f_fbeam(is:ie)
  case(1)      ! direct albedo
    dumfbeam=1.
  case(2)      ! diffuse albedo
    dumfbeam=0.
end select

effbldheight=max(f_bldheight(is:ie)-6.*f_zovegc(is:ie),0.2)/f_bldheight(is:ie)  ! MJT suggestion for tall vegetation
effhwratio=f_hwratio(is:ie)*effbldheight                                        ! MJT suggestion for tall vegetation
  
! canyon
snowdelta=road%snow(is:ie)/(road%snow(is:ie)+maxrdsn)
call getswcoeff(sg_roof,sg_vegr,sg_road,sg_walle,sg_wallw,sg_vegc,sg_rfsn,sg_rdsn,wallpsi,roadpsi,effhwratio, &
                f_vangle(is:ie),f_hangle(is:ie),dumfbeam,f_sigmavegc(is:ie),f_roadalpha(is:ie),f_vegalphac(is:ie), &
                f_wallalpha(is:ie),road%alpha(is:ie),snowdelta)
sg_walle=sg_walle*effbldheight
sg_wallw=sg_wallw*effbldheight
albu=1.-(f_hwratio(is:ie)*(sg_walle+sg_wallw)*(1.-f_wallalpha(is:ie))+snowdelta*sg_rdsn*(1.-road%alpha(is:ie)) &
    +(1.-snowdelta)*((1.-f_sigmavegc(is:ie))*sg_road*(1.-f_roadalpha(is:ie))+f_sigmavegc(is:ie)*sg_vegc*(1.-f_vegalphac(is:ie))))

! roof
snowdelta=roof%snow(is:ie)/(roof%snow(is:ie)+maxrfsn)
albr=(1.-snowdelta)*((1.-f_sigmavegr(is:ie))*sg_roof*f_roofalpha(is:ie)+f_sigmavegr(is:ie)*sg_vegr*f_vegalphar(is:ie)) &
    +snowdelta*sg_rfsn*roof%alpha(is:ie)

! net
alb=f_sigmabld(is:ie)*albr+(1.-f_sigmabld(is:ie))*albu

return
end subroutine atebalbcalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine stores the zenith angle and the solar azimuth angle
! (single grid point)

subroutine atebnewangle1(is,ifin,cosin,azimuthin,ctimein)

implicit none

integer, intent(in) :: is,ifin
integer ifinish,ucount,ib,ie
real, dimension(ifin), intent(in) :: cosin     ! cosine of zenith angle
real, dimension(ifin), intent(in) :: azimuthin ! azimuthal angle
real, dimension(ifin), intent(in) :: ctimein   ! local hour (0<=ctime<=1)

if (ufull==0) return

ifinish=is+ifin-1
ucount=count(upack(is:ifinish))
if (ucount==0) return

ib=count(upack(1:is-1))+1
ie=ucount+ib-1

f_hangle(ib:ie)=0.5*pi-pack(azimuthin,upack(is:ifinish))
f_vangle(ib:ie)=acos(pack(cosin,upack(is:ifinish)))
f_ctime(ib:ie)=pack(ctimein,upack(is:ifinish))

return
end subroutine atebnewangle1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This version of tebnewangle is for CCAM and TAPM
!

subroutine atebccangle(is,ifin,cosin,rlon,rlat,fjd,slag,dt,sdlt)

implicit none

integer, intent(in) :: is, ifin
integer ifinish, ucount, ib, ie
real, intent(in) :: fjd, slag, dt, sdlt
real cdlt
real, dimension(ifin), intent(in) :: cosin, rlon, rlat
real, dimension(ufull) :: hloc, x, y, lattmp

! cosin = cosine of zenith angle
! rlon = longitude
! rlat = latitude
! fjd = day of year
! slag = sun lag angle
! sdlt = sin declination of sun

if ( ufull == 0 ) return

ifinish = is + ifin - 1
ucount = count(upack(is:ifinish))
if ( ucount == 0 ) return

ib = count(upack(1:is-1)) + 1
ie = ucount + ib - 1

cdlt = sqrt(min(max(1.-sdlt*sdlt,0.),1.))

lattmp(ib:ie) = pack(rlat,upack(is:ifinish))

f_vangle(ib:ie) = acos(pack(cosin,upack(is:ifinish)))

! from CCAM zenith.f
hloc(ib:ie) = 2.*pi*fjd + slag + pi + pack(rlon,upack(is:ifinish)) + dt*pi/86400.
f_ctime(ib:ie) = min(max(mod(0.5*hloc(ib:ie)/pi-0.5,1.),0.),1.)

! estimate azimuth angle
x(ib:ie) = sin(-hloc(ib:ie))*cdlt
y(ib:ie) = -cos(-hloc(ib:ie))*cdlt*sin(lattmp(ib:ie))+cos(lattmp(ib:ie))*sdlt
!azimuth=atan2(x,y)
f_hangle(ib:ie) = 0.5*pi - atan2(x(ib:ie),y(ib:ie))


return
end subroutine atebccangle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calcuates screen level diagnostics
!

subroutine atebscrnout(tscrn,qscrn,uscrn,u10,diag,raw)

implicit none

integer, intent(in) :: diag
real, dimension(ifull), intent(inout) :: tscrn,qscrn,uscrn,u10
real, dimension(ufull) :: tmp
logical, intent(in), optional :: raw
logical mode

if (ufull==0) return

mode=.false.
if (present(raw)) mode=raw

if (mode) then
  tscrn=unpack(p_tscrn,upack,tscrn)
  qscrn=unpack(p_qscrn,upack,qscrn)
  uscrn=unpack(p_uscrn,upack,uscrn)
  u10  =unpack(p_u10,  upack,u10  )
else
  tmp=pack(tscrn,upack)
  tmp=sigmau*p_tscrn+(1.-sigmau)*tmp
  tscrn=unpack(tmp,upack,tscrn)
  tmp=pack(qscrn,upack)
  tmp=sigmau*p_qscrn+(1.-sigmau)*tmp
  qscrn=unpack(tmp,upack,qscrn)
  tmp=pack(uscrn,upack)
  tmp=sigmau*p_uscrn+(1.-sigmau)*tmp
  uscrn=unpack(tmp,upack,uscrn)
  tmp=pack(u10,upack)
  tmp=sigmau*p_u10+(1.-sigmau)*tmp
  u10=unpack(tmp,upack,u10)
end if

return
end subroutine atebscrnout

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Extract urban fraction
subroutine atebsigmau(sigu,diag)

implicit none

integer, intent(in) :: diag
real, dimension(ifull), intent(out) :: sigu

sigu=0.
if (ufull==0) return
sigu=unpack(sigmau,upack,sigu)

return
end subroutine atebsigmau

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

subroutine atebcalc(ofg,oeg,ots,owf,orn,dt,zmin,sg,rg,rnd,snd,rho,temp,mixr,ps,uu,vv,umin,diag,raw)

implicit none

integer, intent(in) :: diag
real, intent(in) :: dt,umin
real, dimension(ifull), intent(in) :: sg,rg,rnd,snd,rho,temp,mixr,ps,uu,vv,zmin
real, dimension(ifull), intent(inout) :: ofg,oeg,ots,owf,orn
real, dimension(ufull) :: tmp
real, dimension(ufull) :: a_sg,a_rg,a_rho,a_temp,a_mixr,a_ps,a_umag,a_udir,a_rnd,a_snd,a_zmin
real, dimension(ufull) :: u_fg,u_eg,u_ts,u_wf,u_rn
logical, intent(in), optional :: raw
logical mode

if (ufull==0) return ! no urban grid points

! mode = .false. implies weight output with urban area cover fraction
! mode = .true. implies no weighting of output with urban area cover fraction (assumes 100% cover)
mode=.false.
if (present(raw)) mode=raw

! Host model meteorological data
a_zmin=pack(zmin,                 upack)
a_sg  =pack(sg,                   upack)
a_rg  =pack(rg,                   upack)
a_rho =pack(rho,                  upack)
a_temp=pack(temp,                 upack)
a_mixr=pack(mixr,                 upack)
a_ps  =pack(ps,                   upack)
a_umag=max(pack(sqrt(uu*uu+vv*vv),upack),umin)
a_udir=pack(atan2(vv,uu),         upack)
a_rnd =pack(rnd-snd,              upack)
a_snd =pack(snd,                  upack)

! Update urban prognostic variables
call atebeval(u_fg,u_eg,u_ts,u_wf,u_rn,dt,a_sg,a_rg,a_rho,a_temp,a_mixr,a_ps,a_umag,a_udir,a_rnd,a_snd,a_zmin,diag)

! export urban fluxes on host grid
if (mode) then
  ofg=unpack(u_fg,upack,ofg)
  oeg=unpack(u_eg,upack,oeg)
  ots=unpack(u_ts,upack,ots)
  owf=unpack(u_wf,upack,owf)
  orn=unpack(u_rn,upack,orn)
else
  tmp=pack(ofg,upack)
  tmp=(1.-sigmau)*tmp+sigmau*u_fg
  ofg=unpack(tmp,upack,ofg)
  tmp=pack(oeg,upack)
  tmp=(1.-sigmau)*tmp+sigmau*u_eg
  oeg=unpack(tmp,upack,oeg)
  tmp=pack(ots,upack)
  tmp=((1.-sigmau)*tmp**4+sigmau*u_ts**4)**0.25
  ots=unpack(tmp,upack,ots)
  tmp=pack(owf,upack)
  tmp=(1.-sigmau)*tmp+sigmau*u_wf
  owf=unpack(tmp,upack,owf)
  tmp=pack(orn,upack)
  tmp=(1.-sigmau)*tmp+sigmau*u_rn
  orn=unpack(tmp,upack,orn)
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

subroutine atebeval(u_fg,u_eg,u_ts,u_wf,u_rn,ddt,a_sg,a_rg,a_rho,a_temp,a_mixr,a_ps,a_umag,a_udir,a_rnd,a_snd,a_zmin,diag)

implicit none

integer, intent(in) :: diag
integer k
real, intent(in) :: ddt
real, dimension(ufull) :: garfsn,gardsn,acflx_roof,acflx_walle,acflx_wallw
real, dimension(ufull) :: rdsntemp,rfsntemp,rdsnmelt,rfsnmelt
real, dimension(ufull) :: wallpsi,roadpsi,fgtop,egtop,qsatr,qsata
real, dimension(ufull) :: cu
real, dimension(ufull) :: we,ww,wr,zolog,a,n,zom,zonet,dis
real, dimension(ufull) :: width,newtemp,roofvegwetfac,roadvegwetfac
real, dimension(ufull) :: z_on_l,pa,dts,dtt,effbldheight,effhwratio
real, dimension(ufull), intent(in) :: a_sg,a_rg,a_rho,a_temp,a_mixr,a_ps,a_umag,a_udir,a_rnd,a_snd,a_zmin
real, dimension(ufull), intent(out) :: u_fg,u_eg,u_ts,u_wf,u_rn
real, dimension(ufull) :: sg_roof,sg_vegr,sg_road,sg_walle,sg_wallw,sg_vegc,sg_rfsn,sg_rdsn
real, dimension(ufull) :: rg_roof,rg_road,rg_walle,rg_wallw,rg_vegc,rg_vegr,rg_rfsn,rg_rdsn
real, dimension(ufull) :: fg_roof,fg_road,fg_walle,fg_wallw,fg_vegc,fg_vegr,fg_rfsn,fg_rdsn
real, dimension(ufull) :: eg_roof,eg_road,eg_vegc,eg_vegr,eg_rfsn,eg_rdsn
real, dimension(ufull) :: acond_roof,acond_road,acond_walle,acond_wallw,acond_vegc,acond_vegr,acond_rfsn,acond_rdsn
real, dimension(ufull) :: condterm_roof,condterm_wall
real, dimension(ufull) :: d_roofdelta,d_roaddelta,d_vegdeltac,d_vegdeltar,d_rfsndelta,d_rdsndelta
real, dimension(ufull) :: d_tempc,d_mixrc,d_tempr,d_mixrr,d_sigd,d_sigr,d_rfdzmin
real, dimension(ufull) :: d_accool,d_canyonrgout,d_roofrgout,d_tranc,d_evapc,d_tranr,d_evapr,d_c1c,d_c1r
real, dimension(ufull) :: d_totdepth,d_netemiss,d_netrad,d_topu
real, dimension(ufull) :: d_cwa,d_cw0,d_cww,d_cwr,d_cra,d_crr,d_crw
real, dimension(ufull) :: d_canyontemp,d_canyonmix,d_acout,d_traf

if ( diag>=1 ) write(6,*) "Evaluating aTEB"

! new snowfall
where ( a_snd>1.e-10 )
  ! update snow density
  roof%den = (roof%snow*roof%den+a_snd*ddt*minsnowden)/(roof%snow+ddt*a_snd)
  road%den = (road%snow*road%den+a_snd*ddt*minsnowden)/(road%snow+ddt*a_snd)
  ! reset snow albedo
  roof%alpha = maxsnowalpha
  road%alpha = maxsnowalpha
end where

! estimate water and snow area cover fractions
d_roofdelta = max(roof%surfwater/maxrfwater,0.)**(2./3.)
d_roaddelta = max(road%surfwater/maxrdwater,0.)**(2./3.)
d_vegdeltac = max(road%leafwater/max(maxvwatf*f_vegrlaic,1.E-8),0.)**(2./3.)
d_vegdeltar = max(roof%leafwater/max(maxvwatf*f_vegrlair,1.E-8),0.)**(2./3.)
d_rfsndelta = roof%snow/(roof%snow+maxrfsn)
d_rdsndelta = road%snow/(road%snow+maxrdsn)

! estimate canyon level air temp and water vapor (displacement height at refheight*building height)
pa      = a_ps*exp(-grav*a_zmin/(rd*a_temp))
d_sigd  = a_ps
d_tempc = a_temp*(d_sigd/pa)**(rd/aircp)
call getqsat(qsatr,d_tempc,d_sigd)
call getqsat(qsata,a_temp,pa)
d_mixrc = a_mixr*qsatr/qsata

! estimate roof level air temperature and water vapor (displacement height at building height)
d_sigr  = a_ps*exp(-grav*f_bldheight*(1.-refheight)/(rd*a_temp))
d_tempr = a_temp*(d_sigr/pa)**(rd/aircp)
call getqsat(qsatr,d_tempr,d_sigr)
d_mixrr = a_mixr*qsatr/qsata

! calculate soil data
d_totdepth = sum(f_roaddepth,2)
call getc1(d_c1c,road%soilwater)
call getc1(d_c1r,roof%soilwater)

! calculate minimum heat pumped into canyon by air conditioning (COP updated in canyonflux)
! (use split form to estimate G_{*,4} flux into room for AC.  newtemp is an estimate of the temperature at tau+1)
! MJT notes - a floating building air temperature will allow us to remove this flux estimate and instead calculate
! the AC flux by calculating the flux required to reduce the building air temperature back to the comfort level 
condterm_roof = 1./(0.5*f_roofdepth(:,4)/f_rooflambda(:,4)+r_si)
condterm_wall = 1./(0.5*f_walldepth(:,4)/f_walllambda(:,4)+r_si)
newtemp  = roof%temp(:,4)-condterm_roof*(roof%temp(:,4)-f_bldtemp)/(f_roofcp(:,4)*f_roofdepth(:,4)/ddt    &
          +condterm_roof)
acflx_roof = (1.-f_sigmavegr)*condterm_roof*(newtemp-f_bldtemp)
newtemp  = walle%temp(:,4)-condterm_wall*(walle%temp(:,4)-f_bldtemp)/(f_wallcp(:,4)*f_walldepth(:,4)/ddt  &
          +condterm_wall)
acflx_walle = condterm_wall*(newtemp-f_bldtemp)
newtemp   = wallw%temp(:,4)-condterm_wall*(wallw%temp(:,4)-f_bldtemp)/(f_wallcp(:,4)*f_walldepth(:,4)/ddt &
           +condterm_wall)
acflx_wallw = condterm_wall*(newtemp-f_bldtemp)
if ( acmeth==1 ) then
  d_acout = max(0.,acflx_roof*f_sigmabld/(1.-f_sigmabld)+f_hwratio*(acflx_walle+acflx_wallw))
  !d_acout = acflx_roof*f_sigmabld/(1.-f_sigmabld)+f_hwratio*(acflx_walle+acflx_wallw) ! test energy conservation
else
  d_acout = 0.
end if

! Here we modify the effective canyon geometry to account for in-canyon vegetation
! MJT notes - move this to initialisation?
effbldheight = max(f_bldheight-6.*f_zovegc,0.2)/f_bldheight    ! MJT suggestion for tall vegetation
effhwratio   = f_hwratio*effbldheight                          ! MJT suggestion for tall vegetation

! calculate shortwave reflections
call getswcoeff(sg_roof,sg_vegr,sg_road,sg_walle,sg_wallw,sg_vegc,sg_rfsn,sg_rdsn,wallpsi,roadpsi,effhwratio,     &
                f_vangle,f_hangle,f_fbeam,f_sigmavegc,f_roadalpha,f_vegalphac,f_wallalpha,road%alpha,d_rdsndelta)
sg_roof  = (1.-f_roofalpha)*sg_roof*a_sg
sg_vegr  = (1.-f_vegalphar)*sg_vegr*a_sg
sg_walle = (1.-f_wallalpha)*sg_walle*a_sg*effbldheight ! shadow due to in-canyon vegetation
sg_wallw = (1.-f_wallalpha)*sg_wallw*a_sg*effbldheight ! shadow due to in-canyon vegetation
sg_road  = (1.-f_roadalpha)*sg_road*a_sg               ! currently no road shadow
sg_vegc  = (1.-f_vegalphac)*sg_vegc*a_sg
sg_rfsn  = (1.-roof%alpha)*sg_rfsn*a_sg
sg_rdsn  = (1.-road%alpha)*sg_rdsn*a_sg

! calculate long wave reflections to nrefl order (pregenerated before canyonflux subroutine)
call getlwcoeff(d_netemiss,d_cwa,d_cra,d_cw0,d_cww,d_crw,d_crr,d_cwr,d_rdsndelta,wallpsi,roadpsi,f_sigmavegc,f_roademiss,  &
                f_vegemissc,f_wallemiss)
p_emiss = d_rfsndelta*snowemiss+(1.-d_rfsndelta)*((1.-f_sigmavegr)*f_roofemiss+f_sigmavegr*f_vegemissr)
p_emiss = f_sigmabld*p_emiss+(1.-f_sigmabld)*(2.*f_wallemiss*effhwratio*d_cwa+d_netemiss*d_cra) ! diagnostic only

! estimate in-canyon surface roughness length
dis   = max(max(max(0.1*effbldheight,zocanyon+0.2),f_zovegc+0.2),zosnow+0.2)
zolog = 1./sqrt(d_rdsndelta/log(dis/zosnow)**2+(1.-d_rdsndelta)*(f_sigmavegc/log(dis/f_zovegc)**2  &
       +(1.-f_sigmavegc)/log(dis/zocanyon)**2))
zonet = dis*exp(-zolog)

! estimate overall bulk urban roughness length
zom = zomratio*f_bldheight
where ( zom*f_sigmabld<zonet*(1.-f_sigmabld) ) ! MJT suggestion
  zom = zonet
end where
n   = road%snow/(road%snow+maxrdsn+0.408*grav*zom)     ! snow cover for urban roughness calc (Douville, et al 1995)
zom = (1.-n)*zom+n*zosnow                              ! blend urban and snow roughness lengths (i.e., snow fills canyon)

! here the first model level is always a_zmin above the displacement height
d_rfdzmin = max(max(abs(a_zmin-f_bldheight*(1.-refheight)),zoroof+0.2),f_zovegr+0.2) ! distance to roof displacement height
p_lzom    = log(a_zmin/zom)
p_cndzmin = a_zmin                                                                   ! distance to canyon displacement height

! calculate canyon wind speed and bulk transfer coefficents
! (i.e., acond = 1/(aerodynamic resistance) )
! some terms are updated when calculating canyon air temperature
select case(resmeth)
  case(0) ! Masson (2000)
    cu=exp(-0.25*effhwratio)
    acond_road =cu ! bulk transfer coefficents are updated in canyonflux
    acond_walle=cu
    acond_wallw=cu
    acond_rdsn =cu
    acond_vegc =cu
  case(1) ! Harman et al (2004)
    we=0. ! for cray compiler
    ww=0. ! for cray compiler
    wr=0. ! for cray compiler
    ! estimate wind speed along canyon surfaces
    call getincanwind(we,ww,wr,a_udir,zonet)
    width=f_bldheight/f_hwratio
    dis=max(0.5*width,zocanyon+0.2)
    zolog=log(dis/zocanyon)
    ! calculate terms for turbulent fluxes
    a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
    acond_walle=a*we                 ! east wall bulk transfer
    acond_wallw=a*ww                 ! west wall bulk transfer
    dis=max(max(max(effbldheight*refheight,zocanyon+0.2),f_zovegc+0.2),zosnow+0.2)
    zolog=log(dis/zocanyon)
    a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
    acond_road=a*wr                  ! road bulk transfer
    zolog=log(dis/f_zovegc)
    a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
    acond_vegc=a*wr
    zolog=log(dis/zosnow)
    a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
    acond_rdsn=a*wr                  ! road snow bulk transfer
  case(2) ! Kusaka et al (2001)
    cu=exp(-0.25*effhwratio)
    acond_road =cu ! bulk transfer coefficents are updated in canyonflux
    acond_walle=cu
    acond_wallw=cu
    acond_rdsn =cu
    acond_vegc =cu
  case(3) ! Harman et al (2004)
    we=0. ! for cray compiler
    ww=0. ! for cray compiler
    wr=0. ! for cray compiler
    call getincanwindb(we,ww,wr,a_udir,zonet)
    width=f_bldheight/f_hwratio
    dis=max(0.5*width,zocanyon+0.2)
    zolog=log(dis/zocanyon)
    a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
    acond_walle=a*we                 ! east wall bulk transfer
    acond_wallw=a*ww                 ! west wall bulk transfer
    dis=max(max(max(effbldheight*refheight,zocanyon+0.2),f_zovegc+0.2),zosnow+0.2)
    zolog=log(dis/zocanyon)
    a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
    acond_road=a*wr                  ! road bulk transfer
    zolog=log(dis/f_zovegc)
    a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
    acond_vegc=a*wr
    zolog=log(dis/zosnow)
    a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
    acond_rdsn=a*wr                  ! road snow bulk transfer
end select

! join two walls into a single wall (testing only)
if ( useonewall==1 ) then
  do k = 1,4
    walle%temp(:,k) = 0.5*(walle%temp(:,k)+wallw%temp(:,k))
    wallw%temp(:,k) = walle%temp(:,k)
  end do
  acond_walle = 0.5*(acond_walle+acond_wallw)
  acond_wallw = acond_walle
  sg_walle    = 0.5*(sg_walle+sg_wallw)
  sg_wallw    = sg_walle
end if

! traffic sensible heat flux
call gettraffic(d_traf,f_ctime,f_trafficfg)
d_traf = d_traf/(1.-f_sigmabld)

! calculate canyon fluxes
call solvecanyon(sg_road,rg_road,fg_road,eg_road,acond_road,                          &
                 sg_walle,rg_walle,fg_walle,acond_walle,                              &
                 sg_wallw,rg_wallw,fg_wallw,acond_wallw,                              &
                 sg_vegc,rg_vegc,fg_vegc,eg_vegc,acond_vegc,                          &
                 sg_rdsn,rg_rdsn,fg_rdsn,eg_rdsn,acond_rdsn,rdsntemp,rdsnmelt,gardsn, &
                 a_umag,a_rho,a_rg,a_rnd,a_snd,                                       &
                 d_canyontemp,d_canyonmix,d_tempc,d_mixrc,d_sigd,d_topu,d_netrad,     &
                 d_roaddelta,d_vegdeltac,d_rdsndelta,d_accool,d_acout,d_traf,         &
                 d_canyonrgout,d_tranc,d_evapc,d_cwa,d_cra,d_cw0,d_cww,d_crw,d_crr,   &
                 d_cwr,d_totdepth,d_c1c,ddt)

! calculate roof fluxes (fg_roof updated in solvetridiag)
eg_roof = 0. ! For cray compiler
call solveroof(sg_rfsn,rg_rfsn,fg_rfsn,eg_rfsn,garfsn,rfsnmelt,rfsntemp,acond_rfsn,d_rfsndelta, &
               sg_vegr,rg_vegr,fg_vegr,eg_vegr,acond_vegr,d_vegdeltar,                          &
               sg_roof,rg_roof,eg_roof,acond_roof,d_roofdelta,                                  &
               a_rg,a_umag,a_rho,a_rnd,a_snd,d_tempr,d_mixrr,d_rfdzmin,d_tranr,d_evapr,d_c1r,   &
               d_sigr,ddt)

! tridiagonal solver coefficents for calculating roof, road and wall temperatures
call solvetridiag(sg_roof,rg_roof,fg_roof,eg_roof,garfsn,acflx_roof,  &
                  sg_walle,rg_walle,fg_walle,acflx_walle,             &
                  sg_wallw,rg_wallw,fg_wallw,acflx_wallw,             &
                  sg_road,rg_road,fg_road,eg_road,gardsn,             &
                  d_rfsndelta,d_rdsndelta,d_tempr,a_rho,acond_roof,   &
                  ddt)

! calculate water/snow budgets for road surface
call updatewater(ddt,road%surfwater,road%soilwater,road%leafwater,road%snow,road%den,road%alpha, &
                 rdsnmelt,a_rnd,a_snd,eg_road,eg_rdsn,d_tranc,d_evapc,d_c1c,d_totdepth,          &
                 f_vegrlaic,wbrelaxc)

! calculate water/snow budgets for roof surface
call updatewater(ddt,roof%surfwater,roof%soilwater,roof%leafwater,roof%snow,roof%den,roof%alpha, &
                 rfsnmelt,a_rnd,a_snd,eg_roof,eg_rfsn,d_tranr,d_evapr,d_c1r,f_vegdepthr,         &
                 f_vegrlair,wbrelaxr)

! calculate runoff (leafwater runoff already accounted for in precip reaching canyon floor)
u_rn = max(roof%surfwater-maxrfwater,0.)*f_sigmabld                                   &
      +max(road%surfwater-maxrdwater,0.)*(1.-d_rdsndelta)*(1.-f_sigmavegc)            &
      +max(roof%snow-maxrfsn,0.)*f_sigmabld                                           &
      +max(road%snow-maxrdsn,0.)*d_rdsndelta                                          &
      +max(road%soilwater-f_ssat,0.)*waterden*d_totdepth*(1.-d_rdsndelta)*f_sigmavegc &
      +max(roof%soilwater-f_ssat,0.)*waterden*f_vegdepthr*f_sigmavegr

! remove round-off problems
road%soilwater(1:ufull) = min(max(road%soilwater(1:ufull),f_swilt),f_ssat)
roof%soilwater(1:ufull) = min(max(roof%soilwater(1:ufull),f_swilt),f_ssat)
roof%surfwater(1:ufull) = min(max(roof%surfwater(1:ufull),0.),maxrfwater)
road%surfwater(1:ufull) = min(max(road%surfwater(1:ufull),0.),maxrdwater)
road%leafwater(1:ufull) = min(max(road%leafwater(1:ufull),0.),maxvwatf*f_vegrlaic)
roof%leafwater(1:ufull) = min(max(roof%leafwater(1:ufull),0.),maxvwatf*f_vegrlair)
roof%snow(1:ufull)      = min(max(roof%snow(1:ufull),0.),maxrfsn)
road%snow(1:ufull)      = min(max(road%snow(1:ufull),0.),maxrdsn)
roof%den(1:ufull)       = min(max(roof%den(1:ufull),minsnowden),maxsnowden)
road%den(1:ufull)       = min(max(road%den(1:ufull),minsnowden),maxsnowden)
roof%alpha(1:ufull)     = min(max(roof%alpha(1:ufull),minsnowalpha),maxsnowalpha)
road%alpha(1:ufull)     = min(max(road%alpha(1:ufull),minsnowalpha),maxsnowalpha)

! combine snow and snow-free tiles for fluxes
d_roofrgout = a_rg - d_rfsndelta*rg_rfsn - (1.-d_rfsndelta)*((1.-f_sigmavegr)*rg_roof+f_sigmavegr*rg_vegr)
fg_roof     = d_rfsndelta*fg_rfsn + (1.-d_rfsndelta)*((1.-f_sigmavegr)*fg_roof+f_sigmavegr*fg_vegr)
eg_roof     = d_rfsndelta*eg_rfsn + (1.-d_rfsndelta)*((1.-f_sigmavegr)*eg_roof+f_sigmavegr*eg_vegr)
fgtop       = d_rdsndelta*fg_rdsn + (1.-d_rdsndelta)*((1.-f_sigmavegc)*fg_road+f_sigmavegc*fg_vegc)   &
             + f_hwratio*(fg_walle+fg_wallw) + d_traf + d_accool
egtop       = d_rdsndelta*eg_rdsn + (1.-d_rdsndelta)*((1.-f_sigmavegc)*eg_road+f_sigmavegc*eg_vegc)

! calculate wetfac for roof and road vegetation (see sflux.f or cable_canopy.f90)
roofvegwetfac = max(min((roof%soilwater-f_swilt)/(f_sfc-f_swilt),1.),0.)
roadvegwetfac = max(min((road%soilwater-f_swilt)/(f_sfc-f_swilt),1.),0.)

! calculate longwave, sensible heat latent heat fluxes
! estimate surface temperature from outgoing longwave radiation
u_ts = ((f_sigmabld*d_roofrgout+(1.-f_sigmabld)*d_canyonrgout)/sbconst)**0.25
u_fg = f_sigmabld*fg_roof + (1.-f_sigmabld)*fgtop + f_industryfg
u_eg = f_sigmabld*eg_roof + (1.-f_sigmabld)*egtop
u_wf = f_sigmabld*(1.-d_rfsndelta)*((1.-f_sigmavegr)*d_roofdelta        &
      + f_sigmavegr*((1.-d_vegdeltar)*roofvegwetfac+d_vegdeltar))       &
      + (1.-f_sigmabld)*(1.-d_rdsndelta)*((1.-f_sigmavegc)*d_roaddelta  &
      + f_sigmavegc*((1.-d_vegdeltac)*roadvegwetfac+d_vegdeltac))

! (re)calculate heat roughness length for MOST (diagnostic only)
call getqsat(a,u_ts,d_sigd)
a   = a*u_wf
dts = u_ts*(1.+0.61*a)
dtt = d_tempc*(1.+0.61*d_mixrc)
select case(zohmeth)
  case(0) ! Use veg formulation
    p_lzoh = 2.3 + p_lzom
    call getinvres(p_cdtq,p_cduv,z_on_l,p_lzoh,p_lzom,p_cndzmin,dts,dtt,a_umag,1)
  case(1) ! Use Kanda parameterisation
    p_lzoh = 2.3 + p_lzom ! replaced in getlna
    call getinvres(p_cdtq,p_cduv,z_on_l,p_lzoh,p_lzom,p_cndzmin,dts,dtt,a_umag,2)
  case(2) ! Use Kanda parameterisation
    p_lzoh = 6. + p_lzom
    call getinvres(p_cdtq,p_cduv,z_on_l,p_lzoh,p_lzom,p_cndzmin,dts,dtt,a_umag,4)
end select

! calculate screen level diagnostics
call scrncalc(a_mixr,a_umag,a_temp,u_ts,d_tempc,d_rdsndelta,d_roaddelta,d_vegdeltac,d_sigd,a,rdsntemp,zonet)

return
end subroutine atebeval

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Tridiagonal solver for temperatures

! This version has an implicit estimate for roof sensible heat flux

! [ ggbX ggcX           ] [ temp ] = [ ggdX ]
! [ ggaX ggbX ggcX      ] [ temp ] = [ ggdX ]
! [      ggaX ggbX ggcX ] [ temp ] = [ ggdX ]
! [           ggaX ggbX ] [ temp ] = [ ggdX ]

subroutine solvetridiag(sg_roof,rg_roof,fg_roof,eg_roof,garfsn,acflx_roof,  &
                        sg_walle,rg_walle,fg_walle,acflx_walle,             &
                        sg_wallw,rg_wallw,fg_wallw,acflx_wallw,             &
                        sg_road,rg_road,fg_road,eg_road,gardsn,             &
                        d_rfsndelta,d_rdsndelta,d_tempr,a_rho,acond_roof,   &
                        ddt)

implicit none

integer k
real, intent(in) :: ddt
real, dimension(ufull,1:4) :: ggaroof,ggawall,ggaroad
real, dimension(ufull,0:4) :: ggbroof,ggbwall,ggbroad
real, dimension(ufull,0:3) :: ggcroof,ggcwall,ggcroad
real, dimension(ufull,0:4) :: ggdroof,ggdwalle,ggdwallw,ggdroad
real, dimension(ufull), intent(out) :: fg_roof
real, dimension(ufull), intent(in) :: sg_roof,rg_roof,eg_roof,garfsn,acflx_roof
real, dimension(ufull), intent(in) :: sg_walle,rg_walle,fg_walle,acflx_walle
real, dimension(ufull), intent(in) :: sg_wallw,rg_wallw,fg_wallw,acflx_wallw
real, dimension(ufull), intent(in) :: sg_road,rg_road,fg_road,eg_road,gardsn
real, dimension(ufull), intent(in) :: d_rfsndelta,d_rdsndelta,d_tempr
real, dimension(ufull), intent(in) :: a_rho
real, dimension(ufull), intent(in) :: acond_roof
real, dimension(ufull) :: n

! Conduction terms
ggaroof(:,1) = -2./(f_roofdepth(:,1)/f_rooflambda(:,1))
ggawall(:,1) = -2./(f_walldepth(:,1)/f_walllambda(:,1))
ggaroad(:,1) = -2./(f_roaddepth(:,1)/f_roadlambda(:,1))
do k = 2,4
  ggaroof(:,k) = -2./(f_roofdepth(:,k-1)/f_rooflambda(:,k-1)+f_roofdepth(:,k)/f_rooflambda(:,k))
  ggawall(:,k) = -2./(f_walldepth(:,k-1)/f_walllambda(:,k-1)+f_walldepth(:,k)/f_walllambda(:,k))
  ggaroad(:,k) = -2./(f_roaddepth(:,k-1)/f_roadlambda(:,k-1)+f_roaddepth(:,k)/f_roadlambda(:,k))
end do
ggbroof(:,0) = 2./(f_roofdepth(:,1)/f_rooflambda(:,1))+(1.-d_rfsndelta)*aircp*a_rho*acond_roof
ggbwall(:,0) = 2./(f_walldepth(:,1)/f_walllambda(:,1))
ggbroad(:,0) = 2./(f_roaddepth(:,1)/f_roadlambda(:,1))
ggbroof(:,1) = 2./(f_roofdepth(:,1)/f_rooflambda(:,1))                                     &
              + 2./(f_roofdepth(:,1)/f_rooflambda(:,1)+f_roofdepth(:,2)/f_rooflambda(:,2)) &
              + f_roofcp(:,1)*f_roofdepth(:,1)/ddt
ggbwall(:,1) = 2./(f_walldepth(:,1)/f_walllambda(:,1))                                     &
              + 2./(f_walldepth(:,1)/f_walllambda(:,1)+f_walldepth(:,2)/f_walllambda(:,2)) &
              + f_wallcp(:,1)*f_walldepth(:,1)/ddt
ggbroad(:,1) = 2./(f_roaddepth(:,1)/f_roadlambda(:,1))                                     &
              + 2./(f_roaddepth(:,1)/f_roadlambda(:,1)+f_roaddepth(:,2)/f_roadlambda(:,2)) &
              + f_roadcp(:,1)*f_roaddepth(:,1)/ddt
do k = 2,3
  ggbroof(:,k) = 2./(f_roofdepth(:,k-1)/f_rooflambda(:,k-1)+f_roofdepth(:,k)/f_rooflambda(:,k))  &
                + 2./(f_roofdepth(:,k)/f_rooflambda(:,k)+f_roofdepth(:,k+1)/f_rooflambda(:,k+1)) &
                + f_roofcp(:,k)*f_roofdepth(:,k)/ddt
  ggbwall(:,k) = 2./(f_walldepth(:,k-1)/f_walllambda(:,k-1)+f_walldepth(:,k)/f_walllambda(:,k))  &
                + 2./(f_walldepth(:,k)/f_walllambda(:,k)+f_walldepth(:,k+1)/f_walllambda(:,k+1)) &
                + f_wallcp(:,k)*f_walldepth(:,k)/ddt
  ggbroad(:,k) = 2./(f_roaddepth(:,k-1)/f_roadlambda(:,k-1)+f_roaddepth(:,k)/f_roadlambda(:,k))  &
                + 2./(f_roaddepth(:,k)/f_roadlambda(:,k)+f_roaddepth(:,k+1)/f_roadlambda(:,k+1)) &
                + f_roadcp(:,k)*f_roaddepth(:,k)/ddt
end do
ggbroof(:,4) = 2./(f_roofdepth(:,3)/f_rooflambda(:,3)+f_roofdepth(:,4)/f_rooflambda(:,4))  &
              + f_roofcp(:,4)*f_roofdepth(:,4)/ddt
ggbwall(:,4) = 2./(f_walldepth(:,3)/f_walllambda(:,3)+f_walldepth(:,4)/f_walllambda(:,4))  &
              + f_wallcp(:,4)*f_walldepth(:,4)/ddt
ggbroad(:,4) = 2./(f_roaddepth(:,3)/f_roadlambda(:,3)+f_roaddepth(:,4)/f_roadlambda(:,4))  &
              + f_roadcp(:,4)*f_roaddepth(:,4)/ddt
ggcroof(:,0) = -2./(f_roofdepth(:,1)/f_rooflambda(:,1))
ggcwall(:,0) = -2./(f_walldepth(:,1)/f_walllambda(:,1))
ggcroad(:,0) = -2./(f_roaddepth(:,1)/f_roadlambda(:,1))
do k = 1,3
  ggcroof(:,k) = -2./(f_roofdepth(:,k)/f_rooflambda(:,k)+f_roofdepth(:,k+1)/f_rooflambda(:,k+1))
  ggcwall(:,k) = -2./(f_walldepth(:,k)/f_walllambda(:,k)+f_walldepth(:,k+1)/f_walllambda(:,k+1))
  ggcroad(:,k) = -2./(f_roaddepth(:,k)/f_roadlambda(:,k)+f_roaddepth(:,k+1)/f_roadlambda(:,k+1))
end do

! remove offset from temperature
do k = 1,4
  roof%temp(:,k) = roof%temp(:,k) - 280.
  walle%temp(:,k) = walle%temp(:,k) - 280.
  wallw%temp(:,k) = wallw%temp(:,k) - 280.
  road%temp(:,k) = road%temp(:,k) - 280.
end do

! surface energy budget, AC and previous temperatures
ggdroof(:,0)  = (1.-d_rfsndelta)*(sg_roof+rg_roof-eg_roof+aircp*a_rho*(d_tempr-280.)*acond_roof)+d_rfsndelta*garfsn
ggdwalle(:,0) = sg_walle+rg_walle-fg_walle
ggdwallw(:,0) = sg_wallw+rg_wallw-fg_wallw
ggdroad(:,0)  = (1.-d_rdsndelta)*(sg_road+rg_road-fg_road-eg_road)+d_rdsndelta*gardsn
do k = 1,3
  ggdroof(:,k)  = roof%temp(:,k)*f_roofcp(:,k)*f_roofdepth(:,k)/ddt
  ggdwalle(:,k) = walle%temp(:,k)*f_wallcp(:,k)*f_walldepth(:,k)/ddt
  ggdwallw(:,k) = wallw%temp(:,k)*f_wallcp(:,k)*f_walldepth(:,k)/ddt
  ggdroad(:,k)  = road%temp(:,k)*f_roadcp(:,k)*f_roaddepth(:,k)/ddt
end do
ggdroof(:,4)  = roof%temp(:,4)*f_roofcp(:,4)*f_roofdepth(:,4)/ddt - acflx_roof   ! acflx_roof is AC flux
ggdwalle(:,4) = walle%temp(:,4)*f_wallcp(:,4)*f_walldepth(:,4)/ddt - acflx_walle ! acflx_walle is AC flux
ggdwallw(:,4) = wallw%temp(:,4)*f_wallcp(:,4)*f_walldepth(:,4)/ddt - acflx_wallw ! acflx_wallw is AC flux
ggdroad(:,4)  = road%temp(:,4)*f_roadcp(:,4)*f_roaddepth(:,4)/ddt
    
! tridiagonal solver (Thomas algorithm) to solve for roof, road and wall temperatures
do k = 1,4
  n = ggaroof(:,k)/ggbroof(:,k-1)
  ggbroof(:,k) = ggbroof(:,k) - n*ggcroof(:,k-1)
  ggdroof(:,k) = ggdroof(:,k) - n*ggdroof(:,k-1)
  n = ggawall(:,k)/ggbwall(:,k-1)
  ggbwall(:,k)  = ggbwall(:,k)  - n*ggcwall(:,k-1)
  ggdwalle(:,k) = ggdwalle(:,k) - n*ggdwalle(:,k-1)
  ggdwallw(:,k) = ggdwallw(:,k) - n*ggdwallw(:,k-1)
  n = ggaroad(:,k)/ggbroad(:,k-1)
  ggbroad(:,k) = ggbroad(:,k) - n*ggcroad(:,k-1)
  ggdroad(:,k) = ggdroad(:,k) - n*ggdroad(:,k-1)
end do
roof%temp(:,4)  = ggdroof(:,4)/ggbroof(:,4)
walle%temp(:,4) = ggdwalle(:,4)/ggbwall(:,4)
wallw%temp(:,4) = ggdwallw(:,4)/ggbwall(:,4)
road%temp(:,4)  = ggdroad(:,4)/ggbroad(:,4)
do k = 3,1,-1
  roof%temp(:,k)  = (ggdroof(:,k) -ggcroof(:,k)*roof%temp(:,k+1) )/ggbroof(:,k)
  walle%temp(:,k) = (ggdwalle(:,k)-ggcwall(:,k)*walle%temp(:,k+1))/ggbwall(:,k)
  wallw%temp(:,k) = (ggdwallw(:,k)-ggcwall(:,k)*wallw%temp(:,k+1))/ggbwall(:,k)
  road%temp(:,k)  = (ggdroad(:,k) -ggcroad(:,k)*road%temp(:,k+1) )/ggbroad(:,k)
end do
p_roofskintemp  = (ggdroof(:,0) -ggcroof(:,0)*roof%temp(:,1) )/ggbroof(:,0)
p_walleskintemp = (ggdwalle(:,0)-ggcwall(:,0)*walle%temp(:,1))/ggbwall(:,0)
p_wallwskintemp = (ggdwallw(:,0)-ggcwall(:,0)*wallw%temp(:,1))/ggbwall(:,0)
p_roadskintemp  = (ggdroad(:,0) -ggcroad(:,0)*road%temp(:,1) )/ggbroad(:,0)

! add offset back to temperature
do k = 1,4
  roof%temp(:,k) = roof%temp(:,k) + 280.
  walle%temp(:,k) = walle%temp(:,k) + 280.
  wallw%temp(:,k) = wallw%temp(:,k) + 280.
  road%temp(:,k) = road%temp(:,k) + 280.
end do
p_roofskintemp = p_roofskintemp + 280.
p_walleskintemp = p_walleskintemp + 280.
p_wallwskintemp = p_wallwskintemp + 280.
p_roadskintemp = p_roadskintemp + 280.

! implicit update for fg_roof to improve stability for thin roof layers
fg_roof = aircp*a_rho*(p_roofskintemp-d_tempr)*acond_roof

return
end subroutine solvetridiag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update water prognostic variables for roads and roofs
                            
subroutine updatewater(ddt,surfwater,soilwater,leafwater,snow,den,alpha, &
                       snmelt,a_rnd,a_snd,eg_surf,eg_snow,d_tran,d_evap, &
                       d_c1,d_totdepth,if_vegrlai,iwbrelax)

implicit none

integer, intent(in) :: iwbrelax
real, intent(in) :: ddt
real, dimension(ufull), intent(inout) :: surfwater,soilwater,leafwater,snow,den,alpha
real, dimension(ufull), intent(in) :: snmelt,a_rnd,a_snd,eg_surf,eg_snow
real, dimension(ufull), intent(in) :: d_tran,d_evap,d_c1,d_totdepth,if_vegrlai
real, dimension(ufull) :: modrnd

modrnd = max( a_rnd-d_evap/lv-max( maxvwatf*if_vegrlai-leafwater, 0. )/ddt, 0. ) ! rainfall reaching the soil under vegetation

! note that since sigmaf=1, then there is no soil evaporation, only transpiration.  Evaporation only occurs from water on leafs.
surfwater = surfwater + ddt*(a_rnd-eg_surf/lv+snmelt)                                         ! surface
soilwater = soilwater + ddt*d_c1*(modrnd+snmelt*den/waterden-d_tran/lv)/(waterden*d_totdepth) ! soil
leafwater = leafwater + ddt*(a_rnd-d_evap/lv)                                                 ! leaf
leafwater = min( max( leafwater, 0. ), maxvwatf*if_vegrlai )

if ( iwbrelax==1 ) then
  ! increase soil moisture for irrigation 
  soilwater = soilwater + max( 0.75*f_swilt+0.25*f_sfc-soilwater, 0. )/(86400./ddt+1.) ! 24h e-fold time
end if

! snow fields
snow  = snow + ddt*(a_snd-eg_snow/lv-snmelt)
den   = den + (maxsnowden-den)/(0.24/(86400.*ddt)+1.)
alpha = alpha + (minsnowalpha-alpha)/(0.24/(86400.*ddt)+1.)

return
end subroutine updatewater

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Estimate saturation mixing ratio

subroutine getqsat_v(qsat,temp,ps)

implicit none

real, dimension(:), intent(out) :: qsat
real, dimension(:), intent(in) :: temp
real, dimension(:), intent(in) :: ps
real, dimension(0:220), save :: table
real, dimension(size(qsat)) :: esatf, tdiff, rx
integer ilen
integer, dimension(size(qsat)) :: ix
logical, save :: first=.true.

if ( first ) then
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
  first = .false.
end if

ilen = size(qsat(:))
tdiff(1:ilen) = min( max( temp(1:ilen)-123.16, 0.), 219. )
rx(1:ilen)    = tdiff(1:ilen) - aint(tdiff(1:ilen))
ix(1:ilen)    = int(tdiff(1:ilen))
esatf(1:ilen) = (1.-rx(1:ilen))*table(ix(1:ilen)) + rx(1:ilen)*table(ix(1:ilen)+1)
qsat(1:ilen)  = 0.622*esatf(1:ilen)/max( ps(1:ilen)-esatf(1:ilen), 0.1 )

return
end subroutine getqsat_v

subroutine getqsat_s(qsat,temp,ps)

implicit none

real, intent(in) :: temp,ps
real, intent(out) :: qsat
real, dimension(3) :: dum

dum(2) = temp
dum(3) = ps
call getqsat_v(dum(1:1),dum(2:2),dum(3:3))
qsat = dum(1)

return
end subroutine getqsat_s

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Interface for calcuating ustar and thetastar

subroutine getinvres_v(invres,cd,z_on_l,olzoh,ilzom,zmin,sthetav,thetav,a_umag,mode)

implicit none

integer, intent(in) :: mode
real, dimension(:), intent(in) :: ilzom
real, dimension(size(ilzom)), intent(in) :: zmin,sthetav,thetav
real, dimension(size(ilzom)), intent(in) :: a_umag
real, dimension(size(ilzom)), intent(out) :: invres,cd,z_on_l
real, dimension(size(ilzom)), intent(inout) :: olzoh
real, dimension(size(ilzom)) :: lna,thetavstar,integralh

lna = olzoh - ilzom
call dyerhicks(integralh,z_on_l,cd,thetavstar,thetav,sthetav,a_umag,zmin,ilzom,lna,mode)
invres = vkar*sqrt(cd)*a_umag/integralh
olzoh = lna + ilzom

return
end subroutine getinvres_v

subroutine getinvres_s(invres,cd,z_on_l,olzoh,ilzom,zmin,sthetav,thetav,a_umag,mode)

integer, intent(in) :: mode
real, intent(in) :: ilzom
real, intent(in) :: zmin,sthetav,thetav
real, intent(in) :: a_umag
real, intent(out) :: invres,cd,z_on_l
real, intent(inout) :: olzoh
real, dimension(9) :: dum

dum(4)=olzoh
dum(5)=ilzom
dum(6)=zmin
dum(7)=sthetav
dum(8)=thetav
dum(9)=a_umag
call getinvres_v(dum(1:1),dum(2:2),dum(3:3),dum(4:4),dum(5:5),dum(6:6),dum(7:7),dum(8:8),dum(9:9),mode)
invres=dum(1)
cd=dum(2)
z_on_l=dum(3)
olzoh=dum(4)

return
end subroutine getinvres_s

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate stability functions using Dyerhicks

subroutine dyerhicks(integralh,z_on_l,cd,thetavstar,thetav,sthetav,umag,zmin,ilzom,lna,mode)

implicit none

integer, intent(in) :: mode
integer ic
real, dimension(:), intent(in) :: thetav
real, dimension(size(thetav)), intent(in) :: sthetav,umag,zmin,ilzom
real, dimension(size(thetav)), intent(inout) :: lna
real, dimension(size(thetav)), intent(out) :: cd,thetavstar
real, dimension(size(thetav)), intent(out) :: integralh,z_on_l
real, dimension(size(thetav)) :: z0_on_l,zt_on_l,olzoh
real, dimension(size(thetav)) :: pm0,ph0,pm1,ph1,integralm
real, parameter :: aa1 = 3.8
real, parameter :: bb1 = 0.5
real, parameter :: cc1 = 0.3

cd=(vkar/ilzom)**2                              ! first guess
call getlna(lna,cd,umag,zmin,ilzom,mode)
olzoh=ilzom+lna
integralh=max(sqrt(cd)*ilzom*olzoh/vkar,1.e-10) ! first guess
thetavstar=vkar*(thetav-sthetav)/integralh      ! first guess

do ic=1,icmax
  z_on_l=vkar*zmin*grav*thetavstar/(thetav*cd*umag**2)
  z_on_l=min(z_on_l,10.)
  z0_on_l  = z_on_l*exp(-ilzom)
  zt_on_l  = z0_on_l*exp(-lna)
  where (z_on_l<0.)
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
  integralm = max(integralm,1.e-10)
  integralh = max(integralh,1.e-10)
  !where (z_on_l<=0.4)
    cd = (max(0.01,min(vkar*umag/integralm,2.))/umag)**2
  !elsewhere
  !  cd = (max(0.01,min(vkar*umag/(aa1*( ( z_on_l**bb1)*(1.0+cc1* z_on_l**(1.-bb1)) &
  !      -(z0_on_l**bb1)*(1.+cc1*z0_on_l**(1.-bb1)) )),2.))/umag)**2
  !endwhere
  thetavstar= vkar*(thetav-sthetav)/integralh
  call getlna(lna,cd,umag,zmin,ilzom,mode)
  olzoh=ilzom+lna
end do

return
end subroutine dyerhicks

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Estimate roughness length for heat
!

subroutine getlna(lna,cd,umag,zmin,ilzom,mode)

implicit none

integer, intent(in) :: mode
real, dimension(:), intent(out) :: lna
real, dimension(size(lna)), intent(in) :: cd,umag,zmin,ilzom
real, dimension(size(lna)) :: re
real, parameter :: nu = 1.461E-5
!real, parameter :: eta0 = 1.827E-5
!real, parameter :: t0 = 291.15
!real, parameter :: c = 120.
!eta=eta0*((t0+c)/(theta+c))*(theta/t0)**(2./3.)
!nu=eta/rho

select case(mode) ! roughness length for heat
  case(1) ! zot=zom/10.
    lna=2.3
  case(2) ! Kanda et al 2007
    re=max(sqrt(cd)*umag*zmin*exp(-ilzom)/nu,10.)
    !lna=2.46*re**0.25-2. !(Brutsaet, 1982)
    lna=1.29*re**0.25-2.  !(Kanda et al, 2007)
  case(3) ! zot=zom (neglect molecular diffusion)
    lna=0.
  case(4) ! user defined
    ! no change
  case DEFAULT
    write(6,*) "ERROR: Unknown getlna mode ",mode
    stop
end select

return
end subroutine getlna

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate shortwave radiation coefficents (modified to include 2nd wall)

subroutine getswcoeff(sg_roof,sg_vegr,sg_road,sg_walle,sg_wallw,sg_vegc,sg_rfsn,sg_rdsn,wallpsi,roadpsi,if_hwratio,if_vangle, &
                      if_hangle,if_fbeam,if_sigmavegc,if_roadalpha,if_vegalphac,if_wallalpha,ird_alpha,rdsndelta)

implicit none

integer k
real, dimension(:), intent(in) :: rdsndelta
real, dimension(size(rdsndelta)), intent(in) :: ird_alpha
real, dimension(size(rdsndelta)), intent(out) :: wallpsi, roadpsi
real, dimension(size(rdsndelta)), intent(in) :: if_hwratio, if_vangle, if_hangle, if_fbeam, if_sigmavegc, if_roadalpha, if_vegalphac
real, dimension(size(rdsndelta)), intent(in) :: if_wallalpha
real, dimension(size(rdsndelta)), intent(out) :: sg_roof, sg_vegr, sg_road, sg_walle, sg_wallw, sg_vegc, sg_rfsn, sg_rdsn
real, dimension(size(rdsndelta)) :: thetazero, walles, wallws, roads, ta, tc, xa, ya, roadnetalpha
real, dimension(size(rdsndelta)) :: nwalles, nwallws, nroads

wallpsi = 0.5*(if_hwratio+1.-sqrt(if_hwratio*if_hwratio+1.))/if_hwratio
roadpsi  =sqrt(if_hwratio*if_hwratio+1.) - if_hwratio

! integrate through 180 deg instead of 360 deg.  Hence paritioning to east and west facing walls
where ( if_vangle >= 0.5*pi )
  walles = 0.
  wallws = 1./if_hwratio
  roads = 0.
elsewhere
  ta = tan(if_vangle)
  thetazero = asin(1./max(if_hwratio*ta,1.))
  tc = 2.*(1.-cos(thetazero))
  xa = min(max(if_hangle-thetazero,0.),pi)-max(if_hangle-pi+thetazero,0.)-min(if_hangle+thetazero,0.)
  ya = cos(max(min(0.,if_hangle),if_hangle-pi))-cos(max(min(thetazero,if_hangle),if_hangle-pi)) &
     + cos(min(0.,-if_hangle))-cos(min(thetazero,-if_hangle)) &
     + cos(max(0.,pi-if_hangle))-cos(max(thetazero,pi-if_hangle))
  ! note that these terms now include the azimuth angle
  walles = if_fbeam*(xa/if_hwratio+ta*ya)/pi + (1.-if_fbeam)*wallpsi
  wallws = if_fbeam*((pi-2.*thetazero-xa)/if_hwratio+ta*(tc-ya))/pi + (1.-if_fbeam)*wallpsi
  roads = if_fbeam*(2.*thetazero-if_hwratio*ta*tc)/pi + (1.-if_fbeam)*roadpsi
end where

! Calculate short wave reflections to nrefl order
roadnetalpha = rdsndelta*ird_alpha+(1.-rdsndelta)*((1.-if_sigmavegc)*if_roadalpha+if_sigmavegc*if_vegalphac)
sg_walle = walles
sg_wallw = wallws
sg_road = roads
do k = 1,nrefl
  nwalles = roadnetalpha*wallpsi*roads + if_wallalpha*(1.-2.*wallpsi)*wallws
  nwallws = roadnetalpha*wallpsi*roads + if_wallalpha*(1.-2.*wallpsi)*walles
  nroads = if_wallalpha*(1.-roadpsi)*0.5*(walles+wallws)
  walles = nwalles
  wallws = nwallws
  roads = nroads
  sg_walle = sg_walle+walles
  sg_wallw = sg_wallw+wallws
  sg_road = sg_road+roads
end do
sg_roof = 1.
sg_vegr = 1.
sg_rfsn = 1.
sg_rdsn = sg_road
sg_vegc = sg_road

return
end subroutine getswcoeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate longwave radiation coefficents (modified to include 2nd wall)

subroutine getlwcoeff(d_netemiss,d_cwa,d_cra,d_cw0,d_cww,d_crw,d_crr,d_cwr,d_rdsndelta,wallpsi,roadpsi,if_sigmavegc, &
                      if_roademiss,if_vegemissc,if_wallemiss)

implicit none

integer k
real, dimension(:), intent(inout) :: d_netemiss
real, dimension(size(d_netemiss)), intent(inout) :: d_cwa,d_cra,d_cw0,d_cww,d_crw,d_crr,d_cwr,d_rdsndelta
real, dimension(size(d_netemiss)), intent(in) :: if_sigmavegc,if_roademiss,if_vegemissc,if_wallemiss
real, dimension(size(d_netemiss)) :: wallpsi,roadpsi
real, dimension(size(d_netemiss)) :: rcwa,rcra,rcw0,rcww,rcrw,rcrr,rcwr
real, dimension(size(d_netemiss)) :: ncwa,ncra,ncwe,ncww,ncrw,ncrr,ncwr


d_netemiss=d_rdsndelta*snowemiss+(1.-d_rdsndelta)*((1.-if_sigmavegc)*if_roademiss+if_sigmavegc*if_vegemissc)
d_cwa=wallpsi
d_cra=roadpsi
d_cw0=0.
d_cww=1.-2.*wallpsi
d_crw=0.5*(1.-roadpsi)
d_crr=0.
d_cwr=wallpsi
rcwa=d_cwa
rcra=d_cra
rcw0=d_cw0
rcww=d_cww
rcrw=d_crw
rcrr=d_crr
rcwr=d_cwr
do k=1,nrefl
  ncwa=(1.-d_netemiss)*wallpsi*rcra+(1.-if_wallemiss)*(1.-2.*wallpsi)*rcwa
  ncra=(1.-if_wallemiss)*(1.-roadpsi)*rcwa
  ncwe=(1.-d_netemiss)*wallpsi*rcrw+(1.-if_wallemiss)*(1.-2.*wallpsi)*rcww
  ncww=(1.-d_netemiss)*wallpsi*rcrw+(1.-if_wallemiss)*(1.-2.*wallpsi)*rcw0
  ncrw=(1.-if_wallemiss)*(1.-roadpsi)*0.5*(rcww+rcw0)  
  ncwr=(1.-d_netemiss)*wallpsi*rcrr+(1.-if_wallemiss)*(1.-2.*wallpsi)*rcwr
  ncrr=(1.-if_wallemiss)*(1.-roadpsi)*rcwr
  rcwa=ncwa
  rcra=ncra
  rcw0=ncwe
  rcww=ncww
  rcrw=ncrw
  rcrr=ncrr
  rcwr=ncwr
  d_cwa=d_cwa+rcwa
  d_cra=d_cra+rcra
  d_cw0=d_cw0+rcw0
  d_cww=d_cww+rcww
  d_crw=d_crw+rcrw
  d_cwr=d_cwr+rcwr
  d_crr=d_crr+rcrr
end do

end subroutine getlwcoeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solve for road snow temperature (includes vegetation canopy temperature and canyon temperature)

subroutine solvecanyon(sg_road,rg_road,fg_road,eg_road,acond_road,                          &
                       sg_walle,rg_walle,fg_walle,acond_walle,                              &
                       sg_wallw,rg_wallw,fg_wallw,acond_wallw,                              &
                       sg_vegc,rg_vegc,fg_vegc,eg_vegc,acond_vegc,                          &
                       sg_rdsn,rg_rdsn,fg_rdsn,eg_rdsn,acond_rdsn,rdsntemp,rdsnmelt,gardsn, &
                       a_umag,a_rho,a_rg,a_rnd,a_snd,                                       &
                       d_canyontemp,d_canyonmix,d_tempc,d_mixrc,d_sigd,d_topu,d_netrad,     &
                       d_roaddelta,d_vegdeltac,d_rdsndelta,d_accool,d_acout,d_traf,         &
                       d_canyonrgout,d_tranc,d_evapc,d_cwa,d_cra,d_cw0,d_cww,d_crw,d_crr,   &
                       d_cwr,d_totdepth,d_c1c,ddt)

implicit none

integer k,l
real, intent(in)    :: ddt
real, dimension(ufull), intent(inout) :: rg_road,fg_road,eg_road,acond_road
real, dimension(ufull), intent(inout) :: rg_walle,fg_walle,acond_walle
real, dimension(ufull), intent(inout) :: rg_wallw,fg_wallw,acond_wallw
real, dimension(ufull), intent(inout) :: rg_vegc,fg_vegc,eg_vegc,acond_vegc
real, dimension(ufull), intent(inout) :: rg_rdsn,fg_rdsn,eg_rdsn,acond_rdsn,rdsntemp,rdsnmelt,gardsn
real, dimension(ufull), intent(in) :: sg_road,sg_walle,sg_wallw,sg_vegc,sg_rdsn
real, dimension(ufull), intent(in) :: a_umag,a_rho,a_rg,a_rnd,a_snd
real, dimension(ufull), intent(inout) :: d_canyontemp,d_canyonmix,d_tempc,d_mixrc,d_sigd,d_topu,d_netrad
real, dimension(ufull), intent(inout) :: d_roaddelta,d_vegdeltac,d_rdsndelta,d_accool,d_acout,d_traf
real, dimension(ufull), intent(inout) :: d_canyonrgout,d_tranc,d_evapc,d_cwa,d_cra,d_cw0,d_cww,d_crw,d_crr,d_cwr
real, dimension(ufull), intent(inout) :: d_totdepth,d_c1c
real, dimension(ufull) :: newval,effbldheight,effhwratio,sndepth,snlambda,ldratio,roadqsat,vegqsat,rdsnqsat
real, dimension(ufull) :: cu,topinvres,dts,dtt,cduv,z_on_l,dumroaddelta,dumvegdelta,res
real, dimension(ufull) :: effwalle,effwallw,effroad,effrdsn,effvegc
real, dimension(ufull) :: aa,bb,cc,dd,ee,ff
real, dimension(ufull,2) :: evct,evctx,oldval

! modify canyon geometry for vegetation shadow
effbldheight = max(f_bldheight-6.*f_zovegc,0.2)/f_bldheight  ! MJT suggestion for tall vegetation
effhwratio   = f_hwratio*effbldheight                        ! MJT suggestion for tall vegetation

! snow conductance
sndepth  = road%snow*waterden/road%den
snlambda = icelambda*(road%den/waterden)**1.88

! first guess for canyon air temperature and water vapor mixing ratio
! also guess for canyon veg and snow temperatures
d_canyontemp    = d_tempc
d_canyonmix     = d_mixrc
p_vegtempc      = d_tempc
rdsntemp        = road%temp(:,1)
rdsnmelt        = 0.
dumvegdelta     = 0. ! cray compiler bug
p_roadskintemp  = road%temp(:,1)
p_walleskintemp = walle%temp(:,1)
p_wallwskintemp = wallw%temp(:,1)
d_netrad=d_rdsndelta*snowemiss*rdsntemp**4+(1.-d_rdsndelta)*((1.-f_sigmavegc)*f_roademiss*p_roadskintemp**4 &
                +f_sigmavegc*f_vegemissc*p_vegtempc**4)

! Solve for canyon air temperature and water vapor mixing ratio
do l = 1,ncyits

  !  solve for aerodynamical resistance between canyon and atmosphere  
  ! assume zoh=zom when coupling to canyon air temperature
  p_lzoh = p_lzom
  dts    = d_canyontemp*(1.+0.61*d_canyonmix)
  dtt    = d_tempc*(1.+0.61*d_mixrc)
  call getinvres(topinvres,cduv,z_on_l,p_lzoh,p_lzom,p_cndzmin,dts,dtt,a_umag,3)
  call gettopu(d_topu,a_umag,z_on_l,f_bldheight,f_hwratio,cduv,p_cndzmin)

  if (resmeth==0) then
    acond_road=(11.8+4.2*sqrt(acond_road**2+cduv*a_umag**2))/(aircp*a_rho) ! From Rowley, et al (1930)
    acond_walle=acond_road
    acond_wallw=acond_road
    acond_rdsn=acond_road
    acond_vegc=acond_road
  else if (resmeth==2) then
    cu=acond_road*d_topu
    where (cu<=5.)
      acond_road=(6.15+4.18*cu)/(aircp*a_rho)
    elsewhere
      acond_road=(7.51*cu**0.78)/(aircp*a_rho)
    end where
    acond_walle=acond_road
    acond_wallw=acond_road
    acond_rdsn=acond_road
    acond_vegc=acond_road
  end if

  ! saturated water vapour mixing ratio for road
  call getqsat(roadqsat,p_roadskintemp,d_sigd)
  
  ! correction for dew
  where (roadqsat<d_canyonmix)
    dumroaddelta=1.
  elsewhere
    dumroaddelta=d_roaddelta
  end where
  
  ! calculate canyon road latent heat flux
  aa=lv*(road%surfwater/ddt+a_rnd+rdsnmelt)
  eg_road=min(lv*a_rho*dumroaddelta*(roadqsat-d_canyonmix)*acond_road*d_topu,aa)
  
  ! calculate road and wall skin temperatures
  ! Write energy budget
  !     Solar_net + Longwave_net - Sensible flux - Latent flux - Conduction = 0
  ! or 
  !     sg + a_rg - f_emiss*sbconst*Tskin**4 - aircp*a_rho*(Tskin-d_tempc) &
  !     -eg - (Tskin-temp(:,1))/ldrratio = 0
  ! as a quartic equation
  !      aa*Tskin^4 + dd*Tskin + ee = 0
  ! and solve for Tskin  
  effwalle=f_wallemiss*(a_rg*d_cwa+sbconst*p_walleskintemp**4*f_wallemiss*d_cw0                   & 
                  +sbconst*p_wallwskintemp**4*f_wallemiss*d_cww+sbconst*d_netrad*d_cwr)
  effwallw=f_wallemiss*(a_rg*d_cwa+sbconst*p_wallwskintemp**4*f_wallemiss*d_cw0                   &
                  +sbconst*p_walleskintemp**4*f_wallemiss*d_cww+sbconst*d_netrad*d_cwr)
  effroad=f_roademiss*(a_rg*d_cra+sbconst*(d_netrad*d_crr-p_roadskintemp**4)                      &
                  +sbconst*f_wallemiss*(p_walleskintemp**4+p_wallwskintemp**4)*d_crw)
  ldratio = 0.5*(f_walldepth(:,1)/f_walllambda(:,1))
  aa = f_wallemiss*sbconst
  dd = aircp*a_rho*acond_walle*d_topu+1./ldratio
  ee = -sg_walle-effwalle-aircp*a_rho*acond_walle*d_topu*d_canyontemp-walle%temp(:,1)/ldratio
  call solvequartic(p_walleskintemp,aa,dd,ee) ! This is an estimate of Tskin to be updated in solvetridiag
  dd = aircp*a_rho*acond_wallw*d_topu+1./ldratio
  ee = -sg_wallw-effwallw-aircp*a_rho*acond_wallw*d_topu*d_canyontemp-wallw%temp(:,1)/ldratio
  call solvequartic(p_wallwskintemp,aa,dd,ee) ! This is an estimate of Tskin to be updated in solvetridiag
  ldratio = 0.5*(f_roaddepth(:,1)/f_roadlambda(:,1))
  aa = f_roademiss*sbconst
  dd = aircp*a_rho*acond_road*d_topu+1./ldratio
  ee = -sg_road-effroad-aircp*a_rho*acond_road*d_topu*d_canyontemp-road%temp(:,1)/ldratio+eg_road
  call solvequartic(p_roadskintemp,aa,dd,ee)  ! This is an estimate of Tskin to be updated in solvetridiag
  
  ! Calculate longwave radiation emitted from the canyon floor
  ! MJT notes - This could be included within the iterative solver for snow and vegetation temperatures.
  ! However, it creates a (weak) coupling between these two variables and therefore could require
  ! a multivariate root finding method (e.g,. Broyden's method). Instead we explicitly solve for d_netrad, 
  ! which allows us to decouple the solutions for snow and vegtation temperatures.
  d_netrad=d_rdsndelta*snowemiss*rdsntemp**4+(1.-d_rdsndelta)*((1.-f_sigmavegc)*f_roademiss*p_roadskintemp**4 &
                  +f_sigmavegc*f_vegemissc*p_vegtempc**4)
  
  ! solve for road snow and canyon veg temperatures -------------------------------
  ldratio  = 0.5*( sndepth/snlambda + f_roaddepth(:,1)/f_roadlambda(:,1) )
  oldval(:,1) = p_vegtempc + 0.5
  oldval(:,2) = rdsntemp + 0.5
  call canyonflux(evct,sg_vegc,rg_vegc,fg_vegc,eg_vegc,acond_vegc,vegqsat,res,dumvegdelta,      &
                  sg_rdsn,rg_rdsn,fg_rdsn,eg_rdsn,acond_rdsn,rdsntemp,gardsn,rdsnmelt,rdsnqsat, &
                  a_rg,a_rho,a_rnd,a_snd,                                                       &
                  d_canyontemp,d_canyonmix,d_sigd,d_topu,d_netrad,d_tranc,d_evapc,              &
                  d_cra,d_crr,d_crw,d_totdepth,d_c1c,d_vegdeltac,d_rdsndelta,                   &
                  effvegc,effrdsn,effbldheight,ldratio,ddt)
  p_vegtempc = p_vegtempc - 0.5
  rdsntemp   = rdsntemp - 0.5
  do k = 1,nfgits ! sectant
    evctx = evct
    call canyonflux(evct,sg_vegc,rg_vegc,fg_vegc,eg_vegc,acond_vegc,vegqsat,res,dumvegdelta,      &
                    sg_rdsn,rg_rdsn,fg_rdsn,eg_rdsn,acond_rdsn,rdsntemp,gardsn,rdsnmelt,rdsnqsat, &
                    a_rg,a_rho,a_rnd,a_snd,                                                       &
                    d_canyontemp,d_canyonmix,d_sigd,d_topu,d_netrad,d_tranc,d_evapc,              &
                    d_cra,d_crr,d_crw,d_totdepth,d_c1c,d_vegdeltac,d_rdsndelta,                   &
                    effvegc,effrdsn,effbldheight,ldratio,ddt)
    evctx = evct-evctx
    where (abs(evctx(:,1))>tol)
      newval      = p_vegtempc-alpha*evct(:,1)*(p_vegtempc-oldval(:,1))/evctx(:,1)
      oldval(:,1) = p_vegtempc
      p_vegtempc  = newval
    end where
    where (abs(evctx(:,2))>tol)
      newval      = rdsntemp-alpha*evct(:,2)*(rdsntemp-oldval(:,2))/evctx(:,2)
      oldval(:,2) = rdsntemp
      rdsntemp    = newval
    end where
  end do
  ! ---------------------------------------------------------------    

  ! balance canyon latent heat budget
  aa = d_rdsndelta*acond_rdsn*d_topu
  bb = (1.-d_rdsndelta)*(1.-f_sigmavegc)*dumroaddelta*acond_road*d_topu
  cc = (1.-d_rdsndelta)*f_sigmavegc*(dumvegdelta*acond_vegc*d_topu+(1.-dumvegdelta)/(1./max(acond_vegc*d_topu,1.e-10)+res))
  dd = topinvres
  d_canyonmix = (aa*rdsnqsat+bb*roadqsat+cc*vegqsat+dd*d_mixrc)/(aa+bb+cc+dd)

  ! update heat pumped into canyon with COP
  d_accool=d_acout*(1.+max(d_canyontemp-f_bldtemp,0.)/f_bldtemp)
  !d_accool=d_acout ! test energy conservation

  ! balance sensible heat flux
  aa = aircp*a_rho*topinvres
  bb = d_rdsndelta*aircp*a_rho*acond_rdsn*d_topu
  cc = (1.-d_rdsndelta)*(1.-f_sigmavegc)*aircp*a_rho*acond_road*d_topu
  dd = (1.-d_rdsndelta)*f_sigmavegc*aircp*a_rho*acond_vegc*d_topu
  ee = effhwratio*aircp*a_rho*acond_walle*d_topu
  ff = effhwratio*aircp*a_rho*acond_wallw*d_topu
  d_canyontemp = (aa*d_tempc+bb*rdsntemp+cc*p_roadskintemp+dd*p_vegtempc+ee*p_walleskintemp+ff*p_wallwskintemp+d_traf+d_accool) &
                /(aa+bb+cc+dd+ee+ff)
end do

! solve for canyon sensible heat flux
fg_walle = aircp*a_rho*(p_walleskintemp-d_canyontemp)*acond_walle*d_topu*effbldheight ! canyon vegetation blocks turblent flux
fg_wallw = aircp*a_rho*(p_wallwskintemp-d_canyontemp)*acond_wallw*d_topu*effbldheight ! canyon vegetation blocks turblent flux
fg_road  = aircp*a_rho*(p_roadskintemp-d_canyontemp)*acond_road*d_topu
fg_vegc  = sg_vegc+rg_vegc-eg_vegc
fg_rdsn  = sg_rdsn+rg_rdsn-eg_rdsn*ls/lv-gardsn

! calculate longwave radiation
effwalle=f_wallemiss*(a_rg*d_cwa+sbconst*p_walleskintemp**4*(f_wallemiss*d_cw0-1.)                      & 
                  +sbconst*p_wallwskintemp**4*f_wallemiss*d_cww+sbconst*d_netrad*d_cwr)
rg_walle=effwalle*effbldheight+sbconst*(d_netrad-f_wallemiss*p_walleskintemp**4)*(1.-effbldheight)
effwallw=f_wallemiss*(a_rg*d_cwa+sbconst*p_wallwskintemp**4*(f_wallemiss*d_cw0-1.)                      &
                  +sbconst*p_walleskintemp**4*f_wallemiss*d_cww+sbconst*d_netrad*d_cwr)
rg_wallw=effwallw*effbldheight+sbconst*(d_netrad-f_wallemiss*p_wallwskintemp**4)*(1.-effbldheight)
effroad=f_roademiss*(a_rg*d_cra+sbconst*(d_netrad*d_crr-p_roadskintemp**4)                              &
                  +sbconst*f_wallemiss*(p_walleskintemp**4+p_wallwskintemp**4)*d_crw)
rg_road=effroad+sbconst*(f_wallemiss*(p_walleskintemp**4+p_wallwskintemp**4)                            &
                  -2.*f_roademiss*p_roadskintemp**4)*f_hwratio*(1.-effbldheight)

! outgoing longwave radiation
! note that eff terms are used for outgoing longwave radiation, whereas rg terms are used for heat conduction
d_canyonrgout=a_rg-d_rdsndelta*effrdsn-(1.-d_rdsndelta)*((1.-f_sigmavegc)*effroad+f_sigmavegc*effvegc) &
                  -effhwratio*(effwalle+effwallw)

return
end subroutine solvecanyon

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solve for canyon veg and snow fluxes
                     
subroutine canyonflux(evct,sg_vegc,rg_vegc,fg_vegc,eg_vegc,acond_vegc,vegqsat,res,dumvegdelta,       &
                      sg_rdsn,rg_rdsn,fg_rdsn,eg_rdsn,acond_rdsn,rdsntemp,gardsn,rdsnmelt,rdsnqsat,  &
                      a_rg,a_rho,a_rnd,a_snd,                                                        &
                      d_canyontemp,d_canyonmix,d_sigd,d_topu,d_netrad,d_tranc,d_evapc,               &
                      d_cra,d_crr,d_crw,d_totdepth,d_c1c,d_vegdeltac,d_rdsndelta,                    &
                      effvegc,effrdsn,effbldheight,ldratio,ddt)

implicit none

real, intent(in) :: ddt
real, dimension(ufull,2), intent(out) :: evct
real, dimension(ufull), intent(inout) :: rg_vegc,fg_vegc,eg_vegc,acond_vegc,vegqsat,res,dumvegdelta
real, dimension(ufull), intent(inout) :: rg_rdsn,fg_rdsn,eg_rdsn,acond_rdsn,rdsntemp,gardsn,rdsnmelt,rdsnqsat
real, dimension(ufull), intent(in) :: sg_vegc,sg_rdsn
real, dimension(ufull), intent(in) :: a_rg,a_rho,a_rnd,a_snd
real, dimension(ufull), intent(in) :: effbldheight,ldratio
real, dimension(ufull), intent(out) :: effvegc,effrdsn
real, dimension(ufull), intent(inout) :: d_canyontemp,d_canyonmix,d_sigd,d_topu,d_netrad,d_tranc,d_evapc
real, dimension(ufull), intent(inout) :: d_cra,d_crr,d_crw,d_totdepth,d_c1c,d_vegdeltac,d_rdsndelta
real, dimension(ufull) :: ff,f1,f2,f3,f4

! estimate mixing ratio for vegetation and snow
call getqsat(vegqsat,p_vegtempc,d_sigd)
call getqsat(rdsnqsat,rdsntemp,d_sigd)

! correction for dew
where (vegqsat<d_canyonmix)
  dumvegdelta=1.
elsewhere
  dumvegdelta=d_vegdeltac
end where
  
! vegetation transpiration terms (developed by Eva in CCAM sflux.f and CSIRO9)
where (f_zovegc<0.5)
  ff=1.1*sg_vegc/max(f_vegrlaic*150.,1.E-8)
elsewhere
  ff=1.1*sg_vegc/max(f_vegrlaic*30.,1.E-8)
end where
f1=(1.+ff)/(ff+f_vegrsminc*f_vegrlaic/5000.)
f2=max(0.5*(f_sfc-f_swilt)/max(road%soilwater-f_swilt,1.E-9),1.)
f3=max(1.-0.00025*(vegqsat-d_canyonmix)*d_sigd/0.622,0.5) ! increased limit from 0.05 to 0.5 following Mk3.6    
f4=max(1.-0.0016*(298.-d_canyontemp)**2,0.05)             ! 0.2 in Mk3.6
res=max(30.,f_vegrsminc*f1*f2/(f3*f4))

! solve for vegetation and snow sensible heat fluxes
fg_vegc=aircp*a_rho*(p_vegtempc-d_canyontemp)*acond_vegc*d_topu
fg_rdsn=aircp*a_rho*(rdsntemp-d_canyontemp)*acond_rdsn*d_topu

! calculate longwave radiation for vegetation and snow
effvegc=f_vegemissc*(a_rg*d_cra+sbconst*(d_netrad*d_crr-p_vegtempc**4)                &
                  +sbconst*f_wallemiss*(p_walleskintemp**4+p_wallwskintemp**4)*d_crw)
rg_vegc=effvegc+sbconst*(f_wallemiss*(p_walleskintemp**4+p_wallwskintemp**4)          &
                  -2.*f_vegemissc*p_vegtempc**4)*f_hwratio*(1.-effbldheight)
effrdsn=snowemiss*(a_rg*d_cra+sbconst*(-rdsntemp**4+d_netrad*d_crr)                   &
                  +sbconst*f_wallemiss*(p_walleskintemp**4+p_wallwskintemp**4)*d_crw)
rg_rdsn=effrdsn+sbconst*(f_wallemiss*(p_walleskintemp**4+p_wallwskintemp**4)          &
                  -2.*snowemiss*rdsntemp**4)*f_hwratio*(1.-effbldheight)

! estimate snow melt
rdsnmelt=d_rdsndelta*max(0.,rdsntemp-273.16)/(icecp*road%den*lf*ddt)

! calculate transpiration and evaporation of in-canyon vegetation
d_tranc=lv*min(max((1.-dumvegdelta)*a_rho*(vegqsat-d_canyonmix)/(1./max(acond_vegc*d_topu,1.e-10)+res),0.), &
               max((road%soilwater-f_swilt)*d_totdepth*waterden/(d_c1c*ddt),0.))
d_evapc=lv*min(dumvegdelta*a_rho*(vegqsat-d_canyonmix)*acond_vegc*d_topu,road%leafwater/ddt+a_rnd)
eg_vegc=d_evapc+d_tranc

! calculate canyon snow latent heat and ground fluxes
eg_rdsn=lv*min(a_rho*d_rdsndelta*max(0.,rdsnqsat-d_canyonmix)*acond_rdsn*d_topu &
             ,road%snow/ddt+a_snd-rdsnmelt)
gardsn=(rdsntemp-p_roadskintemp)/ldratio ! use road temperature to represent canyon bottom surface temperature
                                        ! (i.e., we have ommited soil under vegetation temperature)

! vegetation energy budget error term
evct(:,1) = sg_vegc+rg_vegc-fg_vegc-eg_vegc

! road snow energy balance error term
evct(:,2) = sg_rdsn+rg_rdsn-fg_rdsn-eg_rdsn*ls/lv-gardsn

return
end subroutine canyonflux

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve for roof fluxes

subroutine solveroof(sg_rfsn,rg_rfsn,fg_rfsn,eg_rfsn,garfsn,rfsnmelt,rfsntemp,acond_rfsn,d_rfsndelta, &
                     sg_vegr,rg_vegr,fg_vegr,eg_vegr,acond_vegr,d_vegdeltar,                          &
                     sg_roof,rg_roof,eg_roof,acond_roof,d_roofdelta,                                  &
                     a_rg,a_umag,a_rho,a_rnd,a_snd,d_tempr,d_mixrr,d_rfdzmin,d_tranr,d_evapr,d_c1r,   &
                     d_sigr,ddt)

implicit none

integer k
real, intent(in) :: ddt
real, dimension(ufull), intent(inout) :: rg_rfsn,fg_rfsn,eg_rfsn,garfsn,rfsnmelt,rfsntemp,acond_rfsn
real, dimension(ufull), intent(inout) :: rg_vegr,fg_vegr,eg_vegr,acond_vegr
real, dimension(ufull), intent(inout) :: rg_roof,eg_roof,acond_roof
real, dimension(ufull), intent(in) :: sg_rfsn,sg_vegr,sg_roof
real, dimension(ufull), intent(in) :: a_rg,a_umag,a_rho,a_rnd,a_snd
real, dimension(ufull), intent(inout) :: d_tempr,d_mixrr,d_rfdzmin,d_tranr,d_evapr,d_c1r,d_sigr
real, dimension(ufull), intent(inout) :: d_rfsndelta,d_vegdeltar,d_roofdelta
real, dimension(ufull) :: lzomroof,lzohroof,qsatr,dts,dtt,cdroof,z_on_l,newval,ldratio
real, dimension(ufull) :: aa,dd,ee
real, dimension(ufull,2) :: oldval,evctx,evctveg

p_roofskintemp = roof%temp(:,1) ! 1st estimate for calculating roof snow temp
                                ! replaced below
lzomroof=log(d_rfdzmin/zoroof)
lzohroof=2.3+lzomroof
call getqsat(qsatr,p_roofskintemp,d_sigr)
dts=p_roofskintemp*(1.+0.61*d_roofdelta*qsatr)
dtt=d_tempr*(1.+0.61*d_mixrr)
! Assume zot=0.1*zom (i.e., Kanda et al 2007, small experiment)
call getinvres(acond_roof,cdroof,z_on_l,lzohroof,lzomroof,d_rfdzmin,dts,dtt,a_umag,1)
! acond_roof is multiplied by a_umag


! update green roof and snow temperature
p_vegtempr=d_tempr
rfsntemp  =p_roofskintemp
rg_vegr = f_roofemiss*(a_rg-sbconst*p_roofskintemp**4) ! 1st guess
rg_rfsn = f_roofemiss*(a_rg-sbconst*p_roofskintemp**4) ! 1st guess
eg_vegr = 0.
eg_rfsn = 0.
rfsnmelt = 0.
garfsn = 0.
acond_vegr = acond_roof
acond_rfsn = acond_roof
if ( any( d_rfsndelta>0. .or. f_sigmavegr>0. ) ) then
  evctveg = 0.
  oldval(:,1)=p_vegtempr+0.5
  oldval(:,2)=rfsntemp+0.5
  call roofflux(evctveg,rfsntemp,rfsnmelt,garfsn,sg_vegr,rg_vegr,fg_vegr,eg_vegr,acond_vegr, &
                sg_rfsn,rg_rfsn,fg_rfsn,eg_rfsn,acond_rfsn,a_rg,a_umag,a_rho,a_rnd,a_snd,    &
                d_tempr,d_mixrr,d_rfdzmin,d_tranr,d_evapr,d_c1r,d_sigr,d_vegdeltar,          &
                d_rfsndelta,ddt)
! turn off roof snow and roof vegetation if they are not needed
  where ( f_sigmavegr>0. )
    p_vegtempr=p_vegtempr-0.5
  end where
  where ( d_rfsndelta>0. )
    rfsntemp  =rfsntemp-0.5
  end where
  do k=1,nfgits
    evctx=evctveg
    call roofflux(evctveg,rfsntemp,rfsnmelt,garfsn,sg_vegr,rg_vegr,fg_vegr,eg_vegr,acond_vegr, &
                  sg_rfsn,rg_rfsn,fg_rfsn,eg_rfsn,acond_rfsn,a_rg,a_umag,a_rho,a_rnd,a_snd,    &
                  d_tempr,d_mixrr,d_rfdzmin,d_tranr,d_evapr,d_c1r,d_sigr,d_vegdeltar,          &
                  d_rfsndelta,ddt)
    evctx=evctveg-evctx
    where ( abs(evctx(:,1))>tol .and. f_sigmavegr>0. )
      newval=p_vegtempr-alpha*evctveg(:,1)*(p_vegtempr-oldval(:,1))/evctx(:,1)
      oldval(:,1)=p_vegtempr
      p_vegtempr=newval
    end where
    where ( abs(evctx(:,2))>tol .and. d_rfsndelta>0. )
      newval=rfsntemp-alpha*evctveg(:,2)*(rfsntemp-oldval(:,2))/evctx(:,2)
      oldval(:,2)=rfsntemp
      rfsntemp=newval
    end where
  end do
end if
fg_vegr=sg_vegr+rg_vegr-eg_vegr
fg_rfsn=sg_rfsn+rg_rfsn-eg_rfsn*ls/lv-garfsn


! estimate roof latent heat flux (approx roof_skintemp with roof%temp(:,1))
where (qsatr<d_mixrr)
  ! dew
  eg_roof=lv*a_rho*(qsatr-d_mixrr)*acond_roof
elsewhere
  ! evaporation
  aa=lv*(roof%surfwater/ddt+a_rnd+rfsnmelt)
  eg_roof=min(lv*a_rho*d_roofdelta*(qsatr-d_mixrr)*acond_roof,aa)
end where

! estimate roof skin temperature
! Write roof energy budget
!     Solar_net + Longwave_net - Sensible flux - Latent flux - Conduction = 0
! or 
!     sg_roof + a_rg - f_roofemiss*sbconst*Tskin**4 - aircp*a_rho*(Tskin-d_tempr) &
!     -eg_roof - (Tskin-roof%temp(:,1))/ldrratio = 0
! as a quartic equation
!      aa*Tskin^4 + dd*Tskin + ee = 0
! and solve for Tskin
ldratio=0.5*(f_roofdepth(:,1)/f_rooflambda(:,1))
aa=f_roofemiss*sbconst
dd=aircp*a_rho*acond_roof+1./ldratio
ee=-sg_roof-f_roofemiss*a_rg-aircp*a_rho*acond_roof*d_tempr-roof%temp(:,1)/ldratio+eg_roof
call solvequartic(p_roofskintemp,aa,dd,ee) ! This the 2nd estimate of Tskin to be updated in solvetridiag

! calculate net roof longwave radiation
! (sensible heat flux will be updated in solvetridiag using an implicit method)
rg_roof=f_roofemiss*(a_rg-sbconst*p_roofskintemp**4)
!fg_roof=aircp*a_rho*(p_roofskintemp-d_tempr)*acond_roof ! explicit version


return
end subroutine solveroof

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve for green roof and snow fluxes

subroutine roofflux(evct,rfsntemp,rfsnmelt,garfsn,sg_vegr,rg_vegr,fg_vegr,eg_vegr,acond_vegr,   &
                    sg_rfsn,rg_rfsn,fg_rfsn,eg_rfsn,acond_rfsn,a_rg,a_umag,a_rho,a_rnd,a_snd,   &
                    d_tempr,d_mixrr,d_rfdzmin,d_tranr,d_evapr,d_c1r,d_sigr,d_vegdeltar,         &
                    d_rfsndelta,ddt)

implicit none

real, intent(in) :: ddt
real, dimension(ufull,2), intent(inout) :: evct
real, dimension(ufull), intent(in) :: rfsntemp,sg_vegr,sg_rfsn
real, dimension(ufull), intent(inout) :: rfsnmelt,garfsn
real, dimension(ufull), intent(inout) :: rg_vegr,fg_vegr,eg_vegr,acond_vegr
real, dimension(ufull), intent(inout) :: rg_rfsn,fg_rfsn,eg_rfsn,acond_rfsn
real, dimension(ufull), intent(in) :: a_rg,a_umag,a_rho,a_rnd,a_snd
real, dimension(ufull), intent(inout) :: d_tempr,d_mixrr,d_rfdzmin,d_tranr,d_evapr,d_c1r,d_sigr
real, dimension(ufull), intent(inout) :: d_vegdeltar,d_rfsndelta
real, dimension(ufull) :: lzomvegr,lzohvegr,vwetfac,dts,dtt,z_on_l,ff,f1,f2,f3,f4,cdvegr
real, dimension(ufull) :: vegqsat,dumvegdelta,res,sndepth,snlambda,ldratio,lzosnow,rfsnqsat,cdrfsn
real, dimension(ufull) :: lzotdum

call getqsat(vegqsat,p_vegtempr,d_sigr)
where (vegqsat<d_mixrr)
  dumvegdelta=1.
elsewhere
  dumvegdelta=d_vegdeltar
end where

! transpiration terms (developed by Eva in CCAM sflux.f and CSIRO9)
where (f_zovegr<0.5)
  ff=1.1*sg_vegr/max(f_vegrlair*150.,1.E-8)
elsewhere
  ff=1.1*sg_vegr/max(f_vegrlair*30.,1.E-8)
end where
f1=(1.+ff)/(ff+f_vegrsminr*f_vegrlair/5000.)
f2=max(0.5*(f_sfc-f_swilt)/max(roof%soilwater-f_swilt,1.E-9),1.)
f3=max(1.-.00025*(vegqsat-d_mixrr)*d_sigr/0.622,0.5)
f4=max(1.-0.0016*(298.-d_tempr)**2,0.05)
res=max(30.,f_vegrsminr*f1*f2/(f3*f4))

vwetfac=max(min((roof%soilwater-f_swilt)/(f_sfc-f_swilt),1.),0.) ! veg wetfac (see sflux.f or cable_canopy.f90)
vwetfac=(1.-dumvegdelta)*vwetfac+dumvegdelta
lzomvegr=log(d_rfdzmin/f_zovegr)
lzohvegr=2.3+lzomvegr
dts=p_vegtempr*(1.+0.61*vegqsat*vwetfac)
dtt=d_tempr*(1.+0.61*d_mixrr)
! Assume zot=0.1*zom (i.e., Kanda et al 2007, small experiment)
call getinvres(acond_vegr,cdvegr,z_on_l,lzohvegr,lzomvegr,d_rfdzmin,dts,dtt,a_umag,1)
! acond_vegr is multiplied by a_umag

where ( f_sigmavegr>0. )
  ! longwave radiation
  rg_vegr=f_vegemissr*(a_rg-sbconst*p_vegtempr**4)
  
  ! sensible heat flux
  fg_vegr=aircp*a_rho*(p_vegtempr-d_tempr)*acond_vegr

  ! calculate transpiration and evaporation of in-canyon vegetation
  d_tranr=lv*min(max((1.-dumvegdelta)*a_rho*(vegqsat-d_mixrr)/(1./(acond_vegr*a_umag)+res),0.), &
                 max((roof%soilwater-f_swilt)*f_vegdepthr*waterden/(d_c1r*ddt),0.))
  d_evapr=lv*min(dumvegdelta*a_rho*(vegqsat-d_mixrr)*acond_vegr*a_umag,roof%leafwater/ddt+a_rnd)
  eg_vegr=d_evapr+d_tranr

  ! balance green roof energy budget
  evct(:,1)=sg_vegr+rg_vegr-fg_vegr-eg_vegr
end where


! snow conductance
sndepth=roof%snow*waterden/roof%den
snlambda=icelambda*(roof%den/waterden)**1.88
ldratio=0.5*(sndepth/snlambda+f_roofdepth(:,1)/f_rooflambda(:,1))
  
! Update roof snow energy budget
call getqsat(rfsnqsat,rfsntemp,d_sigr)
lzosnow=log(d_rfdzmin/zosnow)
lzotdum=2.3+lzosnow
dts=rfsntemp*(1.+0.61*rfsnqsat)
!dtt already defined above
call getinvres(acond_rfsn,cdrfsn,z_on_l,lzotdum,lzosnow,d_rfdzmin,dts,dtt,a_umag,1)
! acond_rfsn is multiplied by a_umag

where ( d_rfsndelta>0. )
  rfsnmelt=d_rfsndelta*max(0.,rfsntemp-273.16)/(icecp*roof%den*lf*ddt) 
  rg_rfsn=snowemiss*(a_rg-sbconst*rfsntemp**4)
  fg_rfsn=aircp*a_rho*(rfsntemp-d_tempr)*acond_rfsn
  eg_rfsn=lv*min(a_rho*d_rfsndelta*max(0.,rfsnqsat-d_mixrr)*acond_rfsn,roof%snow/ddt+a_snd-rfsnmelt)
  garfsn=(rfsntemp-p_roofskintemp)/ldratio

  ! balance snow energy budget
  evct(:,2)=sg_rfsn+rg_rfsn-fg_rfsn-eg_rfsn*ls/lv-garfsn
end where

return
end subroutine roofflux

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define traffic flux weights during the diurnal cycle

subroutine gettraffic(trafficout,if_ctime,if_trafficfg)

implicit none

real, dimension(:), intent(out) :: trafficout
real, dimension(size(trafficout)), intent(in) :: if_ctime,if_trafficfg
real, dimension(size(trafficout)) :: rp
integer, dimension(size(trafficout)) :: ip

! traffic diurnal cycle weights approximated from Coutts et al (2007)
real, dimension(25), parameter :: trafficcycle = (/ 0.1, 0.1, 0.1, 0.1, 0.5, 1.0, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, &
                                                    1.5, 1.5, 1.5, 1.5, 1.5, 1.4, 1.2,  1., 0.8, 0.6, 0.4, 0.2, &
                                                    0.1 /) 

ip=int(24.*if_ctime)
rp=24.*if_ctime-real(ip)
where (ip<1) ip=ip+24
trafficout=if_trafficfg*((1.-rp)*trafficcycle(ip)+rp*trafficcycle(ip+1))

return
end subroutine gettraffic


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate in-canyon wind speed for walls and road
! This version allows the eddy size to change with canyon orientation
! which requires a numerical solution to the integral

subroutine getincanwind(ueast,uwest,ufloor,a_udir,z0)

implicit none

real, dimension(ufull), intent(out) :: ueast,uwest,ufloor
real, dimension(ufull), intent(in) :: z0
real, dimension(ufull) :: a,b,wsuma,wsumb,fsum
real, dimension(ufull) :: theta1,wdir,h,w,effbldheight
real, dimension(ufull), intent(in) :: a_udir

! rotate wind direction so that all cases are between 0 and pi
! walls are fliped at the end of the subroutine to account for additional pi rotation
where (a_udir>=0.)
  wdir=a_udir
elsewhere
  wdir=a_udir+pi
endwhere

effbldheight=max(f_bldheight-6.*f_zovegc,0.2)/f_bldheight  ! MJT suggestion for tall vegetation

h=f_bldheight*effbldheight
w=f_bldheight/f_hwratio

theta1=asin(min(w/(3.*h),1.))
wsuma=0.
wsumb=0.
fsum=0.  ! floor

! integrate jet on road, venting side (A)
a=0.
b=max(0.,wdir-pi+theta1)
call integratewind(wsuma,wsumb,fsum,a,b,h,w,wdir,z0,0)

! integrate jet on wall, venting side
a=max(0.,wdir-pi+theta1)
b=max(0.,wdir-theta1)
call integratewind(wsuma,wsumb,fsum,a,b,h,w,wdir,z0,1)

! integrate jet on road, venting side (B)
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

! Correct for rotation of winds at start of subroutine
! 0.5 to adjust for factor of 2 in gettopu
where (a_udir>=0.)
  ueast=0.5*wsuma
  uwest=0.5*wsumb
elsewhere
  ueast=0.5*wsumb
  uwest=0.5*wsuma
end where
ufloor=0.5*fsum      ! floor

return
end subroutine getincanwind

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate in-canyon wind speed for walls and road
! This version fixes the eddy size to the canyon width which allows
! for an analytic solution to the integral

subroutine getincanwindb(ueast,uwest,ufloor,a_udir,z0)

implicit none

real, dimension(ufull), intent(out) :: ueast,uwest,ufloor
real, dimension(ufull), intent(in) :: z0
real, dimension(ufull) :: wsuma,wsumb,fsum
real, dimension(ufull) :: theta1,wdir,h,w
real, dimension(ufull) :: dufa,dura,duva,ntheta
real, dimension(ufull) :: dufb,durb,duvb,effbldheight
real, dimension(ufull), intent(in) :: a_udir

! rotate wind direction so that all cases are between 0 and pi
! walls are fliped at the end of the subroutine to account for additional pi rotation
where (a_udir>=0.)
  wdir=a_udir
elsewhere
  wdir=a_udir+pi
endwhere

effbldheight=max(f_bldheight-6.*f_zovegc,0.2)/f_bldheight  ! MJT suggestion for tall vegetation

h=f_bldheight*effbldheight
w=f_bldheight/f_hwratio

theta1=acos(min(w/(3.*h),1.))

call winda(dufa,dura,duva,h,w,z0) ! jet on road
call windb(dufb,durb,duvb,h,w,z0) ! jet on wall
ntheta=2. ! i.e., int_0^pi sin(theta) dtheta = 2.)
where (wdir<theta1.or.wdir>pi-theta1) ! jet on wall
  wsuma=duvb*ntheta
  wsumb=durb*ntheta
  fsum=dufb*ntheta
elsewhere                                   ! jet on road
  wsuma=dura*ntheta
  wsumb=duva*ntheta
  fsum=dufa*ntheta
end where

! Correct for rotation of winds at start of subroutine
! 0.5 to adjust for factor of 2 in gettopu
where (a_udir>=0.)
  ueast=0.5*wsuma
  uwest=0.5*wsumb
elsewhere
  ueast=0.5*wsumb
  uwest=0.5*wsuma
end where
ufloor=0.5*fsum      ! floor

return
end subroutine getincanwindb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! integrate winds

subroutine integratewind(wsuma,wsumb,fsum,a,b,h,w,wdir,z0,mode)

implicit none

integer, intent(in) :: mode
integer n
!integer, parameter :: ntot=45
integer, parameter :: ntot=1 ! simplified method
real, dimension(ufull), intent(in) :: a,b,h,w,wdir,z0
real, dimension(ufull), intent(inout) :: wsuma,wsumb,fsum
real, dimension(ufull) :: theta,dtheta,st,nw
real, dimension(ufull) :: duf,dur,duv

dtheta=(b-a)/real(ntot)
if (any(dtheta>0.)) then
  select case(mode)
    case(0) ! jet on road
      do n=1,ntot
        theta=dtheta*(real(n)-0.5)+a
        st=abs(sin(theta-wdir))
        nw=max(w/max(st,1.E-9),3.*h)
        call winda(duf,dur,duv,h,nw,z0)
        wsuma=wsuma+dur*st*dtheta
        wsumb=wsumb+duv*st*dtheta
        fsum=fsum+duf*st*dtheta
      end do
    case(1) ! jet on wall
      do n=1,ntot
        theta=dtheta*(real(n)-0.5)+a
        st=abs(sin(theta-wdir))
        nw=min(w/max(st,1.E-9),3.*h)
        call windb(duf,dur,duv,h,nw,z0)
        wsuma=wsuma+dur*st*dtheta
        wsumb=wsumb+duv*st*dtheta
        fsum=fsum+duf*st*dtheta
      end do
    case DEFAULT
      write(6,*) "ERROR: Unknown ateb.f90 integratewind mode ",mode
      stop
  end select
end if

return
end subroutine integratewind

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate canyon wind speeds, jet on road

subroutine winda(uf,ur,uv,h,w,z0)

implicit none

real, dimension(ufull), intent(out) :: uf,ur,uv
real, dimension(ufull), intent(in) :: h,w,z0
real, dimension(ufull) :: a,u0,cuven,zolog

a=0.15*max(1.,3.*h/(2.*w))
u0=exp(-0.9*sqrt(13./4.))

zolog=log(max(h,z0+0.2)/z0)
cuven=log(max(refheight*h,z0+0.2)/z0)/log(max(h,z0+0.2)/z0)
cuven=max(cuven*max(1.-3.*h/w,0.),(u0/a)*(h/w)*(1.-exp(-a*max(w/h-3.,0.))))
uf=(u0/a)*(h/w)*(1.-exp(-3.*a))+cuven
!uf=(u0/a)*(h/w)*(2.-exp(-a*3.)-exp(-a*(w/h-3.)))
ur=(u0/a)*exp(-a*3.)*(1.-exp(-a))
! MJT suggestion
cuven=1.-1./zolog
uv=(u0/a)*exp(-a*max(w/h-3.,0.))*(1.-exp(-a))
uv=max(cuven,uv)
!uv=(u0/a)*exp(-a*(w/h-3.))*(1.-exp(-a))

return
end subroutine winda

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate canyon wind speeds, jet on wall

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

zolog=log(max(h,z0+0.2)/z0)
! MJT suggestion (cuven is multipled by dh to avoid divide by zero)
cuven=h-(h-dh)*log(max(h-dh,z0+0.2)/z0)/zolog-dh/zolog
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
subroutine gettopu(d_topu,a_umag,z_on_l,if_bldheight,if_hwratio,ip_cduv,ip_cndzmin)
      
implicit none

real, dimension(ufull), intent(in) :: z_on_l
real, dimension(ufull) :: z0_on_l,bldheight
real, dimension(ufull) :: pm0,pm1,integralm
real, dimension(ufull) :: ustar,neutral
real, dimension(ufull), intent(inout) :: d_topu
real, dimension(ufull), intent(in) :: a_umag
real, dimension(ufull), intent(in) :: if_bldheight,if_hwratio
real, dimension(ufull), intent(inout) :: ip_cduv,ip_cndzmin

bldheight=if_bldheight*(1.-refheight)
ustar=sqrt(ip_cduv)*a_umag

z0_on_l=min(bldheight,ip_cndzmin)*z_on_l/ip_cndzmin ! calculate at canyon top
z0_on_l=min(z0_on_l,10.)
neutral = log(ip_cndzmin/min(bldheight,ip_cndzmin))
where (z_on_l<0.)
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
end where
where (bldheight<ip_cndzmin)
  d_topu=(2./pi)*(a_umag-ustar*integralm/vkar)
elsewhere ! within canyon
  d_topu=(2./pi)*a_umag*exp(0.5*if_hwratio*(1.-ip_cndzmin/bldheight))
end where
d_topu=max(d_topu,0.1)

return
end subroutine gettopu

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Estimate c1 factor for soil moisture availability

subroutine getc1(dc1,moist)

implicit none

real, dimension(ufull), intent(out) :: dc1
real, dimension(ufull), intent(in) :: moist
!real, dimension(ufull) :: n

!n=min(max(moist/f_ssat,0.218),1.)
!dc1=(1.78*n+0.253)/(2.96*n-0.581)

dc1=1.478 ! simplify water conservation

return
end subroutine getc1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate screen diagnostics

subroutine scrncalc(a_mixr,a_umag,a_temp,u_ts,d_tempc,d_rdsndelta,d_roaddelta,d_vegdeltac,d_sigd,smixr,rdsntemp,zonet)
      
implicit none

real, dimension(ufull), intent(in) :: smixr,rdsntemp,zonet
real, dimension(ufull) :: cd,thetav,sthetav
real, dimension(ufull) :: thetavstar,z_on_l,z0_on_l
real, dimension(ufull) :: pm0,ph0,pm1,ph1,integralm,integralh
real, dimension(ufull) :: ustar,qstar,z10_on_l
real, dimension(ufull) :: neutral,neutral10,pm10
real, dimension(ufull) :: integralm10,tts,tetp
real, dimension(ufull) :: tstar,lna
real, dimension(ufull) :: utop,ttop,qtop,wf,tsurf,qsurf,n
real, dimension(ufull), intent(in) :: a_mixr,a_umag,a_temp
real, dimension(ufull), intent(in) :: u_ts
real, dimension(ufull), intent(in) :: d_tempc,d_rdsndelta,d_roaddelta,d_vegdeltac,d_sigd
real, parameter :: z0  = 1.5
real, parameter :: z10 = 10.

select case(scrnmeth)
  case(0) ! estimate screen diagnostics (slab at displacement height approach)
    thetav=d_tempc*(1.+0.61*a_mixr)
    sthetav=u_ts*(1.+0.61*smixr)
    lna=p_lzoh-p_lzom
    call dyerhicks(integralh,z_on_l,cd,thetavstar,thetav,sthetav,a_umag,p_cndzmin,p_lzom,lna,4)
    ustar=sqrt(cd)*a_umag
    qstar=vkar*(a_mixr-smixr)/integralh  
    tstar=vkar*(a_temp-u_ts)/integralh
    
    z0_on_l  = z0*z_on_l/p_cndzmin
    z10_on_l = z10*z_on_l/p_cndzmin
    z0_on_l  = min(z0_on_l,10.)
    z10_on_l = min(z10_on_l,10.)
    neutral   = log(p_cndzmin/z0)
    neutral10 = log(p_cndzmin/z10)
    where (z_on_l<0.)
      ph0     = (1.-16.*z0_on_l)**(-0.50)
      ph1     = (1.-16.*z_on_l)**(-0.50)
      pm0     = (1.-16.*z0_on_l)**(-0.25)
      pm10    = (1.-16.*z10_on_l)**(-0.25)
      pm1     = (1.-16.*z_on_l)**(-0.25)
      integralh   = neutral-2.*log((1.+1./ph1)/(1.+1./ph0))
      integralm   = neutral-2.*log((1.+1./pm1)/(1.+1./pm0)) &
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
      integralh   = neutral-(ph1-ph0)
      integralm   = neutral-(pm1-pm0)
      integralm10 = neutral10-(pm1-pm10)
    endwhere
    p_tscrn = a_temp-tstar*integralh/vkar
    p_qscrn = a_mixr-qstar*integralh/vkar
    p_uscrn = max(a_umag-ustar*integralm/vkar,0.)
    p_u10   = max(a_umag-ustar*integralm10/vkar,0.)
    
  case(1) ! estimate screen diagnostics (two step canopy approach)
    thetav=d_tempc*(1.+0.61*a_mixr)
    sthetav=u_ts*(1.+0.61*smixr)
    lna=p_lzoh-p_lzom
    call dyerhicks(integralh,z_on_l,cd,thetavstar,thetav,sthetav,a_umag,p_cndzmin,p_lzom,lna,4)
    ustar=sqrt(cd)*a_umag
    qstar=vkar*(a_mixr-smixr)/integralh
    tts=vkar*(thetav-sthetav)/integralh
    tstar=vkar*(a_temp-u_ts)/integralh
    
    z0_on_l  = f_bldheight*(1.-refheight)*z_on_l/p_cndzmin ! calculate at canyon top
    z10_on_l = max(z10-f_bldheight*refheight,1.)*z_on_l/p_cndzmin
    z0_on_l  = min(z0_on_l,10.)
    z10_on_l = min(z10_on_l,10.)
    neutral   = log(p_cndzmin/(f_bldheight*(1.-refheight)))
    neutral10 = log(p_cndzmin/max(z10-f_bldheight*refheight,1.))
    where (z_on_l<0.)
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
      integralh   = neutral-(ph1-ph0)
      integralm   = neutral-(pm1-pm0)
      integralm10 = neutral10-(pm1-pm10)
    endwhere
    ttop = thetav-tts*integralh/vkar
    tetp = a_temp-tstar*integralh/vkar
    qtop = a_mixr-qstar*integralh/vkar
    utop = a_umag-ustar*integralm/vkar

    where (f_bldheight<=z10) ! above canyon
      p_u10=max(a_umag-ustar*integralm10/vkar,0.)
    end where

    ! assume standard stability functions hold for urban canyon (needs more work)
    tsurf=d_rdsndelta*rdsntemp+(1.-d_rdsndelta)*((1.-f_sigmavegc)*p_roadskintemp+f_sigmavegc*p_vegtempc)
    n=max(min((road%soilwater-f_swilt)/(f_sfc-f_swilt),1.),0.)
    wf=(1.-d_rdsndelta)*((1.-f_sigmavegc)*d_roaddelta+f_sigmavegc*((1.-d_vegdeltac)*n+d_vegdeltac))
    call getqsat(qsurf,tsurf,d_sigd)
    qsurf=qsurf*wf
    n=log(f_bldheight/zonet)
    
    thetav=ttop*(1.+0.61*qtop)
    sthetav=tsurf*(1.+0.61*qsurf)
    lna=2.3
    call dyerhicks(integralh,z_on_l,cd,thetavstar,thetav,sthetav,utop,f_bldheight,n,lna,1)
    ustar=sqrt(cd)*utop
    tstar=vkar*(tetp-tsurf)/integralh
    qstar=vkar*(qtop-qsurf)/integralh
    
    z0_on_l   = z0*z_on_l/f_bldheight
    z10_on_l  = max(z10,f_bldheight)*z_on_l/f_bldheight
    z0_on_l   = min(z0_on_l,10.)
    z10_on_l  = min(z10_on_l,10.)
    neutral   = log(f_bldheight/z0)
    neutral10 = log(f_bldheight/max(z10,f_bldheight))
    where (z_on_l<0.)
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
      integralh   = neutral-(ph1-ph0)
      integralm   = neutral-(pm1-pm0)
      integralm10 = neutral10-(pm1-pm10)
    endwhere

    p_tscrn = tetp-tstar*integralh/vkar
    p_qscrn = qtop-qstar*integralh/vkar
    p_uscrn = max(utop-ustar*integralm/vkar,0.)
    where (f_bldheight>z10) ! within canyon
      p_u10 = max(utop-ustar*integralm10/vkar,0.)
    end where

  case(2) ! calculate screen diagnostics from canyon only
    tsurf=d_rdsndelta*rdsntemp+(1.-d_rdsndelta)*((1.-f_sigmavegc)*p_roadskintemp+f_sigmavegc*p_vegtempc)
    n=max(min((road%soilwater-f_swilt)/(f_sfc-f_swilt),1.),0.)
    wf=(1.-d_rdsndelta)*((1.-f_sigmavegc)*d_roaddelta+f_sigmavegc*((1.-d_vegdeltac)*n+d_vegdeltac))
    call getqsat(qsurf,tsurf,d_sigd)
    qsurf=qsurf*wf
    n=log(f_bldheight/zonet)

    thetav=d_tempc*(1.+0.61*a_mixr)
    sthetav=tsurf*(1.+0.61*qsurf)
    lna=2.3
    call dyerhicks(integralh,z_on_l,cd,thetavstar,thetav,sthetav,a_umag,p_cndzmin,n,lna,1)
    ustar=sqrt(cd)*a_umag
    qstar=vkar*(a_mixr-smixr)/integralh
    tstar=vkar*(a_temp-tsurf)/integralh
    
    z0_on_l  = z0*z_on_l/p_cndzmin
    z10_on_l = z10*z_on_l/p_cndzmin
    z0_on_l  = min(z0_on_l,10.)
    z10_on_l = min(z10_on_l,10.)
    neutral   = log(p_cndzmin/z0)
    neutral10 = log(p_cndzmin/z10)
    where (z_on_l<0.)
      ph0     = (1.-16.*z0_on_l)**(-0.50)
      ph1     = (1.-16.*z_on_l)**(-0.50)
      pm0     = (1.-16.*z0_on_l)**(-0.25)
      pm10    = (1.-16.*z10_on_l)**(-0.25)
      pm1     = (1.-16.*z_on_l)**(-0.25)
      integralh   = neutral-2.*log((1.+1./ph1)/(1.+1./ph0))
      integralm   = neutral-2.*log((1.+1./pm1)/(1.+1./pm0)) &
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
      integralh   = neutral-(ph1-ph0)
      integralm   = neutral-(pm1-pm0)
      integralm10 = neutral10-(pm1-pm10)
    endwhere
    p_tscrn = a_temp-tstar*integralh/vkar
    p_qscrn = a_mixr-qstar*integralh/vkar
    p_uscrn = max(a_umag-ustar*integralm/vkar,0.)
    p_u10   = max(a_umag-ustar*integralm10/vkar,0.)
    
end select
p_qscrn       = max(p_qscrn,1.E-4)
      
return
end subroutine scrncalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solves Quartic equation of the form a*x^4+d*x+e=0

subroutine solvequartic(x,a,d,e)

implicit none

real, dimension(:), intent(out) :: x
real, dimension(size(x)), intent(in) :: a,d,e
real, dimension(size(x)) :: t1,q,s,qq,d0,d1

d0 = 12.*a*e
d1 = 27.*a*d*d
qq = (0.5*(d1+sqrt(d1**2-4.*d0**3)))**(1./3.)
s = 0.5*sqrt((qq+d0/qq)/(3.*a))
q = d/a
t1 = -s+0.5*sqrt(-4*s**2+q/s) ! only need first solution
!t2 = -s-0.5*sqrt(-4*s**2+q/s)
!t3 = s+0.5*sqrt(-4*s**2-q/s)
!t4 = s-0.5*sqrt(-4*s**2-q/s)

x=t1

return
end subroutine solvequartic

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Disables aTEB so subroutine calls have no effect

subroutine atebdisable(diag)

implicit none

integer, intent(in) :: diag

if (diag>=1) write(6,*) "Disable aTEB"
ufull=0

return
end subroutine atebdisable

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module ateb
