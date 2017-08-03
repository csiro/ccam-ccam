! aTEB urban canopy model
    
! Copyright 2015-2017 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
! This file is part of the aTEB urban canopy model
!
! aTEB is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! aTEB is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with aTEB.  If not, see <http://www.gnu.org/licenses/>.

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
       atebnewangle1,atebccangle,atebdisable,atebloadm,atebsavem,atebcd,               &
       atebdwn,atebscrnout,atebfbeam,atebspitter,atebsigmau,energyrecord,atebdeftype,  &
       atebhydro,atebenergy,atebloadd,atebsaved,atebcalc_thread
public atebnmlfile,urbtemp,energytol,resmeth,useonewall,zohmeth,acmeth,nrefl,vegmode,  &
       soilunder,conductmeth,scrnmeth,wbrelaxc,wbrelaxr,lweff,ncyits,nfgits,tol,alpha, &
       zosnow,snowemiss,maxsnowalpha,minsnowalpha,maxsnowden,minsnowden,refheight,     &
       zomratio,zocanyon,zoroof,maxrfwater,maxrdwater,maxrfsn,maxrdsn,maxvwatf,r_si,   &
       intairtmeth,intmassmeth,ac_cap

#ifdef CCAM
public sigmau_g,upack_g,ufull_g,nl
public f_industryfg,f_bldheight,f_bldwidth,f_coeffbldheight,f_ctime,f_effhwratio
public f_fbeam,f_hangle,f_hwratio,f_intgains_flr,f_intm,f_intmassn,f_rfvegdepth,f_road
public f_roof,f_sfc,f_sigmabld,f_slab,f_ssat,f_swilt,f_trafficfg,f_vangle,f_wall,f_ach
public f_tempcool,f_tempheat,f_bldairtemp
public p_bldheat,p_bldcool,p_traf,p_intgains_full,p_snowmelt,p_cndzmin,p_lzom,p_lzoh
public p_cdtq,p_cduv,p_atmoserr,p_surferr,p_qscrn,p_tscrn,p_u10,p_uscrn
public facetparams,facetdata,hydrodata,vegdata,intm,p_emiss,rdhyd,rfhyd,rfveg
public road,roof,room,slab,walle,wallw,cnveg,int_psi,int_viewf
#endif

! state arrays
integer, save :: ufull_g, ifull, iqut
logical, dimension(:), allocatable, save :: upack_g
real, dimension(:), allocatable, save :: sigmau_g
real, dimension(:,:), allocatable, save :: atebdwn ! These variables are for CCAM onthefly.f
real, dimension(:), allocatable, save :: f_hwratio,f_coeffbldheight,f_effhwratio,f_sigmabld
real, dimension(:), allocatable, save :: f_industryfg,f_intgains_flr,f_trafficfg,f_bldheight,f_bldwidth
real, dimension(:), allocatable, save :: f_ctime,f_vangle,f_hangle,f_fbeam
real, dimension(:), allocatable, save :: f_bldairtemp,p_bldheat,p_bldcool,p_traf,p_intgains_full
real, dimension(:), allocatable, save :: f_swilt,f_sfc,f_ssat,f_rfvegdepth
real, dimension(:), allocatable, save :: f_ach,f_tempheat,f_tempcool
real, dimension(:), allocatable, save :: p_lzom,p_lzoh,p_cndzmin,p_cduv,p_cdtq
real, dimension(:), allocatable, save :: p_tscrn,p_qscrn,p_uscrn,p_u10,p_emiss,p_snowmelt
real, dimension(0:220), save :: table
real(kind=8), dimension(:), allocatable, save :: p_surferr,p_atmoserr,p_surferr_bias,p_atmoserr_bias
real(kind=8), dimension(:), allocatable, save :: p_storagetot_net
real(kind=8), dimension(:,:), allocatable, save :: p_storagetot_road, p_storagetot_walle, p_storagetot_wallw, p_storagetot_roof
real(kind=8), save, allocatable, dimension(:,:,:) :: int_psi, int_viewf         ! internal radiation
integer, dimension(:), allocatable, save :: f_intmassn

type facetdata
  real, dimension(:,:), allocatable :: nodetemp        ! Temperature of node (prognostic)         [K]
  real(kind=8), dimension(:,:), allocatable :: storage ! Facet energy storage (diagnostic)
end type facetdata

type facetparams
  real, dimension(:,:), allocatable :: depth         ! Layer depth                              [m]
  real, dimension(:,:), allocatable :: volcp         ! Layer volumetric heat capacity           [J m^-3 K-1]
  real, dimension(:,:), allocatable :: lambda        ! Layer conductivity                       [W m^-1 K^-1]
  real, dimension(:),   allocatable :: alpha         ! Facet albedo (internal & external)
  real, dimension(:),   allocatable :: emiss         ! Facet emissivity (internal & external)
end type facetparams

type hydrodata
  real, dimension(:), allocatable   :: surfwater
  real, dimension(:), allocatable   :: leafwater
  real, dimension(:), allocatable   :: soilwater
  real, dimension(:), allocatable   :: snow
  real, dimension(:), allocatable   :: den
  real, dimension(:), allocatable   :: snowalpha
end type hydrodata

type vegdata
  real, dimension(:), allocatable :: temp          ! Temperature of veg (prognostic)  [K]
  real, dimension(:), allocatable :: sigma         ! Fraction of veg on roof/canyon
  real, dimension(:), allocatable :: alpha         ! Albedo of veg
  real, dimension(:), allocatable :: emiss         ! Emissivity of veg
  real, dimension(:), allocatable :: lai           ! Leaf area index of veg
  real, dimension(:), allocatable :: zo            ! Roughness of veg
  real, dimension(:), allocatable :: rsmin         ! Minimum stomatal resistance of veg
end type vegdata

type(facetdata),   save :: roof, road, walle, wallw, slab, intm, room
type(facetparams), save :: f_roof, f_road, f_wall, f_slab, f_intm
type(hydrodata),   save :: rfhyd, rdhyd
type(vegdata),     save :: cnveg, rfveg


! model parameters
integer, save      :: atebnmlfile=11       ! Read configuration from nml file (0=off, >0 unit number (default=11))
integer, save      :: resmeth=1            ! Canyon sensible heat transfer (0=Masson, 1=Harman (varying width), 2=Kusaka,
                                           ! 3=Harman (fixed width))
integer, save      :: useonewall=0         ! Combine both wall energy budgets into a single wall (0=two walls, 1=single wall) 
integer, save      :: zohmeth=1            ! Urban roughness length for heat (0=0.1*zom, 1=Kanda, 2=0.003*zom)
integer, save      :: acmeth=1             ! AC heat pump into canyon (0=Off, 1=On, 2=Reversible, COP of 1.0)
integer, save      :: intairtmeth=1        ! Internal air temperature (0=fixed, 1=implicit varying)
integer, save      :: intmassmeth=2        ! Internal thermal mass (0=none, 1=one floor, 2=dynamic floor number)
integer, save      :: nrefl=3              ! Number of canyon reflections for radiation (default=3)
integer, save      :: vegmode=2            ! In-canyon vegetation mode (0=50%/50%, 1=100%/0%, 2=0%/100%, where out/in=X/Y.
                                           ! Negative values are X=abs(vegmode))
integer, save      :: soilunder=1          ! Modify road heat capacity to extend under
                                           ! (0=road only, 1=canveg, 2=bld, 3=canveg & bld)
integer, save      :: conductmeth=1        ! Conduction method (0=half-layer, 1=interface)
integer, save      :: scrnmeth=1           ! Screen diagnostic method (0=Slab, 1=Hybrid, 2=Canyon)
integer, save      :: wbrelaxc=0           ! Relax canyon soil moisture for irrigation (0=Off, 1=On)
integer, save      :: wbrelaxr=0           ! Relax roof soil moisture for irrigation (0=Off, 1=On)
integer, save      :: lweff=2              ! Modification of LW flux for effective canyon height (0=insulated, 1=coupled, 2=full)
integer, parameter :: nl=4                 ! Number of layers (default 4, must be factors of 4)
integer, save      :: iqt=314              ! Diagnostic point (in terms of host grid)
real, save         :: ac_cap=6.            ! capacity of ac in W/m^3
#ifndef CCAM
integer, parameter :: ntiles=1             ! Emulate OMP
#endif
! sectant solver parameters
integer, save      :: ncyits=6             ! Number of iterations for balancing canyon sensible and latent heat fluxes (default=6)
integer, save      :: nfgits=3             ! Number of iterations for balancing veg and snow energy budgets (default=3)
real, save         :: tol=0.001            ! Sectant method tolarance for sensible heat flux (default=0.001)
real, save         :: alpha=1.             ! Weighting for determining the rate of convergence when calculating canyon temperatures
real(kind=8), save :: energytol=0.005_8    ! Tolerance for acceptable energy closure in each timestep
real, save         :: urbtemp=290.         ! reference temperature to improve precision
! physical parameters
real, parameter    :: waterden=1000.       ! water density (kg m^-3)
real, parameter    :: icelambda=2.22       ! conductance of ice (W m^-1 K^-1)
real, parameter    :: aircp=1004.64        ! Heat capapcity of dry air (J kg^-1 K^-1)
real, parameter    :: icecp=2100.          ! Heat capacity of ice (J kg^-1 K^-1)
real, parameter    :: grav=9.80616         ! gravity (m s^-2)
real, parameter    :: vkar=0.4             ! von Karman constant
real, parameter    :: lv=2.501e6           ! Latent heat of vaporisation (J kg^-1)
real, parameter    :: lf=3.337e5           ! Latent heat of fusion (J kg^-1)
real, parameter    :: ls=lv+lf             ! Latent heat of sublimation (J kg^-1)
real, parameter    :: pi=3.14159265        ! pi (must be rounded down for shortwave)
real, parameter    :: rd=287.04            ! Gas constant for dry air
real, parameter    :: rv=461.5             ! Gas constant for water vapor
real, parameter    :: sbconst=5.67e-8      ! Stefan-Boltzmann constant
! snow parameters
real, save         :: zosnow=0.001         ! Roughness length for snow (m)
real, save         :: snowemiss=1.         ! snow emissitivity
real, save         :: maxsnowalpha=0.85    ! max snow albedo
real, save         :: minsnowalpha=0.5     ! min snow albedo
real, save         :: maxsnowden=300.      ! max snow density (kg m^-3)
real, save         :: minsnowden=100.      ! min snow density (kg m^-3)
! generic urban parameters
real, save         :: refheight=0.6        ! Displacement height as a fraction of building height (Kanda et al 2007)
real, save         :: zomratio=0.10        ! Ratio of roughness length to building height (default=0.1 or 10%)
real, save         :: zocanyon=0.01        ! Roughness length of in-canyon surfaces (m)
real, save         :: zoroof=0.01          ! Roughness length of roof surfaces (m)
real, save         :: maxrfwater=1.        ! Maximum roof water (kg m^-2)
real, save         :: maxrdwater=1.        ! Maximum road water (kg m^-2)
real, save         :: maxrfsn=1.           ! Maximum roof snow (kg m^-2)
real, save         :: maxrdsn=1.           ! Maximum road snow (kg m^-2)
real, save         :: maxvwatf=0.1         ! Factor multiplied to LAI to predict maximum leaf water (kg m^-2)
real, save         :: r_si=0.13            ! Building interior surface resistance (W^-1 m^2 K)
real, save         :: acfactor=1.          ! Air conditioning inefficiency factor
! atmosphere stability parameters
integer, save      :: icmax=5              ! number of iterations for stability functions (default=5)
real, save         :: a_1=1.
real, save         :: b_1=2./3.
real, save         :: c_1=5.
real, save         :: d_1=0.35

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
allocate(upack_g(ifull))
upack_g=sigu>0.
ufull_g=count(upack_g)
if (ufull_g==0) then
  deallocate(upack_g)
  return
end if

allocate(f_roof%depth(ufull_g,nl),f_roof%lambda(ufull_g,nl),f_roof%volcp(ufull_g,nl))
allocate(f_wall%depth(ufull_g,nl),f_wall%lambda(ufull_g,nl),f_wall%volcp(ufull_g,nl))
allocate(f_road%depth(ufull_g,nl),f_road%lambda(ufull_g,nl),f_road%volcp(ufull_g,nl))
allocate(f_slab%depth(ufull_g,nl),f_slab%lambda(ufull_g,nl),f_slab%volcp(ufull_g,nl))
allocate(f_intm%depth(ufull_g,nl),f_intm%lambda(ufull_g,nl),f_intm%volcp(ufull_g,nl))
allocate(f_roof%emiss(ufull_g),f_roof%alpha(ufull_g))
allocate(f_wall%emiss(ufull_g),f_wall%alpha(ufull_g))
allocate(f_road%emiss(ufull_g),f_road%alpha(ufull_g))
allocate(f_slab%emiss(ufull_g),f_ach(ufull_g))
allocate(roof%nodetemp(ufull_g,0:nl),road%nodetemp(ufull_g,0:nl),walle%nodetemp(ufull_g,0:nl))
allocate(wallw%nodetemp(ufull_g,0:nl))
allocate(slab%nodetemp(ufull_g,0:nl),intm%nodetemp(ufull_g,0:nl),room%nodetemp(ufull_g,1))
allocate(road%storage(ufull_g,nl),roof%storage(ufull_g,nl),walle%storage(ufull_g,nl),wallw%storage(ufull_g,nl))
allocate(slab%storage(ufull_g,nl),intm%storage(ufull_g,nl),room%storage(ufull_g,1))
allocate(cnveg%emiss(ufull_g),cnveg%sigma(ufull_g),cnveg%alpha(ufull_g))
allocate(rfveg%emiss(ufull_g),rfveg%sigma(ufull_g),rfveg%alpha(ufull_g))
allocate(cnveg%zo(ufull_g),cnveg%lai(ufull_g),cnveg%rsmin(ufull_g))
allocate(rfveg%zo(ufull_g),rfveg%lai(ufull_g),rfveg%rsmin(ufull_g))
allocate(f_rfvegdepth(ufull_g))
allocate(f_ctime(ufull_g),f_bldairtemp(ufull_g))
allocate(f_hangle(ufull_g),f_vangle(ufull_g),f_fbeam(ufull_g))
allocate(f_hwratio(ufull_g),f_coeffbldheight(ufull_g),f_effhwratio(ufull_g),f_bldheight(ufull_g))
allocate(f_sigmabld(ufull_g),f_industryfg(ufull_g),f_intgains_flr(ufull_g),f_trafficfg(ufull_g))
allocate(f_swilt(ufull_g),f_sfc(ufull_g),f_ssat(ufull_g))
allocate(p_lzom(ufull_g),p_lzoh(ufull_g),p_cndzmin(ufull_g),p_cduv(ufull_g),p_cdtq(ufull_g),cnveg%temp(ufull_g))
allocate(rfveg%temp(ufull_g))
allocate(p_tscrn(ufull_g),p_qscrn(ufull_g),p_uscrn(ufull_g),p_u10(ufull_g),p_emiss(ufull_g))
allocate(p_bldheat(ufull_g),p_bldcool(ufull_g),p_traf(ufull_g),p_intgains_full(ufull_g))
allocate(p_surferr(ufull_g),p_atmoserr(ufull_g),p_surferr_bias(ufull_g))
allocate(p_atmoserr_bias(ufull_g))
allocate(rfhyd%surfwater(ufull_g),rfhyd%snow(ufull_g),rfhyd%den(ufull_g),rfhyd%snowalpha(ufull_g))
allocate(rdhyd%surfwater(ufull_g),rdhyd%snow(ufull_g),rdhyd%den(ufull_g),rdhyd%snowalpha(ufull_g))
allocate(rdhyd%leafwater(ufull_g),rdhyd%soilwater(ufull_g),rfhyd%leafwater(ufull_g),rfhyd%soilwater(ufull_g))
allocate(sigmau_g(ufull_g),f_tempheat(ufull_g),f_tempcool(ufull_g))
allocate(int_viewf(ufull_g,4,4), int_psi(ufull_g,4,4),f_intmassn(ufull_g), f_bldwidth(ufull_g))
allocate(p_snowmelt(ufull_g))

! define grid arrays
sigmau_g = pack(sigu,upack_g)

iqu=0
iqut=0
do iq=1,ifull
  if (upack_g(iq)) then
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
roof%nodetemp=1.  ! + urbtemp
road%nodetemp=1.  ! + urbtemp
walle%nodetemp=1. ! + urbtemp
wallw%nodetemp=1. ! + urbtemp
rfhyd%surfwater=0.
rfhyd%snow=0.
rfhyd%den=minsnowden
rfhyd%snowalpha=maxsnowalpha
rfhyd%leafwater=0.
rdhyd%surfwater=0.
rdhyd%snow=0.
rdhyd%den=minsnowden
rdhyd%snowalpha=maxsnowalpha
rdhyd%leafwater=0.
rfhyd%soilwater=0.
rdhyd%soilwater=0.25

f_roof%depth=0.1
f_wall%depth=0.1
f_road%depth=0.1
f_rfvegdepth=0.1
f_roof%volcp=2.E6
f_wall%volcp=2.E6
f_road%volcp=2.E6
f_roof%lambda=2.
f_wall%lambda=2.
f_road%lambda=2.
f_hwratio=1.
f_sigmabld=0.5
cnveg%sigma=0.5
rfveg%sigma=0.
f_industryfg=0.
f_intgains_flr=0.
f_trafficfg=0.
f_bldheight=10.
f_roof%alpha=0.2
f_wall%alpha=0.2
f_road%alpha=0.2
cnveg%alpha=0.2
rfveg%alpha=0.2
f_roof%emiss=0.97
f_wall%emiss=0.97
f_road%emiss=0.97
cnveg%emiss=0.97
rfveg%emiss=0.97
f_bldairtemp=1. ! + urbtemp
f_vangle=0.
f_hangle=0.
f_ctime=0.
f_fbeam=1.
cnveg%zo=0.1
cnveg%lai=1.
cnveg%rsmin=200.
rfveg%zo=0.1
rfveg%lai=1.
rfveg%rsmin=200.
f_swilt=0.
f_sfc=0.5
f_ssat=1.
f_ach=0.5

slab%nodetemp=1. ! + urbtemp
intm%nodetemp=1. ! + urbtemp
room%nodetemp=1.  ! + urbtemp
f_slab%depth=0.1
f_slab%volcp=2.E6
f_slab%lambda=2.
f_slab%emiss=0.97
f_intm%depth=0.1
f_intm%lambda=2.
f_intm%volcp=2.E6
slab%storage=0._8
intm%storage=0._8
room%storage=0._8

utype=1 ! default urban
call atebtype(utype,diag)
call init_internal
call init_lwcoeff

p_cndzmin=max(10.,0.1*f_bldheight+2.)   ! updated in atebcalc
p_lzom=log(p_cndzmin/(0.1*f_bldheight)) ! updated in atebcalc
p_lzoh=6.+p_lzom ! (Kanda et al 2005)   ! updated in atebcalc
p_cduv=(vkar/p_lzom)**2                 ! updated in atebcalc
p_cdtq=vkar**2/(p_lzom*p_lzoh)          ! updated in atebcalc
cnveg%temp=1. ! + urbtemp               ! updated in atebcalc
rfveg%temp=1. ! + urbtemp               ! updated in atebcalc
p_tscrn=1.    ! + urbtemp               ! updated in atebcalc
p_qscrn=0.                              ! updated in atebcalc
p_uscrn=0.                              ! updated in atebcalc
p_u10=0.                                ! updated in atebcalc
p_emiss=0.97                            ! updated in atebcalc
p_bldheat=0._8
p_bldcool=0._8
p_traf=0._8
p_intgains_full=0._8
roof%storage =0._8
road%storage =0._8
walle%storage=0._8
wallw%storage=0._8
p_surferr=0._8
p_atmoserr=0._8
p_surferr_bias=0._8
p_atmoserr_bias=0._8

! for getqsat
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

return
end subroutine atebinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine deallocates arrays used by the TEB scheme

subroutine atebend(diag)

implicit none

integer, intent(in) :: diag

if (diag>=1) write(6,*) "Deallocating aTEB arrays"
if (ufull_g==0) return

deallocate(upack_g)
deallocate(f_roof%depth,f_wall%depth,f_road%depth,f_slab%depth,f_intm%depth)
deallocate(f_roof%volcp,f_wall%volcp,f_road%volcp,f_slab%volcp,f_intm%volcp)
deallocate(f_roof%lambda,f_wall%lambda,f_road%lambda,f_slab%lambda,f_intm%lambda)
deallocate(f_sigmabld,f_hwratio,f_bldheight,f_coeffbldheight,f_effhwratio)
deallocate(f_industryfg,f_intgains_flr,f_trafficfg,f_vangle,f_ctime,f_hangle,f_fbeam)
deallocate(f_roof%alpha,f_wall%alpha,f_road%alpha)
deallocate(f_roof%emiss,f_wall%emiss,f_road%emiss)
deallocate(f_slab%emiss)
deallocate(f_bldairtemp,cnveg%sigma,cnveg%alpha)
deallocate(cnveg%emiss,rfveg%sigma,f_rfvegdepth,rfveg%alpha,rfveg%emiss)
deallocate(cnveg%zo,cnveg%lai,cnveg%rsmin,rfveg%zo,rfveg%lai,rfveg%rsmin)
deallocate(f_swilt,f_sfc,f_ssat)
deallocate(p_lzom,p_lzoh,p_cndzmin,p_cduv,p_cdtq,cnveg%temp,rfveg%temp)
deallocate(p_tscrn,p_qscrn,p_uscrn,p_u10,p_emiss)
deallocate(p_surferr,p_atmoserr,p_surferr_bias,p_atmoserr_bias)
deallocate(p_bldheat,p_bldcool,p_traf,p_intgains_full)
deallocate(rfhyd%surfwater,rfhyd%snow,rfhyd%den,rfhyd%snowalpha)
deallocate(rdhyd%surfwater,rdhyd%snow,rdhyd%den,rdhyd%snowalpha)
deallocate(roof%nodetemp,road%nodetemp,walle%nodetemp,wallw%nodetemp)
deallocate(slab%nodetemp,intm%nodetemp,room%nodetemp)
deallocate(sigmau_g,rdhyd%leafwater,rdhyd%soilwater,rfhyd%leafwater,rfhyd%soilwater)
deallocate(int_viewf,int_psi)
deallocate(road%storage,roof%storage,walle%storage,wallw%storage)
deallocate(slab%storage,intm%storage,room%storage)
deallocate(f_intmassn,f_ach,f_tempheat,f_tempcool)
deallocate(p_snowmelt)

return
end subroutine atebend

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine loads aTEB state arrays (not compulsory)

subroutine atebload(urban,diag)

implicit none

integer, intent(in) :: diag
integer ii
real, dimension(ifull,4*nl+16), intent(in) :: urban

if (diag>=1) write(6,*) "Load aTEB state arrays"
if (ufull_g==0) return

do ii = 0,nl
  roof%nodetemp(:,ii) =pack(urban(:,0*nl+ii+1),   upack_g)
  where ( roof%nodetemp(:,ii)>100. )
    roof%nodetemp(:,ii) = roof%nodetemp(:,ii) - urbtemp
  end where
  walle%nodetemp(:,ii)=pack(urban(:,1*nl+ii+2), upack_g)
  where ( walle%nodetemp(:,ii)>100. )
    walle%nodetemp(:,ii) = walle%nodetemp(:,ii) - urbtemp
  end where
  wallw%nodetemp(:,ii)=pack(urban(:,2*nl+ii+3), upack_g)
  where ( wallw%nodetemp(:,ii)>100. )
    wallw%nodetemp(:,ii) = wallw%nodetemp(:,ii) - urbtemp
  end where
  road%nodetemp(:,ii) =pack(urban(:,3*nl+ii+4),upack_g)
  where ( road%nodetemp(:,ii)>100. )
    road%nodetemp(:,ii) = road%nodetemp(:,ii) - urbtemp
  end where
end do
rdhyd%soilwater=pack(urban(:,4*nl+5),upack_g)
rfhyd%soilwater=pack(urban(:,4*nl+6),upack_g)
rfhyd%surfwater=pack(urban(:,4*nl+7),upack_g)
rdhyd%surfwater=pack(urban(:,4*nl+8),upack_g)
rdhyd%leafwater=pack(urban(:,4*nl+9),upack_g)
rfhyd%leafwater=pack(urban(:,4*nl+10),upack_g)
rfhyd%snow     =pack(urban(:,4*nl+11),upack_g)
rdhyd%snow     =pack(urban(:,4*nl+12),upack_g)
rfhyd%den      =pack(urban(:,4*nl+13),upack_g)
rdhyd%den      =pack(urban(:,4*nl+14),upack_g)
rfhyd%snowalpha    =pack(urban(:,4*nl+15),upack_g)
rdhyd%snowalpha    =pack(urban(:,4*nl+16),upack_g)

return
end subroutine atebload

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! temperature only version of tebload

subroutine atebloadm(urban,moist,diag)

implicit none

integer, intent(in) :: diag
integer ii
real, dimension(ifull,4*nl+4), intent(in) :: urban
real, dimension(ifull,2), intent(in) :: moist

if (diag>=1) write(6,*) "Load aTEB state arrays"
if (ufull_g==0) return

do ii = 0,nl
  roof%nodetemp(:,ii) =pack(urban(:,0*nl+ii+1),   upack_g)
  where ( roof%nodetemp(:,ii)>100. )
    roof%nodetemp(:,ii) = roof%nodetemp(:,ii) - urbtemp
  end where
  walle%nodetemp(:,ii)=pack(urban(:,1*nl+ii+2), upack_g)
  where ( walle%nodetemp(:,ii)>100. )
    walle%nodetemp(:,ii) = walle%nodetemp(:,ii) - urbtemp
  end where
  wallw%nodetemp(:,ii)=pack(urban(:,2*nl+ii+3), upack_g)
  where ( wallw%nodetemp(:,ii)>100. )
    wallw%nodetemp(:,ii) = wallw%nodetemp(:,ii) - urbtemp
  end where
  road%nodetemp(:,ii) =pack(urban(:,3*nl+ii+4),upack_g)
  where ( road%nodetemp(:,ii)>100. )
    road%nodetemp(:,ii) = road%nodetemp(:,ii) - urbtemp
  end where
end do
rdhyd%soilwater=pack(moist(:,1),upack_g)
rfhyd%soilwater=pack(moist(:,2),upack_g)

return
end subroutine atebloadm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! general version of tebload

subroutine atebloadd(urban,mode,diag)

implicit none

integer, intent(in) :: diag
integer ii
real, dimension(ifull), intent(in) :: urban
character(len=*), intent(in) :: mode
character(len=10) :: teststr

if (diag>=1) write(6,*) "Load aTEB state array"
if (ufull_g==0) return

do ii = 0,nl
  write(teststr,'("rooftemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    roof%nodetemp(:,ii)=pack(urban,upack_g)
    where ( roof%nodetemp(:,ii)>100. )
      roof%nodetemp(:,ii) = roof%nodetemp(:,ii) - urbtemp  
    end where    
    return
  end if
  write(teststr,'("walletemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    walle%nodetemp(:,ii)=pack(urban,upack_g)
    where ( walle%nodetemp(:,ii)>100. )
      walle%nodetemp(:,ii) = walle%nodetemp(:,ii) - urbtemp  
    end where  
    return
  end if
  write(teststr,'("wallwtemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    wallw%nodetemp(:,ii)=pack(urban,upack_g)
    where ( wallw%nodetemp(:,ii)>100. )
      wallw%nodetemp(:,ii) = wallw%nodetemp(:,ii) - urbtemp  
    end where  
    return
  end if
  write(teststr,'("roadtemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    road%nodetemp(:,ii)=pack(urban,upack_g)
    where ( road%nodetemp(:,ii)>100. )
      road%nodetemp(:,ii) = road%nodetemp(:,ii) - urbtemp  
    end where  
    return
  end if
  write(teststr,'("slabtemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    slab%nodetemp(:,ii)=pack(urban,upack_g)
    where ( slab%nodetemp(:,ii)>100. )
      slab%nodetemp(:,ii) = slab%nodetemp(:,ii) - urbtemp  
    end where  
    return
  end if  
  write(teststr,'("intmtemp",I1.1)') ii+1
  if ( trim(teststr)==trim(mode) ) then
    intm%nodetemp(:,ii)=pack(urban,upack_g)
    where ( intm%nodetemp(:,ii)>100. )
      intm%nodetemp(:,ii) = intm%nodetemp(:,ii) - urbtemp  
    end where  
    return
  end if   
end do  
  
select case(mode)
  case("canyonsoilmoisture")
    rdhyd%soilwater=pack(urban,upack_g)
    return
  case("roofsoilmoisture")
    rfhyd%soilwater=pack(urban,upack_g)
    return
  case("roadsurfacewater")
    rdhyd%surfwater=pack(urban,upack_g)  
    return
  case("roofsurfacewater")
    rfhyd%surfwater=pack(urban,upack_g)  
    return
  case("canyonleafwater")
    rdhyd%leafwater=pack(urban,upack_g)  
    return
  case("roofleafwater")
    rfhyd%leafwater=pack(urban,upack_g)  
    return
  case("roadsnowdepth")
    rdhyd%snow=pack(urban,upack_g)  
    return
  case("roofsnowdepth")
    rfhyd%snow=pack(urban,upack_g)  
    return
  case("roadsnowdensity")
    rdhyd%den=pack(urban,upack_g)  
    return
  case("roofsnowdensity")
    rfhyd%den=pack(urban,upack_g)  
    return
  case("roadsnowalbedo")
    rdhyd%snowalpha=pack(urban,upack_g)  
    return
  case("roofsnowalbedo")
    rfhyd%snowalpha=pack(urban,upack_g)  
    return
end select
  
write(6,*) "ERROR: Unknown mode for atebloadd ",trim(mode)
stop
    
return
end subroutine atebloadd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine loads aTEB type arrays (not compulsory)

subroutine atebtype(itype,diag)

implicit none

integer, intent(in) :: diag
integer ii,j,ierr,nlp
integer, dimension(ifull), intent(in) :: itype
integer, dimension(ufull_g) :: itmp
integer, parameter :: maxtype = 8
real x
real, dimension(ufull_g) :: tsigveg,tsigmabld
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
! Internal gains sensible heat flux [floor] (W m^-2)
real, dimension(maxtype) :: cintgains_flr=(/ 5.,   5.,   5.,   5.,   5.,   5.,   5.,   5. /)
! Daily averaged traffic sensible heat flux (W m^-2)
real, dimension(maxtype) ::  ctrafficfg=(/  1.5,  1.5,  1.5,  1.5,  1.5,  1.5,  1.5,  1.5 /)
! Comfort temperature (K)
real, dimension(maxtype) :: cbldtemp=(/ 291.16, 291.16, 291.16, 291.16, 291.16, 291.16, 291.16, 291.16 /)
! Roof albedo
real, dimension(maxtype) ::  croofalpha=(/ 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20 /)    ! (Fortuniak 08) Masson = 0.15
! Wall albedo
real, dimension(maxtype) ::  cwallalpha=(/ 0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30 /)    ! (Fortuniak 08) Masson = 0.25
! Road albedo
real, dimension(maxtype) ::  croadalpha=(/ 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10 /)    ! (Fortuniak 08) Masson = 0.08
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
! Slab emissitivity
real, dimension(maxtype) ::  cslabemiss=(/ 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90 /) 
! Canyon veg emissitivity
real, dimension(maxtype) ::  cvegemissc=(/ 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96 /)
! Roof veg emissitivity
real, dimension(maxtype) ::  cvegemissr=(/ 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96 /)
! Green roof soil depth
real, dimension(maxtype) ::   cvegdeptr=(/ 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10 /)
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
! Air volume changes per hour (m^3 m^-3)
real, dimension(maxtype) ::       cach=(/  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,  0.50 /)
! Comfort temperature for heating [k]
real, dimension(maxtype) ::  ctempheat=(/  288.,  288.,  288.,  288.,  288.,  0.,  0.,  0. /)
! Comfort temperature for cooling [k]
real, dimension(maxtype) ::  ctempcool=(/  296.,  296.,  296.,  296.,  296.,  999.,  999.,  999. /)

real, dimension(maxtype,nl) :: croofdepth
real, dimension(maxtype,nl) :: cwalldepth
real, dimension(maxtype,nl) :: croaddepth
real, dimension(maxtype,nl) :: cslabdepth
real, dimension(maxtype,nl) :: croofcp
real, dimension(maxtype,nl) :: cwallcp
real, dimension(maxtype,nl) :: croadcp
real, dimension(maxtype,nl) :: cslabcp
real, dimension(maxtype,nl) :: crooflambda
real, dimension(maxtype,nl) :: cwalllambda
real, dimension(maxtype,nl) :: croadlambda
real, dimension(maxtype,nl) :: cslablambda

namelist /atebnml/  resmeth,useonewall,zohmeth,acmeth,intairtmeth,intmassmeth,nrefl,vegmode,soilunder, &
                    conductmeth,scrnmeth,wbrelaxc,wbrelaxr,lweff,iqt
namelist /atebsnow/ zosnow,snowemiss,maxsnowalpha,minsnowalpha,maxsnowden,minsnowden
namelist /atebgen/  refheight,zomratio,zocanyon,zoroof,maxrfwater,maxrdwater,maxrfsn,maxrdsn,maxvwatf, &
                    acfactor
namelist /atebtile/ czovegc,cvegrlaic,cvegrsminc,czovegr,cvegrlair,cvegrsminr,cswilt,csfc,cssat,       &
                    cvegemissc,cvegemissr,cvegdeptr,cvegalphac,cvegalphar,csigvegc,csigvegr,           &
                    csigmabld,cbldheight,chwratio,cindustryfg,cintgains_flr,ctrafficfg,cbldtemp,       &
                    croofalpha,cwallalpha,croadalpha,croofemiss,cwallemiss,croademiss,croofdepth,      &
                    cwalldepth,croaddepth,croofcp,cwallcp,croadcp,crooflambda,cwalllambda,croadlambda, &
                    cslabdepth,cslabcp,cslablambda,cach,ctempheat,ctempcool

! facet array where: rows=maxtypes (landtypes) and columns=nl (material layers)
nlp=nl/4 ! number of layers in each material segment (over 4 material segments)
! depths (m)
croofdepth= reshape((/ ((0.01/nlp, ii=1,maxtype),j=1,nlp),    &
                       ((0.09/nlp, ii=1,maxtype),j=1,nlp),    & 
                       ((0.40/nlp, ii=1,maxtype),j=1,nlp),    &
                       ((0.10/nlp, ii=1,maxtype),j=1,nlp) /), &
                       (/maxtype,nl/))
cwalldepth= reshape((/ ((0.01/nlp, ii=1,maxtype),j=1,nlp),    &
                       ((0.04/nlp, ii=1,maxtype),j=1,nlp),    & 
                       ((0.10/nlp, ii=1,maxtype),j=1,nlp),    &
                       ((0.05/nlp, ii=1,maxtype),j=1,nlp) /), &
                       (/maxtype,nl/))
croaddepth= reshape((/ ((0.01/nlp, ii=1,maxtype),j=1,nlp),    &
                       ((0.04/nlp, ii=1,maxtype),j=1,nlp),    & 
                       ((0.45/nlp, ii=1,maxtype),j=1,nlp),    &
                       ((3.50/nlp, ii=1,maxtype),j=1,nlp) /), &
                       (/maxtype,nl/))
cslabdepth=reshape((/  ((0.05/nlp, ii=1,maxtype),j=1,nlp),    &
                       ((0.05/nlp, ii=1,maxtype),j=1,nlp),    & 
                       ((0.05/nlp, ii=1,maxtype),j=1,nlp),    &
                       ((0.05/nlp, ii=1,maxtype),j=1,nlp) /), &
                       (/maxtype,nl/))
! heat capacity (J m^-3 K^-1)
croofcp =   reshape((/ ((2.11E6, ii=1,maxtype),j=1,nlp),    & ! dense concrete (Oke 87)
                       ((2.11E6, ii=1,maxtype),j=1,nlp),    & ! dense concrete (Oke 87)
                       ((0.28E6, ii=1,maxtype),j=1,nlp),    & ! light concrete (Oke 87)
                       ((0.29E6, ii=1,maxtype),j=1,nlp) /), & ! insulation (Oke 87)
                       (/maxtype,nl/))
cwallcp =   reshape((/ ((1.55E6, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                       ((1.55E6, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                       ((1.55E6, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                       ((0.29E6, ii=1,maxtype),j=1,nlp) /), & ! insulation (Oke 87)
                       (/maxtype,nl/))
croadcp =   reshape((/ ((1.94E6, ii=1,maxtype),j=1,nlp),    & ! asphalt (Mills 93)
                       ((1.94E6, ii=1,maxtype),j=1,nlp),    & ! asphalt (Mills 93)
                       ((1.28E6, ii=1,maxtype),j=1,nlp),    & ! dry soil (Mills 93)
                       ((1.28E6, ii=1,maxtype),j=1,nlp) /), & ! dry soil (Mills 93)
                       (/maxtype,nl/))
cslabcp=   reshape((/  ((1.55E6, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                       ((1.55E6, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                       ((1.55E6, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                       ((1.55E6, ii=1,maxtype),j=1,nlp) /), & ! concrete (Mills 93)
                       (/maxtype,nl/))
! heat conductivity (W m^-1 K^-1)
crooflambda=reshape((/ ((1.5100, ii=1,maxtype),j=1,nlp),    & ! dense concrete (Oke 87)
                       ((1.5100, ii=1,maxtype),j=1,nlp),    & ! dense concrete (Oke 87)
                       ((0.0800, ii=1,maxtype),j=1,nlp),    & ! light concrete (Oke 87)
                       ((0.0500, ii=1,maxtype),j=1,nlp) /), & ! insulation (Oke 87)
                       (/maxtype,nl/))
cwalllambda=reshape((/ ((0.9338, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                       ((0.9338, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                       ((0.9338, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                       ((0.0500, ii=1,maxtype),j=1,nlp) /), & ! insulation (Oke 87)
                       (/maxtype,nl/))
croadlambda=reshape((/ ((0.7454, ii=1,maxtype),j=1,nlp),    & ! asphalt (Mills 93)
                       ((0.7454, ii=1,maxtype),j=1,nlp),    & ! asphalt (Mills 93)
                       ((0.2513, ii=1,maxtype),j=1,nlp),    & ! dry soil (Mills 93)
                       ((0.2513, ii=1,maxtype),j=1,nlp) /), & ! dry soil (Mills 93)
                       (/maxtype,nl/))
cslablambda=reshape((/ ((0.9338, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                       ((0.9338, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                       ((0.9338, ii=1,maxtype),j=1,nlp),    & ! concrete (Mills 93)
                       ((0.9338, ii=1,maxtype),j=1,nlp) /), & ! concrete (Mills 93)
                       (/maxtype,nl/))

if (diag>=1) write(6,*) "Load aTEB building properties"
if (ufull_g==0) return

itmp=pack(itype,upack_g)
if ((minval(itmp)<1).or.(maxval(itmp)>maxtype)) then
  write(6,*) "ERROR: Urban type is out of range"
  stop
end if

if (atebnmlfile/=0) then
  open(unit=atebnmlfile,file='ateb.nml',action="read",iostat=ierr)
  if (ierr==0) then
    write(6,*) "Reading ateb.nml"
    read(atebnmlfile,nml=atebnml)
    read(atebnmlfile,nml=atebsnow)
    read(atebnmlfile,nml=atebgen)
    read(atebnmlfile,nml=atebtile)
    close(atebnmlfile)  
  end if
end if

select case(vegmode)
  case(0)
    tsigveg=0.5*csigvegc(itmp)/(1.-0.5*csigvegc(itmp))
    tsigmabld=csigmabld(itmp)/(1.-0.5*csigvegc(itmp))
    sigmau_g=sigmau_g*(1.-0.5*csigvegc(itmp))
  case(1)
    tsigveg=0.
    tsigmabld=csigmabld(itmp)/(1.-csigvegc(itmp))
    sigmau_g=sigmau_g*(1.-csigvegc(itmp))
  case(2)
    tsigveg=csigvegc(itmp)
    tsigmabld=csigmabld(itmp)
  case DEFAULT
    if (vegmode<0) then
      x=real(abs(vegmode))/100.
      x=max(min(x,1.),0.)
      tsigveg=x*csigvegc(itmp)/(1.-(1.-x)*csigvegc(itmp))
      tsigmabld=csigmabld(itmp)/(1.-(1.-x)*csigvegc(itmp))
      sigmau_g=sigmau_g*(1.-(1.-x)*csigvegc(itmp))
    else
      write(6,*) "ERROR: Unsupported vegmode ",vegmode
      stop
    end if
end select
cnveg%sigma=max(min(tsigveg/(1.-tsigmabld),1.),0.)
rfveg%sigma=max(min(csigvegr(itmp),1.),0.)
f_sigmabld=max(min(tsigmabld,1.),0.)
!f_hwratio=chwratio(itmp)*f_sigmabld/(1.-f_sigmabld) ! MJT suggested new definition
f_hwratio=chwratio(itmp)          ! MJL simple definition

f_industryfg=cindustryfg(itmp)
f_intgains_flr=cintgains_flr(itmp)
f_trafficfg=ctrafficfg(itmp)
f_bldheight=cbldheight(itmp)
f_roof%alpha=croofalpha(itmp)
f_wall%alpha=cwallalpha(itmp)
f_road%alpha=croadalpha(itmp)
cnveg%alpha=cvegalphac(itmp)
rfveg%alpha=cvegalphar(itmp)
f_roof%emiss=croofemiss(itmp)
f_wall%emiss=cwallemiss(itmp)
f_road%emiss=croademiss(itmp)
cnveg%emiss=cvegemissc(itmp)
rfveg%emiss=cvegemissr(itmp)
f_bldairtemp=cbldtemp(itmp) - urbtemp
! room%nodetemp(:,1)=cbldtemp(itmp) - urbtemp
f_rfvegdepth=cvegdeptr(itmp)
do ii=1,nl
  f_roof%depth(:,ii)=croofdepth(itmp,ii)
  f_wall%depth(:,ii)=cwalldepth(itmp,ii)
  f_road%depth(:,ii)=croaddepth(itmp,ii)
  f_roof%lambda(:,ii)=crooflambda(itmp,ii)
  f_wall%lambda(:,ii)=cwalllambda(itmp,ii)
  f_road%lambda(:,ii)=croadlambda(itmp,ii)
  f_roof%volcp(:,ii)=croofcp(itmp,ii)
  f_wall%volcp(:,ii)=cwallcp(itmp,ii)
  select case(soilunder)
    case(0) ! storage under road only
      f_road%volcp(:,ii)=croadcp(itmp,ii)
    case(1) ! storage under road and canveg
      f_road%volcp(:,ii)=croadcp(itmp,ii)/(1.-cnveg%sigma)
    case(2) ! storage under road and bld
      f_road%volcp(:,ii)=croadcp(itmp,ii)*(1./(1.-cnveg%sigma)*(1./(1.-f_sigmabld)-1.) +1.)
    case(3) ! storage under road and canveg and bld (100% of grid point)
      f_road%volcp(:,ii)=croadcp(itmp,ii)/(1.-cnveg%sigma)/(1.-f_sigmabld)
    case DEFAULT
      write(6,*) "ERROR: Unknown soilunder mode ",soilunder
      stop
  end select
end do
cnveg%zo=czovegc(itmp)
cnveg%lai=cvegrlaic(itmp)
cnveg%rsmin=cvegrsminc(itmp)/max(cnveg%lai,1.E-8)
rfveg%zo=czovegr(itmp)
rfveg%lai=cvegrlair(itmp)
rfveg%rsmin=cvegrsminr(itmp)/max(rfveg%lai,1.E-8)
f_swilt=cswilt(itmp)
f_sfc=csfc(itmp)
f_ssat=cssat(itmp)

! for varying internal temperature
f_slab%emiss=cslabemiss(itmp)
f_ach = cach(itmp)
f_tempheat = ctempheat(itmp)
f_tempcool = ctempcool(itmp)
do ii=1,nl
  f_slab%depth(:,ii)=cslabdepth(itmp,ii)
  f_intm%depth(:,ii)=cslabdepth(itmp,ii)
  f_slab%lambda(:,ii)=cslablambda(itmp,ii)
  f_intm%lambda(:,ii)=cslablambda(itmp,ii)
  f_slab%volcp(:,ii)=cslabcp(itmp,ii)
  f_intm%volcp(:,ii)=cslabcp(itmp,ii)
end do

! Here we modify the effective canyon geometry to account for in-canyon vegetation tall vegetation
f_coeffbldheight = max(f_bldheight-6.*cnveg%zo,0.2)/f_bldheight
f_effhwratio   = f_hwratio*f_coeffbldheight

call init_internal
call init_lwcoeff

if ( diag>0 ) then
  write(6,*) 'hwratio, eff',f_hwratio, f_effhwratio
  write(6,*) 'bldheight, eff',f_bldheight, f_coeffbldheight
  write(6,*) 'sigmabld, sigmavegc', f_sigmabld, cnveg%sigma
  write(6,*) 'roadcp multiple for soilunder:', soilunder,f_road%volcp(itmp,1)/croadcp(itmp,1)
end if

return
end subroutine atebtype

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine atebdeftype(paramdata,typedata,paramname,diag)

implicit none

integer, parameter :: maxtype = 8
integer, intent(in) :: diag
integer, dimension(ifull), intent(in) :: typedata
integer, dimension(ufull_g) :: itmp
real, dimension(maxtype), intent(in) :: paramdata
character(len=*), intent(in) :: paramname

if ( diag>=1 ) write(6,*) "Load aTEB parameters ",trim(paramname)
if ( ufull_g==0 ) return

itmp = pack(typedata,upack_g)
if ( minval(itmp)<1 .or. maxval(itmp)>maxtype ) then
  write(6,*) "ERROR: Urban type is out of range"
  stop
end if

select case(paramname)
  case('bldheight')
    f_bldheight = paramdata(itmp)  
  case('hwratio')
    f_hwratio=paramdata(itmp)  
  case('sigvegc')
    cnveg%sigma=paramdata(itmp)/(1.-f_sigmabld)  
  case('sigmabld')
    f_sigmabld=paramdata(itmp)  
  case('industryfg')
    f_industryfg=paramdata(itmp)  
  case('trafficfg')
    f_trafficfg=paramdata(itmp)  
  case('roofalpha')
    f_roof%alpha=paramdata(itmp)  
  case('wallalpha')
    f_wall%alpha=paramdata(itmp) 
  case('roadalpha')
    f_road%alpha=paramdata(itmp)  
  case('vegalphac')
    cnveg%alpha=paramdata(itmp)  
  case('zovegc')
    cnveg%zo=paramdata(itmp)  
  case default
    write(6,*) "ERROR: Unknown aTEB parameter name ",trim(paramname)
    stop
end select

call init_internal
call init_lwcoeff

return
end subroutine atebdeftype

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine specifies the urban properties for each grid point
!

subroutine atebfndef(ifn,diag)

implicit none

integer, intent(in) :: diag
integer ii
real, dimension(ifull,9*nl+27), intent(in) :: ifn

if (diag>=1) write(6,*) "Load aTEB building properties"
if (ufull_g==0) return

f_hwratio    = pack(ifn(:,1),upack_g)
f_sigmabld   = pack(ifn(:,2),upack_g)
cnveg%sigma  = pack(ifn(:,3)/(1.-ifn(:,2)),upack_g)
rfveg%sigma  = pack(ifn(:,4),upack_g)
f_industryfg = pack(ifn(:,5),upack_g)
f_trafficfg  = pack(ifn(:,6),upack_g)
f_bldheight  = pack(ifn(:,7),upack_g)
f_roof%alpha = pack(ifn(:,8),upack_g)
f_wall%alpha = pack(ifn(:,9),upack_g)
f_road%alpha = pack(ifn(:,10),upack_g)
cnveg%alpha  = pack(ifn(:,11),upack_g)
cnveg%alpha  = pack(ifn(:,12),upack_g)
f_roof%emiss = pack(ifn(:,13),upack_g)
f_wall%emiss = pack(ifn(:,14),upack_g)
f_road%emiss = pack(ifn(:,15),upack_g)
cnveg%emiss  = pack(ifn(:,16),upack_g)
rfveg%emiss  = pack(ifn(:,17),upack_g)
f_bldairtemp = pack(ifn(:,18)-urbtemp,upack_g)

do ii=1,nl
  f_roof%depth(:,ii)   = pack(ifn(:,0*nl+ii+18),upack_g)
  f_wall%depth(:,ii)  = pack(ifn(:,1*nl+ii+18),upack_g)
  f_road%depth(:,ii)   = pack(ifn(:,2*nl+ii+18),upack_g)
  f_roof%volcp(:,ii)   = pack(ifn(:,3*nl+ii+18),upack_g)
  f_wall%volcp(:,ii)  = pack(ifn(:,4*nl+ii+18),upack_g)
  f_road%volcp(:,ii)   = pack(ifn(:,5*nl+ii+18),upack_g)
  f_roof%lambda(:,ii)  = pack(ifn(:,6*nl+ii+18),upack_g)
  f_wall%lambda(:,ii) = pack(ifn(:,7*nl+ii+18),upack_g)
  f_road%lambda(:,ii)  = pack(ifn(:,8*nl+ii+18),upack_g)
end do
cnveg%zo    = pack(ifn(:,9*nl+19),upack_g)
cnveg%lai   = pack(ifn(:,9*nl+20),upack_g)
cnveg%rsmin = pack(ifn(:,9*nl+21),upack_g)
rfveg%zo    = pack(ifn(:,9*nl+22),upack_g)
rfveg%lai   = pack(ifn(:,9*nl+23),upack_g)
rfveg%rsmin = pack(ifn(:,9*nl+24),upack_g)
f_swilt     = pack(ifn(:,9*nl+25),upack_g)
f_sfc       = pack(ifn(:,9*nl+26),upack_g)
f_ssat      = pack(ifn(:,9*nl+27),upack_g)

call init_internal
call init_lwcoeff

return
end subroutine atebfndef

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine saves aTEB state arrays (not compulsory)

subroutine atebsave(urban,diag,rawtemp)

implicit none

integer, intent(in) :: diag
integer ii
real, dimension(ifull,9*nl+27), intent(inout) :: urban
logical, intent(in), optional :: rawtemp
logical rawmode

if ( diag>=1 ) write(6,*) "Save aTEB state arrays"
if ( ufull_g==0 ) return

rawmode = .false.
if ( present(rawtemp) ) then
  rawmode = rawtemp
end if

if ( rawmode ) then                                                                ! if nl=4 then index:
  do ii=0,nl
    urban(:,0*nl+ii+1)=unpack(roof%nodetemp(:,ii),upack_g,urban(:,0*nl+ii+1))        ! 1:5
    urban(:,1*nl+ii+2)=unpack(walle%nodetemp(:,ii),upack_g,urban(:,1*nl+ii+2))       ! 6:10
    urban(:,2*nl+ii+3)=unpack(wallw%nodetemp(:,ii),upack_g,urban(:,2*nl+ii+3))       ! 11:15
    urban(:,3*nl+ii+4)=unpack(road%nodetemp(:,ii),upack_g,urban(:,3*nl+ii+4))        ! 16:20
  end do
else
  do ii=0,nl
    urban(:,0*nl+ii+1)=unpack(roof%nodetemp(:,ii)+urbtemp,upack_g,urban(:,0*nl+ii+1))    ! 1:5
    urban(:,1*nl+ii+2)=unpack(walle%nodetemp(:,ii)+urbtemp,upack_g,urban(:,1*nl+ii+2))   ! 6:10
    urban(:,2*nl+ii+3)=unpack(wallw%nodetemp(:,ii)+urbtemp,upack_g,urban(:,2*nl+ii+3))   ! 11:15
    urban(:,3*nl+ii+4)=unpack(road%nodetemp(:,ii)+urbtemp,upack_g,urban(:,3*nl+ii+4))    ! 16:20
  end do
end if
urban(:,4*nl+5)=unpack(rdhyd%soilwater(:),upack_g,urban(:,4*nl+5))            ! 21
urban(:,4*nl+6)=unpack(rfhyd%soilwater(:),upack_g,urban(:,4*nl+6))            ! 22
urban(:,4*nl+7)=unpack(rfhyd%surfwater(:),upack_g,urban(:,4*nl+7))            ! 23
urban(:,4*nl+8)=unpack(rdhyd%surfwater(:),upack_g,urban(:,4*nl+8))            ! 24
urban(:,4*nl+9)=unpack(rdhyd%leafwater(:),upack_g,urban(:,4*nl+9))            ! 25
urban(:,4*nl+10)=unpack(rfhyd%leafwater(:),upack_g,urban(:,4*nl+10))          ! 26
urban(:,4*nl+11)=unpack(rfhyd%snow(:), upack_g,urban(:,4*nl+11))              ! 27
urban(:,4*nl+12)=unpack(rdhyd%snow(:), upack_g,urban(:,4*nl+12))              ! 28
urban(:,4*nl+13)=unpack(rfhyd%den(:),  upack_g,urban(:,4*nl+13))              ! 29
urban(:,4*nl+14)=unpack(rdhyd%den(:),  upack_g,urban(:,4*nl+14))              ! 30
urban(:,4*nl+15)=unpack(rfhyd%snowalpha(:),upack_g,urban(:,4*nl+15))          ! 31
urban(:,4*nl+16)=unpack(rdhyd%snowalpha(:),upack_g,urban(:,4*nl+16))          ! 32

return
end subroutine atebsave

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! temperature only version of tebsave

subroutine atebsavem(urban,moist,diag,rawtemp)

implicit none

integer, intent(in) :: diag
integer ii
real, dimension(ifull,4*nl+4), intent(inout) :: urban
real, dimension(ifull,2), intent(inout) :: moist
logical, intent(in), optional :: rawtemp
logical rawmode

if ( diag>=1 ) write(6,*) "Save aTEB state arrays"
if ( ufull_g==0 ) return

rawmode = .false.
if ( present(rawtemp) ) then
  rawmode = rawtemp
end if

if ( rawmode ) then                                                                ! if nl=4 then index:
  do ii=0,nl
    urban(:,0*nl+ii+1)=unpack(roof%nodetemp(:,ii),upack_g,urban(:,0*nl+ii+1))            ! 1:5
    urban(:,1*nl+ii+2)=unpack(walle%nodetemp(:,ii),upack_g,urban(:,1*nl+ii+2))           ! 6:10
    urban(:,2*nl+ii+3)=unpack(wallw%nodetemp(:,ii),upack_g,urban(:,2*nl+ii+3))           ! 11:15
    urban(:,3*nl+ii+4)=unpack(road%nodetemp(:,ii),upack_g,urban(:,3*nl+ii+4))            ! 16:20
  end do
else
  do ii=0,nl
    urban(:,0*nl+ii+1)=unpack(roof%nodetemp(:,ii)+urbtemp,upack_g,urban(:,0*nl+ii+1))    ! 1:5
    urban(:,1*nl+ii+2)=unpack(walle%nodetemp(:,ii)+urbtemp,upack_g,urban(:,1*nl+ii+2))   ! 6:10
    urban(:,2*nl+ii+3)=unpack(wallw%nodetemp(:,ii)+urbtemp,upack_g,urban(:,2*nl+ii+3))   ! 11:15
    urban(:,3*nl+ii+4)=unpack(road%nodetemp(:,ii)+urbtemp,upack_g,urban(:,3*nl+ii+4))    ! 16:20
  end do
end if
moist(:,1)=unpack(rdhyd%soilwater(:),upack_g,moist(:,1))
moist(:,2)=unpack(rfhyd%soilwater(:),upack_g,moist(:,2))

return
end subroutine atebsavem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! general version of atebsave

subroutine atebsaved(urban,mode,diag,rawtemp)

implicit none

integer, intent(in) :: diag
integer ii
real, dimension(ifull), intent(out) :: urban
logical, intent(in), optional :: rawtemp
logical rawmode
character(len=*), intent(in) :: mode
character(len=10) :: teststr

if (diag>=1) write(6,*) "Load aTEB state array"
if (ufull_g==0) return

rawmode = .false.
if ( present(rawtemp) ) then
  rawmode = rawtemp
end if

if ( rawmode ) then
  do ii = 0,nl
    write(teststr,'("rooftemp",I1.1)') ii+1
    if ( trim(teststr)==trim(mode) ) then
      urban=unpack(roof%nodetemp(:,ii),upack_g,urban)
      return
    end if
    write(teststr,'("walletemp",I1.1)') ii+1
    if ( trim(teststr)==trim(mode) ) then
      urban=unpack(walle%nodetemp(:,ii),upack_g,urban)
      return
    end if
    write(teststr,'("wallwtemp",I1.1)') ii+1
    if ( trim(teststr)==trim(mode) ) then
      urban=unpack(wallw%nodetemp(:,ii),upack_g,urban)
      return
    end if
    write(teststr,'("roadtemp",I1.1)') ii+1
    if ( trim(teststr)==trim(mode) ) then
      urban=unpack(road%nodetemp(:,ii),upack_g,urban)
      return
    end if
    write(teststr,'("slabtemp",I1.1)') ii+1
    if ( trim(teststr)==trim(mode) ) then
      urban=unpack(slab%nodetemp(:,ii),upack_g,urban)
      return
    end if  
    write(teststr,'("intmtemp",I1.1)') ii+1
    if ( trim(teststr)==trim(mode) ) then
      urban=unpack(intm%nodetemp(:,ii),upack_g,urban)
      return
    end if   
  end do  
else
  do ii = 0,nl
    write(teststr,'("rooftemp",I1.1)') ii+1
    if ( trim(teststr)==trim(mode) ) then
      urban=unpack(roof%nodetemp(:,ii)+urbtemp,upack_g,urban)
      return
    end if
    write(teststr,'("walletemp",I1.1)') ii+1
    if ( trim(teststr)==trim(mode) ) then
      urban=unpack(walle%nodetemp(:,ii)+urbtemp,upack_g,urban)
      return
    end if
    write(teststr,'("wallwtemp",I1.1)') ii+1
    if ( trim(teststr)==trim(mode) ) then
      urban=unpack(wallw%nodetemp(:,ii)+urbtemp,upack_g,urban)
      return
    end if
    write(teststr,'("roadtemp",I1.1)') ii+1
    if ( trim(teststr)==trim(mode) ) then
      urban=unpack(road%nodetemp(:,ii)+urbtemp,upack_g,urban)
      return
    end if
    write(teststr,'("slabtemp",I1.1)') ii+1
    if ( trim(teststr)==trim(mode) ) then
      urban=unpack(slab%nodetemp(:,ii)+urbtemp,upack_g,urban)
      return
    end if  
    write(teststr,'("intmtemp",I1.1)') ii+1
    if ( trim(teststr)==trim(mode) ) then
      urban=unpack(intm%nodetemp(:,ii)+urbtemp,upack_g,urban)
      return
    end if   
  end do  
end if

select case(mode)
  case("canyonsoilmoisture")
    urban=unpack(rdhyd%soilwater,upack_g,urban)  
    return
  case("roofsoilmoisture")
    urban=unpack(rfhyd%soilwater,upack_g,urban)    
    return
  case("roadsurfacewater")
    urban=unpack(rdhyd%surfwater,upack_g,urban)    
    return
  case("roofsurfacewater")
    urban=unpack(rfhyd%surfwater,upack_g,urban)    
    return
  case("canyonleafwater")
    urban=unpack(rdhyd%leafwater,upack_g,urban)    
    return
  case("roofleafwater")
    urban=unpack(rfhyd%leafwater,upack_g,urban)    
    return
  case("roadsnowdepth")
    urban=unpack(rdhyd%snow,upack_g,urban)    
    return
  case("roofsnowdepth")
    urban=unpack(rfhyd%snow,upack_g,urban)    
    return
  case("roadsnowdensity")
    urban=unpack(rdhyd%den,upack_g,urban)    
    return
  case("roofsnowdensity")
    urban=unpack(rfhyd%den,upack_g,urban)    
    return
  case("roadsnowalbedo")
    urban=unpack(rdhyd%snowalpha,upack_g,urban)    
    return
  case("roofsnowalbedo")
    urban=unpack(rfhyd%snowalpha,upack_g,urban)    
    return
end select

write(6,*) "ERROR: Unknown mode for atebsaved ",trim(mode)
stop

return
end subroutine atebsaved

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine collects and passes energy closure information to atebwrap

subroutine energyrecord(o_atmoserr,o_atmoserr_bias,o_surferr,o_surferr_bias, &
                        o_heating,o_cooling,o_intgains,o_traf,o_bldtemp)

implicit none

real, dimension(ufull_g), intent(out) :: o_atmoserr,o_atmoserr_bias,o_surferr,o_surferr_bias
real, dimension(ufull_g), intent(out) :: o_heating,o_cooling,o_intgains,o_traf,o_bldtemp

if ( ufull_g==0 ) return

p_atmoserr_bias = p_atmoserr_bias + p_atmoserr
p_surferr_bias = p_surferr_bias + p_surferr

o_atmoserr      = real(pack(p_atmoserr,upack_g))
o_surferr       = real(pack(p_surferr,upack_g))
o_atmoserr_bias = real(pack(p_atmoserr_bias,upack_g))
o_surferr_bias  = real(pack(p_surferr_bias,upack_g))
o_heating       = pack(p_bldheat,upack_g)
o_cooling       = pack(p_bldcool,upack_g)
o_intgains      = pack(p_intgains_full,upack_g)
o_traf          = pack(p_traf,upack_g)
o_bldtemp       = pack(room%nodetemp(:,1)+urbtemp,upack_g)

return
end subroutine energyrecord

subroutine atebenergy(o_data,mode,if_industryfg,p_bldheat,p_bldcool,p_traf,p_intgains_full,sigmau,upack,ufull,diag)

#ifdef CCAM
use cc_omp                         ! CC OpenMP routines
#endif

implicit none

integer, intent(in) :: ufull, diag
real, dimension(:), intent(inout) :: o_data
real, dimension(ufull) :: ctmp, dtmp
character(len=*), intent(in) :: mode
real, dimension(ufull), intent(in) :: if_industryfg
real, dimension(ufull), intent(in) :: p_bldheat, p_bldcool, p_traf, p_intgains_full
real, dimension(ufull), intent(in) :: sigmau
logical, dimension(size(o_data)), intent(in) :: upack

if ( diag>=1 .and. ntiles==1 ) write(6,*) "Extract energy output"
if ( ufull==0 ) return

select case(mode)
  case("anthropogenic")
    ctmp = pack(o_data, upack)
    dtmp = p_bldheat + p_bldcool + p_traf + if_industryfg + p_intgains_full
    ctmp = (1.-sigmau)*ctmp + sigmau*dtmp
    o_data = unpack(ctmp, upack, o_data)
  case default
    write(6,*) "ERROR: Unknown atebenergy mode ",trim(mode)
    stop
end select    

return
end subroutine atebenergy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine blends urban momentum and heat roughness lengths
! (This version neglects the displacement height (e.g., for CCAM))
!

subroutine atebzo(zom,zoh,zoq,p_cndzmin,p_lzom,p_lzoh,sigmau,upack,ufull,diag,raw)

#ifdef CCAM
use cc_omp                         ! CC OpenMP routines
#endif

implicit none

integer, intent(in) :: ufull, diag
real, dimension(:), intent(inout) :: zom, zoh, zoq
real, dimension(ufull) :: workb,workc,workd,zmtmp,zhtmp,zqtmp
real, parameter :: zr=1.e-15 ! limits minimum roughness length for heat
logical, intent(in), optional :: raw
logical mode
real, dimension(ufull), intent(in) :: p_cndzmin, p_lzom, p_lzoh
real, dimension(ufull), intent(in) :: sigmau
logical, dimension(size(zom)), intent(in) :: upack

if ( diag>=1 .and. ntiles==1 ) write(6,*) "Calculate urban roughness lengths"
if ( ufull==0 ) return

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

subroutine atebcd(cduv,cdtq,p_cdtq,p_cduv,sigmau,upack,ufull,diag,raw)
 
#ifdef CCAM
use cc_omp                         ! CC OpenMP routines
#endif

implicit none
 
integer, intent(in) :: ufull, diag
real, dimension(:), intent(inout) :: cduv, cdtq
real, dimension(ufull) :: ctmp
logical, intent(in), optional :: raw
logical outmode
real, dimension(ufull), intent(in) :: p_cdtq, p_cduv
real, dimension(ufull), intent(in) :: sigmau
logical, dimension(size(cduv)), intent(in) :: upack
 
if ( diag>=1 .and. ntiles==1 ) write(6,*) "Calculate urban drag coeff"
if ( ufull==0 ) return
 
outmode=.false.
if (present(raw)) outmode=raw
 
ctmp=pack(cduv,upack)
if ( outmode ) then
  ctmp=p_cduv 
else
  ctmp=(1.-sigmau)*ctmp+sigmau*p_cduv
end if
cduv=unpack(ctmp,upack,cduv)
 
ctmp=pack(cdtq,upack)
if ( outmode ) then
  ctmp=p_cdtq 
else
  ctmp=(1.-sigmau)*ctmp+sigmau*p_cdtq
end if
cdtq=unpack(ctmp,upack,cdtq)
 
return
end subroutine atebcd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine is for hydrological outputs
!
 
subroutine atebhydro(hydroout,mode,p_snowmelt,sigmau,upack,ufull,diag)

#ifdef CCAM
use cc_omp                         ! CC OpenMP routines
#endif
 
implicit none
 
integer, intent(in) :: ufull, diag
real, dimension(:), intent(inout) :: hydroout
real, dimension(ufull) :: ctmp
character(len=*), intent(in) :: mode
real, dimension(ufull), intent(in) :: p_snowmelt
real, dimension(ufull), intent(in) :: sigmau
logical, dimension(size(hydroout)), intent(in) :: upack
 
if ( diag>=1 .and. ntiles==1 ) write(6,*) "Calculate hydrological outputs"
if ( ufull==0 ) return
 
select case(mode)
  case("snowmelt")
    ctmp=pack(hydroout,upack)
    ctmp=(1.-sigmau)*ctmp+sigmau*p_snowmelt
    hydroout=unpack(ctmp,upack,hydroout)
  case default
    write(6,*) "ERROR: Unknown atebhydro mode ",trim(mode)
    stop
end select
 
return
end subroutine atebhydro

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Store fraction of direct radiation
!

subroutine atebfbeam(is,ifin,fbeam,diag)

implicit none

integer, intent(in) :: is,ifin,diag
integer ifinish,ib,ie,ucount
real, dimension(ifin), intent(in) :: fbeam

if ( diag>=1 ) write(6,*) "Assign urban direct beam ratio"
if ( ufull_g==0 ) return

ifinish=is+ifin-1
ucount=count(upack_g(is:ifinish))
if (ucount==0) return

ib=count(upack_g(1:is-1))+1
ie=ucount+ib-1
f_fbeam(ib:ie)=pack(fbeam,upack_g(is:ifinish))

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
real, dimension(ufull_g) :: tmpr,tmpk,tmprat
real, dimension(ufull_g) :: lsg,lcosin
real, intent(in) :: fjd
real, parameter :: solcon = 1370.

if ( diag>=1 ) write(6,*) "Diagnose urban direct beam ratio"
if ( ufull_g==0 ) return

ifinish=is+ifin-1
ucount=count(upack_g(is:ifinish))
if (ucount==0) return

ib=count(upack_g(1:is-1))+1
ie=ucount+ib-1

lsg(ib:ie)   =pack(sg,upack_g(is:ifinish))
lcosin(ib:ie)=pack(cosin,upack_g(is:ifinish))

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
real, dimension(ufull_g) :: ualb,utmp
logical, intent(in), optional :: raw
logical outmode

if ( diag>=1 ) write(6,*) "Calculate urban albedo (broad)"
if ( ufull_g==0 ) return

outmode=.false.
if (present(raw)) outmode=raw

albmode=0 ! net albedo
if (present(split)) albmode=split

ifinish=is+ifin-1
ucount=count(upack_g(is:ifinish))
if (ucount==0) return

ib=count(upack_g(1:is-1))+1
ie=ucount+ib-1
call atebalbcalc(ib,ucount,ualb(ib:ie),albmode,diag)

if (outmode) then
  alb(:)=unpack(ualb(ib:ie),upack_g(is:ifinish),alb)
else
  utmp(ib:ie)=pack(alb,upack_g(is:ifinish))
  utmp(ib:ie)=(1.-sigmau_g(ib:ie))*utmp(ib:ie)+sigmau_g(ib:ie)*ualb(ib:ie)
  alb(:)=unpack(utmp(ib:ie),upack_g(is:ifinish),alb)
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
real, dimension(ifin) :: snowdeltac, snowdeltar
real, dimension(ifin) :: wallpsi,roadpsi
real, dimension(ifin) :: sg_roof,sg_vegr,sg_road,sg_walle,sg_wallw,sg_vegc,sg_rfsn,sg_rdsn
real, dimension(ifin) :: dumfbeam

if ( diag>=1 ) write(6,*) "Calculate urban albedo"

ie=ifin+is-1

select case(albmode)
  case default ! net albedo
    dumfbeam=f_fbeam(is:ie)
  case(1)      ! direct albedo
    dumfbeam=1.
  case(2)      ! diffuse albedo
    dumfbeam=0.
  end select

! roof
snowdeltar=rfhyd%snow(is:ie)/(rfhyd%snow(is:ie)+maxrfsn)
  
! canyon
snowdeltac=rdhyd%snow(is:ie)/(rdhyd%snow(is:ie)+maxrdsn)
call getswcoeff(sg_roof,sg_vegr,sg_road,sg_walle,sg_wallw,sg_vegc,sg_rfsn,sg_rdsn,wallpsi,roadpsi,f_effhwratio,    &
                f_vangle(is:ie),f_hangle(is:ie),dumfbeam,cnveg%sigma(is:ie),f_road%alpha(is:ie),cnveg%alpha(is:ie), &
                f_wall%alpha(is:ie),rdhyd%snowalpha(is:ie),snowdeltac)
sg_walle=sg_walle*f_coeffbldheight(is:ie)
sg_wallw=sg_wallw*f_coeffbldheight(is:ie)

call getnetalbedo(alb,sg_roof,sg_vegr,sg_road,sg_walle,sg_wallw,sg_vegc,sg_rfsn,sg_rdsn,                       &
                  f_hwratio(is:ie),f_sigmabld(is:ie),rfveg%sigma(is:ie),f_roof%alpha(is:ie),rfveg%alpha(is:ie), &
                  cnveg%sigma(is:ie),f_road%alpha(is:ie),f_wall%alpha(is:ie),cnveg%alpha(is:ie),                 &
                  rfhyd%snowalpha(is:ie),rdhyd%snowalpha(is:ie),snowdeltar,snowdeltac)

return
end subroutine atebalbcalc

subroutine getnetalbedo(alb,sg_roof,sg_vegr,sg_road,sg_walle,sg_wallw,sg_vegc,sg_rfsn,sg_rdsn,  &
                        if_hwratio,if_sigmabld,if_vegsigmar,if_roofalpha,if_vegalphar,          &
                        if_vegsigmac,if_roadalpha,if_wallalpha,if_vegalphac,                    &
                        roofalpha,roadalpha,snowdeltar,snowdeltac)

implicit none

real, dimension(:), intent(out) :: alb
real, dimension(size(alb)), intent(in) :: sg_roof, sg_vegr, sg_road, sg_walle, sg_wallw, sg_vegc
real, dimension(size(alb)), intent(in) :: sg_rfsn, sg_rdsn
real, dimension(size(alb)), intent(in) :: if_hwratio, if_sigmabld
real, dimension(size(alb)), intent(in) :: if_vegsigmar, if_roofalpha, if_vegalphar
real, dimension(size(alb)), intent(in) :: if_vegsigmac, if_roadalpha, if_vegalphac, if_wallalpha
real, dimension(size(alb)), intent(in) :: roofalpha, roadalpha, snowdeltar, snowdeltac
real, dimension(size(alb)) :: albu, albr

! canyon
albu=1.-(if_hwratio*(sg_walle+sg_wallw)*(1.-if_wallalpha)+snowdeltac*sg_rdsn*(1.-roadalpha)                 &
    +(1.-snowdeltac)*((1.-if_vegsigmac)*sg_road*(1.-if_roadalpha)+if_vegsigmac*sg_vegc*(1.-if_vegalphac)))

! roof
albr=(1.-snowdeltar)*((1.-if_vegsigmar)*sg_roof*if_roofalpha+if_vegsigmar*sg_vegr*if_vegalphar) &
    +snowdeltar*sg_rfsn*roofalpha

! net
alb=if_sigmabld*albr+(1.-if_sigmabld)*albu

return
end subroutine getnetalbedo

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

if (ufull_g==0) return

ifinish=is+ifin-1
ucount=count(upack_g(is:ifinish))
if (ucount==0) return

ib=count(upack_g(1:is-1))+1
ie=ucount+ib-1

f_hangle(ib:ie)=0.5*pi-pack(azimuthin,upack_g(is:ifinish))
f_vangle(ib:ie)=acos(pack(cosin,upack_g(is:ifinish)))
f_ctime(ib:ie)=pack(ctimein,upack_g(is:ifinish))

return
end subroutine atebnewangle1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This version of tebnewangle is for CCAM and TAPM
!

subroutine atebccangle(is,ifin,cosin,rlon,rlat,fjd,slag,dt,sdlt)

implicit none

integer, intent(in) :: is,ifin
integer ifinish,ucount,ib,ie
real, intent(in) :: fjd,slag,dt,sdlt
real cdlt
real, dimension(ifin), intent(in) :: cosin,rlon,rlat
real, dimension(ufull_g) :: hloc,x,y,lattmp

! cosin = cosine of zenith angle
! rlon = longitude
! rlat = latitude
! fjd = day of year
! slag = sun lag angle
! sdlt = sin declination of sun

if ( ufull_g==0 ) return

ifinish=is+ifin-1
ucount=count(upack_g(is:ifinish))
if (ucount==0) return

ib=count(upack_g(1:is-1))+1
ie=ucount+ib-1

cdlt=sqrt(min(max(1.-sdlt*sdlt,0.),1.))

lattmp(ib:ie)=pack(rlat,upack_g(is:ifinish))

! from CCAM zenith.f
hloc(ib:ie)=2.*pi*fjd+slag+pi+pack(rlon,upack_g(is:ifinish))+dt*pi/86400.
! estimate azimuth angle
x(ib:ie)=sin(-hloc(ib:ie))*cdlt
y(ib:ie)=-cos(-hloc(ib:ie))*cdlt*sin(lattmp(ib:ie))+cos(lattmp(ib:ie))*sdlt
!azimuth=atan2(x,y)
f_hangle(ib:ie)=0.5*pi-atan2(x(ib:ie),y(ib:ie))
f_vangle(ib:ie)=acos(pack(cosin,upack_g(is:ifinish)))
f_ctime(ib:ie)=min(max(mod(0.5*hloc(ib:ie)/pi-0.5,1.),0.),1.)

return
end subroutine atebccangle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calcuates screen level diagnostics
!

subroutine atebscrnout(tscrn,qscrn,uscrn,u10,diag,raw)

implicit none

integer, intent(in) :: diag
real, dimension(ifull), intent(inout) :: tscrn,qscrn,uscrn,u10
real, dimension(ufull_g) :: tmp
logical, intent(in), optional :: raw
logical mode

if (diag>=1) write(6,*) "Calculate urban 2m diagnostics"
if (ufull_g==0) return

mode=.false.
if (present(raw)) mode=raw

if (mode) then
  tscrn=unpack(p_tscrn+urbtemp,upack_g,tscrn)
  qscrn=unpack(p_qscrn,upack_g,qscrn)
  uscrn=unpack(p_uscrn,upack_g,uscrn)
  u10  =unpack(p_u10,  upack_g,u10  )
else
  tmp=pack(tscrn,upack_g)
  tmp=sigmau_g*(p_tscrn+urbtemp)+(1.-sigmau_g)*tmp
  tscrn=unpack(tmp,upack_g,tscrn)
  tmp=pack(qscrn,upack_g)
  tmp=sigmau_g*p_qscrn+(1.-sigmau_g)*tmp
  qscrn=unpack(tmp,upack_g,qscrn)
  tmp=pack(uscrn,upack_g)
  tmp=sigmau_g*p_uscrn+(1.-sigmau_g)*tmp
  uscrn=unpack(tmp,upack_g,uscrn)
  tmp=pack(u10,upack_g)
  tmp=sigmau_g*p_u10+(1.-sigmau_g)*tmp
  u10=unpack(tmp,upack_g,u10)
end if

return
end subroutine atebscrnout

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Extract urban fraction
subroutine atebsigmau(sigu,diag)

implicit none

integer, intent(in) :: diag
real, dimension(ifull), intent(out) :: sigu

if (diag>=1) write(6,*) "Calculate urban cover fraction"
sigu=0.
if (ufull_g==0) return
sigu=unpack(sigmau_g,upack_g,sigu)

return
end subroutine atebsigmau

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main routine for calculating urban flux contrabution

! ifull = number of horizontal grid points
! dt = model time step (sec)
! zmin = first model level height (m)
! sg = incoming short wave radiation (W/m^2)
! rg = incoming long wave radiation (W/m^2)
! rnd = incoming rainfall/snowfall rate (kg/(m^2 s))
! rho = atmospheric density at first model level (kg/m^3)
! temp = atmospheric temperature at first model level (K)
! mixr = atmospheric mixing ratio at first model level (kg/kg)
! ps = surface pressure (Pa)
! pa = pressure at first model level (Pa)
! uu = U component of wind speed at first model level (m/s)
! vv = V component of wind speed at first model level (m/s)
! umin = minimum wind speed (m/s)
! ofg = Input/Output sensible heat flux (W/m^2)
! oeg = Input/Output latent heat flux (W/m^2)
! ots = Input/Output radiative/skin temperature (K)
! owf = Input/Output wetness fraction/surface water (%)
! diag = diagnostic message mode (0=off, 1=basic messages, 2=more detailed messages, etc)

subroutine atebcalc(ofg,oeg,ots,owf,orn,dt,zmin,sg,rg,rnd,snd,rho,temp,mixr,ps,uu,vv,umin,diag,raw)

implicit none

integer, intent(in) :: diag
real, intent(in) :: dt,umin
real, dimension(ifull), intent(in) :: sg,rg,rnd,snd,rho,temp,mixr,ps,uu,vv,zmin
real, dimension(ifull), intent(inout) :: ofg,oeg,ots,owf,orn
real, dimension(ufull_g) :: tmp
real, dimension(ufull_g) :: a_sg,a_rg,a_rho,a_temp,a_mixr,a_ps,a_umag,a_udir,a_rnd,a_snd,a_zmin
real, dimension(ufull_g) :: u_fg,u_eg,u_ts,u_wf,u_rn
logical, intent(in), optional :: raw
logical mode

if ( ufull_g==0 ) return ! no urban grid points

! mode = .false. implies weight output with urban area cover fraction
! mode = .true. implies no weighting of output with urban area cover fraction (assumes 100% cover)
mode=.false.
if (present(raw)) mode=raw

! Host model meteorological data
a_zmin=pack(zmin,                 upack_g)
a_sg  =pack(sg,                   upack_g)
a_rg  =pack(rg,                   upack_g)
a_rho =pack(rho,                  upack_g)
a_temp=pack(temp-urbtemp,         upack_g)
a_mixr=pack(mixr,                 upack_g)
a_ps  =pack(ps,                   upack_g)
a_umag=max(pack(sqrt(uu*uu+vv*vv),upack_g),umin)
a_udir=pack(atan2(vv,uu),         upack_g)
a_rnd =pack(rnd-snd,              upack_g)
a_snd =pack(snd,                  upack_g)

! Update urban prognostic variables
call atebeval(u_fg,u_eg,u_ts,u_wf,u_rn,dt,a_sg,a_rg,a_rho,a_temp,a_mixr,a_ps,a_umag,a_udir,a_rnd,a_snd,a_zmin,      &
              f_bldheight,f_bldwidth,f_coeffbldheight,f_ctime,f_effhwratio,f_fbeam,f_hangle,f_hwratio,              &
              f_industryfg,f_intgains_flr,f_intm,f_intmassn,f_rfvegdepth,f_road,f_roof,f_sfc,f_sigmabld,f_slab,     &
              f_ssat,f_swilt,f_trafficfg,f_vangle,f_wall,intm,p_cdtq,p_cduv,p_cndzmin,p_emiss,p_intgains_full,      &
              p_lzoh,p_lzom,p_snowmelt,p_traf,rdhyd,rfhyd,rfveg,road,roof,room,slab,walle,wallw,cnveg,p_atmoserr,   &
              p_bldcool,p_bldheat,p_surferr,int_psi,int_viewf,p_qscrn,p_tscrn,p_u10,p_uscrn,f_ach,f_tempcool,       &
              f_tempheat,f_bldairtemp,ufull_g,diag)

! export urban fluxes on host grid
if (mode) then
  ofg=unpack(u_fg,upack_g,ofg)
  oeg=unpack(u_eg,upack_g,oeg)
  ots=unpack(u_ts+urbtemp,upack_g,ots)
  owf=unpack(u_wf,upack_g,owf)
  orn=unpack(u_rn,upack_g,orn)
else
  tmp=pack(ofg,upack_g)
  tmp=(1.-sigmau_g)*tmp+sigmau_g*u_fg
  ofg=unpack(tmp,upack_g,ofg)
  tmp=pack(oeg,upack_g)
  tmp=(1.-sigmau_g)*tmp+sigmau_g*u_eg
  oeg=unpack(tmp,upack_g,oeg)
  tmp=pack(ots,upack_g)
  tmp=((1.-sigmau_g)*tmp**4+sigmau_g*(u_ts+urbtemp)**4)**0.25
  ots=unpack(tmp,upack_g,ots)
  tmp=pack(owf,upack_g)
  tmp=(1.-sigmau_g)*tmp+sigmau_g*u_wf
  owf=unpack(tmp,upack_g,owf)
  tmp=pack(orn,upack_g)
  tmp=(1.-sigmau_g)*tmp+sigmau_g*u_rn
  orn=unpack(tmp,upack_g,orn)
end if

return
end subroutine atebcalc

subroutine atebcalc_thread(ofg,oeg,ots,owf,orn,dt,zmin,sg,rg,rnd,snd,rho,temp,mixr,ps,uu,vv,    &
                    umin,sigmau,if_bldheight,if_bldwidth,if_coeffbldheight,if_ctime,            &
                    if_effhwratio,if_fbeam,if_hangle,if_hwratio,if_industryfg,if_intgains_flr,  &
                    if_intm,if_intmassn,if_rfvegdepth,if_road,if_roof,if_sfc,if_sigmabld,       &
                    if_slab,if_ssat,if_swilt,if_trafficfg,if_vangle,if_wall,intm,p_cdtq,p_cduv, &
                    p_cndzmin,p_emiss,p_intgains_full,p_lzoh,p_lzom,p_snowmelt,p_traf,rdhyd,    &
                    rfhyd,rfveg,road,roof,room,slab,walle,wallw,cnveg,p_atmoserr,p_bldcool,     &
                    p_bldheat,p_surferr,int_psi,int_viewf,p_qscrn,p_tscrn,p_u10,p_uscrn,if_ach, &
                    if_tempcool,if_tempheat,if_bldairtemp,upack,ufull,diag,raw)

implicit none

integer, intent(in) :: ufull, diag
real, intent(in) :: dt, umin
real, dimension(:), intent(in) :: sg,rg,rnd,snd,rho,temp,mixr,ps,uu,vv,zmin
real, dimension(:), intent(inout) :: ofg,oeg,ots,owf,orn
real, dimension(ufull) :: tmp
real, dimension(ufull) :: a_sg,a_rg,a_rho,a_temp,a_mixr,a_ps,a_umag,a_udir,a_rnd,a_snd,a_zmin
real, dimension(ufull) :: u_fg,u_eg,u_ts,u_wf,u_rn
logical, intent(in), optional :: raw
logical mode
real, dimension(ufull), intent(in) :: sigmau
logical, dimension(size(sg)), intent(in) :: upack
real, dimension(ufull), intent(in) :: if_bldheight, if_bldwidth, if_coeffbldheight, if_ctime
real, dimension(ufull), intent(in) :: if_effhwratio, if_fbeam, if_hangle, if_hwratio, if_industryfg
real, dimension(ufull), intent(in) :: if_intgains_flr, if_rfvegdepth, if_sfc, if_sigmabld, if_ssat
real, dimension(ufull), intent(in) :: if_swilt, if_trafficfg, if_vangle, if_ach, if_tempcool, if_tempheat
real, dimension(ufull), intent(in) :: if_bldairtemp
integer, dimension(ufull), intent(in) :: if_intmassn
type(facetparams), intent(in) :: if_intm, if_road, if_roof, if_slab, if_wall
real, dimension(ufull), intent(inout) :: p_cdtq, p_cduv, p_cndzmin, p_emiss, p_intgains_full
real, dimension(ufull), intent(inout) :: p_lzoh, p_lzom, p_snowmelt, p_traf, p_bldcool, p_bldheat
real, dimension(ufull), intent(inout) :: p_qscrn, p_tscrn, p_u10, p_uscrn
real(kind=8), dimension(ufull), intent(inout) :: p_atmoserr, p_surferr
type(hydrodata), intent(inout) :: rdhyd, rfhyd
type(vegdata), intent(inout) :: rfveg
type(facetdata), intent(inout) :: road, roof, room, slab, walle, wallw, intm
type(vegdata), intent(inout) :: cnveg
real(kind=8), dimension(ufull,4,4), intent(in) :: int_psi, int_viewf

if ( ufull==0 ) return ! no urban grid points

! mode = .false. implies weight output with urban area cover fraction
! mode = .true. implies no weighting of output with urban area cover fraction (assumes 100% cover)
mode=.false.
if (present(raw)) mode=raw

! Host model meteorological data
a_zmin=pack(zmin,                 upack)
a_sg  =pack(sg,                   upack)
a_rg  =pack(rg,                   upack)
a_rho =pack(rho,                  upack)
a_temp=pack(temp-urbtemp,         upack)
a_mixr=pack(mixr,                 upack)
a_ps  =pack(ps,                   upack)
a_umag=max(pack(sqrt(uu*uu+vv*vv),upack),umin)
a_udir=pack(atan2(vv,uu),         upack)
a_rnd =pack(rnd-snd,              upack)
a_snd =pack(snd,                  upack)

! Update urban prognostic variables
call atebeval(u_fg,u_eg,u_ts,u_wf,u_rn,dt,a_sg,a_rg,a_rho,a_temp,a_mixr,a_ps,a_umag,a_udir,a_rnd,a_snd,a_zmin,       &
              if_bldheight,if_bldwidth,if_coeffbldheight,if_ctime,if_effhwratio,if_fbeam,if_hangle,if_hwratio,       &
              if_industryfg,if_intgains_flr,if_intm,if_intmassn,if_rfvegdepth,if_road,if_roof,if_sfc,if_sigmabld,    &
              if_slab,if_ssat,if_swilt,if_trafficfg,if_vangle,if_wall,intm,p_cdtq,p_cduv,p_cndzmin,p_emiss,          &
              p_intgains_full,p_lzoh,p_lzom,p_snowmelt,p_traf,rdhyd,rfhyd,rfveg,road,roof,room,slab,walle,wallw,     &
              cnveg,p_atmoserr,p_bldcool,p_bldheat,p_surferr,int_psi,int_viewf,p_qscrn,p_tscrn,p_u10,p_uscrn,        &
              if_ach,if_tempcool,if_tempheat,if_bldairtemp,ufull,diag)

! export urban fluxes on host grid
if (mode) then
  ofg=unpack(u_fg,upack,ofg)
  oeg=unpack(u_eg,upack,oeg)
  ots=unpack(u_ts+urbtemp,upack,ots)
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
  tmp=((1.-sigmau)*tmp**4+sigmau*(u_ts+urbtemp)**4)**0.25
  ots=unpack(tmp,upack,ots)
  tmp=pack(owf,upack)
  tmp=(1.-sigmau)*tmp+sigmau*u_wf
  owf=unpack(tmp,upack,owf)
  tmp=pack(orn,upack)
  tmp=(1.-sigmau)*tmp+sigmau*u_rn
  orn=unpack(tmp,upack,orn)
end if

return
end subroutine atebcalc_thread

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

subroutine atebeval(u_fg,u_eg,u_ts,u_wf,u_rn,ddt,a_sg,a_rg,a_rho,a_temp,a_mixr,a_ps,a_umag,a_udir,a_rnd,a_snd,a_zmin, &
                    if_bldheight,if_bldwidth,if_coeffbldheight,if_ctime,if_effhwratio,if_fbeam,if_hangle,if_hwratio,  &
                    if_industryfg,if_intgains_flr,if_intm,if_intmassn,if_rfvegdepth,if_road,if_roof,if_sfc,           &
                    if_sigmabld,if_slab,if_ssat,if_swilt,if_trafficfg,if_vangle,if_wall,intm,p_cdtq,p_cduv,p_cndzmin, &
                    p_emiss,p_intgains_full,p_lzoh,p_lzom,p_snowmelt,p_traf,rdhyd,rfhyd,rfveg,road,roof,room,slab,    &
                    walle,wallw,cnveg,p_atmoserr,p_bldcool,p_bldheat,p_surferr,int_psi,int_viewf,p_qscrn,p_tscrn,     &
                    p_u10,p_uscrn,if_ach,if_tempcool,if_tempheat,if_bldairtemp,ufull,diag)

#ifdef CCAM
use cc_omp                         ! CC OpenMP routines
#endif

implicit none

integer, intent(in) :: ufull
integer, intent(in) :: diag
integer k
real, intent(in) :: ddt
real, dimension(ufull), intent(in) :: a_sg,a_rg,a_rho,a_temp,a_mixr,a_ps,a_umag,a_udir,a_rnd,a_snd,a_zmin
real, dimension(ufull), intent(out) :: u_fg,u_eg,u_ts,u_wf,u_rn
real, dimension(ufull) :: ggint_roof,ggint_walle,ggint_wallw,ggint_road,ggint_slab,ggint_intm2
real, dimension(ufull) :: rdsntemp,rfsntemp,rdsnmelt,rfsnmelt,garfsn,gardsn
real, dimension(ufull) :: wallpsi,roadpsi,fgtop,egtop,qsatr,qsata
real, dimension(ufull) :: cu,fgrooftop,egrooftop
real, dimension(ufull) :: we,ww,wr,zolog,a,n,zom,zonet,dis
real, dimension(ufull) :: roofvegwetfac,roadvegwetfac
real, dimension(ufull) :: z_on_l,pa,dts,dtt
real, dimension(ufull) :: u_alb, u_melt
real, dimension(ufull) :: sg_roof,sg_vegr,sg_road,sg_walle,sg_wallw,sg_vegc,sg_rfsn,sg_rdsn
real, dimension(ufull) :: rg_roof,rg_road,rg_walle,rg_wallw,rg_vegc,rg_vegr,rg_rfsn,rg_rdsn
real, dimension(ufull) :: rgint_roof,rgint_walle,rgint_wallw,rgint_slab,rgint_zero
real, dimension(ufull) :: fg_roof,fg_road,fg_walle,fg_wallw,fg_vegc,fg_vegr,fg_rfsn,fg_rdsn
real, dimension(ufull) :: eg_roof,eg_road,eg_vegc,eg_vegr,eg_rfsn,eg_rdsn
real, dimension(ufull) :: acond_roof,acond_road,acond_walle,acond_wallw
real, dimension(ufull) :: acond_vegc,acond_vegr,acond_rfsn,acond_rdsn
real, dimension(ufull) :: abase_road,abase_walle,abase_wallw,abase_vegc,abase_rdsn
real, dimension(ufull) :: d_roofdelta,d_roaddelta,d_vegdeltac,d_vegdeltar,d_rfsndelta,d_rdsndelta
real, dimension(ufull) :: d_tempc,d_mixrc,d_tempr,d_mixrr,d_sigd,d_sigr,d_rfdzmin
real, dimension(ufull) :: d_ac_outside,d_canyonrgout,d_roofrgout,d_tranc,d_evapc,d_tranr,d_evapr,d_c1c,d_c1r
real, dimension(ufull) :: d_totdepth,d_netemiss,d_netrad,d_topu
real, dimension(ufull) :: d_cwa,d_cw0,d_cww,d_cwr,d_cra,d_crr,d_crw
real, dimension(ufull) :: d_canyontemp,d_canyonmix,d_traf
real, dimension(ufull) :: ggext_roof,ggext_walle,ggext_wallw,ggext_road,ggext_slab,ggint_intm1,ggext_impl
real, dimension(ufull) :: int_newairtemp, d_ac_inside, d_intgains_bld, int_infilflux
real, dimension(ufull) :: cyc_traffic,cyc_basedemand,cyc_proportion,cyc_translation
real, dimension(ufull) :: ggint_intm1_temp
real, dimension(ufull) :: int_infilfg
real, dimension(ufull,nl) :: depth_cp, depth_lambda 
real, dimension(ufull), intent(in) :: if_bldheight, if_bldwidth, if_coeffbldheight, if_ctime, if_effhwratio
real, dimension(ufull), intent(in) :: if_fbeam, if_hangle, if_hwratio, if_industryfg, if_intgains_flr, if_rfvegdepth
real, dimension(ufull), intent(in) :: if_sfc, if_sigmabld, if_ssat, if_swilt, if_trafficfg, if_vangle, if_ach, if_tempcool
real, dimension(ufull), intent(in) :: if_tempheat, if_bldairtemp
integer, dimension(ufull), intent(in) :: if_intmassn
type(facetparams), intent(in) :: if_intm, if_road, if_roof, if_slab, if_wall
real, dimension(ufull), intent(inout) :: p_cdtq, p_cduv, p_cndzmin, p_emiss, p_intgains_full, p_lzoh, p_lzom, p_snowmelt
real, dimension(ufull), intent(inout) :: p_traf, p_bldcool, p_bldheat, p_qscrn, p_tscrn, p_u10, p_uscrn
real(kind=8), dimension(ufull), intent(inout) :: p_atmoserr, p_surferr
type(hydrodata), intent(inout) :: rdhyd, rfhyd
type(vegdata), intent(inout) :: rfveg
type(facetdata), intent(inout) :: road, roof, room, slab, walle, wallw, intm
type(vegdata), intent(inout) :: cnveg
real(kind=8), dimension(ufull,4,4), intent(in) :: int_psi, int_viewf

if ( diag>=1 .and. ntiles==1 ) write(6,*) "Evaluating aTEB"

! new snowfall
where ( a_snd>1.e-10 )
  ! update snow density
  rfhyd%den = (rfhyd%snow*rfhyd%den+a_snd*ddt*minsnowden)/(rfhyd%snow+ddt*a_snd)
  rdhyd%den = (rdhyd%snow*rdhyd%den+a_snd*ddt*minsnowden)/(rdhyd%snow+ddt*a_snd)
  ! reset snow albedo
  rfhyd%snowalpha = maxsnowalpha
  rdhyd%snowalpha = maxsnowalpha
end where

! calculate water and snow area cover fractions
d_roofdelta = max(rfhyd%surfwater/maxrfwater,0.)**(2./3.)
d_roaddelta = max(rdhyd%surfwater/maxrdwater,0.)**(2./3.)
d_vegdeltac = max(rdhyd%leafwater/max(maxvwatf*cnveg%lai,1.E-8),0.)**(2./3.)
d_vegdeltar = max(rfhyd%leafwater/max(maxvwatf*rfveg%lai,1.E-8),0.)**(2./3.)
d_rfsndelta = rfhyd%snow/(rfhyd%snow+maxrfsn)
d_rdsndelta = rdhyd%snow/(rdhyd%snow+maxrdsn)

! canyon level air temp and water vapor (displacement height at refheight*building height)
pa      = a_ps*exp(-grav*a_zmin/(rd*(a_temp+urbtemp)))
d_sigd  = a_ps
a       = (d_sigd/pa)**(rd/aircp)
d_tempc = a_temp*a + urbtemp*(a-1.)
call getqsat(qsatr,d_tempc,d_sigd)
call getqsat(qsata,a_temp,pa)
d_mixrc = a_mixr*qsatr/qsata

! roof level air temperature and water vapor (displacement height at building height)
d_sigr  = a_ps*exp(-grav*if_bldheight*(1.-refheight)/(rd*(a_temp+urbtemp)))
a       = (d_sigr/pa)**(rd/aircp)
d_tempr = a_temp*a + urbtemp*(a-1.)
call getqsat(qsatr,d_tempr,d_sigr)
d_mixrr = a_mixr*qsatr/qsata

! calculate soil data
d_totdepth = sum(if_road%depth,2)
call getc1(d_c1c,ufull)
call getc1(d_c1r,ufull)

! calculate shortwave reflections
! Here we modify the effective canyon geometry to account for in-canyon vegetation
call getswcoeff(sg_roof,sg_vegr,sg_road,sg_walle,sg_wallw,sg_vegc,sg_rfsn,sg_rdsn,wallpsi,roadpsi,if_effhwratio,  &
                if_vangle,if_hangle,if_fbeam,cnveg%sigma,if_road%alpha,cnveg%alpha,if_wall%alpha,rdhyd%snowalpha, &
                d_rdsndelta)
sg_walle = sg_walle*if_coeffbldheight ! shadow due to in-canyon vegetation
sg_wallw = sg_wallw*if_coeffbldheight ! shadow due to in-canyon vegetation
call getnetalbedo(u_alb,sg_roof,sg_vegr,sg_road,sg_walle,sg_wallw,sg_vegc,sg_rfsn,sg_rdsn,  &
                  if_hwratio,if_sigmabld,rfveg%sigma,if_roof%alpha,rfveg%alpha,             &
                  cnveg%sigma,if_road%alpha,if_wall%alpha,cnveg%alpha,                      &
                  rfhyd%snowalpha,rdhyd%snowalpha,d_rfsndelta,d_rdsndelta)
sg_roof  = (1.-if_roof%alpha)*sg_roof*a_sg
sg_vegr  = (1.-rfveg%alpha)*sg_vegr*a_sg
sg_walle = (1.-if_wall%alpha)*sg_walle*a_sg
sg_wallw = (1.-if_wall%alpha)*sg_wallw*a_sg
sg_road  = (1.-if_road%alpha)*sg_road*a_sg
sg_vegc  = (1.-cnveg%alpha)*sg_vegc*a_sg
sg_rfsn  = (1.-rfhyd%snowalpha)*sg_rfsn*a_sg
sg_rdsn  = (1.-rdhyd%snowalpha)*sg_rdsn*a_sg

! calculate long wave reflections to nrefl order (pregenerated before canyonflux subroutine)
call getlwcoeff(d_netemiss,d_cwa,d_cra,d_cw0,d_cww,d_crw,d_crr,d_cwr,d_rdsndelta,wallpsi,roadpsi,cnveg%sigma,if_road%emiss,  &
                cnveg%emiss,if_wall%emiss)
p_emiss = d_rfsndelta*snowemiss+(1.-d_rfsndelta)*((1.-rfveg%sigma)*if_roof%emiss+rfveg%sigma*rfveg%emiss)
p_emiss = if_sigmabld*p_emiss+(1.-if_sigmabld)*(2.*if_wall%emiss*if_effhwratio*d_cwa+d_netemiss*d_cra) ! diagnostic only

! estimate bulk in-canyon surface roughness length
dis   = max(max(max(0.1*if_coeffbldheight*if_bldheight,zocanyon+0.2),cnveg%zo+0.2),zosnow+0.2)
zolog = 1./sqrt(d_rdsndelta/log(dis/zosnow)**2+(1.-d_rdsndelta)*(cnveg%sigma/log(dis/cnveg%zo)**2  &
       +(1.-cnveg%sigma)/log(dis/zocanyon)**2))
zonet = dis*exp(-zolog)

! estimate overall urban roughness length
zom = zomratio*if_bldheight
where ( zom*if_sigmabld<zonet*(1.-if_sigmabld) ) ! MJT suggestion
  zom = zonet
end where
n   = rdhyd%snow/(rdhyd%snow+maxrdsn+0.408*grav*zom)     ! snow cover for urban roughness calc (Douville, et al 1995)
zom = (1.-n)*zom + n*zosnow                            ! blend urban and snow roughness lengths (i.e., snow fills canyon)

! Calculate distance from atmosphere to displacement height
d_rfdzmin = max(a_zmin-if_bldheight,zoroof+0.2,rfveg%zo+0.2) ! distance to roof displacement height
p_cndzmin = max(a_zmin-refheight*if_bldheight,1.5,zom+0.2)   ! distance to canyon displacement height
p_lzom    = log(p_cndzmin/zom)

! calculate canyon wind speed and bulk transfer coefficents
! (i.e., acond = 1/(aerodynamic resistance) )
! some terms are updated when calculating canyon air temperature
select case(resmeth)
  case(0) ! Masson (2000)
    cu=exp(-0.25*if_effhwratio)
    abase_road =cu ! bulk transfer coefficents are updated in canyonflux
    abase_walle=cu
    abase_wallw=cu
    abase_rdsn =cu
    abase_vegc =cu
  case(1) ! Harman et al (2004)
    we=0. ! for cray compiler
    ww=0. ! for cray compiler
    wr=0. ! for cray compiler
    ! estimate wind speed along canyon surfaces
    call getincanwind(we,ww,wr,a_udir,zonet,if_bldheight,if_coeffbldheight,if_hwratio,ufull)
    dis=max(0.1*if_coeffbldheight*if_bldheight,zocanyon+0.2)
    zolog=log(dis/zocanyon)
    ! calculate terms for turbulent fluxes
    a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
    abase_walle=a*we                 ! east wall bulk transfer
    abase_wallw=a*ww                 ! west wall bulk transfer
    dis=max(0.1*if_coeffbldheight*if_bldheight,zocanyon+0.2,cnveg%zo+0.2,zosnow+0.2)
    zolog=log(dis/zocanyon)
    a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
    abase_road=a*wr                  ! road bulk transfer
    zolog=log(dis/cnveg%zo)
    a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
    abase_vegc=a*wr
    zolog=log(dis/zosnow)
    a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
    abase_rdsn=a*wr                  ! road snow bulk transfer
  case(2) ! Kusaka et al (2001)
    cu=exp(-0.386*if_effhwratio)
    abase_road =cu ! bulk transfer coefficents are updated in canyonflux
    abase_walle=cu
    abase_wallw=cu
    abase_rdsn =cu
    abase_vegc =cu
  case(3) ! Harman et al (2004)
    we=0. ! for cray compiler
    ww=0. ! for cray compiler
    wr=0. ! for cray compiler
    call getincanwindb(we,ww,wr,a_udir,zonet,if_bldheight,if_coeffbldheight,if_hwratio,ufull)
    dis=max(0.1*if_coeffbldheight*if_bldheight,zocanyon+0.2)
    zolog=log(dis/zocanyon)
    a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
    abase_walle=a*we                 ! east wall bulk transfer
    abase_wallw=a*ww                 ! west wall bulk transfer
    dis=max(0.1*if_coeffbldheight*if_bldheight,zocanyon+0.2,cnveg%zo+0.2,zosnow+0.2)
    zolog=log(dis/zocanyon)
    a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
    abase_road=a*wr                  ! road bulk transfer
    zolog=log(dis/cnveg%zo)
    a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
    abase_vegc=a*wr
    zolog=log(dis/zosnow)
    a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
    abase_rdsn=a*wr                  ! road snow bulk transfer
end select
  
! join two walls into a single wall (testing only)
if ( useonewall==1 ) then
  do k = 1,nl
    walle%nodetemp(:,k) = 0.5*(walle%nodetemp(:,k)+wallw%nodetemp(:,k))
    wallw%nodetemp(:,k) = walle%nodetemp(:,k)
  end do
  abase_walle = 0.5*(abase_walle+abase_wallw)
  abase_wallw = abase_walle
  sg_walle    = 0.5*(sg_walle+sg_wallw)
  sg_wallw    = sg_walle
end if

call getdiurnal(if_ctime,cyc_traffic,cyc_basedemand,cyc_proportion,cyc_translation)
! cyc_basedemand=1.
! cyc_proportion=1.
! cyc_translation=0.
! traffic sensible heat flux
p_traf = if_trafficfg*cyc_traffic
d_traf = p_traf/(1.-if_sigmabld)
! internal gains sensible heat flux
d_intgains_bld = (if_intmassn+1.)*if_intgains_flr*cyc_basedemand ! building internal gains 
p_intgains_full= if_sigmabld*d_intgains_bld                     ! full domain internal gains

! calculate canyon fluxes
call solvecanyon(sg_road,rg_road,fg_road,eg_road,acond_road,abase_road,                          &
                 sg_walle,rg_walle,fg_walle,acond_walle,abase_walle,                             &
                 sg_wallw,rg_wallw,fg_wallw,acond_wallw,abase_wallw,                             & 
                 sg_vegc,rg_vegc,fg_vegc,eg_vegc,acond_vegc,abase_vegc,                          &
                 sg_rdsn,rg_rdsn,fg_rdsn,eg_rdsn,acond_rdsn,abase_rdsn,rdsntemp,rdsnmelt,gardsn, &
                 a_umag,a_rho,a_rg,a_rnd,a_snd,                                                  &
                 d_canyontemp,d_canyonmix,d_tempc,d_mixrc,d_sigd,d_topu,d_netrad,                &
                 d_roaddelta,d_vegdeltac,d_rdsndelta,d_ac_outside,d_traf,d_ac_inside,            &
                 d_canyonrgout,d_tranc,d_evapc,d_cwa,d_cra,d_cw0,d_cww,d_crw,d_crr,              &
                 d_cwr,d_totdepth,d_c1c,d_intgains_bld,fgtop,egtop,int_infilflux,int_newairtemp, &
                 int_infilfg,ggint_roof,ggint_walle,ggint_wallw,ggint_road,ggint_slab,           &
                 ggint_intm1,ggint_intm2,cyc_translation,cyc_proportion,ddt,                     &
                 cnveg,if_ach,if_bldairtemp,if_bldheight,if_bldwidth,if_coeffbldheight,          &
                 if_effhwratio,if_hwratio,if_intm,if_intmassn,if_road,if_roof,if_sigmabld,       &
                 if_slab,if_tempcool,if_tempheat,if_wall,intm,p_bldcool,p_bldheat,p_cndzmin,     &
                 p_lzoh,p_lzom,rdhyd,rfveg,road,roof,room,slab,walle,wallw,if_sfc,if_swilt,ufull)

! calculate roof fluxes (fg_roof updated in solvetridiag)
eg_roof = 0. ! For cray compiler
call solveroof(sg_rfsn,rg_rfsn,fg_rfsn,eg_rfsn,garfsn,rfsnmelt,rfsntemp,acond_rfsn,d_rfsndelta, &
               sg_vegr,rg_vegr,fg_vegr,eg_vegr,acond_vegr,d_vegdeltar,                          &
               sg_roof,rg_roof,eg_roof,acond_roof,d_roofdelta,                                  &
               a_rg,a_umag,a_rho,a_rnd,a_snd,d_tempr,d_mixrr,d_rfdzmin,d_tranr,d_evapr,d_c1r,   &
               d_sigr,ddt,if_roof,rfhyd,rfveg,roof,if_rfvegdepth,if_sfc,if_swilt,ufull)

rgint_zero = 0.
! first internal temperature estimation - used for ggint calculation
select case(intairtmeth)
  case(0) ! fixed internal air temperature
    rgint_roof         = 0.
    rgint_walle        = 0.
    rgint_wallw        = 0.
    rgint_slab         = 0.
    
  case(1) ! floating internal air temperature
    call internal_lwflux(rgint_slab,rgint_wallw,rgint_roof,rgint_walle, &
                         if_bldheight,if_bldwidth,int_psi,int_viewf,roof,slab,walle,wallw,ufull)
                
  case DEFAULT
    write(6,*) "ERROR: Unknown intairtmeth mode ",intairtmeth
    stop
end select

! energy balance at facet surfaces
ggext_roof = (1.-d_rfsndelta)*(sg_roof+rg_roof-eg_roof+aircp*a_rho*d_tempr*acond_roof) &
              +d_rfsndelta*garfsn
ggext_walle= sg_walle+rg_walle+aircp*a_rho*d_canyontemp*acond_walle*if_coeffbldheight
ggext_wallw= sg_wallw+rg_wallw+aircp*a_rho*d_canyontemp*acond_wallw*if_coeffbldheight
ggext_road = (1.-d_rdsndelta)*(sg_road+rg_road-eg_road+aircp*a_rho*d_canyontemp*acond_road) &
             +d_rdsndelta*gardsn
             

! tridiagonal solver coefficents for calculating roof, road and wall temperatures
ggext_impl = (1.-d_rfsndelta)*aircp*a_rho*acond_roof  ! later update fg_roof with final roof skin T
depth_cp = if_roof%depth*if_roof%volcp
depth_lambda = if_roof%depth/if_roof%lambda
call solvetridiag(ggext_roof,ggint_roof,rgint_roof,ggext_impl,roof%nodetemp,ddt,     &
                  depth_cp, depth_lambda,ufull)
ggext_impl = aircp*a_rho*acond_walle*if_coeffbldheight ! later update fg_walle with final walle skin T
depth_cp = if_wall%depth*if_wall%volcp
depth_lambda = if_wall%depth/if_wall%lambda
call solvetridiag(ggext_walle,ggint_walle,rgint_walle,ggext_impl,walle%nodetemp,ddt,  &
                  depth_cp,depth_lambda,ufull)
ggext_impl = aircp*a_rho*acond_wallw*if_coeffbldheight ! later update fg_wallw with final wallw skin T
depth_cp = if_wall%depth*if_wall%volcp
depth_lambda = if_wall%depth/if_wall%lambda 
call solvetridiag(ggext_wallw,ggint_wallw,rgint_wallw,ggext_impl,wallw%nodetemp,ddt,  &
                  depth_cp,depth_lambda,ufull)
! rgint_road=0
ggext_impl = (1.-d_rdsndelta)*aircp*a_rho*acond_road ! later update fg_road with final road skin T
depth_cp = if_road%depth*if_road%volcp
depth_lambda = if_road%depth/if_road%lambda 
call solvetridiag(ggext_road,ggint_road,rgint_zero,ggext_impl,road%nodetemp,ddt,      &
                  depth_cp,depth_lambda,ufull)

! implicit update for fg to improve stability for thin layers
fg_roof = aircp*a_rho*(roof%nodetemp(:,0)-d_tempr)*acond_roof
fg_walle = aircp*a_rho*(walle%nodetemp(:,0)-d_canyontemp)*acond_walle*if_coeffbldheight
fg_wallw = aircp*a_rho*(wallw%nodetemp(:,0)-d_canyontemp)*acond_wallw*if_coeffbldheight
fg_road = aircp*a_rho*(road%nodetemp(:,0)-d_canyontemp)*acond_road

! update canyon flux
fgtop = if_hwratio*(fg_walle+fg_wallw) + (1.-d_rdsndelta)*(1.-cnveg%sigma)*fg_road &
      + (1.-d_rdsndelta)*cnveg%sigma*fg_vegc + d_rdsndelta*fg_rdsn                 &
      + d_traf + d_ac_outside - int_infilfg

! calculate internal facet conduction and temperature
ggext_impl = 0.
if ( intairtmeth==1 ) then
  depth_cp = if_slab%depth*if_slab%volcp
  depth_lambda = if_slab%depth/if_slab%lambda
  call solvetridiag(ggext_slab,ggint_slab,rgint_slab,ggext_impl,slab%nodetemp,ddt,     &
                    depth_cp,depth_lambda,ufull)
  if ( intmassmeth/=0 ) then
    ! rgint_intm=0
    ! negative ggint_intm1 (as both ggext and ggint are inside surfaces)
    depth_cp = if_intm%depth*if_intm%volcp
    depth_lambda = if_intm%depth/if_intm%lambda
    ggint_intm1_temp = -ggint_intm1
    call solvetridiag(ggint_intm1_temp,ggint_intm2,rgint_zero,ggext_impl,intm%nodetemp,ddt, &
                      depth_cp,depth_lambda,ufull)
  end if

  ! per m^2
  room%nodetemp(:,1) = room%nodetemp(:,1) + ddt/(a_rho*aircp*if_bldheight) *            & 
                  ((if_bldheight/if_bldwidth)*(ggint_walle + ggint_wallw)               &
                  + ggint_roof + ggint_slab + if_intmassn*(ggint_intm2 + ggint_intm1)   &
                  + int_infilflux + d_ac_inside + d_intgains_bld)
end if

! calculate water/snow budgets for road surface
call updatewater(ddt,rdhyd%surfwater,rdhyd%soilwater,rdhyd%leafwater,rdhyd%snow,    &
                     rdhyd%den,rdhyd%snowalpha,rdsnmelt,a_rnd,a_snd,eg_road,        &
                     eg_rdsn,d_tranc,d_evapc,d_c1c,d_totdepth, cnveg%lai,wbrelaxc,  &
                     if_sfc,if_swilt,ufull)

! calculate water/snow budgets for roof surface
call updatewater(ddt,rfhyd%surfwater,rfhyd%soilwater,rfhyd%leafwater,rfhyd%snow,     &
                     rfhyd%den,rfhyd%snowalpha,rfsnmelt,a_rnd,a_snd,eg_roof,         &
                     eg_rfsn,d_tranr,d_evapr,d_c1r,if_rfvegdepth,rfveg%lai,wbrelaxr, &
                     if_sfc,if_swilt,ufull)

! calculate runoff (leafwater runoff already accounted for in precip reaching canyon floor)
u_rn = max(rfhyd%surfwater-maxrfwater,0.)*if_sigmabld*(1.-rfveg%sigma)                   &
      +max(rdhyd%surfwater-maxrdwater,0.)*(1.-if_sigmabld)*(1.-cnveg%sigma)              &
      +max(rfhyd%snow-maxrfsn,0.)*if_sigmabld                                            &
      +max(rdhyd%snow-maxrdsn,0.)*(1.-if_sigmabld)                                       &
      +max(rfhyd%soilwater-if_ssat,0.)*waterden*if_rfvegdepth*rfveg%sigma*if_sigmabld    &
      +max(rdhyd%soilwater-if_ssat,0.)*waterden*d_totdepth*cnveg%sigma*(1.-if_sigmabld)

! remove round-off problems
rdhyd%soilwater(1:ufull) = min(max(rdhyd%soilwater(1:ufull),if_swilt),if_ssat)
rfhyd%soilwater(1:ufull) = min(max(rfhyd%soilwater(1:ufull),if_swilt),if_ssat)
rfhyd%surfwater(1:ufull) = min(max(rfhyd%surfwater(1:ufull),0.),maxrfwater)
rdhyd%surfwater(1:ufull) = min(max(rdhyd%surfwater(1:ufull),0.),maxrdwater)
rdhyd%leafwater(1:ufull) = min(max(rdhyd%leafwater(1:ufull),0.),maxvwatf*cnveg%lai)
rfhyd%leafwater(1:ufull) = min(max(rfhyd%leafwater(1:ufull),0.),maxvwatf*rfveg%lai)
rfhyd%snow(1:ufull)      = min(max(rfhyd%snow(1:ufull),0.),maxrfsn)
rdhyd%snow(1:ufull)      = min(max(rdhyd%snow(1:ufull),0.),maxrdsn)
rfhyd%den(1:ufull)       = min(max(rfhyd%den(1:ufull),minsnowden),maxsnowden)
rdhyd%den(1:ufull)       = min(max(rdhyd%den(1:ufull),minsnowden),maxsnowden)
rfhyd%snowalpha(1:ufull) = min(max(rfhyd%snowalpha(1:ufull),minsnowalpha),maxsnowalpha)
rdhyd%snowalpha(1:ufull) = min(max(rdhyd%snowalpha(1:ufull),minsnowalpha),maxsnowalpha)

! combine snow and snow-free tiles for fluxes
d_roofrgout = a_rg-d_rfsndelta*rg_rfsn-(1.-d_rfsndelta)*((1.-rfveg%sigma)*rg_roof+rfveg%sigma*rg_vegr)
fgrooftop   = d_rfsndelta*fg_rfsn+(1.-d_rfsndelta)*((1.-rfveg%sigma)*fg_roof+rfveg%sigma*fg_vegr)
egrooftop   = d_rfsndelta*eg_rfsn+(1.-d_rfsndelta)*((1.-rfveg%sigma)*eg_roof+rfveg%sigma*eg_vegr)
!fgtop       = d_rdsndelta*fg_rdsn+(1.-d_rdsndelta)*((1.-cnveg%sigma)*fg_road+cnveg%sigma*fg_vegc)   &
!             +if_hwratio*(fg_walle+fg_wallw)+d_traf+d_ac_outside
!egtop       = d_rdsndelta*eg_rdsn+(1.-d_rdsndelta)*((1.-cnveg%sigma)*eg_road+cnveg%sigma*eg_vegc)

! calculate wetfac for roof and road vegetation (see sflux.f or cable_canopy.f90)
roofvegwetfac = max(min((rfhyd%soilwater-if_swilt)/(if_sfc-if_swilt),1.),0.)
roadvegwetfac = max(min((rdhyd%soilwater-if_swilt)/(if_sfc-if_swilt),1.),0.)

! calculate longwave, sensible heat latent heat outputs
! estimate surface temp from outgoing longwave radiation
u_ts = ((if_sigmabld*d_roofrgout+(1.-if_sigmabld)*d_canyonrgout)/sbconst)**0.25 - urbtemp
u_fg = if_sigmabld*fgrooftop+(1.-if_sigmabld)*fgtop+if_industryfg
u_eg = if_sigmabld*egrooftop+(1.-if_sigmabld)*egtop
u_wf = if_sigmabld*(1.-d_rfsndelta)*((1.-rfveg%sigma)*d_roofdelta       &
      +rfveg%sigma*((1.-d_vegdeltar)*roofvegwetfac+d_vegdeltar))       &
      +(1.-if_sigmabld)*(1.-d_rdsndelta)*((1.-cnveg%sigma)*d_roaddelta  &
      +cnveg%sigma*((1.-d_vegdeltac)*roadvegwetfac+d_vegdeltac))

p_snowmelt = if_sigmabld*rfsnmelt + (1.-if_sigmabld)*rdsnmelt
u_melt = lf*(if_sigmabld*d_rfsndelta*rfsnmelt + (1.-if_sigmabld)*d_rdsndelta*rdsnmelt)

! (re)calculate heat roughness length for MOST (diagnostic only)
call getqsat(a,u_ts,d_sigd)
dts = u_ts + (u_ts+urbtemp)*0.61*a*u_wf
dtt = d_tempc + (d_tempc+urbtemp)*0.61*d_mixrc
select case(zohmeth)
  case(0) ! Use veg formulation
    p_lzoh = 2.3+p_lzom
    call getinvres(p_cdtq,p_cduv,z_on_l,p_lzoh,p_lzom,p_cndzmin,dts,dtt,a_umag,1)
  case(1) ! Use Kanda parameterisation
    p_lzoh = 2.3+p_lzom ! replaced in getlna
    call getinvres(p_cdtq,p_cduv,z_on_l,p_lzoh,p_lzom,p_cndzmin,dts,dtt,a_umag,2)
  case(2) ! Use Kanda parameterisation
    p_lzoh = 6.+p_lzom
    call getinvres(p_cdtq,p_cduv,z_on_l,p_lzoh,p_lzom,p_cndzmin,dts,dtt,a_umag,4)
end select

! calculate screen level diagnostics
call scrncalc(a_mixr,a_umag,a_temp,u_ts,d_tempc,d_rdsndelta,d_roaddelta,d_vegdeltac,d_sigd,a,rdsntemp,zonet, &
              cnveg,if_bldheight,if_sfc,if_swilt,p_cndzmin,p_lzoh,p_lzom,p_qscrn,p_tscrn,p_u10,p_uscrn,rdhyd,road,ufull)

call energyclosure(sg_roof,rg_roof,fg_roof,sg_walle,rg_walle,fg_walle,     &
                   sg_road,rg_road,fg_road,sg_wallw,rg_wallw,fg_wallw,     &
                   rgint_roof,rgint_walle,rgint_wallw,rgint_slab,          &
                   eg_roof,eg_road,garfsn,gardsn,d_rfsndelta,d_rdsndelta,  &
                   a_sg,a_rg,u_ts,u_fg,u_eg,u_alb,u_melt,a_rho,            &
                   ggint_roof,ggint_road,ggint_walle,ggint_wallw,          &
                   ggint_intm1,ggint_slab,ggint_intm2,d_intgains_bld,      &
                   int_infilflux,d_ac_inside,if_bldwidth,ddt,              &
                   cnveg,if_bldheight,if_hwratio,if_industryfg,if_intm,    &
                   if_intmassn,if_road,if_roof,if_sigmabld,if_slab,        &
                   if_wall,intm,p_atmoserr,p_bldcool,p_bldheat,            &
                   p_intgains_full,p_surferr,p_traf,rfveg,road,roof,room,  &
                   slab,walle,wallw,ufull)

return
end subroutine atebeval

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculates flux from facet into room using estimates of next timestep temperatures
! This taylor expansion is necassary for stability where intairtmeth=0 (no internal model)
! and with internal model when internal wall/roof layer has low heat capacity (insulation)
! May be depreciated in future.

subroutine calc_ggint(depth,volcp,lambda,skintemp,newairtemp,cvcoeff,ddt,ggint,ufull)

implicit none

integer, intent(in) :: ufull
real, intent(in)                    :: ddt
real, intent(in), dimension(ufull)  :: depth,volcp,lambda,cvcoeff
real, intent(in), dimension(ufull)  :: skintemp, newairtemp
real, intent(out), dimension(ufull) :: ggint
real, dimension(ufull) :: condterm, newskintemp

select case(conductmeth)
  case(0) ! half-layer conduction
    condterm = 1./(0.5*depth/lambda +1./cvcoeff)
    newskintemp  = skintemp-condterm*(skintemp-newairtemp) &
                    /(volcp*depth/ddt+condterm)

  case(1) ! interface conduction
    condterm = cvcoeff
    newskintemp  = skintemp-condterm*(skintemp-newairtemp) &
                /(0.5*volcp*depth/ddt+condterm)
end select

ggint = condterm*(newskintemp-newairtemp)

end subroutine calc_ggint

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Tridiagonal solver for temperatures

! This version has an implicit estimate for roof sensible heat flux

! [ ggB ggC         ] [ temp ] = [ ggD ]
! [ ggA ggB ggC     ] [ temp ] = [ ggD ]
! [     ggA ggB ggC ] [ temp ] = [ ggD ]
! [         ggA ggB ] [ temp ] = [ ggD ]

subroutine solvetridiag(ggext,ggint,rgint,ggimpl,nodetemp,ddt,cap,res,ufull)

implicit none

integer, intent(in) :: ufull
real, dimension(ufull),     intent(in)    :: ggext,ggint,rgint  ! surface energy fluxes
real, dimension(ufull),     intent(in)    :: ggimpl             ! implicit update for roof only
real, dimension(ufull,0:nl),intent(inout) :: nodetemp           ! temperature of each node
real, dimension(ufull,nl),  intent(in)    :: cap,res            ! layer capacitance & resistance
real, dimension(ufull,0:nl)               :: ggA,ggB,ggC,ggD    ! tridiagonal matrices
real, dimension(ufull)                    :: ggX                ! tridiagonal coefficient
real, intent(in)                          :: ddt                ! timestep
integer k

select case(conductmeth)
  case(0) !!!!!!!!! half-layer conduction !!!!!!!!!!!
    ggA(:,1)      =-2./res(:,1)
    ggA(:,2:nl)   =-2./(res(:,1:nl-1) +res(:,2:nl))
    ggB(:,0)      = 2./res(:,1) + ggimpl
    ggB(:,1)      = 2./res(:,1) +2./(res(:,1)+res(:,2)) + cap(:,1)/ddt
    ggB(:,2:nl-1) = 2./(res(:,1:nl-2) +res(:,2:nl-1)) +2./(res(:,2:nl-1) +res(:,3:nl)) +cap(:,2:nl-1)/ddt
    ggB(:,nl)     = 2./(res(:,nl-1)+res(:,nl)) + cap(:,nl)/ddt
    ggC(:,0)      =-2./res(:,1)
    ggC(:,1:nl-1) =-2./(res(:,1:nl-1)+res(:,2:nl))
    ggD(:,0)      = ggext
    ggD(:,1:nl-1) = nodetemp(:,1:nl-1)*cap(:,1:nl-1)/ddt
    ggD(:,nl)     = nodetemp(:,nl)*cap(:,nl)/ddt - ggint - rgint
  case(1) !!!!!!!!! interface conduction !!!!!!!!!!!
    ggA(:,1:nl)   = -1./res(:,1:nl)
    ggB(:,0)      =  1./res(:,1) +0.5*cap(:,1)/ddt + ggimpl
    ggB(:,1:nl-1) =  1./res(:,1:nl-1) +1./res(:,2:nl) +0.5*(cap(:,1:nl-1) +cap(:,2:nl))/ddt
    ggB(:,nl)     =  1./res(:,nl) + 0.5*cap(:,nl)/ddt
    ggC(:,0:nl-1) = -1./res(:,1:nl)
    ggD(:,0)      = nodetemp(:,0)*0.5*cap(:,1)/ddt + ggext
    ggD(:,1:nl-1) = nodetemp(:,1:nl-1)*0.5*(cap(:,1:nl-1)+cap(:,2:nl))/ddt
    ggD(:,nl)     = nodetemp(:,nl)*0.5*cap(:,nl)/ddt - ggint - rgint
end select
! tridiagonal solver (Thomas algorithm) to solve node temperatures
do k=1,nl
  ggX(:)   = ggA(:,k)/ggB(:,k-1)
  ggB(:,k) = ggB(:,k)-ggX(:)*ggC(:,k-1)
  ggD(:,k) = ggD(:,k)-ggX(:)*ggD(:,k-1)
end do
nodetemp(:,nl) = ggD(:,nl)/ggB(:,nl)
do k=nl-1,0,-1
  nodetemp(:,k) = (ggD(:,k) - ggC(:,k)*nodetemp(:,k+1))/ggB(:,k)
end do

end subroutine solvetridiag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Conservation of energy check

subroutine energyclosure(sg_roof,rg_roof,fg_roof,sg_walle,rg_walle,fg_walle,     &
                         sg_road,rg_road,fg_road,sg_wallw,rg_wallw,fg_wallw,     &
                         rgint_roof,rgint_walle,rgint_wallw,rgint_slab,          &
                         eg_roof,eg_road,garfsn,gardsn,d_rfsndelta,d_rdsndelta,  &
                         a_sg,a_rg,u_ts,u_fg,u_eg,u_alb,u_melt,a_rho,            &
                         ggint_roof,ggint_road,ggint_walle,ggint_wallw,          &
                         ggint_intm1,ggint_slab,ggint_intm2,d_intgains_bld,      &
                         int_infilflux,d_ac_inside,if_bldwidth,ddt,              &
                         cnveg,if_bldheight,if_hwratio,if_industryfg,if_intm,    &
                         if_intmassn,if_road,if_roof,if_sigmabld,if_slab,        &
                         if_wall,intm,p_atmoserr,p_bldcool,p_bldheat,            &
                         p_intgains_full,p_surferr,p_traf,rfveg,road,roof,room,  &
                         slab,walle,wallw,ufull)

implicit none

integer, intent(in) :: ufull
real, intent(in) :: ddt
real, dimension(ufull), intent(in) :: sg_roof,rg_roof,fg_roof,sg_walle,rg_walle,fg_walle
real, dimension(ufull), intent(in) :: sg_road,rg_road,fg_road,sg_wallw,rg_wallw,fg_wallw
real, dimension(ufull), intent(in) :: rgint_roof,rgint_walle,rgint_wallw,rgint_slab
real, dimension(ufull), intent(in) :: eg_roof,eg_road,garfsn,gardsn,d_rfsndelta,d_rdsndelta
real, dimension(ufull), intent(in) :: a_sg,a_rg,u_ts,u_fg,u_eg,u_alb,u_melt,a_rho
real, dimension(ufull), intent(in) :: ggint_roof,ggint_road,ggint_walle,ggint_wallw
real, dimension(ufull), intent(in) :: ggint_intm1,ggint_slab,ggint_intm2,d_intgains_bld
real, dimension(ufull), intent(in) :: int_infilflux,d_ac_inside,if_bldwidth
real(kind=8), dimension(ufull) :: d_roofflux,d_walleflux,d_wallwflux,d_roadflux,d_slabflux,d_intmflux,d_roomflux 
real(kind=8), dimension(ufull) :: d_roofstor,d_wallestor,d_wallwstor,d_roadstor,d_slabstor,d_intmstor,d_roomstor
real(kind=8), dimension(ufull) :: d_faceterr
real(kind=8), dimension(ufull) :: d_storageflux,d_atmosflux
real(kind=8), dimension(ufull,nl) :: roadstorage_prev, roofstorage_prev, wallestorage_prev, wallwstorage_prev
real(kind=8), dimension(ufull,nl) :: slabstorage_prev, intmstorage_prev
real(kind=8), dimension(ufull,1) :: roomstorage_prev
!global
real, dimension(ufull), intent(in) :: if_bldheight, if_hwratio, if_industryfg, if_sigmabld
integer, dimension(ufull), intent(in) :: if_intmassn
type(facetparams), intent(in) :: if_intm, if_road, if_roof, if_slab, if_wall
type(facetdata), intent(inout) :: intm
real, dimension(ufull), intent(in) :: p_bldcool, p_bldheat, p_intgains_full, p_traf
real(kind=8), dimension(ufull), intent(inout) :: p_atmoserr, p_surferr
type(vegdata), intent(in) :: cnveg, rfveg
type(facetdata), intent(inout) :: road, roof, room, slab, walle, wallw
!

! Store previous calculation to determine flux
roofstorage_prev(:,:)  = roof%storage(:,:)
roadstorage_prev(:,:)  = road%storage(:,:)
wallestorage_prev(:,:) = walle%storage(:,:)
wallwstorage_prev(:,:) = wallw%storage(:,:)
slabstorage_prev(:,:)  = slab%storage(:,:)
intmstorage_prev(:,:)  = intm%storage(:,:)
roomstorage_prev(:,:)  = room%storage(:,:)
p_surferr = 0.


room%storage(:,1) = real(if_bldheight(:),8)*real(a_rho(:),8)*real(aircp,8)*real(room%nodetemp(:,1),8)
! Sum heat stored in urban materials from layer 1 to nl
select case(conductmeth)
  case(0) ! half-layer conduction
    roof%storage(:,:) = real(if_roof%depth(:,:),8)*real(if_roof%volcp(:,:),8)*real(roof%nodetemp(:,1:nl),8)
    road%storage(:,:) = real(if_road%depth(:,:),8)*real(if_road%volcp(:,:),8)*real(road%nodetemp(:,1:nl),8)
    walle%storage(:,:)= real(if_wall%depth(:,:),8)*real(if_wall%volcp(:,:),8)*real(walle%nodetemp(:,1:nl),8)
    wallw%storage(:,:)= real(if_wall%depth(:,:),8)*real(if_wall%volcp(:,:),8)*real(wallw%nodetemp(:,1:nl),8)
    slab%storage(:,:) = real(if_slab%depth(:,:),8)*real(if_slab%volcp(:,:),8)*real(slab%nodetemp(:,1:nl),8)
    intm%storage(:,:) = real(if_intm%depth(:,:),8)*real(if_intm%volcp(:,:),8)*real(intm%nodetemp(:,1:nl),8)
  case(1) ! interface conduction
    roof%storage(:,:)  = 0.5_8*real(if_roof%depth(:,:),8)*real(if_roof%volcp(:,:),8)                        & 
                            *(real(roof%nodetemp(:,0:nl-1),8)+real(roof%nodetemp(:,1:nl),8))
    road%storage(:,:)  = 0.5_8*real(if_road%depth(:,:),8)*real(if_road%volcp(:,:),8)                        & 
                            *(real(road%nodetemp(:,0:nl-1),8)+real(road%nodetemp(:,1:nl),8))
    walle%storage(:,:) = 0.5_8*real(if_wall%depth(:,:),8)*real(if_wall%volcp(:,:),8)                        & 
                            *(real(walle%nodetemp(:,0:nl-1),8)+real(walle%nodetemp(:,1:nl),8))
    wallw%storage(:,:) = 0.5_8*real(if_wall%depth(:,:),8)*real(if_wall%volcp(:,:),8)                        & 
                            *(real(wallw%nodetemp(:,0:nl-1),8)+real(wallw%nodetemp(:,1:nl),8))
    slab%storage(:,:)  = 0.5_8*real(if_slab%depth(:,:),8)*real(if_slab%volcp(:,:),8)                        & 
                            *(real(slab%nodetemp(:,0:nl-1),8)+real(slab%nodetemp(:,1:nl),8))
    intm%storage(:,:)  = 0.5_8*real(if_intm%depth(:,:),8)*real(if_intm%volcp(:,:),8)                        & 
                            *(real(intm%nodetemp(:,0:nl-1),8)+real(intm%nodetemp(:,1:nl),8))
end select

if ( all(roofstorage_prev==0._8) ) return
  
d_roofstor = sum(roof%storage-roofstorage_prev,dim=2)/real(ddt,8)
d_roofflux = (1._8-real(d_rfsndelta,8))*(real(sg_roof,8)+real(rg_roof,8)-real(fg_roof,8)-real(eg_roof,8))  &
           + real(d_rfsndelta,8)*real(garfsn,8) - real(ggint_roof,8) - real(rgint_roof,8)
d_faceterr  = d_roofstor - d_roofflux
p_surferr = p_surferr + d_faceterr
if (any(abs(d_faceterr)>=energytol)) write(6,*) "aTEB roof facet closure error:", maxval(abs(d_faceterr))
d_roadstor = sum(road%storage-roadstorage_prev,dim=2)/real(ddt,8)
d_roadflux = (1._8-real(d_rdsndelta,8))*(real(sg_road,8)+real(rg_road,8)-real(fg_road,8)-real(eg_road,8)) &
           + real(d_rdsndelta,8)*real(gardsn,8) - real(ggint_road,8)
d_faceterr  = d_roadstor - d_roadflux
p_surferr = p_surferr + d_faceterr
if (any(abs(d_faceterr)>=energytol)) write(6,*) "aTEB road facet closure error:", maxval(abs(d_faceterr))
d_wallestor= sum(walle%storage-wallestorage_prev,dim=2)/real(ddt,8)
d_walleflux= real(sg_walle,8)+real(rg_walle,8)-real(fg_walle,8) - real(ggint_walle,8) - real(rgint_walle,8)
d_faceterr = d_wallestor - d_walleflux
p_surferr = p_surferr + d_faceterr
if (any(abs(d_faceterr)>=energytol)) write(6,*) "aTEB walle facet closure error:", maxval(abs(d_faceterr))
d_wallwstor= sum(wallw%storage-wallwstorage_prev,dim=2)/real(ddt,8)
d_wallwflux= real(sg_wallw,8)+real(rg_wallw,8)-real(fg_wallw,8) - real(ggint_wallw,8) - real(rgint_wallw,8)
d_faceterr = d_wallwstor - d_wallwflux
p_surferr = p_surferr + d_faceterr
if (any(abs(d_faceterr)>=energytol)) write(6,*) "aTEB wallw facet closure error:", maxval(abs(d_faceterr))
if (intairtmeth==1) then
  d_slabstor = sum(slab%storage-slabstorage_prev,dim=2)/real(ddt,8)
  d_slabflux = -real(ggint_slab,8) - real(rgint_slab,8)
  d_faceterr = d_slabstor - d_slabflux
  p_surferr = p_surferr + d_faceterr
  if (any(abs(d_faceterr)>=energytol)) write(6,*) "aTEB slab facet closure error:", maxval(abs(d_faceterr))
  d_intmstor = sum(intm%storage-intmstorage_prev,dim=2)/real(ddt,8)
  d_intmflux = -real(ggint_intm1,8) - real(ggint_intm2,8)
  d_faceterr = d_intmstor - d_intmflux
  p_surferr = p_surferr + d_faceterr
  if (any(abs(d_faceterr)>=energytol)) write(6,*) "aTEB intm facet closure error:", maxval(abs(d_faceterr))
  d_roomstor = (room%storage(:,1)-roomstorage_prev(:,1))/real(ddt,8)
  d_roomflux = real(ggint_roof,8)+real(ggint_slab,8)-real(if_intmassn,8)*real(d_intmflux,8)            & 
            + (real(if_bldheight,8)/real(if_bldwidth,8))*(real(ggint_walle,8) + real(ggint_wallw,8))   &
            + real(int_infilflux,8) + real(d_ac_inside,8) + real(d_intgains_bld,8)
  d_faceterr = d_roomstor - d_roomflux
  p_surferr = p_surferr + d_faceterr
  if (any(abs(d_faceterr)>=energytol)) write(6,*) "aTEB room volume closure error:", maxval(abs(d_faceterr))
else
  d_slabstor = 0._8
  d_intmstor = 0._8
  d_roomstor = 0._8
end if

d_storageflux = d_roofstor*real(if_sigmabld,8)*(1._8-real(rfveg%sigma,8))           &
              + d_roadstor*(1._8-real(if_sigmabld,8))*(1._8-real(cnveg%sigma,8))    &
              + d_wallestor*(1._8-real(if_sigmabld,8))*real(if_hwratio,8)           &
              + d_wallwstor*(1._8-real(if_sigmabld,8))*real(if_hwratio,8)           &
              + d_slabstor*real(if_sigmabld,8)                                      &
              + d_intmstor*real(if_sigmabld,8)*real(if_intmassn,8)                  &
              + d_roomstor*real(if_sigmabld,8)

! print *, 'd_storageflux',d_storageflux
! print *, 'roof  Qs' ,real(d_roofstor,8)*real(if_sigmabld,8)*(1-real(rfveg%sigma,8))      
! print *, 'road  Qs' ,real(d_roadstor,8)*(1-real(if_sigmabld,8))*(1-real(cnveg%sigma,8))  
! print *, 'walle Qs' ,real(d_wallestor,8)*(1-real(if_sigmabld,8))*real(if_hwratio,8)       
! print *, 'wallw Qs' ,real(d_wallwstor,8)*(1-real(if_sigmabld,8))*real(if_hwratio,8)       
! print *, 'slab  Qs' ,real(d_slabstor,8)*real(if_sigmabld,8)                              
! print *, 'intm  Qs' ,real(d_intmstor,8)*real(if_sigmabld,8)*real(if_intmassn,8)           
! print *, 'room  Qs' ,real(d_roomstor,8)*real(if_sigmabld,8)
! print *, 'room/slab', (d_roomstor*if_sigmabld)/(d_slabstor*if_sigmabld)
! print *, 'infil', real(int_infilflux,8)*real(if_sigmabld,8)

! atmosphere energy flux = (SWdown-SWup) + (LWdown-LWup) - Turbulent + Anthropogenic
d_atmosflux = (real(a_sg,8)-real(a_sg,8)*real(u_alb,8)) + (real(a_rg,8)-real(sbconst,8)*(real(u_ts,8)+urbtemp)**4) &
            - (real(u_fg,8)+real(u_eg,8)+real(u_melt,8)) + real(p_bldheat,8) + real(p_bldcool,8) + real(p_traf,8)  & 
            + real(if_industryfg,8) + real(p_intgains_full,8)
p_atmoserr = d_storageflux - d_atmosflux

if ( any(abs(p_atmoserr)>=energytol) ) then
  write(6,*) "aTEB energy not conserved! Atmos. error:", maxval(abs(p_atmoserr))
end if
! print *, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

return
end subroutine energyclosure

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update water prognostic variables for roads and roofs
                            
subroutine updatewater(ddt,surfwater,soilwater,leafwater,snow,den,alpha, &
                       snmelt,a_rnd,a_snd,eg_surf,eg_snow,d_tran,d_evap, &
                       d_c1,d_totdepth,if_vegrlai,iwbrelax,              &
                       if_sfc,if_swilt,ufull)

implicit none

integer, intent(in) :: ufull
integer, intent(in) :: iwbrelax
real, intent(in) :: ddt
real, dimension(ufull), intent(inout) :: surfwater,soilwater,leafwater,snow,den,alpha
real, dimension(ufull), intent(in) :: snmelt,a_rnd,a_snd,eg_surf,eg_snow
real, dimension(ufull), intent(in) :: d_tran,d_evap,d_c1,d_totdepth,if_vegrlai
real, dimension(ufull) :: modrnd
real, dimension(ufull), intent(in) :: if_sfc, if_swilt

modrnd = max(a_rnd-d_evap/lv-max(maxvwatf*if_vegrlai-leafwater,0.)/ddt,0.) ! rainfall reaching the soil under vegetation

! note that since sigmaf=1, then there is no soil evaporation, only transpiration.  Evaporation only occurs from water on leafs.
surfwater = surfwater+ddt*(a_rnd-eg_surf/lv+snmelt)                                         ! surface
soilwater = soilwater+ddt*d_c1*(modrnd+snmelt*den/waterden-d_tran/lv)/(waterden*d_totdepth) ! soil
leafwater = leafwater+ddt*(a_rnd-d_evap/lv)                                                 ! leaf
leafwater = min(max(leafwater,0.),maxvwatf*if_vegrlai)

if (iwbrelax==1) then
  ! increase soil moisture for irrigation 
  soilwater=soilwater+max(0.75*if_swilt+0.25*if_sfc-soilwater,0.)/(86400./ddt+1.) ! 24h e-fold time
end if

! snow fields
snow  = snow + ddt*(a_snd-eg_snow/lv-snmelt)
den   = den + (maxsnowden-den)/(0.24/(86400.*ddt)+1.)
alpha = alpha + (minsnowalpha-alpha)/(0.24/(86400.*ddt)+1.)

return
end subroutine updatewater

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Estimate saturation mixing ratio

subroutine getqsat(qsat,temp,ps)

implicit none

real, dimension(:), intent(in) :: temp
real, dimension(size(temp)), intent(in) :: ps
real, dimension(size(temp)), intent(out) :: qsat
real, dimension(size(temp)) :: esatf,tdiff,rx
integer, dimension(size(temp)) :: ix

tdiff=min(max( temp+(urbtemp-123.16), 0.), 219.)
rx=tdiff-aint(tdiff)
ix=int(tdiff)
esatf=(1.-rx)*table(ix)+ rx*table(ix+1)
qsat=0.622*esatf/max(ps-esatf,0.1)

return
end subroutine getqsat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Interface for calcuating ustar and thetastar

subroutine getinvres(invres,cd,z_on_l,olzoh,ilzom,zmin,sthetav,thetav,a_umag,mode)

implicit none

integer, intent(in) :: mode
real, dimension(:), intent(in) :: ilzom
real, dimension(size(ilzom)), intent(in) :: zmin,sthetav,thetav
real, dimension(size(ilzom)), intent(in) :: a_umag
real, dimension(size(ilzom)), intent(out) :: invres,cd,z_on_l
real, dimension(size(ilzom)), intent(inout) :: olzoh
real, dimension(size(ilzom)) :: lna,thetavstar,integralh

lna=olzoh-ilzom
call dyerhicks(integralh,z_on_l,cd,thetavstar,thetav,sthetav,a_umag,zmin,ilzom,lna,mode)
invres=vkar*sqrt(cd)*a_umag/integralh
olzoh=lna+ilzom

return
end subroutine getinvres

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate stability functions using Dyerhicks

subroutine dyerhicks(integralh,z_on_l,cd,thetavstar,thetav,sthetav,umagin,zmin,ilzom,lna,mode)

implicit none

integer, intent(in) :: mode
integer ic
real, dimension(:), intent(in) :: thetav
real, dimension(size(thetav)), intent(in) :: sthetav,umagin,zmin,ilzom
real, dimension(size(thetav)), intent(inout) :: lna
real, dimension(size(thetav)), intent(out) :: cd,thetavstar
real, dimension(size(thetav)), intent(out) :: integralh,z_on_l
real, dimension(size(thetav)) :: z0_on_l,zt_on_l,olzoh,umag
real, dimension(size(thetav)) :: pm0,ph0,pm1,ph1,integralm
!real, parameter :: aa1 = 3.8
!real, parameter :: bb1 = 0.5
!real, parameter :: cc1 = 0.3

umag = max(umagin, 0.01)
cd=(vkar/ilzom)**2                         ! first guess
call getlna(lna,cd,umag,zmin,ilzom,mode)
olzoh=ilzom+lna
integralh=sqrt(cd)*ilzom*olzoh/vkar        ! first guess
thetavstar=vkar*(thetav-sthetav)/integralh ! first guess

do ic=1,icmax
  z_on_l=vkar*zmin*grav*thetavstar/((thetav+urbtemp)*cd*umag**2)
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
  integralm = max( integralm, 1.e-10 )
  integralh = max( integralh, 1.e-10 )
  !where (z_on_l<=0.4)
    cd = (max(0.01,min(vkar*umag/integralm,2.))/umag)**2
  !elsewhere
  !  cd = (max(0.01,min(vkar*umag/(aa1*( ( z_on_l**bb1)*(1.0+cc1* z_on_l**(1.-bb1)) &
  !      -(z0_on_l**bb1)*(1.+cc1*z0_on_l**(1.-bb1)) )),2.))/umag)**2
  !endwhere
  thetavstar= vkar*(thetav-sthetav)/integralh
  call getlna(lna,cd,umag,zmin,ilzom,mode)
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

subroutine getswcoeff(sg_roof,sg_vegr,sg_road,sg_walle,sg_wallw,sg_vegc,sg_rfsn,sg_rdsn,wallpsi,roadpsi,if_hwratio, &
                      if_vangle,if_hangle,if_fbeam,if_vegsigmac,if_roadalpha,if_vegalphac,if_wallalpha,ird_alpha,   &
                      rdsndelta)

implicit none

integer k
real, dimension(:), intent(in) :: rdsndelta
real, dimension(size(rdsndelta)), intent(in) :: ird_alpha
real, dimension(size(rdsndelta)), intent(out) :: wallpsi,roadpsi
real, dimension(size(rdsndelta)), intent(in) :: if_hwratio
real, dimension(size(rdsndelta)), intent(in) :: if_vangle,if_hangle,if_fbeam,if_vegsigmac,if_roadalpha,if_vegalphac
real, dimension(size(rdsndelta)), intent(in) :: if_wallalpha
real, dimension(size(rdsndelta)), intent(out) :: sg_roof,sg_vegr,sg_road,sg_walle,sg_wallw,sg_vegc,sg_rfsn,sg_rdsn
real, dimension(size(rdsndelta)) :: thetazero,walles,wallws,roads,ta,tc,xa,ya,roadnetalpha
real, dimension(size(rdsndelta)) :: nwalles,nwallws,nroads

wallpsi=0.5*(if_hwratio+1.-sqrt(if_hwratio*if_hwratio+1.))/if_hwratio
roadpsi=sqrt(if_hwratio*if_hwratio+1.)-if_hwratio

! integrate through 180 deg instead of 360 deg.  Hence paritioning to east and west facing walls
where (if_vangle>=0.5*pi)
  walles=0.
  wallws=1./if_hwratio
  roads=0.
elsewhere
  ta=tan(if_vangle)
  thetazero=asin(1./max(if_hwratio*ta,1.))
  tc=2.*(1.-cos(thetazero))
  xa=min(max(if_hangle-thetazero,0.),pi)-max(if_hangle-pi+thetazero,0.)-min(if_hangle+thetazero,0.)
  ya=cos(max(min(0.,if_hangle),if_hangle-pi))-cos(max(min(thetazero,if_hangle),if_hangle-pi)) &
    +cos(min(0.,-if_hangle))-cos(min(thetazero,-if_hangle)) &
    +cos(max(0.,pi-if_hangle))-cos(max(thetazero,pi-if_hangle))
  ! note that these terms now include the azimuth angle
  walles=if_fbeam*(xa/if_hwratio+ta*ya)/pi+(1.-if_fbeam)*wallpsi
  wallws=if_fbeam*((pi-2.*thetazero-xa)/if_hwratio+ta*(tc-ya))/pi+(1.-if_fbeam)*wallpsi
  roads=if_fbeam*(2.*thetazero-if_hwratio*ta*tc)/pi+(1.-if_fbeam)*roadpsi
end where

! Calculate short wave reflections to nrefl order
roadnetalpha=rdsndelta*ird_alpha+(1.-rdsndelta)*((1.-if_vegsigmac)*if_roadalpha+if_vegsigmac*if_vegalphac)
sg_walle=walles
sg_wallw=wallws
sg_road=roads
do k=1,nrefl
  nwalles=roadnetalpha*wallpsi*roads+if_wallalpha*(1.-2.*wallpsi)*wallws
  nwallws=roadnetalpha*wallpsi*roads+if_wallalpha*(1.-2.*wallpsi)*walles
  nroads=if_wallalpha*(1.-roadpsi)*0.5*(walles+wallws)
  walles=nwalles
  wallws=nwallws
  roads=nroads
  sg_walle=sg_walle+walles
  sg_wallw=sg_wallw+wallws
  sg_road=sg_road+roads
end do
sg_roof=1.
sg_vegr=1.
sg_rfsn=1.
sg_rdsn=sg_road
sg_vegc=sg_road

return
end subroutine getswcoeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate longwave radiation coefficents (modified to include 2nd wall)

subroutine getlwcoeff(d_netemiss,d_cwa,d_cra,d_cw0,d_cww,d_crw,d_crr,d_cwr,d_rdsndelta,wallpsi,roadpsi,if_vegsigmac, &
                      if_roademiss,if_vegemissc,if_wallemiss)

implicit none

integer k
real, dimension(:), intent(inout) :: d_netemiss
real, dimension(size(d_netemiss)), intent(inout) :: d_cwa,d_cra,d_cw0,d_cww,d_crw,d_crr,d_cwr,d_rdsndelta
real, dimension(size(d_netemiss)), intent(in) :: if_vegsigmac,if_roademiss,if_vegemissc,if_wallemiss
real, dimension(size(d_netemiss)), intent(in) :: wallpsi,roadpsi
real, dimension(size(d_netemiss)) :: rcwa,rcra,rcwe,rcww,rcrw,rcrr,rcwr
real, dimension(size(d_netemiss)) :: ncwa,ncra,ncwe,ncww,ncrw,ncrr,ncwr


d_netemiss=d_rdsndelta*snowemiss+(1.-d_rdsndelta)*((1.-if_vegsigmac)*if_roademiss+if_vegsigmac*if_vegemissc)
d_cwa=wallpsi
d_cra=roadpsi
d_cw0=0.
d_cww=1.-2.*wallpsi
d_crw=0.5*(1.-roadpsi)
d_crr=0.
d_cwr=wallpsi
rcwa=d_cwa
rcra=d_cra
rcwe=d_cw0
rcww=d_cww
rcrw=d_crw
rcrr=d_crr
rcwr=d_cwr
do k=1,nrefl
  ncwa=(1.-d_netemiss)*wallpsi*rcra+(1.-if_wallemiss)*(1.-2.*wallpsi)*rcwa
  ncra=(1.-if_wallemiss)*(1.-roadpsi)*rcwa
  ncwe=(1.-d_netemiss)*wallpsi*rcrw+(1.-if_wallemiss)*(1.-2.*wallpsi)*rcww
  ncww=(1.-d_netemiss)*wallpsi*rcrw+(1.-if_wallemiss)*(1.-2.*wallpsi)*rcwe
  ncrw=(1.-if_wallemiss)*(1.-roadpsi)*0.5*(rcww+rcwe)  
  ncwr=(1.-d_netemiss)*wallpsi*rcrr+(1.-if_wallemiss)*(1.-2.*wallpsi)*rcwr
  ncrr=(1.-if_wallemiss)*(1.-roadpsi)*rcwr
  rcwa=ncwa
  rcra=ncra
  rcwe=ncwe
  rcww=ncww
  rcrw=ncrw
  rcrr=ncrr
  rcwr=ncwr
  d_cwa=d_cwa+rcwa
  d_cra=d_cra+rcra
  d_cw0=d_cw0+rcwe
  d_cww=d_cww+rcww
  d_crw=d_crw+rcrw
  d_cwr=d_cwr+rcwr
  d_crr=d_crr+rcrr
end do

end subroutine getlwcoeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solve for road snow temperature (includes vegetation canopy temperature and canyon temperature)

subroutine solvecanyon(sg_road,rg_road,fg_road,eg_road,acond_road,abase_road,                          &
                       sg_walle,rg_walle,fg_walle,acond_walle,abase_walle,                             &
                       sg_wallw,rg_wallw,fg_wallw,acond_wallw,abase_wallw,                             &
                       sg_vegc,rg_vegc,fg_vegc,eg_vegc,acond_vegc,abase_vegc,                          &
                       sg_rdsn,rg_rdsn,fg_rdsn,eg_rdsn,acond_rdsn,abase_rdsn,rdsntemp,rdsnmelt,gardsn, &
                       a_umag,a_rho,a_rg,a_rnd,a_snd,                                                  &
                       d_canyontemp,d_canyonmix,d_tempc,d_mixrc,d_sigd,d_topu,d_netrad,                &
                       d_roaddelta,d_vegdeltac,d_rdsndelta,d_ac_outside,d_traf,d_ac_inside,            &
                       d_canyonrgout,d_tranc,d_evapc,d_cwa,d_cra,d_cw0,d_cww,d_crw,d_crr,              &
                       d_cwr,d_totdepth,d_c1c,d_intgains_bld,fgtop,egtop,int_infilflux,int_newairtemp, &
                       int_infilfg,ggint_roof,ggint_walle,ggint_wallw,ggint_road,ggint_slab,           &
                       ggint_intm1,ggint_intm2,cyc_translation,cyc_proportion,ddt,                     &
                       cnveg,if_ach,if_bldairtemp,if_bldheight,if_bldwidth,if_coeffbldheight,          &
                       if_effhwratio,if_hwratio,if_intm,if_intmassn,if_road,if_roof,if_sigmabld,       &
                       if_slab,if_tempcool,if_tempheat,if_wall,intm,p_bldcool,p_bldheat,p_cndzmin,     &
                       p_lzoh,p_lzom,rdhyd,rfveg,road,roof,room,slab,walle,wallw,if_sfc,if_swilt,ufull)
implicit none

integer, intent(in) :: ufull
integer k,l
real, intent(in)    :: ddt
real, dimension(ufull), intent(inout) :: rg_road,fg_road,eg_road,abase_road
real, dimension(ufull), intent(inout) :: rg_walle,fg_walle,abase_walle
real, dimension(ufull), intent(inout) :: rg_wallw,fg_wallw,abase_wallw
real, dimension(ufull), intent(inout) :: rg_vegc,fg_vegc,eg_vegc,abase_vegc
real, dimension(ufull), intent(inout) :: rg_rdsn,fg_rdsn,eg_rdsn,abase_rdsn,rdsntemp,rdsnmelt,gardsn
real, dimension(ufull), intent(in) :: sg_road,sg_walle,sg_wallw,sg_vegc,sg_rdsn
real, dimension(ufull), intent(in) :: a_umag,a_rho,a_rg,a_rnd,a_snd
real, dimension(ufull), intent(out) :: int_newairtemp
real, dimension(ufull), intent(out) :: d_ac_inside
real, dimension(ufull), intent(inout) :: d_canyontemp,d_canyonmix,d_tempc,d_mixrc,d_sigd,d_topu,d_netrad
real, dimension(ufull), intent(inout) :: d_roaddelta,d_vegdeltac,d_rdsndelta,d_ac_outside,d_traf
real, dimension(ufull), intent(inout) :: d_canyonrgout,d_tranc,d_evapc,d_cwa,d_cra,d_cw0,d_cww,d_crw,d_crr,d_cwr
real, dimension(ufull), intent(inout) :: d_totdepth,d_c1c,d_intgains_bld
real, dimension(ufull), intent(out) :: fgtop,egtop,int_infilflux
real, dimension(ufull), intent(out) :: acond_road,acond_walle,acond_wallw,acond_vegc,acond_rdsn
real, dimension(ufull), intent(out) :: int_infilfg
real, dimension(ufull), intent(out) :: ggint_roof, ggint_walle, ggint_wallw, ggint_road
real, dimension(ufull), intent(out) :: ggint_slab, ggint_intm1, ggint_intm2
real, dimension(ufull) :: newval,sndepth,snlambda,ldratio,roadqsat,vegqsat,rdsnqsat
real, dimension(ufull) :: cu,topinvres,dts,dtt,cduv,z_on_l,dumroaddelta,dumvegdelta,res
real, dimension(ufull) :: effwalle,effwallw,effroad,effrdsn,effvegc
real, dimension(ufull) :: aa,bb,cc,dd,ee,ff
real, dimension(ufull) :: lwflux_walle_road, lwflux_wallw_road, lwflux_walle_rdsn, lwflux_wallw_rdsn
real, dimension(ufull) :: lwflux_walle_vegc, lwflux_wallw_vegc
real, dimension(ufull) :: skintemp, ac_coeff
real, dimension(ufull) :: ac_load,cyc_translation,cyc_proportion
real, dimension(ufull) :: cvcoeff_roof,cvcoeff_walle,cvcoeff_wallw,cvcoeff_slab,cvcoeff_intm1,cvcoeff_intm2
real, dimension(ufull,2) :: evct,evctx,oldval
real, dimension(ufull), intent(in) :: if_ach, if_bldairtemp, if_bldheight, if_bldwidth, if_coeffbldheight
real, dimension(ufull), intent(in) :: if_effhwratio, if_hwratio, if_sigmabld, if_tempcool, if_tempheat
real, dimension(ufull), intent(in) :: if_sfc, if_swilt
integer, dimension(ufull), intent(in) :: if_intmassn
type(facetparams), intent(in) :: if_intm, if_road, if_roof, if_slab, if_wall
type(facetdata), intent(in) :: intm
real, dimension(ufull), intent(inout) :: p_bldcool, p_bldheat, p_cndzmin, p_lzoh
real, dimension(ufull), intent(in) :: p_lzom
type(hydrodata), intent(in) :: rdhyd
type(vegdata), intent(inout) :: cnveg, rfveg
type(facetdata), intent(in) :: roof, slab
type(facetdata), intent(inout) :: road, room, walle, wallw

! snow conductance
sndepth  = rdhyd%snow*waterden/rdhyd%den
snlambda = icelambda*(rdhyd%den/waterden)**1.88

! first guess for canyon air temperature and water vapor mixing ratio
! also guess for canyon veg and snow temperatures
d_canyontemp    = d_tempc
d_canyonmix     = d_mixrc
cnveg%temp      = d_tempc
rdsntemp        = road%nodetemp(:,1)
rdsnmelt        = 0.
dumvegdelta     = 0. ! cray compiler bug
if ( conductmeth==0 ) then
  road%nodetemp(:,0)  = road%nodetemp(:,1)
  walle%nodetemp(:,0) = walle%nodetemp(:,1)
  wallw%nodetemp(:,0) = wallw%nodetemp(:,1)
end if
d_netrad=sbconst*(d_rdsndelta*snowemiss*(rdsntemp+urbtemp)**4                            &
        +(1.-d_rdsndelta)*(1.-cnveg%sigma)*if_road%emiss*(road%nodetemp(:,0)+urbtemp)**4 &
        +(1.-d_rdsndelta)*cnveg%sigma*cnveg%emiss*(cnveg%temp+urbtemp)**4)

! Solve for canyon air temperature and water vapor mixing ratio
do l = 1,ncyits

  ! interior model
  ggint_road = 0.
  ggint_slab = 0.
  ggint_intm1 = 0.
  ggint_intm2 = 0.
  ! first internal temperature estimation - used for ggint calculation
  select case(intairtmeth)
    case(0) ! fixed internal air temperature
      room%nodetemp(:,1) = if_bldairtemp
      call calc_convcoeff(cvcoeff_roof,cvcoeff_walle,cvcoeff_wallw,cvcoeff_slab,  & 
                          cvcoeff_intm1,cvcoeff_intm2,roof,room,slab,ufull)
      ! (use split form to estimate G_{*,4} flux into room for AC.  newtemp is an estimate of the temperature at tau+1)
      call calc_ggint(if_roof%depth(:,nl),if_roof%volcp(:,nl),if_roof%lambda(:,nl),roof%nodetemp(:,nl),  &
                      if_bldairtemp,cvcoeff_roof, ddt, ggint_roof,ufull)
      call calc_ggint(if_wall%depth(:,nl),if_wall%volcp(:,nl),if_wall%lambda(:,nl),walle%nodetemp(:,nl), &
                      if_bldairtemp,cvcoeff_walle, ddt, ggint_walle,ufull)
      call calc_ggint(if_wall%depth(:,nl),if_wall%volcp(:,nl),if_wall%lambda(:,nl),wallw%nodetemp(:,nl), &
                      if_bldairtemp,cvcoeff_wallw, ddt, ggint_wallw,ufull)

      ! flux into room potentially pumped out into canyon (depends on AC method)
      d_ac_inside = -(1.-rfveg%sigma)*ggint_roof - ggint_slab                 & 
                  - (ggint_intm1+ggint_intm2)*if_intmassn                     &
                  - (ggint_walle+ggint_wallw)*(if_bldheight/if_bldwidth)      &
                  - d_intgains_bld
    
    case(1) ! floating internal air temperature
      ! estimate internal surface convection coefficients
      call calc_convcoeff(cvcoeff_roof,cvcoeff_walle,cvcoeff_wallw,cvcoeff_slab,       & 
                          cvcoeff_intm1,cvcoeff_intm2,roof,room,slab,ufull)
      ! estimate new internal air temperature
      call calc_newairtemp(int_newairtemp,a_rho,d_canyontemp,d_intgains_bld,           &
                           cvcoeff_roof,cvcoeff_walle,cvcoeff_wallw,cvcoeff_slab,      &
                           cvcoeff_intm1,cvcoeff_intm2,ddt,                            &
                           if_ach,if_bldheight,if_bldwidth,if_intmassn,intm,roof,      &
                           room,slab,walle,wallw,ufull)
  
      d_ac_inside=0.
      where (int_newairtemp>if_tempcool+cyc_translation-urbtemp)
        ac_load = -(a_rho*aircp*if_bldheight/ddt)*(int_newairtemp+(urbtemp-if_tempcool))
        d_ac_inside = max(-ac_cap*if_bldheight,ac_load)*cyc_proportion
      end where
      where (int_newairtemp<if_tempheat+cyc_translation-urbtemp)
        ac_load = -(a_rho*aircp*if_bldheight/ddt)*(int_newairtemp+(urbtemp-if_tempheat))
        d_ac_inside = min(ac_cap*if_bldheight,ac_load)*cyc_proportion
      end where
    
      call calc_ggint(if_roof%depth(:,nl),if_roof%volcp(:,nl),if_roof%lambda(:,nl),roof%nodetemp(:,nl),   &
                      int_newairtemp,cvcoeff_roof, ddt, ggint_roof,ufull)
      call calc_ggint(if_wall%depth(:,nl),if_wall%volcp(:,nl),if_wall%lambda(:,nl),walle%nodetemp(:,nl),  &
                      int_newairtemp,cvcoeff_walle, ddt, ggint_walle,ufull)
      call calc_ggint(if_wall%depth(:,nl),if_wall%volcp(:,nl),if_wall%lambda(:,nl),wallw%nodetemp(:,nl),  &
                      int_newairtemp,cvcoeff_wallw, ddt, ggint_wallw,ufull)

      call calc_ggint(if_slab%depth(:,nl),if_slab%volcp(:,nl),if_slab%lambda(:,nl),slab%nodetemp(:,nl),   &
                      int_newairtemp,cvcoeff_slab, ddt, ggint_slab,ufull)
      if (intmassmeth/=0) then
        call calc_ggint(if_intm%depth(:,1),if_intm%volcp(:,1),if_intm%lambda(:,1),intm%nodetemp(:,0),     &
                        int_newairtemp,cvcoeff_intm1, ddt, ggint_intm1,ufull)  
        call calc_ggint(if_intm%depth(:,nl),if_intm%volcp(:,nl),if_intm%lambda(:,nl),intm%nodetemp(:,nl), &
                        int_newairtemp,cvcoeff_intm2, ddt, ggint_intm2,ufull)
      end if
                
    case DEFAULT
      write(6,*) "ERROR: Unknown intairtmeth mode ",intairtmeth
      stop
  end select

  !  solve for aerodynamical resistance between canyon and atmosphere  
  ! assume zoh=zom when coupling to canyon air temperature
  p_lzoh = p_lzom
  dts    = d_canyontemp + (d_canyontemp+urbtemp)*0.61*d_canyonmix
  dtt    = d_tempc + (d_tempc+urbtemp)*0.61*d_mixrc
  call getinvres(topinvres,cduv,z_on_l,p_lzoh,p_lzom,p_cndzmin,dts,dtt,a_umag,3)
  call gettopu(d_topu,a_umag,z_on_l,if_bldheight,cduv,p_cndzmin,if_hwratio,ufull)

  if ( resmeth==0 ) then
    acond_road  = (11.8+4.2*sqrt((d_topu*abase_road)**2+cduv*a_umag**2))/(aircp*a_rho)  ! From Rowley, et al (1930)
    acond_walle = acond_road
    acond_wallw = acond_road
    acond_rdsn  = acond_road
    acond_vegc  = acond_road
  else if ( resmeth==2 ) then
    cu = abase_road*d_topu
    where (cu<=5.)
      acond_road = (6.15+4.18*cu)/(aircp*a_rho)
    elsewhere
      acond_road = (7.51*cu**0.78)/(aircp*a_rho)
    end where
    acond_walle = acond_road
    acond_wallw = acond_road
    acond_rdsn  = acond_road
    acond_vegc  = acond_road
  else
    acond_road  = d_topu*abase_road  
    acond_walle = d_topu*abase_walle
    acond_wallw = d_topu*abase_wallw
    acond_vegc  = d_topu*abase_vegc
    acond_rdsn  = d_topu*abase_rdsn
  end if

  ! saturated mixing ratio for road
  call getqsat(roadqsat,road%nodetemp(:,0),d_sigd)   ! evaluate using pressure at displacement height
  
  ! correction for dew
  where (roadqsat<d_canyonmix)
    dumroaddelta=1.
  elsewhere
    dumroaddelta=d_roaddelta
  end where
  
  ! calculate canyon road latent heat flux
  aa=rdhyd%surfwater/ddt+a_rnd+rdsnmelt
  eg_road=lv*min(a_rho*d_roaddelta*(roadqsat-d_canyonmix)*acond_road,aa)
  
  if ( conductmeth==0 ) then ! half-layer diagnostic skin temperature estimate
    ! calculate road and wall skin temperatures
    ! Write energy budget
    !     Solar_net + Longwave_net - Sensible flux - Latent flux - Conduction = 0
    ! or 
    !     sg + a_rg - f_emiss*sbconst*Tskin**4 - aircp*a_rho*(Tskin-d_tempc) &
    !     -eg - (Tskin-temp(:,1))/ldrratio = 0
    ! as a quartic equation
    !      aa*Tskin^4 + dd*Tskin + ee = 0
    ! and solve for Tskin  
    effwalle=if_wall%emiss*(a_rg*d_cwa+sbconst*(walle%nodetemp(:,0)+urbtemp)**4*if_wall%emiss*d_cw0                  & 
                    +sbconst*(wallw%nodetemp(:,0)+urbtemp)**4*if_wall%emiss*d_cww+d_netrad*d_cwr)
    effwallw=if_wall%emiss*(a_rg*d_cwa+sbconst*(wallw%nodetemp(:,0)+urbtemp)**4*if_wall%emiss*d_cw0                  &
                    +sbconst*(walle%nodetemp(:,0)+urbtemp)**4*if_wall%emiss*d_cww+d_netrad*d_cwr)
    effroad=if_road%emiss*(a_rg*d_cra+(d_netrad*d_crr-sbconst*(road%nodetemp(:,0)+urbtemp)**4)                       &
                    +sbconst*if_wall%emiss*((walle%nodetemp(:,0)+urbtemp)**4+(wallw%nodetemp(:,0)+urbtemp)**4)*d_crw)
    ldratio = 0.5*if_wall%depth(:,1)/if_wall%lambda(:,1)
    aa = if_wall%emiss*sbconst
    dd = aircp*a_rho*acond_walle+1./ldratio
    ee = -sg_walle-effwalle-aircp*a_rho*acond_walle*(d_canyontemp+urbtemp)-(walle%nodetemp(:,1)+urbtemp)/ldratio
    call solvequartic(skintemp,aa,dd,ee) ! This is an estimate of Tskin to be updated in solvetridiag
    walle%nodetemp(:,0) = skintemp - urbtemp
    dd = aircp*a_rho*acond_wallw+1./ldratio
    ee = -sg_wallw-effwallw-aircp*a_rho*acond_wallw*(d_canyontemp+urbtemp)-(wallw%nodetemp(:,1)+urbtemp)/ldratio
    call solvequartic(skintemp,aa,dd,ee) ! This is an estimate of Tskin to be updated in solvetridiag
    wallw%nodetemp(:,0) = skintemp - urbtemp
    ldratio = 0.5*if_road%depth(:,1)/if_road%lambda(:,1)
    aa = if_road%emiss*sbconst
    dd = aircp*a_rho*acond_road+1./ldratio
    ee = -sg_road-effroad-aircp*a_rho*acond_road*(d_canyontemp+urbtemp)-(road%nodetemp(:,1)+urbtemp)/ldratio+eg_road
    call solvequartic(skintemp,aa,dd,ee)  ! This is an estimate of Tskin to be updated in solvetridiag
    road%nodetemp(:,0) = skintemp - urbtemp
  end if
  ! Calculate longwave radiation emitted from the canyon floor
  ! MJT notes - This could be included within the iterative solver for snow and vegetation temperatures.
  ! However, it creates a (weak) coupling between these two variables and therefore could require
  ! a multivariate root finding method (e.g,. Broyden's method). Instead we explicitly solve for d_netrad, 
  ! which allows us to decouple the solutions for snow and vegtation temperatures.
  d_netrad=sbconst*(d_rdsndelta*snowemiss*(rdsntemp+urbtemp)**4                             &
          +(1.-d_rdsndelta)*((1.-cnveg%sigma)*if_road%emiss*(road%nodetemp(:,0)+urbtemp)**4 &
          +cnveg%sigma*cnveg%emiss*(cnveg%temp+urbtemp)**4))
  
  if ( lweff/=1 ) then
    lwflux_walle_road = 0.
    lwflux_wallw_road = 0.
    lwflux_walle_rdsn = 0.
    lwflux_wallw_rdsn = 0.
    lwflux_walle_vegc = 0.
    lwflux_wallw_vegc = 0.
  else
    lwflux_walle_road = sbconst*(if_road%emiss*(road%nodetemp(:,0)+urbtemp)**4         &
                       -if_wall%emiss*(walle%nodetemp(:,0)+urbtemp)**4)*(1.-if_coeffbldheight)
    lwflux_wallw_road = sbconst*(if_road%emiss*(road%nodetemp(:,0)+urbtemp)**4         &
                       -if_wall%emiss*(wallw%nodetemp(:,0)+urbtemp)**4)*(1.-if_coeffbldheight)
    lwflux_walle_rdsn = sbconst*(snowemiss*(rdsntemp+urbtemp)**4                       &
                       -if_wall%emiss*(walle%nodetemp(:,0)+urbtemp)**4)*(1.-if_coeffbldheight)
    lwflux_wallw_rdsn = sbconst*(snowemiss*(rdsntemp+urbtemp)**4                       &
                       -if_wall%emiss*(wallw%nodetemp(:,0)+urbtemp)**4)*(1.-if_coeffbldheight)
    lwflux_walle_vegc = sbconst*(cnveg%emiss*(cnveg%temp+urbtemp)**4                   &
                       -if_wall%emiss*(walle%nodetemp(:,0)+urbtemp)**4)*(1.-if_coeffbldheight)
    lwflux_wallw_vegc = sbconst*(cnveg%emiss*(cnveg%temp+urbtemp)**4                   &
                       -if_wall%emiss*(wallw%nodetemp(:,0)+urbtemp)**4)*(1.-if_coeffbldheight)
  end if
  
  ! solve for road snow and canyon veg temperatures -------------------------------
  ldratio  = 0.5*( sndepth/snlambda + if_road%depth(:,1)/if_road%lambda(:,1) )
  oldval(:,1) = cnveg%temp + 0.5
  oldval(:,2) = rdsntemp + 0.5
  call canyonflux(evct,sg_vegc,rg_vegc,fg_vegc,eg_vegc,acond_vegc,vegqsat,res,dumvegdelta,      &
                  sg_rdsn,rg_rdsn,fg_rdsn,eg_rdsn,acond_rdsn,rdsntemp,gardsn,rdsnmelt,rdsnqsat, &
                  a_rg,a_rho,a_rnd,a_snd,                                                       &
                  d_canyontemp,d_canyonmix,d_sigd,d_netrad,d_tranc,d_evapc,                     &
                  d_cra,d_crr,d_crw,d_totdepth,d_c1c,d_vegdeltac,                               &
                  effvegc,effrdsn,ldratio,lwflux_walle_rdsn,lwflux_wallw_rdsn,                  &
                  lwflux_walle_vegc,lwflux_wallw_vegc,ddt,                                      &
                  cnveg,if_sfc,if_swilt,if_wall,rdhyd,road,walle,wallw,ufull)
  cnveg%temp = cnveg%temp - 0.5
  rdsntemp   = rdsntemp - 0.5
  do k = 1,nfgits ! sectant
    evctx = evct
    call canyonflux(evct,sg_vegc,rg_vegc,fg_vegc,eg_vegc,acond_vegc,vegqsat,res,dumvegdelta,      &
                    sg_rdsn,rg_rdsn,fg_rdsn,eg_rdsn,acond_rdsn,rdsntemp,gardsn,rdsnmelt,rdsnqsat, &
                    a_rg,a_rho,a_rnd,a_snd,                                                       &
                    d_canyontemp,d_canyonmix,d_sigd,d_netrad,d_tranc,d_evapc,                     &
                    d_cra,d_crr,d_crw,d_totdepth,d_c1c,d_vegdeltac,                               &
                    effvegc,effrdsn,ldratio,lwflux_walle_rdsn,lwflux_wallw_rdsn,                  &
                    lwflux_walle_vegc,lwflux_wallw_vegc,ddt,                                      &
                    cnveg,if_sfc,if_swilt,if_wall,rdhyd,road,walle,wallw,ufull)
    evctx = evct-evctx
    where (abs(evctx(:,1))>tol)
      newval      = max(min(cnveg%temp-alpha*evct(:,1)*(cnveg%temp-oldval(:,1))/evctx(:,1),400.-urbtemp),200.-urbtemp)
      oldval(:,1) = cnveg%temp
      cnveg%temp  = newval
    end where
    where (abs(evctx(:,2))>tol)
      newval      = max(min(rdsntemp-alpha*evct(:,2)*(rdsntemp-oldval(:,2))/evctx(:,2), 300.-urbtemp),100.-urbtemp)
      oldval(:,2) = rdsntemp
      rdsntemp    = newval
    end where
  end do
  ! ---------------------------------------------------------------    

  ! balance canyon latent heat budget
  aa = d_rdsndelta*acond_rdsn
  bb = (1.-d_rdsndelta)*(1.-cnveg%sigma)*dumroaddelta*acond_road
  cc = (1.-d_rdsndelta)*cnveg%sigma*(dumvegdelta*acond_vegc+(1.-dumvegdelta)/(1./max(acond_vegc,1.e-10)+res))
  dd = topinvres
  d_canyonmix = (aa*rdsnqsat+bb*roadqsat+cc*vegqsat+dd*d_mixrc)/(aa+bb+cc+dd)
  
  ac_coeff = acfactor*max(d_canyontemp-room%nodetemp(:,1),0.)/(room%nodetemp(:,1)+urbtemp)    ! T&H Eq. 10
  ! update heat pumped into canyon
  select case(acmeth) ! AC heat pump into canyon (0=Off, 1=On, 2=Reversible, COP of 1.0)
    case(0) ! unrealistic cooling (buildings act as heat sink)
      d_ac_outside  = 0.
      p_bldheat = max(0.,d_ac_inside*if_sigmabld)
      p_bldcool = max(0.,d_ac_inside*if_sigmabld)
    case(1) ! d_ac_outside pumps conducted heat + ac waste heat back into canyon
      d_ac_outside = max(0.,-d_ac_inside*(1.+ac_coeff)*if_sigmabld/(1.-if_sigmabld))  ! canyon domain W/m/m
      p_bldheat = max(0.,d_ac_inside*if_sigmabld)                                    ! entire domain W/m/m
      p_bldcool = max(0.,-d_ac_inside*ac_coeff*if_sigmabld)                          ! entire domain W/m/m
    case(2) ! reversible heating and cooling (for testing energy conservation)
      d_ac_outside  = -d_ac_inside*if_sigmabld/(1.-if_sigmabld)
      p_bldheat = 0.
      p_bldcool = 0.
    case DEFAULT
      write(6,*) "ERROR: Unknown acmeth mode ",acmeth
      stop
  end select
  ! update infiltration between canyon and room
  select case(intairtmeth)
    case(0)  ! fixed internal temperature
      int_infilflux = 0.
      int_infilfg = 0.
    case(1)
      int_infilflux = (if_ach/3600.)*aircp*a_rho*if_bldheight*(d_canyontemp-int_newairtemp)
      int_infilfg = int_infilflux*if_sigmabld/(1.-if_sigmabld)
    case DEFAULT
      write(6,*) "ERROR: Unknown intairtmeth ",intairtmeth
      stop
  end select

  ! balance sensible heat flux
  aa = aircp*a_rho*topinvres
  bb = d_rdsndelta*aircp*a_rho*acond_rdsn
  cc = (1.-d_rdsndelta)*(1.-cnveg%sigma)*aircp*a_rho*acond_road
  dd = (1.-d_rdsndelta)*cnveg%sigma*aircp*a_rho*acond_vegc
  ee = if_effhwratio*aircp*a_rho*acond_walle
  ff = if_effhwratio*aircp*a_rho*acond_wallw
  !!!!!! infiltration on !!!!!!!
  d_canyontemp = (aa*d_tempc+bb*rdsntemp+cc*road%nodetemp(:,0)+dd*cnveg%temp+ee*walle%nodetemp(:,0) & 
                +ff*wallw%nodetemp(:,0)+d_traf+d_ac_outside-int_infilfg)/(aa+bb+cc+dd+ee+ff)
end do
! solve for canyon sensible heat flux
fg_walle = aircp*a_rho*(walle%nodetemp(:,0)-d_canyontemp)*acond_walle*if_coeffbldheight ! canyon vegetation blocks turblent flux
fg_wallw = aircp*a_rho*(wallw%nodetemp(:,0)-d_canyontemp)*acond_wallw*if_coeffbldheight ! canyon vegetation blocks turblent flux
fg_road  = aircp*a_rho*(road%nodetemp(:,0)-d_canyontemp)*acond_road
fg_vegc  = sg_vegc+rg_vegc-eg_vegc
fg_rdsn  = sg_rdsn+rg_rdsn-eg_rdsn-lf*rdsnmelt-gardsn*(1.-cnveg%sigma)
fgtop = if_hwratio*(fg_walle+fg_wallw) + (1.-d_rdsndelta)*(1.-cnveg%sigma)*fg_road &
      + (1.-d_rdsndelta)*cnveg%sigma*fg_vegc + d_rdsndelta*fg_rdsn                 &
      + d_traf + d_ac_outside - int_infilfg

! solve for canyon latent heat flux
egtop = (1.-d_rdsndelta)*(1.-cnveg%sigma)*eg_road + (1.-d_rdsndelta)*cnveg%sigma*eg_vegc &
      + d_rdsndelta*eg_rdsn

! calculate longwave radiation
if ( lweff/=2 ) then
  effwalle=if_wall%emiss*(a_rg*d_cwa+sbconst*(walle%nodetemp(:,0)+urbtemp)**4*(if_wall%emiss*d_cw0-1.)                & 
                                  +sbconst*(wallw%nodetemp(:,0)+urbtemp)**4*if_wall%emiss*d_cww+d_netrad*d_cwr)
  rg_walle=effwalle*if_coeffbldheight+lwflux_walle_road*(1.-d_rdsndelta)*(1.-cnveg%sigma)/if_hwratio                  &
                                  +lwflux_walle_vegc*(1.-d_rdsndelta)*cnveg%sigma/if_hwratio                          &
                                  +lwflux_walle_rdsn*d_rdsndelta/if_hwratio
  effwallw=if_wall%emiss*(a_rg*d_cwa+sbconst*(wallw%nodetemp(:,0)+urbtemp)**4*(if_wall%emiss*d_cw0-1.)                &
                                  +sbconst*(walle%nodetemp(:,0)+urbtemp)**4*if_wall%emiss*d_cww+d_netrad*d_cwr)
  rg_wallw=effwallw*if_coeffbldheight+lwflux_wallw_road*(1.-d_rdsndelta)*(1.-cnveg%sigma)/if_hwratio                  &
                                  +lwflux_wallw_vegc*(1.-d_rdsndelta)*cnveg%sigma/if_hwratio                          &
                                  +lwflux_wallw_rdsn*d_rdsndelta/if_hwratio
  effroad=if_road%emiss*(a_rg*d_cra+(d_netrad*d_crr-sbconst*(road%nodetemp(:,0)+urbtemp)**4)                          &
                    +sbconst*if_wall%emiss*((walle%nodetemp(:,0)+urbtemp)**4+(wallw%nodetemp(:,0)+urbtemp)**4)*d_crw)
  rg_road=effroad-lwflux_walle_road-lwflux_wallw_road
else
  effwalle=if_wall%emiss*(a_rg*d_cwa+sbconst*(walle%nodetemp(:,0)+urbtemp)**4*(if_wall%emiss*d_cw0-1.)                & 
                                  +sbconst*(wallw%nodetemp(:,0)+urbtemp)**4*if_wall%emiss*d_cww+d_netrad*d_cwr)
  rg_walle=effwalle
  effwallw=if_wall%emiss*(a_rg*d_cwa+sbconst*(wallw%nodetemp(:,0)+urbtemp)**4*(if_wall%emiss*d_cw0-1.)                &
                                  +sbconst*(walle%nodetemp(:,0)+urbtemp)**4*if_wall%emiss*d_cww+d_netrad*d_cwr)
  rg_wallw=effwallw
  effroad=if_road%emiss*(a_rg*d_cra+(d_netrad*d_crr-sbconst*(road%nodetemp(:,0)+urbtemp)**4)                          &
                    +sbconst*if_wall%emiss*((walle%nodetemp(:,0)+urbtemp)**4+(wallw%nodetemp(:,0)+urbtemp)**4)*d_crw)
  rg_road=effroad
end if

! outgoing longwave radiation
! note that eff terms are used for outgoing longwave radiation, whereas rg terms are used for heat conduction
if ( lweff/=2 ) then
  d_canyonrgout=a_rg-d_rdsndelta*effrdsn-(1.-d_rdsndelta)*((1.-cnveg%sigma)*effroad+cnveg%sigma*effvegc)            &
                    -if_hwratio*if_coeffbldheight*(effwalle+effwallw)
else
  d_canyonrgout=a_rg-d_rdsndelta*effrdsn-(1.-d_rdsndelta)*((1.-cnveg%sigma)*effroad+cnveg%sigma*effvegc)            &
                    -if_hwratio*(effwalle+effwallw)
end if
!0. = d_rdsndelta*(lwflux_walle_rdsn+lwflux_wallw_rdsn)                            &
!     +(1.-d_rdsndelta)*((1.-cnveg%sigma)*(lwflux_walle_road+lwflux_wallw_road)    &
!    +cnveg%sigma*(lwflux_walle_vegc+lwflux_wallw_vegc))                           &
!    -if_hwratio*(lwflux_walle_road*(1.-d_rdsndelta)*(1.-cnveg%sigma)/if_hwratio   &
!    +lwflux_walle_vegc*(1.-d_rdsndelta)*cnveg%sigma/if_hwratio                    &
!     +lwflux_walle_rdsn*d_rdsndelta/if_hwratio)                                   &
!    - if_hwratio*(lwflux_wallw_road*(1.-d_rdsndelta)*(1.-cnveg%sigma)/if_hwratio  &
!    +lwflux_wallw_vegc*(1.-d_rdsndelta)*cnveg%sigma/if_hwratio                    &
!    +lwflux_wallw_rdsn*d_rdsndelta/if_hwratio)

!write(6,*) 'd_canyontemp, room%nodetemp',d_canyontemp, room%nodetemp

return
end subroutine solvecanyon

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solve for canyon veg and snow fluxes
                     
subroutine canyonflux(evct,sg_vegc,rg_vegc,fg_vegc,eg_vegc,acond_vegc,vegqsat,res,dumvegdelta,       &
                      sg_rdsn,rg_rdsn,fg_rdsn,eg_rdsn,acond_rdsn,rdsntemp,gardsn,rdsnmelt,rdsnqsat,  &
                      a_rg,a_rho,a_rnd,a_snd,                                                        &
                      d_canyontemp,d_canyonmix,d_sigd,d_netrad,d_tranc,d_evapc,                      &
                      d_cra,d_crr,d_crw,d_totdepth,d_c1c,d_vegdeltac,                                &
                      effvegc,effrdsn,ldratio,lwflux_walle_rdsn,lwflux_wallw_rdsn,                   &
                      lwflux_walle_vegc,lwflux_wallw_vegc,ddt, &
                      cnveg,if_sfc,if_swilt,if_wall,rdhyd,road,walle,wallw,ufull)

implicit none

integer, intent(in) :: ufull
real, intent(in) :: ddt
real, dimension(ufull,2), intent(out) :: evct
real, dimension(ufull), intent(inout) :: rg_vegc,fg_vegc,eg_vegc,acond_vegc,vegqsat,res,dumvegdelta
real, dimension(ufull), intent(inout) :: rg_rdsn,fg_rdsn,eg_rdsn,acond_rdsn,rdsntemp,gardsn,rdsnmelt,rdsnqsat
real, dimension(ufull), intent(in) :: sg_vegc,sg_rdsn
real, dimension(ufull), intent(in) :: a_rg,a_rho,a_rnd,a_snd
real, dimension(ufull), intent(in) :: ldratio,lwflux_walle_rdsn,lwflux_wallw_rdsn,lwflux_walle_vegc,lwflux_wallw_vegc
real, dimension(ufull), intent(out) :: effvegc,effrdsn
real, dimension(ufull), intent(inout) :: d_canyontemp,d_canyonmix,d_sigd,d_netrad,d_tranc,d_evapc
real, dimension(ufull), intent(inout) :: d_cra,d_crr,d_crw,d_totdepth,d_c1c,d_vegdeltac
real, dimension(ufull) :: ff,f1,f2,f3,f4
real, dimension(ufull) :: snevap
type(vegdata), intent(in) :: cnveg
real, dimension(ufull), intent(in) :: if_sfc, if_swilt
type(facetparams), intent(in) :: if_wall
type(hydrodata), intent(in) :: rdhyd
type(facetdata), intent(in) :: road, walle, wallw

! estimate mixing ratio for vegetation and snow
call getqsat(vegqsat,cnveg%temp,d_sigd)
call getqsat(rdsnqsat,rdsntemp,d_sigd)

! correction for dew
where (vegqsat<d_canyonmix)
  dumvegdelta=1.
elsewhere
  dumvegdelta=d_vegdeltac
end where
  
! vegetation transpiration terms (developed by Eva in CCAM sflux.f and CSIRO9)
where (cnveg%zo<0.5)
  ff=1.1*sg_vegc/max(cnveg%lai*150.,1.E-8)
elsewhere
  ff=1.1*sg_vegc/max(cnveg%lai*30.,1.E-8)
end where
f1=(1.+ff)/(ff+cnveg%rsmin*cnveg%lai/5000.)
f2=max(0.5*(if_sfc-if_swilt)/max(rdhyd%soilwater-if_swilt,1.E-9),1.)
f3=max(1.-0.00025*(vegqsat-d_canyonmix)*d_sigd/0.622,0.5) ! increased limit from 0.05 to 0.5 following Mk3.6    
f4=max(1.-0.0016*(298.-urbtemp-d_canyontemp)**2,0.05)     ! 0.2 in Mk3.6
res=max(30.,cnveg%rsmin*f1*f2/(f3*f4))

! solve for vegetation and snow sensible heat fluxes
fg_vegc=aircp*a_rho*(cnveg%temp-d_canyontemp)*acond_vegc
fg_rdsn=aircp*a_rho*(rdsntemp-d_canyontemp)*acond_rdsn

! calculate longwave radiation for vegetation and snow
effvegc=cnveg%emiss*(a_rg*d_cra+(d_netrad*d_crr-sbconst*(cnveg%temp+urbtemp)**4)                          &
                  +sbconst*if_wall%emiss*((walle%nodetemp(:,0)+urbtemp)**4+(wallw%nodetemp(:,0)+urbtemp)**4)*d_crw)
rg_vegc=effvegc-lwflux_walle_vegc-lwflux_wallw_vegc
effrdsn=snowemiss*(a_rg*d_cra+(d_netrad*d_crr-sbconst*(rdsntemp+urbtemp)**4)                              &
                  +sbconst*if_wall%emiss*((walle%nodetemp(:,0)+urbtemp)**4+(wallw%nodetemp(:,0)+urbtemp)**4)*d_crw)
rg_rdsn=effrdsn-lwflux_walle_rdsn-lwflux_wallw_rdsn

! estimate snow melt
rdsnmelt=min(max(0.,rdsntemp+(urbtemp-273.16))*icecp*rdhyd%snow/(ddt*lf),rdhyd%snow/ddt)

! calculate transpiration and evaporation of in-canyon vegetation
d_tranc=lv*min(max((1.-dumvegdelta)*a_rho*(vegqsat-d_canyonmix)/(1./max(acond_vegc,1.e-10)+res),0.), &
               max((rdhyd%soilwater-if_swilt)*d_totdepth*waterden/(d_c1c*ddt),0.))
d_evapc=lv*min(dumvegdelta*a_rho*(vegqsat-d_canyonmix)*acond_vegc,rdhyd%leafwater/ddt+a_rnd)
eg_vegc=d_evapc+d_tranc

! calculate canyon snow latent heat and ground fluxes
snevap=min(a_rho*max(0.,rdsnqsat-d_canyonmix)*acond_rdsn,rdhyd%snow/ddt+a_snd-rdsnmelt)
eg_rdsn=lv*snevap
rdsnmelt=rdsnmelt+snevap
gardsn=(rdsntemp-road%nodetemp(:,0))/ldratio ! use road temperature to represent canyon bottom surface temperature
                                             ! (i.e., we have ommited soil under vegetation temperature)

! vegetation energy budget error term
evct(:,1) = sg_vegc+rg_vegc-fg_vegc-eg_vegc

! road snow energy balance error term
evct(:,2) = sg_rdsn+rg_rdsn-fg_rdsn-eg_rdsn-lf*rdsnmelt-gardsn*(1.-cnveg%sigma)

return
end subroutine canyonflux

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve for roof fluxes

subroutine solveroof(sg_rfsn,rg_rfsn,fg_rfsn,eg_rfsn,garfsn,rfsnmelt,rfsntemp,acond_rfsn,d_rfsndelta, &
                     sg_vegr,rg_vegr,fg_vegr,eg_vegr,acond_vegr,d_vegdeltar,                          &
                     sg_roof,rg_roof,eg_roof,acond_roof,d_roofdelta,                                  &
                     a_rg,a_umag,a_rho,a_rnd,a_snd,d_tempr,d_mixrr,d_rfdzmin,d_tranr,d_evapr,d_c1r,   &
                     d_sigr,ddt,if_roof,rfhyd,rfveg,roof,if_rfvegdepth,if_sfc,if_swilt,ufull)

implicit none

integer, intent(in) :: ufull
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
real, dimension(ufull) :: skintemp
real, dimension(ufull,2) :: oldval,evctx,evctveg
type(facetparams), intent(in) :: if_roof
type(hydrodata), intent(in) :: rfhyd
type(vegdata), intent(inout) :: rfveg
type(facetdata), intent(inout) :: roof
real, dimension(ufull), intent(in) :: if_rfvegdepth, if_sfc, if_swilt

if ( conductmeth==0 ) then
  roof%nodetemp(:,0) = roof%nodetemp(:,1) ! 1st estimate for calculating roof snow temp
end if

lzomroof=log(d_rfdzmin/zoroof)
lzohroof=2.3+lzomroof
call getqsat(qsatr,roof%nodetemp(:,0),d_sigr)
dts=roof%nodetemp(:,0) + (roof%nodetemp(:,0)+urbtemp)*0.61*d_roofdelta*qsatr
dtt=d_tempr + (d_tempr+urbtemp)*0.61*d_mixrr
! Assume zot=0.1*zom (i.e., Kanda et al 2007, small experiment)
call getinvres(acond_roof,cdroof,z_on_l,lzohroof,lzomroof,d_rfdzmin,dts,dtt,a_umag,1)

! update green roof and snow temperature
rfveg%temp=d_tempr
rfsntemp  =roof%nodetemp(:,0)
rg_vegr = if_roof%emiss*(a_rg-sbconst*(roof%nodetemp(:,0)+urbtemp)**4) ! 1st guess
rg_rfsn = if_roof%emiss*(a_rg-sbconst*(roof%nodetemp(:,0)+urbtemp)**4) ! 1st guess
eg_vegr = 0.
eg_rfsn = 0.
rfsnmelt = 0.
garfsn = 0.
d_tranr = 0.
d_evapr = 0.
acond_vegr = acond_roof
acond_rfsn = acond_roof
if ( any( d_rfsndelta>0. .or. rfveg%sigma>0. ) ) then
  evctveg = 0.
  oldval(:,1)=rfveg%temp+0.5
  oldval(:,2)=rfsntemp+0.5
  call roofflux(evctveg,rfsntemp,rfsnmelt,garfsn,sg_vegr,rg_vegr,fg_vegr,eg_vegr,acond_vegr, &
                sg_rfsn,rg_rfsn,fg_rfsn,eg_rfsn,acond_rfsn,a_rg,a_umag,a_rho,a_rnd,a_snd,    &
                d_tempr,d_mixrr,d_rfdzmin,d_tranr,d_evapr,d_c1r,d_sigr,d_vegdeltar,          &
                d_rfsndelta,ddt,if_rfvegdepth,if_roof,if_sfc,if_swilt,rfhyd,rfveg,roof,ufull)
  ! turn off roof snow and roof vegetation if they are not needed
  where ( rfveg%sigma>0. )
    rfveg%temp=rfveg%temp-0.5
  end where
  where ( d_rfsndelta>0. )
    rfsntemp  =rfsntemp-0.5
  end where
  do k=1,nfgits
    evctx=evctveg
    call roofflux(evctveg,rfsntemp,rfsnmelt,garfsn,sg_vegr,rg_vegr,fg_vegr,eg_vegr,acond_vegr, &
                  sg_rfsn,rg_rfsn,fg_rfsn,eg_rfsn,acond_rfsn,a_rg,a_umag,a_rho,a_rnd,a_snd,    &
                  d_tempr,d_mixrr,d_rfdzmin,d_tranr,d_evapr,d_c1r,d_sigr,d_vegdeltar,          &
                  d_rfsndelta,ddt,if_rfvegdepth,if_roof,if_sfc,if_swilt,rfhyd,rfveg,roof,ufull)
    evctx=evctveg-evctx
    where ( abs(evctx(:,1))>tol .and. rfveg%sigma>0. )
      newval=rfveg%temp-alpha*evctveg(:,1)*(rfveg%temp-oldval(:,1))/evctx(:,1)
      oldval(:,1)=rfveg%temp
      rfveg%temp=newval
    end where
    where ( abs(evctx(:,2))>tol .and. d_rfsndelta>0. )
      newval=min(rfsntemp-alpha*evctveg(:,2)*(rfsntemp-oldval(:,2))/evctx(:,2), 300.-urbtemp)
      oldval(:,2)=rfsntemp
      rfsntemp=newval
    end where
  end do
end if
fg_vegr=sg_vegr+rg_vegr-eg_vegr
fg_rfsn=sg_rfsn+rg_rfsn-eg_rfsn-lf*rfsnmelt-garfsn*(1.-rfveg%sigma)

! estimate roof latent heat flux (approx roof_skintemp with roof%nodetemp(:,1))
where ( qsatr<d_mixrr )
  ! dew
  eg_roof=lv*a_rho*(qsatr-d_mixrr)*acond_roof
elsewhere
  ! evaporation
  aa=rfhyd%surfwater/ddt+a_rnd+rfsnmelt
  eg_roof=lv*min(a_rho*d_roofdelta*(qsatr-d_mixrr)*acond_roof,aa)
end where

if ( conductmeth==0 ) then     
  ! estimate roof skin temperature
  ! Write roof energy budget
  !     Solar_net + Longwave_net - Sensible flux - Latent flux - Conduction = 0
  ! or 
  !     sg_roof + a_rg - if_roof%emiss*sbconst*Tskin**4 - aircp*a_rho*(Tskin-d_tempr) &
  !     -eg_roof - (Tskin-roof%nodetemp(:,1))/ldrratio = 0
  ! as a quartic equation
  !      aa*Tskin^4 + dd*Tskin + ee = 0
  ! and solve for Tskin
  ldratio=0.5*(if_roof%depth(:,1)/if_roof%lambda(:,1))
  aa=if_roof%emiss*sbconst
  dd=aircp*a_rho*acond_roof+1./ldratio
  ee=-sg_roof-if_roof%emiss*a_rg-aircp*a_rho*acond_roof*(d_tempr+urbtemp)-(roof%nodetemp(:,1)+urbtemp)/ldratio+eg_roof
  call solvequartic(skintemp,aa,dd,ee) ! This the 2nd estimate of Tskin to be updated in solvetridiag
  roof%nodetemp(:,0) = skintemp - urbtemp
end if

! calculate net roof longwave radiation
! (sensible heat flux will be updated in solvetridiag)
rg_roof=if_roof%emiss*(a_rg-sbconst*(roof%nodetemp(:,0)+urbtemp)**4)

return
end subroutine solveroof

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve for green roof and snow fluxes

subroutine roofflux(evct,rfsntemp,rfsnmelt,garfsn,sg_vegr,rg_vegr,fg_vegr,eg_vegr,acond_vegr,   &
                    sg_rfsn,rg_rfsn,fg_rfsn,eg_rfsn,acond_rfsn,a_rg,a_umag,a_rho,a_rnd,a_snd,   &
                    d_tempr,d_mixrr,d_rfdzmin,d_tranr,d_evapr,d_c1r,d_sigr,d_vegdeltar,         &
                    d_rfsndelta,ddt,if_rfvegdepth,if_roof,if_sfc,if_swilt,rfhyd,rfveg,roof,ufull)

implicit none

integer, intent(in) :: ufull
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
real, dimension(ufull) :: lzotdum, snevap
real, dimension(ufull), intent(in) :: if_rfvegdepth, if_sfc, if_swilt
type(facetparams), intent(in) :: if_roof
type(hydrodata), intent(in) :: rfhyd
type(vegdata), intent(in) :: rfveg
type(facetdata), intent(in) :: roof

call getqsat(vegqsat,rfveg%temp,d_sigr)
where ( vegqsat<d_mixrr )
  dumvegdelta = 1.
elsewhere
  dumvegdelta = d_vegdeltar
end where

! transpiration terms (developed by Eva in CCAM sflux.f and CSIRO9)
where ( rfveg%zo<0.5 )
  ff = 1.1*sg_vegr/max(rfveg%lai*150.,1.E-8)
elsewhere
  ff = 1.1*sg_vegr/max(rfveg%lai*30.,1.E-8)
end where
f1 = (1.+ff)/(ff+rfveg%rsmin*rfveg%lai/5000.)
f2 = max(0.5*(if_sfc-if_swilt)/max(rfhyd%soilwater-if_swilt,1.E-9),1.)
f3 = max(1.-.00025*(vegqsat-d_mixrr)*d_sigr/0.622,0.5)
f4 = max(1.-0.0016*((298.-urbtemp)-d_tempr)**2,0.05)
res = max(30.,rfveg%rsmin*f1*f2/(f3*f4))

vwetfac = max(min((rfhyd%soilwater-if_swilt)/(if_sfc-if_swilt),1.),0.) ! veg wetfac (see sflux.f or cable_canopy.f90)
vwetfac = (1.-dumvegdelta)*vwetfac+dumvegdelta
lzomvegr = log(d_rfdzmin/rfveg%zo)
! xe is a dummy variable for lzohvegr
lzohvegr = 2.3+lzomvegr
dts = rfveg%temp + (rfveg%temp+urbtemp)*0.61*vegqsat*vwetfac
dtt = d_tempr + (d_tempr+urbtemp)*0.61*d_mixrr
! Assume zot=0.1*zom (i.e., Kanda et al 2007, small experiment)
call getinvres(acond_vegr,cdvegr,z_on_l,lzohvegr,lzomvegr,d_rfdzmin,dts,dtt,a_umag,1)
! acond_vegr is multiplied by a_umag

where ( rfveg%sigma>0. )
  ! longwave radiation    
  rg_vegr=rfveg%emiss*(a_rg-sbconst*(rfveg%temp+urbtemp)**4)
  
  ! sensible heat flux
  fg_vegr=aircp*a_rho*(rfveg%temp-d_tempr)*acond_vegr

  ! calculate transpiration and evaporation of in-canyon vegetation
  d_tranr=lv*min(max((1.-dumvegdelta)*a_rho*(vegqsat-d_mixrr)/(1./acond_vegr+res),0.), &
                 max((rfhyd%soilwater-if_swilt)*if_rfvegdepth*waterden/(d_c1r*ddt),0.))
  d_evapr=lv*min(dumvegdelta*a_rho*(vegqsat-d_mixrr)*acond_vegr,rfhyd%leafwater/ddt+a_rnd)
  eg_vegr=d_evapr+d_tranr
  
  ! balance green roof energy budget
  evct(:,1)=sg_vegr+rg_vegr-fg_vegr-eg_vegr
end where


! snow conductance
sndepth=rfhyd%snow*waterden/rfhyd%den
snlambda=icelambda*(rfhyd%den/waterden)**1.88
ldratio=0.5*(sndepth/snlambda+if_roof%depth(:,1)/if_roof%lambda(:,1))

! Update roof snow energy budget
lzosnow=log(d_rfdzmin/zosnow)
call getqsat(rfsnqsat,rfsntemp,d_sigr)
lzotdum=2.3+lzosnow
dts=rfsntemp + (rfsntemp+urbtemp)*0.61*rfsnqsat
call getinvres(acond_rfsn,cdrfsn,z_on_l,lzotdum,lzosnow,d_rfdzmin,dts,dtt,a_umag,1)
! acond_rfsn is multiplied by a_umag

where ( d_rfsndelta>0. )
  rfsnmelt=min(max(0.,rfsntemp+(urbtemp-273.16))*icecp*rfhyd%snow/(ddt*lf),rfhyd%snow/ddt)
  rg_rfsn=snowemiss*(a_rg-sbconst*(rfsntemp+urbtemp)**4)
  fg_rfsn=aircp*a_rho*(rfsntemp-d_tempr)*acond_rfsn
  snevap=min(a_rho*max(0.,rfsnqsat-d_mixrr)*acond_rfsn,rfhyd%snow/ddt+a_snd-rfsnmelt)
  eg_rfsn=lv*snevap
  rfsnmelt=rfsnmelt+snevap
  garfsn=(rfsntemp-roof%nodetemp(:,0))/ldratio
  
  ! balance snow energy budget
  evct(:,2)=sg_rfsn+rg_rfsn-fg_rfsn-eg_rfsn-lf*rfsnmelt-garfsn*(1.-rfveg%sigma)
end where

return
end subroutine roofflux

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define traffic flux weights during the diurnal cycle

subroutine getdiurnal(if_ctime,icyc_traffic,icyc_basedemand,icyc_proportion,icyc_translation)                  

implicit none

real, dimension(:), intent(in) :: if_ctime
real, dimension(:), intent(out) :: icyc_traffic,icyc_basedemand,icyc_proportion,icyc_translation
real, dimension(size(if_ctime)) :: real_p
integer, dimension(size(if_ctime)) :: int_p

! traffic diurnal cycle weights approximated from Coutts et al (2007)
real, dimension(25), parameter :: trafcycle = (/ 0.2, 0.1, 0.1, 0.2, 0.5, 1.1, 1.7, 1.4, 1.1, 1.1, 1.2, 1.3, &
                                                 1.3, 1.4, 1.6, 1.8, 2.0, 1.6, 1.2, 0.9, 0.8, 0.7, 0.5, 0.2, &
                                                 0.2 /)
! base electricity demand cycle weights approximated from Thatcher (2006) for Victoria
real, dimension(25), parameter :: basecycle = (/ 0.9, 0.9, 0.9, 0.9, 0.8, 0.8, 0.9, 1.0, 1.1, 1.1, 1.1, 1.1, &
                                                 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.0, 1.0, 1.0, 1.0, 0.9, &
                                                 0.9 /)
! proportion of heating/cooling appliances in use  approximated from Thatcher (2006)
real, dimension(25), parameter :: propcycle = (/ 0.7, 0.7, 0.6, 0.6, 0.5, 0.5, 0.6, 0.8, 1.0, 1.1, 1.0, 1.0, &
                                                 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6, 1.6, 1.4, 1.3, 1.1, 0.9, &
                                                 0.7 /)
! base temperature translation cycle approximated from Thatcher (2006)
real, dimension(25), parameter :: trancycle = (/ -1.1, -0.9, -1.8, -2.5, -3.0, -2.8, -1.4, 0.3, 1.6, 2.2, 2.5, 2.4, &
                                                  2.1, 1.7, 1.1, 0.4, 0.1, 1.1, 1.4, 0.6, -0.3, -1.0, -1.2, -1.1,   &
                                                 -1.1 /)

int_p=int(24.*if_ctime)
real_p=24.*if_ctime-real(int_p)
where (int_p<1) int_p=int_p+24

icyc_traffic     = ((1.-real_p)*trafcycle(int_p)+real_p*trafcycle(int_p+1))
icyc_basedemand  = ((1.-real_p)*basecycle(int_p)+real_p*basecycle(int_p+1))
icyc_proportion  = ((1.-real_p)*propcycle(int_p)+real_p*propcycle(int_p+1))
icyc_translation = ((1.-real_p)*trancycle(int_p)+real_p*trancycle(int_p+1))

return
end subroutine getdiurnal


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate in-canyon wind speed for walls and road
! This version allows the eddy size to change with canyon orientation
! which requires a numerical solution to the integral

subroutine getincanwind(ueast,uwest,ufloor,a_udir,z0,if_bldheight,if_coeffbldheight,if_hwratio,ufull)

implicit none

integer, intent(in) :: ufull
real, dimension(ufull), intent(out) :: ueast,uwest,ufloor
real, dimension(ufull), intent(in) :: z0
real, dimension(ufull) :: a,b,wsuma,wsumb,fsum
real, dimension(ufull) :: theta1,wdir,h,w
real, dimension(ufull), intent(in) :: a_udir
!global
real, dimension(ufull), intent(in) :: if_bldheight
real, dimension(ufull), intent(in) :: if_coeffbldheight
real, dimension(ufull), intent(in) :: if_hwratio
!

! rotate wind direction so that all cases are between 0 and pi
! walls are fliped at the end of the subroutine to account for additional pi rotation
where (a_udir>=0.)
  wdir=a_udir
elsewhere
  wdir=a_udir+pi
endwhere

h=if_bldheight*if_coeffbldheight
w=if_bldheight/if_hwratio

theta1=asin(min(w/(3.*h),1.))
wsuma=0.
wsumb=0.
fsum=0.  ! floor

! integrate jet on road, venting side (A)
a=0.
b=max(0.,wdir-pi+theta1)
call integratewind(wsuma,wsumb,fsum,a,b,h,w,wdir,z0,0,ufull)

! integrate jet on wall, venting side
a=max(0.,wdir-pi+theta1)
b=max(0.,wdir-theta1)
call integratewind(wsuma,wsumb,fsum,a,b,h,w,wdir,z0,1,ufull)

! integrate jet on road, venting side (B)
a=max(0.,wdir-theta1)
b=wdir
call integratewind(wsuma,wsumb,fsum,a,b,h,w,wdir,z0,0,ufull)

! integrate jet on road, recirculation side (A)
a=wdir
b=min(pi,wdir+theta1)
call integratewind(wsumb,wsuma,fsum,a,b,h,w,wdir,z0,0,ufull)

! integrate jet on wall, recirculation side
a=min(pi,wdir+theta1)
b=min(pi,wdir+pi-theta1)
call integratewind(wsumb,wsuma,fsum,a,b,h,w,wdir,z0,1,ufull)

! integrate jet on road, recirculation side (B)
a=min(pi,wdir+pi-theta1)
b=pi
call integratewind(wsumb,wsuma,fsum,a,b,h,w,wdir,z0,0,ufull)

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

subroutine getincanwindb(ueast,uwest,ufloor,a_udir,z0,if_bldheight,if_coeffbldheight,if_hwratio,ufull)

implicit none

integer, intent(in) :: ufull
real, dimension(ufull), intent(out) :: ueast,uwest,ufloor
real, dimension(ufull), intent(in) :: z0
real, dimension(ufull) :: wsuma,wsumb,fsum
real, dimension(ufull) :: theta1,wdir,h,w
real, dimension(ufull) :: dufa,dura,duva,ntheta
real, dimension(ufull) :: dufb,durb,duvb
real, dimension(ufull), intent(in) :: a_udir
!global
real, dimension(ufull), intent(in) :: if_bldheight
real, dimension(ufull), intent(in) :: if_coeffbldheight
real, dimension(ufull), intent(in) :: if_hwratio
!

! rotate wind direction so that all cases are between 0 and pi
! walls are fliped at the end of the subroutine to account for additional pi rotation
where (a_udir>=0.)
  wdir=a_udir
elsewhere
  wdir=a_udir+pi
endwhere

h=if_bldheight*if_coeffbldheight
w=if_bldheight/if_hwratio

theta1=acos(min(w/(3.*h),1.))

call winda(dufa,dura,duva,h,w,z0,ufull) ! jet on road
call windb(dufb,durb,duvb,h,w,z0,ufull) ! jet on wall
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

subroutine integratewind(wsuma,wsumb,fsum,a,b,h,w,wdir,z0,mode,ufull)

implicit none

integer, intent(in) :: ufull
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
        call winda(duf,dur,duv,h,nw,z0,ufull)
        wsuma=wsuma+dur*st*dtheta
        wsumb=wsumb+duv*st*dtheta
        fsum=fsum+duf*st*dtheta
      end do
    case(1) ! jet on wall
      do n=1,ntot
        theta=dtheta*(real(n)-0.5)+a
        st=abs(sin(theta-wdir))
        nw=min(w/max(st,1.E-9),3.*h)
        call windb(duf,dur,duv,h,nw,z0,ufull)
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

subroutine winda(uf,ur,uv,h,w,z0,ufull)

implicit none

integer, intent(in) :: ufull
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

subroutine windb(uf,ur,uv,h,win,z0,ufull)

implicit none

integer, intent(in) :: ufull
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
subroutine gettopu(d_topu,a_umag,z_on_l,if_bldheight,ip_cduv,ip_cndzmin,if_hwratio,ufull)
      
implicit none

integer, intent(in) :: ufull
real, dimension(ufull), intent(in) :: z_on_l
real, dimension(ufull) :: z0_on_l,bldheight
real, dimension(ufull) :: pm0,pm1,integralm
real, dimension(ufull) :: ustar,neutral
real, dimension(ufull), intent(inout) :: d_topu
real, dimension(ufull), intent(in) :: a_umag
real, dimension(ufull), intent(in) :: if_bldheight
real, dimension(ufull), intent(inout) :: ip_cduv,ip_cndzmin
!global
real, dimension(ufull), intent(in) :: if_hwratio
!

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

subroutine getc1(dc1,ufull)

implicit none

integer, intent(in) :: ufull
real, dimension(ufull), intent(out) :: dc1

!n=min(max(moist/if_ssat,0.218),1.)
!dc1=(1.78*n+0.253)/(2.96*n-0.581)

dc1=1.478 ! simplify water conservation

return
end subroutine getc1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate screen diagnostics

subroutine scrncalc(a_mixr,a_umag,a_temp,u_ts,d_tempc,d_rdsndelta,d_roaddelta,d_vegdeltac,d_sigd,smixr,rdsntemp,zonet, &
                    cnveg,if_bldheight,if_sfc,if_swilt,p_cndzmin,p_lzoh,p_lzom,p_qscrn,p_tscrn,p_u10,p_uscrn,rdhyd,    &
                    road,ufull)
      
implicit none

integer, intent(in) :: ufull
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
real, dimension(ufull), intent(in) :: if_bldheight, if_sfc, if_swilt
real, dimension(ufull), intent(in) :: p_cndzmin, p_lzoh, p_lzom
real, dimension(ufull), intent(inout) :: p_qscrn, p_tscrn, p_u10, p_uscrn
type(vegdata), intent(in) :: cnveg
type(hydrodata), intent(in) :: rdhyd
type(facetdata), intent(in) :: road

real, parameter :: z0  = 1.5
real, parameter :: z10 = 10.

select case(scrnmeth)
  case(0) ! estimate screen diagnostics (slab at displacement height approach)
    thetav=d_tempc + (d_tempc+urbtemp)*0.61*a_mixr
    sthetav=u_ts + (u_ts+urbtemp)*0.61*smixr
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
                    -log((1.+1./pm1**2)/(1.+1./pm0**2))     &
                    +2.*(atan(1./pm1)-atan(1./pm0))
      integralm10 = neutral10-2.*log((1.+1./pm1)/(1.+1./pm10)) &
                    -log((1.+1./pm1**2)/(1.+1./pm10**2))       &
                    +2.*(atan(1./pm1)-atan(1./pm10))     
    elsewhere
      !-------Beljaars and Holtslag (1991) heat function
      ph0  = -((1.+(2./3.)*a_1*z0_on_l)**1.5 &
             +b_1*(z0_on_l-(c_1/d_1))        &
             *exp(-d_1*z0_on_l)+b_1*c_1/d_1-1.)
      ph1  = -((1.+(2./3.)*a_1*z_on_l)**1.5 &
             +b_1*(z_on_l-(c_1/d_1))        &
             *exp(-d_1*z_on_l)+b_1*c_1/d_1-1.)
      pm0 = -(a_1*z0_on_l+b_1*(z0_on_l-(c_1/d_1))*exp(-d_1*z0_on_l)+b_1*c_1/d_1)
      pm10 = -(a_1*z10_on_l+b_1*(z10_on_l-(c_1/d_1))*exp(-d_1*z10_on_l)+b_1*c_1/d_1)
      pm1  = -(a_1*z_on_l+b_1*(z_on_l-(c_1/d_1))*exp(-d_1*z_on_l)+b_1*c_1/d_1)
      integralh   = neutral-(ph1-ph0)
      integralm   = neutral-(pm1-pm0)
      integralm10 = neutral10-(pm1-pm10)
    endwhere
    p_tscrn = a_temp - tstar*integralh/vkar
    p_qscrn = a_mixr - qstar*integralh/vkar
    p_uscrn = max(a_umag-ustar*integralm/vkar,0.)
    p_u10   = max(a_umag-ustar*integralm10/vkar,0.)
    
  case(1) ! estimate screen diagnostics (two step canopy approach)
    thetav=d_tempc + (d_tempc+urbtemp)*0.61*a_mixr
    sthetav=u_ts + (u_ts+urbtemp)*0.61*smixr
    lna=p_lzoh-p_lzom
    call dyerhicks(integralh,z_on_l,cd,thetavstar,thetav,sthetav,a_umag,p_cndzmin,p_lzom,lna,4)
    ustar=sqrt(cd)*a_umag
    qstar=vkar*(a_mixr-smixr)/integralh
    tts=vkar*(thetav-sthetav)/integralh
    tstar=vkar*(a_temp-u_ts)/integralh
    
    z0_on_l  = if_bldheight*(1.-refheight)*z_on_l/p_cndzmin ! calculate at canyon top
    z10_on_l = max(z10-if_bldheight*refheight,1.)*z_on_l/p_cndzmin
    z0_on_l  = min(z0_on_l,10.)
    z10_on_l = min(z10_on_l,10.)
    neutral   = log(p_cndzmin/(if_bldheight*(1.-refheight)))
    neutral10 = log(p_cndzmin/max(z10-if_bldheight*refheight,1.))
    where (z_on_l<0.)
      ph0     = (1.-16.*z0_on_l)**(-0.50)
      ph1     = (1.-16.*z_on_l)**(-0.50)
      pm0     = (1.-16.*z0_on_l)**(-0.25)
      pm10    = (1.-16.*z10_on_l)**(-0.25)
      pm1     = (1.-16.*z_on_l)**(-0.25)
      integralh = neutral-2.*log((1.+1./ph1)/(1.+1./ph0))
      integralm = neutral-2.*log((1.+1./pm1)/(1.+1./pm0)) &
                    -log((1.+1./pm1**2)/(1.+1./pm0**2))   &
                    +2.*(atan(1./pm1)-atan(1./pm0))
      integralm10 = neutral10-2.*log((1.+1./pm1)/(1.+1./pm10)) &
                    -log((1.+1./pm1**2)/(1.+1./pm10**2))       &
                    +2.*(atan(1./pm1)-atan(1./pm10))     
    elsewhere
      !-------Beljaars and Holtslag (1991) heat function
      ph0  = -((1.+(2./3.)*a_1*z0_on_l)**1.5 &
             +b_1*(z0_on_l-(c_1/d_1))        &
             *exp(-d_1*z0_on_l)+b_1*c_1/d_1-1.)
      ph1  = -((1.+(2./3.)*a_1*z_on_l)**1.5 &
             +b_1*(z_on_l-(c_1/d_1))        &
             *exp(-d_1*z_on_l)+b_1*c_1/d_1-1.)
      pm0 = -(a_1*z0_on_l+b_1*(z0_on_l-(c_1/d_1))*exp(-d_1*z0_on_l)+b_1*c_1/d_1)
      pm10 = -(a_1*z10_on_l+b_1*(z10_on_l-(c_1/d_1))*exp(-d_1*z10_on_l)+b_1*c_1/d_1)
      pm1  = -(a_1*z_on_l+b_1*(z_on_l-(c_1/d_1))*exp(-d_1*z_on_l)+b_1*c_1/d_1)
      integralh   = neutral-(ph1-ph0)
      integralm   = neutral-(pm1-pm0)
      integralm10 = neutral10-(pm1-pm10)
    endwhere
    ttop = thetav - tts*integralh/vkar
    tetp = a_temp - tstar*integralh/vkar
    qtop = a_mixr - qstar*integralh/vkar
    utop = a_umag - ustar*integralm/vkar

    where (if_bldheight<=z10) ! above canyon
      p_u10=max(a_umag-ustar*integralm10/vkar,0.)
    end where

    ! assume standard stability functions hold for urban canyon (needs more work)
    tsurf = d_rdsndelta*rdsntemp+(1.-d_rdsndelta)*((1.-cnveg%sigma)*road%nodetemp(:,0)+cnveg%sigma*cnveg%temp)
    n=max(min((rdhyd%soilwater-if_swilt)/(if_sfc-if_swilt),1.),0.)
    wf = (1.-d_rdsndelta)*((1.-cnveg%sigma)*d_roaddelta+cnveg%sigma*((1.-d_vegdeltac)*n+d_vegdeltac))
    call getqsat(qsurf,tsurf,d_sigd)
    qsurf=qsurf*wf
    n=log(if_bldheight/zonet)
    
    thetav=ttop + (ttop+urbtemp)*0.61*qtop
    sthetav=tsurf + (tsurf+urbtemp)*0.61*qsurf
    lna=2.3
    call dyerhicks(integralh,z_on_l,cd,thetavstar,thetav,sthetav,utop,if_bldheight,n,lna,1)
    ustar=sqrt(cd)*utop
    tstar=vkar*(tetp-tsurf)/integralh
    qstar=vkar*(qtop-qsurf)/integralh
    
    z0_on_l   = z0*z_on_l/if_bldheight
    z10_on_l  = max(z10,if_bldheight)*z_on_l/if_bldheight
    z0_on_l   = min(z0_on_l,10.)
    z10_on_l  = min(z10_on_l,10.)
    neutral   = log(if_bldheight/z0)
    neutral10 = log(if_bldheight/max(z10,if_bldheight))
    where (z_on_l<0.)
      ph0     = (1.-16.*z0_on_l)**(-0.50)
      ph1     = (1.-16.*z_on_l)**(-0.50)
      pm0     = (1.-16.*z0_on_l)**(-0.25)
      pm10    = (1.-16.*z10_on_l)**(-0.25)
      pm1     = (1.-16.*z_on_l)**(-0.25)
      integralh = neutral-2.*log((1.+1./ph1)/(1.+1./ph0))
      integralm = neutral-2.*log((1.+1./pm1)/(1.+1./pm0)) &
                    -log((1.+1./pm1**2)/(1.+1./pm0**2))   &
                    +2.*(atan(1./pm1)-atan(1./pm0))
      integralm10 = neutral10-2.*log((1.+1./pm1)/(1.+1./pm10)) &
                    -log((1.+1./pm1**2)/(1.+1./pm10**2))       &
                    +2.*(atan(1./pm1)-atan(1./pm10))     
    elsewhere
      !-------Beljaars and Holtslag (1991) heat function
      ph0  = -((1.+(2./3.)*a_1*z0_on_l)**1.5 &
             +b_1*(z0_on_l-(c_1/d_1))        &
             *exp(-d_1*z0_on_l)+b_1*c_1/d_1-1.)
      ph1  = -((1.+(2./3.)*a_1*z_on_l)**1.5 &
             +b_1*(z_on_l-(c_1/d_1))        &
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
    where (if_bldheight>z10) ! within canyon
      p_u10 = max(utop-ustar*integralm10/vkar,0.)
    end where

  case(2) ! calculate screen diagnostics from canyon only
    tsurf=d_rdsndelta*rdsntemp+(1.-d_rdsndelta)*((1.-cnveg%sigma)*road%nodetemp(:,0)+cnveg%sigma*cnveg%temp)
    n=max(min((rdhyd%soilwater-if_swilt)/(if_sfc-if_swilt),1.),0.)
    wf=(1.-d_rdsndelta)*((1.-cnveg%sigma)*d_roaddelta+cnveg%sigma*((1.-d_vegdeltac)*n+d_vegdeltac))
    call getqsat(qsurf,tsurf,d_sigd)
    qsurf=qsurf*wf
    n=log(if_bldheight/zonet)

    thetav=d_tempc + (d_tempc+urbtemp)*0.61*a_mixr
    sthetav=tsurf + (tsurf+urbtemp)*0.61*qsurf
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
                    -log((1.+1./pm1**2)/(1.+1./pm0**2))     &
                    +2.*(atan(1./pm1)-atan(1./pm0))
      integralm10 = neutral10-2.*log((1.+1./pm1)/(1.+1./pm10)) &
                    -log((1.+1./pm1**2)/(1.+1./pm10**2))       &
                    +2.*(atan(1./pm1)-atan(1./pm10))     
    elsewhere
      !-------Beljaars and Holtslag (1991) heat function
      ph0  = -((1.+(2./3.)*a_1*z0_on_l)**1.5 &
             +b_1*(z0_on_l-(c_1/d_1))        &
             *exp(-d_1*z0_on_l)+b_1*c_1/d_1-1.)
      ph1  = -((1.+(2./3.)*a_1*z_on_l)**1.5 &
             +b_1*(z_on_l-(c_1/d_1))        &
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

d0=12.*a*e
d1=27.*a*d**2
qq=(0.5*(d1+sqrt(d1**2-4.*d0**3)))**(1./3.)
s=0.5*sqrt((qq+d0/qq)/(3.*a))
q=d/a
t1=-s+0.5*sqrt(-4*s**2+q/s)
!t2=-s-0.5*sqrt(-4*s**2+q/s)
!t3=s+0.5*sqrt(-4*s**2-q/s)
!t4=s-0.5*sqrt(-4*s**2-q/s)

x=t1

return
end subroutine solvequartic


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Disables aTEB so subroutine calls have no effect

subroutine atebdisable(diag)

implicit none

integer, intent(in) :: diag

if ( diag>=1 ) write(6,*) "Disable aTEB"

ufull_g = 0

return
end subroutine atebdisable

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The following subroutines are required for internal conditions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This subroutine calculates net longwave radiation flux (flux_rg) at each surface
! longwave flux is temperature dependent, so this subroutine should be run at each timestep
subroutine internal_lwflux(rgint_slab,rgint_wallw,rgint_roof,rgint_walle, &
                           if_bldheight,if_bldwidth,int_psi,int_viewf,    &
                           roof,slab,walle,wallw,ufull)

implicit none
integer, intent(in) :: ufull
real(kind=8), dimension(ufull,4) :: skintemp  ! floor, wall, ceiling, wall temperature array
real(kind=8), dimension(ufull,4) :: epsil     ! floor, wall, ceiling, wall emissivity array
real(kind=8), dimension(ufull,4) :: radnet    ! net flux density on ith surface
real(kind=8), dimension(ufull,4) :: rad       ! net leaving flux density (B) on ith surface
real(kind=8), dimension(ufull)   :: radtot    ! net leaving flux density (B) on ith surface
real(kind=8), dimension(ufull)   :: sum_int_viewf_rad
real, dimension(ufull), intent(out) :: rgint_slab,rgint_wallw,rgint_roof,rgint_walle
integer :: j
real, dimension(ufull), intent(in) :: if_bldheight, if_bldwidth
real(kind=8), dimension(ufull,4,4), intent(in) :: int_psi, int_viewf
type(facetdata), intent(in) :: roof, slab, walle, wallw
!

rad = 0.
radnet = 0.

!epsil = reshape((/(if_slab%emiss,if_wall%emiss,if_roof%emiss,if_wall%emiss, & 
!                    i=1,ufull)/), (/ufull,4/))
epsil = 0.9

skintemp = reshape((/ slab%nodetemp(:,nl),    &
                      wallw%nodetemp(:,nl),   &
                      roof%nodetemp(:,nl),    &
                      walle%nodetemp(:,nl)    &
                   /),(/ufull,4/)) + urbtemp
do j = 1,4
  rad(:,j) = sum(int_psi(:,j,:)*epsil(:,:)*sbconst*skintemp(:,:)**4,dim=2)
end do

do j = 1,4
  sum_int_viewf_rad(:) = sum( int_viewf(:,j,:)*rad(:,:),dim=2)
  ! Harman et al. (2004) Eq. (12)
  where ( abs(epsil(:,j)-1._8)<1.e-20_8 ) ! for black body surface (no reflection)
    radnet(:,j) = epsil(:,j)*sbconst*skintemp(:,j)**4 - sum_int_viewf_rad(:)
  elsewhere               ! for grey body calculation (infinite reflection) 
    radnet(:,j) = (epsil(:,j)*sbconst*skintemp(:,j)**4 - epsil(:,j)*rad(:,j))/(1.-epsil(:,j))
  end where
end do

radtot(:) = abs(if_bldwidth(:)*(radnet(:,1)+radnet(:,3)) + if_bldheight*(radnet(:,2)+radnet(:,4)))

do j = 1,ufull
  if ( radtot(j)>1E-8 ) write(6,*) "error: radiation energy non-closure: ", radtot(j)
end do

rgint_slab  = real(radnet(:,1))
rgint_wallw = real(radnet(:,2))
rgint_roof  = real(radnet(:,3))
rgint_walle = real(radnet(:,4))

return
end subroutine internal_lwflux

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init_lwcoeff
! This subroutine calculates longwave reflection coefficients (int_psi) at each surface
! longwave coefficients do not change, so this subroutine should only be run once
! Infinite reflections per Harman et al., (2004) "Radiative Exchange in Urban Street Canyons"
! Per method in "Radiation Heat Transfer, Sparrow & Cess 1978, Ch 3-3"
! array surface order is: (1) floor; (2) wallw; (3) ceiling; (4) walle

! local variables
real(kind=8), dimension(ufull_g,nl,nl) :: chi
real(kind=8), dimension(nl,nl)       :: krondelta
real(kind=8), dimension(ufull_g)       :: h, w
real(kind=8), dimension(ufull_g,nl)    :: epsil   ! floor, wall, ceiling, wall emissivity array
integer :: i, j
integer :: ierr       ! inverse matrix error flag

krondelta = 0.
chi = 0.
int_psi = 0.
h = f_bldheight
w = f_sigmabld*(f_bldheight/f_hwratio)/(1.-f_sigmabld)

! set int_vfactors
int_viewf(:,1,1) = 0.                                  ! floor to self
int_viewf(:,1,2) = 0.5*(1.+(h/w)-sqrt(1.+(h/w)**2))    ! floor to wallw
int_viewf(:,1,3) = sqrt(1.+(h/w)**2)-(h/w)             ! floor to ceiling
int_viewf(:,1,4) = int_viewf(:,1,2)                    ! floor to walle
int_viewf(:,2,1) = 0.5*(1.+(w/h)-sqrt(1.+(w/h)**2))    ! wallw to floor
int_viewf(:,2,2) = 0.                                  ! wallw to self
int_viewf(:,2,3) = int_viewf(:,2,1)                    ! wallw to ceiling
int_viewf(:,2,4) = sqrt(1.+(w/h)**2)-(w/h)             ! wallw to walle
int_viewf(:,3,1) = int_viewf(:,1,3)                    ! ceiling to floor
int_viewf(:,3,2) = int_viewf(:,1,2)                    ! ceiling to wallw
int_viewf(:,3,3) = 0.                                  ! ceiling to self
int_viewf(:,3,4) = int_viewf(:,1,2)                    ! ceiling walle
int_viewf(:,4,1) = int_viewf(:,2,1)                    ! walle to floor
int_viewf(:,4,2) = int_viewf(:,2,4)                    ! walle to wallw
int_viewf(:,4,3) = int_viewf(:,2,1)                    ! walle to ceiling
int_viewf(:,4,4) = 0                                   ! walle to self

!epsil = reshape((/(if_slab%emiss,if_wall%emiss,if_roof%emiss,if_wall%emiss, & 
!                    i=1,ufull_g)/), (/ufull_g,4/))
epsil = 0.9
do i = 1,nl
  krondelta(i,i) = 1.
end do
do j = 1,nl
  do i = 1,nl
    chi(:,i,j) = krondelta(i,j) - (1.-epsil(:,i))*int_viewf(:,i,j)
  end do
end do

! invert matrix
int_psi = chi
call minverse(int_psi,ierr)

end subroutine init_lwcoeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine init_internal

implicit none

f_bldwidth = f_sigmabld*(f_bldheight/f_hwratio)/(1.-f_sigmabld)
! define number of internal mass floors (based on building height)
select case(intmassmeth)
  case(0) ! no internal mass
    f_intmassn = 0
  case(1) ! one floor of internal mass
    f_intmassn = 1
  case(2) ! dynamic floors of internal mass
    f_intmassn = max((nint(f_bldheight/3.)-1),1)
end select

!write(6,*) 'building height: ', f_bldheight
!write(6,*) 'building width: ' , f_bldwidth
!write(6,*) '# internal floors: ', f_intmassn

end subroutine init_internal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine sets internal surface convective heat transfer coefficients
! Compares temperature of innermost layer temperature with air temperature
! Considers horizontal and vertical orientation
! Based on EnergyPlus: Simple Natural Convection Algorithm [W m^-2 K^-1]

subroutine calc_convcoeff(cvcoeff_roof,cvcoeff_walle,cvcoeff_wallw,cvcoeff_slab, & 
                          cvcoeff_intm1,cvcoeff_intm2,roof,room,slab,ufull)
implicit none

integer, intent(in) :: ufull
real, dimension(ufull), intent(out) :: cvcoeff_roof,cvcoeff_walle,cvcoeff_wallw,cvcoeff_slab
real, dimension(ufull), intent(out) :: cvcoeff_intm1,cvcoeff_intm2
!global
type(facetdata), intent(in) :: roof
type(facetdata), intent(in) :: room
type(facetdata), intent(in) :: slab
! 

cvcoeff_walle = 3.067               ! vertical surface coefficient constant
cvcoeff_wallw = 3.067               ! vertical surface coefficient constant
where ( roof%nodetemp(:,nl)>=room%nodetemp(:,1) )
  cvcoeff_roof(:)=0.948  ! reduced convection
elsewhere
  cvcoeff_roof(:)=4.040  ! enhanced convection  
end where    
cvcoeff_intm1 = 3.076               ! vertical surface coefficient constant
cvcoeff_intm2 = 3.076               ! vertical surface coefficient constant
where (slab%nodetemp(:,nl)<=room%nodetemp(:,1))
  cvcoeff_slab(:)=0.948  ! reduced convection
elsewhere
  cvcoeff_slab(:)=4.040  ! enhanced convection
end where

end subroutine calc_convcoeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calc_newairtemp(int_newairtemp,a_rho,d_canyontemp,d_intgains_bld,      &
                           cvcoeff_roof,cvcoeff_walle,cvcoeff_wallw,cvcoeff_slab, &
                           cvcoeff_intm1,cvcoeff_intm2,ddt,                       &
                           if_ach,if_bldheight,if_bldwidth,if_intmassn,intm,roof, &
                           room,slab,walle,wallw,ufull)

implicit none

integer, intent(in) :: ufull
real, intent(in)    :: ddt
real, dimension(ufull), intent(in) :: cvcoeff_roof,cvcoeff_walle,cvcoeff_wallw,cvcoeff_slab
real, dimension(ufull), intent(in) :: cvcoeff_intm1,cvcoeff_intm2
real, dimension(ufull), intent(in) :: a_rho,d_canyontemp,d_intgains_bld
real, dimension(ufull), intent(out) :: int_newairtemp
real, dimension(ufull) :: rm,rf,we,ww,sl,im1,im2,infl
real, dimension(ufull), intent(in) :: if_ach, if_bldheight, if_bldwidth
integer, dimension(ufull), intent(in) :: if_intmassn
type(facetdata), intent(in) :: intm, roof, room, slab, walle, wallw

rm = a_rho*aircp*if_bldheight/ddt
rf = cvcoeff_roof
we = (if_bldheight/if_bldwidth)*cvcoeff_walle
ww = (if_bldheight/if_bldwidth)*cvcoeff_wallw
sl = cvcoeff_slab
im1 = cvcoeff_intm1*if_intmassn
im2 = cvcoeff_intm2*if_intmassn
infl = if_ach*a_rho*aircp*if_bldheight/3600.

int_newairtemp = (rm*room%nodetemp(:,1)     & ! room temperature
                + rf*roof%nodetemp(:,nl)    & ! roof conduction
                + we*walle%nodetemp(:,nl)   & ! wall conduction east
                + ww*wallw%nodetemp(:,nl)   & ! wall conduction west
                + sl*slab%nodetemp(:,nl)    & ! slab conduction
                + im1*intm%nodetemp(:,0)    & ! mass conduction side 1
                + im2*intm%nodetemp(:,nl)   & ! mass conduction side 2
                + infl*d_canyontemp         & ! infiltration
                + d_intgains_bld            & ! internal gains
                )/ (rm +rf +we +ww +sl +im1 + im2 +infl)

return
end subroutine calc_newairtemp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine calculates the inverse of a NxN matrix
! input/output = a, size=s, error flag=ier
 
subroutine minverse(a,ierr)
 
implicit none
 
real(kind=8), dimension(:,:,:), intent(inout) :: a
real(kind=8), dimension(ufull_g) :: det, d, amax
real(kind=8), dimension(size(a,2)) :: x
real(kind=8) :: y
integer, intent(out)  :: ierr
integer s, ns, iq, i, j
integer, dimension(ufull_g,size(a,2)) :: row, col
integer, dimension(ufull_g) :: prow, pcol
logical, dimension(ufull_g,size(a,2)) :: notpiv

s = size(a,2)

det = 0.
d = 1.
notpiv = .TRUE.
 
do ns = 1,s
 
  amax(:) = 0.
  do j = 1,s
    do i = 1,s
      where ( notpiv(:,j) .and. notpiv(:,i) .and. amax(:)<abs(a(:,i,j)) )
        amax(:) = abs(a(:,i,j))
        prow(:) = i
        pcol(:) = j
      end where  
    end do
  end do
 
  if ( any(amax<0.) ) then
    ierr=1
    return
  end if

  do iq = 1,ufull_g
    notpiv(iq,pcol(iq)) = .FALSE.
    if ( prow(iq)/=pcol(iq) ) then
      d(iq) = -d(iq)
      x(1:s) = a(iq,prow(iq),1:s)
      a(iq,prow(iq),1:s) = a(iq,pcol(iq),1:s)
      a(iq,pcol(iq),1:s) = x(1:s)
    end if
  end do  
 
  row(:,ns) = prow(:)
  col(:,ns) = pcol(:)
  do iq = 1,ufull_g
    amax(iq) = a(iq,pcol(iq),pcol(iq))
  end do  
  d(:) = d(:)*amax(:)
 
  if ( any(abs(d)<=0.) ) then
    ierr=1
    return
  end if
 
  amax(:) = 1./amax(:)
  do iq = 1,ufull_g
    a(iq,pcol(iq),pcol(iq))=1.
    a(iq,pcol(iq),1:s) = a(iq,pcol(iq),1:s)*amax(iq)
  end do  
 
  do i=1,s
    do iq = 1,ufull_g  
      if ( i/=pcol(iq) ) then
        y = a(iq,i,pcol(iq))
        a(iq,i,pcol(iq)) = 0.
        do j = 1,s
          a(iq,i,j)=a(iq,i,j)-y*a(iq,pcol(iq),j)
        end do
      end if
    end do  
  end do
 
end do
 
det(:) = d(:)
 
do ns = s,1,-1
  prow(:) = row(:,ns)
  pcol(:) = col(:,ns)
  do iq = 1,ufull_g
    if ( prow(iq)/=pcol(iq) ) then
      do i = 1,s
        y = a(iq,i,prow(iq))
        a(iq,i,prow(iq)) = a(iq,i,pcol(iq))
        a(iq,i,pcol(iq)) = y
      end do
    end if
  end do  
end do
 
ierr = 0
 
return
end subroutine minverse

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module ateb
