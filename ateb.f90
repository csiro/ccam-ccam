
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
!   call atebsavem    ! to save current state arrays (for use by tebloadm)
!   call atebend      ! to deallocate memory before quiting

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
       atebdwn,atebscrnout,atebfbeam,atebspitter

! state arrays
integer, save :: ufull,ifull,iqut
logical, dimension(:), allocatable, save :: upack
real, dimension(:), allocatable, save :: sigmau
real, dimension(:,:), allocatable, save :: atebdwn ! These variables are for CCAM onthefly.f
real, dimension(:,:), allocatable, save :: rf_temp
real, dimension(:), allocatable, save :: rf_water,rf_snow,rf_den,rf_alpha
real, dimension(:,:), allocatable, save :: rd_temp
real, dimension(:), allocatable, save :: rd_water,rd_snow,rd_den,rd_alpha
real, dimension(:,:), allocatable, save :: we_temp,ww_temp
real, dimension(:), allocatable, save :: v_watrc,v_moistc,v_watrr,v_moistr
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

! model parameters
integer, parameter :: nmlfile=0       ! Read configuration from nml file (0=off, >0 unit number (default=11))
integer, save :: resmeth=1            ! Canyon sensible heat transfer (0=Masson, 1=Harman (varying width), 2=Kusaka, 3=Harman (fixed width))
integer, save :: useonewall=0         ! Combine both wall energy budgets into a single wall (0=two walls, 1=single wall) 
integer, save :: zohmeth=1            ! Urban roughness length for heat (0=0.1*zom, 1=Kanda, 2=0.003*zom)
integer, save :: acmeth=1             ! AC heat pump into canyon (0=Off, 1=On)
integer, save :: nrefl=3              ! Number of canyon reflections (default=3)
integer, save :: vegmode=2            ! In-canyon vegetation mode (0=50%/50%, 1=100%/0%, 2=0%/100%, where out/in=X/Y)
integer, save :: scrnmeth=1           ! Screen diagnostic method (0=Slab, 1=Hybrid, 2=Canyon)
integer, save :: wbrelaxc=0           ! Relax canyon soil moisture for irrigation (0=Off, 1=On)
integer, save :: wbrelaxr=0           ! Relax roof soil moisture for irrigation (0=Off, 1=On)
integer, save :: iqt=314              ! Diagnostic point (in terms of host grid)
! sectant solver parameters
integer, save :: nfgits=6             ! Maximum number of iterations for calculating sensible heat flux (default=4)
integer, save :: negits=3             ! Maximum number of iterations for calculating latent heat flux (default=2)
real, save    :: tol=0.001            ! Sectant method tolarance for sensible heat flux
real, save    :: tolqg=1.E-8          ! Sectant method tolarance for latent heat flux
real, save    :: alpha=0.7            ! Weighting for determining the rate of convergence when calculating canyon temperatures
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
real, parameter :: pi=3.1415927       ! pi
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
real, save :: zomratio=0.1            ! Ratio of roughness length to building height (default=10%)
real, save :: zocanyon=0.01           ! Roughness length of in-canyon surfaces (m)
real, save :: zoroof=0.01             ! Roughness length of roof surfaces (m)
real, save :: maxrfwater=1.           ! Maximum roof water (kg m^-2)
real, save :: maxrdwater=1.           ! Maximum road water (kg m^-2)
real, save :: maxrfsn=1.              ! Maximum roof snow (kg m^-2)
real, save :: maxrdsn=1.              ! Maximum road snow (kg m^-2)
real, save :: maxvwatf=0.1            ! Factor multiplied to LAI to predict maximum leaf water (kg m^-2)
! stability function parameters
integer, save :: icmax=5              ! number of iterations for stability functions (default=5)
real, save    :: a_1=1.
real, save    :: b_1=2./3.
real, save    :: c_1=5.
real, save    :: d_1=0.35


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine prepare the arrays used by the aTEB scheme
! This is a compulsory subroutine that must be called during
! model initalisation

subroutine atebinit(ifin,sigu,diag)

implicit none

integer, intent(in) :: ifin,diag
integer, dimension(ifin) :: utype
integer iqu,iq,ii
real, dimension(ifin), intent(in) :: sigu

if (diag.ge.1) write(6,*) "Initialising aTEB"

ifull=ifin
allocate(upack(ifull))
upack=sigu.gt.0.
ufull=count(upack)
if (ufull.eq.0) then
  deallocate(upack)
  return
end if

allocate(f_roofdepth(ufull,3),f_walldepth(ufull,3),f_roaddepth(ufull,3))
allocate(f_roofcp(ufull,3),f_wallcp(ufull,3),f_roadcp(ufull,3),f_vegcp(ufull,3))
allocate(f_rooflambda(ufull,3),f_walllambda(ufull,3),f_roadlambda(ufull,3),f_veglambda(ufull,3))
allocate(f_hwratio(ufull),f_sigmabld(ufull),f_industryfg(ufull),f_trafficfg(ufull),f_bldheight(ufull),f_vangle(ufull))
allocate(f_hangle(ufull),f_fbeam(ufull),f_roofalpha(ufull),f_wallalpha(ufull),f_roadalpha(ufull),f_ctime(ufull))
allocate(f_roofemiss(ufull),f_wallemiss(ufull),f_roademiss(ufull),f_bldtemp(ufull),f_sigmavegc(ufull),f_vegalphac(ufull))
allocate(f_vegemissc(ufull),f_sigmavegr(ufull),f_vegdepthr(ufull),f_vegalphar(ufull),f_vegemissr(ufull))
allocate(f_zovegc(ufull),f_vegrlaic(ufull),f_vegrsminc(ufull),f_zovegr(ufull),f_vegrlair(ufull),f_vegrsminr(ufull))
allocate(f_swilt(ufull),f_sfc(ufull),f_ssat(ufull))
allocate(p_lzom(ufull),p_lzoh(ufull),p_cndzmin(ufull),p_cduv(ufull),p_cdtq(ufull),p_vegtempc(ufull),p_vegtempr(ufull))
allocate(p_tscrn(ufull),p_qscrn(ufull),p_uscrn(ufull),p_u10(ufull),p_emiss(ufull))
allocate(rf_temp(ufull,3),rf_water(ufull),rf_snow(ufull))
allocate(rf_den(ufull),rf_alpha(ufull))
allocate(rd_temp(ufull,3),rd_water(ufull),rd_snow(ufull))
allocate(rd_den(ufull),rd_alpha(ufull))
allocate(we_temp(ufull,3),ww_temp(ufull,3))
allocate(sigmau(ufull),v_watrc(ufull),v_moistc(ufull),v_watrr(ufull),v_moistr(ufull))

! define grid arrays
sigmau=pack(sigu,upack)
iqu=0
iqut=0
do iq=1,ifull
  if (sigu(iq).gt.0.) then
    iqu=iqu+1
    if (iq.ge.iqt) then
      iqut=iqu
      exit
    end if
  end if
end do

if (iqut.eq.0) then
  !write(6,*) "WARN: Cannot located aTEB diagnostic point.  iqut=1"
  iqut=1
end if

! Initialise state variables
rf_temp=291.
rd_temp=291.
we_temp=291.
ww_temp=291.
rf_water=0.
rf_snow=0.
rf_den=minsnowden
rf_alpha=maxsnowalpha
rd_water=0.
rd_snow=0.
rd_den=minsnowden
rd_alpha=maxsnowalpha

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
f_sigmavegr=0.5
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

p_cndzmin=max(40.,0.1*f_bldheight+1.)   ! updated in atebcalc
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

v_moistc=0.5*(f_ssat+f_swilt)
v_watrc=0.
v_moistr=f_swilt
v_watrr=0.

return
end subroutine atebinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine deallocates arrays used by the TEB scheme

subroutine atebend(diag)

implicit none

integer, intent(in) :: diag

if (ufull.eq.0) return
if (diag.ge.1) write(6,*) "Deallocating aTEB arrays"

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
deallocate(rf_temp,rf_water,rf_snow)
deallocate(rf_den,rf_alpha)
deallocate(rd_temp,rd_water,rd_snow)
deallocate(rd_den,rd_alpha)
deallocate(we_temp,ww_temp)
deallocate(sigmau,v_watrc,v_moistc,v_watrr,v_moistr)

return
end subroutine atebend

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine loads aTEB state arrays (not compulsory)

subroutine atebload(urban,diag)

implicit none

integer, intent(in) :: diag
integer ii
real, dimension(ifull,24), intent(in) :: urban

if (ufull.eq.0) return
if (diag.ge.1) write(6,*) "Load aTEB state arrays"

do ii=1,3
  rf_temp(:,ii)=pack(urban(:,ii),upack)
  we_temp(:,ii)=pack(urban(:,ii+3),upack)
  ww_temp(:,ii)=pack(urban(:,ii+6),upack)
  rd_temp(:,ii)=pack(urban(:,ii+9),upack)
end do
v_moistc=pack(urban(:,13),upack)
v_moistr=pack(urban(:,14),upack)
rf_water=pack(urban(:,15),upack)
rd_water=pack(urban(:,16),upack)
v_watrc=pack(urban(:,17),upack)
v_watrr=pack(urban(:,18),upack)
rf_snow=pack(urban(:,19),upack)
rd_snow=pack(urban(:,20),upack)
rf_den=pack(urban(:,21),upack)
rd_den=pack(urban(:,22),upack)
rf_alpha=pack(urban(:,23),upack)
rd_alpha=pack(urban(:,24),upack)

return
end subroutine atebload

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! temperature only version of tebsave

subroutine atebloadm(urban,moist,diag)

implicit none

integer, intent(in) :: diag
integer ii
real, dimension(ifull,12), intent(in) :: urban
real, dimension(ifull,2), intent(in) :: moist

if (ufull.eq.0) return
if (diag.ge.1) write(6,*) "Load aTEB state arrays"

do ii=1,3
  rf_temp(:,ii)=pack(urban(:,ii),upack)
  we_temp(:,ii)=pack(urban(:,ii+3),upack)
  ww_temp(:,ii)=pack(urban(:,ii+6),upack)
  rd_temp(:,ii)=pack(urban(:,ii+9),upack)
end do
v_moistc=pack(moist(:,1),upack)
v_moistr=pack(moist(:,2),upack)

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
real, dimension(maxtype,3) :: croofdepth=reshape((/ 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, &
                                                    0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 0.40, &
                                                    0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10 /), &
                                                    (/maxtype,3/))
! Wall depths (m)
real, dimension(maxtype,3) :: cwalldepth=reshape((/ 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, &
                                                    0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, &
                                                    0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05 /), &
                                                    (/maxtype,3/))
! Road depths (m)
real, dimension(maxtype,3) :: croaddepth=reshape((/ 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, &
                                                    0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, &
                                                    4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00 /), &
                                                    (/maxtype,3/))
! Roof heat capacity (J m^-3 K^-1)
real, dimension(maxtype,3) :: croofcp=reshape((/ 2.11E6, 2.11E6, 2.11E6, 2.11E6, 2.11E6, 2.11E6, 2.11E6, 2.11E6, &
                                                 0.28E6, 0.28E6, 0.28E6, 0.28E6, 0.28E6, 0.28E6, 0.28E6, 0.28E6, &
                                                 0.29E6, 0.29E6, 0.29E6, 0.29E6, 0.29E6, 0.29E6, 0.29E6, 0.29E6 /), &
                                                 (/maxtype,3/))
! Wall heat capacity (J m^-3 K^-1)
real, dimension(maxtype,3) :: cwallcp=reshape((/ 1.55E6, 1.55E6, 1.55E6, 1.55E6, 1.55E6, 1.55E6, 1.55E6, 1.55E6, &
                                                 1.55E6, 1.55E6, 1.55E6, 1.55E6, 1.55E6, 1.55E6, 1.55E6, 1.55E6, &
                                                 0.29E6, 0.29E6, 0.29E6, 0.29E6, 0.29E6, 0.29E6, 0.29E6, 0.29E6 /), &
                                                 (/maxtype,3/))
! Road heat capacity (J m^-3 K^-1)
real, dimension(maxtype,3) :: croadcp=reshape((/ 1.94E6, 1.94E6, 1.94E6, 1.94E6, 1.94E6, 1.94E6, 1.94E6, 1.94E6, &
                                                 1.28E6, 1.28E6, 1.28E6, 1.28E6, 1.28E6, 1.28E6, 1.28E6, 1.28E6, &
                                                 1.28E6, 1.28E6, 1.28E6, 1.28E6, 1.28E6, 1.28E6, 1.28E6, 1.28E6 /), &
                                                 (/maxtype,3/))
! Roof conductance (W m^-1 K^-1)
real, dimension(maxtype,3) :: crooflambda=reshape((/ 1.5100, 1.5100, 1.5100, 1.5100, 1.5100, 1.5100, 1.5100, 1.5100, &
                                                     0.0800, 0.0800, 0.0800, 0.0800, 0.0800, 0.0800, 0.0800, 0.0800, &
                                                     0.0500, 0.0500, 0.0500, 0.0500, 0.0500, 0.0500, 0.0500, 0.0500 /), &
                                                     (/maxtype,3/))
! Wall conductance (W m^-1 K^-1)
real, dimension(maxtype,3) :: cwalllambda=reshape((/ 0.9338, 0.9338, 0.9338, 0.9338, 0.9338, 0.9338, 0.9338, 0.9338, &
                                                     0.9338, 0.9338, 0.9338, 0.9338, 0.9338, 0.9338, 0.9338, 0.9338, &
                                                     0.0500, 0.0500, 0.0500, 0.0500, 0.0500, 0.0500, 0.0500, 0.0500 /), &
                                                     (/maxtype,3/))
! Road conductance (W m^-1 K^-1)
real, dimension(maxtype,3) :: croadlambda=reshape((/ 0.7454, 0.7454, 0.7454, 0.7454, 0.7454, 0.7454, 0.7454, 0.7454, &
                                                     0.2513, 0.2513, 0.2513, 0.2513, 0.2513, 0.2513, 0.2513, 0.2513, &
                                                     0.2513, 0.2513, 0.2513, 0.2513, 0.2513, 0.2513, 0.2513, 0.2513 /), &
                                                     (/maxtype,3/))
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
                                                                
if (ufull.eq.0) return
if (diag.ge.1) write(6,*) "Load aTEB building properties"

itmp=pack(itype,upack)
if ((minval(itmp).lt.1).or.(maxval(itmp).gt.maxtype)) then
  write(6,*) "ERROR: Urban type is out of range"
  stop
end if

if (nmlfile.ne.0) then
  open(unit=nmlfile,file='ateb.nml',action="read",iostat=ierr)
  if (ierr.eq.0) then
    write(6,*) "Reading ateb.nml"
    read(nmlfile,nml=atebnml)
    read(nmlfile,nml=atebsnow)
    read(nmlfile,nml=atebgen)
    read(nmlfile,nml=atebtile)
    close(nmlfile)  
  end if
end if

select case(vegmode)
  case(0)
    tsigveg=1./(2./csigvegc(itmp)-1.)
    tsigmabld=csigmabld(itmp)*(1.-tsigveg)/(1.-csigvegc(itmp))
  case(1)
    tsigveg=0.
    tsigmabld=csigmabld(itmp)/(1.-csigvegc(itmp))
  case(2)
    tsigveg=csigvegc(itmp)
    tsigmabld=csigmabld(itmp)
  case DEFAULT
    write(6,*) "ERROR: Unsupported vegmode ",vegmode
    stop
end select
f_sigmavegc=tsigveg/(1.-tsigmabld)
f_sigmavegr=csigvegr(itmp)
f_sigmabld=tsigmabld
f_hwratio=chwratio(itmp)
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
do ii=1,3
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
f_vegrsminc=cvegrsminc(itmp)/f_vegrlaic
f_zovegr=czovegr(itmp)
f_vegrlair=cvegrlair(itmp)
f_vegrsminr=cvegrsminr(itmp)/f_vegrlair
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
real, dimension(ifull,54), intent(in) :: ifn

if (ufull.eq.0) return
if (diag.ge.1) write(6,*) "Load aTEB building properties"

f_hwratio=pack(ifn(:,1),upack)
f_sigmabld=pack(ifn(:,2),upack)
f_sigmavegc=pack(ifn(:,3),upack)/(1.-pack(ifn(:,2),upack))
f_sigmavegr=pack(ifn(:,4),upack)
f_industryfg=pack(ifn(:,5),upack)
f_trafficfg=pack(ifn(:,6),upack)
f_bldheight=pack(ifn(:,7),upack)
f_roofalpha=pack(ifn(:,8),upack)
f_wallalpha=pack(ifn(:,9),upack)
f_roadalpha=pack(ifn(:,10),upack)
f_vegalphac=pack(ifn(:,11),upack)
f_vegalphac=pack(ifn(:,12),upack)
f_roofemiss=pack(ifn(:,13),upack)
f_wallemiss=pack(ifn(:,14),upack)
f_roademiss=pack(ifn(:,15),upack)
f_vegemissc=pack(ifn(:,16),upack)
f_vegemissr=pack(ifn(:,17),upack)
f_bldtemp=pack(ifn(:,18),upack)
do ii=1,3
  f_roofdepth(:,ii)=pack(ifn(:,18+ii),upack)
  f_walldepth(:,ii)=pack(ifn(:,21+ii),upack)
  f_roaddepth(:,ii)=pack(ifn(:,24+ii),upack)
  f_roofcp(:,ii)=pack(ifn(:,27+ii),upack)
  f_wallcp(:,ii)=pack(ifn(:,30+ii),upack)
  f_roadcp(:,ii)=pack(ifn(:,33+ii),upack)
  f_rooflambda(:,ii)=pack(ifn(:,36+ii),upack)
  f_walllambda(:,ii)=pack(ifn(:,39+ii),upack)
  f_roadlambda(:,ii)=pack(ifn(:,42+ii),upack)
end do
f_zovegc=pack(ifn(:,46),upack)
f_vegrlaic=pack(ifn(:,47),upack)
f_vegrsminc=pack(ifn(:,48),upack)
f_zovegr=pack(ifn(:,49),upack)
f_vegrlair=pack(ifn(:,50),upack)
f_vegrsminr=pack(ifn(:,51),upack)
f_swilt=pack(ifn(:,52),upack)
f_sfc=pack(ifn(:,53),upack)
f_ssat=pack(ifn(:,54),upack)

return
end subroutine atebfndef

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine saves aTEB state arrays (not compulsory)

subroutine atebsave(urban,diag)

implicit none

integer, intent(in) :: diag
integer ii
real, dimension(ifull,24), intent(inout) :: urban

if (ufull.eq.0) return
if (diag.ge.1) write(6,*) "Save aTEB state arrays"

do ii=1,3
  urban(:,ii)=unpack(rf_temp(:,ii),upack,urban(:,ii))
  urban(:,ii+3)=unpack(we_temp(:,ii),upack,urban(:,ii+3))
  urban(:,ii+6)=unpack(ww_temp(:,ii),upack,urban(:,ii+6))
  urban(:,ii+9)=unpack(rd_temp(:,ii),upack,urban(:,ii+9))
end do
urban(:,13)=unpack(v_moistc,upack,urban(:,13))
urban(:,14)=unpack(v_moistr,upack,urban(:,14))
urban(:,15)=unpack(rf_water,upack,urban(:,15))
urban(:,16)=unpack(rd_water,upack,urban(:,16))
urban(:,17)=unpack(v_watrc,upack,urban(:,17))
urban(:,18)=unpack(v_watrr,upack,urban(:,18))
urban(:,19)=unpack(rf_snow,upack,urban(:,19))
urban(:,20)=unpack(rd_snow,upack,urban(:,20))
urban(:,21)=unpack(rf_den,upack,urban(:,21))
urban(:,22)=unpack(rd_den,upack,urban(:,22))
urban(:,23)=unpack(rf_alpha,upack,urban(:,23))
urban(:,24)=unpack(rd_alpha,upack,urban(:,24))

return
end subroutine atebsave

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! temperature only version of tebsave

subroutine atebsavem(urban,moist,diag)

implicit none

integer, intent(in) :: diag
integer ii
real, dimension(ifull,12), intent(inout) :: urban
real, dimension(ifull,2), intent(inout) :: moist

if (ufull.eq.0) return
if (diag.ge.1) write(6,*) "Save aTEB state arrays"

do ii=1,3
  urban(:,ii)=unpack(rf_temp(:,ii),upack,urban(:,ii))
  urban(:,ii+3)=unpack(we_temp(:,ii),upack,urban(:,ii+3))
  urban(:,ii+6)=unpack(ww_temp(:,ii),upack,urban(:,ii+6))
  urban(:,ii+9)=unpack(rd_temp(:,ii),upack,urban(:,ii+9))
end do
moist(:,1)=unpack(v_moistc,upack,moist(:,1))
moist(:,2)=unpack(v_moistr,upack,moist(:,2))

return
end subroutine atebsavem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine blends urban momentum and heat roughness lengths
! (This version neglects the displacement height (e.g., for CCAM))
!

subroutine atebzo(zom,zoh,diag,raw)

implicit none

integer, intent(in) :: diag
real, dimension(ifull), intent(inout) :: zom,zoh
real, dimension(ufull) :: workb,workc,zmtmp,zhtmp
real, parameter :: zr=1.e-15 ! limits minimum roughness length for heat
logical, intent(in), optional :: raw
logical mode

if (ufull.eq.0) return
if (diag.ge.1) write(6,*) "Calculate urban roughness lengths"

mode=.false.
if (present(raw)) mode=raw

if (mode) then
  zom=unpack(p_cndzmin*exp(-p_lzom),upack,zom)
  zoh=unpack(p_cndzmin*exp(-p_lzoh),upack,zoh)
else 
  ! evaluate at canyon displacement height (really the atmospheric model should provide a displacement height)
  zmtmp=pack(zom,upack)
  zhtmp=pack(zoh,upack)
  workb=sqrt((1.-sigmau)/log(p_cndzmin/zmtmp)**2+sigmau/p_lzom**2)
  workc=(1.-sigmau)/(log(p_cndzmin/zmtmp)*log(p_cndzmin/zhtmp))+sigmau/(p_lzom*p_lzoh)
  workc=workc/workb
  workb=p_cndzmin*exp(-1./workb)
  workc=max(p_cndzmin*exp(-1./workc),zr)
  zom=unpack(workb,upack,zom)
  zoh=unpack(workc,upack,zoh)
  if (minval(workc).le.zr) write(6,*) "WARN: minimum zoh reached"
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

if (ufull.eq.0) return

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

if (ufull.eq.0) return

ifinish=is+ifin-1
ucount=count(upack(is:ifinish))
if (ucount.eq.0) return

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

if (ufull.eq.0) return

ifinish=is+ifin-1
ucount=count(upack(is:ifinish))
if (ucount.eq.0) return

ib=count(upack(1:is-1))+1
ie=ucount+ib-1

f_fbeam(ib:ie)=0.
lsg(ib:ie)=pack(sg,upack(is:ifinish))
lcosin(ib:ie)=pack(cosin,upack(is:ifinish))

tmpr(ib:ie)=0.847+lcosin(ib:ie)*(1.04*lcosin(ib:ie)-1.61)
tmpk(ib:ie)=(1.47-tmpr(ib:ie))/1.66
tmprat(ib:ie)=0.
where (lcosin(ib:ie).gt.1.0e-10 .and. lsg(ib:ie).gt.10.)
  tmprat(ib:ie)=lsg(ib:ie)/(solcon*(1.+0.033*cos(2.*pi*(fjd-10.)/365.))*lcosin(ib:ie))
end where
where (tmprat(ib:ie).gt.0.22)
  f_fbeam(ib:ie)=6.4*(tmprat(ib:ie)-0.22)**2
elsewhere (tmprat(ib:ie).gt.0.35)
  f_fbeam(ib:ie)=min(1.66*tmprat(ib:ie)-0.4728,1.)
elsewhere (tmprat(ib:ie).gt.tmpk(ib:ie))
  f_fbeam(ib:ie)=max(1.-tmpr(ib:ie),0.)
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
integer i,ucount,ib,ie,ifinish,albmode
integer, intent(in), optional :: split
real, dimension(ifin), intent(inout) :: alb
real, dimension(ufull) :: ualb,utmp
logical, intent(in), optional :: raw
logical outmode

if (ufull.eq.0) return

outmode=.false.
if (present(raw)) outmode=raw

albmode=0
if (present(split)) albmode=split

ifinish=is+ifin-1
ucount=count(upack(is:ifinish))
if (ucount.eq.0) return

ib=count(upack(1:is-1))+1
ie=ucount+ib-1
call atebalbcalc(ib,ucount,ualb(ib:ie),albmode,diag)

if (outmode) then
  alb=unpack(ualb(ib:ie),upack(is:ifinish),alb)
else
  utmp(ib:ie)=pack(alb,upack(is:ifinish))
  utmp(ib:ie)=(1.-sigmau(ib:ie))*utmp(ib:ie)+sigmau(ib:ie)*ualb(ib:ie)
  alb=unpack(utmp(ib:ie),upack(is:ifinish),alb)
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

! canyon
effbldheight=max(f_bldheight(is:ie)-6.*f_zovegc(is:ie),0.1)/f_bldheight(is:ie)  ! MJT suggestion for tall vegetation
effhwratio=f_hwratio(is:ie)*effbldheight                                        ! MJT suggestion for tall vegetation
snowdelta=rd_snow(is:ie)/(rd_snow(is:ie)+maxrdsn)
call getswcoeff(ifin,sg_roof,sg_vegr,sg_road,sg_walle,sg_wallw,sg_vegc,sg_rfsn,sg_rdsn,wallpsi,roadpsi,effhwratio, &
                f_vangle(is:ie),f_hangle(is:ie),dumfbeam,f_sigmavegc(is:ie),f_roadalpha(is:ie),f_vegalphac(is:ie), &
                f_wallalpha(is:ie),rd_alpha(is:ie),snowdelta)
albu=1.-(f_hwratio(is:ie)*(sg_walle+sg_wallw)*(1.-f_wallalpha(is:ie))+snowdelta*sg_rdsn*(1.-rd_alpha(is:ie)) &
    +(1.-snowdelta)*((1.-f_sigmavegc(is:ie))*sg_road*(1.-f_roadalpha(is:ie))+f_sigmavegc(is:ie)*sg_vegc*(1.-f_vegalphac(is:ie))))

! roof
snowdelta=rf_snow(is:ie)/(rf_snow(is:ie)+maxrfsn)
albr=(1.-snowdelta)*((1.-f_sigmavegr(is:ie))*sg_roof*f_roofalpha(is:ie)+f_sigmavegr(is:ie)*sg_vegr*f_vegalphar(is:ie)) &
    +snowdelta*sg_rfsn*rf_alpha(is:ie)

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

if (ufull.eq.0) return

ifinish=is+ifin-1
ucount=count(upack(is:ifinish))
if (ucount.eq.0) return

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

integer, intent(in) :: is,ifin
integer ifinish,ucount,ib,ie
real, intent(in) :: fjd,slag,dt,sdlt
real cdlt
real, dimension(ifin), intent(in) :: cosin,rlon,rlat
real, dimension(ufull) :: hloc,x,y,lattmp

! cosin = cosine of zenith angle
! rlon = longitude
! rlat = latitude
! fjd = day of year
! slag = sun lag angle
! sdlt = sin declination of sun

if (ufull.eq.0) return

ifinish=is+ifin-1
ucount=count(upack(is:ifinish))
if (ucount.eq.0) return

ib=count(upack(1:is-1))+1
ie=ucount+ib-1

cdlt=sqrt(min(max(1.-sdlt*sdlt,0.),1.))

lattmp(ib:ie)=pack(rlat,upack(is:ifinish))

! from CCAM zenith.f
hloc(ib:ie)=2.*pi*fjd+slag+pi+pack(rlon,upack(is:ifinish))+dt*pi/86400.
! estimate azimuth angle
x(ib:ie)=sin(-hloc(ib:ie))*cdlt
y(ib:ie)=-cos(-hloc(ib:ie))*cdlt*sin(lattmp(ib:ie))+cos(lattmp(ib:ie))*sdlt
!azimuth=atan2(x,y)
f_hangle(ib:ie)=0.5*pi-atan2(x(ib:ie),y(ib:ie))
f_vangle(ib:ie)=acos(pack(cosin,upack(is:ifinish)))
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
real, dimension(ufull) :: tmp
logical, intent(in), optional :: raw
logical mode

if (ufull.eq.0) return

mode=.false.
if (present(raw)) mode=raw

if (mode) then
  tscrn=unpack(p_tscrn,upack,tscrn)
  qscrn=unpack(p_qscrn,upack,qscrn)
  uscrn=unpack(p_uscrn,upack,uscrn)
  u10=unpack(p_u10,upack,u10)
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

if (ufull.eq.0) return

mode=.false.
if (present(raw)) mode=raw

a_zmin=pack(zmin,upack)
a_sg=pack(sg,upack)
a_rg=pack(rg,upack)
a_rho=pack(rho,upack)
a_temp=pack(temp,upack)
a_mixr=pack(mixr,upack)
a_ps=pack(ps,upack)
a_umag=max(sqrt(pack(uu,upack)**2+pack(vv,upack)**2),umin)
a_udir=atan2(pack(vv,upack),pack(uu,upack))
a_rnd=pack(rnd-snd,upack)
a_snd=pack(snd,upack)

call atebeval(u_fg,u_eg,u_ts,u_wf,u_rn,dt,a_sg,a_rg,a_rho,a_temp,a_mixr,a_ps,a_umag,a_udir,a_rnd,a_snd,a_zmin,diag)

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
integer k,ii,cns
logical, dimension(ufull) :: pgs
real, intent(in) :: ddt
real, dimension(ufull,2:3) :: ggaroof,ggawall,ggaroad
real, dimension(ufull,1:3) :: ggbroof,ggbwall,ggbroad
real, dimension(ufull,1:2) :: ggcroof,ggcwall,ggcroad
real, dimension(ufull,3) :: garoof,gawalle,gawallw,garoad
real, dimension(ufull) :: garfsn,gardsn,gflxroof,gflxwalle,gflxwallw
real, dimension(ufull) :: rdsntemp,rfsntemp,rdsnmelt,rfsnmelt
real, dimension(ufull) :: wallpsi,roadpsi,fgtop,egtop,qsatr
real, dimension(ufull) :: oldval,newval,cu,ctmax,ctmin,evctx,evct
real, dimension(ufull) :: ln,rn,we,ww,wr,zolog,a,xe,xw,cuven,n,zom,zonet,dis
real, dimension(ufull) :: p_wallpsi,p_roadpsi
real, dimension(ufull) :: p_sntemp,p_gasn,p_snmelt
real, dimension(ufull) :: z_on_l,pa,dts,dtt,effbldheight,effhwratio
real, dimension(ufull), intent(in) :: a_sg,a_rg,a_rho,a_temp,a_mixr,a_ps,a_umag,a_udir,a_rnd,a_snd,a_zmin
real, dimension(ufull) :: p_a_umag,p_a_rho,p_a_rg,p_a_rnd,p_a_snd
real, dimension(ufull), intent(out) :: u_fg,u_eg,u_ts,u_wf,u_rn
real, dimension(ufull) :: sg_roof,sg_vegr,sg_road,sg_walle,sg_wallw,sg_vegc,sg_rfsn,sg_rdsn
real, dimension(ufull) :: rg_roof,rg_road,rg_walle,rg_wallw,rg_vegc,rg_vegr,rg_rfsn,rg_rdsn
real, dimension(ufull) :: fg_roof,fg_road,fg_walle,fg_wallw,fg_vegc,fg_vegr,fg_rfsn,fg_rdsn
real, dimension(ufull) :: eg_roof,eg_road,eg_vegc,eg_vegr,eg_rfsn,eg_rdsn
real, dimension(ufull) :: acond_roof,acond_road,acond_walle,acond_wallw,acond_vegc,acond_vegr,acond_rfsn,acond_rdsn
real, dimension(ufull) :: p_sg_vegc,p_sg_rs
real, dimension(ufull) :: p_rg_ro,p_rg_walle,p_rg_wallw,p_rg_vegc,p_rg_rs
real, dimension(ufull) :: p_fg_ro,p_fg_walle,p_fg_wallw,p_fg_vegc,p_fg_rs
real, dimension(ufull) :: p_eg_ro,p_eg_vegc,p_eg_rs
real, dimension(ufull) :: p_acond_ro,p_acond_walle,p_acond_wallw,p_acond_vegc,p_acond_rs
real, dimension(ufull) :: d_roofdelta,d_roaddelta,d_vegdeltac,d_vegdeltar,d_rfsndelta,d_rdsndelta
real, dimension(ufull) :: d_tempc,d_mixrc,d_tempr,d_mixrr,d_sigd,d_sigr,d_rfdzmin
real, dimension(ufull) :: d_accool,d_canyonrgout,d_roofrgout,d_tranc,d_evapc,d_tranr,d_evapr,d_c1c,d_c1r
real, dimension(ufull) :: d_totdepth,d_netemiss,d_netrad,d_topu
real, dimension(ufull) :: d_cwa,d_cwe,d_cww,d_cwr,d_cra,d_crr,d_crw
real, dimension(ufull) :: d_canyontemp,d_canyonmix,d_acout,d_traf
real, dimension(ufull) :: p_d_dzmin,p_d_sig,p_d_temp,p_d_mixr,p_d_rsdelta,p_d_canyontemp,p_d_canyonmix,p_d_topu,p_d_roaddelta
real, dimension(ufull) :: p_d_vegdeltac,p_d_netrad,p_d_netemiss,p_d_tranc,p_d_evapc,p_d_accool,p_d_totdepth,p_d_c1c
real, dimension(ufull) :: p_d_cwa,p_d_cra,p_d_cwe,p_d_cww,p_d_crw,p_d_crr,p_d_cwr,p_d_acout,p_d_traf,p_d_canyonrgout
real, dimension(ufull) :: p_ro_temp,p_ro_water,p_ro_snow,p_ro_den
real, dimension(ufull) :: p_v_watrc,p_v_moistc
real, dimension(ufull) :: p_we_temp,p_ww_temp
real, dimension(ufull) :: p_f_bldheight,p_f_hwratio,p_f_sigmavegc,p_f_roademiss,p_f_vegemissc,p_f_wallemiss
real, dimension(ufull) :: p_f_ctime,p_f_trafficfg,p_f_rodepth,p_f_rolambda,p_f_sigmabld,p_f_bldtemp
real, dimension(ufull) :: p_f_zovegc,p_f_vegrlaic,p_f_vegrsminc,p_f_swilt,p_f_sfc
real, dimension(ufull) :: p_p_vegtempc,p_p_lzom,p_p_lzoh,p_p_cndzmin

if (diag.ge.1) write(6,*) "Evaluating aTEB"

! limit prognostic state variables
do ii=1,3
  rf_temp(:,ii)=min(max(rf_temp(:,ii),200.),400.)
  we_temp(:,ii)=min(max(we_temp(:,ii),200.),400.)
  ww_temp(:,ii)=min(max(ww_temp(:,ii),200.),400.)
  rd_temp(:,ii)=min(max(rd_temp(:,ii),200.),400.)
end do
v_moistc=min(max(v_moistc,f_swilt),f_ssat)
v_moistr=min(max(v_moistr,f_swilt),f_ssat)
rf_water=min(max(rf_water,0.),maxrfwater)
rd_water=min(max(rd_water,0.),maxrdwater)
v_watrc=min(max(v_watrc,0.),maxvwatf*f_vegrlaic)
v_watrr=min(max(v_watrr,0.),maxvwatf*f_vegrlair)
rf_snow=min(max(rf_snow,0.),maxrfsn)
rd_snow=min(max(rd_snow,0.),maxrdsn)
rf_den=min(max(rf_den,minsnowden),maxsnowden)
rd_den=min(max(rd_den,minsnowden),maxsnowden)
rf_alpha=min(max(rf_alpha,minsnowalpha),maxsnowalpha)
rd_alpha=min(max(rd_alpha,minsnowalpha),maxsnowalpha)

! new snowfall
where (a_snd.gt.0.)
  rf_den=(rf_snow*rf_den+a_snd*ddt*minsnowden)/(rf_snow+ddt*a_snd)
  rd_den=(rd_snow*rd_den+a_snd*ddt*minsnowden)/(rd_snow+ddt*a_snd)
  rf_alpha=maxsnowalpha
  rd_alpha=maxsnowalpha
end where
    
! water and snow cover fractions
d_roofdelta=(rf_water/maxrfwater)**(2./3.)
d_roaddelta=(rd_water/maxrdwater)**(2./3.)
d_vegdeltac=(v_watrc/(maxvwatf*f_vegrlaic))**(2./3.)
d_vegdeltar=(v_watrr/(maxvwatf*f_vegrlair))**(2./3.)
d_rfsndelta=rf_snow/(rf_snow+maxrfsn)
d_rdsndelta=rd_snow/(rd_snow+maxrdsn)

! canyon (displacement height at refheight*building height)
d_sigd=a_ps*exp(-grav*f_bldheight*refheight/(rd*a_temp))
pa=a_ps*exp(-grav*a_zmin/(rd*a_temp))
d_tempc=a_temp*(d_sigd/pa)**(rd/aircp)
call getqsat(ufull,qsatr,d_tempc,d_sigd)
call getqsat(ufull,a,a_temp,pa)
d_mixrc=a_mixr*qsatr/a ! a=qsata

! roof (displacement height at building height)
d_sigr=a_ps*exp(-grav*f_bldheight/(rd*a_temp))
d_tempr=a_temp*(d_sigr/pa)**(rd/aircp)
call getqsat(ufull,qsatr,d_tempr,d_sigr)
d_mixrr=a_mixr*qsatr/a ! a=qsata

! total soil depth (for in-canyon vegetation)
d_totdepth=sum(f_roaddepth,2)

! calculate c1 for soil (for vegetation)
call getc1(d_c1c,v_moistc)
call getc1(d_c1r,v_moistr)

! calculate heat pumped into canyon by air conditioning (COP updated in solvecanyon)
! (use split form to estimate G_{*,3} flux into room for AC.  n is an estimate of the temperature at tau+1)
n=rf_temp(:,3)-2.*f_rooflambda(:,3)*(rf_temp(:,3)-f_bldtemp)/(f_roofcp(:,3)*f_roofdepth(:,3)*f_roofdepth(:,3)/ddt+2.*f_rooflambda(:,3))
gflxroof=(1.-f_sigmavegr)*2.*f_rooflambda(:,3)*(n-f_bldtemp)/f_roofdepth(:,3)
n=we_temp(:,3)-2.*f_walllambda(:,3)*(we_temp(:,3)-f_bldtemp)/(f_wallcp(:,3)*f_walldepth(:,3)*f_walldepth(:,3)/ddt+2.*f_walllambda(:,3))
gflxwalle=2.*f_walllambda(:,3)*(n-f_bldtemp)/f_walldepth(:,3)
n=ww_temp(:,3)-2.*f_walllambda(:,3)*(ww_temp(:,3)-f_bldtemp)/(f_wallcp(:,3)*f_walldepth(:,3)*f_walldepth(:,3)/ddt+2.*f_walllambda(:,3))
gflxwallw=2.*f_walllambda(:,3)*(n-f_bldtemp)/f_walldepth(:,3)
if (acmeth.eq.1) then
  d_acout=max(0.,gflxroof*f_sigmabld/(1.-f_sigmabld)+f_hwratio*(gflxwalle+gflxwallw))
  !d_acout=gflxroof*f_sigmabld/(1.-f_sigmabld)+f_hwratio*(gflxwalle+gflxwallw) ! test energy conservation
else
  d_acout=0.
end if

! calculate shortwave reflections
effbldheight=max(f_bldheight-6.*f_zovegc,0.1)/f_bldheight  ! MJT suggestion for tall vegetation
effhwratio=f_hwratio*effbldheight                          ! MJT suggestion for tall vegetation
call getswcoeff(ufull,sg_roof,sg_vegr,sg_road,sg_walle,sg_wallw,sg_vegc,sg_rfsn,sg_rdsn,wallpsi,roadpsi,effhwratio, &
                f_vangle,f_hangle,f_fbeam,f_sigmavegc,f_roadalpha,f_vegalphac,f_wallalpha,rd_alpha,d_rdsndelta)
sg_roof =(1.-f_roofalpha)*sg_roof*a_sg
sg_vegr =(1.-f_vegalphar)*sg_vegr*a_sg
sg_walle=(1.-f_wallalpha)*sg_walle*a_sg*effbldheight ! vegc shadow
sg_wallw=(1.-f_wallalpha)*sg_wallw*a_sg*effbldheight ! vegc shadow
sg_road =(1.-f_roadalpha)*sg_road*a_sg
sg_vegc =(1.-f_vegalphac)*sg_vegc*a_sg
sg_rfsn =(1.-rf_alpha)*sg_rfsn*a_sg
sg_rdsn =(1.-rd_alpha)*sg_rdsn*a_sg

! calculate long wave reflections to nrefl order (pregenerated before solvecanyon subroutine)
call getlwcoeff(ufull,d_netemiss,d_cwa,d_cra,d_cwe,d_cww,d_crw,d_crr,d_cwr,d_rdsndelta,wallpsi,roadpsi,f_sigmavegc,f_roademiss, &
                f_vegemissc,f_wallemiss)
n=d_rfsndelta*snowemiss+(1.-d_rfsndelta)*((1.-f_sigmavegr)*f_roofemiss+f_sigmavegr*f_vegemissr)
p_emiss=f_sigmabld*n+(1.-f_sigmabld)*(2.*f_wallemiss*f_hwratio*d_cwa+d_netemiss*d_cra) ! diagnostic only

! estimate in-canyon surface roughness length
dis=max(max(max(0.1*f_bldheight,zocanyon+0.1),f_zovegc+0.1),zosnow+0.1)
zolog=1./sqrt(d_rdsndelta/log(dis/zosnow)**2             &
     +(1.-d_rdsndelta)*(f_sigmavegc/log(dis/f_zovegc)**2 &
     +(1.-f_sigmavegc)/log(dis/zocanyon)**2))
zonet=dis*exp(-zolog)

! estimate overall urban roughness length
zom=zomratio*f_bldheight
where (zom*f_sigmabld.lt.zonet*(1.-f_sigmabld)) ! MJT suggestion
  zom=zonet
end where
n=rd_snow/(rd_snow+maxrdsn+0.408*grav*zom)       ! snow cover for urban roughness calc (Douville, et al 1995)
zom=(1.-n)*zom+n*zosnow                          ! blend urban and snow roughness lengths (i.e., snow fills canyon)
d_rfdzmin=max(max(abs(a_zmin-f_bldheight),zoroof+0.1),f_zovegr+0.1) ! distance to roof displacement height
where (a_zmin.ge.f_bldheight)
  p_lzom=log(a_zmin/zom)
elsewhere ! lowest atmospheric model level is within the canopy.  Need to interact with boundary layer scheme.
  p_lzom=log(f_bldheight/zom)*exp(0.5*f_hwratio*(1.-a_zmin/f_bldheight))
end where
p_cndzmin=zom*exp(p_lzom)                        ! distance to canyon displacement height

! calculate canyon wind speed and bulk transfer coefficents
! (i.e., acond = 1/(aerodynamic resistance) )
select case(resmeth)
  case(0) ! Masson (2000)
    cu=exp(-0.25*f_hwratio)
    acond_road=cu ! bulk transfer coefficents are updated in solvecanyon
    acond_walle=cu
    acond_wallw=cu
    acond_rdsn=cu
    acond_vegc=cu
  case(1,3) ! Harman et al (2004)
    if (resmeth.eq.1) then
      call getincanwind(we,ww,wr,a_udir,zonet)
    else
      call getincanwindb(we,ww,wr,a_udir,zonet)
    end if
    dis=max(max(max(0.1*f_bldheight,zocanyon+0.1),f_zovegc+0.1),zosnow+0.1)
    zolog=log(dis/zocanyon)
    a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
    acond_road=a*wr                  ! road bulk transfer
    acond_walle=a*we                 ! east wall bulk transfer
    acond_wallw=a*ww                 ! west wall bulk transfer
    zolog=log(dis/f_zovegc)
    a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
    acond_vegc=a*wr
    zolog=log(dis/zosnow)
    a=vkar*vkar/(zolog*(2.3+zolog))  ! Assume zot=zom/10.
    acond_rdsn=a*wr                  ! road snow bulk transfer
  case(2) ! Kusaka et al (2001)
    cu=exp(-0.25*f_hwratio)
    acond_road=cu ! bulk transfer coefficents are updated in solvecanyon
    acond_walle=cu
    acond_wallw=cu
    acond_rdsn=cu
    acond_vegc=cu
end select

! join two walls into a single wall (testing only)
if (useonewall.eq.1) then
  do k=1,3
    we_temp(:,k)=0.5*(we_temp(:,k)+ww_temp(:,k))
    ww_temp(:,k)=we_temp(:,k)
  end do
  acond_walle=0.5*(acond_walle+acond_wallw)
  acond_wallw=acond_walle
  sg_walle=0.5*(sg_walle+sg_wallw)
  sg_wallw=sg_walle
end if

! solve for road snow temperature -------------------------------
! includes solution to vegetation canopy temperature, canopy temperature, sensible heat flux and longwave radiation
pgs=d_rdsndelta.gt.0.
cns=count(pgs)
if (cns.gt.0) then ! road snow
  p_ro_temp(1:cns)=pack(rd_temp(:,1),pgs)
  p_ro_water(1:cns)=pack(rd_water,pgs)
  p_ro_snow(1:cns)=pack(rd_snow,pgs)
  p_ro_den(1:cns)=pack(rd_den,pgs)
  p_we_temp(1:cns)=pack(we_temp(:,1),pgs)
  p_ww_temp(1:cns)=pack(ww_temp(:,1),pgs)
  p_v_watrc(1:cns)=pack(v_watrc,pgs)
  p_v_moistc(1:cns)=pack(v_moistc,pgs)
  p_d_sig(1:cns)=pack(d_sigd,pgs)
  p_d_temp(1:cns)=pack(d_tempc,pgs)
  p_d_mixr(1:cns)=pack(d_mixrc,pgs)
  p_d_roaddelta(1:cns)=pack(d_roaddelta,pgs)
  p_d_vegdeltac(1:cns)=pack(d_vegdeltac,pgs)
  p_d_rsdelta(1:cns)=pack(d_rdsndelta,pgs)
  p_d_netemiss(1:cns)=pack(d_netemiss,pgs)
  p_d_cwa(1:cns)=pack(d_cwa,pgs)
  p_d_cra(1:cns)=pack(d_cra,pgs)
  p_d_cwe(1:cns)=pack(d_cwe,pgs)
  p_d_cww(1:cns)=pack(d_cww,pgs)
  p_d_crw(1:cns)=pack(d_crw,pgs)
  p_d_crr(1:cns)=pack(d_crr,pgs)
  p_d_cwr(1:cns)=pack(d_cwr,pgs)
  p_d_totdepth(1:cns)=pack(d_totdepth,pgs)
  p_d_c1c(1:cns)=pack(d_c1c,pgs)
  p_d_acout(1:cns)=pack(d_acout,pgs)
  p_sg_vegc(1:cns)=pack(sg_vegc,pgs)
  p_sg_rs(1:cns)=pack(sg_rdsn,pgs)
  p_rg_ro(1:cns)=pack(rg_road,pgs)
  p_rg_walle(1:cns)=pack(rg_walle,pgs)
  p_rg_wallw(1:cns)=pack(rg_wallw,pgs)
  p_rg_vegc(1:cns)=pack(rg_vegc,pgs)
  p_rg_rs(1:cns)=pack(rg_rdsn,pgs)
  p_fg_ro(1:cns)=pack(fg_road,pgs)
  p_fg_walle(1:cns)=pack(fg_walle,pgs)
  p_fg_wallw(1:cns)=pack(fg_wallw,pgs)
  p_fg_vegc(1:cns)=pack(fg_vegc,pgs)
  p_fg_rs(1:cns)=pack(fg_rdsn,pgs)
  p_eg_ro(1:cns)=pack(eg_road,pgs)
  p_eg_vegc(1:cns)=pack(eg_vegc,pgs)
  p_eg_rs(1:cns)=pack(eg_rdsn,pgs)
  p_a_umag(1:cns)=pack(a_umag,pgs)
  p_a_rho(1:cns)=pack(a_rho,pgs)
  p_a_rg(1:cns)=pack(a_rg,pgs)
  p_a_rnd(1:cns)=pack(a_rnd,pgs)
  p_a_snd(1:cns)=pack(a_snd,pgs)
  p_acond_ro(1:cns)=pack(acond_road,pgs)
  p_acond_walle(1:cns)=pack(acond_walle,pgs)
  p_acond_wallw(1:cns)=pack(acond_wallw,pgs)
  p_acond_vegc(1:cns)=pack(acond_vegc,pgs)
  p_acond_rs(1:cns)=pack(acond_rdsn,pgs)
  p_wallpsi(1:cns)=pack(wallpsi,pgs)
  p_roadpsi(1:cns)=pack(roadpsi,pgs)
  p_f_bldheight(1:cns)=pack(f_bldheight,pgs)
  p_f_hwratio(1:cns)=pack(f_hwratio,pgs)
  p_f_sigmavegc(1:cns)=pack(f_sigmavegc,pgs)
  p_f_roademiss(1:cns)=pack(f_roademiss,pgs)
  p_f_vegemissc(1:cns)=pack(f_vegemissc,pgs)
  p_f_wallemiss(1:cns)=pack(f_wallemiss,pgs)
  p_f_ctime(1:cns)=pack(f_ctime,pgs)
  p_f_trafficfg(1:cns)=pack(f_trafficfg,pgs)
  p_f_rodepth(1:cns)=pack(f_roaddepth(:,1),pgs)
  p_f_rolambda(1:cns)=pack(f_roadlambda(:,1),pgs)
  p_f_sigmabld(1:cns)=pack(f_sigmabld,pgs)
  p_f_bldtemp(1:cns)=pack(f_bldtemp,pgs)
  p_f_zovegc(1:cns)=pack(f_zovegc,pgs)
  p_f_vegrlaic(1:cns)=pack(f_vegrlaic,pgs)
  p_f_vegrsminc(1:cns)=pack(f_vegrsminc,pgs)
  p_f_swilt(1:cns)=pack(f_swilt,pgs)
  p_f_sfc(1:cns)=pack(f_sfc,pgs)
  p_p_vegtempc(1:cns)=pack(p_vegtempc,pgs)
  p_p_lzom(1:cns)=pack(p_lzom,pgs)
  p_p_lzoh(1:cns)=pack(p_lzoh,pgs)
  p_p_cndzmin(1:cns)=pack(p_cndzmin,pgs)
  ctmax(1:cns)=max(p_d_temp(1:cns),p_ro_temp(1:cns),p_we_temp(1:cns),p_ww_temp(1:cns),p_p_vegtempc(1:cns))+10. ! max road snow temp
  ctmin(1:cns)=min(p_d_temp(1:cns),p_ro_temp(1:cns),p_we_temp(1:cns),p_ww_temp(1:cns),p_p_vegtempc(1:cns))-10. ! min road snow temp
  p_sntemp(1:cns)=0.75*ctmax(1:cns)+0.25*ctmin(1:cns)
  call solverdsn(cns,evct(1:cns),p_rg_ro(1:cns),p_rg_walle(1:cns),p_rg_wallw(1:cns),p_rg_vegc(1:cns),p_rg_rs(1:cns),     &
                 p_fg_ro(1:cns),p_fg_walle(1:cns),p_fg_wallw(1:cns),p_fg_vegc(1:cns),p_fg_rs(1:cns),                     &
                 p_eg_ro(1:cns),p_eg_vegc(1:cns),p_eg_rs(1:cns),p_gasn(1:cns),p_snmelt(1:cns),p_sntemp(1:cns),           &
                 p_ro_temp(1:cns),p_ro_water(1:cns),p_ro_snow(1:cns),p_ro_den(1:cns),p_we_temp(1:cns),                   &
                 p_ww_temp(1:cns),p_v_watrc(1:cns),p_v_moistc(1:cns),p_d_canyontemp(1:cns),p_d_canyonmix(1:cns),         &
                 p_d_sig(1:cns),p_d_temp(1:cns),p_d_mixr(1:cns),p_d_topu(1:cns),p_d_roaddelta(1:cns),                    &
                 p_d_vegdeltac(1:cns),p_d_rsdelta(1:cns),p_d_netrad(1:cns),p_d_netemiss(1:cns),p_d_tranc(1:cns),         &
                 p_d_evapc(1:cns),p_d_cwa(1:cns),p_d_cra(1:cns),p_d_cwe(1:cns),p_d_cww(1:cns),p_d_crw(1:cns),            &
                 p_d_crr(1:cns),p_d_cwr(1:cns),p_d_accool(1:cns),p_d_totdepth(1:cns),                                    &
                 p_d_c1c(1:cns),p_d_acout(1:cns),p_d_traf(1:cns),p_d_canyonrgout(1:cns),p_sg_vegc(1:cns),p_sg_rs(1:cns), &
                 p_a_umag(1:cns),p_a_rho(1:cns),p_a_rg(1:cns),p_a_rnd(1:cns),p_a_snd(1:cns),ddt,p_acond_ro(1:cns),       &
                 p_acond_walle(1:cns),p_acond_wallw(1:cns),p_acond_vegc(1:cns),p_acond_rs(1:cns),p_wallpsi(1:cns),       &
                 p_roadpsi(1:cns),p_f_bldheight(1:cns),p_f_hwratio(1:cns),p_f_sigmavegc(1:cns),p_f_roademiss(1:cns),     &
                 p_f_vegemissc(1:cns),p_f_wallemiss(1:cns),p_f_ctime(1:cns),p_f_trafficfg(1:cns),p_f_rodepth(1:cns),     &
                 p_f_rolambda(1:cns),p_f_sigmabld(1:cns),p_f_bldtemp(1:cns),p_f_zovegc(1:cns),p_f_vegrlaic(1:cns),       &
                 p_f_vegrsminc(1:cns),p_f_swilt(1:cns),p_f_sfc(1:cns),p_p_vegtempc(1:cns),p_p_lzom(1:cns),               &
                 p_p_lzoh(1:cns),p_p_cndzmin(1:cns))
  oldval(1:cns)=p_sntemp(1:cns)
  p_sntemp(1:cns)=0.25*ctmax(1:cns)+0.75*ctmin(1:cns)
  do k=1,nfgits ! sectant
    evctx(1:cns)=evct(1:cns)
    call solverdsn(cns,evct(1:cns),p_rg_ro(1:cns),p_rg_walle(1:cns),p_rg_wallw(1:cns),p_rg_vegc(1:cns),p_rg_rs(1:cns),     &
                   p_fg_ro(1:cns),p_fg_walle(1:cns),p_fg_wallw(1:cns),p_fg_vegc(1:cns),p_fg_rs(1:cns),                     &
                   p_eg_ro(1:cns),p_eg_vegc(1:cns),p_eg_rs(1:cns),p_gasn(1:cns),p_snmelt(1:cns),p_sntemp(1:cns),           &
                   p_ro_temp(1:cns),p_ro_water(1:cns),p_ro_snow(1:cns),p_ro_den(1:cns),p_we_temp(1:cns),                   &
                   p_ww_temp(1:cns),p_v_watrc(1:cns),p_v_moistc(1:cns),p_d_canyontemp(1:cns),p_d_canyonmix(1:cns),         &
                   p_d_sig(1:cns),p_d_temp(1:cns),p_d_mixr(1:cns),p_d_topu(1:cns),p_d_roaddelta(1:cns),                    &
                   p_d_vegdeltac(1:cns),p_d_rsdelta(1:cns),p_d_netrad(1:cns),p_d_netemiss(1:cns),p_d_tranc(1:cns),         &
                   p_d_evapc(1:cns),p_d_cwa(1:cns),p_d_cra(1:cns),p_d_cwe(1:cns),p_d_cww(1:cns),p_d_crw(1:cns),            &
                   p_d_crr(1:cns),p_d_cwr(1:cns),p_d_accool(1:cns),p_d_totdepth(1:cns),                                    &
                   p_d_c1c(1:cns),p_d_acout(1:cns),p_d_traf(1:cns),p_d_canyonrgout(1:cns),p_sg_vegc(1:cns),p_sg_rs(1:cns), &
                   p_a_umag(1:cns),p_a_rho(1:cns),p_a_rg(1:cns),p_a_rnd(1:cns),p_a_snd(1:cns),ddt,p_acond_ro(1:cns),       &
                   p_acond_walle(1:cns),p_acond_wallw(1:cns),p_acond_vegc(1:cns),p_acond_rs(1:cns),p_wallpsi(1:cns),       &
                   p_roadpsi(1:cns),p_f_bldheight(1:cns),p_f_hwratio(1:cns),p_f_sigmavegc(1:cns),p_f_roademiss(1:cns),     &
                   p_f_vegemissc(1:cns),p_f_wallemiss(1:cns),p_f_ctime(1:cns),p_f_trafficfg(1:cns),p_f_rodepth(1:cns),     &
                   p_f_rolambda(1:cns),p_f_sigmabld(1:cns),p_f_bldtemp(1:cns),p_f_zovegc(1:cns),p_f_vegrlaic(1:cns),       &
                   p_f_vegrsminc(1:cns),p_f_swilt(1:cns),p_f_sfc(1:cns),p_p_vegtempc(1:cns),p_p_lzom(1:cns),               &
                   p_p_lzoh(1:cns),p_p_cndzmin(1:cns))
    evctx(1:cns)=evct(1:cns)-evctx(1:cns)
    where (abs(evctx(1:cns)).gt.tol)
      newval(1:cns)=p_sntemp(1:cns)-alpha*evct(1:cns)*(p_sntemp(1:cns)-oldval(1:cns))/evctx(1:cns)
      oldval(1:cns)=p_sntemp(1:cns)
      p_sntemp(1:cns)=newval(1:cns)
    end where
    p_sntemp(1:cns)=min(max(p_sntemp(1:cns),ctmin(1:cns)),ctmax(1:cns))      
  end do
  p_fg_rs(1:cns)=p_sg_rs(1:cns)+p_rg_rs(1:cns)-p_eg_rs(1:cns)-p_gasn(1:cns) !balance snow energy budget
  rg_road=unpack(p_rg_ro(1:cns),pgs,rg_road)
  rg_walle=unpack(p_rg_walle(1:cns),pgs,rg_walle)
  rg_wallw=unpack(p_rg_wallw(1:cns),pgs,rg_wallw)
  rg_vegc=unpack(p_rg_vegc(1:cns),pgs,rg_vegc)
  rg_rdsn=unpack(p_rg_rs(1:cns),pgs,rg_rdsn)
  fg_road=unpack(p_fg_ro(1:cns),pgs,fg_road)
  fg_walle=unpack(p_fg_walle(1:cns),pgs,fg_walle)
  fg_wallw=unpack(p_fg_wallw(1:cns),pgs,fg_wallw)
  fg_vegc=unpack(p_fg_vegc(1:cns),pgs,fg_vegc)
  fg_rdsn=unpack(p_fg_rs(1:cns),pgs,fg_rdsn)
  eg_road=unpack(p_eg_ro(1:cns),pgs,eg_road)
  eg_vegc=unpack(p_eg_vegc(1:cns),pgs,eg_vegc)
  eg_rdsn=unpack(p_eg_rs(1:cns),pgs,eg_rdsn)
  p_vegtempc=unpack(p_p_vegtempc(1:cns),pgs,p_vegtempc)
  p_lzom=unpack(p_p_lzom(1:cns),pgs,p_lzom)
  p_lzoh=unpack(p_p_lzoh(1:cns),pgs,p_lzoh)
  p_cndzmin=unpack(p_p_cndzmin(1:cns),pgs,p_cndzmin)
  gardsn=unpack(p_gasn(1:cns),pgs,gardsn)
  rdsnmelt=unpack(p_snmelt(1:cns),pgs,rdsnmelt)
  rdsntemp=unpack(p_sntemp(1:cns),pgs,rdsntemp)
  acond_road=unpack(p_acond_ro(1:cns),pgs,acond_road)
  acond_walle=unpack(p_acond_walle(1:cns),pgs,acond_walle)
  acond_wallw=unpack(p_acond_wallw(1:cns),pgs,acond_wallw)
  acond_vegc=unpack(p_acond_vegc(1:cns),pgs,acond_vegc)
  acond_rdsn=unpack(p_acond_rs(1:cns),pgs,acond_rdsn)
  d_canyontemp=unpack(p_d_canyontemp(1:cns),pgs,d_canyontemp)
  d_canyonmix=unpack(p_d_canyonmix(1:cns),pgs,d_canyonmix)
  d_topu=unpack(p_d_topu(1:cns),pgs,d_topu)
  d_netrad=unpack(p_d_netrad(1:cns),pgs,d_netrad)
  d_tranc=unpack(p_d_tranc(1:cns),pgs,d_tranc)
  d_evapc=unpack(p_d_evapc(1:cns),pgs,d_evapc)
  d_accool=unpack(p_d_accool(1:cns),pgs,d_accool)
  d_traf=unpack(p_d_traf(1:cns),pgs,d_traf)
  d_canyonrgout=unpack(p_d_canyonrgout(1:cns),pgs,d_canyonrgout)
end if
pgs=.not.pgs
cns=count(pgs)
if (cns.gt.0) then ! no road snow
  p_ro_temp(1:cns)=pack(rd_temp(:,1),pgs)
  p_ro_water(1:cns)=pack(rd_water,pgs)  
  p_ro_snow(1:cns)=pack(rd_snow,pgs)
  p_ro_den(1:cns)=pack(rd_den,pgs)
  p_we_temp(1:cns)=pack(we_temp(:,1),pgs)
  p_ww_temp(1:cns)=pack(ww_temp(:,1),pgs)
  p_v_watrc(1:cns)=pack(v_watrc,pgs)
  p_v_moistc(1:cns)=pack(v_moistc,pgs)
  p_d_sig(1:cns)=pack(d_sigd,pgs)
  p_d_temp(1:cns)=pack(d_tempc,pgs)
  p_d_mixr(1:cns)=pack(d_mixrc,pgs)
  p_d_roaddelta(1:cns)=pack(d_roaddelta,pgs)
  p_d_vegdeltac(1:cns)=pack(d_vegdeltac,pgs)
  p_d_rsdelta(1:cns)=pack(d_rdsndelta,pgs)
  p_d_netemiss(1:cns)=pack(d_netemiss,pgs)
  p_d_cwa(1:cns)=pack(d_cwa,pgs)
  p_d_cra(1:cns)=pack(d_cra,pgs)
  p_d_cwe(1:cns)=pack(d_cwe,pgs)
  p_d_cww(1:cns)=pack(d_cww,pgs)
  p_d_crw(1:cns)=pack(d_crw,pgs)
  p_d_crr(1:cns)=pack(d_crr,pgs)
  p_d_cwr(1:cns)=pack(d_cwr,pgs)
  p_d_totdepth(1:cns)=pack(d_totdepth,pgs)
  p_d_c1c(1:cns)=pack(d_c1c,pgs)
  p_d_acout(1:cns)=pack(d_acout,pgs)  
  p_sg_vegc(1:cns)=pack(sg_vegc,pgs)
  p_sg_rs(1:cns)=pack(sg_rdsn,pgs)
  p_rg_ro(1:cns)=pack(rg_road,pgs)
  p_rg_walle(1:cns)=pack(rg_walle,pgs)
  p_rg_wallw(1:cns)=pack(rg_wallw,pgs)
  p_rg_vegc(1:cns)=pack(rg_vegc,pgs)
  p_rg_rs(1:cns)=pack(rg_rdsn,pgs)
  p_fg_ro(1:cns)=pack(fg_road,pgs)
  p_fg_walle(1:cns)=pack(fg_walle,pgs)
  p_fg_wallw(1:cns)=pack(fg_wallw,pgs)
  p_fg_vegc(1:cns)=pack(fg_vegc,pgs)
  p_fg_rs(1:cns)=pack(fg_rdsn,pgs)
  p_eg_ro(1:cns)=pack(eg_road,pgs)
  p_eg_vegc(1:cns)=pack(eg_vegc,pgs)
  p_eg_rs(1:cns)=pack(eg_rdsn,pgs)
  p_a_umag(1:cns)=pack(a_umag,pgs)
  p_a_rho(1:cns)=pack(a_rho,pgs)
  p_a_rg(1:cns)=pack(a_rg,pgs)
  p_a_rnd(1:cns)=pack(a_rnd,pgs)
  p_a_snd(1:cns)=pack(a_snd,pgs)
  p_acond_ro(1:cns)=pack(acond_road,pgs)
  p_acond_walle(1:cns)=pack(acond_walle,pgs)
  p_acond_wallw(1:cns)=pack(acond_wallw,pgs)
  p_acond_vegc(1:cns)=pack(acond_vegc,pgs)
  p_acond_rs(1:cns)=pack(acond_rdsn,pgs)
  p_wallpsi(1:cns)=pack(wallpsi,pgs)
  p_roadpsi(1:cns)=pack(roadpsi,pgs)
  p_f_bldheight(1:cns)=pack(f_bldheight,pgs)
  p_f_hwratio(1:cns)=pack(f_hwratio,pgs)
  p_f_sigmavegc(1:cns)=pack(f_sigmavegc,pgs)
  p_f_roademiss(1:cns)=pack(f_roademiss,pgs)
  p_f_vegemissc(1:cns)=pack(f_vegemissc,pgs)
  p_f_wallemiss(1:cns)=pack(f_wallemiss,pgs)
  p_f_ctime(1:cns)=pack(f_ctime,pgs)
  p_f_trafficfg(1:cns)=pack(f_trafficfg,pgs)
  p_f_rodepth(1:cns)=pack(f_roaddepth(:,1),pgs)
  p_f_rolambda(1:cns)=pack(f_roadlambda(:,1),pgs)
  p_f_sigmabld(1:cns)=pack(f_sigmabld,pgs)
  p_f_bldtemp(1:cns)=pack(f_bldtemp,pgs)
  p_f_zovegc(1:cns)=pack(f_zovegc,pgs)
  p_f_vegrlaic(1:cns)=pack(f_vegrlaic,pgs)
  p_f_vegrsminc(1:cns)=pack(f_vegrsminc,pgs)
  p_f_swilt(1:cns)=pack(f_swilt,pgs)
  p_f_sfc(1:cns)=pack(f_sfc,pgs)
  p_p_vegtempc(1:cns)=pack(p_vegtempc,pgs)
  p_p_lzom(1:cns)=pack(p_lzom,pgs)
  p_p_lzoh(1:cns)=pack(p_lzoh,pgs)
  p_p_cndzmin(1:cns)=pack(p_cndzmin,pgs)
  p_sntemp(1:cns)=p_ro_temp(1:cns)
  call solverdsn(cns,evct(1:cns),p_rg_ro(1:cns),p_rg_walle(1:cns),p_rg_wallw(1:cns),p_rg_vegc(1:cns),p_rg_rs(1:cns),     &
                 p_fg_ro(1:cns),p_fg_walle(1:cns),p_fg_wallw(1:cns),p_fg_vegc(1:cns),p_fg_rs(1:cns),                     &
                 p_eg_ro(1:cns),p_eg_vegc(1:cns),p_eg_rs(1:cns),p_gasn(1:cns),p_snmelt(1:cns),p_sntemp(1:cns),           &
                 p_ro_temp(1:cns),p_ro_water(1:cns),p_ro_snow(1:cns),p_ro_den(1:cns),p_we_temp(1:cns),                   &
                 p_ww_temp(1:cns),p_v_watrc(1:cns),p_v_moistc(1:cns),p_d_canyontemp(1:cns),p_d_canyonmix(1:cns),         &
                 p_d_sig(1:cns),p_d_temp(1:cns),p_d_mixr(1:cns),p_d_topu(1:cns),p_d_roaddelta(1:cns),                    &
                 p_d_vegdeltac(1:cns),p_d_rsdelta(1:cns),p_d_netrad(1:cns),p_d_netemiss(1:cns),p_d_tranc(1:cns),         &
                 p_d_evapc(1:cns),p_d_cwa(1:cns),p_d_cra(1:cns),p_d_cwe(1:cns),p_d_cww(1:cns),p_d_crw(1:cns),            &
                 p_d_crr(1:cns),p_d_cwr(1:cns),p_d_accool(1:cns),p_d_totdepth(1:cns),                                    &
                 p_d_c1c(1:cns),p_d_acout(1:cns),p_d_traf(1:cns),p_d_canyonrgout(1:cns),p_sg_vegc(1:cns),p_sg_rs(1:cns), &
                 p_a_umag(1:cns),p_a_rho(1:cns),p_a_rg(1:cns),p_a_rnd(1:cns),p_a_snd(1:cns),ddt,p_acond_ro(1:cns),       &
                 p_acond_walle(1:cns),p_acond_wallw(1:cns),p_acond_vegc(1:cns),p_acond_rs(1:cns),p_wallpsi(1:cns),       &
                 p_roadpsi(1:cns),p_f_bldheight(1:cns),p_f_hwratio(1:cns),p_f_sigmavegc(1:cns),p_f_roademiss(1:cns),     &
                 p_f_vegemissc(1:cns),p_f_wallemiss(1:cns),p_f_ctime(1:cns),p_f_trafficfg(1:cns),p_f_rodepth(1:cns),     &
                 p_f_rolambda(1:cns),p_f_sigmabld(1:cns),p_f_bldtemp(1:cns),p_f_zovegc(1:cns),p_f_vegrlaic(1:cns),       &
                 p_f_vegrsminc(1:cns),p_f_swilt(1:cns),p_f_sfc(1:cns),p_p_vegtempc(1:cns),p_p_lzom(1:cns),               &
                 p_p_lzoh(1:cns),p_p_cndzmin(1:cns))
  p_fg_rs(1:cns)=p_sg_rs(1:cns)+p_rg_rs(1:cns)-p_eg_rs(1:cns)-p_gasn(1:cns) !balance snow energy budget
  rg_road=unpack(p_rg_ro(1:cns),pgs,rg_road)
  rg_walle=unpack(p_rg_walle(1:cns),pgs,rg_walle)
  rg_wallw=unpack(p_rg_wallw(1:cns),pgs,rg_wallw)
  rg_vegc=unpack(p_rg_vegc(1:cns),pgs,rg_vegc)
  rg_rdsn=unpack(p_rg_rs(1:cns),pgs,rg_rdsn)
  fg_road=unpack(p_fg_ro(1:cns),pgs,fg_road)
  fg_walle=unpack(p_fg_walle(1:cns),pgs,fg_walle)
  fg_wallw=unpack(p_fg_wallw(1:cns),pgs,fg_wallw)
  fg_vegc=unpack(p_fg_vegc(1:cns),pgs,fg_vegc)
  fg_rdsn=unpack(p_fg_rs(1:cns),pgs,fg_rdsn)
  eg_road=unpack(p_eg_ro(1:cns),pgs,eg_road)
  eg_vegc=unpack(p_eg_vegc(1:cns),pgs,eg_vegc)
  eg_rdsn=unpack(p_eg_rs(1:cns),pgs,eg_rdsn)
  p_vegtempc=unpack(p_p_vegtempc(1:cns),pgs,p_vegtempc)
  p_lzom=unpack(p_p_lzom(1:cns),pgs,p_lzom)
  p_lzoh=unpack(p_p_lzoh(1:cns),pgs,p_lzoh)
  p_cndzmin=unpack(p_p_cndzmin(1:cns),pgs,p_cndzmin)
  acond_road=unpack(p_acond_ro(1:cns),pgs,acond_road)
  acond_walle=unpack(p_acond_walle(1:cns),pgs,acond_walle)
  acond_wallw=unpack(p_acond_wallw(1:cns),pgs,acond_wallw)
  acond_vegc=unpack(p_acond_vegc(1:cns),pgs,acond_vegc)
  acond_rdsn=unpack(p_acond_rs(1:cns),pgs,acond_rdsn)
  d_canyontemp=unpack(p_d_canyontemp(1:cns),pgs,d_canyontemp)
  d_canyonmix=unpack(p_d_canyonmix(1:cns),pgs,d_canyonmix)
  d_topu=unpack(p_d_topu(1:cns),pgs,d_topu)
  d_netrad=unpack(p_d_netrad(1:cns),pgs,d_netrad)
  d_tranc=unpack(p_d_tranc(1:cns),pgs,d_tranc)
  d_evapc=unpack(p_d_evapc(1:cns),pgs,d_evapc)
  d_accool=unpack(p_d_accool(1:cns),pgs,d_accool)
  d_traf=unpack(p_d_traf(1:cns),pgs,d_traf)
  d_canyonrgout=unpack(p_d_canyonrgout(1:cns),pgs,d_canyonrgout)
  where (pgs)
    gardsn=0.
    rdsnmelt=0.
    rdsntemp=rd_temp(:,1)
  end where
end if
! ---------------------------------------------------------------    

! solve for roof snow temperature -------------------------------
pgs=d_rfsndelta.gt.0.
cns=count(pgs)
rg_rfsn=0.
fg_rfsn=0.
eg_rfsn=0.
garfsn=0.
rfsnmelt=0.
rfsntemp=rf_temp(:,1)
acond_rfsn=0.
if (cns.gt.0) then ! roof snow
  p_ro_temp(1:cns)=pack(rf_temp(:,1),pgs)
  p_ro_snow(1:cns)=pack(rf_snow,pgs)
  p_ro_den(1:cns)=pack(rf_den,pgs)
  p_d_dzmin(1:cns)=pack(d_rfdzmin,pgs)
  p_d_sig(1:cns)=pack(d_sigr,pgs)
  p_d_temp(1:cns)=pack(d_tempr,pgs)
  p_d_mixr(1:cns)=pack(d_mixrr,pgs)
  p_d_rsdelta(1:cns)=pack(d_rfsndelta,pgs)
  p_sg_rs(1:cns)=pack(sg_rfsn,pgs)
  p_rg_rs(1:cns)=pack(rg_rfsn,pgs)
  p_fg_rs(1:cns)=pack(fg_rfsn,pgs)
  p_eg_rs(1:cns)=pack(eg_rfsn,pgs)
  p_a_umag(1:cns)=pack(a_umag,pgs)
  p_a_rho(1:cns)=pack(a_rho,pgs)
  p_a_rg(1:cns)=pack(a_rg,pgs)
  p_a_snd(1:cns)=pack(a_snd,pgs)
  p_acond_rs(1:cns)=pack(acond_rfsn,pgs)
  p_f_rodepth(1:cns)=pack(f_roofdepth(:,1),pgs)
  p_f_rolambda(1:cns)=pack(f_rooflambda(:,1),pgs)
  ctmax(1:cns)=max(p_d_temp(1:cns),p_ro_temp(1:cns))+10. ! max roof snow temp
  ctmin(1:cns)=min(p_d_temp(1:cns),p_ro_temp(1:cns))-10. ! min roof snow temp
  p_sntemp(1:cns)=0.75*ctmax(1:cns)+0.25*ctmin(1:cns)
  call solverfsn(cns,evct(1:cns),p_rg_rs(1:cns),p_fg_rs(1:cns),p_eg_rs(1:cns),p_gasn(1:cns),p_snmelt(1:cns),       &
                 p_sntemp(1:cns),p_ro_temp(1:cns),p_ro_snow(1:cns),p_ro_den(1:cns),p_d_dzmin(1:cns),               &
                 p_d_sig(1:cns),p_d_temp(1:cns),p_d_mixr(1:cns),p_d_rsdelta(1:cns),p_sg_rs(1:cns),p_a_umag(1:cns), &
                 p_a_rho(1:cns),p_a_rg(1:cns),p_a_snd(1:cns),p_acond_rs(1:cns),p_f_rodepth(1:cns),                 &
                 p_f_rolambda(1:cns),ddt)
  oldval(1:cns)=p_sntemp(1:cns)
  p_sntemp(1:cns)=0.25*ctmax(1:cns)+0.75*ctmin(1:cns)
  do k=1,nfgits ! sectant
    evctx(1:cns)=evct(1:cns)
    call solverfsn(cns,evct(1:cns),p_rg_rs(1:cns),p_fg_rs(1:cns),p_eg_rs(1:cns),p_gasn(1:cns),p_snmelt(1:cns),       &
                   p_sntemp(1:cns),p_ro_temp(1:cns),p_ro_snow(1:cns),p_ro_den(1:cns),p_d_dzmin(1:cns),               &
                   p_d_sig(1:cns),p_d_temp(1:cns),p_d_mixr(1:cns),p_d_rsdelta(1:cns),p_sg_rs(1:cns),p_a_umag(1:cns), &
                   p_a_rho(1:cns),p_a_rg(1:cns),p_a_snd(1:cns),p_acond_rs(1:cns),p_f_rodepth(1:cns),                 &
                   p_f_rolambda(1:cns),ddt)
    evctx(1:cns)=evct(1:cns)-evctx(1:cns)
    where (abs(evctx(1:cns)).gt.tol)
      newval(1:cns)=p_sntemp(1:cns)-alpha*evct(1:cns)*(p_sntemp(1:cns)-oldval(1:cns))/evctx(1:cns)
      oldval(1:cns)=p_sntemp(1:cns)
      p_sntemp(1:cns)=newval(1:cns)
    end where
    p_sntemp(1:cns)=min(max(p_sntemp(1:cns),ctmin(1:cns)),ctmax(1:cns))      
  end do
  p_fg_rs(1:cns)=p_sg_rs(1:cns)+p_rg_rs(1:cns)-p_eg_rs(1:cns)-p_gasn(1:cns) !balance snow energy budget
  rg_rfsn=unpack(p_rg_rs(1:cns),pgs,rg_rfsn)
  fg_rfsn=unpack(p_fg_rs(1:cns),pgs,fg_rfsn)
  eg_rfsn=unpack(p_eg_rs(1:cns),pgs,eg_rfsn)
  garfsn=unpack(p_gasn(1:cns),pgs,garfsn)
  rfsnmelt=unpack(p_snmelt(1:cns),pgs,rfsnmelt)
  acond_rfsn=unpack(p_acond_rs(1:cns),pgs,acond_rfsn)
end if
!---------------------------------------------------------------- 

! calculate green roof sensible and latent heat fluxes
call solverfveg(rg_vegr,fg_vegr,eg_vegr,acond_vegr,a_rg,a_umag,a_rho,a_rnd,sg_vegr,d_tempr,d_mixrr,d_rfdzmin, &
                d_tranr,d_evapr,d_c1r,d_sigr,d_vegdeltar,ddt)

! calculate roof sensible and latent heat fluxes (without snow)
rg_roof=f_roofemiss*(a_rg-sbconst*rf_temp(:,1)**4)
! a is a dummy variable for lzomroof
a=log(d_rfdzmin/zoroof)
! xe is a dummy variable for lzohroof
xe=2.3+a
xw=d_roofdelta*qsatr ! roof surface mixing ratio
dts=rf_temp(:,1)*(1.+0.61*xw)
dtt=d_tempr*(1.+0.61*d_mixrr)
! Assume zot=0.1*zom (i.e., Kanda et al 2007, small experiment)
! n is a dummy variable for cd
call getinvres(ufull,acond_roof,n,z_on_l,xe,a,d_rfdzmin,dts,dtt,a_umag,1)
acond_roof=acond_roof/a_umag
fg_roof=aircp*a_rho*(rf_temp(:,1)-d_tempr)*acond_roof*a_umag
call getqsat(ufull,qsatr,dts,d_sigr)
where (qsatr.lt.d_mixrr)
  eg_roof=lv*a_rho*(qsatr-d_mixrr)*acond_roof*a_umag
elsewhere
  eg_roof=lv*min(a_rho*d_roofdelta*(qsatr-d_mixrr)*acond_roof*a_umag,rf_water/ddt+a_rnd+rfsnmelt)
end where

! tridiagonal solver coefficents for calculating roof, road and wall temperatures
do k=2,3
  ggaroof(:,k)=-2./(f_roofdepth(:,k-1)/f_rooflambda(:,k-1)+f_roofdepth(:,k)/f_rooflambda(:,k))
  ggawall(:,k)=-2./(f_walldepth(:,k-1)/f_walllambda(:,k-1)+f_walldepth(:,k)/f_walllambda(:,k))
  ggaroad(:,k)=-2./(f_roaddepth(:,k-1)/f_roadlambda(:,k-1)+f_roaddepth(:,k)/f_roadlambda(:,k))
end do
ggbroof(:,1)=2./(f_roofdepth(:,1)/f_rooflambda(:,1)+f_roofdepth(:,2)/f_rooflambda(:,2))+f_roofcp(:,1)*f_roofdepth(:,1)/ddt
ggbwall(:,1)=2./(f_walldepth(:,1)/f_walllambda(:,1)+f_walldepth(:,2)/f_walllambda(:,2))+f_wallcp(:,1)*f_walldepth(:,1)/ddt
ggbroad(:,1)=2./(f_roaddepth(:,1)/f_roadlambda(:,1)+f_roaddepth(:,2)/f_roadlambda(:,2))+f_roadcp(:,1)*f_roaddepth(:,1)/ddt
ggbroof(:,2)=2./(f_roofdepth(:,1)/f_rooflambda(:,1)+f_roofdepth(:,2)/f_rooflambda(:,2)) &
             +2./(f_roofdepth(:,2)/f_rooflambda(:,2)+f_roofdepth(:,3)/f_rooflambda(:,3)) &
             +f_roofcp(:,2)*f_roofdepth(:,2)/ddt
ggbwall(:,2)=2./(f_walldepth(:,1)/f_walllambda(:,1)+f_walldepth(:,2)/f_walllambda(:,2)) &
             +2./(f_walldepth(:,2)/f_walllambda(:,2)+f_walldepth(:,3)/f_walllambda(:,3)) &
             +f_wallcp(:,2)*f_walldepth(:,2)/ddt
ggbroad(:,2)=2./(f_roaddepth(:,1)/f_roadlambda(:,1)+f_roaddepth(:,2)/f_roadlambda(:,2)) &
             +2./(f_roaddepth(:,2)/f_roadlambda(:,2)+f_roaddepth(:,3)/f_roadlambda(:,3)) &
             +f_roadcp(:,2)*f_roaddepth(:,2)/ddt
ggbroof(:,3)=2./(f_roofdepth(:,2)/f_rooflambda(:,2)+f_roofdepth(:,3)/f_rooflambda(:,3)) &
             +f_roofcp(:,3)*f_roofdepth(:,3)/ddt
ggbwall(:,3)=2./(f_walldepth(:,2)/f_walllambda(:,2)+f_walldepth(:,3)/f_walllambda(:,3)) &
             +f_wallcp(:,3)*f_walldepth(:,3)/ddt
ggbroad(:,3)=2./(f_roaddepth(:,2)/f_roadlambda(:,2)+f_roaddepth(:,3)/f_roadlambda(:,3))+f_roadcp(:,3)*f_roaddepth(:,3)/ddt
do k=1,2
  ggcroof(:,k)=-2./(f_roofdepth(:,k)/f_rooflambda(:,k)+f_roofdepth(:,k+1)/f_rooflambda(:,k+1))
  ggcwall(:,k)=-2./(f_walldepth(:,k)/f_walllambda(:,k)+f_walldepth(:,k+1)/f_walllambda(:,k+1))
  ggcroad(:,k)=-2./(f_roaddepth(:,k)/f_roadlambda(:,k)+f_roaddepth(:,k+1)/f_roadlambda(:,k+1))
end do
garoof(:,1)=(1.-d_rfsndelta)*(sg_roof+rg_roof-fg_roof-eg_roof) &
            +d_rfsndelta*garfsn+rf_temp(:,1)*f_roofcp(:,1)*f_roofdepth(:,1)/ddt
gawalle(:,1)=sg_walle+rg_walle-fg_walle+we_temp(:,1)*f_wallcp(:,1)*f_walldepth(:,1)/ddt
gawallw(:,1)=sg_wallw+rg_wallw-fg_wallw+ww_temp(:,1)*f_wallcp(:,1)*f_walldepth(:,1)/ddt
garoad(:,1)=(1.-d_rdsndelta)*(sg_road+rg_road-fg_road-eg_road) &
            +d_rdsndelta*gardsn+rd_temp(:,1)*f_roadcp(:,1)*f_roaddepth(:,1)/ddt
garoof(:,2)=rf_temp(:,2)*f_roofcp(:,2)*f_roofdepth(:,2)/ddt
gawalle(:,2)=we_temp(:,2)*f_wallcp(:,2)*f_walldepth(:,2)/ddt
gawallw(:,2)=ww_temp(:,2)*f_wallcp(:,2)*f_walldepth(:,2)/ddt
garoad(:,2)=rd_temp(:,2)*f_roadcp(:,2)*f_roaddepth(:,2)/ddt
garoof(:,3)=rf_temp(:,3)*f_roofcp(:,3)*f_roofdepth(:,3)/ddt-gflxroof
gawalle(:,3)=we_temp(:,3)*f_wallcp(:,3)*f_walldepth(:,3)/ddt-gflxwalle
gawallw(:,3)=ww_temp(:,3)*f_wallcp(:,3)*f_walldepth(:,3)/ddt-gflxwallw
garoad(:,3)=rd_temp(:,3)*f_roadcp(:,3)*f_roaddepth(:,3)/ddt
    
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
rf_temp(:,3)=garoof(:,3)/ggbroof(:,3)
we_temp(:,3)=gawalle(:,3)/ggbwall(:,3)
ww_temp(:,3)=gawallw(:,3)/ggbwall(:,3)
rd_temp(:,3)=garoad(:,3)/ggbroad(:,3)
do ii=2,1,-1
  rf_temp(:,ii)=(garoof(:,ii)-ggcroof(:,ii)*rf_temp(:,ii+1))/ggbroof(:,ii)
  we_temp(:,ii)=(gawalle(:,ii)-ggcwall(:,ii)*we_temp(:,ii+1))/ggbwall(:,ii)
  ww_temp(:,ii)=(gawallw(:,ii)-ggcwall(:,ii)*ww_temp(:,ii+1))/ggbwall(:,ii)
  rd_temp(:,ii)=(garoad(:,ii)-ggcroad(:,ii)*rd_temp(:,ii+1))/ggbroad(:,ii)
end do

if (any(f_sigmavegc.gt.0.)) then ! in-canyon vegetation
  ! calculate water for canyon surfaces
  n=max(a_rnd-d_evapc/lv-max(maxvwatf*f_vegrlaic-v_watrc,0.)/ddt,0.) ! rainfall reaching the soil under vegetation
  ! note that since sigmaf=1, then there is no soil evaporation, only transpiration.  Evaporation only occurs from water on leafs.
  v_moistc=v_moistc+ddt*(n+rdsnmelt*rd_den/waterden-d_tranc/lv)/(waterden*d_totdepth) ! soil
  v_watrc=v_watrc+ddt*(a_rnd-d_evapc/lv)                                              ! leaf
  v_watrc=min(max(v_watrc,0.),maxvwatf*f_vegrlaic)
  if (wbrelaxc.eq.1) then
    ! increase soil moisture for irrigation 
    v_moistc=v_moistc+max(0.75*f_swilt+0.25*f_sfc-v_moistc,0.)/(86400./ddt+1.) ! 24h e-fold time
  end if
else
  v_moistc=f_swilt
  v_watrc=0.
end if
if (any(f_sigmavegr.gt.0.)) then ! green roof
  ! calculate water for roof surface
  n=max(a_rnd-d_evapr/lv-max(maxvwatf*f_vegrlair-v_watrr,0.)/ddt,0.) ! rainfall reaching the soil under vegetation
  ! note that since sigmaf=1, then there is no soil evaporation, only transpiration.  Evaporation only occurs from water on leafs.
  v_moistr=v_moistr+ddt*(n+rfsnmelt*rf_den/waterden-d_tranr/lv)/(waterden*f_vegdepthr) ! soil
  v_watrr=v_watrr+ddt*(a_rnd-d_evapr/lv)                                               ! leaf
  v_watrr=min(max(v_watrr,0.),maxvwatf*f_vegrlair)
  if (wbrelaxr.eq.1) then
    ! increase soil moisture for irrigation 
    v_moistr=v_moistr+max(0.75*f_swilt+0.25*f_sfc-v_moistr,0.)/(86400./ddt+1.) ! 24h e-fold time
  end if
else
  v_moistr=f_swilt
  v_watrr=0.
end if
rf_water=rf_water+ddt*(a_rnd-eg_roof/lv+rfsnmelt)
rd_water=rd_water+ddt*(a_rnd-eg_road/lv+rdsnmelt)

! calculate snow
rf_snow=rf_snow+ddt*(a_snd-eg_rfsn/lv-rfsnmelt)
rd_snow=rd_snow+ddt*(a_snd-eg_rdsn/lv-rdsnmelt)
rf_den=rf_den+(maxsnowden-rf_den)/(0.24/(86400.*ddt)+1.)
rd_den=rd_den+(maxsnowden-rd_den)/(0.24/(86400.*ddt)+1.)
rf_alpha=rf_alpha+(minsnowalpha-rf_alpha)/(0.24/(86400.*ddt)+1.)
rd_alpha=rd_alpha+(minsnowalpha-rd_alpha)/(0.24/(86400.*ddt)+1.)

! calculate runoff (v_watrc runoff already accounted for in precip reaching canyon floor)
u_rn=0.
u_rn=u_rn+max(rf_water-maxrfwater,0.)*f_sigmabld
u_rn=u_rn+max(rd_water-maxrdwater,0.)*(1.-d_rdsndelta)*(1.-f_sigmavegc)
u_rn=u_rn+max(rf_snow-maxrfsn,0.)*f_sigmabld
u_rn=u_rn+max(rd_snow-maxrdsn,0.)*d_rdsndelta
u_rn=u_rn+max(v_moistc-f_ssat,0.)*waterden*d_totdepth*(1.-d_rdsndelta)*f_sigmavegc
u_rn=u_rn+max(v_moistr-f_ssat,0.)*waterden*f_vegdepthr*f_sigmavegr

! limit temperatures to sensible values
do ii=1,3
  rf_temp(:,ii)=min(max(rf_temp(:,ii),200.),400.)
  we_temp(:,ii)=min(max(we_temp(:,ii),200.),400.)
  ww_temp(:,ii)=min(max(ww_temp(:,ii),200.),400.)
  rd_temp(:,ii)=min(max(rd_temp(:,ii),200.),400.)
end do
v_moistc=min(max(v_moistc,f_swilt),f_ssat)
v_moistr=min(max(v_moistr,f_swilt),f_ssat)
rf_water=min(max(rf_water,0.),maxrfwater)
rd_water=min(max(rd_water,0.),maxrdwater)
v_watrc=min(max(v_watrc,0.),maxvwatf*f_vegrlaic)
v_watrr=min(max(v_watrr,0.),maxvwatf*f_vegrlair)
rf_snow=min(max(rf_snow,0.),maxrfsn)
rd_snow=min(max(rd_snow,0.),maxrdsn)
rf_den=min(max(rf_den,minsnowden),maxsnowden)
rd_den=min(max(rd_den,minsnowden),maxsnowden)
rf_alpha=min(max(rf_alpha,minsnowalpha),maxsnowalpha)
rd_alpha=min(max(rd_alpha,minsnowalpha),maxsnowalpha)

! combine snow and snow-free tiles
d_roofrgout=a_rg-d_rfsndelta*rg_rfsn-(1.-d_rfsndelta)*((1.-f_sigmavegr)*rg_roof+f_sigmavegr*rg_vegr)
fg_roof=d_rfsndelta*fg_rfsn+(1.-d_rfsndelta)*((1.-f_sigmavegr)*fg_roof+f_sigmavegr*fg_vegr)
eg_roof=d_rfsndelta*eg_rfsn+(1.-d_rfsndelta)*((1.-f_sigmavegr)*eg_roof+f_sigmavegr*eg_vegr)
fgtop=d_rdsndelta*fg_rdsn+(1.-d_rdsndelta)*((1.-f_sigmavegc)*fg_road+f_sigmavegc*fg_vegc) &
      +f_hwratio*(fg_walle+fg_wallw)+d_traf+d_accool
egtop=d_rdsndelta*eg_rdsn+(1.-d_rdsndelta)*((1.-f_sigmavegc)*eg_road+f_sigmavegc*eg_vegc)

! calculate longwave, sensible heat latent heat outputs
! estimate surface temp from outgoing longwave radiation
u_ts=((f_sigmabld*d_roofrgout+(1.-f_sigmabld)*d_canyonrgout)/sbconst)**0.25
u_fg=f_sigmabld*fg_roof+(1.-f_sigmabld)*fgtop+f_industryfg
u_eg=f_sigmabld*eg_roof+(1.-f_sigmabld)*egtop

! calculate surface water (i.e., bc for RH) output
n=max(min((v_moistr-f_swilt)/(f_sfc-f_swilt),1.),0.) ! veg wetfac (see sflux.f or cable_canopy.f90)
u_wf=f_sigmabld*(1.-d_rfsndelta)*((1.-f_sigmavegr)*d_roofdelta+f_sigmavegr*((1.-d_vegdeltar)*n+d_vegdeltar))
n=max(min((v_moistc-f_swilt)/(f_sfc-f_swilt),1.),0.) ! veg wetfac (see sflux.f or cable_canopy.f90)
u_wf=u_wf+(1.-f_sigmabld)*(1.-d_rdsndelta)*((1.-f_sigmavegc)*d_roaddelta+f_sigmavegc*((1.-d_vegdeltac)*n+d_vegdeltac))

! (re)calculate heat roughness length for MOST (diagnostic only)
call getqsat(ufull,a,u_ts,d_sigd)
a=a*u_wf
dts=u_ts*(1.+0.61*a)
dtt=d_tempc*(1.+0.61*d_mixrc)
select case(zohmeth)
  case(0) ! Use veg formulation
    p_lzoh=2.3+p_lzom
    call getinvres(ufull,p_cdtq,p_cduv,z_on_l,p_lzoh,p_lzom,p_cndzmin,dts,dtt,a_umag,1)
  case(1) ! Use Kanda parameterisation
    p_lzoh=2.3+p_lzom ! replaced in getlna
    call getinvres(ufull,p_cdtq,p_cduv,z_on_l,p_lzoh,p_lzom,p_cndzmin,dts,dtt,a_umag,2)
  case(2) ! Use Kanda parameterisation
    p_lzoh=6.+p_lzom
    call getinvres(ufull,p_cdtq,p_cduv,z_on_l,p_lzoh,p_lzom,p_cndzmin,dts,dtt,a_umag,4)
end select

! calculate screen level diagnostics
call scrncalc(a_mixr,a_umag,a_temp,u_ts,d_tempc,d_mixrc,d_rdsndelta,d_roaddelta,d_vegdeltac,d_sigd,a,rdsntemp,zonet)

!write(61,'(6F10.4)') fg_walle(iqut),fg_wallw(iqut),fg_road(iqut),fg_vegc(iqut),d_traf(iqut),d_accool(iqut)
!write(61,'(6F10.6)') a_umag(iqut),acond_road(iqut),acond_walle(iqut),acond_wallw(iqut),acond_veg(iqut)

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
! Estimate saturation mixing ratio

subroutine getqsat(ifull,qsat,temp,ps)

implicit none

integer, intent(in) :: ifull
real, dimension(ifull), intent(in) :: temp,ps
real, dimension(ifull), intent(out) :: qsat
real, dimension(0:220) :: table
real, dimension(ifull) :: esatf,tdiff,rx
integer, dimension(ifull) :: ix

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
qsat=0.622*esatf/max(ps-esatf,0.1)

return
end subroutine getqsat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Interface for calcuating ustar and thetastar

subroutine getinvres(cn,invres,cd,z_on_l,olzoh,ilzom,zmin,sthetav,thetav,a_umag,mode)

implicit none

integer, intent(in) :: cn,mode
real, dimension(cn), intent(in) :: ilzom,zmin,sthetav,thetav
real, dimension(cn), intent(out) :: invres,cd,z_on_l
real, dimension(cn), intent(inout) :: olzoh
real, dimension(cn) :: lna,thetavstar,integralh
real, dimension(cn), intent(in) :: a_umag

lna=olzoh-ilzom
call dyerhicks(cn,integralh,z_on_l,cd,thetavstar,thetav,sthetav,a_umag,zmin,ilzom,lna,mode)
invres=vkar*sqrt(cd)*a_umag/integralh
olzoh=lna+ilzom

return
end subroutine getinvres

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate stability functions using Dyerhicks

subroutine dyerhicks(cn,integralh,z_on_l,cd,thetavstar,thetav,sthetav,umag,zmin,ilzom,lna,mode)

implicit none

integer, intent(in) :: cn,mode
integer ic
real, dimension(cn), intent(in) :: thetav,sthetav,umag,zmin,ilzom
real, dimension(cn), intent(inout) :: lna
real, dimension(cn), intent(out) :: cd,thetavstar
real, dimension(cn), intent(out) :: integralh,z_on_l
real, dimension(cn) :: z0_on_l,zt_on_l,olzoh
real, dimension(cn) :: pm0,ph0,pm1,ph1,integralm
real, parameter :: aa1 = 3.8
real, parameter :: bb1 = 0.5
real, parameter :: cc1 = 0.3

cd=(vkar/ilzom)**2                         ! first guess
call getlna(cn,lna,cd,umag,zmin,ilzom,mode)
olzoh=ilzom+lna
integralh=sqrt(cd)*ilzom*olzoh/vkar        ! first guess
thetavstar=vkar*(thetav-sthetav)/integralh ! first guess

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
  call getlna(cn,lna,cd,umag,zmin,ilzom,mode)
end do

return
end subroutine dyerhicks

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Estimate roughness length for heat
!

subroutine getlna(cn,lna,cd,umag,zmin,ilzom,mode)

implicit none

integer, intent(in) :: cn,mode
real, dimension(cn), intent(out) :: lna
real, dimension(cn), intent(in) :: cd,umag,zmin,ilzom
real, dimension(cn) :: re
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

subroutine getswcoeff(cn,sg_roof,sg_vegr,sg_road,sg_walle,sg_wallw,sg_vegc,sg_rfsn,sg_rdsn,wallpsi,roadpsi,if_hwratio,if_vangle, &
                      if_hangle,if_fbeam,if_sigmavegc,if_roadalpha,if_vegalphac,if_wallalpha,ird_alpha,rdsndelta)

implicit none

integer, intent(in) :: cn
integer k
real, dimension(cn), intent(in) :: rdsndelta,ird_alpha
real, dimension(cn), intent(out) :: wallpsi,roadpsi
real, dimension(cn) :: thetazero,walles,wallws,roads,ta,tc,xa,ya,roadnetalpha
real, dimension(cn) :: nwalles,nwallws,nroads
real, dimension(cn), intent(in) :: if_hwratio,if_vangle,if_hangle,if_fbeam,if_sigmavegc,if_roadalpha,if_vegalphac,if_wallalpha
real, dimension(cn), intent(out) :: sg_roof,sg_vegr,sg_road,sg_walle,sg_wallw,sg_vegc,sg_rfsn,sg_rdsn

wallpsi=0.5*(if_hwratio+1.-sqrt(if_hwratio*if_hwratio+1.))/if_hwratio
roadpsi=sqrt(if_hwratio*if_hwratio+1.)-if_hwratio

! integrate through 180 instead of 360
where (if_vangle.ge.0.5*pi)
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
roadnetalpha=rdsndelta*ird_alpha+(1.-rdsndelta)*((1.-if_sigmavegc)*if_roadalpha+if_sigmavegc*if_vegalphac)
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

subroutine getlwcoeff(ifull,d_netemiss,d_cwa,d_cra,d_cwe,d_cww,d_crw,d_crr,d_cwr,d_rdsndelta,wallpsi,roadpsi,if_sigmavegc, &
                      if_roademiss,if_vegemissc,if_wallemiss)

implicit none

integer, intent(in) :: ifull
integer k
real, dimension(ifull) :: wallpsi,roadpsi
real, dimension(ifull) :: rcwa,rcra,rcwe,rcww,rcrw,rcrr,rcwr
real, dimension(ifull) :: ncwa,ncra,ncwe,ncww,ncrw,ncrr,ncwr
real, dimension(ifull), intent(inout) :: d_netemiss,d_cwa,d_cra,d_cwe,d_cww,d_crw,d_crr,d_cwr,d_rdsndelta
real, dimension(ifull), intent(in) :: if_sigmavegc,if_roademiss,if_vegemissc,if_wallemiss

d_netemiss=d_rdsndelta*snowemiss+(1.-d_rdsndelta)*((1.-if_sigmavegc)*if_roademiss+if_sigmavegc*if_vegemissc)
d_cwa=wallpsi
d_cra=roadpsi
d_cwe=0.
d_cww=1.-2.*wallpsi
d_crw=0.5*(1.-roadpsi)
d_crr=0.
d_cwr=wallpsi
rcwa=d_cwa
rcra=d_cra
rcwe=d_cwe
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
  d_cwe=d_cwe+rcwe
  d_cww=d_cww+rcww
  d_crw=d_crw+rcrw
  d_cwr=d_cwr+rcwr
  d_crr=d_crr+rcrr
end do

end subroutine getlwcoeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solve for roof snow temperature

subroutine solverfsn(cn,evct,rg_rfsn,fg_rfsn,eg_rfsn,garfsn,rfsnmelt,rfsntemp,irf_temp,irf_snow,irf_den,d_rfdzmin,d_sigr,d_tempr, &
                     d_mixrr,d_rfsndelta,sg_rfsn,a_umag,a_rg,a_rho,a_snd,acond_rfsn,if_roofdepth,if_rooflambda,ddt)

implicit none

integer, intent(in) :: cn
real, intent(in) :: ddt
real, dimension(cn), intent(out) :: evct,rfsnmelt,garfsn
real, dimension(cn), intent(in) :: rfsntemp
real, dimension(cn) :: cd,rfsnqsat,lzotdum,lzosnow,sndepth,snlambda,ldratio,z_on_l
real, dimension(cn) :: dts,dtt
real, dimension(cn), intent(in) :: sg_rfsn
real, dimension(cn), intent(inout) :: rg_rfsn,fg_rfsn,eg_rfsn,acond_rfsn
real, dimension(cn), intent(in) :: a_umag,a_rg,a_rho,a_snd
real, dimension(cn), intent(in) :: d_rfdzmin,d_sigr,d_tempr,d_mixrr,d_rfsndelta
real, dimension(cn), intent(in) :: irf_temp,irf_snow,irf_den
real, dimension(cn), intent(in) :: if_roofdepth,if_rooflambda

! snow conductance
sndepth=irf_snow*waterden/irf_den
snlambda=icelambda*(irf_den/waterden)**1.88
ldratio=0.5*(sndepth/snlambda+if_roofdepth/if_rooflambda)

! Update roof snow energy budget
lzosnow=log(d_rfdzmin/zosnow)
call getqsat(cn,rfsnqsat,rfsntemp,d_sigr)
lzotdum=2.3+lzosnow
dts=rfsntemp*(1.+0.61*rfsnqsat)
dtt=d_tempr*(1.+0.61*d_mixrr)
call getinvres(cn,acond_rfsn,cd,z_on_l,lzotdum,lzosnow,d_rfdzmin,dts,dtt,a_umag,1)
acond_rfsn=acond_rfsn/a_umag
rfsnmelt=d_rfsndelta*max(0.,rfsntemp-273.16)/(icecp*irf_den*lf*ddt) 
rg_rfsn=snowemiss*(a_rg-sbconst*rfsntemp**4)
fg_rfsn=aircp*a_rho*(rfsntemp-d_tempr)*acond_rfsn*a_umag
eg_rfsn=lv*min(a_rho*d_rfsndelta*max(0.,rfsnqsat-d_mixrr)*acond_rfsn*a_umag,irf_snow/ddt+a_snd-rfsnmelt)
garfsn=(rfsntemp-irf_temp)/ldratio
evct=sg_rfsn+rg_rfsn-fg_rfsn-eg_rfsn*ls/lv-garfsn

return
end subroutine solverfsn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solve for road snow temperature (includes vegetation canopy temperature and canyon temperature)

subroutine solverdsn(cn,evct,rg_road,rg_walle,rg_wallw,rg_vegc,rg_rdsn,fg_road,fg_walle,fg_wallw,fg_vegc,fg_rdsn,eg_road,        &
                     eg_vegc,eg_rdsn,gardsn,rdsnmelt,rdsntemp,ird_temp,ird_water,ird_snow,ird_den,iwe_temp,iww_temp,iv_watrc,    &
                     iv_moistc,d_canyontemp,d_canyonmix,d_sigd,d_tempc,d_mixrc,d_topu,d_roaddelta,d_vegdeltac,d_rdsndelta,       &
                     d_netrad,d_netemiss,d_tranc,d_evapc,d_cwa,d_cra,d_cwe,d_cww,d_crw,d_crr,d_cwr,d_accool,                     &
                     d_totdepth,d_c1c,d_acout,d_traf,d_canyonrgout,sg_vegc,sg_rdsn,a_umag,a_rho,a_rg,a_rnd,a_snd,ddt,acond_road, &
                     acond_walle,acond_wallw,acond_vegc,acond_rdsn,wallpsi,roadpsi,if_bldheight,if_hwratio,if_sigmavegc,         &
                     if_roademiss,if_vegemissc,if_wallemiss,if_ctime,if_trafficfg,if_roaddepth,if_roadlambda,if_sigmabld,        &
                     if_bldtemp,if_zovegc,if_vegrlaic,if_vegrsminc,if_swilt,if_sfc,ip_vegtempc,ip_lzom,ip_lzoh,ip_cndzmin)

implicit none

integer, intent(in) :: cn
integer k,ic
real, intent(in) :: ddt
real, dimension(cn), intent(out) :: evct,gardsn,rdsnmelt
real, dimension(cn), intent(in) :: rdsntemp,wallpsi,roadpsi
real, dimension(cn) :: ctmax,ctmin,topinvres
real, dimension(cn) :: ldratio,sndepth,snlambda
real, dimension(cn) :: oldval,newval,evctveg,evctx
real, dimension(cn), intent(in) :: sg_vegc,sg_rdsn
real, dimension(cn), intent(inout) :: rg_road,rg_walle,rg_wallw,rg_vegc,rg_rdsn
real, dimension(cn), intent(inout) :: fg_road,fg_walle,fg_wallw,fg_vegc,fg_rdsn
real, dimension(cn), intent(inout) :: eg_road,eg_vegc,eg_rdsn
real, dimension(cn), intent(inout) :: acond_road,acond_walle,acond_wallw,acond_vegc,acond_rdsn
real, dimension(cn), intent(in) :: a_umag,a_rho,a_rg,a_rnd,a_snd
real, dimension(cn), intent(inout) :: d_canyontemp,d_canyonmix,d_sigd,d_tempc,d_mixrc,d_topu,d_roaddelta,d_vegdeltac,d_rdsndelta
real, dimension(cn), intent(inout) :: d_netrad,d_netemiss,d_tranc,d_evapc,d_cwa,d_cra,d_cwe,d_cww,d_crw,d_crr,d_cwr,d_accool
real, dimension(cn), intent(inout) :: d_totdepth,d_c1c,d_acout,d_traf,d_canyonrgout
real, dimension(cn), intent(in) :: ird_temp,ird_water,ird_snow,ird_den
real, dimension(cn), intent(in) :: iwe_temp,iww_temp
real, dimension(cn), intent(in) :: iv_watrc,iv_moistc
real, dimension(cn), intent(in) :: if_bldheight,if_hwratio,if_sigmavegc,if_roademiss,if_vegemissc,if_wallemiss
real, dimension(cn), intent(in) :: if_ctime,if_trafficfg,if_roaddepth,if_roadlambda,if_sigmabld,if_bldtemp
real, dimension(cn), intent(in) :: if_zovegc,if_vegrlaic,if_vegrsminc,if_swilt,if_sfc
real, dimension(cn), intent(inout) :: ip_vegtempc,ip_lzom,ip_lzoh,ip_cndzmin

! traffic sensible heat flux
call gettraffic(cn,d_traf,if_ctime,if_trafficfg)
d_traf=d_traf/(1.-if_sigmabld)

! snow conductance
sndepth=ird_snow*waterden/ird_den
snlambda=icelambda*(ird_den/waterden)**1.88
ldratio=0.5*(sndepth/snlambda+if_roaddepth/if_roadlambda)

! solve for vegetation canopy and canyon temperature
! (we previously employed a multi-variable Newton-Raphson method here.  However, it is difficult to handle the
!  dependence of the wind speed at the canyon top (i.e., affecting aerodynamical resistances)
!  when the wind speed at the canyon top depends on the stability functions.  This multiply
!  nested version is slower, but more robust).
ctmax=max(d_tempc,ird_temp,iwe_temp,iww_temp,rdsntemp)+10. ! max temp
ctmin=min(d_tempc,ird_temp,iwe_temp,iww_temp,rdsntemp)-10. ! min temp
if (any(if_sigmavegc.gt.0.)) then ! in-canyon vegetation
  ip_vegtempc=0.6*ctmax+0.4*ctmin
  call solverdvg(cn,evctveg,rg_road,rg_walle,rg_wallw,rg_vegc,rg_rdsn,fg_road,fg_walle,fg_wallw,fg_vegc,fg_rdsn,eg_road,eg_vegc, &
                 eg_rdsn,gardsn,rdsnmelt,topinvres,rdsntemp,ird_temp,ird_water,ird_snow,ird_den,iwe_temp,iww_temp,iv_watrc,      &
                 iv_moistc,d_canyontemp,d_canyonmix,d_sigd,d_tempc,d_mixrc,d_topu,d_roaddelta,d_vegdeltac,d_rdsndelta,d_netrad,  &
                 d_netemiss,d_tranc,d_evapc,d_cwa,d_cra,d_cwe,d_cww,d_crw,d_crr,d_cwr,d_accool,d_totdepth,d_c1c,                 &
                 d_acout,d_canyonrgout,sg_vegc,a_umag,a_rho,a_rg,a_rnd,a_snd,ddt,acond_road,acond_walle,acond_wallw,acond_vegc,  &
                 acond_rdsn,wallpsi,roadpsi,if_bldheight,if_hwratio,if_sigmavegc,if_roademiss,if_vegemissc,if_wallemiss,         &
                 if_bldtemp,if_zovegc,if_vegrlaic,if_vegrsminc,if_swilt,if_sfc,ip_vegtempc,ip_lzom,ip_lzoh,ip_cndzmin,           &
                 d_traf,ldratio)
  oldval=ip_vegtempc
  ip_vegtempc=0.4*ctmax+0.6*ctmin
  do k=1,nfgits ! Sectant
    evctx=evctveg
    call solverdvg(cn,evctveg,rg_road,rg_walle,rg_wallw,rg_vegc,rg_rdsn,fg_road,fg_walle,fg_wallw,fg_vegc,fg_rdsn,eg_road,eg_vegc, &
                   eg_rdsn,gardsn,rdsnmelt,topinvres,rdsntemp,ird_temp,ird_water,ird_snow,ird_den,iwe_temp,iww_temp,iv_watrc,      &
                   iv_moistc,d_canyontemp,d_canyonmix,d_sigd,d_tempc,d_mixrc,d_topu,d_roaddelta,d_vegdeltac,d_rdsndelta,d_netrad,  &
                   d_netemiss,d_tranc,d_evapc,d_cwa,d_cra,d_cwe,d_cww,d_crw,d_crr,d_cwr,d_accool,d_totdepth,d_c1c,                 &
                   d_acout,d_canyonrgout,sg_vegc,a_umag,a_rho,a_rg,a_rnd,a_snd,ddt,acond_road,acond_walle,acond_wallw,acond_vegc,  &
                   acond_rdsn,wallpsi,roadpsi,if_bldheight,if_hwratio,if_sigmavegc,if_roademiss,if_vegemissc,if_wallemiss,         &
                   if_bldtemp,if_zovegc,if_vegrlaic,if_vegrsminc,if_swilt,if_sfc,ip_vegtempc,ip_lzom,ip_lzoh,ip_cndzmin,           &
                   d_traf,ldratio)
    evctx=evctveg-evctx
    where (abs(evctx).gt.tol)  
      newval=ip_vegtempc-alpha*evctveg*(ip_vegtempc-oldval)/evctx
      oldval=ip_vegtempc
      ip_vegtempc=newval
    end where
    ip_vegtempc=min(max(ip_vegtempc,ctmin),ctmax)
  end do
else ! no in-canyon vegetation
  ip_vegtempc=0.5*ctmax+0.5*ctmin
  call solverdvg(cn,evctveg,rg_road,rg_walle,rg_wallw,rg_vegc,rg_rdsn,fg_road,fg_walle,fg_wallw,fg_vegc,fg_rdsn,eg_road,eg_vegc, &
                 eg_rdsn,gardsn,rdsnmelt,topinvres,rdsntemp,ird_temp,ird_water,ird_snow,ird_den,iwe_temp,iww_temp,iv_watrc,      &
                 iv_moistc,d_canyontemp,d_canyonmix,d_sigd,d_tempc,d_mixrc,d_topu,d_roaddelta,d_vegdeltac,d_rdsndelta,d_netrad,  &
                 d_netemiss,d_tranc,d_evapc,d_cwa,d_cra,d_cwe,d_cww,d_crw,d_crr,d_cwr,d_accool,d_totdepth,d_c1c,                 &
                 d_acout,d_canyonrgout,sg_vegc,a_umag,a_rho,a_rg,a_rnd,a_snd,ddt,acond_road,acond_walle,acond_wallw,acond_vegc,  &
                 acond_rdsn,wallpsi,roadpsi,if_bldheight,if_hwratio,if_sigmavegc,if_roademiss,if_vegemissc,if_wallemiss,         &
                 if_bldtemp,if_zovegc,if_vegrlaic,if_vegrsminc,if_swilt,if_sfc,ip_vegtempc,ip_lzom,ip_lzoh,ip_cndzmin,           &
                 d_traf,ldratio)
end if
! ---------------------------------------------------------------    

! balance vegetation energy budget
fg_vegc=sg_vegc+rg_vegc-eg_vegc

! road snow energy balance
evct=sg_rdsn+rg_rdsn-fg_rdsn-eg_rdsn*ls/lv-gardsn

return
end subroutine solverdsn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solve for vegetation canopy temperature

subroutine solverdvg(cn,evct,rg_road,rg_walle,rg_wallw,rg_vegc,rg_rdsn,fg_road,fg_walle,fg_wallw,fg_vegc,fg_rdsn,eg_road,eg_vegc,  &
                     eg_rdsn,gardsn,rdsnmelt,topinvres,rdsntemp,ird_temp,ird_water,ird_snow,ird_den,iwe_temp,iww_temp,iv_watrc,    &
                     iv_moistc,d_canyontemp,d_canyonmix,d_sigd,d_tempc,d_mixrc,d_topu,d_roaddelta,d_vegdeltac,d_rdsndelta,         &
                     d_netrad,d_netemiss,d_tranc,d_evapc,d_cwa,d_cra,d_cwe,d_cww,d_crw,d_crr,d_cwr,d_accool,                       &
                     d_totdepth,d_c1c,d_acout,d_canyonrgout,sg_vegc,a_umag,a_rho,a_rg,a_rnd,a_snd,ddt,acond_road,acond_walle,      &
                     acond_wallw,acond_vegc,acond_rdsn,wallpsi,roadpsi,if_bldheight,if_hwratio,if_sigmavegc,if_roademiss,          &
                     if_vegemissc,if_wallemiss,if_bldtemp,if_zovegc,if_vegrlaic,if_vegrsminc,if_swilt,if_sfc,ip_vegtempc,ip_lzom,  &
                     ip_lzoh,ip_cndzmin,d_traf,ldratio)

implicit none

integer, intent(in) :: cn
integer k
real, intent(in) :: ddt
real, dimension(cn), intent(out) :: evct,gardsn,rdsnmelt,topinvres
real, dimension(cn), intent(in) :: rdsntemp,wallpsi,roadpsi
real, dimension(cn), intent(in) :: ldratio
real, dimension(cn) :: oldval,evctveg,evctx,newval,ctmax,ctmin
real, dimension(cn), intent(in) :: sg_vegc
real, dimension(cn), intent(inout) :: rg_road,rg_walle,rg_wallw,rg_vegc,rg_rdsn
real, dimension(cn), intent(inout) :: fg_road,fg_walle,fg_wallw,fg_vegc,fg_rdsn
real, dimension(cn), intent(inout) :: eg_road,eg_vegc,eg_rdsn
real, dimension(cn), intent(inout) :: acond_road,acond_walle,acond_wallw,acond_vegc,acond_rdsn
real, dimension(cn), intent(in) :: a_umag,a_rho,a_rg,a_rnd,a_snd
real, dimension(cn), intent(inout) :: d_canyontemp,d_canyonmix,d_sigd,d_tempc,d_mixrc,d_topu,d_roaddelta,d_vegdeltac,d_rdsndelta
real, dimension(cn), intent(inout) :: d_netrad,d_netemiss,d_tranc,d_evapc,d_cwa,d_cra,d_cwe,d_cww,d_crw,d_crr,d_cwr,d_accool
real, dimension(cn), intent(inout) :: d_totdepth,d_c1c,d_acout,d_traf,d_canyonrgout
real, dimension(cn), intent(in) :: ird_temp,ird_water,ird_snow,ird_den
real, dimension(cn), intent(in) :: iwe_temp,iww_temp
real, dimension(cn), intent(in) :: iv_watrc,iv_moistc
real, dimension(cn), intent(in) :: if_bldheight,if_hwratio,if_sigmavegc,if_roademiss,if_vegemissc,if_wallemiss,if_bldtemp
real, dimension(cn), intent(in) :: if_zovegc,if_vegrlaic,if_vegrsminc,if_swilt,if_sfc
real, dimension(cn), intent(inout) :: ip_vegtempc,ip_lzom,ip_lzoh,ip_cndzmin

! update canyon temperature
ctmax=max(d_tempc,ird_temp,iwe_temp,iww_temp,rdsntemp)+10. ! max temp
ctmin=min(d_tempc,ird_temp,iwe_temp,iww_temp,rdsntemp)-10. ! min temp
d_canyontemp=0.4*ctmax+0.6*ctmin
call solvecanyon(cn,evctveg,rg_road,rg_walle,rg_wallw,rg_vegc,rg_rdsn,fg_road,fg_walle,fg_wallw,fg_vegc,fg_rdsn,eg_road,eg_vegc, &
                 eg_rdsn,gardsn,rdsnmelt,topinvres,rdsntemp,ird_temp,ird_water,ird_snow,ird_den,iwe_temp,iww_temp,iv_watrc,      &
                 iv_moistc,d_canyontemp,d_canyonmix,d_sigd,d_tempc,d_mixrc,d_topu,d_roaddelta,d_vegdeltac,d_rdsndelta,d_netrad,  &
                 d_netemiss,d_tranc,d_evapc,d_cwa,d_cra,d_cwe,d_cww,d_crw,d_crr,d_cwr,d_accool,d_totdepth,d_c1c,                 &
                 d_acout,d_canyonrgout,sg_vegc,a_umag,a_rho,a_rg,a_rnd,a_snd,ddt,acond_road,acond_walle,acond_wallw,acond_vegc,  &
                 acond_rdsn,wallpsi,roadpsi,if_bldheight,if_hwratio,if_sigmavegc,if_roademiss,if_vegemissc,if_wallemiss,         &
                 if_bldtemp,if_zovegc,if_vegrlaic,if_vegrsminc,if_swilt,if_sfc,ip_vegtempc,ip_lzom,ip_lzoh,ip_cndzmin,   &
                 d_traf,ldratio)
oldval=d_canyontemp
d_canyontemp=0.6*ctmax+0.4*ctmin
do k=1,nfgits
  evctx=evctveg
  call solvecanyon(cn,evctveg,rg_road,rg_walle,rg_wallw,rg_vegc,rg_rdsn,fg_road,fg_walle,fg_wallw,fg_vegc,fg_rdsn,eg_road,eg_vegc, &
                   eg_rdsn,gardsn,rdsnmelt,topinvres,rdsntemp,ird_temp,ird_water,ird_snow,ird_den,iwe_temp,iww_temp,iv_watrc,      &
                   iv_moistc,d_canyontemp,d_canyonmix,d_sigd,d_tempc,d_mixrc,d_topu,d_roaddelta,d_vegdeltac,d_rdsndelta,d_netrad,  &
                   d_netemiss,d_tranc,d_evapc,d_cwa,d_cra,d_cwe,d_cww,d_crw,d_crr,d_cwr,d_accool,d_totdepth,d_c1c,                 &
                   d_acout,d_canyonrgout,sg_vegc,a_umag,a_rho,a_rg,a_rnd,a_snd,ddt,acond_road,acond_walle,acond_wallw,acond_vegc,  &
                   acond_rdsn,wallpsi,roadpsi,if_bldheight,if_hwratio,if_sigmavegc,if_roademiss,if_vegemissc,if_wallemiss,         &
                   if_bldtemp,if_zovegc,if_vegrlaic,if_vegrsminc,if_swilt,if_sfc,ip_vegtempc,ip_lzom,ip_lzoh,ip_cndzmin,   &
                   d_traf,ldratio)
  evctx=evctveg-evctx
  where (abs(evctx).gt.tol) 
    newval=d_canyontemp-alpha*evctveg*(d_canyontemp-oldval)/evctx
    oldval=d_canyontemp
    d_canyontemp=newval
  end where
  d_canyontemp=min(max(d_canyontemp,ctmin),ctmax)
end do

! vegetation energy budget
evct=sg_vegc+rg_vegc-fg_vegc-eg_vegc

return
end subroutine solverdvg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solve for canyon temperature

subroutine solvecanyon(cn,evct,rg_road,rg_walle,rg_wallw,rg_vegc,rg_rdsn,fg_road,fg_walle,fg_wallw,fg_vegc,fg_rdsn,eg_road,eg_vegc, &
                       eg_rdsn,gardsn,rdsnmelt,topinvres,rdsntemp,ird_temp,ird_water,ird_snow,ird_den,iwe_temp,iww_temp,            &
                       iv_watrc,iv_moistc,d_canyontemp,d_canyonmix,d_sigd,d_tempc,d_mixrc,d_topu,d_roaddelta,d_vegdeltac,           &
                       d_rdsndelta,d_netrad,d_netemiss,d_tranc,d_evapc,d_cwa,d_cra,d_cwe,d_cww,d_crw,d_crr,d_cwr,d_accool,          &
                       d_totdepth,d_c1c,d_acout,d_canyonrgout,sg_vegc,a_umag,a_rho,a_rg,a_rnd,a_snd,ddt,acond_road,                 &
                       acond_walle,acond_wallw,acond_vegc,acond_rdsn,wallpsi,roadpsi,if_bldheight,if_hwratio,if_sigmavegc,          &
                       if_roademiss,if_vegemissc,if_wallemiss,if_bldtemp,if_zovegc,if_vegrlaic,if_vegrsminc,if_swilt,if_sfc,        &
                       ip_vegtempc,ip_lzom,ip_lzoh,ip_cndzmin,d_traf,ldratio)

implicit none

integer, intent(in) :: cn
integer ic
real, intent(in) :: ddt
real, dimension(cn), intent(out) :: gardsn,rdsnmelt,topinvres
real, dimension(cn), intent(out) :: evct
real, dimension(cn), intent(in) :: rdsntemp,wallpsi,roadpsi,ldratio
real, dimension(cn) :: roadqsat,rdsnqsat,fgtop
real, dimension(cn) :: vegqsat,res,f1,f2,f3,f4,ff
real, dimension(cn) :: dumroaddelta,dumvegdelta,cu,cduv
real, dimension(cn) :: newcanyonmix,oldcanyonmix,evctmx,evctmxdif
real, dimension(cn) :: z_on_l,dts,dtt,effbldheight,effhwratio,effwalle,effwallw,effvegc,effroad,effrdsn
real, dimension(cn), intent(in) :: sg_vegc
real, dimension(cn), intent(inout) :: rg_road,rg_walle,rg_wallw,rg_vegc,rg_rdsn
real, dimension(cn), intent(inout) :: fg_road,fg_walle,fg_wallw,fg_vegc,fg_rdsn
real, dimension(cn), intent(inout) :: eg_road,eg_vegc,eg_rdsn
real, dimension(cn), intent(inout) :: acond_road,acond_walle,acond_wallw,acond_vegc,acond_rdsn
real, dimension(cn), intent(in) :: a_umag,a_rho,a_rg,a_rnd,a_snd
real, dimension(cn), intent(inout) :: d_canyontemp,d_canyonmix,d_sigd,d_tempc,d_mixrc,d_topu,d_roaddelta,d_vegdeltac,d_rdsndelta
real, dimension(cn), intent(inout) :: d_netrad,d_netemiss,d_tranc,d_evapc,d_cwa,d_cra,d_cwe,d_cww,d_crw,d_crr,d_cwr,d_accool
real, dimension(cn), intent(inout) :: d_totdepth,d_c1c,d_acout,d_traf,d_canyonrgout
real, dimension(cn), intent(in) :: ird_temp,ird_snow,ird_water,ird_den
real, dimension(cn), intent(in) :: iwe_temp,iww_temp
real, dimension(cn), intent(in) :: iv_watrc,iv_moistc
real, dimension(cn), intent(in) :: if_bldheight,if_hwratio,if_sigmavegc,if_roademiss,if_vegemissc,if_wallemiss,if_bldtemp
real, dimension(cn), intent(in) :: if_zovegc,if_vegrlaic,if_vegrsminc,if_swilt,if_sfc
real, dimension(cn), intent(inout) :: ip_vegtempc,ip_lzom,ip_lzoh,ip_cndzmin

effbldheight=max(if_bldheight-6.*if_zovegc,0.1)/if_bldheight  ! MJT suggestion for tall vegetation
effhwratio=if_hwratio*effbldheight                            ! MJT suggestion for tall vegetation

! transpiration terms (developed by Eva in CCAM sflux.f and CSIRO9)
where (if_zovegc.lt.0.5)
  ff=1.1*sg_vegc/(if_vegrlaic*150.)
elsewhere
  ff=1.1*sg_vegc/(if_vegrlaic*30.)
end where
f1=(1.+ff)/(ff+if_vegrsminc*if_vegrlaic/5000.)
f2=max(0.5*(if_sfc-if_swilt)/max(iv_moistc-if_swilt,1.E-9),1.)
f4=max(1.-0.0016*(298.-d_canyontemp)**2,0.05)

! estimate mixing ratio at canyon surfaces
call getqsat(cn,roadqsat,ird_temp,d_sigd)   ! evaluate using pressure at displacement height
call getqsat(cn,vegqsat,ip_vegtempc,d_sigd) ! evaluate using pressure at displacement height
call getqsat(cn,rdsnqsat,rdsntemp,d_sigd)

! Solver for canyonmix
d_canyonmix=d_mixrc
ip_lzoh=ip_lzom
dtt=d_tempc*(1.+0.61*d_mixrc)
do ic=1,negits
  ! solve for aerodynamical resistance between canyon and atmosphere (neglect molecular diffusion)
  dts=d_canyontemp*(1.+0.61*d_canyonmix)
  call getinvres(cn,topinvres,cduv,z_on_l,ip_lzoh,ip_lzom,ip_cndzmin,dts,dtt,a_umag,3)
  call gettopu(cn,d_topu,a_umag,z_on_l,if_bldheight,if_hwratio,cduv,ip_cndzmin)

  if (resmeth.eq.0) then
    acond_road=(11.8+4.2*sqrt(acond_road**2+cduv*a_umag**2))/(aircp*a_rho) ! From Rowley, et al (1930)
    acond_walle=acond_road
    acond_wallw=acond_road
    acond_rdsn=acond_road
    acond_vegc=acond_road
  else if (resmeth.eq.2) then
    cu=acond_road*d_topu
    where (cu.le.5.)
      acond_road=(6.15+4.18*cu)/(aircp*a_rho)
    elsewhere
      acond_road=(7.51*cu**0.78)/(aircp*a_rho)
    end where
    acond_walle=acond_road
    acond_wallw=acond_road
    acond_rdsn=acond_road
    acond_vegc=acond_road
  end if

  ! correction for dew
  where (roadqsat.lt.d_canyonmix)
    dumroaddelta=1.
  elsewhere
    dumroaddelta=d_roaddelta
  endwhere
  where (vegqsat.lt.d_canyonmix)
    dumvegdelta=1.
  elsewhere
    dumvegdelta=d_vegdeltac
  endwhere
  ! remaining transpiration terms
  f3=max(1.-.00025*(vegqsat-d_canyonmix)*d_sigd/0.622,0.05)
  res=max(30.,if_vegrsminc*f1*f2/(f3*f4))
  ! balance canyon latent heat budget
  d_canyonmix=(d_rdsndelta*rdsnqsat*acond_rdsn*d_topu                                                          &
         +(1.-d_rdsndelta)*((1.-f_sigmavegc)*dumroaddelta*roadqsat*acond_road*d_topu                           &
         +if_sigmavegc*vegqsat*(dumvegdelta*acond_vegc*d_topu                                                  &
         +(1.-dumvegdelta)/(1./(acond_vegc*d_topu)+res)))+d_mixrc*topinvres)/                                  &
         (d_rdsndelta*acond_rdsn*d_topu+(1.-d_rdsndelta)*((1.-f_sigmavegc)*dumroaddelta*acond_road*d_topu      &
         +if_sigmavegc*(dumvegdelta*acond_vegc*d_topu+(1.-dumvegdelta)/(1./(acond_vegc*d_topu)+res)))+topinvres)
end do

! solve for canyon sensible heat flux ----------------------------------
fg_walle=aircp*a_rho*(iwe_temp-d_canyontemp)*acond_walle*d_topu
fg_wallw=aircp*a_rho*(iww_temp-d_canyontemp)*acond_wallw*d_topu
fg_road=aircp*a_rho*(ird_temp-d_canyontemp)*acond_road*d_topu
fg_vegc=aircp*a_rho*(ip_vegtempc-d_canyontemp)*acond_vegc*d_topu
where (d_rdsndelta.gt.0.)
  fg_rdsn=aircp*a_rho*(rdsntemp-d_canyontemp)*acond_rdsn*d_topu
elsewhere
  fg_rdsn=0.
end where
fgtop=aircp*a_rho*(d_canyontemp-d_tempc)*topinvres
! ---------------------------------------------------------------   

! calculate longwave radiation
d_netrad=d_rdsndelta*snowemiss*rdsntemp**4+(1.-d_rdsndelta)*((1.-if_sigmavegc)*if_roademiss*ird_temp**4 &
                  +if_sigmavegc*if_vegemissc*ip_vegtempc**4)
effwalle=if_wallemiss*(a_rg*d_cwa+sbconst*iwe_temp**4*(if_wallemiss*d_cwe-1.) & 
                  +sbconst*iww_temp**4*if_wallemiss*d_cww+sbconst*d_netrad*d_cwr)
rg_walle=effwalle*effbldheight+sbconst*(d_netrad-if_wallemiss*iwe_temp**4)*(1.-effbldheight)
effwallw=if_wallemiss*(a_rg*d_cwa+sbconst*iww_temp**4*(if_wallemiss*d_cwe-1.) &
                  +sbconst*iwe_temp**4*if_wallemiss*d_cww+sbconst*d_netrad*d_cwr)
rg_wallw=effwallw*effbldheight+sbconst*(d_netrad-if_wallemiss*iww_temp**4)*(1.-effbldheight)
effroad=if_roademiss*(a_rg*d_cra+sbconst*(d_netrad*d_crr-ird_temp**4) &
                  +sbconst*if_wallemiss*(iwe_temp**4+iww_temp**4)*d_crw)
rg_road=effroad+sbconst*(if_wallemiss*(iwe_temp**4+iww_temp**4) &
                  -2.*if_roademiss*ird_temp**4)*if_hwratio*(1.-effbldheight)
effvegc=if_vegemissc*(a_rg*d_cra+sbconst*(d_netrad*d_crr-ip_vegtempc**4) &
                  +sbconst*if_wallemiss*(iwe_temp**4+iww_temp**4)*d_crw)
rg_vegc=effvegc+sbconst*(if_wallemiss*(iwe_temp**4+iww_temp**4) &
                  -2.*if_vegemissc*ip_vegtempc**4)*if_hwratio*(1.-effbldheight)
where (d_rdsndelta.gt.0.)
  effrdsn=snowemiss*(a_rg*d_cra+sbconst*(-rdsntemp**4+d_netrad*d_crr) &
                    +sbconst*if_wallemiss*(iwe_temp**4+iww_temp**4)*d_crw)
  rg_rdsn=effrdsn+sbconst*(if_wallemiss*(iwe_temp**4+iww_temp**4) &
                    -2.*snowemiss*rdsntemp**4)*if_hwratio*(1.-effbldheight)
elsewhere
  effrdsn=0.
  rg_rdsn=0.
end where
! outgoing longwave radiation
d_canyonrgout=a_rg-d_rdsndelta*effrdsn-(1.-d_rdsndelta)*((1.-if_sigmavegc)*effroad+if_sigmavegc*effvegc) &
                  -effhwratio*(effwalle+effwallw)

! estimate snow melt
rdsnmelt=d_rdsndelta*max(0.,rdsntemp-273.16)/(icecp*ird_den*lf*ddt)

! calculate transpiration and evaporation of in-canyon vegetation
d_tranc=lv*min(max((1.-dumvegdelta)*a_rho*(vegqsat-d_canyonmix)/(1./(acond_vegc*d_topu)+res),0.), &
               max((iv_moistc-if_swilt)*d_totdepth*waterden/(d_c1c*ddt),0.))
d_evapc=lv*min(dumvegdelta*a_rho*(vegqsat-d_canyonmix)*acond_vegc*d_topu,iv_watrc/ddt+a_rnd)
eg_vegc=d_evapc+d_tranc

! calculate canyon latent heat fluxes
eg_road=lv*min(a_rho*dumroaddelta*(roadqsat-d_canyonmix)*acond_road*d_topu &
             ,ird_water/ddt+a_rnd+(1.-if_sigmavegc)*rdsnmelt)
where (d_rdsndelta.gt.0.)
  eg_rdsn=lv*min(a_rho*d_rdsndelta*max(0.,rdsnqsat-d_canyonmix)*acond_rdsn*d_topu &
                ,ird_snow/ddt+a_snd-rdsnmelt)
  gardsn=(rdsntemp-ird_temp)/ldratio ! use road temperature to represent canyon bottom surface temperature
                                     ! (i.e., we have ommited soil under vegetation temperature)
elsewhere
  eg_rdsn=0.
  gardsn=0.
end where

! update heat pumped into canyon with COP
d_accool=d_acout*(1.+max(d_canyontemp-if_bldtemp,0.)/if_bldtemp)
!d_accool=d_acout ! test energy conservation

! balance sensible heat flux
evct=fgtop-(d_rdsndelta*fg_rdsn+(1.-d_rdsndelta)*((1.-if_sigmavegc)*fg_road &
     +if_sigmavegc*fg_vegc)+if_hwratio*(fg_walle+fg_wallw)+d_traf+d_accool)

return
end subroutine solvecanyon

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve for green roof temperature

subroutine solverfveg(rg_vegr,fg_vegr,eg_vegr,acond_vegr,a_rg,a_umag,a_rho,a_rnd,sg_vegr,d_tempr,d_mixrr,d_rfdzmin, &
                      d_tranr,d_evapr,d_c1r,d_sigr,d_vegdeltar,ddt)

implicit none

integer k
real, intent(in) :: ddt
real, dimension(ufull), intent(out) :: rg_vegr,fg_vegr,eg_vegr,acond_vegr
real, dimension(ufull), intent(in) :: a_rg,a_umag,a_rho,a_rnd,sg_vegr
real, dimension(ufull), intent(inout) :: d_tempr,d_mixrr,d_rfdzmin,d_tranr,d_evapr,d_c1r,d_sigr,d_vegdeltar
real, dimension(ufull) :: ctmax,ctmin,newval,oldval,evctx,evctveg

if (any(f_sigmavegr.gt.0.)) then
  ! update green roof temperature
  ctmax=d_tempr+10. ! max temp
  ctmin=d_tempr-10. ! min temp
  p_vegtempr=0.4*ctmax+0.6*ctmin
  call solvegreenroof(evctveg,rg_vegr,fg_vegr,eg_vegr,acond_vegr,a_rg,a_umag,a_rho,a_rnd,sg_vegr,d_tempr,d_mixrr,d_rfdzmin, &
                      d_tranr,d_evapr,d_c1r,d_sigr,d_vegdeltar,ddt)
  oldval=p_vegtempr
  p_vegtempr=0.6*ctmax+0.4*ctmin
  do k=1,nfgits
    evctx=evctveg
    call solvegreenroof(evctveg,rg_vegr,fg_vegr,eg_vegr,acond_vegr,a_rg,a_umag,a_rho,a_rnd,sg_vegr,d_tempr,d_mixrr,d_rfdzmin, &
                        d_tranr,d_evapr,d_c1r,d_sigr,d_vegdeltar,ddt)
    evctx=evctveg-evctx
    where (abs(evctx).gt.tol) 
      newval=p_vegtempr-alpha*evctveg*(p_vegtempr-oldval)/evctx
      oldval=p_vegtempr
      p_vegtempr=newval
    end where
    p_vegtempr=min(max(p_vegtempr,ctmin),ctmax)
  end do
else
  p_vegtempr=d_tempr
  rg_vegr=0.
  fg_vegr=0.
  eg_vegr=0.
  acond_vegr=0.
  d_tranr=0.
  d_evapr=0.
end if

return
end subroutine solverfveg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve for green roof fluxes

subroutine solvegreenroof(evct,rg_vegr,fg_vegr,eg_vegr,acond_vegr,a_rg,a_umag,a_rho,a_rnd,sg_vegr,d_tempr,d_mixrr,d_rfdzmin, &
                          d_tranr,d_evapr,d_c1r,d_sigr,d_vegdeltar,ddt)

implicit none

real, intent(in) :: ddt
real, dimension(ufull), intent(out) :: evct,rg_vegr,fg_vegr,eg_vegr,acond_vegr
real, dimension(ufull), intent(in) :: a_rg,a_umag,a_rho,a_rnd,sg_vegr
real, dimension(ufull), intent(inout) :: d_tempr,d_mixrr,d_rfdzmin,d_tranr,d_evapr,d_c1r,d_sigr,d_vegdeltar
real, dimension(ufull) :: a,n,xe,xw,dts,dtt,z_on_l,ff,f1,f2,f3,f4
real, dimension(ufull) :: vegqsat,dumvegdelta,res

call getqsat(ufull,vegqsat,p_vegtempr,d_sigr)
where (vegqsat.lt.d_mixrr)
  dumvegdelta=1.
elsewhere
  dumvegdelta=d_vegdeltar
endwhere

! transpiration terms (developed by Eva in CCAM sflux.f and CSIRO9)
where (f_zovegr.lt.0.5)
  ff=1.1*sg_vegr/(f_vegrlair*150.)
elsewhere
  ff=1.1*sg_vegr/(f_vegrlair*30.)
end where
f1=(1.+ff)/(ff+f_vegrsminr*f_vegrlair/5000.)
f2=max(0.5*(f_sfc-f_swilt)/max(v_moistr-f_swilt,1.E-9),1.)
f3=max(1.-.00025*(vegqsat-d_mixrr)*d_sigr/0.622,0.05)
f4=max(1.-0.0016*(298.-d_tempr)**2,0.05)
res=max(30.,f_vegrsminr*f1*f2/(f3*f4))

n=max(min((v_moistr-f_swilt)/(f_sfc-f_swilt),1.),0.) ! veg wetfac (see sflux.f or cable_canopy.f90)
xw=(1.-dumvegdelta)*n+dumvegdelta

! calculate green roof sensible and latent heat fluxes
rg_vegr=f_vegemissr*(a_rg-sbconst*p_vegtempr**4)
! a is a dummy variable for lzomvegr
a=log(d_rfdzmin/f_zovegr)
! xe is a dummy variable for lzohvegr
xe=2.3+a
xw=vegqsat*xw ! green roof surface mixing ratio
dts=p_vegtempr*(1.+0.61*xw)
dtt=d_tempr*(1.+0.61*d_mixrr)
! Assume zot=0.1*zom (i.e., Kanda et al 2007, small experiment)
! n is a dummy variable for cd
call getinvres(ufull,acond_vegr,n,z_on_l,xe,a,d_rfdzmin,dts,dtt,a_umag,1)
acond_vegr=acond_vegr/a_umag
fg_vegr=aircp*a_rho*(p_vegtempr-d_tempr)*acond_vegr*a_umag

! calculate transpiration and evaporation of in-canyon vegetation
d_tranr=lv*min(max((1.-dumvegdelta)*a_rho*(vegqsat-d_mixrr)/(1./(acond_vegr*a_umag)+res),0.), &
               max((v_moistr-f_swilt)*f_vegdepthr*waterden/(d_c1r*ddt),0.))
d_evapr=lv*min(dumvegdelta*a_rho*(vegqsat-d_mixrr)*acond_vegr*a_umag,v_watrr/ddt+a_rnd)
eg_vegr=d_evapr+d_tranr

! balance green roof energy budget
evct=sg_vegr+rg_vegr-fg_vegr-eg_vegr

return
end subroutine solvegreenroof

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define traffic flux weights during the diurnal cycle

subroutine gettraffic(cn,trafficout,if_ctime,if_trafficfg)

implicit none

integer, intent(in) :: cn
integer, dimension(cn) :: ip
real, dimension(cn), intent(out) :: trafficout
real, dimension(cn) :: rp
real, dimension(cn), intent(in) :: if_ctime,if_trafficfg
! traffic diurnal cycle weights approximated from Coutts et al (2007)
real, dimension(25), parameter :: trafficcycle = (/ 0.1, 0.1, 0.1, 0.1, 0.5, 1.0, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, &
                                                    1.5, 1.5, 1.5, 1.5, 1.5, 1.4, 1.2,  1., 0.8, 0.6, 0.4, 0.2, &
                                                    0.1 /) 

ip=int(24.*if_ctime)
rp=24.*if_ctime-real(ip)
where (ip.lt.1) ip=ip+24
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
real, dimension(ufull) :: theta1,wdir,h,w
real, dimension(ufull), intent(in) :: a_udir

! rotate wind direction so that all cases are between 0 and pi
! walls are fliped at the end of the subroutine to account for additional pi rotation
where (a_udir.ge.0.)
  wdir=a_udir
elsewhere
  wdir=a_udir+pi
endwhere

h=f_bldheight
w=f_bldheight/f_hwratio

theta1=asin(min(w/(3.*h),1.))
wsuma=0.
wsumb=0.
fsum=0.  ! floor

! integrate jet on road, venting side (A)
a=0.
b=max(0.,wdir-pi+theta1)
if (any(b.gt.a)) call integratewind(wsuma,wsumb,fsum,a,b,h,w,wdir,z0,0)

! integrate jet on wall, venting side
a=max(0.,wdir-pi+theta1)
b=max(0.,wdir-theta1)
if (any(b.gt.a)) call integratewind(wsuma,wsumb,fsum,a,b,h,w,wdir,z0,1)

! integrate jet on road, venting side (B)
a=max(0.,wdir-theta1)
b=wdir
if (any(b.gt.a)) call integratewind(wsuma,wsumb,fsum,a,b,h,w,wdir,z0,0)

! integrate jet on road, recirculation side (A)
a=wdir
b=min(pi,wdir+theta1)
if (any(b.gt.a)) call integratewind(wsumb,wsuma,fsum,a,b,h,w,wdir,z0,0)

! integrate jet on wall, recirculation side
a=min(pi,wdir+theta1)
b=min(pi,wdir+pi-theta1)
if (any(b.gt.a)) call integratewind(wsumb,wsuma,fsum,a,b,h,w,wdir,z0,1)

! integrate jet on road, recirculation side (B)
a=min(pi,wdir+pi-theta1)
b=pi
if (any(b.gt.a)) call integratewind(wsumb,wsuma,fsum,a,b,h,w,wdir,z0,0)

! Correct for rotation of winds at start of subroutine
! 0.5 to adjust for factor of 2 in gettopu
where (a_udir.ge.0.)
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
real, dimension(ufull) :: a,b,wsuma,wsumb,fsum
real, dimension(ufull) :: theta1,wdir,h,w
real, dimension(ufull) :: dufa,dura,duva,ntheta
real, dimension(ufull) :: dufb,durb,duvb
real, dimension(ufull), intent(in) :: a_udir

! rotate wind direction so that all cases are between 0 and pi
! walls are fliped at the end of the subroutine to account for additional pi rotation
where (a_udir.ge.0.)
  wdir=a_udir
elsewhere
  wdir=a_udir+pi
endwhere

h=f_bldheight
w=f_bldheight/f_hwratio

theta1=acos(min(w/(3.*h),1.))

call winda(dufa,dura,duva,h,w,z0) ! jet on road
call windb(dufb,durb,duvb,h,w,z0) ! jet on wall
ntheta=2. ! i.e., int_0^pi sin(theta) dtheta = 2.)
where (wdir.lt.theta1.or.wdir.gt.pi-theta1) ! jet on wall
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
where (a_udir.ge.0.)
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

zolog=log(max(h,z0+0.1)/z0)
cuven=log(max(0.1*h,z0+0.1)/z0)/log(max(h,z0+0.1)/z0)
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
subroutine gettopu(cn,d_topu,a_umag,z_on_l,if_bldheight,if_hwratio,ip_cduv,ip_cndzmin)
      
implicit none

integer, intent(in) :: cn
real, dimension(cn), intent(in) :: z_on_l
real, dimension(cn) :: z0_on_l,bldheight
real, dimension(cn) :: pm0,pm1,integralm
real, dimension(cn) :: ustar,neutral
real, dimension(cn), intent(inout) :: d_topu
real, dimension(cn), intent(in) :: a_umag
real, dimension(cn), intent(in) :: if_bldheight,if_hwratio
real, dimension(cn), intent(inout) :: ip_cduv,ip_cndzmin

bldheight=if_bldheight*(1.-refheight)
ustar=sqrt(ip_cduv)*a_umag

z0_on_l=min(bldheight,ip_cndzmin)*z_on_l/ip_cndzmin ! calculate at canyon top
z0_on_l=min(z0_on_l,10.)
neutral = log(ip_cndzmin/min(bldheight,ip_cndzmin))
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
where (bldheight.lt.ip_cndzmin)
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
real, dimension(ufull) :: n

n=min(max(moist/f_ssat,0.218),1.)
dc1=(1.78*n+0.253)/(2.96*n-0.581)
where (n.gt.0.218)
  dc1=(1.78*(2.96*n-0.581)+2.96*(1.78*n+0.253))/((2.96*n-0.581)**2*f_ssat)
elsewhere
  dc1=0.
end where

return
end subroutine getc1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate screen diagnostics

subroutine scrncalc(a_mixr,a_umag,a_temp,u_ts,d_tempc,d_mixrc,d_rdsndelta,d_roaddelta,d_vegdeltac,d_sigd,smixr,rdsntemp,zonet)
      
implicit none

integer ic
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
real, dimension(ufull), intent(in) :: d_tempc,d_mixrc,d_rdsndelta,d_roaddelta,d_vegdeltac,d_sigd
real, parameter :: z0  = 1.5
real, parameter :: z10 = 10.

select case(scrnmeth)
  case(0) ! estimate screen diagnostics (slab at displacement height approach)
    thetav=d_tempc*(1.+0.61*a_mixr)
    sthetav=u_ts*(1.+0.61*smixr)
    lna=p_lzoh-p_lzom
    call dyerhicks(ufull,integralh,z_on_l,cd,thetavstar,thetav,sthetav,a_umag,p_cndzmin,p_lzom,lna,4)
    ustar=sqrt(cd)*a_umag
    qstar=vkar*(a_mixr-smixr)/integralh  
    tstar=vkar*(a_temp-u_ts)/integralh
    
    z0_on_l  = z0*z_on_l/p_cndzmin
    z10_on_l = z10*z_on_l/p_cndzmin
    z0_on_l  = min(z0_on_l,10.)
    z10_on_l = min(z10_on_l,10.)
    neutral   = log(p_cndzmin/z0)
    neutral10 = log(p_cndzmin/z10)
    where (z_on_l.lt.0.)
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
    call dyerhicks(ufull,integralh,z_on_l,cd,thetavstar,thetav,sthetav,a_umag,p_cndzmin,p_lzom,lna,4)
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
      integralh   = neutral-(ph1-ph0)
      integralm   = neutral-(pm1-pm0)
      integralm10 = neutral10-(pm1-pm10)
    endwhere
    ttop = thetav-tts*integralh/vkar
    tetp = a_temp-tstar*integralh/vkar
    qtop = a_mixr-qstar*integralh/vkar
    utop = a_umag-ustar*integralm/vkar

    where (f_bldheight.le.z10) ! above canyon
      p_u10=max(a_umag-ustar*integralm10/vkar,0.)
    end where

    ! assume standard stability functions hold for urban canyon (needs more work)
    tsurf=d_rdsndelta*rdsntemp+(1.-d_rdsndelta)*((1.-f_sigmavegc)*rd_temp(:,1)+f_sigmavegc*p_vegtempc)
    n=max(min((v_moistc-f_swilt)/(f_sfc-f_swilt),1.),0.)
    wf=(1.-d_rdsndelta)*((1.-f_sigmavegc)*d_roaddelta+f_sigmavegc*((1.-d_vegdeltac)*n+d_vegdeltac))
    call getqsat(ufull,qsurf,tsurf,d_sigd)
    qsurf=qsurf*wf
    n=log(f_bldheight/zonet)
    
    thetav=ttop*(1.+0.61*qtop)
    sthetav=tsurf*(1.+0.61*qsurf)
    lna=2.3
    call dyerhicks(ufull,integralh,z_on_l,cd,thetavstar,thetav,sthetav,utop,f_bldheight,n,lna,1)
    ustar=sqrt(cd)*utop
    tstar=vkar*(tetp-tsurf)/integralh
    qstar=vkar*(qtop-qsurf)/integralh
    
    z0_on_l   = z0*z_on_l/f_bldheight
    z10_on_l  = max(z10,f_bldheight)*z_on_l/f_bldheight
    z0_on_l   = min(z0_on_l,10.)
    z10_on_l  = min(z10_on_l,10.)
    neutral   = log(f_bldheight/z0)
    neutral10 = log(f_bldheight/max(z10,f_bldheight))
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
      integralh   = neutral-(ph1-ph0)
      integralm   = neutral-(pm1-pm0)
      integralm10 = neutral10-(pm1-pm10)
    endwhere

    p_tscrn = tetp-tstar*integralh/vkar
    p_qscrn = qtop-qstar*integralh/vkar
    p_uscrn = max(utop-ustar*integralm/vkar,0.)
    where (f_bldheight.gt.z10) ! within canyon
      p_u10 = max(utop-ustar*integralm10/vkar,0.)
    end where

  case(2) ! calculate screen diagnostics from canyon only
    tsurf=d_rdsndelta*rdsntemp+(1.-d_rdsndelta)*((1.-f_sigmavegc)*rd_temp(:,1)+f_sigmavegc*p_vegtempc)
    n=max(min((v_moistc-f_swilt)/(f_sfc-f_swilt),1.),0.)
    wf=(1.-d_rdsndelta)*((1.-f_sigmavegc)*d_roaddelta+f_sigmavegc*((1.-d_vegdeltac)*n+d_vegdeltac))
    call getqsat(ufull,qsurf,tsurf,d_sigd)
    qsurf=qsurf*wf
    n=log(f_bldheight/zonet)

    thetav=d_tempc*(1.+0.61*a_mixr)
    sthetav=tsurf*(1.+0.61*qsurf)
    lna=2.3
    call dyerhicks(ufull,integralh,z_on_l,cd,thetavstar,thetav,sthetav,a_umag,p_cndzmin,n,lna,1)
    ustar=sqrt(cd)*a_umag
    qstar=vkar*(a_mixr-smixr)/integralh
    tstar=vkar*(a_temp-tsurf)/integralh
    
    z0_on_l  = z0*z_on_l/p_cndzmin
    z10_on_l = z10*z_on_l/p_cndzmin
    z0_on_l  = min(z0_on_l,10.)
    z10_on_l = min(z10_on_l,10.)
    neutral   = log(p_cndzmin/z0)
    neutral10 = log(p_cndzmin/z10)
    where (z_on_l.lt.0.)
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
