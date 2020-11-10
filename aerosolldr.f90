! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2020 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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

! This is a prognostic aerosol model for CCAM based on the LDR scheme used in Mk3.6 (Rotstayn and Lohmann 2002)
    
! Plan to replace diagnosed sea-salt with prognostic version based on CTM.  Eventually employ GLOMAP for modes.

module aerosolldr

implicit none

private
public aldrcalc,aldrinit,aldrend,aldrloademiss,aldrloaderod,cldrop,convscav
public xtg,xtgsav,xtosav,naero
public itracdu,ndust
public dustdd,dustwd,duste,dust_burden
public itracbc,bce,bcdd,bcwd,bc_burden
public itracoc,oce,ocdd,ocwd,oc_burden
public itracdms,itracso2,itracso4
public dmse,dmsso2o,so2e,so2so4o,so2dd,so2wd,so4e,so4dd,so4wd
public dms_burden,so2_burden,so4_burden
public itracsa,nsalt,salte,saltdd,saltwd,salt_burden
public Ch_dust,zvolcemi,ticeu,aeroindir,so4mtn,carbmtn,saltsmallmtn,saltlargemtn
public dustden,dustreff,saltden,saltreff
public xtg_solub,zoxidant_g,erod,ndcls,emissfield,vso2

integer, save :: jk2,jk3,jk4,jk5,jk6,jk8,jk9                ! levels for injection
real, dimension(:,:,:), allocatable, save :: xtg            ! prognostic aerosols (see indexing below)
real, dimension(:,:,:), allocatable, save :: xtgsav         ! save for mass conservation in semi-Lagrangian models
real, dimension(:,:,:), allocatable, save :: xtosav         ! aerosol mixing ratio outside convective cloud
real, dimension(:,:,:), allocatable, save :: xtg_solub      ! aerosol mixing ratio that is dissolved in rain
real, dimension(:,:), allocatable, save :: erod             ! sand, clay and silt fraction that can erode
real, dimension(:,:), allocatable, save :: emissfield       ! non-volcanic emissions
real, dimension(:,:,:), allocatable, save :: zoxidant_g     ! oxidant fields
real, dimension(:), allocatable, save :: vso2               ! volcanic emissions
real, dimension(:,:), allocatable, save :: duste            ! Diagnostic - dust emissions
real, dimension(:,:), allocatable, save :: dustdd           ! Diagnostic - dust dry deposition
real, dimension(:,:), allocatable, save :: dustwd           ! Diagnostic - dust wet deposition
real, dimension(:,:), allocatable, save :: dust_burden      ! Diagnostic - dust burden
real, dimension(:), allocatable, save :: bce                ! Diagnostic - black carbon emissions
real, dimension(:), allocatable, save :: bcdd               ! Diagnostic - black carbon dry deposition
real, dimension(:), allocatable, save :: bcwd               ! Diagnostic - black carbon wet deposition
real, dimension(:), allocatable, save :: bc_burden          ! Diagnostic - black carbon burden
real, dimension(:), allocatable, save :: oce                ! Diagnostic - organic carbon emissions
real, dimension(:), allocatable, save :: ocdd               ! Diagnostic - organic carbon dry deposition
real, dimension(:), allocatable, save :: ocwd               ! Diagnostic - organic carbon wet deposition
real, dimension(:), allocatable, save :: oc_burden          ! Diagnostic - organic carbon burden
real, dimension(:), allocatable, save :: dmse               ! Diagnostic - DMS emissions
real, dimension(:), allocatable, save :: dmsso2o            ! Diagnostic - DMS->so2 oxidation
real, dimension(:), allocatable, save :: so2e               ! Diagnostic - so2 emissions
real, dimension(:), allocatable, save :: so2so4o            ! Diagnostic - so2->so4 oxidation
real, dimension(:), allocatable, save :: so2dd              ! Diagnostic - so2 dry deposition
real, dimension(:), allocatable, save :: so2wd              ! Diagnostic - so2 wet deposition
real, dimension(:), allocatable, save :: so4e               ! Diagnostic - so4 emissions
real, dimension(:), allocatable, save :: so4dd              ! Diagnostic - so4 dry deposition
real, dimension(:), allocatable, save :: so4wd              ! Diagnostic - so4 wet deposition
real, dimension(:), allocatable, save :: salte              ! Diagnostic - salt emission
real, dimension(:), allocatable, save :: saltdd             ! Diagnostic - salt dry deposition
real, dimension(:), allocatable, save :: saltwd             ! Diagnostic - salt wet deposition
real, dimension(:), allocatable, save :: dms_burden         ! Diagnostic - DMS burden
real, dimension(:), allocatable, save :: so2_burden         ! Diagnostic - so2 burden
real, dimension(:), allocatable, save :: so4_burden         ! Diagnostic - so4 burden
real, dimension(:), allocatable, save :: salt_burden        ! Diagnostic - salt burden
!$acc declare create(jk2,jk3,jk4,jk5,jk6,jk8,jk9)

! tracers
integer, parameter :: nsulf = 3
integer, parameter :: ncarb = 4
integer, parameter :: ndust = 4
integer, parameter :: nsalt = 2
integer, parameter :: naero = nsulf+ncarb+ndust+nsalt ! Tracers: DMS, SO2, SO4, BCO, BCI, OCO, OCI, DUST(4), SALT(2)
integer, parameter :: itracdms = 1                  ! Index for DMS tracer
integer, parameter :: itracso2 = 2                  ! Index for SO2 tracer
integer, parameter :: itracso4 = 3                  ! Index for SO4 tracer
integer, parameter :: itracbc = nsulf+1             ! Index for BC tracer (hydrophobic, hydrophillic)
integer, parameter :: itracoc = nsulf+3             ! Index for OC tracer (hydrophobic, hydrophillic)
integer, parameter :: itracdu = nsulf+ncarb+1       ! Index for dust tracer
integer, parameter :: itracsa = nsulf+ncarb+ndust+1 ! Index for salt tracer
integer, parameter :: ndcls = 3                     ! Number of dust emission classes (sand, silt, clay)

! emission indices
integer, parameter :: iso2a1 = 1      ! SO2/SO4 Anthropogenic surface
integer, parameter :: iso2a2 = 2      ! SO2/SO4 Anthropogenic upper level
integer, parameter :: ibca1  = 3      ! BC Anthropogenic surface
integer, parameter :: ibca2  = 4      ! BC Anthropogenic upper level
integer, parameter :: ioca1  = 5      ! OC Anthropogenic surface
integer, parameter :: ioca2  = 6      ! OC Anthropogenic upper level
integer, parameter :: iso2b1 = 7      ! SO2/SO4 BiomassBurning surface
integer, parameter :: iso2b2 = 8      ! SO2/SO4 BiomassBurning upper level
integer, parameter :: ibcb1  = 9      ! BC BiomassBurning surface
integer, parameter :: ibcb2  = 10     ! BC BiomassBurning upper level
integer, parameter :: iocb1  = 11     ! OC BiomassBurning surface
integer, parameter :: iocb2  = 12     ! OC BiomassBurning upper level
integer, parameter :: idmso  = 13     ! DMS ocean
integer, parameter :: idmst  = 14     ! DMS land
integer, parameter :: iocna  = 15     ! Natural organic

! options
integer, save :: enhanceu10 = 0                 ! Modify 10m wind speed for emissions (0=none, 1=quadrature, 2=linear)
integer, save :: aeroindir  = 0                 ! Indirect effect (0=SO4+Carbon+salt, 1=SO4, 2=None)
real, parameter :: zmin     = 1.e-20            ! Minimum concentration tolerance
!$acc declare create(enhanceu10,aeroindir)

! physical constants
real, parameter :: grav      = 9.80616          ! Gravitation constant
real, parameter :: rdry      = 287.04           ! Specific gas const for dry air
real, parameter :: cp        = 1004.64          ! Heat capacity of air
real, parameter :: hl        = 2.5104e6         ! Latent heat of vaporisation
real, parameter :: vkar      = 0.4              ! von Karman constant
real, parameter :: rhos      = 100.             ! Assumed density of snow in kg/m^3

! emission and deposition constants
real, save :: zvolcemi       = 8.               ! Total emission from volcanoes (TgS/yr)
real, save :: Ch_dust        = 1.e-9            ! Transfer coeff for type natural source (kg*s2/m5)
!$acc declare create(zvolcemi,ch_dust)

! Indirect effect coefficients
! converts from mass (kg/m3) to number concentration (1/m3) for dist'n
real, save :: so4mtn = 1.24e17                  ! Penner et al (1998)
real, save :: carbmtn = 1.25e17                 ! Penner et al (1998)
real, save :: saltsmallmtn = 1.89e16            ! Nillson et al. number mode radius = 0.1 um, sd=2
real, save :: saltlargemtn = 1.1e14             ! Nillson et al. number mode radius = 0.5 um, sd=2
!real, save :: so4mtn = 1.69e17                 ! IPCC (2001) Table 5.1
!real, save :: carbmtn = 1.21e17                ! IPCC (2001) Table 5.1
!real, save :: carbmtn = 2.30e17                ! counts Aitken mode as well as accumulation mode carb aerosols
!real, save :: saltsmallmtn = 3.79e17           ! Herzog number mode radius = 0.035 um, sd=1.92, rho=2.165 g/cm3
!real, save :: saltlargemtn = 7.25e14           ! Herzog number mode radius = 0.35 um, sd=1.7, rho=2.165
!$acc declare create(so4mtn,carbmtn,saltsmallmtn,saltlargemtn)

! Dust coefficients
real, dimension(ndust), parameter :: dustden = (/ 2500., 2650., 2650., 2650. /)    ! Density of dust (kg/m3)
                                                                                   ! (Clay, small silt, small slit, small silt)
real, dimension(ndust), parameter :: dustreff = (/ 0.73e-6,1.4e-6,2.4e-6,4.5e-6 /) ! Main effective radius (m)
                                                                                   ! (Clay, small silt, small slit, small silt)

! Salt coefficients
real, dimension(2), parameter :: saltden = (/ 2165., 2165. /)     ! density of salt
real, dimension(2), parameter :: saltreff = (/ 0.1e-6, 0.5e-6 /)  ! radius of salt (um)

! convective scavenging coefficients
real, parameter :: ticeu     = 263.16           ! Temperature for freezing in convective updraft

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialisation

subroutine aldrinit(ifull,iextra,kl,sig)

implicit none

integer, intent(in) :: ifull,iextra,kl
integer pos(1)
real, dimension(kl), intent(in) :: sig

allocate(xtg(ifull+iextra,kl,naero),xtgsav(ifull,kl,naero))
allocate(xtosav(ifull,kl,naero),vso2(ifull))
allocate(emissfield(ifull,15))
allocate(zoxidant_g(ifull,kl,4),erod(ifull,ndcls))
allocate(duste(ifull,ndust),dustdd(ifull,ndust),dustwd(ifull,ndust),dust_burden(ifull,ndust))
allocate(bce(ifull),bcdd(ifull),bcwd(ifull))
allocate(bc_burden(ifull))
allocate(oce(ifull),ocdd(ifull),ocwd(ifull))
allocate(oc_burden(ifull))
allocate(dmse(ifull),dmsso2o(ifull))
allocate(so2e(ifull),so2so4o(ifull),so2dd(ifull),so2wd(ifull))
allocate(so4e(ifull),so4dd(ifull),so4wd(ifull))
allocate(dms_burden(ifull),so2_burden(ifull),so4_burden(ifull))
allocate(salte(ifull),saltdd(ifull),saltwd(ifull),salt_burden(ifull))

xtg=0.
xtgsav=0.
xtosav=0.
vso2=0.
emissfield=0.
zoxidant_g=0.
erod=0.
duste=0.
dustdd=0.
dustwd=0.
dust_burden=0.
bce=0.
bcdd=0.
bcwd=0.
bc_burden=0.
oce=0.
ocdd=0.
ocwd=0.
oc_burden=0.
dmse=0.
dmsso2o=0.
so2e=0.
so2so4o=0.
so2dd=0.
so2wd=0.
so4e=0.
so4dd=0.
so4wd=0.
salte=0.
saltdd=0.
saltwd=0.
dms_burden=0.
so2_burden=0.
so4_burden=0.
salt_burden=0.

! MJT - define injection levels

! scale to CSIRO9 18 levels
! jk2 is top of level=1, bottom of level=2
!pos=maxloc(sig,sig<=0.993) ! 65m
!jk2=max(pos(1),2)
jk2=2
! jk3 is top of level=2, bottom of level=3
pos=maxloc(sig,sig<=0.967) ! 300m
jk3=max(pos(1),jk2+1)
! jk4 is top of level=3, bottom of level=4
pos=maxloc(sig,sig<=0.930) ! 600m
jk4=max(pos(1),jk3+1)
! jk5 is top of level=4, bottom of level=5
pos=maxloc(sig,sig<=0.882) ! 1,000m
jk5=max(pos(1),jk4+1)
! jk6 is top of level=5, bottom of level=6
pos=maxloc(sig,sig<=0.800) ! 1,800m
jk6=max(pos(1),jk5+1)
! jk8 is top of level=7, bottom of level=8
pos=maxloc(sig,sig<=0.530) ! 5,000m
jk8=max(pos(1),jk6+1)
! jk9 is top of level=8, bottom of level=9
pos=maxloc(sig,sig<=0.350) ! 8,000m
jk9=max(pos(1),jk8+1)

!$acc update device(jk2,jk3,jk4,jk5,jk6,jk8,jk9)
!$acc update device(enhanceu10,aeroindir)
!$acc update device(zvolcemi,ch_dust)
!$acc update device(so4mtn,carbmtn,saltsmallmtn,saltlargemtn)

return
end subroutine aldrinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! End

subroutine aldrend

implicit none

deallocate(xtg,xtgsav,xtosav)
deallocate(vso2)
deallocate(emissfield)
deallocate(zoxidant_g,erod)
deallocate(duste,dustdd,dustwd,dust_burden)
deallocate(bce,bcdd,bcwd)
deallocate(bc_burden)
deallocate(oce,ocdd,ocwd)
deallocate(oc_burden)
deallocate(dmse,dmsso2o)
deallocate(so2e,so2so4o,so2dd,so2wd)
deallocate(so4e,so4dd,so4wd)
deallocate(dms_burden,so2_burden,so4_burden)
deallocate(salte,saltdd,saltwd,salt_burden)

return
end subroutine aldrend

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Load emission arrays

subroutine aldrloademiss(index,aa)

implicit none

integer, intent(in) :: index
real, dimension(:), intent(in) :: aa

if ( index<16 ) then
  emissfield(:,index)=aa(1:size(emissfield,1)) ! Then follow SO2, BC and OC from anthro (a) and biomass-burning (b) levels 1 and 2
elseif ( index==16 ) then
  vso2(:)=aa(1:size(vso2))             ! volcanic
else
  write(6,*) "ERROR: index out-of-range for aldrloademiss"
  stop
end if

! already rescaled in aeroemiss
!if (index==iso2b1.or.index==iso2b2.or.index==iso2a1.or.index==iso2a2) then
!  ! convert SO2 emissions to kgS/m2/s
!  emissfield(:,index)=0.5*emissfield(:,index)
!end if

!if (index==iocna) then
! Default yield for natural organics is about 13%, or 16.4 TgC p.a., which may be
! a gross underestimate e.g., Tsigaridis & Kanakidou, Atmos. Chem. Phys. (2003) and
! Kanakidou et al., Atmos. Chem. Phys. (2005). Try 35 TgC for now.
!  emissfield(:,index)=(35./16.4)*emissfield(:,index)
!end if

return
end subroutine aldrloademiss

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Load soil data

subroutine aldrloaderod(inda,aa)

implicit none

integer, intent(in) :: inda
real, dimension(:), intent(in) :: aa

! EROD is the soil fraction of Sand (inda=1), Silt (inda=2) and Clay (inda=3) that can erode
erod(:,inda) = aa(1:size(erod,1))

return
end subroutine aldrloaderod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main routine

subroutine aldrcalc(dt,sig,dz,wg,pblh,prf,ts,ttg,condc,snowd,taudar,fg,eg,v10m,                    &
                    ustar,zo,land,fracice,tsigmf,qvg,qlg,qfg,stratcloud,clcon,cldcon,pccw,rhoa,vt, &
                    pfprec,pfmelt,pfsnow,pfsubl,plambs,pmrate,pmaccr,pfstayice,                    &
                    pqfsedice,prscav,prfreeze,zdayfac,kbsav,xtg,duste,dustdd,xtosav,               &
                    dmsso2o,so2so4o,dust_burden,bc_burden,oc_burden,dms_burden,                    &
                    so2_burden,so4_burden,erod,zoxidant,so2wd,so4wd,bcwd,ocwd,dustwd,              &
                    emissfield,vso2,dmse,so2e,so4e,bce,oce,so2dd,so4dd,bcdd,ocdd,salte,saltdd,     &
                    saltwd,salt_burden,dustden,dustreff,saltden,saltreff,locean,imax,kl)
!$acc routine vector

implicit none

integer, intent(in) :: imax, kl
integer, dimension(imax), intent(in) :: kbsav  ! Bottom of convective cloud
real, intent(in) :: dt                         ! Time step
real, dimension(kl), intent(in) :: sig         ! Sigma levels
real, dimension(imax), intent(in) :: wg        ! Soil moisture fraction of field capacity
real, dimension(imax), intent(in) :: prf       ! Surface pressure
real, dimension(imax), intent(in) :: ts        ! Surface temperture
real, dimension(imax), intent(in) :: pblh      ! Boundary layer height
real, dimension(imax), intent(in) :: v10m      ! 10m wind speed
real, dimension(imax), intent(in) :: condc     ! Convective rainfall
real, dimension(imax), intent(in) :: snowd     ! Snow depth
real, dimension(imax), intent(in) :: taudar    ! Fraction of time sunlit
real, dimension(imax), intent(in) :: fg        ! Sensible heat flux
real, dimension(imax), intent(in) :: eg        ! Latent heat flux
real, dimension(imax), intent(in) :: ustar     ! Friction velocity
real, dimension(imax), intent(in) :: zo        ! Roughness length
real, dimension(imax), intent(in) :: fracice   ! Sea-ice fraction
real, dimension(imax), intent(in) :: tsigmf    ! Vegetation fraction
real, dimension(imax), intent(in) :: vt        ! transfer velocity
real, dimension(imax), intent(in) :: zdayfac   ! scale factor for day length
real, dimension(imax,kl), intent(in) :: dz     ! thickness of vertical levels (m)
real, dimension(imax,kl), intent(in) :: ttg    ! Air temperature
real, dimension(imax,kl), intent(in) :: qvg    ! liquid water mixing ratio
real, dimension(imax,kl), intent(in) :: qlg    ! liquid water mixing ratio
real, dimension(imax,kl), intent(in) :: qfg    ! frozen water mixing ratio
real, dimension(imax,kl), intent(in) :: stratcloud  ! stratiform cloud fraction
real, dimension(imax,kl), intent(in) :: clcon  ! convective cloud fraction
real, dimension(imax), intent(in) :: cldcon    ! Convective rainfall area fraction
real, dimension(imax,kl), intent(in) :: pccw
real, dimension(imax,kl), intent(in) :: rhoa   ! density of air (kg/m3)
real, dimension(imax,kl), intent(in) :: pfprec, pfmelt, pfsnow         ! from LDR prog cloud
real, dimension(imax,kl), intent(in) :: pfsubl, plambs, pmrate         ! from LDR prog cloud
real, dimension(imax,kl), intent(in) :: pmaccr, pqfsedice, prscav      ! from LDR prog cloud
real, dimension(imax,kl), intent(in) :: prfreeze                       ! from LDR prog cloud
real, dimension(imax,kl), intent(in) :: pfstayice                      ! from LDR prog cloud
logical, dimension(imax), intent(in) :: land   ! land/water mask (t=land).  Water includes lakes and ocean
logical, dimension(imax), intent(in) :: locean ! sea mask without lakes (t=ocean)
real, dimension(imax,kl,naero), intent(inout) :: xtg
real, dimension(imax,ndust), intent(inout) :: duste
real, dimension(imax,ndust), intent(inout) :: dustdd
real, dimension(imax,kl,naero), intent(in) :: xtosav
real, dimension(imax), intent(inout) :: dmsso2o
real, dimension(imax), intent(inout) :: so2so4o
real, dimension(imax,ndust), intent(inout) :: dust_burden
real, dimension(imax), intent(inout) :: bc_burden
real, dimension(imax), intent(inout) :: oc_burden
real, dimension(imax), intent(inout) :: dms_burden
real, dimension(imax), intent(inout) :: so2_burden
real, dimension(imax), intent(inout) :: so4_burden
real, dimension(imax,ndcls), intent(in) :: erod
real, dimension(imax,kl,4), intent(in) :: zoxidant
real, dimension(imax), intent(inout) :: so2wd
real, dimension(imax), intent(inout) :: so4wd
real, dimension(imax), intent(inout) :: bcwd
real, dimension(imax), intent(inout) :: ocwd
real, dimension(imax,ndust), intent(inout) :: dustwd
real, dimension(imax,15), intent(in) :: emissfield
real, dimension(imax), intent(in) :: vso2
real, dimension(imax), intent(inout) :: dmse
real, dimension(imax), intent(inout) :: so2e
real, dimension(imax), intent(inout) :: so4e
real, dimension(imax), intent(inout) :: bce
real, dimension(imax), intent(inout) :: oce
real, dimension(imax), intent(inout) :: so2dd
real, dimension(imax), intent(inout) :: so4dd
real, dimension(imax), intent(inout) :: bcdd
real, dimension(imax), intent(inout) :: ocdd
real, dimension(imax), intent(inout) :: salte
real, dimension(imax), intent(inout) :: saltdd
real, dimension(imax), intent(inout) :: saltwd
real, dimension(imax), intent(inout) :: salt_burden
real, dimension(ndust), intent(in) :: dustden, dustreff
real, dimension(2), intent(in) :: saltden, saltreff
real, dimension(imax,naero) :: xtem
real, dimension(imax,kl,naero) :: xte,xtu,xtm1
real, dimension(imax,kl) :: aphp1
real, dimension(imax,kl) :: pclcon
real, dimension(imax,kl) :: prhop1,ptp1
real, dimension(imax,kl) :: pclcover,pcfcover,pmlwc,pmiwc,pfconv
real, dimension(imax) :: bbem,fracc
real, dimension(imax) :: so2oh,so2h2,so2o3,dmsoh,dmsn3
real, dimension(imax) :: cgssnowd
real, dimension(imax) :: veff,vefn
real, dimension(imax) :: qtot
real, dimension(imax) :: rrate,Wstar3,Vgust_free,Vgust_deep
real, dimension(imax) :: v10n,thetav,burden
real, dimension(imax,ndust) :: dcola,dcolb
real, dimension(imax,ndust) :: oldduste
real, dimension(imax) :: oldsalte
real, parameter :: beta = 0.65
integer nt,k

#ifdef debugaero
if ( maxval(xtg(1:imax,:,:))>6.5e-6 ) then
  write(6,*) "xtg out-of-range at start of aldrcalc"
  write(6,*) "xtg maxval,maxloc ",maxval(xtg(1:imax,:,:)),maxloc(xtg(1:imax,:,:))
end if
#endif

cgssnowd(:) = 1.E-3*snowd

! Calculate sub-grid Vgust
v10n(:) = ustar*log(10./zo)/vkar ! neutral wind speed
! Mesoscale enhancement follows Redelsperger et al. (2000), J. Climate 13, 402-421.
! Equation numbers follow Fairall et al. 1996, JGR 101, 3747-3764.

! Calculate convective scaling velocity (Eq.17) and gustiness velocity (Eq.16)
thetav = ttg(1:imax,1)*(1.+0.61*qvg(1:imax,1))
Wstar3 = max(0.,(grav*pblh/thetav)*(fg/cp+0.61*ttg(1:imax,1)*eg/hl)/rhoa(:,1))
Vgust_free = beta*Wstar3**(1./3.)
! Calculate the Redelsperger-based Vgust_deep if deep convection is present.
! Note that Redelspreger gives two other parameterizations, based on
! the updraft or downdraught mass fluxes respectively.
rrate = 8640.*condc/dt   !Rainfall rate in cm/day
Vgust_deep = (19.8*rrate*rrate/(1.5+rrate+rrate*rrate))**0.4
! Calculate effective 10m wind (Eq. 15)
! These can plausibly be added in quadrature, or linearly, the latter giving a much larger
! effect (Lunt & Valdes, JGR, 2002).
select case(enhanceu10)
  case(0)
    veff = v10m
    vefn = v10n
  case(1)
    veff = sqrt( v10m*v10m + Vgust_free*Vgust_free + Vgust_deep*Vgust_deep )
    vefn = sqrt( v10n*v10n + Vgust_free*Vgust_free + Vgust_deep*Vgust_deep )
  case(2)
    veff = v10m + Vgust_free + Vgust_deep
    vefn = v10n + Vgust_free + Vgust_deep
end select

! Emission and dry deposition (sulfur cycle and carbonaceous aerosols)
call xtemiss(dt, rhoa, ts, fracice, vefn, land, tsigmf, cgssnowd, wg, dz,      & !Inputs
             xte, xtem, bbem,                                                  & !Outputs
             emissfield,vso2,dmse,so2e,so4e,bce,oce,xtg,so2dd,so4dd,bcdd,ocdd, &
             imax,kl)                                                            !Inputs
do nt = 1,naero
  do k = 1,kl
    xtg(:,k,nt) = max( xtg(:,k,nt)+xte(:,k,nt)*dt, 0. )
  end do
end do

#ifdef debugaero
if ( maxval(xtg(1:imax,:,:))>6.5e-6 ) then
  write(6,*) "xtg out-of-range after xtemiss"
  write(6,*) "xtg maxval,maxloc ",maxval(xtg(1:imax,:,:)),maxloc(xtg(1:imax,:,:))
end if
#endif

! Emission and dry deposition of dust
do k = 1,kl
  ! calculate air pressure
  aphp1(:,k) = prf(:)*sig(k)*0.01 ! hPa
end do
! Calculate integrated column dust loading before settling and deposition
do k = 1,ndust
  oldduste(:,k) = duste(:,k) ! duste is cumulative dust emissions
  dcola(:,k) = sum( rhoa(:,:)*xtg(1:imax,:,itracdu+k-1)*dz(:,:), dim=2 )
end do  
! Calculate the settling of large dust particles
call dsettling(dt,rhoa,ttg,dz,aphp1,xtg,dustden,dustreff,imax,kl)

#ifdef debugaero
if ( maxval(xtg(1:imax,:,:))>6.5e-6 ) then
  write(6,*) "xtg out-of-range after dsettling"
  write(6,*) "xtg maxval,maxloc ",maxval(xtg(1:imax,:,:)),maxloc(xtg(1:imax,:,:))
end if
#endif

! Calculate dust emission and turbulent dry deposition at the surface
call dustem(dt,rhoa(:,1),wg,veff,dz(:,1),vt,snowd,erod,duste,xtg, &
            dustden,dustreff,imax,kl)
do k = 1,ndust
  ! Calculate integrated column dust after settling
  dcolb(:,k) = sum( rhoa(:,:)*xtg(1:imax,:,itracdu+k-1)*dz(:,:), dim=2 )
  ! Calculate deposition flux to surface
  dustdd(:,k) = dustdd(:,k) + (dcola(:,k)-dcolb(:,k))/dt + duste(:,k) - oldduste(:,k)  
end do  

#ifdef debugaero
if ( maxval(xtg(1:imax,:,:))>6.5e-6 ) then
  write(6,*) "xtg out-of-range after dustem"
  write(6,*) "xtg maxval,maxloc ",maxval(xtg(1:imax,:,:)),maxloc(xtg(1:imax,:,:))
end if
#endif

! Decay of hydrophobic black and organic carbon into hydrophilic forms
call xtsink(dt,xtg,imax,kl)

#ifdef debugaero
if ( maxval(xtg(1:imax,:,:))>6.5e-6 ) then
  write(6,*) "xtg out-of-range after xtsink"
  write(6,*) "xtg maxval,maxloc ",maxval(xtg(1:imax,:,:)),maxloc(xtg(1:imax,:,:))
end if
#endif

oldsalte = salte ! salte is cumulative salt emissions
dcola(:,1) = 0.
do k = 1,nsalt
  dcola(:,1) = dcola(:,1) + sum( rhoa(:,:)*xtg(1:imax,:,itracsa+k-1)*dz(:,:), dim=2 )
end do
  
! Calculate the settling of large salt particles
call ssettling(dt,rhoa,ttg,dz,aphp1,xtg,saltden,saltreff,imax,kl)

#ifdef debugaero
if ( maxval(xtg(1:imax,:,:))>6.5e-6 ) then
  write(6,*) "xtg out-of-range after ssettling"
  write(6,*) "xtg maxval,maxloc ",maxval(xtg(1:imax,:,:)),maxloc(xtg(1:imax,:,:))
end if
#endif

! Calculate salt emission and turbulent dry deposition at the surface
call seasaltem(dt,veff,vt,rhoa(:,1),dz(:,1),salte,xtg,saltreff,locean,imax,kl)

#ifdef debugaero
if ( maxval(xtg(1:imax,:,:))>6.5e-6 ) then
  write(6,*) "xtg out-of-range after seasaltem"
  write(6,*) "xtg maxval,maxloc ",maxval(xtg(1:imax,:,:)),maxloc(xtg(1:imax,:,:))
end if
#endif

dcolb(:,1) = 0.
do k = 1,nsalt
  ! Calculate integrated column dust after settling
  dcolb(:,1) = dcolb(:,1) + sum( rhoa(:,:)*xtg(1:imax,:,itracsa+k-1)*dz(:,:), dim=2 )
end do  
! Calculate deposition flux to surface
saltdd = saltdd + (dcola(:,1)-dcolb(:,1))/dt + salte - oldsalte  

#ifdef debugaero
if ( maxval(xtg(1:imax,:,:))>6.5e-6 ) then
  write(6,*) "xtg out-of-range after seasalt"
  write(6,*) "xtg maxval,maxloc ",maxval(xtg(1:imax,:,:)),maxloc(xtg(1:imax,:,:))
end if
#endif

! Aerosol chemistry and wet deposition
! Need to invert vertical levels for ECHAM code... Don't you hate that?
do nt = 1,naero
  do k = 1,kl
    xtm1(:,kl+1-k,nt) = xtg(1:imax,k,nt)
    ! Convert from aerosol concentration outside convective cloud (used by CCAM)
    ! to aerosol concentration inside convective cloud
    xtu(:,kl+1-k,nt) = max(xtg(1:imax,k,nt)-(1.-clcon(:,k))*xtosav(:,k,nt),0.)/max(clcon(:,k),1.E-8)
  end do
end do
do k = 1,kl
  aphp1(:,kl+1-k)  = rhoa(:,k)*dz(:,k)                          ! density * thickness
  prhop1(:,kl+1-k) = rhoa(:,k)                                  ! air density
  ptp1(:,kl+1-k)   = ttg(1:imax,k)                              ! air temperature
  pclcon(:,kl+1-k) = min(max(clcon(:,k),0.),1.)                 ! convective cloud fraction
  qtot = qlg(1:imax,k) + qfg(1:imax,k)                          ! total liquid and ice mixing ratio
  pclcover(:,kl+1-k) = stratcloud(:,k)*qlg(:,k)/max(qtot,1.E-8) ! Liquid-cloud fraction
  pcfcover(:,kl+1-k) = stratcloud(:,k)*qfg(:,k)/max(qtot,1.E-8) ! Ice-cloud fraction
  pmlwc(:,kl+1-k) = qlg(:,k)
  pmiwc(:,kl+1-k) = qfg(:,k)
  where ( k<=kbsav )
    pfconv(:,kl+1-k) = condc(:)/dt
  elsewhere
    pfconv(:,kl+1-k) = 0.
  end where
end do
!fracc = 0.1   ! LDR suggestion (0.1 to 0.3)
fracc = cldcon ! MJT suggestion (use NCAR scheme)
call xtchemie(2, dt, zdayfac, aphp1, pmrate, pfprec,                    & !Inputs
              pclcover, pmlwc, prhop1, ptp1, taudar, xtm1,              & !Inputs
              pfsnow,pfsubl,pcfcover,pmiwc,pmaccr,pfmelt,pfstayice,     & !Inputs
              pqfsedice,plambs,prscav,prfreeze,pclcon,fracc,            & !Inputs
              pccw,pfconv,xtu,                                          & !Inputs
              xte, so2oh, so2h2, so2o3, dmsoh, dmsn3,                   & !Output
              zoxidant,so2wd,so4wd,bcwd,ocwd,dustwd,saltwd,             &
              imax,kl)                                                    !Inputs
do nt = 1,naero
  do k = 1,kl
    xtg(1:imax,k,nt) = max( xtg(1:imax,k,nt)+xte(:,kl+1-k,nt)*dt, 0. )
  end do
enddo
dmsso2o(:) = dmsso2o(:) + dmsoh(:) + dmsn3(:)             ! oxidation of DMS to SO2
so2so4o(:) = so2so4o(:) + so2oh(:) + so2h2(:) + so2o3(:)  ! oxidation of SO2 to SO4

#ifdef debugaero
if ( maxval(xtg(1:imax,:,:))>6.5e-6 ) then
  write(6,*) "xtg out-of-range after xtchemie"
  write(6,*) "xtg maxval,maxloc ",maxval(xtg(1:imax,:,:)),maxloc(xtg(1:imax,:,:))
end if
#endif

do nt = 1,ndust
  burden(:) = sum( xtg(1:imax,:,nt+itracdu-1)*rhoa(:,:)*dz(:,:), dim=2 )
  dust_burden(:,nt) = dust_burden(:,nt) + burden(:)
end do

burden(:) = 0.
do nt = itracbc,itracbc+1
  burden(:) = burden(:) + sum( xtg(1:imax,:,nt)*rhoa(:,:)*dz(:,:), dim=2 )
end do
bc_burden(:) = bc_burden(:) + burden(:)

burden(:) = 0.
do nt = itracoc,itracoc+1
  burden(:) = burden(:) + sum( xtg(1:imax,:,nt)*rhoa(:,:)*dz(:,:), dim=2 )
end do
oc_burden(:) = oc_burden(:) + burden(:)

burden(:) = sum( xtg(1:imax,:,itracdms)*rhoa(:,:)*dz(:,:), dim=2 )
dms_burden(:) = dms_burden(:) + burden(:)

burden(:) = sum( xtg(1:imax,:,itracso2)*rhoa(:,:)*dz(:,:), dim=2 )
so2_burden(:) = so2_burden(:) + burden(:)

burden(:) = sum( xtg(1:imax,:,itracso4)*rhoa(:,:)*dz(:,:), dim=2 )
so4_burden(:) = so4_burden(:) + burden(:)

burden(:) = 0.
do nt = itracsa,itracsa+nsalt-1
  burden(:) = burden(:) + sum( xtg(1:imax,:,nt)*rhoa(:,:)*dz(:,:), dim=2 )
end do
salt_burden(:) = salt_burden(:) + burden(:)

return
end subroutine aldrcalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! xt emiss

SUBROUTINE XTEMISS(ztmst, rhoa, TSM1M, SEAICEM, ZZSPEED,                         & !Inputs
                   LOLAND, PFOREST, PSNOW, WSM1M, dz,                            & !Inputs
                   XTE, PXTEMS, bbem,                                            & !Outputs
                   emissfield,vso2,dmse,so2e,so4e,bce,oce,xtg,so2dd,so4dd,bcdd,ocdd, &
                   imax,kl)                                                     !Inputs
!$acc routine vector

!
!    THIS ROUTINE CALCULATES THE LOWER BOUNDARY CONDITIONS
!    FOR VDIFF DEPENDING ON THE SURFACE EMISSION AND THE
!    DRY DEPOSITION FLUX.
!
!    JOHANN FEICHTER          UNI/HAMBURG         08/91
!    MODIFIED  U. SCHLESE    DKRZ-HAMBURG        JAN-95
!    Adapted for CSIRO GCM by Leon Rotstayn, 12/99
!    Adapted for CCAM by Marcus Thatcher 2012
!
!    PURPOSE
!   ---------
!    THE LOWER BOUNDARY CONDITION FOR CALCULATING THE
!    TURBULENT EXCHANGE IN THE BOUNDARY LAYER IS
!    DETERMINED BY THE EMISSION AND THE DRY DEPOSITION FLUX.
!

implicit none

! Argument list
integer, intent(in) :: imax, kl
real, intent(in) :: ztmst                           !Timestep [s]
real, dimension(imax,kl), intent(in) :: rhoa        !Density of air
real, dimension(imax), intent(in) :: TSM1M          !Surface temp
real, dimension(imax), intent(in) :: SEAICEM        !Sea-ice fraction
real, dimension(imax), intent(in) :: ZZSPEED        !10m wind (corrected to neutral for Nightingale scheme)
real, dimension(imax,kl), intent(in) :: dz          ! layer thickness [m]
real, dimension(imax), intent(in) :: PFOREST        !Fractional vegetation cover
real, dimension(imax), intent(in) :: PSNOW          !Snow depth [m]
! Land-surface details needed to specify dry deposition velocity
real, dimension(imax), intent(in) :: WSM1M          !surface wetness [vol fraction for CSIRO GCM, not m]
real, dimension(imax,kl,naero), intent(out) :: XTE  !Tracer tendencies (kg/kg/s)
real, dimension(imax,naero), intent(out) :: PXTEMS  !Sfc. flux of tracer passed to vertical mixing [kg/m2/s]
logical, dimension(imax), intent(in) :: LOLAND      !Land flag
! Some diagnostics
real, dimension(imax), intent(out) :: bbem

integer jk,jt

real, dimension(imax,2) :: ZVDRD
real, dimension(imax) :: gdp, zdmsemiss
real, dimension(imax) :: zhilbco, zhilbcy, zhiloco, zhilocy
real, dimension(imax) :: zhilso2, zhilso4
real, dimension(imax) :: zdmscon, ZSST, ScDMS, zVdms, wtliss
real, dimension(imax) :: VpCO2, VpCO2liss
real, dimension(imax) :: zvd2ice, zvd4ice, zvd2nof, zvd4nof
real, dimension(imax,15), intent(in) :: emissfield
real, dimension(imax), intent(in) :: vso2
real, dimension(imax), intent(inout) :: dmse
real, dimension(imax), intent(inout) :: so2e
real, dimension(imax), intent(inout) :: so4e
real, dimension(imax), intent(inout) :: bce
real, dimension(imax), intent(inout) :: oce
real, dimension(imax,kl,naero), intent(in) :: xtg
real, dimension(imax), intent(inout) :: so2dd
real, dimension(imax), intent(inout) :: so4dd
real, dimension(imax), intent(inout) :: bcdd
real, dimension(imax), intent(inout) :: ocdd
!

!     M WATER EQUIVALENT  CRITICAL SNOW HEIGHT (FROM *SURF*)
real, parameter :: ZSNCRI = 0.025
!     COEFFICIENTS FOR ZVDRD = FUNCTION OF SOIL MOISTURE
real, parameter :: ZVWC2 = (0.8E-2 - 0.2E-2)/(1. - 0.9)
real, parameter :: ZVW02 = ZVWC2-0.8E-2
real, parameter :: ZVWC4 = (0.2E-2 - 0.025E-2)/(1. - 0.9)
real, parameter :: ZVW04 = ZVWC4-0.2E-2
real, parameter :: tmelt = 273.05
!     Dry deposition
real, parameter :: ZVDPHOBIC = 0.025E-2
!     DMS emissions
real, parameter :: ScCO2     = 600.
!real, parameter :: a_vpco2   = 0.222 ! nightingale (2000)
!real, parameter :: b_vpco2   = 0.333 ! nightingale (2000)
real, parameter :: a_vpco2  = 0.166 ! approx Liss and Merlivat (see nightingale 2000)
real, parameter :: b_vpco2  = 0.133 ! approx Liss and Merlivat (see nightingale 2000)

! Start code : ----------------------------------------------------------

pxtems(:,:) = 0.
xte(:,:,:) = 0.

! --------------------------------------------------------------
!
!*     1.   SURFACE EMISSION.
!           ------- --------
!
!   CALCULATE DMS EMISSIONS FOLLOWING LISS+MERLIVAT
!   DMS SEAWATER CONC. FROM KETTLE ET AL.
ZDMSCON(:) = EMISSFIELD(:,idmso)*(1.-SEAICEM(:))**2
ZSST(:) = min( TSM1M(:)-273.15, 45. )   ! Even Saltzman Sc formula has trouble over 45 deg C
! The formula for ScDMS from Saltzman et al (1993) is given by Kettle & Andreae (ref below)
ScDMS(:) = 2674. - 147.12*ZSST(:) + 3.726*ZSST(:)**2 - 0.038*ZSST(:)**3 !Sc for DMS (Saltzman et al.)
! Nightingale (2000) scheme (J. Biogeochem. Cycles, 14, 373-387)
! For this scheme, zzspeed is the 10m wind adjusted to neutral stability.
VpCO2(:) = a_vpco2*zzspeed(:)*zzspeed(:) + b_vpco2*zzspeed(:) !Nightingale et al
!  ZZSPEED:  10-M WINDS
where ( ZZSPEED(:)<3.6 )
  zVdms(:) = VpCO2(:)*(ScCO2/ScDMS(:))**(2./3.)
elsewhere ( zzspeed(:)<20. )
  ! Phase in Liss & Merlivat from 13 to 18 m/s, since Nightingale is doubtful for high windspeeds,
  ! due to limited data.
  VpCO2liss(:) = 5.9*ZZSPEED(:) - 49.3
  wtliss(:) = min( max( (zzspeed(:)-13.)/5., 0. ), 1. )
  VpCO2(:) = wtliss(:)*VpCO2liss(:) + (1.-wtliss(:))*VpCO2(:)        
  zVdms(:) = VpCO2(:)*sqrt(ScCO2/ScDMS(:))
elsewhere
  ! limit wind speed to 20 m/s for emissions - MJT suggestion  
  VpCO2liss(:) = 5.9*20. - 49.3
  wtliss(:) = 1.
  VpCO2(:) = VpCO2liss(:)
  zVdms(:) = VpCO2(:)*sqrt(ScCO2/ScDMS(:))
end where
where ( loland(:) )
  !zdmsemiss(:) = emissfield(:,idmst) !kg/m2/s
  zdmsemiss(:) = (1./1.938)*emissfield(:,idmst) !kgS/m2/s
elsewhere
  zdmsemiss(:) = ZDMSCON(:)*ZVDMS(:)*32.06e-11/3600.
  ! NANOMOL/LTR*CM/HOUR --> KG/M**2/SEC
end where
jk = 1
gdp(:) = 1./(rhoa(:,jk)*dz(:,jk))
xte(:,jk,itracdms) = xte(:,jk,itracdms) + zdmsemiss(:)*gdp(:)

! Other biomass emissions of SO2 are done below (with the non-surface S emissions)
PXTEMS(:,ITRACSO2)  =(EMISSFIELD(:,iso2a1)+EMISSFIELD(:,iso2b1))*0.97
PXTEMS(:,ITRACSO4)  =(EMISSFIELD(:,iso2a1)+EMISSFIELD(:,iso2b1))*0.03
! Apply these here as a tendency (XTE), rather than as a surface flux (PXTEMS) via vertmix.
jk=1
gdp(:)=1./(rhoa(:,jk)*dz(:,jk))
xte(:,jk,itracso2)  =xte(:,jk,itracso2)  +pxtems(:,itracso2)*gdp
xte(:,jk,itracso4)  =xte(:,jk,itracso4)  +pxtems(:,itracso4)*gdp

!  EMISSION OF ANTHROPOGENIC SO2 IN THE NEXT HIGHER LEVEL PLUS BIOMASS BURNING
do jk=jk2,jk3-1
  gdp=1./(rhoa(:,jk)*dz(:,jk))/real(jk3-jk2)
  XTE(:,JK,ITRACSO2)  =XTE(:,JK,ITRACSO2)  +0.97*EMISSFIELD(:,iso2a2)*gdp !100% of the "above 100m" SO2 emission
  XTE(:,JK,ITRACSO4)  =XTE(:,JK,ITRACSO4)  +0.03*EMISSFIELD(:,iso2a2)*gdp !100% of the "above 100m" SO4 emission
end do
do jk=jk3,jk4-1
  gdp=1./(rhoa(:,jk)*dz(:,jk))/real(jk4-jk3)
  xte(:,jk,ITRACSO2)=xte(:,jk,ITRACSO2)+0.3*emissfield(:,iso2b2)*gdp
end do
do jk=jk4,jk5-1
  gdp=1./(rhoa(:,jk)*dz(:,jk))/real(jk5-jk4)
  xte(:,jk,ITRACSO2)=xte(:,jk,ITRACSO2)+0.4*emissfield(:,iso2b2)*gdp
end do
do jk=jk5,jk6-1
  gdp=1./(rhoa(:,jk)*dz(:,jk))/real(jk6-jk5)
  xte(:,jk,ITRACSO2)=xte(:,jk,ITRACSO2)+0.3*emissfield(:,iso2b2)*gdp
end do
  
!    VOLCANIC BACKGROUND EMISSIONS 
!
!   3 EMISSION LEVELS: 
!    1. PRE-INTRA ERUPTION IN LEVEL IVOLC-HEIGHT (=TOP OF VOLCANO)
!    2. POST-EXTRA ERUPTION IN LEVEL 15 -16 (CA 550-1736M)
!    3. EXPLOSIVE ERUPTION IN LEVEL 10 - 11 (CA 5000-7900M)
jk=1
gdp=1./(rhoa(:,jk)*dz(:,jk))
XTE(:,jk,ITRACSO2)=XTE(:,jk,ITRACSO2)+ZVOLCEMI*0.36*vso2*gdp
do jk=jk4,jk6-1
  gdp=1./(rhoa(:,jk)*dz(:,jk))/real(jk6-jk4)
  XTE(:,jk,ITRACSO2)=XTE(:,jk,ITRACSO2)+ZVOLCEMI*0.36*vso2*gdp
end do
do jk=jk8,jk9-1
  gdp=1./(rhoa(:,jk)*dz(:,jk))/real(jk9-jk8)
  XTE(:,jk,ITRACSO2)=XTE(:,jk,ITRACSO2)+ZVOLCEMI*0.28*vso2*gdp
end do


!Do carbonaceous aerosols
! Inject the low-level fossil-fuel and natural SOA emissions into layer 1
! Assume BC 80% hydrophobic, OC 50%.
PXTEMS(:,ITRACBC)  =0.8*EMISSFIELD(:,ibca1)
PXTEMS(:,ITRACBC+1)=0.2*EMISSFIELD(:,ibca1)
PXTEMS(:,ITRACOC)  =0.5*(EMISSFIELD(:,ioca1)+EMISSFIELD(:,iocna))
PXTEMS(:,ITRACOC+1)=0.5*(EMISSFIELD(:,ioca1)+EMISSFIELD(:,iocna))
! Apply these here as a tendency (XTE), rather than as a surface flux (PXTEMS) via vertmix.
jk=1
gdp=1./(rhoa(:,jk)*dz(:,jk))
xte(:,jk,itracbc)  =xte(:,jk,itracbc)  +pxtems(:,itracbc)*gdp
xte(:,jk,itracbc+1)=xte(:,jk,itracbc+1)+pxtems(:,itracbc+1)*gdp
xte(:,jk,itracoc)  =xte(:,jk,itracoc)  +pxtems(:,itracoc)*gdp
xte(:,jk,itracoc+1)=xte(:,jk,itracoc+1)+pxtems(:,itracoc+1)*gdp
! Inject the upper-level fossil-fuel emissions into layer 2
! Assume BC 80% hydrophobic, OC 50%.
PXTEMS(:,ITRACBC)  =0.8*EMISSFIELD(:,ibca2)
PXTEMS(:,ITRACBC+1)=0.2*EMISSFIELD(:,ibca2)
PXTEMS(:,ITRACOC)  =0.5*EMISSFIELD(:,ioca2)
PXTEMS(:,ITRACOC+1)=0.5*EMISSFIELD(:,ioca2)
! Apply these here as a tendency (XTE), rather than as a surface flux (PXTEMS) via vertmix.
do jk=jk2,jk3-1
  gdp=1./(rhoa(:,jk)*dz(:,jk))/real(jk3-jk2)
  xte(:,jk,itracbc)  =xte(:,jk,itracbc)  +pxtems(:,itracbc)*gdp
  xte(:,jk,itracbc+1)=xte(:,jk,itracbc+1)+pxtems(:,itracbc+1)*gdp
  xte(:,jk,itracoc)  =xte(:,jk,itracoc)  +pxtems(:,itracoc)*gdp
  xte(:,jk,itracoc+1)=xte(:,jk,itracoc+1)+pxtems(:,itracoc+1)*gdp
end do
! Inject the lower-level biomass emissions into layer 2 (NB: Doesn't include biofuel any more)
! Assume BC and OC both 50% hydrophobic.
PXTEMS(:,ITRACBC)  =0.5*EMISSFIELD(:,ibcb1)
PXTEMS(:,ITRACBC+1)=0.5*EMISSFIELD(:,ibcb1)
PXTEMS(:,ITRACOC)  =0.5*EMISSFIELD(:,iocb1)
PXTEMS(:,ITRACOC+1)=0.5*EMISSFIELD(:,iocb1)
! Apply these here as a tendency (XTE)
do jk=jk2,jk3-1
  gdp=1./(rhoa(:,jk)*dz(:,jk))/real(jk3-jk2)
  xte(:,jk,itracbc)  =xte(:,jk,itracbc)  +pxtems(:,itracbc)*gdp
  xte(:,jk,itracbc+1)=xte(:,jk,itracbc+1)+pxtems(:,itracbc+1)*gdp
  xte(:,jk,itracoc)  =xte(:,jk,itracoc)  +pxtems(:,itracoc)*gdp
  xte(:,jk,itracoc+1)=xte(:,jk,itracoc+1)+pxtems(:,itracoc+1)*gdp
end do
! Inject the upper-level biomass emissions into layers 3-5 (30%, 40%, 30%)
! Assume BC and OC both 50% hydrophobic.
PXTEMS(:,ITRACBC)  =0.5*EMISSFIELD(:,ibcb2)
PXTEMS(:,ITRACBC+1)=0.5*EMISSFIELD(:,ibcb2)
PXTEMS(:,ITRACOC)  =0.5*EMISSFIELD(:,iocb2)
PXTEMS(:,ITRACOC+1)=0.5*EMISSFIELD(:,iocb2)
! Apply these here as a tendency (XTE)
do jk=jk3,jk4-1
  gdp=1./(rhoa(:,jk)*dz(:,jk))/real(jk4-jk3)
  xte(:,jk,itracbc)  =xte(:,jk,itracbc)  +0.3*pxtems(:,itracbc)*gdp
  xte(:,jk,itracbc+1)=xte(:,jk,itracbc+1)+0.3*pxtems(:,itracbc+1)*gdp
  xte(:,jk,itracoc)  =xte(:,jk,itracoc)  +0.3*pxtems(:,itracoc)*gdp
  xte(:,jk,itracoc+1)=xte(:,jk,itracoc+1)+0.3*pxtems(:,itracoc+1)*gdp
end do
do jk=jk4,jk5-1
  gdp=1./(rhoa(:,jk)*dz(:,jk))/real(jk5-jk4)
  xte(:,jk,itracbc)  =xte(:,jk,itracbc)  +0.4*pxtems(:,itracbc)*gdp
  xte(:,jk,itracbc+1)=xte(:,jk,itracbc+1)+0.4*pxtems(:,itracbc+1)*gdp
  xte(:,jk,itracoc)  =xte(:,jk,itracoc)  +0.4*pxtems(:,itracoc)*gdp
  xte(:,jk,itracoc+1)=xte(:,jk,itracoc+1)+0.4*pxtems(:,itracoc+1)*gdp
end do
do jk=jk5,jk6-1
  gdp=1./(rhoa(:,jk)*dz(:,jk))/real(jk6-jk5)
  xte(:,jk,itracbc)  =xte(:,jk,itracbc)  +0.3*pxtems(:,itracbc)*gdp
  xte(:,jk,itracbc+1)=xte(:,jk,itracbc+1)+0.3*pxtems(:,itracbc+1)*gdp
  xte(:,jk,itracoc)  =xte(:,jk,itracoc)  +0.3*pxtems(:,itracoc)*gdp
  xte(:,jk,itracoc+1)=xte(:,jk,itracoc+1)+0.3*pxtems(:,itracoc+1)*gdp
end do

!   --------------------------------------------------------------
!
!*      2.    DRY DEPOSITION.
!             --- ----------

!      DRY DEPOSITION OF SO2, SO4

!           - MELTING/NOT MELTING SEAICE-
where ( tsm1m>=(tmelt-0.1) )
  zvd2ice = 0.8E-2
  zvd4ice = 0.2E-2
elsewhere
  zvd2ice = 0.1E-2
  zvd4ice = 0.025E-2
end where

!         -  SNOW/NO SNOW -
where ( PSNOW>ZSNCRI .and. tsm1m>=tmelt )
!            - MELTING/NOT MELTING SNOW -
  ZVD2NOF=0.8E-2
  ZVD4NOF=0.2E-2
elsewhere ( PSNOW>ZSNCRI )
  ZVD2NOF=0.1E-2
  ZVD4NOF=0.025E-2
elsewhere ( tsm1m<=tmelt )
!           -  FROZEN SOIL -
  ZVD2NOF=0.2E-2
  ZVD4NOF=0.025E-2
elsewhere
!           - PARTLY WET -
  ZVD2NOF=max(min(ZVWC2*WSM1M-ZVW02,0.8E-2),0.2E-2)
  ZVD4NOF=max(min(ZVWC4*WSM1M-ZVW04,0.2E-2),0.025E-2)
end where
    
!     -  SEA -
where (.NOT.LOLAND)
!         - SEA ICE -
  ZVDRD(:,1)=(1.-SEAICEM(:))*0.8E-2+SEAICEM(:)*ZVD2ICE(:) !So leads agree with ocean
  ZVDRD(:,2)=(1.-SEAICEM(:))*0.2E-2+SEAICEM(:)*ZVD4ICE(:)
elsewhere
!      - LAND -
  ZVDRD(:,1)=PFOREST(:)*0.8E-2+(1.-PFOREST(:))*ZVD2NOF(:)
  ZVDRD(:,2)=PFOREST(:)*0.2E-2+(1.-PFOREST(:))*ZVD4NOF(:)
end where


! Sulfur emission diagnostic (hard-coded for 3 sulfur variables)
do jk=1,kl
  dmse=dmse+xte(:,jk,ITRACDMS)*rhoa(:,jk)*dz(:,jk)   !Above surface
  so2e=so2e+xte(:,jk,ITRACSO2)*rhoa(:,jk)*dz(:,jk)   !Above surface
  so4e=so4e+xte(:,jk,ITRACSO4)*rhoa(:,jk)*dz(:,jk)   !Above surface
enddo

! Assume that BC and OC emissions are passed in through xte()
do jt=ITRACBC,ITRACBC+1
  do jk=1,kl
    bce=bce+xte(:,jk,jt)*rhoa(:,jk)*dz(:,jk)
  enddo
enddo
do jt=ITRACOC,ITRACOC+1
  do jk=1,kl
    oce=oce+xte(:,jk,jt)*rhoa(:,jk)*dz(:,jk)
  enddo
enddo

! Total biomass burning primary emissions (note 1.3 for organic carbon)
bbem=emissfield(:,ibcb1)+emissfield(:,ibcb2)+1.3*(emissfield(:,iocb1)+emissfield(:,iocb2))

! ZVDRD   DRY DEPOSITION VELOCITY IN M/S
! ZVDRD(JL,1)  FOR SO2 GAS
! ZVDRD(JL,2)  FOR AEROSOLS
gdp=1./(rhoa(:,1)*dz(:,1))

zhilso2=(xtg(1:imax,1,itracso2)+xte(1:imax,1,itracso2)*ztmst)   &
       *(1.-exp(-ztmst*zvdrd(:,1)/dz(:,1)))/(ztmst*gdp)
xte(:,1,ITRACSO2)  =xte(:,1,ITRACSO2)  -zhilso2*gdp
  
zhilso4=(xtg(1:imax,1,itracso4)+xte(1:imax,1,itracso4)*ztmst)   &
       *(1.-exp(-ztmst*zvdrd(:,2)/dz(:,1)))/(ztmst*gdp)
xte(:,1,ITRACSO4)  =xte(:,1,ITRACSO4)  -zhilso4*gdp

ZHILBCO=(xtg(1:imax,1,ITRACBC)+xte(1:imax,1,itracbc)*ztmst)     &
       *(1.-exp(-ztmst*ZVDPHOBIC/dz(:,1)))/(ztmst*gdp)
xte(:,1,itracbc)  =xte(:,1,itracbc)    -zhilbco*gdp

ZHILBCY=(xtg(1:imax,1,ITRACBC+1)+xte(1:imax,1,itracbc+1)*ztmst) &
       *(1.-exp(-ztmst*ZVDRD(:,2)/dz(:,1)))/(ztmst*gdp)
xte(:,1,itracbc+1)=xte(:,1,itracbc+1)  -zhilbcy*gdp

ZHILOCO=(xtg(1:imax,1,ITRACOC)+xte(1:imax,1,itracoc)*ztmst)     &
       *(1.-exp(-ztmst*ZVDPHOBIC/dz(:,1)))/(ztmst*gdp)
xte(:,1,itracoc)  =xte(:,1,itracoc)    -zhiloco*gdp

ZHILOCY=(xtg(1:imax,1,ITRACOC+1)+xte(1:imax,1,itracoc+1)*ztmst) &
       *(1.-exp(-ztmst*ZVDRD(:,2)/dz(:,1)))/(ztmst*gdp)
xte(:,1,itracoc+1)=xte(:,1,itracoc+1)  -zhilocy*gdp

so2dd=so2dd+zhilso2
so4dd=so4dd+zhilso4
bcdd=bcdd+zhilbco+zhilbcy
ocdd=ocdd+zhiloco+zhilocy

return
end subroutine xtemiss

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! xt sink

SUBROUTINE XTSINK(PTMST,xtg,imax,kl)
!$acc routine vector

!
!   *XTSINK*  CALCULATES THE DECREASE OF TRACER CONCENTRATION
!             FOR  A GIVEN HALF-LIFE-TIME
!
!   JOHANN FEICHTER               UNI-HAMBURG    08-91
!
!   PURPOSE
!  ---------------
!   THE MASS MIXING-RATIO OF TRACERS IS MULTIPLIED WITH
!   EXP(ALOG(0.5)*TIME-STEP/HALF-LIFE-TIME).
!

implicit none

integer, intent(in) :: imax, kl
REAL, intent(in) :: PTMST
real, dimension(imax,kl,naero), intent(inout) :: xtg
real zdxtdt
real pqtmst,zfac,zdecay
integer jk, jl

! Start code : ----------------------------------------------------------

PQTMST=1./PTMST
ZFAC=ALOG(0.5)*PTMST

ZDECAY=EXP(ZFAC/86400.) ! 1 day
!!$acc parallel loop collapse(2) copy(xtg(:,:,itracbc:itracbc+1))
DO CONCURRENT (JK=1:kl)
  DO CONCURRENT(JL=1:imax)
    ZDXTDT=xtg(JL,JK,ITRACBC)*(ZDECAY-1.)
    xtg(JL,JK,ITRACBC)   = xtg(JL,JK,ITRACBC)   + ZDXTDT
    xtg(JL,JK,ITRACBC+1) = xtg(JL,JK,ITRACBC+1) - ZDXTDT
  end do
end do
!!$acc end parallel loop

ZDECAY=EXP(ZFAC/86400.) ! 1 day
!!$acc parallel loop collapse(2) copy(xtg(:,:,itracoc:itracoc+1))
DO CONCURRENT (JK=1:kl)
  DO CONCURRENT (JL=1:imax)
    ZDXTDT=xtg(JL,JK,ITRACOC)*(ZDECAY-1.)
    xtg(JL,JK,ITRACOC)   = xtg(JL,JK,ITRACOC)   + ZDXTDT
    xtg(JL,JK,ITRACOC+1) = xtg(JL,JK,ITRACOC+1) - ZDXTDT
  end do
end do
!!$acc end parallel loop

RETURN
END subroutine xtsink

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! xt chemie

SUBROUTINE XTCHEMIE(KTOP, PTMST,zdayfac,rhodz, PMRATEP, PFPREC,                      & !Inputs
                    PCLCOVER, PMLWC, PRHOP1, PTP1, taudar, xtm1,                     & !Inputs
                    pfsnow,pfsubl,pcfcover,pmiwc,pmaccr,pfmelt,pfstayice,            & !Inputs
                    pqfsedice,plambs,prscav,prfreeze,pclcon,fracc,pccw,pfconv,xtu,   & !Inputs
                    xte,so2oh,so2h2,so2o3,dmsoh,dmsn3,                               & !Outputs
                    zoxidant,so2wd,so4wd,bcwd,ocwd,dustwd,saltwd,                    &
                    imax,kl)                                                           !Inputs
!$acc routine vector

! Inputs
! ktop: top level for aerosol processes (set to 1, counting downwards from top)
! ptmst: timestep (seconds; tdt in main program)
! rhodz: density * thickness (si units)
! pmratep: precip formation rate (kg/kg/s)
! pfprec: rainfall flux (entering from above) (kg/m2/s)
! pclcover: liquid-water cloud fraction (input; don't pass in cfrac though)
! pmlwc: liquid-water mixing ratio (kg/kg)
! prhop1: density of air (kg/m3)
! ptp1: temperature (k)
! taudar: fraction of time sunlit (used to determine if daytime)
! xtm1: tracer mixing ratio (kg/kg)
! pfsnow: snowfall flux (entering from above) (kg/m2/s)
! pfsubl: snowfall flux evaporating in layer k (kg/m2/s)
! pcfcover: ice cloud fraction
! pmiwc: ice mixing ratio (kg/kg)
! pmaccr: accretion rate (kg/kg/s)
! pfmelt: snowfall flux melting in layer k (kg/m2/s)
! pfstayice: snowfall flux staying in layer k (kg/m2/s)
! pqfsedice: fractional ice sedimentation in timestep
! plambs: slope (lambda) for snow crystal size distribution (m**-1)
! prscav: fractional rain scavenging rate in time step (needs to be mult. by coll. eff.)
! pclcon: convective cloud fraction
! fracc: Convective rain fraction
! pccw: convective cloud water mixing ratio (kg/kg)
! pfconv: convective rainfall flux (kg/m2/s)
! xtu: tracer mixing ratio in convective updraught (kg/kg)

! Outputs
! xte: tracer tendency (kg/kg/s)
! so2oh: oxidation of SO2 by OH (diagnostic)
! so2h2: oxidation of SO2 by H2O2 (diagnostic)
! so2o3: oxidation of SO2 by O3 (diagnostic)
! dmsoh: oxidation of DMS by OH (diagnostic)
! dmsn3: oxidation of DMS by NO3 (diagnostic)

!**** *XTCHEMIE*  CALCULATES DRY AND WET CHEMISTRY
!
!      J. FEICHTER             UNI HAMBURG    30/06/92
!
!      PURPOSE
!      ---------
!      THIS ROUTINE COMPUTES THE OXIDATION AND THE WET SCAVENGING
!      OF CHEMICAL SPECIES.
!
!**    INTERFACE
!      -----------
!      *XTCHEMIE*   IS CALLED FROM  progcld in CSIRO GCM
!
!      EXTERNALS
!      ------------
!          *XTWETDEP*  CALCULATES THE WET DEPOSITION

implicit none

! Argument list
integer, intent(in) :: imax, kl
integer, intent(in) :: KTOP
real, intent(in) :: PTMST
real, dimension(imax,kl,naero) :: XTM1
real, dimension(imax,kl,naero) :: xtu
REAL rhodz(imax,kl)
REAL PMRATEP(imax,kl)
REAL PFPREC(imax,kl)
REAL PCLCOVER(imax,kl)
REAL PMLWC(imax,kl)
REAL PRHOP1(imax,kl)
REAL PTP1(imax,kl)
real pfsnow(imax,kl)
real pfconv(imax,kl)
real pfsubl(imax,kl)
real pcfcover(imax,kl)
real pmiwc(imax,kl)
real pmaccr(imax,kl)
real pfmelt(imax,kl)
real pfstayice(imax,kl)
real pqfsedice(imax,kl)
real plambs(imax,kl)
real prscav(imax,kl)
real prfreeze(imax,kl)
real pclcon(imax,kl)
real pccw(imax,kl)
real, dimension(imax) :: taudar
real, dimension(imax) :: fracc
real, dimension(imax), intent(in) :: zdayfac
real, dimension(imax,kl,naero), intent(out) :: xte
real, dimension(imax), intent(out) :: dmsoh, dmsn3, so2oh, so2h2, so2o3 !Diagnostic output

! Local work arrays and variables
integer, dimension(imax) :: ZRDAYL
integer jt,jk,jn
integer jl
real, dimension(imax,kl,naero) :: xto
real, dimension(imax,kl) :: so2oh3d, dmsoh3d, dmsn33d
real, dimension(imax,kl) :: ZXTP10, ZXTP1C, ZHENRY, ZSO4, ZSO4i, ZSO4C, ZHENRYC, ZXTP1CON, zsolub
real, dimension(imax,kl) :: ZZOH, ZZH2O2, ZZO3, ZZNO2
real, dimension(imax,kl) :: zlwcic, ziwcic
real, dimension(imax,kl,4), intent(in) :: zoxidant
real, dimension(imax), intent(inout) :: so2wd
real, dimension(imax), intent(inout) :: so4wd
real, dimension(imax), intent(inout) :: bcwd
real, dimension(imax), intent(inout) :: ocwd
real, dimension(imax,ndust), intent(inout) :: dustwd
real, dimension(imax), intent(inout) :: saltwd
real x,pqtmst
real zlwcl, zlwcv, zhp, zqtp1, zrk, zrke, zxtp1
real zh_so2, zpfac, zp_so2, zf_so2, zh_h2o2, zp_h2o2, zf_h2o2
real ZRKH2O2
real ze1,ze2,ze3,zfac1,zrkfac
real zza,za21,za22,zph_o3,zf_o3,zdt
real zh2o2m,zso2m,zso4m,zsumh2o2,zsumo3
real zq,zso2mh,zdso2h,zso2l,zso4l
real zzb,zzp,zzq,zzp2,zqhp,za2,zheneff
real zrko3,zso2mo,zdso2o,zdso2tot,zfac
real zxtp1dms,zso2,ztk23b
real zhil,zexp,zm,zdms,t,ztk1,zqt,zqt3
real zrhoair,zkno2o3,zkn2o5aq,zrx1,zrx12
real zkno2no3,ztk3,ztk2,zkn2o5
real zno3,zxtp1so2


!    REACTION RATE SO2-OH
real, parameter :: ZK2I=2.0E-12
real, parameter :: ZK2=4.0E-31
real, parameter :: ZK2F=0.45

!   REACTION RATE DMS-NO3
real, parameter :: ZK3=1.9E-13

!   MOLECULAR WEIGHTS IN G
real, parameter :: ZMOLGS=32.064
real, parameter :: ZMOLGAIR=28.84

real, parameter :: ZHPBASE=2.5E-06
real, parameter :: ZE1K=1.1E-02
real, parameter :: ZE1H=2300.
real, parameter :: ZE2K=1.23
real, parameter :: ZE2H=3020.
real, parameter :: ZE3K=1.2E-02
real, parameter :: ZE3H=2010.
real, parameter :: ZQ298=1./298.
real, parameter :: ZRGAS=8.2E-02

real, parameter :: ZAVO=6.022E+23
real, parameter :: ZNAMAIR=1.E-03*ZAVO/ZMOLGAIR

! Start code : ----------------------------------------------------------
dmsoh(:)=0.
dmsn3(:)=0.
so2oh(:)=0.
so2h2(:)=0.
so2o3(:)=0.
do jk = 1,kl
  so2oh3d(:,jk)=0.
  dmsoh3d(:,jk)=0.
  dmsn33d(:,jk)=0.
end do
do jt = 1,naero
  do jk = 1,kl
    xte(:,jk,jt)=0.
  end do
end do
where ( taudar(:)>0.5 )
  zrdayl(:)=1
elsewhere
  zrdayl(:)=0  
end where

! Calculate xto, tracer mixing ratio outside convective updraughts
! Assumes pclcon < 1, but this shouldn't be a problem.
do jt = 1,naero
  do jk = 1,kl
    xto(:,jk,jt)=(xtm1(:,jk,jt)-pclcon(:,jk)*xtu(:,jk,jt))/(1.-pclcon(:,jk))
    xto(:,jk,jt)=max(0.,xto(:,jk,jt))
  end do
end do

#ifdef debugaero
if ( maxval(xtm1(1:imax,:,:)+xte(1:imax,:,:)*PTMST)>6.5e-6 ) then
  write(6,*) "xtg is out-of-range at start of xtchemie"
  write(6,*) "xtg maxval,maxloc ",maxval(xtm1(1:imax,:,:)+xte(1:imax,:,:)*PTMST), &
                                  maxloc(xtm1(1:imax,:,:)+xte(1:imax,:,:)*PTMST)
end if
#endif

!   CALCULATE THE ZRDAYL (=0 --> NIGHT; =1 --> DAY) AND
!                 ZAMUO  =  ZENITH ANGLE

!    CONSTANTS
PQTMST=1./PTMST

do jk = 1,kl

  ! Calculate in-cloud ql
  where ( pclcover(:,jk)>1.e-8 )
    zlwcic(:,jk)=pmlwc(:,jk)/pclcover(:,jk)
  elsewhere
    zlwcic(:,jk)=0.
  end where
  where ( pcfcover(:,jk)>1.e-8 )
    ziwcic(:,jk)=pmiwc(:,jk)/pcfcover(:,jk)
  elsewhere
    ziwcic(:,jk)=0.
  end where

  !  OXIDANT CONCENTRATIONS IN MOLECULE/CM**3
  ! -- levels are already inverted --
  ZZOH(:,jk)   = ZOXIDANT(:,jk,1)
  ZZH2O2(:,jk) = ZOXIDANT(:,jk,2)*PRHOP1(:,jk)*1.e-3
  ZZO3(:,jk)   = ZOXIDANT(:,jk,3)*PRHOP1(:,jk)*1.e-3
  ZZNO2(:,jk)  = ZOXIDANT(:,jk,4)*PRHOP1(:,jk)*1.e-3
  
end do

do jk = 1,kl
  zhenry(:,jk)=0.
  zhenryc(:,jk)=0.
end do

 !   PROCESSES WHICH ARE DIFERENT INSIDE AND OUTSIDE OF CLOUDS
do jk = 1,kl
  ZSO4(:,jk)=amax1(XTO(:,jk,ITRACSO4),0.)
end do

!!$acc data create(prhop1,ptp1,zzh2o2,zzo3,rhodz)
!!$acc update device(prhop1,ptp1,zzh2o2,zzo3,rhodz)

!!$acc parallel loop collapse(2) copy(zhenry,zxtp10,zxtp1c,zso4,so2h2,so2o3) &
!!$acc   copyin(zlwcic,xto(:,:,itracso2),pclcover) present(prhop1,ptp1,zzh2o2,zzo3,rhodz)
DO CONCURRENT (JK=KTOP:kl)
  DO CONCURRENT (JL=1:imax)
    !   CALCULATE THE REACTION-RATES FOR SO2-H2O2
    if ( zlwcic(jl,jk)>zmin ) then
      ZLWCL=ZLWCIC(JL,JK)*PRHOP1(JL,JK)*1.E-06
      ZLWCV=ZLWCIC(JL,JK)*PRHOP1(JL,JK)*1.E-03
      ZHP=ZHPBASE+ZSO4(JL,JK)*1000./(ZLWCIC(JL,JK)*ZMOLGS)
      ZQTP1=1./PTP1(JL,JK)-ZQ298
      ZRK=8.E+04*EXP(-3650.*ZQTP1)/(0.1+ZHP)
      ZRKE=ZRK/(ZLWCL*ZAVO)

      ZH_SO2=ZE2K*EXP(ZE2H*ZQTP1)
      ZPFAC=ZRGAS*ZLWCV*PTP1(JL,JK)
      ZP_SO2=ZH_SO2*ZPFAC
      ZF_SO2=ZP_SO2/(1.+ZP_SO2)

      ZH_H2O2=9.7E+04*EXP(6600.*ZQTP1)
      ZP_H2O2=ZH_H2O2*ZPFAC
      ZF_H2O2=ZP_H2O2/(1.+ZP_H2O2)

      ZRKH2O2=ZRKE*ZF_SO2*ZF_H2O2
    else
      ZRKH2O2=0.
    end if

    !   HETEROGENEOUS CHEMISTRY
    ZXTP1         = XTO(JL,JK,ITRACSO2)
    ZXTP10(JL,JK) = XTO(JL,JK,ITRACSO2)
    ZXTP1C(JL,JK) = XTO(JL,JK,ITRACSO2)
    IF ( ZXTP1>ZMIN .AND. ZLWCIC(JL,JK)>ZMIN ) THEN

      ZQTP1=1./PTP1(JL,JK)-ZQ298
      ZE1=ZE1K*EXP(ZE1H*ZQTP1)
      ZE2=ZE2K*EXP(ZE2H*ZQTP1)
      ZE3=ZE3K*EXP(ZE3H*ZQTP1)

      ZLWCL=ZLWCIC(JL,JK)*PRHOP1(JL,JK)*1.E-06
      !    ZLWCL = LWC IN L/CM**3
      ZLWCV=ZLWCIC(JL,JK)*PRHOP1(JL,JK)*1.E-03
      !   ZLWCV = LWC IN VOL/VOL
      ZFAC1=1./(ZLWCL*ZAVO)
      !   ZFAC1 CALCULATES MOLECULES PER CM**3 TO MOLE PER LTR H2O
      ZRKFAC=ZRGAS*PTP1(JL,JK)*ZLWCV
      !   ZRKFAC CALCULATES DIMENSIONLESS HENRY-COEFF.
      ZZA=ZE2*ZRKFAC
      ZA21=4.39E+11*EXP(-4131./PTP1(JL,JK))
      ZA22=2.56E+03*EXP(-966./PTP1(JL,JK)) !926 corrected to 966 here
      ZPH_O3=ZE1*ZRKFAC
      ZF_O3=ZPH_O3/(1.+ZPH_O3)
      ZDT=PTMST/5.

      ZH2O2M=ZZH2O2(JL,JK)
      ZSO2M=ZXTP1*PRHOP1(JL,JK)*6.022E+20/ZMOLGS
      ZSO4M=ZSO4(JL,JK)*PRHOP1(JL,JK)*6.022E+20/ZMOLGS

      ZSUMH2O2=0.
      ZSUMO3=0.

      DO JN=1,5
        ZQ=ZRKH2O2*ZH2O2M
        ZSO2MH=ZSO2M*EXP(-ZQ*ZDT)

        ZDSO2H=ZSO2M-ZSO2MH
        ZH2O2M=ZH2O2M-ZDSO2H
        ZH2O2M=AMAX1(0.,ZH2O2M)
        ZSUMH2O2=ZSUMH2O2+ZDSO2H

        ZSO4M=ZSO4M+ZDSO2H
        !   CALCULATE THE PH VALUE
        ZSO2L=ZSO2MH*ZFAC1
        ZSO4L=ZSO4M*ZFAC1
        ZZB=ZHPBASE+ZSO4L
        ZZP=(ZZA*ZE3-ZZB-ZZA*ZZB)/(1.+ZZA)
        ZZQ=-ZZA*ZE3*(ZZB+ZSO2L)/(1.+ZZA)
        ZZP=0.5*ZZP
        ZZP2=ZZP*ZZP
        ZHP=-ZZP+SQRT(ZZP2-ZZQ)
        ZQHP=1./ZHP

        !   CALCULATE THE REACTION RATE FOR SO2-O3
        ZA2=(ZA21+ZA22*ZQHP)*ZFAC1
        ZHENEFF=1.+ZE3*ZQHP
        ZP_SO2=ZZA*ZHENEFF
        ZF_SO2=ZP_SO2/(1.+ZP_SO2)
        ZRKO3=ZA2*ZF_O3*ZF_SO2

        ZQ=ZZO3(JL,JK)*ZRKO3
        ZSO2MO=ZSO2MH*EXP(-ZQ*ZDT)
        ZDSO2O=ZSO2MH-ZSO2MO
        ZSO4M=ZSO4M+ZDSO2O
        ZSO2M=ZSO2MO
        ZSUMO3=ZSUMO3+ZDSO2O
      end do  !End of iteration loop

      ZDSO2TOT=ZXTP1-ZSO2M*ZMOLGS/(6.022E+20*PRHOP1(JL,JK))
      ZDSO2TOT=AMIN1(ZDSO2TOT,ZXTP1)
      ZXTP1C(JL,JK)=ZXTP1-ZDSO2TOT
      ZSO4(JL,JK)=ZSO4(JL,JK)+ZDSO2TOT

      ZHENRY(JL,JK)=ZF_SO2
      ! Diagnostic only...
      ZFAC=PQTMST*PCLCOVER(JL,JK)*ZMOLGS/(6.022E+20*PRHOP1(JL,JK))
      ZFAC1=ZFAC*rhodz(JL,JK)
      so2h2(JL)=so2h2(JL)+ZSUMH2O2*ZFAC1
      so2o3(JL)=so2o3(JL)+ZSUMO3*ZFAC1
    END IF
  end do
end do
!!$acc end parallel loop


! Repeat the aqueous oxidation calculation for ice clouds.
do jk = 1,kl
  ZSO4i(:,jk)=amax1(XTO(:,jk,ITRACSO4),0.)
end do

!******************************************************************************
!   CALCULATE THE REACTION-RATES FOR SO2-H2O2
!DO JK=KTOP,KL
!  DO JL=1,imax
!    IF(ziwcic(JL,JK).GT.ZMIN) THEN
!      ZLWCL(jl)=ziwcic(JL,JK)*PRHOP1(JL,JK)*1.E-06
!      ZLWCV(jl)=ziwcic(JL,JK)*PRHOP1(JL,JK)*1.E-03
!      ZHP(jl)=ZHPBASE+ZSO4i(JL,JK)*1000./(ziwcic(JL,JK)*ZMOLGS)
!      ZQTP1(jl)=1./PTP1(JL,JK)-ZQ298
!      ZRK(jl)=8.E+04*EXP(-3650.*ZQTP1(jl))/(0.1+ZHP(jl))
!      ZRKE(jl)=ZRK(jl)/(ZLWCL(jl)*ZAVO)
!
!      ZH_SO2(jl)=ZE2*EXP(ZE2H*ZQTP1(jl))
!      ZPFAC(jl)=ZRGAS*ZLWCV(jl)*PTP1(JL,JK)
!      ZP_SO2(jl)=ZH_SO2(jl)*ZPFAC(jl)
!      ZF_SO2(jl)=ZP_SO2(jl)/(1.+ZP_SO2(jl))
!
!      ZH_H2O2(jl)=9.7E+04*EXP(6600.*ZQTP1(jl))
!      ZP_H2O2(jl)=ZH_H2O2(jl)*ZPFAC(jl)
!      ZF_H2O2(jl)=ZP_H2O2(jl)/(1.+ZP_H2O2(jl))
!
!      ZRKH2O2(JL,JK)=ZRKE(jl)*ZF_SO2(jl)*ZF_H2O2(jl)
!    ELSE
!      ZRKH2O2(JL,JK)=0.
!    ENDIF
!  ENDDO
!
!!   HETEROGENEOUS CHEMISTRY
!  DO JL=1,imax
!    ZXTP1(jl)=XTO(JL,JK,ITRACSO2)
!    IF(ZXTP1(jl)>ZMIN.AND.ziwcic(JL,JK)>ZMIN) THEN
!      X=PRHOP1(JL,JK)
!
!      ZQTP1(jl)=1./PTP1(JL,JK)-ZQ298
!      ZE1=ZE1K*EXP(ZE1H*ZQTP1(jl))
!      ZE2=ZE2K*EXP(ZE2H*ZQTP1(jl))
!      ZE3=ZE3K*EXP(ZE3H*ZQTP1(jl))
!
!      ZLWCL(jl)=ziwcic(JL,JK)*PRHOP1(JL,JK)*1.E-06
!!    ZLWCL = LWC IN L/CM**3
!      ZLWCV(jl)=ziwcic(JL,JK)*PRHOP1(JL,JK)*1.E-03
!!   ZLWCV = LWC IN VOL/VOL
!      ZFAC1=1./(ZLWCL(jl)*ZAVO)
!!   ZFAC1 CALCULATES MOLECULES PER CM**3 TO MOLE PER LTR H2O
!      ZRKFAC=ZRGAS*PTP1(JL,JK)*ZLWCV(jl)
!!   ZRKFAC CALCULATES DIMENSIONLESS HENRY-COEFF.
!      ZZA=ZE2*ZRKFAC
!      ZA21=4.39E+11*EXP(-4131./PTP1(JL,JK))
!      ZA22=2.56E+03*EXP(-966./PTP1(JL,JK)) !926 corrected to 966 here
!      ZPH_O3=ZE1*ZRKFAC
!      ZF_O3=ZPH_O3/(1.+ZPH_O3)
!      ZDT=PTMST/5.
!
!      ZH2O2M=ZZH2O2(JL,JK)
!      ZSO2M=ZXTP1(jl)*X*6.022E+20/ZMOLGS
!      ZSO4M=ZSO4i(JL,JK)*X*6.022E+20/ZMOLGS
!
!      ZSUMH2O2=0.
!      ZSUMO3=0.
!
!      DO JN=1,5
!        ZQ=ZRKH2O2(JL,JK)*ZH2O2M
!        ZSO2MH=ZSO2M*EXP(-ZQ*ZDT)
!
!        ZDSO2H=ZSO2M-ZSO2MH
!        ZH2O2M=ZH2O2M-ZDSO2H
!        ZH2O2M=AMAX1(0.,ZH2O2M)
!        ZSUMH2O2=ZSUMH2O2+ZDSO2H
!
!        ZSO4M=ZSO4M+ZDSO2H
!!   CALCULATE THE PH VALUE
!        ZSO2L=ZSO2MH*ZFAC1
!        ZSO4L=ZSO4M*ZFAC1
!        ZZB=ZHPBASE+ZSO4L
!       ZZP=(ZZA*ZE3-ZZB-ZZA*ZZB)/(1.+ZZA)
!        ZZQ=-ZZA*ZE3*(ZZB+ZSO2L)/(1.+ZZA)
!        ZZP=0.5*ZZP
!        ZZP2=ZZP*ZZP
!        ZHP(jl)=-ZZP+SQRT(ZZP2-ZZQ)
!        ZQHP=1./ZHP(jl)
!
!!   CALCULATE THE REACTION RATE FOR SO2-O3
!        ZA2=(ZA21+ZA22*ZQHP)*ZFAC1
!        ZHENEFF=1.+ZE3*ZQHP
!        ZP_SO2(jl)=ZZA*ZHENEFF
!        ZF_SO2(jl)=ZP_SO2(jl)/(1.+ZP_SO2(jl))
!        ZRKO3=ZA2*ZF_O3*ZF_SO2(jl)
!
!        ZQ=ZZO3(JL,JK)*ZRKO3
!        ZSO2MO=ZSO2MH*EXP(-ZQ*ZDT)
!        ZDSO2O=ZSO2MH-ZSO2MO
!        ZSO4M=ZSO4M+ZDSO2O
!        ZSO2M=ZSO2MO
!        ZSUMO3=ZSUMO3+ZDSO2O
!      ENDDO  !End of iteration loop
!
!      ZDSO2TOT=ZXTP1(jl)-ZSO2M*ZMOLGS/(6.022E+20*X)
!      ZDSO2TOT=AMIN1(ZDSO2TOT,ZXTP1(jl))
!
!      ZXTP10(JL,JK)=ZXTP1(jl)-ZDSO2TOT*pcfcover(jl,jk)/(1.-pclcover(jl,jk))
!      ZSO4i(JL,JK)=ZSO4i(JL,JK)+ZDSO2TOT*pcfcover(jl,jk)/(1.-pclcover(jl,jk))
!      ZHENRY(JL,JK)=ZF_SO2(jl)
!! Diagnostic only...
!      ZFAC=PQTMST*pcfcover(jl,jk)*ZMOLGS/(6.022E+20*X)
!      ZFAC1=ZFAC*rhodz(JL,JK)
!      so2h2(JL)=so2h2(JL)+ZSUMH2O2*ZFAC1
!      so2o3(JL)=so2o3(JL)+ZSUMO3*ZFAC1
!    ENDIF
!  ENDDO
!ENDDO
!******************************************************************************


! Repeat the aqueous oxidation calculation for convective clouds.
do jk = 1,kl
  ZXTP1CON(:,jk)=amax1(XTU(:,jk,ITRACSO2),0.)
  ZSO4C(:,jk)   =amax1(XTU(:,jk,ITRACSO4),0.)
end do

!!$acc parallel loop collapse(2) copy(zhenryc,zxtp1con,zso4c,so2h2,so2o3) &
!!$acc   copyin(pccw,xtu(:,:,itracso2),pclcon) present(prhop1,ptp1,zzh2o2,zzo3,rhodz)
DO CONCURRENT (JK=KTOP:kl)
  DO CONCURRENT (JL=1:imax)
    !   CALCULATE THE REACTION-RATES FOR SO2-H2O2
    if ( PCCW(JL,JK)>ZMIN ) then
      ZLWCL=PCCW(JL,JK)*PRHOP1(JL,JK)*1.E-06
      ZLWCV=PCCW(JL,JK)*PRHOP1(JL,JK)*1.E-03
      ZHP=ZHPBASE+ZSO4C(JL,JK)*1000./(PCCW(JL,JK)*ZMOLGS)
      ZQTP1=1./PTP1(JL,JK)-ZQ298
      ZRK=8.E+04*EXP(-3650.*ZQTP1)/(0.1+ZHP)
      ZRKE=ZRK/(ZLWCL*ZAVO)

      ZH_SO2=ZE2K*EXP(ZE2H*ZQTP1)
      ZPFAC=ZRGAS*ZLWCV*PTP1(JL,JK)
      ZP_SO2=ZH_SO2*ZPFAC
      ZF_SO2=ZP_SO2/(1.+ZP_SO2)

      ZH_H2O2=9.7E+04*EXP(6600.*ZQTP1)
      ZP_H2O2=ZH_H2O2*ZPFAC
      ZF_H2O2=ZP_H2O2/(1.+ZP_H2O2)

      ZRKH2O2=ZRKE*ZF_SO2*ZF_H2O2
    else
      ZRKH2O2=0.
    end if

    !   HETEROGENEOUS CHEMISTRY
    ZXTP1=XTU(JL,JK,ITRACSO2)
    IF(ZXTP1>ZMIN.AND.PCCW(JL,JK)>ZMIN) THEN

      ZQTP1=1./PTP1(JL,JK)-ZQ298
      ZE1=ZE1K*EXP(ZE1H*ZQTP1)
      ZE2=ZE2K*EXP(ZE2H*ZQTP1)
      ZE3=ZE3K*EXP(ZE3H*ZQTP1)

      ZLWCL=PCCW(JL,JK)*PRHOP1(JL,JK)*1.E-06
      !    ZLWCL = LWC IN L/CM**3
      ZLWCV=PCCW(JL,JK)*PRHOP1(JL,JK)*1.E-03
      !   ZLWCV = LWC IN VOL/VOL
      ZFAC1=1./(ZLWCL*ZAVO)
      !   ZFAC1 CALCULATES MOLECULES PER CM**3 TO MOLE PER LTR H2O
      ZRKFAC=ZRGAS*PTP1(JL,JK)*ZLWCV
      !   ZRKFAC CALCULATES DIMENSIONLESS HENRY-COEFF.
      ZZA=ZE2*ZRKFAC
      ZA21=4.39E+11*EXP(-4131./PTP1(JL,JK))
      ZA22=2.56E+03*EXP(-966./PTP1(JL,JK)) !926 corrected to 966 here
      ZPH_O3=ZE1*ZRKFAC
      ZF_O3=ZPH_O3/(1.+ZPH_O3)
      ZDT=PTMST/5.

      ZH2O2M=ZZH2O2(JL,JK)
      ZSO2M=ZXTP1*PRHOP1(JL,JK)*6.022E+20/ZMOLGS
      ZSO4M=ZSO4C(JL,JK)*PRHOP1(JL,JK)*6.022E+20/ZMOLGS

      ZSUMH2O2=0.
      ZSUMO3=0.

      DO JN=1,5
        ZQ=ZRKH2O2*ZH2O2M
        ZSO2MH=ZSO2M*EXP(-ZQ*ZDT)

        ZDSO2H=ZSO2M-ZSO2MH
        ZH2O2M=ZH2O2M-ZDSO2H
        ZH2O2M=AMAX1(0.,ZH2O2M)
        ZSUMH2O2=ZSUMH2O2+ZDSO2H

        ZSO4M=ZSO4M+ZDSO2H
!   CALCULATE THE PH VALUE
        ZSO2L=ZSO2MH*ZFAC1
        ZSO4L=ZSO4M*ZFAC1
        ZZB=ZHPBASE+ZSO4L
        ZZP=(ZZA*ZE3-ZZB-ZZA*ZZB)/(1.+ZZA)
        ZZQ=-ZZA*ZE3*(ZZB+ZSO2L)/(1.+ZZA)
        ZZP=0.5*ZZP
        ZZP2=ZZP*ZZP
        ZHP=-ZZP+SQRT(ZZP2-ZZQ)
        ZQHP=1./ZHP

!   CALCULATE THE REACTION RATE FOR SO2-O3
        ZA2=(ZA21+ZA22*ZQHP)*ZFAC1
        ZHENEFF=1.+ZE3*ZQHP
        ZP_SO2=ZZA*ZHENEFF
        ZF_SO2=ZP_SO2/(1.+ZP_SO2)
        ZRKO3=ZA2*ZF_O3*ZF_SO2
!
        ZQ=ZZO3(JL,JK)*ZRKO3
        ZSO2MO=ZSO2MH*EXP(-ZQ*ZDT)
        ZDSO2O=ZSO2MH-ZSO2MO
        ZSO4M=ZSO4M+ZDSO2O
        ZSO2M=ZSO2MO
        ZSUMO3=ZSUMO3+ZDSO2O
      ENDDO  !End of iteration loop

      ZDSO2TOT=ZXTP1-ZSO2M*ZMOLGS/(6.022E+20*PRHOP1(JL,JK))
      ZDSO2TOT=AMIN1(ZDSO2TOT,ZXTP1)
      ZXTP1CON(JL,JK)=ZXTP1CON(JL,JK)-ZDSO2TOT
      ZSO4C(JL,JK)=ZSO4C(JL,JK)+ZDSO2TOT
      ZHENRYC(JL,JK)=ZF_SO2
      ! Diagnostic only...
      ZFAC=PQTMST*pclcon(jl,jk)*ZMOLGS/(6.022E+20*PRHOP1(JL,JK))
      ZFAC1=ZFAC*rhodz(JL,JK)
      so2h2(JL)=so2h2(JL)+ZSUMH2O2*ZFAC1
      so2o3(JL)=so2o3(JL)+ZSUMO3*ZFAC1
    ENDIF
  ENDDO
ENDDO
!!$acc end parallel loop

!!$acc end data

!*******************************************************************************
!
!    CALCULATE THE WET DEPOSITION
!    (True for all except DMS)
!

!!$acc data create(pclcover,pclcon,rhodz,pmratep,pfprec,pmlwc,ptp1,pfsnow,pfsubl,pcfcover,pmaccr, &
!!$acc   pfmelt,pfstayice,pqfsedice,plambs,prscav,prfreeze,pfconv,fracc)
!!$acc update device(pclcover,pclcon,rhodz,pmratep,pfprec,pmlwc,ptp1,pfsnow,pfsubl,pcfcover,pmaccr, &
!!$acc   pfmelt,pfstayice,pqfsedice,plambs,prscav,prfreeze,pfconv,fracc)

JT=ITRACSO2
zsolub(:,:)=zhenry(:,:)
CALL XTWETDEP( JT,                                   &
               PTMST,                                &
               rhodz,                                &
               PMRATEP, PFPREC,                      &
               PCLCOVER, zsolub, pmlwc, ptp1,        &
               pfsnow,pfsubl,pcfcover,pmaccr,pfmelt, &
               pfstayice,pqfsedice,plambs,           &
               prscav,prfreeze,pfconv,pclcon,        & 
               fracc,                                & !Inputs
               ZXTP10, ZXTP1C, ZXTP1CON,             &
               xtm1(:,:,jt),xte(:,:,jt),so2wd,imax,kl)

JT=ITRACSO4
zxtp10(:,:)=zso4i(:,:)
zxtp1c(:,:)=zso4(:,:)    
zxtp1con(:,:)=zso4c(:,:)
zsolub (:,:)=0.6
CALL XTWETDEP( JT,                                   &
               PTMST,                                &
               rhodz,                                &
               PMRATEP, PFPREC,                      &
               PCLCOVER, zsolub, pmlwc, ptp1,        &
               pfsnow,pfsubl,pcfcover,pmaccr,pfmelt, &
               pfstayice,pqfsedice,plambs,           &
               prscav,prfreeze,pfconv,pclcon,        & 
               fracc,                                & !Inputs
               ZXTP10, ZXTP1C, ZXTP1CON,             &
               xtm1(:,:,jt),xte(:,:,jt),so4wd,imax,kl)

JT=ITRACbc
zxtp10(:,:)=xto(:,:,jt)
zxtp1c(:,:)=xto(:,:,jt)
zxtp1con(:,:)=xtu(:,:,jt)
zsolub(:,:)=0.
CALL XTWETDEP( JT,                                   &
               PTMST,                                &
               rhodz,                                &
               PMRATEP, PFPREC,                      &
               PCLCOVER, zsolub, pmlwc, ptp1,        &
               pfsnow,pfsubl,pcfcover,pmaccr,pfmelt, &
               pfstayice,pqfsedice,plambs,           &
               prscav,prfreeze,pfconv,pclcon,        & 
               fracc,                                & !Inputs
               ZXTP10, ZXTP1C, ZXTP1CON,             &
               xtm1(:,:,jt),xte(:,:,jt),bcwd,imax,kl)

JT=ITRACBC+1
zxtp10(:,:)=xto(:,:,jt)
zxtp1c(:,:)=xto(:,:,jt)
zxtp1con(:,:)=xtu(:,:,jt)
zsolub(:,:)=0.2
CALL XTWETDEP( JT,                                   &
               PTMST,                                &
               rhodz,                                &
               PMRATEP, PFPREC,                      &
               PCLCOVER, zsolub, pmlwc, ptp1,        &
               pfsnow,pfsubl,pcfcover,pmaccr,pfmelt, &
               pfstayice,pqfsedice,plambs,           &
               prscav,prfreeze,pfconv,pclcon,        & 
               fracc,                                & !Inputs
               ZXTP10, ZXTP1C, ZXTP1CON,             &
               xtm1(:,:,jt),xte(:,:,jt),bcwd,imax,kl)

JT=ITRACOC
zxtp10(:,:)=xto(:,:,jt)
zxtp1c(:,:)=xto(:,:,jt)
zxtp1con(:,:)=xtu(:,:,jt)
zsolub(:,:)=0.
CALL XTWETDEP( JT,                                   &
               PTMST,                                &
               rhodz,                                &
               PMRATEP, PFPREC,                      &
               PCLCOVER, zsolub, pmlwc, ptp1,        &
               pfsnow,pfsubl,pcfcover,pmaccr,pfmelt, &
               pfstayice,pqfsedice,plambs,           &
               prscav,prfreeze,pfconv,pclcon,        & 
               fracc,                                & !Inputs
               ZXTP10, ZXTP1C, ZXTP1CON,             &
               xtm1(:,:,jt),xte(:,:,jt),ocwd,imax,kl)

JT=ITRACOC+1
zxtp10(:,:)=xto(:,:,jt)
zxtp1c(:,:)=xto(:,:,jt)
zxtp1con(:,:)=xtu(:,:,jt)
zsolub(:,:)=0.2
CALL XTWETDEP( JT,                                   &
               PTMST,                                &
               rhodz,                                &
               PMRATEP, PFPREC,                      &
               PCLCOVER, zsolub, pmlwc, ptp1,        &
               pfsnow,pfsubl,pcfcover,pmaccr,pfmelt, &
               pfstayice,pqfsedice,plambs,           &
               prscav,prfreeze,pfconv,pclcon,        & 
               fracc,                                & !Inputs
               ZXTP10, ZXTP1C, ZXTP1CON,             &
               xtm1(:,:,jt),xte(:,:,jt),ocwd,imax,kl)

DO JT=ITRACDU,ITRACDU+NDUST-1
  zxtp10(:,:)=xto(:,:,jt)
  zxtp1c(:,:)=xto(:,:,jt)
  zxtp1con(:,:)=xtu(:,:,jt)
  zsolub(:,:)=0.05
  CALL XTWETDEP( JT,                                   &
                 PTMST,                                &
                 rhodz,                                &
                 PMRATEP, PFPREC,                      &
                 PCLCOVER, zsolub, pmlwc, ptp1,        &
                 pfsnow,pfsubl,pcfcover,pmaccr,pfmelt, &
                 pfstayice,pqfsedice,plambs,           &
                 prscav,prfreeze,pfconv,pclcon,        & 
                 fracc,                                & !Inputs
                 ZXTP10, ZXTP1C, ZXTP1CON,             &
                 xtm1(:,:,jt),xte(:,:,jt),             &
                 dustwd(:,jt-itracdu+1),imax,kl)
end do

DO JT=ITRACSA,ITRACSA+NSALT-1
  zxtp10(:,:)=xto(:,:,jt)
  zxtp1c(:,:)=xto(:,:,jt)
  zxtp1con(:,:)=xtu(:,:,jt)
  zsolub(:,:)=0.05
  CALL XTWETDEP( JT,                                   &
                 PTMST,                                &
                 rhodz,                                &
                 PMRATEP, PFPREC,                      &
                 PCLCOVER, zsolub, pmlwc, ptp1,        &
                 pfsnow,pfsubl,pcfcover,pmaccr,pfmelt, &
                 pfstayice,pqfsedice,plambs,           &
                 prscav,prfreeze,pfconv,pclcon,        & 
                 fracc,                                & !Inputs
                 ZXTP10, ZXTP1C, ZXTP1CON,             &
                 xtm1(:,:,jt),xte(:,:,jt),saltwd,imax,kl)
end do

!!$acc end data

#ifdef debugaero
if ( maxval(xtm1(1:imax,:,:)+xte(1:imax,:,:)*PTMST)>6.5e-6 ) then
  write(6,*) "xtg out-of-range after xtwepdep"
  write(6,*) "xtg maxval,maxloc ",maxval(xtm1(1:imax,:,:)+xte(1:imax,:,:)*PTMST), &
                                  maxloc(xtm1(1:imax,:,:)+xte(1:imax,:,:)*PTMST)
end if
#endif

!!$acc parallel loop collapse(2) copy(xte(:,:,itracdms:itracso4),so2oh3d,dmsoh3d,dmsn33d) &
!!$acc   copyin(zrdayl,zdayfac,xtm1(:,:,itracdms:itracso2),zzoh,zzno2,zzo3,ptp1,prhop1)
DO CONCURRENT (JK=1:kl)
  DO CONCURRENT (JL=1:imax)
    X=PRHOP1(JL,JK)      
    IF(ZRDAYL(JL)==1) THEN
      !   DAY-TIME CHEMISTRY        
      ZXTP1SO2=XTM1(JL,JK,ITRACSO2)+XTE(JL,JK,ITRACSO2)*PTMST
      ZTK2=ZK2*(PTP1(JL,JK)/300.)**(-3.3)
      ZM=X*ZNAMAIR
      ZHIL=ZTK2*ZM/ZK2I
      ZEXP=ALOG10(ZHIL)
      ZEXP=1./(1.+ZEXP*ZEXP)
      ZTK23B=ZTK2*ZM/(1.+ZHIL)*ZK2F**ZEXP
      ZSO2=ZXTP1SO2*ZZOH(JL,JK)*ZTK23B*ZDAYFAC(jl)
      ZSO2=AMIN1(ZSO2,ZXTP1SO2*PQTMST)
      ZSO2=AMAX1(ZSO2,0.)
      XTE(JL,JK,ITRACSO2)=XTE(JL,JK,ITRACSO2)-ZSO2
      XTE(JL,JK,ITRACSO4)=XTE(JL,JK,ITRACSO4)+ZSO2
      so2oh3d(jl,jk)=zso2

      ZXTP1DMS=XTM1(JL,JK,ITRACDMS)+XTE(JL,JK,ITRACDMS)*PTMST
      T=PTP1(JL,JK)
!     ZTK1=(T*EXP(-234./T)+8.46E-10*EXP(7230./T)+ &
!          2.68E-10*EXP(7810./T))/(1.04E+11*T+88.1*EXP(7460./T)) !Original
      ztk1=1.646e-10-1.850e-12*t+8.151e-15*t**2-1.253e-17*t**3 !Cubic fit good enough
      ztk1=max(ztk1,5.e-12) !Because cubic falls away for T > 300 K
      ! MJT notes - LDR employs 1.5 here, but CCAM's DMS burden is too high
      ! so we revert to the original factor of 2 used for rotstayn and lohmann 2002
      ztk1=2.*ztk1          !This is the fudge factor to account for other oxidants
      !ztk1=1.5*ztk1        !This is the fudge factor to account for other oxidants
      ZDMS=ZXTP1DMS*ZZOH(JL,JK)*ZTK1*ZDAYFAC(jl)
      ZDMS=AMIN1(ZDMS,ZXTP1DMS*PQTMST)
      XTE(JL,JK,ITRACDMS)=XTE(JL,JK,ITRACDMS)-ZDMS
      XTE(JL,JK,ITRACSO2)=XTE(JL,JK,ITRACSO2)+ZDMS
      dmsoh3d(jl,jk)=zdms
    ELSE

      !   NIGHT-TIME CHEMISTRY
      ZXTP1DMS=XTM1(JL,JK,ITRACDMS)+XTE(JL,JK,ITRACDMS)*PTMST
      ZTK3=ZK3*EXP(520./PTP1(JL,JK))
      !    CALCULATE THE STEADY STATE CONCENTRATION OF NO3
      ZQT=1./PTP1(JL,JK)
      ZQT3=300.*ZQT
      ZRHOAIR=PRHOP1(JL,JK)*ZNAMAIR
      ZKNO2O3=1.2E-13*EXP(-2450.*ZQT)
      ZKN2O5AQ=0.1E-04
      ZRX1=2.2E-30*ZQT3**3.9*ZRHOAIR
      !ZRX2=1.5E-12*ZQT3**0.7
      ZRX12=1.467e-18*ZQT3**3.2*ZRHOAIR !=ZRX1/ZRX2
      ZKNO2NO3=ZRX1/(1.+ZRX12)*0.6**(1./(1.+(ALOG10(ZRX12))**2))
      !ZEQN2O5=4.E-27*EXP(10930.*ZQT)
      !ZKN2O5=ZKNO2NO3/ZEQN2O5
      ZKN2O5=5.5E-4*ZQT3**3.9*ZRHOAIR*EXP(-10930.*ZQT)/(1.+ZRX12)*0.6**(1./(1.+(ALOG10(ZRX12))**2))

      ZNO3=ZKNO2O3*(ZKN2O5+ZKN2O5AQ)*ZZNO2(JL,JK)*ZZO3(JL,JK)
      ZZQ=ZKNO2NO3*ZKN2O5AQ*ZZNO2(JL,JK)+(ZKN2O5+ZKN2O5AQ)*ZTK3*ZXTP1DMS*X*6.022E+20/ZMOLGS
      IF(ZZQ>0.) THEN
        ZNO3=ZNO3/ZZQ
      ELSE
        ZNO3=0.
      ENDIF
      ZDMS=ZXTP1DMS*ZNO3*ZTK3
      ZDMS=AMIN1(ZDMS,ZXTP1DMS*PQTMST)
      XTE(JL,JK,ITRACDMS)=XTE(JL,JK,ITRACDMS)-ZDMS
      XTE(JL,JK,ITRACSO2)=XTE(JL,JK,ITRACSO2)+ZDMS
      dmsn33d(jl,jk)=zdms
    ENDIF
  end do
end do
!!$acc end parallel loop

#ifdef debugaero
if ( maxval(xtm1(1:imax,:,:)+xte(1:imax,:,:)*PTMST)>6.5e-6 ) then
  write(6,*) "xtg out-of-range after day/night chemistry"
  write(6,*) "xtg maxval,maxloc ",maxval(xtm1(1:imax,:,:)+xte(1:imax,:,:)*PTMST), &
                                  maxloc(xtm1(1:imax,:,:)+xte(1:imax,:,:)*PTMST)
end if
#endif


! Calculate tendency of SO2 due to oxidation by OH (diagnostic) and ox. tendencies of DMS
so2oh(:) = so2oh(:) + sum( so2oh3d(:,:)*rhodz(:,:), dim=2 )
dmsoh(:) = dmsoh(:) + sum( dmsoh3d(:,:)*rhodz(:,:), dim=2 )
dmsn3(:) = dmsn3(:) + sum( dmsn33d(:,:)*rhodz(:,:), dim=2 )

RETURN
END subroutine xtchemie

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! xt wetdep

SUBROUTINE XTWETDEP(KTRAC,                                                 &
                    PTMST,                                                 &
                    rhodz,                                                 &
                    PMRATEP, PFPREC,                                       &
                    PCLCOVER, PSOLUB, pmlwc, ptp1,                         &
                    pfsnow,pfsubl,pcfcover,pmaccr,pfmelt,pfstayice,        &
                    pqfsedice,plambs,prscav,prfreeze,pfconv,pclcon,fracc,  & !Inputs
                    PXTP10, PXTP1C, PXTP1CON, xtm1, xte, wd, imax, kl)       !In & Out
!$acc routine vector

!
!   *XTWETDEP* CALCULATES THE WET DEPOSITION OF TRACE GASES OR AEROSOLS
!
!   JOHANN FEICHTER              UNI HAMBURG            08-91
!
!   PURPOSE
!  ---------
!   TO CALCULATE THE WET SCAVENGING OF GASES OR AEROSOLS IN CLOUDS
!
!   INTERFACE
!  -------------
!   THIS ROUTINE IS CALLED FROM *XTCHEM*
!
!  METHOD
!  -------
!
!   NO EXTERNALS
!---------------
!

implicit none

! Argument list
integer, intent(in) :: imax, kl
integer, intent(in) :: KTRAC
real, intent(in) :: PTMST
real, dimension(imax,kl), intent(inout) :: PXTP10   !Tracer m.r. outside liquid-water cloud (clear air/ice cloud)
real, dimension(imax,kl), intent(inout) :: PXTP1C   !Tracer m.r.  inside liquid-water cloud
real, dimension(imax,kl), intent(in) :: PXTP1CON
real, dimension(imax,kl), intent(in) :: rhodz
real, dimension(imax,kl), intent(in) :: PMRATEP
real, dimension(imax,kl), intent(in) :: PFPREC
real, dimension(imax,kl) :: PDEP3D
real, dimension(imax,kl), intent(in) :: PCLCOVER
real, dimension(imax,kl), intent(in) :: PSOLUB
real, dimension(imax,kl), intent(in) :: pmlwc
real, dimension(imax,kl), intent(in) :: ptp1  !temperature
real, dimension(imax,kl), intent(in) :: pfsnow
real, dimension(imax,kl), intent(in) :: pfconv
real, dimension(imax,kl), intent(in) :: pclcon
real, dimension(imax), intent(in) :: fracc    !Convective rain fraction (originially set to 0.1)
real, dimension(imax,kl), intent(in) :: pfsubl
real, dimension(imax,kl), intent(in) :: pcfcover
real, dimension(imax,kl), intent(in) :: pmaccr
real, dimension(imax,kl), intent(in) :: pfmelt
real, dimension(imax,kl), intent(in) :: pfstayice
real, dimension(imax,kl), intent(in) :: pqfsedice
real, dimension(imax,kl), intent(in) :: plambs
real, dimension(imax,kl), intent(in) :: prscav
real, dimension(imax,kl), intent(in) :: prfreeze
real, dimension(imax,kl), intent(in) :: xtm1
real, dimension(imax,kl), intent(inout) :: xte
real, dimension(imax), intent(inout) :: wd

! Local work arrays and variables
real, dimension(imax) :: ZDEPS, ZDEPR
real ZMTOF, ZCLR0, zcollefc
real zilcscav, ziicscav,xdep,plambda,zbcscav,xbcscav,zstay_t,xstay,frc
real zmelt,xmelt,zicscav,xicscav
real xfreeze, zfreeze
real ecols_k, rcoeff_k, zcollefs_k, zcollefr_k
real ZDXTE, zxtp1
logical lmask

integer jk,i
real pqtmst

integer, parameter :: ktop = 2    !Top level for wet deposition (counting from top)
logical, parameter :: assume_convliq = .true. ! assume convective rainfall is liquid

!Below-cloud collection eff. for rain
!real, dimension(naero), parameter :: zcollefr = (/0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.10,0.20,0.40,0.05,0.10/)
real, dimension(naero), parameter :: zcollefr = (/0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.05,0.10,0.20,0.40,0.05,0.10/)
!Below-cloud collection eff. for snow
!real, dimension(naero), parameter :: zcollefs = (/0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.02,0.04,0.08,0.01,0.02/)
real, dimension(naero), parameter :: zcollefs = (/0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.01,0.02,0.04,0.08,0.01,0.02/)
!Retention coeff. on riming
real, dimension(naero), parameter :: Rcoeff = (/1.00,0.62,1.00,0.00,1.00,0.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00/)
! Allow in-cloud scavenging in ice clouds for hydrophobic BC and OC, and dust
real, dimension(naero), parameter :: Ecols = (/0.00,0.00,0.00,0.05,0.00,0.05,0.00,0.05,0.05,0.05,0.05,0.05,0.05/)
! wet deposition coefficients
!!Relative re-evaporation rate
!real, dimension(naero), parameter :: Evfac = (/0.25,1.00,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25/)

! Start code : ----------------------------------------------------------

PQTMST = 1./PTMST
ecols_k = ecols(ktrac)
rcoeff_k = rcoeff(ktrac)
zcollefs_k = zcollefs(ktrac)
zcollefr_k = zcollefr(ktrac)

!!$acc enter data create(zdepr,zdeps,zmtof,zclr0,zcollefc,psolub,pxtp10,pxtp1c,pxtp1con,pdep3d,wd)
!!$acc update device(psolub,pxtp10,pxtp1c,pxtp1con,wd)

!!$acc parallel loop present(zdepr,zdeps)
do concurrent (i=1:imax)
  zdepr(i) = 0.
  zdeps(i) = 0.
end do
!!$acc end parallel loop

!!$acc parallel loop collapse(2) present(pdep3d)
do concurrent (jk=1:kl)
  do concurrent (i=1:imax)
    pdep3d(i,jk) = 0.
  end do
end do
!!$acc end parallel loop

!     BEGIN OF VERTICAL LOOP
do JK = KTOP,kl

  !!$acc parallel loop present(pxtp1c,pxtp10,pclcover,pcfcover,pclcon,pqfsedice,pdep3d,             &
  !!$acc   zdeps,zdepr,pmlwc,psolub,pmaccr,plambs,pfsnow,pfsubl,pfstayice,pfmelt,pmratep,prscav,    &
  !!$acc   prfreeze,pfprec,rhodz)
  do concurrent (i=1:imax)

    !ZCLEAR = 1. - PCLCOVER(i,JK) - pcfcover(i,jk) - pclcon(i,jk)
    ZCLR0 = 1. - PCLCOVER(i,JK) - pclcon(i,jk) !Clear air or ice cloud (applies to pxtp10)
    ZMTOF = rhodz(i,jk)*pqtmst
    PXTP1C(i,JK) = AMAX1( 0., PXTP1C(i,JK) )
    PXTP10(i,JK) = AMAX1( 0., PXTP10(i,JK) )

    ! In-cloud ice scavenging (including vertical redistribution when snow
    ! evaporates or falls into a layer). Include accretion of ql by snow.
    if ( zclr0>zmin ) then
      ziicscav = Ecols_k*pqfsedice(i,jk) !qfsedice is the fractional sedimentation in dt
      ziicscav = max( min( ziicscav, 1. ), 0. )
      xdep = pxtp10(i,jk)*ziicscav
      pdep3d(i,jk) = pdep3d(i,jk) + xdep*pcfcover(i,jk)
      !pxtp10(i,jk) = pxtp10(i,jk)*(zclear+(1.-ziicscav)*pcfcover(:,jk))/(1.-pclcover(:,jk))
      pxtp10(i,jk) = pxtp10(i,jk) - xdep*pcfcover(i,jk)/zclr0 ! MJT suggestion
      zdeps(i) = zdeps(i) + xdep*pcfcover(i,jk)*zmtof
    end if

    ! This loop does riming (accretion of liquid water by falling snow)
    if ( pmlwc(i,jk)>zmin ) then
      zilcscav = Rcoeff_k*psolub(i,jk)*pmaccr(i,jk)*ptmst/pmlwc(i,jk)
      zilcscav = max( min( zilcscav, 1. ), 0. )
      xdep = pxtp1c(i,jk)*zilcscav
      pdep3d(i,jk) = pdep3d(i,jk) + xdep*pclcover(i,jk)
      pxtp1c(i,jk) = pxtp1c(i,jk) - xdep
      zdeps(i) = zdeps(i) + xdep*pclcover(i,jk)*zmtof
    end if

    ! Below-cloud scavenging by snow
    !plambs(:,jk) = 1.6e3*10**(-0.023*(ttg(1:imax,k)-tfrz)) ! for ice
    plambda = min( plambs(i,jk), 8.e3 ) !Cut it off at about -30 deg. C
    zbcscav = zcollefs_k*plambda*pfsnow(i,jk)*ptmst/(2.*rhos)
    zbcscav = max( min( 1., zbcscav/(1.+0.5*zbcscav) ), 0. ) !Time-centred
    xbcscav = zbcscav*pxtp10(i,jk)
    pdep3d(i,jk) = pdep3d(i,jk) + xbcscav*zclr0
    pxtp10(i,jk) = pxtp10(i,jk) - xbcscav
    zdeps(i) = zdeps(i) + xbcscav*zclr0*zmtof

    ! Redistribution by snow that evaporates or stays in layer
    if ( pfsnow(i,jk)>zmin .and. zclr0>zmin ) then
      zstay_t = (pfsubl(i,jk)+pfstayice(i,jk))/pfsnow(i,jk)
      zstay_t = max( min( 1., zstay_t ), 0. )
      xstay = zdeps(i)*zstay_t/zmtof
      pdep3d(i,jk) = pdep3d(i,jk) - xstay
      pxtp10(i,jk) = pxtp10(i,jk) + xstay/zclr0
      zdeps(i) = zdeps(i) - xstay*zmtof
      zdeps(i) = max( 0., zdeps(i) )
    end if

    ! Melting of snow... 
    zmelt = pfmelt(i,jk)/max(pfsnow(i,jk)+pfmelt(i,jk),zmin) 
    zmelt = max( min( 1., zmelt ), 0. )
    xmelt = zmelt*zdeps(i)
    zdepr(i) = zdepr(i) + xmelt
    zdeps(i) = zdeps(i) - xmelt
    zdeps(i) = max( 0., zdeps(i) )
  
    !  In-cloud scavenging by warm-rain processes (autoconversion and collection)
    if ( pmlwc(i,jk)>zmin ) then ! MJT suggestion
      zicscav = psolub(i,jk)*pmratep(i,jk)*ptmst/pmlwc(i,jk)
      zicscav = max( min( zicscav, 1. ), 0. )
      xicscav = pxtp1c(i,jk)*zicscav
      pdep3d(i,jk) = pdep3d(i,jk) + xicscav*pclcover(i,jk)
      pxtp1c(i,jk) = pxtp1c(i,jk) - xicscav
      zdepr(i) = zdepr(i) + xicscav*pclcover(i,jk)*zmtof
    end if
 
    ! Below-cloud scavenging by stratiform rain (conv done below)
    zbcscav = zcollefr_k*prscav(i,jk)
    zbcscav = max( min( 1., zbcscav/(1.+0.5*zbcscav) ), 0. ) !Time-centred
    xbcscav = zbcscav*pxtp10(i,jk)
    pdep3d(i,jk) = pdep3d(i,jk) + xbcscav*zclr0
    pxtp10(i,jk) = pxtp10(i,jk) - xbcscav 
    zdepr(i) = zdepr(i) + xbcscav*zclr0*zmtof
  
    ! Freezing of rain... 
    zfreeze = prfreeze(i,jk)/max(pfprec(i,jk)+prfreeze(i,jk),zmin) 
    zfreeze = max( min( 1., zfreeze ), 0. )
    xfreeze = zfreeze*zdepr(i)
    zdeps(i) = zdeps(i) + xfreeze
    zdepr(i) = zdepr(i) - xfreeze
    zdepr(i) = max( 0., zdepr(i) )

  end do
  !!$acc end parallel loop
  
end do !   END OF VERTICAL LOOP

! Now do the convective below-cloud bit...
! In-cloud convective bit was done in convjlm.

! Search for convective cloud base
!kbase(:) = kl+1
!do jk = ktop,kl
!  where ( pclcon(:,jk)>zmin )
!    kbase(:) = k
!  end where
!enddo

do jk = ktop,kl
  !!$acc parallel loop present(rhodz,pclcover,pclcon,ptp1,zcollefc,pfconv,fracc,pxtp10,pdep3d)
  do concurrent (i=1:imax)
    zmtof = rhodz(i,jk)*pqtmst
    zclr0 = 1. - pclcover(i,jk) - pclcon(i,jk)

    ! Use collection efficiencies for rain below melting level, snow above

    ! MJT notes - Assume rain for JLM convection
    if ( ptp1(i,jk)>273.15 .or. assume_convliq ) then
      zcollefc = zcollefr_k
    else
      zcollefc = zcollefs_k  
    end if

    ! Below-cloud scavenging by convective precipitation
    if ( fracc(i)>zmin ) then
      Frc = max( 0., pfconv(i,jk-1)/fracc(i) )
      zbcscav = zcollefc*fracc(i)*0.24*ptmst*sqrt(Frc*sqrt(Frc))
      !zbcscav = min( 1., zbcscav/(1.+0.5*zbcscav) ) !Time-centred
      zbcscav = max( min( 1., zbcscav ), 0. ) ! MJT suggestion
      xbcscav = zbcscav*pxtp10(i,jk)
      pdep3d(i,jk) = pdep3d(i,jk) + xbcscav*zclr0
      pxtp10(i,jk) = pxtp10(i,jk) - xbcscav
      !conwd(i,ktrac) = conwd(i,ktrac) + xbcscav*zclr0*zmtof
    end if

    ! Below-cloud reevaporation of convective rain
    ! This never triggers for JLM convection because pcevap=0.
    ! lmask(:) = jk>kbase(:) .and. pfconv(:,jk-1)>zmin .and. zclr0(:)>zmin
    ! where ( lmask(:) )
    !   pcevap = pfconv(:,jk-1) - pfconv(:,jk)
    !   zevap = pcevap/pfconv(:,jk-1)
    ! elsewhere
    !   zevap(:)=0.
    ! end where
    ! where ( lmask(:) .and. zevap<1. )
    !   zevap = Evfac(ktrac)*zevap
    ! end where
    ! where ( lmask(:) )
    !   zevap = max( 0., min( 1., zevap ) )
    !   xevap = conwd(:,ktrac)*zevap/zmtof(:) !xevap is the grid-box-mean m.r. change
    !   pdep3d(:,jk) = pdep3d(:,jk) - xevap
    !   pxtp10(:,jk) = pxtp10(:,jk) + xevap/zclr0(:)
    !   conwd(:,ktrac) = conwd(:,ktrac) - xevap*zmtof(:)
    !   conwd(:,ktrac) = max( 0., conwd(:,ktrac) )
    ! end where

  end do
  !!$acc end parallel loop
  
end do

!!$acc parallel loop collapse(2) copy(xte) copyin(xtm1) present(pxtp10,pxtp1c,pxtp1con,pclcover,pclcon)
do concurrent (JK=KTOP:kl)
  do concurrent (i=1:imax)
    ZXTP1=(1.-pclcover(i,jk)-pclcon(i,jk))*PXTP10(i,JK)+ &
              PCLCOVER(i,JK)*PXTP1C(i,JK)+               &
              pclcon(i,jk)*pxtp1con(i,jk)
    zxtp1=max(zxtp1,0.)
    ZDXTE=(ZXTP1-XTM1(i,JK))*PQTMST  !Total tendency (Dep + chem)
    !    CHANGE THE TOTAL TENDENCIES
    xte(i,jk) = xte(i,jk) + zdxte
  end do
end do
!!$acc end parallel loop

do concurrent (jk=1:kl)
  !!$acc parallel loop present(pdep3d,rhodz,wd)
  do concurrent (i=1:imax)
    wd(i) = wd(i) + pqtmst*pdep3d(i,jk)*rhodz(i,jk)
  end do
  !!$acc end parallel loop
end do

!!$acc update device(wd)
!!$acc exit data delete(zdepr,zdeps,zmtof,zclr0,zcollefc,psolub,pxtp10,pxtp1c,pxtp1con,pdep3d,wd)

return
end subroutine xtwetdep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Dust settling

subroutine dsettling(tdt,rhoa,tmp,delz,prf,xtg,dustden,dustreff,imax,kl)
!$acc routine vector

implicit none

!     Inputs:
integer, intent(in) :: imax, kl
real, intent(in) :: tdt                  !timestep (s)
real, dimension(imax,kl), intent(in) :: rhoa !air density (kg/m3)
real, dimension(imax,kl), intent(in) :: tmp  !temperature (K)
real, dimension(imax,kl), intent(in) :: delz !Layer thickness (m)
real, dimension(imax,kl), intent(in) :: prf  !Pressure (hPa)
real, dimension(imax,kl,naero), intent(inout) :: xtg
real, dimension(ndust), intent(in) :: dustden, dustreff

! Local work arrays and variables
real, dimension(imax) :: c_stokes, corr, c_cun
real, dimension(imax) :: newxtg, b, dfall
real, dimension(imax) :: vd_cor
integer nt,k

! Start code : ----------------------------------------------------------

do nt = 1, NDUST
  ! Settling velocity (m/s) for each soil classes (Stokes Law)
  ! DUSTDEN     soil class density             (kg/m3)
  ! DUSTREFF    effective radius of soil class (m)
  ! grav        gravity                        (m/s2)
  ! 0.5         upper limit with temp correction (already incorporated with dzmin_gbl - MJT)
  
  ! Solve at the model top
  ! Dynamic viscosity
  C_Stokes = 1.458E-6*TMP(1:imax,kl)**1.5/(TMP(1:imax,kl)+110.4) 
  ! Cuningham correction
  Corr = 6.6E-8*prf(:,kl)/1013.*TMP(1:imax,kl)/293.15
  C_Cun = 1. + 1.249*corr/dustreff(nt)
  ! Settling velocity
  Vd_cor(:) = 2./9.*grav*dustden(nt)*dustreff(nt)**2/C_Stokes*C_Cun
  
  ! Update mixing ratio
  b = tdt*VD_cor(:)/DELZ(:,kl)
  newxtg = xtg(1:imax,kl,nt+itracdu-1)*exp(-b)
  newxtg = max( newxtg, 0. )
  dfall = max( xtg(1:imax,kl,nt+itracdu-1) - newxtg, 0. )
  xtg(1:imax,kl,nt+itracdu-1) = newxtg
  
  ! Solve each vertical layer successively (layer k)
  do k = kl-1,1,-1
    ! Dynamic viscosity
    C_Stokes = 1.458E-6*TMP(1:imax,k)**1.5/(TMP(1:imax,k)+110.4) 
    ! Cuningham correction
    Corr = 6.6E-8*prf(:,k)/1013.*TMP(1:imax,k)/293.15
    C_Cun = 1. + 1.249*corr/dustreff(nt)
    ! Settling velocity
    Vd_cor(:) = 2./9.*grav*dustden(nt)*dustreff(nt)**2/C_Stokes*C_Cun
      
    ! Update mixing ratio
    b = tdt*Vd_cor(:)/DELZ(:,k)
    dfall = dfall * delz(:,k+1)*rhoa(:,k+1)/(delz(:,k)*rhoa(:,k))
    ! Fout = 1.-exp(-b)
    ! Fthru = 1.-Fout/b
    newxtg = xtg(1:imax,k,nt+itracdu-1)*exp(-b) + dfall*(1.-exp(-b))/b
    newxtg = max( newxtg, 0. )
    dfall = max( xtg(1:imax,k,nt+itracdu-1) + dfall - newxtg, 0. )
    xtg(1:imax,k,nt+itracdu-1) = newxtg
  end do
  
end do

return
end subroutine dsettling

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Dust emissions

subroutine dustem(tdt,rhoa,wg,w10m,dz1,vt,snowd,erod,duste,xtg, &
                  dustden,dustreff,imax,kl)
!$acc routine vector

implicit none

!     Inputs:
integer, intent(in) :: imax, kl
real, intent(in) :: tdt                         !Leapfrog timestep (s) (substep and long step)
real, dimension(imax), intent(in) :: rhoa       !air density (kg/m3)
real, dimension(imax), intent(in) :: wg         !ground wetness (fraction of field capacity)
real, dimension(imax), intent(in) :: w10m       !10m windspeed (m/s)
real, dimension(imax), intent(in) :: dz1        !Lowest layer thickness (m)
real, dimension(imax), intent(in) :: vt         !Transfer velocity at surface for dry deposition (m/s)
real, dimension(imax), intent(in) :: snowd      !Snow depth (mm equivalent water)
real, dimension(imax) :: snowa     !Estimated snow areal coverage
real, dimension(imax) :: u_ts0,u_ts,veff
real, dimension(imax) :: srce,dsrc,airmas
real, dimension(imax) :: a,b
real, dimension(imax) :: airden
real, dimension(imax,ndcls), intent(in) :: erod
real, dimension(imax,ndust), intent(inout) :: duste
real, dimension(imax,kl,naero), intent(inout) :: xtg
real, dimension(ndust), intent(in) :: dustden, dustreff
real g,den,diam
integer n,m

real, parameter :: w_dust = 15. ! maxmium wind speed for dust emissions

! This frac_s array gives fraction of source in each size bin.
! Does not quite add to 1, because larger sizes (omitted) account for some too.
! All source is in first four bins, even when using eight bins, since next four are hydrophilic.
real, dimension(ndust), parameter :: frac_s = (/ 0.1, 0.25, 0.25, 0.25 /)
integer, dimension(ndust), parameter :: ipoint = (/ 3, 2, 2, 2 /)                  ! Pointer used for dust classes
                                                                                   ! (sand=1, silt=2, clay=3)

! Start code : ----------------------------------------------------------

g = grav*1.e2
airden = rhoa*1.e-3

! Convert snow depth to estimated areal coverage (see Zender et al 2003, JGR 108, D14, 4416)
! Must convert from mm to m, then adjust by rho_l/rho_s=10.
! 0.05 m is the geometrical snow thickness for 100% areal coverage.
!hsnow = snowd*0.01 !Geometrical snow thickness in metres
snowa = min( 1., snowd/5. )
airmas = dz1 * rhoa  ! kg/m2

do n = 1, ndust
  ! Threshold velocity as a function of the dust density and the diameter
  den = dustden(n)*1.e-3
  diam = 2.*dustreff(n)*1.e2
  ! Pointer to the 3 classes considered in the source data files
  m = ipoint(n)
  
  ! Following is from Ginoux et al (2004) Env. Modelling & Software.
  u_ts0 = 0.13*1.e-2*sqrt(den*g*diam/airden)*sqrt(1.+0.006/den/g/diam**2.5)/ &
          sqrt(1.928*(1331.*diam**1.56+0.38)**0.092-1.)
  
  ! Case of surface dry enough to erode
  where ( wg<0.1 )
    !Tuning suggested for Asian source by P. Ginoux
    u_ts = max( 0., u_ts0*(1.2+0.2*alog10(max( 1.e-3, wg ))) )
  elsewhere
    ! Case of wet surface, no erosion
    u_ts = 100.
  end where
  
  ! MJT notes - erod should be zero for ocean points
    
  !srce = frac_s(n)*erod(i,m)*dxy(i) ! (m2)
  srce = frac_s(n)*erod(:,m) ! (fraction) - MJT suggestion
  where ( w10m < w_dust )
    dsrc = (1.-snowa)*Ch_dust*srce*W10m*W10m*(W10m-u_ts) ! (kg/s/m2)
  elsewhere
    ! limit maximum wind speed to w_dust m/s for emissions - MJT sugestion  
    dsrc = (1.-snowa)*Ch_dust*srce*w_dust*w_dust*(w_dust-u_ts) ! (kg/s/m2)  
  end where
  dsrc = max( 0., dsrc )

  ! Calculate dust mixing ratio tendency at first model level.
  a = dsrc / max(airmas,0.1)
  duste(:,n) = duste(:,n) + dsrc ! Diagnostic
      
  ! Calculate turbulent dry deposition at surface
  ! Use full layer thickness for CSIRO model (should be correct if Vt is relative to mid-layer)
  veff = max( Vt*(wg+(1.-wg)*exp(-max( 0., w10m-u_ts0 ))), 0. )
  b = Veff / dz1

  ! Update mixing ratio
  ! Write in form dx/dt = a - b*x (a = source term, b = drydep term)
  ! solution is x = a/b + (X0-a/b)*exp(-b*tdt).  However, in split form
  ! x = X0 + a*tdt, and x = X0*exp(-b*tdt), or combined
  ! x = (X0 + a*tdt)*exp(-b*tdt)
  xtg(1:imax,1,n+itracdu-1) = (xtg(1:imax,1,n+itracdu-1)+a*tdt)*exp(-min(b*tdt,50.))
  xtg(1:imax,1,n+itracdu-1) = max( 0., xtg(1:imax,1,n+itracdu-1) )

end do

return
end subroutine dustem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! dustage

!subroutine dustage(tdt,rhg) !Inputs
!
!implicit none
!
!!     Inputs:
!real tdt              !Leapfrog timestep
!real rhg(imax,kl)    !RH (fraction 0 to 1)
!
!! Local work arrays and variables
!real rrate(imax,kl)
!real vvso2,rk,xd,dx
!integer k,mg,nt
!
!! Start code : ----------------------------------------------------------
!
!! Reference for this scheme is
!! Fan et al. (2004): Impact of air pollution on wet deposition of mineral dust aerosols, GRL 31,
!! L02104, doi:10.1029/2003GL018501.
!
!! Loop over 4 dust size bins
!! Reaction rate is a linear function of SO2 concentration
!
!do k=1,kl
!  do mg=1,imax
!    vvso2 = xtg(mg,k,itracso2)*29./64. !Vol. mixing ratio of SO2
!    rk = 0.01 * max (0.1, dim(rhg(mg,k),0.5))
!    xd = 0.
!    do nt = itracdu, itracdu+ndust-1
!      xd = xd + xtg(mg,k,nt)
!    enddo
!    rrate(mg,k) = rk * vvso2 / max (1.e-12, xd)
!  enddo
!enddo
!
!! Reduce hydrophobic dust and increase hydrophilic dust
!do nt = itracdu, itracdu+ndust-1
!  do k=1,kl
!    do mg=1,imax
!      dx = min (xtg(mg,k,nt), rrate(mg,k)*tdt*xtg(mg,k,nt))
!      xtg(mg,k,nt) = xtg(mg,k,nt) - dx
!      xtg(mg,k,nt+ndust) = xtg(mg,k,nt+ndust) + dx
!    enddo
!  enddo
!enddo
!
!return
!end subroutine dustage

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! A simple diagnostic treatment of seasalt aerosol (LDR 3/02)
!
!subroutine seasalt(land,fracice,zmid,pblh,v10m,ssn,imax,kl) !Inputs
!!$acc routine vector
!
!implicit none
!
!integer, intent(in) :: imax, kl
!logical, dimension(imax), intent(in) :: land   !True for land points
!real, dimension(imax), intent(in) :: fracice   !Sea-ice fraction
!real, dimension(imax,kl), intent(in) :: zmid   !Height of full level (m)
!real, dimension(imax), intent(in) :: pblh      !PBL height (m)
!real, dimension(imax), intent(in) :: v10m      !10m windpseed, including effect of sub-grid gustiness (m/s)
!real, dimension(imax,2) :: ssn_base
!real, dimension(imax,kl,2), intent(inout) :: ssn
!integer k
!
!! Calculate number and mass concentration of seasalt within the marine BL.
!! Set seasalt conc. to zero elsewhere.
!! Height of BL taken from ncarpbl scheme, so needs this turned on (although we 
!! set it to 2000m in hvertmx if ncarpbl=F)
!! The first mode is the "film-drop" mode, and the second is the "jet-drop" mode.
!! Smaller number mode radii are given by Nilsson et al. (2001) cf. O'Dowd's.
!!
!! References: O'Dowd et al. (1997) Atmos. Environ. 31, 73-80
!!             Jones et al. (2001) JGR 106, 20293-20310.
!!             Nilsson et al. (2001) JGR 106, 32139-32154.
!
!! Jones et al. give different windspeed relations for v10m < 2, 2 <= v10m <= 17.5,
!! and v10m > 17.5.
!where ( v10m<2. )
!  ssn_base(:,1)=3.856e6*(1.-exp(-0.736*v10m))
!  ssn_base(:,2)=0.671e6*(1.-exp(-1.351*v10m))
!elsewhere ( v10m<17.5 )
!  ssn_base(:,1)=10.**(0.0950*v10m+6.2830) 
!  ssn_base(:,2)=10.**(0.0422*v10m+5.7122) 
!elsewhere
!  ssn_base(:,1)=1.5e8*(1.-97.874*exp(-0.313*v10m))
!  ssn_base(:,2)=3.6e6*(1.-103.926*exp(-0.353*v10m))
!end where
!  
!! Number-to-mass conversion factors are derived from the parameters of the two
!! lognormal modes given by Nilsson, together with rhosalt=2.0e3 kg/m3.
!do k=1,kl
!  where ( .not.land .and. zmid(:,k)<pblh )
!    ssn(:,k,1) = ssn_base(:,1)
!    ssn(:,k,2) = ssn_base(:,2)
!  elsewhere ( .not.land )
!    ssn(:,k,1) = 1.e7*exp(-zmid(:,k)/3000.)
!    ssn(:,k,2) = 1.e6*exp(-zmid(:,k)/1500.)
!  elsewhere
!    ssn(:,k,1) = 0.
!    ssn(:,k,2) = 0.
!  end where
!  ! Reduce over sea ice...
!  ssn(:,k,1) = (1.-fracice)*ssn(:,k,1)
!  ssn(:,k,2) = (1.-fracice)*ssn(:,k,2)
!end do
!
!! These relations give ssm in kg/m3 based on ssn in m^{-3}...
!!ssm(:,:,1)=ssn(:,:,1)/saltsmallmtn
!!ssm(:,:,2)=ssn(:,:,2)/saltlargemtn
!
!return
!end subroutine seasalt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! sea salt settling
subroutine ssettling(tdt,rhoa,tmp,delz,prf,xtg,saltden,saltreff,imax,kl)
!$acc routine vector

implicit none

integer, intent(in) :: imax, kl
real, intent(in) :: tdt                  !timestep (s)
real, dimension(imax,kl), intent(in) :: rhoa !air density (kg/m3)
real, dimension(imax,kl), intent(in) :: tmp  !temperature (K)
real, dimension(imax,kl), intent(in) :: delz !Layer thickness (m)
real, dimension(imax,kl), intent(in) :: prf  !Pressure (hPa)
real, dimension(imax,kl,naero), intent(inout) :: xtg
real, dimension(nsalt), intent(in) :: saltden, saltreff
real, dimension(imax) :: c_stokes, corr, c_cun
real, dimension(imax) :: newxtg, b, dfall
real, dimension(imax) :: vd_cor
integer nt, k

do nt = 1,nsalt
  
  ! Solve at the model top
  ! Dynamic viscosity
  C_Stokes = 1.458E-6*TMP(1:imax,kl)**1.5/(TMP(1:imax,kl)+110.4) 
  ! Cuningham correction
  Corr = 6.6E-8*prf(:,kl)/1013.*TMP(1:imax,kl)/293.15
  C_Cun = 1. + 1.249*corr/saltreff(nt)
  ! Settling velocity
  Vd_cor(:) = 2./9.*grav*saltden(nt)*saltreff(nt)**2/C_Stokes*C_Cun
  
  ! Update mixing ratio
  b = tdt*VD_cor(:)/DELZ(:,kl)
  newxtg = xtg(1:imax,kl,nt+itracsa-1)*exp(-b)
  newxtg = max( newxtg, 0. )
  dfall = max( xtg(1:imax,kl,nt+itracsa-1) - newxtg, 0. )
  xtg(1:imax,kl,nt+itracsa-1) = newxtg
  
  ! Solve each vertical layer successively (layer k)
  do k = kl-1,1,-1
    ! Dynamic viscosity
    C_Stokes = 1.458E-6*TMP(1:imax,k)**1.5/(TMP(1:imax,k)+110.4) 
    ! Cuningham correction
    Corr = 6.6E-8*prf(:,k)/1013.*TMP(1:imax,k)/293.15
    C_Cun = 1. + 1.249*corr/saltreff(nt)
    ! Settling velocity
    Vd_cor(:) = 2./9.*grav*saltden(nt)*saltreff(nt)**2/C_Stokes*C_Cun
      
    ! Update mixing ratio
    b = tdt*Vd_cor(:)/DELZ(:,k)
    dfall = dfall * delz(:,k+1)*rhoa(:,k+1)/(delz(:,k)*rhoa(:,k))
    ! Fout = 1.-exp(-b)
    ! Fthru = 1.-Fout/b
    newxtg = xtg(1:imax,k,nt+itracsa-1)*exp(-b) + dfall*(1.-exp(-b))/b
    newxtg = max( newxtg, 0. )
    dfall = max( xtg(1:imax,k,nt+itracsa-1) + dfall - newxtg, 0. )
    xtg(1:imax,k,nt+itracsa-1) = newxtg
  end do
  
end do

return
end subroutine ssettling

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! sea salt emissions
subroutine seasaltem(tdt,v10m,vt,rhoa,dz1,salte,xtg,saltreff,locean,imax,kl)
!$acc routine vector

implicit none

integer, intent(in) :: imax, kl
integer n
real, intent(in) :: tdt
real, dimension(imax), intent(in) :: v10m    ! 10m wind speed
real, dimension(imax), intent(in) :: vt      ! transfer velocity
real, dimension(imax), intent(in) :: dz1     ! layer thickness
real, dimension(imax), intent(in) :: rhoa    ! air density
real, dimension(imax), intent(inout) :: salte
real, dimension(imax,kl,naero), intent(inout) :: xtg
real, dimension(imax) :: wu10, dfdd, df0dd, a, b, veff, diam
real, dimension(imax) :: ftw, fsw
real, dimension(nsalt) :: mtnfactor
real, dimension(nsalt), intent(in) :: saltreff
real, dimension(nsalt), parameter :: saltrange = (/ 0.4e-6, 3.5e-6 /) ! 0.1-0.5um and 0.5-4um
logical, dimension(imax), intent(in) :: locean

! Follows Sofiev et al 2011 for emissions
wu10 = 3.84e-6*v10m**3.41

mtnfactor(1) = saltsmallmtn
mtnfactor(2) = saltlargemtn

do n = 1,nsalt
  ! Calculate emissions  
  diam = 2.*saltreff(n)*1.e6
    
  !df0dr = 3.6e5*(1.+0.057*diam**1.05)/(diam**3)*10**(1.19*exp(-((0.38-log(diam))/0.65)**2))
  
  ! Follows Sofiev et al 2011 for emissions (accounts for <0.1 microns)
  df0dd = 1.e6*exp(-0.09/(diam+3.e-3))/(2.+exp(-5./diam))*(1.+0.05*diam**1.05)/(diam**3) &
         *10**(1.05*exp(-((0.27-log(diam))/1.1)**2))
    
  ftw = 0.48*diam**(-0.36)   ! 15C
  !ftw = 0.15*diam**(-0.88)  ! 5C
  !ftw = 0.092*diam**(-0.96) ! -2C
  
  fsw = 0.12*diam**(-0.71)
  !fsw = 5.85*e-5*diam**-1.7 ! sal = 0.
  
  dfdd = wu10*df0dd*ftw*fsw  ! number/micro/m^2/s
  
  where ( locean )
    a = dfdd*2.*saltrange(n)*1.e6/(rhoa*dz1) ! number/kg/s
  elsewhere
    a = 0.
  end where 
  
  a = a/mtnfactor(n) ! kg/kg/s  
  
  salte = salte + a*rhoa*dz1 ! Diagnostic
  
  ! Calculate turbulent dry deposition at surface
  veff = Vt
  b = Veff / dz1

  ! Update mixing ratio
  ! Write in form dx/dt = a - b*x (a = source term, b = drydep term)
  ! solution is x = a/b + (X0-a/b)*exp(-b*tdt).  However, in split form
  ! x = X0 + a*tdt, and x = X0*exp(-b*tdt), or combined
  ! x = (X0 + a*tdt)*exp(-b*tdt)
  xtg(1:imax,1,n+itracsa-1) = (xtg(1:imax,1,n+itracsa-1)+a*tdt)*exp(-min(b*tdt,50.))
  xtg(1:imax,1,n+itracsa-1) = max( 0., xtg(1:imax,1,n+itracsa-1) )
end do    

return
end subroutine seasaltem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! cloud droplet concentration

subroutine cldrop(istart,cdn,rhoa,convmode)

implicit none

integer, intent(in) :: istart
integer k,is,ie,imax,kl
real, dimension(:,:), intent(in) :: rhoa
real, dimension(:,:), intent(out) :: cdn
real, dimension(size(cdn,1),size(cdn,2)) :: xtgso4,xtgbc,xtgoc,xtgsa1,xtgsa2
real, dimension(size(cdn,1)) :: so4_n,cphil_n,salt_n,Atot
real, dimension(size(cdn,1)) :: so4mk
logical, intent(in) :: convmode

imax = size(cdn,1)
kl = size(cdn,2)

is = istart
ie = istart + imax - 1

if ( convmode ) then
  ! total grid-box
  xtgso4 = xtg(is:ie,:,itracso4)
  xtgbc  = xtg(is:ie,:,itracbc+1)
  xtgoc  = xtg(is:ie,:,itracoc+1)
  xtgsa1 = xtg(is:ie,:,itracsa)
  xtgsa2 = xtg(is:ie,:,itracsa+1)
else
  ! outside convective fraction of grid-box
  xtgso4 = xtosav(is:ie,:,itracso4)
  xtgbc  = xtosav(is:ie,:,itracbc+1)
  xtgoc  = xtosav(is:ie,:,itracoc+1)
  xtgsa1 = xtosav(is:ie,:,itracsa)
  xtgsa2 = xtosav(is:ie,:,itracsa+1)
end if

select case(aeroindir)
  case(0)
    do k = 1,kl
      ! Factor of 132.14/32.06 converts from sulfur to ammmonium sulfate
      so4_n = so4mtn * (132.14/32.06) * rhoa(:,k) * xtgso4(:,k)
      ! Factor of 1.3 converts from OC to organic matter (OM) 
      cphil_n = carbmtn * rhoa(:,k) * (xtgbc(:,k)+1.3*xtgoc(:,k))
      salt_n = saltsmallmtn*rhoa(:,k)*xtgsa1(:,k) + saltlargemtn*rhoa(:,k)*xtgsa2(:,k)
      ! Jones et al., modified to account for hydrophilic carb aerosols as well
      Atot = max( so4_n + cphil_n + salt_n, 0. )
      cdn(:,k) = max( 1.e7, 3.75e8*(1.-exp(-2.5e-9*Atot)) )
    end do

  case(1)
    ! Use ECHAM SO4 to get cdn_strat.
    do k = 1,kl
      so4mk = max( 1.e-5, 3.e9*rhoa(:,k)*xtgso4(:,k) ) ! x 3 to convert to ug/m3 SO4
      cdn(:,k) = max( 2.e7, 1.62e8*so4mk**0.41 )       !Combined land/ocean.
    end do
    
  case default
    write(6,*) "ERROR: Invaild aeroindir option ",aeroindir
    stop
    
end select

return
end subroutine cldrop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Aerosol scavenging fraction for convective clouds

pure subroutine convscav(fscav,xpkp1,xpold,tt,xs,rho,ntr,kx)
!$acc routine vector

implicit none

integer, intent(in) :: ntr, kx
real, dimension(kx), intent(out) :: fscav ! scavenging fraction
real, dimension(kx), intent(in) :: xpkp1 ! cloud liquid water after precipitation
real, dimension(kx), intent(in) :: xpold ! cloud liquid water before precipitation
real, dimension(kx), intent(in) :: tt    ! parcel temperature
real, dimension(kx), intent(in) :: xs    ! xtg(:,k,3) = so4
real, dimension(kx), intent(in) :: rho   ! air density
real, dimension(kx) :: f_so2,scav_eff
real, dimension(kx) :: zqtp1,ze2,ze3,zfac,zso4l,zso2l,zqhp
real, dimension(kx) :: zza,zzb,zzp,zzq,zzp2,zhp,zheneff,p_so2
logical, dimension(kx) :: bwkp1 

! In-cloud scavenging efficiency for liquid and frozen convective clouds follows.
! Note that value for SO2 (index 2) is overwritten by Henry coefficient f_so2 below.
!real, parameter, dimension(naero) :: scav_effl = (/0.00,1.00,0.90,0.00,0.30,0.00,0.30,0.05,0.05,0.05,0.05,0.05,0.05/) ! liquid
real, parameter, dimension(naero) :: scav_effl = (/0.0,1.0,0.9,0.0,0.3,0.0,0.3,0.3,0.3,0.3,0.3,0.05,0.05/) ! liquid
!real, parameter, dimension(naero) :: scav_effi = (/0.00,0.00,0.00,0.05,0.00,0.05,0.00,0.05,0.05,0.05,0.05,0.05,0.05/) ! ice

!bwkp1(:) = tt(:)>=ticeu ! condensate in parcel is liquid (true) or ice (false)
bwkp1(:) = .true.        ! assume liquid for JLM convection

if ( ntr==itracso2 ) then
  !where ( bwkp1 )
    ! CALCULATE THE SOLUBILITY OF SO2
    ! TOTAL SULFATE  IS ONLY USED TO CALCULATE THE PH OF CLOUD WATER
    ZQTP1 = 1./tt - 1./298.
    ZE2  =1.23*EXP(3020.*ZQTP1)
    ZE3 = 1.2E-02*EXP(2010.*ZQTP1)
    ZFAC = 1000./(max(xpold,1.E-20)*32.064)
    ZSO4L = xs*ZFAC
    ZSO4L = AMAX1(ZSO4L,0.)
    ZSO2L = xs*ZFAC
    ZSO2L = AMAX1(ZSO2L,0.)
    ZZA = ZE2*8.2E-02*tt*max(xpold,1.E-20)*rho*1.E-03
    ZZB = 2.5E-06+ZSO4L
    ZZP = (ZZA*ZE3-ZZB-ZZA*ZZB)/(1.+ZZA)
    ZZQ = -ZZA*ZE3*(ZZB+ZSO2L)/(1.+ZZA)
    ZZP = 0.5*ZZP
    ZZP2 = ZZP*ZZP
    ZHP = -ZZP + SQRT(max(ZZP2-ZZQ,0.))
    ZQHP = 1./ZHP
    ZHENEFF = 1. + ZE3*ZQHP
    P_SO2 = ZZA*ZHENEFF
    F_SO2 = P_SO2/(1.+P_SO2)
    F_SO2 = min(max(0.,F_SO2),1.)
    scav_eff = f_so2
  !elsewhere
  !  scav_eff = scav_effi(ntr)
  !end where
else
  !where ( bwkp1 )
    scav_eff = scav_effl(ntr)
  !elsewhere
  !  scav_eff = scav_effi(ntr)
  !end where
end if

! Wet deposition scavenging fraction
fscav(:) = scav_eff(:)*min(max(xpold(:)-xpkp1(:),0.)/max(xpold(:),1.E-20),1.)
fscav(:) = min( fscav(:), 1. )

return
end subroutine convscav

end module aerosolldr
