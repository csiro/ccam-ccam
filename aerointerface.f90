! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2019 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
! This module interfaces aerosol schemes with CCAM

! Currently, only LDR aerosols are supported.  However, the xtg arrays should be moved to
! here in the future, so that there is a common interface for advection and outcdf.

! - This subroutine assumes only one month at a time is integrated in RCM mode.

module aerointerface

implicit none

private
public load_aerosolldr, aerocalc, aerodrop
public ppfprec, ppfmelt, ppfsnow, ppfevap, ppfsubl, pplambs, ppmrate
public ppmaccr, ppfstayice, ppfstayliq, ppqfsedice, pprscav, pprfreeze
public opticaldepth, updateoxidant, oxidant_timer

integer, save :: ilon, ilat, ilev
integer, save :: oxidant_timer = -9999
integer, parameter :: naerofamilies = 7      ! Number of aerosol families for optical depth
integer, parameter :: updateoxidant = 1440   ! update prescribed oxidant fields once per day
real, dimension(:,:,:), allocatable, save :: oxidantprev_g
real, dimension(:,:,:), allocatable, save :: oxidantnow_g
real, dimension(:,:,:), allocatable, save :: oxidantnext_g
real, dimension(:,:,:), allocatable, save :: opticaldepth
real, dimension(:,:), allocatable, save :: ppfprec, ppfmelt, ppfsnow           ! data saved from LDR cloud scheme
real, dimension(:,:), allocatable, save :: ppfevap, ppfsubl, pplambs, ppmrate  ! data saved from LDR cloud scheme
real, dimension(:,:), allocatable, save :: ppmaccr, ppqfsedice, pprscav        ! data saved from LDR cloud scheme
real, dimension(:,:), allocatable, save :: pprfreeze                           ! data saved from LDR cloud scheme
real, dimension(:,:), allocatable, save :: ppfstayice, ppfstayliq              ! data saved from LDR cloud scheme
real, dimension(:), allocatable, save :: rlev
real, dimension(:), allocatable, save :: zdayfac
real, parameter :: wlc = 0.2e-3         ! LWC of deep conv cloud (kg/m**3)

contains

subroutine aerocalc(oxidant_update,mins)

use aerosolldr           ! LDR prognostic aerosols
use arrays_m             ! Atmosphere dyamics prognostic arrays
use cc_omp               ! CC OpenMP routines
use cfrac_m              ! Cloud fraction
use extraout_m           ! Additional diagnostics
use infile               ! Input file routines
use kuocomb_m            ! JLM convection
use latlong_m            ! Lat/lon coordinates
use liqwpar_m            ! Cloud water mixing ratios
use morepbl_m            ! Additional boundary layer diagnostics
use newmpar_m            ! Grid parameters
use nharrs_m             ! Non-hydrostatic atmosphere arrays
use nsibd_m              ! Land-surface arrays
use ozoneread            ! Ozone input routines
use pbl_m                ! Boundary layer arrays
use screen_m             ! Screen level diagnostics
use sigs_m               ! Atmosphere sigma levels
use soil_m               ! Soil and surface data
use soilsnow_m           ! Soil, snow and surface data
use soilv_m              ! Soil parameters
use vegpar_m             ! Vegetation arrays
use work2_m              ! Diagnostic arrays
use zenith_m             ! Astronomy routines

implicit none

integer, intent(in) :: mins
integer :: tile, is, ie
real, dimension(imax,ilev,4) :: loxidantprev, loxidantnow, loxidantnext
real, dimension(imax,kl,naero) :: lxtg, lxtosav, lxtg_solub
real, dimension(imax,kl,4) :: lzoxidant
real, dimension(imax,kl,2) :: lssn
real, dimension(imax,kl) :: lt, lqg, lqlg, lqfg, lstratcloud
real, dimension(imax,kl) :: lppfprec, lppfmelt, lppfsnow, lppfevap, lppfsubl, lpplambs
real, dimension(imax,kl) :: lppmrate, lppmaccr, lppfstayice, lppfstayliq, lppqfsedice
real, dimension(imax,kl) :: lpprscav, lpprfreeze
real, dimension(imax,ndust) :: lduste, ldustdd, ldust_burden, ldustwd
real, dimension(imax,ndcls) :: lerod
real, dimension(imax,15) :: lemissfield
logical, intent(in) :: oxidant_update

!$omp do schedule(static) private(is,ie),                                                              &
!$omp private(loxidantprev,loxidantnow,loxidantnext,lt,lqg,lqlg,lqfg),                                 &
!$omp private(lstratcloud,lppfprec,lppfmelt,lppfsnow,lppfevap,lppfsubl,lpplambs,lppmrate,lppmaccr),    &
!$omp private(lppfstayice,lppfstayliq,lppqfsedice,lpprscav,lpprfreeze,lxtg,lzoxidant,lduste,ldustdd),  &
!$omp private(lxtosav,lxtg_solub,ldust_burden,lerod,lssn,ldustwd,lemissfield)
do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  
  loxidantprev(1:imax,1:ilev,1:4) = oxidantprev_g(is:ie,1:ilev,1:4)
  loxidantnow(1:imax,1:ilev,1:4)  = oxidantnow_g(is:ie,1:ilev,1:4)
  loxidantnext(1:imax,1:ilev,1:4) = oxidantnext_g(is:ie,1:ilev,1:4)
  lzoxidant(1:imax,1:kl,1:4)      = zoxidant_g(is:ie,1:kl,1:4)
  lxtg                            = xtg(is:ie,:,:)
  lxtosav                         = xtosav(is:ie,:,:)
  lssn                            = ssn(is:ie,:,:)
  lduste                          = duste(is:ie,:)
  ldustdd                         = dustdd(is:ie,:)
  ldustwd                         = dustwd(is:ie,:)
  ldust_burden                    = dust_burden(is:ie,:)
  lemissfield                     = emissfield(is:ie,:)
  lerod                           = erod(is:ie,:)
  lt                              = t(is:ie,:)
  lqg                             = qg(is:ie,:)
  lqlg                            = qlg(is:ie,:)
  lqfg                            = qfg(is:ie,:)
  lstratcloud                     = stratcloud(is:ie,:)
  lppfprec                        = ppfprec(is:ie,:)
  lppfmelt                        = ppfmelt(is:ie,:)
  lppfsnow                        = ppfsnow(is:ie,:)
  lppfevap                        = ppfevap(is:ie,:)
  lppfsubl                        = ppfsubl(is:ie,:)
  lpplambs                        = pplambs(is:ie,:)
  lppmrate                        = ppmrate(is:ie,:)
  lppmaccr                        = ppmaccr(is:ie,:)
  lppfstayice                     = ppfstayice(is:ie,:)
  lppfstayliq                     = ppfstayliq(is:ie,:)
  lppqfsedice                     = ppqfsedice(is:ie,:)
  lpprscav                        = pprscav(is:ie,:)
  lpprfreeze                      = pprfreeze(is:ie,:)  

  call aerocalc_work(loxidantprev,loxidantnow,loxidantnext,ps(is:ie),zdayfac(is:ie),rlatt(is:ie),rlongg(is:ie),       &
                     lt,kbsav(is:ie),ktsav(is:ie),                                                                    &
                     wetfac(is:ie),pblh(is:ie),tss(is:ie),condc(is:ie),snowd(is:ie),fg(is:ie),eg(is:ie),              &
                     u10(is:ie),ustar(is:ie),zo(is:ie),land(is:ie),fracice(is:ie),sigmf(is:ie),lqg,lqlg,lqfg,         &
                     lstratcloud,cdtq(is:ie),lppfprec,lppfmelt,lppfsnow,lppfevap,lppfsubl,lpplambs,lppmrate,lppmaccr, &
                     lppfstayice,lppfstayliq,lppqfsedice,lpprscav,lpprfreeze,so4t(is:ie),lxtg,lzoxidant,lduste,       &
                     ldustdd,                                                                                         &
                     lxtosav,lxtg_solub,dmsso2o(is:ie),so2so4o(is:ie),ldust_burden,bc_burden(is:ie),                  &
                     oc_burden(is:ie),dms_burden(is:ie),                                                              &
                     so2_burden(is:ie),so4_burden(is:ie),lerod,lssn,so2wd(is:ie),so4wd(is:ie),bcwd(is:ie),            &
                     ocwd(is:ie),ldustwd,lemissfield,vso2(is:ie),dmse(is:ie),                                         &
                     so2e(is:ie),so4e(is:ie),bce(is:ie),oce(is:ie),so2dd(is:ie),so4dd(is:ie),bcdd(is:ie),             &
                     ocdd(is:ie),mins,oxidant_update)

  zoxidant_g(is:ie,1:kl,1:4) = lzoxidant(1:imax,1:kl,1:4)
  xtg(is:ie,1:kl,1:naero)    = lxtg(1:imax,1:kl,1:naero)
  ssn(is:ie,:,:)             = lssn
  duste(is:ie,:)             = lduste
  dustdd(is:ie,:)            = ldustdd
  dustwd(is:ie,:)            = ldustwd
  dust_burden(is:ie,:)       = ldust_burden
  
end do
!$omp end do nowait

return
end subroutine aerocalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update prognostic aerosols
subroutine aerocalc_work(oxidantprev,oxidantnow,oxidantnext,ps,zdayfac,rlatt,rlongg,t,kbsav,ktsav,          &
                         wetfac,pblh,tss,condc,snowd,fg,eg,u10,ustar,zo,land,fracice,sigmf,qg,qlg,qfg,      &
                         stratcloud,cdtq,ppfprec,ppfmelt,ppfsnow,ppfevap,ppfsubl,pplambs,ppmrate,ppmaccr,   &
                         ppfstayice,ppfstayliq,ppqfsedice,pprscav,pprfreeze,so4t,xtg,zoxidant,duste,dustdd, &
                         xtosav,xtg_solub,dmsso2o,so2so4o,dust_burden,bc_burden,oc_burden,dms_burden,       &
                         so2_burden,so4_burden,erod,ssn,so2wd,so4wd,bcwd,ocwd,dustwd,emissfield,vso2,dmse,  &
                         so2e,so4e,bce,oce,so2dd,so4dd,bcdd,ocdd,mins,oxidant_update)

use aerosolldr, only : naero,ndcls,aldrcalc,ndust                  ! LDR prognostic aerosols
use cc_mpi                                                         ! CC MPI routines
use cc_omp, only : imax, ntiles                                    ! CC OpenMP routines
use cloudmod, only : convectivecloudfrac                           ! Prognostic strat cloud
use const_phys                                                     ! Physical constants
use newmpar_m                                                      ! Grid parameters
use ozoneread, only : fieldinterpolate                             ! Ozone input routines
use parm_m                                                         ! Model configuration
use sigs_m                                                         ! Atmosphere sigma levels
use zenith_m, only : solargh,zenith                                ! Astronomy routines

implicit none

include 'kuocom.h'      ! Convection parameters

integer, intent(in) :: mins
integer smins
integer j,k,tt,ttx,kinv
integer, dimension(imax), intent(in) :: kbsav, ktsav
real dhr,fjd,r1,dlt,alp,slag
real, dimension(imax,ilev,4), intent(in) :: oxidantprev, oxidantnow, oxidantnext
real, dimension(imax,kl,naero), intent(inout) :: xtg, xtg_solub
real, dimension(imax,kl,naero), intent(in) :: xtosav
real, dimension(imax,kl,4), intent(inout) :: zoxidant
real, dimension(imax,kl,2), intent(inout) :: ssn
real, dimension(imax,kl), intent(in) :: t, qg, qlg, qfg, stratcloud
real, dimension(imax,kl), intent(in) :: ppfprec
real, dimension(imax,kl), intent(in) :: ppfmelt
real, dimension(imax,kl), intent(in) :: ppfsnow
real, dimension(imax,kl), intent(in) :: ppfevap
real, dimension(imax,kl), intent(in) :: ppfsubl
real, dimension(imax,kl), intent(in) :: pplambs
real, dimension(imax,kl), intent(in) :: ppmrate
real, dimension(imax,kl), intent(in) :: ppmaccr
real, dimension(imax,kl), intent(in) :: ppfstayice
real, dimension(imax,kl), intent(in) :: ppfstayliq
real, dimension(imax,kl), intent(in) :: ppqfsedice
real, dimension(imax,kl), intent(in) :: pprscav
real, dimension(imax,kl), intent(in) :: pprfreeze
real, dimension(imax,ndust), intent(inout) :: duste, dustdd, dust_burden, dustwd
real, dimension(imax,ndcls), intent(in) :: erod
real, dimension(imax,15), intent(in) :: emissfield
real, dimension(imax), intent(inout) :: zdayfac
real, dimension(imax), intent(in) :: ps, rlatt, rlongg, wetfac, pblh, tss, condc, snowd
real, dimension(imax), intent(in) :: fg, eg, u10, ustar, zo, fracice, sigmf, cdtq
real, dimension(imax), intent(in) :: vso2
real, dimension(imax), intent(inout) :: dmsso2o, so2so4o, bc_burden, oc_burden, dms_burden
real, dimension(imax), intent(inout) :: so2_burden, so4_burden, so2wd, so4wd, bcwd, ocwd
real, dimension(imax), intent(inout) :: dmse, so2e, so4e, bce, oce, so2dd, so4dd, bcdd, ocdd
real, dimension(imax), intent(out) :: so4t
real, dimension(imax,kl) :: zg,clcon,pccw,rhoa
real, dimension(imax,kl) :: dz
real, dimension(imax) :: coszro,taudar
real, dimension(imax) :: cldcon,wg
logical, dimension(imax), intent(in) :: land
logical, intent(in) :: oxidant_update

! update prescribed oxidant fields
dhr = dt/3600.
if ( oxidant_update ) then
  do j = 1,4 
    ! note levels are inverted by fieldinterpolate
    call fieldinterpolate(zoxidant(:,:,j),oxidantprev(:,:,j),oxidantnow(:,:,j),oxidantnext(:,:,j), &
                          rlev,imax,kl,ilev,mins,sig,ps,interpmeth=0,meanmeth=1)
  end do
  ! estimate day length (presumably to preturb day-time OH levels)
  ttx = nint(86400./dt)
  zdayfac(:) = 0.
  do tt = ttx,1,-1 ! we seem to get a different answer if dhr=24. and ttx=1.
    smins = int(real(tt-1)*dt/60.) + mins
    fjd = float(mod( smins, 525600 ))/1440.  ! 525600 = 1440*365
    call solargh(fjd,bpyear,r1,dlt,alp,slag)
    call zenith(fjd,r1,dlt,slag,rlatt,rlongg,dhr,imax,coszro,taudar)
    where ( taudar>0.5 )
      zdayfac(:) = zdayfac(:) + 1.
    end where
  end do
  where ( zdayfac>0.5 )
    zdayfac(:) = real(ttx)/zdayfac(:)
  end where
else
  fjd = float(mod( mins, 525600 ))/1440.  ! 525600 = 1440*365
  call solargh(fjd,bpyear,r1,dlt,alp,slag)
  call zenith(fjd,r1,dlt,slag,rlatt,rlongg,dhr,imax,coszro,taudar)
  ! taudar is for current timestep - used to indicate sunlit
end if

! set-up input data fields ------------------------------------------------

zg(:,1) = bet(1)*t(1:imax,1)/grav
do k = 2,kl
  zg(:,k) = zg(:,k-1) + (bet(k)*t(1:imax,k)+betm(k)*t(1:imax,k-1))/grav ! height above surface in meters
end do
do k = 1,kl
  dz(:,k) = -rdry*dsig(k)*t(1:imax,k)/(grav*sig(k))
  dz(:,k) = min( max( dz(:,k), 1. ), 2.e4 )
  rhoa(:,k) = ps(1:imax)*sig(k)/(rdry*t(1:imax,k)) ! density of air (kg/m**3)
end do

! estimate convective cloud fraction from leoncld.f
call convectivecloudfrac(clcon,kbsav,ktsav,condc,cldcon=cldcon)
pccw(:,:) = 0.
do k = 1,kl
  kinv = kl + 1 - k  
  ! MJT notes - Assume rain for JLM convection
  !where ( k>kbsav .and. k<=ktsav .and. t(1:imax,k)>ticeu )
  !  pccw(:,kl+1-k) = 0.
  where ( k>kbsav(1:imax) .and. k<=ktsav(1:imax) )
    pccw(1:imax,kinv) = wlc/rhoa(1:imax,k)
  end where
end do

! Water converage at surface
wg(:) = min( max( wetfac, 0. ), 1. )

! MJT notes - We have an option to update the aerosols before the vertical mixing
! or after the vertical mixing.  Updating aerosols before the vertical mixing
! ensures that we can split the convective and non-convective aerosol
! concentrations.  However, updating aerosols after vertical mixing provides a
! better estimate of u10 and pblh.

! update prognostic aerosols
call aldrcalc(dt,sig,zg,dz,wg,pblh,ps,tss,                    &
              t,condc,snowd,taudar,fg,eg,u10,ustar,zo,        &
              land,fracice,sigmf,qg,qlg,qfg,stratcloud,clcon, &
              cldcon,pccw,rhoa,cdtq,ppfprec,ppfmelt,          &
              ppfsnow,ppfevap,ppfsubl,pplambs,ppmrate,        &
              ppmaccr,ppfstayice,ppfstayliq,ppqfsedice,       &
              pprscav,pprfreeze,zdayfac,kbsav,xtg,duste,      &
              dustdd,xtosav,xtg_solub,dmsso2o,so2so4o,        &
              dust_burden,bc_burden,oc_burden,dms_burden,     &
              so2_burden,so4_burden,erod,ssn,zoxidant,        &
              so2wd,so4wd,bcwd,ocwd,dustwd,emissfield,        &
              vso2,dmse,so2e,so4e,bce,oce,so2dd,so4dd,        &
              bcdd,ocdd,imax)
              

! store sulfate for LH+SF radiation scheme.  SEA-ESF radiation scheme imports prognostic aerosols in seaesfrad.f90.
! Factor 1.e3 to convert to gS/m2, x 3 to get sulfate from sulfur
so4t(:) = 0.
do k = 1,kl
  so4t(:) = so4t(:) + 3.e3*xtg(1:imax,k,3)*rhoa(:,k)*dz(:,k)
enddo

if ( diag .and. mydiag .and. ntiles==1 ) then
  write(6,*) "tdiag ",t(idjd,:)
  write(6,*) "qgdiag ",qg(idjd,:)
  write(6,*) "qlgdiag ",qlg(idjd,:)
  write(6,*) "qfgdiag ",qfg(idjd,:)
  write(6,*) "u10diag ",u10(idjd)
  write(6,*) "pblhdiag ",pblh(idjd)
  write(6,*) "fracicediag ",fracice(idjd)
  write(6,*) "DMSdiag ",xtg(idjd,:,1)
  write(6,*) "SO2diag ",xtg(idjd,:,2)
  write(6,*) "SO4diag ",xtg(idjd,:,3)
  write(6,*) "BCphobdiag ",xtg(idjd,:,4)
  write(6,*) "BCphildiag ",xtg(idjd,:,5)
  write(6,*) "OCphobdiag ",xtg(idjd,:,6)
  write(6,*) "OCphildiag ",xtg(idjd,:,7)
  write(6,*) "dust0.8diag ",xtg(idjd,:,8)
  write(6,*) "dust1.0diag ",xtg(idjd,:,9)
  write(6,*) "dust2.0diag ",xtg(idjd,:,10)
  write(6,*) "dust4.0diag ",xtg(idjd,:,11)
  write(6,*) "saltfilmdiag ",ssn(idjd,:,1)
  write(6,*) "saltjetdiag  ",ssn(idjd,:,2)
end if

return
end subroutine aerocalc_work

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Estimate cloud droplet size
subroutine aerodrop(istart,imax,cdn,rhoa,outconv)

use aerosolldr              ! LDR prognostic aerosols
use const_phys              ! Physical constants
use latlong_m, only : rlatt ! Lat/lon coordinates
use newmpar_m               ! Grid parameters
use parm_m                  ! Model configuration
use soil_m, only : land     ! Soil and surface data

implicit none

integer, intent(in) :: istart,imax
integer k,indirmode,iend
real, dimension(imax,kl), intent(out) :: cdn
real, dimension(imax,kl), intent(in) :: rhoa
real, parameter :: cdrops_nh=1.e8, cdropl_nh=3.e8 !Cloud droplet conc sea/land nh
!real, parameter :: cdrops_sh=1.e8, cdropl_sh=3.e8 !Cloud droplet conc sea/land sh
real, parameter :: cdrops_sh=.5e8, cdropl_sh=1.e8 !Cloud droplet conc sea/land sh
logical, intent(in), optional :: outconv
logical convmode

convmode = .true.
if ( present(outconv) ) then
  convmode = .not.outconv
end if

indirmode = abs(iaero)
if ( aeroindir==2 ) then
  indirmode = 0 ! option for no indirect effects
end if

select case( indirmode )
  case( 2, 3 )
    ! prognostic aerosols for indirect effects
    call cldrop(istart,imax,cdn,rhoa,convmode)
  case default
    ! diagnosed for prescribed aerosol indirect effects
    iend = istart + imax - 1
    where ( land(istart:iend).and.rlatt(istart:iend)>0. )
      cdn(:,1) = cdropl_nh
    elsewhere ( land(istart:iend) )
      cdn(:,1) = cdropl_sh
    elsewhere ( rlatt(istart:iend)>0. )
      cdn(:,1) = cdrops_nh
    elsewhere
      cdn(:,1) = cdrops_sh
    end where
    do k = 2,kl
      cdn(:,k) = cdn(:,1)
    end do
end select

return
end subroutine aerodrop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Load aerosols emissions from netcdf
subroutine load_aerosolldr(aerofile, oxidantfile, kdatein)
      
use aerosolldr          ! LDR prognostic aerosols
use cc_mpi              ! CC MPI routines
use cc_omp              ! CC OpenMP routines
use infile              ! Input file routines
use newmpar_m           ! Grid parameters
use ozoneread           ! Ozone input routines
use parmgeom_m          ! Coordinate data
use sigs_m              ! Atmosphere sigma levels
      
implicit none

integer, intent(in) :: kdatein
integer ncstatus, ncid, i, j, varid, tilg
integer jyear, jmonth
integer premonth, nxtmonth
integer, dimension(2) :: spos, npos
integer, dimension(3) :: idum
integer, dimension(4) :: sposs, nposs
real, dimension(:,:), allocatable, save :: dumg
real, dimension(ifull,16) :: duma
real, dimension(:,:,:,:), allocatable, save :: oxidantdum
real, dimension(:), allocatable, save :: rlon, rlat
real, dimension(:), allocatable, save :: rpack
real tlat, tlon, tschmidt
real, parameter :: iotol = 1.E-5 ! tolarance for iotest
character(len=*), intent(in) :: aerofile, oxidantfile
logical tst

if ( myid==0 ) write(6,*) "Initialise LDR prognostic aerosols"

allocate( ppfprec(ifull,kl), ppfmelt(ifull,kl) )
allocate( ppfsnow(ifull,kl) )
allocate( ppfevap(ifull,kl), ppfsubl(ifull,kl) )
allocate( pplambs(ifull,kl), ppmrate(ifull,kl) )
allocate( ppmaccr(ifull,kl) )
allocate( ppqfsedice(ifull,kl), pprscav(ifull,kl), pprfreeze(ifull,kl) )
allocate( ppfstayice(ifull,kl), ppfstayliq(ifull,kl) )
allocate( zdayfac(ifull) )
allocate( opticaldepth(ifull,naerofamilies,3) )
ppfprec = 0.
ppfmelt = 0.
ppfsnow = 0.
ppfevap = 0.
ppfsubl = 0.
pplambs = 0.
ppmrate = 0.
ppmaccr = 0.
ppfstayice = 0.
ppfstayliq = 0.
ppqfsedice = 0.
pprscav = 0.
pprfreeze = 0.
zdayfac = 0.
opticaldepth = 0.

call aldrinit(ifull,iextra,kl,sig)

if ( myid==0 ) then
  allocate( dumg(ifull_g,16) )
  write(6,*) "Opening emissions file ",trim(aerofile)
  call ccnf_open(aerofile,ncid,ncstatus)
  call ncmsg('Aerosol emissions',ncstatus)
  ! check dimensions and location
  call ccnf_get_attg(ncid,'lat0',tlat)
  call ccnf_get_attg(ncid,'lon0',tlon)
  call ccnf_get_attg(ncid,'schmidt0',tschmidt)
  if ( abs(rlong0-tlon) > iotol .or. abs(rlat0-tlat) > iotol .or. abs(schmidt-tschmidt) > iotol ) then
    write(6,*) "ERROR: Grid mismatch for ",trim(aerofile)
    write(6,*) "rlong0,rlat0,schmidt ",rlong0,rlat0,schmidt
    write(6,*) "tlon,tlat,tschmidt   ",tlon,tlat,tschmidt
    call ccmpi_abort(-1)
  end if
  call ccnf_inq_dimlen(ncid,'longitude',tilg)
  if ( tilg /= il_g ) then
    write (6,*) "ERROR: Grid mismatch for ",trim(aerofile)
    write (6,*) "il_g,tilg ",il_g,tilg
    call ccmpi_abort(-1)
  end if
  ! load emission fields
  spos = 1
  npos(1) = il_g
  npos(2) = il_g*6
  write(6,*) "Loading emissions for SO2 anth l1"
  call ccnf_inq_varid(ncid,'so2a1',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate so2a1"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,1))
  write(6,*) "Loading emissions for SO2 anth l2"
  call ccnf_inq_varid(ncid,'so2a2',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate so2a2"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,2))
  write(6,*) "Loading emissions for BC anth l1"
  call ccnf_inq_varid(ncid,'bca1',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate bca1"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,3))
  write(6,*) "Loading emissions for BC anth l2"
  call ccnf_inq_varid(ncid,'bca2',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate bca2"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,4))
  write(6,*) "Loading emissions for OC anth l1"
  call ccnf_inq_varid(ncid,'oca1',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate oca1"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,5))
  write(6,*) "Loading emissions for OC anth l2"
  call ccnf_inq_varid(ncid,'oca2',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate oca2"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,6))
  write(6,*) "Loading emissions for SO2 bio l1"
  call ccnf_inq_varid(ncid,'so2b1',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate so2b1"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,7))
  write(6,*) "Loading emissions for SO2 bio l2"
  call ccnf_inq_varid(ncid,'so2b2',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate so2b2"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,8))
  write(6,*) "Loading emissions for BC bio l1"
  call ccnf_inq_varid(ncid,'bcb1',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate bcb1"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,9))
  write(6,*) "Loading emissions for BC bio l2"
  call ccnf_inq_varid(ncid,'bcb2',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate bcb2"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,10))
  write(6,*) "Loading emissions for OC bio l1"
  call ccnf_inq_varid(ncid,'ocb1',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate ocb1"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,11))
  write(6,*) "Loading emissions for OC bio l2"
  call ccnf_inq_varid(ncid,'ocb2',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate ocb2"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,12))
  write(6,*) "Loading emissions for DMS ocean"
  call ccnf_inq_varid(ncid,'dmso',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate dmso"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,13))
  write(6,*) "Loading emissions for DMS land"
  call ccnf_inq_varid(ncid,'dmst',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate dmst"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,14))
  write(6,*) "Loading emissions for natural organic"
  call ccnf_inq_varid(ncid,'ocna',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate ocna"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,15))
  write(6,*) "Loading emissions for Volcanic SO2"
  call ccnf_inq_varid(ncid,'vso2',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate vso2"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,16))
  call ccmpi_distribute(duma(:,1:16),dumg(:,1:16))
  do i = 1,16
    call aldrloademiss(i,duma(:,i))
  end do
  ! load dust fields
  write(6,*) "Loading emissions for dust (sand)"
  call ccnf_inq_varid(ncid,'sandem',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate sandem"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,1))
  write(6,*) "Loading emissions for dust (slit)"
  call ccnf_inq_varid(ncid,'siltem',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate siltem"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,2))
  write(6,*) "Loading emissions for dust (clay)"
  call ccnf_inq_varid(ncid,'clayem',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate clayem"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg(:,3))
  call ccmpi_distribute(duma(:,1:3),dumg(:,1:3))
  call ccnf_close(ncid)
  do i=1,3
    call aldrloaderod(i,duma(:,i))
  end do
  deallocate( dumg )
  ! load oxidant fields
  write(6,*) "Opening oxidants file ",trim(oxidantfile)
  call ccnf_open(oxidantfile,ncid,ncstatus)
  call ncmsg('Oxidants',ncstatus)
  ! check dimensions and location
  call ccnf_inq_dimlen(ncid,'lon',ilon)
  call ccnf_inq_dimlen(ncid,'lat',ilat)
  call ccnf_inq_dimlen(ncid,'lev',ilev)
  idum(1) = ilon
  idum(2) = ilat
  idum(3) = ilev
  call ccmpi_bcast(idum(1:3),0,comm_world)
  allocate( oxidantprev_g(ifull,ilev,4) )
  allocate( oxidantnow_g(ifull,ilev,4) )
  allocate( oxidantnext_g(ifull,ilev,4) )
  allocate( rlon(ilon), rlat(ilat), rlev(ilev) )
  allocate( oxidantdum(ilon,ilat,ilev,3) )
  sposs = 1
  nposs(1) = ilon
  nposs(2) = ilat
  nposs(3) = ilev
  nposs(4) = 1
  ! use kdate_s as kdate has not yet been defined
  jyear = kdatein/10000
  jmonth = (kdatein-jyear*10000)/100
  premonth = jmonth - 1
  if ( premonth < 1 ) premonth = 12
  nxtmonth = jmonth + 1
  if ( nxtmonth > 12 ) nxtmonth = 1
  write(6,*) "Processing oxidant file for month ",jmonth
  call ccnf_inq_varid(ncid,'lon',varid,tst)
  if (tst) then
    write(6,*) "ERROR: Cannot locate lon"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,sposs(1:1),nposs(1:1),rlon)
  call ccnf_inq_varid(ncid,'lat',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate lat"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,sposs(2:2),nposs(2:2),rlat) ! input latitudes (deg)
  call ccnf_inq_varid(ncid,'lev',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate lev"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,sposs(3:3),nposs(3:3),rlev) ! input vertical levels
  allocate( rpack(ilon+ilat+ilev) )
  rpack(1:ilon) = rlon(1:ilon)
  rpack(ilon+1:ilon+ilat) = rlat(1:ilat)
  rpack(ilon+ilat+1:ilon+ilat+ilev) = rlev(1:ilev)
  call ccmpi_bcast(rpack,0,comm_world)
  deallocate( rpack )
  write(6,*) "Reading OH"
  call ccnf_inq_varid(ncid,'OH',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate OH"
    call ccmpi_abort(-1)
  end if
  sposs(4) = premonth
  call ccnf_get_vara(ncid,varid,sposs,nposs,oxidantdum(:,:,:,1))
  sposs(4) = jmonth
  call ccnf_get_vara(ncid,varid,sposs,nposs,oxidantdum(:,:,:,2))
  sposs(4) = nxtmonth
  call ccnf_get_vara(ncid,varid,sposs,nposs,oxidantdum(:,:,:,3))
  call ccmpi_bcast(oxidantdum(:,:,:,1:3),0,comm_world)
  call o3regrid(oxidantprev_g(:,:,1),oxidantnow_g(:,:,1),oxidantnext_g(:,:,1),oxidantdum,rlon,rlat)
  write(6,*) "Reading H2O2"
  call ccnf_inq_varid(ncid,'H2O2',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate H2O2"
    call ccmpi_abort(-1)
  end if
  sposs(4) = premonth
  call ccnf_get_vara(ncid,varid,sposs,nposs,oxidantdum(:,:,:,1))
  sposs(4) = jmonth
  call ccnf_get_vara(ncid,varid,sposs,nposs,oxidantdum(:,:,:,2))
  sposs(4) = nxtmonth
  call ccnf_get_vara(ncid,varid,sposs,nposs,oxidantdum(:,:,:,3))
  call ccmpi_bcast(oxidantdum(:,:,:,1:3),0,comm_world)
  call o3regrid(oxidantprev_g(:,:,2),oxidantnow_g(:,:,2),oxidantnext_g(:,:,2),oxidantdum,rlon,rlat)
  write(6,*) "Reading O3"
  call ccnf_inq_varid(ncid,'O3',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate O3"
    call ccmpi_abort(-1)
  end if
  sposs(4) = premonth
  call ccnf_get_vara(ncid,varid,sposs,nposs,oxidantdum(:,:,:,1))
  sposs(4) = jmonth
  call ccnf_get_vara(ncid,varid,sposs,nposs,oxidantdum(:,:,:,2))
  sposs(4) = nxtmonth
  call ccnf_get_vara(ncid,varid,sposs,nposs,oxidantdum(:,:,:,3))
  call ccmpi_bcast(oxidantdum(:,:,:,1:3),0,comm_world)
  call o3regrid(oxidantprev_g(:,:,3),oxidantnow_g(:,:,3),oxidantnext_g(:,:,3),oxidantdum,rlon,rlat)
  write(6,*) "Reading NO2"
  call ccnf_inq_varid(ncid,'NO2',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate NO2"
    call ccmpi_abort(-1)
  end if
  sposs(4) = premonth
  call ccnf_get_vara(ncid,varid,sposs,nposs,oxidantdum(:,:,:,1))
  sposs(4) = jmonth
  call ccnf_get_vara(ncid,varid,sposs,nposs,oxidantdum(:,:,:,2))
  sposs(4) = nxtmonth
  call ccnf_get_vara(ncid,varid,sposs,nposs,oxidantdum(:,:,:,3))
  call ccmpi_bcast(oxidantdum(:,:,:,1:3),0,comm_world)
  call o3regrid(oxidantprev_g(:,:,4),oxidantnow_g(:,:,4),oxidantnext_g(:,:,4),oxidantdum,rlon,rlat)
  call ccnf_close(ncid)
  deallocate(oxidantdum,rlat,rlon)
else
  ! load emission fields
  call ccmpi_distribute(duma(:,1:16))  
  do i = 1,16
    call aldrloademiss(i,duma(:,i))
  end do
  ! load dust fields (sand, silt, clay)
  call ccmpi_distribute(duma(:,1:3))  
  do i = 1,3
    call aldrloaderod(i,duma(:,i))
  end do
  ! load oxidant fields
  call ccmpi_bcast(idum(1:3),0,comm_world)
  ilon = idum(1)
  ilat = idum(2)
  ilev = idum(3)
  allocate( oxidantprev_g(ifull,ilev,4) )
  allocate( oxidantnow_g(ifull,ilev,4) )
  allocate( oxidantnext_g(ifull,ilev,4) )
  allocate( rlon(ilon), rlat(ilat), rlev(ilev) )
  allocate( oxidantdum(ilon,ilat,ilev,3) )
  allocate( rpack(ilon+ilat+ilev) )
  call ccmpi_bcast(rpack,0,comm_world)
  rlon(1:ilon) = rpack(1:ilon)
  rlat(1:ilat) = rpack(ilon+1:ilon+ilat)
  rlev(1:ilev) = rpack(ilon+ilat+1:ilon+ilat+ilev)
  deallocate( rpack )
  do j = 1,4
    call ccmpi_bcast(oxidantdum(:,:,:,1:3),0,comm_world)
    call o3regrid(oxidantprev_g(:,:,j),oxidantnow_g(:,:,j),oxidantnext_g(:,:,j),oxidantdum,rlon,rlat)    
  end do
  deallocate(oxidantdum,rlat,rlon)
end if

return
end subroutine load_aerosolldr    
    
end module aerointerface
