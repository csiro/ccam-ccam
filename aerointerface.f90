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
real, dimension(:,:,:), allocatable, save :: oxidantnow_g
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

use aerosolldr                              ! LDR prognostic aerosols
use arrays_m                                ! Atmosphere dyamics prognostic arrays
use cc_mpi, only : mydiag                   ! CC MPI routines
use cc_omp                                  ! CC OpenMP routines
use cfrac_m                                 ! Cloud fraction
use cloudmod, only : convectivecloudfrac    ! Prognostic strat cloud
use const_phys                              ! Physical constants
use extraout_m                              ! Additional diagnostics
use infile                                  ! Input file routines
use kuocomb_m                               ! JLM convection
use latlong_m                               ! Lat/lon coordinates
use liqwpar_m                               ! Cloud water mixing ratios
use morepbl_m                               ! Additional boundary layer diagnostics
use newmpar_m                               ! Grid parameters
use nharrs_m                                ! Non-hydrostatic atmosphere arrays
use nsibd_m                                 ! Land-surface arrays
use ozoneread                               ! Ozone input routines
use parm_m                                  ! Model configuration
use pbl_m                                   ! Boundary layer arrays
use screen_m                                ! Screen level diagnostics
use sigs_m                                  ! Atmosphere sigma levels
use soil_m                                  ! Soil and surface data
use soilsnow_m                              ! Soil, snow and surface data
use soilv_m                                 ! Soil parameters
use vegpar_m                                ! Vegetation arrays
use work2_m                                 ! Diagnostic arrays
use zenith_m                                ! Astronomy routines

implicit none

include 'kuocom.h'                          ! Convection parameters

integer, intent(in) :: mins
integer tile, is, ie, idjd_t
integer k, j, tt, ttx, kinv, smins
real, dimension(imax,ilev) :: loxidantnow
real, dimension(imax,kl,naero) :: lxtg, lxtosav
real, dimension(imax,kl,4) :: lzoxidant
real, dimension(imax,kl) :: lt, lqg, lqlg, lqfg, lstratcloud
real, dimension(imax,kl) :: lppfprec, lppfmelt, lppfsnow, lppfsubl, lpplambs
real, dimension(imax,kl) :: lppmrate, lppmaccr, lppfstayice, lppqfsedice
real, dimension(imax,kl) :: lpprscav, lpprfreeze
real, dimension(imax,kl) :: lclcon
real, dimension(imax,kl) :: dz, rhoa, pccw
real, dimension(imax,ndust) :: lduste, ldustdd, ldust_burden, ldustwd
real, dimension(imax,ndcls) :: lerod
real, dimension(imax,15) :: lemissfield
real, dimension(imax) :: coszro, wg
real, dimension(ifull,kl) :: clcon
real, dimension(ifull) :: taudar, cldcon
real dhr, fjd, r1, dlt, alp, slag
logical, intent(in) :: oxidant_update
logical mydiag_t
logical, dimension(imax) :: locean

! update prescribed oxidant fields
if ( oxidant_update ) then
!$omp do schedule(static) private(is,ie),                       &
!$omp private(loxidantnow,lzoxidant),                           &
!$omp private(dhr,ttx,j,tt,smins,fjd,r1,dlt,alp,slag,coszro)
  do tile = 1,ntiles
    is = (tile-1)*imax + 1
    ie = tile*imax
    do j = 1,4
      loxidantnow(:,1:ilev)  = oxidantnow_g(is:ie,1:ilev,j)
      ! note levels are inverted by fieldinterpolate
      call fieldinterpolate(lzoxidant(:,:,j),loxidantnow,rlev,imax,kl,ilev,sig,ps(is:ie),interpmeth=0)
      zoxidant_g(is:ie,:,j) = lzoxidant(:,:,j)
    end do
    ! estimate day length (presumably to preturb day-time OH levels)
    ttx = nint(86400./dt)
    dhr = dt/3600.
    zdayfac(is:ie) = 0.
    do tt = ttx,1,-1 ! we seem to get a different answer if dhr=24. and ttx=1.
      smins = int(real(tt-1)*dt/60.) + mins
      fjd = float(mod( smins, 525600 ))/1440.  ! 525600 = 1440*365
      call solargh(fjd,bpyear,r1,dlt,alp,slag)
      call zenith(fjd,r1,dlt,slag,rlatt(is:ie),rlongg(is:ie),dhr,imax,coszro,taudar(is:ie))
      where ( taudar(is:ie)>0.5 )
        zdayfac(is:ie) = zdayfac(is:ie) + 1.
      end where
    end do
    where ( zdayfac(is:ie)>0.5 )
      zdayfac(is:ie) = real(ttx)/zdayfac(is:ie)
    end where
  end do
!$omp end do nowait
else
!$omp do schedule(static) private(is,ie),                       &
!$omp private(dhr,fjd,coszro,r1,dlt,alp,slag)
  do tile = 1,ntiles
    is = (tile-1)*imax + 1
    ie = tile*imax
    dhr = dt/3600.  
    fjd = float(mod( mins, 525600 ))/1440.  ! 525600 = 1440*365
    call solargh(fjd,bpyear,r1,dlt,alp,slag)
    call zenith(fjd,r1,dlt,slag,rlatt(is:ie),rlongg(is:ie),dhr,imax,coszro,taudar(is:ie))
    ! taudar is for current timestep - used to indicate sunlit
  end do
!$omp end do nowait
end if
    
! estimate convective cloud fraction from leoncld.f
!$omp do schedule(static) private(is,ie,lclcon)
do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  call convectivecloudfrac(lclcon,kbsav(is:ie),ktsav(is:ie),condc(is:ie),cldcon=cldcon(is:ie))
  clcon(is:ie,:) = lclcon    
end do
!$omp end do nowait    
    
!$omp do schedule(static) private(is,ie,idjd_t,mydiag_t),                                              &
!$omp private(k,dz,rhoa,wg,pccw,kinv,lt,lqg,lqlg,lqfg),                                                &
!$omp private(lstratcloud,lppfprec,lppfmelt,lppfsnow,lppfsubl,lpplambs,lppmrate,lppmaccr),             &
!$omp private(lppfstayice,lppqfsedice,lpprscav,lpprfreeze,lxtg,lzoxidant,lduste,ldustdd),              &
!$omp private(lxtosav,ldust_burden,lerod,ldustwd,lemissfield,lclcon,locean)
!$acc parallel copy(xtg,duste,dustdd,dustwd,dust_burden,so4t,dmsso2o,so2so4o,bc_burden,oc_burden,     &
!$acc   dms_burden,so2_burden,so4_burden,so2wd,so4wd,bcwd,ocwd,dmse,so2e,so4e,bce,oce,so2dd,so4dd,    &
!$acc   bcdd,ocdd,salte,saltdd,saltwd,salt_burden)                                                    &
!$acc copyin(zoxidant_g,xtosav,emissfield,erod,t,qg,qlg,qfg,stratcloud,ppfprec,ppfmelt,ppfsnow,       &
!$acc   ppfsubl,pplambs,ppmrate,ppmaccr,ppfstayice,ppqfsedice,pprscav,pprfreeze,                      &
!$acc   clcon,bet,betm,dsig,sig,ps,kbsav,ktsav,wetfac,pblh,tss,condc,snowd,taudar,fg,eg,u10,ustar,    &
!$acc   zo,land,fracice,sigmf,cldcon,cdtq,zdayfac,vso2,isoilm_in,dustden,dustreff,saltden,saltreff)
!$acc loop gang private(lzoxidant,lxtg,lxtosav,lduste,ldustdd,ldustwd,ldust_burden,lemissfield,       &
!$acc   lerod,lt,lqg,lqlg,lqfg,lstratcloud,lppfprec,lppfmelt,lppfsnow,lppfsubl,lpplambs,              &
!$acc   lppmrate,lppmaccr,lppfstayice,lppqfsedice,lpprscav,lpprfreeze,lclcon,dz,rhoa,                 &
!$acc   wg,pccw,locean)
do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  
  idjd_t = mod(idjd-1,imax)+1
  mydiag_t = ((idjd-1)/imax==tile-1).and.mydiag 

  lzoxidant(:,:,1:4) = zoxidant_g(is:ie,:,1:4)
  lxtg               = xtg(is:ie,:,:)
  lxtosav            = xtosav(is:ie,:,:)
  lduste             = duste(is:ie,:)
  ldustdd            = dustdd(is:ie,:)
  ldustwd            = dustwd(is:ie,:)
  ldust_burden       = dust_burden(is:ie,:)
  lemissfield        = emissfield(is:ie,:)
  lerod              = erod(is:ie,:)
  lt                 = t(is:ie,:)
  lqg                = qg(is:ie,:)
  lqlg               = qlg(is:ie,:)
  lqfg               = qfg(is:ie,:)
  lstratcloud        = stratcloud(is:ie,:)
  lppfprec           = ppfprec(is:ie,:)
  lppfmelt           = ppfmelt(is:ie,:)
  lppfsnow           = ppfsnow(is:ie,:)
  lppfsubl           = ppfsubl(is:ie,:)
  lpplambs           = pplambs(is:ie,:)
  lppmrate           = ppmrate(is:ie,:)
  lppmaccr           = ppmaccr(is:ie,:)
  lppfstayice        = ppfstayice(is:ie,:)
  lppqfsedice        = ppqfsedice(is:ie,:)
  lpprscav           = pprscav(is:ie,:)
  lpprfreeze         = pprfreeze(is:ie,:)
  lclcon             = clcon(is:ie,:)

  !zg(:,1) = bet(1)*lt(:,1)/grav
  !do k = 2,kl
  !  zg(:,k) = zg(:,k-1) + (bet(k)*lt(:,k)+betm(k)*lt(:,k-1))/grav ! height above surface in meters
  !end do
  do k = 1,kl
    dz(:,k) = -rdry*dsig(k)*lt(:,k)/(grav*sig(k))
    dz(:,k) = min( max( dz(:,k), 1. ), 2.e4 )
    rhoa(:,k) = ps(is:ie)*sig(k)/(rdry*lt(:,k)) ! density of air (kg/m**3)
  end do

  pccw(:,:) = 0.
  do k = 1,kl
    kinv = kl + 1 - k  
    ! MJT notes - Assume rain for JLM convection
    !where ( k>kbsav(is:ie) .and. k<=ktsav(is:ie) .and. lt(:,k)>ticeu )
    !  pccw(:,kl+1-k) = 0.
    where ( k>kbsav(is:ie) .and. k<=ktsav(is:ie) )
      pccw(:,kinv) = wlc/rhoa(:,k)
    end where
  end do

  ! Water converage at surface
  wg(:) = min( max( wetfac(is:ie), 0. ), 1. )
  
  locean = isoilm_in(is:ie) == 0 ! excludes lakes


  ! MJT notes - We have an option to update the aerosols before the vertical mixing
  ! or after the vertical mixing.  Updating aerosols before the vertical mixing
  ! ensures that we can split the convective and non-convective aerosol
  ! concentrations.  However, updating aerosols after vertical mixing provides a
  ! better estimate of u10 and pblh.

  ! update prognostic aerosols
  call aldrcalc(dt,sig,dz,wg,pblh(is:ie),ps(is:ie),tss(is:ie),         &
                lt,condc(is:ie),snowd(is:ie),taudar(is:ie),fg(is:ie),  &
                eg(is:ie),u10(is:ie),ustar(is:ie),zo(is:ie),           &
                land(is:ie),fracice(is:ie),sigmf(is:ie),lqg,lqlg,lqfg, &
                lstratcloud,lclcon,cldcon(is:ie),pccw,rhoa,            &
                cdtq(is:ie),lppfprec,lppfmelt,lppfsnow,                &
                lppfsubl,lpplambs,lppmrate,lppmaccr,lppfstayice,       &
                lppqfsedice,lpprscav,lpprfreeze,                       &
                zdayfac(is:ie),kbsav(is:ie),lxtg,lduste,ldustdd,       &
                lxtosav,dmsso2o(is:ie),so2so4o(is:ie),                 &
                ldust_burden,bc_burden(is:ie),oc_burden(is:ie),        &
                dms_burden(is:ie),so2_burden(is:ie),so4_burden(is:ie), &
                lerod,lzoxidant,so2wd(is:ie),so4wd(is:ie),             &
                bcwd(is:ie),ocwd(is:ie),ldustwd,lemissfield,           &
                vso2(is:ie),dmse(is:ie),so2e(is:ie),so4e(is:ie),       &
                bce(is:ie),oce(is:ie),so2dd(is:ie),so4dd(is:ie),       &
                bcdd(is:ie),ocdd(is:ie),salte(is:ie),saltdd(is:ie),    &
                saltwd(is:ie),salt_burden(is:ie),dustden,dustreff,     &
                saltden,saltreff,locean,imax,kl)

  ! MJT notes - passing dustden, dustreff, saltden and saltref due to issues with pgi compiler
  
  ! store sulfate for LH+SF radiation scheme.  SEA-ESF radiation scheme imports prognostic aerosols in seaesfrad.f90.
  ! Factor 1.e3 to convert to gS/m2, x 3 to get sulfate from sulfur
  so4t(is:ie) = 0.
  do k = 1,kl
    so4t(is:ie) = so4t(is:ie) + 3.e3*lxtg(:,k,3)*rhoa(:,k)*dz(:,k)
  enddo
   
  xtg(is:ie,:,:)       = lxtg
  duste(is:ie,:)       = lduste
  dustdd(is:ie,:)      = ldustdd
  dustwd(is:ie,:)      = ldustwd
  dust_burden(is:ie,:) = ldust_burden
  
#ifndef GPU  
  if ( diag .and. mydiag_t ) then
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
    write(6,*) "saltfilmdiag ",xtg(idjd,:,12)
    write(6,*) "saltjetdiag  ",xtg(idjd,:,13)
  end if
#endif
  
end do
!$acc end parallel    
!$omp end do nowait

return
end subroutine aerocalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Estimate cloud droplet size
subroutine aerodrop(istart,cdn,rhoa,outconv)

use aerosolldr              ! LDR prognostic aerosols
use const_phys              ! Physical constants
use latlong_m, only : rlatt ! Lat/lon coordinates
use newmpar_m               ! Grid parameters
use parm_m                  ! Model configuration
use soil_m, only : land     ! Soil and surface data

implicit none

integer, intent(in) :: istart
integer k,indirmode,iend,imax
real, dimension(:,:), intent(out) :: cdn
real, dimension(:,:), intent(in) :: rhoa
real, parameter :: cdrops_nh=1.e8, cdropl_nh=3.e8 !Cloud droplet conc sea/land nh
!real, parameter :: cdrops_sh=1.e8, cdropl_sh=3.e8 !Cloud droplet conc sea/land sh
real, parameter :: cdrops_sh=.5e8, cdropl_sh=1.e8 !Cloud droplet conc sea/land sh
logical, intent(in), optional :: outconv
logical convmode

imax = size(cdn,1)

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
    call cldrop(istart,cdn,rhoa,convmode)
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
integer, dimension(2) :: spos, npos
integer, dimension(3) :: idum
integer, dimension(4) :: sposs, nposs
real, dimension(:,:), allocatable, save :: dumg
real, dimension(ifull,16) :: duma
real, dimension(:,:,:), allocatable, save :: oxidantdum
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
  allocate( oxidantnow_g(ifull,ilev,4) )
  allocate( rlon(ilon), rlat(ilat), rlev(ilev) )
  allocate( oxidantdum(ilon,ilat,ilev) )
  sposs = 1
  nposs(1) = ilon
  nposs(2) = ilat
  nposs(3) = ilev
  nposs(4) = 1
  ! use kdate_s as kdate has not yet been defined
  jyear = kdatein/10000
  jmonth = (kdatein-jyear*10000)/100
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
  sposs(4) = jmonth
  call ccnf_get_vara(ncid,varid,sposs,nposs,oxidantdum)
  call ccmpi_bcast(oxidantdum,0,comm_world)
  call o3regrid(oxidantnow_g(:,:,1),oxidantdum,rlon,rlat)
  write(6,*) "Reading H2O2"
  call ccnf_inq_varid(ncid,'H2O2',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate H2O2"
    call ccmpi_abort(-1)
  end if
  sposs(4) = jmonth
  call ccnf_get_vara(ncid,varid,sposs,nposs,oxidantdum)
  call ccmpi_bcast(oxidantdum,0,comm_world)
  call o3regrid(oxidantnow_g(:,:,2),oxidantdum,rlon,rlat)
  write(6,*) "Reading O3"
  call ccnf_inq_varid(ncid,'O3',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate O3"
    call ccmpi_abort(-1)
  end if
  sposs(4) = jmonth
  call ccnf_get_vara(ncid,varid,sposs,nposs,oxidantdum)
  call ccmpi_bcast(oxidantdum,0,comm_world)
  call o3regrid(oxidantnow_g(:,:,3),oxidantdum,rlon,rlat)
  write(6,*) "Reading NO2"
  call ccnf_inq_varid(ncid,'NO2',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate NO2"
    call ccmpi_abort(-1)
  end if
  sposs(4) = jmonth
  call ccnf_get_vara(ncid,varid,sposs,nposs,oxidantdum)
  call ccmpi_bcast(oxidantdum,0,comm_world)
  call o3regrid(oxidantnow_g(:,:,4),oxidantdum,rlon,rlat)
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
  allocate( oxidantnow_g(ifull,ilev,4) )
  allocate( rlon(ilon), rlat(ilat), rlev(ilev) )
  allocate( oxidantdum(ilon,ilat,ilev) )
  allocate( rpack(ilon+ilat+ilev) )
  call ccmpi_bcast(rpack,0,comm_world)
  rlon(1:ilon) = rpack(1:ilon)
  rlat(1:ilat) = rpack(ilon+1:ilon+ilat)
  rlev(1:ilev) = rpack(ilon+ilat+1:ilon+ilat+ilev)
  deallocate( rpack )
  do j = 1,4
    call ccmpi_bcast(oxidantdum,0,comm_world)
    call o3regrid(oxidantnow_g(:,:,j),oxidantdum,rlon,rlat)    
  end do
  deallocate(oxidantdum,rlat,rlon)
end if

return
end subroutine load_aerosolldr    
    
end module aerointerface
