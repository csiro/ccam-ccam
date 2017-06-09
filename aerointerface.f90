! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2017 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
public opticaldepth

integer, save :: ilon, ilat, ilev
integer, parameter :: naerofamilies = 7      ! Number of aerosol families for optical depth
real, dimension(:,:,:), allocatable, save :: oxidantprev
real, dimension(:,:,:), allocatable, save :: oxidantnow
real, dimension(:,:,:), allocatable, save :: oxidantnext
real, dimension(:,:,:), allocatable, save :: opticaldepth
real, dimension(:,:), allocatable, save :: ppfprec, ppfmelt, ppfsnow           ! data saved from LDR cloud scheme
real, dimension(:,:), allocatable, save :: ppfevap, ppfsubl, pplambs, ppmrate  ! data saved from LDR cloud scheme
real, dimension(:,:), allocatable, save :: ppmaccr, ppqfsedice, pprscav        ! data saved from LDR cloud scheme
real, dimension(:,:), allocatable, save :: pprfreeze                           ! data saved from LDR cloud scheme
real, dimension(:,:), allocatable, save :: ppfstayice, ppfstayliq              ! data saved from LDR cloud scheme
real, dimension(:), allocatable, save :: rlev, zdayfac
real, parameter :: wlc = 0.2e-3         ! LWC of deep conv cloud (kg/m**3)
integer, save :: imax
integer, dimension(:), allocatable, save :: sday

contains

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

if ( myid==0 ) write(6,*) "Initialising prognostic aerosols"

imax=ifull/ntiles

allocate(sday(ntiles))
sday=-9999

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
  write(6,*) "Reading ",trim(aerofile)
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
  write(6,*) "Reading ",trim(oxidantfile)
  call ccnf_open(oxidantfile,ncid,ncstatus)
  call ncmsg('Oxidants',ncstatus)
  ! check dimensions and location
  call ccnf_inq_dimlen(ncid,'lon',ilon)
  call ccnf_inq_dimlen(ncid,'lat',ilat)
  call ccnf_inq_dimlen(ncid,'lev',ilev)
  write(6,*) "Found oxidant dimensions ",ilon,ilat,ilev
  idum(1) = ilon
  idum(2) = ilat
  idum(3) = ilev
  call ccmpi_bcast(idum(1:3),0,comm_world)
  allocate( oxidantprev(ifull,ilev,4) )
  allocate( oxidantnow(ifull,ilev,4) )
  allocate( oxidantnext(ifull,ilev,4) )
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
  call o3regrid(oxidantprev(:,:,1),oxidantnow(:,:,1),oxidantnext(:,:,1),oxidantdum,rlon,rlat)
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
  call o3regrid(oxidantprev(:,:,2),oxidantnow(:,:,2),oxidantnext(:,:,2),oxidantdum,rlon,rlat)
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
  call o3regrid(oxidantprev(:,:,3),oxidantnow(:,:,3),oxidantnext(:,:,3),oxidantdum,rlon,rlat)
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
  call o3regrid(oxidantprev(:,:,4),oxidantnow(:,:,4),oxidantnext(:,:,4),oxidantdum,rlon,rlat)
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
  allocate( oxidantprev(ifull,ilev,4) )
  allocate( oxidantnow(ifull,ilev,4) )
  allocate( oxidantnext(ifull,ilev,4) )
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
    call o3regrid(oxidantprev(:,:,j),oxidantnow(:,:,j),oxidantnext(:,:,j),oxidantdum,rlon,rlat)    
  end do
  deallocate(oxidantdum,rlat,rlon)
end if

if ( myid==0 ) write(6,*) "Finished initialising prognostic aerosols"

return
end subroutine load_aerosolldr

subroutine aerocalc
use aerosolldr           ! LDR prognostic aerosols
use arrays_m             ! Atmosphere dyamics prognostic arrays
use cc_omp
use cfrac_m              ! Cloud fraction
use cloudmod             ! Prognostic strat cloud
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
integer :: tile, is, ie
real, dimension(imax,ilev,4) :: loxidantprev
real, dimension(imax,ilev,4) :: loxidantnow
real, dimension(imax,ilev,4) :: loxidantnext
real, dimension(imax) :: lps
real, dimension(imax) :: lzdayfac
real, dimension(imax) :: lrlatt
real, dimension(imax) :: lrlongg
real, dimension(imax,kl) :: lphi_nh
real, dimension(imax,kl) :: lt
integer, dimension(imax) :: lkbsav
integer, dimension(imax) :: lktsav
real, dimension(imax) :: lwetfac
real, dimension(imax) :: lpblh
real, dimension(imax) :: ltss
real, dimension(imax) :: lcondc
real, dimension(imax) :: lsnowd
real, dimension(imax) :: lfg
real, dimension(imax) :: leg
real, dimension(imax) :: lu10
real, dimension(imax) :: lustar
real, dimension(imax) :: lzo
logical, dimension(imax) :: lland
real, dimension(imax) :: lfracice
real, dimension(imax) :: lsigmf
real, dimension(imax,kl) :: lqg
real, dimension(imax,kl) :: lqlg
real, dimension(imax,kl) :: lqfg
real, dimension(imax,kl) :: lcfrac
real, dimension(imax) :: lcdtq
real, dimension(imax,kl) :: lppfprec
real, dimension(imax,kl) :: lppfmelt
real, dimension(imax,kl) :: lppfsnow
real, dimension(imax,kl) :: lppfevap
real, dimension(imax,kl) :: lppfsubl
real, dimension(imax,kl) :: lpplambs
real, dimension(imax,kl) :: lppmrate
real, dimension(imax,kl) :: lppmaccr
real, dimension(imax,kl) :: lppfstayice
real, dimension(imax,kl) :: lppfstayliq
real, dimension(imax,kl) :: lppqfsedice
real, dimension(imax,kl) :: lpprscav
real, dimension(imax,kl) :: lpprfreeze
real, dimension(imax) :: lso4t
real, dimension(imax,kl,naero) :: lxtg
real, dimension(imax,4*kl) :: lzoxidant
real, dimension(imax) :: lduste
real, dimension(imax) :: ldustdd
real, dimension(imax,kl,naero) :: lxtosav
real, dimension(imax,kl,naero) :: lxtg_solub
real, dimension(imax) :: ldmsso2o
real, dimension(imax) :: lso2so4o
real, dimension(imax) :: ldust_burden
real, dimension(imax) :: lbc_burden
real, dimension(imax) :: loc_burden
real, dimension(imax) :: ldms_burden
real, dimension(imax) :: lso2_burden
real, dimension(imax) :: lso4_burden
real, dimension(imax,ndcls) :: lerod
real, dimension(imax,kl,2) :: lssn
real, dimension(imax) :: lso2wd
real, dimension(imax) :: lso4wd
real, dimension(imax) :: lbcwd
real, dimension(imax) :: locwd
real, dimension(imax) :: ldustwd
real, dimension(imax,15) :: lemissfield
real, dimension(imax) :: lvso2
real, dimension(imax) :: ldmse
real, dimension(imax) :: lso2e
real, dimension(imax) :: lso4e
real, dimension(imax) :: lbce
real, dimension(imax) :: loce
real, dimension(imax) :: lso2dd
real, dimension(imax) :: lso4dd
real, dimension(imax) :: lbcdd
real, dimension(imax) :: locdd

!$omp parallel do private(is,ie), &
!$omp private(loxidantprev,loxidantnow,loxidantnext,lps,lzdayfac,lrlatt,lrlongg,lphi_nh,lt,lkbsav,lktsav), &
!$omp private(lwetfac,lpblh,ltss,lcondc,lsnowd,lfg,leg,lu10,lustar,lzo,lland,lfracice,lsigmf,lqg,lqlg,lqfg), &
!$omp private(lcfrac,lcdtq,lppfprec,lppfmelt,lppfsnow,lppfevap,lppfsubl,lpplambs,lppmrate,lppmaccr), &
!$omp private(lppfstayice,lppfstayliq,lppqfsedice,lpprscav,lpprfreeze,lso4t,lxtg,lzoxidant,lduste), &
!$omp private(ldustdd,lxtosav,lxtg_solub,ldmsso2o,lso2so4o,ldust_burden,lbc_burden,loc_burden), &
!$omp private(ldms_burden,lso2_burden,lso4_burden,lerod,lssn,lso2wd,lso4wd,lbcwd,locwd,ldustwd), &
!$omp private(lemissfield,lvso2,ldmse,lso2e,lso4e,lbce,loce,lso2dd,lso4dd,lbcdd,locdd)
do tile=1,ntiles
  is=(tile-1)*imax+1
  ie=tile*imax

  loxidantprev=oxidantprev(is:ie,:,:)
  loxidantnow=oxidantnow(is:ie,:,:)
  loxidantnext=oxidantnext(is:ie,:,:)
  lps=ps(is:ie)
  lzdayfac=zdayfac(is:ie)
  lrlatt=rlatt(is:ie)
  lrlongg=rlongg(is:ie)
  lphi_nh=phi_nh(is:ie,:)
  lt=t(is:ie,:)
  lkbsav=kbsav(is:ie)
  lktsav=ktsav(is:ie)
  lwetfac=wetfac(is:ie)
  lpblh=pblh(is:ie)
  ltss=tss(is:ie)
  lcondc=condc(is:ie)
  lsnowd=snowd(is:ie)
  lfg=fg(is:ie)
  leg=eg(is:ie)
  lu10=u10(is:ie)
  lustar=ustar(is:ie)
  lzo=zo(is:ie)
  lland=land(is:ie)
  lfracice=fracice(is:ie)
  lsigmf=sigmf(is:ie)
  lqg=qg(is:ie,:)
  lqlg=qlg(is:ie,:)
  lqfg=qfg(is:ie,:)
  lcfrac=cfrac(is:ie,:)
  lcdtq=cdtq(is:ie)
  lppfprec=ppfprec(is:ie,:)
  lppfmelt=ppfmelt(is:ie,:)
  lppfsnow=ppfsnow(is:ie,:)
  lppfevap=ppfevap(is:ie,:)
  lppfsubl=ppfsubl(is:ie,:)
  lpplambs=pplambs(is:ie,:)
  lppmrate=ppmrate(is:ie,:)
  lppmaccr=ppmaccr(is:ie,:)
  lppfstayice=ppfstayice(is:ie,:)
  lppfstayliq=ppfstayliq(is:ie,:)
  lppqfsedice=ppqfsedice(is:ie,:)
  lpprscav=pprscav(is:ie,:)
  lpprfreeze=pprfreeze(is:ie,:)
  lso4t=so4t(is:ie)
  lxtg=xtg(is:ie,:,:)
  lzoxidant=zoxidant(is:ie,:)
  lduste=duste(is:ie)
  ldustdd=dustdd(is:ie)
  lxtosav=xtosav(is:ie,:,:)
  if ( aeromode>=1 ) then
    lxtg_solub=xtg_solub(is:ie,:,:)
  end if
  ldmsso2o=dmsso2o(is:ie)
  lso2so4o=so2so4o(is:ie)
  ldust_burden=dust_burden(is:ie)
  lbc_burden=bc_burden(is:ie)
  loc_burden=oc_burden(is:ie)
  ldms_burden=dms_burden(is:ie)
  lso2_burden=so2_burden(is:ie)
  lso4_burden=so4_burden(is:ie)
  lerod=erod(is:ie,:)
  lssn=ssn(is:ie,:,:)
  lso2wd=so2wd(is:ie)
  lso4wd=so4wd(is:ie)
  lbcwd=bcwd(is:ie)
  locwd=ocwd(is:ie)
  ldustwd=dustwd(is:ie)
  lemissfield=emissfield(is:ie,:)
  lvso2=vso2(is:ie)
  ldmse=dmse(is:ie)
  lso2e=so2e(is:ie)
  lso4e=so4e(is:ie)
  lbce=bce(is:ie)
  loce=oce(is:ie)
  lso2dd=so2dd(is:ie)
  lso4dd=so4dd(is:ie)
  lbcdd=bcdd(is:ie)
  locdd=ocdd(is:ie)

  call aerocalc_work(loxidantprev,loxidantnow,loxidantnext,lps,lzdayfac,lrlatt,lrlongg,lphi_nh,lt,lkbsav,lktsav,   &
                     lwetfac,lpblh,ltss,lcondc,lsnowd,lfg,leg,lu10,lustar,lzo,lland,lfracice,lsigmf,lqg,lqlg,lqfg, &
                     lcfrac,lcdtq,lppfprec,lppfmelt,lppfsnow,lppfevap,lppfsubl,lpplambs,lppmrate,lppmaccr,         &
                     lppfstayice,lppfstayliq,lppqfsedice,lpprscav,lpprfreeze,lso4t,lxtg,lzoxidant,lduste,ldustdd,  &
                     lxtosav,lxtg_solub,ldmsso2o,lso2so4o,ldust_burden,lbc_burden,loc_burden,ldms_burden,          &
                     lso2_burden,lso4_burden,lerod,lssn,lso2wd,lso4wd,lbcwd,locwd,ldustwd,lemissfield,lvso2,ldmse, &
                     lso2e,lso4e,lbce,loce,lso2dd,lso4dd,lbcdd,locdd,tile,imax)

  zdayfac(is:ie)=lzdayfac
  so4t(is:ie)=lso4t
  xtg(is:ie,:,:)=lxtg
  zoxidant(is:ie,:)=lzoxidant
  duste(is:ie)=lduste
  dustdd(is:ie)=ldustdd
  if ( aeromode>=1 ) then
    xtg_solub(is:ie,:,:)=lxtg_solub
  end if
  dmsso2o(is:ie)=ldmsso2o
  so2so4o(is:ie)=lso2so4o
  dust_burden(is:ie)=ldust_burden
  bc_burden(is:ie)=lbc_burden
  oc_burden(is:ie)=loc_burden
  dms_burden(is:ie)=ldms_burden
  so2_burden(is:ie)=lso2_burden
  so4_burden(is:ie)=lso4_burden
  ssn(is:ie,:,:)=lssn
  so2wd(is:ie)=lso2wd
  so4wd(is:ie)=lso4wd
  bcwd(is:ie)=lbcwd
  ocwd(is:ie)=locwd
  dustwd(is:ie)=ldustwd
  dmse(is:ie)=ldmse
  so2e(is:ie)=lso2e
  so4e(is:ie)=lso4e
  bce(is:ie)=lbce
  oce(is:ie)=loce
  so2dd(is:ie)=lso2dd
  so4dd(is:ie)=lso4dd
  bcdd(is:ie)=lbcdd
  ocdd(is:ie)=locdd

end do

end subroutine aerocalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update prognostic aerosols
subroutine aerocalc_work(oxidantprev,oxidantnow,oxidantnext,ps,zdayfac,rlatt,rlongg,phi_nh,t,kbsav,ktsav,   &
                         wetfac,pblh,tss,condc,snowd,fg,eg,u10,ustar,zo,land,fracice,sigmf,qg,qlg,qfg,      &
                         cfrac,cdtq,ppfprec,ppfmelt,ppfsnow,ppfevap,ppfsubl,pplambs,ppmrate,ppmaccr,        &
                         ppfstayice,ppfstayliq,ppqfsedice,pprscav,pprfreeze,so4t,xtg,zoxidant,duste,dustdd, &
                         xtosav,xtg_solub,dmsso2o,so2so4o,dust_burden,bc_burden,oc_burden,dms_burden,       &
                         so2_burden,so4_burden,erod,ssn,so2wd,so4wd,bcwd,ocwd,dustwd,emissfield,vso2,dmse,  &
                         so2e,so4e,bce,oce,so2dd,so4dd,bcdd,ocdd,tile,imax)

use aerosolldr, only : naero,ndcls,aldrloadoxidant,aldrcalc           ! LDR prognostic aerosols
use cc_mpi               ! CC MPI routines
use cc_omp               ! CC OpenMP routines
use cloudmod, only : convectivecloudfrac             ! Prognostic strat cloud
use const_phys           ! Physical constants
use infile, only : getzinp               ! Input file routines
use newmpar_m            ! Grid parameters
use ozoneread, only : fieldinterpolate            ! Ozone input routines
use parm_m               ! Model configuration
use sigs_m               ! Atmosphere sigma levels
use zenith_m, only : solargh,zenith             ! Astronomy routines

implicit none

include 'kuocom.h'      ! Convection parameters

integer, intent(in) :: tile,imax
integer jyear,jmonth,jday,jhour,jmin,mins,smins
integer j,k,tt,ttx
integer, parameter :: updateoxidant = 1440 ! update prescribed oxidant fields once per day
real dhr,fjd,r1,dlt,alp,slag
!global
real, dimension(imax,ilev,4), intent(in) :: oxidantprev
real, dimension(imax,ilev,4), intent(in) :: oxidantnow
real, dimension(imax,ilev,4), intent(in) :: oxidantnext
real, dimension(imax), intent(in) :: ps
real, dimension(imax), intent(inout) :: zdayfac
real, dimension(imax), intent(in) :: rlatt
real, dimension(imax), intent(in) :: rlongg
real, dimension(imax,kl), intent(in) :: phi_nh
real, dimension(imax,kl), intent(in) :: t
integer, dimension(imax), intent(in) :: kbsav
integer, dimension(imax), intent(in) :: ktsav
real, dimension(imax), intent(in) :: wetfac
real, dimension(imax), intent(in) :: pblh
real, dimension(imax), intent(in) :: tss
real, dimension(imax), intent(in) :: condc
real, dimension(imax), intent(in) :: snowd
real, dimension(imax), intent(in) :: fg
real, dimension(imax), intent(in) :: eg
real, dimension(imax), intent(in) :: u10
real, dimension(imax), intent(in) :: ustar
real, dimension(imax), intent(in) :: zo
logical, dimension(imax), intent(in) :: land
real, dimension(imax), intent(in) :: fracice
real, dimension(imax), intent(in) :: sigmf
real, dimension(imax,kl), intent(in) :: qg
real, dimension(imax,kl), intent(in) :: qlg
real, dimension(imax,kl), intent(in) :: qfg
real, dimension(imax,kl), intent(in) :: cfrac
real, dimension(imax), intent(in) :: cdtq
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
real, dimension(imax), intent(inout) :: so4t
real, dimension(imax,kl,naero), intent(inout) :: xtg
real, dimension(imax,4*kl), intent(inout) :: zoxidant
real, dimension(imax), intent(inout) :: duste
real, dimension(imax), intent(inout) :: dustdd
real, dimension(imax,kl,naero), intent(in) :: xtosav
real, dimension(imax,kl,naero), intent(inout) :: xtg_solub
real, dimension(imax), intent(inout) :: dmsso2o
real, dimension(imax), intent(inout) :: so2so4o
real, dimension(imax), intent(inout) :: dust_burden
real, dimension(imax), intent(inout) :: bc_burden
real, dimension(imax), intent(inout) :: oc_burden
real, dimension(imax), intent(inout) :: dms_burden
real, dimension(imax), intent(inout) :: so2_burden
real, dimension(imax), intent(inout) :: so4_burden
real, dimension(imax,ndcls), intent(in) :: erod
real, dimension(imax,kl,2), intent(inout) :: ssn
real, dimension(imax), intent(inout) :: so2wd
real, dimension(imax), intent(inout) :: so4wd
real, dimension(imax), intent(inout) :: bcwd
real, dimension(imax), intent(inout) :: ocwd
real, dimension(imax), intent(inout) :: dustwd
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
!
real, dimension(imax,kl) :: oxout,zg,clcon,pccw,rhoa
real, dimension(imax,kl) :: tnhs,dz
real, dimension(imax) :: coszro,taudar
real, dimension(imax) :: cldcon,wg
real, dimension(kl+1) :: sigh
integer :: nthreads

nthreads=ccomp_get_num_threads()

! timer calculations
call getzinp(jyear,jmonth,jday,jhour,jmin,mins)
! update prescribed oxidant fields
dhr = dt/3600.
if ( sday(tile)<=mins-updateoxidant ) then
  sday(tile) = mins
  do j = 1,4 
    ! note levels are inverted by fieldinterpolate
    call fieldinterpolate(oxout,oxidantprev(:,:,j),oxidantnow(:,:,j),oxidantnext(:,:,j), &
                          rlev,imax,kl,ilev,mins,sig,ps,interpmeth=0)
    do k = 1,kl
      call aldrloadoxidant(k+(j-1)*kl,oxout(:,k),zoxidant,imax)
    end do
  end do
  ! estimate day length (presumably to preturb day-time OH levels)
  ttx = nint(86400./dt)
  zdayfac(:) = 0.
  do tt = ttx,1,-1 ! we seem to get a different answer if dhr=24. and ttx=1.
    smins = int(real(tt-1)*dt/60.)+mins
    fjd = float(mod( smins, 525600 ))/1440.  ! 525600 = 1440*365
    call solargh(fjd,bpyear,r1,dlt,alp,slag)
    call zenith(fjd,r1,dlt,slag,rlatt,rlongg,dhr,imax,coszro,taudar)
    where ( taudar>0.5 )
      zdayfac(:) = zdayfac(:) + 1.
    end where
  end do
  ! final taudar is for current timestep - used to indicate sunlit
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
sigh(1:kl) = sigmh(1:kl) ! store half-levels
sigh(kl+1) = 0.

! Non-hydrostatic terms
tnhs(:,1) = phi_nh(:,1)/bet(1)
zg(:,1) = bet(1)*t(1:imax,1)/grav
do k = 2,kl
  ! representing non-hydrostatic term as a correction to air temperature
  tnhs(:,k) = (phi_nh(:,k)-phi_nh(:,k-1)-betm(k)*tnhs(:,k-1))/bet(k)
  zg(:,k) = zg(:,k-1) + (bet(k)*t(1:imax,k)+betm(k)*t(1:imax,k-1))/grav ! height above surface in meters
end do
do k = 1,kl
  zg(:,k) = zg(:,k) + phi_nh(:,k)/grav
  dz(:,k) = -rdry*dsig(k)*(t(1:imax,k)+tnhs(:,k))/(grav*sig(k))
  rhoa(:,k) = ps(1:imax)*sig(k)/(rdry*t(1:imax,k)) ! density of air (kg/m**3)
end do

! estimate convective cloud fraction from leoncld.f
call convectivecloudfrac(clcon,kbsav,ktsav,condc,imax,cldcon=cldcon)
do k = 1,kl
  ! MJT notes - Assume rain for JLM convection
  !where ( k>kbsav .and. k<=ktsav .and. t(1:imax,k)>ticeu )
  !  pccw(:,kl+1-k) = 0.
  where ( k>kbsav .and. k<=ktsav )
    pccw(:,kl+1-k) = wlc/rhoa(:,k)
  elsewhere
    pccw(:,kl+1-k) = 0.
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
call aldrcalc(dt,sig,zg,dz,wg,pblh,ps,tss,                 &
              t,condc,snowd,taudar,fg,eg,u10,ustar,zo,     &
              land,fracice,sigmf,qg,qlg,qfg,cfrac,clcon,   &
              cldcon,pccw,rhoa,cdtq,ppfprec,ppfmelt,       &
              ppfsnow,ppfevap,ppfsubl,pplambs,ppmrate,     &
              ppmaccr,ppfstayice,ppfstayliq,ppqfsedice,    &
              pprscav,pprfreeze,zdayfac,kbsav,xtg,duste,   &
              dustdd,xtosav,xtg_solub,dmsso2o,so2so4o,     &
              dust_burden,bc_burden,oc_burden,dms_burden,  &
              so2_burden,so4_burden,erod,ssn,zoxidant,     &
              so2wd,so4wd,bcwd,ocwd,dustwd,emissfield,     &
              vso2,dmse,so2e,so4e,bce,oce,so2dd,so4dd,     &
              bcdd,ocdd,imax)
              

! store sulfate for LH+SF radiation scheme.  SEA-ESF radiation scheme imports prognostic aerosols in seaesfrad.f90.
! Factor 1.e3 to convert to gS/m2, x 3 to get sulfate from sulfur
so4t(:) = 0.
do k = 1,kl
  so4t(:) = so4t(:) + 3.e3*xtg(1:imax,k,3)*rhoa(:,k)*dz(:,k)
enddo

if ( diag .and. mydiag .and. nthreads==1 ) then
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

select case(indirmode)
  case(2,3)
    ! prognostic aerosols for indirect effects
    call cldrop(istart,imax,cdn,rhoa,convmode)
  case default
    ! diagnosed for prescribed aerosol indirect effects
    iend = istart + imax - 1
    where (land(istart:iend).and.rlatt(istart:iend)>0.)
      cdn(:,1)=cdropl_nh
    elsewhere (land(istart:iend))
      cdn(:,1)=cdropl_sh
    elsewhere (rlatt(istart:iend)>0.)
      cdn(:,1)=cdrops_nh
    elsewhere
      cdn(:,1)=cdrops_sh
    end where
    do k = 2,kl
      cdn(:,k)=cdn(:,1)
    end do
end select

return
end subroutine aerodrop

end module aerointerface
