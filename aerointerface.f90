! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2016 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
integer, save :: nb,imax
integer, dimension(:), allocatable, save :: sday

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Load aerosols emissions from netcdf
subroutine load_aerosolldr(aerofile, oxidantfile, kdatein, nbin)
      
use aerodata_m, only : aerodata_init
use aerosolldr          ! LDR prognostic aerosols
use cc_mpi              ! CC MPI routines
use infile              ! Input file routines
use newmpar_m           ! Grid parameters
use ozoneread           ! Ozone input routines
use parmgeom_m          ! Coordinate data
use sigs_m              ! Atmosphere sigma levels
      
implicit none

integer, intent(in) :: kdatein,nbin
integer ncstatus, ncid, i, j, varid, tilg
integer jyear, jmonth
integer premonth, nxtmonth
integer, dimension(2) :: spos, npos
integer, dimension(3) :: idum
integer, dimension(4) :: sposs, nposs
real, dimension(:), allocatable, save :: dumg
real, dimension(ifull) :: duma
real, dimension(:,:,:,:), allocatable, save :: oxidantdum
real, dimension(:), allocatable, save :: rlon, rlat
real, dimension(:), allocatable, save :: rpack
real tlat, tlon, tschmidt
real, parameter :: iotol = 1.E-5 ! tolarance for iotest
character(len=*), intent(in) :: aerofile, oxidantfile
logical tst

if ( myid==0 ) write(6,*) "Initialising prognostic aerosols"

nb=nbin
imax=ifull/nb

allocate(sday(nb))
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

if ( myid == 0 ) then
  allocate( dumg(ifull_g) )
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
  call ccnf_get_vara(ncid,varid,spos,npos,dumg)
  call ccmpi_distribute(duma,dumg)
  call aldrloademiss(1,duma)
  write(6,*) "Loading emissions for SO2 anth l2"
  call ccnf_inq_varid(ncid,'so2a2',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate so2a2"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg)
  call ccmpi_distribute(duma,dumg)
  call aldrloademiss(2,duma)
  write(6,*) "Loading emissions for BC anth l1"
  call ccnf_inq_varid(ncid,'bca1',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate bca1"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg)
  call ccmpi_distribute(duma,dumg)
  call aldrloademiss(3,duma)
  write(6,*) "Loading emissions for BC anth l2"
  call ccnf_inq_varid(ncid,'bca2',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate bca2"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg)
  call ccmpi_distribute(duma,dumg)
  call aldrloademiss(4,duma)
  write(6,*) "Loading emissions for OC anth l1"
  call ccnf_inq_varid(ncid,'oca1',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate oca1"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg)
  call ccmpi_distribute(duma,dumg)
  call aldrloademiss(5,duma)
  write(6,*) "Loading emissions for OC anth l2"
  call ccnf_inq_varid(ncid,'oca2',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate oca2"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg)
  call ccmpi_distribute(duma,dumg)
  call aldrloademiss(6,duma)
  write(6,*) "Loading emissions for SO2 bio l1"
  call ccnf_inq_varid(ncid,'so2b1',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate so2b1"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg)
  call ccmpi_distribute(duma,dumg)
  call aldrloademiss(7,duma)
  write(6,*) "Loading emissions for SO2 bio l2"
  call ccnf_inq_varid(ncid,'so2b2',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate so2b2"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg)
  call ccmpi_distribute(duma,dumg)
  call aldrloademiss(8,duma)
  write(6,*) "Loading emissions for BC bio l1"
  call ccnf_inq_varid(ncid,'bcb1',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate bcb1"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg)
  call ccmpi_distribute(duma,dumg)
  call aldrloademiss(9,duma)
  write(6,*) "Loading emissions for BC bio l2"
  call ccnf_inq_varid(ncid,'bcb2',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate bcb2"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg)
  call ccmpi_distribute(duma,dumg)
  call aldrloademiss(10,duma)
  write(6,*) "Loading emissions for OC bio l1"
  call ccnf_inq_varid(ncid,'ocb1',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate ocb1"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg)
  call ccmpi_distribute(duma,dumg)
  call aldrloademiss(11,duma)
  write(6,*) "Loading emissions for OC bio l2"
  call ccnf_inq_varid(ncid,'ocb2',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate ocb2"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg)
  call ccmpi_distribute(duma,dumg)
  call aldrloademiss(12,duma)
  write(6,*) "Loading emissions for DMS ocean"
  call ccnf_inq_varid(ncid,'dmso',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate dmso"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg)
  call ccmpi_distribute(duma,dumg)
  call aldrloademiss(13,duma)
  write(6,*) "Loading emissions for DMS land"
  call ccnf_inq_varid(ncid,'dmst',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate dmst"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg)
  call ccmpi_distribute(duma,dumg)
  call aldrloademiss(14,duma)
  write(6,*) "Loading emissions for natural organic"
  call ccnf_inq_varid(ncid,'ocna',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate ocna"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg)
  call ccmpi_distribute(duma,dumg)
  call aldrloademiss(15,duma)
  write(6,*) "Loading emissions for Volcanic SO2"
  call ccnf_inq_varid(ncid,'vso2',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate vso2"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg)
  call ccmpi_distribute(duma,dumg)
  call aldrloademiss(16,duma)
  ! load dust fields
  write(6,*) "Loading emissions for dust (sand)"
  call ccnf_inq_varid(ncid,'sandem',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate sandem"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg)
  call ccmpi_distribute(duma,dumg)
  call aldrloaderod(1,duma)
  write(6,*) "Loading emissions for dust (slit)"
  call ccnf_inq_varid(ncid,'siltem',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate siltem"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg)
  call ccmpi_distribute(duma,dumg)
  call aldrloaderod(2,duma)
  write(6,*) "Loading emissions for dust (clay)"
  call ccnf_inq_varid(ncid,'clayem',varid,tst)
  if ( tst ) then
    write(6,*) "ERROR: Cannot locate clayem"
    call ccmpi_abort(-1)
  end if
  call ccnf_get_vara(ncid,varid,spos,npos,dumg)
  call ccmpi_distribute(duma,dumg)
  call aldrloaderod(3,duma)
  call ccnf_close(ncid)
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
  do i = 1,16
    call ccmpi_distribute(duma)
    call aldrloademiss(i,duma)
  end do
  ! load dust fields (sand, silt, clay)
  do i = 1,3
    call ccmpi_distribute(duma)
    call aldrloaderod(i,duma)
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

call aerodata_init(nb,imax,ilev,kl,naero,aeromode,ndcls)

if ( myid==0 ) write(6,*) "Finished initialising prognostic aerosols"

return
end subroutine load_aerosolldr

subroutine aerocalc
use cc_mpi, only : start_log,end_log,aero_begin,aero_end
use newmpar_m, only : kl
use aerodata_m
use aerosolldr, only : naero,aeromode,zoxidant,xtg,duste,dustdd,xtosav,xtg_solub,dmsso2o,so2so4o,dust_burden, &
                       bc_burden,dms_burden,so2_burden,so4_burden,EMISSFIELD,vso2,dmse,so2e,so4e,bce,oce,so2dd, &
                       so4dd,bcdd,ocdd,so2wd,so4wd,bcwd,ocwd,dustwd,ndcls,ssn,erod,oc_burden
use kuocomb_m, only : kbsav,ktsav
use arrays_m, only : ps,t,qg
use latlong_m, only : rlatt,rlongg
use nharrs_m, only : phi_nh
use morepbl_m, only : pblh,condc,fg,eg
use work2_m, only : wetfac,zo
use pbl_m, only : tss,cdtq
use soilsnow_m, only : snowd,fracice
use screen_m, only : u10
use extraout_m, only : ustar
use soil_m, only : land,so4t
use nsibd_m, only : sigmf
use liqwpar_m, only : qlg,qfg
use cfrac_m, only : cfrac

implicit none
integer :: i,is,ie

do i=1,nb
  is=(i-1)*imax+1
  ie=i*imax
  b_land(i)%data=land(is:ie)

  b_kbsav(i)%data=kbsav(is:ie)
  b_ktsav(i)%data=ktsav(is:ie)

  b_duste(i)%data=duste(is:ie)
  b_dustdd(i)%data=dustdd(is:ie)
  b_dmsso2o(i)%data=dmsso2o(is:ie)
  b_so2so4o(i)%data=so2so4o(is:ie)
  b_dust_burden(i)%data=dust_burden(is:ie)
  b_bc_burden(i)%data=bc_burden(is:ie)
  b_oc_burden(i)%data=oc_burden(is:ie)
  b_dms_burden(i)%data=dms_burden(is:ie)
  b_so2_burden(i)%data=so2_burden(is:ie)
  b_so4_burden(i)%data=so4_burden(is:ie)
  b_vso2(i)%data=vso2(is:ie)
  b_dmse(i)%data=dmse(is:ie)
  b_so2e(i)%data=so2e(is:ie)
  b_so4e(i)%data=so4e(is:ie)
  b_bce(i)%data=bce(is:ie)
  b_oce(i)%data=oce(is:ie)
  b_so2dd(i)%data=so2dd(is:ie)
  b_so4dd(i)%data=so4dd(is:ie)
  b_bcdd(i)%data=bcdd(is:ie)
  b_ocdd(i)%data=ocdd(is:ie)
  b_so2wd(i)%data=so2wd(is:ie)
  b_so4wd(i)%data=so4wd(is:ie)
  b_bcwd(i)%data=bcwd(is:ie)
  b_ocwd(i)%data=ocwd(is:ie)
  b_dustwd(i)%data=dustwd(is:ie)
  b_ps(i)%data=ps(is:ie)
  b_zdayfac(i)%data=zdayfac(is:ie)
  b_rlatt(i)%data=rlatt(is:ie)
  b_rlongg(i)%data=rlongg(is:ie)
  b_condc(i)%data=condc(is:ie)
  b_wetfac(i)%data=wetfac(is:ie)
  b_pblh(i)%data=pblh(is:ie)
  b_tss(i)%data=tss(is:ie)
  b_snowd(i)%data=snowd(is:ie)
  b_fg(i)%data=fg(is:ie)
  b_eg(i)%data=eg(is:ie)
  b_u10(i)%data=u10(is:ie)
  b_ustar(i)%data=ustar(is:ie)
  b_zo(i)%data=zo(is:ie)
  b_fracice(i)%data=fracice(is:ie)
  b_sigmf(i)%data=sigmf(is:ie)
  b_cdtq(i)%data=cdtq(is:ie)
  b_so4t(i)%data=so4t(is:ie)

  b_zoxidant(i)%data=zoxidant(is:ie,:)
  b_EMISSFIELD(i)%data=EMISSFIELD(is:ie,:)
  b_erod(i)%data=erod(is:ie,:)
  b_phi_nh(i)%data=phi_nh(is:ie,:)
  b_t(i)%data=t(is:ie,:)
  b_qg(i)%data=qg(is:ie,:)
  b_qlg(i)%data=qlg(is:ie,:)
  b_qfg(i)%data=qfg(is:ie,:)
  b_cfrac(i)%data=cfrac(is:ie,:)
  b_ppfprec(i)%data=ppfprec(is:ie,:)
  b_ppfmelt(i)%data=ppfmelt(is:ie,:)
  b_ppfsnow(i)%data=ppfsnow(is:ie,:)
  b_ppfevap(i)%data=ppfevap(is:ie,:)
  b_ppfsubl(i)%data=ppfsubl(is:ie,:)
  b_pplambs(i)%data=pplambs(is:ie,:)
  b_ppmrate(i)%data=ppmrate(is:ie,:)
  b_ppmaccr(i)%data=ppmaccr(is:ie,:)
  b_ppqfsedice(i)%data=ppqfsedice(is:ie,:)
  b_pprscav(i)%data=pprscav(is:ie,:)
  b_pprfreeze(i)%data=pprfreeze(is:ie,:)
  b_ppfstayice(i)%data=ppfstayice(is:ie,:)
  b_ppfstayliq(i)%data=ppfstayliq(is:ie,:)

  b_oxidantprev(i)%data=oxidantprev(is:ie,:,:)
  b_oxidantnow(i)%data=oxidantnow(is:ie,:,:)
  b_oxidantnext(i)%data=oxidantnext(is:ie,:,:)
  b_xtg(i)%data=xtg(is:ie,:,:)
  b_xtosav(i)%data=xtosav(is:ie,:,:)
  if ( aeromode>=1 ) then
    b_xtg_solub(i)%data=xtg_solub(is:ie,:,:)
  end if
  b_ssn(i)%data=ssn(is:ie,:,:)
end do

call start_log(aero_begin)
!$omp parallel do
do i=1,nb
  call aerocalc_work(i,imax)
end do
call end_log(aero_end)

do i=1,nb
  is=(i-1)*imax+1
  ie=i*imax
  duste(is:ie)=b_duste(i)%data
  dustdd(is:ie)=b_dustdd(i)%data
  dmsso2o(is:ie)=b_dmsso2o(i)%data
  so2so4o(is:ie)=b_so2so4o(i)%data
  dust_burden(is:ie)=b_dust_burden(i)%data
  bc_burden(is:ie)=b_bc_burden(i)%data
  oc_burden(is:ie)=b_oc_burden(i)%data
  dms_burden(is:ie)=b_dms_burden(i)%data
  so2_burden(is:ie)=b_so2_burden(i)%data
  so4_burden(is:ie)=b_so4_burden(i)%data
  dmse(is:ie)=b_dmse(i)%data
  so2e(is:ie)=b_so2e(i)%data
  so4e(is:ie)=b_so4e(i)%data
  bce(is:ie)=b_bce(i)%data
  oce(is:ie)=b_oce(i)%data
  so2dd(is:ie)=b_so2dd(i)%data
  so4dd(is:ie)=b_so4dd(i)%data
  bcdd(is:ie)=b_bcdd(i)%data
  ocdd(is:ie)=b_ocdd(i)%data
  so2wd(is:ie)=b_so2wd(i)%data
  so4wd(is:ie)=b_so4wd(i)%data
  bcwd(is:ie)=b_bcwd(i)%data
  ocwd(is:ie)=b_ocwd(i)%data
  dustwd(is:ie)=b_dustwd(i)%data
  zdayfac(is:ie)=b_zdayfac(i)%data
  so4t(is:ie)=b_so4t(i)%data

  zoxidant(is:ie,:)=b_zoxidant(i)%data

  xtg(is:ie,:,:)=b_xtg(i)%data
  if ( aeromode>=1 ) then
    xtg_solub(is:ie,:,:)=b_xtg_solub(i)%data
  end if
  ssn(is:ie,:,:)=b_ssn(i)%data
end do

end subroutine aerocalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update prognostic aerosols
subroutine aerocalc_work(tile,imax)

use aerosolldr, only : xtg,ssn,aldrcalc,aldrloadoxidant        ! LDR prognostic aerosols
use arrays_m, only : t,ps,qg                                   ! Atmosphere dyamics prognostic arrays
use cc_mpi, only : mydiag                                      ! CC MPI routines
use cfrac_m, only : cfrac                                      ! Cloud fraction
use cloudmod, only : convectivecloudfrac                       ! Prognostic strat cloud
use const_phys, only : grav,rdry                               ! Physical constants
use extraout_m, only : ustar                                   ! Additional diagnostics
use infile, only : getzinp                                     ! Input file routines
use kuocomb_m, only : kbsav,ktsav                              ! JLM convection
use latlong_m, only : rlatt,rlongg                             ! Lat/lon coordinates
use liqwpar_m, only : qlg,qfg                                  ! Cloud water mixing ratios
use morepbl_m, only : pblh,condc,fg,eg                         ! Additional boundary layer diagnostics
use newmpar_m, only : ifull,kl                                 ! Grid parameters
use nharrs_m, only : phi_nh                                    ! Non-hydrostatic atmosphere arrays
use nsibd_m, only : sigmf                                      ! Land-surface arrays
use ozoneread, only : fieldinterpolate                         ! Ozone input routines
use parm_m, only : dt,bpyear,diag, lidjd => idjd               ! Model configuration
use pbl_m, only : tss,cdtq                                     ! Boundary layer arrays
use screen_m, only : u10                                       ! Screen level diagnostics
use sigs_m, only : sig,sigmh,bet,betm,dsig                     ! Atmosphere sigma levels
use soil_m, only : land,so4t                                   ! Soil and surface data
use soilsnow_m, only : snowd,fracice                           ! Soil, snow and surface data
!use soilv_m                                                   ! Soil parameters
!use vegpar_m                                                  ! Vegetation arrays
use work2_m, only : wetfac,zo                                  ! Diagnostic arrays
use zenith_m, only : solargh,zenith                            ! Astronomy routines
use aerodata_m

implicit none

include 'kuocom.h'      ! Convection parameters

integer, intent(in) :: tile,imax
integer jyear,jmonth,jday,jhour,jmin,mins,smins
integer j,k,tt,ttx
integer, parameter :: updateoxidant = 1440 ! update prescribed oxidant fields once per day
real dhr,fjd,sfjd,r1,dlt,alp,slag
real, dimension(imax,kl) :: oxout,zg,clcon,pccw,rhoa
real, dimension(imax,kl) :: tnhs,dz
real, dimension(imax) :: coszro,taudar
real, dimension(imax) :: cldcon,wg
real, dimension(kl+1) :: sigh
logical :: have_idjd
integer :: is,ie,idjd
real, dimension(imax,ilev) :: duma,dumb,dumc

is=(tile-1)*imax+1
ie=tile*imax
if ( lidjd>=(tile-1)*imax+1 .and. lidjd<=tile*imax ) then
  have_idjd=.true.
  idjd=mod(lidjd-1,imax)+1
else
  have_idjd=.false.
  idjd=huge(1)
end if

! timer calculations
call getzinp(fjd,jyear,jmonth,jday,jhour,jmin,mins)
! update prescribed oxidant fields
dhr = dt/3600.
if ( sday(tile)<=mins-updateoxidant ) then
  sday(tile) = mins
  do j = 1,4 
    !we do this as with o3set - I'd prefer to pass tile & imax /is & ie
    duma=b_oxidantprev(tile)%data(:,:,j)
    dumb=b_oxidantnow(tile)%data(:,:,j)
    dumc=b_oxidantnext(tile)%data(:,:,j)
    ! note levels are inverted by fieldinterpolate
    call fieldinterpolate(oxout,duma,dumb,dumc, &
                          rlev,imax,kl,ilev,mins,sig,b_ps(tile)%data(:),interpmeth=0)
    do k = 1,kl
      call aldrloadoxidant(k+(j-1)*kl,oxout(:,k),tile,imax)
    end do
  end do
  ! estimate day length (presumably to preturb day-time OH levels)
  ttx = nint(86400./dt)
  b_zdayfac(tile)%data(:) = 0.
  do tt = ttx,1,-1 ! we seem to get a different answer if dhr=24. and ttx=1.
    smins = int(real(tt-1)*dt/60.)+mins
    sfjd = float(mod( smins, 525600 ))/1440.  ! 525600 = 1440*365
    call solargh(sfjd,bpyear,r1,dlt,alp,slag)
    call zenith(sfjd,r1,dlt,slag,b_rlatt(tile)%data(:),b_rlongg(tile)%data(:),dhr,imax,coszro,taudar)
    where ( taudar>0.5 )
      b_zdayfac(tile)%data(:) = b_zdayfac(tile)%data(:) + 1.
    end where
  end do
  ! final taudar is for current timestep - used to indicate sunlit
  where ( b_zdayfac(tile)%data(:)>0.5 )
    b_zdayfac(tile)%data(:) = real(ttx)/b_zdayfac(tile)%data(:)
  end where
else
  sfjd = float(mod( mins, 525600 ))/1440.  ! 525600 = 1440*365
  call solargh(sfjd,bpyear,r1,dlt,alp,slag)
  call zenith(sfjd,r1,dlt,slag,b_rlatt(tile)%data(:),b_rlongg(tile)%data(:),dhr,imax,coszro,taudar)
  ! taudar is for current timestep - used to indicate sunlit
end if

! set-up input data fields ------------------------------------------------
sigh(1:kl) = sigmh(1:kl) ! store half-levels
sigh(kl+1) = 0.

! Non-hydrostatic terms
tnhs(:,1) = b_phi_nh(tile)%data(:,1)/bet(1)
zg(:,1) = bet(1)*b_t(tile)%data(:,1)/grav
do k = 2,kl
  ! representing non-hydrostatic term as a correction to air temperature
  tnhs(:,k) = (b_phi_nh(tile)%data(:,k)-b_phi_nh(tile)%data(:,k-1)-betm(k)*tnhs(:,k-1))/bet(k)
  zg(:,k) = zg(:,k-1) + (bet(k)*b_t(tile)%data(:,k)+betm(k)*b_t(tile)%data(:,k-1))/grav ! height above surface in meters
end do
do k = 1,kl
  zg(:,k) = zg(:,k) + b_phi_nh(tile)%data(:,k)/grav
  dz(:,k) = -rdry*dsig(k)*(b_t(tile)%data(:,k)+tnhs(:,k))/(grav*sig(k))
  rhoa(:,k) = b_ps(tile)%data(:)*sig(k)/(rdry*b_t(tile)%data(:,k)) ! density of air (kg/m**3)
end do

! estimate convective cloud fraction from leoncld.f
call convectivecloudfrac(clcon,b_condc(tile)%data(:),b_kbsav(tile)%data(:),b_ktsav(tile)%data(:),imax,cldcon=cldcon)
do k = 1,kl
  ! MJT notes - Assume rain for JLM convection
  !where ( k>b_kbsav(tile)%data(:) .and. k<=b_ktsav(tile)%data(:) .and. b_t(tile)%data(:,k)>ticeu )
  !  pccw(:,kl+1-k) = 0.
  where ( k>b_kbsav(tile)%data(:) .and. k<=b_ktsav(tile)%data(:) )
    pccw(:,kl+1-k) = wlc/rhoa(:,k)
  elsewhere
    pccw(:,kl+1-k) = 0.
  end where
end do

! Water converage at surface
wg(:) = min( max( b_wetfac(tile)%data(:), 0. ), 1. )

! MJT notes - We have an option to update the aerosols before the vertical mixing
! or after the vertical mixing.  Updating aerosols before the vertical mixing
! ensures that we can split the convective and non-convective aerosol
! concentrations.  However, updating aerosols after vertical mixing provides a
! better estimate of u10 and pblh.

! update prognostic aerosols
call aldrcalc(dt,sig,zg,dz,wg,b_pblh(tile)%data(:),b_ps(tile)%data(:),b_tss(tile)%data(:),                 &
              b_t(tile)%data(:,:),b_condc(tile)%data(:),b_snowd(tile)%data(:),taudar,b_fg(tile)%data(:),b_eg(tile)%data(:),b_u10(tile)%data(:),b_ustar(tile)%data(:),b_zo(tile)%data(:),     &
              b_land(tile)%data(:),b_fracice(tile)%data(:),b_sigmf(tile)%data(:),b_qg(tile)%data(:,:),b_qlg(tile)%data(:,:),b_qfg(tile)%data(:,:),b_cfrac(tile)%data(:,:),clcon,   &
              cldcon,pccw,rhoa,b_cdtq(tile)%data(:),b_ppfprec(tile)%data(:,:),b_ppfmelt(tile)%data(:,:),       &
              b_ppfsnow(tile)%data(:,:),b_ppfevap(tile)%data(:,:),b_ppfsubl(tile)%data(:,:),b_pplambs(tile)%data(:,:),b_ppmrate(tile)%data(:,:),     &
              b_ppmaccr(tile)%data(:,:),b_ppfstayice(tile)%data(:,:),b_ppfstayliq(tile)%data(:,:),b_ppqfsedice(tile)%data(:,:),    &
              b_pprscav(tile)%data(:,:),b_pprfreeze(tile)%data(:,:),b_zdayfac(tile)%data(:),b_kbsav(tile)%data(:),tile,imax)
              

! store sulfate for LH+SF radiation scheme.  SEA-ESF radiation scheme imports prognostic aerosols in seaesfrad.f90.
! Factor 1.e3 to convert to gS/m2, x 3 to get sulfate from sulfur
b_so4t(tile)%data(:) = 0.
do k = 1,kl
  b_so4t(tile)%data(:) = b_so4t(tile)%data(:) + 3.e3*b_xtg(tile)%data(:,k,3)*rhoa(:,k)*dz(:,k)
enddo

if ( diag .and. mydiag .and. have_idjd ) then
  write(6,*) "tdiag ",b_t(tile)%data(idjd,:)
  write(6,*) "qgdiag ",b_qg(tile)%data(idjd,:)
  write(6,*) "qlgdiag ",b_qlg(tile)%data(idjd,:)
  write(6,*) "qfgdiag ",b_qfg(tile)%data(idjd,:)
  write(6,*) "u10diag ",b_u10(tile)%data(idjd)
  write(6,*) "pblhdiag ",b_pblh(tile)%data(idjd)
  write(6,*) "fracicediag ",b_fracice(tile)%data(idjd)
  write(6,*) "DMSdiag ",b_xtg(tile)%data(idjd,:,1)
  write(6,*) "SO2diag ",b_xtg(tile)%data(idjd,:,2)
  write(6,*) "SO4diag ",b_xtg(tile)%data(idjd,:,3)
  write(6,*) "BCphobdiag ",b_xtg(tile)%data(idjd,:,4)
  write(6,*) "BCphildiag ",b_xtg(tile)%data(idjd,:,5)
  write(6,*) "OCphobdiag ",b_xtg(tile)%data(idjd,:,6)
  write(6,*) "OCphildiag ",b_xtg(tile)%data(idjd,:,7)
  write(6,*) "dust0.8diag ",b_xtg(tile)%data(idjd,:,8)
  write(6,*) "dust1.0diag ",b_xtg(tile)%data(idjd,:,9)
  write(6,*) "dust2.0diag ",b_xtg(tile)%data(idjd,:,10)
  write(6,*) "dust4.0diag ",b_xtg(tile)%data(idjd,:,11)
  write(6,*) "saltfilmdiag ",b_ssn(tile)%data(idjd,:,1)
  write(6,*) "saltjetdiag  ",b_ssn(tile)%data(idjd,:,2)
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
