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
    
! This module interfaces aerosol schemes with CCAM

! Currently, only LDR aerosols are supported.  However, the xtg arrays should be moved to
! here in the future, so that there is a common interface for advection and outcdf.

! - This subroutine assumes only one month at a time is integrated in RCM mode.

module aerointerface

implicit none

private
public load_aerosolldr, aerocalc, aerodrop
public ppfprec, ppfmelt, ppfsnow, ppfevap, ppfsubl, pplambs, ppmrate
public ppmaccr, ppfstayice, ppfstayliq, ppqfsed, pprscav
public opticaldepth

integer, save :: ilon, ilat, ilev
integer, parameter :: naerofamilies = 7      ! Number of aerosol families for optical depth
real, dimension(:,:,:), allocatable, save :: oxidantprev
real, dimension(:,:,:), allocatable, save :: oxidantnow
real, dimension(:,:,:), allocatable, save :: oxidantnext
real, dimension(:,:,:), allocatable, save :: opticaldepth
real, dimension(:,:), allocatable, save :: ppfprec, ppfmelt, ppfsnow           ! data saved from LDR cloud scheme
real, dimension(:,:), allocatable, save :: ppfevap, ppfsubl, pplambs, ppmrate  ! data saved from LDR cloud scheme
real, dimension(:,:), allocatable, save :: ppmaccr, ppqfsed, pprscav           ! data saved from LDR cloud scheme
real, dimension(:,:), allocatable, save :: ppfstayice, ppfstayliq              ! data saved from LDR cloud scheme
real, dimension(:), allocatable, save :: rlev, zdayfac
real, parameter :: wlc = 0.2e-3         ! LWC of deep conv cloud (kg/m**3)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Load aerosols emissions from netcdf
subroutine load_aerosolldr(aerofile, oxidantfile, kdatein)
      
use aerosolldr          ! LDR prognostic aerosols
use cc_mpi              ! CC MPI routines
use infile              ! Input file routines
use ozoneread           ! Ozone input routines
use sigs_m              ! Atmosphere sigma levels
      
implicit none

include 'newmpar.h'     ! Grid parameters
include 'parmgeom.h'    ! Coordinate data
      
integer, intent(in) :: kdatein
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

if ( myid == 0 ) write(6,*) "Initialising prognostic aerosols"

allocate( ppfprec(ifull,kl), ppfmelt(ifull,kl) )
allocate( ppfsnow(ifull,kl) )
allocate( ppfevap(ifull,kl), ppfsubl(ifull,kl) )
allocate( pplambs(ifull,kl), ppmrate(ifull,kl) )
allocate( ppmaccr(ifull,kl) )
allocate( ppqfsed(ifull,kl), pprscav(ifull,kl) )
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
ppqfsed = 0.
pprscav = 0.
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

return
end subroutine load_aerosolldr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update prognostic aerosols
subroutine aerocalc

use aerosolldr           ! LDR prognostic aerosols
use arrays_m             ! Atmosphere dyamics prognostic arrays
use cc_mpi               ! CC MPI routines
use cfrac_m              ! Cloud fraction
use cloudmod             ! Prognostic strat cloud
use extraout_m           ! Additional diagnostics
use infile               ! Input file routines
use kuocomb_m            ! JLM convection
use latlong_m            ! Lat/lon coordinates
use liqwpar_m            ! Cloud water mixing ratios
use morepbl_m            ! Additional boundary layer diagnostics
use nharrs_m             ! Non-hydrostatic atmosphere arrays
use nsibd_m              ! Land-surface arrays
use ozoneread            ! Ozone input routines
use pbl_m                ! Boundary layer arrays
use screen_m             ! Screen level diagnostics
use sigs_m               ! Atmosphere sigma levels
use soil_m               ! Soil and surface data
use soilsnow_m           ! Soil, snow and surface data
use tkeeps, only : zidry ! TKE-EPS boundary layer
use vegpar_m             ! Vegetation arrays
use work2_m              ! Diagnostic arrays
use zenith_m             ! Astronomy routines

implicit none

include 'newmpar.h'     ! Grid parameters
include 'const_phys.h'  ! Physical constants
include 'kuocom.h'      ! Convection parameters
include 'parm.h'        ! Model configuration
include 'soilv.h'       ! Soil parameters

integer jyear,jmonth,jday,jhour,jmin,mins,smins
integer j,k,tt,ttx
integer, save :: sday=-9999
integer, parameter :: updateoxidant = 1440 ! update prescribed oxidant fields once per day
real dhr,fjd,sfjd,r1,dlt,alp,slag
real, dimension(ifull,kl) :: oxout,zg,clcon,pccw,rhoa
real, dimension(ifull,kl) :: tnhs,dz,tv
real, dimension(ifull) :: coszro,taudar
real, dimension(ifull) :: cldcon,wg,pblx
real, dimension(kl+1) :: sigh

! timer calculations
call getzinp(fjd,jyear,jmonth,jday,jhour,jmin,mins)
! update prescribed oxidant fields
dhr = dt/3600.
if ( sday<=mins-updateoxidant ) then
  sday = mins
  do j = 1,4 
    ! note levels are inverted
    call fieldinterpolate(oxout,oxidantprev(:,:,j),oxidantnow(:,:,j),oxidantnext(:,:,j), &
                          rlev,ifull,kl,ilev,mins,sig,ps,interpmeth=0)
    do k = 1,kl
      call aldrloadoxidant(k+(j-1)*kl,oxout(:,k))
    end do
  end do
  ! estimate day length (presumably to preturb day-time OH levels)
  ttx = nint(86400./dt)
  zdayfac(:) = 0.
  do tt = ttx,1,-1 ! we seem to get a different answer if dhr=24. and ttx=1.
    smins = int(real(tt-1)*dt/60.)+mins
    sfjd = float(mod( smins, 525600 ))/1440.  ! 525600 = 1440*365
    call solargh(sfjd,bpyear,r1,dlt,alp,slag)
    call zenith(sfjd,r1,dlt,slag,rlatt,rlongg,dhr,ifull,coszro,taudar)
    where ( taudar>0.5 )
      zdayfac(:) = zdayfac(:) + 1.
    end where
  end do
  ! final taudar is for current timestep - used to indicate sunlit
  where ( zdayfac>1.e-20 )
    zdayfac(:) = real(ttx)/zdayfac(:)
  end where
else
  sfjd = float(mod( mins, 525600 ))/1440.  ! 525600 = 1440*365
  call solargh(sfjd,bpyear,r1,dlt,alp,slag)
  call zenith(sfjd,r1,dlt,slag,rlatt,rlongg,dhr,ifull,coszro,taudar)
  ! taudar is for current timestep - used to indicate sunlit
end if

! set-up input data fields ------------------------------------------------
sigh(1:kl) = sigmh(1:kl) ! store half-levels
sigh(kl+1) = 0.

! Non-hydrostatic terms
tv(:,:) = t(1:ifull,:)*(1.+0.61*qg(1:ifull,:)-qlg(1:ifull,:)-qfg(1:ifull,:) &
                       -qrg(1:ifull,:)-qsng(1:ifull,:)-qgrg(1:ifull,:))
tnhs(:,1) = phi_nh(:,1)/bet(1)
zg(:,1) = bet(1)*tv(:,1)/grav
do k = 2,kl
  ! representing non-hydrostatic term as a correction to air temperature
  tnhs(:,k) = (phi_nh(:,k)-phi_nh(:,k-1)-betm(k)*tnhs(:,k-1))/bet(k)
  zg(:,k) = zg(:,k-1) + (bet(k)*tv(:,k)+betm(k)*tv(:,k-1))/grav ! height above surface in meters
end do
do k = 1,kl
  zg(:,k) = zg(:,k) + phi_nh(:,k)/grav
  dz(:,k) = -rdry*dsig(k)*(tv(:,k)+tnhs(:,k))/(grav*sig(k))
  rhoa(:,k) = ps(1:ifull)*sig(k)/(rdry*tv(1:ifull,k)) ! density of air (kg/m**3)
end do

! estimate convective cloud fraction from leoncld.f
call convectivecloudfrac(clcon,cldcon=cldcon)
do k = 1,kl
  ! MJT notes - Assume rain for JLM convection
  !where ( k>kbsav .and. k<=ktsave .and. t(1:ifull,k)>ticeu )
  !  pccw(:,kl+1-k) = 0.
  where ( k>kbsav .and. k<=ktsav )
    pccw(:,kl+1-k) = wlc/rhoa(:,k)
  elsewhere
    pccw(:,kl+1-k) = 0.
  end where
end do

! Water converage at surface
wg(:) = min( max( wetfac, 0. ), 1. )

! Define boundary layer height used by diagnostic sea-salt
if ( nvmix==6 .and. nlocal==7 ) then
  ! TKE
  pblx(:) = zidry(:) ! Dry convective boundary layer height
else
  ! local Ri
  pblx(:) = pblh(:)
end if

! MJT notes - We have an option to update the aerosols before the vertical mixing
! or after the vertical mixing.  Updating aerosols before the vertical mixing
! ensures that we can split the convective and non-convective aerosol
! concentrations.  However, updating aerosols after vertical mixing provides a
! better estimate of u10 and pblh.

! update prognostic aerosols
call aldrcalc(dt,sig,zg,dz,wg,pblx,ps,tss,                 &
              t,condc,snowd,taudar,fg,eg,u10,ustar,zo,     &
              land,fracice,sigmf,qg,qlg,qfg,cfrac,clcon,   &
              cldcon,pccw,rhoa,cdtq,ppfprec,ppfmelt,       &
              ppfsnow,ppfevap,ppfsubl,pplambs,ppmrate,     &
              ppmaccr,ppfstayice,ppfstayliq,ppqfsed,       &
              pprscav,zdayfac,kbsav)
              

! store sulfate for LH+SF radiation scheme.  SEA-ESF radiation scheme imports prognostic aerosols in seaesfrad.f90.
! Factor 1.e3 to convert to gS/m2, x 3 to get sulfate from sulfur
so4t(:) = 0.
do k = 1,kl
  so4t(:) = so4t(:) + 3.e3*xtg(1:ifull,k,3)*rhoa(:,k)*dz(:,k)
enddo

if ( diag .and. mydiag ) then
  write(6,*) "tdiag ",t(idjd,:)
  write(6,*) "qgdiag ",qg(idjd,:)
  write(6,*) "qlgdiag ",qlg(idjd,:)
  write(6,*) "qfgdiag ",qfg(idjd,:)
  write(6,*) "u10diag ",u10(idjd)
  write(6,*) "pblxdiag ",pblx(idjd)
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
end subroutine aerocalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Estimate cloud droplet size
subroutine aerodrop(istart,imax,cdn,rhoa,outconv)

use aerosolldr              ! LDR prognostic aerosols
use latlong_m, only : rlatt ! Lat/lon coordinates
use soil_m, only : land     ! Soil and surface data

implicit none

include 'newmpar.h'         ! Grid parameters
include 'const_phys.h'      ! Physical constants
include 'cparams.h'         ! Input cloud scheme parameters
include 'parm.h'            ! Model configuration

integer, intent(in) :: istart,imax
integer k,indirmode,iend
real, dimension(imax,kl), intent(out) :: cdn
real, dimension(imax,kl), intent(in) :: rhoa
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
  case(2)
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
