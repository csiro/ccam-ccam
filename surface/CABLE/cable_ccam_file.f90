! Conformal Cubic Atmospheric Model
    
! Copyright 2026 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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

! File subroutines for CABLE-CCAM interface

module cable_ccam_file

use cable_ccam_common

implicit none

private
public loadtile, defaulttile, savetiledef, savetile

contains

! *************************************************************************************  
! This subroutine loads CABLE tile data
subroutine defaulttile

use cc_mpi
use darcdf_m
use infile
use newmpar_m
use parm_m
use pbl_m
use soil_m
use soilsnow_m
use vegpar_m
  
integer k, tile, js, je
real, dimension(imax) :: dummy_pack

do tile = 1,ntiles
  js = (tile-1)*imax + 1
  je = tile*imax
  if ( tdata(tile)%mp>0 ) then
      
    do k = 1,cbm_ms
      call cable_pack(tgg(:,k),ssnow(tile)%tgg(:,k),tile)
      call cable_pack(wb(:,k),ssnow(tile)%wb(:,k),tile)
      call cable_pack(wbice(:,k),ssnow(tile)%wbice(:,k),tile)
    end do
    do k = 1,3
      call cable_pack(tggsn(:,k),ssnow(tile)%tggsn(:,k),tile)
      call cable_pack(smass(:,k),ssnow(tile)%smass(:,k),tile)
      call cable_pack(ssdn(:,k),ssnow(tile)%ssdn(:,k),tile)
      dummy_pack = smass(js:je,k)/ssdn(js:je,k)
      call cable_pack(dummy_pack,ssnow(tile)%sdepth(:,k),tile)
      ssnow(tile)%sconds(:,k) = 0.3_8
    end do      
    call cable_pack(tss(:),rad(tile)%trad(:),tile)
    call cable_pack(ssdnn(:),ssnow(tile)%ssdnn,tile)
    call cable_pack(isflag(:),ssnow(tile)%isflag,tile)
    call cable_pack(snowd(:),ssnow(tile)%snowd,tile)
    call cable_pack(snage(:),ssnow(tile)%snage,tile)
    ssnow(tile)%Qrecharge = 0._8
    canopy(tile)%sublayer_dz = 0._8
    ssnow(tile)%rtevap_sat = 0._8
    ssnow(tile)%rtevap_unsat = 0._8
    ssnow(tile)%satfrac = 0.5_8
    ssnow(tile)%wbliq = ssnow(tile)%wb - ssnow(tile)%wbice
    ssnow(tile)%GWwb = 0.5_8*soil(tile)%ssat 
    ssnow(tile)%wtd = 20000._8
    dummy_pack = real(1-isflag(js:je))*tgg(js:je,1) &
               + real(isflag(js:je))*tggsn(js:je,1) - 273.15
    call cable_pack(dummy_pack,ssnow(tile)%tsurface,tile)
    ssnow(tile)%rtsoil = 50._8
    canopy(tile)%cansto = 0._8
    canopy(tile)%ga = 0._8
    canopy(tile)%us = 0.01_8
    ssnow(tile)%pudsto = 0._8
    ssnow(tile)%wetfac = 0._8
    ssnow(tile)%osnowd = ssnow(tile)%snowd
    canopy(tile)%fhs_cor = 0._8
    canopy(tile)%fes_cor = 0._8
    canopy(tile)%fns_cor = 0._8
    canopy(tile)%ga_cor = 0._8
  
    ! default value for fwsoil.  Recaculated by cable_canopy or by SLI
    canopy(tile)%fwsoil = max( 1.e-9_8, sum( veg(tile)%froot*max(1.e-9_8, &
        min(1._8,ssnow(tile)%wb-spread(soil(tile)%swilt,2,cbm_ms))),2)    &
        / ( soil(tile)%sfc-soil(tile)%swilt ) )

  end if ! mp>0
end do   ! tile 
  
if ( mp_global>0 ) then  
  call defaulttile_sli
  call fixtile
end if

return
end subroutine defaulttile

subroutine defaulttile_sli

use cc_mpi
use darcdf_m
use infile
use newmpar_m
use parm_m
use soil_m
use soilsnow_m
use vegpar_m
  
integer k, tile

if ( soil_struc==1 ) then  

  do tile = 1,ntiles
    if ( tdata(tile)%mp>0 ) then
    
      ssnow(tile)%h0 = 0._8      ! pond height
      ssnow(tile)%snowliq = 0._8 ! liquid snow
      ssnow(tile)%Tsoil = ssnow(tile)%tgg - 273.15_8
      ssnow(tile)%thetai = ssnow(tile)%wbice
      do k = 1,cbm_ms
        ssnow(tile)%S(:,k) = ssnow(tile)%wb(:,k)/soil(tile)%ssat
      end do
      where ( ssnow(tile)%snowd>0. )
        ssnow(tile)%nsnow = 1
        ssnow(tile)%sdepth(:,1) = ssnow(tile)%snowd
        ssnow(tile)%smass(:,1) = ssnow(tile)%snowd*ssnow(tile)%ssdn(:,1)
      elsewhere
        ssnow(tile)%nsnow = 0
        ssnow(tile)%sdepth(:,1) = 0._8
        ssnow(tile)%smass(:,1) = 0._8
      end where
      !ssnow(tile)%nsnow = 0
      !ssnow(tile)%snowd = 0._8
      
    end if ! mp>0
  end do   ! tile

  !if ( mp_global>0 ) then
  !  call fixtile
  !end if

end if ! soil_struc==-1

return
end subroutine defaulttile_sli

subroutine loadtile(usedefault)

use cc_mpi
use darcdf_m
use infile
use newmpar_m
use nsibd_m, only : sigmf, carb_plant, carb_litter, carb_soil
use parm_m
use pbl_m
use soil_m
use soilsnow_m
use vegpar_m
  
logical, intent(in), optional :: usedefault
integer k, n, ierr, idv, ierr_casa, ierr_sli, ierr_pop, ierr_svs, ierr_cvc
integer jyear,jmonth,jday,jhour,jmin,mins, ll, cc, hh, dd
integer np_pop, iq, m
integer tile
integer, dimension(6) :: ierr_check
integer, dimension(ifull) :: dati
integer, dimension(ifull,maxtile) :: nmp
integer, dimension(:), allocatable :: dati_out
real, dimension(ifull) :: datr
real(kind=8), dimension(mp_max) :: dummy_unpack
real(kind=8), dimension(ifull) :: dat
real(kind=8), dimension(ifull,cbm_ms) :: datms
real(kind=8), dimension(ifull,3) :: dat3
real(kind=8), dimension(ifull,mplant) :: datmplant
real(kind=8), dimension(ifull,mlitter) :: datmlitter
real(kind=8), dimension(ifull,msoil) :: datmsoil
real(kind=8), dimension(:), allocatable :: dat_out
real(kind=8), dimension(:,:), allocatable :: datpatch
real(kind=8), dimension(:,:), allocatable :: datage
real(kind=8), dimension(:,:,:), allocatable :: datpc
logical tst
logical defaultmode
character(len=80) vname
character(len=21) testname

if ( myid==0 ) write(6,*) 'Read CABLE and CASA initial conditions'

! force CABLE to use generic input for all tiles
! if usedefault = defaultmode = .true.
defaultmode = .false.
if ( present(usedefault) ) then
  defaultmode = usedefault
end if

! check that CABLE data exists in restart file
! and communicate the result to all processors
! as not all processors are assigned an input file
ierr = 1
ierr_casa = 1
ierr_sli = 1
ierr_pop = 1
ierr_svs = 1
ierr_cvc = 1

! io_in==1 ensures no interpolation is required
if ( io_in==1 .and. .not.defaultmode ) then
  if ( myid==0 .or. pfall ) then
    write(testname,'("t",I1.1,"_tgg1")') maxtile  
    call ccnf_inq_varid(ncid,testname,idv,tst)
    if ( .not.tst ) then
      ierr = 0
    end if
    write(testname,'("t",I1.1,"_cplant1")') maxtile  
    call ccnf_inq_varid(ncid,testname,idv,tst)
    if ( .not.tst ) then
      ierr_casa = 0
    end if
    write(testname,'("t",I1.1,"_hzero")') maxtile  
    call ccnf_inq_varid(ncid,testname,idv,tst)
    if ( .not.tst ) then
      ierr_sli = 0
    end if
    write(testname,'("t",I1.1,"_pop_grid_cmass_sum")') maxtile  
    call ccnf_inq_varid(ncid,testname,idv,tst)
    if ( .not.tst ) then
      ierr_pop = 0
    end if
    write(testname,'("t",I1.1,"_svs")') maxtile  
    call ccnf_inq_varid(ncid,testname,idv,tst)
    if ( .not.tst ) then
      ierr_svs = 0
    end if
    write(testname,'("t",I1.1,"_cvc")') maxtile  
    call ccnf_inq_varid(ncid,testname,idv,tst)
    if ( .not.tst ) then
      ierr_cvc = 0
    end if
  end if
end if

! Communicate with processes if not all processes are reading the input file.
if ( .not.pfall ) then
  ierr_check(1) = ierr
  ierr_check(2) = ierr_casa
  ierr_check(3) = ierr_sli
  ierr_check(4) = ierr_pop
  ierr_check(5) = ierr_svs
  ierr_check(6) = ierr_cvc
  call ccmpi_bcast(ierr_check(1:6),0,comm_world)
  ierr       = ierr_check(1)
  ierr_casa  = ierr_check(2)
  ierr_sli   = ierr_check(3)
  ierr_pop   = ierr_check(4)
  ierr_svs   = ierr_check(5)
  ierr_cvc   = ierr_check(6)
end if

if ( myid==0 ) then
  write(6,*) "-> Found ierr,ierr_casa_ierr_sli ",ierr,ierr_casa,ierr_sli
  write(6,*) "->    ierr_pop,ierr_svs,ierr_cvc ",ierr_pop,ierr_svs,ierr_cvc
end if
  
call defaulttile ! initially use default values before overwriting

! default
do tile = 1,ntiles
  if ( tdata(tile)%mp>0 ) then
    tdata(tile)%old_sv = tdata(tile)%sv
    tdata(tile)%old_cveg = tdata(tile)%cveg
  end if
end do

! check for changes
if ( ierr_cvc==0 ) then
  do n = 1,maxtile      
    write(vname,'("t",I1.1,"_cvc")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    dati = nint(dat)
    do tile = 1,ntiles
      call cable_pack(dati,tdata(tile)%old_cveg,tile,nb=n)
    end do  
  end do
end if
call create_new_tile_map(nmp)

if ( ierr/=0 ) then
  ! Cannot locate tile data, use diagnostic data instead    
  if ( myid==0 ) write(6,*) "-> Use gridbox averaged data to initialise CABLE"
else
  ! read tile data
  if ( myid==0 ) write(6,*) "-> Use tiled data to initialise CABLE"  
  do n = 1,maxtile
    if ( ierr_svs == 0 ) then
      write(vname,'("t",I1.1,"_svs")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      datr = real( dat )
      do tile = 1,ntiles
        call cable_pack(datr,tdata(tile)%old_sv,tile,nmp(:,n))
      end do  
    end if        
    write(vname,'("t",I1.1,"_tgg")') n
    call histrd(iarchi-1,ierr,vname,datms(:,1:cbm_ms),ifull)
    do k = 1,cbm_ms
      do iq = 1,ifull
        if ( land(iq) .and. (datms(iq,k)<100._8.or.datms(iq,k)>400._8) ) then
          ! change in land-sea mask?
          datms(iq,k) = tss(iq) ! use surface temperature
        end if
      end do  
      do tile = 1,ntiles
        call cable_pack(datms(:,k),ssnow(tile)%tgg(:,k),tile,nmp(:,n))
      end do  
    end do
    write(vname,'("t",I1.1,"_wb")') n
    call histrd(iarchi-1,ierr,vname,datms(:,1:cbm_ms),ifull)
    do k = 1,cbm_ms
      do tile = 1,ntiles
        call cable_pack(datms(:,k),ssnow(tile)%wb(:,k),tile,nmp(:,n))
      end do  
    end do
    write(vname,'("t",I1.1,"_wbice")') n
    call histrd(iarchi-1,ierr,vname,datms(:,1:cbm_ms),ifull)
    do k = 1,cbm_ms
      do tile = 1,ntiles
        call cable_pack(datms(:,k),ssnow(tile)%wbice(:,k),tile,nmp(:,n))
      end do  
    end do
    write(vname,'("t",I1.1,"_tggsn")') n
    call histrd(iarchi-1,ierr,vname,dat3(:,1:3),ifull)
    do k = 1,3
      do tile = 1,ntiles
        call cable_pack(dat3(:,k),ssnow(tile)%tggsn(:,k),tile,nmp(:,n))
      end do  
    end do
    write(vname,'("t",I1.1,"_smass")') n
    call histrd(iarchi-1,ierr,vname,dat3(:,1:3),ifull)
    do k = 1,3
      do tile = 1,ntiles
        call cable_pack(dat3(:,k),ssnow(tile)%smass(:,k),tile,nmp(:,n))
      end do  
    end do
    write(vname,'("t",I1.1,"_ssdn")') n
    call histrd(iarchi-1,ierr,vname,dat3(:,1:3),ifull)
    do k = 1,3
      do tile = 1,ntiles
        call cable_pack(dat3(:,k),ssnow(tile)%ssdn(:,k),tile,nmp(:,n))
      end do  
    end do
    write(vname,'("t",I1.1,"_sdepth",I1.1)') n
    call histrd(iarchi-1,ierr,vname,dat3(:,1:3),ifull)
    do k = 1,3
      do tile = 1,ntiles
        call cable_pack(dat3(:,k),ssnow(tile)%sdepth(:,k),tile,nmp(:,n))
      end do
    end do
    write(vname,'("t",I1.1,"_sconds")') n
    call histrd(iarchi-1,ierr,vname,dat3(:,1:3),ifull)
    do k = 1,3
      do tile = 1,ntiles
        call cable_pack(dat3(:,k),ssnow(tile)%sconds(:,k),tile,nmp(:,n))
      end do  
    end do
    write(vname,'("t",I1.1,"_ssdnn")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    do tile = 1,ntiles
      call cable_pack(dat,ssnow(tile)%ssdnn(:),tile,nmp(:,n))
    end do  
    write(vname,'("t",I1.1,"_sflag")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    dati = nint(dat)
    do tile = 1,ntiles
      call cable_pack(dati,ssnow(tile)%isflag(:),tile,nmp(:,n))
    end do  
    write(vname,'("t",I1.1,"_snd")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    do tile = 1,ntiles
      call cable_pack(dat,ssnow(tile)%snowd(:),tile,nmp(:,n))
    end do  
    write(vname,'("t",I1.1,"_osnd")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    do tile = 1,ntiles
      call cable_pack(dat,ssnow(tile)%osnowd(:),tile,nmp(:,n))
    end do  
    write(vname,'("t",I1.1,"_snage")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    do tile = 1,ntiles
      call cable_pack(dat,ssnow(tile)%snage(:),tile,nmp(:,n))
    end do
    write(vname,'("t",I1.1,"_rtsoil")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    do tile = 1,ntiles
      call cable_pack(dat,ssnow(tile)%rtsoil(:),tile,nmp(:,n))
    end do  
    write(vname,'("t",I1.1,"_GWwb")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    do tile = 1,ntiles
      call cable_pack(dat,ssnow(tile)%GWwb(:),tile,nmp(:,n))
    end do  
    write(vname,'("t",I1.1,"_wtd")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    do tile = 1,ntiles
      call cable_pack(dat,ssnow(tile)%wtd(:),tile,nmp(:,n))
    end do  
    write(vname,'("t",I1.1,"_cansto")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    do tile = 1,ntiles
      call cable_pack(dat,canopy(tile)%cansto(:),tile,nmp(:,n))
    end do  
    write(vname,'("t",I1.1,"_us")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    do tile = 1,ntiles
      call cable_pack(dat,canopy(tile)%us(:),tile,nmp(:,n))
    end do  
    write(vname,'("t",I1.1,"_pudsto")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    do tile = 1,ntiles
      call cable_pack(dat,ssnow(tile)%pudsto(:),tile,nmp(:,n))
    end do
    write(vname,'("t",I1.1,"_wetfac")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    do tile = 1,ntiles
      call cable_pack(dat,ssnow(tile)%wetfac(:),tile,nmp(:,n))
    end do  
    write(vname,'("t",I1.1,"_ga")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    do tile = 1,ntiles
      call cable_pack(dat,canopy(tile)%ga(:),tile,nmp(:,n))
    end do  
  end do
  
  ! soil temperature check
  do tile = 1,ntiles
    if ( tdata(tile)%mp>0 ) then
      if ( any(ssnow(tile)%tgg>400.) ) then
        write(6,*) "ERROR: Invalid CABLE temperature when reading tile"
        write(6,*) "ssnow%tgg ",maxval(ssnow(tile)%tgg)
        stop -1
      end if
    end if
  end do  

end if ! ierr/=0 ..else..
  
if ( soil_struc==1 ) then
  if ( ierr_sli/=0 ) then
    if ( myid==0 ) write(6,*) "-> Use gridbox averaged data to initialise SLI"
  else 
    if ( myid==0 ) write(6,*) "-> Use tiled data to initialise SLI"  
    do n = 1,maxtile
      write(vname,'("t",I1.1,"_hzero")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        call cable_pack(dat,ssnow(tile)%h0(:),tile,nmp(:,n))
      end do  
      write(vname,'("t",I1.1,"_s")') n
      call histrd(iarchi-1,ierr,vname,datms(:,1:cbm_ms),ifull)
      do k = 1,cbm_ms
        do tile = 1,ntiles
          call cable_pack(datms(:,k),ssnow(tile)%S(:,k),tile,nmp(:,n))
        end do  
      end do
      write(vname,'("t",I1.1,"_tsoil")') n
      call histrd(iarchi-1,ierr,vname,datms(:,1:cbm_ms),ifull)
      do k = 1,cbm_ms
        do tile = 1,ntiles 
          call cable_pack(datms(:,k),ssnow(tile)%tsoil(:,k),tile,nmp(:,n))
        end do  
      end do
      write(vname,'("t",I1.1,"_thetai")') n
      call histrd(iarchi-1,ierr,vname,datms(:,1:cbm_ms),ifull)
      do k = 1,cbm_ms
        do tile = 1,ntiles  
          call cable_pack(datms(:,k),ssnow(tile)%thetai(:,k),tile,nmp(:,n))
        end do  
      end do
      write(vname,'("t",I1.1,"_snowliq",I1.1)') n,1
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        call cable_pack(dat,ssnow(tile)%snowliq(:,1),tile,nmp(:,n)) ! currently nsnow_max=1
      end do  
      write(vname,'("t",I1.1,"_tsurface")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        call cable_pack(dat,ssnow(tile)%tsurface(:),tile,nmp(:,n))
      end do  
      write(vname,'("t",I1.1,"_nsnow")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      dati = nint(dat)
      do tile = 1,ntiles
        call cable_pack(dati,ssnow(tile)%nsnow(:),tile,nmp(:,n))
      end do  
      write(vname,'("t",I1.1,"_fwsoil")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        call cable_pack(dat,canopy(tile)%fwsoil(:),tile,nmp(:,n))
      end do
    end do
  end if ! ierr_sli/=0 ..else..
end if   ! soil_struc==1
  
if ( ccycle/=0 ) then
  if ( ierr_casa/=0 ) then
    if ( .not.allocated( carb_plant ) ) then
      if ( myid==0 ) write(6,*) "-> Use default data to initialise CASA-CNP"  
    else  
      if ( myid==0 ) write(6,*) "-> Use interpolated tiled data to initialise CASA-CNP"
      do m = 1,10
        do tile = 1,ntiles
          call loadtile_carbonpools(carb_plant(:,:,1,m),casapool(tile)%cplant(:,:),tile,m)
          call loadtile_carbonpools(carb_plant(:,:,2,m),casapool(tile)%nplant(:,:),tile,m)
          call loadtile_carbonpools(carb_plant(:,:,3,m),casapool(tile)%pplant(:,:),tile,m)
          call loadtile_carbonpools(carb_litter(:,:,1,m),casapool(tile)%clitter(:,:),tile,m)
          call loadtile_carbonpools(carb_litter(:,:,2,m),casapool(tile)%nlitter(:,:),tile,m)
          call loadtile_carbonpools(carb_litter(:,:,3,m),casapool(tile)%plitter(:,:),tile,m)
          call loadtile_carbonpools(carb_soil(:,:,1,m),casapool(tile)%csoil(:,:),tile,m)
          call loadtile_carbonpools(carb_soil(:,:,2,m),casapool(tile)%nsoil(:,:),tile,m)
          call loadtile_carbonpools(carb_soil(:,:,3,m),casapool(tile)%psoil(:,:),tile,m)
        end do  
      end do   ! mm = 1,10
      m = 14 ! use index 11 to store pft 14
      do tile = 1,ntiles
        call loadtile_carbonpools(carb_plant(:,:,1,11),casapool(tile)%cplant(:,:),tile,m)
        call loadtile_carbonpools(carb_plant(:,:,2,11),casapool(tile)%nplant(:,:),tile,m)
        call loadtile_carbonpools(carb_plant(:,:,3,11),casapool(tile)%pplant(:,:),tile,m)
        call loadtile_carbonpools(carb_litter(:,:,1,11),casapool(tile)%clitter(:,:),tile,m)
        call loadtile_carbonpools(carb_litter(:,:,2,11),casapool(tile)%nlitter(:,:),tile,m)
        call loadtile_carbonpools(carb_litter(:,:,3,11),casapool(tile)%plitter(:,:),tile,m)
        call loadtile_carbonpools(carb_soil(:,:,1,11),casapool(tile)%csoil(:,:),tile,m)
        call loadtile_carbonpools(carb_soil(:,:,2,11),casapool(tile)%nsoil(:,:),tile,m)
        call loadtile_carbonpools(carb_soil(:,:,3,11),casapool(tile)%psoil(:,:),tile,m)
      end do  
      deallocate( carb_plant, carb_litter, carb_soil )
    end if  
  else
    if ( myid==0 ) write(6,*) "-> Use tiled data to initialise CASA-CNP"  
    do n = 1,maxtile
      write(vname,'("t",I1.1,"_cplant")') n
      call histrd(iarchi-1,ierr,vname,datmplant(:,1:mplant),ifull)
      do k = 1,mplant
        do tile = 1,ntiles 
          call cable_pack(datmplant(:,k),casapool(tile)%cplant(:,k),tile,nmp(:,n))
        end do
      end do
      write(vname,'("t",I1.1,"_nplant")') n
      call histrd(iarchi-1,ierr,vname,datmplant(:,1:mplant),ifull)
      do k = 1,mplant
        do tile = 1,ntiles 
          call cable_pack(datmplant(:,k),casapool(tile)%nplant(:,k),tile,nmp(:,n))
        end do  
      end do
      write(vname,'("t",I1.1,"_pplant")') n
      call histrd(iarchi-1,ierr,vname,datmplant(:,1:mplant),ifull)
      do k = 1,mplant
        do tile = 1,ntiles  
          call cable_pack(datmplant(:,k),casapool(tile)%pplant(:,k),tile,nmp(:,n))
        end do  
      end do
      write(vname,'("t",I1.1,"_clitter")') n
      call histrd(iarchi-1,ierr,vname,datmlitter(:,1:mlitter),ifull)
      do k = 1,mlitter
        do tile = 1,ntiles  
          call cable_pack(datmlitter(:,k),casapool(tile)%clitter(:,k),tile,nmp(:,n))
        end do  
      end do
      write(vname,'("t",I1.1,"_nlitter")') n
      call histrd(iarchi-1,ierr,vname,datmlitter(:,1:mlitter),ifull)
      do k = 1,mlitter
        do tile = 1,ntiles  
          call cable_pack(datmlitter(:,k),casapool(tile)%nlitter(:,k),tile,nmp(:,n))
        end do  
      end do
      write(vname,'("t",I1.1,"_plitter")') n
      call histrd(iarchi-1,ierr,vname,datmlitter(:,1:mlitter),ifull)
      do k = 1,mlitter
        do tile = 1,ntiles  
          call cable_pack(datmlitter(:,k),casapool(tile)%plitter(:,k),tile,nmp(:,n))
        end do  
      end do
      write(vname,'("t",I1.1,"_csoil")') n
      call histrd(iarchi-1,ierr,vname,datmsoil(:,1:msoil),ifull)
      do k = 1,msoil
        do tile = 1,ntiles  
          call cable_pack(datmsoil(:,k),casapool(tile)%csoil(:,k),tile,nmp(:,n))
        end do  
      end do
      write(vname,'("t",I1.1,"_nsoil")') n
      call histrd(iarchi-1,ierr,vname,datmsoil(:,1:msoil),ifull)
      do k = 1,msoil
        do tile = 1,ntiles  
          call cable_pack(datmsoil(:,k),casapool(tile)%nsoil(:,k),tile,nmp(:,n))
        end do  
      end do
      write(vname,'("t",I1.1,"_psoil")') n
      call histrd(iarchi-1,ierr,vname,datmsoil(:,1:msoil),ifull)
      do k = 1,msoil
        do tile = 1,ntiles  
          call cable_pack(datmsoil(:,k),casapool(tile)%psoil(:,k),tile,nmp(:,n))
        end do  
      end do
      write(vname,'("t",I1.1,"_glai")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        call cable_pack(dat,casamet(tile)%glai,tile,nmp(:,n))
      end do  
      write(vname,'("t",I1.1,"_phen")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        call cable_pack(dat,phen(tile)%phen,tile,nmp(:,n))
      end do  
      write(vname,'("t",I1.1,"_aphen")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        call cable_pack(dat,phen(tile)%aphen,tile,nmp(:,n))
      end do  
      write(vname,'("t",I1.1,"_phenphase")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      dati = nint(dat)
      do tile = 1,ntiles
        call cable_pack(dati,phen(tile)%phase,tile,nmp(:,n))
      end do
      write(vname,'("t",I1.1,"_doyphase3")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      dati = nint(dat)
      do tile = 1,ntiles
        call cable_pack(dati,phen(tile)%doyphase(:,3),tile,nmp(:,n))
      end do  
      write(vname,'("t",I1.1,"_clabile")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        call cable_pack(dat,casapool(tile)%clabile,tile,nmp(:,n))
      end do  
      write(vname,'("t",I1.1,"_nsoilmin")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        call cable_pack(dat,casapool(tile)%nsoilmin,tile,nmp(:,n))
      end do  
      write(vname,'("t",I1.1,"_psoillab")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        call cable_pack(dat,casapool(tile)%psoillab,tile,nmp(:,n))
      end do  
      write(vname,'("t",I1.1,"_psoilsorb")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        call cable_pack(dat,casapool(tile)%psoilsorb,tile,nmp(:,n))
      end do  
      write(vname,'("t",I1.1,"_psoilocc")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        call cable_pack(dat,casapool(tile)%psoilocc,tile,nmp(:,n))
      end do  
      write(vname,'("t",I1.1,"_crmplant")') n
      call histrd(iarchi-1,ierr,vname,datmplant(:,1:mplant),ifull)
      do k = 1,mplant
        do tile = 1,ntiles
          call cable_pack(datmplant(:,k),casaflux(tile)%crmplant(:,k),tile,nmp(:,n))
        end do  
      end do
      write(vname,'("t",I1.1,"_fracsapwood")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        call cable_pack(dat,casaflux(tile)%frac_sapwood,tile,nmp(:,n))
      end do  
      write(vname,'("t",I1.1,"_sapwoodarea")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        call cable_pack(dat,casaflux(tile)%sapwood_area,tile,nmp(:,n))
      end do  
      write(vname,'("t",I1.1,"_crsoil")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        call cable_pack(dat,casaflux(tile)%crsoil,tile,nmp(:,n))
      end do  
      write(vname,'("t",I1.1,"_cnpp")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        call cable_pack(dat,casaflux(tile)%cnpp,tile,nmp(:,n))
      end do
      write(vname,'("t",I1.1,"_clabloss")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        call cable_pack(dat,casaflux(tile)%clabloss,tile,nmp(:,n))
      end do
      write(vname,'("t",I1.1,"_crgplant")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        call cable_pack(dat,casaflux(tile)%crgplant,tile,nmp(:,n))
      end do
      write(vname,'("t",I1.1,"_stemnpp")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        call cable_pack(dat,casaflux(tile)%stemnpp,tile,nmp(:,n))
      end do
      write(vname,'("t",I1.1,"_LAImax")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        call cable_pack(dat,casabal(tile)%laimax,tile,nmp(:,n))
      end do
      write(vname,'("t",I1.1,"_Cleafmean")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        call cable_pack(dat,casabal(tile)%cleafmean,tile,nmp(:,n))
      end do
      write(vname,'("t",I1.1,"_Crootmean")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        call cable_pack(dat,casabal(tile)%crootmean,tile,nmp(:,n))
      end do
      write(vname,'("t",I1.1,"_fpn")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        call cable_pack(dat,canopy(tile)%fpn,tile,nmp(:,n))
      end do
      write(vname,'("t",I1.1,"_frday")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        call cable_pack(dat,canopy(tile)%frday,tile,nmp(:,n))
      end do
    end do
  end if ! ierr_casa/=0 ..else..
end if   ! ccycle/=0

if ( cable_pop==1 ) then
  if ( ierr_pop/=0 ) then
    if ( myid==0 ) write(6,*) "-> Use default data to initialise POP"
  else
    if ( myid==0 ) write(6,*) "-> Use tiled data to initialise POP"    
    allocate( datpatch(ifull,POP_NPATCH) )  
    allocate( datage(ifull,POP_AGEMAX) )  
    allocate( datpc(ifull,POP_NPATCH,POP_NCOHORT) )
    datpatch = 0._8
    datage = 0._8
    datpc = 0._8
    allocate( dat_out(mp_pop_max), dati_out(mp_pop_max) )
    dat_out = 0._8
    dati_out = 0
    do n = 1,maxtile  
      write(vname,'("t",I1.1,"_pop_grid_cmass_sum")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        np_pop = tdata(tile)%np  
        call pop_pack(dat,dat_out(1:np_pop),tile,nmp(:,n))
        pop(tile)%pop_grid(:)%cmass_sum = dat_out(1:np_pop)
      end do  
      write(vname,'("t",I1.1,"_pop_grid_cmass_sum_old")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        np_pop = tdata(tile)%np  
        call pop_pack(dat,dat_out(1:np_pop),tile,nmp(:,n))
        pop(tile)%pop_grid(:)%cmass_sum_old = dat_out(1:np_pop)
      end do  
      write(vname,'("t",I1.1,"_pop_grid_cheartwood_sum")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        np_pop = tdata(tile)%np  
        call pop_pack(dat,dat_out(1:np_pop),tile,nmp(:,n))
        pop(tile)%pop_grid(:)%cheartwood_sum = dat_out(1:np_pop)
      end do  
      write(vname,'("t",I1.1,"_pop_grid_csapwood_sum")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        np_pop = tdata(tile)%np  
        call pop_pack(dat,dat_out(1:np_pop),tile,nmp(:,n))
        pop(tile)%pop_grid(:)%csapwood_sum = dat_out(1:np_pop)
      end do  
      write(vname,'("t",I1.1,"_pop_grid_csapwood_sum_old")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        np_pop = tdata(tile)%np  
        call pop_pack(dat,dat_out(1:np_pop),tile,nmp(:,n))
        pop(tile)%pop_grid(:)%csapwood_sum_old = dat_out(1:np_pop)
      end do  
      write(vname,'("t",I1.1,"_pop_grid_densindiv")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        np_pop = tdata(tile)%np  
        call pop_pack(dat,dat_out(1:np_pop),tile,nmp(:,n))
        pop(tile)%pop_grid(:)%densindiv = dat_out(1:np_pop)
      end do  
      write(vname,'("t",I1.1,"_pop_grid_height_mean")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        np_pop = tdata(tile)%np  
        call pop_pack(dat,dat_out(1:np_pop),tile,nmp(:,n))
        pop(tile)%pop_grid(:)%height_mean = dat_out(1:np_pop)
      end do  
      write(vname,'("t",I1.1,"_pop_grid_height_max")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        np_pop = tdata(tile)%np  
        call pop_pack(dat,dat_out(1:np_pop),tile,nmp(:,n))
        pop(tile)%pop_grid(:)%height_max = dat_out(1:np_pop) 
      end do  
      write(vname,'("t",I1.1,"_pop_grid_basal_area")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        np_pop = tdata(tile)%np 
        call pop_pack(dat,dat_out(1:np_pop),tile,nmp(:,n))
        pop(tile)%pop_grid(:)%basal_area = dat_out(1:np_pop)
      end do  
      write(vname,'("t",I1.1,"_pop_grid_sapwood_loss")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        np_pop = tdata(tile)%np 
        call pop_pack(dat,dat_out(1:np_pop),tile,nmp(:,n))
        pop(tile)%pop_grid(:)%sapwood_loss = dat_out(1:np_pop)
      end do  
      write(vname,'("t",I1.1,"_pop_grid_sapwood_area_loss")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        np_pop = tdata(tile)%np 
        call pop_pack(dat,dat_out(1:np_pop),tile,nmp(:,n))
        pop(tile)%pop_grid(:)%sapwood_area_loss = dat_out(1:np_pop)
      end do  
      write(vname,'("t",I1.1,"_pop_grid_stress_mortality")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        np_pop = tdata(tile)%np 
        call pop_pack(dat,dat_out(1:np_pop),tile,nmp(:,n))
        pop(tile)%pop_grid(:)%stress_mortality = dat_out(1:np_pop)
      end do  
      write(vname,'("t",I1.1,"_pop_grid_crowding_mortality")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        np_pop = tdata(tile)%np 
        call pop_pack(dat,dat_out(1:np_pop),tile,nmp(:,n))
        pop(tile)%pop_grid(:)%crowding_mortality = dat_out(1:np_pop)
      end do  
      write(vname,'("t",I1.1,"_pop_grid_fire_mortality")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        np_pop = tdata(tile)%np 
        call pop_pack(dat,dat_out(1:np_pop),tile,nmp(:,n))
        pop(tile)%pop_grid(:)%fire_mortality = dat_out(1:np_pop)
      end do  
      write(vname,'("t",I1.1,"_pop_grid_cat_mortality")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        np_pop = tdata(tile)%np 
        call pop_pack(dat,dat_out(1:np_pop),tile,nmp(:,n))
        pop(tile)%pop_grid(:)%cat_mortality = dat_out(1:np_pop)
      end do  
      write(vname,'("t",I1.1,"_pop_grid_res_mortality")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        np_pop = tdata(tile)%np 
        call pop_pack(dat,dat_out(1:np_pop),tile,nmp(:,n))
        pop(tile)%pop_grid(:)%res_mortality = dat_out(1:np_pop)
      end do  
      write(vname,'("t",I1.1,"_pop_grid_growth")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        np_pop = tdata(tile)%np 
        call pop_pack(dat,dat_out(1:np_pop),tile,nmp(:,n))
        pop(tile)%pop_grid(:)%growth = dat_out(1:np_pop)
      end do  
      write(vname,'("t",I1.1,"_pop_grid_area_growth")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        np_pop = tdata(tile)%np 
        call pop_pack(dat,dat_out(1:np_pop),tile,nmp(:,n))
        pop(tile)%pop_grid(:)%area_growth = dat_out(1:np_pop)
      end do  
      write(vname,'("t",I1.1,"_pop_grid_crown_cover")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        np_pop = tdata(tile)%np 
        call pop_pack(dat,dat_out(1:np_pop),tile,nmp(:,n))
        pop(tile)%pop_grid(:)%crown_cover = dat_out(1:np_pop)
      end do  
      write(vname,'("t",I1.1,"_pop_grid_crown_area")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        np_pop = tdata(tile)%np 
        call pop_pack(dat,dat_out(1:np_pop),tile,nmp(:,n))
        pop(tile)%pop_grid(:)%crown_area = dat_out(1:np_pop)
      end do  
      write(vname,'("t",I1.1,"_pop_grid_crown_volume")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        np_pop = tdata(tile)%np 
        call pop_pack(dat,dat_out(1:np_pop),tile,nmp(:,n))
        pop(tile)%pop_grid(:)%crown_volume = dat_out(1:np_pop)
      end do  
      write(vname,'("t",I1.1,"_pop_grid_sapwood_area")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        np_pop = tdata(tile)%np 
        call pop_pack(dat,dat_out(1:np_pop),tile,nmp(:,n))
        pop(tile)%pop_grid(:)%sapwood_area = dat_out(1:np_pop)
      end do  
      write(vname,'("t",I1.1,"_pop_grid_sapwood_area_old")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        np_pop = tdata(tile)%np 
        call pop_pack(dat,dat_out(1:np_pop),tile,nmp(:,n))
        pop(tile)%pop_grid(:)%sapwood_area_old = dat_out(1:np_pop)
      end do  
      write(vname,'("t",I1.1,"_pop_grid_KClump")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      do tile = 1,ntiles
        np_pop = tdata(tile)%np 
        call pop_pack(dat,dat_out(1:np_pop),tile,nmp(:,n))
        pop(tile)%pop_grid(:)%KClump = dat_out(1:np_pop)
      end do  
      write(vname,'("t",I1.1,"_pop_grid_freq_age")') n
      call histrd(iarchi-1,ierr,vname,datage,ifull)
      do k = 1,POP_AGEMAX
        do tile = 1,ntiles
          np_pop = tdata(tile)%np 
          call pop_pack(datage(:,k),dat_out(1:np_pop),tile,nmp(:,n))
          pop(tile)%pop_grid(:)%freq_age(k) = dat_out(1:np_pop)
        end do  
      end do
      write(vname,'("t",I1.1,"_pop_grid_biomass_age")') n
      call histrd(iarchi-1,ierr,vname,datage,ifull)
      do k = 1,POP_AGEMAX
        do tile = 1,ntiles
          np_pop = tdata(tile)%np 
          call pop_pack(datage(:,k),dat_out(1:np_pop),tile,nmp(:,n))
          pop(tile)%pop_grid(:)%biomass_age(k) = dat_out(1:np_pop)
        end do
      end do
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_biomass",I1.1)') n,ll  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        do tile = 1,ntiles
          np_pop = tdata(tile)%np 
          call pop_pack(dat,dat_out(1:np_pop),tile,nmp(:,n))
          pop(tile)%pop_grid(:)%biomass(ll) = dat_out(1:np_pop)
        end do  
      end do
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_density",I1.1)') n,ll
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        do tile = 1,ntiles
          np_pop = tdata(tile)%np 
          call pop_pack(dat,dat_out(1:np_pop),tile,nmp(:,n))
          pop(tile)%pop_grid(:)%density(ll) = dat_out(1:np_pop)
        end do  
      end do  
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_hmean",I1.1)') n,ll  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        do tile = 1,ntiles
          np_pop = tdata(tile)%np 
          call pop_pack(dat,dat_out(1:np_pop),tile,nmp(:,n))
          pop(tile)%pop_grid(:)%hmean(ll) = dat_out(1:np_pop)
        end do  
      end do  
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_hmax",I1.1)') n,ll
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        do tile = 1,ntiles
          np_pop = tdata(tile)%np 
          call pop_pack(dat,dat_out(1:np_pop),tile,nmp(:,n))
          pop(tile)%pop_grid(:)%hmax(ll) = dat_out(1:np_pop)
        end do  
      end do
      do hh = 1,POP_HEIGHT_BINS
        write(vname,'("t",I1.1,"_pop_grid_cmass_stem_bin",I2.2)') n,hh  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        do tile = 1,ntiles
          np_pop = tdata(tile)%np 
          call pop_pack(dat,dat_out(1:np_pop),tile,nmp(:,n)) 
          pop(tile)%pop_grid(:)%cmass_stem_bin(hh) = dat_out(1:np_pop)
        end do  
      end do  
      do hh = 1,POP_HEIGHT_BINS
        write(vname,'("t",I1.1,"_pop_grid_densindiv_bin",I2.2)') n,hh
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        do tile = 1,ntiles
          np_pop = tdata(tile)%np 
          call pop_pack(dat,dat_out(1:np_pop),tile,nmp(:,n))
          pop(tile)%pop_grid(:)%densindiv_bin(hh) = dat_out(1:np_pop)
        end do  
      end do
      do hh = 1,POP_HEIGHT_BINS
        write(vname,'("t",I1.1,"_pop_grid_height_bin",I2.2)') n,hh  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        do tile = 1,ntiles
          np_pop = tdata(tile)%np 
          call pop_pack(dat,dat_out(1:np_pop),tile,nmp(:,n))
          pop(tile)%pop_grid(:)%height_bin(hh) = dat_out(1:np_pop)
        end do  
      end do
      do hh = 1,POP_HEIGHT_BINS
        write(vname,'("t",I1.1,"_pop_grid_diameter_bin",I2.2)') n,hh
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        do tile = 1,ntiles
          np_pop = tdata(tile)%np 
          call pop_pack(dat,dat_out(1:np_pop),tile,nmp(:,n))
          pop(tile)%pop_grid(:)%diameter_bin(hh) = dat_out(1:np_pop)
        end do  
      end do
      do dd = 1,POP_NDISTURB
        write(vname,'("t",I1.1,"_pop_grid_n_age",I1.1)') n,dd  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        dati = nint(dat)
        do tile = 1,ntiles
          np_pop = tdata(tile)%np 
          call pop_pack(dati,dati_out(1:np_pop),tile,nmp(:,n))
          pop(tile)%pop_grid(:)%n_age(dd) = dati_out(1:np_pop)
        end do  
      end do  
      write(vname,'("t",I1.1,"_pop_grid_patch_id")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        dati = nint(datpatch(:,k))  
        do tile = 1,ntiles
          np_pop = tdata(tile)%np        
          call pop_pack(dati,dati_out(1:np_pop),tile,nmp(:,n))
          pop(tile)%pop_grid(:)%patch(k)%id = dati_out(1:np_pop)
        end do  
      end do
      write(vname,'("t",I1.1,"_pop_grid_freq")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles
          np_pop = tdata(tile)%np
          call pop_pack(datpatch(:,k),dat_out(1:np_pop),tile,nmp(:,n))
          pop(tile)%pop_grid(:)%freq(k) = dat_out(1:np_pop)
        end do  
      end do
      write(vname,'("t",I1.1,"_pop_grid_freq_old")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles
          np_pop = tdata(tile)%np
          call pop_pack(datpatch(:,k),dat_out(1:np_pop),tile,nmp(:,n))
          pop(tile)%pop_grid(:)%freq_old(k) = dat_out(1:np_pop)
        end do  
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_factor_recruit")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles
          np_pop = tdata(tile)%np
          call pop_pack(datpatch(:,k),dat_out(1:np_pop),tile,nmp(:,n))
          pop(tile)%pop_grid(:)%patch(k)%factor_recruit = dat_out(1:np_pop)
        end do  
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_pgap")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles
          np_pop = tdata(tile)%np
          call pop_pack(datpatch(:,k),dat_out(1:np_pop),tile,nmp(:,n))
          pop(tile)%pop_grid(:)%patch(k)%pgap = dat_out(1:np_pop)
        end do  
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_lai")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles
          np_pop = tdata(tile)%np
          call pop_pack(datpatch(:,k),dat_out(1:np_pop),tile,nmp(:,n))
          pop(tile)%pop_grid(:)%patch(k)%lai = dat_out(1:np_pop)
        end do
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_biomass")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles
          np_pop = tdata(tile)%np          
          call pop_pack(datpatch(:,k),dat_out(1:np_pop),tile,nmp(:,n))
          pop(tile)%pop_grid(:)%patch(k)%biomass = dat_out(1:np_pop)
        end do  
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_biomass_old")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles
          np_pop = tdata(tile)%np          
          call pop_pack(datpatch(:,k),dat_out(1:np_pop),tile,nmp(:,n))
          pop(tile)%pop_grid(:)%patch(k)%biomass_old = dat_out(1:np_pop)
        end do  
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles
          np_pop = tdata(tile)%np
          call pop_pack(datpatch(:,k),dat_out(1:np_pop),tile,nmp(:,n))
          pop(tile)%pop_grid(:)%patch(k)%sapwood = dat_out(1:np_pop)
        end do  
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_heartwood")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles
          np_pop = tdata(tile)%np          
          call pop_pack(datpatch(:,k),dat_out(1:np_pop),tile,nmp(:,n))
          pop(tile)%pop_grid(:)%patch(k)%heartwood = dat_out(1:np_pop)
        end do  
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_old")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles
          np_pop = tdata(tile)%np
          call pop_pack(datpatch(:,k),dat_out(1:np_pop),tile,nmp(:,n))
          pop(tile)%pop_grid(:)%patch(k)%sapwood_old = dat_out(1:np_pop)
        end do  
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_area")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles
          np_pop = tdata(tile)%np
          call pop_pack(datpatch(:,k),dat_out(1:np_pop),tile,nmp(:,n))
          pop(tile)%pop_grid(:)%patch(k)%sapwood_area = dat_out(1:np_pop)
        end do  
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_area_old")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles
          np_pop = tdata(tile)%np
          call pop_pack(datpatch(:,k),dat_out(1:np_pop),tile,nmp(:,n))
          pop(tile)%pop_grid(:)%patch(k)%sapwood_area_old = dat_out(1:np_pop)
        end do
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_stress_mortality")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles
          np_pop = tdata(tile)%np
          call pop_pack(datpatch(:,k),dat_out(1:np_pop),tile,nmp(:,n))
          pop(tile)%pop_grid(:)%patch(k)%stress_mortality = dat_out(1:np_pop)
        end do  
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_fire_mortality")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles
          np_pop = tdata(tile)%np
          call pop_pack(datpatch(:,k),dat_out(1:np_pop),tile,nmp(:,n))
          pop(tile)%pop_grid(:)%patch(k)%fire_mortality = dat_out(1:np_pop)
        end do  
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_cat_mortality")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles
          np_pop = tdata(tile)%np
          call pop_pack(datpatch(:,k),dat_out(1:np_pop),tile,nmp(:,n))
          pop(tile)%pop_grid(:)%patch(k)%cat_mortality = dat_out(1:np_pop)
        end do  
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_crowding_mortality")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles
          np_pop = tdata(tile)%np
          call pop_pack(datpatch(:,k),dat_out(1:np_pop),tile,nmp(:,n))
          pop(tile)%pop_grid(:)%patch(k)%crowding_mortality = dat_out(1:np_pop)
        end do  
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_cpc")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles
          np_pop = tdata(tile)%np
          call pop_pack(datpatch(:,k),dat_out(1:np_pop),tile,nmp(:,n))
          pop(tile)%pop_grid(:)%patch(k)%cpc = dat_out(1:np_pop)
        end do  
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_mortality")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles
          np_pop = tdata(tile)%np          
          call pop_pack(datpatch(:,k),dat_out(1:np_pop),tile,nmp(:,n))
          pop(tile)%pop_grid(:)%patch(k)%mortality = dat_out(1:np_pop)
        end do  
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_loss")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles
          np_pop = tdata(tile)%np
          call pop_pack(datpatch(:,k),dat_out(1:np_pop),tile,nmp(:,n))
          pop(tile)%pop_grid(:)%patch(k)%sapwood_loss = dat_out(1:np_pop)
        end do  
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_area_loss")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles
          np_pop = tdata(tile)%np          
          call pop_pack(datpatch(:,k),dat_out(1:np_pop),tile,nmp(:,n))
          pop(tile)%pop_grid(:)%patch(k)%sapwood_area_loss = dat_out(1:np_pop)
        end do  
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_growth")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles
          np_pop = tdata(tile)%np
          call pop_pack(datpatch(:,k),dat_out(1:np_pop),tile,nmp(:,n))
          pop(tile)%pop_grid(:)%patch(k)%growth = dat_out(1:np_pop)
        end do  
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_area_growth")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles
          np_pop = tdata(tile)%np          
          call pop_pack(datpatch(:,k),dat_out(1:np_pop),tile,nmp(:,n))
          pop(tile)%pop_grid(:)%patch(k)%area_growth = dat_out(1:np_pop)
        end do  
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_frac_NPP")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles
          np_pop = tdata(tile)%np          
          call pop_pack(datpatch(:,k),dat_out(1:np_pop),tile,nmp(:,n))
          pop(tile)%pop_grid(:)%patch(k)%frac_NPP = dat_out(1:np_pop)
        end do  
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_frac_respiration")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles
          np_pop = tdata(tile)%np
          call pop_pack(datpatch(:,k),dat_out(1:np_pop),tile,nmp(:,n))
          pop(tile)%pop_grid(:)%patch(k)%frac_respiration = dat_out(1:np_pop)
        end do  
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_frac_light_uptake")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles
          np_pop = tdata(tile)%np
          call pop_pack(datpatch(:,k),dat_out(1:np_pop),tile,nmp(:,n))
          pop(tile)%pop_grid(:)%patch(k)%frac_light_uptake = dat_out(1:np_pop)
        end do  
      end do
      do dd = 1,POP_NDISTURB
        write(vname,'("t",I1.1,"_pop_grid_patch_disturbance_interval",I1.1)') n,dd
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        do k = 1,POP_NPATCH
          dati = nint(datpatch(:,k))  
          do tile = 1,ntiles
            np_pop = tdata(tile)%np
            call pop_pack(dati,dati_out(1:np_pop),tile,nmp(:,n))
            pop(tile)%pop_grid(:)%patch(k)%disturbance_interval(dd) = dati_out(1:np_pop)
          end do  
        end do
      end do  
      do dd = 1,POP_NDISTURB
        write(vname,'("t",I1.1,"_pop_grid_patch_first_disturbance_year",I1.1)') n,dd
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        do k = 1,POP_NPATCH
          dati = nint(datpatch(:,k))   
          do tile = 1,ntiles
            np_pop = tdata(tile)%np
            call pop_pack(dati,dati_out(1:np_pop),tile,nmp(:,n))
            pop(tile)%pop_grid(:)%patch(k)%first_disturbance_year(dd) = dati_out(1:np_pop)
          end do  
        end do
      end do  
      do dd = 1,POP_NDISTURB
        write(vname,'("t",I1.1,"_pop_grid_patch_age",I1.1)') n,dd
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        do k = 1,POP_NPATCH
          dati = nint(datpatch(:,k)) 
          do tile = 1,ntiles
            np_pop = tdata(tile)%np
            call pop_pack(dati,dati_out(1:np_pop),tile,nmp(:,n))
            pop(tile)%pop_grid(:)%patch(k)%age(dd) = dati_out(1:np_pop)
          end do  
        end do
      end do  
      do dd = 1,POP_NDISTURB
        write(vname,'("t",I1.1,"_pop_grid_ranked_age_unique",I1.1)') n,dd
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        do k = 1,POP_NPATCH
          dati = nint(datpatch(:,k))  
          do tile = 1,ntiles
            np_pop = tdata(tile)%np
            call pop_pack(dati,dati_out(1:np_pop),tile,nmp(:,n))
            pop(tile)%pop_grid(:)%ranked_age_unique(k,dd) = dati_out(1:np_pop)
          end do  
        end do
      end do  
      do dd = 1,POP_NDISTURB
        write(vname,'("t",I1.1,"_pop_grid_freq_ranked_age_unique",I1.1)') n,dd
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        do k = 1,POP_NPATCH
          do tile = 1,ntiles
            np_pop = tdata(tile)%np
            call pop_pack(datpatch(:,k),dat_out(1:np_pop),tile,nmp(:,n))
            pop(tile)%pop_grid(:)%freq_ranked_age_unique(k,dd) = dat_out(1:np_pop)
          end do  
        end do
      end do  
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_ncohort")') n,ll
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        do k = 1,POP_NPATCH  
          dati = nint(datpatch(:,k))  
          do tile = 1,ntiles
            np_pop = tdata(tile)%np
            call pop_pack(dati,dati_out(1:np_pop),tile,nmp(:,n))
            pop(tile)%pop_grid(:)%patch(k)%layer(ll)%ncohort = dati_out(1:np_pop)
          end do  
        end do
      end do  
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_biomass")') n,ll
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        do k = 1,POP_NPATCH  
          do tile = 1,ntiles
            np_pop = tdata(tile)%np
            call pop_pack(datpatch(:,k),dat_out(1:np_pop),tile,nmp(:,n))
            pop(tile)%pop_grid(:)%patch(k)%layer(ll)%biomass = dat_out(1:np_pop)
          end do  
        end do
      end do  
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_density")') n,ll
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        do k = 1,POP_NPATCH  
          do tile = 1,ntiles
            np_pop = tdata(tile)%np
            call pop_pack(datpatch(:,k),dat_out(1:np_pop),tile,nmp(:,n))
            pop(tile)%pop_grid(:)%patch(k)%layer(ll)%density = dat_out(1:np_pop)
          end do  
        end do
      end do  
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_hmean")') n,ll
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        do k = 1,POP_NPATCH  
          do tile = 1,ntiles
            np_pop = tdata(tile)%np
            call pop_pack(datpatch(:,k),dat_out(1:np_pop),tile,nmp(:,n))
            pop(tile)%pop_grid(:)%patch(k)%layer(ll)%hmean = dat_out(1:np_pop)
          end do  
        end do
      end do  
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_hmax")') n,ll
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        do k = 1,POP_NPATCH  
          do tile = 1,ntiles
            np_pop = tdata(tile)%np
            call pop_pack(datpatch(:,k),dat_out(1:np_pop),tile,nmp(:,n))
            pop(tile)%pop_grid(:)%patch(k)%layer(ll)%hmax = dat_out(1:np_pop)
          end do  
        end do
      end do  
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_age")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT                  
          do k = 1,POP_NPATCH
            dati = nint(datpc(:,k,cc))  
            do tile = 1,ntiles
              np_pop = tdata(tile)%np
              call pop_pack(dati,dati_out(1:np_pop),tile,nmp(:,n))
              pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%age = dati_out(1:np_pop)
            end do  
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_id")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH
            dati = nint(datpc(:,k,cc))  
            do tile = 1,ntiles
              np_pop = tdata(tile)%np
              call pop_pack(dati,dati_out(1:np_pop),tile,nmp(:,n))
              pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%id = dati_out(1:np_pop)
            end do  
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_biomass")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH  
            do tile = 1,ntiles
              np_pop = tdata(tile)%np
              call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),tile,nmp(:,n))
              pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%biomass = dat_out(1:np_pop)
            end do  
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_density")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH  
            do tile = 1,ntiles
              np_pop = tdata(tile)%np
              call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),tile,nmp(:,n))
              pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%density = dat_out(1:np_pop)
            end do  
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_resource_uptake")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH  
            do tile = 1,ntiles
              np_pop = tdata(tile)%np
              call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),tile,nmp(:,n))
              pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_resource_uptake = dat_out(1:np_pop)
            end do  
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_light_uptake")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH  
            do tile = 1,ntiles
              np_pop = tdata(tile)%np
              call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),tile,nmp(:,n))
              pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_light_uptake = dat_out(1:np_pop)
            end do  
          end do  
        end do
      end do        
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_interception")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            do tile = 1,ntiles
              np_pop = tdata(tile)%np
              call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),tile,nmp(:,n))
              pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_interception = dat_out(1:np_pop)
            end do  
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_respiration")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            do tile = 1,ntiles
              np_pop = tdata(tile)%np
              call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),tile,nmp(:,n))
              pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_respiration = dat_out(1:np_pop)
            end do  
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_NPP")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            do tile = 1,ntiles
              np_pop = tdata(tile)%np
              call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),tile,nmp(:,n))
              pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_NPP = dat_out(1:np_pop)
            end do  
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_respiration_scalar")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            do tile = 1,ntiles
              np_pop = tdata(tile)%np              
              call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),tile,nmp(:,n))
              pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%respiration_scalar = dat_out(1:np_pop)
            end do  
          end do  
        end do
      end do         
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_crown_area")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            do tile = 1,ntiles
              np_pop = tdata(tile)%np
              call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),tile,nmp(:,n))
              pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%crown_area = dat_out(1:np_pop)
            end do  
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Pgap")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            do tile = 1,ntiles
              np_pop = tdata(tile)%np
              call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),tile,nmp(:,n))
              pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%Pgap = dat_out(1:np_pop)
            end do  
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_height")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            do tile = 1,ntiles
              np_pop = tdata(tile)%np
              call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),tile,nmp(:,n))
              pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%height = dat_out(1:np_pop)
            end do  
          end do  
        end do
      end do
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_diameter")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            do tile = 1,ntiles
              np_pop = tdata(tile)%np
              call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),tile,nmp(:,n))
              pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%diameter = dat_out(1:np_pop)
            end do  
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_sapwood")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            do tile = 1,ntiles
              np_pop = tdata(tile)%np
              call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),tile,nmp(:,n))
              pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%sapwood = dat_out(1:np_pop)
            end do  
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_heartwood")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            do tile = 1,ntiles
              np_pop = tdata(tile)%np
              call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),tile,nmp(:,n))
              pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%heartwood = dat_out(1:np_pop)
            end do  
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_sapwood_area")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            do tile = 1,ntiles
              np_pop = tdata(tile)%np              
              call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),tile,nmp(:,n))
              pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%sapwood_area = dat_out(1:np_pop)
            end do  
          end do  
        end do
      end do         
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_basal_area")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            do tile = 1,ntiles
              np_pop = tdata(tile)%np
              call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),tile,nmp(:,n))
              pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%basal_area = dat_out(1:np_pop)
            end do  
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_LAI")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            do tile = 1,ntiles
              np_pop = tdata(tile)%np
              call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),tile,nmp(:,n))
              pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%LAI = dat_out(1:np_pop)
            end do  
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Cleaf")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            do tile = 1,ntiles
              np_pop = tdata(tile)%np
              call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),tile,nmp(:,n))
              pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%Cleaf = dat_out(1:np_pop)
            end do  
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Croot")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            do tile = 1,ntiles
              np_pop = tdata(tile)%np
              call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),tile,nmp(:,n))
              pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%Croot = dat_out(1:np_pop)
            end do  
          end do  
        end do
      end do 
    end do
    deallocate( datpatch )
    deallocate( datage )
    deallocate( datpc )
    deallocate( dat_out, dati_out )
  end if ! ierr_pop/=0
end if   ! cable_pop==1  
  
if ( ierr==0 ) then
  ! albvisdir, albvisdif, albnirdir, albnirdif are used when nrad=5
  vname = 'albvisdir'
  call histrd(iarchi-1,ierr,vname,albvisdir,ifull)
  vname = 'albvisdif'
  call histrd(iarchi-1,ierr,vname,albvisdif,ifull)
  vname = 'albnirdir'
  call histrd(iarchi-1,ierr,vname,albnirdir,ifull)
  vname = 'albnirdif'
  call histrd(iarchi-1,ierr,vname,albnirdif,ifull)
  ! albvissav and albnirsav are used when nrad=4
  vname = 'albvis'
  call histrd(iarchi-1,ierr,vname,albvissav,ifull)
  vname = 'albnir'
  call histrd(iarchi-1,ierr,vname,albnirsav,ifull)
end if ! ierr==0  
  
call redistribute_tile
call fixtile

! Calculate LAI and veg fraction diagnostics
vlai(:) = 0.
sigmf(:) = 0.
do tile = 1,ntiles
  if ( tdata(tile)%mp>0 ) then
    call setlai(tdata(tile)%sv,tdata(tile)%vl2,casamet(tile),veg(tile),tdata(tile)%mp)
    call cable_update(vlai,veg(tile)%vlai,tile)
    dummy_unpack(1:tdata(tile)%mp) = 1._8-exp(-0.4_8*veg(tile)%vlai(:))
    call cable_update(sigmf,dummy_unpack(1:tdata(tile)%mp),tile)
  end if
end do  

return
end subroutine loadtile

subroutine loadtile_carbonpools(dat_in,dat_out,tile,mm)

use newmpar_m

implicit none

integer, intent(in) :: mm, tile
integer nb, k, msize, js, je
real, dimension(:,:), intent(in) :: dat_in
real(kind=8), dimension(:,:), intent(inout) :: dat_out
real(kind=8), dimension(size(dat_out,1)) :: dummy_unpack

msize = size(dat_in,2)
js = (tile-1)*imax + 1
je = tile*imax

do k = 1,msize
  dummy_unpack(:) = 0._8  
  call cable_pack(dat_in(js:je,k),dummy_unpack(:),tile)
  where ( veg(tile)%iveg(:)==mm .and. dummy_unpack(:)>1.e-8_8 )
    dat_out(:,k) = dummy_unpack(:)
  end where
end do

return
end subroutine loadtile_carbonpools

subroutine fixtile

!use carbpools_m
use cc_mpi
use darcdf_m
use infile
use newmpar_m
use parm_m
use soil_m
use soilsnow_m
use vegpar_m
  
integer k, tile
real totdepth
 
! Some fixes for rounding errors
do tile = 1,ntiles

  if ( tdata(tile)%mp>0 ) then

    totdepth = 0.
    do k = 1,cbm_ms
      totdepth = totdepth + real(soil(tile)%zse(k))*100.
    enddo
  
    ssnow(tile)%tgg = max(ssnow(tile)%tgg, 200._8)
    ssnow(tile)%tggsn = max(ssnow(tile)%tggsn, 200._8)

    ssnow(tile)%wb = max(ssnow(tile)%wb,0._8)
    ssnow(tile)%wbice = min( max(ssnow(tile)%wbice, 0._8), ssnow(tile)%wb )
    ssnow(tile)%smass = max(ssnow(tile)%smass,0._8)
    ssnow(tile)%rtsoil = max(ssnow(tile)%rtsoil,0._8)
    ssnow(tile)%snowd = max(ssnow(tile)%snowd,0._8)
    ssnow(tile)%osnowd = max(ssnow(tile)%osnowd,0._8)
    ssnow(tile)%wetfac = min(max(ssnow(tile)%wetfac,0._8),1._8)
    canopy(tile)%cansto = max(canopy(tile)%cansto,0._8)
    ssnow(tile)%wbliq = ssnow(tile)%wb - ssnow(tile)%wbice

    ssnow(tile)%wbtot = 0._8
    ssnow(tile)%wbtot1 = 0._8
    ssnow(tile)%wbtot2 = 0._8
    ssnow(tile)%tggav = 0._8
    do k = 1,cbm_ms
      ssnow(tile)%wbtot = ssnow(tile)%wbtot+ssnow(tile)%wb(:,k)*1000._8*soil(tile)%zse(k)
      ssnow(tile)%tggav = ssnow(tile)%tggav+soil(tile)%zse(k)*ssnow(tile)%tgg(:,k)/real(totdepth/100.,8)
      ssnow(tile)%gammzz(:,k) = max((1._8-soil(tile)%ssat)*soil(tile)%css* soil(tile)%rhosoil                         &
          + real(ssnow(tile)%wb(:,k)-ssnow(tile)%wbice(:,k))*4.218e3_8* 1000._8                                       &
          + real(ssnow(tile)%wbice(:,k))*2.100e3_8*1000._8*0.9_8,soil(tile)%css*soil(tile)%rhosoil)*soil(tile)%zse(k) &
          + (1.-ssnow(tile)%isflag)*2090._8*ssnow(tile)%snowd
    end do

    if ( ccycle > 0 ) then
      casapool(tile)%cplant     = max(0._8,casapool(tile)%cplant)
      casapool(tile)%clitter    = max(0._8,casapool(tile)%clitter)
      casapool(tile)%csoil      = max(0._8,casapool(tile)%csoil)
      casabal(tile)%cplantlast  = casapool(tile)%cplant
      casabal(tile)%clitterlast = casapool(tile)%clitter
      casabal(tile)%csoillast   = casapool(tile)%csoil
      casabal(tile)%clabilelast = casapool(tile)%clabile
      casabal(tile)%sumcbal      = 0._8
      casabal(tile)%FCgppyear    = 0._8
      casabal(tile)%FCrpyear     = 0
      casabal(tile)%FCrmleafyear = 0._8
      casabal(tile)%FCrmwoodyear = 0._8
      casabal(tile)%FCrmrootyear = 0._8
      casabal(tile)%FCrgrowyear  = 0._8
      casabal(tile)%FCnppyear    = 0._8
      casabal(tile)%FCrsyear     = 0._8
      casabal(tile)%FCneeyear    = 0._8
      casabal(tile)%dCdtyear     = 0._8
      casapool(tile)%nplant      = max(1.e-6_8,casapool(tile)%nplant)
      casapool(tile)%nlitter     = max(1.e-6_8,casapool(tile)%nlitter)
      casapool(tile)%nsoil       = max(1.e-6_8,casapool(tile)%nsoil)
      casapool(tile)%nsoilmin    = max(1.e-6_8,casapool(tile)%nsoilmin)
      casabal(tile)%nplantlast   = casapool(tile)%nplant
      casabal(tile)%nlitterlast  = casapool(tile)%nlitter
      casabal(tile)%nsoillast    = casapool(tile)%nsoil
      casabal(tile)%nsoilminlast = casapool(tile)%nsoilmin
      casabal(tile)%sumnbal      = 0._8
      casabal(tile)%FNdepyear    = 0._8
      casabal(tile)%FNfixyear    = 0._8
      casabal(tile)%FNsnetyear   = 0._8
      casabal(tile)%FNupyear     = 0._8
      casabal(tile)%FNleachyear  = 0._8
      casabal(tile)%FNlossyear   = 0._8
      casapool(tile)%pplant      = max(1.0e-7_8,casapool(tile)%pplant)
      casapool(tile)%plitter     = max(1.0e-7_8,casapool(tile)%plitter)
      casapool(tile)%psoil       = max(1.0e-7_8,casapool(tile)%psoil)
      casapool(tile)%Psoillab    = max(1.0e-7_8,casapool(tile)%psoillab)
      casapool(tile)%psoilsorb   = max(1.0e-7_8,casapool(tile)%psoilsorb)
      casapool(tile)%psoilocc    = max(1.0e-7_8,casapool(tile)%psoilocc)
      casabal(tile)%pplantlast   = casapool(tile)%pplant
      casabal(tile)%plitterlast  = casapool(tile)%plitter
      casabal(tile)%psoillast    = casapool(tile)%psoil
      casabal(tile)%psoillablast = casapool(tile)%psoillab
      casabal(tile)%psoilsorblast= casapool(tile)%psoilsorb
      casabal(tile)%psoilocclast = casapool(tile)%psoilocc
      casabal(tile)%sumpbal      = 0._8
      casabal(tile)%FPweayear    = 0._8
      casabal(tile)%FPdustyear   = 0._8
      casabal(tile)%FPsnetyear   = 0._8
      casabal(tile)%FPupyear     = 0._8
      casabal(tile)%FPleachyear  = 0._8
      casabal(tile)%FPlossyear   = 0._8
      !casamet(tile)%glai         = max(min( casamet(tile)%glai, &
      !                casabiome(tile)%glaimax(veg(tile)%iveg)), & 
      !                casabiome(tile)%glaimin(veg(tile)%iveg))
    
      where ( .not.( casamet(tile)%iveg2==forest.or.casamet(tile)%iveg2==shrub ) )
        casapool(tile)%cplant(:,wood)  = 0._8
        casapool(tile)%clitter(:,cwd)  = 0._8
        casapool(tile)%nplant(:,wood)  = 0._8
        casapool(tile)%nlitter(:,cwd)  = 0._8
        casapool(tile)%pplant(:,wood)  = 0._8
        casapool(tile)%plitter(:,cwd)  = 0._8
      end where

      ! initializing glai in case not reading pool file (eg. during spin)
      casapool(tile)%ratioNClitter = casapool(tile)%nlitter/(casapool(tile)%clitter+1.0e-10_8)
      casapool(tile)%ratioNPlitter = casapool(tile)%nlitter/(casapool(tile)%plitter+1.0e-10_8)
      casapool(tile)%ratioPClitter = casapool(tile)%plitter/(casapool(tile)%clitter+1.0e-10_8)

      if ( ccycle<2 ) then
        casapool(tile)%Nplant = casapool(tile)%Cplant*casapool(tile)%ratioNCplant
        casapool(tile)%Nsoil  = casapool(tile)%ratioNCsoil*casapool(tile)%Csoil
      else if ( ccycle<3 ) then
        casapool(tile)%Psoil  = casapool(tile)%Nsoil/casapool(tile)%ratioNPsoil
      end if
    
      where ( veg(tile)%iveg==6 .or. veg(tile)%iveg==7 .or. veg(tile)%iveg==8 .or. &
              veg(tile)%iveg==9 .or. veg(tile)%iveg==10 )
        casapool(tile)%cplant(:,wood) = 0._8
        casapool(tile)%clitter(:,cwd) = 0._8
        casapool(tile)%nplant(:,wood) = 0._8
        casapool(tile)%nlitter(:,cwd) = 0._8  
        casapool(tile)%pplant(:,wood) = 0._8
        casapool(tile)%plitter(:,cwd) = 0._8   
      end where
      where ( veg(tile)%iveg==11 .or. veg(tile)%iveg==12 .or. veg(tile)%iveg==13 .or. &
              veg(tile)%iveg==15 .or. veg(tile)%iveg==16 .or. veg(tile)%iveg==17 )
        casapool(tile)%cplant(:,leaf) = 0._8
        casapool(tile)%cplant(:,wood) = 0._8
        casapool(tile)%cplant(:,froot) = 0._8
        casapool(tile)%clitter(:,metb) = 0._8
        casapool(tile)%clitter(:,str) = 0._8
        casapool(tile)%clitter(:,cwd) = 0._8
        casapool(tile)%csoil(:,mic) = 0._8
        casapool(tile)%csoil(:,slow) = 0._8
        casapool(tile)%csoil(:,pass) = 0._8
        casapool(tile)%nplant(:,leaf) = 0._8
        casapool(tile)%nplant(:,wood) = 0._8
        casapool(tile)%nplant(:,froot) = 0._8
        casapool(tile)%nlitter(:,metb) = 0._8
        casapool(tile)%nlitter(:,str) = 0._8
        casapool(tile)%nlitter(:,cwd) = 0._8
        casapool(tile)%nsoil(:,mic) = 0._8
        casapool(tile)%nsoil(:,slow) = 0._8
        casapool(tile)%nsoil(:,pass) = 0._8
        casapool(tile)%pplant(:,leaf) = 0._8
        casapool(tile)%pplant(:,wood) = 0._8
        casapool(tile)%pplant(:,froot) = 0._8
        casapool(tile)%plitter(:,metb) = 0._8
        casapool(tile)%plitter(:,str) = 0._8
        casapool(tile)%plitter(:,cwd) = 0._8
        casapool(tile)%psoil(:,mic) = 0._8
        casapool(tile)%psoil(:,slow) = 0._8
        casapool(tile)%psoil(:,pass) = 0._8
      end where
      where ( veg(tile)%iveg==14 )        
        casapool(tile)%cplant(:,leaf) = 0._8
        casapool(tile)%cplant(:,wood) = 0._8
        casapool(tile)%cplant(:,froot) = 0._8
        casapool(tile)%pplant(:,leaf) = 0._8
        casapool(tile)%pplant(:,wood) = 0._8
        casapool(tile)%pplant(:,froot) = 0._8
      end where
    
    end if ! ccycle>0
  end if   ! mp>0
  
end do ! tile
  
return
end subroutine fixtile

! if the vegetation class has changed for a tile, then attempt to find the tile associated
! with the equivilent vegetation class in the restart data
subroutine create_new_tile_map(nmp)

use cc_mpi
use newmpar_m, only : ifull, imax, ntiles

integer n, m, iq, store_n
integer tile, js, je, i
integer, dimension(ifull,maxtile), intent(inout) :: nmp
integer, dimension(imax,maxtile) :: oldv_up, newv_up
integer, dimension(maxtile) :: oldv_v, newv_v
logical cv_test

if ( myid==0 ) write(6,*) "-> Create map to vegetation types"

do tile = 1,ntiles
  js = (tile-1)*imax + 1
  je = tile*imax
  
  if ( tdata(tile)%mp>0 ) then
    cv_test = any( tdata(tile)%cveg/=tdata(tile)%old_cveg )
  else
    cv_test = .false.
  end if

  oldv_up = 0
  newv_up = 0

  ! default map  
  do iq = js,je
    do n = 1,maxtile  
      nmp(iq,n) = n  
    end do  
  end do

  if ( tdata(tile)%mp>0 ) then
    if ( cv_test ) then
    
      do n = 1,maxtile
        call cable_unpack(tdata(tile)%cveg,newv_up(:,n),tile,nb=n)
        call cable_unpack(tdata(tile)%old_cveg,oldv_up(:,n),tile,nb=n)
      end do
    
      do i = 1,imax
        iq = i + js - 1
        newv_v = newv_up(i,:) ! current vegetation class
        oldv_v = oldv_up(i,:) ! vegetation class in restart file
        ! check each tile for mismatch
        do n = 1,maxtile
          if ( oldv_v(n)/=newv_v(n) ) then
            ! search over all tiles for replacement
            do m = 1,maxtile
              if ( oldv_v(m)==newv_v(n) ) then
                ! shuffle tiles to keep area fraction sum equal to 1.  
                store_n   = nmp(iq,m)
                nmp(iq,m) = nmp(iq,n) ! found required tile
                nmp(iq,n) = store_n
                store_n   = oldv_v(m)
                oldv_v(m) = oldv_v(n)
                oldv_v(n) = store_n !=newv_v(n)   
                exit
              end if  
            end do
          end if  
        end do
      end do
    
    end if ! cv_test
  end if   ! mp>0

end do ! tile
  
return
end subroutine create_new_tile_map


! redistribute temperature and water with a gridbox if tile area fraction changes
subroutine redistribute_tile

use cc_mpi
use newmpar_m, only : imax, ntiles

integer k, tile, js, je
logical sv_test

do tile = 1,ntiles
  js = (tile-1)*imax + 1
  je = tile*imax

  if ( tdata(tile)%mp>0 ) then
    ! check if any point have land-cover change  
    sv_test = any( abs(tdata(tile)%sv(:)-tdata(tile)%old_sv(:))>1.e-8 )  
    
    if ( sv_test ) then

      ! check for errors prior to redistribution
      if ( any(ssnow(tile)%tgg>400.) ) then
        write(6,*) "ERROR: Invalid input temperature for CABLE redistribute_tile"
        write(6,*) "ssnow%tgg ",maxval(ssnow(tile)%tgg)
        call ccmpi_abort(-1)
      end if     
    
      ! assume common soil texture and soil heat capacity
      do k = 1,cbm_ms
        call redistribute_work(ssnow(tile)%tgg(:,k),tile)
        call redistribute_work(ssnow(tile)%wb(:,k),tile)
        call redistribute_work(ssnow(tile)%wbice(:,k),tile)
      end do
      call redistribute_work(ssnow(tile)%GWwb,tile)
      if ( soil_struc==1 ) then
        do k = 1,cbm_ms
          call redistribute_work(ssnow(tile)%tsoil(:,k),tile)
        end do
      end if
    
      ! Do we need a special treatment for snow?
      do k = 1,3
        call redistribute_work(ssnow(tile)%tggsn(:,k),tile)
        call redistribute_work(ssnow(tile)%smass(:,k),tile)
        call redistribute_work(ssnow(tile)%ssdn(:,k),tile)
        call redistribute_work(ssnow(tile)%sdepth(:,k),tile)
        call redistribute_work(ssnow(tile)%sconds(:,k),tile)
      end do
      call redistribute_work(ssnow(tile)%ssdnn(:),tile)
      call redistribute_work(ssnow(tile)%snowd(:),tile)
      call redistribute_work(ssnow(tile)%osnowd(:),tile)
      call redistribute_work(ssnow(tile)%snage(:),tile)
      call redistribute_work(ssnow(tile)%rtsoil(:),tile)
      call redistribute_work(ssnow(tile)%GWwb(:),tile)
      call redistribute_work(ssnow(tile)%wtd(:),tile)
      call redistribute_work(canopy(tile)%cansto(:),tile)
      call redistribute_work(ssnow(tile)%pudsto(:),tile)
      call redistribute_work(ssnow(tile)%wetfac(:),tile)
    
      ! also treatment for carbon
      if ( ccycle/=0 ) then
        do k = 1,3    
          call redistribute_work(casapool(tile)%cplant(:,k),tile)  
          call redistribute_work(casapool(tile)%clitter(:,k),tile)
          call redistribute_work(casapool(tile)%csoil(:,k),tile)
          call redistribute_work(casapool(tile)%nplant(:,k),tile)
          call redistribute_work(casapool(tile)%nlitter(:,k),tile)
          call redistribute_work(casapool(tile)%nsoil(:,k),tile)
          call redistribute_work(casapool(tile)%pplant(:,k),tile)
          call redistribute_work(casapool(tile)%plitter(:,k),tile)
          call redistribute_work(casapool(tile)%psoil(:,k),tile)
        end do    
        !call redistribute_work(casamet(tile)%glai,tile)
        !call redistribute_work(phen(tile)%phen,tile)
        !call redistribute_work(phen(tile)%aphen,tile)
      end if

      ! check for errors after redistribution
      if ( any(ssnow(tile)%tgg>400.) ) then
        write(6,*) "ERROR: Invalid output temperature for CABLE redistribute_tile"
        write(6,*) "ssnow%tgg ",maxval(ssnow(tile)%tgg)
        call ccmpi_abort(-1)
      end if  
    
    end if ! svn_test
  end if   ! mp>0
  
end do ! tile  
  
return
end subroutine redistribute_tile

subroutine redistribute_work(vdata,tile)

use newmpar_m, only : imax, ntiles
use soil_m

integer, intent(in) :: tile
integer nb, iq, is, ie, i
real, dimension(imax,maxtile) :: up_new_svs, up_old_svs
real, dimension(imax) :: svs_sum
real, dimension(maxtile) :: adj_pos_frac, adj_neg_frac, new_svs, old_svs
real(kind=8) ave_neg_vdata, adj_neg_frac_sum
real(kind=8), dimension(:), intent(inout) :: vdata
real(kind=8), dimension(imax,maxtile) :: up_vdata
real(kind=8), dimension(maxtile) :: old_vdata

is = (tile-1)*imax + 1
is = tile*imax

! update data
up_vdata = 0._8
up_new_svs = 0.
up_old_svs = 0.
do nb = 1,tdata(tile)%maxnb
  call cable_unpack(tdata(tile)%sv,up_new_svs(:,nb),tile,nb=nb)
  call cable_unpack(tdata(tile)%old_sv,up_old_svs(:,nb),tile,nb=nb)
  call cable_unpack(vdata,up_vdata(:,nb),tile,nb=nb)
end do  

svs_sum = sum(up_new_svs,dim=2)
do nb = 1,tdata(tile)%maxnb
  where ( land(is:ie) )  
    up_new_svs(:,nb) = up_new_svs(:,nb)/max(svs_sum(:),1.e-10)
  end where  
end do
svs_sum = sum(up_old_svs,dim=2)
do nb = 1,tdata(tile)%maxnb
  where ( land(is:ie) )  
    up_old_svs(:,nb) = up_old_svs(:,nb)/max(svs_sum(:),1.e-10)
  end where  
end do

do i = 1,imax
  iq = i + is - 1  
  do nb = 1,tdata(tile)%maxnb
    new_svs(nb) = up_new_svs(i,nb)
    old_svs(nb) = up_old_svs(i,nb)      
    adj_pos_frac(nb) = max( new_svs(nb)-old_svs(nb), 0. ) ! tiles with increasing fraction
    adj_neg_frac(nb) = max( old_svs(nb)-new_svs(nb), 0. ) ! tiles with decreasing fraction
  end do  
  adj_neg_frac_sum = sum( adj_neg_frac(1:tdata(tile)%maxnb) )   
  ! check adj_neg_frac_sum>0 which can occur if water has changed to land
  if ( land(iq) .and. any( abs(new_svs(1:tdata(tile)%maxnb)-old_svs(1:tdata(tile)%maxnb))>1.e-8 ) .and. &
       adj_neg_frac_sum>1.e-8 ) then 
    do nb = 1,tdata(tile)%maxnb  
      old_vdata(nb) = up_vdata(iq,nb)
    end do  
    ! summarise tiles with decreasing area fraction with an average value
    ave_neg_vdata = sum( adj_neg_frac(1:tdata(tile)%maxnb)*old_vdata(1:tdata(tile)%maxnb) ) &
                   /adj_neg_frac_sum
    ! Only change tiles that are increasing in area fraction
    do nb = 1,tdata(tile)%maxnb
      if ( adj_pos_frac(nb)>1.e-8 ) then  
        up_vdata(i,nb) = (old_vdata(nb)*old_svs(nb) + ave_neg_vdata*adj_pos_frac(nb)) &
                        /(old_svs(nb) + adj_pos_frac(nb))
      end if  
    end do
  end if    
end do  

do nb = 1,tdata(tile)%maxnb
  call cable_pack(up_vdata(:,nb),vdata,tile,nb=nb)
end do  

return
end subroutine redistribute_work    
    
! *************************************************************************************
! This subroutine saves CABLE tile data
subroutine savetiledef(idnc,local,jdim,jsize,cdim,csize,itype)

!use carbpools_m
use cc_mpi, only : myid
use infile
use newmpar_m
use parm_m, only : diaglevel_pop
  
integer, intent(in) :: idnc, jsize
integer k,n
integer ll,dd,hh
integer, dimension(2), intent(in) :: csize
integer, dimension(jsize), intent(in) :: jdim  
integer, dimension(6,4), intent(in) :: cdim
character(len=80) vname
character(len=80) lname
logical, intent(in) :: local
integer, intent(in) :: itype
  
if (myid==0.or.local) then
  if (myid==0) write(6,*) "-> define CABLE tile data"
  if ( itype==-1 ) then !just for restart file
    do n = 1,maxtile  
      write(lname,'("Veg type tile ",I1.1)') n
      write(vname,'("t",I1.1,"_cvc")') n
      call attrib(idnc,jdim,jsize,vname,lname,'none',0.,1.,any_m,point_m,land_m,double_m)
    end do      
    do n = 1,maxtile
      write(lname,'("Veg fraction tile ",I1.1)') n
      write(vname,'("t",I1.1,"_svs")') n
      call attrib(idnc,jdim,jsize,vname,lname,'none',0.,1.,any_m,point_m,land_m,double_m)
      do k = 1,cbm_ms
        write(lname,'("Soil temperature tile ",I1.1," lev ",I1.1)') n,k
        write(vname,'("t",I1.1,"_tgg",I1.1)') n,k
        call attrib(idnc,jdim,jsize,vname,lname,'K',100.,400.,any_m,point_m,land_m,double_m)
      end do
      do k = 1,cbm_ms
        write(lname,'("Soil moisture tile ",I1.1," lev ",I1.1)') n,k
        write(vname,'("t",I1.1,"_wb",I1.1)') n,k 
        call attrib(idnc,jdim,jsize,vname,lname,'m3/m3',0.,2.6,any_m,point_m,land_m,double_m)
      end do
      do k = 1,cbm_ms
        write(lname,'("Soil ice tile ",I1.1," lev ",I1.1)') n,k
        write(vname,'("t",I1.1,"_wbice",I1.1)') n,k 
        call attrib(idnc,jdim,jsize,vname,lname,'m3/m3',0.,2.6,any_m,point_m,land_m,double_m)
      end do
      do k = 1,3
        write(lname,'("Snow temperature tile ",I1.1," lev ",I1.1)') n,k
        write(vname,'("t",I1.1,"_tggsn",I1.1)') n,k
        call attrib(idnc,jdim,jsize,vname,lname,'K',100.,400.,any_m,point_m,land_m,double_m)
      end do
      do k = 1,3
        write(lname,'("Snow mass tile ",I1.1," lev ",I1.1)') n,k
        write(vname,'("t",I1.1,"_smass",I1.1)') n,k 
        call attrib(idnc,jdim,jsize,vname,lname,'K',0.,650.,any_m,point_m,land_m,double_m)
      end do
      do k = 1,3
        write(lname,'("Snow density tile ",I1.1," lev ",I1.1)') n,k
        write(vname,'("t",I1.1,"_ssdn",I1.1)') n,k 
        call attrib(idnc,jdim,jsize,vname,lname,'kg/m3',0.,650.,any_m,point_m,land_m,double_m)
      end do
      do k = 1,3
        write(lname,'("Snow depth tile ",I1.1," lev ",I1.1)') n,k
        write(vname,'("t",I1.1,"_sdepth",I1.1)') n,k 
        call attrib(idnc,jdim,jsize,vname,lname,'mm',0.,6500.,any_m,point_m,land_m,double_m)
      end do
      do k = 1,3
        write(lname,'("Snow sconds tile ",I1.1," lev ",I1.1)') n,k
        write(vname,'("t",I1.1,"_sconds",I1.1)') n,k 
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6.5,any_m,point_m,land_m,double_m)
      end do
      write(lname,'("Snow ssdnn tile ",I1.1)') n
      write(vname,'("t",I1.1,"_ssdnn")') n
      call attrib(idnc,jdim,jsize,vname,lname,'kg/m3',0.,650.,any_m,point_m,land_m,double_m)
      write(lname,'("Snow flag tile ",I1.1)') n
      write(vname,'("t",I1.1,"_sflag")') n
      call attrib(idnc,jdim,jsize,vname,lname,'mm',0.,6.5,any_m,point_m,land_m,double_m)
      write(lname,'("Snow depth tile ",I1.1)') n
      write(vname,'("t",I1.1,"_snd")') n
      call attrib(idnc,jdim,jsize,vname,lname,'mm',0.,6500.,any_m,point_m,land_m,double_m)
      write(lname,'("Old snow depth tile ",I1.1)') n
      write(vname,'("t",I1.1,"_osnd")') n
      call attrib(idnc,jdim,jsize,vname,lname,'mm',0.,6500.,any_m,point_m,land_m,double_m)
      write(lname,'("Snow age tile ",I1.1)') n
      write(vname,'("t",I1.1,"_snage")') n
      call attrib(idnc,jdim,jsize,vname,lname,'none',0.,26.,any_m,point_m,land_m,double_m)
      write(lname,'("Soil turbulent resistance tile ",I1.1)') n
      write(vname,'("t",I1.1,"_rtsoil")') n
      call attrib(idnc,jdim,jsize,vname,lname,'none',0.,1.3e5,any_m,point_m,land_m,double_m)
      write(lname,'("Aquifer moisture content ",I1.1)') n
      write(vname,'("t",I1.1,"_GWwb")') n
      call attrib(idnc,jdim,jsize,vname,lname,'m3/m3',0.,1.3e5,any_m,point_m,land_m,double_m)
      write(lname,'("Water table depth ",I1.1)') n
      write(vname,'("t",I1.1,"_wtd")') n
      call attrib(idnc,jdim,jsize,vname,lname,'m',0.,6.5e4,any_m,point_m,land_m,double_m)
      write(lname,'("cansto tile ",I1.1)') n
      write(vname,'("t",I1.1,"_cansto")') n
      call attrib(idnc,jdim,jsize,vname,lname,'none',0.,13.,any_m,point_m,land_m,double_m)
      write(lname,'("us tile ",I1.1)') n
      write(vname,'("t",I1.1,"_us")') n
      call attrib(idnc,jdim,jsize,vname,lname,'m/s',0.,13.,any_m,point_m,land_m,double_m) 
      write(lname,'("pudsto tile ",I1.1)') n
      write(vname,'("t",I1.1,"_pudsto")') n
      call attrib(idnc,jdim,jsize,vname,lname,'none',0.,13.,any_m,point_m,land_m,double_m)
      write(lname,'("wetfac tile ",I1.1)') n
      write(vname,'("t",I1.1,"_wetfac")') n
      call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6.5,any_m,point_m,land_m,double_m)
      write(lname,'("ga tile ",I1.1)') n
      write(vname,'("t",I1.1,"_ga")') n
      call attrib(idnc,jdim,jsize,vname,lname,'none',-6500.,6500.,any_m,point_m,land_m,double_m)
    end do
    if ( soil_struc==1 ) then
      do n = 1,maxtile  
        write(lname,'("hzero tile ",I1.1)') n
        write(vname,'("t",I1.1,"_hzero")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        do k = 1,cbm_ms
          write(lname,'("S tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_s",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        end do
        do k = 1,cbm_ms
          write(lname,'("tsoil tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_tsoil",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        end do
        do k = 1,cbm_ms
          write(lname,'("thetai tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_thetai",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        end do  
        write(lname,'("snowliq tile ",I1.1," lev ",I1.1)') n,1
        write(vname,'("t",I1.1,"_snowliq",I1.1)') n,1
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        write(lname,'("tsurface tile ",I1.1)') n
        write(vname,'("t",I1.1,"_tsurface")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        write(lname,'("nsnow tile ",I1.1)') n
        write(vname,'("t",I1.1,"_nsnow")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        write(lname,'("fwsoil tile ",I1.1)') n
        write(vname,'("t",I1.1,"_fwsoil")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
      end do
    end if
    if ( ccycle/=0 ) then
      do n = 1,maxtile  
        do k = 1,mplant
          write(lname,'("C leaf pool tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_cplant",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        end do
        do k = 1,mplant
          write(lname,'("N leaf pool tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_nplant",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        end do
        do k = 1,mplant
          write(lname,'("P leaf pool tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_pplant",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        end do
        do k = 1,mlitter
          write(lname,'("C litter pool tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_clitter",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        end do
        do k = 1,mlitter
          write(lname,'("N litter pool tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_nlitter",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        end do
        do k = 1,mlitter
          write(lname,'("P litter pool tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_plitter",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        end do
        do k = 1,msoil
          write(lname,'("C soil pool tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_csoil",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        end do
        do k = 1,msoil
          write(lname,'("N soil pool tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_nsoil",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        end do
        do k = 1,msoil
          write(lname,'("P soil pool tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_psoil",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        end do
        write(lname,'("Prognostic LAI tile ",I1.1)') n
        write(vname,'("t",I1.1,"_glai")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        write(lname,'("Leaf phenology phen tile ",I1.1)') n
        write(vname,'("t",I1.1,"_phen")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        write(lname,'("Leaf phenology rainfall ",I1.1)') n
        write(vname,'("t",I1.1,"_aphen")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        write(lname,'("Leaf phenology phase tile ",I1.1)') n
        write(vname,'("t",I1.1,"_phenphase")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        write(lname,'("Leaf phenology doyphase3 tile ",I1.1)') n
        write(vname,'("t",I1.1,"_doyphase3")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        write(lname,'("C labile tile ",I1.1)') n
        write(vname,'("t",I1.1,"_clabile")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        write(lname,'("N soilmin tile ",I1.1)') n
        write(vname,'("t",I1.1,"_nsoilmin")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        write(lname,'("P soillab tile ",I1.1)') n
        write(vname,'("t",I1.1,"_psoillab")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        write(lname,'("P soilsorb tile ",I1.1)') n
        write(vname,'("t",I1.1,"_psoilsorb")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        write(lname,'("P soilocc tile ",I1.1)') n
        write(vname,'("t",I1.1,"_psoilocc")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        do k = 1,mplant
          write(lname,'("crmplant tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_crmplant",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        end do
        write(lname,'("frac_sapwood tile ",I1.1)') n
        write(vname,'("t",I1.1,"_fracsapwood")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        write(lname,'("sapwoodarea tile ",I1.1)') n
        write(vname,'("t",I1.1,"_sapwoodarea")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        write(lname,'("crsoil tile ",I1.1)') n
        write(vname,'("t",I1.1,"_crsoil")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        write(lname,'("cnpp tile ",I1.1)') n
        write(vname,'("t",I1.1,"_cnpp")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        write(lname,'("clabloss tile ",I1.1)') n
        write(vname,'("t",I1.1,"_clabloss")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        write(lname,'("crgplant tile ",I1.1)') n
        write(vname,'("t",I1.1,"_crgplant")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        write(lname,'("stemnpp tile ",I1.1)') n
        write(vname,'("t",I1.1,"_stemnpp")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        write(lname,'("LAImax tile ",I1.1)') n
        write(vname,'("t",I1.1,"_LAImax")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        write(lname,'("Cleafmean tile ",I1.1)') n
        write(vname,'("t",I1.1,"_Cleafmean")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        write(lname,'("Crootmean tile ",I1.1)') n
        write(vname,'("t",I1.1,"_Crootmean")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        write(lname,'("fpn ",I1.1)') n
        write(vname,'("t",I1.1,"_fpn")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        write(lname,'("frday ",I1.1)') n
        write(vname,'("t",I1.1,"_frday")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
      end do
    end if
  end if ! itype==-1
  if ( cable_pop==1 ) then
    do n = 1,maxtile  
      ! Convention for POP variables are of the form t<n>_pop_grid_<......>
      ! so that ppc2hist can interpolate them
      if ( diaglevel_pop > 0 ) then
        write(lname,'("t",I1.1,"_pop_grid_cmass_sum")') n
        write(vname,'("t",I1.1,"_pop_grid_cmass_sum")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
      end if
      if ( itype==-1 ) then !just for restart file
        write(lname,'("t",I1.1,"_pop_grid_cmass_sum_old")') n
        write(vname,'("t",I1.1,"_pop_grid_cmass_sum_old")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_cheartwood_sum")') n
        write(vname,'("t",I1.1,"_pop_grid_cheartwood_sum")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_csapwood_sum")') n
        write(vname,'("t",I1.1,"_pop_grid_csapwood_sum")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_csapwood_sum_old")') n
        write(vname,'("t",I1.1,"_pop_grid_csapwood_sum_old")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_densindiv")') n
        write(vname,'("t",I1.1,"_pop_grid_densindiv")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_height_mean")') n
        write(vname,'("t",I1.1,"_pop_grid_height_mean")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_height_max")') n
        write(vname,'("t",I1.1,"_pop_grid_height_max")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_basal_area")') n
        write(vname,'("t",I1.1,"_pop_grid_basal_area")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_sapwood_loss")') n
        write(vname,'("t",I1.1,"_pop_grid_sapwood_loss")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_sapwood_area_loss")') n
        write(vname,'("t",I1.1,"_pop_grid_sapwood_area_loss")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_stress_mortality")') n
        write(vname,'("t",I1.1,"_pop_grid_stress_mortality")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_crowding_mortality")') n
        write(vname,'("t",I1.1,"_pop_grid_crowding_mortality")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_fire_mortality")') n
        write(vname,'("t",I1.1,"_pop_grid_fire_mortality")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_cat_mortality")') n
        write(vname,'("t",I1.1,"_pop_grid_cat_mortality")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_res_mortality")') n
        write(vname,'("t",I1.1,"_pop_grid_res_mortality")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_growth")') n
        write(vname,'("t",I1.1,"_pop_grid_growth")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_area_growth")') n
        write(vname,'("t",I1.1,"_pop_grid_area_growth")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_crown_cover")') n
        write(vname,'("t",I1.1,"_pop_grid_crown_cover")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_crown_area")') n
        write(vname,'("t",I1.1,"_pop_grid_crown_area")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_crown_volume")') n
        write(vname,'("t",I1.1,"_pop_grid_crown_volume")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_sapwood_area")') n
        write(vname,'("t",I1.1,"_pop_grid_sapwood_area")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_sapwood_area_old")') n
        write(vname,'("t",I1.1,"_pop_grid_sapwood_area_old")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_KClump")') n
        write(vname,'("t",I1.1,"_pop_grid_KClump")') n
        call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_freq_age")') n
        write(vname,'("t",I1.1,"_pop_grid_freq_age")') n
        call attrib(idnc,cdim(:,3),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_biomass_age")') n
        write(vname,'("t",I1.1,"_pop_grid_biomass_age")') n
        call attrib(idnc,cdim(:,3),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        do ll = 1,POP_NLAYER
          write(lname,'("t",I1.1,"_pop_grid_biomass",I1.1)') n,ll
          write(vname,'("t",I1.1,"_pop_grid_biomass",I1.1)') n,ll
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do      
        do ll = 1,POP_NLAYER
          write(lname,'("t",I1.1,"_pop_grid_density",I1.1)') n,ll
          write(vname,'("t",I1.1,"_pop_grid_density",I1.1)') n,ll
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do
        do ll = 1,POP_NLAYER
          write(lname,'("t",I1.1,"_pop_grid_hmean",I1.1)') n,ll
          write(vname,'("t",I1.1,"_pop_grid_hmean",I1.1)') n,ll
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do         
        do ll = 1,POP_NLAYER
          write(lname,'("t",I1.1,"_pop_grid_hmax",I1.1)') n,ll
          write(vname,'("t",I1.1,"_pop_grid_hmax",I1.1)') n,ll
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do   
        do hh = 1,POP_HEIGHT_BINS
          write(lname,'("t",I1.1,"_pop_grid_cmass_stem_bin",I2.2)') n,hh
          write(vname,'("t",I1.1,"_pop_grid_cmass_stem_bin",I2.2)') n,hh
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do   
        do hh = 1,POP_HEIGHT_BINS
          write(lname,'("t",I1.1,"_pop_grid_densindiv_bin",I2.2)') n,hh
          write(vname,'("t",I1.1,"_pop_grid_densindiv_bin",I2.2)') n,hh
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do         
        do hh = 1,POP_HEIGHT_BINS
          write(lname,'("t",I1.1,"_pop_grid_height_bin",I2.2)') n,hh
          write(vname,'("t",I1.1,"_pop_grid_height_bin",I2.2)') n,hh
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do         
        do hh = 1,POP_HEIGHT_BINS
          write(lname,'("t",I1.1,"_pop_grid_diameter_bin",I2.2)') n,hh
          write(vname,'("t",I1.1,"_pop_grid_diameter_bin",I2.2)') n,hh
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do         
        do dd = 1,POP_NDISTURB
          write(lname,'("t",I1.1,"_pop_grid_n_age",I1.1)') n,dd
          write(vname,'("t",I1.1,"_pop_grid_n_age",I1.1)') n,dd
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do         
      end if
      ! Convention for POP variables are of the form t<n>_pop_grid_<......>
      ! so that ppc2hist can interpolate them
      if ( itype==-1 .or. diaglevel_pop>=9 ) then !just for restart file
        write(lname,'("t",I1.1,"_pop_grid_patch_id")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_id")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_freq")') n
        write(vname,'("t",I1.1,"_pop_grid_freq")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_freq_old")') n
        write(vname,'("t",I1.1,"_pop_grid_freq_old")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_factor_recruit")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_factor_recruit")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_pgap")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_pgap")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_lai")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_lai")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_biomass")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_biomass")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_biomass_old")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_biomass_old")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_sapwood")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_sapwood")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_heartwood")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_heartwood")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_sapwood_old")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_old")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_sapwood_area")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_area")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_sapwood_area_old")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_area_old")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_stress_mortality")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_stress_mortality")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_fire_mortality")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_fire_mortality")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_cat_mortality")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_cat_mortality")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_crowding_mortality")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_crowding_mortality")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_cpc")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_cpc")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_mortality")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_mortality")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_sapwood_loss")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_loss")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_sapwood_area_loss")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_area_loss")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_growth")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_growth")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_area_growth")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_area_growth")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_frac_NPP")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_frac_NPP")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_frac_respiration")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_frac_respiration")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_patch_frac_light_uptake")') n
        write(vname,'("t",I1.1,"_pop_grid_patch_frac_light_uptake")') n
        call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        do dd = 1,POP_NDISTURB  
          write(lname,'("t",I1.1,"_pop_grid_patch_disturbance_interval",I1.1)') n,dd
          write(vname,'("t",I1.1,"_pop_grid_patch_disturbance_interval",I1.1)') n,dd
          call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do  
        do dd = 1,POP_NDISTURB  
          write(lname,'("t",I1.1,"_pop_grid_patch_first_disturbance_year",I1.1)') n,dd
          write(vname,'("t",I1.1,"_pop_grid_patch_first_disturbance_year",I1.1)') n,dd
          call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do  
        do dd = 1,POP_NDISTURB  
          write(lname,'("t",I1.1,"_pop_grid_patch_age",I1.1)') n,dd
          write(vname,'("t",I1.1,"_pop_grid_patch_age",I1.1)') n,dd
          call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do
        do dd = 1,POP_NDISTURB  
          write(lname,'("t",I1.1,"_pop_grid_ranked_age_unique",I1.1)') n,dd
          write(vname,'("t",I1.1,"_pop_grid_ranked_age_unique",I1.1)') n,dd
          call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do
        do dd = 1,POP_NDISTURB  
          write(lname,'("t",I1.1,"_pop_grid_freq_ranked_age_unique",I1.1)') n,dd
          write(vname,'("t",I1.1,"_pop_grid_freq_ranked_age_unique",I1.1)') n,dd
          call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_ncohort")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_ncohort")') n,ll
          call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_biomass")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_biomass")') n,ll
          call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_density")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_density")') n,ll
          call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_hmean")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_hmean")') n,ll
          call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_hmax")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_hmax")') n,ll
          call attrib(idnc,cdim(:,1),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_age")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_age")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_id")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_id")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_biomass")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_biomass")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_density")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_density")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_resource_uptake")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_resource_uptake")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_light_uptake")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_light_uptake")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_interception")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_interception")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_respiration")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_respiration")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_NPP")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_NPP")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_respiration_scalar")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_respiration_scalar")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_crown_area")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_crown_area")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Pgap")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Pgap")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_height")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_height")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_diameter")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_diameter")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_sapwood")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_sapwood")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_heartwood")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_heartwood")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_sapwood_area")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_sapwood_area")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_basal_area")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_basal_area")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_LAI")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_LAI")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Cleaf")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Cleaf")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do
        do ll = 1,POP_NLAYER  
          write(lname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Croot")') n,ll
          write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Croot")') n,ll
          call attrib(idnc,cdim(:,2),csize(2),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        end do
      end if
    end do  
  end if   
  if ( itype==-1 ) then !just for restart file
    !do n=1,maxtile
      !write(lname,'("Sensible correction term ",I1.1)') n
      !write(vname,'("t",I1.1,"_fhscor")') n
      !call attrib(idnc,jdim,jsize,vname,lname,'W/m2',-3000.,3000.,any_m,point_m,land_m,double_m)
      !write(lname,'("Latent correction term ",I1.1)') n
      !write(vname,'("t",I1.1,"_fescor")') n
      !call attrib(idnc,jdim,jsize,vname,lname,'W/m2',-3000.,3000.,any_m,point_m,land_m,double_m)
    !end do
    lname='DIR VIS albedo'
    vname='albvisdir'
    call attrib(idnc,jdim,jsize,vname,lname,'none',0.,1.3,any_m,point_m,land_m,double_m)
    lname='DIF VIS albedo'
    vname='albvisdif'
    call attrib(idnc,jdim,jsize,vname,lname,'none',0.,1.3,any_m,point_m,land_m,double_m)
    lname='DIR NIR albedo'
    vname='albnirdir'
    call attrib(idnc,jdim,jsize,vname,lname,'none',0.,1.3,any_m,point_m,land_m,double_m)
    lname='DIF NIR albedo'
    vname='albnirdif'
    call attrib(idnc,jdim,jsize,vname,lname,'none',0.,1.3,any_m,point_m,land_m,double_m)
    lname='VIS albedo'
    vname='albvis'
    call attrib(idnc,jdim,jsize,vname,lname,'none',0.,1.3,any_m,point_m,land_m,double_m)
    lname='NIR albedo'
    vname='albnir'
    call attrib(idnc,jdim,jsize,vname,lname,'none',0.,1.3,any_m,point_m,land_m,double_m)
  end if
end if
  
return
end subroutine savetiledef

! *************************************************************************************
! This subroutine saves CABLE tile data
subroutine savetile(idnc,local,iarch,itype)

!use carbpools_m
use cc_mpi, only : ccmpi_abort
use infile
use newmpar_m
use parm_m, only : diaglevel_pop
use soil_m
use soilsnow_m
use vegpar_m
  
integer, intent(in) :: idnc,iarch,itype
integer k,n,np_pop,iq
integer cc,dd,hh,ll
integer tile
integer, dimension(ifull) :: dati
real, dimension(ifull) :: datr
real(kind=8), dimension(ifull) :: dat
real(kind=8), dimension(mp_max) :: dummy_unpack
real(kind=8), dimension(:), allocatable :: dat_in
real(kind=8), dimension(:,:), allocatable :: datpatch
real(kind=8), dimension(:,:), allocatable :: datage
real(kind=8), dimension(:,:,:), allocatable :: datpc
character(len=80) vname
logical, intent(in) :: local
  
if ( itype==-1 ) then !just for restart file
  do n = 1,maxtile
    dat = 0
    do tile = 1,ntiles  
      if ( tdata(tile)%mp>0 ) then  
        dummy_unpack(1:tdata(tile)%mp) = real(tdata(tile)%cveg,8)    
        call cable_unpack(dummy_unpack(1:tdata(tile)%mp),dat,tile,nb=n)
      end if  
    end do  
    write(vname,'("t",I1.1,"_cvc")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
  end do      
  
  ! soil temperature check
  do tile = 1,ntiles
    if ( tdata(tile)%mp>0 ) then
      if ( any(ssnow(tile)%tgg>425.) ) then
        write(6,*) "ERROR: Invalid CABLE temperature when writing tile"
        write(6,*) "ssnow%tgg ",maxval(ssnow(tile)%tgg)
        stop -1
      end if
    end if
  end do  
  
  do n = 1,maxtile  ! tile
    datr = 0.
    do tile = 1,ntiles  
      if ( tdata(tile)%mp>0 ) then  
         call cable_unpack(tdata(tile)%sv,datr,tile,nb=n)
      end if
    end do
    dat = real( datr, 8 )
    write(vname,'("t",I1.1,"_svs")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
    do k = 1,cbm_ms     ! soil layer
      dat = real(tgg(:,k),8)
      do tile = 1,ntiles  
        if ( tdata(tile)%mp>0 ) then  
          call cable_unpack(ssnow(tile)%tgg(:,k),dat,tile,nb=n)
        end if
      end do
      write(vname,'("t",I1.1,"_tgg",I1.1)') n,k
      call histwrt(dat,vname,idnc,iarch,local,.true.)
    end do
    do k = 1,cbm_ms
      dat = real(wb(:,k),8)
      do tile = 1,ntiles  
        if ( tdata(tile)%mp>0 ) then  
          call cable_unpack(ssnow(tile)%wb(:,k),dat,tile,nb=n)
        end if
      end do
      write(vname,'("t",I1.1,"_wb",I1.1)') n,k
      call histwrt(dat,vname,idnc,iarch,local,.true.)
    end do
    do k = 1,cbm_ms
      dat = real(wbice(:,k),8)
      do tile = 1,ntiles  
        if ( tdata(tile)%mp>0 ) then  
          call cable_unpack(ssnow(tile)%wbice(:,k),dat,tile,nb=n)
        end if
      end do
      write(vname,'("t",I1.1,"_wbice",I1.1)') n,k
      call histwrt(dat,vname,idnc,iarch,local,.true.)
    end do
    do k = 1,3 ! snow layer
      dat = real(tggsn(:,k),8)
      do tile = 1,ntiles  
        if ( tdata(tile)%mp>0 ) then  
          call cable_unpack(ssnow(tile)%tggsn(:,k),dat,tile,nb=n)
        end if
      end do
      write(vname,'("t",I1.1,"_tggsn",I1.1)') n,k
      call histwrt(dat,vname,idnc,iarch,local,.true.)
    end do
    do k = 1,3
      dat = real(smass(:,k),8)
      do tile = 1,ntiles  
        if ( tdata(tile)%mp>0 ) then  
          call cable_unpack(ssnow(tile)%smass(:,k),dat,tile,nb=n)
        end if
      end do
      write(vname,'("t",I1.1,"_smass",I1.1)') n,k
      call histwrt(dat,vname,idnc,iarch,local,.true.)
    end do
    do k = 1,3
      dat = real(ssdn(:,k),8)
      do tile = 1,ntiles  
        if ( tdata(tile)%mp>0 ) then  
          call cable_unpack(ssnow(tile)%ssdn(:,k),dat,tile,nb=n)
        end if
      end do
      write(vname,'("t",I1.1,"_ssdn",I1.1)') n,k
      call histwrt(dat,vname,idnc,iarch,local,.true.)
    end do  
    do k = 1,3
      dat = real(snowd/3.,8)
      do tile = 1,ntiles  
        if ( tdata(tile)%mp>0 ) then
          call cable_unpack(ssnow(tile)%sdepth(:,k),dat,tile,nb=n)
        end if
      end do
      write(vname,'("t",I1.1,"_sdepth",I1.1)') n,k
      call histwrt(dat,vname,idnc,iarch,local,.true.)
    end do
    do k = 1,3
      dat = 0.2_8
      do tile = 1,ntiles  
        if ( tdata(tile)%mp>0 ) then
          call cable_unpack(ssnow(tile)%sconds(:,k),dat,tile,nb=n)
        end if
      end do
      write(vname,'("t",I1.1,"_sconds",I1.1)') n,k
      call histwrt(dat,vname,idnc,iarch,local,.true.)
    end do
    dat = real(ssdnn,8)
    do tile = 1,ntiles  
      if ( tdata(tile)%mp>0 ) then
        call cable_unpack(ssnow(tile)%ssdnn,dat,tile,nb=n)
      end if
    end do
    write(vname,'("t",I1.1,"_ssdnn")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
    dat = real(isflag,8)
    do tile = 1,ntiles  
      if ( tdata(tile)%mp>0 ) then
        dummy_unpack(1:tdata(tile)%mp) = real(ssnow(tile)%isflag,8)    
        call cable_unpack(dummy_unpack(1:tdata(tile)%mp),dat,tile,nb=n)
      end if    
    end do  
    write(vname,'("t",I1.1,"_sflag")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
    dat = real(snowd,8)
    do tile = 1,ntiles  
      if ( tdata(tile)%mp>0 ) then
        call cable_unpack(ssnow(tile)%snowd,dat,tile,nb=n)
      end if
    end do
    write(vname,'("t",I1.1,"_snd")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
    dat = real(snowd,8)
    do tile = 1,ntiles  
      if ( tdata(tile)%mp>0 ) then
        call cable_unpack(ssnow(tile)%osnowd,dat,tile,nb=n)
      end if
    end do
    write(vname,'("t",I1.1,"_osnd")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
    dat = real(snage,8)
    do tile = 1,ntiles  
      if ( tdata(tile)%mp>0 ) then
        call cable_unpack(ssnow(tile)%snage,dat,tile,nb=n)
      end if
    end do
    write(vname,'("t",I1.1,"_snage")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
    dat = 100._8
    do tile = 1,ntiles  
      if ( tdata(tile)%mp>0 ) then
        call cable_unpack(ssnow(tile)%rtsoil,dat,tile,nb=n)
      end if
    end do  
    write(vname,'("t",I1.1,"_rtsoil")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
    dat = 0._8
    do tile = 1,ntiles  
      if ( tdata(tile)%mp>0 ) then
        call cable_unpack(ssnow(tile)%GWwb,dat,tile,nb=n)
      end if
    end do
    write(vname,'("t",I1.1,"_GWwb")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
    dat = 0._8
    do tile = 1,ntiles  
      if ( tdata(tile)%mp>0 ) then
        call cable_unpack(ssnow(tile)%wtd,dat,tile,nb=n)
      end if   
    end do  
    write(vname,'("t",I1.1,"_wtd")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
    dat=0._8
    do tile = 1,ntiles  
      if ( tdata(tile)%mp>0 ) then
        call cable_unpack(canopy(tile)%cansto,dat,tile,nb=n)
      end if
    end do
    write(vname,'("t",I1.1,"_cansto")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
    dat=0.01_8 ! ustar
    do tile = 1,ntiles  
      if ( tdata(tile)%mp>0 ) then
        call cable_unpack(canopy(tile)%us,dat,tile,nb=n)
      end if
    end do
    write(vname,'("t",I1.1,"_us")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)  
    dat=0._8
    do tile = 1,ntiles  
      if ( tdata(tile)%mp>0 ) then
        call cable_unpack(ssnow(tile)%pudsto,dat,tile,nb=n)
      end if
    end do
    write(vname,'("t",I1.1,"_pudsto")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
    dat=0._8
    do tile = 1,ntiles  
      if ( tdata(tile)%mp>0 ) then
        call cable_unpack(ssnow(tile)%wetfac,dat,tile,nb=n)
      end if
    end do  
    write(vname,'("t",I1.1,"_wetfac")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
    dat=0._8
    do tile = 1,ntiles  
      if ( tdata(tile)%mp>0 ) then
        call cable_unpack(canopy(tile)%ga,dat,tile,nb=n)
      end if
    end do  
    write(vname,'("t",I1.1,"_ga")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
  end do
  if ( soil_struc==1 ) then
    do n = 1,maxtile  ! tile
      dat=0._8
      do tile = 1,ntiles  
        if ( tdata(tile)%mp>0 ) then
          call cable_unpack(ssnow(tile)%h0,dat,tile,nb=n)
        end if
      end do  
      write(vname,'("t",I1.1,"_hzero")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)   
      do k = 1,cbm_ms     ! soil layer
        dat=0._8
        do tile = 1,ntiles  
          if ( tdata(tile)%mp>0 ) then
            call cable_unpack(ssnow(tile)%S(:,k),dat,tile,nb=n)
          end if
        end do  
        write(vname,'("t",I1.1,"_s",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do k = 1,cbm_ms
        dat=0._8
        do tile = 1,ntiles  
          if ( tdata(tile)%mp>0 ) then
            call cable_unpack(ssnow(tile)%tsoil(:,k),dat,tile,nb=n)
          end if
        end do
        write(vname,'("t",I1.1,"_tsoil",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do k = 1,cbm_ms
        dat=0._8
        do tile = 1,ntiles  
          if ( tdata(tile)%mp>0 ) then        
            call cable_unpack(ssnow(tile)%thetai(:,k),dat,tile,nb=n)
          end if
        end do  
        write(vname,'("t",I1.1,"_thetai",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      dat=0._8
      do tile = 1,ntiles  
        if ( tdata(tile)%mp>0 ) then
          call cable_unpack(ssnow(tile)%snowliq(:,1),dat,tile,nb=n)
        end if
      end do  
      write(vname,'("t",I1.1,"_snowliq",I1.1)') n,1
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      do tile = 1,ntiles  
        if ( tdata(tile)%mp>0 ) then
          call cable_unpack(ssnow(tile)%tsurface,dat,tile,nb=n)
        end if
      end do  
      write(vname,'("t",I1.1,"_tsurface")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      do tile = 1,ntiles  
        if ( tdata(tile)%mp>0 ) then
          dummy_unpack(1:tdata(tile)%mp) = real(ssnow(tile)%nsnow,8)  
          call cable_unpack(dummy_unpack(1:tdata(tile)%mp),dat,tile,nb=n)
        end if  
      end do  
      write(vname,'("t",I1.1,"_nsnow")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)   
      dat=0._8
      do tile = 1,ntiles  
        if ( tdata(tile)%mp>0 ) then
          call cable_unpack(canopy(tile)%fwsoil,dat,tile,nb=n)
        end if
      end do  
      write(vname,'("t",I1.1,"_fwsoil")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.) 
    end do
  end if
  if ( ccycle/=0 ) then
    do n = 1,maxtile  ! tile
      do k = 1,mplant     
        dat=0._8
        do tile = 1,ntiles  
          if ( tdata(tile)%mp>0 ) then
            call cable_unpack(casapool(tile)%cplant(:,k),dat,tile,nb=n)
          end if
        end do  
        write(vname,'("t",I1.1,"_cplant",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do k = 1,mplant
        dat=0._8
        do tile = 1,ntiles  
          if ( tdata(tile)%mp>0 ) then        
            call cable_unpack(casapool(tile)%nplant(:,k),dat,tile,nb=n)
          end if
        end do
        write(vname,'("t",I1.1,"_nplant",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do k = 1,mplant
        dat=0._8
        do tile = 1,ntiles  
          if ( tdata(tile)%mp>0 ) then
            call cable_unpack(casapool(tile)%pplant(:,k),dat,tile,nb=n)
          end if
        end do  
        write(vname,'("t",I1.1,"_pplant",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do k = 1,mlitter
        dat=0._8
        do tile = 1,ntiles  
          if ( tdata(tile)%mp>0 ) then
            call cable_unpack(casapool(tile)%clitter(:,k),dat,tile,nb=n)
          end if
        end do
        write(vname,'("t",I1.1,"_clitter",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do k = 1,mlitter
        dat=0._8
        do tile = 1,ntiles  
          if ( tdata(tile)%mp>0 ) then
            call cable_unpack(casapool(tile)%nlitter(:,k),dat,tile,nb=n)
          end if
        end do  
        write(vname,'("t",I1.1,"_nlitter",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do  
      do k = 1,mlitter
        dat=0._8
        do tile = 1,ntiles  
          if ( tdata(tile)%mp>0 ) then
            call cable_unpack(casapool(tile)%plitter(:,k),dat,tile,nb=n)
          end if
        end do  
        write(vname,'("t",I1.1,"_plitter",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do k = 1,msoil
        dat=0._8
        do tile = 1,ntiles  
          if ( tdata(tile)%mp>0 ) then
            call cable_unpack(casapool(tile)%csoil(:,k),dat,tile,nb=n)
          end if
        end do  
        write(vname,'("t",I1.1,"_csoil",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do k = 1,msoil
        dat=0._8
        do tile = 1,ntiles  
          if ( tdata(tile)%mp>0 ) then        
            call cable_unpack(casapool(tile)%nsoil(:,k),dat,tile,nb=n)
          end if
        end do
        write(vname,'("t",I1.1,"_nsoil",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do k = 1,msoil
        dat=0._8
        do tile = 1,ntiles  
          if ( tdata(tile)%mp>0 ) then        
            call cable_unpack(casapool(tile)%psoil(:,k),dat,tile,nb=n)
          end if
        end do  
        write(vname,'("t",I1.1,"_psoil",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      dat=0._8
      do tile = 1,ntiles  
        if ( tdata(tile)%mp>0 ) then
          call cable_unpack(casamet(tile)%glai(:),dat,tile,nb=n)
        end if
      end do  
      write(vname,'("t",I1.1,"_glai")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      do tile = 1,ntiles  
        if ( tdata(tile)%mp>0 ) then
          call cable_unpack(phen(tile)%phen(:),dat,tile,nb=n)
        end if
      end do  
      write(vname,'("t",I1.1,"_phen")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      do tile = 1,ntiles  
        if ( tdata(tile)%mp>0 ) then
          call cable_unpack(phen(tile)%aphen(:),dat,tile,nb=n)
        end if
      end do  
      write(vname,'("t",I1.1,"_aphen")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      do tile = 1,ntiles  
        if ( tdata(tile)%mp>0 ) then
          dummy_unpack(1:tdata(tile)%mp) = real(phen(tile)%phase(:),8)   
          call cable_unpack(dummy_unpack(1:tdata(tile)%mp),dat,tile,nb=n)
        end if  
      end do  
      write(vname,'("t",I1.1,"_phenphase")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      do tile = 1,ntiles  
        if ( tdata(tile)%mp>0 ) then
          dummy_unpack(1:tdata(tile)%mp) = real(phen(tile)%doyphase(:,3),8)    
          call cable_unpack(dummy_unpack(1:tdata(tile)%mp),dat,tile,nb=n)
        end if  
      end do  
      write(vname,'("t",I1.1,"_doyphase3")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      do tile = 1,ntiles  
        if ( tdata(tile)%mp>0 ) then
          call cable_unpack(casapool(tile)%clabile(:),dat,tile,nb=n)
        end if
      end do  
      write(vname,'("t",I1.1,"_clabile")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      do tile = 1,ntiles  
        if ( tdata(tile)%mp>0 ) then
          call cable_unpack(casapool(tile)%nsoilmin(:),dat,tile,nb=n)
        end if
      end do  
      write(vname,'("t",I1.1,"_nsoilmin")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      do tile = 1,ntiles  
        if ( tdata(tile)%mp>0 ) then      
          call cable_unpack(casapool(tile)%psoillab(:),dat,tile,nb=n)
        end if
      end do  
      write(vname,'("t",I1.1,"_psoillab")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      do tile = 1,ntiles  
        if ( tdata(tile)%mp>0 ) then
          call cable_unpack(casapool(tile)%psoilsorb(:),dat,tile,nb=n)
        end if
      end do  
      write(vname,'("t",I1.1,"_psoilsorb")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      do tile = 1,ntiles  
        if ( tdata(tile)%mp>0 ) then
          call cable_unpack(casapool(tile)%psoilocc(:),dat,tile,nb=n)
        end if
      end do  
      write(vname,'("t",I1.1,"_psoilocc")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      do k = 1,mplant     
        dat=0._8
        do tile = 1,ntiles  
          if ( tdata(tile)%mp>0 ) then
            call cable_unpack(casaflux(tile)%crmplant(:,k),dat,tile,nb=n)
          end if
        end do  
        write(vname,'("t",I1.1,"_crmplant",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      dat=0._8
      do tile = 1,ntiles  
        if ( tdata(tile)%mp>0 ) then
          call cable_unpack(casaflux(tile)%frac_sapwood(:),dat,tile,nb=n)
        end if
      end do  
      write(vname,'("t",I1.1,"_fracsapwood")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      do tile = 1,ntiles  
        if ( tdata(tile)%mp>0 ) then      
          call cable_unpack(casaflux(tile)%sapwood_area(:),dat,tile,nb=n)
        end if
      end do  
      write(vname,'("t",I1.1,"_sapwoodarea")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      do tile = 1,ntiles  
        if ( tdata(tile)%mp>0 ) then      
          call cable_unpack(casaflux(tile)%Crsoil(:),dat,tile,nb=n)
        end if
      end do  
      write(vname,'("t",I1.1,"_crsoil")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      do tile = 1,ntiles  
        if ( tdata(tile)%mp>0 ) then
          call cable_unpack(casaflux(tile)%cnpp(:),dat,tile,nb=n)
        end if
      end do  
      write(vname,'("t",I1.1,"_cnpp")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      do tile = 1,ntiles  
        if ( tdata(tile)%mp>0 ) then      
          call cable_unpack(casaflux(tile)%clabloss(:),dat,tile,nb=n)
        end if
      end do  
      write(vname,'("t",I1.1,"_clabloss")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      do tile = 1,ntiles  
        if ( tdata(tile)%mp>0 ) then
          call cable_unpack(casaflux(tile)%crgplant(:),dat,tile,nb=n)
        end if
      end do  
      write(vname,'("t",I1.1,"_crgplant")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      do tile = 1,ntiles  
        if ( tdata(tile)%mp>0 ) then      
          call cable_unpack(casaflux(tile)%stemnpp(:),dat,tile,nb=n)
        end if
      end do  
      write(vname,'("t",I1.1,"_stemnpp")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      do tile = 1,ntiles  
        if ( tdata(tile)%mp>0 ) then      
          call cable_unpack(casabal(tile)%laimax(:),dat,tile,nb=n)
        end if
      end do  
      write(vname,'("t",I1.1,"_LAImax")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      do tile = 1,ntiles  
        if ( tdata(tile)%mp>0 ) then
          call cable_unpack(casabal(tile)%cleafmean(:),dat,tile,nb=n)
        end if
      end do  
      write(vname,'("t",I1.1,"_Cleafmean")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      do tile = 1,ntiles  
        if ( tdata(tile)%mp>0 ) then      
          call cable_unpack(casabal(tile)%crootmean(:),dat,tile,nb=n)
        end if
      end do  
      write(vname,'("t",I1.1,"_Crootmean")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      do tile = 1,ntiles  
        if ( tdata(tile)%mp>0 ) then
          call cable_unpack(canopy(tile)%fpn(:),dat,tile,nb=n)
        end if
      end do  
      write(vname,'("t",I1.1,"_fpn")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      do tile = 1,ntiles  
        if ( tdata(tile)%mp>0 ) then      
          call cable_unpack(canopy(tile)%frday(:),dat,tile,nb=n)
        end if
      end do  
      write(vname,'("t",I1.1,"_frday")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
    end do
  end if
end if ! itype==-1
if ( cable_pop==1 ) then
  allocate( datpatch(ifull,POP_NPATCH) )  
  allocate( datage(ifull,POP_AGEMAX) )  
  allocate( datpc(ifull,POP_NPATCH,POP_NCOHORT) )
  !np_pop = size( pop%pop_grid(:) )
  allocate( dat_in(mp_pop_max) )
  do n = 1,maxtile
    datpatch = 0._8
    datage = 0._8
    datpc = 0._8
    dat = 0._8
    dat_in = 0._8
    if ( diaglevel_pop > 0 ) then
      do tile = 1,ntiles  
        np_pop = tdata(tile)%np  
        if ( np_pop>0 ) then        
          dat_in(1:np_pop) = pop(tile)%pop_grid(:)%cmass_sum 
          call pop_unpack(dat_in(1:np_pop),dat,tile,nb=n)
        end if  
      end do  
      write(vname,'("t",I1.1,"_pop_grid_cmass_sum")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
    end if
    if ( itype==-1 ) then !just for restart file
      do tile = 1,ntiles  
        np_pop = tdata(tile)%np  
        if ( np_pop>0 ) then        
          dat_in(1:np_pop) = pop(tile)%pop_grid(:)%cmass_sum_old  
          call pop_unpack(dat_in(1:np_pop),dat,tile,nb=n)
        end if  
      end do  
      write(vname,'("t",I1.1,"_pop_grid_cmass_sum_old")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      do tile = 1,ntiles  
        np_pop = tdata(tile)%np  
        if ( np_pop>0 ) then        
          dat_in(1:np_pop) = pop(tile)%pop_grid(:)%cheartwood_sum 
          call pop_unpack(dat_in(1:np_pop),dat,tile,nb=n)
        end if  
      end do  
      write(vname,'("t",I1.1,"_pop_grid_cheartwood_sum")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      do tile = 1,ntiles  
        np_pop = tdata(tile)%np  
        if ( np_pop>0 ) then        
          dat_in(1:np_pop) = pop(tile)%pop_grid(:)%csapwood_sum
          call pop_unpack(dat_in(1:np_pop),dat,tile,nb=n)
        end if  
      end do  
      write(vname,'("t",I1.1,"_pop_grid_csapwood_sum")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      do tile = 1,ntiles  
        np_pop = tdata(tile)%np  
        if ( np_pop>0 ) then        
          dat_in(1:np_pop) = pop(tile)%pop_grid(:)%csapwood_sum_old  
          call pop_unpack(dat_in(1:np_pop),dat,tile,nb=n)
        end if  
      end do  
      write(vname,'("t",I1.1,"_pop_grid_csapwood_sum_old")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      do tile = 1,ntiles  
        np_pop = tdata(tile)%np  
        if ( np_pop>0 ) then        
          dat_in(1:np_pop) = pop(tile)%pop_grid(:)%densindiv  
          call pop_unpack(dat_in(1:np_pop),dat,tile,nb=n)
        end if  
      end do  
      write(vname,'("t",I1.1,"_pop_grid_densindiv")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      do tile = 1,ntiles  
        np_pop = tdata(tile)%np  
        if ( np_pop>0 ) then
          dat_in(1:np_pop) = pop(tile)%pop_grid(:)%height_mean  
          call pop_unpack(dat_in(1:np_pop),dat,tile,nb=n)
        end if  
      end do  
      write(vname,'("t",I1.1,"_pop_grid_height_mean")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      do tile = 1,ntiles  
        np_pop = tdata(tile)%np  
        if ( np_pop>0 ) then
          dat_in(1:np_pop) = pop(tile)%pop_grid(:)%height_max  
          call pop_unpack(dat_in(1:np_pop),dat,tile,nb=n)
        end if  
      end do  
      write(vname,'("t",I1.1,"_pop_grid_height_max")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      do tile = 1,ntiles  
        np_pop = tdata(tile)%np  
        if ( np_pop>0 ) then
          dat_in(1:np_pop) = pop(tile)%pop_grid(:)%basal_area  
          call pop_unpack(dat_in(1:np_pop),dat,tile,nb=n)
        end if  
      end do  
      write(vname,'("t",I1.1,"_pop_grid_basal_area")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      do tile = 1,ntiles  
        np_pop = tdata(tile)%np  
        if ( np_pop>0 ) then
          dat_in(1:np_pop) = pop(tile)%pop_grid(:)%sapwood_loss  
          call pop_unpack(dat_in(1:np_pop),dat,tile,nb=n)
        end if  
      end do  
      write(vname,'("t",I1.1,"_pop_grid_sapwood_loss")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      do tile = 1,ntiles  
        np_pop = tdata(tile)%np  
        if ( np_pop>0 ) then
          dat_in(1:np_pop) = pop(tile)%pop_grid(:)%sapwood_area_loss  
          call pop_unpack(dat_in(1:np_pop),dat,tile,nb=n)
        end if
      end do  
      write(vname,'("t",I1.1,"_pop_grid_sapwood_area_loss")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      do tile = 1,ntiles  
        np_pop = tdata(tile)%np  
        if ( np_pop>0 ) then
          dat_in(1:np_pop) = pop(tile)%pop_grid(:)%stress_mortality  
          call pop_unpack(dat_in(1:np_pop),dat,tile,nb=n)
        end if  
      end do  
      write(vname,'("t",I1.1,"_pop_grid_stress_mortality")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      do tile = 1,ntiles  
        np_pop = tdata(tile)%np  
        if ( np_pop>0 ) then
          dat_in(1:np_pop) = pop(tile)%pop_grid(:)%crowding_mortality  
          call pop_unpack(dat_in(1:np_pop),dat,tile,nb=n)
        end if
      end do  
      write(vname,'("t",I1.1,"_pop_grid_crowding_mortality")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      do tile = 1,ntiles  
        np_pop = tdata(tile)%np  
        if ( np_pop>0 ) then
          dat_in(1:np_pop) = pop(tile)%pop_grid(:)%fire_mortality  
          call pop_unpack(dat_in(1:np_pop),dat,tile,nb=n)
        end if
      end do  
      write(vname,'("t",I1.1,"_pop_grid_fire_mortality")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      do tile = 1,ntiles  
        np_pop = tdata(tile)%np  
        if ( np_pop>0 ) then
          dat_in(1:np_pop) = pop(tile)%pop_grid(:)%cat_mortality  
          call pop_unpack(dat_in(1:np_pop),dat,tile,nb=n)
        end if  
      end do  
      write(vname,'("t",I1.1,"_pop_grid_cat_mortality")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      do tile = 1,ntiles  
        np_pop = tdata(tile)%np  
        if ( np_pop>0 ) then
          dat_in(1:np_pop) = pop(tile)%pop_grid(:)%res_mortality  
          call pop_unpack(dat_in(1:np_pop),dat,tile,nb=n)
        end if  
      end do  
      write(vname,'("t",I1.1,"_pop_grid_res_mortality")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      do tile = 1,ntiles  
        np_pop = tdata(tile)%np  
        if ( np_pop>0 ) then
          dat_in(1:np_pop) = pop(tile)%pop_grid(:)%growth  
          call pop_unpack(dat_in(1:np_pop),dat,tile,nb=n)
        end if  
      end do  
      write(vname,'("t",I1.1,"_pop_grid_growth")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      do tile = 1,ntiles  
        np_pop = tdata(tile)%np  
        if ( np_pop>0 ) then
          dat_in(1:np_pop) = pop(tile)%pop_grid(:)%area_growth  
          call pop_unpack(dat_in(1:np_pop),dat,tile,nb=n)
        end if  
      end do  
      write(vname,'("t",I1.1,"_pop_grid_area_growth")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      do tile = 1,ntiles  
        np_pop = tdata(tile)%np  
        if ( np_pop>0 ) then
          dat_in(1:np_pop) = pop(tile)%pop_grid(:)%crown_cover  
          call pop_unpack(dat_in(1:np_pop),dat,tile,nb=n)
        end if  
      end do  
      write(vname,'("t",I1.1,"_pop_grid_crown_cover")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      do tile = 1,ntiles  
        np_pop = tdata(tile)%np  
        if ( np_pop>0 ) then
          dat_in(1:np_pop) = pop(tile)%pop_grid(:)%crown_area  
          call pop_unpack(dat_in(1:np_pop),dat,tile,nb=n)
        end if  
      end do  
      write(vname,'("t",I1.1,"_pop_grid_crown_area")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      do tile = 1,ntiles  
        np_pop = tdata(tile)%np  
        if ( np_pop>0 ) then
          dat_in(1:np_pop) = pop(tile)%pop_grid(:)%crown_volume  
          call pop_unpack(dat_in(1:np_pop),dat,tile,nb=n)
        end if  
      end do  
      write(vname,'("t",I1.1,"_pop_grid_crown_volume")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      do tile = 1,ntiles  
        np_pop = tdata(tile)%np  
        if ( np_pop>0 ) then
          dat_in(1:np_pop) = pop(tile)%pop_grid(:)%sapwood_area  
          call pop_unpack(dat_in(1:np_pop),dat,tile,nb=n)
        end if  
      end do  
      write(vname,'("t",I1.1,"_pop_grid_sapwood_area")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      do tile = 1,ntiles  
        np_pop = tdata(tile)%np  
        if ( np_pop>0 ) then
          dat_in(1:np_pop) = pop(tile)%pop_grid(:)%sapwood_area_old  
          call pop_unpack(dat_in(1:np_pop),dat,tile,nb=n)
        end if  
      end do  
      write(vname,'("t",I1.1,"_pop_grid_sapwood_area_old")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      do tile = 1,ntiles  
        np_pop = tdata(tile)%np  
        if ( np_pop>0 ) then
          dat_in(1:np_pop) = pop(tile)%pop_grid(:)%KClump  
          call pop_unpack(dat_in(1:np_pop),dat,tile,nb=n)
        end if  
      end do  
      write(vname,'("t",I1.1,"_pop_grid_KClump")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      do k = 1,POP_AGEMAX
        do tile = 1,ntiles  
          np_pop = tdata(tile)%np  
          if ( np_pop>0 ) then          
            dat_in(1:np_pop) = pop(tile)%pop_grid(:)%freq_age(k)  
            call pop_unpack(dat_in(1:np_pop),datage(:,k),tile,nb=n)
          end if  
        end do
      end do      
      write(vname,'("t",I1.1,"_pop_grid_freq_age")') n
      call histwrt(datage,vname,idnc,iarch,local,.true.)
      do k = 1,POP_AGEMAX
        do tile = 1,ntiles  
          np_pop = tdata(tile)%np  
          if ( np_pop>0 ) then
            dat_in(1:np_pop) = pop(tile)%pop_grid(:)%biomass_age(k)  
            call pop_unpack(dat_in(1:np_pop),datage(:,k),tile,n)
          end if  
        end do
      end do      
      write(vname,'("t",I1.1,"_pop_grid_biomass_age")') n
      call histwrt(datage,vname,idnc,iarch,local,.true.)
      do ll = 1,POP_NLAYER
        do tile = 1,ntiles  
          np_pop = tdata(tile)%np  
          if ( np_pop>0 ) then
            dat_in(1:np_pop) = pop(tile)%pop_grid(:)%biomass(ll)
            call pop_unpack(dat_in(1:np_pop),dat,tile,nb=n)
          end if  
        end do  
        write(vname,'("t",I1.1,"_pop_grid_biomass",I1.1)') n,ll  
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do tile = 1,ntiles  
          np_pop = tdata(tile)%np  
          if ( np_pop>0 ) then
            dat_in(1:np_pop) = pop(tile)%pop_grid(:)%density(ll)  
            call pop_unpack(dat_in(1:np_pop),dat,tile,nb=n)
          end if  
        end do  
        write(vname,'("t",I1.1,"_pop_grid_density",I1.1)') n,ll  
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do tile = 1,ntiles  
          np_pop = tdata(tile)%np  
          if ( np_pop>0 ) then
            dat_in(1:np_pop) = pop(tile)%pop_grid(:)%hmean(ll)  
            call pop_unpack(dat_in(1:np_pop),dat,tile,nb=n)
          end if  
        end do  
        write(vname,'("t",I1.1,"_pop_grid_hmean",I1.1)') n,ll  
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do tile = 1,ntiles  
          np_pop = tdata(tile)%np  
          if ( np_pop>0 ) then
            dat_in(1:np_pop) = pop(tile)%pop_grid(:)%hmax(ll)  
            call pop_unpack(dat_in(1:np_pop),dat,tile,nb=n)
          end if  
        end do  
        write(vname,'("t",I1.1,"_pop_grid_hmax",I1.1)') n,ll  
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do hh = 1,POP_HEIGHT_BINS
        do tile = 1,ntiles  
          np_pop = tdata(tile)%np  
          if ( np_pop>0 ) then
            dat_in(1:np_pop) = pop(tile)%pop_grid(:)%cmass_stem_bin(hh)  
            call pop_unpack(dat_in(1:np_pop),dat,tile,nb=n)
          end if
        end do  
        write(vname,'("t",I1.1,"_pop_grid_cmass_stem_bin",I2.2)') n,hh  
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do hh = 1,POP_HEIGHT_BINS
        do tile = 1,ntiles  
          np_pop = tdata(tile)%np  
          if ( np_pop>0 ) then
            dat_in(1:np_pop) = pop(tile)%pop_grid(:)%densindiv_bin(hh)  
            call pop_unpack(dat_in(1:np_pop),dat,tile,nb=n)
          end if  
        end do  
        write(vname,'("t",I1.1,"_pop_grid_densindiv_bin",I2.2)') n,hh  
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do hh = 1,POP_HEIGHT_BINS
        do tile = 1,ntiles  
          np_pop = tdata(tile)%np  
          if ( np_pop>0 ) then
            dat_in(1:np_pop) = pop(tile)%pop_grid(:)%height_bin(hh)  
            call pop_unpack(dat_in(1:np_pop),dat,tile,nb=n)
          end if  
        end do  
        write(vname,'("t",I1.1,"_pop_grid_height_bin",I2.2)') n,hh
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do hh = 1,POP_HEIGHT_BINS
        do tile = 1,ntiles  
          np_pop = tdata(tile)%np  
          if ( np_pop>0 ) then
            dat_in(1:np_pop) = pop(tile)%pop_grid(:)%diameter_bin(hh)  
            call pop_unpack(dat_in(1:np_pop),dat,tile,nb=n)
          end if  
        end do  
        write(vname,'("t",I1.1,"_pop_grid_diameter_bin",I2.2)') n,hh 
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do dd = 1,POP_NDISTURB
        do tile = 1,ntiles  
          np_pop = tdata(tile)%np  
          if ( np_pop>0 ) then
            dat_in(1:np_pop) = real(pop(tile)%pop_grid(:)%n_age(dd),8)
            call pop_unpack(dat_in(1:np_pop),dat,tile,nb=n)
          end if  
        end do  
        write(vname,'("t",I1.1,"_pop_grid_n_age",I1.1)') n,dd 
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
    end if
    if ( itype==-1 .or. diaglevel_pop>=9 ) then !just for restart file
      do k = 1,POP_NPATCH
        do tile = 1,ntiles  
          np_pop = tdata(tile)%np  
          if ( np_pop>0 ) then
            dat_in(1:np_pop) = real(pop(tile)%pop_grid(:)%patch(k)%id,8)
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),tile,nb=n)
          end if  
        end do
      end do 
      write(vname,'("t",I1.1,"_pop_grid_patch_id")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles  
          np_pop = tdata(tile)%np  
          if ( np_pop>0 ) then
            dat_in(1:np_pop) = pop(tile)%pop_grid(:)%freq(k)  
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),tile,nb=n)
          end if  
        end do
      end do
      write(vname,'("t",I1.1,"_pop_grid_freq")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles  
          np_pop = tdata(tile)%np  
          if ( np_pop>0 ) then
            dat_in(1:np_pop) = pop(tile)%pop_grid(:)%freq_old(k)   
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),tile,nb=n)
          end if  
        end do
      end do
      write(vname,'("t",I1.1,"_pop_grid_freq_old")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles  
          np_pop = tdata(tile)%np  
          if ( np_pop>0 ) then
            dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%factor_recruit  
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),tile,nb=n)
          end if  
        end do
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_factor_recruit")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles  
          np_pop = tdata(tile)%np  
          if ( np_pop>0 ) then
            dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%pgap  
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),tile,nb=n)
          end if  
        end do
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_pgap")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles  
          np_pop = tdata(tile)%np  
          if ( np_pop>0 ) then          
            dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%lai  
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),tile,nb=n)
          end if  
        end do
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_lai")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles  
          np_pop = tdata(tile)%np  
          if ( np_pop>0 ) then
            dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%biomass  
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),tile,nb=n)
          end if  
        end do
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_biomass")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles  
          np_pop = tdata(tile)%np  
          if ( np_pop>0 ) then          
            dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%biomass_old  
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),tile,nb=n)
          end if  
        end do
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_biomass_old")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)    
      do k = 1,POP_NPATCH
        do tile = 1,ntiles  
          np_pop = tdata(tile)%np  
          if ( np_pop>0 ) then
            dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%sapwood  
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),tile,nb=n)
          end if  
        end do
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles  
          np_pop = tdata(tile)%np  
          if ( np_pop>0 ) then          
            dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%heartwood  
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),tile,nb=n)
          end if  
        end do
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_heartwood")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles  
          np_pop = tdata(tile)%np  
          if ( np_pop>0 ) then
            dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%sapwood_old  
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),tile,nb=n)
          end if  
        end do
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_old")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles  
          np_pop = tdata(tile)%np  
          if ( np_pop>0 ) then
            dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%sapwood_area  
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),tile,nb=n)
          end if  
        end do
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_area")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles  
          np_pop = tdata(tile)%np  
          if ( np_pop>0 ) then
            dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%sapwood_area_old  
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),tile,nb=n)
          end if  
        end do
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_area_old")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles  
          np_pop = tdata(tile)%np  
          if ( np_pop>0 ) then
            dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%stress_mortality   
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),tile,nb=n)
          end if
        end do
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_stress_mortality")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles  
          np_pop = tdata(tile)%np  
          if ( np_pop>0 ) then
            dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%fire_mortality  
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),tile,nb=n)
          end if  
        end do
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_fire_mortality")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles  
          np_pop = tdata(tile)%np  
          if ( np_pop>0 ) then
            dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%cat_mortality  
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),tile,nb=n)
          end if  
        end do
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_cat_mortality")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles  
          np_pop = tdata(tile)%np  
          if ( np_pop>0 ) then
            dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%crowding_mortality 
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),tile,nb=n)
          end if  
        end do
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_crowding_mortality")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles  
          np_pop = tdata(tile)%np  
          if ( np_pop>0 ) then
            dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%cpc  
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),tile,nb=n)
          end if  
        end do
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_cpc")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles  
          np_pop = tdata(tile)%np  
          if ( np_pop>0 ) then
            dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%mortality  
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),tile,nb=n)
          end if  
        end do
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_mortality")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles  
          np_pop = tdata(tile)%np  
          if ( np_pop>0 ) then          
            dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%sapwood_loss  
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),tile,nb=n)
          end if  
        end do
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_loss")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles  
          np_pop = tdata(tile)%np  
          if ( np_pop>0 ) then
            dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%sapwood_area_loss  
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),tile,nb=n)
          end if  
        end do
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_area_loss")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles  
          np_pop = tdata(tile)%np  
          if ( np_pop>0 ) then          
            dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%growth  
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),tile,nb=n)
          end if  
        end do
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_growth")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles  
          np_pop = tdata(tile)%np  
          if ( np_pop>0 ) then          
            dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%area_growth  
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),tile,nb=n)
          end if  
        end do
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_area_growth")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles  
          np_pop = tdata(tile)%np  
          if ( np_pop>0 ) then
            dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%frac_NPP  
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),tile,nb=n)
          end if  
        end do
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_frac_NPP")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles  
          np_pop = tdata(tile)%np  
          if ( np_pop>0 ) then
            dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%frac_respiration  
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),tile,nb=n)
          end if  
        end do
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_frac_respiration")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do k = 1,POP_NPATCH
        do tile = 1,ntiles  
          np_pop = tdata(tile)%np  
          if ( np_pop>0 ) then
            dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%frac_light_uptake  
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),tile,nb=n)
          end if  
        end do
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_frac_light_uptake")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do dd = 1,POP_NDISTURB  
        do k = 1,POP_NPATCH  
          do tile = 1,ntiles  
            np_pop = tdata(tile)%np  
            if ( np_pop>0 ) then
              dat_in(1:np_pop) = real(pop(tile)%pop_grid(:)%patch(k)%disturbance_interval(dd),8)
              call pop_unpack(dat_in(1:np_pop),datpatch(:,k),tile,nb=n)
            end if  
          end do  
        end do
        write(vname,'("t",I1.1,"_pop_grid_patch_disturbance_interval",I1.1)') n,dd
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do dd = 1,POP_NDISTURB 
        do k = 1,POP_NPATCH  
          do tile = 1,ntiles  
            np_pop = tdata(tile)%np  
            if ( np_pop>0 ) then
              dat_in(1:np_pop) = real(pop(tile)%pop_grid(:)%patch(k)%first_disturbance_year(dd),8)
              call pop_unpack(dat_in(1:np_pop),datpatch(:,k),tile,nb=n)
            end if  
          end do  
        end do
        write(vname,'("t",I1.1,"_pop_grid_patch_first_disturbance_year",I1.1)') n,dd
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do dd = 1,POP_NDISTURB 
        do k = 1,POP_NPATCH  
          do tile = 1,ntiles  
            np_pop = tdata(tile)%np  
            if ( np_pop>0 ) then
              dat_in(1:np_pop) = real(pop(tile)%pop_grid(:)%patch(k)%age(dd),8)
              call pop_unpack(dat_in(1:np_pop),datpatch(:,k),tile,nb=n)
            end if  
          end do  
        end do
        write(vname,'("t",I1.1,"_pop_grid_patch_age",I1.1)') n,dd
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do dd = 1,POP_NDISTURB  
        do k = 1,POP_NPATCH  
          do tile = 1,ntiles  
            np_pop = tdata(tile)%np  
            if ( np_pop>0 ) then
              dat_in(1:np_pop) = real(pop(tile)%pop_grid(:)%ranked_age_unique(k,dd),8)  
              call pop_unpack(dat_in(1:np_pop),datpatch(:,k),tile,nb=n)
            end if  
          end do  
        end do
        write(vname,'("t",I1.1,"_pop_grid_ranked_age_unique",I1.1)') n,dd
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do dd = 1,POP_NDISTURB
        do k = 1,POP_NPATCH  
          do tile = 1,ntiles  
            np_pop = tdata(tile)%np  
            if ( np_pop>0 ) then
              dat_in(1:np_pop) = pop(tile)%pop_grid(:)%freq_ranked_age_unique(k,dd)  
              call pop_unpack(dat_in(1:np_pop),datpatch(:,k),tile,nb=n)
            end if  
          end do  
        end do
        write(vname,'("t",I1.1,"_pop_grid_freq_ranked_age_unique",I1.1)') n,dd
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do k = 1,POP_NPATCH  
          do tile = 1,ntiles  
            np_pop = tdata(tile)%np  
            if ( np_pop>0 ) then            
              dat_in(1:np_pop) = real(pop(tile)%pop_grid(:)%patch(k)%layer(ll)%ncohort,8)
              call pop_unpack(dat_in(1:np_pop),datpatch(:,k),tile,nb=n)
            end if  
          end do
        end do
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_ncohort")') n,ll
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do k = 1,POP_NPATCH  
          do tile = 1,ntiles  
            np_pop = tdata(tile)%np  
            if ( np_pop>0 ) then
              dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%layer(ll)%biomass  
              call pop_unpack(dat_in(1:np_pop),datpatch(:,k),tile,nb=n)
            end if  
          end do
        end do
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_biomass")') n,ll
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do k = 1,POP_NPATCH  
          do tile = 1,ntiles  
            np_pop = tdata(tile)%np  
            if ( np_pop>0 ) then
              dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%layer(ll)%density  
              call pop_unpack(dat_in(1:np_pop),datpatch(:,k),tile,nb=n)
            end if  
          end do
        end do
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_density")') n,ll
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do k = 1,POP_NPATCH  
          do tile = 1,ntiles  
            np_pop = tdata(tile)%np  
            if ( np_pop>0 ) then
              dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%layer(ll)%hmean 
              call pop_unpack(dat_in(1:np_pop),datpatch(:,k),tile,nb=n)
            end if  
          end do
        end do
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_hmean")') n,ll
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do k = 1,POP_NPATCH  
          do tile = 1,ntiles  
            np_pop = tdata(tile)%np  
            if ( np_pop>0 ) then
              dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%layer(ll)%hmax  
              call pop_unpack(dat_in(1:np_pop),datpatch(:,k),tile,nb=n)
            end if  
          end do
        end do
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_hmax")') n,ll
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT            
          do k = 1,POP_NPATCH  
            do tile = 1,ntiles  
              np_pop = tdata(tile)%np  
              if ( np_pop>0 ) then
                dat_in(1:np_pop) = real(pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%age,8)
                call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),tile,nb=n)
              end if  
            end do  
          end do  
        end do
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_age")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT 
          do k = 1,POP_NPATCH  
            do tile = 1,ntiles  
              np_pop = tdata(tile)%np  
              if ( np_pop>0 ) then
                dat_in(1:np_pop) = real(pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%id,8)
                call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),tile,nb=n)
              end if  
            end do  
          end do  
        end do
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_id")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT            
          do k = 1,POP_NPATCH  
            do tile = 1,ntiles  
              np_pop = tdata(tile)%np  
              if ( np_pop>0 ) then
                dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%biomass  
                call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),tile,nb=n)
              end if  
            end do  
          end do  
        end do
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_biomass")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH  
            do tile = 1,ntiles  
              np_pop = tdata(tile)%np  
              if ( np_pop>0 ) then
                dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%density  
                call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),tile,nb=n)
              end if  
            end do  
          end do  
        end do
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_density")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT  
          do k = 1,POP_NPATCH  
            do tile = 1,ntiles  
              np_pop = tdata(tile)%np  
              if ( np_pop>0 ) then
                dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_resource_uptake  
                call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),tile,nb=n)
              end if  
            end do  
          end do  
        end do
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_resource_uptake")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH  
            do tile = 1,ntiles  
              np_pop = tdata(tile)%np  
              if ( np_pop>0 ) then
                dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_light_uptake  
                call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),tile,nb=n)
              end if  
            end do  
          end do  
        end do
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_light_uptake")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT            
          do k = 1,POP_NPATCH  
            do tile = 1,ntiles  
              np_pop = tdata(tile)%np  
              if ( np_pop>0 ) then
                dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_interception  
                call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),tile,nb=n)
              end if  
            end do  
          end do  
        end do
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_interception")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH  
            do tile = 1,ntiles  
              np_pop = tdata(tile)%np  
              if ( np_pop>0 ) then
                dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_respiration 
                call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),tile,nb=n)
              end if  
            end do  
          end do  
        end do
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_respiration")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH  
            do tile = 1,ntiles  
              np_pop = tdata(tile)%np  
              if ( np_pop>0 ) then
                dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_NPP  
                call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),tile,nb=n)
              end if  
            end do
          end do  
        end do
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_NPP")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH
            do tile = 1,ntiles  
              np_pop = tdata(tile)%np  
              if ( np_pop>0 ) then
                dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%respiration_scalar  
                call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),tile,nb=n)
              end if  
            end do  
          end do  
        end do
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_respiration_scalar")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT            
          do k = 1,POP_NPATCH  
            do tile = 1,ntiles  
              np_pop = tdata(tile)%np  
              if ( np_pop>0 ) then
                dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%crown_area  
                call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),tile,nb=n)
              end if  
            end do  
          end do  
        end do
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_crown_area")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH  
            do tile = 1,ntiles  
              np_pop = tdata(tile)%np  
              if ( np_pop>0 ) then
                dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%Pgap  
                call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),tile,nb=n)
              end if  
            end do  
          end do  
        end do
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Pgap")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT  
          do k = 1,POP_NPATCH  
            do tile = 1,ntiles  
              np_pop = tdata(tile)%np  
              if ( np_pop>0 ) then
                dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%height  
                call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),tile,nb=n)
              end if  
            end do  
          end do  
        end do
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_height")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH 
            do tile = 1,ntiles  
              np_pop = tdata(tile)%np  
              if ( np_pop>0 ) then
                dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%diameter  
                call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),tile,nb=n)
              end if  
            end do 
          end do  
        end do
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_diameter")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH  
            do tile = 1,ntiles  
              np_pop = tdata(tile)%np  
              if ( np_pop>0 ) then
                dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%sapwood  
                call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),tile,nb=n)
              end if  
            end do  
          end do  
        end do
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_sapwood")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH  
            do tile = 1,ntiles  
              np_pop = tdata(tile)%np  
              if ( np_pop>0 ) then
                dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%heartwood  
                call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),tile,nb=n)
              end if  
            end do  
          end do  
        end do
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_heartwood")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT            
          do k = 1,POP_NPATCH
            do tile = 1,ntiles  
              np_pop = tdata(tile)%np  
              if ( np_pop>0 ) then
                dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%sapwood_area  
                call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),tile,nb=n)
              end if  
            end do  
          end do  
        end do
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_sapwood_area")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH  
            do tile = 1,ntiles  
              np_pop = tdata(tile)%np  
              if ( np_pop>0 ) then
                dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%basal_area  
                call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),tile,nb=n)
              end if  
            end do  
          end do  
        end do
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_basal_area")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH  
            do tile = 1,ntiles  
              np_pop = tdata(tile)%np  
              if ( np_pop>0 ) then
                dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%LAI  
                call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),tile,nb=n)
              end if  
            end do  
          end do  
        end do
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_LAI")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH 
            do tile = 1,ntiles  
              np_pop = tdata(tile)%np  
              if ( np_pop>0 ) then
                dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%Cleaf  
                call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),tile,nb=n)
              end if  
            end do  
          end do  
        end do
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Cleaf")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH  
            do tile = 1,ntiles  
              np_pop = tdata(tile)%np  
              if ( np_pop>0 ) then              
                dat_in(1:np_pop) = pop(tile)%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%Croot  
                call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),tile,nb=n)
              end if  
            end do  
          end do  
        end do  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Croot")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
    end if
  end do
  deallocate( datpatch )
  deallocate( datage )
  deallocate( datpc )
  deallocate( dat_in )
end if    
if ( itype==-1 ) then !just for restart file
  !do n = 1,maxtile  ! tile
    !dat=0._8
    !if (n<=maxnb) call cable_unpack(canopy%fhs_cor(:),dat,n)
    !write(vname,'("t",I1.1,"_fhscor")') n
    !call histwrt(dat,vname,idnc,iarch,local,.true.)
    !dat=0._8
    !if (n<=maxnb) call cable_unpack(canopy%fes_cor(:),dat,n)
    !write(vname,'("t",I1.1,"_fescor")') n
    !call histwrt(dat,vname,idnc,iarch,local,.true.)
    !dat=0._8
    !if (n<=maxnb) call cable_unpack(canopy%fns_cor(:),dat,n)
    !write(vname,'("t",I1.1,"_fnscor")') n
    !call histwrt(dat,vname,idnc,iarch,local,.true.)
    !dat=0._8
    !if (n<=maxnb) call cable_unpack(canopy%ga_cor(:),dat,n)
    !write(vname,'("t",I1.1,"_gacor")') n
    !call histwrt(dat,vname,idnc,iarch,local,.true.)
  !end do
  vname='albvisdir'
  call histwrt(albvisdir,vname,idnc,iarch,local,.true.)
  vname='albvisdif'
  call histwrt(albvisdif,vname,idnc,iarch,local,.true.)
  vname='albnirdir'
  call histwrt(albnirdir,vname,idnc,iarch,local,.true.)
  vname='albnirdif'
  call histwrt(albnirdif,vname,idnc,iarch,local,.true.)
  vname='albvis'
  call histwrt(albvissav,vname,idnc,iarch,local,.true.)
  vname='albnir'
  call histwrt(albnirsav,vname,idnc,iarch,local,.true.)
end if
  
return
end subroutine savetile 

end module cable_ccam_file

