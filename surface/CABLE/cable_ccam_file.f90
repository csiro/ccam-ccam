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

use carbpools_m
use cc_mpi
use darcdf_m
use infile
use newmpar_m
use parm_m
use pbl_m
use soil_m
use soilsnow_m
use vegpar_m
  
integer k
real, dimension(ifull) :: dummy_pack

if ( mp_global>0 ) then
  do k = 1,ms
    call cable_pack(tgg(:,k),ssnow%tgg(:,k))
    call cable_pack(wb(:,k),ssnow%wb(:,k))
    call cable_pack(wbice(:,k),ssnow%wbice(:,k))
  end do
  do k = 1,3
    call cable_pack(tggsn(:,k),ssnow%tggsn(:,k))
    call cable_pack(smass(:,k),ssnow%smass(:,k))
    call cable_pack(ssdn(:,k),ssnow%ssdn(:,k))
    dummy_pack = smass(:,k)/ssdn(:,k)
    call cable_pack(dummy_pack,ssnow%sdepth(:,k))
    ssnow%sconds(:,k) = 0.3_8
  end do      
  call cable_pack(tss,rad%trad(:))
  call cable_pack(ssdnn,ssnow%ssdnn)
  call cable_pack(isflag,ssnow%isflag)
  call cable_pack(snowd,ssnow%snowd)
  call cable_pack(snage,ssnow%snage)
  ssnow%Qrecharge = 0._8
  canopy%sublayer_dz = 0._8
  ssnow%rtevap_sat = 0._8
  ssnow%rtevap_unsat = 0._8
  ssnow%satfrac = 0.5_8
  ssnow%wbliq = ssnow%wb - ssnow%wbice
  ssnow%GWwb = 0.5_8*soil%ssat 
  ssnow%wtd = 20000._8
  dummy_pack = real(1-isflag)*tgg(:,1) + real(isflag)*tggsn(:,1) - 273.15
  call cable_pack(dummy_pack,ssnow%tsurface)
  ssnow%rtsoil = 50._8
  canopy%cansto = 0._8
  canopy%ga = 0._8
  canopy%us = 0.01_8
  ssnow%pudsto = 0._8
  ssnow%wetfac = 0._8
  ssnow%osnowd = ssnow%snowd
  canopy%fhs_cor = 0._8
  canopy%fes_cor = 0._8
  canopy%fns_cor = 0._8
  canopy%ga_cor = 0._8
  
  ! default value for fwsoil.  Recaculated by cable_canopy or by SLI
  canopy%fwsoil = max( 1.e-9_8, sum( veg%froot*max(1.e-9_8,min(1._8,ssnow%wb-spread(soil%swilt,2,ms))),2) &
      / ( soil%sfc-soil%swilt ) )
  
  call defaulttile_sli
  call fixtile
  
end if

return
end subroutine defaulttile

subroutine defaulttile_sli

use carbpools_m
use cc_mpi
use darcdf_m
use infile
use newmpar_m
use parm_m
use pbl_m
use soil_m
use soilsnow_m
use vegpar_m
  
integer k

if ( mp_global>0 ) then
    
  if ( soil_struc==1 ) then  
    
    ssnow%h0 = 0._8      ! pond height
    ssnow%snowliq = 0._8 ! liquid snow
    ssnow%Tsoil = ssnow%tgg - 273.15_8
    ssnow%thetai = ssnow%wbice
    do k = 1,ms
      ssnow%S(:,k) = ssnow%wb(:,k)/soil%ssat
    end do
    where ( ssnow%snowd>0. )
      ssnow%nsnow = 1
      ssnow%sdepth(:,1) = ssnow%snowd
      ssnow%smass(:,1) = ssnow%snowd*ssnow%ssdn(:,1)
    elsewhere
      ssnow%nsnow = 0
      ssnow%sdepth(:,1) = 0._8
      ssnow%smass(:,1) = 0._8
    end where
    !ssnow%nsnow = 0
    !ssnow%snowd = 0._8
  
  end if  
    
  call fixtile
  
end if

return
end subroutine defaulttile_sli

subroutine loadtile(usedefault)

use carbpools_m
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
integer, dimension(6) :: ierr_check
integer, dimension(ifull) :: dati
integer, dimension(mp_global) :: old_cv
integer, dimension(ifull,maxtile) :: nmp
integer, dimension(:), allocatable :: dati_out
real, dimension(ifull) :: datr
real, dimension(mp_global) :: dummy_unpack, old_sv
real(kind=8), dimension(ifull) :: dat
real(kind=8), dimension(ifull,ms) :: datms
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
if ( mp_global>0 ) then
  old_sv(:) = sv(:)
  old_cv(:) = cveg(:)
end if  

! check for changes
if ( ierr_cvc==0 ) then
  do n = 1,maxtile      
    write(vname,'("t",I1.1,"_cvc")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    dati = nint(dat)  
    call cable_pack(dati,old_cv,n)
  end do
end if
call create_new_tile_map(old_cv,nmp)

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
      call cable_pack(datr,old_sv,nmp(:,n))
    end if        
    write(vname,'("t",I1.1,"_tgg")') n
    call histrd(iarchi-1,ierr,vname,datms(:,1:ms),ifull)
    do k = 1,ms
      do iq = 1,ifull
        if ( land(iq) .and. (datms(iq,k)<100._8.or.datms(iq,k)>400._8) ) then
          ! change in land-sea mask?
          datms(iq,k) = tss(iq) ! use surface temperature
        end if
      end do  
      call cable_pack(datms(:,k),ssnow%tgg(:,k),nmp(:,n))
    end do
    write(vname,'("t",I1.1,"_wb")') n
    call histrd(iarchi-1,ierr,vname,datms(:,1:ms),ifull)
    do k = 1,ms
      call cable_pack(datms(:,k),ssnow%wb(:,k),nmp(:,n))
    end do
    write(vname,'("t",I1.1,"_wbice")') n
    call histrd(iarchi-1,ierr,vname,datms(:,1:ms),ifull)
    do k = 1,ms
      call cable_pack(datms(:,k),ssnow%wbice(:,k),nmp(:,n))
    end do
    write(vname,'("t",I1.1,"_tggsn")') n
    call histrd(iarchi-1,ierr,vname,dat3(:,1:3),ifull)
    do k = 1,3
      call cable_pack(dat3(:,k),ssnow%tggsn(:,k),nmp(:,n))
    end do
    write(vname,'("t",I1.1,"_smass")') n
    call histrd(iarchi-1,ierr,vname,dat3(:,1:3),ifull)
    do k = 1,3
      call cable_pack(dat3(:,k),ssnow%smass(:,k),nmp(:,n))
    end do
    write(vname,'("t",I1.1,"_ssdn")') n
    call histrd(iarchi-1,ierr,vname,dat3(:,1:3),ifull)
    do k = 1,3
      call cable_pack(dat3(:,k),ssnow%ssdn(:,k),nmp(:,n))
    end do
    write(vname,'("t",I1.1,"_sdepth",I1.1)') n
    call histrd(iarchi-1,ierr,vname,dat3(:,1:3),ifull)
    do k = 1,3
      call cable_pack(dat3(:,k),ssnow%sdepth(:,k),nmp(:,n))
    end do
    write(vname,'("t",I1.1,"_sconds")') n
    call histrd(iarchi-1,ierr,vname,dat3(:,1:3),ifull)
    do k = 1,3
      call cable_pack(dat3(:,k),ssnow%sconds(:,k),nmp(:,n))
    end do
    write(vname,'("t",I1.1,"_ssdnn")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    call cable_pack(dat,ssnow%ssdnn(:),nmp(:,n))
    write(vname,'("t",I1.1,"_sflag")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    dati = nint(dat)
    call cable_pack(dati,ssnow%isflag(:),nmp(:,n))
    write(vname,'("t",I1.1,"_snd")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    call cable_pack(dat,ssnow%snowd(:),nmp(:,n))
    write(vname,'("t",I1.1,"_osnd")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    call cable_pack(dat,ssnow%osnowd(:),nmp(:,n))
    write(vname,'("t",I1.1,"_snage")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    call cable_pack(dat,ssnow%snage(:),nmp(:,n))
    write(vname,'("t",I1.1,"_rtsoil")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    call cable_pack(dat,ssnow%rtsoil(:),nmp(:,n))
    write(vname,'("t",I1.1,"_GWwb")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    call cable_pack(dat,ssnow%GWwb(:),nmp(:,n))
    write(vname,'("t",I1.1,"_wtd")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    call cable_pack(dat,ssnow%wtd(:),nmp(:,n))
    write(vname,'("t",I1.1,"_cansto")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    call cable_pack(dat,canopy%cansto(:),nmp(:,n))
    write(vname,'("t",I1.1,"_us")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    call cable_pack(dat,canopy%us(:),nmp(:,n))
    write(vname,'("t",I1.1,"_pudsto")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    call cable_pack(dat,ssnow%pudsto(:),nmp(:,n))
    write(vname,'("t",I1.1,"_wetfac")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    call cable_pack(dat,ssnow%wetfac(:),nmp(:,n))
    write(vname,'("t",I1.1,"_ga")') n
    call histrd(iarchi-1,ierr,vname,dat,ifull)
    call cable_pack(dat,canopy%ga(:),nmp(:,n))
  end do
  
  ! soil temperature check
  if ( mp_global>0 ) then
    if ( any(ssnow%tgg>400.) ) then
      write(6,*) "ERROR: Invalid CABLE temperature when reading tile"
      write(6,*) "ssnow%tgg ",maxval(ssnow%tgg)
      stop -1
    end if
  end if 

end if ! ierr/=0 ..else..
  
if ( soil_struc==1 ) then
  if ( ierr_sli/=0 ) then
    if ( myid==0 ) write(6,*) "-> Use gridbox averaged data to initialise SLI"
  else 
    if ( myid==0 ) write(6,*) "-> Use tiled data to initialise SLI"  
    do n = 1,maxtile
      write(vname,'("t",I1.1,"_hzero")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,ssnow%h0(:),nmp(:,n))
      write(vname,'("t",I1.1,"_s")') n
      call histrd(iarchi-1,ierr,vname,datms(:,1:ms),ifull)
      do k = 1,ms
        call cable_pack(datms(:,k),ssnow%S(:,k),nmp(:,n))
      end do
      write(vname,'("t",I1.1,"_tsoil")') n
      call histrd(iarchi-1,ierr,vname,datms(:,1:ms),ifull)
      do k = 1,ms
        call cable_pack(datms(:,k),ssnow%tsoil(:,k),nmp(:,n))
      end do
      write(vname,'("t",I1.1,"_thetai")') n
      call histrd(iarchi-1,ierr,vname,datms(:,1:ms),ifull)
      do k = 1,ms
        call cable_pack(datms(:,k),ssnow%thetai(:,k),nmp(:,n))
      end do
      write(vname,'("t",I1.1,"_snowliq",I1.1)') n,1
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,ssnow%snowliq(:,1),nmp(:,n)) ! currently nsnow_max=1
      write(vname,'("t",I1.1,"_tsurface")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,ssnow%tsurface(:),nmp(:,n))
      write(vname,'("t",I1.1,"_nsnow")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      dati = nint(dat)
      call cable_pack(dati,ssnow%nsnow(:),nmp(:,n))
      write(vname,'("t",I1.1,"_fwsoil")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,canopy%fwsoil(:),nmp(:,n))
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
        call loadtile_carbonpools(carb_plant(:,:,1,m),casapool%cplant(:,:),m)
        call loadtile_carbonpools(carb_plant(:,:,2,m),casapool%nplant(:,:),m)
        call loadtile_carbonpools(carb_plant(:,:,3,m),casapool%pplant(:,:),m)
        call loadtile_carbonpools(carb_litter(:,:,1,m),casapool%clitter(:,:),m)
        call loadtile_carbonpools(carb_litter(:,:,2,m),casapool%nlitter(:,:),m)
        call loadtile_carbonpools(carb_litter(:,:,3,m),casapool%plitter(:,:),m)
        call loadtile_carbonpools(carb_soil(:,:,1,m),casapool%csoil(:,:),m)
        call loadtile_carbonpools(carb_soil(:,:,2,m),casapool%nsoil(:,:),m)
        call loadtile_carbonpools(carb_soil(:,:,3,m),casapool%psoil(:,:),m)
      end do   ! mm = 1,10
      m = 14 ! use index 11 to store pft 14
      call loadtile_carbonpools(carb_plant(:,:,1,11),casapool%cplant(:,:),m)
      call loadtile_carbonpools(carb_plant(:,:,2,11),casapool%nplant(:,:),m)
      call loadtile_carbonpools(carb_plant(:,:,3,11),casapool%pplant(:,:),m)
      call loadtile_carbonpools(carb_litter(:,:,1,11),casapool%clitter(:,:),m)
      call loadtile_carbonpools(carb_litter(:,:,2,11),casapool%nlitter(:,:),m)
      call loadtile_carbonpools(carb_litter(:,:,3,11),casapool%plitter(:,:),m)
      call loadtile_carbonpools(carb_soil(:,:,1,11),casapool%csoil(:,:),m)
      call loadtile_carbonpools(carb_soil(:,:,2,11),casapool%nsoil(:,:),m)
      call loadtile_carbonpools(carb_soil(:,:,3,11),casapool%psoil(:,:),m)
      deallocate( carb_plant, carb_litter, carb_soil )
    end if  
  else
    if ( myid==0 ) write(6,*) "-> Use tiled data to initialise CASA-CNP"  
    do n = 1,maxtile
      write(vname,'("t",I1.1,"_cplant")') n
      call histrd(iarchi-1,ierr,vname,datmplant(:,1:mplant),ifull)
      do k = 1,mplant
        call cable_pack(datmplant(:,k),casapool%cplant(:,k),nmp(:,n))
      end do
      write(vname,'("t",I1.1,"_nplant")') n
      call histrd(iarchi-1,ierr,vname,datmplant(:,1:mplant),ifull)
      do k = 1,mplant
        call cable_pack(datmplant(:,k),casapool%nplant(:,k),nmp(:,n))
      end do
      write(vname,'("t",I1.1,"_pplant")') n
      call histrd(iarchi-1,ierr,vname,datmplant(:,1:mplant),ifull)
      do k = 1,mplant
        call cable_pack(datmplant(:,k),casapool%pplant(:,k),nmp(:,n))
      end do
      write(vname,'("t",I1.1,"_clitter")') n
      call histrd(iarchi-1,ierr,vname,datmlitter(:,1:mlitter),ifull)
      do k = 1,mlitter
        call cable_pack(datmlitter(:,k),casapool%clitter(:,k),nmp(:,n))
      end do
      write(vname,'("t",I1.1,"_nlitter")') n
      call histrd(iarchi-1,ierr,vname,datmlitter(:,1:mlitter),ifull)
      do k = 1,mlitter
        call cable_pack(datmlitter(:,k),casapool%nlitter(:,k),nmp(:,n))
      end do
      write(vname,'("t",I1.1,"_plitter")') n
      call histrd(iarchi-1,ierr,vname,datmlitter(:,1:mlitter),ifull)
      do k = 1,mlitter
        call cable_pack(datmlitter(:,k),casapool%plitter(:,k),nmp(:,n))
      end do
      write(vname,'("t",I1.1,"_csoil")') n
      call histrd(iarchi-1,ierr,vname,datmsoil(:,1:msoil),ifull)
      do k = 1,msoil
        call cable_pack(datmsoil(:,k),casapool%csoil(:,k),nmp(:,n))
      end do
      write(vname,'("t",I1.1,"_nsoil")') n
      call histrd(iarchi-1,ierr,vname,datmsoil(:,1:msoil),ifull)
      do k = 1,msoil
        call cable_pack(datmsoil(:,k),casapool%nsoil(:,k),nmp(:,n))
      end do
      write(vname,'("t",I1.1,"_psoil")') n
      call histrd(iarchi-1,ierr,vname,datmsoil(:,1:msoil),ifull)
      do k = 1,msoil
        call cable_pack(datmsoil(:,k),casapool%psoil(:,k),nmp(:,n))
      end do
      write(vname,'("t",I1.1,"_glai")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,casamet%glai,nmp(:,n))
      write(vname,'("t",I1.1,"_phen")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,phen%phen,nmp(:,n))
      write(vname,'("t",I1.1,"_aphen")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,phen%aphen,nmp(:,n))
      write(vname,'("t",I1.1,"_phenphase")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      dati = nint(dat)
      call cable_pack(dati,phen%phase,nmp(:,n))
      write(vname,'("t",I1.1,"_doyphase3")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      dati = nint(dat)
      call cable_pack(dati,phen%doyphase(:,3),nmp(:,n))
      write(vname,'("t",I1.1,"_clabile")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,casapool%clabile,nmp(:,n))
      write(vname,'("t",I1.1,"_nsoilmin")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,casapool%nsoilmin,nmp(:,n))
      write(vname,'("t",I1.1,"_psoillab")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,casapool%psoillab,nmp(:,n))
      write(vname,'("t",I1.1,"_psoilsorb")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,casapool%psoilsorb,nmp(:,n))
      write(vname,'("t",I1.1,"_psoilocc")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,casapool%psoilocc,nmp(:,n))
      write(vname,'("t",I1.1,"_crmplant")') n
      call histrd(iarchi-1,ierr,vname,datmplant(:,1:mplant),ifull)
      do k = 1,mplant
        call cable_pack(datmplant(:,k),casaflux%crmplant(:,k),nmp(:,n))
      end do
      write(vname,'("t",I1.1,"_fracsapwood")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,casaflux%frac_sapwood,nmp(:,n))
      write(vname,'("t",I1.1,"_sapwoodarea")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,casaflux%sapwood_area,nmp(:,n))
      write(vname,'("t",I1.1,"_crsoil")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,casaflux%crsoil,nmp(:,n))
      write(vname,'("t",I1.1,"_cnpp")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,casaflux%cnpp,nmp(:,n))
      write(vname,'("t",I1.1,"_clabloss")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,casaflux%clabloss,nmp(:,n))
      write(vname,'("t",I1.1,"_crgplant")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,casaflux%crgplant,nmp(:,n))
      write(vname,'("t",I1.1,"_stemnpp")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,casaflux%stemnpp,nmp(:,n))
      write(vname,'("t",I1.1,"_LAImax")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,casabal%laimax,nmp(:,n))
      write(vname,'("t",I1.1,"_Cleafmean")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,casabal%cleafmean,nmp(:,n))
      write(vname,'("t",I1.1,"_Crootmean")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,casabal%crootmean,nmp(:,n))
      write(vname,'("t",I1.1,"_fpn")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,canopy%fpn,nmp(:,n))
      write(vname,'("t",I1.1,"_frday")') n
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call cable_pack(dat,canopy%frday,nmp(:,n))
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
    np_pop = size(pop%pop_grid)
    allocate( dat_out(np_pop), dati_out(np_pop) )
    dat_out = 0._8
    dati_out = 0
    do n = 1,maxtile  
      write(vname,'("t",I1.1,"_pop_grid_cmass_sum")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%cmass_sum = dat_out(1:np_pop)
      write(vname,'("t",I1.1,"_pop_grid_cmass_sum_old")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%cmass_sum_old = dat_out(1:np_pop)
      write(vname,'("t",I1.1,"_pop_grid_cheartwood_sum")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%cheartwood_sum = dat_out(1:np_pop)
      write(vname,'("t",I1.1,"_pop_grid_csapwood_sum")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%csapwood_sum = dat_out(1:np_pop)
      write(vname,'("t",I1.1,"_pop_grid_csapwood_sum_old")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%csapwood_sum_old = dat_out(1:np_pop)
      write(vname,'("t",I1.1,"_pop_grid_densindiv")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%densindiv = dat_out(1:np_pop)
      write(vname,'("t",I1.1,"_pop_grid_height_mean")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%height_mean = dat_out(1:np_pop)
      write(vname,'("t",I1.1,"_pop_grid_height_max")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%height_max = dat_out(1:np_pop) 
      write(vname,'("t",I1.1,"_pop_grid_basal_area")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%basal_area = dat_out(1:np_pop) 
      write(vname,'("t",I1.1,"_pop_grid_sapwood_loss")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%sapwood_loss = dat_out(1:np_pop)
      write(vname,'("t",I1.1,"_pop_grid_sapwood_area_loss")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%sapwood_area_loss = dat_out(1:np_pop)
      write(vname,'("t",I1.1,"_pop_grid_stress_mortality")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%stress_mortality = dat_out(1:np_pop)
      write(vname,'("t",I1.1,"_pop_grid_crowding_mortality")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%crowding_mortality = dat_out(1:np_pop)
      write(vname,'("t",I1.1,"_pop_grid_fire_mortality")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%fire_mortality = dat_out(1:np_pop)
      write(vname,'("t",I1.1,"_pop_grid_cat_mortality")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%cat_mortality = dat_out(1:np_pop)
      write(vname,'("t",I1.1,"_pop_grid_res_mortality")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%res_mortality = dat_out(1:np_pop)
      write(vname,'("t",I1.1,"_pop_grid_growth")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%growth = dat_out(1:np_pop) 
      write(vname,'("t",I1.1,"_pop_grid_area_growth")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%area_growth = dat_out(1:np_pop) 
      write(vname,'("t",I1.1,"_pop_grid_crown_cover")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%crown_cover = dat_out(1:np_pop)
      write(vname,'("t",I1.1,"_pop_grid_crown_area")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%crown_area = dat_out(1:np_pop)
      write(vname,'("t",I1.1,"_pop_grid_crown_volume")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%crown_volume = dat_out(1:np_pop)
      write(vname,'("t",I1.1,"_pop_grid_sapwood_area")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%sapwood_area = dat_out(1:np_pop) 
      write(vname,'("t",I1.1,"_pop_grid_sapwood_area_old")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%sapwood_area_old = dat_out(1:np_pop)
      write(vname,'("t",I1.1,"_pop_grid_KClump")') n  
      call histrd(iarchi-1,ierr,vname,dat,ifull)
      call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
      pop%pop_grid(:)%KClump = dat_out(1:np_pop)
      write(vname,'("t",I1.1,"_pop_grid_freq_age")') n
      call histrd(iarchi-1,ierr,vname,datage,ifull)
      do k = 1,POP_AGEMAX
        call pop_pack(datage(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%freq_age(k) = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_biomass_age")') n
      call histrd(iarchi-1,ierr,vname,datage,ifull)
      do k = 1,POP_AGEMAX
        call pop_pack(datage(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%biomass_age(k) = dat_out(1:np_pop)
      end do
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_biomass",I1.1)') n,ll  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%biomass(ll) = dat_out(1:np_pop)
      end do
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_density",I1.1)') n,ll
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%density(ll) = dat_out(1:np_pop)
      end do  
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_hmean",I1.1)') n,ll  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%hmean(ll) = dat_out(1:np_pop)
      end do  
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_hmax",I1.1)') n,ll
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%hmax(ll) = dat_out(1:np_pop)
      end do
      do hh = 1,POP_HEIGHT_BINS
        write(vname,'("t",I1.1,"_pop_grid_cmass_stem_bin",I2.2)') n,hh  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        call pop_pack(dat,dat_out(1:np_pop),nmp(:,n)) 
        pop%pop_grid(:)%cmass_stem_bin(hh) = dat_out(1:np_pop)
      end do  
      do hh = 1,POP_HEIGHT_BINS
        write(vname,'("t",I1.1,"_pop_grid_densindiv_bin",I2.2)') n,hh
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%densindiv_bin(hh) = dat_out(1:np_pop)
      end do
      do hh = 1,POP_HEIGHT_BINS
        write(vname,'("t",I1.1,"_pop_grid_height_bin",I2.2)') n,hh  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%height_bin(hh) = dat_out(1:np_pop)
      end do
      do hh = 1,POP_HEIGHT_BINS
        write(vname,'("t",I1.1,"_pop_grid_diameter_bin",I2.2)') n,hh
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        call pop_pack(dat,dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%diameter_bin(hh) = dat_out(1:np_pop)
      end do
      do dd = 1,POP_NDISTURB
        write(vname,'("t",I1.1,"_pop_grid_n_age",I1.1)') n,dd  
        call histrd(iarchi-1,ierr,vname,dat,ifull)
        dati = nint(dat)  
        call pop_pack(dati,dati_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%n_age(dd) = dati_out(1:np_pop)
      end do  
      write(vname,'("t",I1.1,"_pop_grid_patch_id")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        dati = nint(datpatch(:,k))  
        call pop_pack(dati,dati_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%id = dati_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_freq")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%freq(k) = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_freq_old")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%freq_old(k) = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_factor_recruit")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%factor_recruit = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_pgap")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%pgap = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_lai")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%lai = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_biomass")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%biomass = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_biomass_old")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%biomass_old = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%sapwood = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_heartwood")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%heartwood = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_old")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%sapwood_old = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_area")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%sapwood_area = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_area_old")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%sapwood_area_old = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_stress_mortality")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%stress_mortality = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_fire_mortality")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%fire_mortality = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_cat_mortality")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%cat_mortality = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_crowding_mortality")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%crowding_mortality = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_cpc")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%cpc = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_mortality")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%mortality = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_loss")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%sapwood_loss = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_area_loss")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%sapwood_area_loss = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_growth")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%growth = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_area_growth")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%area_growth = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_frac_NPP")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%frac_NPP = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_frac_respiration")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%frac_respiration = dat_out(1:np_pop)
      end do
      write(vname,'("t",I1.1,"_pop_grid_patch_frac_light_uptake")') n
      call histrd(iarchi-1,ierr,vname,datpatch,ifull)
      do k = 1,POP_NPATCH
        call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
        pop%pop_grid(:)%patch(k)%frac_light_uptake = dat_out(1:np_pop)
      end do
      do dd = 1,POP_NDISTURB
        write(vname,'("t",I1.1,"_pop_grid_patch_disturbance_interval",I1.1)') n,dd
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        do k = 1,POP_NPATCH
          dati = nint(datpatch(:,k))  
          call pop_pack(dati,dati_out(1:np_pop),nmp(:,n))
          pop%pop_grid(:)%patch(k)%disturbance_interval(dd) = dati_out(1:np_pop)
        end do
      end do  
      do dd = 1,POP_NDISTURB
        write(vname,'("t",I1.1,"_pop_grid_patch_first_disturbance_year",I1.1)') n,dd
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        do k = 1,POP_NPATCH
          dati = nint(datpatch(:,k))   
          call pop_pack(dati,dati_out(1:np_pop),nmp(:,n))
          pop%pop_grid(:)%patch(k)%first_disturbance_year(dd) = dati_out(1:np_pop)
        end do
      end do  
      do dd = 1,POP_NDISTURB
        write(vname,'("t",I1.1,"_pop_grid_patch_age",I1.1)') n,dd
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        do k = 1,POP_NPATCH
          dati = nint(datpatch(:,k)) 
          call pop_pack(dati,dati_out(1:np_pop),nmp(:,n))
          pop%pop_grid(:)%patch(k)%age(dd) = dati_out(1:np_pop)
        end do
      end do  
      do dd = 1,POP_NDISTURB
        write(vname,'("t",I1.1,"_pop_grid_ranked_age_unique",I1.1)') n,dd
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        do k = 1,POP_NPATCH
          dati = nint(datpatch(:,k))  
          call pop_pack(dati,dati_out(1:np_pop),nmp(:,n))
          pop%pop_grid(:)%ranked_age_unique(k,dd) = dati_out(1:np_pop)
        end do
      end do  
      do dd = 1,POP_NDISTURB
        write(vname,'("t",I1.1,"_pop_grid_freq_ranked_age_unique",I1.1)') n,dd
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        do k = 1,POP_NPATCH
          call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
          pop%pop_grid(:)%freq_ranked_age_unique(k,dd) = dat_out(1:np_pop)
        end do
      end do  
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_ncohort")') n,ll
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        do k = 1,POP_NPATCH  
          dati = nint(datpatch(:,k))  
          call pop_pack(dati,dati_out(1:np_pop),nmp(:,n))
          pop%pop_grid(:)%patch(k)%layer(ll)%ncohort = dati_out(1:np_pop)
        end do
      end do  
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_biomass")') n,ll
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        do k = 1,POP_NPATCH  
          call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
          pop%pop_grid(:)%patch(k)%layer(ll)%biomass = dat_out(1:np_pop)
        end do
      end do  
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_density")') n,ll
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        do k = 1,POP_NPATCH  
          call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
          pop%pop_grid(:)%patch(k)%layer(ll)%density = dat_out(1:np_pop)
        end do
      end do  
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_hmean")') n,ll
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        do k = 1,POP_NPATCH  
          call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
          pop%pop_grid(:)%patch(k)%layer(ll)%hmean = dat_out(1:np_pop)
        end do
      end do  
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_hmax")') n,ll
        call histrd(iarchi-1,ierr,vname,datpatch,ifull)
        do k = 1,POP_NPATCH  
          call pop_pack(datpatch(:,k),dat_out(1:np_pop),nmp(:,n))
          pop%pop_grid(:)%patch(k)%layer(ll)%hmax = dat_out(1:np_pop)
        end do
      end do  
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_age")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT                  
          do k = 1,POP_NPATCH
            dati = nint(datpc(:,k,cc))  
            call pop_pack(dati,dati_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%age = dati_out(1:np_pop)
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_id")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH
            dati = nint(datpc(:,k,cc))  
            call pop_pack(dati,dati_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%id = dati_out(1:np_pop)
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_biomass")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH  
            call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%biomass = dat_out(1:np_pop)
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_density")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH  
            call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%density = dat_out(1:np_pop)
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_resource_uptake")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH  
            call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_resource_uptake = dat_out(1:np_pop)
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_light_uptake")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT    
          do k = 1,POP_NPATCH  
            call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_light_uptake = dat_out(1:np_pop)
          end do  
        end do
      end do        
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_interception")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_interception = dat_out(1:np_pop)
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_respiration")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_respiration = dat_out(1:np_pop)
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_NPP")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_NPP = dat_out(1:np_pop)
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_respiration_scalar")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%respiration_scalar = dat_out(1:np_pop)
          end do  
        end do
      end do         
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_crown_area")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%crown_area = dat_out(1:np_pop)
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Pgap")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%Pgap = dat_out(1:np_pop)
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_height")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%height = dat_out(1:np_pop)
          end do  
        end do
      end do
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_diameter")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%diameter = dat_out(1:np_pop)
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_sapwood")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%sapwood = dat_out(1:np_pop)
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_heartwood")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%heartwood = dat_out(1:np_pop)
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_sapwood_area")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%sapwood_area = dat_out(1:np_pop)
          end do  
        end do
      end do         
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_basal_area")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%basal_area = dat_out(1:np_pop)
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_LAI")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%LAI = dat_out(1:np_pop)
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Cleaf")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%Cleaf = dat_out(1:np_pop)
          end do  
        end do
      end do 
      do ll = 1,POP_NLAYER
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Croot")') n,ll
        call histrd(iarchi-1,ierr,vname,datpc,ifull)
        do cc = 1,POP_NCOHORT
          do k = 1,POP_NPATCH  
            call pop_pack(datpc(:,k,cc),dat_out(1:np_pop),nmp(:,n))
            pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%Croot = dat_out(1:np_pop)
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
  
call redistribute_tile(old_sv)
call fixtile

! Calculate LAI and veg fraction diagnostics
vlai(:) = 0.
sigmf(:) = 0.
if ( mp_global>0 ) then
  call getzinp(jyear,jmonth,jday,jhour,jmin,mins)
  call setlai(sigmf,jmonth,jday,jhour,jmin,mp_global,sv,vl2,casamet,veg,ifull)
  dummy_unpack = sv*real(veg%vlai)
  call cable_unpack(dummy_unpack,vlai)
end if

return
end subroutine loadtile

subroutine loadtile_carbonpools(dat_in,dat_out,mm)

use newmpar_m

implicit none

integer, intent(in) :: mm
integer nb, k, msize
real, dimension(:,:), intent(in) :: dat_in
real(kind=8), dimension(:,:), intent(inout) :: dat_out
real(kind=8), dimension(size(dat_out,1)) :: dummy_unpack

msize = size(dat_in,2)

do k = 1,msize
  do nb = 1,maxnb
    dummy_unpack(:) = 0._8  
    call cable_pack(dat_in(:,k),dummy_unpack(:),nb)
  end do
  where ( veg%iveg(:)==mm .and. dummy_unpack(:)>1.e-8_8 )
    dat_out(:,k) = dummy_unpack(:)
  end where
end do

return
end subroutine loadtile_carbonpools

subroutine fixtile

use carbpools_m
use cc_mpi
use darcdf_m
use infile
use newmpar_m
use parm_m
use soil_m
use soilsnow_m
use vegpar_m
  
integer k
real totdepth
 
! Some fixes for rounding errors
if ( mp_global>0 ) then

  totdepth = 0.
  do k = 1,ms
    totdepth = totdepth + real(soil%zse(k))*100.
  enddo
  
  ssnow%tgg = max(ssnow%tgg, 200._8)
  ssnow%tggsn = max(ssnow%tggsn, 200._8)

  ssnow%wb = max(ssnow%wb,0._8)
  ssnow%wbice = min( max(ssnow%wbice, 0._8), ssnow%wb )
  ssnow%smass = max(ssnow%smass,0._8)
  ssnow%rtsoil = max(ssnow%rtsoil,0._8)
  ssnow%snowd = max(ssnow%snowd,0._8)
  ssnow%osnowd = max(ssnow%osnowd,0._8)
  ssnow%wetfac = min(max(ssnow%wetfac,0._8),1._8)
  canopy%cansto = max(canopy%cansto,0._8)
  ssnow%wbliq = ssnow%wb - ssnow%wbice

  ssnow%wbtot = 0._8
  ssnow%wbtot1 = 0._8
  ssnow%wbtot2 = 0._8
  ssnow%tggav = 0._8
  do k = 1,ms
    ssnow%wbtot = ssnow%wbtot+ssnow%wb(:,k)*1000._8*soil%zse(k)
    ssnow%tggav = ssnow%tggav+soil%zse(k)*ssnow%tgg(:,k)/real(totdepth/100.,8)
    ssnow%gammzz(:,k) = max((1._8-soil%ssat)*soil%css* soil%rhosoil                         &
        + real(ssnow%wb(:,k)-ssnow%wbice(:,k))*4.218e3_8* 1000._8                           &
        + real(ssnow%wbice(:,k))*2.100e3_8*1000._8*0.9_8,soil%css*soil%rhosoil)*soil%zse(k) &
        + (1.-ssnow%isflag)*2090._8*ssnow%snowd
  end do

  if ( ccycle > 0 ) then
    casapool%cplant     = max(0._8,casapool%cplant)
    casapool%clitter    = max(0._8,casapool%clitter)
    casapool%csoil      = max(0._8,casapool%csoil)
    casabal%cplantlast  = casapool%cplant
    casabal%clitterlast = casapool%clitter
    casabal%csoillast   = casapool%csoil
    casabal%clabilelast = casapool%clabile
    casabal%sumcbal      = 0._8
    casabal%FCgppyear    = 0._8
    casabal%FCrpyear     = 0
    casabal%FCrmleafyear = 0._8
    casabal%FCrmwoodyear = 0._8
    casabal%FCrmrootyear = 0._8
    casabal%FCrgrowyear  = 0._8
    casabal%FCnppyear    = 0._8
    casabal%FCrsyear     = 0._8
    casabal%FCneeyear    = 0._8
    casabal%dCdtyear     = 0._8
    casapool%nplant      = max(1.e-6_8,casapool%nplant)
    casapool%nlitter     = max(1.e-6_8,casapool%nlitter)
    casapool%nsoil       = max(1.e-6_8,casapool%nsoil)
    casapool%nsoilmin    = max(1.e-6_8,casapool%nsoilmin)
    casabal%nplantlast   = casapool%nplant
    casabal%nlitterlast  = casapool%nlitter
    casabal%nsoillast    = casapool%nsoil
    casabal%nsoilminlast = casapool%nsoilmin
    casabal%sumnbal      = 0._8
    casabal%FNdepyear    = 0._8
    casabal%FNfixyear    = 0._8
    casabal%FNsnetyear   = 0._8
    casabal%FNupyear     = 0._8
    casabal%FNleachyear  = 0._8
    casabal%FNlossyear   = 0._8
    casapool%pplant      = max(1.0e-7_8,casapool%pplant)
    casapool%plitter     = max(1.0e-7_8,casapool%plitter)
    casapool%psoil       = max(1.0e-7_8,casapool%psoil)
    casapool%Psoillab    = max(1.0e-7_8,casapool%psoillab)
    casapool%psoilsorb   = max(1.0e-7_8,casapool%psoilsorb)
    casapool%psoilocc    = max(1.0e-7_8,casapool%psoilocc)
    casabal%pplantlast   = casapool%pplant
    casabal%plitterlast  = casapool%plitter
    casabal%psoillast    = casapool%psoil
    casabal%psoillablast = casapool%psoillab
    casabal%psoilsorblast= casapool%psoilsorb
    casabal%psoilocclast = casapool%psoilocc
    casabal%sumpbal      = 0._8
    casabal%FPweayear    = 0._8
    casabal%FPdustyear   = 0._8
    casabal%FPsnetyear   = 0._8
    casabal%FPupyear     = 0._8
    casabal%FPleachyear  = 0._8
    casabal%FPlossyear   = 0._8
    !casamet%glai         = max(min( casamet%glai, casabiome%glaimax(veg%iveg)), casabiome%glaimin(veg%iveg))
    
    where ( .not.( casamet%iveg2==forest.or.casamet%iveg2==shrub ) )
      casapool%cplant(:,wood)  = 0._8
      casapool%clitter(:,cwd)  = 0._8
      casapool%nplant(:,wood)  = 0._8
      casapool%nlitter(:,cwd)  = 0._8
      casapool%pplant(:,wood)  = 0._8
      casapool%plitter(:,cwd)  = 0._8
    end where

    ! initializing glai in case not reading pool file (eg. during spin)
    casapool%ratioNClitter = casapool%nlitter/(casapool%clitter+1.0e-10_8)
    casapool%ratioNPlitter = casapool%nlitter/(casapool%plitter+1.0e-10_8)
    casapool%ratioPClitter = casapool%plitter/(casapool%clitter+1.0e-10_8)

    if ( ccycle<2 ) then
      casapool%Nplant = casapool%Cplant*casapool%ratioNCplant
      casapool%Nsoil  = casapool%ratioNCsoil*casapool%Csoil
    else if ( ccycle<3 ) then
      casapool%Psoil  = casapool%Nsoil/casapool%ratioNPsoil
    end if
    
    where ( veg%iveg==6 .or. veg%iveg==7 .or. veg%iveg==8 .or. veg%iveg==9 .or. veg%iveg==10 )
      casapool%cplant(:,wood) = 0._8
      casapool%clitter(:,cwd) = 0._8
      casapool%nplant(:,wood) = 0._8
      casapool%nlitter(:,cwd) = 0._8  
      casapool%pplant(:,wood) = 0._8
      casapool%plitter(:,cwd) = 0._8   
    end where
    where ( veg%iveg==11 .or. veg%iveg==12 .or. veg%iveg==13 .or. veg%iveg==15 .or. veg%iveg==16 .or. &
            veg%iveg==17 )
      casapool%cplant(:,leaf) = 0._8
      casapool%cplant(:,wood) = 0._8
      casapool%cplant(:,froot) = 0._8
      casapool%clitter(:,metb) = 0._8
      casapool%clitter(:,str) = 0._8
      casapool%clitter(:,cwd) = 0._8
      casapool%csoil(:,mic) = 0._8
      casapool%csoil(:,slow) = 0._8
      casapool%csoil(:,pass) = 0._8
      casapool%nplant(:,leaf) = 0._8
      casapool%nplant(:,wood) = 0._8
      casapool%nplant(:,froot) = 0._8
      casapool%nlitter(:,metb) = 0._8
      casapool%nlitter(:,str) = 0._8
      casapool%nlitter(:,cwd) = 0._8
      casapool%nsoil(:,mic) = 0._8
      casapool%nsoil(:,slow) = 0._8
      casapool%nsoil(:,pass) = 0._8
      casapool%pplant(:,leaf) = 0._8
      casapool%pplant(:,wood) = 0._8
      casapool%pplant(:,froot) = 0._8
      casapool%plitter(:,metb) = 0._8
      casapool%plitter(:,str) = 0._8
      casapool%plitter(:,cwd) = 0._8
      casapool%psoil(:,mic) = 0._8
      casapool%psoil(:,slow) = 0._8
      casapool%psoil(:,pass) = 0._8
    end where
    where ( veg%iveg==14 )        
      casapool%cplant(:,leaf) = 0._8
      casapool%cplant(:,wood) = 0._8
      casapool%cplant(:,froot) = 0._8
      casapool%pplant(:,leaf) = 0._8
      casapool%pplant(:,wood) = 0._8
      casapool%pplant(:,froot) = 0._8
    end where
    
  end if ! ccycle>0
end if   ! mp_global>0
  
return
end subroutine fixtile

! if the vegetation class has changed for a tile, then attempt to find the tile associated
! with the equivilent vegetation class in the restart data
subroutine create_new_tile_map(old_cv,nmp)

use cc_mpi
use newmpar_m, only : ifull

integer n, m, iq, store_n
integer, dimension(mp_global), intent(in) :: old_cv
integer, dimension(ifull,maxtile), intent(inout) :: nmp
integer, dimension(ifull,maxtile) :: oldv_up, newv_up
integer, dimension(maxtile) :: oldv_v, newv_v
logical :: cv_test

if ( myid==0 ) write(6,*) "-> Create map to vegetation types"
if ( mp_global>0 ) then
  cv_test = any( cveg/=old_cv )
else
  cv_test = .false.
end if

oldv_up = 0
newv_up = 0

! default map  
do iq = 1,ifull
  do n = 1,maxtile  
    nmp(iq,n) = n  
  end do  
end do  

if ( mp_global>0 ) then
  if ( cv_test ) then
    
    do n = 1,maxtile
      call cable_unpack(cveg,newv_up(:,n),n)
      call cable_unpack(old_cv,oldv_up(:,n),n)
    end do     
    
    do iq = 1,ifull
      newv_v = newv_up(iq,:) ! current vegetation class
      oldv_v = oldv_up(iq,:) ! vegetation class in restart file
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
    
  end if
end if

return
end subroutine create_new_tile_map


! redistribute temperature and water with a gridbox if tile area fraction changes
subroutine redistribute_tile(old_sv)

use cc_mpi

integer k
real, dimension(mp_global), intent(in) :: old_sv
logical :: sv_test

! check if any point have land-cover change
if ( mp_global>0 ) then
  sv_test = any( abs(sv-old_sv)>1.e-8 )
else
  sv_test = .false. 
end if

if ( mp_global>0 ) then
  if ( sv_test ) then

    ! check for errors prior to redistribution
    if ( any(ssnow%tgg>400.) ) then
      write(6,*) "ERROR: Invalid input temperature for CABLE redistribute_tile"
      write(6,*) "ssnow%tgg ",maxval(ssnow%tgg)
      call ccmpi_abort(-1)
    end if     
    
    ! assume common soil texture and soil heat capacity
    do k = 1,ms
      call redistribute_work(old_sv,ssnow%tgg(:,k))
      call redistribute_work(old_sv,ssnow%wb(:,k))
      call redistribute_work(old_sv,ssnow%wbice(:,k))
    end do
    call redistribute_work(old_sv,ssnow%GWwb)
    if ( soil_struc==1 ) then
      do k = 1,ms
        call redistribute_work(old_sv,ssnow%tsoil(:,k))
      end do
    end if
    
    ! Do we need a special treatment for snow?
    do k = 1,3
      call redistribute_work(old_sv,ssnow%tggsn(:,k))
      call redistribute_work(old_sv,ssnow%smass(:,k))
      call redistribute_work(old_sv,ssnow%ssdn(:,k))
      call redistribute_work(old_sv,ssnow%sdepth(:,k))
      call redistribute_work(old_sv,ssnow%sconds(:,k))
    end do
    call redistribute_work(old_sv,ssnow%ssdnn(:))
    call redistribute_work(old_sv,ssnow%snowd(:))
    call redistribute_work(old_sv,ssnow%osnowd(:))
    call redistribute_work(old_sv,ssnow%snage(:))
    call redistribute_work(old_sv,ssnow%rtsoil(:))
    call redistribute_work(old_sv,ssnow%GWwb(:))
    call redistribute_work(old_sv,ssnow%wtd(:))
    call redistribute_work(old_sv,canopy%cansto(:))
    call redistribute_work(old_sv,ssnow%pudsto(:))
    call redistribute_work(old_sv,ssnow%wetfac(:))
    
    ! also treatment for carbon
    if ( ccycle/=0 ) then
      do k = 1,3    
        call redistribute_work(old_sv,casapool%cplant(:,k))  
        call redistribute_work(old_sv,casapool%clitter(:,k))
        call redistribute_work(old_sv,casapool%csoil(:,k))
        call redistribute_work(old_sv,casapool%nplant(:,k))
        call redistribute_work(old_sv,casapool%nlitter(:,k))
        call redistribute_work(old_sv,casapool%nsoil(:,k))
        call redistribute_work(old_sv,casapool%pplant(:,k))
        call redistribute_work(old_sv,casapool%plitter(:,k))
        call redistribute_work(old_sv,casapool%psoil(:,k))
      end do    
      !call redistribute_work(old_sv,casamet%glai)
      !call redistribute_work(old_sv,phen%phen)
      !call redistribute_work(old_sv,phen%aphen)
    end if

    ! check for errors after redistribution
    if ( any(ssnow%tgg>400.) ) then
      write(6,*) "ERROR: Invalid output temperature for CABLE redistribute_tile"
      write(6,*) "ssnow%tgg ",maxval(ssnow%tgg)
      call ccmpi_abort(-1)
    end if  
    
  end if
end if
  
return
end subroutine redistribute_tile

subroutine redistribute_work(old_sv,vdata)

use newmpar_m
use soil_m

integer tile, nb, iq, is, ie
real, dimension(mp_global), intent(in) :: old_sv
real, dimension(ifull,maxnb) :: up_new_svs, up_old_svs
real, dimension(ifull) :: svs_sum
real, dimension(maxnb) :: adj_pos_frac, adj_neg_frac, new_svs, old_svs
real(kind=8) :: ave_neg_vdata, adj_neg_frac_sum
real(kind=8), dimension(mp_global), intent(inout) :: vdata
real(kind=8), dimension(ifull,maxnb) :: up_vdata
real(kind=8), dimension(maxnb) :: old_vdata

! update data
up_vdata = 0._8
up_new_svs = 0.
up_old_svs = 0.
do nb = 1,maxnb
  call cable_unpack(sv,up_new_svs(:,nb),nb)
  call cable_unpack(old_sv,up_old_svs(:,nb),nb)
  call cable_unpack(vdata,up_vdata(:,nb),nb)
end do  

svs_sum = sum(up_new_svs,dim=2)
do nb = 1,maxnb
  where ( land(1:ifull) )  
    up_new_svs(:,nb) = up_new_svs(:,nb)/max(svs_sum(:),1.e-10)
  end where  
end do
svs_sum = sum(up_old_svs,dim=2)
do nb = 1,maxnb
  where ( land(1:ifull) )  
    up_old_svs(:,nb) = up_old_svs(:,nb)/max(svs_sum(:),1.e-10)
  end where  
end do

do tile = 1,ntiles
  is = (tile-1)*imax + 1
  ie = tile*imax
  do iq = is,ie
    new_svs(:) = up_new_svs(iq,:)
    old_svs(:) = up_old_svs(iq,:)      
    adj_pos_frac(:) = max( new_svs(:)-old_svs(:), 0. ) ! tiles with increasing fraction
    adj_neg_frac(:) = max( old_svs(:)-new_svs(:), 0. ) ! tiles with decreasing fraction
    adj_neg_frac_sum = sum( adj_neg_frac )   
    ! check adj_neg_frac_sum>0 which can occur if water has changed to land
    if ( land(iq) .and. any( abs(new_svs-old_svs)>1.e-8 ) .and. adj_neg_frac_sum>1.e-8 ) then 
      old_vdata(:) = up_vdata(iq,:)
      ! summarise tiles with decreasing area fraction with an average value
      ave_neg_vdata = sum( adj_neg_frac(:)*old_vdata(:) )/adj_neg_frac_sum
      ! Only change tiles that are increasing in area fraction
      do nb = 1,maxnb
        if ( adj_pos_frac(nb)>1.e-8 ) then  
          up_vdata(iq,nb) = (old_vdata(nb)*old_svs(nb) + ave_neg_vdata*adj_pos_frac(nb)) &
                        /(old_svs(nb) + adj_pos_frac(nb))
        end if  
      end do
    end if    
  end do  
end do

do nb = 1,maxnb
  call cable_pack(up_vdata(:,nb),vdata,nb)
end do  

return
end subroutine redistribute_work    
    
! *************************************************************************************
! This subroutine saves CABLE tile data
subroutine savetiledef(idnc,local,jdim,jsize,cdim,csize,itype)

use carbpools_m
use cc_mpi, only : myid
use infile
use newmpar_m
use parm_m, only : diaglevel_pop
  
integer, intent(in) :: idnc, jsize
integer k,n
integer ll,dd,hh
integer, dimension(2), intent(in) :: csize
integer, dimension(jsize), intent(in) :: jdim  
integer, dimension(6,7), intent(in) :: cdim
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
      do k = 1,ms
        write(lname,'("Soil temperature tile ",I1.1," lev ",I1.1)') n,k
        write(vname,'("t",I1.1,"_tgg",I1.1)') n,k
        call attrib(idnc,jdim,jsize,vname,lname,'K',100.,400.,any_m,point_m,land_m,double_m)
      end do
      do k = 1,ms
        write(lname,'("Soil moisture tile ",I1.1," lev ",I1.1)') n,k
        write(vname,'("t",I1.1,"_wb",I1.1)') n,k 
        call attrib(idnc,jdim,jsize,vname,lname,'m3/m3',0.,2.6,any_m,point_m,land_m,double_m)
      end do
      do k = 1,ms
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
        do k = 1,ms
          write(lname,'("S tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_s",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        end do
        do k = 1,ms
          write(lname,'("tsoil tile ",I1.1," lev ",I1.1)') n,k
          write(vname,'("t",I1.1,"_tsoil",I1.1)') n,k
          call attrib(idnc,jdim,jsize,vname,lname,'none',0.,65000.,any_m,point_m,land_m,double_m)
        end do
        do k = 1,ms
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
        call attrib(idnc,cdim(:,7),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
        write(lname,'("t",I1.1,"_pop_grid_biomass_age")') n
        write(vname,'("t",I1.1,"_pop_grid_biomass_age")') n
        call attrib(idnc,cdim(:,7),csize(1),vname,lname,'none',0.,6500.,any_m,point_m,land_m,double_m)
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

use carbpools_m
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
integer, dimension(ifull) :: dati
real, dimension(ifull) :: datr
real(kind=8), dimension(ifull) :: dat
real(kind=8), dimension(mp_global) :: dummy_unpack
real(kind=8), dimension(:), allocatable :: dat_in
real(kind=8), dimension(:,:), allocatable :: datpatch
real(kind=8), dimension(:,:), allocatable :: datage
real(kind=8), dimension(:,:,:), allocatable :: datpc
character(len=80) vname
logical, intent(in) :: local
  
if ( itype==-1 ) then !just for restart file
  do n = 1,maxtile
    dat = 0
    if ( n<=maxnb ) then
      dummy_unpack = real(cveg,8)    
      call cable_unpack(dummy_unpack,dat,n)
    end if
    write(vname,'("t",I1.1,"_cvc")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
  end do      
  
  ! soil temperature check
  if ( mp_global>0 ) then
    if ( any(ssnow%tgg>400.) ) then
      write(6,*) "ERROR: Invalid CABLE temperature when writing tile"
      write(6,*) "ssnow%tgg ",maxval(ssnow%tgg)
      stop -1
    end if
  end if
  
  do n = 1,maxtile  ! tile
    datr = 0.
    if ( n<=maxnb ) call cable_unpack(sv,datr,n)
    dat = real( datr, 8 )
    write(vname,'("t",I1.1,"_svs")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
    do k = 1,ms     ! soil layer
      dat = real(tgg(:,k),8)
      if ( n<=maxnb ) call cable_unpack(ssnow%tgg(:,k),dat,n)
      write(vname,'("t",I1.1,"_tgg",I1.1)') n,k
      call histwrt(dat,vname,idnc,iarch,local,.true.)
    end do
    do k = 1,ms
      dat = real(wb(:,k),8)
      if ( n<=maxnb ) call cable_unpack(ssnow%wb(:,k),dat,n)
      write(vname,'("t",I1.1,"_wb",I1.1)') n,k
      call histwrt(dat,vname,idnc,iarch,local,.true.)
    end do
    do k = 1,ms
      dat = real(wbice(:,k),8)
      if ( n<=maxnb ) call cable_unpack(ssnow%wbice(:,k),dat,n)
      write(vname,'("t",I1.1,"_wbice",I1.1)') n,k
      call histwrt(dat,vname,idnc,iarch,local,.true.)
    end do
    do k = 1,3 ! snow layer
      dat = real(tggsn(:,k),8)
      if ( n<=maxnb ) call cable_unpack(ssnow%tggsn(:,k),dat,n)
      write(vname,'("t",I1.1,"_tggsn",I1.1)') n,k
      call histwrt(dat,vname,idnc,iarch,local,.true.)
    end do
    do k = 1,3
      dat = real(smass(:,k),8)
      if ( n<=maxnb ) call cable_unpack(ssnow%smass(:,k),dat,n)
      write(vname,'("t",I1.1,"_smass",I1.1)') n,k
      call histwrt(dat,vname,idnc,iarch,local,.true.)
    end do
    do k = 1,3
      dat = real(ssdn(:,k),8)
      if ( n<=maxnb ) call cable_unpack(ssnow%ssdn(:,k),dat,n)
      write(vname,'("t",I1.1,"_ssdn",I1.1)') n,k
      call histwrt(dat,vname,idnc,iarch,local,.true.)
    end do  
    do k = 1,3
      dat = real(snowd/3.,8)
      if ( n<=maxnb ) call cable_unpack(ssnow%sdepth(:,k),dat,n)
      write(vname,'("t",I1.1,"_sdepth",I1.1)') n,k
      call histwrt(dat,vname,idnc,iarch,local,.true.)
    end do
    do k = 1,3
      dat = 0.2_8
      if ( n<=maxnb ) call cable_unpack(ssnow%sconds(:,k),dat,n)
      write(vname,'("t",I1.1,"_sconds",I1.1)') n,k
      call histwrt(dat,vname,idnc,iarch,local,.true.)
    end do
    dat = real(ssdnn,8)
    if ( n<=maxnb ) call cable_unpack(ssnow%ssdnn,dat,n)
    write(vname,'("t",I1.1,"_ssdnn")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
    dat = real(isflag,8)
    if ( n<=maxnb ) then
      dummy_unpack = real(ssnow%isflag,8)    
      call cable_unpack(dummy_unpack,dat,n)
    end if    
    write(vname,'("t",I1.1,"_sflag")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
    dat = real(snowd,8)
    if ( n<=maxnb ) call cable_unpack(ssnow%snowd,dat,n)
    write(vname,'("t",I1.1,"_snd")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
    dat = real(snowd,8)
    if ( n<=maxnb ) call cable_unpack(ssnow%osnowd,dat,n)
    write(vname,'("t",I1.1,"_osnd")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
    dat = real(snage,8)
    if ( n<=maxnb ) call cable_unpack(ssnow%snage,dat,n)
    write(vname,'("t",I1.1,"_snage")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
    dat = 100._8
    if ( n<=maxnb ) call cable_unpack(ssnow%rtsoil,dat,n)
    write(vname,'("t",I1.1,"_rtsoil")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
    dat = 0._8
    if ( n<=maxnb ) call cable_unpack(ssnow%GWwb,dat,n)
    write(vname,'("t",I1.1,"_GWwb")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
    dat = 0._8
    if ( n<=maxnb ) call cable_unpack(ssnow%wtd,dat,n)
    write(vname,'("t",I1.1,"_wtd")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
    dat=0._8
    if (n<=maxnb) call cable_unpack(canopy%cansto,dat,n)
    write(vname,'("t",I1.1,"_cansto")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
    dat=0.01_8 ! ustar
    if (n<=maxnb) call cable_unpack(canopy%us,dat,n)
    write(vname,'("t",I1.1,"_us")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)  
    dat=0._8
    if (n<=maxnb) call cable_unpack(ssnow%pudsto,dat,n)
    write(vname,'("t",I1.1,"_pudsto")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
    dat=0._8
    if (n<=maxnb) call cable_unpack(ssnow%wetfac,dat,n)
    write(vname,'("t",I1.1,"_wetfac")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
    dat=0._8
    if (n<=maxnb) call cable_unpack(canopy%ga,dat,n)
    write(vname,'("t",I1.1,"_ga")') n
    call histwrt(dat,vname,idnc,iarch,local,.true.)
  end do
  if ( soil_struc==1 ) then
    do n = 1,maxtile  ! tile
      dat=0._8
      if (n<=maxnb) call cable_unpack(ssnow%h0,dat,n)
      write(vname,'("t",I1.1,"_hzero")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)   
      do k = 1,ms     ! soil layer
        dat=0._8
        if (n<=maxnb) call cable_unpack(ssnow%S(:,k),dat,n)
        write(vname,'("t",I1.1,"_s",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do k = 1,ms
        dat=0._8
        if (n<=maxnb) call cable_unpack(ssnow%tsoil(:,k),dat,n)
        write(vname,'("t",I1.1,"_tsoil",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do k = 1,ms
        dat=0._8
        if (n<=maxnb) call cable_unpack(ssnow%thetai(:,k),dat,n)
        write(vname,'("t",I1.1,"_thetai",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      dat=0._8
      if (n<=maxnb) call cable_unpack(ssnow%snowliq(:,1),dat,n)
      write(vname,'("t",I1.1,"_snowliq",I1.1)') n,1
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(ssnow%tsurface,dat,n)
      write(vname,'("t",I1.1,"_tsurface")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) then
        dummy_unpack = real(ssnow%nsnow,8)  
        call cable_unpack(dummy_unpack,dat,n)
      end if  
      write(vname,'("t",I1.1,"_nsnow")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)   
      dat=0._8
      if (n<=maxnb) call cable_unpack(canopy%fwsoil,dat,n)
      write(vname,'("t",I1.1,"_fwsoil")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.) 
    end do
  end if
  if ( ccycle/=0 ) then
    do n = 1,maxtile  ! tile
      do k = 1,mplant     
        dat=0._8
        if (n<=maxnb) call cable_unpack(casapool%cplant(:,k),dat,n)
        write(vname,'("t",I1.1,"_cplant",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do k = 1,mplant
        dat=0._8
        if (n<=maxnb) call cable_unpack(casapool%nplant(:,k),dat,n)
        write(vname,'("t",I1.1,"_nplant",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do k = 1,mplant
        dat=0._8
        if (n<=maxnb) call cable_unpack(casapool%pplant(:,k),dat,n)
        write(vname,'("t",I1.1,"_pplant",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do k = 1,mlitter
        dat=0._8
        if (n<=maxnb) call cable_unpack(casapool%clitter(:,k),dat,n)
        write(vname,'("t",I1.1,"_clitter",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do k = 1,mlitter
        dat=0._8
        if (n<=maxnb) call cable_unpack(casapool%nlitter(:,k),dat,n)
        write(vname,'("t",I1.1,"_nlitter",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do  
      do k = 1,mlitter
        dat=0._8
        if (n<=maxnb) call cable_unpack(casapool%plitter(:,k),dat,n)
        write(vname,'("t",I1.1,"_plitter",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do k = 1,msoil
        dat=0._8
        if (n<=maxnb) call cable_unpack(casapool%csoil(:,k),dat,n)
        write(vname,'("t",I1.1,"_csoil",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do k = 1,msoil
        dat=0._8
        if (n<=maxnb) call cable_unpack(casapool%nsoil(:,k),dat,n)
        write(vname,'("t",I1.1,"_nsoil",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do k = 1,msoil
        dat=0._8
        if (n<=maxnb) call cable_unpack(casapool%psoil(:,k),dat,n)
        write(vname,'("t",I1.1,"_psoil",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      dat=0._8
      if (n<=maxnb) call cable_unpack(casamet%glai(:),dat,n)
      write(vname,'("t",I1.1,"_glai")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(phen%phen(:),dat,n)
      write(vname,'("t",I1.1,"_phen")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(phen%aphen(:),dat,n)
      write(vname,'("t",I1.1,"_aphen")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) then
        dummy_unpack = real(phen%phase(:),8)   
        call cable_unpack(dummy_unpack,dat,n)
      end if  
      write(vname,'("t",I1.1,"_phenphase")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) then
        dummy_unpack = real(phen%doyphase(:,3),8)    
        call cable_unpack(dummy_unpack,dat,n)
      end if  
      write(vname,'("t",I1.1,"_doyphase3")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(casapool%clabile(:),dat,n)
      write(vname,'("t",I1.1,"_clabile")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(casapool%nsoilmin(:),dat,n)
      write(vname,'("t",I1.1,"_nsoilmin")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(casapool%psoillab(:),dat,n)
      write(vname,'("t",I1.1,"_psoillab")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(casapool%psoilsorb(:),dat,n)
      write(vname,'("t",I1.1,"_psoilsorb")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(casapool%psoilocc(:),dat,n)
      write(vname,'("t",I1.1,"_psoilocc")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      do k = 1,mplant     
        dat=0._8
        if (n<=maxnb) call cable_unpack(casaflux%crmplant(:,k),dat,n)
        write(vname,'("t",I1.1,"_crmplant",I1.1)') n,k
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      dat=0._8
      if (n<=maxnb) call cable_unpack(casaflux%frac_sapwood(:),dat,n)
      write(vname,'("t",I1.1,"_fracsapwood")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(casaflux%sapwood_area(:),dat,n)
      write(vname,'("t",I1.1,"_sapwoodarea")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(casaflux%Crsoil(:),dat,n)
      write(vname,'("t",I1.1,"_crsoil")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(casaflux%cnpp(:),dat,n)
      write(vname,'("t",I1.1,"_cnpp")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(casaflux%clabloss(:),dat,n)
      write(vname,'("t",I1.1,"_clabloss")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(casaflux%crgplant(:),dat,n)
      write(vname,'("t",I1.1,"_crgplant")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(casaflux%stemnpp(:),dat,n)
      write(vname,'("t",I1.1,"_stemnpp")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(casabal%laimax(:),dat,n)
      write(vname,'("t",I1.1,"_LAImax")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(casabal%cleafmean(:),dat,n)
      write(vname,'("t",I1.1,"_Cleafmean")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(casabal%crootmean(:),dat,n)
      write(vname,'("t",I1.1,"_Crootmean")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(canopy%fpn(:),dat,n)
      write(vname,'("t",I1.1,"_fpn")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      dat=0._8
      if (n<=maxnb) call cable_unpack(canopy%frday(:),dat,n)
      write(vname,'("t",I1.1,"_frday")') n
      call histwrt(dat,vname,idnc,iarch,local,.true.)
    end do
  end if
end if ! itype==-1
if ( cable_pop==1 ) then
  allocate( datpatch(ifull,POP_NPATCH) )  
  allocate( datage(ifull,POP_AGEMAX) )  
  allocate( datpc(ifull,POP_NPATCH,POP_NCOHORT) )
  np_pop = size( pop%pop_grid(:) )
  allocate( dat_in(np_pop) )
  do n = 1,maxtile
    datpatch = 0._8
    datage = 0._8
    datpc = 0._8
    dat = 0._8
    dat_in = 0._8
    if ( diaglevel_pop > 0 ) then
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%cmass_sum 
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_cmass_sum")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
    end if
    if ( itype==-1 ) then !just for restart file
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%cmass_sum_old  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_cmass_sum_old")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%cheartwood_sum 
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_cheartwood_sum")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%csapwood_sum
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_csapwood_sum")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%csapwood_sum_old  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_csapwood_sum_old")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%densindiv  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_densindiv")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%height_mean  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_height_mean")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%height_max  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_height_max")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%basal_area  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_basal_area")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%sapwood_loss  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_sapwood_loss")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%sapwood_area_loss  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if
      write(vname,'("t",I1.1,"_pop_grid_sapwood_area_loss")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%stress_mortality  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_stress_mortality")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%crowding_mortality  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if
      write(vname,'("t",I1.1,"_pop_grid_crowding_mortality")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%fire_mortality  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if
      write(vname,'("t",I1.1,"_pop_grid_fire_mortality")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%cat_mortality  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_cat_mortality")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%res_mortality  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_res_mortality")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%growth  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_growth")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%area_growth  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_area_growth")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%crown_cover  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_crown_cover")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%crown_area  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_crown_area")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%crown_volume  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_crown_volume")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%sapwood_area  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_sapwood_area")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%sapwood_area_old  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_sapwood_area_old")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        dat_in(1:np_pop) = pop%pop_grid(:)%KClump  
        call pop_unpack(dat_in(1:np_pop),dat,n)
      end if  
      write(vname,'("t",I1.1,"_pop_grid_KClump")') n  
      call histwrt(dat,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then      
        do k = 1,POP_AGEMAX
          dat_in(1:np_pop) = pop%pop_grid(:)%freq_age(k)  
          call pop_unpack(dat_in(1:np_pop),datage(:,k),n)
        end do
      end if          
      write(vname,'("t",I1.1,"_pop_grid_freq_age")') n
      call histwrt(datage,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then      
        do k = 1,POP_AGEMAX
          dat_in(1:np_pop) = pop%pop_grid(:)%biomass_age(k)  
          call pop_unpack(dat_in(1:np_pop),datage(:,k),n)
        end do
      end if          
      write(vname,'("t",I1.1,"_pop_grid_biomass_age")') n
      call histwrt(datage,vname,idnc,iarch,local,.true.)
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then
          dat_in(1:np_pop) = pop%pop_grid(:)%biomass(ll)
          call pop_unpack(dat_in(1:np_pop),dat,n)
        end if  
        write(vname,'("t",I1.1,"_pop_grid_biomass",I1.1)') n,ll  
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then
          dat_in(1:np_pop) = pop%pop_grid(:)%density(ll)  
          call pop_unpack(dat_in(1:np_pop),dat,n)
        end if  
        write(vname,'("t",I1.1,"_pop_grid_density",I1.1)') n,ll  
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then
          dat_in(1:np_pop) = pop%pop_grid(:)%hmean(ll)  
          call pop_unpack(dat_in(1:np_pop),dat,n)
        end if  
        write(vname,'("t",I1.1,"_pop_grid_hmean",I1.1)') n,ll  
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then
          dat_in(1:np_pop) = pop%pop_grid(:)%hmax(ll)  
          call pop_unpack(dat_in(1:np_pop),dat,n)
        end if  
        write(vname,'("t",I1.1,"_pop_grid_hmax",I1.1)') n,ll  
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do hh = 1,POP_HEIGHT_BINS
        if ( n<=maxnb ) then
          dat_in(1:np_pop) = pop%pop_grid(:)%cmass_stem_bin(hh)  
          call pop_unpack(dat_in(1:np_pop),dat,n)
        end if
        write(vname,'("t",I1.1,"_pop_grid_cmass_stem_bin",I2.2)') n,hh  
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do hh = 1,POP_HEIGHT_BINS
        if ( n<=maxnb ) then
          dat_in(1:np_pop) = pop%pop_grid(:)%densindiv_bin(hh)  
          call pop_unpack(dat_in(1:np_pop),dat,n)
        end if  
        write(vname,'("t",I1.1,"_pop_grid_densindiv_bin",I2.2)') n,hh  
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do hh = 1,POP_HEIGHT_BINS
        if ( n<=maxnb ) then
          dat_in(1:np_pop) = pop%pop_grid(:)%height_bin(hh)  
          call pop_unpack(dat_in(1:np_pop),dat,n)
        end if  
        write(vname,'("t",I1.1,"_pop_grid_height_bin",I2.2)') n,hh
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do hh = 1,POP_HEIGHT_BINS
        if ( n<=maxnb ) then
          dat_in(1:np_pop) = pop%pop_grid(:)%diameter_bin(hh)  
          call pop_unpack(dat_in(1:np_pop),dat,n)
        end if  
        write(vname,'("t",I1.1,"_pop_grid_diameter_bin",I2.2)') n,hh 
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
      do dd = 1,POP_NDISTURB
        if ( n<=maxnb ) then
          dat_in(1:np_pop) = real(pop%pop_grid(:)%n_age(dd),8)
          call pop_unpack(dat_in(1:np_pop),dat,n)
        end if  
        write(vname,'("t",I1.1,"_pop_grid_n_age",I1.1)') n,dd 
        call histwrt(dat,vname,idnc,iarch,local,.true.)
      end do
    end if
    if ( itype==-1 .or. diaglevel_pop>=9 ) then !just for restart file
      if ( n<=maxnb ) then  
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = real(pop%pop_grid(:)%patch(k)%id,8)
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_id")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%freq(k)  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_freq")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%freq_old(k)   
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if
      write(vname,'("t",I1.1,"_pop_grid_freq_old")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%factor_recruit  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_factor_recruit")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%pgap  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_pgap")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%lai  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_lai")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%biomass  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_biomass")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%biomass_old  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_biomass_old")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)    
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%sapwood  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%heartwood  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_heartwood")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%sapwood_old  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_old")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%sapwood_area  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_area")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%sapwood_area_old  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_area_old")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%stress_mortality   
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_stress_mortality")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%fire_mortality  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_fire_mortality")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%cat_mortality  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_cat_mortality")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%crowding_mortality 
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_crowding_mortality")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%cpc  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_cpc")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%mortality  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_mortality")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%sapwood_loss  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_loss")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%sapwood_area_loss  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_sapwood_area_loss")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%growth  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_growth")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%area_growth  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_area_growth")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%frac_NPP  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_frac_NPP")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%frac_respiration  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_frac_respiration")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      if ( n<=maxnb ) then
        do k = 1,POP_NPATCH
          dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%frac_light_uptake  
          call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
        end do
      end if  
      write(vname,'("t",I1.1,"_pop_grid_patch_frac_light_uptake")') n
      call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      do dd = 1,POP_NDISTURB  
        if ( n<=maxnb ) then  
          do k = 1,POP_NPATCH  
            dat_in(1:np_pop) = real(pop%pop_grid(:)%patch(k)%disturbance_interval(dd),8)
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_disturbance_interval",I1.1)') n,dd
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do dd = 1,POP_NDISTURB 
        if ( n<=maxnb ) then  
          do k = 1,POP_NPATCH  
            dat_in(1:np_pop) = real(pop%pop_grid(:)%patch(k)%first_disturbance_year(dd),8)
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_first_disturbance_year",I1.1)') n,dd
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do dd = 1,POP_NDISTURB 
        if ( n<=maxnb ) then  
          do k = 1,POP_NPATCH  
            dat_in(1:np_pop) = real(pop%pop_grid(:)%patch(k)%age(dd),8)
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_age",I1.1)') n,dd
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do dd = 1,POP_NDISTURB  
        if ( n<=maxnb ) then  
          do k = 1,POP_NPATCH  
            dat_in(1:np_pop) = real(pop%pop_grid(:)%ranked_age_unique(k,dd),8)  
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_ranked_age_unique",I1.1)') n,dd
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do dd = 1,POP_NDISTURB
        if ( n<=maxnb ) then  
          do k = 1,POP_NPATCH  
            dat_in(1:np_pop) = pop%pop_grid(:)%freq_ranked_age_unique(k,dd)  
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_freq_ranked_age_unique",I1.1)') n,dd
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do k = 1,POP_NPATCH  
            dat_in(1:np_pop) = real(pop%pop_grid(:)%patch(k)%layer(ll)%ncohort,8)
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
          end do
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_ncohort")') n,ll
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do k = 1,POP_NPATCH  
            dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%biomass  
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
          end do
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_biomass")') n,ll
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do k = 1,POP_NPATCH  
            dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%density  
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
          end do
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_density")') n,ll
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do k = 1,POP_NPATCH  
            dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%hmean 
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
          end do
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_hmean")') n,ll
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do k = 1,POP_NPATCH  
            dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%hmax  
            call pop_unpack(dat_in(1:np_pop),datpatch(:,k),n)
          end do
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_hmax")') n,ll
        call histwrt(datpatch,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do cc = 1,POP_NCOHORT            
            do k = 1,POP_NPATCH  
              dat_in(1:np_pop) = real(pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%age,8)
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do  
          end do  
        end if
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_age")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do cc = 1,POP_NCOHORT 
            do k = 1,POP_NPATCH  
              dat_in(1:np_pop) = real(pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%id,8)
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do  
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_id")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do cc = 1,POP_NCOHORT            
            do k = 1,POP_NPATCH  
              dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%biomass  
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do  
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_biomass")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do cc = 1,POP_NCOHORT    
            do k = 1,POP_NPATCH  
              dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%density  
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do  
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_density")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do cc = 1,POP_NCOHORT  
            do k = 1,POP_NPATCH  
              dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_resource_uptake  
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do  
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_resource_uptake")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do cc = 1,POP_NCOHORT    
            do k = 1,POP_NPATCH  
              dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_light_uptake  
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do  
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_light_uptake")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do cc = 1,POP_NCOHORT            
            do k = 1,POP_NPATCH  
              dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_interception  
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do  
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_interception")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do cc = 1,POP_NCOHORT    
            do k = 1,POP_NPATCH  
              dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_respiration 
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do  
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_respiration")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do cc = 1,POP_NCOHORT    
            do k = 1,POP_NPATCH  
              dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%frac_NPP  
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_frac_NPP")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do cc = 1,POP_NCOHORT    
            do k = 1,POP_NPATCH
              dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%respiration_scalar  
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do  
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_respiration_scalar")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do cc = 1,POP_NCOHORT            
            do k = 1,POP_NPATCH  
              dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%crown_area  
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do  
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_crown_area")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then
          do cc = 1,POP_NCOHORT    
            do k = 1,POP_NPATCH  
              dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%Pgap  
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do  
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Pgap")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do cc = 1,POP_NCOHORT  
            do k = 1,POP_NPATCH  
              dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%height  
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do  
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_height")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do cc = 1,POP_NCOHORT    
            do k = 1,POP_NPATCH 
              dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%diameter  
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do 
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_diameter")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then 
          do cc = 1,POP_NCOHORT    
            do k = 1,POP_NPATCH  
              dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%sapwood  
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do  
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_sapwood")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do cc = 1,POP_NCOHORT    
            do k = 1,POP_NPATCH  
              dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%heartwood  
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do  
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_heartwood")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do cc = 1,POP_NCOHORT            
            do k = 1,POP_NPATCH
              dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%sapwood_area  
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do  
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_sapwood_area")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do cc = 1,POP_NCOHORT    
            do k = 1,POP_NPATCH  
              dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%basal_area  
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do  
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_basal_area")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do cc = 1,POP_NCOHORT    
            do k = 1,POP_NPATCH  
              dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%LAI  
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do  
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_LAI")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do cc = 1,POP_NCOHORT    
            do k = 1,POP_NPATCH 
              dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%Cleaf  
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do  
          end do  
        end if  
        write(vname,'("t",I1.1,"_pop_grid_patch_layer",I1.1,"_cohort_Cleaf")') n,ll
        call histwrt(datpc,vname,idnc,iarch,local,.true.)
      end do
      do ll = 1,POP_NLAYER
        if ( n<=maxnb ) then  
          do cc = 1,POP_NCOHORT    
            do k = 1,POP_NPATCH  
              dat_in(1:np_pop) = pop%pop_grid(:)%patch(k)%layer(ll)%cohort(cc)%Croot  
              call pop_unpack(dat_in(1:np_pop),datpc(:,k,cc),n)
            end do  
          end do  
        end if  
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

