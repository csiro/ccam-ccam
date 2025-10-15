! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2025 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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

! Main NetCDF input routines.  Host grid is automatically
! interpolated to nested model grid.  Three options are
!   nested=0  Initial conditions
!   nested=1  Nudging fields
!   nested=2  Surface data recycling
!   nested=3  Ensemble fields
!   nested=4  Surface data recycling (without biosphere)
      
! This version supports the parallel file routines contained
! in infile.f90.  Hence, restart files do not require any
! gathers and scatters.
    
! Thanks to Paul Ryan for optimising NetCDF routines
    
module onthefly_m
    
implicit none

private
public onthefly, retopo, onthefly_exit
    
integer, save :: ik, jk, kk, ok, nsibx                        ! input grid size
integer fwsize                                                ! size of temporary arrays
integer, save :: nemi = -1                                    ! land-sea mask method (3=soilt, 2=zht, 1=ocndepth, -1=fail)
integer, save :: fill_land = 0                                ! number of iterations required for land fill
integer, save :: fill_floor = 0                               ! number of iterations required for floor fill (3d ocean)
integer, save :: fill_floorlake = 0                           ! number of iterations required for floor+lake fill (3d ocean)
integer, save :: fill_sea = 0                                 ! number of iterations required for ocean fill
integer, save :: fill_nourban = 0                             ! number of iterations required for urban fill
integer, save :: native_ccam = 0                              ! is host CCAM (native_ccam=1) or cdfivdar (native_ccam=0)
integer, save :: xx4_win, yy4_win                             ! shared memory windows
integer, save :: xy4_count = 0                                ! counter for new allocations of memory windows
integer, dimension(3), save :: xx4save_win, yy4save_win       ! store old memory windows for later deallocation
integer, dimension(:,:), allocatable, save :: nface4          ! interpolation panel index
integer, dimension(:,:,:), allocatable, save :: pijn2pw       ! array for unpacking buffer
integer, dimension(:,:,:), allocatable, save :: pijn2pn       ! array for unpacking buffer
real, save :: rlong0x, rlat0x, schmidtx                       ! input grid coordinates
real, dimension(3,3), save :: rotpoles, rotpole               ! vector rotation data
real, dimension(:,:), allocatable, save :: xg4, yg4           ! interpolation coordinate indices
real, dimension(:), allocatable, save :: axs_a, ays_a, azs_a  ! vector rotation data (single file)
real, dimension(:), allocatable, save :: bxs_a, bys_a, bzs_a  ! vector rotation data (single file)
real, dimension(:), allocatable, save :: axs_w, ays_w, azs_w  ! vector rotation data (multiple files)
real, dimension(:), allocatable, save :: bxs_w, bys_w, bzs_w  ! vector rotation data (multiple files)
real, dimension(:), allocatable, save :: sigin                ! input vertical coordinates
real, dimension(:), allocatable, save :: gosig_1              ! input ocean reference levels
real, dimension(:), allocatable, save :: gosig_h              ! input ocean reference half levels
real, dimension(:,:), allocatable, save :: gosig_3            ! input ocean 3d levels
real(kind=8), dimension(:,:), pointer, save :: xx4, yy4       ! shared arrays used for interpolation
logical iotest, newfile                                       ! tests for interpolation and new metadata
logical allowtrivialfill                                      ! special case where trivial data is allowed

interface fill_cc4
  module procedure fill_cc4_3d, fill_cc4_1
end interface

contains

! *****************************************************************************
! Main interface for input data that reads grid metadata
    
subroutine onthefly(nested,kdate_r,ktime_r,psl,zss,tss,sicedep,fracice,t,u,v,qg,tgg,wb,wbice,snowd,qfg, &
                    qlg,qrg,qsng,qgrg,ni,nr,ns,tggsn,smass,ssdn,ssdnn,snage,isflag,mlodwn,ocndwn,xtgdwn)

use cc_mpi           ! CC MPI routines
use darcdf_m         ! Netcdf data
use infile           ! Input file routines
use mlo_ctrl         ! Ocean physics control layer
use newmpar_m        ! Grid parameters
use parm_m           ! Model configuration
use soil_m           ! Soil and surface data
use stime_m          ! File date data

integer, parameter :: nihead = 54
integer, parameter :: nrhead = 14

integer, intent(in) :: nested
integer, intent(out) :: kdate_r, ktime_r
integer, save :: maxarchi
integer ierx, idvtime
integer kdate_rsav, ktime_rsav
integer, dimension(nihead) :: nahead
integer, dimension(:), intent(out) :: isflag
integer(kind=8) mtimer
real(kind=8) timer
real, dimension(:,:,:), intent(out) :: mlodwn
real, dimension(:,:,:), intent(out) :: xtgdwn
real, dimension(:,:), intent(out) :: wb, wbice, tgg
real, dimension(:,:), intent(out) :: tggsn, smass, ssdn
real, dimension(:,:), intent(out) :: ocndwn
real, dimension(:,:), intent(out) :: t, u, v, qg, qfg, qlg, qrg, qsng, qgrg
real, dimension(:,:), intent(out) :: ni, nr, ns
real, dimension(:), intent(out) :: psl, zss, tss, fracice, snowd
real, dimension(:), intent(out) :: sicedep, ssdnn, snage
real, dimension(nrhead) :: ahead
real, dimension(11) :: rdum
logical, save :: firstcall = .true.
logical ltest
character(len=80) datestring
character(len=120) :: versionstring = ' '

call START_LOG(onthefly_begin)

if ( myid==0 ) then
  write(6,*) 'Entering onthefly for nested,ktau = ',nested,ktau
end if  

!--------------------------------------------------------------------
! pfall indicates all processors have a parallel input file and there
! is no need to broadcast metadata (see infile.f90).  Otherwise read
! metadata on myid=0 and broadcast that data to all processors.
if ( myid==0 .or. pfall ) then
    
  ! Locate new file and read grid metadata --------------------------
  if ( ncid/=ncidold ) then
    if ( myid==0 ) then
      write(6,*) '-> Reading new file metadata'
    end if  
    iarchi = 1    ! default time index for input file
    maxarchi = 0  ! default number of timesteps in input file
    ik = pil_g    ! grid size
    jk = pjl_g    ! grid size
    kk = pka_g    ! atmosphere vertical levels
    ok = pko_g    ! ocean vertical levels
    call ccnf_get_attg(ncid,'rlong0',rlong0x,ierr=ierx)
    if ( ierx==0 ) then
      ! New global attribute method
      call ccnf_get_attg(ncid,'rlat0',rlat0x)
      call ccnf_get_attg(ncid,'schmidt',schmidtx)
      call ccnf_get_attg(ncid,'nsib',nsibx)
    else
      ! Old int_header and real_header method      
      call ccnf_get_attg(ncid,'int_header',nahead)
      call ccnf_get_attg(ncid,'real_header',ahead)
      nsibx    = nahead(44) ! land-surface parameterisation
      rlong0x  = ahead(5)   ! longitude
      rlat0x   = ahead(6)   ! latitude
      schmidtx = ahead(7)   ! schmidt factor
      if ( schmidtx<=0. .or. schmidtx>1. ) then
        ! backwards compatibility option
        rlong0x  = ahead(6)
        rlat0x   = ahead(7)
        schmidtx = ahead(8)
        if ( schmidtx<=0. .or. schmidtx>1. ) then
          write(6,*) "ERROR: Invalid schmidt in input file"
          call ccmpi_abort(-1)
        end if  
      endif  ! (schmidtx<=0..or.schmidtx>1.)
    end if   ! ierx==0 ..else..
    call ccnf_get_attg(ncid,'version',versionstring,ierr=ierx)
    if ( ierx==0 ) then
      native_ccam = 1 ! found CCAM host data
    else
      native_ccam = 0 ! non-native CCAM (possibly cdfvidar)
    end if
    call ccnf_inq_dimlen(ncid,'time',maxarchi)
    if ( myid==0 ) then
      write(6,*) "-> Found ik,jk,kk,ok ",ik,jk,kk,ok
      write(6,*) "->      maxarchi ",maxarchi
      write(6,*) "->      rlong0x,rlat0x,schmidtx ",rlong0x,rlat0x,schmidtx
      write(6,*) "->      native_ccam ",native_ccam
    end if
  end if
  
  ! search for required date ----------------------------------------
  if ( (nrungcm==-14.or.nrungcm==-17) .and. (nested==2.or.nested==4) ) then
    iarchi = 1
    ltest = .true.
    kdate_r = kdate_s 
    ktime_r = ktime_s
    if ( myid==0 ) then
      write(6,*) 'Reading climatology for iarch=1'
    end if
  else  
    if ( myid==0 ) then
      write(6,*) 'Search for kdate_s,ktime_s >= ',kdate_s,ktime_s
    end if
    ltest = .true.       ! flag indicates that the date is not yet found
    iarchi = iarchi - 1  ! move time index back one step to check current position in file
    call ccnf_inq_varid(ncid,'time',idvtime)
    call ccnf_get_att(ncid,idvtime,'units',datestring)
    call processdatestring(datestring,kdate_rsav,ktime_rsav)
    ! start search for required date/time
    do while( ltest .and. iarchi<maxarchi )
      ! could read this as one array, but we only usually need to advance 1 step
      iarchi = iarchi + 1
      kdate_r = kdate_rsav
      ktime_r = ktime_rsav
      call ccnf_get_vara(ncid,idvtime,iarchi,timer)
      mtimer = nint(timer,8)
      call datefix(kdate_r,ktime_r,mtimer)
      ! ltest = .false. when correct date is found
      ltest = (2400*(kdate_r-kdate_s)-1200*nsemble+ktime_r-ktime_s)<0
    end do
    if ( nsemble/=0 ) then
      kdate_r = kdate_s
      ktime_r = ktime_s
    end if
    if ( ltest ) then
      ! ran out of file before correct date was located
      ktime_r = -1
      if ( myid==0 ) then
        write(6,*) 'Search failed with ltest,iarchi =',ltest, iarchi
        write(6,*) '                kdate_r,ktime_r =',kdate_r, ktime_r 
      end if
    else
      ! valid date found  
      if ( nested==1 .and. firstcall ) then
        ! move back to previous time-step  
        firstcall = .false.
        iarchi = max( iarchi-1, 1 )
        kdate_r = kdate_rsav
        ktime_r = ktime_rsav
        call ccnf_get_vara(ncid,idvtime,iarchi,timer)
        mtimer = nint(timer,8)
        call datefix(kdate_r,ktime_r,mtimer)
      end if  
      if ( myid==0 ) then
        write(6,*) 'Search succeeded with ltest,iarchi =',ltest, iarchi
        write(6,*) '                   kdate_r,ktime_r =',kdate_r, ktime_r
      end if
    end if
  end if  

endif  ! ( myid==0 .or. pfall )

! if metadata is not read by all processors, then broadcast ---------
if ( .not.pfall ) then
  rdum(1) = rlong0x
  rdum(2) = rlat0x
  rdum(3) = schmidtx
  ! kdate_r is too large to represent as a single real, so
  ! we split kdate_r into year, month and day
  rdum(4) = real(kdate_r/10000)
  rdum(5) = real(kdate_r/100-nint(rdum(4))*100)
  rdum(6) = real(kdate_r-nint(rdum(4))*10000-nint(rdum(5))*100)
  rdum(7) = real(ktime_r)
  if ( ncid/=ncidold ) then
    rdum(8) = 1.
  else
    rdum(8) = 0.
  end if
  rdum(9)  = real(iarchi)
  rdum(10) = real(nsibx)
  rdum(11) = real(native_ccam)
  call ccmpi_bcast(rdum(1:11),0,comm_world)
  rlong0x     = rdum(1)
  rlat0x      = rdum(2)
  schmidtx    = rdum(3)
  kdate_r     = nint(rdum(4))*10000+nint(rdum(5))*100+nint(rdum(6))
  ktime_r     = nint(rdum(7))
  newfile     = nint(rdum(8))==1
  iarchi      = nint(rdum(9))
  nsibx       = nint(rdum(10))
  native_ccam = nint(rdum(11))
  ik       = pil_g      ! grid size
  jk       = pjl_g      ! grid size
  kk       = pka_g      ! atmosphere vertical levels
  ok       = pko_g      ! ocean vertical levels
else
  newfile = (ncid/=ncidold)
end if

! mark current file as read for metadata
ncidold = ncid

! trap error if correct date/time is not located --------------------
if ( ktime_r<0 ) then
  if ( nested==2 .or. nested==4 ) then
    if ( myid==0 ) write(6,*) "WARN: Cannot locate date/time in input file"
    return
  else
    write(6,*) "ERROR: Cannot locate date/time in input file"
    write(6,*) "myid,pfall,ktime_r ",myid,pfall,ktime_r
    call ccmpi_abort(-1)
  end if
end if
!--------------------------------------------------------------------
      
! Here we call ontheflyx
   
! Note that if histrd fails to find a variable, it returns zero in
! the output array

call onthefly_work(nested,kdate_r,ktime_r,psl,zss,tss,sicedep,fracice,t,u,v,qg,tgg,wb,wbice, &
                   snowd,qfg,qlg,qrg,qsng,qgrg,ni,nr,ns,tggsn,smass,ssdn,ssdnn,snage,isflag, &
                   mlodwn,ocndwn,xtgdwn)

if ( myid==0 ) write(6,*) "Leaving onthefly"

call END_LOG(onthefly_end)

return
end subroutine onthefly


! *****************************************************************************
! Read data from netcdf file
      
! Input usually consists of either a single input file that is
! scattered across processes, or multiple input files that are read
! by many processes.  In the case of restart files, then there is
! no need for message passing.
subroutine onthefly_work(nested,kdate_r,ktime_r,psl,zss,tss,sicedep,fracice,t,u,v,qg,tgg,wb,wbice, &
                         snowd,qfg,qlg,qrg,qsng,qgrg,ni,nr,ns,tggsn,smass,ssdn,ssdnn,snage,isflag, &
                         mlodwn,ocndwn,xtgdwn)
      
use aerointerface, only : opticaldepth         ! Aerosol interface          
use cc_mpi                                     ! CC MPI routines
use cfrac_m                                    ! Cloud fraction
use const_phys                                 ! Physical constants
use darcdf_m                                   ! Netcdf data
use extraout_m                                 ! Additional diagnostics      
use histave_m, only : cbas_ave,ctop_ave,      &     
    wb_ave,tscr_ave,rhscr_ave                  ! Time average arrays
use infile                                     ! Input file routines
use kuocom_m                                   ! JLM convection
use latlong_m                                  ! Lat/lon coordinates
use latltoij_m                                 ! Lat/Lon to cubic ij conversion
use mlodynamics                                ! Ocean dynamics
use morepbl_m                                  ! Additional boundary layer diagnostics
use newmpar_m                                  ! Grid parameters
use nharrs_m, only : phi_nh,lrestart,         &
    lrestart_radiation, lrestart_tracer,      &
    lrestart_convection                        ! Non-hydrostatic atmosphere arrays
use nsibd_m, only : isoilm,isoilm_in,rsmin,   &
    carb_plant,carb_litter,carb_soil           ! Land-surface arrays
use parm_m                                     ! Model configuration
use parmdyn_m                                  ! Dynamics parmaters
use parmgeom_m                                 ! Coordinate data
use prec_m, only : precip,precc                ! Precipitation
use raddiag_m                                  ! Radiation diagnostic
use riverarrays_m                              ! River data
use savuvt_m                                   ! Saved dynamic arrays
use savuv1_m                                   ! Saved dynamic arrays
use screen_m                                   ! Screen level diagnostics
use setxyz_m                                   ! Define CCAM grid
use sflux_m                                    ! Surface flux routines
use sigs_m                                     ! Atmosphere sigma levels
use soil_m                                     ! Soil and surface data
use soilv_m                                    ! Soil parameters
use stime_m                                    ! File date data
use tkeeps, only : tke,eps,u_ema,v_ema,w_ema, &
    thetal_ema,qv_ema,ql_ema,qf_ema,cf_ema,   &
    tke_ema                                    ! TKE-EPS boundary layer
use tracers_m                                  ! Tracer data
use utilities                                  ! Grid utilities
use vecsuv_m                                   ! Map to cartesian coordinates
use vvel_m, only : dpsldt,sdot                 ! Additional vertical velocity
use workglob_m                                 ! Additional grid interpolation
use work2_m                                    ! Diagnostic arrays
use xarrs_m, only : pslx                       ! Saved dynamic arrays

real, parameter :: iotol = 1.E-5               ! tolarance for iotest grid matching
real, parameter :: aerosol_tol = 1.e-4         ! tolarance for aerosol data
      
integer, intent(in) :: nested, kdate_r, ktime_r
integer idv, retopo_test, ktest
integer levk, levkin, ier, igas
integer i, j, k, mm, iq, ifrac, n
integer, dimension(:), intent(out) :: isflag
integer, dimension(11+3*ms) :: ierc
integer, dimension(12), save :: iers
integer, dimension(2) :: shsize
real mxd_o, x_o, y_o, al_o, bt_o, depth_hl_xo, depth_hl_yo
real, dimension(:,:,:), intent(out) :: mlodwn
real, dimension(:,:,:), intent(out) :: xtgdwn
real, dimension(:,:), intent(out) :: ocndwn
real, dimension(:,:), intent(out) :: wb, wbice, tgg
real, dimension(:,:), intent(out) :: tggsn, smass, ssdn
real, dimension(:,:), intent(out) :: t, u, v, qg, qfg, qlg, qrg, qsng, qgrg
real, dimension(:,:), intent(out) :: ni, nr, ns
real, dimension(:), intent(out) :: psl, zss, tss, fracice
real, dimension(:), intent(out) :: snowd, sicedep, ssdnn, snage
real, dimension(ifull) :: dum6, tss_l, tss_s, pmsl, depth !, opldep
real, dimension(ifull) :: duma
real, dimension(ifull,6) :: udum6
real, dimension(ifull,wlev) :: oo
real, dimension(:,:), allocatable :: ucc6
real, dimension(:), allocatable :: ucc
real, dimension(:), allocatable :: fracice_a, sicedep_a
real, dimension(:), allocatable :: tss_l_a, tss_s_a, tss_a
real, dimension(:), allocatable :: t_a_lev, psl_a
real, dimension(:), allocatable, save :: zss_a, ocndep_a !, opldep_a
real, dimension(:), allocatable, save :: wts_a  ! not used here or defined in call setxyz
real, dimension(:,:), allocatable, save :: gosig3_a
real, dimension(kk+ok+13) :: dumr
real(kind=8), dimension(:), pointer, save :: z_a, x_a, y_a
real(kind=8), dimension(ifull) :: dum6r8
character(len=20) vname
character(len=4) trnum
logical, dimension(ms) :: tgg_found, wetfrac_found, wb_found
logical tss_test, tst
logical mixr_found, siced_found, fracice_found, soilt_found
logical u10_found, carbon_found, mlo_found, mlo2_found, mloice_found
logical zht_needed, zht_found, urban1_found, urban2_found
logical aero_found, nllp_found, carbon2_found !, mlo3_found
logical, dimension(:), allocatable, save :: land_a, landlake_a, sea_a, nourban_a
logical, dimension(:,:), allocatable, save :: land_3d
logical, dimension(:,:), allocatable, save :: landlake_3d

! iotest      indicates no interpolation required
! ptest       indicates the grid decomposition of the mesonest file is the same as the model, including the same number of processes
! tss_test    indicates that iotest is true, as well as seaice fraction and seaice depth are present in the input file
! retopo_test indicates that topography adjustment for t, q and psl is required
! fnresid     is the number of processes reading input files.
! fncount     is the number of files read on a process.
! fnproc      is fnresid*fncount or the total number of input files to be read.  fnproc=1 indicates a single input file
! fwsize      is the size of the array for reading input data.  fwsize>0 implies this process id is reading data
! allowtrivialfill is for when the input file has no land-sea mask (e.g., output from one.f)

! default surface height
zss = 0.

! memory needed to read input files
fwsize = pipan*pjpan*pnpan*mynproc 

! default trivial fill is off
allowtrivialfill = .false.
      
! test if retopo fields are required
if ( nud_p==0 .and. nud_t==0 .and. nud_q==0 ) then
  retopo_test = 0
else
  retopo_test = 1
end if
      
! Determine if interpolation is required
iotest = 6*ik*ik==ifull_g .and. abs(rlong0x-rlong0)<iotol .and. abs(rlat0x-rlat0)<iotol .and. &
         abs(schmidtx-schmidt)<iotol .and. (nsib==nsibx.or.nested==1.or.nested==3) .and.      &
         native_ccam==1
if ( abs(nmlo)>=1 .and. abs(nmlo)<=9 ) then
  iotest = iotest .and. wlev==ok
end if

if ( iotest ) then
  io_in = 1   ! no interpolation
  if ( .not.(iotest.and.ptest) ) then
    ! this is a special case, such as when the number of processes changes during an experiment  
    if ( myid==0 ) then
      write(6,*) "Redistribution is required with iotest,ptest,io_in =",iotest, ptest, io_in
    end if
  else    
    if ( myid==0 ) then
      write(6,*) "Interpolation is not required with iotest,io_in =",iotest, io_in
    end if
  end if  
else
  io_in = -1  ! interpolation
  if ( myid==0 ) then
    write(6,*) "Interpolation is required with iotest,io_in =",iotest, io_in
  end if  
end if
  

!--------------------------------------------------------------------
! Allocate interpolation, vertical level and mask arrays
if ( .not.allocated(nface4) ) then
  allocate( nface4(ifull,4), xg4(ifull,4), yg4(ifull,4) )
end if
if ( newfile ) then
  if ( allocated(sigin) ) then
    deallocate( sigin, gosig_1, land_a, landlake_a, land_3d, landlake_3d, sea_a, nourban_a )
    deallocate( gosig_3, gosig_h )
  end if
  allocate( sigin(kk), gosig_1(ok), land_a(fwsize), landlake_a(fwsize), land_3d(fwsize,ok) )
  allocate( landlake_3d(fwsize,ok), sea_a(fwsize), nourban_a(fwsize) )
  allocate( gosig_3(ifull,ok), gosig_h(0:ok) )
  sigin       = 1.
  gosig_1     = 1.
  gosig_3     = 1.
  land_a      = .false.
  landlake_a  = .false.
  land_3d     = .false.
  landlake_3d = .false.
  sea_a       = .true.
  nourban_a   = .false.
  ! reset fill counters
  fill_land      = 0
  fill_floor     = 0
  fill_floorlake = 0
  fill_sea       = 0
  fill_nourban   = 0
  ! reset land-sea mask search
  nemi = -1
end if
      
!--------------------------------------------------------------------
! Determine input grid coordinates and interpolation arrays
if ( newfile .and. .not.iotest ) then
  
#ifdef share_ifullg
  shsize(1) = 1 + 4*ik
  shsize(2) = 1 + 4*ik
  call ccmpi_allocshdatar8(xx4,shsize(1:2),xx4_win)
  call ccmpi_allocshdatar8(yy4,shsize(1:2),yy4_win)
  xy4_count = xy4_count + 1
#else
  allocate( xx4(1+4*ik,1+4*ik), yy4(1+4*ik,1+4*ik) )
#endif

  if ( m_fly==1 ) then
    rlong4_l(:,1) = rlongg(:)*180./pi
    rlat4_l(:,1)  = rlatt(:)*180./pi
  end if
          
  if ( myid==0 ) then
    write(6,*) "Defining input file grid"
    if ( allocated(axs_a) ) then
      deallocate( axs_a, ays_a, azs_a )
      deallocate( bxs_a, bys_a, bzs_a )          
    end if
    allocate( axs_a(ik*ik*6), ays_a(ik*ik*6), azs_a(ik*ik*6) )
    allocate( bxs_a(ik*ik*6), bys_a(ik*ik*6), bzs_a(ik*ik*6) )
    allocate( x_a(ik*ik*6), y_a(ik*ik*6), z_a(ik*ik*6) )
    allocate( wts_a(ik*ik*6) )
    !   following setxyz call is for source data geom    ****   
    do iq = 1,ik*ik*6
      axs_a(iq) = real(iq)
      ays_a(iq) = real(iq)
      azs_a(iq) = real(iq)
    end do 
    call setxyz(ik,rlong0x,rlat0x,-schmidtx,x_a,y_a,z_a,wts_a,axs_a,ays_a,azs_a,bxs_a,bys_a,bzs_a,xx4,yy4, &
                id,jd,ktau,ds)
    deallocate( x_a, y_a, z_a )
    deallocate( wts_a )
    nullify( x_a, y_a, z_a )
  end if ! (myid==0)
  
#ifdef share_ifullg
  if ( node_myid==0 ) then
    call ccmpi_bcastr8(xx4,0,comm_nodecaptain)
    call ccmpi_bcastr8(yy4,0,comm_nodecaptain)
  end if  
  call ccmpi_barrier(comm_node)
#else
  call ccmpi_bcastr8(xx4,0,comm_world)
  call ccmpi_bcastr8(yy4,0,comm_world)
#endif
  
  ! calculate the rotated coords for host and model grid
  rotpoles = calc_rotpole(rlong0x,rlat0x)
  rotpole  = calc_rotpole(rlong0,rlat0)
  if ( myid==0 .and. nmaxpr==1 ) then
    write(6,*)'m_fly,nord ',m_fly,3
    write(6,*)'kdate_r,ktime_r,ktau,ds',kdate_r,ktime_r,ktau,ds
    write(6,*)'rotpoles:'
    do i = 1,3
      write(6,'(3x,2i1,5x,2i1,5x,2i1,5x,3f8.4)') (i,j,j=1,3),(rotpoles(i,j),j=1,3)
    enddo
    write(6,*)'in onthefly rotpole:'
    do i = 1,3
      write(6,'(3x,2i1,5x,2i1,5x,2i1,5x,3f8.4)') (i,j,j=1,3),(rotpole(i,j),j=1,3)
    enddo
    write(6,*)'xx4,yy4 ',xx4(id,jd),yy4(id,jd)
    write(6,*)'before latltoij for id,jd: ',id,jd
    write(6,*)'rlong0x,rlat0x,schmidtx ',rlong0x,rlat0x,schmidtx
  end if                  ! (myid==0)

  ! setup interpolation arrays
  do mm = 1,m_fly  !  was 4, now may be set to 1 in namelist
    call latltoij(rlong4_l(:,mm),rlat4_l(:,mm),      & !input
                  rlong0x,rlat0x,schmidtx,           & !input
                  xg4(:,mm),yg4(:,mm),nface4(:,mm),  & !output (source)
                  xx4,yy4,ik)
  end do

#ifdef share_ifullg
  ! Previously deallocating memory triggered bugs in the MPI library
  ! Here we will leave the memory allocated for now and deallocate
  ! just before mpi_finalize  
  if ( xy4_count <= size(xx4save_win) ) then
    xx4save_win(xy4_count) = xx4_win  
    yy4save_win(xy4_count) = yy4_win
  end if    
  !call ccmpi_freeshdata(xx4_win)
  !call ccmpi_freeshdata(yy4_win)
#else
  deallocate( xx4, yy4 )
#endif
  nullify( xx4, yy4 )
  
  ! Define filemap for multi-file method
  ! Also deallocates axs_a, ays_a, azs_a, bxs_a, bys_a, bzs_a arrays
  call file_wininit
      
end if ! newfile .and. .not.iotest

! allocate working arrays
allocate( ucc(fwsize), tss_a(fwsize), ucc6(fwsize,6) )

! -------------------------------------------------------------------
! read time invariant data when file is first opened
if ( newfile ) then

  if ( myid==0 ) write(6,*) "Reading time invariant fields"  
    
  ! read vertical levels and missing data checks
  if ( myid==0 .or. pfall ) then
    if ( kk>1 ) then
      call ccnf_inq_varid(ncid,'lev',idv,tst)
      if ( tst ) call ccnf_inq_varid(ncid,'layer',idv,tst)
      if ( tst ) call ccnf_inq_varid(ncid,'sigma',idv,tst)
      if ( tst ) then
        write(6,*) "ERORR: multiple levels expected but no sigma data found ",kk
        call ccmpi_abort(-1)
      else
        call ccnf_get_vara(ncid,idv,1,kk,sigin)
        if ( myid==0 .and. nmaxpr==1 ) write(6,'(" sigin=",(9f7.4))') (sigin(k),k=1,kk)
      end if
    else
      sigin(:) = 1.       
    end if  
    if ( ok>0 ) then
      call ccnf_inq_varid(ncid,'olev',idv,tst)
      if ( tst ) then
        ! default for old code without olev  
        mxd_o = 5000.
        x_o = real(wlev)
        al_o = mxd_o*(x_o/mxd_o-1.)/(x_o-x_o*x_o*x_o)         ! sigma levels
        bt_o = mxd_o*(x_o*x_o*x_o/mxd_o-1.)/(x_o*x_o*x_o-x_o) ! sigma levels 
        do k = 1,wlev
          x_o = real(k-1)
          y_o = real(k)
          depth_hl_xo = (al_o*x_o**3+bt_o*x_o)/mxd_o ! ii is for half level ii-0.5
          depth_hl_yo = (al_o*y_o**3+bt_o*y_o)/mxd_o ! ii+1 is for half level ii+0.5
          gosig_1(k) = 0.5*(depth_hl_xo+depth_hl_yo)
        end do
      else
        ! usual  
        call ccnf_get_vara(ncid,idv,1,ok,gosig_1)
      end if
    end if
    ! check for missing data
    iers(1:12) = 0
    call ccnf_inq_varid(ncid,'mixr',idv,tst)
    if ( tst ) iers(1) = -1
    call ccnf_inq_varid(ncid,'siced',idv,tst)
    if ( tst ) iers(2) = -1
    call ccnf_inq_varid(ncid,'fracice',idv,tst)
    if ( tst ) iers(3) = -1
    call ccnf_inq_varid(ncid,'soilt',idv,tst)
    if ( tst ) iers(4) = -1
    call ccnf_inq_varid(ncid,'ocndepth',idv,tst)
    if ( tst ) iers(5) = -1
    call ccnf_inq_varid(ncid,'thetao',idv,tst)
    if ( tst ) iers(6) = -1
    call ccnf_inq_varid(ncid,'t1_intmtgg1',idv,tst)
    if ( tst ) iers(7) = -1
    call ccnf_inq_varid(ncid,'intmtgg1',idv,tst)
    if ( tst ) iers(8) = -1
    call ccnf_inq_varid(ncid,'uic',idv,tst)
    if ( tst ) iers(9) = -1
    call ccnf_inq_varid(ncid,'zht',idv,tst)
    if ( tst ) iers(10) = -1
    call ccnf_inq_varid(ncid,'dms',idv,tst)
    if ( tst ) iers(11) = -1
    call ccnf_inq_varid(ncid,'del_lat',idv,tst)
    if ( tst ) iers(12) = -1
    !call ccnf_inq_varid(ncid,'opldepth',idv,tst)
    !if ( tst ) iers(13) = -1
    call ccnf_inq_varid(ncid,'tsu',idv,tst)
    if ( tst ) then
      write(6,*) "ERROR: Cannot locate tsu in input file"
      call ccmpi_abort(-1)
    end if
  end if
  
  ! bcast data to all processors unless all processes are reading input files
  if ( .not.pfall ) then
    dumr(1:kk) = sigin(1:kk)
    dumr(kk+1:kk+12) = real(iers(1:12))
    if ( ok>0 ) dumr(kk+13:kk+ok+12) = gosig_1(1:ok)
    call ccmpi_bcast(dumr(1:kk+ok+12),0,comm_world)
    sigin(1:kk) = dumr(1:kk)
    iers(1:12) = nint(dumr(kk+1:kk+12))
    if ( ok>0 ) gosig_1(1:ok) = dumr(kk+13:kk+ok+12)
  end if
  
  mixr_found    = iers(1)==0
  siced_found   = iers(2)==0
  fracice_found = iers(3)==0
  soilt_found   = iers(4)==0
  mlo_found     = iers(5)==0
  mlo2_found    = iers(6)==0
  urban1_found  = iers(7)==0
  urban2_found  = iers(8)==0
  mloice_found  = iers(9)==0
  zht_found     = iers(10)==0
  aero_found    = iers(11)==0
  nllp_found    = iers(12)==0
  !mlo3_found    = iers(13)==0  
  
  ! determine whether zht needs to be read
  zht_needed = nested==0 .or. (nested==1.and.retopo_test/=0) .or.          &
      nested==3 .or. .not.(soilt_found.or.mlo_found)
  allowtrivialfill = zht_needed .and. .not.zht_found .and.                 &
      .not.(nested==0 .or. (nested==1.and.retopo_test/=0) .or. nested==3)
  if ( myid==0 ) then
    if ( zht_needed ) then
      write(6,*) "-> Surface height is required with zht_needed     =",zht_needed
    else  
      write(6,*) "-> Surface height is not required with zht_needed =",zht_needed
    end if
    write(6,'(A,2I4)') " -> nested,retopo_test                             =",nested,retopo_test
    write(6,*) "-> soilt_found,mlo_found,mlo2_found               =",soilt_found,mlo_found,mlo2_found
    write(6,*) "-> zht_found,mixr_found,aero_found,mloice_found   =",zht_found,mixr_found,aero_found,mloice_found
    write(6,*) "-> urban1_found,urban2_found,allowtrivialfill     =",urban1_found,urban2_found,allowtrivialfill
    write(6,*) "-> nllp_found                                     =",nllp_found
    if ( zht_needed .and. .not.zht_found .and. .not.allowtrivialfill ) then
      write(6,*) "ERROR: Surface height is required but not found in input file"
      call ccmpi_abort(-1)
    end if     
  end if  
      
  ! determine whether surface temperature needs to be interpolated (tss_test=.false.)
  tss_test = siced_found .and. fracice_found .and. iotest
  if ( myid==0 ) then
    if ( tss_test ) then
      write(6,*) "-> Surface temperature does not require interpolation"
      write(6,*) "-> tss_test,siced_found,fracice_found,iotest,ptest =",tss_test,siced_found,fracice_found,iotest,ptest
    else
      write(6,*) "-> Surface temperature requires interpolation"
      write(6,*) "-> tss_test,siced_found,fracice_found,iotest,ptest =",tss_test,siced_found,fracice_found,iotest,ptest
    end if
  end if
  
  ! read zht
  if ( allocated(zss_a) ) deallocate(zss_a)
  ! read zht for initial conditions or nudging or land-sea mask
  if ( zht_needed ) then
    if ( zht_found ) then  
      if ( tss_test ) then ! tss_test includes iotest=.true.
        allocate( zss_a(ifull) )
        call histrd(iarchi,ier,'zht',zss_a,ifull)
      else     
        allocate( zss_a(fwsize) )
        call histrd(iarchi,ier,'zht',zss_a,6*ik*ik)
        if ( fwsize>0 ) then
          nemi = 2  
          land_a = zss_a>0. ! 2nd guess for land-sea mask
          landlake_a = land_a
        end if
      end if
    else
      if ( tss_test ) then ! tss_test includes iotest=.true.
        allocate( zss_a(ifull) )
        zss_a = 0. ! ocean everywhere
      else     
        allocate( zss_a(fwsize) )
        if ( fwsize>0 ) then
          zss_a = 0.  ! ocean everywhere
          nemi = 2  
          land_a = zss_a>0. ! 2nd guess for land-sea mask
          landlake_a = land_a
        end if
      end if
    end if
  end if
  
  ! read soilt
  if ( soilt_found ) then
    ! read soilt for land-sea mask  
    if ( .not.tss_test ) then ! tss_test includes iotest=.true.
      call histrd(iarchi,ier,'soilt',ucc,6*ik*ik)
      if ( fwsize>0 ) then
        nemi = 3
        land_a = nint(ucc)>0 ! 1st guess for land-sea mask
        landlake_a = nint(ucc)/=0
      end if  
    end if
  end if  
  
  ! read host ocean bathymetry data
  if ( allocated(ocndep_a) ) deallocate( ocndep_a )
  if ( mlo_found ) then
    if ( tss_test ) then ! tss_test includes iotest=.true.
      allocate( ocndep_a(ifull) )
      call histrd(iarchi,ier,'ocndepth',ocndep_a,ifull)
    else     
      allocate( ocndep_a(fwsize) )
      call histrd(iarchi,ier,'ocndepth',ocndep_a,6*ik*ik)
      if ( fwsize>0 ) then
        if ( nemi==-1 ) then
          nemi = 2  
          land_a = ocndep_a<0.1 ! 3rd guess for land-sea mask
          landlake_a = land_a
        end if  
      end if
    end if
    ! calculate ocean half levels
    gosig_h(0) = 0.
    do k = 1,ok
      ! 0.5*(gosig_h(k-1)+gosig_h(k)) = gosig_1(k)  
      gosig_h(k) = 2.*gosig_1(k) - gosig_h(k-1) 
    end do  
    ! default 3d ocean depths
    do k = 1,ok
      gosig_3(:,k) = gosig_1(k)
    end do
  end if
  
  ! read partial depths
  !if ( allocated(opldep_a) ) deallocate( opldep_a )
  !if ( mlo3_found ) then
  !  if ( tss_test ) then ! tss_test includes iotest=.true.
  !    allocate( opldep_a(ifull) )
  !    call histrd(iarchi,ier,'opldepth',opldep_a,ifull)
  !  else
  !    allocate( opldep_a(fwsize) )
  !    call histrd(iarchi,ier,'opldepth',opldep_a,6*ik*ik)  
  !  end if    
  !end if    
  
  ! set-up 3d-depth, land_3d mask for z* ocean
  ! and set-up sea_a array
  if ( fwsize>0 ) then
    if ( mlo_found ) then
      ! land_3d mask
      if ( any(gosig_1>1.) ) then
        ! found z* ocean levels  
        land_3d(:,:) = .false.  
        do k = 1,ok
          ! include +0.1 to avoid rounding issues  
          land_3d(:,k) = ( land_a .or. gosig_1(k)+0.1>=ocndep_a ) 
        end do
      else
        ! found sigma ocean levels - to be depreciated
        do k = 1,ok
          land_3d(:,k) = land_a 
        end do
      end if  ! mlo_found ..else..      
      !if ( mlo3_found ) then
      ! allocate( gosig3_a(fwsize,ok) )  
      !  ! 3d-depth  
      !  do k = 1,ok
      !    gosig3_a(:,k) = gosig_1(:)
      !  end do
      !  do iq = 1,fwsize
      !    if ( opldep_a(iq)>1.e-4 ) then
      !      ! find ocean floor  
      !      ktest = 1  
      !      do k = 1,ok
      !        if ( gosig_h(k-1)<opldep_a(iq) ) then
      !          ktest = k
      !        else
      !          exit
      !        end if
      !      end do
      !      ! update 3d-depth with partial depth
      !      gosig3_a(iq,ktest) = opldep_a(iq)
      !    end if
      !  end do
      !  ! update land_3d mask
      !  if ( any(gosig_1>1.) ) then
      !    ! found z* ocean levels  
      !    land_3d = .false.  
      !    do k = 1,ok
      !      land_3d(:,k) = ( land_a .or. gosig3_a(:,k)+0.1>=ocndep_a ) 
      !    end do
      !  end if
      !  deallocate( gosig3_a )
      !end if ! mlo3_found  
    end if   ! mlo_found
    sea_a = .not.land_a
    do k = 1,ok
      landlake_3d(:,k) = land_3d(:,k) .or. landlake_a
    end do  
  end if ! fwsize>0
  
  ! check that land-sea mask is definied
  if ( fwsize>0 .and. .not.tss_test ) then ! tss_test includes iotest=.true.
    if ( nemi==-1 ) then
      write(6,*) "ERROR: Cannot determine land-sea mask"
      write(6,*) "CCAM requires zht or soilt or ocndepth in input file"
      call ccmpi_abort(-1)
    end if  
    if ( myid==0 .and. nmaxpr==1 ) then
      write(6,*) "-> Land-sea mask using nemi = ",nemi
    end if
  end if  
  
  ! read urban data mask
  ! read urban mask for urban and initial conditions and interpolation
  if ( nurban/=0 .and. nested/=1 .and. nested/=3 .and. .not.iotest .and. nhstest>=0 ) then
    if ( myid==0 .and. nmaxpr==1 ) then
      write(6,*) "-> Determine urban mask"
    end if  
    if ( urban1_found ) then
      call histrd(iarchi,ier,'t1_intmtgg1',ucc,6*ik*ik)  
    else if ( urban2_found ) then
      call histrd(iarchi,ier,'intmtgg1',ucc,6*ik*ik)
    else
      if ( fwsize>0 ) then
        ! will use tsu for urban temperatures so all points are valid  
        ucc = 0.
      end if
    end if  
    if ( fwsize>0 ) then
      nourban_a = ucc>=399.
    end if
  end if  
    
else
    
  ! use saved metadata  
  mixr_found    = iers(1)==0
  siced_found   = iers(2)==0
  fracice_found = iers(3)==0
  soilt_found   = iers(4)==0
  mlo_found     = iers(5)==0
  mlo2_found    = iers(6)==0
  urban1_found  = iers(7)==0
  urban2_found  = iers(8)==0
  mloice_found  = iers(9)==0
  zht_found     = iers(10)==0
  aero_found    = iers(11)==0
  nllp_found    = iers(12)==0
  !mlo3_found    = iers(13)==0
  zht_needed    = nested==0 .or. (nested==1.and.retopo_test/=0) .or.      &
      nested==3 .or. .not.(soilt_found.or.mlo_found)
  allowtrivialfill = zht_needed .and. .not.zht_found .and.                &
      .not.(nested==0 .or. (nested==1.and.retopo_test/=0) .or. nested==3)
  tss_test      = siced_found .and. fracice_found .and. iotest
  
end if ! newfile ..else..

! -------------------------------------------------------------------
! detemine the reference level below sig=0.9 (used to calculate psl)
levk = 0
levkin = 0
if ( nested==0 .or. (nested==1.and.retopo_test/=0) .or. nested==3 ) then
  do while( sig(levk+1)>0.9 ) ! nested grid
    levk = levk + 1
  end do
  do while( sigin(levkin+1)>0.9 ) ! host grid
    levkin = levkin + 1
  end do
  if ( levkin==0 ) then
    write(6,*) "ERROR: Invalid sigma levels in input file"
    write(6,*) "sigin = ",sigin
    call ccmpi_abort(-1)
  end if
end if

if ( myid==0 ) write(6,*) "Reading prognostic fields"

!--------------------------------------------------------------------
! Read surface pressure
! psf read when nested=0 or nested=1.and.retopo_test/=0 or nested=3
psl(1:ifull) = 0.
if ( nested==0 .or. (nested==1.and.retopo_test/=0) .or. nested==3 ) then
  if ( iotest ) then
    call histrd(iarchi,ier,'psf',psl,ifull)
    if ( minval(psl(1:ifull))<-1.6 .or. maxval(psl(1:ifull))>0.4 ) then
      write(6,*) "ERROR: Input psf is out of range"
      write(6,*) "minval,maxval ",minval(psl(1:ifull)),maxval(psl(1:ifull))
      call ccmpi_abort(-1)
    end if
  else
    allocate( psl_a(fwsize) )
    psl_a(:) = 0.
    call histrd(iarchi,ier,'psf',psl_a,6*ik*ik)
    if ( minval(psl_a)<-1.6 .or. maxval(psl_a)>0.4 ) then
      write(6,*) "ERROR: Input psf is out of range"
      write(6,*) "minval,maxval ",minval(psl_a),maxval(psl_a)
      call ccmpi_abort(-1)
    end if    
  end if
endif

! -------------------------------------------------------------------
! Read surface temperature 
! read global tss to diagnose sea-ice or land-sea mask
if ( tss_test ) then ! tss_test includes iotest=.true.
  call histrd(iarchi,ier,'tsu',tss,ifull)
  tss = abs(tss)
  tss = min( max( tss, 100. ), 425. )
else
  call histrd(iarchi,ier,'tsu',tss_a,6*ik*ik)
  tss_a = abs(tss_a)
  tss_a = min( max( tss_a, 100. ), 425. )
end if ! (tss_test) ..else..

 
!--------------------------------------------------------------
! Read ocean data for nudging (sea-ice is read below)
! read when nested=0 or nested=1.and.nud/=0 or nested=2 .or nested=4
if ( abs(nmlo)>=1 .and. abs(nmlo)<=9 .and. nested/=3 ) then
  ! defalt values for surface height
  ocndwn(1:ifull,2) = 0.
  if ( mlo_found ) then
    ! read water surface height
    if ( (nested/=1.or.nud_sfh/=0) .and. ok>0 ) then
      call fillhist1('ocheight',ocndwn(:,2),land_a,fill_land)
      where ( land(1:ifull) )
        ocndwn(1:ifull,2) = 0.
      end where  
    end if
  end if
end if
!--------------------------------------------------------------

!--------------------------------------------------------------
! read sea ice here for prescribed SSTs configuration and for
! mixed-layer-ocean
if ( tss_test ) then ! tss_test includes iotest=.true.

  call histrd(iarchi,ier,'siced',sicedep,ifull)
  call histrd(iarchi,ier,'fracice',fracice,ifull)
  ! check for invalid sea-ice data
  if ( any(fracice(1:ifull)>1.1) ) then
    write(6,*) "ERROR: Invalid fracice in input file"
    write(6,*) "Fracice should be between 0 and 1"
    write(6,*) "maximum fracice ",maxval(fracice(1:ifull))
    call ccmpi_abort(-1)
  end if
  ! fix rounding errors
  fracice(1:ifull) = min( fracice(1:ifull), 1. )
  ! update surface height and ocean depth if required
  if ( zht_needed ) then
    zss(1:ifull) = zss_a(1:ifull) ! use saved zss arrays
  end if
  if ( abs(nmlo)>0 .and. abs(nmlo)<=9 ) then
    if ( mlo_found ) then
      ocndwn(1:ifull,1) = ocndep_a(1:ifull)
      !opldep(1:ifull) = 0.      
    end if
    !if ( mlo3_found ) then
    !  opldep(1:ifull) = opldep_a(1:ifull)
    !end if  
  end if    

else

  allocate( fracice_a(fwsize), sicedep_a(fwsize) )  
  allocate( tss_l_a(fwsize), tss_s_a(fwsize) )
    
  call histrd(iarchi,ier,'siced',sicedep_a,6*ik*ik)
  call histrd(iarchi,ier,'fracice',fracice_a,6*ik*ik)
        
  if ( fwsize>0 ) then

    ! diagnose sea-ice if required      
    if ( any(fracice_a(1:fwsize)>1.1) ) then
      write(6,*) "ERROR: Invalid fracice in input file"
      write(6,*) "Fracice should be between 0 and 1"
      write(6,*) "maximum fracice ",maxval(fracice_a(1:fwsize))
      call ccmpi_abort(-1)
    end if
    ! fix rounding errors
    fracice_a(1:fwsize) = min( fracice_a(1:fwsize), 1. )

    if ( siced_found .and. fracice_found ) then
      ! do nothing
    else if ( siced_found ) then       ! i.e. sicedep read in, fracice not read in
      where ( sicedep_a(1:fwsize)>0. )
        fracice_a(1:fwsize) = 1.
      end where
    else if ( fracice_found ) then     ! i.e. only fracice read in;  sicedep not read in
      where ( fracice_a(1:fwsize)>0.01 )
        sicedep_a(1:fwsize) = 2.
      elsewhere
        sicedep_a(1:fwsize) = 0.
        fracice_a(1:fwsize) = 0.
      end where
    else
      ! neither sicedep nor fracice read in
      sicedep_a(1:fwsize) = 0.  ! Oct 08
      fracice_a(1:fwsize) = 0.
      if ( myid==0 .and. nmaxpr==1 ) write(6,*) 'pre-setting siced in onthefly from tss'
      where ( abs(tss_a(1:fwsize))<=271.6 ) ! for ERA-Interim
        sicedep_a(1:fwsize) = 1.  ! Oct 08  ! previously 271.2
        fracice_a(1:fwsize) = 1.
      end where
    end if  ! siced_found .and. fracice_found ..else..

    ! fill surface temperature and sea-ice
    tss_l_a(1:fwsize) = abs(tss_a(1:fwsize))
    tss_s_a(1:fwsize) = abs(tss_a(1:fwsize))
    call fill_cc1(tss_l_a,sea_a,fill_sea)
    !if ( mlo3_found ) then
    !  ucc6(:,1) = tss_s_a
    !  ucc6(:,2) = sicedep_a
    !  ucc6(:,3) = fracice_a        
    !  ucc6(:,4) = ocndep_a
    !  ucc6(:,5) = opldep_a
    !  call fill_cc4(ucc6(:,1:5),land_a,fill_land)
    ! tss_s_a   = ucc6(:,1)
    !  sicedep_a = ucc6(:,2)
    !  fracice_a = ucc6(:,3)
    !  ocndep_a = ucc6(:,4)
    !  opldep_a = ucc6(:,5)
    ! else ...
    if ( mlo_found ) then
      ucc6(:,1) = tss_s_a
      ucc6(:,2) = sicedep_a
      ucc6(:,3) = fracice_a        
      ucc6(:,4) = ocndep_a
      call fill_cc4(ucc6(:,1:4),land_a,fill_land)
      tss_s_a   = ucc6(:,1)
      sicedep_a = ucc6(:,2)
      fracice_a = ucc6(:,3)
      ocndep_a = ucc6(:,4)
    else    
      ucc6(:,1) = tss_s_a
      ucc6(:,2) = sicedep_a
      ucc6(:,3) = fracice_a        
      call fill_cc4(ucc6(:,1:3),land_a,fill_land)
      tss_s_a   = ucc6(:,1)
      sicedep_a = ucc6(:,2)
      fracice_a = ucc6(:,3)
    end if  
  end if ! fwsize>0

  if ( fwsize>0 ) then
    ucc6(:,1:6) = 0.
    if ( zht_needed ) ucc6(:,1) = zss_a
    if ( mlo_found )  ucc6(:,2) = ocndep_a
    ucc6(:,3) = tss_l_a
    ucc6(:,4) = tss_s_a
    ucc6(:,5) = sicedep_a
    ucc6(:,6) = fracice_a
    !if ( mlo3_found ) ucc6(:,7) = opldep_a
  end if          
  call doints4(ucc6(:,1:6),udum6(:,1:6))
  zss      = udum6(:,1)
  if ( abs(nmlo)>0 .and. abs(nmlo)<=9 ) then
    ocndwn(1:ifull,1) = udum6(1:ifull,2)
  end if  
  tss_l    = udum6(:,3)
  tss_s    = udum6(:,4)
  sicedep  = udum6(:,5)
  fracice  = udum6(:,6)
  !opldep   = udum6(:,7)

  !   incorporate other target land mask effects
  where ( land(1:ifull) )
    sicedep(1:ifull) = 0.
    fracice(1:ifull) = 0.
    tss(1:ifull) = tss_l(1:ifull)
  elsewhere
    tss(1:ifull) = tss_s(1:ifull)
  end where
  where ( sicedep(1:ifull)<0.05 )
    sicedep(1:ifull) = 0.
    fracice(1:ifull) = 0.
  end where
  
  if ( any(tss(1:ifull)>900.) ) then
    write(6,*) "ERROR: Unable to interpolate surface temperature"
    write(6,*) "Possible problem with land-sea mask in input file"
    call ccmpi_abort(-1)
  end if

  deallocate( fracice_a, sicedep_a )
  deallocate( tss_l_a, tss_s_a )
  
end if ! (tss_test ) ..else..


! update 3d ocean depths with partial step depths
!if ( newfile ) then
!  if ( mlo3_found ) then
!    do iq = 1,ifull
!      if ( opldep(iq)>1.e-4 ) then        
!        ktest = 1  
!        do k = 1,ok
!          if ( gosig_h(k-1)<opldep(iq) ) then
!            ktest = k
!          else
!            exit  
!          end if
!        end do
!        gosig_3(iq,ktest) = opldep(iq)
!      end if
!    end do
!  end if ! mlo3_found
!end if   ! newfile


! -------------------------------------------------------------------
! read atmospheric fields for nested=0 or nested=1.and.nud/=0 or nested=3

! air temperature
! read for nested=0 or nested=1.and.retopo_test/=0 or nested=3
if ( nested==0 .or. (nested==1.and.retopo_test/=0) .or. nested==3 ) then
  allocate( t_a_lev(fwsize) )  
  call gethist4a('temp',t,1,levkin=levkin,t_a_lev=t_a_lev)
else
  t(1:ifull,1:kl) = 300.    
end if ! (nested==0.or.(nested==1.and.retopo_test/=0).or.nested==3)

! winds
! read for nested=0 or nested=1.and.nud_uv/=0 or nested=3
if ( nested==0 .or. (nested==1.and.nud_uv/=0) .or. nested==3 ) then
  call gethistuv4a('u','v',u,v,3,4)
  if ( any(u(1:ifull,1:kl)/=u(1:ifull,1:kl)) ) then
    write(6,*) "ERROR: Invalid NaN u-wind in onthefly after gethistuv4a"
    call ccmpi_abort(-1)
  end if  
  if ( any(abs(u(1:ifull,1:kl))>350.) ) then
    write(6,*) "ERROR: Invalid high u-wind in onthefly after gethistuv4a"
    write(6,*) "u max,min ",maxval(u(1:ifull,1:kl)),minval(u(1:ifull,1:kl))
    call ccmpi_abort(-1)
  end if  
  if ( any(v(1:ifull,1:kl)/=v(1:ifull,1:kl)) ) then
    write(6,*) "ERROR: Invalid NaN v-wind in onthefly after gethistuv4a"
    call ccmpi_abort(-1)
  end if  
  if ( any(abs(v(1:ifull,1:kl))>350.) ) then
    write(6,*) "ERROR: Invalid high v-wind in onthefly after gethistuv4a"
    write(6,*) "v max,min ",maxval(v(1:ifull,1:kl)),minval(v(1:ifull,1:kl))
    call ccmpi_abort(-1)
  end if  
else
  u(1:ifull,1:kl) = 0.
  v(1:ifull,1:kl) = 0.
end if ! (nested==0.or.(nested==1.and.nud_uv/=0).or.nested==3)

! mixing ratio
! read for nested=0 or nested=1.and.nud_q/=0 or nested=3
if ( nested==0 .or. (nested==1.and.nud_q/=0) .or. nested==3 ) then
  if ( mixr_found ) then
    call gethist4a('mixr',qg,2)      !     mixing ratio
  else
    call gethist4a('q',qg,2)         !     mixing ratio
  end if
  qg(1:ifull,1:kl) = max( qg(1:ifull,1:kl), 0. )
else
  qg(1:ifull,1:kl) = qgmin
end if ! (nested==0.or.(nested==1.and.nud_q/=0).or.nested==3)

if ( abs(nmlo)>=1 .and. abs(nmlo)<=9 .and. nested/=3 ) then
  ! defalt values
  do k = 1,wlev
    call mloexpdep(0,depth,k,0)
    ! This polynomial fit is from MOM3, based on Levitus
    where (depth(:)<2000.)
      mlodwn(1:ifull,k,1) = 18.4231944       &
        - 0.43030662E-1*depth(1:ifull)       &
        + 0.607121504E-4*depth(1:ifull)**2   &
        - 0.523806281E-7*depth(1:ifull)**3   &
        + 0.272989082E-10*depth(1:ifull)**4  &
        - 0.833224666E-14*depth(1:ifull)**5  &
        + 0.136974583E-17*depth(1:ifull)**6  &
        - 0.935923382E-22*depth(1:ifull)**7
      mlodwn(1:ifull,k,1) = mlodwn(1:ifull,k,1) - wrtemp + tss(1:ifull) - 18.4231944
      mlodwn(1:ifull,k,1) = max( mlodwn(1:ifull,k,1), 275.16-wrtemp )
    elsewhere
      mlodwn(1:ifull,k,1) = 275.16 - wrtemp
    end where
    where ( isoilm_in == 0 )
      mlodwn(1:ifull,k,2) = 34.72 ! sal
    elsewhere
      mlodwn(1:ifull,k,2) = 0.    ! sal (freshwater)  
    end where    
    mlodwn(1:ifull,k,3) = 0.      ! uoc
    mlodwn(1:ifull,k,4) = 0.      ! voc
    mlodwn(1:ifull,k,5) = omink   ! tke
    mlodwn(1:ifull,k,6) = omineps ! eps
  end do  
  if ( mlo_found ) then
    ! ocean potential temperature
    ! ocean temperature and soil temperature use the same arrays
    ! as no fractional land or sea cover is allowed in CCAM
    if ( (nested/=1.or.nud_sst/=0) .and. ok>0 ) then
      if ( mlo2_found ) then
        call fillhist4o('thetao',mlodwn(:,:,1),land_3d,fill_floor,ocndwn(:,1))
      else
        call fillhist4o('tgg',mlodwn(:,:,1),land_3d,fill_floor,ocndwn(:,1))
      end if  
      where ( mlodwn(1:ifull,1:wlev,1)>100. )
        mlodwn(1:ifull,1:wlev,1) = mlodwn(1:ifull,1:wlev,1) - wrtemp ! remove temperature offset for precision
      end where
    end if ! (nestesd/=1.or.nud_sst/=0) ..else..
    ! ocean salinity
    if ( (nested/=1.or.nud_sss/=0) .and. ok>0 ) then
      if ( mlo2_found ) then
        call fillhist4o('so',mlodwn(:,:,2),land_3d,fill_floorlake,ocndwn(:,1))  
      else    
        call fillhist4o('sal',mlodwn(:,:,2),land_3d,fill_floorlake,ocndwn(:,1))
      end if  
      do k = 1,wlev
        where ( isoilm_in == -1 )
          mlodwn(1:ifull,k,2) = 0. ! freshwater lakes
        end where
      end do  
      mlodwn(1:ifull,1:wlev,2) = max( mlodwn(1:ifull,1:wlev,2), 0. )
    end if ! (nestesd/=1.or.nud_sss/=0) ..else..
    ! ocean currents
    if ( (nested/=1.or.nud_ouv/=0) .and. ok>0 ) then
      if ( mlo2_found ) then
        call fillhistuv4o('uo','vo',mlodwn(:,:,3),mlodwn(:,:,4),land_3d,fill_floor,ocndwn(:,1))  
      else    
        call fillhistuv4o('uoc','voc',mlodwn(:,:,3),mlodwn(:,:,4),land_3d,fill_floor,ocndwn(:,1))
      end if  
    end if ! (nestesd/=1.or.nud_ouv/=0) ..else..
  end if   ! mlo_found
end if     ! abs(nmlo)>=1 .and. abs(nmlo)<=9 .and. nested/=3

!------------------------------------------------------------
! Aerosol data
if ( abs(iaero)>=2 .and. nested/=3 ) then
  if ( nested/=1 .or. nud_aero/=0 ) then
    if ( aero_found ) then  
      call gethist4a('dms',  xtgdwn(:,:,1), 5)
      call gethist4a('so2',  xtgdwn(:,:,2), 5)
      call gethist4a('so4',  xtgdwn(:,:,3), 5)
      call gethist4a('bco',  xtgdwn(:,:,4), 5)
      call gethist4a('bci',  xtgdwn(:,:,5), 5)
      call gethist4a('oco',  xtgdwn(:,:,6), 5)
      call gethist4a('oci',  xtgdwn(:,:,7), 5)
      call gethist4a('dust1',xtgdwn(:,:,8), 5)
      call gethist4a('dust2',xtgdwn(:,:,9), 5)
      call gethist4a('dust3',xtgdwn(:,:,10),5)
      call gethist4a('dust4',xtgdwn(:,:,11),5)
      call gethist4a('salt1',xtgdwn(:,:,12),5)
      call gethist4a('salt2',xtgdwn(:,:,13),5)
      do i = 1,13
        if ( any(xtgdwn(:,:,i)>aerosol_tol) ) then
          write(6,*) "ERROR: Bad aerosol data in host"
          write(6,*) "Index ",i
          write(6,*) "Maxval ",maxval(xtgdwn(:,:,i))
          call ccmpi_abort(-1)
        end if  
      end do  
      xtgdwn(:,:,:) = max( xtgdwn(:,:,:), 0. )
    else
      xtgdwn(:,:,:) = 0.  
    end if    
  end if  
end if

!------------------------------------------------------------
! re-grid surface pressure by mapping to MSLP, interpolating and then map to surface pressure
! requires psl_a, zss, zss_a, t and t_a_lev
if ( nested==0 .or. (nested==1.and.retopo_test/=0) .or. nested==3 ) then
  if ( .not.iotest ) then
    if ( fwsize>0 ) then
      ! ucc holds pmsl_a
      call mslpx(ucc,psl_a,zss_a,t_a_lev,sigin(levkin))  ! needs pmsl (preferred)
    end if
    call doints1(ucc,pmsl)
    ! invert pmsl to get psl
    call to_pslx(pmsl,psl,zss,t(:,levk),sig(levk))  ! on target grid
    deallocate( psl_a )
  end if ! .not.iotest
  deallocate( t_a_lev )
end if


if ( abs(iaero)>=2 .and. nested/=3 ) then
  if ( nested/=1.or.nud_aero/=0 ) then
    ! Factor 1.e3 to convert to g/m2, x 3 to get sulfate from sulfur
    so4t(:) = 0.
    do k = 1,kl
      so4t(1:ifull) = so4t(1:ifull) + 3.e3*xtgdwn(1:ifull,k,3)*(-1.e5*exp(psl(1:ifull))*dsig(k))/grav
    end do
  end if  
end if


!**************************************************************
! This is the end of reading the nudging arrays
!**************************************************************


!--------------------------------------------------------------
! The following data is only read for initial conditions
if ( nested/=1 .and. nested/=3 ) then

  ierc(:) = 0  ! flag for located variables
    
  !------------------------------------------------------------------
  ! check soil variables
  if ( myid==0 .or. pfall ) then
    if ( ccycle/=0 ) then
      call ccnf_inq_varid(ncid,'nplant1',idv,tst)
      if ( .not.tst ) ierc(9) = 1
      call ccnf_inq_varid(ncid,'p01_nplant1',idv,tst)
      if ( .not.tst ) ierc(10) = 1
    end if ! ccycle/=0
    do k = 1,ms
      write(vname,'("tgg",I1.1)') k
      call ccnf_inq_varid(ncid,vname,idv,tst)
      if ( .not.tst ) ierc(11+k) = 1
      write(vname,'("wetfrac",I1.1)') k
      call ccnf_inq_varid(ncid,vname,idv,tst)
      if ( .not.tst ) ierc(11+ms+k) = 1
      write(vname,'("wb",I1.1)') k
      call ccnf_inq_varid(ncid,vname,idv,tst)
      if ( .not.tst ) ierc(11+2*ms+k) = 1
    end do
  end if
  
  ! -----------------------------------------------------------------
  ! verify if input is a restart file
  if ( nested==0 ) then
    if ( myid==0 .or. pfall ) then
      if ( kk==kl .and. iotest ) then
        lrestart = .true.
        call ccnf_inq_varid(ncid,'dpsldt',idv,tst)
        if ( tst ) lrestart = .false.
        call ccnf_inq_varid(ncid,'zgnhs',idv,tst)
        if ( tst ) lrestart = .false.
        call ccnf_inq_varid(ncid,'sdot',idv,tst)
        if ( tst ) lrestart = .false.
        call ccnf_inq_varid(ncid,'pslx',idv,tst)
        if ( tst ) lrestart = .false.
        call ccnf_inq_varid(ncid,'savu',idv,tst)
        if ( tst ) lrestart = .false.
        call ccnf_inq_varid(ncid,'savv',idv,tst)
        if ( tst ) lrestart = .false.
        call ccnf_inq_varid(ncid,'savu1',idv,tst)
        if ( tst ) lrestart = .false.
        call ccnf_inq_varid(ncid,'savv1',idv,tst)
        if ( tst ) lrestart = .false.
        call ccnf_inq_varid(ncid,'savu2',idv,tst)
        if ( tst ) lrestart = .false.
        call ccnf_inq_varid(ncid,'savv2',idv,tst)
        if ( tst ) lrestart = .false.
        call ccnf_inq_varid(ncid,'nstag',idv,tst)
        if ( tst ) then
          lrestart = .false.
        else 
          call ccnf_get_vara(ncid,idv,iarchi,ierc(5))
        end if
        call ccnf_inq_varid(ncid,'nstagu',idv,tst)
        if ( tst ) then
          lrestart = .false.
        else 
          call ccnf_get_vara(ncid,idv,iarchi,ierc(6))
        end if
        call ccnf_inq_varid(ncid,'nstagoff',idv,tst)
        if ( tst ) then
          lrestart = .false.
        else 
          call ccnf_get_vara(ncid,idv,iarchi,ierc(7))
        end if
        if ( abs(nmlo)>=3 .and. abs(nmlo)<=9 ) then
          if ( ok==wlev ) then
            call ccnf_inq_varid(ncid,'old1_uo',idv,tst)
            if ( tst ) lrestart = .false.
            call ccnf_inq_varid(ncid,'old1_vo',idv,tst)
            if ( tst ) lrestart = .false.
            call ccnf_inq_varid(ncid,'old2_uo',idv,tst)
            if ( tst ) lrestart = .false.
            call ccnf_inq_varid(ncid,'old2_vo',idv,tst)
            if ( tst ) lrestart = .false.                
            call ccnf_inq_varid(ncid,'ipice',idv,tst)
            if ( tst ) lrestart = .false.
            call ccnf_inq_varid(ncid,'nstagoffmlo',idv,tst)
            if ( tst ) then
              lrestart = .false.
            else
              call ccnf_get_vara(ncid,idv,iarchi,ierc(8))
            end if
          else
            lrestart = .false.
          end if
        end if
        if ( nmlo/=0 .and. abs(nmlo)<=9 ) then
          if ( ok==wlev ) then
            call ccnf_inq_varid(ncid,'old1_uotop',idv,tst)
            if ( tst ) lrestart = .false.
            call ccnf_inq_varid(ncid,'old1_votop',idv,tst)
            if ( tst ) lrestart = .false.
            call ccnf_inq_varid(ncid,'old1_uobot',idv,tst)
            if ( tst ) lrestart = .false.
            call ccnf_inq_varid(ncid,'old1_vobot',idv,tst)
            if ( tst ) lrestart = .false.
          else
            lrestart = .false.
          end if
        end if
        lrestart_radiation = .true.
        call ccnf_inq_varid(ncid,'sgsave',idv,tst)
        if ( tst ) lrestart_radiation = .false. 
        call ccnf_inq_varid(ncid,'rgsave',idv,tst)
        if ( tst ) lrestart_radiation = .false.
        call ccnf_inq_varid(ncid,'rtu',idv,tst)
        if ( tst ) lrestart_radiation = .false. 
        call ccnf_inq_varid(ncid,'rtc',idv,tst)
        if ( tst ) lrestart_radiation = .false.
        call ccnf_inq_varid(ncid,'rgdn',idv,tst)
        if ( tst ) lrestart_radiation = .false. 
        call ccnf_inq_varid(ncid,'rgn',idv,tst)
        if ( tst ) lrestart_radiation = .false.
        call ccnf_inq_varid(ncid,'rgc',idv,tst)
        if ( tst ) lrestart_radiation = .false. 
        call ccnf_inq_varid(ncid,'sint_amp',idv,tst)
        if ( tst ) lrestart_radiation = .false.
        call ccnf_inq_varid(ncid,'sout_amp',idv,tst)
        if ( tst ) lrestart_radiation = .false. 
        call ccnf_inq_varid(ncid,'soutclr_amp',idv,tst)
        if ( tst ) lrestart_radiation = .false.
        call ccnf_inq_varid(ncid,'sgdn_amp',idv,tst)
        if ( tst ) lrestart_radiation = .false. 
        call ccnf_inq_varid(ncid,'sgdndir_amp',idv,tst)
        if ( tst ) lrestart_radiation = .false.
        call ccnf_inq_varid(ncid,'sgn_amp',idv,tst)
        if ( tst ) lrestart_radiation = .false.
        call ccnf_inq_varid(ncid,'sgclr_amp',idv,tst)
        if ( tst ) lrestart_radiation = .false.
        call ccnf_inq_varid(ncid,'sgdclr_amp',idv,tst)
        if ( tst ) lrestart_radiation = .false.
        call ccnf_inq_varid(ncid,'dni_amp',idv,tst)
        if ( tst ) lrestart_radiation = .false.
        call ccnf_inq_varid(ncid,'fbeamvis',idv,tst)
        if ( tst ) lrestart_radiation = .false.
        call ccnf_inq_varid(ncid,'fbeamnir',idv,tst)
        if ( tst ) lrestart_radiation = .false. 
        call ccnf_inq_varid(ncid,'swrsave',idv,tst)
        if ( tst ) lrestart_radiation = .false.
        call ccnf_inq_varid(ncid,'cloudlo',idv,tst)
        if ( tst ) lrestart_radiation = .false.
        call ccnf_inq_varid(ncid,'cloudmi',idv,tst)
        if ( tst ) lrestart_radiation = .false. 
        call ccnf_inq_varid(ncid,'cloudhi',idv,tst)
        if ( tst ) lrestart_radiation = .false.
        call ccnf_inq_varid(ncid,'sw_tend_amp',idv,tst)
        if ( tst ) lrestart_radiation = .false. 
        call ccnf_inq_varid(ncid,'lw_tend',idv,tst)
        if ( tst ) lrestart_radiation = .false.
        lrestart_convection = .true.
        call ccnf_inq_varid(ncid,'dmsedt_rad',idv,tst)
        if ( tst ) lrestart_convection = .false.
        call ccnf_inq_varid(ncid,'dmsedt_pbl',idv,tst)
        if ( tst ) lrestart_convection = .false.
        lrestart_tracer = .true.
        call ccnf_inq_varid(ncid,'tr0001',idv,tst)
        if ( tst ) lrestart_tracer = .false.
      else
        lrestart = .false.
        lrestart_radiation = .false.
        lrestart_convection = .false.
        lrestart_tracer = .false.
      end if ! kk=kl .and. iotest
      ierc(1:5) = 0
      if ( lrestart ) ierc(1) = 1
      if ( lrestart_radiation ) ierc(2) = 1
      if ( lrestart_tracer ) ierc(3) = 1
      call ccnf_inq_varid(ncid,'u10',idv,tst)
      if ( .not.tst ) ierc(4) = 1
      if ( lrestart_convection ) ierc(11) = 1
    end if ! myid==0 .or. pfall
  end if   ! nested==0  
    
  if ( .not.pfall ) then
    call ccmpi_bcast(ierc(1:11+3*ms),0,comm_world)
  end if
  
  lrestart  = (ierc(1)==1)
  lrestart_radiation = (ierc(2)==1)
  lrestart_tracer = (ierc(3)==1)  
  u10_found = (ierc(4)==1)
  lrestart_convection = (ierc(11)==1)
  if ( lrestart ) then
    nstag       = ierc(5)
    nstagu      = ierc(6)
    nstagoff    = ierc(7)
    nstagoffmlo = ierc(8)
    if ( myid==0 .and. nmaxpr==1 ) then
      write(6,*) "Continue staggering from"
      write(6,*) "nstag,nstagu,nstagoff ",nstag,nstagu,nstagoff
      if ( abs(nmlo)>=3 .and. abs(nmlo)<=9 ) then
        write(6,*) "nstagoffmlo ",nstagoffmlo
      end if
    end if
  end if
  carbon_found        = ierc(9)==1
  carbon2_found       = ierc(10)==1
  tgg_found(1:ms)     = ierc(12:11+ms)==1
  wetfrac_found(1:ms) = ierc(12+ms:11+2*ms)==1
  wb_found(1:ms)      = ierc(12+2*ms:11+3*ms)==1
  
  if ( myid==0 ) then
    write(6,*) "-> Found lrestart,lrestart_radiation,lrestart_tracer = ",lrestart,lrestart_radiation,lrestart_tracer
    write(6,*) "->       lrestart_convection,u10_found               = ",lrestart_convection,u10_found
    write(6,*) "->       carbon_found,carbon2_found                  = ",carbon_found,carbon2_found
    write(6,*) "->       tgg_found,wetfrac_found,wb_found            = ",any(tgg_found),any(wetfrac_found),any(wb_found)
  end if

  
  !------------------------------------------------------------------
  ! Read basic fields to be avaliable at zero time-step
  ! These values might need to be written at the zeroth time-step or
  ! used before the physics has updated. We then need to load this
  ! data from the restart file.
  if ( nested==0 .and. lrestart ) then
    if ( nsib==6.or.nsib==7 .and. nhstest>=0 ) then
      call gethist1('rs',rsmin)  
    end if
    call gethist1('zolnd',zo)    
    call gethist1('rnd',duma)
    precip(:) = real(duma/real(nperday),8)
    call gethist1('rnc',duma)
    precc(:) = real(duma/real(nperday),8)
    call gethist1('cll',duma)    
    cll_ave(:) = real(duma,8)
    call gethist1('clm',duma)    
    clm_ave(:) = real(duma,8)
    call gethist1('clh',duma)    
    clh_ave(:) = real(duma,8)
    call gethist1('cld',duma)
    cld_ave(:) = real(duma,8)
    call gethist1('sgdn_ave',duma)
    sgdn_ave(:) = real(duma,8)
    call gethist1('sgdc_ave',duma)
    sgdc_ave(:) = real(duma,8)
  end if
  
  !------------------------------------------------------------------
  ! Read snow and soil tempertaure
  call gethist1('snd',snowd)
  where ( .not.land(1:ifull) .and. (sicedep(1:ifull)<1.e-20 .or. nmlo==0) )
    snowd(1:ifull) = 0.
  end where
  if ( nested/=4 ) then
    if ( all(tgg_found(1:ms)) ) then
      call fillhist4('tgg',tgg,sea_a,fill_sea)
    else
      do k = 1,ms 
        if ( tgg_found(k) ) then
          write(vname,'("tgg",I1.1)') k
        else if ( k<=3 .and. tgg_found(2) ) then
          vname="tgg2"
        else if ( k<=3 ) then
          vname="tb3"
        else if ( tgg_found(6) ) then
          vname="tgg6"
        else
          vname="tb2"
        end if
        if ( iotest ) then
          if ( k==1 .and. .not.tgg_found(1) ) then
            tgg(1:ifull,k) = tss(1:ifull)
          else
            call histrd(iarchi,ier,vname,tgg(:,k),ifull)
          end if
        else
          if ( k==1 .and. .not.tgg_found(1) ) then
            ucc(1:fwsize) = tss_a(1:fwsize)
          else
            call histrd(iarchi,ier,vname,ucc,6*ik*ik)
          end if
          call fill_cc1(ucc,sea_a,fill_sea)
          call doints1(ucc,tgg(:,k))
        end if
      end do
    end if
    do k = 1,ms
      where ( tgg(1:ifull,k)<100. )
        tgg(1:ifull,k) = tgg(1:ifull,k) + wrtemp ! adjust range of soil temp for compressed history file
      end where
    end do  
    if ( .not.iotest ) then
      where ( snowd(1:ifull)>0. .and. land(1:ifull) )
        tgg(1:ifull,1) = min( tgg(1:ifull,1), 270.1 )
      endwhere
      do k = 1,ms
        tgg(1:ifull,k) = max( min( tgg(1:ifull,k), 400. ), 200. )
      end do  
    end if
  end if ! nested/=4  

  !--------------------------------------------------
  ! Read MLO sea-ice data
  if ( abs(nmlo)>=1 .and. abs(nmlo)<=9 ) then
    if ( .not.allocated(micdwn) ) allocate( micdwn(ifull,10) )
    micdwn(1:ifull,1) = 270.
    micdwn(1:ifull,2) = 270.
    micdwn(1:ifull,3) = 270.
    micdwn(1:ifull,4) = 270.
    micdwn(1:ifull,5) = fracice(1:ifull) ! read above with nudging arrays
    micdwn(1:ifull,6) = sicedep(1:ifull) ! read above with nudging arrays
    micdwn(1:ifull,7) = snowd(1:ifull)*1.e-3
    micdwn(1:ifull,8) = 0.  ! sto
    micdwn(1:ifull,9) = 0.  ! uic
    micdwn(1:ifull,10) = 0. ! vic
    if ( mloice_found ) then
      call fillhist4('tggsn',micdwn(:,1:4),land_a,fill_land)
      call fillhist1('sto',micdwn(:,8),land_a,fill_land)
      call fillhistuv1o('uic','vic',micdwn(:,9),micdwn(:,10),land_a,fill_land)
    end if
  end if
  
  !------------------------------------------------------------------
  ! Read river data
  call gethist1('swater',watbdy)

  !------------------------------------------------------------------
  ! Read soil moisture
  if ( nested/=4 ) then
    wb(1:ifull,1:ms) = 20.5
    if ( all(wetfrac_found(1:ms)) ) then
      call fillhist4('wetfrac',wb,sea_a,fill_sea)
      wb(1:ifull,1:ms) = wb(1:ifull,1:ms) + 20. ! flag for fraction of field capacity
    else
      do k = 1,ms
        if ( wetfrac_found(k) ) then
          write(vname,'("wetfrac",I1.1)') k
        else if ( wb_found(k) ) then
          write(vname,'("wb",I1.1)') k
        else if ( k<2 .and. wb_found(2) ) then
          vname = "wb2"
        else if ( k<2 ) then
          vname = "wfg"
        else if ( wb_found(6) ) then
          vname = "wb6"
        else
          vname = "wfb"
        end if
        if ( iotest ) then
          call histrd(iarchi,ier,vname,wb(:,k),ifull)
          if ( wetfrac_found(k) ) then
            wb(1:ifull,k) = wb(1:ifull,k) + 20. ! flag for fraction of field capacity
          end if
        else
          call histrd(iarchi,ier,vname,ucc,6*ik*ik)
          if ( wetfrac_found(k) ) then
            ucc(:) = ucc(:) + 20.   ! flag for fraction of field capacity
          end if
          call fill_cc1(ucc,sea_a,fill_sea)
          call doints1(ucc,wb(:,k))
        end if ! iotest
      end do
    end if
    !unpack field capacity into volumetric soil moisture
    if ( any(wb(1:ifull,1:ms)>10.) ) then
      wb(1:ifull,1:ms) = wb(1:ifull,1:ms) - 20.
      do k = 1,ms
        wb(1:ifull,k) = (1.-wb(1:ifull,k))*swilt(isoilm(1:ifull)) + wb(1:ifull,k)*sfc(isoilm(1:ifull))
        wb(1:ifull,k) = max( wb(1:ifull,k), 0.5*swilt(isoilm(1:ifull)) )
      end do
    end if
    call fillhist1('wetfac',wetfac,sea_a,fill_sea)
    where ( .not.land(1:ifull) )
      wetfac(:) = 1.
    end where
  end if ! nested/=4

  !------------------------------------------------------------------
  ! Read 10m wind speeds for ocean roughness length calculations
  if ( nested==0 ) then
    if ( u10_found ) then
      call gethist1('u10',u10)
    else
      u10 = sqrt(u(1:ifull,1)**2+v(1:ifull,1)**2)*log(10./0.001)/log(zmin/0.001)
    end if
  end if

  !------------------------------------------------------------------
  ! Average fields
  if ( nested==0 .and. lrestart ) then
    call gethist1('tscr_ave',tscr_ave)
    call gethist1('rhscr_ave',rhscr_ave)
    call gethist1('cbas_ave',cbas_ave)
    call gethist1('ctop_ave',ctop_ave)
    call gethist1('wb1_ave',wb_ave(:,1))
    call gethist1('wb2_ave',wb_ave(:,2))
    call gethist1('wb3_ave',wb_ave(:,3))
    call gethist1('wb4_ave',wb_ave(:,4))
    call gethist1('wb5_ave',wb_ave(:,5))
    call gethist1('wb6_ave',wb_ave(:,6))
  end if
  
  !------------------------------------------------------------------
  ! Read diagnostics and fluxes for zeroth time-step output
  if ( nested==0 .and. lrestart ) then
    call gethist1('tscrn',tscrn)
    call gethist1('qgscrn',qgscrn)
    call gethist1('eg',eg)
    call gethist1('fg',fg)
    call gethist1('taux',taux)
    call gethist1('tauy',tauy)
    call gethist1('ustar',ustar)
  end if
  
  !------------------------------------------------------------------
  ! Read boundary layer height for TKE-eps mixing and aerosols
  if ( nested==0 ) then
    call gethist1('pblh',pblh)
    pblh(:) = max(pblh(:), 1.)
  end if
  
  !------------------------------------------------------------------
  ! Read aerosol optical depth
  if ( abs(iaero)>=2 .and. nrad==5 ) then
    if ( nested==0 ) then
      call gethist1('sdust_vis',opticaldepth(:,1,1))
      call gethist1('ldust_vis',opticaldepth(:,2,1))
      call gethist1('so4_vis',opticaldepth(:,3,1))
      call gethist1('bc_vis',opticaldepth(:,5,1))
      call gethist1('oc_vis',opticaldepth(:,6,1))
      call gethist1('ssalt_vis',opticaldepth(:,7,1))
      call gethist1('od550aer',opticaldepth(:,4,1))
    end if    
  end if    
  
  !------------------------------------------------------------------
  ! Read carbon cycle data

  !if ( (nsib==6.or.nsib==7) .and. nhstest>=0 ) then
  !  if ( ccycle/=0 .and. carbon_found ) then
  !    call fillhist4('cplant',cplant,sea_a,fill_sea)
  !    call fillhist4('nplant',niplant,sea_a,fill_sea)
  !    call fillhist4('pplant',pplant,sea_a,fill_sea)
  !    call fillhist4('clitter',clitter,sea_a,fill_sea)
  !    call fillhist4('nlitter',nilitter,sea_a,fill_sea)
  !    call fillhist4('plitter',plitter,sea_a,fill_sea)
  !    call fillhist4('csoil',csoil,sea_a,fill_sea)
  !    call fillhist4('nsoil',nisoil,sea_a,fill_sea)
  !    call fillhist4('psoil',psoil,sea_a,fill_sea)
  !  end if
  !end if

  if ( (nsib==6.or.nsib==7) .and. nhstest>=0 ) then
    if ( ccycle/=0 .and. carbon2_found ) then
      if ( myid==0 ) then
        write(6,*) "-> Read CASA carbon tiles for each PFT"
      end if  
      ! allocated here to be used in cable_ccam2.f90 where it is deallocated
      allocate( carb_plant(ifull,mplant,3,11) )
      allocate( carb_litter(ifull,mlitter,3,11) )
      allocate( carb_soil(ifull,msoil,3,11) )
      carb_plant(:,:,:,:) = 0.
      carb_litter(:,:,:,:) = 0.
      carb_soil(:,:,:,:) = 0.
      do n = 1,10
        write(vname,'("p",I2.2,"_cplant")') n
        call fillhist4(vname,carb_plant(:,:,1,n),sea_a,fill_sea)
        write(vname,'("p",I2.2,"_nplant")') n
        call fillhist4(vname,carb_plant(:,:,2,n),sea_a,fill_sea)
        write(vname,'("p",I2.2,"_pplant")') n
        call fillhist4(vname,carb_plant(:,:,3,n),sea_a,fill_sea)
        write(vname,'("p",I2.2,"_clitter")') n
        call fillhist4(vname,carb_litter(:,:,1,n),sea_a,fill_sea)
        write(vname,'("p",I2.2,"_nlitter")') n
        call fillhist4(vname,carb_litter(:,:,2,n),sea_a,fill_sea)
        write(vname,'("p",I2.2,"_plitter")') n
        call fillhist4(vname,carb_litter(:,:,3,n),sea_a,fill_sea)
        write(vname,'("p",I2.2,"_csoil")') n
        call fillhist4(vname,carb_soil(:,:,1,n),sea_a,fill_sea)
        write(vname,'("p",I2.2,"_nsoil")') n
        call fillhist4(vname,carb_soil(:,:,2,n),sea_a,fill_sea)
        write(vname,'("p",I2.2,"_psoil")') n
        call fillhist4(vname,carb_soil(:,:,3,n),sea_a,fill_sea)
      end do
      n = 14 ! use index 11 to store pft 14
      write(vname,'("p",I2.2,"_cplant")') n
      call fillhist4(vname,carb_plant(:,:,1,11),sea_a,fill_sea)
      write(vname,'("p",I2.2,"_nplant")') n
      call fillhist4(vname,carb_plant(:,:,2,11),sea_a,fill_sea)
      write(vname,'("p",I2.2,"_pplant")') n
      call fillhist4(vname,carb_plant(:,:,3,1),sea_a,fill_sea)
      write(vname,'("p",I2.2,"_clitter")') n
      call fillhist4(vname,carb_litter(:,:,1,11),sea_a,fill_sea)
      write(vname,'("p",I2.2,"_nlitter")') n
      call fillhist4(vname,carb_litter(:,:,2,11),sea_a,fill_sea)
      write(vname,'("p",I2.2,"_plitter")') n
      call fillhist4(vname,carb_litter(:,:,3,11),sea_a,fill_sea)
      write(vname,'("p",I2.2,"_csoil")') n
      call fillhist4(vname,carb_soil(:,:,1,11),sea_a,fill_sea)
      write(vname,'("p",I2.2,"_nsoil")') n
      call fillhist4(vname,carb_soil(:,:,2,11),sea_a,fill_sea)
      write(vname,'("p",I2.2,"_psoil")') n
      call fillhist4(vname,carb_soil(:,:,3,11),sea_a,fill_sea)
    end if
  end if
  
  !------------------------------------------------------------------
  ! Read urban data
  if ( nurban/=0 .and. nhstest>=0 ) then
    if ( urban1_found ) then
      ! restart  
      do ifrac = 1,nfrac
        write(vname,'("t",I1.1,"_rooftgg")') ifrac   
        call fillhist4u(vname,"rooftemp",ifrac,nourban_a,fill_nourban,0.,1)
        write(vname,'("t",I1.1,"_waletgg")') ifrac
        call fillhist4u(vname,"walletemp",ifrac,nourban_a,fill_nourban,0.,1)
        write(vname,'("t",I1.1,"_walwtgg")') ifrac
        call fillhist4u(vname,"wallwtemp",ifrac,nourban_a,fill_nourban,0.,1)
        write(vname,'("t",I1.1,"_roadtgg")') ifrac
        call fillhist4u(vname,"roadtemp",ifrac,nourban_a,fill_nourban,0.,1)
        write(vname,'("t",I1.1,"_slabtgg")') ifrac
        call fillhist4u(vname,"slabtemp",ifrac,nourban_a,fill_nourban,0.,1)
        write(vname,'("t",I1.1,"_intmtgg")') ifrac
        call fillhist4u(vname,"intmtemp",ifrac,nourban_a,fill_nourban,0.,1)
        write(vname,'("t",I1.1,"_roomtgg1")') ifrac
        call fillhist1u(vname,"roomtemp",ifrac,nourban_a,fill_nourban,0.,1)
        write(vname,'("t",I1.1,"_urbnsmc")') ifrac
        call fillhist1u(vname,"canyonsoilmoisture",ifrac,nourban_a,fill_nourban,0.,0)
        write(vname,'("t",I1.1,"_urbnsmr")') ifrac
        call fillhist1u(vname,"roofsoilmoisture",ifrac,nourban_a,fill_nourban,0.,0)
        write(vname,'("t",I1.1,"_roofwtr")') ifrac
        call fillhist1u(vname,"roofsurfacewater",ifrac,nourban_a,fill_nourban,0.,0)
        write(vname,'("t",I1.1,"_roadwtr")') ifrac
        call fillhist1u(vname,"roadsurfacewater",ifrac,nourban_a,fill_nourban,0.,0)
        write(vname,'("t",I1.1,"_urbwtrc")') ifrac
        call fillhist1u(vname,"canyonleafwater",ifrac,nourban_a,fill_nourban,0.,0)
        write(vname,'("t",I1.1,"_urbwtrr")') ifrac
        call fillhist1u(vname,"roofleafwater",ifrac,nourban_a,fill_nourban,0.,0)
        write(vname,'("t",I1.1,"_roofsnd")') ifrac
        call fillhist1u(vname,"roofsnowdepth",ifrac,nourban_a,fill_nourban,0.,0)
        write(vname,'("t",I1.1,"_roadsnd")') ifrac
        call fillhist1u(vname,"roadsnowdepth",ifrac,nourban_a,fill_nourban,0.,0)
        write(vname,'("t",I1.1,"_roofden")') ifrac
        call fillhist1u(vname,"roofsnowdensity",ifrac,nourban_a,fill_nourban,100.,2)
        write(vname,'("t",I1.1,"_roadden")') ifrac
        call fillhist1u(vname,"roadsnowdensity",ifrac,nourban_a,fill_nourban,100.,2)
        write(vname,'("t",I1.1,"_roofsna")') ifrac
        call fillhist1u(vname,"roofsnowalbedo",ifrac,nourban_a,fill_nourban,0.85,2)
        write(vname,'("t",I1.1,"_roadsna")') ifrac
        call fillhist1u(vname,"roadsnowalbedo",ifrac,nourban_a,fill_nourban,0.85,2)
      end do
    else if ( urban2_found ) then
      ! nested with urban data  
      do ifrac = 1,nfrac  
        call fillhist4u("rooftgg","rooftemp",ifrac,nourban_a,fill_nourban,0.,1)
        call fillhist4u("waletgg","walletemp",ifrac,nourban_a,fill_nourban,0.,1)
        call fillhist4u("walwtgg","wallwtemp",ifrac,nourban_a,fill_nourban,0.,1)
        call fillhist4u("roadtgg","roadtemp",ifrac,nourban_a,fill_nourban,0.,1)
        call fillhist4u("slabtgg","slabtemp",ifrac,nourban_a,fill_nourban,0.,1)
        call fillhist4u("intmtgg","intmtemp",ifrac,nourban_a,fill_nourban,0.,1)
      end do  
    else
      ! nested without urban data
      if ( myid==0 .and. nmaxpr==1 ) then
        write(6,*) "Use tsu for urban data"  
      end if
      call gethist1("tsu",dum6)
      dum6 = abs(dum6)
      dum6 = min( max( dum6, 170. ), 380. )
      where ( dum6>150. )  
        dum6 = dum6 - urbtemp
      end where
      dum6r8 = real(dum6,8)
      do ifrac = 1,nfrac
        do k = 1,5
          write(vname,'("rooftemp",I1.1)') k
          call uclem_loadd(dum6r8,vname,ifrac,0)
          write(vname,'("walletemp",I1.1)') k
          call uclem_loadd(dum6r8,vname,ifrac,0)
          write(vname,'("wallwtemp",I1.1)') k
          call uclem_loadd(dum6r8,vname,ifrac,0)
          write(vname,'("roadtemp",I1.1)') k
          call uclem_loadd(dum6r8,vname,ifrac,0)
          write(vname,'("slabtemp",I1.1)') k
          call uclem_loadd(dum6r8,vname,ifrac,0)
          write(vname,'("intmtemp",I1.1)') k
          call uclem_loadd(dum6r8,vname,ifrac,0)
        end do  
      end do        
    end if    
  end if
  
  ! special nllp data
  if ( nested==0 .and. nextout>=4 .and. nllp==3 ) then
    if ( nllp_found ) then  
      call gethist4a('del_lat',tr(:,:,ngas+1),7)  
      call gethist4a('del_lon',tr(:,:,ngas+2),7)  
      call gethist4a('del_p',tr(:,:,ngas+3),7)
    end if
  end if

  ! k-eps data
  if ( nested==0 .and. abs(nmlo)>=1 .and. abs(nmlo)<=9 ) then
    mlodwn(:,:,5:6) = 0.
    if ( mlo2_found ) then
      if ( oclosure==1 ) then
        call fillhist4o('tkeo',mlodwn(:,:,5),land_3d,fill_floor,ocndwn(:,1))  
        call fillhist4o('epso',mlodwn(:,:,6),land_3d,fill_floor,ocndwn(:,1))
        call fillhist4o('uo_ema',oo,land_3d,fill_floor,ocndwn(:,1))
        do k = 1,wlev
          call mloimport('u_ema',oo(:,k),k,0)
        end do  
        call fillhist4o('vo_ema',oo,land_3d,fill_floor,ocndwn(:,1))
        do k = 1,wlev
          call mloimport('v_ema',oo(:,k),k,0)
        end do  
        call fillhist4o('wo_ema',oo,land_3d,fill_floor,ocndwn(:,1))
        do k = 1,wlev
          call mloimport('w_ema',oo(:,k),k,0)
        end do         
        call fillhist4o('temp_ema',oo,land_3d,fill_floor,ocndwn(:,1))
        do k = 1,wlev
          call mloimport('temp_ema',oo(:,k),k,0)
        end do  
        call fillhist4o('sal_ema',oo,land_3d,fill_floor,ocndwn(:,1))
        do k = 1,wlev
          call mloimport('sal_ema',oo(:,k),k,0)
        end do  
      end if  
    end if  
  end if
  
  ! -----------------------------------------------------------------
  ! Read cloud fields
  if ( nested==0 ) then
    call gethist4a('qfg',qfg,5)               ! CLOUD FROZEN WATER
    qfg(1:ifull,1:kl) = max( qfg(1:ifull,1:kl), 0. )
    call gethist4a('qlg',qlg,5)               ! CLOUD LIQUID WATER
    qlg(1:ifull,1:kl) = max( qlg(1:ifull,1:kl), 0. )
    if ( ncloud>=2 ) then
      call gethist4a('qrg',qrg,5)             ! RAIN
      qrg(1:ifull,1:kl) = max( qrg(1:ifull,1:kl), 0. )
    end if
    if ( ncloud>=3 ) then
      call gethist4a('qsng',qsng,5)           ! SNOW
      qsng(1:ifull,1:kl) = max( qsng(1:ifull,1:kl), 0. )
      call gethist4a('qgrg',qgrg,5)           ! GRAUPEL
      qgrg(1:ifull,1:kl) = max( qgrg(1:ifull,1:kl), 0. )
    end if
    call gethist4a('cfrac',cfrac,5)           ! CLOUD FRACTION
    cfrac(1:ifull,1:kl) = max( cfrac(1:ifull,1:kl), 0. )
    cfrac(1:ifull,1:kl) = min( cfrac(1:ifull,1:kl), 1. )
    call gethist4a('stratcf',stratcloud,5)
    stratcloud(1:ifull,1:kl) = max( stratcloud(1:ifull,1:kl), 0. )
    stratcloud(1:ifull,1:kl) = min( stratcloud(1:ifull,1:kl), 1. )
    if ( ncloud>=2 ) then
      call gethist4a('rfrac',rfrac,5)         ! RAIN FRACTION
      rfrac(1:ifull,1:kl) = max( rfrac(1:ifull,1:kl), 0. )
      rfrac(1:ifull,1:kl) = min( rfrac(1:ifull,1:kl), 1. )
    end if
    if ( ncloud>=3 ) then
      call gethist4a('sfrac',sfrac,5)         ! SNOW FRACTION
      sfrac(1:ifull,1:kl) = max( sfrac(1:ifull,1:kl), 0. )
      sfrac(1:ifull,1:kl) = min( sfrac(1:ifull,1:kl), 1. )
      call gethist4a('gfrac',gfrac,5)         ! GRAUPEL FRACTION
      gfrac(1:ifull,1:kl) = max( gfrac(1:ifull,1:kl), 0. )
      gfrac(1:ifull,1:kl) = min( gfrac(1:ifull,1:kl), 1. )
    end if
    call gethist4a('strat_rt',rad_tend,5)    ! STRAT RAD TENDENCY
    call gethist4a('strat_tt',trb_tend,5)    ! STRAT TURB TENDENCY
    call gethist4a('strat_tq',trb_qend,5)    ! STRAT TURB TENDENCY
    if ( ncloud>=100 .and. ncloud<=120 ) then
      call gethist4a('ni',ni,5)               ! Ice number concentration
      call gethist4a('nr',nr,5)               ! Rain number concentration
      call gethist4a('ns',ni,5)               ! Snow number concentration
      ni(1:ifull,1:kl) = max( ni(1:ifull,1:kl), 0. )
      nr(1:ifull,1:kl) = max( nr(1:ifull,1:kl), 0. )
      ns(1:ifull,1:kl) = max( ns(1:ifull,1:kl), 0. )
    end if    
  end if   ! (nested==0)

  !------------------------------------------------------------------
  ! TKE-eps data
  if ( nested==0 .and. (nvmix==6.or.nvmix==9) ) then
    call gethist4a('tke',tke,5)
    if ( all(tke(1:ifull,:)<1.e-20) ) tke(1:ifull,:)=1.5E-4
    call gethist4a('eps',eps,5)
    if  (all(eps(1:ifull,:)<1.e-20) ) eps(1:ifull,:)=1.E-7
    call gethist4a('u_ema',u_ema,5)
    call gethist4a('v_ema',v_ema,5)
    call gethist4a('w_ema',w_ema,5)
    call gethist4a('thetal_ema',thetal_ema,5)
    call gethist4a('qv_ema',qv_ema,5)
    call gethist4a('ql_ema',ql_ema,5)
    call gethist4a('qf_ema',qf_ema,5)
    call gethist4a('cf_ema',cf_ema,5)    
    call gethist4a('tke_ema',tke_ema,5)
  end if
  if ( nested==0 ) then
    call gethist4a('rkm',rkmsave,5)  
    call gethist4a('rkh',rkhsave,5)
  end if

  !------------------------------------------------------------------
  ! Tracer data
  if ( nested==0 .and. ngas>0 ) then
    if ( lrestart_tracer ) then  
      do igas = 1,ngas
        write(trnum,'(i4.4)') igas
        call gethist4a('tr'//trnum,tr(:,:,igas),7)
      end do
    end if  
  end if
  
  ! -----------------------------------------------------------------
  ! Restart fields
  if ( nested==0 ) then
    ! ATMOSPHERE DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    dpsldt(:,:) = -999.
    phi_nh(:,:) = 0.
    sdot(:,:) = -999.
    pslx(:,:) = -999.
    savu(:,:) = -999.
    savv(:,:) = -999.
    savu1(:,:) = -999.
    savv1(:,:) = -999.
    savu2(:,:) = -999.
    savv2(:,:) = -999.
    if ( lrestart ) then
      call gethist4('dpsldt',dpsldt)
      call gethist4('zgnhs',phi_nh)
      sdot(:,1) = 0.
      call gethist4('sdot',sdot(:,2:kk+1))
      call gethist4('pslx',pslx)
      call gethist4('savu',savu)
      call gethist4('savv',savv)
      call gethist4('savu1',savu1)
      call gethist4('savv1',savv1)
      call gethist4('savu2',savu2)
      call gethist4('savv2',savv2)
    end if

    ! OCEAN DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if ( abs(nmlo)>=3 .and. abs(nmlo)<=9 ) then
      oldu1(:,:) = 0.
      oldu2(:,:) = 0.
      oldv1(:,:) = 0.
      oldv2(:,:) = 0.
      ipice(:) = 0.
      if ( lrestart ) then
        call gethist4('old1_uo',oldu1)
        call gethist4('old1_vo',oldv1)
        call gethist4('old2_uo',oldu2)
        call gethist4('old2_vo',oldv2)            
        call gethist1('ipice',ipice)
      end if
    end if
    if ( nmlo/=0 .and. abs(nmlo)<=9 ) then
      ocndwn(:,3:6) = 0.
      if ( lrestart ) then
        call gethist1('old1_uotop',ocndwn(:,3))
        call gethist1('old1_votop',ocndwn(:,4))
        call gethist1('old1_uobot',ocndwn(:,5))
        call gethist1('old1_vobot',ocndwn(:,6))
      end if  
    end if    
       
  end if ! (nested==0)

  ! -----------------------------------------------------------------
  ! soil ice and snow data
  if ( nested/=4 ) then
    call gethist4('wbice',wbice) ! SOIL ICE
    call gethist4('tggsn',tggsn)
    if ( all(tggsn(1:ifull,1:3)<1.e-20) ) tggsn(1:ifull,1:3) = 270.
    call gethist4('smass',smass)
    call gethist4('ssdn',ssdn)
    do k = 1,3
      if ( all(ssdn(1:ifull,k)<1.e-20) ) then
        where ( snowd(1:ifull)>100. )
          ssdn(1:ifull,k)=240.
        elsewhere
          ssdn(1:ifull,k)=140.
        end where
      end if
    end do
    ssdnn(1:ifull) = ssdn(1:ifull,1)
    call gethist1('snage',snage)
    call gethist1('sflag',dum6)
    isflag(1:ifull) = nint(dum6(1:ifull))
  end if ! nested/=4  

  ! -----------------------------------------------------------------
  ! Misc fields
  ! sgsave is needed for convection
  if ( nested==0 ) then
    call gethist1('sgsave',sgsave)
    if ( lrestart_radiation ) then
      call gethist1('rgsave',rgsave)
      call gethist1('rtu',rt)
      call gethist1('rgdn',rgdn)
      call gethist1('rgn',rgn)
      call gethist1('rgc',rgclr)
      call gethist1('sint_amp',sint_amp)
      call gethist1('sout_amp',sout_amp)
      call gethist1('soutclr_amp',soutclr_amp)
      call gethist1('sgdn_amp',sgdn_amp)
      call gethist1('sgdndir_amp',sgdndir_amp)
      call gethist1('sgn_amp',sgn_amp)
      call gethist1('sgclr_amp',sgclr_amp)
      call gethist1('sgdclr_amp',sgdclr_amp)
      call gethist1('dni_amp',dni_amp)
      call gethist1('fbeamvis',fbeamvis)
      call gethist1('fbeamnir',fbeamnir)
      call gethist1('swrsave',swrsave)
      call gethist1('cloudlo',cloudlo)
      call gethist1('cloudmi',cloudmi)
      call gethist1('cloudhi',cloudhi)
      call gethist4('sw_tend_amp',sw_tend_amp)
      call gethist4('lw_tend',lw_tend)
    end if
    if ( lrestart_convection ) then
      call gethist4('dmsedt_rad',dmsedt_rad)
      call gethist4('dmsedt_pbl',dmsedt_pbl)
    end if
  end if
        
endif    ! (nested/=1.and.nested/=3)

!**************************************************************
! This is the end of reading the initial arrays
!**************************************************************  

deallocate( ucc, tss_a, ucc6 )

! -------------------------------------------------------------------
! tgg holds file surface temperature when there is no MLO
if ( nmlo==0 .or. abs(nmlo)>9 ) then
  where ( .not.land(1:ifull) )
    tgg(1:ifull,1) = tss(1:ifull)
  end where
end if

! -------------------------------------------------------------------
! set-up for next read of file
iarchi  = iarchi + 1
kdate_s = kdate_r
ktime_s = ktime_r + 1

if ( myid==0 .and. nested==0 ) then
  write(6,*) "Final lrestart, lrestart_radiation, lrestart_tracer = ",lrestart,lrestart_radiation,lrestart_tracer
end if

return
end subroutine onthefly_work


! *****************************************************************************
! INTERPOLATION ROUTINES                         

! Main interface

subroutine doints1(s,sout)
      
use cc_mpi                 ! CC MPI routines
use infile                 ! Input file routines
use newmpar_m              ! Grid parameters
use parm_m                 ! Model configuration

integer n, iq, mm, idel, jdel
integer iin, jjn, idel_l, jdel_l, no, w, i, j, nn
real, dimension(fwsize), intent(in) :: s
real, dimension(ifull), intent(inout) :: sout
real, dimension(-1:pipan+2,-1:pjpan+2,pnpan,size(filemap_req)) :: abuf
real xxg, yyg, cmin, cmax
real dmul_2, dmul_3, cmul_1, cmul_2, cmul_3, cmul_4
real emul_1, emul_2, emul_3, emul_4, rmul_1, rmul_2, rmul_3, rmul_4
real sx_0m,sx_1m,sx_m0,sx_00,sx_10,sx_20,sx_m1,sx_01,sx_11,sx_21,sx_02,sx_12

call START_LOG(otf_ints_begin)

! This version distributes mutli-file data
call ccmpi_filewinget(abuf,s)

do iq = 1,ifull
  sout(iq) = 0.
end do
do mm = 1,m_fly       !  was 4, now may be 1
  do iq = 1,ifull
    ! target point
    n = nface4(iq,mm)
    idel = int(xg4(iq,mm))
    xxg = xg4(iq,mm) - real(idel)
    jdel = int(yg4(iq,mm))
    yyg = yg4(iq,mm) - real(jdel)
    ! grid index conversion
    iin = (idel-1)/pipan
    jjn = (jdel-1)/pjpan
    idel_l = idel - iin*pipan
    jdel_l = jdel - jjn*pjpan
    w = pijn2pw(iin,jjn,n)
    no = pijn2pn(iin,jjn,n)
    ! bi-cubic
    cmul_1 = (1.-xxg)*(2.-xxg)*(-xxg)/6.
    cmul_2 = (1.-xxg)*(2.-xxg)*(1.+xxg)/2.
    cmul_3 = xxg*(1.+xxg)*(2.-xxg)/2.
    cmul_4 = (1.-xxg)*(-xxg)*(1.+xxg)/6.
    dmul_2 = (1.-xxg)
    dmul_3 = xxg
    emul_1 = (1.-yyg)*(2.-yyg)*(-yyg)/6.
    emul_2 = (1.-yyg)*(2.-yyg)*(1.+yyg)/2.
    emul_3 = yyg*(1.+yyg)*(2.-yyg)/2.
    emul_4 = (1.-yyg)*(-yyg)*(1.+yyg)/6.
    sx_0m = abuf(idel_l,  jdel_l-1,no,w)
    sx_1m = abuf(idel_l+1,jdel_l-1,no,w)
    sx_m0 = abuf(idel_l-1,jdel_l,  no,w)
    sx_00 = abuf(idel_l,  jdel_l,  no,w)
    sx_10 = abuf(idel_l+1,jdel_l,  no,w)
    sx_20 = abuf(idel_l+2,jdel_l,  no,w)
    sx_m1 = abuf(idel_l-1,jdel_l+1,no,w)
    sx_01 = abuf(idel_l,  jdel_l+1,no,w)
    sx_11 = abuf(idel_l+1,jdel_l+1,no,w)
    sx_21 = abuf(idel_l+2,jdel_l+1,no,w)
    sx_02 = abuf(idel_l,  jdel_l+2,no,w)
    sx_12 = abuf(idel_l+1,jdel_l+2,no,w)
    cmin = min(sx_00,sx_01,sx_10,sx_11)
    cmax = max(sx_00,sx_01,sx_10,sx_11)
    rmul_1 = sx_0m*dmul_2 + sx_1m*dmul_3
    rmul_2 = sx_m0*cmul_1 + sx_00*cmul_2 + &
             sx_10*cmul_3 + sx_20*cmul_4
    rmul_3 = sx_m1*cmul_1 + sx_01*cmul_2 + &
             sx_11*cmul_3 + sx_21*cmul_4
    rmul_4 = sx_02*dmul_2 + sx_12*dmul_3
    sout(iq) = sout(iq) + min( max( cmin, rmul_1*emul_1 + rmul_2*emul_2    &
                                        + rmul_3*emul_3 + rmul_4*emul_4 ), &
                                    cmax )/real(m_fly) ! Bermejo & Staniforth
  end do    ! iq loop
end do      ! mm loop

call END_LOG(otf_ints_end)

return
end subroutine doints1

subroutine doints4(s,sout)
      
use cc_mpi                 ! CC MPI routines
use infile                 ! Input file routines
use newmpar_m              ! Grid parameters
use parm_m                 ! Model configuration

integer k, kx, kb, ke, kn, n, iq, mm, idel, jdel
integer iin, jjn, idel_l, jdel_l, no, w, i, j, nn
real, dimension(:,:), intent(in) :: s
real, dimension(:,:), intent(inout) :: sout
real, dimension(-1:pipan+2,-1:pjpan+2,pnpan,size(filemap_req),kblock) :: abuf
real xxg, yyg, cmin, cmax
real dmul_2, dmul_3, cmul_1, cmul_2, cmul_3, cmul_4
real emul_1, emul_2, emul_3, emul_4, rmul_1, rmul_2, rmul_3, rmul_4
real sx_0m,sx_1m,sx_m0,sx_00,sx_10,sx_20,sx_m1,sx_01,sx_11,sx_21,sx_02,sx_12

call START_LOG(otf_ints_begin)

kx = size(sout,2)
    
do kb = 1,kx,kblock
  ke = min(kb+kblock-1, kx)
  kn = ke - kb + 1
  ! This version distributes multi-file data
  call ccmpi_filewinget(abuf(:,:,:,:,1:kn),s(:,kb:ke))
  do k = 1,kn
    do iq = 1,ifull
      sout(iq,k+kb-1) = 0.
    end do
  end do
  do mm = 1,m_fly       !  was 4, now may be 1
    do k = 1,kn
      do iq = 1,ifull
        ! target point
        n = nface4(iq,mm)
        idel = int(xg4(iq,mm))
        xxg = xg4(iq,mm) - real(idel)
        jdel = int(yg4(iq,mm))
        yyg = yg4(iq,mm) - real(jdel)
        ! grid index conversion
        iin = (idel-1)/pipan
        jjn = (jdel-1)/pjpan
        idel_l = idel - iin*pipan
        jdel_l = jdel - jjn*pjpan
        w = pijn2pw(iin,jjn,n)
        no = pijn2pn(iin,jjn,n)
        ! bi-cubic
        cmul_1 = (1.-xxg)*(2.-xxg)*(-xxg)/6.
        cmul_2 = (1.-xxg)*(2.-xxg)*(1.+xxg)/2.
        cmul_3 = xxg*(1.+xxg)*(2.-xxg)/2.
        cmul_4 = (1.-xxg)*(-xxg)*(1.+xxg)/6.
        dmul_2 = (1.-xxg)
        dmul_3 = xxg
        emul_1 = (1.-yyg)*(2.-yyg)*(-yyg)/6.
        emul_2 = (1.-yyg)*(2.-yyg)*(1.+yyg)/2.
        emul_3 = yyg*(1.+yyg)*(2.-yyg)/2.
        emul_4 = (1.-yyg)*(-yyg)*(1.+yyg)/6.
        sx_0m = abuf(idel_l,  jdel_l-1,no,w,k)
        sx_1m = abuf(idel_l+1,jdel_l-1,no,w,k)
        sx_m0 = abuf(idel_l-1,jdel_l,  no,w,k)
        sx_00 = abuf(idel_l,  jdel_l,  no,w,k)
        sx_10 = abuf(idel_l+1,jdel_l,  no,w,k)
        sx_20 = abuf(idel_l+2,jdel_l,  no,w,k)
        sx_m1 = abuf(idel_l-1,jdel_l+1,no,w,k)
        sx_01 = abuf(idel_l,  jdel_l+1,no,w,k)
        sx_11 = abuf(idel_l+1,jdel_l+1,no,w,k)
        sx_21 = abuf(idel_l+2,jdel_l+1,no,w,k)
        sx_02 = abuf(idel_l,  jdel_l+2,no,w,k)
        sx_12 = abuf(idel_l+1,jdel_l+2,no,w,k)
        cmin = min(sx_00,sx_01,sx_10,sx_11)
        cmax = max(sx_00,sx_01,sx_10,sx_11)
        rmul_1 = sx_0m*dmul_2 + sx_1m*dmul_3
        rmul_2 = sx_m0*cmul_1 + sx_00*cmul_2 + sx_10*cmul_3 + sx_20*cmul_4
        rmul_3 = sx_m1*cmul_1 + sx_01*cmul_2 + sx_11*cmul_3 + sx_21*cmul_4
        rmul_4 = sx_02*dmul_2 + sx_12*dmul_3
        sout(iq,k+kb-1) = sout(iq,k+kb-1) + min( max( cmin, rmul_1*emul_1 + rmul_2*emul_2    &
                                                          + rmul_3*emul_3 + rmul_4*emul_4 ), &
                                                      cmax )/real(m_fly) ! Bermejo & Staniforth
      end do    ! iq loop
    end do      ! k loop
  end do        ! mm loop
end do ! kb

call END_LOG(otf_ints_end)

return
end subroutine doints4

! *****************************************************************************
! FILL ROUTINES

subroutine fill_cc1(a_io,land_a,fill_count)
      
! routine fills in interior of an array which has undefined points
! this version is for multiple input files

use cc_mpi          ! CC MPI routines
use infile          ! Input file routines

integer, intent(inout) :: fill_count
integer nrem, j, n
integer ncount, cc, ipf, local_count
integer i
real, parameter :: value=999.       ! missing value flag
real, dimension(fwsize), intent(inout) :: a_io
real, dimension(0:pipan+1,0:pjpan+1,pnpan,mynproc) :: c_io
real csum, ccount
logical, dimension(fwsize), intent(in) :: land_a

! only perform fill on processors reading input files
if ( fwsize==0 ) return

! ignore fill if land_a is trivial
if ( allowtrivialfill ) then
  ncount = count( .not.land_a(1:fwsize) )
  call ccmpi_allreduce(ncount,nrem,'sum',comm_ip)
  if ( nrem==0 .or. nrem==6*ik*ik ) return
end if

call START_LOG(otf_fill_begin)

where ( land_a(1:fwsize) )
  a_io(1:fwsize) = value
end where

nrem = 1
local_count = 0
c_io = value

do while ( nrem>0 )
  do ipf = 1,mynproc
    do n = 1,pnpan
      do j = 1,pjpan
        do i = 1,pipan
          cc = i + (j-1)*pipan + (n-1)*pipan*pjpan + (ipf-1)*pipan*pjpan*pnpan
          c_io(i,j,n,ipf) = a_io(cc)
        end do
      end do
    end do
  end do  
  ncount = count( abs(a_io(1:fwsize)-value)<1.E-20 )
  call ccmpi_filebounds(c_io,comm_ip,corner=.true.)
  ! update body
  if ( ncount>0 ) then
    do ipf = 1,mynproc
      do n = 1,pnpan
        do j = 1,pjpan
          do i = 1,pipan
            cc = (j-1)*pipan + (n-1)*pipan*pjpan + (ipf-1)*pipan*pjpan*pnpan
            if ( abs(a_io(cc+i)-value)<1.e-20 ) then
              csum = 0.
              ccount = 0.
              if ( abs(c_io(i-1,j-1,n,ipf)-value)>=1.e-20 ) then
                csum = csum + c_io(i-1,j-1,n,ipf)
                ccount = ccount + 1.
              end if
              if ( abs(c_io(i,j-1,n,ipf)-value)>=1.e-20 ) then
                csum = csum + c_io(i,j-1,n,ipf)
                ccount = ccount + 1.
              end if
              if ( abs(c_io(i+1,j-1,n,ipf)-value)>=1.e-20 ) then
                csum = csum + c_io(i+1,j-1,n,ipf)
                ccount = ccount + 1.
              end if
              if ( abs(c_io(i-1,j,n,ipf)-value)>=1.e-20 ) then
                csum = csum + c_io(i-1,j,n,ipf)
                ccount = ccount + 1.
              end if
              if ( abs(c_io(i+1,j,n,ipf)-value)>=1.e-20 ) then
                csum = csum + c_io(i+1,j,n,ipf)
                ccount = ccount + 1.
              end if
              if ( abs(c_io(i-1,j+1,n,ipf)-value)>=1.e-20 ) then
                csum = csum + c_io(i-1,j+1,n,ipf)
                ccount = ccount + 1.
              end if
              if ( abs(c_io(i,j+1,n,ipf)-value)>=1.e-20 ) then
                csum = csum + c_io(i,j+1,n,ipf)
                ccount = ccount + 1.
              end if
              if ( abs(c_io(i+1,j+1,n,ipf)-value)>=1.e-20 ) then
                csum = csum + c_io(i+1,j+1,n,ipf)
                ccount = ccount + 1.
              end if
              if ( ccount>0. ) then        
                a_io(cc+i) = csum/ccount
              end if
            end if
          end do
        end do
      end do
    end do
    ncount = count( abs(a_io(1:fwsize)-value)<1.E-6 )  
  end if  
  ! test for convergence
  local_count = local_count + 1
  if ( local_count==fill_count ) then
    nrem = 0
  else if ( local_count>fill_count ) then
    call ccmpi_allreduce(ncount,nrem,'sum',comm_ip)
    if ( nrem==6*ik*ik ) then
      ! Cannot perform fill as all points are trivial    
      nrem = 0
    end if
  end if  
end do
      
fill_count = local_count

call END_LOG(otf_fill_end)

return
end subroutine fill_cc1

subroutine fill_cc4_3d(a_io,land_3d,fill_count)
      
! routine fills in interior of an array which has undefined points
! this version is distributed over processes with input files

use cc_mpi          ! CC MPI routines
use infile          ! Input file routines

real, dimension(:,:), intent(inout) :: a_io
integer, intent(inout) :: fill_count
integer j, n, k, kx
integer cc, ipf, local_count
integer i
integer, dimension(size(a_io,2)) :: ncount, nrem
real, parameter :: value=999.       ! missing value flag
real, dimension(0:pipan+1,0:pjpan+1,pnpan,mynproc,size(a_io,2)) :: c_io
real csum, ccount
logical, dimension(:,:), intent(in) :: land_3d

kx = size(a_io,2)

! only perform fill on processors reading input files
if ( fwsize==0 ) return

! ignore fill if land_a is trivial
if ( allowtrivialfill ) then
  do k = 1,kx  
    ncount(k) = count( .not.land_3d(1:fwsize,k) )
  end do  
  call ccmpi_allreduce(ncount(1:kx),nrem(1:kx),'sum',comm_ip)
  if ( all(nrem(:)==0) .or. all(nrem(:)==6*ik*ik) ) return
end if

call START_LOG(otf_fill_begin)


do k = 1,kx
  where ( land_3d(1:fwsize,k) )
    a_io(1:fwsize,k) = value
  end where
end do

nrem(:) = 1
local_count = 0
c_io = value

do while ( any(nrem(:)>0) )
  do k = 1,kx
    do ipf = 1,mynproc
      do n = 1,pnpan
        do j = 1,pjpan
          do i = 1,pipan
            cc = i + (j-1)*pipan + (n-1)*pipan*pjpan + (ipf-1)*pipan*pjpan*pnpan
            c_io(i,j,n,ipf,k) = a_io(cc,k)
          end do  
        end do
      end do
    end do
  end do  
  call ccmpi_filebounds(c_io,comm_ip,corner=.true.)
  ncount(1:kx) = count( abs(a_io(1:fwsize,1:kx)-value)<1.E-6, dim=1 )
  ! update body
  do k = 1,kx
    if ( ncount(k)>0 ) then
      do ipf = 1,mynproc
        do n = 1,pnpan
          do j = 1,pjpan
            do i = 1,pipan
              cc = (j-1)*pipan + (n-1)*pipan*pjpan + (ipf-1)*pipan*pjpan*pnpan
              if ( abs(a_io(cc+i,k)-value)<1.e-20 ) then
                csum = 0.
                ccount = 0.
                if ( abs(c_io(i-1,j-1,n,ipf,k)-value)>=1.e-20 ) then
                  csum = csum + c_io(i-1,j-1,n,ipf,k)
                  ccount = ccount + 1.
                end if
                if ( abs(c_io(i,j-1,n,ipf,k)-value)>=1.e-20 ) then
                  csum = csum + c_io(i,j-1,n,ipf,k)
                  ccount = ccount + 1.
                end if
                if ( abs(c_io(i+1,j-1,n,ipf,k)-value)>=1.e-20 ) then
                  csum = csum + c_io(i+1,j-1,n,ipf,k)
                  ccount = ccount + 1.
                end if
                if ( abs(c_io(i-1,j,n,ipf,k)-value)>=1.e-20 ) then
                  csum = csum + c_io(i-1,j,n,ipf,k)
                  ccount = ccount + 1.
                end if
                if ( abs(c_io(i+1,j,n,ipf,k)-value)>=1.e-20 ) then
                  csum = csum + c_io(i+1,j,n,ipf,k)
                  ccount = ccount + 1.
                end if
                if ( abs(c_io(i-1,j+1,n,ipf,k)-value)>=1.e-20 ) then
                  csum = csum + c_io(i-1,j+1,n,ipf,k)
                  ccount = ccount + 1.
                end if
                if ( abs(c_io(i,j+1,n,ipf,k)-value)>=1.e-20 ) then
                  csum = csum + c_io(i,j+1,n,ipf,k)
                  ccount = ccount + 1.
                end if
                if ( abs(c_io(i+1,j+1,n,ipf,k)-value)>=1.e-20 ) then
                  csum = csum + c_io(i+1,j+1,n,ipf,k)
                  ccount = ccount + 1.
                end if
                if ( ccount>0. ) then        
                  a_io(cc+i,k) = csum/ccount
                end if
              end if
            end do
          end do
        end do
      end do
      ncount(k) = count( abs(a_io(1:fwsize,k)-value)<1.E-6 )
    end if
  end do
  ! test for convergence
  local_count = local_count + 1
  if ( local_count==fill_count ) then
    nrem(:) = 0  
  else if ( local_count>fill_count ) then
    call ccmpi_allreduce(ncount(1:kx),nrem(1:kx),'sum',comm_ip)
    do k = 1,kx
      if ( nrem(k)==6*ik*ik ) then
        ! Cannot perform fill as all points are trivial  
        nrem(k) = 0
      end if
    end do
  end if  
end do

fill_count = local_count

call END_LOG(otf_fill_end)

return
end subroutine fill_cc4_3d

subroutine fill_cc4_1(a_io,land_a,fill_count)

integer, intent(inout) :: fill_count
integer k
real, dimension(:,:), intent(inout) :: a_io
logical, dimension(:), intent(in) :: land_a
logical, dimension(size(a_io,1),size(a_io,2)) :: land_3d

do k = 1,size(a_io,2)
  land_3d(:,k) = land_a
end do
call fill_cc4_3d(a_io,land_3d,fill_count)

return
end subroutine fill_cc4_1

! *****************************************************************************
! OROGRAPHIC ADJUSTMENT ROUTINES

subroutine mslpx(pmsl,psl,zss,t,siglev)

use const_phys             ! Physical constants
use newmpar_m              ! Grid parameters
use parm_m                 ! Model configuration
use sigs_m                 ! Atmosphere sigma levels
      
!     this one will ignore negative zss (i.e. over the ocean)

real, intent(in) :: siglev
real, dimension(fwsize), intent(in) :: psl, zss, t
real, dimension(fwsize), intent(out) :: pmsl
real, dimension(fwsize) :: dlnps, phi1, tav, tsurf

phi1(1:fwsize)  = t(1:fwsize)*rdry*(1.-siglev)/siglev ! phi of sig(lev) above sfce
tsurf(1:fwsize) = t(1:fwsize)+phi1(1:fwsize)*stdlapse/grav
tav(1:fwsize)   = tsurf(1:fwsize)+max(0.,zss(1:fwsize))*.5*stdlapse/grav
dlnps(1:fwsize) = max(0.,zss(1:fwsize))/(rdry*tav(1:fwsize))
pmsl(1:fwsize)  = 1.e5*exp(psl(1:fwsize)+dlnps(1:fwsize))

return
end subroutine mslpx
      
subroutine to_pslx(pmsl,psl,zss,t,siglev)

use cc_mpi, only : mydiag  ! CC MPI routines
use const_phys             ! Physical constants
use newmpar_m              ! Grid parameters
use parm_m                 ! Model configuration
use sigs_m                 ! Atmosphere sigma levels
      
real, intent(in) :: siglev
real, dimension(ifull), intent(in) :: pmsl, zss, t
real, dimension(ifull), intent(out) :: psl
real, dimension(ifull) :: dlnps, phi1, tav, tsurf

phi1(1:ifull)  = t(1:ifull)*rdry*(1.-siglev)/siglev ! phi of sig(levk) above sfce
tsurf(1:ifull) = t(1:ifull) + phi1(1:ifull)*stdlapse/grav
tav(1:ifull)   = tsurf(1:ifull) + max( 0., zss(1:ifull) )*.5*stdlapse/grav
dlnps(1:ifull) = max( 0., zss(1:ifull))/(rdry*tav(1:ifull) )
psl(1:ifull)   = log(1.e-5*pmsl(1:ifull)) - dlnps(1:ifull)

#ifdef debug
if ( nmaxpr==1 .and. mydiag ) then
  write(6,*)'zs,t_lev,psl,pmsl ',zss(idjd),t(idjd),psl(idjd),pmsl(idjd)
end if
#endif

return
end subroutine to_pslx

subroutine retopo(psl,zsold,zss,t,qg)
!     in Jan 07, renamed recalc of ps from here, to reduce confusion     
!     this version (Aug 2003) allows -ve zsold (from spectral model),
!     but assumes new zss is positive for atmospheric purposes
!     this routine redefines psl, t to compensate for zsold going to zss
!     (but does not overwrite zss, ps themselves here)
!     called by indata and nestin for newtop>=1
!     nowadays just for ps and atmospheric fields Mon  08-23-1999
use cc_mpi, only : mydiag, ccmpi_abort
use const_phys
use diag_m
use newmpar_m
use parm_m
use sigs_m

real, dimension(ifull), intent(inout) :: psl
real, dimension(ifull), intent(in) :: zsold, zss
real, dimension(:,:), intent(inout) :: t, qg
real, dimension(ifull) :: psnew, psold, pslold
real, dimension(kl) :: told, qgold
real sig2
integer iq, k, kk, kold

pslold(1:ifull) = psl(1:ifull)
psold(1:ifull)  = 1.e5*exp(psl(1:ifull))
psl(1:ifull)    = psl(1:ifull) + (zsold(1:ifull)-zss(1:ifull))/(rdry*t(1:ifull,1))
psnew(1:ifull)  = 1.e5*exp(psl(1:ifull))

!     now alter temperatures to compensate for new topography
if ( ktau<100 .and. mydiag ) then
  write(6,*) 'retopo: zsold,zs,psold,psnew ',zsold(idjd),zss(idjd),psold(idjd),psnew(idjd)
  write(6,*) 'retopo: old t ',(t(idjd,k),k=1,kl)
  write(6,*) 'retopo: old qg ',(qg(idjd,k),k=1,kl)
end if  ! (ktau.lt.100)

do iq = 1,ifull
  qgold(1:kl) = qg(iq,1:kl)
  told(1:kl)  = t(iq,1:kl)
  kold = 2
  do k = 1,kl ! MJT suggestion
    sig2 = sig(k)*psnew(iq)/psold(iq)
    if ( sig2>=sig(1) ) then
      ! assume 6.5 deg/km, with dsig=.1 corresponding to 1 km
      t(iq,k) = told(1) + (sig2-sig(1))*6.5/0.1  
    else
      do kk = kold,kl-1
        if ( sig2>sig(kk) ) exit
      end do
      kold = kk
      t(iq,k)  = (told(kk)*(sig(kk-1)-sig2)+told(kk-1)*(sig2-sig(kk)))/(sig(kk-1)-sig(kk))
      qg(iq,k) = (qgold(kk)*(sig(kk-1)-sig2)+qgold(kk-1)*(sig2-sig(kk)))/(sig(kk-1)-sig(kk))
    end if
  end do  ! k loop
end do    ! iq loop

if ( ktau<100 .and. mydiag ) then
  write(6,*) 'retopo: new t ',(t(idjd,k),k=1,kl)
  write(6,*) 'retopo: new qg ',(qg(idjd,k),k=1,kl)
end if  ! (ktau.lt.100)

return
end subroutine retopo

! *****************************************************************************
! VECTOR INTERPOLATION ROUTINES

subroutine interpwind4(uct,vct,ucc,vcc)
      
use cc_mpi           ! CC MPI routines
use newmpar_m        ! Grid parameters
use vecsuv_m         ! Map to cartesian coordinates
      
integer k
real, dimension(fwsize,kk), intent(inout) :: ucc, vcc
real, dimension(fwsize,kk) :: wcc
real, dimension(fwsize) :: uc, vc, wc
real, dimension(ifull,kk), intent(out) :: uct, vct
real, dimension(ifull,kk) :: wct
real, dimension(ifull) :: newu, newv, neww

if ( fwsize>0 ) then
  do k = 1,kk
    ! first set up winds in Cartesian "source" coords            
    uc(1:fwsize) = axs_w(1:fwsize)*ucc(1:fwsize,k) + bxs_w(1:fwsize)*vcc(1:fwsize,k)
    vc(1:fwsize) = ays_w(1:fwsize)*ucc(1:fwsize,k) + bys_w(1:fwsize)*vcc(1:fwsize,k)
    wc(1:fwsize) = azs_w(1:fwsize)*ucc(1:fwsize,k) + bzs_w(1:fwsize)*vcc(1:fwsize,k)
    ! now convert to winds in "absolute" Cartesian components
    ucc(1:fwsize,k) = uc(1:fwsize)*rotpoles(1,1) + vc(1:fwsize)*rotpoles(1,2) + wc(1:fwsize)*rotpoles(1,3)
    vcc(1:fwsize,k) = uc(1:fwsize)*rotpoles(2,1) + vc(1:fwsize)*rotpoles(2,2) + wc(1:fwsize)*rotpoles(2,3)
    wcc(1:fwsize,k) = uc(1:fwsize)*rotpoles(3,1) + vc(1:fwsize)*rotpoles(3,2) + wc(1:fwsize)*rotpoles(3,3)
  end do    ! k loop
end if      ! fwsize>0
! interpolate all required arrays to new C-C positions
! do not need to do map factors and Coriolis on target grid
call doints4(ucc, uct)
call doints4(vcc, vct)
call doints4(wcc, wct)
  
do k = 1,kk
  ! now convert to "target" Cartesian components (transpose used)
  newu(1:ifull) = uct(1:ifull,k)*rotpole(1,1) + vct(1:ifull,k)*rotpole(2,1) + wct(1:ifull,k)*rotpole(3,1)
  newv(1:ifull) = uct(1:ifull,k)*rotpole(1,2) + vct(1:ifull,k)*rotpole(2,2) + wct(1:ifull,k)*rotpole(3,2)
  neww(1:ifull) = uct(1:ifull,k)*rotpole(1,3) + vct(1:ifull,k)*rotpole(2,3) + wct(1:ifull,k)*rotpole(3,3)
  ! then finally to "target" local x-y components
  uct(1:ifull,k) = ax(1:ifull)*newu(1:ifull) + ay(1:ifull)*newv(1:ifull) + az(1:ifull)*neww(1:ifull)
  vct(1:ifull,k) = bx(1:ifull)*newu(1:ifull) + by(1:ifull)*newv(1:ifull) + bz(1:ifull)*neww(1:ifull)
end do  ! k loop
  
return
end subroutine interpwind4

subroutine interpcurrent1(uct,vct,ucc,vcc,mask_a,fill_count)
      
use cc_mpi           ! CC MPI routines
use newmpar_m        ! Grid parameters
use vecsuv_m         ! Map to cartesian coordinates
      
integer, intent(inout) :: fill_count      
real, dimension(fwsize), intent(inout) :: ucc, vcc
real, dimension(fwsize) :: wcc
real, dimension(fwsize) :: uc, vc, wc
real, dimension(fwsize,3) :: uvwcc
real, dimension(ifull), intent(out) :: uct, vct
real, dimension(ifull) :: wct
real, dimension(ifull) :: newu, newv, neww
real, dimension(ifull,3) :: newuvw
logical, dimension(fwsize), intent(in) :: mask_a

if ( fwsize>0 ) then
  uc(1:fwsize) = axs_w(1:fwsize)*ucc(1:fwsize) + bxs_w(1:fwsize)*vcc(1:fwsize)
  vc(1:fwsize) = ays_w(1:fwsize)*ucc(1:fwsize) + bys_w(1:fwsize)*vcc(1:fwsize)
  wc(1:fwsize) = azs_w(1:fwsize)*ucc(1:fwsize) + bzs_w(1:fwsize)*vcc(1:fwsize)
  ! now convert to winds in "absolute" Cartesian components
  ucc(1:fwsize) = uc(1:fwsize)*rotpoles(1,1) + vc(1:fwsize)*rotpoles(1,2) + wc(1:fwsize)*rotpoles(1,3)
  vcc(1:fwsize) = uc(1:fwsize)*rotpoles(2,1) + vc(1:fwsize)*rotpoles(2,2) + wc(1:fwsize)*rotpoles(2,3)
  wcc(1:fwsize) = uc(1:fwsize)*rotpoles(3,1) + vc(1:fwsize)*rotpoles(3,2) + wc(1:fwsize)*rotpoles(3,3)
  ! interpolate all required arrays to new C-C positions
  ! do not need to do map factors and Coriolis on target grid
  uvwcc(:,1) = ucc
  uvwcc(:,2) = vcc
  uvwcc(:,3) = wcc
  call fill_cc4(uvwcc(:,1:3), mask_a, fill_count)
end if
call doints4(uvwcc(:,1:3), newuvw(:,1:3))
  
! now convert to "target" Cartesian components (transpose used)
uct(1:ifull) = newuvw(1:ifull,1)
vct(1:ifull) = newuvw(1:ifull,2)
wct(1:ifull) = newuvw(1:ifull,3)
newu(1:ifull) = uct(1:ifull)*rotpole(1,1) + vct(1:ifull)*rotpole(2,1) + wct(1:ifull)*rotpole(3,1)
newv(1:ifull) = uct(1:ifull)*rotpole(1,2) + vct(1:ifull)*rotpole(2,2) + wct(1:ifull)*rotpole(3,2)
neww(1:ifull) = uct(1:ifull)*rotpole(1,3) + vct(1:ifull)*rotpole(2,3) + wct(1:ifull)*rotpole(3,3)
! then finally to "target" local x-y components
uct(1:ifull) = ax(1:ifull)*newu(1:ifull) + ay(1:ifull)*newv(1:ifull) + az(1:ifull)*neww(1:ifull)
vct(1:ifull) = bx(1:ifull)*newu(1:ifull) + by(1:ifull)*newv(1:ifull) + bz(1:ifull)*neww(1:ifull)

return
end subroutine interpcurrent1

subroutine interpcurrent4(uct,vct,ucc,vcc,mask_3d,fill_count)
      
use cc_mpi           ! CC MPI routines
use newmpar_m        ! Grid parameters
use vecsuv_m         ! Map to cartesian coordinates
      
integer, intent(inout) :: fill_count
integer k
real, dimension(fwsize,ok), intent(inout) :: ucc, vcc
real, dimension(fwsize,ok) :: wcc
real, dimension(ifull,ok), intent(out) :: uct, vct
real, dimension(ifull,ok) :: wct
real, dimension(fwsize) :: uc, vc, wc
real, dimension(ifull) :: newu, newv, neww
logical, dimension(fwsize,ok), intent(in) :: mask_3d

if ( fwsize>0 ) then
  do k = 1,ok
    ! first set up currents in Cartesian "source" coords            
    uc(1:fwsize) = axs_w(1:fwsize)*ucc(1:fwsize,k) + bxs_w(1:fwsize)*vcc(1:fwsize,k)
    vc(1:fwsize) = ays_w(1:fwsize)*ucc(1:fwsize,k) + bys_w(1:fwsize)*vcc(1:fwsize,k)
    wc(1:fwsize) = azs_w(1:fwsize)*ucc(1:fwsize,k) + bzs_w(1:fwsize)*vcc(1:fwsize,k)
    ! now convert to winds in "absolute" Cartesian components
    ucc(1:fwsize,k) = uc(1:fwsize)*rotpoles(1,1) + vc(1:fwsize)*rotpoles(1,2) + wc(1:fwsize)*rotpoles(1,3)
    vcc(1:fwsize,k) = uc(1:fwsize)*rotpoles(2,1) + vc(1:fwsize)*rotpoles(2,2) + wc(1:fwsize)*rotpoles(2,3)
    wcc(1:fwsize,k) = uc(1:fwsize)*rotpoles(3,1) + vc(1:fwsize)*rotpoles(3,2) + wc(1:fwsize)*rotpoles(3,3)
  end do  ! k loop  
  ! interpolate all required arrays to new C-C positions
  ! do not need to do map factors and Coriolis on target grid
  call fill_cc4(ucc, mask_3d, fill_count)
  call fill_cc4(vcc, mask_3d, fill_count)
  call fill_cc4(wcc, mask_3d, fill_count)
end if
call doints4(ucc, uct)
call doints4(vcc, vct)
call doints4(wcc, wct)
  
do k = 1,ok
  ! now convert to "target" Cartesian components (transpose used)  
  newu(1:ifull) = uct(1:ifull,k)*rotpole(1,1) + vct(1:ifull,k)*rotpole(2,1) + wct(1:ifull,k)*rotpole(3,1)
  newv(1:ifull) = uct(1:ifull,k)*rotpole(1,2) + vct(1:ifull,k)*rotpole(2,2) + wct(1:ifull,k)*rotpole(3,2)
  neww(1:ifull) = uct(1:ifull,k)*rotpole(1,3) + vct(1:ifull,k)*rotpole(2,3) + wct(1:ifull,k)*rotpole(3,3)
  ! then finally to "target" local x-y components
  uct(1:ifull,k) = ax(1:ifull)*newu(1:ifull) + ay(1:ifull)*newv(1:ifull) + az(1:ifull)*neww(1:ifull)
  vct(1:ifull,k) = bx(1:ifull)*newu(1:ifull) + by(1:ifull)*newv(1:ifull) + bz(1:ifull)*neww(1:ifull)
end do  ! k loop
  
return
end subroutine interpcurrent4

! *****************************************************************************
! FILE IO ROUTINES

! This version reads and interpolates a surface field
subroutine gethist1(vname,varout)

use cc_mpi             ! CC MPI routines
use darcdf_m           ! Netcdf data
use infile             ! Input file routines
use newmpar_m          ! Grid parameters
      
integer ier
real, dimension(ifull), intent(out) :: varout
real, dimension(fwsize) :: ucc
character(len=*), intent(in) :: vname
      
if ( iotest ) then
  ! read without interpolation or redistribution
  call histrd(iarchi,ier,vname,varout,ifull)
else
  ! for multiple input files
  call histrd(iarchi,ier,vname,ucc,6*ik*ik)
  call doints1(ucc, varout)
end if ! iotest

return
end subroutine gethist1

! This version reads, fills and interpolates a surface field
subroutine fillhist1(vname,varout,mask_a,fill_count)
      
use cc_mpi             ! CC MPI routines
use darcdf_m           ! Netcdf data
use infile             ! Input file routines
use newmpar_m          ! Grid parameters
      
integer, intent(inout) :: fill_count
integer ier
real, dimension(ifull), intent(out) :: varout
real, dimension(fwsize) :: ucc
logical, dimension(fwsize), intent(in) :: mask_a
character(len=*), intent(in) :: vname
      
if ( iotest ) then
  ! read without interpolation or redistribution
  call histrd(iarchi,ier,vname,varout,ifull)
else
  ! for multiple input files
  call histrd(iarchi,ier,vname,ucc,6*ik*ik)
  call fill_cc1(ucc,mask_a,fill_count)
  call doints1(ucc, varout)
end if ! iotest
      
return
end subroutine fillhist1

! This version reads, fills and interpolates a surface velocity field
subroutine fillhistuv1o(uname,vname,uarout,varout,mask_a,fill_count)
   
use cc_mpi             ! CC MPI routines
use darcdf_m           ! Netcdf data
use infile             ! Input file routines
use newmpar_m          ! Grid parameters
      
integer, intent(inout) :: fill_count
integer ier
real, dimension(ifull), intent(out) :: uarout, varout
real, dimension(fwsize) :: ucc, vcc
logical, dimension(fwsize), intent(in) :: mask_a
character(len=*), intent(in) :: uname, vname
      
if ( iotest ) then
  ! read without interpolation or redistribution
  call histrd(iarchi,ier,uname,uarout,ifull)
  call histrd(iarchi,ier,vname,varout,ifull)
else
  ! for multiple input files
  call histrd(iarchi,ier,uname,ucc,6*ik*ik)
  call histrd(iarchi,ier,vname,vcc,6*ik*ik)
  call interpcurrent1(uarout,varout,ucc,vcc,mask_a,fill_count)
end if ! iotest
      
return
end subroutine fillhistuv1o

subroutine fillhist1u(vname,uname,ifrac,mask_a,fill_count,filldefault,fillmode)

use cc_mpi                          ! CC MPI routines
use darcdf_m                        ! Netcdf data
use infile                          ! Input file routines
use newmpar_m                       ! Grid parameters
use uclem_ctrl, only : urbtemp,    &
    uclem_loadd                     ! Urban

integer, intent(inout) :: fill_count
integer, intent(in) :: ifrac, fillmode
integer ierr
real, intent(in) :: filldefault
real, dimension(ifull) :: varoutr4
real(kind=8), dimension(ifull) :: varout
logical, dimension(fwsize), intent(in) :: mask_a
character(len=*), intent(in) :: vname
character(len=*), intent(in) :: uname

if ( iotest ) then
  call histrd(iarchi,ierr,vname,varout,ifull) 
else
  call fillhist1(vname,varoutr4,mask_a,fill_count)
  varout = real(varoutr4,8)
end if
where ( varout>900._8 ) ! missing
  varout = real(filldefault,8)  
end where
select case(fillmode)
  case(0)
    ! do nothing  
  case(1)
    where ( varout>150._8 )  
      varout = min( max( varout, 170._8), 380._8 ) - real(urbtemp,8)
    end where
  case(2)
    where( varout<1.e-20_8 )
      varout = real(filldefault,8)
    end where
  case default
    write(6,*) "ERROR: Unknown fillmode in fillhist1u"
    call ccmpi_abort(-1)
end select
call uclem_loadd(varout,uname,ifrac,0)

return
end subroutine fillhist1u

! This version reads 3D fields
subroutine gethist4(vname,varout)

use cc_mpi             ! CC MPI routines
use darcdf_m           ! Netcdf data
use infile             ! Input file routines
use newmpar_m          ! Grid parameters
      
integer ier
real, dimension(:,:), intent(out) :: varout
real, dimension(fwsize,size(varout,2)) :: ucc
character(len=*), intent(in) :: vname

if ( iotest ) then
  ! read without interpolation or redistribution
  call histrd(iarchi,ier,vname,varout,ifull)
else
  ! for multiple input files
  call histrd(iarchi,ier,vname,ucc,6*ik*ik)
  call doints4(ucc,varout)
end if ! iotest
      
return
end subroutine gethist4   

! This version reads and interpolates 3D atmospheric fields
subroutine gethist4a(vname,varout,vmode,levkin,t_a_lev)
      
use cc_mpi               ! CC MPI routines
use darcdf_m             ! Netcdf data
use infile               ! Input file routines
use newmpar_m            ! Grid parameters
      
integer, intent(in) :: vmode
integer, intent(in), optional :: levkin
integer ier
real, dimension(:,:), intent(out) :: varout
real, dimension(fwsize), intent(out), optional :: t_a_lev
real, dimension(fwsize,kk) :: ucc
real, dimension(ifull,kk) :: u_k
character(len=*), intent(in) :: vname

if ( iotest ) then
  ! read without interpolation or redistribution
  call histrd(iarchi,ier,vname,u_k,ifull)
else
  ! for multiple input files
  call histrd(iarchi,ier,vname,ucc,6*ik*ik)
  if ( fwsize>0.and.present(levkin).and.present(t_a_lev) ) then
    if ( levkin<1 .or. levkin>kk ) then
      write(6,*) "ERROR: Invalid choice of levkin in gethist4a - ",trim(vname)
      call ccmpi_abort(-1)
    end if
    t_a_lev(1:fwsize) = ucc(1:fwsize,levkin)   ! store for psl calculation  
  end if
  call doints4(ucc, u_k)      
end if ! iotest

! vertical interpolation
call vertint(u_k,varout,vmode,sigin)
      
return
end subroutine gethist4a  

! This version reads and interpolates 3D atmospheric winds
subroutine gethistuv4a(uname,vname,uarout,varout,umode,vmode)

use cc_mpi             ! CC MPI routines
use darcdf_m           ! Netcdf data
use infile             ! Input file routines
use newmpar_m          ! Grid parameters

integer, intent(in) :: umode, vmode
integer ier
real, dimension(:,:), intent(out) :: uarout, varout
real, dimension(fwsize,kk) :: ucc, vcc
real, dimension(ifull,kk) :: u_k, v_k
character(len=*), intent(in) :: uname, vname

if ( iotest ) then
  ! read without interpolation or redistribution
  call histrd(iarchi,ier,uname,u_k,ifull)
  call histrd(iarchi,ier,vname,v_k,ifull)
else
  ! for multiple input files
  call histrd(iarchi,ier,uname,ucc,6*ik*ik)
  call histrd(iarchi,ier,vname,vcc,6*ik*ik)
  call interpwind4(u_k,v_k,ucc,vcc)
end if ! iotest

! vertical interpolation
call vertint(u_k,uarout,umode,sigin)
call vertint(v_k,varout,vmode,sigin)
      
return
end subroutine gethistuv4a  

! This version reads, fills a 3D field for the ocean
subroutine fillhist4(vname,varout,mask_a,fill_count)
  
use cc_mpi             ! CC MPI routines
use darcdf_m           ! Netcdf data
use infile             ! Input file routines
use newmpar_m          ! Grid parameters
      
integer, intent(inout) :: fill_count
integer ier
real, dimension(:,:), intent(out) :: varout
real, dimension(fwsize,size(varout,2)) :: ucc
logical, dimension(fwsize), intent(in) :: mask_a
character(len=*), intent(in) :: vname

if ( iotest ) then
  ! read without interpolation or redistribution
  call histrd(iarchi,ier,vname,varout,ifull)
else
  ! for multiple input files
  call histrd(iarchi,ier,vname,ucc,6*ik*ik)
  call fill_cc4(ucc,mask_a,fill_count)
  call doints4(ucc, varout)
end if ! iotest

return
end subroutine fillhist4

! This version reads, fills and interpolates a 3D field for the ocean
subroutine fillhist4o(vname,varout,mask_3d,fill_count,bath)
   
use cc_mpi             ! CC MPI routines
use darcdf_m           ! Netcdf data
use infile             ! Input file routines
use mlo_ctrl           ! Ocean physics control layer
use newmpar_m          ! Grid parameters
      
integer, intent(inout) :: fill_count
integer ier
real, dimension(:,:), intent(out) :: varout
real, dimension(fwsize,ok) :: ucc
real, dimension(ifull,ok) :: u_k
real, dimension(ifull), intent(in) :: bath
logical, dimension(fwsize,ok), intent(in) :: mask_3d
character(len=*), intent(in) :: vname

if ( iotest ) then
  ! read without interpolation or redistribution
  call histrd(iarchi,ier,vname,u_k,ifull)
else
  ! for multiple input files
  call histrd(iarchi,ier,vname,ucc,6*ik*ik)
  call fill_cc4(ucc,mask_3d,fill_count)
  call doints4(ucc,u_k)
end if ! iotest

! vertical interpolation
varout(:,:) = 0.5*(minval(u_k)+maxval(u_k)) ! easier for debuging
call mloregrid(ok,gosig_3,bath,u_k,varout,0) ! should use mode 3 or 4?

return
end subroutine fillhist4o

! This version reads, fills and interpolates ocean currents
subroutine fillhistuv4o(uname,vname,uarout,varout,mask_a,fill_count,bath)
  
use cc_mpi             ! CC MPI routines
use darcdf_m           ! Netcdf data
use infile             ! Input file routines
use mlo_ctrl           ! Ocean physics control layer
use newmpar_m          ! Grid parameters
      
integer, intent(inout) :: fill_count
integer ier
real, dimension(:,:), intent(out) :: uarout, varout
real, dimension(fwsize,ok) :: ucc, vcc
real, dimension(ifull,ok) :: u_k, v_k
real, dimension(ifull), intent(in) :: bath
logical, dimension(fwsize,ok), intent(in) :: mask_a
character(len=*), intent(in) :: uname, vname

if ( iotest ) then
  ! read without interpolation or redistribution
  call histrd(iarchi,ier,uname,u_k,ifull)
  call histrd(iarchi,ier,vname,v_k,ifull)
else
  ! for multiple input files
  call histrd(iarchi,ier,uname,ucc,6*ik*ik)
  call histrd(iarchi,ier,vname,vcc,6*ik*ik)
  call interpcurrent4(u_k,v_k,ucc,vcc,mask_a,fill_count)
end if ! iotest

! vertical interpolation
uarout = 0.
varout = 0.
call mloregrid(ok,gosig_3,bath,u_k,uarout,0)
call mloregrid(ok,gosig_3,bath,v_k,varout,0)

return
end subroutine fillhistuv4o

subroutine fillhist4u(vname,uname,ifrac,mask_a,fill_count,filldefault,fillmode)

use cc_mpi                          ! CC MPI routines
use darcdf_m                        ! Netcdf data
use infile                          ! Input file routines
use newmpar_m                       ! Grid parameters
use uclem_ctrl, only : urbtemp,    &
    uclem_loadd                     ! Urban

integer, intent(inout) :: fill_count
integer, intent(in) :: ifrac, fillmode
integer k, ierr
real, intent(in) :: filldefault
real, dimension(ifull,5) :: varoutr4
real(kind=8), dimension(ifull,5) :: varout
logical, dimension(fwsize), intent(in) :: mask_a
character(len=*), intent(in) :: vname
character(len=*), intent(in) :: uname
character(len=20) lname

if ( iotest ) then
  call histrd(iarchi,ierr,vname,varout,ifull) 
else
  call fillhist4(vname,varoutr4,mask_a,fill_count)
  varout = real(varoutr4,8)
end if
where ( varout>900._8 ) ! missing
  varout = real(filldefault,8)  
end where
select case(fillmode)
  case(0)
    ! do nothing  
  case(1)
    where ( varout>150._8 )  
      varout = min( max( varout, 170._8), 380._8 ) - real(urbtemp,8)
    end where
  case(2)
    where( varout<1.e-20_8 )
      varout = real(filldefault,8)
    end where
  case default
    write(6,*) "ERROR: Unknown fillmode in fillhist4u"
    call ccmpi_abort(-1)
end select
do k = 1,5
  write(lname,'(A,I1.1)') uname,k
  call uclem_loadd(varout(:,k),lname,ifrac,0)
end do
        
return
end subroutine fillhist4u

! *****************************************************************************
! FILE DATA MESSAGE PASSING ROUTINES

! Define mapping for distributing file data to processors
subroutine file_wininit

use cc_mpi             ! CC MPI routines
use infile             ! Input file routines
use newmpar_m          ! Grid parameters
use parm_m             ! Model configuration

integer n, ipf
integer mm, iq, idel, jdel
integer ncount, w
integer ncount_a, ncount_b
integer iin, jjn, no, ip
logical, dimension(0:fnproc-1) :: lfile
integer, dimension(nproc*fncount) :: tempmap_send, tempmap_smod
integer, dimension(nproc*fncount) :: tempmap_recv, tempmap_rmod
logical, dimension(0:nproc-1) :: lproc, lproc_t

if ( allocated(filemap_req) ) then
  write(6,*) "ERROR: Close input file before opening a new file"
  call ccmpi_abort(-1)
end if
if ( allocated(filemap_recv) ) then
  write(6,*) "ERROR: Close input file before opening a new file"
  call ccmpi_abort(-1)
end if
if ( allocated(filemap_send) ) then
  write(6,*) "ERROR: Close input file before opening a new file"
  call ccmpi_abort(-1)
end if
if ( allocated(axs_w) ) then
  deallocate( axs_w, ays_w, azs_w )
  deallocate( bxs_w, bys_w, bzs_w )
end if
if ( allocated(pijn2pw) ) then
  deallocate( pijn2pw, pijn2pn )  
end if

if ( myid==0 ) then
  write(6,*) "Create map for input file data"
end if

! calculate which grid points and input files are needed by this processor
lfile(:) = .false.
do mm = 1,m_fly
  do iq = 1,ifull
    idel = int(xg4(iq,mm))
    jdel = int(yg4(iq,mm))
    n = nface4(iq,mm)
    ! search stencil of bi-cubic interpolation
    lfile(procarray(idel,  jdel+2,n)) = .true.
    lfile(procarray(idel+1,jdel+2,n)) = .true.
    lfile(procarray(idel-1,jdel+1,n)) = .true.
    lfile(procarray(idel  ,jdel+1,n)) = .true.
    lfile(procarray(idel+1,jdel+1,n)) = .true.
    lfile(procarray(idel+2,jdel+1,n)) = .true.
    lfile(procarray(idel-1,jdel,  n)) = .true.
    lfile(procarray(idel  ,jdel,  n)) = .true.
    lfile(procarray(idel+1,jdel,  n)) = .true.
    lfile(procarray(idel+2,jdel,  n)) = .true.
    lfile(procarray(idel,  jdel-1,n)) = .true.
    lfile(procarray(idel+1,jdel-1,n)) = .true.   
  end do
end do

! Construct a map of files to be accessed by this process
if ( myid==0 ) then
  write(6,*) "-> Create map of files required by this process"
end if
ncount = count(lfile(0:fnproc-1))
allocate( filemap_req(ncount), filemap_qmod(ncount) )
ncount = 0
do w = 0,fnproc-1
  if ( lfile(w) ) then
    ncount = ncount + 1
    filemap_req(ncount) = mod( w, fnresid )
    filemap_qmod(ncount) = w/fnresid
  end if
end do

! Construct a map of processes that need this file
if ( myid==0 ) then
  write(6,*) "-> Create map for communication between processes"  
end if
tempmap_send(:) = -1
tempmap_smod(:) = -1
tempmap_recv(:) = -1
tempmap_rmod(:) = -1
allocate( filemap_indx(0:nproc-1,0:fncount-1) )
filemap_indx = 0
filemap_indxlen = 0
ncount_a = 0
ncount_b = 0
do ipf = 0,fncount-1
  lproc(:) = .false.
  do w = 1,size(filemap_req)
    if ( filemap_qmod(w) == ipf ) then
      lproc(filemap_req(w)) = .true.
    end if
  end do 
  call ccmpi_allreduce(lproc,lproc_t,"or",comm_node)
  ncount = 0
  do w = 0,nproc-1
    if ( lproc_t(w) ) then
      filemap_indxlen = filemap_indxlen + 1
      filemap_indx(w,ipf) = filemap_indxlen
      ncount = ncount + 1
      if ( mod(ncount,node_nproc)/=node_myid ) then
        lproc_t(w) = .false.
      end if
    end if
    if ( lproc_t(w) ) then
      ncount_b = ncount_b + 1
      tempmap_recv(ncount_b) = w
      tempmap_rmod(ncount_b) = ipf
    end if
  end do
  call ccmpi_alltoall(lproc_t,comm_world) ! global transpose
  do w = 0,nproc-1
    if ( lproc_t(w) ) then
      ncount_a = ncount_a + 1
      tempmap_send(ncount_a) = w
      tempmap_smod(ncount_a) = ipf
    end if
  end do  
end do
allocate( filemap_send(ncount_a), filemap_smod(ncount_a) )
allocate( filemap_recv(ncount_b), filemap_rmod(ncount_b) )
filemap_send(1:ncount_a) = tempmap_send(1:ncount_a)
filemap_smod(1:ncount_a) = tempmap_smod(1:ncount_a)
filemap_recv(1:ncount_b) = tempmap_recv(1:ncount_b)
filemap_rmod(1:ncount_b) = tempmap_rmod(1:ncount_b)

call ccmpi_filewininit(kblock)

! Define halo indices for ccmpi_filebounds
if ( myid==0 ) then
  write(6,*) "-> Setup bounds function for processors reading input files"  
end if
call ccmpi_filebounds_setup(comm_ip)

! Distribute fields for vector rotation
if ( myid==0 ) then
  write(6,*) "-> Distribute vector rotation data to processors reading input files"
end if
allocate(axs_w(fwsize), ays_w(fwsize), azs_w(fwsize))
allocate(bxs_w(fwsize), bys_w(fwsize), bzs_w(fwsize))
if ( myid==0 ) then
  call ccmpi_filedistribute(axs_w,axs_a,comm_ip)
  call ccmpi_filedistribute(ays_w,ays_a,comm_ip)
  call ccmpi_filedistribute(azs_w,azs_a,comm_ip)
  call ccmpi_filedistribute(bxs_w,bxs_a,comm_ip)
  call ccmpi_filedistribute(bys_w,bys_a,comm_ip)
  call ccmpi_filedistribute(bzs_w,bzs_a,comm_ip)
  deallocate( axs_a, ays_a, azs_a )
  deallocate( bxs_a, bys_a, bzs_a )
else if ( fwsize>0 ) then
  call ccmpi_filedistribute(axs_w,comm_ip)
  call ccmpi_filedistribute(ays_w,comm_ip)
  call ccmpi_filedistribute(azs_w,comm_ip)
  call ccmpi_filedistribute(bxs_w,comm_ip)
  call ccmpi_filedistribute(bys_w,comm_ip)
  call ccmpi_filedistribute(bzs_w,comm_ip)
end if

! Calculate indices for unpacking
if ( myid==0 ) then
  write(6,*) "-> Calculate unpacking indices"
end if
iin = (pil_g-1)/pipan
jjn = (pil_g-1)/pjpan
allocate( pijn2pw(0:iin,0:jjn,0:npanels), pijn2pn(0:iin,0:jjn,0:npanels) )
pijn2pw = -1 ! missing
pijn2pn = -1 ! missing
do w = 1,size(filemap_req)
   ip = filemap_req(w) + filemap_qmod(w)*fnresid
   do n = 1,pnpan
      no = n - pnoff(ip)
      idel = pioff(ip,no) + 1
      jdel = pjoff(ip,no) + 1
      iin = (idel-1)/pipan
      jjn = (jdel-1)/pjpan
      pijn2pw(iin,jjn,no) = w
      pijn2pn(iin,jjn,no) = n
   end do
end do

if ( myid==0 ) then
  write(6,*) "-> Finished creating control data for input file data"
end if

return
end subroutine file_wininit

subroutine onthefly_exit

use cc_mpi

integer i

#ifdef share_ifullg
do i = 1,xy4_count
  call ccmpi_freeshdata(xx4save_win(i))
  call ccmpi_freeshdata(yy4save_win(i))
end do
#endif

end subroutine onthefly_exit

end module onthefly_m
