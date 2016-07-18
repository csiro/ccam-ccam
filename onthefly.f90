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

! Main NetCDF input routines.  Host grid is automatically
! interpolated to nested model grid.  Three options are
!   nested=0  Initial conditions
!   nested=1  Nudging fields
!   nested=2  Surface data recycling
      
! This version supports the parallel file routines contained
! in infile.f90.  Hence, restart files do not require any
! gathers and scatters.

! In the case where the grid needs to be interpolated, RMA
! is used to distribute host data from processes, which
! reduces the amount of message passing.
    
! When -Dusempi3 is enabled, then host data arrays are
! shared between processes on a node.  The node captian
! is then responsible for obtaining interpolation data
! for all processes on a node.
    
! Thanks to Paul Ryan for optimising NetCDF routines
    
module onthefly_m
    
implicit none

private
public onthefly, retopo
    
integer, parameter :: nord = 3                                ! 1 for bilinear, 3 for bicubic interpolation
integer, save :: ik, jk, kk, ok, nsibx                        ! input grid size
integer dk, fwsize                                            ! size of temporary arrays
integer, dimension(:,:), allocatable, save :: nface4          ! interpolation panel index
integer, dimension(0:5), save :: comm_face                    ! communicator for processes requiring a input panel
real, save :: rlong0x, rlat0x, schmidtx                       ! input grid coordinates
real, dimension(3,3), save :: rotpoles, rotpole               ! vector rotation data
real, dimension(:,:), allocatable, save :: xg4, yg4           ! interpolation coordinate indices
real, dimension(:), allocatable, save :: axs_a, ays_a, azs_a  ! vector rotation data
real, dimension(:), allocatable, save :: bxs_a, bys_a, bzs_a  ! vector rotation data 
real, dimension(:), allocatable, save :: axs_w, ays_w, azs_w  ! vector rotation data
real, dimension(:), allocatable, save :: bxs_w, bys_w, bzs_w  ! vector rotation data
real, dimension(:), allocatable, save :: sigin                ! input vertical coordinates
logical iotest, newfile                                       ! tests for interpolation and new metadata
logical, dimension(0:5), save :: nfacereq = .false.           ! list of panels required for interpolation
logical, save :: bcst_allocated = .false.                     ! Bcast communicator groups have been defined

contains

! *****************************************************************************
! Main interface for input data that reads grid metadata
    
subroutine onthefly(nested,kdate_r,ktime_r,psl,zss,tss,sicedep,fracice,t,u,v,qg,tgg,wb,wbice,snowd,qfg, &
                    qlg,qrg,qsng,qgrg,tggsn,smass,ssdn,ssdnn,snage,isflag,mlodwn,ocndwn,xtgdwn)

use aerosolldr       ! LDR prognostic aerosols
use cc_mpi           ! CC MPI routines
use darcdf_m         ! Netcdf data
use infile           ! Input file routines
use mlo              ! Ocean physics and prognostic arrays
use newmpar_m        ! Grid parameters
use parm_m           ! Model configuration
use soil_m           ! Soil and surface data
use stime_m          ! File date data

implicit none

integer, parameter :: nihead = 54
integer, parameter :: nrhead = 14

integer, intent(in) :: nested
integer, intent(out) :: kdate_r, ktime_r
integer, save :: maxarchi
integer mtimer, ierx, idvkd, idvkt, idvmt, idvtime
integer kdate_rsav, ktime_rsav
integer, dimension(nihead) :: nahead
integer, dimension(ifull), intent(out) :: isflag
real timer
real, dimension(ifull,wlev,4), intent(out) :: mlodwn
real, dimension(ifull,kl,naero), intent(out) :: xtgdwn
real, dimension(ifull,ms), intent(out) :: wb, wbice, tgg
real, dimension(ifull,3), intent(out) :: tggsn, smass, ssdn
real, dimension(ifull,2), intent(out) :: ocndwn
real, dimension(:,:), intent(out) :: t, u, v, qg, qfg, qlg, qrg, qsng, qgrg
real, dimension(ifull), intent(out) :: psl, zss, tss, fracice, snowd
real, dimension(ifull), intent(out) :: sicedep, ssdnn, snage
real, dimension(nrhead) :: ahead
real, dimension(10) :: rdum
logical ltest, tst
character(len=80) datestring

call START_LOG(onthefly_begin)

if ( myid==0 ) write(6,*) 'Entering onthefly for nested,ktau = ',nested,ktau

!--------------------------------------------------------------------
! pfall indicates all processors have a parallel input file and there
! is no need to broadcast metadata (see infile.f90).  Otherwise read
! metadata on myid=0 and broadcast that data to all processors.
if ( myid==0 .or. pfall ) then
    
    ! Locate new file and read grid metadata --------------------------
  if ( ncid/=ncidold ) then
    if ( myid==0 ) write(6,*) 'Reading new file metadata'
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
      endif  ! (schmidtx<=0..or.schmidtx>1.)
    end if   ! ierx==0 ..else..
    call ccnf_inq_dimlen(ncid,'time',maxarchi)
    if ( myid==0 ) then
      write(6,*) "Found ik,jk,kk,ok ",ik,jk,kk,ok
      write(6,*) "      maxarchi ",maxarchi
      write(6,*) "      rlong0x,rlat0x,schmidtx ",rlong0x,rlat0x,schmidtx
    end if
  end if
  
  ! search for required date ----------------------------------------
  if ( myid==0 ) write(6,*)'Search for kdate_s,ktime_s >= ',kdate_s,ktime_s
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
    mtimer = nint(timer)
    call datefix(kdate_r,ktime_r,mtimer)
    ! ltest = .false. when correct date is found
    ltest = (2400*(kdate_r-kdate_s)-1200*nsemble+(ktime_r-ktime_s))<0
  end do
  if ( nsemble/=0 ) then
    kdate_r = kdate_s
    ktime_r = ktime_s
  end if
  if ( ltest ) then
    ! ran out of file before correct date was located
    ktime_r = -1
  end if
  if ( myid==0 ) then
    if ( ltest ) then
      write(6,*) 'Search failed with ltest,iarchi =',ltest, iarchi
      write(6,*) '                kdate_r,ktime_r =',kdate_r, ktime_r
    else
      write(6,*) 'Search was sucessful with ltest,iarchi =',ltest, iarchi
      write(6,*) '                       kdate_r,ktime_r =',kdate_r, ktime_r
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
  call ccmpi_bcast(rdum(1:10),0,comm_world)
  rlong0x  = rdum(1)
  rlat0x   = rdum(2)
  schmidtx = rdum(3)
  kdate_r  = nint(rdum(4))*10000+nint(rdum(5))*100+nint(rdum(6))
  ktime_r  = nint(rdum(7))
  newfile  = (nint(rdum(8))==1)
  iarchi   = nint(rdum(9))
  nsibx    = nint(rdum(10))
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
  if ( nested==2 ) then
    if ( myid==0 ) write(6,*) "WARN: Cannot locate date/time in input file"
    return
  else
    write(6,*) "ERROR: Cannot locate date/time in input file"
    call ccmpi_abort(-1)
  end if
end if
!--------------------------------------------------------------------
      
! Here we call ontheflyx with different automatic array
! sizes.  dk is used for global arrays that are defined
! on myid==0.  fwsize is used for MPI RMA in loading 
! files over multiple processors.
      
! Note that if histrd fails to find a variable, it returns
! zero in the output array
      
if ( myid==0 ) then
  dk = ik ! non-zero automatic array size in onthefly_work
else
  dk = 0  ! zero automatic array size in onthefly_work
end if

! memory needed to read input files
fwsize = pil*pjl*pnpan*mynproc 

call onthefly_work(nested,kdate_r,ktime_r,psl,zss,tss,sicedep,fracice,t,u,v,qg,tgg,wb,wbice, &
                   snowd,qfg,qlg,qrg,qsng,qgrg,tggsn,smass,ssdn,ssdnn,snage,isflag,mlodwn,   &
                   ocndwn,xtgdwn)

if ( myid==0 ) write(6,*) "Leaving onthefly"

call END_LOG(onthefly_end)

return
end subroutine onthefly


! *****************************************************************************
! Read data from netcdf file
      
! Input usually consists of either a single input file that is
! scattered across processes, or multiple input files that are
! read by many processes and shared by RMA.  In the case of
! restart files, then there is no need for message passing.
subroutine onthefly_work(nested,kdate_r,ktime_r,psl,zss,tss,sicedep,fracice,t,u,v,qg,tgg,wb,wbice, &
                         snowd,qfg,qlg,qrg,qsng,qgrg,tggsn,smass,ssdn,ssdnn,snage,isflag,mlodwn,   &
                         ocndwn,xtgdwn)
      
use aerosolldr, only : ssn,naero               ! LDR aerosol scheme
use ateb, only : atebdwn, urbtemp              ! Urban
use cable_def_types_mod, only : ncs, ncp       ! CABLE dimensions
use casadimension, only : mplant,mlitter,msoil ! CASA dimensions
use carbpools_m                                ! Carbon pools
use cc_mpi                                     ! CC MPI routines
use cfrac_m                                    ! Cloud fraction
use cloudmod                                   ! Prognostic strat cloud
use darcdf_m                                   ! Netcdf data
use extraout_m                                 ! Additional diagnostics      
use histave_m, only : cbas_ave,ctop_ave, &     ! Time average arrays
    wb_ave,tscr_ave
use infile                                     ! Input file routines
use latlong_m                                  ! Lat/lon coordinates
use latltoij_m                                 ! Lat/Lon to cubic ij conversion
use mlo, only : wlev,micdwn,mloregrid,wrtemp   ! Ocean physics and prognostic arrays
use mlodynamics                                ! Ocean dynamics
use mlodynamicsarrays_m                        ! Ocean dynamics data
use morepbl_m                                  ! Additional boundary layer diagnostics
use newmpar_m                                  ! Grid parameters
use nharrs_m, only : phi_nh,lrestart           ! Non-hydrostatic atmosphere arrays
use nsibd_m, only : isoilm,rsmin               ! Land-surface arrays
use parm_m                                     ! Model configuration
use parmdyn_m                                  ! Dynamics parmaters
use parmgeom_m                                 ! Coordinate data
use prec_m, only : precip,precc                ! Precipitation
use riverarrays_m                              ! River data
use savuvt_m                                   ! Saved dynamic arrays
use savuv1_m                                   ! Saved dynamic arrays
use screen_m                                   ! Screen level diagnostics
use setxyz_m                                   ! Define CCAM grid
use sigs_m                                     ! Atmosphere sigma levels
use soil_m                                     ! Soil and surface data
use soilv_m                                    ! Soil parameters
use stime_m                                    ! File date data
use tkeeps, only : tke,eps                     ! TKE-EPS boundary layer
use tracers_m                                  ! Tracer data
use utilities                                  ! Grid utilities
use vecsuv_m                                   ! Map to cartesian coordinates
use vvel_m, only : dpsldt,sdot                 ! Additional vertical velocity
use xarrs_m, only : pslx                       ! Saved dynamic arrays
use workglob_m                                 ! Additional grid interpolation
use work2_m                                    ! Diagnostic arrays

implicit none

include 'const_phys.h'                         ! Physical constants
include 'kuocom.h'                             ! Convection parameters

real, parameter :: iotol = 1.E-5      ! tolarance for iotest grid matching
      
integer, intent(in) :: nested, kdate_r, ktime_r
integer idv, isoil, nud_test
integer levk, levkin, ier, igas, nemi
integer i, j, k, n, mm, iq, numneg
integer, dimension(fwsize) :: isoilm_a
integer, dimension(ifull), intent(out) :: isflag
integer, dimension(7+3*ms) :: ierc
integer, dimension(4), save :: iers
real, dimension(ifull,wlev,4), intent(out) :: mlodwn
real, dimension(ifull,kl,naero), intent(out) :: xtgdwn
real, dimension(ifull,2), intent(out) :: ocndwn
real, dimension(ifull,ms), intent(out) :: wb, wbice, tgg
real, dimension(ifull,3), intent(out) :: tggsn, smass, ssdn
real, dimension(:,:), intent(out) :: t, u, v, qg, qfg, qlg, qrg, qsng, qgrg
real, dimension(ifull), intent(out) :: psl, zss, tss, fracice
real, dimension(ifull), intent(out) :: snowd, sicedep, ssdnn, snage
real, dimension(ifull) :: dum6, tss_l, tss_s, pmsl
real, dimension(fwsize) :: ucc
real, dimension(fwsize) :: fracice_a, sicedep_a
real, dimension(fwsize) :: tss_l_a, tss_s_a, tss_a
real, dimension(fwsize) :: t_a_lev, psl_a
real, dimension(:), allocatable, save :: zss_a, ocndep_l
real, dimension(kk+4) :: dumr
character(len=8) vname
character(len=3) trnum
logical, dimension(ms) :: tgg_found, wetfrac_found, wb_found
logical tsstest, tst
logical mixr_found, siced_found, fracice_found, soilt_found
logical u10_found, carbon_found
logical, dimension(:), allocatable, save :: land_a, sea_a

integer, dimension(3) :: shsize
integer, save :: xx4_win, yy4_win
real, dimension(:), allocatable :: wts_a  ! not used here or defined in call setxyz
real(kind=8), dimension(:,:), pointer :: xx4, yy4
real(kind=8), dimension(:,:), allocatable, target :: xx4_dummy, yy4_dummy
real(kind=8), dimension(:), allocatable, target :: z_a_dummy, x_a_dummy, y_a_dummy
real(kind=8), dimension(:), pointer :: z_a, x_a, y_a

! land-sea mask method (nemi=3 use soilt, nemi=2 use tgg, nemi=1 use zs)
nemi = 3
      
! test if retopo fields are required
if ( nud_p==0 .and. nud_t==0 .and. nud_q==0 ) then
  nud_test = 0
else
  nud_test = 1
end if
      
! Determine if interpolation is required
iotest = 6*ik*ik==ifull_g .and. abs(rlong0x-rlong0)<iotol .and. abs(rlat0x-rlat0)<iotol .and. &
         abs(schmidtx-schmidt)<iotol .and. nsib==nsibx
if ( iotest ) then
  io_in = 1   ! no interpolation
  if ( myid==0 ) write(6,*) "Interpolation is required with iotest,io_in =",iotest, io_in
else
  io_in = -1  ! interpolation
  if ( myid==0 ) write(6,*) "Interpolation is not required with iotest,io_in =",iotest, io_in
end if

!--------------------------------------------------------------------
! Allocate interpolation, vertical level and mask arrays
! dk is only non-zero on myid==0
if ( .not.allocated(nface4) ) then
  allocate( nface4(ifull,4), xg4(ifull,4), yg4(ifull,4) )
end if
if ( newfile ) then
  if ( allocated(sigin) ) then
    deallocate( sigin, land_a, sea_a )
    deallocate( axs_a, ays_a, azs_a )
    deallocate( bxs_a, bys_a, bzs_a )          
  end if
  allocate( sigin(kk), land_a(fwsize), sea_a(fwsize) )
  allocate( axs_a(dk*dk*6), ays_a(dk*dk*6), azs_a(dk*dk*6) )
  allocate( bxs_a(dk*dk*6), bys_a(dk*dk*6), bzs_a(dk*dk*6) )
end if
      
!--------------------------------------------------------------------
! Determine input grid coordinates and interpolation arrays
if ( newfile .and. .not.iotest ) then
#ifdef usempi3
  shsize(1) = 1 + 4*ik
  shsize(2) = 1 + 4*ik
  call ccmpi_allocshdatar8(xx4,shsize(1:2),xx4_win)
  call ccmpi_allocshdatar8(yy4,shsize(1:2),yy4_win)
#else
  allocate( xx4_dummy(1+4*ik,1+4*ik), yy4_dummy(1+4*ik,1+4*ik) )
  xx4 => xx4_dummy
  yy4 => yy4_dummy
#endif

  if ( m_fly==1 ) then
    rlong4_l(:,1) = rlongg(:)*180./pi
    rlat4_l(:,1)  = rlatt(:)*180./pi
  end if
          
  if ( myid==0 ) then
    write(6,*) "Defining input file grid"
!   following setxyz call is for source data geom    ****   
    do iq = 1,ik*ik*6
      axs_a(iq) = real(iq)
      ays_a(iq) = real(iq)
      azs_a(iq) = real(iq)
    end do 
    allocate(x_a_dummy(ik*ik*6),y_a_dummy(ik*ik*6),z_a_dummy(ik*ik*6))
    allocate(wts_a(ik*ik*6))
    x_a => x_a_dummy
    y_a => y_a_dummy
    z_a => z_a_dummy
    call setxyz(ik,rlong0x,rlat0x,-schmidtx,x_a,y_a,z_a,wts_a,axs_a,ays_a,azs_a,bxs_a,bys_a,bzs_a,xx4,yy4)
    nullify(x_a, y_a, z_a)
    deallocate(x_a_dummy,y_a_dummy,z_a_dummy)
    deallocate(wts_a)
  end if ! (myid==0)
  
#ifdef usempi3
  call ccmpi_shepoch(xx4_win) ! also yy4_win
  if ( node_myid==0 ) then
    call ccmpi_bcastr8(xx4,0,comm_nodecaptian)
    call ccmpi_bcastr8(yy4,0,comm_nodecaptian)
  end if
  call ccmpi_shepoch(xx4_win) ! also yy4_win
#else
  call ccmpi_bcastr8(xx4,0,comm_world)
  call ccmpi_bcastr8(yy4,0,comm_world)
#endif
  
  ! calculate the rotated coords for host and model grid
  rotpoles = calc_rotpole(rlong0x,rlat0x)
  rotpole  = calc_rotpole(rlong0,rlat0)
  if ( myid==0 ) then
    write(6,*)'m_fly,nord ',m_fly,nord
    write(6,*)'kdate_r,ktime_r,ktau,ds',kdate_r,ktime_r,ktau,ds
    write(6,*)'rotpoles:'
    do i = 1,3
      write(6,'(3x,2i1,5x,2i1,5x,2i1,5x,3f8.4)') (i,j,j=1,3),(rotpoles(i,j),j=1,3)
    enddo
    if ( nmaxpr==1 ) then
      write(6,*)'in onthefly rotpole:'
      do i = 1,3
        write(6,'(3x,2i1,5x,2i1,5x,2i1,5x,3f8.4)') (i,j,j=1,3),(rotpole(i,j),j=1,3)
      enddo
      write(6,*)'xx4,yy4 ',xx4(id,jd),yy4(id,jd)
      write(6,*)'before latltoij for id,jd: ',id,jd
      write(6,*)'rlong0x,rlat0x,schmidtx ',rlong0x,rlat0x,schmidtx
    end if                ! (nmaxpr==1)
  end if                  ! (myid==0)

  ! setup interpolation arrays
  do mm = 1,m_fly  !  was 4, now may be set to 1 in namelist
    do iq = 1,ifull
      call latltoij(rlong4_l(iq,mm),rlat4_l(iq,mm),       & !input
                    rlong0x,rlat0x,schmidtx,              & !input
                    xg4(iq,mm),yg4(iq,mm),nface4(iq,mm),  & !output (source)
                    xx4,yy4,ik)
    end do
  end do

  nullify( xx4, yy4 )
#ifdef usempi3
  call ccmpi_freeshdata(xx4_win)
  call ccmpi_freeshdata(yy4_win)
#else
  deallocate( xx4_dummy, yy4_dummy )  
#endif

  ! Identify cubic panels to be processed
  if ( myid==0 ) then
    nfacereq(:) = .true. ! this is the host processor for bcast
  else
    nfacereq(:) = .false.
    do n = 0,npanels
      nfacereq(n) = any( nface4(:,:)==n )
    end do
  end if
  
  ! Define filemap for MPI RMA method
  call file_wininit
  
  ! Define comm_face for MPI IBcast method
  call splitface
       
end if ! newfile .and. .not.iotest
      
! -------------------------------------------------------------------
! read time invariant data when file is first opened
! need global zss_a for (potentially) landsea mask and psl interpolation
! need global isoilm_a for (potentially) landsea mask
if ( newfile ) then

  if ( myid==0 ) write(6,*) "Reading time invariant fields"  
    
  ! read vertical levels and missing data checks
  if ( myid==0 .or. pfall ) then
    call ccnf_inq_varid(ncid,'lev',idv,tst)
    if ( tst ) call ccnf_inq_varid(ncid,'layer',idv,tst)
    if ( tst ) call ccnf_inq_varid(ncid,'sigma',idv,tst)
    if ( tst ) then
      if ( myid==0 ) write(6,*) "No sigma level data found in input file"
      if ( kk>1 ) then
        write(6,*) "ERORR: multiple levels expected but no sigma data found ",kk
        call ccmpi_abort(-1)
      end if
      sigin(:) = 1.
    else
      call ccnf_get_vara(ncid,idv,1,kk,sigin)
      if ( myid==0 ) write(6,'(" sigin=",(9f7.4))') (sigin(k),k=1,kk)
    end if
    ! check for missing data
    iers(1:4) = 0
    call ccnf_inq_varid(ncid,'mixr',idv,tst)
    if ( tst ) iers(1) = -1
    call ccnf_inq_varid(ncid,'siced',idv,tst)
    if ( tst ) iers(2) = -1
    call ccnf_inq_varid(ncid,'fracice',idv,tst)
    if ( tst ) iers(3) = -1
    call ccnf_inq_varid(ncid,'soilt',idv,tst)
    if ( tst ) iers(4) = -1
  end if
  
  ! bcast data to all processors unless all processes are reading input files
  if ( .not.pfall ) then
    dumr(1:kk)      = sigin(1:kk)
    dumr(kk+1:kk+4) = real(iers(1:4))
    call ccmpi_bcast(dumr(1:kk+4),0,comm_world)
    sigin(1:kk) = dumr(1:kk)
    iers(1:4)   = nint(dumr(kk+1:kk+4))
  end if
  
  mixr_found    = (iers(1)==0)
  siced_found   = (iers(2)==0)
  fracice_found = (iers(3)==0)
  soilt_found   = (iers(4)==0)

  ! determine whether surface temperature needs to be interpolated (tsstest=.false.)
  tsstest = siced_found .and. fracice_found .and. iotest
  if ( myid==0 ) then
    if ( tsstest ) then
      write(6,*) "Surface temperature does not require interpolation"
      write(6,*) "tsstest,siced_found,fracice_found,iotest =",tsstest,siced_found,fracice_found,iotest
    else
      write(6,*) "Surface temperature requires interpolation"
      write(6,*) "tsstest,siced_found,fracice_found,iotest =",tsstest,siced_found,fracice_found,iotest
    end if
  end if
  if ( allocated(zss_a) ) deallocate(zss_a)
  if ( tsstest ) then
    ! load local surface temperature
    allocate( zss_a(ifull) )
    call histrd1(iarchi,ier,'zht',ik,zss_a,ifull)
  else if ( fnresid==1 ) then
    ! load global surface temperature using gather
    allocate( zss_a(6*dk*dk) )
    call histrd1(iarchi,ier,'zht',  ik,zss_a,6*ik*ik,nogather=.false.)
    call histrd1(iarchi,ier,'soilt',ik,ucc  ,6*ik*ik,nogather=.false.)
    if ( myid==0 ) then
      isoilm_a(:) = nint(ucc(:))
      if ( .not.soilt_found ) isoilm_a(:) = -1 ! missing value flag
    end if
  else
    ! load global surface temperature using RMA
    allocate( zss_a(fwsize) )
    call histrd1(iarchi,ier,'zht',  ik,zss_a,6*ik*ik,nogather=.true.)
    call histrd1(iarchi,ier,'soilt',ik,ucc  ,6*ik*ik,nogather=.true.)
    if ( fwsize>0 ) then
      isoilm_a(:) = nint(ucc(:))
      if ( .not.soilt_found ) isoilm_a(:) = -1 ! missing value flag
    end if
  end if
  
  ! read host ocean bathymetry data
  if ( nmlo/=0 .and. abs(nmlo)<=9 ) then
    if ( .not.allocated(ocndep_l) ) allocate(ocndep_l(ifull))
    call gethist1('ocndepth',ocndep_l)
  end if
  
  if ( myid==0 ) write(6,*) "Finished reading invariant fields"
  
else
    
  ! use saved metadata  
  mixr_found    = (iers(1)==0)
  siced_found   = (iers(2)==0)
  fracice_found = (iers(3)==0)
  soilt_found   = (iers(4)==0)
  tsstest = siced_found .and. fracice_found .and. iotest
  
end if ! newfile ..else..

! -------------------------------------------------------------------
! detemine the reference level below sig=0.9 (used to calculate psl)
levk = 0
levkin = 0
if ( nested==0 .or. ( nested==1.and.nud_test/=0 ) ) then
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

!--------------------------------------------------------------------
! Read surface pressure
! psf read when nested=0 or nested=1.and.nud_p/=0
psl_a(:) = 0.
psl(:)   = 0.
if ( nested==0 .or. ( nested==1 .and. nud_test/=0 ) ) then
  if ( iotest ) then
    call histrd1(iarchi,ier,'psf',ik,psl,ifull)
  else if ( fnresid==1 ) then
    call histrd1(iarchi,ier,'psf',ik,psl_a,6*ik*ik,nogather=.false.)
  else
    call histrd1(iarchi,ier,'psf',ik,psl_a,6*ik*ik,nogather=.true.)
  end if
endif

! -------------------------------------------------------------------
! Read surface temperature 
! read global tss to diagnose sea-ice or land-sea mask
if ( tsstest ) then
  call histrd1(iarchi,ier,'tsu',ik,tss,ifull)
  zss(:) = zss_a(:) ! use saved zss arrays
else
  if ( fnresid==1 ) then
    call histrd1(iarchi,ier,'tsu',ik,tss_a,6*ik*ik,nogather=.false.)
  else
    call histrd1(iarchi,ier,'tsu',ik,tss_a,6*ik*ik,nogather=.true.)
  end if
      
  ! set up land-sea mask from either soilt, tss or zss
  if ( newfile .and. fwsize>0 ) then
    if ( nemi==3 ) then 
      land_a(:) = isoilm_a(:)>0
      numneg = count( .not.land_a(:) )
      if ( any(isoilm_a(:)<0) ) nemi = 2
    end if ! (nemi==3)
    if ( nemi==2 ) then
      numneg = 0
      do iq = 1,fwsize
        if ( tss_a(iq)>0. ) then ! over land
          land_a(iq) = .true.
        else                     ! over sea
          land_a(iq) = .false.
          numneg = numneg + 1
        end if               ! (tss(iq)>0) .. else ..
      end do
      if ( numneg==0 ) nemi = 1  ! should be using zss in that case
    end if !  (nemi==2)
    tss_a(:) = abs(tss_a(:))
    if ( nemi==1 ) then
      land_a(:) = zss_a(:)>0.
      numneg = count(.not.land_a)
    end if ! (nemi==1)
    if ( myid==0 ) then
      write(6,*)'Land-sea mask using nemi = ',nemi
    end if
    sea_a(:) = .not.land_a(:)
  end if ! (newfile.and.fwsize>0)
end if ! (tsstest) ..else..

      
!--------------------------------------------------------------
! Read ocean data for nudging (sea-ice is read below)
! read when nested=0 or nested==1.and.nud/=0 or nested=2
if ( abs(nmlo)>=1 .and. abs(nmlo)<=9 ) then
  ! fixed ocean depth
  ocndwn(:,1) = ocndep_l(:)
  ! ocean potential temperature
  ! ocean temperature and soil temperature use the same arrays
  ! as no fractional land or sea cover is allowed in CCAM
  if ( ( nested/=1 .or. nud_sst/=0 ) .and. ok>0 ) then
    call fillhist4o('tgg',mlodwn(:,:,1),land_a,ocndwn(:,1))
    where ( mlodwn(:,:,1)>100. )
      mlodwn(:,:,1) = mlodwn(:,:,1) - wrtemp ! remove temperature offset for precision
    end where
  else
    mlodwn(:,:,1) = 293. - wrtemp
  end if ! (nestesd/=1.or.nud_sst/=0) ..else..
  ! ocean salinity
  if ( ( nested/=1 .or. nud_sss/=0 ) .and. ok>0 ) then
    call fillhist4o('sal',mlodwn(:,:,2),land_a,ocndwn(:,1))
    mlodwn(:,:,2) = max( mlodwn(:,:,2), 0. )
  else
    mlodwn(:,:,2) = 34.72   
  end if ! (nestesd/=1.or.nud_sss/=0) ..else..
  ! ocean currents
  if ( ( nested/=1 .or. nud_ouv/=0 ) .and. ok>0 ) then
    call fillhistuv4o('uoc','voc',mlodwn(:,:,3),mlodwn(:,:,4),land_a,ocndwn(:,1))
  else
    mlodwn(:,:,3:4) = 0.               
  end if ! (nestesd/=1.or.nud_ouv/=0) ..else..
  ! water surface height
  if ( nested/=1 .or. nud_sfh/=0 ) then
    call fillhist1('ocheight',ocndwn(:,2),land_a)
  else
    ocndwn(:,2) = 0.
  end if ! (nested/=1.or.nud_sfh/=0) ..else..
end if
!--------------------------------------------------------------


!--------------------------------------------------------------
! read sea ice here for prescribed SSTs configuration and for
! mixed-layer-ocean
if ( tsstest ) then
  call histrd1(iarchi,ier,'siced',  ik,sicedep,ifull)
  call histrd1(iarchi,ier,'fracice',ik,fracice,ifull)
  if ( any(fracice>1.) ) then
    write(6,*) "ERROR: Invalid fracice in input file"
    write(6,*) "Fracice should be between 0 and 1"
    write(6,*) "maximum fracice ",maxval(fracice)
    call ccmpi_abort(-1)
  end if
else
  if ( fnresid==1 ) then
    call histrd1(iarchi,ier,'siced',  ik,sicedep_a,6*ik*ik,nogather=.false.)
    call histrd1(iarchi,ier,'fracice',ik,fracice_a,6*ik*ik,nogather=.false.)
  else
    call histrd1(iarchi,ier,'siced',  ik,sicedep_a,6*ik*ik,nogather=.true.)
    call histrd1(iarchi,ier,'fracice',ik,fracice_a,6*ik*ik,nogather=.true.)
  end if
  if ( myid<fnresid ) then
    if ( any(fracice_a>1.) ) then
      write(6,*) "ERROR: Invalid fracice in input file"
      write(6,*) "Fracice should be between 0 and 1"
      write(6,*) "maximum fracice ",maxval(fracice_a)
      call ccmpi_abort(-1)
    end if
  end if
        
  ! diagnose sea-ice if required
  if ( fwsize>0 ) then
    if ( siced_found ) then  ! i.e. sicedep read in 
      if ( .not.fracice_found ) then ! i.e. sicedep read in; fracice not read in
        where ( sicedep_a(:)>0. )
          fracice_a(:) = 1.
        end where
      end if  ! (ierr/=0)  fracice
    else      ! sicedep not read in
      if ( .not.fracice_found ) then  ! neither sicedep nor fracice read in
        sicedep_a(:) = 0.  ! Oct 08
        fracice_a(:) = 0.
        if ( myid==0 ) then
          write(6,*) 'pre-setting siced in onthefly from tss'
        end if
        where ( abs(tss_a(:))<=271.6 ) ! for ERA-Interim
          sicedep_a(:) = 1.  ! Oct 08   ! previously 271.2
          fracice_a(:) = 1.
        end where
      else  ! i.e. only fracice read in;  done in indata, nestin
            ! but needed here for onthefly (different dims) 28/8/08        
        where ( fracice_a(:)>0.01 )
          sicedep_a(:) = 2.
        elsewhere
          sicedep_a(:) = 0.
          fracice_a(:) = 0.
        end where
      end if  ! .not.fracice_found ..else..
    end if    ! siced_found .. else ..    for sicedep

    ! fill surface temperature and sea-ice
    tss_l_a(:) = abs(tss_a(:))
    tss_s_a(:) = abs(tss_a(:))
    if ( fnresid==1 ) then
      call fill_cc1_gather(tss_l_a,sea_a)
      call fill_cc1_gather(tss_s_a,land_a)
      call fill_cc1_gather(sicedep_a,land_a)
      call fill_cc1_gather(fracice_a,land_a)
    else
      call fill_cc1_nogather(tss_l_a,sea_a)
      call fill_cc1_nogather(tss_s_a,land_a)
      call fill_cc1_nogather(sicedep_a,land_a)
      call fill_cc1_nogather(fracice_a,land_a)
    end if
  end if ! myid==0

  if ( iotest ) then
    ! This case occurs for missing sea-ice data
    if ( myid==0 ) then
      call ccmpi_distribute(zss,zss_a)
      call ccmpi_distribute(tss_l,tss_l_a)
      call ccmpi_distribute(tss_s,tss_s_a)
      call ccmpi_distribute(sicedep,sicedep_a)
      call ccmpi_distribute(fracice,fracice_a)
    else
      call ccmpi_distribute(zss)
      call ccmpi_distribute(tss_l)
      call ccmpi_distribute(tss_s)
      call ccmpi_distribute(sicedep)
      call ccmpi_distribute(fracice)
    end if
!   incorporate other target land mask effects
    where ( land(1:ifull) )
      sicedep = 0.
      fracice = 0.
      tss = tss_l
    elsewhere
      tss = tss_s
    end where
  else
!   The routine doints1 does the gather, calls ints4 and redistributes
    if ( fnresid==1 ) then
      call doints1_gather(zss_a,     zss)
      call doints1_gather(tss_l_a,   tss_l)
      call doints1_gather(tss_s_a,   tss_s)
      call doints1_gather(fracice_a, fracice)
      call doints1_gather(sicedep_a, sicedep)
    else
      call doints1_nogather(zss_a,     zss)
      call doints1_nogather(tss_l_a,   tss_l)
      call doints1_nogather(tss_s_a,   tss_s)
      call doints1_nogather(fracice_a, fracice)
      call doints1_nogather(sicedep_a, sicedep)
    end if
!   incorporate other target land mask effects
    where ( land(1:ifull) )
      tss(1:ifull) = tss_l(1:ifull)
    elsewhere
      tss(1:ifull) = tss_s(1:ifull)   ! no sign switch in CCAM
    end where
    where ( land(1:ifull) .or. sicedep(1:ifull)<0.05 ) ! for sflux
      sicedep(1:ifull) = 0.
      fracice(1:ifull) = 0.
    end where
  end if ! iotest
end if ! (tsstest) ..else..

! to be depeciated !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!if (nspecial==44.or.nspecial==46) then
!  do iq=1,ifull
!    rlongd=rlongg(iq)*180./pi
!    rlatd=rlatt(iq)*180./pi
!    if (rlatd>=-43..and.rlatd<=-30.) then
!      if (rlongd>=155..and.rlongd<=170.) then
!        tss(iq)=tss(iq)+1.
!      end if
!    end if
!  end do
!end if
!if (nspecial==45.or.nspecial==46) then
!  do iq=1,ifull
!    rlongd=rlongg(iq)*180./pi
!    rlatd=rlatt(iq)*180./pi
!    if (rlatd>=-15..and.rlatd<=-5.) then
!      if (rlongd>=150..and.rlongd<=170.) then
!        tss(iq)=tss(iq)+1.
!      end if
!    end if
!  end do
!end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! -------------------------------------------------------------------
! read atmospheric fields for nested=0 or nested=1.and.nud/=0

! air temperature
! read for nested=0 or nested=1.and.(nud_t/=0.or.nud_p/=0)
if ( nested==0 .or. ( nested==1.and.nud_test/=0 ) ) then
  call gethist4a('temp',t,2,levkin=levkin,t_a_lev=t_a_lev)
else
  t(1:ifull,:) = 300.    
end if ! (nested==0.or.(nested==1.and.nud_test/=0))

! winds
! read for nested=0 or nested=1.and.nud_uv/=0
if ( nested==0 .or. ( nested==1.and.nud_uv/=0 ) ) then
  call gethistuv4a('u','v',u,v,3,4)
else
  u(1:ifull,:) = 0.
  v(1:ifull,:) = 0.
end if ! (nested==0.or.(nested==1.and.nud_uv/=0))

! mixing ratio
! read for nested=0 or nested=1.and.nud_q/=0
if ( nested==0 .or. ( nested==1.and.nud_q/=0 ) ) then
  if ( mixr_found ) then
    call gethist4a('mixr',qg,2)      !     mixing ratio
  else
    call gethist4a('q',qg,2)         !     mixing ratio
  end if
else
  qg(1:ifull,:) = qgmin
end if ! (nested==0.or.(nested==1.and.nud_q/=0))

!------------------------------------------------------------
! Aerosol data
if ( abs(iaero)>=2 .and. ( nested/=1.or.nud_aero/=0 ) ) then
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
end if

!------------------------------------------------------------
! re-grid surface pressure by mapping to MSLP, interpolating and then map to surface pressure
! requires psl_a, zss, zss_a, t and t_a_lev
if ( nested==0 .or. ( nested==1.and.nud_test/=0 ) ) then
  if ( .not.iotest ) then
    if ( fwsize>0 ) then
      ! ucc holds pmsl_a
      call mslpx(ucc,psl_a,zss_a,t_a_lev,sigin(levkin))  ! needs pmsl (preferred)
    end if
    if ( fnresid==1 ) then
      call doints1_gather(ucc,pmsl)
    else
      call doints1_nogather(ucc,pmsl)
    end if
    ! invert pmsl to get psl
    call to_pslx(pmsl,psl,zss,t(:,levk),levk)  ! on target grid
  end if ! .not.iotest
end if


!**************************************************************
! This is the end of reading the nudging arrays
!**************************************************************


!--------------------------------------------------------------
! The following data is only read for initial conditions
if ( nested/=1 ) then

  ierc(:) = 0  ! flag for located variables
    
  !------------------------------------------------------------------
  ! check soil variables
  if ( myid==0 .or. pfall ) then
    if ( ccycle==0 ) then
      !call ccnf_inq_varid(ncid,'cplant1',idv,tst)
      !if ( tst ) ierc(7)=-1
      ierc(7) = 0
    else
      call ccnf_inq_varid(ncid,'glai',idv,tst)
      if ( .not.tst ) ierc(7) = 1
    end if
    do k = 1,ms
      write(vname,'("tgg",I1.1)') k
      call ccnf_inq_varid(ncid,vname,idv,tst)
      if ( .not.tst ) ierc(7+k) = 1
      write(vname,'("wetfrac",I1.1)') k
      call ccnf_inq_varid(ncid,vname,idv,tst)
      if ( .not.tst ) ierc(7+ms+k) = 1
      write(vname,'("wb",I1.1)') k
      call ccnf_inq_varid(ncid,vname,idv,tst)
      if ( .not.tst ) ierc(7+2*ms+k) = 1
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
          call ccnf_get_vara(ncid,idv,iarchi,ierc(3))
        end if
        call ccnf_inq_varid(ncid,'nstagu',idv,tst)
        if ( tst ) then
          lrestart = .false.
        else 
          call ccnf_get_vara(ncid,idv,iarchi,ierc(4))
        end if
        call ccnf_inq_varid(ncid,'nstagoff',idv,tst)
        if ( tst ) then
          lrestart = .false.
        else 
          call ccnf_get_vara(ncid,idv,iarchi,ierc(5))
        end if
        if ( abs(nmlo)>=3 .and. abs(nmlo)<=9 ) then
          if ( ok==wlev ) then
            call ccnf_inq_varid(ncid,'oldu101',idv,tst)
            if ( tst ) lrestart = .false.
            call ccnf_inq_varid(ncid,'oldv101',idv,tst)
            if ( tst ) lrestart = .false.
            call ccnf_inq_varid(ncid,'oldu201',idv,tst)
            if ( tst ) lrestart = .false.
            call ccnf_inq_varid(ncid,'oldv201',idv,tst)
            if ( tst ) lrestart = .false.
            call ccnf_inq_varid(ncid,'ipice',idv,tst)
            if ( tst ) lrestart = .false.
            call ccnf_inq_varid(ncid,'nstagoffmlo',idv,tst)
            if ( tst ) then
              lrestart = .false.
            else
              call ccnf_get_vara(ncid,idv,iarchi,ierc(6))
            end if
          else
            lrestart = .false.
          end if
        end if
      else
        lrestart = .false.
      end if ! kk=kl .and. iotest
      ierc(1:2) = 0
      if ( lrestart ) ierc(1) = 1
      call ccnf_inq_varid(ncid,'u10',idv,tst)
      if ( .not.tst ) ierc(2) = 1
    end if ! myid==0 .or. pfall
    
  end if   ! nested==0  
    
  if ( .not.pfall ) then
    call ccmpi_bcast(ierc(1:7+3*ms),0,comm_world)
  end if
  
  lrestart  = (ierc(1)==1)
  u10_found = (ierc(2)==1)
  if ( lrestart ) then
    nstag       = ierc(3)
    nstagu      = ierc(4)
    nstagoff    = ierc(5)
    nstagoffmlo = ierc(6)
    if ( myid==0 ) then
      write(6,*) "Continue staggering from"
      write(6,*) "nstag,nstagu,nstagoff ",nstag,nstagu,nstagoff
      if ( abs(nmlo)>=3 .and. abs(nmlo)<=9 ) then
        write(6,*) "nstagoffmlo ",nstagoffmlo
      end if
    end if
  end if
  carbon_found        = (ierc(7)==1)
  tgg_found(1:ms)     = (ierc(8:7+ms)==1)
  wetfrac_found(1:ms) = (ierc(8+ms:7+2*ms)==1)
  wb_found(1:ms)      = (ierc(8+2*ms:7+3*ms)==1)
        
  !------------------------------------------------------------------
  ! Read basic fields
  if ( nested==0 ) then
    if ( nsib==6 .or. nsib==7 ) then
      call gethist1('rs',rsmin)  
      call gethist1('zolnd',zo)
    end if
    call gethist1('rnd',precip)
    precip(:) = precip(:)/real(nperday)
    call gethist1('rnc',precc)
    precc(:) = precc(:)/real(nperday)
  end if
  
  !------------------------------------------------------------------
  ! Read snow and soil tempertaure
  call gethist1('snd',snowd)
  where ( .not.land(1:ifull) .and. (sicedep==0. .or. nmlo==0) )
    snowd = 0.
  end where
  if ( all(tgg_found(1:ms)) ) then
    call fillhist4('tgg',tgg,ms,sea_a)
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
          tgg(:,k) = tss(:)
        else
          call histrd1(iarchi,ier,vname,ik,tgg(:,k),ifull)
        end if
      else if ( fnresid==1 ) then
        if ( k==1 .and. .not.tgg_found(1) ) then
          ucc(1:dk*dk*6) = tss_a(1:dk*dk*6)
        else
          call histrd1(iarchi,ier,vname,ik,ucc,6*ik*ik,nogather=.false.)
        end if
        call fill_cc1_gather(ucc,sea_a)
        call doints1_gather(ucc,tgg(:,k))
      else
        if ( k==1 .and. .not.tgg_found(1) ) then
          ucc(1:fwsize) = tss_a(1:fwsize)
        else
          call histrd1(iarchi,ier,vname,ik,ucc,6*ik*ik,nogather=.true.)
        end if
        call fill_cc1_nogather(ucc,sea_a)
        call doints1_nogather(ucc,tgg(:,k))
      end if
    end do
  end if
  do k = 1,ms
    where ( tgg(:,k)<100. )
      tgg(:,k) = tgg(:,k) + wrtemp ! adjust range of soil temp for compressed history file
    end where
  end do  
  if ( .not.iotest ) then
    where ( snowd>0. .and. land(1:ifull) )
      tgg(:,1) = min( tgg(:,1), 270.1 )
    endwhere
  end if

  !--------------------------------------------------
  ! Read MLO sea-ice data
  if ( abs(nmlo)>=1 .and. abs(nmlo)<=9 ) then
    if ( .not.allocated(micdwn) ) allocate( micdwn(ifull,11) )
    call fillhist4('tggsn',micdwn(:,1:4),4,land_a)
    if ( all(micdwn(:,1)==0.) ) micdwn(:,1:4) = 270.
    micdwn(:,5) = fracice ! read above with nudging arrays
    micdwn(:,6) = sicedep ! read above with nudging arrays
    micdwn(:,7) = snowd*1.e-3
    call fillhist1('sto',micdwn(:,8),land_a)
    call fillhistuv1o('uic','vic',micdwn(:,9),micdwn(:,10),land_a)
    call fillhist1('icesal',micdwn(:,11),land_a)
  end if
  
  !------------------------------------------------------------------
  ! Read river data
  if ( (abs(nmlo)>=2.and.abs(nmlo)<=9) .or. nriver==1 ) then
    call gethist1('swater',watbdy)
  end if

  !------------------------------------------------------------------
  ! Read soil moisture
  wb(:,:) = 20.5
  if ( all(wetfrac_found(1:ms)) ) then
    call fillhist4('wetfrac',wb,ms,sea_a)
    wb(:,:) = wb(:,:) + 20. ! flag for fraction of field capacity
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
        call histrd1(iarchi,ier,vname,ik,wb(:,k),ifull)
        if ( wetfrac_found(k) ) then
          wb(:,k) = wb(:,k) + 20. ! flag for fraction of field capacity
        end if
      else if ( fnresid==1 ) then
        call histrd1(iarchi,ier,vname,ik,ucc,6*ik*ik,nogather=.false.)
        if ( wetfrac_found(k) ) then
          ucc(:) = ucc(:) + 20.   ! flag for fraction of field capacity
        end if
        call fill_cc1_gather(ucc,sea_a)
        call doints1_gather(ucc,wb(:,k))
      else
        call histrd1(iarchi,ier,vname,ik,ucc,6*ik*ik,nogather=.true.)
        if ( wetfrac_found(k) ) then
          ucc(:) = ucc(:) + 20.   ! flag for fraction of field capacity
        end if
        call fill_cc1_nogather(ucc,sea_a)
        call doints1_nogather(ucc,wb(:,k))
      end if ! iotest
    end do
  end if
  !unpack field capacity into volumetric soil moisture
  if ( any(wb(:,:)>10.) ) then
    if ( mydiag ) write(6,*) "Unpacking wetfrac to wb",wb(idjd,1)
    wb(:,:) = wb(:,:) - 20.
    do iq = 1,ifull
      isoil = isoilm(iq)
      wb(iq,:) = (1.-wb(iq,:))*swilt(isoil) + wb(iq,:)*sfc(isoil)
    end do
    if ( mydiag ) write(6,*) "giving wb",wb(idjd,1)
  end if
  call fillhist1('wetfac',wetfac,sea_a)
  where ( .not.land(1:ifull) )
    wetfac(:) = 1.
  end where

  !------------------------------------------------------------------
  ! Read 10m wind speeds for special sea roughness length calculations
  if ( nested==0 ) then
    if ( u10_found ) then
      call gethist1('u10',u10)
    else
      u10 = sqrt(u(1:ifull,1)**2+v(1:ifull,1)**2)*log(10./0.001)/log(zmin/0.001)
    end if
  end if

  !------------------------------------------------------------------
  ! Average fields
  if ( nested==0 ) then
    call gethist1('tscr_ave',tscr_ave)
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
  if ( nested==0 ) then
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
  ! Read CABLE/CASA aggregate C+N+P pools
  if ( nsib>=6 ) then
    if ( ccycle==0 ) then
      !if ( carbon_found ) then
      !  do k=1,ncp
      !    write(vname,'("cplant",I1.1)') k
      !    call fillhist1(vname,cplant(:,k),sea_a)
      !  end do
      !  do k=1,ncs
      !    write(vname,'("csoil",I1.1)') k
      !    call fillhist1(vname,csoil(:,k),sea_a)
      !  end do
      !end if
    else
      if ( carbon_found ) then
        call fillhist4('cplant',cplant,mplant,sea_a)
        call fillhist4('nplant',niplant,mplant,sea_a)
        call fillhist4('pplant',pplant,mplant,sea_a)
        call fillhist4('clitter',clitter,mlitter,sea_a)
        call fillhist4('nlitter',nilitter,mlitter,sea_a)
        call fillhist4('plitter',plitter,mlitter,sea_a)
        call fillhist4('csoil',csoil,msoil,sea_a)
        call fillhist4('nsoil',nisoil,msoil,sea_a)
        call fillhist4('psoil',psoil,msoil,sea_a)
        call fillhist1('glai',glai,sea_a)
      end if ! carbon_found
    end if   ! ccycle==0 ..else..
  end if     ! if nsib==6.or.nsib==7

  !------------------------------------------------------------------
  ! Read urban data
  if ( nurban/=0 ) then
    if ( .not.allocated(atebdwn) ) allocate(atebdwn(ifull,32))
    call fillhist4('rooftgg',atebdwn(:,1:5),  5,sea_a,filllimit=399.)
    call fillhist4('waletgg',atebdwn(:,6:10), 5,sea_a,filllimit=399.)
    call fillhist4('walwtgg',atebdwn(:,11:15),5,sea_a,filllimit=399.)
    call fillhist4('roadtgg',atebdwn(:,16:20),5,sea_a,filllimit=399.)
    call fillhist1('urbnsmc',atebdwn(:,21),sea_a,filllimit=399.)
    call fillhist1('urbnsmr',atebdwn(:,22),sea_a,filllimit=399.)
    call fillhist1('roofwtr',atebdwn(:,23),sea_a,filllimit=399.)
    call fillhist1('roadwtr',atebdwn(:,24),sea_a,filllimit=399.)
    call fillhist1('urbwtrc',atebdwn(:,25),sea_a,filllimit=399.)
    call fillhist1('urbwtrr',atebdwn(:,26),sea_a,filllimit=399.)
    call fillhist1('roofsnd',atebdwn(:,27),sea_a,filllimit=399.)
    call fillhist1('roadsnd',atebdwn(:,28),sea_a,filllimit=399.)
    call fillhist1('roofden',atebdwn(:,29),sea_a,filllimit=399.)
    if ( all(atebdwn(:,29)==0.) ) atebdwn(:,29)=100.
    call fillhist1('roadden',atebdwn(:,30),sea_a,filllimit=399.)
    if ( all(atebdwn(:,30)==0.) ) atebdwn(:,30)=100.
    call fillhist1('roofsna',atebdwn(:,31),sea_a,filllimit=399.)
    if ( all(atebdwn(:,31)==0.) ) atebdwn(:,31)=0.85
    call fillhist1('roadsna',atebdwn(:,32),sea_a,filllimit=399.)
    if ( all(atebdwn(:,32)==0.) ) atebdwn(:,32)=0.85
    do k = 1,20
      where ( atebdwn(:,k)>100. )
        atebdwn(:,k) = atebdwn(:,k) - urbtemp
      end where
    end do
  end if

  ! -----------------------------------------------------------------
  ! Read cloud fields
  if ( nested==0 .and. ldr/=0 ) then
    call gethist4a('qfg',qfg,5)               ! CLOUD FROZEN WATER
    call gethist4a('qlg',qlg,5)               ! CLOUD LIQUID WATER
    if ( ncloud>=2 ) then
      call gethist4a('qrg',qrg,5)             ! RAIN
    end if
    if ( ncloud>=3 ) then
      call gethist4a('qsng',qsng,5)           ! SNOW
      call gethist4a('qgrg',qgrg,5)           ! GRAUPEL
    end if
    call gethist4a('cfrac',cfrac,5)           ! CLOUD FRACTION
    if ( ncloud>=2 ) then
      call gethist4a('rfrac',rfrac,5)         ! RAIN FRACTION
    end if
    if ( ncloud>=3 ) then
      call gethist4a('sfrac',sfrac,5)         ! SNOW FRACTION
      call gethist4a('gfrac',gfrac,5)         ! GRAUPEL FRACTION
    end if
    if ( ncloud>=4 ) then
      call gethist4a('stratcf',stratcloud,5)  ! STRAT CLOUD FRACTION
      call gethist4a('strat_nt',nettend,5)    ! STRAT NET TENDENCY
    end if ! (ncloud>=4)
  end if   ! (nested==0)

  !------------------------------------------------------------------
  ! TKE-eps data
  if ( nested==0 .and. nvmix==6 ) then
    call gethist4a('tke',tke,5)
    if ( all(tke(1:ifull,:)==0.) ) tke(1:ifull,:)=1.5E-4
    call gethist4a('eps',eps,5)
    if  (all(eps(1:ifull,:)==0.) ) eps(1:ifull,:)=1.E-7
  end if

  !------------------------------------------------------------------
  ! Tracer data
  if ( nested==0 .and. ngas>0 ) then              
    do igas = 1,ngas              
      write(trnum,'(i3.3)') igas
      call gethist4a('tr'//trnum,tr(:,:,igas),7)
    end do
  end if

  !------------------------------------------------------------------
  ! Aerosol data ( non-nudged or diagnostic )
  if ( abs(iaero)>=2 ) then
    call gethist4a('seasalt1',ssn(:,:,1),5)
    call gethist4a('seasalt2',ssn(:,:,2),5)
    ! Factor 1.e3 to convert to g/m2, x 3 to get sulfate from sulfur
    so4t(:) = 0.
    do k = 1,kl
      so4t(:) = so4t(:) + 3.e3*xtgdwn(:,k,3)*(-1.e5*exp(psl(1:ifull))*dsig(k))/grav
    enddo
  end if
  
  ! -----------------------------------------------------------------
  ! Restart fields
  if ( nested==0 ) then
    ! DPSLDT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    dpsldt(:,:) = -999.
    if ( lrestart ) then
      call histrd4(iarchi,ier,'dpsldt',ik,kk,dpsldt,ifull)
    end if

    ! ZGNHS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    phi_nh(:,:) = 0.
    if ( lrestart ) then
      call histrd4(iarchi,ier,'zgnhs',ik,kk,phi_nh,ifull)
    end if
    
    ! SDOT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    sdot(:,:) = -999.
    if ( lrestart ) then
      sdot(:,1) = 0.
      call histrd4(iarchi,ier,'sdot',ik,kk,sdot(:,2:kk+1),ifull)
    end if

    ! PSLX !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    pslx(:,:) = -999.
    if ( lrestart ) then
      call histrd4(iarchi,ier,'pslx',ik,kk,pslx,ifull)
    end if
          
    ! SAVU !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    savu(:,:) = -999.
    if ( lrestart ) then
      call histrd4(iarchi,ier,'savu',ik,kk,savu,ifull)
    end if
          
    ! SAVV !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    savv(:,:) = -999.
    if ( lrestart ) then
      call histrd4(iarchi,ier,'savv',ik,kk,savv,ifull)
    end if

    ! SAVU1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    savu1(:,:) = -999.
    if ( lrestart ) then
      call histrd4(iarchi,ier,'savu1',ik,kk,savu1,ifull)
    end if
          
    ! SAVV1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    savv1(:,:) = -999.
    if ( lrestart ) then
      call histrd4(iarchi,ier,'savv1',ik,kk,savv1,ifull)
    end if

    ! SAVU2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    savu2(:,:) = -999.
    if ( lrestart ) then
      call histrd4(iarchi,ier,'savu2',ik,kk,savu2,ifull)
    end if
          
    ! SAVV2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    savv2(:,:) = -999.
    if ( lrestart ) then
      call histrd4(iarchi,ier,'savv2',ik,kk,savv2,ifull)
    end if

    ! OCEAN DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if ( abs(nmlo)>=3 .and. abs(nmlo)<=9 ) then
      oldu1(:,:) = 0.
      oldu2(:,:) = 0.
      oldv1(:,:) = 0.
      oldv2(:,:) = 0.
      ipice(:) = 0.
      if ( lrestart ) then
        call histrd4(iarchi,ier,'oldu1',ik,ok,oldu1,ifull)
        call histrd4(iarchi,ier,'oldv1',ik,ok,oldv1,ifull)
        call histrd4(iarchi,ier,'oldu2',ik,ok,oldu2,ifull)
        call histrd4(iarchi,ier,'oldv2',ik,ok,oldv2,ifull)
        call histrd1(iarchi,ier,'ipice',ik,ipice,ifull)
      end if
    end if
       
  end if ! (nested==0)

  ! -----------------------------------------------------------------
  ! soil ice and snow data
  call gethist4('wbice',wbice,ms) ! SOIL ICE
  call gethist4('tggsn',tggsn,3)
  if ( all(tggsn(:,:)==0.) ) tggsn(:,:) = 270.
  call gethist4('smass',smass,3)
  call gethist4('ssdn',ssdn,3)
  do k = 1,3
    if ( all(ssdn(:,k)==0.) ) then
      where ( snowd(:)>100. )
        ssdn(:,k)=240.
      elsewhere
        ssdn(:,k)=140.
      end where
    end if
  end do
  ssdnn(:) = ssdn(:,1)
  call gethist1('snage',snage)
  call gethist1('sflag',dum6)
  isflag(:) = nint(dum6(:))

  ! -----------------------------------------------------------------
  ! Misc fields
  ! sgsave is needed for convection
  if ( nested==0 ) then
    call gethist1('sgsave',sgsave)
  end if
        
endif    ! (nested/=1)

!**************************************************************
! This is the end of reading the initial arrays
!**************************************************************         

! -------------------------------------------------------------------
! tgg holds file surface temperature when no MLO
if ( nmlo==0 .or. abs(nmlo)>9 ) then
  where ( .not.land(1:ifull) )
    tgg(:,1) = tss
  end where
end if

! -------------------------------------------------------------------
! set-up for next read of file
iarchi  = iarchi + 1
kdate_s = kdate_r
ktime_s = ktime_r + 1

if ( myid==0 .and. nested==0 ) then
  write(6,*) "Final lrestart ",lrestart
end if

return
end subroutine onthefly_work


! *****************************************************************************
! INTERPOLATION ROUTINES                         

! Main interface
! Note that sx is a global array for all processors

subroutine doints1_nogather(s,sout)
      
use cc_mpi                 ! CC MPI routines
use infile                 ! Input file routines
use newmpar_m              ! Grid parameters
use parm_m                 ! Model configuration

implicit none
      
integer mm
real, dimension(:), intent(in) :: s
real, dimension(:), intent(inout) :: sout
real, dimension(ifull,m_fly) :: wrk
real, dimension(pil*pjl*pnpan,size(filemap),fncount) :: abuf
real, dimension(ik+4,ik+4,0:npanels) :: sx

call START_LOG(otf_ints1_begin)

if ( .not.allocated(filemap) ) then
  write(6,*) "ERROR: Mapping for RMA file windows has not been defined"
  call ccmpi_abort(-1)
end if

! This version uses MPI RMA to distribute data
call ccmpi_filewinget(abuf,s)

sx(1:ik+4,1:ik+4,0:npanels) = 0.
call ccmpi_filewinunpack(sx,abuf)
call sxpanelbounds(sx)

if ( nord==1 ) then   ! bilinear
  do mm = 1,m_fly     !  was 4, now may be 1
    call ints_blb(sx,wrk(:,mm),nface4(:,mm),xg4(:,mm),yg4(:,mm))
  end do
else                  ! bicubic
  do mm = 1,m_fly     !  was 4, now may be 1
    call intsb(sx,wrk(:,mm),nface4(:,mm),xg4(:,mm),yg4(:,mm))
  end do
end if   ! (nord==1)  .. else ..
sout(1:ifull) = sum(wrk(:,:), dim=2)/real(m_fly)

call END_LOG(otf_ints1_end)

return
end subroutine doints1_nogather

subroutine doints1_gather(s,sout)
      
use cc_mpi                 ! CC MPI routines
use infile                 ! Input file routines
use newmpar_m              ! Grid parameters
use parm_m                 ! Model configuration

implicit none
      
integer mm, n, ik2
real, dimension(:), intent(in) :: s
real, dimension(:), intent(inout) :: sout
real, dimension(ifull,m_fly) :: wrk
real, dimension(ik+4,ik+4,0:npanels) :: sx

call START_LOG(otf_ints1_begin)

if ( .not.bcst_allocated ) then
  write(6,*) "ERROR: Bcst commuicators have not been defined"
  call ccmpi_abort(-1)
end if

! This version uses MPI_Bcast to distribute data
sx(1:ik+4,1:ik+4,0:npanels) = 0.
if ( dk>0 ) then
  ik2 = ik*ik
  sx(3:ik+2,3:ik+2,0:npanels) = reshape( s(1:(npanels+1)*ik2), (/ ik, ik, npanels+1 /) )
  call sxpanelbounds(sx)
end if
do n = 0,npanels
  ! send each face of the host dataset to processes that require it
  if ( nfacereq(n) ) then
    call ccmpi_bcast(sx(:,:,n),0,comm_face(n))
  end if
end do  ! n loop

if ( nord==1 ) then   ! bilinear
  do mm = 1,m_fly     !  was 4, now may be 1
    call ints_blb(sx,wrk(:,mm),nface4(:,mm),xg4(:,mm),yg4(:,mm))
  end do
else                  ! bicubic
  do mm = 1,m_fly     !  was 4, now may be 1
    call intsb(sx,wrk(:,mm),nface4(:,mm),xg4(:,mm),yg4(:,mm))
  end do
end if   ! (nord==1)  .. else ..
sout(1:ifull) = sum(wrk(:,:), dim=2)/real(m_fly)

call END_LOG(otf_ints1_end)

return
end subroutine doints1_gather

subroutine doints4_nogather(s,sout)
      
use cc_mpi                 ! CC MPI routines
use infile                 ! Input file routines
use newmpar_m              ! Grid parameters
use parm_m                 ! Model configuration

implicit none
      
integer mm, k, kx, kb, ke, kn
real, dimension(:,:), intent(in) :: s
real, dimension(:,:), intent(inout) :: sout
real, dimension(ifull,m_fly) :: wrk
real, dimension(pil*pjl*pnpan,size(filemap),fncount,kblock) :: abuf
real, dimension(ik+4,ik+4,0:npanels) :: sx

call START_LOG(otf_ints4_begin)

kx = size(sout, 2)

if ( .not.allocated(filemap) ) then
  write(6,*) "ERROR: Mapping for RMA file windows has not been defined"
  call ccmpi_abort(-1)
end if

do kb = 1,kx,kblock
  ke = min(kb+kblock-1, kx)
  kn = ke - kb + 1

  ! This version uses MPI RMA to distribute data
  call ccmpi_filewinget(abuf(:,:,:,1:kn),s(:,kb:ke))
    
  ! MJT notes - sx can be made into a shared memory array,
  ! although this requires a MPI_Fence when the abuf
  ! arrays are unpacked for each level.
  
  if ( nord==1 ) then   ! bilinear
    do k = 1,kn
      sx(1:ik+4,1:ik+4,0:npanels) = 0.
      call ccmpi_filewinunpack(sx,abuf(:,:,:,k))
      call sxpanelbounds(sx)
      do mm = 1,m_fly     !  was 4, now may be 1
        call ints_blb(sx,wrk(:,mm),nface4(:,mm),xg4(:,mm),yg4(:,mm))
      end do
      sout(1:ifull,k+kb-1) = sum(wrk(:,:), dim=2)/real(m_fly)
    end do
  else                  ! bicubic
    do k = 1,kn
      sx(1:ik+4,1:ik+4,0:npanels) = 0.
      call ccmpi_filewinunpack(sx,abuf(:,:,:,k))
      call sxpanelbounds(sx)
      do mm = 1,m_fly     !  was 4, now may be 1
        call intsb(sx,wrk(:,mm),nface4(:,mm),xg4(:,mm),yg4(:,mm))
      end do
      sout(1:ifull,k+kb-1) = sum(wrk(:,:), dim=2)/real(m_fly)
    end do
  end if   ! (nord==1)  .. else ..

end do
  
call END_LOG(otf_ints4_end)

return
end subroutine doints4_nogather

subroutine doints4_gather(s,sout)
      
use cc_mpi                 ! CC MPI routines
use infile                 ! Input file routines
use newmpar_m              ! Grid parameters
use parm_m                 ! Model configuration

implicit none
      
integer mm, n, k, kx, ik2
real, dimension(:,:), intent(in) :: s
real, dimension(:,:), intent(inout) :: sout
real, dimension(ifull,m_fly) :: wrk
real, dimension(ik+4,ik+4,size(sout,2),0:npanels) :: sx
real, dimension(ik+4,ik+4,0:npanels) :: sy

call START_LOG(otf_ints4_begin)

kx = size(sout,2)

if ( .not.bcst_allocated ) then
  write(6,*) "ERROR: Bcst commuicators have not been defined"
  call ccmpi_abort(-1)
end if

! This version uses MPI_Bcast to distribute data
sx(1:ik+4,1:ik+4,1:kx,0:npanels) = 0.
if ( dk>0 ) then
  ik2 = ik*ik
  !     first extend s arrays into sx - this one -1:il+2 & -1:il+2
  do k = 1,kx
    sy(3:ik+2,3:ik+2,0:npanels) = reshape( s(1:(npanels+1)*ik2,k), (/ ik, ik, npanels+1 /) )
    call sxpanelbounds(sy(:,:,:))
    sx(:,:,k,:) = sy(:,:,:)
  end do
end if
do n = 0,npanels
  if ( nfacereq(n) ) then
    call ccmpi_bcast(sx(:,:,:,n),0,comm_face(n))
   end if
end do

do k = 1,kx
  sy(:,:,:) = sx(:,:,k,:)
  if ( nord==1 ) then   ! bilinear
    do mm = 1,m_fly     !  was 4, now may be 1
      call ints_blb(sy(:,:,:),wrk(:,mm),nface4(:,mm),xg4(:,mm),yg4(:,mm))
    end do
  else                  ! bicubic
    do mm = 1,m_fly     !  was 4, now may be 1
      call intsb(sy(:,:,:),wrk(:,mm),nface4(:,mm),xg4(:,mm),yg4(:,mm))
    end do
  end if   ! (nord==1)  .. else ..
  sout(1:ifull,k) = sum( wrk(:,:), dim=2 )/real(m_fly)
end do
  
call END_LOG(otf_ints4_end)

return
end subroutine doints4_gather

subroutine sxpanelbounds(sx_l)

use newmpar_m

implicit none

integer i, n, n_w, n_e, n_n, n_s
real, dimension(-1:ik+2,-1:ik+2,0:npanels), intent(inout) :: sx_l

do n = 0,npanels
  if ( nfacereq(n) ) then
    if ( mod(n,2)==0 ) then
      n_w = mod(n+5, 6)
      n_e = mod(n+2, 6)
      n_n = mod(n+1, 6)
      n_s = mod(n+4, 6)
      do i = 1,ik
        sx_l(0,i,n)    = sx_l(ik,i,n_w)
        sx_l(-1,i,n)   = sx_l(ik-1,i,n_w)
        sx_l(ik+1,i,n) = sx_l(ik+1-i,1,n_e)
        sx_l(ik+2,i,n) = sx_l(ik+1-i,2,n_e)
        sx_l(i,ik+1,n) = sx_l(i,1,n_n)
        sx_l(i,ik+2,n) = sx_l(i,2,n_n)
        sx_l(i,0,n)    = sx_l(ik,ik+1-i,n_s)
        sx_l(i,-1,n)   = sx_l(ik-1,ik+1-i,n_s)
      end do ! i
      sx_l(-1,0,n)      = sx_l(ik,2,n_w)        ! wws
      sx_l(0,-1,n)      = sx_l(ik,ik-1,n_s)     ! wss
      sx_l(0,0,n)       = sx_l(ik,1,n_w)        ! ws
      sx_l(ik+1,0,n)    = sx_l(ik,1,n_e)        ! es  
      sx_l(ik+2,0,n)    = sx_l(ik-1,1,n_e)      ! ees 
      sx_l(-1,ik+1,n)   = sx_l(ik,ik-1,n_w)     ! wwn
      sx_l(0,ik+2,n)    = sx_l(ik-1,ik,n_w)     ! wnn
      sx_l(ik+2,ik+1,n) = sx_l(2,1,n_e)         ! een  
      sx_l(ik+1,ik+2,n) = sx_l(1,2,n_e)         ! enn  
      sx_l(0,ik+1,n)    = sx_l(ik,ik,n_w)       ! wn  
      sx_l(ik+1,ik+1,n) = sx_l(1,1,n_e)         ! en  
      sx_l(ik+1,-1,n)   = sx_l(ik,2,n_e)        ! ess        
    else
      n_w = mod(n+4, 6)
      n_e = mod(n+1, 6)
      n_n = mod(n+2, 6)
      n_s = mod(n+5, 6)
      do i = 1,ik
        sx_l(0,i,n)    = sx_l(ik+1-i,ik,n_w)
        sx_l(-1,i,n)   = sx_l(ik+1-i,ik-1,n_w)
        sx_l(ik+1,i,n) = sx_l(1,i,n_e)
        sx_l(ik+2,i,n) = sx_l(2,i,n_e)
        sx_l(i,ik+1,n) = sx_l(1,ik+1-i,n_n)
        sx_l(i,ik+2,n) = sx_l(2,ik+1-i,n_n)
        sx_l(i,0,n)    = sx_l(i,ik,n_s)
        sx_l(i,-1,n)   = sx_l(i,ik-1,n_s)
      end do ! i
      sx_l(-1,0,n)      = sx_l(ik-1,ik,n_w)    ! wws
      sx_l(0,-1,n)      = sx_l(2,ik,n_s)       ! wss
      sx_l(0,0,n)       = sx_l(ik,ik,n_w)      ! ws
      sx_l(ik+1,0,n)    = sx_l(1,1,n_e)        ! es
      sx_l(ik+2,0,n)    = sx_l(1,2,n_e)        ! ees
      sx_l(-1,ik+1,n)   = sx_l(2,ik,n_w)       ! wwn   
      sx_l(0,ik+2,n)    = sx_l(1,ik-1,n_w)     ! wnn  
      sx_l(ik+2,ik+1,n) = sx_l(1,ik-1,n_e)     ! een  
      sx_l(ik+1,ik+2,n) = sx_l(2,ik,n_e)       ! enn  
      sx_l(0,ik+1,n)    = sx_l(1,ik,n_w)       ! wn  
      sx_l(ik+1,ik+1,n) = sx_l(1,ik,n_e)       ! en  
      sx_l(ik+1,-1,n)   = sx_l(2,1,n_e)        ! ess         
    end if   ! mod(n,2)==0 ..else..
  end if     ! nfacereq(n)
end do       ! n loop

return
end subroutine sxpanelbounds

subroutine intsb(sx_l,sout,nface_l,xg_l,yg_l)
      
!     same as subr ints, but with sout passed back and no B-S      
!     s is input; sout is output array
!     later may wish to save idel etc between array calls
!     this one does linear interp in x on outer y sides
!     doing x-interpolation before y-interpolation
!     This is a global routine 

use newmpar_m              ! Grid parameters
use parm_m                 ! Model configuration

implicit none

integer, dimension(ifull), intent(in) :: nface_l
integer :: idel, jdel
integer :: n, iq
real, dimension(ifull), intent(inout) :: sout
real xxg, yyg
real cmin, cmax
real, intent(in), dimension(ifull) :: xg_l, yg_l
real, dimension(-1:ik+2,-1:ik+2,0:npanels), intent(in) :: sx_l
real, dimension(2:3) :: dmul
real, dimension(1:4) :: cmul, emul, rmul

do iq = 1,ifull   ! runs through list of target points
  n = nface_l(iq)
  idel = int(xg_l(iq))
  xxg = xg_l(iq) - real(idel)
  jdel = int(yg_l(iq))
  yyg = yg_l(iq) - real(jdel)

  ! bi-cubic
  cmul(1) = (1.-xxg)*(2.-xxg)*(-xxg)/6.
  cmul(2) = (1.-xxg)*(2.-xxg)*(1.+xxg)/2.
  cmul(3) = xxg*(1.+xxg)*(2.-xxg)/2.
  cmul(4) = (1.-xxg)*(-xxg)*(1.+xxg)/6.
  dmul(2) = (1.-xxg)
  dmul(3) = xxg
  emul(1) = (1.-yyg)*(2.-yyg)*(-yyg)/6.
  emul(2) = (1.-yyg)*(2.-yyg)*(1.+yyg)/2.
  emul(3) = yyg*(1.+yyg)*(2.-yyg)/2.
  emul(4) = (1.-yyg)*(-yyg)*(1.+yyg)/6.
  cmin = minval( sx_l(idel:idel+1,jdel:jdel+1,n) )
  cmax = maxval( sx_l(idel:idel+1,jdel:jdel+1,n) )  
  rmul(1) = sum( sx_l(idel:idel+1,  jdel-1,n)*dmul(2:3) )
  rmul(2) = sum( sx_l(idel-1:idel+2,jdel,  n)*cmul(1:4) )
  rmul(3) = sum( sx_l(idel-1:idel+2,jdel+1,n)*cmul(1:4) )
  rmul(4) = sum( sx_l(idel:idel+1,  jdel+2,n)*dmul(2:3) )
  
  sout(iq) = min( max( cmin, sum( rmul(1:4)*emul(1:4) ) ), cmax ) ! Bermejo & Staniforth
end do    ! iq loop

return
end subroutine intsb

subroutine ints_blb(sx_l,sout,nface_l,xg_l,yg_l) 
      
!     this one does bi-linear interpolation only

use newmpar_m              ! Grid parameters
use parm_m                 ! Model configuration

implicit none

integer :: n, iq, idel, jdel
integer, intent(in), dimension(ifull) :: nface_l
real, dimension(ifull), intent(inout) :: sout
real, intent(in), dimension(ifull) :: xg_l, yg_l
real, dimension(-1:ik+2,-1:ik+2,0:npanels), intent(in) :: sx_l
real :: xxg, yyg

do iq = 1,ifull  ! runs through list of target points
  n = nface_l(iq)
  idel = int(xg_l(iq))
  xxg = xg_l(iq) - real(idel)
  jdel = int(yg_l(iq))
  yyg = yg_l(iq) - real(jdel)
  sout(iq) = yyg*(xxg*sx_l(idel+1,jdel+1,n) + (1.-xxg)*sx_l(idel,jdel+1,n)) + &
          (1.-yyg)*(xxg*sx_l(idel+1,jdel,n) + (1.-xxg)*sx_l(idel,jdel,n))
enddo    ! iq loop

return
end subroutine ints_blb

! *****************************************************************************
! FILL ROUTINES

subroutine fill_cc1_nogather(a_io,land_a)
      
! routine fills in interior of an array which has undefined points
! this version is for multiple input files

use cc_mpi          ! CC MPI routines
use infile          ! Input file routines

implicit none

integer nrem, j, n
integer ncount, cc, ipf
integer, dimension(pil) :: neighc
real, parameter :: value=999.       ! missing value flag
real, dimension(:), intent(inout) :: a_io
real, dimension(0:pil+1,0:pjl+1,pnpan,mynproc) :: c_io
real, dimension(pil,4) :: c
logical, dimension(:), intent(in) :: land_a
logical, dimension(pil,4) :: maskc

! only perform fill on processors reading input files
if ( fwsize==0 ) return
  
where ( land_a(1:fwsize) )
  a_io(1:fwsize) = value
end where
ncount = count( abs(a_io(1:fwsize)-value)<1.E-6 )
call ccmpi_allreduce(ncount,nrem,'sum',comm_ip)
if ( nrem==6*ik*ik ) return
 
do while ( nrem>0 )
  c_io(1:pil,1:pjl,1:pnpan,1:mynproc) = reshape( a_io(1:fwsize), (/ pil, pjl, pnpan, mynproc /) )
  call ccmpi_filebounds(c_io,comm_ip)
  do ipf = 1,mynproc
    do n = 1,pnpan
      do j = 1,pjl
        c(1:pil,1) = c_io(1:pil,j+1,n,ipf)
        c(1:pil,2) = c_io(1:pil,j-1,n,ipf)
        c(1:pil,3) = c_io(2:pil+1,j,n,ipf)
        c(1:pil,4) = c_io(0:pil-1,j,n,ipf)
        maskc(1:pil,1:4) = c(1:pil,1:4)/=value
        neighc(1:pil) = count( maskc(1:pil,1:4), dim=2 )
        cc = (j-1)*pil + (n-1)*pil*pjl + (ipf-1)*pil*pjl*pnpan
        where ( neighc(1:pil)>0 .and. c_io(1:pil,j,n,ipf)==value )
          a_io(1+cc:pil+cc) = sum( c(1:pil,1:4), mask=maskc(1:pil,1:4), dim=2 )/real(neighc(1:pil))
        end where
      end do
    end do
  end do
  ! test for convergence
  ncount = count( abs(a_io(1:fwsize)-value)<1.E-6 )
  call ccmpi_allreduce(ncount,nrem,'sum',comm_ip)
end do
      
return
end subroutine fill_cc1_nogather

subroutine fill_cc1_gather(a_io,land_a)
      
! routine fills in interior of an array which has undefined points
! this version is for a single input file

use cc_mpi          ! CC MPI routines
use infile          ! Input file routines

implicit none

integer nrem, i, iq, j, n
integer iminb, imaxb, jminb, jmaxb
integer is, ie, js, je
integer, dimension(0:5) :: imin, imax, jmin, jmax
integer, dimension(dk) :: neighb
integer, parameter, dimension(0:5) :: npann=(/1,103,3,105,5,101/)
integer, parameter, dimension(0:5) :: npane=(/102,2,104,4,100,0/)
integer, parameter, dimension(0:5) :: npanw=(/5,105,1,101,3,103/)
integer, parameter, dimension(0:5) :: npans=(/104,0,100,2,102,4/)
real, parameter :: value=999.       ! missing value flag
real, dimension(:), intent(inout) :: a_io
real, dimension(6*dk*dk) :: b_io
real, dimension(0:dk+1) :: a
real, dimension(dk) :: b_north, b_south, b_east, b_west
real, dimension(dk,4) :: b
logical, dimension(:), intent(in) :: land_a
logical, dimension(dk,4) :: mask
logical lflag

! only perform fill on myid==0
if ( dk==0 ) return

where ( land_a(1:6*dk*dk) )
  a_io(1:6*dk*dk) = value
end where
if ( all(abs(a_io(1:6*dk*dk)-value)<1.E-6) ) return

imin(0:5) = 1
imax(0:5) = dk
jmin(0:5) = 1
jmax(0:5) = dk
          
nrem = 1    ! Just for first iteration
do while ( nrem>0 )
  nrem = 0
  b_io(1:6*dk*dk) = a_io(1:6*dk*dk)
  ! MJT restricted fill
  do n = 0,5
    
    iminb = dk
    imaxb = 1
    jminb = dk
    jmaxb = 1
    
    ! north
    if (npann(n)<100) then
      do i = 1,dk
        iq=i+npann(n)*dk*dk
        b_north(i) = b_io(iq)
      end do
    else
      do i = 1,dk
        iq=1+(dk-i)*dk+(npann(n)-100)*dk*dk
        b_north(i) = b_io(iq)
      end do
    end if
    ! south
    if (npans(n)<100) then
      do i = 1,dk
        iq=i+(dk-1)*dk+npans(n)*dk*dk
        b_south(i) = b_io(iq)
      end do
    else
      do i = 1,dk
        iq=dk+(dk-i)*dk+(npans(n)-100)*dk*dk
        b_south(i) = b_io(iq)
      end do
    end if
    ! east
    if (npane(n)<100) then
      do j = 1,dk
        iq=1+(j-1)*dk+npane(n)*dk*dk
        b_east(j) = b_io(iq)
      end do
    else
      do j = 1,dk
        iq=dk+1-j+(npane(n)-100)*dk*dk
        b_east(j) = b_io(iq)
      end do
    end if
    ! west
    if (npanw(n)<100) then
      do j = 1,dk
        iq=dk+(j-1)*dk+npanw(n)*dk*dk
        b_west(j) = b_io(iq)
      end do
    else
      do j = 1,dk
        iq=dk+1-j+(dk-1)*dk+(npanw(n)-100)*dk*dk
        b_west(j) = b_io(iq)
      end do
    end if

    is = imin(n)
    ie = imax(n)
    js = jmin(n)
    je = jmax(n)
    
    if ( js==1 ) then
      ! j = 1
      a(0)     = b_west(1)
      a(dk+1)  = b_east(1)
      a(max(is-1,1))  = b_io(max(is-1,1)+n*dk*dk)
      a(min(ie+1,dk)) = b_io(min(ie+1,dk)+n*dk*dk)
      a(is:ie) = b_io(is+n*dk*dk:ie+n*dk*dk)
      b(is:ie,1) = b_io(is+dk+n*dk*dk:ie+dk+n*dk*dk) ! north
      b(is:ie,2) = b_south(is:ie)                    ! south
      b(is:ie,3) = a(is+1:ie+1)                      ! east
      b(is:ie,4) = a(is-1:ie-1)                      ! west
      mask(is:ie,1:4) = b(is:ie,1:4)/=value
      neighb(is:ie) = count( mask(is:ie,1:4), dim=2)
      where ( neighb(is:ie)>0 .and. a(is:ie)==value )
        a_io(is+n*dk*dk:ie+n*dk*dk) = sum( b(is:ie,1:4), mask=mask(is:ie,1:4), dim=2)/real(neighb(is:ie))
      end where
      lflag = .false.
      do i = is,ie
        if ( neighb(i)==0 ) then
          nrem = nrem + 1 ! current number of points without a neighbour
          iminb = min(i, iminb)
          imaxb = max(i, imaxb)
          lflag = .true.
        end if
      end do
      if ( lflag ) then
        jminb = min(1, jminb)
        jmaxb = max(1, jmaxb)
      end if
    end if
    do j = max(js,2),min(je,dk-1)
      a(0)     = b_west(j)
      a(dk+1)  = b_east(j)
      a(max(is-1,1))  = b_io(max(is-1,1)+(j-1)*dk+n*dk*dk)
      a(min(ie+1,dk)) = b_io(min(ie+1,dk)+(j-1)*dk+n*dk*dk)
      a(is:ie) = b_io(is+(j-1)*dk+n*dk*dk:ie+(j-1)*dk+n*dk*dk)
      b(is:ie,1) = b_io(is+j*dk+n*dk*dk:ie+j*dk+n*dk*dk)         ! north
      b(is:ie,2) = b_io(is+(j-2)*dk+n*dk*dk:ie+(j-2)*dk+n*dk*dk) ! south
      b(is:ie,3) = a(is+1:ie+1)                                  ! east
      b(is:ie,4) = a(is-1:ie-1)                                  ! west
      mask(is:ie,1:4) = b(is:ie,1:4)/=value
      neighb(is:ie) = count( mask(is:ie,1:4), dim=2)
      where ( neighb(is:ie)>0 .and. a(is:ie)==value )
        a_io(is+(j-1)*dk+n*dk*dk:ie+(j-1)*dk+n*dk*dk) = sum( b(is:ie,1:4), mask=mask(is:ie,1:4), dim=2)/real(neighb(is:ie))
      end where
      lflag = .false.
      do i = is,ie
        if ( neighb(i)==0 ) then
          nrem = nrem + 1 ! current number of points without a neighbour
          iminb = min(i, iminb)
          imaxb = max(i, imaxb)
          lflag = .true.
        end if
      end do
      if ( lflag ) then
        jminb = min(j, jminb)
        jmaxb = max(j, jmaxb)
      end if
    end do
    if ( je==dk ) then
      ! j = dk
      a(0)     = b_west(dk)
      a(dk+1)  = b_east(dk)
      a(max(is-1,1))  = b_io(max(is-1,1)-dk+(n+1)*dk*dk)
      a(min(ie+1,dk)) = b_io(min(ie+1,dk)-dk+(n+1)*dk*dk)
      a(is:ie) = b_io(is-dk+(n+1)*dk*dk:ie-dk+(n+1)*dk*dk)
      b(is:ie,1) = b_north(is:ie)                                ! north
      b(is:ie,2) = b_io(is-2*dk+(n+1)*dk*dk:ie-2*dk+(n+1)*dk*dk) ! south
      b(is:ie,3) = a(is+1:ie+1)                                  ! east
      b(is:ie,4) = a(is-1:ie-1)                                  ! west
      mask(is:ie,1:4) = b(is:ie,1:4)/=value
      neighb(is:ie) = count( mask(is:ie,1:4), dim=2)
      where ( neighb(is:ie)>0 .and. a(is:ie)==value )
        a_io(is-dk+(n+1)*dk*dk:ie-dk+(n+1)*dk*dk) = sum( b(is:ie,1:4), mask=mask(is:ie,1:4), dim=2)/real(neighb(is:ie))
      end where
      lflag = .false.
      do i = is,ie
        if ( neighb(i)==0 ) then
          nrem = nrem + 1 ! current number of points without a neighbour
          iminb = min(i, iminb)
          imaxb = max(i, imaxb)
          lflag = .true.
        end if
      end do
      if ( lflag ) then
        jminb = min(dk, jminb)
        jmaxb = max(dk, jmaxb)
      end if
    end if
    
    imin(n) = iminb
    imax(n) = imaxb
    jmin(n) = jminb
    jmax(n) = jmaxb
  end do
end do
  
return
end subroutine fill_cc1_gather

subroutine fill_cc4_nogather(a_io,land_a)
      
! routine fills in interior of an array which has undefined points
! this version is distributed over processes with input files

use cc_mpi          ! CC MPI routines
use infile          ! Input file routines

implicit none

integer nrem, j, n, k, kx
integer ncount, cc, ipf
integer, dimension(pil) :: neighc
real, parameter :: value=999.       ! missing value flag
real, dimension(:,:), intent(inout) :: a_io
real, dimension(0:pil+1,0:pjl+1,pnpan,mynproc,size(a_io,2)) :: c_io
real, dimension(pil,4) :: c
logical, dimension(:), intent(in) :: land_a
logical, dimension(pil,4) :: maskc

kx = size(a_io,2)

! only perform fill on processors reading input files
if ( fwsize==0 ) return

do k = 1,kx
  where ( land_a(1:fwsize) )
    a_io(1:fwsize,k) = value
  end where
end do
ncount = count( abs(a_io(1:fwsize,kx)-value)<1.E-6 )
call ccmpi_allreduce(ncount,nrem,'sum',comm_ip)
if ( nrem==6*ik*ik ) return
 
do while ( nrem > 0 )
  c_io(1:pil,1:pjl,1:pnpan,1:mynproc,1:kx) = reshape( a_io(1:fwsize,1:kx), (/ pil, pjl, pnpan, mynproc, kx /) )
  call ccmpi_filebounds(c_io,comm_ip)
  do k = 1,kx
    do ipf = 1,mynproc
      do n = 1,pnpan
        do j = 1,pjl
          c(1:pil,1) = c_io(1:pil,j+1,n,ipf,k)
          c(1:pil,2) = c_io(1:pil,j-1,n,ipf,k)
          c(1:pil,3) = c_io(2:pil+1,j,n,ipf,k)
          c(1:pil,4) = c_io(0:pil-1,j,n,ipf,k)
          maskc(1:pil,1:4) = c(1:pil,1:4)/=value
          neighc(1:pil) = count( maskc(1:pil,1:4), dim=2)
          cc = (j-1)*pil + (n-1)*pil*pjl + (ipf-1)*pil*pjl*pnpan
          where ( neighc(1:pil)>0 .and. c_io(1:pil,j,n,ipf,k)==value )
            a_io(1+cc:pil+cc,k) = sum( c(1:pil,1:4), mask=maskc(1:pil,1:4), dim=2)/real(neighc(1:pil))
          end where
        end do
      end do
    end do
  end do
  ! test for convergence
  ncount = count( abs(a_io(1:fwsize,kx)-value)<1.E-6 )
  call ccmpi_allreduce(ncount,nrem,'sum',comm_ip)
end do
      
return
end subroutine fill_cc4_nogather

subroutine fill_cc4_gather(a_io,land_a)
      
! routine fills in interior of an array which has undefined points
! this version is for a single input file

use cc_mpi          ! CC MPI routines
use infile          ! Input file routines

implicit none

integer nrem, i, iq, j, n, k, kx
integer iminb, imaxb, jminb, jmaxb
integer is, ie, js, je
integer, dimension(0:5) :: imin, imax, jmin, jmax
integer, dimension(dk) :: neighb
integer, parameter, dimension(0:5) :: npann=(/1,103,3,105,5,101/)
integer, parameter, dimension(0:5) :: npane=(/102,2,104,4,100,0/)
integer, parameter, dimension(0:5) :: npanw=(/5,105,1,101,3,103/)
integer, parameter, dimension(0:5) :: npans=(/104,0,100,2,102,4/)
real, parameter :: value=999.       ! missing value flag
real, dimension(:,:), intent(inout) :: a_io
real, dimension(6*dk*dk,size(a_io,2)) :: b_io
real, dimension(0:dk+1) :: a
real, dimension(dk,size(a_io,2)) :: b_north, b_south, b_east, b_west
real, dimension(dk,4) :: b
logical, dimension(:), intent(in) :: land_a
logical, dimension(dk,4) :: mask
logical lflag

kx = size(a_io,2)

! only perform fill on myid==0
if ( dk==0 ) return

do k = 1,kx
  where ( land_a(1:6*dk*dk) )
    a_io(1:6*dk*dk,k)=value
  end where
end do
if ( all(abs(a_io(1:6*dk*dk,kx)-value)<1.E-6) ) return

imin(0:5) = 1
imax(0:5) = dk
jmin(0:5) = 1
jmax(0:5) = dk
          
nrem = 1    ! Just for first iteration
do while ( nrem>0 )
  nrem = 0
  b_io(1:6*dk*dk,1:kx) = a_io(1:6*dk*dk,1:kx)
  ! MJT restricted fill
  do n = 0,5
       
    iminb = dk
    imaxb = 1
    jminb = dk
    jmaxb = 1
    
    ! north
    if (npann(n)<100) then
      do k = 1,kx
        do i = 1,dk
          iq=i+npann(n)*dk*dk
          b_north(i,k) = b_io(iq,k)
        end do
      end do
    else
      do k = 1,kx
        do i = 1,dk
          iq=1+(dk-i)*dk+(npann(n)-100)*dk*dk
          b_north(i,k) = b_io(iq,k)
        end do
      end do
    end if
    ! south
    if (npans(n)<100) then
      do k = 1,kx
        do i = 1,dk
          iq=i+(dk-1)*dk+npans(n)*dk*dk
          b_south(i,k) = b_io(iq,k)
        end do
      end do
    else
      do k = 1,kx
        do i = 1,dk
          iq=dk+(dk-i)*dk+(npans(n)-100)*dk*dk
          b_south(i,k) = b_io(iq,k)
        end do
      end do
    end if
    ! east
    if (npane(n)<100) then
      do k = 1,kx
        do j = 1,dk
          iq=1+(j-1)*dk+npane(n)*dk*dk
          b_east(j,k) = b_io(iq,k)
        end do
      end do
    else
      do k = 1,kx
        do j = 1,dk
          iq=dk+1-j+(npane(n)-100)*dk*dk
          b_east(j,k) = b_io(iq,k)
        end do
      end do
    end if
    ! west
    if (npanw(n)<100) then
      do k = 1,kx
        do j = 1,dk
          iq=dk+(j-1)*dk+npanw(n)*dk*dk
          b_west(j,k) = b_io(iq,k)
        end do
      end do
    else
      do k = 1,kx
        do j = 1,dk
          iq=dk+1-j+(dk-1)*dk+(npanw(n)-100)*dk*dk
          b_west(j,k) = b_io(iq,k)
        end do
      end do
    end if

    is = imin(n)
    ie = imax(n)
    js = jmin(n)
    je = jmax(n)
    
    if ( js==1 ) then
      ! j = 1
      do k = 1,kx
        a(0)     = b_west(1,k)
        a(dk+1)  = b_east(1,k)
        a(max(is-1,1))  = b_io(max(is-1,1)+n*dk*dk,k)
        a(min(ie+1,dk)) = b_io(min(ie+1,dk)+n*dk*dk,k)
        a(is:ie) = b_io(is+n*dk*dk:ie+n*dk*dk,k)
        b(is:ie,1) = b_io(is+dk+n*dk*dk:ie+dk+n*dk*dk,k) ! north
        b(is:ie,2) = b_south(is:ie,k)                    ! south
        b(is:ie,3) = a(is+1:ie+1)                        ! east
        b(is:ie,4) = a(is-1:ie-1)                        ! west
        mask(is:ie,1:4) = b(is:ie,1:4)/=value
        neighb(is:ie) = count( mask(is:ie,1:4), dim=2 )
        where ( neighb(is:ie)>0 .and. a(is:ie)==value )
          a_io(is+n*dk*dk:ie+n*dk*dk,k) = sum( b(is:ie,1:4), mask=mask(is:ie,1:4), dim=2 )/real(neighb(is:ie))
        end where
      end do
      lflag = .false.
      do i = is,ie
        if ( neighb(i)==0 ) then
          nrem = nrem + 1 ! current number of points without a neighbour
          iminb = min(i, iminb)
          imaxb = max(i, imaxb)
          lflag = .true.
        end if
      end do
      if ( lflag ) then
        jminb = min(1, jminb)
        jmaxb = max(1, jmaxb)
      end if
    end if
    do j = max(js,2),min(je,dk-1)
      do k = 1,kx
        a(0)     = b_west(j,k)
        a(dk+1)  = b_east(j,k)
        a(max(is-1,1))  = b_io(max(is-1,1)+(j-1)*dk+n*dk*dk,k)
        a(min(ie+1,dk)) = b_io(min(ie+1,dk)+(j-1)*dk+n*dk*dk,k)
        a(is:ie) = b_io(is+(j-1)*dk+n*dk*dk:ie+(j-1)*dk+n*dk*dk,k)
        b(is:ie,1) = b_io(is+j*dk+n*dk*dk:ie+j*dk+n*dk*dk,k)         ! north
        b(is:ie,2) = b_io(is+(j-2)*dk+n*dk*dk:ie+(j-2)*dk+n*dk*dk,k) ! south
        b(is:ie,3) = a(is+1:ie+1)                                    ! east
        b(is:ie,4) = a(is-1:ie-1)                                    ! west
        mask(is:ie,1:4) = b(is:ie,1:4)/=value
        neighb(is:ie) = count( mask(is:ie,1:4), dim=2)
        where ( neighb(is:ie)>0 .and. a(is:ie)==value )
          a_io(is+(j-1)*dk+n*dk*dk:ie+(j-1)*dk+n*dk*dk,k) = sum( b(is:ie,1:4), mask=mask(is:ie,1:4), dim=2)/real(neighb(is:ie))
        end where
      end do
      lflag = .false.
      do i = is,ie
        if ( neighb(i)==0 ) then
          nrem = nrem + 1 ! current number of points without a neighbour
          iminb = min(i, iminb)
          imaxb = max(i, imaxb)
          lflag = .true.
        end if
      end do
      if ( lflag ) then
        jminb = min(j, jminb)
        jmaxb = max(j, jmaxb)
      end if
    end do
    if ( je==dk ) then
      ! j = dk
      do k = 1,kx
        a(0)     = b_west(dk,k)
        a(dk+1)  = b_east(dk,k)
        a(max(is-1,1))  = b_io(max(is-1,1)-dk+(n+1)*dk*dk,k)
        a(min(ie+1,dk)) = b_io(min(ie+1,dk)-dk+(n+1)*dk*dk,k)
        a(is:ie) = b_io(is-dk+(n+1)*dk*dk:ie-dk+(n+1)*dk*dk,k)
        b(is:ie,1) = b_north(is:ie,k)                                ! north
        b(is:ie,2) = b_io(is-2*dk+(n+1)*dk*dk:ie-2*dk+(n+1)*dk*dk,k) ! south
        b(is:ie,3) = a(is+1:ie+1)                                    ! east
        b(is:ie,4) = a(is-1:ie-1)                                    ! west
        mask(is:ie,1:4) = b(is:ie,1:4)/=value
        neighb(is:ie) = count( mask(is:ie,1:4), dim=2)
        where ( neighb(is:ie)>0 .and. a(is:ie)==value )
          a_io(is-dk+(n+1)*dk*dk:ie-dk+(n+1)*dk*dk,k) = sum( b(is:ie,1:4), mask=mask(is:ie,1:4), dim=2)/real(neighb(is:ie))
        end where
      end do
      lflag = .false.
      do i = is,ie
        if ( neighb(i)==0 ) then
          nrem = nrem + 1 ! current number of points without a neighbour
          iminb = min(i, iminb)
          imaxb = max(i, imaxb)
          lflag = .true.
        end if
      end do
      if ( lflag ) then
        jminb = min(dk, jminb)
        jmaxb = max(dk, jmaxb)
      end if
    end if
    
    imin(n) = iminb
    imax(n) = imaxb
    jmin(n) = jminb
    jmax(n) = jmaxb
  end do
end do
      
return
end subroutine fill_cc4_gather

! *****************************************************************************
! OROGRAPHIC ADJUSTMENT ROUTINES

subroutine mslpx(pmsl,psl,zs,t,siglev)

use newmpar_m              ! Grid parameters
use parm_m                 ! Model configuration
use sigs_m                 ! Atmosphere sigma levels
      
!     this one will ignore negative zs (i.e. over the ocean)

implicit none
      
include 'const_phys.h'    ! Physical constants

integer nfull
real siglev, c, con, conr
real, dimension(:), intent(inout) :: pmsl, psl, zs, t
real, dimension(size(pmsl)) :: dlnps, phi1, tav, tsurf

nfull = size(pmsl)
c     = grav/stdlapse
conr  = c/rdry
con   = siglev**(rdry/c)/c

phi1(1:nfull)  = t(1:nfull)*rdry*(1.-siglev)/siglev ! phi of sig(lev) above sfce
tsurf(1:nfull) = t(1:nfull)+phi1(1:nfull)*stdlapse/grav
tav(1:nfull)   = tsurf(1:nfull)+max(0.,zs(1:nfull))*.5*stdlapse/grav
dlnps(1:nfull) = max(0.,zs(1:nfull))/(rdry*tav(1:nfull))
pmsl(1:nfull)  = 1.e5*exp(psl(1:nfull)+dlnps(1:nfull))

return
end subroutine mslpx
      
subroutine to_pslx(pmsl,psl,zs,t,levk)

use newmpar_m              ! Grid parameters
use parm_m                 ! Model configuration
use sigs_m                 ! Atmosphere sigma levels
      
implicit none
      
include 'const_phys.h'     ! Physical constants
      
integer levk
real, dimension(ifull) :: pmsl, psl, zs, t
real, dimension(ifull) :: dlnps, phi1, tav, tsurf

phi1(:)  = t(:)*rdry*(1.-sig(levk))/sig(levk) ! phi of sig(levk) above sfce
tsurf(:) = t(:) + phi1(:)*stdlapse/grav
tav(:)   = tsurf(:) + max( 0., zs(:) )*.5*stdlapse/grav
dlnps(:) = max( 0., zs(:))/(rdry*tav(:) )
psl(:)   = log(1.e-5*pmsl(:)) - dlnps(:)

#ifdef debug
if ( nmaxpr==1 .and. mydiag ) then
  write(6,*)'to_psl levk,sig(levk) ',levk,sig(levk)
  write(6,*)'zs,t_lev,psl,pmsl ',zs(idjd),t(idjd),psl(idjd),pmsl(idjd)
end if
#endif

return
end subroutine to_pslx

subroutine retopo(psl,zsold,zs,t,qg)
!     in Jan 07, renamed recalc of ps from here, to reduce confusion     
!     this version (Aug 2003) allows -ve zsold (from spectral model),
!     but assumes new zs is positive for atmospheric purposes
!     this routine redefines psl, t to compensate for zsold going to zs
!     (but does not overwrite zs, ps themselves here)
!     called by indata and nestin for newtop>=1
!     nowadays just for ps and atmospheric fields Mon  08-23-1999
use cc_mpi, only : mydiag
use diag_m
use newmpar_m
use parm_m
use sigs_m

implicit none

include 'const_phys.h'

real, dimension(:), intent(inout) :: psl
real, dimension(:), intent(in) :: zsold, zs
real, dimension(:,:), intent(inout) :: t, qg
real, dimension(size(psl)) :: psnew, psold, pslold
real, dimension(kl) :: told, qgold
real sig2
integer iq, k, kkk, kold

pslold(1:ifull) = psl(1:ifull)
psold(1:ifull)  = 1.e5*exp(psl(1:ifull))
psl(1:ifull)    = psl(1:ifull) + (zsold(1:ifull)-zs(1:ifull))/(rdry*t(1:ifull,1))
psnew(1:ifull)  = 1.e5*exp(psl(1:ifull))

!     now alter temperatures to compensate for new topography
if ( ktau<100 .and. mydiag ) then
  write(6,*) 'retopo: zsold,zs,psold,psnew ',zsold(idjd),zs(idjd),psold(idjd),psnew(idjd)
  write(6,*) 'retopo: old t ',(t(idjd,k),k=1,kl)
  write(6,*) 'retopo: old qg ',(qg(idjd,k),k=1,kl)
end if  ! (ktau.lt.100)

do iq = 1,ifull
  qgold(1:kl) = qg(iq,1:kl)
  told(1:kl)  = t(iq,1:kl)
  kold = 2
  !do k = 1,kl-1
  do k = 1,kl ! MJT suggestion
    sig2 = sig(k)*psnew(iq)/psold(iq)
    if ( sig2>=sig(1) ) then
      ! assume 6.5 deg/km, with dsig=.1 corresponding to 1 km
      t(iq,k) = told(1) + (sig2-sig(1))*6.5/0.1  
    else
      do kkk = kold,kl-1
        if ( sig2>sig(kkk) ) exit
      end do
      kold = kkk
      t(iq,k)  = (told(kkk)*(sig(kkk-1)-sig2)+told(kkk-1)*(sig2-sig(kkk)))/(sig(kkk-1)-sig(kkk))
      qg(iq,k) = (qgold(kkk)*(sig(kkk-1)-sig2)+qgold(kkk-1)*(sig2-sig(kkk)))/(sig(kkk-1)-sig(kkk))
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

subroutine interpwind4(uct,vct,ucc,vcc,nogather)
      
use cc_mpi           ! CC MPI routines
use newmpar_m        ! Grid parameters
use vecsuv_m         ! Map to cartesian coordinates
      
implicit none
      
integer k
real, dimension(:,:), intent(inout) :: ucc, vcc
real, dimension(size(ucc,1),kk) :: wcc
real, dimension(size(ucc,1)) :: uc, vc, wc
real, dimension(ifull,kk), intent(out) :: uct, vct
real, dimension(ifull,kk) :: wct
real, dimension(ifull) :: newu, newv, neww
logical, intent(in), optional :: nogather
logical ngflag

ngflag = .false.
if ( present(nogather) ) then
  ngflag = nogather
end if

if ( ngflag ) then
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
  call doints4_nogather(ucc, uct)
  call doints4_nogather(vcc, vct)
  call doints4_nogather(wcc, wct)
else
  ! dk is only non-zero on myid==0
  if ( dk>0 ) then
    do k = 1,kk
      ! first set up winds in Cartesian "source" coords            
      uc(1:6*dk*dk) = axs_a(1:6*dk*dk)*ucc(1:6*dk*dk,k) + bxs_a(1:6*dk*dk)*vcc(1:6*dk*dk,k)
      vc(1:6*dk*dk) = ays_a(1:6*dk*dk)*ucc(1:6*dk*dk,k) + bys_a(1:6*dk*dk)*vcc(1:6*dk*dk,k)
      wc(1:6*dk*dk) = azs_a(1:6*dk*dk)*ucc(1:6*dk*dk,k) + bzs_a(1:6*dk*dk)*vcc(1:6*dk*dk,k)
      ! now convert to winds in "absolute" Cartesian components
      ucc(1:6*dk*dk,k) = uc(1:6*dk*dk)*rotpoles(1,1) + vc(1:6*dk*dk)*rotpoles(1,2) + wc(1:6*dk*dk)*rotpoles(1,3)
      vcc(1:6*dk*dk,k) = uc(1:6*dk*dk)*rotpoles(2,1) + vc(1:6*dk*dk)*rotpoles(2,2) + wc(1:6*dk*dk)*rotpoles(2,3)
      wcc(1:6*dk*dk,k) = uc(1:6*dk*dk)*rotpoles(3,1) + vc(1:6*dk*dk)*rotpoles(3,2) + wc(1:6*dk*dk)*rotpoles(3,3)
    end do    ! k loop
  end if      ! dk>0
  ! interpolate all required arrays to new C-C positions
  ! do not need to do map factors and Coriolis on target grid
  call doints4_gather(ucc, uct)
  call doints4_gather(vcc, vct)
  call doints4_gather(wcc, wct)
end if
  
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

subroutine interpcurrent1(uct,vct,ucc,vcc,mask_a,nogather)
      
use cc_mpi           ! CC MPI routines
use newmpar_m        ! Grid parameters
use vecsuv_m         ! Map to cartesian coordinates
      
implicit none
      
real, dimension(:), intent(inout) :: ucc, vcc
real, dimension(size(ucc)) :: wcc
real, dimension(size(ucc)) :: uc, vc, wc
real, dimension(ifull), intent(out) :: uct, vct
real, dimension(ifull) :: wct
real, dimension(ifull) :: newu, newv, neww
logical, dimension(:), intent(in) :: mask_a
logical, intent(in), optional :: nogather
logical ngflag

ngflag = .false.
if ( present(nogather) ) then
  ngflag = nogather
end if

if ( ngflag ) then
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
    call fill_cc1_nogather(ucc, mask_a)
    call fill_cc1_nogather(vcc, mask_a)
    call fill_cc1_nogather(wcc, mask_a)
  end if
  call doints1_nogather(ucc, uct)
  call doints1_nogather(vcc, vct)
  call doints1_nogather(wcc, wct)
else
  ! dk is only non-zero on myid==0
  if ( dk>0 ) then
    ! first set up currents in Cartesian "source" coords            
    uc(1:6*dk*dk) = axs_a(1:6*dk*dk)*ucc(1:6*dk*dk) + bxs_a(1:6*dk*dk)*vcc(1:6*dk*dk)
    vc(1:6*dk*dk) = ays_a(1:6*dk*dk)*ucc(1:6*dk*dk) + bys_a(1:6*dk*dk)*vcc(1:6*dk*dk)
    wc(1:6*dk*dk) = azs_a(1:6*dk*dk)*ucc(1:6*dk*dk) + bzs_a(1:6*dk*dk)*vcc(1:6*dk*dk)
    ! now convert to winds in "absolute" Cartesian components
    ucc(1:6*dk*dk) = uc(1:6*dk*dk)*rotpoles(1,1) + vc(1:6*dk*dk)*rotpoles(1,2) + wc(1:6*dk*dk)*rotpoles(1,3)
    vcc(1:6*dk*dk) = uc(1:6*dk*dk)*rotpoles(2,1) + vc(1:6*dk*dk)*rotpoles(2,2) + wc(1:6*dk*dk)*rotpoles(2,3)
    wcc(1:6*dk*dk) = uc(1:6*dk*dk)*rotpoles(3,1) + vc(1:6*dk*dk)*rotpoles(3,2) + wc(1:6*dk*dk)*rotpoles(3,3)
    ! interpolate all required arrays to new C-C positions
    ! do not need to do map factors and Coriolis on target grid
    call fill_cc1_gather(ucc, mask_a)
    call fill_cc1_gather(vcc, mask_a)
    call fill_cc1_gather(wcc, mask_a)
  end if ! dk>0
  call doints1_gather(ucc, uct)
  call doints1_gather(vcc, vct)
  call doints1_gather(wcc, wct)
end if
  
! now convert to "target" Cartesian components (transpose used)
newu(1:ifull) = uct(1:ifull)*rotpole(1,1) + vct(1:ifull)*rotpole(2,1) + wct(1:ifull)*rotpole(3,1)
newv(1:ifull) = uct(1:ifull)*rotpole(1,2) + vct(1:ifull)*rotpole(2,2) + wct(1:ifull)*rotpole(3,2)
neww(1:ifull) = uct(1:ifull)*rotpole(1,3) + vct(1:ifull)*rotpole(2,3) + wct(1:ifull)*rotpole(3,3)
! then finally to "target" local x-y components
uct(1:ifull) = ax(1:ifull)*newu(1:ifull) + ay(1:ifull)*newv(1:ifull) + az(1:ifull)*neww(1:ifull)
vct(1:ifull) = bx(1:ifull)*newu(1:ifull) + by(1:ifull)*newv(1:ifull) + bz(1:ifull)*neww(1:ifull)

return
end subroutine interpcurrent1

subroutine interpcurrent4(uct,vct,ucc,vcc,mask_a,nogather)
      
use cc_mpi           ! CC MPI routines
use newmpar_m        ! Grid parameters
use vecsuv_m         ! Map to cartesian coordinates
      
implicit none
      
integer k
real, dimension(:,:), intent(inout) :: ucc, vcc
real, dimension(size(ucc,1),ok) :: wcc
real, dimension(ifull,ok), intent(out) :: uct, vct
real, dimension(ifull,ok) :: wct
real, dimension(size(ucc,1)) :: uc, vc, wc
real, dimension(ifull) :: newu, newv, neww
logical, dimension(:), intent(in) :: mask_a
logical, intent(in), optional :: nogather
logical ngflag

ngflag = .false.
if ( present(nogather) ) then
  ngflag = nogather
end if

if ( ngflag ) then
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
    call fill_cc4_nogather(ucc, mask_a)
    call fill_cc4_nogather(vcc, mask_a)
    call fill_cc4_nogather(wcc, mask_a)
  end if
  call doints4_nogather(ucc, uct)
  call doints4_nogather(vcc, vct)
  call doints4_nogather(wcc, wct)
else
  ! dk is only non-zero on myid==0
  if ( dk>0 ) then
    do k = 1,ok
      ! first set up currents in Cartesian "source" coords            
      uc(1:6*dk*dk) = axs_a(1:6*dk*dk)*ucc(1:6*dk*dk,k) + bxs_a(1:6*dk*dk)*vcc(1:6*dk*dk,k)
      vc(1:6*dk*dk) = ays_a(1:6*dk*dk)*ucc(1:6*dk*dk,k) + bys_a(1:6*dk*dk)*vcc(1:6*dk*dk,k)
      wc(1:6*dk*dk) = azs_a(1:6*dk*dk)*ucc(1:6*dk*dk,k) + bzs_a(1:6*dk*dk)*vcc(1:6*dk*dk,k)
      ! now convert to winds in "absolute" Cartesian components
      ucc(1:6*dk*dk,k) = uc(1:6*dk*dk)*rotpoles(1,1) + vc(1:6*dk*dk)*rotpoles(1,2) + wc(1:6*dk*dk)*rotpoles(1,3)
      vcc(1:6*dk*dk,k) = uc(1:6*dk*dk)*rotpoles(2,1) + vc(1:6*dk*dk)*rotpoles(2,2) + wc(1:6*dk*dk)*rotpoles(2,3)
      wcc(1:6*dk*dk,k) = uc(1:6*dk*dk)*rotpoles(3,1) + vc(1:6*dk*dk)*rotpoles(3,2) + wc(1:6*dk*dk)*rotpoles(3,3)
    end do  ! k loop  
    ! interpolate all required arrays to new C-C positions
    ! do not need to do map factors and Coriolis on target grid
    call fill_cc4_gather(ucc, mask_a)
    call fill_cc4_gather(vcc, mask_a)
    call fill_cc4_gather(wcc, mask_a)
  end if    ! dk>0  
  call doints4_gather(ucc, uct)
  call doints4_gather(vcc, vct)
  call doints4_gather(wcc, wct)
end if
  
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
      
implicit none

integer ier
real, dimension(:), intent(out) :: varout
real, dimension(fwsize) :: ucc
character(len=*), intent(in) :: vname
      
if ( iotest ) then
  ! read without interpolation or redistribution
  call histrd1(iarchi,ier,vname,ik,varout,ifull)
else if ( fnresid==1 ) then
  ! use bcast method for single input file
  ! requires interpolation and redistribution
  call histrd1(iarchi,ier,vname,ik,ucc,6*ik*ik,nogather=.false.)
  call doints1_gather(ucc, varout)
else
  ! use RMA method for multiple input files
  ! requires interpolation and redistribution
  call histrd1(iarchi,ier,vname,ik,ucc,6*ik*ik,nogather=.true.)
  call doints1_nogather(ucc, varout)
end if ! iotest

return
end subroutine gethist1

! This version reads, fills and interpolates a surface field
subroutine fillhist1(vname,varout,mask_a,filllimit)
      
use cc_mpi             ! CC MPI routines
use darcdf_m           ! Netcdf data
use infile             ! Input file routines
use newmpar_m          ! Grid parameters
      
implicit none
      
integer ier
real, intent(in), optional :: filllimit
real, dimension(:), intent(out) :: varout
real, dimension(fwsize) :: ucc
logical, dimension(:), intent(in) :: mask_a
character(len=*), intent(in) :: vname
      
if ( iotest ) then
  ! read without interpolation or redistribution
  call histrd1(iarchi,ier,vname,ik,varout,ifull)
else if ( fnresid==1 ) then
  ! use bcast method for single input file
  ! requires interpolation and redistribution
  call histrd1(iarchi,ier,vname,ik,ucc,6*ik*ik,nogather=.false.)
  if ( present(filllimit) ) then
    where ( ucc(:)>=filllimit )
      ucc(:) = 999.
    end where
  end if  
  call fill_cc1_gather(ucc,mask_a)
  call doints1_gather(ucc, varout)
else
  ! use RMA method for multiple input files
  ! requires interpolation and redistribution
  call histrd1(iarchi,ier,vname,ik,ucc,6*ik*ik,nogather=.true.)
  if ( present(filllimit) ) then
    where ( ucc(:)>=filllimit )
      ucc(:) = 999.
    end where
  end if  
  call fill_cc1_nogather(ucc,mask_a)
  call doints1_nogather(ucc, varout)
end if ! iotest
      
return
end subroutine fillhist1

! This version reads, fills and interpolates a surface velocity field
subroutine fillhistuv1o(uname,vname,uarout,varout,mask_a)
   
use cc_mpi             ! CC MPI routines
use darcdf_m           ! Netcdf data
use infile             ! Input file routines
use newmpar_m          ! Grid parameters
      
implicit none
      
integer ier
real, dimension(:), intent(out) :: uarout, varout
real, dimension(fwsize) :: ucc, vcc
logical, dimension(:), intent(in) :: mask_a
character(len=*), intent(in) :: uname, vname
      
if ( iotest ) then
  ! read without interpolation or redistribution
  call histrd1(iarchi,ier,uname,ik,uarout,ifull)
  call histrd1(iarchi,ier,vname,ik,varout,ifull)
else if ( fnresid==1 ) then
  ! use bcast method for single input file
  ! requires interpolation and redistribution
  call histrd1(iarchi,ier,uname,ik,ucc,6*ik*ik,nogather=.false.)
  call histrd1(iarchi,ier,vname,ik,vcc,6*ik*ik,nogather=.false.)
  call interpcurrent1(uarout,varout,ucc,vcc,mask_a,nogather=.false.)
else
  ! use RMA method for multiple input files
  ! requires interpolation and redistribution
  call histrd1(iarchi,ier,uname,ik,ucc,6*ik*ik,nogather=.true.)
  call histrd1(iarchi,ier,vname,ik,vcc,6*ik*ik,nogather=.true.)
  call interpcurrent1(uarout,varout,ucc,vcc,mask_a,nogather=.true.)
end if ! iotest
      
return
end subroutine fillhistuv1o

! This version reads 3D fields
subroutine gethist4(vname,varout,kx)

use cc_mpi             ! CC MPI routines
use darcdf_m           ! Netcdf data
use infile             ! Input file routines
use newmpar_m          ! Grid parameters
      
implicit none

integer, intent(in) :: kx
integer ier
real, dimension(:,:), intent(out) :: varout
real, dimension(fwsize,kx) :: ucc
character(len=*), intent(in) :: vname

if ( iotest ) then
  ! read without interpolation or redistribution
  call histrd4(iarchi,ier,vname,ik,kx,varout,ifull)
else if ( fnresid==1 ) then
  ! use bcast method for single input file
  ! requires interpolation and redistribution
  call histrd4(iarchi,ier,vname,ik,kx,ucc,6*ik*ik,nogather=.false.)
  call doints4_gather(ucc, varout)
else
  ! use RMA method for multiple input files
  ! requires interpolation and redistribution
  call histrd4(iarchi,ier,vname,ik,kx,ucc,6*ik*ik,nogather=.true.)
  call doints4_nogather(ucc,varout)
end if ! iotest
      
return
end subroutine gethist4   

! This version reads and interpolates 3D atmospheric fields
subroutine gethist4a(vname,varout,vmode,levkin,t_a_lev)
      
use cc_mpi               ! CC MPI routines
use darcdf_m             ! Netcdf data
use infile               ! Input file routines
use newmpar_m            ! Grid parameters
      
implicit none

integer, intent(in) :: vmode
integer, intent(in), optional :: levkin
integer ier
real, dimension(:,:), intent(out) :: varout
real, dimension(:), intent(out), optional :: t_a_lev
real, dimension(fwsize,kk) :: ucc
real, dimension(ifull,kk) :: u_k
character(len=*), intent(in) :: vname

if ( iotest ) then
  ! read without interpolation or redistribution
  call histrd4(iarchi,ier,vname,ik,kk,u_k,ifull)
else
  if ( fnresid==1 ) then
    ! use bcast method for single input file
    ! requires interpolation and redistribution
    call histrd4(iarchi,ier,vname,ik,kk,ucc,6*ik*ik,nogather=.false.)
    if ( fwsize>0.and.present(levkin).and.present(t_a_lev) ) then
      t_a_lev(:) = ucc(:,levkin)   ! store for psl calculation
    end if
    call doints4_gather(ucc, u_k)
  else
    ! use RMA method for multiple input files
    ! requires interpolation and redistribution
    call histrd4(iarchi,ier,vname,ik,kk,ucc,6*ik*ik,nogather=.true.)
    if ( fwsize>0.and.present(levkin).and.present(t_a_lev) ) then
      t_a_lev(:) = ucc(:,levkin)   ! store for psl calculation  
    end if
    call doints4_nogather(ucc, u_k)      
  end if
end if ! iotest

! vertical interpolation
call vertint(u_k,varout,vmode,kk,sigin)
      
return
end subroutine gethist4a  

! This version reads and interpolates 3D atmospheric winds
subroutine gethistuv4a(uname,vname,uarout,varout,umode,vmode)

use cc_mpi             ! CC MPI routines
use darcdf_m           ! Netcdf data
use infile             ! Input file routines
use newmpar_m          ! Grid parameters

implicit none

integer, intent(in) :: umode, vmode
integer ier
real, dimension(:,:), intent(out) :: uarout, varout
real, dimension(fwsize,kk) :: ucc, vcc
real, dimension(ifull,kk) :: u_k, v_k
character(len=*), intent(in) :: uname, vname

if ( iotest ) then
  ! read without interpolation or redistribution
  call histrd4(iarchi,ier,uname,ik,kk,u_k,ifull)
  call histrd4(iarchi,ier,vname,ik,kk,v_k,ifull)
else if ( fnresid==1 ) then
  ! use bcast method for single input file
  ! requires interpolation and redistribution
  call histrd4(iarchi,ier,uname,ik,kk,ucc,6*ik*ik,nogather=.false.)
  call histrd4(iarchi,ier,vname,ik,kk,vcc,6*ik*ik,nogather=.false.)
  call interpwind4(u_k,v_k,ucc,vcc,nogather=.false.)
else
  ! use RMA method for multiple input files
  ! requires interpolation and redistribution
  call histrd4(iarchi,ier,uname,ik,kk,ucc,6*ik*ik,nogather=.true.)
  call histrd4(iarchi,ier,vname,ik,kk,vcc,6*ik*ik,nogather=.true.)
  call interpwind4(u_k,v_k,ucc,vcc,nogather=.true.)
end if ! iotest

! vertical interpolation
call vertint(u_k,uarout,umode,kk,sigin)
call vertint(v_k,varout,vmode,kk,sigin)
      
return
end subroutine gethistuv4a  

! This version reads, fills a 3D field for the ocean
subroutine fillhist4(vname,varout,kx,mask_a,filllimit)
  
use cc_mpi             ! CC MPI routines
use darcdf_m           ! Netcdf data
use infile             ! Input file routines
use newmpar_m          ! Grid parameters
      
implicit none
      
integer, intent(in) :: kx
integer ier
real, intent(in), optional :: filllimit
real, dimension(:,:), intent(out) :: varout
real, dimension(fwsize,kx) :: ucc
logical, dimension(:), intent(in) :: mask_a
character(len=*), intent(in) :: vname

if ( iotest ) then
  ! read without interpolation or redistribution
  call histrd4(iarchi,ier,vname,ik,kx,varout,ifull)
else if ( fnresid==1 ) then
  ! use bcast method for single input file
  ! requires interpolation and redistribution
  call histrd4(iarchi,ier,vname,ik,kx,ucc,6*ik*ik,nogather=.false.)
  if ( present(filllimit) ) then
    where ( ucc(:,:)>=filllimit )
      ucc(:,:) = 999.
    end where
  end if
  call fill_cc4_gather(ucc,mask_a)
  call doints4_gather(ucc, varout)
else
  ! use RMA method for multiple input files
  ! requires interpolation and redistribution
  call histrd4(iarchi,ier,vname,ik,kx,ucc,6*ik*ik,nogather=.true.)
  if ( present(filllimit) ) then
    where ( ucc(:,:)>=filllimit )
      ucc(:,:) = 999.
    end where
  end if
  call fill_cc4_nogather(ucc,mask_a)
  call doints4_nogather(ucc, varout)
end if ! iotest

return
end subroutine fillhist4

! This version reads, fills and interpolates a 3D field for the ocean
subroutine fillhist4o(vname,varout,mask_a,bath)
   
use cc_mpi             ! CC MPI routines
use darcdf_m           ! Netcdf data
use infile             ! Input file routines
use mlo                ! Ocean physics and prognostic arrays
use newmpar_m          ! Grid parameters
      
implicit none
      
integer ier
real, dimension(:,:), intent(out) :: varout
real, dimension(fwsize,ok) :: ucc
real, dimension(ifull,ok) :: u_k
real, dimension(ifull), intent(in) :: bath
logical, dimension(:), intent(in) :: mask_a
character(len=*), intent(in) :: vname

if ( iotest ) then
  ! read without interpolation or redistribution
  call histrd4(iarchi,ier,vname,ik,ok,u_k,ifull)
else if ( fnresid==1 ) then
  ! use bcast method for single input file
  ! requires interpolation and redistribution
  call histrd4(iarchi,ier,vname,ik,ok,ucc,6*ik*ik,nogather=.false.)
  call fill_cc4_gather(ucc,mask_a)
  call doints4_gather(ucc,u_k)
else
  ! use RMA method for multiple input files
  ! requires interpolation and redistribution
  call histrd4(iarchi,ier,vname,ik,ok,ucc,6*ik*ik,nogather=.true.)
  call fill_cc4_nogather(ucc,mask_a)
  call doints4_nogather(ucc,u_k)
end if ! iotest

! vertical interpolation
varout=0.
call mloregrid(ok,bath,u_k,varout,0)

return
end subroutine fillhist4o

! This version reads, fills and interpolates ocean currents
subroutine fillhistuv4o(uname,vname,uarout,varout,mask_a,bath)
  
use cc_mpi             ! CC MPI routines
use darcdf_m           ! Netcdf data
use infile             ! Input file routines
use mlo                ! Ocean physics and prognostic arrays
use newmpar_m          ! Grid parameters
      
implicit none
      
integer ier
real, dimension(:,:), intent(out) :: uarout, varout
real, dimension(fwsize,ok) :: ucc, vcc
real, dimension(ifull,ok) :: u_k, v_k
real, dimension(ifull), intent(in) :: bath
logical, dimension(:), intent(in) :: mask_a
character(len=*), intent(in) :: uname, vname

if ( iotest ) then
  ! read without interpolation or redistribution
  call histrd4(iarchi,ier,uname,ik,ok,u_k,ifull)
  call histrd4(iarchi,ier,vname,ik,ok,v_k,ifull)
else if ( fnresid==1 ) then
  ! use bcast method for single input file
  ! requires interpolation and redistribution
  call histrd4(iarchi,ier,uname,ik,ok,ucc,6*ik*ik,nogather=.false.)
  call histrd4(iarchi,ier,vname,ik,ok,vcc,6*ik*ik,nogather=.false.)
  call interpcurrent4(u_k,v_k,ucc,vcc,mask_a,nogather=.false.)
else
  ! use RMA method for multiple input files
  ! requires interpolation and redistribution
  call histrd4(iarchi,ier,uname,ik,ok,ucc,6*ik*ik,nogather=.true.)
  call histrd4(iarchi,ier,vname,ik,ok,vcc,6*ik*ik,nogather=.true.)
  call interpcurrent4(u_k,v_k,ucc,vcc,mask_a,nogather=.true.)
end if ! iotest

! vertical interpolation
call mloregrid(ok,bath,u_k,uarout,0)
call mloregrid(ok,bath,v_k,varout,0)

return
end subroutine fillhistuv4o

! *****************************************************************************
! FILE DATA MESSAGE PASSING ROUTINES

! Define RMA windows for distributing file data to processors
subroutine file_wininit

use cc_mpi             ! CC MPI routines
use infile             ! Input file routines
use newmpar_m          ! Grid parameters

implicit none

integer, dimension(:,:,:,:), pointer :: procarray
integer, dimension(:,:,:,:), allocatable, target :: procarray_dummy
integer, dimension(4) :: shsize
integer procarray_win

if ( allocated(filemap) ) then
  deallocate( filemap )
end if
if ( allocated(axs_w) ) then
  deallocate( axs_w, ays_w, azs_w )
  deallocate( bxs_w, bys_w, bzs_w )
end if

! No RMA window for single input file
if ( fnresid<=1 ) return

if ( myid==0 ) then
  write(6,*) "Create map for file RMA windows"
end if

#ifdef usempi3
shsize(1) = ik + 4
shsize(2) = ik + 4
shsize(3) = npanels + 1
shsize(4) = 2
call ccmpi_allocshdata(procarray,shsize(1:4),procarray_win)
call ccmpi_shepoch(procarray_win)
if ( node_myid==0 ) then
  call file_wininit_defineprocarray(procarray)
end if
call ccmpi_shepoch(procarray_win)
#else
allocate( procarray_dummy(ik+4,ik+4,npanels+1,2) )
procarray => procarray_dummy
call file_wininit_defineprocarray(procarray)
#endif

call file_wininit_definefilemap(procarray)

! Distribute fields for vector rotation
if ( myid==0 ) then
  write(6,*) "Distribute vector rotation data to processors reading input files"
end if
    
allocate(axs_w(fwsize), ays_w(fwsize), azs_w(fwsize))
allocate(bxs_w(fwsize), bys_w(fwsize), bzs_w(fwsize))
if ( myid==0 ) then
  call file_distribute(axs_w,axs_a)
  call file_distribute(ays_w,ays_a)
  call file_distribute(azs_w,azs_a)
  call file_distribute(bxs_w,bxs_a)
  call file_distribute(bys_w,bys_a)
  call file_distribute(bzs_w,bzs_a)
else if ( fwsize>0 ) then
  call file_distribute(axs_w)
  call file_distribute(ays_w)
  call file_distribute(azs_w)
  call file_distribute(bxs_w)
  call file_distribute(bys_w)
  call file_distribute(bzs_w)
end if

! Define halo indices for ccmpi_filebounds
if ( myid==0 ) then
  write(6,*) "Setup bounds function for processors reading input files"
end if

call ccmpi_filebounds_setup(procarray,comm_ip,ik)

#ifdef usempi3
call ccmpi_freeshdata(procarray_win)
#else
deallocate( procarray_dummy )
#endif

if ( myid==0 ) then
  write(6,*) "Finished creating control data for file RMA windows"
end if

return
end subroutine file_wininit

subroutine file_wininit_defineprocarray(procarray)

use cc_mpi             ! CC MPI routines
use infile             ! Input file routines
use newmpar_m          ! Grid parameters

implicit none

integer i, n
integer n_n, n_e, n_s, n_w
integer ip, ipf, jpf, no, ca, cb
integer, dimension(-1:ik+2,-1:ik+2,0:npanels,2), intent(inout) :: procarray ! can be a pointer

! define host process of each input file gridpoint
procarray(-1:ik+2,-1:ik+2,0:npanels,1:2) = -1
do ipf = 0,fnproc/fnresid-1
  do jpf = 1,fnresid
    ip = ipf*fnresid + jpf - 1
    do n = 0,pnpan-1
      no = n - pnoff(ip) + 1
      ca = pioff(ip,no)
      cb = pjoff(ip,no)
      procarray(1+ca:pil+ca,1+cb:pjl+cb,no,1) = jpf - 1 ! processor rank
      procarray(1+ca:pil+ca,1+cb:pjl+cb,no,2) = ipf + 1 ! file rank
    end do
  end do
end do

! update boundaries
do n = 0,npanels
  if ( mod(n,2)==0 ) then
    n_w = mod(n+5, 6)
    n_e = mod(n+2, 6)
    n_n = mod(n+1, 6)
    n_s = mod(n+4, 6)
    do i = 1,ik
      procarray(0,i,n,:)    = procarray(ik,i,n_w,:)
      procarray(-1,i,n,:)   = procarray(ik-1,i,n_w,:)
      procarray(ik+1,i,n,:) = procarray(ik+1-i,1,n_e,:)
      procarray(ik+2,i,n,:) = procarray(ik+1-i,2,n_e,:)
      procarray(i,ik+1,n,:) = procarray(i,1,n_n,:)
      procarray(i,ik+2,n,:) = procarray(i,2,n_n,:)
      procarray(i,0,n,:)    = procarray(ik,ik+1-i,n_s,:)
      procarray(i,-1,n,:)   = procarray(ik-1,ik+1-i,n_s,:)
    end do ! i
    procarray(-1,0,n,:)      = procarray(ik,2,n_w,:)        ! wws
    procarray(0,-1,n,:)      = procarray(ik,ik-1,n_s,:)     ! wss
    procarray(0,0,n,:)       = procarray(ik,1,n_w,:)        ! ws
    procarray(ik+1,0,n,:)    = procarray(ik,1,n_e,:)        ! es  
    procarray(ik+2,0,n,:)    = procarray(ik-1,1,n_e,:)      ! ees 
    procarray(-1,ik+1,n,:)   = procarray(ik,ik-1,n_w,:)     ! wwn
    procarray(0,ik+2,n,:)    = procarray(ik-1,ik,n_w,:)     ! wnn
    procarray(ik+2,ik+1,n,:) = procarray(2,1,n_e,:)         ! een  
    procarray(ik+1,ik+2,n,:) = procarray(1,2,n_e,:)         ! enn  
    procarray(0,ik+1,n,:)    = procarray(ik,ik,n_w,:)       ! wn  
    procarray(ik+1,ik+1,n,:) = procarray(1,1,n_e,:)         ! en  
    procarray(ik+1,-1,n,:)   = procarray(ik,2,n_e,:)        ! ess  
  else
    n_w = mod(n+4, 6)
    n_e = mod(n+1, 6)
    n_n = mod(n+2, 6)
    n_s = mod(n+5, 6)
    do i = 1,ik
      procarray(0,i,n,:)    = procarray(ik+1-i,ik,n_w,:)
      procarray(-1,i,n,:)   = procarray(ik+1-i,ik-1,n_w,:)
      procarray(ik+1,i,n,:) = procarray(1,i,n_e,:)
      procarray(ik+2,i,n,:) = procarray(2,i,n_e,:)
      procarray(i,ik+1,n,:) = procarray(1,ik+1-i,n_n,:)
      procarray(i,ik+2,n,:) = procarray(2,ik+1-i,n_n,:)
      procarray(i,0,n,:)    = procarray(i,ik,n_s,:)
      procarray(i,-1,n,:)   = procarray(i,ik-1,n_s,:)
    end do ! i
    procarray(-1,0,n,:)      = procarray(ik-1,ik,n_w,:)    ! wws
    procarray(0,-1,n,:)      = procarray(2,ik,n_s,:)       ! wss
    procarray(0,0,n,:)       = procarray(ik,ik,n_w,:)      ! ws
    procarray(ik+1,0,n,:)    = procarray(1,1,n_e,:)        ! es
    procarray(ik+2,0,n,:)    = procarray(1,2,n_e,:)        ! ees
    procarray(-1,ik+1,n,:)   = procarray(2,ik,n_w,:)       ! wwn   
    procarray(0,ik+2,n,:)    = procarray(1,ik-1,n_w,:)     ! wnn  
    procarray(ik+2,ik+1,n,:) = procarray(1,ik-1,n_e,:)     ! een  
    procarray(ik+1,ik+2,n,:) = procarray(2,ik,n_e,:)       ! enn  
    procarray(0,ik+1,n,:)    = procarray(1,ik,n_w,:)       ! wn  
    procarray(ik+1,ik+1,n,:) = procarray(1,ik,n_e,:)       ! en  
    procarray(ik+1,-1,n,:)   = procarray(2,1,n_e,:)        ! ess          
  end if     ! if mod(n,2)==0 ..else..
end do       ! n

return
end subroutine file_wininit_defineprocarray

subroutine file_wininit_definefilemap(procarray)

use cc_mpi             ! CC MPI routines
use newmpar_m          ! Grid parameters
use parm_m             ! Model configuration

implicit none

integer mm, iq, idel, jdel, n
integer ncount, iproc, rproc
integer, dimension(-1:ik+2,-1:ik+2,0:npanels,2), intent(in) :: procarray
logical, dimension(-1:nproc-1) :: lproc

! calculate which grid points and input files are needed by this processor
lproc(-1:nproc-1) = .false.
do mm = 1,m_fly
  do iq = 1,ifull
    idel = int(xg4(iq,mm))
    jdel = int(yg4(iq,mm))
    n = nface4(iq,mm)
    ! search stencil of bi-cubic interpolation
    lproc(procarray(idel,  jdel+2,n,1)) = .true.
    lproc(procarray(idel+1,jdel+2,n,1)) = .true.
    lproc(procarray(idel-1,jdel+1,n,1)) = .true.
    lproc(procarray(idel  ,jdel+1,n,1)) = .true.
    lproc(procarray(idel+1,jdel+1,n,1)) = .true.
    lproc(procarray(idel+2,jdel+1,n,1)) = .true.
    lproc(procarray(idel-1,jdel,  n,1)) = .true.
    lproc(procarray(idel  ,jdel,  n,1)) = .true.
    lproc(procarray(idel+1,jdel,  n,1)) = .true.
    lproc(procarray(idel+2,jdel,  n,1)) = .true.
    lproc(procarray(idel,  jdel-1,n,1)) = .true.
    lproc(procarray(idel+1,jdel-1,n,1)) = .true.
  end do
end do
if ( lproc(-1) ) then
  write(6,*) "ERROR: Internal error in file_wininit"
  call ccmpi_abort(-1)
end if

! Construct a map of files to be accessed by MPI_Get
ncount = count(lproc)
allocate( filemap(ncount) )
ncount = 0
do iproc = 0,nproc-1
  ! stagger reading of windows - does this make any difference with active RMA?
  rproc = modulo( myid+iproc, nproc )
  if ( lproc(rproc) ) then
    ncount = ncount + 1
    filemap(ncount) = rproc
  end if
end do

return
end subroutine file_wininit_definefilemap

! Define commuication group for broadcasting file panel data
subroutine splitface

use cc_mpi            ! CC MPI routines
use newmpar_m         ! Grid parameters

implicit none

integer n, colour

! Free any existing comm_face
if ( bcst_allocated ) then
  do n = 0,npanels
    call ccmpi_commfree(comm_face(n))
  end do
  bcst_allocated = .false.
end if

! No split face for multiple input files
if ( fnresid>1 ) return

if ( myid == 0 ) then
  write(6,*) "Create communication groups for Bcast method in onthefly"
end if

do n = 0,npanels
  if ( nfacereq(n) ) then
    colour = 1
  else
    colour = -1 ! undefined
  end if
  call ccmpi_commsplit(comm_face(n),comm_world,colour,myid)
end do
bcst_allocated = .true.

if ( myid==0 ) then
  write(6,*) "Finished initalising Bcast method for onthefly"
end if

return
end subroutine splitface

subroutine processdatestring(datestring,kdate_rsav,ktime_rsav)

use cc_mpi            ! CC MPI routines

implicit none

integer, intent(out) :: kdate_rsav, ktime_rsav
integer iposa, iposb, ierx
integer yyyy, mm, dd, hh, mt
character(len=*), intent(in) :: datestring

if ( datestring(1:7)/='minutes' ) then
  write(6,*) "ERROR: Time units expected to be minutes"
  write(6,*) "Found ",trim(datestring)
  call ccmpi_abort(-1)
end if

! process year
iposa = index(trim(datestring),'since')
iposa = iposa + 5 ! skip 'since'
iposb = index(trim(datestring(iposa:)),'-')
iposb = iposa + iposb - 2 ! remove '-'
read(datestring(iposa:iposb),FMT=*,iostat=ierx) yyyy
if ( ierx/=0 ) then
  write(6,*) "ERROR reading time units.  Expecting year but found ",datestring(iposa:iposb)
  call ccmpi_abort(-1)
end if

! process month
iposa = iposb + 2 ! skip '-'
iposb = index(trim(datestring(iposa:)),'-')
iposb = iposa + iposb - 2 ! remove '-'
read(datestring(iposa:iposb),FMT=*,iostat=ierx) mm
if ( ierx/=0 ) then
  write(6,*) "ERROR reading time units.  Expecting month but found ",datestring(iposa:iposb)
  call ccmpi_abort(-1)
end if

! process day
iposa = iposb + 2 ! skip '-'
iposb = index(trim(datestring(iposa:)),' ')
iposb = iposa + iposb - 2 ! remove ' '
read(datestring(iposa:iposb),FMT=*,iostat=ierx) dd
if ( ierx/=0 ) then
  write(6,*) "ERROR reading time units.  Expecting day but found ",datestring(iposa:iposb)
  call ccmpi_abort(-1)
end if

! process hour
iposa = iposb + 2 ! skip ' '
iposb = index(trim(datestring(iposa:)),':')
iposb = iposa + iposb - 2 ! remove ':'
read(datestring(iposa:iposb),FMT=*,iostat=ierx) hh
if ( ierx/=0 ) then
  write(6,*) "ERROR reading time units.  Expecting hour but found ",datestring(iposa:iposb)
  call ccmpi_abort(-1)
end if

! process mins
iposa = iposb + 2 ! skip ':'
iposb = index(trim(datestring(iposa:)),':')
iposb = iposa + iposb - 2 ! remove ':'
read(datestring(iposa:iposb),FMT=*,iostat=ierx) mt
if ( ierx/=0 ) then
  write(6,*) "ERROR reading time units.  Expecting minutes but found ",datestring(iposa:iposb)
  call ccmpi_abort(-1)
end if

! final date and time
kdate_rsav = yyyy*10000 + mm*100 + dd
ktime_rsav = hh*100 + mt

return
end subroutine processdatestring

end module onthefly_m
