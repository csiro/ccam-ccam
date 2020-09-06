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
    
module cc_mpi

! This module manages all MPI communications between processors.  The system was originally developed
! by Martin Dix and subsequently modified by Marcus Thatcher.  Thanks to Aaron McDonough for developing
! the Vampir trace routines and upgrading the timer calls.  Thanks to Paul Ryan for the design of the
! shared memory arrays.


! Preprocessor directives:

! -Dusempi_mod is for users that want to link to MPI using Fortran 90 bindings
! -Dusempi3    exploits MPI-3 shared memory that is currently used for sharing global arrays within a node to reduce
!              the size of the CCAM memory footprint
! -Dvampir     is for coupling with VAMPIR for tracers
! -Di8r8       is for running in double precision mode
! -Dscm        is for Single Column Model which disables MPI

   use cc_omp
#ifndef scm
#ifdef usempi_mod
   use mpi
#else
   use mpif_m
#endif
   use newmpar_m

   implicit none

   private

   integer, save, public :: comm_world                                     ! global communication group
   integer, save, public :: myid                                           ! processor rank for comm_world
   integer, save, public :: ipan, jpan                                     ! grid size on processor
   integer, save, public :: ioff, joff, noff                               ! offset of processor grid relative to global grid
   integer, save, public :: nxproc, nyproc                                 ! number of processors in the x and y directions
   integer, parameter, private :: nagg = 3                                 ! maximum number of levels to aggregate
   integer, parameter, public :: maxcolour = 2                             ! maximum number of colours for iterative solvers
   
   integer, save, private :: maxbuflen, maxvertlen                         ! bounds buffer size   
   logical, save, public :: uniform_decomp                                 ! uniform decomposition flag
   logical, save, public :: mydiag                                         ! true if diagnostic point id, jd is in my region
   !$acc declare create(mydiag,ipan,jpan)
   
   integer, save, public :: comm_node, node_myid, node_nproc               ! node communicator
   integer, save, public :: comm_nodecaptain, nodecaptain_myid, &
                            nodecaptain_nproc                              ! node captain communicator
   integer, save, public :: node_captainid                                 ! rank of the node captain in the comm_nodecaptain group
   
   integer, save, public :: comm_proc, comm_rows, comm_cols                ! comm groups for scale-selective filter
   integer, save, public :: hproc, mproc, npta, pprocn, pprocx             ! decomposition parameters for scale-selective filter

   integer, save, public :: comm_vnode, vnode_nproc, vnode_myid            ! procformat communicator for node
   integer, save, public :: comm_vleader, vleader_nproc, vleader_myid      ! procformat communicator for node captain group   
   integer, save, public :: vnode_vleaderid                                ! rank of the procformat node captain

   integer(kind=4), save, private :: nreq, rreq                            ! number of messages requested and to be received
   integer(kind=4), allocatable, dimension(:), save, private :: ireq       ! requested message index
   integer, allocatable, dimension(:), save, private :: rlist              ! map of processor index from requested message index
   
   integer, allocatable, dimension(:), save, public :: neighlist           ! list of neighbour processors
   integer, allocatable, dimension(:), save, private :: neighmap           ! map of processor to neighbour index
   integer, save, public :: neighnum                                       ! number of neighbours

   integer(kind=4), allocatable, dimension(:), save, public ::              &
      specmap_recv                                                         ! gather map recived for spectral filter
   integer(kind=4), allocatable, dimension(:), save, public ::              &
      specmap_send                                                         ! gather map sent for spectral filter
   integer, allocatable, dimension(:), save, public :: specmap_ext         ! gather map for spectral filter (includes filter final
                                                                           ! pass for sparse arrays)
   type globalpack_info
     real, allocatable, dimension(:,:,:) :: localdata
   end type globalpack_info
   type(globalpack_info), allocatable, dimension(:,:,:), save, private ::   & 
      globalpack                                                           ! store sparse global arrays for spectral filter

   integer, save, public :: pipan, pjpan, pnpan                            ! decomposition parameters file window
   integer, save, public :: pil_g, pjl_g, pka_g, pko_g                     ! decomposition parameters file window
   integer, save, public :: fnproc, fnresid, fncount, mynproc              ! number and decomposition of input files
   integer, allocatable, dimension(:), save, public :: pnoff               ! file window panel offset
   integer, allocatable, dimension(:,:), save, public :: pioff, pjoff      ! file window coordinate offset
   integer(kind=4), allocatable, dimension(:), save, public ::              &
      filemap_recv, filemap_rmod                                           ! file map received for onthefly
   integer(kind=4), allocatable, dimension(:), save, public ::              &
      filemap_send, filemap_smod                                           ! file map sent for onthefly
   integer, allocatable, dimension(:), save, public :: filemap_facecomm    ! communicator for face bcast
   integer, allocatable, dimension(:), save, public :: filemap_rinv        ! inverse map for filemap_recv/rmod
   integer, save, public :: fnproc_bcast_max = 6                           ! max number of input files for bcast method
   
   integer, allocatable, dimension(:), save, private :: fileneighlist      ! list of file neighbour processors
   integer, save, public :: fileneighnum                                   ! number of file neighbours
   
   public :: ccmpi_setup, ccmpi_distribute, ccmpi_gather, ccmpi_gatherr8,   &
             ccmpi_distributer8, ccmpi_gatherall, bounds, boundsuv,         &
             deptsync, intssync_send, intssync_recv, start_log, end_log,    &
             log_on, log_off, log_flush, log_setup,                         &
             ccglobal_posneg, readglobvar, writeglobvar, ccmpi_reduce,      &
             ccmpi_reducer8, ccmpi_allreduce, ccmpi_abort, ccmpi_bcast,     &
             ccmpi_bcastr8, ccmpi_barrier, ccmpi_gatherx, ccmpi_gatherxr8,  &
             ccmpi_scatterx, ccmpi_allgatherx, ccmpi_init,                  &
             ccmpi_finalize, ccmpi_commsplit, ccmpi_commfree,               &
             bounds_colour_send, bounds_colour_recv, boundsr8,              &
             ccmpi_reinit, ccmpi_alltoall, ccmpi_procformat_init
   public :: mgbounds, mgcollect, mgbcast, mgbcastxn, mgbcasta, mg_index,   &
             mg_fproc, mg_fproc_1
   public :: ind, indx, indp, indg, iq2iqg, indv_mpi, indglobal, fproc,     &
             face_set, uniform_set, dix_set
   public :: allocateglobalpack, copyglobalpack,                            &
             ccmpi_gathermap_send2, ccmpi_gathermap_recv2,                  &
             ccmpi_gathermap_send3, ccmpi_gathermap_recv3, getglobalpack_v, &
             setglobalpack_v
   public :: ccmpi_filewinget, ccmpi_filewinunpack, ccmpi_filebounds_setup, &
             ccmpi_filebounds_send, ccmpi_filebounds_recv,                  &
             ccmpi_filegather, ccmpi_filedistribute, procarray
#ifdef usempi3
   public :: ccmpi_allocshdata, ccmpi_allocshdatar8, ccmpi_remap
   public :: ccmpi_shepoch, ccmpi_freeshdata
#endif
   
   interface ccmpi_gather
      module procedure host_gather2, host_gather3, host_gather4
      module procedure proc_gather2, proc_gather3, proc_gather4
   end interface
   interface ccmpi_gatherr8
      module procedure host_gather2r8, host_gather3r8, host_gather4r8
      module procedure proc_gather2r8, proc_gather3r8, proc_gather4r8
   end interface
   interface ccmpi_distribute
      module procedure host_distribute2, host_distribute2i
      module procedure proc_distribute2, proc_distribute2i
      module procedure host_distribute3, host_distribute3i
      module procedure proc_distribute3, proc_distribute3i
      module procedure host_distribute4
      module procedure proc_distribute4
   end interface
   interface ccmpi_distributer8
      module procedure host_distribute2r8, host_distribute3r8, host_distribute4r8
      module procedure proc_distribute2r8, proc_distribute3r8, proc_distribute4r8
   end interface
   interface ccmpi_gatherall
      module procedure ccmpi_gatherall2, ccmpi_gatherall3
   end interface
   interface bounds
      module procedure bounds2, bounds3, bounds4
   end interface
   interface boundsuv
      module procedure boundsuv2, boundsuv3
   end interface
   interface boundsr8
      module procedure bounds3r8
   end interface
   interface intssync_send
      module procedure intssync_send3, intssync_send4
   end interface
   interface intssync_recv
      module procedure intssync_recv3, intssync_recv4
   end interface
   interface ccglobal_posneg
      module procedure ccglobal_posneg2, ccglobal_posneg3, ccglobal_posneg4
      module procedure ccglobal_posneg3o, ccglobal_posneg4o
   end interface
   interface readglobvar
      module procedure readglobvar2, readglobvar3, readglobvar2i
   end interface
   interface writeglobvar
      module procedure writeglobvar2, writeglobvar3
   end interface
   interface ccmpi_reduce
      module procedure ccmpi_reduce2i, ccmpi_reduce1r, ccmpi_reduce2r, ccmpi_reduce3r
      module procedure ccmpi_reduce2c, ccmpi_reduce2l
   end interface
   interface ccmpi_reducer8
     module procedure ccmpi_reduce1rr8
   end interface
   interface ccmpi_allreduce
      module procedure ccmpi_allreduce1i, ccmpi_allreduce2i, ccmpi_allreduce1r, ccmpi_allreduce2r
      module procedure ccmpi_allreduce3r, ccmpi_allreduce1c, ccmpi_allreduce2c
   end interface
   interface ccmpi_bcast
      module procedure ccmpi_bcast1i, ccmpi_bcast2i, ccmpi_bcast3i, ccmpi_bcast1r, ccmpi_bcast2r
      module procedure ccmpi_bcast3r, ccmpi_bcast4r, ccmpi_bcast5r, ccmpi_bcast1s
   end interface
   interface ccmpi_bcastr8
      module procedure ccmpi_bcast2r8, ccmpi_bcast3r8, ccmpi_bcast4r8
   end interface
   interface ccmpi_gatherx
      module procedure ccmpi_gatherx2r,  ccmpi_gatherx3r,  ccmpi_gatherx4r
      module procedure ccmpi_gatherx23r, ccmpi_gatherx34r, ccmpi_gatherx45r
      module procedure ccmpi_gatherx2i,  ccmpi_gatherx3i
      module procedure ccmpi_gatherx23l
   end interface
   interface ccmpi_gatherxr8
     module procedure ccmpi_gatherx2rr8,  ccmpi_gatherx3rr8
     module procedure ccmpi_gatherx23rr8, ccmpi_gatherx34rr8, ccmpi_gatherx45rr8
   end interface
   interface ccmpi_scatterx
      module procedure ccmpi_scatterx2r, ccmpi_scatterx32r, ccmpi_scatterx3r
      module procedure ccmpi_scatterx32l
   end interface
   interface ccmpi_allgatherx
      module procedure ccmpi_allgatherx2i, ccmpi_allgatherx2r
   end interface
   interface ccmpi_alltoall
      module procedure ccmpi_alltoall2l
   end interface
   interface ccmpi_filebounds_send
      module procedure ccmpi_filebounds_send2, ccmpi_filebounds_send3
   end interface ccmpi_filebounds_send
   interface ccmpi_filebounds_recv
      module procedure ccmpi_filebounds_recv2, ccmpi_filebounds_recv3
   end interface ccmpi_filebounds_recv
   interface ccmpi_filegather
      module procedure host_filegather2, proc_filegather2
      module procedure host_filegather3, proc_filegather3
   end interface ccmpi_filegather
   interface ccmpi_filedistribute
      module procedure host_filedistribute2, proc_filedistribute2
      module procedure host_filedistribute3, proc_filedistribute3
   end interface ccmpi_filedistribute
   interface ccmpi_filewinget
     module procedure ccmpi_filewinget2, ccmpi_filewinget3
   end interface
   interface mgbounds
      module procedure mgbounds2, mgbounds3
   end interface
   interface mgbcast
      module procedure mgbcast2, mgbcast3
   end interface
   interface mgbcasta
      module procedure mgbcasta2, mgbcasta3
   end interface
   interface mgcollect
      module procedure mgcollect1, mgcollectreduce, mgcollectxn
   end interface
#ifdef usempi3
   interface ccmpi_allocshdata
      module procedure ccmpi_allocshdata2r, ccmpi_allocshdata3r, ccmpi_allocshdata4r
      module procedure ccmpi_allocshdata5r
      module procedure ccmpi_allocshdata2i, ccmpi_allocshdata3i, ccmpi_allocshdata5i
   end interface
   interface ccmpi_allocshdatar8
      module procedure ccmpi_allocshdata2_r8, ccmpi_allocshdata3_r8, ccmpi_allocshdata4_r8
   end interface
#endif
   
   ! Do directions need to be swapped
   logical, parameter, private, dimension(0:npanels) ::                      &
           swap_e = (/ .true., .false., .true., .false., .true., .false. /), &
           swap_w = (/ .false., .true., .false., .true., .false., .true. /), &
           swap_n = (/ .false., .true., .false., .true., .false., .true. /), &
           swap_s = (/ .true., .false., .true., .false., .true., .false. /)

   type bounds_info
      ! Buffer arrays 
      real, dimension(:), allocatable :: sbuf, rbuf
      ! Flag for whether u and v need to be reversed
      real, dimension(:), allocatable :: send_neg, uv_neg
      ! Buffer arrays
      real(kind=8), dimension(:), allocatable :: s8buf, r8buf
      ! Index arrays
      integer, dimension(:), allocatable :: request_list, send_list, unpack_list
      integer, dimension(:), allocatable :: request_list_uv, send_list_uv, unpack_list_uv
      ! Flag for whether u and v need to be swapped
      logical, dimension(:), allocatable :: uv_swap, send_swap
      ! Buffer length
      integer :: len, sbuflen, rbuflen
      ! Number of points for each processor. Also double row versions.
      ! lenx is first row plux corner points.  lenh is just the ne side.
      integer :: slen2, rlen2
      integer, dimension(maxcolour) :: rlenh_bg, rlenh_fn
      integer, dimension(maxcolour) :: rlen_bg, rlen_fn
      integer, dimension(maxcolour) :: rlenx_bg, rlenx_fn
      integer, dimension(maxcolour) :: slenh_bg, slenh_fn
      integer, dimension(maxcolour) :: slen_bg, slen_fn
      integer, dimension(maxcolour) :: slenx_bg, slenx_fn
      ! Vector groups
      integer :: rlen_su_bg, rlen_ev_fn, rlen_nnu_bg, rlen_eev_fn
      integer :: rlen_sv_bg, rlen_wu_fn, rlen_nv_bg, rlen_eu_fn
      integer :: rlen_ssv_bg, rlen_wwu_fn, rlen_nnv_bg, rlen_eeu_fn
      integer :: slen_su_bg, slen_ev_fn, slen_nnu_bg, slen_eev_fn
      integer :: slen_sv_bg, slen_wu_fn, slen_nv_bg, slen_eu_fn
      integer :: slen_ssv_bg, slen_wwu_fn, slen_nnv_bg, slen_eeu_fn
   end type bounds_info

   ! bounds data
   type(bounds_info), allocatable, dimension(:), save :: bnds

   ! partition boundary indices into colours
   integer, dimension(:,:), allocatable, save, public :: iqx, iqn, iqe, iqw, iqs
   integer, dimension(:,:), allocatable, save, public :: iqne, iqen, iqse, iqes, iqnw, iqwn, iqsw, iqws
   integer, public, save :: ifull_maxcolour
   integer, dimension(maxcolour), public, save :: ifull_colour, ifull_colour_border

   ! flag whether processor region edge is a face edge.
   logical, public, save :: edge_w, edge_n, edge_s, edge_e
   
   type dpoints_info
      real, dimension(:,:), allocatable :: a
      real, dimension(:), allocatable :: b
   end type dpoints_info
   type dindex_info
      integer, dimension(:,:), allocatable :: a
   end type dindex_info
   type sextra_info
      real, dimension(:), allocatable :: a
   end type sextra_info
   
   ! Off processor departure points
   type(dpoints_info), allocatable, dimension(:), public, save :: dpoints ! request list from other proc
   type(dpoints_info), allocatable, dimension(:), private, save :: dbuf   ! recv buffer
   type(dindex_info), allocatable, dimension(:), private, save :: dindex  ! request list for my proc
   type(sextra_info), allocatable, dimension(:), public, save :: sextra   ! send buffer
   ! Number of points for each processor.
   integer, dimension(:), allocatable, private, save :: dslen
   integer, dimension(:), allocatable, public, save :: drlen

   ! Multi-grid arrays
   type mgtype
      integer :: ifull, iextra, ifull_fine, ifull_coarse
      integer :: merge_len, merge_row, ipan, merge_pos, nmax
      integer :: comm_merge, neighnum, npanx
      integer, dimension(:), allocatable :: merge_list
      integer, dimension(:), allocatable :: in, ie, is, iw
      integer, dimension(:), allocatable :: ine, ien, ise, ies, inw, iwn, isw, iws
      integer, dimension(:), allocatable :: coarse_a, coarse_b, coarse_c
      integer, dimension(:), allocatable :: fine, fine_n, fine_e, fine_ne
      integer, dimension(:), allocatable :: neighlist
      integer, dimension(:), allocatable :: procmap
      real, dimension(:), allocatable :: zzn, zze, zzs, zzw, zz
   end type mgtype

   type mgbndtype
      integer :: len, rlen, slen, rlenx, slenx
      integer, dimension(:), allocatable :: send_list
      integer, dimension(:), allocatable :: unpack_list
      integer, dimension(:), allocatable :: request_list
   end type mgbndtype

   ! Multi-grid levels and buffers
   type(mgtype), dimension(:), allocatable, save, public :: mg
   type(mgbndtype), dimension(:,:), allocatable, save, public :: mg_bnds
   integer, save, public :: mg_maxlevel, mg_maxlevel_local
   integer, save, public :: mg_ifull_maxcolour
   integer, dimension(:,:), allocatable, save, public :: col_iq, col_iqn, col_iqe, col_iqs, col_iqw
   integer, dimension(:,:), allocatable, save, public :: col_iqne, col_iqen, col_iqse, col_iqes
   integer, dimension(:,:), allocatable, save, public :: col_iqnw, col_iqwn, col_iqsw, col_iqws

   ! File IO
   type filebounds_info
      integer, dimension(:,:), allocatable :: send_list
      integer, dimension(:,:), allocatable :: request_list
      integer, dimension(:,:), allocatable :: unpack_list
      integer :: len, rlen, slen, rlenx, slenx
   end type filebounds_info
   
   type(filebounds_info), allocatable, dimension(:), save :: filebnds
   
   ! Timer
   integer, public, save :: ints_begin, ints_end
   integer, public, save :: nonlin_begin, nonlin_end
   integer, public, save :: helm_begin, helm_end
   integer, public, save :: adjust_begin, adjust_end
   integer, public, save :: upglobal_begin, upglobal_end
   integer, public, save :: hordifg_begin, hordifg_end
   integer, public, save :: vadv_begin, vadv_end
   integer, public, save :: depts_begin, depts_end
   integer, public, save :: stag_begin, stag_end
   integer, public, save :: ocnstag_begin, ocnstag_end
   integer, public, save :: mfix_begin, mfix_end
   integer, public, save :: phys_begin, phys_end
   integer, public, save :: outfile_begin, outfile_end
   integer, public, save :: onthefly_begin, onthefly_end
   integer, public, save :: otf_fill_begin, otf_fill_end
   integer, public, save :: otf_ints1_begin, otf_ints1_end
   integer, public, save :: otf_ints4_begin, otf_ints4_end
   integer, public, save :: histrd3_begin, histrd3_end
   integer, public, save :: histrd4_begin, histrd4_end
   integer, public, save :: histrd5_begin, histrd5_end
   integer, public, save :: indata_begin, indata_end
   integer, public, save :: nestin_begin, nestin_end
   integer, public, save :: nestOTF_begin, nestOTF_end
   integer, public, save :: nestWIN_begin, nestWIN_end
   integer, public, save :: nestcalc_begin, nestcalc_end
   integer, public, save :: nestcomm_begin, nestcomm_end
   integer, public, save :: ensemble_begin, ensemble_end
   integer, public, save :: amipsst_begin, amipsst_end
   integer, public, save :: gwdrag_begin, gwdrag_end
   integer, public, save :: convection_begin, convection_end
   integer, public, save :: cloud_begin, cloud_end
   integer, public, save :: radnet_begin, radnet_end
   integer, public, save :: radinit_begin, radinit_end
   integer, public, save :: sfluxnet_begin, sfluxnet_end
   integer, public, save :: sfluxwater_begin, sfluxwater_end
   integer, public, save :: sfluxland_begin, sfluxland_end
   integer, public, save :: sfluxurban_begin, sfluxurban_end
   integer, public, save :: vertmix_begin, vertmix_end
   integer, public, save :: aerosol_begin, aerosol_end
   integer, public, save :: maincalc_begin, maincalc_end
   integer, public, save :: gatherfile_begin, gatherfile_end
   integer, public, save :: distributefile_begin, distributefile_end
   integer, public, save :: allreducepn_begin, allreducepn_end   
   integer, public, save :: precon_begin, precon_end
   integer, public, save :: waterdynamics_begin, waterdynamics_end
   integer, public, save :: waterunpack_begin, waterunpack_end
   integer, public, save :: watermfix_begin, watermfix_end
   integer, public, save :: waterdeps_begin, waterdeps_end
   integer, public, save :: watereos_begin, watereos_end
   integer, public, save :: waterints_begin, waterints_end
   integer, public, save :: watervadv_begin, watervadv_end
   integer, public, save :: waterhelm_begin, waterhelm_end
   integer, public, save :: wateriadv_begin, wateriadv_end
   integer, public, save :: waterdiff_begin, waterdiff_end
   integer, public, save :: waterpack_begin, waterpack_end
   integer, public, save :: river_begin, river_end
   integer, public, save :: bcast_begin, bcast_end
   integer, public, save :: alltoall_begin, alltoall_end
   integer, public, save :: allgather_begin, allgather_end
   integer, public, save :: allgathercc_begin, allgathercc_end
   integer, public, save :: gather_begin, gather_end
   integer, public, save :: gathercc_begin, gathercc_end
   integer, public, save :: scatter_begin, scatter_end
   integer, public, save :: scattercc_begin, scattercc_end
   integer, public, save :: scatterfile_begin, scatterfile_end
   integer, public, save :: reduce_begin, reduce_end
   integer, public, save :: allreduce_begin, allreduce_end
   integer, public, save :: mpiwait_begin, mpiwait_end
   integer, public, save :: mpiwaitmap_begin, mpiwaitmap_end
   integer, public, save :: mpiwaitdep_begin, mpiwaitdep_end
   integer, public, save :: mpiwaitmg_begin, mpiwaitmg_end
   integer, public, save :: mpiwaitfile_begin, mpiwaitfile_end
   integer, public, save :: mpiwaitmapfile_begin, mpiwaitmapfile_end
   integer, public, save :: gathermg_begin, gathermg_end
   integer, public, save :: bcastmg_begin, bcastmg_end
   integer, public, save :: mgsetup_begin, mgsetup_end
   integer, public, save :: mgfine_begin, mgfine_end
   integer, public, save :: mgup_begin, mgup_end
   integer, public, save :: mgcoarse_begin, mgcoarse_end
   integer, public, save :: mgdown_begin, mgdown_end
   integer, public, save :: p1_begin, p1_end
   integer, public, save :: p2_begin, p2_end
   integer, public, save :: p3_begin, p3_end
   integer, public, save :: p4_begin, p4_end
   integer, public, save :: p5_begin, p5_end
   integer, public, save :: p6_begin, p6_end
   integer, parameter :: nevents = 86
#ifdef simple_timer
   public :: simple_timer_finalize
   !real(kind=8), dimension(nevents), save :: tot_time = 0._8, start_time
   integer(kind=8), dimension(nevents), save :: tot_time = 0_8, start_time
#endif
   character(len=15), dimension(nevents), save :: event_name
   real, save, public :: mpiinit_time, total_time

#ifdef vampir
#include "vt_user.inc"
#endif


contains

   subroutine ccmpi_setup(id,jd,idjd,dt)
      use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer
      use indices_m
      use latlong_m
      use map_m
      use sumdd_m
      use vecsuv_m
      use workglob_m
      use xyzinfo_m
      integer, intent(in) :: id, jd
      integer, intent(out) :: idjd
      integer iproc, dproc, iq, iqg, i, j, n
      integer(kind=4) :: ierr, colour, rank, lcommin, lcommout
      integer, dimension(ifull) :: colourmask
      real, intent(in) :: dt
      real, dimension(:,:), allocatable :: dum, dum_g
      real, dimension(:,:), allocatable :: dumu, dumv
      real(kind=8), dimension(:,:), allocatable :: dumr8, dumr8_g
      logical(kind=4) :: ltrue

      nreq = 0
      allocate( bnds(0:nproc-1) )
      
      
      ! Decompose grid over processes
      if ( uniform_decomp ) then
         call proc_setup_uniform(id,jd,idjd)
         ! may require two boundries from the same process
         maxbuflen = (max(ipan,jpan)+4)*2*2*2*npan
      else
         call proc_setup(id,jd,idjd)
         if ( nproc < npanels+1 ) then
            ! possible to have two boundaries from the same process 
            maxbuflen = (il_g+4)*2*2*2*npan + 4
         else
            ! only one boundary can be sent from a process 
            maxbuflen = (max(ipan,jpan)+4)*2*2 + 4
         end if    
      end if
      maxvertlen = max( kl, ol, 15 )
!$acc update device(mydiag,ipan,jpan)
      
      
      ! Distribute global arrays over processes
      if ( myid == 0 ) then
         allocate( dum_g(ifull_g,15), dum(ifull,15) ) 
         dum_g(1:ifull_g,1)  = wts_g(1:ifull_g)
         dum_g(1:ifull_g,2)  = em_g(1:ifull_g)
         dum_g(1:ifull_g,3)  = emu_g(1:ifull_g)
         dum_g(1:ifull_g,4)  = emv_g(1:ifull_g)
         dum_g(1:ifull_g,5)  = ax_g(1:ifull_g)
         dum_g(1:ifull_g,6)  = ay_g(1:ifull_g)
         dum_g(1:ifull_g,7)  = az_g(1:ifull_g)
         dum_g(1:ifull_g,8)  = bx_g(1:ifull_g)
         dum_g(1:ifull_g,9)  = by_g(1:ifull_g)
         dum_g(1:ifull_g,10) = bz_g(1:ifull_g)
         dum_g(1:ifull_g,11) = f_g(1:ifull_g)
         dum_g(1:ifull_g,12) = fu_g(1:ifull_g)
         dum_g(1:ifull_g,13) = fv_g(1:ifull_g)
         dum_g(1:ifull_g,14) = rlatt_g(1:ifull_g)
         dum_g(1:ifull_g,15) = rlongg_g(1:ifull_g)
         call ccmpi_distribute(dum(:,1:15),dum_g(:,1:15)) 
         wts(1:ifull)    = dum(1:ifull,1)
         em(1:ifull)     = dum(1:ifull,2)
         emu(1:ifull)    = dum(1:ifull,3)
         emv(1:ifull)    = dum(1:ifull,4)
         ax(1:ifull)     = dum(1:ifull,5)
         ay(1:ifull)     = dum(1:ifull,6)
         az(1:ifull)     = dum(1:ifull,7)
         bx(1:ifull)     = dum(1:ifull,8)
         by(1:ifull)     = dum(1:ifull,9)
         bz(1:ifull)     = dum(1:ifull,10)
         f(1:ifull)      = dum(1:ifull,11)
         fu(1:ifull)     = dum(1:ifull,12)
         fv(1:ifull)     = dum(1:ifull,13)
         rlatt(1:ifull)  = dum(1:ifull,14)
         rlongg(1:ifull) = dum(1:ifull,15)
         dum_g(1:ifull_g,1:4) = rlat4(1:ifull_g,1:4)
         dum_g(1:ifull_g,5:8) = rlong4(1:ifull_g,1:4)
         call ccmpi_distribute(dum(:,1:8),dum_g(:,1:8))
         rlat4_l(1:ifull,1:4)  = dum(1:ifull,1:4)
         rlong4_l(1:ifull,1:4) = dum(1:ifull,5:8)
         deallocate( dum_g, dum )
         allocate( dumr8_g(ifull_g,3), dumr8(ifull,3) )
         dumr8_g(1:ifull_g,1) = x_g(1:ifull_g)
         dumr8_g(1:ifull_g,2) = y_g(1:ifull_g)
         dumr8_g(1:ifull_g,3) = z_g(1:ifull_g)
         call ccmpi_distributer8(dumr8(:,1:3),dumr8_g(:,1:3))
         x(1:ifull) = dumr8(1:ifull,1)
         y(1:ifull) = dumr8(1:ifull,2)
         z(1:ifull) = dumr8(1:ifull,3)
         deallocate( dumr8_g, dumr8 )
      else
         allocate( dum(ifull,15) )
         call ccmpi_distribute(dum(:,1:15))
         wts(1:ifull)    = dum(1:ifull,1)
         em(1:ifull)     = dum(1:ifull,2)
         emu(1:ifull)    = dum(1:ifull,3)
         emv(1:ifull)    = dum(1:ifull,4)
         ax(1:ifull)     = dum(1:ifull,5)
         ay(1:ifull)     = dum(1:ifull,6)
         az(1:ifull)     = dum(1:ifull,7)
         bx(1:ifull)     = dum(1:ifull,8)
         by(1:ifull)     = dum(1:ifull,9)
         bz(1:ifull)     = dum(1:ifull,10)
         f(1:ifull)      = dum(1:ifull,11)
         fu(1:ifull)     = dum(1:ifull,12)
         fv(1:ifull)     = dum(1:ifull,13)
         rlatt(1:ifull)  = dum(1:ifull,14)
         rlongg(1:ifull) = dum(1:ifull,15)
         call ccmpi_distribute(dum(:,1:8))
         rlat4_l(1:ifull,1:4)  = dum(1:ifull,1:4)
         rlong4_l(1:ifull,1:4) = dum(1:ifull,5:8)
         deallocate( dum )         
         allocate( dumr8(ifull,3) )
         call ccmpi_distributer8(dumr8(:,1:3))
         x(1:ifull) = dumr8(1:ifull,1)
         y(1:ifull) = dumr8(1:ifull,2)
         z(1:ifull) = dumr8(1:ifull,3)
         deallocate( dumr8 )
      end if

      
      ! Configure halos      
      call bounds_setup(dt)
      allocate( dum(1:ifull+iextra,2) )
      dum = 0.
      dum(1:ifull,1) = em(1:ifull)
      dum(1:ifull,2) = f(1:ifull)
      call bounds(dum(:,1:2),corner=.true.)
      em(ifull+1:ifull+iextra) = dum(ifull+1:ifull+iextra,1)
      f(ifull+1:ifull+iextra) = dum(ifull+1:ifull+iextra,2)
      deallocate( dum )
      allocate( dumu(ifull+iextra,4), dumv(ifull+iextra,4) )
      dumu = 0.
      dumv = 0.
      dumu(1:ifull,1) = emu(1:ifull)
      dumv(1:ifull,1) = emv(1:ifull)
      dumu(1:ifull,2) = ax(1:ifull)
      dumv(1:ifull,2) = bx(1:ifull)
      dumu(1:ifull,3) = ay(1:ifull)
      dumv(1:ifull,3) = by(1:ifull)
      dumu(1:ifull,4) = az(1:ifull)
      dumv(1:ifull,4) = bz(1:ifull)
      call boundsuv(dumu(:,1:4),dumv(:,1:4))
      emu(ifull+1:ifull+iextra) = dumu(ifull+1:ifull+iextra,1)     
      emv(ifull+1:ifull+iextra) = dumv(ifull+1:ifull+iextra,1)
      ax(ifull+1:ifull+iextra)  = dumu(ifull+1:ifull+iextra,2)     
      bx(ifull+1:ifull+iextra)  = dumv(ifull+1:ifull+iextra,2)
      ay(ifull+1:ifull+iextra)  = dumu(ifull+1:ifull+iextra,3)     
      by(ifull+1:ifull+iextra)  = dumv(ifull+1:ifull+iextra,3)
      az(ifull+1:ifull+iextra)  = dumu(ifull+1:ifull+iextra,4)     
      bz(ifull+1:ifull+iextra)  = dumv(ifull+1:ifull+iextra,4)
      deallocate( dumu, dumv )
      
      
      ! Off processor departure points
      allocate( dpoints(neighnum) )
      allocate( dbuf(0:neighnum) )
      allocate( dindex(0:neighnum) )
      allocate( sextra(neighnum) )
      allocate( dslen(0:neighnum), drlen(0:neighnum) )
      dslen(:) = 0
      drlen(:) = 0
      do dproc = 1,neighnum
        iproc = neighlist(dproc)
        allocate( dpoints(dproc)%a(4,bnds(iproc)%len) )
        allocate( dbuf(dproc)%a(4,bnds(iproc)%len) )
        allocate( dbuf(dproc)%b(nagg*bnds(iproc)%len) )
        allocate( dindex(dproc)%a(bnds(iproc)%len,2) )
        allocate( sextra(dproc)%a(nagg*bnds(iproc)%len) )
      end do
      ! store invalid points in dproc=0
      allocate( dbuf(0)%a(4,1) )
      allocate( dbuf(0)%b(1) )
      allocate( dindex(0)%a(1,2) )


      ! Pack colour indices
      do n = 1,npan
         do j = 1,jpan
            do i = 1,ipan
               iq  = indp(i,j,n)  ! Local
               iqg = indg(i,j,n)  ! Global
               colourmask(iq) = findcolour(iqg)
            end do
         end do
      end do

      do n = 1,maxcolour
         ifull_colour(n) = count( colourmask == n )
      end do
      ifull_maxcolour = maxval( ifull_colour )
      allocate( iqx(ifull_maxcolour,maxcolour) )
      allocate( iqn(ifull_maxcolour,maxcolour), iqe(ifull_maxcolour,maxcolour) )
      allocate( iqw(ifull_maxcolour,maxcolour), iqs(ifull_maxcolour,maxcolour) )
      allocate( iqne(ifull_maxcolour,maxcolour), iqen(ifull_maxcolour,maxcolour) )
      allocate( iqse(ifull_maxcolour,maxcolour), iqes(ifull_maxcolour,maxcolour) )
      allocate( iqnw(ifull_maxcolour,maxcolour), iqwn(ifull_maxcolour,maxcolour) )
      allocate( iqsw(ifull_maxcolour,maxcolour), iqws(ifull_maxcolour,maxcolour) )
      ifull_colour = 0
      ! first process border
      do n = 1,npan
        j = 1
        do i = 1,ipan
          iq = indp(i,j,n)
          ifull_colour(colourmask(iq)) = ifull_colour(colourmask(iq)) + 1
          iqx(ifull_colour(colourmask(iq)),colourmask(iq)) = iq
          iqn(ifull_colour(colourmask(iq)),colourmask(iq)) = in(iq)
          iqe(ifull_colour(colourmask(iq)),colourmask(iq)) = ie(iq)
          iqw(ifull_colour(colourmask(iq)),colourmask(iq)) = iw(iq)
          iqs(ifull_colour(colourmask(iq)),colourmask(iq)) = is(iq)
          iqne(ifull_colour(colourmask(iq)),colourmask(iq)) = ine(iq)
          iqen(ifull_colour(colourmask(iq)),colourmask(iq)) = ien(iq)
          iqse(ifull_colour(colourmask(iq)),colourmask(iq)) = ise(iq)
          iqes(ifull_colour(colourmask(iq)),colourmask(iq)) = ies(iq)
          iqnw(ifull_colour(colourmask(iq)),colourmask(iq)) = inw(iq)
          iqwn(ifull_colour(colourmask(iq)),colourmask(iq)) = iwn(iq)
          iqsw(ifull_colour(colourmask(iq)),colourmask(iq)) = isw(iq)
          iqws(ifull_colour(colourmask(iq)),colourmask(iq)) = iws(iq)
        end do
        j = jpan
        do i = 1,ipan
          iq = indp(i,j,n)
          ifull_colour(colourmask(iq)) = ifull_colour(colourmask(iq)) + 1
          iqx(ifull_colour(colourmask(iq)),colourmask(iq)) = iq
          iqn(ifull_colour(colourmask(iq)),colourmask(iq)) = in(iq)
          iqe(ifull_colour(colourmask(iq)),colourmask(iq)) = ie(iq)
          iqw(ifull_colour(colourmask(iq)),colourmask(iq)) = iw(iq)
          iqs(ifull_colour(colourmask(iq)),colourmask(iq)) = is(iq)
          iqne(ifull_colour(colourmask(iq)),colourmask(iq)) = ine(iq)
          iqen(ifull_colour(colourmask(iq)),colourmask(iq)) = ien(iq)
          iqse(ifull_colour(colourmask(iq)),colourmask(iq)) = ise(iq)
          iqes(ifull_colour(colourmask(iq)),colourmask(iq)) = ies(iq)
          iqnw(ifull_colour(colourmask(iq)),colourmask(iq)) = inw(iq)
          iqwn(ifull_colour(colourmask(iq)),colourmask(iq)) = iwn(iq)
          iqsw(ifull_colour(colourmask(iq)),colourmask(iq)) = isw(iq)
          iqws(ifull_colour(colourmask(iq)),colourmask(iq)) = iws(iq)
        end do
        i = 1
        do j = 2,jpan-1
          iq = indp(i,j,n)
          ifull_colour(colourmask(iq)) = ifull_colour(colourmask(iq)) + 1
          iqx(ifull_colour(colourmask(iq)),colourmask(iq)) = iq
          iqn(ifull_colour(colourmask(iq)),colourmask(iq)) = in(iq)
          iqe(ifull_colour(colourmask(iq)),colourmask(iq)) = ie(iq)
          iqw(ifull_colour(colourmask(iq)),colourmask(iq)) = iw(iq)
          iqs(ifull_colour(colourmask(iq)),colourmask(iq)) = is(iq)
          iqne(ifull_colour(colourmask(iq)),colourmask(iq)) = ine(iq)
          iqen(ifull_colour(colourmask(iq)),colourmask(iq)) = ien(iq)
          iqse(ifull_colour(colourmask(iq)),colourmask(iq)) = ise(iq)
          iqes(ifull_colour(colourmask(iq)),colourmask(iq)) = ies(iq)
          iqnw(ifull_colour(colourmask(iq)),colourmask(iq)) = inw(iq)
          iqwn(ifull_colour(colourmask(iq)),colourmask(iq)) = iwn(iq)
          iqsw(ifull_colour(colourmask(iq)),colourmask(iq)) = isw(iq)
          iqws(ifull_colour(colourmask(iq)),colourmask(iq)) = iws(iq)
        end do
        i = ipan
        do j = 2,jpan-1
          iq = indp(i,j,n)
          ifull_colour(colourmask(iq)) = ifull_colour(colourmask(iq)) + 1
          iqx(ifull_colour(colourmask(iq)),colourmask(iq)) = iq
          iqn(ifull_colour(colourmask(iq)),colourmask(iq)) = in(iq)
          iqe(ifull_colour(colourmask(iq)),colourmask(iq)) = ie(iq)
          iqw(ifull_colour(colourmask(iq)),colourmask(iq)) = iw(iq)
          iqs(ifull_colour(colourmask(iq)),colourmask(iq)) = is(iq)
          iqne(ifull_colour(colourmask(iq)),colourmask(iq)) = ine(iq)
          iqen(ifull_colour(colourmask(iq)),colourmask(iq)) = ien(iq)
          iqse(ifull_colour(colourmask(iq)),colourmask(iq)) = ise(iq)
          iqes(ifull_colour(colourmask(iq)),colourmask(iq)) = ies(iq)
          iqnw(ifull_colour(colourmask(iq)),colourmask(iq)) = inw(iq)
          iqwn(ifull_colour(colourmask(iq)),colourmask(iq)) = iwn(iq)
          iqsw(ifull_colour(colourmask(iq)),colourmask(iq)) = isw(iq)
          iqws(ifull_colour(colourmask(iq)),colourmask(iq)) = iws(iq)
        end do
      end do
      ifull_colour_border(1:maxcolour) = ifull_colour(1:maxcolour)
      ! next process interior
      do n = 1,npan
        do j = 2,jpan-1
          do i = 2,ipan-1
            iq = indp(i,j,n)
            ifull_colour(colourmask(iq)) = ifull_colour(colourmask(iq)) + 1
            iqx(ifull_colour(colourmask(iq)),colourmask(iq)) = iq
            iqn(ifull_colour(colourmask(iq)),colourmask(iq)) = in(iq)
            iqe(ifull_colour(colourmask(iq)),colourmask(iq)) = ie(iq)
            iqw(ifull_colour(colourmask(iq)),colourmask(iq)) = iw(iq)
            iqs(ifull_colour(colourmask(iq)),colourmask(iq)) = is(iq)
            iqne(ifull_colour(colourmask(iq)),colourmask(iq)) = ine(iq)
            iqen(ifull_colour(colourmask(iq)),colourmask(iq)) = ien(iq)
            iqse(ifull_colour(colourmask(iq)),colourmask(iq)) = ise(iq)
            iqes(ifull_colour(colourmask(iq)),colourmask(iq)) = ies(iq)
            iqnw(ifull_colour(colourmask(iq)),colourmask(iq)) = inw(iq)
            iqwn(ifull_colour(colourmask(iq)),colourmask(iq)) = iwn(iq)
            iqsw(ifull_colour(colourmask(iq)),colourmask(iq)) = isw(iq)
            iqws(ifull_colour(colourmask(iq)),colourmask(iq)) = iws(iq)
          end do
        end do
      end do


      ! Create MPI_SUMDR for calculating global sums with high precision
      ltrue = .true. 
      ! Operator MPI_SUMDR is created based on an external function DRPDR.
      call MPI_OP_CREATE( DRPDR, ltrue, MPI_SUMDR, ierr )
      
      
      ! prepare comm groups - used by scale-selective filter
      if ( uniform_decomp ) then
         npta = 6                     ! number of panels per processor
         mproc = nproc                ! number of processors per panel
         pprocn = 0                   ! start panel
         pprocx = 5                   ! end panel
         hproc = 0                    ! host processor for panel
      else
         npta = max( 6/nproc, 1 )     ! number of panels per processor
         mproc = max( nproc/6, 1 )    ! number of processors per panel
         pprocn = myid*npta/mproc     ! start panel
         pprocx = pprocn + npta - 1   ! end panel
         hproc = pprocn*mproc/npta    ! host processor for panel
      end if

      ! comm between work groups with captain hproc
      colour = hproc
      rank = myid - hproc
      lcommin = comm_world
      call MPI_Comm_Split( lcommin, colour, rank, lcommout, ierr )
      comm_proc = lcommout
      
      ! comm between columns in work group
      colour = ioff
      rank = joff/jpan
      lcommin = comm_proc
      call MPI_Comm_Split( lcommin, colour, rank, lcommout, ierr )
      comm_cols = lcommout
      
      ! comm between rows in work group      
      colour = joff
      rank = ioff/ipan
      lcommin = comm_proc
      call MPI_Comm_Split( lcommin, colour, rank, lcommout, ierr )
      comm_rows = lcommout
            
   return
   end subroutine ccmpi_setup
   
   subroutine host_distribute2(af,a1)
      ! Convert standard 1D arrays to face form and distribute to processors
      real, dimension(ifull), intent(out) :: af
      real, dimension(ifull_g), intent(in) :: a1
      real, dimension(ifull,0:nproc-1) :: sbuf
      integer :: j, n, iq, iproc
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      integer(kind=4) :: ierr, lsize, lcomm
      integer :: npoff, ipoff, jpoff ! Offsets for target
      integer :: slen
      
      ! map array in order of processor rank
      if ( uniform_decomp ) then
         do iproc = 0,nproc-1
            call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,ipan,jpan)
            do n = 1,npan
               do j = 1,jpan
                  iq = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g
                  slen = (j-1)*ipan + (n-1)*ipan*jpan
                  sbuf(slen+1:slen+ipan,iproc) = a1(iq+1:iq+ipan)
               end do
            end do
         end do
      else
         do iproc = 0,nproc-1
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan)
            do n = 1,npan
               do j = 1,jpan
                  iq = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g
                  slen = (j-1)*ipan + (n-1)*ipan*jpan
                  sbuf(slen+1:slen+ipan,iproc) = a1(iq+1:iq+ipan)
               end do
            end do
         end do
      end if

      lsize = ifull
      lcomm = comm_world
      call START_LOG(scattercc_begin)
      call MPI_Scatter( sbuf, lsize, ltype, af, lsize, ltype, 0_4, lcomm, ierr )
      call END_LOG(scattercc_end)
      
   end subroutine host_distribute2

   subroutine proc_distribute2(af)
      ! Convert standard 1D arrays to face form and distribute to processors
      real, dimension(ifull), intent(out) :: af
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      integer(kind=4) :: ierr, lsize, lcomm
      real, dimension(1,1) :: sbuf

      if ( myid == 0 ) then
         write(6,*) "Error: ccmpi_distribute argument required on proc 0"
         call ccmpi_abort(-1)
      end if
      
      lsize = ifull
      lcomm = comm_world
      call START_LOG(scattercc_begin)      
      call MPI_Scatter( sbuf, lsize, ltype, af, lsize, ltype, 0_4, lcomm, ierr )
      call END_LOG(scattercc_end)
      
   end subroutine proc_distribute2

   subroutine host_distribute2r8(af,a1)
      ! Convert standard 1D arrays to face form and distribute to processors
      real(kind=8), dimension(ifull), intent(out) :: af
      real(kind=8), dimension(ifull_g), intent(in) :: a1
      integer :: j, n, iq, iproc
      integer(kind=4) :: ierr, lsize, lcomm
      real(kind=8), dimension(ifull,0:nproc-1) :: sbuf
      integer :: npoff, ipoff, jpoff ! Offsets for target
      integer :: slen
      
      ! map array in order of processor rank
      if ( uniform_decomp ) then
         do iproc = 0,nproc-1
            call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,ipan,jpan)
            do n = 1,npan
               do j = 1,jpan
                  iq = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g
                  slen = (j-1)*ipan + (n-1)*ipan*jpan
                  sbuf(slen+1:slen+ipan,iproc) = a1(iq+1:iq+ipan)
               end do
            end do
         end do
      else
         do iproc = 0,nproc-1
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan)
            do n = 1,npan
               do j = 1,jpan
                  iq = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g
                  slen = (j-1)*ipan + (n-1)*ipan*jpan
                  sbuf(slen+1:slen+ipan,iproc) = a1(iq+1:iq+ipan)
               end do
            end do
         end do
      end if

      lsize = ifull
      lcomm = comm_world
      call START_LOG(scattercc_begin)
      call MPI_Scatter( sbuf, lsize, MPI_DOUBLE_PRECISION, af, lsize, MPI_DOUBLE_PRECISION, 0_4, lcomm, ierr )
      call END_LOG(scattercc_end)
      
   end subroutine host_distribute2r8
   
   subroutine proc_distribute2r8(af)
      ! Convert standard 1D arrays to face form and distribute to processors
      real(kind=8), dimension(ifull), intent(out) :: af
      integer(kind=4) :: ierr, lsize, lcomm
      real(kind=8), dimension(1,1) :: sbuf
      
      if ( myid == 0 ) then
         write(6,*) "Error: ccmpi_distribute argument required on proc 0"
         call ccmpi_abort(-1)
      end if
      
      lsize = ifull
      lcomm = comm_world
      call START_LOG(scattercc_begin)
      call MPI_Scatter( sbuf, lsize, MPI_DOUBLE_PRECISION, af, lsize, MPI_DOUBLE_PRECISION, 0_4, lcomm, ierr )
      call END_LOG(scattercc_end)
      
   end subroutine proc_distribute2r8

   subroutine host_distribute2i(af,a1)
      ! Convert standard 1D arrays to face form and distribute to processors
      integer, dimension(ifull), intent(out) :: af
      integer, dimension(ifull_g), intent(in) :: a1
      integer :: j, n, iq, iproc
#ifdef i8r8
      integer(kind=4),parameter :: ltype = MPI_INTEGER8
#else
      integer(kind=4),parameter :: ltype = MPI_INTEGER
#endif
      integer(kind=4) :: ierr, lsize, lcomm
      integer, dimension(ifull,0:nproc-1) :: sbuf
      integer :: npoff, ipoff, jpoff ! Offsets for target
      integer :: slen
      
      ! map array in order of processor rank
      if ( uniform_decomp ) then
         do iproc = 0,nproc-1
            call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,ipan,jpan)
            do n = 1,npan
               do j = 1,jpan
                  iq = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g
                  slen = (j-1)*ipan + (n-1)*ipan*jpan
                  sbuf(slen+1:slen+ipan,iproc) = a1(iq+1:iq+ipan)
               end do
            end do
         end do
      else
         do iproc = 0,nproc-1
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan)
            do n = 1,npan
               do j = 1,jpan
                  iq = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g
                  slen = (j-1)*ipan + (n-1)*ipan*jpan
                  sbuf(slen+1:slen+ipan,iproc) = a1(iq+1:iq+ipan)
               end do
            end do
         end do
      end if

      lsize = ifull
      lcomm = comm_world
      call START_LOG(scattercc_begin)
      call MPI_Scatter( sbuf, lsize, ltype, af, lsize, ltype, 0_4, lcomm, ierr )
      call END_LOG(scattercc_end)
      
   end subroutine host_distribute2i

   subroutine proc_distribute2i(af)
      ! Convert standard 1D arrays to face form and distribute to processors
      integer, dimension(ifull), intent(out) :: af
#ifdef i8r8
      integer(kind=4),parameter :: ltype = MPI_INTEGER8
#else
      integer(kind=4),parameter :: ltype = MPI_INTEGER
#endif
      integer(kind=4) :: ierr, lsize, lcomm
      integer, dimension(1,1) :: sbuf
      
      if ( myid == 0 ) then
         write(6,*) "Error: ccmpi_distribute argument required on proc 0"
         call ccmpi_abort(-1)
      end if
      
      lsize = ifull
      lcomm = comm_world
      call START_LOG(scattercc_begin)
      call MPI_Scatter( sbuf, lsize, ltype, af, lsize, ltype, 0_4, lcomm, ierr )
      call END_LOG(scattercc_end)

   end subroutine proc_distribute2i

   subroutine host_distribute3(af,a1)
      ! Convert standard 1D arrays to face form and distribute to processors
      ! This is also used for tracers, so second dimension is not necessarily
      ! the number of levels
      real, dimension(:,:), intent(out) :: af
      real, dimension(:,:), intent(in) :: a1
      integer :: j, n, k, iq, iproc
#ifdef i8r8
      integer(kind=4),parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4),parameter :: ltype = MPI_REAL
#endif
      integer(kind=4) :: ierr, lsize, lcomm
      real, dimension(ifull,size(af,2),0:nproc-1) :: sbuf
      real, dimension(ifull,size(af,2)) :: aftemp
      integer :: npoff, ipoff, jpoff ! Offsets for target
      integer :: slen, kx
      
      kx = size(af,2)

      ! map array in order of processor rank
      if ( uniform_decomp ) then
         do k = 1,kx
            do iproc = 0,nproc-1
               call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,ipan,jpan)
               do n = 1,npan
                  do j = 1,jpan
                     iq = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g
                     slen = (j-1)*ipan + (n-1)*ipan*jpan
                     sbuf(slen+1:slen+ipan,k,iproc) = a1(iq+1:iq+ipan,k)
                  end do
               end do
            end do
         end do
      else
         do k = 1,kx
            do iproc = 0,nproc-1
               call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan)
               do n = 1,npan
                  do j = 1,jpan
                     iq = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g
                     slen = (j-1)*ipan + (n-1)*ipan*jpan
                     sbuf(slen+1:slen+ipan,k,iproc) = a1(iq+1:iq+ipan,k)
                  end do
               end do
            end do
         end do
      end if

      lsize = ifull*kx
      lcomm = comm_world
      if ( size(af,1) == ifull ) then
         call START_LOG(scattercc_begin)
         call MPI_Scatter( sbuf, lsize, ltype, af, lsize, ltype, 0_4, lcomm, ierr )
         call END_LOG(scattercc_end)
      else
         call START_LOG(scattercc_begin)
         call MPI_Scatter( sbuf, lsize, ltype, aftemp, lsize, ltype, 0_4, lcomm, ierr )  
         call END_LOG(scattercc_end)
         af(1:ifull,1:kx) = aftemp(1:ifull,1:kx)
      end if   

   end subroutine host_distribute3

   subroutine proc_distribute3(af)
      ! Convert standard 1D arrays to face form and distribute to processors
      ! This is also used for tracers, so second dimension is not necessarily
      ! the number of levels
      real, dimension(:,:), intent(out) :: af
#ifdef i8r8
      integer(kind=4),parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4),parameter :: ltype = MPI_REAL
#endif
      integer(kind=4) :: ierr, lsize, lcomm
      real, dimension(1,1,1) :: sbuf
      real, dimension(ifull,size(af,2)) :: aftemp
      integer :: kx

      if ( myid == 0 ) then
         write(6,*) "Error: ccmpi_distribute argument required on proc 0"
         call ccmpi_abort(-1)
      end if
      
      kx = size(af,2)
      lsize = ifull*kx
      lcomm = comm_world
      if ( size(af,1) == ifull ) then
         call START_LOG(scattercc_begin)
         call MPI_Scatter( sbuf, lsize, ltype, af, lsize, ltype, 0_4, lcomm, ierr )
         call END_LOG(scattercc_end)
      else
         call START_LOG(scattercc_begin)
         call MPI_Scatter( sbuf, lsize, ltype, aftemp, lsize, ltype, 0_4, lcomm, ierr ) 
         call END_LOG(scattercc_end)
         af(1:ifull,1:kx) = aftemp(1:ifull,1:kx)
      end if   

   end subroutine proc_distribute3

   subroutine host_distribute3r8(af,a1)
      ! Convert standard 1D arrays to face form and distribute to processors
      real(kind=8), dimension(:,:), intent(out) :: af
      real(kind=8), dimension(:,:), intent(in) :: a1
      integer :: j, n, iq, iproc
      integer(kind=4) :: ierr, lsize, lcomm
      real(kind=8), dimension(ifull,size(af,2),0:nproc-1) :: sbuf
      real(kind=8), dimension(ifull,size(af,2)) :: aftemp
      integer :: npoff, ipoff, jpoff ! Offsets for target
      integer :: slen, kx, k

      kx = size(af,2)
      
      ! map array in order of processor rank
      if ( uniform_decomp ) then
         do k = 1,kx 
            do iproc = 0,nproc-1
               call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,ipan,jpan)
               do n = 1,npan
                  do j = 1,jpan
                     iq = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g
                     slen = (j-1)*ipan + (n-1)*ipan*jpan
                     sbuf(slen+1:slen+ipan,k,iproc) = a1(iq+1:iq+ipan,k)
                  end do
               end do
            end do
         end do
      else
         do k = 1,kx 
            do iproc = 0,nproc-1
               call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan)
               do n = 1,npan
                  do j = 1,jpan
                     iq = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g
                     slen = (j-1)*ipan + (n-1)*ipan*jpan
                     sbuf(slen+1:slen+ipan,k,iproc) = a1(iq+1:iq+ipan,k)
                  end do
               end do
            end do
         end do
      end if

      lsize = ifull*kx
      lcomm = comm_world
      if ( size(af,1) == ifull ) then
         call START_LOG(scattercc_begin)
         call MPI_Scatter( sbuf, lsize, MPI_DOUBLE_PRECISION, af, lsize, MPI_DOUBLE_PRECISION, 0_4, lcomm, ierr )
         call END_LOG(scattercc_end)
      else
         call START_LOG(scattercc_begin) 
         call MPI_Scatter( sbuf, lsize, MPI_DOUBLE_PRECISION, aftemp, lsize, MPI_DOUBLE_PRECISION, 0_4, lcomm, ierr )
         call END_LOG(scattercc_end)
         af(1:ifull,1:kx) = aftemp(1:ifull,1:kx)
      end if    

   end subroutine host_distribute3r8
   
   subroutine proc_distribute3r8(af)
      ! Convert standard 1D arrays to face form and distribute to processors
      real(kind=8), dimension(:,:), intent(out) :: af
      integer :: kx
      integer(kind=4) :: ierr, lsize, lcomm
      real(kind=8), dimension(1,1,1) :: sbuf
      real(kind=8), dimension(ifull,size(af,2)) :: aftemp

      if ( myid == 0 ) then
         write(6,*) "Error: ccmpi_distribute argument required on proc 0"
         call ccmpi_abort(-1)
      end if
      
      kx = size(af,2)
      lsize = ifull*kx
      lcomm = comm_world
      if ( size(af,1) == ifull ) then
         call START_LOG(scattercc_begin)
         call MPI_Scatter( sbuf, lsize, MPI_DOUBLE_PRECISION, af, lsize, MPI_DOUBLE_PRECISION, 0_4, lcomm, ierr )
         call END_LOG(scattercc_end)
      else
         call START_LOG(scattercc_begin) 
         call MPI_Scatter( sbuf, lsize, MPI_DOUBLE_PRECISION, aftemp, lsize, MPI_DOUBLE_PRECISION, 0_4, lcomm, ierr ) 
         call END_LOG(scattercc_end)
         af(1:ifull,1:kx) = aftemp(1:ifull,1:kx)
      end if   

   end subroutine proc_distribute3r8

   subroutine host_distribute3i(af,a1)
      ! Convert standard 1D arrays to face form and distribute to processors
      ! This is also used for tracers, so second dimension is not necessarily
      ! the number of levels
      integer, dimension(:,:), intent(out) :: af
      integer, dimension(:,:), intent(in) :: a1
      integer :: j, n, k, iq, iproc
      integer :: npoff, ipoff, jpoff ! Offsets for target
      integer :: slen, kx
      integer, dimension(ifull,size(af,2),0:nproc-1) :: sbuf
      integer, dimension(ifull,size(af,2)) :: aftemp
#ifdef i8r8
      integer(kind=4),parameter :: ltype = MPI_INTEGER8
#else
      integer(kind=4),parameter :: ltype = MPI_INTEGER
#endif
      integer(kind=4) :: ierr, lsize, lcomm

      kx = size(af,2)

      ! map array in order of processor rank
      if ( uniform_decomp ) then
         do k = 1,kx
            do iproc = 0,nproc-1
               call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,ipan,jpan)
               do n = 1,npan
                  do j = 1,jpan
                     iq = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g
                     slen = (j-1)*ipan + (n-1)*ipan*jpan
                     sbuf(slen+1:slen+ipan,k,iproc) = a1(iq+1:iq+ipan,k)
                  end do
               end do
            end do
         end do
      else
         do k = 1,kx
            do iproc = 0,nproc-1
               call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan)
               do n = 1,npan
                  do j = 1,jpan
                     iq = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g
                     slen = (j-1)*ipan + (n-1)*ipan*jpan
                     sbuf(slen+1:slen+ipan,k,iproc) = a1(iq+1:iq+ipan,k)
                  end do
               end do
            end do
         end do
      end if

      lsize = ifull*kx
      lcomm = comm_world
      if ( size(af,1) == ifull ) then
         call START_LOG(scattercc_begin)
         call MPI_Scatter( sbuf, lsize, ltype, af, lsize, ltype, 0_4, lcomm, ierr )
         call END_LOG(scattercc_end)
      else
         call START_LOG(scattercc_begin) 
         call MPI_Scatter( sbuf, lsize, ltype, aftemp, lsize, ltype, 0_4, lcomm, ierr )
         call END_LOG(scattercc_end)
         af(1:ifull,1:kx) = aftemp(1:ifull,1:kx)
      end if   

   end subroutine host_distribute3i

   subroutine proc_distribute3i(af)
      ! Convert standard 1D arrays to face form and distribute to processors
      ! This is also used for tracers, so second dimension is not necessarily
      ! the number of levels
      integer, dimension(:,:), intent(out) :: af
      integer, dimension(1,1,1) :: sbuf
      integer, dimension(ifull,size(af,2)) :: aftemp
      integer :: kx
#ifdef i8r8
      integer(kind=4),parameter :: ltype = MPI_INTEGER8
#else
      integer(kind=4),parameter :: ltype = MPI_INTEGER
#endif
      integer(kind=4) :: ierr, lsize, lcomm
      
      if ( myid == 0 ) then
         write(6,*) "Error: ccmpi_distribute argument required on proc 0"
         call ccmpi_abort(-1)
      end if
      
      kx = size(af,2)
      lsize = ifull*kx
      lcomm = comm_world
      if ( size(af,1) == ifull ) then
         call START_LOG(scattercc_begin)
         call MPI_Scatter( sbuf, lsize, ltype, af, lsize, ltype, 0_4, lcomm, ierr )
         call END_LOG(scattercc_end)
      else
         call START_LOG(scattercc_begin) 
         call MPI_Scatter( sbuf, lsize, ltype, aftemp, lsize, ltype, 0_4, lcomm, ierr )
         call END_LOG(scattercc_end)
         af(1:ifull,1:kx) = aftemp(1:ifull,1:kx)
      end if   

   end subroutine proc_distribute3i   
   
   subroutine host_distribute4(af,a1)
      ! Convert standard 1D arrays to face form and distribute to processors
      ! This is also used for tracers, so second dimension is not necessarily
      ! the number of levels
      real, dimension(:,:,:), intent(out) :: af
      real, dimension(:,:,:), intent(in) :: a1
      real, dimension(ifull,size(af,2),size(af,3),0:nproc-1) :: sbuf
      real, dimension(ifull,size(af,2),size(af,3)) :: aftemp
      integer :: j, n, k, l, iq, iproc
      integer :: npoff, ipoff, jpoff ! Offsets for target
      integer :: slen, kx, lx
#ifdef i8r8
      integer(kind=4),parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4),parameter :: ltype = MPI_REAL
#endif
      integer(kind=4) :: ierr, lsize, lcomm

      kx = size(af,2)
      lx = size(af,3)

      ! map array in order of processor rank
      if ( uniform_decomp ) then
         do l = 1,lx 
            do k = 1,kx
               do iproc = 0,nproc-1
                  call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,ipan,jpan)
                  do n = 1,npan
                     do j = 1,jpan
                        iq = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g
                        slen = (j-1)*ipan + (n-1)*ipan*jpan
                        sbuf(slen+1:slen+ipan,k,l,iproc) = a1(iq+1:iq+ipan,k,l)
                     end do   
                  end do
               end do
            end do
         end do
      else
         do l = 1,lx 
            do k = 1,kx
               do iproc = 0,nproc-1
                  call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan)
                  do n = 1,npan
                     do j = 1,jpan
                        iq = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g
                        slen = (j-1)*ipan + (n-1)*ipan*jpan
                        sbuf(slen+1:slen+ipan,k,l,iproc) = a1(iq+1:iq+ipan,k,l)
                     end do   
                  end do
               end do
            end do
         end do
      end if

      lsize = ifull*kx*lx
      lcomm = comm_world
      if ( size(af,1) == ifull ) then
         call START_LOG(scattercc_begin)
         call MPI_Scatter( sbuf, lsize, ltype, af, lsize, ltype, 0_4, lcomm, ierr )
         call END_LOG(scattercc_end)
      else
         call START_LOG(scattercc_begin) 
         call MPI_Scatter( sbuf, lsize, ltype, aftemp, lsize, ltype, 0_4, lcomm, ierr )
         call END_LOG(scattercc_end)
         af(1:ifull,1:kx,1:lx) = aftemp(1:ifull,1:kx,1:lx)
      end if   

   end subroutine host_distribute4

   subroutine proc_distribute4(af)
      ! Convert standard 1D arrays to face form and distribute to processors
      ! This is also used for tracers, so second dimension is not necessarily
      ! the number of levels
      real, dimension(:,:,:), intent(out) :: af
      real, dimension(1,1,1,1) :: sbuf
      real, dimension(ifull,size(af,2),size(af,3)) :: aftemp
      integer :: kx, lx
#ifdef i8r8
      integer(kind=4),parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4),parameter :: ltype = MPI_REAL
#endif
      integer(kind=4) :: ierr, lsize, lcomm
      
      if ( myid == 0 ) then
         write(6,*) "Error: ccmpi_distribute argument required on proc 0"
         call ccmpi_abort(-1)
      end if

      kx = size(af,2)
      lx = size(af,3)
      lsize = ifull*kx*lx
      lcomm = comm_world
      if ( size(af,1) == ifull ) then
         call START_LOG(scattercc_begin)
         call MPI_Scatter( sbuf, lsize, ltype, af, lsize, ltype, 0_4, lcomm, ierr )
         call END_LOG(scattercc_end)
      else
         call START_LOG(scattercc_begin) 
         call MPI_Scatter( sbuf, lsize, ltype, aftemp, lsize, ltype, 0_4, lcomm, ierr ) 
         call END_LOG(scattercc_end)
         af(1:ifull,1:kx,1:lx) = aftemp(1:ifull,1:kx,1:lx)
      end if   

   end subroutine proc_distribute4

   subroutine host_distribute4r8(af,a1)
      ! Convert standard 1D arrays to face form and distribute to processors
      real(kind=8), dimension(:,:,:), intent(out) :: af
      real(kind=8), dimension(:,:,:), intent(in) :: a1
      real(kind=8), dimension(ifull,size(af,2),size(af,3),0:nproc-1) :: sbuf
      real(kind=8), dimension(ifull,size(af,2),size(af,3)) :: aftemp
      integer :: j, n, iq, iproc
      integer :: npoff, ipoff, jpoff ! Offsets for target
      integer :: slen, kx, k, lx, l
      integer(kind=4) :: ierr, lsize, lcomm
      
      kx = size(af,2)
      lx = size(af,3)
      
      ! map array in order of processor rank
      if ( uniform_decomp ) then
         do l = 1,lx  
            do k = 1,kx 
               do iproc = 0,nproc-1
                  call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,ipan,jpan)
                  do n = 1,npan
                     do j = 1,jpan
                        iq = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g
                        slen = (j-1)*ipan + (n-1)*ipan*jpan
                        sbuf(slen+1:slen+ipan,k,l,iproc) = a1(iq+1:iq+ipan,k,l)
                     end do   
                  end do
               end do
            end do
         end do
      else
         do l = 1,lx  
            do k = 1,kx 
               do iproc = 0,nproc-1
                  call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan)
                  do n = 1,npan
                     do j = 1,jpan
                        iq = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g
                        slen = (j-1)*ipan + (n-1)*ipan*jpan
                        sbuf(slen+1:slen+ipan,k,l,iproc) = a1(iq+1:iq+ipan,k,l)
                     end do   
                  end do
               end do
            end do
         end do
      end if

      lsize = ifull*kx*lx
      lcomm = comm_world
      if ( size(af,1) == ifull ) then
         call START_LOG(scattercc_begin)
         call MPI_Scatter( sbuf, lsize, MPI_DOUBLE_PRECISION, af, lsize, MPI_DOUBLE_PRECISION, 0_4, lcomm, ierr )
         call END_LOG(scattercc_end)
      else
         call START_LOG(scattercc_begin) 
         call MPI_Scatter( sbuf, lsize, MPI_DOUBLE_PRECISION, aftemp, lsize, MPI_DOUBLE_PRECISION, 0_4, lcomm, ierr )
         call END_LOG(scattercc_end)
         af(1:ifull,1:kx,1:lx) = aftemp(1:ifull,1:kx,1:lx)
      end if

   end subroutine host_distribute4r8
   
   subroutine proc_distribute4r8(af)
      ! Convert standard 1D arrays to face form and distribute to processors
      real(kind=8), dimension(:,:,:), intent(out) :: af
      real(kind=8), dimension(1,1,1,1) :: sbuf
      real(kind=8), dimension(ifull,size(af,2),size(af,3)) :: aftemp
      integer :: kx, lx
      integer(kind=4) :: ierr, lsize, lcomm
      
      if ( myid == 0 ) then
         write(6,*) "Error: ccmpi_distribute argument required on proc 0"
         call ccmpi_abort(-1)
      end if

      kx = size(af,2)
      lx = size(af,3)
      lsize = ifull*kx*lx
      lcomm = comm_world
      if ( size(af,1) == ifull ) then
         call START_LOG(scattercc_begin)
         call MPI_Scatter( sbuf, lsize, MPI_DOUBLE_PRECISION, af, lsize, MPI_DOUBLE_PRECISION, 0_4, lcomm, ierr )
         call END_LOG(scattercc_end)
      else
         call START_LOG(scattercc_begin) 
         call MPI_Scatter( sbuf, lsize, MPI_DOUBLE_PRECISION, aftemp, lsize, MPI_DOUBLE_PRECISION, 0_4, lcomm, ierr ) 
         call END_LOG(scattercc_end)
         af(1:ifull,1:kx,1:lx) = aftemp(1:ifull,1:kx,1:lx)
      end if   

   end subroutine proc_distribute4r8

   subroutine host_gather2(a,ag)
      ! Collect global arrays.
      real, dimension(ifull), intent(in) :: a
      real, dimension(ifull_g), intent(out) :: ag
      integer :: iproc
#ifdef i8r8
      integer(kind=4),parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4),parameter :: ltype = MPI_REAL
#endif
      integer(kind=4) :: ierr, lsize, lcomm
      real, dimension(ifull,0:nproc-1) :: abuf
      integer :: ipoff, jpoff, npoff
      integer :: j, n, iq, iqg
      
      lsize = ifull
      lcomm = comm_world
      call START_LOG(gathercc_begin)
      call MPI_Gather( a, lsize, ltype, abuf, lsize, ltype, 0_4, lcomm, ierr )
      call END_LOG(gathercc_end)

      ! map array in order of processor rank
      if ( uniform_decomp ) then
         do iproc = 0,nproc-1
            call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,ipan,jpan)
            do n = 1,npan
               do j = 1,jpan
                  ! Global indices are i+ipoff, j+jpoff, n-npoff
                  iqg = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g ! True global 1D index
                  iq = (j-1)*ipan + (n-1)*ipan*jpan
                  ag(iqg+1:iqg+ipan) = abuf(iq+1:iq+ipan,iproc)
               end do
            end do
         end do
      else
         do iproc = 0,nproc-1
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan)
            do n = 1,npan
               do j = 1,jpan
                  ! Global indices are i+ipoff, j+jpoff, n-npoff
                  iqg = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g ! True global 1D index
                  iq = (j-1)*ipan + (n-1)*ipan*jpan
                  ag(iqg+1:iqg+ipan) = abuf(iq+1:iq+ipan,iproc)
               end do
            end do
         end do
      end if

   end subroutine host_gather2
   
   subroutine proc_gather2(a)
      ! Collect global arrays.
      real, dimension(ifull), intent(in) :: a
#ifdef i8r8
      integer(kind=4),parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4),parameter :: ltype = MPI_REAL
#endif
      integer(kind=4) :: ierr, lsize, lcomm
      real, dimension(1,1) :: abuf

      if ( myid == 0 ) then
         write(6,*) "Error: ccmpi_gather argument required on proc 0"
         call ccmpi_abort(-1)
      end if
      
      lsize = ifull
      lcomm = comm_world
      call START_LOG(gathercc_begin)
      call MPI_Gather( a, lsize, ltype, abuf, lsize, ltype, 0_4, lcomm, ierr )
      call END_LOG(gathercc_end)

   end subroutine proc_gather2

   subroutine host_gather2r8(a,ag)
      ! Collect global arrays.
      real(kind=8), dimension(ifull), intent(in) :: a
      real(kind=8), dimension(ifull_g), intent(out) :: ag
      integer :: iproc
      integer(kind=4),parameter :: ltype = MPI_DOUBLE_PRECISION
      integer(kind=4) :: ierr, lsize, lcomm
      real(kind=8), dimension(ifull,0:nproc-1) :: abuf
      integer :: ipoff, jpoff, npoff
      integer :: j, n, iq, iqg
      
      lsize = ifull
      lcomm = comm_world
      call START_LOG(gathercc_begin)
      call MPI_Gather( a, lsize, ltype, abuf, lsize, ltype, 0_4, lcomm, ierr )
      call END_LOG(gathercc_end)

      ! map array in order of processor rank
      if ( uniform_decomp ) then
         do iproc = 0,nproc-1
            call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,ipan,jpan)
            do n = 1,npan
               do j = 1,jpan
                  ! Global indices are i+ipoff, j+jpoff, n-npoff
                  iqg = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g ! True global 1D index
                  iq = (j-1)*ipan + (n-1)*ipan*jpan
                  ag(iqg+1:iqg+ipan) = abuf(iq+1:iq+ipan,iproc)
               end do
            end do
         end do
      else
         do iproc = 0,nproc-1
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan)
            do n = 1,npan
               do j = 1,jpan
                  ! Global indices are i+ipoff, j+jpoff, n-npoff
                  iqg = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g ! True global 1D index
                  iq = (j-1)*ipan + (n-1)*ipan*jpan
                  ag(iqg+1:iqg+ipan) = abuf(iq+1:iq+ipan,iproc)
               end do
            end do
         end do
      end if
      
   end subroutine host_gather2r8
   
   subroutine proc_gather2r8(a)
      ! Collect global arrays.
      real(kind=8), dimension(ifull), intent(in) :: a
      integer(kind=4),parameter :: ltype = MPI_DOUBLE_PRECISION
      integer(kind=4) :: ierr, lsize, lcomm
      real(kind=8), dimension(1,1) :: abuf

      if ( myid == 0 ) then
         write(6,*) "Error: ccmpi_gather argument required on proc 0"
         call ccmpi_abort(-1)
      end if

      lsize = ifull
      lcomm = comm_world
      call START_LOG(gathercc_begin)
      call MPI_Gather( a, lsize, ltype, abuf, lsize, ltype, 0_4, lcomm, ierr )
      call END_LOG(gathercc_end)
      
   end subroutine proc_gather2r8
   
   subroutine host_gather3(a,ag)
      ! Collect global arrays.
      real, dimension(:,:), intent(in) :: a
      real, dimension(:,:), intent(out) :: ag
      integer :: iproc
      integer :: ipoff, jpoff, npoff
      integer :: j, n, k, iq, iqg, kx
#ifdef i8r8
      integer(kind=4),parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4),parameter :: ltype = MPI_REAL
#endif
      integer(kind=4) :: ierr, lsize, lcomm
      real, dimension(ifull,size(a,2),0:nproc-1) :: abuf
      real, dimension(ifull,size(a,2)) :: atemp
      
      kx = size(a,2)
      lsize = ifull*kx
      lcomm = comm_world
      if ( size(a,1) == ifull ) then
         call START_LOG(gathercc_begin)
         call MPI_Gather( a, lsize, ltype, abuf, lsize, ltype, 0_4, lcomm, ierr )
         call END_LOG(gathercc_end)
      else    
         atemp(1:ifull,1:kx) = a(1:ifull,1:kx)
         call START_LOG(gathercc_begin)
         call MPI_Gather( atemp, lsize, ltype, abuf, lsize, ltype, 0_4, lcomm, ierr )
         call END_LOG(gathercc_end)
      end if   

      ! map array in order of processor rank
      if ( uniform_decomp ) then
         do iproc = 0,nproc-1
            call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,ipan,jpan)
            do k = 1,kx
               do n = 1,npan
                  do j = 1,jpan
                     ! Global indices are i+ipoff, j+jpoff, n-npoff
                     iqg = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g ! True global 1D index
                     iq = (j-1)*ipan + (n-1)*ipan*jpan
                     ag(iqg+1:iqg+ipan,k) = abuf(iq+1:iq+ipan,k,iproc)
                  end do
               end do
            end do
         end do
      else
         do iproc = 0,nproc-1
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan)
            do k = 1,kx
               do n = 1,npan
                  do j = 1,jpan
                     ! Global indices are i+ipoff, j+jpoff, n-npoff
                     iqg = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g ! True global 1D index
                     iq = (j-1)*ipan + (n-1)*ipan*jpan
                     ag(iqg+1:iqg+ipan,k) = abuf(iq+1:iq+ipan,k,iproc)
                  end do
               end do
            end do
         end do
      end if

   end subroutine host_gather3
   
   subroutine proc_gather3(a)
      ! Collect global arrays.
      real, dimension(:,:), intent(in) :: a
#ifdef i8r8
      integer(kind=4),parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4),parameter :: ltype = MPI_REAL
#endif
      integer(kind=4) :: ierr, lsize, lcomm
      real, dimension(1,1,1) :: abuf
      real, dimension(ifull,size(a,2)) :: atemp
      integer :: kx

      if ( myid == 0 ) then
         write(6,*) "Error: ccmpi_gather argument required on proc 0"
         call ccmpi_abort(-1)
      end if

      kx = size(a,2)
      lsize = ifull*kx
      lcomm = comm_world
      if ( size(a,1) == ifull ) then
         call START_LOG(gathercc_begin)
         call MPI_Gather( a, lsize, ltype, abuf, lsize, ltype, 0_4, lcomm, ierr )
         call END_LOG(gathercc_end)
      else    
         atemp(1:ifull,1:kx) = a(1:ifull,1:kx)
         call START_LOG(gathercc_begin)
         call MPI_Gather( atemp, lsize, ltype, abuf, lsize, ltype, 0_4, lcomm, ierr )
         call END_LOG(gathercc_end)
      end if   

   end subroutine proc_gather3

   subroutine host_gather3r8(a,ag)
      ! Collect global arrays.
      real(kind=8), dimension(:,:), intent(in) :: a
      real(kind=8), dimension(:,:), intent(out) :: ag
      integer :: iproc
      integer :: ipoff, jpoff, npoff
      integer :: j, n, k, iq, iqg, kx
      integer(kind=4),parameter :: ltype = MPI_DOUBLE_PRECISION
      integer(kind=4) :: ierr, lsize, lcomm
      real(kind=8), dimension(ifull,size(a,2),0:nproc-1) :: abuf
      real(kind=8), dimension(ifull,size(a,2)) :: atemp
      
      kx = size(a,2)
      lsize = ifull*kx
      lcomm = comm_world
      if ( size(a,1) == ifull ) then
         call START_LOG(gathercc_begin)
         call MPI_Gather( a, lsize, ltype, abuf, lsize, ltype, 0_4, lcomm, ierr )
         call END_LOG(gathercc_end)
      else    
         atemp(1:ifull,1:kx) = a(1:ifull,1:kx)
         call START_LOG(gathercc_begin)
         call MPI_Gather( atemp, lsize, ltype, abuf, lsize, ltype, 0_4, lcomm, ierr )
         call END_LOG(gathercc_end)
      end if  

      ! map array in order of processor rank
      if ( uniform_decomp ) then
         do iproc = 0,nproc-1
            call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,ipan,jpan)
            do k = 1,kx
               do n = 1,npan
                  do j = 1,jpan
                     ! Global indices are i+ipoff, j+jpoff, n-npoff
                     iqg = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g ! True global 1D index
                     iq = (j-1)*ipan + (n-1)*ipan*jpan
                     ag(iqg+1:iqg+ipan,k) = abuf(iq+1:iq+ipan,k,iproc)
                  end do
               end do
            end do
         end do
      else
         do iproc = 0,nproc-1
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan)
            do k = 1,kx
               do n = 1,npan
                  do j = 1,jpan
                     ! Global indices are i+ipoff, j+jpoff, n-npoff
                     iqg = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g ! True global 1D index
                     iq = (j-1)*ipan + (n-1)*ipan*jpan
                     ag(iqg+1:iqg+ipan,k) = abuf(iq+1:iq+ipan,k,iproc)
                  end do
               end do
            end do
         end do
      end if

   end subroutine host_gather3r8
   
   subroutine proc_gather3r8(a)
      ! Collect global arrays.
      real(kind=8), dimension(:,:), intent(in) :: a
      integer(kind=4),parameter :: ltype = MPI_DOUBLE_PRECISION
      integer(kind=4) :: ierr, lsize, lcomm
      real(kind=8), dimension(1,1,1) :: abuf
      real(kind=8), dimension(ifull,size(a,2)) :: atemp
      integer :: kx

      if ( myid == 0 ) then
         write(6,*) "Error: ccmpi_gather argument required on proc 0"
         call ccmpi_abort(-1)
      end if
      
      kx = size(a,2)
      lsize = ifull*kx
      lcomm = comm_world
      if ( size(a,1) == ifull ) then
         call START_LOG(gathercc_begin)
         call MPI_Gather( a, lsize, ltype, abuf, lsize, ltype, 0_4, lcomm, ierr ) 
         call END_LOG(gathercc_end)
      else    
         atemp(1:ifull,1:kx) = a(1:ifull,1:kx)
         call START_LOG(gathercc_begin)
         call MPI_Gather( atemp, lsize, ltype, abuf, lsize, ltype, 0_4, lcomm, ierr )
         call END_LOG(gathercc_end)
      end if   
      
   end subroutine proc_gather3r8
   
   subroutine host_gather4(a,ag)
      ! Collect global arrays.
      real, dimension(:,:,:), intent(in) :: a
      real, dimension(:,:,:), intent(out) :: ag
      integer :: iproc
      integer :: ipoff, jpoff, npoff
      integer :: j, n, k, l, iq, iqg, kx, lx
#ifdef i8r8
      integer(kind=4),parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4),parameter :: ltype = MPI_REAL
#endif
      integer(kind=4) :: ierr, lsize, lcomm
      real, dimension(ifull,size(a,2),size(a,3),0:nproc-1) :: abuf
      real, dimension(ifull,size(a,2),size(a,3)) :: atemp

      kx = size(a,2)
      lx = size(a,3)
      lsize = ifull*kx*lx
      lcomm = comm_world
      if ( size(a,1) == ifull ) then
         call START_LOG(gathercc_begin)
         call MPI_Gather( a, lsize, ltype, abuf, lsize, ltype, 0_4, lcomm, ierr )
         call END_LOG(gathercc_end)
      else    
         atemp(1:ifull,1:kx,1:lx) = a(1:ifull,1:kx,1:lx)
         call START_LOG(gathercc_begin)
         call MPI_Gather( atemp, lsize, ltype, abuf, lsize, ltype, 0_4, lcomm, ierr )
         call END_LOG(gathercc_end)
      end if   

      ! map array in order of processor rank
      if ( uniform_decomp ) then
         do iproc = 0,nproc-1
            call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,ipan,jpan)
            do l = 1,lx
               do k = 1,kx
                  do n = 1,npan
                     do j = 1,jpan
                        ! Global indices are i+ipoff, j+jpoff, n-npoff
                        iqg = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g ! True global 1D index
                        iq = (j-1)*ipan + (n-1)*ipan*jpan
                        ag(iqg+1:iqg+ipan,k,l) = abuf(iq+1:iq+ipan,k,l,iproc)
                     end do  
                  end do
               end do
            end do
         end do
      else
         do iproc = 0,nproc-1
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan)
            do l = 1,lx
               do k = 1,kx
                  do n = 1,npan
                     do j = 1,jpan
                        ! Global indices are i+ipoff, j+jpoff, n-npoff
                        iqg = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g ! True global 1D index
                        iq = (j-1)*ipan + (n-1)*ipan*jpan
                        ag(iqg+1:iqg+ipan,k,l) = abuf(iq+1:iq+ipan,k,l,iproc)
                     end do   
                  end do
               end do
            end do
         end do
      end if

   end subroutine host_gather4
   
   subroutine proc_gather4(a)
      ! Collect global arrays.
      real, dimension(:,:,:), intent(in) :: a
#ifdef i8r8
      integer(kind=4),parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4),parameter :: ltype = MPI_REAL
#endif
      integer(kind=4) :: ierr, lsize, lcomm
      real, dimension(1,1,1,1) :: abuf
      real, dimension(ifull,size(a,2),size(a,3)) :: atemp
      integer :: kx, lx

      if ( myid == 0 ) then
         write(6,*) "Error: ccmpi_gather argument required on proc 0"
         call ccmpi_abort(-1)
      end if

      kx = size(a,2)
      lx = size(a,3)
      lsize = ifull*kx*lx
      lcomm = comm_world
      if ( size(a,1) == ifull ) then
         call START_LOG(gathercc_begin)
         call MPI_Gather( a, lsize, ltype, abuf, lsize, ltype, 0_4, lcomm, ierr )
         call END_LOG(gathercc_end)
      else    
         atemp(1:ifull,1:kx,1:lx) = a(1:ifull,1:kx,1:lx)
         call START_LOG(gathercc_begin)
         call MPI_Gather( atemp, lsize, ltype, abuf, lsize, ltype, 0_4, lcomm, ierr )
         call END_LOG(gathercc_end)
      end if   

   end subroutine proc_gather4

   subroutine host_gather4r8(a,ag)
      ! Collect global arrays.
      real(kind=8), dimension(:,:,:), intent(in) :: a
      real(kind=8), dimension(:,:,:), intent(out) :: ag
      integer :: iproc
      integer :: ipoff, jpoff, npoff
      integer :: j, n, k, l, iq, iqg, kx, lx
      integer(kind=4),parameter :: ltype = MPI_DOUBLE_PRECISION
      integer(kind=4) :: ierr, lsize, lcomm
      real(kind=8), dimension(ifull,size(a,2),size(a,3),0:nproc-1) :: abuf
      real(kind=8), dimension(ifull,size(a,2),size(a,3)) :: atemp

      kx = size(a,2)
      lx = size(a,3)
      lsize = ifull*kx*lx
      lcomm = comm_world
      if ( size(a,1) == ifull ) then
         call START_LOG(gathercc_begin)
         call MPI_Gather( a, lsize, ltype, abuf, lsize, ltype, 0_4, lcomm, ierr )
         call END_LOG(gathercc_end)
      else    
         atemp(1:ifull,1:kx,1:lx) = a(1:ifull,1:kx,1:lx)
         call START_LOG(gathercc_begin)
         call MPI_Gather( atemp, lsize, ltype, abuf, lsize, ltype, 0_4, lcomm, ierr )
         call END_LOG(gathercc_end)
      end if   

      ! map array in order of processor rank
      if ( uniform_decomp ) then
         do iproc = 0,nproc-1
            call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,ipan,jpan)
            do l = 1,lx
               do k = 1,kx
                  do n = 1,npan
                     do j = 1,jpan
                        ! Global indices are i+ipoff, j+jpoff, n-npoff
                        iqg = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g ! True global 1D index
                        iq = (j-1)*ipan + (n-1)*ipan*jpan
                        ag(iqg+1:iqg+ipan,k,l) = abuf(iq+1:iq+ipan,k,l,iproc)
                     end do   
                  end do
               end do
            end do
         end do
      else
         do iproc = 0,nproc-1
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan)
            do l = 1,lx
               do k = 1,kx
                  do n = 1,npan
                     do j = 1,jpan
                        ! Global indices are i+ipoff, j+jpoff, n-npoff
                        iqg = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g ! True global 1D index
                        iq = (j-1)*ipan + (n-1)*ipan*jpan
                        ag(iqg+1:iqg+ipan,k,l) = abuf(iq+1:iq+ipan,k,l,iproc)
                     end do   
                  end do
               end do
            end do
         end do
      end if

   end subroutine host_gather4r8
   
   subroutine proc_gather4r8(a)
      ! Collect global arrays.
      real(kind=8), dimension(:,:,:), intent(in) :: a
      integer(kind=4),parameter :: ltype = MPI_DOUBLE_PRECISION
      integer(kind=4) :: ierr, lsize, lcomm
      real(kind=8), dimension(1,1,1,1) :: abuf
      real(kind=8), dimension(ifull,size(a,2),size(a,3)) :: atemp
      integer :: kx, lx

      if ( myid == 0 ) then
         write(6,*) "Error: ccmpi_gather argument required on proc 0"
         call ccmpi_abort(-1)
      end if
      
      kx = size(a,2)
      lx = size(a,3)
      lsize = ifull*kx*lx
      lcomm = comm_world
      if ( size(a,1) == ifull ) then
         call START_LOG(gathercc_begin)
         call MPI_Gather( a, lsize, ltype, abuf, lsize, ltype, 0_4, lcomm, ierr )
         call END_LOG(gathercc_end)
      else    
         atemp(1:ifull,1:kx,1:lx) = a(1:ifull,1:kx,1:lx)
         call START_LOG(gathercc_begin)
         call MPI_Gather( atemp, lsize, ltype, abuf, lsize, ltype, 0_4, lcomm, ierr )
         call END_LOG(gathercc_end)
      end if   

   end subroutine proc_gather4r8
   
   subroutine ccmpi_gatherall2(a,ag)
      ! Collect global arrays.
      real, dimension(ifull), intent(in) :: a
      real, dimension(ifull_g), intent(out) :: ag
      real, dimension(ifull,0:nproc-1) :: abuf
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      integer(kind=4) :: ierr, lsize, lcomm
      integer :: ipoff, jpoff, npoff
      integer :: j, n, iq, iqg, iproc

      lsize = ifull
      lcomm = comm_world
      call START_LOG(allgathercc_begin)
      call MPI_AllGather( a, lsize, ltype, abuf, lsize, ltype, lcomm, ierr )
      call END_LOG(allgathercc_end)

      ! map array in order of processor rank
      if ( uniform_decomp ) then
         do iproc = 0,nproc-1
            call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,ipan,jpan)
            do n = 1,npan
               do j = 1,jpan
                  ! Global indices are i+ipoff, j+jpoff, n-npoff
                  iqg = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g ! True global 1D index
                  iq = (j-1)*ipan + (n-1)*ipan*jpan
                  ag(iqg+1:iqg+ipan) = abuf(iq+1:iq+ipan,iproc)
               end do
            end do
         end do
      else
         do iproc = 0,nproc-1
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan)
            do n = 1,npan
               do j = 1,jpan
                  ! Global indices are i+ipoff, j+jpoff, n-npoff
                  iqg = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g ! True global 1D index
                  iq = (j-1)*ipan + (n-1)*ipan*jpan
                  ag(iqg+1:iqg+ipan) = abuf(iq+1:iq+ipan,iproc)
               end do
            end do
         end do
      end if

   end subroutine ccmpi_gatherall2
   
   subroutine ccmpi_gatherall3(a,ag)
      ! Collect global arrays.
      real, dimension(:,:), intent(in) :: a
      real, dimension(:,:), intent(out) :: ag
#ifdef i8r8
      integer(kind=4),parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4),parameter :: ltype = MPI_REAL
#endif
      integer(kind=4) :: ierr, lsize, lcomm
      real, dimension(ifull,size(a,2),0:nproc-1) :: abuf
      real, dimension(ifull,size(a,2)) :: atemp
      integer :: ipoff, jpoff, npoff
      integer :: j, n, k, iq, iqg, kx, iproc

      kx = size(a,2)
      lsize = ifull*kx
      lcomm = comm_world
      if ( size(a,1) == ifull ) then
         call START_LOG(allgathercc_begin)
         call MPI_AllGather( a, lsize, ltype, abuf, lsize, ltype, lcomm, ierr ) 
         call END_LOG(allgathercc_end)
      else    
         atemp(:,:) = a(1:ifull,1:kx)
         call START_LOG(allgathercc_begin)
         call MPI_AllGather( atemp, lsize, ltype, abuf, lsize, ltype, lcomm, ierr )
         call END_LOG(allgathercc_end)
      end if   

      ! map array in order of processor rank
      if ( uniform_decomp ) then
         do iproc = 0,nproc-1
            call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,ipan,jpan)
            do k = 1,kx
               do n = 1,npan
                  do j = 1,jpan
                     ! Global indices are i+ipoff, j+jpoff, n-npoff
                     iqg = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g ! True global 1D index
                     iq = (j-1)*ipan + (n-1)*ipan*jpan
                     ag(iqg+1:iqg+ipan,k) = abuf(iq+1:iq+ipan,k,iproc)
                  end do
               end do
            end do
         end do
      else
         do iproc = 0,nproc-1
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan)
            do k = 1,kx
               do n = 1,npan
                  do j = 1,jpan
                     ! Global indices are i+ipoff, j+jpoff, n-npoff
                     iqg = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g ! True global 1D index
                     iq = (j-1)*ipan + (n-1)*ipan*jpan
                     ag(iqg+1:iqg+ipan,k) = abuf(iq+1:iq+ipan,k,iproc)
                  end do
               end do
            end do
         end do
      end if

   end subroutine ccmpi_gatherall3

   subroutine ccmpi_gathermap_send2(a)

      real, dimension(ifull), intent(in) :: a
      integer :: ncount, w
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      integer(kind=4) :: ierr, lsize
      integer(kind=4) :: lcomm
      integer(kind=4) :: itag = 52
        
      ncount = size(specmap_recv)
      lsize = ifull
      lcomm = comm_world

      ! Set up the buffers to recv
      nreq = 0
      do w = 1,ncount
         nreq = nreq + 1
         rlist(nreq) = w
         call MPI_IRecv( bnds(specmap_recv(w))%rbuf, lsize, ltype, specmap_recv(w), itag, lcomm, ireq(nreq), ierr )
      end do
      rreq = nreq

      ! Set up the buffers to send
      bnds(myid)%sbuf(1:ifull) = a
      do w = 1,size(specmap_send)
         nreq = nreq + 1
         call MPI_ISend( bnds(myid)%sbuf, lsize, ltype, specmap_send(w), itag, lcomm, ireq(nreq), ierr )
      end do

   end subroutine ccmpi_gathermap_send2

   subroutine ccmpi_gathermap_recv2(kref)

      integer, intent(in) :: kref
      integer :: ncount, w, iproc, n, iq
      integer :: ipoff, jpoff, npoff
      integer :: ipak, jpak
      integer :: sreq, rcount, jproc
      integer(kind=4) :: ierr, lsize
      integer(kind=4) :: lcomm, ldone
      integer(kind=4), dimension(size(specmap_recv)) :: donelist
      integer(kind=4), dimension(MPI_STATUS_SIZE,size(ireq)) :: status
      
      ncount = size(specmap_recv)
      lsize = ifull
      lcomm = comm_world

      ! Unpack incomming messages
      rcount = rreq
      if ( uniform_decomp ) then
         do while ( rcount > 0 )
            call START_LOG(mpiwaitmap_begin) 
            call MPI_Waitsome( rreq, ireq, ldone, donelist, status, ierr )
            call END_LOG(mpiwaitmap_end)
            rcount = rcount - ldone
            do jproc = 1,ldone
               w = rlist(donelist(jproc))
               iproc = specmap_recv(w)
               call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,ipan,jpan)
               ipak = ipoff/ipan
               jpak = jpoff/jpan
               do n = 1,npan
                  ! Global indices are i+ipoff, j+jpoff, n-npoff
                  iq = (n-1)*ipan*jpan
                  globalpack(ipak,jpak,n-npoff)%localdata(:,:,kref+1) = &
                     reshape( bnds(iproc)%rbuf(iq+1:iq+ipan*jpan), (/ ipan, jpan /) )
               end do
            end do
         end do   
      else
         do while ( rcount > 0 )
            call START_LOG(mpiwaitmap_begin) 
            call MPI_Waitsome( rreq, ireq, ldone, donelist, status, ierr )
            call END_LOG(mpiwaitmap_end)
            rcount = rcount - ldone
            do jproc = 1,ldone
               w = rlist(donelist(jproc))
               iproc = specmap_recv(w)
               call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan)
               ipak = ipoff/ipan
               jpak = jpoff/jpan
               do n = 1,npan
                  ! Global indices are i+ipoff, j+jpoff, n-npoff
                  iq = (n-1)*ipan*jpan
                  globalpack(ipak,jpak,n-npoff)%localdata(:,:,kref+1) = &
                     reshape( bnds(iproc)%rbuf(iq+1:iq+ipan*jpan), (/ ipan, jpan /) )
               end do
            end do
         end do
      end if
      
      sreq = nreq - rreq
      if ( sreq > 0 ) then
         call START_LOG(mpiwaitmap_begin) 
         call MPI_Waitall( sreq, ireq(rreq+1:nreq), status, ierr )
         call END_LOG(mpiwaitmap_end)
      end if   

   end subroutine ccmpi_gathermap_recv2

    subroutine ccmpi_gathermap_send3(a)

      real, dimension(:,:), intent(in) :: a
      integer :: ncount, w, kx
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      integer(kind=4) :: ierr, lsize
      integer(kind=4) :: lcomm
      integer(kind=4) :: itag = 52

      kx = size(a,2)
      ncount = size(specmap_recv)
      lsize = ifull*kx
      lcomm = comm_world

      bnds(myid)%sbuf(1:ifull*kx) = reshape( a(1:ifull,1:kx), (/ ifull*kx /) )      
      
      ! Set up the buffers to recv
      nreq = 0
      do w = 1,ncount
         nreq = nreq + 1
         rlist(nreq) = w
         call MPI_IRecv( bnds(specmap_recv(w))%rbuf, lsize, ltype, specmap_recv(w), itag, lcomm, ireq(nreq), ierr )
      end do
      rreq = nreq
      
      ! Set up the buffers to send
      do w = 1,size(specmap_send)
         nreq = nreq + 1
         call MPI_ISend( bnds(myid)%sbuf, lsize, ltype, specmap_send(w), itag, lcomm, ireq(nreq), ierr )
      end do

    end subroutine ccmpi_gathermap_send3
  
   subroutine ccmpi_gathermap_recv3(kx,kref)

      integer, intent(in) :: kx, kref
      integer :: ncount, w, iproc, k, n, iq
      integer :: ipoff, jpoff, npoff
      integer :: ipak, jpak
      integer :: sreq, rcount, jproc
      integer(kind=4) :: ierr, lsize
      integer(kind=4) :: lcomm, ldone
      integer(kind=4), dimension(size(specmap_recv)) :: donelist
      integer(kind=4), dimension(MPI_STATUS_SIZE,size(ireq)) :: status

      ncount = size(specmap_recv)
      lsize = ifull*kx
      lcomm = comm_world

      ! Unpack incomming messages
      rcount = rreq
      if ( uniform_decomp ) then
         do while ( rcount > 0 )
            call START_LOG(mpiwaitmap_begin) 
            call MPI_Waitsome( rreq, ireq, ldone, donelist, status, ierr )
            call END_LOG(mpiwaitmap_end)
            rcount = rcount - ldone
            do jproc = 1,ldone
               w = rlist(donelist(jproc))        
               iproc = specmap_recv(w)
               call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,ipan,jpan)
               ipak = ipoff/ipan
               jpak = jpoff/jpan
               do k = 1,kx
                  do n = 1,npan
                     ! Global indices are i+ipoff, j+jpoff, n-npoff
                     iq = (n-1)*ipan*jpan + (k-1)*ifull
                     globalpack(ipak,jpak,n-npoff)%localdata(:,:,kref+k) = &
                        reshape( bnds(iproc)%rbuf(iq+1:iq+ipan*jpan), (/ ipan, jpan /) )
                  end do
               end do
            end do
         end do   
      else
         do while ( rcount > 0 )
            call START_LOG(mpiwaitmap_begin) 
            call MPI_Waitsome( rreq, ireq, ldone, donelist, status, ierr )
            call END_LOG(mpiwaitmap_end)
            rcount = rcount - ldone
            do jproc = 1,ldone
               w = rlist(donelist(jproc))
               iproc = specmap_recv(w)
               call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan)
               ipak = ipoff/ipan
               jpak = jpoff/jpan
               do k = 1,kx
                  do n = 1,npan
                     ! Global indices are i+ipoff, j+jpoff, n-npoff
                     iq = (n-1)*ipan*jpan + (k-1)*ifull
                     globalpack(ipak,jpak,n-npoff)%localdata(:,:,kref+k) = &
                       reshape( bnds(iproc)%rbuf(iq+1:iq+ipan*jpan), (/ ipan, jpan /) )
                  end do
               end do
            end do
         end do
      end if
      
      sreq = nreq - rreq
      if ( sreq > 0 ) then
         call START_LOG(mpiwaitmap_begin) 
         call MPI_Waitall( sreq, ireq(rreq+1:nreq), status, ierr )
         call END_LOG(mpiwaitmap_end)
      end if   

   end subroutine ccmpi_gathermap_recv3
    
   subroutine setglobalpack_v(datain,ibeg,iend,k)
   
      ! This subroutine assigns a value to a gridpoint
      ! in the global sparse array
   
      integer, intent(in) :: ibeg, iend, k
      integer :: il2, iqg, im1, jm1, ilen, jlen, iq
      integer :: b_n, b_ipak, b_jpak, b_iloc, b_jloc
      integer :: e_n, e_ipak, e_jpak, e_iloc, e_jloc
      integer :: s_ipak, s_jpak, s_iloc, s_jloc
      integer :: c_ipak, c_jpak
      real, dimension(:), intent(in) :: datain
      
      il2 = il_g*il_g
      
      iqg = ibeg - 1
      b_n = iqg/il2
      iqg = iqg - b_n*il2
      jm1 = iqg/il_g
      im1 = iqg - jm1*il_g
      b_ipak = im1/ipan
      b_jpak = jm1/jpan
      b_iloc = im1 + 1 - b_ipak*ipan
      b_jloc = jm1 + 1 - b_jpak*jpan

      iqg = iend - 1
      e_n = iqg/il2
      iqg = iqg - e_n*il2
      jm1 = iqg/il_g
      im1 = iqg - jm1*il_g
      e_ipak = im1/ipan
      e_jpak = jm1/jpan
      e_iloc = im1 + 1 - e_ipak*ipan
      e_jloc = jm1 + 1 - e_jpak*jpan
      
      if ( e_jpak >= b_jpak) then
         s_jpak = 1
      else
         s_jpak = -1
      end if
      if ( e_ipak >= b_ipak) then
         s_ipak = 1
      else
         s_ipak = -1
      end if
      if ( e_jloc >= b_jloc) then
         s_jloc = 1
      else
         s_jloc = -1
      end if
      if ( e_iloc >= b_iloc) then
         s_iloc = 1
      else
         s_iloc = -1
      end if

      ilen = abs(e_iloc-b_iloc) + 1
      jlen = abs(e_jloc-b_jloc) + 1
            
      if ( jlen > ilen ) then
         iq = 0
         do c_jpak = b_jpak,e_jpak,s_jpak
            globalpack(b_ipak,c_jpak,b_n)%localdata(b_iloc,b_jloc:e_jloc:s_jloc,k) = datain(iq+1:iq+jlen)
            iq = iq + jlen
         end do
      else
         iq = 0
         do c_ipak = b_ipak,e_ipak,s_ipak
            globalpack(c_ipak,b_jpak,b_n)%localdata(b_iloc:e_iloc:s_iloc,b_jloc,k) = datain(iq+1:iq+ilen)
            iq = iq + ilen
         end do
      end if
   
   end subroutine setglobalpack_v
   
   subroutine getglobalpack_v(dataout,ibeg,iend,k)
   
      ! This subroutine returns a value from a gridpoint
      ! in the global sparse array

      integer, intent(in) :: ibeg, iend, k
      integer :: il2, iqg, im1, jm1, ilen, jlen, iq
      integer :: b_n, b_ipak, b_jpak, b_iloc, b_jloc
      integer :: e_n, e_ipak, e_jpak, e_iloc, e_jloc
      integer :: s_ipak, s_jpak, s_iloc, s_jloc
      integer :: c_ipak, c_jpak
      real, dimension(:), intent(out) :: dataout
      
      il2 = il_g*il_g
      
      iqg = ibeg - 1
      b_n = iqg/il2
      iqg = iqg - b_n*il2
      jm1 = iqg/il_g
      im1 = iqg - jm1*il_g
      b_ipak = im1/ipan
      b_jpak = jm1/jpan
      b_iloc = im1 + 1 - b_ipak*ipan
      b_jloc = jm1 + 1 - b_jpak*jpan

      iqg = iend - 1
      e_n = iqg/il2
      iqg = iqg - e_n*il2
      jm1 = iqg/il_g
      im1 = iqg - jm1*il_g
      e_ipak = im1/ipan
      e_jpak = jm1/jpan
      e_iloc = im1 + 1 - e_ipak*ipan
      e_jloc = jm1 + 1 - e_jpak*jpan
            
      if ( e_jpak >= b_jpak) then
         s_jpak = 1
      else
         s_jpak = -1
      end if
      if ( e_ipak >= b_ipak) then
         s_ipak = 1
      else
         s_ipak = -1
      end if
      if ( e_jloc >= b_jloc) then
         s_jloc = 1
      else
         s_jloc = -1
      end if
      if ( e_iloc >= b_iloc) then
         s_iloc = 1
      else
         s_iloc = -1
      end if

      ilen = abs(e_iloc-b_iloc) + 1
      jlen = abs(e_jloc-b_jloc) + 1
      
      if ( jlen > ilen ) then
         iq = 0
         do c_jpak = b_jpak,e_jpak,s_jpak
            dataout(iq+1:iq+jlen) = globalpack(b_ipak,c_jpak,b_n)%localdata(b_iloc,b_jloc:e_jloc:s_jloc,k)
            iq = iq + jlen
         end do
      else
         iq = 0
         do c_ipak = b_ipak,e_ipak,s_ipak
            dataout(iq+1:iq+ilen) = globalpack(c_ipak,b_jpak,b_n)%localdata(b_iloc:e_iloc:s_iloc,b_jloc,k)
            iq = iq + ilen
         end do
      end if
      
   end subroutine getglobalpack_v
   
   subroutine copyglobalpack(krefin,krefout,kx)

      ! This routine copies one section of the global sparse array
      ! to another section.  Note that it only copies the memory
      ! assigned by gathermap.  specmap needs to be replaced with
      ! spectmapext to copy all parts of the global sparse array.
   
      integer, intent(in) :: krefin, krefout, kx
      integer :: w, n, ncount, iproc, ipak, jpak
      integer :: ipoff, jpoff, npoff
   
      ncount = size(specmap_recv)
      if ( uniform_decomp ) then
         do w = 1,ncount
            iproc = specmap_recv(w)
            call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,ipan,jpan)
            ipak = ipoff/ipan
            jpak = jpoff/jpan
            do n = 1,npan
               ! Global indices are i+ipoff, j+jpoff, n-npoff
               globalpack(ipak,jpak,n-npoff)%localdata(:,:,krefout+1:krefout+kx) = &
                  globalpack(ipak,jpak,n-npoff)%localdata(:,:,krefin+1:krefin+kx)
            end do
         end do
      else
         do w = 1,ncount
            iproc = specmap_recv(w)
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan)
            ipak = ipoff/ipan
            jpak = jpoff/jpan
            do n = 1,npan
               ! Global indices are i+ipoff, j+jpoff, n-npoff
               globalpack(ipak,jpak,n-npoff)%localdata(:,:,krefout+1:krefout+kx) = &
                  globalpack(ipak,jpak,n-npoff)%localdata(:,:,krefin+1:krefin+kx)
            end do
         end do
      end if
   
   end subroutine copyglobalpack

   subroutine allocateglobalpack(kx,ky)
   
      ! This allocates global sparse arrays for the digital filter.
      ! Usually this is 1:kl or 1:ol in size, but for some configurations
      ! we need to store the original values and hence use 1:2*kl or 1:2*ol.
      ! Also, the 0 index is to store the sum term for the digital filter.
   
      integer, intent(in) :: kx, ky
      integer :: ncount, w, ipak, jpak, n, iproc
      integer :: ipoff, jpoff, npoff
      integer :: xlen
      
      xlen = size(specmap_recv) + size(specmap_send)
      if ( size(ireq)<xlen ) then
         deallocate( ireq )
         allocate( ireq(xlen) )
      end if
      xlen = size(specmap_recv)
      if ( size(rlist)<xlen ) then
         deallocate( rlist )
         allocate( rlist(xlen) )
      end if
      xlen = ifull*ky
      do w = 1,size(specmap_recv)
         iproc = specmap_recv(w) 
         if ( bnds(iproc)%rbuflen < xlen ) then
            if ( allocated(bnds(iproc)%rbuf) ) then
               deallocate( bnds(iproc)%rbuf )
               deallocate( bnds(iproc)%r8buf )
            end if
            allocate( bnds(iproc)%rbuf(xlen) )
            allocate( bnds(iproc)%r8buf(xlen) )
            bnds(iproc)%rbuflen = xlen
         end if
      end do   
      iproc = myid
      if ( bnds(iproc)%sbuflen < xlen ) then
         if ( allocated(bnds(iproc)%sbuf) ) then
            deallocate( bnds(iproc)%sbuf )
            deallocate( bnds(iproc)%s8buf )
         end if
         allocate( bnds(iproc)%sbuf(xlen) )
         allocate( bnds(iproc)%s8buf(xlen) )
         bnds(iproc)%sbuflen = xlen
      end if 
      
      ! allocate globalpack arrays for 1D scale-selective filter
      allocate(globalpack(0:nxproc-1,0:nyproc-1,0:5))
      ncount = size(specmap_ext)
      if ( uniform_decomp ) then
         do w = 1,ncount
            iproc = specmap_ext(w)
            call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,ipan,jpan)
            ! Global indices are i+ipoff, j+jpoff, n-npoff
            ipak = ipoff/ipan
            jpak = jpoff/jpan
            do n = 1,npan
               allocate(globalpack(ipak,jpak,n-npoff)%localdata(ipan,jpan,0:kx))
               globalpack(ipak,jpak,n-npoff)%localdata = 0.
            end do
         end do
      else
         do w = 1,ncount
            iproc = specmap_ext(w)
            call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan)
            ! Global indices are i+ipoff, j+jpoff, n-npoff
            ipak = ipoff/ipan
            jpak = jpoff/jpan
            do n = 1,npan
               allocate(globalpack(ipak,jpak,n-npoff)%localdata(ipan,jpan,0:kx))
               globalpack(ipak,jpak,n-npoff)%localdata = 0.
            end do
         end do
      end if
      
      deallocate(specmap_ext) ! not needed after allocation of global sparse arrays
   
   end subroutine allocateglobalpack
   
   subroutine ccglobal_posneg2(array, delpos, delneg)
      ! Calculate global sums of positive and negative values of array
      use sumdd_m
      use xyzinfo_m       

      real, intent(in), dimension(ifull) :: array
      real, intent(out) :: delpos, delneg
      integer(kind=4) :: ierr, lcomm
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_COMPLEX
#else
      integer(kind=4), parameter :: ltype = MPI_COMPLEX
#endif   
      complex, dimension(2) :: local_sum, global_sum
      real, dimension(ifull) :: tmparr 

      local_sum(1:2) = cmplx(0., 0.)
      tmparr(1:ifull)  = max(0., array(1:ifull)*wts(1:ifull))
      call drpdr_local(tmparr, local_sum(1))
      tmparr(1:ifull)  = min(0., array(1:ifull)*wts(1:ifull))
      call drpdr_local(tmparr, local_sum(2))
      global_sum(1:2) = cmplx(0., 0.)
      lcomm = comm_world
      call START_LOG(allreducepn_begin)
      call MPI_Allreduce( local_sum, global_sum, 2_4, ltype, MPI_SUMDR, lcomm, ierr )
      call END_LOG(allreducepn_end)
      delpos = real(global_sum(1))
      delneg = real(global_sum(2))

   end subroutine ccglobal_posneg2
    
   subroutine ccglobal_posneg3(array, delpos, delneg, dsig)
      ! Calculate global sums of positive and negative values of array
      use sumdd_m
      use xyzinfo_m
      real, intent(in), dimension(:,:) :: array
      real, intent(in), dimension(:) :: dsig
      real, intent(out) :: delpos, delneg
      integer :: k, kx
      integer(kind=4) :: ierr, lcomm
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_COMPLEX
#else
      integer(kind=4), parameter :: ltype = MPI_COMPLEX
#endif   
      complex, dimension(2) :: local_sum, global_sum
      real, dimension(ifull) :: tmparr
       
      kx = size(array,2)   
      local_sum(1:2) = cmplx(0., 0.)
      do k = 1,kx
         tmparr(1:ifull) = max(0., abs(dsig(k))*array(1:ifull,k)*wts(1:ifull))
         call drpdr_local(tmparr, local_sum(1))
         tmparr(1:ifull) = min(0., abs(dsig(k))*array(1:ifull,k)*wts(1:ifull))
         call drpdr_local(tmparr, local_sum(2))
      end do ! k loop
      global_sum(1:2) = cmplx(0., 0.)
      lcomm = comm_world
      call START_LOG(allreducepn_begin)
      call MPI_Allreduce( local_sum, global_sum, 2_4, ltype, MPI_SUMDR, lcomm, ierr )
      call END_LOG(allreducepn_end)
      delpos = real(global_sum(1))
      delneg = real(global_sum(2))

   end subroutine ccglobal_posneg3
   
   subroutine ccglobal_posneg3o(array, delpos, delneg, dsig)
      ! Calculate global sums of positive and negative values of array
      use sumdd_m
      use xyzinfo_m
      real, intent(in), dimension(:,:) :: array
      real, intent(in), dimension(:,:) :: dsig
      real, intent(out) :: delpos, delneg
      integer :: k, kx
      integer(kind=4) :: ierr, lcomm
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_COMPLEX
#else
      integer(kind=4), parameter :: ltype = MPI_COMPLEX
#endif   
      complex, dimension(2) :: local_sum, global_sum
      real, dimension(ifull) :: tmparr
       
      kx = size(array,2)   
      local_sum(1:2) = cmplx(0., 0.)
      do k = 1,kx
         tmparr(1:ifull) = max(0., abs(dsig(1:ifull,k))*array(1:ifull,k)*wts(1:ifull))
         call drpdr_local(tmparr, local_sum(1))
         tmparr(1:ifull) = min(0., abs(dsig(1:ifull,k))*array(1:ifull,k)*wts(1:ifull))
         call drpdr_local(tmparr, local_sum(2))
      end do ! k loop
      global_sum(1:2) = cmplx(0., 0.)
      lcomm = comm_world
      call START_LOG(allreducepn_begin)
      call MPI_Allreduce( local_sum, global_sum, 2_4, ltype, MPI_SUMDR, lcomm, ierr )
      call END_LOG(allreducepn_end)
      delpos = real(global_sum(1))
      delneg = real(global_sum(2))

   end subroutine ccglobal_posneg3o

   subroutine ccglobal_posneg4 (array, delpos, delneg, dsig)
      ! Calculate global sums of positive and negative values of array
      use sumdd_m
      use xyzinfo_m
      real, intent(in), dimension(:,:,:) :: array
      real, intent(in), dimension(:) :: dsig
      real, intent(out), dimension(:) :: delpos, delneg
      integer :: i, k, kx, ntr
      integer(kind=4) :: ierr, mnum, lcomm
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_COMPLEX
#else
      integer(kind=4), parameter :: ltype = MPI_COMPLEX
#endif   
      complex, dimension(2*size(array,3)) :: local_sum, global_sum
      complex, dimension(2*size(array,2)) :: local_sum_t
      real, dimension(ifull,2*size(array,2)) :: tmparr

      kx  = size(array,2)
      ntr = size(array,3)
      do i = 1,ntr
         do k = 1,kx
            tmparr(1:ifull,k) = max(0.,abs(dsig(k))*array(1:ifull,k,i)*wts(1:ifull))
            tmparr(1:ifull,k+kx) = min(0.,abs(dsig(k))*array(1:ifull,k,i)*wts(1:ifull))
         end do ! k loop
         local_sum_t(1:2*kx) = cmplx(0.,0.)
         call drpdr_local_v(tmparr, local_sum_t)
         local_sum(i) = cmplx(0.,0.)
         local_sum(i+ntr) = cmplx(0.,0.)
         do k = 1,kx
            call drpdr(local_sum_t(k:k), local_sum(i:i), 1, ltype) 
            call drpdr(local_sum_t(k+kx:k+kx), local_sum(i+ntr:i+ntr), 1, ltype) 
         end do
      end do
      mnum = 2*ntr
      global_sum(1:2*ntr) = cmplx(0.,0.)
      lcomm = comm_world
      call START_LOG(allreducepn_begin)
      call MPI_Allreduce( local_sum, global_sum, mnum, ltype, MPI_SUMDR, lcomm, ierr )
      call END_LOG(allreducepn_end)
      delpos(1:ntr) = real(global_sum(1:ntr))
      delneg(1:ntr) = real(global_sum(ntr+1:2*ntr))

   end subroutine ccglobal_posneg4
   
   subroutine ccglobal_posneg4o(array, delpos, delneg, dsig)
      ! Calculate global sums of positive and negative values of array
      use sumdd_m
      use xyzinfo_m
      real, intent(in), dimension(:,:,:) :: array
      real, intent(in), dimension(:,:) :: dsig
      real, intent(out), dimension(:) :: delpos, delneg
      integer :: i, k, kx, ntr
      integer(kind=4) :: ierr, mnum, lcomm
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_COMPLEX
#else
      integer(kind=4), parameter :: ltype = MPI_COMPLEX
#endif   
      complex, dimension(2*size(array,3)) :: local_sum, global_sum
      real, dimension(ifull) :: tmparr

      kx  = size(array,2)
      ntr = size(array,3)
      local_sum(1:2*ntr) = cmplx(0.,0.)
      do i = 1,ntr
         do k = 1,kx
            tmparr(1:ifull) = max(0.,abs(dsig(1:ifull,k))*array(1:ifull,k,i)*wts(1:ifull))
            call drpdr_local(tmparr, local_sum(i))
            tmparr(1:ifull) = min(0.,abs(dsig(1:ifull,k))*array(1:ifull,k,i)*wts(1:ifull))
            call drpdr_local(tmparr, local_sum(i+ntr))
         end do ! k loop
      end do
      mnum = 2*ntr
      global_sum(1:2*ntr) = cmplx(0.,0.)
      lcomm = comm_world
      call START_LOG(allreducepn_begin)
      call MPI_Allreduce( local_sum, global_sum, mnum, ltype, MPI_SUMDR, lcomm, ierr )
      call END_LOG(allreducepn_end)
      delpos(1:ntr) = real(global_sum(1:ntr))
      delneg(1:ntr) = real(global_sum(ntr+1:2*ntr))

   end subroutine ccglobal_posneg4o
   
   subroutine bounds_setup(dt)

      use const_phys, only : rearth
      use indices_m
      
      integer :: n, i, j, iq, iqq, mycol, ncount
      integer :: iproc, rproc, sproc
      integer :: iqg, iql, iloc, jloc, nloc, icol
      integer :: iext, iextu, iextv
      integer, dimension(:), allocatable :: dumi
      integer, dimension(:,:), allocatable :: dums, dumr
      integer(kind=4) :: ierr, itag=0
      integer(kind=4) :: llen, lproc, lcomm
      integer(kind=4), dimension(:,:), allocatable :: status
#ifdef i8r8
      integer(kind=4), parameter :: ltypei = MPI_INTEGER8
      integer(kind=4), parameter :: ltyper = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltypei = MPI_INTEGER
      integer(kind=4), parameter :: ltyper = MPI_REAL
#endif
      logical :: swap
      logical, dimension(:), allocatable :: neigharray_g
      logical(kind=4), dimension(:,:), allocatable :: dumsl, dumrl
      real :: maxdis
      real, intent(in) :: dt

      ! Just set values that point to values within own processors region.
      ! Other values are set up later
      ! Set initial values to make sure missing values are caught properly.
      in = huge(1)
      is = huge(1)
      iw = huge(1)
      ie = huge(1)
      inn = huge(1)
      iss = huge(1)
      iww = huge(1)
      iee = huge(1)
      ien = huge(1)
      ine = huge(1)
      ise = huge(1)
      iwn = huge(1)
      inw = huge(1)
      isw = huge(1)
      ies = huge(1)
      iws = huge(1)
      ieu = huge(1)
      iwu = huge(1)
      inv = huge(1)
      isv = huge(1)
      iev = huge(1)
      iwv = huge(1)
      inu = huge(1)
      isu = huge(1)
      ieeu = huge(1)
      iwwu = huge(1)
      innv = huge(1)
      issv = huge(1)
      innu = huge(1)
      ieev = huge(1)
      lwws = huge(1)
      lwss = huge(1)
      lees = huge(1)
      less = huge(1)
      lwwn = huge(1)
      lwnn = huge(1)
      leen = huge(1)
      lenn = huge(1)
      lsww = huge(1)
      lssw = huge(1)
      lsee = huge(1)
      lsse = huge(1)
      lnww = huge(1)
      lnnw = huge(1)
      lnee = huge(1)
      lnne = huge(1)
      do n=1,npan
         do j=1,jpan
            do i=1,ipan
               iq = indp(i,j,n)   ! Local
               iqg = indg(i,j,n)  ! Global

               iqq = in_g(iqg)    ! Global neighbour index
               rproc = qproc(iqq) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  call indv_mpi(iqq,iloc,jloc,nloc)
                  in(iq) = indp(iloc,jloc,nloc)
               end if
               iqq = inn_g(iqg)   ! Global neighbour index
               rproc = qproc(iqq) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  call indv_mpi(iqq,iloc,jloc,nloc)
                  inn(iq) = indp(iloc,jloc,nloc)
               end if

               iqq = is_g(iqg)    ! Global neighbour index
               rproc = qproc(iqq) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  call indv_mpi(iqq,iloc,jloc,nloc)
                  is(iq) = indp(iloc,jloc,nloc)
               end if
               iqq = iss_g(iqg)    ! Global neighbour index
               rproc = qproc(iqq) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  call indv_mpi(iqq,iloc,jloc,nloc)
                  iss(iq) = indp(iloc,jloc,nloc)
               end if

               iqq = ie_g(iqg)    ! Global neighbour index
               rproc = qproc(iqq) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  call indv_mpi(iqq,iloc,jloc,nloc)
                  ie(iq) = indp(iloc,jloc,nloc)
               end if
               iqq = iee_g(iqg)    ! Global neighbour index
               rproc = qproc(iqq) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  call indv_mpi(iqq,iloc,jloc,nloc)
                  iee(iq) = indp(iloc,jloc,nloc)
               end if

               iqq = iw_g(iqg)    ! Global neighbour index
               rproc = qproc(iqq) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  call indv_mpi(iqq,iloc,jloc,nloc)
                  iw(iq) = indp(iloc,jloc,nloc)
               end if
               iqq = iww_g(iqg)    ! Global neighbour index
               rproc = qproc(iqq) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  call indv_mpi(iqq,iloc,jloc,nloc)
                  iww(iq) = indp(iloc,jloc,nloc)
               end if

               ! Note that the model only needs a limited set of the diagonal
               ! index arrays
               iqq = ine_g(iqg)    ! Global neighbour index
               rproc = qproc(iqq) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  call indv_mpi(iqq,iloc,jloc,nloc)
                  ine(iq) = indp(iloc,jloc,nloc)
               end if

               iqq = ise_g(iqg)    ! Global neighbour index
               rproc = qproc(iqq) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  call indv_mpi(iqq,iloc,jloc,nloc)
                  ise(iq) = indp(iloc,jloc,nloc)
               end if

               iqq = ien_g(iqg)    ! Global neighbour index
               rproc = qproc(iqq) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  call indv_mpi(iqq,iloc,jloc,nloc)
                  ien(iq) = indp(iloc,jloc,nloc)
               end if

               iqq = iwn_g(iqg)    ! Global neighbour index
               rproc = qproc(iqq) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  call indv_mpi(iqq,iloc,jloc,nloc)
                  iwn(iq) = indp(iloc,jloc,nloc)
               end if

               iqq = inw_g(iqg)    ! Global neighbour index
               rproc = qproc(iqq) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  call indv_mpi(iqq,iloc,jloc,nloc)
                  inw(iq) = indp(iloc,jloc,nloc)
               end if

               iqq = isw_g(iqg)    ! Global neighbour index
               rproc = qproc(iqq) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  call indv_mpi(iqq,iloc,jloc,nloc)
                  isw(iq) = indp(iloc,jloc,nloc)
               end if

               iqq = ies_g(iqg)    ! Global neighbour index
               rproc = qproc(iqq) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  call indv_mpi(iqq,iloc,jloc,nloc)
                  ies(iq) = indp(iloc,jloc,nloc)
               end if

               iqq = iws_g(iqg)    ! Global neighbour index
               rproc = qproc(iqq) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  call indv_mpi(iqq,iloc,jloc,nloc)
                  iws(iq) = indp(iloc,jloc,nloc)
               end if

            end do
         end do
      end do

      ! Correct within the same face only (not necessarily the same
      ! processor, but will be corrected later).
      ieu = ie
      iwu = iw
      inv = in
      isv = is
      iev = ie
      iwv = iw
      inu = in
      isu = is

      ! Initialise the edge variables
      edge_w = ioff == 0
      edge_s = joff == 0
      edge_n = joff == il_g - jpan
      edge_e = ioff == il_g - ipan

      do n = 0,nproc-1
         bnds(n)%len = 0
         bnds(n)%rbuflen = 0
         bnds(n)%sbuflen = 0
         bnds(n)%rlenh_bg(1:maxcolour) = 0
         bnds(n)%rlenh_fn(1:maxcolour) = 0
         bnds(n)%slenh_bg(1:maxcolour) = 0
         bnds(n)%slenh_fn(1:maxcolour) = 0
         bnds(n)%rlen_bg(1:maxcolour) = 0
         bnds(n)%rlen_fn(1:maxcolour) = 0
         bnds(n)%slen_bg(1:maxcolour) = 0
         bnds(n)%slen_fn(1:maxcolour) = 0
         bnds(n)%rlenx_bg(1:maxcolour) = 0
         bnds(n)%rlenx_fn(1:maxcolour) = 0
         bnds(n)%slenx_bg(1:maxcolour) = 0
         bnds(n)%slenx_fn(1:maxcolour) = 0
         bnds(n)%rlen2 = 0
         bnds(n)%slen2 = 0
         bnds(n)%rlen_su_bg = 0
         bnds(n)%rlen_ev_fn = 0
         bnds(n)%rlen_nnu_bg = 0
         bnds(n)%rlen_eev_fn = 0
         bnds(n)%rlen_sv_bg = 0
         bnds(n)%rlen_wu_fn = 0
         bnds(n)%rlen_nv_bg = 0
         bnds(n)%rlen_eu_fn = 0
         bnds(n)%rlen_ssv_bg = 0
         bnds(n)%rlen_wwu_fn = 0
         bnds(n)%rlen_nnv_bg = 0
         bnds(n)%rlen_eu_fn = 0
         bnds(n)%slen_su_bg = 0
         bnds(n)%slen_ev_fn = 0
         bnds(n)%slen_nnu_bg = 0
         bnds(n)%slen_eev_fn = 0
         bnds(n)%slen_sv_bg = 0
         bnds(n)%slen_wu_fn = 0
         bnds(n)%slen_nv_bg = 0
         bnds(n)%slen_eu_fn = 0
         bnds(n)%slen_ssv_bg = 0
         bnds(n)%slen_wwu_fn = 0
         bnds(n)%slen_nnv_bg = 0
         bnds(n)%slen_eeu_fn = 0     
      end do   

!     In the first pass through, set up list of points to be requested from
!     other processors. These points are placed in the "iextra" region at the
!     end of arrays. The indirect indices are updated to point into this 
!     region.
      iext = 0
 
      bnds(:)%rlenh_bg(1) = 1
      
      do icol = 1,maxcolour
      
         do n = 1,npan

            !     Start with N edge
            j = jpan
            do i = 1,ipan
               iq = indg(i,j,n)
               iqq = in_g(iq)
               ! Which processor has this point
               rproc = qproc(iqq)
               if ( rproc /= myid ) then ! Don't add points already on this proc.
                  mycol = findcolour(iqq)
                  if ( mycol == icol ) then
                     iql = indp(i,j,n)  !  Local index
                     ! Add this point to request list
                     bnds(rproc)%rlenh_fn(icol) = bnds(rproc)%rlenh_fn(icol) + 1
                     call check_bnds_alloc(rproc, iext)
                     bnds(rproc)%request_list(bnds(rproc)%rlenh_fn(icol)) = iqq
                     ! Increment extended region index
                     iext = iext + 1
                     bnds(rproc)%unpack_list(bnds(rproc)%rlenh_fn(icol)) = iext
                     in(iql) = ifull + iext
                  end if
               end if   
            end do

            !     E edge
            i = ipan
            do j = 1,jpan
               iq = indg(i,j,n)
               iqq = ie_g(iq)
               ! Which processor has this point
               rproc = qproc(iqq)
               if ( rproc /= myid ) then ! Don't add points already on this proc.
                  mycol = findcolour(iqq)
                  if ( mycol == icol ) then
                     iql = indp(i,j,n)  !  Local index
                     ! Add this point to request list
                     bnds(rproc)%rlenh_fn(icol) = bnds(rproc)%rlenh_fn(icol) + 1
                     call check_bnds_alloc(rproc, iext)
                     bnds(rproc)%request_list(bnds(rproc)%rlenh_fn(icol)) = iqq
                     ! Increment extended region index
                     iext = iext + 1
                     bnds(rproc)%unpack_list(bnds(rproc)%rlenh_fn(icol)) = iext
                     ie(iql) = ifull + iext
                  end if
               end if   
            end do
         end do ! n=1,npan
         
         if ( icol < maxcolour ) then
            bnds(:)%rlenh_bg(icol+1) = bnds(:)%rlenh_fn(icol) + 1
            bnds(:)%rlenh_fn(icol+1) = bnds(:)%rlenh_fn(icol)
         end if
      
      end do ! icol=1,maxcolour
      
      bnds(:)%rlen_bg(1) = bnds(:)%rlenh_fn(maxcolour) + 1
      bnds(:)%rlen_fn(1) = bnds(:)%rlenh_fn(maxcolour)
      
      do icol = 1,maxcolour
      
         do n = 1,npan

            !     W edge
            i = 1
            do j = 1,jpan
               ! 1D code takes care of corners separately at the end so only goes
               ! over 1,jpan here.
               iq = indg(i,j,n)
               iqq = iw_g(iq) ! iqq is the global index of the required neighbouring point.
               ! Which processor has this point
               rproc = qproc(iqq)
               if ( rproc /= myid ) then ! Don't add points already on this proc.
                  mycol = findcolour(iqq)
                  if ( mycol == icol ) then
                     iql = indp(i,j,n)  !  Local index
                     ! Add this point to request list
                     bnds(rproc)%rlen_fn(icol) = bnds(rproc)%rlen_fn(icol) + 1
                     call check_bnds_alloc(rproc, iext)
                     bnds(rproc)%request_list(bnds(rproc)%rlen_fn(icol)) = iqq
                     ! Increment extended region index
                     iext = iext + 1
                     bnds(rproc)%unpack_list(bnds(rproc)%rlen_fn(icol)) = iext
                     iw(iql) = ifull + iext
                  end if
               end if   
            end do

            !     S edge
            j = 1
            do i = 1,ipan
               iq = indg(i,j,n)
               iqq = is_g(iq)
               ! Which processor has this point
               rproc = qproc(iqq)
               if ( rproc /= myid ) then ! Don't add points already on this proc.
                  mycol = findcolour(iqq)
                  if ( mycol == icol ) then
                     iql = indp(i,j,n)  !  Local index
                     ! Add this point to request list
                     bnds(rproc)%rlen_fn(icol) = bnds(rproc)%rlen_fn(icol) + 1                     
                     call check_bnds_alloc(rproc, iext)
                     bnds(rproc)%request_list(bnds(rproc)%rlen_fn(icol)) = iqq
                     ! Increment extended region index
                     iext = iext + 1
                     bnds(rproc)%unpack_list(bnds(rproc)%rlen_fn(icol)) = iext
                     is(iql) = ifull + iext
                  end if
               end if   
            end do
         end do ! n=1,npan

         if ( icol < maxcolour ) then
            bnds(:)%rlen_bg(icol+1) = bnds(:)%rlen_fn(icol) + 1
            bnds(:)%rlen_fn(icol+1) = bnds(:)%rlen_fn(icol)
         end if
         
      end do ! icol=1,maxcolour
      
      bnds(:)%rlenx_bg(1) = bnds(:)%rlen_fn(maxcolour) + 1
      bnds(:)%rlenx_fn(1) = bnds(:)%rlen_fn(maxcolour)
      
      do icol = 1,maxcolour
      
         do n = 1,npan
            ! NE, EN
            iql = indp(ipan,jpan,n)
            iqg = indg(ipan,jpan,n)
            iqq = ine_g(iqg)
            ! Which processor has this point
            rproc = qproc(iqq)
            if ( rproc /= myid ) then ! Don't add points already on this proc.
               mycol = findcolour(iqq)
               if ( mycol == icol ) then
                  ! Add this point to request list
                  bnds(rproc)%rlenx_fn(icol) = bnds(rproc)%rlenx_fn(icol) + 1
                  call check_bnds_alloc(rproc, iext)
                  bnds(rproc)%request_list(bnds(rproc)%rlenx_fn(icol)) = iqq
                  ! Increment extended region index
                  iext = iext + 1
                  bnds(rproc)%unpack_list(bnds(rproc)%rlenx_fn(icol)) = iext
                  ine(iql) = ifull + iext
                end if
             end if
            if ( ien_g(iqg) == ine_g(iqg) ) then
               ien(iql) = ine(iql)
            else
               iqq = ien_g(iqg)
               ! Which processor has this point
               rproc = qproc(iqq)
               if ( rproc /= myid ) then ! Add to list
                  bnds(rproc)%rlenx_fn(icol) = bnds(rproc)%rlenx_fn(icol) + 1
                  call check_bnds_alloc(rproc, iext)
                  bnds(rproc)%request_list(bnds(rproc)%rlenx_fn(icol)) = iqq
                  ! Increment extended region index
                  iext = iext + 1
                  bnds(rproc)%unpack_list(bnds(rproc)%rlenx_fn(icol)) = iext
                  ien(iql) = ifull + iext
               end if
            end if
         end do

         do n = 1,npan
            ! SE, ES
            iql = indp(ipan,1,n)
            iqg = indg(ipan,1,n)
            iqq = ise_g(iqg)
            ! Which processor has this point
            rproc = qproc(iqq)
            if ( rproc /= myid ) then ! Don't add points already on this proc.
               mycol = findcolour(iqq)
               if ( mycol == icol ) then
                  ! Add this point to request list
                  bnds(rproc)%rlenx_fn(icol) = bnds(rproc)%rlenx_fn(icol) + 1
                  call check_bnds_alloc(rproc, iext)
                  bnds(rproc)%request_list(bnds(rproc)%rlenx_fn(icol)) = iqq
                  ! Increment extended region index
                  iext = iext + 1
                  bnds(rproc)%unpack_list(bnds(rproc)%rlenx_fn(icol)) = iext
                  ise(iql) = ifull + iext
                end if
             end if
            if ( ies_g(iqg) == ise_g(iqg) ) then
               ies(iql) = ise(iql)
            else
               iqq = ies_g(iqg)
               ! Which processor has this point
               rproc = qproc(iqq)
               if ( rproc /= myid ) then ! Add to list
                  bnds(rproc)%rlenx_fn(icol) = bnds(rproc)%rlenx_fn(icol) + 1
                  call check_bnds_alloc(rproc, iext)
                  bnds(rproc)%request_list(bnds(rproc)%rlenx_fn(icol)) = iqq
                  ! Increment extended region index
                  iext = iext + 1
                  bnds(rproc)%unpack_list(bnds(rproc)%rlenx_fn(icol)) = iext
                  ies(iql) = ifull + iext
               end if
            end if
         end do
         
         do n = 1,npan
            ! NW, WN
            iql = indp(1,jpan,n)
            iqg = indg(1,jpan,n)
            iqq = inw_g(iqg)
            ! Which processor has this point
            rproc = qproc(iqq)
            if ( rproc /= myid ) then ! Don't add points already on this proc.
               mycol = findcolour(iqq)
               if ( mycol == icol ) then
                  ! Add this point to request list
                  bnds(rproc)%rlenx_fn(icol) = bnds(rproc)%rlenx_fn(icol) + 1
                  call check_bnds_alloc(rproc, iext)
                  bnds(rproc)%request_list(bnds(rproc)%rlenx_fn(icol)) = iqq
                  ! Increment extended region index
                  iext = iext + 1
                  bnds(rproc)%unpack_list(bnds(rproc)%rlenx_fn(icol)) = iext
                  inw(iql) = ifull + iext
               end if
            end if
            if ( iwn_g(iqg) == inw_g(iqg) ) then
               iwn(iql) = inw(iql)
            else
               iqq = iwn_g(iqg)
               ! Which processor has this point
               rproc = qproc(iqq)
               if ( rproc /= myid ) then ! Add to list
                  bnds(rproc)%rlenx_fn(icol) = bnds(rproc)%rlenx_fn(icol) + 1
                  call check_bnds_alloc(rproc, iext)
                  bnds(rproc)%request_list(bnds(rproc)%rlenx_fn(icol)) = iqq
                  ! Increment extended region index
                  iext = iext + 1
                  bnds(rproc)%unpack_list(bnds(rproc)%rlenx_fn(icol)) = iext
                  iwn(iql) = ifull + iext
               end if
            end if
         end do

         do n = 1,npan
            ! SW, WS
            iql = indp(1,1,n)
            iqg = indg(1,1,n)
            iqq = isw_g(iqg)
            ! Which processor has this point
            rproc = qproc(iqq)
            if ( rproc /= myid ) then ! Don't add points already on this proc.
               mycol = findcolour(iqq)
               if ( mycol == icol ) then
                  ! Add this point to request list
                  bnds(rproc)%rlenx_fn(icol) = bnds(rproc)%rlenx_fn(icol) + 1
                  call check_bnds_alloc(rproc, iext)
                  bnds(rproc)%request_list(bnds(rproc)%rlenx_fn(icol)) = iqq
                  ! Increment extended region index
                  iext = iext + 1
                  bnds(rproc)%unpack_list(bnds(rproc)%rlenx_fn(icol)) = iext
                  isw(iql) = ifull + iext
               end if
            end if
            if ( iws_g(iqg) == isw_g(iqg) ) then
               iws(iql) = isw(iql)
            else
               iqq = iws_g(iqg)
               ! Which processor has this point
               rproc = qproc(iqq)
               if ( rproc /= myid ) then ! Add to list
                  bnds(rproc)%rlenx_fn(icol) = bnds(rproc)%rlenx_fn(icol) + 1
                  call check_bnds_alloc(rproc, iext)
                  bnds(rproc)%request_list(bnds(rproc)%rlenx_fn(icol)) = iqq
                  ! Increment extended region index
                  iext = iext + 1
                  bnds(rproc)%unpack_list(bnds(rproc)%rlenx_fn(icol)) = iext
                  iws(iql) = ifull + iext
               end if
            end if
         end do
         
         if ( icol < maxcolour ) then
            bnds(:)%rlenx_bg(icol+1) = bnds(:)%rlenx_fn(icol) + 1
            bnds(:)%rlenx_fn(icol+1) = bnds(:)%rlenx_fn(icol)
         end if
         
      end do ! icol=1,maxcolour

      ! Now set up the second row
      bnds(:)%rlen2 = bnds(:)%rlenx_fn(maxcolour)  ! so that they're appended.

      do n = 1,npan

         !     Start with W edge
         i = 1
         do j = 1,jpan
            iq = indg(i,j,n)
            iqq = iww_g(iq)
            ! Which processor has this point
            rproc = qproc(iqq)
            if ( rproc /= myid ) then ! Don't add points already on this proc.
               ! Add this point to request list
               bnds(rproc)%rlen2 = bnds(rproc)%rlen2 + 1    
               call check_bnds_alloc(rproc, iext)
               bnds(rproc)%request_list(bnds(rproc)%rlen2) = iqq
               ! Increment extended region index
               iext = iext + 1
               bnds(rproc)%unpack_list(bnds(rproc)%rlen2) = iext
               iql = indp(i,j,n)  !  Local index
               iww(iql) = ifull + iext
            end if   
         end do

         !     N edge
         j = jpan
         do i = 1,ipan
            iq = indg(i,j,n)
            iqq = inn_g(iq)
            ! Which processor has this point
            rproc = qproc(iqq)
            if ( rproc /= myid ) then ! Don't add points already on this proc.
               ! Add this point to request list
               bnds(rproc)%rlen2 = bnds(rproc)%rlen2 + 1 
               call check_bnds_alloc(rproc, iext)
               bnds(rproc)%request_list(bnds(rproc)%rlen2) = iqq
               ! Increment extended region index
               iext = iext + 1
               bnds(rproc)%unpack_list(bnds(rproc)%rlen2) = iext
               iql = indp(i,j,n)  !  Local index
               inn(iql) = ifull + iext
            end if   
         end do

         !     E edge
         i = ipan
         do j = 1,jpan
            iq = indg(i,j,n)
            iqq = iee_g(iq)
            ! Which processor has this point
            rproc = qproc(iqq)
            if ( rproc /= myid ) then ! Don't add points already on this proc.
               ! Add this point to request list
               bnds(rproc)%rlen2 = bnds(rproc)%rlen2 + 1 
               call check_bnds_alloc(rproc, iext)
               bnds(rproc)%request_list(bnds(rproc)%rlen2) = iqq
               ! Increment extended region index
               iext = iext + 1
               bnds(rproc)%unpack_list(bnds(rproc)%rlen2) = iext
               iql = indp(i,j,n)  !  Local index
               iee(iql) = ifull + iext
            end if   
         end do

         !     S edge
         j = 1
         do i = 1,ipan
            iq = indg(i,j,n)
            iqq = iss_g(iq)
            ! Which processor has this point
            rproc = qproc(iqq)
            if ( rproc /= myid ) then ! Don't add points already on this proc.
               ! Add this point to request list
               bnds(rproc)%rlen2 = bnds(rproc)%rlen2 + 1 
               call check_bnds_alloc(rproc, iext)
               bnds(rproc)%request_list(bnds(rproc)%rlen2) = iqq
               ! Increment extended region index
               iext = iext + 1
               bnds(rproc)%unpack_list(bnds(rproc)%rlen2) = iext
               iql = indp(i,j,n)  !  Local index
               iss(iql) = ifull + iext
            end if   
         end do
      end do ! n=1,npan

!     Now handle the second order special corner values
      do n = 1,npan
         ! NE corner, lnee, leen, lenn, lnne
         iqg = indg(ipan,jpan,n)
         if ( edge_n .and. edge_e ) then
            call fix_index2(lnee_g(n-noff),lnee,n,bnds,iext)
            call fix_index2(leen_g(n-noff),leen,n,bnds,iext)
            call fix_index2(lenn_g(n-noff),lenn,n,bnds,iext)
            call fix_index2(lnne_g(n-noff),lnne,n,bnds,iext)
         else if ( edge_e ) then
            ! Use in first because it'll be on same face.
            call fix_index2(iee_g(in_g(iqg)),lnee,n,bnds,iext)
            leen(n) = lnee(n)
            call fix_index2(ie_g(inn_g(iqg)),lenn,n,bnds,iext)
            lnne(n) = lenn(n)
         else
            ! Use ie first
            call fix_index2(in_g(iee_g(iqg)),lnee,n,bnds,iext)
            leen(n) = lnee(n)
            call fix_index2(inn_g(ie_g(iqg)),lenn,n,bnds,iext)
            lnne(n) = lenn(n)
         end if

         ! NW corner, lnww, lwwn, lwnn, lnnw
         iqg = indg(1,jpan,n)
         if ( edge_n .and. edge_w ) then
            call fix_index2(lnww_g(n-noff),lnww,n,bnds,iext)
            call fix_index2(lwwn_g(n-noff),lwwn,n,bnds,iext)
            call fix_index2(lwnn_g(n-noff),lwnn,n,bnds,iext)
            call fix_index2(lnnw_g(n-noff),lnnw,n,bnds,iext)
         else if ( edge_w ) then
            ! Use in first because it'll be on same face.
            call fix_index2(iww_g(in_g(iqg)),lnww,n,bnds,iext)
            lwwn(n) = lnww(n)
            call fix_index2(iw_g(inn_g(iqg)),lwnn,n,bnds,iext)
            lnnw(n) = lwnn(n)
         else
            ! Use iw first
            call fix_index2(in_g(iww_g(iqg)),lnww,n,bnds,iext)
            lwwn(n) = lnww(n)
            call fix_index2(inn_g(iw_g(iqg)),lwnn,n,bnds,iext)
            lnnw(n) = lwnn(n)
         end if

         ! SE corner, lsee, lees, less, lsse
         iqg = indg(ipan,1,n)
         if ( edge_s .and. edge_e ) then
            call fix_index2(lsee_g(n-noff),lsee,n,bnds,iext)
            call fix_index2(lees_g(n-noff),lees,n,bnds,iext)
            call fix_index2(less_g(n-noff),less,n,bnds,iext)
            call fix_index2(lsse_g(n-noff),lsse,n,bnds,iext)
         else if ( edge_e ) then
            ! Use is first because it will be on same face.
            call fix_index2(iee_g(is_g(iqg)),lsee,n,bnds,iext)
            lees(n) = lsee(n)
            call fix_index2(ie_g(iss_g(iqg)),less,n,bnds,iext)
            lsse(n) = less(n)
         else
            ! Use ie first
            call fix_index2(is_g(iee_g(iqg)),lsee,n,bnds,iext)
            lees(n) = lsee(n)
            call fix_index2(iss_g(ie_g(iqg)),less,n,bnds,iext)
            lsse(n) = less(n)
         end if

         ! SW corner, lsww, lwws, lwss, lssw
         iqg = indg(1,1,n)
         if ( edge_s .and. edge_w ) then
            call fix_index2(lsww_g(n-noff),lsww,n,bnds,iext)
            call fix_index2(lwws_g(n-noff),lwws,n,bnds,iext)
            call fix_index2(lwss_g(n-noff),lwss,n,bnds,iext)
            call fix_index2(lssw_g(n-noff),lssw,n,bnds,iext)
         else if ( edge_w ) then
            ! Use is first because it will be on same face.
            call fix_index2(iww_g(is_g(iqg)),lsww,n,bnds,iext)
            lwws(n) = lsww(n)
            call fix_index2(iw_g(iss_g(iqg)),lwss,n,bnds,iext)
            lssw(n) = lwss(n)
         else
            ! Use iw first
            call fix_index2(is_g(iww_g(iqg)),lsww,n,bnds,iext)
            lwws(n) = lsww(n)
            call fix_index2(iss_g(iw_g(iqg)),lwss,n,bnds,iext)
            lssw(n) = lwss(n)
         end if

      end do

      ! Indices that are missed above (should be a better way to get these)
      do n = 1,npan
         do j = 1,jpan
            iww(indp(2,j,n)) = iw(indp(1,j,n))
            iee(indp(ipan-1,j,n)) = ie(indp(ipan,j,n))
         end do
         do i = 1,ipan
            iss(indp(i,2,n)) = is(indp(i,1,n))
            inn(indp(i,jpan-1,n)) = in(indp(i,jpan,n))
         end do
      end do

      !     Set up the diagonal index arrays. Most of the points here will have
      !     already been added to copy lists above. The corners are handled
      !     separately here. This means some points may be copied twice but it's
      !     a small overhead.
      do n = 1,npan
         do j = 1,jpan
            do i = 1,ipan
               iq = indp(i,j,n)   ! Local
               ! Except at corners, ien = ine etc.
               if ( i > 1 ) then
                  inw(iq) = in(iw(iq))
                  isw(iq) = is(iw(iq))
               else
                  if ( j < jpan ) inw(iq) = iw(in(iq))
                  if ( j > 1 )    isw(iq) = iw(is(iq))
               end if
               if ( i < ipan ) then
                  ! ie will be defined
                  ine(iq) = in(ie(iq))
                  ise(iq) = is(ie(iq))
               else
                  ! i = ipan, ie will have been remapped
                  if ( j > 1 )    ise(iq) = ie(is(iq))
                  if ( j < jpan ) ine(iq) = ie(in(iq))
               end if
               if ( j > 1 ) then
                  ies(iq) = ie(is(iq))
                  iws(iq) = iw(is(iq))
               else
                  if ( i < ipan ) ies(iq)=is(ie(iq))
                  if ( i > 1 )    iws(iq)=is(iw(iq))
               end if
               if ( j < jpan ) then
                  ien(iq) = ie(in(iq))
                  iwn(iq) = iw(in(iq))
               else
                  if ( i < ipan) ien(iq) = in(ie(iq))
                  if ( i > 1 )   iwn(iq) = in(iw(iq))
               end if
            end do
         end do
      end do

      if ( iext > iextra ) then
         write(6,*) "IEXT too large", iext, iextra
         call ccmpi_abort(-1)
      end if

      ! determine neighbour processes
      allocate( neigharray_g(0:nproc-1) )
      ! default neighbour list
      neigharray_g(:) = bnds(:)%rlen2 > 0
      ! estimate maximum distance for departure points
      ! assume wind speed is less than 350 m/s
      maxdis = 350.*dt/rearth ! unit sphere
      do n = 1,npan
         j = 1
         do i = 1,ipan
            iqq = indg(i,j,n)
            call checkdistance(iqq,maxdis,neigharray_g)
         end do
         j = jpan
         do i = 1,ipan
            iqq = indg(i,j,n)
            call checkdistance(iqq,maxdis,neigharray_g)
         end do
         i = 1
         do j = 1,jpan
            iqq = indg(i,j,n)
            call checkdistance(iqq,maxdis,neigharray_g)
         end do
         i = ipan
         do j = 1,jpan
            iqq = indg(i,j,n)
            call checkdistance(iqq,maxdis,neigharray_g)
         end do
      end do
      neigharray_g(myid) = .false.
      neighnum = count( neigharray_g )
      where( neigharray_g )
         bnds(:)%len = max( bnds(:)%len, maxbuflen*maxvertlen )
      end where

      ! set-up neighbour lists
      allocate ( neighlist(neighnum) )
      allocate ( neighmap(0:nproc-1) )
      ncount = 0
      neighmap(:) = 0 ! missing
      do iproc = 1,nproc-1
         rproc = modulo(myid+iproc,nproc)
         if ( neigharray_g(rproc) ) then
            ncount = ncount + 1
            neighlist(ncount) = rproc
            neighmap(rproc) = ncount
         end if
      end do
      deallocate( neigharray_g )
      if ( ncount /= neighnum ) then
         write(6,*) "ERROR: neighnum mismatch"
         write(6,*) "neighnum, ncount ",neighnum, ncount
         call ccmpi_abort(-1)
      end if

      
      ! allocate arrays that depend on neighnum
      ! ireq needs 1 point for the MPI_Waitall which can use ireq(rreq+1)
      allocate( ireq(max(2*neighnum,1)), rlist(max(neighnum,1)) )
      allocate( dums(7,neighnum), dumr(7,neighnum) )
      allocate( status(MPI_STATUS_SIZE,size(ireq)) )

      
      ! Communicate lengths for rlenh, rlen, rlenx and rlen2
      lcomm = comm_world
      nreq = 0
      do iproc = 1,neighnum
         rproc = neighlist(iproc) ! Recv from
         if ( bnds(rproc)%rlen2 > 0 ) then
            nreq = nreq + 1
            lproc = rproc
            call MPI_IRecv( dumr(:,iproc), 7_4, ltypei, lproc, &
                            itag, lcomm, ireq(nreq), ierr )
         end if
      end do
      do iproc = neighnum,1,-1
         sproc = neighlist(iproc)  ! Send to
         if ( bnds(sproc)%rlen2 > 0 ) then
            nreq = nreq + 1
            dums(1,iproc) = bnds(sproc)%rlen2
            dums(2,iproc) = bnds(sproc)%rlenh_fn(1)
            dums(3,iproc) = bnds(sproc)%rlenh_fn(2)
            dums(4,iproc) = bnds(sproc)%rlen_fn(1)
            dums(5,iproc) = bnds(sproc)%rlen_fn(2)
            dums(6,iproc) = bnds(sproc)%rlenx_fn(1)
            dums(7,iproc) = bnds(sproc)%rlenx_fn(2)
            lproc = sproc
            call MPI_ISend( dums(:,iproc), 7_4, ltypei, lproc, &
                            itag, lcomm, ireq(nreq), ierr )
         end if
      end do
      call MPI_Waitall( nreq, ireq, status, ierr )
      do iproc = 1,neighnum
         rproc = neighlist(iproc)
         if ( bnds(rproc)%rlen2 > 0 ) then
            bnds(rproc)%slen2       = dumr(1,iproc)
            bnds(rproc)%slenh_fn(1) = dumr(2,iproc)
            bnds(rproc)%slenh_fn(2) = dumr(3,iproc)
            bnds(rproc)%slen_fn(1)  = dumr(4,iproc)
            bnds(rproc)%slen_fn(2)  = dumr(5,iproc)
            bnds(rproc)%slenx_fn(1) = dumr(6,iproc)
            bnds(rproc)%slenx_fn(2) = dumr(7,iproc)
         end if
      end do
      bnds(:)%slenh_bg(1) = 1
      bnds(:)%slenh_bg(2) = bnds(:)%slenh_fn(1) + 1
      bnds(:)%slen_bg(1)  = bnds(:)%slenh_fn(2) + 1
      bnds(:)%slen_bg(2)  = bnds(:)%slen_fn(1) + 1
      bnds(:)%slenx_bg(1) = bnds(:)%slen_fn(2) + 1
      bnds(:)%slenx_bg(2) = bnds(:)%slenx_fn(1) + 1
      
      deallocate( dums, dumr )

      
      ! Now, for each processor send the list of points I want.
      ! The state of being a neighbour is reflexive so only expect to
      ! recv from those processors I send to (are there grid arrangements for
      ! which this would not be true?)
      lcomm = comm_world
      nreq = 0
      do iproc = 1,neighnum
         rproc = neighlist(iproc)  ! Recv from
         llen = bnds(rproc)%slen2
         if ( llen > 0 ) then
            allocate( bnds(rproc)%send_list(llen) ) 
            nreq = nreq + 1         
            lproc = rproc
            call MPI_IRecv( bnds(rproc)%send_list, llen, ltypei, lproc, &
                            itag, lcomm, ireq(nreq), ierr )
         end if
      end do
      do iproc = 1,neighnum
         sproc = neighlist(iproc)  ! Send to
         ! Send list of requests
         llen = bnds(sproc)%rlen2
         if ( llen > 0 ) then
            lproc = sproc
            nreq = nreq + 1
            call MPI_ISend( bnds(sproc)%request_list, llen, &
                 ltypei, lproc, itag, lcomm, ireq(nreq), ierr )
         end if
      end do      
      call MPI_Waitall( nreq, ireq, status, ierr )

!     Now deallocate arrays that are no longer needed
      allocate( dumi(maxbuflen) )
      do iproc = 1,neighnum
         rproc = neighlist(iproc)  ! Recv from
         if ( bnds(rproc)%rlen2 > 0 ) then
            deallocate( bnds(rproc)%request_list )
            dumi(1:bnds(rproc)%rlen2) = bnds(rproc)%unpack_list(1:bnds(rproc)%rlen2)
            deallocate( bnds(rproc)%unpack_list )
            allocate( bnds(rproc)%unpack_list(bnds(rproc)%rlen2) )
            bnds(rproc)%unpack_list(1:bnds(rproc)%rlen2) = dumi(1:bnds(rproc)%rlen2)
         end if   
      end do
      deallocate( dumi )
      
      
!     Start of UV section

!     In the first pass through, set up list of points to be requested from
!     other processors. In the 1D code values on the same processor are
!     copied only if they have to be swapped.
!     This only makes a difference on 1, 2 or 3 processors.

      iextu = 0
      iextv = 0

      ! save start of isv indices
      bnds(:)%rlen_sv_bg = 1

      !     S edge, V
      j = 1
      do n = 1,npan
         do i = 1,ipan
            iqg = indg(i,j,n)
            iqq = is_g(iqg)
            rproc = qproc(iqq)
            swap = edge_s .and. swap_s(n-noff)
            if ( rproc /= myid .or. swap ) then
               bnds(rproc)%rlen_wu_fn = bnds(rproc)%rlen_wu_fn + 1                
               call check_bnds_alloc(rproc, iextv)
               bnds(rproc)%request_list_uv(bnds(rproc)%rlen_wu_fn) = -iqq
               ! Increment extended region index
               iextv = iextv + 1
               bnds(rproc)%unpack_list_uv(bnds(rproc)%rlen_wu_fn) = -iextv
               iql = indp(i,j,n)  !  Local index
               isv(iql) = ifull + iextv
               bnds(rproc)%uv_swap(bnds(rproc)%rlen_wu_fn) = swap
               bnds(rproc)%uv_neg(bnds(rproc)%rlen_wu_fn) = 1.
            end if   
         end do
      end do
         
      !     Start with W edge, U values
      i = 1
      do n = 1,npan
         do j = 1,jpan
            iqg = indg(i,j,n)
            iqq = iw_g(iqg)
            ! Which processor has this point
            rproc = qproc(iqq)
            ! Only need to add to bounds region if it's on another processor
            ! or if it's on this processor and needs to be swapped.
            ! Decide if u/v need to be swapped. My face is n-noff
            swap = edge_w .and. swap_w(n-noff)
            if ( rproc /= myid .or. swap ) then
               bnds(rproc)%rlen_wu_fn = bnds(rproc)%rlen_wu_fn + 1 
               call check_bnds_alloc(rproc, iextu)
               bnds(rproc)%request_list_uv(bnds(rproc)%rlen_wu_fn) = iqq
               ! Increment extended region index
               iextu = iextu + 1
               bnds(rproc)%unpack_list_uv(bnds(rproc)%rlen_wu_fn) = iextu
               iql = indp(i,j,n)  !  Local index
               iwu(iql) = ifull + iextu
               bnds(rproc)%uv_swap(bnds(rproc)%rlen_wu_fn) = swap
               bnds(rproc)%uv_neg(bnds(rproc)%rlen_wu_fn) = 1.
            end if  
         end do
      end do

      ! save start of inv indices
      bnds(:)%rlen_nv_bg = bnds(:)%rlen_wu_fn + 1
      bnds(:)%rlen_eu_fn = bnds(:)%rlen_wu_fn

      !     N edge (V)
      j = jpan
      do n = 1,npan
         do i = 1,ipan
            iqg = indg(i,j,n)
            iqq = in_g(iqg)
            rproc = qproc(iqq)
            swap = edge_n .and. swap_n(n-noff)
            if ( rproc /= myid .or. swap ) then
               ! Add this point to request list
               bnds(rproc)%rlen_eu_fn = bnds(rproc)%rlen_eu_fn + 1
               call check_bnds_alloc(rproc, iextv)
               ! to show that this is v rather than u, flip sign
               bnds(rproc)%request_list_uv(bnds(rproc)%rlen_eu_fn) = -iqq
               ! Increment extended region index
               iextv = iextv + 1
               bnds(rproc)%unpack_list_uv(bnds(rproc)%rlen_eu_fn) = -iextv
               iql = indp(i,j,n)  !  Local index
               inv(iql) = ifull + iextv
               bnds(rproc)%uv_swap(bnds(rproc)%rlen_eu_fn) = swap
               bnds(rproc)%uv_neg(bnds(rproc)%rlen_eu_fn) = 1.
            end if
         end do
      end do

      !     E edge, U
      i = ipan
      do n = 1,npan
         do j = 1,jpan
            iqg = indg(i,j,n)
            iqq = ie_g(iqg)
            rproc = qproc(iqq)
            swap = edge_e .and. swap_e(n-noff)
            if ( rproc /= myid .or. swap ) then
               ! Add this point to request list
               bnds(rproc)%rlen_eu_fn = bnds(rproc)%rlen_eu_fn + 1 
               call check_bnds_alloc(rproc, iextu)
               bnds(rproc)%request_list_uv(bnds(rproc)%rlen_eu_fn) = iqq
               ! Increment extended region index
               iextu = iextu + 1
               bnds(rproc)%unpack_list_uv(bnds(rproc)%rlen_eu_fn) = iextu
               iql = indp(i,j,n)  !  Local index
               ieu(iql) = ifull + iextu
               bnds(rproc)%uv_swap(bnds(rproc)%rlen_eu_fn) = swap
               bnds(rproc)%uv_neg(bnds(rproc)%rlen_eu_fn) = 1.
            end if   
         end do
      end do

      ! Second pass
      ieeu = iee
      iwwu = iww
      innv = inn
      issv = iss
      innu = inn
      ieev = iee

      ! save start of issv indices
      bnds(:)%rlen_ssv_bg = bnds(:)%rlen_eu_fn + 1
      bnds(:)%rlen_wwu_fn = bnds(:)%rlen_eu_fn

      !     SS edge, V
      j = 1
      do n = 1,npan
         do i = 1,ipan
            iqg = indg(i,j,n)
            iqq = iss_g(iqg)
            rproc = qproc(iqq)
            swap = edge_s .and. swap_s(n-noff)
            if ( rproc /= myid .or. swap ) then
               bnds(rproc)%rlen_wwu_fn = bnds(rproc)%rlen_wwu_fn + 1
               call check_bnds_alloc(rproc, iextv)
               bnds(rproc)%request_list_uv(bnds(rproc)%rlen_wwu_fn) = -iqq
               ! Increment extended region index
               iextv = iextv + 1
               bnds(rproc)%unpack_list_uv(bnds(rproc)%rlen_wwu_fn) = -iextv
               iql = indp(i,j,n)  !  Local index
               issv(iql) = ifull + iextv
               bnds(rproc)%uv_swap(bnds(rproc)%rlen_wwu_fn) = swap
               bnds(rproc)%uv_neg(bnds(rproc)%rlen_wwu_fn) = 1.
            end if   
         end do
      end do ! n=1,npan

      !     Start with WW edge, U values
      i = 1
      do n = 1,npan
         do j = 1,jpan
            iqg = indg(i,j,n)
            iqq = iww_g(iqg)
            ! Which processor has this point
            rproc = qproc(iqq)
            swap = edge_w .and. swap_w(n-noff)
            if ( rproc /= myid .or. swap ) then
               ! Add this point to request list
               bnds(rproc)%rlen_wwu_fn = bnds(rproc)%rlen_wwu_fn + 1
               call check_bnds_alloc(rproc, iextu)
               bnds(rproc)%request_list_uv(bnds(rproc)%rlen_wwu_fn) = iqq
               ! Increment extended region index
               iextu = iextu + 1
               bnds(rproc)%unpack_list_uv(bnds(rproc)%rlen_wwu_fn) = iextu
               iql = indp(i,j,n)  !  Local index
               iwwu(iql) = ifull + iextu
               ! Decide if u/v need to be swapped. My face is n-noff
               bnds(rproc)%uv_swap(bnds(rproc)%rlen_wwu_fn) = swap
               bnds(rproc)%uv_neg(bnds(rproc)%rlen_wwu_fn) = 1.
            end if   
         end do
      end do

      ! save start of innv indices
      bnds(:)%rlen_nnv_bg = bnds(:)%rlen_wwu_fn + 1
      bnds(:)%rlen_eeu_fn = bnds(:)%rlen_wwu_fn

      !     NN edge (V)
      j = jpan
      do n = 1,npan
         do i = 1,ipan
            iqg = indg(i,j,n)
            iqq = inn_g(iqg)
            rproc = qproc(iqq)
            swap = edge_n .and. swap_n(n-noff)
            if ( rproc /= myid .or. swap ) then
               ! Add this point to request list
               bnds(rproc)%rlen_eeu_fn = bnds(rproc)%rlen_eeu_fn + 1                
               call check_bnds_alloc(rproc, iextv)
               ! to show that this is v rather than u, flip sign
               bnds(rproc)%request_list_uv(bnds(rproc)%rlen_eeu_fn) = -iqq
               ! Increment extended region index
               iextv = iextv + 1
               bnds(rproc)%unpack_list_uv(bnds(rproc)%rlen_eeu_fn) = -iextv
               iql = indp(i,j,n)  !  Local index
               innv(iql) = ifull + iextv
               bnds(rproc)%uv_swap(bnds(rproc)%rlen_eeu_fn) = swap
               bnds(rproc)%uv_neg(bnds(rproc)%rlen_eeu_fn) = 1.
            end if   
         end do
      end do

      !     EE edge, U
      i = ipan
      do n = 1,npan
         do j = 1,jpan
            iqg = indg(i,j,n)
            iqq = iee_g(iqg)
            rproc = qproc(iqq)
            swap = edge_e .and. swap_e(n-noff)
            if ( rproc /= myid .or. swap ) then
               ! Add this point to request list
               bnds(rproc)%rlen_eeu_fn = bnds(rproc)%rlen_eeu_fn + 1
               call check_bnds_alloc(rproc, iextu)
               bnds(rproc)%request_list_uv(bnds(rproc)%rlen_eeu_fn) = iqq
               ! Increment extended region index
               iextu = iextu + 1
               bnds(rproc)%unpack_list_uv(bnds(rproc)%rlen_eeu_fn) = iextu
               iql = indp(i,j,n)  !  Local index
               ieeu(iql) = ifull + iextu
               bnds(rproc)%uv_swap(bnds(rproc)%rlen_eeu_fn) = swap
               bnds(rproc)%uv_neg(bnds(rproc)%rlen_eeu_fn) = 1.
            end if   
         end do
      end do

      ! Third pass
      ! save start of isu indices
      bnds(:)%rlen_su_bg = bnds(:)%rlen_eeu_fn + 1
      bnds(:)%rlen_ev_fn = bnds(:)%rlen_eeu_fn

      !     S edge, U
      j = 1
      do n = 1,npan
         do i = 1,ipan
            iqg = indg(i,j,n)
            iqq = is_g(iqg)
            rproc = qproc(iqq)
            swap = edge_s .and. swap_s(n-noff)
            if ( rproc /= myid .or. swap ) then
               bnds(rproc)%rlen_ev_fn = bnds(rproc)%rlen_ev_fn + 1 
               call check_bnds_alloc(rproc, iextu)
               bnds(rproc)%request_list_uv(bnds(rproc)%rlen_ev_fn) = iqq
               ! Increment extended region index
               iextu = iextu + 1
               bnds(rproc)%unpack_list_uv(bnds(rproc)%rlen_ev_fn) = iextu
               iql = indp(i,j,n)  !  Local index
               isu(iql) = ifull + iextu
               bnds(rproc)%uv_swap(bnds(rproc)%rlen_ev_fn) = swap
               if ( swap ) then
                  bnds(rproc)%uv_neg(bnds(rproc)%rlen_ev_fn) = -1.
               else
                  bnds(rproc)%uv_neg(bnds(rproc)%rlen_ev_fn) = 1. 
               end if   
            end if   
         end do
      end do ! n=1,npan

      !     Start with W edge, V values
      i = 1
      do n = 1,npan
         do j = 1,jpan
            iqg = indg(i,j,n)
            iqq = iw_g(iqg)
            ! Which processor has this point
            rproc = qproc(iqq)
            ! Only need to add to bounds region if it's on another processor
            ! or if it's on this processor and needs to be swapped.
            ! Decide if u/v need to be swapped. My face is n-noff
            swap = edge_w .and. swap_w(n-noff)
            if ( rproc /= myid .or. swap ) then
               bnds(rproc)%rlen_ev_fn = bnds(rproc)%rlen_ev_fn + 1 
               call check_bnds_alloc(rproc, iextv)
               ! to show that this is v rather than u, flip sign
               bnds(rproc)%request_list_uv(bnds(rproc)%rlen_ev_fn) = -iqq
               ! Increment extended region index
               iextv = iextv + 1
               bnds(rproc)%unpack_list_uv(bnds(rproc)%rlen_ev_fn) = -iextv
               iql = indp(i,j,n)  !  Local index
               iwv(iql) = ifull + iextv
               bnds(rproc)%uv_swap(bnds(rproc)%rlen_ev_fn) = swap
               if ( swap ) then
                  bnds(rproc)%uv_neg(bnds(rproc)%rlen_ev_fn) = -1.
               else
                  bnds(rproc)%uv_neg(bnds(rproc)%rlen_ev_fn) = 1. 
               end if   
            end if   
         end do
      end do

      !     N edge (U)
      j = jpan
      do n = 1,npan
         do i = 1,ipan
            iqg = indg(i,j,n)
            iqq = in_g(iqg)
            rproc = qproc(iqq)
            swap = edge_n .and. swap_n(n-noff)
            if ( rproc /= myid .or. swap ) then
               ! Add this point to request list
               bnds(rproc)%rlen_ev_fn = bnds(rproc)%rlen_ev_fn + 1
               call check_bnds_alloc(rproc, iextu)
               bnds(rproc)%request_list_uv(bnds(rproc)%rlen_ev_fn) = iqq
               ! Increment extended region index
               iextu = iextu + 1
               bnds(rproc)%unpack_list_uv(bnds(rproc)%rlen_ev_fn) = iextu
               iql = indp(i,j,n)  !  Local index
               inu(iql) = ifull + iextu
               bnds(rproc)%uv_swap(bnds(rproc)%rlen_ev_fn) = swap
               if ( swap ) then
                  bnds(rproc)%uv_neg(bnds(rproc)%rlen_ev_fn) = -1.
               else
                  bnds(rproc)%uv_neg(bnds(rproc)%rlen_ev_fn) = 1. 
               end if   
            end if   
         end do
      end do

      !     E edge, V
      i = ipan
      do n = 1,npan
         do j = 1,jpan
            iqg = indg(i,j,n)
            iqq = ie_g(iqg)
            rproc = qproc(iqq)
            swap = edge_e .and. swap_e(n-noff)
            if ( rproc /= myid .or. swap ) then
               ! Add this point to request list
               bnds(rproc)%rlen_ev_fn = bnds(rproc)%rlen_ev_fn + 1 
               call check_bnds_alloc(rproc, iextv)
               ! to show that this is v rather than u, flip sign
               bnds(rproc)%request_list_uv(bnds(rproc)%rlen_ev_fn) = -iqq
               ! Increment extended region index
               iextv = iextv + 1
               bnds(rproc)%unpack_list_uv(bnds(rproc)%rlen_ev_fn) = -iextv
               iql = indp(i,j,n)  !  Local index
               iev(iql) = ifull + iextv
               bnds(rproc)%uv_swap(bnds(rproc)%rlen_ev_fn) = swap
               if ( swap ) then
                  bnds(rproc)%uv_neg(bnds(rproc)%rlen_ev_fn) = -1.
               else
                  bnds(rproc)%uv_neg(bnds(rproc)%rlen_ev_fn) = 1. 
               end if   
            end if   
         end do
      end do

      ! save start of innu indices
      bnds(:)%rlen_nnu_bg = bnds(:)%rlen_ev_fn + 1
      bnds(:)%rlen_eev_fn = bnds(:)%rlen_ev_fn
      
      !     NN edge (U)
      j = jpan
      do n = 1,npan
         do i = 1,ipan
            iqg = indg(i,j,n)
            iqq = inn_g(iqg)
            rproc = qproc(iqq)
            swap = edge_n .and. swap_n(n-noff)
            if ( rproc /= myid .or. swap ) then
               ! Add this point to request list
               bnds(rproc)%rlen_eev_fn = bnds(rproc)%rlen_eev_fn + 1 
               call check_bnds_alloc(rproc, iextu)
               bnds(rproc)%request_list_uv(bnds(rproc)%rlen_eev_fn) = iqq
               ! Increment extended region index
               iextu = iextu + 1
               bnds(rproc)%unpack_list_uv(bnds(rproc)%rlen_eev_fn) = iextu
               iql = indp(i,j,n)  !  Local index
               innu(iql) = ifull + iextu
               bnds(rproc)%uv_swap(bnds(rproc)%rlen_eev_fn) = swap
               if ( swap ) then
                  bnds(rproc)%uv_neg(bnds(rproc)%rlen_eev_fn) = -1.
               else
                  bnds(rproc)%uv_neg(bnds(rproc)%rlen_eev_fn) = 1. 
               end if   
            end if   
         end do
      end do

      !     EE edge, V
      i = ipan
      do n = 1,npan
         do j = 1,jpan
            iqg = indg(i,j,n)
            iqq = iee_g(iqg)
            rproc = qproc(iqq)
            swap = edge_e .and. swap_e(n-noff)
            if ( rproc /= myid .or. swap ) then
               ! Add this point to request list
               bnds(rproc)%rlen_eev_fn = bnds(rproc)%rlen_eev_fn + 1
               call check_bnds_alloc(rproc, iextv)
               ! to show that this is v rather than u, flip sign
               bnds(rproc)%request_list_uv(bnds(rproc)%rlen_eev_fn) = -iqq
               ! Increment extended region index
               iextv = iextv + 1
               bnds(rproc)%unpack_list_uv(bnds(rproc)%rlen_eev_fn) = -iextv
               iql = indp(i,j,n)  !  Local index
               ieev(iql) = ifull + iextv
               bnds(rproc)%uv_swap(bnds(rproc)%rlen_eev_fn) = swap
               if ( swap ) then
                  bnds(rproc)%uv_neg(bnds(rproc)%rlen_eev_fn) = -1.
               else
                  bnds(rproc)%uv_neg(bnds(rproc)%rlen_eev_fn) = 1. 
               end if   
            end if   
         end do
      end do
      
      if ( iextu > iextra ) then
         write(6,*) "IEXTU too large", iextu, iextra
         call ccmpi_abort(-1)
      end if

      if ( iextv > iextra ) then
         write(6,*) "IEXTV too large", iextv, iextra
         call ccmpi_abort(-1)
      end if

      
      allocate( dums(6,neighnum), dumr(6,neighnum) )
      
      ! Communicate lengths for rlen_uv and rlenx_uv, etc
      lcomm = comm_world
      nreq = 0
      do iproc = 1,neighnum
         rproc = neighlist(iproc)  ! Recv from
         if ( bnds(rproc)%rlen_eev_fn > 0 ) then
            nreq = nreq + 1
            lproc = rproc
            call MPI_IRecv( dumr(:,iproc), 6_4, ltypei, lproc, &
                 itag, lcomm, ireq(nreq), ierr )
         end if
      end do
      do iproc = neighnum,1,-1
         sproc = neighlist(iproc)  ! Send to
         if ( bnds(sproc)%rlen_eev_fn > 0 ) then
            nreq = nreq + 1
            dums(1,iproc) = bnds(sproc)%rlen_wu_fn
            dums(2,iproc) = bnds(sproc)%rlen_eu_fn
            dums(3,iproc) = bnds(sproc)%rlen_wwu_fn
            dums(4,iproc) = bnds(sproc)%rlen_eeu_fn
            dums(5,iproc) = bnds(sproc)%rlen_ev_fn
            dums(6,iproc) = bnds(sproc)%rlen_eev_fn
            lproc = sproc
            call MPI_ISend( dums(:,iproc), 6_4, ltypei, lproc, &
                 itag, lcomm, ireq(nreq), ierr )
         end if
      end do
      call MPI_Waitall( nreq, ireq, status, ierr )
      do iproc = 1,neighnum
         rproc = neighlist(iproc)
         if ( bnds(rproc)%rlen_eev_fn > 0 ) then
            bnds(rproc)%slen_wu_fn  = dumr(1,iproc)
            bnds(rproc)%slen_eu_fn  = dumr(2,iproc)
            bnds(rproc)%slen_wwu_fn = dumr(3,iproc)
            bnds(rproc)%slen_eeu_fn = dumr(4,iproc)
            bnds(rproc)%slen_ev_fn  = dumr(5,iproc)
            bnds(rproc)%slen_eev_fn = dumr(6,iproc)
         end if
      end do
      bnds(:)%slen_sv_bg  = 1
      bnds(:)%slen_nv_bg  = bnds(:)%slen_wu_fn  + 1
      bnds(:)%slen_ssv_bg = bnds(:)%slen_eu_fn  + 1
      bnds(:)%slen_nnv_bg = bnds(:)%slen_wwu_fn + 1
      bnds(:)%slen_su_bg  = bnds(:)%slen_eeu_fn + 1
      bnds(:)%slen_nnu_bg = bnds(:)%slen_ev_fn  + 1
      
      deallocate( dumr, dums )

      
      ! Now, for each processor send the list of points I want.
      ! Also have to send the swap list
      lcomm = comm_world
      nreq = 0
      do iproc = 1,neighnum
         rproc = neighlist(iproc)
         llen = bnds(rproc)%slen_eev_fn
         if ( llen > 0 ) then
            allocate( bnds(rproc)%send_list_uv(llen) ) 
            nreq = nreq + 1
            lproc = rproc
            call MPI_IRecv( bnds(rproc)%send_list_uv, llen, ltypei, lproc, &
                  itag, lcomm, ireq(nreq), ierr )
         end if
      end do
      do iproc = neighnum,1,-1
         sproc = neighlist(iproc)
         if ( bnds(sproc)%rlen_eev_fn > 0 ) then
            ! Send list of requests
            nreq = nreq + 1
            llen = bnds(sproc)%rlen_eev_fn
            lproc = sproc
            call MPI_ISend( bnds(sproc)%request_list_uv, llen, ltypei, lproc, &
                 itag, lcomm, ireq(nreq), ierr )
         end if
      end do
      call MPI_Waitall( nreq, ireq, status, ierr )

      ! Deallocate arrays that are no longer needed
      ! Note that the request_list for myid is still allocated
      allocate( dumi(maxbuflen) )
      do iproc = 1,neighnum
         rproc = neighlist(iproc)
         if ( bnds(rproc)%rlen_eev_fn > 0 ) then
            deallocate( bnds(rproc)%request_list_uv ) 
            dumi(1:bnds(rproc)%rlen_eev_fn) = bnds(rproc)%unpack_list_uv(1:bnds(rproc)%rlen_eev_fn)
            deallocate( bnds(rproc)%unpack_list_uv )
            allocate( bnds(rproc)%unpack_list_uv(bnds(rproc)%rlen_eev_fn) )
            bnds(rproc)%unpack_list_uv(1:bnds(rproc)%rlen_eev_fn) = dumi(1:bnds(rproc)%rlen_eev_fn)
         end if
      end do
      deallocate( dumi )

      
      allocate( dumsl(maxbuflen,neighnum), dumrl(maxbuflen,neighnum) )
      
      ! Only send the swap list once
      lcomm = comm_world
      nreq = 0
      do iproc = 1,neighnum
         rproc = neighlist(iproc) ! Recv from
         llen = bnds(rproc)%slen_eev_fn
         if ( llen > 0 ) then
            nreq = nreq + 1
            lproc = rproc
            call MPI_IRecv( dumrl(:,iproc), llen, MPI_LOGICAL, &
                  lproc, itag, lcomm, ireq(nreq), ierr )
         end if
      end do
      do iproc = neighnum,1,-1
         sproc = neighlist(iproc) ! Send to
         llen = bnds(sproc)%rlen_eev_fn
         if ( llen > 0 ) then
            ! Send list of requests
            nreq = nreq + 1
            lproc = sproc
            dumsl(1:bnds(sproc)%rlen_eev_fn,iproc) = bnds(sproc)%uv_swap(1:bnds(sproc)%rlen_eev_fn)
            call MPI_ISend( dumsl(:,iproc), llen, MPI_LOGICAL, &
                  lproc, itag, lcomm, ireq(nreq), ierr )
         end if
      end do
      call MPI_Waitall( nreq, ireq, status, ierr )
      
      do iproc = 1,neighnum
         rproc = neighlist(iproc)
         if ( bnds(rproc)%slen_eev_fn > 0 ) then
            allocate( bnds(rproc)%send_swap(bnds(rproc)%slen_eev_fn) ) 
            bnds(rproc)%send_swap(1:bnds(rproc)%slen_eev_fn) = dumrl(1:bnds(rproc)%slen_eev_fn,iproc)
         end if
         if ( bnds(rproc)%rlen_eev_fn > 0 ) then
            deallocate( bnds(rproc)%uv_swap )
         end if    
      end do
      
      deallocate( dumrl, dumsl )
      

      ! Only send the neg list once
      lcomm = comm_world
      nreq = 0
      do iproc = 1,neighnum
         rproc = neighlist(iproc)
         llen = bnds(rproc)%slen_eev_fn
         if ( llen > 0 ) then
            allocate( bnds(rproc)%send_neg(llen) ) 
            nreq = nreq + 1
            lproc = rproc
            call MPI_IRecv( bnds(rproc)%send_neg, llen, ltyper, &
                  lproc, itag, lcomm, ireq(nreq), ierr )
         end if
      end do
      do iproc = neighnum,1,-1
         sproc = neighlist(iproc)
         llen = bnds(sproc)%rlen_eev_fn
         if ( llen > 0 ) then
            ! Send list of requests
            nreq = nreq + 1
            lproc = sproc
            call MPI_ISend( bnds(sproc)%uv_neg, llen, ltyper,&
                  lproc, itag, lcomm, ireq(nreq), ierr )
         end if
      end do
      call MPI_Waitall( nreq, ireq, status, ierr )
      
      
      deallocate( status )
      
      
      do iproc = 1,neighnum
         rproc = neighlist(iproc)
         if ( bnds(rproc)%rlen_eev_fn > 0 ) then
           deallocate( bnds(rproc)%uv_neg )
         end if
      end do
      
      
      ! Flag that all messages have been cleared
      nreq = 0
      rreq = 0
      
      
      ! Indices that are missed above (should be a better way to get these)
      do n = 1,npan
         do j = 1,jpan
            iwwu(indp(2,j,n)) = iwu(indp(1,j,n))
            ieeu(indp(ipan-1,j,n)) = ieu(indp(ipan,j,n))
            ieev(indp(ipan-1,j,n)) = iev(indp(ipan,j,n))
         end do
         do i = 1,ipan
            issv(indp(i,2,n)) = isv(indp(i,1,n))
            innv(indp(i,jpan-1,n)) = inv(indp(i,jpan,n))
            innu(indp(i,jpan-1,n)) = inu(indp(i,jpan,n))
         end do
      end do

      !  At the moment send_lists use global indices. Convert these to local.
      do iproc = neighnum,1,-1
         sproc = neighlist(iproc)  ! Send to
         do iq = 1,bnds(sproc)%slen2
            ! send_list(iq) is global point index, i, j, n are local
            iqq = bnds(sproc)%send_list(iq)
            call indv_mpi(iqq,i,j,n)
            bnds(sproc)%send_list(iq) = indp(i,j,n)
         end do
         do iq = 1,bnds(sproc)%slen_eev_fn
            ! send_list(iq) is global point index, i, j, n are local
            ! Use abs because sign is used as u/v flag
            iqq = abs(bnds(sproc)%send_list_uv(iq))
            call indv_mpi(iqq,i,j,n)
            bnds(sproc)%send_list_uv(iq) = sign(indp(i,j,n),bnds(sproc)%send_list_uv(iq))
         end do
      end do
      if ( bnds(myid)%rlen2 /= 0 ) then
         write(6,*) "ERROR: Invalid rlen2 for myid"
         call ccmpi_abort(-1)
      end if
      do iq = 1,bnds(myid)%rlen_eev_fn
         iqq = abs(bnds(myid)%request_list_uv(iq))
         call indv_mpi(iqq,i,j,n)
         bnds(myid)%request_list_uv(iq) = sign(indp(i,j,n),bnds(myid)%request_list_uv(iq))
      end do

      
      ! Allocate buffer arrays
      do iproc = 1,neighnum
         rproc = neighlist(iproc)
         if ( bnds(rproc)%len > 0 ) then
            allocate( bnds(rproc)%rbuf(nagg*bnds(rproc)%len) )
            allocate( bnds(rproc)%sbuf(nagg*bnds(rproc)%len) )
            bnds(rproc)%rbuflen = nagg*bnds(rproc)%len
            allocate( bnds(rproc)%r8buf(nagg*bnds(rproc)%len) )
            allocate( bnds(rproc)%s8buf(nagg*bnds(rproc)%len) )
            bnds(rproc)%sbuflen = nagg*bnds(rproc)%len     
         end if
      end do
      

      ! Final check for values that haven't been set properly
      do n = 1,npan
         do j = 1,jpan
            do i = 1,ipan
               iq = indp(i,j,n)
               call check_set( in(iq), "IN", i, j, n, iq)
               call check_set( is(iq), "IS", i, j, n, iq)
               call check_set( iw(iq), "IW", i, j, n, iq)
               call check_set( ie(iq), "IE", i, j, n, iq)
               call check_set( inn(iq), "INN", i, j, n, iq)
               call check_set( iss(iq), "ISS", i, j, n, iq)
               call check_set( iww(iq), "IWW", i, j, n, iq)
               call check_set( iee(iq), "IEE", i, j, n, iq)
               call check_set( ien(iq), "IEN", i, j, n, iq)
               call check_set( ine(iq), "INE", i, j, n, iq)
               call check_set( ise(iq), "ISE", i, j, n, iq)
               call check_set( iwn(iq), "IWN", i, j, n, iq)
               call check_set( inw(iq), "INW", i, j, n, iq)
               call check_set( isw(iq), "ISW", i, j, n, iq)
               call check_set( ies(iq), "IES", i, j, n, iq)
               call check_set( iws(iq), "IWS", i, j, n, iq)
               call check_set( ieu(iq), "IEU", i, j, n, iq)
               call check_set( iwu(iq), "IWU", i, j, n, iq)
               call check_set( inv(iq), "INV", i, j, n, iq)
               call check_set( isv(iq), "ISV", i, j, n, iq)
               call check_set( iev(iq), "IEV", i, j, n, iq)
               call check_set( iwv(iq), "IWV", i, j, n, iq)
               call check_set( inu(iq), "INU", i, j, n, iq)
               call check_set( isu(iq), "ISU", i, j, n, iq)
               call check_set( ieeu(iq), "IEEU", i, j, n, iq)
               call check_set( iwwu(iq), "IWWU", i, j, n, iq)
               call check_set( innv(iq), "INNV", i, j, n, iq)
               call check_set( issv(iq), "ISSV", i, j, n, iq)
               call check_set( innu(iq), "INNU", i, j, n, iq)
               call check_set( ieev(iq), "IEEV", i, j, n, iq)
            end do
         end do
         call check_set( lwws(n), "LWWS", 1, 1, n, 1)
         call check_set( lwss(n), "LWSS", 1, 1, n, 1)
         call check_set( lees(n), "LEES", 1, 1, n, 1)
         call check_set( less(n), "LESS", 1, 1, n, 1)
         call check_set( lwwn(n), "LWWN", 1, 1, n, 1)
         call check_set( lwnn(n), "LWNN", 1, 1, n, 1)
         call check_set( leen(n), "LEEN", 1, 1, n, 1)
         call check_set( lenn(n), "LENN", 1, 1, n, 1)
         call check_set( lsww(n), "LSWW", 1, 1, n, 1)
         call check_set( lssw(n), "LSSW", 1, 1, n, 1)
         call check_set( lsee(n), "LSEE", 1, 1, n, 1)
         call check_set( lsse(n), "LSSE", 1, 1, n, 1)
         call check_set( lnww(n), "LNWW", 1, 1, n, 1)
         call check_set( lnnw(n), "LNNW", 1, 1, n, 1)
         call check_set( lnee(n), "LNEE", 1, 1, n, 1)
         call check_set( lnne(n), "LNNE", 1, 1, n, 1)
      end do

   end subroutine bounds_setup

   subroutine check_bnds_alloc(rproc, iext)
      integer, intent(in) :: rproc
      integer, intent(in) :: iext
      integer :: testlen

!     Allocate the components of the bnds array. It's too much work to
!     get the exact sizes, so allocate a fixed size for each case where
!     there's an interaction.
      if ( bnds(rproc)%len == 0 ) then
         ! Not allocated yet.
         if (rproc /= myid) then
            allocate ( bnds(rproc)%request_list(maxbuflen) )
            allocate ( bnds(rproc)%unpack_list(maxbuflen) )
         end if
         allocate( bnds(rproc)%request_list_uv(maxbuflen) )
         allocate( bnds(rproc)%unpack_list_uv(maxbuflen) )
         allocate( bnds(rproc)%uv_swap(maxbuflen) )
         allocate( bnds(rproc)%uv_neg(maxbuflen) )
         bnds(rproc)%uv_neg = 1.
         bnds(rproc)%len = maxbuflen*maxvertlen
      else
         ! Just check length
         testlen = max( bnds(rproc)%rlenh_fn(1), bnds(rproc)%rlenh_fn(2), bnds(rproc)%rlen_fn(1), bnds(rproc)%rlen_fn(2), &
                        bnds(rproc)%rlenx_fn(1), bnds(rproc)%rlenx_fn(2), bnds(rproc)%rlen2 )*maxvertlen
         if ( testlen >  bnds(rproc)%len ) then
            write(6,*) "Error, maximum length error in check_bnds_alloc"
            write(6,*) myid, rproc, testlen,  bnds(rproc)%len, maxvertlen
            write(6,*) bnds(rproc)%rlenh_fn(1), bnds(rproc)%rlenh_fn(2), bnds(rproc)%rlen_fn(1), bnds(rproc)%rlen_fn(2)
            write(6,*) bnds(rproc)%rlenx_fn(1), bnds(rproc)%rlenx_fn(2), bnds(rproc)%rlen2
            call ccmpi_abort(-1)
         end if
         if ( iext >= iextra ) then
            write(6,*) "Error, iext maximum length error in check_bnds_alloc"
            write(6,*) myid, iext, iextra
            call ccmpi_abort(-1)
         end if
      end if
      
   end subroutine check_bnds_alloc
   
   subroutine checkdistance(iqq,maxdis,neigharray_g)

      use xyzinfo_m, only : x_g, y_g, z_g
   
      integer, intent(in) :: iqq
      integer :: iqg, i_g, j_g, iproc, ioff_g, joff_g
      real, intent(in) :: maxdis
      real :: disarray_g
      logical, dimension(0:nproc-1), intent(inout) :: neigharray_g
   
      do joff_g = 0,jl_g/jpan-1
         do ioff_g = 0,il_g/ipan-1
            j_g = 1
            do i_g = 1,ipan
               iqg = i_g + ioff_g*ipan + (j_g + joff_g*jpan - 1)*il_g
               disarray_g = real(x_g(iqq)*x_g(iqg) + y_g(iqq)*y_g(iqg) + z_g(iqq)*z_g(iqg))
               disarray_g = acos( max( min( disarray_g, 1. ), -1. ) )
               if ( disarray_g < maxdis ) then
                  iproc = qproc(iqg)
                  neigharray_g(iproc) = .true.
               end if
            end do
            j_g = jpan
            do i_g = 1,ipan
               iqg = i_g + ioff_g*ipan + (j_g + joff_g*jpan - 1)*il_g
               disarray_g = real(x_g(iqq)*x_g(iqg) + y_g(iqq)*y_g(iqg) + z_g(iqq)*z_g(iqg))
               disarray_g = acos( max( min( disarray_g, 1. ), -1. ) )
               if ( disarray_g < maxdis ) then
                  iproc = qproc(iqg)
                  neigharray_g(iproc) = .true.
               end if
            end do
            i_g = 1
            do j_g = 1,jpan
               iqg = i_g + ioff_g*ipan + (j_g + joff_g*jpan - 1)*il_g
               disarray_g = real(x_g(iqq)*x_g(iqg) + y_g(iqq)*y_g(iqg) + z_g(iqq)*z_g(iqg))
               disarray_g = acos( max( min( disarray_g, 1. ), -1. ) )
               if ( disarray_g < maxdis ) then
                  iproc = qproc(iqg)
                  neigharray_g(iproc) = .true.
               end if
            end do
            i_g = ipan
            do j_g = 1,jpan
               iqg = i_g + ioff_g*ipan + (j_g + joff_g*jpan - 1)*il_g
               disarray_g = real(x_g(iqq)*x_g(iqg) + y_g(iqq)*y_g(iqg) + z_g(iqq)*z_g(iqg))
               disarray_g = acos( max( min( disarray_g, 1. ), -1. ) )
               if ( disarray_g < maxdis ) then
                  iproc = qproc(iqg)
                  neigharray_g(iproc) = .true.
               end if
            end do
         end do
      end do   
            
   end subroutine checkdistance            
   
   subroutine check_set(ind,str,i,j,n,iq)
      integer, intent(in) :: ind,i,j,n,iq
      character(len=*) :: str
      if ( ind == huge(1) ) then
         write(6,*) str, " not set", myid, i, j, n, iq
         call ccmpi_abort(-1)
      end if
   end subroutine check_set

   subroutine bounds2(t, nrows, corner, nehalf)
      ! Copy the boundary regions
      real, dimension(ifull+iextra), intent(inout) :: t
      integer, intent(in), optional :: nrows
      logical, intent(in), optional :: corner
      logical, intent(in), optional :: nehalf
      logical :: extra, single, double
      integer :: iproc, send_len, recv_len
      integer :: rcount, jproc, mproc, iq
      integer, dimension(neighnum) :: rslen, sslen
      integer(kind=4) :: ierr, itag=1, llen, sreq, lproc
      integer(kind=4) :: ldone, lcomm
      integer(kind=4), dimension(neighnum) :: donelist
      integer(kind=4), dimension(MPI_STATUS_SIZE,size(ireq)) :: status
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif   

      double = .false.
      extra = .false.
      single = .true.
      if ( present(nrows) ) then
         if ( nrows == 2 ) then
            double = .true.
         end if
      end if
      if ( .not. double ) then
         if ( present(corner) ) then
            extra = corner
         end if
         if ( .not. extra ) then
            if ( present(nehalf) ) then
               single = .not. nehalf
            end if
         end if
      end if
      
      ! Split messages into corner and non-corner processors
      if ( double ) then
         rslen(:) = bnds(neighlist(:))%rlen2
         sslen(:) = bnds(neighlist(:))%slen2
      else if ( extra ) then
         rslen(:) = bnds(neighlist(:))%rlenx_fn(maxcolour)
         sslen(:) = bnds(neighlist(:))%slenx_fn(maxcolour)
      else if ( single ) then
         rslen(:) = bnds(neighlist(:))%rlen_fn(maxcolour)
         sslen(:) = bnds(neighlist(:))%slen_fn(maxcolour)
      else
         rslen(:) = bnds(neighlist(:))%rlenh_fn(maxcolour)
         sslen(:) = bnds(neighlist(:))%slenh_fn(maxcolour)
      end if

      ! Set up the buffers to recv
      lcomm = comm_world
      nreq = 0
      do iproc = 1,neighnum
         recv_len = rslen(iproc)
         if ( recv_len > 0 ) then
            lproc = neighlist(iproc)  ! Recv from
            nreq  = nreq + 1
            rlist(nreq) = iproc
            llen  = recv_len
            call MPI_IRecv( bnds(lproc)%rbuf, llen, ltype, lproc, &
                   itag, lcomm, ireq(nreq), ierr )
         end if
      end do
      rreq = nreq
      !     Set up the buffers to send
      do iproc = neighnum,1,-1
         ! Build up list of points
         send_len = sslen(iproc)
         if ( send_len > 0 ) then
            lproc = neighlist(iproc)  ! Send to
            do concurrent (iq = 1:send_len)
               bnds(lproc)%sbuf(iq) = t(bnds(lproc)%send_list(iq))
            end do   
            nreq  = nreq + 1
            llen  = send_len
            call MPI_ISend( bnds(lproc)%sbuf, llen, ltype, lproc, &
                 itag, lcomm, ireq(nreq), ierr )
         end if
      end do
      
      ! Unpack incomming messages
      rcount = rreq
      do while ( rcount > 0 )
         call START_LOG(mpiwait_begin)
         call MPI_Waitsome( rreq, ireq(1:rreq), ldone, donelist, status, ierr )
         call END_LOG(mpiwait_end)
         rcount = rcount - ldone
         do jproc = 1,ldone
            mproc = donelist(jproc)
            iproc = rlist(mproc)  ! Recv from
            lproc = neighlist(iproc)
            ! unpack_list(iq) is index into extended region
            do concurrent (iq = 1:rslen(iproc))
               t(ifull+bnds(lproc)%unpack_list(iq)) = bnds(lproc)%rbuf(iq)
            end do   
         end do
      end do

      ! Clear any remaining messages
      sreq = nreq - rreq
      if ( sreq > 0 ) then
         call START_LOG(mpiwait_begin)
         call MPI_Waitall( sreq, ireq(rreq+1:nreq), status, ierr)
         call END_LOG(mpiwait_end)
      end if   

   end subroutine bounds2

   subroutine bounds3(t, nrows, klim, corner, nehalf)
      ! Copy the boundary regions. Only this routine requires the extra klim
      ! argument (for helmsol).
      real, dimension(:,:), intent(inout) :: t
      integer, intent(in), optional :: nrows, klim
      logical, intent(in), optional :: corner
      logical, intent(in), optional :: nehalf
      logical :: extra, single, double
      integer :: iproc, kx, send_len, recv_len
      integer :: rcount, jproc, mproc, iq, k
      integer, dimension(neighnum) :: rslen, sslen
      integer(kind=4) :: ierr, itag=2, llen, sreq, lproc
      integer(kind=4) :: ldone, lcomm
      integer(kind=4), dimension(neighnum) :: donelist
      integer(kind=4), dimension(MPI_STATUS_SIZE,size(ireq)) :: status
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif  
      
      kx = size(t,2)
      double = .false.
      extra  = .false.
      single = .true.
      if ( present(klim) ) then
         kx = klim
      end if
      if ( present(nrows) ) then
         if ( nrows == 2 ) then
            double = .true.
         end if
      end if
      if ( .not. double ) then
         if ( present(corner) ) then
            extra = corner
         end if
         if ( .not. extra ) then
            if ( present(nehalf) ) then
               single = .not. nehalf
            end if
         end if
      end if

      ! Split messages into corner and non-corner processors
      if ( double ) then
         rslen = bnds(neighlist)%rlen2
         sslen = bnds(neighlist)%slen2
      else if ( extra ) then
         rslen = bnds(neighlist)%rlenx_fn(maxcolour)
         sslen = bnds(neighlist)%slenx_fn(maxcolour)
      else if ( single ) then
         rslen = bnds(neighlist)%rlen_fn(maxcolour)
         sslen = bnds(neighlist)%slen_fn(maxcolour)
      else
         rslen = bnds(neighlist)%rlenh_fn(maxcolour)
         sslen = bnds(neighlist)%slenh_fn(maxcolour)
      end if

      ! Set up the buffers to send and recv
      lcomm = comm_world
      nreq = 0
      do iproc = 1,neighnum
         recv_len = rslen(iproc)
         if ( recv_len > 0 ) then
            lproc = neighlist(iproc)  ! Recv from
            nreq = nreq + 1
            rlist(nreq) = iproc
            llen = recv_len*kx
            call MPI_IRecv( bnds(lproc)%rbuf, llen, ltype, lproc, &
                 itag, lcomm, ireq(nreq), ierr )
         end if
      end do
      rreq = nreq
      do iproc = neighnum,1,-1
         send_len = sslen(iproc)
         if ( send_len > 0 ) then
            lproc = neighlist(iproc)  ! Send to
            do k = 1, kx
               do concurrent (iq = 1:send_len) 
                  bnds(lproc)%sbuf(iq+(k-1)*send_len) = t(bnds(lproc)%send_list(iq),k)
               end do
            end do   
            nreq = nreq + 1
            llen = send_len*kx
            call MPI_ISend( bnds(lproc)%sbuf, llen, ltype, lproc, &
                 itag, lcomm, ireq(nreq), ierr )
         end if
      end do

      ! Unpack incomming messages
      rcount = rreq
      do while ( rcount > 0 )
         call START_LOG(mpiwait_begin)
         call MPI_Waitsome( rreq, ireq(1:rreq), ldone, donelist, status, ierr )
         call END_LOG(mpiwait_end)
         rcount = rcount - ldone
         do jproc = 1,ldone
            mproc = donelist(jproc)
            iproc = rlist(mproc)  ! Recv from
            lproc = neighlist(iproc)
            do k = 1,kx
               do concurrent (iq = 1:rslen(iproc))
                  t(ifull+bnds(lproc)%unpack_list(iq),k) &
                      = bnds(lproc)%rbuf(iq+(k-1)*rslen(iproc))
               end do
            end do   
         end do
      end do

      ! Clear any remaining messages
      sreq = nreq - rreq
      if ( sreq > 0 ) then
         call START_LOG(mpiwait_begin)
         call MPI_Waitall( sreq, ireq(rreq+1:nreq), status, ierr)
         call END_LOG(mpiwait_end)
      end if

   end subroutine bounds3

   subroutine bounds3r8(t, nrows, klim, corner, nehalf)
      ! Copy the boundary regions. Only this routine requires the extra klim
      ! argument (for helmsol).
      real(kind=8), dimension(:,:), intent(inout) :: t
      integer, intent(in), optional :: nrows, klim
      logical, intent(in), optional :: corner
      logical, intent(in), optional :: nehalf
      logical :: extra, single, double
      integer :: iproc, kx, send_len, recv_len
      integer :: rcount, jproc, mproc, iq, k
      integer, dimension(neighnum) :: rslen, sslen
      integer(kind=4) :: ierr, itag = 2, llen, sreq, lproc
      integer(kind=4) :: ldone, lcomm
      integer(kind=4), dimension(neighnum) :: donelist
      integer(kind=4), dimension(MPI_STATUS_SIZE,size(ireq)) :: status
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
      
      kx = size(t,2)
      double = .false.
      extra  = .false.
      single = .true.
      if ( present(klim) ) then
         kx = klim
      end if
      if ( present(nrows) ) then
         if ( nrows == 2 ) then
            double = .true.
         end if
      end if
      if ( .not. double ) then
         if ( present(corner) ) then
            extra = corner
         end if
         if ( .not. extra ) then
            if ( present(nehalf) ) then
               single = .not. nehalf
            end if
         end if
      end if

      ! Split messages into corner and non-corner processors
      if ( double ) then
         rslen = bnds(neighlist)%rlen2
         sslen = bnds(neighlist)%slen2
      else if ( extra ) then
         rslen = bnds(neighlist)%rlenx_fn(maxcolour)
         sslen = bnds(neighlist)%slenx_fn(maxcolour)
      else if ( single ) then
         rslen = bnds(neighlist)%rlen_fn(maxcolour)
         sslen = bnds(neighlist)%slen_fn(maxcolour)
      else
         rslen = bnds(neighlist)%rlenh_fn(maxcolour)
         sslen = bnds(neighlist)%slenh_fn(maxcolour)
      end if

      ! Set up the buffers to send and recv
      lcomm = comm_world
      nreq = 0
      do iproc = 1,neighnum
         recv_len = rslen(iproc)
         if ( recv_len > 0 ) then
            lproc = neighlist(iproc)  ! Recv from
            nreq = nreq + 1
            rlist(nreq) = iproc
            llen = recv_len*kx
            call MPI_IRecv( bnds(lproc)%r8buf, llen, ltype, lproc, &
                 itag, lcomm, ireq(nreq), ierr )
         end if
      end do
      rreq = nreq
      do iproc = neighnum,1,-1
         send_len = sslen(iproc)
         if ( send_len > 0 ) then
            lproc = neighlist(iproc)  ! Send to
            do k = 1,kx
               do concurrent (iq = 1:send_len)
                  bnds(lproc)%s8buf(iq+(k-1)*send_len) = t(bnds(lproc)%send_list(iq),k)
               end do
            end do   
            nreq = nreq + 1
            llen = send_len*kx
            call MPI_ISend( bnds(lproc)%s8buf, llen, ltype, lproc, &
                 itag, lcomm, ireq(nreq), ierr )
         end if
      end do
      
      ! Unpack incomming messages
      rcount = rreq
      do while ( rcount>0 )
         call START_LOG(mpiwait_begin)
         call MPI_Waitsome( rreq, ireq(1:rreq), ldone, donelist, status, ierr )
         call END_LOG(mpiwait_end)
         rcount = rcount - ldone
         do jproc = 1,ldone
            mproc = donelist(jproc)
            iproc = rlist(mproc)  ! Recv from
            lproc = neighlist(iproc)
            do k = 1,kx
               do concurrent (iq = 1:rslen(iproc))
                  t(ifull+bnds(lproc)%unpack_list(iq),k) &
                      = bnds(lproc)%r8buf(iq+(k-1)*rslen(iproc))
               end do
            end do   
         end do
      end do

      ! Clear any remaining messages
      sreq = nreq - rreq
      if ( sreq > 0 ) then
         call START_LOG(mpiwait_begin)
         call MPI_Waitall( sreq, ireq(rreq+1:nreq), status, ierr )
         call END_LOG(mpiwait_end)
      end if   

   end subroutine bounds3r8
   
   subroutine bounds4(t, nrows, corner, nehalf)
      ! Copy the boundary regions.
      real, dimension(:,:,:), intent(inout) :: t
      integer, intent(in), optional :: nrows
      logical, intent(in), optional :: corner
      logical, intent(in), optional :: nehalf
      logical :: extra, single, double
      integer :: iproc, kx, send_len, recv_len
      integer :: rcount, jproc, mproc, ntr, iq, k, l
      integer :: nstart, nend, ntot
      integer, dimension(neighnum) :: rslen, sslen
      integer(kind=4) :: ierr, itag=3, llen, sreq, lproc
      integer(kind=4) :: ldone, lcomm
      integer(kind=4), dimension(neighnum) :: donelist
      integer(kind=4), dimension(MPI_STATUS_SIZE,size(ireq)) :: status
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif  
      
      kx = size(t,2)
      ntr = size(t,3)
      double = .false.
      extra  = .false.
      single = .true.
      if ( present(nrows) ) then
         if ( nrows == 2 ) then
            double = .true.
         end if
      end if
      if ( .not. double ) then
         if ( present(corner) ) then
            extra = corner
         end if
         if ( .not. extra ) then
            if ( present(nehalf) ) then
               single = .not. nehalf
            end if
         end if
      end if

      ! Split messages into corner and non-corner processors
      if ( double ) then
         rslen = bnds(neighlist)%rlen2
         sslen = bnds(neighlist)%slen2
      else if ( extra ) then
         rslen = bnds(neighlist)%rlenx_fn(maxcolour)
         sslen = bnds(neighlist)%slenx_fn(maxcolour)
      else if ( single ) then
         rslen = bnds(neighlist)%rlen_fn(maxcolour)
         sslen = bnds(neighlist)%slen_fn(maxcolour)
      else
         rslen = bnds(neighlist)%rlenh_fn(maxcolour)
         sslen = bnds(neighlist)%slenh_fn(maxcolour)
      end if

      lcomm = comm_world
      
      do nstart = 1,ntr,nagg
         nend = min(nstart+nagg-1,ntr)
         ntot = nend - nstart + 1
      
         ! Set up the buffers to send
         nreq = 0
         do iproc = 1,neighnum
            recv_len = rslen(iproc)
            if ( recv_len > 0 ) then
               lproc = neighlist(iproc)  ! Recv from
               nreq = nreq + 1
               rlist(nreq) = iproc
               llen = recv_len*kx*ntot
               call MPI_IRecv( bnds(lproc)%rbuf, llen, ltype, lproc, &
                    itag, lcomm, ireq(nreq), ierr )
            end if
         end do
         rreq = nreq
         do iproc = neighnum,1,-1
            send_len = sslen(iproc)
            if ( send_len > 0 ) then
               lproc = neighlist(iproc)  ! Send to
               do l = 1,ntot
                  do k = 1,kx 
                     do concurrent (iq = 1:send_len)
                        bnds(lproc)%sbuf(iq+(k-1)*send_len+(l-1)*send_len*kx) = &
                          t(bnds(lproc)%send_list(iq),k,l+nstart-1)
                     end do
                  end do
               end do   
               nreq = nreq + 1
               llen = send_len*kx*ntot
               call MPI_ISend( bnds(lproc)%sbuf, llen, ltype, lproc, &
                    itag, lcomm, ireq(nreq), ierr )
            end if
         end do

         ! Unpack incomming messages
         rcount = rreq
         do while ( rcount > 0 )
            call START_LOG(mpiwait_begin)
            call MPI_Waitsome( rreq, ireq(1:rreq), ldone, donelist, status, ierr )
            call END_LOG(mpiwait_end)
            rcount = rcount - ldone
            do jproc = 1,ldone
               mproc = donelist(jproc)
               iproc = rlist(mproc)  ! Recv from
               lproc = neighlist(iproc)
               do l = 1,ntot
                  do k = 1,kx 
                     do concurrent (iq = 1:rslen(iproc))  
                        t(ifull+bnds(lproc)%unpack_list(iq),k,l+nstart-1)              &
                             = bnds(lproc)%rbuf(iq+(k-1)*rslen(iproc)+(l-1)*rslen(iproc)*kx)
                     end do
                  end do
               end do   
            end do
         end do

         ! Clear any remaining messages
         sreq = nreq - rreq
         if ( sreq > 0 ) then
            call START_LOG(mpiwait_begin)
            call MPI_Waitall( sreq, ireq(rreq+1:nreq), status, ierr )
            call END_LOG(mpiwait_end)
         end if

      end do   

   end subroutine bounds4

   subroutine bounds_colour_send(t, lcolour, klim, corner)
      ! Copy the boundary regions. This version allows supports updating
      ! different gridpoint colours
      real, dimension(:,:), intent(in) :: t
      integer, intent(in) :: lcolour
      integer, intent(in), optional :: klim
      logical, intent(in), optional :: corner
      logical :: extra
      integer :: iproc, kx, recv_len, iqq, ibeg, iend, iq, k
      integer(kind=4) :: ierr, itag=4, llen, lproc, lcomm
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif 
      
      if ( lcolour<1 .or. lcolour>maxcolour ) then
         write(6,*) "ERROR: Invalid colour for bounds_colour_send"
         call ccmpi_abort(-1)
      end if   
      
      kx = size(t, 2)
      extra = .false.
      if ( present(klim) ) then
         kx = klim
      end if
      if ( present(corner) ) then
         extra = corner
      end if

!     Set up the buffers to send and recv
      lcomm = comm_world
      nreq = 0
      do iproc = 1,neighnum
         lproc = neighlist(iproc)  ! Recv from
         recv_len = (bnds(lproc)%rlenh_fn(lcolour)-bnds(lproc)%rlenh_bg(lcolour)+1   &
                    +bnds(lproc)%rlen_fn(lcolour)-bnds(lproc)%rlen_bg(lcolour)+1)*kx
         if ( extra ) then
            recv_len = recv_len + (bnds(lproc)%rlenx_fn(lcolour)-bnds(lproc)%rlenx_bg(lcolour)+1)*kx 
         end if    
         if ( recv_len > 0 ) then
            nreq = nreq + 1
            rlist(nreq) = iproc
            llen = recv_len
            call MPI_IRecv( bnds(lproc)%rbuf, llen, ltype, lproc, &
                            itag, lcomm, ireq(nreq), ierr )
         end if
      end do
      rreq = nreq
      do iproc = neighnum,1,-1
         lproc = neighlist(iproc)  ! Send to
         iqq = 0
         ibeg = bnds(lproc)%slenh_bg(lcolour)
         iend = bnds(lproc)%slenh_fn(lcolour)
         if ( iend >= ibeg ) then
            do k = 1,kx 
               do concurrent (iq = 1:iend-ibeg+1)    
                  bnds(lproc)%sbuf(iqq+iq+(k-1)*(iend-ibeg+1))  &
                      = t(bnds(lproc)%send_list(iq+ibeg-1),k)
               end do
            end do   
         end if
         iqq = iqq + (iend-ibeg+1)*kx
         ibeg = bnds(lproc)%slen_bg(lcolour)
         iend = bnds(lproc)%slen_fn(lcolour)
         if ( iend >= ibeg ) then
            do k = 1,kx 
               do concurrent (iq = 1:iend-ibeg+1)    
                  bnds(lproc)%sbuf(iqq+iq+(k-1)*(iend-ibeg+1))  &
                      = t(bnds(lproc)%send_list(iq+ibeg-1),k)
               end do
            end do   
         end if   
         iqq = iqq + (iend-ibeg+1)*kx
         if ( extra ) then
            ibeg = bnds(lproc)%slenx_bg(lcolour)
            iend = bnds(lproc)%slenx_fn(lcolour)
            if ( iend >= ibeg ) then
               do k = 1,kx 
                  do concurrent (iq = 1:iend-ibeg+1)    
                     bnds(lproc)%sbuf(iqq+iq+(k-1)*(iend-ibeg+1))  &
                         = t(bnds(lproc)%send_list(iq+ibeg-1),k)
                  end do
               end do   
            end if   
            iqq = iqq + (iend-ibeg+1)*kx
         end if  
         if ( iqq > 0 ) then
            nreq = nreq + 1
            llen = iqq
            call MPI_ISend( bnds(lproc)%sbuf, llen, ltype, lproc, &
                            itag, lcomm, ireq(nreq), ierr )
         end if
      end do

   end subroutine bounds_colour_send
   
   subroutine bounds_colour_recv(t, lcolour, klim, corner)
      ! Copy the boundary regions. This version allows supports updating
      ! different gridpoint colours
      real, dimension(:,:), intent(inout) :: t
      integer, intent(in) :: lcolour
      integer, intent(in), optional :: klim
      logical, intent(in), optional :: corner
      logical :: extra
      integer :: iproc, kx, iqq, ibeg, iend
      integer :: rcount, jproc, iq, k
      integer(kind=4) :: ierr, sreq, lproc, ldone
      integer(kind=4), dimension(neighnum) :: donelist
      integer(kind=4), dimension(MPI_STATUS_SIZE,size(ireq)) :: status

      if ( lcolour<1 .or. lcolour>maxcolour ) then
         write(6,*) "ERROR: Invalid colour for bounds_colour_recv"
         call ccmpi_abort(-1)
      end if   
      
      kx = size(t, 2)
      extra = .false.
      if ( present(klim) ) then
         kx = klim
      end if
      if ( present(corner) ) then
         extra = corner
      end if
      
      ! Unpack incomming messages
      rcount = rreq
      do while ( rcount > 0 )
         call START_LOG(mpiwait_begin)
         call MPI_Waitsome( rreq, ireq, ldone, donelist, status, ierr )
         call END_LOG(mpiwait_end)
         rcount = rcount - ldone
         do jproc = 1,ldone
            iproc = rlist(donelist(jproc))  ! Recv from
            lproc = neighlist(iproc)
            iqq = 0
            ibeg = bnds(lproc)%rlenh_bg(lcolour)
            iend = bnds(lproc)%rlenh_fn(lcolour)
            do k = 1,kx
               do concurrent (iq = 1:iend-ibeg+1)
                  t(ifull+bnds(lproc)%unpack_list(iq+ibeg-1),k)  &
                      = bnds(lproc)%rbuf(iqq+iq+(k-1)*(iend-ibeg+1))
               end do
            end do   
            iqq = iqq + (iend-ibeg+1)*kx
            ibeg = bnds(lproc)%rlen_bg(lcolour)
            iend = bnds(lproc)%rlen_fn(lcolour)
            do k = 1,kx
               do concurrent (iq = 1:iend-ibeg+1)
                  t(ifull+bnds(lproc)%unpack_list(iq+ibeg-1),k)  &
                      = bnds(lproc)%rbuf(iqq+iq+(k-1)*(iend-ibeg+1))
               end do
            end do   
            iqq = iqq + (iend-ibeg+1)*kx
            if ( extra ) then
               ibeg = bnds(lproc)%rlenx_bg(lcolour)
               iend = bnds(lproc)%rlenx_fn(lcolour)
               do k = 1,kx
                  do concurrent (iq = 1:iend-ibeg+1)
                     t(ifull+bnds(lproc)%unpack_list(iq+ibeg-1),k)  &
                         = bnds(lproc)%rbuf(iqq+iq+(k-1)*(iend-ibeg+1))
                  end do
               end do
               iqq = iqq + (iend-ibeg+1)*kx
            end if   
         end do
      end do

      ! Clear any remaining messages
      sreq = nreq - rreq
      if ( sreq > 0 ) then
         call START_LOG(mpiwait_begin)
         call MPI_Waitall( sreq, ireq(rreq+1:nreq), status, ierr )
         call END_LOG(mpiwait_end)
      end if   

   end subroutine bounds_colour_recv
   
   subroutine boundsuv2(u, v, nrows, stag, allvec)
      ! Copy the boundary regions of u and v. This doesn't require the
      ! diagonal points like (0,0), but does have to take care of the
      ! direction changes.
      real, dimension(ifull+iextra), intent(inout) :: u, v
      integer, intent(in), optional :: nrows, stag
      logical, intent(in), optional :: allvec
      logical :: double, extra
      logical :: fsvwu, fnveu, fssvwwu, fnnveeu
      logical :: fsuev, fnnueev
      integer :: iq, iqz, iproc, iqq, recv_len
      integer :: rcount, myrlen, jproc, mproc, stagmode
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif   
      integer(kind=4) :: ierr, itag=5, llen, sreq, lproc
      integer(kind=4) :: ldone, lcomm
      integer(kind=4), dimension(neighnum) :: donelist
      integer(kind=4), dimension(MPI_STATUS_SIZE,size(ireq)) :: status

      double = .false.
      extra = .false.
      stagmode = 0
      if ( present(nrows)) then
         if ( nrows == 2 ) then
            double = .true.
         end if
      end if
      if ( present(stag) ) then
         stagmode = stag
      end if
      if ( present(allvec) ) then
         extra = allvec
      end if
      
      if ( extra .and. double ) then
         fsvwu = .true.
         fnveu = .true.
         fssvwwu = .true.
         fnnveeu = .true.
         fsuev = .true.
         fnnueev = .true.     
      else if ( extra ) then
         fsvwu = .true.
         fnveu = .true.
         fssvwwu = .false.
         fnnveeu = .false.
         fsuev = .true.
         fnnueev = .false.
      else if ( double ) then
         fsvwu = .true.
         fnveu = .true.
         fssvwwu = .true.
         fnnveeu = .true.
         fsuev = .false.
         fnnueev = .false.
      else if ( stagmode == 1 ) then
         fsvwu = .false.
         fnveu = .true.
         fssvwwu = .false.
         fnnveeu = .true.
         fsuev = .false.
         fnnueev = .false.
      else if ( stagmode == 2 ) then
         fsvwu = .true.
         fnveu = .true.
         fssvwwu = .false.
         fnnveeu = .true. ! fnnveeu requires fnveu
         fsuev = .false.
         fnnueev = .false.
      else if ( stagmode == 3 ) then
         fsvwu = .true.
         fnveu = .true.
         fssvwwu = .true. ! fssvwwu requires fsvwu
         fnnveeu = .false.
         fsuev = .false.
         fnnueev = .false.
      else if ( stagmode == 5 ) then
         fsvwu = .true.
         fnveu = .false.
         fssvwwu = .true.
         fnnveeu = .false.
         fsuev = .false.
         fnnueev = .false.
      else if ( stagmode == -9 ) then
         fsvwu = .true.
         fnveu = .false.
         fssvwwu = .false.
         fnnveeu = .false.
         fsuev = .false.
         fnnueev = .false.
      else if ( stagmode == -10 ) then
         fsvwu = .false.
         fnveu = .true.
         fssvwwu = .false.
         fnnveeu = .false.
         fsuev = .false.
         fnnueev = .false.
      else
         fsvwu = .true.
         fnveu = .true.
         fssvwwu = .false.
         fnnveeu = .false.
         fsuev = .false.
         fnnueev = .false.
      end if
      myrlen = bnds(myid)%rlen_eev_fn

!     Set up the buffers to send
      lcomm = comm_world
      nreq = 0
      do iproc = 1,neighnum
         lproc = neighlist(iproc)  ! Recv from
         recv_len = 0
         if ( fsvwu ) then
            recv_len = recv_len + bnds(lproc)%rlen_wu_fn - bnds(lproc)%rlen_sv_bg + 1
         end if
         if ( fnveu ) then
            recv_len = recv_len + bnds(lproc)%rlen_eu_fn - bnds(lproc)%rlen_nv_bg + 1
         end if         
         if ( fssvwwu ) then
            recv_len = recv_len + bnds(lproc)%rlen_wwu_fn - bnds(lproc)%rlen_ssv_bg + 1
         end if         
         if ( fnnveeu ) then
            recv_len = recv_len + bnds(lproc)%rlen_eeu_fn - bnds(lproc)%rlen_nnv_bg + 1
         end if
         if ( fsuev ) then
            recv_len = recv_len + bnds(lproc)%rlen_ev_fn - bnds(lproc)%rlen_su_bg + 1 
         end if
         if ( fnnueev ) then
            recv_len = recv_len + bnds(lproc)%rlen_eev_fn - bnds(lproc)%rlen_nnu_bg + 1 
         end if
         if ( recv_len > 0 ) then 
            nreq = nreq + 1
            rlist(nreq) = iproc
            llen = recv_len
            call MPI_IRecv( bnds(lproc)%rbuf, llen, ltype, lproc, &
                 itag, lcomm, ireq(nreq), ierr )
         end if
      end do
      rreq = nreq
      do iproc = neighnum,1,-1
         lproc = neighlist(iproc)  ! Send to
         ! Build up list of points
         iqq = 0
         if ( fsvwu ) then
            iqz = iqq - bnds(lproc)%slen_sv_bg + 1 
            do concurrent (iq = bnds(lproc)%slen_sv_bg:bnds(lproc)%slen_wu_fn)
               ! Use abs because sign is used as u/v flag
               if ( bnds(lproc)%send_list_uv(iq)>0 .neqv. bnds(lproc)%send_swap(iq) ) then
                  bnds(lproc)%sbuf(iqz+iq) = bnds(lproc)%send_neg(iq)*u(abs(bnds(lproc)%send_list_uv(iq)))
               else
                  bnds(lproc)%sbuf(iqz+iq) = bnds(lproc)%send_neg(iq)*v(abs(bnds(lproc)%send_list_uv(iq)))
               end if 
            end do
            iqq = iqq + bnds(lproc)%slen_wu_fn - bnds(lproc)%slen_sv_bg + 1
         end if
         if ( fnveu ) then
            iqz = iqq - bnds(lproc)%slen_nv_bg + 1 
            do concurrent (iq = bnds(lproc)%slen_nv_bg:bnds(lproc)%slen_eu_fn)
               ! Use abs because sign is used as u/v flag
               if ( bnds(lproc)%send_list_uv(iq)>0 .neqv. bnds(lproc)%send_swap(iq) ) then
                  bnds(lproc)%sbuf(iqz+iq) = bnds(lproc)%send_neg(iq)*u(abs(bnds(lproc)%send_list_uv(iq)))
               else
                  bnds(lproc)%sbuf(iqz+iq) = bnds(lproc)%send_neg(iq)*v(abs(bnds(lproc)%send_list_uv(iq)))
               end if 
            end do
            iqq = iqq + bnds(lproc)%slen_eu_fn - bnds(lproc)%slen_nv_bg + 1
         end if
         if ( fssvwwu ) then
            iqz = iqq - bnds(lproc)%slen_ssv_bg + 1 
            do concurrent (iq = bnds(lproc)%slen_ssv_bg:bnds(lproc)%slen_wwu_fn)
               ! Use abs because sign is used as u/v flag
               if ( bnds(lproc)%send_list_uv(iq)>0 .neqv. bnds(lproc)%send_swap(iq) ) then
                  bnds(lproc)%sbuf(iqz+iq) = bnds(lproc)%send_neg(iq)*u(abs(bnds(lproc)%send_list_uv(iq)))
               else
                  bnds(lproc)%sbuf(iqz+iq) = bnds(lproc)%send_neg(iq)*v(abs(bnds(lproc)%send_list_uv(iq)))
               end if 
            end do
            iqq = iqq + bnds(lproc)%slen_wwu_fn - bnds(lproc)%slen_ssv_bg + 1
         end if
         if ( fnnveeu ) then
            iqz = iqq - bnds(lproc)%slen_nnv_bg + 1 
            do concurrent (iq = bnds(lproc)%slen_nnv_bg:bnds(lproc)%slen_eeu_fn)
               ! Use abs because sign is used as u/v flag
               if ( bnds(lproc)%send_list_uv(iq)>0 .neqv. bnds(lproc)%send_swap(iq) ) then
                  bnds(lproc)%sbuf(iqz+iq) = bnds(lproc)%send_neg(iq)*u(abs(bnds(lproc)%send_list_uv(iq)))
               else
                  bnds(lproc)%sbuf(iqz+iq) = bnds(lproc)%send_neg(iq)*v(abs(bnds(lproc)%send_list_uv(iq)))
               end if 
            end do
            iqq = iqq + bnds(lproc)%slen_eeu_fn - bnds(lproc)%slen_nnv_bg + 1
         end if
         if ( fsuev ) then
            iqz = iqq - bnds(lproc)%slen_su_bg + 1 
            do concurrent (iq = bnds(lproc)%slen_su_bg:bnds(lproc)%slen_ev_fn)
               ! Use abs because sign is used as u/v flag
               if ( bnds(lproc)%send_list_uv(iq)>0 .neqv. bnds(lproc)%send_swap(iq) ) then
                  bnds(lproc)%sbuf(iqz+iq) = bnds(lproc)%send_neg(iq)*u(abs(bnds(lproc)%send_list_uv(iq)))
               else
                  bnds(lproc)%sbuf(iqz+iq) = bnds(lproc)%send_neg(iq)*v(abs(bnds(lproc)%send_list_uv(iq)))
               end if 
            end do
            iqq = iqq + bnds(lproc)%slen_ev_fn - bnds(lproc)%slen_su_bg + 1             
         end if
         if ( fnnueev ) then
            iqz = iqq - bnds(lproc)%slen_nnu_bg + 1 
            do concurrent (iq = bnds(lproc)%slen_nnu_bg:bnds(lproc)%slen_eev_fn)
               ! Use abs because sign is used as u/v flag
               if ( bnds(lproc)%send_list_uv(iq)>0 .neqv. bnds(lproc)%send_swap(iq) ) then
                  bnds(lproc)%sbuf(iqz+iq) = bnds(lproc)%send_neg(iq)*u(abs(bnds(lproc)%send_list_uv(iq)))
               else
                  bnds(lproc)%sbuf(iqz+iq) = bnds(lproc)%send_neg(iq)*v(abs(bnds(lproc)%send_list_uv(iq)))
               end if 
            end do
            iqq = iqq + bnds(lproc)%slen_eev_fn - bnds(lproc)%slen_nnu_bg + 1
         end if
         if ( iqq > 0 ) then
            nreq = nreq + 1
            llen = iqq
            call MPI_ISend( bnds(lproc)%sbuf, llen, ltype, lproc, &
                 itag, lcomm, ireq(nreq), ierr )
         end if
      end do

      ! See if there are any points on my own processor that need
      ! to be fixed up. This will only be in the case when nproc < npanels.
      do concurrent (iq = 1:myrlen)
         ! request_list is same as send_list in this case
         ! unpack_list(iq) is index into extended region 
         if ( bnds(myid)%unpack_list_uv(iq) > 0 .and. (bnds(myid)%request_list_uv(iq)>0 .neqv. bnds(myid)%uv_swap(iq)) ) then
            u(ifull+bnds(myid)%unpack_list_uv(iq)) = bnds(myid)%uv_neg(iq)*u(abs(bnds(myid)%request_list_uv(iq)))
         else if ( bnds(myid)%request_list_uv(iq)>0 .neqv. bnds(myid)%uv_swap(iq) ) then
            v(ifull-bnds(myid)%unpack_list_uv(iq)) = bnds(myid)%uv_neg(iq)*u(abs(bnds(myid)%request_list_uv(iq)))
         else if ( bnds(myid)%unpack_list_uv(iq) > 0 ) then
            u(ifull+bnds(myid)%unpack_list_uv(iq)) = bnds(myid)%uv_neg(iq)*v(abs(bnds(myid)%request_list_uv(iq)))
         else
            v(ifull-bnds(myid)%unpack_list_uv(iq)) = bnds(myid)%uv_neg(iq)*v(abs(bnds(myid)%request_list_uv(iq)))
         end if
      end do
      
      ! Unpack incomming messages
      rcount = rreq
      do while ( rcount > 0 )
         call START_LOG(mpiwait_begin)
         call MPI_Waitsome( rreq, ireq, ldone, donelist, status, ierr )
         call END_LOG(mpiwait_end)
         rcount = rcount - ldone
         do jproc = 1,ldone
            mproc = donelist(jproc)
            iproc = rlist(mproc)  ! Recv from
            lproc = neighlist(iproc)
            iqq = 0
            if ( fsvwu ) then
               iqz = iqq - bnds(lproc)%rlen_sv_bg + 1 
               do concurrent (iq = bnds(lproc)%rlen_sv_bg:bnds(lproc)%rlen_wu_fn)
                  ! unpack_list(iq) is index into extended region
                  if ( bnds(lproc)%unpack_list_uv(iq) > 0 ) then
                     u(ifull+bnds(lproc)%unpack_list_uv(iq)) = bnds(lproc)%rbuf(iqz+iq)
                  else
                     v(ifull-bnds(lproc)%unpack_list_uv(iq)) = bnds(lproc)%rbuf(iqz+iq)
                  end if
               end do
               iqq = iqq + bnds(lproc)%rlen_wu_fn - bnds(lproc)%rlen_sv_bg + 1
            end if
            if ( fnveu ) then
               iqz = iqq - bnds(lproc)%rlen_nv_bg + 1 
               do concurrent (iq = bnds(lproc)%rlen_nv_bg:bnds(lproc)%rlen_eu_fn)
                  ! unpack_list(iq) is index into extended region
                  if ( bnds(lproc)%unpack_list_uv(iq) > 0 ) then
                     u(ifull+bnds(lproc)%unpack_list_uv(iq)) = bnds(lproc)%rbuf(iqz+iq)
                  else
                     v(ifull-bnds(lproc)%unpack_list_uv(iq)) = bnds(lproc)%rbuf(iqz+iq)
                  end if
               end do
               iqq = iqq + bnds(lproc)%rlen_eu_fn - bnds(lproc)%rlen_nv_bg + 1
            end if
            if ( fssvwwu ) then
               iqz = iqq - bnds(lproc)%rlen_ssv_bg + 1 
               do concurrent (iq = bnds(lproc)%rlen_ssv_bg:bnds(lproc)%rlen_wwu_fn)
                  ! unpack_list(iq) is index into extended region
                  if ( bnds(lproc)%unpack_list_uv(iq) > 0 ) then
                     u(ifull+bnds(lproc)%unpack_list_uv(iq)) = bnds(lproc)%rbuf(iqz+iq)
                  else
                     v(ifull-bnds(lproc)%unpack_list_uv(iq)) = bnds(lproc)%rbuf(iqz+iq)
                  end if
               end do
               iqq = iqq + bnds(lproc)%rlen_wwu_fn - bnds(lproc)%rlen_ssv_bg + 1
            end if
            if ( fnnveeu ) then
               iqz = iqq - bnds(lproc)%rlen_nnv_bg + 1 
               do concurrent (iq = bnds(lproc)%rlen_nnv_bg:bnds(lproc)%rlen_eeu_fn)
                  ! unpack_list(iq) is index into extended region
                  if ( bnds(lproc)%unpack_list_uv(iq) > 0 ) then
                     u(ifull+bnds(lproc)%unpack_list_uv(iq)) = bnds(lproc)%rbuf(iqz+iq)
                  else
                     v(ifull-bnds(lproc)%unpack_list_uv(iq)) = bnds(lproc)%rbuf(iqz+iq)
                  end if
               end do
               iqq = iqq + bnds(lproc)%rlen_eeu_fn - bnds(lproc)%rlen_nnv_bg + 1
            end if
            if ( fsuev ) then
               iqz = iqq - bnds(lproc)%rlen_su_bg + 1 
               do concurrent (iq = bnds(lproc)%rlen_su_bg:bnds(lproc)%rlen_ev_fn)
                  ! unpack_list(iq) is index into extended region
                  if ( bnds(lproc)%unpack_list_uv(iq) > 0 ) then
                     u(ifull+bnds(lproc)%unpack_list_uv(iq)) = bnds(lproc)%rbuf(iqz+iq)
                  else
                     v(ifull-bnds(lproc)%unpack_list_uv(iq)) = bnds(lproc)%rbuf(iqz+iq)
                  end if
               end do
               iqq = iqq + bnds(lproc)%rlen_ev_fn - bnds(lproc)%rlen_su_bg + 1
            end if
            if ( fnnueev ) then
               iqz = iqq - bnds(lproc)%rlen_nnu_bg + 1 
               do concurrent (iq = bnds(lproc)%rlen_nnu_bg:bnds(lproc)%rlen_eev_fn)
                  ! unpack_list(iq) is index into extended region
                  if ( bnds(lproc)%unpack_list_uv(iq) > 0 ) then
                     u(ifull+bnds(lproc)%unpack_list_uv(iq)) = bnds(lproc)%rbuf(iqz+iq)
                  else
                     v(ifull-bnds(lproc)%unpack_list_uv(iq)) = bnds(lproc)%rbuf(iqz+iq)
                  end if
               end do
               iqq = iqq + bnds(lproc)%rlen_eev_fn - bnds(lproc)%rlen_nnu_bg + 1 
            end if
         end do
      end do

      ! Clear any remaining messages
      sreq = nreq - rreq
      if ( sreq > 0 ) then
         call START_LOG(mpiwait_begin)
         call MPI_Waitall( sreq, ireq(rreq+1:nreq), status, ierr )
         call END_LOG(mpiwait_end)
      end if   

   end subroutine boundsuv2

   subroutine boundsuv3(u, v, nrows, stag, allvec)
      ! Copy the boundary regions of u and v. This doesn't require the
      ! diagonal points like (0,0), but does have to take care of the
      ! direction changes.
      real, dimension(:,:), intent(inout) :: u, v
      integer, intent(in), optional :: nrows
      integer, intent(in), optional :: stag
      logical, intent(in), optional :: allvec
      logical :: double, extra
      logical :: fsvwu, fnveu, fssvwwu, fnnveeu
      logical :: fsuev, fnnueev
      integer :: iq, iqz, iproc, kx, rproc, sproc, iqq, recv_len
      integer :: rcount, myrlen, jproc, mproc, stagmode, k
      integer(kind=4) :: ierr, itag=6, llen, sreq, lproc
      integer(kind=4) :: ldone, lcomm
      integer(kind=4), dimension(neighnum) :: donelist
      integer(kind=4), dimension(MPI_STATUS_SIZE,size(ireq)) :: status
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif   
      real, dimension(maxbuflen*maxvertlen,neighnum) :: sbuf, rbuf

      kx = size(u, 2)
      double = .false.
      extra = .false.
      stagmode = 0
      if ( present(nrows) ) then
         if ( nrows == 2 ) then
            double = .true.
         end if
      end if
      if ( present(stag) ) then
         stagmode = stag
      end if
      if ( present(allvec) ) then
         extra = allvec
      end if

      if ( extra .and. double ) then
         fsvwu = .true.
         fnveu = .true.
         fssvwwu = .true.
         fnnveeu = .true.
         fsuev = .true.
         fnnueev = .true.
      else if ( extra ) then
         fsvwu = .true.
         fnveu = .true.
         fssvwwu = .false.
         fnnveeu = .false.
         fsuev = .true.
         fnnueev = .false.
      else if ( double ) then
         fsvwu = .true.
         fnveu = .true.
         fssvwwu = .true.
         fnnveeu = .true.
         fsuev = .false.
         fnnueev = .false.
      else if ( stagmode == 1 ) then
         fsvwu = .false.
         fnveu = .true.
         fssvwwu = .false.
         fnnveeu = .true.
         fsuev = .false.
         fnnueev = .false.
      else if ( stagmode == 2 ) then
         fsvwu = .true.
         fnveu = .true.
         fssvwwu = .false.
         fnnveeu = .true. ! fnnveeu requires fnveu
         fsuev = .false.
         fnnueev = .false.
      else if ( stagmode == 3 ) then
         fsvwu = .true.
         fnveu = .true.
         fssvwwu = .true. ! fssvwwu requires fsvwu
         fnnveeu = .false.
         fsuev = .false.
         fnnueev = .false.
      else if ( stagmode == 5 ) then
         fsvwu = .true.
         fnveu = .false.
         fssvwwu = .true.
         fnnveeu = .false.
         fsuev = .false.
         fnnueev = .false.
      else if ( stagmode == -9 ) then
         fsvwu = .true.
         fnveu = .false.
         fssvwwu = .false.
         fnnveeu = .false.
         fsuev = .false.
         fnnueev = .false.
      else if ( stagmode == -10 ) then
         fsvwu = .false.
         fnveu = .true.
         fssvwwu = .false.
         fnnveeu = .false.
         fsuev = .false.
         fnnueev = .false.
      else
         fsvwu = .true.
         fnveu = .true.
         fssvwwu = .false.
         fnnveeu = .false.
         fsuev = .false.
         fnnueev = .false.
      end if
      myrlen = bnds(myid)%rlen_eev_fn

!     Set up the buffers to send and recv
      lcomm = comm_world
      nreq = 0
      do iproc = 1,neighnum
         rproc = neighlist(iproc)  ! Recv from
         recv_len = 0
         if ( fsvwu ) then
            recv_len = recv_len + (bnds(rproc)%rlen_wu_fn - bnds(rproc)%rlen_sv_bg + 1)*kx
         end if
         if ( fnveu ) then
            recv_len = recv_len + (bnds(rproc)%rlen_eu_fn - bnds(rproc)%rlen_nv_bg + 1)*kx
         end if         
         if ( fssvwwu ) then
            recv_len = recv_len + (bnds(rproc)%rlen_wwu_fn - bnds(rproc)%rlen_ssv_bg + 1)*kx
         end if         
         if ( fnnveeu ) then
            recv_len = recv_len + (bnds(rproc)%rlen_eeu_fn - bnds(rproc)%rlen_nnv_bg + 1)*kx
         end if
         if ( fsuev ) then
            recv_len = recv_len + (bnds(rproc)%rlen_ev_fn - bnds(rproc)%rlen_su_bg + 1)*kx
         end if
         if ( fnnueev ) then
            recv_len = recv_len + (bnds(rproc)%rlen_eev_fn - bnds(rproc)%rlen_nnu_bg + 1)*kx
         end if
         if ( recv_len > 0 ) then 
            nreq = nreq + 1
            rlist(nreq) = iproc
            llen = recv_len
            lproc = rproc
            call MPI_IRecv( rbuf(:,iproc), llen, ltype, lproc, &
                 itag, lcomm, ireq(nreq), ierr )
         end if
      end do
      rreq = nreq
      do iproc = neighnum,1,-1
         sproc = neighlist(iproc)  ! Send to
         ! Build up list of points
         iqq = 0
         if ( fsvwu ) then
            do k = 1,kx
               iqz = iqq - bnds(sproc)%slen_sv_bg + 1 
               do concurrent (iq = bnds(sproc)%slen_sv_bg:bnds(sproc)%slen_wu_fn)
                  ! send_list_uv(iq) is point index.
                  ! Use abs because sign is used as u/v flag
                  if ( bnds(sproc)%send_list_uv(iq)>0 .neqv. bnds(sproc)%send_swap(iq) ) then
                     sbuf(iqz+iq,iproc) = bnds(sproc)%send_neg(iq)*u(abs(bnds(sproc)%send_list_uv(iq)),k)
                  else
                     sbuf(iqz+iq,iproc) = bnds(sproc)%send_neg(iq)*v(abs(bnds(sproc)%send_list_uv(iq)),k)
                 end if
               end do
               iqq = iqq + bnds(sproc)%slen_wu_fn - bnds(sproc)%slen_sv_bg + 1
            end do   
         end if
         if ( fnveu ) then
            do k = 1,kx 
               iqz = iqq - bnds(sproc)%slen_nv_bg + 1 
               do concurrent (iq = bnds(sproc)%slen_nv_bg:bnds(sproc)%slen_eu_fn)
                  ! send_list_uv(iq) is point index.
                  ! Use abs because sign is used as u/v flag
                  if ( bnds(sproc)%send_list_uv(iq)>0 .neqv. bnds(sproc)%send_swap(iq) ) then
                     sbuf(iqz+iq,iproc) = bnds(sproc)%send_neg(iq)*u(abs(bnds(sproc)%send_list_uv(iq)),k)
                  else
                     sbuf(iqz+iq,iproc) = bnds(sproc)%send_neg(iq)*v(abs(bnds(sproc)%send_list_uv(iq)),k)
                  end if 
               end do
               iqq = iqq + bnds(sproc)%slen_eu_fn - bnds(sproc)%slen_nv_bg + 1
            end do   
         end if
         if ( fssvwwu ) then
            do k = 1,kx 
               iqz = iqq - bnds(sproc)%slen_ssv_bg + 1 
               do concurrent (iq = bnds(sproc)%slen_ssv_bg:bnds(sproc)%slen_wwu_fn)
                  ! send_list_uv(iq) is point index.
                  ! Use abs because sign is used as u/v flag
                  if ( bnds(sproc)%send_list_uv(iq)>0 .neqv. bnds(sproc)%send_swap(iq) ) then
                     sbuf(iqz+iq,iproc) = bnds(sproc)%send_neg(iq)*u(abs(bnds(sproc)%send_list_uv(iq)),k)
                  else
                     sbuf(iqz+iq,iproc) = bnds(sproc)%send_neg(iq)*v(abs(bnds(sproc)%send_list_uv(iq)),k)
                  end if 
               end do
               iqq = iqq + bnds(sproc)%slen_wwu_fn - bnds(sproc)%slen_ssv_bg + 1
            end do   
         end if
         if ( fnnveeu ) then
            do k = 1,kx
               iqz = iqq - bnds(sproc)%slen_nnv_bg + 1 
               do concurrent (iq = bnds(sproc)%slen_nnv_bg:bnds(sproc)%slen_eeu_fn)
                  ! send_list_uv(iq) is point index.
                  ! Use abs because sign is used as u/v flag
                  if ( bnds(sproc)%send_list_uv(iq)>0 .neqv. bnds(sproc)%send_swap(iq) ) then
                     sbuf(iqz+iq,iproc) = bnds(sproc)%send_neg(iq)*u(abs(bnds(sproc)%send_list_uv(iq)),k)
                  else
                     sbuf(iqz+iq,iproc) = bnds(sproc)%send_neg(iq)*v(abs(bnds(sproc)%send_list_uv(iq)),k)
                  end if
               end do
               iqq = iqq + bnds(sproc)%slen_eeu_fn - bnds(sproc)%slen_nnv_bg + 1
            end do   
         end if
         if ( fsuev ) then
            do k = 1,kx 
               iqz = iqq - bnds(sproc)%slen_su_bg + 1 
               do concurrent (iq = bnds(sproc)%slen_su_bg:bnds(sproc)%slen_ev_fn)
                  ! send_list_uv(iq) is point index.
                  ! Use abs because sign is used as u/v flag
                  if ( bnds(sproc)%send_list_uv(iq)>0 .neqv. bnds(sproc)%send_swap(iq) ) then
                     sbuf(iqz+iq,iproc) = bnds(sproc)%send_neg(iq)*u(abs(bnds(sproc)%send_list_uv(iq)),k)
                  else
                     sbuf(iqz+iq,iproc) = bnds(sproc)%send_neg(iq)*v(abs(bnds(sproc)%send_list_uv(iq)),k)
                 end if
              end do
              iqq = iqq + bnds(sproc)%slen_ev_fn - bnds(sproc)%slen_su_bg + 1
            end do  
         end if
         if ( fnnueev ) then
            do k = 1,kx 
               iqz = iqq - bnds(sproc)%slen_nnu_bg + 1 
               do concurrent (iq = bnds(sproc)%slen_nnu_bg:bnds(sproc)%slen_eev_fn)
                  ! send_list_uv(iq) is point index.
                  ! Use abs because sign is used as u/v flag
                  if ( bnds(sproc)%send_list_uv(iq)>0 .neqv. bnds(sproc)%send_swap(iq) ) then
                     sbuf(iqz+iq,iproc) = bnds(sproc)%send_neg(iq)*u(abs(bnds(sproc)%send_list_uv(iq)),k)
                  else
                     sbuf(iqz+iq,iproc) = bnds(sproc)%send_neg(iq)*v(abs(bnds(sproc)%send_list_uv(iq)),k)
                  end if
               end do
               iqq = iqq + bnds(sproc)%slen_eev_fn - bnds(sproc)%slen_nnu_bg + 1
            end do   
         end if
         if ( iqq > 0 ) then
            nreq = nreq + 1
            llen = iqq
            lproc = sproc
            call MPI_ISend( sbuf(:,iproc), llen, ltype, lproc, &
                 itag, lcomm, ireq(nreq), ierr )
         end if
      end do

      ! See if there are any points on my own processor that need
      ! to be fixed up. This will only be in the case when nproc < npanels.
      do k = 1,kx
         do concurrent (iq = 1:myrlen)
            ! request_list is same as send_list in this case
            ! unpack_list(iq) is index into extended region
            if ( bnds(myid)%unpack_list_uv(iq)>0 .and. (bnds(myid)%request_list_uv(iq)>0 .neqv. bnds(myid)%uv_swap(iq)) ) then
               u(ifull+bnds(myid)%unpack_list_uv(iq),k) = bnds(myid)%uv_neg(iq)*u(abs(bnds(myid)%request_list_uv(iq)),k)
            else if ( bnds(myid)%request_list_uv(iq)>0 .neqv. bnds(myid)%uv_swap(iq) ) then
               v(ifull-bnds(myid)%unpack_list_uv(iq),k) = bnds(myid)%uv_neg(iq)*u(abs(bnds(myid)%request_list_uv(iq)),k)
            else if ( bnds(myid)%unpack_list_uv(iq)>0 ) then
               u(ifull+bnds(myid)%unpack_list_uv(iq),k) = bnds(myid)%uv_neg(iq)*v(abs(bnds(myid)%request_list_uv(iq)),k)
            else
               v(ifull-bnds(myid)%unpack_list_uv(iq),k) = bnds(myid)%uv_neg(iq)*v(abs(bnds(myid)%request_list_uv(iq)),k)
            end if
         end do   
      end do
      
      ! Unpack incomming messages
      rcount = rreq
      do while ( rcount > 0 )
         call START_LOG(mpiwait_begin)
         call MPI_Waitsome( rreq, ireq, ldone, donelist, status, ierr )
         call END_LOG(mpiwait_end)
         rcount = rcount - ldone
         do jproc = 1,ldone
            mproc = donelist(jproc)
            iproc = rlist(mproc)  ! Recv from
            rproc = neighlist(iproc)
            iqq = 0
            if ( fsvwu ) then
               do k = 1,kx  
                  iqz = iqq - bnds(rproc)%rlen_sv_bg + 1 
                  do concurrent (iq = bnds(rproc)%rlen_sv_bg:bnds(rproc)%rlen_wu_fn)
                     ! unpack_list(iq) is index into extended region
                     if ( bnds(rproc)%unpack_list_uv(iq) > 0 ) then
                        u(ifull+bnds(rproc)%unpack_list_uv(iq),k) = rbuf(iqz+iq,iproc)
                     else
                        v(ifull-bnds(rproc)%unpack_list_uv(iq),k) = rbuf(iqz+iq,iproc)
                     end if
                  end do
                  iqq = iqq + bnds(rproc)%rlen_wu_fn - bnds(rproc)%rlen_sv_bg + 1
               end do   
            end if
            if ( fnveu ) then
               do k = 1,kx 
                  iqz = iqq - bnds(rproc)%rlen_nv_bg + 1 
                  do concurrent (iq = bnds(rproc)%rlen_nv_bg:bnds(rproc)%rlen_eu_fn)
                     ! unpack_list(iq) is index into extended region
                     if ( bnds(rproc)%unpack_list_uv(iq) > 0 ) then
                        u(ifull+bnds(rproc)%unpack_list_uv(iq),k) = rbuf(iqz+iq,iproc)
                     else
                        v(ifull-bnds(rproc)%unpack_list_uv(iq),k) = rbuf(iqz+iq,iproc)
                     end if
                  end do
                  iqq = iqq + bnds(rproc)%rlen_eu_fn - bnds(rproc)%rlen_nv_bg + 1
               end do   
            end if         
            if ( fssvwwu ) then
               do k = 1,kx 
                  iqz = iqq - bnds(rproc)%rlen_ssv_bg + 1 
                  do concurrent (iq = bnds(rproc)%rlen_ssv_bg:bnds(rproc)%rlen_wwu_fn)
                     ! unpack_list(iq) is index into extended region
                     if ( bnds(rproc)%unpack_list_uv(iq) > 0 ) then
                        u(ifull+bnds(rproc)%unpack_list_uv(iq),k) = rbuf(iqz+iq,iproc)
                     else
                        v(ifull-bnds(rproc)%unpack_list_uv(iq),k) = rbuf(iqz+iq,iproc)
                     end if
                  end do
                  iqq = iqq + bnds(rproc)%rlen_wwu_fn - bnds(rproc)%rlen_ssv_bg + 1
               end do  
            end if         
            if ( fnnveeu ) then
               do k = 1,kx 
                  iqz = iqq - bnds(rproc)%rlen_nnv_bg + 1 
                  do concurrent (iq = bnds(rproc)%rlen_nnv_bg:bnds(rproc)%rlen_eeu_fn)
                     ! unpack_list(iq) is index into extended region
                     if ( bnds(rproc)%unpack_list_uv(iq) > 0 ) then
                        u(ifull+bnds(rproc)%unpack_list_uv(iq),k) = rbuf(iqz+iq,iproc)
                     else
                        v(ifull-bnds(rproc)%unpack_list_uv(iq),k) = rbuf(iqz+iq,iproc)
                     end if
                  end do
                  iqq = iqq + bnds(rproc)%rlen_eeu_fn - bnds(rproc)%rlen_nnv_bg + 1
               end do   
            end if     
            if ( fsuev ) then
               do k = 1,kx 
                  iqz = iqq - bnds(rproc)%rlen_su_bg + 1 
                  do concurrent (iq = bnds(rproc)%rlen_su_bg:bnds(rproc)%rlen_ev_fn)
                     ! unpack_list(iq) is index into extended region
                     if ( bnds(rproc)%unpack_list_uv(iq) > 0 ) then
                        u(ifull+bnds(rproc)%unpack_list_uv(iq),k) = rbuf(iqz+iq,iproc)
                     else
                        v(ifull-bnds(rproc)%unpack_list_uv(iq),k) = rbuf(iqz+iq,iproc)
                     end if
                  end do
                  iqq = iqq + bnds(rproc)%rlen_ev_fn - bnds(rproc)%rlen_su_bg + 1
               end do
            end if
            if ( fnnueev ) then
               do k = 1,kx 
                  iqz = iqq - bnds(rproc)%rlen_nnu_bg + 1 
                  do concurrent (iq = bnds(rproc)%rlen_nnu_bg:bnds(rproc)%rlen_eev_fn)
                     ! unpack_list(iq) is index into extended region
                     if ( bnds(rproc)%unpack_list_uv(iq) > 0 ) then
                        u(ifull+bnds(rproc)%unpack_list_uv(iq),k) = rbuf(iqz+iq,iproc)
                     else
                        v(ifull-bnds(rproc)%unpack_list_uv(iq),k) = rbuf(iqz+iq,iproc)
                     end if
                  end do
                  iqq = iqq + bnds(rproc)%rlen_eev_fn - bnds(rproc)%rlen_nnu_bg + 1
               end do   
            end if
         end do
      end do

      ! Clear any remaining messages
      sreq = nreq - rreq
      if ( sreq > 0 ) then
         call START_LOG(mpiwait_begin)
         call MPI_Waitall( sreq, ireq(rreq+1:nreq), status, ierr )
         call END_LOG(mpiwait_end)
      end if

   end subroutine boundsuv3
   
   subroutine deptsync(nface,xg,yg)
      ! Different levels will have different winds, so the list of points is
      ! different on each level.
      ! xg ranges from 0.5 to il+0.5 on a face. A given processors range
      ! is 0.5+ioff to 0.5+ipan+ioff
      ! Assignment of points to processors needs to match what ints does
      ! in case there's an exact edge point.
      ! Because of the boundary region, the range [0:ipan+1) can be handled.
      ! Need floor(xxg) in range [0:ipan]
      integer, dimension(:,:), intent(in) :: nface
      real, dimension(:,:), intent(in) :: xg, yg
      integer :: iproc, jproc, dproc
      integer :: ip, jp, xn, kx
      integer :: iq, k, idel, jdel, nf
      integer :: rcount
      integer(kind=4) :: itag=99, ierr, llen, ncount, sreq, lproc
      integer(kind=4) :: ldone, lcomm
      integer(kind=4), dimension(MPI_STATUS_SIZE,size(ireq)) :: status
      integer(kind=4), dimension(neighnum) :: donelist
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif   

      ! This does nothing in the one processor case
      if ( neighnum < 1 ) return
      
      kx = size(nface,2)
      dslen(:) = 0
      drlen(:) = 0
      dproc = 0
      lcomm = comm_world
      
      ! In this case the length of each buffer is unknown and will not
      ! be symmetric between processors. Therefore need to get the length
      ! from the message status
      do iproc = 1,neighnum
         lproc = neighlist(iproc)  ! Recv from
         ! Use the maximum size in the recv call.
         llen = 4*bnds(lproc)%len
         call MPI_IRecv( dpoints(iproc)%a, llen, ltype, lproc, &
                      itag, lcomm, ireq(iproc), ierr )
      end do
      nreq = neighnum
      
      ! Calculate request list
      do k = 1,kx
         do iq = 1,ifull
            nf = nface(iq,k) + noff ! Make this a local index
            idel = int(xg(iq,k)) - ioff
            jdel = int(yg(iq,k)) - joff
            if ( idel<0 .or. idel>ipan .or. jdel<0 .or. jdel>jpan .or. nf<1 .or. nf>npan ) then
               ! If point is on a different processor, add to a list 
               ip = min( il_g, max( 1, nint(xg(iq,k)) ) )
               jp = min( il_g, max( 1, nint(yg(iq,k)) ) )
               iproc = fproc(ip,jp,nface(iq,k)) ! processor that owns global grid point
               dproc = neighmap(iproc) ! returns 0 if not in neighlist
               ! Add this point to the list of requests I need to send to iproc
               dslen(dproc) = dslen(dproc) + 1
               ! Limit request index to valid range to avoid seg fault
               xn = max( min( dslen(dproc), bnds(iproc)%len ), 1 )
               ! Since nface is a small integer it can be exactly represented by a
               ! real. It is simpler to send like this than use a proper structure.
               dbuf(dproc)%a(:,xn) = (/ real(nface(iq,k)), xg(iq,k), yg(iq,k), real(k) /)
               dindex(dproc)%a(xn,1) = iq
               dindex(dproc)%a(xn,2) = k
            end if
         end do
      end do
 
      ! Check for errors
      if ( dslen(0) > 0 ) then
         ip = min( il_g, max( 1, nint(dbuf(0)%a(2,1)) ) )
         jp = min( il_g, max( 1, nint(dbuf(0)%a(3,1)) ) )
         iproc = fproc(ip,jp,nint(dbuf(0)%a(1,1)))
         write(6,*) "myid,dslen,len ",myid,dslen(0),0
         write(6,*) "Example error iq,k,iproc ",dindex(0)%a(1,1),dindex(0)%a(1,2),iproc
         write(6,*) "dbuf ",dbuf(0)%a(:,1)
         write(6,*) "neighlist ",neighlist
         write(6,*) "ERROR: Wind speed is very large and the departure point is"
         write(6,*) "       further away than the neighbouring processes. This"
         write(6,*) "       error could be due to an excessively large"
         write(6,*) "       time-step or alternatively caused by a NaN originating"
         write(6,*) "       from earlier in the simulation."
         call checksize( dslen(0), 0, "Deptsync" )
      end if
      do dproc = 1,neighnum
         iproc = neighlist(dproc)
         if ( dslen(dproc) > bnds(iproc)%len ) then
            write(6,*) "myid,dslen,len ",myid,dslen(dproc),bnds(iproc)%len
            write(6,*) "Example error iq,k,iproc ",dindex(dproc)%a(1,1),dindex(dproc)%a(1,2),iproc
            write(6,*) "dbuf ",dbuf(dproc)%a(:,1)
            write(6,*) "neighlist ",neighlist
            write(6,*) "ERROR: Wind speed is very large and the departure point is"
            write(6,*) "       further away than the neighbouring processes. This"
            write(6,*) "       error could be due to an excessively large"
            write(6,*) "       time-step or alternatively caused by a NaN originating"
            write(6,*) "       from earlier in the simulation."
            call checksize( dslen(dproc), bnds(iproc)%len, "Deptsync" )
         end if
      end do

      ! Send request list
      rreq = nreq
      do iproc = neighnum,1,-1
         lproc = neighlist(iproc)  ! Send to
         ! Send, even if length is zero
         nreq = nreq + 1
         llen = 4*dslen(iproc)
         call MPI_ISend( dbuf(iproc)%a, llen, ltype, lproc, &
                 itag, lcomm, ireq(nreq), ierr )
      end do
      
      ! Unpack incomming messages
      rcount = rreq
      do while ( rcount > 0 )
         call START_LOG(mpiwaitdep_begin)
         call MPI_Waitsome( rreq, ireq, ldone, donelist, status, ierr )
         call END_LOG(mpiwaitdep_end)
         rcount = rcount - ldone
         do jproc = 1,ldone
            ! Now get the actual sizes from the status
            call MPI_Get_count( status(:,jproc), ltype, ncount, ierr )
            drlen(donelist(jproc)) = ncount/4
         end do
      end do

      ! Clear any remaining message requests
      sreq = nreq - rreq
      if ( sreq > 0 ) then
         call START_LOG(mpiwaitdep_begin)
         call MPI_Waitall( sreq, ireq(rreq+1:nreq), status, ierr )
         call END_LOG(mpiwaitdep_end)
      end if

   end subroutine deptsync

   subroutine intssync_send3
      
      call intssync_send4(1)

   end subroutine intssync_send3
   
   subroutine intssync_send4(ntr)
      integer, intent(in) :: ntr
      integer :: iproc
      integer(kind=4) :: itag=98, ierr, llen, lproc, lcomm
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      
      if ( ntr > nagg ) then
         write(6,*) "ERROR: Internal error in intssync_send.  ntr > nagg"
         call ccmpi_abort(-1)
      end if

      ! When sending the results, roles of dslen and drlen are reversed
      lcomm = comm_world
      nreq = 0
      do iproc = 1,neighnum
         if ( dslen(iproc) > 0 ) then
            nreq = nreq + 1
            rlist(nreq) = iproc
            llen = dslen(iproc)*ntr
            lproc = neighlist(iproc)  ! Recv from
            call MPI_IRecv( dbuf(iproc)%b, llen, ltype, lproc, itag, &
                            lcomm, ireq(nreq), ierr )
         end if
      end do
      rreq = nreq
      do iproc = neighnum,1,-1
         if ( drlen(iproc) > 0 ) then
            nreq = nreq + 1
            llen = drlen(iproc)*ntr
            lproc = neighlist(iproc)  ! Send to
            call MPI_ISend( sextra(iproc)%a, llen, ltype, lproc, itag, & 
                            lcomm, ireq(nreq), ierr )
         end if
      end do

   end subroutine intssync_send4

   subroutine intssync_recv3(s)
      real, dimension(:,:), intent(inout) :: s
      integer :: iproc, iq, jproc
      integer :: rcount
      integer(kind=4) :: ierr, ldone, sreq
      integer(kind=4), dimension(neighnum) :: donelist
      integer(kind=4), dimension(MPI_STATUS_SIZE,size(ireq)) :: status
      
      ! Unpack incomming messages
      rcount = rreq
      do while ( rcount > 0 )
         call START_LOG(mpiwaitdep_begin)
         call MPI_Waitsome( rreq, ireq, ldone, donelist, status, ierr )
         call END_LOG(mpiwaitdep_end)
         rcount = rcount - ldone
         do jproc = 1,ldone
            iproc = rlist(donelist(jproc))
            do concurrent (iq = 1:dslen(iproc))
               s(dindex(iproc)%a(iq,1),dindex(iproc)%a(iq,2)) = dbuf(iproc)%b(iq)
            end do   
         end do
      end do

      ! Clear any remaining messages
      sreq = nreq - rreq
      if ( sreq > 0 ) then
         call START_LOG(mpiwaitdep_begin)
         call MPI_Waitall( sreq, ireq(rreq+1:nreq), status, ierr )
         call END_LOG(mpiwaitdep_end)
      end if

   end subroutine intssync_recv3

   subroutine intssync_recv4(s)
      real, dimension(:,:,:), intent(inout) :: s
      integer :: iproc, iq, jproc
      integer :: rcount, ntr, l
      integer(kind=4) :: ierr, ldone, sreq
      integer(kind=4), dimension(neighnum) :: donelist
      integer(kind=4), dimension(MPI_STATUS_SIZE,size(ireq)) :: status

      ntr = size(s,3)
      
      ! Unpack incomming messages
      rcount = rreq
      do while ( rcount > 0 )
         call START_LOG(mpiwaitdep_begin)
         call MPI_Waitsome( rreq, ireq, ldone, donelist, status, ierr )
         call END_LOG(mpiwaitdep_end)
         rcount = rcount - ldone
         do jproc = 1,ldone
            iproc = rlist(donelist(jproc))
            do l = 1,ntr
               do concurrent (iq = 1:dslen(iproc))
                  s(dindex(iproc)%a(iq,1),dindex(iproc)%a(iq,2),l) = dbuf(iproc)%b(iq+(l-1)*dslen(iproc))
               end do
            end do   
         end do
      end do

      ! Clear any remaining messages
      sreq = nreq - rreq
      if ( sreq > 0 ) then
         call START_LOG(mpiwaitdep_begin)
         call MPI_Waitall( sreq, ireq(rreq+1:nreq), status, ierr )
         call END_LOG(mpiwaitdep_end)
      end if

   end subroutine intssync_recv4
   
   subroutine indv_mpi(iq, i, j, n)
      integer , intent(in) :: iq
      integer , intent(out) :: i
      integer , intent(out) :: j
      integer , intent(out) :: n

      ! Calculate local i, j, n from global iq

      ! Global i, j, n
      n = (iq - 1)/(il_g*il_g)
      j = 1 + (iq - n*il_g*il_g - 1)/il_g
      i = iq - (j - 1)*il_g - n*il_g*il_g
      if ( fproc(i,j,n) /= myid ) then
         write(*,"(a,5i5)") "Consistency failure in indv_mpi", myid, iq, i, j, n
         call ccmpi_abort(-1)
      end if
      ! Reduced to values on my processor
      j = j - joff
      i = i - ioff
      n = n + noff      
   end subroutine indv_mpi

   pure function indglobal(i,j,n) result(iq)
      integer, intent(in) :: i, j, n
      integer :: iq

      ! Calculate a 1D global index from the global indices
      ! n in range 0:npanels
      iq = i + (j-1)*il_g + n*il_g*il_g
   end function indglobal

   pure function indg(i,j,n) result(iq)
      integer, intent(in) :: i, j, n
      integer :: iq

      ! Calculate a 1D global index from the local processors indices
      ! n in range 1..npan
      iq = i+ioff + (j+joff-1)*il_g + (n-noff)*il_g*il_g
   end function indg

   pure function indp(i,j,n) result(iq)
      integer, intent(in) :: i, j, n
      integer :: iq

      ! Calculate a 1D local index from the local processors indices
      ! Note that face number runs from 1 here.
      iq = i + (j-1)*ipan + (n-1)*ipan*jpan
   end function indp

  pure function iq2iqg(iq) result(iqg)
      integer, intent(in) :: iq
      integer :: iqg
      integer :: i, j, n

      ! Calculate global iqg from local iq

      ! MJT bug fix (should be ipan and jpan, not il and jl)
      ! Calculate local i, j, n
      n = 1 + (iq-1)/(ipan*jpan)  ! In range 1 .. npan
      j = 1 + ( iq - (n-1)*(ipan*jpan) - 1) / ipan
      i = iq - (j-1)*ipan - (n-1)*(ipan*jpan)
      iqg = indg(i,j,n)

   end function iq2iqg

   pure function fproc(i,j,n) result(fpout)
      ! locates processor that owns a global grid point
      integer, intent(in) :: i, j, n
      integer :: fpout
      integer :: ip, jp

      ip = (i-1)/ipan
      jp = (j-1)/jpan
      if ( uniform_decomp ) then
         fpout = ip + jp*nxproc
      else
         fpout = ip + jp*nxproc + n*nxproc*nyproc/npan
      end if
   
   end function fproc

   pure function qproc(iqg) result(qpout)
      ! locates processor that owns a global grid point
      integer, intent(in) :: iqg
      integer :: qpout
      integer :: i, j, n

      n = (iqg - 1) / (il_g*il_g)
      j = 1 + (iqg - n*il_g*il_g - 1)/il_g
      i = iqg - (j - 1)*il_g - n*il_g*il_g

      qpout = fproc(i,j,n)
   
   end function qproc

   pure function indx_indv(iq_g, mil_g, mipan, mjpan, mioff, mjoff, mnoff) result(iq)
      ! converts a global index to a local index
      integer, intent(in) :: iq_g, mil_g, mipan, mjpan, mioff, mjoff, mnoff
      integer iq
      integer i, j, n
      
      n = (iq_g - 1)/(mil_g*mil_g)
      j = 1 + (iq_g - n*mil_g*mil_g - 1)/mil_g
      i = iq_g - (j - 1)*mil_g - n*mil_g*mil_g

      n = n + mnoff  
      j = j - mjoff
      i = i - mioff
      
      iq = i + (j-1)*mipan + (n-1)*mipan*mjpan
   
   end function indx_indv

   pure function indx(i,j,n,il,jl) result(iq)
      ! more general version of ind function

      integer, intent(in) :: i, j, n, il, jl
      integer iq

      iq = i + (j-1)*il + n*il*jl

   end function indx
   
   pure function ind(i,j,n) result(iq)

      integer, intent(in) :: i, j, n
      integer iq

      iq = i + (j-1)*ipan + (n-1)*ipan*jpan

   end function ind
   
   pure function findcolour(iqg) result(icol)
   
      integer, intent(in) :: iqg
      integer icol
      integer ig, jg, ng, tg

      icol = -1 ! for gfortran
      
      ! calculate global i,j,n
      tg = iqg - 1
      ng = tg/(il_g*il_g)
      tg = tg - ng*il_g*il_g
      jg = tg/il_g
      tg = tg - jg*il_g
      ig = tg
      ig = ig + 1
      jg = jg + 1
   
      ! MJT notes - we use two colours for both
      ! uniform_decomp and face_decomp to ensure
      ! the results are bit-reproducible
      icol = mod( ig + jg + ng*il_g, 2 ) + 1
   
   end function findcolour
   
   subroutine checksize(len, msize, mesg)
   
      integer, intent(in) :: len
      integer, intent(in) :: msize
      character(len=*), intent(in) :: mesg
      
      if ( len > msize ) then
         write(6,*) "Error, maxsize exceeded in ", mesg
         call ccmpi_abort(-1)
      end if
      
   end subroutine checksize

   subroutine fix_index2(iqq,larray,n,bnds,iext)
      integer, intent(in) :: iqq, n
      integer, dimension(:), intent(out) :: larray
      integer, intent(inout) :: iext
      type(bounds_info), dimension(0:), intent(inout) :: bnds
      integer :: rproc
      integer :: iloc,jloc,nloc

      ! Which processor has this point
      rproc = qproc(iqq)
      if ( rproc /= myid ) then ! Add to list
         bnds(rproc)%rlen2 = bnds(rproc)%rlen2 + 1 
         call check_bnds_alloc(rproc, iext)
         bnds(rproc)%request_list(bnds(rproc)%rlen2) = iqq
         ! Increment extended region index
         iext = iext + 1
         bnds(rproc)%unpack_list(bnds(rproc)%rlen2) = iext
         larray(n) = ifull + iext
      else
         ! If it's on this processor, just use the local index
         call indv_mpi(iqq,iloc,jloc,nloc)
         larray(n) = indp(iloc,jloc,nloc)
      end if
      
   end subroutine fix_index2

   subroutine proc_setup(id,jd,idjd)
!     Routine to set up offsets etc.
      integer :: i, j, n, nd, jdf, idjd_g
      integer, intent(in) :: id, jd
      integer, intent(out) :: idjd
      integer, dimension(0:npanels) :: ipoff, jpoff

      call face_set( ipan, jpan, noff, ipoff, jpoff, npan, il_g, myid, nproc, nxproc, nyproc )
      ioff = ipoff(0)
      joff = jpoff(0)

      ! Convert standard jd to a face index
      nd = (jd-1)/il_g ! 0: to match fproc
      jdf = jd - nd*il_g
      mydiag = ( myid == fproc(id,jdf,nd) )
      ! Convert global indices to ones on this processors region
      idjd_g = id + (jd-1)*il_g
      if ( mydiag ) then
         call indv_mpi(idjd_g,i,j,n)
         idjd = indp(i,j,n)
      else
         ! This should never be used so set a value that will give a bounds error
         idjd = huge(1)
      end if

   end subroutine proc_setup

   subroutine proc_setup_uniform(id,jd,idjd)
!     Routine to set up offsets etc for the uniform decomposition
      integer :: i, j, n, nd, jdf, idjd_g
      integer, intent(in) :: id, jd
      integer, intent(out) :: idjd
      integer, dimension(0:npanels) :: ipoff, jpoff

      call dix_set( ipan, jpan, noff, ipoff, jpoff, npan, il_g, myid, nproc, nxproc, nyproc)
      ioff = ipoff(0)
      joff = jpoff(0)

!     Check that the values calculated here match those set as parameters
      if ( ipan /= il ) then
         write(6,*) "Error, parameter mismatch, ipan /= il", ipan, il
         call ccmpi_abort(-1)
      end if
      if ( jpan*npan /= jl ) then
         write(6,*) "Error, parameter mismatch, jpan*npan /= jl", jpan, npan, jl
         call ccmpi_abort(-1)
      end if

      ! Convert standard jd to a face index
      nd = (jd-1)/il_g ! 0: to match fproc
      jdf = jd - nd*il_g
      mydiag = ( myid == fproc(id,jdf,nd) )
      ! Convert global indices to ones on this processors region
      idjd_g = id + (jd-1)*il_g
      if ( mydiag ) then
         call indv_mpi(idjd_g,i,j,n)
         idjd = indp(i,j,n)
      else
         ! This should never be used so set a value that will give a bounds error
         idjd = huge(1)
      end if

   end subroutine proc_setup_uniform

   subroutine face_set(ipan_l, jpan_l, noff_l, ioff_l, joff_l, npan_l, il_gx, myid_l, nproc_l, nxproc_l, nyproc_l)
      integer, intent(in) :: myid_l, nproc_l, npan_l, il_gx
      integer, intent(out) :: ipan_l, jpan_l, noff_l, nxproc_l, nyproc_l
      integer, dimension(0:npanels), intent(out) :: ioff_l, joff_l 
      integer :: n

      !  Processor allocation
      !  if  nproc_l <= npanels+1, then each gets a number of full panels
      if ( nproc_l<=npanels+1 ) then
         if ( modulo( npanels+1, nproc_l )/=0 ) then
            write(6,*) "Error, number of processors must divide number of panels"
            call ccmpi_abort(-1)
         end if
!        npan_l = (npanels+1)/nproc_l
         ipan_l = il_gx
         jpan_l = il_gx
         noff_l = 1 - myid_l*npan_l
         ioff_l(:) = 0
         joff_l(:) = 0
         nxproc_l = 1
         nyproc_l = 1
      else  ! nproc_l >= npanels+1
         if ( modulo( nproc_l, npanels+1 )/=0 ) then
            write(6,*) "Error, number of processors must be a multiple of number of panels"
            call ccmpi_abort(-1)
         end if
!        npan_l = 1
         n = nproc_l/(npanels+1)
         !  n is the number of processors on each face
         !  Try to factor this into two values are close as possible.
         !  nxproc is the smaller of the 2.
         nxproc_l = nint(sqrt(real(n)))
         nyproc_l = n / nxproc_l
         do nxproc_l = nint(sqrt(real(n))), 1, -1
            nyproc_l = n / nxproc_l
            if ( modulo( il_gx, nxproc_l )==0 .and. modulo( il_gx, nyproc_l )==0 .and. &
                 nxproc_l*nyproc_l==n ) exit
         end do
         if ( nxproc_l*nyproc_l/=n ) then
            write(6,*) "Error in splitting up faces"
            call ccmpi_abort(-1)
         end if

         ! Still need to check that the processor distribution is compatible
         ! with the grid.
         if ( modulo( il_gx, nxproc_l )/=0 ) then
            write(6,*) "Error, il not a multiple of nxproc", il_gx, nxproc_l
            call ccmpi_abort(-1)
         end if
         if ( modulo( il_gx, nyproc_l )/=0 ) then
            write(6,*) "Error, il not a multiple of nyproc", il_gx, nyproc_l
            call ccmpi_abort(-1)
         end if
         ipan_l = il_gx/nxproc_l
         jpan_l = il_gx/nyproc_l

         ! Set offsets for this processor
         call proc_region_face(myid_l,ioff_l(0),joff_l(0),noff_l,nxproc_l,nyproc_l,ipan_l,jpan_l,npan_l)
         ioff_l(1:npanels) = ioff_l(0)
         joff_l(1:npanels) = joff_l(0)
      end if
   
   end subroutine face_set

   subroutine dix_set(ipan_l,jpan_l,noff_l,ioff_l,joff_l,npan_l,il_gx,myid_l,nproc_l,nxproc_l,nyproc_l)
      integer, intent(in) :: myid_l, nproc_l, npan_l, il_gx
      integer, intent(out) :: ipan_l, jpan_l, noff_l, nxproc_l, nyproc_l
      integer, dimension(0:npanels), intent(out) :: ioff_l, joff_l 
      
      if ( npan_l /= npanels+1 ) then
         write(6,*) "Error: inconsistency in proc_setup_uniform"
         call ccmpi_abort(-1)
      end if
      !  Processor allocation: each processor gets a part of each panel
      !  Try to factor nproc into two values are close as possible.
      !  nxproc is the smaller of the 2.
      nxproc_l = nint(sqrt(real(nproc_l)))
      do nxproc_l = nint(sqrt(real(nproc_l))),1,-1
         ! This will always exit eventually because it's trivially true 
         ! for nxproc=1
         nyproc_l = nproc_l/nxproc_l
         if ( modulo( nproc_l, nxproc_l )==0 .and. &
              modulo( il_gx, nxproc_l )==0  .and.  &
              modulo( il_gx, nyproc_l )==0 ) exit
      end do
      nyproc_l = nproc_l/nxproc_l
      if ( nxproc_l*nyproc_l/=nproc_l ) then
         write(6,*) "Error in splitting up faces"
         call ccmpi_abort(-1)
      end if

      ! Still need to check that the processor distribution is compatible
      ! with the grid.
      if ( modulo( il_gx, nxproc_l )/=0 ) then
         write(6,*) "Error, il not a multiple of nxproc", il_gx, nxproc_l
         call ccmpi_abort(-1)
      end if
      if ( modulo( il_gx, nyproc_l )/=0 ) then
         write(6,*) "Error, il not a multiple of nyproc", il_gx, nyproc_l
         call ccmpi_abort(-1)
      end if
      ipan_l = il_gx/nxproc_l
      jpan_l = il_gx/nyproc_l

      ! Set offsets for this processor
      call proc_region_dix(myid_l,ioff_l(0),joff_l(0),noff_l,nxproc_l,ipan_l,jpan_l)
      ioff_l(1:npanels)=ioff_l(0)
      joff_l(1:npanels)=joff_l(0)

   end subroutine dix_set

   subroutine uniform_set(ipan_l,jpan_l,noff_l,ioff_l,joff_l,npan_l,il_gx,myid_l,nproc_l,nxproc_l,nyproc_l)
      ! backwards compatibility
      integer, intent(in) :: myid_l, nproc_l, npan_l, il_gx
      integer, intent(out) :: ipan_l, jpan_l, noff_l, nxproc_l, nyproc_l
      integer, dimension(0:npanels), intent(out) :: ioff_l, joff_l 
      integer n
      
      if ( npan_l /= npanels+1 ) then
         write(6,*) "Error: inconsistency in proc_setup_uniform"
         call ccmpi_abort(-1)
      end if
      !  Processor allocation: each processor gets a part of each panel
      !  Try to factor nproc into two values are close as possible.
      !  nxproc is the smaller of the 2.
      nxproc_l = nint(sqrt(real(nproc_l)))
      do nxproc_l = nint(sqrt(real(nproc_l))), 1, -1
         ! This will always exit eventually because it's trivially true 
         ! for nxproc=1
         nyproc_l = nproc_l/nxproc_l
         if ( modulo( nproc_l, nxproc_l )==0 .and. &
              modulo( il_gx, nxproc_l )==0  .and.  &
              modulo( il_gx, nyproc_l )==0 ) exit
      end do
      nyproc_l = nproc_l/nxproc_l
      if ( nxproc_l*nyproc_l/=nproc_l ) then
         write(6,*) "Error in splitting up faces"
         call ccmpi_abort(-1)
      end if

      ! Still need to check that the processor distribution is compatible
      ! with the grid.
      if ( modulo( il_gx, nxproc_l )/=0 ) then
         write(6,*) "Error, il not a multiple of nxproc", il_gx, nxproc_l
         call ccmpi_abort(-1)
      end if
      if ( modulo( il_gx, nyproc_l )/=0 ) then
         write(6,*) "Error, il not a multiple of nyproc", il_gx, nyproc_l
         call ccmpi_abort(-1)
      end if
      ipan_l = il_gx/nxproc_l
      jpan_l = il_gx/nyproc_l

      ! Set offsets for this processor
      do n = 0,npanels
         call proc_region_uniform(myid_l,n,ioff_l(n),joff_l(n),noff_l,nxproc_l,nyproc_l,ipan_l,jpan_l)
      end do

   end subroutine uniform_set

   pure subroutine proc_region_face(procid,ipoff,jpoff,npoff,nxproc_l,nyproc_l,ipan_l,jpan_l,npan_l)
      ! Calculate the offsets for a given processor
      integer, intent(in) :: procid, nxproc_l, nyproc_l, ipan_l, jpan_l, npan_l
      integer, intent(out) :: ipoff, jpoff, npoff
      integer :: myface, mtmp, nproc_l

      nproc_l = nxproc_l*nyproc_l*(npanels+1)/npan_l
      if ( nproc_l <= npanels+1 ) then
         npoff = 1 - procid*npan_l
         ipoff = 0
         jpoff = 0
      else
         myface = procid/(nxproc_l*nyproc_l)
         npoff = 1 - myface
         ! mtmp is the processor index on this face, 0:(nxprox*nyproc-1)
         mtmp = procid - myface*nxproc_l*nyproc_l
         jpoff = (mtmp/nxproc_l)*jpan_l
         ipoff = modulo( mtmp, nxproc_l )*ipan_l
      end if
     
   end subroutine proc_region_face

   pure subroutine proc_region_dix(procid,ipoff,jpoff,npoff,nxproc_l,ipan_l,jpan_l)
      ! Calculate the offsets for a given processor
      integer, intent(in) :: procid, nxproc_l, ipan_l, jpan_l
      integer, intent(out) :: ipoff, jpoff, npoff

      ! Original Dix uniform decomposition
      ! Set offsets for this processor (same on all faces)
      npoff = 1
      jpoff = (procid/nxproc_l) * jpan_l
      ipoff = modulo(procid,nxproc_l)*ipan_l
     
   end subroutine proc_region_dix
   
   pure subroutine proc_region_uniform(procid,panid,ipoff,jpoff,npoff,nxproc_l,nyproc_l,ipan_l,jpan_l)
      ! Calculate the offsets for a given processor
      integer, intent(in) :: procid, panid, nxproc_l, nyproc_l, ipan_l, jpan_l
      integer, intent(out) :: ipoff, jpoff, npoff
      integer il_gx

      ! MJT suggested decomposition to improve load balance
      il_gx = ipan_l*nxproc_l
      npoff = 1
      select case(panid)
         case(0)
            jpoff = (procid/nxproc_l) * jpan_l
            ipoff = modulo( procid, nxproc_l ) * ipan_l
            if ( jpoff >= il_gx/2 ) then
               ipoff = il_gx - ipoff - ipan_l
            end if
         case(1)
            jpoff = (procid/nxproc_l) * jpan_l
            ipoff = modulo( procid, nxproc_l ) * ipan_l
         case(2)
            jpoff = (procid/nxproc_l) * jpan_l
            ipoff = modulo( procid, nxproc_l ) * ipan_l
            if ( ipoff >= il_gx/2 ) then
               jpoff = il_gx - jpoff - jpan_l
            end if
         case(3)
            jpoff = modulo( procid, nyproc_l ) * jpan_l
            ipoff = (procid/nyproc_l) * ipan_l
            if ( ipoff >= il_gx/2 ) then
               jpoff = il_gx - jpoff - jpan_l
            end if
         case(4)
            jpoff = modulo( procid, nyproc_l ) * jpan_l
            ipoff = (procid/nyproc_l) * ipan_l
         case(5)
            jpoff = modulo( procid, nyproc_l ) * jpan_l
            ipoff = (procid/nyproc_l) * ipan_l
            if ( jpoff >= il_gx/2 ) then
               ipoff = il_gx - ipoff - ipan_l
            end if        
      end select
     
   end subroutine proc_region_uniform

   subroutine start_log ( event )
      integer, intent(in) :: event
      if ( ccomp_get_thread_num() /= 0 ) return
#ifdef vampir
      VT_USER_START(event_name(event))
#endif
#ifdef simple_timer
      !start_time(event) = MPI_Wtime()
      call system_clock( start_time(event) ) 
#endif 
   end subroutine start_log

   subroutine end_log ( event )
      integer, intent(in) :: event
      integer(kind=8) :: end_time
      if ( ccomp_get_thread_num() /= 0 ) return
#ifdef vampir
      VT_USER_END(event_name(event))
#endif
#ifdef simple_timer
      !tot_time(event) = tot_time(event) + MPI_Wtime() - start_time(event)
      call system_clock( end_time )
      tot_time(event) = tot_time(event) + (end_time-start_time(event))
#endif 
   end subroutine end_log

   subroutine log_off()
#ifdef vampir
       VT_OFF()
#endif
   end subroutine log_off
   
   subroutine log_on()
#ifdef vampir
      VT_ON()
#endif
   end subroutine log_on
   
   subroutine log_flush()
#ifdef vampir   
    VT_BUFFER_FLUSH()
#endif
   end subroutine log_flush

   subroutine log_setup()
#ifdef vampir
#ifdef simple_timer
      write(6,*) "ERROR: vampir and simple_timer should not be compiled together"
      call ccmpi_abort(-1)
#endif
#endif

      call add_event(maincalc_begin,      maincalc_end,      "MainCalc")
      call add_event(indata_begin,        indata_end,        "Indata")
      call add_event(phys_begin,          phys_end,          "Phys")
      call add_event(nonlin_begin,        nonlin_end,        "Nonlin")
      call add_event(upglobal_begin,      upglobal_end,      "Upglobal")
      call add_event(adjust_begin,        adjust_end,        "Adjust")
      call add_event(hordifg_begin,       hordifg_end,       "Hordifg")
      call add_event(helm_begin,          helm_end,          "Adv_Helm")
      call add_event(vadv_begin,          vadv_end,          "Adv_Vadv")
      call add_event(depts_begin,         depts_end,         "Adv_Depts")
      call add_event(ints_begin,          ints_end,          "Adv_Ints")
      call add_event(mfix_begin,          mfix_end,          "Adv_Mfix")
      call add_event(stag_begin,          stag_end,          "Adv_Stag")
      call add_event(waterdynamics_begin, waterdynamics_end, "Waterdynamics")
      call add_event(waterdiff_begin,     waterdiff_end,     "Whordifg")
      call add_event(waterunpack_begin,   waterunpack_end,   "Water_Unpack")
      call add_event(watereos_begin,      watereos_end,      "Water_EOS")
      call add_event(waterhelm_begin,     waterhelm_end,     "Water_Helm")
      call add_event(watervadv_begin,     watervadv_end,     "Water_Vadv")
      call add_event(waterdeps_begin,     waterdeps_end,     "Water_Deps")
      call add_event(waterints_begin,     waterints_end,     "Water_Ints")
      call add_event(wateriadv_begin,     wateriadv_end,     "Water_Iadv")
      call add_event(watermfix_begin,     watermfix_end,     "Water_Mfix")
      call add_event(ocnstag_begin,       ocnstag_end,       "Water_Stag")
      call add_event(waterpack_begin,     waterpack_end,     "Water_Pack")
      call add_event(river_begin,         river_end,         "River")
      call add_event(outfile_begin,       outfile_end,       "Outfile")
      call add_event(onthefly_begin,      onthefly_end,      "Onthefly")
      call add_event(otf_fill_begin,      otf_fill_end,      "OTF_Fill")
      call add_event(otf_ints1_begin,     otf_ints1_end,     "OTF_Doints1")
      call add_event(otf_ints4_begin,     otf_ints4_end,     "OTF_Doints4")
      call add_event(histrd3_begin,       histrd3_end,       "OTF_HistRd3")
      call add_event(histrd4_begin,       histrd4_end,       "OTF_HistRd4")
      call add_event(histrd5_begin,       histrd5_end,       "OTF_HistRd5")
      call add_event(nestin_begin,        nestin_end,        "Nestin")
      call add_event(nestotf_begin,       nestotf_end,       "Nest_OTF")
      call add_event(nestwin_begin,       nestwin_end,       "Nest_MAP")
      call add_event(nestcalc_begin,      nestcalc_end,      "Nest_calc")
      call add_event(nestcomm_begin,      nestcomm_end,      "Nest_comm")
      call add_event(ensemble_begin,      ensemble_end,      "Ensemble")
      call add_event(amipsst_begin,       amipsst_end,       "AMIPSST")
      call add_event(gwdrag_begin,        gwdrag_end,        "GWdrag")
      call add_event(convection_begin,    convection_end,    "Convection")
      call add_event(cloud_begin,         cloud_end,         "Cloud")
      call add_event(radnet_begin,        radnet_end,        "Rad")
      call add_event(radinit_begin,       radinit_end,       "Rad_init")
      call add_event(sfluxnet_begin,      sfluxnet_end,      "Sflux")
      call add_event(sfluxwater_begin,    sfluxwater_end,    "Sflux_water")
      call add_event(sfluxland_begin,     sfluxland_end,     "Sflux_land")
      call add_event(sfluxurban_begin,    sfluxurban_end,    "Sflux_urban")
      call add_event(aerosol_begin,       aerosol_end,       "Aerosol")
      call add_event(vertmix_begin,       vertmix_end,       "Vertmix")
      call add_event(precon_begin,        precon_end,        "Precon")
      call add_event(mgsetup_begin,       mgsetup_end,       "MG_Setup")
      call add_event(mgfine_begin,        mgfine_end,        "MG_Fine")
      call add_event(mgup_begin,          mgup_end,          "MG_Up")
      call add_event(mgcoarse_begin,      mgcoarse_end,      "MG_Coarse")
      call add_event(mgdown_begin,        mgdown_end,        "MG_Down")
      call add_event(bcast_begin,         bcast_end,         "MPI_Bcast")
      call add_event(bcastmg_begin,       bcastmg_end,       "MPI_BcastMG")
      call add_event(alltoall_begin,      alltoall_end,      "MPI_AlltoAll")
      call add_event(allgather_begin,     allgather_end,     "MPI_AllGather")
      call add_event(allgathercc_begin,   allgathercc_end,   "MPI_AllGatherCC")
      call add_event(gather_begin,        gather_end,        "MPI_Gather")
      call add_event(gathercc_begin,      gathercc_end,      "MPI_GatherCC")
      call add_event(gathermg_begin,      gathermg_end,      "MPI_GatherMG")
      call add_event(gatherfile_begin,    gatherfile_end,    "MPI_GatherFILE")
      call add_event(scatter_begin,       scatter_end,       "MPI_Scatter")
      call add_event(scattercc_begin,     scattercc_end,     "MPI_ScatterCC")
      call add_event(scatterfile_begin,   scatterfile_end,   "MPI_ScatterFILE")
      call add_event(allreduce_begin,     allreduce_end,     "MPI_AllReduce")
      call add_event(allreducepn_begin,   allreducepn_end,   "MPI_AllReducePN")
      call add_event(reduce_begin,        reduce_end,        "MPI_Reduce")
      call add_event(mpiwait_begin,       mpiwait_end,       "MPI_Wait")
      call add_event(mpiwaitmap_begin,    mpiwaitmap_end,    "MPI_WaitMAP")
      call add_event(mpiwaitdep_begin,    mpiwaitdep_end,    "MPI_WaitDEP")
      call add_event(mpiwaitmg_begin,     mpiwaitmg_end,     "MPI_WaitMG")
      call add_event(mpiwaitfile_begin,   mpiwaitfile_end,   "MPI_WaitFILE")
      call add_event(mpiwaitmapfile_begin,mpiwaitmapfile_end,"MPI_WaitMAPFILE")
      call add_event(p1_begin,            p1_end,            "Probe1")
      call add_event(p2_begin,            p2_end,            "Probe2")
      call add_event(p3_begin,            p3_end,            "Probe3")
      call add_event(p4_begin,            p4_end,            "Probe4")
      call add_event(p5_begin,            p5_end,            "Probe5")
      call add_event(p6_begin,            p6_end,            "Probe6")
      
   end subroutine log_setup

   subroutine add_event(e_begin, e_end, e_name)
      integer, intent(out) :: e_begin, e_end
      character(len=*), intent(in) :: e_name
      integer, save :: e_counter = 0

      e_counter = e_counter + 1

      if ( e_counter > nevents ) then
         write(6,*) "ERROR: nevents is incorrectly specified"
         write(6,*) e_counter, nevents
         stop
      end if

      e_begin = e_counter
      e_end = e_counter
      event_name(e_counter) = e_name
        
   end subroutine add_event

#ifdef simple_timer
   subroutine simple_timer_finalize
      ! Calculate the mean, min and max times for each case
      integer :: i
      integer(kind=4) :: ierr, llen, lcomm
      real(kind=8), dimension(nevents) :: emean, emax, emin
      real(kind=8), dimension(nevents) :: tot_time_r
      real, dimension(2) :: time_l, time_mean, time_max, time_min
      
      llen = nevents
      lcomm = comm_world
      emean = 0._8
      emax = 0._8
      emin = 0._8
      tot_time_r = real(tot_time,8)/1000000.
      call MPI_Reduce(tot_time_r, emean, llen, MPI_DOUBLE_PRECISION, &
                      MPI_SUM, 0_4, lcomm, ierr )
      call MPI_Reduce(tot_time_r, emax, llen, MPI_DOUBLE_PRECISION,  &
                      MPI_MAX, 0_4, lcomm, ierr )
      call MPI_Reduce(tot_time_r, emin, llen, MPI_DOUBLE_PRECISION,  &
                      MPI_MIN, 0_4, lcomm, ierr )
      if ( myid == 0 ) then
         write(6,*) "==============================================="
         write(6,*) "  Times over all processes"
         write(6,*) "  Routine        Mean time  Min time  Max time"
         do i = 1,nevents
            if ( emean(i) > 0. ) then
               ! This stops boundsa, b getting written when they're not used.
               write(*,"(a,3f10.3)") event_name(i), emean(i)/nproc, emin(i), emax(i)
            end if
         end do
      end if
      
      llen = 2
      time_mean = 0.
      time_max = 0.
      time_min = 0.
      time_l(1:2) = (/ mpiinit_time, total_time /)
      call MPI_Reduce(time_l, time_mean, llen, MPI_REAL, &
                      MPI_SUM, 0_4, lcomm, ierr )
      time_mean = time_mean/real(nproc)
      call MPI_Reduce(time_l, time_max, llen, MPI_REAL,  &
                      MPI_MAX, 0_4, lcomm, ierr )
      call MPI_Reduce(time_l, time_min, llen, MPI_REAL,  &
                      MPI_MIN, 0_4, lcomm, ierr )
      if ( myid == 0 ) then
         write(*,"(a,3f10.3)") "MPI_Initialise ",time_mean(1),time_min(1),time_max(1)
         write(*,"(a,3f10.3)") "Total_Time     ",time_mean(2),time_min(2),time_max(2)
      end if   
        
   end subroutine simple_timer_finalize
#endif

   ! Read and distribute a global variable
   ! Optional arguments for format and to skip over records
   subroutine readglobvar2(un,var,skip,fmt)
      integer, intent(in) :: un
      real, dimension(:), intent(out) :: var
      logical, intent(in), optional :: skip
      character(len=*), intent(in), optional :: fmt
      real, dimension(ifull_g) :: varg
      logical :: doskip

      doskip = .false.
      if ( present(skip) ) doskip = skip

      if ( doskip ) then
         if ( myid == 0 ) then
            read(un)
         end if
      else
         if ( myid == 0 ) then
            if ( present(fmt) ) then
               if ( fmt == "*" ) then
                  read(un,*) varg
               else
                  read(un,fmt) varg
               end if
            else
               read(un) varg
            end if
            ! Use explicit ranges here because some arguments might be extended.
            call ccmpi_distribute(var(1:ifull),varg)
         else
            call ccmpi_distribute(var(1:ifull))
         end if
      end if

   end subroutine readglobvar2

   subroutine readglobvar2i(un,var,skip,fmt)
      integer, intent(in) :: un
      integer, dimension(:), intent(out) :: var
      logical, intent(in), optional :: skip
      character(len=*), intent(in), optional :: fmt
      integer, dimension(ifull_g) :: varg
      logical :: doskip

      doskip = .false.
      if ( present(skip) ) doskip = skip

      if ( doskip ) then
         if ( myid == 0 ) then
            read(un)
         end if
      else
         if ( myid == 0 ) then
            if ( present(fmt) ) then
               if ( fmt == "*" ) then
                  read(un,*) varg
               else
                  read(un,fmt) varg
               end if
            else
               read(un) varg
            end if
            ! Use explicit ranges here because some arguments might be extended.
            call ccmpi_distribute(var(1:ifull),varg)
         else
            call ccmpi_distribute(var(1:ifull))
         end if
      end if

   end subroutine readglobvar2i

   subroutine readglobvar3(un,var,skip,fmt)
      integer, intent(in) :: un
      real, dimension(:,:), intent(out) :: var
      logical, intent(in), optional :: skip
      real, dimension(ifull_g,size(var,2)) :: varg
      character(len=*), intent(in), optional :: fmt
      integer :: k, kk
      logical :: doskip

      doskip = .false.
      if ( present(skip) ) doskip = skip

      kk = size(var,2)

      if ( doskip ) then
         if ( myid == 0 ) then
            read(un)
         end if
      else
         if ( myid == 0 ) then
            if ( present(fmt) ) then
               if ( fmt == "*" ) then
                  read(un,*) varg
               else
                  read(un,fmt) varg
               end if
            else
               read(un) varg
            end if
            call ccmpi_distribute(var(1:ifull,:),varg)
         else
            call ccmpi_distribute(var(1:ifull,:))
         end if
      end if

   end subroutine readglobvar3

   ! Gather and write a global variable
   ! Optional argument for format
   subroutine writeglobvar2(un,var,fmt)
      integer, intent(in) :: un
      real, dimension(:), intent(in) :: var
      character(len=*), intent(in), optional :: fmt
      real, dimension(ifull_g) :: varg

      if ( myid == 0 ) then
         ! Use explicit ranges here because some arguments might be extended.
         call ccmpi_gather(var(1:ifull),varg)
         if ( present(fmt) ) then
            if ( fmt == "*" ) then
               write(un,*) varg
            else
               write(un,fmt) varg
            end if
         else
            write(un) varg
         end if
      else
         call ccmpi_gather(var(1:ifull))
      end if

   end subroutine writeglobvar2

   subroutine writeglobvar3(un,var,fmt)
      integer, intent(in) :: un
      real, dimension(:,:), intent(in) :: var
      real, dimension(ifull_g,size(var,2)) :: varg
      character(len=*), intent(in), optional :: fmt
      integer :: k, kk

      kk = size(var,2)

      if ( myid == 0 ) then
         call ccmpi_gather(var(1:ifull,:),varg)
         if ( present(fmt) ) then
            if ( fmt == "*" ) then
               write(un,*) varg
            else
               write(un,fmt) varg
            end if
         else
            write(un) varg
         end if
      else
         call ccmpi_gather(var(1:ifull,:))
      end if

   end subroutine writeglobvar3

   subroutine ccmpi_reduce2i(ldat,gdat,op,host,comm)
   
      integer, intent(in) :: host,comm
      integer(kind=4) :: lop, lcomm, ierr, lsize, lhost
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_INTEGER8
#else
      integer(kind=4), parameter :: ltype = MPI_INTEGER
#endif 
      integer, dimension(:), intent(in) :: ldat
      integer, dimension(:), intent(out) :: gdat
      character(len=*), intent(in) :: op

      lhost = host
      lcomm = comm
      lsize = size(ldat)
      
      select case( op )
         case( "max" )
            lop = MPI_MAX
         case( "min" )
            lop = MPI_MIN
         case( "sum" )
            lop = MPI_SUM
         case default
            write(6,*) "ERROR: Unknown option for ccmpi_reduce ",op
            call ccmpi_abort(-1)
      end select
      
      call START_LOG(reduce_begin)
      call MPI_Reduce(ldat, gdat, lsize, ltype, lop, lhost, lcomm, ierr )
      call END_LOG(reduce_end)
   
   end subroutine ccmpi_reduce2i

   subroutine ccmpi_reduce1r(ldat,gdat,op,host,comm)
   
      integer, intent(in) :: host, comm
      integer(kind=4) :: ltype, lop, lcomm, lerr, lhost
      real, intent(in) :: ldat
      real, intent(out) :: gdat
      character(len=*), intent(in) :: op

      lhost = host
      lcomm = comm
      
      select case( op )
         case( "max" )
            lop = MPI_MAX
#ifdef i8r8
            ltype = MPI_DOUBLE_PRECISION
#else
            ltype = MPI_REAL
#endif 
         case( "min" )
            lop = MPI_MIN
#ifdef i8r8
            ltype = MPI_DOUBLE_PRECISION
#else
            ltype = MPI_REAL
#endif 
         case( "sum" )
            lop = MPI_SUM
#ifdef i8r8
            ltype = MPI_DOUBLE_PRECISION
#else
            ltype = MPI_REAL
#endif 
         case default
            write(6,*) "ERROR: Unknown option for ccmpi_reduce ",op
            call ccmpi_abort(-1)
      end select

      call START_LOG(reduce_begin)
      call MPI_Reduce(ldat, gdat, 1_4, ltype, lop, lhost, lcomm, lerr )
      call END_LOG(reduce_end)
   
   end subroutine ccmpi_reduce1r
   
   subroutine ccmpi_reduce2r(ldat,gdat,op,host,comm)
   
      integer, intent(in) :: host, comm
      integer(kind=4) :: ltype, lop, lcomm, lerr, lsize, lhost
      real, dimension(:), intent(in) :: ldat
      real, dimension(:), intent(out) :: gdat
      character(len=*), intent(in) :: op
      
      lhost = host
      lcomm = comm
      lsize = size(ldat)
            
      select case( op )
         case( "max" )
            lop = MPI_MAX
#ifdef i8r8
            ltype = MPI_DOUBLE_PRECISION
#else
            ltype = MPI_REAL
#endif 
         case( "min" )
            lop = MPI_MIN
#ifdef i8r8
            ltype = MPI_DOUBLE_PRECISION
#else
            ltype = MPI_REAL
#endif 
         case( "sum" )
            lop = MPI_SUM
#ifdef i8r8
            ltype = MPI_DOUBLE_PRECISION
#else
            ltype = MPI_REAL
#endif 
         case( "maxloc" )
            lop = MPI_MAXLOC
            lsize = lsize/2
#ifdef i8r8
            ltype = MPI_2DOUBLE_PRECISION
#else
            ltype = MPI_2REAL
#endif 
         case( "minloc" )
            lop = MPI_MINLOC
            lsize = lsize/2
#ifdef i8r8
            ltype = MPI_2DOUBLE_PRECISION
#else
            ltype = MPI_2REAL
#endif 
         case default
            write(6,*) "ERROR: Unknown option for ccmpi_reduce ",op
            call ccmpi_abort(-1)
      end select

      call START_LOG(reduce_begin)
      call MPI_Reduce(ldat, gdat, lsize, ltype, lop, lhost, lcomm, lerr )
      call END_LOG(reduce_end)
   
   end subroutine ccmpi_reduce2r

   subroutine ccmpi_reduce3r(ldat,gdat,op,host,comm)
   
      integer, intent(in) :: host, comm
      integer(kind=4) ltype, lop, lcomm, lerr, lsize, lhost
      real, dimension(:,:), intent(in) :: ldat
      real, dimension(:,:), intent(out) :: gdat
      character(len=*), intent(in) :: op

      lhost = host
      lcomm = comm
      lsize = size(ldat)
      
      select case( op )
         case( "max" )
            lop = MPI_MAX
#ifdef i8r8
            ltype = MPI_DOUBLE_PRECISION
#else
            ltype = MPI_REAL
#endif 
         case( "min" )
            lop = MPI_MIN
#ifdef i8r8
            ltype = MPI_DOUBLE_PRECISION
#else
            ltype = MPI_REAL
#endif 
         case( "sum" )
            lop = MPI_SUM
#ifdef i8r8
            ltype = MPI_DOUBLE_PRECISION
#else
            ltype = MPI_REAL
#endif 
         case( "maxloc" )
            lop = MPI_MAXLOC
            lsize = lsize/2
#ifdef i8r8
            ltype = MPI_2DOUBLE_PRECISION
#else
            ltype = MPI_2REAL
#endif 
         case( "minloc" )
            lop = MPI_MINLOC
            lsize = lsize/2
#ifdef i8r8
            ltype = MPI_2DOUBLE_PRECISION
#else
            ltype = MPI_2REAL
#endif 
         case default
            write(6,*) "ERROR: Unknown option for ccmpi_reduce ",op
            call ccmpi_abort(-1)

      end select
      
      call START_LOG(reduce_begin)
      call MPI_Reduce(ldat, gdat, lsize, ltype, lop, lhost, lcomm, lerr )
      call END_LOG(reduce_end)
   
   end subroutine ccmpi_reduce3r

   subroutine ccmpi_reduce2c(ldat,gdat,op,host,comm)
   
      use sumdd_m
   
      integer, intent(in) :: host,comm
      integer(kind=4) :: lop, lcomm, lerr, lsize, lhost
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_COMPLEX
#else
      integer(kind=4), parameter :: ltype = MPI_COMPLEX
#endif 
      complex, dimension(:), intent(in) :: ldat
      complex, dimension(:), intent(out) :: gdat
      character(len=*), intent(in) :: op

      lhost = host
      lcomm = comm
      lsize = size(ldat)
      
      select case( op )
         case( "max" )
            lop = MPI_MAX
         case( "min" )
            lop = MPI_MIN
         case( "sum" )
            lop = MPI_SUM
         case( "sumdr" )
            lop = MPI_SUMDR
         case default
            write(6,*) "ERROR: Unknown option for ccmpi_reduce ",op
            call ccmpi_abort(-1)
      end select
      
      call START_LOG(reduce_begin)
      call MPI_Reduce(ldat, gdat, lsize, ltype, lop, lhost, lcomm, lerr )
      call END_LOG(reduce_end)
   
   end subroutine ccmpi_reduce2c

   subroutine ccmpi_reduce2l(ldat,gdat,op,host,comm)
   
      integer, intent(in) :: host,comm
      integer(kind=4) :: lop, lcomm, ierr, lsize, lhost
      integer(kind=4), parameter :: ltype = MPI_LOGICAL
      logical, dimension(:), intent(in) :: ldat
      logical, dimension(:), intent(out) :: gdat
      character(len=*), intent(in) :: op

      lhost = host
      lcomm = comm
      lsize = size(ldat)
      
      select case( op )
         case( "or" )
            lop = MPI_LOR
         case( "and" )
            lop = MPI_LAND
         case default
            write(6,*) "ERROR: Unknown option for ccmpi_reduce ",op
            call ccmpi_abort(-1)
      end select
      
      call START_LOG(reduce_begin)
      call MPI_Reduce(ldat, gdat, lsize, ltype, lop, lhost, lcomm, ierr )
      call END_LOG(reduce_end)
   
   end subroutine ccmpi_reduce2l

   subroutine ccmpi_reduce1rr8(ldat,gdat,op,host,comm)
   
      integer, intent(in) :: host, comm
      integer(kind=4) :: ltype, lop, lcomm, lerr, lhost
      real(kind=8), intent(in) :: ldat
      real(kind=8), intent(out) :: gdat
      character(len=*), intent(in) :: op
      
      lhost = host
      lcomm = comm
            
      select case( op )
         case( "max" )
            lop = MPI_MAX
            ltype = MPI_DOUBLE_PRECISION
         case( "min" )
            lop = MPI_MIN
            ltype = MPI_DOUBLE_PRECISION
         case( "sum" )
            lop = MPI_SUM
            ltype = MPI_DOUBLE_PRECISION
         case default
            write(6,*) "ERROR: Unknown option for ccmpi_reduce ",op
            call ccmpi_abort(-1)
      end select

      call START_LOG(reduce_begin)
      call MPI_Reduce(ldat, gdat, 1_4, ltype, lop, lhost, lcomm, lerr )
      call END_LOG(reduce_end)
   
   end subroutine ccmpi_reduce1rr8
   
   subroutine ccmpi_allreduce1i(ldat,gdat,op,comm)
   
      integer, intent(in) :: comm
      integer(kind=4) :: lop, lcomm, lerr
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_INTEGER8
#else
      integer(kind=4), parameter :: ltype = MPI_INTEGER
#endif 
      integer, intent(in) :: ldat
      integer, intent(out) :: gdat
      character(len=*), intent(in) :: op

      lcomm = comm
      
      select case( op )
         case( "max" )
            lop = MPI_MAX
         case( "min" )
            lop = MPI_MIN
         case( "sum" )
            lop = MPI_SUM
         case default
            write(6,*) "ERROR: Unknown option for ccmpi_allreduce ",op
            call ccmpi_abort(-1)
      end select
      
      call START_LOG(allreduce_begin)
      call MPI_AllReduce(ldat, gdat, 1_4, ltype, lop, lcomm, lerr )
      call END_LOG(allreduce_end)
   
   end subroutine ccmpi_allreduce1i
   
   subroutine ccmpi_allreduce2i(ldat,gdat,op,comm)
   
      integer, intent(in) :: comm
      integer(kind=4) :: lop, lcomm, lerr, lsize
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_INTEGER8
#else
      integer(kind=4), parameter :: ltype = MPI_INTEGER
#endif 
      integer, dimension(:), intent(in) :: ldat
      integer, dimension(:), intent(out) :: gdat
      character(len=*), intent(in) :: op

      lcomm = comm
      lsize = size(ldat)
      
      select case( op )
         case( "max" )
            lop = MPI_MAX
         case( "min" )
            lop = MPI_MIN
         case( "sum" )
            lop = MPI_SUM
         case default
            write(6,*) "ERROR: Unknown option for ccmpi_allreduce ",op
            call ccmpi_abort(-1)
      end select
      
      call START_LOG(allreduce_begin)
      call MPI_AllReduce(ldat, gdat, lsize, ltype, lop, lcomm, lerr )
      call END_LOG(allreduce_end)
   
   end subroutine ccmpi_allreduce2i

   subroutine ccmpi_allreduce1r(ldat,gdat,op,comm)
   
      integer, intent(in) :: comm
      integer(kind=4) lop, lcomm, lerr
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      real, intent(in) :: ldat
      real, intent(out) :: gdat
      character(len=*), intent(in) :: op

      lcomm = comm
      
      select case( op )
         case( "max" )
            lop = MPI_MAX
         case( "min" )
            lop = MPI_MIN
         case( "sum" )
            lop = MPI_SUM
         case default
            write(6,*) "ERROR: Unknown option for ccmpi_allreduce ",op
            call ccmpi_abort(-1)
      end select

      call START_LOG(allreduce_begin)
      call MPI_AllReduce(ldat, gdat, 1_4, ltype, lop, lcomm, lerr )
      call END_LOG(allreduce_end)
   
   end subroutine ccmpi_allreduce1r
   
   subroutine ccmpi_allreduce2r(ldat,gdat,op,comm)
   
      integer, intent(in) :: comm
      integer(kind=4) lop, lcomm, lerr, lsize
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      real, dimension(:), intent(in) :: ldat
      real, dimension(:), intent(out) :: gdat
      character(len=*), intent(in) :: op

      lcomm = comm
      lsize = size(ldat)
      
      select case( op )
         case( "max" )
            lop = MPI_MAX
         case( "min" )
            lop = MPI_MIN
         case( "sum" )
            lop = MPI_SUM
         case default
            write(6,*) "ERROR: Unknown option for ccmpi_allreduce ",op
            call ccmpi_abort(-1)
      end select

      call START_LOG(allreduce_begin)
      call MPI_AllReduce(ldat, gdat, lsize, ltype, lop, lcomm, lerr )
      call END_LOG(allreduce_end)
   
   end subroutine ccmpi_allreduce2r
  
   subroutine ccmpi_allreduce3r(ldat,gdat,op,comm)
   
      use sumdd_m
   
      integer, intent(in) :: comm
      integer(kind=4) ltype, lop, lcomm, lerr, lsize
      real, dimension(:,:), intent(in) :: ldat
      real, dimension(:,:), intent(out) :: gdat
      character(len=*), intent(in) :: op

      lcomm = comm
      lsize = size(ldat)
      
      select case( op )
         case( "max" )
            lop = MPI_MAX
#ifdef i8r8
            ltype = MPI_DOUBLE_PRECISION
#else
            ltype = MPI_REAL
#endif 
         case( "min" )
            lop = MPI_MIN
#ifdef i8r8
            ltype = MPI_DOUBLE_PRECISION
#else
            ltype = MPI_REAL
#endif 
         case( "sum" )
            lop = MPI_SUM
#ifdef i8r8
            ltype = MPI_DOUBLE_PRECISION
#else
            ltype = MPI_REAL
#endif 
         case( "maxloc" )
            lop = MPI_MAXLOC
            lsize = lsize/2
#ifdef i8r8
            ltype = MPI_2DOUBLE_PRECISION
#else
            ltype = MPI_2REAL
#endif 
         case( "minloc" )
            lop = MPI_MINLOC
            lsize = lsize/2
#ifdef i8r8
            ltype = MPI_2DOUBLE_PRECISION
#else
            ltype = MPI_2REAL
#endif 
         case default
            write(6,*) "ERROR: Unknown option for ccmpi_allreduce ",op
            call ccmpi_abort(-1)
      end select
         
      call START_LOG(allreduce_begin)
      call MPI_AllReduce(ldat, gdat, lsize, ltype, lop, lcomm, lerr )
      call END_LOG(allreduce_end)
   
   end subroutine ccmpi_allreduce3r

   subroutine ccmpi_allreduce1c(ldat,gdat,op,comm)
   
      use sumdd_m
   
      integer, intent(in) :: comm
      integer(kind=4) :: lop, lcomm, lerr
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_COMPLEX
#else
      integer(kind=4), parameter :: ltype = MPI_COMPLEX
#endif 
      complex, intent(in) :: ldat
      complex, intent(out) :: gdat
      character(len=*), intent(in) :: op

      lcomm = comm
      
      select case( op )
         case( "max" )
            lop = MPI_MAX
         case( "min" )
            lop = MPI_MIN
         case( "sum" )
            lop = MPI_SUM
         case( "sumdr" )
            lop = MPI_SUMDR
         case default
            write(6,*) "ERROR: Unknown option for ccmpi_allreduce ",op
            call ccmpi_abort(-1)
      end select
      
      call START_LOG(allreduce_begin)
      call MPI_AllReduce(ldat, gdat, 1_4, ltype, lop, lcomm, lerr )
      call END_LOG(allreduce_end)
   
   end subroutine ccmpi_allreduce1c
   
   subroutine ccmpi_allreduce2c(ldat,gdat,op,comm)
   
      use sumdd_m
   
      integer, intent(in) :: comm
      integer(kind=4) :: lop, lcomm, lerr, lsize
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_COMPLEX
#else
      integer(kind=4), parameter :: ltype = MPI_COMPLEX
#endif 
      complex, dimension(:), intent(in) :: ldat
      complex, dimension(:), intent(out) :: gdat
      character(len=*), intent(in) :: op

      lcomm = comm
      lsize = size(ldat)
      
      select case( op )
         case( "max" )
            lop = MPI_MAX
         case( "min" )
            lop = MPI_MIN
         case( "sum" )
            lop = MPI_SUM
         case( "sumdr" )
            lop = MPI_SUMDR
         case default
            write(6,*) "ERROR: Unknown option for ccmpi_allreduce ",op
            call ccmpi_abort(-1)
      end select
      
      call START_LOG(allreduce_begin)
      call MPI_AllReduce(ldat, gdat, lsize, ltype, lop, lcomm, lerr )
      call END_LOG(allreduce_end)
   
   end subroutine ccmpi_allreduce2c
   
   subroutine ccmpi_abort(ierrin)
   
      integer, intent(in) :: ierrin
      integer(kind=4) :: lerrin, ierr
      
      if ( myid==0 ) then
        call finishbanner
      end if
      lerrin = ierrin
      call MPI_Abort(MPI_COMM_WORLD, lerrin ,ierr)
   
   end subroutine ccmpi_abort

   subroutine ccmpi_bcast1i(ldat,host,comm)

      integer, intent(in) :: host, comm
      integer(kind=4) :: lcomm, lhost, lerr
#ifdef i8r8
      integer(kind=4) :: ltype = MPI_INTEGER8
#else
      integer(kind=4) :: ltype = MPI_INTEGER
#endif
      integer, intent(inout) :: ldat

      lhost = host
      lcomm = comm
      call START_LOG(bcast_begin)
      call MPI_Bcast(ldat,1_4,ltype,lhost,lcomm,lerr)
      call END_LOG(bcast_end)
         
   end subroutine ccmpi_bcast1i

   subroutine ccmpi_bcast2i(ldat,host,comm)

      integer, intent(in) :: host, comm
      integer(kind=4) :: lcomm, lhost, lerr, lsize
#ifdef i8r8
      integer(kind=4) :: ltype = MPI_INTEGER8
#else
      integer(kind=4) :: ltype = MPI_INTEGER
#endif
      integer, dimension(:), intent(inout) :: ldat

      lhost = host
      lcomm = comm
      lsize = size(ldat)
      call START_LOG(bcast_begin)
      call MPI_Bcast(ldat,lsize,ltype,lhost,lcomm,lerr)
      call END_LOG(bcast_end)
         
   end subroutine ccmpi_bcast2i

   subroutine ccmpi_bcast3i(ldat,host,comm)

      integer, intent(in) :: host, comm
      integer(kind=4) :: lcomm, lhost, lerr, lsize
#ifdef i8r8
      integer(kind=4) :: ltype = MPI_INTEGER8
#else
      integer(kind=4) :: ltype = MPI_INTEGER
#endif
      integer, dimension(:,:), intent(inout) :: ldat

      lhost = host
      lcomm = comm
      lsize = size(ldat)
      call START_LOG(bcast_begin)
      call MPI_Bcast(ldat,lsize,ltype,lhost,lcomm,lerr)
      call END_LOG(bcast_end)
      
   end subroutine ccmpi_bcast3i

   subroutine ccmpi_bcast1r(ldat,host,comm)
   
      integer, intent(in) :: host, comm
      integer(kind=4) :: lcomm, lhost, lerr
#ifdef i8r8
      integer(kind=4) :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4) :: ltype = MPI_REAL
#endif
      real, intent(inout) :: ldat

      lhost = host
      lcomm = comm
      call START_LOG(bcast_begin)
      call MPI_Bcast(ldat,1_4,ltype,lhost,lcomm,lerr)
      call END_LOG(bcast_end)
   
   end subroutine ccmpi_bcast1r
   
   subroutine ccmpi_bcast2r(ldat,host,comm)
   
      integer, intent(in) :: host, comm
      integer(kind=4) :: lcomm, lhost, lerr, lsize
#ifdef i8r8
      integer(kind=4) :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4) :: ltype = MPI_REAL
#endif
      real, dimension(:), intent(inout) :: ldat

      lhost = host
      lcomm = comm
      lsize = size(ldat)
      call START_LOG(bcast_begin)
      call MPI_Bcast(ldat,lsize,ltype,lhost,lcomm,lerr)
      call END_LOG(bcast_end)
   
   end subroutine ccmpi_bcast2r

   subroutine ccmpi_bcast3r(ldat,host,comm)
   
      integer, intent(in) :: host, comm
      integer(kind=4) :: lcomm, lhost, lerr, lsize
#ifdef i8r8
      integer(kind=4) :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4) :: ltype = MPI_REAL
#endif
      real, dimension(:,:), intent(inout) :: ldat

      lhost = host
      lcomm = comm
      lsize = size(ldat)
      call START_LOG(bcast_begin)
      call MPI_Bcast(ldat,lsize,ltype,lhost,lcomm,lerr)
      call END_LOG(bcast_end)
   
   end subroutine ccmpi_bcast3r

   subroutine ccmpi_bcast4r(ldat,host,comm)
   
      integer, intent(in) :: host, comm
      integer(kind=4) :: lcomm, lhost, lerr, lsize
#ifdef i8r8
      integer(kind=4) :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4) :: ltype = MPI_REAL
#endif
      real, dimension(:,:,:), intent(inout) :: ldat

      lhost = host
      lcomm = comm
      lsize = size(ldat)
      call START_LOG(bcast_begin)
      call MPI_Bcast(ldat, lsize, ltype, lhost, lcomm, lerr)
      call END_LOG(bcast_end)
   
   end subroutine ccmpi_bcast4r

   subroutine ccmpi_bcast5r(ldat,host,comm)
   
      integer, intent(in) :: host, comm
      integer(kind=4) :: lcomm, lhost, lerr, lsize
#ifdef i8r8
      integer(kind=4) :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4) :: ltype = MPI_REAL
#endif
      real, dimension(:,:,:,:), intent(inout) :: ldat

      lhost = host
      lcomm = comm
      lsize = size(ldat)
      call START_LOG(bcast_begin)
      call MPI_Bcast(ldat,lsize,ltype,lhost,lcomm,lerr)
      call END_LOG(bcast_end)
   
   end subroutine ccmpi_bcast5r

   subroutine ccmpi_bcast1s(ldat,host,comm)
   
      integer, intent(in) :: host, comm
      integer(kind=4) :: lcomm, lhost, lerr, lsize
      character(len=*), intent(inout) :: ldat
      integer, parameter :: maxdummysize = 1024
      integer i
      integer(kind=1), dimension(maxdummysize) :: dummy

      ! MJT notes - MS Windows MPI_CHARACTER seems broken

      lhost = host
      lcomm = comm
      lsize = len(ldat)
      if ( lsize > maxdummysize ) then
        write(6,*) "ERROR: Dummy array too small in ccmpi_bcast1s"
        call ccmpi_abort(-1)
      end if
      do i = 1,lsize
         dummy(i) = int(iachar(ldat(i:i)),1)
      end do
      call START_LOG(bcast_begin)
      call MPI_Bcast(dummy,lsize,MPI_BYTE,lhost,lcomm,lerr)
      call END_LOG(bcast_end)
      do i = 1,lsize
         ldat(i:i) = achar(dummy(i))
      end do
   
   end subroutine ccmpi_bcast1s
   
   subroutine ccmpi_bcast2r8(ldat,host,comm)
   
      integer, intent(in) :: host, comm
      integer(kind=4) :: lcomm, lhost, ierr, lsize
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
      real(kind=8), dimension(:), intent(inout) :: ldat
      
      lhost = host
      lcomm = comm
      lsize = size(ldat)
      call START_LOG(bcast_end)
      call MPI_Bcast(ldat,lsize,ltype,lhost,lcomm,ierr)
      call END_LOG(bcast_end)
   
   end subroutine ccmpi_bcast2r8
   
   subroutine ccmpi_bcast3r8(ldat,host,comm)
   
      integer, intent(in) :: host,comm
      integer(kind=4) :: lcomm, lhost, ierr, lsize
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
      real(kind=8), dimension(:,:), intent(inout) :: ldat
   
      lhost = host
      lcomm = comm
      lsize = size(ldat)
      call START_LOG(bcast_begin)
      call MPI_Bcast(ldat, lsize, ltype, lhost, lcomm, ierr)
      call END_LOG(bcast_end)
   
   end subroutine ccmpi_bcast3r8   
   
   subroutine ccmpi_bcast4r8(ldat,host,comm)
   
      integer, intent(in) :: host,comm
      integer(kind=4) :: lcomm, lhost, ierr, lsize
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
      real(kind=8), dimension(:,:,:), intent(inout) :: ldat

      lhost = host
      lcomm = comm
      lsize = size(ldat)
      call START_LOG(bcast_begin)
      call MPI_Bcast(ldat,lsize,ltype,lhost,lcomm,ierr)
      call END_LOG(bcast_end)
   
   end subroutine ccmpi_bcast4r8

   subroutine ccmpi_barrier(comm)
   
      integer, intent(in) :: comm
      integer(kind=4) :: lcomm, ierr
      
      lcomm = comm
      call MPI_Barrier( lcomm, ierr )
   
   end subroutine ccmpi_barrier
   
   subroutine ccmpi_gatherx2r(gdat,ldat,host,comm)

      integer, intent(in) :: host, comm
      integer(kind=4) :: lsize, lhost, lcomm, lerr
#ifdef i8r8
      integer(kind=4) :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4) :: ltype = MPI_REAL
#endif
      real, dimension(:), intent(out) :: gdat
      real, dimension(:), intent(in) :: ldat

      lcomm = comm
      lhost = host
      lsize = size(ldat)
      call START_LOG(gather_begin)
      call MPI_Gather(ldat,lsize,ltype,gdat,lsize,ltype,lhost,lcomm,lerr)
      call END_LOG(gather_end)
      
   end subroutine ccmpi_gatherx2r
   
   subroutine ccmpi_gatherx3r(gdat,ldat,host,comm)

      integer, intent(in) :: host, comm
      integer(kind=4) :: lsize, lhost, lcomm, lerr
#ifdef i8r8
      integer(kind=4) :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4) :: ltype = MPI_REAL
#endif
      real, dimension(:,:), intent(out) :: gdat
      real, dimension(:,:), intent(in) :: ldat

      lcomm = comm
      lhost = host
      lsize = size(ldat)
      call START_LOG(gather_begin)
      call MPI_Gather(ldat,lsize,ltype,gdat,lsize,ltype,lhost,lcomm,lerr)
      call END_LOG(gather_end)
      
   end subroutine ccmpi_gatherx3r

    subroutine ccmpi_gatherx4r(gdat,ldat,host,comm)

      integer, intent(in) :: host, comm
      integer(kind=4) :: lsize, lhost, lcomm, lerr
#ifdef i8r8
      integer(kind=4) :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4) :: ltype = MPI_REAL
#endif
      real, dimension(:,:,:), intent(out) :: gdat
      real, dimension(:,:,:), intent(in) :: ldat
      
      lcomm = comm
      lhost = host
      lsize = size(ldat)
      call START_LOG(gather_begin)
      call MPI_Gather(ldat,lsize,ltype,gdat,lsize,ltype,lhost,lcomm,lerr)
      call END_LOG(gather_end)
      
   end subroutine ccmpi_gatherx4r
  
   subroutine ccmpi_gatherx23r(gdat,ldat,host,comm)

      integer, intent(in) :: host, comm
      integer(kind=4) :: lsize, lhost, lcomm, lerr
#ifdef i8r8
      integer(kind=4) :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4) :: ltype = MPI_REAL
#endif
      real, dimension(:,:), intent(out) :: gdat
      real, dimension(:), intent(in) :: ldat

      lcomm = comm
      lhost = host
      lsize = size(ldat)
      call START_LOG(gather_begin)
      call MPI_Gather(ldat,lsize,ltype,gdat,lsize,ltype,lhost,lcomm,lerr)
      call END_LOG(gather_end)
      
   end subroutine ccmpi_gatherx23r
   
   subroutine ccmpi_gatherx34r(gdat,ldat,host,comm)

      integer, intent(in) :: host, comm
      integer(kind=4) :: lsize, lhost, lcomm, lerr
#ifdef i8r8
      integer(kind=4) :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4) :: ltype = MPI_REAL
#endif
      real, dimension(:,:,:), intent(out) :: gdat
      real, dimension(:,:), intent(in) :: ldat

      lcomm = comm
      lhost = host
      lsize = size(ldat)
      call START_LOG(gather_begin)
      call MPI_Gather(ldat,lsize,ltype,gdat,lsize,ltype,lhost,lcomm,lerr)
      call END_LOG(gather_end)
      
   end subroutine ccmpi_gatherx34r

   subroutine ccmpi_gatherx45r(gdat,ldat,host,comm)

      integer, intent(in) :: host, comm
      integer(kind=4) :: lsize, lhost, lcomm, lerr
#ifdef i8r8
      integer(kind=4) :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4) :: ltype = MPI_REAL
#endif
      real, dimension(:,:,:,:), intent(out) :: gdat
      real, dimension(:,:,:), intent(in) :: ldat

      lcomm = comm
      lhost = host
      lsize = size(ldat)
      call START_LOG(gather_begin)
      call MPI_Gather(ldat,lsize,ltype,gdat,lsize,ltype,lhost,lcomm,lerr)
      call END_LOG(gather_end)
      
   end subroutine ccmpi_gatherx45r
   
   subroutine ccmpi_gatherx2i(gdat,ldat,host,comm)

      integer, intent(in) :: host, comm
      integer(kind=4) :: lsize, lhost, lcomm, lerr
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_INTEGER8
#else
      integer(kind=4), parameter :: ltype = MPI_INTEGER
#endif 
      integer, dimension(:), intent(out) :: gdat
      integer, dimension(:), intent(in) :: ldat

      lcomm = comm
      lhost = host
      lsize = size(ldat)
      call START_LOG(gather_begin)
      call MPI_Gather(ldat,lsize,ltype,gdat,lsize,ltype,lhost,lcomm,lerr)
      call END_LOG(gather_end)
      
   end subroutine ccmpi_gatherx2i
   
   subroutine ccmpi_gatherx3i(gdat,ldat,host,comm)

      integer, intent(in) :: host, comm
      integer(kind=4) :: lsize, lhost, lcomm, lerr
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_INTEGER8
#else
      integer(kind=4), parameter :: ltype = MPI_INTEGER
#endif 
      integer, dimension(:,:), intent(out) :: gdat
      integer, dimension(:,:), intent(in) :: ldat

      lcomm = comm
      lhost = host
      lsize = size(ldat)
      call START_LOG(gather_begin)
      call MPI_Gather(ldat,lsize,ltype,gdat,lsize,ltype,lhost,lcomm,lerr)
      call END_LOG(gather_end)
      
   end subroutine ccmpi_gatherx3i
   
   subroutine ccmpi_gatherx23l(gdat,ldat,host,comm)

      integer, intent(in) :: host, comm
      integer(kind=4) :: lsize, lhost, lcomm, lerr
      logical, dimension(:,:), intent(out) :: gdat
      logical, dimension(:), intent(in) :: ldat

      lcomm = comm
      lhost = host
      lsize = size(ldat)
      call START_LOG(gather_begin)
      call MPI_Gather(ldat,lsize,MPI_LOGICAL,gdat,lsize,MPI_LOGICAL,lhost,lcomm,lerr)
      call END_LOG(gather_end)
      
   end subroutine ccmpi_gatherx23l
   
   subroutine ccmpi_gatherx2rr8(gdat,ldat,host,comm)

      integer, intent(in) :: host, comm
      integer(kind=4) :: lsize, lhost, lcomm, lerr
      integer(kind=4) :: ltype = MPI_DOUBLE_PRECISION
      real(kind=8), dimension(:), intent(out) :: gdat
      real(kind=8), dimension(:), intent(in) :: ldat

      lcomm = comm
      lhost = host
      lsize = size(ldat)
      call START_LOG(gather_begin)
      call MPI_Gather(ldat,lsize,ltype,gdat,lsize,ltype,lhost,lcomm,lerr)
      call END_LOG(gather_end)
      
   end subroutine ccmpi_gatherx2rr8
   
   subroutine ccmpi_gatherx3rr8(gdat,ldat,host,comm)

      integer, intent(in) :: host, comm
      integer(kind=4) :: lsize, lhost, lcomm, lerr
      integer(kind=4) :: ltype = MPI_DOUBLE_PRECISION
      real(kind=8), dimension(:,:), intent(out) :: gdat
      real(kind=8), dimension(:,:), intent(in) :: ldat

      lcomm = comm
      lhost = host
      lsize = size(ldat)
      call START_LOG(gather_begin)
      call MPI_Gather(ldat,lsize,ltype,gdat,lsize,ltype,lhost,lcomm,lerr)
      call END_LOG(gather_end)
      
   end subroutine ccmpi_gatherx3rr8
   
   subroutine ccmpi_gatherx23rr8(gdat,ldat,host,comm)

      integer, intent(in) :: host, comm
      integer(kind=4) :: lsize, lhost, lcomm, lerr
      integer(kind=4) :: ltype = MPI_DOUBLE_PRECISION
      real(kind=8), dimension(:,:), intent(out) :: gdat
      real(kind=8), dimension(:), intent(in) :: ldat

      lcomm = comm
      lhost = host
      lsize = size(ldat)
      call START_LOG(gather_begin)
      call MPI_Gather(ldat,lsize,ltype,gdat,lsize,ltype,lhost,lcomm,lerr)
      call END_LOG(gather_end)
      
   end subroutine ccmpi_gatherx23rr8
   
   subroutine ccmpi_gatherx34rr8(gdat,ldat,host,comm)

      integer, intent(in) :: host, comm
      integer(kind=4) :: lsize, lhost, lcomm, lerr
      integer(kind=4) :: ltype = MPI_DOUBLE_PRECISION
      real(kind=8), dimension(:,:,:), intent(out) :: gdat
      real(kind=8), dimension(:,:), intent(in) :: ldat
      
      lcomm = comm
      lhost = host
      lsize = size(ldat)
      call START_LOG(gather_begin)
      call MPI_Gather(ldat,lsize,ltype,gdat,lsize,ltype,lhost,lcomm,lerr)
      call END_LOG(gather_end)
      
   end subroutine ccmpi_gatherx34rr8

   subroutine ccmpi_gatherx45rr8(gdat,ldat,host,comm)

      integer, intent(in) :: host, comm
      integer(kind=4) :: lsize, lhost, lcomm, lerr
      integer(kind=4) :: ltype = MPI_DOUBLE_PRECISION
      real(kind=8), dimension(:,:,:,:), intent(out) :: gdat
      real(kind=8), dimension(:,:,:), intent(in) :: ldat

      lcomm = comm
      lhost = host
      lsize = size(ldat)
      call START_LOG(gather_begin)
      call MPI_Gather(ldat,lsize,ltype,gdat,lsize,ltype,lhost,lcomm,lerr)
      call END_LOG(gather_end)
      
   end subroutine ccmpi_gatherx45rr8
   
   subroutine ccmpi_scatterx2r(gdat,ldat,host,comm)
   
      integer, intent(in) :: host, comm
      integer(kind=4) :: lsize, lhost, lcomm, lerr
#ifdef i8r8
      integer(kind=4) :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4) :: ltype = MPI_REAL
#endif
      real, dimension(:), intent(in) :: gdat
      real, dimension(:), intent(out) :: ldat
     
      lcomm = comm
      lhost = host
      lsize = size(ldat)
      call START_LOG(scatter_begin)
      call MPI_Scatter(gdat,lsize,ltype,ldat,lsize,ltype,lhost,lcomm,lerr)
      call END_LOG(scatter_end)
      
   end subroutine ccmpi_scatterx2r

   subroutine ccmpi_scatterx32r(gdat,ldat,host,comm)
   
      integer, intent(in) :: host, comm
      integer(kind=4) :: lsize, lhost, lcomm, lerr
#ifdef i8r8
      integer(kind=4) :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4) :: ltype = MPI_REAL
#endif
      real, dimension(:,:), intent(in) :: gdat
      real, dimension(:), intent(out) :: ldat
      
      lcomm = comm
      lhost = host
      lsize = size(ldat)
      call START_LOG(scatter_begin)
      call MPI_Scatter(gdat,lsize,ltype,ldat,lsize,ltype,lhost,lcomm,lerr)
      call END_LOG(scatter_end)
      
   end subroutine ccmpi_scatterx32r

   subroutine ccmpi_scatterx3r(gdat,ldat,host,comm)
   
      integer, intent(in) :: host, comm
      integer(kind=4) :: lsize, lhost, lcomm, lerr
#ifdef i8r8
      integer(kind=4) :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4) :: ltype = MPI_REAL
#endif
      real, dimension(:,:), intent(in) :: gdat
      real, dimension(:,:), intent(out) :: ldat
      
      lcomm = comm
      lhost = host
      lsize = size(ldat)
      call START_LOG(scatter_begin)
      call MPI_Scatter(gdat,lsize,ltype,ldat,lsize,ltype,lhost,lcomm,lerr)
      call END_LOG(scatter_end)
      
   end subroutine ccmpi_scatterx3r
   
   subroutine ccmpi_scatterx32l(gdat,ldat,host,comm)
   
      integer, intent(in) :: host, comm
      integer(kind=4) :: lsize, lhost, lcomm, lerr
#ifdef i8r8
      integer(kind=4) :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4) :: ltype = MPI_REAL
#endif
      logical, dimension(:,:), intent(in) :: gdat
      logical, dimension(:), intent(out) :: ldat
      
      lcomm = comm
      lhost = host
      lsize = size(ldat)
      call START_LOG(scatter_begin)
      call MPI_Scatter(gdat,lsize,ltype,ldat,lsize,ltype,lhost,lcomm,lerr)
      call END_LOG(scatter_end)
      
   end subroutine ccmpi_scatterx32l
   
   subroutine ccmpi_allgatherx2i(gdat,ldat,comm)
   
      integer, intent(in) :: comm
      integer(kind=4) lsize, lcomm, lerr
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_INTEGER8
#else
      integer(kind=4), parameter :: ltype = MPI_INTEGER
#endif  
      integer, dimension(:), intent(in) :: ldat
      integer, dimension(:), intent(out) :: gdat
      
      lcomm = comm
      lsize = size(ldat)
      call START_LOG(allgather_begin)
      call MPI_AllGather(ldat,lsize,ltype,gdat,lsize,ltype,lcomm,lerr)
      call END_LOG(allgather_end)
      
   end subroutine ccmpi_allgatherx2i

   subroutine ccmpi_allgatherx2r(gdat,ldat,comm)
   
      integer, intent(in) :: comm
      integer(kind=4) lsize, lcomm, lerr
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif  
      real, dimension(:), intent(in) :: ldat
      real, dimension(:), intent(out) :: gdat
   
      lcomm = comm
      lsize = size(ldat)
      call START_LOG(allgather_begin)
      call MPI_AllGather(ldat,lsize,ltype,gdat,lsize,ltype,lcomm,lerr)
      call END_LOG(allgather_end)
      
   end subroutine ccmpi_allgatherx2r
   
   subroutine ccmpi_alltoall2l(gdat,comm)
   
      integer, intent(in) :: comm
      integer(kind=4) :: lsize, lcomm, lerr
      logical, dimension(:), intent(inout) :: gdat
      logical, dimension(size(gdat)) :: rdat
      
      lcomm = comm
      lsize = size(gdat)/nproc
      call START_LOG(alltoall_begin)
      call MPI_AlltoAll( gdat, lsize, MPI_LOGICAL, rdat, lsize, MPI_LOGICAL, lcomm, lerr )
      call END_LOG(alltoall_end)
      gdat = rdat
   
   end subroutine ccmpi_alltoall2l
   
   subroutine ccmpi_commsplit(commout,comm,colour,rank)
   
      integer, intent(out) :: commout
      integer, intent(in) :: comm, colour, rank
      integer(kind=4) :: lcomm, lcommout, lerr, lrank, lcolour
   
      lcomm = comm
      lrank = rank
      if ( colour >= 0 ) then
        lcolour = colour
      else
        lcolour = MPI_UNDEFINED
      end if
      call MPI_Comm_Split(lcomm, lcolour, lrank, lcommout, lerr)
      commout = lcommout
   
   end subroutine ccmpi_commsplit
   
   subroutine ccmpi_commfree(comm)
   
      integer, intent(in) :: comm
      integer(kind=4) :: lcomm, lerr
      
      lcomm = comm
      if ( lcomm /= MPI_COMM_NULL ) then
        call MPI_Comm_Free(lcomm,lerr)
      end if
   
   end subroutine ccmpi_commfree
   
   subroutine ccmpi_init

      integer(kind=4) :: lerr, lproc, lid
      integer, dimension(8) :: times_a, times_b
#ifdef _OPENMP
      integer(kind=4) :: lprovided
#endif

      call date_and_time(values=times_a)

      ! Global communicator
#ifdef _OPENMP
      call MPI_Init_Thread(MPI_THREAD_FUNNELED, lprovided, lerr)
      if ( lprovided < MPI_THREAD_FUNNELED ) then
         write(6,*) "ERROR: MPI does not support MPI_THREAD_FUNNELED"
         call ccmpi_abort(-1)
      end if
#else
      call MPI_Init(lerr)
#endif
      call MPI_Comm_size(MPI_COMM_WORLD, lproc, lerr) ! Find number of processes
      call MPI_Comm_rank(MPI_COMM_WORLD, lid, lerr)   ! Find local processor id
      nproc      = lproc
      myid       = lid
      comm_world = MPI_COMM_WORLD
      
      call date_and_time(values=times_b)
      mpiinit_time = sum( real(times_b(5:8) - times_a(5:8))*(/ 3600., 60., 1., 0.001 /) )
      if ( mpiinit_time < 0. ) mpiinit_time = mpiinit_time + 86400.

   end subroutine ccmpi_init
   
   subroutine ccmpi_reinit(newnproc)
   
      integer, intent(in) :: newnproc
      integer(kind=4) :: lerr, lid
      integer(kind=4) :: lcommout, lcommin
      integer(kind=4) :: lcolour
#ifdef usempi3
      integer(kind=4) :: lproc
#endif
      
      if ( newnproc < nproc ) then
         if ( myid == 0 ) then
            write(6,*) "Reducing number of processes from ",nproc," to ",newnproc  
         end if
         if ( myid < newnproc ) then
            lcolour = 1
         else
            lcolour = MPI_UNDEFINED
         end if
         lcommin = comm_world
         lid = myid
         call MPI_Comm_Split(lcommin, lcolour, lid, lcommout, lerr) ! redefine comm_world
         comm_world = lcommout
         nproc = newnproc
      end if
      
      if ( myid < nproc ) then
   
#ifdef usempi3
         if ( nproc > 1 ) then

            ! Intra-node communicator 
            lid = myid
            lcommin = comm_world
            call MPI_Comm_split_type(lcommin, MPI_COMM_TYPE_SHARED, lid, MPI_INFO_NULL, lcommout, lerr)
            call MPI_Comm_size(lcommout, lproc, lerr) ! Find number of processes on node
            call MPI_Comm_rank(lcommout, lid, lerr)   ! Find local processor id on node
            comm_node  = lcommout
            node_nproc = lproc
            node_myid  = lid

            ! Inter-node commuicator
            lcolour = node_myid
            lid = myid
            lcommin = comm_world
            call MPI_Comm_Split(lcommin, lcolour, lid, lcommout, lerr)
            call MPI_Comm_size(lcommout, lproc, lerr) ! Find number of nodes that have rank node_myid
            call MPI_Comm_rank(lcommout, lid, lerr)   ! Node id for process rank node_myid
            comm_nodecaptain  = lcommout
            nodecaptain_nproc = lproc
            nodecaptain_myid  = lid

            ! Communicate node id 
            lid = nodecaptain_myid
            lcommout = comm_node
            call MPI_Bcast(lid, 1_4, MPI_INTEGER, 0_4, lcommout, lerr)
            node_captainid = lid
            
         else
             
            comm_node  = comm_world
            node_nproc = nproc
            node_myid  = myid
         
            comm_nodecaptain  = comm_world
            nodecaptain_nproc = nproc
            nodecaptain_myid  = myid
         
            node_captainid = myid
            
         end if
      
         if ( myid==0 .and. (node_myid/=0.or.nodecaptain_myid/=0) ) then
            write(6,*) "ERROR: Intra-node communicator failed"
            write(6,*) "myid, node_myid, nodecaptain_myid ",myid,node_myid,nodecaptain_myid
            call ccmpi_abort(-1)
         end if
#else
         ! each process is treated as a node
         lcommin = comm_world
         lcolour = myid
         lid = 0
         call MPI_Comm_Split(lcommin, lcolour, lid, lcommout, lerr)
         comm_node  = lcommout
         node_nproc = 1
         node_myid  = 0
      
         comm_nodecaptain = comm_world
         nodecaptain_nproc = nproc
         nodecaptain_myid = myid
      
         node_captainid = myid
#endif

      else
         call ccmpi_finalize
         stop
      end if

   end subroutine ccmpi_reinit
   
#ifdef usempi3
   subroutine ccmpi_remap
   
      integer :: node_nx, node_ny, node_dx, node_dy
      integer :: oldrank, ty, cy, tx, cx
      integer :: new_node_nproc
      integer :: testid, newid
      integer(kind=4) :: ref_nodecaptain_nproc
      integer(kind=4) :: lerr, lid, lcommin, lcommout
      
      lcommin = comm_world
      node_nx = 0
      node_ny = 0
      node_dx = nxp
      node_dy = nyp
      
      ! communicate number of nodes to all processes
      ref_nodecaptain_nproc = nodecaptain_nproc
      call MPI_Bcast(ref_nodecaptain_nproc, 1_4, MPI_INTEGER, 0_4, lcommin, lerr )
   
      if ( mod(nproc, ref_nodecaptain_nproc)==0 .and. ref_nodecaptain_nproc>1 ) then
         new_node_nproc = nproc / ref_nodecaptain_nproc
         ! calculate virtual node decomposition
         node_nx = max( int(sqrt(real(new_node_nproc))), 1 )
         node_ny = new_node_nproc / node_nx
         node_dx = nxp/node_nx
         node_dy = nyp/node_ny
         do while ( (node_nx*node_dx/=nxp.or.node_ny*node_dy/=nyp.or.node_nx*node_ny/=new_node_nproc) .and. node_nx>0 )
            node_nx = node_nx - 1
            node_ny = new_node_nproc/max( node_nx, 1 )
            node_dx = nxp/max( node_nx, 1 )
            node_dy = nyp/node_ny
         end do   

        ! remap ranks if a valid decomposition has been found
         if ( myid == 0 ) then
            write(6,*) "Remapping ranks using node_nx,node_ny ",node_nx,node_ny
            write(6,*) "node_dx,node_dy                       ",node_dx,node_dy
         end if
         do testid = 0,nproc-1
            oldrank = testid 
            ty = oldrank/(node_ny*node_dx*node_nx) ! node position in y
            oldrank = oldrank - ty*node_ny*node_dx*node_nx
            cy = oldrank/(node_dx*node_nx)         ! y-row position in node
            oldrank = oldrank - cy*node_dx*node_nx
            tx = oldrank/node_nx                   ! node position in x
            oldrank = oldrank - tx*node_nx
            cx = oldrank                           ! x-column position in node
            newid = ty*node_dx*node_nx*node_ny + tx*node_nx*node_ny + cy*node_nx + cx
            if ( newid==myid ) then
               lid = testid
               exit
            end if   
         end do
         lcommin = comm_world
         call MPI_Comm_Split(lcommin, 0_4, lid, lcommout, lerr) ! redefine comm_world
         call MPI_Comm_rank(lcommout, lid, lerr)                ! find local processor id
         comm_world = lcommout
         myid = lid
         if ( lcommin/=MPI_COMM_WORLD .and. lcommin/=MPI_COMM_NULL ) then
            call MPI_Comm_Free(lcommin, lerr) 
         end if
      end if
   
   end subroutine ccmpi_remap
#endif   
   
   subroutine ccmpi_finalize
   
      integer(kind=4) :: lerr
      
      !lcomm = MPI_COMM_WORLD
      !call MPI_Barrier( lcomm, lerr ) ! work around for Pearcey bug
      call MPI_Finalize( lerr )
   
   end subroutine ccmpi_finalize
   
   subroutine ccmpi_procformat_init(localhist,procmode)
   
      integer, intent(inout) :: procmode
      logical, intent(in) :: localhist
#ifdef usempi3
      integer(kind=4) :: lcolour, lcomm, lrank
      integer(kind=4) :: lcommout, lerr, lsize
#endif

      if ( localhist ) then
#ifdef usempi3
         ! configure procmode
         if ( procmode == 0 ) then
            ! procmode=0 uses existing nodes, even if they have different numbers of processes
            if ( myid == 0 ) then
               write(6,*) "Configure procformat output with nodes=",nodecaptain_nproc
            end if
            ! Intra-procmode communicator 
            comm_vnode  = comm_node
            vnode_nproc = node_nproc
            vnode_myid  = node_myid
            ! Inter-procmode communicator
            comm_vleader  = comm_nodecaptain
            vleader_nproc = nodecaptain_nproc
            vleader_myid  = nodecaptain_myid
            ! Communicate procmode id
            vnode_vleaderid = node_captainid
            procmode = node_nproc ! can be different on different nodes
         else
            ! user specified procmode>0 
            procmode = max(procmode, 1) 
            do while ( mod(node_nproc, procmode)/=0 )
               procmode = procmode - 1 ! can be different on different nodes
            end do  
            if ( myid == 0 ) then
               write(6,*) "Configure procformat output with procmode=",procmode
            end if
            ! Intra-procmode communicator
            lcolour = node_myid/procmode
            lcomm = comm_node
            lrank = node_myid
            call MPI_Comm_Split(lcomm, lcolour, lrank, lcommout, lerr)
            comm_vnode = lcommout
            call MPI_Comm_Size(lcommout, lsize, lerr)
            vnode_nproc = lsize
            call MPI_Comm_Rank(lcommout, lrank, lerr)
            vnode_myid = lrank
            ! Inter-procmode communicator
            lcolour = vnode_myid
            lcomm = comm_world
            lrank = myid
            call MPI_Comm_Split(lcomm, lcolour, lrank, lcommout, lerr)
            comm_vleader = lcommout
            call MPI_Comm_Size(lcommout, lsize, lerr)
            vleader_nproc = lsize
            call MPI_Comm_Rank(lcommout, lrank, lerr)
            vleader_myid = lrank
            vnode_vleaderid = vleader_myid
            ! Communicate procmode id
            lcomm = comm_vnode
            lrank = vnode_vleaderid
            call MPI_Bcast( lrank, 1_4, MPI_INTEGER, 0_4, lcomm, lerr )
            vnode_vleaderid = lrank
         end if    
         !call ccmpi_node_leader ! setup comm_vleader and comm_reordered with myid2
         !call ccmpi_node_ioreaders
#else
         if ( myid == 0 ) then  
            write(6,*) "Set procmode=1 as CCAM was compiled without -Dusempi3"
         end if  
         procmode = 1
         ! Intra-procmode communicator
         comm_vnode = comm_node
         vnode_nproc = node_nproc ! =1
         vnode_myid = node_myid   ! =0
         ! Inter-procmode communicator
         comm_vleader = comm_world
         vleader_nproc = nproc
         vleader_myid = myid
         vnode_vleaderid = vleader_myid
#endif
      else
         ! localhist = .false.
         procmode = nproc
         comm_vnode = comm_world
         vnode_nproc = nproc
         vnode_myid = myid
         comm_vleader = comm_world
         vleader_nproc = nproc
         vleader_myid = myid
         vnode_vleaderid = 0
      end if

   end subroutine ccmpi_procformat_init
   
   ! This routine allows multi-grid bounds updates
   ! The code is based on cc_mpi bounds routines, but
   ! accomodates the g-th multi-grid
   ! This subroutine merges datasets when upscaling with the multi-grid solver
   subroutine mgcollectreduce(g,vdat,dsolmax,klim)

      integer, intent(in) :: g
      integer, intent(in), optional :: klim
      integer :: kx, msg_len, ipanx, jpanx
      real, dimension(:,:), intent(inout) :: vdat
      real, dimension(:), intent(inout) :: dsolmax

      ! merge length
      if ( mg(g)%merge_len <= 1 ) return

      if ( present(klim) ) then
         kx = klim
      else
         kx = size(vdat,2)      
      end if

      msg_len = mg(g)%ifull/(mg(g)%merge_len*mg(g)%npanx) ! message unit size
      ipanx = mg(g)%ipan
      jpanx = mg(g)%ifull/(mg(g)%ipan*mg(g)%npanx)
      call mgcollectreduce_work( g, vdat, dsolmax, kx, mg(g)%nmax, msg_len, mg(g)%npanx, ipanx, jpanx )
  
   return
   end subroutine mgcollectreduce

   subroutine mgcollectreduce_work(g,vdat,dsolmax,kx,nmax,msg_len,npanx,ipanx,jpanx)

      integer, intent(in) :: g, kx, nmax, msg_len, npanx, ipanx, jpanx
      integer :: nrow, ncol
      integer :: yproc, ir, ic, is, js, k, n, j, iq, iqq
      integer :: nrm1, hoz_len
      integer(kind=4) :: ierr, ilen, lcomm
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      real, dimension(:,:), intent(inout) :: vdat
      real, dimension(:), intent(inout) :: dsolmax
      real, dimension((msg_len*npanx+1)*kx) :: tdat
      real, dimension((msg_len*npanx+1)*kx,nmax) :: tdat_g

      ! prep data for sending around the merge
      nrow    = mg(g)%ipan/mg(g)%merge_row  ! number of points along a row per processor
      ncol    = msg_len/nrow                ! number of points along a col per processor
      nrm1    = nrow - 1
      hoz_len = msg_len*npanx
      ilen = (hoz_len+1)*kx

      ! pack contiguous buffer
      tdat(1:hoz_len*kx) = reshape( vdat(1:hoz_len,1:kx), (/ hoz_len*kx /) )
      tdat(hoz_len*kx+1:ilen) = dsolmax(1:kx)

      lcomm = mg(g)%comm_merge
      call START_LOG(gathermg_begin)
      call MPI_Gather( tdat, ilen, ltype, tdat_g, ilen, ltype, 0_4, lcomm, ierr )      
      call END_LOG(gathermg_end)

      ! unpack buffers (nmax is zero unless this is the host processor)
      if ( nmax > 0 ) then
         ! general case
         do yproc = 1,nmax
            ir = mod(yproc-1,mg(g)%merge_row) + 1   ! index for proc row
            ic = (yproc-1)/mg(g)%merge_row + 1      ! index for proc col
            is = (ir-1)*nrow + 1
            js = (ic-1)*ncol + 1
            do k = 1,kx
               do n = 1,npanx
                  do j = js,js+ncol-1
                     iq = (j-1)*ipanx + (n-1)*ipanx*jpanx
                     iqq = 1 - is + (j-js)*nrow + (n-1)*nrow*ncol + (k-1)*nrow*ncol*npanx
                     vdat(iq+is:iq+is+nrow-1,k) = tdat_g(iqq+is:iqq+is+nrow-1,yproc)
                  end do
               end do
            end do   
         end do
         dsolmax(1:kx) = maxval( tdat_g(hoz_len*kx+1:ilen,1:nmax), dim=2 )
      end if
  
   return
   end subroutine mgcollectreduce_work

   ! This routing collects data from other processors without a reduction (max or min) array
   subroutine mgcollect1(g,vdat,klim)

      integer, intent(in) :: g
      integer, intent(in), optional :: klim
      integer :: kx, msg_len, ipanx, jpanx
      real, dimension(:,:), intent(inout) :: vdat

      ! merge length
      if ( mg(g)%merge_len <= 1 ) return
      
      if (present(klim)) then
         kx = klim
      else
         kx = size(vdat,2)
      end if

      msg_len = mg(g)%ifull/(mg(g)%merge_len*mg(g)%npanx) ! message unit size
      ipanx = mg(g)%ipan
      jpanx = mg(g)%ifull/(mg(g)%ipan*mg(g)%npanx)
      call mgcollect_work( g, vdat, kx, mg(g)%nmax, msg_len, mg(g)%npanx, ipanx, jpanx )
  
   return
   end subroutine mgcollect1

   subroutine mgcollect_work(g,vdat,kx,nmax,msg_len,npanx,ipanx,jpanx)

      integer, intent(in) :: g, kx, nmax, msg_len, npanx, ipanx, jpanx
      integer :: nrow, ncol
      integer :: yproc, ir, ic, is, js, k, n, j, iq, iqq
      integer :: nrm1, hoz_len
      integer(kind=4) :: ierr, ilen, lcomm
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      real, dimension(:,:), intent(inout) :: vdat
      real, dimension(msg_len*npanx*kx) :: tdat
      real, dimension(msg_len*npanx*kx,nmax) :: tdat_g

      ! prep data for sending around the merge
      nrow    = mg(g)%ipan/mg(g)%merge_row       ! number of points along a row per processor
      ncol    = msg_len/nrow                     ! number of points along a col per processor
      nrm1    = nrow - 1
      hoz_len = msg_len*npanx
      ilen = hoz_len*kx

      ! pack contiguous buffer
      tdat(1:ilen) = reshape( vdat(1:hoz_len,1:kx), (/ hoz_len*kx /) )

      lcomm = mg(g)%comm_merge
      call START_LOG(gathermg_begin)
      call MPI_Gather( tdat, ilen, ltype, tdat_g, ilen, ltype, 0_4, lcomm, ierr )
      call END_LOG(gathermg_end)

      ! unpack buffers (nmax is zero unless this is the host processor)
      if ( nmax > 0 ) then
         do yproc = 1,nmax
            ir = mod(yproc-1,mg(g)%merge_row) + 1   ! index for proc row
            ic = (yproc-1)/mg(g)%merge_row + 1      ! index for proc col
            is = (ir-1)*nrow + 1
            js = (ic-1)*ncol + 1
            do k = 1,kx
               do n = 1,npanx
                  do j = js,js+ncol-1
                     iq = (j-1)*ipanx + (n-1)*ipanx*jpanx
                     iqq = 1 - is + (j-js)*nrow + (n-1)*nrow*ncol + (k-1)*nrow*ncol*npanx
                     vdat(iq+is:iq+is+nrow-1,k) = tdat_g(iqq+is:iqq+is+nrow-1,yproc)
                  end do
               end do
            end do   
         end do
      end if

   return
   end subroutine mgcollect_work

   ! This version of mgcollect also performs a max and min reduction
   subroutine mgcollectxn(g,vdat,smaxmin,klim)

      integer, intent(in) :: g
      integer, intent(in), optional :: klim
      integer kx, msg_len, ipanx, jpanx
      real, dimension(:,:), intent(inout) :: vdat
      real, dimension(:,:), intent(inout) :: smaxmin

      ! merge length
      if ( mg(g)%merge_len <= 1 ) return
      
      if (present(klim)) then
         kx = klim
      else
         kx = size(vdat,2)
      end if

      msg_len = mg(g)%ifull/(mg(g)%merge_len*mg(g)%npanx) ! message unit size
      ipanx = mg(g)%ipan
      jpanx = mg(g)%ifull/(mg(g)%ipan*mg(g)%npanx)
      call mgcollectxn_work( g, vdat, smaxmin, kx, mg(g)%nmax, msg_len, mg(g)%npanx, ipanx, jpanx )
  
   return
   end subroutine mgcollectxn

   subroutine mgcollectxn_work(g,vdat,smaxmin,kx,nmax,msg_len,npanx,ipanx,jpanx)

      integer, intent(in) :: g, kx, nmax, msg_len, npanx, ipanx, jpanx
      integer :: nrow, ncol
      integer :: yproc, ir, ic, is, js, k, n, j, iq, iqq
      integer :: nrm1, hoz_len
      integer(kind=4) :: ierr, ilen, lcomm
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      real, dimension(:,:), intent(inout) :: vdat
      real, dimension(:,:), intent(inout) :: smaxmin
      real, dimension((msg_len*npanx+2)*kx) :: tdat
      real, dimension((msg_len*npanx+2)*kx,nmax) :: tdat_g

      ! prep data for sending around the merge
      nrow    = mg(g)%ipan/mg(g)%merge_row  ! number of points along a row per processor
      ncol    = msg_len/nrow                ! number of points along a col per processor
      nrm1    = nrow - 1
      hoz_len = msg_len*npanx
      ilen = (hoz_len+2)*kx

      ! pack contiguous buffer
      tdat(1:hoz_len*kx) = reshape( vdat(1:hoz_len,1:kx), (/ hoz_len*kx /) )
      tdat(hoz_len*kx+1:ilen) = reshape( smaxmin(1:kx,1:2), (/ kx*2 /) )

      lcomm = mg(g)%comm_merge
      call START_LOG(gathermg_begin)
      call MPI_Gather( tdat, ilen, ltype, tdat_g, ilen, ltype, 0_4, lcomm, ierr )
      call END_LOG(gathermg_end)

      ! unpack buffers (nmax is zero unless this is the host processor)
      if ( nmax > 0 ) then
         do yproc = 1,nmax
            ir = mod(yproc-1,mg(g)%merge_row) + 1   ! index for proc row
            ic = (yproc-1)/mg(g)%merge_row + 1      ! index for proc col
            is = (ir-1)*nrow + 1
            js = (ic-1)*ncol + 1
            do k = 1,kx
               do n = 1,npanx
                  do j = js,js+ncol-1
                     iq = (j-1)*ipanx + (n-1)*ipanx*jpanx
                     iqq = 1 - is + (j-js)*nrow + (n-1)*nrow*ncol + (k-1)*nrow*ncol*npanx
                     vdat(iq+is:iq+is+nrow-1,k) = tdat_g(iqq+is:iqq+is+nrow-1,yproc)
                  end do
               end do
            end do   
         end do
         smaxmin(1:kx,1) = maxval( tdat_g(hoz_len*kx+1:(hoz_len+1)*kx,1:nmax), dim=2 )
         smaxmin(1:kx,2) = minval( tdat_g((hoz_len+1)*kx+1:ilen,1:nmax), dim=2 )
      end if   
  
   return
   end subroutine mgcollectxn_work

   ! This subroutine merges datasets when upscaling with the multi-grid solver
   ! This version also updates the halo

   subroutine mgbcast2(g,vdat,dsolmax,nobounds)
   
      integer, intent(in) :: g
      integer :: out_len
      integer(kind=4) :: ierr, ilen, lcomm
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      real, dimension(:), intent(inout) :: vdat
      real, dimension(:), intent(inout) :: dsolmax
      real, dimension(mg(g)%ifull+mg(g)%iextra+size(dsolmax)) :: tdat
      logical, intent(in), optional :: nobounds
      logical :: nbflag
      
      if ( mg(g)%merge_len <= 1 ) return
      
      nbflag = .false.
      if ( present(nobounds) ) then
         nbflag = nobounds
      end if   
         
      if ( nbflag ) then
         out_len = mg(g)%ifull
      else   
         out_len = mg(g)%ifull + mg(g)%iextra
      end if  
      ilen = out_len + size(dsolmax)
      
      ! pack contiguous buffer
      tdat(1:out_len) = vdat(1:out_len)
      tdat(out_len+1:ilen) = dsolmax(:)
      
      lcomm = mg(g)%comm_merge
      call START_LOG(bcastmg_begin)
      call MPI_Bcast( tdat, ilen, ltype, 0_4, lcomm, ierr )
      call END_LOG(bcastmg_end)

      ! extract data from Bcast
      vdat(1:out_len) = tdat(1:out_len)
      dsolmax(:) = tdat(out_len+1:ilen)
   
   return
   end subroutine mgbcast2
   
   subroutine mgbcast3(g,vdat,dsolmax,klim,nobounds)
   
      integer, intent(in) :: g
      integer, intent(in), optional :: klim
      integer :: kx, out_len
      integer(kind=4) :: ierr, ilen, lcomm
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      real, dimension(:,:), intent(inout) :: vdat
      real, dimension(:), intent(inout) :: dsolmax
      real, dimension((mg(g)%ifull+mg(g)%iextra+1)*size(vdat,2)) :: tdat
      logical, intent(in), optional :: nobounds
      logical :: nbflag
      
      if ( mg(g)%merge_len <= 1 ) return

      if ( present(klim) ) then
         kx = klim
      else
         kx = size(vdat,2)
      end if
      
      nbflag = .false.
      if ( present(nobounds) ) then
         nbflag = nobounds
      end if   
         
      if ( nbflag ) then
         out_len = mg(g)%ifull
      else   
         out_len = mg(g)%ifull + mg(g)%iextra
      end if  
      ilen = (out_len+1)*kx
      
      ! pack contiguous buffer
      tdat(1:out_len*kx) = reshape( vdat(1:out_len,1:kx), (/ out_len*kx /) )
      tdat(out_len*kx+1:ilen) = dsolmax(1:kx)

      lcomm = mg(g)%comm_merge
      call START_LOG(bcastmg_begin)
      call MPI_Bcast( tdat, ilen, ltype, 0_4, lcomm, ierr )
      call END_LOG(bcastmg_end)

      ! extract data from Bcast
      vdat(1:out_len,1:kx) = reshape( tdat(1:out_len*kx), (/ out_len, kx /) )
      dsolmax(1:kx) = tdat(out_len*kx+1:ilen)
   
   return
   end subroutine mgbcast3

   subroutine mgbcasta2(g,vdat,nobounds)
   
      integer, intent(in) :: g
      integer :: out_len
      integer(kind=4) :: ierr, ilen, lcomm
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      real, dimension(:), intent(inout) :: vdat
      logical, intent(in), optional :: nobounds
      logical :: nbflag
      
      if ( mg(g)%merge_len <= 1 ) return
      
      nbflag = .false.
      if ( present(nobounds) ) then
         nbflag = nobounds
      end if   
         
      if ( nbflag ) then
         out_len = mg(g)%ifull
      else   
         out_len = mg(g)%ifull + mg(g)%iextra
      end if  

      ilen = out_len
      lcomm = mg(g)%comm_merge
      call START_LOG(bcastmg_begin)
      call MPI_Bcast( vdat, ilen, ltype, 0_4, lcomm, ierr )      
      call END_LOG(bcastmg_end)
   
   return
   end subroutine mgbcasta2
   
   subroutine mgbcasta3(g,vdat,nobounds)
   
      integer, intent(in) :: g
      integer :: kx, out_len
      integer(kind=4) :: ierr, ilen, lcomm
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      real, dimension(:,:), intent(inout) :: vdat
      real, dimension((mg(g)%ifull + mg(g)%iextra)*size(vdat,2)) :: tdat
      logical, intent(in), optional :: nobounds
      logical :: nbflag

      if ( mg(g)%merge_len <= 1 ) return

      kx = size(vdat,2)
      
      nbflag = .false.
      if ( present(nobounds) ) then
         nbflag = nobounds
      end if   
         
      if ( nbflag ) then
         out_len = mg(g)%ifull
      else   
         out_len = mg(g)%ifull + mg(g)%iextra
      end if  
      ilen = out_len*kx

      ! pack contiguous buffer
      tdat(1:ilen) = reshape( vdat(1:out_len,1:kx), (/ out_len*kx /) )
      
      lcomm = mg(g)%comm_merge
      call START_LOG(bcastmg_begin)
      call MPI_Bcast( tdat, ilen, ltype, 0_4, lcomm, ierr )      
      call END_LOG(bcastmg_end)

      ! extract data from Bcast
      vdat(1:out_len,1:kx) = reshape( tdat(1:ilen), (/ out_len, kx /) )
   
   return
   end subroutine mgbcasta3
   
   subroutine mgbcastxn(g,vdat,smaxmin,klim,nobounds)
   
      integer, intent(in) :: g
      integer, intent(in), optional :: klim
      integer :: kx, out_len
      integer(kind=4) :: ierr, ilen, lcomm
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      real, dimension(:,:), intent(inout) :: vdat
      real, dimension(:,:), intent(inout) :: smaxmin
      real, dimension((mg(g)%ifull+mg(g)%iextra+2)*size(vdat,2)) :: tdat
      logical, intent(in), optional :: nobounds
      logical :: nbflag
      
      if ( mg(g)%merge_len <= 1 ) return

      if (present(klim)) then
         kx = klim
      else
         kx = size(vdat,2)
      end if
      
      nbflag = .false.
      if ( present(nobounds) ) then
         nbflag = nobounds
      end if   
         
      if ( nbflag ) then
         out_len = mg(g)%ifull
      else   
         out_len = mg(g)%ifull + mg(g)%iextra
      end if   
      ilen = (out_len+2)*kx

      ! pack contiguous buffer
      tdat(1:out_len*kx) = reshape( vdat(1:out_len,1:kx), (/ out_len*kx /) )
      tdat(out_len*kx+1:ilen) = reshape( smaxmin(1:kx,1:2), (/ kx*2 /) )
      
      lcomm = mg(g)%comm_merge
      call START_LOG(bcastmg_begin)
      call MPI_Bcast( tdat, ilen, ltype, 0_4, lcomm, ierr )      
      call END_LOG(bcastmg_end)

      ! extract data from Bcast
      vdat(1:out_len,1:kx) = reshape( tdat(1:out_len*kx), (/ out_len, kx /) )
      smaxmin(1:kx,1:2) = reshape( tdat(out_len*kx+1:ilen), (/ kx, 2 /) )
   
   return
   end subroutine mgbcastxn
   
   ! Set up the indices required for the multigrid scheme.
   subroutine mg_index(g,mil_g,mipan,mjpan)

      use indices_m
   
      integer, intent(in) :: g, mil_g, mipan, mjpan
      integer, dimension(:), allocatable :: mg_colourmask
      integer, dimension(2*(mipan+mjpan+2)*(npanels+1)) :: dum
      integer, dimension(3) :: mg_ifull_colour
      integer, dimension(:,:), allocatable :: dums, dumr
      integer :: mioff, mjoff
      integer :: i, j, n, iq, iqq, iqg, mfull_g
      integer :: iext, iproc, xlen, jx, nc, xlev, rproc, sproc
      integer :: ntest, ncount
      integer(kind=4) :: itag=22, lproc, ierr, llen, lcomm
      integer(kind=4), dimension(:,:), allocatable :: status
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_INTEGER8
#else
      integer(kind=4), parameter :: ltype = MPI_INTEGER
#endif
      logical lflag, lglob      

      ! size of this grid
      mfull_g = 6*mil_g*mil_g


      ! calculate processor map in iq coordinates
      lglob = .true.
      lflag = .true.
      mioff = 0
      mjoff = 0
      do n = 0,npanels
         do j = 1,mil_g
            do i = 1,mil_g
               iq = indx(i,j,n,mil_g,mil_g)
               if ( mg_qproc(iq,mil_g,g) /= myid ) then
                  ! found grid point that does not belong to myid
                  lglob = .false.
               else if ( lflag ) then
                  ! first grid point for myid
                  mioff = i - 1
                  mjoff = j - 1
                  lflag = .false.
               end if
            end do
         end do
      end do
      
      if ( lflag ) then
         write(6,*) "ERROR: Cannot find myid in mg_proc"
         write(6,*) "myid,g ",myid,g
         call ccmpi_abort(-1)
      end if
      

      mg_bnds(:,g)%len = 0
      mg_bnds(:,g)%rlen = 0
      mg_bnds(:,g)%slen = 0
      mg_bnds(:,g)%rlenx = 0
      mg_bnds(:,g)%slenx = 0

      ! Calculate local indices on this processor
      if ( lglob ) then
         do iq = 1,mfull_g
            mg(g)%in(iq) = jn_g(iq,mil_g)
            mg(g)%is(iq) = js_g(iq,mil_g)
            mg(g)%ie(iq) = je_g(iq,mil_g)
            mg(g)%iw(iq) = jw_g(iq,mil_g)
            mg(g)%ine(iq) = jne_g(iq,mil_g)
            mg(g)%ien(iq) = jen_g(iq,mil_g)
            mg(g)%ise(iq) = jse_g(iq,mil_g)
            mg(g)%ies(iq) = jes_g(iq,mil_g)
            mg(g)%inw(iq) = jnw_g(iq,mil_g)
            mg(g)%iwn(iq) = jwn_g(iq,mil_g)
            mg(g)%isw(iq) = jsw_g(iq,mil_g)
            mg(g)%iws(iq) = jws_g(iq,mil_g)
         end do
         mg(g)%iextra = 0
         mg(g)%neighnum = 0
         allocate ( mg(g)%neighlist(mg(g)%neighnum) )
      else
         mg(g)%iextra = 2*(mipan+mjpan+4)*npan ! first guess

         ! This only occurs with grids prior to globgath.  So npan and noff are still valid.
         do n = 1,npan
            do j = 1,mjpan
               do i = 1,mipan
                  iq = indx(i,j,n-1,mipan,mjpan)                 ! Local
                  iqg = indx(i+mioff,j+mjoff,n-noff,mil_g,mil_g) ! Global

                  iqq = jn_g(iqg,mil_g)         ! Global neighbour index
                  rproc = mg_qproc(iqq,mil_g,g) ! Processor that has this point
                  if ( rproc == myid ) then ! Just copy the value
                     ! Convert global iqq to local value
                     mg(g)%in(iq) = indx_indv(iqq,mil_g,mipan,mjpan,mioff,mjoff,noff)
                  end if

                  iqq = js_g(iqg,mil_g)         ! Global neighbour index
                  rproc = mg_qproc(iqq,mil_g,g) ! Processor that has this point
                  if ( rproc == myid ) then ! Just copy the value
                     ! Convert global iqq to local value
                     mg(g)%is(iq) = indx_indv(iqq,mil_g,mipan,mjpan,mioff,mjoff,noff)
                  end if

                  iqq = je_g(iqg,mil_g)         ! Global neighbour index
                  rproc = mg_qproc(iqq,mil_g,g) ! Processor that has this point
                  if ( rproc == myid ) then ! Just copy the value
                     ! Convert global iqq to local value
                     mg(g)%ie(iq) = indx_indv(iqq,mil_g,mipan,mjpan,mioff,mjoff,noff)
                  end if

                  iqq = jw_g(iqg,mil_g)         ! Global neighbour index
                  rproc = mg_qproc(iqq,mil_g,g) ! Processor that has this point
                  if ( rproc == myid ) then ! Just copy the value
                     ! Convert global iqq to local value
                     mg(g)%iw(iq) = indx_indv(iqq,mil_g,mipan,mjpan,mioff,mjoff,noff)
                  end if
            
                  iqq = jne_g(iqg,mil_g)        ! Global neighbour index
                  rproc = mg_qproc(iqq,mil_g,g) ! Processor that has this point
                  if ( rproc == myid ) then ! Just copy the value
                     ! Convert global iqq to local value
                     mg(g)%ine(iq) = indx_indv(iqq,mil_g,mipan,mjpan,mioff,mjoff,noff)
                  end if      

                  iqq = jen_g(iqg,mil_g)        ! Global neighbour index
                  rproc = mg_qproc(iqq,mil_g,g) ! Processor that has this point
                  if ( rproc == myid ) then ! Just copy the value
                     ! Convert global iqq to local value
                     mg(g)%ien(iq) = indx_indv(iqq,mil_g,mipan,mjpan,mioff,mjoff,noff)
                  end if

                  iqq = jnw_g(iqg,mil_g)        ! Global neighbour index
                  rproc = mg_qproc(iqq,mil_g,g) ! Processor that has this point
                  if ( rproc == myid ) then ! Just copy the value
                     ! Convert global iqq to local value
                     mg(g)%inw(iq) = indx_indv(iqq,mil_g,mipan,mjpan,mioff,mjoff,noff)
                  end if
 
                  iqq = jwn_g(iqg,mil_g)        ! Global neighbour index
                  rproc = mg_qproc(iqq,mil_g,g) ! Processor that has this point
                  if ( rproc == myid ) then ! Just copy the value
                     ! Convert global iqq to local value
                     mg(g)%iwn(iq) = indx_indv(iqq,mil_g,mipan,mjpan,mioff,mjoff,noff)
                  end if

                  iqq = jse_g(iqg,mil_g)        ! Global neighbour index
                  rproc = mg_qproc(iqq,mil_g,g) ! Processor that has this point
                  if ( rproc == myid ) then ! Just copy the value
                     ! Convert global iqq to local value
                     mg(g)%ise(iq) = indx_indv(iqq,mil_g,mipan,mjpan,mioff,mjoff,noff)
                  end if

                  iqq = jes_g(iqg,mil_g)        ! Global neighbour index
                  rproc = mg_qproc(iqq,mil_g,g) ! Processor that has this point
                  if ( rproc == myid ) then ! Just copy the value
                     ! Convert global iqq to local value
                     mg(g)%ies(iq) = indx_indv(iqq,mil_g,mipan,mjpan,mioff,mjoff,noff)
                  end if

                  iqq = jsw_g(iqg,mil_g)        ! Global neighbour index
                  rproc = mg_qproc(iqq,mil_g,g) ! Processor that has this point
                  if ( rproc == myid ) then ! Just copy the value
                     ! Convert global iqq to local value
                     mg(g)%isw(iq) = indx_indv(iqq,mil_g,mipan,mjpan,mioff,mjoff,noff)
                  end if
 
                  iqq = jws_g(iqg,mil_g)        ! Global neighbour index
                  rproc = mg_qproc(iqq,mil_g,g) ! Processor that has this point
                  if ( rproc == myid ) then ! Just copy the value
                     ! Convert global iqq to local value
                     mg(g)%iws(iq) = indx_indv(iqq,mil_g,mipan,mjpan,mioff,mjoff,noff)
                  end if
                  
               end do
            end do
         end do


         ! Calculate local indices in halo
         iext = 0
         do n = 1,npan

            !     Start with N edge
            j = mjpan
            do i = 1,mipan
               iq = indx(i,j,n-1,mipan,mjpan)  !  Local index 
               iqg = indx(i+mioff,j+mjoff,n-noff,mil_g,mil_g)
               iqq = jn_g(iqg,mil_g)
               ! Which processor has this point
               rproc = mg_qproc(iqq,mil_g,g)
               if ( rproc /= myid ) then ! Don't add points already on this proc.
                  ! Add this point to request list
                  call mgcheck_bnds_alloc(g, rproc, iext)
                  mg_bnds(rproc,g)%rlen = mg_bnds(rproc,g)%rlen + 1
                  mg_bnds(rproc,g)%request_list(mg_bnds(rproc,g)%rlen) = iqq
                  ! Increment extended region index
                  iext = iext + 1
                  mg_bnds(rproc,g)%unpack_list(mg_bnds(rproc,g)%rlen) = iext
                  mg(g)%in(iq) = mg(g)%ifull + iext
               end if   
            end do

            !     E edge
            i = mipan
            do j = 1,mjpan
               iq = indx(i,j,n-1,mipan,mjpan)  !  Local index 
               iqg = indx(i+mioff,j+mjoff,n-noff,mil_g,mil_g)
               iqq = je_g(iqg,mil_g)
               ! Which processor has this point
               rproc = mg_qproc(iqq,mil_g,g)
               if ( rproc /= myid ) then ! Don't add points already on this proc.
                  ! Add this point to request list
                  call mgcheck_bnds_alloc(g, rproc, iext)
                  mg_bnds(rproc,g)%rlen = mg_bnds(rproc,g)%rlen + 1
                  mg_bnds(rproc,g)%request_list(mg_bnds(rproc,g)%rlen) = iqq
                  ! Increment extended region index
                  iext = iext + 1
                  mg_bnds(rproc,g)%unpack_list(mg_bnds(rproc,g)%rlen) = iext
                  mg(g)%ie(iq) = mg(g)%ifull + iext
               end if   
            end do

            !     W edge
            i = 1
            do j = 1,mjpan
               iq = indx(i,j,n-1,mipan,mjpan)  !  Local index  
               iqg = indx(i+mioff,j+mjoff,n-noff,mil_g,mil_g)
               iqq = jw_g(iqg,mil_g)
               ! Which processor has this point
               rproc = mg_qproc(iqq,mil_g,g)
               if ( rproc /= myid ) then ! Don't add points already on this proc.
                  ! Add this point to request list
                  call mgcheck_bnds_alloc(g, rproc, iext)
                  mg_bnds(rproc,g)%rlen = mg_bnds(rproc,g)%rlen + 1
                  mg_bnds(rproc,g)%request_list(mg_bnds(rproc,g)%rlen) = iqq
                  ! Increment extended region index
                  iext = iext + 1
                  mg_bnds(rproc,g)%unpack_list(mg_bnds(rproc,g)%rlen) = iext
                  mg(g)%iw(iq) = mg(g)%ifull + iext
               end if   
            end do

            !     S edge
            j = 1
            do i = 1,mipan
               iq = indx(i,j,n-1,mipan,mjpan)  !  Local index
               iqg = indx(i+mioff,j+mjoff,n-noff,mil_g,mil_g)
               iqq = js_g(iqg,mil_g)
               ! Which processor has this point
               rproc = mg_qproc(iqq,mil_g,g)
               if ( rproc /= myid ) then ! Don't add points already on this proc.
                  ! Add this point to request list
                  call mgcheck_bnds_alloc(g, rproc, iext)
                  mg_bnds(rproc,g)%rlen = mg_bnds(rproc,g)%rlen + 1
                  mg_bnds(rproc,g)%request_list(mg_bnds(rproc,g)%rlen) = iqq
                  ! Increment extended region index
                  iext = iext + 1
                  mg_bnds(rproc,g)%unpack_list(mg_bnds(rproc,g)%rlen) = iext
                  mg(g)%is(iq) = mg(g)%ifull + iext
                end if  
            end do
            
         end do ! n=1,npan
         
         mg_bnds(:,g)%rlenx = mg_bnds(:,g)%rlen
         
         do n = 1,npan
            
            ! NE, EN
            iq = indx(mipan,mjpan,n-1,mipan,mjpan)  !  Local index
            iqg = indx(mipan+mioff,mjpan+mjoff,n-noff,mil_g,mil_g)
            iqq = jne_g(iqg,mil_g)
            ! Which processor has this point
            rproc = mg_qproc(iqq,mil_g,g)
            if ( rproc /= myid ) then ! Don't add points already on this proc.
               ! Add this point to request list
               call mgcheck_bnds_alloc(g, rproc, iext)
               mg_bnds(rproc,g)%rlenx = mg_bnds(rproc,g)%rlenx + 1
               mg_bnds(rproc,g)%request_list(mg_bnds(rproc,g)%rlenx) = iqq
               ! Increment extended region index
               iext = iext + 1
               mg_bnds(rproc,g)%unpack_list(mg_bnds(rproc,g)%rlenx) = iext
               mg(g)%ine(iq) = mg(g)%ifull + iext
            end if  
            if ( jen_g(iqg,mil_g) == jne_g(iqg,mil_g) ) then
               mg(g)%ien(iq) = mg(g)%ine(iq)
            else
               iqq = jen_g(iqg,mil_g)
               ! Which processor has this point
               rproc = mg_qproc(iqq,mil_g,g)
               if ( rproc /= myid ) then ! Add to list
                 ! Add this point to request list
                 call mgcheck_bnds_alloc(g, rproc, iext)
                 mg_bnds(rproc,g)%rlenx = mg_bnds(rproc,g)%rlenx + 1
                 mg_bnds(rproc,g)%request_list(mg_bnds(rproc,g)%rlenx) = iqq
                 ! Increment extended region index
                 iext = iext + 1
                 mg_bnds(rproc,g)%unpack_list(mg_bnds(rproc,g)%rlenx) = iext
                 mg(g)%ien(iq) = mg(g)%ifull + iext
               end if
            end if

            ! SE, ES
            iq = indx(mipan,1,n-1,mipan,mjpan)  !  Local index
            iqg = indx(mipan+mioff,1+mjoff,n-noff,mil_g,mil_g)
            iqq = jse_g(iqg,mil_g)
            ! Which processor has this point
            rproc = mg_qproc(iqq,mil_g,g)
            if ( rproc /= myid ) then ! Don't add points already on this proc.
               ! Add this point to request list
               call mgcheck_bnds_alloc(g, rproc, iext)
               mg_bnds(rproc,g)%rlenx = mg_bnds(rproc,g)%rlenx + 1
               mg_bnds(rproc,g)%request_list(mg_bnds(rproc,g)%rlenx) = iqq
               ! Increment extended region index
               iext = iext + 1
               mg_bnds(rproc,g)%unpack_list(mg_bnds(rproc,g)%rlenx) = iext
               mg(g)%ise(iq) = mg(g)%ifull + iext
            end if  
            if ( jes_g(iqg,mil_g) == jse_g(iqg,mil_g) ) then
               mg(g)%ies(iq) = mg(g)%ise(iq)
            else
               iqq = jes_g(iqg,mil_g)
               ! Which processor has this point
               rproc = mg_qproc(iqq,mil_g,g)
               if ( rproc /= myid ) then ! Add to list
                 ! Add this point to request list
                 call mgcheck_bnds_alloc(g, rproc, iext)
                 mg_bnds(rproc,g)%rlenx = mg_bnds(rproc,g)%rlenx + 1
                 mg_bnds(rproc,g)%request_list(mg_bnds(rproc,g)%rlenx) = iqq
                 ! Increment extended region index
                 iext = iext + 1
                 mg_bnds(rproc,g)%unpack_list(mg_bnds(rproc,g)%rlenx) = iext
                 mg(g)%ies(iq) = mg(g)%ifull + iext
               end if
            end if

            ! NW, WN
            iq = indx(1,mjpan,n-1,mipan,mjpan)  !  Local index
            iqg = indx(1+mioff,mjpan+mjoff,n-noff,mil_g,mil_g)
            iqq = jnw_g(iqg,mil_g)
            ! Which processor has this point
            rproc = mg_qproc(iqq,mil_g,g)
            if ( rproc /= myid ) then ! Don't add points already on this proc.
               ! Add this point to request list
               call mgcheck_bnds_alloc(g, rproc, iext)
               mg_bnds(rproc,g)%rlenx = mg_bnds(rproc,g)%rlenx + 1
               mg_bnds(rproc,g)%request_list(mg_bnds(rproc,g)%rlenx) = iqq
               ! Increment extended region index
               iext = iext + 1
               mg_bnds(rproc,g)%unpack_list(mg_bnds(rproc,g)%rlenx) = iext
               mg(g)%inw(iq) = mg(g)%ifull + iext
            end if  
            if ( jwn_g(iqg,mil_g) == jnw_g(iqg,mil_g) ) then
               mg(g)%iwn(iq) = mg(g)%inw(iq)
            else
               iqq = jwn_g(iqg,mil_g)
               ! Which processor has this point
               rproc = mg_qproc(iqq,mil_g,g)
               if ( rproc /= myid ) then ! Add to list
                 ! Add this point to request list
                 call mgcheck_bnds_alloc(g, rproc, iext)
                 mg_bnds(rproc,g)%rlenx = mg_bnds(rproc,g)%rlenx + 1
                 mg_bnds(rproc,g)%request_list(mg_bnds(rproc,g)%rlenx) = iqq
                 ! Increment extended region index
                 iext = iext + 1
                 mg_bnds(rproc,g)%unpack_list(mg_bnds(rproc,g)%rlenx) = iext
                 mg(g)%iwn(iq) = mg(g)%ifull + iext
               end if
            end if

            ! SW, WS
            iq = indx(1,1,n-1,mipan,mjpan)  !  Local index
            iqg = indx(1+mioff,1+mjoff,n-noff,mil_g,mil_g)
            iqq = jsw_g(iqg,mil_g)
            ! Which processor has this point
            rproc = mg_qproc(iqq,mil_g,g)
            if ( rproc /= myid ) then ! Don't add points already on this proc.
               ! Add this point to request list
               call mgcheck_bnds_alloc(g, rproc, iext)
               mg_bnds(rproc,g)%rlenx = mg_bnds(rproc,g)%rlenx + 1
               mg_bnds(rproc,g)%request_list(mg_bnds(rproc,g)%rlenx) = iqq
               ! Increment extended region index
               iext = iext + 1
               mg_bnds(rproc,g)%unpack_list(mg_bnds(rproc,g)%rlenx) = iext
               mg(g)%isw(iq) = mg(g)%ifull + iext
            end if  
            if ( jws_g(iqg,mil_g) == jsw_g(iqg,mil_g) ) then
               mg(g)%iws(iq) = mg(g)%isw(iq)
            else
               iqq = jws_g(iqg,mil_g)
               ! Which processor has this point
               rproc = mg_qproc(iqq,mil_g,g)
               if ( rproc /= myid ) then ! Add to list
                 ! Add this point to request list
                 call mgcheck_bnds_alloc(g, rproc, iext)
                 mg_bnds(rproc,g)%rlenx = mg_bnds(rproc,g)%rlenx + 1
                 mg_bnds(rproc,g)%request_list(mg_bnds(rproc,g)%rlenx) = iqq
                 ! Increment extended region index
                 iext = iext + 1
                 mg_bnds(rproc,g)%unpack_list(mg_bnds(rproc,g)%rlenx) = iext
                 mg(g)%iws(iq) = mg(g)%ifull + iext
               end if
            end if
            
         end do ! n=1,npan

         mg(g)%iextra = iext

        
         ! Set up the diagonal index arrays
         do n = 1,npan
            do j = 1,mjpan
               do i = 1,mipan
                  iq = indx(i,j,n-1,mipan,mjpan)   ! Local
                  ! Except at corners, ien = ine etc.
                  if ( i > 1 ) then
                     mg(g)%inw(iq) = mg(g)%in(mg(g)%iw(iq))
                     mg(g)%isw(iq) = mg(g)%is(mg(g)%iw(iq))
                  else
                     if ( j < mjpan ) mg(g)%inw(iq) = mg(g)%iw(mg(g)%in(iq))
                     if ( j > 1 )     mg(g)%isw(iq) = mg(g)%iw(mg(g)%is(iq))
                  end if
                  if ( i < mipan ) then
                     ! ie will be defined
                     mg(g)%ine(iq) = mg(g)%in(mg(g)%ie(iq))
                     mg(g)%ise(iq) = mg(g)%is(mg(g)%ie(iq))
                  else
                     ! i = ipan, ie will have been remapped
                     if ( j > 1 )     mg(g)%ise(iq) = mg(g)%ie(mg(g)%is(iq))
                     if ( j < mjpan ) mg(g)%ine(iq) = mg(g)%ie(mg(g)%in(iq))
                  end if
                  if ( j > 1 ) then
                     mg(g)%ies(iq) = mg(g)%ie(mg(g)%is(iq))
                     mg(g)%iws(iq) = mg(g)%iw(mg(g)%is(iq))
                  else
                     if ( i < mipan ) mg(g)%ies(iq)=mg(g)%is(mg(g)%ie(iq))
                     if ( i > 1 )     mg(g)%iws(iq)=mg(g)%is(mg(g)%iw(iq))
                  end if
                  if ( j < mjpan ) then
                     mg(g)%ien(iq) = mg(g)%ie(mg(g)%in(iq))
                     mg(g)%iwn(iq) = mg(g)%iw(mg(g)%in(iq))
                  else
                     if ( i < mipan) mg(g)%ien(iq) = mg(g)%in(mg(g)%ie(iq))
                     if ( i > 1 )    mg(g)%iwn(iq) = mg(g)%in(mg(g)%iw(iq))
                  end if
               end do
            end do
         end do
      
         
         ! check neighbours
         mg(g)%neighnum = count( mg_bnds(:,g)%rlenx > 0 )

         ! increase size of request list if needed
         ntest = mg(g)%neighnum
         if ( 2*ntest > size(ireq) ) then
            deallocate( ireq, rlist )
            allocate( ireq(max(2*ntest,1)) )
            allocate( rlist(max(ntest,1)) )
         end if
         allocate( status(MPI_STATUS_SIZE,size(ireq)) )
  
         ! set up neighbour lists
         allocate ( mg(g)%neighlist(mg(g)%neighnum) )
         ncount = 0
         do iproc = 1,nproc-1
            rproc = modulo( myid+iproc, nproc )
            if ( mg_bnds(rproc,g)%rlenx > 0 ) then
               ncount = ncount + 1
               mg(g)%neighlist(ncount) = rproc
            end if
         end do
      
         if ( ncount/=mg(g)%neighnum ) then
            write(6,*) "ERROR: Multi-grid neighnum mismatch"
            write(6,*) "neighnum, ncount ",mg(g)%neighnum, ncount
            call ccmpi_abort(-1)
         end if

         
         ! Now, for each processor send the length of points I want.
         allocate( dums(2,mg(g)%neighnum), dumr(2,mg(g)%neighnum) )
         lcomm = comm_world
         nreq = 0
         do iproc = 1,mg(g)%neighnum
            rproc = mg(g)%neighlist(iproc)  ! Recv from
            if ( mg_bnds(rproc,g)%rlenx > 0 ) then
               nreq = nreq + 1
               lproc = rproc
               call MPI_IRecv( dumr(:,iproc), 2_4, ltype, lproc, itag, lcomm, ireq(nreq), ierr )
            end if
         end do
         do iproc = mg(g)%neighnum,1,-1
            sproc = mg(g)%neighlist(iproc)  ! Send to
            if ( mg_bnds(sproc,g)%rlenx > 0 ) then
               nreq = nreq + 1
               dums(1,iproc) = mg_bnds(sproc,g)%rlenx
               dums(2,iproc) = mg_bnds(sproc,g)%rlen
               lproc = sproc
               call MPI_ISend( dums(:,iproc), 2_4, ltype, lproc, itag, lcomm, ireq(nreq), ierr )
            end if
         end do
         call MPI_Waitall( nreq, ireq, status, ierr )
         do iproc = 1,mg(g)%neighnum
            rproc = mg(g)%neighlist(iproc)
            if ( mg_bnds(rproc,g)%rlenx > 0 ) then
               mg_bnds(rproc,g)%slenx = dumr(1,iproc)
               mg_bnds(rproc,g)%slen  = dumr(2,iproc)
            end if
         end do   
         nreq = 0
         rreq = 0
         deallocate( dums, dumr )

  
         ! Now start sending messages  
         lcomm = comm_world
         nreq = 0
         do iproc = 1,mg(g)%neighnum
            lproc = mg(g)%neighlist(iproc)  ! Recv from
            allocate( mg_bnds(lproc,g)%send_list(mg_bnds(lproc,g)%slenx) )
            nreq = nreq + 1
            ! Use the maximum size in the recv call.
            llen = mg_bnds(lproc,g)%slenx
            call MPI_IRecv( mg_bnds(lproc,g)%send_list(1), llen, ltype, lproc, &
                            itag, lcomm, ireq(nreq), ierr )
         end do
         do iproc = mg(g)%neighnum,1,-1
            lproc = mg(g)%neighlist(iproc)  ! Send to
            ! Send list of requests
            nreq = nreq + 1
            llen = mg_bnds(lproc,g)%rlenx
            call MPI_ISend( mg_bnds(lproc,g)%request_list(1), llen, ltype, lproc, &
                            itag, lcomm, ireq(nreq), ierr )
         end do      
         call MPI_Waitall( nreq, ireq, status, ierr )
         nreq = 0
         rreq = 0
        
         deallocate( status )

         ! At the moment send_lists use global indices. Convert these to local.
         do iproc = mg(g)%neighnum,1,-1
            sproc = mg(g)%neighlist(iproc)  ! Send to
            do iq = 1,mg_bnds(sproc,g)%slenx
               ! send_list(iq) is global point index, i, j, n are local
               iqq = mg_bnds(sproc,g)%send_list(iq)
               mg_bnds(sproc,g)%send_list(iq) = indx_indv(iqq,mil_g,mipan,mjpan,mioff,mjoff,noff)
            end do
         end do
         if ( mg_bnds(myid,g)%rlenx /= 0 ) then
            write(6,*) "ERROR: Invalid rlenx in myid"
            call ccmpi_abort(-1)
         end if   

         ! reduce array size where possible
         do iproc = 0,nproc-1
            xlen = mg_bnds(iproc,g)%rlenx
            if ( mg_bnds(iproc,g)%len > xlen ) then
               deallocate( mg_bnds(iproc,g)%request_list )
               dum(1:xlen) = mg_bnds(iproc,g)%unpack_list(1:xlen)
               deallocate( mg_bnds(iproc,g)%unpack_list )
               if ( xlen > 0 ) then
                 allocate( mg_bnds(iproc,g)%unpack_list(xlen) )
                 mg_bnds(iproc,g)%unpack_list(1:xlen) = dum(1:xlen)
               end if
               mg_bnds(iproc,g)%len = xlen
            end if

            ! set up buffers
            xlev = max( kl, 1 ) ! ol is not required
            xlen = xlev*mg_bnds(iproc,g)%rlenx
            if ( bnds(iproc)%rbuflen < xlen ) then
               if ( allocated(bnds(iproc)%rbuf) ) then
                  deallocate( bnds(iproc)%rbuf )
                  deallocate( bnds(iproc)%r8buf )
               end if
               allocate( bnds(iproc)%rbuf(xlen) )
               allocate( bnds(iproc)%r8buf(xlen) )
               bnds(iproc)%rbuflen = xlen
            end if
            xlen = xlev*mg_bnds(iproc,g)%slenx
            if ( bnds(iproc)%sbuflen < xlen ) then
               if ( allocated(bnds(iproc)%sbuf) ) then
                  deallocate( bnds(iproc)%sbuf )
                  deallocate( bnds(iproc)%s8buf )
               end if
               allocate( bnds(iproc)%sbuf(xlen) )
               allocate( bnds(iproc)%s8buf(xlen) )
               bnds(iproc)%sbuflen = xlen
            end if
         end do

      end if

      
      ! calculate colours
      if ( g == mg_maxlevel ) then
          
         if ( myid == 0 ) then
  
            allocate( mg_colourmask(6*mil_g*mil_g) ) 
          
            ! always a three colour mask for coarse grid
            do n = 0,npanels
               do j = 1,mil_g
                  do i = 1,mil_g
                     iq = indx(i,j,n,mil_g,mil_g)

                     jx = mod( i+j+n*mil_g, 2 )
                     select case( n+jx*(npanels+1) )
                        case( 0, 1, 3, 4 )
                           mg_colourmask(iq) = 1
                        case( 2, 5, 6, 9 )
                           mg_colourmask(iq) = 2
                        case( 7, 8, 10, 11 )
                           mg_colourmask(iq) = 3
                     end select
                  end do
               end do
            end do
  
            mg_ifull_maxcolour = count( mg_colourmask == 1 )
            if ( mg_ifull_maxcolour /= count( mg_colourmask == 2 ) .or.  &
                 mg_ifull_maxcolour /= count( mg_colourmask == 3 ) ) then
               write(6,*) "ERROR: Unbalanced MG colours"
               call ccmpi_abort(-1)
            end if
         
            allocate( col_iq(mg_ifull_maxcolour,3), col_iqn(mg_ifull_maxcolour,3), col_iqe(mg_ifull_maxcolour,3) )
            allocate( col_iqs(mg_ifull_maxcolour,3), col_iqw(mg_ifull_maxcolour,3) )
            allocate( col_iqne(mg_ifull_maxcolour,3), col_iqen(mg_ifull_maxcolour,3) )
            allocate( col_iqse(mg_ifull_maxcolour,3), col_iqes(mg_ifull_maxcolour,3) )
            allocate( col_iqnw(mg_ifull_maxcolour,3), col_iqwn(mg_ifull_maxcolour,3) )
            allocate( col_iqsw(mg_ifull_maxcolour,3), col_iqws(mg_ifull_maxcolour,3) )
  
            mg_ifull_colour = 0
            col_iq(:,:) = 0
            col_iqn(:,:) = 0
            col_iqe(:,:) = 0
            col_iqs(:,:) = 0
            col_iqw(:,:) = 0
            col_iqne(:,:) = 0
            col_iqen(:,:) = 0
            col_iqse(:,:) = 0
            col_iqes(:,:) = 0
            col_iqnw(:,:) = 0
            col_iqwn(:,:) = 0
            col_iqsw(:,:) = 0
            col_iqws(:,:) = 0
            do iq = 1,mg(g)%ifull
               nc = mg_colourmask(iq)
               mg_ifull_colour(nc) = mg_ifull_colour(nc) + 1
               iqq = mg_ifull_colour(nc)
               col_iq(iqq,nc) = iq
               col_iqn(iqq,nc) = mg(g)%in(iq)
               col_iqe(iqq,nc) = mg(g)%ie(iq)
               col_iqs(iqq,nc) = mg(g)%is(iq)
               col_iqw(iqq,nc) = mg(g)%iw(iq)
               col_iqne(iqq,nc) = mg(g)%ine(iq)
               col_iqen(iqq,nc) = mg(g)%ien(iq)
               col_iqse(iqq,nc) = mg(g)%ise(iq)
               col_iqes(iqq,nc) = mg(g)%ies(iq)
               col_iqnw(iqq,nc) = mg(g)%inw(iq)
               col_iqwn(iqq,nc) = mg(g)%iwn(iq)
               col_iqsw(iqq,nc) = mg(g)%isw(iq)
               col_iqws(iqq,nc) = mg(g)%iws(iq)
            end do
         
            deallocate( mg_colourmask )

         end if

      end if
       
   return
   end subroutine mg_index

   subroutine mgcheck_bnds_alloc(g,iproc,iext)

      integer, intent(in) :: iproc
      integer, intent(in) :: g, iext

      if ( mg_bnds(iproc,g)%len <= 0 ) then
         allocate( mg_bnds(iproc,g)%request_list(mg(g)%iextra) )
         allocate( mg_bnds(iproc,g)%unpack_list(mg(g)%iextra) )
         mg_bnds(iproc,g)%len = mg(g)%iextra
      else
         if ( iext>mg(g)%iextra ) then
            write(6,*) "ERROR: MG grid undersized in mgcheck_bnds_alloc"
            write(6,*) "iext,iextra,g,iproc,myid ",iext,mg(g)%iextra,g,iproc,myid
            call ccmpi_abort(-1)
         end if
      end if

   return
   end subroutine mgcheck_bnds_alloc
   
   subroutine mgbounds2(g,vdat,corner)

      integer, intent(in) :: g
      logical, intent(in), optional :: corner
      logical :: extra
      integer :: iproc, recv_len, send_len
      integer :: rcount, jproc, iq
      integer, dimension(mg(g)%neighnum) :: rslen, sslen
      integer(kind=4) :: ierr, itag=20, llen, sreq, lproc
      integer(kind=4) :: ldone, lcomm
      integer(kind=4), dimension(size(ireq)) :: donelist
      integer(kind=4), dimension(MPI_STATUS_SIZE,size(ireq)) :: status
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      real, dimension(:), intent(inout) :: vdat
      
      extra = .false.
      if ( present(corner) ) then
         extra = corner
      end if   
      
      if ( extra ) then
         rslen(1:mg(g)%neighnum) = mg_bnds(mg(g)%neighlist,g)%rlenx
         sslen(1:mg(g)%neighnum) = mg_bnds(mg(g)%neighlist,g)%slenx
      else    
         rslen(1:mg(g)%neighnum) = mg_bnds(mg(g)%neighlist,g)%rlen
         sslen(1:mg(g)%neighnum) = mg_bnds(mg(g)%neighlist,g)%slen
      end if   

      !     Set up the buffers to send and recv
      lcomm = comm_world
      nreq = 0
      do iproc = 1,mg(g)%neighnum
         recv_len = rslen(iproc)
         if ( recv_len > 0 ) then
            lproc = mg(g)%neighlist(iproc)  ! Recv from
            nreq = nreq + 1
            rlist(nreq) = iproc
            llen = recv_len
            call MPI_IRecv( bnds(lproc)%rbuf, llen, ltype, lproc, &
                            itag, lcomm, ireq(nreq), ierr )
         end if
      end do
      rreq = nreq
      do iproc = mg(g)%neighnum,1,-1
         send_len = sslen(iproc)
         if ( send_len > 0 ) then
            lproc = mg(g)%neighlist(iproc)  ! Send to
            do concurrent (iq = 1:send_len)
               bnds(lproc)%sbuf(iq) = vdat(mg_bnds(lproc,g)%send_list(iq))
            end do   
            nreq = nreq + 1
            llen = send_len
            call MPI_ISend( bnds(lproc)%sbuf, llen, ltype, lproc, &
                            itag, lcomm, ireq(nreq), ierr )
         end if
      end do

      rcount = rreq
      do while ( rcount > 0 )
         call START_LOG(mpiwaitmg_begin)
         call MPI_Waitsome( rreq, ireq, ldone, donelist, status, ierr )
         call END_LOG(mpiwaitmg_end)
         rcount = rcount - ldone
         do jproc = 1,ldone
            iproc = rlist(donelist(jproc))  ! Recv from
            lproc = mg(g)%neighlist(iproc)
            recv_len = rslen(iproc)
            do concurrent (iq = 1:recv_len)
               vdat(mg(g)%ifull+mg_bnds(lproc,g)%unpack_list(iq)) &
                   = bnds(lproc)%rbuf(iq)
            end do    
         end do
      end do

      ! clear any remaining messages
      sreq = nreq - rreq
      if ( sreq > 0 ) then
         call START_LOG(mpiwaitmg_begin)
         call MPI_Waitall( sreq, ireq(rreq+1:nreq), status, ierr )
         call END_LOG(mpiwaitmg_end)
      end if   

   return
   end subroutine mgbounds2
 
   subroutine mgbounds3(g,vdat,klim,corner)

      integer, intent(in) :: g
      integer, intent(in), optional :: klim
      logical, intent(in), optional :: corner
      logical :: extra
      integer :: kx, iproc, recv_len, send_len
      integer :: rcount, jproc, iq, k
      integer, dimension(mg(g)%neighnum) :: rslen, sslen
      integer(kind=4) :: ierr, itag=20, llen, sreq, lproc
      integer(kind=4) :: ldone, lcomm
      integer(kind=4), dimension(size(ireq)) :: donelist
      integer(kind=4), dimension(MPI_STATUS_SIZE,size(ireq)) :: status
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      real, dimension(:,:), intent(inout) :: vdat

      kx = size(vdat,2)
      extra = .false.
      if (present(klim)) then
         kx = klim
      end if
      if ( present(corner) ) then
         extra = corner
      end if   
      
      if ( extra ) then
         rslen(1:mg(g)%neighnum) = mg_bnds(mg(g)%neighlist,g)%rlenx
         sslen(1:mg(g)%neighnum) = mg_bnds(mg(g)%neighlist,g)%slenx
      else    
         rslen(1:mg(g)%neighnum) = mg_bnds(mg(g)%neighlist,g)%rlen
         sslen(1:mg(g)%neighnum) = mg_bnds(mg(g)%neighlist,g)%slen
      end if 
      
      rslen(1:mg(g)%neighnum) = mg_bnds(mg(g)%neighlist,g)%rlen
      sslen(1:mg(g)%neighnum) = mg_bnds(mg(g)%neighlist,g)%slen

      !     Set up the buffers to send and recv
      lcomm = comm_world
      nreq = 0
      do iproc = 1,mg(g)%neighnum
         recv_len = rslen(iproc)
         if ( recv_len > 0 ) then
            lproc = mg(g)%neighlist(iproc)  ! Recv from
            nreq = nreq + 1
            rlist(nreq) = iproc
            llen = recv_len*kx
            call MPI_IRecv( bnds(lproc)%rbuf, llen, ltype, lproc, &
                            itag, lcomm, ireq(nreq), ierr )
         end if
      end do
      rreq = nreq
      do iproc = mg(g)%neighnum,1,-1
         send_len = sslen(iproc)
         if ( send_len > 0 ) then
            lproc = mg(g)%neighlist(iproc)  ! Send to
            do k = 1,kx
               do concurrent (iq = 1:send_len)
                  bnds(lproc)%sbuf(iq+(k-1)*send_len) = vdat(mg_bnds(lproc,g)%send_list(iq),k)
               end do
            end do   
            nreq = nreq + 1
            llen = send_len*kx
            call MPI_ISend( bnds(lproc)%sbuf, llen, ltype, lproc, &
                            itag, lcomm, ireq(nreq), ierr )
         end if
      end do

      rcount = rreq
      do while ( rcount > 0 )
         call START_LOG(mpiwaitmg_begin)
         call MPI_Waitsome( rreq, ireq, ldone, donelist, status, ierr )
         call END_LOG(mpiwaitmg_end)
         rcount = rcount - ldone
         do jproc = 1,ldone
            iproc = rlist(donelist(jproc))  ! Recv from
            lproc = mg(g)%neighlist(iproc)
            recv_len = rslen(iproc)
            do k = 1,kx
               do concurrent (iq = 1:recv_len)
                  vdat(mg(g)%ifull+mg_bnds(lproc,g)%unpack_list(iq),k) &
                      = bnds(lproc)%rbuf(iq+(k-1)*recv_len)
               end do
           end do    
         end do
      end do

      ! clear any remaining messages
      sreq = nreq - rreq
      if ( sreq > 0 ) then
         call START_LOG(mpiwaitmg_begin)
         call MPI_Waitall( sreq, ireq(rreq+1:nreq), status, ierr )
         call END_LOG(mpiwaitmg_end)
      end if   

   return
   end subroutine mgbounds3

   pure function mg_fproc_1(g,i,j,n) result(mg_fpout)
     ! locates processor that owns a global grid point
     integer, intent(in) :: i, j, n, g
     integer mg_fpout
     integer g_l, i_l, j_l
     
     i_l = i
     j_l = j
     do g_l = g-1,1,-1
        i_l = (i_l-1)*2 + 1 
        j_l = (j_l-1)*2 + 1 
     end do
     
     mg_fpout = fproc(i_l,j_l,n)
     
   end function mg_fproc_1
   
   pure function mg_fproc(g,i,j,n) result(mg_fpout)
     ! locates processor that owns a global grid point
     integer, intent(in) :: i, j, n, g
     integer mg_fpout
     integer fp_l
     
     fp_l = mg_fproc_1(g,i,j,n)
     mg_fpout = mg(g)%procmap(fp_l)
     
   end function mg_fproc   
   
   pure function mg_qproc(iqg,mil_g,g) result(mg_qpout)
      ! locates processor that owns a global grid point
      integer, intent(in) :: iqg, mil_g, g
      integer :: mg_qpout
      integer :: i, j, n

      n = (iqg - 1) / (mil_g*mil_g)
      j = 1 + (iqg - n*mil_g*mil_g - 1)/mil_g
      i = iqg - (j - 1)*mil_g - n*mil_g*mil_g

      mg_qpout = mg_fproc(g,i,j,n)
   
   end function mg_qproc

   subroutine ccmpi_filewinget2(abuf,sinp)
   
      integer :: w, ncount, nlen, cc, ipf
      integer :: n, is, ie, ipin
      integer(kind=4) :: lsize, ierr
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      integer(kind=4) :: lcomm
      integer(kind=4) :: itag = 50
      integer(kind=4), dimension(size(filemap_recv)+size(filemap_send)) :: i_req
      integer(kind=4), dimension(MPI_STATUS_SIZE,size(filemap_recv)+size(filemap_send)) :: status
      real, dimension(:), intent(in) :: sinp
      real, dimension(pipan*pjpan*pnpan,size(filemap_recv)), intent(out) :: abuf 

      nlen = pipan*pjpan*pnpan
      lsize = nlen
      
      if ( fnproc <= fnproc_bcast_max ) then
        
         ! use Bcast for single input panel
         abuf = 0.
         do ipin = 0,fnproc-1
            lcomm = filemap_facecomm(ipin) 
            n = filemap_rinv(ipin)
            if ( myid == mod(ipin,fnresid) ) then
               ipf = ipin/fnresid
               cc = nlen*ipf
               call START_LOG(bcast_begin)
               call MPI_Bcast( sinp(1+cc:nlen+cc), lsize, ltype, 0_4, lcomm, ierr )
               call END_LOG(bcast_end)
               if ( n > 0 ) then
                  abuf(1:nlen,n) = sinp(1+cc:nlen+cc)
               end if   
            else
               if ( n > 0 ) then 
                  call START_LOG(bcast_begin)
                  call MPI_Bcast( abuf(1:nlen,n), lsize, ltype, 0_4, lcomm, ierr )
                  call END_LOG(bcast_end)
               end if   
            end if   
         end do
          
      else    

        ! MJT notes - We could use RMA or MPI_Alltoallv for this problem.  RMA fence forces synchronisation
        ! across processes.  Some implentations of MPI_Alltoallv also can employ synchronisation as they
        ! are optimsied for most sendcounts and recvcounts are non-zero.  Splitting communicators and using
        ! MPI_Bcast could work, but is expensive to initialise.  
          
         ncount = size(filemap_recv)
         lcomm = comm_world
          
         !     Set up the buffers to recv
         nreq = 0
         do w = 1,ncount
            ipf = filemap_rmod(w)
            itag = 50 + ipf
            nreq  = nreq + 1
            call MPI_IRecv( abuf(:,w), lsize, ltype, filemap_recv(w), itag, lcomm, i_req(nreq), ierr )
         end do
      
         !     Set up the buffers to send
         do w = 1,size(filemap_send)
            ipf = filemap_smod(w)
            itag = 50 + ipf
            cc = nlen*ipf
            nreq  = nreq + 1
            call MPI_ISend( sinp(1+cc:nlen+cc), lsize, ltype, filemap_send(w), itag, lcomm, i_req(nreq), ierr )
         end do
      
         call START_LOG(mpiwaitmapfile_begin) 
         call MPI_Waitall( nreq, i_req, status, ierr )
         call END_LOG(mpiwaitmapfile_end)
         
      end if   
      
   end subroutine ccmpi_filewinget2

   subroutine ccmpi_filewinget3(abuf,sinp)
   
      integer :: n, w, ncount, nlen, kx
      integer :: cc, ipf, ipin
      integer :: rcount, jproc
      integer :: is, ie
      integer(kind=4) :: lsize, ierr
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      integer(kind=4) :: lcomm, sreq, ldone
      integer(kind=4) :: itag = 51
      integer(kind=4), dimension(size(filemap_recv)+size(filemap_send)) :: i_req
      integer(kind=4), dimension(size(filemap_recv)) :: i_list, donelist
      integer(kind=4), dimension(MPI_STATUS_SIZE,size(filemap_recv)+size(filemap_send)) :: status
      real, dimension(:,:), intent(in) :: sinp
      real, dimension(:,:,:), intent(out) :: abuf
      real, dimension(pipan*pjpan*pnpan,size(sinp,2),size(filemap_recv)) :: bbuf
      real, dimension(pipan*pjpan*pnpan,size(sinp,2),mynproc) :: cbuf

      kx = size(sinp,2)
      nlen = pipan*pjpan*pnpan
      lsize = nlen*kx
      
      if ( fnproc <= fnproc_bcast_max ) then
      
         ! use Bcast for single input file  
         abuf = 0.
         do ipin = 0,fnproc-1
            lcomm = filemap_facecomm(ipin) 
            n = filemap_rinv(ipin)
            if ( myid == mod(ipin,fnresid) ) then
               ipf = ipin/fnresid
               cc = nlen*ipf             
               cbuf(:,:,ipf+1) = sinp(1+cc:nlen+cc,1:kx) 
               call START_LOG(bcast_begin)
               call MPI_Bcast( cbuf(:,:,ipf+1), lsize, ltype, 0_4, lcomm, ierr )
               call END_LOG(bcast_end)
               if ( n > 0 ) then
                  abuf(1:nlen,n,1:kx) = cbuf(1:nlen,1:kx,ipf+1) 
               end if   
            else
               if ( n > 0 ) then
                  call START_LOG(bcast_begin) 
                  call MPI_Bcast( bbuf(:,:,n), lsize, ltype, 0_4, lcomm, ierr )
                  call END_LOG(bcast_end)
                  abuf(1:nlen,n,1:kx) = bbuf(1:nlen,1:kx,n)
               end if
            end if
         end do      
         
      else   

         ncount = size(filemap_recv)
         lcomm = comm_world
          
         do ipf = 0,mynproc-1
            cc = nlen*ipf 
            cbuf(1:nlen,1:kx,ipf+1) = sinp(1+cc:nlen+cc,1:kx)
         end do  
      
         !     Set up the buffers to recv
         nreq = 0
         do w = 1,ncount
            ipf = filemap_rmod(w)
            itag = 51 + ipf
            nreq  = nreq + 1
            i_list(nreq) = w
            call MPI_IRecv( bbuf(:,:,w), lsize, ltype, filemap_recv(w), itag, lcomm, i_req(nreq), ierr )
         end do
         rreq = nreq
      
         !     Set up the buffers to send
         do w = 1,size(filemap_send)
            ipf = filemap_smod(w)
            itag = 51 + ipf
            nreq  = nreq + 1
            call MPI_ISend( cbuf(:,:,ipf+1), lsize, ltype, filemap_send(w), itag, lcomm, i_req(nreq), ierr )  
         end do
      
         ! Unpack incomming messages
         rcount = rreq
         do while ( rcount > 0 )
            call START_LOG(mpiwaitmapfile_begin) 
            call MPI_Waitsome( rreq, i_req, ldone, donelist, status, ierr )
            call END_LOG(mpiwaitmapfile_end)
            rcount = rcount - ldone
            do jproc = 1,ldone
               w = i_list(donelist(jproc))
               abuf(1:nlen,w,1:kx) = bbuf(1:nlen,1:kx,w)
            end do
         end do
      
         sreq = nreq - rreq
         if ( sreq > 0 ) then
            call START_LOG(mpiwaitmapfile_begin) 
            call MPI_Waitall( sreq, i_req(rreq+1:nreq), status, ierr )
            call END_LOG(mpiwaitmapfile_end)
         end if
         
      end if   
      
   end subroutine ccmpi_filewinget3
   
   subroutine ccmpi_filewinunpack(sout,abuf)

      integer :: ncount, w, ip, n, no, ca, cb, cc
      real, dimension(-1:,-1:,0:), intent(inout) :: sout
      real, dimension(:,:), intent(in) :: abuf

      ncount = size(filemap_recv)
      
      do w = 1,ncount
         ip = filemap_recv(w) + filemap_rmod(w)*fnresid
         do n = 0,pnpan-1
            no = n - pnoff(ip) + 1
            ca = pioff(ip,no)
            cb = pjoff(ip,no)
            cc = n*pipan*pjpan
            sout(1+ca:pipan+ca,1+cb:pjpan+cb,no) = reshape( abuf(1+cc:pipan*pjpan+cc,w), (/ pipan, pjpan /) )
         end do
      end do
      
   end subroutine ccmpi_filewinunpack
   
   subroutine ccmpi_filebounds_setup(comm)

      integer, intent(in) :: comm
      integer, dimension(:,:), allocatable, save :: dumi
      integer :: ipf, n, i, j, iq, ncount, ca, cb, no, ip
      integer :: filemaxbuflen, xlen
      integer :: iproc, jproc, iloc, jloc, nloc, floc
      integer(kind=4) :: lproc, lcomm, llen, ierr
      integer(kind=4) :: itag=42
      integer(kind=4), dimension(:,:), allocatable :: status
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_INTEGER8
#else
      integer(kind=4), parameter :: ltype = MPI_INTEGER
#endif
      logical, save :: fileallocate = .false.
     
      if ( myid >= fnresid ) return

      lcomm = comm
      filemaxbuflen = 2*(pipan+pjpan+2)*pnpan*fncount
      maxvertlen = max( maxvertlen, pka_g, pko_g )
      
      ! allocate memory to filebnds
      if ( fileallocate ) then
         do iproc = 0,size(filebnds)-1
            if ( filebnds(iproc)%slenx > 0 ) then
               deallocate( filebnds(iproc)%send_list )
            end if
            if ( filebnds(iproc)%rlenx > 0 ) then
               deallocate( filebnds(iproc)%unpack_list )
            end if
         end do
         if ( filebnds(myid)%rlenx > 0 ) then
           deallocate( filebnds(myid)%request_list )
         end if  
         deallocate( filebnds )
      end if
      allocate( filebnds(0:fnresid-1) )
      fileallocate = .true.

      ! calculate recv message length
      filebnds(:)%len = 0
      filebnds(:)%rlen = 0
      filebnds(:)%slen = 0
      filebnds(:)%rlenx = 0
      filebnds(:)%slenx = 0
      do ipf = 0,fncount-1  ! fncount=fnproc/fnresid
         ip = ipf*fnresid + myid
         do n = 1,pnpan
            no = n - pnoff(ip)
            ca = pioff(ip,no)
            cb = pjoff(ip,no)
            do i = 1,pipan
               iproc = mod(procarray(i+ca,cb,no), fnresid) ! South
               floc  = procarray(i+ca,cb,no)/fnresid + 1
               filebnds(iproc)%rlenx = filebnds(iproc)%rlenx + 1
               call check_filebnds_alloc(iproc,filemaxbuflen)
               ! store global index
               filebnds(iproc)%request_list(filebnds(iproc)%rlenx,1) = i + ca
               filebnds(iproc)%request_list(filebnds(iproc)%rlenx,2) = cb
               filebnds(iproc)%request_list(filebnds(iproc)%rlenx,3) = no
               filebnds(iproc)%request_list(filebnds(iproc)%rlenx,4) = floc
               ! store local index
               filebnds(iproc)%unpack_list(filebnds(iproc)%rlenx,1) = i
               filebnds(iproc)%unpack_list(filebnds(iproc)%rlenx,2) = 0
               filebnds(iproc)%unpack_list(filebnds(iproc)%rlenx,3) = n
               filebnds(iproc)%unpack_list(filebnds(iproc)%rlenx,4) = ipf + 1
               iproc = mod(procarray(i+ca,pjpan+1+cb,no), fnresid) ! North
               floc  = procarray(i+ca,pjpan+1+cb,no)/fnresid + 1
               filebnds(iproc)%rlenx = filebnds(iproc)%rlenx + 1
               call check_filebnds_alloc(iproc,filemaxbuflen)
               ! store global index
               filebnds(iproc)%request_list(filebnds(iproc)%rlenx,1) = i + ca
               filebnds(iproc)%request_list(filebnds(iproc)%rlenx,2) = pjpan + 1 + cb
               filebnds(iproc)%request_list(filebnds(iproc)%rlenx,3) = no
               filebnds(iproc)%request_list(filebnds(iproc)%rlenx,4) = floc
               ! store local index
               filebnds(iproc)%unpack_list(filebnds(iproc)%rlenx,1) = i
               filebnds(iproc)%unpack_list(filebnds(iproc)%rlenx,2) = pjpan + 1
               filebnds(iproc)%unpack_list(filebnds(iproc)%rlenx,3) = n
               filebnds(iproc)%unpack_list(filebnds(iproc)%rlenx,4) = ipf + 1
            end do
            do j = 1,pjpan
               iproc = mod(procarray(ca,j+cb,no), fnresid) ! West
               floc  = procarray(ca,j+cb,no)/fnresid + 1
               filebnds(iproc)%rlenx = filebnds(iproc)%rlenx + 1
               call check_filebnds_alloc(iproc,filemaxbuflen)
               ! store global index
               filebnds(iproc)%request_list(filebnds(iproc)%rlenx,1) = ca
               filebnds(iproc)%request_list(filebnds(iproc)%rlenx,2) = j + cb
               filebnds(iproc)%request_list(filebnds(iproc)%rlenx,3) = no
               filebnds(iproc)%request_list(filebnds(iproc)%rlenx,4) = floc
               ! store local index
               filebnds(iproc)%unpack_list(filebnds(iproc)%rlenx,1) = 0
               filebnds(iproc)%unpack_list(filebnds(iproc)%rlenx,2) = j
               filebnds(iproc)%unpack_list(filebnds(iproc)%rlenx,3) = n
               filebnds(iproc)%unpack_list(filebnds(iproc)%rlenx,4) = ipf + 1
               iproc = mod(procarray(pipan+1+ca,j+cb,no), fnresid) ! East
               floc  = procarray(pipan+1+ca,j+cb,no)/fnresid + 1
               filebnds(iproc)%rlenx = filebnds(iproc)%rlenx + 1
               call check_filebnds_alloc(iproc,filemaxbuflen)
               ! store global index
               filebnds(iproc)%request_list(filebnds(iproc)%rlenx,1) = pipan + 1 + ca
               filebnds(iproc)%request_list(filebnds(iproc)%rlenx,2) = j + cb
               filebnds(iproc)%request_list(filebnds(iproc)%rlenx,3) = no
               filebnds(iproc)%request_list(filebnds(iproc)%rlenx,4) = floc
               ! store local index
               filebnds(iproc)%unpack_list(filebnds(iproc)%rlenx,1) = pipan + 1
               filebnds(iproc)%unpack_list(filebnds(iproc)%rlenx,2) = j
               filebnds(iproc)%unpack_list(filebnds(iproc)%rlenx,3) = n
               filebnds(iproc)%unpack_list(filebnds(iproc)%rlenx,4) = ipf + 1
            end do
         end do
      end do
      
      filebnds(:)%rlen = filebnds(:)%rlenx
      do ipf = 0,fncount-1  ! fncount=fnproc/fnresid
         ip = ipf*fnresid + myid
         do n = 1,pnpan
            no = n - pnoff(ip)
            ca = pioff(ip,no)
            cb = pjoff(ip,no)
            iproc = mod(procarray(ca,cb,no), fnresid) ! South - West
            floc  = procarray(ca,cb,no)/fnresid + 1
            filebnds(iproc)%rlenx = filebnds(iproc)%rlenx + 1
            call check_filebnds_alloc(iproc,filemaxbuflen)
            ! store global index
            filebnds(iproc)%request_list(filebnds(iproc)%rlenx,1) = ca
            filebnds(iproc)%request_list(filebnds(iproc)%rlenx,2) = cb
            filebnds(iproc)%request_list(filebnds(iproc)%rlenx,3) = no
            filebnds(iproc)%request_list(filebnds(iproc)%rlenx,4) = floc
            ! store local index
            filebnds(iproc)%unpack_list(filebnds(iproc)%rlenx,1) = 0
            filebnds(iproc)%unpack_list(filebnds(iproc)%rlenx,2) = 0
            filebnds(iproc)%unpack_list(filebnds(iproc)%rlenx,3) = n
            filebnds(iproc)%unpack_list(filebnds(iproc)%rlenx,4) = ipf + 1
            iproc = mod(procarray(ca,pjpan+1+cb,no), fnresid) ! North - West
            floc  = procarray(ca,pjpan+1+cb,no)/fnresid + 1
            filebnds(iproc)%rlenx = filebnds(iproc)%rlenx + 1
            call check_filebnds_alloc(iproc,filemaxbuflen)
            ! store global index
            filebnds(iproc)%request_list(filebnds(iproc)%rlenx,1) = ca
            filebnds(iproc)%request_list(filebnds(iproc)%rlenx,2) = pjpan + 1 + cb
            filebnds(iproc)%request_list(filebnds(iproc)%rlenx,3) = no
            filebnds(iproc)%request_list(filebnds(iproc)%rlenx,4) = floc
            ! store local index
            filebnds(iproc)%unpack_list(filebnds(iproc)%rlenx,1) = 0
            filebnds(iproc)%unpack_list(filebnds(iproc)%rlenx,2) = pjpan + 1
            filebnds(iproc)%unpack_list(filebnds(iproc)%rlenx,3) = n
            filebnds(iproc)%unpack_list(filebnds(iproc)%rlenx,4) = ipf + 1
            iproc = mod(procarray(pipan+1+ca,cb,no), fnresid) ! South - East
            floc  = procarray(pipan+1+ca,cb,no)/fnresid + 1
            filebnds(iproc)%rlenx = filebnds(iproc)%rlenx + 1
            call check_filebnds_alloc(iproc,filemaxbuflen)
            ! store global index
            filebnds(iproc)%request_list(filebnds(iproc)%rlenx,1) = pipan + 1 + ca
            filebnds(iproc)%request_list(filebnds(iproc)%rlenx,2) = cb
            filebnds(iproc)%request_list(filebnds(iproc)%rlenx,3) = no
            filebnds(iproc)%request_list(filebnds(iproc)%rlenx,4) = floc
            ! store local index
            filebnds(iproc)%unpack_list(filebnds(iproc)%rlenx,1) = pipan + 1
            filebnds(iproc)%unpack_list(filebnds(iproc)%rlenx,2) = 0
            filebnds(iproc)%unpack_list(filebnds(iproc)%rlenx,3) = n
            filebnds(iproc)%unpack_list(filebnds(iproc)%rlenx,4) = ipf + 1
            iproc = mod(procarray(pipan+1+ca,pjpan+1+cb,no), fnresid) ! North - East
            floc  = procarray(pipan+1+ca,pjpan+1+cb,no)/fnresid + 1
            filebnds(iproc)%rlenx = filebnds(iproc)%rlenx + 1
            call check_filebnds_alloc(iproc,filemaxbuflen)
            ! store global index
            filebnds(iproc)%request_list(filebnds(iproc)%rlenx,1) = pipan + 1 + ca
            filebnds(iproc)%request_list(filebnds(iproc)%rlenx,2) = pjpan + 1 + cb
            filebnds(iproc)%request_list(filebnds(iproc)%rlenx,3) = no
            filebnds(iproc)%request_list(filebnds(iproc)%rlenx,4) = floc
            ! store local index
            filebnds(iproc)%unpack_list(filebnds(iproc)%rlenx,1) = pipan + 1
            filebnds(iproc)%unpack_list(filebnds(iproc)%rlenx,2) = pjpan + 1
            filebnds(iproc)%unpack_list(filebnds(iproc)%rlenx,3) = n
            filebnds(iproc)%unpack_list(filebnds(iproc)%rlenx,4) = ipf + 1
         end do
      end do
 
      ! identify neighbour processors
      fileneighnum = count( filebnds(:)%rlenx > 0 )
      if ( filebnds(myid)%rlenx > 0 ) then
        fileneighnum = fileneighnum - 1
      end if
      if ( allocated(fileneighlist) ) then
         deallocate(fileneighlist)
      end if
      allocate(fileneighlist(fileneighnum))
      ncount = 0
      do jproc = 1,fnresid-1
         iproc = modulo(myid+jproc,fnresid)
          if ( filebnds(iproc)%rlenx > 0 ) then
            ncount = ncount + 1
            fileneighlist(ncount) = iproc
         end if
      end do

      ! increase size of request list if needed
      if ( 2*fileneighnum > size(ireq) ) then
         deallocate( ireq, rlist )
         allocate( ireq(max(2*fileneighnum,1)) )
         allocate( rlist(max(fileneighnum,1)) )
      end if
      allocate( status(MPI_STATUS_SIZE,size(ireq)) )

      ! Now, for each processor send the length of points I want.
      allocate( dumi(4,fileneighnum) )
      nreq = 0
      do jproc = 1,fileneighnum
         iproc = fileneighlist(jproc) ! Recv from
         nreq = nreq + 1
         lproc = iproc
         call MPI_IRecv( dumi(1:2,jproc), 2_4, ltype, lproc, itag, lcomm, ireq(nreq), ierr )
      end do
      do jproc = fileneighnum,1,-1
         iproc = fileneighlist(jproc)
         nreq = nreq + 1
         lproc = iproc
         dumi(3,jproc) = filebnds(iproc)%rlen
         dumi(4,jproc) = filebnds(iproc)%rlenx
         call MPI_ISend( dumi(3:4,jproc), 2_4, ltype, lproc, itag, lcomm, ireq(nreq), ierr )
      end do
      call MPI_Waitall( nreq, ireq, status, ierr )
      do jproc = 1,fileneighnum
         iproc = fileneighlist(jproc)
         filebnds(iproc)%slen = dumi(1,jproc)
         filebnds(iproc)%slenx = dumi(2,jproc)
      end do
      deallocate( dumi )
      
      ! send unpack_list to file neighbours
      do i = 1,4
         nreq = 0
         do jproc = 1,fileneighnum
            iproc = fileneighlist(jproc)  ! Recv from
            if ( i == 1 ) then
               allocate(filebnds(iproc)%send_list(filebnds(iproc)%slenx,4))
            end if   
            nreq = nreq + 1
            ! Recv list of requests
            llen = filebnds(iproc)%slenx
            lproc = iproc
            call MPI_IRecv( filebnds(iproc)%send_list(:,i), llen, &
                 ltype, lproc, itag, lcomm, ireq(nreq), ierr )
         end do
         do jproc = fileneighnum,1,-1
            iproc = fileneighlist(jproc)  ! Send to
            ! Send list of requests
            nreq = nreq + 1
            llen = filebnds(iproc)%rlenx
            lproc = iproc
            call MPI_ISend( filebnds(iproc)%request_list(:,i), llen, &
                 ltype, lproc, itag, lcomm, ireq(nreq), ierr )
         end do      
         call MPI_Waitall( nreq, ireq, status, ierr )
      end do
      
      deallocate( status )
      

      ! convert send_list and unpack_list to local indices
      allocate( dumi(filemaxbuflen,4) )
      do jproc = 1,fileneighnum
         iproc = fileneighlist(jproc)  ! Send to
         deallocate( filebnds(iproc)%request_list )
         dumi(1:filebnds(iproc)%rlenx,1:4) = filebnds(iproc)%unpack_list(1:filebnds(iproc)%rlenx,1:4)
         deallocate( filebnds(iproc)%unpack_list )
         allocate( filebnds(iproc)%unpack_list(filebnds(iproc)%rlenx,4) )
         filebnds(iproc)%unpack_list(1:filebnds(iproc)%rlenx,1:4) = dumi(1:filebnds(iproc)%rlenx,1:4)
         do iq = 1,filebnds(iproc)%slenx
            iloc = filebnds(iproc)%send_list(iq,1)
            jloc = filebnds(iproc)%send_list(iq,2)
            nloc = filebnds(iproc)%send_list(iq,3)
            floc = filebnds(iproc)%send_list(iq,4)
            call file_ijnpg2ijnp(iloc,jloc,nloc,floc,myid,pil_g)
            filebnds(iproc)%send_list(iq,1) = iloc
            filebnds(iproc)%send_list(iq,2) = jloc
            filebnds(iproc)%send_list(iq,3) = nloc
            filebnds(iproc)%send_list(iq,4) = floc
         end do
      end do
      deallocate( dumi )
      do iq = 1,filebnds(myid)%rlenx
         iloc = filebnds(myid)%request_list(iq,1)
         jloc = filebnds(myid)%request_list(iq,2)
         nloc = filebnds(myid)%request_list(iq,3)
         floc = filebnds(myid)%request_list(iq,4)
         call file_ijnpg2ijnp(iloc,jloc,nloc,floc,myid,pil_g)
         filebnds(myid)%request_list(iq,1) = iloc
         filebnds(myid)%request_list(iq,2) = jloc
         filebnds(myid)%request_list(iq,3) = nloc
         filebnds(myid)%request_list(iq,4) = floc
      end do
   
      ! set up buffers
      do jproc = 1,fileneighnum
         iproc = fileneighlist(jproc)
         xlen = maxvertlen*filebnds(iproc)%rlenx
         if ( bnds(iproc)%rbuflen < xlen ) then
            if ( allocated(bnds(iproc)%rbuf) ) then
               deallocate( bnds(iproc)%rbuf )
               deallocate( bnds(iproc)%r8buf )
            end if
            allocate( bnds(iproc)%rbuf(xlen) )
            allocate( bnds(iproc)%r8buf(xlen) )
            bnds(iproc)%rbuflen = xlen
         end if
         xlen = maxvertlen*filebnds(iproc)%slenx
         if ( bnds(iproc)%sbuflen < xlen ) then
            if ( allocated(bnds(iproc)%sbuf) ) then
               deallocate( bnds(iproc)%sbuf )
               deallocate( bnds(iproc)%s8buf )
            end if
            allocate( bnds(iproc)%sbuf(xlen) )
            allocate( bnds(iproc)%s8buf(xlen) )
            bnds(iproc)%sbuflen = xlen
         end if
      end do
      
   end subroutine ccmpi_filebounds_setup
   
   subroutine check_filebnds_alloc(iproc,filemaxbuflen)
   
      integer, intent(in) :: iproc, filemaxbuflen
      
      if ( filebnds(iproc)%len<=0 ) then
         filebnds(iproc)%len = filemaxbuflen
         allocate(filebnds(iproc)%request_list(filemaxbuflen,4))
         allocate(filebnds(iproc)%unpack_list(filemaxbuflen,4))
      else
         if ( filebnds(iproc)%rlenx > filemaxbuflen ) then
            write(6,*) "ERROR: file grid is undersized in check_filebnds_alloc"
            write(6,*) "myid,iproc,rlenx,filemaxbuflen ",myid,iproc,filebnds(iproc)%rlenx,filemaxbuflen
            call ccmpi_abort(-1)
         end if
      end if
   
   end subroutine check_filebnds_alloc

   subroutine ccmpi_filebounds_send2(sdat,comm,corner)

      integer, intent(in) :: comm
      integer :: iproc, iq, send_len
      integer, dimension(fileneighnum) :: rslen, sslen
      integer(kind=4) :: llen, lproc, ierr, lcomm
      integer(kind=4) :: itag=40
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      real, dimension(0:pipan+1,0:pjpan+1,pnpan,1:fncount), intent(inout) :: sdat
      logical, intent(in), optional :: corner
      logical :: extra
      
      lcomm = comm
      
      if ( present(corner) ) then
         extra = corner
      else
         extra = .false.
      end if   
          
      if ( extra ) then
         rslen(:) = filebnds(fileneighlist)%rlenx
         sslen(:) = filebnds(fileneighlist)%slenx
      else    
         rslen(:) = filebnds(fileneighlist)%rlen
         sslen(:) = filebnds(fileneighlist)%slen
      end if 

      !     Set up the buffers to recv
      nreq = 0
      do iproc = 1,fileneighnum
         if ( rslen(iproc) > 0 ) then 
            llen = rslen(iproc)
            lproc = fileneighlist(iproc)  ! Recv from
            nreq  = nreq + 1
            rlist(nreq) = iproc
            call MPI_IRecv( bnds(lproc)%rbuf, llen, ltype, lproc, itag, lcomm, ireq(nreq), ierr )
         end if   
      end do
      rreq = nreq
      !     Set up the buffers to send
      do iproc = fileneighnum,1,-1
         if ( sslen(iproc) > 0 ) then 
            ! Build up list of points
            send_len = sslen(iproc)
            llen = send_len
            lproc = fileneighlist(iproc)  ! Send to
            do concurrent (iq = 1:send_len)
               bnds(lproc)%sbuf(iq) = sdat(filebnds(lproc)%send_list(iq,1),filebnds(lproc)%send_list(iq,2), &
                                           filebnds(lproc)%send_list(iq,3),filebnds(lproc)%send_list(iq,4))
            end do
            nreq  = nreq + 1
            call MPI_ISend( bnds(lproc)%sbuf, llen, ltype, lproc, itag, lcomm, ireq(nreq), ierr )
         end if   
      end do

   end subroutine ccmpi_filebounds_send2

   subroutine ccmpi_filebounds_recv2(sdat,comm,corner)

      integer, intent(in) :: comm
      integer :: myrlen, iproc, jproc, mproc, iq, rcount
      integer, dimension(fileneighnum) :: rslen
      integer(kind=4) :: lproc, ierr, ldone, sreq, lcomm
      integer(kind=4), dimension(fileneighnum) :: donelist
      integer(kind=4), dimension(MPI_STATUS_SIZE,size(ireq)) :: status
      real, dimension(0:pipan+1,0:pjpan+1,pnpan,1:fncount), intent(inout) :: sdat
      logical, intent(in), optional :: corner
      logical :: extra

      lcomm = comm
      
      if ( present(corner) ) then
         extra = corner
      else
         extra = .false.
      end if   
          
      if ( extra ) then
         rslen(:) = filebnds(fileneighlist)%rlenx
      else    
         rslen(:) = filebnds(fileneighlist)%rlen
      end if 
      myrlen = filebnds(myid)%rlenx

      ! See if there are any points on my own processor that need
      ! to be fixed up.
      do concurrent (iq = 1:myrlen)
         ! request_list is same as send_list in this case
         sdat(filebnds(myid)%unpack_list(iq,1),filebnds(myid)%unpack_list(iq,2),   &
              filebnds(myid)%unpack_list(iq,3),filebnds(myid)%unpack_list(iq,4)) = &
         sdat(filebnds(myid)%request_list(iq,1),filebnds(myid)%request_list(iq,2), &
              filebnds(myid)%request_list(iq,3),filebnds(myid)%request_list(iq,4))
      end do

      ! Unpack incomming messages
      rcount = rreq
      do while ( rcount > 0 )
         call START_LOG(mpiwaitfile_begin)
         call MPI_Waitsome( rreq, ireq, ldone, donelist, status, ierr )
         call END_LOG(mpiwaitfile_end)
         rcount = rcount - ldone
         do jproc = 1,ldone
            mproc = donelist(jproc)
            iproc = rlist(mproc)  ! Recv from
            lproc = fileneighlist(iproc)
            ! unpack_list(iq) is index into extended region
            do concurrent (iq = 1:rslen(iproc))
               sdat(filebnds(lproc)%unpack_list(iq,1),filebnds(lproc)%unpack_list(iq,2),   &
                    filebnds(lproc)%unpack_list(iq,3),filebnds(lproc)%unpack_list(iq,4)) = &
               bnds(lproc)%rbuf(iq)
            end do
         end do
      end do

      ! Clear any remaining messages
      sreq = nreq - rreq
      if ( sreq > 0 ) then
         call START_LOG(mpiwaitfile_begin)
         call MPI_Waitall( sreq, ireq(rreq+1:nreq), status, ierr )
         call END_LOG(mpiwaitfile_end)
      end if   

   end subroutine ccmpi_filebounds_recv2
   
   subroutine ccmpi_filebounds_send3(sdat,comm,corner)

      integer, intent(in) :: comm
      integer :: iproc, iq, kx, send_len, k
      integer, dimension(fileneighnum) :: rslen, sslen
      integer(kind=4) :: llen, lproc, ierr, lcomm
      integer(kind=4) :: itag=41
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      real, dimension(0:,0:,1:,1:,1:), intent(inout) :: sdat
      logical, intent(in), optional :: corner
      logical :: extra
      
      kx = size(sdat,5)
      lcomm = comm

      if ( present(corner) ) then
         extra = corner
      else
         extra = .false.
      end if   
          
      if ( extra ) then
         rslen(:) = filebnds(fileneighlist)%rlenx
         sslen(:) = filebnds(fileneighlist)%slenx
      else    
         rslen(:) = filebnds(fileneighlist)%rlen
         sslen(:) = filebnds(fileneighlist)%slen
      end if 

      !     Set up the buffers to send and recv
      nreq = 0
      do iproc = 1,fileneighnum
         if ( rslen(iproc)>0 ) then
            llen = rslen(iproc)*kx
            lproc = fileneighlist(iproc)  ! Recv from
            nreq = nreq + 1
            rlist(nreq) = iproc
            call MPI_IRecv( bnds(lproc)%rbuf, llen, ltype, lproc, itag, lcomm, ireq(nreq), ierr )
         end if   
      end do
      rreq = nreq
      do iproc = fileneighnum,1,-1
         if ( sslen(iproc) > 0 ) then 
            send_len = sslen(iproc)
            llen = send_len*kx
            lproc = fileneighlist(iproc)  ! Send to
            do k = 1,kx
               do concurrent (iq = 1:send_len)
                  bnds(lproc)%sbuf(iq+(k-1)*send_len) = sdat(filebnds(lproc)%send_list(iq,1),filebnds(lproc)%send_list(iq,2),      &
                                                             filebnds(lproc)%send_list(iq,3),filebnds(lproc)%send_list(iq,4),k)
               end do
            end do   
            nreq = nreq + 1
            call MPI_ISend( bnds(lproc)%sbuf, llen, ltype, lproc, itag, lcomm, ireq(nreq), ierr )
         end if   
      end do

   end subroutine ccmpi_filebounds_send3

   subroutine ccmpi_filebounds_recv3(sdat,comm,corner)

      integer, intent(in) :: comm
      integer :: myrlen, iproc, jproc, mproc, iq, rcount, kx
      integer :: k
      integer, dimension(fileneighnum) :: rslen
      integer(kind=4) :: lproc, ierr, ldone, sreq, lcomm
      integer(kind=4), dimension(fileneighnum) :: donelist
      integer(kind=4), dimension(MPI_STATUS_SIZE,size(ireq)) :: status
      real, dimension(0:,0:,1:,1:,1:), intent(inout) :: sdat
      logical, intent(in), optional :: corner
      logical :: extra

      kx = size(sdat,5)
      lcomm = comm

      if ( present(corner) ) then
         extra = corner
      else
         extra = .false.
      end if   
          
      if ( extra ) then
         rslen(:) = filebnds(fileneighlist)%rlenx
      else    
         rslen(:) = filebnds(fileneighlist)%rlen
      end if 
      myrlen = filebnds(myid)%rlenx

      ! See if there are any points on my own processor that need
      ! to be fixed up.
      do k = 1,kx
         do concurrent (iq = 1:myrlen)
            ! request_list is same as send_list in this case
            sdat(filebnds(myid)%unpack_list(iq,1),filebnds(myid)%unpack_list(iq,2),     &
                 filebnds(myid)%unpack_list(iq,3),filebnds(myid)%unpack_list(iq,4),k) = &
            sdat(filebnds(myid)%request_list(iq,1),filebnds(myid)%request_list(iq,2),   &
                 filebnds(myid)%request_list(iq,3),filebnds(myid)%request_list(iq,4),k)
         end do   
      end do

      ! Unpack incomming messages
      rcount = rreq
      do while ( rcount > 0 )
         call START_LOG(mpiwaitfile_begin)
         call MPI_Waitsome( rreq, ireq, ldone, donelist, status, ierr )
         call END_LOG(mpiwaitfile_end)
         rcount = rcount - ldone
         do jproc = 1,ldone
            mproc = donelist(jproc)
            iproc = rlist(mproc)  ! Recv from
            lproc = fileneighlist(iproc)
            do k = 1,kx
               do concurrent (iq = 1:rslen(iproc))
                  sdat(filebnds(lproc)%unpack_list(iq,1),filebnds(lproc)%unpack_list(iq,2),     &
                       filebnds(lproc)%unpack_list(iq,3),filebnds(lproc)%unpack_list(iq,4),k) = &
                  bnds(lproc)%rbuf(iq+(k-1)*rslen(iproc))
               end do
            end do   
         end do
      end do

      ! Clear any remaining messages
      sreq = nreq - rreq
      if ( sreq > 0 ) then
         call START_LOG(mpiwaitfile_begin)
         call MPI_Waitall( sreq, ireq(rreq+1:nreq), status, ierr )
         call END_LOG(mpiwaitfile_end)
      end if   
      
   end subroutine ccmpi_filebounds_recv3
   
   subroutine host_filegather2(a,ag,comm)
      ! Collect global arrays.
      integer, intent(in) :: comm
      integer :: iproc, ipf, ip, n, j
      integer :: ca, cb
#ifdef i8r8
      integer(kind=4),parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4),parameter :: ltype = MPI_REAL
#endif
      integer(kind=4) :: ierr, lsize, lcomm
      real, dimension(pipan*pjpan*pnpan*fncount), intent(in) :: a
      real, dimension(pil_g*pjl_g), intent(out) :: ag
      real, dimension(pipan*pjpan*pnpan*fncount,0:fnresid-1) :: abuf

      lsize = pipan*pjpan*pnpan*fncount
      lcomm = comm
      call START_LOG(gatherfile_begin)
      call MPI_Gather(a,lsize,ltype,abuf,lsize,ltype,0_4,lcomm,ierr)
      call END_LOG(gatherfile_end)

      ! map array in order of processor rank
      do iproc = 0,fnresid-1
         do ipf = 0,fncount-1
            ip = iproc + ipf*fnresid
            do n = 0,pnpan-1
               do j = 1,pjpan
                  ca = (j-1)*pipan + n*pipan*pjpan + ipf*pipan*pjpan*pnpan
                  cb = pioff(ip,n) + (j+pjoff(ip,n)-1)*pil_g + (n-pnoff(ip)+1)*pil_g*pil_g
                  ag(1+cb:pipan+cb) = abuf(1+ca:pipan+ca,iproc)
               end do
            end do
         end do
      end do   

   end subroutine host_filegather2
   
   subroutine proc_filegather2(a,comm)
      ! Collect global arrays.
      integer, intent(in) :: comm
#ifdef i8r8
      integer(kind=4),parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4),parameter :: ltype = MPI_REAL
#endif
      integer(kind=4) :: ierr, lsize, lcomm
      real, dimension(pipan*pjpan*pnpan*fncount), intent(in) :: a
      real, dimension(1,1) :: abuf

      if ( myid == 0 ) then
         write(6,*) "Error: ccmpi_filegather argument is required on proc 0"
         call ccmpi_abort(-1)
      end if
      
      lsize = pipan*pjpan*pnpan*fncount
      lcomm = comm
      call START_LOG(gatherfile_begin)
      call MPI_Gather(a,lsize,ltype,abuf,lsize,ltype,0_4,lcomm,ierr)
      call END_LOG(gatherfile_end)

   end subroutine proc_filegather2

   subroutine host_filegather3(a,ag,comm)
      ! Collect global arrays.
      integer, intent(in) :: comm
      integer :: iproc, k, ipf, ip, n, j
      integer :: kx, ca, cb
#ifdef i8r8
      integer(kind=4),parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4),parameter :: ltype = MPI_REAL
#endif
      integer(kind=4) :: ierr, lsize, lcomm
      real, dimension(:,:), intent(in) :: a
      real, dimension(:,:), intent(out) :: ag
      real, dimension(size(a,1),size(a,2),0:fnresid-1) :: abuf

      kx = size(a,2)
      lsize = pipan*pjpan*pnpan*fncount*kx
      lcomm = comm
      call START_LOG(gatherfile_begin)
      call MPI_Gather(a,lsize,ltype,abuf,lsize,ltype,0_4,lcomm,ierr)
      call END_LOG(gatherfile_end)

      ! map array in order of processor rank
      do iproc = 0,fnresid-1
         do k = 1,kx 
            do ipf = 0,fncount-1
               ip = iproc + ipf*fnresid
               do n = 0,pnpan-1
                  do j = 1,pjpan
                     ca = (j-1)*pipan + n*pipan*pjpan + ipf*pipan*pjpan*pnpan
                     cb = pioff(ip,n) + (j+pjoff(ip,n)-1)*pil_g + (n-pnoff(ip)+1)*pil_g*pil_g
                     ag(1+cb:pipan+cb,k) = abuf(1+ca:pipan+ca,k,iproc)
                  end do   
               end do
            end do
         end do
      end do   

   end subroutine host_filegather3
   
   subroutine proc_filegather3(a,comm)
      ! Collect global arrays.
      integer, intent(in) :: comm
      integer :: kx
#ifdef i8r8
      integer(kind=4),parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4),parameter :: ltype = MPI_REAL
#endif
      integer(kind=4) :: ierr, lsize, lcomm
      real, dimension(:,:), intent(in) :: a
      real, dimension(1,1,1) :: abuf

      if ( myid == 0 ) then
         write(6,*) "Error: ccmpi_filegather argument is required on proc 0"
         call ccmpi_abort(-1)
      end if
      
      kx = size(a,2)
      lsize = pipan*pjpan*pnpan*fncount*kx
      lcomm = comm
      call START_LOG(gatherfile_begin)
      call MPI_Gather(a,lsize,ltype,abuf,lsize,ltype,0_4,lcomm,ierr)
      call END_LOG(gatherfile_end)

   end subroutine proc_filegather3
   
   subroutine host_filedistribute2(af,a1,comm)
      ! Convert standard 1D arrays to face form and distribute to processors
      integer, intent(in) :: comm
      integer :: j, n, iproc, ipf, ip
      integer :: ca, cb
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      integer(kind=4) :: ierr, lsize, lcomm
      real, dimension(pipan*pjpan*pnpan*fncount), intent(out) :: af
      real, dimension(pil_g*pjl_g), intent(in) :: a1
      real, dimension(pipan*pjpan*pnpan*fncount,0:fnresid-1) :: sbuf

      ! map array in order of processor rank
      do iproc = 0,fnresid-1
         do ipf = 0,fncount-1
            ip = iproc + ipf*fnresid
            do n = 0,pnpan-1
               do j = 1,pjpan
                  ca = (j-1)*pipan + n*pipan*pjpan + ipf*pipan*pjpan*pnpan
                  cb = pioff(ip,n) + (j+pjoff(ip,n)-1)*pil_g + (n-pnoff(ip)+1)*pil_g*pil_g
                  sbuf(1+ca:pipan+ca,iproc) = a1(1+cb:pipan+cb)
               end do
            end do
         end do
      end do 

      lsize = pipan*pjpan*pnpan*fncount
      lcomm = comm
      call START_LOG(scatterfile_begin)
      call MPI_Scatter( sbuf, lsize, ltype, af, lsize, ltype, 0_4, lcomm, ierr )
      call END_LOG(scatterfile_end)
      
   end subroutine host_filedistribute2

   subroutine proc_filedistribute2(af,comm)
      ! Convert standard 1D arrays to face form and distribute to processors
      integer, intent(in) :: comm
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      integer(kind=4) :: ierr, lsize, lcomm
      real, dimension(pipan*pjpan*pnpan*fncount), intent(out) :: af
      real, dimension(1,1) :: sbuf
      
      if ( myid == 0 ) then
         write(6,*) "Error: ccmpi_distribute argument required on proc 0"
         call ccmpi_abort(-1)
      end if
      
      lsize = pipan*pjpan*pnpan*fncount
      lcomm = comm
      call START_LOG(scatterfile_begin)
      call MPI_Scatter(sbuf,lsize,ltype,af,lsize,ltype,0_4,lcomm,ierr)
      call END_LOG(scatterfile_end)
      
   end subroutine proc_filedistribute2

   subroutine host_filedistribute3(af,a1,comm)
      ! Convert standard 1D arrays to face form and distribute to processors
      integer, intent(in) :: comm
      integer :: j, n, iproc, k, ipf, ip
      integer :: kx, ca, cb
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      integer(kind=4) :: ierr, lsize, lcomm
      real, dimension(:,:), intent(out) :: af
      real, dimension(:,:), intent(in) :: a1
      real, dimension(size(af,1),size(af,2),0:fnresid-1) :: sbuf
      
      kx = size(af,2)
      
      ! map array in order of processor rank
      do iproc = 0,fnresid-1
         do k = 1,kx 
            do ipf = 0,fncount-1
               ip = iproc + ipf*fnresid
               do n = 0,pnpan-1
                  do j = 1,pjpan
                     ca = (j-1)*pipan + n*pipan*pjpan + ipf*pipan*pjpan*pnpan
                     cb = pioff(ip,n) + (j+pjoff(ip,n)-1)*pil_g + (n-pnoff(ip)+1)*pil_g*pil_g
                     sbuf(1+ca:pipan+ca,k,iproc) = a1(1+cb:pipan+cb,k)
                  end do   
               end do
            end do
         end do
      end do 

      lsize = pipan*pjpan*pnpan*fncount*kx
      lcomm = comm
      call START_LOG(scatterfile_begin)
      call MPI_Scatter( sbuf, lsize, ltype, af, lsize, ltype, 0_4, lcomm, ierr )
      call END_LOG(scatterfile_end)
      
   end subroutine host_filedistribute3

   subroutine proc_filedistribute3(af,comm)
      ! Convert standard 1D arrays to face form and distribute to processors
      integer, intent(in) :: comm
      integer :: kx
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      integer(kind=4) :: ierr, lsize, lcomm
      real, dimension(:,:), intent(out) :: af
      real, dimension(1,1,1) :: sbuf

      if ( myid == 0 ) then
         write(6,*) "Error: ccmpi_distribute argument required on proc 0"
         call ccmpi_abort(-1)
      end if
      
      kx = size(af,2)
      lsize = pipan*pjpan*pnpan*fncount*kx
      lcomm = comm
      call START_LOG(scatterfile_begin)
      call MPI_Scatter(sbuf,lsize,ltype,af,lsize,ltype,0_4,lcomm,ierr)
      call END_LOG(scatterfile_end)
      
   end subroutine proc_filedistribute3
   
   pure subroutine procarray_ijn(i,j,n)
      ! converts (-2:il+2,-2:il+2,0:5) coordinate to (1:il,1:il,0:5)
      integer, intent(inout) :: i, j, n
      integer :: i_in, j_in, n_in, n_n, n_s, n_e, n_w
      
      i_in = i
      j_in = j
      n_in = n

      ! adjust grid index for halo points
      if ( mod(n_in,2) == 0 ) then
         n_w = mod(n_in+5, 6)
         n_e = mod(n_in+2, 6)
         n_n = mod(n_in+1, 6)
         n_s = mod(n_in+4, 6)
         if ( i_in == -1 ) then
            if ( j_in == 0 ) then ! wws
               i = pil_g
               j = 2
               n = n_w
            else if ( j_in == pil_g+1 ) then ! wwn
               i = pil_g
               j = pil_g - 1
               n = n_w
            else
               i = pil_g - 1
               j = j_in
               n = n_w
            end if
         else if ( i_in == 0 ) then
            if ( j_in == -1 ) then ! wss 
               i = pil_g
               j = pil_g - 1
               n = n_s
            else if ( j_in == 0 ) then ! ws
               i = pil_g
               j = 1
               n = n_w
            else if ( j_in == pil_g+1 ) then ! wn
               i = pil_g
               j = pil_g
               n = n_w
            else if ( j_in == pil_g+2 ) then ! wwn
               i = pil_g - 1
               j = pil_g
               n = n_w
            else
               i = pil_g
               j = j_in
               n = n_w
            end if
         else if ( i_in == pil_g+1 ) then
            if ( j_in == -1 ) then ! ess 
               i = pil_g
               j = 2
               n = n_e
            else if ( j_in == 0 ) then ! es
               i = pil_g
               j = 1
               n = n_e
            else if ( j_in == pil_g+1 ) then ! en
               i = 1
               j = 1
               n = n_e
            else if ( j_in == pil_g+2 ) then ! enn
               i = 1
               j = 2
               n = n_e
            else
               i = pil_g + 1 - j_in
               j = 1
               n = n_e
            end if
         else if ( i_in == pil_g+2 ) then
            if ( j_in == 0 ) then ! ees
               i = pil_g - 1
               j = 1
               n = n_e
            else if ( j_in == pil_g+1 ) then ! een
               i = 2
               j = 1
               n = n_e
            else
               i = pil_g + 1 - j_in
               j = 2
               n = n_e
            end if
         else
            if ( j_in == -1 ) then 
               i = pil_g - 1
               j = pil_g + 1 - i_in
               n = n_s
            else if ( j_in == 0 ) then
               i = pil_g
               j = pil_g + 1 - i_in
               n = n_s
            else if ( j_in == pil_g+1 ) then
               i = i_in
               j = 1
               n = n_n
            else if ( j_in == pil_g+2 ) then
               i = i_in
               j = 2
               n = n_n
            else ! interior
               i = i_in
               j = j_in
               n = n_in
            end if
         end if   
      else
         n_w = mod(n_in+4, 6)
         n_e = mod(n_in+1, 6)
         n_n = mod(n_in+2, 6)
         n_s = mod(n_in+5, 6)          
         if ( i_in == -1 ) then
            if ( j_in == 0 ) then ! wws
               i = pil_g - 1
               j = pil_g
               n = n_w
            else if ( j_in == pil_g+1 ) then ! wwn
               i = 2
               j = pil_g
               n = n_w
            else
               i = pil_g + 1 - j_in
               j = pil_g - 1
               n = n_w
            end if
         else if ( i_in == 0 ) then
            if ( j_in == -1 ) then ! wss
               i = 2
               j = pil_g
               n = n_s
            else if ( j_in == 0 ) then ! ws
               i = pil_g
               j = pil_g
               n = n_w
            else if ( j_in == pil_g+1 ) then ! wn
               i = 1
               j = pil_g
               n = n_w
            else if ( j_in == pil_g+2 ) then ! wnn
               i = 1
               j = pil_g - 1
               n = n_w
            else
               i = pil_g + 1 - j_in
               j = pil_g
               n = n_w
            end if
         else if ( i_in == pil_g+1 ) then
            if ( j_in == -1 ) then ! ess
               i = 2
               j = 1
               n = n_e
            else if ( j_in == 0 ) then ! es
               i = 1
               j = 1
               n = n_e
            else if ( j_in == pil_g+1 ) then ! en
               i = 1
               j = pil_g
               n = n_e
            else if ( j_in == pil_g+2 ) then ! enn
               i = 2
               j = pil_g
               n = n_e
            else
               i = 1
               j = j_in
               n = n_e
            end if
         else if ( i_in == pil_g+2 ) then
            if ( j_in == 0 ) then ! ees
               i = 1
               j = 2
               n = n_e
            else if ( j_in == pil_g+1 ) then ! een
               i = 1
               j = pil_g - 1
               n = n_e
            else
               i = 2
               j = j_in
               n = n_e
            end if
          else
            if ( j_in == -1 ) then 
               i = i_in
               j = pil_g - 1
               n = n_s
            else if ( j_in == 0 ) then
               i = i_in
               j = pil_g
               n = n_s
            else if ( j_in == pil_g+1 ) then
               i = 1
               j = pil_g + 1 - i_in
               n = n_n
            else if ( j_in == pil_g+2 ) then
               i = 2
               j = pil_g + 1 - i_in
               n = n_n
            else ! interior
               i = i_in
               j = j_in
               n = n_in
            end if
         end if   
      end if
  
   end subroutine procarray_ijn
   
   pure function procarray(i_in,j_in,n_in) result(proc_out)
      ! returns the process id the holds the file corresponding to the input grid point
      integer, intent(in) :: i_in, j_in, n_in
      integer :: proc_out
      integer :: ip, ca, cb, cc, i, j, n
      
      proc_out = -1 ! missing
      i = i_in
      j = j_in
      n = n_in
      call procarray_ijn(i,j,n)
      
      ! determine processor that owns the grid point
      do ip = 0,fnproc-1
         cc = n + pnoff(ip) - 1
         if ( cc>=0 .and. cc<=pnpan-1 ) then
            cb = pjoff(ip,n)
            if ( j>=cb .and. j<=cb+pjpan ) then
               ca = pioff(ip,n)
               if ( i>=ca .and. i<=ca+pipan ) then
                  proc_out = ip 
                  exit
               end if
            end if
         end if   
      end do   
  
   end function procarray
   
   pure subroutine file_ijnpg2ijnp(iloc,jloc,nloc,floc,iproc,ik)
      ! converts file global index to local index
      integer, intent(inout) :: iloc, jloc, nloc
      integer, intent(in) :: floc, iproc,ik
      integer :: i, j, n
      integer :: ip, ca, cb
  
      i = iloc ! global i
      j = jloc ! global j
      n = nloc ! global panel
      
      ! fix up out-of-bounds indices
      if ( iloc == 0 ) then
         if ( jloc == 0 ) then
            if ( mod(n,2) == 0 ) then
               i = ik
               j = 1
               n = mod(nloc+5,6) ! n_w
            else
               i = ik
               j = ik
               n = mod(nloc+4,6) ! n_w
            end if    
         else if ( jloc == ik+1 ) then
            if ( mod(n,2) == 0 ) then
               i = ik
               j = ik
               n = mod(nloc+5,6) ! n_w
            else
               i = 1
               j = ik
               n = mod(nloc+4,6) ! n_w
            end if 
         else    
            if ( mod(nloc,2) == 0 ) then          
               i = ik
               j = jloc
               n = mod(nloc+5,6) ! n_w
            else
               i = ik + 1 - jloc
               j = ik
               n = mod(nloc+4,6) ! n_w
            end if
         end if   
      else if ( iloc == ik+1 ) then
         if ( jloc == 0 ) then
            if ( mod(n,2) == 0 ) then
               i = ik
               j = 1
               n = mod(nloc+2,6) ! n_e
            else
               i = 1
               j = 1
               n = mod(nloc+1,6) ! n_e
            end if 
         else if ( jloc == ik+1 ) then
            if ( mod(n,2) == 0 ) then
               i = 1
               j = 1
               n = mod(nloc+2,6) ! n_e
            else
               i = 1
               j = ik
               n = mod(nloc+1,6) ! n_e
            end if 
         else    
            if ( mod(nloc,2) == 0 ) then          
               i = ik + 1 - jloc
               j = 1
               n = mod(nloc+2,6) ! n_e
            else
               i = 1
               j = jloc
               n = mod(nloc+1,6) ! n_e
            end if
         end if   
      else
         if ( jloc == 0 ) then
            if ( mod(nloc,2) == 0 ) then              
               i = ik
               j = ik + 1 - iloc
               n = mod(nloc+4,6) ! n_s
            else
               i = iloc
               j = ik
               n = mod(nloc+5,6) ! n_s
            end if
         else if ( jloc == ik+1 ) then
            if ( mod(nloc,2) == 0 ) then              
               i = iloc
               j = 1
               n = mod(nloc+1,6) ! n_n
            else
               i = 1
               j = ik + 1 - iloc
               n = mod(nloc+2,6) ! n_n
            end if
         end if   
      end if
         
      ip = (floc-1)*fnresid + iproc
      ca = pioff(ip,n)
      cb = pjoff(ip,n)
      iloc = i - ca        ! local i
      jloc = j - cb        ! local j
      nloc = n + pnoff(ip) ! local n
   
   end subroutine file_ijnpg2ijnp

#ifdef usempi3   
   subroutine ccmpi_allocshdata2r(pdata,sshape,win,comm_in,myid_in)
      use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer

      integer, intent(out) :: win
      integer, intent(in), optional :: comm_in, myid_in
      integer :: lmyid
      integer, dimension(1), intent(in) :: sshape
      integer(kind=MPI_ADDRESS_KIND) :: qsize, lsize
      integer(kind=4) :: disp_unit, ierr, tsize
      integer(kind=4) :: lcomm, lwin
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      real, pointer, dimension(:) :: pdata 
      type(c_ptr) :: baseptr

!     allocted a single shared memory region on each node
      if ( present(comm_in) ) then
         lcomm = comm_in
      else
         lcomm = comm_node
      end if
      if ( present(myid_in) ) then
         lmyid = myid_in
      else
         lmyid = node_myid
      end if
      call MPI_Type_size( ltype, tsize, ierr )
      if ( lmyid == 0 ) then
         lsize = sshape(1)*tsize
      else
         lsize = 0_4
      end if
      call MPI_Win_allocate_shared( lsize, 1_4, MPI_INFO_NULL, lcomm, baseptr, lwin, ierr )
      if ( lmyid /= 0 ) then
         call MPI_Win_shared_query( lwin, 0_4, qsize, disp_unit, baseptr, ierr )
      end if
      call c_f_pointer( baseptr, pdata, sshape )
      win = lwin

   end subroutine ccmpi_allocshdata2r 

   subroutine ccmpi_allocshdata3r(pdata,sshape,win,comm_in,myid_in)
      use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer

      integer, intent(out) :: win
      integer, intent(in), optional :: comm_in, myid_in
      integer :: lmyid
      integer, dimension(2), intent(in) :: sshape
      integer(kind=MPI_ADDRESS_KIND) :: qsize, lsize
      integer(kind=4) :: disp_unit, ierr, tsize
      integer(kind=4) :: lcomm, lwin
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      real, pointer, dimension(:,:) :: pdata 
      type(c_ptr) :: baseptr

!     allocted a single shared memory region on each node
      if ( present(comm_in) ) then
         lcomm = comm_in
      else
         lcomm = comm_node
      end if
      if ( present(myid_in) ) then
         lmyid = myid_in
      else
         lmyid = node_myid
      end if
      call MPI_Type_size( ltype, tsize, ierr )
      if ( lmyid == 0 ) then
         lsize = sshape(1)*sshape(2)*tsize
      else
         lsize = 0_4
      end if
      call MPI_Win_allocate_shared( lsize, 1_4, MPI_INFO_NULL, lcomm, baseptr, lwin, ierr )
      if ( lmyid /= 0 ) then
         call MPI_Win_shared_query( lwin, 0_4, qsize, disp_unit, baseptr, ierr )
      end if
      call c_f_pointer( baseptr, pdata, sshape )
      win = lwin

   end subroutine ccmpi_allocshdata3r 
   
   subroutine ccmpi_allocshdata4r(pdata,sshape,win,comm_in,myid_in)
      use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer

      integer, dimension(3), intent(in) :: sshape
      integer, intent(out) :: win
      integer, intent(in), optional :: comm_in, myid_in
      integer :: lmyid
      integer(kind=MPI_ADDRESS_KIND) :: qsize, lsize
      integer(kind=4) :: disp_unit, ierr, tsize
      integer(kind=4) :: lcomm, lwin
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      real, pointer, dimension(:,:,:) :: pdata 
      type(c_ptr) :: baseptr

!     allocted a single shared memory region on each node
      if ( present(comm_in) ) then
         lcomm = comm_in
      else
         lcomm = comm_node
      end if
      if ( present(myid_in) ) then
         lmyid = myid_in
      else
         lmyid = node_myid
      end if
      call MPI_Type_size(ltype, tsize, ierr)
      if ( lmyid == 0 ) then
         lsize = sshape(1)*sshape(2)*sshape(3)*tsize
      else
         lsize = 0_4
      end if
      call MPI_Win_allocate_shared(lsize, 1_4, MPI_INFO_NULL, lcomm, baseptr, lwin, ierr)
      if ( lmyid /= 0 ) then
         call MPI_Win_shared_query(lwin, 0_4, qsize, disp_unit, baseptr, ierr)
      end if
      call c_f_pointer(baseptr, pdata, sshape)
      win = lwin

   end subroutine ccmpi_allocshdata4r 

   subroutine ccmpi_allocshdata5r(pdata,sshape,win,comm_in,myid_in)
      use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer

      integer, intent(out) :: win
      integer, intent(in), optional :: comm_in, myid_in
      integer :: lmyid
      integer, dimension(4), intent(in) :: sshape
      integer(kind=MPI_ADDRESS_KIND) :: qsize, lsize
      integer(kind=4) :: disp_unit, ierr, tsize
      integer(kind=4) :: lcomm, lwin
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      real, pointer, dimension(:,:,:,:) :: pdata 
      type(c_ptr) :: baseptr

!     allocted a single shared memory region on each node
      if ( present(comm_in) ) then
         lcomm = comm_in
      else
         lcomm = comm_node
      end if
      if ( present(myid_in) ) then
         lmyid = myid_in
      else
         lmyid = node_myid
      end if
      call MPI_Type_size( ltype, tsize, ierr )
      if ( lmyid == 0 ) then
         lsize = sshape(1)*sshape(2)*sshape(3)*sshape(4)*tsize
      else
         lsize = 0_4
      end if
      call MPI_Win_allocate_shared( lsize, 1_4, MPI_INFO_NULL, lcomm, baseptr, lwin, ierr )
      if ( lmyid /= 0 ) then
         call MPI_Win_shared_query( lwin, 0_4, qsize, disp_unit, baseptr, ierr )
      end if
      call c_f_pointer( baseptr, pdata, sshape )
      win = lwin

   end subroutine ccmpi_allocshdata5r 

   subroutine ccmpi_allocshdata2i(pdata,sshape,win,comm_in,myid_in)
      use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer

      integer, pointer, dimension(:) :: pdata 
      integer, intent(out) :: win
      integer, intent(in), optional :: comm_in, myid_in
      integer :: lmyid
      integer, dimension(1), intent(in) :: sshape
      integer(kind=MPI_ADDRESS_KIND) :: qsize, lsize
      integer(kind=4) :: disp_unit, ierr, tsize
      integer(kind=4) :: lcomm, lwin
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_INTEGER8
#else
      integer(kind=4), parameter :: ltype = MPI_INTEGER
#endif
      type(c_ptr) :: baseptr

!     allocted a single shared memory region on each node
      if ( present(comm_in) ) then
         lcomm = comm_in
      else
         lcomm = comm_node
      end if
      if ( present(myid_in) ) then
         lmyid = myid_in
      else
         lmyid = node_myid
      end if
      call MPI_Type_size( ltype, tsize, ierr )
      if ( lmyid == 0 ) then
         lsize = sshape(1)*tsize
      else
         lsize = 0_4
      end if
      call MPI_Win_allocate_shared( lsize, 1_4, MPI_INFO_NULL, lcomm, baseptr, lwin, ierr )
      if ( lmyid /= 0 ) then
         call MPI_Win_shared_query( lwin, 0_4, qsize, disp_unit, baseptr, ierr )
      end if
      call c_f_pointer( baseptr, pdata, sshape )
      win = lwin

   end subroutine ccmpi_allocshdata2i

   subroutine ccmpi_allocshdata3i(pdata,sshape,win,comm_in,myid_in)
      use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer

      integer, pointer, dimension(:,:) :: pdata 
      integer, intent(out) :: win
      integer, intent(in), optional :: comm_in, myid_in
      integer :: lmyid
      integer, dimension(2), intent(in) :: sshape
      integer(kind=MPI_ADDRESS_KIND) :: qsize, lsize
      integer(kind=4) :: disp_unit, ierr, tsize
      integer(kind=4) :: lcomm, lwin
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_INTEGER8
#else
      integer(kind=4), parameter :: ltype = MPI_INTEGER
#endif
      type(c_ptr) :: baseptr

!     allocted a single shared memory region on each node
      if ( present(comm_in) ) then
         lcomm = comm_in
      else
         lcomm = comm_node
      end if
      if ( present(myid_in) ) then
         lmyid = myid_in
      else
         lmyid = node_myid
      end if
      call MPI_Type_size( ltype, tsize, ierr )
      if ( lmyid == 0 ) then
         lsize = sshape(1)*sshape(2)*tsize
      else
         lsize = 0_4
      end if
      call MPI_Win_allocate_shared( lsize, 1_4, MPI_INFO_NULL, lcomm, baseptr, lwin, ierr )
      if ( lmyid /= 0 ) then
         call MPI_Win_shared_query( lwin, 0_4, qsize, disp_unit, baseptr, ierr )
      end if
      call c_f_pointer( baseptr, pdata, sshape )
      win = lwin

   end subroutine ccmpi_allocshdata3i
   
   subroutine ccmpi_allocshdata5i(pdata,sshape,win,comm_in,myid_in)
      use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer

      integer, pointer, dimension(:,:,:,:) :: pdata 
      integer, intent(out) :: win
      integer, intent(in), optional :: comm_in, myid_in
      integer :: lmyid
      integer, dimension(4), intent(in) :: sshape
      integer(kind=MPI_ADDRESS_KIND) :: qsize, lsize
      integer(kind=4) :: disp_unit, ierr, tsize
      integer(kind=4) :: lcomm, lwin
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_INTEGER8
#else
      integer(kind=4), parameter :: ltype = MPI_INTEGER
#endif
      type(c_ptr) :: baseptr

!     allocted a single shared memory region on each node
      if ( present(comm_in) ) then
         lcomm = comm_in
      else
         lcomm = comm_node
      end if
      if ( present(myid_in) ) then
         lmyid = myid_in
      else
         lmyid = node_myid
      end if
      call MPI_Type_size( ltype, tsize, ierr )
      if ( lmyid == 0 ) then
         lsize = sshape(1)*sshape(2)*sshape(3)*sshape(4)*tsize
      else
         lsize = 0_4
      end if
      call MPI_Win_allocate_shared( lsize, 1_4, MPI_INFO_NULL, lcomm, baseptr, lwin, ierr )
      if ( lmyid /= 0 ) then
         call MPI_Win_shared_query( lwin, 0_4, qsize, disp_unit, baseptr, ierr )
      end if
      call c_f_pointer( baseptr, pdata, sshape )
      win = lwin

   end subroutine ccmpi_allocshdata5i
   
   subroutine ccmpi_allocshdata2_r8(pdata,sshape,win,comm_in,myid_in)
      use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer

      integer, intent(out) :: win
      integer, intent(in), optional :: comm_in, myid_in
      integer :: lmyid
      integer, dimension(1), intent(in) :: sshape
      integer(kind=MPI_ADDRESS_KIND) :: qsize, lsize
      integer(kind=4) :: disp_unit, ierr, tsize
      integer(kind=4) :: lcomm, lwin
      real(kind=8), pointer, dimension(:) :: pdata 
      type(c_ptr) :: baseptr

!     allocted a single shared memory region on each node
      if ( present(comm_in) ) then
         lcomm = comm_in
      else
         lcomm = comm_node
      end if
      if ( present(myid_in) ) then
         lmyid = myid_in
      else
         lmyid = node_myid
      end if
      call MPI_Type_size( MPI_DOUBLE_PRECISION, tsize, ierr )
      if ( lmyid == 0 ) then
         lsize = sshape(1)*tsize
      else
         lsize = 0_4
      end if
      call MPI_Win_allocate_shared( lsize, 1_4, MPI_INFO_NULL, lcomm, baseptr, lwin, ierr )
      if ( lmyid /= 0 ) then
         call MPI_Win_shared_query( lwin, 0_4, qsize, disp_unit, baseptr, ierr )
      end if
      call c_f_pointer( baseptr, pdata, sshape )
      win = lwin

   end subroutine ccmpi_allocshdata2_r8
   
   subroutine ccmpi_allocshdata3_r8(pdata,sshape,win,comm_in,myid_in)
      use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer

      integer, intent(out) :: win
      integer, intent(in), optional :: comm_in, myid_in
      integer :: lmyid
      integer, dimension(2), intent(in) :: sshape
      integer(kind=MPI_ADDRESS_KIND) :: qsize, lsize
      integer(kind=4) :: disp_unit, ierr, tsize
      integer(kind=4) :: lcomm, lwin
      real(kind=8), pointer, dimension(:,:) :: pdata 
      type(c_ptr) :: baseptr

!     allocted a single shared memory region on each node
      if ( present(comm_in) ) then
         lcomm = comm_in
      else
         lcomm = comm_node
      end if
      if ( present(myid_in) ) then
         lmyid = myid_in
      else
         lmyid = node_myid
      end if
      call MPI_Type_size( MPI_DOUBLE_PRECISION, tsize, ierr )
      if ( lmyid == 0 ) then
         lsize = product(sshape)*tsize
      else
         lsize = 0_4
      end if
      call MPI_Win_allocate_shared( lsize, 1_4, MPI_INFO_NULL, lcomm, baseptr, lwin, ierr )
      if ( lmyid /= 0 ) then
         call MPI_Win_shared_query( lwin, 0_4, qsize, disp_unit, baseptr, ierr )
      end if
      call c_f_pointer( baseptr, pdata, sshape )
      win = lwin

   end subroutine ccmpi_allocshdata3_r8 

   subroutine ccmpi_allocshdata4_r8(pdata,sshape,win,comm_in,myid_in)
      use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer

      integer, intent(out) :: win
      integer, intent(in), optional :: comm_in, myid_in
      integer :: lmyid
      integer, dimension(3), intent(in) :: sshape
      integer(kind=MPI_ADDRESS_KIND) :: qsize, lsize
      integer(kind=4) :: disp_unit, ierr, tsize
      integer(kind=4) :: lcomm, lwin
      real(kind=8), pointer, dimension(:,:,:) :: pdata 
      type(c_ptr) :: baseptr

!     allocted a single shared memory region on each node
      if ( present(comm_in) ) then
         lcomm = comm_in
      else
         lcomm = comm_node
      end if
      if ( present(myid_in) ) then
         lmyid = myid_in
      else
         lmyid = node_myid
      end if
      call MPI_Type_size( MPI_DOUBLE_PRECISION, tsize, ierr )
      if ( lmyid == 0 ) then
         lsize = product(sshape)*tsize
      else
         lsize = 0_4
      end if
      call MPI_Win_allocate_shared( lsize, 1_4, MPI_INFO_NULL, lcomm, baseptr, lwin, ierr )
      if ( lmyid /= 0 ) then
         call MPI_Win_shared_query( lwin, 0_4, qsize, disp_unit, baseptr, ierr )
      end if
      call c_f_pointer( baseptr, pdata, sshape )
      win = lwin

   end subroutine ccmpi_allocshdata4_r8 

   subroutine ccmpi_shepoch(win,assert)
   
       integer, intent(in) :: win
       integer(kind=4) :: lwin, lerr, lassert
       character(len=*), intent(in), optional :: assert
       
       lassert = 0_4
       if ( present(assert) ) then
          select case(assert)
             case('noprecede')
                lassert = MPI_MODE_NOPRECEDE 
             case('nosucceed')
                lassert = MPI_MODE_NOSUCCEED
          end select
       end if
       
       lwin = win
       call MPI_Win_fence( lassert, lwin, lerr )

   end subroutine ccmpi_shepoch
   
   subroutine ccmpi_freeshdata(win)
   
      integer, intent(in) :: win
      integer(kind=4) :: lwin, lerr
      
      lwin = win
      call MPI_win_free( lwin, lerr ) 
      
   end subroutine ccmpi_freeshdata
#endif
#endif

#ifdef scm
   use newmpar_m
   implicit none
   private

   integer, save, public :: comm_world                                     ! global communication group
   integer, save, public :: myid                                           ! processor rank for comm_world
   
   integer, parameter, public :: nagg = 3                                  ! maximum number of levels to aggregate for message
                                                                           ! passing
   logical, save, public :: uniform_decomp                                 ! uniform decomposition flag
   integer, save, public :: vnode_nproc, comm_vnode, vnode_myid
   integer, save, public :: vnode_vleaderid
   integer, save, public :: vleader_nproc, comm_vleader, vleader_myid
      
   logical, public, save :: mydiag ! True if diagnostic point id, jd is in my region
  
   integer, save, public :: pipan, pjpan, pnpan                            ! decomposition parameters file window
   integer, save, public :: pil_g, pjl_g, pka_g, pko_g                     ! decomposition parameters file window
   integer, save, public :: fnproc, fnresid, fncount, mynproc              ! number and decomposition of input files
   integer, allocatable, dimension(:), save, public :: pnoff               ! file window panel offset
   integer, allocatable, dimension(:,:), save, public :: pioff, pjoff      ! file window coordinate offset
   integer(kind=4), allocatable, dimension(:), save, public ::              &
      filemap_recv, filemap_rmod                                           ! file map received for onthefly
   integer(kind=4), allocatable, dimension(:), save, public ::              &
      filemap_send, filemap_smod                                           ! file map sent for onthefly
   integer, allocatable, dimension(:), save, public :: filemap_facecomm    ! communicator for face bcast
   integer, allocatable, dimension(:), save, public :: filemap_rinv        ! inverse map for filemap_recv/rmod
   integer, save, public :: fnproc_bcast_max = 6                           ! max number of input files for bcast method

   integer, public, save :: histrd3_begin, histrd3_end
   integer, public, save :: histrd4_begin, histrd4_end
   integer, public, save :: histrd5_begin, histrd5_end
   integer, public, save :: distribute_begin, distribute_end
   integer, public, save :: radinit_begin,radinit_end
   integer, public, save :: radmisc_begin,radmisc_end
   integer, public, save :: radsw_begin, radsw_end
   integer, public, save :: radlw_begin, radlw_end
   integer, public, save :: sfluxnet_begin, sfluxnet_end
   integer, public, save :: sfluxwater_begin, sfluxwater_end
   integer, public, save :: sfluxland_begin, sfluxland_end
   integer, public, save :: sfluxurban_begin, sfluxurban_end

   public :: ccmpi_bcast, ccmpi_abort, ccmpi_barrier, ccmpi_scatterx, ccmpi_gather
   public :: ccmpi_filewinfree, ccmpi_commfree, ccmpi_commsplit
   public :: ccmpi_gatherx, ccmpi_distribute, ccmpi_reduce, ccmpi_bcastr8
   public :: ccmpi_reducer8, ccmpi_distributer8, ccmpi_gatherxr8, ccmpi_gatherr8
   public :: dix_set, face_set, uniform_set
   public :: start_log, end_log
   
   interface ccmpi_bcast
      module procedure ccmpi_bcast1i, ccmpi_bcast2i, ccmpi_bcast3i, ccmpi_bcast1r, ccmpi_bcast2r, &
                       ccmpi_bcast3r, ccmpi_bcast4r, ccmpi_bcast5r, ccmpi_bcast1s
   end interface
   interface ccmpi_scatterx
      module procedure ccmpi_scatterx2r, ccmpi_scatterx32r, ccmpi_scatterx3r
   end interface
   interface ccmpi_gather
      module procedure host_gather2, host_gather3, host_gather4
   end interface
   interface ccmpi_gatherr8
      module procedure host_gather2r8, host_gather3r8, host_gather4r8
   end interface
   interface ccmpi_gatherx
      module procedure ccmpi_gatherx2r, ccmpi_gatherx3r, ccmpi_gatherx4r
      module procedure ccmpi_gatherx23r, ccmpi_gatherx34r, ccmpi_gatherx45r
      module procedure ccmpi_gatherx3i
   end interface
   interface ccmpi_gatherxr8
     module procedure ccmpi_gatherx2rr8, ccmpi_gatherx3rr8
     module procedure ccmpi_gatherx23rr8, ccmpi_gatherx34rr8, ccmpi_gatherx45rr8
   end interface
   interface ccmpi_distribute
      module procedure host_distribute2, host_distribute2i,  &    
                       host_distribute3, host_distribute3i
      module procedure host_distribute4
   end interface
   interface ccmpi_distributer8
      module procedure host_distribute2r8, host_distribute3r8, host_distribute4r8
   end interface
   interface ccmpi_reduce
      module procedure ccmpi_reduce2i, ccmpi_reduce1r, ccmpi_reduce2r, ccmpi_reduce3r, &
                       ccmpi_reduce2c, ccmpi_reduce2l
   end interface
   interface ccmpi_reducer8
     module procedure ccmpi_reduce1rr8
   end interface
   interface ccmpi_bcastr8
      module procedure ccmpi_bcast2r8, ccmpi_bcast3r8, ccmpi_bcast4r8
   end interface   
   
   contains
    
   subroutine ccmpi_bcast1i(ldat,host,comm)
      integer, intent(in) :: host, comm
      integer, intent(inout) :: ldat
   end subroutine ccmpi_bcast1i

   subroutine ccmpi_bcast2i(ldat,host,comm)
      integer, intent(in) :: host, comm
      integer, dimension(:), intent(inout) :: ldat
   end subroutine ccmpi_bcast2i

   subroutine ccmpi_bcast3i(ldat,host,comm)
      integer, intent(in) :: host, comm
      integer, dimension(:,:), intent(inout) :: ldat
   end subroutine ccmpi_bcast3i

   subroutine ccmpi_bcast1r(ldat,host,comm)
      integer, intent(in) :: host, comm
      real, intent(inout) :: ldat
   end subroutine ccmpi_bcast1r
   
   subroutine ccmpi_bcast2r(ldat,host,comm)
      integer, intent(in) :: host, comm
      real, dimension(:), intent(inout) :: ldat
   end subroutine ccmpi_bcast2r

   subroutine ccmpi_bcast3r(ldat,host,comm)
      integer, intent(in) :: host, comm
      real, dimension(:,:), intent(inout) :: ldat
   end subroutine ccmpi_bcast3r

   subroutine ccmpi_bcast4r(ldat,host,comm)
      integer, intent(in) :: host, comm
      real, dimension(:,:,:), intent(inout) :: ldat
   end subroutine ccmpi_bcast4r

   subroutine ccmpi_bcast5r(ldat,host,comm)
      integer, intent(in) :: host, comm
      real, dimension(:,:,:,:), intent(inout) :: ldat
   end subroutine ccmpi_bcast5r

   subroutine ccmpi_bcast1s(ldat,host,comm)
      integer, intent(in) :: host, comm
      character(len=*), intent(inout) :: ldat
   end subroutine ccmpi_bcast1s
    
   subroutine ccmpi_abort(ierr)
      integer, intent(in) :: ierr
      call finishbanner
      stop -1
   end subroutine ccmpi_abort
   
   subroutine ccmpi_barrier(comm)
      integer, intent(in) :: comm
   end subroutine ccmpi_barrier

   subroutine ccmpi_scatterx2r(gdat,ldat,host,comm)
      integer, intent(in) :: host, comm
      real, dimension(:), intent(in) :: gdat
      real, dimension(:), intent(out) :: ldat
      ldat = gdat
   end subroutine ccmpi_scatterx2r

   subroutine ccmpi_scatterx32r(gdat,ldat,host,comm)
      integer, intent(in) :: host, comm
      real, dimension(:,:), intent(in) :: gdat
      real, dimension(:), intent(out) :: ldat
      ldat = reshape(gdat,(/size(ldat)/))
   end subroutine ccmpi_scatterx32r

   subroutine ccmpi_scatterx3r(gdat,ldat,host,comm)
      integer, intent(in) :: host, comm
      real, dimension(:,:), intent(in) :: gdat
      real, dimension(:,:), intent(out) :: ldat
      ldat = gdat
   end subroutine ccmpi_scatterx3r
   
   subroutine host_gather2(a,ag)
      real, dimension(ifull), intent(in) :: a
      real, dimension(ifull_g), intent(out), optional :: ag
      if ( present(ag) ) then
        ag = a
      end if  
   end subroutine host_gather2

   subroutine host_gather3(a,ag)
      real, dimension(:,:), intent(in) :: a
      real, dimension(:,:), intent(out), optional :: ag
      if ( present(ag) ) then
        ag = a
      end if  
   end subroutine host_gather3

   subroutine host_gather4(a,ag)
      real, dimension(:,:,:), intent(in) :: a
      real, dimension(:,:,:), intent(out), optional :: ag
      if ( present(ag) ) then
        ag = a
      end if  
   end subroutine host_gather4
   
   subroutine host_gather2r8(a,ag)
      real(kind=8), dimension(ifull), intent(in) :: a
      real(kind=8), dimension(ifull_g), intent(out), optional :: ag
      if ( present(ag) ) then
        ag = a
      end if  
   end subroutine host_gather2r8

   subroutine host_gather3r8(a,ag)
      real(kind=8), dimension(:,:), intent(in) :: a
      real(kind=8), dimension(:,:), intent(out), optional :: ag
      if ( present(ag) ) then
        ag = a
      end if  
   end subroutine host_gather3r8

   subroutine host_gather4r8(a,ag)
      real(kind=8), dimension(:,:,:), intent(in) :: a
      real(kind=8), dimension(:,:,:), intent(out), optional :: ag
      if ( present(ag) ) then
        ag = a
      end if  
   end subroutine host_gather4r8
   
   subroutine ccmpi_filewinfree
   end subroutine ccmpi_filewinfree
   
   subroutine ccmpi_commfree(comm)
      integer, intent(in) :: comm
   end subroutine ccmpi_commfree
   
   subroutine ccmpi_commsplit(commout,comm,colour,rank)
      integer, intent(out) :: commout
      integer, intent(in) :: comm, colour, rank
         commout = comm
   end subroutine ccmpi_commsplit
   
   subroutine dix_set(ipan_l,jpan_l,noff_l,ioff_l,joff_l,npan_l,il_gx,myid_l,nproc_l,nxproc_l,nyproc_l)
      integer, intent(in) :: myid_l, nproc_l, npan_l, il_gx
      integer, intent(out) :: ipan_l, jpan_l, noff_l, nxproc_l, nyproc_l
      integer, dimension(0:npanels), intent(out) :: ioff_l, joff_l
         ipan_l = 1
         jpan_l = 1
         noff_l = 1
         nxproc_l = 1
         nyproc_l = 1
         ioff_l(:) = 0
         joff_l(:) = 0 
   end subroutine dix_set
   
   subroutine face_set(ipan_l, jpan_l, noff_l, ioff_l, joff_l, npan_l, il_gx, myid_l, nproc_l, nxproc_l, nyproc_l)
      integer, intent(in) :: myid_l, nproc_l, npan_l, il_gx
      integer, intent(out) :: ipan_l, jpan_l, noff_l, nxproc_l, nyproc_l
      integer, dimension(0:npanels), intent(out) :: ioff_l, joff_l
         ipan_l = 1
         jpan_l = 1
         noff_l = 1
         nxproc_l = 1
         nyproc_l = 1
         ioff_l(:) = 0
         joff_l(:) = 0
   end subroutine face_set

   subroutine uniform_set(ipan_l,jpan_l,noff_l,ioff_l,joff_l,npan_l,il_gx,myid_l,nproc_l,nxproc_l,nyproc_l)
      integer, intent(in) :: myid_l, nproc_l, npan_l, il_gx
      integer, intent(out) :: ipan_l, jpan_l, noff_l, nxproc_l, nyproc_l
      integer, dimension(0:npanels), intent(out) :: ioff_l, joff_l 
         ipan_l = 1
         jpan_l = 1
         noff_l = 1
         nxproc_l = 1
         nyproc_l = 1
         ioff_l(:) = 0
         joff_l(:) = 0
   end subroutine uniform_set
   
   subroutine ccmpi_gatherx2r(gdat,ldat,host,comm)
      integer, intent(in) :: host, comm
      real, dimension(:), intent(out) :: gdat
      real, dimension(:), intent(in) :: ldat
      gdat = ldat
   end subroutine ccmpi_gatherx2r
   
   subroutine ccmpi_gatherx3r(gdat,ldat,host,comm)
      integer, intent(in) :: host, comm
      real, dimension(:,:), intent(out) :: gdat
      real, dimension(:,:), intent(in) :: ldat
      gdat = ldat
   end subroutine ccmpi_gatherx3r

    subroutine ccmpi_gatherx4r(gdat,ldat,host,comm)
      integer, intent(in) :: host, comm
      real, dimension(:,:,:), intent(out) :: gdat
      real, dimension(:,:,:), intent(in) :: ldat
      gdat = ldat
   end subroutine ccmpi_gatherx4r
  
   subroutine ccmpi_gatherx23r(gdat,ldat,host,comm)
      integer, intent(in) :: host, comm
      real, dimension(:,:), intent(out) :: gdat
      real, dimension(:), intent(in) :: ldat
      gdat = reshape(ldat,(/size(gdat,1),size(gdat,2)/))
   end subroutine ccmpi_gatherx23r
   
   subroutine ccmpi_gatherx34r(gdat,ldat,host,comm)
      integer, intent(in) :: host, comm
      real, dimension(:,:,:), intent(out) :: gdat
      real, dimension(:,:), intent(in) :: ldat
      gdat = reshape(ldat,(/size(gdat,1),size(gdat,2),size(gdat,3)/))
   end subroutine ccmpi_gatherx34r

   subroutine ccmpi_gatherx45r(gdat,ldat,host,comm)
      integer, intent(in) :: host, comm
      real, dimension(:,:,:,:), intent(out) :: gdat
      real, dimension(:,:,:), intent(in) :: ldat
      gdat = reshape(ldat,(/size(gdat,1),size(gdat,2),size(gdat,3),size(gdat,4)/))
   end subroutine ccmpi_gatherx45r
   
   subroutine ccmpi_gatherx3i(gdat,ldat,host,comm)
      integer, intent(in) :: host, comm
      integer, dimension(:,:), intent(out) :: gdat
      integer, dimension(:,:), intent(in) :: ldat
      gdat = ldat
   end subroutine ccmpi_gatherx3i
 
   subroutine ccmpi_gatherx2rr8(gdat,ldat,host,comm)
      integer, intent(in) :: host, comm
      real(kind=8), dimension(:), intent(out) :: gdat
      real(kind=8), dimension(:), intent(in) :: ldat
      gdat = ldat
   end subroutine ccmpi_gatherx2rr8

   subroutine ccmpi_gatherx3rr8(gdat,ldat,host,comm)
      integer, intent(in) :: host, comm
      real(kind=8), dimension(:,:), intent(out) :: gdat
      real(kind=8), dimension(:,:), intent(in) :: ldat
      gdat = ldat
   end subroutine ccmpi_gatherx3rr8
   
   subroutine ccmpi_gatherx23rr8(gdat,ldat,host,comm)
      integer, intent(in) :: host, comm
      real(kind=8), dimension(:,:), intent(out) :: gdat
      real(kind=8), dimension(:), intent(in) :: ldat
      gdat = reshape(ldat,(/size(gdat,1),size(gdat,2)/))
   end subroutine ccmpi_gatherx23rr8
   
   subroutine ccmpi_gatherx34rr8(gdat,ldat,host,comm)
      integer, intent(in) :: host, comm
      real(kind=8), dimension(:,:,:), intent(out) :: gdat
      real(kind=8), dimension(:,:), intent(in) :: ldat
      gdat = reshape(ldat,(/size(gdat,1),size(gdat,2),size(gdat,3)/))
   end subroutine ccmpi_gatherx34rr8

   subroutine ccmpi_gatherx45rr8(gdat,ldat,host,comm)
      integer, intent(in) :: host, comm
      real(kind=8), dimension(:,:,:,:), intent(out) :: gdat
      real(kind=8), dimension(:,:,:), intent(in) :: ldat
      gdat = reshape(ldat,(/size(gdat,1),size(gdat,2),size(gdat,3),size(gdat,4)/))
   end subroutine ccmpi_gatherx45rr8
   
   subroutine host_distribute2(af,a1)
      real, dimension(ifull), intent(out) :: af
      real, dimension(ifull_g), intent(in), optional :: a1
      if ( present(a1) ) then
        af = a1
      end if
   end subroutine host_distribute2

   subroutine host_distribute2i(af,a1)
      integer, dimension(ifull), intent(out) :: af
      integer, dimension(ifull_g), intent(in), optional :: a1
      if ( present(a1) ) then
        af = a1
      end if
   end subroutine host_distribute2i

   subroutine host_distribute3(af,a1)
      real, dimension(:,:), intent(out) :: af
      real, dimension(:,:), intent(in), optional :: a1
      if ( present(a1) ) then
        af = a1
      end if
   end subroutine host_distribute3

   subroutine host_distribute4(af,a1)
      real, dimension(:,:,:), intent(out) :: af
      real, dimension(:,:,:), intent(in), optional :: a1
      if ( present(a1) ) then
        af = a1
      end if
   end subroutine host_distribute4
   
   subroutine host_distribute3i(af,a1)
      integer, dimension(:,:), intent(out) :: af
      integer, dimension(:,:), intent(in), optional :: a1
      if ( present(a1) ) then
        af = a1
      end if
   end subroutine host_distribute3i

   subroutine host_distribute2r8(af,a1)
      real(kind=8), dimension(ifull), intent(out) :: af
      real(kind=8), dimension(ifull_g), intent(in), optional :: a1
      if ( present(a1) ) then
        af = a1
      end if
   end subroutine host_distribute2r8
   
   subroutine host_distribute3r8(af,a1)
      real(kind=8), dimension(:,:), intent(out) :: af
      real(kind=8), dimension(:,:), intent(in), optional :: a1
      if ( present(a1) ) then
        af = a1
      end if
   end subroutine host_distribute3r8
 
   subroutine host_distribute4r8(af,a1)
      real(kind=8), dimension(:,:,:), intent(out) :: af
      real(kind=8), dimension(:,:,:), intent(in), optional :: a1
      if ( present(a1) ) then
        af = a1
      end if
   end subroutine host_distribute4r8
   
   subroutine ccmpi_reduce2i(ldat,gdat,op,host,comm)
      integer, intent(in) :: host,comm
      integer, dimension(:), intent(in) :: ldat
      integer, dimension(:), intent(out) :: gdat
      character(len=*), intent(in) :: op
      gdat = ldat
   end subroutine ccmpi_reduce2i

   subroutine ccmpi_reduce1r(ldat,gdat,op,host,comm)
      integer, intent(in) :: host, comm
      real, intent(in) :: ldat
      real, intent(out) :: gdat
      character(len=*), intent(in) :: op
      gdat = ldat
   end subroutine ccmpi_reduce1r
   
   subroutine ccmpi_reduce2r(ldat,gdat,op,host,comm)
      integer, intent(in) :: host, comm
      real, dimension(:), intent(in) :: ldat
      real, dimension(:), intent(out) :: gdat
      character(len=*), intent(in) :: op
      gdat = ldat
   end subroutine ccmpi_reduce2r

   subroutine ccmpi_reduce3r(ldat,gdat,op,host,comm)
      integer, intent(in) :: host, comm
      real, dimension(:,:), intent(in) :: ldat
      real, dimension(:,:), intent(out) :: gdat
      character(len=*), intent(in) :: op
      gdat = ldat
   end subroutine ccmpi_reduce3r

   subroutine ccmpi_reduce2c(ldat,gdat,op,host,comm)
      integer, intent(in) :: host,comm
      complex, dimension(:), intent(in) :: ldat
      complex, dimension(:), intent(out) :: gdat
      character(len=*), intent(in) :: op
      gdat = ldat
   end subroutine ccmpi_reduce2c

   subroutine ccmpi_reduce2l(ldat,gdat,op,host,comm)
      integer, intent(in) :: host,comm
      logical, dimension(:), intent(in) :: ldat
      logical, dimension(:), intent(out) :: gdat
      character(len=*), intent(in) :: op
      gdat = ldat
   end subroutine ccmpi_reduce2l

   subroutine ccmpi_reduce1rr8(ldat,gdat,op,host,comm)
      integer, intent(in) :: host,comm
      real(kind=8), intent(in) :: ldat
      real(kind=8), intent(out) :: gdat
      character(len=*), intent(in) :: op
      gdat = ldat
   end subroutine ccmpi_reduce1rr8 
   
   subroutine ccmpi_bcast2r8(ldat,host,comm)
      integer, intent(in) :: host, comm
      real(kind=8), dimension(:), intent(inout) :: ldat
   end subroutine ccmpi_bcast2r8
   
   subroutine ccmpi_bcast3r8(ldat,host,comm)
      integer, intent(in) :: host,comm
      real(kind=8), dimension(:,:), intent(inout) :: ldat
   end subroutine ccmpi_bcast3r8   
   
   subroutine ccmpi_bcast4r8(ldat,host,comm)
      integer, intent(in) :: host,comm
      real(kind=8), dimension(:,:,:), intent(inout) :: ldat
   end subroutine ccmpi_bcast4r8
   
   subroutine start_log ( event )
      integer, intent(in) :: event
   end subroutine start_log
   
   subroutine end_log( event )
      integer, intent(in) :: event
   end subroutine end_log
#endif

end module cc_mpi

