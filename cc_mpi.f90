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
    
module cc_mpi

! This module manages all MPI communications between processors.  The system was originally developed
! by Martin Dix and subsequently modified by Marcus Thatcher.  Thanks to Aaron McDonough for developing
! the Vampir trace routines and upgrading the timer calls.  Thanks to Paul Ryan for the design of the
! shared memory arrays.

#ifdef usempi_mod
   use mpi
#else
   use mpif_m
#endif
   implicit none
   private
   include 'newmpar.h'

   integer, save, public :: comm_world                                     ! global communication group
   integer, save, public :: myid                                           ! processor rank for comm_world
   integer, save, public :: ipan, jpan                                     ! grid size on processor
   integer, save, public :: ioff, joff, noff                               ! offset of processor grid relative to global grid
   integer, save, public :: nxproc, nyproc                                 ! number of processors in the x and y directions
   integer, save, public :: nagg                                           ! maximum number of levels to aggregate for message
                                                                           ! passing

#ifdef usempi3
   integer, save, public :: comm_node                                      ! node communication group
   integer, save, public :: comm_nodecaptian                               ! node captian communication group
   integer, save, public :: node_myid                                      ! processor rank for comm_node
   integer, save, public :: node_nproc                                     ! number of processors on a node
   integer, save, public :: nodecaptian_myid                               ! node rank (with captian)
   integer, save, public :: nodecaptian_nproc                              ! number of nodes (with captians)
#endif
   
   integer, save, public :: comm_proc, comm_rows, comm_cols                ! comm groups for scale-selective filter
   integer, save, public :: hproc, mproc, npta, pprocn, pprocx             ! decomposition parameters for scale-selective filter

   integer(kind=4), save, private :: nreq, rreq                            ! number of messages requested and to be received
   integer(kind=4), allocatable, dimension(:), save, private :: ireq       ! requested message index
   integer, allocatable, dimension(:), save, private :: rlist              ! map of processor index from requested message index
   
   integer, allocatable, dimension(:), save, public :: neighlist           ! list of neighbour processors
   integer, allocatable, dimension(:), save, private :: neighmap           ! map of processor to neighbour index
   integer, save, public :: neighnum                                       ! number of neighbours
   
   integer(kind=4), save, private :: localwin                              ! local window handle for spectral filter
   integer(kind=4), allocatable, dimension(:), save, public :: specmap     ! gather map for spectral filter
   integer, allocatable, dimension(:), save, public :: specmapext          ! gather map for spectral filter (includes filter final
                                                                           ! pass for sparse arrays)
   real, allocatable, dimension(:,:), save, private :: specstore           ! window for gather map
   type globalpack_info
     real, allocatable, dimension(:,:,:) :: localdata
   end type globalpack_info
   ! store sparse global arrays for spectral filter
   type(globalpack_info), allocatable, dimension(:,:,:), save, private :: globalpack                                            

   integer, save, public :: pil, pjl, pnpan                                ! decomposition parameters file window
   integer, save, public :: fnproc, fnresid, fncount                       ! number and decomposition of input files
   integer, allocatable, dimension(:), save, public :: pnoff               ! file window panel offset
   integer, allocatable, dimension(:,:), save, public :: pioff, pjoff      ! file window coordinate offset
   integer(kind=4), save, private :: filewin                               ! local window handle for onthefly 
   integer(kind=4), allocatable, dimension(:), save, public :: filemap     ! file map for onthefly
   real, allocatable, dimension(:,:), save, private :: filestore           ! window for file map
   
   integer, allocatable, dimension(:), save, private :: fileneighlist      ! list of file neighbour processors
   integer, save, public :: fileneighnum                                   ! number of file neighbours
   
   public :: ccmpi_setup, ccmpi_distribute, ccmpi_gather,                   &
             ccmpi_distributer8, ccmpi_gatherall, bounds, boundsuv,         &
             deptsync, intssync_send, intssync_recv, start_log, end_log,    &
             log_on, log_off, log_setup, phys_loadbal, ccglobal_posneg,     &
             ccglobal_sum, readglobvar, writeglobvar, ccmpi_reduce,         &
             ccmpi_allreduce, ccmpi_abort, ccmpi_bcast, ccmpi_bcastr8,      &
             ccmpi_barrier, ccmpi_gatherx, ccmpi_scatterx,                  &
             ccmpi_allgatherx, ccmpi_init, ccmpi_finalize, ccmpi_commsplit, &
             ccmpi_commfree, bounds_colour_send, bounds_colour_recv,        &
             boundsuv_allvec, boundsr8
   public :: mgbounds, mgcollect, mgbcast, mgbcastxn, mgbcasta, mg_index,   &
             mg_fproc, mg_fproc_1
   public :: ind, indx, indp, indg, iq2iqg, indv_mpi, indglobal, fproc,     &
             proc_region, proc_region_face, proc_region_dix, face_set,      &
             uniform_set, dix_set
   public :: mgbndtype, dpoints_t, dindex_t, sextra_t, bnds
   public :: getglobalpack, setglobalpack, allocateglobalpack,              &
             copyglobalpack, ccmpi_gathermap
   public :: ccmpi_filewincreate, ccmpi_filewinfree, ccmpi_filewinget,      &
             ccmpi_filewinunpack, ccmpi_filebounds_setup, ccmpi_filebounds
   private :: ccmpi_distribute2, ccmpi_distribute2i, ccmpi_distribute2r8,   &
              ccmpi_distribute3, ccmpi_distribute3i, ccmpi_gather2,         &
              ccmpi_gather3, checksize, ccglobal_posneg2, ccglobal_posneg3, &
              ccglobal_posneg4, ccglobal_sum2, ccglobal_sum3
#ifdef usempi3
   public :: ccmpi_allocshdata, ccmpi_allocshdatar8
   public :: ccmpi_shepoch, ccmpi_freeshdata
   
   interface ccmpi_allocshdata
      module procedure ccmpi_allocshdata2r, ccmpi_allocshdata4r, ccmpi_allocshdata5r, &
                       ccmpi_allocshdata2i, ccmpi_allocshdata5i
   end interface
   interface ccmpi_allocshdatar8
      module procedure ccmpi_allocshdata2_r8, ccmpi_allocshdata3_r8, ccmpi_allocshdata4_r8
   end interface
#endif
   
   interface ccmpi_gather
      module procedure ccmpi_gather2, ccmpi_gather3
   end interface
   interface ccmpi_distribute
      module procedure ccmpi_distribute2, ccmpi_distribute2i,  &    
                       ccmpi_distribute3, ccmpi_distribute3i
   end interface
   interface ccmpi_distributer8
      module procedure ccmpi_distribute2r8
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
   interface ccglobal_posneg
      module procedure ccglobal_posneg2, ccglobal_posneg3, ccglobal_posneg4
   end interface
   interface ccglobal_sum
      module procedure ccglobal_sum2, ccglobal_sum3
   end interface
   interface readglobvar
      module procedure readglobvar2, readglobvar3, readglobvar2i
   end interface
   interface writeglobvar
      module procedure writeglobvar2, writeglobvar3
   end interface
   interface ccmpi_reduce
      module procedure ccmpi_reduce2i, ccmpi_reduce1r, ccmpi_reduce2r, ccmpi_reduce3r, &
                       ccmpi_reduce2c, ccmpi_reduce2l
   end interface
   interface ccmpi_allreduce
      module procedure ccmpi_allreduce1i, ccmpi_allreduce2i, ccmpi_allreduce2r, ccmpi_allreduce3r, &
                       ccmpi_allreduce2c
   end interface
   interface ccmpi_bcast
      module procedure ccmpi_bcast1i, ccmpi_bcast2i, ccmpi_bcast3i, ccmpi_bcast1r, ccmpi_bcast2r, &
                       ccmpi_bcast3r, ccmpi_bcast4r, ccmpi_bcast5r, ccmpi_bcast1s
   end interface
   interface ccmpi_bcastr8
      module procedure ccmpi_bcast2r8, ccmpi_bcast3r8, ccmpi_bcast4r8
   end interface
   interface ccmpi_gatherx
      module procedure ccmpi_gatherx2r, ccmpi_gatherx3r, ccmpi_gatherx4r
      module procedure ccmpi_gatherx23r, ccmpi_gatherx34r
      module procedure ccmpi_gatherx3i
   end interface
   interface ccmpi_scatterx
      module procedure ccmpi_scatterx2r, ccmpi_scatterx32r, ccmpi_scatterx3r
   end interface
   interface ccmpi_allgatherx
      module procedure ccmpi_allgatherx2i, ccmpi_allgatherx2r
   end interface
   interface ccmpi_gathermap
      module procedure ccmpi_gathermap2, ccmpi_gathermap3
   end interface
   interface ccmpi_filebounds
      module procedure ccmpi_filebounds2, ccmpi_filebounds3
   end interface ccmpi_filebounds
   interface mgcollect
      module procedure mgcollect1, mgcollectreduce, mgcollectxn
   end interface
   interface ccmpi_filewinget
     module procedure ccmpi_filewinget2, ccmpi_filewinget3
   end interface

   ! Do directions need to be swapped
   logical, parameter, private, dimension(0:npanels) ::    &
           swap_e = (/ .true., .false., .true., .false., .true., .false. /), &
           swap_w = (/ .false., .true., .false., .true., .false., .true. /), &
           swap_n = (/ .false., .true., .false., .true., .false., .true. /), &
           swap_s = (/ .true., .false., .true., .false., .true., .false. /)

   type bounds_info
      real, dimension(:), allocatable :: sbuf, rbuf
      real, dimension(:), allocatable :: send_neg
      real(kind=8), dimension(:), allocatable :: s8buf, r8buf
      integer, dimension(:), allocatable :: request_list
      integer, dimension(:), allocatable :: send_list
      integer, dimension(:), allocatable :: unpack_list
      integer, dimension(:), allocatable :: request_list_uv
      integer, dimension(:), allocatable :: send_list_uv
      integer, dimension(:), allocatable :: unpack_list_uv
      ! Flag for whether u and v need to be swapped
      logical, dimension(:), allocatable :: uv_swap, send_swap, uv_neg
      ! Number of points for each processor. Also double row versions.
      ! lenx is first row plux corner points.  lenh is just the ne side.
      integer :: slen, rlen, slenx, rlenx, slen2, rlen2
      integer :: slenh, rlenh
      ! Number of points for each processor. lenx is for nu, su, ev and wv
      integer :: slen_uv, rlen_uv, slen2_uv, rlen2_uv
      integer :: slenx_uv, rlenx_uv
      integer :: len, sbuflen, rbuflen
   end type bounds_info
   
   type dpoints_t
      real, dimension(:,:), allocatable :: a
      real, dimension(:), allocatable :: b
   end type dpoints_t
   type dindex_t
      integer, dimension(:,:), allocatable :: a
   end type dindex_t
   type sextra_t
      real, dimension(:), allocatable :: a
   end type sextra_t

   type boundsplit
      integer :: isubg, ievfn
      integer :: isvbg, iwufn, invbg, ieufn
      integer :: issvbg, iwwufn, innvbg, ieeufn
   end type boundsplit
   type coloursplit
      integer, dimension(3) :: ihbg, ihfn, ifbg, iffn
   end type coloursplit

   type(bounds_info), allocatable, dimension(:), save :: bnds

   ! partition boundary indices into geometric groups or colours
   type(boundsplit), allocatable, dimension(:), save, private :: rsplit
   type(boundsplit), allocatable, dimension(:), save, private :: ssplit
   type(coloursplit), allocatable, dimension(:), save, private :: rcolsp
   type(coloursplit), allocatable, dimension(:), save, private :: scolsp
   integer, dimension(:,:), allocatable, save, public :: iqx, iqn, iqe, iqw, iqs
   integer, save, public :: ifullmaxcol
#ifdef uniform_decomp
   integer, parameter, public :: maxcolour = 3
#else
   integer, parameter, public :: maxcolour = 2
#endif
   integer, public, save, dimension(maxcolour) :: ifullcol, ifullcol_border

   integer, public, save :: maxbuflen

   ! Flag whether processor region edge is a face edge.
   logical, public, save :: edge_w, edge_n, edge_s, edge_e

   ! Off processor departure points
   type(dpoints_t), allocatable, dimension(:), public, save :: dpoints
   type(dpoints_t), allocatable, dimension(:), private, save :: dbuf
   type(dindex_t), allocatable, dimension(:), public, save :: dindex
   type(sextra_t), allocatable, dimension(:), public, save :: sextra
   ! Number of points for each processor.
   integer, dimension(:), allocatable, public, save :: dslen, drlen

   logical, public, save :: mydiag ! True if diagnostic point id, jd is in my region

   ! Multi-grid arrays
   type mgtype
      integer :: ifull, iextra, ixlen, ifull_fine, ifull_coarse
      integer :: merge_len, merge_row, ipan, merge_pos, nmax
      integer :: comm_merge, neighnum, npanx
      !integer, dimension(:,:,:), allocatable :: fproc
      integer, dimension(:), allocatable :: merge_list
      integer, dimension(:), allocatable :: in, ie, is, iw, ine, inw, ise, isw
      integer, dimension(:), allocatable :: coarse_a, coarse_b, coarse_c, coarse_d
      integer, dimension(:), allocatable :: fine, fine_n, fine_e, fine_ne
      integer, dimension(:), allocatable :: neighlist
      integer, dimension(:), allocatable :: procmap
      real, dimension(:), allocatable :: zzn, zze, zzs, zzw, zz
      real, dimension(:), allocatable :: wgt_a, wgt_bc, wgt_d
   end type mgtype

   type mgbndtype
      integer :: len, rlen, rlenx, slen, slenx
      integer, dimension(:), allocatable :: send_list
      integer, dimension(:), allocatable :: unpack_list
      integer, dimension(:), allocatable :: request_list
   end type mgbndtype

   ! Multi-grid levels and buffers
   type(mgtype), dimension(:), allocatable, save, public :: mg
   type(mgbndtype), dimension(:,:), allocatable, save, public :: mg_bnds
   integer, save, public :: mg_maxlevel, mg_maxlevel_local
   integer, save, public :: mg_ifullmaxcol
   integer, dimension(:,:), allocatable, save, public :: col_iq, col_iqn, col_iqe, col_iqs, col_iqw

   ! File IO
   type filebounds_info
      integer, dimension(:,:), allocatable :: send_list
      integer, dimension(:,:), allocatable :: request_list
      integer, dimension(:,:), allocatable :: unpack_list
      integer :: len, rlen, slen
   end type filebounds_info
   
   type(filebounds_info), allocatable, dimension(:), save :: filebnds
   
   ! Timer
   integer, public, save :: bounds_begin, bounds_end
   integer, public, save :: boundsuv_begin, boundsuv_end
   integer, public, save :: ints_begin, ints_end
   integer, public, save :: nonlin_begin, nonlin_end
   integer, public, save :: helm_begin, helm_end
   integer, public, save :: adjust_begin, adjust_end
   integer, public, save :: upglobal_begin, upglobal_end
   integer, public, save :: hordifg_begin, hordifg_end
   integer, public, save :: vadv_begin, vadv_end
   integer, public, save :: depts_begin, depts_end
   integer, public, save :: deptsync_begin, deptsync_end
   integer, public, save :: intssync_begin, intssync_end
   integer, public, save :: stag_begin, stag_end
   integer, public, save :: ocnstag_begin, ocnstag_end
   integer, public, save :: toij_begin, toij_end
   integer, public, save :: physloadbal_begin, physloadbal_end
   integer, public, save :: phys_begin, phys_end
   integer, public, save :: outfile_begin, outfile_end
   integer, public, save :: onthefly_begin, onthefly_end
   integer, public, save :: otf_ints1_begin, otf_ints1_end
   integer, public, save :: otf_ints4_begin, otf_ints4_end
   integer, public, save :: histrd1_begin, histrd1_end
   integer, public, save :: histrd4_begin, histrd4_end
   integer, public, save :: indata_begin, indata_end
   integer, public, save :: nestin_begin, nestin_end
   integer, public, save :: nestpack_begin, nestpack_end
   integer, public, save :: nestcalc_begin, nestcalc_end
   integer, public, save :: nestcomm_begin, nestcomm_end
   integer, public, save :: nestunpack_begin, nestunpack_end
   integer, public, save :: gwdrag_begin, gwdrag_end
   integer, public, save :: convection_begin, convection_end
   integer, public, save :: cloud_begin, cloud_end
   integer, public, save :: radnet_begin,radnet_end
   integer, public, save :: radmisc_begin,radmisc_end
   integer, public, save :: radsw_begin, radsw_end
   integer, public, save :: radlw_begin, radlw_end   
   integer, public, save :: sfluxnet_begin, sfluxnet_end
   integer, public, save :: sfluxwater_begin, sfluxwater_end
   integer, public, save :: sfluxland_begin, sfluxland_end
   integer, public, save :: sfluxurban_begin, sfluxurban_end
   integer, public, save :: vertmix_begin, vertmix_end
   integer, public, save :: aerosol_begin, aerosol_end
   integer, public, save :: model_begin, model_end
   integer, public, save :: maincalc_begin, maincalc_end
   integer, public, save :: gatherrma_begin, gatherrma_end
   integer, public, save :: gather_begin, gather_end
   integer, public, save :: distribute_begin, distribute_end
   integer, public, save :: globsum_begin, globsum_end
   integer, public, save :: posneg_begin, posneg_end   
   integer, public, save :: precon_begin, precon_end
   integer, public, save :: waterdynamics_begin, waterdynamics_end
   integer, public, save :: watermisc_begin, watermisc_end
   integer, public, save :: waterdeps_begin, waterdeps_end
   integer, public, save :: watereos_begin, watereos_end
   integer, public, save :: waterhadv_begin, waterhadv_end
   integer, public, save :: watervadv_begin, watervadv_end
   integer, public, save :: waterhelm_begin, waterhelm_end
   integer, public, save :: wateriadv_begin, wateriadv_end
   integer, public, save :: waterdiff_begin, waterdiff_end
   integer, public, save :: river_begin, river_end
   integer, public, save :: bcast_begin, bcast_end   
   integer, public, save :: allgatherx_begin, allgatherx_end
   integer, public, save :: gatherx_begin, gatherx_end
   integer, public, save :: scatterx_begin, scatterx_end
   integer, public, save :: reduce_begin, reduce_end
   integer, public, save :: mpiwait_begin, mpiwait_end
   integer, public, save :: mpiwaituv_begin, mpiwaituv_end
   integer, public, save :: mpiwaituvtile_begin, mpiwaituvtile_end
   integer, public, save :: mpiwaitdep_begin, mpiwaitdep_end
   integer, public, save :: mpiwaitmg_begin, mpiwaitmg_end
   integer, public, save :: mgbounds_begin, mgbounds_end
   integer, public, save :: mgcollect_begin, mgcollect_end
   integer, public, save :: mgbcast_begin, mgbcast_end
   integer, public, save :: mgsetup_begin, mgsetup_end
   integer, public, save :: mgdecomp_begin, mgdecomp_end
   integer, public, save :: mgfine_begin, mgfine_end
   integer, public, save :: mgup_begin, mgup_end
   integer, public, save :: mgcoarse_begin, mgcoarse_end
   integer, public, save :: mgdown_begin, mgdown_end
   integer, public, save :: mgmlosetup_begin, mgmlosetup_end
   integer, public, save :: mgmlodecomp_begin, mgmlodecomp_end
   integer, public, save :: mgmlofine_begin, mgmlofine_end
   integer, public, save :: mgmloup_begin, mgmloup_end
   integer, public, save :: mgmlocoarse_begin, mgmlocoarse_end
   integer, public, save :: mgmlodown_begin, mgmlodown_end
   integer, parameter :: nevents = 85
#ifdef simple_timer
   public :: simple_timer_finalize
   real(kind=8), dimension(nevents), save :: tot_time = 0., start_time
#endif
   character(len=15), dimension(nevents), save :: event_name

#ifdef vampir
#include "vt_user.inc"
#endif


contains

   subroutine ccmpi_setup(kx)
      use indices_m
      use latlong_m
      use map_m
      use sumdd_m
      use vecsuv_m
      use xyzinfo_m
      integer, intent(in) :: kx
      integer iproc, dproc, iq, iqg, i, j, n
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      integer(kind=4) ierr
      integer(kind=4) colour, rank, lcommin, lcommout
      integer(kind=4) asize
      integer(kind=MPI_ADDRESS_KIND) wsize
      integer, dimension(ifull) :: colourmask
      logical(kind=4) :: ltrue

      nreq = 0

      allocate( bnds(0:nproc-1) )
      allocate( rsplit(0:nproc-1), ssplit(0:nproc-1) )
      allocate( rcolsp(0:nproc-1), scolsp(0:nproc-1) )
      
#ifdef uniform_decomp
      call proc_setup_uniform
      ! Faces may not line up properly so need extra factor here
      maxbuflen = (max(ipan,jpan)+4)*3*max(nagg*kl,3*ol)*8*2  !*3 for extra vector row (e.g., inu,isu,iev,iwv)
#else
      call proc_setup
      if ( nproc<npanels+1 ) then
         ! This is the maximum size, each face has 4 edges
         maxbuflen = npan*4*(il_g+4)*3*max(nagg*kl,3*ol)    !*3 for extra vector row (e.g., inu,isu,iev,iwv)
      else
         maxbuflen = (max(ipan,jpan)+4)*3*max(nagg*kl,3*ol) !*3 for extra vector row (e.g., inu,isu,iev,iwv)
      end if
#endif


      if ( myid == 0 ) then
         call ccmpi_distribute(wts,wts_g)
         call ccmpi_distribute(em,em_g)
         call ccmpi_distribute(emu,emu_g)
         call ccmpi_distribute(emv,emv_g)
         call ccmpi_distribute(ax,ax_g)
         call ccmpi_distribute(bx,bx_g)
         call ccmpi_distribute(ay,ay_g)
         call ccmpi_distribute(by,by_g)
         call ccmpi_distribute(az,az_g)
         call ccmpi_distribute(bz,bz_g)
         call ccmpi_distribute(f,f_g)
         call ccmpi_distribute(fu,fu_g)
         call ccmpi_distribute(fv,fv_g)
         call ccmpi_distribute(dmdx,dmdx_g)
         call ccmpi_distribute(dmdy,dmdy_g)
         call ccmpi_distributer8(x,x_g)
         call ccmpi_distributer8(y,y_g)
         call ccmpi_distributer8(z,z_g)
         call ccmpi_distribute(rlatt,rlatt_g)
         call ccmpi_distribute(rlongg,rlongg_g)
      else
         call ccmpi_distribute(wts)
         call ccmpi_distribute(em)
         call ccmpi_distribute(emu)
         call ccmpi_distribute(emv)
         call ccmpi_distribute(ax)
         call ccmpi_distribute(bx)
         call ccmpi_distribute(ay)
         call ccmpi_distribute(by)
         call ccmpi_distribute(az)
         call ccmpi_distribute(bz)
         call ccmpi_distribute(f)
         call ccmpi_distribute(fu)
         call ccmpi_distribute(fv)
         call ccmpi_distribute(dmdx)
         call ccmpi_distribute(dmdy)
         call ccmpi_distributer8(x)
         call ccmpi_distributer8(y)
         call ccmpi_distributer8(z)
         call ccmpi_distribute(rlatt)
         call ccmpi_distribute(rlongg)
      end if

      call bounds_setup
      call bounds(em,corner=.true.)
      call boundsuv(emu,emv)
      call boundsuv(ax,bx)
      call boundsuv(ay,by)
      call boundsuv(az,bz)
      call bounds(f,corner=.true.)

      allocate( dpoints(neighnum) )
      allocate( dbuf(0:neighnum) )
      allocate( dindex(0:neighnum) )
      allocate( sextra(neighnum) )
      allocate( dslen(0:neighnum), drlen(0:neighnum) )
      dslen(:) = 0
      drlen(:) = 0

      ! Off processor departure points
      do dproc = 1,neighnum
        iproc = neighlist(dproc)
        allocate( dpoints(dproc)%a(4,bnds(iproc)%len) )
        allocate( dbuf(dproc)%a(4,bnds(iproc)%len) )
        allocate( dbuf(dproc)%b(bnds(iproc)%len) )
        allocate( dindex(dproc)%a(2,bnds(iproc)%len) )
        allocate( sextra(dproc)%a(bnds(iproc)%len) )
      end do
      ! store invalid points in dproc=0
      allocate( dbuf(0)%a(4,1) )
      allocate( dbuf(0)%b(1) )
      allocate( dindex(0)%a(2,1) )


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
         ifullcol(n) = count( colourmask == n )
      end do
      ifullmaxcol = maxval( ifullcol )
      allocate ( iqx(ifullmaxcol,maxcolour) )
      allocate ( iqn(ifullmaxcol,maxcolour), iqe(ifullmaxcol,maxcolour) )
      allocate ( iqw(ifullmaxcol,maxcolour), iqs(ifullmaxcol,maxcolour) )
      ifullcol = 0
      do n = 1,npan
        j = 1
        do i = 1,ipan
          iq = indp(i,j,n)
          ifullcol(colourmask(iq)) = ifullcol(colourmask(iq)) + 1
          iqx(ifullcol(colourmask(iq)),colourmask(iq)) = iq
          iqn(ifullcol(colourmask(iq)),colourmask(iq)) = in(iq)
          iqe(ifullcol(colourmask(iq)),colourmask(iq)) = ie(iq)
          iqw(ifullcol(colourmask(iq)),colourmask(iq)) = iw(iq)
          iqs(ifullcol(colourmask(iq)),colourmask(iq)) = is(iq)
        end do
        j = jpan
        do i = 1,ipan
          iq = indp(i,j,n)
          ifullcol(colourmask(iq)) = ifullcol(colourmask(iq)) + 1
          iqx(ifullcol(colourmask(iq)),colourmask(iq)) = iq
          iqn(ifullcol(colourmask(iq)),colourmask(iq)) = in(iq)
          iqe(ifullcol(colourmask(iq)),colourmask(iq)) = ie(iq)
          iqw(ifullcol(colourmask(iq)),colourmask(iq)) = iw(iq)
          iqs(ifullcol(colourmask(iq)),colourmask(iq)) = is(iq)
        end do
        i = 1
        do j = 2,jpan-1
          iq = indp(i,j,n)
          ifullcol(colourmask(iq)) = ifullcol(colourmask(iq)) + 1
          iqx(ifullcol(colourmask(iq)),colourmask(iq)) = iq
          iqn(ifullcol(colourmask(iq)),colourmask(iq)) = in(iq)
          iqe(ifullcol(colourmask(iq)),colourmask(iq)) = ie(iq)
          iqw(ifullcol(colourmask(iq)),colourmask(iq)) = iw(iq)
          iqs(ifullcol(colourmask(iq)),colourmask(iq)) = is(iq)
        end do
        i = ipan
        do j = 2,jpan-1
          iq = indp(i,j,n)
          ifullcol(colourmask(iq)) = ifullcol(colourmask(iq)) + 1
          iqx(ifullcol(colourmask(iq)),colourmask(iq)) = iq
          iqn(ifullcol(colourmask(iq)),colourmask(iq)) = in(iq)
          iqe(ifullcol(colourmask(iq)),colourmask(iq)) = ie(iq)
          iqw(ifullcol(colourmask(iq)),colourmask(iq)) = iw(iq)
          iqs(ifullcol(colourmask(iq)),colourmask(iq)) = is(iq)
        end do
      end do
      ifullcol_border = ifullcol
      do n = 1,npan
        do j = 2,jpan-1
          do i = 2,ipan-1
            iq = indp(i,j,n)
            ifullcol(colourmask(iq)) = ifullcol(colourmask(iq)) + 1
            iqx(ifullcol(colourmask(iq)),colourmask(iq)) = iq
            iqn(ifullcol(colourmask(iq)),colourmask(iq)) = in(iq)
            iqe(ifullcol(colourmask(iq)),colourmask(iq)) = ie(iq)
            iqw(ifullcol(colourmask(iq)),colourmask(iq)) = iw(iq)
            iqs(ifullcol(colourmask(iq)),colourmask(iq)) = is(iq)
          end do
        end do
      end do

      ltrue = .true. 
!     operator MPI_SUMDR is created based on an external function DRPDR.
      call MPI_OP_CREATE (DRPDR,  ltrue, MPI_SUMDR,  ierr)

      ! prepare comm groups - used by scale-selective filter
#ifdef uniform_decomp
      npta = 6                            ! number of panels per processor
      mproc = nproc                       ! number of processors per panel
      pprocn = 0                          ! start panel
      pprocx = 5                          ! end panel
      hproc = 0                           ! host processor for panel
#else
      npta = max( 6/nproc, 1 )            ! number of panels per processor
      mproc = max( nproc/6, 1 )           ! number of processors per panel
      pprocn = myid*npta/mproc            ! start panel
      pprocx = pprocn + npta - 1          ! end panel
      hproc = pprocn*mproc/npta           ! host processor for panel
#endif

      ! comm between work groups with captain hproc
      colour = hproc
      rank = myid - hproc
      call MPI_Comm_Split(MPI_COMM_WORLD,colour,rank,lcommout,ierr)
      comm_proc = lcommout
      
      ! comm between columns in work group
      colour = ioff
      rank = joff/jpan
      lcommin = comm_proc
      call MPI_Comm_Split(lcommin,colour,rank,lcommout,ierr)
      comm_cols = lcommout
      
      ! comm between rows in work group      
      colour = joff
      rank = ioff/ipan
      call MPI_Comm_Split(lcommin,colour,rank,lcommout,ierr)
      comm_rows = lcommout
      
      if ( myid==hproc ) then
         if ( ioff/=0 .or. joff/=0 ) then
            write(6,*) "ERROR: hproc incorrectly assigned"
            call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
         end if
      end if
      
      ! prep RMA windows for gathermap
      if ( nproc>1 ) then
         allocate(specstore(ifull,kx))
         !call MPI_Info_create(info,ierr)
         !call MPI_Info_set(info,"no_locks","true",ierr)
         !call MPI_Info_set(info,"same_size","true",ierr)
         !call MPI_Info_set(info,"same_disp_unit","true",ierr)
         call MPI_Type_size(ltype,asize,ierr)
         wsize = asize*ifull*kx
         call MPI_Win_create(specstore,wsize,asize,MPI_INFO_NULL,MPI_COMM_WORLD,localwin,ierr)
         !call MPI_Info_free(info,ierr)
      end if
      
   return
   end subroutine ccmpi_setup
   
   subroutine ccmpi_distribute2(af,a1)
      ! Convert standard 1D arrays to face form and distribute to processors
      real, dimension(ifull), intent(out) :: af
      real, dimension(ifull_g), intent(in), optional :: a1
      integer(kind=4) :: ierr

      call START_LOG(distribute_begin)

      ! Copy internal region
      if ( myid == 0 ) then
         if ( .not. present(a1) ) then
            write(6,*) "Error: ccmpi_distribute argument required on proc 0"
            call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
         end if
         call host_distribute2(af,a1)
      else
         call proc_distribute2(af)
      end if

      call END_LOG(distribute_end)

   end subroutine ccmpi_distribute2

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
      integer(kind=4) :: ierr, lsize
      integer :: npoff, ipoff, jpoff ! Offsets for target
      integer :: slen

      ! map array in order of processor rank
      do iproc = 0,nproc-1
#ifdef uniform_decomp
         call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,ipan,jpan)
#else
         call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan)
#endif
         do n = 1,npan
            do j = 1,jpan
               iq = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g
               slen = (j-1)*ipan + (n-1)*ipan*jpan
               sbuf(slen+1:slen+ipan,iproc) = a1(iq+1:iq+ipan)
            end do
         end do
      end do

      lsize = ifull
      call MPI_Scatter(sbuf,lsize,ltype,af,lsize,ltype,0_4,MPI_COMM_WORLD,ierr)

   end subroutine host_distribute2

   subroutine proc_distribute2(af)
      ! Convert standard 1D arrays to face form and distribute to processors
      real, dimension(ifull), intent(out) :: af
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      integer(kind=4) :: ierr, lsize
      real, dimension(0,0) :: sbuf

      lsize = ifull
      call MPI_Scatter(sbuf,lsize,ltype,af,lsize,ltype,0_4,MPI_COMM_WORLD,ierr)

   end subroutine proc_distribute2

   subroutine ccmpi_distribute2r8(af,a1)
      ! Convert standard 1D arrays to face form and distribute to processors
      real(kind=8), dimension(ifull), intent(out) :: af
      real(kind=8), dimension(ifull_g), intent(in), optional :: a1
      integer(kind=4) :: ierr

      call START_LOG(distribute_begin)

      if ( myid == 0 ) then
         if ( .not. present(a1) ) then
            write(6,*) "Error: ccmpi_distribute argument required on proc 0"
            call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
         end if
         call host_distribute2r8(af,a1)
      else
         call proc_distribute2r8(af)
      end if

      call END_LOG(distribute_begin)
      
   end subroutine ccmpi_distribute2r8

   subroutine host_distribute2r8(af,a1)
      ! Convert standard 1D arrays to face form and distribute to processors
      real(kind=8), dimension(ifull), intent(out) :: af
      real(kind=8), dimension(ifull_g), intent(in) :: a1
      integer :: j, n, iq, iproc
      integer(kind=4) :: ierr, lsize
      real(kind=8), dimension(ifull,0:nproc-1) :: sbuf
      integer :: npoff, ipoff, jpoff ! Offsets for target
      integer :: slen

      ! map array in order of processor rank
      do iproc = 0,nproc-1
#ifdef uniform_decomp
         call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,ipan,jpan)
#else
         call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan)
#endif
         do n = 1,npan
            do j = 1,jpan
               iq = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g
               slen = (j-1)*ipan + (n-1)*ipan*jpan
               sbuf(slen+1:slen+ipan,iproc) = a1(iq+1:iq+ipan)
            end do
         end do
      end do

      lsize = ifull
      call MPI_Scatter(sbuf,lsize,MPI_DOUBLE_PRECISION,af,lsize,MPI_DOUBLE_PRECISION,0_4,MPI_COMM_WORLD,ierr)

   end subroutine host_distribute2r8
   
   subroutine proc_distribute2r8(af)
      ! Convert standard 1D arrays to face form and distribute to processors
      real(kind=8), dimension(ifull), intent(out) :: af
      integer(kind=4) :: ierr, lsize
      real(kind=8), dimension(0,0) :: sbuf

      lsize = ifull
      call MPI_Scatter(sbuf,lsize,MPI_DOUBLE_PRECISION,af,lsize,MPI_DOUBLE_PRECISION,0_4,MPI_COMM_WORLD,ierr)

   end subroutine proc_distribute2r8

   subroutine ccmpi_distribute2i(af,a1)
      ! Convert standard 1D arrays to face form and distribute to processors
      integer, dimension(ifull), intent(out) :: af
      integer, dimension(ifull_g), intent(in), optional :: a1
      integer(kind=4) :: ierr

      call START_LOG(distribute_begin)

      if ( myid == 0 ) then
         if ( .not. present(a1) ) then
            write(6,*) "Error: ccmpi_distribute argument required on proc 0"
            call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
         end if
         call host_distribute2i(af,a1)
      else
         call proc_distribute2i(af)
      end if

      call END_LOG(distribute_end)
      
   end subroutine ccmpi_distribute2i

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
      integer(kind=4) :: ierr, lsize
      integer, dimension(ifull,0:nproc-1) :: sbuf
      integer :: npoff, ipoff, jpoff ! Offsets for target
      integer :: slen

      ! map array in order of processor rank
      do iproc = 0,nproc-1
#ifdef uniform_decomp
         call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,ipan,jpan)
#else
         call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan)
#endif
         do n = 1,npan
            do j = 1,jpan
               iq = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g
               slen = (j-1)*ipan + (n-1)*ipan*jpan
               sbuf(slen+1:slen+ipan,iproc) = a1(iq+1:iq+ipan)
            end do
         end do
      end do

      lsize = ifull
      call MPI_Scatter(sbuf,lsize,ltype,af,lsize,ltype,0_4,MPI_COMM_WORLD,ierr)
 
   end subroutine host_distribute2i

   subroutine proc_distribute2i(af)
      ! Convert standard 1D arrays to face form and distribute to processors
      integer, dimension(ifull), intent(out) :: af
#ifdef i8r8
      integer(kind=4),parameter :: ltype = MPI_INTEGER8
#else
      integer(kind=4),parameter :: ltype = MPI_INTEGER
#endif
      integer(kind=4) :: ierr, lsize
      integer, dimension(0,0) :: sbuf

      lsize = ifull
      call MPI_Scatter(sbuf,lsize,ltype,af,lsize,ltype,0_4,MPI_COMM_WORLD,ierr)
 
   end subroutine proc_distribute2i

   subroutine ccmpi_distribute3(af,a1)
      ! Convert standard 1D arrays to face form and distribute to processors
      ! This is also used for tracers, so second dimension is not necessarily
      ! the number of levels
      real, dimension(:,:), intent(out) :: af
      real, dimension(:,:), intent(in), optional :: a1
      integer(kind=4) :: ierr

      call START_LOG(distribute_begin)

      if ( myid == 0 ) then
         if ( .not. present(a1) ) then
            write(6,*) "Error: ccmpi_distribute argument required on proc 0"
            call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
         end if
         call host_distribute3(af,a1)
      else
         call proc_distribute3(af)
      end if

      call END_LOG(distribute_end)
      
   end subroutine ccmpi_distribute3

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
      integer(kind=4) :: ierr, lsize
      real, dimension(ifull,size(af,2),0:nproc-1) :: sbuf
      real, dimension(ifull,size(af,2)) :: aftemp
      integer :: npoff, ipoff, jpoff ! Offsets for target
      integer :: slen, kx

      kx = size(af,2)

      ! map array in order of processor rank
      do iproc = 0,nproc-1
#ifdef uniform_decomp
         call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,ipan,jpan)
#else
         call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan)
#endif
         do k = 1,kx
            do n = 1,npan
               do j = 1,jpan
                  iq = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g
                  slen = (j-1)*ipan + (n-1)*ipan*jpan
                  sbuf(slen+1:slen+ipan,k,iproc) = a1(iq+1:iq+ipan,k)
               end do
            end do
         end do
      end do

      lsize = ifull*kx
      call MPI_Scatter(sbuf,lsize,ltype,aftemp,lsize,ltype,0_4,MPI_COMM_WORLD,ierr)
      af(1:ifull,1:kx) = aftemp(:,:)

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
      integer(kind=4) :: ierr, lsize
      real, dimension(0,0,0) :: sbuf
      real, dimension(ifull,size(af,2)) :: aftemp
      integer :: kx

      kx = size(af,2)
      lsize = ifull*kx
      call MPI_Scatter(sbuf,lsize,ltype,aftemp,lsize,ltype,0_4,MPI_COMM_WORLD,ierr)
      af(1:ifull,1:kx) = aftemp(:,:)

   end subroutine proc_distribute3

   subroutine ccmpi_distribute3i(af,a1)
      ! Convert standard 1D arrays to face form and distribute to processors
      ! This is also used for tracers, so second dimension is not necessarily
      ! the number of levels
      integer, dimension(:,:), intent(out) :: af
      integer, dimension(:,:), intent(in), optional :: a1
      integer(kind=4) :: ierr

      call START_LOG(distribute_begin)
      
      if ( myid == 0 ) then
         if ( .not. present(a1) ) then
            write(6,*) "Error: ccmpi_distribute argument required on proc 0"
            call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
         end if
         call host_distribute3i(af,a1)
      else
         call proc_distribute3i(af)
      end if

      call END_LOG(distribute_end)
      
   end subroutine ccmpi_distribute3i

   subroutine host_distribute3i(af,a1)
      ! Convert standard 1D arrays to face form and distribute to processors
      ! This is also used for tracers, so second dimension is not necessarily
      ! the number of levels
      integer, dimension(:,:), intent(out) :: af
      integer, dimension(:,:), intent(in) :: a1
#ifdef i8r8
      integer(kind=4),parameter :: ltype = MPI_INTEGER8
#else
      integer(kind=4),parameter :: ltype = MPI_INTEGER
#endif
      integer :: j, n, k, iq, iproc
      integer(kind=4) :: ierr, lsize
      integer, dimension(ifull,size(af,2),0:nproc-1) :: sbuf
      integer, dimension(ifull,size(af,2)) :: aftemp
      integer :: npoff, ipoff, jpoff ! Offsets for target
      integer :: slen, kx

      kx = size(af,2)

      ! map array in order of processor rank
      do iproc = 0,nproc-1
#ifdef uniform_decomp
         call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,ipan,jpan)
#else
         call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan)
#endif
         do k = 1,kx
            do n = 1,npan
               do j = 1,jpan
                  iq = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g
                  slen = (j-1)*ipan + (n-1)*ipan*jpan
                  sbuf(slen+1:slen+ipan,k,iproc) = a1(iq+1:iq+ipan,k)
               end do
            end do
         end do
      end do

      lsize = ifull*kx
      call MPI_Scatter(sbuf,lsize,ltype,aftemp,lsize,ltype,0_4,MPI_COMM_WORLD,ierr)      
      af(1:ifull,1:kx) = aftemp(:,:)

   end subroutine host_distribute3i

   subroutine proc_distribute3i(af)
      ! Convert standard 1D arrays to face form and distribute to processors
      ! This is also used for tracers, so second dimension is not necessarily
      ! the number of levels
      integer, dimension(:,:), intent(out) :: af
#ifdef i8r8
      integer(kind=4),parameter :: ltype = MPI_INTEGER8
#else
      integer(kind=4),parameter :: ltype = MPI_INTEGER
#endif
      integer(kind=4) :: ierr, lsize
      integer, dimension(0,0,0) :: sbuf
      integer, dimension(ifull,size(af,2)) :: aftemp
      integer :: kx

      kx = size(af,2)
      lsize = ifull*kx
      call MPI_Scatter(sbuf,lsize,ltype,aftemp,lsize,ltype,0_4,MPI_COMM_WORLD,ierr)      
      af(1:ifull,1:kx) = aftemp(:,:)

   end subroutine proc_distribute3i   

   subroutine ccmpi_gather2(a,ag)
      ! Collect global arrays.

      real, dimension(ifull), intent(in) :: a
      real, dimension(ifull_g), intent(out), optional :: ag
      integer(kind=4) :: ierr

      call START_LOG(gather_begin)

      if ( myid == 0 ) then
         if ( .not. present(ag) ) then
            write(6,*) "Error: ccmpi_gather argument required on proc 0"
            call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
         end if
         call host_gather2(a,ag)
      else
         call proc_gather2(a)
      end if

      call END_LOG(gather_end)

   end subroutine ccmpi_gather2

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
      integer(kind=4) :: ierr, lsize
      real, dimension(ifull,0:nproc-1) :: abuf
      integer :: ipoff, jpoff, npoff
      integer :: j, n, iq, iqg

      lsize = ifull
      call MPI_Gather(a,lsize,ltype,abuf,lsize,ltype,0_4,MPI_COMM_WORLD,ierr)

      ! map array in order of processor rank
      do iproc = 0,nproc-1
#ifdef uniform_decomp
         call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,ipan,jpan)
#else
         call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan)
#endif
         do n = 1,npan
            ! Use the face indices for unpacking
            do j = 1,jpan
               ! Global indices are i+ipoff, j+jpoff, n-npoff
               iqg = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g ! True global 1D index
               iq = (j-1)*ipan + (n-1)*ipan*jpan
               ag(iqg+1:iqg+ipan) = abuf(iq+1:iq+ipan,iproc)
            end do
         end do
      end do

   end subroutine host_gather2
   
   subroutine proc_gather2(a)
      real, dimension(ifull), intent(in) :: a
#ifdef i8r8
      integer(kind=4),parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4),parameter :: ltype = MPI_REAL
#endif
      integer(kind=4) :: ierr, lsize
      real, dimension(0,0) :: abuf

      lsize = ifull
      call MPI_Gather(a,lsize,ltype,abuf,lsize,ltype,0_4,MPI_COMM_WORLD,ierr)

   end subroutine proc_gather2

   subroutine ccmpi_gather3(a,ag)
      ! Collect global arrays.

      real, dimension(:,:), intent(in) :: a
      real, dimension(:,:), intent(out), optional :: ag
      integer(kind=4) :: ierr

      call START_LOG(gather_begin)

      if ( myid == 0 ) then
         if ( .not. present(ag) ) then
            write(6,*) "Error: ccmpi_gather argument required on proc 0"
            call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
         end if
         call host_gather3(a,ag)
      else
         call proc_gather3(a)
      end if
      
      call END_LOG(gather_end)

   end subroutine ccmpi_gather3

   subroutine host_gather3(a,ag)
      real, dimension(:,:), intent(in) :: a
      real, dimension(:,:), intent(out) :: ag
      integer :: iproc
#ifdef i8r8
      integer(kind=4),parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4),parameter :: ltype = MPI_REAL
#endif
      integer(kind=4) :: ierr, lsize
      real, dimension(ifull,size(a,2),0:nproc-1) :: abuf
      real, dimension(ifull,size(a,2)) :: atemp
      integer :: ipoff, jpoff, npoff
      integer :: j, n, k, iq, iqg, kx

      kx = size(a,2)
      lsize = ifull*kx
      atemp(:,:) = a(1:ifull,1:kx)
      call MPI_Gather(atemp,lsize,ltype,abuf,lsize,ltype,0_4,MPI_COMM_WORLD,ierr)

      ! map array in order of processor rank
      do iproc = 0,nproc-1
#ifdef uniform_decomp
         call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,ipan,jpan)
#else
         call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan)
#endif
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

   end subroutine host_gather3
   
   subroutine proc_gather3(a)
      real, dimension(:,:), intent(in) :: a
#ifdef i8r8
      integer(kind=4),parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4),parameter :: ltype = MPI_REAL
#endif
      integer(kind=4) :: ierr, lsize
      real, dimension(0,0,0) :: abuf
      real, dimension(ifull,size(a,2)) :: atemp
      integer :: kx

      kx = size(a,2)
      lsize = ifull*kx
      atemp(:,:) = a(1:ifull,1:kx)
      call MPI_Gather(atemp,lsize,ltype,abuf,lsize,ltype,0_4,MPI_COMM_WORLD,ierr)

   end subroutine proc_gather3

   subroutine ccmpi_gatherall2(a,ag)
      ! Collect global arrays.

      real, dimension(:), intent(in) :: a
      real, dimension(:), intent(out) :: ag
      real, dimension(ifull,0:nproc-1) :: abuf
      integer :: iproc
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      integer(kind=4) :: ierr, lsize
      integer :: ipoff, jpoff, npoff
      integer :: j, n, iq, iqg

      call START_LOG(gather_begin)

      lsize = ifull
      call MPI_AllGather(a,lsize,ltype,abuf,lsize,ltype,MPI_COMM_WORLD,ierr)

      ! map array in order of processor rank
      do iproc = 0,nproc-1
#ifdef uniform_decomp
         call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,ipan,jpan)
#else
         call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan)
#endif
         do n = 1,npan
            do j = 1,jpan
               ! Global indices are i+ipoff, j+jpoff, n-npoff
               iqg = ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g ! True global 1D index
               iq = (j-1)*ipan + (n-1)*ipan*jpan
               ag(iqg+1:iqg+ipan) = abuf(iq+1:iq+ipan,iproc)
            end do
         end do
      end do

      call END_LOG(gather_end)

   end subroutine ccmpi_gatherall2
   
   subroutine ccmpi_gatherall3(a,ag)
      ! Collect global arrays.

      real, dimension(:,:), intent(in) :: a
      real, dimension(:,:), intent(out) :: ag
      integer :: iproc
#ifdef i8r8
      integer(kind=4),parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4),parameter :: ltype = MPI_REAL
#endif
      integer(kind=4) :: ierr, lsize
      real, dimension(ifull,size(a,2),0:nproc-1) :: abuf
      real, dimension(ifull,size(a,2)) :: atemp
      integer :: ipoff, jpoff, npoff
      integer :: j, n, k, iq, iqg, kx

      call START_LOG(gather_begin)

      kx = size(a,2)
      lsize = ifull*kx
      atemp(:,:) = a(1:ifull,1:kx)
      call MPI_AllGather(atemp,lsize,ltype,abuf,lsize,ltype,MPI_COMM_WORLD,ierr)

      ! map array in order of processor rank
      do iproc = 0,nproc-1
#ifdef uniform_decomp
         call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,ipan,jpan)
#else
         call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan)
#endif
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

      call END_LOG(gather_end)

   end subroutine ccmpi_gatherall3

   subroutine ccmpi_gathermap2(a,kref)

      integer, intent(in) :: kref
      real, dimension(ifull), intent(in) :: a
      real, dimension(ifull,size(specmap)) :: abuf 
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      integer(kind=4) :: ierr, lsize, itest
      integer(kind=MPI_ADDRESS_KIND) :: displ
      integer :: ncount, w, iproc, n, iq
      integer :: ipoff, jpoff, npoff
      integer :: ipak, jpak
      
      if ( nproc==1 ) then
         do n = 1,npan
            iq = (n-1)*ipan*jpan
            globalpack(0,0,n-noff)%localdata(:,:,kref+1) = &
               reshape( a(iq+1:iq+ipan*jpan), (/ ipan, jpan /) )
         end do
         return
      end if
   
      call START_LOG(gatherrma_begin)
   
      ncount = size(specmap)
      specstore(1:ifull,1) = a(1:ifull)
   
      lsize = ifull
      displ = 0
      call MPI_Win_fence( MPI_MODE_NOPRECEDE, localwin, ierr )
      do w = 1,ncount
         call MPI_Get(abuf(:,w),lsize,ltype,specmap(w),displ,lsize,ltype,localwin,ierr)
      end do
      itest = ior( MPI_MODE_NOSUCCEED, MPI_MODE_NOPUT )
      call MPI_Win_fence( itest, localwin, ierr )
   
      do w = 1,ncount
         iproc = specmap(w)
#ifdef uniform_decomp
         call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,ipan,jpan)
#else
         call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan)
#endif
         ipak = ipoff/ipan
         jpak = jpoff/jpan
         do n = 1,npan
            ! Global indices are i+ipoff, j+jpoff, n-npoff
            iq = (n-1)*ipan*jpan
            globalpack(ipak,jpak,n-npoff)%localdata(:,:,kref+1) = &
               reshape( abuf(iq+1:iq+ipan*jpan,w), (/ ipan, jpan /) )
         end do
      end do
      
      call END_LOG(gatherrma_end)
   
   end subroutine ccmpi_gathermap2

   subroutine ccmpi_gathermap3(a,kref)

      integer, intent(in) :: kref
      real, dimension(:,:), intent(in) :: a
      real, dimension(ifull*size(a,2),size(specmap)) :: abuf 
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      integer(kind=4) :: ierr, lsize, itest
      integer(kind=MPI_ADDRESS_KIND) :: displ
      integer :: ncount, w, iproc, k, n, iq, kx
      integer :: ipoff, jpoff, npoff
      integer :: ipak, jpak
      
      kx = size(a,2)
      
      if ( nproc==1 ) then
         do k = 1,kx
            do n = 1,npan
               iq = (n-1)*ipan*jpan
               globalpack(0,0,n-noff)%localdata(:,:,kref+k) = &
                  reshape( a(iq+1:iq+ipan*jpan,k), (/ ipan, jpan /) )
            end do
         end do
         return
      end if
   
      call START_LOG(gatherrma_begin)
   
      ncount = size(specmap)
      
      if ( kx>size(specstore,2) ) then
         write(6,*) "ERROR: gathermap array is too big for window buffer"
         call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
      end if
      
      specstore(1:ifull,1:kx) = a(1:ifull,:)

      lsize = ifull*kx
      displ = 0   
      call MPI_Win_fence(MPI_MODE_NOPRECEDE,localwin,ierr)
      do w = 1,ncount
         call MPI_Get(abuf(:,w),lsize,ltype,specmap(w),displ,lsize,ltype,localwin,ierr)
      end do
      itest = ior(MPI_MODE_NOSUCCEED,MPI_MODE_NOPUT)
      call MPI_Win_fence(itest,localwin,ierr)
   
      do w = 1,ncount
         iproc = specmap(w)
#ifdef uniform_decomp
         call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,ipan,jpan)
#else
         call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan)
#endif
         ipak = ipoff/ipan
         jpak = jpoff/jpan
         do k = 1,kx
            do n = 1,npan
               ! Global indices are i+ipoff, j+jpoff, n-npoff
               iq = (n-1)*ipan*jpan + (k-1)*ifull
               globalpack(ipak,jpak,n-npoff)%localdata(:,:,kref+k) = &
                  reshape( abuf(iq+1:iq+ipan*jpan,w), (/ ipan, jpan /) )
            end do
         end do
      end do
      
      call END_LOG(gatherrma_end)
   
   end subroutine ccmpi_gathermap3
   
   subroutine setglobalpack(datain,jbeg,ibeg,iend,k,trans)
   
      ! This subroutine assigns a value to a gridpoint
      ! in the global sparse array
   
      integer, intent(in) :: jbeg, ibeg, iend, k
      integer :: il2, iqg, im1, jm1, ilen, jlen, ijlen, iq
      integer :: b_n, b_ipak, b_jpak, b_iloc, b_jloc
      integer :: e_n, e_ipak, e_jpak, e_iloc, e_jloc
      integer :: s_ipak, s_jpak, s_iloc, s_jloc
      integer :: c_ipak, c_jpak
      integer(kind=4) :: ierr
      real, dimension(:), intent(in) :: datain
      logical :: transp
      logical, intent(in), optional :: trans
      
      transp = .false.
      if ( present(trans) ) then
         transp = trans
      end if
      
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
      
      if ( b_n/=e_n ) then
         write(6,*) "ERROR: getglobalpack requires ibeg and iend to belong to the same face"
         call MPI_abort(MPI_COMM_WORLD,-1,ierr)
      end if

      if ( e_jpak>=b_jpak) then
         s_jpak = 1
      else
         s_jpak = -1
      end if
      if ( e_ipak>=b_ipak) then
         s_ipak = 1
      else
         s_ipak = -1
      end if
      if ( e_jloc>=b_jloc) then
         s_jloc = 1
      else
         s_jloc = -1
      end if
      if ( e_iloc>=b_iloc) then
         s_iloc = 1
      else
         s_iloc = -1
      end if

      ilen = abs(e_iloc-b_iloc) + 1
      jlen = abs(e_jloc-b_jloc) + 1
      ijlen = ilen*jlen
      
      if ( transp ) then
         iq = jbeg - 1
         do c_ipak = b_ipak,e_ipak,s_ipak
            do c_jpak = b_jpak,e_jpak,s_jpak
               globalpack(c_ipak,c_jpak,b_n)%localdata(b_iloc:e_iloc:s_iloc,b_jloc:e_jloc:s_jloc,k) = &
                 transpose( reshape( datain(iq+1:iq+ijlen), (/ jlen, ilen /) ) )
               iq = iq + ijlen
            end do
         end do
      else
         iq = jbeg - 1
         do c_jpak = b_jpak,e_jpak,s_jpak
            do c_ipak = b_ipak,e_ipak,s_ipak
               globalpack(c_ipak,c_jpak,b_n)%localdata(b_iloc:e_iloc:s_iloc,b_jloc:e_jloc:s_jloc,k) = &
                 reshape( datain(iq+1:iq+ijlen), (/ ilen, jlen /) )
               iq = iq + ijlen
            end do
         end do
      end if
   
   end subroutine setglobalpack
   
   subroutine getglobalpack(dataout,jbeg,ibeg,iend,k,trans)
   
      ! This subroutine returns a value from a gridpoint
      ! in the global sparse array

      integer, intent(in) :: jbeg, ibeg, iend, k
      integer :: il2, iqg, im1, jm1, ilen, jlen, ijlen, iq
      integer :: b_n, b_ipak, b_jpak, b_iloc, b_jloc
      integer :: e_n, e_ipak, e_jpak, e_iloc, e_jloc
      integer :: s_ipak, s_jpak, s_iloc, s_jloc
      integer :: c_ipak, c_jpak
      integer(kind=4) :: ierr
      real, dimension(:), intent(out) :: dataout
      logical :: transp
      logical, intent(in), optional :: trans
      
      transp = .false.
      if ( present(trans) ) then
         transp = trans
      end if

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
      
      if ( b_n/=e_n ) then
         write(6,*) "ERROR: getglobalpack requires ibeg and iend to belong to the same face"
         call MPI_abort(MPI_COMM_WORLD,-1,ierr)
      end if
      
      if ( e_jpak>=b_jpak) then
         s_jpak = 1
      else
         s_jpak = -1
      end if
      if ( e_ipak>=b_ipak) then
         s_ipak = 1
      else
         s_ipak = -1
      end if
      if ( e_jloc>=b_jloc) then
         s_jloc = 1
      else
         s_jloc = -1
      end if
      if ( e_iloc>=b_iloc) then
         s_iloc = 1
      else
         s_iloc = -1
      end if

      ilen = abs(e_iloc-b_iloc) + 1
      jlen = abs(e_jloc-b_jloc) + 1
      ijlen = ilen*jlen
      
      if ( transp ) then
         iq = jbeg - 1
         do c_ipak = b_ipak,e_ipak,s_ipak
            do c_jpak = b_jpak,e_jpak,s_jpak
               dataout(iq+1:iq+ijlen) = reshape( transpose(                                            &
                 globalpack(c_ipak,c_jpak,b_n)%localdata(b_iloc:e_iloc:s_iloc,b_jloc:e_jloc:s_jloc,k)  &
                 ), (/ ijlen /) )
               iq = iq + ijlen
            end do
         end do
      else
         iq = jbeg - 1
         do c_jpak = b_jpak,e_jpak,s_jpak
            do c_ipak = b_ipak,e_ipak,s_ipak
               dataout(iq+1:iq+ijlen) = reshape(                                                       &
                 globalpack(c_ipak,c_jpak,b_n)%localdata(b_iloc:e_iloc:s_iloc,b_jloc:e_jloc:s_jloc,k), &
                 (/ ijlen /) )
               iq = iq + ijlen
            end do
         end do
      end if
      
   end subroutine getglobalpack
   
   subroutine copyglobalpack(krefin,krefout,kx)

      ! This routine copies one section of the global sparse array
      ! to another section.  Note that it only copies the memory
      ! assigned by gathermap.  specmap needs to be replaced with
      ! spectmapext to copy all parts of the global sparse array.
   
      integer, intent(in) :: krefin,krefout,kx
      integer :: w, n, ncount, iproc, ipak, jpak
      integer :: ipoff, jpoff, npoff
   
      ncount = size(specmap)
      do w = 1,ncount
         iproc = specmap(w)
#ifdef uniform_decomp
         call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,ipan,jpan)
#else
         call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan)
#endif
         ipak = ipoff/ipan
         jpak = jpoff/jpan
         do n = 1,npan
            ! Global indices are i+ipoff, j+jpoff, n-npoff
            globalpack(ipak,jpak,n-npoff)%localdata(:,:,krefout+1:krefout+kx) = &
               globalpack(ipak,jpak,n-npoff)%localdata(:,:,krefin+1:krefin+kx)
         end do
      end do
   
   end subroutine copyglobalpack

   subroutine allocateglobalpack(kx)
   
      ! This allocates global sparse arrays for the digital filter.
      ! Usually this is 1:kl or 1:ol in size, but for some configurations
      ! we need to store the original values and hence use 1:2*kl or 1:2*ol.
      ! Also, the 0 index is to store the sum term for the digital filter.
   
      integer, intent(in) :: kx
      integer :: ncount, w, ipak, jpak, n, iproc
      integer :: ipoff, jpoff, npoff
      
      ! allocate globalpack arrays for 1D scale-selective filter
      allocate(globalpack(0:nxproc-1,0:nyproc-1,0:5))
      ncount = size(specmapext)
      do w = 1,ncount
         iproc = specmapext(w)
#ifdef uniform_decomp
         call proc_region_dix(iproc,ipoff,jpoff,npoff,nxproc,ipan,jpan)
#else
         call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan)
#endif
         ! Global indices are i+ipoff, j+jpoff, n-npoff
         ipak = ipoff/ipan
         jpak = jpoff/jpan
         do n = 1,npan
            allocate(globalpack(ipak,jpak,n-npoff)%localdata(ipan,jpan,0:kx))
            globalpack(ipak,jpak,n-npoff)%localdata = 0.
         end do
      end do
      
      deallocate(specmapext) ! not needed after allocation of global sparse arrays
   
   end subroutine allocateglobalpack
   
   subroutine ccmpi_filewincreate(kx)
   
      integer, intent(in) :: kx
      integer(kind=4) :: asize, ierr
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      integer(kind=MPI_ADDRESS_KIND) :: wsize
      
      if ( nproc>1 ) then
         if ( myid<fnresid ) then
            allocate( filestore(pil*pjl*pnpan,kx) )
         else
            allocate( filestore(0,0) )
         end if
         !call MPI_Info_create(info,ierr)
         !call MPI_Info_set(info,"no_locks","true",ierr)
         call MPI_Type_size(ltype, asize, ierr)
         wsize = asize*pil*pjl*pnpan*kx
         call MPI_Win_create(filestore, wsize, asize, MPI_INFO_NULL, MPI_COMM_WORLD, filewin, ierr)
         !call MPI_Info_free(info,ierr)
      end if
   
   end subroutine ccmpi_filewincreate
   
   subroutine ccmpi_filewinfree
   
      integer(kind=4) ierr
   
      if ( nproc>1 ) then
         deallocate( filestore )
         call MPI_Win_Free( filewin, ierr )
      end if
   
   end subroutine ccmpi_filewinfree
   
   subroutine ccmpi_filewinget2(abuf,sinp)
   
      integer n, w, ncount, nlen
      integer ca, cc, ipf
      integer(kind=4) :: lsize, ierr, itest
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      integer(kind=MPI_ADDRESS_KIND) :: displ
      real, dimension(:), intent(in) :: sinp
      real, dimension(:,:,:), intent(out) :: abuf 
   
      if ( nproc==1 ) then
         do ipf = 0,fncount-1
            do n = 0,pnpan-1
               ca = n*pil*pjl
               cc = n*pil*pjl + pil*pjl*pnpan*ipf
               abuf(1+ca:pil*pjl+ca,1,ipf+1) = sinp(1+cc:pil*pjl+cc)
            end do
         end do
         return
      end if
   
      call START_LOG(gatherrma_begin)
   
      ncount = size(filemap)
      nlen = pil*pjl*pnpan
      lsize = nlen
      displ = 0_4
      itest = ior( MPI_MODE_NOSUCCEED, MPI_MODE_NOPUT )
      
      do ipf = 0,fncount-1
          
         if ( myid<fnresid ) then
            cc = nlen*ipf             
            filestore(1:nlen,1) = sinp(1+cc:nlen+cc)
         end if
   
         call MPI_Win_fence(MPI_MODE_NOPRECEDE, filewin, ierr)
         do w = 1,ncount
            call MPI_Get(abuf(:,w,ipf+1), lsize, ltype, filemap(w), displ, lsize, ltype, filewin, ierr)
         end do
         call MPI_Win_fence(itest, filewin, ierr)
   
      end do
      
      call END_LOG(gatherrma_end)
      
   end subroutine ccmpi_filewinget2

   subroutine ccmpi_filewinget3(abuf,sinp)
   
      integer n, w, ncount, nlen, kx, k
      integer ca, cc, ipf
      integer(kind=4) :: lsize, ierr, itest
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      integer(kind=MPI_ADDRESS_KIND) :: displ
      real, dimension(:,:), intent(in) :: sinp
      real, dimension(:,:,:,:), intent(out) :: abuf
      real, dimension(pil*pjl*pnpan,size(sinp,2),size(filemap)) :: bbuf
   
      kx = size(sinp,2)
      
      if ( nproc==1 ) then
         do ipf = 0,fncount-1
            do k = 1,kx
               do n = 0,pnpan-1
                  ca = n*pil*pjl
                  cc = n*pil*pjl + pil*pjl*pnpan*ipf
                  abuf(1+ca:pil*pjl+ca,1,ipf+1,k) = sinp(1+cc:pil*pjl+cc,k)
               end do
            end do
         end do
         return
      end if
   
      call START_LOG(gatherrma_begin)

      if ( kx>size(filestore,2) .and. myid<fnresid ) then
         write(6,*) "ERROR: Size of file window is too small to support input array size"
         write(6,*) "Window levels ",size(filestore,2)
         write(6,*) "Input levels ",kx
         call ccmpi_abort(-1)
      end if
      
      ncount = size(filemap)
      nlen = pil*pjl*pnpan
      lsize = nlen*kx
      displ = 0_4
      itest = ior( MPI_MODE_NOSUCCEED, MPI_MODE_NOPUT )
      
      do ipf = 0,fncount-1
          
         if ( myid<fnresid ) then
            cc = nlen*ipf             
            filestore(1:nlen,1:kx) = sinp(1+cc:nlen+cc,1:kx)
         end if
   
         call MPI_Win_fence(MPI_MODE_NOPRECEDE, filewin, ierr)
         do w = 1,ncount
            call MPI_Get(bbuf(:,:,w), lsize, ltype, filemap(w), displ, lsize, ltype, filewin, ierr)
         end do
         call MPI_Win_fence(itest, filewin, ierr)
         do w = 1,ncount
            abuf(1:nlen,w,ipf+1,1:kx) = bbuf(1:nlen,1:kx,w)
         end do

      end do
      
      call END_LOG(gatherrma_end)
      
   end subroutine ccmpi_filewinget3
   
   subroutine ccmpi_filewinunpack(sout,abuf)

      integer :: ncount, ipf, w, ip, n, no, ca, cb, cc
      real, dimension(-1:,-1:,0:), intent(inout) :: sout
      real, dimension(pil*pjl*pnpan,size(filemap),fncount), intent(in) :: abuf

      ncount = size(filemap)

      do ipf = 1,fncount ! fncount=fnproc/fnresid
         do w = 1,ncount
            ip = filemap(w) + (ipf-1)*fnresid
            do n = 0,pnpan-1
               no = n - pnoff(ip) + 1
               ca = pioff(ip,no)
               cb = pjoff(ip,no)
               cc = n*pil*pjl
               sout(1+ca:pil+ca,1+cb:pjl+cb,no) = reshape( abuf(1+cc:pil*pjl+cc,w,ipf), (/ pil, pjl /) )
            end do
         end do
      end do
      
   end subroutine ccmpi_filewinunpack
   
   subroutine bounds_setup

      use indices_m
      
      integer :: n, i, j, iq, iqq, mycol, ncount
      integer :: iproc, rproc, sproc
      integer(kind=4), dimension(:,:), allocatable :: status
      integer(kind=4) :: ierr, itag=0, lcount
      integer(kind=4) :: llen, lproc
      integer, dimension(:,:), allocatable :: dums, dumr
      integer, dimension(:,:), allocatable :: dumsb, dumrb
      integer :: iqg, iql, iloc, jloc, nloc, icol
      integer :: iext, iextu, iextv
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_INTEGER8
#else
      integer(kind=4), parameter :: ltype = MPI_INTEGER
#endif
      logical :: swap
      logical(kind=4), dimension(:,:), allocatable :: dumsl, dumrl

      ! Just set values that point to values within own processors region.
      ! Other values are set up later
      ! Set initial values to make sure missing values are caught properly.
!cdir iexpand(indp, indg, indv_mpi)
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

      bnds(:)%len = 0
      bnds(:)%rlenh = 0
      bnds(:)%slenh = 0
      bnds(:)%rlen = 0
      bnds(:)%slen = 0
      bnds(:)%rlenx = 0
      bnds(:)%slenx = 0
      bnds(:)%rlen2 = 0
      bnds(:)%slen2 = 0
      bnds(:)%rlen_uv = 0
      bnds(:)%slen_uv = 0
      bnds(:)%rlenx_uv = 0
      bnds(:)%slenx_uv = 0
      bnds(:)%rlen2_uv = 0
      bnds(:)%slen2_uv = 0

!     In the first pass through, set up list of points to be requested from
!     other processors. These points are placed in the "iextra" region at the
!     end of arrays. The indirect indices are updated to point into this 
!     region.
      iext = 0
 
      rcolsp(:)%ihbg(1) = 1
      rcolsp(:)%ihfn(1) = 0
      
      do icol=1,maxcolour
      
         do n=1,npan

            !     Start with N edge
            j=jpan
            do i=1,ipan
               iq = indg(i,j,n)
               iqq = in_g(iq)
               ! Which processor has this point
               rproc = qproc(iqq)
               if ( rproc == myid ) cycle ! Don't add points already on this proc.
               mycol = findcolour(iqq)
               if ( mycol /= icol ) cycle
               iql = indp(i,j,n)  !  Local index
               ! Add this point to request list
               call check_bnds_alloc(rproc, iext)
               bnds(rproc)%rlenh = bnds(rproc)%rlenh + 1
               rcolsp(rproc)%ihfn(icol) = rcolsp(rproc)%ihfn(icol) + 1
               bnds(rproc)%request_list(bnds(rproc)%rlenh) = iqq
               ! Increment extended region index
               iext = iext + 1
               bnds(rproc)%unpack_list(bnds(rproc)%rlenh) = iext
               in(iql) = ifull+iext
            end do

            !     E edge
            i = ipan
            do j=1,jpan
               iq = indg(i,j,n)
               iqq = ie_g(iq)
               ! Which processor has this point
               rproc = qproc(iqq)
               if ( rproc == myid ) cycle ! Don't add points already on this proc.
               mycol = findcolour(iqq)
               if ( mycol /= icol ) cycle
               iql = indp(i,j,n)  !  Local index
               ! Add this point to request list
               call check_bnds_alloc(rproc, iext)
               bnds(rproc)%rlenh = bnds(rproc)%rlenh + 1
               rcolsp(rproc)%ihfn(icol) = rcolsp(rproc)%ihfn(icol) + 1
               bnds(rproc)%request_list(bnds(rproc)%rlenh) = iqq
               ! Increment extended region index
               iext = iext + 1
               bnds(rproc)%unpack_list(bnds(rproc)%rlenh) = iext
               ie(iql) = ifull+iext
            end do
         end do ! n=1,npan
         
         if ( icol < 3 ) then
            rcolsp(:)%ihbg(icol+1) = rcolsp(:)%ihfn(icol) + 1
            rcolsp(:)%ihfn(icol+1) = rcolsp(:)%ihfn(icol)
         end if
      
      end do ! icol=1,maxcolour
      
      bnds(:)%rlen = bnds(:)%rlenh  ! so that they are appended
      rcolsp(:)%ifbg(1) = rcolsp(:)%ihfn(3) + 1
      rcolsp(:)%iffn(1) = rcolsp(:)%ihfn(3)
      
      do icol=1,maxcolour
      
         do n=1,npan

            !     W edge
            i = 1
            do j=1,jpan
               ! 1D code takes care of corners separately at the end so only goes
               ! over 1,jpan here.
               iq = indg(i,j,n)
               iqq = iw_g(iq) ! iqq is the global index of the required neighbouring point.
               ! Which processor has this point
               rproc = qproc(iqq)
               if ( rproc == myid ) cycle ! Don't add points already on this proc.
               mycol = findcolour(iqq)
               if ( mycol /= icol ) cycle
               iql = indp(i,j,n)  !  Local index
               ! Add this point to request list
               call check_bnds_alloc(rproc, iext)
               bnds(rproc)%rlen = bnds(rproc)%rlen + 1
               rcolsp(rproc)%iffn(icol) = rcolsp(rproc)%iffn(icol) + 1
               bnds(rproc)%request_list(bnds(rproc)%rlen) = iqq
               ! Increment extended region index
               iext = iext + 1
               bnds(rproc)%unpack_list(bnds(rproc)%rlen) = iext
               iw(iql) = ifull+iext
            end do

            !     S edge
            j=1
            do i=1,ipan
               iq = indg(i,j,n)
               iqq = is_g(iq)
               ! Which processor has this point
               rproc = qproc(iqq)
               if ( rproc == myid ) cycle ! Don't add points already on this proc.
               mycol = findcolour(iqq)
               if ( mycol /= icol ) cycle
               iql = indp(i,j,n)  !  Local index
               ! Add this point to request list
               call check_bnds_alloc(rproc, iext)
               bnds(rproc)%rlen = bnds(rproc)%rlen + 1
               rcolsp(rproc)%iffn(icol) = rcolsp(rproc)%iffn(icol) + 1
               bnds(rproc)%request_list(bnds(rproc)%rlen) = iqq
               ! Increment extended region index
               iext = iext + 1
               bnds(rproc)%unpack_list(bnds(rproc)%rlen) = iext
               is(iql) = ifull+iext
            end do
         end do ! n=1,npan

         if ( icol < 3 ) then
            rcolsp(:)%ifbg(icol+1) = rcolsp(:)%iffn(icol) + 1
            rcolsp(:)%iffn(icol+1) = rcolsp(:)%iffn(icol)
         end if
         
      end do ! icol=1,maxcolour

      bnds(:)%rlenx = bnds(:)%rlen  ! so that they're appended.
      
!     Now handle the special corner values that need to be remapped
!     This adds to rlen, so needs to come before the _XX stuff.
      do n=1,npan
         ! NE, EN
         iq = indp(ipan,jpan,n)
         iqg = indg(ipan,jpan,n)
         iqq = ine_g(iqg)
         ! Which processor has this point
         rproc = qproc(iqq)
         if ( rproc /= myid ) then ! Add to list
            call check_bnds_alloc(rproc, iext)
            bnds(rproc)%rlenx = bnds(rproc)%rlenx + 1
            bnds(rproc)%request_list(bnds(rproc)%rlenx) = iqq
            ! Increment extended region index
            iext = iext + 1
            bnds(rproc)%unpack_list(bnds(rproc)%rlenx) = iext
            ine(iq) = ifull+iext
         end if

         if (ien_g(iqg)==ine_g(iqg)) then
            ien(iq)=ine(iq)
         else
            iqq = ien_g(iqg)
            ! Which processor has this point
            rproc = qproc(iqq)
            if ( rproc /= myid ) then ! Add to list
               call check_bnds_alloc(rproc, iext)
               bnds(rproc)%rlenx = bnds(rproc)%rlenx + 1
               bnds(rproc)%request_list(bnds(rproc)%rlenx) = iqq
               ! Increment extended region index
               iext = iext + 1
               bnds(rproc)%unpack_list(bnds(rproc)%rlenx) = iext
               ien(iq) = ifull+iext
            end if
         end if

         ! SE, ES
         iq = indp(ipan,1,n)
         iqg = indg(ipan,1,n)
         iqq = ise_g(iqg)
         ! Which processor has this point
         rproc = qproc(iqq)
         if ( rproc /= myid ) then ! Add to list
            call check_bnds_alloc(rproc, iext)
            bnds(rproc)%rlenx = bnds(rproc)%rlenx + 1
            bnds(rproc)%request_list(bnds(rproc)%rlenx) = iqq
            ! Increment extended region index
            iext = iext + 1
            bnds(rproc)%unpack_list(bnds(rproc)%rlenx) = iext
            ise(iq) = ifull+iext
         end if

         if (ies_g(iqg)==ise_g(iqg)) then
            ies(iq)=ise(iq)
         else
            iqq = ies_g(iqg)
            ! Which processor has this point
            rproc = qproc(iqq)
            if ( rproc /= myid ) then ! Add to list
               call check_bnds_alloc(rproc, iext)
               bnds(rproc)%rlenx = bnds(rproc)%rlenx + 1
               bnds(rproc)%request_list(bnds(rproc)%rlenx) = iqq
               ! Increment extended region index
               iext = iext + 1
               bnds(rproc)%unpack_list(bnds(rproc)%rlenx) = iext
               ies(iq) = ifull+iext
            end if
         end if

         ! WN, NW
         iq = indp(1,jpan,n)
         iqg = indg(1,jpan,n)
         iqq = iwn_g(iqg)
         ! Which processor has this point
         rproc = qproc(iqq)
         if ( rproc /= myid ) then ! Add to list
            call check_bnds_alloc(rproc, iext)
            bnds(rproc)%rlenx = bnds(rproc)%rlenx + 1
            bnds(rproc)%request_list(bnds(rproc)%rlenx) = iqq
            ! Increment extended region index
            iext = iext + 1
            bnds(rproc)%unpack_list(bnds(rproc)%rlenx) = iext
            iwn(iq) = ifull+iext
         end if

         if (inw_g(iqg)==iwn_g(iqg)) then
            inw(iq)=iwn(iq)
         else
            iqq = inw_g(iqg)
            ! Which processor has this point
            rproc = qproc(iqq)
            if ( rproc /= myid ) then ! Add to list
               call check_bnds_alloc(rproc, iext)
               bnds(rproc)%rlenx = bnds(rproc)%rlenx + 1
               bnds(rproc)%request_list(bnds(rproc)%rlenx) = iqq
               ! Increment extended region index
               iext = iext + 1
               bnds(rproc)%unpack_list(bnds(rproc)%rlenx) = iext
               inw(iq) = ifull+iext
            end if
         end if

         ! SW, WS
         iq = indp(1,1,n)
         iqg = indg(1,1,n)
         iqq = isw_g(iqg)
         ! Which processor has this point
         rproc = qproc(iqq)
         if ( rproc /= myid ) then ! Add to list
            call check_bnds_alloc(rproc, iext)
            bnds(rproc)%rlenx = bnds(rproc)%rlenx + 1
            bnds(rproc)%request_list(bnds(rproc)%rlenx) = iqq
            ! Increment extended region index
            iext = iext + 1
            bnds(rproc)%unpack_list(bnds(rproc)%rlenx) = iext
            isw(iq) = ifull+iext
         end if

         if (iws_g(iqg)==isw_g(iqg)) then
            iws(iq)=isw(iq)
         else
            iqq = iws_g(iqg)
            ! Which processor has this point
            rproc = qproc(iqq)
            if ( rproc /= myid ) then ! Add to list
               call check_bnds_alloc(rproc, iext)
               bnds(rproc)%rlenx = bnds(rproc)%rlenx + 1
               bnds(rproc)%request_list(bnds(rproc)%rlenx) = iqq
               ! Increment extended region index
               iext = iext + 1
               bnds(rproc)%unpack_list(bnds(rproc)%rlenx) = iext
               iws(iq) = ifull+iext
            end if
         end if

      end do

!     Now set up the second row
      bnds(:)%rlen2 = bnds(:)%rlenx  ! so that they're appended.
      
      do n=1,npan

         !     Start with W edge
         i = 1
         do j=1,jpan
            iq = indg(i,j,n)
            iqq = iww_g(iq)

            ! Which processor has this point
            rproc = qproc(iqq)
            if ( rproc == myid ) cycle ! Don't add points already on this proc.
            ! Add this point to request list
            call check_bnds_alloc(rproc, iext)
            bnds(rproc)%rlen2 = bnds(rproc)%rlen2 + 1
            bnds(rproc)%request_list(bnds(rproc)%rlen2) = iqq
            ! Increment extended region index
            iext = iext + 1
            bnds(rproc)%unpack_list(bnds(rproc)%rlen2) = iext
            iql = indp(i,j,n)  !  Local index
            iww(iql) = ifull+iext
         end do

         !     N edge
         j=jpan
         do i=1,ipan
            iq = indg(i,j,n)
            iqq = inn_g(iq)

            ! Which processor has this point
            rproc = qproc(iqq)
            if ( rproc == myid ) cycle ! Don't add points already on this proc.
            ! Add this point to request list
            call check_bnds_alloc(rproc, iext)
            bnds(rproc)%rlen2 = bnds(rproc)%rlen2 + 1
            bnds(rproc)%request_list(bnds(rproc)%rlen2) = iqq
            ! Increment extended region index
            iext = iext + 1
            bnds(rproc)%unpack_list(bnds(rproc)%rlen2) = iext
            iql = indp(i,j,n)  !  Local index
            inn(iql) = ifull+iext
         end do

         !     E edge
         i = ipan
         do j=1,jpan
            iq = indg(i,j,n)
            iqq = iee_g(iq)

            ! Which processor has this point
            rproc = qproc(iqq)
            if ( rproc == myid ) cycle ! Don't add points already on this proc.
            ! Add this point to request list
            call check_bnds_alloc(rproc, iext)
            bnds(rproc)%rlen2 = bnds(rproc)%rlen2 + 1
            bnds(rproc)%request_list(bnds(rproc)%rlen2) = iqq
            ! Increment extended region index
            iext = iext + 1
            bnds(rproc)%unpack_list(bnds(rproc)%rlen2) = iext
            iql = indp(i,j,n)  !  Local index
            iee(iql) = ifull+iext
         end do

         !     S edge
         j=1
         do i=1,ipan
            iq = indg(i,j,n)
            iqq = iss_g(iq)

            ! Which processor has this point
            rproc = qproc(iqq)
            if ( rproc == myid ) cycle ! Don't add points already on this proc.
            ! Add this point to request list
            call check_bnds_alloc(rproc, iext)
            bnds(rproc)%rlen2 = bnds(rproc)%rlen2 + 1
            bnds(rproc)%request_list(bnds(rproc)%rlen2) = iqq
            ! Increment extended region index
            iext = iext + 1
            bnds(rproc)%unpack_list(bnds(rproc)%rlen2) = iext
            iql = indp(i,j,n)  !  Local index
            iss(iql) = ifull+iext
         end do
      end do ! n=1,npan

!     Now handle the second order special corner values
      do n=1,npan
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
            ! Use is first because it'll be on same face.
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
            ! Use is first because it'll be on same face.
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

!     Indices that are missed above (should be a better way to get these)
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
         call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
      end if

      neighnum = count( bnds(:)%rlen2 > 0 )
      ! ireq needs 1 point for the MPI_Waitall which can use ireq(rreq+1)
      allocate( ireq(max(2*neighnum,1)) )
      allocate( rlist(neighnum) )
      allocate( status(MPI_STATUS_SIZE,2*neighnum ))
      allocate( dums(7,neighnum),dumr(7,neighnum) )
      allocate( dumsb(9,neighnum),dumrb(9,neighnum) )
      allocate( dumsl(maxbuflen,neighnum),dumrl(maxbuflen,neighnum) )

      ! set up neighbour lists
      allocate ( neighlist(neighnum) )
      allocate ( neighmap(0:nproc-1) )
      ncount = 0
      neighmap = 0 ! missing
      do iproc = 1,nproc-1
         rproc = modulo(myid+iproc,nproc)
         if ( bnds(rproc)%rlen2 > 0 ) then
            ncount = ncount + 1
            neighlist(ncount) = rproc
            neighmap(rproc) = ncount
         end if
      end do
      
      if ( ncount /= neighnum ) then
         write(6,*) "ERROR: neighnum mismatch"
         write(6,*) "neighnum, ncount ",neighnum, ncount
         call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
      end if
   
      
!     Now, for each processor send the list of points I want.
!     The state of being a neighbour is reflexive so only expect to
!     recv from those processors I send to (are there grid arrangements for
!     which this would not be true?)
!     Get the complete request lists by using rlen2
      nreq = 0
      do iproc = 1,neighnum
         rproc = neighlist(iproc)  ! Recv from
         nreq = nreq + 1
         ! Use the maximum size in the recv call.
         llen = bnds(rproc)%len
         lproc = rproc
         call MPI_IRecv( bnds(rproc)%send_list(1), llen, &
              ltype, lproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
      end do
      do iproc = 1,neighnum
         sproc = neighlist(iproc)  ! Send to
         ! Send list of requests
         nreq = nreq + 1
         llen = bnds(sproc)%rlen2
         lproc = sproc
         call MPI_ISend( bnds(sproc)%request_list(1), llen, &
              ltype, lproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
      end do      
      call MPI_Waitall( nreq, ireq, status, ierr )

!     Now get the actual sizes from the status
      nreq = 0
      do iproc = 1,neighnum
         rproc = neighlist(iproc)  ! Recv from
         ! First half of nreq are recv
         nreq = nreq + 1
         call MPI_Get_count( status(:,nreq), ltype, lcount, ierr )
         ! This the number of points I have to send to rproc.
         bnds(rproc)%slen2 = lcount
      end do

      
!     For rlen and rlen2, just communicate the lengths. The indices have 
!     already been taken care of.
      scolsp(:)%ihfn(1) = 0
      scolsp(:)%ihfn(2) = 0
      scolsp(:)%ihfn(3) = 0
      scolsp(:)%iffn(1) = 0
      scolsp(:)%iffn(2) = 0
      scolsp(:)%iffn(3) = 0
      nreq = 0
      do iproc = 1,neighnum
         rproc = neighlist(iproc) ! Recv from
         nreq = nreq + 1
         lproc = rproc
         call MPI_IRecv( dumrb(:,iproc), 9_4, ltype, lproc, &
              itag, MPI_COMM_WORLD, ireq(nreq), ierr )
      end do
      do iproc = neighnum,1,-1
         ! Send and recv from same proc
         sproc = neighlist(iproc)  ! Send to
         nreq = nreq + 1
         dumsb(1,iproc) = bnds(sproc)%rlenh
         dumsb(2,iproc) = bnds(sproc)%rlen
         dumsb(3,iproc) = bnds(sproc)%rlenx
         dumsb(4,iproc) = rcolsp(sproc)%ihfn(1)
         dumsb(5,iproc) = rcolsp(sproc)%ihfn(2)
         dumsb(6,iproc) = rcolsp(sproc)%ihfn(3)
         dumsb(7,iproc) = rcolsp(sproc)%iffn(1)
         dumsb(8,iproc) = rcolsp(sproc)%iffn(2)
         dumsb(9,iproc) = rcolsp(sproc)%iffn(3)
         lproc = sproc
         call MPI_ISend( dumsb(:,iproc), 9_4, ltype, lproc, &
              itag, MPI_COMM_WORLD, ireq(nreq), ierr )
      end do
      call MPI_Waitall(nreq,ireq,MPI_STATUSES_IGNORE,ierr)
      do iproc = 1,neighnum
         rproc = neighlist(iproc)
         bnds(rproc)%slenh     = dumrb(1,iproc)
         bnds(rproc)%slen      = dumrb(2,iproc)
         bnds(rproc)%slenx     = dumrb(3,iproc)
         scolsp(rproc)%ihfn(1) = dumrb(4,iproc)
         scolsp(rproc)%ihfn(2) = dumrb(5,iproc)
         scolsp(rproc)%ihfn(3) = dumrb(6,iproc)
         scolsp(rproc)%iffn(1) = dumrb(7,iproc)
         scolsp(rproc)%iffn(2) = dumrb(8,iproc)
         scolsp(rproc)%iffn(3) = dumrb(9,iproc)
      end do
      scolsp(:)%ihbg(1) = 1
      scolsp(:)%ihbg(2) = scolsp(:)%ihfn(1) + 1
      scolsp(:)%ihbg(3) = scolsp(:)%ihfn(2) + 1
      scolsp(:)%ifbg(1) = scolsp(:)%ihfn(3) + 1
      scolsp(:)%ifbg(2) = scolsp(:)%iffn(1) + 1
      scolsp(:)%ifbg(3) = scolsp(:)%iffn(2) + 1

!     Start of UV section

!     In the first pass through, set up list of points to be requested from
!     other processors. In the 1D code values on the same processor are
!     copied only if they have to be swapped.
!     This only makes a difference on 1, 2 or 3 processors.

      iextu = 0
      iextv = 0

      ! save start of isv indices
      rsplit(:)%isvbg = 1
      rsplit(:)%iwufn = 0

      !     S edge, V
      j=1
      do n=1,npan
         do i=1,ipan
            iqg = indg(i,j,n)
            iqq = is_g(iqg)
            rproc = qproc(iqq)
            swap = edge_s .and. swap_s(n-noff)
            if ( rproc == myid .and. .not. swap ) cycle
            call check_bnds_alloc(rproc, iextv)
            bnds(rproc)%rlen_uv = bnds(rproc)%rlen_uv + 1
            bnds(rproc)%request_list_uv(bnds(rproc)%rlen_uv) = -iqq
            rsplit(rproc)%iwufn = rsplit(rproc)%iwufn + 1
            ! Increment extended region index
            iextv = iextv + 1
            bnds(rproc)%unpack_list_uv(bnds(rproc)%rlen_uv) = -iextv
            iql = indp(i,j,n)  !  Local index
            isv(iql) = ifull+iextv
            bnds(rproc)%uv_swap(bnds(rproc)%rlen_uv) = swap
            bnds(rproc)%uv_neg(bnds(rproc)%rlen_uv) = .false.
         end do
      end do
         
      !     Start with W edge, U values
      i = 1
      do n=1,npan
         do j=1,jpan
            iqg = indg(i,j,n)
            iqq = iw_g(iqg)
            ! Which processor has this point
            rproc = qproc(iqq)
            ! Only need to add to bounds region if it's on another processor
            ! or if it's on this processor and needs to be swapped.
            ! Decide if u/v need to be swapped. My face is n-noff
            swap = edge_w .and. swap_w(n-noff)
            if ( rproc == myid .and. .not. swap ) cycle
            call check_bnds_alloc(rproc, iextu)
            bnds(rproc)%rlen_uv = bnds(rproc)%rlen_uv + 1
            bnds(rproc)%request_list_uv(bnds(rproc)%rlen_uv) = iqq
            rsplit(rproc)%iwufn = rsplit(rproc)%iwufn + 1
            ! Increment extended region index
            iextu = iextu + 1
            bnds(rproc)%unpack_list_uv(bnds(rproc)%rlen_uv) = iextu
            iql = indp(i,j,n)  !  Local index
            iwu(iql) = ifull+iextu
            bnds(rproc)%uv_swap(bnds(rproc)%rlen_uv) = swap
            bnds(rproc)%uv_neg(bnds(rproc)%rlen_uv) = .false.
         end do
      end do

      ! save start of inv indices
      rsplit(:)%invbg = rsplit(:)%iwufn + 1
      rsplit(:)%ieufn = rsplit(:)%iwufn

      !     N edge (V)
      j=jpan
      do n=1,npan
         do i=1,ipan
            iqg = indg(i,j,n)
            iqq = in_g(iqg)
            rproc = qproc(iqq)
            swap = edge_n .and. swap_n(n-noff)
            if ( rproc == myid .and. .not. swap ) cycle
            ! Add this point to request list
            call check_bnds_alloc(rproc, iextv)
            bnds(rproc)%rlen_uv = bnds(rproc)%rlen_uv + 1
            ! to show that this is v rather than u, flip sign
            bnds(rproc)%request_list_uv(bnds(rproc)%rlen_uv) = -iqq
            rsplit(rproc)%ieufn = rsplit(rproc)%ieufn + 1
            ! Increment extended region index
            iextv = iextv + 1
            bnds(rproc)%unpack_list_uv(bnds(rproc)%rlen_uv) = -iextv
            iql = indp(i,j,n)  !  Local index
            inv(iql) = ifull+iextv
            bnds(rproc)%uv_swap(bnds(rproc)%rlen_uv) = swap
            bnds(rproc)%uv_neg(bnds(rproc)%rlen_uv) = .false.
         end do
      end do

      !     E edge, U
      i = ipan
      do n=1,npan
         do j=1,jpan
            iqg = indg(i,j,n)
            iqq = ie_g(iqg)
            rproc = qproc(iqq)
            swap = edge_e .and. swap_e(n-noff)
            if ( rproc == myid .and. .not. swap ) cycle
            ! Add this point to request list
            call check_bnds_alloc(rproc, iextu)
            bnds(rproc)%rlen_uv = bnds(rproc)%rlen_uv + 1
            bnds(rproc)%request_list_uv(bnds(rproc)%rlen_uv) = iqq
            rsplit(rproc)%ieufn = rsplit(rproc)%ieufn + 1
            ! Increment extended region index
            iextu = iextu + 1
            bnds(rproc)%unpack_list_uv(bnds(rproc)%rlen_uv) = iextu
            iql = indp(i,j,n)  !  Local index
            ieu(iql) = ifull+iextu
            bnds(rproc)%uv_swap(bnds(rproc)%rlen_uv) = swap
            bnds(rproc)%uv_neg(bnds(rproc)%rlen_uv) = .false.
         end do
      end do

!     Second pass
      bnds(:)%rlen2_uv = bnds(:)%rlen_uv
      bnds(:)%slen2_uv = bnds(:)%slen_uv
      ieeu = iee
      iwwu = iww
      innv = inn
      issv = iss

      ! save start of issv indices
      rsplit(:)%issvbg = rsplit(:)%ieufn + 1
      rsplit(:)%iwwufn = rsplit(:)%ieufn

      !     S edge, V
      j=1
      do n=1,npan
         do i=1,ipan
            iqg = indg(i,j,n)
            iqq = iss_g(iqg)
            rproc = qproc(iqq)
            swap = edge_s .and. swap_s(n-noff)
            if ( rproc == myid .and. .not. swap ) cycle
            call check_bnds_alloc(rproc, iextv)
            bnds(rproc)%rlen2_uv = bnds(rproc)%rlen2_uv + 1
            bnds(rproc)%request_list_uv(bnds(rproc)%rlen2_uv) = -iqq
            rsplit(rproc)%iwwufn = rsplit(rproc)%iwwufn + 1
            ! Increment extended region index
            iextv = iextv + 1
            bnds(rproc)%unpack_list_uv(bnds(rproc)%rlen2_uv) = -iextv
            iql = indp(i,j,n)  !  Local index
            issv(iql) = ifull+iextv
            bnds(rproc)%uv_swap(bnds(rproc)%rlen2_uv) = swap
            bnds(rproc)%uv_neg(bnds(rproc)%rlen2_uv) = .false.
         end do
      end do ! n=1,npan

      !     Start with W edge, U values
      i = 1
      do n=1,npan
         do j=1,jpan
            iqg = indg(i,j,n)
            iqq = iww_g(iqg)
            ! Which processor has this point
            rproc = qproc(iqq)
            swap = edge_w .and. swap_w(n-noff)
            if ( rproc == myid .and. .not. swap ) cycle
            ! Add this point to request list
            call check_bnds_alloc(rproc, iextu)
            bnds(rproc)%rlen2_uv = bnds(rproc)%rlen2_uv + 1
            bnds(rproc)%request_list_uv(bnds(rproc)%rlen2_uv) = iqq
            rsplit(rproc)%iwwufn = rsplit(rproc)%iwwufn + 1
            ! Increment extended region index
            iextu = iextu + 1
            bnds(rproc)%unpack_list_uv(bnds(rproc)%rlen2_uv) = iextu
            iql = indp(i,j,n)  !  Local index
            iwwu(iql) = ifull+iextu
            ! Decide if u/v need to be swapped. My face is n-noff
            bnds(rproc)%uv_swap(bnds(rproc)%rlen2_uv) = swap
            bnds(rproc)%uv_neg(bnds(rproc)%rlen2_uv) = .false.
         end do
      end do

      ! save start of innv indices
      rsplit(:)%innvbg = rsplit(:)%iwwufn + 1
      rsplit(:)%ieeufn = rsplit(:)%iwwufn

      !     N edge (V)
      j=jpan
      do n=1,npan
         do i=1,ipan
            iqg = indg(i,j,n)
            iqq = inn_g(iqg)
            rproc = qproc(iqq)
            swap = edge_n .and. swap_n(n-noff)
            if ( rproc == myid .and. .not. swap ) cycle
            ! Add this point to request list
            call check_bnds_alloc(rproc, iextv)
            bnds(rproc)%rlen2_uv = bnds(rproc)%rlen2_uv + 1
            ! to show that this is v rather than u, flip sign
            bnds(rproc)%request_list_uv(bnds(rproc)%rlen2_uv) = -iqq
            rsplit(rproc)%ieeufn = rsplit(rproc)%ieeufn + 1
            ! Increment extended region index
            iextv = iextv + 1
            bnds(rproc)%unpack_list_uv(bnds(rproc)%rlen2_uv) = -iextv
            iql = indp(i,j,n)  !  Local index
            innv(iql) = ifull+iextv
            bnds(rproc)%uv_swap(bnds(rproc)%rlen2_uv) = swap
            bnds(rproc)%uv_neg(bnds(rproc)%rlen2_uv) = .false.
         end do
      end do

      !     E edge, U
      i = ipan
      do n=1,npan
         do j=1,jpan
            iqg = indg(i,j,n)
            iqq = iee_g(iqg)
            rproc = qproc(iqq)
            swap = edge_e .and. swap_e(n-noff)
            if ( rproc == myid .and. .not. swap ) cycle
            ! Add this point to request list
            call check_bnds_alloc(rproc, iextu)
            bnds(rproc)%rlen2_uv = bnds(rproc)%rlen2_uv + 1
            bnds(rproc)%request_list_uv(bnds(rproc)%rlen2_uv) = iqq
            rsplit(rproc)%ieeufn = rsplit(rproc)%ieeufn + 1
            ! Increment extended region index
            iextu = iextu + 1
            bnds(rproc)%unpack_list_uv(bnds(rproc)%rlen2_uv) = iextu
            iql = indp(i,j,n)  !  Local index
            ieeu(iql) = ifull+iextu
            bnds(rproc)%uv_swap(bnds(rproc)%rlen2_uv) = swap
            bnds(rproc)%uv_neg(bnds(rproc)%rlen2_uv) = .false.
         end do
      end do

      ! Third pass
      bnds(:)%rlenx_uv = bnds(:)%rlen2_uv
      bnds(:)%slenx_uv = bnds(:)%slen2_uv

      ! save start of isu indices
      rsplit(:)%isubg = rsplit(:)%ieeufn + 1
      rsplit(:)%ievfn = rsplit(:)%ieeufn

      !     S edge, U
      j=1
      do n=1,npan
         do i=1,ipan
            iqg = indg(i,j,n)
            iqq = is_g(iqg)
            rproc = qproc(iqq)
            swap = edge_s .and. swap_s(n-noff)
            if ( rproc == myid .and. .not. swap ) cycle
            call check_bnds_alloc(rproc, iextu)
            bnds(rproc)%rlenx_uv = bnds(rproc)%rlenx_uv + 1
            bnds(rproc)%request_list_uv(bnds(rproc)%rlenx_uv) = iqq
            rsplit(rproc)%ievfn = rsplit(rproc)%ievfn + 1
            ! Increment extended region index
            iextu = iextu + 1
            bnds(rproc)%unpack_list_uv(bnds(rproc)%rlenx_uv) = iextu
            iql = indp(i,j,n)  !  Local index
            isu(iql) = ifull+iextu
            bnds(rproc)%uv_swap(bnds(rproc)%rlenx_uv) = swap
            bnds(rproc)%uv_neg(bnds(rproc)%rlenx_uv) = swap
         end do
      end do ! n=1,npan

      !     Start with W edge, V values
      i = 1
      do n=1,npan
         do j=1,jpan
            iqg = indg(i,j,n)
            iqq = iw_g(iqg)
            ! Which processor has this point
            rproc = qproc(iqq)
            ! Only need to add to bounds region if it's on another processor
            ! or if it's on this processor and needs to be swapped.
            ! Decide if u/v need to be swapped. My face is n-noff
            swap = edge_w .and. swap_w(n-noff)
            if ( rproc == myid .and. .not. swap ) cycle
            call check_bnds_alloc(rproc, iextv)
            bnds(rproc)%rlenx_uv = bnds(rproc)%rlenx_uv + 1
            ! to show that this is v rather than u, flip sign
            bnds(rproc)%request_list_uv(bnds(rproc)%rlenx_uv) = -iqq
            rsplit(rproc)%ievfn = rsplit(rproc)%ievfn + 1
            ! Increment extended region index
            iextv = iextv + 1
            bnds(rproc)%unpack_list_uv(bnds(rproc)%rlenx_uv) = -iextv
            iql = indp(i,j,n)  !  Local index
            iwv(iql) = ifull+iextv
            bnds(rproc)%uv_swap(bnds(rproc)%rlenx_uv) = swap
            bnds(rproc)%uv_neg(bnds(rproc)%rlenx_uv) = swap
         end do
      end do

      !     N edge (U)
      j=jpan
      do n=1,npan
         do i=1,ipan
            iqg = indg(i,j,n)
            iqq = in_g(iqg)
            rproc = qproc(iqq)
            swap = edge_n .and. swap_n(n-noff)
            if ( rproc == myid .and. .not. swap ) cycle
            ! Add this point to request list
            call check_bnds_alloc(rproc, iextu)
            bnds(rproc)%rlenx_uv = bnds(rproc)%rlenx_uv + 1
            bnds(rproc)%request_list_uv(bnds(rproc)%rlenx_uv) = iqq
            rsplit(rproc)%ievfn = rsplit(rproc)%ievfn + 1
            ! Increment extended region index
            iextu = iextu + 1
            bnds(rproc)%unpack_list_uv(bnds(rproc)%rlenx_uv) = iextu
            iql = indp(i,j,n)  !  Local index
            inu(iql) = ifull+iextu
            bnds(rproc)%uv_swap(bnds(rproc)%rlenx_uv) = swap
            bnds(rproc)%uv_neg(bnds(rproc)%rlenx_uv) = swap
         end do
      end do

      !     E edge, V
      i = ipan
      do n=1,npan
         do j=1,jpan
            iqg = indg(i,j,n)
            iqq = ie_g(iqg)
            rproc = qproc(iqq)
            swap = edge_e .and. swap_e(n-noff)
            if ( rproc == myid .and. .not. swap ) cycle
            ! Add this point to request list
            call check_bnds_alloc(rproc, iextv)
            bnds(rproc)%rlenx_uv = bnds(rproc)%rlenx_uv + 1
            ! to show that this is v rather than u, flip sign
            bnds(rproc)%request_list_uv(bnds(rproc)%rlenx_uv) = -iqq
            rsplit(rproc)%ievfn = rsplit(rproc)%ievfn + 1
            ! Increment extended region index
            iextv = iextv + 1
            bnds(rproc)%unpack_list_uv(bnds(rproc)%rlenx_uv) = -iextv
            iql = indp(i,j,n)  !  Local index
            iev(iql) = ifull+iextv
            bnds(rproc)%uv_swap(bnds(rproc)%rlenx_uv) = swap
            bnds(rproc)%uv_neg(bnds(rproc)%rlenx_uv) = swap
         end do
      end do

      if ( iextu > iextra ) then
         write(6,*) "IEXTU too large", iextu, iextra
         call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
      end if

      if ( iextv > iextra ) then
         write(6,*) "IEXTV too large", iextv, iextra
         call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
      end if

!     Now, for each processor send the list of points I want.
!     Send all possible pairs and don't assume anything about the symmetry.
!     For completeness, send zero length messages too.
!     Get the length from the message status
!     Also have to send the swap list

      nreq = 0
      do iproc = 1,neighnum
         rproc = neighlist(iproc)
         if ( bnds(rproc)%rlenx_uv > 0 ) then
            nreq = nreq + 1
            ! Use the maximum size in the recv call.
            llen = bnds(rproc)%len
            lproc = rproc
            call MPI_IRecv( bnds(rproc)%send_list_uv(1), llen, ltype, lproc, &
                  itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do
      do iproc = neighnum,1,-1
         sproc = neighlist(iproc)
         if ( bnds(sproc)%rlenx_uv > 0 ) then
            ! Send list of requests
            nreq = nreq + 1
            llen = bnds(sproc)%rlenx_uv
            lproc = sproc
            call MPI_ISend( bnds(sproc)%request_list_uv(1), llen, ltype, lproc, &
                 itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do
      call MPI_Waitall( nreq, ireq, status, ierr )

!     Now get the actual sizes from the status
      nreq = 0
      do iproc = 1,neighnum
         rproc = neighlist(iproc)
         if ( bnds(rproc)%rlenx_uv > 0 ) then
            ! First half of nreq are recv
            nreq = nreq + 1
            call MPI_Get_count( status(:,nreq), ltype, lcount, ierr )
            ! This the number of points I have to send to rproc.
            bnds(rproc)%slenx_uv = lcount
         end if
      end do

!     For rlen_uv and rlenx_uv, etc, just communicate the lengths. The indices have 
!     already been taken care of.
      ssplit(:)%iwufn  = 0
      ssplit(:)%ieufn  = 0
      ssplit(:)%iwwufn = 0
      ssplit(:)%ieeufn = 0
      ssplit(:)%ievfn  = 0
      nreq = 0
      do iproc = 1,neighnum
         rproc = neighlist(iproc)  ! Recv from
         if ( bnds(rproc)%rlenx_uv > 0 ) then
            nreq = nreq + 1
            lproc = rproc
            call MPI_IRecv( dumr(:,iproc), 7_4, ltype, lproc, &
                 itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do
      do iproc = neighnum,1,-1
         sproc = neighlist(iproc)  ! Send to
         if ( bnds(sproc)%rlenx_uv > 0 ) then
            nreq = nreq + 1
            dums(1,iproc) = bnds(sproc)%rlen_uv
            dums(2,iproc) = bnds(sproc)%rlen2_uv
            dums(3,iproc) = rsplit(sproc)%iwufn
            dums(4,iproc) = rsplit(sproc)%ieufn
            dums(5,iproc) = rsplit(sproc)%iwwufn
            dums(6,iproc) = rsplit(sproc)%ieeufn
            dums(7,iproc) = rsplit(sproc)%ievfn
            lproc = sproc
            call MPI_ISend( dums(:,iproc), 7_4, ltype, lproc, &
                 itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do
      call MPI_Waitall(nreq,ireq,MPI_STATUSES_IGNORE,ierr)
      do iproc = 1,neighnum
         rproc = neighlist(iproc)
         if ( bnds(rproc)%rlenx_uv > 0 ) then
            bnds(rproc)%slen_uv  = dumr(1,iproc)
            bnds(rproc)%slen2_uv = dumr(2,iproc)
            ssplit(rproc)%iwufn  = dumr(3,iproc)
            ssplit(rproc)%ieufn  = dumr(4,iproc)
            ssplit(rproc)%iwwufn = dumr(5,iproc)
            ssplit(rproc)%ieeufn = dumr(6,iproc)
            ssplit(rproc)%ievfn  = dumr(7,iproc)
         end if
      end do
      ssplit(:)%isvbg  = 1
      ssplit(:)%invbg  = ssplit(:)%iwufn  + 1
      ssplit(:)%issvbg = ssplit(:)%ieufn  + 1
      ssplit(:)%innvbg = ssplit(:)%iwwufn + 1
      ssplit(:)%isubg  = ssplit(:)%ieeufn + 1

      ! Only send the swap list once
      nreq = 0
      do iproc = 1,neighnum
         rproc = neighlist(iproc) ! Recv from
         if ( bnds(rproc)%rlenx_uv > 0 ) then
            nreq = nreq + 1
            llen = bnds(rproc)%slenx_uv
            lproc = rproc
            call MPI_IRecv( dumrl(:,iproc), llen, MPI_LOGICAL, &
                  lproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do
      do iproc = neighnum,1,-1
         sproc = neighlist(iproc) ! Send to
         if ( bnds(sproc)%rlenx_uv > 0 ) then
            ! Send list of requests
            nreq = nreq + 1
            llen = bnds(sproc)%rlenx_uv
            lproc = sproc
            dumsl(1:bnds(sproc)%rlenx_uv,iproc) = bnds(sproc)%uv_swap(1:bnds(sproc)%rlenx_uv)
            call MPI_ISend( dumsl(:,iproc), llen, MPI_LOGICAL, &
                  lproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do
      call MPI_Waitall( nreq, ireq, MPI_STATUSES_IGNORE, ierr )
      do iproc = 1,neighnum
         rproc = neighlist(iproc)
         if ( bnds(rproc)%slenx_uv > 0 ) then
            bnds(rproc)%send_swap(1:bnds(rproc)%slenx_uv) = dumrl(1:bnds(rproc)%slenx_uv,iproc)
         end if
      end do

      ! Only send the neg list once
      nreq = 0
      do iproc = 1,neighnum
         rproc = neighlist(iproc)
         if ( bnds(rproc)%rlenx_uv > 0 ) then
            nreq = nreq + 1
            llen = bnds(rproc)%slenx_uv
            lproc = rproc
            call MPI_IRecv( dumrl(:,iproc), llen, MPI_LOGICAL, &
                  lproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do
      do iproc = neighnum,1,-1
         sproc = neighlist(iproc)
         if ( bnds(sproc)%rlenx_uv > 0 ) then
            ! Send list of requests
            nreq = nreq + 1
            llen = bnds(sproc)%rlenx_uv
            lproc = sproc
            dumsl(1:bnds(sproc)%rlenx_uv,iproc) = bnds(sproc)%uv_neg(1:bnds(sproc)%rlenx_uv)
            call MPI_ISend( dumsl(:,iproc), llen, MPI_LOGICAL,&
                  lproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do
      call MPI_Waitall( nreq, ireq, MPI_STATUSES_IGNORE, ierr )
      do iproc=1,neighnum
         rproc=neighlist(iproc)
         if ( bnds(rproc)%slenx_uv > 0 ) then
           where ( dumrl(1:bnds(rproc)%slenx_uv,iproc) )
              bnds(rproc)%send_neg(1:bnds(rproc)%slenx_uv) = -1.
           elsewhere
              bnds(rproc)%send_neg(1:bnds(rproc)%slenx_uv) =  1.
           end where
         end if
      end do
      
      ! Flag that all messages have been cleared
      nreq = 0
      rreq = 0

!     Indices that are missed above (should be a better way to get these)
      do n=1,npan
         do j=1,jpan
            iwwu(indp(2,j,n)) = iwu(indp(1,j,n))
            ieeu(indp(ipan-1,j,n)) = ieu(indp(ipan,j,n))
         end do
         do i=1,ipan
            issv(indp(i,2,n)) = isv(indp(i,1,n))
            innv(indp(i,jpan-1,n)) = inv(indp(i,jpan,n))
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
         do iq = 1,bnds(sproc)%slenx_uv
            ! send_list(iq) is global point index, i, j, n are local
            ! Use abs because sign is used as u/v flag
            iqq = abs(bnds(sproc)%send_list_uv(iq))
            call indv_mpi(iqq,i,j,n)
            bnds(sproc)%send_list_uv(iq) = sign(indp(i,j,n),bnds(sproc)%send_list_uv(iq))
         end do
      end do
      do iq = 1,bnds(myid)%rlen2
         iqq = bnds(myid)%request_list(iq)
         call indv_mpi(iqq,i,j,n)
         bnds(myid)%request_list(iq) = indp(i,j,n)
      end do
      do iq = 1,bnds(myid)%rlenx_uv
         iqq = abs(bnds(myid)%request_list_uv(iq))
         call indv_mpi(iqq,i,j,n)
         bnds(myid)%request_list_uv(iq) = sign(indp(i,j,n),bnds(myid)%request_list_uv(iq))
      end do


      ! resize arrays
      deallocate( status )
      deallocate( dumr, dums )
      deallocate( dumrb, dumsb )
      deallocate( dumrl, dumsl )
      call reducealloc 
      do iproc = 0,nproc-1
         bnds(iproc)%sbuflen = bnds(iproc)%len
         bnds(iproc)%rbuflen = bnds(iproc)%len
      end do
      do iproc = 1,neighnum
         rproc = neighlist(iproc)
         allocate ( bnds(rproc)%rbuf(bnds(rproc)%len) )
         allocate ( bnds(rproc)%sbuf(bnds(rproc)%len) )
         allocate ( bnds(rproc)%r8buf(bnds(rproc)%len) )
         allocate ( bnds(rproc)%s8buf(bnds(rproc)%len) )
      end do
      

!  Final check for values that haven't been set properly
      do n=1,npan
         do j=1,jpan
            do i=1,ipan
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
   
   subroutine reducealloc
      ! free memory   
      integer iproc,nlen
      integer(kind=4) ierr
      real, dimension(maxbuflen) :: rdum
      integer, dimension(maxbuflen) :: idum
      logical, dimension(maxbuflen) :: ldum
   
      do iproc = 0,nproc-1
         nlen = max(nagg*kl,3*ol)*max(bnds(iproc)%rlen2,bnds(iproc)%rlenx_uv,bnds(iproc)%slen2,bnds(iproc)%slenx_uv)
         if ( nlen < bnds(iproc)%len ) then
            bnds(iproc)%len = nlen
            if ( iproc /= myid ) then
               idum(1:bnds(iproc)%len) = bnds(iproc)%request_list(1:bnds(iproc)%len)
               deallocate ( bnds(iproc)%request_list )
               allocate ( bnds(iproc)%request_list(bnds(iproc)%len) )
               bnds(iproc)%request_list(1:bnds(iproc)%len) = idum(1:bnds(iproc)%len)
               idum(1:bnds(iproc)%len) = bnds(iproc)%send_list(1:bnds(iproc)%len)
               deallocate ( bnds(iproc)%send_list )
               allocate ( bnds(iproc)%send_list(bnds(iproc)%len) )
               bnds(iproc)%send_list(1:bnds(iproc)%len) = idum(1:bnds(iproc)%len)
               idum(1:bnds(iproc)%len) = bnds(iproc)%unpack_list(1:bnds(iproc)%len)
               deallocate ( bnds(iproc)%unpack_list )
               allocate ( bnds(iproc)%unpack_list(bnds(iproc)%len) )
               bnds(iproc)%unpack_list(1:bnds(iproc)%len) = idum(1:bnds(iproc)%len)
            end if
            idum(1:bnds(iproc)%len) = bnds(iproc)%request_list_uv(1:bnds(iproc)%len)
            deallocate ( bnds(iproc)%request_list_uv )
            allocate ( bnds(iproc)%request_list_uv(bnds(iproc)%len) )
            bnds(iproc)%request_list_uv(1:bnds(iproc)%len) = idum(1:bnds(iproc)%len)
            idum(1:bnds(iproc)%len) = bnds(iproc)%send_list_uv(1:bnds(iproc)%len)
            deallocate ( bnds(iproc)%send_list_uv )
            allocate ( bnds(iproc)%send_list_uv(bnds(iproc)%len) )
            bnds(iproc)%send_list_uv(1:bnds(iproc)%len) = idum(1:bnds(iproc)%len)
            idum(1:bnds(iproc)%len) = bnds(iproc)%unpack_list_uv(1:bnds(iproc)%len)
            deallocate ( bnds(iproc)%unpack_list_uv )
            allocate ( bnds(iproc)%unpack_list_uv(bnds(iproc)%len) )
            bnds(iproc)%unpack_list_uv(1:bnds(iproc)%len) = idum(1:bnds(iproc)%len)
            ldum(1:bnds(iproc)%len) = bnds(iproc)%uv_swap(1:bnds(iproc)%len)
            deallocate ( bnds(iproc)%uv_swap )
            allocate ( bnds(iproc)%uv_swap(bnds(iproc)%len) )
            bnds(iproc)%uv_swap(1:bnds(iproc)%len) = ldum(1:bnds(iproc)%len)
            ldum(1:bnds(iproc)%len) = bnds(iproc)%send_swap(1:bnds(iproc)%len)
            deallocate ( bnds(iproc)%send_swap )
            allocate ( bnds(iproc)%send_swap(bnds(iproc)%len) )
            bnds(iproc)%send_swap(1:bnds(iproc)%len) = ldum(1:bnds(iproc)%len)
            ldum(1:bnds(iproc)%len) = bnds(iproc)%uv_neg(1:bnds(iproc)%len)
            deallocate ( bnds(iproc)%uv_neg )
            allocate ( bnds(iproc)%uv_neg(bnds(iproc)%len) )
            bnds(iproc)%uv_neg(1:bnds(iproc)%len) = ldum(1:bnds(iproc)%len)
            rdum(1:bnds(iproc)%len) = bnds(iproc)%send_neg(1:bnds(iproc)%len)
            deallocate ( bnds(iproc)%send_neg )
            allocate ( bnds(iproc)%send_neg(bnds(iproc)%len) )
            bnds(iproc)%send_neg(1:bnds(iproc)%len) = rdum(1:bnds(iproc)%len)
         else if ( nlen > bnds(iproc)%len ) then
            write(6,*) "ERROR reducing array size"
            write(6,*) "myid,iproc,nlen,len ",myid,iproc,nlen,bnds(iproc)%len
            write(6,*) "maxbuflen ",maxbuflen
            call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
         end if
      end do
   end subroutine reducealloc

   subroutine check_set(ind,str,i,j,n,iq)
      integer, intent(in) :: ind,i,j,n,iq
      integer(kind=4) :: ierr
      character(len=*) :: str
      if ( ind == huge(1) ) then
         write(6,*) str, " not set", myid, i, j, n, iq
         call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
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
      integer :: rcount, myrlen, jproc, mproc
      integer, dimension(neighnum) :: rslen, sslen
      integer(kind=4) :: ierr, itag=1, llen, sreq, lproc
      integer(kind=4) :: ldone
      integer(kind=4), dimension(neighnum) :: donelist
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
         rslen  = bnds(neighlist)%rlen2
         sslen  = bnds(neighlist)%slen2
         myrlen = bnds(myid)%rlen2
      else if ( extra ) then
         rslen  = bnds(neighlist)%rlenx
         sslen  = bnds(neighlist)%slenx
         myrlen = bnds(myid)%rlenx
      else if ( single ) then
         rslen  = bnds(neighlist)%rlen
         sslen  = bnds(neighlist)%slen
         myrlen = bnds(myid)%rlen
      else
         rslen  = bnds(neighlist)%rlenh
         sslen  = bnds(neighlist)%slenh
         myrlen = bnds(myid)%rlenh
      end if

      call START_LOG(bounds_begin)

!     Set up the buffers to recv
      nreq = 0
      do iproc = 1,neighnum
         recv_len = rslen(iproc)
         if ( recv_len > 0 ) then
            lproc = neighlist(iproc)  ! Recv from
            nreq  = nreq + 1
            rlist(nreq) = iproc
            llen  = recv_len
            call MPI_IRecv( bnds(lproc)%rbuf(1), llen, ltype, lproc, &
                   itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do
      rreq = nreq
      !     Set up the buffers to send
      do iproc = neighnum,1,-1
         ! Build up list of points
         send_len = sslen(iproc)
         if ( send_len > 0 ) then
            lproc = neighlist(iproc)  ! Send to
            bnds(lproc)%sbuf(1:send_len) = t(bnds(lproc)%send_list(1:send_len))
            nreq  = nreq + 1
            llen  = send_len
            call MPI_ISend( bnds(lproc)%sbuf(1), llen, ltype, lproc, &
                 itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do

      ! Finally see if there are any points on my own processor that need
      ! to be fixed up. This will only be in the case when nproc < npanels.
      if ( myrlen > 0 ) then
         ! request_list is same as send_list in this case
         t(ifull+bnds(myid)%unpack_list(1:myrlen)) = t(bnds(myid)%request_list(1:myrlen))
      end if

      ! Unpack incomming messages
      rcount = rreq
      do while ( rcount > 0 )
         call START_LOG(mpiwait_begin)
         call MPI_Waitsome( rreq, ireq, ldone, donelist, MPI_STATUSES_IGNORE, ierr )
         call END_LOG(mpiwait_end)
         rcount = rcount - ldone
         do jproc = 1,ldone
            mproc = donelist(jproc)
            iproc = rlist(mproc)  ! Recv from
            lproc = neighlist(iproc)
            ! unpack_list(iq) is index into extended region
            t(ifull+bnds(lproc)%unpack_list(1:rslen(iproc))) = bnds(lproc)%rbuf(1:rslen(iproc))
         end do
      end do

      ! Clear any remaining messages
      ! MJT notes - We could also call Waitall at the start of the next
      ! call to bounds or boundsuv.  However, this assumes MPI will
      ! process the message in the background which is not always the case.
      ! Instead we call Waitall here to ensure the MPI messages are
      ! progressed.
      sreq = nreq - rreq
      call START_LOG(mpiwait_begin)
      call MPI_Waitall(sreq,ireq(rreq+1:nreq),MPI_STATUSES_IGNORE,ierr)
      call END_LOG(mpiwait_end)

      call END_LOG(bounds_end)

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
      integer :: rcount, myrlen, jproc, mproc
      integer, dimension(neighnum) :: rslen, sslen
      integer(kind=4) :: ierr, itag = 2, llen, sreq, lproc
      integer(kind=4) :: ldone
      integer(kind=4), dimension(neighnum) :: donelist
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
         myrlen = bnds(myid)%rlen2
      else if ( extra ) then
         rslen = bnds(neighlist)%rlenx
         sslen = bnds(neighlist)%slenx
         myrlen = bnds(myid)%rlenx
      else if ( single ) then
         rslen = bnds(neighlist)%rlen
         sslen = bnds(neighlist)%slen
         myrlen = bnds(myid)%rlen
      else
         rslen  = bnds(neighlist)%rlenh
         sslen  = bnds(neighlist)%slenh
         myrlen = bnds(myid)%rlenh
      end if

      call START_LOG(bounds_begin)

!     Set up the buffers to send and recv
      nreq = 0
      do iproc = 1,neighnum
         recv_len = rslen(iproc)
         if ( recv_len > 0 ) then
            lproc = neighlist(iproc)  ! Recv from
            nreq = nreq + 1
            rlist(nreq) = iproc
            llen = recv_len*kx
            call MPI_IRecv( bnds(lproc)%rbuf(1), llen, ltype, lproc, &
                 itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do
      rreq = nreq
      do iproc = neighnum,1,-1
         send_len = sslen(iproc)
         if ( send_len > 0 ) then
            lproc = neighlist(iproc)  ! Send to
            bnds(lproc)%sbuf(1:send_len*kx) = reshape( t(bnds(lproc)%send_list(1:send_len),1:kx), (/ send_len*kx /) )
            nreq = nreq + 1
            llen = send_len*kx
            call MPI_ISend( bnds(lproc)%sbuf(1), llen, ltype, lproc, &
                 itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do

      ! Finally see if there are any points on my own processor that need
      ! to be fixed up. This will only be in the case when nproc < npanels.
      if ( myrlen>0 ) then
         ! request_list is same as send_list in this case
         t(ifull+bnds(myid)%unpack_list(1:myrlen),1:kx) = t(bnds(myid)%request_list(1:myrlen),1:kx)
      end if

      ! Unpack incomming messages
      rcount = rreq
      do while ( rcount > 0 )
         call START_LOG(mpiwait_begin)
         call MPI_Waitsome( rreq, ireq, ldone, donelist, MPI_STATUSES_IGNORE, ierr )
         call END_LOG(mpiwait_end)
         rcount = rcount - ldone
         do jproc = 1,ldone
            mproc = donelist(jproc)
            iproc = rlist(mproc)  ! Recv from
            lproc = neighlist(iproc)
            t(ifull+bnds(lproc)%unpack_list(1:rslen(iproc)),1:kx) &
                = reshape( bnds(lproc)%rbuf(1:rslen(iproc)*kx), (/ rslen(iproc), kx /) )
         end do
      end do

      ! Clear any remaining messages
      sreq = nreq - rreq
      call START_LOG(mpiwait_begin)
      call MPI_Waitall(sreq,ireq(rreq+1:nreq),MPI_STATUSES_IGNORE,ierr)
      call END_LOG(mpiwait_end)

      call END_LOG(bounds_end)

   end subroutine bounds3

   subroutine bounds4(t, nrows, corner, nehalf)
      ! Copy the boundary regions.
      real, dimension(:,:,:), intent(inout) :: t
      integer, intent(in), optional :: nrows
      logical, intent(in), optional :: corner
      logical, intent(in), optional :: nehalf
      logical :: extra, single, double
      integer :: iproc, kx, send_len, recv_len
      integer :: rcount, myrlen, jproc, mproc, ntr
      integer, dimension(neighnum) :: rslen, sslen
      integer(kind=4) :: ierr, itag=3, llen, sreq, lproc
      integer(kind=4) :: ldone
      integer(kind=4), dimension(neighnum) :: donelist
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
         myrlen = bnds(myid)%rlen2
      else if ( extra ) then
         rslen = bnds(neighlist)%rlenx
         sslen = bnds(neighlist)%slenx
         myrlen = bnds(myid)%rlenx
      else if ( single ) then
         rslen = bnds(neighlist)%rlen
         sslen = bnds(neighlist)%slen
         myrlen = bnds(myid)%rlen
      else
         rslen  = bnds(neighlist)%rlenh
         sslen  = bnds(neighlist)%slenh
         myrlen = bnds(myid)%rlenh
      end if

      call START_LOG(bounds_begin)

!     Set up the buffers to send
      nreq = 0
      do iproc = 1,neighnum
         recv_len = rslen(iproc)
         if ( recv_len > 0 ) then
            lproc = neighlist(iproc)  ! Recv from
            nreq = nreq + 1
            rlist(nreq) = iproc
            llen = recv_len*kx*ntr
            call MPI_IRecv( bnds(lproc)%rbuf(1), llen, ltype, lproc, &
                 itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do
      rreq = nreq
      do iproc = neighnum,1,-1
         send_len = sslen(iproc)
         if ( send_len > 0 ) then
            lproc = neighlist(iproc)  ! Send to
            bnds(lproc)%sbuf(1:send_len*kx*ntr) = reshape( t(bnds(lproc)%send_list(1:send_len),1:kx,1:ntr), (/ send_len*kx*ntr /) )
            nreq = nreq + 1
            llen = send_len*kx*ntr
            call MPI_ISend( bnds(lproc)%sbuf(1), llen, ltype, lproc, &
                 itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do

      ! Finally see if there are any points on my own processor that need
      ! to be fixed up. This will only be in the case when nproc < npanels.
      if ( myrlen>0 ) then
         ! request_list is same as send_list in this case
         t(ifull+bnds(myid)%unpack_list(1:myrlen),1:kx,1:ntr) = t(bnds(myid)%request_list(1:myrlen),1:kx,1:ntr)
      end if

      ! Unpack incomming messages
      rcount = rreq
      do while ( rcount > 0 )
         call START_LOG(mpiwait_begin)
         call MPI_Waitsome( rreq, ireq, ldone, donelist, MPI_STATUSES_IGNORE, ierr )
         call END_LOG(mpiwait_end)
         rcount = rcount - ldone
         do jproc = 1,ldone
            mproc = donelist(jproc)
            iproc = rlist(mproc)  ! Recv from
            lproc = neighlist(iproc)
            t(ifull+bnds(lproc)%unpack_list(1:rslen(iproc)),1:kx,1:ntr)                           &
                 = reshape( bnds(lproc)%rbuf(1:rslen(iproc)*kx*ntr), (/ rslen(iproc), kx, ntr /) )
         end do
      end do

      ! Clear any remaining messages
      sreq = nreq - rreq
      call START_LOG(mpiwait_begin)
      call MPI_Waitall(sreq,ireq(rreq+1:nreq),MPI_STATUSES_IGNORE,ierr)
      call END_LOG(mpiwait_end)

      call END_LOG(bounds_end)

   end subroutine bounds4

   subroutine bounds_colour_send(t, lcolour, klim)
      ! Copy the boundary regions. This version allows supports updating
      ! different gridpoint colours
      real, dimension(:,:), intent(in) :: t
      integer, intent(in) :: lcolour
      integer, intent(in), optional :: klim
      integer :: iproc, kx, recv_len, iqq, ibeg, iend
      integer(kind=4) :: ierr, itag=4, llen, lproc
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif 

      call START_LOG(bounds_begin)

      if ( present(klim) ) then
         kx = klim
      else
         kx = size(t,2)
      end if

!     Set up the buffers to send and recv
      nreq = 0
      do iproc = 1,neighnum
         lproc = neighlist(iproc)  ! Recv from
         recv_len = rcolsp(lproc)%ihfn(lcolour)-rcolsp(lproc)%ihbg(lcolour)   &
                   +rcolsp(lproc)%iffn(lcolour)-rcolsp(lproc)%ifbg(lcolour)+2
         if ( recv_len > 0 ) then
            nreq = nreq + 1
            rlist(nreq) = iproc
            llen = recv_len*kx
            call MPI_IRecv( bnds(lproc)%rbuf(1), llen, ltype, lproc, &
                            itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do
      rreq = nreq
      do iproc = neighnum,1,-1
         lproc = neighlist(iproc)  ! Send to
         iqq = 0
         ibeg = scolsp(lproc)%ihbg(lcolour)
         iend = scolsp(lproc)%ihfn(lcolour)
         bnds(lproc)%sbuf(iqq+1:iqq+(iend-ibeg+1)*kx)  &
             = reshape( t(bnds(lproc)%send_list(ibeg:iend),1:kx), (/ (iend-ibeg+1)*kx /) )
         iqq = iqq + (iend-ibeg+1)*kx
         ibeg = scolsp(lproc)%ifbg(lcolour)
         iend = scolsp(lproc)%iffn(lcolour)
         bnds(lproc)%sbuf(iqq+1:iqq+(iend-ibeg+1)*kx)  &
             = reshape( t(bnds(lproc)%send_list(ibeg:iend),1:kx), (/ (iend-ibeg+1)*kx /) )
         iqq = iqq + (iend-ibeg+1)*kx
         if ( iqq > 0 ) then
            nreq = nreq + 1
            llen = iqq
            call MPI_ISend( bnds(lproc)%sbuf(1), llen, ltype, lproc, &
                            itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do

      call END_LOG(bounds_end)

   end subroutine bounds_colour_send
   
   subroutine bounds_colour_recv(t, lcolour, klim)
      ! Copy the boundary regions. This version allows supports updating
      ! different gridpoint colours
      real, dimension(:,:), intent(inout) :: t
      integer, intent(in) :: lcolour
      integer, intent(in), optional :: klim
      integer :: iproc, kx, iqq, ibeg, iend
      integer :: rcount, jproc, myrlen
      integer(kind=4) :: ierr, sreq, lproc, ldone
      integer(kind=4), dimension(neighnum) :: donelist

      call START_LOG(bounds_begin)

      if ( present(klim) ) then
         kx = klim
      else
         kx = size(t,2)
      end if
      myrlen = bnds(myid)%rlen
      
      ! Finally see if there are any points on my own processor that need
      ! to be fixed up. This will only be in the case when nproc < npanels.
      ! request_list is same as send_list in this case
      if ( myrlen>0 ) then
         t(ifull+bnds(myid)%unpack_list(1:myrlen),1:kx) = t(bnds(myid)%request_list(1:myrlen),1:kx)
      end if
      
      ! Unpack incomming messages
      rcount = rreq
      do while ( rcount>0 )
         call START_LOG(mpiwait_begin)
         call MPI_Waitsome( rreq, ireq, ldone, donelist, MPI_STATUSES_IGNORE, ierr )
         call END_LOG(mpiwait_end)
         rcount = rcount - ldone
         do jproc = 1,ldone
            iproc = rlist(donelist(jproc))  ! Recv from
            lproc = neighlist(iproc)
            iqq = 0
            ibeg = rcolsp(lproc)%ihbg(lcolour)
            iend = rcolsp(lproc)%ihfn(lcolour)
            t(ifull+bnds(lproc)%unpack_list(ibeg:iend),1:kx)  &
                = reshape( bnds(lproc)%rbuf(iqq+1:iqq+(iend-ibeg+1)*kx), (/ iend-ibeg+1, kx /) )
            iqq = iqq + (iend-ibeg+1)*kx
            ibeg = rcolsp(lproc)%ifbg(lcolour)
            iend = rcolsp(lproc)%iffn(lcolour)
            t(ifull+bnds(lproc)%unpack_list(ibeg:iend),1:kx)  &
                = reshape( bnds(lproc)%rbuf(iqq+1:iqq+(iend-ibeg+1)*kx), (/ iend-ibeg+1, kx /) )
         end do
      end do

      ! Clear any remaining messages
      sreq = nreq - rreq
      call START_LOG(mpiwait_begin)
      call MPI_Waitall(sreq,ireq(rreq+1:nreq),MPI_STATUSES_IGNORE,ierr)
      call END_LOG(mpiwait_end)

      call END_LOG(bounds_end)

   end subroutine bounds_colour_recv
   
   subroutine boundsuv2(u, v, nrows, stag)
      ! Copy the boundary regions of u and v. This doesn't require the
      ! diagonal points like (0,0), but does have to take care of the
      ! direction changes.
      real, dimension(ifull+iextra), intent(inout) :: u, v
      integer, intent(in), optional :: nrows, stag
      logical :: double
      logical :: fsvwu, fnveu, fssvwwu, fnnveeu
      integer :: iq, iqz, iqt, iproc, iqq, recv_len
      integer :: rcount, myrlen, jproc, mproc, stagmode
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif   
      integer(kind=4) :: ierr, itag=5, llen, sreq, lproc
      integer(kind=4) :: ldone
      integer(kind=4), dimension(neighnum) :: donelist
      real :: tmp

      double = .false.
      stagmode = 0
      if ( present(nrows)) then
         if ( nrows == 2 ) double = .true.
      end if
      if ( present(stag) ) then
         stagmode = stag
      end if

      call START_LOG(boundsuv_begin)
      
      if ( double ) then
         fsvwu = .true.
         fnveu = .true.
         fssvwwu = .true.
         fnnveeu = .true.
         myrlen = bnds(myid)%rlen2_uv
      else if ( stagmode == 1 ) then
         fsvwu = .false.
         fnveu = .true.
         fssvwwu = .false.
         fnnveeu = .true.
         myrlen = bnds(myid)%rlen2_uv
      else if ( stagmode == 2 ) then
         fsvwu = .true.
         fnveu = .true.
         fssvwwu = .false.
         fnnveeu = .true. ! fnnveeu requires fnveu
         myrlen = bnds(myid)%rlen2_uv
      else if ( stagmode == 3 ) then
         fsvwu = .true.
         fnveu = .true.
         fssvwwu = .true. ! fssvwwu requires fsvwu
         fnnveeu = .false.
         myrlen = bnds(myid)%rlen2_uv
      else if ( stagmode == 5 ) then
         fsvwu = .true.
         fnveu = .false.
         fssvwwu = .true.
         fnnveeu = .false.
         myrlen = bnds(myid)%rlen2_uv
      else if ( stagmode == -9 ) then
         fsvwu = .true.
         fnveu = .false.
         fssvwwu = .false.
         fnnveeu = .false.
         myrlen = bnds(myid)%rlen_uv
      else if ( stagmode == -10 ) then
         fsvwu = .false.
         fnveu = .true.
         fssvwwu = .false.
         fnnveeu = .false.
         myrlen = bnds(myid)%rlen_uv
      else
         fsvwu = .true.
         fnveu = .true.
         fssvwwu = .false.
         fnnveeu = .false.
         myrlen = bnds(myid)%rlen_uv
      end if

!     Set up the buffers to send
      nreq = 0
      do iproc = 1,neighnum
         lproc = neighlist(iproc)  ! Recv from
         recv_len = 0
         if ( fsvwu ) then
            recv_len = recv_len+rsplit(lproc)%iwufn-rsplit(lproc)%isvbg+1
         end if
         if ( fnveu ) then
            recv_len = recv_len+rsplit(lproc)%ieufn-rsplit(lproc)%invbg+1
         end if         
         if ( fssvwwu ) then
            recv_len = recv_len+rsplit(lproc)%iwwufn-rsplit(lproc)%issvbg+1
         end if         
         if ( fnnveeu ) then
            recv_len = recv_len+rsplit(lproc)%ieeufn-rsplit(lproc)%innvbg+1
         end if
         if ( recv_len > 0 ) then 
            nreq = nreq + 1
            rlist(nreq) = iproc
            llen = recv_len
            call MPI_IRecv( bnds(lproc)%rbuf(1), llen, ltype, lproc, &
                 itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do
      rreq = nreq
      do iproc = neighnum,1,-1
         lproc = neighlist(iproc)  ! Send to
         ! Build up list of points
         iqq = 0
         if ( fsvwu ) then
!cdir nodep
            do iq=ssplit(lproc)%isvbg,ssplit(lproc)%iwufn
               ! Use abs because sign is used as u/v flag
               iqz = iqq+iq-ssplit(lproc)%isvbg+1
               iqt = bnds(lproc)%send_list_uv(iq)
               if ( (iqt > 0) .neqv. bnds(lproc)%send_swap(iq) ) then
                  bnds(lproc)%sbuf(iqz) = bnds(lproc)%send_neg(iq)*u(abs(iqt))
               else
                  bnds(lproc)%sbuf(iqz) = bnds(lproc)%send_neg(iq)*v(abs(iqt))
               end if 
            end do
            iqq = iqq+ssplit(lproc)%iwufn-ssplit(lproc)%isvbg+1
         end if
         if ( fnveu ) then
!cdir nodep
            do iq=ssplit(lproc)%invbg,ssplit(lproc)%ieufn
               ! Use abs because sign is used as u/v flag
               iqz = iqq+iq-ssplit(lproc)%invbg+1
               iqt = bnds(lproc)%send_list_uv(iq)
               if ( (iqt > 0) .neqv. bnds(lproc)%send_swap(iq) ) then
                  bnds(lproc)%sbuf(iqz) = bnds(lproc)%send_neg(iq)*u(abs(iqt))
               else
                  bnds(lproc)%sbuf(iqz) = bnds(lproc)%send_neg(iq)*v(abs(iqt))
               end if 
            end do
            iqq = iqq+ssplit(lproc)%ieufn-ssplit(lproc)%invbg+1
         end if
         if ( fssvwwu ) then
!cdir nodep
            do iq=ssplit(lproc)%issvbg,ssplit(lproc)%iwwufn
               ! Use abs because sign is used as u/v flag
               iqz = iqq+iq-ssplit(lproc)%issvbg+1
               iqt = bnds(lproc)%send_list_uv(iq)
               if ( (iqt > 0) .neqv. bnds(lproc)%send_swap(iq) ) then
                  bnds(lproc)%sbuf(iqz) = bnds(lproc)%send_neg(iq)*u(abs(iqt))
               else
                  bnds(lproc)%sbuf(iqz) = bnds(lproc)%send_neg(iq)*v(abs(iqt))
               end if 
            end do
            iqq = iqq+ssplit(lproc)%iwwufn-ssplit(lproc)%issvbg+1
         end if
         if ( fnnveeu ) then
!cdir nodep
            do iq=ssplit(lproc)%innvbg,ssplit(lproc)%ieeufn
               ! Use abs because sign is used as u/v flag
               iqz = iqq+iq-ssplit(lproc)%innvbg+1
               iqt = bnds(lproc)%send_list_uv(iq)
               if ( iqt > 0 .neqv. bnds(lproc)%send_swap(iq) ) then
                  bnds(lproc)%sbuf(iqz) = bnds(lproc)%send_neg(iq)*u(abs(iqt))
               else
                  bnds(lproc)%sbuf(iqz) = bnds(lproc)%send_neg(iq)*v(abs(iqt))
               end if 
            end do
            iqq = iqq+ssplit(lproc)%ieeufn-ssplit(lproc)%innvbg+1
         end if
         if ( iqq > 0 ) then
            nreq = nreq + 1
            llen = iqq
            call MPI_ISend( bnds(lproc)%sbuf(1), llen, ltype, lproc, &
                 itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do

      ! Finally see if there are any points on my own processor that need
      ! to be fixed up. This will only be in the case when nproc < npanels.
!cdir nodep
      do iq = 1,myrlen
         ! request_list is same as send_list in this case
         if ( ( bnds(myid)%request_list_uv(iq) > 0) .neqv. &
                   bnds(myid)%uv_swap(iq) ) then  ! haven't copied to send_swap yet
            tmp = u(abs(bnds(myid)%request_list_uv(iq)))
         else
            tmp = v(abs(bnds(myid)%request_list_uv(iq)))
         end if
         if ( bnds(myid)%uv_neg(iq) ) tmp = -tmp
         ! unpack_list(iq) is index into extended region
         if ( bnds(myid)%unpack_list_uv(iq) > 0 ) then
            u(ifull+bnds(myid)%unpack_list_uv(iq)) = tmp
         else
            v(ifull-bnds(myid)%unpack_list_uv(iq)) = tmp
         end if
      end do

      ! Unpack incomming messages
      rcount = rreq
      do while ( rcount > 0 )
         call START_LOG(mpiwaituv_begin)
         call MPI_Waitsome( rreq, ireq, ldone, donelist, MPI_STATUSES_IGNORE, ierr )
         call END_LOG(mpiwaituv_end)
         rcount = rcount - ldone
         do jproc = 1,ldone
            mproc = donelist(jproc)
            iproc = rlist(mproc)  ! Recv from
            lproc = neighlist(iproc)
            iqq = 0
            if ( fsvwu ) then
!cdir nodep
               do iq=rsplit(lproc)%isvbg,rsplit(lproc)%iwufn
                  ! unpack_list(iq) is index into extended region
                  iqz = iqq+iq-rsplit(lproc)%isvbg+1
                  iqt = bnds(lproc)%unpack_list_uv(iq) 
                  if ( iqt > 0 ) then
                     u(ifull+iqt) = bnds(lproc)%rbuf(iqz)
                  else
                     v(ifull-iqt) = bnds(lproc)%rbuf(iqz)
                  end if
               end do
               iqq = iqq+rsplit(lproc)%iwufn-rsplit(lproc)%isvbg+1
            end if
            if ( fnveu ) then
!cdir nodep
               do iq=rsplit(lproc)%invbg,rsplit(lproc)%ieufn
                  ! unpack_list(iq) is index into extended region
                  iqz = iqq+iq-rsplit(lproc)%invbg+1
                  iqt = bnds(lproc)%unpack_list_uv(iq) 
                  if ( iqt > 0 ) then
                     u(ifull+iqt) = bnds(lproc)%rbuf(iqz)
                  else
                     v(ifull-iqt) = bnds(lproc)%rbuf(iqz)
                  end if
               end do
               iqq = iqq+rsplit(lproc)%ieufn-rsplit(lproc)%invbg+1
            end if
            if ( fssvwwu ) then
!cdir nodep
               do iq=rsplit(lproc)%issvbg,rsplit(lproc)%iwwufn
                  ! unpack_list(iq) is index into extended region
                  iqz = iqq+iq-rsplit(lproc)%issvbg+1
                  iqt = bnds(lproc)%unpack_list_uv(iq)
                  if ( iqt > 0 ) then
                     u(ifull+iqt) = bnds(lproc)%rbuf(iqz)
                  else
                     v(ifull-iqt) = bnds(lproc)%rbuf(iqz)
                  end if
               end do
               iqq = iqq+rsplit(lproc)%iwwufn-rsplit(lproc)%issvbg+1
            end if
            if ( fnnveeu ) then
!cdir nodep
               do iq=rsplit(lproc)%innvbg,rsplit(lproc)%ieeufn
                  ! unpack_list(iq) is index into extended region
                  iqz = iqq+iq-rsplit(lproc)%innvbg+1
                  iqt = bnds(lproc)%unpack_list_uv(iq) 
                  if ( iqt > 0 ) then
                     u(ifull+iqt) = bnds(lproc)%rbuf(iqz)
                  else
                     v(ifull-iqt) = bnds(lproc)%rbuf(iqz)
                  end if
               end do
               iqq = iqq+rsplit(lproc)%ieeufn-rsplit(lproc)%innvbg+1
            end if
         end do
      end do

      ! Clear any remaining messages
      sreq = nreq - rreq
      call START_LOG(mpiwaituv_begin)
      call MPI_Waitall(sreq,ireq(rreq+1:nreq),MPI_STATUSES_IGNORE,ierr)
      call END_LOG(mpiwaituv_end)

      call END_LOG(boundsuv_end)

   end subroutine boundsuv2

   subroutine boundsuv3(u, v, nrows, stag)
      ! Copy the boundary regions of u and v. This doesn't require the
      ! diagonal points like (0,0), but does have to take care of the
      ! direction changes.
      real, dimension(:,:), intent(inout) :: u, v
      integer, intent(in), optional :: nrows
      integer, intent(in), optional :: stag
      logical :: double
      logical :: fsvwu, fnveu, fssvwwu, fnnveeu
      integer :: iq, iqz, iqt, iproc, kx, rproc, sproc, iqq, recv_len
      integer :: rcount, myrlen, jproc, mproc, stagmode
      integer(kind=4) :: ierr, itag=6, llen, sreq, lproc
      integer(kind=4) :: ldone
      integer(kind=4), dimension(neighnum) :: donelist
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif   
      real, dimension(size(u,2)) :: tmp
      
      kx = size(u,2)
      double = .false.
      stagmode = 0
      if ( present(nrows) ) then
         if ( nrows == 2 ) then
            double = .true.
         end if
      end if
      if ( present(stag) ) then
         stagmode = stag
      end if

      call START_LOG(boundsuv_begin)
      
      if ( double ) then
         fsvwu = .true.
         fnveu = .true.
         fssvwwu = .true.
         fnnveeu = .true.
         myrlen = bnds(myid)%rlen2_uv
      else if ( stagmode == 1 ) then
         fsvwu = .false.
         fnveu = .true.
         fssvwwu = .false.
         fnnveeu = .true.
         myrlen = bnds(myid)%rlen2_uv
      else if ( stagmode == 2 ) then
         fsvwu = .true.
         fnveu = .true.
         fssvwwu = .false.
         fnnveeu = .true. ! fnnveeu requires fnveu
         myrlen = bnds(myid)%rlen2_uv
      else if ( stagmode == 3 ) then
         fsvwu = .true.
         fnveu = .true.
         fssvwwu = .true. ! fssvwwu requires fsvwu
         fnnveeu = .false.
         myrlen = bnds(myid)%rlen2_uv
      else if ( stagmode == 5 ) then
         fsvwu = .true.
         fnveu = .false.
         fssvwwu = .true.
         fnnveeu = .false.
         myrlen = bnds(myid)%rlen2_uv
      else if ( stagmode == -9 ) then
         fsvwu = .true.
         fnveu = .false.
         fssvwwu = .false.
         fnnveeu = .false.
         myrlen = bnds(myid)%rlen_uv
      else if ( stagmode == -10 ) then
         fsvwu = .false.
         fnveu = .true.
         fssvwwu = .false.
         fnnveeu = .false.
         myrlen = bnds(myid)%rlen_uv
      else
         fsvwu = .true.
         fnveu = .true.
         fssvwwu = .false.
         fnnveeu = .false.
         myrlen = bnds(myid)%rlen_uv
      end if

!     Set up the buffers to send and recv
      nreq = 0
      do iproc = 1,neighnum
         rproc = neighlist(iproc)  ! Recv from
         recv_len = 0
         if ( fsvwu ) then
            recv_len = recv_len+rsplit(rproc)%iwufn-rsplit(rproc)%isvbg+1
         end if
         if ( fnveu ) then
            recv_len = recv_len+rsplit(rproc)%ieufn-rsplit(rproc)%invbg+1
         end if         
         if ( fssvwwu ) then
            recv_len = recv_len+rsplit(rproc)%iwwufn-rsplit(rproc)%issvbg+1
         end if         
         if ( fnnveeu ) then
            recv_len = recv_len+rsplit(rproc)%ieeufn-rsplit(rproc)%innvbg+1
         end if
         if ( recv_len > 0 ) then 
            nreq = nreq + 1
            rlist(nreq) = iproc
            llen = recv_len*kx
            lproc = rproc
            call MPI_IRecv( bnds(rproc)%rbuf(1), llen, ltype, lproc, &
                 itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do
      rreq = nreq
      do iproc = neighnum,1,-1
         sproc = neighlist(iproc)  ! Send to
         ! Build up list of points
         iqq = 0
         if ( fsvwu ) then
!cdir nodep
            do iq = ssplit(sproc)%isvbg,ssplit(sproc)%iwufn
               ! send_list_uv(iq) is point index.
               ! Use abs because sign is used as u/v flag
               iqz = iqq + iq - ssplit(sproc)%isvbg + 1
               iqt = bnds(sproc)%send_list_uv(iq)
               if ( (iqt > 0) .neqv. bnds(sproc)%send_swap(iq) ) then
                  bnds(sproc)%sbuf(1+(iqz-1)*kx:iqz*kx) = bnds(sproc)%send_neg(iq)*u(abs(iqt),1:kx)
               else
                  bnds(sproc)%sbuf(1+(iqz-1)*kx:iqz*kx) = bnds(sproc)%send_neg(iq)*v(abs(iqt),1:kx)
               end if
            end do
            iqq = iqq+ssplit(sproc)%iwufn-ssplit(sproc)%isvbg+1
         end if
         if ( fnveu ) then
!cdir nodep
            do iq = ssplit(sproc)%invbg,ssplit(sproc)%ieufn
               ! send_list_uv(iq) is point index.
               ! Use abs because sign is used as u/v flag
               iqz = iqq + iq - ssplit(sproc)%invbg + 1
               iqt = bnds(sproc)%send_list_uv(iq)
               if ( (iqt > 0) .neqv. bnds(sproc)%send_swap(iq) ) then
                  bnds(sproc)%sbuf(1+(iqz-1)*kx:iqz*kx) = bnds(sproc)%send_neg(iq)*u(abs(iqt),1:kx)
               else
                  bnds(sproc)%sbuf(1+(iqz-1)*kx:iqz*kx) = bnds(sproc)%send_neg(iq)*v(abs(iqt),1:kx)
               end if 
            end do
            iqq = iqq+ssplit(sproc)%ieufn-ssplit(sproc)%invbg+1
         end if
         if ( fssvwwu ) then
!cdir nodep
            do iq = ssplit(sproc)%issvbg,ssplit(sproc)%iwwufn
               ! send_list_uv(iq) is point index.
               ! Use abs because sign is used as u/v flag
               iqz = iqq + iq - ssplit(sproc)%issvbg + 1
               iqt = bnds(sproc)%send_list_uv(iq)
               if ( (iqt > 0) .neqv. bnds(sproc)%send_swap(iq) ) then
                  bnds(sproc)%sbuf(1+(iqz-1)*kx:iqz*kx) = bnds(sproc)%send_neg(iq)*u(abs(iqt),1:kx)
               else
                  bnds(sproc)%sbuf(1+(iqz-1)*kx:iqz*kx) = bnds(sproc)%send_neg(iq)*v(abs(iqt),1:kx)
               end if 
            end do
            iqq = iqq+ssplit(sproc)%iwwufn-ssplit(sproc)%issvbg+1
         end if
         if ( fnnveeu ) then
!cdir nodep
            do iq = ssplit(sproc)%innvbg,ssplit(sproc)%ieeufn
               ! send_list_uv(iq) is point index.
               ! Use abs because sign is used as u/v flag
               iqz = iqq + iq - ssplit(sproc)%innvbg + 1
               iqt = bnds(sproc)%send_list_uv(iq)
               if ( (iqt > 0) .neqv. bnds(sproc)%send_swap(iq) ) then
                  bnds(sproc)%sbuf(1+(iqz-1)*kx:iqz*kx) = bnds(sproc)%send_neg(iq)*u(abs(iqt),1:kx)
               else
                  bnds(sproc)%sbuf(1+(iqz-1)*kx:iqz*kx) = bnds(sproc)%send_neg(iq)*v(abs(iqt),1:kx)
               end if
            end do
            iqq = iqq + ssplit(sproc)%ieeufn - ssplit(sproc)%innvbg + 1
         end if
         if ( iqq>0 ) then
            nreq = nreq + 1
            llen = iqq*kx
            lproc = sproc
            call MPI_ISend( bnds(sproc)%sbuf(1), llen, ltype, lproc, &
                 itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do

      ! Finally see if there are any points on my own processor that need
      ! to be fixed up. This will only be in the case when nproc < npanels.
!cdir nodep
      do iq = 1,myrlen
         ! request_list is same as send_list in this case
         if ( (bnds(myid)%request_list_uv(iq) > 0) .neqv. &
                  bnds(myid)%uv_swap(iq) ) then  ! haven't copied to send_swap yet
            tmp(1:kx) = u(abs(bnds(myid)%request_list_uv(iq)),1:kx)
         else
            tmp(1:kx) = v(abs(bnds(myid)%request_list_uv(iq)),1:kx)
         end if
         if ( bnds(myid)%uv_neg(iq) ) tmp(1:kx) = -tmp(1:kx)

         ! unpack_list(iq) is index into extended region
         if ( bnds(myid)%unpack_list_uv(iq) > 0 ) then
            u(ifull+bnds(myid)%unpack_list_uv(iq),1:kx) = tmp(1:kx)
         else
            v(ifull-bnds(myid)%unpack_list_uv(iq),1:kx) = tmp(1:kx)
         end if
      end do

      ! Unpack incomming messages
      rcount = rreq
      do while ( rcount > 0 )
         call START_LOG(mpiwaituv_begin)
         call MPI_Waitsome( rreq, ireq, ldone, donelist, MPI_STATUSES_IGNORE, ierr )
         call END_LOG(mpiwaituv_end)
         rcount = rcount - ldone
         do jproc = 1,ldone
            mproc = donelist(jproc)
            iproc = rlist(mproc)  ! Recv from
            rproc = neighlist(iproc)
            iqq = 0
            if ( fsvwu ) then
!cdir nodep
               do iq = rsplit(rproc)%isvbg,rsplit(rproc)%iwufn
                  ! unpack_list(iq) is index into extended region
                  iqz = iqq + iq - rsplit(rproc)%isvbg + 1
                  iqt = bnds(rproc)%unpack_list_uv(iq) 
                  if ( iqt > 0 ) then
                     u(ifull+iqt,1:kx) = bnds(rproc)%rbuf(1+(iqz-1)*kx:iqz*kx)
                  else
                     v(ifull-iqt,1:kx) = bnds(rproc)%rbuf(1+(iqz-1)*kx:iqz*kx)
                  end if
               end do
               iqq = iqq+rsplit(rproc)%iwufn-rsplit(rproc)%isvbg+1
            end if
            if ( fnveu ) then
!cdir nodep
               do iq = rsplit(rproc)%invbg,rsplit(rproc)%ieufn
                  ! unpack_list(iq) is index into extended region
                  iqz = iqq + iq - rsplit(rproc)%invbg + 1
                  iqt = bnds(rproc)%unpack_list_uv(iq)
                  if ( iqt > 0 ) then
                     u(ifull+iqt,1:kx) = bnds(rproc)%rbuf(1+(iqz-1)*kx:iqz*kx)
                  else
                     v(ifull-iqt,1:kx) = bnds(rproc)%rbuf(1+(iqz-1)*kx:iqz*kx)
                  end if
               end do
               iqq = iqq+rsplit(rproc)%ieufn-rsplit(rproc)%invbg+1
            end if         
            if ( fssvwwu ) then
!cdir nodep
               do iq = rsplit(rproc)%issvbg,rsplit(rproc)%iwwufn
                  ! unpack_list(iq) is index into extended region
                  iqz = iqq + iq - rsplit(rproc)%issvbg + 1
                  iqt = bnds(rproc)%unpack_list_uv(iq)
                  if ( iqt > 0 ) then
                     u(ifull+iqt,1:kx) = bnds(rproc)%rbuf(1+(iqz-1)*kx:iqz*kx)
                  else
                     v(ifull-iqt,1:kx) = bnds(rproc)%rbuf(1+(iqz-1)*kx:iqz*kx)
                  end if
               end do
               iqq = iqq+rsplit(rproc)%iwwufn-rsplit(rproc)%issvbg+1
            end if         
            if ( fnnveeu ) then
!cdir nodep
               do iq = rsplit(rproc)%innvbg,rsplit(rproc)%ieeufn
                  ! unpack_list(iq) is index into extended region
                  iqz = iqq + iq - rsplit(rproc)%innvbg + 1
                  iqt = bnds(rproc)%unpack_list_uv(iq)
                  if ( iqt > 0 ) then
                     u(ifull+iqt,1:kx) = bnds(rproc)%rbuf(1+(iqz-1)*kx:iqz*kx)
                  else
                     v(ifull-iqt,1:kx) = bnds(rproc)%rbuf(1+(iqz-1)*kx:iqz*kx)
                  end if
               end do
               iqq = iqq+rsplit(rproc)%ieeufn-rsplit(rproc)%innvbg+1
            end if         
         end do
      end do

      ! Clear any remaining messages
      sreq = nreq - rreq
      call START_LOG(mpiwaituv_begin)
      call MPI_Waitall(sreq,ireq(rreq+1:nreq),MPI_STATUSES_IGNORE,ierr)
      call END_LOG(mpiwaituv_end)

      call END_LOG(boundsuv_end)

   end subroutine boundsuv3
   
   subroutine boundsuv_allvec(u, v)
      ! Copy the boundary regions of u and v. This doesn't require the
      ! diagonal points like (0,0), but does have to take care of the
      ! direction changes.
      real, dimension(:,:), intent(inout) :: u, v
      integer :: iq, iqz, iqt, iproc, kx, rproc, sproc, iqq, recv_len
      integer :: rcount, myrlen, jproc, mproc
      integer(kind=4) :: ierr, itag=7, llen, sreq, lproc
      integer(kind=4) :: ldone
      integer(kind=4), dimension(neighnum) :: donelist
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif   
      real, dimension(1:size(u,2)) :: tmp
      
      kx = size(u,2)
      myrlen = bnds(myid)%rlenx_uv

      call START_LOG(boundsuv_begin)

!     Set up the buffers to recv and send
      nreq = 0
      do iproc = 1,neighnum
         rproc = neighlist(iproc)  ! Recv from
         recv_len = rsplit(rproc)%ieufn - rsplit(rproc)%isvbg + 1
         recv_len = recv_len + rsplit(rproc)%ievfn - rsplit(rproc)%isubg + 1
         if ( recv_len > 0 ) then 
            nreq = nreq + 1
            rlist(nreq) = iproc
            llen = recv_len*kx
            lproc = rproc
            call MPI_IRecv( bnds(rproc)%rbuf(1), llen, ltype, lproc, &
                 itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do
      rreq = nreq
      do iproc = neighnum,1,-1
         sproc = neighlist(iproc)  ! Send to
         ! Build up list of points
         do iq = ssplit(sproc)%isvbg,ssplit(sproc)%ieufn
            ! send_list_uv(iq) is point index.
            ! Use abs because sign is used as u/v flag
            iqz = iq - ssplit(sproc)%isvbg+1
            iqt = bnds(sproc)%send_list_uv(iq)
            if ( (iqt > 0) .neqv. bnds(sproc)%send_swap(iq) ) then
               bnds(sproc)%sbuf(1+(iqz-1)*kx:iqz*kx) = bnds(sproc)%send_neg(iq)*u(abs(iqt),1:kx)
            else
               bnds(sproc)%sbuf(1+(iqz-1)*kx:iqz*kx) = bnds(sproc)%send_neg(iq)*v(abs(iqt),1:kx)
            end if
         end do
         iqq = ssplit(sproc)%ieufn-ssplit(sproc)%isvbg+1
         do iq = ssplit(sproc)%isubg,ssplit(sproc)%ievfn
            ! send_list_uv(iq) is point index.
            ! Use abs because sign is used as u/v flag
            iqz = iqq + iq - ssplit(sproc)%isubg + 1
            iqt = bnds(sproc)%send_list_uv(iq)
            if ( (iqt > 0) .neqv. bnds(sproc)%send_swap(iq) ) then
               bnds(sproc)%sbuf(1+(iqz-1)*kx:iqz*kx) = bnds(sproc)%send_neg(iq)*u(abs(iqt),1:kx)
            else
               bnds(sproc)%sbuf(1+(iqz-1)*kx:iqz*kx) = bnds(sproc)%send_neg(iq)*v(abs(iqt),1:kx)
            end if 
         end do
         iqq = iqq + ssplit(sproc)%ievfn - ssplit(sproc)%isubg + 1
         if ( iqq>0 ) then
            nreq = nreq + 1
            llen = iqq*kx
            lproc = sproc
            call MPI_ISend( bnds(sproc)%sbuf(1), llen, ltype, lproc, &
                 itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do

      ! Finally see if there are any points on my own processor that need
      ! to be fixed up. This will only be in the case when nproc < npanels.
      do iq = 1,myrlen
         ! request_list is same as send_list in this case
         if ( (bnds(myid)%request_list_uv(iq) > 0) .neqv. bnds(myid)%uv_swap(iq) ) then  ! haven't copied to send_swap yet
            tmp(1:kx) = u(abs(bnds(myid)%request_list_uv(iq)),1:kx)
         else
            tmp(1:kx) = v(abs(bnds(myid)%request_list_uv(iq)),1:kx)
         end if
         if ( bnds(myid)%uv_neg(iq) ) tmp(1:kx) = -tmp(1:kx)

         ! unpack_list(iq) is index into extended region
         if ( bnds(myid)%unpack_list_uv(iq) > 0 ) then
            u(ifull+bnds(myid)%unpack_list_uv(iq),1:kx) = tmp(1:kx)
         else
            v(ifull-bnds(myid)%unpack_list_uv(iq),1:kx) = tmp(1:kx)
         end if
      end do

      ! Unpack incomming messages
      rcount = rreq
      do while ( rcount > 0 )
         call START_LOG(mpiwaituv_begin)
         call MPI_Waitsome( rreq, ireq, ldone, donelist, MPI_STATUSES_IGNORE, ierr )
         call END_LOG(mpiwaituv_end)
         rcount = rcount - ldone
         do jproc = 1,ldone
            mproc = donelist(jproc)
            iproc = rlist(mproc)  ! Recv from
            rproc = neighlist(iproc)
            do iq = rsplit(rproc)%isvbg,rsplit(rproc)%ieufn
               ! unpack_list(iq) is index into extended region
               iqz = iq - rsplit(rproc)%isvbg + 1
               iqt = bnds(rproc)%unpack_list_uv(iq) 
               if ( iqt > 0 ) then
                  u(ifull+iqt,1:kx) = bnds(rproc)%rbuf(1+(iqz-1)*kx:iqz*kx)
               else
                  v(ifull-iqt,1:kx) = bnds(rproc)%rbuf(1+(iqz-1)*kx:iqz*kx)
               end if
            end do
            iqq = rsplit(rproc)%ieufn - rsplit(rproc)%isvbg + 1
            do iq = rsplit(rproc)%isubg,rsplit(rproc)%ievfn
               ! unpack_list(iq) is index into extended region
               iqz = iqq + iq - rsplit(rproc)%isubg + 1
               iqt = bnds(rproc)%unpack_list_uv(iq)
               if ( iqt > 0 ) then
                  u(ifull+iqt,1:kx) = bnds(rproc)%rbuf(1+(iqz-1)*kx:iqz*kx)
               else
                  v(ifull-iqt,1:kx) = bnds(rproc)%rbuf(1+(iqz-1)*kx:iqz*kx)
               end if
            end do
            iqq = iqq + rsplit(rproc)%ievfn - rsplit(rproc)%isubg + 1
         end do
      end do

      ! Clear any remaining messages
      sreq = nreq - rreq
      call START_LOG(mpiwaituv_begin)
      call MPI_Waitall(sreq,ireq(rreq+1:nreq),MPI_STATUSES_IGNORE,ierr)
      call END_LOG(mpiwaituv_end)

      call END_LOG(boundsuv_end)

   end subroutine boundsuv_allvec

   subroutine bounds3r8(t, nrows, klim, corner, nehalf)
      ! Copy the boundary regions. Only this routine requires the extra klim
      ! argument (for helmsol).
      real(kind=8), dimension(:,:), intent(inout) :: t
      integer, intent(in), optional :: nrows, klim
      logical, intent(in), optional :: corner
      logical, intent(in), optional :: nehalf
      logical :: extra, single, double
      integer :: iproc, kx, send_len, recv_len
      integer :: rcount, myrlen, jproc, mproc
      integer, dimension(neighnum) :: rslen, sslen
      integer(kind=4) :: ierr, itag = 2, llen, sreq, lproc
      integer(kind=4) :: ldone
      integer(kind=4), dimension(neighnum) :: donelist
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
         myrlen = bnds(myid)%rlen2
      else if ( extra ) then
         rslen = bnds(neighlist)%rlenx
         sslen = bnds(neighlist)%slenx
         myrlen = bnds(myid)%rlenx
      else if ( single ) then
         rslen = bnds(neighlist)%rlen
         sslen = bnds(neighlist)%slen
         myrlen = bnds(myid)%rlen
      else
         rslen  = bnds(neighlist)%rlenh
         sslen  = bnds(neighlist)%slenh
         myrlen = bnds(myid)%rlenh
      end if

      call START_LOG(bounds_begin)

!     Set up the buffers to send and recv
      nreq = 0
      do iproc = 1,neighnum
         recv_len = rslen(iproc)
         if ( recv_len > 0 ) then
            lproc = neighlist(iproc)  ! Recv from
            nreq = nreq + 1
            rlist(nreq) = iproc
            llen = recv_len*kx
            call MPI_IRecv( bnds(lproc)%r8buf(1), llen, ltype, lproc, &
                 itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do
      rreq = nreq
      do iproc = neighnum,1,-1
         send_len = sslen(iproc)
         if ( send_len > 0 ) then
            lproc = neighlist(iproc)  ! Send to
            bnds(lproc)%s8buf(1:send_len*kx) = reshape( t(bnds(lproc)%send_list(1:send_len),1:kx), (/ send_len*kx /) )
            nreq = nreq + 1
            llen = send_len*kx
            call MPI_ISend( bnds(lproc)%s8buf(1), llen, ltype, lproc, &
                 itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do

      ! Finally see if there are any points on my own processor that need
      ! to be fixed up. This will only be in the case when nproc < npanels.
      if ( myrlen>0 ) then
         ! request_list is same as send_list in this case
         t(ifull+bnds(myid)%unpack_list(1:myrlen),1:kx) = t(bnds(myid)%request_list(1:myrlen),1:kx)
      end if

      ! Unpack incomming messages
      rcount = rreq
      do while ( rcount>0 )
         call START_LOG(mpiwait_begin)
         call MPI_Waitsome( rreq, ireq, ldone, donelist, MPI_STATUSES_IGNORE, ierr )
         call END_LOG(mpiwait_end)
         rcount = rcount - ldone
         do jproc = 1,ldone
            mproc = donelist(jproc)
            iproc = rlist(mproc)  ! Recv from
            lproc = neighlist(iproc)
            t(ifull+bnds(lproc)%unpack_list(1:rslen(iproc)),1:kx) &
                = reshape( bnds(lproc)%r8buf(1:rslen(iproc)*kx), (/ rslen(iproc), kx /) )
         end do
      end do

      ! Clear any remaining messages
      sreq = nreq - rreq
      call START_LOG(mpiwait_begin)
      call MPI_Waitall(sreq,ireq(rreq+1:nreq),MPI_STATUSES_IGNORE,ierr)
      call END_LOG(mpiwait_end)

      call END_LOG(bounds_end)

   end subroutine bounds3r8
   
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
      real :: xf, yf
      integer :: iproc, jproc, dproc
      integer :: ip, jp, xn, kx
      integer :: iq, k, idel, jdel, nf, gf
      integer :: rcount
      integer(kind=4) :: itag=99, ierr, llen, ncount, sreq, lproc
      integer(kind=4) :: ldone
      integer(kind=4), dimension(MPI_STATUS_SIZE,neighnum) :: status
      integer(kind=4), dimension(neighnum) :: donelist
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif   

      ! This does nothing in the one processor case
      if ( neighnum < 1 ) return

      call START_LOG(deptsync_begin)      

      kx = size(nface,2)
      dslen(:) = 0
      drlen(:) = 0
      dproc = 0
      
!     In this case the length of each buffer is unknown and will not
!     be symmetric between processors. Therefore need to get the length
!     from the message status
      do iproc = 1,neighnum
         lproc = neighlist(iproc)  ! Recv from
         ! Use the maximum size in the recv call.
         llen = 4*bnds(lproc)%len/nagg
         call MPI_IRecv( dpoints(iproc)%a, llen, ltype, lproc, &
                      itag, MPI_COMM_WORLD, ireq(iproc), ierr )
      end do
      nreq = neighnum
      
      ! Calculate request list
      do k = 1,kx
         do iq = 1,ifull
            gf = nface(iq,k)
            xf = xg(iq,k)
            yf = yg(iq,k)
            nf = gf + noff ! Make this a local index
            idel = int(xf) - ioff
            jdel = int(yf) - joff
            if ( idel<0 .or. idel>ipan .or. jdel<0 .or. jdel>jpan .or. nf<1 .or. nf>npan ) then
               ! If point is on a different processor, add to a list 
               ip = min( il_g, max( 1, nint(xf) ) )
               jp = min( il_g, max( 1, nint(yf) ) )
               iproc = fproc(ip,jp,gf) ! processor that owns global grid point
               dproc = neighmap(iproc) ! returns 0 if not in neighlist
               ! Add this point to the list of requests I need to send to iproc
               dslen(dproc) = dslen(dproc) + 1
               ! Limit request index to valid range to avoid seg fault
               xn = max( min( dslen(dproc), bnds(iproc)%len/nagg ), 1 )
               ! Since nface is a small integer it can be exactly represented by a
               ! real. It is simpler to send like this than use a proper structure.
               dbuf(dproc)%a(:,xn) = (/ real(gf), xf, yf, real(k) /)
               dindex(dproc)%a(:,xn) = (/ iq, k /)
            end if
         end do
      end do
 
      ! Check for errors
      if ( dslen(0) > 0 ) then
         write(6,*) "myid,dslen(0) ",myid,dslen(0)
         gf = nint(dbuf(dproc)%a(1,1))
         ip = min( il_g, max( 1, nint(dbuf(dproc)%a(2,1)) ) )
         jp = min( il_g, max( 1, nint(dbuf(dproc)%a(3,1)) ) )
         iproc = fproc(ip,jp,gf)
         write(6,*) "Example error iq,k,iproc ",dindex(0)%a(:,1),iproc
         write(6,*) "dbuf ", dbuf(0)%a(:,1)
         write(6,*) "neighlist ",neighlist
         call checksize( dslen(0), 0, "Deptsync" )
      end if
      do dproc = 1,neighnum
         iproc = neighlist(dproc)
         if ( dslen(dproc) > bnds(iproc)%len/nagg ) then
            write(6,*) "myid,iproc,dslen,len ",myid,iproc,dslen(dproc),bnds(iproc)%len/nagg
            write(6,*) "Example error iq,k, ",dindex(dproc)%a(:,1)
            write(6,*) "dbuf ",dbuf(dproc)%a(:,1)
            write(6,*) "neighlist ",neighlist
            call checksize( dslen(dproc), bnds(iproc)%len/nagg, "Deptsync" )
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
                 itag, MPI_COMM_WORLD, ireq(nreq), ierr )
      end do
      
      ! Unpack incomming messages
      rcount = rreq
      do while ( rcount > 0 )
         call START_LOG(mpiwaitdep_begin)
         call MPI_Waitsome( rreq, ireq, ldone, donelist, status, ierr )
         call END_LOG(mpiwaitdep_end)
         rcount = rcount - ldone
         do jproc = 1,ldone
!           Now get the actual sizes from the status
            call MPI_Get_count( status(:,jproc), ltype, ncount, ierr )
            drlen(donelist(jproc)) = ncount/4
         end do
      end do

      ! Clear any remaining message requests
      sreq = nreq - rreq
      call START_LOG(mpiwaitdep_begin)
      call MPI_Waitall( sreq, ireq(rreq+1:nreq), MPI_STATUSES_IGNORE, ierr )
      call END_LOG(mpiwaitdep_end)
      
      call END_LOG(deptsync_end)

   end subroutine deptsync

   subroutine intssync_send(ntr)
      integer, intent(in) :: ntr
      integer :: iproc
      integer(kind=4) :: itag=98, ierr, llen, lproc
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif

      call START_LOG(intssync_begin)
      
      ! When sending the results, roles of dslen and drlen are reversed
      nreq = 0
      do iproc = 1,neighnum
         if ( dslen(iproc) > 0 ) then
            nreq = nreq + 1
            rlist(nreq) = iproc
            llen = dslen(iproc)*ntr
            lproc = neighlist(iproc)  ! Recv from
            call MPI_IRecv( dbuf(iproc)%b, llen, ltype, lproc, itag, &
                            MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do
      ! MJT notes - we could further split this subroutine into two subroutines,
      ! so that work can be done between the IRecv and ISend.  However, it is
      ! not clear if there is enough work to ensure IRecv are posted early given
      ! the load imbalance arising from the physics.
      rreq = nreq
      do iproc = neighnum,1,-1
         if ( drlen(iproc) > 0 ) then
            nreq = nreq + 1
            llen = drlen(iproc)*ntr
            lproc = neighlist(iproc)  ! Send to
            call MPI_ISend( sextra(iproc)%a, llen, ltype, lproc, itag, & 
                            MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do
      
      call END_LOG(intssync_end)

   end subroutine intssync_send

   subroutine intssync_recv(s)
      real, dimension(:,:,:), intent(inout) :: s
      integer :: iproc, iq, jproc
      integer :: rcount, ntr
      integer(kind=4) :: ierr, ldone, sreq
      integer(kind=4), dimension(neighnum) :: donelist

      call START_LOG(intssync_begin)
      
      ntr = size( s, 3)
      
      ! Unpack incomming messages
      rcount = rreq
      do while ( rcount > 0 )
         call START_LOG(mpiwaitdep_begin)
         call MPI_Waitsome( rreq, ireq, ldone, donelist, MPI_STATUSES_IGNORE, ierr )
         call END_LOG(mpiwaitdep_end)
         rcount = rcount - ldone
         do jproc = 1,ldone
            iproc = rlist(donelist(jproc))
            do iq = 1,dslen(iproc)
               s(dindex(iproc)%a(1,iq),dindex(iproc)%a(2,iq),1:ntr) = dbuf(iproc)%b(1+(iq-1)*ntr:iq*ntr)
            end do
         end do
      end do

      ! Clear any remaining messages
      sreq = nreq - rreq
      call START_LOG(mpiwaitdep_begin)
      call MPI_Waitall(sreq,ireq(rreq+1:nreq),MPI_STATUSES_IGNORE,ierr)
      call END_LOG(mpiwaitdep_end)

      call END_LOG(intssync_end)

   end subroutine intssync_recv

   subroutine indv_mpi(iq, i, j, n)
      integer , intent(in) :: iq
      integer , intent(out) :: i
      integer , intent(out) :: j
      integer , intent(out) :: n
      integer(kind=4) :: ierr

      ! Calculate local i, j, n from global iq

      ! Global i, j, n
      n = (iq - 1)/(il_g*il_g)
      j = 1 + (iq - n*il_g*il_g - 1)/il_g
      i = iq - (j - 1)*il_g - n*il_g*il_g
      if ( fproc(i,j,n) /= myid ) then
         write(*,"(a,5i5)") "Consistency failure in indv_mpi", myid, iq, i, j, n
         call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
      end if
      ! Reduced to values on my processor
      j = j - joff
      i = i - ioff
      n = n + noff      
   end subroutine indv_mpi

   function indglobal(i,j,n) result(iq)
      integer, intent(in) :: i, j, n
      integer :: iq

      ! Calculate a 1D global index from the global indices
      ! n in range 0:npanels
      iq = i + (j-1)*il_g + n*il_g*il_g
   end function indglobal

   function indg(i,j,n) result(iq)
      integer, intent(in) :: i, j, n
      integer :: iq

      ! Calculate a 1D global index from the local processors indices
      ! n in range 1..npan
      iq = i+ioff + (j+joff-1)*il_g + (n-noff)*il_g*il_g
   end function indg

   function indp(i,j,n) result(iq)
      integer, intent(in) :: i, j, n
      integer :: iq

      ! Calculate a 1D local index from the local processors indices
      ! Note that face number runs from 1 here.
      iq = i + (j-1)*ipan + (n-1)*ipan*jpan
   end function indp

   function iq2iqg(iq) result(iqg)
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

   function fproc(i,j,n) result(fpout)
      ! locates processor that owns a global grid point
      integer, intent(in) :: i, j, n
      integer :: fpout
      integer :: ip, jp

      ip = (i-1)/ipan
      jp = (j-1)/jpan
#ifdef uniform_decomp
      fpout = ip + jp*nxproc
#else
      fpout = ip + jp*nxproc + n*nxproc*nyproc/npan
#endif
   
   end function fproc

   function qproc(iqg) result(qpout)
      ! locates processor that owns a global grid point
      integer, intent(in) :: iqg
      integer :: qpout
      integer :: i, j, n

      n = (iqg - 1) / (il_g*il_g)
      j = 1 + (iqg - n*il_g*il_g - 1)/il_g
      i = iqg - (j - 1)*il_g - n*il_g*il_g

      qpout = fproc(i,j,n)
   
   end function qproc

   subroutine checksize(len, msize, mesg)
      integer, intent(in) :: len
      integer, intent(in) :: msize
      character(len=*), intent(in) :: mesg
      integer(kind=4) :: ierr
      if ( len > msize ) then
         write(6,*) "Error, maxsize exceeded in ", mesg
         call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
      end if
   end subroutine checksize

   subroutine check_bnds_alloc(rproc, iext)
      integer, intent(in) :: rproc
      integer, intent(in) :: iext
      integer :: len
      integer(kind=4) ierr

!     Allocate the components of the bnds array. It's too much work to
!     get the exact sizes, so allocate a fixed size for each case where
!     there's an interaction.
      if ( bnds(rproc)%len == 0 ) then
         ! Not allocated yet.
         len = maxbuflen
         if (rproc /= myid) then
            allocate ( bnds(rproc)%request_list(len) )
            allocate ( bnds(rproc)%send_list(len) )
            allocate ( bnds(rproc)%unpack_list(len) )
         end if
         allocate ( bnds(rproc)%request_list_uv(len) )
         allocate ( bnds(rproc)%send_list_uv(len) )
         allocate ( bnds(rproc)%unpack_list_uv(len) )
         allocate ( bnds(rproc)%uv_swap(len), bnds(rproc)%send_swap(len) )
         allocate ( bnds(rproc)%uv_neg(len), bnds(rproc)%send_neg(len) )
         bnds(rproc)%uv_neg = .false.
         bnds(rproc)%send_neg = 1.
         bnds(rproc)%len = len
      else
         ! Just check length
         if ( max(kl,ol)*bnds(rproc)%rlen >=  bnds(rproc)%len ) then
            write(6,*) "Error, maximum length error in check_bnds_alloc"
            write(6,*) myid, rproc, bnds(rproc)%rlen,  bnds(rproc)%len, max(kl,ol)
            call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
         end if
         if ( iext >= iextra ) then
            write(6,*) "Error, iext maximum length error in check_bnds_alloc"
            write(6,*) myid, iext, iextra
            call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
         end if
      end if
   end subroutine check_bnds_alloc

   !subroutine fix_index(iqq, larray, n, bnds, iext)
   !   integer, intent(in) :: iqq, n
   !   integer, dimension(:), intent(out) :: larray
   !   integer, intent(inout) :: iext
   !   type(bounds_info), dimension(0:), intent(inout) :: bnds
   !   integer :: rproc
   !   integer :: iloc,jloc,nloc
   !
   !   ! This processes extra corner points, so adds to rlenx
   !   ! Which processor has this point
   !   rproc = qproc(iqq)
   !   if ( rproc /= myid ) then ! Add to list
   !      call check_bnds_alloc(rproc, iext)
   !      bnds(rproc)%rlenx = bnds(rproc)%rlenx + 1
   !      bnds(rproc)%request_list(bnds(rproc)%rlenx) = iqq
   !      ! Increment extended region index
   !      iext = iext + 1
   !      bnds(rproc)%unpack_list(bnds(rproc)%rlenx) = iext
   !      larray(n) = ifull+iext
   !   else
   !      ! If it's on this processor, just use the local index
   !      call indv_mpi(iqq,iloc,jloc,nloc)
   !      larray(n) = indp(iloc,jloc,nloc)
   !   end if
   !end subroutine fix_index

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
         call check_bnds_alloc(rproc, iext)
         bnds(rproc)%rlen2 = bnds(rproc)%rlen2 + 1
         bnds(rproc)%request_list(bnds(rproc)%rlen2) = iqq
         ! Increment extended region index
         iext = iext + 1
         bnds(rproc)%unpack_list(bnds(rproc)%rlen2) = iext
         larray(n) = ifull+iext
      else
         ! If it's on this processor, just use the local index
         call indv_mpi(iqq,iloc,jloc,nloc)
         larray(n) = indp(iloc,jloc,nloc)
      end if
   end subroutine fix_index2

   subroutine proc_setup
      include 'parm.h'
!     Routine to set up offsets etc.
      integer :: i, j, n, nd, jdf, idjd_g
      integer, dimension(0:npanels) :: ipoff, jpoff

      call face_set( ipan, jpan, noff, ipoff, jpoff, npan, il_g, myid, nproc, nxproc, nyproc )
      ioff = ipoff(0)
      joff = jpoff(0)

!      ipfull = ipan*jpan*npan
!      iextra = 4*npan*(ipan+jpan+8)

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

#ifdef uniform_decomp
   subroutine proc_setup_uniform
      include 'parm.h'
!     Routine to set up offsets etc for the uniform decomposition
      integer :: i, j, n, nd, jdf, idjd_g
      integer, dimension(0:npanels) :: ipoff, jpoff
      integer(kind=4) ierr

      call dix_set( ipan, jpan, noff, ipoff, jpoff, npan, il_g, myid, nproc, nxproc, nyproc)
      ioff = ipoff(0)
      joff = jpoff(0)

!     Check that the values calculated here match those set as parameters
      if ( ipan /= il ) then
         write(6,*) "Error, parameter mismatch, ipan /= il", ipan, il
         call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
      end if
      if ( jpan*npan /= jl ) then
         write(6,*) "Error, parameter mismatch, jpan*npan /= jl", jpan, npan, jl
         call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
      end if

!      ipfull = ipan*jpan*npan
!      iextra = 4*npan*(ipan+jpan+8)

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
#endif

   subroutine face_set(ipan_l, jpan_l, noff_l, ioff_l, joff_l, npan_l, il_gx, myid_l, nproc_l, nxproc_l, nyproc_l)
      integer, intent(in) :: myid_l, nproc_l, npan_l, il_gx
      integer, intent(out) :: ipan_l, jpan_l, noff_l, nxproc_l, nyproc_l
      integer, dimension(0:npanels), intent(out) :: ioff_l, joff_l 
      integer n
      integer(kind=4) ierr

      !  Processor allocation
      !  if  nproc_l <= npanels+1, then each gets a number of full panels
      if ( nproc_l<=npanels+1 ) then
         if ( modulo( npanels+1, nproc_l )/=0 ) then
            write(6,*) "Error, number of processors must divide number of panels"
            call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
         end if
!         npan_l = (npanels+1)/nproc_l
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
            call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
         end if
!         npan_l = 1
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
            call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
         end if

         ! Still need to check that the processor distribution is compatible
         ! with the grid.
         if ( modulo( il_gx, nxproc_l )/=0 ) then
            write(6,*) "Error, il not a multiple of nxproc", il_gx, nxproc_l
            call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
         end if
         if ( modulo( il_gx, nyproc_l )/=0 ) then
            write(6,*) "Error, il not a multiple of nyproc", il_gx, nyproc_l
            call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
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
      integer(kind=4) ierr
      
      if ( npan_l /= npanels+1 ) then
         write(6,*) "Error: inconsistency in proc_setup_uniform"
         call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
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
         call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
      end if

      ! Still need to check that the processor distribution is compatible
      ! with the grid.
      if ( modulo( il_gx, nxproc_l )/=0 ) then
         write(6,*) "Error, il not a multiple of nxproc", il_gx, nxproc_l
         call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
      end if
      if ( modulo( il_gx, nyproc_l )/=0 ) then
         write(6,*) "Error, il not a multiple of nyproc", il_gx, nyproc_l
         call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
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
      integer(kind=4) ierr
      
      if ( npan_l /= npanels+1 ) then
         write(6,*) "Error: inconsistency in proc_setup_uniform"
         call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
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
         call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
      end if

      ! Still need to check that the processor distribution is compatible
      ! with the grid.
      if ( modulo( il_gx, nxproc_l )/=0 ) then
         write(6,*) "Error, il not a multiple of nxproc", il_gx, nxproc_l
         call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
      end if
      if ( modulo( il_gx, nyproc_l )/=0 ) then
         write(6,*) "Error, il not a multiple of nyproc", il_gx, nyproc_l
         call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
      end if
      ipan_l = il_gx/nxproc_l
      jpan_l = il_gx/nyproc_l

      ! Set offsets for this processor
      do n = 0,npanels
         call proc_region_uniform(myid_l,n,ioff_l(n),joff_l(n),noff_l,nxproc_l,nyproc_l,ipan_l,jpan_l)
      end do

   end subroutine uniform_set

   subroutine proc_region(procid,panid,ipoff,jpoff,npoff,nxproc_l,nyproc_l,ipan_l,jpan_l,npan_l,dmode)
      ! Calculate the offsets for a given processor
      integer, intent(in) :: procid, panid, nxproc_l, nyproc_l, ipan_l, jpan_l, npan_l
      integer, intent(in) :: dmode
      integer, intent(out) :: ipoff, jpoff, npoff

      select case(dmode)
         case(0) ! Face
            call proc_region_face(procid,ipoff,jpoff,npoff,nxproc_l,nyproc_l,ipan_l,jpan_l,npan_l)
      
         case(1) ! Old uniform
            call proc_region_uniform(procid,panid,ipoff,jpoff,npoff,nxproc_l,nyproc_l,ipan_l,jpan_l)
            
         case(2) ! New uniform
            call proc_region_dix(procid,ipoff,jpoff,npoff,nxproc_l,ipan_l,jpan_l)

#ifdef debug            
         case default
            write(6,*) "ERROR: Invalid decomposition ",dmode
            call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
#endif
      end select
     
   end subroutine proc_region

   subroutine proc_region_face(procid,ipoff,jpoff,npoff,nxproc_l,nyproc_l,ipan_l,jpan_l,npan_l)
      ! Calculate the offsets for a given processor
      integer, intent(in) :: procid, nxproc_l, nyproc_l, ipan_l, jpan_l, npan_l
      integer, intent(out) :: ipoff, jpoff, npoff
      integer :: myface, mtmp, nproc_l

      nproc_l = nxproc_l*nyproc_l*(npanels+1)/npan_l
      if ( nproc_l<=npanels+1 ) then
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

   subroutine proc_region_uniform(procid,panid,ipoff,jpoff,npoff,nxproc_l,nyproc_l,ipan_l,jpan_l)
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
            ipoff = modulo(procid,nxproc_l) * ipan_l
            if (jpoff>=il_gx/2) then
               ipoff=il_gx-ipoff-ipan_l
            end if
         case(1)
            jpoff = (procid/nxproc_l) * jpan_l
            ipoff = modulo(procid,nxproc_l) * ipan_l
         case(2)
            jpoff = (procid/nxproc_l) * jpan_l
            ipoff = modulo(procid,nxproc_l) * ipan_l
            if (ipoff>=il_gx/2) then
               jpoff=il_gx-jpoff-jpan_l
            end if
         case(3)
            jpoff = modulo(procid,nyproc_l) * jpan_l
            ipoff = (procid/nyproc_l) * ipan_l
            if (ipoff>=il_gx/2) then
               jpoff=il_gx-jpoff-jpan_l
            end if
         case(4)
            jpoff = modulo(procid,nyproc_l) * jpan_l
            ipoff = (procid/nyproc_l) * ipan_l
         case(5)
            jpoff = modulo(procid,nyproc_l) * jpan_l
            ipoff = (procid/nyproc_l) * ipan_l
            if (jpoff>=il_gx/2) then
               ipoff=il_gx-ipoff-ipan_l
            end if        
      end select
     
   end subroutine proc_region_uniform

   subroutine proc_region_dix(procid,ipoff,jpoff,npoff,nxproc_l,ipan_l,jpan_l)
      ! Calculate the offsets for a given processor
      integer, intent(in) :: procid, nxproc_l, ipan_l, jpan_l
      integer, intent(out) :: ipoff, jpoff, npoff

      ! Original Dix uniform decomposition
      ! Set offsets for this processor (same on all faces)
      npoff = 1
      jpoff = (procid/nxproc_l) * jpan_l
      ipoff = modulo(procid,nxproc_l)*ipan_l
     
   end subroutine proc_region_dix

   subroutine start_log ( event )
      integer, intent(in) :: event
#ifdef vampir
      VT_USER_START(event_name(event))
#endif
#ifdef simple_timer
      start_time(event) = MPI_Wtime()
#endif 
   end subroutine start_log

   subroutine end_log ( event )
      integer, intent(in) :: event
#ifdef vampir
      VT_USER_END(event_name(event))
#endif

#ifdef simple_timer
      tot_time(event) = tot_time(event) + MPI_Wtime() - start_time(event)
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

   subroutine log_setup()
#ifdef vampir
#ifdef simple_timer
      write(6,*) "ERROR: vampir and simple_timer should not be compiled together"
      stop
#endif
#endif

      model_begin = 1
      model_end =  model_begin
      event_name(model_begin) = "Model"

      maincalc_begin = 2
      maincalc_end =  maincalc_begin
      event_name(maincalc_begin) = "MainCalc"

      phys_begin = 3
      phys_end =  phys_begin
      event_name(phys_begin) = "Phys"

      physloadbal_begin = 4
      physloadbal_end =  physloadbal_begin
      event_name(physloadbal_begin) = "PhysLoadBal"

      nonlin_begin = 5
      nonlin_end = nonlin_begin 
      event_name(nonlin_begin) = "Nonlin"

      helm_begin = 6
      helm_end = helm_begin
      event_name(helm_begin) = "Helm"

      adjust_begin = 7
      adjust_end = adjust_begin
      event_name(Adjust_begin) = "Adjust"

      upglobal_begin = 8
      upglobal_end = upglobal_begin
      event_name(upglobal_begin) = "Upglobal"

      hordifg_begin = 9
      hordifg_end = hordifg_begin
      event_name(hordifg_begin) = "Hordifg"

      vadv_begin = 10
      vadv_end = vadv_begin
      event_name(vadv_begin) = "Vadv"

      depts_begin = 11
      depts_end = depts_begin
      event_name(depts_begin) = "Depts"

      ints_begin = 12
      ints_end = ints_begin 
      event_name(ints_begin) = "Ints"
      
      stag_begin = 13
      stag_end = stag_begin
      event_name(stag_begin) = "Stag"

      toij_begin = 14
      toij_end =  toij_begin
      event_name(toij_begin) = "Toij"

      outfile_begin = 15
      outfile_end = outfile_begin
      event_name(outfile_begin) = "Outfile"

      onthefly_begin = 16
      onthefly_end = onthefly_begin
      event_name(onthefly_begin) = "Onthefly"

      otf_ints1_begin = 17
      otf_ints1_end = otf_ints1_begin
      event_name(otf_ints1_begin) = "OTF_ints1"      

      otf_ints4_begin = 18
      otf_ints4_end = otf_ints4_begin
      event_name(otf_ints4_begin) = "OTF_ints4"       
      
      histrd1_begin = 19
      histrd1_end = histrd1_begin
      event_name(histrd1_begin) = "HistRd1"
      
      histrd4_begin = 20
      histrd4_end = histrd4_begin
      event_name(histrd4_begin) = "HistRd4"

       bounds_begin = 21
      bounds_end = bounds_begin
      event_name(bounds_begin) = "Bounds"

      boundsuv_begin = 22
      boundsuv_end = boundsuv_begin
      event_name(boundsuv_begin) = "BoundsUV"

      deptsync_begin = 23
      deptsync_end = deptsync_begin
      event_name(deptsync_begin) = "Deptsync"

      intssync_begin = 24
      intssync_end = intssync_begin
      event_name(intssync_begin) = "Intssync"

      gatherrma_begin = 25
      gatherrma_end = gatherrma_begin
      event_name(gatherrma_begin) = "GatherRMA"      
      
      gather_begin = 26
      gather_end = gather_begin
      event_name(gather_begin) = "Gather"

      distribute_begin = 27
      distribute_end = distribute_begin
      event_name(distribute_begin) = "Distribute"

      posneg_begin = 28
      posneg_end = posneg_begin
      event_name(posneg_begin) = "Posneg"

      globsum_begin = 29
      globsum_end = globsum_begin
      event_name(globsum_begin) = "Globsum"
      
      precon_begin = 30
      precon_end = precon_begin
      event_name(precon_begin) = "Precon"

      indata_begin = 31
      indata_end =  indata_begin
      event_name(indata_begin) = "Indata"

      nestin_begin = 32
      nestin_end =  nestin_begin
      event_name(nestin_begin) = "Nestin"
      
      nestpack_begin = 33
      nestpack_end =  nestpack_begin
      event_name(nestpack_begin) = "Nest_pack"

      nestcalc_begin = 34
      nestcalc_end =  nestcalc_begin
      event_name(nestcalc_begin) = "Nest_calc"

      nestcomm_begin = 35
      nestcomm_end =  nestcomm_begin
      event_name(nestcomm_begin) = "Nest_comm"
      
      nestunpack_begin = 36
      nestunpack_end =  nestunpack_begin
      event_name(nestunpack_begin) = "Nest_unpack"
      
      gwdrag_begin = 37
      gwdrag_end =  gwdrag_begin
      event_name(gwdrag_begin) = "GWdrag"

      convection_begin = 38
      convection_end =  convection_begin
      event_name(convection_begin) = "Convection"

      cloud_begin = 39
      cloud_end =  cloud_begin
      event_name(cloud_begin) = "Cloud"

      radnet_begin = 40
      radnet_end =  radnet_begin
      event_name(radnet_begin) = "Rad_net"

      radmisc_begin = 41
      radmisc_end =  radmisc_begin
      event_name(radmisc_begin) = "Rad_misc"
      
      radsw_begin = 42
      radsw_end =  radsw_begin
      event_name(radsw_begin) = "Rad_SW"

      radlw_begin = 43
      radlw_end =  radlw_begin
      event_name(radlw_begin) = "Rad_LW"      

      sfluxnet_begin = 44
      sfluxnet_end =  sfluxnet_begin
      event_name(sfluxnet_begin) = "Sflux_net"
      
      sfluxwater_begin = 45
      sfluxwater_end =  sfluxwater_begin
      event_name(sfluxwater_begin) = "Sflux_water"

      sfluxland_begin = 46
      sfluxland_end =  sfluxland_begin
      event_name(sfluxland_begin) = "Sflux_land"

      sfluxurban_begin = 47
      sfluxurban_end =  sfluxurban_begin
      event_name(sfluxurban_begin) = "Sflux_urban"

      vertmix_begin = 48
      vertmix_end =  vertmix_begin
      event_name(vertmix_begin) = "Vertmix"

      aerosol_begin = 49
      aerosol_end =  aerosol_begin
      event_name(aerosol_begin) = "Aerosol"

      waterdynamics_begin = 50
      waterdynamics_end =  waterdynamics_begin
      event_name(waterdynamics_begin) = "Waterdynamics"

      watermisc_begin = 51
      watermisc_end =  watermisc_begin
      event_name(watermisc_begin) = "Water_misc"

      waterdeps_begin = 52
      waterdeps_end =  waterdeps_begin
      event_name(waterdeps_begin) = "Water_deps"

      watereos_begin = 53
      watereos_end =  watereos_begin
      event_name(watereos_begin) = "Water_EOS"

      waterhadv_begin = 54
      waterhadv_end =  waterhadv_begin
      event_name(waterhadv_begin) = "Water_Hadv"

      watervadv_begin = 55
      watervadv_end =  watervadv_begin
      event_name(watervadv_begin) = "Water_Vadv"

      waterhelm_begin = 56
      waterhelm_end =  waterhelm_begin
      event_name(waterhelm_begin) = "Water_helm"

      wateriadv_begin = 57
      wateriadv_end =  wateriadv_begin
      event_name(wateriadv_begin) = "Water_Iadv"
      
      ocnstag_begin = 58
      ocnstag_end = ocnstag_begin
      event_name(ocnstag_begin) = "Water_stag"      

      waterdiff_begin = 59
      waterdiff_end =  waterdiff_begin
      event_name(waterdiff_begin) = "Waterdiff"

      river_begin = 60
      river_end =  river_begin
      event_name(river_begin) = "River"

      mgsetup_begin = 61
      mgsetup_end =  mgsetup_begin
      event_name(mgsetup_begin) = "MG_Setup"

      mgdecomp_begin = 62
      mgdecomp_end =  mgdecomp_begin
      event_name(mgdecomp_begin) = "MG_Decomp"
      
      mgfine_begin = 63
      mgfine_end =  mgfine_begin
      event_name(mgfine_begin) = "MG_Fine"

      mgup_begin = 64
      mgup_end =  mgup_begin
      event_name(mgup_begin) = "MG_Up"

      mgcoarse_begin = 65
      mgcoarse_end =  mgcoarse_begin
      event_name(mgcoarse_begin) = "MG_Coarse"

      mgdown_begin = 66
      mgdown_end = mgdown_begin
      event_name(mgdown_begin) = "MG_Down"

      mgmlosetup_begin = 67
      mgmlosetup_end = mgmlosetup_begin
      event_name(mgmlosetup_begin) = "MGMLO_Setup"

      mgmlodecomp_begin = 68
      mgmlodecomp_end = mgmlodecomp_begin
      event_name(mgmlodecomp_begin) = "MGMLO_Decomp"      
      
      mgmlofine_begin = 69
      mgmlofine_end = mgmlofine_begin
      event_name(mgmlofine_begin) = "MGMLO_Fine"

      mgmloup_begin = 70
      mgmloup_end = mgmloup_begin
      event_name(mgmloup_begin) = "MGMLO_Up"

      mgmlocoarse_begin = 71
      mgmlocoarse_end = mgmlocoarse_begin
      event_name(mgmlocoarse_begin) = "MGMLO_Coarse"

      mgmlodown_begin = 72
      mgmlodown_end = mgmlodown_begin
      event_name(mgmlodown_begin) = "MGMLO_Down"

      mgbounds_begin = 73
      mgbounds_end = mgbounds_begin
      event_name(mgbounds_begin) = "MG_bounds"
      
      mgcollect_begin = 74
      mgcollect_end = mgcollect_begin
      event_name(mgcollect_begin) = "MG_collect"      

      mgbcast_begin = 75
      mgbcast_end = mgbcast_begin
      event_name(mgbcast_begin) = "MG_bcast"   

      bcast_begin = 76
      bcast_end = bcast_begin
      event_name(bcast_begin) = "MPI_Bcast"

      allgatherx_begin = 77
      allgatherx_end = allgatherx_begin
      event_name(allgatherx_begin) = "MPI_AllGather" 
      
      gatherx_begin = 78
      gatherx_end = gatherx_begin
      event_name(gatherx_begin) = "MPI_Gather"

      scatterx_begin = 79
      scatterx_end = scatterx_begin
      event_name(scatterx_begin) = "MPI_Scatter"
      
      reduce_begin = 80
      reduce_end = reduce_begin
      event_name(reduce_begin) = "MPI_Reduce"
      
      mpiwait_begin = 81
      mpiwait_end = mpiwait_begin
      event_name(mpiwait_begin) = "MPI_Wait"

      mpiwaituv_begin = 82
      mpiwaituv_end = mpiwaituv_begin
      event_name(mpiwaituv_begin) = "MPI_WaitUV"

      mpiwaituvtile_begin = 83
      mpiwaituvtile_end = mpiwaituvtile_begin
      event_name(mpiwaituvtile_begin) = "MPI_WaitUV_Tile"

      mpiwaitdep_begin = 84
      mpiwaitdep_end = mpiwaitdep_begin
      event_name(mpiwaitdep_begin) = "MPI_WaitDEP"

      mpiwaitmg_begin = 85
      mpiwaitmg_end = mpiwaitmg_begin
      event_name(mpiwaitmg_begin) = "MPI_WaitMG"
     
   end subroutine log_setup
   
   subroutine phys_loadbal()
!     This forces a sychronisation to make the physics load imbalance overhead
!     explicit. 
      integer(kind=4) :: ierr

      call START_LOG(physloadbal_begin)
      call MPI_Barrier( MPI_COMM_WORLD, ierr )
      call END_LOG(physloadbal_end)
   end subroutine phys_loadbal

#ifdef simple_timer
   subroutine simple_timer_finalize()
      ! Calculate the mean, min and max times for each case
      integer :: i
      integer(kind=4) :: ierr, llen
      real(kind=8), dimension(nevents) :: emean, emax, emin
      llen=nevents
      call MPI_Reduce(tot_time, emean, llen, MPI_DOUBLE_PRECISION, &
                      MPI_SUM, 0_4, MPI_COMM_WORLD, ierr )
      call MPI_Reduce(tot_time, emax, llen, MPI_DOUBLE_PRECISION, &
                      MPI_MAX, 0_4, MPI_COMM_WORLD, ierr )
      call MPI_Reduce(tot_time, emin, llen, MPI_DOUBLE_PRECISION, &
                      MPI_MIN, 0_4, MPI_COMM_WORLD, ierr )
      if ( myid == 0 ) then
         write(6,*) "==============================================="
         write(6,*) "  Times over all processes"
         write(6,*) "  Routine        Mean time  Min time  Max time"
         do i=1,nevents
            if ( emean(i) > 0. ) then
               ! This stops boundsa, b getting written when they're not used.
               write(*,"(a,3f10.3)") event_name(i), emean(i)/nproc, emin(i), emax(i)
            end if
         end do
      end if
         
   end subroutine simple_timer_finalize
#endif

    subroutine ccglobal_posneg2 (array, delpos, delneg)
       ! Calculate global sums of positive and negative values of array
       use sumdd_m
       use xyzinfo_m       

       real, intent(in), dimension(ifull) :: array
       real, intent(out) :: delpos, delneg
       integer(kind=4) :: ierr
#ifdef i8r8
       integer(kind=4), parameter :: ltype = MPI_DOUBLE_COMPLEX
#else
       integer(kind=4), parameter :: ltype = MPI_COMPLEX
#endif   
       complex, dimension(2) :: local_sum, global_sum
!      Temporary array for the drpdr_local function
       real, dimension(ifull) :: tmparr 

       call START_LOG(posneg_begin)

       local_sum(1:2) = cmplx(0., 0.)
       tmparr(1:ifull)  = max(0., array(1:ifull)*wts(1:ifull))
       call drpdr_local(tmparr, local_sum(1))
       tmparr(1:ifull)  = min(0., array(1:ifull)*wts(1:ifull))
       call drpdr_local(tmparr, local_sum(2))
       global_sum(1:2) = cmplx(0., 0.)
       call MPI_Allreduce ( local_sum, global_sum, 2_4, ltype, MPI_SUMDR, MPI_COMM_WORLD, ierr )
       delpos = real(global_sum(1))
       delneg = real(global_sum(2))

       call END_LOG(posneg_end)

    end subroutine ccglobal_posneg2
    
    subroutine ccglobal_posneg3 (array, delpos, delneg, dsigin)
       ! Calculate global sums of positive and negative values of array
       use sigs_m
       use sumdd_m
       use xyzinfo_m
       real, intent(in), dimension(:,:) :: array
       real, intent(in), dimension(:), optional :: dsigin
       real, intent(out) :: delpos, delneg
       real, dimension(size(array,2)) :: dsigx
       integer :: k, kx
       integer(kind=4) :: ierr
#ifdef i8r8
       integer(kind=4), parameter :: ltype = MPI_DOUBLE_COMPLEX
#else
       integer(kind=4), parameter :: ltype = MPI_COMPLEX
#endif   
       complex, dimension(2) :: local_sum, global_sum
!      Temporary array for the drpdr_local function
       real, dimension(ifull) :: tmparr

       call START_LOG(posneg_begin)

       kx = size(array,2)
       if (present(dsigin)) then
         dsigx(1:kx) = -dsigin(1:kx)
       else
         dsigx(1:kx) = dsig(1:kx)
       end if
       
       local_sum(1:2) = cmplx(0., 0.)
       do k=1,kx
          tmparr(1:ifull) = max(0., -dsigx(k)*array(1:ifull,k)*wts(1:ifull))
          call drpdr_local(tmparr, local_sum(1))
          tmparr(1:ifull) = min(0., -dsigx(k)*array(1:ifull,k)*wts(1:ifull))
          call drpdr_local(tmparr, local_sum(2))
       end do ! k loop
       global_sum(1:2) = cmplx(0., 0.)
       call MPI_Allreduce ( local_sum, global_sum, 2_4, ltype, MPI_SUMDR, MPI_COMM_WORLD, ierr )
       delpos = real(global_sum(1))
       delneg = real(global_sum(2))

       call END_LOG(posneg_end)

    end subroutine ccglobal_posneg3

    subroutine ccglobal_posneg4 (array, delpos, delneg, dsigin)
       ! Calculate global sums of positive and negative values of array
       use sigs_m
       use sumdd_m
       use xyzinfo_m
       real, intent(in), dimension(:,:,:) :: array
       real, intent(in), dimension(:), optional :: dsigin
       real, intent(out), dimension(:) :: delpos, delneg
       real, dimension(size(array,2)) :: dsigx
       integer :: i, k, kx, ntr
       integer(kind=4) :: ierr, mnum
#ifdef i8r8
       integer(kind=4), parameter :: ltype = MPI_DOUBLE_COMPLEX
#else
       integer(kind=4), parameter :: ltype = MPI_COMPLEX
#endif   
       complex, dimension(2*size(array,3)) :: local_sum, global_sum
!      Temporary array for the drpdr_local function
       real, dimension(ifull) :: tmparr

       call START_LOG(posneg_begin)

       kx  = size(array,2)
       ntr = size(array,3)
       if (present(dsigin)) then
         dsigx(1:kx) = -dsigin(1:kx)
       else
         dsigx(1:kx) = dsig(1:kx)
       end if

       local_sum(1:2*ntr) = cmplx(0.,0.)
       do i = 1,ntr
          do k=1,kx
             tmparr(1:ifull) = max(0.,-dsigx(k)*array(1:ifull,k,i)*wts(1:ifull))
             call drpdr_local(tmparr, local_sum(i))
             tmparr(1:ifull) = min(0.,-dsigx(k)*array(1:ifull,k,i)*wts(1:ifull))
             call drpdr_local(tmparr, local_sum(i+ntr))
          end do ! k loop
       end do
       mnum = 2*ntr
       global_sum(1:2*ntr) = cmplx(0.,0.)
       call MPI_Allreduce ( local_sum, global_sum, mnum, ltype, MPI_SUMDR, MPI_COMM_WORLD, ierr )
       delpos(1:ntr) = real(global_sum(1:ntr))
       delneg(1:ntr) = real(global_sum(ntr+1:2*ntr))

       call END_LOG(posneg_end)

    end subroutine ccglobal_posneg4

    subroutine ccglobal_sum2 (array, result)
       ! Calculate global sum of an array
       use sumdd_m
       use xyzinfo_m
       real, intent(in), dimension(ifull) :: array
       real, intent(out) :: result
       real :: result_l
       integer :: iq
       integer(kind=4) :: ierr
#ifdef i8r8
       integer(kind=4), parameter :: ltype = MPI_DOUBLE_COMPLEX
#else
       integer(kind=4), parameter :: ltype = MPI_COMPLEX
#endif   
       complex :: local_sum, global_sum
!      Temporary array for the drpdr_local function
       real, dimension(ifull) :: tmparr

       call START_LOG(globsum_begin)

       result_l = 0.
       do iq = 1,ifull
          tmparr(iq)  = array(iq)*wts(iq)
       enddo
       local_sum = cmplx(0.,0.)
       call drpdr_local(tmparr, local_sum)
       call MPI_Allreduce ( local_sum, global_sum, 1_4, ltype, MPI_SUMDR, MPI_COMM_WORLD, ierr )
       result = real(global_sum)

       call END_LOG(globsum_end)

    end subroutine ccglobal_sum2

    subroutine ccglobal_sum3 (array, result, dsigin)
       ! Calculate global sum of 3D array, appyling vertical weighting
       use sigs_m
       use sumdd_m
       use xyzinfo_m
       real, intent(in), dimension(:,:) :: array
       real, intent(in), dimension(:), optional :: dsigin
       real, intent(out) :: result
       real, dimension(size(array,2)) :: dsigx
       real :: result_l
       integer :: k, iq, kx
       integer(kind=4) ierr
#ifdef i8r8
       integer(kind=4), parameter :: ltype = MPI_DOUBLE_COMPLEX
#else
       integer(kind=4), parameter :: ltype = MPI_COMPLEX
#endif   
       complex :: local_sum, global_sum
!      Temporary array for the drpdr_local function
       real, dimension(ifull) :: tmparr

       call START_LOG(globsum_begin)

       kx = size(array,2)
       if (present(dsigin)) then
         dsigx(1:kx) = -dsigin(1:kx)
       else
         dsigx(1:kx) = dsig(1:kx)
       end if
       
       result_l = 0.
       local_sum = cmplx(0.,0.)
       do k = 1,kx
          do iq = 1,ifull
             tmparr(iq)  = -dsigx(k)*array(iq,k)*wts(iq)
          enddo
          call drpdr_local(tmparr, local_sum)
       end do ! k
       call MPI_Allreduce ( local_sum, global_sum, 1_4, ltype, MPI_SUMDR, MPI_COMM_WORLD, ierr )
       result = real(global_sum)

       call END_LOG(globsum_end)

    end subroutine ccglobal_sum3

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
            ! Use explicit ranges here because some arguments might be extended.
            ! ccmpi_distribute3 expects kl, it's not general
            if ( kk == kl ) then
               call ccmpi_distribute(var(1:ifull,:),varg)
            else
               do k=1,kk
                  call ccmpi_distribute(var(1:ifull,k),varg(:,k))
               end do
            end if
         else
            if ( kk == kl ) then
               call ccmpi_distribute(var(1:ifull,:))
            else
               do k=1,kk
                  call ccmpi_distribute(var(1:ifull,k))
               end do
            end if
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
         ! Use explicit ranges here because some arguments might be extended.
         ! ccmpi_gather3 expects kl, it's not general
         if ( kk == kl ) then
            call ccmpi_gather(var(1:ifull,:),varg)
         else
            do k=1,kk
               call ccmpi_gather(var(1:ifull,k),varg(:,k))
            end do
         end if
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
         if ( kk == kl ) then
            call ccmpi_gather(var(1:ifull,:))
         else
            do k=1,kk
               call ccmpi_gather(var(1:ifull,k))
            end do
         end if
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
      
      call START_LOG(reduce_begin)
      
      select case( op )
         case( "max" )
            lop = MPI_MAX
         case( "min" )
            lop = MPI_MIN
         case( "sum" )
            lop = MPI_SUM
         case default
            write(6,*) "ERROR: Unknown option for ccmpi_reduce ",op
            call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
      end select
      
      lhost = host
      lcomm = comm
      lsize = size(ldat)

      call MPI_Reduce(ldat, gdat, lsize, ltype, lop, lhost, lcomm, ierr )
 
      call END_LOG(reduce_end)
   
   end subroutine ccmpi_reduce2i

   subroutine ccmpi_reduce1r(ldat,gdat,op,host,comm)
   
      integer, intent(in) :: host, comm
      integer(kind=4) :: ltype, lop, lcomm, lerr, lhost
      real, intent(in) :: ldat
      real, intent(out) :: gdat
      character(len=*), intent(in) :: op
      
      call START_LOG(reduce_begin)
      
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
            call MPI_Abort(MPI_COMM_WORLD,-1_4,lerr)
      end select
     
      call MPI_Reduce(ldat, gdat, 1_4, ltype, lop, lhost, lcomm, lerr )
   
      call END_LOG(reduce_end)
   
   end subroutine ccmpi_reduce1r
   
   subroutine ccmpi_reduce2r(ldat,gdat,op,host,comm)
   
      integer, intent(in) :: host, comm
      integer(kind=4) :: ltype, lop, lcomm, lerr, lsize, lhost
      real, dimension(:), intent(in) :: ldat
      real, dimension(:), intent(out) :: gdat
      character(len=*), intent(in) :: op
      
      call START_LOG(reduce_begin)
      
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
            call MPI_Abort(MPI_COMM_WORLD,-1_4,lerr)
      end select
     
      call MPI_Reduce(ldat, gdat, lsize, ltype, lop, lhost, lcomm, lerr )
   
      call END_LOG(reduce_end)
   
   end subroutine ccmpi_reduce2r

   subroutine ccmpi_reduce3r(ldat,gdat,op,host,comm)
   
      integer, intent(in) :: host, comm
      integer(kind=4) ltype, lop, lcomm, lerr, lsize, lhost
      real, dimension(:,:), intent(in) :: ldat
      real, dimension(:,:), intent(out) :: gdat
      character(len=*), intent(in) :: op
      
      call START_LOG(reduce_begin)

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
            call MPI_Abort(MPI_COMM_WORLD,-1_4,lerr)
      end select
      
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
      
      call START_LOG(reduce_begin)
      
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
            call MPI_Abort(MPI_COMM_WORLD,-1_4,lerr)
      end select
      
      lhost = host
      lcomm = comm
      lsize = size(ldat)
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
      
      call START_LOG(reduce_begin)
      
      select case( op )
         case( "or" )
            lop = MPI_LOR
         case( "and" )
            lop = MPI_LAND
         case default
            write(6,*) "ERROR: Unknown option for ccmpi_reduce ",op
            call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
      end select
      
      lhost = host
      lcomm = comm
      lsize = size(ldat)

      call MPI_Reduce(ldat, gdat, lsize, ltype, lop, lhost, lcomm, ierr )
 
      call END_LOG(reduce_end)
   
   end subroutine ccmpi_reduce2l
   
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
      
      call START_LOG(reduce_begin)
      
      select case( op )
         case( "max" )
            lop = MPI_MAX
         case( "min" )
            lop = MPI_MIN
         case( "sum" )
            lop = MPI_SUM
         case default
            write(6,*) "ERROR: Unknown option for ccmpi_allreduce ",op
            call MPI_Abort(MPI_COMM_WORLD,-1_4,lerr)
      end select
      
      lcomm = comm
      call MPI_AllReduce(ldat, gdat, 1_4, ltype, lop, lcomm, lerr )
 
      call END_LOG(reduce_end)
   
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
      
      call START_LOG(reduce_begin)
      
      select case( op )
         case( "max" )
            lop = MPI_MAX
         case( "min" )
            lop = MPI_MIN
         case( "sum" )
            lop = MPI_SUM
         case default
            write(6,*) "ERROR: Unknown option for ccmpi_allreduce ",op
            call MPI_Abort(MPI_COMM_WORLD,-1_4,lerr)
      end select
      
      lcomm = comm
      lsize = size(ldat)
      call MPI_AllReduce(ldat, gdat, lsize, ltype, lop, lcomm, lerr )
 
      call END_LOG(reduce_end)
   
   end subroutine ccmpi_allreduce2i

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
      
      call START_LOG(reduce_begin)
      
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
            call MPI_Abort(MPI_COMM_WORLD,-1_4,lerr)
      end select
     
      call MPI_AllReduce(ldat, gdat, lsize, ltype, lop, lcomm, lerr )
   
      call END_LOG(reduce_end)
   
   end subroutine ccmpi_allreduce2r
  
   subroutine ccmpi_allreduce3r(ldat,gdat,op,comm)
   
      use sumdd_m
   
      integer, intent(in) :: comm
      integer(kind=4) ltype, lop, lcomm, lerr, lsize
      real, dimension(:,:), intent(in) :: ldat
      real, dimension(:,:), intent(out) :: gdat
      character(len=*), intent(in) :: op
      
      call START_LOG(reduce_begin)
      
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
            call MPI_Abort(MPI_COMM_WORLD,-1_4,lerr)
      end select
      
      call MPI_AllReduce(ldat, gdat, lsize, ltype, lop, lcomm, lerr )
   
      call END_LOG(reduce_end)
   
   end subroutine ccmpi_allreduce3r
   
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
      
      call START_LOG(reduce_begin)
      
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
            call MPI_Abort(MPI_COMM_WORLD,-1_4,lerr)
      end select
      
      lcomm = comm
      lsize = size(ldat)
      call MPI_AllReduce(ldat, gdat, lsize, ltype, lop, lcomm, lerr )
   
      call END_LOG(reduce_end)
   
   end subroutine ccmpi_allreduce2c
   
   subroutine ccmpi_abort(ierrin)
   
      integer, intent(in) :: ierrin
      integer(kind=4) :: lerrin, ierr
      
      lerrin = ierrin
      call MPI_Abort(MPI_COMM_WORLD,lerrin,ierr)
   
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

      call START_LOG(bcast_begin)

      lhost = host
      lcomm = comm
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

      call START_LOG(bcast_begin)

      lhost = host
      lcomm = comm
      lsize = size(ldat)
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

      call START_LOG(bcast_begin)

      lhost = host
      lcomm = comm
      lsize = size(ldat)
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

      call START_LOG(bcast_begin)

      lhost = host
      lcomm = comm
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

      call START_LOG(bcast_begin)

      lhost = host
      lcomm = comm
      lsize = size(ldat)
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

      call START_LOG(bcast_begin)

      lhost = host
      lcomm = comm
      lsize = size(ldat)
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

      call START_LOG(bcast_begin)

      lhost = host
      lcomm = comm
      lsize = size(ldat)
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

      call START_LOG(bcast_begin)

      lhost = host
      lcomm = comm
      lsize = size(ldat)
      call MPI_Bcast(ldat,lsize,ltype,lhost,lcomm,lerr)
   
      call END_LOG(bcast_end)
   
   end subroutine ccmpi_bcast5r

   subroutine ccmpi_bcast1s(ldat,host,comm)
   
      integer, intent(in) :: host, comm
      integer(kind=4) :: lcomm, lhost, lerr, lsize
      character(len=*), intent(inout) :: ldat
      integer, parameter :: maxdummysize = 256
      integer i
      integer(kind=1), dimension(maxdummysize) :: dummy

      ! MJT notes - Windows MPI_CHARACTER seems broken
      
      call START_LOG(bcast_begin)

      lhost = host
      lcomm = comm
      lsize = len(ldat)
      if ( lsize > maxdummysize ) then
        write(6,*) "ERROR: Dummy array too small in ccmpi_bcast1s"
        call mpi_abort(MPI_COMM_WORLD,-1_4,lerr)
      end if
      do i = 1,lsize
         dummy(i) = int(iachar(ldat(i:i)),1)
      end do
      call MPI_Bcast(dummy,lsize,MPI_BYTE,lhost,lcomm,lerr)
      do i = 1,lsize
         ldat(i:i) = achar(dummy(i))
      end do
   
      call END_LOG(bcast_end)
   
   end subroutine ccmpi_bcast1s
   
   subroutine ccmpi_bcast2r8(ldat,host,comm)
   
      integer, intent(in) :: host, comm
      integer(kind=4) :: lcomm, lhost, ierr, lsize
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
      real(kind=8), dimension(:), intent(inout) :: ldat
   
      call START_LOG(bcast_end)
      
      lhost = host
      lcomm = comm
      lsize = size(ldat)
      call MPI_Bcast(ldat,lsize,ltype,lhost,lcomm,ierr)
   
      call END_LOG(bcast_end)
   
   end subroutine ccmpi_bcast2r8
   
   subroutine ccmpi_bcast3r8(ldat,host,comm)
   
      integer, intent(in) :: host,comm
      integer(kind=4) :: lcomm, lhost, ierr, lsize
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
      real(kind=8), dimension(:,:), intent(inout) :: ldat
   
      call START_LOG(bcast_begin)
      
      lhost = host
      lcomm = comm
      lsize = size(ldat)
      call MPI_Bcast(ldat, lsize, ltype, lhost, lcomm, ierr)
   
      call END_LOG(bcast_end)
   
   end subroutine ccmpi_bcast3r8   
   
   subroutine ccmpi_bcast4r8(ldat,host,comm)
   
      integer, intent(in) :: host,comm
      integer(kind=4) :: lcomm, lhost, ierr, lsize
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
      real(kind=8), dimension(:,:,:), intent(inout) :: ldat

      call START_LOG(bcast_begin)

      lhost = host
      lcomm = comm
      lsize = size(ldat)
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

      call START_LOG(gatherx_begin)
      
      lcomm = comm
      lhost = host
      lsize = size(ldat)
      call MPI_Gather(ldat,lsize,ltype,gdat,lsize,ltype,lhost,lcomm,lerr)
   
      call END_LOG(gatherx_end)
      
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

      call START_LOG(gatherx_begin)
      
      lcomm = comm
      lhost = host
      lsize = size(ldat)
      call MPI_Gather(ldat,lsize,ltype,gdat,lsize,ltype,lhost,lcomm,lerr)
   
      call END_LOG(gatherx_end)
      
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

      call START_LOG(gatherx_begin)
      
      lcomm = comm
      lhost = host
      lsize = size(ldat)
      call MPI_Gather(ldat,lsize,ltype,gdat,lsize,ltype,lhost,lcomm,lerr)
   
      call END_LOG(gatherx_end)
      
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

      call START_LOG(gatherx_begin)
      
      lcomm = comm
      lhost = host
      lsize = size(ldat)
      call MPI_Gather(ldat,lsize,ltype,gdat,lsize,ltype,lhost,lcomm,lerr)
   
      call END_LOG(gatherx_end)
      
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

      call START_LOG(gatherx_begin)
      
      lcomm = comm
      lhost = host
      lsize = size(ldat)
      call MPI_Gather(ldat,lsize,ltype,gdat,lsize,ltype,lhost,lcomm,lerr)
   
      call END_LOG(gatherx_end)
      
   end subroutine ccmpi_gatherx34r

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

      call START_LOG(gatherx_begin)

      lcomm = comm
      lhost = host
      lsize = size(ldat)
      call MPI_Gather(ldat,lsize,ltype,gdat,lsize,ltype,lhost,lcomm,lerr)

      call END_LOG(gatherx_end)
      
   end subroutine ccmpi_gatherx3i
   
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

      call START_LOG(scatterx_begin)
     
      lcomm = comm
      lhost = host
      lsize = size(ldat)
      call MPI_Scatter(gdat,lsize,ltype,ldat,lsize,ltype,lhost,lcomm,lerr)
   
      call END_LOG(scatterx_end)
      
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

      call START_LOG(scatterx_begin)
      
      lcomm = comm
      lhost = host
      lsize = size(ldat)
      call MPI_Scatter(gdat,lsize,ltype,ldat,lsize,ltype,lhost,lcomm,lerr)
   
      call END_LOG(scatterx_end)
      
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

      call START_LOG(scatterx_begin)
      
      lcomm = comm
      lhost = host
      lsize = size(ldat)
      call MPI_Scatter(gdat,lsize,ltype,ldat,lsize,ltype,lhost,lcomm,lerr)
   
      call END_LOG(scatterx_end)
      
   end subroutine ccmpi_scatterx3r
   
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
   
      call START_LOG(allgatherx_begin)
      
      lcomm = comm
      lsize = size(ldat)
      call MPI_AllGather(ldat,lsize,ltype,gdat,lsize,ltype,lcomm,lerr)
      
      call END_LOG(allgatherx_end)
      
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
   
      call START_LOG(allgatherx_begin)
      
      lcomm = comm
      lsize = size(ldat)
      call MPI_AllGather(ldat,lsize,ltype,gdat,lsize,ltype,lcomm,lerr)
      
      call END_LOG(allgatherx_end)
      
   end subroutine ccmpi_allgatherx2r
   
   subroutine ccmpi_init

      integer(kind=4) :: lerr, lproc, lid
#ifdef usempi3
      integer(kind=4) :: lcommout
      integer(kind=4) :: lcolour
#endif

      ! Global communicator
      call MPI_Init(lerr)
      call MPI_Comm_size(MPI_COMM_WORLD, lproc, lerr) ! Find number of processes
      call MPI_Comm_rank(MPI_COMM_WORLD, lid, lerr)   ! Find local processor id
      nproc      = lproc
      myid       = lid
      comm_world = MPI_COMM_WORLD
      
#ifdef usempi3
      ! Intra-node communicator
      lid = myid
      call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, lid, MPI_INFO_NULL, lcommout, lerr)
      call MPI_Comm_size(lcommout, lproc, lerr) ! Find number of processes on node
      call MPI_Comm_rank(lcommout, lid, lerr)   ! Find local processor id on node
      comm_node  = lcommout
      node_nproc = lproc
      node_myid  = lid
      
      ! Inter-node commuicator
      lcolour = node_myid
      lid = myid
      call MPI_Comm_Split(MPI_COMM_WORLD, lcolour, lid, lcommout, lerr)
      call MPI_Comm_size(lcommout, lproc, lerr) ! Find number of processes on node
      call MPI_Comm_rank(lcommout, lid, lerr)   ! Find local processor id on node
      comm_nodecaptian  = lcommout
      nodecaptian_nproc = lproc
      nodecaptian_myid  = lid
      
      if ( myid==0 .and. (node_myid/=0.or.nodecaptian_myid/=0) ) then
         write(6,*) "ERROR: Intra-node communicator failed"
         write(6,*) "myid, node_myid, nodecaptian_myid ",myid,node_myid,nodecaptian_myid
         call MPI_ABORT(MPI_COMM_WORLD, -1_4, lerr)
      end if
#endif

   end subroutine ccmpi_init
   
   subroutine ccmpi_finalize
   
      integer(kind=4) :: lerr
   
      if ( nproc>1 ) then
         call MPI_Win_free(localwin, lerr)
      end if
      call MPI_Finalize(lerr)
   
   end subroutine ccmpi_finalize

   subroutine ccmpi_commsplit(commout,comm,colour,rank)
   
      integer, intent(out) :: commout
      integer, intent(in) :: comm, colour, rank
      integer(kind=4) :: lcomm, lcommout, lerr, lrank, lcolour
   
      lcomm = comm
      lrank = rank
      if ( colour>=0 ) then
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

   ! This routine allows multi-grid bounds updates
   ! The code is based on cc_mpi bounds routines, but
   ! accomodates the g-th multi-grid
   
   subroutine mgbounds(g,vdat,klim,corner)

      integer, intent(in) :: g
      integer, intent(in), optional :: klim
      integer :: kx
      integer :: iproc, recv_len, send_len
      integer :: rcount, myrlen, jproc
      integer, dimension(mg(g)%neighnum) :: rslen, sslen
      integer(kind=4) :: ierr, itag=20, llen, sreq, lproc
      integer(kind=4) :: ldone
      integer(kind=4), dimension(size(ireq)) :: donelist
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      real, dimension(:,:), intent(inout) :: vdat
      logical, intent(in), optional :: corner
      logical extra

      call START_LOG(mgbounds_begin)
      
      if (present(klim)) then
         kx = klim
      else
         kx = size(vdat,2)
      end if
      if (present(corner)) then
         extra = corner
      else
         extra = .false.
      end if

      if ( extra ) then
         rslen  = mg_bnds(mg(g)%neighlist,g)%rlenx
         sslen  = mg_bnds(mg(g)%neighlist,g)%slenx
         myrlen = mg_bnds(myid,g)%rlenx
      else
         rslen  = mg_bnds(mg(g)%neighlist,g)%rlen
         sslen  = mg_bnds(mg(g)%neighlist,g)%slen
         myrlen = mg_bnds(myid,g)%rlen
      end if

      !     Set up the buffers to send and recv
      nreq = 0
      do iproc = 1,mg(g)%neighnum
         recv_len = rslen(iproc)
         if ( recv_len>0 ) then
            lproc = mg(g)%neighlist(iproc)  ! Recv from
            nreq = nreq + 1
            rlist(nreq) = iproc
            llen = recv_len*kx
            call MPI_IRecv( bnds(lproc)%rbuf(1), llen, ltype, lproc, &
                            itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do
      rreq = nreq
      do iproc = mg(g)%neighnum,1,-1
         send_len = sslen(iproc)
         if ( send_len>0 ) then
            lproc = mg(g)%neighlist(iproc)  ! Send to
            bnds(lproc)%sbuf(1:send_len*kx) = reshape( vdat(mg_bnds(lproc,g)%send_list(1:send_len),1:kx), (/ send_len*kx /) )
            nreq = nreq + 1
            llen = send_len*kx
            call MPI_ISend( bnds(lproc)%sbuf(1), llen, ltype, lproc, &
                            itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do

      ! Finally see if there are any points on my own processor that need
      ! to be fixed up. This will only be in the case when nproc < npanels.
      if ( myrlen>0 ) then
        vdat(mg(g)%ifull+mg_bnds(myid,g)%unpack_list(1:myrlen),1:kx) = vdat(mg_bnds(myid,g)%request_list(1:myrlen),1:kx)
      end if

      rcount = rreq
      do while ( rcount>0 )
         call START_LOG(mpiwaitmg_begin)
         call MPI_Waitsome( rreq, ireq, ldone, donelist, MPI_STATUSES_IGNORE, ierr )
         call END_LOG(mpiwaitmg_end)
         rcount = rcount - ldone
         do jproc = 1,ldone
            iproc = rlist(donelist(jproc))  ! Recv from
            lproc = mg(g)%neighlist(iproc)
            vdat(mg(g)%ifull+mg_bnds(lproc,g)%unpack_list(1:rslen(iproc)),1:kx) &
                = reshape( bnds(lproc)%rbuf(1:rslen(iproc)*kx), (/ rslen(iproc), kx /) )
         end do
      end do

      ! clear any remaining messages
      sreq = nreq - rreq
      call START_LOG(mpiwaitmg_begin)
      call MPI_Waitall(sreq,ireq(rreq+1:nreq),MPI_STATUSES_IGNORE,ierr)
      call END_LOG(mpiwaitmg_end)

      call END_LOG(mgbounds_end)

   return
   end subroutine mgbounds

   ! This subroutine merges datasets when upscaling with the multi-grid solver
   subroutine mgcollectreduce(g,vdat,dsolmax,klim)

      integer, intent(in) :: g
      integer, intent(in), optional :: klim
      integer :: kx, msg_len
      real, dimension(:,:), intent(inout) :: vdat
      real, dimension(:), intent(inout) :: dsolmax

      ! merge length
      if ( mg(g)%merge_len<=1 ) return

      call START_LOG(mgcollect_begin)

      if ( present(klim) ) then
         kx = klim
      else
         kx = size(vdat,2)      
      end if

      msg_len = mg(g)%ifull/(mg(g)%merge_len*mg(g)%npanx) ! message unit size
      call mgcollectreduce_work( g, vdat, dsolmax, kx, mg(g)%nmax, msg_len, mg(g)%npanx )

      call END_LOG(mgcollect_end)
  
   return
   end subroutine mgcollectreduce

   subroutine mgcollectreduce_work(g,vdat,dsolmax,kx,nmax,msg_len,npanx)

      integer, intent(in) :: g, kx, nmax, msg_len, npanx
      integer n, iq_a, iq_c
      integer nrow, ncol, na, nb
      integer yproc, ir, ic, is, js, je, jj, k
      integer nrm1, hoz_len
      integer(kind=4) :: ierr, ilen, lcomm
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      real, dimension(:,:), intent(inout) :: vdat
      real, dimension(:), intent(inout) :: dsolmax
      real, dimension(msg_len*npanx+1,kx) :: tdat
      real, dimension(msg_len*npanx+1,kx,nmax) :: tdat_g

      ! prep data for sending around the merge
      nrow    = mg(g)%ipan/mg(g)%merge_row  ! number of points along a row per processor
      ncol    = msg_len/nrow                ! number of points along a col per processor
      nrm1    = nrow - 1
      hoz_len = msg_len*npanx

      ! pack contiguous buffer
      tdat(1:hoz_len,1:kx) = vdat(1:hoz_len,1:kx)
      tdat(hoz_len+1,1:kx) = dsolmax(1:kx)

      ilen = (hoz_len+1)*kx
      lcomm = mg(g)%comm_merge
      call MPI_Gather( tdat, ilen, ltype, tdat_g, ilen, ltype, 0_4, lcomm, ierr )      

      ! unpack buffers (nmax is zero unless this is the host processor)
      if ( nmax > 0 ) then
         if ( hoz_len == 1 ) then
            ! usual case
            do yproc = 1,nmax
               ir = mod(yproc-1,mg(g)%merge_row) + 1   ! index for proc row
               ic = (yproc-1)/mg(g)%merge_row + 1      ! index for proc col
               iq_a = ir + (ic-1)*mg(g)%ipan
               vdat(iq_a,1:kx) = tdat_g(1,1:kx,yproc)
            end do
         else
            ! general case
            do yproc = 1,nmax
               ir = mod(yproc-1,mg(g)%merge_row) + 1   ! index for proc row
               ic = (yproc-1)/mg(g)%merge_row + 1      ! index for proc col
               is = (ir-1)*nrow + 1
               js = (ic-1)*ncol + 1
               je = ic*ncol
               do k = 1,kx
                  do n = 1,npanx
                     na = is + (n-1)*msg_len*nmax
                     nb =  1 + (n-1)*msg_len
                     do jj = js,je
                        iq_a = na + (jj-1)*mg(g)%ipan
                        iq_c = nb + (jj-js)*nrow
                        vdat(iq_a:iq_a+nrm1,k) = tdat_g(iq_c:iq_c+nrm1,k,yproc)
                     end do
                  end do
               end do
            end do
         end if
         dsolmax(1:kx) = maxval( tdat_g(hoz_len+1,1:kx,:), dim=2 )
      end if
  
   return
   end subroutine mgcollectreduce_work

   ! This routing collects data from other processors without a reduction (max or min) array
   subroutine mgcollect1(g,vdat,klim)

      integer, intent(in) :: g
      integer, intent(in), optional :: klim
      integer kx, msg_len
      real, dimension(:,:), intent(inout) :: vdat

      ! merge length
      if ( mg(g)%merge_len <= 1 ) return

      call START_LOG(mgcollect_begin)

      if (present(klim)) then
         kx = klim
      else
         kx = size(vdat,2)
      end if

      msg_len = mg(g)%ifull/(mg(g)%merge_len*mg(g)%npanx) ! message unit size
      call mgcollect_work( g, vdat, kx, mg(g)%nmax, msg_len, mg(g)%npanx )

      call END_LOG(mgcollect_end)
  
   return
   end subroutine mgcollect1

   subroutine mgcollect_work(g,vdat,kx,nmax,msg_len,npanx)

      integer, intent(in) :: g, kx, nmax, msg_len, npanx
      integer n, iq_a, iq_c
      integer nrow, ncol, na, nb
      integer yproc, ir, ic, is, js, je, jj, k
      integer nrm1, hoz_len
      integer(kind=4) :: ierr, ilen, lcomm
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      real, dimension(:,:), intent(inout) :: vdat
      real, dimension(msg_len*npanx,kx) :: tdat
      real, dimension(msg_len*npanx,kx,nmax) :: tdat_g

      ! prep data for sending around the merge
      nrow    = mg(g)%ipan/mg(g)%merge_row       ! number of points along a row per processor
      ncol    = msg_len/nrow                     ! number of points along a col per processor
      nrm1    = nrow - 1
      hoz_len = msg_len*npanx

      ! pack contiguous buffer
      tdat(1:hoz_len,1:kx) = vdat(1:hoz_len,1:kx)

      ilen = hoz_len*kx
      lcomm = mg(g)%comm_merge
      call MPI_Gather( tdat, ilen, ltype, tdat_g, ilen, ltype, 0_4, lcomm, ierr )

      ! unpack buffers (nmax is zero unless this is the host processor)
      if ( nmax>0 ) then
         if ( hoz_len == 1 ) then
            ! usual case
            do yproc = 1,nmax
               ir = mod(yproc-1,mg(g)%merge_row) + 1   ! index for proc row
               ic = (yproc-1)/mg(g)%merge_row + 1      ! index for proc col
               iq_a = ir + (ic-1)*mg(g)%ipan
               vdat(iq_a,1:kx) = tdat_g(1,1:kx,yproc)
            end do
         else
            ! general case      
            do yproc = 1,nmax
               ir = mod(yproc-1,mg(g)%merge_row) + 1   ! index for proc row
               ic = (yproc-1)/mg(g)%merge_row + 1      ! index for proc col
               is = (ir-1)*nrow + 1
               js = (ic-1)*ncol + 1
               je = ic*ncol
               do k = 1,kx
                  do n = 1,npanx
                     na = is + (n-1)*msg_len*nmax
                     nb =  1 + (n-1)*msg_len
                     do jj = js,je
                        iq_a = na + (jj-1)*mg(g)%ipan
                        iq_c = nb + (jj-js)*nrow
                        vdat(iq_a:iq_a+nrm1,k) = tdat_g(iq_c:iq_c+nrm1,k,yproc)
                     end do
                  end do
               end do
            end do
	 end if
      end if
  
   return
   end subroutine mgcollect_work

   ! This version of mgcollect also performs a max and min reduction
   subroutine mgcollectxn(g,vdat,smaxmin,klim)

      integer, intent(in) :: g
      integer, intent(in), optional :: klim
      integer kx, msg_len
      real, dimension(:,:), intent(inout) :: vdat
      real, dimension(:,:), intent(inout) :: smaxmin

      ! merge length
      if ( mg(g)%merge_len <= 1 ) return

      call START_LOG(mgcollect_begin)

      if (present(klim)) then
         kx = klim
      else
         kx = size(vdat,2)
      end if

      msg_len = mg(g)%ifull/(mg(g)%merge_len*mg(g)%npanx) ! message unit size
      call mgcollectxn_work( g, vdat, smaxmin, kx, mg(g)%nmax, msg_len, mg(g)%npanx )

      call END_LOG(mgcollect_end)
  
   return
   end subroutine mgcollectxn

   subroutine mgcollectxn_work(g,vdat,smaxmin,kx,nmax,msg_len,npanx)

      integer, intent(in) :: g, kx, nmax, msg_len, npanx
      integer :: n, iq_a, iq_c
      integer :: nrow, ncol, na, nb
      integer :: yproc, ir, ic, is, js, je, jj, k 
      integer :: nrm1, hoz_len
      integer(kind=4) :: ierr, ilen, lcomm
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      real, dimension(:,:), intent(inout) :: vdat
      real, dimension(:,:), intent(inout) :: smaxmin
      real, dimension(msg_len*npanx+2,kx) :: tdat
      real, dimension(msg_len*npanx+2,kx,nmax) :: tdat_g

      ! prep data for sending around the merge
      nrow    = mg(g)%ipan/mg(g)%merge_row  ! number of points along a row per processor
      ncol    = msg_len/nrow                ! number of points along a col per processor
      nrm1    = nrow - 1
      hoz_len = msg_len*npanx

      ! pack contiguous buffer
      tdat(1:hoz_len,1:kx) = vdat(1:hoz_len,1:kx) 
      tdat(hoz_len+1,1:kx) = smaxmin(1:kx,1)
      tdat(hoz_len+2,1:kx) = smaxmin(1:kx,2)

      ilen = (hoz_len+2)*kx
      lcomm = mg(g)%comm_merge
      call MPI_Gather( tdat, ilen, ltype, tdat_g, ilen, ltype, 0_4, lcomm, ierr )

      ! unpack buffers (nmax is zero unless this is the host processor)
      if ( nmax > 0 ) then
         if ( hoz_len == 1 ) then
            ! usual case
            do yproc = 1,nmax
               ir = mod(yproc-1,mg(g)%merge_row) + 1   ! index for proc row
               ic = (yproc-1)/mg(g)%merge_row + 1      ! index for proc col
               iq_a = ir + (ic-1)*mg(g)%ipan
               vdat(iq_a,1:kx) = tdat_g(1,1:kx,yproc)
            end do
         else
            ! general case          
            do yproc = 1,nmax
               ir = mod(yproc-1,mg(g)%merge_row) + 1   ! index for proc row
               ic = (yproc-1)/mg(g)%merge_row + 1      ! index for proc col
               is = (ir-1)*nrow + 1
               js = (ic-1)*ncol + 1
               je = ic*ncol
               do k = 1,kx
                  do n = 1,npanx
                     na = is + (n-1)*msg_len*nmax
                     nb =  1 + (n-1)*msg_len
                     do jj = js,je
                        iq_a = na + (jj-1)*mg(g)%ipan
                        iq_c = nb + (jj-js)*nrow
                        vdat(iq_a:iq_a+nrm1,k) = tdat_g(iq_c:iq_c+nrm1,k,yproc)
                     end do
                  end do
               end do
            end do
         end if
         smaxmin(1:kx,1) = maxval( tdat_g(hoz_len+1,1:kx,:), dim=2 )
         smaxmin(1:kx,2) = minval( tdat_g(hoz_len+2,1:kx,:), dim=2 )
      end if
  
   return
   end subroutine mgcollectxn_work

   ! This subroutine merges datasets when upscaling with the multi-grid solver
   ! This version also updates the halo
   subroutine mgbcast(g,vdat,dsolmax,klim)
   
      integer, intent(in) :: g
      integer, intent(in), optional :: klim
      integer :: kx, out_len
      real, dimension(:,:), intent(inout) :: vdat
      real, dimension(:), intent(inout) :: dsolmax

      ! merge length
      if ( mg(g)%merge_len <= 1 ) return
   
      call START_LOG(mgbcast_begin)
   
      if ( present(klim) ) then
         kx = klim
      else
         kx = size(vdat,2)
      end if
   
      out_len = mg(g)%ifull + mg(g)%iextra
      call mgbcast_work( g, vdat, dsolmax, kx, out_len )
   
      call END_LOG(mgbcast_end)
   
   return
   end subroutine mgbcast
   
   subroutine mgbcast_work(g,vdat,dsolmax,kx,out_len)
   
      integer, intent(in) :: g, kx, out_len
      integer(kind=4) :: ierr, ilen, lcomm
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      real, dimension(:,:), intent(inout) :: vdat
      real, dimension(:), intent(inout) :: dsolmax
      real, dimension(out_len+1,kx) :: tdat

      ! pack contiguous buffer
      tdat(1:out_len,1:kx) = vdat(1:out_len,1:kx)
      tdat(out_len+1,1:kx) = dsolmax(1:kx)

      ilen = (out_len+1)*kx
      lcomm = mg(g)%comm_merge
      call MPI_Bcast( tdat, ilen, ltype, 0_4, lcomm, ierr )

      ! extract data from Bcast
      vdat(1:out_len,1:kx) = tdat(1:out_len,1:kx)
      dsolmax(1:kx) = tdat(out_len+1,1:kx)

   return
   end subroutine mgbcast_work

   subroutine mgbcasta(g,vdat)
   
      integer, intent(in) :: g
      integer :: kx, out_len
      real, dimension(:,:), intent(inout) :: vdat

      ! merge length
      if ( mg(g)%merge_len <= 1 ) return
   
      call START_LOG(mgbcast_begin)
   
      kx = size(vdat,2)
      out_len = mg(g)%ifull + mg(g)%iextra
      call mgbcasta_work( g, vdat, kx, out_len )
   
      call END_LOG(mgbcast_end)
   
   return
   end subroutine mgbcasta
   
   subroutine mgbcasta_work(g,vdat,kx,out_len)
   
      integer, intent(in) :: g, kx, out_len
      integer(kind=4) :: ierr, ilen, lcomm
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      real, dimension(:,:), intent(inout) :: vdat
      real, dimension(out_len,kx) :: tdat

      ! pack contiguous buffer
      tdat(1:out_len,1:kx) = vdat(1:out_len,1:kx)
      
      ilen = out_len*kx
      lcomm = mg(g)%comm_merge
      call MPI_Bcast( tdat, ilen, ltype, 0_4, lcomm, ierr )      

      ! extract data from Bcast
      vdat(1:out_len,1:kx) = tdat(1:out_len,1:kx)
   
   return
   end subroutine mgbcasta_work

   subroutine mgbcastxn(g,vdat,smaxmin,klim)
   
      integer, intent(in) :: g
      integer, intent(in), optional :: klim
      integer :: kx, out_len
      real, dimension(:,:), intent(inout) :: vdat
      real, dimension(:,:), intent(inout) :: smaxmin

      ! merge length
      if ( mg(g)%merge_len <= 1 ) return
   
      call START_LOG(mgbcast_begin)
   
      if (present(klim)) then
         kx = klim
      else
         kx = size(vdat,2)
      end if
   
      out_len = mg(g)%ifull + mg(g)%iextra
      call mgbcastxn_work( g, vdat, smaxmin, kx, out_len )
   
      call END_LOG(mgbcast_end)
   
   return
   end subroutine mgbcastxn
   
   subroutine mgbcastxn_work(g,vdat,smaxmin,kx,out_len)
   
      integer, intent(in) :: g, kx, out_len
      integer(kind=4) :: ierr, ilen, lcomm
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      real, dimension(:,:), intent(inout) :: vdat
      real, dimension(:,:), intent(inout) :: smaxmin
      real, dimension(out_len+2,1:kx) :: tdat

      ! pack contiguous buffer
      tdat(1:out_len,1:kx) = vdat(1:out_len,1:kx)
      tdat(out_len+1,1:kx) = smaxmin(1:kx,1)
      tdat(out_len+2,1:kx) = smaxmin(1:kx,2)
      
      ilen = (out_len+2)*kx
      lcomm = mg(g)%comm_merge
      call MPI_Bcast( tdat, ilen, ltype, 0_4, lcomm, ierr )      

      ! extract data from Bcast
      vdat(1:out_len,1:kx) = tdat(1:out_len,1:kx)
      smaxmin(1:kx,1) = tdat(out_len+1,1:kx)
      smaxmin(1:kx,2) = tdat(out_len+2,1:kx)
   
   return
   end subroutine mgbcastxn_work

   ! Set up the indices required for the multigrid scheme.
   subroutine mg_index(g,mil_g,mipan,mjpan)

      use indices_m
   
      integer, intent(in) :: g, mil_g, mipan, mjpan
      integer, dimension(:), allocatable :: mg_colourmask
      integer, dimension(2*(mipan+mjpan+2)*(npanels+1)) :: dum
      integer, dimension(2,0:nproc-1) :: sdum, rdum
      integer, dimension(3) :: mg_ifullcol
      integer mioff, mjoff
      integer i, j, n, iq, iqq, iqg, iql, iqb, iqtmp, mfull_g
      integer iloc, jloc, nloc
      integer iext, iproc, xlen, jx, nc, xlev, rproc, sproc
      integer ntest, ncount
      integer(kind=4) :: itag=22, lproc, ierr, llen
      ! 13 is the maximum number of possibe neigbours (i.e., uniform decomposition).
      integer(kind=4), dimension(26) :: dreq
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
         call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
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
            mg(g)%inw(iq) = jnw_g(iq,mil_g)
            mg(g)%ise(iq) = jse_g(iq,mil_g)
            mg(g)%isw(iq) = jsw_g(iq,mil_g)
         end do
         mg(g)%ixlen = 0
         mg(g)%iextra = 0
         mg(g)%neighnum = 0
         allocate ( mg(g)%neighlist(mg(g)%neighnum) )
      else
         mg(g)%iextra = 2*(mipan+mjpan+2)*npan ! first guess

         ! This only occurs with grids prior to globgath.  So npan and noff are still valid.
         do n = 1,npan
            do j = 1,mjpan
               do i = 1,mipan
                  iq = indx(i,j,n-1,mipan,mjpan) ! Local
                  iqg = indx(i+mioff,j+mjoff,n-noff,mil_g,mil_g) ! Global

                  iqq = jn_g(iqg,mil_g)       ! Global neighbour index
                  rproc = mg_qproc(iqq,mil_g,g) ! Processor that has this point
                  if ( rproc == myid ) then ! Just copy the value
                     ! Convert global iqq to local value
                     call indv_mpix(iqq,iloc,jloc,nloc,mil_g,mioff,mjoff,noff)
                     mg(g)%in(iq) = indx(iloc,jloc,nloc-1,mipan,mjpan)
                  end if

                  iqq = js_g(iqg,mil_g)       ! Global neighbour index
                  rproc = mg_qproc(iqq,mil_g,g) ! Processor that has this point
                  if ( rproc == myid ) then ! Just copy the value
                     ! Convert global iqq to local value
                     call indv_mpix(iqq,iloc,jloc,nloc,mil_g,mioff,mjoff,noff)
                     mg(g)%is(iq) = indx(iloc,jloc,nloc-1,mipan,mjpan)
                  end if

                  iqq = je_g(iqg,mil_g)       ! Global neighbour index
                  rproc = mg_qproc(iqq,mil_g,g) ! Processor that has this point
                  if ( rproc == myid ) then ! Just copy the value
                     ! Convert global iqq to local value
                     call indv_mpix(iqq,iloc,jloc,nloc,mil_g,mioff,mjoff,noff)
                     mg(g)%ie(iq) = indx(iloc,jloc,nloc-1,mipan,mjpan)
                  end if

                  iqq = jw_g(iqg,mil_g)    ! Global neighbour index
                  rproc = mg_qproc(iqq,mil_g,g) ! Processor that has this point
                  if ( rproc == myid ) then ! Just copy the value
                     ! Convert global iqq to local value
                     call indv_mpix(iqq,iloc,jloc,nloc,mil_g,mioff,mjoff,noff)
                     mg(g)%iw(iq) = indx(iloc,jloc,nloc-1,mipan,mjpan)
                  end if

                  ! Note that the model only needs a limited set of the diagonal
                  ! index arrays
                  iqq = jne_g(iqg,mil_g)    ! Global neighbour index
                  rproc = mg_qproc(iqq,mil_g,g) ! Processor that has this point
                  if ( rproc == myid ) then ! Just copy the value
                     ! Convert global iqq to local value
                     call indv_mpix(iqq,iloc,jloc,nloc,mil_g,mioff,mjoff,noff)
                     mg(g)%ine(iq) = indx(iloc,jloc,nloc-1,mipan,mjpan)
                  end if

                  iqq = jse_g(iqg,mil_g)    ! Global neighbour index
                  rproc = mg_qproc(iqq,mil_g,g) ! Processor that has this point
                  if ( rproc == myid ) then ! Just copy the value
                     ! Convert global iqq to local value
                     call indv_mpix(iqq,iloc,jloc,nloc,mil_g,mioff,mjoff,noff)
                     mg(g)%ise(iq) = indx(iloc,jloc,nloc-1,mipan,mjpan)
                  end if

                  iqq = jnw_g(iqg,mil_g)    ! Global neighbour index
                  rproc = mg_qproc(iqq,mil_g,g) ! Processor that has this point
                  if ( rproc == myid ) then ! Just copy the value
                     ! Convert global iqq to local value
                     call indv_mpix(iqq,iloc,jloc,nloc,mil_g,mioff,mjoff,noff)
                     mg(g)%inw(iq) = indx(iloc,jloc,nloc-1,mipan,mjpan)
                  end if

                  iqq = jsw_g(iqg,mil_g)    ! Global neighbour index
                  rproc = mg_qproc(iqq,mil_g,g) ! Processor that has this point
                  if ( rproc == myid ) then ! Just copy the value
                     ! Convert global iqq to local value
                     call indv_mpix(iqq,iloc,jloc,nloc,mil_g,mioff,mjoff,noff)
                     mg(g)%isw(iq) = indx(iloc,jloc,nloc-1,mipan,mjpan)
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
               iq = indx(i+mioff,j+mjoff,n-noff,mil_g,mil_g)
               iqq = jn_g(iq,mil_g)
               ! Which processor has this point
               rproc = mg_qproc(iqq,mil_g,g)
               if ( rproc == myid ) cycle ! Don't add points already on this proc.
               iql = indx(i,j,n-1,mipan,mjpan)  !  Local index
               ! Add this point to request list
               call mgcheck_bnds_alloc(g, rproc, iext)
               iqtmp = -1
               do iqb = 1,mg_bnds(rproc,g)%rlen
                  if ( mg_bnds(rproc,g)%request_list(iqb) == iqq ) then
                     iqtmp = mg_bnds(rproc,g)%unpack_list(iqb)
                     exit
                  end if
               end do
               if ( iqtmp > 0 ) then
                  mg(g)%in(iql) = mg(g)%ifull+iqtmp
               else
                  mg_bnds(rproc,g)%rlen = mg_bnds(rproc,g)%rlen + 1
                  mg_bnds(rproc,g)%request_list(mg_bnds(rproc,g)%rlen) = iqq
                  ! Increment extended region index
                  iext = iext + 1
                  mg_bnds(rproc,g)%unpack_list(mg_bnds(rproc,g)%rlen) = iext
                  mg(g)%in(iql) = mg(g)%ifull+iext
               end if
            end do

            !     E edge
            i = mipan
            do j=1,mjpan
               iq = indx(i+mioff,j+mjoff,n-noff,mil_g,mil_g)
               iqq = je_g(iq,mil_g)
               ! Which processor has this point
               rproc = mg_qproc(iqq,mil_g,g)
               if ( rproc == myid ) cycle ! Don't add points already on this proc.
               iql = indx(i,j,n-1,mipan,mjpan)  !  Local index
               ! Add this point to request list
               call mgcheck_bnds_alloc(g, rproc, iext)
               iqtmp = -1
               do iqb = 1,mg_bnds(rproc,g)%rlen
                  if ( mg_bnds(rproc,g)%request_list(iqb) == iqq ) then
                     iqtmp = mg_bnds(rproc,g)%unpack_list(iqb)
                     exit
                  end if
               end do
               if ( iqtmp > 0 ) then
                  mg(g)%ie(iql) = mg(g)%ifull+iqtmp
               else
                  mg_bnds(rproc,g)%rlen = mg_bnds(rproc,g)%rlen + 1
                  mg_bnds(rproc,g)%request_list(mg_bnds(rproc,g)%rlen) = iqq
                  ! Increment extended region index
                  iext = iext + 1
                  mg_bnds(rproc,g)%unpack_list(mg_bnds(rproc,g)%rlen) = iext
                  mg(g)%ie(iql) = mg(g)%ifull+iext
               end if
            end do

            !     W edge
            i = 1
            do j = 1,mjpan
               iq = indx(i+mioff,j+mjoff,n-noff,mil_g,mil_g)
               iqq = jw_g(iq,mil_g)
               ! Which processor has this point
               rproc = mg_qproc(iqq,mil_g,g)
               if ( rproc == myid ) cycle ! Don't add points already on this proc.
               iql = indx(i,j,n-1,mipan,mjpan)  !  Local index
               ! Add this point to request list
               call mgcheck_bnds_alloc(g, rproc, iext)
               iqtmp = -1
               do iqb = 1,mg_bnds(rproc,g)%rlen
                  if ( mg_bnds(rproc,g)%request_list(iqb) == iqq ) then
                     iqtmp = mg_bnds(rproc,g)%unpack_list(iqb)
                     exit
                  end if
               end do
               if ( iqtmp > 0 ) then
                  mg(g)%iw(iql) = mg(g)%ifull+iqtmp
               else
                  mg_bnds(rproc,g)%rlen = mg_bnds(rproc,g)%rlen + 1
                  mg_bnds(rproc,g)%request_list(mg_bnds(rproc,g)%rlen) = iqq
                  ! Increment extended region index
                  iext = iext + 1
                  mg_bnds(rproc,g)%unpack_list(mg_bnds(rproc,g)%rlen) = iext
                  mg(g)%iw(iql) = mg(g)%ifull+iext
               end if
            end do

            !     S edge
            j = 1
            do i = 1,mipan
               iq = indx(i+mioff,j+mjoff,n-noff,mil_g,mil_g)
               iqq = js_g(iq,mil_g)
               ! Which processor has this point
               rproc = mg_qproc(iqq,mil_g,g)
               if ( rproc == myid ) cycle ! Don't add points already on this proc.
               iql = indx(i,j,n-1,mipan,mjpan)  !  Local index
               ! Add this point to request list
               call mgcheck_bnds_alloc(g, rproc, iext)
               iqtmp = -1
               do iqb = 1,mg_bnds(rproc,g)%rlen
                  if ( mg_bnds(rproc,g)%request_list(iqb) == iqq ) then
                     iqtmp = mg_bnds(rproc,g)%unpack_list(iqb)
                     exit
                  end if
               end do
               if ( iqtmp > 0 ) then
                  mg(g)%is(iql) = mg(g)%ifull+iqtmp
               else
                  mg_bnds(rproc,g)%rlen = mg_bnds(rproc,g)%rlen + 1
                  mg_bnds(rproc,g)%request_list(mg_bnds(rproc,g)%rlen) = iqq
                  ! Increment extended region index
                  iext = iext + 1
                  mg_bnds(rproc,g)%unpack_list(mg_bnds(rproc,g)%rlen) = iext
                  mg(g)%is(iql) = mg(g)%ifull+iext
               end if
            end do
         end do ! n=1,npan

         mg(g)%ixlen = iext
         mg_bnds(:,g)%rlenx = mg_bnds(:,g)%rlen  ! so that they're appended.
      
         ! Now handle the special corner values that need to be remapped
         do n = 1,npan
            ! NE
            iq = indx(mipan,mjpan,n-1,mipan,mjpan)
            iqg = indx(mipan+mioff,mjpan+mjoff,n-noff,mil_g,mil_g)
            iqq = jne_g(iqg,mil_g)
            ! Which processor has this point
            rproc = mg_qproc(iqq,mil_g,g)
            if ( rproc /= myid ) then ! Add to list
               call mgcheck_bnds_alloc(g, rproc, iext)
               iqtmp=-1
               do iqb=1,mg_bnds(rproc,g)%rlenx
                  if ( mg_bnds(rproc,g)%request_list(iqb) == iqq ) then
                     iqtmp = mg_bnds(rproc,g)%unpack_list(iqb)
                     exit
                  end if
               end do
               if ( iqtmp > 0 ) then
                  mg(g)%ine(iq) = mg(g)%ifull+iqtmp
               else
                  mg_bnds(rproc,g)%rlenx = mg_bnds(rproc,g)%rlenx + 1
                  mg_bnds(rproc,g)%request_list(mg_bnds(rproc,g)%rlenx) = iqq
                  ! Increment extended region index
                  iext = iext + 1
                  mg_bnds(rproc,g)%unpack_list(mg_bnds(rproc,g)%rlenx) = iext
                  mg(g)%ine(iq) = mg(g)%ifull+iext
               end if
            end if

            ! SE
            iq = indx(mipan,1,n-1,mipan,mjpan)
            iqg = indx(mipan+mioff,1+mjoff,n-noff,mil_g,mil_g)
            iqq = jse_g(iqg,mil_g)
            ! Which processor has this point
            rproc = mg_qproc(iqq,mil_g,g)
            if ( rproc /= myid ) then ! Add to list
               call mgcheck_bnds_alloc(g, rproc, iext)
               iqtmp=-1
               do iqb=1,mg_bnds(rproc,g)%rlenx
                  if ( mg_bnds(rproc,g)%request_list(iqb) == iqq ) then
                     iqtmp = mg_bnds(rproc,g)%unpack_list(iqb)
                     exit
                  end if
               end do
               if ( iqtmp > 0 ) then
                  mg(g)%ise(iq) = mg(g)%ifull+iqtmp
               else
                  mg_bnds(rproc,g)%rlenx = mg_bnds(rproc,g)%rlenx + 1
                  mg_bnds(rproc,g)%request_list(mg_bnds(rproc,g)%rlenx) = iqq
                  ! Increment extended region index
                  iext = iext + 1
                  mg_bnds(rproc,g)%unpack_list(mg_bnds(rproc,g)%rlenx) = iext
                  mg(g)%ise(iq) = mg(g)%ifull+iext
               end if
            end if

            ! WN
            iq = indx(1,mjpan,n-1,mipan,mjpan)
            iqg = indx(1+mioff,mjpan+mjoff,n-noff,mil_g,mil_g)
            iqq = jnw_g(iqg,mil_g)
            ! Which processor has this point
            rproc = mg_qproc(iqq,mil_g,g)
            if ( rproc /= myid ) then ! Add to list
               call mgcheck_bnds_alloc(g, rproc, iext)
               iqtmp=-1
               do iqb=1,mg_bnds(rproc,g)%rlenx
                  if ( mg_bnds(rproc,g)%request_list(iqb) == iqq ) then
                     iqtmp = mg_bnds(rproc,g)%unpack_list(iqb)
                     exit
                  end if
               end do
               if ( iqtmp > 0 ) then
                  mg(g)%inw(iq) = mg(g)%ifull+iqtmp
               else
                  mg_bnds(rproc,g)%rlenx = mg_bnds(rproc,g)%rlenx + 1
                  mg_bnds(rproc,g)%request_list(mg_bnds(rproc,g)%rlenx) = iqq
                  ! Increment extended region index
                  iext = iext + 1
                  mg_bnds(rproc,g)%unpack_list(mg_bnds(rproc,g)%rlenx) = iext
                  mg(g)%inw(iq) = mg(g)%ifull+iext
               end if
            end if

            ! SW
            iq = indx(1,1,n-1,mipan,mjpan)
            iqg = indx(1+mioff,1+mjoff,n-noff,mil_g,mil_g)
            iqq = jsw_g(iqg,mil_g)
            ! Which processor has this point
            rproc = mg_qproc(iqq,mil_g,g)
            if ( rproc /= myid ) then ! Add to list
               call mgcheck_bnds_alloc(g, rproc, iext)
               iqtmp=-1
               do iqb=1,mg_bnds(rproc,g)%rlenx
                  if ( mg_bnds(rproc,g)%request_list(iqb) == iqq ) then
                     iqtmp = mg_bnds(rproc,g)%unpack_list(iqb)
                     exit
                  end if
               end do
               if ( iqtmp > 0 ) then
                  mg(g)%isw(iq) = mg(g)%ifull+iqtmp
               else
                  mg_bnds(rproc,g)%rlenx = mg_bnds(rproc,g)%rlenx + 1
                  mg_bnds(rproc,g)%request_list(mg_bnds(rproc,g)%rlenx) = iqq
                  ! Increment extended region index
                  iext = iext + 1
                  mg_bnds(rproc,g)%unpack_list(mg_bnds(rproc,g)%rlenx) = iext
                  mg(g)%isw(iq) = mg(g)%ifull+iext
               end if
            end if

         end do
         mg(g)%iextra = iext

         ! Set up the diagonal index arrays.
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
               end do
            end do
         end do


         mg(g)%neighnum = count( mg_bnds(:,g)%rlenx > 0 )
         if ( mg(g)%neighnum > 13 ) then
            write(6,*) "ERROR: More than 13 MG neighbours at level ",g
            call MPI_Abort( MPI_COMM_WORLD, -1_4, ierr )
         end if

         ! Now, for each processor send the length of points I want.
         nreq = 0
         llen = 2
         do iproc = 1,nproc-1
            rproc = modulo( myid+iproc, nproc )
            if ( mg_bnds(rproc,g)%rlenx>0 ) then
               nreq = nreq + 1
               lproc = rproc
               call MPI_IRecv( rdum(:,rproc), llen, ltype, lproc, itag, MPI_COMM_WORLD, dreq(nreq), ierr )
            end if
         end do
         do iproc = 1,nproc-1
            sproc = modulo( myid-iproc, nproc )  ! Send to
            if ( mg_bnds(sproc,g)%rlenx>0 ) then
               nreq = nreq + 1
               sdum(1,sproc) = mg_bnds(sproc,g)%rlenx
               sdum(2,sproc) = mg_bnds(sproc,g)%rlen
               lproc = sproc
               call MPI_ISend( sdum(:,sproc), llen, ltype, lproc, itag, MPI_COMM_WORLD, dreq(nreq), ierr )
            end if
         end do
         call MPI_Waitall(nreq,dreq,MPI_STATUSES_IGNORE,ierr)
         nreq = 0
         rreq = 0

         do iproc = 1,nproc-1
            rproc = modulo( myid+iproc, nproc )
            if ( mg_bnds(rproc,g)%rlenx>0 ) then
               mg_bnds(rproc,g)%slenx = rdum(1,rproc)
               mg_bnds(rproc,g)%slen  = rdum(2,rproc)
            end if
         end do

         ! increase size of request list if needed
         ntest = mg(g)%neighnum
         if ( 2*ntest > size(ireq) ) then
            deallocate( ireq, rlist )
            allocate( ireq(2*ntest) )
            allocate( rlist(ntest) )
         end if
  
         ! set up neighbour lists
         allocate ( mg(g)%neighlist(mg(g)%neighnum) )
         ncount = 0
         do iproc = 1,nproc-1
            rproc = modulo( myid+iproc, nproc )
            if ( mg_bnds(rproc,g)%rlenx>0 ) then
               ncount = ncount + 1
               mg(g)%neighlist(ncount) = rproc
            end if
         end do
      
         if ( ncount/=mg(g)%neighnum ) then
            write(6,*) "ERROR: Multi-grid neighnum mismatch"
            write(6,*) "neighnum, ncount ",mg(g)%neighnum, ncount
            call MPI_Abort( MPI_COMM_WORLD, -1_4, ierr )
         end if
  
         ! Now start sending messages  
         nreq = 0
         do iproc = 1,mg(g)%neighnum
            lproc = mg(g)%neighlist(iproc)  ! Recv from
            allocate( mg_bnds(lproc,g)%send_list(mg_bnds(lproc,g)%slenx) )
            nreq = nreq + 1
            ! Use the maximum size in the recv call.
            llen = mg_bnds(lproc,g)%slenx
            call MPI_IRecv( mg_bnds(lproc,g)%send_list(1), llen, ltype, lproc, &
                            itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end do
         do iproc = mg(g)%neighnum,1,-1
            lproc = mg(g)%neighlist(iproc)  ! Send to
            ! Send list of requests
            nreq = nreq + 1
            llen = mg_bnds(lproc,g)%rlenx
            call MPI_ISend( mg_bnds(lproc,g)%request_list(1), llen, ltype, lproc, &
                            itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end do      
         call MPI_Waitall(nreq,ireq,MPI_STATUSES_IGNORE,ierr)
         nreq = 0
         rreq = 0

         ! At the moment send_lists use global indices. Convert these to local.
         do iproc = mg(g)%neighnum,1,-1
            sproc = mg(g)%neighlist(iproc)  ! Send to
            do iq = 1,mg_bnds(sproc,g)%slenx
               ! send_list(iq) is global point index, i, j, n are local
               iqq = mg_bnds(sproc,g)%send_list(iq)
               call indv_mpix(iqq,i,j,n,mil_g,mioff,mjoff,noff)
               mg_bnds(sproc,g)%send_list(iq) = indx(i,j,n-1,mipan,mjpan)
            end do
         end do
         do iq = 1,mg_bnds(myid,g)%rlenx
            iqq = mg_bnds(myid,g)%request_list(iq)
            call indv_mpix(iqq,i,j,n,mil_g,mioff,mjoff,noff)
            mg_bnds(myid,g)%request_list(iq) = indx(i,j,n-1,mipan,mjpan)
         end do

         ! reduce array size where possible
         do iproc = 0,nproc-1
            xlen = mg_bnds(iproc,g)%rlenx
            if ( mg_bnds(iproc,g)%len > xlen ) then
               dum(1:xlen) = mg_bnds(iproc,g)%request_list(1:xlen)
               deallocate( mg_bnds(iproc,g)%request_list )
               if ( xlen > 0 ) then
                 allocate( mg_bnds(iproc,g)%request_list(xlen) )
                 mg_bnds(iproc,g)%request_list(1:xlen) = dum(1:xlen)
               end if
               dum(1:xlen) = mg_bnds(iproc,g)%unpack_list(1:xlen)
               deallocate( mg_bnds(iproc,g)%unpack_list )
               if ( xlen > 0 ) then
                 allocate( mg_bnds(iproc,g)%unpack_list(xlen) )
                 mg_bnds(iproc,g)%unpack_list(1:xlen) = dum(1:xlen)
               end if
               mg_bnds(iproc,g)%len = xlen
            end if

            ! set up buffers
            xlev = max(kl,ol)
            xlen = xlev*mg_bnds(iproc,g)%rlenx
            if ( bnds(iproc)%rbuflen < xlen ) then
               if ( bnds(iproc)%rbuflen > 0 ) then
                  deallocate( bnds(iproc)%rbuf )
                  deallocate( bnds(iproc)%r8buf )
               end if
               allocate( bnds(iproc)%rbuf(xlen) )
               allocate( bnds(iproc)%r8buf(xlen) )
               bnds(iproc)%rbuflen = xlen
            end if
            xlen = xlev*mg_bnds(iproc,g)%slenx
            if ( bnds(iproc)%sbuflen < xlen ) then
               if ( bnds(iproc)%sbuflen > 0 ) then
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
  
         mg_ifullmaxcol = count( mg_colourmask == 1 )
         if ( mg_ifullmaxcol /= count( mg_colourmask == 2 ) .or. mg_ifullmaxcol /= count( mg_colourmask == 3 ) ) then
           write(6,*) "ERROR: Unbalanced MG colours"
           call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
         end if
         allocate( col_iq(mg_ifullmaxcol,3),  col_iqn(mg_ifullmaxcol,3), col_iqe(mg_ifullmaxcol,3) )
         allocate( col_iqs(mg_ifullmaxcol,3), col_iqw(mg_ifullmaxcol,3) )
  
         mg_ifullcol = 0
         col_iq = 0
         col_iqn = 0
         col_iqe = 0
         col_iqs = 0
         col_iqw = 0
         do iq = 1,mg(g)%ifull
            nc = mg_colourmask(iq)
            mg_ifullcol(nc) = mg_ifullcol(nc) + 1
            iqq = mg_ifullcol(nc)
            col_iq(iqq,nc) = iq
            col_iqn(iqq,nc) = mg(g)%in(iq)
            col_iqe(iqq,nc) = mg(g)%ie(iq)
            col_iqs(iqq,nc) = mg(g)%is(iq)
            col_iqw(iqq,nc) = mg(g)%iw(iq)
         end do
         
         deallocate( mg_colourmask )

      end if

   return
   end subroutine mg_index

   subroutine mgcheck_bnds_alloc(g,iproc,iext)

      integer, intent(in) :: iproc
      integer, intent(in) :: g, iext
      integer(kind=4) :: ierr

      if ( mg_bnds(iproc,g)%len <= 0 ) then
         allocate( mg_bnds(iproc,g)%request_list(mg(g)%iextra) )
         allocate( mg_bnds(iproc,g)%unpack_list(mg(g)%iextra) )
         mg_bnds(iproc,g)%len = mg(g)%iextra
      else
         if ( iext>mg(g)%iextra ) then
            write(6,*) "ERROR: MG grid undersized in mgcheck_bnds_alloc"
            write(6,*) "iext,iextra,g,iproc,myid ",iext,mg(g)%iextra,g,iproc,myid
            call MPI_Abort(MPI_COMM_WORLD,-1_4,ierr)
         end if
      end if

   return
   end subroutine mgcheck_bnds_alloc

   subroutine indv_mpix(iq, i, j, n, mil_g, mioff, mjoff, mnoff)

      integer, intent(in) :: iq, mil_g, mnoff
      integer, intent(in) :: mioff, mjoff
      integer, intent(out) :: i
      integer, intent(out) :: j
      integer, intent(out) :: n

      ! Calculate local i, j, n from global iq

      ! Global i, j, n
      n = (iq - 1)/(mil_g*mil_g)
      j = 1 + (iq - n*mil_g*mil_g - 1)/mil_g
      i = iq - (j - 1)*mil_g - n*mil_g*mil_g

      ! Reduced to values on my processor
      n = n + mnoff  
      j = j - mjoff
      i = i - mioff

   return
   end subroutine indv_mpix

   function mg_fproc_1(g,i,j,n) result(mg_fpout)
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
     
   return
   end function mg_fproc_1
   
   function mg_fproc(g,i,j,n) result(mg_fpout)
     ! locates processor that owns a global grid point
     integer, intent(in) :: i, j, n, g
     integer mg_fpout
     integer fp_l
     
     fp_l = mg_fproc_1(g,i,j,n)
     mg_fpout = mg(g)%procmap(fp_l)
     
   return
   end function mg_fproc   
   
   function mg_qproc(iqg,mil_g,g) result(mg_qpout)
      ! locates processor that owns a global grid point
      integer, intent(in) :: iqg, mil_g, g
      integer :: mg_qpout
      integer :: i, j, n

      n = (iqg - 1) / (mil_g*mil_g)
      j = 1 + (iqg - n*mil_g*mil_g - 1)/mil_g
      i = iqg - (j - 1)*mil_g - n*mil_g*mil_g

      mg_qpout = mg_fproc(g,i,j,n)
   
   return
   end function mg_qproc
   
   function indx(i,j,n,il,jl) result(iq)
      ! more general version of ind function

      integer, intent(in) :: i, j, n, il, jl
      integer iq

      iq = i + (j-1)*il + n*il*jl

   return
   end function indx
   
   function ind(i,j,n) result(iq)

      integer, intent(in) :: i, j, n
      integer iq

      iq = i + (j-1)*ipan + (n-1)*ipan*jpan

   return
   end function ind
   
   function findcolour(iqg) result(icol)
   
      integer, intent(in) :: iqg
      integer icol
      integer ig, jg, ng, tg, jx

      ! calculate global i,j,n
      tg = iqg - 1
      ng = tg/(il_g*il_g)
      tg = tg - ng*il_g*il_g
      jg = tg/il_g
      tg = tg - jg*il_g
      ig = tg
      ig = ig + 1
      jg = jg + 1
   
#ifdef uniform_decomp
      ! three colour mask
      jx = mod( ig + jg + ng*il_g, 2 )
      select case( ng + jx*(npanels+1) )
         case( 0, 1, 3, 4 )
            icol = 1
         case( 2, 5, 6, 9 )
            icol = 2
         case( 7, 8, 10, 11 )
            icol = 3
      end select
#else
      ! two colour mask
      jx = mod( ig + jg + ng*il_g, 2 )
      if (jx == 0 ) then
         icol = 1
      else
         icol = 2
      end if
#endif
   
   return
   end function findcolour

   subroutine ccmpi_filebounds_setup(procarray,comm,ik)

      integer, intent(in) :: comm, ik
      integer, dimension(-1:ik+2,-1:ik+2,0:npanels,1:2), intent(in) :: procarray
      integer, dimension(:,:), allocatable :: dummy
      integer :: ipf, n, i, j, iq, ncount, ca, cb, no, ip
      integer :: filemaxbuflen, xlen, xlev
      integer :: iproc, jproc, iloc, jloc, nloc, floc
      integer(kind=4) :: lproc, lcomm, llen, ierr
      integer(kind=4) :: itag=42
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_INTEGER8
#else
      integer(kind=4), parameter :: ltype = MPI_INTEGER
#endif
      logical, save :: fileallocate = .false.
     
      if ( myid >= fnresid ) return

      lcomm = comm
      filemaxbuflen = 2*(pil+pjl)*pnpan*fncount
      
      ! allocate memory to filebnds
      if ( fileallocate ) then
         do iproc = 0,size(filebnds)-1
            if ( filebnds(iproc)%slen > 0 ) then
               deallocate(filebnds(iproc)%send_list)
            end if
            if ( filebnds(iproc)%rlen > 0 ) then
               deallocate(filebnds(iproc)%request_list)
               deallocate(filebnds(iproc)%unpack_list)
            end if
         end do
         deallocate(filebnds)
      end if
      allocate(filebnds(0:fnresid-1))
      fileallocate = .true.

      ! calculate recv message length
      filebnds(:)%len = 0
      filebnds(:)%rlen = 0
      filebnds(:)%slen = 0
      do ipf = 0,fncount-1
         ip = ipf*fnresid + myid
         do n = 1,pnpan
            no = n - pnoff(ip)
            ca = pioff(ip,no)
            cb = pjoff(ip,no)
            do i = 1,pil
               iproc = procarray(i+ca,cb,no,1)
               floc = procarray(i+ca,cb,no,2)
               filebnds(iproc)%rlen = filebnds(iproc)%rlen + 1
               call check_filebnds_alloc(iproc,filemaxbuflen)
               ! store global index
               filebnds(iproc)%request_list(filebnds(iproc)%rlen,1) = i + ca
               filebnds(iproc)%request_list(filebnds(iproc)%rlen,2) = cb
               filebnds(iproc)%request_list(filebnds(iproc)%rlen,3) = no
               filebnds(iproc)%request_list(filebnds(iproc)%rlen,4) = floc
               ! store local index
               filebnds(iproc)%unpack_list(filebnds(iproc)%rlen,1) = i
               filebnds(iproc)%unpack_list(filebnds(iproc)%rlen,2) = 0
               filebnds(iproc)%unpack_list(filebnds(iproc)%rlen,3) = n
               filebnds(iproc)%unpack_list(filebnds(iproc)%rlen,4) = ipf + 1
               iproc = procarray(i+ca,pjl+1+cb,no,1)
               floc = procarray(i+ca,pjl+1+cb,no,2)
               filebnds(iproc)%rlen = filebnds(iproc)%rlen + 1
               call check_filebnds_alloc(iproc,filemaxbuflen)
               ! store global index
               filebnds(iproc)%request_list(filebnds(iproc)%rlen,1) = i + ca
               filebnds(iproc)%request_list(filebnds(iproc)%rlen,2) = pjl + 1 + cb
               filebnds(iproc)%request_list(filebnds(iproc)%rlen,3) = no
               filebnds(iproc)%request_list(filebnds(iproc)%rlen,4) = floc
               ! store local index
               filebnds(iproc)%unpack_list(filebnds(iproc)%rlen,1) = i
               filebnds(iproc)%unpack_list(filebnds(iproc)%rlen,2) = pjl + 1
               filebnds(iproc)%unpack_list(filebnds(iproc)%rlen,3) = n
               filebnds(iproc)%unpack_list(filebnds(iproc)%rlen,4) = ipf + 1
            end do
            do j = 1,pjl
               iproc = procarray(ca,j+cb,no,1)
               floc = procarray(ca,j+cb,no,2)
               filebnds(iproc)%rlen = filebnds(iproc)%rlen + 1
               call check_filebnds_alloc(iproc,filemaxbuflen)
               ! store global index
               filebnds(iproc)%request_list(filebnds(iproc)%rlen,1) = ca
               filebnds(iproc)%request_list(filebnds(iproc)%rlen,2) = j + cb
               filebnds(iproc)%request_list(filebnds(iproc)%rlen,3) = no
               filebnds(iproc)%request_list(filebnds(iproc)%rlen,4) = floc
               ! store local index
               filebnds(iproc)%unpack_list(filebnds(iproc)%rlen,1) = 0
               filebnds(iproc)%unpack_list(filebnds(iproc)%rlen,2) = j
               filebnds(iproc)%unpack_list(filebnds(iproc)%rlen,3) = n
               filebnds(iproc)%unpack_list(filebnds(iproc)%rlen,4) = ipf + 1
               iproc = procarray(pil+1+ca,j+cb,no,1)
               floc = procarray(pil+1+ca,j+cb,no,2)
               filebnds(iproc)%rlen = filebnds(iproc)%rlen + 1
               call check_filebnds_alloc(iproc,filemaxbuflen)
               ! store global index
               filebnds(iproc)%request_list(filebnds(iproc)%rlen,1) = pil + 1 + ca
               filebnds(iproc)%request_list(filebnds(iproc)%rlen,2) = j + cb
               filebnds(iproc)%request_list(filebnds(iproc)%rlen,3) = no
               filebnds(iproc)%request_list(filebnds(iproc)%rlen,4) = floc
               ! store local index
               filebnds(iproc)%unpack_list(filebnds(iproc)%rlen,1) = pil + 1
               filebnds(iproc)%unpack_list(filebnds(iproc)%rlen,2) = j
               filebnds(iproc)%unpack_list(filebnds(iproc)%rlen,3) = n
               filebnds(iproc)%unpack_list(filebnds(iproc)%rlen,4) = ipf + 1
            end do
         end do
      end do
 
      ! identify neighbour processors
      fileneighnum = count( filebnds(:)%rlen > 0 )
      if ( filebnds(myid)%rlen > 0 ) then
        fileneighnum = fileneighnum - 1
      end if
      if ( allocated(fileneighlist) ) then
         deallocate(fileneighlist)
      end if
      allocate(fileneighlist(fileneighnum))
      ncount = 0
      do jproc = 1,fnresid-1
         iproc = modulo(myid+jproc,fnresid)
          if ( filebnds(iproc)%rlen > 0 ) then
            ncount = ncount + 1
            fileneighlist(ncount) = iproc
         end if
      end do

      ! reduce array size where possible
      allocate(dummy(filemaxbuflen,4))
      do jproc = 1,fileneighnum
         iproc = fileneighlist(jproc)
         xlen = filebnds(iproc)%rlen
         if ( filebnds(iproc)%len > xlen ) then
            dummy(1:xlen,1:4) = filebnds(iproc)%request_list(1:xlen,1:4)
            deallocate(filebnds(iproc)%request_list)
            allocate(filebnds(iproc)%request_list(xlen,4))
            filebnds(iproc)%request_list(1:xlen,1:4) = dummy(1:xlen,1:4)
            dummy(1:xlen,1:4) = filebnds(iproc)%unpack_list(1:xlen,1:4)
            deallocate(filebnds(iproc)%unpack_list)
            allocate(filebnds(iproc)%unpack_list(xlen,4))
            filebnds(iproc)%unpack_list(1:xlen,1:4) = dummy(1:xlen,1:4)
         end if
      end do
      if ( filebnds(myid)%rlen > 0 ) then
         xlen = filebnds(myid)%rlen
         if ( filebnds(myid)%len > xlen ) then
            dummy(1:xlen,1:4) = filebnds(myid)%request_list(1:xlen,1:4)
            deallocate(filebnds(myid)%request_list)
            allocate(filebnds(myid)%request_list(xlen,4))
            filebnds(myid)%request_list(1:xlen,1:4) = dummy(1:xlen,1:4)
            dummy(1:xlen,1:4) = filebnds(myid)%unpack_list(1:xlen,1:4)
            deallocate(filebnds(myid)%unpack_list)
            allocate(filebnds(myid)%unpack_list(xlen,4))
            filebnds(myid)%unpack_list(1:xlen,1:4) = dummy(1:xlen,1:4)
         end if
      end if
      deallocate(dummy)

      ! increase size of request list if needed
      if ( 2*fileneighnum > size(ireq) ) then
         deallocate( ireq, rlist )
         allocate( ireq(2*fileneighnum) )
         allocate( rlist(fileneighnum) )
      end if

      ! Now, for each processor send the length of points I want.
      nreq = 0
      do jproc = 1,fileneighnum
         iproc = fileneighlist(jproc) ! Recv from
         nreq = nreq + 1
         lproc = iproc
         call MPI_IRecv( filebnds(iproc)%slen, 1_4, ltype, lproc, itag, lcomm, ireq(nreq), ierr )
      end do
      do jproc = fileneighnum,1,-1
         iproc = fileneighlist(jproc)
         nreq = nreq + 1
         lproc = iproc
         call MPI_ISend( filebnds(iproc)%rlen, 1_4, ltype, lproc, itag, lcomm, ireq(nreq), ierr )
      end do
      call MPI_Waitall( nreq, ireq, MPI_STATUSES_IGNORE, ierr )
      
      ! send unpack_list to file neighbours
      nreq = 0
      do jproc = 1,fileneighnum
         iproc = fileneighlist(jproc)  ! Recv from
         allocate(filebnds(iproc)%send_list(filebnds(iproc)%slen,4))
         nreq = nreq + 1
         ! Use the maximum size in the recv call.
         llen = 4*filebnds(iproc)%slen
         lproc = iproc
         call MPI_IRecv( filebnds(iproc)%send_list(1,1), llen, &
              ltype, lproc, itag, lcomm, ireq(nreq), ierr )
      end do
      do jproc = fileneighnum,1,-1
         iproc = fileneighlist(jproc)  ! Send to
         ! Send list of requests
         nreq = nreq + 1
         llen = 4*filebnds(iproc)%rlen
         lproc = iproc
         call MPI_ISend( filebnds(iproc)%request_list(1,1), llen, &
              ltype, lproc, itag, lcomm, ireq(nreq), ierr )
      end do      
      call MPI_Waitall( nreq, ireq, MPI_STATUSES_IGNORE, ierr )

      ! convert send_list and unpack_list to local indices
      do jproc = 1,fileneighnum
         iproc = fileneighlist(jproc)  ! Send to
         do iq = 1,filebnds(iproc)%slen
            iloc = filebnds(iproc)%send_list(iq,1)
            jloc = filebnds(iproc)%send_list(iq,2)
            nloc = filebnds(iproc)%send_list(iq,3)
            floc = filebnds(iproc)%send_list(iq,4)
            call file_ijnpg2ijnp(iloc,jloc,nloc,floc,myid,ik)
            filebnds(iproc)%send_list(iq,1) = iloc
            filebnds(iproc)%send_list(iq,2) = jloc
            filebnds(iproc)%send_list(iq,3) = nloc
            filebnds(iproc)%send_list(iq,4) = floc
         end do
      end do
      do iq = 1,filebnds(myid)%rlen
         iloc = filebnds(myid)%request_list(iq,1)
         jloc = filebnds(myid)%request_list(iq,2)
         nloc = filebnds(myid)%request_list(iq,3)
         floc = filebnds(myid)%request_list(iq,4)
         call file_ijnpg2ijnp(iloc,jloc,nloc,floc,myid,ik)
         filebnds(myid)%request_list(iq,1) = iloc
         filebnds(myid)%request_list(iq,2) = jloc
         filebnds(myid)%request_list(iq,3) = nloc
         filebnds(myid)%request_list(iq,4) = floc
      end do
   
      ! set up buffers
      xlev = max( kl, ol )
      do jproc = 1,fileneighnum
         iproc = fileneighlist(jproc)
         xlen = xlev*filebnds(iproc)%rlen
         if ( bnds(iproc)%rbuflen < xlen ) then
            if ( bnds(iproc)%rbuflen > 0 ) then
               deallocate( bnds(iproc)%rbuf )
               deallocate( bnds(iproc)%r8buf )
            end if
            allocate( bnds(iproc)%rbuf(xlen) )
            allocate( bnds(iproc)%r8buf(xlen) )
            bnds(iproc)%rbuflen = xlen
         end if
         xlen = xlev*filebnds(iproc)%slen
         if ( bnds(iproc)%sbuflen < xlen ) then
            if ( bnds(iproc)%sbuflen > 0 ) then
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
      end if
   
   end subroutine check_filebnds_alloc
   
   subroutine file_ijnpg2ijnp(iloc,jloc,nloc,floc,iproc,ik)
   
      integer, intent(inout) :: iloc, jloc, nloc
      integer, intent(in) :: floc, iproc,ik
      integer :: i, j, n
      integer :: ip, ca, cb
  
      i = iloc ! global i
      j = jloc ! global j
      n = nloc ! global panel
      
      ! fix up out-of-bounds indices
      if ( iloc == 0 ) then
         if ( mod(nloc,2) == 0 ) then          
            i = ik
            j = jloc
            n = mod(nloc+5,6) ! n_w
         else
            i = ik + 1 - jloc
            j = ik
            n = mod(nloc+4,6) ! n_w
         end if
      else if ( iloc == ik+1 ) then
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
         
      ip = (floc-1)*fnresid + iproc
      ca = pioff(ip,n)
      cb = pjoff(ip,n)
      iloc = i - ca        ! local i
      jloc = j - cb        ! local j
      nloc = n + pnoff(ip) ! local n
   
   end subroutine file_ijnpg2ijnp
   
   subroutine ccmpi_filebounds2(sdat,comm)

      integer, intent(in) :: comm
      integer :: myrlen, iproc, jproc, mproc, iq, rcount
      integer :: send_len
      integer, dimension(fileneighnum) :: rslen, sslen
      integer(kind=4) :: llen, lproc, ierr, ldone, sreq, lcomm
      integer(kind=4) :: itag=40
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      integer(kind=4), dimension(fileneighnum) :: donelist
      real, dimension(0:,0:,1:,1:), intent(inout) :: sdat

      call START_LOG(bounds_begin)
      
      lcomm = comm
      
      rslen(:) = filebnds(fileneighlist)%rlen
      sslen(:) = filebnds(fileneighlist)%slen
      myrlen = filebnds(myid)%rlen

      !     Set up the buffers to recv
      nreq = 0
      do iproc = 1,fileneighnum
         llen = rslen(iproc)
         lproc = fileneighlist(iproc)  ! Recv from
         nreq  = nreq + 1
         rlist(nreq) = iproc
         call MPI_IRecv( bnds(lproc)%rbuf(1), llen, ltype, lproc, itag, lcomm, ireq(nreq), ierr )
      end do
      rreq = nreq
      !     Set up the buffers to send
      do iproc = fileneighnum,1,-1
         ! Build up list of points
         send_len = sslen(iproc)
         llen = send_len
         lproc = fileneighlist(iproc)  ! Send to
         do iq = 1,send_len
            bnds(lproc)%sbuf(iq) = sdat(filebnds(lproc)%send_list(iq,1),filebnds(lproc)%send_list(iq,2), &
                                        filebnds(lproc)%send_list(iq,3),filebnds(lproc)%send_list(iq,4))
         end do
         nreq  = nreq + 1
         call MPI_ISend( bnds(lproc)%sbuf(1), llen, ltype, lproc, itag, lcomm, ireq(nreq), ierr )
      end do

      ! Finally see if there are any points on my own processor that need
      ! to be fixed up.
      do iq = 1,myrlen
         ! request_list is same as send_list in this case
         sdat(filebnds(myid)%unpack_list(iq,1),filebnds(myid)%unpack_list(iq,2),   &
              filebnds(myid)%unpack_list(iq,3),filebnds(myid)%unpack_list(iq,4)) = &
         sdat(filebnds(myid)%request_list(iq,1),filebnds(myid)%request_list(iq,2), &
              filebnds(myid)%request_list(iq,3),filebnds(myid)%request_list(iq,4))
      end do

      ! Unpack incomming messages
      rcount = rreq
      do while ( rcount > 0 )
         call START_LOG(mpiwait_begin)
         call MPI_Waitsome( rreq, ireq, ldone, donelist, MPI_STATUSES_IGNORE, ierr )
         call END_LOG(mpiwait_end)
         rcount = rcount - ldone
         do jproc = 1,ldone
            mproc = donelist(jproc)
            iproc = rlist(mproc)  ! Recv from
            lproc = fileneighlist(iproc)
            ! unpack_list(iq) is index into extended region
!cdir nodep
            do iq = 1,rslen(iproc)
               sdat(filebnds(lproc)%unpack_list(iq,1),filebnds(lproc)%unpack_list(iq,2),   &
                    filebnds(lproc)%unpack_list(iq,3),filebnds(lproc)%unpack_list(iq,4)) = &
               bnds(lproc)%rbuf(iq)
            end do
         end do
      end do

      ! Clear any remaining messages
      sreq = nreq - rreq
      call START_LOG(mpiwait_begin)
      call MPI_Waitall( sreq, ireq(rreq+1:nreq), MPI_STATUSES_IGNORE, ierr )
      call END_LOG(mpiwait_end)

      call END_LOG(bounds_end)

   end subroutine ccmpi_filebounds2

   subroutine ccmpi_filebounds3(sdat,comm)

      integer, intent(in) :: comm
      integer :: myrlen, iproc, jproc, mproc, iq, rcount, kx
      integer :: send_len
      integer, dimension(fileneighnum) :: rslen, sslen
      integer(kind=4) :: llen, lproc, ierr, ldone, sreq, lcomm
      integer(kind=4) :: itag=41
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      integer(kind=4), dimension(fileneighnum) :: donelist
      real, dimension(0:,0:,1:,1:,1:), intent(inout) :: sdat

      call START_LOG(bounds_begin)
      
      kx = size(sdat,5)
      lcomm = comm

      rslen(:) = filebnds(fileneighlist)%rlen
      sslen(:) = filebnds(fileneighlist)%slen
      myrlen = filebnds(myid)%rlen

      !     Set up the buffers to send and recv
      nreq = 0
      do iproc = 1,fileneighnum
         llen = rslen(iproc)*kx
         lproc = fileneighlist(iproc)  ! Recv from
         nreq = nreq + 1
         rlist(nreq) = iproc
         call MPI_IRecv( bnds(lproc)%rbuf(1), llen, ltype, lproc, itag, lcomm, ireq(nreq), ierr )
      end do
      rreq = nreq
      do iproc = fileneighnum,1,-1
         send_len = sslen(iproc)
         llen = send_len*kx
         lproc = fileneighlist(iproc)  ! Send to
         do iq = 1,send_len
            bnds(lproc)%sbuf(1+(iq-1)*kx:iq*kx) = sdat(filebnds(lproc)%send_list(iq,1),filebnds(lproc)%send_list(iq,2),      &
                                                       filebnds(lproc)%send_list(iq,3),filebnds(lproc)%send_list(iq,4),1:kx)
         end do
         nreq = nreq + 1
         call MPI_ISend( bnds(lproc)%sbuf(1), llen, ltype, lproc, itag, lcomm, ireq(nreq), ierr )
      end do

      ! Finally see if there are any points on my own processor that need
      ! to be fixed up.
      do iq = 1,myrlen
         ! request_list is same as send_list in this case
         sdat(filebnds(myid)%unpack_list(iq,1),filebnds(myid)%unpack_list(iq,2),        &
              filebnds(myid)%unpack_list(iq,3),filebnds(myid)%unpack_list(iq,4),1:kx) = &
         sdat(filebnds(myid)%request_list(iq,1),filebnds(myid)%request_list(iq,2),      &
              filebnds(myid)%request_list(iq,3),filebnds(myid)%request_list(iq,4),1:kx)
      end do

      ! Unpack incomming messages
      rcount = rreq
      do while ( rcount > 0 )
         call START_LOG(mpiwait_begin)
         call MPI_Waitsome( rreq, ireq, ldone, donelist, MPI_STATUSES_IGNORE, ierr )
         call END_LOG(mpiwait_end)
         rcount = rcount - ldone
         do jproc = 1,ldone
            mproc = donelist(jproc)
            iproc = rlist(mproc)  ! Recv from
            lproc = fileneighlist(iproc)
            do iq = 1,rslen(iproc)
               sdat(filebnds(lproc)%unpack_list(iq,1),filebnds(lproc)%unpack_list(iq,2),        &
                    filebnds(lproc)%unpack_list(iq,3),filebnds(lproc)%unpack_list(iq,4),1:kx) = &
               bnds(lproc)%rbuf(1+(iq-1)*kx:iq*kx)
            end do
         end do
      end do

      ! Clear any remaining messages
      sreq = nreq - rreq
      call START_LOG(mpiwait_begin)
      call MPI_Waitall( sreq, ireq(rreq+1:nreq), MPI_STATUSES_IGNORE, ierr )
      call END_LOG(mpiwait_end)

      call END_LOG(bounds_end)
      
   end subroutine ccmpi_filebounds3

#ifdef usempi3   
   subroutine ccmpi_allocshdata2r(pdata,sshape,win)
      use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer

      real, pointer, dimension(:), intent(inout) :: pdata 
      integer, intent(out) :: win
      integer, dimension(1), intent(in) :: sshape
      integer(kind=MPI_ADDRESS_KIND) :: qsize, lsize
      integer(kind=4) :: disp_unit, ierr, tsize
      integer(kind=4) :: lcomm, lwin
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      type(c_ptr) :: baseptr

!     allocted a single shared memory region on each node
      lcomm = comm_node
      call MPI_Type_size( ltype, tsize, ierr )
      if ( node_myid==0 ) then
         lsize = sshape(1)*tsize
      else
         lsize = 0_4
      end if
      call MPI_Win_allocate_shared( lsize, 1_4, MPI_INFO_NULL, lcomm, baseptr, lwin, ierr )
      if ( node_myid/=0 ) then
         call MPI_Win_shared_query( lwin, 0_4, qsize, disp_unit, baseptr, ierr )
      end if
      call c_f_pointer( baseptr, pdata, sshape )
      win = lwin

   end subroutine ccmpi_allocshdata2r 

   subroutine ccmpi_allocshdata4r(pdata,sshape,win)
      use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer

      real, pointer, dimension(:,:,:), intent(inout) :: pdata 
      integer, dimension(3), intent(in) :: sshape
      integer, intent(out) :: win
      integer(kind=MPI_ADDRESS_KIND) :: qsize, lsize
      integer(kind=4) :: disp_unit, ierr, tsize
      integer(kind=4) :: lcomm, lwin
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      type(c_ptr) :: baseptr

!     allocted a single shared memory region on each node
      lcomm = comm_node
      call MPI_Type_size(ltype, tsize, ierr)
      if ( node_myid==0 ) then
         lsize = sshape(1)*sshape(2)*sshape(3)*tsize
      else
         lsize = 0_4
      end if
      call MPI_Win_allocate_shared(lsize, 1_4, MPI_INFO_NULL, lcomm, baseptr, lwin, ierr)
      if ( node_myid/=0 ) then
         call MPI_Win_shared_query(lwin, 0_4, qsize, disp_unit, baseptr, ierr)
      end if
      call c_f_pointer(baseptr, pdata, sshape)
      win = lwin

   end subroutine ccmpi_allocshdata4r 

   subroutine ccmpi_allocshdata5r(pdata,sshape,win)
      use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer

      real, pointer, dimension(:,:,:,:), intent(inout) :: pdata 
      integer, intent(out) :: win
      integer, dimension(4), intent(in) :: sshape
      integer(kind=MPI_ADDRESS_KIND) :: qsize, lsize
      integer(kind=4) :: disp_unit, ierr, tsize
      integer(kind=4) :: lcomm, lwin
#ifdef i8r8
      integer(kind=4), parameter :: ltype = MPI_DOUBLE_PRECISION
#else
      integer(kind=4), parameter :: ltype = MPI_REAL
#endif
      type(c_ptr) :: baseptr

!     allocted a single shared memory region on each node
      lcomm = comm_node
      call MPI_Type_size( ltype, tsize, ierr )
      if ( node_myid==0 ) then
         lsize = sshape(1)*sshape(2)*sshape(3)*sshape(4)*tsize
      else
         lsize = 0_4
      end if
      call MPI_Win_allocate_shared( lsize, 1_4, MPI_INFO_NULL, lcomm, baseptr, lwin, ierr )
      if ( node_myid/=0 ) then
         call MPI_Win_shared_query( lwin, 0_4, qsize, disp_unit, baseptr, ierr )
      end if
      call c_f_pointer( baseptr, pdata, sshape )
      win = lwin

   end subroutine ccmpi_allocshdata5r 

   subroutine ccmpi_allocshdata2i(pdata,sshape,win)
      use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer

      integer, pointer, dimension(:), intent(inout) :: pdata 
      integer, intent(out) :: win
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
      lcomm = comm_node
      call MPI_Type_size( ltype, tsize, ierr )
      if ( node_myid==0 ) then
         lsize = sshape(1)*tsize
      else
         lsize = 0_4
      end if
      call MPI_Win_allocate_shared( lsize, 1_4, MPI_INFO_NULL, lcomm, baseptr, lwin, ierr )
      if ( node_myid/=0 ) then
         call MPI_Win_shared_query( lwin, 0_4, qsize, disp_unit, baseptr, ierr )
      end if
      call c_f_pointer( baseptr, pdata, sshape )
      win = lwin

   end subroutine ccmpi_allocshdata2i
   
   subroutine ccmpi_allocshdata5i(pdata,sshape,win)
      use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer

      integer, pointer, dimension(:,:,:,:), intent(inout) :: pdata 
      integer, intent(out) :: win
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
      lcomm = comm_node
      call MPI_Type_size( ltype, tsize, ierr )
      if ( node_myid==0 ) then
         lsize = sshape(1)*sshape(2)*sshape(3)*sshape(4)*tsize
      else
         lsize = 0_4
      end if
      call MPI_Win_allocate_shared( lsize, 1_4, MPI_INFO_NULL, lcomm, baseptr, lwin, ierr )
      if ( node_myid/=0 ) then
         call MPI_Win_shared_query( lwin, 0_4, qsize, disp_unit, baseptr, ierr )
      end if
      call c_f_pointer( baseptr, pdata, sshape )
      win = lwin

   end subroutine ccmpi_allocshdata5i
   
   subroutine ccmpi_allocshdata2_r8(pdata,sshape,win)
      use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer

      real(kind=8), pointer, dimension(:), intent(inout) :: pdata 
      integer, intent(out) :: win
      integer, dimension(1), intent(in) :: sshape
      integer(kind=MPI_ADDRESS_KIND) :: qsize, lsize
      integer(kind=4) :: disp_unit, ierr, tsize
      integer(kind=4) :: lcomm, lwin
      type(c_ptr) :: baseptr

!     allocted a single shared memory region on each node
      lcomm = comm_node
      call MPI_Type_size( MPI_DOUBLE_PRECISION, tsize, ierr )
      if ( node_myid==0 ) then
         lsize = sshape(1)*tsize
      else
         lsize = 0_4
      end if
      call MPI_Win_allocate_shared( lsize, 1_4, MPI_INFO_NULL, lcomm, baseptr, lwin, ierr )
      if ( node_myid/=0 ) then
         call MPI_Win_shared_query( lwin, 0_4, qsize, disp_unit, baseptr, ierr )
      end if
      call c_f_pointer( baseptr, pdata, sshape )
      win = lwin

   end subroutine ccmpi_allocshdata2_r8
   
   subroutine ccmpi_allocshdata3_r8(pdata,sshape,win)
      use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer

      real(kind=8), pointer, dimension(:,:), intent(inout) :: pdata 
      integer, intent(out) :: win
      integer, dimension(2), intent(in) :: sshape
      integer(kind=MPI_ADDRESS_KIND) :: qsize, lsize
      integer(kind=4) :: disp_unit, ierr, tsize
      integer(kind=4) :: lcomm, lwin
      type(c_ptr) :: baseptr

!     allocted a single shared memory region on each node
      lcomm = comm_node
      call MPI_Type_size( MPI_DOUBLE_PRECISION, tsize, ierr )
      if ( node_myid==0 ) then
         lsize = sshape(1)*sshape(2)*tsize
      else
         lsize = 0_4
      end if
      call MPI_Win_allocate_shared( lsize, 1_4, MPI_INFO_NULL, lcomm, baseptr, lwin, ierr )
      if ( node_myid/=0 ) then
         call MPI_Win_shared_query( lwin, 0_4, qsize, disp_unit, baseptr, ierr )
      end if
      call c_f_pointer( baseptr, pdata, sshape )
      win = lwin

   end subroutine ccmpi_allocshdata3_r8 

   subroutine ccmpi_allocshdata4_r8(pdata,sshape,win)
      use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer

      real(kind=8), pointer, dimension(:,:,:), intent(inout) :: pdata 
      integer, intent(out) :: win
      integer, dimension(3), intent(in) :: sshape
      integer(kind=MPI_ADDRESS_KIND) :: qsize, lsize
      integer(kind=4) :: disp_unit, ierr, tsize
      integer(kind=4) :: lcomm, lwin
      type(c_ptr) :: baseptr

!     allocted a single shared memory region on each node
      lcomm = comm_node
      call MPI_Type_size( MPI_DOUBLE_PRECISION, tsize, ierr )
      if ( node_myid==0 ) then
         lsize = sshape(1)*sshape(2)*sshape(3)*tsize
      else
         lsize = 0_4
      end if
      call MPI_Win_allocate_shared( lsize, 1_4, MPI_INFO_NULL, lcomm, baseptr, lwin, ierr )
      if ( node_myid/=0 ) then
         call MPI_Win_shared_query( lwin, 0_4, qsize, disp_unit, baseptr, ierr )
      end if
      call c_f_pointer( baseptr, pdata, sshape )
      win = lwin

   end subroutine ccmpi_allocshdata4_r8 

   subroutine ccmpi_shepoch(win)
   
       integer, intent(in) :: win
       integer(kind=4) :: lwin, lerr
       
       lwin = win
       call MPI_Win_fence( 0_4, lwin, lerr )

   end subroutine ccmpi_shepoch
   
   subroutine ccmpi_freeshdata(win)
   
      integer, intent(in) :: win
      integer(kind=4) :: lwin, lerr
      
      lwin = win
      call MPI_win_free( lwin, lerr ) 
      
   end subroutine ccmpi_freeshdata
#endif
   
end module cc_mpi

