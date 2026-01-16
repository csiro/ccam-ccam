! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2026 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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

! This module manages all MPI communications between processes.  The system was originally developed
! by Martin Dix and subsequently modified by Marcus Thatcher.  Thanks to Aaron McDonough for developing
! the Vampir trace routines and upgrading the timer calls.  Thanks to Paul Ryan for the design of the
! shared memory arrays.

! Preprocess directives:

! -Dusempi3        exploits MPI-3 shared memory to optmise communication between nodes
! -Dusempimod      uses f90 mod interface for MPI.
! -Dusempimod_f09  uses f09 mod interface for MPI.  Requires -Dusempimod
! -Dvampir         is for coupling with VAMPIR for tracers
! -Di8r8           is for running in double precision mode

   use cc_mpi_collective
   use cc_mpi_common
   use cc_mpi_file
   use cc_mpi_multigrid
   use cc_mpi_point
   
   private
   
   integer, allocatable, dimension(:), save, public :: neighlist           ! list of neighbour processes
   integer, allocatable, dimension(:), save, private :: neighmap           ! map of process to neighbour index
   integer, save, public :: neighnum                                       ! number of neighbours

   public :: ccmpi_setup, ccmpi_init, ccmpi_reinit, ccmpi_finalize,         &
             ccmpi_procformat_init, ccmpi_remap 

   ! Imported from cc_mpi_collective
   public :: ccmpi_distribute, ccmpi_gather, ccmpi_gatherr8,                &
             ccmpi_distributer8, ccmpi_gatherall
   public :: allocateglobalpack, copyglobalpack,                            &
             ccmpi_gathermap_send, ccmpi_gathermap_recv, getglobalpack_v,   &
             setglobalpack_v, ccmpi_gathermap_wait, deallocateglobalpack
   public :: ccglobal_posneg   
   public :: readglobvar, writeglobvar
   public :: specmap_req, specmap_recv, specmap_send, specmap_indx,         &
             specmap_indxlen, specmap_ext
   
   ! Imported from cc_mpi_common
   public :: comm_world
   public :: comm_node, node_myid, node_nproc
   public :: comm_nodecaptain, nodecaptain_myid, nodecaptain_nproc,         &
             nodecaptainid
   public :: comm_proc, comm_rows, comm_cols, hproc, mproc, npta, pprocn,   &
             pprocx
   public :: comm_vnode, vnode_nproc, vnode_myid
   public :: comm_vleader, vleader_nproc, vleader_myid, vnode_vleaderid
   public :: myid, ipan, jpan, ioff, joff, noff, nxproc, nyproc, node_nx,   &
             node_ny, nagg, maxcolour, mydiag
   public :: indp, indg, iq2iqg, iqg2iq, indv_mpi, fproc, face_set,         &
             uniform_set, dix_set
   public :: start_log, end_log, log_on, log_off, log_flush, log_setup
   public :: ccmpi_reduce, ccmpi_reducer8, ccmpi_allreduce, ccmpi_abort,    &
             ccmpi_bcast, ccmpi_bcastr8, ccmpi_barrier, ccmpi_gatherx,      &
             ccmpi_gatherxr8, ccmpi_scatterx, ccmpi_allgatherx,             &
             ccmpi_commsplit, ccmpi_commfree, ccmpi_alltoall
   public :: edge_w, edge_n, edge_s, edge_e
   public :: ifull_maxcolour, iqx, ifull_colour, ifull_colour_border
   public :: simple_timer_finalize      
   public :: ints_begin, ints_end
   public :: nonlin_begin, nonlin_end
   public :: helm_begin, helm_end
   public :: adjust_begin, adjust_end
   public :: upglobal_begin, upglobal_end
   public :: hordifg_begin, hordifg_end
   public :: vadv_begin, vadv_end
   public :: depts_begin, depts_end
   public :: stag_begin, stag_end
   public :: ocnstag_begin, ocnstag_end
   public :: mfix_begin, mfix_end
   public :: phys_begin, phys_end
   public :: outfile_begin, outfile_end
   public :: onthefly_begin, onthefly_end
   public :: otf_fill_begin, otf_fill_end
   public :: otf_ints_begin, otf_ints_end
   public :: histrd_begin, histrd_end
   public :: indata_begin, indata_end
   public :: nestin_begin, nestin_end
   public :: nestOTF_begin, nestOTF_end
   public :: nestWIN_begin, nestWIN_end
   public :: nestcalc_begin, nestcalc_end
   public :: nestcomm_begin, nestcomm_end
   public :: ensemble_begin, ensemble_end
   public :: amipsst_begin, amipsst_end
   public :: gwdrag_begin, gwdrag_end
   public :: convection_begin, convection_end
   public :: cloud_begin, cloud_end
   public :: radnet_begin, radnet_end
   public :: radinit_begin, radinit_end
   public :: radSW_begin, radSW_end
   public :: radLW_begin, radLW_end
   public :: sfluxnet_begin, sfluxnet_end
   public :: sfluxwater_begin, sfluxwater_end
   public :: sfluxland_begin, sfluxland_end
   public :: sfluxurban_begin, sfluxurban_end
   public :: vertmix_begin, vertmix_end
   public :: aerosol_begin, aerosol_end
   public :: cape_begin, cape_end
   public :: maincalc_begin, maincalc_end
   public :: precon_begin, precon_end
   public :: waterdynamics_begin, waterdynamics_end
   public :: watermfix_begin, watermfix_end
   public :: waterdeps_begin, waterdeps_end
   public :: watereos_begin, watereos_end
   public :: waterints_begin, waterints_end
   public :: watervadv_begin, watervadv_end
   public :: waterhelm_begin, waterhelm_end
   public :: wateriadv_begin, wateriadv_end
   public :: waterdiff_begin, waterdiff_end
   public :: river_begin, river_end
   public :: bcast_begin, bcast_end
   public :: alltoall_begin, alltoall_end
   public :: allgather_begin, allgather_end
   public :: gather_begin, gather_end
   public :: scatter_begin, scatter_end
   public :: reduce_begin, reduce_end
   public :: allreduce_begin, allreduce_end
   public :: mpiwait_begin, mpiwait_end
   public :: mpibarrier_begin, mpibarrier_end
   public :: mgfine_begin, mgfine_end
   public :: mgup_begin, mgup_end
   public :: mgcoarse_begin, mgcoarse_end
   public :: mgdown_begin, mgdown_end
   public :: p1_begin, p1_end
   public :: p2_begin, p2_end
   public :: p3_begin, p3_end
   public :: p4_begin, p4_end
   public :: p5_begin, p5_end
   public :: p6_begin, p6_end
   public :: p7_begin, p7_end
   public :: p8_begin, p8_end
   public :: p9_begin, p9_end
   public :: p10_begin, p10_end
   public :: p11_begin, p11_end
   public :: p12_begin, p12_end
   public :: p13_begin, p13_end
   public :: p14_begin, p14_end
   public :: p15_begin, p15_end
   public :: p16_begin, p16_end
   public :: p17_begin, p17_end   
   public :: mpiinit_time, total_time
#ifdef usempi3
   public :: ccmpi_allocshdata, ccmpi_allocshdatar8
   public :: ccmpi_shepoch, ccmpi_freeshdata
#endif

   ! Imported from cc_mpi_multigrid
   public :: mgbounds, mgcollect, mgbcast, mg_index, mg_fproc, mg_fproc_1,  &
             mgbounds_colour
   public :: mg, mg_bnds, mg_maxlevel_local, mg_maxlevel

   ! Imported from cc_mpi_file
   public :: ccmpi_filewinget, ccmpi_filebounds_setup,                      &
             ccmpi_filebounds, ccmpi_filedistribute, procarray,             &
             ccmpi_filewininit, ccmpi_filewinfinalize,                      &
             ccmpi_filewinfinalize_exit
   public :: pipan, pjpan, pnpan
   public :: pil_g, pjl_g, pka_g, pko_g
   public :: fnproc, fnresid, fncount, mynproc
   public :: pnoff, pioff, pjoff
   public :: filemap_req, filemap_qmod, filemap_recv, filemap_rmod,         &
             filemap_send, filemap_smod, filemap_indx, filemap_indxlen,     &
             nodefile_win, nodefilesave_win, nodefile_count
   public :: fileneighnum

   ! Imported from cc_mpi_point
   public :: bounds, boundsuv, bounds_colour_send, bounds_colour_recv,      &
             bounds_send, bounds_recv, boundsr8, deptsync, intssync_send,   &
             intssync_recv
   public :: drlen, dpoints, sextra

contains

   subroutine ccmpi_setup(id,jd,idjd,dt)
      use indices_m
      use latlong_m
      use map_m
      use sumdd_m
      use vecsuv_m
      use workglob_m
      use xyzinfo_m
      integer, intent(in) :: id, jd
      integer, intent(out) :: idjd
      integer :: iproc, dproc, iq, iqg, i, j, n
      integer :: colourmask
      integer(kind=4) :: ierr, colour, rank, lcommin, lcommout
      real, intent(in) :: dt
      real, dimension(:,:), allocatable :: dum
      real, dimension(:,:), allocatable :: dumu, dumv
      real, dimension(:,:), allocatable :: dumr, dumr_g
      real(kind=8), dimension(:,:), allocatable :: dumr8, dumr8_g
      logical(kind=4) :: ltrue

      
      if ( maxcolour<1 .and. maxcolour>3 ) then
         write(6,*) "ERROR: Invalid choice for maxcolour"
         call ccmpi_abort(-1)
      end if
      
      
      ! Allocate arrays for boundary information
      nreq = 0
      allocate( bnds(0:nproc-1) )
      ! index=0 is for all coloured grid points
      do n = 0,nproc-1
         allocate( bnds(n)%rlenh_bg(0:maxcolour) )
         allocate( bnds(n)%rlenh_fn(0:maxcolour) )
         allocate( bnds(n)%slenh_bg(0:maxcolour) )
         allocate( bnds(n)%slenh_fn(0:maxcolour) )
         allocate( bnds(n)%rlen_bg(0:maxcolour) )
         allocate( bnds(n)%rlen_fn(0:maxcolour) )
         allocate( bnds(n)%slen_bg(0:maxcolour) )
         allocate( bnds(n)%slen_fn(0:maxcolour) )
         allocate( bnds(n)%rlenx_bg(0:maxcolour) )
         allocate( bnds(n)%rlenx_fn(0:maxcolour) )
         allocate( bnds(n)%slenx_bg(0:maxcolour) )
         allocate( bnds(n)%slenx_fn(0:maxcolour) )
      end do  


      ! Decompose grid over processes
      call proc_setup(id,jd,idjd)
      if ( nproc < npanels+1 ) then
         ! possible to have two boundaries from the same process 
         maxbuflen = (il_g+4)*3*2*2*npan + 4
      else
         ! only one boundary can be sent from a process 
         maxbuflen = (max(ipan,jpan)+4)*3*2 + 4
      end if    
      maxvertlen = max( kl, ol, 15 )

      !$acc update device(ipan,jpan)
      
      ! Distribute global arrays over processes
      if ( myid==0 ) then
         write(6,*) "-> Distribute global arrays"
      end if
      allocate( dumr(ifull,23) )
      ! combine 2d arrays into a single 3d arrary for distribute
      if ( myid == 0 ) then
         allocate( dumr_g(ifull_g,23) ) 
         dumr_g(1:ifull_g,1)  = wts_g(1:ifull_g)
         dumr_g(1:ifull_g,2)  = em_g(1:ifull_g)
         dumr_g(1:ifull_g,3)  = emu_g(1:ifull_g)
         dumr_g(1:ifull_g,4)  = emv_g(1:ifull_g)
         dumr_g(1:ifull_g,5)  = ax_g(1:ifull_g)
         dumr_g(1:ifull_g,6)  = ay_g(1:ifull_g)
         dumr_g(1:ifull_g,7)  = az_g(1:ifull_g)
         dumr_g(1:ifull_g,8)  = bx_g(1:ifull_g)
         dumr_g(1:ifull_g,9)  = by_g(1:ifull_g)
         dumr_g(1:ifull_g,10) = bz_g(1:ifull_g)
         dumr_g(1:ifull_g,11) = f_g(1:ifull_g)
         dumr_g(1:ifull_g,12) = fu_g(1:ifull_g)
         dumr_g(1:ifull_g,13) = fv_g(1:ifull_g)
         dumr_g(1:ifull_g,14) = rlatt_g(1:ifull_g)
         dumr_g(1:ifull_g,15) = rlongg_g(1:ifull_g)
         dumr_g(1:ifull_g,16:19) = rlat4(1:ifull_g,1:4)
         dumr_g(1:ifull_g,20:23) = rlong4(1:ifull_g,1:4)
         call ccmpi_distribute(dumr(:,1:23),dumr_g(:,1:23)) 
         deallocate( dumr_g )
      else
         call ccmpi_distribute(dumr(:,1:23))
      end if
      ! update 3d array into 2d arrays
      wts(1:ifull)    = dumr(1:ifull,1)
      em(1:ifull)     = dumr(1:ifull,2)
      emu(1:ifull)    = dumr(1:ifull,3)
      emv(1:ifull)    = dumr(1:ifull,4)
      ax(1:ifull)     = dumr(1:ifull,5)
      ay(1:ifull)     = dumr(1:ifull,6)
      az(1:ifull)     = dumr(1:ifull,7)
      bx(1:ifull)     = dumr(1:ifull,8)
      by(1:ifull)     = dumr(1:ifull,9)
      bz(1:ifull)     = dumr(1:ifull,10)
      f(1:ifull)      = dumr(1:ifull,11)
      fu(1:ifull)     = dumr(1:ifull,12)
      fv(1:ifull)     = dumr(1:ifull,13)
      rlatt(1:ifull)  = dumr(1:ifull,14)
      rlongg(1:ifull) = dumr(1:ifull,15)
      rlat4_l(1:ifull,1:4)  = dumr(1:ifull,16:19)
      rlong4_l(1:ifull,1:4) = dumr(1:ifull,20:23)
      deallocate( dumr )
      allocate( dumr8(ifull,3) )
      ! combine 2d arrays into a single 3d arrary for distribute
      if ( myid == 0 ) then
         allocate( dumr8_g(ifull_g,3) ) 
         dumr8_g(1:ifull_g,1) = x_g(1:ifull_g)
         dumr8_g(1:ifull_g,2) = y_g(1:ifull_g)
         dumr8_g(1:ifull_g,3) = z_g(1:ifull_g)
         call ccmpi_distributer8(dumr8(:,1:3),dumr8_g(:,1:3)) 
         deallocate( dumr8_g )
      else
         call ccmpi_distributer8(dumr8(:,1:3))
      end if
      ! update 3d array into 2d arrays
      x(1:ifull) = dumr8(1:ifull,1)
      y(1:ifull) = dumr8(1:ifull,2)
      z(1:ifull) = dumr8(1:ifull,3)
      deallocate( dumr8 )

      
      ! Configure halos
      call bounds_setup(dt)
      if ( myid==0 ) then
         write(6,*) "-> Update halos for constant arrays"
      end if
      allocate( dum(1:ifull+iextra,2) )
      dum = 0.
      dum(1:ifull,1) = em(1:ifull)
      dum(1:ifull,2) = f(1:ifull)
      call bounds(dum(:,1:2),corner=.true.)
      em(ifull+1:ifull+iextra) = dum(ifull+1:ifull+iextra,1)
      f(ifull+1:ifull+iextra) = dum(ifull+1:ifull+iextra,2)
      deallocate( dum )
      allocate( dumr8(ifull+iextra,3) )
      dumr8(1:ifull,1) = x(1:ifull)
      dumr8(1:ifull,2) = y(1:ifull)
      dumr8(1:ifull,3) = z(1:ifull)
      call boundsr8(dumr8(:,1:3),corner=.true.)
      x(ifull+1:ifull+iextra) = dumr8(ifull+1:ifull+iextra,1)
      y(ifull+1:ifull+iextra) = dumr8(ifull+1:ifull+iextra,2)
      z(ifull+1:ifull+iextra) = dumr8(ifull+1:ifull+iextra,3)
      deallocate( dumr8 )
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
      
      
      ! Off process departure points
      if ( myid==0 ) then
         write(6,*) "-> Allocate memory for departure points"
      end if
      allocate( dpoints(neighnum) )
      allocate( dbuf(0:neighnum) )
      allocate( dindex(0:neighnum) )
      allocate( sextra(neighnum) )
      allocate( dslen(0:neighnum), drlen(0:neighnum) )
      dslen(:) = 0
      drlen(:) = 0
      do dproc = 1,neighnum
        iproc = neighlist(dproc)
        allocate( dpoints(dproc)%a(bnds(iproc)%len,4) )
        allocate( dbuf(dproc)%a(bnds(iproc)%len,4) )
        allocate( dbuf(dproc)%b(nagg*bnds(iproc)%len) )
        allocate( dindex(dproc)%a(bnds(iproc)%len,2) )
        allocate( sextra(dproc)%a(nagg*bnds(iproc)%len) )
      end do
      ! store invalid points in dproc=0
      allocate( dbuf(0)%a(4,1) )
      allocate( dbuf(0)%b(1) )
      allocate( dindex(0)%a(1,2) )


      ! Pack colour indices
      if ( myid==0 ) then
         write(6,*) "-> Define colour indices"
      end if
      allocate( ifull_colour(maxcolour) )
      ifull_colour(:) = 0
      do n = 1,npan
         do j = 1,jpan
            do i = 1,ipan
               iq  = indp(i,j,n)  ! Local
               iqg = indg(i,j,n)  ! Global
               colourmask = findcolour(iqg,il_g)
               ifull_colour(colourmask) = ifull_colour(colourmask) + 1 
            end do
         end do
      end do
      ifull_maxcolour = maxval( ifull_colour(:) )

      ! order points to allow border only updating
      allocate( ifull_colour_border(maxcolour), iqx(ifull_maxcolour,maxcolour) )
      ifull_colour(:) = 0
      ! first process border
      do n = 1,npan
         j = 1
         do i = 1,ipan
            iq  = indp(i,j,n)  ! Local
            iqg = indg(i,j,n)  ! Global
            colourmask = findcolour(iqg,il_g)
            ifull_colour(colourmask) = ifull_colour(colourmask) + 1
            iqx(ifull_colour(colourmask),colourmask) = iq
         end do
         do j = 2,jpan-1
            i = 1  
            iq  = indp(i,j,n)  ! Local
            iqg = indg(i,j,n)  ! Global
            colourmask = findcolour(iqg,il_g)
            ifull_colour(colourmask) = ifull_colour(colourmask) + 1
            iqx(ifull_colour(colourmask),colourmask) = iq
            i = ipan
            iq  = indp(i,j,n)  ! Local
            iqg = indg(i,j,n)  ! Global
            colourmask = findcolour(iqg,il_g)
            ifull_colour(colourmask) = ifull_colour(colourmask) + 1
            iqx(ifull_colour(colourmask),colourmask) = iq
         end do
         j = jpan
         do i = 1,ipan
            iq  = indp(i,j,n)  ! Local
            iqg = indg(i,j,n)  ! Global
            colourmask = findcolour(iqg,il_g)
            ifull_colour(colourmask) = ifull_colour(colourmask) + 1
            iqx(ifull_colour(colourmask),colourmask) = iq
         end do
      end do
      ifull_colour_border(1:maxcolour) = ifull_colour(1:maxcolour)
      ! next process interior
      do n = 1,npan
         do j = 2,jpan-1
            do i = 2,ipan-1
               iq  = indp(i,j,n)  ! Local
               iqg = indg(i,j,n)  ! Global
               colourmask = findcolour(iqg,il_g)
               ifull_colour(colourmask) = ifull_colour(colourmask) + 1
               iqx(ifull_colour(colourmask),colourmask) = iq
            end do
         end do
      end do


      ! Create MPI_SUMDR for calculating global sums with high precision
      if ( myid==0 ) then
         write(6,*) "-> Create MPI_SUMDR"
      end if
      ltrue = .true. 
      ! Operator MPI_SUMDR is created based on an external function DRPDR.
      call MPI_OP_CREATE( DRPDR, ltrue, MPI_SUMDR, ierr )
      
      
      ! prepare comm groups - used by scale-selective filter
      if ( myid==0 ) then
         write(6,*) "-> Create comm groups"
      end if
      npta = max( 6/nproc, 1 )     ! number of panels per process
      mproc = max( nproc/6, 1 )    ! number of processes per panel
      pprocn = myid*npta/mproc     ! start panel
      pprocx = pprocn + npta - 1   ! end panel
      hproc = pprocn*mproc/npta    ! host process for panel

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
   
   subroutine ccmpi_init
      integer(kind=4) :: lerr, lproc, lid
      integer(kind=8) :: begin_time, end_time, count_rate, count_max
#ifdef _OPENMP
      integer(kind=4) :: lprovided
#endif

      call system_clock( begin_time, count_rate, count_max )

      ! Global communicator
      call MPI_Init(lerr)
      call MPI_Comm_size(MPI_COMM_WORLD, lproc, lerr) ! Find number of processes
      call MPI_Comm_rank(MPI_COMM_WORLD, lid, lerr)   ! Find local process id
      nproc      = lproc
      myid       = lid
      comm_world = MPI_COMM_WORLD
      
      call system_clock( end_time, count_rate, count_max )
      mpiinit_time = real( end_time - begin_time, 8)/real(count_rate,8)

   end subroutine ccmpi_init
   
   subroutine ccmpi_reinit(newnproc)
      ! redefine comm_world to drop unused processes
      integer, intent(in) :: newnproc
      integer(kind=4) :: lerr, lid
      integer(kind=4) :: lcommout, lcommin
      integer(kind=4) :: colour
#ifdef usempi3
      integer(kind=4) :: lproc
#endif
      
      if ( newnproc < nproc ) then
         if ( myid == 0 ) then
            write(6,*) "Reducing number of processes from ",nproc," to ",newnproc  
         end if
         if ( myid < newnproc ) then
            colour = 1
         else
            colour = MPI_UNDEFINED
         end if
         lcommin = comm_world
         lid = myid
         call MPI_Comm_Split(lcommin, colour, lid, lcommout, lerr) ! redefine comm_world
         comm_world = lcommout
         nproc = newnproc
         if ( myid >= nproc ) then
            call ccmpi_finalize
            stop
         end if    
      end if
      
#ifdef usempi3
      ! Intra-node communicator 
      lid = myid
      lcommin = comm_world
      call MPI_Comm_split_type(lcommin, MPI_COMM_TYPE_SHARED, lid, MPI_INFO_NULL, lcommout, lerr)
      call MPI_Comm_size(lcommout, lproc, lerr) ! Find number of processes on node
      call MPI_Comm_rank(lcommout, lid, lerr)   ! Find local process id on node
      comm_node  = lcommout
      node_nproc = lproc
      node_myid  = lid

      ! Inter-node commuicator
      colour = node_myid
      lid = myid
      lcommin = comm_world
      call MPI_Comm_Split(lcommin, colour, lid, lcommout, lerr)
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
      
      if ( myid==0 .and. (node_myid/=0.or.nodecaptain_myid/=0) ) then
         write(6,*) "ERROR: Intra-node communicator failed"
         write(6,*) "myid, node_myid, nodecaptain_myid ",myid,node_myid,nodecaptain_myid
         call ccmpi_abort(-1)
      end if
#else
      ! each process is treated as a node
      lcommin = comm_world
      colour = myid
      lid = 0
      call MPI_Comm_Split(lcommin, colour, lid, lcommout, lerr)
      comm_node  = lcommout
      node_nproc = 1
      node_myid  = 0
      
      comm_nodecaptain = comm_world
      nodecaptain_nproc = nproc
      nodecaptain_myid = myid
      
      node_captainid = myid
#endif

   end subroutine ccmpi_reinit
   
   subroutine ccmpi_remap
      ! redefine comm_world for optimal internode communicsation
      integer :: node_dx, node_dy
      integer :: oldrank, ty, cy, tx, cx
      integer :: new_node_nproc
      integer :: vsize
      integer(kind=4) :: ref_nodecaptain_nproc
      integer(kind=4) :: lerr, lid, lcommin, lcommout
      integer(kind=4), dimension(:), allocatable :: nodesize
      integer(kind=4), dimension(1) :: lnode

      node_nx = 1
      node_ny = 1
      node_dx = nxp
      node_dy = nyp
      
#ifdef usempi3
      new_node_nproc = -1

      ! communicate number of nodes to all processes
      ref_nodecaptain_nproc = nodecaptain_nproc
      lcommin = comm_world      
      call MPI_Bcast(ref_nodecaptain_nproc, 1_4, MPI_INTEGER, 0_4, lcommin, lerr )
   
      if ( mod(nproc, ref_nodecaptain_nproc)==0 .and. ref_nodecaptain_nproc>1 ) then
         ! fully loaded nodes method
         new_node_nproc = nproc / ref_nodecaptain_nproc
      else
         ! virtual nodes method
         allocate( nodesize(ref_nodecaptain_nproc) )
         if ( node_myid==0 ) then
            lnode(1) = node_nproc
            lcommin = comm_nodecaptain
            call START_LOG(gather_begin)
            call MPI_Gather(lnode,1_4,MPI_INTEGER,nodesize,1_4,MPI_INTEGER,0_4,lcommin,lerr)
            call END_LOG(gather_end)
         end if
         lcommin = comm_world
         call MPI_Bcast(nodesize, ref_nodecaptain_nproc, MPI_INTEGER, 0_4, lcommin, lerr )
         do vsize = maxval(nodesize),1,-1 
            if ( all( mod(nodesize(:),vsize)==0 ) ) then
               exit 
            end if
         end do
         deallocate( nodesize )
         if ( vsize > 2 ) then
            new_node_nproc = vsize
         end if
      end if

      if ( new_node_nproc >= 4 ) then

         ! calculate virtual node decomposition
         node_nx = max( int(sqrt(real(new_node_nproc))), 1 )
         node_ny = new_node_nproc / node_nx
         node_dx = nxp / node_nx
         node_dy = nyp / node_ny
         do while ( (node_nx*node_dx/=nxp.or.node_ny*node_dy/=nyp.or.node_nx*node_ny/=new_node_nproc) .and. node_nx>0 )
            node_nx = node_nx - 1
            node_ny = new_node_nproc / max( node_nx, 1 )
            node_dx = nxp / max( node_nx, 1 )
            node_dy = nyp / node_ny
         end do   

        ! remap ranks if a valid decomposition has been found
         if ( myid == 0 ) then
            write(6,*) "Remapping ranks using node_nx,node_ny ",node_nx,node_ny
            write(6,*) "with node_dx,node_dy                  ",node_dx,node_dy
         end if
         !myid = ty*node_ny*node_dx*node_nx + tx*node_dx*node_nx + cy*node_nx + cx
         ty = myid/(node_dx*node_nx*node_ny)
         tx = (myid-ty*node_dx*node_nx*node_ny)/(node_nx*node_ny)
         cy = (myid-ty*node_dx*node_nx*node_ny-tx*node_nx*node_ny)/node_nx
         cx = myid-ty*node_dx*node_nx*node_ny-tx*node_nx*node_ny-cy*node_nx
         lid = ty*node_ny*node_dx*node_nx + tx*node_nx + cy*node_dx*node_nx + cx
         lcommin = comm_world
         call MPI_Comm_Split(lcommin, 0_4, lid, lcommout, lerr) ! redefine comm_world
         ! The next line is not needed since lid is already the correct rank
         call MPI_Comm_rank(lcommout, lid, lerr)                ! find local process id
         comm_world = lcommout
         myid = lid
         if ( lcommin/=MPI_COMM_WORLD .and. lcommin/=MPI_COMM_NULL ) then
            call MPI_Comm_Free(lcommin, lerr) 
         end if

      end if
#endif   
   
   end subroutine ccmpi_remap
   
   subroutine ccmpi_finalize   
      integer(kind=4) :: lerr
      
      call MPI_Finalize( lerr )
   
   end subroutine ccmpi_finalize
   
   subroutine ccmpi_procformat_init(localhist,procmode)
      ! define commuication handle to support procformat parallel output
      integer, intent(inout) :: procmode
      logical, intent(in) :: localhist
#ifdef usempi3
      integer(kind=4) :: colour, lcomm, lrank
      integer(kind=4) :: lcommout, lerr, lsize
#endif

      if ( localhist ) then
#ifdef usempi3
         ! configure procmode
         if ( procmode == 0 ) then
            ! procmode=0 uses existing nodes, even if they have different numbers of processes
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
            if ( myid == 0 ) then
               write(6,*) "Configure procformat output with procmode=",procmode
            end if
         else
            ! user specified procmode>0 
            procmode = max(procmode, 1)
            if ( procmode >= node_nx ) then
               do while ( mod(node_nproc,procmode)/=0 .and. mod(procmode,node_nx)/=0 .and. procmode>=node_nx )
                  procmode = procmode - 1 ! can be different on different nodes
               end do  
            end if    
            if ( procmode < node_nx ) then
               do while ( mod(node_nproc,procmode)/=0 .and. mod(node_nx,procmode)/=0 .and. procmode>0 )
                  procmode = procmode - 1 ! can be different on different nodes
               end do  
            end if
            if ( myid == 0 ) then
               write(6,*) "Configure procformat output with procmode=",procmode
            end if
            ! Intra-procmode communicator
            colour = node_myid/procmode
            lcomm = comm_node
            lrank = node_myid
            call MPI_Comm_Split(lcomm, colour, lrank, lcommout, lerr)
            comm_vnode = lcommout
            call MPI_Comm_Size(lcommout, lsize, lerr)
            vnode_nproc = lsize ! should equal procmode
            call MPI_Comm_Rank(lcommout, lrank, lerr)
            vnode_myid = lrank
            ! Inter-procmode communicator
            colour = vnode_myid
            lcomm = comm_world
            lrank = myid
            call MPI_Comm_Split(lcomm, colour, lrank, lcommout, lerr)
            comm_vleader = lcommout
            call MPI_Comm_Size(lcommout, lsize, lerr)
            vleader_nproc = lsize
            call MPI_Comm_Rank(lcommout, lrank, lerr)
            vleader_myid = lrank
            ! Communicate procmode id
            lcomm = comm_vnode
            lrank = vleader_myid
            call MPI_Bcast( lrank, 1_4, MPI_INTEGER, 0_4, lcomm, lerr )
            vnode_vleaderid = lrank
         end if
         ! just an internal check to make sure there are no errors
         if ( myid==0 .and. vnode_myid/=0 ) then
            write(6,*) "ERROR: vnode_myid/=0 with myid==0"
            call ccmpi_abort(-1)
         end if   
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

end module cc_mpi

