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

! Common subroutines and variables for MPI communication
    
module cc_mpi_common

#ifdef usempimod
#ifdef usempimod_f08
   use mpi_f08
#else   
   use mpi
#endif   
#endif
   use newmpar_m

   implicit none

   public

#ifndef usempimod
   include 'mpif.h'
#endif
   
   integer, save :: comm_world                                     ! global communication group
   integer, save :: myid                                           ! process rank for comm_world
   integer, save :: ipan, jpan                                     ! grid size on process
   integer, save :: ioff, joff, noff                               ! offset of process grid relative to global grid
   integer, save :: nxproc, nyproc                                 ! number of processes in the x and y directions
                                                                   ! for the global grid
   integer, save :: node_nx, node_ny                               ! number of processes in the x and y directions
                                                                   ! on a node
   integer, save :: nagg = 8                                       ! maximum number of levels to aggregate
   integer, save :: maxcolour = 2                                  ! maximum number of colours for iterative solvers
   
   integer, save :: maxbuflen, maxvertlen                          ! bounds buffer size   
   logical, save :: mydiag                                         ! true if diagnostic point id, jd is in my region
   
   integer, save :: comm_node, node_myid, node_nproc               ! node communicator
   integer, save :: comm_nodecaptain, nodecaptain_myid, &
                    nodecaptain_nproc                              ! node captain communicator
   integer, save :: node_captainid                                 ! rank of the node captain in the comm_nodecaptain group
   
   integer, save :: comm_proc, comm_rows, comm_cols                ! comm groups for scale-selective filter
   integer, save :: hproc, mproc, npta, pprocn, pprocx             ! decomposition parameters for scale-selective filter

   integer, save :: comm_vnode, vnode_nproc, vnode_myid            ! procformat communicator for node
   integer, save :: comm_vleader, vleader_nproc, vleader_myid      ! procformat communicator for node captain group   
   integer, save :: vnode_vleaderid                                ! rank of the procformat node captain

   integer(kind=4), save :: nreq, rreq                             ! number of messages requested and to be received
   integer(kind=4), allocatable, dimension(:), save :: ireq        ! requested message index
   integer, allocatable, dimension(:), save :: rlist               ! map of process index from requested message index
   
   interface ccmpi_reduce
      module procedure ccmpi_reduce1r, ccmpi_reduce2r, ccmpi_reduce3r
      module procedure ccmpi_reduce1c, ccmpi_reduce2c
   end interface
   interface ccmpi_reducer8
     module procedure ccmpi_reduce1rr8
   end interface
   interface ccmpi_allreduce
      module procedure ccmpi_allreduce1r, ccmpi_allreduce2r
      module procedure ccmpi_allreduce1c, ccmpi_allreduce2c
      module procedure ccmpi_allreduce1i, ccmpi_allreduce2i
      module procedure ccmpi_allreduce1l, ccmpi_allreduce2l
   end interface
   interface ccmpi_bcast
      module procedure ccmpi_bcast1i, ccmpi_bcast2i, ccmpi_bcast3i
      module procedure ccmpi_bcast1r, ccmpi_bcast2r, ccmpi_bcast3r, ccmpi_bcast4r
      module procedure ccmpi_bcast1s, ccmpi_bcast2s
   end interface
   interface ccmpi_bcastr8
      module procedure ccmpi_bcast1r8, ccmpi_bcast2r8, ccmpi_bcast3r8
   end interface
   interface ccmpi_gatherx
      module procedure ccmpi_gatherx23r, ccmpi_gatherx34r, ccmpi_gatherx45r
      module procedure ccmpi_gatherx12i, ccmpi_gatherx23i
   end interface
   interface ccmpi_gatherxr8
     module procedure ccmpi_gatherx23rr8, ccmpi_gatherx34rr8, ccmpi_gatherx45rr8
   end interface
   interface ccmpi_allgatherx
      module procedure ccmpi_allgatherx32r
   end interface
   interface ccmpi_scatterx
      module procedure ccmpi_scatterx32r
   end interface
   interface ccmpi_alltoall
      module procedure ccmpi_alltoall2l
   end interface
#ifdef usempi3
   interface ccmpi_allocshdata
      module procedure ccmpi_allocshdata2r, ccmpi_allocshdata3r, ccmpi_allocshdata4r, ccmpi_allocshdata5r, ccmpi_allocshdata6r 
   end interface
   interface ccmpi_allocshdatar8
      module procedure ccmpi_allocshdata2_r8, ccmpi_allocshdata3_r8
   end interface
#endif
   
   ! Do directions need to be swapped
   logical, parameter, dimension(0:npanels) ::                      &
           swap_e = (/ .true., .false., .true., .false., .true., .false. /), &
           swap_w = (/ .false., .true., .false., .true., .false., .true. /), &
           swap_n = (/ .false., .true., .false., .true., .false., .true. /), &
           swap_s = (/ .true., .false., .true., .false., .true., .false. /)

   ! partition indices into colours
   integer, dimension(:,:), allocatable, save :: iqx
   integer, save :: ifull_maxcolour
   integer, dimension(:), allocatable, save :: ifull_colour, ifull_colour_border

   ! flag whether process region edge is a face edge.
   logical, save :: edge_w, edge_n, edge_s, edge_e
   
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
      ! Number of points for each process. Also double row versions.
      ! lenx is first row plux corner points.  lenh is just the ne side.
      integer :: slen2, rlen2
      integer, dimension(:), allocatable :: rlenh_bg, rlenh_fn, slenh_bg, slenh_fn
      integer, dimension(:), allocatable :: rlen_bg, rlen_fn, slen_bg, slen_fn
      integer, dimension(:), allocatable :: rlenx_bg, rlenx_fn, slenx_bg, slenx_fn
      ! Vector groups
      integer :: rlen_su_bg, rlen_ev_fn
      integer :: rlen_sv_bg, rlen_wu_fn, rlen_nv_bg, rlen_eu_fn
      integer :: rlen_ssv_bg, rlen_wwu_fn, rlen_nnv_bg, rlen_eeu_fn
      integer :: slen_su_bg, slen_ev_fn
      integer :: slen_sv_bg, slen_wu_fn, slen_nv_bg, slen_eu_fn
      integer :: slen_ssv_bg, slen_wwu_fn, slen_nnv_bg, slen_eeu_fn
   end type bounds_info

   ! bounds data
   type(bounds_info), allocatable, dimension(:), save :: bnds
   
   ! Timer
   integer, save :: ints_begin, ints_end
   integer, save :: nonlin_begin, nonlin_end
   integer, save :: helm_begin, helm_end
   integer, save :: adjust_begin, adjust_end
   integer, save :: upglobal_begin, upglobal_end
   integer, save :: hordifg_begin, hordifg_end
   integer, save :: vadv_begin, vadv_end
   integer, save :: depts_begin, depts_end
   integer, save :: stag_begin, stag_end
   integer, save :: ocnstag_begin, ocnstag_end
   integer, save :: mfix_begin, mfix_end
   integer, save :: phys_begin, phys_end
   integer, save :: outfile_begin, outfile_end
   integer, save :: onthefly_begin, onthefly_end
   integer, save :: otf_fill_begin, otf_fill_end
   integer, save :: otf_ints_begin, otf_ints_end
   integer, save :: histrd_begin, histrd_end
   integer, save :: indata_begin, indata_end
   integer, save :: nestin_begin, nestin_end
   integer, save :: nestOTF_begin, nestOTF_end
   integer, save :: nestWIN_begin, nestWIN_end
   integer, save :: nestcalc_begin, nestcalc_end
   integer, save :: nestcomm_begin, nestcomm_end
   integer, save :: ensemble_begin, ensemble_end
   integer, save :: amipsst_begin, amipsst_end
   integer, save :: gwdrag_begin, gwdrag_end
   integer, save :: convection_begin, convection_end
   integer, save :: cloud_begin, cloud_end
   integer, save :: radnet_begin, radnet_end
   integer, save :: radinit_begin, radinit_end
   integer, save :: radSW_begin, radSW_end
   integer, save :: radLW_begin, radLW_end
   integer, save :: sfluxnet_begin, sfluxnet_end
   integer, save :: sfluxwater_begin, sfluxwater_end
   integer, save :: sfluxland_begin, sfluxland_end
   integer, save :: sfluxurban_begin, sfluxurban_end
   integer, save :: vertmix_begin, vertmix_end
   integer, save :: aerosol_begin, aerosol_end
   integer, save :: cape_begin, cape_end
   integer, save :: maincalc_begin, maincalc_end
   integer, save :: precon_begin, precon_end
   integer, save :: waterdynamics_begin, waterdynamics_end
   integer, save :: watermfix_begin, watermfix_end
   integer, save :: waterdeps_begin, waterdeps_end
   integer, save :: watereos_begin, watereos_end
   integer, save :: waterints_begin, waterints_end
   integer, save :: watervadv_begin, watervadv_end
   integer, save :: waterhelm_begin, waterhelm_end
   integer, save :: wateriadv_begin, wateriadv_end
   integer, save :: waterdiff_begin, waterdiff_end
   integer, save :: river_begin, river_end
   integer, save :: bcast_begin, bcast_end
   integer, save :: alltoall_begin, alltoall_end
   integer, save :: allgather_begin, allgather_end
   integer, save :: gather_begin, gather_end
   integer, save :: scatter_begin, scatter_end
   integer, save :: reduce_begin, reduce_end
   integer, save :: allreduce_begin, allreduce_end
   integer, save :: mpiwait_begin, mpiwait_end
   integer, save :: mpibarrier_begin, mpibarrier_end
   integer, save :: mgfine_begin, mgfine_end
   integer, save :: mgup_begin, mgup_end
   integer, save :: mgcoarse_begin, mgcoarse_end
   integer, save :: mgdown_begin, mgdown_end
   integer, save :: p1_begin, p1_end
   integer, save :: p2_begin, p2_end
   integer, save :: p3_begin, p3_end
   integer, save :: p4_begin, p4_end
   integer, save :: p5_begin, p5_end
   integer, save :: p6_begin, p6_end
   integer, save :: p7_begin, p7_end
   integer, save :: p8_begin, p8_end
   integer, save :: p9_begin, p9_end
   integer, save :: p10_begin, p10_end
   integer, save :: p11_begin, p11_end
   integer, save :: p12_begin, p12_end
   integer, save :: p13_begin, p13_end
   integer, save :: p14_begin, p14_end
   integer, save :: p15_begin, p15_end
   integer, save :: p16_begin, p16_end
   integer, save :: p17_begin, p17_end
   integer, parameter :: nevents = 83
   real(kind=8), dimension(nevents), save, private :: tot_time = 0._8, start_time
   real(kind=8), save :: mpiinit_time, total_time
   character(len=15), dimension(nevents), save, private :: event_name
   
#ifdef vampir
#include "vt_user.inc"
#endif

!$acc declare create(ipan,jpan)

contains
   
   pure subroutine indv_mpi(iq, i, j, n)
      ! Calculate local i, j, n from global iq
      integer , intent(in) :: iq
      integer , intent(out) :: i, j, n

      ! Global i, j, n
      n = (iq - 1)/(il_g**2)
      j = 1 + (iq - n*il_g**2 - 1)/il_g
      i = iq - (j - 1)*il_g - n*il_g**2
      ! Reduced to values on my process
      j = j - joff
      i = i - ioff
      n = n + noff      
   end subroutine indv_mpi

   pure function indg(i,j,n) result(iqg)
      ! Calculate a 1D global index from the local process indices
      ! n in range 1..npan
      integer, intent(in) :: i, j, n
      integer :: iqg

      iqg = i+ioff + (j+joff-1)*il_g + (n-noff)*il_g*il_g

   end function indg

   pure function indp(i,j,n) result(iq)
      ! Calculate a 1D local index from the local process indices
      ! Note that face number runs from 1 here.
      integer, intent(in) :: i, j, n
      integer :: iq

      iq = i + (j-1)*ipan + (n-1)*ipan*jpan

   end function indp

  pure function iq2iqg(iq) result(iqg)
      ! Calculate global iqg from local iq
      integer, intent(in) :: iq
      integer :: iqg
      integer :: i, j, n

      n = 1 + (iq-1)/(ipan*jpan)  ! In range 1 .. npan
      j = 1 + ( iq - (n-1)*ipan*jpan - 1) / ipan
      i = iq - (j-1)*ipan - (n-1)*ipan*jpan
      iqg = i+ioff + (j+joff-1)*il_g + (n-noff)*il_g**2

  end function iq2iqg
  
  pure function iqg2iq(iqg) result(iq)
      ! Calculate local iq from global iqg
      integer, intent(in) :: iqg
      integer :: iq
      integer :: i, j, n

      n = (iqg - 1)/(il_g**2)
      j = 1 + (iqg - n*il_g**2 - 1)/il_g
      i = iqg - (j - 1)*il_g - n*il_g**2
      iq = i-ioff + (j-joff-1)*ipan + (n+noff-1)*ipan*jpan

   end function iqg2iq

   pure function fproc(i,j,n) result(fpout)
      ! locate process that owns a global grid point
      integer, intent(in) :: i, j, n
      integer :: fpout

      fpout = (i-1)/ipan + ((j-1)/jpan)*nxproc + n*nxproc*nyproc/npan
   
   end function fproc

   pure function qproc(iqg) result(qpout)
      ! locates process that owns a global grid point
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
      
      n = (iq_g - 1)/(mil_g**2)
      j = 1 + (iq_g - n*mil_g**2 - 1)/mil_g
      i = iq_g - (j - 1)*mil_g - n*mil_g**2
      iq = i-mioff + (j-mjoff-1)*mipan + (n + mnoff-1)*mipan*mjpan
   
   end function indx_indv
   
   pure function findcolour(iqg,mil_g) result(icol)
      ! find colour of a global grid box for grid size mil_g
      integer, intent(in) :: iqg, mil_g
      integer icol
      integer ig, jg, ng, tg
      integer jx

      icol = -1 ! for gfortran
      
      ! calculate global i,j,n
      tg = iqg - 1
      ng = tg/(mil_g**2)
      tg = tg - ng*mil_g**2
      jg = tg/il_g
      tg = tg - jg*mil_g
      ig = tg + 1
      ig = ig + 1
      jg = jg + 1
      
      if ( maxcolour==3 ) then
         ! three colour mask
         jx = mod( ig + jg + ng*mil_g, 2 )
         select case( ng + jx*(npanels+1) )
            case( 0, 1, 3, 4 )
               icol = 1
            case( 2, 5, 6, 9 )
               icol = 2
            case( 7, 8, 10, 11 )
               icol = 3
         end select

      else if ( maxcolour==2 ) then
         ! two colour mask
         icol = mod( ig + jg + ng*mil_g, 2 ) + 1
         
      else if ( maxcolour==1 ) then
         ! single colour 
         icol = 1

      end if
   
   end function findcolour

   subroutine proc_setup(id,jd,idjd)
!     Routine to set up offsets etc.
      integer :: nd, jdf, idjd_g
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
      ! Convert global indices to ones on this process region
      idjd_g = id + (jd-1)*il_g
      if ( mydiag ) then
         idjd = iqg2iq(idjd_g)
      else
         ! This should never be used so set a value that will give a bounds error
         idjd = huge(1)
      end if

   end subroutine proc_setup

   subroutine face_set(ipan_l, jpan_l, noff_l, ioff_l, joff_l, npan_l, il_gx, myid_l, nproc_l, nxproc_l, nyproc_l)
      integer, intent(in) :: myid_l, nproc_l, npan_l, il_gx
      integer, intent(out) :: ipan_l, jpan_l, noff_l, nxproc_l, nyproc_l
      integer, dimension(0:npanels), intent(out) :: ioff_l, joff_l 
      integer :: n

      !  Process allocation
      !  if  nproc_l <= npanels+1, then each gets a number of full panels
      if ( nproc_l<=npanels+1 ) then
         if ( modulo( npanels+1, nproc_l )/=0 ) then
            write(6,*) "Error, number of processes must divide number of panels"
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
            write(6,*) "Error, number of processes must be a multiple of number of panels"
            call ccmpi_abort(-1)
         end if
!        npan_l = 1
         n = nproc_l/(npanels+1)
         !  n is the number of processes on each face
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

         ! Still need to check that the process distribution is compatible
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

         ! Set offsets for this process
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
      !  Process allocation: each process gets a part of each panel
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

      ! Still need to check that the process distribution is compatible
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

      ! Set offsets for this process
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
      !  Process allocation: each process gets a part of each panel
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

      ! Still need to check that the process distribution is compatible
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

      ! Set offsets for this process
      do n = 0,npanels
         call proc_region_uniform(myid_l,n,ioff_l(n),joff_l(n),noff_l,nxproc_l,nyproc_l,ipan_l,jpan_l)
      end do

   end subroutine uniform_set

   pure subroutine proc_region_face(procid,ipoff,jpoff,npoff,nxproc_l,nyproc_l,ipan_l,jpan_l,npan_l)
      ! Calculate the offsets for a given process
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
         ! mtmp is the process index on this face, 0:(nxprox*nyproc-1)
         mtmp = procid - myface*nxproc_l*nyproc_l
         jpoff = (mtmp/nxproc_l)*jpan_l
         ipoff = modulo( mtmp, nxproc_l )*ipan_l
      end if
     
   end subroutine proc_region_face

   pure subroutine proc_region_dix(procid,ipoff,jpoff,npoff,nxproc_l,ipan_l,jpan_l)
      ! Calculate the offsets for a given process
      integer, intent(in) :: procid, nxproc_l, ipan_l, jpan_l
      integer, intent(out) :: ipoff, jpoff, npoff

      ! Original Dix uniform decomposition
      ! Set offsets for this process (same on all faces)
      npoff = 1
      jpoff = (procid/nxproc_l) * jpan_l
      ipoff = modulo(procid,nxproc_l)*ipan_l
     
   end subroutine proc_region_dix
   
   pure subroutine proc_region_uniform(procid,panid,ipoff,jpoff,npoff,nxproc_l,nyproc_l,ipan_l,jpan_l)
      ! Calculate the offsets for a given process
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
      integer(kind=8) :: begin_time, count_rate, count_max
#ifdef vampir
      VT_USER_START(event_name(event))
#endif
      call system_clock( begin_time, count_rate, count_max )
      start_time(event) = real(begin_time,8)/real(count_rate,8)
   end subroutine start_log

   subroutine end_log ( event )
      integer, intent(in) :: event
      integer(kind=8) :: end_time, count_rate, count_max
#ifdef vampir
      VT_USER_END(event_name(event))
#endif
      call system_clock( end_time, count_rate, count_max )
      tot_time(event) = tot_time(event) + ( real(end_time,8)/real(count_rate,8) - start_time(event) )
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
      ! define events for log
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
      call add_event(watereos_begin,      watereos_end,      "Water_EOS")
      call add_event(waterhelm_begin,     waterhelm_end,     "Water_Helm")
      call add_event(watervadv_begin,     watervadv_end,     "Water_Vadv")
      call add_event(waterdeps_begin,     waterdeps_end,     "Water_Deps")
      call add_event(waterints_begin,     waterints_end,     "Water_Ints")
      call add_event(wateriadv_begin,     wateriadv_end,     "Water_Iadv")
      call add_event(watermfix_begin,     watermfix_end,     "Water_Mfix")
      call add_event(ocnstag_begin,       ocnstag_end,       "Water_Stag")
      call add_event(river_begin,         river_end,         "River")
      call add_event(mgfine_begin,        mgfine_end,        "MG_Fine")
      call add_event(mgup_begin,          mgup_end,          "MG_Up")
      call add_event(mgcoarse_begin,      mgcoarse_end,      "MG_Coarse")
      call add_event(mgdown_begin,        mgdown_end,        "MG_Down")
      call add_event(outfile_begin,       outfile_end,       "Outfile")
      call add_event(onthefly_begin,      onthefly_end,      "Onthefly")
      call add_event(otf_fill_begin,      otf_fill_end,      "OTF_Fill")
      call add_event(otf_ints_begin,      otf_ints_end,      "OTF_Doints")
      call add_event(histrd_begin,        histrd_end,        "OTF_HistRd")
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
      call add_event(radsw_begin,         radsw_end,         "Rad_SW")
      call add_event(radlw_begin,         radlw_end,         "Rad_LW")
      call add_event(sfluxnet_begin,      sfluxnet_end,      "Sflux")
      call add_event(sfluxwater_begin,    sfluxwater_end,    "Sflux_water")
      call add_event(sfluxland_begin,     sfluxland_end,     "Sflux_land")
      call add_event(sfluxurban_begin,    sfluxurban_end,    "Sflux_urban")
      call add_event(vertmix_begin,       vertmix_end,       "Vertmix")
      call add_event(aerosol_begin,       aerosol_end,       "Aerosol")
      call add_event(cape_begin,          cape_end,          "CAPE")
      call add_event(precon_begin,        precon_end,        "Precon")
      call add_event(bcast_begin,         bcast_end,         "MPI_Bcast")
      call add_event(alltoall_begin,      alltoall_end,      "MPI_AlltoAll")
      call add_event(allgather_begin,     allgather_end,     "MPI_AllGather")
      call add_event(gather_begin,        gather_end,        "MPI_Gather")
      call add_event(scatter_begin,       scatter_end,       "MPI_Scatter")
      call add_event(allreduce_begin,     allreduce_end,     "MPI_AllReduce")
      call add_event(reduce_begin,        reduce_end,        "MPI_Reduce")
      call add_event(mpiwait_begin,       mpiwait_end,       "MPI_Wait")
      call add_event(mpibarrier_begin,    mpibarrier_end,    "MPI_Barrier")
      call add_event(p1_begin,            p1_end,            "Probe1")
      call add_event(p2_begin,            p2_end,            "Probe2")
      call add_event(p3_begin,            p3_end,            "Probe3")
      call add_event(p4_begin,            p4_end,            "Probe4")
      call add_event(p5_begin,            p5_end,            "Probe5")
      call add_event(p6_begin,            p6_end,            "Probe6")
      call add_event(p7_begin,            p7_end,            "Probe7")
      call add_event(p8_begin,            p8_end,            "Probe8")
      call add_event(p9_begin,            p9_end,            "Probe9")
      call add_event(p10_begin,           p10_end,           "Probe10")
      call add_event(p11_begin,           p11_end,           "Probe11")
      call add_event(p12_begin,           p12_end,           "Probe12")
      call add_event(p13_begin,           p13_end,           "Probe13")
      call add_event(p14_begin,           p14_end,           "Probe14")
      call add_event(p15_begin,           p15_end,           "Probe15")
      call add_event(p16_begin,           p16_end,           "Probe16")
      call add_event(p17_begin,           p17_end,           "Probe17")
      
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

   subroutine simple_timer_finalize
      ! Calculate the mean, min and max times for each case
      integer :: i
      integer(kind=4) :: ierr, llen, lcomm
      real(kind=8), dimension(nevents) :: emean, emax, emin
      real(kind=8), dimension(2) :: time_l, time_mean, time_max, time_min
      
      llen = nevents
      lcomm = comm_world
      emean = 0._8
      emax = 0._8
      emin = 0._8
      call MPI_Reduce(tot_time, emean, llen, MPI_DOUBLE_PRECISION, &
                      MPI_SUM, 0_4, lcomm, ierr )
      call MPI_Reduce(tot_time, emax, llen, MPI_DOUBLE_PRECISION,  &
                      MPI_MAX, 0_4, lcomm, ierr )
      call MPI_Reduce(tot_time, emin, llen, MPI_DOUBLE_PRECISION,  &
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
      call MPI_Reduce(time_l, time_mean, llen, MPI_DOUBLE_PRECISION, &
                      MPI_SUM, 0_4, lcomm, ierr )
      time_mean = time_mean/real(nproc)
      call MPI_Reduce(time_l, time_max, llen, MPI_DOUBLE_PRECISION,  &
                      MPI_MAX, 0_4, lcomm, ierr )
      call MPI_Reduce(time_l, time_min, llen, MPI_DOUBLE_PRECISION,  &
                      MPI_MIN, 0_4, lcomm, ierr )
      if ( myid == 0 ) then
         write(*,"(a,3f10.3)") "MPI_Initialise ",time_mean(1),time_min(1),time_max(1)
         write(*,"(a,3f10.3)") "Total_Time     ",time_mean(2),time_min(2),time_max(2)
      end if   
        
   end subroutine simple_timer_finalize

   subroutine ccmpi_reduce1r(ldat,gdat,op,host,comm)
      integer, intent(in) :: host, comm
      real, intent(in) :: ldat
      real, intent(out) :: gdat
      character(len=*), intent(in) :: op
      real, dimension(1,1) :: ldat_l
      real, dimension(1,1) :: gdat_l

      ldat_l(1,1) = ldat
      gdat_l(:,:) = 0.
      call ccmpi_reduce3r(ldat_l,gdat_l,op,host,comm)
      gdat = gdat_l(1,1)
      
   end subroutine ccmpi_reduce1r
   
   subroutine ccmpi_reduce2r(ldat,gdat,op,host,comm)
   
      integer, intent(in) :: host, comm
      real, dimension(:), intent(in) :: ldat
      real, dimension(:), intent(out) :: gdat
      character(len=*), intent(in) :: op
      real, dimension(size(ldat),1) :: ldat_l
      real, dimension(size(gdat),1) :: gdat_l
 
      ldat_l(:,1) = ldat(:)
      gdat_l(:,:) = 0.
      call ccmpi_reduce3r(ldat_l,gdat_l,op,host,comm)
      gdat(:) = gdat_l(:,1)
   
   end subroutine ccmpi_reduce2r

   subroutine ccmpi_reduce3r(ldat,gdat,op,host,comm)
   
      integer, intent(in) :: host, comm
      real, dimension(:,:), intent(in) :: ldat
      real, dimension(:,:), intent(out) :: gdat
      character(len=*), intent(in) :: op
      integer(kind=4) :: lcomm, lerr, lsize, lhost

      call START_LOG(reduce_begin)
      
      lhost = host
      lcomm = comm
      lsize = size(ldat)
      
      select case( op )
         case( "max" )
#ifdef i8r8
            call MPI_Reduce(ldat, gdat, lsize, MPI_DOUBLE_PRECISION, MPI_MAX, lhost, lcomm, lerr )
#else
            call MPI_Reduce(ldat, gdat, lsize, MPI_REAL, MPI_MAX, lhost, lcomm, lerr )
#endif

         case( "min" )
#ifdef i8r8
            call MPI_Reduce(ldat, gdat, lsize, MPI_DOUBLE_PRECISION, MPI_MIN, lhost, lcomm, lerr )
#else
            call MPI_Reduce(ldat, gdat, lsize, MPI_REAL, MPI_MIN, lhost, lcomm, lerr )
#endif 
         case( "sum" )
#ifdef i8r8
            call MPI_Reduce(ldat, gdat, lsize, MPI_DOUBLE_PRECISION, MPI_SUM, lhost, lcomm, lerr )
#else
            call MPI_Reduce(ldat, gdat, lsize, MPI_REAL, MPI_SUM, lhost, lcomm, lerr )
#endif 
         case( "maxloc" )
            lsize = lsize/2
#ifdef i8r8
            call MPI_Reduce(ldat, gdat, lsize, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, lhost, lcomm, lerr )
#else
            call MPI_Reduce(ldat, gdat, lsize, MPI_2REAL, MPI_MAXLOC, lhost, lcomm, lerr )
#endif 
         case( "minloc" )
            lsize = lsize/2
#ifdef i8r8
            call MPI_Reduce(ldat, gdat, lsize, MPI_2DOUBLE_PRECISION, MPI_MINLOC, lhost, lcomm, lerr )
#else
            call MPI_Reduce(ldat, gdat, lsize, MPI_2REAL, MPI_MINLOC, lhost, lcomm, lerr )
#endif 
         case default
            write(6,*) "ERROR: Unknown option for ccmpi_reduce ",op
            call ccmpi_abort(-1)
      end select
            
      call END_LOG(reduce_end)
   
   end subroutine ccmpi_reduce3r

   subroutine ccmpi_reduce1c(ldat,gdat,op,host,comm)
      integer, intent(in) :: host,comm
      complex, intent(in) :: ldat
      complex, intent(out) :: gdat
      character(len=*), intent(in) :: op
      complex, dimension(1) :: ldat_l
      complex, dimension(1) :: gdat_l

      ldat_l(1) = ldat
      gdat_l(1) = cmplx( 0., 0. )
      call ccmpi_reduce2c(ldat_l,gdat_l,op,host,comm)
      gdat = gdat_l(1)

   end subroutine ccmpi_reduce1c

   subroutine ccmpi_reduce2c(ldat,gdat,op,host,comm)
      use sumdd_m   
      integer, intent(in) :: host,comm
      integer(kind=4) :: lcomm, lerr, lsize, lhost
      complex, dimension(:), intent(in) :: ldat
      complex, dimension(:), intent(out) :: gdat
      character(len=*), intent(in) :: op

      call START_LOG(reduce_begin)      
      
      lhost = host
      lcomm = comm
      lsize = size(ldat)
      
      select case( op )
         case( "max" )
#ifdef i8r8
            call MPI_Reduce(ldat, gdat, lsize, MPI_DOUBLE_COMPLEX, MPI_MAX, lhost, lcomm, lerr )
#else
            call MPI_Reduce(ldat, gdat, lsize, MPI_COMPLEX, MPI_MAX, lhost, lcomm, lerr )
#endif
         case( "min" )
#ifdef i8r8
            call MPI_Reduce(ldat, gdat, lsize, MPI_DOUBLE_COMPLEX, MPI_MIN, lhost, lcomm, lerr )
#else
            call MPI_Reduce(ldat, gdat, lsize, MPI_COMPLEX, MPI_MIN, lhost, lcomm, lerr )
#endif
         case( "sum" )
#ifdef i8r8
            call MPI_Reduce(ldat, gdat, lsize, MPI_DOUBLE_COMPLEX, MPI_SUM, lhost, lcomm, lerr )
#else
            call MPI_Reduce(ldat, gdat, lsize, MPI_COMPLEX, MPI_SUM, lhost, lcomm, lerr )
#endif
         case( "sumdr" )
#ifdef i8r8
            call MPI_Reduce(ldat, gdat, lsize, MPI_DOUBLE_COMPLEX, MPI_SUMDR, lhost, lcomm, lerr )
#else
            call MPI_Reduce(ldat, gdat, lsize, MPI_COMPLEX, MPI_SUMDR, lhost, lcomm, lerr )
#endif
         case default
            write(6,*) "ERROR: Unknown option for ccmpi_reduce ",op
            call ccmpi_abort(-1)
      end select
      
      call END_LOG(reduce_end)
   
   end subroutine ccmpi_reduce2c

   subroutine ccmpi_reduce1rr8(ldat,gdat,op,host,comm)
   
      integer, intent(in) :: host, comm
      integer(kind=4) :: lcomm, lerr, lhost
      real(kind=8), intent(in) :: ldat
      real(kind=8), intent(out) :: gdat
      character(len=*), intent(in) :: op

      call START_LOG(reduce_begin)
      
      lhost = host
      lcomm = comm
            
      select case( op )
         case( "max" )
            call MPI_Reduce(ldat, gdat, 1_4, MPI_DOUBLE_PRECISION, MPI_MAX, lhost, lcomm, lerr )
         case( "min" )
            call MPI_Reduce(ldat, gdat, 1_4, MPI_DOUBLE_PRECISION, MPI_MIN, lhost, lcomm, lerr )
         case( "sum" )
            call MPI_Reduce(ldat, gdat, 1_4, MPI_DOUBLE_PRECISION, MPI_SUM, lhost, lcomm, lerr )
         case default
            write(6,*) "ERROR: Unknown option for ccmpi_reduce ",op
            call ccmpi_abort(-1)
      end select

      call END_LOG(reduce_end)
   
   end subroutine ccmpi_reduce1rr8
   
   subroutine ccmpi_allreduce1i(ldat,gdat,op,comm)
      integer, intent(in) :: comm
      integer, intent(in) :: ldat
      integer, intent(out) :: gdat
      character(len=*), intent(in) :: op
      integer, dimension(1) :: ldat_l
      integer, dimension(1) :: gdat_l
      
      ldat_l(1) = ldat
      gdat_l(:) = 0
      call ccmpi_allreduce2i(ldat_l,gdat_l,op,comm)
      gdat = gdat_l(1)
   
   end subroutine ccmpi_allreduce1i
   
   subroutine ccmpi_allreduce2i(ldat,gdat,op,comm)
      integer, intent(in) :: comm
      integer(kind=4) :: lcomm, lerr, lsize
      integer, dimension(:), intent(in) :: ldat
      integer, dimension(:), intent(out) :: gdat
      character(len=*), intent(in) :: op

      call START_LOG(allreduce_begin)
      
      lcomm = comm
      lsize = size(ldat)
      
      select case( op )
         case( "max" )
#ifdef i8r8
            call MPI_AllReduce(ldat, gdat, lsize, MPI_INTEGER8, MPI_MAX, lcomm, lerr )
#else
            call MPI_AllReduce(ldat, gdat, lsize, MPI_INTEGER, MPI_MAX, lcomm, lerr )
#endif 
         case( "min" )
#ifdef i8r8
            call MPI_AllReduce(ldat, gdat, lsize, MPI_INTEGER8, MPI_MIN, lcomm, lerr )
#else
            call MPI_AllReduce(ldat, gdat, lsize, MPI_INTEGER, MPI_MIN, lcomm, lerr )
#endif 
         case( "sum" )
#ifdef i8r8
            call MPI_AllReduce(ldat, gdat, lsize, MPI_INTEGER8, MPI_SUM, lcomm, lerr )
#else
            call MPI_AllReduce(ldat, gdat, lsize, MPI_INTEGER, MPI_SUM, lcomm, lerr )
#endif 
         case default
            write(6,*) "ERROR: Unknown option for ccmpi_allreduce ",op
            call ccmpi_abort(-1)
      end select
      
      call END_LOG(allreduce_end)
   
   end subroutine ccmpi_allreduce2i

   subroutine ccmpi_allreduce1r(ldat,gdat,op,comm)
      integer, intent(in) :: comm
      real, intent(in) :: ldat
      real, intent(out) :: gdat
      character(len=*), intent(in) :: op
      real, dimension(1) :: ldat_l
      real, dimension(1) :: gdat_l

      ldat_l(1) = ldat
      gdat_l(:) = 0.
      call ccmpi_allreduce2r(ldat_l,gdat_l,op,comm)
      gdat = gdat_l(1)
   
   end subroutine ccmpi_allreduce1r
   
   subroutine ccmpi_allreduce2r(ldat,gdat,op,comm)
      integer, intent(in) :: comm
      integer(kind=4) lcomm, lerr, lsize
      real, dimension(:), intent(in) :: ldat
      real, dimension(:), intent(out) :: gdat
      character(len=*), intent(in) :: op

      call START_LOG(allreduce_begin)
      
      lcomm = comm
      lsize = size(ldat)
      
      select case( op )
         case( "max" )
#ifdef i8r8
            call MPI_AllReduce(ldat, gdat, lsize, MPI_DOUBLE_PRECISION, MPI_MAX, lcomm, lerr )
#else
            call MPI_AllReduce(ldat, gdat, lsize, MPI_REAL, MPI_MAX, lcomm, lerr )
#endif
         case( "min" )
#ifdef i8r8
            call MPI_AllReduce(ldat, gdat, lsize, MPI_DOUBLE_PRECISION, MPI_MIN, lcomm, lerr )
#else
            call MPI_AllReduce(ldat, gdat, lsize, MPI_REAL, MPI_MIN, lcomm, lerr )
#endif
         case( "sum" )
#ifdef i8r8
            call MPI_AllReduce(ldat, gdat, lsize, MPI_DOUBLE_PRECISION, MPI_SUM, lcomm, lerr )
#else
            call MPI_AllReduce(ldat, gdat, lsize, MPI_REAL, MPI_SUM, lcomm, lerr )
#endif
         case default
            write(6,*) "ERROR: Unknown option for ccmpi_allreduce ",op
            call ccmpi_abort(-1)
      end select

      call END_LOG(allreduce_end)
   
   end subroutine ccmpi_allreduce2r
  
   subroutine ccmpi_allreduce1c(ldat,gdat,op,comm)
      integer, intent(in) :: comm
      integer(kind=4) :: lcomm, lerr
      complex, intent(in) :: ldat
      complex, intent(out) :: gdat
      character(len=*), intent(in) :: op
      complex, dimension(1) :: ldat_l
      complex, dimension(1) :: gdat_l
      
      ldat_l(1) = ldat
      gdat_l(:) = cmplx( 0., 0. )
      call ccmpi_allreduce2c(ldat_l,gdat_l,op,comm)
      gdat = gdat_l(1)
   
   end subroutine ccmpi_allreduce1c
   
   subroutine ccmpi_allreduce2c(ldat,gdat,op,comm)
      use sumdd_m
      integer, intent(in) :: comm
      integer(kind=4) :: lcomm, lerr, lsize
      complex, dimension(:), intent(in) :: ldat
      complex, dimension(:), intent(out) :: gdat
      character(len=*), intent(in) :: op

      call START_LOG(allreduce_begin)
      
      lcomm = comm
      lsize = size(ldat)
      
      select case( op )
         case( "max" )
#ifdef i8r8
            call MPI_AllReduce(ldat, gdat, lsize, MPI_DOUBLE_COMPLEX, MPI_MAX, lcomm, lerr )
#else
            call MPI_AllReduce(ldat, gdat, lsize, MPI_COMPLEX, MPI_MAX, lcomm, lerr )
#endif
         case( "min" )
#ifdef i8r8
            call MPI_AllReduce(ldat, gdat, lsize, MPI_DOUBLE_COMPLEX, MPI_MIN, lcomm, lerr )
#else
            call MPI_AllReduce(ldat, gdat, lsize, MPI_COMPLEX, MPI_MIN, lcomm, lerr )
#endif
         case( "sum" )
#ifdef i8r8
            call MPI_AllReduce(ldat, gdat, lsize, MPI_DOUBLE_COMPLEX, MPI_SUM, lcomm, lerr )
#else
            call MPI_AllReduce(ldat, gdat, lsize, MPI_COMPLEX, MPI_SUM, lcomm, lerr )
#endif
         case( "sumdr" )
#ifdef i8r8
            call MPI_AllReduce(ldat, gdat, lsize, MPI_DOUBLE_COMPLEX, MPI_SUMDR, lcomm, lerr )
#else
            call MPI_AllReduce(ldat, gdat, lsize, MPI_COMPLEX, MPI_SUMDR, lcomm, lerr )
#endif
         case default
            write(6,*) "ERROR: Unknown option for ccmpi_allreduce ",op
            call ccmpi_abort(-1)
      end select      
      
      call END_LOG(allreduce_end)
   
   end subroutine ccmpi_allreduce2c
   
   subroutine ccmpi_allreduce1l(ldat,gdat,op,comm)
      integer, intent(in) :: comm
      integer(kind=4) :: lop, lcomm, ierr
      logical, intent(in) :: ldat
      logical, intent(out) :: gdat
      character(len=*), intent(in) :: op
      logical, dimension(1) :: ldat_l
      logical, dimension(1) :: gdat_l

      ldat_l(1) = ldat
      gdat_l(1) = .false.
      call ccmpi_allreduce2l(ldat_l,gdat_l,op,comm)
      gdat = gdat_l(1)
      
   end subroutine ccmpi_allreduce1l
   
   subroutine ccmpi_allreduce2l(ldat_in,gdat_out,op,comm)
      integer, intent(in) :: comm
      integer(kind=4) :: lcomm, ierr, lsize
      logical, dimension(:), intent(in) :: ldat_in
      logical, dimension(:), intent(out) :: gdat_out
      logical(kind=4), dimension(size(ldat_in)) :: ldat
      logical(kind=4), dimension(size(gdat_out)) :: gdat
      character(len=*), intent(in) :: op

      call START_LOG(allreduce_begin)
      
      lcomm = comm
      lsize = size(ldat_in)
      ldat = ldat_in

      select case( op )
         case( "or" )
            call MPI_AllReduce(ldat, gdat, lsize, MPI_LOGICAL, MPI_LOR, lcomm, ierr )
         case( "and" )
            call MPI_AllReduce(ldat, gdat, lsize, MPI_LOGICAL, MPI_LAND, lcomm, ierr )
         case default
            write(6,*) "ERROR: Unknown option for ccmpi_reduce ",op
            call ccmpi_abort(-1)
      end select
      
      gdat_out = gdat

      call END_LOG(allreduce_end)
   
   end subroutine ccmpi_allreduce2l
   
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
      integer, intent(inout) :: ldat
      integer, dimension(1,1) :: ldat_l

      ldat_l(1,1) = ldat
      call ccmpi_bcast3i(ldat_l,host,comm)
      ldat = ldat_l(1,1)
      
   end subroutine ccmpi_bcast1i

   subroutine ccmpi_bcast2i(ldat,host,comm)
      integer, intent(in) :: host, comm
      integer, dimension(:), intent(inout) :: ldat
      integer, dimension(size(ldat),1) :: ldat_l

      ldat_l(:,1) = ldat(:)
      call ccmpi_bcast3i(ldat_l,host,comm)
      ldat(:) = ldat_l(:,1)
         
   end subroutine ccmpi_bcast2i

   subroutine ccmpi_bcast3i(ldat,host,comm)
      integer, intent(in) :: host, comm
      integer(kind=4) :: lcomm, lhost, lerr, lsize
      integer, dimension(:,:), intent(inout) :: ldat

      call START_LOG(bcast_begin)
      
      lhost = host
      lcomm = comm
      lsize = size(ldat)
#ifdef i8r8
      call MPI_Bcast(ldat,lsize,MPI_INTEGER8,lhost,lcomm,lerr)
#else
      call MPI_Bcast(ldat,lsize,MPI_INTEGER,lhost,lcomm,lerr)
#endif

      call END_LOG(bcast_end)
      
   end subroutine ccmpi_bcast3i

   subroutine ccmpi_bcast1r(ldat,host,comm)
      integer, intent(in) :: host, comm
      real, intent(inout) :: ldat
      real, dimension(1,1,1) :: ldat_l

      ldat_l(1,1,1) = ldat
      call ccmpi_bcast4r(ldat_l,host,comm)
      ldat = ldat_l(1,1,1)
   
   end subroutine ccmpi_bcast1r
   
   subroutine ccmpi_bcast2r(ldat,host,comm)   
      integer, intent(in) :: host, comm
      real, dimension(:), intent(inout) :: ldat
      real, dimension(size(ldat),1,1) :: ldat_l

      ldat_l(:,1,1) = ldat(:)
      call ccmpi_bcast4r(ldat_l,host,comm)
      ldat(:) = ldat_l(:,1,1)
   
   end subroutine ccmpi_bcast2r

   subroutine ccmpi_bcast3r(ldat,host,comm)   
      integer, intent(in) :: host, comm
      real, dimension(:,:), intent(inout) :: ldat
      real, dimension(size(ldat,1),size(ldat,2),1) :: ldat_l

      ldat_l(:,:,1) = ldat(:,:)
      call ccmpi_bcast4r(ldat_l,host,comm)
      ldat(:,:) = ldat_l(:,:,1)
   
   end subroutine ccmpi_bcast3r

   subroutine ccmpi_bcast4r(ldat,host,comm)
      integer, intent(in) :: host, comm
      integer(kind=4) :: lcomm, lhost, lerr, lsize
      real, dimension(:,:,:), intent(inout) :: ldat

      call START_LOG(bcast_begin)
      
      lhost = host
      lcomm = comm
      lsize = size(ldat)
#ifdef i8r8
      call MPI_Bcast(ldat, lsize, MPI_DOUBLE_PRECISION, lhost, lcomm, lerr)
#else
      call MPI_Bcast(ldat, lsize, MPI_REAL, lhost, lcomm, lerr)
#endif

      call END_LOG(bcast_end)
   
   end subroutine ccmpi_bcast4r

   subroutine ccmpi_bcast1s(ldat,host,comm)
      integer, intent(in) :: host, comm
      character(len=*), intent(inout) :: ldat
      character(len=len(ldat)), dimension(1) :: ldat_l

      ldat_l(1) = ldat
      call ccmpi_bcast2s(ldat_l,host,comm)
      ldat = ldat_l(1)
      
   end subroutine ccmpi_bcast1s

   subroutine ccmpi_bcast2s(ldat,host,comm)
      integer, intent(in) :: host, comm
      integer(kind=4) :: lcomm, lhost, lerr, lsize, llen, lnum
      character(len=*), dimension(:), intent(inout) :: ldat
      integer :: i, j
      integer(kind=1), dimension(:,:), allocatable :: dummy

      ! MJT notes - MS Windows MPI_CHARACTER seems broken

      call START_LOG(bcast_begin)

      lhost = host
      lcomm = comm
      llen = len(ldat(1))
      lnum = size(ldat)
      lsize = llen*lnum
      allocate( dummy(llen,lnum) )
      do j = 1,lnum
         do i = 1,llen
            dummy(i,j) = int(iachar(ldat(j)(i:i)),1)
         end do   
      end do
      call MPI_Bcast(dummy,lsize,MPI_BYTE,lhost,lcomm,lerr)
      do j = 1,lnum
         do i = 1,llen
            ldat(j)(i:i) = achar(dummy(i,j))
         end do   
      end do
      deallocate( dummy )

      call END_LOG(bcast_end)
      
   end subroutine ccmpi_bcast2s
   
   subroutine ccmpi_bcast1r8(ldat,host,comm)
      integer, intent(in) :: host, comm
      real(kind=8), intent(inout) :: ldat
      real(kind=8), dimension(1,1) :: ldat_l
      
      ldat_l(1,1) = ldat
      call ccmpi_bcast3r8(ldat_l,host,comm)
      ldat = ldat_l(1,1)      
   
   end subroutine ccmpi_bcast1r8
   
   subroutine ccmpi_bcast2r8(ldat,host,comm)
      integer, intent(in) :: host, comm
      real(kind=8), dimension(:), intent(inout) :: ldat
      real(kind=8), dimension(size(ldat),1) :: ldat_l
      
      ldat_l(:,1) = ldat(:)
      call ccmpi_bcast3r8(ldat_l,host,comm)
      ldat(:) = ldat_l(:,1)
   
   end subroutine ccmpi_bcast2r8
   
   subroutine ccmpi_bcast3r8(ldat,host,comm)
      integer, intent(in) :: host,comm
      integer(kind=4) :: lcomm, lhost, ierr, lsize
      real(kind=8), dimension(:,:), intent(inout) :: ldat
   
      lhost = host
      lcomm = comm
      lsize = size(ldat)
      call START_LOG(bcast_begin)
      call MPI_Bcast(ldat, lsize, MPI_DOUBLE_PRECISION, lhost, lcomm, ierr)
      call END_LOG(bcast_end)
   
   end subroutine ccmpi_bcast3r8

   subroutine ccmpi_barrier(comm)
      integer, intent(in) :: comm
      integer(kind=4) :: lcomm, ierr
      
      lcomm = comm
      call START_LOG(mpibarrier_begin)
      call MPI_Barrier( lcomm, ierr )
      call END_LOG(mpibarrier_end)
   
   end subroutine ccmpi_barrier
     
   subroutine ccmpi_gatherx23r(gdat,ldat,host,comm)
      integer, intent(in) :: host, comm
      real, dimension(:,:), intent(out) :: gdat
      real, dimension(:), intent(in) :: ldat
      real, dimension(size(gdat,1),size(gdat,2),1,1) :: gdat_l
      real, dimension(size(ldat,1),1,1) :: ldat_l
      
      ldat_l(:,1,1) = ldat(:)
      gdat_l(:,:,:,:) = 0.
      call ccmpi_gatherx45r(gdat_l,ldat_l,host,comm)
      gdat(:,:) = gdat_l(:,:,1,1)
      
   end subroutine ccmpi_gatherx23r
   
   subroutine ccmpi_gatherx34r(gdat,ldat,host,comm)
      integer, intent(in) :: host, comm
      real, dimension(:,:,:), intent(out) :: gdat
      real, dimension(:,:), intent(in) :: ldat
      real, dimension(size(gdat,1),size(gdat,2),size(gdat,3),1) :: gdat_l
      real, dimension(size(ldat,1),size(ldat,2),1) :: ldat_l
      
      ldat_l(:,:,1) = ldat(:,:)
      gdat_l(:,:,:,:) = 0.
      call ccmpi_gatherx45r(gdat_l,ldat_l,host,comm)
      gdat(:,:,:) = gdat_l(:,:,:,1)
      
   end subroutine ccmpi_gatherx34r

   subroutine ccmpi_gatherx45r(gdat,ldat,host,comm)
      integer, intent(in) :: host, comm
      integer(kind=4) :: lsize, lhost, lcomm, lerr
      real, dimension(:,:,:,:), intent(out) :: gdat
      real, dimension(:,:,:), intent(in) :: ldat

      call START_LOG(gather_begin)
      
      lcomm = comm
      lhost = host
      lsize = size(ldat)

#ifdef i8r8
      call MPI_Gather(ldat,lsize,MPI_DOUBLE_PRECISION,gdat,lsize,MPI_DOUBLE_PRECISION,lhost,lcomm,lerr)
#else
      call MPI_Gather(ldat,lsize,MPI_REAL,gdat,lsize,MPI_REAL,lhost,lcomm,lerr)
#endif

      call END_LOG(gather_end)
      
   end subroutine ccmpi_gatherx45r
   
   subroutine ccmpi_gatherx12i(gdat,ldat,host,comm)
      integer, intent(in) :: host, comm
      integer, dimension(:), intent(out) :: gdat
      integer, intent(in) :: ldat
      integer, dimension(size(gdat,1),1) :: gdat_l
      integer, dimension(1) :: ldat_l
      
      ldat_l(1) = ldat
      gdat_l(:,:) = 0
      call ccmpi_gatherx23i(gdat_l,ldat_l,host,comm)
      gdat(:) = gdat_l(:,1)
      
   end subroutine ccmpi_gatherx12i

   subroutine ccmpi_gatherx23i(gdat,ldat,host,comm)
      integer, intent(in) :: host, comm
      integer(kind=4) :: lsize, lhost, lcomm, lerr
      integer, dimension(:,:), intent(out) :: gdat
      integer, dimension(:), intent(in) :: ldat

      call START_LOG(gather_begin)
      
      lcomm = comm
      lhost = host
      lsize = size(ldat)

#ifdef i8r8
      call MPI_Gather(ldat,lsize,MPI_INTEGER8,gdat,lsize,MPI_INTEGER8,lhost,lcomm,lerr)
#else
      call MPI_Gather(ldat,lsize,MPI_INTEGER,gdat,lsize,MPI_INTEGER,lhost,lcomm,lerr)
#endif

      call END_LOG(gather_end)
      
   end subroutine ccmpi_gatherx23i
      
   subroutine ccmpi_gatherx23rr8(gdat,ldat,host,comm)
      integer, intent(in) :: host, comm
      integer(kind=4) :: lsize, lhost, lcomm, lerr
      real(kind=8), dimension(:,:), intent(out) :: gdat
      real(kind=8), dimension(:), intent(in) :: ldat
      real(kind=8), dimension(size(gdat,1),size(gdat,2),1,1) :: gdat_l
      real(kind=8), dimension(size(ldat,1),1,1) :: ldat_l
      
      ldat_l(:,1,1) = ldat(:)
      gdat_l(:,:,:,:) = 0._8
      call ccmpi_gatherx45rr8(gdat_l,ldat_l,host,comm)
      gdat(:,:) = gdat_l(:,:,1,1)
      
   end subroutine ccmpi_gatherx23rr8
   
   subroutine ccmpi_gatherx34rr8(gdat,ldat,host,comm)
      integer, intent(in) :: host, comm
      real(kind=8), dimension(:,:,:), intent(out) :: gdat
      real(kind=8), dimension(:,:), intent(in) :: ldat
      real(kind=8), dimension(size(gdat,1),size(gdat,2),size(gdat,3),1) :: gdat_l
      real(kind=8), dimension(size(ldat,1),size(ldat,2),1) :: ldat_l
      
      ldat_l(:,:,1) = ldat(:,:)
      gdat_l(:,:,:,:) = 0._8
      call ccmpi_gatherx45rr8(gdat_l,ldat_l,host,comm)
      gdat(:,:,:) = gdat_l(:,:,:,1)
      
   end subroutine ccmpi_gatherx34rr8

   subroutine ccmpi_gatherx45rr8(gdat,ldat,host,comm)
      integer, intent(in) :: host, comm
      integer(kind=4) :: lsize, lhost, lcomm, lerr
      real(kind=8), dimension(:,:,:,:), intent(out) :: gdat
      real(kind=8), dimension(:,:,:), intent(in) :: ldat

      call START_LOG(gather_begin)

      lcomm = comm
      lhost = host
      lsize = size(ldat)

      call MPI_Gather(ldat,lsize,MPI_DOUBLE_PRECISION,gdat,lsize,MPI_DOUBLE_PRECISION,lhost,lcomm,lerr)

      call END_LOG(gather_end)
   
   end subroutine ccmpi_gatherx45rr8

   subroutine ccmpi_scatterx32r(gdat,ldat,host,comm)
      integer, intent(in) :: host, comm
      integer(kind=4) :: lsize, lhost, lcomm, lerr
      real, dimension(:,:), intent(in) :: gdat
      real, dimension(:), intent(out) :: ldat

      call START_LOG(scatter_begin)

      lcomm = comm
      lhost = host
      lsize = size(ldat)

#ifdef i8r8
      call MPI_Scatter(gdat,lsize,MPI_DOUBLE_PRECISION,ldat,lsize,MPI_DOUBLE_PRECISION,lhost,lcomm,lerr)
#else
      call MPI_Scatter(gdat,lsize,MPI_REAL,ldat,lsize,MPI_REAL,lhost,lcomm,lerr)
#endif

      call END_LOG(scatter_end)
      
   end subroutine ccmpi_scatterx32r
      
   subroutine ccmpi_allgatherx32r(gdat,ldat,comm)
      integer, intent(in) :: comm
      integer(kind=4) lsize, lcomm, lerr
      real, dimension(:,:), intent(in) :: ldat
      real, dimension(:), intent(out) :: gdat

      call START_LOG(allgather_begin)
      
      lcomm = comm
      lsize = size(ldat)

#ifdef i8r8
      call MPI_AllGather(ldat,lsize,MPI_DOUBLE_PRECISION,gdat,lsize,MPI_DOUBLE_PRECISION,lcomm,lerr)
#else
      call MPI_AllGather(ldat,lsize,MPI_REAL,gdat,lsize,MPI_REAL,lcomm,lerr)
#endif

      call END_LOG(allgather_end)
      
   end subroutine ccmpi_allgatherx32r
   
   subroutine ccmpi_alltoall2l(gdat,comm)
      integer, intent(in) :: comm
      integer(kind=4) :: lsize, lcomm, lerr
      logical, dimension(:), intent(inout) :: gdat
      logical(kind=4), dimension(size(gdat)) :: rdat, qdat
      
      lcomm = comm
      lsize = size(gdat)/nproc
      call START_LOG(alltoall_begin)
      qdat(:) = gdat(:)
      call MPI_AlltoAll( qdat, lsize, MPI_LOGICAL, rdat, lsize, MPI_LOGICAL, lcomm, lerr )
      gdat(:) = rdat(:)
      call END_LOG(alltoall_end)
      
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
#ifdef i8r8
      call MPI_Type_size( MPI_DOUBLE_PRECISION, tsize, ierr )
#else
      call MPI_Type_size( MPI_REAL, tsize, ierr )
#endif
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
#ifdef i8r8
      call MPI_Type_size( MPI_DOUBLE_PRECISION, tsize, ierr )
#else
      call MPI_Type_size( MPI_REAL, tsize, ierr )
#endif
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
#ifdef i8r8
      call MPI_Type_size( MPI_DOUBLE_PRECISION, tsize, ierr )
#else
      call MPI_Type_size( MPI_REAL, tsize, ierr )
#endif
      if ( lmyid == 0 ) then
         lsize = product(sshape)*tsize
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
#ifdef i8r8
      call MPI_Type_size( MPI_DOUBLE_PRECISION, tsize, ierr )
#else
      call MPI_Type_size( MPI_REAL, tsize, ierr )
#endif
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

   end subroutine ccmpi_allocshdata5r 

   subroutine ccmpi_allocshdata6r(pdata,sshape,win,comm_in,myid_in)
      use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer
      integer, intent(out) :: win
      integer, intent(in), optional :: comm_in, myid_in
      integer :: lmyid
      integer, dimension(5), intent(in) :: sshape
      integer(kind=MPI_ADDRESS_KIND) :: qsize, lsize
      integer(kind=4) :: disp_unit, ierr, tsize
      integer(kind=4) :: lcomm, lwin
      real, pointer, dimension(:,:,:,:,:) :: pdata 
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
#ifdef i8r8
      call MPI_Type_size( MPI_DOUBLE_PRECISION, tsize, ierr )
#else
      call MPI_Type_size( MPI_REAL, tsize, ierr )
#endif
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

   end subroutine ccmpi_allocshdata6r 
   
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
   
end module cc_mpi_common

