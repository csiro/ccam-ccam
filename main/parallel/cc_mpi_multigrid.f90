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

! Subroutines and variables for the geometric multigrid solver    
    
module cc_mpi_multigrid

   use cc_mpi_common

   implicit none

   private

   public :: mgbounds, mgcollect, mgbcast, mg_index, mg_fproc, mg_fproc_1,  &
             mgbounds_colour
   
   interface mgbounds
      module procedure mgbounds2, mgbounds3
   end interface
   interface mgbounds_colour
      module procedure mgbounds_colour3
   end interface
   interface mgbcast
      module procedure mgbcast2, mgbcast3
      module procedure mgbcasta2, mgbcasta3
      module procedure mgbcastxn3
   end interface
   interface mgcollect
      module procedure mgcollect1, mgcollectreduce, mgcollectxn
   end interface

   ! Multi-grid arrays
   type mgtype
      integer :: ifull, iextra, ifull_fine, ifull_coarse
      integer :: merge_len, merge_row, ipan, merge_pos, nmax
      integer :: comm_merge, neighnum, npanx
      integer :: ifull_maxcolour
      integer, dimension(:), allocatable :: merge_list
      integer, dimension(:), allocatable :: in, ie, is, iw
      integer, dimension(:), allocatable :: ine, ien, ise, ies, inw, iwn, isw, iws
      integer, dimension(:), allocatable :: coarse_a, coarse_b, coarse_c
      integer, dimension(:), allocatable :: fine, fine_n, fine_e, fine_ne
      integer, dimension(:), allocatable :: neighlist
      integer, dimension(:), allocatable :: procmap
      integer, dimension(:), allocatable :: ifull_colour
      integer, dimension(:,:), allocatable :: iqx
      real, dimension(:), allocatable :: zzn, zze, zzs, zzw, zz
   end type mgtype

   type mgbndtype
      integer :: len
      integer, dimension(:), allocatable :: rlen_bg, rlen_fn, slen_bg, slen_fn
      integer, dimension(:), allocatable :: rlenx_bg, rlenx_fn, slenx_bg, slenx_fn      
      integer, dimension(:), allocatable :: send_list
      integer, dimension(:), allocatable :: unpack_list
      integer, dimension(:), allocatable :: request_list
   end type mgbndtype

   ! Multi-grid levels and buffers
   type(mgtype), dimension(:), allocatable, save, public :: mg
   type(mgbndtype), dimension(:,:), allocatable, save, public :: mg_bnds
   integer, save, public :: mg_maxlevel, mg_maxlevel_local

contains
   
   subroutine mgcollectreduce(g,vdat,dsolmax,klim)
      ! This subroutine merges datasets when upscaling with the multi-grid solver
      integer, intent(in) :: g
      integer, intent(in), optional :: klim
      integer :: kx, msg_len, ipanx, jpanx
      real, dimension(:,:), intent(inout) :: vdat
      real, dimension(:), intent(inout) :: dsolmax

      ! merge length
      if ( mg(g)%merge_len <= 1 ) return

      kx = size(vdat, 2)
      
      if ( present(klim) ) then
         kx = klim
      end if

      msg_len = mg(g)%ifull/(mg(g)%merge_len*mg(g)%npanx) ! message unit size
      ipanx = mg(g)%ipan
      jpanx = mg(g)%ifull/(mg(g)%ipan*mg(g)%npanx)
      call mgcollectreduce_work( g, vdat, dsolmax, kx, mg(g)%nmax, msg_len, mg(g)%npanx, ipanx, jpanx )

   return
   end subroutine mgcollectreduce

   subroutine mgcollectreduce_work(g,vdat,dsolmax,kx,nmax,msg_len,npanx,ipanx,jpanx)
      integer, intent(in) :: g, kx, nmax, msg_len, npanx, ipanx, jpanx
      real, dimension(:,:), intent(inout) :: vdat
      real, dimension(:), intent(inout) :: dsolmax
      integer :: nrow, ncol
      integer :: yproc, ir, ic, is, js, k, n, j, iq, iqq
      integer :: nrm1, hoz_len
      integer(kind=4) :: ierr, ilen, lcomm
      real, dimension((msg_len*npanx+1)*kx) :: tdat
      real, dimension((msg_len*npanx+1)*kx,nmax) :: tdat_g

      ! prep data for sending around the merge
      nrow    = mg(g)%ipan/mg(g)%merge_row  ! number of points along a row per process
      ncol    = msg_len/nrow                ! number of points along a col per process
      nrm1    = nrow - 1
      hoz_len = msg_len*npanx
      ilen = (hoz_len+1)*kx

      ! pack contiguous buffer
      tdat(1:hoz_len*kx) = reshape( vdat(1:hoz_len,1:kx), (/ hoz_len*kx /) )
      tdat(hoz_len*kx+1:ilen) = dsolmax(1:kx)

      lcomm = mg(g)%comm_merge
      call START_LOG(gather_begin)
#ifdef i8r8
      call MPI_Gather( tdat, ilen, MPI_DOUBLE_PRECISION, tdat_g, ilen, MPI_DOUBLE_PRECISION, 0_4, lcomm, ierr )
#else
      call MPI_Gather( tdat, ilen, MPI_REAL, tdat_g, ilen, MPI_REAL, 0_4, lcomm, ierr )
#endif
      call END_LOG(gather_end)

      ! unpack buffers (nmax is zero unless this is the host process)
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

   end subroutine mgcollectreduce_work

   subroutine mgcollect1(g,vdat,klim)
      ! This routing collects data from other processes without a reduction (max or min) array
      integer, intent(in) :: g
      integer, intent(in), optional :: klim
      integer :: kx, msg_len, ipanx, jpanx
      real, dimension(:,:), intent(inout) :: vdat

      ! merge length
      if ( mg(g)%merge_len <= 1 ) return
      
      kx = size(vdat, 2)
      
      if (present(klim)) then
         kx = klim
      end if

      msg_len = mg(g)%ifull/(mg(g)%merge_len*mg(g)%npanx) ! message unit size
      ipanx = mg(g)%ipan
      jpanx = mg(g)%ifull/(mg(g)%ipan*mg(g)%npanx)
      call mgcollect_work( g, vdat, kx, mg(g)%nmax, msg_len, mg(g)%npanx, ipanx, jpanx )
  
   end subroutine mgcollect1

   subroutine mgcollect_work(g,vdat,kx,nmax,msg_len,npanx,ipanx,jpanx)
      integer, intent(in) :: g, kx, nmax, msg_len, npanx, ipanx, jpanx
      real, dimension(:,:), intent(inout) :: vdat
      integer :: nrow, ncol
      integer :: yproc, ir, ic, is, js, k, n, j, iq, iqq
      integer :: nrm1, hoz_len
      integer(kind=4) :: ierr, ilen, lcomm
      real, dimension(msg_len*npanx*kx) :: tdat
      real, dimension(msg_len*npanx*kx,nmax) :: tdat_g

      ! prep data for sending around the merge
      nrow    = mg(g)%ipan/mg(g)%merge_row       ! number of points along a row per process
      ncol    = msg_len/nrow                     ! number of points along a col per process
      nrm1    = nrow - 1
      hoz_len = msg_len*npanx
      ilen = hoz_len*kx

      ! pack contiguous buffer
      tdat(1:ilen) = reshape( vdat(1:hoz_len,1:kx), (/ hoz_len*kx /) )

      lcomm = mg(g)%comm_merge
      call START_LOG(gather_begin)
#ifdef i8r8
      call MPI_Gather( tdat, ilen, MPI_DOUBLE_PRECISION, tdat_g, ilen, MPI_DOUBLE_PRECISION, 0_4, lcomm, ierr )
#else
      call MPI_Gather( tdat, ilen, MPI_REAL, tdat_g, ilen, MPI_REAL, 0_4, lcomm, ierr )
#endif
      call END_LOG(gather_end)

      ! unpack buffers (nmax is zero unless this is the host process)
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

   end subroutine mgcollect_work

   subroutine mgcollectxn(g,vdat,smaxmin,klim)
      ! This version of mgcollect also performs a max and min reduction
      integer, intent(in) :: g
      integer, intent(in), optional :: klim
      integer kx, msg_len, ipanx, jpanx
      real, dimension(:,:), intent(inout) :: vdat
      real, dimension(:,:), intent(inout) :: smaxmin

      ! merge length
      if ( mg(g)%merge_len <= 1 ) return
      
      kx = size(vdat, 2)
      
      if (present(klim)) then
         kx = klim
      end if

      msg_len = mg(g)%ifull/(mg(g)%merge_len*mg(g)%npanx) ! message unit size
      ipanx = mg(g)%ipan
      jpanx = mg(g)%ifull/(mg(g)%ipan*mg(g)%npanx)
      call mgcollectxn_work( g, vdat, smaxmin, kx, mg(g)%nmax, msg_len, mg(g)%npanx, ipanx, jpanx )
  
   end subroutine mgcollectxn

   subroutine mgcollectxn_work(g,vdat,smaxmin,kx,nmax,msg_len,npanx,ipanx,jpanx)
      integer, intent(in) :: g, kx, nmax, msg_len, npanx, ipanx, jpanx
      integer :: nrow, ncol
      integer :: yproc, ir, ic, is, js, k, n, j, iq, iqq
      integer :: nrm1, hoz_len
      integer(kind=4) :: ierr, ilen, lcomm
      real, dimension(:,:), intent(inout) :: vdat
      real, dimension(:,:), intent(inout) :: smaxmin
      real, dimension((msg_len*npanx+2)*kx) :: tdat
      real, dimension((msg_len*npanx+2)*kx,nmax) :: tdat_g

      ! prep data for sending around the merge
      nrow    = mg(g)%ipan/mg(g)%merge_row  ! number of points along a row per process
      ncol    = msg_len/nrow                ! number of points along a col per process
      nrm1    = nrow - 1
      hoz_len = msg_len*npanx
      ilen = (hoz_len+2)*kx

      ! pack contiguous buffer
      tdat(1:hoz_len*kx) = reshape( vdat(1:hoz_len,1:kx), (/ hoz_len*kx /) )
      tdat(hoz_len*kx+1:ilen) = reshape( smaxmin(1:kx,1:2), (/ kx*2 /) )

      lcomm = mg(g)%comm_merge
      call START_LOG(gather_begin)
#ifdef i8r8
      call MPI_Gather( tdat, ilen, MPI_DOUBLE_PRECISION, tdat_g, ilen, MPI_DOUBLE_PRECISION, 0_4, lcomm, ierr )
#else
      call MPI_Gather( tdat, ilen, MPI_REAL, tdat_g, ilen, MPI_REAL, 0_4, lcomm, ierr )
#endif
      call END_LOG(gather_end)

      ! unpack buffers (nmax is zero unless this is the host process)
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
  
   end subroutine mgcollectxn_work

   subroutine mgbcast2(g,vdat,dsolmax,nobounds)
      integer, intent(in) :: g
      real, dimension(:), intent(inout) :: vdat
      real, dimension(:), intent(inout) :: dsolmax
      logical, intent(in), optional :: nobounds
      real, dimension(size(vdat),1) :: vdat_l
      logical :: nobounds_l
      
      nobounds_l = .false.
      if ( present(nobounds) ) then
         nobounds_l = nobounds
      end if
      
      vdat_l(:,1) = vdat(:)
      call mgbcast3( g, vdat_l, dsolmax, nobounds=nobounds_l )
      vdat(:) = vdat_l(:,1)

   end subroutine mgbcast2
   
   subroutine mgbcast3(g,vdat,dsolmax,klim,nobounds)
      ! This subroutine broadcasts data when downscaling with the multi-grid solver
      ! This version also updates the convergence (dsolmax)
      integer, intent(in) :: g
      integer, intent(in), optional :: klim
      integer :: kx, out_len, dx
      integer(kind=4) :: ierr, ilen, lcomm
      real, dimension(:,:), intent(inout) :: vdat
      real, dimension(:), intent(inout) :: dsolmax
      real, dimension((mg(g)%ifull+mg(g)%iextra+1)*size(vdat,2)) :: tdat
      logical, intent(in), optional :: nobounds
      logical :: nbflag
      
      if ( mg(g)%merge_len <= 1 ) return

      kx = size(vdat,2)
      dx = size(dsolmax)
      if ( present(klim) ) then
         kx = klim
         dx = min( dx, kx )
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
      ilen = out_len*kx + dx
      
      ! pack contiguous buffer
      tdat(1:out_len*kx) = reshape( vdat(1:out_len,1:kx), (/ out_len*kx /) )
      tdat(out_len*kx+1:ilen) = dsolmax(1:dx)

      lcomm = mg(g)%comm_merge
      call START_LOG(bcast_begin)
#ifdef i8r8
      call MPI_Bcast( tdat, ilen, MPI_DOUBLE_PRECISION, 0_4, lcomm, ierr )
#else
      call MPI_Bcast( tdat, ilen, MPI_REAL, 0_4, lcomm, ierr )
#endif
      call END_LOG(bcast_end)

      ! extract data from Bcast
      vdat(1:out_len,1:kx) = reshape( tdat(1:out_len*kx), (/ out_len, kx /) )
      dsolmax(1:dx) = tdat(out_len*kx+1:ilen)
   
   end subroutine mgbcast3

   subroutine mgbcasta2(g,vdat,nobounds)
      integer, intent(in) :: g
      real, dimension(:), intent(inout) :: vdat
      logical, intent(in), optional :: nobounds
      real, dimension(size(vdat),1) :: vdat_l
      logical :: nobounds_l
      
      nobounds_l = .false.
      if ( present(nobounds) ) then
         nobounds_l = nobounds
      end if
      
      vdat_l(:,1) = vdat(:)
      call mgbcasta3( g, vdat_l, nobounds=nobounds_l )
      vdat(:) = vdat_l(:,1)
   
   end subroutine mgbcasta2
   
   subroutine mgbcasta3(g,vdat,nobounds)
      ! This subroutine broadcasts data when downscaling with the multi-grid solver
      integer, intent(in) :: g
      integer :: kx, out_len
      integer(kind=4) :: ierr, ilen, lcomm
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
      call START_LOG(bcast_begin)
#ifdef i8r8
      call MPI_Bcast( tdat, ilen, MPI_DOUBLE_PRECISION, 0_4, lcomm, ierr )
#else
      call MPI_Bcast( tdat, ilen, MPI_REAL, 0_4, lcomm, ierr )
#endif
      call END_LOG(bcast_end)

      ! extract data from Bcast
      vdat(1:out_len,1:kx) = reshape( tdat(1:ilen), (/ out_len, kx /) )
   
   end subroutine mgbcasta3
   
   subroutine mgbcastxn3(g,vdat,smaxmin,klim,nobounds)
      ! This subroutine broadcasts data when downscaling with the multi-grid solver
      ! This version also updates the range (smaxmin)
      integer, intent(in) :: g
      integer, intent(in), optional :: klim
      integer :: kx, out_len
      integer(kind=4) :: ierr, ilen, lcomm
      real, dimension(:,:), intent(inout) :: vdat
      real, dimension(:,:), intent(inout) :: smaxmin
      real, dimension((mg(g)%ifull+mg(g)%iextra+2)*size(vdat,2)) :: tdat
      logical, intent(in), optional :: nobounds
      logical :: nbflag
      
      if ( mg(g)%merge_len <= 1 ) return

      kx = size(vdat, 2)

      nbflag = .false.
      if (present(klim)) then
         kx = klim
      end if
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
      call START_LOG(bcast_begin)
#ifdef i8r8
      call MPI_Bcast( tdat, ilen, MPI_DOUBLE_PRECISION, 0_4, lcomm, ierr )
#else
      call MPI_Bcast( tdat, ilen, MPI_REAL, 0_4, lcomm, ierr )
#endif
      call END_LOG(bcast_end)

      ! extract data from Bcast
      vdat(1:out_len,1:kx) = reshape( tdat(1:out_len*kx), (/ out_len, kx /) )
      smaxmin(1:kx,1:2) = reshape( tdat(out_len*kx+1:ilen), (/ kx, 2 /) )
   
   end subroutine mgbcastxn3
   
   ! Set up the indices required for the multigrid scheme.
   subroutine mg_index(g,mil_g,mipan,mjpan)
      use indices_m
      integer, intent(in) :: g, mil_g, mipan, mjpan
      integer, dimension(2*(mipan+mjpan+2)*(npanels+1)) :: dum
      integer, dimension(:,:), allocatable :: dums, dumr
      integer :: mioff, mjoff
      integer :: i, j, n, iq, iqq, iqg, mfull_g
      integer :: iext, iproc, xlen, jx, nc, xlev, rproc, sproc
      integer :: ntest, ncount, icol, mycol
      integer(kind=4) :: itag=22, lproc, ierr, llen, lcomm
      logical lflag, lglob      

      ! size of this grid
      mfull_g = 6*mil_g**2


      ! calculate process map in iq coordinates
      lglob = .true.
      lflag = .true.
      mioff = 0
      mjoff = 0
      do n = 0,npanels
         do j = 1,mil_g
            do i = 1,mil_g
               iq = i + (j-1)*mil_g + n*mil_g**2 
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
      
      ! index=0 is for all colours
      do n = 0,nproc-1
         mg_bnds(n,g)%len = 0
         allocate( mg_bnds(n,g)%rlen_bg(0:maxcolour), mg_bnds(n,g)%rlen_fn(0:maxcolour) )
         allocate( mg_bnds(n,g)%slen_bg(0:maxcolour), mg_bnds(n,g)%slen_fn(0:maxcolour) )
         allocate( mg_bnds(n,g)%rlenx_bg(0:maxcolour), mg_bnds(n,g)%rlenx_fn(0:maxcolour) )
         allocate( mg_bnds(n,g)%slenx_bg(0:maxcolour), mg_bnds(n,g)%slenx_fn(0:maxcolour) )
         mg_bnds(n,g)%rlen_bg(0:maxcolour) = 0
         mg_bnds(n,g)%rlen_fn(0:maxcolour) = 0
         mg_bnds(n,g)%slen_bg(0:maxcolour) = 0
         mg_bnds(n,g)%slen_fn(0:maxcolour) = 0
         mg_bnds(n,g)%rlenx_bg(0:maxcolour) = 0
         mg_bnds(n,g)%rlenx_fn(0:maxcolour) = 0
         mg_bnds(n,g)%slenx_bg(0:maxcolour) = 0
         mg_bnds(n,g)%slenx_fn(0:maxcolour) = 0
      end do   
         
      ! Calculate local indices on this process
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
         
         ! calculate colours for global level
         allocate( mg(g)%ifull_colour(maxcolour) )
         mg(g)%ifull_colour(:) = 0
         do iq = 1,mfull_g
            mycol = findcolour(iq,mil_g) 
            mg(g)%ifull_colour(mycol) = mg(g)%ifull_colour(mycol) + 1
         end do
         mg(g)%ifull_maxcolour = maxval( mg(g)%ifull_colour(:) )
         allocate( mg(g)%iqx(mg(g)%ifull_maxcolour,maxcolour) )
         mg(g)%ifull_colour(:) = 0
         do iq = 1,mfull_g
            mycol = findcolour(iq,mil_g) 
            mg(g)%ifull_colour(mycol) = mg(g)%ifull_colour(mycol) + 1
            mg(g)%iqx(mg(g)%ifull_colour(mycol),mycol) = iq
         end do
         
      else

         mg(g)%iextra = 2*(mipan+mjpan+4)*npan ! first guess
          
         ! This only occurs with grids prior to globgath.  So npan and noff are still valid.
         do n = 1,npan
            do j = 1,mjpan
               do i = 1,mipan
                  iq = i + (j-1)*mipan + (n-1)*mipan*mjpan              ! Local
                  iqg = i+mioff + (j+mjoff-1)*mil_g + (n-noff)*mil_g**2 ! Global

                  iqq = jn_g(iqg,mil_g)         ! Global neighbour index
                  if ( mg_qproc(iqq,mil_g,g) == myid ) then ! Just copy the value
                     ! Convert global iqq to local value
                     mg(g)%in(iq) = indx_indv(iqq,mil_g,mipan,mjpan,mioff,mjoff,noff)
                  end if

                  iqq = js_g(iqg,mil_g)         ! Global neighbour index
                  if ( mg_qproc(iqq,mil_g,g) == myid ) then ! Just copy the value
                     ! Convert global iqq to local value
                     mg(g)%is(iq) = indx_indv(iqq,mil_g,mipan,mjpan,mioff,mjoff,noff)
                  end if

                  iqq = je_g(iqg,mil_g)         ! Global neighbour index
                  if ( mg_qproc(iqq,mil_g,g) == myid ) then ! Just copy the value
                     ! Convert global iqq to local value
                     mg(g)%ie(iq) = indx_indv(iqq,mil_g,mipan,mjpan,mioff,mjoff,noff)
                  end if

                  iqq = jw_g(iqg,mil_g)         ! Global neighbour index
                  if ( mg_qproc(iqq,mil_g,g) == myid ) then ! Just copy the value
                     ! Convert global iqq to local value
                     mg(g)%iw(iq) = indx_indv(iqq,mil_g,mipan,mjpan,mioff,mjoff,noff)
                  end if
            
                  iqq = jne_g(iqg,mil_g)        ! Global neighbour index
                  if ( mg_qproc(iqq,mil_g,g) == myid ) then ! Just copy the value
                     ! Convert global iqq to local value
                     mg(g)%ine(iq) = indx_indv(iqq,mil_g,mipan,mjpan,mioff,mjoff,noff)
                  end if      

                  iqq = jen_g(iqg,mil_g)        ! Global neighbour index
                  if ( mg_qproc(iqq,mil_g,g) == myid ) then ! Just copy the value
                     ! Convert global iqq to local value
                     mg(g)%ien(iq) = indx_indv(iqq,mil_g,mipan,mjpan,mioff,mjoff,noff)
                  end if

                  iqq = jnw_g(iqg,mil_g)        ! Global neighbour index
                  if ( mg_qproc(iqq,mil_g,g) == myid ) then ! Just copy the value
                     ! Convert global iqq to local value
                     mg(g)%inw(iq) = indx_indv(iqq,mil_g,mipan,mjpan,mioff,mjoff,noff)
                  end if
 
                  iqq = jwn_g(iqg,mil_g)        ! Global neighbour index
                  if ( mg_qproc(iqq,mil_g,g) == myid ) then ! Just copy the value
                     ! Convert global iqq to local value
                     mg(g)%iwn(iq) = indx_indv(iqq,mil_g,mipan,mjpan,mioff,mjoff,noff)
                  end if

                  iqq = jse_g(iqg,mil_g)        ! Global neighbour index
                  if ( mg_qproc(iqq,mil_g,g) == myid ) then ! Just copy the value
                     ! Convert global iqq to local value
                     mg(g)%ise(iq) = indx_indv(iqq,mil_g,mipan,mjpan,mioff,mjoff,noff)
                  end if

                  iqq = jes_g(iqg,mil_g)        ! Global neighbour index
                  if ( mg_qproc(iqq,mil_g,g) == myid ) then ! Just copy the value
                     ! Convert global iqq to local value
                     mg(g)%ies(iq) = indx_indv(iqq,mil_g,mipan,mjpan,mioff,mjoff,noff)
                  end if

                  iqq = jsw_g(iqg,mil_g)        ! Global neighbour index
                  if ( mg_qproc(iqq,mil_g,g) == myid ) then ! Just copy the value
                     ! Convert global iqq to local value
                     mg(g)%isw(iq) = indx_indv(iqq,mil_g,mipan,mjpan,mioff,mjoff,noff)
                  end if
 
                  iqq = jws_g(iqg,mil_g)        ! Global neighbour index
                  if ( mg_qproc(iqq,mil_g,g) == myid ) then ! Just copy the value
                     ! Convert global iqq to local value
                     mg(g)%iws(iq) = indx_indv(iqq,mil_g,mipan,mjpan,mioff,mjoff,noff)
                  end if
                  
               end do
            end do
         end do

         
         ! Calculate local indices in halo
         iext = 0
         
         do n = 0,nproc-1
            mg_bnds(n,g)%rlen_bg(1) = 1
         end do
         
         do icol = 1,maxcolour
         
            do n = 1,npan

               !     Start with N edge
               j = mjpan
               do i = 1,mipan
                  iq = i + (j-1)*mipan + (n-1)*mipan*mjpan              !  Local index 
                  iqg = i+mioff + (j+mjoff-1)*mil_g + (n-noff)*mil_g**2 !  Global index
                  iqq = jn_g(iqg,mil_g)
                  ! Which process has this point
                  rproc = mg_qproc(iqq,mil_g,g)
                  if ( rproc /= myid ) then ! Don't add points already on this proc.
                     mycol = findcolour(iqq,mil_g)
                     if ( mycol == icol ) then
                        ! Add this point to request list
                        mg_bnds(rproc,g)%rlen_fn(icol) = mg_bnds(rproc,g)%rlen_fn(icol) + 1
                        call mgcheck_bnds_alloc(g, rproc, iext)
                        mg_bnds(rproc,g)%request_list(mg_bnds(rproc,g)%rlen_fn(icol)) = iqq
                        ! Increment extended region index
                        iext = iext + 1
                        mg_bnds(rproc,g)%unpack_list(mg_bnds(rproc,g)%rlen_fn(icol)) = iext
                        mg(g)%in(iq) = mg(g)%ifull + iext
                     end if   
                  end if   
               end do

               !     E edge
               i = mipan
               do j = 1,mjpan
                  iq = i + (j-1)*mipan + (n-1)*mipan*mjpan              !  Local index 
                  iqg = i+mioff + (j+mjoff-1)*mil_g + (n-noff)*mil_g**2 !  Global index
                  iqq = je_g(iqg,mil_g)
                  ! Which process has this point
                  rproc = mg_qproc(iqq,mil_g,g)
                  if ( rproc /= myid ) then ! Don't add points already on this proc.
                     mycol = findcolour(iqq,mil_g)
                     if ( mycol == icol ) then
                        ! Add this point to request list
                        mg_bnds(rproc,g)%rlen_fn(icol) = mg_bnds(rproc,g)%rlen_fn(icol) + 1
                        call mgcheck_bnds_alloc(g, rproc, iext)
                        mg_bnds(rproc,g)%request_list(mg_bnds(rproc,g)%rlen_fn(icol)) = iqq
                        ! Increment extended region index
                        iext = iext + 1
                        mg_bnds(rproc,g)%unpack_list(mg_bnds(rproc,g)%rlen_fn(icol)) = iext
                        mg(g)%ie(iq) = mg(g)%ifull + iext
                     end if   
                  end if   
               end do

               !     W edge
               i = 1
               do j = 1,mjpan
                  iq = i + (j-1)*mipan + (n-1)*mipan*mjpan              !  Local index 
                  iqg = i+mioff + (j+mjoff-1)*mil_g + (n-noff)*mil_g**2 !  Global index
                  iqq = jw_g(iqg,mil_g)
                  ! Which process has this point
                  rproc = mg_qproc(iqq,mil_g,g)
                  if ( rproc /= myid ) then ! Don't add points already on this proc.
                     mycol = findcolour(iqq,mil_g)
                     if ( mycol == icol ) then
                        ! Add this point to request list
                        mg_bnds(rproc,g)%rlen_fn(icol) = mg_bnds(rproc,g)%rlen_fn(icol) + 1
                        call mgcheck_bnds_alloc(g, rproc, iext)
                        mg_bnds(rproc,g)%request_list(mg_bnds(rproc,g)%rlen_fn(icol)) = iqq
                        ! Increment extended region index
                        iext = iext + 1
                        mg_bnds(rproc,g)%unpack_list(mg_bnds(rproc,g)%rlen_fn(icol)) = iext
                        mg(g)%iw(iq) = mg(g)%ifull + iext
                     end if   
                  end if   
               end do

               !     S edge
               j = 1
               do i = 1,mipan
                  iq = i + (j-1)*mipan + (n-1)*mipan*mjpan              !  Local index 
                  iqg = i+mioff + (j+mjoff-1)*mil_g + (n-noff)*mil_g**2 !  Global index
                  iqq = js_g(iqg,mil_g)
                  ! Which process has this point
                  rproc = mg_qproc(iqq,mil_g,g)
                  if ( rproc /= myid ) then ! Don't add points already on this proc.
                     mycol = findcolour(iqq,mil_g)
                     if ( mycol == icol ) then
                        ! Add this point to request list
                        mg_bnds(rproc,g)%rlen_fn(icol) = mg_bnds(rproc,g)%rlen_fn(icol) + 1    
                        call mgcheck_bnds_alloc(g, rproc, iext)
                        mg_bnds(rproc,g)%request_list(mg_bnds(rproc,g)%rlen_fn(icol)) = iqq
                        ! Increment extended region index
                        iext = iext + 1
                        mg_bnds(rproc,g)%unpack_list(mg_bnds(rproc,g)%rlen_fn(icol)) = iext
                        mg(g)%is(iq) = mg(g)%ifull + iext
                     end if   
                  end if  
               end do
            
            end do ! n=1,npan
            
            if ( icol < maxcolour ) then
               do n = 0,nproc-1
                  mg_bnds(n,g)%rlen_bg(min(icol+1,maxcolour)) = mg_bnds(n,g)%rlen_fn(icol) + 1
                  mg_bnds(n,g)%rlen_fn(min(icol+1,maxcolour)) = mg_bnds(n,g)%rlen_fn(icol)
               end do
            end if
            
         end do ! icol = 1,maxcolour   
         
         do n = 0,nproc-1
            mg_bnds(n,g)%rlenx_bg(1) = mg_bnds(n,g)%rlen_fn(maxcolour) + 1
            mg_bnds(n,g)%rlenx_fn(1) = mg_bnds(n,g)%rlen_fn(maxcolour)
         end do
         
         do icol = 1,maxcolour
         
            do n = 1,npan
            
               ! NE, EN
               iq = mipan + (mjpan-1)*mipan + (n-1)*mipan*mjpan              !  Local index 
               iqg = mipan+mioff + (mjpan+mjoff-1)*mil_g + (n-noff)*mil_g**2 !  Global index
               iqq = jne_g(iqg,mil_g)
               ! Which process has this point
               rproc = mg_qproc(iqq,mil_g,g)
               if ( rproc /= myid ) then ! Don't add points already on this proc.
                  mycol = findcolour(iqq,mil_g)
                  if ( mycol == icol ) then
                     ! Add this point to request list
                     mg_bnds(rproc,g)%rlenx_fn(icol) = mg_bnds(rproc,g)%rlenx_fn(icol) + 1 
                     call mgcheck_bnds_alloc(g, rproc, iext)
                     mg_bnds(rproc,g)%request_list(mg_bnds(rproc,g)%rlenx_fn(icol)) = iqq
                     ! Increment extended region index
                     iext = iext + 1
                     mg_bnds(rproc,g)%unpack_list(mg_bnds(rproc,g)%rlenx_fn(icol)) = iext
                     mg(g)%ine(iq) = mg(g)%ifull + iext
                  end if
               end if
               if ( jen_g(iqg,mil_g) == jne_g(iqg,mil_g) ) then
                  mg(g)%ien(iq) = mg(g)%ine(iq)
               else
                  iqq = jen_g(iqg,mil_g)
                  ! Which process has this point
                  rproc = mg_qproc(iqq,mil_g,g)
                  if ( rproc /= myid ) then ! Add to list
                     mycol = findcolour(iqq,mil_g)
                     if ( mycol == icol ) then
                        ! Add this point to request list
                        mg_bnds(rproc,g)%rlenx_fn(icol) = mg_bnds(rproc,g)%rlenx_fn(icol) + 1 
                        call mgcheck_bnds_alloc(g, rproc, iext)
                        mg_bnds(rproc,g)%request_list(mg_bnds(rproc,g)%rlenx_fn(icol)) = iqq
                        ! Increment extended region index
                        iext = iext + 1
                        mg_bnds(rproc,g)%unpack_list(mg_bnds(rproc,g)%rlenx_fn(icol)) = iext
                        mg(g)%ien(iq) = mg(g)%ifull + iext
                     end if   
                  end if
               end if

               ! SE, ES
               iq = mipan + (n-1)*mipan*mjpan                      !  Local index
               iqg = mipan+mioff + mjoff*mil_g + (n-noff)*mil_g**2 !  Global index
               iqq = jse_g(iqg,mil_g)
               ! Which process has this point
               rproc = mg_qproc(iqq,mil_g,g)
               if ( rproc /= myid ) then ! Don't add points already on this proc.
                  mycol = findcolour(iqq,mil_g)
                  if ( mycol == icol ) then
                     ! Add this point to request list
                     mg_bnds(rproc,g)%rlenx_fn(icol) = mg_bnds(rproc,g)%rlenx_fn(icol) + 1 
                     call mgcheck_bnds_alloc(g, rproc, iext)
                     mg_bnds(rproc,g)%request_list(mg_bnds(rproc,g)%rlenx_fn(icol)) = iqq
                     ! Increment extended region index
                     iext = iext + 1
                     mg_bnds(rproc,g)%unpack_list(mg_bnds(rproc,g)%rlenx_fn(icol)) = iext
                     mg(g)%ise(iq) = mg(g)%ifull + iext
                  end if   
               end if  
               if ( jes_g(iqg,mil_g) == jse_g(iqg,mil_g) ) then
                  mg(g)%ies(iq) = mg(g)%ise(iq)
               else
                  iqq = jes_g(iqg,mil_g)
                  ! Which process has this point
                  rproc = mg_qproc(iqq,mil_g,g)
                  if ( rproc /= myid ) then ! Add to list
                     mycol = findcolour(iqq,mil_g)
                     if ( mycol == icol ) then
                        ! Add this point to request list
                        mg_bnds(rproc,g)%rlenx_fn(icol) = mg_bnds(rproc,g)%rlenx_fn(icol) + 1    
                        call mgcheck_bnds_alloc(g, rproc, iext)
                        mg_bnds(rproc,g)%request_list(mg_bnds(rproc,g)%rlenx_fn(icol)) = iqq
                        ! Increment extended region index
                        iext = iext + 1
                        mg_bnds(rproc,g)%unpack_list(mg_bnds(rproc,g)%rlenx_fn(icol)) = iext
                        mg(g)%ies(iq) = mg(g)%ifull + iext
                     end if   
                  end if
               end if

               ! NW, WN
               iq = 1 + (mjpan-1)*mipan + (n-1)*mipan*mjpan              !  Local index
               iqg = 1+mioff + (mjpan+mjoff-1)*mil_g + (n-noff)*mil_g**2 !  Global index
               iqq = jnw_g(iqg,mil_g)
               ! Which process has this point
               rproc = mg_qproc(iqq,mil_g,g)
               if ( rproc /= myid ) then ! Don't add points already on this proc.
                  mycol = findcolour(iqq,mil_g)
                  if ( mycol == icol ) then
                     ! Add this point to request list
                     mg_bnds(rproc,g)%rlenx_fn(icol) = mg_bnds(rproc,g)%rlenx_fn(icol) + 1 
                     call mgcheck_bnds_alloc(g, rproc, iext)
                     mg_bnds(rproc,g)%request_list(mg_bnds(rproc,g)%rlenx_fn(icol)) = iqq
                     ! Increment extended region index
                     iext = iext + 1
                     mg_bnds(rproc,g)%unpack_list(mg_bnds(rproc,g)%rlenx_fn(icol)) = iext
                     mg(g)%inw(iq) = mg(g)%ifull + iext
                  end if   
               end if  
               if ( jwn_g(iqg,mil_g) == jnw_g(iqg,mil_g) ) then
                  mg(g)%iwn(iq) = mg(g)%inw(iq)
               else
                  iqq = jwn_g(iqg,mil_g)
                  ! Which process has this point
                  rproc = mg_qproc(iqq,mil_g,g)
                  if ( rproc /= myid ) then ! Add to list
                     mycol = findcolour(iqq,mil_g)
                     if ( mycol == icol ) then
                        ! Add this point to request list
                        mg_bnds(rproc,g)%rlenx_fn(icol) = mg_bnds(rproc,g)%rlenx_fn(icol) + 1
                        call mgcheck_bnds_alloc(g, rproc, iext)
                        mg_bnds(rproc,g)%request_list(mg_bnds(rproc,g)%rlenx_fn(icol)) = iqq
                        ! Increment extended region index
                        iext = iext + 1
                        mg_bnds(rproc,g)%unpack_list(mg_bnds(rproc,g)%rlenx_fn(icol)) = iext
                        mg(g)%iwn(iq) = mg(g)%ifull + iext
                     end if   
                  end if
               end if

               ! SW, WS
               iq = 1 + (n-1)*mipan*mjpan                        !  Local index
               iqg = 1+mioff + (mjoff)*mil_g + (n-noff)*mil_g**2 !  Global index
               iqq = jsw_g(iqg,mil_g)
               ! Which process has this point
               rproc = mg_qproc(iqq,mil_g,g)
               if ( rproc /= myid ) then ! Don't add points already on this proc.
                  mycol = findcolour(iqq,mil_g)
                  if ( mycol == icol ) then
                     ! Add this point to request list
                     mg_bnds(rproc,g)%rlenx_fn(icol) = mg_bnds(rproc,g)%rlenx_fn(icol) + 1    
                     call mgcheck_bnds_alloc(g, rproc, iext)
                     mg_bnds(rproc,g)%request_list(mg_bnds(rproc,g)%rlenx_fn(icol)) = iqq
                     ! Increment extended region index
                     iext = iext + 1
                     mg_bnds(rproc,g)%unpack_list(mg_bnds(rproc,g)%rlenx_fn(icol)) = iext
                     mg(g)%isw(iq) = mg(g)%ifull + iext
                  end if   
               end if  
               if ( jws_g(iqg,mil_g) == jsw_g(iqg,mil_g) ) then
                  mg(g)%iws(iq) = mg(g)%isw(iq)
               else
                  iqq = jws_g(iqg,mil_g)
                  ! Which process has this point
                  rproc = mg_qproc(iqq,mil_g,g)
                  if ( rproc /= myid ) then ! Add to list
                     mycol = findcolour(iqq,mil_g)
                     if ( mycol == icol ) then
                        ! Add this point to request list
                        mg_bnds(rproc,g)%rlenx_fn(icol) = mg_bnds(rproc,g)%rlenx_fn(icol) + 1 
                        call mgcheck_bnds_alloc(g, rproc, iext)
                        mg_bnds(rproc,g)%request_list(mg_bnds(rproc,g)%rlenx_fn(icol)) = iqq
                        ! Increment extended region index
                        iext = iext + 1
                        mg_bnds(rproc,g)%unpack_list(mg_bnds(rproc,g)%rlenx_fn(icol)) = iext
                        mg(g)%iws(iq) = mg(g)%ifull + iext
                     end if
                  end if   
               end if
            
            end do ! n=1,npan
               
            if ( icol < maxcolour ) then
               do n = 0,nproc-1
                  mg_bnds(n,g)%rlenx_bg(min(icol+1,maxcolour)) = mg_bnds(n,g)%rlenx_fn(icol) + 1
                  mg_bnds(n,g)%rlenx_fn(min(icol+1,maxcolour)) = mg_bnds(n,g)%rlenx_fn(icol)
               end do
            end if
            
         end do ! icol = 1,maxcolour   


         mg(g)%iextra = iext

        
         ! Set up the diagonal index arrays
         do n = 1,npan
            do j = 1,mjpan
               do i = 1,mipan
                  iq = i + (j-1)*mipan + (n-1)*mipan*mjpan !  Local index
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
         mg(g)%neighnum = 0
         do iproc = 1,nproc-1
            rproc = modulo( myid+iproc, nproc ) 
            if ( mg_bnds(rproc,g)%rlenx_fn(maxcolour) > 0 ) then 
              mg(g)%neighnum = mg(g)%neighnum + 1
            end if
         end do

         ! increase size of request list if needed
         ntest = mg(g)%neighnum
         if ( 2*ntest > size(ireq) ) then
            deallocate( ireq, rlist )
            allocate( ireq(max(2*ntest,1)) )
            allocate( rlist(max(ntest,1)) )
         end if
  
         ! set up neighbour lists
         allocate ( mg(g)%neighlist(mg(g)%neighnum) )
         ncount = 0
         do iproc = 1,nproc-1
            rproc = modulo( myid+iproc, nproc )
            if ( mg_bnds(rproc,g)%rlenx_fn(maxcolour) > 0 ) then
               ncount = ncount + 1
               mg(g)%neighlist(ncount) = rproc
            end if
         end do
      
         if ( ncount/=mg(g)%neighnum ) then
            write(6,*) "ERROR: Multi-grid neighnum mismatch"
            write(6,*) "myid, neighnum, ncount ",myid,mg(g)%neighnum, ncount
            do n = 0,nproc-1
              if ( mg_bnds(n,g)%rlenx_fn(maxcolour) > 0 ) then
                 write(6,*) "myid,n,rlenx ",myid,n,mg_bnds(n,g)%rlenx_fn(maxcolour) 
              end if
            end do  
            call ccmpi_abort(-1)
         end if

         
         ! Now, for each process send the length of points I want.
         allocate( dums(2*maxcolour,mg(g)%neighnum), dumr(2*maxcolour,mg(g)%neighnum) )
         lcomm = comm_world
         nreq = 0
         do iproc = 1,mg(g)%neighnum
            rproc = mg(g)%neighlist(iproc)  ! Recv from
            if ( mg_bnds(rproc,g)%rlenx_fn(maxcolour) > 0 ) then
               nreq = nreq + 1
               lproc = rproc
#ifdef i8r8
               call MPI_IRecv( dumr(:,iproc), int(2*maxcolour,4), MPI_INTEGER8, lproc, itag, lcomm, ireq(nreq), ierr )
#else
               call MPI_IRecv( dumr(:,iproc), int(2*maxcolour,4), MPI_INTEGER, lproc, itag, lcomm, ireq(nreq), ierr )
#endif
            end if
         end do
         do iproc = mg(g)%neighnum,1,-1
            sproc = mg(g)%neighlist(iproc)  ! Send to
            if ( mg_bnds(sproc,g)%rlenx_fn(maxcolour) > 0 ) then
               nreq = nreq + 1
               dums(1:maxcolour,iproc) = mg_bnds(sproc,g)%rlen_fn(1:maxcolour)
               dums(maxcolour+1:2*maxcolour,iproc) = mg_bnds(sproc,g)%rlenx_fn(1:maxcolour)
               lproc = sproc
#ifdef i8r8
               call MPI_ISend( dums(:,iproc), int(2*maxcolour,4), MPI_INTEGER8, lproc, itag, lcomm, ireq(nreq), ierr )
#else
               call MPI_ISend( dums(:,iproc), int(2*maxcolour,4), MPI_INTEGER, lproc, itag, lcomm, ireq(nreq), ierr )
#endif
            end if
         end do
         call MPI_Waitall( nreq, ireq, MPI_STATUSES_IGNORE, ierr )
         do iproc = 1,mg(g)%neighnum
            rproc = mg(g)%neighlist(iproc)
            if ( mg_bnds(rproc,g)%rlenx_fn(maxcolour) > 0 ) then
               mg_bnds(rproc,g)%slen_fn(1:maxcolour) = dumr(1:maxcolour,iproc)
               mg_bnds(rproc,g)%slenx_fn(1:maxcolour) = dumr(maxcolour+1:2*maxcolour,iproc)
            end if
         end do   
         do n = 0,nproc-1
            mg_bnds(n,g)%slen_bg(1) = 1
            do i = 2,maxcolour
               mg_bnds(n,g)%slen_bg(i) = mg_bnds(n,g)%slen_fn(i-1) + 1
            end do
            mg_bnds(n,g)%slenx_bg(1) = mg_bnds(n,g)%slen_fn(maxcolour) + 1
            do i = 2,maxcolour
               mg_bnds(n,g)%slenx_bg(i) = mg_bnds(n,g)%slenx_fn(i-1) + 1
            end do
         end do
         
         ! create index=0 for all colours
         do n = 0,nproc-1
            mg_bnds(n,g)%rlen_bg(0) = mg_bnds(n,g)%rlen_bg(1)
            mg_bnds(n,g)%rlen_fn(0) = mg_bnds(n,g)%rlen_fn(maxcolour)
            mg_bnds(n,g)%slen_bg(0) = mg_bnds(n,g)%slen_bg(1)
            mg_bnds(n,g)%slen_fn(0) = mg_bnds(n,g)%slen_fn(maxcolour)
            mg_bnds(n,g)%rlenx_bg(0) = mg_bnds(n,g)%rlenx_bg(1)
            mg_bnds(n,g)%rlenx_fn(0) = mg_bnds(n,g)%rlenx_fn(maxcolour)
            mg_bnds(n,g)%slenx_bg(0) = mg_bnds(n,g)%slenx_bg(1)
            mg_bnds(n,g)%slenx_fn(0) = mg_bnds(n,g)%slenx_fn(maxcolour)            
         end do
         
         nreq = 0
         rreq = 0
         deallocate( dums, dumr )

  
         ! Now start sending messages  
         lcomm = comm_world
         nreq = 0
         do iproc = 1,mg(g)%neighnum
            lproc = mg(g)%neighlist(iproc)  ! Recv from
            allocate( mg_bnds(lproc,g)%send_list(mg_bnds(lproc,g)%slenx_fn(maxcolour)) )
            nreq = nreq + 1
            ! Use the maximum size in the recv call.
            llen = mg_bnds(lproc,g)%slenx_fn(maxcolour)
#ifdef i8r8
            call MPI_IRecv( mg_bnds(lproc,g)%send_list(1), llen, MPI_INTEGER8, lproc, &
                            itag, lcomm, ireq(nreq), ierr )
#else
            call MPI_IRecv( mg_bnds(lproc,g)%send_list(1), llen, MPI_INTEGER, lproc, &
                            itag, lcomm, ireq(nreq), ierr )
#endif
         end do
         do iproc = mg(g)%neighnum,1,-1
            lproc = mg(g)%neighlist(iproc)  ! Send to
            ! Send list of requests
            nreq = nreq + 1
            llen = mg_bnds(lproc,g)%rlenx_fn(maxcolour)
#ifdef i8r8
            call MPI_ISend( mg_bnds(lproc,g)%request_list(1), llen, MPI_INTEGER8, lproc, &
                            itag, lcomm, ireq(nreq), ierr )
#else
            call MPI_ISend( mg_bnds(lproc,g)%request_list(1), llen, MPI_INTEGER, lproc, &
                            itag, lcomm, ireq(nreq), ierr )
#endif
         end do      
         call MPI_Waitall( nreq, ireq, MPI_STATUSES_IGNORE, ierr )
         nreq = 0
         rreq = 0

         
         ! At the moment send_lists use global indices. Convert these to local.
         do iproc = mg(g)%neighnum,1,-1
            sproc = mg(g)%neighlist(iproc)  ! Send to
            do iq = 1,mg_bnds(sproc,g)%slenx_fn(maxcolour)
               ! send_list(iq) is global point index, i, j, n are local
               iqq = mg_bnds(sproc,g)%send_list(iq)
               mg_bnds(sproc,g)%send_list(iq) = indx_indv(iqq,mil_g,mipan,mjpan,mioff,mjoff,noff)
            end do
         end do
         if ( mg_bnds(myid,g)%rlenx_fn(maxcolour) /= 0 ) then
            write(6,*) "ERROR: Invalid rlenx in myid"
            call ccmpi_abort(-1)
         end if   

         ! reduce array size where possible
         do iproc = 0,nproc-1
            xlen = mg_bnds(iproc,g)%rlenx_fn(maxcolour)
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
            xlen = xlev*mg_bnds(iproc,g)%rlenx_fn(maxcolour)
            if ( bnds(iproc)%rbuflen < xlen ) then
               if ( allocated(bnds(iproc)%rbuf) ) then
                  deallocate( bnds(iproc)%rbuf )
                  deallocate( bnds(iproc)%r8buf )
               end if
               allocate( bnds(iproc)%rbuf(xlen) )
               allocate( bnds(iproc)%r8buf(xlen) )
               bnds(iproc)%rbuflen = xlen
            end if
            xlen = xlev*mg_bnds(iproc,g)%slenx_fn(maxcolour)
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
      
         ! calculate colours per level
         allocate( mg(g)%ifull_colour(maxcolour) )
         mg(g)%ifull_colour(:) = 0
         do n = 1,npan
            do j = 1,mjpan
               do i = 1,mipan
                  iq = i + (j-1)*mipan + (n-1)*mipan*mjpan              ! Local index
                  iqg = i+mioff + (j+mjoff-1)*mil_g + (n-noff)*mil_g**2 ! Global index
                  mycol = findcolour(iqg,mil_g)
                  mg(g)%ifull_colour(mycol) = mg(g)%ifull_colour(mycol) + 1
               end do
            end do
         end do
         mg(g)%ifull_maxcolour = maxval( mg(g)%ifull_colour(:) )
         allocate( mg(g)%iqx(mg(g)%ifull_maxcolour,maxcolour) )
         mg(g)%ifull_colour(:) = 0
         do n = 1,npan
            do j = 1,mjpan
               do i = 1,mipan
                  iq = i + (j-1)*mipan + (n-1)*mipan*mjpan              ! Local index
                  iqg = i+mioff + (j+mjoff-1)*mil_g + (n-noff)*mil_g**2 ! Global index
                  mycol = findcolour(iqg,mil_g)
                  mg(g)%ifull_colour(mycol) = mg(g)%ifull_colour(mycol) + 1
                  mg(g)%iqx(mg(g)%ifull_colour(mycol),mycol) = iq
               end do
            end do
         end do   
      
      end if   

   end subroutine mg_index

   subroutine mgcheck_bnds_alloc(g,iproc,iext)
      integer, intent(in) :: iproc
      integer, intent(in) :: g, iext
      integer :: testlen, i

      if ( mg_bnds(iproc,g)%len <= 0 ) then
         allocate( mg_bnds(iproc,g)%request_list(mg(g)%iextra) )
         allocate( mg_bnds(iproc,g)%unpack_list(mg(g)%iextra) )
         mg_bnds(iproc,g)%len = mg(g)%iextra
      else
         ! Just check length 
         testlen = 0
         do i = 0,maxcolour
            testlen = max( mg_bnds(iproc,g)%rlen_fn(i), mg_bnds(iproc,g)%rlenx_fn(i), &
                           testlen )
         end do
         if ( testlen > mg_bnds(iproc,g)%len ) then
            write(6,*) "ERROR: MG grid undersized in mgcheck_bnds_alloc" 
            call ccmpi_abort(-1)
         end if
         if ( iext>mg(g)%iextra ) then
            write(6,*) "ERROR: MG grid undersized in mgcheck_bnds_alloc"
            write(6,*) "iext,iextra,g,iproc,myid ",iext,mg(g)%iextra,g,iproc,myid
            call ccmpi_abort(-1)
         end if
      end if

   end subroutine mgcheck_bnds_alloc
   
   subroutine mgbounds2(g,vdat,corner)
      integer, intent(in) :: g
      real, dimension(:), intent(inout) :: vdat
      logical, intent(in), optional :: corner
      real, dimension(size(vdat),1) :: vdat_l
      logical :: corner_l
      
      corner_l = .true.
      if ( present(corner) ) then
         corner_l = corner
      end if

      ! colour=0 is for all grid points
      vdat_l(:,1) = vdat(:)
      call mgbounds_colour3( g, vdat_l, 0, corner=corner_l )
      vdat(:) = vdat_l(:,1)

   end subroutine mgbounds2
 
   subroutine mgbounds3(g,vdat,klim,corner)
      integer, intent(in) :: g
      real, dimension(:,:), intent(inout) :: vdat
      integer, intent(in), optional :: klim
      logical, intent(in), optional :: corner
      integer :: klim_l
      logical :: corner_l
      
      klim_l = size(vdat, 2)
      if ( present(klim) ) then
         klim_l = klim
      end if
      corner_l = .false.
      if ( present(corner) ) then
         corner_l = corner
      end if
      
      ! colour=0 is for all grid points
      call mgbounds_colour3( g, vdat, 0, klim=klim_l, corner=corner_l )

   end subroutine mgbounds3

   subroutine mgbounds_colour3(g,vdat,colour,klim,corner)
      ! update halo for specified colour with multi-grid level g
      integer, intent(in) :: g, colour
      real, dimension(:,:), intent(inout) :: vdat
      integer, intent(in), optional :: klim
      logical, intent(in), optional :: corner
      logical :: extra
      integer :: kx, iproc, recv_len
      integer :: rcount, jproc, mproc, iq, k, iqq, ibeg, iend
      integer(kind=4) :: ierr, itag=20, llen, lproc
      integer(kind=4) :: ldone, lcomm
      integer(kind=4), dimension(2*mg(g)%neighnum) :: donelist

      kx = size(vdat,2)
      extra = .false.
      if (present(klim)) then
         kx = klim
      end if
      if ( present(corner) ) then
         extra = corner
      end if
      
      if ( colour<0 .or. colour>maxcolour ) then
         write(6,*) "ERROR: Invalid colour for mgbounds_colour"
         call ccmpi_abort(-1)
      end if

      !     Set up the buffers to send and recv
      lcomm = comm_world
      nreq = 0
      do iproc = 1,mg(g)%neighnum
         lproc = mg(g)%neighlist(iproc)  ! Recv from
         recv_len = 0
         ibeg = mg_bnds(lproc,g)%rlen_bg(colour)
         iend = mg_bnds(lproc,g)%rlen_fn(colour)
         if ( iend >= ibeg ) then
            recv_len = recv_len + (iend-ibeg+1)*kx
         end if   
         if ( extra ) then
            ibeg = mg_bnds(lproc,g)%rlenx_bg(colour)
            iend = mg_bnds(lproc,g)%rlenx_fn(colour)
            if ( iend >= ibeg ) then
               recv_len = recv_len + (iend-ibeg+1)*kx
            end if   
         end if
         if ( recv_len > 0 ) then
            nreq = nreq + 1
            rlist(nreq) = iproc
            llen = recv_len*kx
#ifdef i8r8
            call MPI_IRecv( bnds(lproc)%rbuf, llen, MPI_DOUBLE_PRECISION, lproc, &
                            itag, lcomm, ireq(nreq), ierr )
#else
            call MPI_IRecv( bnds(lproc)%rbuf, llen, MPI_REAL, lproc, &
                            itag, lcomm, ireq(nreq), ierr )
#endif
         end if
      end do
      rreq = nreq
      do iproc = mg(g)%neighnum,1,-1
         lproc = mg(g)%neighlist(iproc)  ! Send to
         iqq = 0
         ibeg = mg_bnds(lproc,g)%slen_bg(colour)
         iend = mg_bnds(lproc,g)%slen_fn(colour)
         if ( iend >= ibeg ) then
            do k = 1,kx
               do iq = 1,iend-ibeg+1
                  bnds(lproc)%sbuf(iqq+iq+(k-1)*(iend-ibeg+1)) = &
                      vdat(mg_bnds(lproc,g)%send_list(iq+ibeg-1),k)
               end do
            end do   
         end if   
         iqq = iqq + (iend-ibeg+1)*kx
         if ( extra ) then
            ibeg = mg_bnds(lproc,g)%slenx_bg(colour)
            iend = mg_bnds(lproc,g)%slenx_fn(colour)
            if ( iend >= ibeg ) then
               do k = 1,kx
                  do iq = 1,iend-ibeg+1
                     bnds(lproc)%sbuf(iqq+iq+(k-1)*(iend-ibeg+1)) = &
                         vdat(mg_bnds(lproc,g)%send_list(iq+ibeg-1),k)
                  end do
               end do   
            end if
            iqq = iqq + (iend-ibeg+1)*kx
         end if   
         if ( iqq > 0 ) then
            nreq = nreq + 1
            llen = iqq
#ifdef i8r8
            call MPI_ISend( bnds(lproc)%sbuf, llen, MPI_DOUBLE_PRECISION, lproc, &
                            itag, lcomm, ireq(nreq), ierr )
#else
            call MPI_ISend( bnds(lproc)%sbuf, llen, MPI_REAL, lproc, &
                            itag, lcomm, ireq(nreq), ierr )
#endif
         end if
      end do

      rcount = nreq
      do while ( rcount > 0 )
         call START_LOG(mpiwait_begin)
         call MPI_Waitsome( nreq, ireq, ldone, donelist, MPI_STATUSES_IGNORE, ierr )
         call END_LOG(mpiwait_end)
         rcount = rcount - ldone
         do jproc = 1,ldone
            mproc = donelist(jproc)
            if ( mproc <= rreq ) then
               iproc = rlist(mproc)  ! Recv from
               lproc = mg(g)%neighlist(iproc)
               iqq = 0
               ibeg = mg_bnds(lproc,g)%rlen_bg(colour)
               iend = mg_bnds(lproc,g)%rlen_fn(colour)               
               do k = 1,kx
                  do iq = 1,iend-ibeg+1
                     vdat(mg(g)%ifull+mg_bnds(lproc,g)%unpack_list(iq+ibeg-1),k) &
                         = bnds(lproc)%rbuf(iqq+iq+(k-1)*(iend-ibeg+1))
                  end do
               end do
               iqq = iqq + (iend-ibeg+1)*kx
               if ( extra ) then
                  ibeg = mg_bnds(lproc,g)%rlenx_bg(colour)
                  iend = mg_bnds(lproc,g)%rlenx_fn(colour)               
                  do k = 1,kx
                     do iq = 1,iend-ibeg+1
                        vdat(mg(g)%ifull+mg_bnds(lproc,g)%unpack_list(iq+ibeg-1),k) &
                            = bnds(lproc)%rbuf(iqq+iq+(k-1)*(iend-ibeg+1))
                     end do
                  end do
                  iqq = iqq + (iend-ibeg+1)*kx 
               end if    
            end if    ! mproc <= rreq
         end do
      end do

   end subroutine mgbounds_colour3
   
   pure function mg_fproc_1(g,i,j,n) result(mg_fpout)
     ! locates process that owns a global grid point
     integer, intent(in) :: i, j, n, g
     integer mg_fpout, g_l, i_l, j_l
     
     i_l = i
     j_l = j
     do g_l = g-1,1,-1
        i_l = (i_l-1)*2 + 1 
        j_l = (j_l-1)*2 + 1 
     end do
     mg_fpout = fproc(i_l,j_l,n)
     
   end function mg_fproc_1
   
   pure function mg_fproc(g,i,j,n) result(mg_fpout)
     ! locates process that owns a global grid point
     integer, intent(in) :: i, j, n, g
     integer mg_fpout, fp_l
     
     fp_l = mg_fproc_1(g,i,j,n)
     mg_fpout = mg(g)%procmap(fp_l)
     
   end function mg_fproc   
   
   pure function mg_qproc(iqg,mil_g,g) result(mg_qpout)
      ! locates process that owns a global grid point
      integer, intent(in) :: iqg, mil_g, g
      integer :: mg_qpout, i, j, n

      n = (iqg - 1) / (mil_g*mil_g)
      j = 1 + (iqg - n*mil_g*mil_g - 1)/mil_g
      i = iqg - (j - 1)*mil_g - n*mil_g*mil_g
      mg_qpout = mg_fproc(g,i,j,n)
   
   end function mg_qproc

end module cc_mpi_multigrid

