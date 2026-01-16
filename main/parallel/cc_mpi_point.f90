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
    
! Subroutines for point-to-point commuication    
    
module cc_mpi_point

   use cc_mpi_common
   
   private
   
   integer, allocatable, dimension(:), save, public :: neighlist           ! list of neighbour processes
   integer, allocatable, dimension(:), save, private :: neighmap           ! map of process to neighbour index
   integer, save, public :: neighnum                                       ! number of neighbours

   public :: bounds_setup 
   public :: bounds, boundsuv, bounds_colour_send, bounds_colour_recv,      &
             bounds_send, bounds_recv, boundsr8, deptsync, intssync_send,   &
             intssync_recv

   interface bounds
      module procedure bounds2, bounds3, bounds4
   end interface
   interface boundsr8
      module procedure bounds2r8, bounds3r8, bounds4r8
   end interface
   interface bounds_send
      module procedure bounds_send3, bounds_send4
   end interface
   interface bounds_recv
      module procedure bounds_recv3, bounds_recv4
   end interface
   interface bounds_colour_send
      module procedure bounds_colour_send3, bounds_colour_send4
   end interface
   interface bounds_colour_recv
      module procedure bounds_colour_recv3, bounds_colour_recv4
   end interface
   interface boundsuv
      module procedure boundsuv2, boundsuv3
   end interface
   interface intssync_send
      module procedure intssync_send3, intssync_send4
   end interface
   interface intssync_recv
      module procedure intssync_recv3, intssync_recv4
   end interface

   type dpoints_info
      real, dimension(:,:), allocatable :: a
   end type dpoints_info
   type dbuf_info
      real, dimension(:,:), allocatable :: a
      real, dimension(:), allocatable :: b
   end type dbuf_info   
   type dindex_info
      integer, dimension(:,:), allocatable :: a
   end type dindex_info
   type sextra_info
      real, dimension(:), allocatable :: a
   end type sextra_info
   
   ! Off process departure points
   type(dpoints_info), allocatable, dimension(:), save, public :: dpoints ! request list from other proc
   type(dbuf_info), allocatable, dimension(:), save, public :: dbuf       ! recv buffer
   type(dindex_info), allocatable, dimension(:), save, public :: dindex   ! request list for my proc
   type(sextra_info), allocatable, dimension(:), save, public :: sextra   ! send buffer
   ! Number of points for each process.
   integer, dimension(:), allocatable, save, public :: dslen
   integer, dimension(:), allocatable, save, public :: drlen

contains

   subroutine bounds_setup(dt)

      use const_phys, only : rearth
      use indices_m
      
      integer :: n, i, j, iq, iqq, mycol, ncount
      integer :: iproc, rproc, sproc
      integer :: iqg, iql, icol
      integer :: iext, iextu, iextv
      integer, dimension(:), allocatable :: dumi
      integer, dimension(:,:), allocatable :: dums, dumr
      integer(kind=4) :: ierr, itag=0
      integer(kind=4) :: llen, lproc, lcomm
      real :: maxdis
      real, intent(in) :: dt
      logical :: swap
      logical, dimension(:), allocatable :: neigharray_g
      logical(kind=4), dimension(:,:), allocatable :: dumsl, dumrl

      ! Just set values that point to values within own process region.
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

      if ( myid==0 ) then
         write(6,*) "-> Calculate local direction index for scalars"
      end if
      do n = 1,npan
         do j = 1,jpan
            do i = 1,ipan
               iq = indp(i,j,n)   ! Local
               iqg = indg(i,j,n)  ! Global

               iqq = in_g(iqg)    ! Global neighbour index
               if ( qproc(iqq) == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  in(iq) = iqg2iq(iqq)
               end if
               iqq = inn_g(iqg)   ! Global neighbour index
               if ( qproc(iqq) == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  inn(iq) = iqg2iq(iqq)
               end if

               iqq = is_g(iqg)    ! Global neighbour index
               if ( qproc(iqq) == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  is(iq) = iqg2iq(iqq)
               end if
               iqq = iss_g(iqg)    ! Global neighbour index
               if ( qproc(iqq) == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  iss(iq) = iqg2iq(iqq)
               end if

               iqq = ie_g(iqg)    ! Global neighbour index
               if ( qproc(iqq) == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  ie(iq) = iqg2iq(iqq)
               end if
               iqq = iee_g(iqg)    ! Global neighbour index
               if ( qproc(iqq) == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  iee(iq) = iqg2iq(iqq)
               end if

               iqq = iw_g(iqg)    ! Global neighbour index
               if ( qproc(iqq) == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  iw(iq) = iqg2iq(iqq)
               end if
               iqq = iww_g(iqg)    ! Global neighbour index
               if ( qproc(iqq) == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  iww(iq) = iqg2iq(iqq)
               end if

               ! Note that the model only needs a limited set of the diagonal
               ! index arrays
               iqq = ine_g(iqg)    ! Global neighbour index
               if ( qproc(iqq) == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  ine(iq) = iqg2iq(iqq)
               end if

               iqq = ise_g(iqg)    ! Global neighbour index
               if ( qproc(iqq) == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  ise(iq) = iqg2iq(iqq)
               end if

               iqq = ien_g(iqg)    ! Global neighbour index
               if ( qproc(iqq) == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  ien(iq) = iqg2iq(iqq)
               end if

               iqq = iwn_g(iqg)    ! Global neighbour index
               if ( qproc(iqq) == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  iwn(iq) = iqg2iq(iqq)
               end if

               iqq = inw_g(iqg)    ! Global neighbour index
               if ( qproc(iqq) == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  inw(iq) = iqg2iq(iqq)
               end if

               iqq = isw_g(iqg)    ! Global neighbour index
               if ( qproc(iqq) == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  isw(iq) = iqg2iq(iqq)
               end if

               iqq = ies_g(iqg)    ! Global neighbour index
               if ( qproc(iqq) == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  ies(iq) = iqg2iq(iqq)
               end if

               iqq = iws_g(iqg)    ! Global neighbour index
               if ( qproc(iqq) == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  iws(iq) = iqg2iq(iqq)
               end if

            end do
         end do
      end do

      ! Correct within the same face only (not necessarily the same
      ! process, but will be corrected later).
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
         ! index=0 is for all coloured grid points
         bnds(n)%rlenh_bg(0:maxcolour) = 0
         bnds(n)%rlenh_fn(0:maxcolour) = 0
         bnds(n)%slenh_bg(0:maxcolour) = 0
         bnds(n)%slenh_fn(0:maxcolour) = 0
         bnds(n)%rlen_bg(0:maxcolour) = 0
         bnds(n)%rlen_fn(0:maxcolour) = 0
         bnds(n)%slen_bg(0:maxcolour) = 0
         bnds(n)%slen_fn(0:maxcolour) = 0
         bnds(n)%rlenx_bg(0:maxcolour) = 0
         bnds(n)%rlenx_fn(0:maxcolour) = 0
         bnds(n)%slenx_bg(0:maxcolour) = 0
         bnds(n)%slenx_fn(0:maxcolour) = 0
         bnds(n)%rlen2 = 0
         bnds(n)%slen2 = 0
         bnds(n)%rlen_su_bg = 0
         bnds(n)%rlen_ev_fn = 0
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
         bnds(n)%slen_sv_bg = 0
         bnds(n)%slen_wu_fn = 0
         bnds(n)%slen_nv_bg = 0
         bnds(n)%slen_eu_fn = 0
         bnds(n)%slen_ssv_bg = 0
         bnds(n)%slen_wwu_fn = 0
         bnds(n)%slen_nnv_bg = 0
         bnds(n)%slen_eeu_fn = 0    
      end do   

      if ( myid==0 ) then
         write(6,*) "-> Calculate neighbour direction index for scalars"
      end if
!     In the first pass through, set up list of points to be requested from
!     other processes. These points are placed in the "iextra" region at the
!     end of arrays. The indirect indices are updated to point into this 
!     region.
      iext = 0
 
      do n = 0,nproc-1
         bnds(n)%rlenh_bg(1) = 1
      end do   
      
      do icol = 1,maxcolour
      
         do n = 1,npan

            !     Start with N edge
            j = jpan
            do i = 1,ipan
               iq = indg(i,j,n)
               iqq = in_g(iq)
               ! Which process has this point
               rproc = qproc(iqq)
               if ( rproc /= myid ) then ! Don't add points already on this proc.
                  mycol = findcolour(iqq,il_g)
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
               ! Which process has this point
               rproc = qproc(iqq)
               if ( rproc /= myid ) then ! Don't add points already on this proc.
                  mycol = findcolour(iqq,il_g)
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
            do n = 0,nproc-1 
               bnds(n)%rlenh_bg(min(icol+1,maxcolour)) = bnds(n)%rlenh_fn(icol) + 1
               bnds(n)%rlenh_fn(min(icol+1,maxcolour)) = bnds(n)%rlenh_fn(icol)
            end do   
         end if
      
      end do ! icol=1,maxcolour

      do n = 0,nproc-1
         bnds(n)%rlen_bg(1) = bnds(n)%rlenh_fn(maxcolour) + 1
         bnds(n)%rlen_fn(1) = bnds(n)%rlenh_fn(maxcolour)
      end do   
      
      do icol = 1,maxcolour
      
         do n = 1,npan

            !     W edge
            i = 1
            do j = 1,jpan
               ! 1D code takes care of corners separately at the end so only goes
               ! over 1,jpan here.
               iq = indg(i,j,n)
               iqq = iw_g(iq) ! iqq is the global index of the required neighbouring point.
               ! Which process has this point
               rproc = qproc(iqq)
               if ( rproc /= myid ) then ! Don't add points already on this proc.
                  mycol = findcolour(iqq,il_g)
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
               ! Which process has this point
               rproc = qproc(iqq)
               if ( rproc /= myid ) then ! Don't add points already on this proc.
                  mycol = findcolour(iqq,il_g)
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
            do n = 0,nproc-1 
               bnds(n)%rlen_bg(min(icol+1,maxcolour)) = bnds(n)%rlen_fn(icol) + 1
               bnds(n)%rlen_fn(min(icol+1,maxcolour)) = bnds(n)%rlen_fn(icol)
            end do   
         end if
         
      end do ! icol=1,maxcolour
      
      do n = 0,nproc-1
         bnds(n)%rlenx_bg(1) = bnds(n)%rlen_fn(maxcolour) + 1
         bnds(n)%rlenx_fn(1) = bnds(n)%rlen_fn(maxcolour)
      end do   
      
      do icol = 1,maxcolour
      
         do n = 1,npan
            ! NE, EN
            iql = indp(ipan,jpan,n)
            iqg = indg(ipan,jpan,n)
            iqq = ine_g(iqg)
            ! Which process has this point
            rproc = qproc(iqq)
            if ( rproc /= myid ) then ! Don't add points already on this proc.
               mycol = findcolour(iqq,il_g)
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
            iqq = ien_g(iqg)
            ! Which process has this point
            rproc = qproc(iqq)
            if ( rproc /= myid ) then ! Don't add points already on this proc.
               mycol = findcolour(iqq,il_g)
               if ( mycol == icol ) then
                  if ( ien_g(iqg) == ine_g(iqg) ) then
                     ! Avoid duplicate 
                     ien(iql) = ine(iql)
                  else
                     ! Add this point to request list
                     bnds(rproc)%rlenx_fn(icol) = bnds(rproc)%rlenx_fn(icol) + 1
                     call check_bnds_alloc(rproc, iext)
                     bnds(rproc)%request_list(bnds(rproc)%rlenx_fn(icol)) = iqq
                     ! Increment extended region index
                     iext = iext + 1
                     bnds(rproc)%unpack_list(bnds(rproc)%rlenx_fn(icol)) = iext
                     ien(iql) = ifull + iext
                  end if   
               end if
            end if
         end do

         do n = 1,npan
            ! SE, ES
            iql = indp(ipan,1,n)
            iqg = indg(ipan,1,n)
            iqq = ise_g(iqg)
            ! Which process has this point
            rproc = qproc(iqq)
            if ( rproc /= myid ) then ! Don't add points already on this proc.
               mycol = findcolour(iqq,il_g)
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
            iqq = ies_g(iqg)
            ! Which process has this point
            rproc = qproc(iqq)
            if ( rproc /= myid ) then ! Don't add points already on this proc.
               mycol = findcolour(iqq,il_g)
               if ( mycol == icol ) then
                  if ( ies_g(iqg) == ise_g(iqg) ) then
                     ! Avoid duplicate
                     ies(iql) = ise(iql)
                  else
                     ! Add this point to request list    
                     bnds(rproc)%rlenx_fn(icol) = bnds(rproc)%rlenx_fn(icol) + 1
                     call check_bnds_alloc(rproc, iext)
                     bnds(rproc)%request_list(bnds(rproc)%rlenx_fn(icol)) = iqq
                     ! Increment extended region index
                     iext = iext + 1
                     bnds(rproc)%unpack_list(bnds(rproc)%rlenx_fn(icol)) = iext
                     ies(iql) = ifull + iext
                  end if   
               end if
            end if
         end do
         
         do n = 1,npan
            ! NW, WN
            iql = indp(1,jpan,n)
            iqg = indg(1,jpan,n)
            iqq = inw_g(iqg)
            ! Which process has this point
            rproc = qproc(iqq)
            if ( rproc /= myid ) then ! Don't add points already on this proc.
               mycol = findcolour(iqq,il_g)
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
            iqq = iwn_g(iqg)
            ! Which process has this point
            rproc = qproc(iqq)
            if ( rproc /= myid ) then ! Don't add points already on this proc.
               mycol = findcolour(iqq,il_g)
               if ( mycol == icol ) then
                  if ( iwn_g(iqg) == inw_g(iqg) ) then
                     ! Avoid duplicate 
                     iwn(iql) = inw(iql)
                  else
                     ! Add this point to request list
                     bnds(rproc)%rlenx_fn(icol) = bnds(rproc)%rlenx_fn(icol) + 1
                     call check_bnds_alloc(rproc, iext)
                     bnds(rproc)%request_list(bnds(rproc)%rlenx_fn(icol)) = iqq
                     ! Increment extended region index
                     iext = iext + 1
                     bnds(rproc)%unpack_list(bnds(rproc)%rlenx_fn(icol)) = iext
                     iwn(iql) = ifull + iext
                  end if   
               end if
            end if
         end do

         do n = 1,npan
            ! SW, WS
            iql = indp(1,1,n)
            iqg = indg(1,1,n)
            iqq = isw_g(iqg)
            ! Which process has this point
            rproc = qproc(iqq)
            if ( rproc /= myid ) then ! Don't add points already on this proc.
               mycol = findcolour(iqq,il_g)
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
            iqq = iws_g(iqg)
            ! Which process has this point
            rproc = qproc(iqq)
            if ( rproc /= myid ) then ! Don't add points already on this proc.
               mycol = findcolour(iqq,il_g)
               if ( mycol == icol ) then
                  if ( iws_g(iqg) == isw_g(iqg) ) then
                     ! Avoid duplicate 
                     iws(iql) = isw(iql)
                  else
                     ! Add this point to request list
                     bnds(rproc)%rlenx_fn(icol) = bnds(rproc)%rlenx_fn(icol) + 1
                     call check_bnds_alloc(rproc, iext)
                     bnds(rproc)%request_list(bnds(rproc)%rlenx_fn(icol)) = iqq
                     ! Increment extended region index
                     iext = iext + 1
                     bnds(rproc)%unpack_list(bnds(rproc)%rlenx_fn(icol)) = iext
                     iws(iql) = ifull + iext
                  end if   
               end if
            end if
         end do
         
         if ( icol < maxcolour ) then
            do n = 0,nproc-1 
               bnds(n)%rlenx_bg(min(icol+1,maxcolour)) = bnds(n)%rlenx_fn(icol) + 1
               bnds(n)%rlenx_fn(min(icol+1,maxcolour)) = bnds(n)%rlenx_fn(icol)
            end do   
         end if
         
      end do ! icol=1,maxcolour

      ! Now set up the second row
      do n = 0,nproc-1
         bnds(n)%rlen2 = bnds(n)%rlenx_fn(maxcolour)  ! so that they're appended.
      end do   

      do n = 1,npan

         !     Start with W edge
         i = 1
         do j = 1,jpan
            iq = indg(i,j,n)
            iqq = iww_g(iq)
            ! Which process has this point
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
            ! Which process has this point
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
            ! Which process has this point
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
            ! Which process has this point
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
      if ( myid==0 ) then
         write(6,*) "-> Construct list of neighbours"
      end if
      allocate( neigharray_g(0:nproc-1) )
      ! default neighbour list
      do n = 0,nproc-1
         neigharray_g(n) = bnds(n)%rlen2 > 0
      end do   
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
      neighnum = count( neigharray_g(:) )
      do n = 0,nproc-1
         if ( neigharray_g(n) ) then
            bnds(n)%len = max( bnds(n)%len, maxbuflen*maxvertlen )
         end if
      end do

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
      allocate( dums(3*maxcolour+1,neighnum), dumr(3*maxcolour+1,neighnum) )

      
      ! Communicate lengths for rlenh, rlen, rlenx and rlen2
      if ( myid==0 ) then
         write(6,*) "-> Communicate halo message length for scalars"
      end if
      lcomm = comm_world
      nreq = 0
      do iproc = 1,neighnum
         rproc = neighlist(iproc) ! Recv from
         if ( bnds(rproc)%rlen2 > 0 ) then
            nreq = nreq + 1
            lproc = rproc
#ifdef i8r8
            call MPI_IRecv( dumr(:,iproc), int(3*maxcolour+1,4), MPI_INTEGER8, lproc, &
                            itag, lcomm, ireq(nreq), ierr )
#else
            call MPI_IRecv( dumr(:,iproc), int(3*maxcolour+1,4), MPI_INTEGER, lproc,  &
                            itag, lcomm, ireq(nreq), ierr )
#endif
         end if
      end do
      do iproc = neighnum,1,-1
         sproc = neighlist(iproc)  ! Send to
         if ( bnds(sproc)%rlen2 > 0 ) then
            nreq = nreq + 1
            dums(1,iproc) = bnds(sproc)%rlen2
            dums(2:maxcolour+1,iproc) = bnds(sproc)%rlenh_fn(1:maxcolour)
            dums(maxcolour+2:2*maxcolour+1,iproc) = bnds(sproc)%rlen_fn(1:maxcolour)
            dums(2*maxcolour+2:3*maxcolour+1,iproc) = bnds(sproc)%rlenx_fn(1:maxcolour)
            lproc = sproc
#ifdef i8r8
            call MPI_ISend( dums(:,iproc), int(3*maxcolour+1,4), MPI_INTEGER8, lproc, &
                            itag, lcomm, ireq(nreq), ierr )
#else
            call MPI_ISend( dums(:,iproc), int(3*maxcolour+1,4), MPI_INTEGER, lproc,  &
                            itag, lcomm, ireq(nreq), ierr )
#endif
         end if
      end do
      call MPI_Waitall( nreq, ireq, MPI_STATUSES_IGNORE, ierr )
      do iproc = 1,neighnum
         rproc = neighlist(iproc)
         if ( bnds(rproc)%rlen2 > 0 ) then
            bnds(rproc)%slen2 = dumr(1,iproc)
            bnds(rproc)%slenh_fn(1:maxcolour) = dumr(2:maxcolour+1,iproc)
            bnds(rproc)%slen_fn(1:maxcolour)  = dumr(maxcolour+2:2*maxcolour+1,iproc)
            bnds(rproc)%slenx_fn(1:maxcolour) = dumr(2*maxcolour+2:3*maxcolour+1,iproc)
         end if
      end do
      do n = 0,nproc-1
         bnds(n)%slenh_bg(1) = 1
         do i = 2,maxcolour
            bnds(n)%slenh_bg(i) = bnds(n)%slenh_fn(i-1) + 1
         end do
         bnds(n)%slen_bg(1)  = bnds(n)%slenh_fn(maxcolour) + 1
         do i = 2,maxcolour
            bnds(n)%slen_bg(i)  = bnds(n)%slen_fn(i-1) + 1
         end do
         bnds(n)%slenx_bg(1) = bnds(n)%slen_fn(maxcolour) + 1
         do i = 2,maxcolour
            bnds(n)%slenx_bg(i) = bnds(n)%slenx_fn(i-1) + 1
         end do
      end do
      
      ! define index=0 for all coloured grid points
      do n = 0,nproc-1
         bnds(n)%rlenh_bg(0) = bnds(n)%rlenh_bg(1)
         bnds(n)%rlenh_fn(0) = bnds(n)%rlenh_fn(maxcolour)
         bnds(n)%slenh_bg(0) = bnds(n)%slenh_bg(1)
         bnds(n)%slenh_fn(0) = bnds(n)%slenh_fn(maxcolour)
         bnds(n)%rlen_bg(0) = bnds(n)%rlen_bg(1)
         bnds(n)%rlen_fn(0) = bnds(n)%rlen_fn(maxcolour)
         bnds(n)%slen_bg(0) = bnds(n)%slen_bg(1)
         bnds(n)%slen_fn(0) = bnds(n)%slen_fn(maxcolour)
         bnds(n)%rlenx_bg(0) = bnds(n)%rlenx_bg(1)
         bnds(n)%rlenx_fn(0) = bnds(n)%rlenx_fn(maxcolour)
         bnds(n)%slenx_bg(0) = bnds(n)%slenx_bg(1)
         bnds(n)%slenx_fn(0) = bnds(n)%slenx_fn(maxcolour)
      end do
      
      deallocate( dums, dumr )

      
      ! Now, for each process send the list of points I want.
      ! The state of being a neighbour is reflexive so only expect to
      ! recv from those processes I send to (are there grid arrangements for
      ! which this would not be true?)
      if ( myid==0 ) then
         write(6,*) "-> Communicate halo request list for scalars"
      end if
      lcomm = comm_world
      nreq = 0
      do iproc = 1,neighnum
         rproc = neighlist(iproc)  ! Recv from
         llen = bnds(rproc)%slen2
         if ( llen > 0 ) then
            allocate( bnds(rproc)%send_list(llen) ) 
            nreq = nreq + 1         
            lproc = rproc
#ifdef i8r8
            call MPI_IRecv( bnds(rproc)%send_list, llen, MPI_INTEGER8, lproc, &
                            itag, lcomm, ireq(nreq), ierr )
#else
            call MPI_IRecv( bnds(rproc)%send_list, llen, MPI_INTEGER, lproc,  &
                            itag, lcomm, ireq(nreq), ierr )
#endif
         end if
      end do
      do iproc = 1,neighnum
         sproc = neighlist(iproc)  ! Send to
         ! Send list of requests
         llen = bnds(sproc)%rlen2
         if ( llen > 0 ) then
            lproc = sproc
            nreq = nreq + 1
#ifdef i8r8            
            call MPI_ISend( bnds(sproc)%request_list, llen, MPI_INTEGER8, &
                            lproc, itag, lcomm, ireq(nreq), ierr )
#else
            call MPI_ISend( bnds(sproc)%request_list, llen, MPI_INTEGER,  &
                            lproc, itag, lcomm, ireq(nreq), ierr )
#endif
         end if
      end do
      call MPI_Waitall( nreq, ireq, MPI_STATUSES_IGNORE, ierr )

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
!     other processes. In the 1D code values on the same process are
!     copied only if they have to be swapped.
!     This only makes a difference on 1, 2 or 3 processes.

      if ( myid==0 ) then
         write(6,*) "-> Calculate neighbour direction index for vectors"
      end if
      iextu = 0
      iextv = 0

      ! save start of isv indices
      do n = 0,nproc-1
         bnds(n)%rlen_sv_bg = 1
      end do   

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
            ! Which process has this point
            rproc = qproc(iqq)
            ! Only need to add to bounds region if it's on another process
            ! or if it's on this process and needs to be swapped.
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
      do n = 0,nproc-1
         bnds(n)%rlen_nv_bg = bnds(n)%rlen_wu_fn + 1
         bnds(n)%rlen_eu_fn = bnds(n)%rlen_wu_fn
      end do   

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

      ! save start of issv indices
      do n = 0,nproc-1
         bnds(n)%rlen_ssv_bg = bnds(n)%rlen_eu_fn + 1
         bnds(n)%rlen_wwu_fn = bnds(n)%rlen_eu_fn
      end do   

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
            ! Which process has this point
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
      do n = 0,nproc-1
         bnds(n)%rlen_nnv_bg = bnds(n)%rlen_wwu_fn + 1
         bnds(n)%rlen_eeu_fn = bnds(n)%rlen_wwu_fn
      end do   

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
      do n = 0,nproc-1
         bnds(n)%rlen_su_bg = bnds(n)%rlen_eeu_fn + 1
         bnds(n)%rlen_ev_fn = bnds(n)%rlen_eeu_fn
      end do   

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
            ! Which process has this point
            rproc = qproc(iqq)
            ! Only need to add to bounds region if it's on another process
            ! or if it's on this process and needs to be swapped.
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
      
      if ( iextu > iextra ) then
         write(6,*) "IEXTU too large", iextu, iextra
         call ccmpi_abort(-1)
      end if

      if ( iextv > iextra ) then
         write(6,*) "IEXTV too large", iextv, iextra
         call ccmpi_abort(-1)
      end if

      
      allocate( dums(5,neighnum), dumr(5,neighnum) )
      
      ! Communicate lengths for rlen_uv and rlenx_uv, etc
      if ( myid==0 ) then
         write(6,*) "-> Communicate halo message length for vectors"
      end if
      lcomm = comm_world
      nreq = 0
      do iproc = 1,neighnum
         rproc = neighlist(iproc)  ! Recv from
         if ( bnds(rproc)%rlen_ev_fn > 0 ) then
            nreq = nreq + 1
            lproc = rproc
#ifdef i8r8
            call MPI_IRecv( dumr(:,iproc), 5_4, MPI_INTEGER8, lproc, &
                 itag, lcomm, ireq(nreq), ierr )
#else
            call MPI_IRecv( dumr(:,iproc), 5_4, MPI_INTEGER, lproc,  &
                 itag, lcomm, ireq(nreq), ierr )
#endif
         end if
      end do
      do iproc = neighnum,1,-1
         sproc = neighlist(iproc)  ! Send to
         if ( bnds(sproc)%rlen_ev_fn > 0 ) then
            nreq = nreq + 1
            dums(1,iproc) = bnds(sproc)%rlen_wu_fn
            dums(2,iproc) = bnds(sproc)%rlen_eu_fn
            dums(3,iproc) = bnds(sproc)%rlen_wwu_fn
            dums(4,iproc) = bnds(sproc)%rlen_eeu_fn
            dums(5,iproc) = bnds(sproc)%rlen_ev_fn
            lproc = sproc
#ifdef i8r8
            call MPI_ISend( dums(:,iproc), 5_4, MPI_INTEGER8, lproc, &
                 itag, lcomm, ireq(nreq), ierr )
#else
            call MPI_ISend( dums(:,iproc), 5_4, MPI_INTEGER, lproc,  &
                 itag, lcomm, ireq(nreq), ierr )
#endif
         end if
      end do
      call MPI_Waitall( nreq, ireq, MPI_STATUSES_IGNORE, ierr )
      do iproc = 1,neighnum
         rproc = neighlist(iproc)
         if ( bnds(rproc)%rlen_ev_fn > 0 ) then
            bnds(rproc)%slen_wu_fn  = dumr(1,iproc)
            bnds(rproc)%slen_eu_fn  = dumr(2,iproc)
            bnds(rproc)%slen_wwu_fn = dumr(3,iproc)
            bnds(rproc)%slen_eeu_fn = dumr(4,iproc)
            bnds(rproc)%slen_ev_fn  = dumr(5,iproc)
         end if
      end do
      do n = 0,nproc-1
         bnds(n)%slen_sv_bg  = 1
         bnds(n)%slen_nv_bg  = bnds(n)%slen_wu_fn  + 1
         bnds(n)%slen_ssv_bg = bnds(n)%slen_eu_fn  + 1
         bnds(n)%slen_nnv_bg = bnds(n)%slen_wwu_fn + 1
         bnds(n)%slen_su_bg  = bnds(n)%slen_eeu_fn + 1
      end do   
      
      deallocate( dumr, dums )

      
      ! Now, for each process send the list of points I want.
      ! Also have to send the swap list
      if ( myid==0 ) then
         write(6,*) "-> Communicate halo request list for vectors"
      end if
      lcomm = comm_world
      nreq = 0
      do iproc = 1,neighnum
         rproc = neighlist(iproc)
         llen = bnds(rproc)%slen_ev_fn
         if ( llen > 0 ) then
            allocate( bnds(rproc)%send_list_uv(llen) ) 
            nreq = nreq + 1
            lproc = rproc
#ifdef i8r8
            call MPI_IRecv( bnds(rproc)%send_list_uv, llen, MPI_INTEGER8, lproc, &
                  itag, lcomm, ireq(nreq), ierr )
#else
            call MPI_IRecv( bnds(rproc)%send_list_uv, llen, MPI_INTEGER, lproc,  &
                  itag, lcomm, ireq(nreq), ierr )
#endif
         end if
      end do
      do iproc = neighnum,1,-1
         sproc = neighlist(iproc)
         if ( bnds(sproc)%rlen_ev_fn > 0 ) then
            ! Send list of requests
            nreq = nreq + 1
            llen = bnds(sproc)%rlen_ev_fn
            lproc = sproc
#ifdef i8r8
            call MPI_ISend( bnds(sproc)%request_list_uv, llen, MPI_INTEGER8, lproc, &
                 itag, lcomm, ireq(nreq), ierr )
#else
            call MPI_ISend( bnds(sproc)%request_list_uv, llen, MPI_INTEGER, lproc,  &
                 itag, lcomm, ireq(nreq), ierr )
#endif
         end if
      end do
      call MPI_Waitall( nreq, ireq, MPI_STATUSES_IGNORE, ierr )

      ! Deallocate arrays that are no longer needed
      ! Note that the request_list for myid is still allocated
      allocate( dumi(maxbuflen) )
      do iproc = 1,neighnum
         rproc = neighlist(iproc)
         if ( bnds(rproc)%rlen_ev_fn > 0 ) then
            deallocate( bnds(rproc)%request_list_uv ) 
            dumi(1:bnds(rproc)%rlen_ev_fn) = bnds(rproc)%unpack_list_uv(1:bnds(rproc)%rlen_ev_fn)
            deallocate( bnds(rproc)%unpack_list_uv )
            allocate( bnds(rproc)%unpack_list_uv(bnds(rproc)%rlen_ev_fn) )
            bnds(rproc)%unpack_list_uv(1:bnds(rproc)%rlen_ev_fn) = dumi(1:bnds(rproc)%rlen_ev_fn)
         end if
      end do
      deallocate( dumi )

      
      allocate( dumsl(maxbuflen,neighnum), dumrl(maxbuflen,neighnum) )
      
      ! Only send the swap list once
      if ( myid==0 ) then
         write(6,*) "-> Communicate halo swap list for vectors"
      end if
      lcomm = comm_world
      nreq = 0
      do iproc = 1,neighnum
         rproc = neighlist(iproc) ! Recv from
         llen = bnds(rproc)%slen_ev_fn
         if ( llen > 0 ) then
            nreq = nreq + 1
            lproc = rproc
            call MPI_IRecv( dumrl(:,iproc), llen, MPI_LOGICAL, &
                  lproc, itag, lcomm, ireq(nreq), ierr )
         end if
      end do
      do iproc = neighnum,1,-1
         sproc = neighlist(iproc) ! Send to
         llen = bnds(sproc)%rlen_ev_fn
         if ( llen > 0 ) then
            ! Send list of requests
            nreq = nreq + 1
            lproc = sproc
            dumsl(1:bnds(sproc)%rlen_ev_fn,iproc) = bnds(sproc)%uv_swap(1:bnds(sproc)%rlen_ev_fn)
            call MPI_ISend( dumsl(:,iproc), llen, MPI_LOGICAL, &
                  lproc, itag, lcomm, ireq(nreq), ierr )
         end if
      end do
      call MPI_Waitall( nreq, ireq, MPI_STATUSES_IGNORE, ierr )
     
      do iproc = 1,neighnum
         rproc = neighlist(iproc)
         if ( bnds(rproc)%slen_ev_fn > 0 ) then
            allocate( bnds(rproc)%send_swap(bnds(rproc)%slen_ev_fn) ) 
            bnds(rproc)%send_swap(1:bnds(rproc)%slen_ev_fn) = dumrl(1:bnds(rproc)%slen_ev_fn,iproc)
         end if
         if ( bnds(rproc)%rlen_ev_fn > 0 ) then
            deallocate( bnds(rproc)%uv_swap )
         end if    
      end do
      
      deallocate( dumrl, dumsl )
      

      ! Only send the neg list once
      if ( myid==0 ) then
         write(6,*) "-> Communicate halo neg list for vectors"
      end if
      lcomm = comm_world
      nreq = 0
      do iproc = 1,neighnum
         rproc = neighlist(iproc)
         llen = bnds(rproc)%slen_ev_fn
         if ( llen > 0 ) then
            allocate( bnds(rproc)%send_neg(llen) ) 
            nreq = nreq + 1
            lproc = rproc
#ifdef i8r8
            call MPI_IRecv( bnds(rproc)%send_neg, llen, MPI_DOUBLE_PRECISION, &
                  lproc, itag, lcomm, ireq(nreq), ierr )
#else
            call MPI_IRecv( bnds(rproc)%send_neg, llen, MPI_REAL, &
                  lproc, itag, lcomm, ireq(nreq), ierr )
#endif
         end if
      end do
      do iproc = neighnum,1,-1
         sproc = neighlist(iproc)
         llen = bnds(sproc)%rlen_ev_fn
         if ( llen > 0 ) then
            ! Send list of requests
            nreq = nreq + 1
            lproc = sproc
#ifdef i8r8
            call MPI_ISend( bnds(sproc)%uv_neg, llen, MPI_DOUBLE_PRECISION, &
                  lproc, itag, lcomm, ireq(nreq), ierr )
#else
            call MPI_ISend( bnds(sproc)%uv_neg, llen, MPI_REAL, &
                  lproc, itag, lcomm, ireq(nreq), ierr )
#endif
         end if
      end do
      call MPI_Waitall( nreq, ireq, MPI_STATUSES_IGNORE, ierr )
      

      do iproc = 1,neighnum
         rproc = neighlist(iproc)
         if ( bnds(rproc)%rlen_ev_fn > 0 ) then
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
         end do
         do i = 1,ipan
            issv(indp(i,2,n)) = isv(indp(i,1,n))
            innv(indp(i,jpan-1,n)) = inv(indp(i,jpan,n))
         end do
      end do

      !  At the moment send_lists use global indices. Convert these to local
      if ( myid==0 ) then
         write(6,*) "-> Convert send list from global to local index"
      end if
      do iproc = neighnum,1,-1
         sproc = neighlist(iproc)  ! Send to
         do iq = 1,bnds(sproc)%slen2
            ! send_list(iq) is global point index, i, j, n are local
            iqq = bnds(sproc)%send_list(iq)
            bnds(sproc)%send_list(iq) = iqg2iq(iqq)
         end do
         do iq = 1,bnds(sproc)%slen_ev_fn
            ! send_list(iq) is global point index, i, j, n are local
            ! Use abs because sign is used as u/v flag
            iqq = abs(bnds(sproc)%send_list_uv(iq))
            bnds(sproc)%send_list_uv(iq) = sign(iqg2iq(iqq),bnds(sproc)%send_list_uv(iq))
         end do
      end do
      if ( bnds(myid)%rlen2 /= 0 ) then
         write(6,*) "ERROR: Invalid rlen2 for myid"
         call ccmpi_abort(-1)
      end if
      do iq = 1,bnds(myid)%rlen_ev_fn
         iqq = abs(bnds(myid)%request_list_uv(iq))
         bnds(myid)%request_list_uv(iq) = sign(iqg2iq(iqq),bnds(myid)%request_list_uv(iq))
      end do

      
      ! Allocate buffer arrays
      if ( myid==0 ) then
         write(6,*) "-> Allocate buffers"
      end if
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
      if ( myid==0 ) then
         write(6,*) "-> Verify halo"
      end if
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
      integer :: testlen, i

!     Allocate the components of the bnds array. It's too much work to
!     get the exact sizes, so allocate a fixed size for each case where
!     there's an interaction.
      if ( bnds(rproc)%len == 0 ) then
         ! Not allocated yet.
         if ( rproc /= myid ) then
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
         testlen = 0 
         do i = 1,maxcolour 
            testlen = max( bnds(rproc)%rlenh_fn(i), bnds(rproc)%rlen_fn(i), &
                           bnds(rproc)%rlenx_fn(i), bnds(rproc)%rlen2,      &
                           testlen )
         end do   
         testlen = testlen*maxvertlen   
         if ( testlen >  bnds(rproc)%len ) then
            write(6,*) "Error, maximum length error in check_bnds_alloc"
            write(6,*) myid, rproc, testlen,  bnds(rproc)%len, maxvertlen
            write(6,*) bnds(rproc)%rlenh_fn(1:maxcolour)
            write(6,*) bnds(rproc)%rlen_fn(1:maxcolour)
            write(6,*) bnds(rproc)%rlenx_fn(1:maxcolour)
            write(6,*) bnds(rproc)%rlen2
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
      real(kind=8) :: disarray_g
      logical, dimension(0:nproc-1), intent(inout) :: neigharray_g
   
      do joff_g = 0,jl_g/jpan-1
         do ioff_g = 0,il_g/ipan-1
            j_g = 1
            do i_g = 1,ipan
               iqg = i_g + ioff_g*ipan + (j_g + joff_g*jpan - 1)*il_g
               disarray_g = x_g(iqq)*x_g(iqg) + y_g(iqq)*y_g(iqg) + z_g(iqq)*z_g(iqg)
               disarray_g = acos( max( min( disarray_g, 1._8 ), -1._8 ) )
               if ( disarray_g < real(maxdis,8) ) then
                  iproc = qproc(iqg)
                  neigharray_g(iproc) = .true.
               end if
            end do
            j_g = jpan
            do i_g = 1,ipan
               iqg = i_g + ioff_g*ipan + (j_g + joff_g*jpan - 1)*il_g
               disarray_g = x_g(iqq)*x_g(iqg) + y_g(iqq)*y_g(iqg) + z_g(iqq)*z_g(iqg)
               disarray_g = acos( max( min( disarray_g, 1._8 ), -1._8 ) )
               if ( disarray_g < real(maxdis,8) ) then
                  iproc = qproc(iqg)
                  neigharray_g(iproc) = .true.
               end if
            end do
            i_g = 1
            do j_g = 1,jpan
               iqg = i_g + ioff_g*ipan + (j_g + joff_g*jpan - 1)*il_g
               disarray_g = x_g(iqq)*x_g(iqg) + y_g(iqq)*y_g(iqg) + z_g(iqq)*z_g(iqg)
               disarray_g = acos( max( min( disarray_g, 1._8 ), -1._8 ) )
               if ( disarray_g < real(maxdis,8) ) then
                  iproc = qproc(iqg)
                  neigharray_g(iproc) = .true.
               end if
            end do
            i_g = ipan
            do j_g = 1,jpan
               iqg = i_g + ioff_g*ipan + (j_g + joff_g*jpan - 1)*il_g
               disarray_g = x_g(iqq)*x_g(iqg) + y_g(iqq)*y_g(iqg) + z_g(iqq)*z_g(iqg)
               disarray_g = acos( max( min( disarray_g, 1._8 ), -1._8 ) )
               if ( disarray_g < real(maxdis,8) ) then
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
      real, dimension(ifull+iextra), intent(inout) :: t
      real, dimension(size(t,1),1,1) :: t_l
      integer, intent(in), optional :: nrows
      logical, intent(in), optional :: corner, nehalf
      integer :: nrows_l
      logical :: corner_l, nehalf_l
      
      nrows_l = 1
      if ( present(nrows) ) then
         nrows_l = nrows 
      end if
      corner_l = .false.
      if ( present(corner) ) then
         corner_l = corner
      end if
      nehalf_l = .false.
      if ( present(nehalf) ) then
         nehalf_l = nehalf 
      end if

      ! colour=0 is for all grid points
      t_l(:,1,1) = t(:)
      call bounds_colour_send4( t_l, 0, nrows=nrows_l, corner=corner_l, nehalf=nehalf_l )
      call bounds_colour_recv4( t_l, 0, nrows=nrows_l, corner=corner_l, nehalf=nehalf_l )
      t(:) = t_l(:,1,1)      

   end subroutine bounds2
   
   subroutine bounds2r8(t, nrows, corner, nehalf)
      real(kind=8), dimension(ifull+iextra), intent(inout) :: t
      real(kind=8), dimension(size(t,1),1,1) :: t_l
      integer, intent(in), optional :: nrows
      logical, intent(in), optional :: corner, nehalf
      integer :: nrows_l
      logical :: corner_l, nehalf_l
      
      nrows_l = 1
      if ( present(nrows) ) then
         nrows_l = nrows 
      end if
      corner_l = .false.
      if ( present(corner) ) then
         corner_l = corner
      end if
      nehalf_l = .false.
      if ( present(nehalf) ) then
         nehalf_l = nehalf 
      end if

      t_l(:,1,1) = t(:)
      call bounds4r8( t_l, nrows=nrows_l, corner=corner_l, nehalf=nehalf_l )
      t(:) = t_l(:,1,1)      

   end subroutine bounds2r8   

   subroutine bounds3(t, nrows, klim, corner, nehalf)
      real, dimension(:,:), intent(inout) :: t
      integer, intent(in), optional :: nrows, klim
      logical, intent(in), optional :: corner, nehalf
      real, dimension(size(t,1),size(t,2),1) :: t_l
      integer :: nrows_l, klim_l
      logical :: corner_l, nehalf_l
      
      klim_l = size(t,2)
      if ( present(klim) ) then
         klim_l = klim
      end if
      nrows_l = 1
      if ( present(nrows) ) then
         nrows_l = nrows 
      end if
      corner_l = .false.
      if ( present(corner) ) then
         corner_l = corner
      end if
      nehalf_l = .false.
      if ( present(nehalf) ) then
         nehalf_l = nehalf 
      end if

      ! colour=0 is for all grid points
      t_l(:,:,1) = t(:,:)
      call bounds_colour_send4( t_l, 0, nrows=nrows_l, klim=klim_l, corner=corner_l, nehalf=nehalf_l )
      call bounds_colour_recv4( t_l, 0, nrows=nrows_l, klim=klim_l, corner=corner_l, nehalf=nehalf_l )
      t(:,:) = t_l(:,:,1)

   end subroutine bounds3

   subroutine bounds3r8(t, nrows, klim, corner, nehalf)
      real(kind=8), dimension(:,:), intent(inout) :: t
      real(kind=8), dimension(size(t,1),size(t,2),1) :: t_l
      integer, intent(in), optional :: nrows, klim
      logical, intent(in), optional :: corner, nehalf
      integer :: nrows_l, klim_l
      logical :: corner_l, nehalf_l
      
      klim_l = size(t,2)
      if ( present(klim) ) then
         klim_l = klim
      end if
      nrows_l = 1
      if ( present(nrows) ) then
         nrows_l = nrows 
      end if
      corner_l = .false.
      if ( present(corner) ) then
         corner_l = corner
      end if
      nehalf_l = .false.
      if ( present(nehalf) ) then
         nehalf_l = nehalf 
      end if

      t_l(:,:,1) = t(:,:)
      call bounds4r8( t_l, nrows=nrows_l, klim=klim_l, corner=corner_l, nehalf=nehalf_l )
      t(:,:) = t_l(:,:,1)

   end subroutine bounds3r8
   
   subroutine bounds4(t, nrows, klim, corner, nehalf)
      real, dimension(:,:,:), intent(inout) :: t
      integer, intent(in), optional :: nrows, klim
      logical, intent(in), optional :: corner, nehalf
      integer :: nrows_l, klim_l
      integer :: istart, iend
      logical :: corner_l, nehalf_l
      
      klim_l = size(t, 2)
      if ( present(klim) ) then
         klim_l = klim
      end if
      nrows_l = 1
      if ( present(nrows) ) then
         nrows_l = nrows 
      end if
      corner_l = .false.
      if ( present(corner) ) then
         corner_l = corner
      end if
      nehalf_l = .false.
      if ( present(nehalf) ) then
         nehalf_l = nehalf 
      end if

      do istart = 1,size(t,3),nagg
         iend = min( istart + nagg - 1, size(t,3) )
         ! colour=0 should send all grid points
         call bounds_colour_send4( t(:,:,istart:iend), 0, nrows=nrows_l, klim=klim_l, corner=corner_l, nehalf=nehalf_l )
         call bounds_colour_recv4( t(:,:,istart:iend), 0, nrows=nrows_l, klim=klim_l, corner=corner_l, nehalf=nehalf_l )
      end do   
      
   end subroutine bounds4
   
   subroutine bounds4r8(t, nrows, klim, corner, nehalf)
      ! Copy the boundary regions.
      real(kind=8), dimension(:,:,:), intent(inout) :: t
      integer, intent(in), optional :: nrows, klim
      logical, intent(in), optional :: corner, nehalf
      logical :: extra, single, double
      integer :: iproc, kx, send_len, recv_len
      integer :: rcount, jproc, mproc, ntr, iq, k, l, n
      integer :: nstart, nend, ntot
      integer, dimension(neighnum) :: rslen, sslen
      integer(kind=4), save :: itag=3
      integer(kind=4) :: ierr, llen, lproc, ldone, lcomm
      integer(kind=4), dimension(2*neighnum) :: donelist

      kx = size(t, 2)
      ntr = size(t, 3)
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

      ! Split messages into corner and non-corner processes
      if ( double ) then
         do n = 1,neighnum  
            rslen(n) = bnds(neighlist(n))%rlen2
            sslen(n) = bnds(neighlist(n))%slen2
         end do
      else if ( extra ) then
         do n = 1,neighnum     
            rslen(n) = bnds(neighlist(n))%rlenx_fn(maxcolour)
            sslen(n) = bnds(neighlist(n))%slenx_fn(maxcolour)
         end do
      else if ( single ) then
         do n = 1,neighnum     
            rslen(n) = bnds(neighlist(n))%rlen_fn(maxcolour)
            sslen(n) = bnds(neighlist(n))%slen_fn(maxcolour)
         end do
      else
         do n = 1,neighnum 
            rslen(n) = bnds(neighlist(n))%rlenh_fn(maxcolour)
            sslen(n) = bnds(neighlist(n))%slenh_fn(maxcolour)
         end do
      end if

      lcomm = comm_world
      
      do nstart = 1,ntr,nagg
         nend = min(nstart+nagg-1,ntr)
         ntot = nend - nstart + 1
         
         itag = mod(itag + 1, 10000)
      
         ! Set up the buffers to send
         nreq = 0
         do iproc = 1,neighnum
            recv_len = rslen(iproc)
            if ( recv_len > 0 ) then
               lproc = neighlist(iproc)  ! Recv from
               nreq = nreq + 1
               rlist(nreq) = iproc
               llen = recv_len*kx*ntot
               call MPI_IRecv( bnds(lproc)%r8buf, llen, MPI_DOUBLE_PRECISION, lproc, &
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
                     do iq = 1,send_len
                        bnds(lproc)%s8buf(iq+(k-1)*send_len+(l-1)*send_len*kx) = &
                          t(bnds(lproc)%send_list(iq),k,l+nstart-1)
                     end do
                  end do
               end do   
               nreq = nreq + 1
               llen = send_len*kx*ntot
               call MPI_ISend( bnds(lproc)%s8buf, llen, MPI_DOUBLE_PRECISION, lproc, &
                    itag, lcomm, ireq(nreq), ierr )
            end if
         end do

         ! Unpack incomming messages
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
                  lproc = neighlist(iproc)
                  do l = 1,ntot
                     do k=1,kx
                        do iq = 1,rslen(iproc)
                           t(ifull+bnds(lproc)%unpack_list(iq),k,l+nstart-1)              &
                                = bnds(lproc)%r8buf(iq+(k-1)*rslen(iproc)+(l-1)*rslen(iproc)*kx)
                        end do
                     end do   
                  end do
               end if
            end do
         end do

      end do  ! nstart 
      
   end subroutine bounds4r8
   
   subroutine bounds_colour_send3(t, colour, nrows, klim, corner, nehalf)
      real, dimension(:,:), intent(in) :: t
      integer, intent(in) :: colour
      integer, intent(in), optional :: nrows, klim
      logical, intent(in), optional :: corner, nehalf
      real, dimension(size(t,1),size(t,2),1) :: t_l
      integer :: nrows_l, klim_l
      logical :: corner_l, nehalf_l
      
      klim_l = size(t, 2)
      if ( present(klim) ) then
         klim_l = klim
      end if
      nrows_l = 1
      if ( present(nrows) ) then
         nrows_l = nrows 
      end if
      corner_l = .false.
      if ( present(corner) ) then
         corner_l = corner
      end if
      nehalf_l = .false.
      if ( present(nehalf) ) then
         nehalf_l = nehalf 
      end if
      
      t_l(:,:,1) = t(:,:)
      call bounds_colour_send4( t_l, colour, nrows=nrows_l, klim=klim_l, corner=corner_l, nehalf=nehalf_l )

   end subroutine bounds_colour_send3
   
   subroutine bounds_colour_send4(t, colour, nrows, klim, corner, nehalf)
      ! Copy the boundary regions. This version allows supports updating
      ! different gridpoint colours
      real, dimension(:,:,:), intent(in) :: t
      integer, intent(in) :: colour
      integer, intent(in), optional :: nrows, klim
      logical, intent(in), optional :: corner, nehalf
      logical :: extra, single, double
      integer :: iproc, kx, recv_len, iqq, ibeg, iend, iq, k
      integer :: l, ntr
      integer(kind=4), save :: itag=4
      integer(kind=4) :: ierr, llen, lproc, lcomm

      if ( colour<0 .or. colour>maxcolour ) then
         write(6,*) "ERROR: Invalid colour for bounds_colour_send"
         call ccmpi_abort(-1)
      end if   
      
      kx = size(t, 2)
      ntr = size(t, 3)
      single = .true.
      extra = .false.
      double = .false.
      if ( present(klim) ) then
         kx = klim
      end if
      if ( present(nehalf) ) then
        single = .not.nehalf
      end if
      if ( present(corner) ) then
         extra = corner
      end if
      if ( present(nrows) ) then
         double = nrows==2
      end if
      
      if ( double ) then
         extra = .true.
         single = .true.
      else if ( extra ) then
         single = .true.
      end if
      
      itag = mod(itag + 1, 10000)
      
      if ( ntr > nagg ) then
         write(6,*) "ERROR: bounds_colour_send4 can only send nagg tracers"
         write(6,*) "ntr, nagg ",ntr, nagg
         call ccmpi_abort(-1)
      end if

!     Set up the buffers to send and recv
      lcomm = comm_world
      nreq = 0
      do iproc = 1,neighnum
         lproc = neighlist(iproc)  ! Recv from
         recv_len = 0
         ibeg = bnds(lproc)%rlenh_bg(colour)
         iend = bnds(lproc)%rlenh_fn(colour)
         if ( iend >= ibeg ) then
            recv_len = recv_len + (iend-ibeg+1)*kx*ntr
         end if   
         if ( single ) then
            ibeg = bnds(lproc)%rlen_bg(colour)
            iend = bnds(lproc)%rlen_fn(colour)       
            if ( iend >= ibeg ) then
               recv_len = recv_len + (iend-ibeg+1)*kx*ntr
            end if   
         end if
         if ( extra ) then
            ibeg = bnds(lproc)%rlenx_bg(colour)
            iend = bnds(lproc)%rlenx_fn(colour)       
            if ( iend >= ibeg ) then
               recv_len = recv_len + (iend-ibeg+1)*kx*ntr
            end if   
         end if
         if ( double ) then
            ibeg = bnds(lproc)%rlenx_fn(maxcolour) + 1
            iend = bnds(lproc)%rlen2       
            if ( iend >= ibeg ) then
               recv_len = recv_len + (iend-ibeg+1)*kx*ntr
            end if   
         end if
         if ( recv_len > 0 ) then
            nreq = nreq + 1
            rlist(nreq) = iproc
            llen = recv_len
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
      do iproc = neighnum,1,-1
         lproc = neighlist(iproc)  ! Send to
         iqq = 0
         ibeg = bnds(lproc)%slenh_bg(colour)
         iend = bnds(lproc)%slenh_fn(colour)
         if ( iend >= ibeg ) then
            do l = 1,ntr
               do k = 1,kx 
                  do iq = 1,iend-ibeg+1
                     bnds(lproc)%sbuf(iqq+iq+(k-1)*(iend-ibeg+1)+(l-1)*(iend-ibeg+1)*kx)  &
                         = t(bnds(lproc)%send_list(iq+ibeg-1),k,l)
                  end do   
               end do
            end do   
         end if
         iqq = iqq + (iend-ibeg+1)*kx*ntr
         if ( single ) then
            ibeg = bnds(lproc)%slen_bg(colour)
            iend = bnds(lproc)%slen_fn(colour)
            if ( iend >= ibeg ) then
               do l = 1,ntr
                  do k = 1,kx
                     do iq = 1,iend-ibeg+1
                        bnds(lproc)%sbuf(iqq+iq+(k-1)*(iend-ibeg+1)+(l-1)*(iend-ibeg+1)*kx)  &
                            = t(bnds(lproc)%send_list(iq+ibeg-1),k,l)
                     end do   
                  end do
               end do   
            end if   
            iqq = iqq + (iend-ibeg+1)*kx*ntr
         end if   
         if ( extra ) then
            ibeg = bnds(lproc)%slenx_bg(colour)
            iend = bnds(lproc)%slenx_fn(colour)
            if ( iend >= ibeg ) then
               do l = 1,ntr
                  do k = 1,kx 
                     do iq = 1,iend-ibeg+1
                        bnds(lproc)%sbuf(iqq+iq+(k-1)*(iend-ibeg+1)+(l-1)*(iend-ibeg+1)*kx)  &
                            = t(bnds(lproc)%send_list(iq+ibeg-1),k,l)
                     end do   
                  end do
               end do   
            end if   
            iqq = iqq + (iend-ibeg+1)*kx*ntr
         end if
         if ( double ) then
            ibeg = bnds(lproc)%slenx_fn(maxcolour) + 1
            iend = bnds(lproc)%slen2
            if ( iend >= ibeg ) then
               do l = 1,ntr
                  do k = 1,kx 
                     do iq = 1,iend-ibeg+1
                        bnds(lproc)%sbuf(iqq+iq+(k-1)*(iend-ibeg+1)+(l-1)*(iend-ibeg+1)*kx)  &
                            = t(bnds(lproc)%send_list(iq+ibeg-1),k,l)
                     end do   
                  end do
               end do   
            end if   
            iqq = iqq + (iend-ibeg+1)*kx*ntr
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

   end subroutine bounds_colour_send4
   
   subroutine bounds_colour_recv3(t, colour, nrows, klim, corner, nehalf)
      real, dimension(:,:), intent(inout) :: t
      integer, intent(in) :: colour
      integer, intent(in), optional :: nrows, klim
      logical, intent(in), optional :: corner, nehalf
      real, dimension(size(t,1),size(t,2),1) :: t_l
      integer :: nrows_l, klim_l
      logical :: corner_l, nehalf_l
      
      klim_l = size(t, 2)
      if ( present(klim) ) then
         klim_l = klim
      end if
      nrows_l = 1
      if ( present(nrows) ) then
         nrows_l = nrows 
      end if
      corner_l = .false.
      if ( present(corner) ) then
         corner_l = corner
      end if
      nehalf_l = .false.
      if ( present(nehalf) ) then
         nehalf_l = nehalf 
      end if
      
      t_l(:,:,1) = t(:,:)
      call bounds_colour_recv4( t_l, colour, nrows=nrows_l, klim=klim_l, corner=corner_l, nehalf=nehalf_l )
      t(:,:) = t_l(:,:,1)

   end subroutine bounds_colour_recv3

   subroutine bounds_colour_recv4(t, colour, nrows, klim, corner, nehalf)
      ! Copy the boundary regions. This version allows supports updating
      ! different gridpoint colours
      real, dimension(:,:,:), intent(inout) :: t
      integer, intent(in) :: colour
      integer, intent(in), optional :: nrows, klim
      logical, intent(in), optional :: corner, nehalf
      logical :: extra, single, double
      integer :: iproc, kx, iqq, ibeg, iend
      integer :: rcount, jproc, mproc, iq, k
      integer :: l, ntr
      integer(kind=4) :: ierr, lproc, ldone
      integer(kind=4), dimension(2*neighnum) :: donelist
      
      if ( colour<0 .or. colour>maxcolour ) then
         write(6,*) "ERROR: Invalid colour for bounds_colour_recv"
         call ccmpi_abort(-1)
      end if   
      
      kx = size(t, 2)
      ntr = size(t, 3)
      single = .true.
      extra = .false.
      double = .false.
      if ( present(klim) ) then
         kx = klim
      end if
      if ( present(nehalf) ) then
         single = .not.nehalf
      end if
      if ( present(corner) ) then
         extra = corner
      end if
      if ( present(nrows) ) then
         double = nrows==2
      end if
      
      if ( double ) then
         extra = .true.
         single = .true.
      else if ( extra ) then
         single = .true.
      end if
      
      if ( ntr > nagg ) then
         write(6,*) "ERROR: bounds_colour_recv4 can only send nagg tracers"
         write(6,*) "ntr, nagg ",ntr, nagg
         call ccmpi_abort(-1)
      end if
      
      ! Unpack incomming messages
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
               lproc = neighlist(iproc)
               iqq = 0
               ibeg = bnds(lproc)%rlenh_bg(colour)
               iend = bnds(lproc)%rlenh_fn(colour)
               if ( iend >= ibeg ) then
                  do l = 1,ntr
                     do k = 1,kx
                        do iq = 1,iend-ibeg+1
                           t(ifull+bnds(lproc)%unpack_list(iq+ibeg-1),k,l)  &
                               = bnds(lproc)%rbuf(iqq+iq+(k-1)*(iend-ibeg+1)+(l-1)*(iend-ibeg+1)*kx)
                        end do   
                     end do   
                  end do
                  iqq = iqq + (iend-ibeg+1)*kx*ntr
               end if   
               if ( single ) then
                  ibeg = bnds(lproc)%rlen_bg(colour)
                  iend = bnds(lproc)%rlen_fn(colour)
                  if ( iend >= ibeg ) then
                     do l = 1,ntr
                        do k = 1,kx
                           do iq = 1,iend-ibeg+1
                              t(ifull+bnds(lproc)%unpack_list(iq+ibeg-1),k,l)  &
                                  = bnds(lproc)%rbuf(iqq+iq+(k-1)*(iend-ibeg+1)+(l-1)*(iend-ibeg+1)*kx)
                           end do   
                        end do   
                     end do
                     iqq = iqq + (iend-ibeg+1)*kx*ntr
                  end if   
               end if   
               if ( extra ) then
                  ibeg = bnds(lproc)%rlenx_bg(colour)
                  iend = bnds(lproc)%rlenx_fn(colour)
                  if ( iend >= ibeg ) then
                     do l = 1,ntr
                        do k = 1,kx
                           do iq = 1,iend-ibeg+1
                              t(ifull+bnds(lproc)%unpack_list(iq+ibeg-1),k,l)  &
                                  = bnds(lproc)%rbuf(iqq+iq+(k-1)*(iend-ibeg+1)+(l-1)*(iend-ibeg+1)*kx)
                           end do   
                        end do   
                     end do
                     iqq = iqq + (iend-ibeg+1)*kx*ntr
                  end if   
               end if
               if ( double ) then
                  ibeg = bnds(lproc)%rlenx_fn(maxcolour) + 1
                  iend = bnds(lproc)%rlen2
                  if ( iend >= ibeg ) then
                     do l = 1,ntr
                        do k = 1,kx
                           do iq = 1,iend-ibeg+1
                              t(ifull+bnds(lproc)%unpack_list(iq+ibeg-1),k,l)  &
                                  = bnds(lproc)%rbuf(iqq+iq+(k-1)*(iend-ibeg+1)+(l-1)*(iend-ibeg+1)*kx)
                           end do   
                        end do   
                     end do
                     iqq = iqq + (iend-ibeg+1)*kx*ntr
                  end if   
               end if
            end if ! mproc <= rreq  
         end do    ! jproc = 1,ldone
      end do       ! while( rcount > 0 )

   end subroutine bounds_colour_recv4
   
   subroutine bounds_send3(t, nrows, klim, corner, nehalf)
      real, dimension(:,:), intent(in) :: t
      integer, intent(in), optional :: nrows, klim
      logical, intent(in), optional :: corner, nehalf
      real, dimension(size(t,1),size(t,2),1) :: t_l
      integer :: nrows_l, klim_l
      logical :: corner_l, nehalf_l
      
      nrows_l = 1
      if ( present(nrows) ) then
         nrows_l = nrows
      end if
      klim_l = size(t,2)
      if ( present(klim) ) then
         klim_l = klim
      end if
      corner_l = .false.
      if ( present(corner) ) then
         corner_l = corner
      end if
      nehalf_l = .false.
      if ( present(nehalf) ) then
         nehalf_l = nehalf
      end if
      
      t_l(:,:,1) = t(:,:)
      call bounds_colour_send4( t_l, 0, nrows=nrows_l, klim=klim_l, corner=corner_l, nehalf=nehalf_l )

   end subroutine bounds_send3

   subroutine bounds_send4(t, nrows, klim, corner, nehalf)
      real, dimension(:,:,:), intent(in) :: t
      integer, intent(in), optional :: nrows, klim
      logical, intent(in), optional :: corner, nehalf
      integer :: nrows_l, klim_l
      logical :: corner_l, nehalf_l
      
      nrows_l = 1
      if ( present(nrows) ) then
         nrows_l = nrows
      end if
      klim_l = size(t,2)
      if ( present(klim) ) then
         klim_l = klim
      end if
      corner_l = .false.
      if ( present(corner) ) then
         corner_l = corner
      end if
      nehalf_l = .false.
      if ( present(nehalf) ) then
         nehalf_l = nehalf
      end if
      
      ! colour=0 sends all grid points
      call bounds_colour_send4( t, 0, nrows=nrows_l, klim=klim_l, corner=corner_l, nehalf=nehalf_l )

   end subroutine bounds_send4
   
   subroutine bounds_recv3(t, nrows, klim, corner, nehalf)
      real, dimension(:,:), intent(inout) :: t
      integer, intent(in), optional :: nrows, klim
      logical, intent(in), optional :: corner, nehalf
      real, dimension(size(t,1),size(t,2),1) :: t_l
      integer :: nrows_l, klim_l
      logical :: corner_l, nehalf_l
      
      nrows_l = 1
      if ( present(nrows) ) then
         nrows_l = nrows
      end if
      klim_l = size(t,2)
      if ( present(klim) ) then
         klim_l = klim
      end if
      corner_l = .false.
      if ( present(corner) ) then
         corner_l = corner
      end if
      nehalf_l = .false.
      if ( present(nehalf) ) then
         nehalf_l = nehalf
      end if
      
      t_l(:,:,1) = t(:,:)
      call bounds_colour_recv4( t_l, 0, nrows=nrows_l, klim=klim_l, corner=corner_l, nehalf=nehalf_l )
      t(:,:) = t_l(:,:,1)

   end subroutine bounds_recv3

   subroutine bounds_recv4(t, nrows, klim, corner, nehalf)
      real, dimension(:,:,:), intent(inout) :: t
      integer, intent(in), optional :: nrows, klim
      logical, intent(in), optional :: corner, nehalf
      integer :: nrows_l, klim_l
      logical :: corner_l, nehalf_l
      
      nrows_l = 1
      if ( present(nrows) ) then
         nrows_l = nrows
      end if
      klim_l = size(t,2)
      if ( present(klim) ) then
         klim_l = klim
      end if
      corner_l = .false.
      if ( present(corner) ) then
         corner_l = corner
      end if
      nehalf_l = .false.
      if ( present(nehalf) ) then
         nehalf_l = nehalf
      end if
      
      ! colour=0 receives all grid points
      call bounds_colour_recv4( t, 0, nrows=nrows_l, klim=klim_l, corner=corner_l, nehalf=nehalf_l )

   end subroutine bounds_recv4
   
   subroutine boundsuv2(u, v, stag, allvec)
      ! Copy the boundary regions of u and v. This doesn't require the
      ! diagonal points like (0,0), but does have to take care of the
      ! direction changes.
      real, dimension(ifull+iextra), intent(inout) :: u, v
      integer, intent(in), optional :: stag
      logical, intent(in), optional :: allvec
      real, dimension(ifull+iextra,1) :: u_l, v_l
      integer :: stag_l
      logical :: allvec_l
      
      stag_l = 0
      allvec_l = .false.
      if ( present(stag) ) then
         stag_l = stag
      end if
      if ( present(allvec) ) then
         allvec_l = allvec
      end if
      u_l(:,1) = u(:)
      v_l(:,1) = v(:)
      call boundsuv3(u_l, v_l, stag=stag_l, allvec=allvec_l)
      u(:) = u_l(:,1)
      v(:) = v_l(:,1)

   end subroutine boundsuv2

   subroutine boundsuv3(u, v, stag, allvec)
      ! Copy the boundary regions of u and v. This doesn't require the
      ! diagonal points like (0,0), but does have to take care of the
      ! direction changes.
      real, dimension(:,:), intent(inout) :: u, v
      integer, intent(in), optional :: stag
      logical, intent(in), optional :: allvec
      logical :: extra, fsvwu, fnveu, fssvwwu, fnnveeu, fsuev
      integer :: iq, iqz, iproc, kx, rproc, sproc, iqq, recv_len
      integer :: rcount, myrlen, jproc, mproc, stagmode, k, iqlen
      integer(kind=4), save :: itag=6
      integer(kind=4) :: ierr, llen, lproc
      integer(kind=4) :: ldone, lcomm
      integer(kind=4), dimension(2*neighnum) :: donelist  

      kx = size(u, 2)
      extra = .false.
      stagmode = 0
      if ( present(stag) ) then
         stagmode = stag
      end if
      if ( present(allvec) ) then
         extra = allvec
      end if

      if ( extra ) then
         fsvwu = .true.
         fnveu = .true.
         fssvwwu = .false.
         fnnveeu = .false.
         fsuev = .true.
      else if ( stagmode == 1 ) then
         fsvwu = .false.
         fnveu = .true.
         fssvwwu = .false.
         fnnveeu = .true.
         fsuev = .false.
      else if ( stagmode == 2 ) then
         fsvwu = .true.
         fnveu = .true.
         fssvwwu = .false.
         fnnveeu = .true. ! fnnveeu requires fnveu
         fsuev = .false.
      else if ( stagmode == 3 ) then
         fsvwu = .true.
         fnveu = .true.
         fssvwwu = .true. ! fssvwwu requires fsvwu
         fnnveeu = .false.
         fsuev = .false.
      else if ( stagmode == 5 ) then
         fsvwu = .true.
         fnveu = .false.
         fssvwwu = .true.
         fnnveeu = .false.
         fsuev = .false.
      else if ( stagmode == -9 ) then
         fsvwu = .true.
         fnveu = .false.
         fssvwwu = .false.
         fnnveeu = .false.
         fsuev = .false.
      else if ( stagmode == -10 ) then
         fsvwu = .false.
         fnveu = .true.
         fssvwwu = .false.
         fnnveeu = .false.
         fsuev = .false.
      else
         fsvwu = .true.
         fnveu = .true.
         fssvwwu = .false.
         fnnveeu = .false.
         fsuev = .false.
      end if
      myrlen = bnds(myid)%rlen_ev_fn
      
      itag = mod(itag + 1, 10000)

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
         if ( recv_len > 0 ) then 
            nreq = nreq + 1
            rlist(nreq) = iproc
            llen = recv_len
            lproc = rproc
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
      do iproc = neighnum,1,-1
         sproc = neighlist(iproc)  ! Send to
         ! Build up list of points
         iqq = 0
         if ( fsvwu ) then
            iqz = iqq - bnds(sproc)%slen_sv_bg + 1 
            iqlen = bnds(sproc)%slen_wu_fn - bnds(sproc)%slen_sv_bg + 1
            do k = 1,kx
               do iq = bnds(sproc)%slen_sv_bg,bnds(sproc)%slen_wu_fn
                  ! send_list_uv(iq) is point index.
                  ! Use abs because sign is used as u/v flag
                  if ( bnds(sproc)%send_list_uv(iq)>0 .neqv. bnds(sproc)%send_swap(iq) ) then
                     bnds(sproc)%sbuf(iqz+iq+(k-1)*iqlen) = bnds(sproc)%send_neg(iq)*u(abs(bnds(sproc)%send_list_uv(iq)),k)
                  else
                     bnds(sproc)%sbuf(iqz+iq+(k-1)*iqlen) = bnds(sproc)%send_neg(iq)*v(abs(bnds(sproc)%send_list_uv(iq)),k)
                 end if
               end do
            end do   
            iqq = iqq + iqlen*kx
         end if
         if ( fnveu ) then
            iqz = iqq - bnds(sproc)%slen_nv_bg + 1 
            iqlen = bnds(sproc)%slen_eu_fn - bnds(sproc)%slen_nv_bg + 1
            do k = 1,kx 
               do iq = bnds(sproc)%slen_nv_bg,bnds(sproc)%slen_eu_fn
                  ! send_list_uv(iq) is point index.
                  ! Use abs because sign is used as u/v flag
                  if ( bnds(sproc)%send_list_uv(iq)>0 .neqv. bnds(sproc)%send_swap(iq) ) then
                     bnds(sproc)%sbuf(iqz+iq+(k-1)*iqlen) = bnds(sproc)%send_neg(iq)*u(abs(bnds(sproc)%send_list_uv(iq)),k)
                  else
                     bnds(sproc)%sbuf(iqz+iq+(k-1)*iqlen) = bnds(sproc)%send_neg(iq)*v(abs(bnds(sproc)%send_list_uv(iq)),k)
                  end if 
               end do
            end do   
            iqq = iqq + iqlen*kx
         end if
         if ( fssvwwu ) then
            iqz = iqq - bnds(sproc)%slen_ssv_bg + 1 
            iqlen = bnds(sproc)%slen_wwu_fn - bnds(sproc)%slen_ssv_bg + 1
            do k = 1,kx 
               do iq = bnds(sproc)%slen_ssv_bg,bnds(sproc)%slen_wwu_fn
                  ! send_list_uv(iq) is point index.
                  ! Use abs because sign is used as u/v flag
                  if ( bnds(sproc)%send_list_uv(iq)>0 .neqv. bnds(sproc)%send_swap(iq) ) then
                     bnds(sproc)%sbuf(iqz+iq+(k-1)*iqlen) = bnds(sproc)%send_neg(iq)*u(abs(bnds(sproc)%send_list_uv(iq)),k)
                  else
                     bnds(sproc)%sbuf(iqz+iq+(k-1)*iqlen) = bnds(sproc)%send_neg(iq)*v(abs(bnds(sproc)%send_list_uv(iq)),k)
                  end if 
               end do
            end do   
            iqq = iqq + iqlen*kx
         end if
         if ( fnnveeu ) then
            iqz = iqq - bnds(sproc)%slen_nnv_bg + 1 
            iqlen = bnds(sproc)%slen_eeu_fn - bnds(sproc)%slen_nnv_bg + 1
            do k = 1,kx
               do iq = bnds(sproc)%slen_nnv_bg,bnds(sproc)%slen_eeu_fn
                  ! send_list_uv(iq) is point index.
                  ! Use abs because sign is used as u/v flag
                  if ( bnds(sproc)%send_list_uv(iq)>0 .neqv. bnds(sproc)%send_swap(iq) ) then
                     bnds(sproc)%sbuf(iqz+iq+(k-1)*iqlen) = bnds(sproc)%send_neg(iq)*u(abs(bnds(sproc)%send_list_uv(iq)),k)
                  else
                     bnds(sproc)%sbuf(iqz+iq+(k-1)*iqlen) = bnds(sproc)%send_neg(iq)*v(abs(bnds(sproc)%send_list_uv(iq)),k)
                  end if
               end do
            end do   
            iqq = iqq + iqlen*kx
         end if
         if ( fsuev ) then
            iqz = iqq - bnds(sproc)%slen_su_bg + 1 
            iqlen = bnds(sproc)%slen_ev_fn - bnds(sproc)%slen_su_bg + 1
            do k = 1,kx 
               do iq = bnds(sproc)%slen_su_bg,bnds(sproc)%slen_ev_fn
                  ! send_list_uv(iq) is point index.
                  ! Use abs because sign is used as u/v flag
                  if ( bnds(sproc)%send_list_uv(iq)>0 .neqv. bnds(sproc)%send_swap(iq) ) then
                     bnds(sproc)%sbuf(iqz+iq+(k-1)*iqlen) = bnds(sproc)%send_neg(iq)*u(abs(bnds(sproc)%send_list_uv(iq)),k)
                  else
                     bnds(sproc)%sbuf(iqz+iq+(k-1)*iqlen) = bnds(sproc)%send_neg(iq)*v(abs(bnds(sproc)%send_list_uv(iq)),k)
                 end if
              end do
            end do  
            iqq = iqq + iqlen*kx
         end if
         if ( iqq > 0 ) then
            nreq = nreq + 1
            llen = iqq
            lproc = sproc
#ifdef i8r8
            call MPI_ISend( bnds(lproc)%sbuf, llen, MPI_DOUBLE_PRECISION, lproc, &
                 itag, lcomm, ireq(nreq), ierr )
#else
            call MPI_ISend( bnds(lproc)%sbuf, llen, MPI_REAL, lproc, &
                 itag, lcomm, ireq(nreq), ierr )
#endif
         end if
      end do

      ! See if there are any points on my own process that need
      ! to be fixed up. This will only be in the case when nproc < npanels.
      do k = 1,kx
         do iq = 1,myrlen
            ! request_list is same as send_list in this case
            ! unpack_list(iq) is index into extended region
            if ( bnds(myid)%request_list_uv(iq)>0 .neqv. bnds(myid)%uv_swap(iq) ) then 
               if ( bnds(myid)%unpack_list_uv(iq)>0 ) then
                  u(ifull+bnds(myid)%unpack_list_uv(iq),k) = bnds(myid)%uv_neg(iq)*u(abs(bnds(myid)%request_list_uv(iq)),k)
               else
                  v(ifull-bnds(myid)%unpack_list_uv(iq),k) = bnds(myid)%uv_neg(iq)*u(abs(bnds(myid)%request_list_uv(iq)),k)
               end if
            else
               if ( bnds(myid)%unpack_list_uv(iq)>0 ) then
                  u(ifull+bnds(myid)%unpack_list_uv(iq),k) = bnds(myid)%uv_neg(iq)*v(abs(bnds(myid)%request_list_uv(iq)),k)
               else
                  v(ifull-bnds(myid)%unpack_list_uv(iq),k) = bnds(myid)%uv_neg(iq)*v(abs(bnds(myid)%request_list_uv(iq)),k)
               end if
            end if   
         end do   
      end do
      
      ! Unpack incomming messages
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
               rproc = neighlist(iproc)
               iqq = 0
               if ( fsvwu ) then
                  iqz = iqq - bnds(rproc)%rlen_sv_bg + 1 
                  iqlen = bnds(rproc)%rlen_wu_fn - bnds(rproc)%rlen_sv_bg + 1
                  do k = 1,kx
                     do iq = bnds(rproc)%rlen_sv_bg,bnds(rproc)%rlen_wu_fn
                        ! unpack_list(iq) is index into extended region
                        if ( bnds(rproc)%unpack_list_uv(iq) > 0 ) then
                           u(ifull+bnds(rproc)%unpack_list_uv(iq),k) = bnds(rproc)%rbuf(iqz+iq+(k-1)*iqlen)
                        else
                           v(ifull-bnds(rproc)%unpack_list_uv(iq),k) = bnds(rproc)%rbuf(iqz+iq+(k-1)*iqlen)
                        end if
                     end do   
                  end do   
                  iqq = iqq + iqlen*kx
               end if
               if ( fnveu ) then
                  iqz = iqq - bnds(rproc)%rlen_nv_bg + 1                 
                  iqlen = bnds(rproc)%rlen_eu_fn - bnds(rproc)%rlen_nv_bg + 1
                  do k = 1,kx
                     do iq = bnds(rproc)%rlen_nv_bg,bnds(rproc)%rlen_eu_fn
                        ! unpack_list(iq) is index into extended region
                        if ( bnds(rproc)%unpack_list_uv(iq) > 0 ) then
                           u(ifull+bnds(rproc)%unpack_list_uv(iq),k) = bnds(rproc)%rbuf(iqz+iq+(k-1)*iqlen)
                        else
                           v(ifull-bnds(rproc)%unpack_list_uv(iq),k) = bnds(rproc)%rbuf(iqz+iq+(k-1)*iqlen)
                        end if
                     end do   
                  end do   
                  iqq = iqq + iqlen*kx
               end if         
               if ( fssvwwu ) then
                  iqz = iqq - bnds(rproc)%rlen_ssv_bg + 1                 
                  iqlen = bnds(rproc)%rlen_wwu_fn - bnds(rproc)%rlen_ssv_bg + 1
                  do k = 1,kx
                     do iq = bnds(rproc)%rlen_ssv_bg,bnds(rproc)%rlen_wwu_fn
                        ! unpack_list(iq) is index into extended region
                        if ( bnds(rproc)%unpack_list_uv(iq) > 0 ) then
                           u(ifull+bnds(rproc)%unpack_list_uv(iq),k) = bnds(rproc)%rbuf(iqz+iq+(k-1)*iqlen)
                        else
                           v(ifull-bnds(rproc)%unpack_list_uv(iq),k) = bnds(rproc)%rbuf(iqz+iq+(k-1)*iqlen)
                        end if
                     end do   
                  end do  
                  iqq = iqq + iqlen*kx
               end if         
               if ( fnnveeu ) then
                  iqz = iqq - bnds(rproc)%rlen_nnv_bg + 1                 
                  iqlen = bnds(rproc)%rlen_eeu_fn - bnds(rproc)%rlen_nnv_bg + 1
                  do k = 1,kx
                     do iq = bnds(rproc)%rlen_nnv_bg,bnds(rproc)%rlen_eeu_fn
                        ! unpack_list(iq) is index into extended region
                        if ( bnds(rproc)%unpack_list_uv(iq) > 0 ) then
                           u(ifull+bnds(rproc)%unpack_list_uv(iq),k) = bnds(rproc)%rbuf(iqz+iq+(k-1)*iqlen)
                        else
                           v(ifull-bnds(rproc)%unpack_list_uv(iq),k) = bnds(rproc)%rbuf(iqz+iq+(k-1)*iqlen)
                        end if
                     end do   
                  end do   
                  iqq = iqq + iqlen*kx
               end if     
               if ( fsuev ) then
                  iqz = iqq - bnds(rproc)%rlen_su_bg + 1                 
                  iqlen = bnds(rproc)%rlen_ev_fn - bnds(rproc)%rlen_su_bg + 1
                  do k = 1,kx
                     do iq = bnds(rproc)%rlen_su_bg,bnds(rproc)%rlen_ev_fn
                        ! unpack_list(iq) is index into extended region
                        if ( bnds(rproc)%unpack_list_uv(iq) > 0 ) then
                           u(ifull+bnds(rproc)%unpack_list_uv(iq),k) = bnds(rproc)%rbuf(iqz+iq+(k-1)*iqlen)
                        else
                           v(ifull-bnds(rproc)%unpack_list_uv(iq),k) = bnds(rproc)%rbuf(iqz+iq+(k-1)*iqlen)
                        end if
                     end do   
                  end do
                  iqq = iqq + iqlen*kx
               end if
            end if ! mproc <= rreq
         end do
      end do

   end subroutine boundsuv3
   
   subroutine deptsync(nface,xg,yg)
      ! Different levels will have different winds, so the list of points is
      ! different on each level.
      ! xg ranges from 0.5 to il+0.5 on a face. A given process range
      ! is 0.5+ioff to 0.5+ipan+ioff
      ! Assignment of points to processes needs to match what ints does
      ! in case there's an exact edge point.
      ! Because of the boundary region, the range [0:ipan+1) can be handled.
      ! Need floor(xxg) in range [0:ipan]
      integer, dimension(:,:), intent(in) :: nface
      real, dimension(:,:), intent(in) :: xg, yg
      integer :: iproc, jproc, dproc
      integer :: ip, jp, xn, kx
      integer :: iq, k
      integer :: rcount
      integer(kind=4), save :: itag=99
      integer(kind=4) :: ierr, llen, ncount, lproc
      integer(kind=4) :: ldone, lcomm
      integer(kind=4), dimension(MPI_STATUS_SIZE,2*neighnum) :: status
      integer(kind=4), dimension(2*neighnum) :: donelist
      real, dimension(4,maxbuflen*maxvertlen,neighnum) :: buf_dpoints, buf_dbuf 

      ! This does nothing in the one process case
      if ( neighnum < 1 ) return
      
      kx = size(nface,2)
      dslen(:) = 0
      drlen(:) = 0
      dproc = 0
      lcomm = comm_world
      itag = mod(itag + 1, 10000)
      
      ! In this case the length of each buffer is unknown and will not
      ! be symmetric between processes. Therefore need to get the length
      ! from the message status
      do iproc = 1,neighnum
         lproc = neighlist(iproc)  ! Recv from
         ! Use the maximum size in the recv call.
         llen = 4*bnds(lproc)%len
#ifdef i8r8
         call MPI_IRecv( buf_dpoints(:,:,iproc), llen, MPI_DOUBLE_PRECISION, lproc, &
                      itag, lcomm, ireq(iproc), ierr )
#else
         call MPI_IRecv( buf_dpoints(:,:,iproc), llen, MPI_REAL, lproc, &
                      itag, lcomm, ireq(iproc), ierr )
#endif
      end do
      
      ! Calculate request list
      do k = 1,kx
         do iq = 1,ifull
            ! If point is on a different process, add to a list 
            ip = min( il_g, max( 1, nint(xg(iq,k)) ) )
            jp = min( il_g, max( 1, nint(yg(iq,k)) ) )
            iproc = fproc(ip,jp,nface(iq,k)) ! process that owns global grid point
            if ( iproc/=myid ) then
               dproc = neighmap(iproc) ! returns 0 if not in neighlist
               ! Add this point to the list of requests I need to send to iproc
               dslen(dproc) = dslen(dproc) + 1
               ! Limit request index to valid range to avoid seg fault
               xn = max( min( dslen(dproc), bnds(iproc)%len ), 1 )
               ! Since nface is a small integer it can be exactly represented by a
               ! real. It is simpler to send like this than use a proper structure.
               dbuf(dproc)%a(xn,1) = real(nface(iq,k))
               dbuf(dproc)%a(xn,2) = xg(iq,k)
               dbuf(dproc)%a(xn,3) = yg(iq,k)
               dbuf(dproc)%a(xn,4) = real(k)
               dindex(dproc)%a(xn,1) = iq
               dindex(dproc)%a(xn,2) = k
            end if
         end do
      end do
 
      ! Check for errors
      if ( dslen(0) > 0 ) then
         ip = min( il_g, max( 1, nint(dbuf(0)%a(1,2)) ) )
         jp = min( il_g, max( 1, nint(dbuf(0)%a(1,3)) ) )
         iproc = fproc(ip,jp,nint(dbuf(0)%a(1,1)))
         write(6,*) "myid,dslen,len ",myid,dslen(0),0
         write(6,*) "Example error iq,k,iproc ",dindex(0)%a(1,1),dindex(0)%a(1,2),iproc
         write(6,*) "dbuf ",dbuf(0)%a(1,:)
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
            write(6,*) "dbuf ",dbuf(dproc)%a(1,:)
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
      do iproc = 1,neighnum
         lproc = neighlist(iproc)  ! Send to
         ! Send, even if length is zero
         llen = 4*dslen(iproc)
         buf_dbuf(1:4,1:dslen(iproc),iproc) = transpose( dbuf(iproc)%a(1:dslen(iproc),1:4) )
#ifdef i8r8
         call MPI_ISend( buf_dbuf(:,:,iproc), llen, MPI_DOUBLE_PRECISION, lproc, &
                 itag, lcomm, ireq(neighnum+iproc), ierr )
#else
         call MPI_ISend( buf_dbuf(:,:,iproc), llen, MPI_REAL, lproc, &
                 itag, lcomm, ireq(neighnum+iproc), ierr )
#endif
      end do
      
      ! Unpack incomming messages
      nreq = 2*neighnum
      rcount = nreq
      do while ( rcount > 0 )
         call START_LOG(mpiwait_begin)
         call MPI_Waitsome( nreq, ireq, ldone, donelist, status, ierr )
         call END_LOG(mpiwait_end)
         rcount = rcount - ldone
         do jproc = 1,ldone
            iproc = donelist(jproc)
            if ( iproc <= neighnum ) then
               ! Now get the actual sizes from the status
#ifdef i8r8
               call MPI_Get_count( status(:,jproc), MPI_DOUBLE_PRECISION, ncount, ierr )
#else
               call MPI_Get_count( status(:,jproc), MPI_REAL, ncount, ierr )
#endif
               drlen(iproc) = ncount/4
               dpoints(iproc)%a(1:drlen(iproc),1:4) = transpose( buf_dpoints(1:4,1:drlen(iproc),iproc) )
            end if
         end do
      end do

   end subroutine deptsync

   subroutine intssync_send3
      
      call intssync_send4(1)

   end subroutine intssync_send3
   
   subroutine intssync_send4(ntr)
      integer, intent(in) :: ntr
      integer :: iproc
      integer(kind=4), save :: itag=98
      integer(kind=4) :: ierr, llen, lproc, lcomm
      
      if ( ntr > nagg ) then
         write(6,*) "ERROR: Internal error in intssync_send.  ntr > nagg"
         call ccmpi_abort(-1)
      end if
      
      itag = mod(itag + 1, 10000)

      ! When sending the results, roles of dslen and drlen are reversed
      lcomm = comm_world
      nreq = 0
      do iproc = 1,neighnum
         if ( dslen(iproc) > 0 ) then
            nreq = nreq + 1
            rlist(nreq) = iproc
            llen = dslen(iproc)*ntr
            lproc = neighlist(iproc)  ! Recv from
#ifdef i8r8
            call MPI_IRecv( dbuf(iproc)%b, llen, MPI_DOUBLE_PRECISION, lproc, itag, &
                            lcomm, ireq(nreq), ierr )
#else
            call MPI_IRecv( dbuf(iproc)%b, llen, MPI_REAL, lproc, itag, &
                            lcomm, ireq(nreq), ierr )
#endif
         end if
      end do
      rreq = nreq
      do iproc = neighnum,1,-1
         if ( drlen(iproc) > 0 ) then
            nreq = nreq + 1
            llen = drlen(iproc)*ntr
            lproc = neighlist(iproc)  ! Send to
#ifdef i8r8
            call MPI_ISend( sextra(iproc)%a, llen, MPI_DOUBLE_PRECISION, lproc, itag, & 
                            lcomm, ireq(nreq), ierr )
#else
            call MPI_ISend( sextra(iproc)%a, llen, MPI_REAL, lproc, itag, & 
                            lcomm, ireq(nreq), ierr )
#endif
         end if
      end do

   end subroutine intssync_send4

   subroutine intssync_recv3(s)
      real, dimension(:,:), intent(inout) :: s
      real, dimension(ifull,size(s,2),1) :: s_l
      
      s_l(1:ifull,:,1) = s(1:ifull,:)
      call intssync_recv4(s_l)
      s(1:ifull,:) = s_l(1:ifull,:,1)

   end subroutine intssync_recv3

   subroutine intssync_recv4(s)
      real, dimension(:,:,:), intent(inout) :: s
      integer :: iproc, mproc, iq, jproc
      integer :: rcount, ntr, l
      integer(kind=4) :: ierr, ldone
      integer(kind=4), dimension(2*neighnum) :: donelist

      ntr = size(s,3)
      if ( ntr > nagg ) then
         write(6,*) "ERROR: Internal error in intssync_recv.  ntr > nagg"
         call ccmpi_abort(-1)
      end if
      
      ! Unpack incomming messages
      rcount = nreq
      do while ( rcount > 0 )
         call START_LOG(mpiwait_begin)
         call MPI_Waitsome( nreq, ireq, ldone, donelist, MPI_STATUSES_IGNORE, ierr )
         call END_LOG(mpiwait_end)
         rcount = rcount - ldone
         do jproc = 1,ldone
            mproc = donelist(jproc)
            if ( mproc <= rreq ) then
               iproc = rlist(mproc)
               do l = 1,ntr
                  do iq = 1,dslen(iproc)
                     s(dindex(iproc)%a(iq,1),dindex(iproc)%a(iq,2),l) &
                         = dbuf(iproc)%b(iq+(l-1)*dslen(iproc))
                  end do   
               end do
            end if   
         end do
      end do

   end subroutine intssync_recv4
   
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

      ! Which process has this point
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
         ! If it's on this process, just use the local index
         larray(n) = iqg2iq(iqq)
      end if
      
   end subroutine fix_index2

end module cc_mpi_point

