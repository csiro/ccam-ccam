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
    
! Collective MPI communication subroutines    
    
module cc_mpi_collective

   use cc_mpi_common
   
   private
   
   integer(kind=4), allocatable, dimension(:), save, public ::              &
      specmap_req                                                          ! gather map required for spectral filter
   integer(kind=4), allocatable, dimension(:), save, public ::              &
      specmap_recv                                                         ! gather map recieved for spectral filter
   integer(kind=4), allocatable, dimension(:), save, public ::              &
      specmap_send                                                         ! gather map sent for spectral filter
   integer, allocatable, dimension(:), save, public :: specmap_indx        ! index for nodepack
   integer, save, public :: specmap_indxlen                                ! number of messages for nodepack
   integer, allocatable, dimension(:), save, public :: specmap_ext         ! gather map for spectral filter (includes filter final
                                                                           ! pass for sparse arrays)
   real, dimension(:,:,:,:,:), pointer, save, private :: nodepack          ! node buffer for spectral filter
   integer, save, private :: nodepack_win
   type globalpack_info
     real, allocatable, dimension(:,:,:) :: localdata
   end type globalpack_info
   type(globalpack_info), allocatable, dimension(:,:,:), save, private ::   & 
      globalpack                                                           ! store sparse global arrays for spectral filter

   public :: ccmpi_distribute, ccmpi_gather, ccmpi_gatherr8,                &
             ccmpi_distributer8, ccmpi_gatherall
   public :: allocateglobalpack, copyglobalpack,                            &
             ccmpi_gathermap_send, ccmpi_gathermap_recv, getglobalpack_v,   &
             setglobalpack_v, ccmpi_gathermap_wait, deallocateglobalpack
   public :: ccglobal_posneg   
   public :: readglobvar, writeglobvar

   interface ccmpi_gather
      module procedure host_gather2, host_gather3, host_gather4
      module procedure proc_gather2, proc_gather3, proc_gather4
   end interface
   interface ccmpi_gatherr8
      module procedure host_gather2r8, host_gather3r8, host_gather4r8
      module procedure proc_gather2r8, proc_gather3r8, proc_gather4r8
   end interface
   interface ccmpi_distribute
      module procedure host_distribute2, host_distribute3, host_distribute4
      module procedure proc_distribute2, proc_distribute3, proc_distribute4      
      module procedure host_distribute2i, host_distribute3i, host_distribute4i
      module procedure proc_distribute2i, proc_distribute3i, proc_distribute4i      
   end interface
   interface ccmpi_distributer8
      module procedure host_distribute2r8, host_distribute3r8, host_distribute4r8
      module procedure proc_distribute2r8, proc_distribute3r8, proc_distribute4r8      
   end interface
   interface ccmpi_gatherall
      module procedure ccmpi_gatherall2, ccmpi_gatherall3
   end interface
   interface ccmpi_gathermap_send
      module procedure ccmpi_gathermap_send2, ccmpi_gathermap_send3
   end interface
   interface ccmpi_gathermap_recv
      module procedure ccmpi_gathermap_recv2, ccmpi_gathermap_recv3
   end interface
   interface ccglobal_posneg
      ! ccglobal_posneg2 and ccglobal_posneg2o are equivilent
      module procedure ccglobal_posneg2
      module procedure ccglobal_posneg3, ccglobal_posneg4
      module procedure ccglobal_posneg3o, ccglobal_posneg4o
   end interface
   interface readglobvar
      module procedure readglobvar2, readglobvar3, readglobvar2i
   end interface
   interface writeglobvar
      module procedure writeglobvar2, writeglobvar3
   end interface

contains
   
   subroutine host_distribute2(af,a1)
      real, dimension(ifull), intent(out) :: af
      real, dimension(ifull_g), intent(in) :: a1
      real, dimension(ifull,1,1) :: af_l
      real, dimension(ifull_g,1,1) :: a1_l
      
      a1_l(1:ifull_g,1,1) = a1(1:ifull_g)
      call host_distribute4(af_l,a1_l)
      af(1:ifull) = af_l(1:ifull,1,1)
      
   end subroutine host_distribute2

   subroutine proc_distribute2(af)
      real, dimension(ifull), intent(out) :: af
      real, dimension(ifull,1,1) :: af_l
      
      call proc_distribute4(af_l)
      af(1:ifull) = af_l(1:ifull,1,1)
      
   end subroutine proc_distribute2   
   
   subroutine host_distribute2r8(af,a1)
      real(kind=8), dimension(ifull), intent(out) :: af
      real(kind=8), dimension(ifull_g), intent(in) :: a1
      real(kind=8), dimension(ifull,1,1) :: af_l
      real(kind=8), dimension(ifull_g,1,1) :: a1_l
      
      a1_l(1:ifull_g,1,1) = a1(1:ifull_g)
      call host_distribute4r8(af_l,a1_l)
      af(1:ifull) = af_l(1:ifull,1,1)
      
   end subroutine host_distribute2r8

   subroutine proc_distribute2r8(af)
      real(kind=8), dimension(ifull), intent(out) :: af
      real(kind=8), dimension(ifull,1,1) :: af_l
      
      call proc_distribute4r8(af_l)
      af(1:ifull) = af_l(1:ifull,1,1)
      
   end subroutine proc_distribute2r8
      
   subroutine host_distribute2i(af,a1)
      integer, dimension(ifull), intent(out) :: af
      integer, dimension(ifull_g), intent(in) :: a1
      integer, dimension(ifull,1,1) :: af_l
      integer, dimension(ifull_g,1,1) :: a1_l

      a1_l(1:ifull_g,1,1) = a1(1:ifull_g)
      call host_distribute4i(af_l,a1_l)
      af(1:ifull) = af_l(1:ifull,1,1)

   end subroutine host_distribute2i
   
   subroutine proc_distribute2i(af)
      integer, dimension(ifull), intent(out) :: af
      integer, dimension(ifull,1,1) :: af_l

      call proc_distribute4i(af_l)
      af(1:ifull) = af_l(1:ifull,1,1)

   end subroutine proc_distribute2i   

   subroutine host_distribute3(af,a1)
      real, dimension(:,:), intent(out) :: af
      real, dimension(:,:), intent(in) :: a1
      real, dimension(ifull,size(af,2),1) :: af_l
      real, dimension(ifull_g,size(af,2),1) :: a1_l
      
      a1_l(1:ifull_g,1:size(af,2),1) = a1(1:ifull_g,1:size(af,2))
      call host_distribute4(af_l,a1_l)
      af(1:ifull,:) = af_l(1:ifull,:,1)  

   end subroutine host_distribute3
   
   subroutine proc_distribute3(af)
      real, dimension(:,:), intent(out) :: af
      real, dimension(ifull,size(af,2),1) :: af_l
      
      call proc_distribute4(af_l)
      af(1:ifull,:) = af_l(1:ifull,:,1)  

   end subroutine proc_distribute3

   subroutine host_distribute3r8(af,a1)
      real(kind=8), dimension(:,:), intent(out) :: af
      real(kind=8), dimension(:,:), intent(in) :: a1
      real(kind=8), dimension(ifull,size(af,2),1) :: af_l
      real(kind=8), dimension(ifull_g,size(af,2),1) :: a1_l
      
      a1_l(1:ifull_g,1:size(af,2),1) = a1(1:ifull_g,1:size(af,2))
      call host_distribute4r8(af_l,a1_l)
      af(1:ifull,:) = af_l(1:ifull,:,1)

   end subroutine host_distribute3r8
   
   subroutine proc_distribute3r8(af)
      real(kind=8), dimension(:,:), intent(out) :: af
      real(kind=8), dimension(ifull,size(af,2),1) :: af_l
      
      call proc_distribute4r8(af_l)
      af(1:ifull,:) = af_l(1:ifull,:,1)

   end subroutine proc_distribute3r8

   subroutine host_distribute3i(af,a1)
      integer, dimension(:,:), intent(out) :: af
      integer, dimension(:,:), intent(in) :: a1
      integer, dimension(ifull,size(af,2),1) :: af_l
      integer, dimension(ifull_g,size(af,2),1) :: a1_l
      
      a1_l(1:ifull_g,1:size(af,2),1) = a1(1:ifull_g,1:size(af,2))
      call host_distribute4i(af_l,a1_l)
      af(1:ifull,:) = af_l(1:ifull,:,1)

   end subroutine host_distribute3i
   
   subroutine proc_distribute3i(af)
      integer, dimension(:,:), intent(out) :: af
      integer, dimension(ifull,size(af,2),1) :: af_l
      
      call proc_distribute4i(af_l)
      af(1:ifull,:) = af_l(1:ifull,:,1)

   end subroutine proc_distribute3i
  
   subroutine host_distribute4(af,a1)
      ! Convert standard 1D arrays to face form and distribute to processes
      ! This is also used for tracers, so second dimension is not necessarily
      ! the number of levels
      real, dimension(:,:,:), intent(out) :: af
      real, dimension(:,:,:), intent(in) :: a1
      real, dimension(ifull,size(af,2),size(af,3),0:nproc-1) :: sbuf
      real, dimension(ifull,size(af,2),size(af,3)) :: aftemp
      integer :: j, n, k, l, iq, iproc
      integer :: npoff, ipoff, jpoff ! Offsets for target
      integer :: slen, kx, lx
      integer(kind=4) :: ierr, lsize, lcomm

      kx = size(af,2)
      lx = size(af,3)

      ! map array in order of process rank
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

      lsize = ifull*kx*lx
      lcomm = comm_world
      call START_LOG(scatter_begin) 
#ifdef i8r8      
      call MPI_Scatter( sbuf, lsize, MPI_DOUBLE_PRECISION, aftemp, lsize, MPI_DOUBLE_PRECISION, 0_4, lcomm, ierr )
#else
      call MPI_Scatter( sbuf, lsize, MPI_REAL, aftemp, lsize, MPI_REAL, 0_4, lcomm, ierr )
#endif
      call END_LOG(scatter_end)
      af(1:ifull,1:kx,1:lx) = aftemp(1:ifull,1:kx,1:lx)

   end subroutine host_distribute4

   subroutine proc_distribute4(af)
      ! Convert standard 1D arrays to face form and distribute to processes
      ! This is also used for tracers, so second dimension is not necessarily
      ! the number of levels
      real, dimension(:,:,:), intent(out) :: af
      real, dimension(1,1,1,1) :: sbuf
      real, dimension(ifull,size(af,2),size(af,3)) :: aftemp
      integer :: kx, lx
      integer(kind=4) :: ierr, lsize, lcomm

      kx = size(af,2)
      lx = size(af,3)
      lsize = ifull*kx*lx
      lcomm = comm_world
      call START_LOG(scatter_begin) 
#ifdef i8r8      
      call MPI_Scatter( sbuf, lsize, MPI_DOUBLE_PRECISION, aftemp, lsize, MPI_DOUBLE_PRECISION, 0_4, lcomm, ierr )
#else
      call MPI_Scatter( sbuf, lsize, MPI_REAL, aftemp, lsize, MPI_REAL, 0_4, lcomm, ierr )
#endif
      call END_LOG(scatter_end)
      af(1:ifull,1:kx,1:lx) = aftemp(1:ifull,1:kx,1:lx) 

   end subroutine proc_distribute4

   subroutine host_distribute4r8(af,a1)
      ! Convert standard 1D arrays to face form and distribute to processes
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
      
      ! map array in order of process rank
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

      lsize = ifull*kx*lx
      lcomm = comm_world
      call START_LOG(scatter_begin) 
      call MPI_Scatter( sbuf, lsize, MPI_DOUBLE_PRECISION, aftemp, lsize, MPI_DOUBLE_PRECISION, 0_4, lcomm, ierr )
      call END_LOG(scatter_end)
      af(1:ifull,1:kx,1:lx) = aftemp(1:ifull,1:kx,1:lx)

   end subroutine host_distribute4r8
   
   subroutine proc_distribute4r8(af)
      ! Convert standard 1D arrays to face form and distribute to processes
      real(kind=8), dimension(:,:,:), intent(out) :: af
      real(kind=8), dimension(1,1,1,1) :: sbuf
      real(kind=8), dimension(ifull,size(af,2),size(af,3)) :: aftemp
      integer :: kx, lx
      integer(kind=4) :: ierr, lsize, lcomm

      kx = size(af,2)
      lx = size(af,3)
      lsize = ifull*kx*lx
      lcomm = comm_world
      call START_LOG(scatter_begin) 
      call MPI_Scatter( sbuf, lsize, MPI_DOUBLE_PRECISION, aftemp, lsize, MPI_DOUBLE_PRECISION, 0_4, lcomm, ierr ) 
      call END_LOG(scatter_end)
      af(1:ifull,1:kx,1:lx) = aftemp(1:ifull,1:kx,1:lx)

   end subroutine proc_distribute4r8
   
   subroutine host_distribute4i(af,a1)
      ! Convert standard 1D arrays to face form and distribute to processes
      ! This is also used for tracers, so second dimension is not necessarily
      ! the number of levels
      integer, dimension(:,:,:), intent(out) :: af
      integer, dimension(:,:,:), intent(in) :: a1
      integer, dimension(ifull,size(af,2),size(af,3),0:nproc-1) :: sbuf
      integer, dimension(ifull,size(af,2),size(af,3)) :: aftemp
      integer :: j, n, k, l, iq, iproc
      integer :: npoff, ipoff, jpoff ! Offsets for target
      integer :: slen, kx, lx
      integer(kind=4) :: ierr, lsize, lcomm

      kx = size(af,2)
      lx = size(af,3)

      ! map array in order of process rank
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

      lsize = ifull*kx*lx
      lcomm = comm_world
      call START_LOG(scatter_begin) 
#ifdef i8r8      
      call MPI_Scatter( sbuf, lsize, MPI_INTEGER8, aftemp, lsize, MPI_INTEGER8, 0_4, lcomm, ierr )
#else
      call MPI_Scatter( sbuf, lsize, MPI_INTEGER, aftemp, lsize, MPI_INTEGER, 0_4, lcomm, ierr )
#endif
      call END_LOG(scatter_end)
      af(1:ifull,1:kx,1:lx) = aftemp(1:ifull,1:kx,1:lx)

   end subroutine host_distribute4i

   subroutine proc_distribute4i(af)
      ! Convert standard 1D arrays to face form and distribute to processes
      ! This is also used for tracers, so second dimension is not necessarily
      ! the number of levels
      integer, dimension(:,:,:), intent(out) :: af
      integer, dimension(1,1,1,1) :: sbuf
      integer, dimension(ifull,size(af,2),size(af,3)) :: aftemp
      integer :: kx, lx
      integer(kind=4) :: ierr, lsize, lcomm

      kx = size(af,2)
      lx = size(af,3)
      lsize = ifull*kx*lx
      lcomm = comm_world
      call START_LOG(scatter_begin) 
#ifdef i8r8      
      call MPI_Scatter( sbuf, lsize, MPI_INTEGER8, aftemp, lsize, MPI_INTEGER8, 0_4, lcomm, ierr )
#else
      call MPI_Scatter( sbuf, lsize, MPI_INTEGER, aftemp, lsize, MPI_INTEGER, 0_4, lcomm, ierr )
#endif
      call END_LOG(scatter_end)
      af(1:ifull,1:kx,1:lx) = aftemp(1:ifull,1:kx,1:lx) 

   end subroutine proc_distribute4i

   subroutine host_gather2(a,ag)
      real, dimension(ifull), intent(in) :: a
      real, dimension(ifull_g), intent(out) :: ag
      real, dimension(ifull,1,1) :: a_l
      real, dimension(ifull_g,1,1) :: ag_l
      
      a_l(1:ifull,1,1) = a(1:ifull)
      call host_gather4(a_l,ag_l)
      ag(1:ifull_g) = ag_l(1:ifull_g,1,1)

   end subroutine host_gather2

   subroutine proc_gather2(a)
      real, dimension(ifull), intent(in) :: a
      real, dimension(ifull,1,1) :: a_l
      
      a_l(1:ifull,1,1) = a(1:ifull)
      call proc_gather4(a_l) 

   end subroutine proc_gather2
   
   subroutine host_gather2r8(a,ag)
      real(kind=8), dimension(ifull), intent(in) :: a
      real(kind=8), dimension(ifull_g), intent(out) :: ag
      real(kind=8), dimension(ifull,1,1) :: a_l
      real(kind=8), dimension(ifull_g,1,1) :: ag_l
      
      a_l(1:ifull,1,1) = a(1:ifull)
      call host_gather4r8(a_l,ag_l)
      ag(1:ifull_g) = ag_l(1:ifull_g,1,1)
      
   end subroutine host_gather2r8
   
   subroutine proc_gather2r8(a)
      real(kind=8), dimension(ifull), intent(in) :: a
      real(kind=8), dimension(ifull,1,1) :: a_l
      
      a_l(1:ifull,1,1) = a(1:ifull)
      call proc_gather4r8(a_l) 
      
   end subroutine proc_gather2r8
   
   subroutine host_gather3(a,ag)
      real, dimension(:,:), intent(in) :: a
      real, dimension(:,:), intent(out) :: ag
      real, dimension(ifull,size(a,2),1) :: a_l
      real, dimension(ifull_g,size(a,2),1) :: ag_l
      
      a_l(1:ifull,:,1) = a(1:ifull,:)
      call host_gather4(a_l,ag_l)
      ag(1:ifull_g,1:size(a,2)) = ag_l(1:ifull_g,1:size(a,2),1)

   end subroutine host_gather3
   
   subroutine proc_gather3(a)
      real, dimension(:,:), intent(in) :: a
      real, dimension(ifull,size(a,2),1) :: a_l
      
      a_l(1:ifull,:,1) = a(1:ifull,:)
      call proc_gather4(a_l) 

   end subroutine proc_gather3   

   subroutine host_gather3r8(a,ag)
      real(kind=8), dimension(:,:), intent(in) :: a
      real(kind=8), dimension(:,:), intent(out) :: ag
      real(kind=8), dimension(ifull,size(a,2),1) :: a_l
      real(kind=8), dimension(ifull_g,size(a,2),1) :: ag_l
      
      a_l(1:ifull,:,1) = a(1:ifull,:)
      call host_gather4r8(a_l,ag_l)
      ag(1:ifull_g,1:size(a,2)) = ag_l(1:ifull_g,1:size(a,2),1)

   end subroutine host_gather3r8
   
   subroutine proc_gather3r8(a)
      real(kind=8), dimension(:,:), intent(in) :: a
      real(kind=8), dimension(ifull,size(a,2),1) :: a_l
      
      a_l(1:ifull,:,1) = a(1:ifull,:)
      call proc_gather4r8(a_l) 

   end subroutine proc_gather3r8
   
   subroutine host_gather4(a,ag)
      ! Collect global arrays.
      real, dimension(:,:,:), intent(in) :: a
      real, dimension(:,:,:), intent(out) :: ag
      real, dimension(ifull,size(a,2),size(a,3),0:nproc-1) :: abuf
      real, dimension(ifull,size(a,2),size(a,3)) :: atemp
      integer :: iproc, ipoff, jpoff, npoff
      integer :: j, n, k, l, iq, iqg, kx, lx
      integer(kind=4) :: ierr, lsize, lcomm

      kx = size(a,2)
      lx = size(a,3)
      lsize = ifull*kx*lx
      lcomm = comm_world
      atemp(1:ifull,1:kx,1:lx) = a(1:ifull,1:kx,1:lx)
      call START_LOG(gather_begin)
#ifdef i8r8
      call MPI_Gather( atemp, lsize, MPI_DOUBLE_PRECISION, abuf, lsize, MPI_DOUBLE_PRECISION, 0_4, lcomm, ierr )
#else
      call MPI_Gather( atemp, lsize, MPI_REAL, abuf, lsize, MPI_REAL, 0_4, lcomm, ierr )
#endif
      call END_LOG(gather_end)

      ! map array in order of process rank
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

   end subroutine host_gather4
   
   subroutine proc_gather4(a)
      ! Collect global arrays.
      real, dimension(:,:,:), intent(in) :: a
      real, dimension(1,1,1,1) :: abuf
      real, dimension(ifull,size(a,2),size(a,3)) :: atemp
      integer :: kx, lx
      integer(kind=4) :: ierr, lsize, lcomm

      kx = size(a,2)
      lx = size(a,3)
      lsize = ifull*kx*lx
      lcomm = comm_world
      atemp(1:ifull,1:kx,1:lx) = a(1:ifull,1:kx,1:lx)
      call START_LOG(gather_begin)
#ifdef i8r8
      call MPI_Gather( atemp, lsize, MPI_DOUBLE_PRECISION, abuf, lsize, MPI_DOUBLE_PRECISION, 0_4, lcomm, ierr )
#else
      call MPI_Gather( atemp, lsize, MPI_REAL, abuf, lsize, MPI_REAL, 0_4, lcomm, ierr )
#endif
      call END_LOG(gather_end)

   end subroutine proc_gather4

   subroutine host_gather4r8(a,ag)
      ! Collect global arrays.
      real(kind=8), dimension(:,:,:), intent(in) :: a
      real(kind=8), dimension(:,:,:), intent(out) :: ag
      real(kind=8), dimension(ifull,size(a,2),size(a,3),0:nproc-1) :: abuf
      real(kind=8), dimension(ifull,size(a,2),size(a,3)) :: atemp
      integer :: iproc, ipoff, jpoff, npoff
      integer :: j, n, k, l, iq, iqg, kx, lx
      integer(kind=4) :: ierr, lsize, lcomm

      kx = size(a,2)
      lx = size(a,3)
      lsize = ifull*kx*lx
      lcomm = comm_world
      atemp(1:ifull,1:kx,1:lx) = a(1:ifull,1:kx,1:lx)
      call START_LOG(gather_begin)
      call MPI_Gather( atemp, lsize, MPI_DOUBLE_PRECISION, abuf, lsize, MPI_DOUBLE_PRECISION, 0_4, lcomm, ierr )
      call END_LOG(gather_end)

      ! map array in order of process rank
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

   end subroutine host_gather4r8
   
   subroutine proc_gather4r8(a)
      ! Collect global arrays.
      real(kind=8), dimension(:,:,:), intent(in) :: a
      real(kind=8), dimension(1,1,1,1) :: abuf
      real(kind=8), dimension(ifull,size(a,2),size(a,3)) :: atemp
      integer :: kx, lx
      integer(kind=4) :: ierr, lsize, lcomm
      
      kx = size(a,2)
      lx = size(a,3)
      lsize = ifull*kx*lx
      lcomm = comm_world
      atemp(1:ifull,1:kx,1:lx) = a(1:ifull,1:kx,1:lx)
      call START_LOG(gather_begin)
      call MPI_Gather( atemp, lsize, MPI_DOUBLE_PRECISION, abuf, lsize, MPI_DOUBLE_PRECISION, 0_4, lcomm, ierr )
      call END_LOG(gather_end)

   end subroutine proc_gather4r8
   
   subroutine ccmpi_gatherall2(a,ag)
      real, dimension(ifull), intent(in) :: a
      real, dimension(ifull_g), intent(out) :: ag
      real, dimension(ifull,1) :: a_l
      real, dimension(ifull_g,1) :: ag_l
      
      a_l(1:ifull,1) = a(1:ifull)
      call ccmpi_gatherall3(a_l,ag_l)
      ag(1:ifull_g) = ag_l(1:ifull_g,1)

   end subroutine ccmpi_gatherall2
   
   subroutine ccmpi_gatherall3(a,ag)
      ! Collect global arrays.
      real, dimension(:,:), intent(in) :: a
      real, dimension(:,:), intent(out) :: ag
      real, dimension(ifull,size(a,2),0:nproc-1) :: abuf
      real, dimension(ifull,size(a,2)) :: atemp
      integer :: ipoff, jpoff, npoff
      integer :: j, n, k, iq, iqg, kx, iproc
      integer(kind=4) :: ierr, lsize, lcomm

      kx = size(a,2)
      lsize = ifull*kx
      lcomm = comm_world
      atemp(:,:) = a(1:ifull,1:kx)
      call START_LOG(allgather_begin)
#ifdef i8r8
      call MPI_AllGather( atemp, lsize, MPI_DOUBLE_PRECISION, abuf, lsize, MPI_DOUBLE_PRECISION, lcomm, ierr )
#else
      call MPI_AllGather( atemp, lsize, MPI_REAL, abuf, lsize, MPI_REAL, lcomm, ierr )
#endif
      call END_LOG(allgather_end)

      ! map array in order of process rank
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

   end subroutine ccmpi_gatherall3

   subroutine ccmpi_gathermap_send2(a)
      real, dimension(ifull), intent(in) :: a
      real, dimension(ifull,1) :: a_l
      
      a_l(1:ifull,1) = a(1:ifull)
      call ccmpi_gathermap_send3(a_l)

   end subroutine ccmpi_gathermap_send2

   subroutine ccmpi_gathermap_send3(a)
      ! Initiate nonblocking broadcast for spectral nudging.  Data will be stored
      ! in sparse array globalpack.  This version uses shared memory for more efficent
      ! transfer of data to processes.
      real, dimension(:,:), intent(in) :: a
      integer :: w, kx
      integer(kind=4) :: ierr, lsize, lcomm, itag=52
      
      kx = size(a,2)
      lsize = ifull*kx
      lcomm = comm_world

      ! Set up the buffers to recv
      nreq = 0
      do w = 1,size(specmap_recv)
         nreq = nreq + 1
         rlist(nreq) = w
#ifdef i8r8
         call MPI_IRecv( bnds(specmap_recv(w))%rbuf, lsize, MPI_DOUBLE_PRECISION, specmap_recv(w), itag, lcomm, &
                         ireq(nreq), ierr )
#else
         call MPI_IRecv( bnds(specmap_recv(w))%rbuf, lsize, MPI_REAL, specmap_recv(w), itag, lcomm, ireq(nreq), &
                         ierr )
#endif
      end do
      rreq = nreq
      
      ! Set up the buffers to send
      bnds(myid)%sbuf(1:ifull*kx) = reshape( a(1:ifull,1:kx), (/ ifull*kx /) )
      do w = 1,size(specmap_send)
         nreq = nreq + 1
#ifdef i8r8
         call MPI_ISend( bnds(myid)%sbuf, lsize, MPI_DOUBLE_PRECISION, specmap_send(w), itag, lcomm, ireq(nreq), &
                         ierr )
#else
         call MPI_ISend( bnds(myid)%sbuf, lsize, MPI_REAL, specmap_send(w), itag, lcomm, ireq(nreq), ierr )
#endif
      end do

   end subroutine ccmpi_gathermap_send3

   subroutine ccmpi_gathermap_wait
      ! wait for nonblocking broadcast for spectral nudging to complete 
      integer(kind=4) :: ierr
      
      if ( nreq > 0 ) then
         call START_LOG(mpiwait_begin) 
         call MPI_Waitall( nreq, ireq, MPI_STATUSES_IGNORE, ierr )
         call END_LOG(mpiwait_end)
         nreq = 0
      end if   

   end subroutine ccmpi_gathermap_wait
   
   subroutine ccmpi_gathermap_recv2(kref)
      integer, intent(in) :: kref
      
      call ccmpi_gathermap_recv3(1,kref)
   
   end subroutine ccmpi_gathermap_recv2

   subroutine ccmpi_gathermap_recv3(kx,kref)
      ! Receive nonblocking broadcast for spectral nudging.  Shared memory is used to make transfer
      ! more efficient for sharing information across a node.  Data is stored in
      ! sparse array globalpack.
      integer, intent(in) :: kx, kref
      integer :: w, iproc, k, n, iq
      integer :: ipoff, jpoff, npoff
      integer :: ipak, jpak
      integer :: rcount, jproc, kproc, mproc
      integer(kind=4) :: ierr, ldone, lcomm
      integer(kind=4), dimension(size(specmap_send)+size(specmap_recv)) :: donelist

      ! kx is the size of the message we expect to recieve
      ! kref is where we want to store the message in sparse array (globalpack)
      
      ! use shared memory to transfer messages within a node
      
      if ( nreq > 0 ) then
          
         ! Unpack incomming messages into shared memory (nodepack)
         rcount = nreq
         do while ( rcount > 0 )
            call START_LOG(mpiwait_begin) 
            call MPI_Waitsome( nreq, ireq, ldone, donelist, MPI_STATUSES_IGNORE, ierr )
            call END_LOG(mpiwait_end)
            rcount = rcount - ldone
            do jproc = 1,ldone
               mproc = donelist(jproc)
               if ( mproc <= rreq ) then
                  w = rlist(mproc)
                  iproc = specmap_recv(w)
                  kproc = specmap_indx(iproc)
                  do k = 1,kx
                     do n = 1,npan
                        ! Global indices are i+ipoff, j+jpoff, n-npoff
                        iq = (n-1)*ipan*jpan + (k-1)*ifull
                        nodepack(:,:,n,k,kproc) = &
                           reshape( bnds(iproc)%rbuf(iq+1:iq+ipan*jpan), (/ ipan, jpan /) )
                     end do
                  end do
               end if
            end do
         end do
         nreq = 0
      
      else

         ! Unpack messages that have already been received 
         ! into shared memory (nodepack) 
         do w = 1,size(specmap_recv)
            iproc = specmap_recv(w)
            kproc = specmap_indx(iproc)
            do k = 1,kx
               do n = 1,npan
                  ! Global indices are i+ipoff, j+jpoff, n-npoff
                  iq = (n-1)*ipan*jpan + (k-1)*ifull
                  nodepack(:,:,n,k,kproc) = &
                    reshape( bnds(iproc)%rbuf(iq+1:iq+ipan*jpan), (/ ipan, jpan /) )
               end do
            end do
         end do
          
      end if ! if nreq>0 ..else..

      ! wait until all messages have been received before extracting data from
      ! shared memory
      call START_LOG(mpibarrier_begin) 
      lcomm = comm_node
      call MPI_Barrier( lcomm, ierr )
      call END_LOG(mpibarrier_end)

      ! unpack for process from shared memory (nodepack) into sparse array (globalpack)
      do w = 1,size(specmap_req)
         iproc = specmap_req(w)
         kproc = specmap_indx(iproc)
         call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan)
         ipak = ipoff/ipan
         jpak = jpoff/jpan
         do k = 1,kx
            do n = 1,npan
               globalpack(ipak,jpak,n-npoff)%localdata(:,:,kref+k) = &
                  nodepack(:,:,n,k,kproc)
            end do
         end do   
      end do 

      ! wait until all processes have extracted data from shared memory before
      ! allowing it to be used for new messages
      call START_LOG(mpibarrier_begin) 
      lcomm = comm_node
      call MPI_Barrier( lcomm, ierr )
      call END_LOG(mpibarrier_end)
      
   end subroutine ccmpi_gathermap_recv3
    
   subroutine setglobalpack_v(datain,ibeg,iend,k)
      ! This subroutine stores data in the sparse array (globalpack)
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
      ! This subroutine retrieves data from the sparse array (globalpack)
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
      ! to another section.  We can then refresh data for spectral
      ! nudging as each panel is processed (see spechost in nesting.f90)
      integer, intent(in) :: krefin, krefout, kx
      integer :: w, n, ncount, iproc, ipak, jpak
      integer :: ipoff, jpoff, npoff
      
      ! specmap_ext defines all allocated sparse memory in globalpack
  
      ncount = size(specmap_ext)
      do w = 1,ncount
         iproc = specmap_ext(w)
         call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan)
         ipak = ipoff/ipan
         jpak = jpoff/jpan
         do n = 1,npan
            ! Global indices are i+ipoff, j+jpoff, n-npoff
            globalpack(ipak,jpak,n-npoff)%localdata(:,:,krefout+1:krefout+kx) = &
               globalpack(ipak,jpak,n-npoff)%localdata(:,:,krefin+1:krefin+kx)
         end do
      end do
   
   end subroutine copyglobalpack

   subroutine allocateglobalpack(kx,ky)
      ! This allocates global sparse arrays for spectral nudging.
      ! Usually this is 1:kl or 1:ol in size, but for some configurations
      ! we need to store twice the data and hence use 1:2*kl or 1:2*ol.
      ! Also, the 0 index is to store the sum term for the digital filter.
      integer, intent(in) :: kx, ky
      integer :: ncount, w, ipak, jpak, n, iproc
      integer :: ipoff, jpoff, npoff
      integer :: xlen
      integer, dimension(5) :: shsize
      
      ! increase request list size
      xlen = size(specmap_recv) + size(specmap_send)
      if ( size(ireq)<xlen ) then
         deallocate( ireq )
         allocate( ireq(xlen) )
      end if
      
      ! increase receive list size
      xlen = size(specmap_recv)
      if ( size(rlist)<xlen ) then
         deallocate( rlist )
         allocate( rlist(xlen) )
      end if
      
      ! increase receive buffers
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
      
      ! increase send buffer
      xlen = ifull*ky
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
   
      ! allocated shared memory for internal node buffer
#ifdef usempi3
      shsize(1) = ipan
      shsize(2) = jpan
      shsize(3) = npan
      shsize(4) = kx
      shsize(5) = specmap_indxlen
      call ccmpi_allocshdata(nodepack,shsize(1:5),nodepack_win)
#else
      allocate( nodepack(ipan,jpan,npan,kx,specmap_indxlen) )
#endif
   
   end subroutine allocateglobalpack
   
   subroutine deallocateglobalpack
      ! deallocate sparse arrays and shared memory
      integer :: ncount, w, iproc, ipak, jpak, n
      integer :: ipoff, jpoff, npoff

      ncount = size(specmap_ext)
      do w = 1,ncount
         iproc = specmap_ext(w)
         call proc_region_face(iproc,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan)
         ! Global indices are i+ipoff, j+jpoff, n-npoff
         ipak = ipoff/ipan
         jpak = jpoff/jpan
         do n = 1,npan
            deallocate(globalpack(ipak,jpak,n-npoff)%localdata)
         end do
      end do
      
      deallocate(globalpack)
   
#ifdef usempi3
      call ccmpi_freeshdata(nodepack_win)
#else
      deallocate( nodepack )
#endif
      nullify( nodepack )
      
      deallocate( specmap_ext, specmap_send, specmap_recv, specmap_req, specmap_indx )
   
   end subroutine deallocateglobalpack
   
   subroutine ccglobal_posneg2(array, delpos, delneg)
      real, intent(in), dimension(ifull) :: array
      real, intent(out) :: delpos, delneg
      real, dimension(ifull,1,1) :: array_l
      real, dimension(1) :: delpos_l, delneg_l
      real, dimension(1) :: dsig
      
      dsig(1) = 1.
      array_l(1:ifull,1,1) = array(1:ifull)
      call ccglobal_posneg4(array_l, delpos_l, delneg_l, dsig)
      delpos = delpos_l(1)
      delneg = delneg_l(1)

   end subroutine ccglobal_posneg2
    
   subroutine ccglobal_posneg3(array, delpos, delneg, dsig)
      real, intent(in), dimension(:,:) :: array
      real, intent(in), dimension(:) :: dsig
      real, intent(out) :: delpos, delneg
      real, dimension(ifull,size(array,2),1) :: array_l
      real, dimension(1) :: delpos_l, delneg_l
      
      array_l(1:ifull,:,1) = array(1:ifull,:)
      call ccglobal_posneg4(array_l, delpos_l, delneg_l, dsig)
      delpos = delpos_l(1)
      delneg = delneg_l(1)

   end subroutine ccglobal_posneg3
   
   subroutine ccglobal_posneg3o(array, delpos, delneg, dsig)
      real, intent(in), dimension(:,:) :: array
      real, intent(in), dimension(:,:) :: dsig
      real, intent(out) :: delpos, delneg
      real, dimension(ifull,size(array,2),1) :: array_l
      real, dimension(1) :: delpos_l, delneg_l
      
      array_l(1:ifull,:,1) = array(1:ifull,:)
      call ccglobal_posneg4o(array_l, delpos_l, delneg_l, dsig)
      delpos = delpos_l(1)
      delneg = delneg_l(1)

   end subroutine ccglobal_posneg3o

   subroutine ccglobal_posneg4(array, delpos, delneg, dsig)
      ! Calculate global sums of positive and negative values of array
      use sumdd_m
      use xyzinfo_m
      real, intent(in), dimension(:,:,:) :: array
      real, intent(in), dimension(:) :: dsig
      real, intent(out), dimension(:) :: delpos, delneg
      real, dimension(ifull,2*size(array,2)*size(array,3)) :: tmparr
      integer :: i, k, kx, ntr
      integer(kind=4) :: ierr, mnum, lcomm
      complex, dimension(2*size(array,3)) :: local_sum, global_sum
      complex, dimension(size(array,2),2*size(array,3)) :: tmparr_k
      complex, dimension(2*size(array,2)*size(array,3)) :: local_sum_k

      kx  = size(array,2)
      ntr = size(array,3)
      local_sum_k(1:2*kx*ntr) = cmplx(0., 0.)
      do i = 1,ntr
         do k = 1,kx
            tmparr(1:ifull,k+2*kx*(i-1)) = max(0.,abs(dsig(k))*array(1:ifull,k,i)*wts(1:ifull))
            tmparr(1:ifull,k+kx+2*kx*(i-1)) = min(0.,abs(dsig(k))*array(1:ifull,k,i)*wts(1:ifull))
         end do ! k loop
      end do
      call drpdr_local_v(tmparr,local_sum_k)
      local_sum(1:2*ntr) = cmplx(0.,0.)      
      do i = 1,ntr
         tmparr_k(1:kx,i) = local_sum_k(1+2*kx*(i-1):kx+2*kx*(i-1))
         tmparr_k(1:kx,i+ntr) = local_sum_k(kx+1+2*kx*(i-1):2*kx+2*kx*(i-1))
      end do
      call drpdr_v(tmparr_k,local_sum)
      mnum = 2*ntr
      global_sum(1:2*ntr) = cmplx(0.,0.)
      lcomm = comm_world
      call START_LOG(allreduce_begin)
#ifdef i8r8
      call MPI_Allreduce( local_sum, global_sum, mnum, MPI_DOUBLE_COMPLEX, MPI_SUMDR, lcomm, ierr )
#else
      call MPI_Allreduce( local_sum, global_sum, mnum, MPI_COMPLEX, MPI_SUMDR, lcomm, ierr )
#endif
      call END_LOG(allreduce_end)
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
      real, dimension(ifull,2*size(array,2)*size(array,3)) :: tmparr
      integer :: i, k, kx, ntr
      integer(kind=4) :: ierr, mnum, lcomm
      complex, dimension(2*size(array,3)) :: local_sum, global_sum
      complex, dimension(size(array,2),2*size(array,3)) :: tmparr_k
      complex, dimension(2*size(array,2)*size(array,3)) :: local_sum_k

      kx  = size(array,2)
      ntr = size(array,3)
      local_sum_k(1:2*kx*ntr) = cmplx(0., 0.)
      do i = 1,ntr
         do k = 1,kx
            tmparr(1:ifull,k+2*kx*(i-1)) = max(0.,abs(dsig(1:ifull,k))*array(1:ifull,k,i)*wts(1:ifull))
            tmparr(1:ifull,k+kx+2*kx*(i-1)) = min(0.,abs(dsig(1:ifull,k))*array(1:ifull,k,i)*wts(1:ifull))
         end do ! k loop
      end do
      call drpdr_local_v(tmparr,local_sum_k)
      local_sum(1:2*ntr) = cmplx(0.,0.)
      do i = 1,ntr
         tmparr_k(1:kx,i) = local_sum_k(1+2*kx*(i-1):kx+2*kx*(i-1))
         tmparr_k(1:kx,i+ntr) = local_sum_k(kx+1+2*kx*(i-1):2*kx+2*kx*(i-1))
      end do
      call drpdr_v(tmparr_k,local_sum)
      mnum = 2*ntr
      global_sum(1:2*ntr) = cmplx(0.,0.)
      lcomm = comm_world
      call START_LOG(allreduce_begin)
#ifdef i8r8
      call MPI_Allreduce( local_sum, global_sum, mnum, MPI_DOUBLE_COMPLEX, MPI_SUMDR, lcomm, ierr )
#else
      call MPI_Allreduce( local_sum, global_sum, mnum, MPI_COMPLEX, MPI_SUMDR, lcomm, ierr )
#endif
      call END_LOG(allreduce_end)
      delpos(1:ntr) = real(global_sum(1:ntr))
      delneg(1:ntr) = real(global_sum(ntr+1:2*ntr))

   end subroutine ccglobal_posneg4o

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
      integer :: kk
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
      integer :: kk

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

end module cc_mpi_collective

