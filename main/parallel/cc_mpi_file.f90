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

! Subroutines for MPI communication between processes with input files    
    
module cc_mpi_file

   use cc_mpi_common

   implicit none

   private

   integer, save, public :: pipan, pjpan, pnpan                            ! decomposition parameters file window
   integer, save, public :: pil_g, pjl_g, pka_g, pko_g                     ! decomposition parameters file window
   integer, save, public :: fnproc, fnresid, fncount, mynproc              ! number and decomposition of input files
   integer, allocatable, dimension(:), save, public :: pnoff               ! file window panel offset
   integer, allocatable, dimension(:,:), save, public :: pioff, pjoff      ! file window coordinate offset
   integer(kind=4), allocatable, dimension(:), save, public ::              &
      filemap_req, filemap_qmod                                            ! file map requested for onthefly
   integer(kind=4), allocatable, dimension(:), save, public ::              &
      filemap_recv, filemap_rmod                                           ! file map received for onthefly
   integer(kind=4), allocatable, dimension(:), save, public ::              &
      filemap_send, filemap_smod                                           ! file map sent for onthefly
   integer, allocatable, dimension(:,:), save, public :: filemap_indx      ! index for file map messages
   integer, save, public :: filemap_indxlen
   real, dimension(:,:,:,:,:), pointer, save, private :: nodefile          ! node buffer for file map
   integer, save, private :: nodefile_win
   integer, dimension(3), save, private :: nodefilesave_win
   integer, save, private :: nodefile_count = 0
   
   integer, allocatable, dimension(:), save, private :: fileneighlist      ! list of file neighbour processes
   integer, save, public :: fileneighnum                                   ! number of file neighbours
  
   public :: ccmpi_filewinget, ccmpi_filebounds_setup,                      &
             ccmpi_filebounds, ccmpi_filedistribute, procarray,             &
             ccmpi_filewininit, ccmpi_filewinfinalize,                      &
             ccmpi_filewinfinalize_exit
   
   interface ccmpi_filebounds
      module procedure ccmpi_filebounds2, ccmpi_filebounds3
   end interface ccmpi_filebounds
   interface ccmpi_filedistribute
      module procedure host_filedistribute2, proc_filedistribute2
      module procedure host_filedistribute3, proc_filedistribute3
   end interface ccmpi_filedistribute
   interface ccmpi_filewinget
     module procedure ccmpi_filewinget2, ccmpi_filewinget3
   end interface

   ! File IO
   type filebounds_info
      integer, dimension(:,:), allocatable :: send_list
      integer, dimension(:,:), allocatable :: request_list
      integer, dimension(:,:), allocatable :: unpack_list
      integer :: len, rlen, slen, rlenx, slenx
   end type filebounds_info
   
   type(filebounds_info), allocatable, dimension(:), save, private :: filebnds

contains

   subroutine ccmpi_filewinget2(abuf,sinp)   
      real, dimension(:), intent(in) :: sinp
      real, dimension(-1:pipan+2,-1:pjpan+2,pnpan,size(filemap_req)), intent(out) :: abuf 
      real, dimension(size(sinp),1) :: sinp_l
      real, dimension(-1:pipan+2,-1:pjpan+2,pnpan,size(filemap_req),1) :: abuf_l

      sinp_l(:,1) = sinp(:)
      call ccmpi_filewinget3(abuf_l,sinp_l)
      abuf(:,:,:,:) = abuf_l(:,:,:,:,1)
      
   end subroutine ccmpi_filewinget2

   subroutine ccmpi_filewinget3(abuf,sinp)
      ! broadcast file data to processes using shared memory for efficent data transfer
      real, dimension(:,:), intent(in) :: sinp
      real, dimension(-1:pipan+2,-1:pjpan+2,1:pnpan,1:size(filemap_req),1:size(sinp,2)), intent(out) :: abuf
      real, dimension(pipan*pjpan*pnpan,size(sinp,2),size(filemap_recv)) :: bbuf
      real, dimension(pipan*pjpan*pnpan,size(sinp,2),mynproc) :: cbuf
      integer :: n, w, nlen, kx, cc, ipf
      integer :: rcount, jproc, kproc
      integer :: k, ip, no, ca, cb
      integer(kind=4) :: lcomm, ldone, lsize, ierr, lproc, mproc
      integer(kind=4) :: itag = 51
      integer(kind=4), dimension(size(filemap_recv)+size(filemap_send)) :: i_req, donelist
      integer(kind=4), dimension(size(filemap_recv)) :: i_list

      kx = size(sinp,2)
      nlen = pipan*pjpan*pnpan
      lsize = nlen*kx
      lcomm = comm_world

      !     Set up the buffers to recv
      nreq = 0
      do w = 1,size(filemap_recv)
         ipf = filemap_rmod(w)
         itag = 51 + ipf
         nreq  = nreq + 1
         i_list(nreq) = w
         lproc = filemap_recv(w)
#ifdef i8r8
         call MPI_IRecv( bbuf(:,:,w), lsize, MPI_DOUBLE_PRECISION, lproc, itag, lcomm, i_req(nreq), ierr )
#else
         call MPI_IRecv( bbuf(:,:,w), lsize, MPI_REAL, lproc, itag, lcomm, i_req(nreq), ierr )
#endif
      end do
      rreq = nreq
      
      !     Set up the buffers to send
      do ipf = 0,mynproc-1
         cc = nlen*ipf 
         cbuf(1:nlen,1:kx,ipf+1) = sinp(1+cc:nlen+cc,1:kx)
      end do
      do w = 1,size(filemap_send)
         ipf = filemap_smod(w)
         itag = 51 + ipf
         nreq  = nreq + 1
         lproc = filemap_send(w)
#ifdef i8r8
         call MPI_ISend( cbuf(:,:,ipf+1), lsize, MPI_DOUBLE_PRECISION, lproc, itag, lcomm, i_req(nreq), ierr )
#else
         call MPI_ISend( cbuf(:,:,ipf+1), lsize, MPI_REAL, lproc, itag, lcomm, i_req(nreq), ierr )
#endif
      end do

      ! Unpack incomming messages
      rcount = nreq
      do while ( rcount > 0 )
         call START_LOG(mpiwait_begin) 
         call MPI_Waitsome( nreq, i_req, ldone, donelist, MPI_STATUSES_IGNORE, ierr )
         call END_LOG(mpiwait_end)
         rcount = rcount - ldone
         do jproc = 1,ldone
            mproc = donelist(jproc)
            if ( mproc <= rreq ) then
               w = i_list(mproc)
               kproc = filemap_indx(filemap_recv(w),filemap_rmod(w))
               do k = 1,kx
                  do n = 0,pnpan-1
                     cc = n*pipan*pjpan
                     nodefile(:,:,n+1,k,kproc) = reshape( bbuf(1+cc:pipan*pjpan+cc,k,w), (/ pipan, pjpan /) )
                   end do  
               end do
            end if
         end do
      end do

      call START_LOG(mpibarrier_begin) 
      lcomm = comm_node
      call MPI_Barrier( lcomm, ierr )
      call END_LOG(mpibarrier_end)

      do w = 1,size(filemap_req)
         kproc = filemap_indx(filemap_req(w),filemap_qmod(w))
         do k = 1,kx
            abuf(1:pipan,1:pjpan,1:pnpan,w,k) = nodefile(1:pipan,1:pjpan,1:pnpan,k,kproc)
         end do
      end do  

      do k = 1,kx
         call abufpanelbounds(abuf(:,:,:,:,k))
      end do    
      
      call START_LOG(mpibarrier_begin) 
      lcomm = comm_node
      call MPI_Barrier( lcomm, ierr )
      call END_LOG(mpibarrier_end)

   end subroutine ccmpi_filewinget3
   
   subroutine abufpanelbounds(abuf)
      real, dimension(-1:pipan+2,-1:pjpan+2,1:pnpan,size(filemap_req)), intent(inout) :: abuf
      real, dimension(-1:pil_g+2,-1:pil_g+2,0:npanels) :: sx_l
      integer :: w, ip, n, no, ca, cb
      integer :: i, n_w, n_e, n_n, n_s
      
      sx_l = 0.
      
      do w = 1,size(filemap_req)
         ip = filemap_req(w) + filemap_qmod(w)*fnresid
         do n = 1,pnpan
            no = n - pnoff(ip)
            ca = pioff(ip,no)
            cb = pjoff(ip,no)
            sx_l(1+ca:pipan+ca,1+cb:pjpan+cb,no) = abuf(1:pipan,1:pjpan,n,w)
         end do
      end do
      
      do n = 0,npanels
         if ( mod(n,2)==0 ) then
            n_w = mod(n+5, 6)
            n_e = mod(n+2, 6)
            n_n = mod(n+1, 6)
            n_s = mod(n+4, 6)
            do i = 1,pil_g
               sx_l(-1,i,n)      = sx_l(pil_g-1,i,n_w)
               sx_l(0,i,n)       = sx_l(pil_g,i,n_w)
               sx_l(pil_g+1,i,n) = sx_l(pil_g+1-i,1,n_e)
               sx_l(pil_g+2,i,n) = sx_l(pil_g+1-i,2,n_e)
               sx_l(i,-1,n)      = sx_l(pil_g-1,pil_g+1-i,n_s)
               sx_l(i,0,n)       = sx_l(pil_g,pil_g+1-i,n_s)
               sx_l(i,pil_g+1,n) = sx_l(i,1,n_n)
               sx_l(i,pil_g+2,n) = sx_l(i,2,n_n)
            end do ! i
            sx_l(0,0,n)             = sx_l(pil_g,1,n_w)        ! ws
            sx_l(-1,0,n)            = sx_l(pil_g,2,n_w)        ! wws
            sx_l(0,-1,n)            = sx_l(pil_g,pil_g-1,n_s)  ! wss
            sx_l(pil_g+1,0,n)       = sx_l(pil_g,1,n_e)        ! es  
            sx_l(pil_g+2,0,n)       = sx_l(pil_g-1,1,n_e)      ! ees 
            sx_l(pil_g+1,-1,n)      = sx_l(pil_g,2,n_e)        ! ess        
            sx_l(0,pil_g+1,n)       = sx_l(pil_g,pil_g,n_w)    ! wn  
            sx_l(-1,pil_g+1,n)      = sx_l(pil_g,pil_g-1,n_w)  ! wwn
            sx_l(0,pil_g+2,n)       = sx_l(pil_g-1,pil_g,n_w)  ! wnn
            sx_l(pil_g+1,pil_g+1,n) = sx_l(1,1,n_e)            ! en  
            sx_l(pil_g+2,pil_g+1,n) = sx_l(2,1,n_e)            ! een  
            sx_l(pil_g+1,pil_g+2,n) = sx_l(1,2,n_e)            ! enn  
         else
            n_w = mod(n+4, 6)
            n_e = mod(n+1, 6)
            n_n = mod(n+2, 6)
            n_s = mod(n+5, 6)
            do i = 1,pil_g
               sx_l(-1,i,n)      = sx_l(pil_g+1-i,pil_g-1,n_w)  
               sx_l(0,i,n)       = sx_l(pil_g+1-i,pil_g,n_w)
               sx_l(pil_g+1,i,n) = sx_l(1,i,n_e)
               sx_l(pil_g+2,i,n) = sx_l(2,i,n_e)
               sx_l(i,-1,n)      = sx_l(i,pil_g-1,n_s)
               sx_l(i,0,n)       = sx_l(i,pil_g,n_s)
               sx_l(i,pil_g+1,n) = sx_l(1,pil_g+1-i,n_n)
               sx_l(i,pil_g+2,n) = sx_l(2,pil_g+1-i,n_n)
            end do ! i
            sx_l(0,0,n)             = sx_l(pil_g,pil_g,n_w)   ! ws
            sx_l(-1,0,n)            = sx_l(pil_g-1,pil_g,n_w) ! wws
            sx_l(0,-1,n)            = sx_l(2,pil_g,n_s)       ! wss
            sx_l(pil_g+1,0,n)       = sx_l(1,1,n_e)           ! es
            sx_l(pil_g+2,0,n)       = sx_l(1,2,n_e)           ! ees
            sx_l(pil_g+1,-1,n)      = sx_l(2,1,n_e)           ! ess
            sx_l(0,pil_g+1,n)       = sx_l(1,pil_g,n_w)       ! wn       
            sx_l(-1,pil_g+1,n)      = sx_l(2,pil_g,n_w)       ! wwn   
            sx_l(0,pil_g+2,n)       = sx_l(1,pil_g-1,n_w)     ! wnn
            sx_l(pil_g+1,pil_g+1,n) = sx_l(1,pil_g,n_e)       ! en  
            sx_l(pil_g+2,pil_g+1,n) = sx_l(1,pil_g-1,n_e)     ! een  
            sx_l(pil_g+1,pil_g+2,n) = sx_l(2,pil_g,n_e)       ! enn  
         end if   ! mod(n,2)==0 ..else..
      end do       ! n loop

      do w = 1,size(filemap_req)
         ip = filemap_req(w) + filemap_qmod(w)*fnresid
         do n = 1,pnpan
            no = n - pnoff(ip)
            ca = pioff(ip,no)
            cb = pjoff(ip,no)
            abuf(-1:pipan+2,-1:pjpan+2,n,w) = sx_l(-1+ca:pipan+2+ca,-1+cb:pjpan+2+cb,no) 
         end do
      end do
      
   end subroutine abufpanelbounds

   subroutine ccmpi_filewininit(kblock)
      integer, intent(in) :: kblock
      integer, dimension(5) :: shsize
      
      ! allocated shared memory for internal node buffer
#ifdef usempi3
      shsize(1) = pipan
      shsize(2) = pjpan
      shsize(3) = pnpan
      shsize(4) = kblock
      shsize(5) = filemap_indxlen
      call ccmpi_allocshdata(nodefile,shsize(1:5),nodefile_win)
      nodefile_count = nodefile_count + 1
#else
      allocate( nodefile(pipan,pjpan,pnpan,kblock,filemap_indxlen) )
#endif
   
   end subroutine ccmpi_filewininit
   
   subroutine ccmpi_filewinfinalize
      
#ifdef usempi3
      ! Previously deallocating memory triggered bugs in the MPI library
      ! Here we will leave the memory allocated for now and deallocate
      ! just before mpi_finalize with ccmpi_filewinfinalize_exit
      if ( nodefile_count > 0 ) then
         if ( nodefile_count <= size(nodefilesave_win) ) then
            nodefilesave_win(nodefile_count) = nodefile_win
         end if   
         !call ccmpi_freeshdata(nodefile_win)         
      end if   
#else
      if ( associated( nodefile ) ) then
         deallocate( nodefile )
      end if   
#endif
      nullify(nodefile)
   
   end subroutine ccmpi_filewinfinalize
   
   subroutine ccmpi_filewinfinalize_exit
      integer :: i
   
#ifdef usempi3
      do i = 1,nodefile_count
         call ccmpi_freeshdata(nodefilesave_win(i))
      end do   
#endif
   
   end subroutine ccmpi_filewinfinalize_exit
   
   subroutine ccmpi_filebounds_setup(comm_ip)
      integer, intent(in) :: comm_ip
      integer, dimension(:,:), allocatable, save :: dumi
      integer :: ipf, n, i, j, iq, ncount, ca, cb, no, ip
      integer :: filemaxbuflen, xlen
      integer :: iproc, jproc, iloc, jloc, nloc, floc
      integer(kind=4) :: lproc, lcomm, llen, ierr
      integer(kind=4) :: itag=42
      logical, save :: fileallocate = .false.
     
      if ( myid>=fnresid ) return
      
      lcomm = comm_ip
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
      do n = 0,fnresid-1
         filebnds(n)%len = 0
         filebnds(n)%rlen = 0
         filebnds(n)%slen = 0
         filebnds(n)%rlenx = 0
         filebnds(n)%slenx = 0
      end do   
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
      
      do n = 0,fnresid-1
         filebnds(n)%rlen = filebnds(n)%rlenx
      end do   
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
 
      ! identify neighbour processes
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

      ! Now, for each process send the length of points I want.
      allocate( dumi(4,fileneighnum) )
      nreq = 0
      do jproc = 1,fileneighnum
         iproc = fileneighlist(jproc) ! Recv from
         nreq = nreq + 1
         lproc = iproc
#ifdef i8r8
         call MPI_IRecv( dumi(1:2,jproc), 2_4, MPI_INTEGER8, lproc, itag, lcomm, ireq(nreq), ierr )
#else
         call MPI_IRecv( dumi(1:2,jproc), 2_4, MPI_INTEGER, lproc, itag, lcomm, ireq(nreq), ierr )
#endif
      end do
      do jproc = fileneighnum,1,-1
         iproc = fileneighlist(jproc)
         nreq = nreq + 1
         lproc = iproc
         dumi(3,jproc) = filebnds(iproc)%rlen
         dumi(4,jproc) = filebnds(iproc)%rlenx
#ifdef i8r8
         call MPI_ISend( dumi(3:4,jproc), 2_4, MPI_INTEGER8, lproc, itag, lcomm, ireq(nreq), ierr )
#else
         call MPI_ISend( dumi(3:4,jproc), 2_4, MPI_INTEGER, lproc, itag, lcomm, ireq(nreq), ierr )
#endif
      end do
      call MPI_Waitall( nreq, ireq, MPI_STATUSES_IGNORE, ierr )
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
#ifdef i8r8
            call MPI_IRecv( filebnds(iproc)%send_list(:,i), llen, &
                 MPI_INTEGER8, lproc, itag, lcomm, ireq(nreq), ierr )
#else
            call MPI_IRecv( filebnds(iproc)%send_list(:,i), llen, &
                 MPI_INTEGER, lproc, itag, lcomm, ireq(nreq), ierr )
#endif
         end do
         do jproc = fileneighnum,1,-1
            iproc = fileneighlist(jproc)  ! Send to
            ! Send list of requests
            nreq = nreq + 1
            llen = filebnds(iproc)%rlenx
            lproc = iproc
#ifdef i8r8
            call MPI_ISend( filebnds(iproc)%request_list(:,i), llen, &
                 MPI_INTEGER8, lproc, itag, lcomm, ireq(nreq), ierr )
#else
            call MPI_ISend( filebnds(iproc)%request_list(:,i), llen, &
                 MPI_INTEGER, lproc, itag, lcomm, ireq(nreq), ierr )
#endif
         end do
         call MPI_Waitall( nreq, ireq, MPI_STATUSES_IGNORE, ierr )
      end do

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
         lproc = iproc
         if ( bnds(lproc)%rbuflen < xlen ) then
            if ( allocated(bnds(lproc)%rbuf) ) then
               deallocate( bnds(lproc)%rbuf )
               deallocate( bnds(lproc)%r8buf )
            end if
            allocate( bnds(lproc)%rbuf(xlen) )
            allocate( bnds(lproc)%r8buf(xlen) )
            bnds(lproc)%rbuflen = xlen
         end if
         xlen = maxvertlen*filebnds(iproc)%slenx
         if ( bnds(lproc)%sbuflen < xlen ) then
            if ( allocated(bnds(lproc)%sbuf) ) then
               deallocate( bnds(lproc)%sbuf )
               deallocate( bnds(lproc)%s8buf )
            end if
            allocate( bnds(lproc)%sbuf(xlen) )
            allocate( bnds(lproc)%s8buf(xlen) )
            bnds(lproc)%sbuflen = xlen
         end if
      end do
      
   end subroutine ccmpi_filebounds_setup
   
   subroutine check_filebnds_alloc(iproc,filemaxbuflen)
      ! allocate memory for bounds with processes reading files
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

   subroutine ccmpi_filebounds2(sdat,comm_ip,corner)
      real, dimension(0:pipan+1,0:pjpan+1,pnpan,1:fncount), intent(inout) :: sdat
      integer, intent(in) :: comm_ip
      logical, intent(in), optional :: corner
      real, dimension(0:pipan+1,0:pjpan+1,pnpan,1:fncount,1) :: sdat_l
      logical :: corner_l
      
      corner_l = .true.
      if ( present(corner) ) then
         corner_l = corner
      end if
      sdat_l(:,:,:,:,1) = sdat(:,:,:,:)
      call ccmpi_filebounds3( sdat_l, comm_ip, corner=corner_l )
      sdat(:,:,:,:) = sdat_l(:,:,:,:,1)

   end subroutine ccmpi_filebounds2
   
   subroutine ccmpi_filebounds3(sdat,comm_ip,corner)
      ! update halo for processes part of the file communicator
      real, dimension(0:,0:,1:,1:,1:), intent(inout) :: sdat
      integer, intent(in) :: comm_ip
      logical, intent(in), optional :: corner
      integer :: iproc, iq, kx, send_len, k
      integer :: myrlen, jproc, mproc, rcount
      integer, dimension(fileneighnum) :: rslen, sslen
      integer(kind=4) :: llen, lproc, ierr, lcomm, ldone
      integer(kind=4) :: itag=41
      integer(kind=4), dimension(2*fileneighnum) :: donelist
      logical :: extra
      
      kx = size(sdat,5)
      lcomm = comm_ip

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
      myrlen = filebnds(myid)%rlenx

      !     Set up the buffers to send and recv
      nreq = 0
      do iproc = 1,fileneighnum
         if ( rslen(iproc)>0 ) then
            llen = rslen(iproc)*kx
            lproc = fileneighlist(iproc)  ! Recv from
            nreq = nreq + 1
            rlist(nreq) = iproc
#ifdef i8r8
            call MPI_IRecv( bnds(lproc)%rbuf, llen, MPI_DOUBLE_PRECISION, lproc, itag, lcomm, ireq(nreq), ierr )
#else
            call MPI_IRecv( bnds(lproc)%rbuf, llen, MPI_REAL, lproc, itag, lcomm, ireq(nreq), ierr )
#endif
         end if   
      end do
      rreq = nreq
      do iproc = fileneighnum,1,-1
         if ( sslen(iproc) > 0 ) then 
            send_len = sslen(iproc)
            llen = send_len*kx
            lproc = fileneighlist(iproc)  ! Send to
            do k = 1,kx
               do iq = 1,send_len
                  bnds(lproc)%sbuf(iq+(k-1)*send_len) =                              &
                     sdat(filebnds(lproc)%send_list(iq,1),filebnds(lproc)%send_list(iq,2),    &
                          filebnds(lproc)%send_list(iq,3),filebnds(lproc)%send_list(iq,4),k)
               end do
            end do   
            nreq = nreq + 1
#ifdef i8r8
            call MPI_ISend( bnds(lproc)%sbuf, llen, MPI_DOUBLE_PRECISION, lproc, itag, lcomm, ireq(nreq), ierr )
#else
            call MPI_ISend( bnds(lproc)%sbuf, llen, MPI_REAL, lproc, itag, lcomm, ireq(nreq), ierr )
#endif
         end if   
      end do

      ! See if there are any points on my own process that need
      ! to be fixed up.
      do k = 1,kx
         do iq = 1,myrlen
            ! request_list is same as send_list in this case
            sdat(filebnds(myid)%unpack_list(iq,1),filebnds(myid)%unpack_list(iq,2),     &
                 filebnds(myid)%unpack_list(iq,3),filebnds(myid)%unpack_list(iq,4),k) = &
            sdat(filebnds(myid)%request_list(iq,1),filebnds(myid)%request_list(iq,2),   &
                 filebnds(myid)%request_list(iq,3),filebnds(myid)%request_list(iq,4),k)
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
               lproc = fileneighlist(iproc)
               do k = 1,kx
                  do iq = 1,rslen(iproc)
                     sdat(filebnds(lproc)%unpack_list(iq,1),filebnds(lproc)%unpack_list(iq,2),     &
                          filebnds(lproc)%unpack_list(iq,3),filebnds(lproc)%unpack_list(iq,4),k) = &
                     bnds(lproc)%rbuf(iq+(k-1)*rslen(iproc))
                  end do
               end do
            end if   
         end do
      end do
      
   end subroutine ccmpi_filebounds3
   
   subroutine host_filedistribute2(af,a1,comm)
      integer, intent(in) :: comm
      real, dimension(pipan*pjpan*pnpan*fncount), intent(out) :: af
      real, dimension(pil_g*pjl_g), intent(in) :: a1
      real, dimension(pipan*pjpan*pnpan*fncount,1) :: af_l
      real, dimension(pil_g*pjl_g,1) :: a1_l
      
      a1_l(:,1) = a1(:)
      call host_filedistribute3(af_l,a1_l,comm)
      af(:) = af_l(:,1)
      
   end subroutine host_filedistribute2

   subroutine proc_filedistribute2(af,comm)
      integer, intent(in) :: comm
      real, dimension(pipan*pjpan*pnpan*fncount), intent(out) :: af
      real, dimension(pipan*pjpan*pnpan*fncount,1) :: af_l
      
      call proc_filedistribute3(af_l,comm)
      af(:) = af_l(:,1)
      
   end subroutine proc_filedistribute2

   subroutine host_filedistribute3(af,a1,comm)
      ! Convert standard 1D arrays to face form and distribute to processes
      integer, intent(in) :: comm
      real, dimension(:,:), intent(out) :: af
      real, dimension(:,:), intent(in) :: a1
      real, dimension(size(af,1),size(af,2),0:fnresid-1) :: sbuf
      integer :: j, n, iproc, k, ipf, ip
      integer :: kx, ca, cb
      integer(kind=4) :: ierr, lsize, lcomm
      
      kx = size(af,2)
      
      ! map array in order of process rank
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
      call START_LOG(scatter_begin)
#ifdef i8r8
      call MPI_Scatter( sbuf, lsize, MPI_DOUBLE_PRECISION, af, lsize, MPI_DOUBLE_PRECISION, 0_4, lcomm, ierr )
#else
      call MPI_Scatter( sbuf, lsize, MPI_REAL, af, lsize, MPI_REAL, 0_4, lcomm, ierr )
#endif
      call END_LOG(scatter_end)
      
   end subroutine host_filedistribute3

   subroutine proc_filedistribute3(af,comm)
      ! Convert standard 1D arrays to face form and distribute to processes
      integer, intent(in) :: comm
      real, dimension(:,:), intent(out) :: af
      real, dimension(1,1,1) :: sbuf
      integer :: kx
      integer(kind=4) :: ierr, lsize, lcomm
      
      kx = size(af,2)
      lsize = pipan*pjpan*pnpan*fncount*kx
      lcomm = comm
      call START_LOG(scatter_begin)
#ifdef i8r8
      call MPI_Scatter( sbuf, lsize, MPI_DOUBLE_PRECISION, af, lsize, MPI_DOUBLE_PRECISION, 0_4, lcomm, ierr )
#else
      call MPI_Scatter( sbuf, lsize, MPI_REAL, af, lsize, MPI_REAL, 0_4, lcomm, ierr )
#endif
      call END_LOG(scatter_end)
      
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
      
      ! determine process that owns the grid point
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
      integer, intent(in) :: floc, iproc, ik
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

end module cc_mpi_file

