module cc_mpi
   implicit none
   private
   include 'newmpar.h'
   include 'mpif.h'

   integer, public :: myid ! Processor number
   ! These are just defined here to minimise changes from the ccshallow code.
   ! These could be parameters
   integer, public :: ipan, jpan !, npan
   integer, public :: ioff, joff, noff
   integer, private :: nxproc, nyproc

   integer, public, save, dimension(il_g,il_g,0:npanels) :: fproc
   integer, public, save, dimension(ifull_g) :: qproc


   ! These only need to be module variables for the overlap case
   ! when they have to keep values between boundsa and boundsb.
   ! For now leave them here.
   integer, dimension(2*nproc), save, private :: ireq
   integer, save, private :: nreq

   public :: bounds, boundsuv, ccmpi_setup, ccmpi_distribute, ccmpi_gather, &
             indp, indg, deptsync, intssync
   private :: indv_mpi, ccmpi_distribute2, ccmpi_distribute3, &
              ccmpi_gather2, ccmpi_gather3, checksize
   interface ccmpi_gather
      module procedure ccmpi_gather2, ccmpi_gather3
   end interface
   interface ccmpi_distribute
      module procedure ccmpi_distribute2, ccmpi_distribute3
   end interface
   interface bounds
      module procedure bounds2, bounds3
   end interface
   interface boundsuv
      module procedure boundsuv2, boundsuv3
   end interface

   ! Define neighbouring faces
   integer, parameter, private, dimension(0:npanels) ::    &
           n_e = (/ 2, 2, 4, 4, 0, 0 /), &
           n_w = (/ 5, 5, 1, 1, 3, 3 /), &
           n_n = (/ 1, 3, 3, 5, 5, 1 /), &
           n_s = (/ 4, 0, 0, 2, 2, 4 /)
   ! Do directions need to be swapped
   logical, parameter, private, dimension(0:npanels) ::    &
           swap_e = (/ .true., .false., .true., .false., .true., .false. /), &
           swap_w = (/ .false., .true., .false., .true., .false., .true. /), &
           swap_n = (/ .false., .true., .false., .true., .false., .true. /), &
           swap_s = (/ .true., .false., .true., .false., .true., .false. /)

   type bounds_info
      real, dimension(:), pointer :: sbuf, rbuf
      integer, dimension(:), pointer :: request_list
      integer, dimension(:), pointer :: send_list
      integer, dimension(:), pointer :: unpack_list
      integer, dimension(:), pointer :: request_list_uv
      integer, dimension(:), pointer :: send_list_uv
      integer, dimension(:), pointer :: unpack_list_uv
      ! Flag for whether u and v need to be swapped
      logical, dimension(:), pointer :: uv_swap, send_swap
      ! Number of points for each processor. Also double row versions.
      integer :: slen, rlen, slen2, rlen2
      ! Number of points for each processor.
      integer :: slen_uv, rlen_uv, slen2_uv, rlen2_uv
      integer :: len
   end type bounds_info

   type(bounds_info), dimension(0:nproc-1), save :: bnds

   ! Flag whether processor region edge is a face edge.
   logical, public, save :: edge_w, edge_n, edge_s, edge_e

   ! Off processor departure points
   integer, private, save :: maxsize
   real, dimension(:,:,:), allocatable, public, save :: dpoints
   integer, dimension(:,:,:), allocatable, public, save :: dindex
   real, dimension(:,:),  allocatable, public, save :: sextra
   ! Number of points for each processor.
   integer, dimension(:), allocatable, public, save :: dslen, drlen

   logical, public, save :: mydiag ! True if diagnostic point id, jd is in my region

#if defined(vampir)  || defined(mpilog)
   integer, public, save :: bounds_begin, bounds_end
   integer, public, save :: boundsa_begin, boundsa_end
   integer, public, save :: boundsb_begin, boundsb_end
   integer, public, save :: boundsuv_begin, boundsuv_end
   integer, public, save :: ints_begin, ints_end
   integer, public, save :: intsbl_begin, intsbl_end
   integer, public, save :: nonlin_begin, nonlin_end
   integer, public, save :: helm_begin, helm_end
   integer, public, save :: adjust_begin, adjust_end
   integer, public, save :: upglobal_begin, upglobal_end
   integer, public, save :: depts_begin, depts_end
   integer, public, save :: deptsync_begin, deptsync_end
   integer, public, save :: intssync_begin, intssync_end
   integer, public, save :: stag_begin, stag_end
   integer, public, save :: toij_begin, toij_end
#endif

#ifdef mpilog
   public :: mpe_log_event, mpe_log_get_event_number, mpe_describe_state
   interface 
      function mpe_log_event (event, idata, string) result(val)
         integer, intent(in) :: event, idata
         character(len=*), intent(in) :: string
         integer :: val
      end function mpe_log_event
      function mpe_log_get_event_number() result(val)
         integer :: val
      end function mpe_log_get_event_number
      function mpe_describe_state(start,finish, name, color) result(val)
         integer, intent(in) :: start, finish
         character(len=*), intent(in) :: name, color
         integer :: val
      end function mpe_describe_state
   end interface
#endif

contains

   subroutine ccmpi_setup()
      include 'xyzinfo.h'
      include 'xyzinfo_g.h'
      include 'map.h'
      include 'map_g.h'
      include 'vecsuv.h'
      include 'vecsuv_g.h'
      include 'latlong.h'
      include 'latlong_g.h'

      call proc_setup(il,npanels,ifull)

      print*, "Grid", npan, ipan, jpan
      print*, "Offsets", myid, ioff, joff, noff

      ! Also do the initialisation for deptsync here
      allocate ( dslen(0:nproc-1), drlen(0:nproc-1) )
      if ( nproc == 1 ) then
         maxsize = 0 ! Not used in this case
      else
         ! 4 rows should be plenty (size is checked anyway).
         maxsize = 4*max(ipan,jpan)*npan*kl
      end if
      ! Off processor departure points
      allocate ( dpoints(4,maxsize,0:nproc-1), sextra(maxsize,0:nproc-1) )
      allocate ( dindex(2,maxsize,0:nproc-1) )
      
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
!!$         call ccmpi_distribute(axu,axu_g)
!!$         call ccmpi_distribute(ayu,ayu_g)
!!$         call ccmpi_distribute(azu,azu_g)
!!$         call ccmpi_distribute(bxv,bxv_g)
!!$         call ccmpi_distribute(byv,byv_g)
!!$         call ccmpi_distribute(bzv,bzv_g)
         call ccmpi_distribute(f,f_g)
         call ccmpi_distribute(fu,fu_g)
         call ccmpi_distribute(fv,fv_g)
         call ccmpi_distribute(x,x_g)
         call ccmpi_distribute(y,y_g)
         call ccmpi_distribute(z,z_g)
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
!!$         call ccmpi_distribute(axu)
!!$         call ccmpi_distribute(ayu)
!!$         call ccmpi_distribute(azu)
!!$         call ccmpi_distribute(bxv)
!!$         call ccmpi_distribute(byv)
!!$         call ccmpi_distribute(bzv)
         call ccmpi_distribute(f)
         call ccmpi_distribute(fu)
         call ccmpi_distribute(fv)
         call ccmpi_distribute(x)
         call ccmpi_distribute(y)
         call ccmpi_distribute(z)
         call ccmpi_distribute(rlatt)
         call ccmpi_distribute(rlongg)
      end if

      call bounds_setup()
      call bounds(em)
      call boundsuv(emu,emv)
      call boundsuv(ax,bx)
      call boundsuv(ay,by)
      call boundsuv(az,bz)
!!$      call bounds(axu)
!!$      call bounds(ayu)
!!$      call bounds(azu)
!!$      call bounds(bxv)
!!$      call bounds(byv)
!!$      call bounds(bzv)

   end subroutine ccmpi_setup

   subroutine ccmpi_distribute2(af,a1)
      ! Convert standard 1D arrays to face form and distribute to processors
      real, dimension(ifull), intent(out) :: af
      real, dimension(ifull_g), intent(in), optional :: a1
      integer :: i, j, n, iq, iqg, n1, n2, iq1, iq2, itag=0, iproc, ierr, count
      integer, dimension(MPI_STATUS_SIZE) :: status
!     Note ipfull = ipan*jpan*npan
      real, dimension(ipan*jpan*npan) :: sbuf
      integer :: npoff, ipoff, jpoff ! Offsets for target
      integer :: slen

      if ( myid == 0 .and. .not. present(a1) ) then
         print*, "Error: ccmpi_distribute argument required on proc 0"
         stop
      end if
      ! Copy internal region
      if ( myid == 0 ) then
         ! First copy own region
         do n=1,npan
            do j=1,jpan
               do i=1,ipan
                  iqg = indg(i,j,n)  ! True global index
                  iq = indp(i,j,n)
                  af(iq) = a1(iqg)
               end do
            end do
         end do
         ! Send appropriate regions to all other processes. In this version
         ! processor regions are no longer necessarily a continuous iq range.
         do iproc=1,nproc-1
            ! Panel range on the target processor
            call proc_region(iproc,ipoff,jpoff,npoff)
!            print*, "TARGET", ipoff, jpoff, npoff
            slen = 0
            do n=1,npan
               do j=1,jpan
                  do i=1,ipan
                     iq = i+ipoff + (j+jpoff-1)*il + (n-npoff)*il*il
                     slen = slen+1
                     sbuf(slen) = a1(iq)
                  end do
               end do
            end do
            call MPI_SSend( sbuf, slen, MPI_REAL, iproc, itag, &
                            MPI_COMM_WORLD, ierr )
         end do
      else ! myid /= 0
         call MPI_Recv( af, ipan*jpan*npan, MPI_REAL, 0, itag, &
                        MPI_COMM_WORLD, status, ierr )
         ! Check that the length is the expected value.
         call MPI_Get_count(status, MPI_REAL, count, ierr)
         if ( count /= ifull ) then
            print*, "Error, wrong length in ccmpi_distribute", myid, ifull, count
            call MPI_Abort(MPI_COMM_WORLD)
         end if

      end if

   end subroutine ccmpi_distribute2

   subroutine ccmpi_distribute3(af,a1)
      ! Convert standard 1D arrays to face form and distribute to processors
      real, dimension(ifull,kl), intent(out) :: af
      real, dimension(ifull_g,kl), intent(in), optional :: a1
      integer :: i, j, n, iq, iqg, n1, n2, iq1, iq2, itag=0, iproc, ierr, count
      integer, dimension(MPI_STATUS_SIZE) :: status
!     Note ipfull = ipan*jpan*npan
      real, dimension(ipan*jpan*npan,kl) :: sbuf
      integer :: npoff, ipoff, jpoff ! Offsets for target
      integer :: slen

      if ( myid == 0 .and. .not. present(a1) ) then
         print*, "Error: ccmpi_distribute argument required on proc 0"
         stop
      end if
      ! Copy internal region
      if ( myid == 0 ) then
         ! First copy own region
         do n=1,npan
            do j=1,jpan
               do i=1,ipan
                  iqg = indg(i,j,n)  ! True global index
                  iq = indp(i,j,n)
                  af(iq,:) = a1(iqg,:)
               end do
            end do
         end do
         ! Send appropriate regions to all other processes. In this version
         ! processor regions are no longer necessarily a continuous iq range.
         do iproc=1,nproc-1
            ! Panel range on the target processor
            call proc_region(iproc,ipoff,jpoff,npoff)
!            print*, "TARGET", ipoff, jpoff, npoff
            slen = 0
            do n=1,npan
               do j=1,jpan
                  do i=1,ipan
                     iq = i+ipoff + (j+jpoff-1)*il + (n-npoff)*il*il
                     slen = slen+kl
                     sbuf(slen,:) = a1(iq,:)
                  end do
               end do
            end do
            call MPI_SSend( sbuf, slen, MPI_REAL, iproc, itag, &
                            MPI_COMM_WORLD, ierr )
         end do
      else ! myid /= 0
         call MPI_Recv( af, ipan*jpan*npan*kl, MPI_REAL, 0, itag, &
                        MPI_COMM_WORLD, status, ierr )
         ! Check that the length is the expected value.
         call MPI_Get_count(status, MPI_REAL, count, ierr)
         if ( count /= ifull*kl ) then
            print*, "Error, wrong length in ccmpi_distribute", myid, ifull*kl, count
            call MPI_Abort(MPI_COMM_WORLD)
         end if

      end if

   end subroutine ccmpi_distribute3

   subroutine ccmpi_gather2(a,ag)
      ! Collect global arrays.

      real, dimension(ifull), intent(in) :: a
      real, dimension(ifull_g), intent(out), optional :: ag
      integer :: ierr, itag = 0, iproc
      integer, dimension(MPI_STATUS_SIZE) :: status
      real, dimension(ifull) :: abuf
      integer :: ipoff, jpoff, npoff
      integer :: i, j, n, iq, iqg

      if ( myid == 0 .and. .not. present(ag) ) then
         print*, "Error: ccmpi_gather argument required on proc 0"
         stop
      end if

      itag = itag + 1
      if ( myid == 0 ) then
         ! Use the face indices for unpacking
         do n=1,npan
            do j=1,jpan
               do i=1,ipan
                  iqg = indg(i,j,n)  ! True global index
                  iq = indp(i,j,n)
                  ag(iqg) = a(iq)
               end do
            end do
         end do
         do iproc=1,nproc-1
            call MPI_Recv( abuf, ipan*jpan*npan, MPI_REAL, iproc, itag, &
                        MPI_COMM_WORLD, status, ierr )
            ! Panel range on the source processor
            call proc_region(iproc,ipoff,jpoff,npoff)
            ! Use the face indices for unpacking
            do n=1,npan
               do j=1,jpan
                  do i=1,ipan
                     ! Global indices are i+ipoff, j+jpoff, n-npoff
                     iqg = ind(i+ipoff,j+jpoff,n-npoff) ! True global 1D index
                     iq = indp(i,j,n)
                     ag(iqg) = abuf(iq)
                  end do
               end do
            end do
         end do
      else
         abuf = a
         call MPI_SSend( abuf, ipan*jpan*npan, MPI_REAL, 0, itag, &
                         MPI_COMM_WORLD, ierr )
      end if

   end subroutine ccmpi_gather2

   subroutine ccmpi_gather3(a,ag)
      ! Collect global arrays.

      real, dimension(ifull,kl), intent(in) :: a
      real, dimension(ifull_g,kl), intent(out), optional :: ag
      integer :: ierr, itag = 0, iproc
      integer, dimension(MPI_STATUS_SIZE) :: status
      real, dimension(ifull,kl) :: abuf
      integer :: ipoff, jpoff, npoff
      integer :: i, j, n, iq, iqg

      if ( myid == 0 .and. .not. present(ag) ) then
         print*, "Error: ccmpi_gather argument required on proc 0"
         stop
      end if

      itag = itag + 1
      if ( myid == 0 ) then
         ! Use the face indices for unpacking
         do n=1,npan
            do j=1,jpan
               do i=1,ipan
                  iqg = indg(i,j,n)  ! True global index
                  iq = indp(i,j,n)
                  ag(iqg,:) = a(iq,:)
               end do
            end do
         end do
         do iproc=1,nproc-1
            call MPI_Recv( abuf, ipan*jpan*npan*kl, MPI_REAL, iproc, itag, &
                        MPI_COMM_WORLD, status, ierr )
            ! Panel range on the source processor
            call proc_region(iproc,ipoff,jpoff,npoff)
            ! Use the face indices for unpacking
            do n=1,npan
               do j=1,jpan
                  do i=1,ipan
                     iqg = ind(i+ipoff,j+jpoff,n-npoff) ! True global 1D index
                     iq = indp(i,j,n)
                     ag(iqg,:) = abuf(iq,:)
                  end do
               end do
            end do
         end do
      else
         abuf = a
         call MPI_SSend( abuf, ipan*jpan*npan*kl, MPI_REAL, 0, itag, &
                         MPI_COMM_WORLD, ierr )
      end if

   end subroutine ccmpi_gather3

   subroutine bounds_setup()

      include 'indices.h'
      include 'indices_g.h'
      integer :: n, ir, nr, i, j, iq, iqx, count
      logical :: double
      integer :: bstart, ierr, itag = 0, iproc, rproc, sproc
      integer, dimension(MPI_STATUS_SIZE,2*nproc) :: status
      real :: tmp
      integer :: i2, j2, n2, iqg
      integer :: iext, iql, iloc, jloc, nloc

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
      ieu = huge(1)
      iwu = huge(1)
      inv = huge(1)
      isv = huge(1)
      lwws = huge(1)
      lws = huge(1)
      lwss = huge(1)
      les = huge(1)
      lees = huge(1)
      less = huge(1)
      lwwn = huge(1)
      lwnn = huge(1)
      leen = huge(1)
      lenn = huge(1)
      lsww = huge(1)
      lsw = huge(1)
      lssw = huge(1)
      lsee = huge(1)
      lsse = huge(1)
      lnww = huge(1)
      lnw = huge(1)
      lnnw = huge(1)
      lnee = huge(1)
      lnne = huge(1)
      do n=1,npan
         do j=1,jpan
            do i=1,ipan
               iq = indp(i,j,n)   ! Local
               iqg = indg(i,j,n)  ! Global

               iqx = in_g(iqg)    ! Global neighbour index
               rproc = qproc(iqx) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqx to local value
                  call indv_mpi(iqx,iloc,jloc,nloc)
                  in(iq) = indp(iloc,jloc,nloc)
               end if
               iqx = inn_g(iqg)    ! Global neighbour index
               rproc = qproc(iqx) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqx to local value
                  call indv_mpi(iqx,iloc,jloc,nloc)
                  inn(iq) = indp(iloc,jloc,nloc)
               end if

               iqx = is_g(iqg)    ! Global neighbour index
               rproc = qproc(iqx) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqx to local value
                  call indv_mpi(iqx,iloc,jloc,nloc)
                  is(iq) = indp(iloc,jloc,nloc)
               end if
               iqx = iss_g(iqg)    ! Global neighbour index
               rproc = qproc(iqx) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqx to local value
                  call indv_mpi(iqx,iloc,jloc,nloc)
                  iss(iq) = indp(iloc,jloc,nloc)
               end if

               iqx = ie_g(iqg)    ! Global neighbour index
               rproc = qproc(iqx) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqx to local value
                  call indv_mpi(iqx,iloc,jloc,nloc)
                  ie(iq) = indp(iloc,jloc,nloc)
               end if
               iqx = iee_g(iqg)    ! Global neighbour index
               rproc = qproc(iqx) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqx to local value
                  call indv_mpi(iqx,iloc,jloc,nloc)
                  iee(iq) = indp(iloc,jloc,nloc)
               end if

               iqx = iw_g(iqg)    ! Global neighbour index
               rproc = qproc(iqx) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqx to local value
                  call indv_mpi(iqx,iloc,jloc,nloc)
                  iw(iq) = indp(iloc,jloc,nloc)
               end if
               iqx = iww_g(iqg)    ! Global neighbour index
               rproc = qproc(iqx) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqx to local value
                  call indv_mpi(iqx,iloc,jloc,nloc)
                  iww(iq) = indp(iloc,jloc,nloc)
               end if

               ! Note that the model only needs a limited set of the diagonal
               ! index arrays
               iqx = ine_g(iqg)    ! Global neighbour index
               rproc = qproc(iqx) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqx to local value
                  call indv_mpi(iqx,iloc,jloc,nloc)
                  ine(iq) = indp(iloc,jloc,nloc)
               end if

               iqx = ise_g(iqg)    ! Global neighbour index
               rproc = qproc(iqx) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqx to local value
                  call indv_mpi(iqx,iloc,jloc,nloc)
                  ise(iq) = indp(iloc,jloc,nloc)
               end if

               iqx = ien_g(iqg)    ! Global neighbour index
               rproc = qproc(iqx) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqx to local value
                  call indv_mpi(iqx,iloc,jloc,nloc)
                  ien(iq) = indp(iloc,jloc,nloc)
               end if

               iqx = iwn_g(iqg)    ! Global neighbour index
               rproc = qproc(iqx) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqx to local value
                  call indv_mpi(iqx,iloc,jloc,nloc)
                  iwn(iq) = indp(iloc,jloc,nloc)
               end if

            end do
         end do
      end do

      ! Correct within the same face only ( not necessarily the same
      ! processor, but will be corrected later.
      ieu = ie
      iwu = iw
      inv = in
      isv = is

      ! Initialise the edge variables
      edge_w = ioff == 0
      edge_s = joff == 0
      edge_n = joff == il-jpan
      edge_e = ioff == il-ipan

      ! Allocate array to hold values for each processor, including self.
!      allocate(bnds(0:nproc-1))
!      allocate(ireq(2*nproc))

      nr = 1

      bnds(:)%len = 0
      bnds(:)%rlen = 0
      bnds(:)%slen = 0
      bnds(:)%rlen2 = 0
      bnds(:)%slen2 = 0


!     In the first pass through, set up list of points to be requested from
!     other processors. This is basically the same in the face and 1D versions,
!     except that the face version calculates an index for unpacking into
!     the boundary region, while the 1D code puts it at the end and
!     updates the index arrays appropriately.

      iext = 0
      do n=1,npan

         !     Start with W edge
         i = 1
         do j=1,jpan
            ! 1D code takes care of corners separately at the end so only goes
            ! over 1,jpan here.
               iq = indg(i,j,n)
               iqx = iw_g(iq)

            ! iqx is the global index of the required neighbouring point.

            ! Which processor has this point
            rproc = qproc(iqx)
            if ( rproc == myid ) cycle ! Don't add points already on this proc.
            ! Add this point to request list
            call check_bnds_alloc(rproc, iext)
            bnds(rproc)%rlen = bnds(rproc)%rlen + 1
            bnds(rproc)%request_list(bnds(rproc)%rlen) = iqx
            ! Increment extended region index
            iext = iext + 1
            bnds(rproc)%unpack_list(bnds(rproc)%rlen) = iext
            iql = indp(i,j,n)  !  Local index
            iw(iql) = ifull+iext
         end do

         !     N edge
         j=jpan
         do i=1,ipan
               iq = indg(i,j,n)
               iqx = in_g(iq)

            ! Which processor has this point
            rproc = qproc(iqx)
            if ( rproc == myid ) cycle ! Don't add points already on this proc.
            ! Add this point to request list
            call check_bnds_alloc(rproc, iext)
            bnds(rproc)%rlen = bnds(rproc)%rlen + 1
            bnds(rproc)%request_list(bnds(rproc)%rlen) = iqx
            ! Increment extended region index
            iext = iext + 1
            bnds(rproc)%unpack_list(bnds(rproc)%rlen) = iext
            iql = indp(i,j,n)  !  Local index
            in(iql) = ifull+iext
         end do

         !     E edge
         i = ipan
         do j=1,jpan
               iq = indg(i,j,n)
               iqx = ie_g(iq)

            ! Which processor has this point
            rproc = qproc(iqx)
            if ( rproc == myid ) cycle ! Don't add points already on this proc.
            ! Add this point to request list
            call check_bnds_alloc(rproc, iext)
            bnds(rproc)%rlen = bnds(rproc)%rlen + 1
            bnds(rproc)%request_list(bnds(rproc)%rlen) = iqx
            ! Increment extended region index
            iext = iext + 1
            bnds(rproc)%unpack_list(bnds(rproc)%rlen) = iext
            iql = indp(i,j,n)  !  Local index
            ie(iql) = ifull+iext
         end do

         !     S edge
         j=1
         do i=1,ipan
               iq = indg(i,j,n)
               iqx = is_g(iq)

            ! Which processor has this point
            rproc = qproc(iqx)
            if ( rproc == myid ) cycle ! Don't add points already on this proc.
            ! Add this point to request list
            call check_bnds_alloc(rproc, iext)
            bnds(rproc)%rlen = bnds(rproc)%rlen + 1
            bnds(rproc)%request_list(bnds(rproc)%rlen) = iqx
            ! Increment extended region index
            iext = iext + 1
            bnds(rproc)%unpack_list(bnds(rproc)%rlen) = iext
            iql = indp(i,j,n)  !  Local index
            is(iql) = ifull+iext
         end do
      end do ! n=1,npan

!     Now handle the special corner values that need to be remapped
!     This adds to rlen, so needs to come before the _XX stuff.
      do n=1,npan
         ! NE, EN
         iq = indp(ipan,jpan,n)
         iqg = indg(ipan,jpan,n)
!        This might matter for leen etc but not here
!         if ( edge_e .and. edge_n ) then
!            ! This is a real vertex
         iqx = ine_g(iqg)
         ! Which processor has this point
         rproc = qproc(iqx)
         if ( rproc /= myid ) then ! Add to list
            call check_bnds_alloc(rproc, iext)
            bnds(rproc)%rlen = bnds(rproc)%rlen + 1
            bnds(rproc)%request_list(bnds(rproc)%rlen) = iqx
            ! Increment extended region index
            iext = iext + 1
            bnds(rproc)%unpack_list(bnds(rproc)%rlen) = iext
            ine(iq) = ifull+iext
         end if

         iqx = ien_g(iqg)
         ! Which processor has this point
         rproc = qproc(iqx)
         if ( rproc /= myid ) then ! Add to list
            call check_bnds_alloc(rproc, iext)
            bnds(rproc)%rlen = bnds(rproc)%rlen + 1
            bnds(rproc)%request_list(bnds(rproc)%rlen) = iqx
            ! Increment extended region index
            iext = iext + 1
            bnds(rproc)%unpack_list(bnds(rproc)%rlen) = iext
            ien(iq) = ifull+iext
         end if

         iq = indp(ipan,1,n)
         iqg = indg(ipan,1,n)
         iqx = ise_g(iqg)
         ! Which processor has this point
         rproc = qproc(iqx)
         if ( rproc /= myid ) then ! Add to list
            call check_bnds_alloc(rproc, iext)
            bnds(rproc)%rlen = bnds(rproc)%rlen + 1
            bnds(rproc)%request_list(bnds(rproc)%rlen) = iqx
            ! Increment extended region index
            iext = iext + 1
            bnds(rproc)%unpack_list(bnds(rproc)%rlen) = iext
            ise(iq) = ifull+iext
         end if

         iq = indp(1,jpan,n)
         iqg = indg(1,jpan,n)
         iqx = iwn_g(iqg)
         ! Which processor has this point
         rproc = qproc(iqx)
         if ( rproc /= myid ) then ! Add to list
            call check_bnds_alloc(rproc, iext)
            bnds(rproc)%rlen = bnds(rproc)%rlen + 1
            bnds(rproc)%request_list(bnds(rproc)%rlen) = iqx
            ! Increment extended region index
            iext = iext + 1
            bnds(rproc)%unpack_list(bnds(rproc)%rlen) = iext
            iwn(iq) = ifull+iext
         end if

!        Special vertex arrays.
!        First order ws, es, sw, nw

         iq = indp(1,1,n)
         iqg = indg(1,1,n)
         if ( edge_w .and. edge_s ) then
            iqx = lws_g(n-noff)  ! This is the true face index
         else if ( edge_w ) then
            ! Middle of W edge, use is first because this is on the same
            ! face
            iqx = iw_g(is_g(iqg))
         else
            iqx = is_g(iw_g(iqg))
         end if
         call fix_index(iqx,lws,n,bnds,iext)

         if ( edge_w .and. edge_s ) then
            iqx = lsw_g(n-noff)  ! This is the true face index
         else if ( edge_w ) then
            ! Middle of W edge, use is first because this is on the same
            ! face
            iqx = iw_g(is_g(iqg))
         else
            iqx = is_g(iw_g(iqg))
         end if
         call fix_index(iqx,lsw,n,bnds,iext)

         iq = indp(ipan,1,n)
         iqg = indg(ipan,1,n)
         if ( edge_e .and. edge_s ) then
            iqx = les_g(n-noff)  ! This is the true face index
         else if ( edge_e ) then
            ! Middle of E edge, use is first because this is on the same
            ! face
            iqx = ie_g(is_g(iqg))
         else
            iqx = is_g(ie_g(iqg))
         end if
         call fix_index(iqx,les,n,bnds,iext)

         iq = indp(1,jpan,n)
         iqg = indg(1,jpan,n)
         if ( edge_w .and. edge_n ) then
            iqx = lnw_g(n-noff)  ! This is the true face index
         else if ( edge_w ) then
            ! Middle of W edge, use in first because this is on the same
            ! face
            iqx = iw_g(in_g(iqg))
         else
            iqx = in_g(iw_g(iqg))
         end if
         call fix_index(iqx,lnw,n,bnds,iext)

      end do



!     Now set up the second row
      bnds(:)%rlen2 =  bnds(:)%rlen  ! so that they're appended.

      do n=1,npan

         !     Start with W edge
         i = 1
         do j=1,jpan
               iq = indg(i,j,n)
               iqx = iww_g(iq)

            ! Which processor has this point
            rproc = qproc(iqx)
            if ( rproc == myid ) cycle ! Don't add points already on this proc.
            ! Add this point to request list
            call check_bnds_alloc(rproc, iext)
            bnds(rproc)%rlen2 = bnds(rproc)%rlen2 + 1
            bnds(rproc)%request_list(bnds(rproc)%rlen2) = iqx
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
               iqx = inn_g(iq)

            ! Which processor has this point
            rproc = qproc(iqx)
            if ( rproc == myid ) cycle ! Don't add points already on this proc.
            ! Add this point to request list
            call check_bnds_alloc(rproc, iext)
            bnds(rproc)%rlen2 = bnds(rproc)%rlen2 + 1
            bnds(rproc)%request_list(bnds(rproc)%rlen2) = iqx
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
               iqx = iee_g(iq)

            ! Which processor has this point
            rproc = qproc(iqx)
            if ( rproc == myid ) cycle ! Don't add points already on this proc.
            ! Add this point to request list
            call check_bnds_alloc(rproc, iext)
            bnds(rproc)%rlen2 = bnds(rproc)%rlen2 + 1
            bnds(rproc)%request_list(bnds(rproc)%rlen2) = iqx
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
               iqx = iss_g(iq)

            ! Which processor has this point
            rproc = qproc(iqx)
            if ( rproc == myid ) cycle ! Don't add points already on this proc.
            ! Add this point to request list
            call check_bnds_alloc(rproc, iext)
            bnds(rproc)%rlen2 = bnds(rproc)%rlen2 + 1
            bnds(rproc)%request_list(bnds(rproc)%rlen2) = iqx
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
      do n=1,npan
         do j=1,jpan
            iww(indp(2,j,n)) = iw(indp(1,j,n))
            iee(indp(ipan-1,j,n)) = ie(indp(ipan,j,n))
         end do
         do i=1,ipan
            iss(indp(i,2,n)) = is(indp(i,1,n))
            inn(indp(i,jpan-1,n)) = in(indp(i,jpan,n))
         end do
      end do

!     Set up the diagonal index arrays. Most of the points here will have
!     already been added to copy lists above. The corners are handled
!     separately here. This means some points may be copied twice but it's
!     a small overhead.
      do n=1,npan
         do j=1,jpan
            do i=1,ipan
               iq = indp(i,j,n)   ! Local
               ! Except at corners, ien = ine etc.
               if ( i < ipan ) then
                  ! ie will be defined
                  ine(iq) = in(ie(iq))
                  ise(iq) = is(ie(iq))
               else
                  ! i = ipan, ie will have been remapped
                  if ( j > 1 )    ise(iq) = ie(is(iq))
                  if ( j < jpan ) ine(iq) = ie(in(iq))
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
         print*, "IEXT too large", iext, iextra
         call MPI_Abort(MPI_COMM_WORLD)
      end if



!     Now, for each processor send the list of points I want.
!     The state of being a neighbour is reflexive so only expect to
!     recv from those processors I send to (are there grid arrangements for
!     which this would not be true?)

      nreq = 0
      do iproc = 1,nproc-1  !
         ! Send and recv from same proc
         sproc = modulo(myid+iproc,nproc)  ! Send to
         if (bnds(sproc)%rlen > 0 ) then
            ! Send list of requests
            nreq = nreq + 1
            ! Using array(1) rather than array is neccessary on the NEC
            ! (only a problem with pointers, not regular arrays).
            call MPI_ISend( bnds(sproc)%request_list(1), bnds(sproc)%rlen, &
                 MPI_INTEGER, sproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
            nreq = nreq + 1
            call MPI_IRecv( bnds(sproc)%send_list(1), bnds(sproc)%len, &
                 MPI_INTEGER, sproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do

      if ( nreq > 0 ) then
         call MPI_Waitall(nreq,ireq,status,ierr)
      end if

!     Now get the actual sizes from the status
      nreq = 0
      do iproc = 1,nproc-1  !
         sproc = modulo(myid+iproc,nproc)
         if (bnds(sproc)%rlen > 0 ) then
            ! To get recv status, advance nreq by 2
            nreq = nreq + 2
            call MPI_Get_count(status(1,nreq), MPI_INTEGER, count, ierr)
            ! This the number of points I have to send to rproc.
            bnds(sproc)%slen = count
         end if
      end do

!     For simplicity do exactly the same thing for rlen2, even though the
!     data are largely repeated.
      nreq = 0
      do iproc = 1,nproc-1  !
         sproc = modulo(myid+iproc,nproc)  ! Send to
         if (bnds(sproc)%rlen > 0 ) then
            ! Send list of requests
            nreq = nreq + 1
            call MPI_ISend( bnds(sproc)%request_list(1), bnds(sproc)%rlen2, &
                 MPI_INTEGER, sproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
            nreq = nreq + 1
            ! Use the maximum size in the recv call.
            call MPI_IRecv( bnds(sproc)%send_list(1), bnds(sproc)%len, &
                 MPI_INTEGER, sproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do
      if ( nreq > 0 ) then
         call MPI_Waitall(nreq,ireq,status,ierr)
      end if

!     Now get the actual sizes from the status
      nreq = 0
      do iproc = 1,nproc-1  !
         sproc = modulo(myid+iproc,nproc)  ! Recv from
         if (bnds(sproc)%rlen > 0 ) then
            ! To get recv status, advance nreq by 2
            nreq = nreq + 2
            call MPI_Get_count(status(1,nreq), MPI_INTEGER, count, ierr)
            ! This the number of points I have to send to rproc.
            bnds(sproc)%slen2 = count
         end if
      end do


!     Start of UV section

      bnds(:)%rlen_uv = 0
      bnds(:)%slen_uv = 0

!     In the first pass through, set up list of points to be requested from
!     other processors. In the 1D code values on the same processor are
!     copied if they have to be swapped.
!     (Actually make this a second stage optimisation).

      iext = 0
      do n=1,npan

         !     Start with W edge, U values
         i = 1
         do j=1,jpan
            iqg = indg(i,j,n)
            iqx = iw_g(iqg)
            ! Which processor has this point
            rproc = qproc(iqx)
            call check_bnds_alloc(rproc, iext)
            bnds(rproc)%rlen_uv = bnds(rproc)%rlen_uv + 1
            bnds(rproc)%request_list_uv(bnds(rproc)%rlen_uv) = iqx
            ! Increment extended region index
            iext = iext + 1
            bnds(rproc)%unpack_list_uv(bnds(rproc)%rlen_uv) = iext
            iql = indp(i,j,n)  !  Local index
            iwu(iql) = ifull+iext
            ! Decide if u/v need to be swapped. My face is n-noff
            bnds(rproc)%uv_swap(bnds(rproc)%rlen_uv) = edge_w .and. swap_w(n-noff)
         end do

         !     N edge (V)
         j=jpan
         do i=1,ipan
            iqg = indg(i,j,n)
            iqx = in_g(iqg)
            rproc = qproc(iqx)
            ! Add this point to request list
            call check_bnds_alloc(rproc, iext)
            bnds(rproc)%rlen_uv = bnds(rproc)%rlen_uv + 1
            ! to show that this is v rather than u, flip sign
            bnds(rproc)%request_list_uv(bnds(rproc)%rlen_uv) = -iqx
            ! Increment extended region index
            iext = iext + 1
            bnds(rproc)%unpack_list_uv(bnds(rproc)%rlen_uv) = -iext
            iql = indp(i,j,n)  !  Local index
            inv(iql) = ifull+iext
            bnds(rproc)%uv_swap(bnds(rproc)%rlen_uv) = edge_n .and. swap_n(n-noff)
         end do

         !     E edge, U
         i = ipan
         do j=1,jpan
            iqg = indg(i,j,n)
            iqx = ie_g(iqg)
            rproc = qproc(iqx)
            ! Add this point to request list
            call check_bnds_alloc(rproc, iext)
            bnds(rproc)%rlen_uv = bnds(rproc)%rlen_uv + 1
            bnds(rproc)%request_list_uv(bnds(rproc)%rlen_uv) = iqx
            ! Increment extended region index
            iext = iext + 1
            bnds(rproc)%unpack_list_uv(bnds(rproc)%rlen_uv) = iext
            iql = indp(i,j,n)  !  Local index
            ieu(iql) = ifull+iext
            bnds(rproc)%uv_swap(bnds(rproc)%rlen_uv) = edge_e .and. swap_e(n-noff)
         end do

         !     S edge, V
         j=1
         do i=1,ipan
            iqg = indg(i,j,n)
            iqx = is_g(iqg)
            rproc = qproc(iqx)
            call check_bnds_alloc(rproc, iext)
            bnds(rproc)%rlen_uv = bnds(rproc)%rlen_uv + 1
            bnds(rproc)%request_list_uv(bnds(rproc)%rlen_uv) = -iqx
            ! Increment extended region index
            iext = iext + 1
            bnds(rproc)%unpack_list_uv(bnds(rproc)%rlen_uv) = -iext
            iql = indp(i,j,n)  !  Local index
            isv(iql) = ifull+iext
            bnds(rproc)%uv_swap(bnds(rproc)%rlen_uv) = edge_s .and. swap_s(n-noff)
         end do
      end do ! n=1,npan

!     Second pass
      bnds(:)%rlen2_uv = bnds(:)%rlen_uv
      bnds(:)%slen2_uv = bnds(:)%slen_uv

      print*, "Bounds", myid, "SLEN", bnds(:)%slen
      print*, "Bounds", myid, "RLEN", bnds(:)%rlen
      print*, "Bounds", myid, "SLEN2", bnds(:)%slen2
      print*, "Bounds", myid, "RLEN2", bnds(:)%rlen2

! 1D code doesn't use double row U/V


!     Now, for each processor send the list of points I want.
!     Send all possible pairs and don't assume anything about the symmetry.
!     For completeness, send zero length messages too.
!     Get the length from the message status
!     Also have to send the swap list

      nreq = 0
      do iproc = 1,nproc-1  !
         sproc = modulo(myid+iproc,nproc)
         if ( bnds(sproc)%rlen_uv > 0 ) then
            ! Send list of requests
            nreq = nreq + 1
            call MPI_ISend( bnds(sproc)%request_list_uv(1), bnds(sproc)%rlen_uv, &
                 MPI_INTEGER, sproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
            nreq = nreq + 1
            ! Use the maximum size in the recv call.
            call MPI_IRecv(  bnds(sproc)%send_list_uv(1),  bnds(sproc)%len, &
                 MPI_INTEGER, sproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do
      if ( nreq > 0 ) then
         call MPI_Waitall(nreq,ireq,status,ierr)
      end if

!     Now get the actual sizes from the status
      nreq = 0
      do iproc = 1,nproc-1  !
         sproc = modulo(myid+iproc,nproc)
         if ( bnds(sproc)%rlen_uv > 0 ) then
            ! To get recv status, advance nreq by 2
            nreq = nreq + 2
            call MPI_Get_count(status(1,nreq), MPI_INTEGER, count, ierr)
            ! This the number of points I have to send to rproc.
            bnds(sproc)%slen_uv = count
         end if
      end do

      ! Now send the second set
      nreq = 0
      do iproc = 1,nproc-1  !
         sproc = modulo(myid+iproc,nproc)
         if ( bnds(sproc)%rlen_uv > 0 ) then
            ! Send list of requests
            nreq = nreq + 1
            call MPI_ISend( bnds(sproc)%request_list_uv(1), bnds(sproc)%rlen2_uv, &
                 MPI_INTEGER, sproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
            nreq = nreq + 1
            ! Use the maximum size in the recv call.
            call MPI_IRecv( bnds(sproc)%send_list_uv(1), bnds(sproc)%len, &
                 MPI_INTEGER, sproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do
      if ( nreq > 0 ) then
         call MPI_Waitall(nreq,ireq,status,ierr)
      end if

!     Now get the actual sizes from the status
      nreq = 0
      do iproc = 1,nproc-1  !
         sproc = modulo(myid+iproc,nproc)
         if ( bnds(sproc)%rlen_uv > 0 ) then
            ! To get recv status, advance nreq by 2
            nreq = nreq + 2
            call MPI_Get_count(status(1,nreq), MPI_INTEGER, count, ierr)
            ! This the number of points I have to send to rproc.
            bnds(sproc)%slen2_uv = count
         end if
      end do

      ! Only send the swap list once
      nreq = 0
      do iproc = 1,nproc-1  !
         sproc = modulo(myid+iproc,nproc)
         if ( bnds(sproc)%rlen_uv > 0 ) then
            ! Send list of requests
            nreq = nreq + 1
            call MPI_ISend( bnds(sproc)%uv_swap(1), bnds(sproc)%rlen2_uv, &
                 MPI_LOGICAL, sproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
            nreq = nreq + 1
            ! Use the maximum size in the recv call.
            call MPI_IRecv( bnds(sproc)%send_swap(1), bnds(sproc)%len, &
                 MPI_LOGICAL, sproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do
      if ( nreq > 0 ) then
         call MPI_Waitall(nreq,ireq,status,ierr)
      end if


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
               call check_set( ieu(iq), "IEU", i, j, n, iq)
               call check_set( iwu(iq), "IWU", i, j, n, iq)
               call check_set( inv(iq), "INV", i, j, n, iq)
               call check_set( isv(iq), "ISV", i, j, n, iq)
            end do
         end do
         call check_set( lwws(n), "LWWS", 1, 1, n, 1)
         call check_set( lws(n),  "LWS",  1, 1, n, 1)
         call check_set( lwss(n), "LWSS", 1, 1, n, 1)
         call check_set( les(n),  "LES",  1, 1, n, 1)
         call check_set( lees(n), "LEES", 1, 1, n, 1)
         call check_set( less(n), "LESS", 1, 1, n, 1)
         call check_set( lwwn(n), "LWWN", 1, 1, n, 1)
         call check_set( lwnn(n), "LWNN", 1, 1, n, 1)
         call check_set( leen(n), "LEEN", 1, 1, n, 1)
         call check_set( lenn(n), "LENN", 1, 1, n, 1)
         call check_set( lsww(n), "LSWW", 1, 1, n, 1)
         call check_set( lsw(n),  "LSW",  1, 1, n, 1)
         call check_set( lssw(n), "LSSW", 1, 1, n, 1)
         call check_set( lsee(n), "LSEE", 1, 1, n, 1)
         call check_set( lsse(n), "LSSE", 1, 1, n, 1)
         call check_set( lnww(n), "LNWW", 1, 1, n, 1)
         call check_set( lnw(n),  "LNW",  1, 1, n, 1)
         call check_set( lnnw(n), "LNNW", 1, 1, n, 1)
         call check_set( lnee(n), "LNEE", 1, 1, n, 1)
         call check_set( lnne(n), "LNNE", 1, 1, n, 1)
      end do


   end subroutine bounds_setup

   subroutine check_set(ind,str,i,j,n,iq)
      integer, intent(in) :: ind,i,j,n,iq
      character(len=*) :: str
      if ( ind == huge(1) ) then
         print*, str, " not set", myid, i, j, n, iq
         call MPI_Abort(MPI_COMM_WORLD)
      end if
   end subroutine check_set

   subroutine bounds2(t, nrows)
      ! Copy the boundary regions
      real, dimension(ifull+iextra), intent(inout) :: t
      integer, intent(in), optional :: nrows
      integer :: n, ir, i, j, iq, iqx, count
      logical :: double
      integer :: bstart, ierr, itag = 0, iproc, rproc, sproc
      integer, dimension(MPI_STATUS_SIZE,2*nproc) :: status
      real :: tmp
      integer :: i2, j2, n2
      integer :: send_len, recv_len

#ifdef mpilog
      ierr = MPE_log_event(bounds_begin,0,"")
#endif
#ifdef vampir
      call vtbegin(bounds_begin, ierr)
#endif

      double = .false.
      if (present(nrows)) then
         if ( nrows == 2 ) then
            double = .true.
         end if
      end if

!     Set up the buffers to send
      nreq = 0
      do iproc = 1,nproc-1  !
         sproc = modulo(myid+iproc,nproc)  ! Send to
         rproc = modulo(myid-iproc,nproc)  ! Recv from
         if ( double ) then
            recv_len = bnds(rproc)%rlen2
            send_len = bnds(sproc)%slen2
         else
            recv_len = bnds(rproc)%rlen
            send_len = bnds(sproc)%slen
         end if
         if ( recv_len /= 0 ) then
            nreq = nreq + 1
            call MPI_IRecv( bnds(rproc)%rbuf(1),  bnds(rproc)%len, &
                 MPI_REAL, rproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
         if ( send_len > 0 ) then
            ! Build up list of points
            do iq=1,send_len
               ! send_list(iq) is global point index, i, j, n are local
               call indv_mpi(bnds(sproc)%send_list(iq),i,j,n)
               bnds(sproc)%sbuf(iq) = t(indp(i,j,n))
            end do
            nreq = nreq + 1
            call MPI_ISend( bnds(sproc)%sbuf(1), send_len, &
                 MPI_REAL, sproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do

      if ( nreq > 0 ) then
         call MPI_Waitall(nreq,ireq,status,ierr)
      end if

      do iproc = 1,nproc-1  !
         rproc = modulo(myid-iproc,nproc)  ! Recv from
         if ( double ) then
            recv_len = bnds(rproc)%rlen2
         else
            recv_len = bnds(rproc)%rlen
         end if
         if ( recv_len > 0 ) then
            do iq=1,recv_len
               ! unpack_list(iq) is index into extended region
               t(ifull+bnds(rproc)%unpack_list(iq)) = bnds(rproc)%rbuf(iq)
            end do
         end if
      end do

      ! Finally see if there are any points on my own processor that need
      ! to be fixed up. This will only be in the case when nproc < npanels.
      if ( double ) then
         recv_len = bnds(myid)%rlen2
      else
         recv_len = bnds(myid)%rlen
      end if
      if ( recv_len > 0 ) then
         do iq=1,recv_len
            ! request_list is same as send_list in this case
            call indv_mpi(bnds(myid)%request_list(iq),i,j,n)
            ! i, j, n are local
            t(ifull+bnds(myid)%unpack_list(iq)) = t(indp(i,j,n))
         end do
      end if

#ifdef mpilog
      ierr = MPE_log_event(bounds_end,0,"")
#endif
#ifdef vampir
      call vtend(bounds_end, ierr)
#endif

   end subroutine bounds2

   subroutine bounds3(t, nrows, klim)
      ! Copy the boundary regions. Only this routine requires the extra klim
      ! argument (for helmsol).
      real, dimension(ifull+iextra,kl), intent(inout) :: t
      integer, intent(in), optional :: nrows, klim
      integer :: n, ir, i, j, iq, iqx, count
      logical :: double
      integer :: bstart, ierr, itag = 0, iproc, rproc, sproc
      integer, dimension(MPI_STATUS_SIZE,2*nproc) :: status
      real :: tmp
      integer :: i2, j2, n2
      integer :: send_len, recv_len, kx


#ifdef mpilog
      ierr = MPE_log_event(bounds_begin,0,"")
#endif
#endif
#ifdef vampir
      call vtbegin(bounds_begin, ierr)
#endif

      double = .false.
      if (present(nrows)) then
         if ( nrows == 2 ) then
            double = .true.
         end if
      end if
      if ( present(klim) ) then
         kx = klim
      else
         kx = kl
      end if

!     Set up the buffers to send
      nreq = 0
      do iproc = 1,nproc-1  !
         sproc = modulo(myid+iproc,nproc)  ! Send to
         rproc = modulo(myid-iproc,nproc)  ! Recv from
         if ( double ) then
            recv_len = bnds(rproc)%rlen2
            send_len = bnds(sproc)%slen2
         else
            recv_len = bnds(rproc)%rlen
            send_len = bnds(sproc)%slen
         end if
         if ( recv_len /= 0 ) then
            nreq = nreq + 1
            call MPI_IRecv( bnds(rproc)%rbuf(1), bnds(rproc)%len, &
                 MPI_REAL, rproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
         if ( send_len > 0 ) then
            ! Build up list of points
            do iq=1,send_len
               ! send_list(iq) is point index.
               call indv_mpi(bnds(sproc)%send_list(iq),i,j,n)
               bnds(sproc)%sbuf(1+(iq-1)*kx:iq*kx) = t(indp(i,j,n),1:kx)
            end do
            nreq = nreq + 1
            call MPI_ISend( bnds(sproc)%sbuf(1), send_len*kx, &
                 MPI_REAL, sproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do

      if ( nreq > 0 ) then
         call MPI_Waitall(nreq,ireq,status,ierr)
      end if

      do iproc = 1,nproc-1  !
         rproc = modulo(myid-iproc,nproc)  ! Recv from
         if ( double ) then
            recv_len = bnds(rproc)%rlen2
         else
            recv_len = bnds(rproc)%rlen
         end if
         if ( recv_len > 0 ) then
            do iq=1,recv_len
               ! i, j, n are local
               t(ifull+bnds(rproc)%unpack_list(iq),1:kx) = bnds(rproc)%rbuf(1+(iq-1)*kx:iq*kx)
            end do
         end if
      end do

      ! Finally see if there are any points on my own processor that need
      ! to be fixed up. This will only be in the case when nproc < npanels.
      if ( double ) then
         recv_len = bnds(myid)%rlen2
      else
         recv_len = bnds(myid)%rlen
      end if
      if ( recv_len > 0 ) then
         do iq=1,recv_len
            ! request_list is same as send_list in this case
            call indv_mpi(bnds(myid)%request_list(iq),i,j,n)
            ! i, j, n are local
            t(ifull+bnds(myid)%unpack_list(iq),:) = t(indp(i,j,n),:)
         end do
      end if

#ifdef mpilog
      ierr = MPE_log_event(bounds_end,0,"")
#endif
#ifdef vampir
      call vtend(bounds_end, ierr)
#endif

   end subroutine bounds3

   subroutine boundsuv2(u, v, nrows)
      ! Copy the boundary regions of u and v. This doesn't require the
      ! diagonal points like (0,0), but does have to take care of the
      ! direction changes.
      real, dimension(ifull+iextra), intent(inout) :: u, v
      integer, intent(in), optional :: nrows
      integer :: n, ir, i, j, iq, iqx, count
      logical :: double
      integer :: bstart, ierr, itag = 0, iproc, rproc, sproc
      integer, dimension(MPI_STATUS_SIZE,2*nproc) :: status
      real :: tmp
      integer :: i2, j2, n2
      integer :: send_len, recv_len

#ifdef mpilog
      ierr = MPE_log_event(boundsuv_begin,0,"")
#endif
#ifdef vampir
      call vtbegin(boundsuv_begin, ierr)
#endif

      double = .false.
      if (present(nrows)) then
         if ( nrows == 2 ) then
            double = .true.
         end if
      end if

      if ( double ) then
         print*, "NROWS=2 not implemented in 1D boundsuv"
         call MPI_Abort(MPI_COMM_WORLD)
      end if

!     Set up the buffers to send
      nreq = 0
      do iproc = 1,nproc-1  !
         sproc = modulo(myid+iproc,nproc)  ! Send to
         rproc = modulo(myid-iproc,nproc)  ! Recv from
         if ( double ) then
            recv_len = bnds(rproc)%rlen2_uv
            send_len = bnds(sproc)%slen2_uv
         else
            recv_len = bnds(rproc)%rlen_uv
            send_len = bnds(sproc)%slen_uv
         end if
         if ( recv_len /= 0 ) then
            nreq = nreq + 1
            call MPI_IRecv( bnds(rproc)%rbuf(1), bnds(rproc)%len, &
                 MPI_REAL, rproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
         if ( send_len > 0 ) then
            ! Build up list of points
            do iq=1,send_len
               ! send_list_uv(iq) is point index.
               ! Use abs because sign is used as u/v flag
               call indv_mpi(abs(bnds(sproc)%send_list_uv(iq)),i,j,n)
               !
               if ( (bnds(sproc)%send_list_uv(iq) > 0) .neqv. &
                     bnds(sproc)%send_swap(iq) ) then
                  bnds(sproc)%sbuf(iq) = u(indp(i,j,n))
               else
                  bnds(sproc)%sbuf(iq) = v(indp(i,j,n))
               end if
            end do
            nreq = nreq + 1
            call MPI_ISend( bnds(sproc)%sbuf(1), send_len, &
                 MPI_REAL, sproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do

      if ( nreq > 0 ) then
         call MPI_Waitall(nreq,ireq,status,ierr)
      end if

      do iproc = 1,nproc-1  !
         rproc = modulo(myid-iproc,nproc)  ! Recv from
         if ( double ) then
            recv_len = bnds(rproc)%rlen2_uv
         else
            recv_len = bnds(rproc)%rlen_uv
         end if
         if ( recv_len > 0 ) then
            do iq=1,recv_len
               ! unpack_list(iq) is index into extended region
               if ( bnds(rproc)%unpack_list_uv(iq) > 0 ) then
                  u(ifull+bnds(rproc)%unpack_list_uv(iq)) = bnds(rproc)%rbuf(iq)
               else
                  v(ifull-bnds(rproc)%unpack_list_uv(iq)) = bnds(rproc)%rbuf(iq)
               end if
            end do
         end if
      end do

      ! Finally see if there are any points on my own processor that need
      ! to be fixed up. This will only be in the case when nproc < npanels.
      if ( bnds(myid)%rlen_uv > 0 ) then
         if ( double ) then
            recv_len = bnds(myid)%rlen2_uv
         else
            recv_len = bnds(myid)%rlen_uv
         end if
         do iq=1,recv_len
            ! request_list is same as send_list in this case
            call indv_mpi(abs( bnds(myid)%request_list_uv(iq)),i,j,n)
            if ( ( bnds(myid)%request_list_uv(iq) > 0) .neqv. &
                      bnds(myid)%uv_swap(iq) ) then  ! haven't copied to send_swap yet
               tmp = u(indp(i,j,n))
            else
               tmp = v(indp(i,j,n))
            end if
            ! unpack_list(iq) is index into extended region
            if ( bnds(myid)%unpack_list_uv(iq) > 0 ) then
               u(ifull+bnds(myid)%unpack_list_uv(iq)) = tmp
            else
               v(ifull-bnds(myid)%unpack_list_uv(iq)) = tmp
            end if
         end do
      end if

#ifdef mpilog
      ierr = MPE_log_event(boundsuv_end,0,"")
#endif
#ifdef vampir
      call vtend(boundsuv_end, ierr)
#endif

   end subroutine boundsuv2

   subroutine boundsuv3(u, v, nrows)
      ! Copy the boundary regions of u and v. This doesn't require the
      ! diagonal points like (0,0), but does have to take care of the
      ! direction changes.
      real, dimension(ifull+iextra,kl), intent(inout) :: u, v
      integer, intent(in), optional :: nrows
      integer :: n, ir, i, j, iq, iqx, count
      logical :: double
      integer :: bstart, ierr, itag = 0, iproc, rproc, sproc
      integer, dimension(0:nproc) :: buflen
      integer, dimension(MPI_STATUS_SIZE,2*nproc) :: status
      real, dimension(kl) :: tmp
      integer :: i2, j2, n2
      integer :: send_len, recv_len

#ifdef mpilog
      ierr = MPE_log_event(boundsuv_begin,0,"")
#endif
#ifdef vampir
      call vtbegin(boundsuv_begin, ierr)
#endif

      double = .false.
      if (present(nrows)) then
         if ( nrows == 2 ) then
            double = .true.
         end if
      end if

      if ( double ) then
         print*, "NROWS=2 not implemented in 1D boundsuv"
         call MPI_Abort(MPI_COMM_WORLD)
      end if

!     Set up the buffers to send
      nreq = 0
      do iproc = 1,nproc-1  !
         sproc = modulo(myid+iproc,nproc)  ! Send to
         rproc = modulo(myid-iproc,nproc)  ! Recv from
         if ( double ) then
            recv_len = bnds(rproc)%rlen2_uv
            send_len = bnds(sproc)%slen2_uv
         else
            recv_len = bnds(rproc)%rlen_uv
            send_len = bnds(sproc)%slen_uv
         end if
         if ( recv_len /= 0 ) then
            nreq = nreq + 1
            call MPI_IRecv( bnds(rproc)%rbuf(1), bnds(rproc)%len, &
                 MPI_REAL, rproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
         if ( send_len > 0 ) then
            ! Build up list of points
            do iq=1,send_len
               ! send_list_uv(iq) is point index.
               ! Use abs because sign is used as u/v flag
               call indv_mpi(abs(bnds(sproc)%send_list_uv(iq)),i,j,n)
               !
               if ( (bnds(sproc)%send_list_uv(iq) > 0) .neqv. &
                     bnds(sproc)%send_swap(iq) ) then
                  bnds(sproc)%sbuf(1+(iq-1)*kl:iq*kl) = u(indp(i,j,n),:)
               else
                  bnds(sproc)%sbuf(1+(iq-1)*kl:iq*kl) = v(indp(i,j,n),:)
               end if
            end do
            nreq = nreq + 1
            call MPI_ISend( bnds(sproc)%sbuf(1), send_len*kl, &
                 MPI_REAL, sproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do

      if ( nreq > 0 ) then
         call MPI_Waitall(nreq,ireq,status,ierr)
      end if

      do iproc = 1,nproc-1  !
         rproc = modulo(myid-iproc,nproc)  ! Recv from
         if ( double ) then
            recv_len = bnds(rproc)%rlen2_uv
         else
            recv_len = bnds(rproc)%rlen_uv
         end if
         if ( recv_len > 0 ) then
            do iq=1,recv_len
               ! unpack_list(iq) is index into extended region
               if ( bnds(rproc)%unpack_list_uv(iq) > 0 ) then
                  u(ifull+bnds(rproc)%unpack_list_uv(iq),:) = bnds(rproc)%rbuf(1+(iq-1)*kl:iq*kl)
               else
                  v(ifull-bnds(rproc)%unpack_list_uv(iq),:) = bnds(rproc)%rbuf(1+(iq-1)*kl:iq*kl)
               end if
            end do
         end if
      end do

      ! Finally see if there are any points on my own processor that need
      ! to be fixed up. This will only be in the case when nproc < npanels.
      if ( bnds(myid)%rlen_uv > 0 ) then
         if ( double ) then
            recv_len = bnds(myid)%rlen2_uv
         else
            recv_len = bnds(myid)%rlen_uv
         end if
         do iq=1,recv_len
            ! request_list is same as send_list in this case
            call indv_mpi(abs(bnds(myid)%request_list_uv(iq)),i,j,n)
            if ( (bnds(myid)%request_list_uv(iq) > 0) .neqv. &
                     bnds(myid)%uv_swap(iq) ) then  ! haven't copied to send_swap yet
               tmp = u(indp(i,j,n),:)
            else
               tmp = v(indp(i,j,n),:)
            end if
            ! unpack_list(iq) is index into extended region
            if ( bnds(myid)%unpack_list_uv(iq) > 0 ) then
               u(ifull+bnds(myid)%unpack_list_uv(iq),:) = tmp
            else
               v(ifull-bnds(myid)%unpack_list_uv(iq),:) = tmp
            end if
         end do
      end if

#ifdef mpilog
      ierr = MPE_log_event(boundsuv_end,0,"")
#endif
#ifdef vampir
      call vtend(boundsuv_end, ierr)
#endif

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
      integer, dimension(ifull,kl), intent(in) :: nface
      real, dimension(ifull,kl), intent(in) :: xg, yg
      integer :: iproc
      real, dimension(4,maxsize,0:nproc) :: buf
      integer :: nreq, itag = 99, ierr, rproc, sproc
      integer, dimension(2*nproc) :: ireq
      integer, dimension(MPI_STATUS_SIZE,2*nproc) :: status
      real, dimension(0:nproc) :: dslen_r, drlen_r
      integer :: count, ip, jp
      integer :: iq, iqk, k, idel, jdel, nf

#ifdef mpilog
      ierr = MPE_log_event(deptsync_begin,0,"")
#endif
#ifdef vampir
      call vtbegin(deptsync_begin, ierr)
#endif
      dslen = 0
      drlen = 0
      dindex = 0
      do k=1,kl
         do iq=1,ifull
            nf = nface(iq,k) + noff ! Make this a local index
            idel = int(xg(iq,k)) - ioff
            jdel = int(yg(iq,k)) - joff
            if ( idel < 0 .or. idel > ipan .or. jdel < 0 .or. jdel > jpan &
                 .or. nf < 1 .or. nf > npan ) then
               ! If point is on a different face, add to a list 
               ip = min(il_g,max(1,nint(xg(iq,k))))
               jp = min(il_g,max(1,nint(yg(iq,k))))
               iproc = fproc(ip,jp,nface(iq,k))
               if ( iproc == myid ) then
                  print*, "Inconsistency in deptsync"
                  stop
               end if
!!$            if ( mydiag .and. iq==1 ) then
!!$               print*, "DEPTS", k, xg(iq,k), xg(iq,k), nface(iq,k), ip,jp, iproc
!!$            end if
               ! Add this point to the list of requests I need to send to iproc
               dslen(iproc) = dslen(iproc) + 1
               call checksize(dslen(iproc),"Deptssync")
               ! Since nface is a small integer it can be exactly represented by a
               ! real. It's simpler to send like this than use a proper structure.
               buf(:,dslen(iproc),iproc) = (/ real(nface(iq,k)), xg(iq,k), yg(iq,k), real(k) /)
               dindex(:,dslen(iproc),iproc) = (/ iq,k /)
            end if
         end do
      end do
!     In this case the length of each buffer is unknown and will not
!     be symmetric between processors. Therefore need to get the length
!     from the message status


      nreq = 0
      do iproc = 1,nproc-1  !
         sproc = modulo(myid+iproc,nproc)  ! Send to
         rproc = modulo(myid-iproc,nproc)  ! Recv from
         ! Send, even if length is zero
         nreq = nreq + 1
         call MPI_ISend( buf(1,1,sproc), 4*dslen(sproc), &
                 MPI_REAL, sproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         nreq = nreq + 1
         ! Use the maximum size in the recv call.
         call MPI_IRecv( dpoints(1,1,rproc), 4*maxsize, &
                         MPI_REAL, rproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
      end do
      if ( nreq > 0 ) then
         call MPI_Waitall(nreq,ireq,status,ierr)
      end if

!     Now get the actual sizes from the status
      nreq = 0
      do iproc = 1,nproc-1  !
         rproc = modulo(myid-iproc,nproc)  ! Recv from
         ! To get recv status, advance nreq by 2
         nreq = nreq + 2
         call MPI_Get_count(status(1,nreq), MPI_REAL, count, ierr)
         drlen(rproc) = count/4
      end do

#ifdef mpilog
      ierr = MPE_log_event(deptsync_end,0,"")
#endif
#ifdef vampir
      call vtend(deptsync_end, ierr)
#endif
!!!      print*, "DEPTSYNC", myid, drlen, dindex(:,:10,:)

   end subroutine deptsync

   subroutine intssync(s)
      real, dimension(:,:), intent(inout) :: s
      integer :: iproc
      real, dimension(maxsize,0:nproc) :: buf
      integer :: nreq, itag = 0, ierr, rproc, sproc
      integer :: iq
      integer, dimension(0:nproc) :: kount
      integer, dimension(2*nproc) :: ireq
      integer, dimension(MPI_STATUS_SIZE,2*nproc) :: status

#ifdef mpilog
      ierr = MPE_log_event(intssync_begin,0,"")
#endif
#ifdef vampir
      call vtbegin(intssync_begin, ierr)
#endif

      ! When sending the results, roles of dslen and drlen are reversed
      nreq = 0
      do iproc = 1,nproc-1  !
         sproc = modulo(myid+iproc,nproc)  ! Send to
         rproc = modulo(myid-iproc,nproc)  ! Recv from
         if ( drlen(sproc) /= 0 ) then
            nreq = nreq + 1
            call MPI_ISend( sextra(1,sproc), drlen(sproc), &
                 MPI_REAL, sproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
         if ( dslen(rproc) /= 0 ) then
            nreq = nreq + 1
            call MPI_IRecv( buf(1,rproc), dslen(rproc), &
                            MPI_REAL, rproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do
      if ( nreq > 0 ) then
         call MPI_Waitall(nreq,ireq,status,ierr)
      end if

      do iproc=0,nproc-1
         if ( iproc == myid ) then
            cycle
         end if
!!!         print*, "INTSSYNC recvd", myid, buf(1:dslen(iproc),iproc)
         do iq=1,dslen(iproc)
!!!            print*, "UNPACKING", myid, dindex(1,iq,iproc),dindex(2,iq,iproc), buf(iq,iproc) 
            s(dindex(1,iq,iproc),dindex(2,iq,iproc)) = buf(iq,iproc)
         end do
      end do
#ifdef mpilog
      ierr = MPE_log_event(intssync_end,0,"")
#endif
#ifdef vampir
      call vtend(intssync_end, ierr)
#endif

   end subroutine intssync

   subroutine indv_mpi(iq, i, j, n)
      integer , intent(in) :: iq
      integer , intent(out) :: i
      integer , intent(out) :: j
      integer , intent(out) :: n

      ! Global i, j, n
      n = (iq - 1)/(il_g*il_g)
      j = 1 + (iq - n*il_g*il_g - 1)/il_g
      i = iq - (j - 1)*il_g - n*il_g*il_g
      if ( fproc(i,j,n) /= myid ) then
         write(*,"(a,5i5)") "Consistency failure in indv_mpi", myid, iq, i, j, n
         call MPI_Abort(MPI_COMM_WORLD)
      end if
      ! Reduced to values on my processor
      n = n + noff
      j = j - joff
      i = i - ioff
   end subroutine indv_mpi

   function ind(i,j,n) result(iq)
      integer, intent(in) :: i, j, n
      integer :: iq

      ! Calculate a global index from the processors indices
      ! n in range 1..npan
      iq = i + (j-1)*il_g + n*il_g*il_g
   end function ind

   function indg(i,j,n) result(iq)
      integer, intent(in) :: i, j, n
      integer :: iq

      ! Calculate a global index from the processors indices
      ! n in range 1..npan
      iq = i+ioff + (j+joff-1)*il_g + (n-noff)*il_g*il_g
   end function indg

   function indp(i,j,n) result(iq)
      integer, intent(in) :: i, j, n
      integer :: iq

      ! Calculate a local index from the processors indices
      ! Note that face number runs from 1 here.
      iq = i + (j-1)*ipan + (n-1)*ipan*jpan
   end function indp

   subroutine checksize(len, mesg)
      integer, intent(in) :: len
      character(len=*), intent(in) :: mesg
      if ( len > maxsize ) then
         print*, "Error, maxsize exceeded in ", mesg
         stop
      end if
   end subroutine checksize

   subroutine check_bnds_alloc(rproc, iext)
      integer, intent(in) :: rproc
      integer, intent(in) :: iext
      integer :: len, ierr

!     Allocate the components of the bnds array. It's too much work to
!     get the exact sizes, so allocate a fixed size for each case where
!     there's an interaction.
      if ( bnds(rproc)%len == 0 ) then
         ! Not allocated yet.
         if ( nproc < npanels+1 ) then
            ! This is the maximum size, each face has 4 edges
            len = npan*4*(il_g+4)*2*kl
         else
            len = (max(ipan,jpan)+4)*2*kl
         end if
         allocate ( bnds(rproc)%rbuf(len) )
         allocate ( bnds(rproc)%sbuf(len) )
         allocate ( bnds(rproc)%request_list(len) )
         allocate ( bnds(rproc)%send_list(len) )
         allocate ( bnds(rproc)%unpack_list(len) )
         allocate ( bnds(rproc)%request_list_uv(len) )
         allocate ( bnds(rproc)%send_list_uv(len) )
         allocate ( bnds(rproc)%unpack_list_uv(len) )
         allocate ( bnds(rproc)%uv_swap(len), bnds(rproc)%send_swap(len) )
         bnds(rproc)%len = len
      else
         ! Just check length
         if ( kl*bnds(rproc)%rlen >=  bnds(rproc)%len ) then
            print*, "Error, maximum length error in check_bnds_alloc"
            print*, myid, rproc, bnds(rproc)%rlen,  bnds(rproc)%len, kl
            call MPI_Abort(MPI_COMM_WORLD,ierr)
         end if
         if ( iext >= iextra ) then
            print*, "Error, iext maximum length error in check_bnds_alloc"
            print*, myid, iext, iextra
            call MPI_Abort(MPI_COMM_WORLD,ierr)
         end if
      end if
   end subroutine check_bnds_alloc

   subroutine fix_index(iqx,larray,n,bnds,iext)
      integer, intent(in) :: iqx, n
      integer, dimension(:), intent(out) :: larray
      integer, intent(inout) :: iext
      type(bounds_info), dimension(0:), intent(inout) :: bnds
      integer :: rproc, iloc,jloc,nloc

      ! Which processor has this point
      rproc = qproc(iqx)
      if ( rproc /= myid ) then ! Add to list
         call check_bnds_alloc(rproc, iext)
         bnds(rproc)%rlen = bnds(rproc)%rlen + 1
         bnds(rproc)%request_list(bnds(rproc)%rlen) = iqx
         ! Increment extended region index
         iext = iext + 1
         bnds(rproc)%unpack_list(bnds(rproc)%rlen) = iext
         larray(n) = ifull+iext
      else
         ! If it's on this processor, just use the local index
         call indv_mpi(iqx,iloc,jloc,nloc)
         larray(n) = indp(iloc,jloc,nloc)
      end if
   end subroutine fix_index

   subroutine fix_index2(iqx,larray,n,bnds,iext)
      integer, intent(in) :: iqx, n
      integer, dimension(:), intent(out) :: larray
      integer, intent(inout) :: iext
      type(bounds_info), dimension(0:), intent(inout) :: bnds
      integer :: rproc, iloc,jloc,nloc

      ! Which processor has this point
      rproc = qproc(iqx)
      if ( rproc /= myid ) then ! Add to list
         call check_bnds_alloc(rproc, iext)
         bnds(rproc)%rlen2 = bnds(rproc)%rlen2 + 1
         bnds(rproc)%request_list(bnds(rproc)%rlen2) = iqx
         ! Increment extended region index
         iext = iext + 1
         bnds(rproc)%unpack_list(bnds(rproc)%rlen2) = iext
         larray(n) = ifull+iext
      else
         ! If it's on this processor, just use the local index
         call indv_mpi(iqx,iloc,jloc,nloc)
         larray(n) = indp(iloc,jloc,nloc)
      end if
   end subroutine fix_index2

   subroutine proc_setup(il,npanels,ifull)
      include 'parm.h'
!     Routine to set up offsets etc.
      integer, intent(in) :: il, npanels, ifull
      integer :: i, j, n, ierr, iproc, nd, jdf, idjd_g

      !  Processor allocation
      !  if  nproc <= npanels+1, then each gets a number of full panels
      if ( nproc <= npanels+1 ) then
         if ( modulo(npanels+1,nproc) /= 0 ) then
            print*, "Error, number of processors must divide number of panels"
            stop
         end if
!         npan = (npanels+1)/nproc
         ipan = il_g
         jpan = il_g
         noff = 1 - myid*npan
         ioff = 0
         joff = 0
         do n=0,npanels
            fproc(:,:,n) = n/npan
            qproc(1+n*il_g**2:(n+1)*il_g**2) = n/npan
         end do
      else  ! nproc >= npanels+1
         if ( modulo (nproc, npanels+1) /= 0 ) then
            print*, "Error, number of processors must be a multiple of number of panels"
            stop
         end if
!         npan = 1
         n = nproc / (npanels+1)
         !     n must be a power of 2
         nxproc = 1 ! Number of processors in x direction
         nyproc = 1 ! Number of processors in y direction
         do
            nxproc = nxproc*2
            n = n / 2
            if ( n == 1 ) exit
            nyproc = nyproc*2
            n = n / 2
            if ( n == 1 ) exit
         end do
         n = nproc / (npanels+1)
         print*, "NX, NY", nxproc, nyproc, n
         if ( nxproc*nyproc /= n ) then
            print*, "Error in splitting up faces"
            call MPI_finalize(ierr)
            stop
         end if
         if ( modulo(il_g,nxproc) /= 0 ) then
            print*, "Error, il not a multiple of nxproc"
            call MPI_finalize(ierr)
            stop
         end if
         if ( modulo(il_g,nyproc) /= 0 ) then
            print*, "Error, il not a multiple of nyproc"
            call MPI_finalize(ierr)
            stop
         end if
         ipan = il_g/nxproc
         jpan = il_g/nyproc

         iproc = 0
         qproc = -9999
         do n=0,npanels
            do j=1,il_g,jpan
               do i=1,il_g,ipan
                  fproc(i:i+ipan-1,j:j+jpan-1,n) = iproc
                  iproc = iproc + 1
               end do
            end do
         end do

         do n=0,npanels
            do j=1,il_g
               do i=1,il_g
                  qproc(ind(i,j,n)) = fproc(i,j,n)
               end do
            end do
         end do

         ! Set offsets for this processor
         call proc_region(myid,ioff,joff,noff)
      end if

!      ipfull = ipan*jpan*npan
!      iextra = 4*npan*(ipan+jpan+8)
!      print*, "ipfull, iextra", ipfull, iextra

      ! Convert standard jd to a face index
      nd = (jd-1)/il_g ! 0: to match fproc
      jdf = jd - nd*il_g
      mydiag = ( myid == fproc(id,jdf,nd) )
      ! Convert global indices to ones on this processors region
      idjd_g = id + (jd-1)*il_g
      ! Use offset to stop ifull going to zero
      idjd = 1 + modulo(idjd_g-1,ifull) ! Correct value on the appropriate processor

   end subroutine proc_setup

   subroutine proc_region(procid,ipoff,jpoff,npoff)
      ! Calculate the offsets for a given processor
      integer, intent(in) :: procid
      integer, intent(out) :: ipoff, jpoff, npoff
      integer :: myface, mtmp

      if ( nproc <= npanels+1 ) then
         npoff = 1 - procid*npan
         ipoff = 0
         jpoff = 0
      else
         myface = procid / (nxproc*nyproc)
         npoff = 1 - myface
         mtmp = procid - myface*(nxproc*nyproc)
         if ( nyproc == 1 ) then
            jpoff = 0
         else
            jpoff = mtmp/nyproc * jpan
         end if
         ipoff = modulo(mtmp,nxproc)*ipan
      end if
   end subroutine proc_region

end module cc_mpi
