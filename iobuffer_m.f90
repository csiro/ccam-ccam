module iobuffer_m
   use cc_mpi
   use mpi
#ifdef usenc3
   use netcdf_m
#else
   use netcdf
#endif
   implicit none

   include 'parm.h'         ! Model configuration

   private

   integer, parameter :: iobuff_max=300
   type, public :: iobuffer_t
      integer :: fid
      integer :: vid
      integer :: ndims
      integer(kind=4), dimension(5) :: start, ncount
      integer(kind=4) :: request
      integer(kind=4), dimension(MPI_STATUS_SIZE) :: status
      real, allocatable, dimension(:,:) :: var
      real, allocatable, dimension(:,:,:) :: gvar
      integer(kind=2), allocatable, dimension(:,:) :: ipack
      integer(kind=2), allocatable, dimension(:,:,:) :: gipack
   end type
   type(iobuffer_t), dimension(iobuff_max), save :: iobuff
   integer, save :: ibc=0
   logical, save :: restart

   public :: init_iobuffer,add_iobuffer
   private :: del_iobuffer,sync_iobuffer,flush_iobuffer,write_iobuffer
   private :: add_iobuffer_r, add_iobuffer_s

   interface add_iobuffer
      module procedure add_iobuffer_r, add_iobuffer_s
   end interface
contains

   subroutine init_iobuffer(fmode)
      integer :: fmode

      if ( fmode==-1 ) then
         restart=.true.
      else
         restart=.false.
      end if
      call del_iobuffer

   end subroutine init_iobuffer

   subroutine del_iobuffer
      integer :: i

      if ( .not.useiobuffer .or. restart ) return

      call sync_iobuffer
      call flush_iobuffer

      do i=1,ibc
         iobuff(i)%fid=0
         iobuff(i)%vid=0
         iobuff(i)%ndims=0
         iobuff(i)%start=0
         iobuff(i)%ncount=0
         iobuff(i)%request=0
         iobuff(i)%status=0
         if ( allocated(iobuff(i)%var) ) deallocate(iobuff(i)%var)
         if ( allocated(iobuff(i)%gvar) ) deallocate(iobuff(i)%gvar)
         if ( allocated(iobuff(i)%ipack) ) deallocate(iobuff(i)%ipack)
         if ( allocated(iobuff(i)%gipack) ) deallocate(iobuff(i)%gipack)
      end do

      ibc=0
   end subroutine del_iobuffer

   subroutine add_iobuffer_r(fid,vid,ndims,ifull,istep,inproc,start,ncount,var)
      integer, intent(in) :: fid,vid,ndims,ifull,istep,inproc
      real, dimension(ifull,istep) :: var
      integer(kind=4), dimension(:) :: start, ncount
      integer :: ierr

      if ( .not.useiobuffer .or. restart ) return

      ibc=ibc+1
      if ( ibc.gt.iobuff_max ) call ccmpi_abort(-1)

      iobuff(ibc)%fid=fid
      iobuff(ibc)%vid=vid
      iobuff(ibc)%ndims=ndims
      iobuff(ibc)%start(1:ndims)=start(1:ndims)
      iobuff(ibc)%ncount(1:ndims)=ncount(1:ndims)

      allocate( iobuff(ibc)%var(ifull,istep) )
      iobuff(ibc)%var=var
      if ( myid_node.eq.0 ) then
         allocate( iobuff(ibc)%gvar(ifull,istep,inproc) )
      else
         allocate( iobuff(ibc)%gvar(0,0,0) )
      end if

      call MPI_Igather(iobuff(ibc)%var,ifull*istep,MPI_REAL,iobuff(ibc)%gvar,ifull*istep,MPI_REAL,0,comm_vnode,iobuff(ibc)%request,ierr)

   end subroutine add_iobuffer_r

   subroutine add_iobuffer_s(fid,vid,ndims,ifull,istep,inproc,start,ncount,var)
      integer, intent(in) :: fid,vid,ndims,ifull,istep,inproc
      integer(kind=2), dimension(ifull,istep) :: var
      integer(kind=4), dimension(:) :: start, ncount
      integer :: ierr

      if ( .not.useiobuffer .or. restart ) return

      ibc=ibc+1
      if ( ibc.gt.iobuff_max ) call ccmpi_abort(-1)

      iobuff(ibc)%fid=fid
      iobuff(ibc)%vid=vid
      iobuff(ibc)%ndims=ndims
      iobuff(ibc)%start(1:ndims)=start(1:ndims)
      iobuff(ibc)%ncount(1:ndims)=ncount(1:ndims)

      allocate( iobuff(ibc)%ipack(ifull,istep) )
      iobuff(ibc)%ipack=var
      if ( myid_node.eq.0 ) then
         allocate( iobuff(ibc)%gipack(ifull,istep,inproc) )
      else
         allocate( iobuff(ibc)%gipack(0,0,0) )
      end if

      call MPI_Igather(iobuff(ibc)%ipack,ifull*istep,MPI_INTEGER2,iobuff(ibc)%gipack,ifull*istep,MPI_INTEGER2,0,comm_vnode,iobuff(ibc)%request,ierr)

   end subroutine add_iobuffer_s

   subroutine sync_iobuffer(idx)
      integer :: i, ierr
      integer, intent(in), optional :: idx

      if ( present(idx) ) then
         call MPI_Wait( iobuff(idx)%request, iobuff(idx)%status, ierr )
      else
         do i=1,ibc
            call MPI_Wait( iobuff(i)%request, iobuff(i)%status, ierr )
         end do
      end if

   end subroutine sync_iobuffer

   subroutine flush_iobuffer(idx)
      integer :: i
      integer, intent(in), optional :: idx

      if (myid_node.eq.0) then
         if ( present(idx) ) then
            call write_iobuffer(idx)
         else
            do i=1,ibc
               call write_iobuffer(i)
            end do
         end if
      end if
   end subroutine flush_iobuffer

   subroutine write_iobuffer(idx)
      integer :: ier, lndims
      integer, intent(in) :: idx

      lndims=iobuff(idx)%ndims
      if ( allocated(iobuff(idx)%gipack) ) then
         ier = nf90_put_var(iobuff(idx)%fid,iobuff(idx)%vid,iobuff(idx)%gipack,iobuff(idx)%start(1:lndims),iobuff(idx)%ncount(1:lndims))
      elseif ( allocated(iobuff(idx)%gvar) ) then
         ier = nf90_put_var(iobuff(idx)%fid,iobuff(idx)%vid,iobuff(idx)%gvar,iobuff(idx)%start(1:lndims),iobuff(idx)%ncount(1:lndims))
      end if

   end subroutine write_iobuffer

end module iobuffer_m
