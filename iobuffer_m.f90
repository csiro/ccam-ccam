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
      type(iobuffer_t), pointer :: next => null()
   end type
   type(iobuffer_t), pointer, save :: head,tail
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
      type(iobuffer_t), pointer :: current

      if ( .not.useiobuffer .or. restart ) return

      call sync_iobuffer
      call flush_iobuffer

      do while ( associated(head) )
         if ( allocated(head%var) ) deallocate(head%var)
         if ( allocated(head%gvar) ) deallocate(head%gvar)
         if ( allocated(head%ipack) ) deallocate(head%ipack)
         if ( allocated(head%gipack) ) deallocate(head%gipack)
         current => head
         if ( associated(current%next) ) then
            head => current%next
         else
            nullify(head)
         end if
         deallocate(current)
         nullify(current)
         ibc=ibc-1
      end do
      nullify(head,tail)

   end subroutine del_iobuffer

   subroutine add_iobuffer_r(fid,vid,ndims,ifull,istep,inproc,start,ncount,var)
      integer, intent(in) :: fid,vid,ndims,ifull,istep,inproc
      real, dimension(ifull,istep) :: var
      integer(kind=4), dimension(:) :: start, ncount
      integer :: ierr
      type(iobuffer_t), pointer :: iobuff

      if ( .not.useiobuffer .or. restart ) return

      allocate(iobuff)
      ibc=ibc+1

      iobuff%fid=fid
      iobuff%vid=vid
      iobuff%ndims=ndims
      iobuff%start(1:ndims)=start(1:ndims)
      iobuff%ncount(1:ndims)=ncount(1:ndims)

      allocate( iobuff%var(ifull,istep) )
      iobuff%var=var
      if ( myid_node.eq.0 ) then
         allocate( iobuff%gvar(ifull,istep,inproc) )
      else
         allocate( iobuff%gvar(0,0,0) )
      end if

      call MPI_Igather(iobuff%var,ifull*istep,MPI_REAL,iobuff%gvar,ifull*istep,MPI_REAL,0,comm_vnode,iobuff%request,ierr)

      if ( associated(head) ) then
         tail%next => iobuff
         tail => iobuff
      else
         head => iobuff
         tail => iobuff
      end if
      nullify(tail%next)

   end subroutine add_iobuffer_r

   subroutine add_iobuffer_s(fid,vid,ndims,ifull,istep,inproc,start,ncount,var)
      integer, intent(in) :: fid,vid,ndims,ifull,istep,inproc
      integer(kind=2), dimension(ifull,istep) :: var
      integer(kind=4), dimension(:) :: start, ncount
      integer :: ierr
      type(iobuffer_t), pointer :: iobuff

      if ( .not.useiobuffer .or. restart ) return

      allocate(iobuff)
      ibc=ibc+1

      iobuff%fid=fid
      iobuff%vid=vid
      iobuff%ndims=ndims
      iobuff%start(1:ndims)=start(1:ndims)
      iobuff%ncount(1:ndims)=ncount(1:ndims)

      allocate( iobuff%ipack(ifull,istep) )
      iobuff%ipack=var
      if ( myid_node.eq.0 ) then
         allocate( iobuff%gipack(ifull,istep,inproc) )
      else
         allocate( iobuff%gipack(0,0,0) )
      end if

      call MPI_Igather(iobuff%ipack,ifull*istep,MPI_INTEGER2,iobuff%gipack,ifull*istep,MPI_INTEGER2,0,comm_vnode,iobuff%request,ierr)

      if ( associated(head) ) then
         tail%next => iobuff
         tail => iobuff
      else
         head => iobuff
         tail => iobuff
      end if
      nullify(tail%next)

   end subroutine add_iobuffer_s

   subroutine sync_iobuffer
      integer :: i, ierr
      type(iobuffer_t), pointer :: current

      current => head
      do while ( associated(current) )
         call MPI_Wait( current%request, current%status, ierr )
         current => current%next
      end do

   end subroutine sync_iobuffer

   subroutine flush_iobuffer
      integer :: i
      type(iobuffer_t), pointer :: current

      if (myid_node.eq.0) then
         current => head
         do while ( associated(current) )
            call write_iobuffer(current)
            current => current%next
         end do
      end if
   end subroutine flush_iobuffer

   subroutine write_iobuffer(item)
      type(iobuffer_t), pointer :: item
      integer :: ier, lndims

      lndims=item%ndims
      if ( allocated(item%gipack) ) then
         ier = nf90_put_var(item%fid,item%vid,item%gipack,item%start(1:lndims),item%ncount(1:lndims))
      elseif ( allocated(item%gvar) ) then
         ier = nf90_put_var(item%fid,item%vid,item%gvar,item%start(1:lndims),item%ncount(1:lndims))
      end if

   end subroutine write_iobuffer

end module iobuffer_m
