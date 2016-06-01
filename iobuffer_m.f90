module iobuffer_m
   use cc_mpi
   use mpi
#ifdef usenc_mod
use netcdf
#else
use netcdf_m
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

   subroutine init_iobuffer(fid,fmode)
      integer, intent(in) :: fid,fmode

      if ( fmode==-1 ) then
         restart=.true.
      else
         restart=.false.
      end if
      call del_iobuffer(fid)

   end subroutine init_iobuffer

   subroutine del_iobuffer(fid)
      integer, intent(in) :: fid
      integer :: i
      type(iobuffer_t), pointer :: current,last

      if ( .not.useiobuffer .or. restart ) return

      call sync_iobuffer(fid)
      call flush_iobuffer(fid)

      current => head
      nullify(last)
      do while ( associated(current) )
         if ( current%fid .eq. fid ) then
            if ( allocated(current%var) ) deallocate(current%var)
            if ( allocated(current%gvar) ) deallocate(current%gvar)
            if ( allocated(current%ipack) ) deallocate(current%ipack)
            if ( allocated(current%gipack) ) deallocate(current%gipack)

            if ( associated(current%next) ) then
               if ( associated(last) ) then
                  last%next => current%next
                  deallocate(current)
                  current => last%next
               else
                  head => current%next
                  deallocate(current)
                  current => head
               end if
            else
               if ( associated(last) ) then
                  deallocate(current)
                  nullify(current,last%next)
                  tail => last
               else
                  deallocate(current)
                  nullify(current,head,tail)
               end if
            end if
            ibc=ibc-1
         else
            last => current
            current => current%next
         end if
      end do

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
      if ( vnode_myid.eq.0 ) then
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
      if ( vnode_myid.eq.0 ) then
         allocate( iobuff%gipack(ifull,istep,inproc) )
      else
         allocate( iobuff%gipack(0,0,0) )
      end if

      call MPI_Igather(iobuff%ipack,ifull*istep,MPI_INTEGER2,iobuff%gipack,ifull*istep,MPI_INTEGER2,0,comm_vnode,iobuff%request,&
                       ierr)

      if ( associated(head) ) then
         tail%next => iobuff
         tail => iobuff
      else
         head => iobuff
         tail => iobuff
      end if
      nullify(tail%next)

   end subroutine add_iobuffer_s

   subroutine sync_iobuffer(fid)
      integer, intent(in) :: fid
      integer :: i, ierr
      type(iobuffer_t), pointer :: current

      current => head
      do while ( associated(current) )
         if ( current%fid .eq. fid ) then
            call MPI_Wait( current%request, current%status, ierr )
         end if
         current => current%next
      end do

   end subroutine sync_iobuffer

   subroutine flush_iobuffer(fid)
      integer, intent(in) :: fid
      integer :: i
      type(iobuffer_t), pointer :: current

      if (vnode_myid.eq.0) then
         current => head
         do while ( associated(current) )
            if ( current%fid .eq. fid ) then
               call write_iobuffer(current)
            end if
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
