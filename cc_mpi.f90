module cc_mpi
   implicit none
   private
   include 'newmpar.h'
   include 'mpif.h'

   integer, save, public :: myid ! Processor number
   integer, save, public :: ipan, jpan
   integer, save, dimension(0:npanels), public :: ioff, joff ! these can change on different panels
   integer, save, public :: noff
   integer, save, public :: comm_world
   integer, save, private :: nxproc, nyproc

   integer, public, allocatable, save, dimension(:,:,:) :: fproc
   integer, public, allocatable, save, dimension(:) :: qproc


   ! These only need to be module variables for the overlap case
   ! when they have to keep values between boundsa and boundsb.
   ! For now leave them here.
   integer(kind=4), allocatable, dimension(:), save, private :: ireq
   integer(kind=4), save, private :: nreq

   public :: bounds, boundsuv, ccmpi_setup, ccmpi_distribute, ccmpi_gather, &
             ccmpi_distributer8, indp, indg, deptsync, intssync, start_log, &
             end_log, log_on, log_off, log_setup, phys_loadbal,             &
             ccglobal_posneg, ccglobal_sum, iq2iqg, indv_mpi, indglobal,    &
             readglobvar, writeglobvar, face_set, uniform_set,              &
             ccmpi_reduce, ccmpi_allreduce, ccmpi_abort, ccmpi_bcast,       &
             ccmpi_bcastr8, ccmpi_barrier, ccmpi_gatherx, ccmpi_scatterx,   &
             ccmpi_allgatherx, ccmpi_recv, ccmpi_ssend, ccmpi_init,         &
             ccmpi_finalize, ccmpi_commsplit, ccmpi_commfree
   public :: dpoints_t,dindex_t,sextra_t,bnds
   private :: ccmpi_distribute2, ccmpi_distribute2i, ccmpi_distribute2r8,   &
              ccmpi_distribute3, ccmpi_distribute3i, ccmpi_gather2,         &
              ccmpi_gather3, checksize, ccglobal_posneg2, ccglobal_posneg3, &
              ccglobal_sum2, ccglobal_sum3
   interface ccmpi_gather
      module procedure ccmpi_gather2, ccmpi_gather3
   end interface
   interface ccmpi_distribute
      module procedure ccmpi_distribute2, ccmpi_distribute2i,  &    
                       ccmpi_distribute3, ccmpi_distribute3i
   end interface
   interface ccmpi_distributer8
      module procedure ccmpi_distribute2r8
   end interface
   interface bounds
      module procedure bounds2, bounds3
   end interface
   interface boundsuv
      module procedure boundsuv2, boundsuv3
   end interface
   interface ccglobal_posneg
      module procedure ccglobal_posneg2, ccglobal_posneg3
   end interface
   interface ccglobal_sum
      module procedure ccglobal_sum2, ccglobal_sum3
   end interface
   interface readglobvar
      module procedure readglobvar2, readglobvar3, readglobvar2i
   end interface
   interface writeglobvar
      module procedure writeglobvar2, writeglobvar3
   end interface
   interface ccmpi_reduce
      module procedure ccmpi_reduce2i, ccmpi_reduce2r, ccmpi_reduce3r, ccmpi_reduce2c
   end interface ccmpi_reduce
   interface ccmpi_allreduce
      module procedure ccmpi_allreduce2i, ccmpi_allreduce2r, ccmpi_allreduce3r, ccmpi_allreduce2c
   end interface ccmpi_allreduce
   interface ccmpi_bcast
      module procedure ccmpi_bcast2i, ccmpi_bcast3i, ccmpi_bcast2r, ccmpi_bcast3r, ccmpi_bcast4r, &
                       ccmpi_bcast5r, ccmpi_bcast2s
   end interface ccmpi_bcast
   interface ccmpi_bcastr8
      module procedure ccmpi_bcast2r8, ccmpi_bcast3r8, ccmpi_bcast4r8
   end interface ccmpi_bcastr8
   interface ccmpi_gatherx
      module procedure ccmpi_gatherx2r
   end interface ccmpi_gatherx
   interface ccmpi_scatterx
      module procedure ccmpi_scatterx2r
   end interface ccmpi_scatterx
   interface ccmpi_allgatherx
      module procedure ccmpi_allgatherx2i
   end interface ccmpi_allgatherx
   interface ccmpi_recv
      module procedure ccmpi_recv2r
   end interface ccmpi_recv
   interface ccmpi_ssend
      module procedure ccmpi_ssend2r
   end interface ccmpi_ssend

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
      logical, dimension(:), pointer :: uv_swap, send_swap,uv_neg,send_neg
      ! Number of points for each processor. Also double row versions.
      ! lenx is first row plux corner points.  lenh is just the ne side.
      integer(kind=4) :: slen, rlen, slenx, rlenx, slen2, rlen2
      integer(kind=4) :: slenh, rlenh
      ! Number of points for each processor. lenx is for nu, su, ev and wv
      integer(kind=4) :: slen_uv, rlen_uv, slen2_uv, rlen2_uv
      integer(kind=4) :: slenx_uv, rlenx_uv
      integer(kind=4) :: len
      ! ocean mask
      integer(kind=4) :: mlomsk
   end type bounds_info
   
   type dpoints_t
      real, dimension(:,:), allocatable :: a
      real, dimension(:), allocatable :: b
   end type dpoints_t
   type dindex_t
      integer, dimension(:,:), allocatable :: a
   end type dindex_t
   type sextra_t
      real, dimension(:), allocatable :: a
   end type sextra_t

   type boundsplit
      integer :: isubg, iwvfn, inubg, ievfn
      integer :: isvbg, iwufn, invbg, ieufn
      integer :: issvbg, iwwufn, innvbg, ieeufn
   end type boundsplit
   type coloursplit
      integer, dimension(3) :: ihbg, ihfn, ifbg, iffn
   end type coloursplit

   type(bounds_info), allocatable, dimension(:), save :: bnds
   integer, save, public :: mlomsk_ne, mlomsk_sw

   type(boundsplit), allocatable, dimension(:), save, private :: rsplit
   type(boundsplit), allocatable, dimension(:), save, private :: ssplit
   type(coloursplit), allocatable, dimension(:), save, private :: rcolsp
   type(coloursplit), allocatable, dimension(:), save, private :: scolsp
   integer, allocatable, dimension(:), save, public :: colourmask
   integer, dimension(3), save, public :: ifullx
   integer, dimension(:,:), allocatable, save, public :: iqx,iqn,iqe,iqw,iqs
   integer, dimension(:,:), allocatable, save, public :: iqwu,iqsv
#ifdef uniform_decomp
   integer, parameter, public :: maxcolour = 3
#else
   integer, parameter, public :: maxcolour = 2
#endif

   integer, public, save :: maxbuflen

   ! Flag whether processor region edge is a face edge.
   logical, dimension(0:npanels), public, save :: edge_w, edge_n, edge_s, edge_e

   ! Off processor departure points
   type(dpoints_t), allocatable, dimension(:), public, save :: dpoints
   type(dpoints_t), allocatable, dimension(:), private, save :: dbuf
   type(dindex_t), allocatable, dimension(:), public, save :: dindex
   type(sextra_t), allocatable, dimension(:), public, save :: sextra
   ! Number of points for each processor.
   integer(kind=4), dimension(:), allocatable, public, save :: dslen, drlen

   ! True if processor is a nearest neighbour
   logical, allocatable, dimension(:), public, save :: neighbour

   logical, public, save :: mydiag ! True if diagnostic point id, jd is in my region

   integer, public, save :: bounds_begin, bounds_end
   integer, public, save :: boundsa_begin, boundsa_end
   integer, public, save :: boundsb_begin, boundsb_end
   integer, public, save :: boundsuv_begin, boundsuv_end
   integer, public, save :: ints_begin, ints_end
   integer, public, save :: nonlin_begin, nonlin_end
   integer, public, save :: helm_begin, helm_end
   integer, public, save :: adjust_begin, adjust_end
   integer, public, save :: upglobal_begin, upglobal_end
   integer, public, save :: hordifg_begin, hordifg_end
   integer, public, save :: vadv_begin, vadv_end
   integer, public, save :: depts_begin, depts_end
   integer, public, save :: deptsync_begin, deptsync_end
   integer, public, save :: intssync_begin, intssync_end
   integer, public, save :: stag_begin, stag_end
   integer, public, save :: ocnstag_begin, ocnstag_end
   integer, public, save :: toij_begin, toij_end
   integer, public, save :: physloadbal_begin, physloadbal_end
   integer, public, save :: phys_begin, phys_end
   integer, public, save :: outfile_begin, outfile_end
   integer, public, save :: onthefly_begin, onthefly_end
   integer, public, save :: indata_begin, indata_end
   integer, public, save :: nestin_begin, nestin_end
   integer, public, save :: gwdrag_begin, gwdrag_end
   integer, public, save :: convection_begin, convection_end
   integer, public, save :: cloud_begin, cloud_end
   integer, public, save :: radnet_begin,radnet_end
   integer, public, save :: radmisc_begin,radmisc_end
   integer, public, save :: radsw_begin, radsw_end
   integer, public, save :: radlw_begin, radlw_end   
   integer, public, save :: sfluxnet_begin, sfluxnet_end
   integer, public, save :: sfluxwater_begin, sfluxwater_end
   integer, public, save :: sfluxland_begin, sfluxland_end
   integer, public, save :: sfluxurban_begin, sfluxurban_end
   integer, public, save :: vertmix_begin, vertmix_end
   integer, public, save :: aerosol_begin, aerosol_end
   integer, public, save :: model_begin, model_end
   integer, public, save :: maincalc_begin, maincalc_end
   integer, public, save :: gather_begin, gather_end
   integer, public, save :: distribute_begin, distribute_end
   integer, public, save :: posneg_begin, posneg_end
   integer, public, save :: globsum_begin, globsum_end
   integer, public, save :: precon_begin, precon_end
   integer, public, save :: waterdynamics_begin, waterdynamics_end
   integer, public, save :: waterdiff_begin, waterdiff_end
   integer, public, save :: river_begin, river_end
   integer, public, save :: mpiwait_begin, mpiwait_end
   integer, public, save :: reduce_begin, reduce_end
   integer, public, save :: bcast_begin, bcast_end
#ifdef simple_timer
   public :: simple_timer_finalize
   integer(kind=4), parameter :: nevents=49
   double precision, dimension(nevents), save :: tot_time = 0., start_time
   character(len=15), dimension(nevents), save :: event_name
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
      use indices_m
      use latlong_m
      use map_m
      use sumdd_m
      use vecsuv_m
      use xyzinfo_m
      integer iproc, iq, iqg, i, j, n, ig, jg, ng, tg, jx, mcc
      integer(kind=4) ierr

      allocate(fproc(il_g,il_g,0:npanels))
      allocate(qproc(ifull_g))
      allocate(ireq(2*nproc))
      allocate(bnds(0:nproc-1))
      allocate(dpoints(0:nproc-1))
      allocate(dbuf(0:nproc-1))
      allocate(dindex(0:nproc-1))
      allocate(sextra(0:nproc-1))
      allocate(neighbour(0:nproc-1))
      
#ifdef uniform_decomp
      call proc_setup_uniform(npanels,ifull)
      ! Faces may not line up properly so need extra factor here
      maxbuflen = (max(ipan,jpan)+4)*3*max(kl+1,ol+1) * 8 * 2  !*3 for extra vector row (e.g., inu,isu,iev,iwv)
                                                               ! kl+1 and ol+1 for 2D + 3D packing
#else
      call proc_setup(npanels,ifull)
      if ( nproc < npanels+1 ) then
         ! This is the maximum size, each face has 4 edges
         maxbuflen = npan*4*(il_g+4)*3*max(kl+1,ol+1) !*3 for extra vector row (e.g., inu,isu,iev,iwv)
                                                      ! kl+1 and ol+1 for 2D + 3D packing
      else
         maxbuflen = (max(ipan,jpan)+4)*3*max(kl+1,ol+1) !*3 for extra vector row (e.g., inu,isu,iev,iwv)
                                                         ! kl+1 and ol+1 for 2D + 3D packing
      end if
#endif

#ifdef DEBUG
      write(6,*) "Grid", npan, ipan, jpan
      write(6,*) "Offsets", myid, ioff(:), joff(:), noff
#endif

      allocate ( rsplit(0:nproc-1), ssplit(0:nproc-1) )
      allocate ( rcolsp(0:nproc-1), scolsp(0:nproc-1) )
      allocate ( colourmask(ifull) )

      ! Also do the initialisation for deptsync here
      allocate ( dslen(0:nproc-1), drlen(0:nproc-1) )
      dslen=0
      drlen=0

#ifdef uniform_decomp
      ! three colour mask
      do n=1,npan
         do j=1,jpan
            do i=1,ipan
               iq = indp(i,j,n)   ! Local
               iqg = indg(i,j,n)  ! Global
           
               ! calculate global i,j,n
               tg = iqg-1
               ng = tg/(il_g*il_g)
               tg = tg-ng*il_g*il_g
               jg = tg/il_g
               tg = tg-jg*il_g
               ig = tg
               ig = ig+1
               jg = jg+1
           
               jx = mod(ig+jg+ng*il_g,2)
               select case( ng+jx*6 )
                  case(0,1,3,4)
                     colourmask(iq)=1
                  case(2,5,6,9)
                     colourmask(iq)=2
                  case(7,8,10,11)
                     colourmask(iq)=3
               end select
            end do
         end do
      end do
#else
      ! two colour mask
      do n=1,npan
         do j=1,jpan
            do i=1,ipan
               iq = indp(i,j,n)   ! Local
               iqg = indg(i,j,n)  ! Global
           
               ! calculate global i,j,n
               tg = iqg-1
               ng = tg/(il_g*il_g)
               tg = tg-ng*il_g*il_g
               jg = tg/il_g
               tg = tg-jg*il_g
               ig = tg
               ig = ig+1
               jg = jg+1
           
               jx = mod(ig+jg+ng*il_g,2)
               if (jx == 0 ) then
                  colourmask(iq)=1
               else
                  colourmask(iq)=2
               end if
            end do
         end do
      end do
#endif
      
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
         call ccmpi_distribute(f,f_g)
         call ccmpi_distribute(fu,fu_g)
         call ccmpi_distribute(fv,fv_g)
         call ccmpi_distribute(dmdx,dmdx_g)
         call ccmpi_distribute(dmdy,dmdy_g)
         call ccmpi_distributer8(x,x_g)
         call ccmpi_distributer8(y,y_g)
         call ccmpi_distributer8(z,z_g)
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
         call ccmpi_distribute(f)
         call ccmpi_distribute(fu)
         call ccmpi_distribute(fv)
         call ccmpi_distribute(dmdx)
         call ccmpi_distribute(dmdy)
         call ccmpi_distributer8(x)
         call ccmpi_distributer8(y)
         call ccmpi_distributer8(z)
         call ccmpi_distribute(rlatt)
         call ccmpi_distribute(rlongg)
      end if

      call bounds_setup()
      call bounds(em)
      call boundsuv(emu,emv)
      call boundsuv(ax,bx)
      call boundsuv(ay,by)
      call boundsuv(az,bz)
      call bounds(f,corner=.true.)

      ! Off processor departure points
      do iproc=0,nproc-1
        if (neighbour(iproc)) then
          allocate(dpoints(iproc)%a(4,bnds(iproc)%len))
          allocate(dbuf(iproc)%a(4,bnds(iproc)%len))
          allocate(dbuf(iproc)%b(bnds(iproc)%len))
          allocate(sextra(iproc)%a(bnds(iproc)%len))
          allocate(dindex(iproc)%a(2,bnds(iproc)%len))
        else
          allocate(dpoints(iproc)%a(4,1))
          allocate(dbuf(iproc)%a(4,1))
          allocate(dbuf(iproc)%b(1))
          allocate(sextra(iproc)%a(1))
          allocate(dindex(iproc)%a(2,1))
        end if
      end do

      ! Pack colour indices
      mcc=max( count( colourmask == 1 ), count( colourmask == 2), count( colourmask == 3) )
      allocate ( iqx(mcc,maxcolour) )
      allocate ( iqn(mcc,maxcolour), iqe(mcc,maxcolour) )
      allocate ( iqw(mcc,maxcolour), iqs(mcc,maxcolour) )
      allocate ( iqwu(mcc,maxcolour), iqsv(mcc,maxcolour) )
      ifullx=0
      iqx=0
      iqn=0
      iqe=0
      iqw=0
      iqs=0
      iqwu=0
      iqsv=0
      do iq=1,ifull
         ifullx(colourmask(iq))=ifullx(colourmask(iq))+1
         iqx(ifullx(colourmask(iq)),colourmask(iq))=iq
         iqn(ifullx(colourmask(iq)),colourmask(iq))=in(iq)
         iqe(ifullx(colourmask(iq)),colourmask(iq))=ie(iq)
         iqw(ifullx(colourmask(iq)),colourmask(iq))=iw(iq)
         iqs(ifullx(colourmask(iq)),colourmask(iq))=is(iq)
         iqwu(ifullx(colourmask(iq)),colourmask(iq))=iwu(iq)
         iqsv(ifullx(colourmask(iq)),colourmask(iq))=isv(iq)
      end do

#ifdef sumdd
!     operator MPI_SUMDR is created based on an external function DRPDR. 
      call MPI_OP_CREATE (DRPDR, .TRUE., MPI_SUMDR, ierr) 
#endif

   end subroutine ccmpi_setup

   subroutine ccmpi_distribute2(af,a1)
      ! Convert standard 1D arrays to face form and distribute to processors
      real, dimension(ifull), intent(out) :: af
      real, dimension(ifull_g), intent(in), optional :: a1
      integer :: i, j, n, iq, iproc
      integer(kind=4) ierr, lsize, ltype
      real, dimension(ifull,0:nproc-1) :: sbuf
      integer :: npoff, ipoff, jpoff ! Offsets for target
      integer :: slen

      call start_log(distribute_begin)
!cdir iexpand(indp, indg)
      if ( myid == 0 .and. .not. present(a1) ) then
         write(6,*) "Error: ccmpi_distribute argument required on proc 0"
         call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
      end if
      ! Copy internal region
      if ( myid == 0 ) then
         ! map array in order of processor rank
         do iproc=0,nproc-1
            slen = 0
#ifdef uniform_decomp
            do n=1,npan
               ! Panel range on the source processor
               call proc_region(iproc,n-1,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan,il_g,nproc,1)
#else
            call proc_region(iproc,0,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan,il_g,nproc,0)
            do n=1,npan
#endif
               do j=1,jpan
                  do i=1,ipan
                     iq = i+ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g
                     slen = slen + 1
                     sbuf(slen,iproc) = a1(iq)
                  end do
               end do
            end do
         end do
      end if
      lsize=ifull
#ifdef r8i8
      ltype=MPI_DOUBLE_PRECISION
#else
      ltype=MPI_REAL
#endif
      call MPI_Scatter(sbuf,lsize,ltype,af,lsize,ltype,0,MPI_COMM_WORLD,ierr)

      call end_log(distribute_end)
   end subroutine ccmpi_distribute2

   subroutine ccmpi_distribute2r8(af,a1)
      ! Convert standard 1D arrays to face form and distribute to processors
      real*8, dimension(ifull), intent(out) :: af
      real*8, dimension(ifull_g), intent(in), optional :: a1
      integer :: i, j, n, iq, iproc
      integer(kind=4) ierr,lsize
      real*8, dimension(ifull,0:nproc-1) :: sbuf
      integer :: npoff, ipoff, jpoff ! Offsets for target
      integer :: slen

      call start_log(distribute_begin)
!cdir iexpand(indp, indg)
      if ( myid == 0 .and. .not. present(a1) ) then
         write(6,*) "Error: ccmpi_distribute argument required on proc 0"
         call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
      end if
      ! Copy internal region
      if ( myid == 0 ) then
         ! map array in order of processor rank
         do iproc=0,nproc-1
            slen = 0
#ifdef uniform_decomp
            do n=1,npan
               ! Panel range on the source processor
               call proc_region(iproc,n-1,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan,il_g,nproc,1)
#else
            call proc_region(iproc,0,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan,il_g,nproc,0)
            do n=1,npan
#endif
               do j=1,jpan
                  do i=1,ipan
                     iq = i+ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g
                     slen = slen + 1
                     sbuf(slen,iproc) = a1(iq)
                  end do
               end do
            end do
         end do
      end if
      lsize=ifull
      call MPI_Scatter(sbuf,lsize,MPI_DOUBLE_PRECISION,af,lsize,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      call end_log(distribute_end)
   end subroutine ccmpi_distribute2r8

   subroutine ccmpi_distribute2i(af,a1)
      ! Convert standard 1D arrays to face form and distribute to processors
      integer, dimension(ifull), intent(out) :: af
      integer, dimension(ifull_g), intent(in), optional :: a1
      integer :: i, j, n, iq, iproc
      integer(kind=4) ierr,lsize,ltype
      integer, dimension(ifull,0:nproc-1) :: sbuf
      integer :: npoff, ipoff, jpoff ! Offsets for target
      integer :: slen

      call start_log(distribute_begin)
!cdir iexpand(indp, indg)
      if ( myid == 0 .and. .not. present(a1) ) then
         write(6,*) "Error: ccmpi_distribute argument required on proc 0"
         call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
      end if
      ! Copy internal region
      if ( myid == 0 ) then
         ! map array in order of processor rank
         do iproc=0,nproc-1
            slen = 0
#ifdef uniform_decomp
            do n=1,npan
               ! Panel range on the source processor
               call proc_region(iproc,n-1,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan,il_g,nproc,1)
#else
            call proc_region(iproc,0,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan,il_g,nproc,0)
            do n=1,npan
#endif
               do j=1,jpan
                  do i=1,ipan
                     iq = i+ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g
                     slen = slen + 1
                     sbuf(slen,iproc) = a1(iq)
                  end do
               end do
            end do
         end do
      end if
      lsize=ifull
#ifdef r8i8
      ltype=MPI_INTEGER8
#else
      ltype=MPI_INTEGER
#endif
      call MPI_Scatter(sbuf,lsize,ltype,af,lsize,ltype,0,MPI_COMM_WORLD,ierr)
 
      call end_log(distribute_end)
   end subroutine ccmpi_distribute2i

   subroutine ccmpi_distribute3(af,a1)
      ! Convert standard 1D arrays to face form and distribute to processors
      ! This is also used for tracers, so second dimension is not necessarily
      ! the number of levels
      real, dimension(:,:), intent(out) :: af
      real, dimension(:,:), intent(in), optional :: a1
      integer :: i, j, n, iq, iproc
      integer(kind=4) ierr,lsize,ltype
      real, dimension(ifull,size(af,2),0:nproc-1) :: sbuf
      integer :: npoff, ipoff, jpoff ! Offsets for target
      integer :: slen, kx

      call start_log(distribute_begin)
!cdir iexpand(indp, indg)
      if ( myid == 0 .and. .not. present(a1) ) then
         write(6,*) "Error: ccmpi_distribute argument required on proc 0"
         call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
      end if
      kx=size(af,2)
      ! Copy internal region
      if ( myid == 0 ) then
         ! map array in order of processor rank
         do iproc=0,nproc-1
            slen = 0
#ifdef uniform_decomp
            do n=1,npan
               ! Panel range on the source processor
               call proc_region(iproc,n-1,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan,il_g,nproc,1)
#else
            call proc_region(iproc,0,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan,il_g,nproc,0)
            do n=1,npan
#endif
               do j=1,jpan
                  do i=1,ipan
                     iq = i+ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g
                     slen = slen + 1
                     sbuf(slen,:,iproc) = a1(iq,:)
                  end do
               end do
            end do
         end do
      end if
      lsize=ifull*kx
#ifdef r8i8
      ltype=MPI_DOUBLE_PRECISION
#else
      ltype=MPI_REAL
#endif
      call MPI_Scatter(sbuf,lsize,ltype,af,lsize,ltype,0,MPI_COMM_WORLD,ierr)      

      call end_log(distribute_end)
   end subroutine ccmpi_distribute3

   subroutine ccmpi_distribute3i(af,a1)
      ! Convert standard 1D arrays to face form and distribute to processors
      ! This is also used for tracers, so second dimension is not necessarily
      ! the number of levels
      integer, dimension(:,:), intent(out) :: af
      integer, dimension(:,:), intent(in), optional :: a1
      integer :: i, j, n, iq, iproc
      integer(kind=4) ierr,lsize,ltype
      integer, dimension(ifull,size(af,2),0:nproc-1) :: sbuf
      integer :: npoff, ipoff, jpoff ! Offsets for target
      integer :: slen, kx

      call start_log(distribute_begin)
!cdir iexpand(indp, indg)
      if ( myid == 0 .and. .not. present(a1) ) then
         write(6,*) "Error: ccmpi_distribute argument required on proc 0"
         call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
      end if
      kx=size(af,2)
      ! Copy internal region
      if ( myid == 0 ) then
         ! map array in order of processor rank
         do iproc=0,nproc-1
            slen = 0
#ifdef uniform_decomp
            do n=1,npan
               ! Panel range on the source processor
               call proc_region(iproc,n-1,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan,il_g,nproc,1)
#else
            call proc_region(iproc,0,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan,il_g,nproc,0)
            do n=1,npan
#endif
               do j=1,jpan
                  do i=1,ipan
                     iq = i+ipoff + (j+jpoff-1)*il_g + (n-npoff)*il_g*il_g
                     slen = slen + 1
                     sbuf(slen,:,iproc) = a1(iq,:)
                  end do
               end do
            end do
         end do
      end if
      lsize=ifull*kx
#ifdef r8i8
      ltype=MPI_INTEGER8
#else
      ltype=MPI_INTEGER
#endif
      call MPI_Scatter(sbuf,lsize,ltype,af,lsize,ltype,0,MPI_COMM_WORLD,ierr)      

      call end_log(distribute_end)
   end subroutine ccmpi_distribute3i

   subroutine ccmpi_gather2(a,ag)
      ! Collect global arrays.

      real, dimension(ifull), intent(in) :: a
      real, dimension(ifull_g), intent(out), optional :: ag
      integer :: iproc
      integer(kind=4) ierr,lsize,ltype
      real, dimension(ifull,0:nproc-1) :: abuf
      integer :: ipoff, jpoff, npoff
      integer :: i, j, n, iq, iqg

      call start_log(gather_begin)
!cdir iexpand(indp, indg, ind)
      if ( myid == 0 .and. .not. present(ag) ) then
         write(6,*) "Error: ccmpi_gather argument required on proc 0"
         call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
      end if
      
      lsize=ifull
#ifdef r8i8
      ltype=MPI_DOUBLE_PRECISION
#else
      ltype=MPI_REAL
#endif
      call MPI_Gather(a,lsize,ltype,abuf,lsize,ltype,0,MPI_COMM_WORLD,ierr)
      if ( myid == 0 ) then
         ! map array in order of processor rank
         do iproc=0,nproc-1
#ifdef uniform_decomp
            do n=1,npan
               ! Panel range on the source processor
               call proc_region(iproc,n-1,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan,il_g,nproc,1)
#else
            call proc_region(iproc,0,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan,il_g,nproc,0)
            do n=1,npan
#endif
               ! Use the face indices for unpacking
               do j=1,jpan
!cdir nodep
                  do i=1,ipan
                     ! Global indices are i+ipoff, j+jpoff, n-npoff
                     iqg = indglobal(i+ipoff,j+jpoff,n-npoff) ! True global 1D index
                     iq = indp(i,j,n)
                     ag(iqg) = abuf(iq,iproc)
                  end do
               end do
            end do
         end do
      end if

      call end_log(gather_end)

   end subroutine ccmpi_gather2

   subroutine ccmpi_gather3(a,ag)
      ! Collect global arrays.

      real, dimension(:,:), intent(in) :: a
      real, dimension(:,:), intent(out), optional :: ag
      integer :: iproc
      integer(kind=4) ierr,lsize,ltype
      real, dimension(ifull,size(a,2),0:nproc-1) :: abuf
      integer :: ipoff, jpoff, npoff
      integer :: i, j, n, iq, iqg, kx

      call start_log(gather_begin)
!cdir iexpand(indp, indg, ind)
      if ( myid == 0 .and. .not. present(ag) ) then
         write(6,*) "Error: ccmpi_gather argument required on proc 0"
         call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
      end if

      kx=size(a,2)
      lsize=ifull*kx
#ifdef r8i8
      ltype=MPI_DOUBLE_PRECISION
#else
      ltype=MPI_REAL
#endif
      call MPI_Gather(a,lsize,ltype,abuf,lsize,ltype,0,MPI_COMM_WORLD,ierr)
      if ( myid == 0 ) then
         ! map array in order of processor rank
         do iproc=0,nproc-1
#ifdef uniform_decomp
            do n=1,npan
               ! Panel range on the source processor
               call proc_region(iproc,n-1,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan,il_g,nproc,1)
#else
            call proc_region(iproc,0,ipoff,jpoff,npoff,nxproc,nyproc,ipan,jpan,npan,il_g,nproc,0)
            do n=1,npan
#endif
               ! Use the face indices for unpacking
               do j=1,jpan
!cdir nodep
                  do i=1,ipan
                     ! Global indices are i+ipoff, j+jpoff, n-npoff
                     iqg = indglobal(i+ipoff,j+jpoff,n-npoff) ! True global 1D index
                     iq = indp(i,j,n)
                     ag(iqg,:) = abuf(iq,:,iproc)
                  end do
               end do
            end do
         end do
      end if
      
      call end_log(gather_end)

   end subroutine ccmpi_gather3

   subroutine bounds_setup()

      use indices_m
      integer :: n, nr, i, j, iq, iqq
      integer :: iproc
      integer(kind=4), dimension(MPI_STATUS_SIZE,2*nproc) :: status
      integer(kind=4) :: ierr, itag=0, count, rproc, sproc
      integer(kind=4), dimension(8,0:nproc-1) :: dums, dumr
      integer(kind=4), dimension(9,0:nproc-1) :: dumsb, dumrb
      integer :: iqg, iql, iloc, jloc, nloc, icol
      integer :: iext, iextu, iextv
      logical :: swap

      ! Just set values that point to values within own processors region.
      ! Other values are set up later
      ! Set initial values to make sure missing values are caught properly.
!cdir iexpand(indp, indg, indv_mpi)
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
      do n=1,npan
         do j=1,jpan
            do i=1,ipan
               iq = indp(i,j,n)   ! Local
               iqg = indg(i,j,n)  ! Global

               iqq = in_g(iqg)    ! Global neighbour index
               rproc = qproc(iqq) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  call indv_mpi(iqq,iloc,jloc,nloc)
                  in(iq) = indp(iloc,jloc,nloc)
               end if
               iqq = inn_g(iqg)    ! Global neighbour index
               rproc = qproc(iqq) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  call indv_mpi(iqq,iloc,jloc,nloc)
                  inn(iq) = indp(iloc,jloc,nloc)
               end if

               iqq = is_g(iqg)    ! Global neighbour index
               rproc = qproc(iqq) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  call indv_mpi(iqq,iloc,jloc,nloc)
                  is(iq) = indp(iloc,jloc,nloc)
               end if
               iqq = iss_g(iqg)    ! Global neighbour index
               rproc = qproc(iqq) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  call indv_mpi(iqq,iloc,jloc,nloc)
                  iss(iq) = indp(iloc,jloc,nloc)
               end if

               iqq = ie_g(iqg)    ! Global neighbour index
               rproc = qproc(iqq) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  call indv_mpi(iqq,iloc,jloc,nloc)
                  ie(iq) = indp(iloc,jloc,nloc)
               end if
               iqq = iee_g(iqg)    ! Global neighbour index
               rproc = qproc(iqq) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  call indv_mpi(iqq,iloc,jloc,nloc)
                  iee(iq) = indp(iloc,jloc,nloc)
               end if

               iqq = iw_g(iqg)    ! Global neighbour index
               rproc = qproc(iqq) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  call indv_mpi(iqq,iloc,jloc,nloc)
                  iw(iq) = indp(iloc,jloc,nloc)
               end if
               iqq = iww_g(iqg)    ! Global neighbour index
               rproc = qproc(iqq) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  call indv_mpi(iqq,iloc,jloc,nloc)
                  iww(iq) = indp(iloc,jloc,nloc)
               end if

               ! Note that the model only needs a limited set of the diagonal
               ! index arrays
               iqq = ine_g(iqg)    ! Global neighbour index
               rproc = qproc(iqq) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  call indv_mpi(iqq,iloc,jloc,nloc)
                  ine(iq) = indp(iloc,jloc,nloc)
               end if

               iqq = ise_g(iqg)    ! Global neighbour index
               rproc = qproc(iqq) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  call indv_mpi(iqq,iloc,jloc,nloc)
                  ise(iq) = indp(iloc,jloc,nloc)
               end if

               iqq = ien_g(iqg)    ! Global neighbour index
               rproc = qproc(iqq) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  call indv_mpi(iqq,iloc,jloc,nloc)
                  ien(iq) = indp(iloc,jloc,nloc)
               end if

               iqq = iwn_g(iqg)    ! Global neighbour index
               rproc = qproc(iqq) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  call indv_mpi(iqq,iloc,jloc,nloc)
                  iwn(iq) = indp(iloc,jloc,nloc)
               end if

               iqq = inw_g(iqg)    ! Global neighbour index
               rproc = qproc(iqq) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  call indv_mpi(iqq,iloc,jloc,nloc)
                  inw(iq) = indp(iloc,jloc,nloc)
               end if

               iqq = isw_g(iqg)    ! Global neighbour index
               rproc = qproc(iqq) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  call indv_mpi(iqq,iloc,jloc,nloc)
                  isw(iq) = indp(iloc,jloc,nloc)
               end if

               iqq = ies_g(iqg)    ! Global neighbour index
               rproc = qproc(iqq) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  call indv_mpi(iqq,iloc,jloc,nloc)
                  ies(iq) = indp(iloc,jloc,nloc)
               end if

               iqq = iws_g(iqg)    ! Global neighbour index
               rproc = qproc(iqq) ! Processor that has this point
               if ( rproc == myid ) then ! Just copy the value
                  ! Convert global iqq to local value
                  call indv_mpi(iqq,iloc,jloc,nloc)
                  iws(iq) = indp(iloc,jloc,nloc)
               end if

            end do
         end do
      end do

      ! Correct within the same face only (not necessarily the same
      ! processor, but will be corrected later).
      ieu = ie
      iwu = iw
      inv = in
      isv = is
      iev = ie
      iwv = iw
      inu = in
      isu = is

      ! Initialise the edge variables
      edge_w(:) = ioff(:) == 0
      edge_s(:) = joff(:) == 0
      edge_n(:) = joff(:) == il_g - jpan
      edge_e(:) = ioff(:) == il_g - ipan

      nr = 1

      bnds(:)%len = 0
      bnds(:)%rlenh = 0
      bnds(:)%slenh = 0
      bnds(:)%rlen = 0
      bnds(:)%slen = 0
      bnds(:)%rlenx = 0
      bnds(:)%slenx = 0
      bnds(:)%rlen2 = 0
      bnds(:)%slen2 = 0
      bnds(:)%rlen_uv = 0
      bnds(:)%slen_uv = 0
      bnds(:)%rlenx_uv = 0
      bnds(:)%slenx_uv = 0
      bnds(:)%rlen2_uv = 0
      bnds(:)%slen2_uv = 0

!     In the first pass through, set up list of points to be requested from
!     other processors. These points are placed in the "iextra" region at the
!     end of arrays. The indirect indices are updated to point into this 
!     region.
      iext = 0
 
      rcolsp(:)%ihbg(1) = 1
      rcolsp(:)%ihfn(1) = 0
      
      do icol=1,maxcolour
      
         do n=1,npan

            !     Start with N edge
            j=jpan
            do i=1,ipan
               iq = indg(i,j,n)
               iqq = in_g(iq)
               ! Which processor has this point
               rproc = qproc(iqq)
               if ( rproc == myid ) cycle ! Don't add points already on this proc.
               iql = indp(i,j,n)  !  Local index
               if ( colourmask(iql) /= icol ) cycle
               ! Add this point to request list
               call check_bnds_alloc(rproc, iext)
               bnds(rproc)%rlenh = bnds(rproc)%rlenh + 1
               rcolsp(rproc)%ihfn(icol) = rcolsp(rproc)%ihfn(icol) + 1
               bnds(rproc)%request_list(bnds(rproc)%rlenh) = iqq
               ! Increment extended region index
               iext = iext + 1
               bnds(rproc)%unpack_list(bnds(rproc)%rlenh) = iext
               in(iql) = ifull+iext
            end do

            !     E edge
            i = ipan
            do j=1,jpan
               iq = indg(i,j,n)
               iqq = ie_g(iq)
               ! Which processor has this point
               rproc = qproc(iqq)
               if ( rproc == myid ) cycle ! Don't add points already on this proc.
               iql = indp(i,j,n)  !  Local index
               if ( colourmask(iql) /= icol ) cycle
               ! Add this point to request list
               call check_bnds_alloc(rproc, iext)
               bnds(rproc)%rlenh = bnds(rproc)%rlenh + 1
               rcolsp(rproc)%ihfn(icol) = rcolsp(rproc)%ihfn(icol) + 1
               bnds(rproc)%request_list(bnds(rproc)%rlenh) = iqq
               ! Increment extended region index
               iext = iext + 1
               bnds(rproc)%unpack_list(bnds(rproc)%rlenh) = iext
               ie(iql) = ifull+iext
            end do
         end do ! n=1,npan
         
         if ( icol < 3 ) then
            rcolsp(:)%ihbg(icol+1) = rcolsp(:)%ihfn(icol) + 1
            rcolsp(:)%ihfn(icol+1) = rcolsp(:)%ihfn(icol)
         end if
      
      end do ! icol=1,maxcolour
      
      bnds(:)%rlen = bnds(:)%rlenh  ! so that they are appended
      rcolsp(:)%ifbg(1) = rcolsp(:)%ihfn(3) + 1
      rcolsp(:)%iffn(1) = rcolsp(:)%ihfn(3)
      
      do icol=1,maxcolour
      
         do n=1,npan

            !     W edge
            i = 1
            do j=1,jpan
               ! 1D code takes care of corners separately at the end so only goes
               ! over 1,jpan here.
               iq = indg(i,j,n)
               iqq = iw_g(iq) ! iqq is the global index of the required neighbouring point.
               ! Which processor has this point
               rproc = qproc(iqq)
               if ( rproc == myid ) cycle ! Don't add points already on this proc.
               iql = indp(i,j,n)  !  Local index
               if ( colourmask(iql) /= icol ) cycle
               ! Add this point to request list
               call check_bnds_alloc(rproc, iext)
               bnds(rproc)%rlen = bnds(rproc)%rlen + 1
               rcolsp(rproc)%iffn(icol) = rcolsp(rproc)%iffn(icol) + 1
               bnds(rproc)%request_list(bnds(rproc)%rlen) = iqq
               ! Increment extended region index
               iext = iext + 1
               bnds(rproc)%unpack_list(bnds(rproc)%rlen) = iext
               iw(iql) = ifull+iext
            end do

            !     S edge
            j=1
            do i=1,ipan
               iq = indg(i,j,n)
               iqq = is_g(iq)
               ! Which processor has this point
               rproc = qproc(iqq)
               if ( rproc == myid ) cycle ! Don't add points already on this proc.
               iql = indp(i,j,n)  !  Local index
               if ( colourmask(iql) /= icol ) cycle
               ! Add this point to request list
               call check_bnds_alloc(rproc, iext)
               bnds(rproc)%rlen = bnds(rproc)%rlen + 1
               rcolsp(rproc)%iffn(icol) = rcolsp(rproc)%iffn(icol) + 1
               bnds(rproc)%request_list(bnds(rproc)%rlen) = iqq
               ! Increment extended region index
               iext = iext + 1
               bnds(rproc)%unpack_list(bnds(rproc)%rlen) = iext
               is(iql) = ifull+iext
            end do
         end do ! n=1,npan

         if ( icol < 3 ) then
            rcolsp(:)%ifbg(icol+1) = rcolsp(:)%iffn(icol) + 1
            rcolsp(:)%iffn(icol+1) = rcolsp(:)%iffn(icol)
         end if
         
      end do ! icol=1,maxcolour

      bnds(:)%rlenx = bnds(:)%rlen  ! so that they're appended.
      
!     Now handle the special corner values that need to be remapped
!     This adds to rlen, so needs to come before the _XX stuff.
      do n=1,npan
         ! NE, EN
         iq = indp(ipan,jpan,n)
         iqg = indg(ipan,jpan,n)
!        This might matter for leen etc but not here
!         if ( edge_e .and. edge_n ) then
!            ! This is a real vertex
         iqq = ine_g(iqg)
         ! Which processor has this point
         rproc = qproc(iqq)
         if ( rproc /= myid ) then ! Add to list
            call check_bnds_alloc(rproc, iext)
            bnds(rproc)%rlenx = bnds(rproc)%rlenx + 1
            bnds(rproc)%request_list(bnds(rproc)%rlenx) = iqq
            ! Increment extended region index
            iext = iext + 1
            bnds(rproc)%unpack_list(bnds(rproc)%rlenx) = iext
            ine(iq) = ifull+iext
         end if

         if (ien_g(iqg)==ine_g(iqg)) then
            ien(iq)=ine(iq)
         else
            iqq = ien_g(iqg)
            ! Which processor has this point
            rproc = qproc(iqq)
            if ( rproc /= myid ) then ! Add to list
               call check_bnds_alloc(rproc, iext)
               bnds(rproc)%rlenx = bnds(rproc)%rlenx + 1
               bnds(rproc)%request_list(bnds(rproc)%rlenx) = iqq
               ! Increment extended region index
               iext = iext + 1
               bnds(rproc)%unpack_list(bnds(rproc)%rlenx) = iext
               ien(iq) = ifull+iext
            end if
         end if

         iq = indp(ipan,1,n)
         iqg = indg(ipan,1,n)
         iqq = ise_g(iqg)
         ! Which processor has this point
         rproc = qproc(iqq)
         if ( rproc /= myid ) then ! Add to list
            call check_bnds_alloc(rproc, iext)
            bnds(rproc)%rlenx = bnds(rproc)%rlenx + 1
            bnds(rproc)%request_list(bnds(rproc)%rlenx) = iqq
            ! Increment extended region index
            iext = iext + 1
            bnds(rproc)%unpack_list(bnds(rproc)%rlenx) = iext
            ise(iq) = ifull+iext
         end if

         if (ies_g(iqg)==ise_g(iqg)) then
            ies(iq)=ise(iq)
         else
            iqq = ies_g(iqg)
            ! Which processor has this point
            rproc = qproc(iqq)
            if ( rproc /= myid ) then ! Add to list
               call check_bnds_alloc(rproc, iext)
               bnds(rproc)%rlenx = bnds(rproc)%rlenx + 1
               bnds(rproc)%request_list(bnds(rproc)%rlenx) = iqq
               ! Increment extended region index
               iext = iext + 1
               bnds(rproc)%unpack_list(bnds(rproc)%rlenx) = iext
               ies(iq) = ifull+iext
            end if
         end if

         iq = indp(1,jpan,n)
         iqg = indg(1,jpan,n)
         iqq = iwn_g(iqg)
         ! Which processor has this point
         rproc = qproc(iqq)
         if ( rproc /= myid ) then ! Add to list
            call check_bnds_alloc(rproc, iext)
            bnds(rproc)%rlenx = bnds(rproc)%rlenx + 1
            bnds(rproc)%request_list(bnds(rproc)%rlenx) = iqq
            ! Increment extended region index
            iext = iext + 1
            bnds(rproc)%unpack_list(bnds(rproc)%rlenx) = iext
            iwn(iq) = ifull+iext
         end if

         if (inw_g(iqg)==iwn_g(iqg)) then
            inw(iq)=iwn(iq)
         else
            iqq = inw_g(iqg)
            ! Which processor has this point
            rproc = qproc(iqq)
            if ( rproc /= myid ) then ! Add to list
               call check_bnds_alloc(rproc, iext)
               bnds(rproc)%rlenx = bnds(rproc)%rlenx + 1
               bnds(rproc)%request_list(bnds(rproc)%rlenx) = iqq
               ! Increment extended region index
               iext = iext + 1
               bnds(rproc)%unpack_list(bnds(rproc)%rlenx) = iext
               inw(iq) = ifull+iext
            end if
         end if

         iq = indp(1,1,n)
         iqg = indg(1,1,n)
         iqq = isw_g(iqg)
         ! Which processor has this point
         rproc = qproc(iqq)
         if ( rproc /= myid ) then ! Add to list
            call check_bnds_alloc(rproc, iext)
            bnds(rproc)%rlenx = bnds(rproc)%rlenx + 1
            bnds(rproc)%request_list(bnds(rproc)%rlenx) = iqq
            ! Increment extended region index
            iext = iext + 1
            bnds(rproc)%unpack_list(bnds(rproc)%rlenx) = iext
            isw(iq) = ifull+iext
         end if

         if (iws_g(iqg)==isw_g(iqg)) then
            iws(iq)=isw(iq)
         else
            iqq = iws_g(iqg)
            ! Which processor has this point
            rproc = qproc(iqq)
            if ( rproc /= myid ) then ! Add to list
               call check_bnds_alloc(rproc, iext)
               bnds(rproc)%rlenx = bnds(rproc)%rlenx + 1
               bnds(rproc)%request_list(bnds(rproc)%rlenx) = iqq
               ! Increment extended region index
               iext = iext + 1
               bnds(rproc)%unpack_list(bnds(rproc)%rlenx) = iext
               iws(iq) = ifull+iext
            end if
         end if

      end do

!     Now set up the second row
      bnds(:)%rlen2 = bnds(:)%rlenx  ! so that they're appended.
      do n=1,npan

         !     Start with W edge
         i = 1
         do j=1,jpan
            iq = indg(i,j,n)
            iqq = iww_g(iq)

            ! Which processor has this point
            rproc = qproc(iqq)
            if ( rproc == myid ) cycle ! Don't add points already on this proc.
            ! Add this point to request list
            call check_bnds_alloc(rproc, iext)
            bnds(rproc)%rlen2 = bnds(rproc)%rlen2 + 1
            bnds(rproc)%request_list(bnds(rproc)%rlen2) = iqq
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
            iqq = inn_g(iq)

            ! Which processor has this point
            rproc = qproc(iqq)
            if ( rproc == myid ) cycle ! Don't add points already on this proc.
            ! Add this point to request list
            call check_bnds_alloc(rproc, iext)
            bnds(rproc)%rlen2 = bnds(rproc)%rlen2 + 1
            bnds(rproc)%request_list(bnds(rproc)%rlen2) = iqq
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
            iqq = iee_g(iq)

            ! Which processor has this point
            rproc = qproc(iqq)
            if ( rproc == myid ) cycle ! Don't add points already on this proc.
            ! Add this point to request list
            call check_bnds_alloc(rproc, iext)
            bnds(rproc)%rlen2 = bnds(rproc)%rlen2 + 1
            bnds(rproc)%request_list(bnds(rproc)%rlen2) = iqq
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
            iqq = iss_g(iq)

            ! Which processor has this point
            rproc = qproc(iqq)
            if ( rproc == myid ) cycle ! Don't add points already on this proc.
            ! Add this point to request list
            call check_bnds_alloc(rproc, iext)
            bnds(rproc)%rlen2 = bnds(rproc)%rlen2 + 1
            bnds(rproc)%request_list(bnds(rproc)%rlen2) = iqq
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
         if ( edge_n(n-noff) .and. edge_e(n-noff) ) then
            call fix_index2(lnee_g(n-noff),lnee,n,bnds,iext)
            call fix_index2(leen_g(n-noff),leen,n,bnds,iext)
            call fix_index2(lenn_g(n-noff),lenn,n,bnds,iext)
            call fix_index2(lnne_g(n-noff),lnne,n,bnds,iext)
         else if ( edge_e(n-noff) ) then
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
         if ( edge_n(n-noff) .and. edge_w(n-noff) ) then
            call fix_index2(lnww_g(n-noff),lnww,n,bnds,iext)
            call fix_index2(lwwn_g(n-noff),lwwn,n,bnds,iext)
            call fix_index2(lwnn_g(n-noff),lwnn,n,bnds,iext)
            call fix_index2(lnnw_g(n-noff),lnnw,n,bnds,iext)
         else if ( edge_w(n-noff) ) then
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
         if ( edge_s(n-noff) .and. edge_e(n-noff) ) then
            call fix_index2(lsee_g(n-noff),lsee,n,bnds,iext)
            call fix_index2(lees_g(n-noff),lees,n,bnds,iext)
            call fix_index2(less_g(n-noff),less,n,bnds,iext)
            call fix_index2(lsse_g(n-noff),lsse,n,bnds,iext)
         else if ( edge_e(n-noff) ) then
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
         if ( edge_s(n-noff) .and. edge_w(n-noff) ) then
            call fix_index2(lsww_g(n-noff),lsww,n,bnds,iext)
            call fix_index2(lwws_g(n-noff),lwws,n,bnds,iext)
            call fix_index2(lwss_g(n-noff),lwss,n,bnds,iext)
            call fix_index2(lssw_g(n-noff),lssw,n,bnds,iext)
         else if ( edge_w(n-noff) ) then
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
         call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
      end if



!     Now, for each processor send the list of points I want.
!     The state of being a neighbour is reflexive so only expect to
!     recv from those processors I send to (are there grid arrangements for
!     which this would not be true?)
!     Get the complete request lists by using rlen2
      nreq = 0
      do iproc = 1,nproc-1  !
         sproc = modulo(myid+iproc,nproc)  ! Send to
         if (bnds(sproc)%rlen2 > 0 ) then
            nreq = nreq + 1
            ! Use the maximum size in the recv call.
            call MPI_IRecv( bnds(sproc)%send_list(1), bnds(sproc)%len, &
                 MPI_INTEGER, sproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do
      do iproc = 1,nproc-1  !
         sproc = modulo(myid+iproc,nproc)  ! Send to
         if (bnds(sproc)%rlen2 > 0 ) then
            ! Send list of requests
            nreq = nreq + 1
            call MPI_ISend( bnds(sproc)%request_list(1), bnds(sproc)%rlen2, &
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
         if (bnds(sproc)%rlen2 > 0 ) then
            ! First half of nreq are recv
            nreq = nreq + 1
            call MPI_Get_count(status(1,nreq), MPI_INTEGER, count, ierr)
            ! This the number of points I have to send to rproc.
            bnds(sproc)%slen2 = count
         end if
      end do

!     For rlen and rlenx, just communicate the lengths. The indices have 
!     already been taken care of.
      scolsp(:)%ihfn(1) = 0
      scolsp(:)%ihfn(2) = 0
      scolsp(:)%ihfn(3) = 0
      scolsp(:)%iffn(1) = 0
      scolsp(:)%iffn(2) = 0
      scolsp(:)%iffn(3) = 0
      nreq = 0
      do iproc = 1,nproc-1  !
         ! Send and recv from same proc
         sproc = modulo(myid+iproc,nproc)  ! Send to
         if (bnds(sproc)%rlen2 > 0 ) then
            nreq = nreq + 1
            call MPI_IRecv( dumrb(:,sproc), 9, &
                 MPI_INTEGER, sproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do
      do iproc = 1,nproc-1  !
         ! Send and recv from same proc
         sproc = modulo(myid+iproc,nproc)  ! Send to
         if (bnds(sproc)%rlen2 > 0 ) then
            nreq = nreq + 1
            dumsb(1,sproc) = bnds(sproc)%rlenh
            dumsb(2,sproc) = bnds(sproc)%rlen
            dumsb(3,sproc) = bnds(sproc)%rlenx
            dumsb(4,sproc) = rcolsp(sproc)%ihfn(1)
            dumsb(5,sproc) = rcolsp(sproc)%ihfn(2)
            dumsb(6,sproc) = rcolsp(sproc)%ihfn(3)
            dumsb(7,sproc) = rcolsp(sproc)%iffn(1)
            dumsb(8,sproc) = rcolsp(sproc)%iffn(2)
            dumsb(9,sproc) = rcolsp(sproc)%iffn(3)
            call MPI_ISend( dumsb(:,sproc), 9, &
                 MPI_INTEGER, sproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do
      if ( nreq > 0 ) then
         call MPI_Waitall(nreq,ireq,status,ierr)
      end if
      do iproc = 1,nproc-1
         sproc = modulo(myid+iproc,nproc)
         if ( bnds(sproc)%rlen2 > 0 ) then
            bnds(sproc)%slenh     = dumrb(1,sproc)
            bnds(sproc)%slen      = dumrb(2,sproc)
            bnds(sproc)%slenx     = dumrb(3,sproc)
            scolsp(sproc)%ihfn(1) = dumrb(4,sproc)
            scolsp(sproc)%ihfn(2) = dumrb(5,sproc)
            scolsp(sproc)%ihfn(3) = dumrb(6,sproc)
            scolsp(sproc)%iffn(1) = dumrb(7,sproc)
            scolsp(sproc)%iffn(2) = dumrb(8,sproc)
            scolsp(sproc)%iffn(3) = dumrb(9,sproc)
         end if
      end do
      scolsp(:)%ihbg(1) = 1
      scolsp(:)%ihbg(2) = scolsp(:)%ihfn(1) + 1
      scolsp(:)%ihbg(3) = scolsp(:)%ihfn(2) + 1
      scolsp(:)%ifbg(1) = scolsp(:)%ihfn(3) + 1
      scolsp(:)%ifbg(2) = scolsp(:)%iffn(1) + 1
      scolsp(:)%ifbg(3) = scolsp(:)%iffn(2) + 1

!     Start of UV section

!     In the first pass through, set up list of points to be requested from
!     other processors. In the 1D code values on the same processor are
!     copied only if they have to be swapped.
!     This only makes a difference on 1, 2 or 3 processors.

      iextu = 0
      iextv = 0

      ! save start of isv indices
      rsplit(:)%isvbg = 1
      rsplit(:)%iwufn = 0

      !     S edge, V
      j=1
      do n=1,npan
         do i=1,ipan
            iqg = indg(i,j,n)
            iqq = is_g(iqg)
            rproc = qproc(iqq)
            swap = edge_s(n-noff) .and. swap_s(n-noff)
            if ( rproc == myid .and. .not. swap ) cycle
            call check_bnds_alloc(rproc, iextv)
            bnds(rproc)%rlen_uv = bnds(rproc)%rlen_uv + 1
            bnds(rproc)%request_list_uv(bnds(rproc)%rlen_uv) = -iqq
            rsplit(rproc)%iwufn = rsplit(rproc)%iwufn + 1
            ! Increment extended region index
            iextv = iextv + 1
            bnds(rproc)%unpack_list_uv(bnds(rproc)%rlen_uv) = -iextv
            iql = indp(i,j,n)  !  Local index
            isv(iql) = ifull+iextv
            bnds(rproc)%uv_swap(bnds(rproc)%rlen_uv) = swap
            bnds(rproc)%uv_neg(bnds(rproc)%rlen_uv) = .false.
         end do
      end do
         
      !     Start with W edge, U values
      i = 1
      do n=1,npan
         do j=1,jpan
            iqg = indg(i,j,n)
            iqq = iw_g(iqg)
            ! Which processor has this point
            rproc = qproc(iqq)
            ! Only need to add to bounds region if it's on another processor
            ! or if it's on this processor and needs to be swapped.
            ! Decide if u/v need to be swapped. My face is n-noff
            swap = edge_w(n-noff) .and. swap_w(n-noff)
            if ( rproc == myid .and. .not. swap ) cycle
            call check_bnds_alloc(rproc, iextu)
            bnds(rproc)%rlen_uv = bnds(rproc)%rlen_uv + 1
            bnds(rproc)%request_list_uv(bnds(rproc)%rlen_uv) = iqq
            rsplit(rproc)%iwufn = rsplit(rproc)%iwufn + 1
            ! Increment extended region index
            iextu = iextu + 1
            bnds(rproc)%unpack_list_uv(bnds(rproc)%rlen_uv) = iextu
            iql = indp(i,j,n)  !  Local index
            iwu(iql) = ifull+iextu
            bnds(rproc)%uv_swap(bnds(rproc)%rlen_uv) = swap
            bnds(rproc)%uv_neg(bnds(rproc)%rlen_uv) = .false.
         end do
      end do

      ! save start of inv indices
      rsplit(:)%invbg = rsplit(:)%iwufn + 1
      rsplit(:)%ieufn = rsplit(:)%iwufn

      !     N edge (V)
      j=jpan
      do n=1,npan
         do i=1,ipan
            iqg = indg(i,j,n)
            iqq = in_g(iqg)
            rproc = qproc(iqq)
            swap = edge_n(n-noff) .and. swap_n(n-noff)
            if ( rproc == myid .and. .not. swap ) cycle
            ! Add this point to request list
            call check_bnds_alloc(rproc, iextv)
            bnds(rproc)%rlen_uv = bnds(rproc)%rlen_uv + 1
            ! to show that this is v rather than u, flip sign
            bnds(rproc)%request_list_uv(bnds(rproc)%rlen_uv) = -iqq
            rsplit(rproc)%ieufn = rsplit(rproc)%ieufn + 1
            ! Increment extended region index
            iextv = iextv + 1
            bnds(rproc)%unpack_list_uv(bnds(rproc)%rlen_uv) = -iextv
            iql = indp(i,j,n)  !  Local index
            inv(iql) = ifull+iextv
            bnds(rproc)%uv_swap(bnds(rproc)%rlen_uv) = swap
            bnds(rproc)%uv_neg(bnds(rproc)%rlen_uv) = .false.
         end do
      end do

      !     E edge, U
      i = ipan
      do n=1,npan
         do j=1,jpan
            iqg = indg(i,j,n)
            iqq = ie_g(iqg)
            rproc = qproc(iqq)
            swap = edge_e(n-noff) .and. swap_e(n-noff)
            if ( rproc == myid .and. .not. swap ) cycle
            ! Add this point to request list
            call check_bnds_alloc(rproc, iextu)
            bnds(rproc)%rlen_uv = bnds(rproc)%rlen_uv + 1
            bnds(rproc)%request_list_uv(bnds(rproc)%rlen_uv) = iqq
            rsplit(rproc)%ieufn = rsplit(rproc)%ieufn + 1
            ! Increment extended region index
            iextu = iextu + 1
            bnds(rproc)%unpack_list_uv(bnds(rproc)%rlen_uv) = iextu
            iql = indp(i,j,n)  !  Local index
            ieu(iql) = ifull+iextu
            bnds(rproc)%uv_swap(bnds(rproc)%rlen_uv) = swap
            bnds(rproc)%uv_neg(bnds(rproc)%rlen_uv) = .false.
         end do
      end do

!     Second pass
      bnds(:)%rlen2_uv = bnds(:)%rlen_uv
      bnds(:)%slen2_uv = bnds(:)%slen_uv
      ieeu = iee
      iwwu = iww
      innv = inn
      issv = iss

      ! save start of issv indices
      rsplit(:)%issvbg = rsplit(:)%ieufn + 1
      rsplit(:)%iwwufn = rsplit(:)%ieufn

      !     S edge, V
      j=1
      do n=1,npan
         do i=1,ipan
            iqg = indg(i,j,n)
            iqq = iss_g(iqg)
            rproc = qproc(iqq)
            swap = edge_s(n-noff) .and. swap_s(n-noff)
            if ( rproc == myid .and. .not. swap ) cycle
            call check_bnds_alloc(rproc, iextv)
            bnds(rproc)%rlen2_uv = bnds(rproc)%rlen2_uv + 1
            bnds(rproc)%request_list_uv(bnds(rproc)%rlen2_uv) = -iqq
            rsplit(rproc)%iwwufn = rsplit(rproc)%iwwufn + 1
            ! Increment extended region index
            iextv = iextv + 1
            bnds(rproc)%unpack_list_uv(bnds(rproc)%rlen2_uv) = -iextv
            iql = indp(i,j,n)  !  Local index
            issv(iql) = ifull+iextv
            bnds(rproc)%uv_swap(bnds(rproc)%rlen2_uv) = swap
            bnds(rproc)%uv_neg(bnds(rproc)%rlen2_uv) = .false.
         end do
      end do ! n=1,npan

      !     Start with W edge, U values
      i = 1
      do n=1,npan
         do j=1,jpan
            iqg = indg(i,j,n)
            iqq = iww_g(iqg)
            ! Which processor has this point
            rproc = qproc(iqq)
            swap = edge_w(n-noff) .and. swap_w(n-noff)
            if ( rproc == myid .and. .not. swap ) cycle
            ! Add this point to request list
            call check_bnds_alloc(rproc, iextu)
            bnds(rproc)%rlen2_uv = bnds(rproc)%rlen2_uv + 1
            bnds(rproc)%request_list_uv(bnds(rproc)%rlen2_uv) = iqq
            rsplit(rproc)%iwwufn = rsplit(rproc)%iwwufn + 1
            ! Increment extended region index
            iextu = iextu + 1
            bnds(rproc)%unpack_list_uv(bnds(rproc)%rlen2_uv) = iextu
            iql = indp(i,j,n)  !  Local index
            iwwu(iql) = ifull+iextu
            ! Decide if u/v need to be swapped. My face is n-noff
            bnds(rproc)%uv_swap(bnds(rproc)%rlen2_uv) = swap
            bnds(rproc)%uv_neg(bnds(rproc)%rlen2_uv) = .false.
         end do
      end do

      ! save start of innv indices
      rsplit(:)%innvbg = rsplit(:)%iwwufn + 1
      rsplit(:)%ieeufn = rsplit(:)%iwwufn

      !     N edge (V)
      j=jpan
      do n=1,npan
         do i=1,ipan
            iqg = indg(i,j,n)
            iqq = inn_g(iqg)
            rproc = qproc(iqq)
            swap = edge_n(n-noff) .and. swap_n(n-noff)
            if ( rproc == myid .and. .not. swap ) cycle
            ! Add this point to request list
            call check_bnds_alloc(rproc, iextv)
            bnds(rproc)%rlen2_uv = bnds(rproc)%rlen2_uv + 1
            ! to show that this is v rather than u, flip sign
            bnds(rproc)%request_list_uv(bnds(rproc)%rlen2_uv) = -iqq
            rsplit(rproc)%ieeufn = rsplit(rproc)%ieeufn + 1
            ! Increment extended region index
            iextv = iextv + 1
            bnds(rproc)%unpack_list_uv(bnds(rproc)%rlen2_uv) = -iextv
            iql = indp(i,j,n)  !  Local index
            innv(iql) = ifull+iextv
            bnds(rproc)%uv_swap(bnds(rproc)%rlen2_uv) = swap
            bnds(rproc)%uv_neg(bnds(rproc)%rlen2_uv) = .false.
         end do
      end do

      !     E edge, U
      i = ipan
      do n=1,npan
         do j=1,jpan
            iqg = indg(i,j,n)
            iqq = iee_g(iqg)
            rproc = qproc(iqq)
            swap = edge_e(n-noff) .and. swap_e(n-noff)
            if ( rproc == myid .and. .not. swap ) cycle
            ! Add this point to request list
            call check_bnds_alloc(rproc, iextu)
            bnds(rproc)%rlen2_uv = bnds(rproc)%rlen2_uv + 1
            bnds(rproc)%request_list_uv(bnds(rproc)%rlen2_uv) = iqq
            rsplit(rproc)%ieeufn = rsplit(rproc)%ieeufn + 1
            ! Increment extended region index
            iextu = iextu + 1
            bnds(rproc)%unpack_list_uv(bnds(rproc)%rlen2_uv) = iextu
            iql = indp(i,j,n)  !  Local index
            ieeu(iql) = ifull+iextu
            bnds(rproc)%uv_swap(bnds(rproc)%rlen2_uv) = swap
            bnds(rproc)%uv_neg(bnds(rproc)%rlen2_uv) = .false.
         end do
      end do

      ! Third pass
      bnds(:)%rlenx_uv = bnds(:)%rlen2_uv
      bnds(:)%slenx_uv = bnds(:)%slen2_uv

      ! save start of isu indices
      rsplit(:)%isubg = rsplit(:)%ieeufn + 1
      rsplit(:)%iwvfn = rsplit(:)%ieeufn

      !     S edge, U
      j=1
      do n=1,npan
         do i=1,ipan
            iqg = indg(i,j,n)
            iqq = is_g(iqg)
            rproc = qproc(iqq)
            swap = edge_s(n-noff) .and. swap_s(n-noff)
            if ( rproc == myid .and. .not. swap ) cycle
            call check_bnds_alloc(rproc, iextu)
            bnds(rproc)%rlenx_uv = bnds(rproc)%rlenx_uv + 1
            bnds(rproc)%request_list_uv(bnds(rproc)%rlenx_uv) = iqq
            rsplit(rproc)%iwvfn = rsplit(rproc)%iwvfn + 1
            ! Increment extended region index
            iextu = iextu + 1
            bnds(rproc)%unpack_list_uv(bnds(rproc)%rlenx_uv) = iextu
            iql = indp(i,j,n)  !  Local index
            isu(iql) = ifull+iextu
            bnds(rproc)%uv_swap(bnds(rproc)%rlenx_uv) = swap
            bnds(rproc)%uv_neg(bnds(rproc)%rlenx_uv) = swap
         end do
      end do ! n=1,npan

      !     Start with W edge, V values
      i = 1
      do n=1,npan
         do j=1,jpan
            iqg = indg(i,j,n)
            iqq = iw_g(iqg)
            ! Which processor has this point
            rproc = qproc(iqq)
            ! Only need to add to bounds region if it's on another processor
            ! or if it's on this processor and needs to be swapped.
            ! Decide if u/v need to be swapped. My face is n-noff
            swap = edge_w(n-noff) .and. swap_w(n-noff)
            if ( rproc == myid .and. .not. swap ) cycle
            call check_bnds_alloc(rproc, iextv)
            bnds(rproc)%rlenx_uv = bnds(rproc)%rlenx_uv + 1
            ! to show that this is v rather than u, flip sign
            bnds(rproc)%request_list_uv(bnds(rproc)%rlenx_uv) = -iqq
            rsplit(rproc)%iwvfn = rsplit(rproc)%iwvfn + 1
            ! Increment extended region index
            iextv = iextv + 1
            bnds(rproc)%unpack_list_uv(bnds(rproc)%rlenx_uv) = -iextv
            iql = indp(i,j,n)  !  Local index
            iwv(iql) = ifull+iextv
            bnds(rproc)%uv_swap(bnds(rproc)%rlenx_uv) = swap
            bnds(rproc)%uv_neg(bnds(rproc)%rlenx_uv) = swap
         end do
      end do

      ! start save of inu indices
      rsplit(:)%inubg = rsplit(:)%iwvfn + 1
      rsplit(:)%ievfn = rsplit(:)%iwvfn

      !     N edge (U)
      j=jpan
      do n=1,npan
         do i=1,ipan
            iqg = indg(i,j,n)
            iqq = in_g(iqg)
            rproc = qproc(iqq)
            swap = edge_n(n-noff) .and. swap_n(n-noff)
            if ( rproc == myid .and. .not. swap ) cycle
            ! Add this point to request list
            call check_bnds_alloc(rproc, iextu)
            bnds(rproc)%rlenx_uv = bnds(rproc)%rlenx_uv + 1
            bnds(rproc)%request_list_uv(bnds(rproc)%rlenx_uv) = iqq
            rsplit(rproc)%ievfn = rsplit(rproc)%ievfn + 1
            ! Increment extended region index
            iextu = iextu + 1
            bnds(rproc)%unpack_list_uv(bnds(rproc)%rlenx_uv) = iextu
            iql = indp(i,j,n)  !  Local index
            inu(iql) = ifull+iextu
            bnds(rproc)%uv_swap(bnds(rproc)%rlenx_uv) = swap
            bnds(rproc)%uv_neg(bnds(rproc)%rlenx_uv) = swap
         end do
      end do

      !     E edge, V
      i = ipan
      do n=1,npan
         do j=1,jpan
            iqg = indg(i,j,n)
            iqq = ie_g(iqg)
            rproc = qproc(iqq)
            swap = edge_e(n-noff) .and. swap_e(n-noff)
            if ( rproc == myid .and. .not. swap ) cycle
            ! Add this point to request list
            call check_bnds_alloc(rproc, iextv)
            bnds(rproc)%rlenx_uv = bnds(rproc)%rlenx_uv + 1
            ! to show that this is v rather than u, flip sign
            bnds(rproc)%request_list_uv(bnds(rproc)%rlenx_uv) = -iqq
            rsplit(rproc)%ievfn = rsplit(rproc)%ievfn + 1
            ! Increment extended region index
            iextv = iextv + 1
            bnds(rproc)%unpack_list_uv(bnds(rproc)%rlenx_uv) = -iextv
            iql = indp(i,j,n)  !  Local index
            iev(iql) = ifull+iextv
            bnds(rproc)%uv_swap(bnds(rproc)%rlenx_uv) = swap
            bnds(rproc)%uv_neg(bnds(rproc)%rlenx_uv) = swap
         end do
      end do

      if ( iextu > iextra ) then
         write(6,*) "IEXTU too large", iextu, iextra
         call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
      end if

      if ( iextv > iextra ) then
         write(6,*) "IEXTV too large", iextv, iextra
         call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
      end if

!     Nearest neighbours are defined as those points which send/recv 
!     boundary information.
      neighbour = bnds(:)%rlen2 > 0

#ifdef DEBUG
      write(6,*) "Bounds", myid, "SLEN ", bnds(:)%slen
      write(6,*) "Bounds", myid, "RLEN ", bnds(:)%rlen
      write(6,*) "Bounds", myid, "SLENX", bnds(:)%slenx
      write(6,*) "Bounds", myid, "RLENX", bnds(:)%rlenx
      write(6,*) "Bounds", myid, "SLEN2", bnds(:)%slen2
      write(6,*) "Bounds", myid, "RLEN2", bnds(:)%rlen2
      write(6,*) "Neighbour", myid, neighbour
#endif

!     Now, for each processor send the list of points I want.
!     Send all possible pairs and don't assume anything about the symmetry.
!     For completeness, send zero length messages too.
!     Get the length from the message status
!     Also have to send the swap list

      nreq = 0
      do iproc = 1,nproc-1  !
         sproc = modulo(myid+iproc,nproc)
         if ( bnds(sproc)%rlenx_uv > 0 ) then
            nreq = nreq + 1
            ! Use the maximum size in the recv call.
            call MPI_IRecv(  bnds(sproc)%send_list_uv(1),  bnds(sproc)%len, &
                 MPI_INTEGER, sproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do
      do iproc = 1,nproc-1  !
         sproc = modulo(myid+iproc,nproc)
         if ( bnds(sproc)%rlenx_uv > 0 ) then
            ! Send list of requests
            nreq = nreq + 1
            call MPI_ISend( bnds(sproc)%request_list_uv(1), bnds(sproc)%rlenx_uv, &
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
         if ( bnds(sproc)%rlenx_uv > 0 ) then
            ! First half of nreq are recv
            nreq = nreq + 1
            call MPI_Get_count(status(1,nreq), MPI_INTEGER, count, ierr)
            ! This the number of points I have to send to rproc.
            bnds(sproc)%slenx_uv = count
         end if
      end do

!     For rlen_uv and rlenx_uv, etc, just communicate the lengths. The indices have 
!     already been taken care of.
      ssplit(:)%iwufn  = 0
      ssplit(:)%ieufn  = 0
      ssplit(:)%iwwufn = 0
      ssplit(:)%ieeufn = 0
      ssplit(:)%iwvfn  = 0
      ssplit(:)%ievfn  = 0
      nreq = 0
      do iproc = 1,nproc-1  !
         ! Send and recv from same proc
         sproc = modulo(myid+iproc,nproc)  ! Send to
         if ( bnds(sproc)%rlenx_uv > 0 ) then
            nreq = nreq + 1
            call MPI_IRecv( dumr(:,sproc), 8, &
                 MPI_INTEGER, sproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do
      do iproc = 1,nproc-1  !
         ! Send and recv from same proc
         sproc = modulo(myid+iproc,nproc)  ! Send to
         if ( bnds(sproc)%rlenx_uv > 0 ) then
            nreq = nreq + 1
            dums(1,sproc) = bnds(sproc)%rlen_uv
            dums(2,sproc) = bnds(sproc)%rlen2_uv
            dums(3,sproc) = rsplit(sproc)%iwufn
            dums(4,sproc) = rsplit(sproc)%ieufn
            dums(5,sproc) = rsplit(sproc)%iwwufn
            dums(6,sproc) = rsplit(sproc)%ieeufn
            dums(7,sproc) = rsplit(sproc)%iwvfn
            dums(8,sproc) = rsplit(sproc)%ievfn
            call MPI_ISend( dums(:,sproc), 8, &
                 MPI_INTEGER, sproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do
      if ( nreq > 0 ) then
         call MPI_Waitall(nreq,ireq,status,ierr)
      end if
      do iproc = 1,nproc-1
         sproc = modulo(myid+iproc,nproc)
         if ( bnds(sproc)%rlenx_uv > 0 ) then
            bnds(sproc)%slen_uv  = dumr(1,sproc)
            bnds(sproc)%slen2_uv = dumr(2,sproc)
            ssplit(sproc)%iwufn  = dumr(3,sproc)
            ssplit(sproc)%ieufn  = dumr(4,sproc)
            ssplit(sproc)%iwwufn = dumr(5,sproc)
            ssplit(sproc)%ieeufn = dumr(6,sproc)
            ssplit(sproc)%iwvfn  = dumr(7,sproc)
            ssplit(sproc)%ievfn  = dumr(8,sproc)
         end if
      end do
      ssplit(:)%isvbg  = 1
      ssplit(:)%invbg  = ssplit(:)%iwufn  + 1
      ssplit(:)%issvbg = ssplit(:)%ieufn  + 1
      ssplit(:)%innvbg = ssplit(:)%iwwufn + 1
      ssplit(:)%isubg  = ssplit(:)%ieeufn + 1
      ssplit(:)%inubg  = ssplit(:)%iwvfn  + 1

      ! Only send the swap list once
      nreq = 0
      do iproc = 1,nproc-1  !
         sproc = modulo(myid+iproc,nproc)
         if ( bnds(sproc)%rlenx_uv > 0 ) then
            nreq = nreq + 1
            ! Use the maximum size in the recv call.
            call MPI_IRecv( bnds(sproc)%send_swap(1), bnds(sproc)%slenx_uv, &
                 MPI_LOGICAL, sproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do
      do iproc = 1,nproc-1  !
         sproc = modulo(myid+iproc,nproc)
         if ( bnds(sproc)%rlenx_uv > 0 ) then
            ! Send list of requests
            nreq = nreq + 1
            call MPI_ISend( bnds(sproc)%uv_swap(1), bnds(sproc)%rlenx_uv, &
                 MPI_LOGICAL, sproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do
      if ( nreq > 0 ) then
         call MPI_Waitall(nreq,ireq,status,ierr)
      end if

      ! Only send the neg list once
      nreq = 0
      do iproc = 1,nproc-1  !
         sproc = modulo(myid+iproc,nproc)
         if ( bnds(sproc)%rlenx_uv > 0 ) then
            nreq = nreq + 1
            ! Use the maximum size in the recv call.
            call MPI_IRecv( bnds(sproc)%send_neg(1), bnds(sproc)%len, &
                 MPI_LOGICAL, sproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do
      do iproc = 1,nproc-1  !
         sproc = modulo(myid+iproc,nproc)
         if ( bnds(sproc)%rlenx_uv > 0 ) then
            ! Send list of requests
            nreq = nreq + 1
            call MPI_ISend( bnds(sproc)%uv_neg(1), bnds(sproc)%rlenx_uv, &
                 MPI_LOGICAL, sproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do
      if ( nreq > 0 ) then
         call MPI_Waitall(nreq,ireq,status,ierr)
      end if

!     Indices that are missed above (should be a better way to get these)
      do n=1,npan
         do j=1,jpan
            iwwu(indp(2,j,n)) = iwu(indp(1,j,n))
            ieeu(indp(ipan-1,j,n)) = ieu(indp(ipan,j,n))
         end do
         do i=1,ipan
            issv(indp(i,2,n)) = isv(indp(i,1,n))
            innv(indp(i,jpan-1,n)) = inv(indp(i,jpan,n))
         end do
      end do

#ifdef DEBUG
      write(6,*) "Bounds", myid, "SLEN_UV ", bnds(:)%slen_uv
      write(6,*) "Bounds", myid, "RLEN_UV ", bnds(:)%rlen_uv
      write(6,*) "Bounds", myid, "SLEN2_UV ", bnds(:)%slen2_uv
      write(6,*) "Bounds", myid, "RLEN2_UV ", bnds(:)%rlen2_uv
#endif

!  At the moment send_lists use global indices. Convert these to local.
      do iproc = 1,nproc-1  !
         sproc = modulo(myid+iproc,nproc)  ! Send to
         do iq=1,bnds(sproc)%slen2
            ! send_list(iq) is global point index, i, j, n are local
            call indv_mpi(bnds(sproc)%send_list(iq),i,j,n)
            bnds(sproc)%send_list(iq) = indp(i,j,n)
         end do
         do iq=1,bnds(sproc)%slenx_uv
            ! send_list(iq) is global point index, i, j, n are local
            ! Use abs because sign is used as u/v flag
            call indv_mpi(abs(bnds(sproc)%send_list_uv(iq)),i,j,n)
            bnds(sproc)%send_list_uv(iq) = sign(indp(i,j,n),bnds(sproc)%send_list_uv(iq))
         end do
      end do
      do iq=1,bnds(myid)%rlen2
         call indv_mpi(bnds(myid)%request_list(iq),i,j,n)
         bnds(myid)%request_list(iq) = indp(i,j,n)
      end do
      do iq=1,bnds(myid)%rlenx_uv
         call indv_mpi(abs(bnds(myid)%request_list_uv(iq)),i,j,n)
         bnds(myid)%request_list_uv(iq) = sign(indp(i,j,n),bnds(myid)%request_list_uv(iq))
      end do

      call reducealloc ! resize arrays

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
   
   ! MJT memory
   subroutine reducealloc
   
   implicit none
   
   integer iproc,nlen
   integer(kind=4) ierr
   integer(kind=4), dimension(maxbuflen) :: idum
   logical, dimension(maxbuflen) :: ldum
   
   do iproc=0,nproc-1
     nlen=max(kl+1,ol+1)*max(bnds(iproc)%rlen2,bnds(iproc)%rlenx_uv,bnds(iproc)%slen2,bnds(iproc)%slenx_uv)
     if (nlen<bnds(iproc)%len) then
       !write(6,*) "Reducing array size.  myid,iproc,nlen,len ",myid,iproc,nlen,bnds(iproc)%len
       bnds(iproc)%len=nlen
       if (iproc/=myid) then
         deallocate ( bnds(iproc)%rbuf )
         allocate ( bnds(iproc)%rbuf(bnds(iproc)%len) )
         deallocate ( bnds(iproc)%sbuf )
         allocate ( bnds(iproc)%sbuf(bnds(iproc)%len) )
         idum(1:bnds(iproc)%len)=bnds(iproc)%request_list(1:bnds(iproc)%len)
         deallocate ( bnds(iproc)%request_list )
         allocate ( bnds(iproc)%request_list(bnds(iproc)%len) )
         bnds(iproc)%request_list(1:bnds(iproc)%len)=idum(1:bnds(iproc)%len)
         idum(1:bnds(iproc)%len)=bnds(iproc)%send_list(1:bnds(iproc)%len)
         deallocate ( bnds(iproc)%send_list )
         allocate ( bnds(iproc)%send_list(bnds(iproc)%len) )
         bnds(iproc)%send_list(1:bnds(iproc)%len)=idum(1:bnds(iproc)%len)
         idum(1:bnds(iproc)%len)=bnds(iproc)%unpack_list(1:bnds(iproc)%len)
         deallocate ( bnds(iproc)%unpack_list )
         allocate ( bnds(iproc)%unpack_list(bnds(iproc)%len) )
         bnds(iproc)%unpack_list(1:bnds(iproc)%len)=idum(1:bnds(iproc)%len)
       end if
       idum(1:bnds(iproc)%len)=bnds(iproc)%request_list_uv(1:bnds(iproc)%len)
       deallocate ( bnds(iproc)%request_list_uv )
       allocate ( bnds(iproc)%request_list_uv(bnds(iproc)%len) )
       bnds(iproc)%request_list_uv(1:bnds(iproc)%len)=idum(1:bnds(iproc)%len)
       idum(1:bnds(iproc)%len)=bnds(iproc)%send_list_uv(1:bnds(iproc)%len)
       deallocate ( bnds(iproc)%send_list_uv )
       allocate ( bnds(iproc)%send_list_uv(bnds(iproc)%len) )
       bnds(iproc)%send_list_uv(1:bnds(iproc)%len)=idum(1:bnds(iproc)%len)
       idum(1:bnds(iproc)%len)=bnds(iproc)%unpack_list_uv(1:bnds(iproc)%len)
       deallocate ( bnds(iproc)%unpack_list_uv )
       allocate ( bnds(iproc)%unpack_list_uv(bnds(iproc)%len) )
       bnds(iproc)%unpack_list_uv(1:bnds(iproc)%len)=idum(1:bnds(iproc)%len)
       ldum(1:bnds(iproc)%len)=bnds(iproc)%uv_swap(1:bnds(iproc)%len)
       deallocate ( bnds(iproc)%uv_swap )
       allocate ( bnds(iproc)%uv_swap(bnds(iproc)%len) )
       bnds(iproc)%uv_swap(1:bnds(iproc)%len)=ldum(1:bnds(iproc)%len)
       ldum(1:bnds(iproc)%len)=bnds(iproc)%send_swap(1:bnds(iproc)%len)
       deallocate ( bnds(iproc)%send_swap )
       allocate ( bnds(iproc)%send_swap(bnds(iproc)%len) )
       bnds(iproc)%send_swap(1:bnds(iproc)%len)=ldum(1:bnds(iproc)%len)
       ldum(1:bnds(iproc)%len)=bnds(iproc)%uv_neg(1:bnds(iproc)%len)
       deallocate ( bnds(iproc)%uv_neg )
       allocate ( bnds(iproc)%uv_neg(bnds(iproc)%len) )
       bnds(iproc)%uv_neg(1:bnds(iproc)%len)=ldum(1:bnds(iproc)%len)
       ldum(1:bnds(iproc)%len)=bnds(iproc)%send_neg(1:bnds(iproc)%len)
       deallocate ( bnds(iproc)%send_neg )
       allocate ( bnds(iproc)%send_neg(bnds(iproc)%len) )
       bnds(iproc)%send_neg(1:bnds(iproc)%len)=ldum(1:bnds(iproc)%len)
     elseif (nlen>bnds(iproc)%len) then
       write(6,*) "ERROR reducing array size"
       write(6,*) "myid,iproc,nlen,len ",myid,iproc,nlen,bnds(iproc)%len
       write(6,*) "maxbuflen ",maxbuflen
       call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
     end if
   end do
   
   return
   end subroutine reducealloc

   subroutine check_set(ind,str,i,j,n,iq)
      integer, intent(in) :: ind,i,j,n,iq
      integer(kind=4) :: ierr
      character(len=*) :: str
      if ( ind == huge(1) ) then
         write(6,*) str, " not set", myid, i, j, n, iq
         call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
      end if
   end subroutine check_set

   subroutine bounds2(t, nrows, corner, nehalf, gmode, colour)
      ! Copy the boundary regions
      ! MJT - modified to restrict bounds calls to ocean processors with gmode=1
      ! MJT - modified to allow NE half updates to reduce message size
      ! MJT - modified to allow colours for SOR solver
      real, dimension(ifull+iextra), intent(inout) :: t
      integer, intent(in), optional :: nrows
      integer, intent(in), optional :: gmode
      logical, intent(in), optional :: corner
      logical, intent(in), optional :: nehalf
      integer, intent(in), optional :: colour
      logical :: double, extra, single
      logical :: fcrh, fcrf, fcgh, fcgf, fcbh, fcbf
      logical :: fextra, fdouble
      integer :: iq, iproc, iqz
      integer(kind=4) :: ierr, itag = 0, rproc, sproc, ltype
      integer(kind=4), dimension(MPI_STATUS_SIZE,2*nproc) :: status
      integer(kind=4) :: send_len, recv_len, iqq
      integer :: lmode, lcolour, msk_ne, msk_sw

      call start_log(bounds_begin)
      
      lmode=0
      lcolour=0 ! all colours
      if (present(gmode)) lmode=gmode
      if (present(colour)) lcolour=colour
#ifdef r8i8
      ltype=MPI_DOUBLE_PRECISION
#else
      ltype=MPI_REAL
#endif   

      if ( lcolour == 0 ) then
         t(ifull+1:ifull+iextra)=9.E9 ! MJT test for bad bounds call
      end if

      double = .false.
      extra = .false.
      single = .true.
      if (present(nrows)) then
         if ( nrows == 2 ) then
            double = .true.
         end if
      end if
      ! corner is irrelevant in double case
      if ( .not. double .and. present(corner) ) then
         extra = corner
      end if
      ! nehalf is irrelevant in double or corner case
      if ( .not. double .and. .not. extra .and. present(nehalf) ) then
         if (nehalf) single = .false.
      end if
      
      fcrh=.true.
      fcgh=.true.
      fcbh=.true.
      fcrf=.false.
      fcgf=.false.
      fcbf=.false.
      fextra=.false.
      fdouble=.false.
      if ( double ) then
         fcrf=.true.
         fcgf=.true.
         fcbf=.true.
         fextra=.true.
         fdouble=.true.
      else if ( extra ) then
         fcrf=.true.
         fcgf=.true.
         fcbf=.true.
         fextra=.true.
      else if ( single ) then
         fcrf=.true.
         fcgf=.true.
         fcbf=.true.
      end if
      select case ( lcolour )
         case( 1 )
           fcgh=.false.
           fcbh=.false.
           fcgf=.false.
           fcbf=.false.
         case( 2 )
           fcrh=.false.
           fcbh=.false.
           fcrf=.false.
           fcbf=.false.
         case( 3 )
           fcrh=.false.
           fcgh=.false.
           fcrf=.false.
           fcgf=.false.
      end select

!     Set up the buffers to send
      nreq = 0
      do iproc = 1,nproc-1  !
         rproc = modulo(myid-iproc,nproc)  ! Recv from
         if ( double ) then
            recv_len = bnds(rproc)%rlen2
         else if ( extra ) then
            recv_len = bnds(rproc)%rlenx
         else if ( single ) then
            recv_len = bnds(rproc)%rlen
         else
            recv_len = bnds(rproc)%rlenh
         end if
         ! processor mask
         msk_ne=1
         msk_sw=1
         if (lmode==1) then
            recv_len=recv_len*bnds(rproc)%mlomsk
            msk_ne=mlomsk_ne
            msk_sw=mlomsk_sw
         end if
         if ( recv_len > 0 ) then
            recv_len=0
            if ( fcrh ) then
               recv_len=recv_len+(rcolsp(rproc)%ihfn(1)-rcolsp(rproc)%ihbg(1)+1)*msk_ne
            end if
            if ( fcgh ) then
               recv_len=recv_len+(rcolsp(rproc)%ihfn(2)-rcolsp(rproc)%ihbg(2)+1)*msk_ne
            end if
            if ( fcbh ) then
               recv_len=recv_len+(rcolsp(rproc)%ihfn(3)-rcolsp(rproc)%ihbg(3)+1)*msk_ne
            end if
            if ( fcrf ) then
               recv_len=recv_len+(rcolsp(rproc)%iffn(1)-rcolsp(rproc)%ifbg(1)+1)*msk_sw
            end if
            if ( fcgf ) then
               recv_len=recv_len+(rcolsp(rproc)%iffn(2)-rcolsp(rproc)%ifbg(2)+1)*msk_sw
            end if
            if ( fcbf ) then
               recv_len=recv_len+(rcolsp(rproc)%iffn(3)-rcolsp(rproc)%ifbg(3)+1)*msk_sw
            end if
            if ( fextra ) then
               recv_len=recv_len+bnds(rproc)%rlenx-rcolsp(rproc)%iffn(3)
            end if
            if ( fdouble ) then
               recv_len=recv_len+bnds(rproc)%rlen2-bnds(rproc)%rlenx
            end if
            if ( recv_len > 0 ) then
               nreq = nreq + 1
               call MPI_IRecv( bnds(rproc)%rbuf(1),  recv_len, &
                    ltype, rproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
            end if
         end if
      end do
      do iproc = 1,nproc-1  !
         sproc = modulo(myid+iproc,nproc)  ! Send to
         if ( double ) then
            send_len = bnds(sproc)%slen2
         else if ( extra ) then
            send_len = bnds(sproc)%slenx
         else if ( single ) then
            send_len = bnds(sproc)%slen
         else
            send_len = bnds(sproc)%slenh
         end if
         ! processor mask
         msk_ne=1
         msk_sw=1
         if (lmode==1) then
            send_len=send_len*bnds(sproc)%mlomsk
            msk_ne=mlomsk_ne
            msk_sw=mlomsk_sw
         end if
         if ( send_len > 0 ) then
            ! Build up list of points
            iqq = 0
            if ( fcrh .and. msk_ne==1 ) then
!cdir nodep
               do iq=scolsp(sproc)%ihbg(1),scolsp(sproc)%ihfn(1)
                  iqz = iqq+iq-scolsp(sproc)%ihbg(1)+1
                  bnds(sproc)%sbuf(iqz) = t(bnds(sproc)%send_list(iq))
               end do
               iqq = iqq+scolsp(sproc)%ihfn(1)-scolsp(sproc)%ihbg(1)+1
            end if
             if ( fcgh .and. msk_ne==1 ) then
!cdir nodep
               do iq=scolsp(sproc)%ihbg(2),scolsp(sproc)%ihfn(2)
                  iqz = iqq+iq-scolsp(sproc)%ihbg(2)+1
                  bnds(sproc)%sbuf(iqz) = t(bnds(sproc)%send_list(iq))
               end do
               iqq = iqq+scolsp(sproc)%ihfn(2)-scolsp(sproc)%ihbg(2)+1
            end if           
            if ( fcbh .and. msk_ne==1 ) then
!cdir nodep
               do iq=scolsp(sproc)%ihbg(3),scolsp(sproc)%ihfn(3)
                  iqz = iqq+iq-scolsp(sproc)%ihbg(3)+1
                  bnds(sproc)%sbuf(iqz) = t(bnds(sproc)%send_list(iq))
               end do
               iqq = iqq+scolsp(sproc)%ihfn(3)-scolsp(sproc)%ihbg(3)+1
            end if            
            if ( fcrf .and. msk_sw==1 ) then
!cdir nodep
               do iq=scolsp(sproc)%ifbg(1),scolsp(sproc)%iffn(1)
                  iqz = iqq+iq-scolsp(sproc)%ifbg(1)+1
                  bnds(sproc)%sbuf(iqz) = t(bnds(sproc)%send_list(iq))
               end do
               iqq = iqq+scolsp(sproc)%iffn(1)-scolsp(sproc)%ifbg(1)+1
            end if
            if ( fcgf .and. msk_sw==1 ) then
!cdir nodep
               do iq=scolsp(sproc)%ifbg(2),scolsp(sproc)%iffn(2)
                  iqz = iqq+iq-scolsp(sproc)%ifbg(2)+1
                  bnds(sproc)%sbuf(iqz) = t(bnds(sproc)%send_list(iq))
               end do
               iqq = iqq+scolsp(sproc)%iffn(2)-scolsp(sproc)%ifbg(2)+1
            end if
            if ( fcbf .and. msk_sw==1 ) then
!cdir nodep
               do iq=scolsp(sproc)%ifbg(3),scolsp(sproc)%iffn(3)
                  iqz = iqq+iq-scolsp(sproc)%ifbg(3)+1
                  bnds(sproc)%sbuf(iqz) = t(bnds(sproc)%send_list(iq))
               end do
               iqq = iqq+scolsp(sproc)%iffn(3)-scolsp(sproc)%ifbg(3)+1
            end if
            if ( fextra ) then
!cdir nodep
                do iq=scolsp(sproc)%iffn(3)+1,bnds(sproc)%slenx
                  iqz = iqq+iq-scolsp(sproc)%iffn(3)
                  bnds(sproc)%sbuf(iqz) = t(bnds(sproc)%send_list(iq))
               end do
               iqq = iqq+bnds(sproc)%slenx-scolsp(sproc)%iffn(3)
            end if
            if ( fdouble ) then
!cdir nodep
                do iq=bnds(sproc)%slenx+1,bnds(sproc)%slen2
                  iqz = iqq+iq-bnds(sproc)%slenx
                  bnds(sproc)%sbuf(iqz) = t(bnds(sproc)%send_list(iq))
               end do
               iqq = iqq+bnds(sproc)%slen2-bnds(sproc)%slenx
            end if
            if ( iqq > 0 ) then
               nreq = nreq + 1
               call MPI_ISend( bnds(sproc)%sbuf(1), iqq, &
                    ltype, sproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
            end if
         end if
      end do

      ! Finally see if there are any points on my own processor that need
      ! to be fixed up. This will only be in the case when nproc < npanels.
      if ( double ) then
         recv_len = bnds(myid)%rlen2
      else if ( extra ) then
         recv_len = bnds(myid)%rlenx
      else if ( single ) then
         recv_len = bnds(myid)%rlen
      else
         recv_len = bnds(myid)%rlenh
      end if
!cdir nodep
      do iq=1,recv_len
         ! request_list is same as send_list in this case
         t(ifull+bnds(myid)%unpack_list(iq)) = t(bnds(myid)%request_list(iq))
      end do

      ! We have experimented with using MPI_Waitany.  The speed-up was
      ! marginal (i.e., roughly 1%), and did not justify the need to
      ! increase MPI_TYPE_MAX for large grids.
      if ( nreq > 0 ) then
         call start_log(mpiwait_begin)
         call MPI_Waitall(nreq,ireq,status,ierr)
         call end_log(mpiwait_end)
      end if

      do iproc = 1,nproc-1
         rproc = modulo(myid-iproc,nproc)  ! Recv from
         if ( double ) then
            recv_len = bnds(rproc)%rlen2
         else if ( extra ) then
            recv_len = bnds(rproc)%rlenx
         else if ( single ) then
            recv_len = bnds(rproc)%rlen
         else
            recv_len = bnds(rproc)%rlenh
         end if
         msk_ne=1
         msk_sw=1
         if (lmode==1) then
            recv_len=recv_len*bnds(rproc)%mlomsk
            msk_ne=mlomsk_ne
            msk_sw=mlomsk_sw
         end if
         if ( recv_len > 0 ) then
            ! unpack_list(iq) is index into extended region
            iqq = 0
            if ( fcrh .and. msk_ne==1 ) then
!cdir nodep
               do iq=rcolsp(rproc)%ihbg(1),rcolsp(rproc)%ihfn(1)
                  iqz = iqq+iq-rcolsp(rproc)%ihbg(1)+1
                  t(ifull+bnds(rproc)%unpack_list(iq)) = bnds(rproc)%rbuf(iqz)
               end do
               iqq = iqq+rcolsp(rproc)%ihfn(1)-rcolsp(rproc)%ihbg(1)+1
            end if
             if ( fcgh .and. msk_ne==1 ) then
!cdir nodep
               do iq=rcolsp(rproc)%ihbg(2),rcolsp(rproc)%ihfn(2)
                  iqz = iqq+iq-rcolsp(rproc)%ihbg(2)+1
                  t(ifull+bnds(rproc)%unpack_list(iq)) = bnds(rproc)%rbuf(iqz)
               end do
               iqq = iqq+rcolsp(rproc)%ihfn(2)-rcolsp(rproc)%ihbg(2)+1
            end if           
            if ( fcbh .and. msk_ne==1 ) then
!cdir nodep
               do iq=rcolsp(rproc)%ihbg(3),rcolsp(rproc)%ihfn(3)
                  iqz = iqq+iq-rcolsp(rproc)%ihbg(3)+1
                  t(ifull+bnds(rproc)%unpack_list(iq)) = bnds(rproc)%rbuf(iqz)
               end do
               iqq = iqq+rcolsp(rproc)%ihfn(3)-rcolsp(rproc)%ihbg(3)+1
            end if            
            if ( fcrf .and. msk_sw==1 ) then
!cdir nodep
               do iq=rcolsp(rproc)%ifbg(1),rcolsp(rproc)%iffn(1)
                  iqz = iqq+iq-rcolsp(rproc)%ifbg(1)+1
                  t(ifull+bnds(rproc)%unpack_list(iq)) = bnds(rproc)%rbuf(iqz)
               end do
               iqq = iqq+rcolsp(rproc)%iffn(1)-rcolsp(rproc)%ifbg(1)+1
            end if
            if ( fcgf .and. msk_sw==1 ) then
!cdir nodep
               do iq=rcolsp(rproc)%ifbg(2),rcolsp(rproc)%iffn(2)
                  iqz = iqq+iq-rcolsp(rproc)%ifbg(2)+1
                  t(ifull+bnds(rproc)%unpack_list(iq)) = bnds(rproc)%rbuf(iqz)
               end do
               iqq = iqq+rcolsp(rproc)%iffn(2)-rcolsp(rproc)%ifbg(2)+1
            end if
            if ( fcbf .and. msk_sw==1 ) then
!cdir nodep
               do iq=rcolsp(rproc)%ifbg(3),rcolsp(rproc)%iffn(3)
                  iqz = iqq+iq-rcolsp(rproc)%ifbg(3)+1
                  t(ifull+bnds(rproc)%unpack_list(iq)) = bnds(rproc)%rbuf(iqz)
               end do
               iqq = iqq+rcolsp(rproc)%iffn(3)-rcolsp(rproc)%ifbg(3)+1
            end if
            if ( fextra ) then
!cdir nodep
                do iq=rcolsp(rproc)%iffn(3)+1,bnds(rproc)%rlenx
                  iqz = iqq+iq-rcolsp(rproc)%iffn(3)
                  t(ifull+bnds(rproc)%unpack_list(iq)) = bnds(rproc)%rbuf(iqz)
               end do
               iqq = iqq+bnds(rproc)%rlenx-rcolsp(rproc)%iffn(3)
            end if
            if ( fdouble ) then
!cdir nodep
                do iq=bnds(rproc)%rlenx+1,bnds(rproc)%rlen2
                  iqz = iqq+iq-bnds(rproc)%rlenx
                  t(ifull+bnds(rproc)%unpack_list(iq)) = bnds(rproc)%rbuf(iqz)
               end do
               iqq = iqq+bnds(rproc)%rlen2-bnds(rproc)%rlenx
            end if
         end if         
      end do

      call end_log(bounds_end)

   end subroutine bounds2

   subroutine bounds3(t, nrows, klim, corner, nehalf, gmode, colour)
      ! Copy the boundary regions. Only this routine requires the extra klim
      ! argument (for helmsol).
      ! MJT - modified to restrict bounds calls to ocean processors with gmode=1
      ! MJT - modified to allow NE half updates to reduce message size
      ! MJT - modified to allow colours for SOR solver
      real, dimension(:,:), intent(inout) :: t
      integer, intent(in), optional :: nrows, klim
      integer, intent(in), optional :: gmode
      logical, intent(in), optional :: corner
      logical, intent(in), optional :: nehalf
      integer, intent(in), optional :: colour
      logical :: double, extra, single
      logical :: fcrh, fcrf, fcgh, fcgf, fcbh, fcbf
      logical :: fextra, fdouble
      integer :: iq, iproc, kx, iqz, iqb, iqe
      integer(kind=4) :: ierr, itag = 0, rproc, sproc, ltype
      integer(kind=4), dimension(MPI_STATUS_SIZE,2*nproc) :: status
      integer(kind=4) :: send_len, recv_len, iqq
      integer :: lmode, lcolour, msk_ne, msk_sw

      call start_log(bounds_begin)
      
      lmode=0
      lcolour=0 ! all colours
      if (present(gmode)) lmode=gmode
      if (present(colour)) lcolour=colour
#ifdef r8i8
      ltype=MPI_DOUBLE_PRECISION
#else
      ltype=MPI_REAL
#endif   

      if ( lcolour == 0 ) then
         t(ifull+1:ifull+iextra,:)=9.E9 ! MJT test for bad bounds call
      end if

      double = .false.
      extra = .false.
      single = .true.
      if (present(nrows)) then
         if ( nrows == 2 ) then
            double = .true.
         end if
      end if
      ! corner is irrelevant in double case
      if ( .not. double .and. present(corner) ) then
         extra = corner
      end if
      ! nehalf is irrelevant in double or corner case
      if ( .not. double .and. .not. extra .and. present(nehalf) ) then
         if (nehalf) single = .false.
      end if
      if ( present(klim) ) then
         kx = klim
      else
         kx = size(t,2)
      end if

      fcrh=.true.
      fcgh=.true.
      fcbh=.true.
      fcrf=.false.
      fcgf=.false.
      fcbf=.false.
      fextra=.false.
      fdouble=.false.
      if ( double ) then
         fcrf=.true.
         fcgf=.true.
         fcbf=.true.
         fextra=.true.
         fdouble=.true.
      else if ( extra ) then
         fcrf=.true.
         fcgf=.true.
         fcbf=.true.
         fextra=.true.
      else if ( single ) then
         fcrf=.true.
         fcgf=.true.
         fcbf=.true.
      end if
      select case ( lcolour )
         case( 1 )
           fcgh=.false.
           fcbh=.false.
           fcgf=.false.
           fcbf=.false.
         case( 2 )
           fcrh=.false.
           fcbh=.false.
           fcrf=.false.
           fcbf=.false.
         case( 3 )
           fcrh=.false.
           fcgh=.false.
           fcrf=.false.
           fcgf=.false.
      end select

!     Set up the buffers to send
      nreq = 0
      do iproc = 1,nproc-1  !
         rproc = modulo(myid-iproc,nproc)  ! Recv from
         if ( double ) then
            recv_len = bnds(rproc)%rlen2
         else if ( extra ) then
            recv_len = bnds(rproc)%rlenx
         else if ( single ) then
            recv_len = bnds(rproc)%rlen
         else
            recv_len = bnds(rproc)%rlenh
         end if
         msk_ne=1
         msk_sw=1
         if (lmode==1) then
            recv_len=recv_len*bnds(rproc)%mlomsk
            msk_ne=mlomsk_ne
            msk_sw=mlomsk_sw
         end if
         if ( recv_len > 0 ) then
            recv_len=0
            if ( fcrh ) then
               recv_len=recv_len+(rcolsp(rproc)%ihfn(1)-rcolsp(rproc)%ihbg(1)+1)*msk_ne
            end if
            if ( fcgh ) then
               recv_len=recv_len+(rcolsp(rproc)%ihfn(2)-rcolsp(rproc)%ihbg(2)+1)*msk_ne
            end if
            if ( fcbh ) then
               recv_len=recv_len+(rcolsp(rproc)%ihfn(3)-rcolsp(rproc)%ihbg(3)+1)*msk_ne
            end if
            if ( fcrf ) then
               recv_len=recv_len+(rcolsp(rproc)%iffn(1)-rcolsp(rproc)%ifbg(1)+1)*msk_sw
            end if
            if ( fcgf ) then
               recv_len=recv_len+(rcolsp(rproc)%iffn(2)-rcolsp(rproc)%ifbg(2)+1)*msk_sw
            end if
            if ( fcbf ) then
               recv_len=recv_len+(rcolsp(rproc)%iffn(3)-rcolsp(rproc)%ifbg(3)+1)*msk_sw
            end if
            if ( fextra ) then
               recv_len=recv_len+bnds(rproc)%rlenx-rcolsp(rproc)%iffn(3)
            end if
            if ( fdouble ) then
               recv_len=recv_len+bnds(rproc)%rlen2-bnds(rproc)%rlenx
            end if
            if ( recv_len > 0 ) then
               nreq = nreq + 1
               recv_len=recv_len*kx
               call MPI_IRecv( bnds(rproc)%rbuf(1), recv_len, &
                    ltype, rproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
            end if
         end if
      end do
      do iproc = 1,nproc-1  !
         sproc = modulo(myid+iproc,nproc)  ! Send to
         if ( double ) then
            send_len = bnds(sproc)%slen2
         else if ( extra ) then
            send_len = bnds(sproc)%slenx
         else if ( single ) then
            send_len = bnds(sproc)%slen
         else
            send_len = bnds(sproc)%slenh
         end if
         msk_ne=1
         msk_sw=1
         if (lmode==1) then
            send_len=send_len*bnds(sproc)%mlomsk
            msk_ne=mlomsk_ne
            msk_sw=mlomsk_sw
         end if
         if ( send_len > 0 ) then
            iqq = 0
            if ( fcrh .and. msk_ne==1 ) then
!cdir nodep
               do iq=scolsp(sproc)%ihbg(1),scolsp(sproc)%ihfn(1)
                  iqz = iqq+iq-scolsp(sproc)%ihbg(1)+1
                  iqb = 1+(iqz-1)*kx
                  iqe = iqz*kx
                  bnds(sproc)%sbuf(iqb:iqe) = t(bnds(sproc)%send_list(iq),:)
               end do
               iqq = iqq+scolsp(sproc)%ihfn(1)-scolsp(sproc)%ihbg(1)+1
            end if
             if ( fcgh .and. msk_ne==1 ) then
!cdir nodep
               do iq=scolsp(sproc)%ihbg(2),scolsp(sproc)%ihfn(2)
                  iqz = iqq+iq-scolsp(sproc)%ihbg(2)+1
                  iqb = 1+(iqz-1)*kx
                  iqe = iqz*kx
                  bnds(sproc)%sbuf(iqb:iqe) = t(bnds(sproc)%send_list(iq),:)
               end do
               iqq = iqq+scolsp(sproc)%ihfn(2)-scolsp(sproc)%ihbg(2)+1
            end if           
            if ( fcbh .and. msk_ne==1 ) then
!cdir nodep
               do iq=scolsp(sproc)%ihbg(3),scolsp(sproc)%ihfn(3)
                  iqz = iqq+iq-scolsp(sproc)%ihbg(3)+1
                  iqb = 1+(iqz-1)*kx
                  iqe = iqz*kx
                  bnds(sproc)%sbuf(iqb:iqe) = t(bnds(sproc)%send_list(iq),:)
               end do
               iqq = iqq+scolsp(sproc)%ihfn(3)-scolsp(sproc)%ihbg(3)+1
            end if            
            if ( fcrf .and. msk_sw==1 ) then
!cdir nodep
               do iq=scolsp(sproc)%ifbg(1),scolsp(sproc)%iffn(1)
                  iqz = iqq+iq-scolsp(sproc)%ifbg(1)+1
                  iqb = 1+(iqz-1)*kx
                  iqe = iqz*kx
                  bnds(sproc)%sbuf(iqb:iqe) = t(bnds(sproc)%send_list(iq),:)
               end do
               iqq = iqq+scolsp(sproc)%iffn(1)-scolsp(sproc)%ifbg(1)+1
            end if
            if ( fcgf .and. msk_sw==1 ) then
!cdir nodep
               do iq=scolsp(sproc)%ifbg(2),scolsp(sproc)%iffn(2)
                  iqz = iqq+iq-scolsp(sproc)%ifbg(2)+1
                  iqb = 1+(iqz-1)*kx
                  iqe = iqz*kx
                  bnds(sproc)%sbuf(iqb:iqe) = t(bnds(sproc)%send_list(iq),:)
               end do
               iqq = iqq+scolsp(sproc)%iffn(2)-scolsp(sproc)%ifbg(2)+1
            end if
            if ( fcbf .and. msk_sw==1 ) then
!cdir nodep
               do iq=scolsp(sproc)%ifbg(3),scolsp(sproc)%iffn(3)
                  iqz = iqq+iq-scolsp(sproc)%ifbg(3)+1
                  iqb = 1+(iqz-1)*kx
                  iqe = iqz*kx
                  bnds(sproc)%sbuf(iqb:iqe) = t(bnds(sproc)%send_list(iq),:)
               end do
               iqq = iqq+scolsp(sproc)%iffn(3)-scolsp(sproc)%ifbg(3)+1
            end if
            if ( fextra ) then
!cdir nodep
                do iq=scolsp(sproc)%iffn(3)+1,bnds(sproc)%slenx
                  iqz = iqq+iq-scolsp(sproc)%iffn(3)
                  iqb = 1+(iqz-1)*kx
                  iqe = iqz*kx
                  bnds(sproc)%sbuf(iqb:iqe) = t(bnds(sproc)%send_list(iq),:)
               end do
               iqq = iqq+bnds(sproc)%slenx-scolsp(sproc)%iffn(3)
            end if
            if ( fdouble ) then
!cdir nodep
                do iq=bnds(sproc)%slenx+1,bnds(sproc)%slen2
                  iqz = iqq+iq-bnds(sproc)%slenx
                  iqb = 1+(iqz-1)*kx
                  iqe = iqz*kx
                  bnds(sproc)%sbuf(iqb:iqe) = t(bnds(sproc)%send_list(iq),:)
               end do
               iqq = iqq+bnds(sproc)%slen2-bnds(sproc)%slenx
            end if
            if ( iqq > 0 ) then
               nreq = nreq + 1
               iqq=iqq*kx
               call MPI_ISend( bnds(sproc)%sbuf(1), iqq, &
                    ltype, sproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
            end if
         end if
      end do

      ! Finally see if there are any points on my own processor that need
      ! to be fixed up. This will only be in the case when nproc < npanels.
      if ( double ) then
         recv_len = bnds(myid)%rlen2
      else if ( extra ) then
         recv_len = bnds(myid)%rlenx
      else if ( single ) then
         recv_len = bnds(myid)%rlen
      else
         recv_len = bnds(myid)%rlenh
      end if
!cdir nodep
      do iq=1,recv_len
         ! request_list is same as send_list in this case
         t(ifull+bnds(myid)%unpack_list(iq),:) = t(bnds(myid)%request_list(iq),:)
      end do

      if ( nreq > 0 ) then
         call start_log(mpiwait_begin)
         call MPI_Waitall(nreq,ireq,status,ierr)
         call end_log(mpiwait_end)
      end if

      do iproc = 1,nproc-1
         rproc = modulo(myid-iproc,nproc)  ! Recv from
         if ( double ) then
            recv_len = bnds(rproc)%rlen2
         else if ( extra ) then
            recv_len = bnds(rproc)%rlenx
         else if ( single ) then
            recv_len = bnds(rproc)%rlen
         else
            recv_len = bnds(rproc)%rlenh
         end if
         msk_ne=1
         msk_sw=1
         if (lmode==1) then
            recv_len=recv_len*bnds(rproc)%mlomsk
            msk_ne=mlomsk_ne
            msk_sw=mlomsk_sw
         end if
         if ( recv_len > 0 ) then
            iqq = 0
            if ( fcrh .and. msk_ne==1 ) then
!cdir nodep
               do iq=rcolsp(rproc)%ihbg(1),rcolsp(rproc)%ihfn(1)
                  iqz = iqq+iq-rcolsp(rproc)%ihbg(1)+1
                  iqb = 1+(iqz-1)*kx
                  iqe = iqz*kx
                  t(ifull+bnds(rproc)%unpack_list(iq),:) = bnds(rproc)%rbuf(iqb:iqe)
               end do
               iqq = iqq+rcolsp(rproc)%ihfn(1)-rcolsp(rproc)%ihbg(1)+1
            end if
             if ( fcgh .and. msk_ne==1 ) then
!cdir nodep
               do iq=rcolsp(rproc)%ihbg(2),rcolsp(rproc)%ihfn(2)
                  iqz = iqq+iq-rcolsp(rproc)%ihbg(2)+1
                  iqb = 1+(iqz-1)*kx
                  iqe = iqz*kx
                  t(ifull+bnds(rproc)%unpack_list(iq),:) = bnds(rproc)%rbuf(iqb:iqe)
               end do
               iqq = iqq+rcolsp(rproc)%ihfn(2)-rcolsp(rproc)%ihbg(2)+1
            end if           
            if ( fcbh .and. msk_ne==1 ) then
!cdir nodep
               do iq=rcolsp(rproc)%ihbg(3),rcolsp(rproc)%ihfn(3)
                  iqz = iqq+iq-rcolsp(rproc)%ihbg(3)+1
                  iqb = 1+(iqz-1)*kx
                  iqe = iqz*kx
                  t(ifull+bnds(rproc)%unpack_list(iq),:) = bnds(rproc)%rbuf(iqb:iqe)
               end do
               iqq = iqq+rcolsp(rproc)%ihfn(3)-rcolsp(rproc)%ihbg(3)+1
            end if            
            if ( fcrf .and. msk_sw==1 ) then
!cdir nodep
               do iq=rcolsp(rproc)%ifbg(1),rcolsp(rproc)%iffn(1)
                  iqz = iqq+iq-rcolsp(rproc)%ifbg(1)+1
                  iqb = 1+(iqz-1)*kx
                  iqe = iqz*kx
                  t(ifull+bnds(rproc)%unpack_list(iq),:) = bnds(rproc)%rbuf(iqb:iqe)
               end do
               iqq = iqq+rcolsp(rproc)%iffn(1)-rcolsp(rproc)%ifbg(1)+1
            end if
            if ( fcgf .and. msk_sw==1 ) then
!cdir nodep
               do iq=rcolsp(rproc)%ifbg(2),rcolsp(rproc)%iffn(2)
                  iqz = iqq+iq-rcolsp(rproc)%ifbg(2)+1
                  iqb = 1+(iqz-1)*kx
                  iqe = iqz*kx
                  t(ifull+bnds(rproc)%unpack_list(iq),:) = bnds(rproc)%rbuf(iqb:iqe)
               end do
               iqq = iqq+rcolsp(rproc)%iffn(2)-rcolsp(rproc)%ifbg(2)+1
            end if
            if ( fcbf .and. msk_sw==1 ) then
!cdir nodep
               do iq=rcolsp(rproc)%ifbg(3),rcolsp(rproc)%iffn(3)
                  iqz = iqq+iq-rcolsp(rproc)%ifbg(3)+1
                  iqb = 1+(iqz-1)*kx
                  iqe = iqz*kx
                  t(ifull+bnds(rproc)%unpack_list(iq),:) = bnds(rproc)%rbuf(iqb:iqe)
               end do
               iqq = iqq+rcolsp(rproc)%iffn(3)-rcolsp(rproc)%ifbg(3)+1
            end if
            if ( fextra ) then
!cdir nodep
                do iq=rcolsp(rproc)%iffn(3)+1,bnds(rproc)%rlenx
                  iqz = iqq+iq-rcolsp(rproc)%iffn(3)
                  iqb = 1+(iqz-1)*kx
                  iqe = iqz*kx
                  t(ifull+bnds(rproc)%unpack_list(iq),:) = bnds(rproc)%rbuf(iqb:iqe)
               end do
               iqq = iqq+bnds(rproc)%rlenx-rcolsp(rproc)%iffn(3)
            end if
            if ( fdouble ) then
!cdir nodep
                do iq=bnds(rproc)%rlenx+1,bnds(rproc)%rlen2
                  iqz = iqq+iq-bnds(rproc)%rlenx
                  iqb = 1+(iqz-1)*kx
                  iqe = iqz*kx
                  t(ifull+bnds(rproc)%unpack_list(iq),:) = bnds(rproc)%rbuf(iqb:iqe)
               end do
               iqq = iqq+bnds(rproc)%rlen2-bnds(rproc)%rlenx
            end if
         end if
      end do

      call end_log(bounds_end)

   end subroutine bounds3

   subroutine boundsuv2(u, v, nrows, stag, allvec, gmode)
      ! Copy the boundary regions of u and v. This doesn't require the
      ! diagonal points like (0,0), but does have to take care of the
      ! direction changes.
      ! MJT - Modified to send smaller message lengths for staguv and diffusion
      ! MJT - modified to restrict bounds calls to ocean processors with gmode=1
      real, dimension(ifull+iextra), intent(inout) :: u, v
      integer, intent(in), optional :: nrows
      integer, intent(in), optional :: stag
      integer, intent(in), optional :: gmode
      logical, intent(in), optional :: allvec
      logical :: double, extra
      logical :: fsvwu, fnveu, fssvwwu, fnnveeu
      logical :: fsuwv, fnuev
      integer :: stagmode, iq, iqz, iproc
      integer(kind=4) :: ierr, itag = 0, rproc, sproc, iqq, ltype
      integer(kind=4), dimension(MPI_STATUS_SIZE,2*nproc) :: status
      integer(kind=4) :: send_len, recv_len
      integer :: lmode, msk_ne, msk_sw
      real :: tmp

      call start_log(boundsuv_begin)
      
      lmode=0
      if (present(gmode)) lmode=gmode
#ifdef r8i8
      ltype=MPI_DOUBLE_PRECISION
#else
      ltype=MPI_REAL
#endif   

      u(ifull+1:ifull+iextra)=9.E9 ! MJT test for bad bounds call
      v(ifull+1:ifull+iextra)=9.E9 ! MJT test for bad bounds call

      double = .false.
      extra = .false.
      stagmode = 0
      if ( present(allvec) ) then
        extra = allvec
      end if
      ! double is irrelevant in extra case
      if ( .not. extra .and. present(nrows)) then
         if ( nrows == 2 ) double = .true.
      end if
      if ( present(stag) ) then
        stagmode = stag
      end if

      fsvwu = .true.
      fnveu = .true.
      fssvwwu = .false.
      fnnveeu = .false.
      fsuwv = .false.
      fnuev = .false.
      if ( double ) then
        fssvwwu = .true.
        fnnveeu = .true.
      else if ( extra ) then
        fsuwv = .true.
        fnuev = .true.      
      else if ( stagmode == 1 ) then
        fsvwu = .false.
        fnnveeu = .true.
      else if ( stagmode == 2 ) then
        fnnveeu = .true. ! fnnveeu requires fnveu
      else if ( stagmode == 3 ) then
        fssvwwu = .true. ! fssvwwu requires fsvwu
      else if ( stagmode == 5 ) then
        fnveu = .false.
        fssvwwu = .true.
      else if ( stagmode == -9 ) then
        fnveu = .false.
      else if ( stagmode == -10 ) then
        fsvwu = .false.
      else if ( stagmode == -19 ) then
        fsvwu = .false.
        fnveu = .false.
        fnuev = .true.
      end if
      
!     Set up the buffers to send
      nreq = 0
      do iproc = 1,nproc-1  !
         rproc = modulo(myid-iproc,nproc)  ! Recv from
         if ( extra .or. stagmode == -19 ) then
            recv_len = bnds(rproc)%rlenx_uv
         else if ( double .or. stagmode > 0 ) then
            recv_len = bnds(rproc)%rlen2_uv
         else
            recv_len = bnds(rproc)%rlen_uv
         end if
         msk_ne=1
         msk_sw=1
         if (lmode==1) then
            recv_len=recv_len*bnds(rproc)%mlomsk
            msk_ne=mlomsk_ne
            msk_sw=mlomsk_sw
         end if
         if ( recv_len /= 0 ) then
            recv_len=0
            if ( fsvwu ) then
               recv_len=recv_len+(rsplit(rproc)%iwufn-rsplit(rproc)%isvbg+1)*msk_sw
            end if
            if ( fnveu ) then
               recv_len=recv_len+(rsplit(rproc)%ieufn-rsplit(rproc)%invbg+1)*msk_ne
            end if         
            if ( fssvwwu ) then
               recv_len=recv_len+(rsplit(rproc)%iwwufn-rsplit(rproc)%issvbg+1)*msk_sw
            end if         
            if ( fnnveeu ) then
               recv_len=recv_len+(rsplit(rproc)%ieeufn-rsplit(rproc)%innvbg+1)*msk_ne
            end if         
            if ( fsuwv ) then
               recv_len=recv_len+(rsplit(rproc)%iwvfn-rsplit(rproc)%isubg+1)*msk_sw
            end if
            if ( fnuev ) then
               recv_len=recv_len+(rsplit(rproc)%ievfn-rsplit(rproc)%inubg+1)*msk_ne
            end if
            if ( recv_len > 0 ) then 
               nreq = nreq + 1
               call MPI_IRecv( bnds(rproc)%rbuf(1), recv_len, &
                    ltype, rproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
            end if
         end if
      end do
      do iproc = 1,nproc-1  !
         sproc = modulo(myid+iproc,nproc)  ! Send to
         if ( extra .or. stagmode == -19 ) then
            send_len = bnds(sproc)%slenx_uv
         else if ( double .or. stagmode > 0 ) then
            send_len = bnds(sproc)%slen2_uv
         else
            send_len = bnds(sproc)%slen_uv
         end if
         msk_ne=1
         msk_sw=1
         if (lmode==1) then
            send_len=send_len*bnds(sproc)%mlomsk
            msk_ne=mlomsk_ne
            msk_sw=mlomsk_sw
         end if
         if ( send_len > 0 ) then
            ! Build up list of points
            iqq = 0
            if ( fsvwu .and. msk_sw==1 ) then
!cdir nodep
               do iq=ssplit(sproc)%isvbg,ssplit(sproc)%iwufn
                  ! Use abs because sign is used as u/v flag
                  iqz = iqq+iq-ssplit(sproc)%isvbg+1
                  if ( (bnds(sproc)%send_list_uv(iq) > 0) .neqv. &
                        bnds(sproc)%send_swap(iq) ) then
                     bnds(sproc)%sbuf(iqz) = u(abs(bnds(sproc)%send_list_uv(iq)))
                  else
                     bnds(sproc)%sbuf(iqz) = v(abs(bnds(sproc)%send_list_uv(iq)))
                  end if 
                  if ( bnds(sproc)%send_neg(iq) ) bnds(sproc)%sbuf(iqz) = -bnds(sproc)%sbuf(iqz)
               end do
               iqq = iqq+ssplit(sproc)%iwufn-ssplit(sproc)%isvbg+1
            end if
            if ( fnveu .and. msk_ne==1 ) then
!cdir nodep
               do iq=ssplit(sproc)%invbg,ssplit(sproc)%ieufn
                  ! Use abs because sign is used as u/v flag
                  iqz = iqq+iq-ssplit(sproc)%invbg+1
                  if ( (bnds(sproc)%send_list_uv(iq) > 0) .neqv. &
                        bnds(sproc)%send_swap(iq) ) then
                     bnds(sproc)%sbuf(iqz) = u(abs(bnds(sproc)%send_list_uv(iq)))
                  else
                     bnds(sproc)%sbuf(iqz) = v(abs(bnds(sproc)%send_list_uv(iq)))
                  end if 
                  if ( bnds(sproc)%send_neg(iq) ) bnds(sproc)%sbuf(iqz) = -bnds(sproc)%sbuf(iqz)
               end do
               iqq = iqq+ssplit(sproc)%ieufn-ssplit(sproc)%invbg+1
            end if
            if ( fssvwwu .and. msk_sw==1 ) then
!cdir nodep
               do iq=ssplit(sproc)%issvbg,ssplit(sproc)%iwwufn
                  ! Use abs because sign is used as u/v flag
                  iqz = iqq+iq-ssplit(sproc)%issvbg+1
                  if ( (bnds(sproc)%send_list_uv(iq) > 0) .neqv. &
                        bnds(sproc)%send_swap(iq) ) then
                     bnds(sproc)%sbuf(iqz) = u(abs(bnds(sproc)%send_list_uv(iq)))
                  else
                     bnds(sproc)%sbuf(iqz) = v(abs(bnds(sproc)%send_list_uv(iq)))
                  end if 
                  if ( bnds(sproc)%send_neg(iq) ) bnds(sproc)%sbuf(iqz) = -bnds(sproc)%sbuf(iqz)
               end do
               iqq = iqq+ssplit(sproc)%iwwufn-ssplit(sproc)%issvbg+1
            end if
            if ( fnnveeu .and. msk_ne==1 ) then
!cdir nodep
               do iq=ssplit(sproc)%innvbg,ssplit(sproc)%ieeufn
                  ! Use abs because sign is used as u/v flag
                  iqz = iqq+iq-ssplit(sproc)%innvbg+1
                  if ( (bnds(sproc)%send_list_uv(iq) > 0) .neqv. &
                        bnds(sproc)%send_swap(iq) ) then
                     bnds(sproc)%sbuf(iqz) = u(abs(bnds(sproc)%send_list_uv(iq)))
                  else
                     bnds(sproc)%sbuf(iqz) = v(abs(bnds(sproc)%send_list_uv(iq)))
                  end if 
                  if ( bnds(sproc)%send_neg(iq) ) bnds(sproc)%sbuf(iqz) = -bnds(sproc)%sbuf(iqz)
               end do
               iqq = iqq+ssplit(sproc)%ieeufn-ssplit(sproc)%innvbg+1
            end if
            if ( fsuwv .and. msk_sw==1 ) then
!cdir nodep
               do iq=ssplit(sproc)%isubg,ssplit(sproc)%iwvfn
                  ! Use abs because sign is used as u/v flag
                  iqz = iqq+iq-ssplit(sproc)%isubg+1
                  if ( (bnds(sproc)%send_list_uv(iq) > 0) .neqv. &
                        bnds(sproc)%send_swap(iq) ) then
                     bnds(sproc)%sbuf(iqz) = u(abs(bnds(sproc)%send_list_uv(iq)))
                  else
                     bnds(sproc)%sbuf(iqz) = v(abs(bnds(sproc)%send_list_uv(iq)))
                  end if 
                  if ( bnds(sproc)%send_neg(iq) ) bnds(sproc)%sbuf(iqz) = -bnds(sproc)%sbuf(iqz)
               end do
               iqq = iqq+ssplit(sproc)%iwvfn-ssplit(sproc)%isubg+1
            end if
            if ( fnuev .and. msk_ne==1 ) then
!cdir nodep
               do iq=ssplit(sproc)%inubg,ssplit(sproc)%ievfn
                  ! Use abs because sign is used as u/v flag
                  iqz = iqq+iq-ssplit(sproc)%inubg+1
                  if ( (bnds(sproc)%send_list_uv(iq) > 0) .neqv. &
                        bnds(sproc)%send_swap(iq) ) then
                     bnds(sproc)%sbuf(iqz) = u(abs(bnds(sproc)%send_list_uv(iq)))
                  else
                     bnds(sproc)%sbuf(iqz) = v(abs(bnds(sproc)%send_list_uv(iq)))
                  end if 
                  if ( bnds(sproc)%send_neg(iq) ) bnds(sproc)%sbuf(iqz) = -bnds(sproc)%sbuf(iqz)
               end do
               iqq = iqq+ssplit(sproc)%ievfn-ssplit(sproc)%inubg+1
            end if
            if ( iqq > 0 ) then
               nreq = nreq + 1
               call MPI_ISend( bnds(sproc)%sbuf(1), iqq, &
                    ltype, sproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
            end if
         end if
      end do

      ! Finally see if there are any points on my own processor that need
      ! to be fixed up. This will only be in the case when nproc < npanels.
      if ( bnds(myid)%rlen_uv > 0 ) then
         if ( extra .or. stagmode == -19 ) then
            recv_len = bnds(myid)%rlenx_uv
         else if ( double .or. stagmode > 0 ) then
            recv_len = bnds(myid)%rlen2_uv
         else
            recv_len = bnds(myid)%rlen_uv
         end if
!cdir nodep
         do iq=1,recv_len
            ! request_list is same as send_list in this case
            if ( ( bnds(myid)%request_list_uv(iq) > 0) .neqv. &
                      bnds(myid)%uv_swap(iq) ) then  ! haven't copied to send_swap yet
               tmp = u(abs(bnds(myid)%request_list_uv(iq)))
            else
               tmp = v(abs(bnds(myid)%request_list_uv(iq)))
            end if
            if ( bnds(myid)%uv_neg(iq) ) tmp = -tmp
            ! unpack_list(iq) is index into extended region
            if ( bnds(myid)%unpack_list_uv(iq) > 0 ) then
               u(ifull+bnds(myid)%unpack_list_uv(iq)) = tmp
            else
               v(ifull-bnds(myid)%unpack_list_uv(iq)) = tmp
            end if
         end do
      end if

      if ( nreq > 0 ) then
         call start_log(mpiwait_begin)
         call MPI_Waitall(nreq,ireq,status,ierr)
         call end_log(mpiwait_end)
      end if

      do iproc = 1,nproc-1  !
         rproc = modulo(myid-iproc,nproc)  ! Recv from
         if ( extra .or. stagmode == -19 ) then
            recv_len = bnds(rproc)%rlenx_uv
         else if ( double .or. stagmode > 0 ) then
            recv_len = bnds(rproc)%rlen2_uv
         else
            recv_len = bnds(rproc)%rlen_uv
         end if
         msk_ne=1
         msk_sw=1
         if (lmode==1) then
            recv_len=recv_len*bnds(rproc)%mlomsk
            msk_ne=mlomsk_ne
            msk_sw=mlomsk_sw
         end if
         iqq = 0
         if ( fsvwu .and. msk_sw==1 ) then
!cdir nodep
            do iq=rsplit(rproc)%isvbg,rsplit(rproc)%iwufn
               ! unpack_list(iq) is index into extended region
               iqz = iqq+iq-rsplit(rproc)%isvbg+1
               if ( bnds(rproc)%unpack_list_uv(iq) > 0 ) then
                  u(ifull+bnds(rproc)%unpack_list_uv(iq)) = bnds(rproc)%rbuf(iqz)
               else
                  v(ifull-bnds(rproc)%unpack_list_uv(iq)) = bnds(rproc)%rbuf(iqz)
               end if
            end do
            iqq = iqq+rsplit(rproc)%iwufn-rsplit(rproc)%isvbg+1
         end if
         if ( fnveu .and. msk_ne==1 ) then
!cdir nodep
            do iq=rsplit(rproc)%invbg,rsplit(rproc)%ieufn
               ! unpack_list(iq) is index into extended region
               iqz = iqq+iq-rsplit(rproc)%invbg+1
               if ( bnds(rproc)%unpack_list_uv(iq) > 0 ) then
                  u(ifull+bnds(rproc)%unpack_list_uv(iq)) = bnds(rproc)%rbuf(iqz)
               else
                  v(ifull-bnds(rproc)%unpack_list_uv(iq)) = bnds(rproc)%rbuf(iqz)
               end if
            end do
            iqq = iqq+rsplit(rproc)%ieufn-rsplit(rproc)%invbg+1
         end if
         if ( fssvwwu .and. msk_sw==1 ) then
!cdir nodep
            do iq=rsplit(rproc)%issvbg,rsplit(rproc)%iwwufn
               ! unpack_list(iq) is index into extended region
               iqz = iqq+iq-rsplit(rproc)%issvbg+1
               if ( bnds(rproc)%unpack_list_uv(iq) > 0 ) then
                  u(ifull+bnds(rproc)%unpack_list_uv(iq)) = bnds(rproc)%rbuf(iqz)
               else
                  v(ifull-bnds(rproc)%unpack_list_uv(iq)) = bnds(rproc)%rbuf(iqz)
               end if
            end do
            iqq = iqq+rsplit(rproc)%iwwufn-rsplit(rproc)%issvbg+1
         end if
         if ( fnnveeu .and. msk_ne==1 ) then
!cdir nodep
            do iq=rsplit(rproc)%innvbg,rsplit(rproc)%ieeufn
               ! unpack_list(iq) is index into extended region
               iqz = iqq+iq-rsplit(rproc)%innvbg+1
               if ( bnds(rproc)%unpack_list_uv(iq) > 0 ) then
                  u(ifull+bnds(rproc)%unpack_list_uv(iq)) = bnds(rproc)%rbuf(iqz)
               else
                  v(ifull-bnds(rproc)%unpack_list_uv(iq)) = bnds(rproc)%rbuf(iqz)
               end if
            end do
            iqq = iqq+rsplit(rproc)%ieeufn-rsplit(rproc)%innvbg+1
         end if
         if ( fsuwv .and. msk_sw==1 ) then
!cdir nodep
            do iq=rsplit(rproc)%isubg,rsplit(rproc)%iwvfn
               ! unpack_list(iq) is index into extended region
               iqz = iqq+iq-rsplit(rproc)%isubg+1
               if ( bnds(rproc)%unpack_list_uv(iq) > 0 ) then
                  u(ifull+bnds(rproc)%unpack_list_uv(iq)) = bnds(rproc)%rbuf(iqz)
               else
                  v(ifull-bnds(rproc)%unpack_list_uv(iq)) = bnds(rproc)%rbuf(iqz)
               end if
            end do
            iqq = iqq+rsplit(rproc)%iwvfn-rsplit(rproc)%isubg+1
         end if
         if ( fnuev .and. msk_ne==1 ) then
!cdir nodep
            do iq=rsplit(rproc)%inubg,rsplit(rproc)%ievfn
               ! unpack_list(iq) is index into extended region
               iqz = iqq+iq-rsplit(rproc)%inubg+1
               if ( bnds(rproc)%unpack_list_uv(iq) > 0 ) then
                  u(ifull+bnds(rproc)%unpack_list_uv(iq)) = bnds(rproc)%rbuf(iqz)
               else
                  v(ifull-bnds(rproc)%unpack_list_uv(iq)) = bnds(rproc)%rbuf(iqz)
               end if
            end do
            iqq = iqq+rsplit(rproc)%ievfn-rsplit(rproc)%inubg+1
         end if
      end do

      call end_log(boundsuv_end)

   end subroutine boundsuv2

   subroutine boundsuv3(u, v, nrows, stag, allvec, gmode)
      ! Copy the boundary regions of u and v. This doesn't require the
      ! diagonal points like (0,0), but does have to take care of the
      ! direction changes.
      ! MJT - Modified to send smaller message lengths for staguv and diffusion
      ! MJT - modified to restrict bounds calls to ocean processors with gmode=1
      real, dimension(:,:), intent(inout) :: u, v
      integer, intent(in), optional :: nrows
      integer, intent(in), optional :: stag
      integer, intent(in), optional :: gmode
      logical, intent(in), optional :: allvec
      logical :: double, extra
      logical :: fsvwu, fnveu, fssvwwu, fnnveeu
      logical :: fsuwv, fnuev
      integer :: stagmode, iq, iqz, iqb, iqe, iproc, kx
      integer(kind=4) :: ierr, itag = 0, rproc, sproc, iqq, ltype
      integer(kind=4), dimension(MPI_STATUS_SIZE,2*nproc) :: status
      integer(kind=4) :: send_len, recv_len
      integer :: lmode, msk_ne, msk_sw
      real, dimension(maxbuflen) :: tmp
      
      call start_log(boundsuv_begin)

      lmode=0
      if (present(gmode)) lmode=gmode
#ifdef r8i8
      ltype=MPI_DOUBLE_PRECISION
#else
      ltype=MPI_REAL
#endif   

      kx=size(u,2)
      u(ifull+1:ifull+iextra,:)=9.E9 ! MJT test for bad bounds call
      v(ifull+1:ifull+iextra,:)=9.E9 ! MJT test for bad bounds call

      double = .false.
      extra = .false.
      stagmode = 0
      if ( present(allvec) ) then
        extra=allvec
      end if 
      ! double is irrelevant in extra case
      if ( .not. extra .and. present(nrows) ) then
         if ( nrows == 2 ) double = .true.
      end if
      if ( present(stag) ) then
        stagmode = stag
      end if
      
      fsvwu = .true.
      fnveu = .true.
      fssvwwu = .false.
      fnnveeu = .false.
      fsuwv = .false.
      fnuev = .false.
      if ( double ) then
        fssvwwu = .true.
        fnnveeu = .true.
      else if ( extra ) then
        fsuwv = .true.
        fnuev = .true.      
      else if ( stagmode == 1 ) then
        fsvwu = .false.
        fnnveeu = .true.
      else if ( stagmode == 2 ) then
        fnnveeu = .true. ! fnnveeu requires fnveu
      else if ( stagmode == 3 ) then
        fssvwwu = .true. ! fssvwwu requires fsvwu
      else if ( stagmode == 5 ) then
        fnveu = .false.
        fssvwwu = .true.
      else if ( stagmode == -9 ) then
        fnveu = .false.
      else if ( stagmode == -10 ) then
        fsvwu = .false.
      else if ( stagmode == -19 ) then
        fsvwu = .false.
        fnveu = .false.
        fnuev = .true.
      end if

!     Set up the buffers to send
      nreq = 0
      do iproc = 1,nproc-1  !
         rproc = modulo(myid-iproc,nproc)  ! Recv from
         if ( extra .or. stagmode == -19 ) then
            recv_len = bnds(rproc)%rlenx_uv
         else if ( double .or. stagmode > 0 ) then
            recv_len = bnds(rproc)%rlen2_uv
         else
            recv_len = bnds(rproc)%rlen_uv
         end if
         msk_ne=1
         msk_sw=1
         if (lmode==1) then
            recv_len=recv_len*bnds(rproc)%mlomsk
            msk_ne=mlomsk_ne
            msk_sw=mlomsk_sw
         end if
         if ( recv_len > 0 ) then
            recv_len=0
            if ( fsvwu ) then
               recv_len=recv_len+(rsplit(rproc)%iwufn-rsplit(rproc)%isvbg+1)*msk_sw
            end if
            if ( fnveu ) then
               recv_len=recv_len+(rsplit(rproc)%ieufn-rsplit(rproc)%invbg+1)*msk_ne
            end if         
            if ( fssvwwu ) then
               recv_len=recv_len+(rsplit(rproc)%iwwufn-rsplit(rproc)%issvbg+1)*msk_sw
            end if         
            if ( fnnveeu ) then
               recv_len=recv_len+(rsplit(rproc)%ieeufn-rsplit(rproc)%innvbg+1)*msk_ne
            end if         
            if ( fsuwv ) then
               recv_len=recv_len+(rsplit(rproc)%iwvfn-rsplit(rproc)%isubg+1)*msk_sw
            end if
            if ( fnuev ) then
               recv_len=recv_len+(rsplit(rproc)%ievfn-rsplit(rproc)%inubg+1)*msk_ne
            end if
            if ( recv_len > 0 ) then 
               nreq = nreq + 1
               recv_len=recv_len*kx
               call MPI_IRecv( bnds(rproc)%rbuf(1), recv_len, &
                    ltype, rproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
            end if
         end if
      end do
      do iproc = 1,nproc-1  !
         sproc = modulo(myid+iproc,nproc)  ! Send to
         if ( extra .or. stagmode == -19 ) then
            send_len = bnds(sproc)%slenx_uv
         else if ( double .or. stagmode > 0 ) then
            send_len = bnds(sproc)%slen2_uv
         else
            send_len = bnds(sproc)%slen_uv
         end if
         msk_ne=1
         msk_sw=1
         if (lmode==1) then
            send_len=send_len*bnds(sproc)%mlomsk
            msk_ne=mlomsk_ne
            msk_sw=mlomsk_sw
         end if
         if ( send_len > 0 ) then
            ! Build up list of points
            iqq = 0
            if ( fsvwu .and. msk_sw==1 ) then
!cdir nodep
               do iq=ssplit(sproc)%isvbg,ssplit(sproc)%iwufn
                  ! send_list_uv(iq) is point index.
                  ! Use abs because sign is used as u/v flag
                  iqz = iqq+iq-ssplit(sproc)%isvbg+1
                  iqb = 1+(iqz-1)*kx
                  iqe = iqz*kx
                  if ( (bnds(sproc)%send_list_uv(iq) > 0) .neqv. &
                        bnds(sproc)%send_swap(iq) ) then
                     bnds(sproc)%sbuf(iqb:iqe) = u(abs(bnds(sproc)%send_list_uv(iq)),:)
                  else
                     bnds(sproc)%sbuf(iqb:iqe) = v(abs(bnds(sproc)%send_list_uv(iq)),:)
                  end if 
                  if ( bnds(sproc)%send_neg(iq) ) bnds(sproc)%sbuf(iqb:iqe) = -bnds(sproc)%sbuf(iqb:iqe)
               end do
               iqq = iqq+ssplit(sproc)%iwufn-ssplit(sproc)%isvbg+1
            end if
            if ( fnveu .and. msk_ne==1 ) then
!cdir nodep
               do iq=ssplit(sproc)%invbg,ssplit(sproc)%ieufn
                  ! send_list_uv(iq) is point index.
                  ! Use abs because sign is used as u/v flag
                  iqz = iqq+iq-ssplit(sproc)%invbg+1
                  iqb = 1+(iqz-1)*kx
                  iqe = iqz*kx
                  if ( (bnds(sproc)%send_list_uv(iq) > 0) .neqv. &
                        bnds(sproc)%send_swap(iq) ) then
                     bnds(sproc)%sbuf(iqb:iqe) = u(abs(bnds(sproc)%send_list_uv(iq)),:)
                  else
                     bnds(sproc)%sbuf(iqb:iqe) = v(abs(bnds(sproc)%send_list_uv(iq)),:)
                  end if 
                  if ( bnds(sproc)%send_neg(iq) ) bnds(sproc)%sbuf(iqb:iqe) = -bnds(sproc)%sbuf(iqb:iqe)
               end do
               iqq = iqq+ssplit(sproc)%ieufn-ssplit(sproc)%invbg+1
            end if
            if ( fssvwwu .and. msk_sw==1 ) then
!cdir nodep
               do iq=ssplit(sproc)%issvbg,ssplit(sproc)%iwwufn
                  ! send_list_uv(iq) is point index.
                  ! Use abs because sign is used as u/v flag
                  iqz = iqq+iq-ssplit(sproc)%issvbg+1
                  iqb = 1+(iqz-1)*kx
                  iqe = iqz*kx
                  if ( (bnds(sproc)%send_list_uv(iq) > 0) .neqv. &
                        bnds(sproc)%send_swap(iq) ) then
                     bnds(sproc)%sbuf(iqb:iqe) = u(abs(bnds(sproc)%send_list_uv(iq)),:)
                  else
                     bnds(sproc)%sbuf(iqb:iqe) = v(abs(bnds(sproc)%send_list_uv(iq)),:)
                  end if 
                  if ( bnds(sproc)%send_neg(iq) ) bnds(sproc)%sbuf(iqb:iqe) = -bnds(sproc)%sbuf(iqb:iqe)
               end do
               iqq = iqq+ssplit(sproc)%iwwufn-ssplit(sproc)%issvbg+1
            end if
            if ( fnnveeu .and. msk_ne==1 ) then
!cdir nodep
               do iq=ssplit(sproc)%innvbg,ssplit(sproc)%ieeufn
                  ! send_list_uv(iq) is point index.
                  ! Use abs because sign is used as u/v flag
                  iqz = iqq+iq-ssplit(sproc)%innvbg+1
                  iqb = 1+(iqz-1)*kx
                  iqe = iqz*kx
                  if ( (bnds(sproc)%send_list_uv(iq) > 0) .neqv. &
                        bnds(sproc)%send_swap(iq) ) then
                     bnds(sproc)%sbuf(iqb:iqe) = u(abs(bnds(sproc)%send_list_uv(iq)),:)
                  else
                     bnds(sproc)%sbuf(iqb:iqe) = v(abs(bnds(sproc)%send_list_uv(iq)),:)
                  end if 
                  if ( bnds(sproc)%send_neg(iq) ) bnds(sproc)%sbuf(iqb:iqe) = -bnds(sproc)%sbuf(iqb:iqe)
               end do
               iqq = iqq+ssplit(sproc)%ieeufn-ssplit(sproc)%innvbg+1
            end if
            if ( fsuwv .and. msk_sw==1 ) then
!cdir nodep
               do iq=ssplit(sproc)%isubg,ssplit(sproc)%iwvfn
                  ! send_list_uv(iq) is point index.
                  ! Use abs because sign is used as u/v flag
                  iqz = iqq+iq-ssplit(sproc)%isubg+1
                  iqb = 1+(iqz-1)*kx
                  iqe = iqz*kx
                  if ( (bnds(sproc)%send_list_uv(iq) > 0) .neqv. &
                        bnds(sproc)%send_swap(iq) ) then
                     bnds(sproc)%sbuf(iqb:iqe) = u(abs(bnds(sproc)%send_list_uv(iq)),:)
                  else
                     bnds(sproc)%sbuf(iqb:iqe) = v(abs(bnds(sproc)%send_list_uv(iq)),:)
                  end if 
                  if ( bnds(sproc)%send_neg(iq) ) bnds(sproc)%sbuf(iqb:iqe) = -bnds(sproc)%sbuf(iqb:iqe)
               end do
               iqq = iqq+ssplit(sproc)%iwvfn-ssplit(sproc)%isubg+1
            end if
            if ( fnuev .and. msk_ne==1 ) then
!cdir nodep
               do iq=ssplit(sproc)%inubg,ssplit(sproc)%ievfn
                  ! send_list_uv(iq) is point index.
                  ! Use abs because sign is used as u/v flag
                  iqz = iqq+iq-ssplit(sproc)%inubg+1
                  iqb = 1+(iqz-1)*kx
                  iqe = iqz*kx
                  if ( (bnds(sproc)%send_list_uv(iq) > 0) .neqv. &
                        bnds(sproc)%send_swap(iq) ) then
                     bnds(sproc)%sbuf(iqb:iqe) = u(abs(bnds(sproc)%send_list_uv(iq)),:)
                  else
                     bnds(sproc)%sbuf(iqb:iqe) = v(abs(bnds(sproc)%send_list_uv(iq)),:)
                  end if 
                  if ( bnds(sproc)%send_neg(iq) ) bnds(sproc)%sbuf(iqb:iqe) = -bnds(sproc)%sbuf(iqb:iqe)
               end do
               iqq = iqq+ssplit(sproc)%ievfn-ssplit(sproc)%inubg+1
            end if
            if ( iqq > 0 ) then
               nreq = nreq + 1
               iqq=iqq*kx
               call MPI_ISend( bnds(sproc)%sbuf(1), iqq, &
                    ltype, sproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
            end if
         end if
      end do

      ! Finally see if there are any points on my own processor that need
      ! to be fixed up. This will only be in the case when nproc < npanels.
      if ( bnds(myid)%rlen_uv > 0 ) then
         if ( extra .or. stagmode == -19 ) then
            recv_len = bnds(myid)%rlenx_uv
         else if ( double .or. stagmode > 0 ) then
            recv_len = bnds(myid)%rlen2_uv
         else
            recv_len = bnds(myid)%rlen_uv
         end if
!cdir nodep
         do iq=1,recv_len
            ! request_list is same as send_list in this case
            iqb=1+(iq-1)*kx
            iqe=iq*kx
            if ( (bnds(myid)%request_list_uv(iq) > 0) .neqv. &
                     bnds(myid)%uv_swap(iq) ) then  ! haven't copied to send_swap yet
               tmp(iqb:iqe) = u(abs(bnds(myid)%request_list_uv(iq)),:)
            else
               tmp(iqb:iqe) = v(abs(bnds(myid)%request_list_uv(iq)),:)
            end if
            if ( bnds(myid)%uv_neg(iq) ) tmp(iqb:iqe) = -tmp(iqb:iqe)

            ! unpack_list(iq) is index into extended region
            if ( bnds(myid)%unpack_list_uv(iq) > 0 ) then
               u(ifull+bnds(myid)%unpack_list_uv(iq),:) = tmp(iqb:iqe)
            else
               v(ifull-bnds(myid)%unpack_list_uv(iq),:) = tmp(iqb:iqe)
            end if
         end do
      end if
      
      if ( nreq > 0 ) then
         call start_log(mpiwait_begin)
         call MPI_Waitall(nreq,ireq,status,ierr)
         call end_log(mpiwait_end)
      end if

      do iproc = 1,nproc-1  !
         rproc = modulo(myid-iproc,nproc)  ! Recv from
         if ( extra .or. stagmode == -19 ) then
            recv_len = bnds(rproc)%rlenx_uv
         else if ( double .or. stagmode > 0 ) then
            recv_len = bnds(rproc)%rlen2_uv
         else
            recv_len = bnds(rproc)%rlen_uv
         end if
         msk_ne=1
         msk_sw=1
         if (lmode==1) then
            recv_len=recv_len*bnds(rproc)%mlomsk
            msk_ne=mlomsk_ne
            msk_sw=mlomsk_sw
         end if
         iqq = 0
         if ( fsvwu .and. msk_sw==1 ) then
!cdir nodep
            do iq=rsplit(rproc)%isvbg,rsplit(rproc)%iwufn
               ! unpack_list(iq) is index into extended region
               iqz = iqq+iq-rsplit(rproc)%isvbg+1
               iqb = 1+(iqz-1)*kx
               iqe = iqz*kx
               if ( bnds(rproc)%unpack_list_uv(iq) > 0 ) then
                  u(ifull+bnds(rproc)%unpack_list_uv(iq),:) = bnds(rproc)%rbuf(iqb:iqe)
               else
                  v(ifull-bnds(rproc)%unpack_list_uv(iq),:) = bnds(rproc)%rbuf(iqb:iqe)
               end if
            end do
            iqq = iqq+rsplit(rproc)%iwufn-rsplit(rproc)%isvbg+1
         end if
         if ( fnveu .and. msk_ne==1 ) then
!cdir nodep
            do iq=rsplit(rproc)%invbg,rsplit(rproc)%ieufn
               ! unpack_list(iq) is index into extended region
               iqz = iqq+iq-rsplit(rproc)%invbg+1
               iqb = 1+(iqz-1)*kx
               iqe = iqz*kx
               if ( bnds(rproc)%unpack_list_uv(iq) > 0 ) then
                  u(ifull+bnds(rproc)%unpack_list_uv(iq),:) = bnds(rproc)%rbuf(iqb:iqe)
               else
                  v(ifull-bnds(rproc)%unpack_list_uv(iq),:) = bnds(rproc)%rbuf(iqb:iqe)
               end if
            end do
            iqq = iqq+rsplit(rproc)%ieufn-rsplit(rproc)%invbg+1
         end if         
         if ( fssvwwu .and. msk_sw==1 ) then
!cdir nodep
            do iq=rsplit(rproc)%issvbg,rsplit(rproc)%iwwufn
               ! unpack_list(iq) is index into extended region
               iqz = iqq+iq-rsplit(rproc)%issvbg+1
               iqb = 1+(iqz-1)*kx
               iqe = iqz*kx
               if ( bnds(rproc)%unpack_list_uv(iq) > 0 ) then
                  u(ifull+bnds(rproc)%unpack_list_uv(iq),:) = bnds(rproc)%rbuf(iqb:iqe)
               else
                  v(ifull-bnds(rproc)%unpack_list_uv(iq),:) = bnds(rproc)%rbuf(iqb:iqe)
               end if
            end do
            iqq = iqq+rsplit(rproc)%iwwufn-rsplit(rproc)%issvbg+1
         end if         
         if ( fnnveeu .and. msk_ne==1 ) then
!cdir nodep
            do iq=rsplit(rproc)%innvbg,rsplit(rproc)%ieeufn
               ! unpack_list(iq) is index into extended region
               iqz = iqq+iq-rsplit(rproc)%innvbg+1
               iqb = 1+(iqz-1)*kx
               iqe = iqz*kx
               if ( bnds(rproc)%unpack_list_uv(iq) > 0 ) then
                  u(ifull+bnds(rproc)%unpack_list_uv(iq),:) = bnds(rproc)%rbuf(iqb:iqe)
               else
                  v(ifull-bnds(rproc)%unpack_list_uv(iq),:) = bnds(rproc)%rbuf(iqb:iqe)
               end if
            end do
            iqq = iqq+rsplit(rproc)%ieeufn-rsplit(rproc)%innvbg+1
         end if         
         if ( fsuwv .and. msk_sw==1 ) then
!cdir nodep
            do iq=rsplit(rproc)%isubg,rsplit(rproc)%iwvfn
               ! unpack_list(iq) is index into extended region
               iqz = iqq+iq-rsplit(rproc)%isubg+1
               iqb = 1+(iqz-1)*kx
               iqe = iqz*kx
               if ( bnds(rproc)%unpack_list_uv(iq) > 0 ) then
                  u(ifull+bnds(rproc)%unpack_list_uv(iq),:) = bnds(rproc)%rbuf(iqb:iqe)
               else
                  v(ifull-bnds(rproc)%unpack_list_uv(iq),:) = bnds(rproc)%rbuf(iqb:iqe)
               end if
            end do
            iqq = iqq+rsplit(rproc)%iwvfn-rsplit(rproc)%isubg+1
         end if 
         if ( fnuev .and. msk_sw==1 ) then
!cdir nodep
            do iq=rsplit(rproc)%inubg,rsplit(rproc)%ievfn
               ! unpack_list(iq) is index into extended region
               iqz = iqq+iq-rsplit(rproc)%inubg+1
               iqb = 1+(iqz-1)*kx
               iqe = iqz*kx
               if ( bnds(rproc)%unpack_list_uv(iq) > 0 ) then
                  u(ifull+bnds(rproc)%unpack_list_uv(iq),:) = bnds(rproc)%rbuf(iqb:iqe)
               else
                  v(ifull-bnds(rproc)%unpack_list_uv(iq),:) = bnds(rproc)%rbuf(iqb:iqe)
               end if
            end do
            iqq = iqq+rsplit(rproc)%ievfn-rsplit(rproc)%inubg+1
         end if        
      end do

      call end_log(boundsuv_end)

   end subroutine boundsuv3

   subroutine deptsync(nface,xg,yg,gmode)
      ! Different levels will have different winds, so the list of points is
      ! different on each level.
      ! xg ranges from 0.5 to il+0.5 on a face. A given processors range
      ! is 0.5+ioff to 0.5+ipan+ioff
      ! Assignment of points to processors needs to match what ints does
      ! in case there's an exact edge point.
      ! Because of the boundary region, the range [0:ipan+1) can be handled.
      ! Need floor(xxg) in range [0:ipan]
      ! MJT - modified to restrict bounds calls to ocean processors with gmode=1
      use arrays_m
      integer, dimension(:,:), intent(in) :: nface
      integer, intent(in), optional :: gmode
      real, dimension(:,:), intent(in) :: xg, yg      
      integer :: iproc
      integer(kind=4) :: nreq, itag = 99, ierr, rproc, sproc, ltype
      integer(kind=4), dimension(MPI_STATUS_SIZE,2*nproc) :: status
      integer, dimension(0:nproc-1) :: binlen,msglen
      integer :: count, ip, jp, xn, kx
      integer :: iq, k, idel, jdel, nf
      integer :: lmode
      logical, dimension(0:nproc-1) :: lneighbour

      ! This does nothing in the one processor case
      if ( nproc == 1 ) return

      lmode=0
      if (present(gmode)) lmode=gmode
#ifdef r8i8
      ltype=MPI_DOUBLE_PRECISION
#else
      ltype=MPI_REAL
#endif   

      call start_log(deptsync_begin)
      dslen = 0
      drlen = 0
      do iproc=0,nproc-1
        if (neighbour(iproc)) then
          dindex(iproc)%a = 0
        end if
      end do
      kx=size(nface,2)
      select case(lmode)
         case(1)
            msglen=bnds(:)%len*bnds(:)%mlomsk
            lneighbour=neighbour.and.(bnds(:)%mlomsk==1)
         case default
            msglen=bnds(:)%len
            lneighbour=neighbour
      end select
      binlen=max(msglen,1)
      do k=1,kx
         do iq=1,ifull
            nf = nface(iq,k) + noff ! Make this a local index
            idel = int(xg(iq,k)) - ioff(nface(iq,k))
            jdel = int(yg(iq,k)) - joff(nface(iq,k))
            if ( idel < 0 .or. idel > ipan .or. jdel < 0 .or. jdel > jpan &
                 .or. nf < 1 .or. nf > npan ) then
               ! If point is on a different face, add to a list 
               ip = min(il_g,max(1,nint(xg(iq,k))))
               jp = min(il_g,max(1,nint(yg(iq,k))))
               iproc = fproc(ip,jp,nface(iq,k))
! This prevents vectorisation
#ifdef debug
               if ( iproc == myid ) then
                  write(6,*) "Inconsistency in deptsync"
                  call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
               end if
#endif
               ! Add this point to the list of requests I need to send to iproc
               dslen(iproc) = dslen(iproc) + 1
#ifdef debug
               call checksize(dslen(iproc),bnds(iproc)%len,"Deptssync")
#endif
               ! Since nface is a small integer it can be exactly represented by a
               ! real. It's simpler to send like this than use a proper structure.
               xn=min(dslen(iproc),binlen(iproc))
               dbuf(iproc)%a(:,xn) = (/ real(nface(iq,k)), xg(iq,k), yg(iq,k), real(k) /)
               dindex(iproc)%a(:,xn) = (/ iq,k /)
            end if
         end do
      end do
      
      ! Error check
      do iproc=0,nproc-1
        if (dslen(iproc)>msglen(iproc)) then
          write(6,*) "myid,iproc,neighbour,dslen,len ",myid,iproc,neighbour(iproc),dslen(iproc),msglen(iproc)
          iq=dindex(iproc)%a(1,1)
          k=dindex(iproc)%a(2,1)
          write(6,*) "Example error iq,k,u,v ",iq,k,u(iq,k),v(iq,k)
          call checksize(dslen(iproc),msglen(iproc),"Deptssync")
        end if
      end do

!     In this case the length of each buffer is unknown and will not
!     be symmetric between processors. Therefore need to get the length
!     from the message status
      nreq = 0
      do iproc = 1,nproc-1  !
         ! Is there any advantage to this ordering here or would send/recv
         ! to the same processor be just as good?
         rproc = modulo(myid-iproc,nproc)  ! Recv from
         if ( lneighbour(rproc) ) then
            nreq = nreq + 1
            ! Use the maximum size in the recv call.
            call MPI_IRecv( dpoints(rproc)%a, 4*bnds(rproc)%len, &
                         ltype, rproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do
      do iproc = 1,nproc-1  !
         ! Is there any advantage to this ordering here or would send/recv
         ! to the same processor be just as good?
         sproc = modulo(myid+iproc,nproc)  ! Send to
         if ( lneighbour(sproc) ) then
            ! Send, even if length is zero
            nreq = nreq + 1
            call MPI_ISend( dbuf(sproc)%a, 4*dslen(sproc), &
                    ltype, sproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         else
            if ( dslen(sproc) > 0 ) then
               write(6,*) "Error, dslen > 0 for non neighbour",      &
                    myid, sproc, dslen(sproc)
               call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
            end if
         end if
      end do
      if ( nreq > 0 ) then
         call start_log(mpiwait_begin)
         call MPI_Waitall(nreq,ireq,status,ierr)
         call end_log(mpiwait_end)
      end if

!     Now get the actual sizes from the status
      nreq = 0
      do iproc = 1,nproc-1  !
         rproc = modulo(myid-iproc,nproc)  ! Recv from
         if ( lneighbour(rproc) ) then
            nreq = nreq + 1
            call MPI_Get_count(status(1,nreq), MPI_REAL, count, ierr)
            drlen(rproc) = count/4
         end if
      end do

      call end_log(deptsync_end)

   end subroutine deptsync

   subroutine intssync(s)
      real, dimension(:,:), intent(inout) :: s
      integer :: iproc
      integer(kind=4) :: nreq, itag = 0, ierr, rproc, sproc, ltype
      integer :: iq
      integer(kind=4), dimension(MPI_STATUS_SIZE,2*nproc) :: status

      call start_log(intssync_begin)

#ifdef r8i8
      ltype=MPI_DOUBLE_PRECISION
#else
      ltype=MPI_REAL
#endif   

      ! When sending the results, roles of dslen and drlen are reversed
      nreq = 0
      do iproc = 1,nproc-1  !
         rproc = modulo(myid-iproc,nproc)  ! Recv from
         if ( dslen(rproc) /= 0 ) then
            nreq = nreq + 1
            call MPI_IRecv( dbuf(rproc)%b, dslen(rproc), &
                            ltype, rproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do
      do iproc = 1,nproc-1  !
         sproc = modulo(myid+iproc,nproc)  ! Send to
         if ( drlen(sproc) /= 0 ) then
            nreq = nreq + 1
            call MPI_ISend( sextra(sproc)%a, drlen(sproc), & 
                 ltype, sproc, itag, MPI_COMM_WORLD, ireq(nreq), ierr )
         end if
      end do
      
      if ( nreq > 0 ) then
         call start_log(mpiwait_begin)
         call MPI_Waitall(nreq,ireq,status,ierr)
         call end_log(mpiwait_end)
      end if

      do iproc=0,nproc-1
         if ( iproc == myid ) then
            cycle
         end if
         do iq=1,dslen(iproc)
            s(dindex(iproc)%a(1,iq),dindex(iproc)%a(2,iq)) = dbuf(iproc)%b(iq)
         end do
      end do

      call end_log(intssync_end)

   end subroutine intssync

   subroutine indv_mpi(iq, i, j, n)
      integer , intent(in) :: iq
      integer , intent(out) :: i
      integer , intent(out) :: j
      integer , intent(out) :: n
      integer :: ierr,ierr2

      ! Calculate local i, j, n from global iq

      ! Global i, j, n
      n = (iq - 1)/(il_g*il_g)
      j = 1 + (iq - n*il_g*il_g - 1)/il_g
      i = iq - (j - 1)*il_g - n*il_g*il_g
      if ( fproc(i,j,n) /= myid ) then
         write(*,"(a,5i5)") "Consistency failure in indv_mpi", myid, iq, i, j, n
         call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
      end if
      ! Reduced to values on my processor
      j = j - joff(n)
      i = i - ioff(n)
      n = n + noff      
   end subroutine indv_mpi

   function indglobal(i,j,n) result(iq)
      integer, intent(in) :: i, j, n
      integer :: iq

      ! Calculate a 1D global index from the global indices
      ! n in range 0:npanels
      iq = i + (j-1)*il_g + n*il_g*il_g
   end function indglobal

   function indg(i,j,n) result(iq)
      integer, intent(in) :: i, j, n
      integer :: iq

      ! Calculate a 1D global index from the local processors indices
      ! n in range 1..npan
      iq = i+ioff(n-noff) + (j+joff(n-noff)-1)*il_g + (n-noff)*il_g*il_g
   end function indg

   function indp(i,j,n) result(iq)
      integer, intent(in) :: i, j, n
      integer :: iq

      ! Calculate a 1D local index from the local processors indices
      ! Note that face number runs from 1 here.
      iq = i + (j-1)*ipan + (n-1)*ipan*jpan
   end function indp

   function iq2iqg(iq) result(iqg)
      integer, intent(in) :: iq
      integer :: iqg
      integer :: i, j, n

      ! Calculate global iqg from local iq

      ! MJT bug fix (should be ipan and jpan, not il and jl)
      ! Calculate local i, j, n
      n = 1 + (iq-1)/(ipan*jpan)  ! In range 1 .. npan
      j = 1 + ( iq - (n-1)*(ipan*jpan) - 1) / ipan
      i = iq - (j-1)*ipan - (n-1)*(ipan*jpan)
      iqg = indg(i,j,n)

   end function iq2iqg

   subroutine checksize(len, msize, mesg)
      integer, intent(in) :: len,msize
      character(len=*), intent(in) :: mesg
      integer(kind=4) :: ierr
      if ( len > msize ) then
         write(6,*) "Error, maxsize exceeded in ", mesg
         call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
      end if
   end subroutine checksize

   subroutine check_bnds_alloc(rproc, iext)
      integer, intent(in) :: rproc
      integer, intent(in) :: iext
      integer :: len
      integer(kind=4) ierr

!     Allocate the components of the bnds array. It's too much work to
!     get the exact sizes, so allocate a fixed size for each case where
!     there's an interaction.
      if ( bnds(rproc)%len == 0 ) then
         ! Not allocated yet.
         len = maxbuflen
         if (rproc /= myid) then
            allocate ( bnds(rproc)%rbuf(len) )
            allocate ( bnds(rproc)%sbuf(len) )
            allocate ( bnds(rproc)%request_list(len) )
            allocate ( bnds(rproc)%send_list(len) )
            allocate ( bnds(rproc)%unpack_list(len) )
         end if
         allocate ( bnds(rproc)%request_list_uv(len) )
         allocate ( bnds(rproc)%send_list_uv(len) )
         allocate ( bnds(rproc)%unpack_list_uv(len) )
         allocate ( bnds(rproc)%uv_swap(len), bnds(rproc)%send_swap(len) )
         allocate ( bnds(rproc)%uv_neg(len), bnds(rproc)%send_neg(len) )
         bnds(rproc)%uv_neg = .false.
         bnds(rproc)%send_neg = .false.
         bnds(rproc)%len = len
      else
         ! Just check length
         if ( max(kl,ol)*bnds(rproc)%rlen >=  bnds(rproc)%len ) then
            write(6,*) "Error, maximum length error in check_bnds_alloc"
            write(6,*) myid, rproc, bnds(rproc)%rlen,  bnds(rproc)%len, max(kl,ol)
            call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
         end if
         if ( iext >= iextra ) then
            write(6,*) "Error, iext maximum length error in check_bnds_alloc"
            write(6,*) myid, iext, iextra
            call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
         end if
      end if
   end subroutine check_bnds_alloc

   subroutine fix_index(iqq,larray,n,bnds,iext)
      integer, intent(in) :: iqq, n
      integer, dimension(:), intent(out) :: larray
      integer, intent(inout) :: iext
      type(bounds_info), dimension(0:), intent(inout) :: bnds
      integer :: rproc, iloc,jloc,nloc

      ! This processes extra corner points, so adds to rlenx
      ! Which processor has this point
      rproc = qproc(iqq)
      if ( rproc /= myid ) then ! Add to list
         call check_bnds_alloc(rproc, iext)
         bnds(rproc)%rlenx = bnds(rproc)%rlenx + 1
         bnds(rproc)%request_list(bnds(rproc)%rlenx) = iqq
         ! Increment extended region index
         iext = iext + 1
         bnds(rproc)%unpack_list(bnds(rproc)%rlenx) = iext
         larray(n) = ifull+iext
      else
         ! If it's on this processor, just use the local index
         call indv_mpi(iqq,iloc,jloc,nloc)
         larray(n) = indp(iloc,jloc,nloc)
      end if
   end subroutine fix_index

   subroutine fix_index2(iqq,larray,n,bnds,iext)
      integer, intent(in) :: iqq, n
      integer, dimension(:), intent(out) :: larray
      integer, intent(inout) :: iext
      type(bounds_info), dimension(0:), intent(inout) :: bnds
      integer :: rproc, iloc,jloc,nloc

      ! Which processor has this point
      rproc = qproc(iqq)
      if ( rproc /= myid ) then ! Add to list
         call check_bnds_alloc(rproc, iext)
         bnds(rproc)%rlen2 = bnds(rproc)%rlen2 + 1
         bnds(rproc)%request_list(bnds(rproc)%rlen2) = iqq
         ! Increment extended region index
         iext = iext + 1
         bnds(rproc)%unpack_list(bnds(rproc)%rlen2) = iext
         larray(n) = ifull+iext
      else
         ! If it's on this processor, just use the local index
         call indv_mpi(iqq,iloc,jloc,nloc)
         larray(n) = indp(iloc,jloc,nloc)
      end if
   end subroutine fix_index2

   subroutine proc_setup(npanels,ifull)
      include 'parm.h'
!     Routine to set up offsets etc.
      integer, intent(in) :: npanels, ifull
      integer :: i, j, n, iproc, nd, jdf, idjd_g

      call face_set(ipan,jpan,noff,ioff,joff,npan,il_g,myid,nproc,nxproc,nyproc)

      if ( nproc <= npanels+1 ) then

         do n=0,npanels
            fproc(:,:,n) = n/npan
            qproc(1+n*il_g**2:(n+1)*il_g**2) = n/npan
         end do

      else  ! nproc >= npanels+1

         iproc = 0
         qproc = -9999 ! Mask value so any points not set are obvious.
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
                  qproc(indglobal(i,j,n)) = fproc(i,j,n)
               end do
            end do
         end do

      end if

!      ipfull = ipan*jpan*npan
!      iextra = 4*npan*(ipan+jpan+8)

      ! Convert standard jd to a face index
      nd = (jd-1)/il_g ! 0: to match fproc
      jdf = jd - nd*il_g
      mydiag = ( myid == fproc(id,jdf,nd) )
      ! Convert global indices to ones on this processors region
      idjd_g = id + (jd-1)*il_g
      if ( mydiag ) then
         call indv_mpi(idjd_g,i,j,n)
         idjd = indp(i,j,n)
      else
         ! This should never be used so set a value that will give a bounds error
         idjd = huge(1)
      end if

   end subroutine proc_setup

   subroutine face_set(ipanx,jpanx,noffx,ioffx,joffx,npanx,il_gx,myidx,nprocx,nxprocx,nyprocx)
      integer, intent(in) :: myidx,nprocx,npanx,il_gx
      integer, intent(out) :: ipanx,jpanx,noffx,nxprocx,nyprocx
      integer, dimension(0:npanels), intent(out) :: ioffx,joffx 
      integer n
      integer(kind=4) ierr

      !  Processor allocation
      !  if  nprocx <= npanels+1, then each gets a number of full panels
      if ( nprocx <= npanels+1 ) then
         if ( modulo(npanels+1,nprocx) /= 0 ) then
            write(6,*) "Error, number of processors must divide number of panels"
            call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
         end if
!         npanx = (npanels+1)/nprocx
         ipanx = il_gx
         jpanx = il_gx
         noffx = 1 - myidx*npanx
         ioffx(:) = 0
         joffx(:) = 0
         nxprocx = 1
         nyprocx = 1
      else  ! nprocx >= npanels+1
         if ( modulo (nprocx, npanels+1) /= 0 ) then
            write(6,*) "Error, number of processors must be a multiple of number of panels"
            call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
         end if
!         npanx = 1
         n = nprocx / (npanels+1)
         !  n is the number of processors on each face
         !  Try to factor this into two values are close as possible.
         !  nxproc is the smaller of the 2.
         nxprocx = nint(sqrt(real(n)))
         nyprocx = n / nxprocx
         do nxprocx = nint(sqrt(real(n))), 1, -1
            nyprocx = n / nxprocx
            if ( modulo(il_gx,nxprocx) == 0 .and. modulo(il_gx,nyprocx) == 0 .and. &
                 nxprocx*nyprocx == n ) exit
         end do
         if ( nxprocx*nyprocx /= n ) then
            write(6,*) "Error in splitting up faces"
            call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
         end if

         ! Still need to check that the processor distribution is compatible
         ! with the grid.
         if ( modulo(il_gx,nxprocx) /= 0 ) then
            write(6,*) "Error, il not a multiple of nxproc", il_gx, nxprocx
            call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
         end if
         if ( modulo(il_gx,nyprocx) /= 0 ) then
            write(6,*) "Error, il not a multiple of nyproc", il_gx, nyprocx
            call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
         end if
         ipanx = il_gx/nxprocx
         jpanx = il_gx/nyprocx

         ! Set offsets for this processor
         call proc_region(myidx,0,ioffx(0),joffx(0),noffx,nxprocx,nyprocx,ipanx,jpanx,npanx,il_gx,nprocx,0)
         ioffx(1:npanels)=ioffx(0)
         joffx(1:npanels)=joffx(0)
      end if
   
   end subroutine face_set

   subroutine proc_setup_uniform(npanels,ifull)
      include 'parm.h'
!     Routine to set up offsets etc for the uniform decomposition
      integer, intent(in) :: npanels, ifull
      integer :: i, j, n, iproc, nd, jdf, idjd_g
      integer(kind=4) ierr

      call uniform_set(ipan,jpan,noff,ioff,joff,npan,il_g,myid,nproc,nxproc,nyproc)

      !iproc = 0
      !qproc = -9999 ! Mask value so any points not set are obvious.
      !do j=1,il_g,jpan
      !   do i=1,il_g,ipan
      !      fproc(i:i+ipan-1,j:j+jpan-1,:) = iproc
      !      iproc = iproc + 1
      !   end do
      !end do
      
      ! MJT suggested decomposition to improve load balance
      ! Here, processors assigned to different faces are
      ! (quasi-)equidistant.  For it to work correctly, we
      ! need to invert ipan and jpan on processors 3 to 5
      ! which wastes memory with x(il,jl) arrays.  Hence
      ! the processor assignement is only quasi-equidistant
      ! when nxproc.ne.nyproc
      iproc=0
      qproc=-9999
      do j=1,il_g,jpan
        do i=1,il_g,ipan
          if (j.le.il_g/2) then
            fproc(i:i+ipan-1,j:j+jpan-1,0)=iproc
          else
            fproc(il_g-i-ipan+2:il_g-i+1,j:j+jpan-1,0)=iproc
          end if
          iproc = iproc + 1
        end do
      end do
      iproc=0
      do j=1,il_g,jpan
        do i=1,il_g,ipan
          fproc(i:i+ipan-1,j:j+jpan-1,1)=iproc
          iproc = iproc + 1
        end do
      end do
      iproc=0
      do j=1,il_g,jpan
        do i=1,il_g,ipan
          if (i.le.il_g/2) then
            fproc(i:i+ipan-1,j:j+jpan-1,2)=iproc
          else
            fproc(i:i+ipan-1,il_g-j-jpan+2:il_g-j+1,2)=iproc
          end if
          iproc = iproc + 1
        end do
      end do
      iproc=0
      do i=1,il_g,ipan
        do j=1,il_g,jpan
          if (i.le.il_g/2) then
            fproc(i:i+ipan-1,j:j+jpan-1,3)=iproc
          else
            fproc(i:i+ipan-1,il_g-j-jpan+2:il_g-j+1,3)=iproc
          end if
          iproc = iproc + 1
        end do
      end do
      iproc=0
      do i=1,il_g,ipan
        do j=1,il_g,jpan
          fproc(i:i+ipan-1,j:j+jpan-1,4)=iproc
          iproc = iproc + 1
        end do
      end do
      iproc=0
      do i=1,il_g,ipan
        do j=1,il_g,jpan
          if (j.le.il_g/2) then
            fproc(i:i+ipan-1,j:j+jpan-1,5)=iproc
          else
            fproc(il_g-i-ipan+2:il_g-i+1,j:j+jpan-1,5)=iproc
          end if
          iproc = iproc + 1
        end do
      end do

      do n=0,npanels
         do j=1,il_g
            do i=1,il_g
               qproc(indglobal(i,j,n)) = fproc(i,j,n)
            end do
         end do
      end do
      
      if (any(qproc<0)) then
         write(6,*) "Error, incorrect assignment of processors to qproc"
         call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
      end if

!     Check that the values calculated here match those set as parameters
      if ( ipan /= il ) then
         write(6,*) "Error, parameter mismatch, ipan /= il", ipan, il
         call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
      end if
      if ( jpan*npan /= jl ) then
         write(6,*) "Error, parameter mismatch, jpan*npan /= jl", jpan, npan, jl
         call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
      end if

!      ipfull = ipan*jpan*npan
!      iextra = 4*npan*(ipan+jpan+8)

      ! Convert standard jd to a face index
      nd = (jd-1)/il_g ! 0: to match fproc
      jdf = jd - nd*il_g
      mydiag = ( myid == fproc(id,jdf,nd) )
      ! Convert global indices to ones on this processors region
      idjd_g = id + (jd-1)*il_g
      if ( mydiag ) then
         call indv_mpi(idjd_g,i,j,n)
         idjd = indp(i,j,n)
      else
         ! This should never be used so set a value that will give a bounds error
         idjd = huge(1)
      end if

   end subroutine proc_setup_uniform

   subroutine uniform_set(ipanx,jpanx,noffx,ioffx,joffx,npanx,il_gx,myidx,nprocx,nxprocx,nyprocx)
      integer, intent(in) :: myidx, nprocx, npanx, il_gx
      integer, intent(out) :: ipanx, jpanx, noffx, nxprocx, nyprocx
      integer, dimension(0:npanels), intent(out) :: ioffx, joffx 
      integer n
      integer(kind=4) ierr
      
      if ( npanx /= npanels+1 ) then
         write(6,*) "Error: inconsistency in proc_setup_uniform"
         stop
      end if
      !  Processor allocation: each processor gets a part of each panel
      !  Try to factor nproc into two values are close as possible.
      !  nxproc is the smaller of the 2.
      nxprocx = nint(sqrt(real(nprocx)))
      do nxprocx = nint(sqrt(real(nprocx))), 1, -1
         ! This will always exit eventually because it's trivially true 
         ! for nxproc=1
         nyprocx = nprocx / nxprocx
         if ( modulo(nprocx,nxprocx) == 0 .and. &
              modulo(il_gx,nxprocx) == 0  .and. &
              modulo(il_gx,nyprocx) == 0 ) exit
      end do
      nyprocx = nprocx / nxprocx
      if ( nxprocx*nyprocx /= nprocx ) then
         write(6,*) "Error in splitting up faces"
         call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
      end if

      ! Still need to check that the processor distribution is compatible
      ! with the grid.
      if ( modulo(il_gx,nxprocx) /= 0 ) then
         write(6,*) "Error, il not a multiple of nxproc", il_gx, nxprocx
         call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
      end if
      if ( modulo(il_gx,nyprocx) /= 0 ) then
         write(6,*) "Error, il not a multiple of nyproc", il_gx, nyprocx
         call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
      end if
      ipanx = il_gx/nxprocx
      jpanx = il_gx/nyprocx

      ! Set offsets for this processor
      do n=0,npanels
         call proc_region(myidx,n,ioffx(n),joffx(n),noffx,nxprocx,nyprocx,ipanx,jpanx,npanx,il_gx,nprocx,1)
      end do

   end subroutine uniform_set

   subroutine proc_region(procid,panid,ipoff,jpoff,npoff,nxprocx,nyprocx,ipanx,jpanx,npanx,il_gx,nprocx,dmode)
      ! Calculate the offsets for a given processor
      integer, intent(in) :: procid, panid, nxprocx, nyprocx, ipanx, jpanx, npanx, il_gx, nprocx, dmode
      integer, intent(out) :: ipoff, jpoff, npoff
      integer :: myface, mtmp
      integer(kind=4) ierr

      select case(dmode)
         case(0)
            if ( nprocx <= npanels+1 ) then
               npoff = 1 - procid*npanx
               ipoff = 0
               jpoff = 0
            else
               myface = procid / (nxprocx*nyprocx)
               npoff = 1 - myface
               ! mtmp is the processor index on this face, 0:(nxprox*nyproc-1)
               mtmp = procid - myface*(nxprocx*nyprocx)
               jpoff = (mtmp/nxprocx) * jpanx
               ipoff = modulo(mtmp,nxprocx) * ipanx
            end if
      
         case(1)
            ! Set offsets for this processor (same on all faces)
            !npoff = 1
            !jpoff = (procid/nxprocx) * jpanx
            !ipoff = modulo(procid,nxprocx)*ipanx
      
            ! MJT suggested decomposition to improve load balance
            npoff=1
            select case(panid)
               case(0)
                  jpoff = (procid/nxprocx) * jpanx
                  ipoff = modulo(procid,nxprocx) * ipanx
                  if (jpoff>=il_gx/2) then
                     ipoff=il_gx-ipoff-ipanx
                  end if
               case(1)
                  jpoff = (procid/nxprocx) * jpanx
                  ipoff = modulo(procid,nxprocx) * ipanx
               case(2)
                  jpoff = (procid/nxprocx) * jpanx
                  ipoff = modulo(procid,nxprocx) * ipanx
                  if (ipoff>=il_gx/2) then
                     jpoff=il_gx-jpoff-jpanx
                  end if
               case(3)
                  jpoff = modulo(procid,nyprocx) * jpanx
                  ipoff = (procid/nyprocx) * ipanx
                  if (ipoff>=il_gx/2) then
                     jpoff=il_gx-jpoff-jpanx
                  end if
               case(4)
                  jpoff = modulo(procid,nyprocx) * jpanx
                  ipoff = (procid/nyprocx) * ipanx
               case(5)
                  jpoff = modulo(procid,nyprocx) * jpanx
                  ipoff = (procid/nyprocx) * ipanx
                  if (jpoff>=il_gx/2) then
                     ipoff=il_gx-ipoff-ipanx
                  end if        
            end select
         case default
            write(6,*) "ERROR: Invalid decomposition ",dmode
            call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
      end select
     
   end subroutine proc_region

   subroutine start_log ( event )
      integer, intent(in) :: event
      integer :: ierr
#ifdef vampir
      include "VT.inc"
#endif
#ifdef mpilog
      ierr = MPE_log_event(event,0,"")
#endif
#ifdef vampir
      call vtenter(event, VT_NOSCL, ierr)
#endif
#ifdef simple_timer
#ifdef scyld
      double precision :: MPI_Wtime
#endif
      start_time(event) = MPI_Wtime()
#endif 
   end subroutine start_log

   subroutine end_log ( event )
      integer, intent(in) :: event
      integer :: ierr
#ifdef vampir
      include "VT.inc"
#endif
#ifdef mpilog
      ierr = MPE_log_event(event,0,"")
#endif
#ifdef vampir
      call vtleave(VT_NOSCL, ierr)
#endif
#ifdef simple_timer
#ifdef scyld
      double precision :: MPI_Wtime
#endif
      tot_time(event) = tot_time(event) + MPI_Wtime() - start_time(event)
#endif 
   end subroutine end_log

   subroutine log_off()
#ifdef vampir
      call vttraceoff()
#endif
   end subroutine log_off
   
   subroutine log_on()
#ifdef vampir
      call vttraceon()
#endif
   end subroutine log_on

   subroutine log_setup()
      integer :: ierr
      integer :: classhandle
#ifdef mpilog
      bounds_begin = MPE_Log_get_event_number()
      bounds_end = MPE_Log_get_event_number()
      ierr = MPE_Describe_state(bounds_begin, bounds_end, "Bounds", "yellow")
      boundsa_begin = MPE_Log_get_event_number()
      boundsa_end = MPE_Log_get_event_number()
      ierr = MPE_Describe_state(boundsa_begin, boundsa_end, "BoundsA", "DarkOrange1")
      boundsb_begin = MPE_Log_get_event_number()
      boundsb_end = MPE_Log_get_event_number()
      ierr = MPE_Describe_state(boundsb_begin, boundsb_end, "BoundsB", "DarkOrange1")
      boundsuv_begin = MPE_Log_get_event_number()
      boundsuv_end = MPE_Log_get_event_number()
      ierr = MPE_Describe_state(boundsuv_begin, boundsuv_end, "BoundsUV", "khaki")
      ints_begin = MPE_Log_get_event_number()
      ints_end = MPE_Log_get_event_number()
      ierr = MPE_Describe_state(ints_begin, ints_end, "Ints","goldenrod")
      nonlin_begin = MPE_Log_get_event_number()
      nonlin_end = MPE_Log_get_event_number()
      ierr = MPE_Describe_state(nonlin_begin, nonlin_end, "Nonlin", "IndianRed")
      helm_begin = MPE_Log_get_event_number()
      helm_end = MPE_Log_get_event_number()
      ierr = MPE_Describe_state(helm_begin, helm_end, "Helm", "magenta1")
      adjust_begin = MPE_Log_get_event_number()
      adjust_end = MPE_Log_get_event_number()
      ierr = MPE_Describe_state(adjust_begin, adjust_end, "Adjust", "blue")
      upglobal_begin = MPE_Log_get_event_number()
      upglobal_end = MPE_Log_get_event_number()
      ierr = MPE_Describe_state(upglobal_begin, upglobal_end, "Upglobal", "ForestGreen")
      hordifg_begin = MPE_Log_get_event_number()
      hordifg_end = MPE_Log_get_event_number()
      ierr = MPE_Describe_state(hordifg_begin, hordifg_end, "Hordifg", "blue")
      vadv_begin = MPE_Log_get_event_number()
      vadv_end = MPE_Log_get_event_number()
      ierr = MPE_Describe_state(vadv_begin, vadv_end, "Vadv", "blue")
      depts_begin = MPE_Log_get_event_number()
      depts_end = MPE_Log_get_event_number()
      ierr = MPE_Describe_state(depts_begin, depts_end, "Depts", "pink1")
      deptsync_begin = MPE_Log_get_event_number()
      deptsync_end = MPE_Log_get_event_number()
      ierr = MPE_Describe_state(deptsync_begin, deptsync_end, "Deptsync", "YellowGreen")
      intssync_begin = MPE_Log_get_event_number()
      intssync_end = MPE_Log_get_event_number()
      ierr = MPE_Describe_state(intssync_begin, intssync_end, "Intssync", "YellowGreen")
      stag_begin = MPE_Log_get_event_number()
      stag_end = MPE_Log_get_event_number()
      ierr = MPE_Describe_state(stag_begin, stag_end, "Stag", "YellowGreen")
      ocnstag_begin = MPE_Log_get_event_number()
      ocnstag_end = MPE_Log_get_event_number()
      ierr = MPE_Describe_state(ocnstag_begin, ocnstag_end, "Ocnstag", "YellowGreen")
      toij_begin = MPE_Log_get_event_number()
      toij_end = MPE_Log_get_event_number()
      ierr = MPE_Describe_state(toij_begin, toij_end, "Toij", "blue")
      physloadbal_begin = MPE_Log_get_event_number()
      physloadbal_end = MPE_Log_get_event_number()
      ierr = MPE_Describe_state(physloadbal_begin, physloadbal_end, "PhysLBbal", "blue")
      phys_begin = MPE_Log_get_event_number()
      phys_end = MPE_Log_get_event_number()
      ierr = MPE_Describe_state(phys_begin, phys_end, "Phys", "Yellow")
      outfile_begin = MPE_Log_get_event_number()
      outfile_end = MPE_Log_get_event_number()
      ierr = MPE_Describe_state(outfile_begin, outfile_end, "Outfile", "Yellow")
      onthefly_begin = MPE_Log_get_event_number()
      onthefly_end = MPE_Log_get_event_number()
      ierr = MPE_Describe_state(onthefly_begin, onthefly_end, "Onthefly", "Yellow")
      indata_begin = MPE_Log_get_event_number()
      indata_end = MPE_Log_get_event_number()
      ierr = MPE_Describe_state(indata_begin, indata_end, "Indata", "Yellow")
      nestin_begin = MPE_Log_get_event_number()
      nestin_end = MPE_Log_get_event_number()
      ierr = MPE_Describe_state(nestin_begin, nestin_end, "Nestin", "Yellow")
      gwdrag_begin = MPE_Log_get_event_number()
      gwdrag_end = MPE_Log_get_event_number()
      ierr = MPE_Describe_state(gwdrag_begin, gwdrag_end, "GWdrag", "Yellow")
      convection_begin = MPE_Log_get_event_number()
      convection_end = MPE_Log_get_event_number()
      ierr = MPE_Describe_state(convection_begin, convection_end, "Convection", "Yellow")
      cloud_begin = MPE_Log_get_event_number()
      cloud_end = MPE_Log_get_event_number()
      ierr = MPE_Describe_state(cloud_begin, cloud_end, "Cloud", "Yellow")
      radnet_begin = MPE_Log_get_event_number()
      radnet_end = MPE_Log_get_event_number()
      ierr = MPE_Describe_state(radnet_begin, radnet_end, "Rad_net", "Yellow")
      radmisc_begin = MPE_Log_get_event_number()
      radmisc_end = MPE_Log_get_event_number()
      ierr = MPE_Describe_state(radmisc_begin, radmisc_end, "Rad_misc", "Yellow")
      radsw_begin = MPE_Log_get_event_number()
      radsw_end = MPE_Log_get_event_number()
      ierr = MPE_Describe_state(radsw_begin, radsw_end, "Rad_SW", "Yellow")      
      radlw_begin = MPE_Log_get_event_number()
      radlw_end = MPE_Log_get_event_number()      
      ierr = MPE_Describe_state(radlw_begin, radlw_end, "Rad_LW", "Yellow")
      sfluxnet_begin = MPE_Log_get_event_number()
      sfluxnet_end = MPE_Log_get_event_number()
      ierr = MPE_Describe_state(sfluxnet_begin, sfluxnet_end, "Sflux_net", "Yellow")
      sfluxwater_begin = MPE_Log_get_event_number()
      sfluxwater_end = MPE_Log_get_event_number()
      ierr = MPE_Describe_state(sfluxwater_begin, sfluxwater_end, "Sflux_water", "Yellow")
      sfluxland_begin = MPE_Log_get_event_number()
      sfluxland_end = MPE_Log_get_event_number()
      ierr = MPE_Describe_state(sfluxland_begin, sfluxland_end, "Sflux_land", "Yellow")
      sfluxurban_begin = MPE_Log_get_event_number()
      sfluxurban_end = MPE_Log_get_event_number()
      ierr = MPE_Describe_state(sfluxurban_begin, sfluxurban_end, "Sflux_urban", "Yellow")
      vertmix_begin = MPE_Log_get_event_number()
      vertmix_end = MPE_Log_get_event_number()
      ierr = MPE_Describe_state(vertmix_begin, vertmix_end, "Vertmix", "Yellow")
      aerosol_begin = MPE_Log_get_event_number()
      aerosol_end = MPE_Log_get_event_number()
      ierr = MPE_Describe_state(aerosol_begin, aerosol_end, "Aerosol", "Yellow")
      waterdynamics_begin = MPE_Log_get_event_number()
      waterdynamics_end = MPE_Log_get_event_number()
      ierr = MPE_Describe_state(waterdynamics_begin, waterdynamics_end, "Waterdynamics", "blue")
      waterdiff_begin = MPE_Log_get_event_number()
      waterdiff_end = MPE_Log_get_event_number()
      ierr = MPE_Describe_state(waterdiff_begin, waterdiff_end, "Waterdiff", "blue")
      river_begin = MPE_Log_get_event_number()
      river_end = MPE_Log_get_event_number()
      ierr = MPE_Describe_state(river_begin, river_end, "River", "Yellow")
#endif
#ifdef vampir
      call vtfuncdef("Bounds", classhandle, bounds_begin, ierr)
      bounds_end = bounds_begin
      call vtfuncdef("BoundsA", classhandle, boundsa_begin, ierr)
      boundsa_end = boundsa_end
      call vtfuncdef("BoundsB", classhandle, boundsb_begin, ierr)
      boundsb_end = boundsb_begin
      call vtfuncdef("BoundsUV", classhandle, boundsuv_begin, ierr)
      boundsuv_end = boundsuv_begin
      call vtfuncdef("Ints", classhandle, ints_begin, ierr)
      ints_end = ints_begin 
      call vtfuncdef("Nonlin", classhandle, nonlin_begin, ierr)
      nonlin_end = nonlin_begin 
      call vtfuncdef("Helm", classhandle, helm_begin, ierr)
      helm_end = helm_begin
      call vtfuncdef("Adjust", classhandle, adjust_begin, ierr)
      adjust_end = adjust_begin
      call vtfuncdef("Upglobal", classhandle, upglobal_begin, ierr)
      upglobal_end = upglobal_begin
      call vtfuncdef("Hordifg", classhandle, hordifg_begin, ierr)
      hordifg_end = hordifg_begin
      call vtfuncdef("Vadv", classhandle, vadv_begin, ierr)
      vadv_end = vadv_begin
      call vtfuncdef("Depts", classhandle, depts_begin, ierr)
      depts_end = depts_begin
      call vtfuncdef("Deptsync", classhandle, deptsync_begin, ierr)
      deptsync_end = deptsync_begin
      call vtfuncdef("Intssync", classhandle, intssync_begin, ierr)
      intssync_end = intssync_begin
      call vtfuncdef("Stag", classhandle, stag_begin, ierr)
      stag_end = stag_begin
      call vtfuncdef("Ocnstag", classhandle, ocnstag_begin, ierr)
      ocnstag_end = ocnstag_begin
      call vtfuncdef("Toij", classhandle, toij_begin, ierr)
      toij_end =  toij_begin
      call vtfuncdef("PhysLBal", classhandle, physloadbal_begin, ierr)
      physloadbal_end =  physloadbal_begin
      call vtfuncdef("Phys", classhandle, phys_begin, ierr)
      phys_end =  physloadbal_begin
      call vtfuncdef("Outfile", classhandle, outfile_begin, ierr)
      outfile_end =  outfile_begin
      call vtfuncdef("Onthefly", classhandle, onthefly_begin, ierr)
      onthefly_end =  onthefly_begin
      call vtfuncdef("Indata", classhandle, indata_begin, ierr)
      indata_end =  indata_begin
      call vtfuncdef("Nestin", classhandle, nestin_begin, ierr)
      nestin_end =  nestin_begin
      call vtfuncdef("GWdrag", classhandle, gwdrag_begin, ierr)
      gwdrag_end =  gwdrag_begin
      call vtfuncdef("Convection", classhandle, convection_begin, ierr)
      convection_end =  convection_begin
      call vtfuncdef("Cloud", classhandle, cloud_begin, ierr)
      cloud_end =  cloud_begin
      call vtfuncdef("Rad_net", classhandle, radnet_begin, ierr)
      radnet_end =  radnet_begin
      call vtfuncdef("Rad_misc", classhandle, radmisc_begin, ierr)
      radmisc_end =  radmisc_begin
      call vtfuncdef("Rad_SW", classhandle, radsw_begin, ierr)
      radsw_end =  radsw_begin
      call vtfuncdef("Rad_LW", classhandle, radlw_begin, ierr)
      radlw_end =  radlw_begin      
      call vtfuncdef("Sflux_net", classhandle, sfluxnet_begin, ierr)
      sfluxnet_end =  sfluxnet_begin
      call vtfuncdef("Sflux_water", classhandle, sfluxwater_begin, ierr)
      sfluxwater_end =  sfluxwater_begin
      call vtfuncdef("Sflux_land", classhandle, sfluxland_begin, ierr)
      sfluxland_end =  sfluxland_begin
      call vtfuncdef("Sflux_urban", classhandle, sfluxurban_begin, ierr)
      sfluxurban_end =  sfluxurban_begin
      call vtfuncdef("Vertmix", classhandle, vertmix_begin, ierr)
      vertmix_end =  vertmix_begin
      call vtfuncdef("Aerosol", classhandle, aerosol_begin, ierr)
      aerosol_end = aerosol_begin
      call vtfuncdef("Waterdynamics", classhandle, waterdynamics_begin, ierr)
      waterdynamics_end =  waterdynamics_begin
      call vtfuncdef("Waterdiff", classhandle, waterdiff_begin, ierr)
      waterdiff_end =  waterdiff_begin
      call vtfuncdef("River", classhandle, river_begin, ierr)
      river_end =  river_begin
#endif
#ifdef simple_timer

      model_begin = 1
      model_end =  model_begin
      event_name(model_begin) = "Whole Model"

      maincalc_begin = 2
      maincalc_end =  maincalc_begin
      event_name(maincalc_begin) = "Main Calc loop"

      phys_begin = 3
      phys_end =  phys_begin
      event_name(phys_begin) = "Phys"

      physloadbal_begin = 4
      physloadbal_end =  physloadbal_begin
      event_name(physloadbal_begin) = "Phys Load Bal"

      ints_begin = 5
      ints_end = ints_begin 
      event_name(ints_begin) = "Ints"

      nonlin_begin = 6
      nonlin_end = nonlin_begin 
      event_name(nonlin_begin) = "Nonlin"

      helm_begin = 7
      helm_end = helm_begin
      event_name(helm_begin) = "Helm"

      adjust_begin = 8
      adjust_end = adjust_begin
      event_name(Adjust_begin) = "Adjust"

      upglobal_begin = 9
      upglobal_end = upglobal_begin
      event_name(upglobal_begin) = "Upglobal"

      hordifg_begin = 10
      hordifg_end = hordifg_begin
      event_name(hordifg_begin) = "Hordifg"

      vadv_begin = 11
      vadv_end = vadv_begin
      event_name(vadv_begin) = "Vadv"

      depts_begin = 12
      depts_end = depts_begin
      event_name(depts_begin) = "Depts"

      stag_begin = 13
      stag_end = stag_begin
      event_name(stag_begin) = "Stag"

      ocnstag_begin = 14
      ocnstag_end = ocnstag_begin
      event_name(ocnstag_begin) = "Ocnstag"

      toij_begin = 15
      toij_end =  toij_begin
      event_name(toij_begin) = "Toij"

      outfile_begin = 16
      outfile_end =  outfile_begin
      event_name(outfile_begin) = "Outfile"

      onthefly_begin = 17
      onthefly_end =  onthefly_begin
      event_name(onthefly_begin) = "Onthefly"

      bounds_begin = 18
      bounds_end = bounds_begin
      event_name(bounds_begin) = "Bounds"

      boundsa_begin = 19
      boundsa_end = boundsa_end
      event_name(boundsa_begin) = "BoundsA"

      boundsb_begin = 20
      boundsb_end = boundsb_begin
      event_name(boundsb_begin) = "BoundsB"

      boundsuv_begin = 21
      boundsuv_end = boundsuv_begin
      event_name(boundsuv_begin) = "BoundsUV"

      deptsync_begin = 22
      deptsync_end = deptsync_begin
      event_name(deptsync_begin) = "Deptsync"

      intssync_begin = 23
      intssync_end = intssync_begin
      event_name(intssync_begin) = "Intssync"

      gather_begin = 24
      gather_end = gather_begin
      event_name(gather_begin) = "Gather"

      distribute_begin = 25
      distribute_end = distribute_begin
      event_name(distribute_begin) = "Distribute"

      posneg_begin = 26
      posneg_end = posneg_begin
      event_name(posneg_begin) = "Posneg"

      globsum_begin = 27
      globsum_end = globsum_begin
      event_name(globsum_begin) = "Globsum"

      reduce_begin = 28
      reduce_end = reduce_begin
      event_name(reduce_begin) = "Reduce"

      bcast_begin = 29
      bcast_end = bcast_begin
      event_name(bcast_begin) = "Bcast"

      mpiwait_begin = 30
      mpiwait_end = mpiwait_begin
      event_name(mpiwait_begin) = "MPI_Wait"

      precon_begin = 31
      precon_end = precon_begin
      event_name(precon_begin) = "Precon"

      indata_begin = 32
      indata_end =  indata_begin
      event_name(indata_begin) = "Indata"

      nestin_begin = 33
      nestin_end =  nestin_begin
      event_name(nestin_begin) = "Nestin"
      
      gwdrag_begin = 34
      gwdrag_end =  gwdrag_begin
      event_name(gwdrag_begin) = "GWdrag"

      convection_begin = 35
      convection_end =  convection_begin
      event_name(convection_begin) = "Convection"

      cloud_begin = 36
      cloud_end =  cloud_begin
      event_name(cloud_begin) = "Cloud"

      radnet_begin = 37
      radnet_end =  radnet_begin
      event_name(radnet_begin) = "Rad_net"

      radmisc_begin = 38
      radmisc_end =  radmisc_begin
      event_name(radmisc_begin) = "Rad_misc"
      
      radsw_begin = 39
      radsw_end =  radsw_begin
      event_name(radsw_begin) = "Rad_SW"

      radlw_begin = 40
      radlw_end =  radlw_begin
      event_name(radlw_begin) = "Rad_LW"      

      sfluxnet_begin = 41
      sfluxnet_end =  sfluxnet_begin
      event_name(sfluxnet_begin) = "Sflux_net"
      
      sfluxwater_begin = 42
      sfluxwater_end =  sfluxwater_begin
      event_name(sfluxwater_begin) = "Sflux_water"

      sfluxland_begin = 43
      sfluxland_end =  sfluxland_begin
      event_name(sfluxland_begin) = "Sflux_land"

      sfluxurban_begin = 44
      sfluxurban_end =  sfluxurban_begin
      event_name(sfluxurban_begin) = "Sflux_urban"

      vertmix_begin = 45
      vertmix_end =  vertmix_begin
      event_name(vertmix_begin) = "Vertmix"

      aerosol_begin = 46
      aerosol_end =  aerosol_begin
      event_name(aerosol_begin) = "Aerosol"

      waterdynamics_begin = 47
      waterdynamics_end =  waterdynamics_begin
      event_name(waterdynamics_begin) = "Waterdynamics"

      waterdiff_begin = 48
      waterdiff_end =  waterdiff_begin
      event_name(waterdiff_begin) = "Waterdiff"

      river_begin = 49
      river_end =  river_begin
      event_name(river_begin) = "River"

#endif
   end subroutine log_setup
   
   subroutine phys_loadbal()
!     This forces a sychronisation to make the physics load imbalance overhead
!     explicit. 
      integer(kind=4) :: ierr
      call start_log(physloadbal_begin)
      call MPI_Barrier( MPI_COMM_WORLD, ierr )
      call end_log(physloadbal_end)
   end subroutine phys_loadbal

#ifdef simple_timer
      subroutine simple_timer_finalize()
         ! Calculate the mean, min and max times for each case
         integer :: i
         integer(kind=4) ierr
         double precision, dimension(nevents) :: emean, emax, emin
         call MPI_Reduce(tot_time, emean, nevents, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, 0, MPI_COMM_WORLD, ierr )
         call MPI_Reduce(tot_time, emax, nevents, MPI_DOUBLE_PRECISION, &
                         MPI_MAX, 0, MPI_COMM_WORLD, ierr )
         call MPI_Reduce(tot_time, emin, nevents, MPI_DOUBLE_PRECISION, &
                         MPI_MIN, 0, MPI_COMM_WORLD, ierr )
         if ( myid == 0 ) then
            write(6,*) "==============================================="
            write(6,*) "  Times over all processes"
            write(6,*) "  Routine        Mean time  Min time  Max time"
            do i=1,nevents
               if ( emean(i) > 0. ) then
                  ! This stops boundsa, b getting written when they're not used.
                  write(*,"(a,3f10.3)") event_name(i), emean(i)/nproc, emin(i), emax(i)
               end if
            end do
         end if
      end subroutine simple_timer_finalize
#endif

    subroutine ccglobal_posneg2 (array, delpos, delneg, comm)
       ! Calculate global sums of positive and negative values of array
       ! MJT - modified to restrict calls to comm for ocean processors
       use sumdd_m
       use xyzinfo_m       
       include 'newmpar.h'

       real, intent(in), dimension(ifull) :: array
       real, intent(out) :: delpos, delneg
       real :: delpos_l, delneg_l
       real, dimension(2) :: delarr, delarr_l
       integer, intent(in), optional :: comm
       integer :: iq
       integer(kind=4) :: lcomm, ierr, ltype
       complex, dimension(2) :: local_sum, global_sum
!      Temporary array for the drpdr_local function
       real, dimension(ifull) :: tmparr, tmparr2 

       call start_log(posneg_begin)

       lcomm=MPI_COMM_WORLD
       if (present(comm)) lcomm=comm

       delpos_l = 0.
       delneg_l = 0.
       do iq=1,ifull
          tmparr(iq)  = max(0.,array(iq)*wts(iq))
          tmparr2(iq) = min(0.,array(iq)*wts(iq))
       enddo
       local_sum = (0.,0.)
       call drpdr_local(tmparr, local_sum(1))
       call drpdr_local(tmparr2, local_sum(2))
       delpos_l=real(local_sum(1))
       delneg_l=real(local_sum(2))
#ifdef sumdd
#ifdef r8i8
      ltype=MPI_COMPLEX8
#else
      ltype=MPI_COMPLEX
#endif   
       call MPI_Allreduce ( local_sum, global_sum, 2, ltype,     &
                            MPI_SUMDR, lcomm, ierr )
       delpos = real(global_sum(1))
       delneg = real(global_sum(2))
#else
#ifdef r8i8
      ltype=MPI_DOUBLE_PRECISION
#else
      ltype=MPI_REAL
#endif   
       delarr_l(1:2) = (/ delpos_l, delneg_l /)
       call MPI_Allreduce ( delarr_l, delarr, 2, ltype, MPI_SUM,    &
                            lcomm, ierr )
       delpos = delarr(1)
       delneg = delarr(2)
#endif

       call end_log(posneg_end)

    end subroutine ccglobal_posneg2
    
    subroutine ccglobal_posneg3 (array, delpos, delneg, dsigin, comm)
       ! Calculate global sums of positive and negative values of array
       ! MJT - modified to restrict calls to comm for ocean processors
       use sigs_m
       use sumdd_m
       use xyzinfo_m
       include 'newmpar.h'
       real, intent(in), dimension(:,:) :: array
       real, intent(in), dimension(:), optional :: dsigin
       real, intent(out) :: delpos, delneg
       integer, intent(in), optional :: comm
       real :: delpos_l, delneg_l
       real, dimension(size(array,2)) :: dsigx
       real, dimension(2) :: delarr, delarr_l
       integer :: k, iq, kx
       integer(kind=4) ierr, lcomm, ltype
       complex, dimension(2) :: local_sum, global_sum
!      Temporary array for the drpdr_local function
       real, dimension(ifull) :: tmparr, tmparr2 

       call start_log(posneg_begin)

       lcomm=MPI_COMM_WORLD
       if (present(comm)) lcomm=comm

       delpos_l = 0.
       delneg_l = 0.
       kx=size(array,2)
       if (present(dsigin)) then
         dsigx=-dsigin
       else
         dsigx=dsig
       end if
       local_sum = (0.,0.)
       do k=1,kx
          do iq=1,ifull
             tmparr(iq)  = max(0.,-dsigx(k)*array(iq,k)*wts(iq))
             tmparr2(iq) = min(0.,-dsigx(k)*array(iq,k)*wts(iq))
          end do
          call drpdr_local(tmparr, local_sum(1))
          call drpdr_local(tmparr2, local_sum(2))
          delpos_l=real(local_sum(1))
          delneg_l=real(local_sum(2))
       end do ! k loop
#ifdef sumdd
#ifdef r8i8
      ltype=MPI_COMPLEX8
#else
      ltype=MPI_COMPLEX
#endif   
       call MPI_Allreduce ( local_sum, global_sum, 2, ltype,     &
                            MPI_SUMDR, lcomm, ierr )
       delpos = real(global_sum(1))
       delneg = real(global_sum(2))
#else
#ifdef r8i8
      ltype=MPI_DOUBLE_PRECISION
#else
      ltype=MPI_REAL
#endif   
       delarr_l(1:2) = (/ delpos_l, delneg_l /)
       call MPI_Allreduce ( delarr_l, delarr, 2,ltype, MPI_SUM,    &
                            lcomm, ierr )
       delpos = delarr(1)
       delneg = delarr(2)
#endif

       call end_log(posneg_end)

    end subroutine ccglobal_posneg3

    subroutine ccglobal_sum2 (array, result)
       ! Calculate global sum of an array
       use sumdd_m
       use xyzinfo_m
       include 'newmpar.h'
       real, intent(in), dimension(ifull) :: array
       real, intent(out) :: result
       real :: result_l
       integer :: iq
       integer(kind=4) ierr
#ifdef sumdd
       complex :: local_sum, global_sum
!      Temporary array for the drpdr_local function
       real, dimension(ifull) :: tmparr
#endif

       call start_log(globsum_begin)

       result_l = 0.
       do iq=1,ifull
#ifdef sumdd         
          tmparr(iq)  = array(iq)*wts(iq)
#else
          result_l = result_l + array(iq)*wts(iq)
#endif
       enddo
#ifdef sumdd
       local_sum = (0.,0.)
       call drpdr_local(tmparr, local_sum)
       call MPI_Allreduce ( local_sum, global_sum, 1, MPI_COMPLEX,     &
                            MPI_SUMDR, MPI_COMM_WORLD, ierr )
       result = real(global_sum)
#else
       call MPI_Allreduce ( result_l, result, 1, MPI_REAL, MPI_SUM,    &
                            MPI_COMM_WORLD, ierr )
#endif

       call end_log(globsum_end)

    end subroutine ccglobal_sum2

    subroutine ccglobal_sum3 (array, result)
       ! Calculate global sum of 3D array, appyling vertical weighting
       use sigs_m
       use sumdd_m
       use xyzinfo_m
       include 'newmpar.h'
       real, intent(in), dimension(ifull,kl) :: array
       real, intent(out) :: result
       real :: result_l
       integer :: k, iq
       integer(kind=4) ierr
#ifdef sumdd
       complex :: local_sum, global_sum
!      Temporary array for the drpdr_local function
       real, dimension(ifull) :: tmparr
#endif

       call start_log(globsum_begin)

       result_l = 0.
#ifdef sumdd
       local_sum = (0.,0.)
#endif
       do k=1,kl
          do iq=1,ifull
#ifdef sumdd         
             tmparr(iq)  = -dsig(k)*array(iq,k)*wts(iq)
#else
             result_l = result_l - dsig(k)*array(iq,k)*wts(iq)
#endif
          enddo
#ifdef sumdd
          call drpdr_local(tmparr, local_sum)
#endif
       end do ! k

#ifdef sumdd
       call MPI_Allreduce ( local_sum, global_sum, 1, MPI_COMPLEX,    &
                            MPI_SUMDR, MPI_COMM_WORLD, ierr )
       result = real(global_sum)
#else
       call MPI_Allreduce ( result_l, result, 1, MPI_REAL, MPI_SUM,   &
                            MPI_COMM_WORLD, ierr )
#endif

       call end_log(globsum_end)

    end subroutine ccglobal_sum3

    ! Read and distribute a global variable
    ! Optional arguments for format and to skip over records
    subroutine readglobvar2(un,var,skip,fmt)
      include 'newmpar.h'
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
      include 'newmpar.h'
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
      include 'newmpar.h'
      integer, intent(in) :: un
      real, dimension(:,:), intent(out) :: var
      logical, intent(in), optional :: skip
      real, dimension(ifull_g,size(var,2)) :: varg
      character(len=*), intent(in), optional :: fmt
      integer :: k, kk
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
            ! Use explicit ranges here because some arguments might be extended.
            ! ccmpi_distribute3 expects kl, it's not general
            if ( kk == kl ) then
               call ccmpi_distribute(var(1:ifull,:),varg)
            else
               do k=1,kk
                  call ccmpi_distribute(var(1:ifull,k),varg(:,k))
               end do
            end if
         else
            if ( kk == kl ) then
               call ccmpi_distribute(var(1:ifull,:))
            else
               do k=1,kk
                  call ccmpi_distribute(var(1:ifull,k))
               end do
            end if
         end if
      end if

   end subroutine readglobvar3

    ! Gather and write a global variable
    ! Optional argument for format
    subroutine writeglobvar2(un,var,fmt)
      include 'newmpar.h'
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
      include 'newmpar.h'
      integer, intent(in) :: un
      real, dimension(:,:), intent(in) :: var
      real, dimension(ifull_g,size(var,2)) :: varg
      character(len=*), intent(in), optional :: fmt
      integer :: k, kk

      kk = size(var,2)

      if ( myid == 0 ) then
         ! Use explicit ranges here because some arguments might be extended.
         ! ccmpi_gather3 expects kl, it's not general
         if ( kk == kl ) then
            call ccmpi_gather(var(1:ifull,:),varg)
         else
            do k=1,kk
               call ccmpi_gather(var(1:ifull,k),varg(:,k))
            end do
         end if
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
         if ( kk == kl ) then
            call ccmpi_gather(var(1:ifull,:))
         else
            do k=1,kk
               call ccmpi_gather(var(1:ifull,k))
            end do
         end if
      end if

   end subroutine writeglobvar3

   subroutine ccmpi_reduce2i(ldat,gdat,op,host,comm)
   
      integer, intent(in) :: host,comm
      integer(kind=4) ltype,lop,lcomm,ierr,lsize,lkind,lhost
      integer, dimension(:), intent(in) :: ldat
      integer, dimension(:), intent(out) :: gdat
      character(len=*), intent(in) :: op
      
      call start_log(reduce_begin)
      
      select case(op)
         case("max")
            lop=MPI_MAX
         case("min")
            lop=MPI_MIN
         case("sum")
            lop=MPI_SUM
         case default
            write(6,*) "ERROR: Unknown option for ccmpi_reduce ",op
            call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
      end select
      
      lhost=host
      lcomm=comm
      lsize=size(ldat)
#ifdef r8i8
      ltype=MPI_INTEGER8
#else
      ltype=MPI_INTEGER
#endif 

      call MPI_Reduce(ldat, gdat, lsize, ltype, lop, lhost, lcomm, ierr )
 
      call end_log(reduce_end)
   
   end subroutine ccmpi_reduce2i

   subroutine ccmpi_reduce2r(ldat,gdat,op,host,comm)
   
      integer, intent(in) :: host,comm
      integer(kind=4) ltype,lop,lcomm,lerr,lsize,lhost
      real, dimension(:), intent(in) :: ldat
      real, dimension(:), intent(out) :: gdat
      character(len=*), intent(in) :: op
      
      call start_log(reduce_begin)
      
      lhost=host
      lcomm=comm
      lsize=size(ldat)
            
      select case(op)
         case("max")
            lop=MPI_MAX
#ifdef r8i8
            ltype=MPI_DOUBLE_PRECISION
#else
            ltype=MPI_REAL
#endif 
         case("min")
            lop=MPI_MIN
#ifdef r8i8
            ltype=MPI_DOUBLE_PRECISION
#else
            ltype=MPI_REAL
#endif 
         case("sum")
            lop=MPI_SUM
#ifdef r8i8
            ltype=MPI_DOUBLE_PRECISION
#else
            ltype=MPI_REAL
#endif 
         case("maxloc")
            lop=MPI_MAXLOC
            lsize=lsize/2
#ifdef r8i8
            ltype=MPI_2DOUBLE_PRECISION
#else
            ltype=MPI_2REAL
#endif 
         case("minloc")
            lop=MPI_MINLOC
            lsize=lsize/2
#ifdef r8i8
            ltype=MPI_2DOUBLE_PRECISION
#else
            ltype=MPI_2REAL
#endif 
         case default
            write(6,*) "ERROR: Unknown option for ccmpi_reduce ",op
            call MPI_Abort(MPI_COMM_WORLD,-1,lerr)
      end select
     
      call MPI_Reduce(ldat, gdat, lsize, ltype, lop, lhost, lcomm, lerr )
   
      call end_log(reduce_end)
   
   end subroutine ccmpi_reduce2r

   subroutine ccmpi_reduce3r(ldat,gdat,op,host,comm)
   
      integer, intent(in) :: host,comm
      integer(kind=4) ltype,lop,lcomm,lerr,lsize,lhost
      real, dimension(:,:), intent(in) :: ldat
      real, dimension(:,:), intent(out) :: gdat
      character(len=*), intent(in) :: op
      
      call start_log(reduce_begin)

      lhost=host
      lcomm=comm
      lsize=size(ldat)
            
      select case(op)
         case("max")
            lop=MPI_MAX
#ifdef r8i8
            ltype=MPI_DOUBLE_PRECISION
#else
            ltype=MPI_REAL
#endif 
         case("min")
            lop=MPI_MIN
#ifdef r8i8
            ltype=MPI_DOUBLE_PRECISION
#else
            ltype=MPI_REAL
#endif 
         case("sum")
            lop=MPI_SUM
#ifdef r8i8
            ltype=MPI_DOUBLE_PRECISION
#else
            ltype=MPI_REAL
#endif 
         case("maxloc")
            lop=MPI_MAXLOC
            lsize=lsize/2
#ifdef r8i8
            ltype=MPI_2DOUBLE_PRECISION
#else
            ltype=MPI_2REAL
#endif 
         case("minloc")
            lop=MPI_MINLOC
            lsize=lsize/2
#ifdef r8i8
            ltype=MPI_2DOUBLE_PRECISION
#else
            ltype=MPI_2REAL
#endif 
         case default
            write(6,*) "ERROR: Unknown option for ccmpi_reduce ",op
            call MPI_Abort(MPI_COMM_WORLD,-1,lerr)
      end select
      
      call MPI_Reduce(ldat, gdat, lsize, ltype, lop, lhost, lcomm, lerr )
   
      call end_log(reduce_end)
   
   end subroutine ccmpi_reduce3r

   subroutine ccmpi_reduce2c(ldat,gdat,op,host,comm)
   
      use sumdd_m
   
      integer, intent(in) :: host,comm
      integer(kind=4) ltype,lop,lcomm,lerr,lsize,lhost
      complex, dimension(:), intent(in) :: ldat
      complex, dimension(:), intent(out) :: gdat
      character(len=*), intent(in) :: op
      
      call start_log(reduce_begin)
      
      select case(op)
         case("max")
            lop=MPI_MAX
         case("min")
            lop=MPI_MIN
         case("sum")
            lop=MPI_SUM
         case("sumdr")
            lop=MPI_SUMDR
         case default
            write(6,*) "ERROR: Unknown option for ccmpi_reduce ",op
            call MPI_Abort(MPI_COMM_WORLD,-1,lerr)
      end select
      
      lhost=host
      lcomm=comm
      lsize=size(ldat)
#ifdef r8i8
      ltype=MPI_COMPLEX8
#else
      ltype=MPI_COMPLEX
#endif 
     
      call MPI_Reduce(ldat, gdat, lsize, ltype, lop, lhost, lcomm, lerr )
   
      call end_log(reduce_end)
   
   end subroutine ccmpi_reduce2c
      
   subroutine ccmpi_allreduce2i(ldat,gdat,op,comm)
   
      integer, intent(in) :: comm
      integer(kind=4) ltype,lop,lcomm,lerr,lsize
      integer, dimension(:), intent(in) :: ldat
      integer, dimension(:), intent(out) :: gdat
      character(len=*), intent(in) :: op
      
      call start_log(reduce_begin)
      
      select case(op)
         case("max")
            lop=MPI_MAX
         case("min")
            lop=MPI_MIN
         case("sum")
            lop=MPI_SUM
         case default
            write(6,*) "ERROR: Unknown option for ccmpi_allreduce ",op
            call MPI_Abort(MPI_COMM_WORLD,-1,lerr)
      end select
      
      lcomm=comm
      lsize=size(ldat)
#ifdef r8i8
      ltype=MPI_INTEGER8
#else
      ltype=MPI_INTEGER
#endif 
     
      call MPI_AllReduce(ldat, gdat, lsize, ltype, lop, lcomm, lerr )
 
      call end_log(reduce_end)
   
   end subroutine ccmpi_allreduce2i

   subroutine ccmpi_allreduce2r(ldat,gdat,op,comm)
   
      integer, intent(in) :: comm
      integer(kind=4) ltype,lop,lcomm,lerr,lsize
      real, dimension(:), intent(in) :: ldat
      real, dimension(:), intent(out) :: gdat
      character(len=*), intent(in) :: op
      
      call start_log(reduce_begin)
      
      lcomm=comm
      lsize=size(ldat)
      
      select case(op)
         case("max")
            lop=MPI_MAX
#ifdef r8i8
            ltype=MPI_DOUBLE_PRECISION
#else
            ltype=MPI_REAL
#endif 
         case("min")
            lop=MPI_MIN
#ifdef r8i8
            ltype=MPI_DOUBLE_PRECISION
#else
            ltype=MPI_REAL
#endif 
         case("sum")
            lop=MPI_SUM
#ifdef r8i8
            ltype=MPI_DOUBLE_PRECISION
#else
            ltype=MPI_REAL
#endif 
         case("maxloc")
            lop=MPI_MAXLOC
            lsize=lsize/2
#ifdef r8i8
            ltype=MPI_2DOUBLE_PRECISION
#else
            ltype=MPI_2REAL
#endif 
         case("minloc")
            lop=MPI_MINLOC
            lsize=lsize/2
#ifdef r8i8
            ltype=MPI_2DOUBLE_PRECISION
#else
            ltype=MPI_2REAL
#endif 
         case default
            write(6,*) "ERROR: Unknown option for ccmpi_allreduce ",op
            call MPI_Abort(MPI_COMM_WORLD,-1,lerr)
      end select
     
      call MPI_AllReduce(ldat, gdat, lsize, ltype, lop, lcomm, lerr )
   
      call end_log(reduce_end)
   
   end subroutine ccmpi_allreduce2r
  
   subroutine ccmpi_allreduce3r(ldat,gdat,op,comm)
   
      integer, intent(in) :: comm
      integer(kind=4) ltype,lop,lcomm,lerr,lsize
      real, dimension(:,:), intent(in) :: ldat
      real, dimension(:,:), intent(out) :: gdat
      character(len=*), intent(in) :: op
      
      call start_log(reduce_begin)
      
      lcomm=comm
      lsize=size(ldat)
      
      select case(op)
         case("max")
            lop=MPI_MAX
#ifdef r8i8
            ltype=MPI_DOUBLE_PRECISION
#else
            ltype=MPI_REAL
#endif 
         case("min")
            lop=MPI_MIN
#ifdef r8i8
            ltype=MPI_DOUBLE_PRECISION
#else
            ltype=MPI_REAL
#endif 
         case("sum")
            lop=MPI_SUM
#ifdef r8i8
            ltype=MPI_DOUBLE_PRECISION
#else
            ltype=MPI_REAL
#endif 
         case("maxloc")
            lop=MPI_MAXLOC
            lsize=lsize/2
#ifdef r8i8
            ltype=MPI_2DOUBLE_PRECISION
#else
            ltype=MPI_2REAL
#endif 
         case("minloc")
            lop=MPI_MINLOC
            lsize=lsize/2
#ifdef r8i8
            ltype=MPI_2DOUBLE_PRECISION
#else
            ltype=MPI_2REAL
#endif 
         case default
            write(6,*) "ERROR: Unknown option for ccmpi_allreduce ",op
            call MPI_Abort(MPI_COMM_WORLD,-1,lerr)
      end select
      
      call MPI_AllReduce(ldat, gdat, lsize, ltype, lop, lcomm, lerr )
   
      call end_log(reduce_end)
   
   end subroutine ccmpi_allreduce3r
   
   subroutine ccmpi_allreduce2c(ldat,gdat,op,comm)
   
      use sumdd_m
   
      integer, intent(in) :: comm
      integer(kind=4) ltype,lop,lcomm,lerr,lsize
      complex, dimension(:), intent(in) :: ldat
      complex, dimension(:), intent(out) :: gdat
      character(len=*), intent(in) :: op
      
      call start_log(reduce_begin)
      
      select case(op)
         case("max")
            lop=MPI_MAX
         case("min")
            lop=MPI_MIN
         case("sum")
            lop=MPI_SUM
         case("sumdr")
            lop=MPI_SUMDR
         case default
            write(6,*) "ERROR: Unknown option for ccmpi_allreduce ",op
            call MPI_Abort(MPI_COMM_WORLD,-1,lerr)
      end select
      
      lcomm=comm
      lsize=size(ldat)
#ifdef r8i8
      ltype=MPI_COMPLEX8
#else
      ltype=MPI_COMPLEX
#endif 
     
      call MPI_AllReduce(ldat, gdat, lsize, ltype, lop, lcomm, lerr )
   
      call end_log(reduce_end)
   
   end subroutine ccmpi_allreduce2c
   
   subroutine ccmpi_abort(ierrin)
   
      integer, intent(in) :: ierrin
      integer(kind=4) lerrin,ierr
      
      lerrin=ierrin
      call MPI_Abort(MPI_COMM_WORLD,lerrin,ierr)
   
   end subroutine ccmpi_abort

   subroutine ccmpi_bcast2i(ldat,host,comm)

      integer, intent(in) :: host,comm
      integer(kind=4) ltype,lcomm,lhost,lerr,lsize
      integer, dimension(:), intent(in) :: ldat

      call start_log(bcast_begin)

      lhost=host
      lcomm=comm
      lsize=size(ldat)
#ifdef r8i8
      ltype=MPI_INTEGER8
#else
      ltype=MPI_INTEGER
#endif 
   
      call MPI_Bcast(ldat,lsize,ltype,lhost,lcomm,lerr)
         
      call end_log(bcast_end)
         
   end subroutine ccmpi_bcast2i

   subroutine ccmpi_bcast3i(ldat,host,comm)

      integer, intent(in) :: host,comm
      integer(kind=4) ltype,lcomm,lhost,lerr,lsize
      integer, dimension(:,:), intent(in) :: ldat

      call start_log(bcast_begin)

      lhost=host
      lcomm=comm
      lsize=size(ldat)
#ifdef r8i8
      ltype=MPI_INTEGER8
#else
      ltype=MPI_INTEGER
#endif   
   
      call MPI_Bcast(ldat,lsize,ltype,lhost,lcomm,lerr)
      
      call end_log(bcast_end)   
      
   end subroutine ccmpi_bcast3i
   
   subroutine ccmpi_bcast2r(ldat,host,comm)
   
      integer, intent(in) :: host,comm
      integer(kind=4) ltype,lcomm,lhost,lerr,lsize
      real, dimension(:), intent(in) :: ldat

      call start_log(bcast_begin)

      lhost=host
      lcomm=comm
      lsize=size(ldat)
#ifdef r8i8
      ltype=MPI_DOUBLE_PRECISION
#else
      ltype=MPI_REAL
#endif 
   
      call MPI_Bcast(ldat,lsize,ltype,lhost,lcomm,lerr)
   
      call end_log(bcast_end)
   
   end subroutine ccmpi_bcast2r

   subroutine ccmpi_bcast3r(ldat,host,comm)
   
      integer, intent(in) :: host,comm
      integer(kind=4) ltype,lcomm,lhost,lerr,lsize
      real, dimension(:,:), intent(in) :: ldat

      call start_log(bcast_begin)

      lhost=host
      lcomm=comm
      lsize=size(ldat)
#ifdef r8i8
      ltype=MPI_DOUBLE_PRECISION
#else
      ltype=MPI_REAL
#endif   
   
      call MPI_Bcast(ldat,lsize,ltype,lhost,lcomm,lerr)
   
      call end_log(bcast_end)
   
   end subroutine ccmpi_bcast3r

   subroutine ccmpi_bcast4r(ldat,host,comm)
   
      integer, intent(in) :: host,comm
      integer(kind=4) ltype,lcomm,lhost,lerr,lsize
      real, dimension(:,:,:), intent(in) :: ldat

      call start_log(bcast_begin)

      lhost=host
      lcomm=comm
      lsize=size(ldat)
#ifdef r8i8
      ltype=MPI_DOUBLE_PRECISION
#else
      ltype=MPI_REAL
#endif  
   
      call MPI_Bcast(ldat,lsize,ltype,lhost,lcomm,lerr)
   
      call end_log(bcast_end)
   
   end subroutine ccmpi_bcast4r

   subroutine ccmpi_bcast5r(ldat,host,comm)
   
      integer, intent(in) :: host,comm
      integer(kind=4) ltype,lcomm,lhost,lerr,lsize
      real, dimension(:,:,:,:), intent(in) :: ldat

      call start_log(bcast_begin)

      lhost=host
      lcomm=comm
      lsize=size(ldat)
#ifdef r8i8
      ltype=MPI_DOUBLE_PRECISION
#else
      ltype=MPI_REAL
#endif
      call MPI_Bcast(ldat,lsize,ltype,lhost,lcomm,lerr)
   
      call end_log(bcast_end)
   
   end subroutine ccmpi_bcast5r

   subroutine ccmpi_bcast2s(ldat,host,comm)
   
      integer, intent(in) :: host,comm
      integer(kind=4) ltype,lcomm,lhost,lerr,lsize
      character(len=*), dimension(:), intent(in) :: ldat

      call start_log(bcast_begin)

      lhost=host
      lcomm=comm
      lsize=size(ldat)*len(ldat(1))
      ltype=MPI_CHARACTER
   
      call MPI_Bcast(ldat,lsize,ltype,lhost,lcomm,lerr)
   
      call end_log(bcast_end)
   
   end subroutine ccmpi_bcast2s
   
   subroutine ccmpi_bcast2r8(ldat,host,comm)
   
      integer, intent(in) :: host,comm
      integer(kind=4) ltype,lcomm,lhost,ierr,lsize
      double precision, dimension(:), intent(in) :: ldat
   
      call start_log(bcast_begin)
      
      lhost=host
      lcomm=comm
      lsize=size(ldat)
      ltype=MPI_DOUBLE_PRECISION
   
      call MPI_Bcast(ldat,lsize,ltype,lhost,lcomm,ierr)
   
      call end_log(bcast_end)
   
   end subroutine ccmpi_bcast2r8
   
   subroutine ccmpi_bcast3r8(ldat,host,comm)
   
      integer, intent(in) :: host,comm
      integer(kind=4) ltype,lcomm,lhost,ierr,lsize
      double precision, dimension(:,:), intent(in) :: ldat
   
      call start_log(bcast_begin)
      
      lhost=host
      lcomm=comm
      lsize=size(ldat)
      ltype=MPI_DOUBLE_PRECISION
   
      call MPI_Bcast(ldat,lsize,ltype,lhost,lcomm,ierr)
   
      call end_log(bcast_end)
   
   end subroutine ccmpi_bcast3r8   
   
   subroutine ccmpi_bcast4r8(ldat,host,comm)
   
      integer, intent(in) :: host,comm
      integer(kind=4) ltype,lcomm,lhost,ierr,lsize
      double precision, dimension(:,:,:), intent(in) :: ldat

      call start_log(bcast_begin)

      lhost=host
      lcomm=comm
      lsize=size(ldat)
      ltype=MPI_DOUBLE_PRECISION
   
      call MPI_Bcast(ldat,lsize,ltype,lhost,lcomm,ierr)
   
      call end_log(bcast_end)
   
   end subroutine ccmpi_bcast4r8

   subroutine ccmpi_barrier(comm)
   
      integer, intent(in) :: comm
      integer(kind=4) lcomm,ierr
      
      lcomm=comm
      call MPI_Barrier( lcomm, ierr )
   
   end subroutine ccmpi_barrier
   
   subroutine ccmpi_gatherx2r(gdat,ldat,host,comm)

      integer, intent(in) :: host,comm
      integer(kind=4) lsize,ltype,lhost,lcomm,lerr
      real, dimension(:), intent(out) :: gdat
      real, dimension(:), intent(in) :: ldat

#ifdef debug        
      if (size(gdat)/=size(ldat)*nproc) then
         write(6,*) "ERROR: Incorrect size for ccmpi_gather"
         call MPI_Abort(MPI_COMM_WORLD,-1,lerr)
      end if
#endif
      
      lcomm=comm
      lhost=host
      lsize=size(ldat)
#ifdef r8i8
      ltype=MPI_DOUBLE_PRECISION
#else
      ltype=MPI_REAL
#endif     
   
      call MPI_Gather(ldat,lsize,ltype,gdat,lsize,ltype,lhost,lcomm,lerr)
   
   end subroutine ccmpi_gatherx2r
   
   subroutine ccmpi_scatterx2r(gdat,ldat,host,comm)
   
      integer, intent(in) :: host,comm
      integer(kind=4) lsize,ltype,lhost,lcomm,lerr
      real, dimension(:), intent(in) :: gdat
      real, dimension(:), intent(out) :: ldat

#ifdef debug        
      if (size(gdat)/=size(ldat)*nproc) then
         write(6,*) "ERROR: Incorrect size for ccmpi_scatter"
         call MPI_Abort(MPI_COMM_WORLD,-1,lerr)
      end if
#endif
      
      lcomm=comm
      lhost=host
      lsize=size(ldat)
#ifdef r8i8
      ltype=MPI_DOUBLE_PRECISION
#else
      ltype=MPI_REAL
#endif     
   
      call MPI_Scatter(gdat,lsize,ltype,ldat,lsize,ltype,lhost,lcomm,lerr)
   
   end subroutine ccmpi_scatterx2r
   
   subroutine ccmpi_allgatherx2i(gdat,ldat,comm)
   
      integer, intent(in) :: comm
      integer(kind=4) lsize,ltype,lcomm,lerr
      integer, dimension(:), intent(in) :: ldat
      integer, dimension(:), intent(out) :: gdat

#ifdef debug      
      if (size(gdat)/=size(ldat)*nproc) then
         write(6,*) "ERROR: Incorrect size for ccmpi_allgather"
         call MPI_Abort(MPI_COMM_WORLD,-1,lerr)
      end if
#endif
   
      lcomm=comm
      lsize=size(ldat)
#ifdef r8i8
      ltype=MPI_INTEGER8
#else
      ltype=MPI_INTEGER
#endif  
      
      call MPI_AllGather(ldat,lsize,ltype,gdat,lsize,ltype,lcomm,lerr)
      
   end subroutine ccmpi_allgatherx2i
   
   subroutine ccmpi_recv2r(ldat,iproc,itag,comm)
   
      integer, intent(in) :: iproc,itag,comm
      integer(kind=4) lproc,ltag,lcomm,lerr,lsize,ltype
      integer(kind=4), dimension(MPI_STATUS_SIZE) :: lstatus
      real, dimension(:), intent(out) :: ldat
   
      lproc=iproc
      ltag=itag
      lcomm=comm
      lsize=size(ldat)      
#ifdef r8i8
      ltype=MPI_DOUBLE_PRECISION
#else
      ltype=MPI_REAL
#endif
      
      call MPI_Recv(ldat,lsize,ltype,lproc,ltag,lcomm,lstatus,lerr)
   
   end subroutine ccmpi_recv2r
   
   subroutine ccmpi_ssend2r(ldat,iproc,itag,comm)
   
      integer, intent(in) :: iproc,itag,comm
      integer(kind=4) lproc,ltag,lcomm,lerr,lsize,ltype
      integer(kind=4), dimension(MPI_STATUS_SIZE) :: lstatus
      real, dimension(:), intent(in) :: ldat

      lproc=iproc
      ltag=itag
      lcomm=comm
      lsize=size(ldat)      
#ifdef r8i8
      ltype=MPI_DOUBLE_PRECISION
#else
      ltype=MPI_REAL
#endif
  
      call MPI_SSend(ldat,lsize,ltype,lproc,ltag,lcomm,lerr)
   
   end subroutine ccmpi_ssend2r

   subroutine ccmpi_init

      integer(kind=4) lerr,lproc,lid

      call MPI_Init(lerr)
      call MPI_Comm_size(MPI_COMM_WORLD, lproc, lerr) ! Find number of processes
      call MPI_Comm_rank(MPI_COMM_WORLD, lid, lerr)   ! Find local processor id
      nproc=lproc
      myid=lid
      comm_world=MPI_COMM_WORLD
   
   end subroutine ccmpi_init
   
   subroutine ccmpi_finalize
   
      integer(kind=4) lerr
   
      call MPI_Finalize(lerr)
   
   end subroutine ccmpi_finalize

   subroutine ccmpi_commsplit(commout,comm,colour,rank)
   
      integer, intent(out) :: commout
      integer, intent(in) :: comm,colour,rank
      integer(kind=4) lcomm,lcommout,lerr,lrank,lcolour
   
      lcomm=comm
      lcolour=colour
      lrank=rank
      call MPI_Comm_Split(lcomm,lcolour,lrank,lcommout,lerr)
      commout=lcommout
   
   end subroutine ccmpi_commsplit
   
   subroutine ccmpi_commfree(comm)
   
      integer, intent(in) :: comm
      integer(kind=4) lcomm,lerr
      
      lcomm=comm
      call MPI_Comm_Free(lcomm,lerr)
   
   end subroutine ccmpi_commfree

end module cc_mpi

