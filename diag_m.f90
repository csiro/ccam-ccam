module diag_m
   implicit none
   public :: printa, maxmin, average
   private :: maxmin1, maxmin2
   interface maxmin
      module procedure maxmin1, maxmin2
   end interface
contains
   subroutine printa(name,a,ktau,level,i1,i2,j1,j2,bias,facti)
      ! printa has to work with arrays dimension both ifull and ifull+iextra
      ! printb has printj entry, and automatic choice of fact if facti=0.
      ! printa is always passed a 1D section. Level isn't used. Really need
      ! overloaded version. But for now, just use 1D
      use cc_mpi
      include 'newmpar.h'
      character(len=*), intent(in) :: name
      real, dimension(:), intent(in) :: a
      integer, intent(in) :: ktau, level, i1, i2, j1, j2
      real, intent(in) :: bias, facti
      integer jrow, k1, k2
      integer i, ia, ib, j, ja, jb, k, n, nlocal
      real fact

      ! The indices i1, i2, j1, j2 are global
      n = (j1-1)/il ! Global n
      j = j1 - n*il ! Value on face
      if ( myid == fproc(i1,j,n) ) then
         nlocal = n + noff
         ja = j1 - n*il - joff
         jb = j2 - n*il - joff
         jb = min(jb,ja+24,il)
         fact=facti
         ! Calculate factor from the middle of the face 
         ! This will be consistent for 1-6 processors at least
         if(facti.eq.0.) fact=10./abs(a(indp(ipan/2,jpan/2,nlocal)))
         print 9 ,name,ktau,level,bias,fact
 9       format(/1x,a4,' ktau =',i7,'  level =',i3,'  addon =',g8.2,   &      
     & '  has been mult by',1pe8.1)
         print 91,(j,j=j1,j2)
91       format(4x,25i6)
         do i=i1-ioff,i2-ioff
            write(unit=*,fmt="(i5)", advance="no") i+ioff
            do j=ja,jb
               write(unit=*,fmt="(f6.2)", advance="no")                &
     &              (a(indp(i,j,nlocal))-bias)*fact
            end do
            write(*,*)
         end do
      end if

!!!c     following entry prints j cross section in standard orientation
!!!      entry printj(name,a,ktau,jrow,i1,i2,j1,j2,k1,k2,bias,facti)
!!!      fact=facti
!!!      if(facti.eq.0.)fact=10./abs(a((i1+i2)/2,(j1+j2)/2,(k1+k2)/2))
!!!      print 93 ,name,ktau,jrow,bias,fact
!!!93    format(/1x,a4,' ktau =',i7,'  for j =',i3,'  addon =',g8.2,
!!!     1 '  has been mult by',1pe8.1)
!!!      ia=i1
!!!32    ib=min(ia+24,i2)
!!!      print 91,(i,i=ia,ib)
!!!      do 38 k=k2,k1,-1
!!!      do 35 i=ia,ib
!!!35    col(i)=(a(i,jrow,k)-bias)*fact
!!!38    print 92,k,(col(i),i=ia,ib)
!!!      if(ib.eq.i2)return
!!!      ia=ib+1
!!!      go to 32

   end subroutine printa

! has more general format & scaling factor  with subroutine average at bottom
   subroutine maxmin2(u,char,ktau,fact,kup)
      use cc_mpi
      include 'newmpar.h'
      include 'mpif.h'
      character(len=2), intent(in) :: char
      integer, intent(in) :: ktau, kup
      real, intent(in) :: fact
      real, dimension(:,:), intent(in) :: u
      real :: umin(kl),umax(kl)
!      integer iumax(kl),jumax(kl),iumin(kl),jumin(kl)
      integer i, j, k, iq, ierr
      real, dimension(kl) :: gumax, gumin


      do k=1,kup
         umax(k) = maxval(u(1:ifull,k))*fact
         umin(k) = minval(u(1:ifull,k))*fact
      end do
#ifdef scyld
      call MPI_Reduce_( umax, gumax, kup, MPI_REAL, MPI_MAX, 0,         &
     &                   MPI_COMM_WORLD, ierr )
      call MPI_Reduce_( umin, gumin, kup, MPI_REAL, MPI_MIN, 0,         &
     &                  MPI_COMM_WORLD, ierr )
#else
      call MPI_Reduce ( umax, gumax, kup, MPI_REAL, MPI_MAX, 0,         &
     &                  MPI_COMM_WORLD, ierr )
      call MPI_Reduce ( umin, gumin, kup, MPI_REAL, MPI_MIN, 0,         &
     &                  MPI_COMM_WORLD, ierr )
#endif


      if ( myid == 0 ) then
      if(kup.eq.1)then
        print 970,ktau,char,gumax(1),char,gumin(1)
970     format(i7,1x,a2,'max ',f8.3,3x,a2,'min ',f8.3)
!!!        print 9705,ktau,iumax(1),jumax(1),iumin(1),jumin(1)
!!!9705    format(i7,'  posij',i4,i4,10x,i3,i4)
        return
      endif   !  (kup.eq.1)

      if(gumax(1).ge.1000.)then   ! for radon
        print 961,ktau,char,(gumax(k),k=1,kup)
961     format(i7,1x,a2,'max ',10f7.1/(14x,10f7.1)/(14x,10f7.1))
!!!        print 977,ktau,(iumax(k),jumax(k),k=1,kup)
        print 962,ktau,char,(gumin(k),k=1,kup)
962     format(i7,1x,a2,'min ',10f7.1/(14x,10f7.1)/(14x,10f7.1))
!!!        print 977,ktau,(iumin(k),jumin(k),k=1,kup)
      elseif(gumax(kup).gt.30.)then  ! format for T, & usually u,v
        print 971,ktau,char,(gumax(k),k=1,kup)
971     format(i7,1x,a2,'max ',10f7.2/(14x,10f7.2)/(14x,10f7.2))
!!!        print 977,ktau,(iumax(k),jumax(k),k=1,kup)
        print 972,ktau,char,(gumin(k),k=1,kup)
972     format(i7,1x,a2,'min ',10f7.2/(14x,10f7.2)/(14x,10f7.2))
!!!        print 977,ktau,(iumin(k),jumin(k),k=1,kup)
977     format(i7,'  posij',10(i3,i4)/(14x,10(i3,i4))/(14x,10(i3,i4)))
      else  ! for qg & sd
        print 981,ktau,char,(gumax(k),k=1,kup)
981     format(i7,1x,a2,'max ',10f7.3/(14x,10f7.3)/(14x,10f7.3))
!!!        print 977,ktau,(iumax(k),jumax(k),k=1,kup)
        print 982,ktau,char,(gumin(k),k=1,kup)
982     format(i7,1x,a2,'min ',10f7.3/(14x,10f7.3)/(14x,10f7.3))
!!!        print 977,ktau,(iumin(k),jumin(k),k=1,kup)
      endif
      end if ! myid == 0
      return
   end subroutine maxmin2

   subroutine maxmin1(u,char,ktau,fact,kup)
      use cc_mpi
      include 'newmpar.h'
      include 'mpif.h'
      character(len=2), intent(in) :: char
      integer, intent(in) :: ktau, kup
      real, intent(in) :: fact
      real, dimension(:), intent(in) :: u
      real :: umin, umax
!      integer iumax(kl),jumax(kl),iumin(kl),jumin(kl)
      integer i, j, k, iq, ierr
      real :: gumax, gumin


      umax = maxval(u(1:ifull))*fact
      umin = minval(u(1:ifull))*fact

#ifdef scyld
      call MPI_Reduce_( umax, gumax, 1, MPI_REAL, MPI_MAX, 0,         &
     &                  MPI_COMM_WORLD, ierr )
      call MPI_Reduce_( umin, gumin, 1, MPI_REAL, MPI_MIN, 0,         &
     &                  MPI_COMM_WORLD, ierr )
#else
      call MPI_Reduce ( umax, gumax, 1, MPI_REAL, MPI_MAX, 0,         &
     &                  MPI_COMM_WORLD, ierr )
      call MPI_Reduce ( umin, gumin, 1, MPI_REAL, MPI_MIN, 0,         &
     &                  MPI_COMM_WORLD, ierr )
#endif


      if ( myid == 0 ) then
        print 970,ktau,char,gumax,char,gumin
970     format(i7,1x,a2,'max ',f8.3,3x,a2,'min ',f8.3)
!!!        print 9705,ktau,iumax(1),jumax(1),iumin(1),jumin(1)
!!!9705    format(i7,'  posij',i4,i4,10x,i3,i4)
      end if ! myid == 0
      return
   end subroutine maxmin1

   subroutine average(speed,spmean_g,spavge_g)
      include 'mpif.h'
      include 'newmpar.h'
      include 'sigs.h'
      include 'xyzinfo.h'  ! wts
      real, dimension(:,: ), intent(in) :: speed
      real, dimension(:), intent(out) :: spmean_g
      real, intent(out) :: spavge_g
      real, dimension(kl) :: spmean
      integer k, iq, ierr

      do k=1,kl
         spmean(k)=0.
         do iq=1,ifull
            spmean(k) = spmean(k)+speed(iq,k)*wts(iq)
         enddo                  !  iq loop
      end do
#ifdef scyld
      call MPI_Reduce_( spmean, spmean_g, kl, MPI_REAL, MPI_SUM, 0,   &
     &                  MPI_COMM_WORLD, ierr )
#else
      call MPI_Reduce ( spmean, spmean_g, kl, MPI_REAL, MPI_SUM, 0,   &
     &                  MPI_COMM_WORLD, ierr )
#endif

      spavge_g = 0.0
      do k=1,kl
         spavge_g = spavge_g-dsig(k)*spmean_g(k) ! dsig is -ve
      end do

   end subroutine average
end module diag_m
