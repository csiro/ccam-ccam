module sumdd_m
   implicit none
   integer(kind=4), save :: MPI_SUMDR, MPI_SUMDRA, MPI_MAXMIN
contains
   subroutine drpdr (dra, drb, len, itype) 
!  Modification of original codes written by David H. Bailey. 
!  This subroutine computes drb(i) = dra(i) + drb(i) 
!  From He and Ding 2001
!  This version is modified to just use double real precision rather
!  than double double. This makes it work with the default real precision
!  on a range of machines.
!  Here we're more concerned with reproducibility rather than accuracy
!  so there's no need for double precision.

!  With ifort, this has to be compiled with -mp 
! (no other optimisation options, not even -g)
!  On APAC it works with default optimisation and on the NEC with -Cvopt.

      integer(kind=4), intent(in) :: len, itype
      complex, dimension(len), intent(in)  :: dra
      complex, dimension(len), intent(inout) :: drb
      real :: e, t1, t2 
      integer :: i

      do i = 1, len 
         !  Compute dra + drb using Knuth's trick. 
         t1 = real(dra(i)) + real(drb(i)) 
         e = t1 - real(dra(i)) 
         t2 = ((real(drb(i)) - e) + (real(dra(i)) - (t1 - e)))  + &
              aimag(dra(i)) + aimag(drb(i)) 
         !    The result is t1 + t2, after normalization. 
         drb(i) = cmplx (t1 + t2, t2 - ((t1 + t2) - t1)) 
      end do
   end subroutine drpdr

   subroutine drpdra (dra, drb, len, itype) 
!  Array version of drpdr

      integer(kind=4), intent(in) :: len, itype
      complex, dimension(2,len/2), intent(in)  :: dra
      complex, dimension(2,len/2), intent(inout) :: drb
      real, dimension(2) :: e, t1, t2 
      integer :: i

      do i = 1, len/2
         !  Compute dra + drb using Knuth's trick. 
         t1(:) = real(dra(:,i)) + real(drb(:,i)) 
         e(:) = t1(:) - real(dra(:,i)) 
         t2(:) = ((real(drb(:,i)) - e(:)) + (real(dra(:,i)) - (t1(:) - e(:))))  + &
                 aimag(dra(:,i)) + aimag(drb(:,i)) 
         !    The result is t1 + t2, after normalization. 
         drb(:,i) = cmplx (t1(:) + t2(:), t2(:) - ((t1(:) + t2(:)) - t1(:))) 
      end do
   end subroutine drpdra

   subroutine drpdr_local (array, local_sum)
   ! This is a local version of drpdr that takes an array of reals on 
   ! one processor and returns the double-real sum
   ! Note that it accumulates into local_sum so this has to be zeroed
   ! before use.
      implicit none 
      real, dimension(:), intent(in)  :: array
      complex, intent(inout) :: local_sum
      real :: e, t1, t2 
      integer :: i
      
      do i=1,size(array)
         t1 = array(i) + real(local_sum) 
         e = t1 - array(i) 
         t2 = ((real(local_sum) - e) + (array(i) - (t1 - e)))  + &
              aimag(local_sum) 
         local_sum = cmplx (t1 + t2, t2 - ((t1 + t2) - t1)) 
      end do
      
   end subroutine drpdr_local

   subroutine maxmin (dra, drb, len, itype) 
!  Combine max and min into a single operation

      integer(kind=4), intent(in) :: len, itype
      real, dimension(2,len/2), intent(in)  :: dra
      real, dimension(2,len/2), intent(inout) :: drb

      drb(1,:)=max(dra(1,:),drb(1,:))
      drb(2,:)=min(dra(2,:),drb(2,:))

   end subroutine maxmin

end module sumdd_m

