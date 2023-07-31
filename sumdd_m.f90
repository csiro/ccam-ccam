! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2023 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
module sumdd_m
   implicit none
   public drpdr, drpdr_local
   private
   integer(kind=4), save, public :: MPI_SUMDR
   interface drpdr_local
     module procedure drpdr_local2r, drpdr_local3r
     module procedure drpdr_local2c
   end interface drpdr_local
contains
   subroutine drpdr(dra, drb, len, itype) 
!  Modification of original codes written by David H. Bailey. 
!  This subroutine computes drb(i) = dra(i) + drb(i) 
!  From He and Ding 2001
!  This version is modified to just use double real precision rather
!  than double double. This makes it work with the default real precision
!  on a range of machines.
!  Here we're more concerned with reproducibility rather than accuracy
!  so there's no need for double precision.

      integer(kind=4), intent(in) :: len
      real :: e, t1, t2 
      integer :: i
      complex, dimension(len), intent(in)  :: dra
      complex, dimension(len), intent(inout) :: drb
      integer(kind=4), intent(in) :: itype

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

   pure subroutine drpdr_local2r(array, local_sum)
   ! This is a local version of drpdr that takes an array of reals on 
   ! one processor and returns the double-real sum
   ! Note that it accumulates into local_sum so this has to be zeroed
   ! before use.
      implicit none 
      real, dimension(:), intent(in)  :: array
      complex, intent(inout) :: local_sum
      real :: e, t1, t2 
      integer :: i
      
      do i = 1,size(array)
         t1 = array(i) + real(local_sum) 
         e = t1 - array(i) 
         t2 = ((real(local_sum) - e) + (array(i) - (t1 - e)))  + aimag(local_sum)
         local_sum = cmplx (t1 + t2, t2 - ((t1 + t2) - t1)) 
      end do
      
   end subroutine drpdr_local2r
   
   pure subroutine drpdr_local3r(array, local_sum)
   ! This is a local version of drpdr that takes an array of reals on 
   ! one processor and returns the double-real sum
   ! Note that it accumulates into local_sum so this has to be zeroed
   ! before use.
      implicit none 
      real, dimension(:,:), intent(in)  :: array
      complex, dimension(:), intent(inout) :: local_sum
      real :: e, t1, t2 
      real, dimension(size(array,2),size(array,1)) :: array_t
      integer :: i, n
      
      array_t(:,:) = transpose(array)

      do i = 1,size(array_t,2)
         do n = 1,size(array_t,1)
            t1 = array_t(n,i) + real(local_sum(n))
            e = t1 - array_t(n,i) 
            t2 = ((real(local_sum(n)) - e) + (array_t(n,i) - (t1 - e)))  + aimag(local_sum(n))
            local_sum(n) = cmplx(t1 + t2, t2 - ((t1 + t2) - t1))
         end do
      end do
      
   end subroutine drpdr_local3r

   pure subroutine drpdr_local2c(array, local_sum)
   ! This is a local version of drpdr that takes an array of reals on 
   ! one processor and returns the double-real sum
   ! Note that it accumulates into local_sum so this has to be zeroed
   ! before use.
      implicit none 
      complex, dimension(:), intent(in)  :: array
      complex, intent(inout) :: local_sum
      real :: e, t1, t2 
      integer :: i
      
      do i = 1,size(array)
         t1 = real(array(i)) + real(local_sum) 
         e = t1 - real(array(i)) 
         t2 = ((real(local_sum) - e) + (real(array(i)) - (t1 - e))) &
              + aimag(array(i)) + aimag(local_sum)
         local_sum = cmplx (t1 + t2, t2 - ((t1 + t2) - t1)) 
      end do
      
   end subroutine drpdr_local2c
   
end module sumdd_m

