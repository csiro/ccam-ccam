! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2024 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
module cc_omp

#ifdef _OPENMP
   use omp_lib, only : omp_get_max_threads, omp_get_thread_num,     &
                       omp_get_default_device, omp_get_num_devices, &
                       omp_set_default_device
#endif

   implicit none
   private

   integer, save, public :: maxthreads, ntiles, imax
#ifdef GPU
   integer, save, public :: maxtilesize = 32 ! suggested value
#else
   integer, save, public :: maxtilesize = 96 ! suggested value
#endif   

   public ::  ccomp_init
   public ::  ccomp_get_thread_num

   contains
    
   function ccomp_get_thread_num() result(threadn)
      integer :: threadn

#ifdef _OPENMP
      threadn = omp_get_thread_num()
#else
      threadn = 0
#endif

   end function ccomp_get_thread_num

   subroutine ccomp_init(myid)
      use newmpar_m, only : ifull
      integer :: i, tmp
      integer, intent(in) :: myid
   
#ifdef _OPENMP
      maxthreads = omp_get_max_threads()
#else
      maxthreads = 1
#endif

      !find imax if maxtilesize isn't already a factor of ifull
      imax = min( max( maxtilesize, 1 ), ifull )
      tmp = imax
      imax = -1 ! missing flag
      ! first attempt to find multiple of 16
      do i = tmp,16,-1
         if ( mod(ifull,i)==0 .and. mod(i,16)==0 ) then
            imax = i
            exit
         end if
      end do
      if ( imax<1 ) then
         ! second attempt if multiple of 16 is not possible
         do i = tmp,1,-1
            if ( mod(ifull,i)==0 ) then
               imax = i
               exit
            end if
         end do
      end if

      !find the number of tiles
      ntiles = ifull/imax
      
      return
   end subroutine ccomp_init
 
end module cc_omp
