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

   integer, save, public :: maxthreads

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
      integer, intent(in) :: myid
   
#ifdef _OPENMP
      maxthreads = omp_get_max_threads()
#else
      maxthreads = 1
#endif
      
      return
   end subroutine ccomp_init
 
end module cc_omp
