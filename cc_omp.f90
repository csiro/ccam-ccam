! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2017 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
   use omp_lib, only : omp_get_num_threads, omp_get_max_threads
#endif

   implicit none
   private

   public ::  ccomp_get_num_threads
   public ::  ccomp_get_max_threads

   contains

   function ccomp_get_num_threads result(nthreads)
      integer :: nthreads

#ifdef _OPENMP
      nthreads=omp_get_num_threads()
#else
      nthreads=1
#endif

   end function ccomp_get_num_threads

   function ccomp_get_max_threads result(nthreads)
      integer :: nthreads

#ifdef _OPENMP
      nthreads=omp_get_max_threads()
#else
      nthreads=1
#endif

   end function ccomp_get_max_threads
end module cc_omp
