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
   use omp_lib, only : omp_get_max_threads
#endif

   implicit none
   private

#ifdef _OPENMP
   logical, parameter, public :: using_omp = .true.
#else
   logical, parameter, public :: using_omp = .false.
#endif
   integer, save, public :: maxthreads, ntiles, imax
   integer, save, public :: maxtilesize = huge(1) ! actually 96 is preffered

   public ::  ccomp_init
   public ::  ccomp_ntiles

   contains

   function ccomp_get_max_threads() result(nthreads)
      integer :: nthreads

#ifdef _OPENMP
      nthreads = omp_get_max_threads()
#else
      nthreads = 1
#endif

   end function ccomp_get_max_threads

   subroutine ccomp_init

      maxthreads = ccomp_get_max_threads()

   end subroutine ccomp_init
      
   subroutine ccomp_ntiles
      use newmpar_m, only : ifull
      integer :: i, tmp

      !find a tiling at least as much as the number of threads 
      ntiles = ifull
      do i = maxthreads,ifull
         if ( mod(ifull,i) == 0 ) then
            ntiles = i
            exit
         end if
      end do

      !find the next biggest maxtilesize if maxtilesize isn't already a factor of ifull
      maxtilesize = min( max( maxtilesize, 1 ), ifull )
      tmp = maxtilesize
      do i = tmp,ifull
         if ( mod(ifull,i) == 0 ) then
            maxtilesize = i
            exit
         end if
      end do

      !increase the number of tiles if the resultant tile size is too big
      if ( ifull/ntiles > maxtilesize ) ntiles = ifull/maxtilesize

      imax = ifull/ntiles
      
      return
   end subroutine ccomp_ntiles
 
end module cc_omp
