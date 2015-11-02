! Conformal Cubic Atmospheric Model
    
! Copyright 2015 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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


! This is a cut down version of the MPICH mpif.h, suitable for CCAM using
! mpidummy.f90 in place of a real MPI library.

!
!     The types MPI_INTEGER1,2,4 and MPI_REAL4,8 are OPTIONAL.
!     Their values are zero if they are not available.  Note that
!     using these reduces the portability of code.
!
      integer, parameter :: MPI_CHARACTER =  1
      integer, parameter :: MPI_LOGICAL   = 25
      integer, parameter :: MPI_REAL      = 26
      integer, parameter :: MPI_2REAL     = 32
      integer, parameter :: MPI_DOUBLE_PRECISION = 27
      integer, parameter :: MPI_2DOUBLE_PRECISION = 1275072547
      integer, parameter :: MPI_INTEGER   = 28
      integer, parameter :: MPI_INTEGER8  = 1275070513
      integer, parameter :: MPI_COMPLEX   = 1275070494
      integer, parameter :: MPI_COMPLEX8  = 1275070504
      integer, parameter :: MPI_COMPLEX16 = 1275072554
      integer, parameter :: MPI_DOUBLE_COMPLEX = z'4c001022'
      integer, parameter :: MPI_PACKED    = z'4c00010f'
      integer, parameter :: MPI_BYTE      = z'4c00010d'
!
!     Collective operations
!
      integer, parameter :: MPI_MAX    = 100
      integer, parameter :: MPI_MIN    = 101
      integer, parameter :: MPI_SUM    = 102
      integer, parameter :: MPI_MINLOC = 110
      integer, parameter :: MPI_MAXLOC = 111
      integer, parameter :: MPI_LAND   = 1476395013
      integer, parameter :: MPI_LOR    = 1476395015
!
!     Communicator
!
      integer, parameter :: MPI_COMM_WORLD  = 91
      integer, parameter :: MPI_COMM_NULL   = 67108864
      integer, parameter :: MPI_COMM_TYPE_SHARED = 1
      integer, parameter :: MPI_UNDEFINED   = -32766
      
      integer, parameter :: MPI_STATUS_SIZE = 6
      integer, dimension(MPI_STATUS_SIZE), parameter :: MPI_STATUS_IGNORE = 0
      integer, dimension(MPI_STATUS_SIZE,1), parameter :: MPI_STATUSES_IGNORE = 0
      integer, parameter :: MPI_ADDRESS_KIND = INT_PTR_KIND()
      integer, parameter :: MPI_INFO_NULL = 469762048
      integer, parameter :: MPI_MODE_NOPUT = 4096
      integer, parameter :: MPI_MODE_NOPRECEDE = 8192
      integer, parameter :: MPI_MODE_NOSUCCEED = 16384

!
!     All other MPI routines are subroutines
!
      double precision mpi_wtime, mpi_wtick
      external         MPI_WTIME, MPI_WTICK
