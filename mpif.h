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
      integer, parameter :: MPI_DOUBLE_PRECISION = 27
      integer, parameter :: MPI_INTEGER   = 28
      integer, parameter :: MPI_2REAL     = 32
!
!     Collective operations
!
      integer, parameter :: MPI_MAX    = 100
      integer, parameter :: MPI_MIN    = 101
      integer, parameter :: MPI_SUM    = 102
      integer, parameter :: MPI_MINLOC = 110
      integer, parameter :: MPI_MAXLOC = 111
!
!     Communicator
!
      integer, parameter :: MPI_COMM_WORLD  = 91
      
      integer, parameter :: MPI_STATUS_SIZE = 6

!
!     All other MPI routines are subroutines
!
      double precision mpi_wtime, mpi_wtick
      external         MPI_WTIME, MPI_WTICK
