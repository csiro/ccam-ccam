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
    
module usage_m
   implicit none
   contains
   subroutine usage()
   use cc_mpi                                 ! CC MPI routines
   if (myid==0) then
      write(*,"(a)") &
"Usage: mpirun -np nproc globpea [-h] [-c input_file]", &
"  globpea -h for full list of options and more information."
   end if
   call ccmpi_barrier(comm_world)
   call ccmpi_abort(-1)
   end subroutine usage

   subroutine help
   use cc_mpi                                 ! CC MPI routines
   if (myid==0) then
      write(*,"(a)") &
"Conformal Cubic Atmospheric Model (CCAM)", &
"", &
"Usage: mpirun -np nproc globpea [-h] [-c input_file]",&
"", &
"Command line options are", &
"", &
" -h for help (this message)", &
"", &
" -c input_file where input_file is the input namelist file", &
"   (Default is input).", &
""
   end if
   
   call ccmpi_barrier(comm_world)
   call ccmpi_finalize
   call finishbanner
   stop

   end subroutine help
end module usage_m
