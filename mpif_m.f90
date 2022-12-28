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
    
module mpif_m

#ifdef usempimod
use mpi
#endif

implicit none
public

#ifndef usempimod
include 'mpif.h'
#endif

#ifdef usempiinterface
interface

subroutine MPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, ierror)
type(*), dimension(..), intent(in) :: sendbuf
type(*), dimension(..), intent(inout) :: recvbuf
integer(kind=4), intent(in) :: sendcount, sendtype, recvcount, recvtype, root, comm
integer(kind=4), intent(out) :: ierror
end subroutine

subroutine MPI_Scatter(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, ierror)
type(*), dimension(..), intent(in) :: sendbuf
type(*), dimension(..), intent(inout) :: recvbuf
integer(kind=4), intent(in) :: sendcount, sendtype, recvcount, recvtype, root, comm
integer(kind=4), intent(out) :: ierror
end subroutine

subroutine MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm, ierror)
type(*), dimension(..), intent(in) :: sendbuf
type(*), dimension(..), intent(inout) :: recvbuf
integer(kind=4), intent(in) :: count, datatype, op, comm
integer(kind=4), intent(out) :: ierror
end subroutine

subroutine MPI_ISend(buf, count, datatype, dest, tag, comm, request, ierror)
type(*), dimension(..), intent(in) :: buf
integer(kind=4), intent(in) :: count, datatype, dest, tag, comm
integer(kind=4), intent(out) :: request, ierror
end subroutine

subroutine MPI_IRecv(buf, count, datatype, source, tag, comm, request, ierror)
type(*), dimension(..), intent(inout) :: buf
integer(kind=4), intent(in) :: count, datatype, source, tag, comm
integer(kind=4), intent(out) :: request, ierror
end subroutine

subroutine MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm, ierror)
type(*), dimension(..), intent(in) :: sendbuf
type(*), dimension(..), intent(inout) :: recvbuf
integer(kind=4), intent(in) :: count, datatype, op, root, comm
integer(kind=4), intent(out) :: ierror
end subroutine

subroutine MPI_Bcast(buf, count, datatype, root, comm, ierror)
type(*), dimension(..), intent(inout) :: buf
integer(kind=4), intent(in) :: count, datatype, root, comm
integer(kind=4), intent(out) :: ierror
end subroutine

subroutine MPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, ierror)
type(*), dimension(..), intent(in) :: sendbuf
type(*), dimension(..), intent(inout) :: recvbuf
integer(kind=4), intent(in) :: sendcount, sendtype, recvcount, recvtype, comm
integer(kind=4), intent(out) :: ierror
end subroutine

end interface
#endif

end module mpif_m
