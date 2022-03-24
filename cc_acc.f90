! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2022 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
module cc_acc

#ifdef _OPENACC
   use openacc, only : acc_get_device_type, acc_get_device_type, &
                       acc_get_num_devices, acc_set_device_num, &
                       acc_get_device_num, acc_device_kind
#endif

   implicit none
   private

   integer, save, private :: gpuid = -1
   integer, save, public :: async_length = 2

   public ::  ccacc_init

   contains

   subroutine ccacc_init(myid,ngpus)

     integer, intent(in) :: myid
     integer, intent(inout) :: ngpus

#ifndef _OPENMP     
#ifdef _OPENACC
     integer :: device_num
     integer(acc_device_kind) :: devicetype

     devicetype = acc_get_device_type()
     ngpus = acc_get_num_devices(devicetype)

     if ( ngpus > 0 ) then
        call acc_set_device_num(mod(myid,ngpus),devicetype)
        gpuid = acc_get_device_num(devicetype)
     end if   
#endif
#endif

   end subroutine ccacc_init
 
end module cc_acc
