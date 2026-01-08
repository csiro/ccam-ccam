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
   public ccacc_init

   integer, save, public :: async_length = 2
   
   contains

   subroutine ccacc_init(myid,ngpus)

     integer, intent(in) :: myid
     integer, intent(inout) :: ngpus
     integer(kind=4) :: lmyid, lngpus
#ifdef _OPENACC
     integer(kind=4) :: device_num
     integer(acc_device_kind) :: devicetype

     lmyid = myid
     devicetype = acc_get_device_type()
     lngpus = acc_get_num_devices(devicetype)
     ngpus = lngpus

     if ( ngpus > 0 ) then
        call acc_set_device_num(mod(lmyid,lngpus),devicetype)
        !lgpuid = acc_get_device_num(devicetype)
     end if   
#endif

   end subroutine ccacc_init
 
end module cc_acc
