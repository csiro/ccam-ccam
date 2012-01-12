#include "cable_directives.h"

!========================================================================!
!=== PURPOSE: to enable accessibility of fundamental vars in global   ===!
!=== model thru out program (mainly CABLE).                           ===!
!=== USE: use module in subroutines (SRs) at top level of global model===!
!--- 2 define cable_timestep_ data with vars peculiar to global model ===!
!=== optionally call alias_ SR to give name smore familiar to cable   ===!
!=== people. then again use this module in  any subsequent SR in which===!
!=== you want to access this data.                                    ===!  
!========================================================================!

!jhan:can be used for MUCH more than this

#ifdef cable_common

module cable_common_module
   implicit none 
   integer, save :: cable_gltimestep_i, cable_gltimestep_tot, &
         cable_gltimestep_width, cable_glnode_i
   integer, save :: ktau_gl, kend_gl, knode_gl

   contains

   
   subroutine alias_cable_timestep()
      implicit none 
      ktau_gl = cable_gltimestep_i
      kend_gl = cable_gltimestep_tot
      knode_gl = cable_glnode_i
      return
   end subroutine alias_cable_timestep

   subroutine cable_switch_status(diag_message, diag_var )
      implicit none
      character(len=*), intent(in) :: diag_message, diag_var 
      character(len=68) :: message
         message=trim(trim('cable switch ::')//trim(diag_message))            &
                  //' '//trim(diag_var)
         !kdcorbin, 09/10 - changed knode_gl test from 1
         if(ktau_gl==1 .and. knode_gl==0 ) then
            print *, message
         endif   
      return
   end subroutine cable_switch_status


end module cable_common_module

#endif 
