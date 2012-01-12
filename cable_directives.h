!=============================================================
!--- this file should be included in any src file which attempts to use 
!--- any of the below defs or else's. include UPDEF here #ifdef CABLE
!=============================================================
#define ONLINE_UM 
!#define ONLINE_Mk3L 
!#define OFFLINE_CABLE 

#define FDIAG_SOIL_RESP "on"
!--- leaf respiration to be included in calc of GPP_tile currently in cable_IMunpack
#define FLEAF_RESPIRATION "off"

!---this is only here for now as jhan_ntiles is scattered thru the code
#define jhan_ntiles

!--- def number soil types (cable_....)
!#  define JHNSOIL 10
! MJT change for CCAM
#  define JHNSOIL 13

#ifdef jhan_ntiles
#  define JHNTILES 17
#  define jhan_soil_layers
#  ifdef jhan_soil_layers
#     define six_soil_layers
#  endif
#endif

#ifndef jhan_ntiles
#  define JHNTILES 9
#endif

!--- cable_common vars
#define cable_common

!=============================================================








!=============================================================
!!--- USE CABLE_DIAG  
!=============================================================
!--- below directives enable generic build/use of diag 
!--- : cable_diag.F90 - includes module in build
!
!#define cable_diag
!#ifdef cable_diag
#define inc_cable_diag
!#endif

