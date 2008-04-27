! eak 16/03/06
! carbon pools
!      real, dimension(ifull) ::carb_lf,carb_wd,carb_rts
!      real, dimension(ifull) ::carb_slf,carb_sls
      integer :: inyear_carb  ! year to initialise carbon pools
      real, dimension(ncp) ::ratecp
      real, dimension(ifull,ncp) ::cplant
      real, dimension(ncs) ::ratecs
      real, dimension(ifull,ncs) ::csoil
      real, dimension(ifull,12) ::csoil01
      real, dimension(ifull) ::csoil02
      real, dimension(mxvt,ncp) ::tcplant
      real, dimension(mxvt,ncs) ::tcsoil
      real, dimension(ifull) ::fnee,fpn,frd,frp,frpw,frpr,frs

!      common/carbpools/ carb_lf,carb_wd,carb_rts,carb_slf,carb_sls
      common/carbpools/ cplant,ratecp
      common/carbpools/ csoil,ratecs
      common/carbpools/ csoil01,csoil02
      common/carbpools/ inyear_carb
      common/carbflux/ fnee,fpn,frd,frp,frpw,frpr,frs

 

