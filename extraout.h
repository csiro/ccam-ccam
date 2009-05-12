      real, dimension(ifull) :: cloudlo, cloudmi, cloudhi, cloudtot,      &
     &                          rgsave, rtsave, sintsave, sgsave,         &
     &                          rtclsave, sgclsave, taux, tauy, ustar,    &
     &                          swrsave                                     ! MJT cable
      real, dimension(ifull,8) :: u10_3hr,v10_3hr,tscr_3hr,rh1_3hr
      common/extraout/cloudlo,cloudmi,cloudhi,cloudtot,rgsave,rtsave,     &
     &     sintsave,sgsave,rtclsave,sgclsave,taux,tauy,ustar,             &
     &     u10_3hr,v10_3hr,tscr_3hr,rh1_3hr,swrsave                         ! MJT cable

