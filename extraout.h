      real, dimension(ifull) :: cloudlo, cloudmi, cloudhi, cloudtot,      &
     &                          rgsave, rtsave, sintsave, sgsave,         &
     &                          rtclsave, sgclsave, taux, tauy, ustar,    &
     &                          swrsave,fbeamvis,fbeamnir                   ! MJT cable ! MJT radiation
      real, dimension(ifull,8) :: u10_3hr,v10_3hr,tscr_3hr,rh1_3hr
      common/extraout/cloudlo,cloudmi,cloudhi,cloudtot,rgsave,rtsave,     &
     &     sintsave,sgsave,rtclsave,sgclsave,taux,tauy,ustar,             &
     &     u10_3hr,v10_3hr,tscr_3hr,rh1_3hr,swrsave,fbeamvis,fbeamnir       ! MJT cable ! MJT radiation

