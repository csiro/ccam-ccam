!  eak version 16/03/06
      real, dimension(ifull) :: cloudlo, cloudmi, cloudhi, cloudtot,
     &                          rgsave, rtsave, sintsave, sgsave,
     &                          rtclsave, sgclsave, taux, tauy, ustar,
     &                          coszros
      real, dimension(ifull,8) :: u10_3hr,v10_3hr
      real, dimension(ifull,4) :: tscr_6hr,rh1_6hr
      common/extraout/cloudlo,cloudmi,cloudhi,cloudtot,rgsave,rtsave,
     &     sintsave,sgsave,rtclsave,sgclsave,taux,tauy,ustar,
     &     u10_3hr,v10_3hr,tscr_6hr,rh1_6hr,coszros

