      real, dimension(ifull) :: eg_ave, fg_ave, ga_ave,epan_ave,dew_ave,  &
     &                 cbas_ave, ctop_ave, rndmax, qscrn_ave,             &
     &                 tmaxscr, tminscr, tscr_ave,                        &
     &                 rhmaxscr, rhminscr,                                &
     &                 riwp_ave,rlwp_ave,u10max,v10max,u10mx,             &
     &                 u1max, v1max, u2max, v2max, capemax, epot_ave
      real, dimension(ifull,ms) :: wb_ave,tgg_ave
      real, dimension(ifull) :: theta_ave,fpn_ave,frday_ave,frp_ave,      &
     &                          rnet_ave,tsu_ave,alb_ave
      common /histave/ eg_ave, fg_ave, ga_ave, epan_ave, dew_ave,         &
     &                 cbas_ave, ctop_ave, rndmax, qscrn_ave,             &
     &                 tmaxscr, tminscr, tscr_ave, riwp_ave, rlwp_ave,    &
     &                 u10max, v10max, rhmaxscr, rhminscr, u10mx,         &
     &                 u1max, v1max, u2max, v2max, capemax, epot_ave,     &
     &                 theta_ave,wb_ave,tgg_ave,fpn_ave,frday_ave,        &
     &                 frp_ave,rnet_ave,tsu_ave,alb_ave
