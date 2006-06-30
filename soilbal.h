! eak 16/03/06
!                                       variables for SEB and Hyd. balance
      real tevap, tprecip, trnoff, totenbal, osnowd0, wbtot0
      common /soilbal/ tevap(ifull)
      common /soilbal/ tprecip(ifull)
      common /soilbal/ trnoff(ifull)
      common /soilbal/ totenbal(ifull)
      common /soilbal/ osnowd0(ifull) ! snow depth at time step 0
      common /soilbal/ wbtot0(ifull)  ! total soil water at time 0 
