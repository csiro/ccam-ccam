! eak version 16/03/06
!     tggsn needs to directly precede tgg
      real tggsn, tgg, wb, wbice, wbtot, smass, ssdn, ssdnn, snowd 
      real  osnowd, snage, sno, gflux, sgflux, snowflx, otgsoil 
      real  runoff, rnof1, rnof2, rtsoil, albvisnir, albsoil, albsoilsn
      real  fracice, sicedep
      integer isflag

      common /soilcom/ tggsn(ifull,3) ! snow temperatures
      common /soilcom/ tgg(ifull,ms)  ! soil temperature (in K)
      common /soilcom/ wb(ifull,ms)   ! soil moisture (volumetric)
      common /soilcom/ wbice(ifull,ms)! soil ice
      common /soilcom/ wbtot(ifull)   ! total soil water (mm)
      common /soilcom/ smass(ifull,3) ! snow masses
      common /soilcom/ ssdn(ifull,3)  ! snow densities
      common /soilcom/ ssdnn(ifull)   ! average snow density
      common /soilcom/ snowd(ifull)   ! snow depth (liquid water)  
      common /soilcom/ osnowd(ifull)  ! snow depth from previous time step
      common /soilcom/ snage(ifull)   ! snow age
      common /soilcom/ sno(ifull)     ! accum. snow in mm since last write 
      common /soilcom/ isflag(ifull)  !
      common /soilcom/ gflux(ifull)   !
      common /soilcom/ sgflux(ifull)  !
      common /soilcom/ snowflx(ifull) ! surface snow melt heat flux
      common /soilcom/ otgsoil(ifull) ! soil or snow surface temp. at the previous dt
      common /soilcom/ runoff(ifull)  ! runoff (mm)
      common /soilcom/ rnof1(ifull)   ! surface runoff (mm)
      common /soilcom/ rnof2(ifull)   ! deep drainage (mm)
      common /soilcom/ rtsoil(ifull)  ! 
      common /soilcom/ albsoil(ifull)      ! soil albedo
      common /soilcom/ albsoilsn(ifull,2)    ! soil+snow albedo
      common /soilcom/ albvisnir(ifull,2)  ! soil+snow+plant albedo
      common /seaice/  fracice(ifull) ! sea-ice fraction
      common /seaice/  sicedep(ifull) ! sea-ice depth


!      common /soilcom/ coszros(ifull) ! cos of zenith angle
!                                       variables for SEB and Hyd. balance
!      common /soilcom/ tevap(ifull)
!      common /soilcom/ tprecip(ifull)
!      common /soilcom/ trnoff(ifull)
!      common /soilcom/ totenbal(ifull)
!      common /soilcom/ osnowd0(ifull) ! snow depth at time step 0
!      common /soilcom/ wbtot0(ifull)  ! total soil water at time 0 
!!      common /soilcom/ zse(ms)        ! soil layers thickenss (m)
!!      common /soilcom/ zshh(ms+1)     ! 
!!      common /soilcom/ froot(13,ms)      ! root distribution 
