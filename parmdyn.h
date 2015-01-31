!     dynamics options (globpe, adjust5, nonlin, upglobal)

!     parameter (mfix_qg=1)   ! 1 "mass" fix for qg
!                               2 "mass" fix for qg and trace gases

!     mfix in namelist:        -1 on pslx in upglobal
!                               0 off
!                               1 cunning in adjust5
!                               2 more-cunning in adjust5

      integer         mex,mfix,mfix_qg,mspec,mup,mfix_tr
      integer         nh,nonl,nritch_t,nrot,mfix_aero
      integer         nstag,nstagu,ntbar,nvsplit,nxmap,precon,helmmeth
      integer         nstagoff
      real            epsp,epsu,epsf,restol
      
      common/paramdyn/mex,mfix,mfix_qg,mspec,mup,                        &
     &                nh,nonl,nritch_t,nrot,                             &
     &                nstag,nstagu,nstagoff,ntbar,nvsplit,nxmap,precon,  &
     &                helmmeth,epsp,epsu,epsf,restol,mfix_tr,mfix_aero

!            (ntbar=0)           ! 0 for standard
!            (ntbar=(kl+1)/2)    ! level# for tbar2d with T set in nonlin

!            nvsplit    0  uses tendencies for vadv, radn & vertmix
!                       1  splits radn, vertmix, gwdrag, conjob (not vadv)
!                       2  splits radn, vertmix, gwdrag, conjob & vadv
!                       3  splits just vadv
!                      -1  splits just vertmix 
!                      N.B. qg always split for vadv
!                      N.B. always split for vadv called from adjust5
