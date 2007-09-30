!     dynamics options (globpe, adjust5, nonlin, upglobal)

!     parameter (mfix_qg=1)   ! 1 "mass" fix for qg
!                               2 "mass" fix for qg and trace gases

!     mfix in namelist:        -1 on pslx in upglobal
!                               0 off
!                               1 cunning in adjust5
!                               2 more-cunning in adjust5

      integer         m,mex,mfix,mfix_qg,mspec,mup
      integer         nh,nonl,npex,nritch_t,nrot
      integer         nstag,nstagu,ntbar,nvsplit,nxmap,precon
      real            epsp,epsu,epsf,restol
      
      common/paramdyn/m,mex,mfix,mfix_qg,mspec,mup,                      &
     &                nh,nonl,npex,nritch_t,nrot,                        &
     &                nstag,nstagu,ntbar,nvsplit,nxmap,precon,           &
     &                epsp,epsu,epsf,restol

!            (ntbar=0)           ! 0 for standard
!            (ntbar=(kl+1)/2)    ! level# for tbar2d with T set in nonlin

!            nvsplit    0  uses tendencies for vadv, radn & vertmix
!                       1  splits radn, vertmix, gwdrag, conjob (not vadv)
!                       2  splits radn, vertmix, gwdrag, conjob & vadv
!                       3  splits just vadv
!                      -1  splits just vertmix 
!                      N.B. qg always split for vadv
!                      N.B. always split for vadv called from adjust5
