!     dynamics options (globpe, adjust5, nonlin, upglobal)

!     parameter (mfix_qg=1)   ! 1 "mass" fix for qg
!                               2 "mass" fix for qg and trace gases

!     mfix in namelist:        -1 on pslx in upglobal
!                               0 off
!                               1 cunning in adjust5
!                               2 more-cunning in adjust5
!
!     morder was in namelist, e.g. 24 ! set morder to 2 or 4 for nonlin
!        4 till Mon  08-31-1998, then 2; 4 again on Tue  09-01-1998
!        currently favouring 24
!
!     morder_r was in namelist ! set morder_r to 2 or 4 for upglobal
!     this gives order for nritch advection options - not used nowadays

      integer m,mex,mfix,mfix_qg,mspec,mup,nonl,nritch,nritch_t,nrot,      &
     &        nstag,nstagu,ntbar,nuvfilt,nvsplit,nxmap,precon,nh,npex
      real    epsp,epsu,epsf,restol
      common/paramdyn/epsp,epsu,epsf,m,mex,mfix,mfix_qg,                   &
     &                mspec,mup,nonl,nritch,nritch_t,nrot,                 &
     &                nstag,nstagu,ntbar,nuvfilt,                          &
     &                nvsplit,nxmap,restol,precon,                         &
     &                nh,npex      ! for future use
!              (ntbar=0)           ! 0 for standard
!              (ntbar=(kl+1)/2)    ! level# for tbar2d with T set in nonlin
!           nspec_us set to 1 by nritch=102 (unstaggering at opposite ktau)

!                 nvsplit    0  uses tendencies for vadv, radn & vertmix
!                            1  splits radn, vertmix, gwdrag, conjob (not vadv)
!                            2  splits radn, vertmix, gwdrag, conjob & vadv
!                            3  splits just vadv
!                           -1  splits just vertmix (always done in DARLAM)
!                          N.B. qg always split for vadv
!                          N.B. always split for vadv called from adjust5
!                          N.B. nkuo=45 also splits radn & vertmix, but
!                               usually uses tn tendencies from conjob
