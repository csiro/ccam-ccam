!     vertical [advection] options (globpe, vadvtvd,vadv30,vadv30in)
      integer nvad, nvadh
      common/paramvrt/nvad,nvadh

      integer nimp, nthub, ntvd
!     for RMIP1 following were (1,1,2); during 2002-3 were (0,2,3)
      parameter (nimp=0)  !  0 for original explicit non-flux TVD term (usual)
!                            1 for implicit non-flux TVD term
      parameter (nthub=2) !  1 equivalent to original TVD higher-order
!                            2 higher-order is Lax-Wendroff TVD (usual)
      parameter (ntvd=3)  !  0 original phi flux-limiter
!                            1 redefined original phi (sames answers as 0)
!                            2 MC phi flux-limiter
!                            3 superbee flux-limiter (usual)
