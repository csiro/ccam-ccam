!     vertical [advection] options (globpe, vadvtvd,vadv30,vadv30in)
      common/paramvrt/nvad,nvadh

!     for RMIP1 following were (1,1,2); during 2002 were (0,2,3)
      parameter (nimp=0)  !  0 for original explicit non-flux TVD term
!                            1 for implicit non-flux TVD term
      parameter (nthub=2) !  0 original TVD higher-order
!                            1 equivalent to original TVD higher-order
!                            2 higher-order is Lax-Wendroff TVD
      parameter (ntvd=3)  !  0 original phi flux-limiter
!                            1 redefined original phi (sames answers as 0)
!                            2 MC phi flux-limiter
!                            3 superbee flux-limiter

!     parameter (n_bs=0)  !  0 off            was  in vadv30in, vadvl_w
!                            1 for Sun-Yeh limiter           (not for spline)
!                            2 for Bermejo-Staniforth limiter
!                            3 all B-S, except S-Y for T     (not for spline)
!                            only used in vadv30in, vadvl_w (not vadvtvd)
!                            N.B. Gadd/L_W has different meaning for n_bs


