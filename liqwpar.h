      integer, parameter :: ifullw=1
c     if really doing liquid water variables use the following instead:
c     integer, parameter :: ifullw=ifull+iextra
      real, dimension(ifullw,kl) :: qfg, qlg
      common/liqw/qfg,qlg
