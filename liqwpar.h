      integer ifullw
      parameter(ifullw=1)
c     if really doing liquid water variables use the following instead:
c     parameter(ifullw=ifull)
      real qfg,qlg
      common/liqw/qfg(ifullw,kl),qlg(ifullw,kl)
