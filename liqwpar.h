!      integer, parameter :: ifullw=1  being removed from 29/1/04
!     if really doing liquid water variables use the following instead:
      integer, parameter :: ifullw=ifull+iextra
      real, dimension(ifullw,kl) :: qfg, qlg
      common/liqw/qfg,qlg
