!     include 'newmpar.h'
      parameter (la=25)

      common / eta / deta(kl), aeta(kl), eta(kl+1)

      common / masks / hbm2(ifull), htm(ifull,kl)
     .               , sm(ifull), klh(ifull), res(ifull)

      parameter ( itb=76, jtb=171 )
      common / phys / ktm, dtq2
     &              , pl,  rdq, rdth, rdp, rdthe
     &              , qs0(jtb), sqs(jtb), the0(itb), sthe(itb)
     &              , ptbl(itb,jtb), ttbl(jtb,itb)
      real thl
      data thl/210./

      common /vrbls/ pd(ifull)
      common /pvrbls/ prec(ifull), cuprec(ifull), cldefi(ifull)
      common/work3e/t(ifull,kl),q(ifull,kl)  ! 28/9/00   jlm

