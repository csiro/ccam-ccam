      integer il, npanels, jl, kl, ksl, ms, ifull, ij, ijk, iquad
!     plan to replace ij by ifull eventually
c     parameter(il=32,npanels=13,jl=il+npanels*il,kl=18,ksl=3,ms=6)
      parameter(il=48,npanels=5,jl=il+npanels*il,kl=18,ksl=3,ms=6)
!      parameter(il=48,npanels=5,jl=il+npanels*il,kl=18,ksl=3,ms=6)
!      parameter(il=20,npanels=5,jl=il+npanels*il,kl=9,ksl=3,ms=6)
      parameter(ifull=il*jl,ij=il*jl,ijk=il*jl*kl)
      parameter( iquad=1+il*((8*npanels)/(npanels+4)) )
!     for     npanels:   0          5        13
!                  jl:   -         6*il     14*il
!                quad:   1         4*il+1   6*il+1
