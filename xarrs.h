      real, dimension(ifull,kl) :: ux, vx
      ! tx and pslx have iextra because they're used in ints call.
      real, dimension(ifull+iextra,kl) :: tx, pslx
      common /xarrs/ ux,vx,tx,pslx
