      real ax(ifull),bx(ifull),ay(ifull),by(ifull)                       &
     &             ,az(ifull),bz(ifull)
      common/vecsuv_g/ax,bx,ay,by,az,bz
      real  ax6(il,il,0:5),ay6(il,il,0:5),az6(il,il,0:5)                 &
     &     ,bx6(il,il,0:5),by6(il,il,0:5),bz6(il,il,0:5)
      equivalence (ax6,ax),(ay6,ay),(az6,az)                             &
     &           ,(bx6,bx),(by6,by),(bz6,bz)
