      real ax, bx, ay, by, az, bz
      common/vecsuv/ax(ifull),bx(ifull),ay(ifull),by(ifull),  ! re-ordered
     &              az(ifull),bz(ifull)    !  ,cosa(ifull)    ! for globpea
      real ax6(il,il,0:5),ay6(il,il,0:5),az6(il,il,0:5),
     &     bx6(il,il,0:5),by6(il,il,0:5),bz6(il,il,0:5)
      equivalence (ax6,ax),(ay6,ay),(az6,az),
     &            (bx6,bx),(by6,by),(bz6,bz)
