      real x, y, z, wts
      common/xyz/x(ifull),y(ifull),z(ifull),wts(ifull) ! rlat and rlong gone
      real x6(il,il,0:5),y6(il,il,0:5),z6(il,il,0:5)
      real x06(il,il,0:5),y06(il,il,0:5),z06(il,il,0:5)
      equivalence (x,x6,x06),(y,y6,y06),(z,z6,z06)
