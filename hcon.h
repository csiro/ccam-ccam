           real amolwt,csubp,diffctr,o3difctr,p0,p0xzp2,p0xzp8,
     &     p0x2,radcon,rgas,rgassp,secpda,ratco2mw,rath2omw,radcon1,
     &     ginv,p0inv,gp0inv 
      real hundred,hninety,sixty,fifty,ten,eight,five, 
     &     four,three,two,one,haf,quartr,zero
      real h83e26,h71e26,h1e15,h1e13,h1e11,h1e8,h4e5,
     &     h165e5,h5725e4,h488e4,h1e4,h24e3,h20788e3,
     &     h2075e3,h1224e3,h5e2,h3082e2,h3e2,h2945e2,
     &     h23e2,h15e2,h35e1,h3p6,h181e1,h18e1,h2p9,h2p8,
     &     h2p5,h1p8,h1p4387,h1p4,h1p25892,hp8,hp518,
     &     hp369,hp1 
      real h44871m2,h559m3,h1m3,h987m4,h285m4,h1m4,
     &     h6938m5,h394m5,h37412m5,h1439m5,h128m5,h1m5,
     &     h7m6,h4999m6,h25452m6,h1m6,h391m7,h1174m7,
     &     h8725m8,h327m8,h257m8,h1m8,h23m10,h14m10, 
     &     h11m10,h1m10,h83m11,h82m11,h8m11,h77m11,
     &     h72m11,h53m11,h48m11,h44m11,h42m11,h37m11,
     &     h35m11,h32m11,h3m11,h28m11,h24m11,h23m11, 
     &     h2m11,h18m11,h15m11,h14m11,h114m11,h11m11,
     &     h1m11,h96m12,h93m12,h77m12,h74m12,h65m12, 
     &     h62m12,h6m12,h45m12,h44m12,h4m12,h38m12,
     &     h37m12,h3m12,h29m12,h28m12,h24m12,h21m12, 
     &     h16m12,h14m12,h12m12,h8m13,h46m13,h36m13, 
     &     h135m13,h12m13,h1m13,h3m14,h15m14,h14m14, 
     &     h1m17,h1m18,h1m19,h1m20,h1m21,h1m22,h1m23,
     &     h1m24,h26m30,h14m30,h25m31,h21m31,h12m31, 
     &     h9m32,h55m32,h45m32,h4m33,h62m34,h1m60
      real hmp575,hm13ez,hm19ez,hm1e1,hm181e1,hm1e2
      real h1e6,h2e6,h1m2,hmp66667,hm6666m2,hp166666,
     &     h41666m2,hmp5,hm2m2,h29316e2,h1226e1,h3116e1, 
     &     h9p94,hp6,h625m2,hp228,hp60241,hm1797e1,
     &     h8121e1,h2e2,hm1ez,h26e2,h44194m2,h1p41819
      real hp219,hp144,hp816,h69766e5,h235m3,hp26, 
     &     h129m2,h75826m4,h1p082,hp805,h1386e2, 
     &     h658m2,h1036e2,h2118m2,h42m2,h323m4,
     &     h67390e2,hp3795,hp5048,h102m5,h451m6
      real h16e1,hm161e1,h161e1,h3m3,h101m16,
     &     hm1597e1,h25e2,hp118666,h15m5,h3p5,h18e3, 
     &     h6p08108,hmp805,hp602409,hp526315,
     &     h28571m2,h1m16
      real h3m4
      real hm8e1 
      real h28e1 

      common/phycon/amolwt,csubp,diffctr,o3difctr,p0,
     &            p0xzp2,p0xzp8,p0x2,radcon,rgas,rgassp,secpda
      common/phycon/ratco2mw,rath2omw 
      common/phycon/radcon1 
      common/phycon/ginv,p0inv,gp0inv 
      common/hcon/hundred,hninety,sixty,fifty,ten,eight,five, 
     &            four,three,two,one,haf,quartr,zero
      common/hcon/h83e26,h71e26,h1e15,h1e13,h1e11,h1e8,h4e5,
     &            h165e5,h5725e4,h488e4,h1e4,h24e3,h20788e3,
     &            h2075e3,h1224e3,h5e2,h3082e2,h3e2,h2945e2,
     &            h23e2,h15e2,h35e1,h3p6,h181e1,h18e1,h2p9,h2p8,
     &            h2p5,h1p8,h1p4387,h1p4,h1p25892,hp8,hp518,
     &            hp369,hp1 
      common/hcon/h44871m2,h559m3,h1m3,h987m4,h285m4,h1m4,
     &            h6938m5,h394m5,h37412m5,h1439m5,h128m5,h1m5,
     &            h7m6,h4999m6,h25452m6,h1m6,h391m7,h1174m7,
     &            h8725m8,h327m8,h257m8,h1m8,h23m10,h14m10, 
     &            h11m10,h1m10,h83m11,h82m11,h8m11,h77m11,
     &            h72m11,h53m11,h48m11,h44m11,h42m11,h37m11,
     &            h35m11,h32m11,h3m11,h28m11,h24m11,h23m11, 
     &            h2m11,h18m11,h15m11,h14m11,h114m11,h11m11,
     &            h1m11,h96m12,h93m12,h77m12,h74m12,h65m12, 
     &            h62m12,h6m12,h45m12,h44m12,h4m12,h38m12,
     &            h37m12,h3m12,h29m12,h28m12,h24m12,h21m12, 
     &            h16m12,h14m12,h12m12,h8m13,h46m13,h36m13, 
     &            h135m13,h12m13,h1m13,h3m14,h15m14,h14m14, 
     &            h1m17,h1m18,h1m19,h1m20,h1m21,h1m22,h1m23,
     &            h1m24,h26m30,h14m30,h25m31,h21m31,h12m31, 
     &            h9m32,h55m32,h45m32,h4m33,h62m34,h1m60
      common/hcon/hmp575,hm13ez,hm19ez,hm1e1,hm181e1,hm1e2
      common/hcon/h1e6,h2e6,h1m2,hmp66667,hm6666m2,hp166666,
     &            h41666m2,hmp5,hm2m2,h29316e2,h1226e1,h3116e1, 
     &            h9p94,hp6,h625m2,hp228,hp60241,hm1797e1,
     &            h8121e1,h2e2,hm1ez,h26e2,h44194m2,h1p41819
      common/hcon/hp219,hp144,hp816,h69766e5,h235m3,hp26, 
     &            h129m2,h75826m4,h1p082,hp805,h1386e2, 
     &            h658m2,h1036e2,h2118m2,h42m2,h323m4,
     &            h67390e2,hp3795,hp5048,h102m5,h451m6
      common/hcon/h16e1,hm161e1,h161e1,h3m3,h101m16,
     &            hm1597e1,h25e2,hp118666,h15m5,h3p5,h18e3, 
     &            h6p08108,hmp805,hp602409,hp526315,
     &            h28571m2,h1m16
      common/hcon/h3m4
      common/hcon/hm8e1 
      common/hcon/h28e1 

