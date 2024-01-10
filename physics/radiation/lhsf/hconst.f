      subroutine hconst 
c 
      include 'hcon.h'
c     subroutine hconst defines variables to represent floating-
c       point constants.
c 
c     comdeck hcon contains the common block for these floating-
c       point constants.
c 
c     the naming conventions for the floating-point variables are 
c       as follows: 
c 
c   1) physical and mathematical constants will be given names
c     relevant to their meaning 
c   2) other constants will be given names relevant to their value
c      and adhering to the following conventions: 
c       a) the first letter will be represented with an 'h' except
c          for i) and j) below
c       b) a decimal point will be represented with a 'p' 
c       c) there will be no embedded '0'(zero); all 0's will
c          be represented with a 'z'
c       d) a minus sign will be represented with an 'm' 
c       e) the decimal point is assumed after the first digit for 
c          numbers with exponents 
c       f) positive exponents are indicated with 'e';negative 
c          exponents with 'm' 
c       g) digits are truncated in order to have no more than 8 
c          characters per name
c       h) numbers less than 0.1 and greater than 10. will be 
c          represented in exponent format (except a few special cases)
c       i) the whole numbers from 0.0 through 10.,and 20.,30.,40.,50.,
c          60.,70.,80.,90.,100.,will be spelled out 
c       j) good judgment will prevail over all conventions
c 
c       examples
c     constant           variable name             convention 
c      600.                 lheatc                  1)
c      680.                 lheats                  1)
c     1.4142               sqroot2                  1)
c     2.0                    two                    2)-(i)
c    -3.0                  hm3pz                    2)-(a,b,d)
c    310.                  c31e2                    2)-(a,e,f,h)
c   -0.7239e-9             hm723m1z                 2)-(a,c,d,e,f,g,h)
c     0.0                   zero                    2)-(i)
c     0.1                   hp1                     2)-(a,b,h)
c     0.01                 h1m2                     2)-(a,e,f,h)
c     30.                  thirty                   2)-(h,i)
c     0.5                  haf                      2)-(j)
c     9.0                  hnine                    2)-(j)
c
c******the following are physical constants*****
c        arranged in alphabetical order

!!      include 'PHYSPARAMS.h'

      amolwt=28.9644
      csubp=1.00484e7 
      diffctr=1.66
!!      localg=grav*100 !grav=9.806 set in PHYSPARAMS.h  -- not for globpe
      ginv=1./980.616
      o3difctr=1.90 
      p0=1013250. 
      p0inv=1./p0 
      gp0inv=ginv*p0inv 
      p0xzp2=202649.902 
      p0xzp8=810600.098 
      p0x2=2.*1013250.
      radcon=8.427
      radcon1=1./8.427
      ratco2mw=1.519449738
      rath2omw=.622 
      rgas=8.3142e7 
      rgassp=8.31432e7
      secpda=8.64e4 
c 
c******the following are mathematical constants*******
c        arranged in decreasing order 
      hundred=100.
      hninety=90. 
      sixty=60. 
      fifty=50. 
      ten=10. 
      eight=8.
      five=5. 
      four=4. 
      three=3.
      two=2.
      one=1.
      haf=0.5 
      quartr=0.25 
      zero=0. 
c 
c******following are positive floating point constants(h's) 
c       arranged in decreasing order
      h83e26=8.3e26 
      h71e26=7.1e26 
      h1e15=1.e15 
      h1e13=1.e13 
      h1e11=1.e11 
      h1e8=1.e8 
      h2e6=2.0e6
      h1e6=1.0e6
      h69766e5=6.97667e5
      h4e5=4.e5 
      h165e5=1.65e5 
      h5725e4=57250.
      h488e4=48800. 
      h1e4=1.e4 
      h24e3=2400. 
      h20788e3=2078.8 
      h2075e3=2075. 
      h18e3=1800. 
      h1224e3=1224. 
      h67390e2=673.9057 
      h5e2=500. 
      h3082e2=308.2 
      h3e2=300. 
      h2945e2=294.5 
      h29316e2=293.16 
      h26e2=260.0 
      h25e2=250.
      h23e2=230.
      h2e2=200.0
      h15e2=150.
      h1386e2=138.6 
      h1036e2=103.6 
      h8121e1=81.21 
      h35e1=35. 
      h3116e1=31.16 
      h28e1=28. 
      h181e1=18.1 
      h18e1=18. 
      h161e1=16.1 
      h16e1=16. 
      h1226e1=12.26 
      h9p94=9.94
      h6p08108=6.081081081
      h3p6=3.6
      h3p5=3.5
      h2p9=2.9
      h2p8=2.8
      h2p5=2.5
      h1p8=1.8
      h1p4387=1.4387
      h1p41819=1.418191 
      h1p4=1.4
      h1p25892=1.258925411
      h1p082=1.082
      hp816=0.816 
      hp805=0.805 
      hp8=0.8 
      hp60241=0.60241 
      hp602409=0.60240964 
      hp6=0.6 
      hp526315=0.52631579 
      hp518=0.518 
      hp5048=0.5048 
      hp3795=0.3795 
      hp369=0.369 
      hp26=0.26 
      hp228=0.228 
      hp219=0.219 
      hp166666=.166666
      hp144=0.144 
      hp118666=0.118666192
      hp1=0.1 
c        (negative exponentials begin here) 
      h658m2=0.0658 
      h625m2=0.0625 
      h44871m2=4.4871e-2
      h44194m2=.044194
      h42m2=0.042 
      h41666m2=0.0416666
      h28571m2=.02857142857 
      h2118m2=0.02118 
      h129m2=0.0129 
      h1m2=.01
      h559m3=5.59e-3
      h3m3=0.003
      h235m3=2.35e-3
      h1m3=1.0e-3 
      h987m4=9.87e-4
      h323m4=0.000323 
      h3m4=0.0003 
      h285m4=2.85e-4
      h1m4=0.0001 
      h75826m4=7.58265e-4 
      h6938m5=6.938e-5
      h394m5=3.94e-5
      h37412m5=3.7412e-5
      h15m5=1.5e-5
      h1439m5=1.439e-5
      h128m5=1.28e-5
      h102m5=1.02e-5
      h1m5=1.0e-5 
      h7m6=7.e-6
      h4999m6=4.999e-6
      h451m6=4.51e-6
      h25452m6=2.5452e-6
      h1m6=1.e-6
      h391m7=3.91e-7
      h1174m7=1.174e-7
      h8725m8=8.725e-8
      h327m8=3.27e-8
      h257m8=2.57e-8
      h1m8=1.0e-8 
      h23m10=2.3e-10
      h14m10=1.4e-10
      h11m10=1.1e-10
      h1m10=1.e-10
      h83m11=8.3e-11
      h82m11=8.2e-11
      h8m11=8.e-11
      h77m11=7.7e-11
      h72m11=7.2e-11
      h53m11=5.3e-11
      h48m11=4.8e-11
      h44m11=4.4e-11
      h42m11=4.2e-11
      h37m11=3.7e-11
      h35m11=3.5e-11
      h32m11=3.2e-11
      h3m11=3.0e-11 
      h28m11=2.8e-11
      h24m11=2.4e-11
      h23m11=2.3e-11
      h2m11=2.e-11
      h18m11=1.8e-11
      h15m11=1.5e-11
      h14m11=1.4e-11
      h114m11=1.14e-11
      h11m11=1.1e-11
      h1m11=1.e-11
      h96m12=9.6e-12
      h93m12=9.3e-12
      h77m12=7.7e-12
      h74m12=7.4e-12
      h65m12=6.5e-12
      h62m12=6.2e-12
      h6m12=6.e-12
      h45m12=4.5e-12
      h44m12=4.4e-12
      h4m12=4.e-12
      h38m12=3.8e-12
      h37m12=3.7e-12
      h3m12=3.e-12
      h29m12=2.9e-12
      h28m12=2.8e-12
      h24m12=2.4e-12
      h21m12=2.1e-12
      h16m12=1.6e-12
      h14m12=1.4e-12
      h12m12=1.2e-12
      h8m13=8.e-13
      h46m13=4.6e-13
      h36m13=3.6e-13
      h135m13=1.35e-13
      h12m13=1.2e-13
      h1m13=1.e-13
      h3m14=3.e-14
      h15m14=1.5e-14
      h14m14=1.4e-14
      h101m16=1.01e-16
      h1m16=1.0e-16 
      h1m17=1.e-17
      h1m18=1.e-18
      h1m19=1.e-19
      h1m20=1.e-20
      h1m21=1.e-21
      h1m22=1.e-22
      h1m23=1.e-23
      h1m24=1.e-24
      h26m30=2.6e-30
      h14m30=1.4e-30
      h25m31=2.5e-31
      h21m31=2.1e-31
      h12m31=1.2e-31
      h9m32=9.e-32
      h55m32=5.5e-32
      h45m32=4.5e-32
      h4m33=4.e-33
      h62m34=6.2e-34
!!      h1m60=1.0e-60 
c 
c******following are negative floating point constants (hm's) 
c          arranged in descending order 
      hm2m2=-.02
      hm6666m2=-.066667 
      hmp5=-0.5 
      hmp575=-0.575 
      hmp66667=-.66667
      hmp805=-0.805 
      hm1ez=-1. 
      hm13ez=-1.3 
      hm19ez=-1.9 
      hm1e1=-10.
      hm1597e1=-15.97469413 
      hm161e1=-16.1 
      hm1797e1=-17.97469413 
      hm181e1=-18.1 
      hm8e1=-80.
      hm1e2=-100. 
c 
      return
      end 

