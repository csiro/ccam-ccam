c     conformal octagon - code of Jim Purser
c     JMcG usage:
c            1) calculate xx4, yy4, zz4 just for NE corner of regular octagon 1
c            2) in setxyz apply to all 14 panels
c            3) in setxyz convert panel values via schmidt
c            4) in setxyz save only xxy4, yy4 (after mapping via schmidta)
      subroutine jimco
      include 'newmpar.h'
      parameter(n=128)
      include 'bigxy4.h' ! common/bigxy4/xx4(iquad,iquad),yy4(iquad,iquad)
c     common/cstoct/a1(40),a2(40),b1(40),b2(40),m1g,mag,a,b,cx,cy
c    *,ssr,n1,rotm(3,3),sc,sci,zac,st,sti,ci,cir,ciq,ciap
      common/work2/zz4(iquad,iquad),work(1200)     ! to agree with setxyz
     .    ,dum2(18*il*jl - iquad*iquad -1200 )
      dimension xe(3),xep(3),dxedx(3),dxedy(3)
     *,gxedx(3),gxedy(3),hxedx(3),hxedy(3)
      dimension ra(0:n),qa(0:n)
     *,jumble(n)

c     call inoct(n1grid,nagrid,n1,clon,clat,sc,ra,qa,n,jumble,work)
c                 24      8    20  0.  90.  1.      128
      call inoct(6*il  ,2*il  ,20, 0. ,90. ,1.,ra,qa,128,jumble,work)

c     print *,'il,jl,iquad,arrsize ',il,jl,iquad,
c    .     18*il*jl - iquad*iquad -1200
      do j=1,iquad
c     do j=45,46
       y=(j-1)/real(iquad-1)
       do i=1,iquad
c      do i=45,46
        x=(i-1)/real(iquad-1)
c       print *,'i,j,n1,iquad,x,y ',i,j,n1,iquad,x,y
        call otoc1(x,y,xe)
c       jlm chooses just a positive quadrant
c       print *,'xea,n1               ',xe,n1
        xx4(i,j)=xe(1)
c       print *,'xeb,n1               ',xe,n1
        yy4(i,j)=xe(2)
c       print *,'xec,n1               ',xe,n1
        zz4(i,j)=xe(3)
c       print *,'xed,n1               ',xe,n1
       enddo
      enddo
      return
      end

c------------------------------------------------------------------------------c
c   r.j.purser, national meteorological center, washington d.c.  1995          c
c                   subroutine inoct                                           c
c   initialize the constants of common/cstoct/ used to describe grid and to    c
c   perform transformations between map coordinates and the sphere or disk     c
c                                                                              c
c --> n1grid    number of grid-spaces from center to side of octagon           c
c --> nagrid    number of grid-spaces from center-line to nearest vertex       c
c --> ntay      number of taylor series coefficients retained (ntay<41)        c
c --> flon0     longitude of center of octagon-1                               c
c --> flat0     latitude of center of octagon-1                                c
c --> scen      scale-enhancement parameter = radius of image of octagon-1 in  c
c               a stereographic projection centered on (flat0,flon0) that maps c
c               the concentric hemisphere to the unit disk                     c
c --- ra,qa     work arrays of dimension n+1 (n defined below)                 c
c --> n         number of the form 2**m (eg 128) and > 3*ntay defining the     c
c               period of fft used in the construction of taylor series        c
c --- jumble,work   work arrays of dimension n                                 c
c------------------------------------------------------------------------------c
      subroutine inoct(n1grid,nagrid,ntay,flon0,flat0,scen
     *,ra,qa,n,jumble,work)
      complex ci,cir,ciq,ciap
      common/cstoct/a1(40),a2(40),b1(40),b2(40),m1g,mag,a,b,cx,cy
     *,ssr,n1,rotm(3,3),sc,sci,zac,st,sti,ci,cir,ciq,ciap
      dimension ra(0:*),qa(0:*),jumble(n),work(n)
      n1=ntay
      if(3*n1.gt.n)
     *stop'n too small in inoct, choose a larger power of 2'
      m1g=n1grid ! grid intervals between center and side
      mag=nagrid ! grid intervals between center and first corner
      sc=scen
      sci=1./sc
      zac=(sci-sc)/(sci+sc)    !   (1-c**2)/(1+c**2)  jlm
      a=float(mag)/m1g ! nice ratios: 1/3, 3/7, 7/17
c     print *,'calling indoct, n1= ',n1
      call indoct(ra,qa,n,jumble,work)
c     print *,'calling ininmap'
      call ininmap(flon0,flat0,rotm)
      return
      end

      subroutine indoct(ra,qa,n,jumble,work)
      complex ci,cir,ciq,cia,ciap,cirs2,z,w,cang,cdang1,cdang2
      common/cstoct/a1(40),a2(40),b1(40),b2(40),m1g,mag,a,b,cx,cy
     *,ssr,n1,rotm(3,3),sc,sci,zac,st,sti,ci,cir,ciq,ciap
      dimension ra(0:*),qa(0:*),jumble(n),work(n)
      st=4./3.
      sti=3./4.
      ci=cmplx(0.,1.)
      cia=ci*a
      ciap=cia+1.
      cir=(-ci)**.75
      ciq=ci**.25
      nh=n/2
      r2=sqrt(2.)
      orc=.5            ! relaxation factor
      sw1=sqrt(1.+a*a)
      sw2=amin1(a*2,r2*(1.-a))
      ssr=(sw1/sw2)**2
      rs1=.95*sw1
      rs2=.95*sw2
      rs3=rs1**4
      rs4=rs2**st

      do i=1,n1
       a1(i)=0.
      enddo

c  first guess for first few taylor coefficients (based on a=1/3)
c  needed to begin the boot-strap procedure:
      a1(1)= 1.05060
      a1(2)= -.172629
      a1(3)=  .125716
      a1(4)=  .00015657
      a1(5)= -.0017866
      a1(6)= -.0031168

      do i=1,n1
       a2(i)=0.
      enddo
      a2(1)= .688411
      a2(2)=-.180924
      a2(3)=-.186593
      a2(4)= .024811
      a2(5)=-.044888
      a2(6)=-.049190
      b=-.147179

      jumble(1)=0
      pio4=atan(1.)
      cdang1=ci*pio4/nh
      cdang2=-3.*cdang1
      cirs2=ci*rs2

      do kit=1,50
       do i=0,nh
        cang=cdang1*i
        z=rs1*cexp(cang)
        call toct(z,w)
        w=w**4
        u=real(w)
        v=aimag(w)
        ra(i)=u
        qa(i)=-v
        ra(n-i)=u
        qa(n-i)=v
       enddo
       qa(0)=0.
       qa(nh)=0.
       call cfft(ra,qa,n,1.,work,jumble)
       rs3p=1.
       do i=1,n1
        rs3p=rs3p*rs3
        a1(i)=a1(i) +orc*(ra(i)/rs3p-a1(i))
       enddo

      do i=0,nh
       cang=cdang2*i
       z=ciap-cirs2*cexp(cang)
       call toct(z,w)
       w=ci*(w-1.)/(w+1.)
       u=real(w)
       v=aimag(w)
       ra(i)=u
       qa(i)=v
       ra(n-i)=u
       qa(n-i)=-v
      enddo
      qa(0)=0.
      qa(nh)=0.
      call cfft(ra,qa,n,1.,work,jumble)
      b=b+orc*(ra(0)-b)
      rs4p=1.
      do i=1,n1
       rs4p=rs4p*rs4
       a2(i)=a2(i) +orc*(ra(i)/rs4p-a2(i))
      enddo

      enddo
      w=b
      w=(w+ci)/(ci-w)
      cx=real(w)
      cy=aimag(w)
      call cinvrt(a1,b1,n1,ra)
      call cinvrt(a2,b2,n1,ra)

c     print'('' image of octagon-vertex on unit-circle:'')'
c     write(6,62)cx,cy
c     print'('' taylor coefficients, series a1,a2,b1,b2:'')'
c     write(6,63)b
c     do i=1,n1
c      write(6,64)i,a1(i),a2(i),b1(i),b2(i)
c     enddo
62    format(2(1x,e12.6))
63    format(4x,'0',13x,4(1x,e12.6))
64    format(1x,i4,4(1x,e12.6))
      return
      end

c------------------------------------------------------------------------------c
c   r.j.purser, national meteorological center, washington d.c.  1995          c
c                   subroutine  otoc1                                           c
c   transform to earth-centered cartesians from octagon-1 (otoc1) or from      c
c   octagon-2 (otoc2).                                                         c
c                                                                              c
c --> x,y   x and y map-coordinate of point in octagon-1 (otoc1) or            c
c           octagon-2 (otoc2)                                                  c
c <-- xe    3-vector of earth-centered cartesian coordinates corresponding     c
c           to map location (x,y).                                             c
c------------------------------------------------------------------------------c
      subroutine otoc1(x,y,xe)
      complex ci,cir,ciq,z,w,ciap
      common/cstoct/a1(40),a2(40),b1(40),b2(40),m1g,mag,a,b,cx,cy
     *,ssr,n1,rotm(3,3),sc,sci,zac,st,sti,ci,cir,ciq,ciap
      dimension xe(3)
      z=cmplx(x,y)
      call toct(z,w)
      w=w*sc
      xa=real(w)
      ya=aimag(w)
      xx=xa**2+ya**2
      s=2./(1.+xx)
      xe(1)=s-1.
      xe(2)=s*xa
      xe(3)=s*ya
      call nmapt(rotm,xe,xe)
      return

      entry otoc2(x,y,xe)
      z=cmplx(-x,y)
      call toct(z,w)
      w=w*sci
      xa=real(w)
      ya=aimag(w)
      xx=xa**2+ya**2
      s=-2./(1.+xx)
      xe(1)=1.+s
      xe(2)=s*xa
      xe(3)=s*ya
      call nmapt(rotm,xe,xe)
      return
      end

c------------------------------------------------------------------------------c
c   r.j.purser, national meteorological center, washington d.c.  1995          c
c                   subroutine  ctoo                                            c
c   transform from earth-centered cartesians to octagon-1 (kmap=1) or to       c
c   octagon-2 (kmap=-1)                                                        c
c                                                                              c
c --> xe    earth-centered coordinates (3-vector) of a point                   c
c <-- x,y   map-coordinates of this point                                      c
c <-- kmap  map indicator (=1 for octagon-1, =-1 for octagon-2)                c
c------------------------------------------------------------------------------c
      subroutine ctoo(xe,x,y,kmap)
      logical kx,ky,kxy,ks
      complex ci,cir,ciq,z,w,ciap
      common/cstoct/a1(40),a2(40),b1(40),b2(40),m1g,mag,a,b,cx,cy
     *,ssr,n1,rotm(3,3),sc,sci,zac,st,sti,ci,cir,ciq,ciap
      dimension xe(3)
      call nmap(rotm,xe,xe)
      xa=xe(2)
      ya=xe(3)
      za=xe(1)
      if(za.gt.zac)then
       kmap=1
       zai=sci/(1.+za)
       x=xa*zai
       y=ya*zai
      else
       kmap=-1
       zai=sc/(1.-za)
       x=xa*zai
       y=-ya*zai
      endif
      kx=x.lt.0.
      if(kx)x=-x
      ky=y.lt.0.
      if(ky)y=-y
      kxy=y.gt.x
      if(kxy)then
       t=x
       x=y
       y=t
      endif
      dd1=x**2+y**2
      dd2=(cx-x)**2+(cy-y)**2
      w=cmplx(x,y)
      ks=dd1.lt.ssr*dd2
      if(ks)then
       w=w**4
       call tay(w,b1,n1,z)
       z=ciq*(-z*ci)**.25
      else
       w=ci*(w-1.)/(w+1.)
       w=w-b
       call tay(w,b2,n1,z)
       z=cir*(ci*z)**sti
       z=ciap-ci*z
      endif
      x=real(z)
      y=aimag(z)
      if(kxy)then
       t=x
       x=y
       y=t
      endif
      if(ky)y=-y
      if(kx)x=-x
      return
      end

c------------------------------------------------------------------------------c
c   r.j.purser, national meteorological center, washington d.c.  1995          c
c                   subroutine  toct                                            c
c   transform from complex-z in the standard unit-octagon to complex-w in the  c
c   unit-circle                                                                c
c------------------------------------------------------------------------------c
      subroutine toct(z,w)
      complex ci,cir,ciq,z,w,zt,ciap
      logical kx,ky,kxy,kx1,kxy1,ks
      common/cstoct/a1(40),a2(40),b1(40),b2(40),m1g,mag,a,b,cx,cy
     *,ssr,n1,rotm(3,3),sc,sci,zac,st,sti,ci,cir,ciq,ciap
      x=real(z)
      y=aimag(z)
      kx=x.lt.0.
      if(kx)x=-x
      ky=y.lt.0.
      if(ky)y=-y
      kxy=y.gt.x
      if(kxy)then
       t=x
       x=y
       y=t
      endif
      kx1=x.gt.1.
      kxy1=x+y.gt.1.+a
      if(kx1)then
       x=2.-x
      elseif(kxy1)then
       t=x
       x=1.+a-y
       y=1.+a-t
      endif
      dd1=x**2+y**2
      dd2=(1.-x)**2+(a-y)**2

      zt=cmplx(x,y)
      ks=dd1.lt.ssr*dd2
      if(ks)then
       zt=zt**4
       call tay(zt,a1,n1,w)
       w=ciq*(-w*ci)**.25
      else
c      jlm: note that DEC alpha has trouble with (0.,0.)**4/3
c      print *,'in toct,n1, st = ',n1,st
c      print *,'zt 1 ',zt
       zt=a+ci*zt-ci
c      print *,'a ',a
c      print *,'ci ',ci
c      print *,'zt 2 ',zt
       if(zt.ne.cmplx(0.,0.))zt=zt**st   !  jlm fix for DEC computers
c      print *,'zt**st ',zt
c      print *,'a2 ',a2
c      print *,'n1,w ',n1,w
       call tay(zt,a2,n1,w)
c      print *,'w,b,ci ',w,b,ci
       w=w+b
       w=(w+ci)/(ci-w)
c      print *,'n1,wa ',n1,w
      endif

      if(kx1.or.kxy1)w=w/cabs(w)**2
c      print *,'n1,wb ',n1,w
      x=real(w)
      y=aimag(w)
      if(kxy)then
       t=x
       x=y
       y=t
      endif
      if(ky)y=-y
      if(kx)x=-x
      w=cmplx(x,y)
c     print *,'n1,wc ',n1,w
      return
      end

c------------------------------------------------------------------------------c
c   r.j.purser, national meteorological center, washington d.c.  1995          c
c                   subroutine  ininmap                                         c
c  initialize the rotation matrix rot3 needed to transform standard            c
c  earth-centered cartesian components to the alternative cartesian frame      c
c  oriented so as to put geographical point (alat0,alon0) on the projection    c
c  axis.                                                                       c
c------------------------------------------------------------------------------c
      subroutine ininmap(alon0,alat0,rot3)
      dimension rot3(3,3)
      dimension x(3),t(3),u(3)
      print *,'entering ininmap, alat0 = ' ,alat0
      pi=4*atan(1.)
      dr=pi/180.
      blon0=dr*alon0
      blat0=dr*alat0
      clon0=cos(blon0)
      slon0=sin(blon0)
      clat0=cos(blat0)
      slat0=sin(blat0)
c     following as in J&J p.122 but for theta use -latitude (jlm)
      rot3(1,1)=clat0*clon0
      rot3(1,2)=clat0*slon0
      rot3(1,3)=slat0
      rot3(2,1)=-slon0
      rot3(2,2)=clon0
      rot3(2,3)=0.
      rot3(3,1)=-slat0*clon0
      rot3(3,2)=-slat0*slon0
      rot3(3,3)=clat0
      print *,'within ininmap, alat0, blat0, clat0, slat0 = '
     .                        ,alat0, blat0, clat0, slat0
      return

      entry nmap(rot3,x,t)
      call mulmm(rot3,x,u,3,3,1,3,3,3)
      do i=1,3
      t(i)=u(i)
      enddo
      return
      entry nmapt(rot3,x,t)
      call multm(rot3,x,u,3,3,1,3,3,3)
      do i=1,3
      t(i)=u(i)
      enddo
      return
      end

c------------------------------------------------------------------------------c
c   r.j.purser, national meteorological center, washington d.c.  1994          c
c                   subroutine  cinvrt                                          c
c  compute the taylor series coefficients z for the functional-inverse of      c
c  the function whose taylor series coefficients are w, for the case where     c
c  the constant coefficients of both z and w are 0                             c
c                                                                              c
c  --> w    taylor coefficients of original function (starting with linear)    c
c  <-- z    taylor coefficients of inverse function                            c
c  --> m    number of taylor series coefficients computed                      c
c  --- work workspace array consisting of at least 3*m elements                c
c------------------------------------------------------------------------------c
      subroutine cinvrt(w,z,m,work)
      dimension w(*),z(*),work(m,3)
      do i=1,m
       work(i,1)=w(i)
       work(i,2)=0.
      enddo
      work(1,2)=1.
      do j=1,m
       zj=work(j,2)/work(j,1)
       z(j)=zj
       do i=j,m
        work(i,2)=work(i,2)-zj*work(i,1)
       enddo
       call conv(work(1,1),w,work(1,3),m)
       do i=1,m
        work(i,1)=work(i,3)
       enddo
      enddo
      return
      end

c------------------------------------------------------------------------------c
c   r.j.purser, national meteorological center, washington d.c.  1994          c
c                   subroutine  conv                                            c
c   convolve double-precision series a with b to form c, up to m terms         c
c   starting with element 1                                                    c
c                                                                              c
c a,b   --> inputs (convolution factors)                                       c
c c     <-- output (convolution product)                                       c
c m     --> number of elements of a, b, c                                      c
c------------------------------------------------------------------------------c
      subroutine conv(a,b,c,m)
      dimension a(*),b(*),c(*)
      do k=1,m
       c(k)=0.
      enddo
      do i=1,m
       do j=1,m-i
        k=i+j
        c(k)=c(k)+a(i)*b(j)
       enddo
      enddo
      return
      end
