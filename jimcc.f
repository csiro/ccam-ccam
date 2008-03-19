      subroutine jimcc(em4,ax4,ay4,az4,xx4,yy4,il,myid)
c     like jim6.f but without stretch option
c     hedra1 data is hardwired
c     xx-->xx4, yy-->yy4, fm-->em4, dxa-->ax4, dxb-->ay4, dxc-->az4
c     include 'newmpar_gx.h'
c     include 'parm.h'    ! ktau
c      parameter(n=4*il ,np=n+1,non2=n/2)    !jlm for quad res. grid
      parameter(ipanel=2,ngrmax=1,ndiagj=0)
c      include 'bigxy4.h' ! common/bigxy4/xx4(iquad,iquad),yy4(iquad,iquad)
      integer :: myid,num
      real em4(4*il+1,4*il+1),
     .     ax4(4*il+1,4*il+1),ay4(4*il+1,4*il+1),az4(4*il+1,4*il+1)
      real xa(4*il+1,4*il+1),xb(4*il+1,4*il+1),xc(4*il+1,4*il+1)
      real*8 xx4(4*il+1,4*il+1),yy4(4*il+1,4*il+1)
c     equivalence (xb,xx4),(xc,yy4)  ! just to save on storage  jlm
      save num
      data num/0/
c     real dya(np,np),dyb(np,np),dyc(np,np)
c     ngr = 1  at unstaggered positions
c         = 2  at i+.5,j      positions
c         = 3  at i,j+.5      positions
c         = 4  at i+.5,j+.5   positions
      np=4*il+1
      CALL INROT

      do ngr=1,ngrmax
       call rgrid(xa,xb,xc,ax4,ay4,az4,            em4,np,ipanel,ngr)
c      these are values on the sphere

      if(num.eq.0.and.myid==0)then
       do j=1,np,np-1
        do i=1,np,np-1
         print *,'in jimcc xa,xb,xc: ',i,j,xa(i,j),xb(i,j),xc(i,j)
        enddo
       enddo
       do j=np/4,np+1-np/4,np+1-np/4-np/4
        do i=np/4,np+1-np/4,np+1-np/4-np/4
         print *,'in jimcc xa,xb,xc: ',i,j,xa(i,j),xb(i,j),xc(i,j)
        enddo
       enddo
       print *,'in jimcc xa,xb,xc: ',1,5,xa(1,5),(1,5),xc(1,5)
       print *,'now imposing 8-fold symmetry (jlm)'
      endif  ! (num.eq.0)

       do i=1,(np+1)/2
         xc(i,1)=max(-1.,xc(i,1)) ! for rounding errors
       enddo
       do j=1,(np+1)/2
        do i=j+1,(np+1/2)  ! rest of bottom LH corner
         xa(j,i)=xa(i,j)
         xb(j,i)=xc(i,j)
         xc(j,i)=xb(i,j)
        enddo
        do i=1,np/2        ! all of bottom RH corner
         xa(np+1-i,j)=xa(i,j)
         xb(np+1-i,j)=-xb(i,j)
         xc(np+1-i,j)=xc(i,j)
        enddo
        do i=1,np          ! all of top half
         xa(i,np+1-j)=xa(i,j)
         xb(i,np+1-j)=xb(i,j)
         xc(i,np+1-j)=-xc(i,j)
        enddo
       enddo

      if(num.eq.0.and.myid==0)then
       num=1
       do j=1,np,np-1
        do i=1,np,np-1
         print *,'in jimcc xa,xb,xc: ',i,j,xa(i,j),xb(i,j),xc(i,j)
        enddo
       enddo
       do j=np/4,np+1-np/4,np+1-np/4-np/4
        do i=np/4,np+1-np/4,np+1-np/4-np/4
         print *,'in jimcc xa,xb,xc: ',i,j,xa(i,j),xb(i,j),xc(i,j)
        enddo
       enddo
       print *,'in jimcc xa,xb,xc: ',1,5,xa(1,5),xb(1,5),xc(1,5)
      endif  ! (num.eq.0)

c      now convert these to just x,y values on the cube
       do j=1,np
        do i=1,np
         xx4(i,j)=real(xb(i,j),8)/real(xa(i,j),8)
         yy4(i,j)=real(xc(i,j),8)/real(xa(i,j),8)
        enddo
       enddo

c      if(ndiag.gt.0)then
c        now print the converted x,y values on the cube
c        print *,'xx4'
c        do j=np,1,-1
c         print *,(xx4(i,j),i=1,np)
c        enddo
c        print *,'yy4'
c        do j=np,1,-1
c         print *,(yy4(i,j),i=1,np)
c        enddo
c      endif

c      nw is dimension to be written out
c      nw=np
c      if(ngr.gt.1)nw=np-1   ! staggered cases
c      write (lu,*) nw,stretch,ngr
c      write (lu,*) ((xx4(i,j),i=1,nw),j=1,nw)
c      write (lu,*) ((yy4(i,j),i=1,nw),j=1,nw)
c      write (lu,*) n,' map factors for top face'
c      write (lu,*) ((em4(i,j),i=1,n),j=1,n)
c      if(ngr.eq.1.or.ngr.eq.3)then
c        write (lu,*) n,' ax4, ay4, az4 vector components for top panel'
c        write (lu,*) ((ax4(i,j),i=1,n),j=1,n)
c        write (lu,*) ((ay4(i,j),i=1,n),j=1,n)
c        write (lu,*) ((az4(i,j),i=1,n),j=1,n)
c      elseif(ngr.eq.4)then   ! for staggered v grid:
c        write (lu,*) n,' dya, dyb, dyc vector components for top panel'
c        write (lu,*) ((dya(i,j),i=1,n),j=1,n)
c        write (lu,*) ((dyb(i,j),i=1,n),j=1,n)
c        write (lu,*) ((dyc(i,j),i=1,n),j=1,n)
c      endif
c      if(ngr.eq.1)then
c        following only printed out for interest
c        i=1
c        j=(n+2)/2
c        edge=sqrt((xa(i+1,j)-xa(i,j))**2+(xb(i+1,j)-xb(i,j))**2
c    .         +(xc(i+1,j)-xc(i,j))**2)
c        i=(n+2)/2
c        j=(n+2)/2
c        cent=sqrt((xa(i+1,j)-xa(i,j))**2+(xb(i+1,j)-xb(i,j))**2
c    .         +(xc(i+1,j)-xc(i,j))**2)
c        print *,'edge,cent,rat ',edge,cent,edge/cent
c      endif  ! ngr.eq.1
      enddo  ! ngr loop
      return
      END

C------------------------------------------------------------------------------C
C   R.J.Purser, National Meteorological Center, Washington D.C.  1994          C
C                   SUBROUTINE  INROT                                          C
C                                                                              C
C   Initialize the rotation matrix needed to set the orientation               C
C   of the cube to user's specification, then call general initialization      C
C   routine INHEDRA to make the transformation code ready for use.             C
C                                                                              C
C  ROT is the user-defined orthogonal matrix determining the orientation of    C
C      the mapping-cube with respect to the standard earth-centered            C
C      coordinate basis defined in the comments below.                         C
C  IG1  is an array of 6 group indices. For designated map-panel IPANEL, the   C
C           element IG1(IPANEL) is associated with the octant of this map      C
C           in contact with the map-origin and that abuts the x-axis           C
C   Note: it is up to the user to ensure that ROT is indeed orthogonal and     C
C   that the 6 elements specified in IG1 do indeed belong to the 6 distinct    C
C   faces of the cube in accordance with the numbering convention adopted      C
C   for the triangulation into "group elements".                               C
C------------------------------------------------------------------------------C
      SUBROUTINE INROT
      DIMENSION ROT(3,3),IG1(6)
      DATA IG1/22,13,40,43,9,45/
      R2=SQRT(2.)
      R3=SQRT(3.)
      R6=R2*R3

C  SPECIFY THE ORIENTATION OF THE TRIANGULATED MAPPING-CUBE IN TERMS OF 3
C  ORTHOGONAL UNIT-VECTORS, REFERRED TO HERE AS p, q, r, DEFINED AS FOLLOWS:
C  VECTOR q POINTS TOWARDS THE MID-POINT OF EDGE-3 WHERE TRIANGULAR ELEMENTS
C  24, 29, 41 AND 36 MEET; VECTOR r POINTS TOWARDS VERTEX-0 WHERE ELEMENTS
C  0, 1, 2, 3, 4 AND 5 MEET; VECTOR p POINTS TOWARDS A POSITION ON THE
C  COMMON BOUNDARY OF ELEMENTS 14 AND 15 (IN FACE-3) SUCH THAT IT IS
C  ORTHOGONAL TO BOTH q AND r, OR EQUIVALENTLY, p IS DEFINED AS THE
C  (RIGHT-HANDED) CROSS-PRODUCT, q*r. THE BASIS-VECTORS USED TO EXPRESS
C  p, q, AND r ARE (1): THE UNIT VECTOR POINTING TO LAT,LONG=0,0; (2): THE
C  VECTOR POINTING TO LAT,LONG=0,90E; (3) THE UNIT VECTOR POINTING TO THE
C  NORTH POLE.

C  DEFINE VECTOR p AND MAKE IT THE FIRST COLUMN OF MATRIX ROT:
      ROT(1,1)=R6/6.
      ROT(2,1)=-R6/6.
      ROT(3,1)=-R6/3.

C  DEFINE VECTOR q AND MAKE IT THE SECOND COLUMN OF MATRIX ROT:
      ROT(1,2)=R2/2.
      ROT(2,2)=R2/2.
      ROT(3,2)=0.

C  DEFINE VECTOR r AND MAKE IT THE THIRD COLUMN OF MATRIX ROT:
      ROT(1,3)=R3/3.
      ROT(2,3)=-R3/3.
      ROT(3,3)=R3/3.

C  CUSTOMIZATION OF THE MAPPING TRANSFORMATION IS COMPLETED BY SPECIFYING,
C  FOR EACH NUMBERED MAP-PANEL (FROM 1 TO 6) THE TRIANGULAR ELEMENT THAT
C  CORRESPONDS TO THE 3 RESTRICTIONS IN THE LOCAL COORDINATES x,y:
C          a: x < .5;
C          b: 0. < y;
C          c: y < x.
C  FOR EACH MAP-PANEL, IP, THIS BASIC ELEMENT IS PRESCRIBED IN IG1(IP).
C  THESE, TOGETHER WITH THE ORTHOGONAL MATRIX ROT MADE UP OF p,q,r, ARE
C  PASSED TO THE GENERAL INITIALIZATION ROUTINE INHEDRA, AFTER WHICH, THE
C  MAP-TRANSFORMATION ROUTINES ARE READY FOR USE.
      CALL INHEDRA(ROT,IG1)
      RETURN
      END
c------------------------------------------------------------------------------c
c   r.j.purser, national meteorological center, washington d.c.  1994          c
c                   subroutine  rgrid                                          c
c                                                                              c
c   set up array of standard earth-centered coordinates of a chosen panel      c
c   of the rancic map. convention used for map coordinates here is that        c
c   each origin is the pole (north or south as appropriate) and the (x,y)      c
c   map coordinates in each panel form a right-handed pair, with x and y both  c
c   between 0 and 1. to choose another panel-corner for origin, or to alter    c
c   chiral convention, just rearrange the table igofkg. for more radical       c
c   change in map-coordinate convention, rewrite this routine!                 c
c                                                                              c
c  <-- xe,ye,ze     earth-centered coordinate of regular map-grid              c
c  --> np           number of points along each grid line (edges included)     c
c  --> ipanel       map-panel index [0,5]                                      c
c------------------------------------------------------------------------------c
c     subroutine rgrid(xe,ye,ze,dxa,dxb,dxc,dya,dyb,dyc,em4,
      subroutine rgrid(xe,ye,ze,dxa,dxb,dxc,            em4,
     . np,ipanel,ngr)
c     ngr = 1  at unstaggered positions
c         = 2  at i+.5,j+.5   positions
c         = 3  at i+.5,j      positions
c         = 4  at i,j+.5      positions
c     for staggered grids, actually only use 1,n part of the array
c set up earth-centered coordinates for points of chosen panel of the
c  rancic map
      parameter (stretch=1.,stretchm=1.-stretch)
      dimension xe(np,np),ye(np,np),ze(np,np)
      real dxa(np,np),dxb(np,np),dxc(np,np),em4(np,np)
c     real dya(np,np),dyb(np,np),dyc(np,np)
      dimension ex(3),xc(3,3),DFDX(2)
      n=np-1
      d=1./n
      xadd=0.
      yadd=0.
      if(ngr.eq.2.or.ngr.eq.3)xadd=.5*d
      if(ngr.eq.2.or.ngr.eq.4)yadd=.5*d
      do j=0,n
       jp=j+1
       y=j*d  + yadd   ! jlm allows staggered v
       y=.5+(y-.5)*(stretch+stretchm*(2.*y-1.)**2)   !jlm
       do i=0,n
        ip=i+1
        x=i*d + xadd   ! jlm allows staggered u
        x=.5+(x-.5)*(stretch+stretchm*(2.*x-1.)**2)   !jlm
        call mtoc(x,y,ipanel,ex)
        xe(ip,jp)=ex(1)
        ye(ip,jp)=ex(2)
        ze(ip,jp)=ex(3)
        call mtocd(x,y,ipanel,xc,em4(ip,jp),dfdx)
c       return dxa etc as unit vectors
        den=sqrt(xc(1,1)**2+xc(2,1)**2+xc(3,1)**2)
        if(den.lt.1.e-6)den=1.
        dxa(ip,jp)=xc(1,1)/den   ! the three components of a vector along dx
        dxb(ip,jp)=xc(2,1)/den
        dxc(ip,jp)=xc(3,1)/den
c       den=sqrt(xc(1,2)**2+xc(2,2)**2+xc(3,2)**2)
c       if(den.lt.1.e-6)den=1.
c       dya(ip,jp)=xc(1,2)/den   ! the three components of a vector along dy
c       dyb(ip,jp)=xc(2,2)/den
c       dyc(ip,jp)=xc(3,2)/den
       enddo
      enddo
      return
      end

C                                               *****************
C                                               *   HED6A.FOR   *
C                                               *  PURSER 1994  *
C                                               *****************
C
C  SUBROUTINES NEEDED TO PERFORM RANCIC-CUBE TRANSFORMATIONS
C
C------------------------------------------------------------------------------C
C   R.J.Purser, National Meteorological Center, Washington D.C.  1994          C
C               BLOCKDATA                                                      C
C   Initialize some of the arrays in common, particularly those associated     C
C   with the implicit representation of the group of symmetries of the cube    C
C------------------------------------------------------------------------------C
      BLOCKDATA
      complex ci,cip4,cip3oss
      common/cstoc/
     *a(30),b(30),ss,third,fourth,ci,cip4,cip3oss
      COMMON/CHEDRA/
     * KDA(0:47),KDB(0:47),KDC(0:47)
     *,KNA(0:47),KNB(0:47),KNC(0:47)
     *,IGOFKG(8,6),IPOFIG(0:47),KGOFIG(0:47)
     *,IG,FLIP8(2,2,8),F(3,0:5),E(3,0:11),TXE(3,3)
     *,ROTG(3,3,0:47)
      DATA KDA/0,1,1,2,2,0, 0,1,1,5,5,0, 3,1,1,2,2,3, 3,1,1,5,5,3
     *,        0,4,4,2,2,0, 0,4,4,5,5,0, 3,4,4,2,2,3, 3,4,4,5,5,3/
     *,KDB/3,10,4,11,5,9, 6,7,1,11,5,0, 3,1,7,8,2,9,  6,4,10,8,2,0
     *,     0,10,4,2,8,6, 9,7,1,2,8,3,  0,1,7,5,11,6, 9,4,10,5,11,3/
     *,KDC/5,11,3,9,4,10, 5,11,6,0,1,7, 2,8,3,9,7,1,  2,8,6,0,10,4
     *,     8,2,0,6,4,10, 8,2,9,3,1,7,  11,5,0,6,7,1, 11,5,9,3,10,4/
     *,KNA/12,25,26,9,10,17,   18,31,32,3,4,23
     *,     0,37,38,21,22,5,   6,43,44,15,16,11
     *,     36,1,2,33,34,41,   42,7,8,27,28,47
     *,     24,13,14,45,46,29, 30,19,20,39,40,35/
     *,KNB/ 5,2,1,4,3,0,       11,8,7,10,9,6
     *,     17,14,13,16,15,12, 23,20,19,22,21,18
     *,     29,26,25,28,27,24, 35,32,31,34,33,30
     *,     41,38,37,40,39,36, 47,44,43,46,45,42/
     *,KNC/ 1,0,3,2,5,4,       7,6,9,8,11,10
     *,     13,12,15,14,17,16, 19,18,21,20,23,22
     *,     25,24,27,26,29,28, 31,30,33,32,35,34
     *,     37,36,39,38,41,40, 43,42,45,44,47,46/
     *,IG/0/
      DATA FLIP8
     */1.,0.,0.,1.,  -1.,0.,0.,1.,  1.,0.,0.,-1.,  -1.,0.,0.,-1.
     *,0.,1.,1.,0.,  0.,-1.,1.,0.,  0.,1.,-1.,0.,  0.,-1.,-1.,0./

c     following originally read from hedra1.dat
      data a/1.47713062600964,-.38183510510174,-.05573058001191
     . ,-.00895883606818,-.00791315785221,-.00486625437708
     . ,-.00329251751279,-.00235481488325,-.00175870527475
     . ,-.00135681133278,-.00107459847699,-.00086944475948
     . ,-.00071607115121,-.00059867100093,-.00050699063239
     . ,-.00043415191279,-.00037541003286,-.00032741060100
     . ,-.00028773091482,-.00025458777519,-.00022664642371
     . ,-.00020289261022,-.00018254510830,-.00016499474460
     . ,-.00014976117167,-.00013646173947,-.00012478875822
     . ,-.00011449267279,-.00010536946150,-.00009725109376/
      data b/.67698819751739,.11847293456554,.05317178134668
     . ,.02965810434052,.01912447304028,.01342565621117,.00998873323180
     . ,.00774868996406,.00620346979888,.00509010874883,.00425981184328
     . ,.00362308956077,.00312341468940,.00272360948942,.00239838086555
     . ,.00213001905118,.00190581316131,.00171644156404,.00155493768255
     . ,.00141600715207,.00129556597754,.00119042140226,.00109804711790
     . ,.00101642216628,.00094391366522,.00087919021224,.00082115710311
     . ,.00076890728775,.00072168382969,.00067885087750/
      END

C------------------------------------------------------------------------------C
C   R.J.Purser, National Meteorological Center, Washington D.C.  1994          C
C                   SUBROUTINE INHEDRA                                         C
C   Initialize variables needed to perform the Rancic-transformation and       C
C   its inverse                                                                C
C------------------------------------------------------------------------------C
      SUBROUTINE INHEDRA(ROT,IG1)
      COMPLEX CI,CIP4,CIP3OSS
      COMMON/CSTOC/
     *A(30),B(30),SS,THIRD,FOURTH,CI,CIP4,CIP3OSS
      COMMON/CHEDRA/
     * KDA(0:47),KDB(0:47),KDC(0:47)
     *,KNA(0:47),KNB(0:47),KNC(0:47)
     *,IGOFKG(8,6),IPOFIG(0:47),KGOFIG(0:47)
     *,IG,FLIP8(2,2,8),F(3,0:5),E(3,0:11),TXE(3,3)
     *,ROTG(3,3,0:47)
      DIMENSION ROT(3,3),IG1(6),IGK(8)
C  SET UP CROSS-REFERENCE TABLES CONNECTING GROUP-ELEMENTS IG TO
C  USER-DEFINED NUMBERING AND ORIENTATIONS OF PANELS:
      DO IP=1,6
       IGK(1)=IG1(IP)
       IGK(2)=KNA(IGK(1))
       IGK(5)=KNC(IGK(1))
       IGK(6)=KNC(IGK(2))
       IGK(7)=KNA(IGK(5))
       IGK(8)=KNA(IGK(6))
       IGK(3)=KNC(IGK(7))
       IGK(4)=KNC(IGK(8))
       DO KG=1,8
        IG=IGK(KG)
        IGOFKG(KG,IP)=IG
        IPOFIG(IG)=IP
        KGOFIG(IG)=KG
       ENDDO
      ENDDO
      R2=SQRT(2.)
      R3=SQRT(3.)
      R6=SQRT(6.)
      R2O2=R2/2.
      R3O2=R3/2.
      R3O3=R3/3.
      R6O3=R6/3.
      R6O6=R6/6.
      R3O6=R3/6.
      SS=R2
      THIRD=1./3
      FOURTH=1./4
      CI=CMPLX(0.,1.)
      CIP4=CI**FOURTH
      CIP3OSS=CI**THIRD/SS
      F(1,0)=-R6O3
      F(2,0)=0.
      F(3,0)=R3O3
      F(1,1)=R6O6
      F(2,1)=-R2O2
      F(3,1)=R3O3
      F(1,2)=R6O6
      F(2,2)=R2O2
      F(3,2)=R3O3
      E(1,0)=R3O3
      E(2,0)=0.
      E(3,0)=R6O3
      E(1,1)=-R3O6
      E(2,1)=.5
      E(3,1)=R6O3
      E(1,2)=-R3O6
      E(2,2)=-.5
      E(3,2)=R6O3
      E(1,3)=0.
      E(2,3)=1.
      E(3,3)=0.
      E(1,4)=-R3O2
      E(2,4)=-.5
      E(3,4)=0.
      E(1,5)=R3O2
      E(2,5)=-.5
      E(3,5)=0.
      DO J=0,2
       K=J+3
       DO I=1,3
        F(I,K)=-F(I,J)
       ENDDO
      ENDDO
      DO J=0,5
       K=J+6
       DO I=1,3
        E(I,K)=-E(I,J)
       ENDDO
      ENDDO
      DO I=1,3
       TXE(1,I)=F(I,0)
       TXE(2,I)=E(I,3)
       TXE(3,I)=E(I,5)
      ENDDO
      CALL INVMM(TXE,TXE,3,3,3)
C  ROTATE THE 6 FACE-VECTORS, F, TO USER-DEFINED ORIENTATION:
      CALL MULMM(ROT,F,ROTG,3,3,6,3,3,9)
      CALL COPM(ROTG,F,3,6,9,3)

C  ROTATE THE 12 EDGE-VECTORS, E, TO USER-DEFINED ORIENTATION:
      CALL MULMM(ROT,E,ROTG,3,3,12,3,3,9)
      CALL COPM(ROTG,E,3,12,9,3)

C  BASED ON THE PRESCRIBED ORIENTATION (DEFINED BY "ROT"),
C  CONSTRUCT THE ROTATION MATRIX ASSOCIATED WITH EACH GROUP ELEMENT LG:

      DO LG=0,47
       DO J=1,3
        DO I=1,3
         ROTG(I,J,LG)=F(I,KDA(LG))*TXE(J,1)
     *               +E(I,KDB(LG))*TXE(J,2)
     *               +E(I,KDC(LG))*TXE(J,3)
        ENDDO
       ENDDO
      ENDDO

c     OPEN(UNIT=1,FILE='hedra1.dat'
c    *,ACCESS='SEQUENTIAL',STATUS='UNKNOWN')
c     DO I=1,30
c     READ(1,100)A(I),B(I)
c100  FORMAT(2(F20.14,1X))
c     ENDDO
c     CLOSE(1)
      RETURN
      END

C------------------------------------------------------------------------------C
C   R.J.Purser, National Meteorological Center, Washington D.C.  1994          C
C                   SUBROUTINE MTOC                                            C
C   Transform from map-coordinates to standard earth-centered cartesians       C
C                                                                              C
C  -->    XX,YY  map-coordinates                                               C
C  -->    IPANEL map-panel index                                               C
C  <--    XC     standard earth-centered cartesian coordinates                 C
C------------------------------------------------------------------------------C
      SUBROUTINE MTOC(XX,YY,IPANEL,XC)
      COMPLEX CI,W,Z,CIP4,CIP3OSS,ARG
      COMMON/CSTOC/
     *A(30),B(30),SS,THIRD,FOURTH,CI,CIP4,CIP3OSS
      COMMON/CHEDRA/
     * KDA(0:47),KDB(0:47),KDC(0:47)
     *,KNA(0:47),KNB(0:47),KNC(0:47)
     *,IGOFKG(8,6),IPOFIG(0:47),KGOFIG(0:47)
     *,IG,FLIP8(2,2,8),F(3,0:5),E(3,0:11),TXE(3,3)
     *,ROTG(3,3,0:47)
      DIMENSION XC(3),XV(3)
      X=XX
      Y=YY
      KG=1
      IF(X.GT..5)THEN
       KG=KG+1
       X=1.-X
      ENDIF
      IF(Y.GT..5)THEN
       KG=KG+2
       Y=1.-Y
      ENDIF
      IF(Y.GT.X)THEN
       KG=KG+4
       T=X
       X=Y
       Y=T
      ENDIF
c     Z=CMPLX(X,Y)**4
      z=cmplx(x,y)*cmplx(x,y)
      z=z*z
      CALL TAY(Z,A,30,W)
      ARG = -CI*W                                      ! mrd
      IF ( CABS(ARG).EQ.0. ) THEN                      ! mrd
         W = (0.,0.)                                   ! mrd
      ELSE                                             ! mrd
         W=CIP3OSS*(-CI*W)**THIRD                      ! mrd
      END IF                                           ! mrd
cjim  W=CIP3OSS*(-CI*W)**THIRD
cjlm  if(w.ne.cmplx(0.,0.))w=cip3oss*(-ci*w)**third    ! jlm
      XW=REAL(W)
      YW=AIMAG(W)
      H=2./(1.+XW*XW+YW*YW)
      XV(1)=XW*H
      XV(2)=YW*H
      XV(3)=H-1.
      IG=IGOFKG(KG,IPANEL)
      CALL MULMM(ROTG(1,1,IG),XV,XC,3,3,1,3,3,3)
      RETURN
      END

C------------------------------------------------------------------------------C
C   R.J.Purser, National Meteorological Center, Washington D.C.  1994          C
C                   SUBROUTINE MTOCD                                           C
C   Transform from map-coordinates to standard earth-centered cartesians       C
C   and simultaneously accumulate the derivative of the transformation in      C
C   order to provide information on map-scaling factor and relative            C
C   orientation                                                                C
C                                                                              C
C  --> XX,YY  map-coordinates (from corner IG, each panel a unit square)       C
C  --> IPANEL map-panel index                                                  C
C  <-- XC   augmented jacobian matrix: first two columns represent the         C
C           derivative with respect to X and Y of the earth-centered           C
C           cartesian coordinates of the point corresponding to map-image      C
C           X,Y. These cartesian coordinates themselves are inserted into      C
C           column-3 of XDC.                                                   C
C  <-- em4   map-factor at this point                                           C
C  <-- DFDX x- and y-derivatives of map-factor here                            C
C------------------------------------------------------------------------------C
      SUBROUTINE MTOCD(XX,YY,IPANEL,XC,em4,DFDX)
      COMPLEX CI,W,Z,CIP4,CIP3OSS,ZU,WU,CD,CDD,ARG
      COMMON/CSTOC/
     *A(30),B(30),SS,THIRD,FOURTH,CI,CIP4,CIP3OSS
      COMMON/CHEDRA/
     * KDA(0:47),KDB(0:47),KDC(0:47)
     *,KNA(0:47),KNB(0:47),KNC(0:47)
     *,IGOFKG(8,6),IPOFIG(0:47),KGOFIG(0:47)
     *,IG,FLIP8(2,2,8),F(3,0:5),E(3,0:11),TXE(3,3)
     *,ROTG(3,3,0:47)
      DIMENSION XDC(3,3),XD(3,2),V1(2),DFDX(2)
      DIMENSION XC(3,3)
      X=XX
      Y=YY
      KG=1
      IF(X.GT..5)THEN
       KG=KG+1
       X=1.-X
      ENDIF
      IF(Y.GT..5)THEN
       KG=KG+2
       Y=1.-Y
      ENDIF
      IF(Y.GT.X)THEN
       KG=KG+4
       T=X
       X=Y
       Y=T
      ENDIF
      ZU=CMPLX(X,Y)
      Z=ZU**4
      CALL TAYDD(Z,A,30,W,CD,CDD)
cjim  WU=CIP3OSS*(-CI*W)**THIRD
cjlm  wu=cmplx(0.,0.)                                  ! jlm
cjlm  if(w.ne.cmplx(0.,0.))wu=cip3oss*(-ci*w)**third   ! jlm
      ARG = -CI*W                                      ! mrd
      IF ( CABS(ARG).EQ.0. ) THEN                      ! mrd
         WU = (0.,0.)                                  ! mrd
      ELSE                                             ! mrd
         WU=CIP3OSS*(-CI*W)**THIRD                     ! mrd
      END IF                                           ! mrd
      XW=REAL(WU)
      YW=AIMAG(WU)
      XWXW=XW*XW
      XWYW=XW*YW
      YWYW=YW*YW
      H=2./(1.+XWXW+YWYW)
      HH=H*H
      XDC(1,3)=XW*H
      XDC(2,3)=YW*H
      XDC(3,3)=H-1.
      XDC(1,1)=H-HH*XWXW
      XDC(2,1)= -HH*XWYW
      XDC(3,1)= -HH*XW
      XDC(1,2)=XDC(2,1)
      XDC(2,2)=H-HH*YWYW
      XDC(3,2)= -HH*YW
      IF(CABS(Z).EQ.0.)THEN
       CD=0.
       CDD=0.
       em4=0.
       V1(1)=0.
       V1(2)=0.
      ELSE
       CD=4.*WU*CD*Z/(3.*W*ZU)
       CDD=3.*CD/ZU-2.*CD*CD/WU+16.*WU*Z**2*CDD/(3.*W*ZU**2)
       RD=REAL(CD)
       QD=AIMAG(CD)
       RDD=REAL(CDD)
       QDD=AIMAG(CDD)
       S=SQRT(RD*RD+QD*QD)
       DSDX=(RDD*RD+QDD*QD)/S
       DSDY=(RDD*QD-QDD*RD)/S
       DHDX=-HH*(XW*RD+YW*QD)
       DHDY=-HH*(-XW*QD+YW*RD)
       em4=H*S
       V1(1)=DHDX*S+H*DSDX
       V1(2)=DHDY*S+H*DSDY
      ENDIF
      RD=REAL(CD)
      QD=AIMAG(CD)
      DO I=1,3
       XD(I,1)= XDC(I,1)*RD+XDC(I,2)*QD
       XD(I,2)=-XDC(I,1)*QD+XDC(I,2)*RD
      ENDDO
      IG=IGOFKG(KG,IPANEL)
      CALL MULMM(V1,FLIP8(1,1,KG),DFDX,1,2,2,1,2,1)
      CALL MULMM(XD,FLIP8(1,1,KG),XDC,3,2,2,3,2,3)
      CALL MULMM(ROTG(1,1,IG),XDC,XC,3,3,3,3,3,3)
      RETURN
      END


C------------------------------------------------------------------------------C
C   R.J.Purser, National Meteorological Center, Washington D.C.  1994          C
C                   SUBROUTINE FLIP                                            C
C   Use standard earth-centered cartesian coordinates of a point to determine  C
C   the group-element IG of the orthogonal transformation needed to transform  C
C   this point into the standard triangular wedge, and the new cartesian       C
C   coordinates of the point once it has indergone this transformation. The    C
C   transformed position puts the point close to the pole and with small non-  C
C   negative components XC(1),XC(2), to ensure that it lies well inside the    C
C   circle of convergence of the Taylor series and with appropriate azimuth    C
C   to ensure intended solution is obtained when fractional power is taken.    C
C  --> XE standard earth-centered cartesian coordinates [3]                    C
C  <-- XC transformed earth-centered cartesian coordinates [3]                 C
C------------------------------------------------------------------------------C
      SUBROUTINE FLIP(XE,XC)
      COMMON/CHEDRA/
     * KDA(0:47),KDB(0:47),KDC(0:47)
     *,KNA(0:47),KNB(0:47),KNC(0:47)
     *,IGOFKG(8,6),IPOFIG(0:47),KGOFIG(0:47)
     *,IG,FLIP8(2,2,8),F(3,0:5),E(3,0:11),TXE(3,3)
     *,ROTG(3,3,0:47)
      DIMENSION XE(3),XC(3)
C  VALIDITY OF CONDITIONS A & B & C UNCERTAIN:
      DAX=DOT(F(1,KDA(IG)),XE,3)
      DBX=DOT(E(1,KDB(IG)),XE,3)
      IF(DAX.LT.0.)THEN
       DAX=-DAX
       IG=KNA(IG)
      ENDIF
      IF(DBX.LT.0.)THEN
       DBX=-DBX
       IG=KNB(IG)
      ENDIF

C  VALIDITY OF CONDITIONS A & B TRUE, C UNCERTAIN:
400   DCX=DOT(E(1,KDC(IG)),XE,3)
      IF(DCX.LT.0.)THEN
       DCX=-DCX
       IG=KNC(IG)

C  VALIDITY OF CONDITIONS A & B UNCERTAIN, C TRUE:
       DAX=DOT(F(1,KDA(IG)),XE,3)
       DBX=DOT(E(1,KDB(IG)),XE,3)
       IF(DAX.LT.0.)THEN
        DAX=-DAX
        IG=KNA(IG)
        IF(DBX.LT.0.)THEN
         DBX=-DBX
         IG=KNB(IG)
        ENDIF
       ELSE
        IF(DBX.GE.0.)GOTO 300
        DBX=-DBX
        IG=KNB(IG)
       ENDIF
       GOTO 400
      ENDIF

C  VALIDITY OF CONDITIONS A & B & C TRUE:
300   DO I=1,3
        XC(I)=TXE(I,1)*DAX+TXE(I,2)*DBX+TXE(I,3)*DCX
      ENDDO
      RETURN
      END

C------------------------------------------------------------------------------C
C   R.J.Purser, National Meteorological Center, Washington D.C.  1994          C
C                   SUBROUTINE STOC                                            C
C   Transform from latitude,longitude, to corner-origin coordinates of one     C
C   of the maps.                                                               C
C   -->   DLAT, DLON    Latitude and longitude (degrees)                       C
C   <--   X,Y           map coordinates scaled to the unit square              C
C   <--   IPANEL        map-panel index                                        C
C------------------------------------------------------------------------------C
      SUBROUTINE STOC(DLAT,DLON,X,Y,IPANEL)
      DIMENSION XE(3)
      CALL STOE(DLAT,DLON,XE)
      CALL CTOM(XE,X,Y,IPANEL)
      RETURN
      END

C------------------------------------------------------------------------------C
C   R.J.Purser, National Meteorological Center, Washington D.C.  1994          C
C                   SUBROUTINE CTOM                                            C
C   Transform from group-oriented cartesians to map-coordinates.               C
C   Use CTTOM to transform first to corner-coordinates. Then use group         C
C   element IG to determine map-panel and perform the conversion from          C
C   corner-coordinates (0<y<x<.5) to full map-panel coordinates (0<x<1, 0<y<1) C
C  -->  XE      earth-centered cartesians                                      C
C  <--  X,Y     map-coordinates                                                C
C  <--  IPANEL  map-panel index                                                C
C------------------------------------------------------------------------------C
      SUBROUTINE CTOM(XE,X,Y,IPANEL)
      COMMON/CHEDRA/
     * KDA(0:47),KDB(0:47),KDC(0:47)
     *,KNA(0:47),KNB(0:47),KNC(0:47)
     *,IGOFKG(8,6),IPOFIG(0:47),KGOFIG(0:47)
     *,IG,FLIP8(2,2,8),F(3,0:5),E(3,0:11),TXE(3,3)
     *,ROTG(3,3,0:47)
      DIMENSION XE(3),XC(3)
      CALL FLIP(XE,XC)
      IPANEL=IPOFIG(IG)
      KG=KGOFIG(IG)
      CALL CTTOM(XC,X,Y)
      IF(KG.GT.4)THEN
       T=X
       X=Y
       Y=T
       KG=KG-4
      ENDIF
      IF(KG.GT.2)THEN
       Y=1.-Y
       KG=KG-2
      ENDIF
      IF(KG.GT.1)THEN
       X=1.-X
      ENDIF
      RETURN
      END

C------------------------------------------------------------------------------C
C   R.J.Purser, National Meteorological Center, Washington D.C.  1994          C
C                   SUBROUTINE CTTOM                                           C
C   Transform from earth-centered cartesians to map-coordinates                C
C   Point must lie in standard wedge near the pole. Transform first to         C
C   complex-position in the stereographic plane. Then apply fractional-        C
C   power of Taylor-series.                                                    C
C  -->  XE  earth-centered cartesians                                          C
C  <--  X,Y corner-origin map-coordinates oriented with respect to principal   C
C           edge.                                                              C
C------------------------------------------------------------------------------C
      SUBROUTINE CTTOM(XE,X,Y)
      COMPLEX CI,W,Z,CIP4,CIP3OSS
      COMMON/CSTOC/
     *A(30),B(30),SS,THIRD,FOURTH,CI,CIP4,CIP3OSS
      DIMENSION XE(3)
C  SS=SQRT(2.) AND REPRESENTS THE SCALE FACTOR NEEDED
C  TO PLACE THE NEIGHBORING 3 VERTICES OF THE NOMINATED POLE OF THE
C      6-HEDRON ON THE UNIT CIRCLE OF THE RESCALED STEREOGRAPHIC MAP
C  CIP4 IS THE COMPLEX 4th-ROOT OF i (i.e., 8th-ROOT OF -1)
      HI=SS/(1.+XE(3))
      W=(CMPLX(XE(1),XE(2))*HI)**3
      CALL TAY(W,B,30,Z)
C  ROTATE AWAY FROM THE BRANCH-CUT THAT MARKS THE NEGATIVE-REAL AXIS:
C  TAKE 4th-ROOT AND ROTATE BACK AGAIN
      Z=CIP4*(-CI*Z)**FOURTH
      X=REAL(Z)
      Y=AIMAG(Z)
      RETURN
      END

C------------------------------------------------------------------------------C
C   R.J.Purser, National Meteorological Center, Washington D.C.  1994          C
C                   SUBROUTINE STOE                                            C
C   Transform latitude and longitude (degrees) to earth-centered cartesian     C
C   coordinates.                                                               C
C  --> DLAT     latitude                                                       C
C  --> DLON     longitude                                                      C
C  <-- XE       three cartesian components.                                    C
C------------------------------------------------------------------------------C
      SUBROUTINE STOE(DLAT,DLON,XE)
      DIMENSION XE(3)
C  FOR 64-BIT PRECISION, USE:
C      DATA DTOR/.01745329251994329577/
C  FOR 32-BIT PRECISION, USE:
      DATA DTOR/.017453293/
      RLAT=DTOR*DLAT
      RLON=DTOR*DLON
      SLA=SIN(RLAT)
      CLA=COS(RLAT)
      SLO=SIN(RLON)
      CLO=COS(RLON)
      XE(1)=CLA*CLO
      XE(2)=CLA*SLO
      XE(3)=SLA
      RETURN
      END

C------------------------------------------------------------------------------C
C   R.J.Purser, National Meteorological Center, Washington D.C.  1994          C
C                   SUBROUTINE  TAY                                            C
C  Evaluate the complex function W of Z whose real                             C
C  Taylor series coefficients are RA.                                          C
C                                                                              C
C  --> Z    function argument (complex)                                        C
C  --> RA   Taylor coefficients (real)                                         C
C  --> N    number of coefficients (starting with the linear term)             C
C  <-- W    Taylor-series approximation of the function (complex)              C
C------------------------------------------------------------------------------C
      SUBROUTINE TAY(Z,RA,N,W)
      COMPLEX Z,W
      DIMENSION RA(*)
      W=0.
      DO I=N,1,-1
c      print *,'in tay i,w,ra(i),z ',i,w,ra(i),z
       W=(W+RA(I))*Z
      ENDDO
      RETURN
      END

C------------------------------------------------------------------------------C
C   R.J.Purser, National Meteorological Center, Washington D.C.  1994          C
C                   SUBROUTINE  TAYD                                           C
C  Evaluate the complex function W of Z whose real                             C
C  Taylor series coefficients are RA, together with its derivative WD.         C
C                                                                              C
C  --> Z    function argument (complex)                                        C
C  --> RA   Taylor coefficients (real)                                         C
C  --> N    number of coefficients (starting with the linear term)             C
C  <-- W    Taylor-series approximation of the function (complex)              C
C  <-- WD   Taylor series approximation of the derivative of the function W    C
C------------------------------------------------------------------------------C
      SUBROUTINE TAYD(Z,RA,N,W,WD)
      COMPLEX Z,W,WD
      DIMENSION RA(*)
      W=0.
      WD=0.
      DO I=N,1,-1
       W=(W+RA(I))*Z
       WD=Z*WD+I*RA(I)
      ENDDO
      RETURN
      END

C------------------------------------------------------------------------------C
C   R.J.Purser, National Meteorological Center, Washington D.C.  1994          C
C                   SUBROUTINE  TAYDD                                          C
C  Evaluate the complex function W of Z whose real                             C
C  Taylor series coefficients are RA, together with its derivative WD.         C
C  and second derivative WDD                                                   C
C                                                                              C
C  --> Z    function argument (complex)                                        C
C  --> RA   Taylor coefficients (real)                                         C
C  --> N    number of coefficients (starting with the linear term)             C
C  <-- W    Taylor-series approximation of the function (complex)              C
C  <-- WD   Taylor series approximation of the derivative of the function W    C
C  <-- WDD  Taylor series approximation of the derivative of the function WD   C
C------------------------------------------------------------------------------C
      SUBROUTINE TAYDD(Z,RA,N,W,WD,WDD)
      COMPLEX Z,W,WD,WDD
      DIMENSION RA(*)
      W=0.
      WD=0.
      WDD=0.
      DO I=N,1,-1
       W=(W+RA(I))*Z
       WD=Z*WD+I*RA(I)
      ENDDO
      DO I=N,2,-1
       WDD=Z*WDD+I*(I-1)*RA(I)
      ENDDO
      RETURN
      END

      FUNCTION DOT(A,B,M)
      DIMENSION A(M),B(M)
      DOT=0.
      DO I=1,M
       DOT=DOT+A(I)*B(I)
      ENDDO
      RETURN
      END

C------------------------------------------------------------------------------C
C   R.J.Purser, National Meteorological Center, Washington D.C.  1994          C
C                   SUBROUTINE COPM etc                                        C
C  Mainly routines to copy matrices                                            C
C------------------------------------------------------------------------------C
      SUBROUTINE COPM(A,B,MI,MJ,NA,NB)
      DIMENSION A(NA,*),B(NB,*)
      ENTRY EQMM(A,B,MI,MJ,NA,NB)
      DO 230 I=1,MI
      DO 230 J=1,MJ
  230 B(I,J)=A(I,J)
      RETURN
      ENTRY NEQMM(A,B,MI,MJ,NA,NB)
      ENTRY CONM(A,B,MI,MJ,NA,NB)
      DO I=1,MI
      DO J=1,MJ
         B(I,J)=-A(I,J)
      ENDDO
      ENDDO
      RETURN
      ENTRY MULMS(A,SS,B,MI,MJ,NA,NB)
      DO I=1,MI
      DO J=1,MJ
         B(I,J)=A(I,J)*SS
      ENDDO
      ENDDO
      RETURN
      ENTRY EQTM(A,B,MI,MJ,NA,NB)
      ENTRY COPT(A,B,MI,MJ,NA,NB)
      DO I=1,MI
      DO J=1,MJ
       B(I,J)=A(J,I)
      ENDDO
      ENDDO
      RETURN
      ENTRY CONT(A,B,MI,MJ,NA,NB)
      DO I=1,MI
      DO J=1,MJ
       B(I,J)=-A(J,I)
      ENDDO
      ENDDO
      RETURN
      ENTRY ZERM(A,MI,MJ,NA)
      DO I=1,MI
      DO J=1,MJ
       A(I,J)=0.
      ENDDO
      ENDDO
      RETURN
      END
      SUBROUTINE MULMM(A,B,C,MI,MJ,MK,NA,NB,NC)
      LOGICAL OMUL
      DIMENSION A(NA,*),B(NB,*),C(NC,*)
      OMUL=.TRUE.
      GOTO 320
      ENTRY MADMM(A,B,C,MI,MJ,MK,NA,NB,NC)
      OMUL=.FALSE.
320   DO 200 I=1,MI
      DO 200 K=1,MK
      IF(OMUL)C(I,K)=0.
      DO 200 J=1,MJ
  200 C(I,K)=C(I,K)+A(I,J)*B(J,K)
      RETURN
      ENTRY MULMT(A,B,C,MI,MJ,MK,NA,NB,NC)
      OMUL=.TRUE.
      GOTO 321
      ENTRY MADMT(A,B,C,MI,MJ,MK,NA,NB,NC)
      OMUL=.FALSE.
321   DO 201 I=1,MI
      DO 201 K=1,MK
      IF(OMUL)C(I,K)=0.
      DO 201 J=1,MJ
  201 C(I,K)=C(I,K)+A(I,J)*B(K,J)
      RETURN
      ENTRY MULTM(A,B,C,MI,MJ,MK,NA,NB,NC)
      OMUL=.TRUE.
      GOTO 322
      ENTRY MADTM(A,B,C,MI,MJ,MK,NA,NB,NC)
      OMUL=.FALSE.
322   DO 202 I=1,MI
      DO 202 K=1,MK
      IF(OMUL)C(I,K)=0.
      DO 202 J=1,MJ
  202 C(I,K)=C(I,K)+A(J,I)*B(J,K)
      RETURN
      ENTRY MULTT(A,B,C,MI,MJ,MK,NA,NB,NC)
      OMUL=.TRUE.
      GOTO 323
      ENTRY MADTT(A,B,C,MI,MJ,MK,NA,NB,NC)
      OMUL=.FALSE.
323   DO 203 I=1,MI
      DO 203 K=1,MK
      IF(OMUL)C(I,K)=0.
      DO 203 J=1,MJ
  203 C(I,K)=C(I,K)+A(J,I)*B(K,J)
      RETURN
      ENTRY MSBMM(A,B,C,MI,MJ,MK,NA,NB,NC)
      DO 206 I=1,MI
      DO 206 K=1,MK
      DO 206 J=1,MJ
206   C(I,K)=C(I,K)-A(I,J)*B(J,K)
      RETURN
      ENTRY MSBMT(A,B,C,MI,MJ,MK,NA,NB,NC)
      DO 207 I=1,MI
      DO 207 K=1,MK
      DO 207 J=1,MJ
207   C(I,K)=C(I,K)-A(I,J)*B(K,J)
      RETURN
      ENTRY MSBTM(A,B,C,MI,MJ,MK,NA,NB,NC)
      DO 208 I=1,MI
      DO 208 K=1,MK
      DO 208 J=1,MJ
208   C(I,K)=C(I,K)-A(J,I)*B(J,K)
      RETURN
      ENTRY MSBTT(A,B,C,MI,MJ,MK,NA,NB,NC)
      DO 209 I=1,MI
      DO 209 K=1,MK
      DO 209 J=1,MJ
209   C(I,K)=C(I,K)-A(J,I)*B(K,J)
      RETURN
      END

C------------------------------------------------------------------------------C
C   R.J.Purser, National Meteorological Center, Washington D.C.  1993          C
C                   SUBROUTINE LUFM                                            C
C  perform l-u decomposition of square matrix a in place with                  C
C  partial pivoting                                                            C
C  For DOUBLE PRECISION version see DLUFM                                      C
C                                                                              C
C  --> a    square matrix to be factorized                                     C
C  <-- ipiv array encoding the pivoting sequence                               C
C  <-- d    indicator for possible sign change of determinant                  C
C  --> m    degree of (active part of) a                                       C
C  --> na   first fortran dimension of a                                       C
C                                                                              C
C------------------------------------------------------------------------------C
      subroutine LUFM(a,ipiv,d,m,na)
      dimension a(na,*),ipiv(*)
      d=1.
      ipiv(m)=m
      do j=1,m-1
       jp=j+1
       abig=abs(a(j,j))
       iquad=j
       do i=jp,m
        aa=abs(a(i,j))
        if(aa.gt.abig)then
         iquad=i
         abig=aa
        endif
       enddo
c  swap rows, recording changed sign of determinant
       ipiv(j)=iquad
       if(iquad.ne.j)then
        d=-d
        do k=1,m
         t=a(j,k)
         a(j,k)=a(iquad,k)
         a(iquad,k)=t
        enddo
       endif
       ajj=a(j,j)
       if(ajj.eq.0.)then
        jm=j-1
        print *,'failure in lufact: matrix singular, rank=',jm
        stop
       endif
       ajji=1./ajj
       do i=jp,m
        aij=ajji*a(i,j)
        a(i,j)=aij
        do k=jp,m
         a(i,k)=a(i,k)-aij*a(j,k)
        enddo
       enddo
      enddo
      return
      end

C------------------------------------------------------------------------------C
C   R.J.Purser, National Meteorological Center, Washington D.C.  1993          C
C                   SUBROUTINE INVMM                                           C
C  invert matrix, possibly in place (a=b), using the l-u decomposition method  C
C  For DOUBLE PRECISION version see DINVMM                                     C
C                                                                              C
C  --> b    square matrix to be inverted                                       C
C  <-- a    inverse of b                                                       C
C  --> m    degree of (active part of) b and a                                 C
C  --> nb   first fortran dimension of b                                       C
C  --> na   first fortran dimension of a                                       C
C                                                                              C
C   LIMITATION:                                                                C
C    ipiv is an index array, internal to this array, encoding the              C
C    pivoting sequence used. It is given a fortran dimension of NN=500         C
C    in the parameter statement below. If the order of the linear system       C
C    exceeds 500, increase this parameter accordingly                          C
C                                                                              C
C------------------------------------------------------------------------------C
      subroutine INVMM(b,a,m,nb,na)
      PARAMETER (NN=500)
      DIMENSION IPIV(NN)
      dimension a(na,*),b(nb,*)
      do j=1,m
      do i=1,m
       a(i,j)=b(i,j)
      enddo
      enddo
      call lufm(a,ipiv,d,m,na)
c  invert u in place:
      do i=1,m
       a(i,i)=1./a(i,i)
      enddo
      do i=1,m-1
       do j=i+1,m
        s=0.
        do k=i,j-1
         s=s-a(i,k)*a(k,j)
        enddo
        a(i,j)=a(j,j)*s
       enddo
      enddo
c  invert l in place assuming implicitly diagonal elements of unity
      do j=1,m-1
       do i=j+1,m
        s=-a(i,j)
        do k=j+1,i-1
         s=s-a(i,k)*a(k,j)
        enddo
        a(i,j)=s
       enddo
      enddo
c  form the product of u**-1 and l**-1 in place
      do j=1,m-1
       do i=1,j
        s=a(i,j)
        do k=j+1,m
         s=s+a(i,k)*a(k,j)
        enddo
        a(i,j)=s
       enddo
       do i=j+1,m
        s=0.
        do k=i,m
         s=s+a(i,k)*a(k,j)
        enddo
        a(i,j)=s
       enddo
      enddo
c  permute columns according to ipiv
      do j=m-1,1,-1
       l=ipiv(j)
       do i=1,m
        t=a(i,j)
        a(i,j)=a(i,l)
        a(i,l)=t
       enddo
      enddo
      return
      end
!     Not in MPI version
!      include 'jimco.f'
!      include 'nfft.f'

