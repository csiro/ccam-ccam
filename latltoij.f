      subroutine latltoij(rlongin,rlatin,xout,yout,nf)
c     given a pair of latitudes and longitudes (in degrees),
c     returns i and j values on the conformal-cubic grid as
c     xout ranging between .5 and   il +.5, and
c     yout ranging between .5 and 6*il +.5

c     modify for Cray; used by plotg.f and topgencc.f
      parameter (ncray=0)    ! 0 for most computers, 1 for Cray
c     contains a version of xytoij
      parameter (ntest=0)
      include 'newmpar.h'
      include 'bigxy4.h' ! common/bigxy4/xx4(iquad,iquad),yy4(iquad,iquad)
      include 'const_phys.h'
      include 'parm.h'
      include 'parmdyn.h'
      real rotpolei(3,3)
      dimension xgx(0:5),xgy(0:5),xgz(0:5),ygx(0:5),ygy(0:5),ygz(0:5)
      data xgx/0., 0., 0., 0., 1., 1./, ygx/0.,-1.,-1., 0., 0., 0./,
     .     xgy/1., 1., 0., 0., 0., 0./, ygy/0., 0., 0.,-1.,-1., 0./,
     .     xgz/0., 0.,-1.,-1., 0., 0./, ygz/1., 0., 0., 0., 0., 1./

!     following used by npanels=13
      dimension npanetab(-1:1,-1:4),acon(-1:1,-1:4),bcon(-1:1,-1:4)
      dimension     xadd(-1:1,-1:4),yadd(-1:1,-1:4)
      data npanetab/9,12,13, 7,2,1, 6,4,3, 6,5,3, 8,10,0, 9,11,13/
      data acon/0,1,1, -1,0,0, -1,0,0, 0,0,-1, 1,1,0, 1,1,0/
      data bcon/1,0,0, 0,-1,-1, 0,-1,-1, -1,-1,0, 0,0,1, 0,0,1/
      data xadd/ -.5,  .5, -.5,   -.5,  .5,  .5,    -.5,-.5,-.5,
     .          -1.5,-1.5, 1.5,    1.5, .5, 3.5,    1.5, .5, 4.5/
      data yadd/ 1.5, 1.5, 1.5,     .5, .5, 1.5,    1.5, .5, 1.5,
     .           -.5,  .5, 2.5,   -2.5,-2.5,-.5,   -3.5,-3.5,-.5/
c        acon:   0    -1    -1       0     1     1     9  7  6   6  8  9
c                1     0     0       0     1     1    12  2  4   5 10 11
c                1     0     0      -1     0     0    13  1  3   3  0 13
c        bcon:   1     0     0      -1     0     0
c                0    -1    -1      -1     0     0
c                0    -1    -1       0     1     1
c        xadd: -0.5  -0.5  -0.5    -1.5   1.5   1.5      correct ones
c               0.5   0.5  -0.5    -1.5   0.5   0.5
c              -0.5   0.5  -0.5     1.5   3.5   4.5
c        yadd:  1.5   0.5   1.5    -0.5  -2.5  -3.5
c               1.5   0.5   0.5     0.5  -2.5  -3.5
c               1.5   1.5   1.5     2.5  -0.5  -0.5

      save num,rotpolei,alf,schmidtp,alfp,schmidtm,alfm
      data nmaploop/3/,numtst/-1/,num/0/
c     print *,'entering latltoij;num,rlongin,rlatin,schmidt,schm13 '
c    .                          ,num,rlongin,rlatin,schmidt,schm13
!     if(num.eq.0)then     ! not with onthefly
        alf=(1.-schmidt**2)/(1.+schmidt**2)
!       schm13=.1  ! originally used .5 ! now comes from setxyz through parm.h
!       alf13=(1.-schm13**2)/(1.+schm13**2)
!       schm13rad=schm13*(1.+alf13)
        schmidtp=schm13/schmidt        ! for z3d above the "equator"
        alfp=(1.-schmidtp**2)/(1.+schmidtp**2)
        schmidtm=1./(schm13*schmidt)   ! for z3d below the "equator"
        alfm=(1.-schmidtm**2)/(1.+schmidtm**2)
c       print *,'latltoij; rlong0,rlat0,schmidt,schm13,alf:',
c    .                     rlong0,rlat0,schmidt,schm13,alf
c       print *,'xx4(iquad,1) ',xx4(iquad,1)
        coslong=cos(rlong0*pi/180.)
        sinlong=sin(rlong0*pi/180.)
        coslat=cos(rlat0*pi/180.)
        sinlat=sin(rlat0*pi/180.)
        rotpolei(1,1)=coslong*sinlat
        rotpolei(2,1)=-sinlong
        rotpolei(3,1)=coslong*coslat
        rotpolei(1,2)=sinlong*sinlat
        rotpolei(2,2)=coslong
        rotpolei(3,2)=sinlong*coslat
        rotpolei(1,3)=-coslat
        rotpolei(2,3)=0.
        rotpolei(3,3)=sinlat
!     endif
      num=num+1
c     numtst=num
      if(num.eq.numtst)print *,'a rlongin,rlatin ',rlongin,rlatin
      xa=cos(rlongin*pi/180.)*cos(rlatin*pi/180.)
      ya=sin(rlongin*pi/180.)*cos(rlatin*pi/180.)
      za=sin(rlatin*pi/180.)
      if(num.eq.numtst)print *,'b xa,ya,za ',xa,ya,za
      x=rotpolei(1,1)*xa+rotpolei(1,2)*ya+rotpolei(1,3)*za
      y=rotpolei(2,1)*xa+rotpolei(2,2)*ya+rotpolei(2,3)*za
      z=rotpolei(3,1)*xa+rotpolei(3,2)*ya+rotpolei(3,3)*za
      if(num.eq.numtst)print *,'c x,y,z ',x,y,z

      if(npanels.eq.5)then

!       if necessary, transform physical (x, y, z) to equivalent coordinates
!       on regular gnomonic panels
        if(schmidt.ne.1.)then
          x1=x

          z1=z
          x=x*(1.-alf)/(schmidt*(1.-alf*z))
          y=y*(1.-alf)/(schmidt*(1.-alf*z))
          z=(z-alf)/(1.-alf*z)
          if(ntest.eq.1.and.z1.gt..82.and.z1.lt..821)then
            print *,'latltoij: rlongin, rlatin ',rlongin, rlatin
            print *,'latltoij: xa,ya,za ',xa,ya,za
            print *,'latltoij: x1,z1 ',x1,z1
            print *,'latltoij: x,y,z ',x,y,z
          endif
        endif         ! (schmidt.ne.1.)

        denxyz=max( abs(x),abs(y),abs(z) )
        xx=x/denxyz
        yy=y/denxyz
        zz=z/denxyz
c       deduce corresponding face
        if(ncray.eq.1)then
c         all these if statements are replaced by the subsequent cunning code
          if(abs(x ).eq.denxyz)then             ! Cray
            if(x .eq.denxyz)then                ! Cray
              nf    =0                          ! Cray
              xgrid =       yy                  ! Cray
              ygrid =       zz                  ! Cray
            else                                ! Cray
              nf    =3                          ! Cray
              xgrid =     -zz                   ! Cray
              ygrid =     -yy                   ! Cray
            endif                               ! Cray
          elseif(abs(z ).eq.denxyz)then         ! Cray
            if(z .eq.denxyz)then                ! Cray
              nf    =1                          ! Cray
              xgrid =      yy                   ! Cray
              ygrid =     -xx                   ! Cray
            else                                ! Cray
              nf    =4                          ! Cray
              xgrid =      xx                   ! Cray
              ygrid =     -yy                   ! Cray
            endif                               ! Cray
          else                                  ! Cray
            if(y .eq.denxyz)then                ! Cray
              nf    =2                          ! Cray
              xgrid =     -zz                   ! Cray
              ygrid =     -xx                   ! Cray
            else                                ! Cray
              nf    =5                          ! Cray
              xgrid =      xx                   ! Cray
              ygrid =      zz                   ! Cray
            endif                               ! Cray
          endif                                 ! Cray
        else  ! e.g. ncray=0
          nf=max( int(xx)*(3*int(xx)-3) ,   ! ** n=0,5 version
     .            int(zz)*(5*int(zz)-3) ,
     .            int(yy)*(7*int(yy)-3) )/2
          xgrid=xgx(nf)*xx+xgy(nf)*yy+xgz(nf)*zz  ! -1 to 1
          ygrid=ygx(nf)*xx+ygy(nf)*yy+ygz(nf)*zz
        endif    ! (ncray.eq.1)

c       convert to grid point numbering
c       the xytoij routine follows

c       use 4* resolution grid il --> 4*il
        xgrid=min(max(-.99999,xgrid),.99999)
        ygrid=min(max(-.99999,ygrid),.99999)
c       first guess for ri, rj and nearest ig,jg
        ri=1.+(1.+xgrid)*2*il
        rj=1.+(1.+ygrid)*2*il
        do loop=1,nmaploop
         ig=nint(ri)
         jg=nint(rj)
         is=sign(1.,ri-ig)
         js=sign(1.,rj-jg)
c        predict new value for ri, rj
         dxx=xx4(ig+is,jg)-xx4(ig,jg)
         dyx=xx4(ig,jg+js)-xx4(ig,jg)
         dxy=yy4(ig+is,jg)-yy4(ig,jg)
         dyy=yy4(ig,jg+js)-yy4(ig,jg)
         den=dxx*dyy-dyx*dxy
         ri=ig+is*((xgrid-xx4(ig,jg))*dyy-(ygrid-yy4(ig,jg))*dyx)/den
         rj=jg+js*((ygrid-yy4(ig,jg))*dxx-(xgrid-xx4(ig,jg))*dxy)/den
        enddo  ! loop loop
        xout=.25*(ri+3.) -.5  ! -.5 for stag; back to normal ri, rj defn
        yout=.25*(rj+3.) -.5  ! -.5 for stag
c       expect xout, yout (at this point) to range between .5 and il+.5

      elseif(npanels.eq.13)then

c       first convert to equivalent of schmidt=.5 grid
        if(z.gt.alf)then
          zsign=1.
          xstr=x*schmidtp*(1.+alfp)/(1.+alfp*z)
          ystr=y*schmidtp*(1.+alfp)/(1.+alfp*z)
        else
          zsign=-1.
          xstr=x*schmidtm*(1.+alfm)/(1.+alfm*z)
          ystr=y*schmidtm*(1.+alfm)/(1.+alfm*z)
        endif     !  (z.gt.alf)
!       could avoid above "if", by first doing 1/schmidt, then schmidt13
!       using abs(z). Extra operations may then increase errors slightly?
c       now remember departure quadrants
        xsign=sign(1.,xstr)
        ysign=sign(1.,ystr)
        if(num.eq.numtst)print *,'d xstr,ystr,xsign,ysign,zsign ',
     .                              xstr,ystr,xsign,ysign,zsign

c       use 4* resolution grid
c       N.B. for toij13, have xx4 and yy4 between 0 and .5 (after schmidtx )
        xgr=abs(xstr)
        ygr=abs(ystr)
c       first guess for ri, rj (1 to 6*il+1) and nearest i,j
        ri=1.+xgr*(iquad-1)/xx4(iquad,1)    ! divide by schm13 "equator" radius
        rj=1.+ygr*(iquad-1)/xx4(iquad,1)
        if(num.eq.numtst)print *,'e xgr,ygr,ri,rj',xgr,ygr,ri,rj
        do loop=1,nmaploop
!        ri=max(1. , min(real(iquad)-.0001,ri) )  !  not needed
!        rj=max(1. , min(real(iquad)-.0001,rj) )  !  not needed
         i=nint(ri)
         j=nint(rj)
         is=sign(1.,ri-i)
         js=sign(1.,rj-j)
c        predict new value for ri, rj
         dxx=xx4(i+is,j)-xx4(i,j)
         dyx=xx4(i,j+js)-xx4(i,j)
         dxy=yy4(i+is,j)-yy4(i,j)
         dyy=yy4(i,j+js)-yy4(i,j)
         den=dxx*dyy-dyx*dxy
         ri=i+is*((xgr-xx4(i,j))*dyy-(ygr-yy4(i,j))*dyx)/den
         rj=j+js*((ygr-yy4(i,j))*dxx-(xgr-xx4(i,j))*dxy)/den
         if(num.eq.numtst)print *,'e1,i,j,is,js,ri,rj ',
     .                                i,j,is,js,ri,rj
        enddo  ! loop loop
!       write ri,rj on a BIG grid (-1.5*il to 1.5*il, -1.5*il to 4.5*il)
!       where the Y variable is wrapping around the globe
!       and the "north pole" is now at (0.,0.)
        ri=.25*(ri-1.)*xsign
        rj=.25*( (rj-1.)*ysign*zsign + (1.-zsign)*real(6*il) )
        if(num.eq.numtst)print *,'e2 bigy  ri,rj ',ri,rj
!       allocate to a box (-1:1, -1:4)
        ibox=max(-1,min(nint(ri/il),1))   ! allows for -1.5 or 1.5
        jbox=max(-1,min(nint(rj/il),4))   ! allows for -1.5 or 4.5
c       convert  xg, yg ( .5 to il+.5)
        if(num.eq.numtst)print *,'f ri,rj,ibox,jbox',ri,rj,ibox,jbox
        if(num.eq.numtst)print *,'f1 nf,xadd,yadd,acon,bcon ',
     .   npanetab(ibox,jbox),
     .   xadd(ibox,jbox),yadd(ibox,jbox),acon(ibox,jbox),bcon(ibox,jbox)
        nf=npanetab(ibox,jbox)
        xout=.5 +xadd(ibox,jbox)*real(il)
     .          +acon(ibox,jbox)*ri -bcon(ibox,jbox)*rj
        yout=.5 +yadd(ibox,jbox)*real(il)
     .          +bcon(ibox,jbox)*ri +acon(ibox,jbox)*rj
        if(num.eq.numtst)print *,'f1 xadd,yadd,acon,bcon ',
     .   xadd(ibox,jbox),yadd(ibox,jbox),acon(ibox,jbox),bcon(ibox,jbox)
        if(num.eq.numtst)print *,'g nf,xout,yout ',nf,xout,yout

      endif  !  (npanels.eq.5) elseif(npanels.eq.13)

      if(ntest.eq.1.and.rlongin.gt.43.1.and.rlongin.lt.49.9)then
        if(rlatin.gt.-24.2.and.rlatin.lt.-23.8)then
          print *,'lat,long,x,y,z,den ',rlatin,rlongin,x,y,z,denxyz
          print *,'xx,yy,zz ',xx,yy,zz
          print *,'nf,xout,yout,ri,rj ',nf,xout,yout,ri,rj
          print *,'youtb ',yout+nf*il
        endif
      endif
      return
      end
