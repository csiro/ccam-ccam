      subroutine latltoij(rlongin,rlatin,xout,yout,nf,xx4,yy4,ik)
c     given a pair of latitudes and longitudes (in degrees),
c     returns i and j values on the conformal-cubic grid as
c     xout ranging between .5 and   il +.5, and
c     yout ranging between .5 and   il +.5
c     Note that the parallel version still returns xout, yout on the 
c     global grid.

c     modify for Cray; used by plotg.f and topgencc.f
      use utilities
      parameter (ncray=1)    ! 0 for most computers, 1 for Cray
c     contains a version of xytoij
      parameter (ntest=1)
      include 'newmpar.h'
      real*8 xx4(1+4*ik,1+4*ik),yy4(1+4*ik,1+4*ik)
      real*8 dxx,dyy,dxy,dyx,denxyz
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
      real*8 alf,one,x,y,z  ! 6/11/07
      data one/1./    ! just to force real*8 calculation
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

      save num,rotpolei,alf
      data nmaploop/3/,numtst/100/,num/0/
!     if(num.eq.0)then     ! not with onthefly
        alf=(one-schmidt**2)/(one+schmidt**2)
        rotpolei = transpose(calc_rotpole(rlong0,rlat0))
!     endif
      num=num+1
c     numtst=num
      if(num<numtst)print *,'a rlongin,rlatin ',rlongin,rlatin
      xa=cos(rlongin*pi/180.)*cos(rlatin*pi/180.)
      ya=sin(rlongin*pi/180.)*cos(rlatin*pi/180.)
      za=sin(rlatin*pi/180.)
      if(num<numtst)print *,'b xa,ya,za ',xa,ya,za
      x=rotpolei(1,1)*xa+rotpolei(1,2)*ya+rotpolei(1,3)*za
      y=rotpolei(2,1)*xa+rotpolei(2,2)*ya+rotpolei(2,3)*za
      z=rotpolei(3,1)*xa+rotpolei(3,2)*ya+rotpolei(3,3)*za
      if(num<numtst)print *,'c x,y,z ',x,y,z

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
!       ri=1.+(1.+xgrid)*2*il_g
!       rj=1.+(1.+ygrid)*2*il_g
        ri=1.+(1.+xgrid)*2*ik
        rj=1.+(1.+ygrid)*2*ik
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

      if(ntest.eq.1.and.rlongin.gt.43.1.and.rlongin.lt.49.9)then
        if(rlatin.gt.-24.2.and.rlatin.lt.-23.8)then
          print *,'lat,long,x,y,z,den ',rlatin,rlongin,x,y,z,denxyz
          print *,'xx,yy,zz ',xx,yy,zz
          print *,'nf,xout,yout,ri,rj ',nf,xout,yout,ri,rj
!         print *,'youtb ',yout+nf*il_g
          print *,'youtb ',yout+nf*ik
        endif
      endif
      return
      end
