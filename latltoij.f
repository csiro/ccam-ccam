      subroutine latltoij(rlongin,rlatin,rlong0,rlat0,schmidt,
     &                    xout,yout,nf,xx4,yy4,ik)
c     given a pair of latitudes and longitudes (in degrees),
c     returns i and j values on the conformal-cubic grid as
c     xout ranging between .5 and   il +.5, and
c     yout ranging between .5 and   il +.5
c     Note that the parallel version still returns xout, yout on the 
c     global grid.

c     modify for Cray; used by plotg.f and topgencc.f
      use utilities
      implicit none
      integer, parameter :: ncray=1    ! 0 for many computers, 1 for Cray
      integer, parameter :: ntest=0, numtst=-1
      include 'newmpar.h'
      include 'const_phys.h'
      include 'parm.h'
      include 'parmdyn.h'
      real, save :: rotpolei(3,3)
      real xgx(0:5),xgy(0:5),xgz(0:5),ygx(0:5),ygy(0:5),ygz(0:5)
      data xgx/0., 0., 0., 0., 1., 1./, ygx/0.,-1.,-1., 0., 0., 0./,
     .     xgy/1., 1., 0., 0., 0., 0./, ygy/0., 0., 0.,-1.,-1., 0./,
     .     xgz/0., 0.,-1.,-1., 0., 0./, ygz/1., 0., 0., 0., 0., 1./

      real*8 xx4(1+4*ik,1+4*ik),yy4(1+4*ik,1+4*ik)
      real*8 dxx,dyy,dxy,dyx,denxyz,x,y,z 
      real*8, save :: alf, one=1.
      real den,ri,rj,xout,yout,xa,ya,za,xgrid,ygrid,xx,yy,zz,x1,z1
      real rlatin,rlongin,rlong0,rlat0,schmidt 
      integer ig,jg,ik,is,js,loop,nf
      integer, save :: nmaploop=3, num=0

!     if(num==0)then     ! not with onthefly
        alf=(one-schmidt**2)/(one+schmidt**2)
        rotpolei = transpose(calc_rotpole(rlong0,rlat0))
!     endif
      num=num+1
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
