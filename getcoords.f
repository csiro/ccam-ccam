      program getcoords

!  Get CC mdoel grid coordinates of specified lat/lon

!  It uses the model routine setxyz

      include 'newmpar.h'
      include 'parm.h'
      parameter(degtorad=3.1415926/180.)

      print*, ' Enter lat/lon (deg) '
      read(*,*) rlat, rlong
      rlat = rlat * degtorad
      rlong = rlong * degtorad
      call setxyz
      call xytoijf(rlat,rlong)

      end

!-------------------------------------------------------------------------
      subroutine xytoijf(rlat,rlong)
!     Taken from plotg.f, modified for the model history
      include 'newmpar.h'
      include 'bigxy4.h' ! common/bigxy4/xx4(4*il+1,4*il+1),yy4(4*il+1,4*il+1)
      dimension xgx(0:5),xgy(0:5),xgz(0:5),ygx(0:5),ygy(0:5),ygz(0:5)
      data xgx/0., 0., 0., 0., 1., 1./, ygx/0.,-1.,-1., 0., 0., 0./,
     .     xgy/1., 1., 0., 0., 0., 0./, ygy/0., 0., 0.,-1.,-1., 0./,
     .     xgz/0., 0.,-1.,-1., 0., 0./, ygz/1., 0., 0., 0., 0., 1./
      data nmaploop/3/
      data pi/3.1415926536/

      if (rlong.gt.pi) rlong = rlong-2.*pi !  -pi< rlong <= pi
!           set corresponding cartesian x,y,z, value on unit sphere
      x=cos(rlong)*cos(rlat)
      y=sin(rlong)*cos(rlat)
      z=sin(rlat)
      den=max( abs(x),abs(y),abs(z) )
      x=x/den
      y=y/den
      z=z/den
!           deduce corresponding face
      nf=max( int(x)*(3*int(x)-3) , ! ** n=0,5 version
     &              int(z)*(5*int(z)-3) ,
     &              int(y)*(7*int(y)-3) )/2
      xgrid=xgx(nf)*x+xgy(nf)*y+xgz(nf)*z ! -1 to 1
      ygrid=ygx(nf)*x+ygy(nf)*y+ygz(nf)*z

!           convert to grid point numbering
!           the xytoij routine follows

!           use 4* resolution grid il --> 4*il
      xgrid=min(max(-.99999,xgrid),.99999)
      ygrid=min(max(-.99999,ygrid),.99999)
!           first guess for ri, rj and nearest ig,jg
      ri=1.+(1.+xgrid)*2*il
      rj=1.+(1.+ygrid)*2*il
      do loop=1,nmaploop
         ig=nint(ri)
         jg=nint(rj)
         is=sign(1.,ri-ig)
         js=sign(1.,rj-jg)
!              predict new value for ri, rj
         dxx=xx4(ig+is,jg)-xx4(ig,jg)
         dyx=xx4(ig,jg+js)-xx4(ig,jg)
         dxy=yy4(ig+is,jg)-yy4(ig,jg)
         dyy=yy4(ig,jg+js)-yy4(ig,jg)
         den=dxx*dyy-dyx*dxy
         ri=ig+is*((xgrid-xx4(ig,jg))*dyy-
     &              (ygrid-yy4(ig,jg))*dyx)/den
         rj=jg+js*((ygrid-yy4(ig,jg))*dxx-
     &              (xgrid-xx4(ig,jg))*dxy)/den
      enddo                     ! loop loop
      xg=.25*(ri+3.)
      yg=.25*(rj+3.)
      print*, ' NFACE ', nf
      print*, ' XG, YG ', xg, yg
      print*, ' RI, RJ ', ri, rj

      return
      end
