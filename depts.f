      subroutine depts1(x3d,y3d,z3d)  ! input ubar,vbar are unstaggered vels for level k
!     3D version
      implicit none
c     modify toij5 for Cray
      integer, parameter :: ntest=0
      include 'newmpar.h'
      include 'constant.h'   ! rearth
      include 'indices.h' ! in,is,iw,ie,inn,iss,iww,iee
      include 'map.h'
      include 'parm.h'
      include 'vecsuv.h'   ! vecsuv info
      include 'xyzinfo.h'  ! x,y,z,wts
      real ubar, vbar
      common/uvbar/ubar(ifull,kl),vbar(ifull,kl)
      integer nface
      real xg, yg
      common/work3f/nface(ifull,kl),xg(ifull,kl),yg(ifull,kl) ! depts, upglobal
!     Work common for these?
      real uc(ifull,kl),vc(ifull,kl),wc(ifull,kl), temp(ifull,kl)
      real x3d(ifull,kl),y3d(ifull,kl),z3d(ifull,kl)   ! upglobal depts 
      integer iq, k, intsch

      do k=1,kl
         do iq=1,ifull
c           departure point x, y, z is called x3d, y3d, z3d
c           first find corresponding cartesian vels
            uc(iq,k)=(ax(iq)*ubar(iq,k) + bx(iq)*vbar(iq,k))*dt/rearth ! unit sphere 
            vc(iq,k)=(ay(iq)*ubar(iq,k) + by(iq)*vbar(iq,k))*dt/rearth ! unit sphere 
            wc(iq,k)=(az(iq)*ubar(iq,k) + bz(iq)*vbar(iq,k))*dt/rearth ! unit sphere 
            x3d(iq,k)=x(iq)-uc(iq,k) ! 1st guess
            y3d(iq,k)=y(iq)-vc(iq,k)
            z3d(iq,k)=z(iq)-wc(iq,k)
         enddo                  ! iq loop
      end do

c     convert to grid point numbering
      do k=1,kl
         call toij5 (k,x3d(1,k),y3d(1,k),z3d(1,k)) ! maybe remove k dependency
      end do

      if(ntest.eq.1)then
        print *,'ubar,vbar ',ubar(idjd,nlv),vbar(idjd,nlv)
        print *,'uc,vc,wc ',uc(idjd,nlv),vc(idjd,nlv),wc(idjd,nlv)
        print *,'1st guess for k = ',nlv
        print *,'x3d,y3d,z3d ',x3d(idjd,nlv),y3d(idjd,nlv),z3d(idjd,nlv)
        print *,'xg,yg,nface ',xg(idjd,nlv),yg(idjd,nlv),nface(idjd,nlv)
      endif

      intsch=mod(ktau,2)
      temp = uc
      call ints(temp,intsch,nface,xg,yg,2)
      do k=1,kl
         x3d(:,k) = x(:) - 0.5*(uc(:,k)+temp(:,k)) ! 2nd guess
      end do
      temp = vc
      call ints(temp,intsch,nface,xg,yg,2)
      do k=1,kl
         y3d(:,k) = y(:) - 0.5*(vc(:,k)+temp(:,k)) ! 2nd guess
      end do
      temp = wc
      call ints(temp,intsch,nface,xg,yg,2)
      do k=1,kl
         z3d(:,k) = z(:) - 0.5*(wc(:,k)+temp(:,k)) ! 2nd guess
      end do

      do k=1,kl
         call toij5 (k,x3d(1,k),y3d(1,k),z3d(1,k)) ! maybe remove k dependency
      end do
      if(ntest.eq.1)then
        print *,'2nd guess for k = ',nlv
        print *,'x3d,y3d,z3d ',x3d(idjd,nlv),y3d(idjd,nlv),z3d(idjd,nlv)
        print *,'xg,yg,nface ',xg(idjd,nlv),yg(idjd,nlv),nface(idjd,nlv)
      endif

      temp = uc
      call ints(temp,intsch,nface,xg,yg,2)
      do k=1,kl
         x3d(:,k) = x(:) - 0.5*(uc(:,k)+temp(:,k)) ! 3rd guess
      end do
      temp = vc
      call ints(temp,intsch,nface,xg,yg,2)
      do k=1,kl
         y3d(:,k) = y(:) - 0.5*(vc(:,k)+temp(:,k)) ! 3rd guess
      end do
      temp = wc
      call ints(temp,intsch,nface,xg,yg,2)
      do k=1,kl
         z3d(:,k) = z(:) - 0.5*(wc(:,k)+temp(:,k)) ! 3rd guess
      end do

      do k=1,kl
         call toij5 (k,x3d(1,k),y3d(1,k),z3d(1,k)) ! maybe remove k dependency
      end do
      if(ntest.eq.1)then
        print *,'3rd guess for k = ',nlv
        print *,'x3d,y3d,z3d ',x3d(idjd,nlv),y3d(idjd,nlv),z3d(idjd,nlv)
        print *,'xg,yg,nface ',xg(idjd,nlv),yg(idjd,nlv),nface(idjd,nlv)
      endif

      return
      end

      subroutine depts(x3d,y3d,z3d)   ! input ubar,vbar are unstaggered vels for level k
!     3D version
c     modify toij5 for Cray
      implicit none
      integer, parameter :: ntest=0
      include 'newmpar.h'
      include 'constant.h'   ! rearth
      include 'indices.h' ! in,is,iw,ie,inn,iss,iww,iee
      include 'map.h'
      include 'parm.h'
      include 'vecsuv.h'   ! vecsuv info
      include 'xyzinfo.h'  ! x,y,z,wts
      real ubar, vbar
      common/uvbar/ubar(ifull,kl),vbar(ifull,kl)
      integer nface
      real xg, yg
      common/work3f/nface(ifull,kl),xg(ifull,kl),yg(ifull,kl) ! depts, upglobal
      ! Is there an appropriate work common for this?
      real gx(ifull,kl),gy(ifull,kl),gz(ifull,kl),
     &     derx(ifull,kl),dery(ifull,kl),derz(ifull,kl)
      real x3d(ifull,kl),y3d(ifull,kl),z3d(ifull,kl)
      integer nit, ndiag
      data nit/3/,ndiag/0/
      integer :: iq, itn, k
      real :: uc, vc, wc

      if(ntest.eq.1)then
        print *,'entering depts'
        print *,'ubar,vbar ',ubar(idjd,nlv),vbar(idjd,nlv)
      endif
      do k=1,kl
         do iq=1,ifull
!           departure point x, y, z is called x3d, y3d, z3d
!           first find corresponding cartesian vels
            uc = ax(iq)*ubar(iq,k) + bx(iq)*vbar(iq,k)
            vc = ay(iq)*ubar(iq,k) + by(iq)*vbar(iq,k)
            wc = az(iq)*ubar(iq,k) + bz(iq)*vbar(iq,k)
!           for gnomic can do first term analytically!
            derx(iq,k) = -uc*dt/rearth ! because working on unit sphere
            dery(iq,k) = -vc*dt/rearth
            derz(iq,k) = -wc*dt/rearth
            x3d(iq,k) = x(iq)+derx(iq,k)
            y3d(iq,k) = y(iq)+dery(iq,k)
            z3d(iq,k) = z(iq)+derz(iq,k)
!           for other terms want vels in units of double grid-points/timestep
            ubar(iq,k) = ubar(iq,k)*dt *.5*em(iq)/ds
            vbar(iq,k) = vbar(iq,k)*dt *.5*em(iq)/ds
         enddo ! iq loop
      end do
      if(ntest.eq.1)then
         print *,'itn,x3d,y3d,z3d,u,v: ',1,x3d(idjd,nlv),y3d(idjd,nlv),
     &        z3d(idjd,nlv),ubar(idjd,nlv),vbar(idjd,nlv)
         print *,'x,ax,bx,derx: ',x(idjd),ax(idjd),bx(idjd),
     &        derx(idjd,nlv)
      endif
!       Need a version that works on 3D arrays
!       if(ndiag.eq.2)then
!         call printp('derx',derx)
!         call printp('dery',dery)
!         call printp('derz',derz)
!         call printp('x3d ',x3d)
!         call printp('y3d ',y3d)
!         call printp('z3d ',z3d)
!       endif

      do itn=2,nit
         do k=1,kl
*cdir nodep
            do iq=1,ifull
               gx(iq,k) = -(ubar(iq,k)*(derx(ie(iq),k)-derx(iw(iq),k))
     &               + vbar(iq,k)*(derx(in(iq),k)-derx(is(iq),k)) )/itn
               gy(iq,k) = -(ubar(iq,k)*(dery(ie(iq),k)-dery(iw(iq),k))
     &               + vbar(iq,k)*(dery(in(iq),k)-dery(is(iq),k)) )/itn
               gz(iq,k) = -(ubar(iq,k)*(derz(ie(iq),k)-derz(iw(iq),k))
     &               + vbar(iq,k)*(derz(in(iq),k)-derz(is(iq),k)) )/itn
            enddo               ! iq loop
         end do

c        if(diag)print *,'depts itn,nit,dt: ',itn,nit,dt

         do k=1,kl
            do iq=1,ifull
               derx(iq,k) = gx(iq,k)
               dery(iq,k) = gy(iq,k)
               derz(iq,k) = gz(iq,k)
               x3d(iq,k) = x3d(iq,k)+derx(iq,k)
               y3d(iq,k) = y3d(iq,k)+dery(iq,k)
               z3d(iq,k) = z3d(iq,k)+derz(iq,k)
            enddo               ! iq loop
         end do

c        if(diag)print *,'itn,x3d,y3d,z3d,u,v: ',itn,x3d(idjd),y3d(idjd)
c    .                    ,z3d(idjd),ubar(idjd,1),vbar(idjd,1)
!         if(ndiag.eq.2)then
!c          call printp('gx   ',gx)
!c          call printp('gy   ',gy)
!c          call printp('gz   ',gz)
!           call printp('derx',derx)
!           call printp('dery',dery)
!           call printp('derz',derz)
!           call printp('x3d ',x3d)
!           call printp('y3d ',y3d)
!           call printp('z3d ',z3d)
!         endif
      enddo                     ! itn=2,nit

c     convert to grid point numbering
      do k=1,kl
         if(npanels.eq.5) call toij5 (k,x3d(1,k),y3d(1,k),z3d(1,k)) ! maybe remove k dependency
         if(npanels.eq.13)call toij13(k,x3d(1,k),y3d(1,k),z3d(1,k)) ! maybe remove k dependency
      end do

      if(ntest.eq.1)then
        print *,'at end of depts for k = ',k
        print *,'x3d,y3d,z3d ',x3d(idjd,nlv),y3d(idjd,nlv),z3d(idjd,nlv)
        print *,'xg,yg,nface ',xg(idjd,nlv),yg(idjd,nlv),nface(idjd,nlv)
      endif
      if(ndiag.eq.2)then
        print *,'after toij5/toij13'
        call printp('xg  ',xg)
        call printp('yg  ',yg)
      endif
      return
      end

      subroutine toij5(k,x3d,y3d,z3d)
c     modify toij5 for Cray
      parameter (ncray=0)    ! 0 for most computers, 1 for Cray
      parameter (ntest=0)
      include 'newmpar.h'
      include 'parm.h'
      include 'bigxy4.h' ! common/bigxy4/xx4(iquad,iquad),yy4(iquad,iquad)
      common/work3f/nface(ifull,kl),xg(ifull,kl),yg(ifull,kl) ! depts, upglobal
      real x3d(ifull),y3d(ifull),z3d(ifull)
      common/work2b/xstr(ifull),ystr(ifull),zstr(ifull)
      include 'xyzinfo.h'  ! x,y,z,wts
      dimension xgx(0:5),xgy(0:5),xgz(0:5),ygx(0:5),ygy(0:5),ygz(0:5)
      data xgx/0., 0., 0., 0., 1., 1./, ygx/0.,-1.,-1., 0., 0., 0./,
     .     xgy/1., 1., 0., 0., 0., 0./, ygy/0., 0., 0.,-1.,-1., 0./,
     .     xgz/0., 0.,-1.,-1., 0., 0./, ygz/1., 0., 0., 0., 0., 1./
      data nmaploop/3/,ndiag/0/,num/0/
      save num
      if(num.eq.0.and.ncray.eq.0)then  ! check if divide by itself is working
	do iq=1,100
	 xstr(iq)=-1.+.0002*iq
	 ystr(iq)=min(abs(xstr(iq)),4.)
	 zstr(iq)=xstr(iq)/ystr(iq)
	enddo 
	do iq=1,100
!       print *,'iq,xstr,ystr,zstr ',iq,xstr(iq),ystr(iq),zstr(iq)
	 if(zstr(iq).ne.-1.)stop 'must use ncray=1 on this PC'
	enddo     
       num=1
      endif

!     if necessary, transform (x3d, y3d, z3d) to equivalent
!     coordinates (xstr, ystr, zstr) on regular gnomonic panels
      if(schmidt.eq.1.)then
        do iq=1,ifull
         xstr(iq)=x3d(iq)
         ystr(iq)=y3d(iq)
         zstr(iq)=z3d(iq)
        enddo   ! iq loop
      else      ! (schmidt.ne.1.)
        alf=(1.-schmidt**2)/(1.+schmidt**2)
        alfonsch=(1.-alf)/schmidt
        do iq=1,ifull
         xstr(iq)=x3d(iq)*alfonsch/(1.-alf*z3d(iq))
         ystr(iq)=y3d(iq)*alfonsch/(1.-alf*z3d(iq))
         zstr(iq)=   (z3d(iq)-alf)/(1.-alf*z3d(iq))
        enddo   ! iq loop
      endif     ! (schmidt.ne.1.)

c      first deduce departure faces
c      instead calculate cubic coordinates
c      The faces are:
c      0: X=1   1: Z=1   2: Y=1   3: X=-1   4: Z=-1   5: Y=-1

       do iq=1,ifull
        denxyz=max( abs(xstr(iq)),abs(ystr(iq)),abs(zstr(iq)) )
        xd=xstr(iq)/denxyz
        yd=ystr(iq)/denxyz
        zd=zstr(iq)/denxyz

        if(ncray.eq.1)then
c         all these if statements are replaced by the subsequent cunning code
          if(abs(xstr(iq)).eq.denxyz)then         ! Cray
             if(xstr(iq).eq.denxyz)then           ! Cray
              nface(iq,k)    =0                   ! Cray
              xg(iq,k) =       yd                 ! Cray
              yg(iq,k) =       zd                 ! Cray
            else                                  ! Cray
              nface(iq,k)    =3                   ! Cray
              xg(iq,k) =     -zd                  ! Cray
              yg(iq,k) =     -yd                  ! Cray
            endif                               ! Cray
          elseif(abs(zstr(iq)).eq.denxyz)then    ! Cray
             if(zstr(iq).eq.denxyz)then           ! Cray
              nface(iq,k)    =1                   ! Cray
              xg(iq,k) =      yd                  ! Cray
              yg(iq,k) =     -xd                  ! Cray
            else                                ! Cray
              nface(iq,k)    =4                   ! Cray
              xg(iq,k) =      xd                  ! Cray
              yg(iq,k) =     -yd                  ! Cray
            endif                               ! Cray
          else                                  ! Cray
            if(ystr(iq).eq.denxyz)then           ! Cray
              nface(iq,k)    =2                   ! Cray
              xg(iq,k) =     -zd                  ! Cray
              yg(iq,k) =     -xd                  ! Cray
            else                                ! Cray
              nface(iq,k)    =5                   ! Cray
              xg(iq,k) =      xd                  ! Cray
              yg(iq,k) =      zd                  ! Cray
            endif                               ! Cray
          endif                                 ! Cray
        else  ! e.g. ncray=0
c         N.B. the Cray copes poorly with the following (sometimes .ne.1),
c              with e.g. division of  .978 by itself giving  .99999.....53453
c         max() allows for 2 of x,y,z being 1.  This is the cunning code:
          nf=max( int(xd)*(3*int(xd)-3) ,   ! ** n=0,5 version
     .            int(zd)*(5*int(zd)-3) ,
     .            int(yd)*(7*int(yd)-3) )/2
          nface(iq,k)=nf
          xg(iq,k)=xgx(nf)*xd+xgy(nf)*yd+xgz(nf)*zd  ! -1 to 1
          yg(iq,k)=ygx(nf)*xd+ygy(nf)*yd+ygz(nf)*zd
        endif    ! (ncray.eq.1)
       enddo   ! iq loop   
	if(ntest.eq.1)then
	 iq=idjd
        denxyz=max( abs(xstr(iq)),abs(ystr(iq)),abs(zstr(iq)) )
        xd=xstr(iq)/denxyz
        yd=ystr(iq)/denxyz
        zd=zstr(iq)/denxyz
	 print *,'xstr,ystr,zstr,denxyz ',
     .           xstr(iq),ystr(iq),zstr(iq),denxyz
        write(6,'("xstr,ystr,zstr = ",3z20)') xstr(iq),ystr(iq),zstr(iq)
        write(6,'("denxyz = ",z20)') denxyz
        write(6,'("xd,yd,zd = ",3z20)') xd,yd,zd
	 print *,'abs(xstr,ystr,zstr) ',
     .           abs(xstr(iq)),abs(ystr(iq)),abs(zstr(iq))
	 print *,'xd,yd,zd,nf ',xd,yd,zd,nf
       endif
       if(ndiag.eq.2)then
         call printp('xcub',xd)  ! need to reinstate as arrays for this diag
         call printp('ycub',yd)
         call printp('zcub',zd)
c       call printn('nfac',nface)
         print *,'before xytoiq'
         call printp('xg  ',xg)
         call printp('yg  ',yg)
       endif

c     use 4* resolution grid il --> 4*il
      do iq=1,ifull
       xg(iq,k)=min(max(-.99999,xg(iq,k)),.99999)
       yg(iq,k)=min(max(-.99999,yg(iq,k)),.99999)
c      first guess for ri, rj and nearest i,j
       ri=1.+(1.+xg(iq,k))*2*il
       rj=1.+(1.+yg(iq,k))*2*il
       do loop=1,nmaploop
        i=nint(ri)
        j=nint(rj)
        is=sign(1.,ri-i)
        js=sign(1.,rj-j)
c       predict new value for ri, rj
        dxx=xx4(i+is,j)-xx4(i,j)
        dyx=xx4(i,j+js)-xx4(i,j)
        dxy=yy4(i+is,j)-yy4(i,j)
        dyy=yy4(i,j+js)-yy4(i,j)
        den=dxx*dyy-dyx*dxy
        ri=i+is*((xg(iq,k)-xx4(i,j))*dyy-(yg(iq,k)-yy4(i,j))*dyx)/den
        rj=j+js*((yg(iq,k)-yy4(i,j))*dxx-(xg(iq,k)-xx4(i,j))*dxy)/den
       enddo  ! loop loop
c      expect xg, yg to range between .5 and il+.5
       xg(iq,k)=.25*(ri+3.) -.5  ! -.5 for stag; back to normal ri, rj defn
       yg(iq,k)=.25*(rj+3.) -.5  ! -.5 for stag
      enddo   ! iq loop
      return
      end

      subroutine toij13(k,x3d,y3d,z3d)  ! no Schmidt
      include 'newmpar.h'
      include 'parm.h'
      include 'bigxy4.h' ! common/bigxy4/xx4(iquad,iquad),yy4(iquad,iquad)
      common/work3f/nface(ifull,kl),xg(ifull,kl),yg(ifull,kl) ! depts, upglobal
      real x3d(ifull),y3d(ifull),z3d(ifull)
      include 'xyzinfo.h'  ! x,y,z,wts
      dimension npanetab(-1:1,-1:4),acon(-1:1,-1:4),bcon(-1:1,-1:4)
      dimension     xadd(-1:1,-1:4),yadd(-1:1,-1:4)
      data nmaploop/3/
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

c     first convert to equivalent of schmidt=.5 (or .1 better) grid
      alf=(1.-schmidt**2)/(1.+schmidt**2)  ! value of z on "equator"
!     schm13=.1  ! originally used .5 ! now comes from setxyz through parm.h
      schmidtp=schm13/schmidt        ! for z3d above the "equator"
      alfp=(1.-schmidtp**2)/(1.+schmidtp**2)
      schmidtm=1./(schm13*schmidt)   ! for z3d below the "equator"
      alfm=(1.-schmidtm**2)/(1.+schmidtm**2)
      do iq=1,ifull
       if(z3d(iq).gt.alf)then
         zsign=1.
         xstr=x3d(iq)*schmidtp*(1.+alfp)/(1.+alfp*z3d(iq))
         ystr=y3d(iq)*schmidtp*(1.+alfp)/(1.+alfp*z3d(iq))
       else
         zsign=-1.
         xstr=x3d(iq)*schmidtm*(1.+alfm)/(1.+alfm*z3d(iq))
         ystr=y3d(iq)*schmidtm*(1.+alfm)/(1.+alfm*z3d(iq))
       endif     !  (z3d(iq).gt.alf)
!      could avoid above "if", by first doing 1/schmidt, then schmidt13
!      using abs(z3d). Extra operations may then increase errors slightly?
c      now remember departure quadrants
       xsign=sign(1.,xstr)
       ysign=sign(1.,ystr)

c      use 4* resolution grid
c      N.B. for toij13, have xx4 and yy4 between 0 and .5 (after schmidtx)
       xgr=abs(xstr)
       ygr=abs(ystr)
c      first guess for ri, rj (1 to 6*il+1) and nearest i,j
       ri=1.+xgr*(iquad-1)/xx4(iquad,1)    ! divide by schm13 "equator" radius
       rj=1.+ygr*(iquad-1)/xx4(iquad,1)
       do loop=1,nmaploop
!       ri=max(1. , min(real(iquad)-.0001,ri) )  !  not needed
!       rj=max(1. , min(real(iquad)-.0001,rj) )  !  not needed
        i=nint(ri)
        j=nint(rj)
        is=sign(1.,ri-i)
        js=sign(1.,rj-j)
c       predict new value for ri, rj
        dxx=xx4(i+is,j)-xx4(i,j)
        dyx=xx4(i,j+js)-xx4(i,j)
        dxy=yy4(i+is,j)-yy4(i,j)
        dyy=yy4(i,j+js)-yy4(i,j)
        den=dxx*dyy-dyx*dxy
        ri=i+is*((xgr-xx4(i,j))*dyy-(ygr-yy4(i,j))*dyx)/den
        rj=j+js*((ygr-yy4(i,j))*dxx-(xgr-xx4(i,j))*dxy)/den
       enddo  ! loop loop
!      write ri,rj on a BIG grid (-1.5*il to 1.5*il, -1.5*il to 4.5*il)
!      where the Y variable is wrapping around the globe
!      and the "north pole" is now at (0.,0.)
       ri=.25*(ri-1.)*xsign
       rj=.25*( (rj-1.)*ysign*zsign + (1.-zsign)*real(6*il) )
!      allocate to a box (-1:1, -1:4)
        ibox=max(-1,min(nint(ri/il),1))   ! allows for -1.5 or 1.5
        jbox=max(-1,min(nint(rj/il),4))   ! allows for -1.5 or 4.5
c      convert  xg, yg ( .5 to il+.5) and nface
       nface(iq,k)=npanetab(ibox,jbox)
       xg(iq,k)=.5 +xadd(ibox,jbox)*real(il)
     .             +acon(ibox,jbox)*ri -bcon(ibox,jbox)*rj
       yg(iq,k)=.5 +yadd(ibox,jbox)*real(il)
     .             +bcon(ibox,jbox)*ri +acon(ibox,jbox)*rj
      enddo   ! iq loop
      return
      end
