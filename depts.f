      subroutine depts1(x3d,y3d,z3d)  ! input ubar,vbar are unstaggered vels for level k
!     3D version
      use cc_mpi
      implicit none
c     modify toij5 for Cray
      integer, parameter :: ntest=0
      include 'newmpar.h'
      include 'const_phys.h'   ! rearth
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
!     temp needs iextra becaue it's used in ints
      real uc(ifull,kl),vc(ifull,kl),wc(ifull,kl), temp(ifull+iextra,kl)
      real*8 x3d(ifull,kl),y3d(ifull,kl),z3d(ifull,kl)   ! upglobal depts 
      integer iq, k, intsch

      call start_log(depts_begin)
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
!     Share off processor departure points.
      call deptsync(nface,xg,yg)

      if(diag.and.mydiag)then
        print *,'ubar,vbar ',ubar(idjd,nlv),vbar(idjd,nlv)
        print *,'uc,vc,wc ',uc(idjd,nlv),vc(idjd,nlv),wc(idjd,nlv)
        print *,'x,y,z ',x(idjd),y(idjd),z(idjd)
        
        print *,'1st guess for k = ',nlv
        print *,'x3d,y3d,z3d ',x3d(idjd,nlv),y3d(idjd,nlv),z3d(idjd,nlv)
        print *,'xg,yg,nface ',xg(idjd,nlv),yg(idjd,nlv),nface(idjd,nlv)
      endif

      intsch=mod(ktau,2)
      temp(1:ifull,:) = uc(1:ifull,:)
      call ints(temp,intsch,nface,xg,yg,2)
      do k=1,kl
         x3d(1:ifull,k) = x(1:ifull) -
     &                    0.5*(uc(1:ifull,k)+temp(1:ifull,k)) ! 2nd guess
      end do
      temp(1:ifull,:) = vc(1:ifull,:)
      call ints(temp,intsch,nface,xg,yg,2)
      do k=1,kl
         y3d(1:ifull,k) = y(1:ifull) -
     &                    0.5*(vc(1:ifull,k)+temp(1:ifull,k)) ! 2nd guess
      end do
      temp(1:ifull,:) = wc(1:ifull,:)
      call ints(temp,intsch,nface,xg,yg,2)
      do k=1,kl
         z3d(1:ifull,k) = z(1:ifull) -
     &                    0.5*(wc(1:ifull,k)+temp(1:ifull,k)) ! 2nd guess
      end do

      do k=1,kl
         call toij5 (k,x3d(1,k),y3d(1,k),z3d(1,k)) ! maybe remove k dependency
      end do
!     Share off processor departure points.
      call deptsync(nface,xg,yg)
      if(diag.and.mydiag)then
        print *,'2nd guess for k = ',nlv
        print *,'x3d,y3d,z3d ',x3d(idjd,nlv),y3d(idjd,nlv),z3d(idjd,nlv)
        print *,'xg,yg,nface ',xg(idjd,nlv),yg(idjd,nlv),nface(idjd,nlv)
      endif

      temp(1:ifull,:) = uc(1:ifull,:)
      call ints(temp,intsch,nface,xg,yg,2)
      do k=1,kl
         x3d(1:ifull,k) = x(1:ifull) -
     &                    0.5*(uc(1:ifull,k)+temp(1:ifull,k)) ! 3rd guess
      end do
      temp(1:ifull,:) = vc(1:ifull,:)
      call ints(temp,intsch,nface,xg,yg,2)
      do k=1,kl
         y3d(1:ifull,k) = y(1:ifull) -
     &                    0.5*(vc(1:ifull,k)+temp(1:ifull,k)) ! 3rd guess
      end do
      temp(1:ifull,:) = wc(1:ifull,:)
      call ints(temp,intsch,nface,xg,yg,2)
      do k=1,kl
         z3d(1:ifull,k) = z(1:ifull) -
     &                    0.5*(wc(1:ifull,k)+temp(1:ifull,k)) ! 3rd guess
      end do

      do k=1,kl
         call toij5 (k,x3d(1,k),y3d(1,k),z3d(1,k)) ! maybe remove k dependency
      end do
      if(diag.and.mydiag)then
        print *,'3rd guess for k = ',nlv
        print *,'x3d,y3d,z3d ',x3d(idjd,nlv),y3d(idjd,nlv),z3d(idjd,nlv)
        print *,'xg,yg,nface ',xg(idjd,nlv),yg(idjd,nlv),nface(idjd,nlv)
      endif

!     Share off processor departure points.
      call deptsync(nface,xg,yg)

      call end_log(depts_end)
      return
      end

      subroutine depts(x3d,y3d,z3d)   ! input ubar,vbar are unstaggered vels for level k
!     3D version
c     modify toij5 for Cray
      use cc_mpi
      implicit none
      integer, parameter :: ntest=0
      include 'newmpar.h'
      include 'const_phys.h'   ! rearth
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
      real, dimension(ifull,kl) :: gx, gy, gz
      real, dimension(ifull+iextra,kl) :: derx, dery, derz
      real*8 x3d(ifull,kl),y3d(ifull,kl),z3d(ifull,kl)
      integer, parameter :: nit=3, ndiag=0
      integer :: iq, itn, k
      real :: uc, vc, wc

      call start_log(depts_begin)
      if(ntest.eq.1.and.mydiag)then
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
         call bounds(derx)
         call bounds(dery)
         call bounds(derz)
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
!         Not implemented in the MPI version.
!         if(npanels.eq.13)call toij13(k,x3d(1,k),y3d(1,k),z3d(1,k)) ! maybe remove k dependency
      end do

      if(ntest.eq.1.and.mydiag)then
        print *,'at end of depts for k = ',k
        print *,'x3d,y3d,z3d ',x3d(idjd,nlv),y3d(idjd,nlv),z3d(idjd,nlv)
        print *,'xg,yg,nface ',xg(idjd,nlv),yg(idjd,nlv),nface(idjd,nlv)
      endif
!     printp not implemented in MPI version.
!      if(ndiag.eq.2)then
!        print *,'after toij5/toij13'
!        call printp('xg  ',xg)
!        call printp('yg  ',yg)
!      endif

!     Share off processor departure points.
      call deptsync(nface,xg,yg)

      call end_log(depts_end)
      return
      end

      subroutine toij5(k,x3d,y3d,z3d)
      use cc_mpi
c     modify toij5 for Cray
      parameter (ncray=1)    ! 0 for most computers, 1 for Cray
      parameter (ntest=0)
      include 'newmpar.h'
      include 'parm.h'
      include 'parmgeom.h'  ! rlong0,rlat0,schmidt  
      include 'bigxy4.h' ! common/bigxy4/xx4(iquad,iquad),yy4(iquad,iquad)
      common/work3f/nface(ifull,kl),xg(ifull,kl),yg(ifull,kl) ! depts, upglobal
      real*8 x3d(ifull),y3d(ifull),z3d(ifull)
      common/work2b/xstr(ifull),ystr(ifull),zstr(ifull)
      real*8 alf,alfonsch,den,one  ! 6/11/07 esp for 200m
      data one/1./             ! to force real*8
      include 'xyzinfo.h'  ! x,y,z,wts
      dimension xgx(0:5),xgy(0:5),xgz(0:5),ygx(0:5),ygy(0:5),ygz(0:5)
      data xgx/0., 0., 0., 0., 1., 1./, ygx/0.,-1.,-1., 0., 0., 0./,
     .     xgy/1., 1., 0., 0., 0., 0./, ygy/0., 0., 0.,-1.,-1., 0./,
     .     xgz/0., 0.,-1.,-1., 0., 0./, ygz/1., 0., 0., 0., 0., 1./
      data nmaploop/3/,ndiag/0/,num/0/
      save num

      call start_log(toij_begin)
      if(num==0)then
        if(mydiag)print *,'checking for ncray = ',ncray
        If(ncray==0)then  ! check if divide by itself is working
          call checkdiv(xstr,ystr,zstr)
        endif
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
        alf=(one-schmidt**2)/(one+schmidt**2)
!       alfonsch=(1.-alf)/schmidt
        alfonsch=2.*schmidt/(one+schmidt**2)  ! same but bit more accurate
        do iq=1,ifull
         den=one-alf*z3d(iq) ! to force real*8
         xstr(iq)=x3d(iq)*(alfonsch/den)
         ystr(iq)=y3d(iq)*(alfonsch/den)
         zstr(iq)=   (z3d(iq)-alf)/den
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
	if(ntest==1.and.k==nlv)then
	 iq=idjd
	 print *,'x3d,y3d,z3d ',x3d(iq),y3d(iq),z3d(iq)
         den=one-alf*z3d(iq) ! to force real*8
         print *,'den ',den
         denxyz=max( abs(xstr(iq)),abs(ystr(iq)),abs(zstr(iq)) )
         xd=xstr(iq)/denxyz
         yd=ystr(iq)/denxyz
         zd=zstr(iq)/denxyz
	 print *,'k,xstr,ystr,zstr,denxyz ',
     .            k,xstr(iq),ystr(iq),zstr(iq),denxyz
	 print *,'abs(xstr,ystr,zstr) ',
     .            abs(xstr(iq)),abs(ystr(iq)),abs(zstr(iq))
	 print *,'xd,yd,zd,nface ',xd,yd,zd,nface(iq,k)
         print *,'alf,alfonsch ',alf,alfonsch
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
c      if(iq<10.and.k==nlv)print *,'iq,xg,yg ',iq,xg(iq,k),yg(iq,k)
       xg(iq,k)=min(max(-.99999,xg(iq,k)),.99999)
       yg(iq,k)=min(max(-.99999,yg(iq,k)),.99999)
c      first guess for ri, rj and nearest i,j
       ri=1.+(1.+xg(iq,k))*2*il_g
       rj=1.+(1.+yg(iq,k))*2*il_g
       if(ntest==1.and.iq==idjd.and.k==nlv)then
         print *,'A: xg,yg,ri,rj ',xg(iq,k),yg(iq,k),ri,rj
       endif
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
c       if(iq<10.and.k==nlv)then
c        print *,'i,j,is,js,ri,rj ',i,j,is,js,ri,rj
c        print *,'xx4,yy4 ',xx4(i,j),yy4(i,j)
c        print *,'xx4a,xx4b ',xx4(i+is,j),xx4(i,j+js)
c        print *,'yy4a,yy4b ',yy4(i+is,j),yy4(i,j+js)
c        print *,'loop,dxx,dyy,dyx,dxy,den ',loop,dxx,dyy,dyx,dxy,den
c       endif
        ri=i+is*((xg(iq,k)-xx4(i,j))*dyy-(yg(iq,k)-yy4(i,j))*dyx)/den
        rj=j+js*((yg(iq,k)-yy4(i,j))*dxx-(xg(iq,k)-xx4(i,j))*dxy)/den
        if(ntest==1.and.iq==idjd.and.k==nlv)then
          print *,'B: xg,yg,ri,rj ',xg(iq,k),yg(iq,k),ri,rj
          print *,'i,j,xx4,yy4 ',i,j,xx4(i,j),yy4(i,j)
        endif
       enddo  ! loop loop
c      expect xg, yg to range between .5 and il+.5
       xg(iq,k)=.25*(ri+3.) -.5  ! -.5 for stag; back to normal ri, rj defn
       yg(iq,k)=.25*(rj+3.) -.5  ! -.5 for stag
      enddo   ! iq loop
      call end_log(toij_end)
      return
      end
      subroutine checkdiv(xstr,ystr,zstr)
!     Check whether optimisation uses multiplication by reciprocal so
!     that x/x /= 1.
      implicit none
      include 'newmpar.h'
      real xstr(ifull),ystr(ifull),zstr(ifull)
      real denxyz
      integer iq
      integer, parameter :: n=100

      call random_number(xstr(1:n))
      ystr(1:n) = 0.9*xstr(1:n)
      zstr(1:n) = 0.8*xstr(1:n)
      ! By construction here, xstr is largest, so xstr(iq)/denxyz should be
      ! 1.
      do iq=1,n
        denxyz=max( abs(xstr(iq)),abs(ystr(iq)),abs(zstr(iq)) )
        xstr(iq) = xstr(iq)/denxyz
        ystr(iq) = ystr(iq)/denxyz
        zstr(iq) = zstr(iq)/denxyz
      end do
      if ( any(xstr(1:n)/=1.0) ) then
         print*, "Error, must use ncray=1 on this machine"
         stop
      end if
      end subroutine checkdiv
