c     numtop > 0 for topography overlay (chooses contour in metres also)
c            = 0 for no overlay
c            = -1 for land-sea overlay
c            = -3 for lat-long supmap
c            =-10 or less for panel 0 and so on
c     nplot=1 does plots of grid structure  
c             (nplot & mplot still being debugged)
c     mplot=0,1 for plots on circular domain (1 with supmap)
c          =11 circular domain of grid structure (as in Helmholtz plot) -not working
c          =2,3 for plots on lat/long domain (3 with supmap)
c          =-2,-4 for plots on circular domain (every 2nd or 4th point)
c          =-102,-104 for alt. plots on circular domain 
c                             (every 2nd or 4th point)
c               e.g. with nplot=-3
!     maplon (0 off), maplat, radmap give viewing window with mplot=3
c     rlong0,rlat0,schmidt can be read in namelist
      include 'jimcc.f'
      include 'setxyz.f'
      program plotg    ! for globpe
c     all contour plotting done via [1:360]x[-89:89] array
      include 'newmpar.h'
      parameter (nmax=4001)
      include 'latlong.h'  ! rlatt,rlongg
      include 'parm.h'
      include 'xyzinfo.h'  ! x,y,z,wts
      common/lunit/luw
      real vecx(3),vecy(3),vecz(3)
      character ofile*40,d*20
      common/xynf/nface(1:360,-89:89),xg(1:360,-89:89),yg(1:360,-89:89)

      namelist/plotnml/ofile,ndeg
     . ,mplot,nplot,rlong0,rlat0,schmidt,viewrad,d
     . ,maplon,maplat,radmap
     . ,id,jd,vecx,vecy
      data ndeg/20/,mplot/3/,nplot/0/,viewrad/0./,maplon/0/
      data vecx/0.,0.,0./
      print *,'Plot routine compiled for il,jl,kl = ',il,jl,kl
      luw=30
      rlong0=0.
      rlat0=90.
      schmidt=1.
      read (5, plotnml)
      write(6, plotnml)
      open(luw,file=ofile,status='unknown')
      call setxyz
      call pltgridx(nplot,mplot,ndeg,viewrad,vecx,vecy,vecz,
     .              maplon,maplat,radmap)
      end

      subroutine pltgridx(nplot,mplot,ndeg,viewrad,vecx,vecy,vecz,
     .                    maplon,maplat,radmap)
      parameter (ntest=0)  ! 0 for off
      include 'newmpar.h'
      include 'bigxy4.h' 
      include 'indices.h' ! in,is,iw,ie,inn,iss,iww,iee
      include 'latlong.h'  ! rlatt,rlongg
      include 'parm.h'
      common/lunit/luw
      include 'xyzinfo.h'  ! x,y,z,wts
      common/work2/rlong(ifull),rlat(ifull),dum(ifull,16)
      real xa(0:7),ya(0:7),za(0:7),xb(0:7),yb(0:7),zb(0:7)
      real rlata(0:7),rlonga(0:7),rlatb(0:7),rlongb(0:7)
      real vecx(3),vecy(3),vecz(3)
      real veca(3),vecb(3)
      real rotpole(3,3)
      data pi/3.1415926536/
      print *,'in pltgridx rlong0,rlat0,schmidt ',rlong0,rlat0,schmidt
      viewradv=0.
      if(viewrad.gt.0.)viewradv=1./viewrad  ! viewrad=0 is infinity
      if(nplot.eq.-2)then    ! unrotated view from above face
        vecx(1)=1.
        vecx(2)=0.
        vecx(3)=0.
        vecy(1)=0.
        vecy(2)=1.
        vecy(3)=0.
        vecz(1)=0.
        vecz(2)=0.
        vecz(3)=1.
      endif
      if(nplot.eq.-3)then    ! view from above vertex
        vecx(1)=1./sqrt(3.)
        vecx(2)=1./sqrt(3.)
        vecx(3)=1./sqrt(3.)
        vecy(1)=-2./sqrt(6.)
        vecy(2)=1./sqrt(6.)
        vecy(3)=1./sqrt(6.)
        vecz(1)=0.
        vecz(2)=-1./sqrt(2.)
        vecz(3)=1./sqrt(2.)
      endif
      if(nplot.eq.-4)then    ! view from above edge
        vecx(1)=1./sqrt(2.)
        vecx(2)=0.
        vecx(3)=1./sqrt(2.)
        vecy(1)=0.
        vecy(2)=1.
        vecy(3)=0.
        vecz(1)=-1./sqrt(2.)
        vecz(2)=0.
        vecz(3)=1./sqrt(2.)
      endif
      if(nplot.eq.-5)then    ! read in rotation vectors
        open(65,file='3dvecs',status='old')
        read (65,*) vecx(1),vecy(1),vecz(1)  ! this one lat/longs
        print *,'lat/longs ',vecx(1),vecy(1),vecz(1)
        read (65,*) vecx(1),vecy(1),vecz(1)
        read (65,*) vecx(2),vecy(2),vecz(2)
        read (65,*) vecx(3),vecy(3),vecz(3)
      endif
      print *,'in plgridx '
      print *,'vecx ',vecx
      print *,'vecy ',vecy
      if(nplot.eq.2)then
!       uses vecx and vecy to rotate view of grid      
!       vecx is thru panel 0, vecy thru panel 2, vecz thru panel 1
!       as default for rlat0=rlong0=0
!       normalize vecx
        call unit(vecx)
!       define vecz=vecx X vecy
        call cross1(vecz,vecx,vecy)
        call unit(vecz)
!       assume vecy is only approximate vector, so correct it
!       recalculate vecy  = vecz X vecx
        call cross1(vecy,vecz,vecx)
      endif  ! nplot.eq.2
      if(nplot.eq.3)then
!       take vecx, vecy as being coords of desired diagonal upper vertices
!       find these using fort.22, then run ccam with approriate id,jd 
!       to determine correpsonding values of x,y,z
        vecz(1)=.5*(vecx(1)+vecy(1))
        vecz(2)=.5*(vecx(2)+vecy(2))
        vecz(3)=.5*(vecx(3)+vecy(3))
        call unit(vecz)
!       let veca be vector joining vertices, vecb = veca X vecz
        veca(1)=vecx(1)-vecy(1)
        veca(2)=vecx(2)-vecy(2)
        veca(3)=vecx(3)-vecy(3)
        call unit(veca)
        call cross1(vecb,veca,vecz)
        vecx(1)=.5*(veca(1)+vecb(1))
        vecx(2)=.5*(veca(2)+vecb(2))
        vecx(3)=.5*(veca(3)+vecb(3))
        call unit(vecx)
        call cross1(vecy,vecz,vecx)
      endif  ! nplot.eq.3

      print *,'vecx ',vecx
      print *,'vecy ',vecy
      print *,'vecz ',vecz
      print *,'x,y,z ',x(idjd),y(idjd),z(idjd)
      print *,'rlongg,rlatt ',rlongg(idjd),rlatt(idjd)
      print *,'rlong,rlat in degrees ',
     .         rlongg(idjd)*180./pi,rlatt(idjd)*180./pi
      if(nplot.eq.2.or.nplot.eq.3)then
!       view plots rotated by vecx, vecy, vecz; x goes to vecx etc
        do iq=1,ifull
         xx=x(iq)
         yy=y(iq)
         zz=z(iq)
c         x(iq)=vecx(1)*xx+vecx(2)*yy+vecx(3)*zz
c         y(iq)=vecy(1)*xx+vecy(2)*yy+vecy(3)*zz
c         z(iq)=vecz(1)*xx+vecz(2)*yy+vecz(3)*zz
         x(iq)=vecx(1)*xx+vecy(1)*yy+vecz(1)*zz
         y(iq)=vecx(2)*xx+vecy(2)*yy+vecz(2)*zz
         z(iq)=vecx(3)*xx+vecy(3)*yy+vecz(3)*zz
        enddo   ! iq loop
        print *,'new x,y,z ',x(idjd),y(idjd),z(idjd)
      endif     ! (nplot.eq.2.or.nplot.eq.3)
!     rotpole(1,) is x-axis of rotated coords in terms of orig Cartesian
!     rotpole(2,) is y-axis of rotated coords in terms of orig Cartesian
!     rotpole(3,) is z-axis of rotated coords in terms of orig Cartesian
      coslong=cos(rlong0*pi/180.)
      sinlong=sin(rlong0*pi/180.)
      coslat=cos(rlat0*pi/180.)
      sinlat=sin(rlat0*pi/180.)
      rotpole(1,1)=coslong*sinlat
      rotpole(1,2)=-sinlong
      rotpole(1,3)=coslong*coslat
      rotpole(2,1)=sinlong*sinlat
      rotpole(2,2)=coslong
      rotpole(2,3)=sinlong*coslat
      rotpole(3,1)=-coslat
      rotpole(3,2)=0.
      rotpole(3,3)=sinlat
      do iq=1,ifull
!      x(), y(z), z() are "local" coords with z out of central panel
!      while xx, yy, zz are "true" Cartesian values
!      xx is new x after rot by rlong0 then rlat0
       xx=rotpole(1,1)*x(iq)+rotpole(1,2)*y(iq)+rotpole(1,3)*z(iq)
       yy=rotpole(2,1)*x(iq)+rotpole(2,2)*y(iq)+rotpole(2,3)*z(iq)
       zz=rotpole(3,1)*x(iq)+rotpole(3,2)*y(iq)+rotpole(3,3)*z(iq)
c      xx=x(iq)
c      yy=y(iq)
c      zz=z(iq)
       rlatt(iq)=asin(zz)
       if(yy.ne.0..or.xx.ne.0.)then
         rlongg(iq)=atan2(yy,xx)                       ! N.B. -pi to pi
         if(rlongg(iq).lt.0.)rlongg(iq)=rlongg(iq)+2.*pi ! 0 to 2*pi  
       else
         rlongg(iq)=0.    ! a default value for NP/SP
       endif
      enddo   ! iq loop
      print *,'re-calculated rlongg, rlatt ',
     .         rlongg(idjd),rlatt(idjd)
      print *,'re-calculated rlong, rlat in degrees',
     .         rlongg(idjd)*180./pi,rlatt(idjd)*180./pi

      nnup=5
      if(npanels.eq.13)nnup=7
      if(mplot.eq.0.or.mplot.eq.1)then   ! plots on circular domain
c       just plot x,y links for +ve z
        if(viewradv.gt.0.)then
          call mapstr('SA',1./viewradv)
          if(mplot.eq.1)call supmap(7,rlat0,rlong0,
     .                             0.,-90.,90.,-90.,90.,1,ndeg,0,0,ierr)
        else
          if(mplot.eq.1)call supmap(2,rlat0,rlong0,
     .                             0.,-90.,90.,-90.,90.,1,ndeg,0,0,ierr)
        endif    !  (viewradv.gt.0.)
        ztest=viewradv    !  i.e. 1/viewing radius, default 0.
c       facta scales x so as to have maximum of 1. (i.e. for z=ztest)
        facta=sqrt(1.-ztest**2)            ! default 1.
        do iq=1,ifull
         if(z(iq).ge.ztest.and.z(ie(iq)).ge.ztest)then
           fact=facta/(1.-ztest*z(iq))     ! default 1.
           call frstpt(x(iq)*fact,y(iq)*fact)
           fact=facta/(1.-ztest*z(ie(iq)))
           call vector(x(ie(iq))*fact,y(ie(iq))*fact)
         endif
         if(z(iq).ge.ztest.and.z(in(iq)).ge.ztest)then
           fact=facta/(1.-ztest*z(iq))
           call frstpt(x(iq)*fact,y(iq)*fact)
           fact=facta/(1.-ztest*z(in(iq)))
           call vector(x(in(iq))*fact,y(in(iq))*fact)
         endif
        enddo  ! iq loop

      elseif(mplot.eq.12)then
!       in this section, xa coord is out of page (like z above)
        ztest=viewradv    !  i.e. 1/viewing radius, default 0.
c       facta scales x so as to have maximum of 1. (i.e. for z=ztest)
        facta=sqrt(1.-ztest**2)            ! default 1.
        print *,'in pltgridx with mplot=11'
        if(viewradv.gt.0.)then
          call mapstr('SA',1./viewradv)
          call supmap(7,rlat0,rlong0,
     .                             0.,-90.,90.,-90.,90.,1,ndeg,0,0,ierr)
        else
          call supmap(2,rlat0,rlong0,
     .                             0.,-90.,90.,-90.,90.,1,ndeg,0,0,ierr)
        endif    !  (viewradv.gt.0.)

!       call set(.05,.95,.05,.95,-1.,1.,-1.,1.,1)  
!       draw outline of grid elements on circular domain
!       work through xx4, yy4 values
        do j=1,iquad,4
         xx=xx4(1,j)
         yy=yy4(1,j)
         call xyz8u(xa,ya,za,xx,yy)
         do i=2,iquad
          xx=xx4(i,j)
          yy=yy4(i,j)
          call xyz8u(xb,yb,zb,xx,yy)
          print *,'i,j,xx,yy: ',i,j,xx,yy
          do nn=0,nnup
           print *,'nn  xb,yb,zb: ',nn,xb(nn),yb(nn),zb(nn)
           if(za(nn).ge.ztest.and.zb(nn).ge.ztest)then
             fact=facta/(1.-ztest*za(nn))
             call frstpt(xa(nn)*fact,ya(nn)*fact)
             fact=facta/(1.-ztest*zb(nn))
             call vector(xb(nn)*fact,yb(nn)*fact)
           if(ntest.ne.0)print *,'join A ',xa(nn),ya(nn),xb(nn),yb(nn)
           endif
           xa(nn)=xb(nn)
           ya(nn)=yb(nn)
           za(nn)=zb(nn)
          enddo  ! nn loop
         enddo   ! i loop
        enddo    ! j loop
        do i=1,iquad,4
         xx=xx4(i,1)
         yy=yy4(i,1)
         call xyz8u(xa,ya,za,xx,yy)
         do j=2,iquad
          xx=xx4(i,j)
          yy=yy4(i,j)
          call xyz8u(xb,yb,zb,xx,yy)
          do nn=0,nnup
           if(za(nn).ge.ztest.and.zb(nn).ge.ztest)then
             fact=facta/(1.-ztest*za(nn))
             call frstpt(xa(nn)*fact,ya(nn)*fact)
             fact=facta/(1.-ztest*zb(nn))
             call vector(xb(nn)*fact,yb(nn)*fact)
           if(ntest.ne.0)print *,'join B ',xa(nn),ya(nn),xb(nn),yb(nn)
           endif
           xa(nn)=xb(nn)
           ya(nn)=yb(nn)
           za(nn)=zb(nn)
          enddo  ! nn loop
         enddo   ! i loop
        enddo    ! j loop
!c       draw a unit circle around the edge
!        call frstpt(1.,0.)
!        do ii=1,360
!         call vector(cos(2*ii*pi/360.),sin(2*ii*pi/360.))
!        enddo

      elseif(mplot.eq.4)then
!       draw outline of grid elements on lat-long map
        call supmap(8,0.,180.,0.,0.,360.,-90.,90.,1,ndeg,0,0,ierr)
c       call set(.05,.95,.275,.725,0.,360.,-90.,90.,1)  ! old supmap
        call set(0. ,1. ,.25 ,.75 ,0.,360.,-90.,90.,1)  ! new ncar
        do j=1,iquad,4
         xx=xx4(1,j)
         yy=yy4(1,j)
         call xyz8(xa,ya,za,xx,yy)
         do nn=0,nnup
          rlata(nn)=asin(za(nn))*180./pi
          if(ya(nn).ne.0..or.xa(nn).ne.0.)then
            rlonga(nn)=atan2(ya(nn),xa(nn))*180./pi
            if(rlonga(nn).lt.0.)rlonga(nn)=rlonga(nn)+360.
          else
            rlonga(nn)=0.
          endif
         enddo  ! nn loop
         do i=2,iquad
          xx=xx4(i,j)
          yy=yy4(i,j)
          call xyz8(xb,yb,zb,xx,yy)
          do nn=0,nnup
           rlatb(nn)=asin(zb(nn))*180./pi
           if(yb(nn).ne.0..or.xb(nn).ne.0.)then
             rlongb(nn)=atan2(yb(nn),xb(nn))*180./pi
             if(rlongb(nn).lt.0.)rlongb(nn)=rlongb(nn)+360.
           else
             rlongb(nn)=0.
           endif
           write(21,'(4f8.3)')rlonga(nn),rlata(nn),rlongb(nn),rlatb(nn)
           if(abs(rlongb(nn)-rlonga(nn)).lt.40..and.
     .        abs(rlatb(nn)-rlata(nn)).lt.40.)then
              call frstpt(rlonga(nn),rlata(nn))
              call vector(rlongb(nn),rlatb(nn))
           endif  ! (abs(rlongb(nn)-rlonga(nn)).lt.40...
           rlata(nn)=rlatb(nn)
           rlonga(nn)=rlongb(nn)
          enddo  ! nn loop
         enddo  ! i loop
        enddo   ! j loop

        do i=1,iquad,4
         xx=xx4(i,1)
         yy=yy4(i,1)
         call xyz8(xa,ya,za,xx,yy)
         do nn=0,nnup
          rlata(nn)=asin(za(nn))*180./pi
          if(ya(nn).ne.0..or.xa(nn).ne.0.)then
            rlonga(nn)=atan2(ya(nn),xa(nn))*180./pi
            if(rlonga(nn).lt.0.)rlonga(nn)=rlonga(nn)+360.
          else
            rlonga(nn)=0.
          endif
         enddo  ! nn loop
         do j=2,iquad
          xx=xx4(i,j)
          yy=yy4(i,j)
          call xyz8(xb,yb,zb,xx,yy)
          do nn=0,nnup
           rlatb(nn)=asin(zb(nn))*180./pi
           if(yb(nn).ne.0..or.xb(nn).ne.0.)then
             rlongb(nn)=atan2(yb(nn),xb(nn))*180./pi
             if(rlongb(nn).lt.0.)rlongb(nn)=rlongb(nn)+360.
           else
             rlongb(nn)=0.
           endif
           write(21,'(4f7.2)')rlonga(nn),rlata(nn),rlongb(nn),rlatb(nn)
           if(abs(rlongb(nn)-rlonga(nn)).lt.40..and.
     .        abs(rlatb(nn)-rlata(nn)).lt.40.)then
              call frstpt(rlonga(nn),rlata(nn))
              call vector(rlongb(nn),rlatb(nn))
           endif  ! (abs(rlongb(nn)-rlonga(nn)).lt.40...
           rlata(nn)=rlatb(nn)
           rlonga(nn)=rlongb(nn)
          enddo  ! nn loop
         enddo  ! j loop
        enddo   ! i loop

      else   ! plots on lat/long domain   mplot=2 or 3
        if(maplon.eq.0)then
	   rlon1=0.
	   rlon2=360.
	   rlat1=-90.
	   rlat2=90.
	   yset1=.25
	   yset2=.75
	 else
	   rlon1=maplon-radmap
	   rlon2=maplon+radmap
	   rlat1=maplat-radmap
	   rlat2=maplat+radmap
	   yset1=0.
	   yset2=1.
	   print *,'non-zero maplon: ',maplon
	   print *,'rlon1,rlon2,rlat1,rlat2 ',rlon1,rlon2,rlat1,rlat2
	 endif
        if(mplot.eq.3)call supmap
     .            (8,0.,180.,0.,rlon1,rlon2,rlat1,rlat2,1,ndeg,0,0,ierr)
        call set(0.,1.,yset1 ,yset2,rlon1,rlon2,rlat1,rlat2,1)  ! new ncar
c       just plot y,z links for +ve x
c       calculate lat/longs on strips  longs 0 to 360
        do iq=1,ifull
            rlat(iq)=rlatt(iq)*180./pi
            rlong(iq)=rlongg(iq)*180./pi
            if(rlong(iq).lt.0.)rlong(iq)=rlong(iq)+360.
        enddo   ! iq loop
        do iq=1,ifull
            if(abs(rlong(ie(iq))-rlong(iq)).lt.40..and.
     .         abs(rlat(ie(iq))-rlat(iq)).lt.40.)then
              call frstpt(rlong(iq),rlat(iq))
              call vector(rlong(ie(iq)),rlat(ie(iq)))
            endif
            if(abs(rlong(in(iq))-rlong(iq)).lt.40..and.
     .         abs(rlat(in(iq))-rlat(iq)).lt.40.)then
              call frstpt(rlong(iq),rlat(iq))
              call vector(rlong(in(iq)),rlat(in(iq)))
            endif
        enddo   ! iq loop
      endif
      call frame
      if(nplot.ne.0)stop 'finished in pltgridx'
      return
      end
      subroutine xyz8(xa,ya,za,xx,yy)
!     for unrotated version, see xyz8u
      include 'newmpar.h'
      include 'parm.h'
      real xa(0:7),ya(0:7),za(0:7)
      real rotpole(3,3)
      save num
      data pi/3.1415926536/,num/0/
c     calculates the 6 (or 8) rotated and stretched values corresponding
c     to input (xx,yy) from the reference panel
c     assumes the order is irrelevant
      if(npanels.eq.5)then
        nnup=5
        xa(0)=1.
        ya(0)=xx
        za(0)=yy
        xa(3)=-1.
        za(3)=-xx
        ya(3)=-yy
        xa(1)=-yy
        ya(1)=xx
        za(1)= 1.
        ya(4)=-yy
        xa(4)=xx
        za(4)=-1.
        xa(2)=-yy
        ya(2)= 1.
        za(2)=-xx
        za(5)=yy
        ya(5)=-1.
        xa(5)=xx
        do n=0,npanels
         call norm(xa(n),ya(n),za(n),den)
!        xa, ya, za are now coords on sphere  -1 to 1
        enddo    ! n loop
      elseif(npanels.eq.13)then
        nnup=7
c       unstretch xx,yy,zz (assuming they were created with schm13=.1)
!       schm13=.1  ! originally used .5 ! now comes through parm.h
        alf=(1.-schm13**2)/(1.+schm13**2)
        alfonsch=(1.-alf)/schm13
        zz=sqrt(1.-xx**2-yy**2)
        x=xx*alfonsch/(1.-alf*zz)
        y=yy*alfonsch/(1.-alf*zz)
        z=(zz-alf)/(1.-alf*zz)
c       calculate xxn,yyn,zzn for upper hemisphere
c       and       xxs,yys,zzs for lower hemisphere assuming schmidt=.1
        xa(0)= x
        ya(0)= y
        za(0)= z
        xa(1)= x
        ya(1)=-y
        za(1)= z
        xa(2)=-x
        ya(2)=-y
        za(2)= z
        xa(3)=-x
        ya(3)= y
        za(3)= z
        xa(4)= x
        ya(4)= y
        za(4)=-z
        xa(5)= x
        ya(5)=-y
        za(5)=-z
        xa(6)=-x
        ya(6)=-y
        za(6)=-z
        xa(7)=-x
        ya(7)= y
        za(7)=-z
      endif      ! (npanels.eq.5) elseif(npanels.eq.13)
c     now perform schmidt etc
      coslong=cos(rlong0*pi/180.)
      sinlong=sin(rlong0*pi/180.)
      coslat=cos(rlat0*pi/180.)
      sinlat=sin(rlat0*pi/180.)
      rotpole(1,1)=coslong*sinlat
      rotpole(1,2)=-sinlong
      rotpole(1,3)=coslong*coslat
      rotpole(2,1)=sinlong*sinlat
      rotpole(2,2)=coslong
      rotpole(2,3)=sinlong*coslat
      rotpole(3,1)=-coslat
      rotpole(3,2)=0.
      rotpole(3,3)=sinlat
      alf=(1.-schmidt**2)/(1.+schmidt**2)
      num=num+1
      do nn=0,nnup
c      if(num.lt.1000)print *,'a x,y,z: ',nn,xa(nn),ya(nn),za(nn)
       zin=za(nn)
       xx=xa(nn)*schmidt*(1.+alf)/(1.+alf*zin)
       yy=ya(nn)*schmidt*(1.+alf)/(1.+alf*zin)
       zz=(alf+zin)/(1.+alf*zin)
       xa(nn)=rotpole(1,1)*xx+rotpole(1,2)*yy+rotpole(1,3)*zz
       ya(nn)=rotpole(2,1)*xx+rotpole(2,2)*yy+rotpole(2,3)*zz
       za(nn)=rotpole(3,1)*xx+rotpole(3,2)*yy+rotpole(3,3)*zz
c      if(num.lt.1000)print *,'b x,y,z: ',nn,xa(nn),ya(nn),za(nn)
      enddo  !  nn loop
      return
      end

      subroutine xyz8u(xa,ya,za,xx,yy)
!     this is unrotated version
      include 'newmpar.h'
      include 'parm.h'
      real xa(0:7),ya(0:7),za(0:7)
      real rotpole(3,3)
      save num
      data pi/3.1415926536/,num/0/
c     calculates the 6 values corresponding to
c     input (xx,yy) from the reference panel
c     assumes the order is irrelevant
        nnup=5
        xa(0)=1.
        ya(0)=xx
        za(0)=yy
        xa(3)=-1.
        za(3)=-xx
        ya(3)=-yy
        xa(1)=-yy
        ya(1)=xx
        za(1)= 1.
        ya(4)=-yy
        xa(4)=xx
        za(4)=-1.
        xa(2)=-yy
        ya(2)= 1.
        za(2)=-xx
        za(5)=yy
        ya(5)=-1.
        xa(5)=xx
        do n=0,npanels
         call norm(xa(n),ya(n),za(n),den)
!        xa, ya, za are now coords on sphere  -1 to 1
        enddo    ! n loop
      alf=(1.-schmidt**2)/(1.+schmidt**2)
      num=num+1
      do nn=0,nnup
c      if(num.lt.1000)print *,'a x,y,z: ',nn,xa(nn),ya(nn),za(nn)
       zin=za(nn)
       xa(nn)=xa(nn)*schmidt*(1.+alf)/(1.+alf*zin)
       ya(nn)=ya(nn)*schmidt*(1.+alf)/(1.+alf*zin)
       za(nn)=(alf+zin)/(1.+alf*zin)
      enddo  !  nn loop
      return
      end
      subroutine cross1(c,a,b)
c     calculate vector components of c = a x b
c     where each RHS component represents 3 vector components
c     and convert c to a unit vector
      real c(3),a(3),b(3)
      c(1)=a(2)*b(3)-b(2)*a(3)
      c(2)=a(3)*b(1)-b(3)*a(1)
      c(3)=a(1)*b(2)-b(1)*a(2)
      return
      end
      subroutine unit(c)
c     convert vector c to unit vector
      real c(3)
      cmag=sqrt(c(1)**2+c(2)**2+c(3)**2)
      c(1)=c(1)/cmag
      c(2)=c(2)/cmag
      c(3)=c(3)/cmag
      return
      end

      subroutine set(a,b,c,d,e,f,g,h,i)
      common/lunit/luw
      write(luw,91)
91    format('SET   ')
      write(luw,92)a,b,c,d,e,f,g,h,i
92    format(4f8.5,1p,4e10.3,i4)
      return
      end
      subroutine frstpt(x,y)
      common/lunit/luw
      write(luw,91)
91    FORMAT('FRSTPT')
      write(luw,92)x,y
92    format(2e13.5)
      return
      end
      subroutine vector(x,y)
      common/lunit/luw
      write(luw,91)
91    FORMAT('VECTOR')
      write(luw,92)x,y
92    format(2e13.5)
      return
      end
      subroutine packit(a,idim,iout,jout)
      include 'newmpar.h'
      common/lunit/luw
      dimension a(idim,jout),ia(360*179)
c     dimension ib(2500),ic(2500)
c  scale array and write out as integers between 0 and 4095 ( = fff hex)
c     for decimal version replace 4095 by 9999, and 26z3 by 20i4
      iijj=iout*jout
      amin=a(1,1)
      amax=amin+1.e-10
      do j=1,jout
       do i=1,iout
        amin=min(a(i,j),amin)
        amax=max(a(i,j), amax)
       enddo
      enddo
      write(luw,91)amin,amax
91    format(1p,2e13.5)
      n=0
      do j=1,jout
       do i=1,iout
        n=n+1
        ia(n)=(a(i,j)-amin)*4095/(amax-amin)
       enddo
      enddo
      write(luw,92)(ia(i),i=1,iijj)
92    format(26z3)
c      this is a base 64 version
c      ib(1;iijj)=ia(1;iijj)/64
c      ic(1;iijj)=ia(1;iijj)-64*ib(1;iijj)+32
c      ib(1;iijj)=ib(1;iijj)+32
c      write(luw,92) (ib(i),ic(i),i=1,iijj)
c92    format(80r1)
      return
      end
      subroutine supmap(i1,a1,a2,a3,a4,a5,a6,a7,i2,i3,i4,i5,ierr)
      common/lunit/luw
      write(luw,91)
91    FORMAT('SUPMAP')
      write(luw,92) i1,a1,a2,a3,a4,a5,a6,a7,i2,i3,i4,i5
92    format(i3,7f7.2,4i3)
      return
      end
      subroutine mapstr(i1,a1)
      common/lunit/luw
      character i1*2
      write(luw,91)
91    FORMAT('MAPSTR')
      write(luw,92) a1
92    format(f8.4)
      return
      end
      subroutine frame
      common/lunit/luw
      data numframe/0/
      numframe=numframe+1
      write(luw,91) numframe
91    FORMAT('FRAME ',i3)
      end
