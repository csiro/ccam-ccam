      subroutine setxyz(ik,rlong0,rlat0,schmidtin,            ! input
     &      x,y,z,wts, ax,ay,az,bx,by,bz,xx4,yy4)             ! output
!     this routine modified sept '06 to accept smaller il_g via abs(ik)
!     essentially to temporarily provide xx4, yy4, ax...bz for onthefly  
!     note that ax6 etc not needed for onthefly     
      use cc_mpi, only : indx
      use indices_m
      use jimcc_m
      use latlong_m
      use map_m
      use utilities
      use workglob_m
      implicit none
c     integer, parameter :: ntang=2  ! always done
                            ! ntang=0 for tang. vectors from Rancic et al.
                            !         not for stretched as vecpanel not ready
                            ! ntang=1 for tang. vectors by finite diffs
                            ! ntang=2 for map factors by finite diffs too
c     schmidt included
c     sets up x, y, z on sphere and unit u,v vectors
c     note that x,y,z have been normalized by rearth, the radius of the earth
c     suffix 6 denotes hex (6)
      include 'newmpar.h'
      include 'const_phys.h'   ! rearth
      include 'parm.h'
      integer, intent(in) :: ik  ! passed as argument. Actual i dimension.
!                                  if negative, suppress calc of rlat4, rlong4, indices,em_g                             
      integer i,j,n,ikk,idjd_g,iq,ii,i0,j0,n0,in0,jn0,nn0
      integer is0,js0,ns0,ie0,je0,ne0,iw0,jw0,nw0,inn0,jnn0
      integer iss0,jss0,nss0,nnn0,iee0,jee0,nee0,iww0,jww0,nww0,m
      integer iq11,iq12,iq13,iq22,iq32,iqcc,iqnn
      integer imin,imax,jmin,jmax,numpts
      integer iqm,iqp, n_n, n_e, n_w, n_s
      integer iquadx
      integer, save :: num = 0
      integer, parameter :: ndiag = 0
      real(kind=8) x(ik*ik*6),y(ik*ik*6),z(ik*ik*6)
      real ax(ik*ik*6),ay(ik*ik*6),az(ik*ik*6),wts(ik*ik*6)
      real bx(ik*ik*6),by(ik*ik*6),bz(ik*ik*6)
      real(kind=8) xx4(1+4*ik,1+4*ik),yy4(1+4*ik,1+4*ik)
      real em4(1+4*ik,1+4*ik)
     .    ,ax4(1+4*ik,1+4*ik),ay4(1+4*ik,1+4*ik)
     &    ,az4(1+4*ik,1+4*ik)
     .    ,axx(ik*ik*6),ayy(ik*ik*6),azz(ik*ik*6)
     .    ,bxx(ik*ik*6),byy(ik*ik*6),bzz(ik*ik*6)
      real rlong0,rlat0,schmidt,schmidtin
      real rotpole(3,3)
      real(kind=8) alf,den1,xx,yy,zz,x4_iq_m,y4_iq_m,z4_iq_m
      real(kind=8), parameter :: one = 1._8
      real dsfact,xin,yin,zin
      real den, dot,eps,dx2,dy2,sumwts,rlat,rlong,ratmin,ratmax,rat
      real rlatdeg,rlondeg
      num=num+1
c     When using the ifull_g notation: in_g, ie_g, iw_g and is_g give the
c     indices for the n, e, w, s neighbours respectively
c     a, b denote unit vectors in the direction of x, y (e & n) respectively
      idjd_g = id+il_g*(jd-1)  ! Global value
      schmidt=abs(schmidtin)
      ikk=abs(ik)
      iquadx=1+ik*((8*npanels)/(npanels+4))
      
      em4=0. ! for cray compiler

      ! MJT notes - indices are now defined in indices_m.f90 as
      ! functions to save memory with global arrays
      
      if(ndiag.eq.3)then
        do n=0,npanels
         do j=1,ikk
          do i=1,ikk
           iq=indx(i,j,n,ikk,ikk)
           call indv(iq,i0,j0,n0)
           call indv(in_g(iq),in0,jn0,nn0)
           call indv(is_g(iq),is0,js0,ns0)
           call indv(ie_g(iq),ie0,je0,ne0)
           call indv(iw_g(iq),iw0,jw0,nw0)
           call indv(inn_g(iq),inn0,jnn0,nnn0)
           call indv(iss_g(iq),iss0,jss0,nss0)
           call indv(iee_g(iq),iee0,jee0,nee0)
           call indv(iww_g(iq),iww0,jww0,nww0)
           print 91,i0,j0,n0,
     .      in0,jn0,nn0,is0,js0,ns0,ie0,je0,ne0,iw0,jw0,nw0,
     .      inn0,jnn0,nnn0,iss0,jss0,nss0,iee0,jee0,nee0,iww0,jww0,nww0
91         format(9(i4,i2,i2))
          enddo   ! i loop
         enddo   ! j loop
        enddo    ! n loop
      endif      ! (ndiag.eq.3)

!----------------------------------------------------------------------------
c     calculate grid information using quadruple resolution grid
      call jimcc(em4,ax4,ay4,az4,xx4,yy4,ikk)
c     call jimcc(em4,ax4,ay4,az4,myid)
      if(ktau<=1)then
        print *,'xx4 first & last ',xx4(1,1),xx4(iquadx,iquadx)
        print *,'xx4 (5,5),(7,7),(9,9) ',xx4(5,5),xx4(7,7),xx4(9,9)
        print *,'yy4 first & last ',yy4(1,1),yy4(iquadx,iquadx)
        print *,'yy4 (5,5),(7,7),(9,9) ',yy4(5,5),yy4(7,7),yy4(9,9)
        print *,'xx4, yy4 central',xx4(2*ikk+1,2*ikk+1),
     &                               yy4(2*ikk+1,2*ikk+1)
      endif  ! (ktau<=1)

!     rotpole(1,) is x-axis of rotated coords in terms of orig Cartesian
!     rotpole(2,) is y-axis of rotated coords in terms of orig Cartesian
!     rotpole(3,) is z-axis of rotated coords in terms of orig Cartesian
      rotpole = calc_rotpole(rlong0,rlat0)

      if(schmidtin<0.)go to 3
!     following just for rlong4, rlat4 (& x4,y4,z4)
      alf=(one-schmidt**2)/(one+schmidt**2)
      do m=1,4
       do j=1,ikk
         do i=1,ikk
          if(m.eq.1)then
            xx=xx4(4*i-1-1,4*j-1-1)
            yy=yy4(4*i-1-1,4*j-1-1)
          endif
          if(m.eq.2)then
            xx=xx4(4*i-1-1,4*j-1+1)
            yy=yy4(4*i-1-1,4*j-1+1)
          endif
          if(m.eq.3)then
            xx=xx4(4*i-1+1,4*j-1-1)
            yy=yy4(4*i-1+1,4*j-1-1)
          endif
          if(m.eq.4)then
            xx=xx4(4*i-1+1,4*j-1+1)
            yy=yy4(4*i-1+1,4*j-1+1)
          endif
c         set up x0, y0, z0 coords on cube -1 to 1
c         avoids earlier equivalencing of x,x0  etc
          x(indx(i,j,0,ikk,ikk))= 1.
          y(indx(i,j,0,ikk,ikk))=xx
          z(indx(i,j,0,ikk,ikk))=yy
          x(indx(i,j,3,ikk,ikk))=-1.
          z(indx(i,j,3,ikk,ikk))=-xx
          y(indx(i,j,3,ikk,ikk))=-yy
          x(indx(i,j,1,ikk,ikk))=-yy
          y(indx(i,j,1,ikk,ikk))=xx
          z(indx(i,j,1,ikk,ikk))= 1.
          y(indx(i,j,4,ikk,ikk))=-yy
          x(indx(i,j,4,ikk,ikk))=xx
          z(indx(i,j,4,ikk,ikk))=-1.
          x(indx(i,j,2,ikk,ikk))=-yy
          y(indx(i,j,2,ikk,ikk))= 1.
          z(indx(i,j,2,ikk,ikk))=-xx
          z(indx(i,j,5,ikk,ikk))=yy
          y(indx(i,j,5,ikk,ikk))=-1.
          x(indx(i,j,5,ikk,ikk))=xx
c         if(i.eq.(ikk+1)/2.and.j.eq.(ikk+1)/2)print *,'n,xx,yy: ',n,xx,yy
         enddo  ! i loop
       enddo   ! j loop
       do iq=1,6*ik*ik
        call norm8(x(iq),y(iq),z(iq),den1) ! x, y, z are coords on sphere  -1 to 1
        x4_iq_m=x(iq)*schmidt*(1.+alf)/(1.+alf*z(iq))
        y4_iq_m=y(iq)*schmidt*(1.+alf)/(1.+alf*z(iq))
        z4_iq_m=(alf+z(iq))/(1.+alf*z(iq)) 
!       here is calculation of rlong4, rlat4
c       also provide latitudes and longitudes (-pi to pi)
        if(rlong0.eq.0..and.rlat0.eq.90.)then
          xx=x4_iq_m
          yy=y4_iq_m
          zz=z4_iq_m
        else
!         x4(), y4(z), z4() are "local" coords with z4 out of central panel
!         while xx, yy, zz are "true" Cartesian values
!         xx is new x after rot by rlong0 then rlat0
          xx=rotpole(1,1)*x4_iq_m+rotpole(1,2)*y4_iq_m+
     .       rotpole(1,3)*z4_iq_m
          yy=rotpole(2,1)*x4_iq_m+rotpole(2,2)*y4_iq_m+
     .       rotpole(2,3)*z4_iq_m
          zz=rotpole(3,1)*x4_iq_m+rotpole(3,2)*y4_iq_m+
     .       rotpole(3,3)*z4_iq_m
        endif
        rlat4(iq,m)=asin(zz)
        if(yy.ne.0..or.xx.ne.0.)then
          rlong4(iq,m)=atan2(yy,xx)                       ! N.B. -pi to pi
          if(rlong4(iq,m).lt.0.)rlong4(iq,m)=rlong4(iq,m)+2.*pi ! 0 to 2*pi  09-25-1997
        else
          rlong4(iq,m)=0.    ! a default value for NP/SP
        endif
!       convert long4 and lat4 (used by cctocc4) to degrees	
        rlat4(iq,m)=rlat4(iq,m)*180./pi
        rlong4(iq,m)=rlong4(iq,m)*180./pi
       enddo   ! iq loop
      enddo    ! m loop

      dsfact=4*ikk/(2.*pi)     ! con-cube
      ds=rearth/dsfact
c     extend em4 to uppermost i and j rows
      do j=1,4*ikk
       em4(iquadx,j)=em4(1,j)
      enddo
      do i=1,4*ikk
       em4(i,iquadx)=em4(i,1)
      enddo
      do j=1,ikk
       do i=1,ikk
        do n=0,5
c        average Purser em is pi/2
         em_g(indx(i,j,n,ikk,ikk))=pi/(2.*em4(4*i-1,4*j-1))
         emu_g(indx(i,j,n,ikk,ikk))=pi/(2.*em4(4*i+1,4*j-1))
         emv_g(indx(i,j,n,ikk,ikk))=pi/(2.*em4(4*i-1,4*j+1))
        enddo ! n loop
        ax(indx(i,j,0,ikk,ikk))=ax4(4*i-1,4*j-1)
        ay(indx(i,j,0,ikk,ikk))=ay4(4*i-1,4*j-1)
        az(indx(i,j,0,ikk,ikk))=az4(4*i-1,4*j-1)
       enddo  ! i loop
      enddo   ! j loop
      if(ktau<=1)then
        print *,'ax6 (1,1,0) & (2,2,0) ',ax(indx(1,1,0,ikk,ikk)),
     &                                   ax(indx(2,2,0,ikk,ikk))
        print *,'ay6 (1,1,0) & (2,2,0) ',ay(indx(1,1,0,ikk,ikk)),
     &                                   ay(indx(2,2,0,ikk,ikk))
        print *,'az6 (1,1,0) & (2,2,0) ',az(indx(1,1,0,ikk,ikk)),
     &                                   az(indx(2,2,0,ikk,ikk))
      endif ! (ktau<=1)
         
3     do j=1,ikk
       do i=1,ikk
        xx=xx4(4*i-1,4*j-1)
        yy=yy4(4*i-1,4*j-1)

c       set up x0, y0, z0 coords on cube -1 to 1
c       avoids earlier equivalencing of x,x0  etc
        x(indx(i,j,0,ikk,ikk))= 1.
        y(indx(i,j,0,ikk,ikk))=xx
        z(indx(i,j,0,ikk,ikk))=yy
        x(indx(i,j,3,ikk,ikk))=-1.
        z(indx(i,j,3,ikk,ikk))=-xx
        y(indx(i,j,3,ikk,ikk))=-yy
        x(indx(i,j,1,ikk,ikk))=-yy
        y(indx(i,j,1,ikk,ikk))=xx
        z(indx(i,j,1,ikk,ikk))= 1.
        y(indx(i,j,4,ikk,ikk))=-yy
        x(indx(i,j,4,ikk,ikk))=xx
        z(indx(i,j,4,ikk,ikk))=-1.
        x(indx(i,j,2,ikk,ikk))=-yy
        y(indx(i,j,2,ikk,ikk))= 1.
        z(indx(i,j,2,ikk,ikk))=-xx
        z(indx(i,j,5,ikk,ikk))=yy
        y(indx(i,j,5,ikk,ikk))=-1.
        x(indx(i,j,5,ikk,ikk))=xx
c       if(i.eq.(ikk+1)/2.and.j.eq.(ikk+1)/2)print *,'xx,yy: ',xx,yy
       enddo  ! i loop
      enddo   ! j loop
      do iq=1,ik*ik*6
       call norm8(x(iq),y(iq),z(iq),den1) ! x, y, z are coords on sphere  -1 to 1
      enddo   ! iq loop
      if(ktau.eq.0)print *,'basic grid length ds =',ds

      if(schmidt.ne.1.)then
        alf=(one-schmidt**2)/(one+schmidt**2)
        print *,'doing schmidt with schmidt,alf: ',schmidt,alf
        do iq=1,ik*ik*6
         xin=x(iq)
         yin=y(iq)
         zin=z(iq)
         x(iq)=xin*schmidt*(1.+alf)/(1.+alf*zin)
         y(iq)=yin*schmidt*(1.+alf)/(1.+alf*zin)
         z(iq)=(alf+zin)/(1.+alf*zin)
         if(schmidtin>0.)em_g(iq)=em_g(iq)*schmidt*(1.+alf*zin)/(1.-alf)
        enddo   ! iq loop

        if(schmidtin>0.)then
         do iq=1,ifull_g
!         with schmidt, for ntang=1 or 2 must average em_g to get emu_g & emv_g
          emu_g(iq)=.5*(em_g(iq)+em_g(ie_g(iq)))
          emv_g(iq)=.5*(em_g(iq)+em_g(in_g(iq)))
         enddo   ! iq loop
        endif    ! (schmidtin>0.)
      endif      ! (schmidt.ne.1.)

      print *,'ktau,ikk,schmidtin ',ktau,ikk,schmidtin
#ifdef debug
      if(ndiag.eq.2)call printp('x   ', x)
      if(ndiag.eq.2)call printp('y   ', y)
      if(ndiag.eq.2)call printp('z   ', z)
#endif
      if(ktau.eq.0.and.schmidtin>0.)then
        print *,'On each panel (ntang=0)_em_g for ',
     .          '(1,1),(1,2),(1,3),(2,2),(3,2),(ic,ic),(ikk,ikk)'
        do n=0,npanels
         iq11=indx(1,1,n,ikk,ikk)
         iq12=indx(1,2,n,ikk,ikk)
         iq13=indx(1,3,n,ikk,ikk)
         iq22=indx(2,2,n,ikk,ikk)
         iq32=indx(3,2,n,ikk,ikk)
         iqcc=indx((ikk+1)/2,(ikk+1)/2,n,ikk,ikk)
         iqnn=indx(ikk,ikk,n,ikk,ikk)
         print '(i3,7f8.3)',n,em_g(iq11),em_g(iq12),em_g(iq13),
     .                        em_g(iq22),em_g(iq32),em_g(iqcc),
     &                        em_g(iqnn)
        enddo
        print *,'On each panel (ntang=0)_emu_g for ',
     .          '(1,1),(1,2),(1,3),(2,2),(3,2),(ic,ic),(ikk,ikk)'
        do n=0,npanels
         iq11=indx(1,1,n,ikk,ikk)
         iq12=indx(1,2,n,ikk,ikk)
         iq13=indx(1,3,n,ikk,ikk)
         iq22=indx(2,2,n,ikk,ikk)
         iq32=indx(3,2,n,ikk,ikk)
         iqcc=indx((ikk+1)/2,(ikk+1)/2,n,ikk,ikk)
         iqnn=indx(ikk,ikk,n,ikk,ikk)
         print '(i3,7f8.3)',n,emu_g(iq11),emu_g(iq12),emu_g(iq13),
     .                        emu_g(iq22),emu_g(iq32),emu_g(iqcc),
     &                        emu_g(iqnn)
        enddo
        print *,'On each panel (ntang=0)_emv_g for ',
     .          '(1,1),(1,2),(1,3),(2,2),(3,2),(ic,ic),(ikk,ikk)'
        do n=0,npanels
         iq11=indx(1,1,n,ikk,ikk)
         iq12=indx(1,2,n,ikk,ikk)
         iq13=indx(1,3,n,ikk,ikk)
         iq22=indx(2,2,n,ikk,ikk)
         iq32=indx(3,2,n,ikk,ikk)
         iqcc=indx((ikk+1)/2,(ikk+1)/2,n,ikk,ikk)
         iqnn=indx(ikk,ikk,n,ikk,ikk)
         print '(i3,7f8.3)',n,emv_g(iq11),emv_g(iq12),emv_g(iq13),
     .                        emv_g(iq22),emv_g(iq32),emv_g(iqcc),
     &                        emv_g(iqnn)
        enddo
      endif  ! (ktau.eq.0.and.schmidtin>0.)

c     set up vectors in direction of u and v
      if(schmidtin<0.)then
c       same as below but avoids recalc in_g, ie_g, iw_g, is_g)  
        do n=0,5,2
         n_w=mod(n+5,6)
         n_e=mod(n+2,6)
         n_n=mod(n+1,6)
         n_s=mod(n+4,6)
         do j=1,ikk
          do i=1,ikk
           iq=indx(i,j,n,ikk,ikk)
           iqp=iq+1
           iqm=iq-1
           if(i==1)iqm=indx(ikk,j,n_w,ikk,ikk)
           if(i==ikk)iqp=indx(ikk+1-j,1,n_e,ikk,ikk)
           ax(iq)=x(iqp)-x(iqm)
           az(iq)=z(iqp)-z(iqm)
           ay(iq)=y(iqp)-y(iqm)
           iqp=iq+ikk
           iqm=iq-ikk
           if(j==1)iqm=indx(ikk,ikk+1-i,n_s,ikk,ikk)
           if(j==ikk)iqp=indx(i,1,n_n,ikk,ikk)
           bx(iq)=x(iqp)-x(iqm)
           by(iq)=y(iqp)-y(iqm)
           bz(iq)=z(iqp)-z(iqm)
          enddo
         enddo
        enddo
        do n=1,5,2
         n_w=mod(n+4,6)
         n_e=mod(n+1,6)
         n_n=mod(n+2,6)
         n_s=mod(n+5,6)
         do j=1,ikk
          do i=1,ikk
           iq=indx(i,j,n,ikk,ikk)
           iqp=iq+1
           iqm=iq-1
           if(i==1)iqm=indx(ikk+1-j,ikk,n_w,ikk,ikk)
           if(i==ikk)iqp=indx(1,j,n_e,ikk,ikk)
           ax(iq)=x(iqp)-x(iqm)
           ay(iq)=y(iqp)-y(iqm)
           az(iq)=z(iqp)-z(iqm)
           iqp=iq+ikk
           iqm=iq-ikk
           if(j==1)iqm=indx(i,ikk,n_s,ikk,ikk)
           if(j==ikk)iqp=indx(1,ikk+1-i,n_n,ikk,ikk)
           bx(iq)=x(iqp)-x(iqm)
           by(iq)=y(iqp)-y(iqm)
           bz(iq)=z(iqp)-z(iqm)
          enddo
         enddo
        enddo
      else  !  usual with (schmidtin>0. (but equiv. to above)
        do iq=1,ifull_g
c        first guess tang vectors by finite differences
         ax(iq)=x(ie_g(iq))-x(iw_g(iq))
         ay(iq)=y(ie_g(iq))-y(iw_g(iq))
         az(iq)=z(ie_g(iq))-z(iw_g(iq))
         bx(iq)=x(in_g(iq))-x(is_g(iq))
         by(iq)=y(in_g(iq))-y(is_g(iq))
         bz(iq)=z(in_g(iq))-z(is_g(iq))
        enddo   ! iq loop
      endif   ! (schmidtin<0.  ... else)
c     do iq=1,ikk*ikk*6
c      print *,'iq,ax,ay,az: ',iq,ax(iq),ay(iq),az(iq)
c      print *,'iq,bx,by,bz: ',iq,bx(iq),by(iq),bz(iq)
c     enddo
c       form axx and bxx tangential to the sphere
        call cross3b(axx,ayy,azz, bx,by,bz, x,y,z, ik*ik*6)
        call cross3(bxx,byy,bzz, x,y,z, ax,ay,az, ik*ik*6)
        do iq=1,ikk*ikk*6
         call norm(axx(iq),ayy(iq),azz(iq),den)
         call norm(bxx(iq),byy(iq),bzz(iq),den)
c        make sure they are perpendicular & normalize
         dot=axx(iq)*bxx(iq)+ayy(iq)*byy(iq)+azz(iq)*bzz(iq)
         eps=-dot/(1.+sqrt(1.-dot*dot))
         ax(iq)=axx(iq)+eps*bxx(iq)
         ay(iq)=ayy(iq)+eps*byy(iq)
         az(iq)=azz(iq)+eps*bzz(iq)
         bx(iq)=bxx(iq)+eps*axx(iq)
         by(iq)=byy(iq)+eps*ayy(iq)
         bz(iq)=bzz(iq)+eps*azz(iq)
         call norm(ax(iq),ay(iq),az(iq),den)
         call norm(bx(iq),by(iq),bz(iq),den)
        enddo   ! iq loop
        if(schmidtin<0.)return  ! finish of stuff needed for onthefly
        
          do iq=1,ifull_g
!          calculate inverse of emu_g & emv_g first
           dx2=(x(ie_g(iq))-x(iq))**2+(y(ie_g(iq))-y(iq))**2
     .                             +(z(ie_g(iq))-z(iq))**2
!          include arc-length corrn using 2*arcsin(theta/2)
           emu_g(iq)=sqrt(dx2)*(1.+dx2/24.) *dsfact
           dy2=(x(in_g(iq))-x(iq))**2+(y(in_g(iq))-y(iq))**2
     .                             +(z(in_g(iq))-z(iq))**2
           emv_g(iq)=sqrt(dy2)*(1.+dy2/24.) *dsfact
          enddo   ! iq loop
          do iq=1,ifull_g   ! based on inverse values of emu_g & emv_g
           if (isv2_g(iq)<1.and.iwu2_g(iq)>ifull_g) then        ! MJT bug fix
           em_g(iq)=4./(emv_g(iwu2_g(iq)-ifull_g)+emu_g(iq)+
     .                  emu_g(isv2_g(iq)+ifull_g)+emv_g(iq))    ! MJT bug fix
           else if (isv2_g(iq)<1) then                          ! MJT bug fix
           em_g(iq)=4./(emu_g(iwu2_g(iq))+emu_g(iq)+
     .                  emu_g(isv2_g(iq)+ifull_g)+emv_g(iq))    ! MJT bug fix
           else if (iwu2_g(iq)>ifull_g) then                    ! MJT bug fix
           em_g(iq)=4./(emv_g(iwu2_g(iq)-ifull_g)+emu_g(iq)+
     .                  emv_g(isv2_g(iq))+emv_g(iq))            ! MJT bug fix
           else                                                 ! MJT bug fix
           em_g(iq)=4./(emu_g(iwu2_g(iq))+emu_g(iq)+
     .                  emv_g(isv2_g(iq))+emv_g(iq))            ! MJT bug fix
           end if                                               ! MJT bug fix
c          experimental option follows - only tiniest difference for ikk=20
c          em_g(iq)=2./sqrt((emu_g(iwu2_g(iq))+emu_g(iq))*
c    .                  (emv_g(isv2_g(iq))+emv_g(iq)))
          enddo   ! iq loop
          do iq=1,ifull_g
           emu_g(iq)=1./emu_g(iq)
           emv_g(iq)=1./emv_g(iq)
          enddo   ! iq loop
      
      if(ktau.eq.0)then
        do iq=ikk-2,ikk
         print *,'iq,em_g,emu_g,emv_g',iq,em_g(iq),emu_g(iq),emv_g(iq)
        enddo   ! iq loop
        if(id.le.ikk.and.jd.le.jl)then
          iq=id+ikk*(jd-1)
          print *,'values at idjd'
          print *,'iq,x,y,z',iq,x(iq),y(iq),z(iq)
          print *,'iq,ax,ay,az',iq,ax(iq),ay(iq),az(iq)
          print *,'iq,bx,by,bz',iq,bx(iq),by(iq),bz(iq)
          print *,'values at in_g(idjd)'
          print *,'iq,x,y,z',in_g(iq),x(in_g(iq)),y(in_g(iq)),
     &                       z(in_g(iq))
          print *,'iq,ax,ay,az',in_g(iq),ax(in_g(iq)),ay(in_g(iq)),
     &                          az(in_g(iq))
          print *,'iq,bx,by,bz',in_g(iq),bx(in_g(iq)),by(in_g(iq)),
     &                          bz(in_g(iq))
          print *,'values at ie_g(idjd)'
          print *,'iq,x,y,z',ie_g(iq),x(ie_g(iq)),y(ie_g(iq)),
     &                       z(ie_g(iq))
          print *,'iq,ax,ay,az',ie_g(iq),ax(ie_g(iq)),ay(ie_g(iq)),
     &                          az(ie_g(iq))
          print *,'iq,bx,by,bz',ie_g(iq),bx(ie_g(iq)),by(ie_g(iq)),
     &                          bz(ie_g(iq))
          print *,'values at iw_g(idjd)'
          print *,'iq,x,y,z',iw_g(iq),x(iw_g(iq)),y(iw_g(iq)),
     &                       z(iw_g(iq))
          print *,'iq,ax,ay,az',iw_g(iq),ax(iw_g(iq)),ay(iw_g(iq)),
     &                          az(iw_g(iq))
          print *,'iq,bx,by,bz',iw_g(iq),bx(iw_g(iq)),by(iw_g(iq)),
     &                          bz(iw_g(iq))
          print *,'values at is_g(idjd)'
          print *,'iq,x,y,z',is_g(iq),x(is_g(iq)),y(is_g(iq)),
     &                       z(is_g(iq))
          print *,'iq,ax,ay,az',is_g(iq),ax(is_g(iq)),ay(is_g(iq)),
     &                          az(is_g(iq))
          print *,'iq,bx,by,bz',is_g(iq),bx(is_g(iq)),by(is_g(iq)),
     &                          bz(is_g(iq))
        endif
      endif  ! (ktau.eq.0)

c     calculate approx areas around each grid point
c     just used for error diagnostics
c     now use 1/(em_g**2) to cope with schmidt, rotated and ocatagon coordinates
      sumwts=0.
      do iq=1,ifull_g
       wts(iq)=1./em_g(iq)**2
       sumwts=sumwts+wts(iq)
!c     cosa is dot product of unit vectors
!c     *** only useful as diagnostic for gnew
!      cosa(iq)=ax(iq)*bx(iq)+ay(iq)*by(iq)+az(iq)*bz(iq)
      enddo   ! iq loop
      if(ktau.eq.0)then
        print *,'sumwts/ifull_g ',sumwts/ifull_g  ! ideally equals 4*pi ??
        print *,'in setxyz rlong0,rlat0,schmidt ',rlong0,rlat0,schmidt
      endif  ! (ktau.eq.0)

c     previously calculates approx areas around each grid point using modulus of
c     cross-products; used for error diagnostics
c     do j=1,ikk
c      do i=1,ikk
c       iq=ind(i,j,0)
c       wts(iq)=.5*(crossmod(iq,4*i-3,4*j+1,4*i+1,4*j+1)  ! nw,ne  with xx4,yy4
c    .             +crossmod(iq,4*i+1,4*j+1,4*i+1,4*j-3)  ! ne,se  with xx4,yy4
c    .             +crossmod(iq,4*i+1,4*j-3,4*i-3,4*j-3)  ! se,sw  with xx4,yy4
c    .             +crossmod(iq,4*i-3,4*j-3,4*i-3,4*j+1)) ! sw,nw  with xx4,yy4
c       sumwts=sumwts+wts(iq)
c       do n=1,5
c        wts(iq+n*ikk*ikk)=wts(iq)
c       enddo  ! n loop
c      enddo   ! i loop
c     enddo    ! j loop

      do iq=1,ifull_g
c      scale wts so sum over globe is 1.
c      wts(iq)=wts(iq)/(6.*sumwts)  ! for old conf-cub defn
       wts(iq)=wts(iq)/sumwts
c      also provide latitudes and longitudes (-pi to pi)
       if(rlong0.eq.0..and.rlat0.eq.90.)then
         xx=x(iq)
         yy=y(iq)
         zz=z(iq)
       else
!        x(), y(z), z() are "local" coords with z out of central panel
!        while xx, yy, zz are "true" Cartesian values
!        xx is new x after rot by rlong0 then rlat0
         xx=rotpole(1,1)*x(iq)+rotpole(1,2)*y(iq)+rotpole(1,3)*z(iq)
         yy=rotpole(2,1)*x(iq)+rotpole(2,2)*y(iq)+rotpole(2,3)*z(iq)
         zz=rotpole(3,1)*x(iq)+rotpole(3,2)*y(iq)+rotpole(3,3)*z(iq)
       endif
c      f_g(iq)=2. *2.*pi *(z(iq)/rdiv) /86400.
       rlatt_g(iq)=asin(zz)
       f_g(iq)=2. *2.*pi *zz /86400.  !  zz along "true" N-S axis
       if(yy.ne.0..or.xx.ne.0.)then
         rlongg_g(iq)=atan2(yy,xx)                       ! N.B. -pi to pi
         if(rlongg_g(iq).lt.0.)rlongg_g(iq)=rlongg_g(iq)+2.*pi ! 0 to 2*pi  09-25-1997
       else
         rlongg_g(iq)=0.    ! a default value for NP/SP
       endif
c      if(iq.eq.idjd)print *,'iq,x,y,z,xx,yy,zz,long,lat ',
c    .  iq,x(iq),y(iq),z(iq),xx,yy,zz,
c    .  rlongg_g(iq)*180./pi,rlatt_g(iq)*180./pi
      enddo   ! iq loop
#ifdef debug
      if(ndiag.eq.2)then
!       do iq=1,ifull_g
!        cosa(iq)=100.*wts(iq)
!       enddo
!       call printp('wts ',cosa)
        call printp('lat ',rlat)
        call printp('long',rlong)
      endif
#endif
      if(ktau.eq.0)then
        print *,'At centre of the faces:'
        do n=0,npanels
         iq=indx((ikk+1)/2,(ikk+1)/2,n,ikk,ikk)
         print '('' n,iq,x,y,z,long,lat,f ''i2,i7,3f7.3,2f8.2,f9.5)',n,
     .     iq,x(iq),y(iq),z(iq),
     .     rlongg_g(iq)*180./pi,rlatt_g(iq)*180./pi,f_g(iq)
        enddo
        print *,'At mid-x along edges:'
        do n=0,npanels
         iq=indx((ikk+1)/2,1,n,ikk,ikk)
         print '('' n,iq,x,y,z,long,lat,f_g ''i2,i7,3f7.3,2f8.2,f9.5)',
     &     n,iq,x(iq),y(iq),z(iq),
     .     rlongg_g(iq)*180./pi,rlatt_g(iq)*180./pi,f_g(iq)
        enddo
        print *,'At mid-y along edges:'
        do n=0,npanels
         iq=indx(1,(ikk+1)/2,n,ikk,ikk)
         print '('' n,iq,x,y,z,long,lat,f_g ''i2,i7,3f7.3,2f8.2,f9.5)',
     &     n,iq,x(iq),y(iq),z(iq),
     .     rlongg_g(iq)*180./pi,rlatt_g(iq)*180./pi,f_g(iq)
        enddo
        print *,'On each panel final_em_g for ',
     .          '(1,1),(1,2),(1,3),(2,2),(3,2),(ic,ic),(ikk,ikk)'
        do n=0,npanels
         iq11=indx(1,1,n,ikk,ikk)
         iq12=indx(1,2,n,ikk,ikk)
         iq13=indx(1,3,n,ikk,ikk)
         iq22=indx(2,2,n,ikk,ikk)
         iq32=indx(3,2,n,ikk,ikk)
         iqcc=indx((ikk+1)/2,(ikk+1)/2,n,ikk,ikk)
         iqnn=indx(ikk,ikk,n,ikk,ikk)
         print '(i3,7f8.3)',n,em_g(iq11),em_g(iq12),em_g(iq13),
     .                        em_g(iq22),em_g(iq32),em_g(iqcc),
     &                        em_g(iqnn)
        enddo
        print *,'On each panel final_emu_g for ',
     .          '(1,1),(1,2),(1,3),(2,2),(3,2),(ic,ic),(ikk,ikk)'
        do n=0,npanels
         iq11=indx(1,1,n,ikk,ikk)
         iq12=indx(1,2,n,ikk,ikk)
         iq13=indx(1,3,n,ikk,ikk)
         iq22=indx(2,2,n,ikk,ikk)
         iq32=indx(3,2,n,ikk,ikk)
         iqcc=indx((ikk+1)/2,(ikk+1)/2,n,ikk,ikk)
         iqnn=indx(ikk,ikk,n,ikk,ikk)
         print '(i3,7f8.3)',n,emu_g(iq11),emu_g(iq12),emu_g(iq13),
     .                        emu_g(iq22),emu_g(iq32),emu_g(iqcc),
     &                        emu_g(iqnn)
        enddo
        print *,'On each panel final_emv_g for ',
     .          '(1,1),(1,2),(1,3),(2,2),(3,2),(ic,ic),(ikk,ikk)'
        do n=0,npanels
         iq11=indx(1,1,n,ikk,ikk)
         iq12=indx(1,2,n,ikk,ikk)
         iq13=indx(1,3,n,ikk,ikk)
         iq22=indx(2,2,n,ikk,ikk)
         iq32=indx(3,2,n,ikk,ikk)
         iqcc=indx((ikk+1)/2,(ikk+1)/2,n,ikk,ikk)
         iqnn=indx(ikk,ikk,n,ikk,ikk)
         print '(i3,7f8.3)',n,emv_g(iq11),emv_g(iq12),emv_g(iq13),
     .                        emv_g(iq22),emv_g(iq32),emv_g(iqcc),
     &                        emv_g(iqnn)
        enddo
      endif  ! (ktau.eq.0)
      do iq=1,ifull_g   ! set up Coriolis
       fu_g(iq)=(f_g(iq)+f_g(ie_g(iq)))*.5
       fv_g(iq)=(f_g(iq)+f_g(in_g(iq)))*.5
      enddo   ! iq loop
      do iq=1,ifull_g   ! average map factor derivs needed for nxmap=1
       dmdx_g(iq)=.5*(em_g(ie_g(iq))-em_g(iw_g(iq)))/ds  
       dmdy_g(iq)=.5*(em_g(in_g(iq))-em_g(is_g(iq)))/ds  
      enddo   ! iq loop

      ratmin=100.
      ratmax=0.
      do n=0,npanels
        do i=1,ikk
         iq=indx(i,ikk/2,n,ikk,ikk)
         rat=em_g(iq)/em_g(ie_g(iq))
         if(rat<ratmin)then
           ratmin=rat
           imin=i
         endif
         if(rat>ratmax)then
           ratmax=rat
           imax=i
         endif
        enddo
        if(num==1)then
          print *,'em_g ratio for j=ikk/2 on npanel ',n
          write (6,"(12f6.3)")
     &          (em_g(indx(i,ikk/2,n,ikk,ikk))
     &          /em_g(ie_g(indx(i,ikk/2,n,ikk,ikk))),
     &          i=1,ikk)
        endif
      enddo
      print *,'for j=ikk/2 & myid=0, ratmin,ratmax = ',ratmin,ratmax
      print *,'with imin,imax ',imin,imax
      ratmin=100.
      ratmax=0.
      do n=0,npanels
       do j=1,ikk
        do i=1,ikk
         iq=indx(i,j,n,ikk,ikk)
         rat=em_g(iq)/em_g(ie_g(iq))
         if(rat<ratmin)then
          ratmin=rat
          imin=i
          jmin=j
         endif
         if(rat>ratmax)then
          ratmax=rat
          imax=i
          jmax=j
         endif
        enddo
       enddo
      enddo
      print *,'for all j & myid=0, ratmin,ratmax = ',ratmin,ratmax
      print *,'with imin,jmin,imax,jmax ',imin,jmin,imax,jmax
      write (6,"('1st 10 ratios',10f6.3)") (em_g(iq)/em_g(ie_g(iq)),
     &                                     iq=1,10)
      numpts=0
      do iq=1,ifull_g
       rlatdeg=rlatt_g(iq)*180./pi
       rlondeg=rlongg_g(iq)*180./pi
       if(rlatdeg>20..and.rlatdeg<60.
     &    .and.rlondeg>230..and.rlatdeg<300.)numpts=numpts+1
      enddo
      print *,'points in SGMIP region ',numpts

      return
      end

      subroutine indv(iq,i,j,n)
c     calculates simple i,j,n indices from supplied iq
      include 'newmpar.h'
      n=(iq-1)/(il_g*il_g)
      j=1+(iq-n*il_g*il_g-1)/il_g
      i=iq-(j-1)*il_g-n*il_g*il_g
      return
      end
      subroutine norm(a,b,c,den)
      den=sqrt(a**2+b**2+c**2)
      a=a/den
      b=b/den
      c=c/den
      return
      end
      subroutine norm8(a,b,c,den)
      real(kind=8) a,b,c,den
      den=sqrt(a**2+b**2+c**2)
      a=a/den
      b=b/den
      c=c/den
      return
      end
      subroutine vecpanel(ax6,ay6,az6)
c     define vectors on panels 1:5 from panel 0
      include 'newmpar.h'
      dimension ax6(il_g,il_g,0:5),ay6(il_g,il_g,0:5),az6(il_g,il_g,0:5)
      do j=1,il_g
       do i=1,il_g
        a1=ax6(i,j,0)
        a2=ay6(i,j,0)
        a3=az6(i,j,0)
        ax6(i,j,1)=-a3
        ay6(i,j,1)=a2
        az6(i,j,1)=a1
        ax6(i,j,2)=-a3
        ay6(i,j,2)=a1
        az6(i,j,2)=-a2
        ax6(i,j,3)=-a1
        ay6(i,j,3)=-a3
        az6(i,j,3)=-a2
        ax6(i,j,4)=a2
        ay6(i,j,4)=-a3
        az6(i,j,4)=-a1
        ax6(i,j,5)=a2
        ay6(i,j,5)=-a1
        az6(i,j,5)=a3
       enddo  ! i loop
      enddo   ! j loop
      return
      end

      subroutine printp(name,s6)
      include 'newmpar.h'
      character(len=4) name
      dimension s6(il_g,il_g,0:5)  ! no longer access s(ifull_g-1), s(ifull_g)
      dimension s1f(0:il_g+1,3*il_g),s2f(0:il_g+1,3*il_g)

c     s1 is Grenwich-NP section i.e.  0-1-3
c     s2 is Oz-SP section i.e.  2-4-5
      call strip2(s6,s6,s1f,s2f)
      print *, name,'  013'
        do j=3*il_g,1,-1
         print 9,j,(s1f(i,j),i=0,il_g+1)
        enddo
9        format(i3,1x,21f6.3)
      print *
      print *, name,'  245'
        do j=3*il_g,1,-1
         print 9,j,(s2f(i,j),i=0,il_g+1)
        enddo
      return
      end

      subroutine strip2(s,s6,s1,s2)
      use cc_mpi, only : indx
      use indices_m
      implicit none
      include 'newmpar.h'
      integer i,j
!      dimension in6(il_g,il_g,0:5),is6(il_g,il_g,0:5),iw6(il_g,il_g,0:5)
!     .         ,ie6(il_g,il_g,0:5)
!      equivalence (in_g,in6),(is_g,is6),(iw_g,iw6),(ie_g,ie6)
      real s(ifull_g),s6(il_g,il_g,0:5)
c     N.B.  s & s6 are equivalenced via the call
c     dimension s1f(0:il_g+1,3*il_g),s2f(0:il_g+1,3*il_g)
      real s1(0:il_g+1,il_g,3),s2(0:il_g+1,il_g,3)  ! equiv to s1f, s2f via call
c     s1 is Grenwich-NP section i.e.  0-1-3
c     s2 is Oz-SP section i.e.  2-4-5
c     for gnewst, these are extended on the sides only (i=0 & il_g+1)
      do j=1,il_g
       do i=1,il_g
        s1(i,j,1)=s6(i,j,0)
        s1(i,j,2)=s6(i,j,1)
        s1(i,j,3)=s6(j,il_g+1-i,3)
        s2(i,j,1)=s6(j,il_g+1-i,2)
        s2(i,j,2)=s6(i,j,4)
        s2(i,j,3)=s6(i,j,5)
       enddo  ! i loop
       
       s1(0,j,1)=s(iw_g(indx(1,j,0,il_g,il_g)))
c      print *,'j,iw6(1,j,0),s1(0,j,1) ',j,iw6(1,j,0),s1(0,j,1)
       s1(0,j,2)=s(iw_g(indx(1,j,1,il_g,il_g)))
       s1(0,j,3)=s(in_g(indx(j,il_g,3,il_g,il_g)))
       s2(0,j,1)=s(in_g(indx(j,il_g,2,il_g,il_g)))
       s2(0,j,2)=s(iw_g(indx(1,j,4,il_g,il_g)))
       s2(0,j,3)=s(iw_g(indx(1,j,5,il_g,il_g)))
       s1(il_g+1,j,1)=s(ie_g(indx(il_g,j,0,il_g,il_g)))
       s1(il_g+1,j,2)=s(ie_g(indx(il_g,j,1,il_g,il_g)))
       s1(il_g+1,j,3)=s(is_g(indx(j,1,3,il_g,il_g)))
       s2(il_g+1,j,1)=s(is_g(indx(j,1,2,il_g,il_g)))
       s2(il_g+1,j,2)=s(ie_g(indx(il_g,j,4,il_g,il_g)))
       s2(il_g+1,j,3)=s(ie_g(indx(il_g,j,5,il_g,il_g)))
      enddo   ! j loop
      return
      end

      subroutine cross3(c1,c2,c3,a1,a2,a3,b1,b2,b3, ifull_g)
c     calculate vector components of c = a x b
c     where each RHS component represents 3 vector components
c     this one need not have contiguous memory in common
c     include 'newmpar_gx.h'
      real(kind=8)    a1(ifull_g),a2(ifull_g),a3(ifull_g)
      dimension b1(ifull_g),b2(ifull_g),b3(ifull_g)
      dimension c1(ifull_g),c2(ifull_g),c3(ifull_g)
      do i=1,ifull_g
       c1(i)=a2(i)*b3(i)-b2(i)*a3(i)
       c2(i)=a3(i)*b1(i)-b3(i)*a1(i)
       c3(i)=a1(i)*b2(i)-b1(i)*a2(i)
      enddo
      return
      end

      subroutine cross3b(c1,c2,c3,a1,a2,a3,b1,b2,b3, ifull_g)
c     calculate vector components of c = a x b
c     where each RHS component represents 3 vector components
c     this one need not have contiguous memory in common
c     include 'newmpar_gx.h'
      dimension a1(ifull_g),a2(ifull_g),a3(ifull_g)
      real(kind=8)    b1(ifull_g),b2(ifull_g),b3(ifull_g)
      dimension c1(ifull_g),c2(ifull_g),c3(ifull_g)
      do i=1,ifull_g
       c1(i)=a2(i)*b3(i)-b2(i)*a3(i)
       c2(i)=a3(i)*b1(i)-b3(i)*a1(i)
       c3(i)=a1(i)*b2(i)-b1(i)*a2(i)
      enddo
      return
      end


