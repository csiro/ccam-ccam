      subroutine onthefly(kdate_r,ktime_r,
     .   psl,zss,tss,wb,wbice,snowd,sicedep,
     .   t,u,v,qg,tgg,
     .   tggsn,smass,ssdn, ssdnn,osnowd,snage,isflag,nested)
!     Target points use values interpolated to their 4 grid "corners";
!     these corner values are then averaged to the grid centres
!     Called by either indata or nestin
!     nested=0  for calls from indata; 1  for calls from nestin     

      parameter (ntest=0)
      parameter (nord=3)   ! 1 for bilinear, 3 for bicubic
!     related to cctocc4                       

!     Note: 1) The arrays are replaced in place
!           2) kl is assumed to be the same for both grids
      include 'newmpar.h'
      include 'bigxy4.h' ! common/bigxy4/xx4(iquad,iquad),yy4(iquad,iquad)
      include 'const_phys.h'
      include 'latlong.h'  ! rlatt,rlongg,
      include 'map.h'  ! zs,land & used for giving info after all setxyz
      include 'parm.h'
      include 'parm_nqg.h'  ! nqg_r,nqg_set
      include 'sigs.h'
      include 'stime.h'   ! kdate_s,ktime_s  sought values for data read
      include 'tracers.h'
      include 'vecsuv.h'   ! vecsuv info
      include 'vvel.h'
      include 'xyzinfo.h'  ! x,y,z,wts
      common/schmidtx/rlong0x,rlat0x,schmidtx ! infile, newin, nestin, indata
      common/sigin/sigin(kl),kk  ! for vertint, infile
      common/work3b/rlong4(ifull,4),rlat4(ifull,4),   ! shared with setxyz
     .              aa(ifull),bb(ifull),
     .              xx4_sav(iquad,iquad),yy4_sav(iquad,iquad),
     .              dum3d(2*ijk-10*ifull-2*iquad*iquad) 
      common/work3f/nface(ifull),xg(ifull),yg(ifull),
     .              uct(ifull),vct(ifull),wct(ifull),
     .              ax_s(ifull),ay_s(ifull),az_s(ifull),
     .              bx_s(ifull),by_s(ifull),bz_s(ifull),
     .              zs_t(ifull),tss_l(ifull),tss_s(ifull),land_t(ifull),
     .              ucc(ifull),vcc(ifull),wcc(ifull),
     .              uc(ifull),vc(ifull),wc(ifull),
     .              pmsl(ifull),dum3b(3*ijk-23*ifull)
      logical land_t
      dimension xg4(ifull,4),yg4(ifull,4),nface4(ifull,4)
      real rotpoles(3,3),rotpole(3,3)

      real psl(ifull),zss(ifull),tss(ifull),
     . wb(ifull,ms),wbice(ifull,ms),snowd(ifull),sicedep(ifull),
     . t(ifull,kl),u(ifull,kl),v(ifull,kl),qg(ifull,kl),
     . tgg(ifull,ms),tggsn(ifull,3),smass(ifull,3),ssdn(ifull,3),
     . ssdnn(ifull),osnowd(ifull),snage(ifull)
      integer isflag(ifull)
c     data pi/3.1415926536/,id1/16/,jd1/79/
      data id1/3/,jd1/60/

      nqg_set=8    
!     save cc target file geometry
      rlong0_t=rlong0
      rlat0_t=rlat0
      schmidt_t=schmidt
      ds_t=ds
      do iq=1,ifull
       zs_t(iq)=zs(iq)
       land_t(iq)=land(iq)  
      enddo

c     start of processing loop 
      nemi=2   !  assume source land-mask based on tss sign first
      if(ktau.lt.3)print *,'search for kdate_s,ktime_s >= ',
     .                                 kdate_s,ktime_s
      id_t=id
      jd_t=jd
      id=id1
      jd=jd1
      idjd=id+il*(jd-1)
      idjd1=idjd

      call infile(io_in2,kdate_r,ktime_r,nem2,
     . timegb,ds,psl,aa,zss,aa,bb,
     . tss,aa,wb,wbice,aa,snowd,sicedep,
     . t,u,v,qg,tgg,
     . tggsn,smass,ssdn, ssdnn,osnowd,snage,isflag,nested)

      do k=1,kl
       sig(k)=sigin(k)
      enddo
     
      rlong0=rlong0x
      rlat0=rlat0x
      schmidt=schmidtx
      call setxyz   ! for source data geometry        ******************
!     rotpole(1,) is x-axis of rotated coords in terms of orig Cartesian
!     rotpole(2,) is y-axis of rotated coords in terms of orig Cartesian
!     rotpole(3,) is z-axis of rotated coords in terms of orig Cartesian
      coslong=cos(rlong0*pi/180.)
      sinlong=sin(rlong0*pi/180.)
      coslat=cos(rlat0*pi/180.)
      sinlat=sin(rlat0*pi/180.)
      rotpoles(1,1)=coslong*sinlat
      rotpoles(1,2)=-sinlong
      rotpoles(1,3)=coslong*coslat
      rotpoles(2,1)=sinlong*sinlat
      rotpoles(2,2)=coslong
      rotpoles(2,3)=sinlong*coslat
      rotpoles(3,1)=-coslat
      rotpoles(3,2)=0.
      rotpoles(3,3)=sinlat
      if(ktau.lt.3)then
        print *,'nfly,nord ',nfly,nord
        print *,'io_in2,kdate_r,ktime_r,ktau,ds',
     .         io_in2,kdate_r,ktime_r,ktau,ds
        print *,'ds,ds_t ',ds,ds_t
        print *,'a zss(idjd1) ',zss(idjd1)
        print *,'rotpoles:'
        do i=1,3
         print 9,(i,j,j=1,3),(rotpoles(i,j),j=1,3)
        enddo
      endif ! (ktau.lt.3)
!     save "source" ax, bx etc  - used in transforming source u & v
      do iq=1,ifull
       ax_s(iq)=ax(iq)
       ay_s(iq)=ay(iq)
       az_s(iq)=az(iq)
       bx_s(iq)=bx(iq)
       by_s(iq)=by(iq)
       bz_s(iq)=bz(iq)
c      rlat_t(iq)=rlatt(iq)*180./pi
c      rlong_t(iq)=rlongg(iq)*180./pi
      enddo
      do j=1,iquad
       do i=1,iquad
        xx4_sav(i,j)=xx4(i,j)
        yy4_sav(i,j)=yy4(i,j)
       enddo
      enddo

!     restore cc target file geometry
      rlong0=rlong0_t
      rlat0=rlat0_t
      schmidt=schmidt_t
      ds_t=ds
      id=id_t
      jd=jd_t
      idjd=id+il*(jd-1)
      call setxyz  ! for target        *********************************
      print *,'after target setxyz'
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
      if(nmaxpr.eq.1)then
        print *,'rotpole:'
        do i=1,3
         print 9,(i,j,j=1,3),(rotpole(i,j),j=1,3)
9        format(3x,2i1,5x,2i1,5x,2i1,5x,3f8.4)
        enddo
        print *,'xx4,xx4_sav,yy4,yy4_sav ',
     .           xx4(id,jd),xx4_sav(id,jd),yy4(id,jd),yy4_sav(id,jd)
      endif  ! (nmaxpr.eq.1)
!     put source values into xx4, yy4; target values into xx4_sav,yy4_sav
      do j=1,iquad
       do i=1,iquad
	 xx_s=xx4_sav(i,j)
	 xx4_sav(i,j)=xx4(i,j)
	 xx4(i,j)=xx_s
	 yy_s=yy4_sav(i,j)
	 yy4_sav(i,j)=yy4(i,j)
	 yy4(i,j)=yy_s
       enddo
      enddo

!     restore cc source file geometry for latltoij
      rlong0=rlong0x
      rlat0=rlat0x
      schmidt=schmidtx
      if(nmaxpr.eq.1)then
        print *,'before latltoij for id,jd: ',id,jd
        print *,'rlong4(1-4) ',(rlong4(idjd,m),m=1,4)
        print *,'rlat4(1-4) ',(rlat4(idjd,m),m=1,4)
        print *,'rlong0x,rlat0x,schmidtx ',rlong0x,rlat0x,schmidtx 
      endif  ! (nmaxpr.eq.1)
      do m=1,4
        do iq=1,ifull
         call latltoij(rlong4(iq,m),rlat4(iq,m),
     .                 xg4(iq,m),yg4(iq,m),nface4(iq,m))
        enddo
       enddo
	id2=nint(xg4(idjd,1))
	jd2=il*nface4(idjd,1)+nint(yg4(idjd,1))
       idjd2=id2+il*(jd2-1)
	if(nmaxpr.eq.1)then
         print *,'after latltoij giving id2,jd2: ',id2,jd2
         print *,'nface4(1-4) ',(nface4(idjd,m),m=1,4)
         print *,'xg4(1-4) ',(xg4(idjd,m),m=1,4)
         print *,'yg4(1-4) ',(yg4(idjd,m),m=1,4)
         write(6,"('wb_s(1)#  ',9f7.3)") 
     .          ((wb(ii+(jj-1)*il,1),ii=id2-1,id2+1),jj=jd2-1,jd2+1)
         write(6,"('wb_s(ms)# ',9f7.3)") 
     .          ((wb(ii+(jj-1)*il,ms),ii=id2-1,id2+1),jj=jd2-1,jd2+1)
       endif
       if(nfly.eq.2)then    ! needs pmsl in this case (preferred)
         call mslp(pmsl,psl,zss,t)  
	endif

        do k=1,kl
         do iq=1,ifull
!         first set up winds in Cartesian "source" coords
          uc(iq)=ax_s(iq)*u(iq,k) + bx_s(iq)*v(iq,k)
          vc(iq)=ay_s(iq)*u(iq,k) + by_s(iq)*v(iq,k)
          wc(iq)=az_s(iq)*u(iq,k) + bz_s(iq)*v(iq,k)
!         now convert to winds in "absolute" Cartesian components
          ucc(iq)=uc(iq)*rotpoles(1,1)+vc(iq)*rotpoles(1,2)
     .                                +wc(iq)*rotpoles(1,3)
          vcc(iq)=uc(iq)*rotpoles(2,1)+vc(iq)*rotpoles(2,2)
     .                                +wc(iq)*rotpoles(2,3)
          wcc(iq)=uc(iq)*rotpoles(3,1)+vc(iq)*rotpoles(3,2)
     .                                +wc(iq)*rotpoles(3,3)
        enddo  ! iq loop
	 if(ktau.lt.3.and.k.eq.1)then
          print *,'uc,vc,wc: ',uc(id),vc(idjd1),wc(idjd1)
          print *,'ucc,vcc,wcc: ',ucc(idjd1),vcc(idjd1),wcc(idjd1)
          print *,'calling ints4 for k= ',k
	 endif

!      interpolate all required arrays to new C-C positions
!      don't need to do map factors and Coriolis on target grid
       np=0  ! controls prints in ints4
       call ints4(t (1,k),nface4,xg4,yg4,nord)
       call ints4(qg(1,k),nface4,xg4,yg4,nord)
       call ints4(ucc,      nface4,xg4,yg4,nord)
       call ints4(vcc,      nface4,xg4,yg4,nord)
       call ints4(wcc,      nface4,xg4,yg4,nord)
c      if(nsd.eq.1)then
c        call ints4(sdot(1,k),   nface4,xg4,yg4,nord)
c      endif

c      ********************** N.B. tracers not ready yet
c      if(iltin.gt.1)then
c        do ntr=1,ntracin
c         call ints4(tr(1,k,ntr),nface4,xg4,yg4,nord)
c        enddo
c      endif
 
       do iq=1,ifull
!       now convert to "target" Cartesian components (transpose used)
        uct(iq)=ucc(iq)*rotpole(1,1)+vcc(iq)*rotpole(2,1)
     .                           +wcc(iq)*rotpole(3,1)
        vct(iq)=ucc(iq)*rotpole(1,2)+vcc(iq)*rotpole(2,2)
     .                           +wcc(iq)*rotpole(3,2)
        wct(iq)=ucc(iq)*rotpole(1,3)+vcc(iq)*rotpole(2,3)
     .                           +wcc(iq)*rotpole(3,3)
!       then finally to "target" local x-y components
        u(iq,k)=ax(iq)*uct(iq) +ay(iq)*vct(iq) +az(iq)*wct(iq)
        v(iq,k)=bx(iq)*uct(iq) +by(iq)*vct(iq) +bz(iq)*wct(iq)
       enddo  ! iq loop
	if(ktau.lt.3.and.k.eq.1)then
         print *,'interp. ucc,vcc,wcc: ',ucc(idjd),vcc(idjd),wcc(idjd)
         print *,'uct,vct,wct: ',uct(idjd),vct(idjd),wct(idjd)
	  print *,'ax,ay,az ',ax(idjd),ay(idjd),az(idjd)
	  print *,'bx,by,bz ',bx(idjd),by(idjd),bz(idjd)
         print *,'final u , v: ',u(idjd,k),v(idjd,k)
	endif
      enddo  ! k loop

!     below we interpolate quantities which may be affected by land-sea mask

!     set up land-sea mask from either tss or zss
      if(nemi.eq.2)then
        numneg=0
        do iq=1,ifull
         if(tss(iq).gt.0)then  ! over land
           land(iq)=.true.
         else                    ! over sea
           land(iq)=.false.
	    numneg=numneg+1
         endif    ! (tss(iq).gt.0) .. else ..
        enddo
        if(numneg.eq.0)nemi=1   ! should be using zss in that case
      endif  !  (nemi.eq.2)
      print *,'using nemi = ',nemi
      if(nemi.eq.1)then
        do iq=1,ifull
         if(zss(iq).gt.0)then  ! over land
           land(iq)=.true.
         else                    ! over sea
           land(iq)=.false.
         endif    ! (zss(iq).gt.0) .. else ..
        enddo
      endif  !  (nemi.eq.1)

      spval=999.
      do iq=1,ifull
         if(land(iq))then  ! over land
           tss_l(iq)=tss(iq)
           tss_s(iq)=spval
           sicedep(iq)=spval
         else                    ! over sea
	    numneg=numneg+1
           tss_s(iq)=abs(tss(iq))
           tss_l(iq)=spval
!          w(iq)=spval
!          w2(iq)=spval
!          ts(iq)=spval
!          ts(iq,2)=spval
           snowd(iq)=spval
           do k=1,ms
            tgg(iq,k)=spval
            wb(iq,k)=spval
           enddo
         endif    ! (tss(iq).gt.0) .. else ..
      enddo
      
      if(nmaxpr.eq.1)then
        print *,'before fill tss ',tss(idjd2)
        print *,'before fill tss_l, tss_s ',tss_l(idjd2),tss_s(idjd2)
        print *,'before fill/ints4 sicedep ',sicedep(idjd2)
        print *,'before fill wb'
        write(6,"('wb_s(1)#  ',9f7.3)") 
     .          ((wb(ii+(jj-1)*il,1),ii=id2-1,id2+1),jj=jd2-1,jd2+1)
      endif  ! (nmaxpr.eq.1)
      call fill_cc(tss_l,spval)
      call fill_cc(tss_s,spval)
      call fill_cc(snowd,spval)
      call fill_cc(sicedep,spval)
      do k=1,ms
        call fill_cc(tgg(1,k),spval)
        call fill_cc(wb(1,k),spval)
      enddo
      if(nmaxpr.eq.1)then
        print *,'after fill tss_l, tss_s ',tss_l(idjd2),tss_s(idjd2)
        print *,'after fill sicedep ',sicedep(idjd2)
        print *,'after fill wb'
        write(6,"('wb_s(1)#  ',9f7.3)") 
     .          ((wb(ii+(jj-1)*il,1),ii=id2-1,id2+1),jj=jd2-1,jd2+1)
        print *,'before ints4 psl(idjd2),zss(idjd2) ',
     .                        psl(idjd),zss(idjd2)
      endif  ! (nmaxpr.eq.1)

      if(nfly.gt.0)then
        norder=1
      else
        norder=nord
      endif
      call ints4(psl ,     nface4,xg4,yg4,norder)
      call ints4(zss ,nface4,xg4,yg4,norder)  
      if(nfly.eq.2)then
        call ints4(pmsl,   nface4,xg4,yg4,nord)
!       invert pmsl to get psl
        call to_psl(pmsl,psl,zss,t)  
      endif  ! (nfly.eq.2)
      call ints4(tss_l ,   nface4,xg4,yg4,nord)
      call ints4(tss_s ,   nface4,xg4,yg4,nord)
c     call ints4(precip,   nface4,xg4,yg4,nord)
      do k=1,ms
       call ints4(tgg(1,k),nface4,xg4,yg4,nord)
       call ints4(wb(1,k) ,nface4,xg4,yg4,nord)
      enddo
      if(nmaxpr.eq.1)then
        print *,'after ints4 idjd,zss(idjd) ',idjd,zss(idjd)
	 print *,'after ints4 psl,pmsl ',psl(idjd),pmsl(idjd)
        print *,'after ints4 wb_t'
        write(6,"('wb_t(1)#  ',9f7.3)") 
     .           ((wb(ii+(jj-1)*il,1),ii=id2-1,id2+1),jj=jd2-1,jd2+1)
      endif  ! (nmaxpr.eq.1)
c     call ints4(alb ,     nface4,xg4,yg4,nord)
c     call ints4(precc,    nface4,xg4,yg4,nord)
      if(nqg.ge.6)then
        call ints4(snowd,  nface4,xg4,yg4,nord)
        call ints4(sicedep,nface4,xg4,yg4,nord)
        print *,'after ints4 sicedep ',sicedep(idjd)
c       call ints4(cloudlo,nface4,xg4,yg4,nord)
c       call ints4(cloudmi,nface4,xg4,yg4,nord)
c       call ints4(cloudhi,nface4,xg4,yg4,nord)
      endif
      if(nqg.ge.7)then
c       call ints4(tscrn,  nface4,xg4,yg4,nord)

c       call ints4(qgscrn, nface4,xg4,yg4,nord)
c       call ints4(u10,    nface4,xg4,yg4,nord)
      endif

c     incorporate target land mask effects, e.g. into surface temperature
      do iq=1,ifull
       tss(iq)=tss_l(iq)
       if(land_t(iq))then
         sicedep(iq)=0.
	else
         tss(iq)=tss_s(iq)  ! no sign switch in CCAM
         snowd(iq)=0.
!        zs(iq)=-.1   ! dont do this
       endif
      enddo   ! iq loop
      if(nmaxpr.eq.1)then
        print *,'after ints tss_l, tss_s ',tss_l(idjd),tss_s(idjd)
        print *,'after ints tss',tss(idjd)
      endif  ! (nmaxpr.eq.1)

!     end of processing loop

!     restore target values into xx4, yy4
      rlong0x=rlong0_t  ! for indata cross-check
      rlat0x=rlat0_t
      schmidtx=schmidt_t
      do j=1,iquad
       do i=1,iquad
	 xx4(i,j)=xx4_sav(i,j)
	 yy4(i,j)=yy4_sav(i,j)
       enddo
      enddo

!     restore target zs and land arrays
      do iq=1,ifull
       zs(iq)=zs_t(iq)
       land(iq)=land_t(iq)  
      enddo

      rlong0=rlong0_t
      rlat0=rlat0_t
      schmidt=schmidt_t
      ds=ds_t
      return
      end

      subroutine ints4(s,nface4 ,xg4 ,yg4,nord)  ! does calls to intsb
      include 'newmpar.h'
      include 'parm.h'
      parameter (ntest=0)
      dimension s(ifull),nface4(ifull,4),xg4(ifull,4),yg4(ifull,4)
      real wrk(ifull,4)
      if(nord.eq.1)then
        do m=1,4
         call ints_blb(s,wrk(1,m),nface4(1,m),xg4(1,m),yg4(1,m))
        enddo
      else
        do m=1,4
         call intsb(s,wrk(1,m),nface4(1,m),xg4(1,m),yg4(1,m))
        enddo
      endif   ! (nord.eq.1)  .. else ..
      if(ntest.gt.0)then
        print *,'in ints4 for id,jd,nord: ',id,jd,nord
        print *,'nface4(1-4) ',(nface4(idjd,m),m=1,4)
        print *,'xg4(1-4) ',(xg4(idjd,m),m=1,4)
        print *,'yg4(1-4) ',(yg4(idjd,m),m=1,4)
        print *,'wrk(1-4) ',(wrk(idjd,m),m=1,4)
      endif
!     average 4 m values to central value
      do iq=1,ifull
       s(iq)=.25*(wrk(iq,1)+wrk(iq,2)+wrk(iq,3)+wrk(iq,4))
      enddo
      return
      end

      subroutine intsb(s,sout,nface,xg,yg)   ! N.B. sout here
!     same as subr ints, ut with sout passed back and no B-S      
!     s is input; sout is output array
c     later may wish to save idel etc between array calls
c     this one does linear interp in x on outer y sides
c     doing x-interpolation before y-interpolation
      include 'newmpar.h'
      include 'parm.h'
      dimension s(ifull),nface(ifull),xg(ifull),yg(ifull)
      common/work2b/sx(-1:il+2,-1:il+2,0:npanels)
     .       ,dum2(3*ifull -(il+4)*(il+4)*(npanels+1) )
c     real sx(-1:il+2,-1:il+2,0:npanels),
      real r(4),sout(ifull)
      include 'indices.h' ! in,is,iw,ie,inn,iss,iww,iee
      dimension in6(il,il,0:npanels),is6(il,il,0:npanels)
     .         ,iw6(il,il,0:npanels),ie6(il,il,0:npanels)
     .         ,inn6(il,il,0:npanels),iss6(il,il,0:npanels)
     .         ,iww6(il,il,0:npanels),iee6(il,il,0:npanels)
     .         ,ine6(il,il,0:npanels),ise6(il,il,0:npanels)
     .         ,ien6(il,il,0:npanels),iwn6(il,il,0:npanels)
      equivalence (in,in6),(is,is6),(iw,iw6),(ie,ie6)
      equivalence (inn,inn6),(iss,iss6),(iww,iww6),(iee,iee6)
      equivalence (ine,ine6),(ise,ise6),(ien,ien6),(iwn,iwn6)
      ind(i,j,n)=i+(j-1)*il+n*il*il  ! *** for n=0,npanels
c     this is intsb           EW interps done first
c     first extend s arrays into sx - this one -1:il+2 & -1:il+2
      do n=0,npanels
       do j=1,il
        do i=1,il
         sx(i,j,n)=s(ind(i,j,n))
        enddo  ! i loop
        sx(0,j,n)=s(iw6(1,j,n))
        sx(-1,j,n)=s(iww6(1,j,n))
        sx(il+1,j,n)=s(ie6(il,j,n))
        sx(il+2,j,n)=s(iee6(il,j,n))
       enddo   ! j loop
       do i=1,il
        sx(i,0,n)=s(is6(i,1,n))
        sx(i,-1,n)=s(iss6(i,1,n))
        sx(i,il+1,n)=s(in6(i,il,n))
        sx(i,il+2,n)=s(inn6(i,il,n))
       enddo  ! i loop
c      for ew interpolation, sometimes need (different from ns):
c          (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
c         (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)
       sx(-1,0,n)=s(lwws(n))
       sx(0,0,n)=s(lws(n))
       sx(0,-1,n)=s(lwss(n))
       sx(il+1,0,n)=s(les(n))
       sx(il+2,0,n)=s(lees(n))
       sx(il+1,-1,n)=s(less(n))
       sx(-1,il+1,n)=s(lwwn(n))
       sx(0,il+2,n)=s(lwnn(n))
       sx(il+2,il+1,n)=s(leen(n))
       sx(il+1,il+2,n)=s(lenn(n))
       sx(0,il+1,n)   =s(iwn6(1,il,n))
       sx(il+1,il+1,n)=s(ien6(il,il,n))
      enddo    ! n loop

      do iq=1,ifull
c       if(iq.eq.idjd)print *,'iq,nface,xg,yg ',
c    .                         iq,nface(iq),xg(iq),yg(iq)
        n=nface(iq)
        idel=int(xg(iq))
        xxg=xg(iq)-idel
c       yg here goes from .5 to il +.5
        jdel=int(yg(iq))
        yyg=yg(iq)-jdel
c       if(iq.eq.idjd)then
c         print *,'iq,idel,xxg,jdel,yyg,n ',
c    .             iq,idel,xxg,jdel,yyg,n
c         print *,'sx nn=1',sx(idel  ,jdel-1,n),sx(idel+1,jdel-1,n)
c         print *,'sx nn=2',sx(idel-1,jdel  ,n),sx(idel  ,jdel  ,n)
c    .                     ,sx(idel+1,jdel  ,n),sx(idel+2,jdel  ,n)
c         print *,'sx nn=3',sx(idel-1,jdel+1,n),sx(idel  ,jdel+1,n)
c    .                     ,sx(idel+1,jdel+1,n),sx(idel+2,jdel+1,n)
c         print *,'sx nn=4',sx(idel  ,jdel+2,n),sx(idel+1,jdel+2,n)
c       endif
        do nn=2,3       ! N.B.
         c1=sx(idel-1,jdel+nn-2,n)
         c2=sx(idel  ,jdel+nn-2,n)
         c3=sx(idel+1,jdel+nn-2,n)
         c4=sx(idel+2,jdel+nn-2,n)
         r(nn)=((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.)
     .     -xxg*(1.+xxg)*c4/3.)
     .     +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
        enddo    ! nn loop
c       r       ={(1-x     )*{(2-x     )*[(1+x     )*c2-x     *c1/3]
c         -x     *(1+x     )*c4/3}
c         +x    *(1+x     )*(2-x     )*c3}/2
        do nn=1,4,3       ! N.B.
         c2=sx(idel  ,jdel+nn-2,n)
         c3=sx(idel+1,jdel+nn-2,n)
         r(nn)=(1.-xxg)*c2 +xxg*c3
        enddo    ! nn loop
c       array(iq)=((1.-yyg)*((2.-yyg)*((1.+yyg)*r(2)-yyg*r(1)/3.)
c    .             -yyg*(1.+yyg)*r(4)/3.)
c    .             +yyg*(1.+yyg)*(2.-yyg)*r(3))/2.
c      following does Bermejo Staniforth
        aaa=((1.-yyg)*((2.-yyg)*((1.+yyg)*r(2)-yyg*r(1)/3.)
     .      -yyg*(1.+yyg)*r(4)/3.)
     .      +yyg*(1.+yyg)*(2.-yyg)*r(3))/2.
c       if(iq.eq.idjd)print *,'before BS aaa sx sx+ sx0+ sx++',aaa,
c    .                         sx(idel,jdel,n),sx(idel+1,jdel,n),
c    .                         sx(idel,jdel+1,n),sx(idel+1,jdel+1,n)
        aaa=min( aaa , max( sx(idel,jdel,n),sx(idel+1,jdel,n),
     .                      sx(idel,jdel+1,n),sx(idel+1,jdel+1,n) ) )
        sout(iq)=max( aaa , min( sx(idel,jdel,n),sx(idel+1,jdel,n),
     .                        sx(idel,jdel+1,n),sx(idel+1,jdel+1,n) ) )
      enddo    ! iq loop
      return
      end

      subroutine ints_blb(s,sout,nface,xg,yg) 
c     this one does bi-linear interpolation only
      include 'newmpar.h'
      include 'parm.h'
      dimension s(ifull),nface(ifull),xg(ifull),yg(ifull)
      common/work2b/sx(-1:il+2,-1:il+2,0:npanels)
     .       ,dum2(3*ifull -(il+4)*(il+4)*(npanels+1) )
c     real sx(-1:il+2,-1:il+2,0:npanels),
      real sout(ifull)
      include 'indices.h' ! in,is,iw,ie,inn,iss,iww,iee
      dimension in6(il,il,0:npanels),is6(il,il,0:npanels)
     .         ,iw6(il,il,0:npanels),ie6(il,il,0:npanels)
     .         ,inn6(il,il,0:npanels),iss6(il,il,0:npanels)
     .         ,iww6(il,il,0:npanels),iee6(il,il,0:npanels)
     .         ,ine6(il,il,0:npanels),ise6(il,il,0:npanels)
     .         ,ien6(il,il,0:npanels),iwn6(il,il,0:npanels)
      equivalence (in,in6),(is,is6),(iw,iw6),(ie,ie6)
      equivalence (inn,inn6),(iss,iss6),(iww,iww6),(iee,iee6)
      equivalence (ine,ine6),(ise,ise6),(ien,ien6),(iwn,iwn6)
      ind(i,j,n)=i+(j-1)*il+n*il*il  ! *** for n=0,npanels
c     first extend s arrays into sx - this one -1:il+2 & -1:il+2
c                    but for bi-linear only need 0:il+1 &  0:il+1
      do n=0,npanels
       do j=1,il
        do i=1,il
         sx(i,j,n)=s(ind(i,j,n))
        enddo  ! i loop
       enddo   ! j loop
       do j=1,il
        sx(0,j,n)=s(iw6(1,j,n))
        sx(il+1,j,n)=s(ie6(il,j,n))
       enddo   ! j loop
       do i=1,il
        sx(i,0,n)=s(is6(i,1,n))
        sx(i,il+1,n)=s(in6(i,il,n))
       enddo  ! i loop

       sx(0,0,n)=s(lws(n))
       sx(il+1,0,n)=s(les(n))
       sx(0,il+1,n)   =s(iwn6(1,il,n))
       sx(il+1,il+1,n)=s(ien6(il,il,n))
      enddo    ! n loop

      do iq=1,ifull
       n=nface(iq)
       idel=int(xg(iq))
       xxg=xg(iq)-idel
       jdel=int(yg(iq))
       yyg=yg(iq)-jdel
c      print *,'iq,idel,jdel,n ',iq,idel,jdel,n
       sout(iq)=yyg*(xxg*sx(idel+1,jdel+1,n)
     .               +(1.-xxg)*sx(idel,jdel+1,n))
     .    +(1.-yyg)*(xxg*sx(idel+1,jdel,n)
     .               +(1.-xxg)*sx(idel,jdel,n))
      enddo    ! iq loop
      return
      end

      subroutine fill_cc(a,value)
c     routine fills in interior of an array which has undefined points
      include 'newmpar.h'
      include 'indices.h'
      real a(ifull)         ! input and output array
      real value            ! array value denoting undefined
      real b(ifull)
      
      num=0
2     nrem=0
      num=num+1
      do 6 iq=1,ifull
      b(iq)=a(iq)
      if(a(iq).eq.value)then
        neighb=0
        av=0.
        if(a(in(iq)).ne.value)then
          neighb=neighb+1
          av=av+a(in(iq))
        endif
        if(a(ie(iq)).ne.value)then
          neighb=neighb+1
          av=av+a(ie(iq))
        endif
        if(a(iw(iq)).ne.value)then
          neighb=neighb+1
          av=av+a(iw(iq))
        endif
        if(a(is(iq)).ne.value)then
          neighb=neighb+1
          av=av+a(is(iq))
        endif
        if(neighb.gt.0)then
          b(iq)=av/neighb
	   avx=av
        else
          nrem=nrem+1    ! current number of points without a neighbour
        endif
      endif
6     continue
      do iq=1,ifull
        a(iq)=b(iq)
      enddo
      if(nrem.gt.0)then
        if(num.le.2)
     .    print *,'after 1/2 time thru fill loop num,nrem,avx = ',
     .                                           num,nrem,avx
        go to 2
      endif  ! (nrem.gt.0)
      return
      end

