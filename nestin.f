      subroutine nestin                ! for globpea
      include 'newmpar.h'
c     ik,jk,kk are array dimensions read in infile - not for globpea
c     int2d code - not used for globpea
      include 'aalat.h'
      include 'arrays.h'
      include 'const_phys.h'
      include 'dates.h'    ! mtimer
      include 'dava.h'
      include 'davb.h'     ! psls,qgg,tt,uu,vv
      include 'map.h'
      include 'parm.h'     ! qgmin
      include 'pbl.h'      ! sice
      include 'soil.h'     ! tss
      include 'soilsnow.h' ! tgg
      include 'stime.h'    ! kdate_s,ktime_s  sought values for data read
      common/nest/ta(ifull,kl),ua(ifull,kl),va(ifull,kl),psla(ifull),
     .            tb(ifull,kl),ub(ifull,kl),vb(ifull,kl),pslb(ifull),
     .            qa(ifull,kl),qb(ifull,kl),tssa(ifull),tssb(ifull)
      common/schmidtx/rlong0x,rlat0x,schmidtx ! infile, newin, nestin, indata
      common/sigin/sigin(kl),kk  ! for vertint, infile
      common/work2/zsb(ifull),dum1(ifull),tssi(ifull),dum2(ifull),
     .  dum3(ifull),dum4(ifull),dum5(ifull),dum6(ifull),dum7(ifull),
     .                        dumxx(ifull,18-9)
      common/work3/dum3a(ifull,kl),dum3b(ifull,kl,4)
      integer num,mtimea,mtimeb
      data num/0/,mtimea/0/,mtimeb/-1/
      save num,mtimea,mtimeb
c     mtimer, mtimeb are in minutes
      if(ktau.lt.100)print *,'in nestin ktau,mtimer,mtimea,mtimeb ',
     .                                  ktau,mtimer,mtimea,mtimeb
      if(ktau.lt.100)print *,'with kdate_s,ktime_s >= ',kdate_s,ktime_s
      if(mtimeb.eq.-1)then
        mtimeb=mtimer  ! zero in fact
        print *,'set nesting fields to those already read in via indata'
        do iq=1,ifull
         pslb(iq)=psl(iq)
         tssb(iq)=tss(iq)
        enddo
        tb(:,:)=t(:,:)
        qb(:,:)=max(qg(:,:),qgmin)
        ub(:,:)=u(:,:)
        vb(:,:)=v(:,:)
        return
      endif       ! (mtimeb.eq.-1)

      if(mtimer.gt.0.and.mtimer.le.mtimeb)go to 6  ! allows for dt<1 minute
c     transfer mtimeb fields to mtimea
      mtimea=mtimeb
      do iq=1,ifull
       psla(iq)=pslb(iq)
       tssa(iq)=tssb(iq)
      enddo
      ta(:,:)=tb(:,:)
      qa(:,:)=qb(:,:)
      ua(:,:)=ub(:,:)
      va(:,:)=vb(:,:)

c     read tb etc  - for globpea, straight into tb etc
      if(io_in.eq.1.or.io_in.eq.3)then
      call infile(meso2,kdate_r,ktime_r,nem2, ! different from DARLAM
     .  timeg_b,ds_r,pslb,dum1,zsb,dum2,dum3, ! timeg_r,ds_r,psl,pmsl,zs,em,f
     .  tssb,dum4,dum3a,dum3a,dum5,dum6,dum7, ! tss,precip,wb,wbice,alb,snowd,sicedep
     .  tb,ub,vb,qb,dum3a,                    ! t,u,v,qg,tgg
     .	 dum3a,dum3a,dum3a, dum2,dum3,dum4,dum5,1)
!                above are:   tggsn,smass,ssdn, ssdnn,osnowd,snage,isflag
      endif   ! (io_in.eq.1.or.io_in.eq.3)

      if(io_in.eq.-1.or.io_in.eq.-3)then
          call onthefly(kdate_r,ktime_r,
     .     pslb,zsb,tssb,dum3a,dum3a,dum6,dum7,
     .     tb,ub,vb,qb,dum3a,                    ! t,u,v,qg,tgg
     .	    dum3a,dum3a,dum3a, dum2,dum3,dum4,dum5,1)
!                above are:   tggsn,smass,ssdn, ssdnn,osnowd,snage,isflag
      endif   ! (io_in.eq.1.or.io_in.eq.3)
      write (6,"('zsb# nestin  ',9f7.1)") 
     .          ((zsb(ii+(jj-1)*il),ii=id-1,id+1),jj=jd-1,jd+1)
      write (6,"('tssb# nestin ',9f7.1)") 
     .          ((tssb(ii+(jj-1)*il),ii=id-1,id+1),jj=jd-1,jd+1)
   
      if(abs(rlong0  -rlong0x).gt..01.or.
     .   abs(rlat0    -rlat0x).gt..01.or.
     .   abs(schmidt-schmidtx).gt..01)stop "grid mismatch in infile"

      kdhour=(ktime_r-ktime)/100                      ! integer hour diff
      kdmin=(ktime_r-100*(ktime_r/100))-(ktime-100*(ktime/100))
      print *,'nesting file has: kdate_r,ktime_r,kdhour,kdmin ',
     .                           kdate_r,ktime_r,kdhour,kdmin
      mtimeb=60*24*(iabsdate(kdate_r,kdate)-iabsdate(kdate,kdate))
     .               +60*kdhour+kdmin
c     mtimeb=60*24*(iabsdate(kdate_r,kdate)-iabsdate(kdate,kdate))
c    .               +nint(60*(ktime_r-ktime)/100.)  ! up till 29/11/02
      print *,'kdate_r,iabsdate ',kdate_r,iabsdate(kdate_r,kdate)
      print *,'giving mtimeb = ',mtimeb

c     print additional information
      print *,' kdate ',kdate,' ktime ',ktime
      print *,'timeg,mtimer,mtimea,mtimeb: ',
     .         timeg,mtimer,mtimea,mtimeb
      print *,'ds,ds_r ',ds,ds_r

!     if(mod(ktau,nmaxpr).eq.0.or.diag)then
c       following is useful if troublesome data is read in
        print *,'following max/min values printed from nestin'
        call maxmin(ub,'ub',ktau,1.,kl)
        call maxmin(vb,'vb',ktau,1.,kl)
        call maxmin(tb,'tb',ktau,1.,kl)
        call maxmin(qb,'qb',ktau,1.e3,kl)
        print *,'following are really psl not ps'
        call maxmin(pslb,'ps',ktau,100.,1)
!     endif

      qb(:,:)=max(qb(:,:),qgmin)

      if(kk.lt.kl)then
c       this section allows for different number of vertical levels
c       presently assume sigin (up to kk) levels are set up as per nsig=6
c       option in eigenv, though original csiro9 levels are sufficiently
c       close for these interpolation purposes.
        if(diag)then
          print *,'kk,sigin ',kk,(sigin(k),k=1,kk)
          print *,'tb before vertint ',(tb(idjd,k),k=1,kk)
        endif
        call vertint(tb,1)  ! transforms tb from kk to kl
        if(diag)then
          print *,'tb after vertint ',(tb(idjd,k),k=1,kk)
          print *,'qb before vertint ',(qb(idjd,k),k=1,kk)
        endif
        call vertint(qb,2)
        if(diag)print *,'qb after vertint ',(qb(idjd,k),k=1,kk)
        call vertint(ub,3)
        call vertint(vb,4)
      endif  ! (kk.lt.kl)

c     N.B. tssb (sea) only altered for newtop=2 (done here now)
      if(newtop.eq.2)then
!       reduce sea tss to mslp      e.g. for QCCA in NCEP GCM
        do iq=1,ifull
         if(tssb(iq).lt.0.)tssb(iq)=
     .                       tssb(iq)-zsb(iq)*stdlapse/grav  ! N.B. -
        enddo
      endif  ! (newtop.eq.2)

      do iq=1,ifull
        tssb(iq) = abs(tssb(iq))
      enddo
!     test code in nestin for modifying SSTs
c	do j=270,288
c	 do i=9,22
c	  tssb(i,j)=tssb(i,j)-.5
c	 enddo
c	enddo      

      if(newtop.ge.1)then
c       in these cases redefine pslb, tb and (effectively) zsb using zs
c       this keeps inner-mesh land mask & zs
c       presently simplest to do whole pslb, tb (& qb) arrays
        if(nmaxpr.eq.1)then
          print *,'zs (idjd) :',zs(idjd)
          print *,'zsb (idjd) :',zsb(idjd)
          print *,'psl (idjd) :',psl(idjd)
          print *,'pslb in(idjd) :',pslb(idjd)
          print *,'now call retopo from nestin'
        endif
        call retopo(pslb,zsb,zs,tb,qb)
        if(nmaxpr.eq.1)then
          print *,'pslb out(idjd) :',pslb(idjd)
          print *,'after pslb print; num= ',num
        endif
      endif   !  newtop.ge.1

      if(num.eq.0)then
        num=1
        call printa('zs  ',zs        ,ktau,0  ,ia,ib,ja,jb,0.,.01)
        call printa('zsb ',zsb       ,ktau,0  ,ia,ib,ja,jb,0.,.01)
        call printa('psl ',psl       ,ktau,0  ,ia,ib,ja,jb,0.,1.e2)
        call printa('pslb',pslb      ,ktau,0  ,ia,ib,ja,jb,0.,1.e2)
        call printa('t   ',t(1,nlv),ktau,nlv,ia,ib,ja,jb,200.,1.)
        call printa('tb  ',tb(1,nlv),ktau,nlv,ia,ib,ja,jb,200.,1.)
        call printa('u   ',u(1,nlv),ktau,nlv,ia,ib,ja,jb,0.,1.)
        call printa('ub  ',ub(1,nlv),ktau,nlv,ia,ib,ja,jb,0.,1.)
        call printa('v   ',v(1,nlv),ktau,nlv,ia,ib,ja,jb,0.,1.)
        call printa('vb  ',vb(1,nlv),ktau,nlv,ia,ib,ja,jb,0.,1.)
        call printa('davt',davt,0,0,ia,ib,ja,jb,0.,10.)
        return
      endif   !  num.eq.0

c     now use tt, uu, vv arrays for time interpolated values
6     timerm=ktau*dt/60.   ! real value in minutes (in case dt < 60 seconds)
      cona=(mtimeb-timerm)/real(mtimeb-mtimea)
      conb=(timerm-mtimea)/real(mtimeb-mtimea)
      psls(:)=cona*psla(:)+conb*pslb(:)
      tt (:,:)=cona*ta(:,:)+conb*tb(:,:)
      qgg(:,:)=cona*qa(:,:)+conb*qb(:,:)
      uu (:,:)=cona*ua(:,:)+conb*ub(:,:)
      vv (:,:)=cona*va(:,:)+conb*vb(:,:)

c     calculate time interpolated tss (into tssi)
      if(namip.eq.0)return     ! namip SSTs/sea-ice take precedence
      tssi(:)=cona*tssa(:)+conb*tssb(:)  
      do iq=1,ifull
       if(.not.land(iq))then
	  if(tssi(iq).gt.273.2)then
	    tss(iq)=tssi(iq)
	    tgg(iq,1)=tssi(iq)
!          all others are for implied sea ice points
!          N.B. fracice, sice etc are updated once daily in sflux
	  elseif(tss(iq).gt.273.2)then  
	    tss(iq)=tssi(iq)
	    tgg(iq,1)=tssi(iq)
!          if already sea ice implied, use present tss, tgg1
         endif	    
       endif  ! (.not.land(iq))
      enddo  
	      
      return
      end
