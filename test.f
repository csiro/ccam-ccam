      subroutine convjlm      ! jlm convective scheme - Version 2
!     with nkuo<0 does large-scale first
!     has fldownn depending on delta sigma      
      include 'newmpar.h'
      parameter (ntest=0)      ! 1 or 2 to turn on
!     parameter (iterconv=1)   ! to kuocom.h
!     parameter (fldown=.6)    ! to kuocom.h
!     parameter (detrain=.05)  ! to kuocom.h
!     parameter (alfcon=.015)  ! e.g. .01 or .015   or 1.1 (over sea)
!     parameter (alflnd=1.1)   ! e.g. 1.1  with nqbot=3
!     parameter (alfsea=1.1)   ! e.g. 1.1  with nqbot=3
!     parameter (nqbot=2)      ! 0: original qq(1), 1: qs1, 2: qbas for land/sea
!                                3:for RH-enhanced base
      parameter (nalfs=0)      ! 1: usual, 0: off, -1: for omgf effect
      parameter (nfluxq=1)     ! 1 from 6/8/02
      parameter (nfluxdsk=2)   ! 1 till 6/8/02
      parameter (nbase=1)      ! 0 shared below, 1 uses only one (originally 0)
!     parameter (methdetr=1)   ! 1 (top only); 2 (top half); 4 top quarter (= kbconv)
!     parameter (methprec=2)   ! 1 (top only); 2 (top half); 4 top quarter (= kbconv)
      parameter (no2layer=0)   ! usually 0, 1 to suppress 2-layer clouds (even wetter!)
      parameter (nuvconv=0)    ! usually 0, 1 to turn on momentum mixing
      parameter (kcl_top=kl-2) !max level for cloud top (convjlm,radrive,vertmix)
!     parameter (dsig2=.1, dsig4=.7)  ! used with detrainx
c     nevapls, nevapcc:  turn off/on evap of ls or cc --- now through parm.h
c        0 off, 1 for Hal's evap, 2 for jlm, 3 for UK (ls only), 4 & 5 newer UK
      parameter (ars=461.,grav=9.806,hl=2.5104e6,hlars=hl/ars)
      include 'arrays.h'
      include 'constant.h'
      include 'dava.h'    ! davt
      include 'kuocom.h'  ! kbsav,ktsav,convfact,convpsav,ndavconv
      include 'map.h'     ! land
      include 'morepbl.h'
      include 'nlin.h'
      include 'parm.h'
      include 'prec.h'
      include 'sigs.h'
      include 'soil.h'
      include 'tracers.h'  ! ngas, nllp, ntrac
      include 'vvel.h'
      include 'xarrs.h'
      common/work2/ktmax(ifull),conrev(ifull),fluxr(ifull)
     .  ,fluxb(ifull),fluxt(ifull),kbsav_b(ifull),sumd(ifull)
     .  ,rnrt(ifull),rnrtc(ifull),sumqs(ifull),sumqt(ifull)
     .  ,sumss(ifull),factdav(ifull),kbsav_ls(ifull)
     .  ,qbase(ifull),fldow(ifull),qxcess(ifull),fluxq(ifull)
      real sbase(ifull),alfqarr(ifull),alfsarr(ifull)
      common/work3/tt(ifull,kl),qq(ifull,kl),qs(ifull,kl)
     .  ,dqsdt(ifull,kl)
      common/work3b/s(ifull,kl),hs(ifull,kl)
      common/work3c/rnrt3(ifull,kl)
      common/work3f/delq(ifull,kl),dels(ifull,kl),rnrt3d(ifull,kl)
      real dsk(kl),phi(ifull,kl),delu(ifull,kl),delv(ifull,kl)
      real revc(ifull,kl)  ! just for evap diags
      real alfss(kl),alfm(kl),qbass(kl)
      real alfqq_l(kl),alfqq_s(kl),alfqq1_l(kl),alfqq1_s(kl)
      equivalence (phi,s),(delu,revc,delq),(delv,dels)
!     data rhcv/.75/                  now in kuocom

      common /es_table/ table(0:220)
c     arithmetic statement functions to replace call to establ.
c     t is temp in kelvin, which should lie between 123.16 and 343.16;
c     tdiff is difference between t and 123.16, subject to 0 <= tdiff <= 220
      tdiff(tm)=min(max(tm-123.16 , 0.) , 220.)
      establ(tm) =(1.-(tdiff(tm)-aint(tdiff(tm))))*table(int(tdiff(tm)))
     &           + (tdiff(tm)-aint(tdiff(tm)))*table(int(tdiff(tm))+1)
      Aev(tm) = 2.008e-9*tm**2 - 1.385e-6*tm + 2.424e-4  !For UKMO evap scheme
      Asb(tm) = max (-5.2e-9*tm**2+2.5332e-6*tm-2.9111e-4,1.e-5) !For UKMO subl

      hlcp=hl/cp
      roncp=r/cp
      do k=1,kl
       dsk(k)=-dsig(k)    !   dsk = delta sigma (positive)
       alfm(k)=min(.5,alfcon/(1.e-10+sig(1)-sig(k)))  ! e.g. alfcon=.01
       alfss(k)=1.-alfm(k)
	if(nalfs.eq.0)alfss(k)=1.  ! for no alf effect
      enddo     ! k loop
      nqbot=3
      if(alflnd.lt.1.)nqbot=2
      if(nqbot.eq.2)then
        do k=1,kl
         alfqq_s(k)=1.-alfsea*sig(k)**8   ! e.g. alfsea=.5; =0. for off
         alfqq_l(k)=1.-alflnd*sig(k)**8   ! e.g. alfsea=.4; =0. for off
         alfqq1_s(k)=1.-alfqq_s(k)
         alfqq1_l(k)=1.-alfqq_l(k)
        enddo     ! k loop
      endif   ! (nqbot.eq.2)
      if(nqbot.eq.3)then
        do k=1,kl
         alfqq_s(k)=alfsea
         alfqq_l(k)=alflnd
         alfqq1_s(k)=0.
         alfqq1_l(k)=0.
        enddo     ! k loop
      endif   ! (nqbot.eq.3)
  
      ecmwa1=5.44e-4
      ecmwa2=5.09
      ecmwcc=sqrt(.05)
      ecmwls=1.

c     convective then L/S rainfall
      do iq=1,ifull
       conrev(iq)=1000.*ps(iq)/(grav*dt) ! factor to conv. precip to g/m2/s
       rnrt(iq)=0.     ! initialize large-scale rainfall array
       ktmax(iq)=kl      ! preset 1 level above current topmost-permitted ktsav
       kbsav_ls(iq)=0    ! for L/S
      enddo  ! iq loop

      tt=t       !3D
      qq=qg      !3D
      rnrt3d=0.  !3D initialize vertically accum. 3d convective rainfall 
      if(nkuo.lt.0)go to 5

!__________________________beginning of convective calculations_____________________

2     do itn=1,iterconv
!     calculate geopotential height 
      do iq=1,ifull
       kbsav(iq)=kl     ! preset value for no convection
       kbsav_b(iq)=kl   ! preset value for no convection
       ktsav(iq)=kl     ! preset value for no convection
       convpsav(iq)=0.
       phi(iq,1)=bet(1)*tt(iq,1)
      enddo     ! iq loop

      do k=2,kl
       do iq=1,ifull
        phi(iq,k)=phi(iq,k-1)+bet(k)*tt(iq,k)+betm(k)*tt(iq,k-1)
       enddo     ! iq loop
      enddo      ! k  loop

      do k=1,kl   
         do iq=1,ifull
          dels(iq,k)=0.
          delq(iq,k)=0.
          rnrt3(iq,k)=0.  ! initialize iter. 3d convective rainfall array
          es=establ(tt(iq,k))
          pk=ps(iq)*sig(k)
          qs(iq,k)=max(.622*es/(pk-es),1.5e-6)  
          dqsdt(iq,k)=qs(iq,k)*pk*hlars/(tt(iq,k)**2*(pk-es))
          s(iq,k)=cp*tt(iq,k)+phi(iq,k)  ! dry static energy
c         calculate hs
          hs(iq,k)=s(iq,k)+hl*qs(iq,k) ! saturated moist static energy
         enddo  ! iq loop
      enddo     ! k loop

      if(ktau.eq.1)then
       print *,'itn,iterconv,nbase,no2layer,nuvconv ',
     .          itn,iterconv,nbase,no2layer,nuvconv
       print *,'methdetr,detrain,detrainx,dsig2,dsig4 ',
     .          methdetr,detrain,detrainx,dsig2,dsig4
       print *,'methprec,fldown ',
     .          methprec,fldown
       print *,'nfluxq,nfluxdsk,nalfs ',
     .          nfluxq,nfluxdsk,nalfs
       print *,'nqbot,alfcon,alflnd,alfsea ',nqbot,alfcon,alflnd,alfsea
       write (6,"('alfqq_l ',19f7.3/(8x,19f7.3))") alfqq_l
       write (6,"('alfqq_s ',19f7.3/(8x,19f7.3))") alfqq_s
       write (6,"('alfss   ',19f7.3/(8x,19f7.3))") alfss
      endif
      if(ntest.ne.0.or.diag)then
        print *,'near beginning of convjlm loop; ktau,itn,id,jd: ',
     .                                           ktau,itn,id,jd
        write (6,"('rh   ',19f7.2/(8x,19f7.2))") 
     .             (100.*qq(idjd,k)/qs(idjd,k),k=1,kl)
        write (6,"('qs   ',19f7.3/(8x,19f7.3))") 
     .             (1000.*qs(idjd,k),k=1,kl)
        write (6,"('qq   ',19f7.3/(8x,19f7.3))") 
     .             (1000.*qq(idjd,k),k=1,kl)

        pwater0=0.   ! in mm     
	 do k=1,kl
	  iq=idjd
         pwater0=pwater0-dsig(k)*qg(idjd,k)*ps(idjd)/g
	  if(land(iq))then
c*	    qbass(k)=alfqq_l(k)*qq(iq,k)+(1.-alfqq_l(k))*qq(iq,1) 
	    qbass(k)=alfqq_l(k)*qq(iq,k)+alfqq1_l(k)*qq(iq,1) 
	  else
c*	    qbass(k)=alfqq_s(k)*qq(iq,k)+(1.-alfqq_s(k))*qq(iq,1) 
	    qbass(k)=alfqq_s(k)*qq(iq,k)+alfqq1_s(k)*qq(iq,1) 
	  endif
	 enddo
        write (6,"('qbas ',19f7.3/(8x,19f7.3))") 
     .             (1000.*qbass(k),k=1,kl)
        write (6,"('tt   ',19f7.2/(8x,19f7.2))") 
     .             (tt(idjd,k),k=1,kl)
        write (6,"('s/cp ',19f7.2/(8x,19f7.2))") (s(idjd,k)/cp,k=1,kl)
        write (6,"('h/cp ',19f7.2/(8x,19f7.2))") 
     .             (s(idjd,k)/cp+hlcp*qq(idjd,k),k=1,kl)
        write (6,"('hb/cp',19f7.2/(8x,19f7.2))") 
     .             (s(idjd,k)/cp+hlcp*qbass(k),k=1,kl)
        write (6,"('hs/cp',19f7.2/(8x,19f7.2))") 
     .             (hs(idjd,k)/cp,k=1,kl)
        write (6,"('k  ',19i7/(3x,19i7))") (k,k=1,kl)
      endif

      do k=kl/2,kuocb+1,-1   ! downwards to find cloud base and kbsav_b
         do iq=1,ifull
c         find tentative cloud base, and bottom of below-cloud layer
          if(qq(iq,k-1).gt.rhcv*qs(iq,k-1))then
!           next line ensures sub-cloud layer moist and contiguous 
            if(kbsav_b(iq).eq.k.and.sig(kbsav(iq)).gt..8)kbsav_b(iq)=k-1 
	     if(land(iq))then
		qbas=alfqq_l(k-1)*qq(iq,k-1)+alfqq1_l(k-1)*qq(iq,1) 
              alfq=alfqq_l(k-1)   
            else
		qbas=alfqq_s(k-1)*qq(iq,k-1)+alfqq1_s(k-1)*qq(iq,1) 
              alfq=alfqq_s(k-1)   
	     endif
	     alfs=alfss(k-1)  ! simple one
!           following  gives alfs=1. for omgf +ve	     
            sbas=alfs*s(iq,k-1)+(1.-alfs)*s(iq,1) ! e.g. alfs=1.
            if(sbas+hl*qbas.gt.hs(iq,k))then 
              if(qbas.gt.max(qs(iq,k),qq(iq,k)))then  ! newer test 
                kbsav(iq)=k-1
	         qbase(iq)=qbas
	         sbase(iq)=sbas
                kbsav_b(iq)=k-1
                ktsav(iq)=k
                alfqarr(iq)=alfq    ! for use in flux calc
                alfsarr(iq)=alfs    ! for use in flux calc
              endif
            endif
c	     if(ntest.eq.1.and.iq.eq.idjd)then   ! needs sopt
c	       print *,'k,(sbas+hl*qbas)/cp,hs(iq,k)/cp ',
c     .	         k,(sbas+hl*qbas)/cp,hs(iq,k)/cp
c		print *,'kbsav,qbas,qs(.,k),qq(.,k) ',
c     .	         kbsav(iq),qbas,qs(iq,k),qq(iq,k)
c	     endif
          endif
         enddo  ! iq loop
      enddo     ! k loop

      do k=kuocb+2,kl   ! upwards to find cloud top
       do iq=1,ifull
        if(ktsav(iq).eq.k-1.and.k.lt.ktmax(iq))then ! ktmax allows for itn
          hbase=sbase(iq)+hl*qbase(iq)
          if(hbase.gt.hs(iq,k))ktsav(iq)=k
        endif
       enddo    ! iq loop
      enddo     ! k loop

      if(no2layer.eq.1)then
        do iq=1,ifull
         if(ktsav(iq).eq.kbsav(iq)+1)then
            ktsav(iq)=kl
         endif
        enddo
      endif

      if(ntest.ne.0.or.diag)then
       iq=idjd
       iq=idjd
       kdown=min(kl,kbsav(iq)+nint(.6+.75*(ktsav(iq)-kbsav(iq))))
       print *,"ktau,kbsav_b,kbsav,ktsav,kdown ",
     .    ktau,kbsav_b(idjd),kbsav(idjd),ktsav(idjd),kdown
       print *,'alfqarr,alfsarr ',alfqarr(idjd),alfsarr(idjd)
      endif

!     calculate dels and delq for all layers for unit base mass flux
!     designed to be OK if q(k)>qs(k)
      do iq=1,ifull
!      dels and delq for top layer
!      following line allows for qq>qs, to avoid drying layer
       qsk=max(qs(iq,ktsav(iq)),qq(iq,ktsav(iq)))  
       qprec=max(0.,qbase(iq)-qsk)               ! tempry for Himalaya problem
       dels(iq,ktsav(iq))=dels(iq,ktsav(iq))+hl*qprec        ! precip. heating
       dels(iq,ktsav(iq))=dels(iq,ktsav(iq))+sbase(iq)       ! s flux
       delq(iq,ktsav(iq))=delq(iq,ktsav(iq))+qsk         ! moistening done later

       kdown=min(kl , kbsav(iq)+nint(.6+.75*(ktsav(iq)-kbsav(iq))))
       qsk=min(qs(iq,kdown),qq(iq,kdown))    ! N.B. min for downdraft subs
       tdown=t(iq,kbsav(iq))+(s(iq,kdown)-s(iq,kbsav(iq))
     .      +hl*(qsk-qs(iq,kbsav(iq))))/(cp+hl*dqsdt(iq,kbsav(iq)))
!      don't want moistening of cloud base layer (can lead to -ve fluxes)     
c      qdown=min( qq(iq,kbsav(iq)) , qs(iq,kbsav(iq))+
c    .           (tdown-t(iq,kbsav(iq)))*dqsdt(iq,kbsav(iq)) )
       qdown= qs(iq,kbsav(iq))+
     .           (tdown-t(iq,kbsav(iq)))*dqsdt(iq,kbsav(iq)) 
       dprec=qdown-qsk                      ! to be mult by fldownn
!      fldownn=fldown   ! up till 27/11/01 fldown=.2 say	
       fldownn=fldown*(sig(kbsav(iq))-sig(ktsav(iq))) ! fldown=.66 say	
	totprec=qprec-fldownn*dprec 
	fldow(iq)=(.5+sign(.5,totprec))*fldownn ! suppr. downdraft for totprec<0
       rnrtc(iq)=qprec-fldow(iq)*dprec         ! already has dsk factor
       if(ntest.eq.1.and.iq.eq.idjd)then
        print *,'kdown,tdown,qdown ',
     .           kdown,tdown,qdown
        print *,'qsk,qprec,dels0 ',qsk,qprec,hl*qprec 
        print *,'fldow,dprec,rnrtc ',
     .           fldow(iq),dprec,rnrtc(iq)
       endif
!      add in downdraft contributions
       dels(iq,kdown)=dels(iq,kdown)-fldow(iq)*s(iq,kdown)
       delq(iq,kdown)=delq(iq,kdown)-fldow(iq)*qsk
!      calculate emergent downdraft properties
       delq(iq,kbsav(iq))=delq(iq,kbsav(iq))+fldow(iq)*qdown
       dels(iq,kbsav(iq))=dels(iq,kbsav(iq))+fldow(iq)*
     .        (s(iq,kbsav(iq))+cp*(tdown-t(iq,kbsav(iq))))
!      add contrib to cloud base layer
       delq(iq,kbsav(iq))=delq(iq,kbsav(iq))-qbase(iq)
       dels(iq,kbsav(iq))=dels(iq,kbsav(iq))-sbase(iq)
      enddo  ! iq loop

c     do k=kl-1,kuocb+1,-1
      do k=kuocb+1,kl-1
       do iq=1,ifull
        kdown=min(kl , kbsav(iq)+nint(.6+.75*(ktsav(iq)-kbsav(iq))))
        if(k.ge.kdown+1.and.k.le.ktsav(iq))then
         savg=.5*(s(iq,k)+s(iq,k-1))
         qavg=.5*(qq(iq,k)+qq(iq,k-1))   
         dels(iq,k)=dels(iq,k)-savg      ! subsidence 
         dels(iq,k-1)=dels(iq,k-1)+savg  ! subsidence into lower layer
         delq(iq,k)=delq(iq,k)-qavg      ! subsidence 
         delq(iq,k-1)=delq(iq,k-1)+qavg  ! subsidence into lower layer
         if(ntest.eq.1.and.iq.eq.idjd)then
           print *,'in top loop; savg,qavg ',savg,qavg
         endif
        endif  ! (k.ge.kdown+1.and.k.le.ktsav(iq))
	 
        flm=1.-fldow(iq)
        if(k.le.kdown.and.k.gt.kbsav(iq))then
         savgb=.5*(s(iq,k)+s(iq,k-1))       
         qavgb=.5*(qq(iq,k)+qq(iq,k-1))   
         dels(iq,k)=dels(iq,k)-flm*savgb       ! subsidence 
         dels(iq,k-1)=dels(iq,k-1)+flm*savgb   ! subsidence into lower layer
         delq(iq,k)=delq(iq,k)-flm*qavgb       ! subsidence 
         delq(iq,k-1)=delq(iq,k-1)+flm*qavgb   ! subsidence into lower layer
        endif
       enddo  ! iq loop
      enddo   ! k loop
      
      if(ntest.ne.0.or.diag)then
        print *,"before diag print of dels,delh "
	 iq=idjd
        nlayersd=max(1,nint((ktsav(iq)-kbsav(iq)-.1)/methdetr)) ! round down
        khalfd=ktsav(iq)+1-nlayersd
        nlayersp=max(1,nint((ktsav(iq)-kbsav(iq)-.1)/methprec)) ! round down
        khalfp=ktsav(iq)+1-nlayersp
	 print *,'kbsav,ktsav,khalfd,nlayersd,khalfp,nlayersp ',
     .           kbsav(iq),ktsav(iq),khalfd,nlayersd,khalfp,nlayersp 
        print *,'dels1 ',(dels(idjd,k),k=1,kl)
        print *,'delh ',(dels(idjd,k)+hl*delq(idjd,k),k=1,kl)
        sum=0.
        do k=kbsav(idjd),ktsav(idjd)
         sum=sum+dels(idjd,k)+hl*delq(idjd,k)
        enddo
        print *,'qbase,sum_delh ',qbase(idjd),sum
      endif

!     calculate actual delq and dels
      do k=1,kl-1   
       do iq=1,ifull
        if(k.gt.kbsav(iq))then
          delq(iq,k)=delq(iq,k)/dsk(k)
          dels(iq,k)=dels(iq,k)/dsk(k)
        else
          if(nbase.eq.0.and.k.ge.kbsav_b(iq))then
!           assume cloud base layer, and those below it, are mixed well
c           delq(iq,k)=delq(iq,kbsav(iq))/(1.-sigmh(kbsav(iq)+1))
c           dels(iq,k)=dels(iq,kbsav(iq))/(1.-sigmh(kbsav(iq)+1))
            delq(iq,k)=delq(iq,kbsav(iq))/
     .                           (sigmh(kbsav_b(iq))-sigmh(kbsav(iq)+1))
            dels(iq,k)=dels(iq,kbsav(iq))/
     .                           (sigmh(kbsav_b(iq))-sigmh(kbsav(iq)+1))
          elseif(k.eq.kbsav(iq))then  ! s and qg only from kb layer
            delq(iq,k)=delq(iq,k)/dsk(k)
            dels(iq,k)=dels(iq,k)/dsk(k)
          endif  ! (nbase.eq.0)
        endif    ! (k.gt.kbsav(iq)) .. else ..
       enddo     ! iq loop
      enddo      ! k loop

      if(ntest.ne.0.or.diag)then
        print *,"before convpsav calc "
        print *,'delq ',(delq(idjd,k),k=1,kl)
        print *,'dels ',(dels(idjd,k),k=1,kl)
      endif

!     calculate base mass flux 
      do iq=1,ifull
       if(ktsav(iq).lt.kl)then 
!        Base limiter: new_qq(kb)=rhcv*new_qs(kb), i.e.
!        [qq+M*delq]_kb=rhcv*[qs+M*dqsdt*dels/cp]_kb 
         fluxb(iq)=max(0.,(rhcv*qs(iq,kbsav(iq))-qq(iq,kbsav(iq)))/ ! with dqsdt term
     .              ( delq(iq,kbsav(iq))
     .               -rhcv*dqsdt(iq,kbsav(iq))*dels(iq,kbsav(iq))/cp ) ) 
         convpsav(iq)=fluxb(iq)
	  if(nfluxq.eq.1)then  ! does little cf flux_dsk
!          fluxq limiter: new_qq(kb)=new_qq(kt)
!          i.e. [qq+M*delq]_kb=[qq+M*delq]_kt
           fluxq(iq)=max(0.,(qq(iq,kbsav(iq))-qq(iq,ktsav(iq)))/ 
     .                      (delq(iq,ktsav(iq)) - delq(iq,kbsav(iq))))
           convpsav(iq)=min(convpsav(iq),fluxq(iq))
	  endif  ! (nfluxq.eq.1)
	  if(nfluxdsk.eq.1)convpsav(iq)=min(convpsav(iq),dsk(kbsav(iq)))
	  if(nfluxdsk.eq.2)convpsav(iq)=
     .                              min(convpsav(iq),.5*dsk(kbsav(iq)))
       endif    ! (ktsav(iq).lt.kl)
      enddo     ! iq loop
      
      factr=dt/max(dt,convtime*3600.) ! was using 1.01 till 29/3/01
      do k=kuocb+1,kl-1
       do iq=1,ifull
        if(k.gt.kbsav(iq).and.k.le.ktsav(iq).and.ktsav(iq).lt.kl)then
!         want: new_hbas>=new_hs(k), i.e. in limiting case:
!         [h+alfsarr*M*dels+alfqarr*M*hl*delq]_base 
!                                          = [hs+M*dels+M*hlcp*dels*dqsdt]_k
          fluxt(iq)=max(0., (sbase(iq)+hl*qbase(iq)-hs(iq,k))/
     .     (dels(iq,k)*(1.+hlcp*dqsdt(iq,k))
     .        -alfsarr(iq)*dels(iq,kbsav(iq))
     .        -alfqarr(iq)*hl*delq(iq,kbsav(iq)) ) )
          convpsav(iq)=min(convpsav(iq),factr*fluxt(iq))
        endif   ! (k.gt.kbsav(iq).and.k.le.ktsav(iq).and.ktsav(iq).lt.kl)
       enddo    ! iq loop
       if(ntest.ne.0)then
	 iq=idjd
        if(k.gt.kbsav(iq).and.k.le.ktsav(iq).and.ktsav(iq).lt.kl)then
	   print *,'k,fluxt,convpsav,top,bott ',k,fluxt(iq),convpsav(iq),
     .     (sbase(iq)+hl*qbase(iq)-hs(iq,k)),
     .     ( dels(iq,k)*(1.+hlcp*dqsdt(iq,k))
     .         -dels(iq,kbsav(iq))-alfqarr(iq)*hl*delq(iq,kbsav(iq)) )
         if(itn.gt.1.and.ktsav(iq)-kbsav(iq).gt.6.
     .                    and.fluxt(iq).eq.0.)then
           print *,'ktau,itn,iq,kbsav,ktsav,ktmax ',
     .              ktau,itn,iq,kbsav(iq),ktsav(iq),ktmax(iq)
           print *,'hsdiff,fluxb ',
     .              hs(iq,kbsav(iq))-hs(iq,ktsav(iq)),fluxb(iq)
         endif
        endif   ! (k.gt.kbsav(iq).and.k.le.ktsav(iq).and.ktsav(iq).lt.kl)
       endif    ! (ntest.ne.0)
      enddo     ! k loop
      
      if(iterconv.gt.1)then 
        if(convfact.le.1.001)then 
          do iq=1,ifull          ! not needed if cloud base rising (19/7/02)
           ktmax(iq)=ktsav(iq)   ! ready for next itn iteration
          enddo  ! iq loop
        else
          do iq=1,ifull          
           convpsav(iq)=convfact*convpsav(iq) ! typically convfact=1.02  
          enddo  ! iq loop
        endif    ! (convfact.le.1.001)
      endif      ! (iterconv.gt.1)
      
!!!!!!!!!!!!!!!!!!!!! moistening & precip distribution now here!!!!!v2!!!!!!      

      if(methprec.eq.1)then     ! expel rain from top layer only
       do iq=1,ifull  
        rnrt3(iq,ktsav(iq))=rnrtc(iq)
       enddo  ! iq loop
      elseif(methprec.eq.9)then ! expel rain from all cloud layers
       do iq=1,ifull  
        nlayers=ktsav(iq)-kbsav(iq)
        do k=kbsav(iq)+1,ktsav(iq)
         rnrt3(iq,k)=rnrtc(iq)/nlayers
        enddo  
       enddo   ! iq loop
      else                      ! expel rain from upper layers - usual
       do k=kuocb+1,kl-1 
        do iq=1,ifull  
c        khalf=min( ktsav(iq) , max( (2-methprec)*ktsav(iq),
c    .         nint((kbsav(iq)+(methprec-1.)*ktsav(iq))/methprec+1.1) ))
c        nlayers=ktsav(iq)+1-khalf
         nlayers=max(1,nint((ktsav(iq)-kbsav(iq)-.1)/methprec)) ! round down
         khalf=ktsav(iq)+1-nlayers
         if(k.ge.khalf.and.k.le.ktsav(iq))then
          rnrt3(iq,k)=rnrtc(iq)/nlayers
         endif
        enddo  ! iq loop
       enddo   ! k loop
      endif    ! (methprec.eq.1)  .. else ..

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      if(ntest.ne.0.or.diag)then
        iq=idjd
        detr=max( detrain,detrain+(dsig4-sig(kbsav(iq))+sig(ktsav(iq)))*
     .                 (detrainx-detrain)/(dsig4-dsig2) )
        detr=min(detr,detrainx)
        print *,'detr,delta_sig,qxcess ',
     .           detr,sig(kbsav(iq))-sig(ktsav(iq)),qxcess(iq)
        flux_dsk=dsk(kbsav(idjd))*convfact
	 if(nfluxdsk.eq.2)flux_dsk=.5*dsk(kbsav(idjd))*convfact
        write(6,"('flux_dsk,fluxq,fluxb,fluxt,convpsav',5f6.3)")
     .       flux_dsk,fluxq(idjd),fluxb(idjd),fluxt(idjd),convpsav(idjd)
        write(6,"('delQ',18f6.3)")
     .	 (1000.*convpsav(idjd)*delq(idjd,k),k=1,kl)
        write(6,"('delt',18f6.3)")
     .	 (convpsav(idjd)*dels(idjd,k)/cp,k=1,kl)
        write(6,"('delQ*dsk',18f6.3)")
     .	 (1000.*convpsav(idjd)*delq(idjd,k)*dsk(k),k=1,kl)
        write(6,"('delt*dsk',18f6.3)")
     .	 (convpsav(idjd)*dels(idjd,k)*dsk(k)/cp,k=1,kl)
	 convmax=0.
	 nums=0
        do iq=1,ifull     
	  if(ktsav(iq)-kbsav(iq).eq.1)then
	    nums=nums+1
	    if(convpsav(iq).gt.convmax)then
	      convmax=convpsav(iq)
             write(6,"('bc  ',2i5,2i3,2x,3f6.3,f7.3,f6.3,' rh:',2f5.2)") 
     .       iq,nums,kbsav(iq),ktsav(iq),
     .       dsk(kbsav(iq)),fluxq(iq),fluxb(iq),fluxt(iq),convpsav(iq),
     .       qq(iq,kbsav(iq))/qs(iq,kbsav(iq)),
     .       qq(iq,ktsav(iq))/qs(iq,ktsav(iq))
           endif
         endif
        enddo  ! iq loop
	 convmax=0.
	 nums=0
        do iq=1,ifull     
	  if(ktsav(iq)-kbsav(iq).eq.2)then
	    nums=nums+1
	    if(convpsav(iq).gt.convmax)then
	      convmax=convpsav(iq)
             write(6,"('bcd ',2i5,2i3,2x,3f6.3,f7.3,f6.3,' rh:',2f5.2)") 
     .       iq,nums,kbsav(iq),ktsav(iq),
     .       dsk(kbsav(iq)),fluxq(iq),fluxb(iq),fluxt(iq),convpsav(iq),
     .       qq(iq,kbsav(iq))/qs(iq,kbsav(iq)),
     .       qq(iq,ktsav(iq))/qs(iq,ktsav(iq))
           endif
         endif
        enddo  ! iq loop
	 convmax=0.
	 nums=0
        do iq=1,ifull     
	  if(ktsav(iq)-kbsav(iq).eq.3)then
	    nums=nums+1
	    if(convpsav(iq).gt.convmax)then
	      convmax=convpsav(iq)
             write(6,"('bcde',2i5,2i3,2x,3f6.3,f7.3,f6.3,' rh:',2f5.2)") 
     .       iq,nums,kbsav(iq),ktsav(iq),
     .       dsk(kbsav(iq)),fluxq(iq),fluxb(iq),fluxt(iq),convpsav(iq),
     .       qq(iq,kbsav(iq))/qs(iq,kbsav(iq)),
     .       qq(iq,ktsav(iq))/qs(iq,ktsav(iq))
           endif
         endif
        enddo  ! iq loop
      endif    ! (ntest.ne.0.or.diag)

      
!     update qq, tt and precip
      do k=1,kl   
       do iq=1,ifull
        qq(iq,k)=qq(iq,k)+convpsav(iq)*delq(iq,k)
        tt(iq,k)=tt(iq,k)+convpsav(iq)*dels(iq,k)/cp
        rnrt3d(iq,k)=rnrt3d(iq,k)+convpsav(iq)*rnrt3(iq,k)*conrev(iq) !  gm/m**2/sec
       enddo    ! iq loop
      enddo     ! k loop

      if(ntest.ne.0.or.diag)then
        iq=idjd
!       N.B. convpsav(iq) is already mult by dt jlm: mass flux is convpsav/dt
        print *,"after convection: ktau,itn,kbsav,ktsav ",
     .                 ktau,itn,kbsav(idjd),ktsav(idjd)
        write (6,"('qgc ',19f7.3/(8x,19f7.3))") 
     .             (1000.*qq(idjd,k),k=1,kl)
        write (6,"('ttc ',19f7.2/(8x,19f7.2))") 
     .             (tt(idjd,k),k=1,kl)
        print *,'rnrtc,fldow,qxcess ',rnrtc(iq),fldow(iq),qxcess(iq)
	 print *,'ktsav,qbase,qs_ktsav,qq_ktsav ',
     .	          ktsav(iq),qbase(iq),qs(iq,ktsav(iq)),qq(iq,ktsav(iq)) 
        write (6,"('rnrt3 ',19f7.5/(8x,19f7.5))") 
     .             (rnrt3(idjd,k),k=1,kl)
        write (6,"('rnrt3d_c',19f7.3/(8x,19f7.3))") 
     .             (rnrt3d(idjd,k),k=1,kl)
	 iq=idjd
        delq_av=0.
        delt_av=0.
        rnd_av=0.
        heatlev=0.
        do k=1,kl
         delq_av=delq_av+dsk(k)*convpsav(iq)*delq(iq,k)
         delt_av=delt_av+dsk(k)*convpsav(iq)*dels(iq,k)/cp
         heatlev=heatlev+sig(k)*dsk(k)*convpsav(iq)*dels(iq,k)/cp
         rnd_av=rnd_av+convpsav(iq)*rnrt3(iq,k)*conrev(iq)
        enddo
        print *,'delq_av,delt_exp,rnd_av,rnd_exp ',
     .         delq_av,-delq_av*hl/cp,rnd_av,-delq_av*conrev(iq)
        if(delt_av.ne.0.)print *,'ktau,itn,kbsav,ktsav,delt_av,heatlev',
     .          ktau,itn,kbsav(idjd),ktsav(idjd),delt_av,heatlev/delt_av
      endif   ! (ntest.ne.0.or.diag)

      if(nuvconv.gt.0)then   !  momentum calculations, similar to above
       do k=1,kl   
        do iq=1,ifull
         delu(iq,k)=0.
         delv(iq,k)=0.
        enddo   ! iq loop
       enddo    ! k loop
       do iq=1,ifull
        k=ktsav(iq)   ! dels and delq for top layer
        delu(iq,k)=delu(iq,k)
     .             +(u(iq,kbsav(iq))-u(iq,k)) 
        delu(iq,k-1)=u(iq,k)    ! subsidence into lower layer
        delv(iq,k)=delv(iq,k)
     .             +(v(iq,kbsav(iq))-v(iq,k))
         delv(iq,k-1)=v(iq,k)    ! subsidence into lower layer
        do k=ktsav(iq)-1,kbsav(iq)+1,-1
         delu(iq,k)=delu(iq,k)-u(iq,k) ! subs.
         delu(iq,k-1)=u(iq,k)          ! subsidence into lower layer
         delv(iq,k)=delv(iq,k)-v(iq,k) ! subs.
         delv(iq,k-1)=v(iq,k)          ! subsidence into lower layer
        enddo
!       assume cloud base layer, and those below it, are mixed well
        delu(iq,kbsav(iq))=delu(iq,kbsav(iq))-u(iq,kbsav(iq))
        delv(iq,kbsav(iq))=delv(iq,kbsav(iq))-v(iq,kbsav(iq))
       enddo  ! iq loop
!      calculate actual delu and delv
       do k=1,kl-1   
        do iq=1,ifull
         if(k.gt.kbsav(iq))then
           delu(iq,k)=delu(iq,k)/dsk(k)
           delv(iq,k)=delv(iq,k)/dsk(k)
         else
           if(nbase.eq.0.and.k.ge.kbsav_b(iq))then
!            assume cloud base layer, and those below it, are mixed well
             delu(iq,k)=delu(iq,kbsav(iq))/
     .                           (sigmh(kbsav_b(iq))-sigmh(kbsav(iq)+1))
             delv(iq,k)=delv(iq,kbsav(iq))/
     .                           (sigmh(kbsav_b(iq))-sigmh(kbsav(iq)+1))
           elseif(k.eq.kbsav(iq))then  ! u and v only from kb layer
             delu(iq,k)=delu(iq,k)/dsk(k)
             delv(iq,k)=delv(iq,k)/dsk(k)
           endif  ! (nbase.eq.0)
         endif    ! (k.gt.kbsav(iq)) .. else ..
        enddo     ! iq loop
       enddo      ! k loop
!      update u, v
       do k=1,kl   
        do iq=1,ifull
         u(iq,k)=u(iq,k)+convpsav(iq)*delu(iq,k)
         v(iq,k)=v(iq,k)+convpsav(iq)*delv(iq,k)
        enddo    ! iq loop
       enddo     ! k loop
      endif     ! (nuvconv.gt.0)
      
      if(ntest.ne.0.or.diag)then
        write (6,"('uuc ',19f6.1/(8x,19f6.1))") 
     .             (u(idjd,k),k=1,kl)
        write (6,"('vvc ',19f6.1/(8x,19f6.1))") 
     .             (v(idjd,k),k=1,kl)
      endif
 
      enddo     ! itn=1,iterconv

!__________________________end of convective calculations_____________________
!     now evap of convective rainfall, and vertical accumulation of rnrt3d
      if(nevapcc.eq.2)then  ! jlm
        do k=kl-2,1,-1
         do iq=1,ifull
          revc(iq,k)=
     .      rhmois*3.*(sig(k)-sig(k+1))*(1.-qq(iq,k)/qs(iq,k)) ! e.g. rhmois=3.
     .      *rnrt3d(iq,k+1)/(conrev(iq)*.1)   ! jlm suggestion
          revc(iq,k)=min(revc(iq,k),rnrt3d(iq,k+1)/(conrev(iq)*dsk(k)))
!         max needed for roundoff
          rnrt3d(iq,k)=rnrt3d(iq,k)
     .        +max(0., rnrt3d(iq,k+1)-revc(iq,k)*dsk(k)*conrev(iq))
          tt(iq,k)=tt(iq,k)-revc(iq,k)*hlcp
          qq(iq,k)=qq(iq,k)+revc(iq,k)
         enddo    ! iq loop
        enddo     ! k loop
      endif       ! (nevapcc.eq.2)

      if(ntest.ne.0.or.diag)then
        print *,"before nevapcc=5 section"
        write (6,"('rnrt3d_d',19f7.3/(8x,19f7.3))") 
     .             (rnrt3d(idjd,k),k=1,kl)
      endif   ! (ntest.ne.0.or.diag)

      if(nevapcc.eq.5)then  !  like nevapls=5
        rKa=2.4e-2
        rvap=461.
        Dva=2.21
        cfls=1.     ! cld frac
        cfls=rhmois ! cld frac
        cflscon=4560.*cfls**.3125
        do k=kl-2,1,-1
         do iq=1,ifull
c         if(k.eq.ktsav(iq))then
c           cfls=1.
c         else
c           cfls=rhmois
c         endif
!         fluxr(iq)=rnrt(iq)*1.e-3*dt ! kg/m2
          fluxc=rnrt3d(iq,k+1)*1.e-3*dt ! kg/m2
          rhodz=ps(iq)*dsk(k)/grav
          qpf=fluxc/rhodz     ! Mix ratio of rain which falls into layer
          pk=ps(iq)*sig(k)
          es=qs(iq,k)*pk/.622
          Apr=hl*hl/(rKa*tt(iq,k)*(rvap*tt(iq,k))-1.)
          Bpr=rvap*tt(iq,k)*pk/(Dva*es)
          Fr=fluxc/(cfls*dt)
c	   if(Fr.lt.0.)then
c	     print *,'negative flux for ktau,iq,k,cfls ',ktau,iq,k,cfls
c           print *,'rnrt3d ',(rnrt3d(iq,kk),kk=1,kl)
c	     print *,'kbsav,ktsav ',kbsav(iq),ktsav(iq)
c	     stop
c	   endif
          rhoa=pk/(r*tt(iq,k))
          dz=pk/(rhoa*grav)
          Vr=max(.01 , 11.3*Fr**(1./9.)/sqrt(rhoa)) ! Actual fall speed
          dtev=dz/Vr
          qr=fluxc/(dt*rhoa*Vr)
          qgdiff=qs(iq,k)-qq(iq,k)
          Cev2=cflscon*qgdiff/(qs(iq,k)*(Apr+Bpr)) ! Ignore rhoa**0.12
          qr2=max(0. , qr**.3125 - .3125*Cev2*dtev)**3.2
          Cevx=(qr-qr2)/dtev  ! i.e. Cev*qgdiff
          alphal=hl*qs(iq,k)/(ars*tt(iq,k)**2)
          blx=qgdiff+Cevx*dt*(1.+hlcp*alphal)  ! i.e. bl*qgdiff
          evapls=cfls*dt*Cevx*qgdiff/blx ! UKMO
          evapls=max(0. , min(evapls,qpf))
          revc(iq,k)=min(evapls , rnrt3d(iq,k+1)/(conrev(iq)*dsk(k)))
!         max needed for roundoff
          rnrt3d(iq,k)=rnrt3d(iq,k)
     .        +max(0., rnrt3d(iq,k+1)-revc(iq,k)*dsk(k)*conrev(iq))
          tt(iq,k)=tt(iq,k)-revc(iq,k)*hlcp
          qq(iq,k)=qq(iq,k)+revc(iq,k)
         enddo    ! iq loop
        enddo     ! k loop
      endif       ! (nevapcc.eq.5)

      if(nevapcc.eq.0)then
        do k=2,kl-1
         do iq=1,ifull
          rnrt3d(iq,1)=rnrt3d(iq,1)+rnrt3d(iq,k)
         enddo    ! iq loop
        enddo     ! k loop
      endif       ! (nevapcc.eq.0) .. else ..
      do iq=1,ifull
       rnrtc(iq)=rnrt3d(iq,1) ! save in convective rainfall array
      enddo  ! iq loop

      if(ntest.ne.0.or.diag)then
        print *,"after  convection & evap"
        write (6,"('qge ',19f7.3/(8x,19f7.3))") 
     .             (1000.*qq(idjd,k),k=1,kl)
        write (6,"('tte ',19f6.1/(8x,19f6.1))") 
     .             (tt(idjd,k),k=1,kl)
        write (6,"('rnrt3d  ',19f7.3/(8x,19f7.3))") 
     .             (rnrt3d(idjd,k),k=1,kl)
        print *,'rnrtc ',rnrtc(idjd)
!       following calc of total heatlev really needs iterconv=1	 
	 do k=kl-2,1,-1
	  delt_av=delt_av-dsk(k)*revc(idjd,k)*hlcp
	  heatlev=heatlev-sig(k)*dsk(k)*revc(idjd,k)*hlcp
	  print *,'k,rh,delt_av, heatlev ',
     .            k,100.*qq(idjd,k)/qs(idjd,k),delt_av,heatlev
	 enddo
	 if(delt_av.ne.0.)print *,'ktau,delt_av-net,heatlev_net ',
     .                            ktau,delt_av,heatlev/delt_av
      endif
      if(nkuo.lt.0)go to 8

!__________________________beginning of large-scale calculations_____________________
c     check for grid-scale rainfall 
5     do k=kl,1,-1    ! top down to end up with proper kbsav_ls
       do iq=1,ifull
        es=establ(tt(iq,k))
        pk=ps(iq)*sig(k)
        qs(iq,k)=max(.622*es/(pk-es),1.5e-6)  ! Sat  10-31-1998
        if(qq(iq,k).gt.rhsat*qs(iq,k))then
          kbsav_ls(iq)=k
          gam=max(hlcp*qs(iq,k)*pk*hlars/(tt(iq,k)**2*(pk-es)),0.) 
          dqrx=(qq(iq,k)-rhsat*qs(iq,k))/(1.+rhsat*gam)
          tt(iq,k)=tt(iq,k)+hlcp*dqrx
          qq(iq,k)=qq(iq,k)-dqrx
          rnrt(iq)=rnrt(iq)+dqrx*dsk(k)
        endif   ! (qq(iq,k).gt.rhsat*qs(iq,k))
       enddo    ! iq loop
      enddo     ! k loop

c     now do evaporation of L/S precip
c      Here rnrt is the rainfall rate in gm/m**2/sec
!      conrev(iq)=1000.*ps(iq)/(grav*dt)     
       rnrt=rnrt*conrev                 !2D
       fluxr=rnrt*1.e-3*dt ! kg/m2      !2D

      if(ntest.ne.0.or.diag)then
        print *,'after large scale rain: kbsav_ls,rnrt ',
     .                                   kbsav_ls(idjd),rnrt(idjd)
        write (6,"('qg ',19f7.3/(8x,19f7.3))") 
     .             (1000.*qq(idjd,k),k=1,kl)
        write (6,"('tt ',19f7.2/(8x,19f7.2))") 
     .             (tt(idjd,k),k=1,kl)
      endif

      if(nevapls.eq.1)then
        do k=2*kl/3,1,-1
         do iq=1,ifull
          if(k.lt.kbsav_ls(iq))then
            revq=dt*ecmwls*ecmwa1*max((qs(iq,k)-qq(iq,k)),0.)
     .                     *sqrt(sqrt(sig(k))*rnrt(iq)/ecmwa2)
            revq=min(revq,rnrt(iq)/(conrev(iq)*dsk(k)))
!           max needed for roundoff
            rnrt(iq)=max(1.e-10,rnrt(iq)-revq*dsk(k)*conrev(iq))
            tt(iq,k)=tt(iq,k)-revq*hlcp
            qq(iq,k)=qq(iq,k)+revq
          endif   !  (k.lt.kbsav_ls(iq))
         enddo    ! iq loop
        enddo     ! k loop
      endif       ! (nevapls.eq.1)

      if(nevapls.eq.2)then
        do k=2*kl/3,1,-1
         do iq=1,ifull
          if(k.lt.kbsav_ls(iq))then
            revq=
     .        3.*(sig(k)-sig(kbsav_ls(iq)))*(1.-qq(iq,k)/qs(iq,k))
     .        *rnrt(iq)/(conrev(iq)*.1)   ! jlm suggestion
            revq=min(revq,rnrt(iq)/(conrev(iq)*dsk(k)))
!           max needed for roundoff
            rnrt(iq)=max(1.e-10,rnrt(iq)-revq*dsk(k)*conrev(iq))
            tt(iq,k)=tt(iq,k)-revq*hlcp
            qq(iq,k)=qq(iq,k)+revq
          endif   !  (k.lt.kbsav_ls(iq))
         enddo    ! iq loop
        enddo     ! k loop
      endif       ! (nevapls.eq.2)

      if(nevapls.eq.3)then ! from ldr newrain.f 29 june 95
        do k=2*kl/3,1,-1
         do iq=1,ifull
          if(k.lt.kbsav_ls(iq))then
            rhodz=ps(iq)*dsk(k)/grav
            qpf=fluxr(iq)/rhodz     !Mix ratio of rain which falls into layer
            Cev=859.*Aev(tt(iq,k))*sqrt(fluxr(iq)/dt)
            alphal=hl*qs(iq,k)/(ars*tt(iq,k)**2)
            bl=1+Cev*dt*(1+hlcp*alphal)
            cfls=1. ! cld frac large scale
            evapls= cfls*dt*(Cev/bl)*(qs(iq,k)-qq(iq,k)) !UKMO
            evapls=max(0.,min(evapls,qpf))
            qpf=qpf-evapls
            fluxr(iq)=fluxr(iq)+rhodz*qpf
            revq=evapls
            revq=min(revq,rnrt(iq)/(conrev(iq)*dsk(k)))
!           max needed for roundoff
            rnrt(iq)=max(1.e-10,rnrt(iq)-revq*dsk(k)*conrev(iq))
            tt(iq,k)=tt(iq,k)-revq*hlcp
            qq(iq,k)=qq(iq,k)+revq
          endif   !  (k.lt.kbsav_ls(iq))
         enddo    ! iq loop
        enddo     ! k loop
      endif       ! (nevapls.eq.3)

      if(nevapls.eq.4)then ! newer UKMO Thu  03-05-1998
        do k=2*kl/3,1,-1
         do iq=1,ifull
          if(k.lt.kbsav_ls(iq))then
            rhodz=ps(iq)*dsk(k)/grav
            qpf=fluxr(iq)/rhodz     !Mix ratio of rain which falls into layer
            if( tt(iq,k).gt.273.15 )then
              Cev=405.*Aev(tt(iq,k))*sqrt(fluxr(iq)/dt)
            else
              Cev=3200.*Asb(tt(iq,k))*sqrt(fluxr(iq)/dt)
            endif
            alphal=hl*qs(iq,k)/(ars*tt(iq,k)**2)
            bl=1+Cev*dt*(1+hlcp*alphal)
            cfls=1. ! cld frac large scale
            evapls= cfls*dt*(Cev/bl)*(qs(iq,k)-qq(iq,k)) !UKMO
            evapls=max(0.,min(evapls,qpf))
            qpf=qpf-evapls
            fluxr(iq)=fluxr(iq)+rhodz*qpf
            revq=evapls
            revq=min(revq,rnrt(iq)/(conrev(iq)*dsk(k)))
!           max needed for roundoff
            rnrt(iq)=max(1.e-10,rnrt(iq)-revq*dsk(k)*conrev(iq))
            tt(iq,k)=tt(iq,k)-revq*hlcp
            qq(iq,k)=qq(iq,k)+revq
          endif   !  (k.lt.kbsav_ls(iq))
         enddo    ! iq loop
        enddo     ! k loop
      endif       ! (nevapls.eq.4)

      if(nevapls.eq.5)then ! even newer UKMO Thu  03-05-1998
        rKa=2.4e-2
        rvap=461.
        Dva=2.21
        cfls=1. ! cld frac7 large scale
        cflscon=4560.*cfls**.3125
        do k=2*kl/3,1,-1
         do iq=1,ifull
          if(k.lt.kbsav_ls(iq))then
            rhodz=ps(iq)*dsk(k)/grav
            qpf=fluxr(iq)/rhodz     ! Mix ratio of rain which falls into layer
            pk=ps(iq)*sig(k)
            es=qs(iq,k)*pk/.622
            Apr=hl*hl/(rKa*tt(iq,k)*(rvap*tt(iq,k))-1.)
            Bpr=rvap*tt(iq,k)*pk/(Dva*es)
            Fr=fluxr(iq)/(cfls*dt)
            rhoa=pk/(r*tt(iq,k))
            dz=pk/(rhoa*grav)
            Vr=max(.01 , 11.3*Fr**(1./9.)/sqrt(rhoa)) ! Actual fall speed
            dtev=dz/Vr
            qr=fluxr(iq)/(dt*rhoa*Vr)
            qgdiff=qs(iq,k)-qq(iq,k)
            Cev2=cflscon*qgdiff/(qs(iq,k)*(Apr+Bpr)) ! Ignore rhoa**0.12
            qr2=max(0. , qr**.3125 - .3125*Cev2*dtev)**3.2
!           Cev=(qr-qr2)/(dtev*qgdiff)
            Cevx=(qr-qr2)/dtev  ! i.e. Cev*qgdiff
            alphal=hl*qs(iq,k)/(ars*tt(iq,k)**2)
!           bl=1.+Cev*dt*(1.+hlcp*alphal)
            blx=qgdiff+Cevx*dt*(1.+hlcp*alphal)  ! i.e. bl*qgdiff
!           evapls= cfls*dt*qgdiff*Cev/bl ! UKMO
            evapls= cfls*dt*Cevx*qgdiff/blx ! UKMO
            evapls=max(0. , min(evapls,qpf))
            qpf=qpf-evapls
c           if(evapls.gt.0.)print *,'iq,k,evapls ',iq,k,evapls
            fluxr(iq)=fluxr(iq)+rhodz*qpf
            revq=evapls
            revq=min(revq , rnrt(iq)/(conrev(iq)*dsk(k)))
!           max needed for roundoff
            rnrt(iq)=max(1.e-10 , rnrt(iq)-revq*dsk(k)*conrev(iq))
            tt(iq,k)=tt(iq,k)-revq*hlcp
            qq(iq,k)=qq(iq,k)+revq
          endif   !  (k.lt.kbsav_ls(iq))
         enddo    ! iq loop
        enddo     ! k loop
      endif       ! (nevapls.eq.5)
      if(nkuo.lt.0)go to 2

!__________________________end of large-scale calculations_____________________


8     if(ntest.ne.0.or.diag)then
        print *,'near end of convjlm: rnrt,rnrtc ',
     .                                  rnrt(idjd),rnrtc(idjd)
        write (6,"('qg ',19f7.3/(8x,19f7.3))") 
     .             (1000.*qq(idjd,k),k=1,kl)
        write (6,"('tt ',19f7.2/(8x,19f7.2))") 
     .             (tt(idjd,k),k=1,kl)
      endif

!     section for convective transport of trace gases (jlm 22/2/01)
      if(ilt.gt.1)then
        if(iterconv.ne.1)stop 'need 1 for trace gases'
        do ntr=1,ntrac
         do k=1,kl
	   do iq=1,ifull
           s(iq,k)=tr(iq,k,ntr)
          enddo    ! iq loop
         enddo     ! k loop
         do iq=1,ifull
          if(ktsav(iq).lt.kl)then
            kb=kbsav(iq)
            kt=ktsav(iq)
            veldt=convpsav(iq)*(1.-fldow(iq))  ! simple treatment
            fluxup=veldt*s(iq,kb)
!           remove gas from cloud base layer
            tr(iq,kb,ntr)=tr(iq,kb,ntr)-fluxup/dsk(kb)
!           put flux of gas into top convective layer
            tr(iq,kt,ntr)=tr(iq,kt,ntr)+fluxup/dsk(kt)
            do k=kb+1,kt
             tr(iq,k,ntr)=tr(iq,k,ntr)-s(iq,k)*veldt/dsk(k)
             tr(iq,k-1,ntr)=tr(iq,k-1,ntr)+s(iq,k)*veldt/dsk(k-1)
            enddo
          endif
         enddo   ! iq loop
        enddo    ! ntr loop    
      endif      ! ilt.gt.1

      if(npanels.eq.0)then     ! for DARLAM
        i1=3       ! for DARLAM
        j1=3       ! for DARLAM
        imax=il-1  ! for DARLAM
        jmax=jl-1  ! for DARLAM
        do j=j1,jmax
         do i=i1,imax
          iq=i+(j-1)*il
          prcon= .001*dt*rnrtc(iq)
          condc(iq)=prcon              ! convective precip for this timestep
          precc(iq)=precc(iq)+prcon
          prl_s= .001*dt*rnrt(iq)
          condx(iq)=prcon+prl_s        ! total precip for this timestep
          precip(iq)=precip(iq)+prcon+prl_s
         enddo  !  i loop
        enddo   !  j loop
        if(ndavconv.eq.0)then  
          factdav=1.  !2D
        else          ! i.e. for ndavconv=1 or 2
          factdav=1.-davt  !2D to reduce heating near boundaries
        endif         ! (ndavconv.eq.0) .. else ..
        if(ndavconv.eq.2)then     ! davies-style for qg near boundaries too
          do k=1,kl
           do j=j1,jmax
            do i=i1,imax
             iq=i+(j-1)*il
             qg(iq,k)=qg(iq,k)+factdav(iq)*(qq(iq,k)-qg(iq,k))
            enddo  !  i loop
           enddo   !  j loop
          enddo    !  k loop
        else         ! i.e. for ndavconv=1 don't reduce moisture changes
          do k=1,kl
           do j=j1,jmax
            do i=i1,imax
             iq=i+(j-1)*il
             qg(iq,k)=qq(iq,k)
            enddo  !  i loop
           enddo   !  j loop
          enddo    !  k loop
        endif      ! (ndavconv.eq.2)
        if(abs(nkuo).eq.23)then  ! split T at end of timestep (tn not used)
          do k=1,kl
           do j=j1,jmax
            do i=i1,imax
             iq=i+(j-1)*il
             t(iq,k)=tt(iq,k)
            enddo  !  i loop
           enddo   !  j loop
          enddo    !  k loop
        else
          do k=1,kl
           do j=j1,jmax
            do i=i1,imax
             iq=i+(j-1)*il
             tn(iq,k)=tn(iq,k)+factdav(iq)*(tt(iq,k)-t(iq,k))/dt
            enddo  !  i loop
           enddo   !  j loop
          enddo    !  k loop
        endif      ! (abs(nkuo).eq.23)
      else         ! usual conformal-cubic
        qg=qq      ! 3D
        do iq=1,ifull
         prcon= .001*dt*rnrtc(iq)
         condc(iq)=prcon              ! convective precip for this timestep
         precc(iq)=precc(iq)+prcon
         prl_s= .001*dt*rnrt(iq)
         condx(iq)=prcon+prl_s        ! total precip for this timestep
         precip(iq)=precip(iq)+prcon+prl_s
        enddo
        if(abs(nkuo).eq.23)then  ! split T at end of timestep (tn not used)
          t=tt
        else
          tn=tn+(tt-t)/dt
        endif
      endif        ! (npanels.eq.0)

      if(ntest.eq.1.or.diag)then
        pwater=0.   ! in mm     
	 do k=1,kl
         pwater=pwater-dsig(k)*qg(idjd,k)*ps(idjd)/g
	 enddo
	 print *,'pwater0,pwater+condx,pwater ',
     .	          pwater0,pwater+condx(idjd),pwater
        print *,'D rnrt,rnrtc,condx ',
     .             rnrt(idjd),rnrtc(idjd),condx(idjd)
        print *,'precc,precip ',
     .           precc(idjd),precip(idjd)
        call maxmin(rnrtc,'rc',ktau,1.,1)
      endif
      return
      end
