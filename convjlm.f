      subroutine convjlm      ! jlm convective scheme - Version v3
!     has fldownn depending on delta sigma      
      use cc_mpi, only : mydiag
      use diag_m
      include 'newmpar.h'
      parameter (ntest=0)      ! 1 or 2 to turn on
!     parameter (iterconv=1)   ! to kuocom.h
!     parameter (fldown=.6)    ! to kuocom.h
!     parameter (detrain=.05)  ! to kuocom.h
      parameter (alfcon=.015)  ! e.g. .01 or .015   or 1.1 (over sea)
!     parameter (alflnd=1.15)  ! e.g. 1.15  with nqbot=3
!     parameter (alfsea=1.05)  ! e.g. 1.05  with nqbot=3
!     parameter (nqbot=3)      ! 0: original qq(1), 1: qs1, 2: qbas for land/sea
!                                3:for RH-enhanced base
      parameter (nalfs=0)      ! 1: was usual, 0: off
      parameter (nfluxq=1)     ! 1 from 6/8/02
      parameter (nfluxdsk=-2)  ! 1 till 6/8/02
      parameter (nbase=1)      ! 0 shared below, 1 uses only one (originally 0)
!     parameter (methdetr=8)   ! 1 (top only); 2 (top half); 4 top quarter (= kbconv)
!     parameter (methprec=2)   ! 1 (top only); 2 (top half); 4 top quarter (= kbconv)
      parameter (no2layer=0)   ! usually 0, 1 to suppress 2-layer clouds (even wetter!)
      parameter (nuvconv=0)    ! usually 0, 1 to turn on momentum mixing
      parameter (kcl_top=kl-2) !max level for cloud top (convjlm,radrive,vertmix)
!     parameter (dsig2=.1, dsig4=.7)  ! used with detrainx
c     nevapls:  turn off/on ls evap - through parm.h; 0 off, 5 newer UK
      include 'arrays.h'
      include 'const_phys.h'
      include 'dava.h'     ! davt
      include 'kuocom.h'   ! kbsav,ktsav,convfact,convpsav,ndavconv
      include 'liqwpar.h'  ! ifullw
      include 'map.h'      ! land
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
     .  ,fluxb(ifull),kbsav_b(ifull),kbsav_ls(ifull)
     .  ,rnrt(ifull),rnrtc(ifull),sumqs(ifull),sumqt(ifull)
     .  ,sumss(ifull),factdav(ifull),qbase(ifull),fldow(ifull)
     .  ,qxcess(ifull),fluxq(ifull),dum2(ifull,2)
      common/work2b/sbase(ifull),alfqarr(ifull),alfsarr(ifull)
      common/work2c/rnrtcn(ifull),kb_sav(ifull),kt_sav(ifull)
      common/work3/tt(ifull,kl),qq(ifull,kl),qs(ifull,kl)
     .  ,dqsdt(ifull,kl),es(ifull,kl)  ! es array mainly for diag prints
      common/work3b/s(ifull,kl),hs(ifull,kl)
      common/work3d/fluxt(ifull,kl)   ! 3d for diags
      common/work3f/delq(ifull,kl),dels(ifull,kl),qliqw(ifull,kl)
      real dsk(kl),phi(ifull,kl),delu(ifull,kl),delv(ifull,kl)
      real revc(ifull,kl)  ! just for evap diags
      real alfss(kl),alfm(kl),qbass(kl)
      real alfqq_l(kl),alfqq_s(kl),alfqq1_l(kl),alfqq1_s(kl)
      equivalence (phi,s),(delu,revc,delq),(delv,dels)
!     data rhcv/.75/                  now in kuocom

      include 'establ.h'
      Aev(tm) = 2.008e-9*tm**2 - 1.385e-6*tm + 2.424e-4  !For UKMO evap scheme
      Asb(tm) = max (-5.2e-9*tm**2+2.5332e-6*tm-2.9111e-4,1.e-5) !For UKMO subl

      do k=1,kl
       dsk(k)=-dsig(k)    !   dsk = delta sigma (positive)
       alfm(k)=min(.5,alfcon/(max(1.e-10,sig(1)-sig(k))))  ! e.g. alfcon=.01
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
      qliqw(:,:)=0.  
      conrev(1:ifull)=1000.*ps(1:ifull)/(grav*dt) ! factor to convert precip to g/m2/s
      rnrt(:)=0.       ! initialize large-scale rainfall array
      rnrtc(:)=0.      ! initialize convective  rainfall array
      ktmax(:)=kl      ! preset 1 level above current topmost-permitted ktsav
      kbsav_ls(:)=0    ! for L/S

      tt(1:ifull,:)=t(1:ifull,:)       
      qq(1:ifull,:)=qg(1:ifull,:)      
      factr=dt/max(dt,convtime*3600.)  ! to re-scale convpsav

!__________________________beginning of convective calculations_____________________
      do itn=1,iterconv
!     calculate geopotential height 
      kb_sav(:)=kl     ! preset value for no convection
      kbsav_b(:)=kl    ! preset value for no convection
      kt_sav(:)=kl     ! preset value for no convection
      rnrtcn(:)=0.     ! initialize convective rainfall array (pre convpsav)
      convpsav(:)=0.
      dels(:,:)=1.e-20
      delq(:,:)=0.
      phi(:,1)=bet(1)*tt(:,1)

      do k=2,kl
       do iq=1,ifull
        phi(iq,k)=phi(iq,k-1)+bet(k)*tt(iq,k)+betm(k)*tt(iq,k-1)
       enddo     ! iq loop
      enddo      ! k  loop

      do k=1,kl   
       do iq=1,ifull
        es(iq,k)=establ(tt(iq,k))
       enddo  ! iq loop
      enddo   ! k loop
      do k=1,kl   
       do iq=1,ifull
        pk=ps(iq)*sig(k)
        qs(iq,k)=max(.622*es(iq,k)/(pk-es(iq,k)),1.5e-6)  
        dqsdt(iq,k)=qs(iq,k)*pk*hlars/(tt(iq,k)**2*(pk-es(iq,k)))
        s(iq,k)=cp*tt(iq,k)+phi(iq,k)  ! dry static energy
c       calculate hs
        hs(iq,k)=s(iq,k)+hl*qs(iq,k)   ! saturated moist static energy
       enddo  ! iq loop
      enddo   ! k loop

      if(ktau.eq.1.and.mydiag)then
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
      if((ntest.ne.0.or.diag).and.mydiag)then
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
         pwater0=pwater0-dsig(k)*qg(idjd,k)*ps(idjd)/grav
         if(land(iq))then
           qbass(k)=alfqq_l(k)*qq(iq,k)+alfqq1_l(k)*qq(iq,1) 
         else
           qbass(k)=alfqq_s(k)*qq(iq,k)+alfqq1_s(k)*qq(iq,1) 
         endif
	  qbass(k)=min(qbass(k),max(qq(iq,k),qs(iq,k)))  ! added 3/6/03
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
        print *,'es ',(es(idjd,k),k=1,kl)
        print *,'dqsdt ',(dqsdt(idjd,k),k=1,kl)
      endif

      sbase(:)=0.
      qbase(:)=0.
      do k=kl/2,kuocb+1,-1   ! downwards to find cloud base and kbsav_b
       do iq=1,ifull
!       find tentative cloud base, and bottom of below-cloud layer
        if(qq(iq,k-1).gt.rhcv*qs(iq,k-1))then
!         next line ensures sub-cloud layer moist and contiguous 
          if(kbsav_b(iq).eq.k.and.sig(kb_sav(iq)).gt..8)kbsav_b(iq)=k-1 
!         if(kbsav_b(iq).eq.k)kbsav_b(iq)=k-1 
          if(land(iq))then
            qbas=alfqq_l(k-1)*qq(iq,k-1)+alfqq1_l(k-1)*qq(iq,1) 
            alfq=alfqq_l(k-1)   
          else
            qbas=alfqq_s(k-1)*qq(iq,k-1)+alfqq1_s(k-1)*qq(iq,1) 
            alfq=alfqq_s(k-1)   
          endif  ! (land(iq)) .. else ..
	   qbas=min(qbas,max(qs(iq,k-1),qq(iq,k-1)))    ! added 3/6/03
          alfs=alfss(k-1)  ! simple one
          sbas=alfs*s(iq,k-1)+(1.-alfs)*s(iq,1) ! e.g. alfs=1.
!         now choose base layer to also have local max of qbase      
          if(sbas+hl*qbas.gt.max(hs(iq,k),sbase(iq)+hl*qbase(iq)))then
            if(qbas.gt.max(qs(iq,k),qq(iq,k)))then  ! newer qbas test 
              kb_sav(iq)=k-1
                qbase(iq)=qbas
                sbase(iq)=sbas
              kbsav_b(iq)=k-1
              kt_sav(iq)=k
              alfqarr(iq)=alfq    ! for use in flux calc
              alfsarr(iq)=alfs    ! for use in flux calc
            endif ! (qbas.gt.max(qs(iq,k),qq(iq,k)))
          endif   ! (sbas+hl*qbas.gt.hs(iq,k))
c           if(ntest.eq.1.and.iq.eq.idjd)then   ! needs sopt
c             print *,'k,(sbas+hl*qbas)/cp,hs(iq,k)/cp ',
c     .          k,(sbas+hl*qbas)/cp,hs(iq,k)/cp
c             print *,'kbsav,qbas,qs(.,k),qq(.,k) ',
c     .          kb_sav(iq),qbas,qs(iq,k),qq(iq,k)
c           endif ! (ntest.eq.1.and.iq.eq.idjd)
        endif   ! (qq(iq,k-1).gt.rhcv*qs(iq,k-1))
       enddo    ! iq loop
      enddo     ! k loop

      do k=kuocb+2,kl   ! upwards to find cloud top
       do iq=1,ifull
        if(kt_sav(iq).eq.k-1.and.k.lt.ktmax(iq))then ! ktmax allows for itn
          hbase=sbase(iq)+hl*qbase(iq)
          if(hbase.gt.hs(iq,k))kt_sav(iq)=k
        endif
       enddo    ! iq loop
      enddo     ! k loop

      if(no2layer.eq.1)then
!       present scheme does not properly handle 1- or 2-layer clouds
        do iq=1,ifull
         if(kt_sav(iq)-kb_sav(iq).le.no2layer)then
            kb_sav(iq)=kl
         endif
        enddo  ! iq loop
      endif

      if(no2layer.lt.0)then
!       present scheme does not properly handle 1- or 2-layer clouds,
!       so leave them to shallow convection scheme
        do iq=1,ifull
         if(sig(kt_sav(iq)).gt.sigksct)then
c        if(sig(kt_sav(iq)).gt..8)then
            kb_sav(iq)=kl
         endif
        enddo  ! iq loop
      endif

      if((ntest.ne.0.or.diag).and.mydiag)then
       iq=idjd
       iq=idjd
       kdown=min(kl,kb_sav(iq)+nint(.6+.75*(kt_sav(iq)-kb_sav(iq))))
       print *,"ktau,itn,kbsav_b,kbsav,ktsav,kdown ",
     .    ktau,itn,kbsav_b(idjd),kb_sav(idjd),kt_sav(idjd),kdown
       print *,'alfqarr,alfsarr ',alfqarr(idjd),alfsarr(idjd)
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     calculate dels and delq for all layers for unit base mass flux
!     designed to be OK if q(k)>qs(k)
      do iq=1,ifull
!      dels and delq for top layer
!      following line allows for qq>qs, to avoid drying layer
       qsk=max(qs(iq,kt_sav(iq)),qq(iq,kt_sav(iq)))  
       qprec=max(0.,qbase(iq)-qsk)             
       dels(iq,kt_sav(iq))=dels(iq,kt_sav(iq))+hl*qprec        ! precip. heating
       dels(iq,kt_sav(iq))=dels(iq,kt_sav(iq))+sbase(iq)   ! s flux
       delq(iq,kt_sav(iq))=delq(iq,kt_sav(iq))+qsk  

       kdown=min(kl , kb_sav(iq)+nint(.6+.75*(kt_sav(iq)-kb_sav(iq))))
       qsk=min(qs(iq,kdown),qq(iq,kdown))    ! N.B. min for downdraft subs
       tdown=t(iq,kb_sav(iq))+(s(iq,kdown)-s(iq,kb_sav(iq))
     .      +hl*(qsk-qs(iq,kb_sav(iq))))/(cp+hl*dqsdt(iq,kb_sav(iq)))
!      don't want moistening of cloud base layer (can lead to -ve fluxes)     
c      qdown=min( qq(iq,kb_sav(iq)) , qs(iq,kb_sav(iq))+
c    .           (tdown-t(iq,kb_sav(iq)))*dqsdt(iq,kb_sav(iq)) )
       qdown= qs(iq,kb_sav(iq))+
     .           (tdown-t(iq,kb_sav(iq)))*dqsdt(iq,kb_sav(iq)) 
       dprec=qdown-qsk                      ! to be mult by fldownn
       fldownn=fldown*(sig(kb_sav(iq))-sig(kt_sav(iq))) ! fldown=.66 say        
       totprec=qprec-fldownn*dprec 
       fldow(iq)=(.5+sign(.5,totprec))*fldownn ! suppr. downdraft for totprec<0
       rnrtcn(iq)=qprec-fldow(iq)*dprec        ! already has dsk factor
       if(ntest.eq.1.and.iq.eq.idjd)then
         print *,'kdown,tdown,qdown ',kdown,tdown,qdown
         print *,'qsk,qprec,dels0 ',qsk,qprec,hl*qprec 
         print *,'fldow,dprec,rnrtcn ',fldow(iq),dprec,rnrtcn(iq)
       endif
!      add in downdraft contributions
       dels(iq,kdown)=dels(iq,kdown)-fldow(iq)*s(iq,kdown)
       delq(iq,kdown)=delq(iq,kdown)-fldow(iq)*qsk
!      calculate emergent downdraft properties
       delq(iq,kb_sav(iq))=delq(iq,kb_sav(iq))+fldow(iq)*qdown
       dels(iq,kb_sav(iq))=dels(iq,kb_sav(iq))+fldow(iq)*
     .        (s(iq,kb_sav(iq))+cp*(tdown-t(iq,kb_sav(iq))))
!      add contrib to cloud base layer (unit flux this one)
       delq(iq,kb_sav(iq))=delq(iq,kb_sav(iq))-qbase(iq)
       dels(iq,kb_sav(iq))=dels(iq,kb_sav(iq))-sbase(iq)
      enddo  ! iq loop
      if((ntest.ne.0.or.diag).and.mydiag)then
       print *,"delsa ",(dels(idjd,k),k=1,kt_sav(idjd))
       print *,"delqa ",(delq(idjd,k),k=1,kt_sav(idjd))
      endif

!     subsidence and (possible) "shallow" detrainment
      do k=kuocb+1,kl-1
       do iq=1,ifull
        if(k.gt.kb_sav(iq).and.k.le.kt_sav(iq))then
c         savg=.5*(s(iq,k)+s(iq,k-1))
c         qavg=.5*(qq(iq,k)+qq(iq,k-1))   
          savg=s(iq,k)   ! 23/5/03
          qavg=qq(iq,k)  ! 23/5/03 to avoid giving -ve qg for the layer
          dels(iq,k)=dels(iq,k)-savg             ! subsidence 
          dels(iq,k-1)=dels(iq,k-1)+savg         ! subsidence into lower layer
          delq(iq,k)=delq(iq,k)-qavg             ! subsidence 
          delq(iq,k-1)=delq(iq,k-1)+qavg         ! subsidence into lower layer
        endif  ! (k.gt.kb_sav(iq).and.k.le.kt_sav(iq))
       enddo   ! iq loop
      enddo    ! k loop
      if((ntest.ne.0.or.diag).and.mydiag)then
       print *,"delsb ",(dels(idjd,k),k=1,kt_sav(idjd))
       print *,"delqb ",(delq(idjd,k),k=1,kt_sav(idjd))
      endif

!     modify calculated subsidence for downdraft effects
      do k=kuocb+1,kl-1
       do iq=1,ifull
        kdown=min(kl , kb_sav(iq)+nint(.6+.75*(kt_sav(iq)-kb_sav(iq))))
        if(k.le.kdown.and.k.gt.kb_sav(iq))then
c        savgb=.5*(s(iq,k)+s(iq,k-1))       
c        qavgb=.5*(qq(iq,k)+qq(iq,k-1))   
         savgb=s(iq,k)   ! 23/5/03
         qavgb=qq(iq,k)  ! 23/5/03
         dels(iq,k)=dels(iq,k)+fldow(iq)*savgb     ! anti-subsidence 
         dels(iq,k-1)=dels(iq,k-1)-fldow(iq)*savgb ! anti-subsidence into l-l
         delq(iq,k)=delq(iq,k)+fldow(iq)*qavgb     ! anti-subsidence 
         delq(iq,k-1)=delq(iq,k-1)-fldow(iq)*qavgb ! anti-subsidence into l-l
        endif
       enddo   ! iq loop
      enddo    ! k loop
      if((ntest.ne.0.or.diag).and.mydiag)then
       print *,"delsc ",(dels(idjd,k),k=1,kt_sav(idjd))
       print *,"delqc ",(delq(idjd,k),k=1,kt_sav(idjd))
      endif

!!!!!!!!! "shallow" detrainment using detrainx, dsig2, dsig4 !!!!!v3!!!!!! 
      if(detrainx.gt.0.)then   
       do iq=1,ifull  
        detr=min(detrainx,max( 0.,
     .  (dsig4-sig(kb_sav(iq))+sig(kt_sav(iq)))*detrainx/(dsig4-dsig2)))
        qxcess(iq)=detr*rnrtcn(iq)              ! e.g. .2* gives 20% detrainment
        rnrtcn(iq)=rnrtcn(iq)-qxcess(iq)
       enddo   ! iq loop
     
!      N.B. pre convpsav & dsk
       if(methdetr.eq.1)then                    ! moisten top layer only 
         do iq=1,ifull  
          k=kt_sav(iq)
          deltaq1=qxcess(iq)                    ! moistening
          delq(iq,k)=delq(iq,k)+deltaq1
          dels(iq,k)=dels(iq,k)-deltaq1*hl
         enddo  ! iq loop
       elseif(methdetr.eq.-1)then               ! moisten bottom layer only 
         do iq=1,ifull  
          k=kb_sav(iq)+1
          deltaq1=qxcess(iq)                    ! moistening
          delq(iq,k)=delq(iq,k)+deltaq1
          dels(iq,k)=dels(iq,k)-deltaq1*hl
         enddo  ! iq loop
       elseif(methdetr.eq.8)then        ! moisten all cloud layers, top most
         do iq=1,ifull  
          nlayers=kt_sav(iq)-kb_sav(iq)
          sum=.5*nlayers*(nlayers+1)
          do k=kb_sav(iq)+1,kt_sav(iq)
          deltaq8=qxcess(iq)*(k-kb_sav(iq))/sum ! moistening
          delq(iq,k)=delq(iq,k)+deltaq8
          dels(iq,k)=dels(iq,k)-deltaq8*hl
          enddo  
         enddo   ! iq loop
       elseif(methdetr.eq.-8)then ! moisten all cloud layers, bottom most
        do iq=1,ifull  
         nlayers=kt_sav(iq)-kb_sav(iq)
         sum=.5*nlayers*(nlayers+1)
         do k=kb_sav(iq)+1,kt_sav(iq)
	   deltaq8=qxcess(iq)*(kt_sav(iq)+1-k)/sum ! moistening
          delq(iq,k)=delq(iq,k)+deltaq8
          dels(iq,k)=dels(iq,k)-deltaq8*hl
         enddo  
        enddo   ! iq loop
       elseif(methdetr.eq.9)then                ! moisten all cloud layers
         do iq=1,ifull  
          nlayers=kt_sav(iq)-kb_sav(iq)
          do k=kb_sav(iq)+1,kt_sav(iq)
           deltaq9=qxcess(iq)/nlayers           ! moistening
           delq(iq,k)=delq(iq,k)+deltaq9
           dels(iq,k)=dels(iq,k)-deltaq9*hl
          enddo  
         enddo   ! iq loop
       endif     ! (methdetr.eq.1)  .. else ..
      endif  ! (detrainx.gt.0.)

      if((ntest.ne.0.or.diag).and.mydiag)then
        print *,"before diag print of dels,delh "
        iq=idjd
        nlayersd=max(1,nint((kt_sav(iq)-kb_sav(iq)-.1)/
     .           max(1,methdetr))) ! round down
        khalfd=kt_sav(iq)+1-nlayersd
        nlayersp=max(1,nint((kt_sav(iq)-kb_sav(iq)-.1)/methprec)) ! round down
        khalfp=kt_sav(iq)+1-nlayersp
        print *,'kbsav,ktsav,khalfd ',
     .           kb_sav(iq),kt_sav(iq),khalfd 
        print *,'nlayersd,khalfp,nlayersp ',
     .           nlayersd,khalfp,nlayersp 
        print *,"delsd ",(dels(idjd,k),k=1,kt_sav(idjd))
        print *,"delqd ",(delq(idjd,k),k=1,kt_sav(idjd))
        print *,'delh ',(dels(idjd,k)+hl*delq(idjd,k),k=1,kl)
        print *,'delhb',(dels(iq,k)+alfqarr(iq)*hl*delq(iq,k),k=1,kl)
        sum=0.
        do k=kb_sav(idjd),kt_sav(idjd)
         sum=sum+dels(idjd,k)+hl*delq(idjd,k)
        enddo
        print *,'qbase,sum_delh ',qbase(idjd),sum
      endif

!     calculate actual delq and dels
      do k=1,kl-1   
       do iq=1,ifull
        if(k.gt.kb_sav(iq))then
          delq(iq,k)=delq(iq,k)/dsk(k)
          dels(iq,k)=dels(iq,k)/dsk(k)
        else
          if(nbase.eq.0.and.k.ge.kbsav_b(iq))then
!           assume cloud base layer, and those below it, are mixed well
c           delq(iq,k)=delq(iq,kb_sav(iq))/(1.-sigmh(kb_sav(iq)+1))
c           dels(iq,k)=dels(iq,kb_sav(iq))/(1.-sigmh(kb_sav(iq)+1))
            delq(iq,k)=delq(iq,kb_sav(iq))/
     .                          (sigmh(kbsav_b(iq))-sigmh(kb_sav(iq)+1))
            dels(iq,k)=dels(iq,kb_sav(iq))/
     .                          (sigmh(kbsav_b(iq))-sigmh(kb_sav(iq)+1))
          elseif(k.eq.kb_sav(iq))then  ! s and qg only from kb layer
            delq(iq,k)=delq(iq,k)/dsk(k)
            dels(iq,k)=dels(iq,k)/dsk(k)
          endif  ! (nbase.eq.0)
        endif    ! (k.gt.kb_sav(iq)) .. else ..
       enddo     ! iq loop
      enddo      ! k loop

      if((ntest.ne.0.or.diag).and.mydiag)then
        print *,"before convpsav calc, after division by dsk"
        print *,"dels ",(dels(idjd,k),k=1,kt_sav(idjd))
        print *,"delq ",(delq(idjd,k),k=1,kt_sav(idjd))
      endif

!     calculate base mass flux 
      do iq=1,ifull
       if(kt_sav(iq).lt.kl)then 
!        Base limiter: new_qq(kb)=rhcv*new_qs(kb), i.e.
!        [qq+M*delq]_kb=rhcv*[qs+M*dqsdt*dels/cp]_kb 
         fluxb(iq)=max(0.,(rhcv*qs(iq,kb_sav(iq))-qq(iq,kb_sav(iq)))/ ! with dqsdt term
     .            ( delq(iq,kb_sav(iq))
     .             -rhcv*dqsdt(iq,kb_sav(iq))*dels(iq,kb_sav(iq))/cp ) ) 
         convpsav(iq)=fluxb(iq)
         if(nfluxq.eq.1)then  ! does little cf flux_dsk
!          fluxq limiter: new_qq(kb)=new_qq(kt)
!          i.e. [qq+M*delq]_kb=[qq+M*delq]_kt
           fluxq(iq)=max(0.,(qq(iq,kb_sav(iq))-qq(iq,kt_sav(iq)))/ 
     .                      (delq(iq,kt_sav(iq)) - delq(iq,kb_sav(iq))))
           convpsav(iq)=min(convpsav(iq),fluxq(iq))
         endif  ! (nfluxq.eq.1)
         if(nfluxdsk.eq.1)convpsav(iq)=
     .                     min(convpsav(iq),dsk(kb_sav(iq)))
         if(nfluxdsk.eq.2)convpsav(iq)=
     .                     min(convpsav(iq),.5*dsk(kb_sav(iq)))
         if(nfluxdsk.eq.-2)convpsav(iq)=
     .                    min(convpsav(iq),dsk(kb_sav(iq))/(2*iterconv))
       endif    ! (kt_sav(iq).lt.kl)
      enddo     ! iq loop
      
      do k=kuocb+1,kl-1
       do iq=1,ifull
!       want: new_hbas>=new_hs(k), i.e. in limiting case:
!       [h+alfsarr*M*dels+alfqarr*M*hl*delq]_base 
!                                          = [hs+M*dels+M*hlcp*dels*dqsdt]_k
        fluxt(iq,k)=(sbase(iq)+hl*qbase(iq)-hs(iq,k))/
     .              (dels(iq,k)*(1.+hlcp*dqsdt(iq,k))
     .              -alfsarr(iq)*dels(iq,kb_sav(iq))
     .              -alfqarr(iq)*hl*delq(iq,kb_sav(iq)) ) 
c        if(k.eq.kt_sav(iq))then
c          convpsav(iq)=min(convpsav(iq),max(0.,fluxt(iq,k)))
c        endif   ! (k.eq.kt_sav(iq))
c        if(k.gt.kb_sav(iq).and.k.lt.kt_sav(iq).and.
        if(k.gt.kb_sav(iq).and.k.le.kt_sav(iq).and.
     .     fluxt(iq,k).gt.0.)then
           convpsav(iq)=min(convpsav(iq),fluxt(iq,k))
        endif   ! (k.gt.kb_sav(iq).and.k.lt.kt_sav(iq).and.fluxt(iq,k).gt.0.)
       enddo    ! iq loop
      enddo     ! k loop

      if(abs(ksc).gt.90)then  ! will do shallow clouds in vertmix
        do iq=1,ifull
c        if(sig(kt_sav(iq)).gt.sigksct)then
         if(sig(kt_sav(iq)).gt..8)then
           convpsav(iq)=0.
         endif
        enddo  ! iq loop
      elseif(ksc.eq.22)then  ! will do shallow clouds in vertmix
        do iq=1,ifull
         if(sig(kt_sav(iq)).gt..01*ksc)then  ! e.g. ksc=80
           convpsav(iq)=0.
         endif
        enddo  ! iq loop
      endif
      
      if(itn.eq.1)then 
        kbsav(:)=kb_sav(:) 
        ktsav(:)=kt_sav(:) 
      endif                              ! (itn.eq.1)
      if(itn.lt.iterconv)then 
        convpsav(:)=convfact*convpsav(:) ! typically convfact=1.02  
      endif                              ! (itn.lt.iterconv)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      if((ntest.ne.0.or.diag).and.mydiag)then
        iq=idjd
        detr=max( detrain,detrain+(dsig4-sig(kb_sav(iq))+
     .            sig(kt_sav(iq)))*(detrainx-detrain)/(dsig4-dsig2) )
        detr=min(detr,detrainx)
        print *,'detr,delta_sig,qxcess ',
     .           detr,sig(kb_sav(iq))-sig(kt_sav(iq)),qxcess(iq)
        flux_dsk=dsk(kb_sav(idjd))*convfact
        if(nfluxdsk.eq.2)flux_dsk=.5*dsk(kb_sav(idjd))*convfact
        write(6,"('flux_dsk,fluxq,fluxb,convpsav',5f8.5)")
     .       flux_dsk,fluxq(idjd),fluxb(idjd),convpsav(idjd)
        write(6,"('fluxt',9f8.5/9f8.5)") (fluxt(iq,k),k=1,kt_sav(iq))
        write(6,"('delQ',18f6.3)")
     .   (1000.*convpsav(idjd)*delq(idjd,k),k=1,kl)
        write(6,"('delt',18f6.3)")
     .   (convpsav(idjd)*dels(idjd,k)/cp,k=1,kl)
        write(6,"('delQ*dsk',18f6.3)")
     .   (1000.*convpsav(idjd)*delq(idjd,k)*dsk(k),k=1,kl)
        write(6,"('delt*dsk',18f6.3)")
     .   (convpsav(idjd)*dels(idjd,k)*dsk(k)/cp,k=1,kl)
        convmax=0.
        nums=0
        do iq=1,ifull     
         if(kt_sav(iq)-kb_sav(iq).eq.1)then
           nums=nums+1
           if(convpsav(iq).gt.convmax)then
c          if(nums.lt.20)then
             convmax=convpsav(iq)
             write(6,"('bc  ',2i5,2i3,2x,3f6.3,f7.3,f6.3,' rh:',2f5.2)") 
     .       iq,nums,kb_sav(iq),kt_sav(iq),
     .       dsk(kb_sav(iq)),fluxq(iq),fluxb(iq),
     .       fluxt(iq,kt_sav(iq)),convpsav(iq),
     .       qq(iq,kb_sav(iq))/qs(iq,kb_sav(iq)),
     .       qq(iq,kt_sav(iq))/qs(iq,kt_sav(iq))
           endif
         endif
        enddo  ! iq loop
        convmax=0.
        nums=0
        do iq=1,ifull     
         if(kt_sav(iq)-kb_sav(iq).eq.2)then
           nums=nums+1
           if(convpsav(iq).gt.convmax)then
c          if(nums.lt.20)then
             convmax=convpsav(iq)
             write(6,"('bcd ',2i5,2i3,2x,3f6.3,f7.3,f6.3,' rh:',2f5.2)") 
     .       iq,nums,kb_sav(iq),kt_sav(iq),
     .       dsk(kb_sav(iq)),fluxq(iq),fluxb(iq),
     .       fluxt(iq,kt_sav(iq)),convpsav(iq),
     .       qq(iq,kb_sav(iq))/qs(iq,kb_sav(iq)),
     .       qq(iq,kt_sav(iq))/qs(iq,kt_sav(iq))
           endif
         endif
        enddo  ! iq loop
        convmax=0.
        nums=0
        do iq=1,ifull     
         if(kt_sav(iq)-kb_sav(iq).eq.3)then
           nums=nums+1
           if(convpsav(iq).gt.convmax)then
c          if(nums.lt.20)then
             convmax=convpsav(iq)
             write(6,"('bcde',2i5,2i3,2x,3f6.3,f7.3,f6.3,' rh:',2f5.2)") 
     .       iq,nums,kb_sav(iq),kt_sav(iq),
     .       dsk(kb_sav(iq)),fluxq(iq),fluxb(iq),
     .       fluxt(iq,kt_sav(iq)),convpsav(iq),
     .       qq(iq,kb_sav(iq))/qs(iq,kb_sav(iq)),
     .       qq(iq,kt_sav(iq))/qs(iq,kt_sav(iq))
           endif
         endif
        enddo  ! iq loop
      endif    ! (ntest.ne.0.or.diag)

!     update qq, tt and precip
      rnrtc(:)=rnrtc(:)+convpsav(:)*rnrtcn(:)*conrev(:) ! g/m**2/s
      do k=1,kl   
       do iq=1,ifull
        qq(iq,k)=qq(iq,k)+convpsav(iq)*delq(iq,k)
        tt(iq,k)=tt(iq,k)+convpsav(iq)*dels(iq,k)/cp
       enddo    ! iq loop
      enddo     ! k loop
      
!!!!!!!!!!!!!!!!!!!!! "deep" detrainment using detrain !!!!!v3!!!!!!    
!     N.B. convpsav has been updated here, but without factr term  
      qxcess(:)=detrain*rnrtcn(:)             ! e.g. .2* gives 20% detrainment
      rnrtcn(:)=rnrtcn(:)-qxcess(:)
      if(methprec.eq.1)then                   ! moisten top layer only 
       do iq=1,ifull  
        deltaq1=convpsav(iq)*qxcess(iq)/dsk(kt_sav(iq))             
        qliqw(iq,kt_sav(iq))=qliqw(iq,kt_sav(iq))+deltaq1
       enddo  ! iq loop
      elseif(methprec.eq.8)then               ! moisten all cloud layers, top most
       do iq=1,ifull  
        nlayers=kt_sav(iq)-kb_sav(iq)
        sum=.5*nlayers*(nlayers+1)
        do k=kb_sav(iq)+1,kt_sav(iq)
           deltaq8=convpsav(iq)*qxcess(iq)*(k-kb_sav(iq))/(sum*dsk(k)) 
         qliqw(iq,k)=qliqw(iq,k)+deltaq8
        enddo  
       enddo   ! iq loop
      elseif(methprec.eq.9)then               ! moisten all cloud layers
       do iq=1,ifull  
        nlayers=kt_sav(iq)-kb_sav(iq)
        do k=kb_sav(iq)+1,kt_sav(iq)
         deltaq9=convpsav(iq)*qxcess(iq)/(nlayers*dsk(k))       
         qliqw(iq,k)=qliqw(iq,k)+deltaq9
        enddo  
       enddo   ! iq loop
      else                                    ! moisten upper layers
       do k=kuocb+1,kl-1 
        do iq=1,ifull  
         nlayers=max(1,nint((kt_sav(iq)-kb_sav(iq)-.1)/methdetr))   ! round down
         khalf=kt_sav(iq)+1-nlayers
         if(k.ge.khalf.and.k.le.kt_sav(iq))then
           deltaq=convpsav(iq)*qxcess(iq)/(nlayers*dsk(k))         
           qliqw(iq,k)=qliqw(iq,k)+deltaq
         endif
        enddo  ! iq loop
       enddo   ! k loop
      endif    ! (methprec.eq.1)  .. else ..

      if((ntest.ne.0.or.diag).and.mydiag)then
        iq=idjd
!       N.B. convpsav(iq) is already mult by dt jlm: mass flux is convpsav/dt
        print *,"after convection: ktau,itn,kbsav,ktsav ",
     .                 ktau,itn,kb_sav(idjd),kt_sav(idjd)
        write (6,"('qgc ',19f7.3/(8x,19f7.3))") 
     .             (1000.*qq(idjd,k),k=1,kl)
        write (6,"('ttc ',19f7.2/(8x,19f7.2))") 
     .             (tt(idjd,k),k=1,kl)
        print *,'rnrtc,fldow,qxcess ',rnrtc(iq),fldow(iq),qxcess(iq)
        print *,'ktsav,qbase,qs_ktsav,qq_ktsav ',
     .          kt_sav(iq),qbase(iq),qs(iq,kt_sav(iq)),qq(iq,kt_sav(iq)) 
        iq=idjd
        delq_av=0.
        delt_av=0.
        rnd_av=rnd_av+rnrtc(iq)
        heatlev=0.
        do k=1,kl
         delq_av=delq_av+dsk(k)*convpsav(iq)*delq(iq,k)
         delt_av=delt_av+dsk(k)*convpsav(iq)*dels(iq,k)/cp
         heatlev=heatlev+sig(k)*dsk(k)*convpsav(iq)*dels(iq,k)/cp
        enddo
        print *,'delq_av,delt_exp,rnd_av,rnd_exp ',
     .         delq_av,-delq_av*hl/cp,rnd_av,-delq_av*conrev(iq)
        if(delt_av.ne.0.)print *,'ktau,itn,kbsav,ktsav,delt_av,heatlev',
     .        ktau,itn,kb_sav(idjd),kt_sav(idjd),delt_av,heatlev/delt_av
      endif   ! (ntest.ne.0.or.diag)

      if(nuvconv.gt.0)then   !  momentum calculations, similar to above
        do k=1,kl   
         do iq=1,ifull
          delu(iq,k)=0.
          delv(iq,k)=0.
         enddo   ! iq loop
        enddo    ! k loop
        do iq=1,ifull
         k=kt_sav(iq)                        ! dels and delq for top layer
         delu(iq,k)=delu(iq,k)+(u(iq,kb_sav(iq))-u(iq,k)) 
         delu(iq,k-1)=u(iq,k)                ! subsidence into lower layer
         delv(iq,k)=delv(iq,k)+(v(iq,kb_sav(iq))-v(iq,k))
         delv(iq,k-1)=v(iq,k)                ! subsidence into lower layer
         do k=kt_sav(iq)-1,kb_sav(iq)+1,-1
          delu(iq,k)=delu(iq,k)-u(iq,k)      ! subs.
          delu(iq,k-1)=u(iq,k)               ! subsidence into lower layer
          delv(iq,k)=delv(iq,k)-v(iq,k)      ! subs.
          delv(iq,k-1)=v(iq,k)               ! subsidence into lower layer
         enddo
!        assume cloud base layer, and those below it, are mixed well
         delu(iq,kb_sav(iq))=delu(iq,kb_sav(iq))-u(iq,kb_sav(iq))
         delv(iq,kb_sav(iq))=delv(iq,kb_sav(iq))-v(iq,kb_sav(iq))
        enddo  ! iq loop
!       calculate actual delu and delv
        do k=1,kl-1   
         do iq=1,ifull
          if(k.gt.kb_sav(iq))then
            delu(iq,k)=delu(iq,k)/dsk(k)
            delv(iq,k)=delv(iq,k)/dsk(k)
          else
            if(nbase.eq.0.and.k.ge.kbsav_b(iq))then
!             assume cloud base layer, and those below it, are mixed well
              delu(iq,k)=delu(iq,kb_sav(iq))/
     .                   (sigmh(kbsav_b(iq))-sigmh(kb_sav(iq)+1))
              delv(iq,k)=delv(iq,kb_sav(iq))/
     .                   (sigmh(kbsav_b(iq))-sigmh(kb_sav(iq)+1))
            elseif(k.eq.kb_sav(iq))then        ! u and v only from kb layer
              delu(iq,k)=delu(iq,k)/dsk(k)
              delv(iq,k)=delv(iq,k)/dsk(k)
            endif  ! (nbase.eq.0)
          endif    ! (k.gt.kb_sav(iq)) .. else ..
         enddo     ! iq loop
        enddo      ! k loop
!       update u, v
        do k=1,kl   
         do iq=1,ifull
          u(iq,k)=u(iq,k)+factr*convpsav(iq)*delu(iq,k)
          v(iq,k)=v(iq,k)+factr*convpsav(iq)*delv(iq,k)
         enddo    ! iq loop
        enddo     ! k loop
      endif     ! (nuvconv.gt.0)

!     section for convective transport of trace gases (jlm 22/2/01)
      if(ilt.gt.1)then
!       if(iterconv.ne.1)stop 'need 1 for trace gases' ! should be OK now
        do ntr=1,ntrac
         do k=1,kl
           do iq=1,ifull
           s(iq,k)=tr(iq,k,ntr)
          enddo    ! iq loop
         enddo     ! k loop
         do iq=1,ifull
          if(kt_sav(iq).lt.kl)then
            kb=kb_sav(iq)
            kt=kt_sav(iq)
            veldt=factr*convpsav(iq)*(1.-fldow(iq))  ! simple treatment
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
      
      if((ntest.ne.0.or.diag).and.mydiag)then
        write (6,"('uuc ',19f6.1/(8x,19f6.1))") (u(idjd,k),k=1,kl)
        write (6,"('vvc ',19f6.1/(8x,19f6.1))") (v(idjd,k),k=1,kl)
      endif
      enddo     ! itn=1,iterconv

      if(factr.lt.1.)then
        rnrtc(:)=factr*rnrtc(:)
        qq(1:ifull,:)=qg(1:ifull,:)+factr*(qq(1:ifull,:)-qg(1:ifull,:))      
        qliqw(1:ifull,:)=factr*qliqw(1:ifull,:)      
        tt(1:ifull,:)= t(1:ifull,:)+factr*(tt(1:ifull,:)- t(1:ifull,:))      
      endif  ! (factr.lt.1.)

!     update qq, tt for evap of qliqw (qliqw arose from moistening)
      if(ifullw.eq.ifull)then
!       Leon's stuff here, e.g.
        do k=1,kl
          do iq=1,ifullw
           qlg(iq,k)=qlg(iq,k)+qliqw(iq,k)
          enddo
         enddo
      else
        qq(1:ifull,:)=qq(1:ifull,:)+qliqw(1:ifull,:)         
        tt(1:ifull,:)=tt(1:ifull,:)-hl*qliqw(1:ifull,:)/cp   
      endif  ! (ifullw.eq.ifull)
!__________________________end of convective calculations_____________________
     
      if((ntest.ne.0.or.diag).and.mydiag)then
        print *,"after convection"
        write (6,"('qge ',19f7.3/(8x,19f7.3))")(1000.*qq(idjd,k),k=1,kl)
        write (6,"('tte ',19f6.1/(8x,19f6.1))")(tt(idjd,k),k=1,kl)
        print *,'rnrtc ',rnrtc(idjd)
!       following calc of total heatlev really needs iterconv=1  
         do k=kl-2,1,-1
          delt_av=delt_av-dsk(k)*revc(idjd,k)*hlcp
          heatlev=heatlev-sig(k)*dsk(k)*revc(idjd,k)*hlcp
          print *,'k,rh,delt_av, heatlev ',
     .           k,100.*qq(idjd,k)/qs(idjd,k),delt_av,heatlev
         enddo
         if(delt_av.ne.0.)print *,'ktau,delt_av-net,heatlev_net ',
     .                             ktau,delt_av,heatlev/delt_av
      endif
      if(ifullw.gt.1)go to 8

!__________________________beginning of large-scale calculations_____________________
!     check for grid-scale rainfall 
5     do k=1,kl   
       do iq=1,ifull
        es(iq,k)=establ(tt(iq,k))
       enddo  ! iq loop
      enddo   ! k loop
      do k=kl,1,-1    ! top down to end up with proper kbsav_ls
       do iq=1,ifull
        pk=ps(iq)*sig(k)
        qs(iq,k)=max(.622*es(iq,k)/(pk-es(iq,k)),1.5e-6)  
        if(qq(iq,k).gt.rhsat*qs(iq,k))then
          kbsav_ls(iq)=k
          gam=max(hlcp*qs(iq,k)*pk*hlars/(tt(iq,k)**2*(pk-es(iq,k))),0.) 
          dqrx=(qq(iq,k)-rhsat*qs(iq,k))/(1.+rhsat*gam)
          tt(iq,k)=tt(iq,k)+hlcp*dqrx
          qq(iq,k)=qq(iq,k)-dqrx
          rnrt(iq)=rnrt(iq)+dqrx*dsk(k)
        endif   ! (qq(iq,k).gt.rhsat*qs(iq,k))
       enddo    ! iq loop
      enddo     ! k loop

!!!!!!!!!!!!!!!!  now do evaporation of L/S precip !!!!!!!!!!!!!!!!!!
!     conrev(iq)=1000.*ps(iq)/(grav*dt)     
      rnrt(:)=rnrt(:)*conrev(:)                 
!     here the rainfall rate rnrt has been converted to g/m**2/sec
      fluxr(:)=rnrt(:)*1.e-3*dt ! kg/m2      

      if((ntest.ne.0.or.diag).and.mydiag)then
        print *,'after large scale rain: kbsav_ls,rnrt ',
     .                                   kbsav_ls(idjd),rnrt(idjd)
        write (6,"('qg ',19f7.3/(8x,19f7.3))") 
     .             (1000.*qq(idjd,k),k=1,kl)
        write (6,"('tt ',19f7.2/(8x,19f7.2))") 
     .             (tt(idjd,k),k=1,kl)
      endif

      if(nevapls.eq.5)then ! even newer UKMO Thu  03-05-1998
        rKa=2.4e-2
        Dva=2.21
        cfls=1. ! cld frac7 large scale
        cflscon=4560.*cfls**.3125
        do k=2*kl/3,1,-1
         do iq=1,ifull
          if(k.lt.kbsav_ls(iq))then
            rhodz=ps(iq)*dsk(k)/grav
            qpf=fluxr(iq)/rhodz     ! Mix ratio of rain which falls into layer
            pk=ps(iq)*sig(k)
c           es(iq,k)=qs(iq,k)*pk/.622
            Apr=hl*hl/(rKa*tt(iq,k)*(rvap*tt(iq,k))-1.)
            Bpr=rvap*tt(iq,k)*pk/(Dva*es(iq,k))
            Fr=fluxr(iq)/(cfls*dt)
            rhoa=pk/(rdry*tt(iq,k))
            dz=pk/(rhoa*grav)
            Vr=max(.01 , 11.3*Fr**(1./9.)/sqrt(rhoa)) ! Actual fall speed
            dtev=dz/Vr
            qr=fluxr(iq)/(dt*rhoa*Vr)
            qgdiff=qs(iq,k)-qq(iq,k)
            Cev2=cflscon*qgdiff/(qs(iq,k)*(Apr+Bpr))  ! Ignore rhoa**0.12
            qr2=max(0. , qr**.3125 - .3125*Cev2*dtev)**3.2
!           Cev=(qr-qr2)/(dtev*qgdiff)
            Cevx=(qr-qr2)/dtev  ! i.e. Cev*qgdiff
            alphal=hl*qs(iq,k)/(ars*tt(iq,k)**2)
!           bl=1.+Cev*dt*(1.+hlcp*alphal)
            blx=qgdiff+Cevx*dt*(1.+hlcp*alphal)  ! i.e. bl*qgdiff
!           evapls= cfls*dt*qgdiff*Cev/bl        ! UKMO
            evapls= cfls*dt*Cevx*qgdiff/blx      ! UKMO
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
          endif   ! (k.lt.kbsav_ls(iq))
         enddo    ! iq loop
        enddo     ! k loop
      endif       ! (nevapls.eq.5)
!__________________________end of large-scale calculations_____________________


8     qg(1:ifull,:)=qq(1:ifull,:)                   
      condc(:)=.001*dt*rnrtc(:)      ! convective precip for this timestep
      precc(:)=precc(:)+condc(:)        
      condx(:)=condc(:)+.001*dt*rnrt(:) ! total precip for this timestep
      precip(:)=precip(:)+condx(:)      
      t(1:ifull,:)=tt(1:ifull,:)             

      if((ntest.eq.1.or.diag))then
        if(mydiag) then
        print *,'at end of convjlm: rnrt,rnrtc ',
     .                              rnrt(idjd),rnrtc(idjd)
        write (6,"('qg ',19f7.3/(8x,19f7.3))") 
     .             (1000.*qq(idjd,k),k=1,kl)
        write (6,"('tt ',19f7.2/(8x,19f7.2))") 
     .             (tt(idjd,k),k=1,kl)
        pwater=0.   ! in mm     
          do k=1,kl
         pwater=pwater-dsig(k)*qg(idjd,k)*ps(idjd)/grav
          enddo
          print *,'pwater0,pwater+condx,pwater ',
     .           pwater0,pwater+condx(idjd),pwater
        print *,'D rnrt,rnrtc,condx ',
     .             rnrt(idjd),rnrtc(idjd),condx(idjd)
        print *,'precc,precip ',
     .           precc(idjd),precip(idjd)
        endif
        call maxmin(rnrtc,'rc',ktau,1.,1)
      endif
      return
      end
