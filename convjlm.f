      subroutine convjlm      ! jlm convective scheme - Version v4
!     the shallow convection options here just for iterconv=1
!     has +ve fldownn depending on delta sigma; -ve uses older abs(fldown)   
!     N.B. nevapcc option has been removed   
      use cc_mpi, only : mydiag
      use diag_m
      implicit none
      integer itn,iq,k,k13,k23,kb,kbsav_b,kbsav_ls,kb_sav,kcl_top
     .       ,kdown,khalf,khalfd,khalfp,kt,kt_sav,ktmax
     .       ,nalfs,nbase,nfluxq,nfluxdsk,nlayers,nlayersd,nlayersp
     .       ,ntest,ntr,nums
      real Aev,Asb,alfcon,alfq,alfs,alfqarr,alfsarr
     .    ,conrev,convmax,delq,dels,deltaq,delq_av,delt_av
     .    ,detr,dprec,dqrx,dqsdt,es
     .    ,factdav,factr,factuv,fldow,fldownn
     .    ,fluxq,fluxr,fluxb,flux_dsk,fluxup,fluxt
     .    ,frac,fraca,fracb,gam,hbase,heatlev,hs
     .    ,pwater,pwater0,qavg,qavgb,qbas,qbase,qdown,qliqw,qprec
     .    ,qq,qs,qsk,qxcess
     .    ,rkmid,rnd_av,rnrt,rnrtc,rnrtcn
     .    ,s,savg,savgb,sbas,sbase,sum,sumqs,sumqt,sumss
     .    ,tm,totprec,tt,tdown,veldt
     
     .    ,rKa,Dva,cfls,cflscon,rhodz,qpf,pk,Apr,Bpr,Fr,rhoa,dz,Vr
     .    ,dtev,qr,qgdiff,Cev2,qr2,Cevx,alphal,blx,evapls,revq

      include 'newmpar.h'
      parameter (ntest=0)      ! 1 or 2 to turn on
!     parameter (iterconv=1)   ! to kuocom.h
!     parameter (fldown=.6)    ! to kuocom.h
!     parameter (detrain=.05)  ! to kuocom.h
      parameter (alfcon=.015)  ! e.g. .01 or .015   or 1.1 (over sea)
!     parameter (alflnd=1.15)  ! e.g. 1.15  
!     parameter (alfsea=1.05)  ! e.g. 1.05  
!     parameter (nqbot=3)      ! 0: original qq(1), 1: qs1, 2: qbas for land/sea
!                                3:for RH-enhanced base
      parameter (nalfs=0)      ! 1: was usual, 0: off
      parameter (nfluxq=0)     ! 1 from 6/8/02  off from 16/10/03
      parameter (nfluxdsk=-2)  ! 1 till 6/8/02
      parameter (nbase=1)      ! 0 shared below, 1 uses only one (originally 0)
!     parameter (methdetr=2)   ! various - applies only to shallow clouds
!     parameter (methprec=8)   ! 1 (top only); 2 (top half); 4 top quarter (= kbconv)
!     parameter (nuvconv=0)    ! usually 0, >0 to turn on momentum mixing
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
     
      common/work2/ktmax(ifull),conrev(ifull),fluxr(ifull)
     .  ,fluxb(ifull),kbsav_b(ifull),kbsav_ls(ifull)
     .  ,rnrt(ifull),rnrtc(ifull),sumqs(ifull),sumqt(ifull)
     .  ,sumss(ifull),factdav(ifull),qbase(ifull),fldow(ifull)
     .  ,qxcess(ifull),fluxq(ifull),tdown(ifull),qdown(ifull)
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
      real h0(kl),q0(kl),t0(kl)  ! for idjd diags
      equivalence (phi,s),(delu,revc,delq),(delv,dels)
!     data rhcv/.75/                  now in kuocom
      real convshal(ifull)
      save convshal
      data convshal/ifull*0./

      include 'establ.h'
      Aev(tm) = 2.008e-9*tm**2 - 1.385e-6*tm + 2.424e-4  !For UKMO evap scheme
      Asb(tm) = max (-5.2e-9*tm**2+2.5332e-6*tm-2.9111e-4,1.e-5) !For UKMO subl

      do k=1,kl
       dsk(k)=-dsig(k)    !   dsk = delta sigma (positive)
       alfm(k)=min(.5,alfcon/(max(1.e-10,sig(1)-sig(k))))  ! e.g. alfcon=.01
       alfss(k)=1.-alfm(k)
       if(nalfs.eq.0)alfss(k)=1.  ! for no alf effect  - usual setting
      enddo     ! k loop
      if(alflnd.lt.1.)then     
        if(alflnd.gt.0.)then  ! common alfsea, alflnd method
          do k=1,kl
           alfqq_s(k)=1.-alfsea*sig(k)**8   ! e.g. alfsea=.5; =0. for off
           alfqq_l(k)=1.-alflnd*sig(k)**8   ! e.g. alflnd=.4; =0. for off
          enddo     ! k loop
	 else   ! this is alflnd < 0, i.e. oldest alfcon method
          do k=1,kl
           alfqq_s(k)=1.-alfm(k)
           alfqq_l(k)=1.-alfm(k)
          enddo     ! k loop
        endif  ! (alflnd.gt.0.)  ..else
        do k=1,kl
         alfqq1_s(k)=1.-alfqq_s(k)
         alfqq1_l(k)=1.-alfqq_l(k)
        enddo     ! k loop
      else   ! usual with alflnd>1 (was called nqbot=3 method)
        do k=1,kl
         alfqq_s(k)=abs(alfsea)
         alfqq_l(k)=alflnd
         alfqq1_s(k)=0.
         alfqq1_l(k)=0.
        enddo     ! k loop
      endif  !   (alflnd.lt.1.)  .. else ..
  
c     convective first, then L/S rainfall
      qliqw(:,:)=0.  
      conrev(1:ifull)=1000.*ps(1:ifull)/(grav*dt) ! factor to convert precip to g/m2/s
      rnrt(:)=0.       ! initialize large-scale rainfall array
      rnrtc(:)=0.      ! initialize convective  rainfall array
      ktmax(:)=kl      ! preset 1 level above current topmost-permitted ktsav
      kbsav_ls(:)=0    ! for L/S
      ktsav(:)=kl      ! preset value to show no deep or shallow convection
      kbsav(:)=kl      ! preset value to show no deep or shallow convection

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
       print *,'itn,iterconv,nbase,nuvconv ',
     .          itn,iterconv,nbase,nuvconv
       print *,'methdetr,detrain,detrainx,dsig2,dsig4 ',
     .          methdetr,detrain,detrainx,dsig2,dsig4
       print *,'methprec,fldown ',
     .          methprec,fldown
       print *,'nfluxq,nfluxdsk,nalfs ',
     .          nfluxq,nfluxdsk,nalfs
       print *,'alfcon,alflnd,alfsea ',alfcon,alflnd,alfsea
       write (6,"('alfqq_l ',12f7.3/(8x,12f7.3))") alfqq_l   ! max is 18f7.3
       write (6,"('alfqq_s ',12f7.3/(8x,12f7.3))") alfqq_s
       write (6,"('alfss   ',12f7.3/(8x,12f7.3))") alfss
      endif
      if( (ntest.ne.0.or.diag.or.(ktau.eq.1.and.nmaxpr.eq.1))
     &     .and. mydiag ) then
        print *,'near beginning of convjlm loop; ktau,itn,id,jd: ',
     .                                           ktau,itn,id,jd
        write (6,"('rh   ',12f7.2/(5x,12f7.2))") 
     .             (100.*qq(idjd,k)/qs(idjd,k),k=1,kl)
        write (6,"('qs   ',12f7.3/(5x,12f7.3))") 
     .             (1000.*qs(idjd,k),k=1,kl)
        write (6,"('qq   ',12f7.3/(5x,12f7.3))") 
     .             (1000.*qq(idjd,k),k=1,kl)
        pwater0=0.   ! in mm     
        do k=1,kl
         iq=idjd
	  h0(k)=s(idjd,k)/cp+hlcp*qq(idjd,k)
	  q0(k)=qq(idjd,k)
	  t0(k)=tt(idjd,k)
         pwater0=pwater0-dsig(k)*qg(idjd,k)*ps(idjd)/grav
         if(land(iq))then
           qbass(k)=alfqq_l(k)*qq(iq,k)+alfqq1_l(k)*qq(iq,1) 
         else
           qbass(k)=alfqq_s(k)*qq(iq,k)+alfqq1_s(k)*qq(iq,1) 
         endif
	  qbass(k)=min(qbass(k),max(qq(iq,k),qs(iq,k)))  ! added 3/6/03
        enddo
        write (6,"('qbas ',12f7.3/(5x,12f7.3))") 
     .             (1000.*qbass(k),k=1,kl)
        write (6,"('tt   ',12f7.2/(5x,12f7.2))") 
     .             (tt(idjd,k),k=1,kl)
        write (6,"('s/cp ',12f7.2/(5x,12f7.2))") (s(idjd,k)/cp,k=1,kl)
        write (6,"('h/cp ',12f7.2/(5x,12f7.2))") 
     .             (s(idjd,k)/cp+hlcp*qq(idjd,k),k=1,kl)
        write (6,"('hb/cp',12f7.2/(5x,12f7.2))") 
     .             (s(idjd,k)/cp+hlcp*qbass(k),k=1,kl)
        write (6,"('hs/cp',12f7.2/(5x,12f7.2))") 
     .             (hs(idjd,k)/cp,k=1,kl)
        write (6,"('k  ',12i7/(3x,12i7))") (k,k=1,kl)
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
c         if(kbsav_b(iq).eq.k.and.sig(kb_sav(iq)).gt..8)kbsav_b(iq)=k-1 
          if(kbsav_b(iq).eq.k)kbsav_b(iq)=k-1 
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
c         kbasep=kb_sav(iq)+1
c	   if((u(iq,k-1)-u(iq,kbasep))**2+(v(iq,k-1)-v(iq,kbasep))**2.
c     .        gt.convshr**2)kt_sav(iq)=k-1   ! possible shear test2 Feb '04
        endif   ! (kt_sav(iq).eq.k-1.and.k.lt.ktmax(iq))
       enddo    ! iq loop
      enddo     ! k loop
      
      if(itn.gt.1)then  ! added April 04
        do iq=1,ifull
	  if(rnrtc(iq).gt.0..and.kb_sav(iq).eq.kbsav(iq).
     .                      and.kt_sav(iq).eq.ktsav(iq))then
	    kt_sav(iq)=max(kbsav(iq)+1,ktsav(iq)-1)
	  endif   !  (rnrtc(iq).gt.0.....)
	 enddo
      endif   ! (itn.gt.1)

c      if(convrh.gt.0.)then
c        do k=kuocb+2,kl   ! upwards to find legal cloud top
c         do iq=1,ifull
c          if(kt_sav(iq).lt.kl.and.k.lt.kt_sav(iq))then 
c!            use lowest such layer, e.g. with convrh=.4
c  	     if(qq(iq,k)/qs(iq,k).lt.convrh)kt_sav(iq)=k  
c          endif
c         enddo    ! iq loop
c        enddo     ! k loop
c      endif       ! (convrh.gt.0.)

      if((ntest.ne.0.or.diag).and.mydiag)then
       iq=idjd
       kdown=min(kl,kb_sav(iq)+nint(.6+.75*(kt_sav(iq)-kb_sav(iq))))
       print *,"ktau,itn,kbsav_b,kb_sav,kt_sav,kdown ",
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
       dels(iq,kt_sav(iq))=dels(iq,kt_sav(iq))+hl*qprec    ! precip. heating
       dels(iq,kt_sav(iq))=dels(iq,kt_sav(iq))+sbase(iq)   ! s flux
       delq(iq,kt_sav(iq))=delq(iq,kt_sav(iq))+qsk  

       kdown=min(kl , kb_sav(iq)+nint(.6+.75*(kt_sav(iq)-kb_sav(iq))))
       qsk=min(qs(iq,kdown),qq(iq,kdown))    ! N.B. min for downdraft subs
       tdown(iq)=t(iq,kb_sav(iq))+(s(iq,kdown)-s(iq,kb_sav(iq))
     .      +hl*(qsk-qs(iq,kb_sav(iq))))/(cp+hl*dqsdt(iq,kb_sav(iq)))
!      don't want moistening of cloud base layer (can lead to -ve fluxes)     
c      qdown(iq)=min( qq(iq,kb_sav(iq)) , qs(iq,kb_sav(iq))+
c    .           (tdown(iq)-t(iq,kb_sav(iq)))*dqsdt(iq,kb_sav(iq)) )
       qdown(iq)= qs(iq,kb_sav(iq))+
     .           (tdown(iq)-t(iq,kb_sav(iq)))*dqsdt(iq,kb_sav(iq)) 
       dprec=qdown(iq)-qsk                      ! to be mult by fldownn
!	typically fldown=.6 or -.2        
       fldownn=max(-fldown,fldown*(sig(kb_sav(iq))-sig(kt_sav(iq))))
       totprec=qprec-fldownn*dprec 
       fldow(iq)=(.5+sign(.5,totprec))*fldownn ! suppr. downdraft for totprec<0
       rnrtcn(iq)=qprec-fldow(iq)*dprec        ! already has dsk factor
       if(ntest.eq.1.and.iq.eq.idjd)then
         print *,'kdown,tdown,qdown ',kdown,tdown(iq),qdown(iq)
         print *,'qsk,qprec,dels0 ',qsk,qprec,hl*qprec 
         print *,'fldow,dprec,rnrtcn ',fldow(iq),dprec,rnrtcn(iq)
       endif
!      add in downdraft contributions
       dels(iq,kdown)=dels(iq,kdown)-fldow(iq)*s(iq,kdown)
       delq(iq,kdown)=delq(iq,kdown)-fldow(iq)*qsk
!      calculate emergent downdraft properties
       delq(iq,kb_sav(iq))=delq(iq,kb_sav(iq))+fldow(iq)*qdown(iq)
       dels(iq,kb_sav(iq))=dels(iq,kb_sav(iq))+fldow(iq)*
     .        (s(iq,kb_sav(iq))+cp*(tdown(iq)-t(iq,kb_sav(iq))))
!      add contrib to cloud base layer (unit flux this one)
       delq(iq,kb_sav(iq))=delq(iq,kb_sav(iq))-qbase(iq)
       dels(iq,kb_sav(iq))=dels(iq,kb_sav(iq))-sbase(iq)
      enddo  ! iq loop
      if((ntest.ne.0.or.diag).and.mydiag)then
       print *,"delsa ",(dels(idjd,k),k=1,kt_sav(idjd))
       print *,"delqa ",(delq(idjd,k),k=1,kt_sav(idjd))
      endif

      if(alfsea.gt.0)then 
        fraca=1.    ! effectively used from 23/5/03
	fracb=0.    ! to avoid giving -ve qg for the layer
      else
        fraca=.5
	fracb=.5
      endif  ! (alfsea.gt.0)

!     subsidence and (possible) "shallow" detrainment
      do k=kuocb+1,kl-1
       do iq=1,ifull
        if(k.gt.kb_sav(iq).and.k.le.kt_sav(iq))then
          savg=fraca*s(iq,k)+fracb*s(iq,k-1)
          qavg=fraca*qq(iq,k)+fracb*qq(iq,k-1)
!         savg=s(iq,k)   ! 23/5/03
!         qavg=qq(iq,k)  ! 23/5/03 to avoid giving -ve qg for the layer
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
         savgb=fraca*s(iq,k)+fracb*s(iq,k-1)
         qavgb=fraca*qq(iq,k)+fracb*qq(iq,k-1)
!        savgb=s(iq,k)   ! 23/5/03
!        qavgb=qq(iq,k)  ! 23/5/03
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
     
!      N.B. pre convpsav & dsk; methdetr only used if detrainx>0.
       if(methdetr.eq.1)then                    ! moisten top layer only 
         do iq=1,ifull  
          k=kt_sav(iq)
          deltaq=qxcess(iq)                    ! moistening
          delq(iq,k)=delq(iq,k)+deltaq
          dels(iq,k)=dels(iq,k)-deltaq*hl
         enddo  ! iq loop
       elseif(methdetr.eq.-1)then               ! moisten bottom layer only 
         do iq=1,ifull  
          k=kb_sav(iq)+1
          deltaq=qxcess(iq)                    ! moistening
          delq(iq,k)=delq(iq,k)+deltaq
          dels(iq,k)=dels(iq,k)-deltaq*hl
         enddo  ! iq loop
       elseif(methdetr.eq.8)then        ! moisten all cloud layers, top most
         do iq=1,ifull  
          nlayers=kt_sav(iq)-kb_sav(iq)
          sum=.5*nlayers*(nlayers+1)
          do k=kb_sav(iq)+1,kt_sav(iq)
          deltaq=qxcess(iq)*(k-kb_sav(iq))/sum ! moistening
          delq(iq,k)=delq(iq,k)+deltaq
          dels(iq,k)=dels(iq,k)-deltaq*hl
          enddo  
         enddo   ! iq loop
       elseif(methdetr.eq.-8)then ! moisten all cloud layers, bottom most
        do iq=1,ifull  
         nlayers=kt_sav(iq)-kb_sav(iq)
         sum=.5*nlayers*(nlayers+1)
         do k=kb_sav(iq)+1,kt_sav(iq)
	   deltaq=qxcess(iq)*(kt_sav(iq)+1-k)/sum ! moistening
          delq(iq,k)=delq(iq,k)+deltaq
          dels(iq,k)=dels(iq,k)-deltaq*hl
         enddo  
        enddo   ! iq loop
       elseif(methdetr.eq.9)then                ! moisten all cloud layers
         do iq=1,ifull  
          nlayers=kt_sav(iq)-kb_sav(iq)
          do k=kb_sav(iq)+1,kt_sav(iq)
           deltaq=qxcess(iq)/nlayers           ! moistening
           delq(iq,k)=delq(iq,k)+deltaq
           dels(iq,k)=dels(iq,k)-deltaq*hl
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
        print *,'kb_sav,kt_sav,khalfd ',
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
!      if(kt_sav(iq).lt.kl)then 
       if(kb_sav(iq).lt.kl)then   ! from Oct 03, same answers, but more effic.
         if(nfluxdsk.eq.1)convpsav(iq)=dsk(kb_sav(iq))
         if(nfluxdsk.eq.2)convpsav(iq)=.5*dsk(kb_sav(iq))
         if(nfluxdsk.eq.-2)convpsav(iq)=dsk(kb_sav(iq))/(2*iterconv)
         if(nfluxq.eq.1)then  ! does little cf flux_dsk - will be removed
!          fluxq limiter: new_qq(kb)=new_qq(kt)
!          i.e. [qq+M*delq]_kb=[qq+M*delq]_kt
           fluxq(iq)=max(0.,(qq(iq,kb_sav(iq))-qq(iq,kt_sav(iq)))/ 
     .                      (delq(iq,kt_sav(iq)) - delq(iq,kb_sav(iq))))
           convpsav(iq)=min(convpsav(iq),fluxq(iq))
         endif  ! (nfluxq.eq.1)
       endif    ! (kb_sav(iq).lt.kl)
      enddo     ! iq loop
      if(rhcv.gt.0.)then  ! usually run with rhcv=0.
        do iq=1,ifull
         if(kb_sav(iq).lt.kl)then  
!          Base limiter: new_qq(kb)=rhcv*new_qs(kb), i.e.
!          [qq+M*delq]_kb=rhcv*[qs+M*dqsdt*dels/cp]_kb 
           fluxb(iq)=max(0.,(rhcv*qs(iq,kb_sav(iq))-qq(iq,kb_sav(iq)))/ 
     .            ( delq(iq,kb_sav(iq))
     .             -rhcv*dqsdt(iq,kb_sav(iq))*dels(iq,kb_sav(iq))/cp ) ) 
           convpsav(iq)=min(convpsav(iq),fluxb(iq))
         endif ! (kb_sav(iq).lt.kl)
        enddo  ! iq loop
      endif    ! (rhcv.gt.0.)
      
      do k=kuocb+1,kl-1
       do iq=1,ifull
!       want: new_hbas>=new_hs(k), i.e. in limiting case:
!       [h+alfsarr*M*dels+alfqarr*M*hl*delq]_base 
!                                          = [hs+M*dels+M*hlcp*dels*dqsdt]_k
        fluxt(iq,k)=(sbase(iq)+hl*qbase(iq)-hs(iq,k))/
     .              (dels(iq,k)*(1.+hlcp*dqsdt(iq,k))
     .              -alfsarr(iq)*dels(iq,kb_sav(iq))
     .              -alfqarr(iq)*hl*delq(iq,kb_sav(iq)) ) 

c        if(k.gt.kb_sav(iq).and.k.le.kt_sav(iq).and.
c     .     fluxt(iq,k).gt.0.)then
c           convpsav(iq)=min(convpsav(iq),fluxt(iq,k))
c        endif   ! (k.gt.kb_sav(iq).and.k.lt.kt_sav(iq).and.fluxt(iq,k).gt.0.)
!       may occasionally get -ve fluxt, for supersaturated layers
        if(k.gt.kb_sav(iq).and.k.le.kt_sav(iq))then
          convpsav(iq)=min(convpsav(iq),max(0.,fluxt(iq,k)))
        endif   ! (k.gt.kb_sav(iq).and.k.lt.kt_sav(iq))

       enddo    ! iq loop
      enddo     ! k loop      

      if(ntest.eq.2)then
        convmax=0.
        do iq=1,ifull
         if(convpsav(iq).gt.convmax.and.kb_sav(iq).eq.2)then
	    print *,'ktau,iq,convpsav,fluxt3,kb_sav2,kt_sav ',
     .            ktau,iq,convpsav(iq),fluxt(iq,3),kb_sav(iq),kt_sav(iq)
           idjd=iq
	    convmax=convpsav(iq)
         endif
        enddo	
      endif   !  (ntest.eq.2)
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
     .             flux_dsk,fluxq(idjd),fluxb(idjd),convpsav(idjd)
        write(6,"('fluxt',9f8.5/(5x,9f8.5))")
     .            (fluxt(iq,k),k=1,kt_sav(iq))
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
             write(6,"('bc  ',2i5,2i3,2x,2f6.3,f7.3,f6.3,' rh:',2f5.2)") 
     .       iq,nums,kb_sav(iq),kt_sav(iq),
     .       dsk(kb_sav(iq)),fluxb(iq),
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
             write(6,"('bcd ',2i5,2i3,2x,2f6.3,f7.3,f6.3,' rh:',2f5.2)") 
     .       iq,nums,kb_sav(iq),kt_sav(iq),
     .       dsk(kb_sav(iq)),fluxb(iq),
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
             write(6,"('bcde',2i5,2i3,2x,2f6.3,f7.3,f6.3,' rh:',2f5.2)") 
     .       iq,nums,kb_sav(iq),kt_sav(iq),
     .       dsk(kb_sav(iq)),fluxb(iq),
     .       fluxt(iq,kt_sav(iq)),convpsav(iq),
     .       qq(iq,kb_sav(iq))/qs(iq,kb_sav(iq)),
     .       qq(iq,kt_sav(iq))/qs(iq,kt_sav(iq))
           endif
         endif
        enddo  ! iq loop
      endif    ! (ntest.ne.0.or.diag)
      
      if(shaltime.gt.0.)then
        do iq=1,ifull
	  if(convshal(iq).lt.shaltime*3600.)convpsav(iq)=0.
	 enddo
      endif     ! (shaltime.gt.0.)	   

!      if(abs(ksc).gt.0)then  ! will do shallow clouds in vertmix
        if(sig_ct.lt.0.)then  ! use abs(sig_ct) as thickness of shallow clouds
          do iq=1,ifull
           if(sigmh(kb_sav(iq)+1)-sigmh(kt_sav(iq)+1).lt.-sig_ct)then  
             convpsav(iq)=0.         ! N.B. will get same result on later itns
	     if(ktsav(iq).eq.kl)then
               kbsav(iq)=kb_sav(iq) 
               ktsav(iq)=kt_sav(iq)  ! for possible use in vertmix
	     endif  ! (ktsav(iq).eq.kl)
           endif
          enddo  ! iq loop
	 else
          do iq=1,ifull
           if(sig(kt_sav(iq)).gt.sig_ct)then  ! typically sig_ct ~ .8
             convpsav(iq)=0.         ! N.B. will get same result on later itns
             if(ktsav(iq).eq.kl)then
               kbsav(iq)=kb_sav(iq) 
               ktsav(iq)=kt_sav(iq)  ! for possible use in vertmix
	     endif  ! (ktsav(iq).eq.kl)
           endif
          enddo  ! iq loop
        endif    !  (sig_ct.lt.0.).. else ..
!     endif
      
      do iq=1,ifull
       if(ktsav(iq).eq.kl.and.convpsav(iq).gt.0.)then
         kbsav(iq)=kb_sav(iq) 
         ktsav(iq)=kt_sav(iq)  
       endif  ! (ktsav(iq).eq.kl.and.convpsav(iq).gt.0.)
      enddo   ! iq loop

      if(itn.lt.iterconv)then 
        convpsav(:)=convfact*convpsav(:) ! typically convfact=1.02  
      endif                              ! (itn.lt.iterconv)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!     following 2 lines moved up (correctly) on 30/3/04
      qxcess(:)=detrain*rnrtcn(:)             ! e.g. .2* gives 20% detrainment
      rnrtcn(:)=rnrtcn(:)-qxcess(:)

!     update qq, tt and precip
      rnrtc(:)=rnrtc(:)+convpsav(:)*rnrtcn(:)*conrev(:) ! g/m**2/s
      do k=1,kl   
       do iq=1,ifull
        qq(iq,k)=qq(iq,k)+convpsav(iq)*delq(iq,k)
        tt(iq,k)=tt(iq,k)+convpsav(iq)*dels(iq,k)/cp
       enddo    ! iq loop
      enddo     ! k loop
      
!!!!!!!!!!!!!!!!!!!!! "deep" detrainment using detrain !!!!!v3!!!!!!    
!     N.B. convpsav has been updated here with convfact, but without factr term  
      if(methprec.eq.1)then                   ! moisten top layer only 
       do iq=1,ifull  
        deltaq=convpsav(iq)*qxcess(iq)/dsk(kt_sav(iq))             
        qliqw(iq,kt_sav(iq))=qliqw(iq,kt_sav(iq))+deltaq
       enddo  ! iq loop
      elseif(methprec.eq.8)then               ! moisten all cloud layers, top most
       do k=kuocb+1,kl-1 
        do iq=1,ifull  
         nlayers=kt_sav(iq)-kb_sav(iq)
         sum=.5*nlayers*(nlayers+1)
         if(k.gt.kb_sav(iq).and.k.le.kt_sav(iq))then
           deltaq=convpsav(iq)*qxcess(iq)*(k-kb_sav(iq))/(sum*dsk(k))    
           qliqw(iq,k)=qliqw(iq,k)+deltaq
         endif
        enddo  ! iq loop
       enddo   ! k loop
      elseif(methprec.eq.9)then               ! moisten all cloud layers
       do iq=1,ifull  
        nlayers=kt_sav(iq)-kb_sav(iq)
        do k=kb_sav(iq)+1,kt_sav(iq)
         deltaq=convpsav(iq)*qxcess(iq)/(nlayers*dsk(k))       
         qliqw(iq,k)=qliqw(iq,k)+deltaq
        enddo  
       enddo   ! iq loop
      elseif(methprec.eq.3)then               ! moisten middle third
       do k=kuocb+1,kl-1 
        do iq=1,ifull  
         k13=nint(kb_sav(iq)+.6 +(kt_sav(iq)-kb_sav(iq))/3.)  ! up by .1
         k23=nint(kb_sav(iq)+.4 +(kt_sav(iq)-kb_sav(iq))/1.5) ! down by .1
         nlayers=k23+1-k13
         if(k.ge.k13.and.k.le.k23)then
           deltaq=convpsav(iq)*qxcess(iq)/(nlayers*dsk(k))         
           qliqw(iq,k)=qliqw(iq,k)+deltaq
         endif
        enddo  ! iq loop
       enddo   ! k loop
      elseif(methprec.eq.-3)then          ! moisten mostly middle third
       do k=kuocb+1,kl-1 
        do iq=1,ifull  
	  sum=.25*(kt_sav(iq)-kb_sav(iq))**2
	  rkmid=.5*(kt_sav(iq)+kb_sav(iq)+1)
c        if(iq.eq.idjd)print *,'iq,kb_sav,kt_sav,sum,rkmid ',
c     .      	          iq,kb_sav(iq),kt_sav(iq),sum,rkmid
         if(k.gt.kb_sav(iq).and.k.le.kt_sav(iq))then
	    frac=kt_sav(iq)-rkmid+.5-abs(rkmid-k)
	    if(abs(rkmid-k).lt..1)frac=frac-.25    ! for central rkmid value
c          if(iq.eq.idjd)print *,'k,frac ',k,frac
           deltaq=convpsav(iq)*qxcess(iq)*frac/(sum*dsk(k))         
           qliqw(iq,k)=qliqw(iq,k)+deltaq
         endif
        enddo  ! iq loop
       enddo   ! k loop
      else                                    ! moisten upper layers
       do k=kuocb+1,kl-1 
        do iq=1,ifull  
         nlayers=max(1,nint((kt_sav(iq)-kb_sav(iq)-.1)/methprec))   ! round down
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
        write (6,"('qgc ',12f7.3/(8x,12f7.3))") 
     .             (1000.*qq(idjd,k),k=1,kl)
        write (6,"('ttc ',12f7.2/(8x,12f7.2))") 
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
        factuv=.1*nuvconv*factr ! full effect for nuvconv=10
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
          u(iq,k)=u(iq,k)+factuv*convpsav(iq)*delu(iq,k)
          v(iq,k)=v(iq,k)+factuv*convpsav(iq)*delv(iq,k)
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
        write (6,"('uuc ',12f6.1/(8x,12f6.1))") (u(idjd,k),k=1,kl)
        write (6,"('vvc ',12f6.1/(8x,12f6.1))") (v(idjd,k),k=1,kl)
      endif

      if(nmaxpr.eq.1.and.mydiag)then
        print *,'ktau,itn,kb_sav,kt_sav,convpsav ',
     .           ktau,itn,kb_sav(idjd),kt_sav(idjd),convpsav(idjd)
      endif

      enddo     ! itn=1,iterconv

      if(factr.lt.1.)then
        rnrtc(:)=factr*rnrtc(:)
        qq(1:ifull,:)=qg(1:ifull,:)+factr*(qq(1:ifull,:)-qg(1:ifull,:))      
        qliqw(1:ifull,:)=factr*qliqw(1:ifull,:)      
        tt(1:ifull,:)= t(1:ifull,:)+factr*(tt(1:ifull,:)- t(1:ifull,:))      
      endif  ! (factr.lt.1.)

!     update qq, tt for evap of qliqw (qliqw arose from moistening)
      if(ldr.ne.0)then
!       Leon's stuff here, e.g.
        do k=1,kl
          do iq=1,ifullw
!          qlg(iq,k)=qlg(iq,k)+qliqw(iq,k)
           if(tt(iq,k).lt.253.16)then   ! i.e. -20C
             qfg(iq,k)=qfg(iq,k)+qliqw(iq,k)
             tt(iq,k)=tt(iq,k)+3.35e5*qliqw(iq,k)/cp   ! fusion heating
	    else
             qlg(iq,k)=qlg(iq,k)+qliqw(iq,k)
	    endif
          enddo
         enddo
      else
        qq(1:ifull,:)=qq(1:ifull,:)+qliqw(1:ifull,:)         
        tt(1:ifull,:)=tt(1:ifull,:)-hl*qliqw(1:ifull,:)/cp   
      endif  ! (ldr.ne.0)
!__________________________end of convective calculations_____________________
     
      if((ntest.ne.0.or.diag).and.mydiag)then
        print *,"after convection"
        write (6,"('qge ',12f7.3/(8x,12f7.3))")(1000.*qq(idjd,k),k=1,kl)
        write (6,"('tte ',12f6.1/(8x,12f6.1))")(tt(idjd,k),k=1,kl)
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
      if(ldr.ne.0)go to 8

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
        write (6,"('qg ',12f7.3/(8x,12f7.3))") 
     .             (1000.*qq(idjd,k),k=1,kl)
        write (6,"('tt ',12f7.2/(8x,12f7.2))") 
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

      if(ntest.ne.0.or.diag.or.(ktau.eq.1.and.nmaxpr.eq.1))then
        if ( mydiag ) then
        print *,'at end of convjlm: rnrt,rnrtc ',rnrt(idjd),rnrtc(idjd)
        iq=idjd
        kdown=min(kl , kb_sav(iq)+nint(.6+.75*(kt_sav(iq)-kb_sav(iq))))
        print *,'kdown,tdown,qdown ',kdown,tdown(iq),qdown(iq)
        phi(iq,1)=bet(1)*tt(iq,1)
        do k=2,kl
         phi(iq,k)=phi(iq,k-1)+bet(k)*tt(iq,k)+betm(k)*tt(iq,k-1)
        enddo      ! k  loop
        do k=1,kl   
         es(iq,k)=establ(tt(iq,k))
         pk=ps(iq)*sig(k)
         qs(iq,k)=max(.622*es(iq,k)/(pk-es(iq,k)),1.5e-6)  
         s(iq,k)=cp*tt(iq,k)+phi(iq,k)  ! dry static energy
         hs(iq,k)=s(iq,k)+hl*qs(iq,k)   ! saturated moist static energy
         if(land(iq))then
           qbass(k)=alfqq_l(k)*qq(iq,k)+alfqq1_l(k)*qq(iq,1) 
         else
           qbass(k)=alfqq_s(k)*qq(iq,k)+alfqq1_s(k)*qq(iq,1) 
         endif
        enddo   ! k loop
        write (6,"('rh   ',12f7.2/(5x,12f7.2))") 
     .             (100.*qq(idjd,k)/qs(idjd,k),k=1,kl)
        write (6,"('qs   ',12f7.3/(5x,12f7.3))") 
     .             (1000.*qs(idjd,k),k=1,kl)
        write (6,"('qq   ',12f7.3/(5x,12f7.3))") 
     .             (1000.*qq(idjd,k),k=1,kl)
        write (6,"('tt   ',12f7.2/(5x,12f7.2))") 
     .             (tt(idjd,k),k=1,kl)
        write (6,"('s/cp ',12f7.2/(5x,12f7.2))") (s(idjd,k)/cp,k=1,kl)
        write (6,"('h/cp ',12f7.2/(5x,12f7.2))") 
     .             (s(idjd,k)/cp+hlcp*qq(idjd,k),k=1,kl)
        write (6,"('hb/cp',12f7.2/(5x,12f7.2))") 
     .             (s(idjd,k)/cp+hlcp*qbass(k),k=1,kl)
        write (6,"('hs/cp',12f7.2/(5x,12f7.2))") 
     .             (hs(idjd,k)/cp,k=1,kl)
        write (6,"('k  ',12i7/(3x,12i7))") (k,k=1,kl)
	 print *,'following are h,q,t changes during timestep'
        write (6,"('delh ',12f7.2/(5x,12f7.2))") 
     .             (s(idjd,k)/cp+hlcp*qq(idjd,k)-h0(k),k=1,kl)
        write (6,"('delq ',12f7.3/(5x,12f7.3))") 
     .             (1000.*(qq(idjd,k)-q0(k)),k=1,kl)
        write (6,"('delt  ',12f7.2/(5x,12f7.2))") 
     .             (tt(idjd,k)-t0(k),k=1,kl)
        flux_dsk=.5*dsk(kb_sav(idjd))  ! may depend on nfluxdsk
        write(6,"('flux_dsk,fluxq,fluxb,convpsav',5f8.5)")
     .             flux_dsk,fluxq(idjd),fluxb(idjd),convpsav(idjd)
        write(6,"('fluxt',9f8.5/(5x,9f8.5))")
     .            (fluxt(iq,k),k=1,kt_sav(iq))
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
      
      if(shaltime.gt.0.)then
        do iq=1,ifull
	  if(ktsav(iq).lt.kl)then
	    convshal(iq)=convshal(iq)+dt
	  else
	    convshal(iq)=0.
	  endif  ! (ktsav(iq).lt.kl)
	 enddo
      endif     ! (shaltime.gt.0.)	      
      return
      end
