      subroutine convjlm      ! jlm convective scheme - Version v3
!     the shallow convection options here just for iterconv=1
!     has +ve fldownn depending on delta sigma; -ve uses older abs(fldown)   
!     N.B. nevapcc option has been removed   
      use cc_mpi, only : mydiag
      use diag_m
      implicit none
      integer itn,iq,k,k13,k23,kb,kbsav_ls,kb_sav,kcl_top
     .       ,khalf,khalfd,khalfp,kt,kt_sav,ktmax
     .       ,nalfs,nfluxq,nfluxdsk,nlayers,nlayersd,nlayersp
     .       ,ntest,ntr,nums
      real Aev,Asb,alfcon,alfq,alfs,alfqarr,alfsarr
     .    ,conrev,convmax,delq,dels,deltaq,delq_av,delt_av
     .    ,den1,den2,den3,detr,dprec,dqrx,dqsdt
     .    ,es,entrainp,entrr,entrradd
     .    ,factdav,factr,facuv,fldow,fldownn
     .    ,fluxq,fluxr,fluxb,flux_dsk,fluxup,fluxt
     .    ,frac,fraca,fracb,gam,hbas,hbase,heatlev,hs
     .    ,pwater,pwater0,qavg,qavgb,qbas,qbase,qdown,qentrr,qliqw,qprec
     .    ,qq,qs,qsk,qxcess
     .    ,rkmid,rnrt,rnrtc,rnrtcn
     .    ,s,savg,savgb,sbas,sbase,sentrr,sum,sumqs,sumqt,sumss
     .    ,tm,totprec,tt,tdown,veldt
     
     .    ,rKa,Dva,cfls,cflscon,rhodz,qpf,pk,Apr,Bpr,Fr,rhoa,dz,Vr
     .    ,dtev,qr,qgdiff,Cev2,qr2,Cevx,alphal,blx,evapls,revq

      include 'newmpar.h'
      parameter (ntest=0)     ! 1 or 2 to turn on; -1 for ldr writes
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
!     parameter (mbase=0)      ! 1 cfrac test; 2 cfrac and omega test 
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
      include 'latlong.h'  ! rlatt
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
      real cfrac,duml,dumm
      common/nonlsav/cfrac(ifull,kl),duml(ifull,kl,2)
      common/work2/ktmax(ifull),conrev(ifull),fluxr(ifull)
     .  ,fluxb(ifull),dumm(ifull),kbsav_ls(ifull)
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
!!    equivalence (phi,s),(delu,revc,delq),(delv,dels)
      equivalence (delu,revc,delq),(delv,dels)
!     data rhcv/.75/                  now in kuocom
      integer kdown(ifull),kscbase,ksctop
      real prcpv(kl),theeb(ifull),thee(ifull,kl)
      equivalence (phi,thee)
      save kscbase,ksctop,prcpv
      real entr(ifull),qentr(ifull),sentr(ifull),factr(ifull)
!     real convshal(ifull)
!     save convshal
!     data convshal/ifull*0./

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
      if(convtime>10.)then                          !  ds-dependent     
!       following uses convtime_eff(hrs)=convtime(~60)/dx(in km)  
        factr(:)=dt/max(dt,convtime*3.6e6*em(1:ifull)/ds) ! to re-scale convpsav
      elseif(convtime<0.)then         
!       following behaves like 2*dt/(-convtime*3600.) for small dt     
        factr(:)=1.-(1.+dt/(convtime*3600.))**2     ! to re-scale convpsav 
      else                                          ! original option
        factr(:)=dt/max(dt,convtime*3600.)          ! to re-scale convpsav
!     testing linear from -15 to 0 to 15:  1 to 4 to 1
c       do iq=1,ifull
c	  convmax=convtime*max(1.,1.+(15.-abs(rlatt(iq)*180./3.14))/5.)
c        factr(iq)=dt/max(dt,convmax*3600.)          ! to re-scale convpsav
c	 enddo
      endif  ! (convtime<0.)
c     following modification was found to be not much use
c     do iq=1,ifull
c       k=kb_sav(iq)
c       if(qq(iq,k)+qlg(iq,k)+qfg(iq,k)>qs(iq,k))factr(iq)=1.
c     enddo

!__________________________beginning of convective calculations_____________________
      do itn=1,iterconv
!     calculate geopotential height 
      kb_sav(:)=kl     ! preset value for no convection
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
       print *,'itn,iterconv,nuvconv ',
     .          itn,iterconv,nuvconv
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
      if((ntest>0.or.diag.or.nmaxpr.eq.1).and.mydiag) then
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
      hs(:,:)=max(hs(:,:),s(:,:)+hl*qq(:,:))   ! added 21/7/04
      do k=kl/2,kuocb+1,-1   ! downwards to find cloud base
       kdown(:)=1
       if(rhcv>0.)then
         do iq=1,ifull
          if(qq(iq,k-1)<=rhcv*qs(iq,k-1))kdown(iq)=0
         enddo
       endif
       if(ktau>1.and.rhcv<0.)then
         do iq=1,ifull
          if(phi(iq,k-1)/9.806>pblh(iq))kdown(iq)=0
         enddo
       endif
       if(mbase==1)then
         do iq=1,ifull
          if(cfrac(iq,k)==0.)kdown(iq)=0
         enddo
       endif
       if(mbase==2)then
         do iq=1,ifull
          if(dpsldt(iq,k)>dsig2.or.cfrac(iq,k)==0.)kdown(iq)=0
         enddo
       endif
       do iq=1,ifull
!       find tentative cloud base, and bottom of below-cloud layer
        if(kdown(iq)==1)then   ! kdown just 0 or 1 here
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
          if(sbas+hl*qbas>max(hs(iq,k),sbase(iq)+hl*qbase(iq)))then
            if(qbas>max(qs(iq,k),qq(iq,k)))then  ! newer qbas test 
c           if(qbas>max(qs(iq,k),qq(iq,k)).and.cfrac(iq,k)>0.)then  ! newer qbas test 
c            if(dpsldt(iq,k)<dsig2)then   ! omega test july '04
              kb_sav(iq)=k-1
              qbase(iq)=qbas
              sbase(iq)=sbas
              kt_sav(iq)=k
              alfqarr(iq)=alfq    ! for use in flux calc
              alfsarr(iq)=alfs    ! for use in flux calc
c            endif   ! (dpsldt(iq,k)<dsig2)
            endif ! (qbas>max(qs(iq,k),qq(iq,k)))
          endif   ! (sbas+hl*qbas>hs(iq,k))
        endif   ! (kdown(iq)==1)
       enddo    ! iq loop
       if(ntest>0.and.mydiag) then
	  iq=idjd
         if(land(iq))then
           qbas=alfqq_l(k-1)*qq(iq,k-1)+alfqq1_l(k-1)*qq(iq,1) 
           alfq=alfqq_l(k-1)   
         else
           qbas=alfqq_s(k-1)*qq(iq,k-1)+alfqq1_s(k-1)*qq(iq,1) 
           alfq=alfqq_s(k-1)   
         endif  ! (land(iq)) .. else ..
	  qbas=min(qbas,max(qs(iq,k-1),qq(iq,k-1)))   
         alfs=alfss(k-1)  ! simple one
         sbas=alfs*s(iq,k-1)+(1.-alfs)*s(iq,1) 
	  print *,'k,kdown,kb_sav,kt_sav ',
     &            k,kdown(iq),kb_sav(iq),kt_sav(iq)
	  print *,'phi(,k-1)/g,pblh ',phi(iq,k-1)/9.806,pblh(iq)
	  print *,'qbas,qs(.,k),qq(.,k) ',qbas,qs(iq,k),qq(iq,k)
         print *,'k,(sbas+hl*qbas)/cp,hs(iq,k)/cp ',
     .            k,(sbas+hl*qbas)/cp,hs(iq,k)/cp
       endif    ! (ntest>0.and.mydiag) 
      enddo     ! k loop

!     subtract contrib to cloud base layer (unit flux this one)
      do iq=1,ifull
       delq(iq,kb_sav(iq))=delq(iq,kb_sav(iq))-qbase(iq)
       dels(iq,kb_sav(iq))=dels(iq,kb_sav(iq))-sbase(iq)
      enddo
      entrainp=1.+abs(entrain)
      entr(:)=0.
      qentr(:)=0.
      sentr(:)=0.
      if(entrain>=0.)then
        do k=kuocb+2,kl   ! upwards to find cloud top
         do iq=1,ifull
          if(kt_sav(iq).eq.k-1.and.k.lt.ktmax(iq))then ! ktmax allows for itn
            hbas=sbase(iq)+hl*qbase(iq)
            entrr=-entrain*dsig(k-1)
            qentrr=qentr(iq)+entrr*qq(iq,k-1)
            sentrr=sentr(iq)+entrr*s(iq,k-1)
            entrradd=entr(iq)+entrr
            hbase=(hbas+sentrr+hl*qentrr)/(1.+entrradd)
            if(hbase.gt.hs(iq,k))then
              kt_sav(iq)=k
              entr(iq)=entrradd
              qentr(iq)=qentrr
              sentr(iq)=sentrr
	     endif
          endif   ! (kt_sav(iq).eq.k-1.and.k.lt.ktmax(iq))
         enddo    ! iq loop
        enddo     ! k loop
!       update top-emerging qbase & sbase to include entrainment
        qbase(:)=(qbase(:)+qentr(:))/(1.+entr(:))
        sbase(:)=(sbase(:)+sentr(:))/(1.+entr(:))
        entr(:)=entrain
      endif  ! (entrain>=0.)
      if(entrain<0.)then
        fluxr(:)=1.  ! just for qbase & sbase denom (allows for 2-layer clouds)
        do k=kuocb+2,kl   ! upwards to find cloud top
         do iq=1,ifull
c	   iq=idjd
          if(kt_sav(iq).eq.k-1.and.k.lt.ktmax(iq))then ! ktmax allows for itn
            hbas=sbase(iq)+hl*qbase(iq)
            entrr=abs(entrain)/(sigmh(kb_sav(iq)+1)-sigmh(k))
            qentrr=qentr(iq)-dsig(k-1)*qq(iq,k-1)
            sentrr=sentr(iq)-dsig(k-1)*s(iq,k-1)
            hbase=(hbas+entrr*(sentrr+hl*qentrr))/entrainp
c	     print *,'k,hbas,entrr,qentrr,sentrr ',
c     .               k,hbas,entrr,qentrr,sentrr 	     
c	     print *,'dsig(k-1),hbase,hs ',dsig(k-1),hbase,hs(iq,k)
            if(hbase.gt.hs(iq,k))then
              kt_sav(iq)=k
              entr(iq)=entrr
              qentr(iq)=qentrr
              sentr(iq)=sentrr
		fluxr(iq)=entrainp
	     endif
          endif   ! (kt_sav(iq).eq.k-1.and.k.lt.ktmax(iq))
         enddo    ! iq loop
        enddo     ! k loop
!       update top-emerging qbase & sbase to include entrainment
        qbase(:)=(qbase(:)+entr(:)*qentr(:))/fluxr(:)
        sbase(:)=(sbase(:)+entr(:)*sentr(:))/fluxr(:)
      endif  ! (entrain<0.)
      
c      if(itn.gt.1)then  ! added April 04
c        do iq=1,ifull
c	  if(rnrtc(iq).gt.0..and.kb_sav(iq).eq.kbsav(iq).
c     .                      and.kt_sav(iq).eq.ktsav(iq))then
c	    kt_sav(iq)=max(kbsav(iq)+1,ktsav(iq)-1)
c	  endif   !  (rnrtc(iq).gt.0.....)
c	 enddo
c      endif   ! (itn.gt.1)

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

      if((ntest>0.or.diag).and.mydiag)then
       iq=idjd
       print *,"sbase,qbase ",sbase(iq),qbase(iq)
       entrradd=1.+entr(iq)*(sigmh(kb_sav(iq)+1)-sigmh(kt_sav(iq)))
       print *,"ktmax,ktsav,entr,entrradd ",
     .          ktmax(iq),ktsav(iq),entr(iq),entrradd
       kdown(iq)=min(kl,kb_sav(iq)+nint(.6+.75*(kt_sav(iq)-kb_sav(iq))))
       print *,"ktau,itn,kb_sav,kt_sav,kdown ",
     .    ktau,itn,kb_sav(iq),kt_sav(iq),kdown(iq)
       print *,'alfqarr,alfsarr ',alfqarr(iq),alfsarr(iq)
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     calculate dels and delq for all layers for unit base mass flux
!     designed to be OK if q(k)>qs(k)
      do iq=1,ifull
!      dels and delq for top layer
!      following line allows for qq>qs, to avoid drying layer
       qsk=max(qs(iq,kt_sav(iq)),qq(iq,kt_sav(iq)))  
!      qprec=max(0.,qbase(iq)-qsk)             
!      dels(iq,kt_sav(iq))=dels(iq,kt_sav(iq))+sbase(iq)   ! s flux
!      delq(iq,kt_sav(iq))=delq(iq,kt_sav(iq))+qsk  
       entrradd=1.+entr(iq)*(sigmh(kb_sav(iq)+1)-sigmh(kt_sav(iq)))
       qprec=entrradd*max(0.,qbase(iq)-qsk)             
       dels(iq,kt_sav(iq))=dels(iq,kt_sav(iq))+entrradd*sbase(iq) ! s flux
       delq(iq,kt_sav(iq))=delq(iq,kt_sav(iq))+entrradd*qsk  

       dels(iq,kt_sav(iq))=dels(iq,kt_sav(iq))+hl*qprec    ! precip. heating
       kdown(iq)=min(kl,kb_sav(iq)+nint(.6+.75*(kt_sav(iq)-kb_sav(iq))))
       qsk=min(qs(iq,kdown(iq)),qq(iq,kdown(iq))) ! N.B. min for downdraft subs
       tdown(iq)=t(iq,kb_sav(iq))+(s(iq,kdown(iq))-s(iq,kb_sav(iq))
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
         print *,'kdown,tdown,qdown ',kdown(iq),tdown(iq),qdown(iq)
         print *,'qsk,qprec,dels0 ',qsk,qprec,hl*qprec 
         print *,'fldow,dprec,rnrtcn ',fldow(iq),dprec,rnrtcn(iq)
       endif
!      add in downdraft contributions
       dels(iq,kdown(iq))=dels(iq,kdown(iq))-fldow(iq)*s(iq,kdown(iq))
       delq(iq,kdown(iq))=delq(iq,kdown(iq))-fldow(iq)*qsk
!      calculate emergent downdraft properties
       delq(iq,kb_sav(iq))=delq(iq,kb_sav(iq))+fldow(iq)*qdown(iq)
       dels(iq,kb_sav(iq))=dels(iq,kb_sav(iq))+fldow(iq)*
     .        (s(iq,kb_sav(iq))+cp*(tdown(iq)-t(iq,kb_sav(iq))))
      enddo  ! iq loop
      if((ntest>0.or.diag).and.mydiag)then
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
        if(k>kb_sav(iq).and.k<=kt_sav(iq))then
          savg=(fraca*s(iq,k)+fracb*s(iq,k-1))
     .         *(1.+(sigmh(kb_sav(iq)+1)-sigmh(k))*entr(iq))
          qavg=(fraca*qq(iq,k)+fracb*qq(iq,k-1))
     .         *(1.+(sigmh(kb_sav(iq)+1)-sigmh(k))*entr(iq))
!         savg=s(iq,k)   ! 23/5/03
!         qavg=qq(iq,k)  ! 23/5/03 to avoid giving -ve qg for the layer
          dels(iq,k)=dels(iq,k)-savg             ! subsidence 
          dels(iq,k-1)=dels(iq,k-1)+savg         ! subsidence into lower layer
          delq(iq,k)=delq(iq,k)-qavg             ! subsidence 
          delq(iq,k-1)=delq(iq,k-1)+qavg         ! subsidence into lower layer
          if(k<kt_sav(iq))then
            dels(iq,k)=dels(iq,k)+dsig(k)*entr(iq)*s(iq,k)  ! entr into updraft
            delq(iq,k)=delq(iq,k)+dsig(k)*entr(iq)*qq(iq,k) ! entr into updraft
	   endif
        endif  ! (k>kb_sav(iq).and.k<=kt_sav(iq))
       enddo   ! iq loop
      enddo    ! k loop
      if((ntest>0.or.diag).and.mydiag)then
       print *,"delsb ",(dels(idjd,k),k=1,kt_sav(idjd))
       print *,"delqb ",(delq(idjd,k),k=1,kt_sav(idjd))
      endif

!     modify calculated subsidence for downdraft effects
      do k=kuocb+1,kl-1
       do iq=1,ifull
        if(k.le.kdown(iq).and.k.gt.kb_sav(iq))then
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
      if((ntest>0.or.diag).and.mydiag)then
       print *,"delsc ",(dels(idjd,k),k=1,kt_sav(idjd))
       print *,"delqc ",(delq(idjd,k),k=1,kt_sav(idjd))
      endif

!!!!!!!!! "shallow" detrainment using detrainx, dsig2, dsig4 !!!!!v3!!!!!! 
      if(detrainx.gt.0.)then   
       do iq=1,ifull  
!       detr=min(detrainx,max( 0.,
!    .  (dsig4-sig(kb_sav(iq))+sig(kt_sav(iq)))*detrainx/(dsig4-dsig2)))
        detr=detrainx
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

      if((ntest>0.or.diag).and.mydiag)then
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
        elseif(k.eq.kb_sav(iq))then  ! s and qg only from kb layer
          delq(iq,k)=delq(iq,k)/dsk(k)
          dels(iq,k)=dels(iq,k)/dsk(k)
        endif    ! (k.gt.kb_sav(iq)) .. else ..
       enddo     ! iq loop
      enddo      ! k loop

      if((ntest>0.or.diag).and.mydiag)then
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
!       pre 21/7/04 occasionally got -ve fluxt, for supersaturated layers
        if(k.gt.kb_sav(iq).and.k.le.kt_sav(iq))then
c	   if(          (dels(iq,k)*(1.+hlcp*dqsdt(iq,k))  ! test for Steve
c     .                -alfsarr(iq)*dels(iq,kb_sav(iq))
c     .                -alfqarr(iq)*hl*delq(iq,kb_sav(iq)) ).eq.0.)then
c	     den1=dels(iq,k)*(1.+hlcp*dqsdt(iq,k))
c	     den2=alfsarr(iq)*dels(iq,kb_sav(iq))
c	     den3=alfqarr(iq)*hl*delq(iq,kb_sav(iq))
c           print *,'k,iq,kb_sav(iq),kt_sav(iq),den1,den2,den3 ',
c     .              k,iq,kb_sav(iq),kt_sav(iq),den1,den2,den3
c	     stop
c         endif 
          fluxt(iq,k)=(sbase(iq)+hl*qbase(iq)-hs(iq,k))/
     .                (dels(iq,k)*(1.+hlcp*dqsdt(iq,k))
     .                -alfsarr(iq)*dels(iq,kb_sav(iq))
     .                -alfqarr(iq)*hl*delq(iq,kb_sav(iq)) ) 
          convpsav(iq)=min(convpsav(iq),max(0.,fluxt(iq,k)))
        endif   ! (k.gt.kb_sav(iq).and.k.lt.kt_sav(iq))
       enddo    ! iq loop
       if((ntest>0.or.diag).and.mydiag)then
        iq=idjd
        if(k.gt.kb_sav(iq).and.k.le.kt_sav(iq))then
	   den1=dels(iq,k)*(1.+hlcp*dqsdt(iq,k))
	   den2=alfsarr(iq)*dels(iq,kb_sav(iq))
	   den3=alfqarr(iq)*hl*delq(iq,kb_sav(iq))
          fluxt(iq,k)=(sbase(iq)+hl*qbase(iq)-hs(iq,k))/
     .                (den1-den2-den3) 
          print *,'k,den1,den2,den3,fluxt ',k,den1,den2,den3,fluxt(iq,k)
        endif   ! (k.gt.kb_sav(iq).and.k.lt.kt_sav(iq))
       endif   ! ((ntest>0.or.diag).and.mydiag)
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
      if((ntest>0.or.diag).and.mydiag)then
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
      endif    ! (ntest>0.or.diag)
      
!     if(shaltime.gt.0.)then
!       do iq=1,ifull
!        if(convshal(iq).lt.shaltime*3600.)convpsav(iq)=0.
!       enddo
!     endif     ! (shaltime.gt.0.)	   

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
      
      do iq=1,ifull
       if(ktsav(iq).eq.kl.and.convpsav(iq).gt.0.)then
         kbsav(iq)=kb_sav(iq) 
         ktsav(iq)=kt_sav(iq)  
       endif  ! (ktsav(iq).eq.kl.and.convpsav(iq).gt.0.)
      enddo   ! iq loop

      if(itn.lt.iterconv)then 
        convpsav(:)=convfact*convpsav(:) ! typically convfact=1.02  
        ktmax(:)=kt_sav(:)    ! added July 04
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
      if(methprec.eq.1)then               ! moisten top layer only 
       do iq=1,ifull  
        deltaq=convpsav(iq)*qxcess(iq)/dsk(kt_sav(iq))             
        qliqw(iq,kt_sav(iq))=qliqw(iq,kt_sav(iq))+deltaq
       enddo  ! iq loop
      elseif(methprec.eq.8)then           ! moisten all cloud layers, top most
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
      elseif(methprec.eq.9)then           ! moisten all cloud layers
       do iq=1,ifull  
        nlayers=kt_sav(iq)-kb_sav(iq)
        do k=kb_sav(iq)+1,kt_sav(iq)
         deltaq=convpsav(iq)*qxcess(iq)/(nlayers*dsk(k))       
         qliqw(iq,k)=qliqw(iq,k)+deltaq
        enddo  
       enddo   ! iq loop
      elseif(methprec.eq.3)then           ! moisten middle third
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
      else                                ! moisten upper layers
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

      if((ntest>0.or.diag).and.mydiag)then
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
        heatlev=0.
        do k=1,kl
         delq_av=delq_av+dsk(k)*convpsav(iq)*delq(iq,k)
         delt_av=delt_av+dsk(k)*convpsav(iq)*dels(iq,k)/cp
         heatlev=heatlev+sig(k)*dsk(k)*convpsav(iq)*dels(iq,k)/cp
        enddo
        print *,'delq_av,delt_exp,rnd_exp ',
     .         delq_av,-delq_av*hl/cp,-delq_av*conrev(iq)
        if(delt_av.ne.0.)print *,'ktau,itn,kbsav,ktsav,delt_av,heatlev',
     .        ktau,itn,kb_sav(idjd),kt_sav(idjd),delt_av,heatlev/delt_av
      endif   ! (ntest>0.or.diag)

      if(nuvconv>0.and.nuvconv<100)then !  momentum calculations, original style
        facuv=.1*nuvconv ! full effect for nuvconv=10
        delu(:,:)=0.
        delv(:,:)=0.
	 do k=kl,1,-1
	  do iq=1,ifull
          if(k==kt_sav(iq))then                 ! delu and delv for top layer
            delu(iq,k)=u(iq,kb_sav(iq))-u(iq,k)
            delv(iq,k)=v(iq,kb_sav(iq))-v(iq,k)
          elseif(k>=kb_sav(iq).and.k<kt_sav(iq))then
            delu(iq,k)=u(iq,k+1)-u(iq,k)       ! subs.
            delv(iq,k)=v(iq,k+1)-v(iq,k)       ! subs.
          endif  ! (k==kb_sav(iq))  .. else ..      
         enddo   ! iq loop
        enddo    ! k loop
!       update u & v using actual delu and delv (i.e. divided by dsk)
        do k=1,kl-1   
	  do iq=1,ifull
          u(iq,k)=u(iq,k)+facuv*factr(iq)*convpsav(iq)*delu(iq,k)/dsk(k)
          v(iq,k)=v(iq,k)+facuv*factr(iq)*convpsav(iq)*delv(iq,k)/dsk(k)
         enddo  ! iq loop
        enddo   ! k loop
      endif     ! (nuvconv>0.and.nuvconv<100)

      if(nuvconv>100.and.nuvconv<200)then   !  momentum, with downdrafts
        facuv=.1*(nuvconv-100) ! full effect for nuvconv=110
        delu(:,:)=0.
        delv(:,:)=0.
	 do k=kl,1,-1
	  do iq=1,ifull
          if(k==kt_sav(iq))then                 ! delu and delv for top layer
            delu(iq,k)=u(iq,kb_sav(iq))-u(iq,k)
            delv(iq,k)=v(iq,kb_sav(iq))-v(iq,k)
          elseif(k>=kdown(iq).and.k<kt_sav(iq))then
            delu(iq,k)=u(iq,k+1)-u(iq,k)       ! subs.
            delv(iq,k)=v(iq,k+1)-v(iq,k)       ! subs.
          elseif(k>kb_sav(iq).and.k<kdown(iq))then
            delu(iq,k)=(1.-fldow(iq))*(u(iq,k+1)-u(iq,k))   ! subs.
            delv(iq,k)=(1.-fldow(iq))*(v(iq,k+1)-v(iq,k))   ! subs.
          elseif(k==kb_sav(iq))then
            delu(iq,k)=(1.-fldow(iq))*u(iq,k+1)+
     .                  fldow(iq)*u(iq,kdown(iq))-u(iq,k)   ! subs.
            delv(iq,k)=(1.-fldow(iq))*v(iq,k+1)+
     .                  fldow(iq)*v(iq,kdown(iq))-v(iq,k)   ! subs.
          endif  ! (k==kb_sav(iq))  .. else ..    
         enddo   ! iq loop
        enddo    ! k loop
!       update u & v using actual delu and delv (i.e. divided by dsk)
        do k=1,kl-1   
	  do iq=1,ifull
          u(iq,k)=u(iq,k)+facuv*factr(iq)*convpsav(iq)*delu(iq,k)/dsk(k)
          v(iq,k)=v(iq,k)+facuv*factr(iq)*convpsav(iq)*delv(iq,k)/dsk(k)
         enddo  ! iq loop
        enddo   ! k loop
      endif     ! (nuvconv>100.and.nuvconv<200)

      if(nuvconv>200.and.nuvconv<300)then !  momentum, general mixing
!       this one assumes horiz. pressure gradients will be generated on parcel      
        facuv=.1*(nuvconv-200) ! full effect for nuvconv=210
        delu(:,:)=0.
        delv(:,:)=0.
	 do k=1,kl-1
	  do iq=1,ifull
          if(k==kt_sav(iq))then                 ! delu and delv for top layer
            delu(iq,k)=u(iq,k-1)-u(iq,k)
            delv(iq,k)=v(iq,k-1)-v(iq,k)
          elseif(k>kb_sav(iq).and.k<kt_sav(iq))then
            delu(iq,k)=u(iq,k+1)+u(iq,k-1)-2.*u(iq,k)       
            delv(iq,k)=v(iq,k+1)+v(iq,k-1)-2.*v(iq,k)      
          elseif(k==kb_sav(iq))then
            delu(iq,k)=u(iq,k+1)-u(iq,k)       
            delv(iq,k)=v(iq,k+1)-v(iq,k)      
          endif  ! (k==kb_sav(iq))  .. else ..      
         enddo   ! iq loop
        enddo    ! k loop
!       update u & v using actual delu and delv (i.e. divided by dsk)
        do k=1,kl-1   
	  do iq=1,ifull
          u(iq,k)=u(iq,k)+facuv*factr(iq)*convpsav(iq)*delu(iq,k)/dsk(k)
          v(iq,k)=v(iq,k)+facuv*factr(iq)*convpsav(iq)*delv(iq,k)/dsk(k)
         enddo  ! iq loop
        enddo   ! k loop
      endif     ! (nuvconv>200.and.nuvconv<300)

      if(nuvconv>300.and.nuvconv<400)then !  mixture of 10 & 210
!       this one assumes horiz. pressure gradients will be generated on parcel      
        facuv=.1*(nuvconv-300) ! full effect for nuvconv=310
        delu(:,:)=0.
        delv(:,:)=0.
	 do k=kl,1,-1
	  do iq=1,ifull
          if(k==kt_sav(iq))then                 ! delu and delv for top layer
            delu(iq,k)=.5*(u(iq,kb_sav(iq))+u(iq,k-1))-u(iq,k)
            delv(iq,k)=.5*(v(iq,kb_sav(iq))+v(iq,k-1))-v(iq,k)
          elseif(k>kb_sav(iq).and.k<kt_sav(iq))then
            delu(iq,k)=.5*(u(iq,k+1)-u(iq,k)       ! subs.
     .                    +u(iq,k+1)+u(iq,k-1)-2.*u(iq,k) )      
            delv(iq,k)=.5*(v(iq,k+1)-v(iq,k)       ! subs.
     .                    +v(iq,k+1)+v(iq,k-1)-2.*v(iq,k) )      
          elseif(k==kb_sav(iq))then
            delu(iq,k)=u(iq,k+1)-u(iq,k)       
            delv(iq,k)=v(iq,k+1)-v(iq,k)      
          endif  ! (k==kb_sav(iq))  .. else ..      
         enddo   ! iq loop
        enddo    ! k loop
!       update u & v using actual delu and delv (i.e. divided by dsk)
        do k=1,kl-1   
	  do iq=1,ifull
          u(iq,k)=u(iq,k)+facuv*factr(iq)*convpsav(iq)*delu(iq,k)/dsk(k)
          v(iq,k)=v(iq,k)+facuv*factr(iq)*convpsav(iq)*delv(iq,k)/dsk(k)
         enddo  ! iq loop
        enddo   ! k loop
      endif     ! (nuvconv>200.and.nuvconv<300)

!     if(nuvconv>400.and.nuvconv<500)then ! pseudo implicit
!       this one assumes horiz. pressure gradients will be generated on parcel      
!       updated u & v using actual delu and delv (i.e. divided by dsk)
!     endif     ! (nuvconv>400.and.nuvconv<500)

      if(nuvconv>500.and.nuvconv<600)then   !  with extended downdrafts
        facuv=.1*(nuvconv-500) ! full effect for nuvconv=510
        delu(:,:)=0.
        delv(:,:)=0.
	 do k=1,kl-1
	  do iq=1,ifull
	  if(kt_sav(iq)<kl)then
          if(k==kt_sav(iq))then                 ! delu and delv for top layer
            delu(iq,k)=u(iq,kb_sav(iq))-u(iq,k)
            delv(iq,k)=v(iq,kb_sav(iq))-v(iq,k)
          elseif(k>=kdown(iq).and.k<kt_sav(iq))then
            delu(iq,k)=u(iq,k+1)-u(iq,k)       ! subs.
            delv(iq,k)=v(iq,k+1)-v(iq,k)       ! subs.
          elseif(k>kb_sav(iq).and.k<kdown(iq))then
            delu(iq,k)=(1.-fldow(iq))*(u(iq,k+1)-u(iq,k))   ! subs.
            delv(iq,k)=(1.-fldow(iq))*(v(iq,k+1)-v(iq,k))   ! subs.
          elseif(k==kb_sav(iq))then
            delu(iq,k)=(1.-fldow(iq))*u(iq,k+1)+            ! subs.
     .                      fldow(iq)*u(iq,k-1) -u(iq,k)    ! ups.
            delv(iq,k)=(1.-fldow(iq))*v(iq,k+1)+            ! subs.
     .                      fldow(iq)*v(iq,k-1) -v(iq,k)    ! ups.
          elseif(k>1.and.k<kb_sav(iq))then
            delu(iq,k)=fldow(iq)*(u(iq,k-1)-u(iq,k))        ! ups.
            delv(iq,k)=fldow(iq)*(v(iq,k-1)-v(iq,k))        ! ups.
          elseif(k==1)then
            delu(iq,k)=fldow(iq)*(u(iq,kdown(iq))-u(iq,k))   
            delv(iq,k)=fldow(iq)*(v(iq,kdown(iq))-v(iq,k))
          endif  ! (k==kb_sav(iq))  .. else ..    
	  endif   ! (kt_sav(iq)<kl)
         enddo   ! iq loop
        enddo    ! k loop
!       update u & v using actual delu and delv (i.e. divided by dsk)
        do k=1,kl-1   
	  do iq=1,ifull
          u(iq,k)=u(iq,k)+facuv*factr(iq)*convpsav(iq)*delu(iq,k)/dsk(k)
          v(iq,k)=v(iq,k)+facuv*factr(iq)*convpsav(iq)*delv(iq,k)/dsk(k)
         enddo  ! iq loop
        enddo   ! k loop
      endif     ! (nuvconv>500.and.nuvconv<600)

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
            veldt=factr(iq)*convpsav(iq)*(1.-fldow(iq))  ! simple treatment
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
      
      if((ntest>0.or.diag).and.mydiag)then
        write (6,"('uuc ',12f6.1/(8x,12f6.1))") (u(idjd,k),k=1,kl)
        write (6,"('vvc ',12f6.1/(8x,12f6.1))") (v(idjd,k),k=1,kl)
      endif

      if(nmaxpr.eq.1.and.mydiag)then
        print *,'ktau,itn,kb_sav,kt_sav,convpsav ',
     .           ktau,itn,kb_sav(idjd),kt_sav(idjd),convpsav(idjd)
      endif

      enddo     ! itn=1,iterconv

      rnrtc(:)=factr(:)*rnrtc(:)
      do k=1,kl
      qq(1:ifull,k)=qg(1:ifull,k)+factr(:)*(qq(1:ifull,k)-qg(1:ifull,k))      
      qliqw(1:ifull,k)=factr(:)*qliqw(1:ifull,k)      
      tt(1:ifull,k)= t(1:ifull,k)+factr(:)*(tt(1:ifull,k)- t(1:ifull,k))      
      enddo

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
     
      if((ntest>0.or.diag).and.mydiag)then
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

      if((ntest>0.or.diag).and.mydiag)then
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

      if(ntest>0.or.diag.or.(ktau.eq.1.and.nmaxpr.eq.1))then
        if ( mydiag ) then
        print *,'at end of convjlm: rnrt,rnrtc ',rnrt(idjd),rnrtc(idjd)
        iq=idjd
        print *,'kdown,tdown,qdown ',kdown(iq),tdown(iq),qdown(iq)
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
      
!     if(shaltime.gt.0.)then
!       do iq=1,ifull
!        if(ktsav(iq).lt.kl)then
!          convshal(iq)=convshal(iq)+dt
!        else
!          convshal(iq)=0.
!        endif  ! (ktsav(iq).lt.kl)
!       enddo
!     endif     ! (shaltime.gt.0.)	      

      if(ksc>70.and.ksc<80)then
       if(ktau.eq.1)then
!       set ksctop for shallow convection
        ksctop=1    ! ksctop will be first level below sigksct
        do while(sig(ksctop+1).gt.sigksct)  !  e.g. sigksct=.75
         ksctop=ksctop+1
        enddo
        kscbase=1  ! kscbase will be first level above sigkcsb
        do while(sig(kscbase).gt.sigkscb.and.sigkscb.gt.0.) ! e.g. sigkscb=.99
         kscbase=kscbase+1
        enddo
        if (mydiag) then
        print *,'For shallow convection, set in convjlm:'
        print *,'ksc,kscbase,ksctop,kscsea ',
     .           ksc,kscbase,ksctop,kscsea
	 write (6,"(' sigkscb,sigksct,tied_con,tied_over,tied_rh:',
     .       5f8.3)")sigkscb,sigksct,tied_con,tied_over,tied_rh
        end if
        do k=1,kl
         prcpv(k)=sig(k)**(-roncp)
        enddo  ! k loop
       endif   ! (ktau.eq.1)
       if(ksc.eq.79)then    ! resembles ksc=99
        kt_sav(:)=ktsav(:)  ! as a working array
        do iq=1,ifull
         theeb(iq)=prcpv(kscbase)*t(iq,kscbase)*
     .                  (t(iq,kscbase) + .5*hlcp*qs(iq,kscbase))
     .                 /(t(iq,kscbase) - .5*hlcp*qs(iq,kscbase))
        enddo    ! iq loop
        do k=kscbase+1,ksctop  ! typically kscbase=3 & ksctop=6
         do iq=1,ifull
          thee(iq,k)=prcpv(k)*t(iq,k)*(t(iq,k) + .5*hlcp*qs(iq,k))
     .                         /(t(iq,k) - .5*hlcp*qs(iq,k))
          if(qg(iq,kscbase).gt.tied_rh*qs(iq,kscbase).
     .      and.thee(iq,k).lt.theeb(iq))then         !  default tied_rh=.75
            if(kt_sav(iq).eq.kl)then
              kbsav(iq)=kscbase
              ktsav(iq)=k
!             rk_shal(iq,k-1)=tied_con                
!             rk_shal(iq,k)=tied_over                 
	     endif
          endif ! (qg(iq,kscbase).gt......)
         enddo  ! iq loop
        enddo   ! end of k=kscbase+1,ksctop loop
       endif       ! (ksc.eq.79)
      endif  ! (ksc>70.and.ksc<80)
      
      if(ntest==-1.and.nproc==1.and.ktau.eq.10)then
        do k=1,kl
	  do iq=il*il+1,2*il*il
	   den1=qfg(iq,k)+qlg(iq,k)
	   if(land(iq).and.den1>0.)then
            write(26,'(4f9.4,2i6)')
     .	     t(iq,k),qfg(iq,k)/den1,1000.*qfg(iq,k),1000.*qlg(iq,k),iq,k
          endif
	  enddo
	 enddo     
      endif  ! (ntest==-1.and.nproc==1)    
      return
      end
