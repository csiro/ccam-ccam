      subroutine conjob      ! globpea & rcsb (non-chen); nkuo=46 only
      use cc_mpi, only : mydiag
      include 'newmpar.h' 
      parameter (itermax=1)  ! originally 2 - not available in vect. version
      parameter (ntest=0)    ! 1 or 2 to turn on; for 3 uncomment !!! lines
      parameter (nrenorm=1)  ! usually 1
      parameter (no2layer=0) ! usually 0, 1 to suppress 2-layer clouds (no good)
      parameter (kcl_top=kl-2) !max level for cloud top (conjob,radrive,vertmix)
!     N.B. usually set nkuo=44, rhcv=.75, rhmois=.6, sigcb=.9      
c     this one has safer code for nevapls=5
c     this one vectorized by jlm Tue  04-07-1998 (N.B. won't do chen)
c     nevapls, nevapcc:  turn off/on evap of ls or cc --- now through parm.h
c        0 off, 1 for Hal's evap, 2 for jlm, 3 for UK (ls only), 4 & 5 newer UK
c     Has dq(1) fix for cloud base drying
c     For 32-bit machine need real*8  cam()
c     Hal's ds() renamed dsh()
      include 'arrays.h'
      include 'const_phys.h'
      include 'dava.h' ! davt
      include 'kuocom.h'   ! also with kbsav,ktsav,convpsav,ndavconv
      include 'morepbl.h'
      include 'nlin.h'
      include 'parm.h'
      include 'prec.h'
      include 'sigs.h'
      include 'soil.h'
      include 'tracers.h'  ! ngas, nllp, ntrac
      common/epst/epst(ifull)
      common/work2/alphac(ifull),conrev(ifull),fluxr(ifull)
     .  ,heatpc(ifull),rhmx(ifull),smin(ifull),sumd(ifull)
     .  ,rnrt(ifull),rnrtc(ifull),sumqs(ifull),sumqt(ifull)
     .  ,sumss(ifull),factdav(ifull),kbsav_ls(ifull),dum2(ifull,4)
      common/work3/tt(ifull,kl),qq(ifull,kl),qsg(ifull,kl),gam(ifull,kl)
     .  ,spare(ifull,kl)
      common/work3b/ss(ifull,kl),uh(ifull,kl)
      common/work3c/dsh(ifull,kl)
      real dq(ifull,kl)
      equivalence (dq,ss)
      real*8 cam(kl-1,kl)
      real algf(kl),delt(kl),dsk(kl),rhs(kl),sigadd(kl),s(kl)
!     set convective relaxation time (convtime [was btcnv], usually 1 hour)
!     data convtime/1./,rhcv/.75/,rhmois/.6/ ! usually 1., .75, .6 now in kuocom
      data epsconv/0./

      include 'establ.h'
      Aev(tm) = 2.008e-9*tm**2 - 1.385e-6*tm + 2.424e-4  !For UKMO evap scheme
      Asb(tm) = max (-5.2e-9*tm**2+2.5332e-6*tm-2.9111e-4,1.e-5) !For UKMO subl
!!!   ntest=3   !!!

      do k=1,kl
       dsk(k)=-dsig(k)    !   dsk = delta sigma (positive)
       algf(k)=log(sig(k))
      enddo ! ! k=1,kl

      do k=2,kl
       delt(k)=.5*roncp*(algf(k-1)-algf(k))
      enddo ! ! k=2,kl
!     prev. used tdt=2*dt for leapfrog scheme, tdt=dt for adjustment scheme
!     convtime (hours) is convective relaxation time, originally 1 h
!     but zero (or <= dt) forces "hard" adjustment
!     bcnv=min( 1. , dt/max(1.e-20,convtime*3600.) )  ! is basically same as next line:
      bcnv=dt/max(dt,convtime*3600.) 
      ecmwa1=5.44e-4
      ecmwa2=5.09
      ecmwcc=sqrt(.05)
      ecmwls=1.

c     loop over all points : L/S rainfall, then conv
      do iq=1,ifull
       conrev(iq)=1000.*ps(iq)/(grav*dt) ! factor to conv. precip to g/m2/s
       rnrt(iq)=0.  ! initialize large-scale rainfall array
       rnrtc(iq)=0. ! initialize convective rainfall array
       kbsav_ls(iq)=0    ! for L/S
      enddo  ! iq loop

      do k=1,kl
       do iq=1,ifull
         tt(iq,k)=t(iq,k)
         qq(iq,k)=qg(iq,k)
         es=establ(tt(iq,k))
         pk=ps(iq)*sig(k)
         qsg(iq,k)=max(.622*es/(pk-es),1.5e-6)  ! Sat  10-31-1998
         gam(iq,k)=max(hlcp*qsg(iq,k)*pk*hlars/(tt(iq,k)**2*(pk-es)),
     .                 0.)   ! because establ poor for stratosphere   29/11/00
       enddo  ! iq loop
      enddo   ! k=1,kl

      if(diag.and.mydiag)then
        print *,'near beginning of conjob; idjd = ',idjd
        print *,'itermax,nrenorm,no2layer ',itermax,nrenorm,no2layer
        write (6,"('rh  ',12f7.2/(8x,12f7.2))") 
     .             (100.*qq(idjd,k)/qsg(idjd,k),k=1,kl)
        write (6,"('qs ',12f7.3/(8x,12f7.3))") 
     .             (1000.*qsg(idjd,k),k=1,kl)
        write (6,"('qg ',12f7.3/(8x,12f7.3))") 
     .             (1000.*qq(idjd,k),k=1,kl)
        write (6,"('tt ',12f7.2/(8x,12f7.2))") 
     .             (tt(idjd,k),k=1,kl)
        print *,'gam ',(gam(idjd,k),k=1,kl)
      endif
      if(rhcv.lt.0.)go to 5

!__________________________beginning of large-scale calculations_____________________
c     check for grid-scale rainfall first.
2     do k=kl,1,-1    ! top down to end up with proper kbsav_ls
       do iq=1,ifull
        if(qq(iq,k).gt.rhsat*qsg(iq,k))then
          kbsav_ls(iq)=k
          dqrx=(qq(iq,k)-rhsat*qsg(iq,k))/(1.+rhsat*gam(iq,k))
          tt(iq,k)=tt(iq,k)+hlcp*dqrx
          qq(iq,k)=qq(iq,k)-dqrx
          rnrt(iq)=rnrt(iq)+dqrx*dsk(k)
        endif   ! (qq(iq,k).gt.rhsat*qsg(iq,k))
       enddo    ! iq loop
      enddo     ! k loop

c     now do evaporation of L/S precip
      do iq=1,ifull
c      Here rnrt is the rainfall rate in gm/m**2/sec
!      conrev=1000.*ps(iq)/(grav*dt)
       rnrt(iq)=rnrt(iq)*conrev(iq) 
       fluxr(iq)=rnrt(iq)*1.e-3*dt ! kg/m2
      enddo  ! iq loop

      if(diag.and.mydiag)then
        print *,'before evap of large scale rain: kbsav_ls,rnrt ',
     .                                   kbsav_ls(idjd),rnrt(idjd)
        write (6,"('qg ',12f7.3/(8x,12f7.3))") 
     .             (1000.*qq(idjd,k),k=1,kl)
        write (6,"('tt ',12f7.2/(8x,12f7.2))") 
     .             (tt(idjd,k),k=1,kl)
      endif

      if(nevapls.eq.1)then
        do k=2*kl/3,1,-1
         do iq=1,ifull
          if(k.lt.kbsav_ls(iq))then
            revq=dt*ecmwls*ecmwa1*max((qsg(iq,k)-qq(iq,k)),0.)
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
     .        3.*(sig(k)-sig(kbsav_ls(iq)))*(1.-qq(iq,k)/qsg(iq,k))
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
            alphal=hl*qsg(iq,k)/(ars*tt(iq,k)**2)
            bl=1+Cev*dt*(1+hlcp*alphal)
            cfls=1. ! cld frac large scale
            evapls= cfls*dt*(Cev/bl)*(qsg(iq,k)-qq(iq,k)) !UKMO
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
            alphal=hl*qsg(iq,k)/(ars*tt(iq,k)**2)
            bl=1+Cev*dt*(1+hlcp*alphal)
            cfls=1. ! cld frac large scale
            evapls= cfls*dt*(Cev/bl)*(qsg(iq,k)-qq(iq,k)) !UKMO
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
        Dva=2.21
        cfls=1. ! cld frac7 large scale
        cflscon=4560.*cfls**.3125
        do k=2*kl/3,1,-1
         do iq=1,ifull
          if(k.lt.kbsav_ls(iq))then
            rhodz=ps(iq)*dsk(k)/grav
            qpf=fluxr(iq)/rhodz     ! Mix ratio of rain which falls into layer
            pk=ps(iq)*sig(k)
            es=qsg(iq,k)*pk/.622
            Apr=hl*hl/(rKa*tt(iq,k)*(rvap*tt(iq,k))-1.)
            Bpr=rvap*tt(iq,k)*pk/(Dva*es)
            Fr=fluxr(iq)/(cfls*dt)
            rhoa=pk/(rdry*tt(iq,k))
            dz=pk/(rhoa*grav)
            Vr=max(.01 , 11.3*Fr**(1./9.)/sqrt(rhoa)) ! Actual fall speed
            dtev=dz/Vr
            qr=fluxr(iq)/(dt*rhoa*Vr)
            qgdiff=qsg(iq,k)-qq(iq,k)
            Cev2=cflscon*qgdiff/(qsg(iq,k)*(Apr+Bpr)) ! Ignore rhoa**0.12
            qr2=max(0. , qr**.3125 - .3125*Cev2*dtev)**3.2
!           Cev=(qr-qr2)/(dtev*qgdiff)
            Cevx=(qr-qr2)/dtev  ! i.e. Cev*qgdiff
            alphal=hl*qsg(iq,k)/(ars*tt(iq,k)**2)
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
      if(rhcv.lt.0.)go to 8

!__________________________end of large-scale calculations_____________________

4     do k=1,kl         ! recompute qsg & gam, in case tt changed
       do iq=1,ifull
         es=establ(tt(iq,k))
         pk=ps(iq)*sig(k)
         qsg(iq,k)=max(.622*es/(pk-es),1.5e-6)  ! Sat  10-31-1998
         gam(iq,k)=max(hlcp*qsg(iq,k)*pk*hlars/(tt(iq,k)**2*(pk-es)),
     .                 0.)   ! because establ poor for stratosphere   29/11/00
       enddo  ! iq loop
      enddo   ! k=1,kl

      if(diag.and.mydiag)then
        print *,'before convection: rnrt,rnrtc ',
     .                                rnrt(idjd),rnrtc(idjd)
        write (6,"('rh   ',12f7.2/(8x,12f7.2))") 
     .             (100.*qq(idjd,k)/qsg(idjd,k),k=1,kl)
        write (6,"('qs   ',12f7.3/(8x,12f7.3))") 
     .             (1000.*qsg(idjd,k),k=1,kl)
        write (6,"('qq   ',12f7.3/(8x,12f7.3))") 
     .             (1000.*qq(idjd,k),k=1,kl)
        write (6,"('tt   ',12f7.2/(8x,12f7.2))") 
     .             (tt(idjd,k),k=1,kl)
      endif
      if(rhcv.lt.0.)go to 2

!__________________________beginning of convective calculations_____________________

5     do iq=1,ifull
       kbsav(iq)=0   ! preset value for no convection
       ktsav(iq)=kcl_top
       convpsav(iq)=0.
       sumss(iq)=0.
       sumqs(iq)=0.
       sumqt(iq)=0.
       dsh(iq,1)=0.
      enddo  ! iq loop

c     Multi-level penetrating convection : H.B.Gordon 1991
c     as used in the CSIRO4 and CSIRO9 spectral models.
c     The scheme is originally a "soft" moist adjustment generating
c     a mass flux. Convective columns are adjusted towards a
c     uniform moist profile.

c     The following cloud base/top as per BFR test model.
c     convection cloud base set at kuocb

c     compute instability parameters for all levels.
      do k=kcl_top,kuocb+1,-1   ! downwards to end just above cloud base
       if(epsconv.eq.0.)then  ! original conjob code
         do iq=1,ifull
c         dry static energy ss() [=(s(k)-s(k-1))/cp]
          ss(iq,k)=max(0.,                                             ! jlm
     .               tt(iq,k)*(1.+delt(k))-tt(iq,k-1)*(1.-delt(k)))
c         moist instability parameters
          uh(iq,k)=max(hlcp*(qq(iq,k-1)-qsg(iq,k))-ss(iq,k),0.)
c         if the rh value is insufficient, reset uh() value to 0.
          if(qq(iq,k-1).lt.abs(rhcv)*qsg(iq,k-1)) uh(iq,k)=0.
!         create an indicator for cloud base (if any)
          if(uh(iq,k).gt.0.)kbsav(iq)=k-1
         enddo  ! iq loop
       else   ! alternative code with different lowest base depending on epst
         do iq=1,ifull
c         dry static energy ss()
          ss(iq,k)=max(0.,                                             ! jlm
     .               tt(iq,k)*(1.+delt(k))-tt(iq,k-1)*(1.-delt(k)))
c         moist instability parameters
          uh(iq,k)=max(hlcp*(qq(iq,k-1)-qsg(iq,k))-ss(iq,k),0.)
c         if the rh value is insufficient, reset uh() value to 0.
          if(qq(iq,k-1).lt.abs(rhcv)*qsg(iq,k-1)) uh(iq,k)=0.
!         create an indicator for cloud base (if any)
c         if(uh(iq,k).gt.0.)kbsav(iq)=k-1
!         special: for min sea-cloud-base 1 layer above land-base
c         if(land(iq))then
c	    if(uh(iq,k).gt.0.)kbsav(iq)=k-1
c	  else  ! for sea:
cc          if(k.gt.kuocb+1.and.uh(iq,k).gt.0.)kbsav(iq)=k-1  ! land1
c	    if(k.gt.kuocb+2.and.uh(iq,k).gt.0.)kbsav(iq)=k-1  ! land2
c	  endif
!	  special: for usual cloud-base 1 layer above high-orog-base
!         to get just land-sea contrast use small epsconv, e.g. .0000001
 	  if(epst(iq).gt.epsconv)then  ! high-orog:
 	    if(uh(iq,k).gt.0.)kbsav(iq)=k-1
 	  else			! non-high-orog, lowest cloud base at kuocb+1:
	    if(k.gt.kuocb+1.and.uh(iq,k).gt.0.)kbsav(iq)=k-1
	   endif
         enddo  ! iq loop
       endif    !(epsconv.eq.0.)
      enddo     ! k loop

c!     special: for min sea-cloud-base 1 layer above land-base
c      do iq=1,ifull
c       if(.not.land(iq).and.kbsav(iq).gt.0)
c     .           kbsav(iq)=max(kbsav(iq),kuocb+1)
c      enddo

      if(diag.and.mydiag)then
        phi=bet(1)*tt(idjd,1)
        do k=2,kuocb
         phi=phi+bet(k)*tt(idjd,k)+betm(k)*tt(idjd,k-1)
	  s(k)=tt(idjd,k)+phi/cp
	 enddo
        do k=kuocb+1,kcl_top
         s(k)=s(k-1)+ss(idjd,k)  ! to give "usual" s
	 enddo
        write (6,"('s/cp ',12f7.2/(8x,12f7.2))") (s(k),k=1,kl)
        write (6,"('h/cp ',12f7.2/(8x,12f7.2))") 
     .             (s(k)+hlcp*qq(idjd,k),k=1,kl)
        write (6,"('hs/cp',12f7.2/(8x,12f7.2))") 
     .             (s(k)+hlcp*qsg(idjd,k),k=1,kl)
        write (6,"('uh a ',12f7.2/(8x,12f7.2))") (uh(idjd,k),k=1,kl)
      endif

c     check for conv at each point from level=kuocb up
c     test for convection: uh(iq,k+1)=(h(kbase)-h*(k+1))/cp >0
c     from computed instability parameters uh(iq,k)
c     set min cloud base for conv to be from level kuocb
c     convection control by rh>rhcv
c     stability parameter : uh(iq,k+1)=(h(kbase)-h(k+1))/cp
c     uh(iq,k+1)=max(hlcp*(qg(k)-qsg(iq,k+1))-ss(iq,k+1),0.)   ss?
c     instability shown by kbsav>0
!     don't need test here for kbsav>0, because values only used later within clouds
      do k=kuocb+1,kcl_top
       do iq=1,ifull
        if(k.gt.kbsav(iq).and.k.le.ktsav(iq))then
c         run up the levels stopping when no instability
          if(hlcp*(qq(iq,kbsav(iq))-qsg(iq,k)).gt.
     .                                       sumss(iq)+ss(iq,k))then
            sumss(iq)=sumss(iq)+ss(iq,k)
            uh(iq,k)=hlcp*(qq(iq,kbsav(iq))-qsg(iq,k))-sumss(iq)
          else   !  "new" uh(iq,k) .le.0., i.e. stable level reached
            ktsav(iq)=k-1
            If(no2layer.eq.1.and.ktsav(iq).eq.kbsav(iq)+1)
!    .                           ktsav(iq)=kbsav(iq)    
     .                           kbsav(iq)=0                ! 6/2/01
          endif ! (hlcp*.....gt.sumss(iq)+ss(iq,k))
        endif   ! (k.gt.kbsav(iq).and.k.le.ktsav(iq))
       enddo    ! iq loop
      enddo     ! k loop
      if(diag.and.mydiag)then
        print *,"after stab test: kbsav,ktsav,sumss ",
     .                            kbsav(idjd),ktsav(idjd),sumss(idjd)
        write (6,"('uh b ',12f7.2/(8x,12f7.2))") (uh(idjd,k),k=1,kl)
        print *,'ss ',(ss(idjd,k),k=1,ktsav(idjd))
        print *,'gam ',(gam(idjd,k),k=1,ktsav(idjd))
      endif

      do k=kuocb,kcl_top
       do iq=1,ifull
        if(k.ge.kbsav(iq).and.k.le.ktsav(iq))then
          sumqt(iq)=sumqt(iq)+qq(iq,k)*dsk(k)
          sumqs(iq)=sumqs(iq)+qsg(iq,k)*dsk(k)
          if(ntest.ne.0.and.iq.eq.idjd.and.mydiag)          
     .       print *,'k,sumqs,sumqt ',k,sumqs(idjd),sumqt(idjd)
        endif   ! (k.ge.kbsav(iq).and.k.le.ktsav(iq))
       enddo    ! iq loop
      enddo     ! k loop
      do iq=1,ifull
c      set up moisture changes
c      general moistening of driest levels : cutoff at rhmois
!      moisten dry levels towards rhmx, dry moist levels towards rhmx (jlm)
!      jlm: larger rhmois gives more moistening and less convective rain
!     don't need following "if" because values only used later within clouds
!      if(kbsav(iq).gt.0)then   
         rhmx(iq)=min(rhmois,.9*sumqt(iq)/sumqs(iq))
         temt=uh(iq,ktsav(iq))/(1.+gam(iq,ktsav(iq)))
         qclt=qsg(iq,ktsav(iq))+(gam(iq,ktsav(iq))*temt/hlcp)
         alphac(iq)=(qq(iq,kbsav(iq))-qclt)/
     .                               (sumqt(iq)-rhmx(iq)*sumqs(iq))
         heatpc(iq)=temt+sumss(iq)   ! this is  total heating
!      endif   ! (kbsav(iq).gt.0)then    ! don't need this "if" 10/9/99
      enddo    ! iq loop

      do k=kuocb,kcl_top
       do iq=1,ifull
!       if(k.ge.kbsav(iq))then          ! can omit this "if"
          dq(iq,k)=alphac(iq)*(qq(iq,k)-rhmx(iq)*qsg(iq,k))
!       endif  ! (k.ge.kbsav(iq))then   ! can omit this "if"
       enddo   ! iq loop
      enddo    ! k loop
      if(diag.and.mydiag)then
!       N.B. heatpc diag only meaningful if kbsav(iq)>0
!       dq & dq_eff are +ve for drying of environment
        print *,'rhmx,alphac,heatpc ',
     .           rhmx(idjd),alphac(idjd),heatpc(idjd)
        print *,'sumqt,sumqs ',sumqt(idjd),sumqs(idjd)
        print *,'dq_eff ',(dq(idjd,k)*dsk(k),k=1,ktsav(idjd))
      endif

c     from here, don't need ss, but need dsh(3d), dq(3d) and uh(3d)
      do iq=1,ifull
       kbase=kbsav(iq)
       if(kbase.gt.0)then
         ktop=ktsav(iq)
c        mc=number of levels involved in convection (>=2)
         mc=ktop+1-kbase
c        First guess heating at cloud base ds1'=0.
c        This determines the approx shape of the profile from
c        kbase+1 to ktsav(iq). If ds2' is found to be less than zero, then
c        recompute with ds1'' which is less than the ds2' .
         if(mc.eq.2)then
c          special case for mc=2 (no iterations)
c          mc=2 : convection from 1 level to next only - special coding
           dsh(iq,kbase+1)=heatpc(iq)/dsk(ktop)
         else
c          mc>=3 : convection for 3 or more levels
c          the scheme requires cam(kbase to ktop-1, kbase to ktop)
c          set up matrix cam(,) for solution by back substitution
c           'tr' elimination terms

           do k=kbase+2,ktop-1   ! N.B. mc=ktop+1-kbase
            rhs(k-1)=-hlcp*dq(iq,kbase)*(uh(iq,k+1)-uh(iq,kbase+1)) ! rhs
            cam(k-1,kbase)=uh(iq,k+1)*(1.+gam(iq,kbase+1)+delt(kbase+1))
     .            -uh(iq,kbase+1)*(delt(kbase+1)+delt(kbase+2))
            cam(k-1,k)=-uh(iq,kbase+1)*(1.+gam(iq,k+1)+delt(k+1))
!           do kic=kbase,k-2     ! not needed in newest code
!            cam(kic,k)=0.       ! not needed in newest code
!           enddo  ! kic loop    ! not needed in newest code
            do kic=k,ktop-1
             cam(kic-1,k-1)=-uh(iq,kbase+1)*(delt(k)+delt(k+1))
            enddo        ! kic loop
           enddo  ! k loop

           cam(kbase,kbase)=
     .            uh(iq,kbase+2)*(1.+gam(iq,kbase+1)+delt(kbase+1))
     .           -uh(iq,kbase+1)*(delt(kbase+1)+delt(kbase+2))
           cam(kbase,kbase+1)=-uh(iq,kbase+1)*
     .            (1.+gam(iq,kbase+2)+delt(kbase+2))
           rhs(kbase)=-hlcp*dq(iq,kbase)*(uh(iq,kbase+2)-uh(iq,kbase+1)) ! rhs
c          energy conservation terms
           do k=kbase+1,ktop
            cam(ktop-1,k-1)=dsk(k)
           enddo   ! k loop
           rhs(ktop-1)=heatpc(iq)
           if(ntest.ge.3.and.ntest.eq.mc)then
c            do k=kbase+2,ktop-1   ! N.B. mc=ktop+1-kbase
c             do kic=kbase,k-2     ! not needed in newest code - diag here
c              cam(kic,k)=0.       ! not needed in newest code - diag here
c             enddo  ! kic loop    ! not needed in newest code - diag here
c            enddo  ! k loop
             print *,'conjob cam_in for iq,kbase,ktop: ',iq,kbase,ktop
             print 93,(rhs(kic),kic=kbase,ktop-1)
             do k=ktop-1,kbase,-1
             print 93,(cam(kic,k),kic=kbase,ktop-1)
93           format(10f8.3)
            enddo
           endif

c          solve this matrix set of equations; rhs is in cam(...,ktop)
           do kic=ktop-1,kbase+1,-1
            x=cam(kic,kic)
            y=cam(kic-1,kic)
            do k=kbase,kic       ! with this, can get rid of zero cam above
             cam(kic-1,k)=cam(kic-1,k)*x-cam(kic,k)*y
            enddo  ! k loop
            rhs(kic-1)=rhs(kic-1)*x-rhs(kic)*y      ! rhs term
           enddo   ! kic loop

c          now substitute into the triangular cam
           dsh(iq,kbase+1)=rhs(kbase)/cam(kbase,kbase)
           do kic=kbase+1,ktop-1
            sum=rhs(kic)
            do k=kbase,kic-1
             sum=sum-cam(kic,k)*dsh(iq,k+1)
            enddo  ! k loop
            dsh(iq,kic+1)=sum/cam(kic,kic)
           enddo   ! kic loop
           if(ntest.ge.3.and.ntest.eq.mc)then
             print *,'conjob cam_out for iq: ',iq
             print 93,(rhs(kic),kic=kbase,ktop-1)
             do k=ktop-1,kbase,-1
             print 93,(cam(kic,k),kic=kbase,ktop-1)
             enddo
             print *,'conjob dsh_out (kbase value not used)'
             print 93,(dsh(iq,k),k=kbase,ktop)
!!!          ntest=ntest+1     !!!
           endif
         endif   ! (mc.eq.2)  ... else ...
       endif  ! (kbase.gt.0)
      enddo   ! iq loop

      if(nrenorm.eq.1)then  ! re-normalize
!       renormalize the profile to ensure >= 0 heating rates for all levels
        do iq=1,ifull
         smin(iq)=0.
         sumd(iq)=0.
        enddo   ! iq loop
        do k=kuocb+1,kcl_top
         do iq=1,ifull
          if(k.gt.kbsav(iq).and.k.le.ktsav(iq))then
            smin(iq)=min( smin(iq) , dsh(iq,k) )
            sumd(iq)=sumd(iq)+dsk(k)
          endif   !  (k.gt.kbsav(iq).and.k.le.ktsav(iq))
         enddo   ! iq loop
        enddo    ! k loop
        if(diag.and.mydiag)print *,
     .            'renorm smin,heatpc,sumd,alxx: ',
     .                    smin(idjd),heatpc(idjd),sumd(idjd),
     .             heatpc(idjd)/(heatpc(idjd)-smin(idjd)*sumd(idjd))
        do k=kuocb+1,kcl_top
         do iq=1,ifull
          if(smin(iq).lt.0..and.kbsav(iq).gt.0.and.k.gt.kbsav(iq))then
            alxx=heatpc(iq)/(heatpc(iq)-smin(iq)*sumd(iq))
            if(ntest.ne.0.and.iq.eq.idjd.and.mydiag)print *,
     .            'renorm k,dsh_in,dsh_out: ',
     .                    k,dsh(iq,k),alxx*(dsh(idjd,k)-smin(idjd))
           dsh(iq,k)=alxx*(dsh(iq,k)-smin(iq))
          endif  ! (smin(iq).lt.0..and.kbsav(iq).gt.0.and.k.gt.kbsav(iq))
         enddo   ! iq loop
        enddo    ! k loop
      endif      ! (nrenorm.eq.1)

c     determine the mass flux : use the uh(iq,kbsav+1) exp decay equation
      do iq=1,ifull
       if(kbsav(iq).gt.0)then
         kbase=kbsav(iq)
         denom=hlcp*dq(iq,kbase)+
     .      (1.+gam(iq,kbsav(iq)+1)+delt(kbsav(iq)+1))*dsh(iq,kbase+1)
!    .     -(1.-delt(kbsav(iq)+1))*dsh(iq,kbase) ! dsh(iq,kbase) stays zero now
         convpsav(iq)=bcnv*uh(iq,kbsav(iq)+1)/denom
c        if(denom.ne.0.)convpsav(iq)=bcnv*uh(iq,kbsav(iq)+1)/denom
       endif   ! (kbsav(iq).gt.0)
      enddo    ! iq loop
!      N.B. convpsav(iq) is already mult by dt jlm: mass flux is convpsav/dt

c     add the convection changes to the temp and moisture fields.
c     calc total precip.
c     Add to determine additive rainfall
      do k=kuocb,kcl_top
       do iq=1,ifull
        if(kbsav(iq).gt.0.and.k.ge.kbsav(iq).and.k.le.ktsav(iq))then
          dsh(iq,k)=convpsav(iq)*dsh(iq,k)
          dq(iq,k)=convpsav(iq)*dq(iq,k)
c         rnrtc(iq)=sum(dq(iq,i)*dsk(k))*100.*pg(mg)/(grav*dt)  ! kgm/m**2/sec
c         rnrtc(iq)=1000.*rnrtc(iq)                           !  gm/m**2/sec
          rnrtc(iq)=rnrtc(iq)+dq(iq,k)*dsk(k)*conrev(iq)      !  gm/m**2/sec
          qq(iq,k)=qq(iq,k)-dq(iq,k)
          tt(iq,k)=tt(iq,k)+dsh(iq,k)
        endif   !  (kbsav(iq).gt.0.and.k.ge.kbsav(iq).and.k.le.ktsav(iq))
       enddo    ! iq loop
      enddo     ! k loop

      if(diag.and.mydiag)then
        print *,"after Hal's convection: ktau,kbsav,ktsav,convpsav ",
     .                       ktau,kbsav(idjd),ktsav(idjd),convpsav(idjd)
        print *,'rhmx,sumss ',rhmx(idjd),sumss(idjd)
        do k=1,kbsav(idjd)-1
         dq(idjd,k)=0.   ! just for clean diagnostic
        enddo
        write (6,"('1000*dq ',12f6.3/(8x,12f6.3))") 
     .             (1000.*dq(idjd,k),k=1,ktsav(idjd))
        write (6,"('qg ',12f7.3/(8x,12f7.3))") 
     .             (1000.*qq(idjd,k),k=1,kl)
        write (6,"('tt      ',12f6.1/(8x,12f6.1))") 
     .             (tt(idjd,k),k=1,kl)
c       print *,'dq ',(dq(idjd,k),k=1,ktsav(idjd))  ! really dq
c       print *,'tt ',(tt(idjd,k),k=1,kl)
        print *,'C rnrt,rnrtc ',
     .             rnrt(idjd),rnrtc(idjd)
	 iq=idjd
        delq_av=0.
        delt_av=0.
	 heatlev=0.
	 do k=kbsav(idjd),ktsav(idjd)
	  delq_av=delq_av-dsk(k)*dq(iq,k)
	  delt_av=delt_av+dsk(k)*dsh(iq,k)
	  heatlev=heatlev+sig(k)*dsk(k)*dsh(iq,k)
	 enddo
	 print *,'delq_av,delt_exp,rnd_av,rnd_exp ',
     .    delq_av,-delq_av*hl/cp,rnrtc(iq),-delq_av*conrev(iq)
	 if(delt_av.ne.0.)print *,'ktau,kbsav,ktsav,delt_av,heatlev ',
     .          ktau,kbsav(idjd),ktsav(idjd),delt_av,heatlev/delt_av
      endif

c     now do evaporation of conv. precip
      if(nevapcc.eq.1)then
c       ECMWF based evaporation of convective rainfall
c       passing through lower layers. Convective rain has an assumed
c       grid coverage of Cc = 5%. The evap rate is controlled by
c       E = Cc.a1.(Qsat-Q) ( sqrt(sigma) * rnrtc / a2 / Cc)**a3
c       where a1=5.44E-04, a2=5.09E-03, and a3 = 0.5777
c       Here rnrtc is the rainfall rate in gm/m**2/sec.
c       See "Research manual 3" from ECMWF.
c       Calculate evap from top down : The evap may become
c       equal to the total rain (integrated down to that level) before
c       reaching the next level (or the surface).
        do k=kl/3,1,-1
         do iq=1,ifull
          if(k.lt.kbsav(iq))then
            revq=dt*ecmwcc*ecmwa1*max((qsg(iq,k)-qq(iq,k)),0.)*
     &         sqrt(sqrt(sig(k))*rnrtc(iq)/ecmwa2)
            revq=min(revq,rnrtc(iq)/(conrev(iq)*dsk(k)))
!           max needed for roundoff
            rnrtc(iq)=max(1.e-10,rnrtc(iq)-revq*dsk(k)*conrev(iq))
            tt(iq,k)=tt(iq,k)-revq*hlcp
            qq(iq,k)=qq(iq,k)+revq
          endif   !  (k.lt.kbsav(iq))
         enddo    ! iq loop
        enddo     ! k loop
      endif       ! (nevapcc.eq.1)

      if(nevapcc.eq.2)then  
        do k=kl/3,1,-1
         do iq=1,ifull
          if(k.lt.kbsav(iq))then
            revq=
     .      3.*(sig(k)-sig(kbsav(iq)))*(1.-qq(iq,k)/qsg(iq,k))
     .           *rnrtc(iq)/(conrev(iq)*.1)   ! jlm suggestion
            revq=min(revq,rnrtc(iq)/(conrev(iq)*dsk(k)))
!           max needed for roundoff
            rnrtc(iq)=max(1.e-10,rnrtc(iq)-revq*dsk(k)*conrev(iq))
            tt(iq,k)=tt(iq,k)-revq*hlcp
            qq(iq,k)=qq(iq,k)+revq
          endif   !  (k.lt.kbsav(iq))
         enddo    ! iq loop
        enddo     ! k loop
      endif       ! (nevapcc.eq.2)

      if(nevapcc.eq.3)then
        frac=.10   ! e.g. 20%   applied below cloud base
        sigadd(1)=dsk(1)
        do k=2,kl
         sigadd(k)=sigadd(k-1)+dsk(k)  ! where dsk=-dsig
        enddo
        do k=1,kl/2
         do iq=1,ifull
          if(k.lt.kbsav(iq))then
            revq=frac*rnrtc(iq)/(sigadd(kbsav(iq)-1)*conrev(iq))
            tt(iq,k)=tt(iq,k)-revq*hlcp
            qq(iq,k)=qq(iq,k)+revq
          endif   !  (k.lt.kbsav(iq))
         enddo    ! iq loop
        enddo     ! k loop
        do iq=1,ifull
         rnrtc(iq)=(1.-frac)*rnrtc(iq)
        enddo     ! iq loop
      endif       ! (nevapcc.eq.3)

      if(nevapcc.eq.4)then
        fract=.50   ! e.g. 50% from half-cloud down
        fracb=.00   ! e.g. 10% below base
        sigadd(1)=dsk(1)
        do k=2,kl
         sigadd(k)=sigadd(k-1)+dsk(k)  ! where dsk=-dsig
        enddo
        do k=1,kl-1
         do iq=1,ifull
          kav=(kbsav(iq)+ktsav(iq))/2
          if(kbsav(iq).gt.0.and.k.le.kav)then
            revq=fract*rnrtc(iq)/(sigadd(kav)*conrev(iq))
            tt(iq,k)=tt(iq,k)-revq*hlcp
            qq(iq,k)=qq(iq,k)+revq
          endif   !  (kbsav(iq).gt.0.and.k.le.ktsav(iq))
          if(k.lt.kbsav(iq))then
            revq=fracb*rnrtc(iq)/(sigadd(kbsav(iq)-1)*conrev(iq))
            tt(iq,k)=tt(iq,k)-revq*hlcp
            qq(iq,k)=qq(iq,k)+revq
          endif   !  (k.lt.kbsav(iq))
         enddo    ! iq loop
        enddo     ! k loop
        do iq=1,ifull
         rnrtc(iq)=(1.-fract-fracb)*rnrtc(iq)
        enddo     ! iq loop
      endif       ! (nevapcc.eq.4)
      if(rhcv.lt.0.)go to 4
!__________________________end of convective calculations_____________________

8     if(diag.and.mydiag)then
        print *,'near end of conjob: rnrt,rnrtc ',
     .                                  rnrt(idjd),rnrtc(idjd)
        write (6,"('qg ',12f7.3/(8x,12f7.3))") 
     .             (1000.*qq(idjd,k),k=1,kl)
        write (6,"('tt ',12f7.2/(8x,12f7.2))") 
     .             (tt(idjd,k),k=1,kl)
      endif

!     section for convective transport of trace gases (jlm 22/2/01)
      if(ilt.gt.1)then
        do ntr=1,ntrac
	  do k=1,klt
	   do iq=1,ilt*jlt
           ss(iq,k)=tr(iq,k,ntr)
          enddo
         enddo
         do iq=1,ilt*jlt
          if(kbsav(iq).gt.0)then
            kb=kbsav(iq)
            kt=ktsav(iq)
            veldt=convpsav(iq)
            fluxup=veldt*ss(iq,kb)
!           remove gas from cloud base layer
            tr(iq,kb,ntr)=tr(iq,kb,ntr)-fluxup/dsk(kb)
!           put flux of gas into top convective layer
            tr(iq,kt,ntr)=tr(iq,kt,ntr)+fluxup/dsk(kt)
            do k=kb+1,kt
             tr(iq,k,ntr)=tr(iq,k,ntr)-ss(iq,k)*veldt/dsk(k)
             tr(iq,k-1,ntr)=tr(iq,k-1,ntr)+ss(iq,k)*veldt/dsk(k-1)
            enddo
          endif
         enddo   ! iq loop
        enddo    ! ntr loop    
      endif      ! ilt.gt.1

!     usual conformal-cubic (DARLAM option removed)
        qg(1:ifull,:)=qq(:,:)
c       print *,'B rnrt,rnrtc,condx ',
c    .             rnrt(idjd),rnrtc(idjd),condx(idjd)
        do iq=1,ifull
         prcon= .001*dt*rnrtc(iq)
         condc(iq)=prcon              ! convective precip for this timestep
         precc(iq)=precc(iq)+prcon
         prl_s= .001*dt*rnrt(iq)
         condx(iq)=prcon+prl_s        ! total precip for this timestep
         precip(iq)=precip(iq)+prcon+prl_s
        enddo
        t(1:ifull,:)=tt(:,:)           ! usually split from June '03

      if(ntest.eq.1.and.mydiag)then
        print *,'D rnrt,rnrtc,condx ',
     .             rnrt(idjd),rnrtc(idjd),condx(idjd)
        print *,'precc,precip ',
     .           precc(idjd),precip(idjd)
      endif
      return
      end
