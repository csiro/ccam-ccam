      subroutine vertmix
!     inputs & outputs: t,u,v,qg
      use cc_mpi, only : mydiag
      use diag_m
      include 'newmpar.h'
      parameter (ntest=0)
c     parameter (ipwr=1)       ! really can use (ipwr=min(1,nlocal))
c     parameter (ilnl=il**ipwr,jlnl=jl**ipwr)
      parameter (kcl_top=kl-2) ! maximum level for cloud top (conjob & vertmix)
      parameter (ndvmod=0)     ! 0 default, 1+ for dvmod tests
!     typically tied_con=6., tied_over=2., tied_rh=.75
!     nlocal in parm.h         ! 0 local scheme, 1 or 2 non-local scheme
!     real t(ifull,kl),u(ifull,kl),v(ifull,kl),qg(ifull,kl),ps(ifull)
      include 'arrays.h'
      include 'dates.h'
      include 'indices.h'
      include 'kuocom.h'   ! also with kbsav,ktsav,convpsav,kscsea,sigksct
      include 'liqwpar.h'  ! ifullw, qfg, qlg
      include 'map.h'      ! em, f, fu, fv, etc  not needed here?
      include 'nlin.h'
      include 'morepbl.h'
      include 'parm.h'
      include 'pbl.h'
      include 'permsurf.h'
      include 'savuvt.h'
      include 'sigs.h'
      include 'soil.h'
      include 'tracers.h'  ! ngas, nllp, ntrac
      include 'xarrs.h'
      include 'const_phys.h'
      common/work3/qs(ifull,kl),delthet(ifull,kl),
     .             thebas(ifull,kl),cu(ifull,kl),dum3(ifull,kl)
      common/work3b/uav(ifull,kl),vav(ifull,kl)   
!     n.b. uav & vav also used by pbldif; all of work3 used by tracervmix
      common/work3c/rhs(ifull,kl)   
      common/work3f/wrk1(ijk),wrk2(ijk),wrk3(ijk) 
      common/work2/csq(ifull),delq(ifull),dvmod(ifull)
     . ,dz(ifull),dzr(ifull),fm(ifull),fh(ifull),ri(ifull),sqmxl(ifull)
     . ,x(ifull),zhv(ifull),theeb(ifull),dum2(ifull,6)
      real sighkap(kl),delons(kl),delh(kl),prcpv(kl)
      real sigkap(kl)
      real alfqq_l(kl),alfqq_s(kl),alfqq1_l(kl),alfqq1_s(kl)
      real at(ifull,kl),au(ifull,kl),ct(ifull,kl)
      real zh(ifull,kl),tmnht(ifull,kl)
      real gt(ifull,kl),guv(ifull,kl)
      real rkm(ifull,kl),rkh(ifull,kl)
!     equivalence (gt,rkh_nl),(guv,rkm_nl)
      equivalence (rkh,wrk1),(rkm,wrk2)
      equivalence (tmnht,at,un),(zh,au,wrk3)
!     equivalence (gamat,ct)
c     set coefficients for Louis scheme
      data bprm/4.7/,cm/7.4/,ch/5.3/,amxlsq/100./,vkar3/.35/,vkar4/.4/
      data bprmj/5./,cmj/5./,chj/2.6/
      save kscbase,ksctop,prcpv
      include 'establ.h'

      rong=rdry/grav
      do k=1,kl-1
       sighkap(k)=sigmh(k+1)**(-roncp)
       delons(k)=rong *((sig(k+1)-sig(k))/sigmh(k+1))
      enddo      ! k loop
      do k=1,kl
       delh(k)=-rong *dsig(k)/sig(k)  ! sign of delh defined so always +ve
       sigkap(k)=sig(k)**(-roncp)
      enddo      ! k loop
      if( (diag.or.ntest.ge.1).and.mydiag)then
        print *,'sig ',sig
        print *,'dsig ',dsig
        print *,'delh ',delh
	 print *,'in vertmix'
        write (6,"('uin ',19f7.2/(8x,19f7.2))") (u(idjd,k),k=1,kl) 
        write (6,"('vin ',19f7.2/(8x,19f7.2))") (v(idjd,k),k=1,kl) 
      endif
      rlogs1=log(sig(1))
      rlogs2=log(sig(2))
      rlogh1=log(sigmh(2))
      rlog12=1./(rlogs1-rlogs2)
      tmnht(1:ifull,1)=(t(1:ifull,2)*rlogs1-t(1:ifull,1)*rlogs2+
     .           (t(1:ifull,1)-t(1:ifull,2))*rlogh1)*rlog12
!     n.b. an approximate zh is quite adequate for this routine
      zh(1:ifull,1)=t(1:ifull,1)*delh(1)
      do k=2,kl-1
       do iq=1,ifull
        zh(iq,k)=zh(iq,k-1)+t(iq,k)*delh(k)
        tmnht(iq,k) =(t(iq,k)+t(iq,k+1))*.5
       enddo     ! iq loop
      enddo      !  k loop
      do k=1,kl
       do iq=1,ifull
        rhs(iq,k)=t(iq,k)*sigkap(k)  ! rhs is theta here
        es=establ(t(iq,k))
        qs(iq,k)=.622*es/(ps(iq)*sig(k)-es)
       enddo     ! iq loop
      enddo      !  k loop
      do k=1,kl-1  ! top level set separately to 0 for gt & guv
       rkh(:,k)=0.
       rkm(:,k)=0.
       delthet(:,k)=rhs(:,k+1)-rhs(:,k)  ! rhs is theta here
      enddo      !  k loop

      if(ktau.eq.1)then
!       set ksctop for shallow convection
        ksctop=1    ! ksctop will be first level below sigkcst
        do while(sig(ksctop+1).gt.sigksct)  !  e.g. sigksct=.8
         ksctop=ksctop+1
        enddo
        kscbase=1  ! kscbase will be first level above sigkcsb
        do while(sig(kscbase).gt.sigkscb)  !  e.g. sigkscb=.99
         kscbase=kscbase+1
        enddo
        if ( myid == 0 ) then
        print *,'For shallow convection:'
        print *,'ksc,kscbase,ksctop,kscsea ',
     .           ksc,kscbase,ksctop,kscsea
	 write (6,"(' sigkscb,sigksct,tied_con,tied_over,tied_rh:',
     .       5f8.3)")sigkscb,sigksct,tied_con,tied_over,tied_rh
        end if
        do k=1,kl
         prcpv(k)=sig(k)**(-rdry/cp)
        enddo    ! k loop
        do k=1,kl
         alfqq_s(k)=alfsea
         alfqq_l(k)=alflnd
         alfqq_s(k)=1./tied_rh
         alfqq_l(k)=1./tied_rh
         alfqq1_s(k)=0.
         alfqq1_l(k)=0.
        enddo     ! k loop
      endif

c     ************ section for Tiedtke shallow convection *******************
      if(ksc.eq.99)then
        do iq=1,ifull
        theeb(iq)=prcpv(kscbase)*t(iq,kscbase)*
     .                  (t(iq,kscbase) + .5*hlcp*qs(iq,kscbase))
     .                 /(t(iq,kscbase) - .5*hlcp*qs(iq,kscbase))
        enddo    ! iq loop
        if(kscsea.eq.1)then  ! Tiedtke done only over sea
          do k=kscbase+1,ksctop
           do ip=ipland+1,ipsea                                      ! sea only
            iq=iperm(ip)                                             ! sea only
            thee=prcpv(k)*t(iq,k)*(t(iq,k) + .5*hlcp*qs(iq,k))
     .                           /(t(iq,k) - .5*hlcp*qs(iq,k))
            if(qg(iq,kscbase).gt.tied_rh*qs(iq,kscbase).
     .       and.thee.lt.theeb(iq))then               !  default tied_rh=.75
              rkh(iq,k-1)=tied_con                !  m**2/sec  6., originally 10.
              rkh(iq,k)=tied_over                 !  m**2/sec
            endif ! (qg(iq,kscbase).gt.rhscon*qs(iq,kscbase).....
           enddo  ! iq loop
          enddo   ! end of k=kscbase+1,ksctop loop
        else      !  i.e. Tiedtke original scheme over land and sea
          do k=kscbase+1,ksctop  ! typically kscbase=3 & ksctop=6
           do iq=1,ifull
            thee=prcpv(k)*t(iq,k)*(t(iq,k) + .5*hlcp*qs(iq,k))
     .                           /(t(iq,k) - .5*hlcp*qs(iq,k))
            if(qg(iq,kscbase).gt.tied_rh*qs(iq,kscbase).
     .       and.thee.lt.theeb(iq))then               !  default tied_rh=.75
              rkh(iq,k-1)=tied_con                !  m**2/sec  6., originally 10.
              rkh(iq,k)=tied_over                 !  m**2/sec
		if(ntest.eq.3.and.k.eq.ksctop)then
		  print *,'ktau,iq,theeb,thee,delthee ',
     .       		    ktau,iq,theeb(iq),thee,theeb(iq)-thee
		endif
            endif ! (qg(iq,kscbase).gt.rhscon*qs(iq,kscbase).....
           enddo  ! iq loop
          enddo   ! end of k=kscbase+1,ksctop loop
        endif     ! (kscsea.eq.1)  .. else ..
      endif       ! (ksc.eq.99)
c     *********** end of Tiedtke shallow convection section *****************

c     ************ Tiedtke_jlm shallow convection 97 ***************
      if(ksc.eq.97)then
        do k=kscbase,ksctop-1
         do iq=1,ifull
	   qbas=min(qg(iq,k)/tied_rh,qs(iq,k))
c         thebas(iq,k)=t(iq,k)+hlcp*qbas  + hght
          thebas(iq,k)=prcpv(k)*t(iq,k)*(t(iq,k) + .5*hlcp*qbas)
     .                                 /(t(iq,k) - .5*hlcp*qbas)
         enddo  ! iq loop
        enddo   ! k loop
	 theeb(:)=thebas(:,kscbase)
        do k=kscbase+1,ksctop+1
         do iq=1,ifull
          thee=prcpv(k)*t(iq,k)*(t(iq,k) + .5*hlcp*qs(iq,k))
     .                         /(t(iq,k) - .5*hlcp*qs(iq,k))
          if(theeb(iq).gt.thee)then
            rkh(iq,k-1)=tied_con                !  m**2/sec  6., originally 10.
            rkh(iq,k)=tied_over                 !  m**2/sec
	   else
	     theeb(iq)=thebas(iq,k) ! ready for next k in k-loop
          endif
         enddo  ! iq loop
        enddo   ! k loop
!       suppress pseudo-deep convection	 
        do k=kscbase,ksctop
         do iq=1,ifull
          if(rkh(iq,ksctop).gt..9*tied_con)then  ! to allow for tied_over=0.
            rkh(iq,k)=0.
          endif
         enddo  ! iq loop
        enddo   ! k loop
      endif     ! (ksc.eq.97)
c     *********** end of Tiedtke_jlm shallow convection 97 *************

c     ************ section for Geleyn shallow convection *******************
      if(ksc.eq.-99)then
        do k=kscbase,ksctop    ! new usage of ksc thu  02-17-2000
         do iq=1,ifull
          delthet(iq,k)=delthet(iq,k)-hlcp*max(0.,
     .                 qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k) )
         enddo  ! iq loop
        enddo   !  k loop
      endif     ! (ksc.eq.-99)
c     *********** end of Geleyn shallow convection section *****************

      if( (ntest.ne.0.or.diag) .and. mydiag)then
        iq=idjd
        print *,'for shallow convection in vertmix '
	 print *,'kbsav,ktsav,theeb: ',kbsav(iq),ktsav(iq),theeb(iq)
        write (6,"('rh   ',19f7.2/(8x,19f7.2))") 
     .             (100.*qg(idjd,k)/qs(idjd,k),k=1,kl)
        write (6,"('qs   ',19f7.3/(8x,19f7.3))") 
     .             (1000.*qs(idjd,k),k=1,kl)
        write (6,"('qg   ',19f7.3/(8x,19f7.3))") 
     .             (1000.*qg(idjd,k),k=1,kl)
        write (6,"('qbas ',19f7.3/(8x,19f7.3))") 
     .             (1000.*qg(idjd,k)/tied_rh,k=1,kl)
        write (6,"('t    ',19f7.2/(8x,19f7.2))") 
     .             (t(idjd,k),k=1,kl)
        write (6,"('thebas',19f7.2/(8x,19f7.2))") 
     .             (thebas(iq,k),k=1,kl)
c        write (6,"('hs',19f7.2/(8x,19f7.2))") 
c     .             (t(idjd,k)+hlcp*qs(idjd,k),k=1,kl)
        write (6,"('thee',19f7.2/(8x,19f7.2))") 
     .            (prcpv(k)*t(iq,k)*(t(iq,k) + .5*hlcp*qs(iq,k))
     .                             /(t(iq,k) - .5*hlcp*qs(iq,k)),k=1,kl)
      endif

!     following now defined in vertmix (don't need to pass from sflux)
      uav(1:ifull,:)=av_vmod*u(1:ifull,:)+(1.-av_vmod)*savu(1:ifull,:)   ! 3D
      vav(1:ifull,:)=av_vmod*v(1:ifull,:)+(1.-av_vmod)*savv(1:ifull,:)   ! 3D
      do k=1,kl-1
       do iq=1,ifull
        dz(iq) =-tmnht(iq,k)*delons(k)  ! this is z(k+1)-z(k)
        dzr(iq)=1./dz(iq)
        zhv(iq)=1./zh(iq,k)
        if(ndvmod.eq.0)then
          dvmod(iq)=sqrt( (uav(iq,k+1)-uav(iq,k))**2
     .                   +(vav(iq,k+1)-vav(iq,k))**2 )
        else
          dvmod(iq)=ndvmod  ! just for tests
        endif    ! (ndvmod.eq.0)

c       x is bulk ri *(dvmod **2)
        x(iq)=grav*dz(iq)*(delthet(iq,k)/
     .   (tmnht(iq,k)*sighkap(k))+.61*(qg(iq,k+1)-qg(iq,k)))

c       fm and fh denote f(Louis style)*dvmod
c       nb. an error exists in the Louis expression for c; this is corrected
c       in the current code
        csq(iq) = zhv(iq)*
     .              (((1.+dz(iq)*zhv(iq))**(1./3.)-1.)*dzr(iq))**3
       enddo     ! iq loop

       if(nvmix.lt.0)then   ! use quasi-neutral mixing for test runs only
         do iq=1,ifull
          fm(iq)=1.
          fh(iq)=1.
         enddo   ! iq loop
         if(k.lt.kl/2)then
           do iq=1,ifull
            if(land(iq))fh(iq)=abs(nvmix)
           enddo   ! iq loop
         endif
       endif     !  (nvmix.lt.0)

       if(nvmix.eq.4)then
c        newest code, stable same as csiro9 here
         do iq=1,ifull
          sqmxl(iq)=(vkar4*zh(iq,k)/(1.+vkar4*zh(iq,k)/amxlsq))**2
          dvmod(iq)=max( dvmod(iq) , 1. )
          ri(iq)=x(iq)/dvmod(iq)**2
          if(ri(iq).lt. 0.)then  ! unstable case
c           first do momentum
            denma=dvmod(iq)+cmj*( 2.*bprmj*sqmxl(iq)*
     .                            sqrt(-x(iq)*csq(iq)) )
            fm(iq)=dvmod(iq)-(2.*bprmj *x(iq))/denma
c           now heat
            denha=dvmod(iq)+chj*( 2.*bprmj*sqmxl(iq)*
     .                            sqrt(-x(iq)*csq(iq)) )
            fh(iq)=dvmod(iq)-(2.*bprmj *x(iq))/denha
          else                     ! stable case
c           the following is the original Louis stable formula
            fm(iq)=dvmod(iq)/(1.+4.7*ri(iq))**2
            fh(iq)=fm(iq)
          endif
         enddo   ! iq loop
       endif     ! (nvmix.eq.4)

c      calculate k's, guv and gt
c      nonlocal scheme usually applied to temperature and moisture, not momentum
c      (i.e. local scheme is applied to momentum for nlocal=0,1)
       if(nvmix.ne.0)then    ! use nvmix=0 only for special test runs
         do iq=1,ifull
          rkm(iq,k)=rkm(iq,k)+fm(iq)*sqmxl(iq)*dzr(iq)
          rkh(iq,k)=rkh(iq,k)+fh(iq)*sqmxl(iq)*dzr(iq)
         enddo   ! iq loop
       endif     ! (nvmix.ne.0)

       if((diag.or.ntest.eq.2).and.mydiag)then
         print *,'k,dt,delsig,sqmxl,dzr ',
     .            k,dt,delsig,sqmxl(idjd),dzr(idjd)
         print *,'k,t,t+,ps ',k,t(idjd,k),t(idjd,k+1),ps(idjd)
         print *,'k,qg,qg+ ',k,qg(idjd,k),qg(idjd,k+1)
         print *,'k,qs,qs+ ',k,qs(idjd,k),qs(idjd,k+1)
         es=establ(t(idjd,k))
         esp=establ(t(idjd,k+1))
         print *,'k,es,es+,delthet ',k,es,esp,delthet(idjd,k)
         print *,'k,fm,dvmod,ri,csq ',
     .            k,fm(idjd),dvmod(idjd),ri(idjd),csq(idjd)
         print *,'x,zh,tmnht,dz,delh,sighkap ',x(idjd),zh(idjd,k),
     .                 tmnht(idjd,k),dz(idjd),delh(k),sighkap(k)
       endif
      enddo      ! end of k loop

      if( (diag.or.ntest.ge.1) .and. mydiag )then
        print *,'before possible call to pbldif in vertmix'
          write (6,"('rkh0 ',16f8.3)") (rkh(idjd,k),k=1,16)
          write (6,"('rkm0 ',16f8.3)") (rkm(idjd,k),k=1,16)
        write (6,"('uav ',19f7.2/(8x,19f7.2))") (uav(idjd,k),k=1,kl) 
        write (6,"('vav ',19f7.2/(8x,19f7.2))") (vav(idjd,k),k=1,kl) 
        write (6,"('thet',19f7.2/(8x,19f7.2))") (rhs(idjd,k),k=1,kl) 
        write (6,"('t   ',19f7.2/(8x,19f7.2))") (t(idjd,k),k=1,kl) 
        write (6,"('qg ',19f7.3/(8x,19f7.3))") 
     .             (1000.*qg(idjd,k),k=1,kl)
        write (6,"('qs ',19f7.3/(8x,19f7.3))") 
     .             (1000.*qs(idjd,k),k=1,kl)
        write (6,"('thee',19f7.2/(8x,19f7.2))") 
     .        (prcpv(k)*t(idjd,k)*(t(idjd,k) + .5*hlcp*qs(idjd,k))
     .                   /(t(idjd,k) - .5*hlcp*qs(idjd,k)),k=1,kl)
      endif

      if(nlocal.ne.0)then
        call pbldif(rhs,rkh,rkm,uav,vav)
!       n.b. *** pbldif partially updates qg and theta (t done during trim)	 
!       ncar info is returned in rkm_nl and rkh_nl arrays
        if( (diag.or.ntest.ge.1) .and. mydiag )then
	   print *,'after pbldif in vertmix'
          write (6,"('rkh1 ',16f8.3)") (rkh(idjd,k),k=1,16)
          write (6,"('rkm1 ',16f8.3)") (rkm(idjd,k),k=1,16)
          write (6,"('thet',19f7.2/(8x,19f7.2))") (rhs(idjd,k),k=1,kl) 
          write (6,"('qg ',19f7.3/(8x,19f7.3))") 
     .              (1000.*qg(idjd,k),k=1,kl)
        endif
        if(diag)then
          call printa('rkh ',rkh,ktau,nlv,ia,ib,ja,jb,0.,1.)
          call printa('cond',condx,ktau,1,ia,ib,ja,jb,0.,1.)
          call printa('zs  ',zs,ktau,1,ia,ib,ja,jb,0.,1.)
        endif
      endif      ! (nlocal.gt.0)

      do k=1,kl-1
        delsig=sig(k+1)-sig(k)
        do iq=1,ifull
          dz(iq) =-tmnht(iq,k)*delons(k)  ! this is z(k+1)-z(k)
          dzr(iq)=1./dz(iq)
          guv(iq,k)= rkm(iq,k)*dt *delsig *dzr(iq)**2
          gt(iq,k)=  rkh(iq,k)*dt *delsig *dzr(iq)**2
         enddo   ! iq loop
      enddo      ! k loop

      guv(:,kl)=0.
      gt(:,kl)=0.
      if(diag)then
        call maxmin(guv,'gu',ktau,1.e3,kl)
        if ( mydiag ) then
           print *,'vertmix guv ',(guv(idjd,k),k=1,kl)
           print *,'vertmix gt ',(gt(idjd,k),k=1,kl)
        end if
        call printa('t   ',t,ktau,nlv,ia,ib,ja,jb,200.,1.)
        call printa('tss ',tss,ktau,nlv,ia,ib,ja,jb,200.,1.)
        call printa('eg  ',eg,ktau,nlv,ia,ib,ja,jb,0.,1.)
        call printa('fg  ',fg,ktau,nlv,ia,ib,ja,jb,0.,1.)
        call printa('thet',rhs,ktau,nlv,ia,ib,ja,jb,200.,1.)
      endif

      if(ncvmix.gt.0)then  ! cumulus mixing of momentum - jlm version
        do k=kuocb+1-ncvmix,kcl_top-1
         do iq=1,ifull
!          for no-conv-points, doing k=1-ncvmix,1 with convpsav=0.
!          1 below for ncvmix=2
           if(kbsav(iq).gt.0.and.k.ge.kbsav(iq)+1-ncvmix
     .                      .and.k.lt.ktsav(iq))
     .       guv(iq,k)=guv(iq,k)-convpsav(iq)*.5      ! with factor of .5
          enddo  ! iq loop
        enddo    ! k loop
        if(diag.and.mydiag)then
          print *,'vertmix after conv; kb,kt,convp'
     .     ,kbsav(idjd),ktsav(idjd),convpsav(idjd)
          print *,'new guv',(guv(idjd,k),k=1,kl)
        endif
      endif      !   (ncvmix.gt.0)

      conflux=grav*dt/dsig(1)
      condrag=grav*dt/(dsig(1)*rdry)
c     first do theta (then convert back to t)
      at(:,1)=0.
      rhs(1:ifull,1)=rhs(1:ifull,1)-(conflux/cp)*fg(1:ifull)/ps(1:ifull)

      if(npanels.gt.0)then
        do k=2,kl
         do iq=1,ifull
          at(iq,k) =-gt(iq,k-1)/dsig(k)
         enddo   ! iq loop
        enddo    !  k loop
        do k=1,kl
*cdir nodep
         do iq=1,ifull
          ct(iq,k) =-gt(iq,k)/dsig(k)
         enddo   ! iq loop
        enddo    !  k loop
      else       ! i.e. npanels=0   darlam
        do k=2,kl
         do iq=1,ifull
          at(iq,k) =-gt(iq,k-1)/dsig(k)
         enddo   ! iq loop
        enddo    !  k loop
        do k=1,kl
*cdir nodep
         do iq=1,ifull
          ct(iq,k) =-gt(iq,k)/dsig(k)
         enddo   ! iq loop
        enddo    !  k loop
      endif      !  (npanels.gt.0) .. else ..
      if((diag.or.ntest.eq.2).and.mydiag)then
        print *,'ktau,fg,tss,ps ',ktau,fg(idjd),tss(idjd),ps(idjd)
        print *,'at ',(at(idjd,k),k=1,kl)
        print *,'ct ',(ct(idjd,k),k=1,kl)
        print *,'rhs ',(rhs(idjd,k),k=1,kl)
      endif      ! (ntest.eq.2)

      call trim(at,ct,rhs,0)   ! for t
      do k=1,kl
        do iq=1,ifull
         t(iq,k)=rhs(iq,k)/sigkap(k)
        enddo  ! iq loop
      enddo      !  k loop
      if(diag)then
         if (mydiag) then
            print *,'vertmix eg,fg,cdtq,land '
     .                  ,eg(idjd),fg(idjd),cdtq(idjd),land(idjd)
            print *,'vertmix theta after trim ',(rhs(idjd,k),k=1,kl)
         end if
        call printa('thet',rhs,ktau,nlv,ia,ib,ja,jb,200.,1.)
      endif

c     now do moisture
      rhs(1:ifull,:)=qg(1:ifull,:)
      rhs(1:ifull,1)=rhs(1:ifull,1)-(conflux/hl)*eg(1:ifull)/ps(1:ifull)
c     could add extra sfce moisture flux term for crank-nicholson
      call trim(at,ct,rhs,0)    ! for qg
      qg(1:ifull,:)=rhs(1:ifull,:)

      if(ifullw.eq.ifull)then
c       now do qfg
	 do k=1,kl
	  do iq=1,ifullw
          rhs(iq,k)=qfg(iq,k)
         enddo
        enddo
        call trim(at,ct,rhs,0)    ! for qfg
	 do k=1,kl
	  do iq=1,ifullw
          qfg(iq,k)=rhs(iq,k)
         enddo
        enddo
c       now do qlg
	 do k=1,kl
	  do iq=1,ifullw
          rhs(iq,k)=qlg(iq,k)
         enddo
        enddo
        call trim(at,ct,rhs,0)    ! for qlg
	 do k=1,kl
	  do iq=1,ifullw
          qlg(iq,k)=rhs(iq,k)
         enddo
        enddo
      endif    ! (ifullw.eq.ifull)

c     now do trace gases
      if(ngas.gt.0)call tracervmix( at, ct, rhs )

!     from here on just u and v stuff
      au(:,1)=cduv(:)*condrag/tss(:)   ! globpea
      do k=2,kl
       do iq=1,ifull
        au(iq,k) =-guv(iq,k-1)/dsig(k)  ! globpea
       enddo   ! iq loop
      enddo    !  k loop
      do k=1,kl
       do iq=1,ifull
        cu(iq,k) =-guv(iq,k)/dsig(k)    ! globpea
       enddo   ! iq loop
      enddo    !  k loop
      if((diag.or.ntest.eq.2).and.mydiag)then
        print *,'au ',(au(idjd,k),k=1,kl)
        print *,'cu ',(cu(idjd,k),k=1,kl)
      endif      ! (ntest.eq.2)

c     first do u
      rhs(1:ifull,:)=u(1:ifull,:)
      call trim(au,cu,rhs,0)
      u(1:ifull,:)=rhs(1:ifull,:)
      if(diag.and.mydiag)then
        print *,'vertmix au ',(au(idjd,k),k=1,kl)
      endif

c     now do v; with properly unstaggered au,cu
      rhs(1:ifull,:)=v(1:ifull,:)
      call trim(au,cu,rhs,0)    ! note now that au, cu unstaggered globpea
      v(1:ifull,:)=rhs(1:ifull,:)
      if((diag.or.ntest.ge.1).and.mydiag)then
        print *,'after trim in vertmix '
        write (6,"('thet',19f7.2/(8x,19f7.2))") 
     .              (sigkap(k)*t(idjd,k),k=1,kl) 
        write (6,"('t   ',19f7.2/(8x,19f7.2))") (t(idjd,k),k=1,kl) 
        write (6,"('qg ',19f7.3/(8x,19f7.3))") 
     .             (1000.*qg(idjd,k),k=1,kl)
        write (6,"('u   ',19f7.2/(8x,19f7.2))") (u(idjd,k),k=1,kl) 
        write (6,"('v   ',19f7.2/(8x,19f7.2))") (v(idjd,k),k=1,kl) 
        print *,'cduv,cduv+1,cduvj+1,tss ',
     .           cduv(idjd),cduv(idjd+1),cduv(idjd+il),tss(idjd)
        write (6,"('au ',19f7.3/(8x,19f7.3))") (au(idjd,k),k=1,kl) 
        write (6,"('cu ',19f7.3/(8x,19f7.3))") (cu(idjd,k),k=1,kl) 
c       call maxmin(u,'uv',ktau,1.,kl)
c       call maxmin(v,'vv',ktau,1.,kl)
      endif
      return
      end

      subroutine tracervmix( at, ct, updtr )
c     this routine does the vertical mixing of tracers
      use cc_mpi, only : mydiag
      use diag_m
      parameter (ntest=0)  ! 1 for diag prints
      include 'newmpar.h'
      include 'const_phys.h'
      include 'dates.h'         ! dt
      include 'nsibd.h'
      include 'parm.h'
      include 'sigs.h'          ! dsig
      include 'soil.h'
      include 'tracers.h'       ! jlm don't need to pass co2fact etc?
      common/work3/vmixarrs(ifull,kl,3),trsrc(ifull,kl),spare(ifull,kl)
      real updtr(ifull,kl),at(ifull,kl),ct(ifull,kl)
      trfact = grav * dt / dsig(1)
      co2fact=1000.*trfact*fair_molm/fc_molm
      o2fact=1000.*trfact*fair_molm/fo2_molm

c     rates of change of mixing ratio m from the surface are computed
c     using the expression
c       dm/dt = phi/(rho*dz)
c     where phi is a mass flux of air with same number of molecules as in
c     the tracer rho is the air density, dz is the height of the lowest box
c       dp/dz = -g*rho
c     so    rho*dz = - dp/g
c       p = ps*sigma, so dp = ps*dsigma
c     then  dm = - phi * g * dt / (ps * dsigma)
c          = - conflux * phi / ps

      if( iradon.ne.0 ) then
	 k2=max(2*iradon,1)/max(iradon,1)  ! 1 for iradon=0	 
        if(ntest.eq.1.and.mydiag)
     &        print *,'ktau,iradon,trfact ',ktau,iradon,trfact
        call radonsflux
        if(ntest.eq.1.and.mydiag)
     &       print *,'after radsflux tr1&2,trsrc ',
     &  tr(idjd,1,max(1,iradon)),tr(idjd,k2,max(1,iradon)),trsrc(idjd,1)
        call radonvmix (updtr, trfact )
        if(ntest.eq.1.and.mydiag)
     &       print *,'after radonvmix tr1&2,updtr ',
     &  tr(idjd,1,max(1,iradon)),tr(idjd,k2,max(1,iradon)),updtr(idjd,1)
        call trimcopy(at,ct,updtr,iradon) ! same as trim but holds bdy vals
        if(ntest.eq.1.and.mydiag)
     &       print *,'after trim tr1&2,updtr ',
     &  tr(idjd,1,max(1,iradon)),tr(idjd,k2,max(1,iradon)),updtr(idjd,1)
      end if

      if( ico2.ne.0 ) then
	 k2=max(2*ico2,1)/max(ico2,1)  ! 1 for ico2=0
        if(ntest.eq.1.and.mydiag)
     &        print *,'ktau,ico2,co2fact ',ktau,ico2,co2fact
        call co2sflux
        if(ntest.eq.1.and.mydiag)
     &       print *,'after co2sflux tr1&2,trsrc ',
     &      tr(idjd,1,max(1,ico2)),tr(idjd,k2,max(1,ico2)),trsrc(idjd,1)
        call co2vmix(updtr, co2fact )
        if(ntest.eq.1.and.mydiag)
     &       print *,'after co2vmix tr1&2,updtr ',
     &      tr(idjd,1,max(1,ico2)),tr(idjd,k2,max(1,ico2)),updtr(idjd,1)
        call trimcopy(at,ct,updtr,ico2)
        if(ntest.eq.1.and.mydiag)
     &       print *,'after trim tr1&2,updtr ',
     &      tr(idjd,1,max(1,ico2)),tr(idjd,k2,max(1,ico2)),updtr(idjd,1)
        if(ntest.eq.1)then
c         print *,'can mel sources ',ktau,trsrc(46,57,1),trsrc(39,52,1)
c         print *,'can mel conc lev1 ',ktau,tr(46,57,1,2),tr(39,52,1,2)
          call printa('co2 ',tr(:,:,k2),ktau,1,ia,ib,ja,jb,357.,1.)
          call printa('rado',tr(:,:,1),ktau,1,ia,ib,ja,jb,0.,1.)
        endif  ! (ntest.eq.1)
      endif

      if( io2.ne.0 ) then
        call o2vmix(updtr, o2fact )
        call trimcopy(at,ct,updtr,io2)
      end if

      if(iso2.gt.0) then
        so2fact(1) = trfact
        do i=2,nso2lev
          so2fact(i) = grav * dt / dsig(i)
        end do  !   i=1,nso2lev
        if( iso2.gt.0.and.nso2lev.lt.1 )stop 'vertmix: nso2lev.lt.1'
      endif

      if( iso2.ne.0 ) then
       call so2sflux
       call so2vmix(updtr )
       call trimcopy(at,ct,updtr,iso2)
	do k=1,klt
	 do iq=1,ilt*jlt
         tr(iq,k,max(1,iso2))=max(-100.,min(tr(iq,k,max(1,iso2)),9000.))
        enddo
       enddo
      endif

      if( iso4.ne.0 ) then
*       call so4sflux
*       call so4vmix( updtr )
*       call trimcopy(at,ct,updtr,iso4)
      end if

      if( ich4.ne.0 ) then
*       call ch4sflux
*       call ch4vmix( updtr )
*       call trimcopy(at,ct,updtr,ich4)
      end if

      if( io2.ne.0 ) then
*       call o2sflux
*       call o2vmix(updtr )
*       call trimcopy(at,ct,updtr,io2)
      end if
       return
      end
*............................................................
c     not needed, except for setting darlam boundary values?
      subroutine trimcopy(at,ct,updtr,itracer)
      include 'const_phys.h'
      include 'newmpar.h'
      include 'parm.h'
      include 'tracers.h'
      real updtr(ifull,kl),at(ifull,kl),ct(ifull,kl)
      call trim(at,ct,updtr,0)
      do k=1,klt
         do iq=1,ilt*jlt   !  alter all values for c-c
          tr(iq,k,itracer)=updtr(iq,k)
         enddo  ! iq loop
      enddo     !  k loop
      return
      end

