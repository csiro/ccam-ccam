      subroutine vertmix       ! globpea (for DARLAM redo equivs)
      include 'newmpar.h'
      parameter (ntest=0)
c     parameter (ipwr=1)       ! really can use (ipwr=min(1,nlocal))
c     parameter (ilnl=il**ipwr,jlnl=jl**ipwr)
      parameter (kcl_top=kl-2) ! maximum level for cloud top (conjob & vertmix)
      parameter (kscmom=0)     ! 0 default, 1 turns off Tiedtke momentum
      parameter (ndvmod=0)     ! 0 default, 1+ for dvmod tests
!     parameter (tied_rh=.75)  ! usually .75
!     parameter (tied_con=6.)  ! use 6., 10., or 25.
c     nlocal in parm.h         ! 0 local scheme, 1 or 2 non-local scheme
      include 'arrays.h'
      include 'dates.h'
      include 'indices.h'
      include 'kuocom.h'   ! also with kbsav,ktsav,convpsav,kscsea
      include 'map.h'      ! em, f, dpsldt, fu, fv, etc  not needed here?
      include 'nlin.h'
      include 'morepbl.h'
      include 'parm.h'
      include 'pbl.h'
      include 'sigs.h'
      include 'soil.h'
      include 'tracers.h'  ! ngas, nllp, ntrac
!     include 'xarrs.h'
      common/permsurf/ipsice,ipsea,ipland,iperm(ifull)
c     following common block passes uav, vav from sflux to vertmix
c     first 3 arrays also used in trim
      common/work3/qs(il,jl,kl),delthet(il,jl,kl),
     .             uav(il,jl,kl),vav(il,jl,kl)
!     N.B. uav & vav also used by pbldif; all of work3 used by tracervmix
      real cu(il,jl,kl)
      equivalence (cu,vav)
      common/work3b/rkh_nl(il,jl,kl),rkm_nl(il,jl,kl)
      common/work3cc/gamat(il,jl,kl)      ! gamat only used in vertmix, pbldif
      common/work3d/gamaq(il,jl,kl)      ! gamaq only used in vertmix, pbldif
      common/work3e/rhs(il,jl,kl),dzk(il,jl,kl) ! rhs also used in pbldif
      common/work3f/e1(ijk),e2(ijk),e4(ijk)
      real gamau(il,jl,kl),gamav(il,jl,kl) ! just for pbldif
      common/work2/csq(il,jl),delq(il,jl),dvmod(il,jl)
     . ,dz(il,jl),dzr(il,jl),fm(il,jl),fh(il,jl),ri(il,jl)
     . ,sqmxl(il,jl),x(il,jl),zhv(il,jl),the1(ifull),dum2(il,jl,6)
      real sighkap(kl),delons(kl),delh(kl),prcpv(kl)
      real sigkap(kl)
      real at(il,jl,kl),au(il,jl,kl),ct(il,jl,kl)
      real zh(il,jl,kl),tmnht(il,jl,kl)
      real gt(il,jl,kl),guv(il,jl,kl)
      real rkm(il,jl,kl),rkh(il,jl,kl)
      equivalence (gt,rkh_nl),(guv,rkm_nl)
      equivalence (rkh,e1),(rkm,e2)
      equivalence (tmnht,at,un),(zh,au,e4)
      equivalence (gamat,ct)
      data r/287./,hl/2.5104e6/,cp/1.00464e3/
c     set coefficients for Louis scheme
      data bprm/4.7/,cm/7.4/,ch/5.3/,amxlsq/100./,vkar3/.35/,vkar4/.4/
      data bprmj/5./,cmj/5./,chj/2.6/
      common /es_table/ table(0:220)
c Arithmetic statement functions to replace call to establ.
c T is temp in Kelvin, which should lie between 123.16 and 343.16;
c TDIFF is difference between T and 123.16, subject to 0 <= TDIFF <= 220
      tdiff(tm)=min(max(tm-123.16 , 0.) , 220.)
      establ(tm) =(1.-(tdiff(tm)-aint(tdiff(tm))))*table(int(tdiff(tm)))
     &           + (tdiff(tm)-aint(tdiff(tm)))*table(int(tdiff(tm))+1)

      hlcp=hl/cp
      roncp=r/cp
      rong=r/9.806
      do k=1,kl-1
       sighkap(k)=sigmh(k+1)**(-roncp)
       delons(k)=rong *((sig(k+1)-sig(k))/sigmh(k+1))
      enddo      ! k loop
      do k=1,kl
       delh(k)=-rong *dsig(k)/sig(k)  ! sign of delh defined so always +ve
       sigkap(k)=sig(k)**(-roncp)
      enddo      ! k loop
      if(diag)then
        print *,'sig ',sig
        print *,'dsig ',dsig
        print *,'delh ',delh
        print *,'vertmix uin ',(u(id,jd,k),k=1,kl)
        print *,'vertmix vin ',(v(id,jd,k),k=1,kl)
      endif
      rlogs1=log(sig(1))
      rlogs2=log(sig(2))
      rlogh1=log(sigmh(2))
      rlog12=1./(rlogs1-rlogs2)
      do iq=1,ifull
       tmnht(iq,1,1)=(t(iq,1,2)*rlogs1-t(iq,1,1)*rlogs2+
     .               (t(iq,1,1)-t(iq,1,2))*rlogh1)*rlog12
c      n.b. an approximate zh is quite adequate for this routine
       zh(iq,1,1)=t(iq,1,1)*delh(1)
      enddo      ! iq loop
      do k=2,kl-1
       do iq=1,ifull
        zh(iq,1,k)=zh(iq,1,k-1)+t(iq,1,k)*delh(k)
        tmnht(iq,1,k) =(t(iq,1,k)+t(iq,1,k+1))*.5
       enddo     ! iq loop
      enddo      !  k loop
      do kq=1,ijk-ifull  ! top level set separately to 0 for gt & guv
       rkh(kq,1,1)=0.
       rkm(kq,1,1)=0.
      enddo
      do k=1,kl
       do iq=1,ifull
        rhs(iq,1,k)=t(iq,1,k)*sigkap(k)  ! rhs is theta here
        es=establ(t(iq,1,k))
        qs(iq,1,k)=.622*es/(ps(iq,1)*sig(k)-es)
       enddo     ! iq loop
      enddo      !  k loop

c     ************ section for tiedtke shallow convection *******************
!     effective difference from bmrc's is mixing of momentum
!     and range of levels
      if(ksc.gt.0)then
        do k=kscbase,ksc
         prcpv(k)=sig(k)**(-r/cp)
        enddo    ! k loop
        do iq=1,ifull
         the1(iq)=prcpv(kscbase)*t(iq,1,kscbase)*
     .                   (t(iq,1,kscbase) + .5*hlcp*qs(iq,1,kscbase))
     .                  /(t(iq,1,kscbase) - .5*hlcp*qs(iq,1,kscbase))
        enddo    ! iq loop

        if(kscsea.eq.1)then   ! Tiedtke done only over sea
          do k=kscbase+1,ksc
           do ip=ipland+1,ipsea                                      ! sea only
            iq=iperm(ip)                                             ! sea only
            the=prcpv(k)*t(iq,1,k)*(t(iq,1,k) + .5*hlcp*qs(iq,1,k))
     .                            /(t(iq,1,k) - .5*hlcp*qs(iq,1,k))
            if(qg(iq,1,kscbase).gt.tied_rh*qs(iq,1,kscbase).
     .       and.the.lt.the1(iq))then               !  default tied_rh=.75
              if(kscmom.eq.0)rkm(iq,1,k-1)=tied_con !  m**2/sec  6., originally 10.
              rkh(iq,1,k-1)=tied_con                !  m**2/sec  6., originally 10.
              if(kscmom.eq.0)rkm(iq,1,k)=2.         !  m**2/sec
              rkh(iq,1,k)=2.                        !  m**2/sec
            endif ! (qg(iq,1,kscbase).gt.rhscon*qs(iq,1,kscbase).....
           enddo  ! iq loop
          enddo   ! end of k=kscbase+1,ksc loop
        else      !  i.e. original scheme over land and sea
          do k=kscbase+1,ksc
           do iq=1,ifull
            the=prcpv(k)*t(iq,1,k)*(t(iq,1,k) + .5*hlcp*qs(iq,1,k))
     .                            /(t(iq,1,k) - .5*hlcp*qs(iq,1,k))
            if(qg(iq,1,kscbase).gt.tied_rh*qs(iq,1,kscbase).
     .       and.the.lt.the1(iq))then               !  default tied_rh=.75
              if(kscmom.eq.0)rkm(iq,1,k-1)=tied_con !  m**2/sec  6., originally 10.
              rkh(iq,1,k-1)=tied_con                !  m**2/sec  6., originally 10.
              if(kscmom.eq.0)rkm(iq,1,k)=2.         !  m**2/sec
              rkh(iq,1,k)=2.                        !  m**2/sec
            endif ! (qg(iq,1,kscbase).gt.rhscon*qs(iq,1,kscbase).....
           enddo  ! iq loop
          enddo   ! end of k=kscbase+1,ksc loop
        endif     ! (kscsea.eq.1)  .. else ..
      endif       ! (ksc.gt.0)
c     *********** end of tiedtke shallow convection section *****************

      do k=1,kl-1
       do iq=1,ifull
        delthet(iq,1,k)=rhs(iq,1,k+1)-rhs(iq,1,k)  ! rhs is theta here
       enddo     ! iq loop
      enddo      !  k loop

c     ************ section for geleyn shallow convection *******************
      if(ksc.lt.0)then
!       do k=-ksc,kl/2
        do k=kscbase,-ksc    ! new usage of ksc Thu  02-17-2000
         do iq=1,ifull
          delthet(iq,1,k)=delthet(iq,1,k)-hlcp*max(0.,
     .                 qs(iq,1,k+1)-qg(iq,1,k+1)-qs(iq,1,k)+qg(iq,1,k) )
         enddo  ! iq loop
        enddo   !  k loop
      endif     ! (ksc.lt.0)
c     *********** end of geleyn shallow convection section *****************

      if(diag.or.ntest.eq.1)then
        Print *,'after shallow convection'
        print *,'rkh ',(rkh(id,jd,k),k=1,kl)
        print *,'rkm ',(rkm(id,jd,k),k=1,kl)
        print *,'t ',(t(id,jd,k),k=1,kl)
        print *,'qg ',(qg(id,jd,k),k=1,kl)
      endif

      do k=1,kl-1
       do iq=1,ifull
        dz(iq,1) =-tmnht(iq,1,k)*delons(k)  ! this is z(k+1)-z(k)
        dzr(iq,1)=1./dz(iq,1)
        zhv(iq,1)=1./zh(iq,1,k)
        if(ndvmod.eq.0)then
          dvmod(iq,1)=sqrt( (uav(iq,1,k+1)-uav(iq,1,k))**2
     .                     +(vav(iq,1,k+1)-vav(iq,1,k))**2 )
        else
          dvmod(iq,1)=ndvmod  ! just for tests
        endif    ! (ndvmod.eq.0)

c       x is bulk ri *(dvmod **2)
        x(iq,1)=9.806*dz(iq,1)*(delthet(iq,1,k)/
     .   (tmnht(iq,1,k)*sighkap(k))+.61*(qg(iq,1,k+1)-qg(iq,1,k)))

c       fm and fh denote f(louis style)*dvmod
c       nb. an error exists in the Louis expression for c; this is corrected
c       in the current code
        csq(iq,1) = zhv(iq,1)*
     .              (((1.+dz(iq,1)*zhv(iq,1))**(1./3.)-1.)*dzr(iq,1))**3
       enddo     ! iq loop

       if(nlocal.gt.0)then
         do iq=1,ifull
          dzk(iq,1,k)=-tmnht(iq,1,k)*delons(k) ! this is z(k+1)-z(k)
         enddo   ! iq loop
       endif     ! (nlocal.gt.0)

       if(nvmix.lt.0)then   ! use quasi-neutral mixing for test runs only
         do iq=1,ifull
          fm(iq,1)=1.
          fh(iq,1)=1.
         enddo   ! iq loop
         if(k.lt.kl/2)then
           do iq=1,ifull
            if(land(iq,1))fh(iq,1)=abs(nvmix)
           enddo   ! iq loop
         endif
       endif     !  (nvmix.lt.0)

       if(nvmix.eq.2)then
c        code as used in early mesoscale version
         ric=0.067 *(dz(1+il/2,1+jl/2)*100.)**0.25
         do iq=1,ifull
          sqmxl(iq,1)=(vkar3*zh(iq,1,k)/(1.+vkar3*zh(iq,1,k)/amxlsq))**2
          dvmod(iq,1)=max(dvmod(iq,1),1.e-4*dz(iq,1))
          ri(iq,1)=x(iq,1)/dvmod(iq,1)**2
          if(ri(iq,1).lt. 0.)then  ! unstable case
c           first do momentum
            denma=dvmod(iq,1)+cm*( 2.*bprm*sqmxl(iq,1)*
     .                            sqrt(-x(iq,1)*csq(iq,1)) )
            fm(iq,1)=dvmod(iq,1)-(2.*bprmj *x(iq,1))/denma
c           now heat
            denha=dvmod(iq,1)+ch*( 2.*bprm*sqmxl(iq,1)*
     .                            sqrt(-x(iq,1)*csq(iq,1)) )
            fh(iq,1)=( dvmod(iq,1)-(2.*bprm *x(iq,1))/denha  )/.74
          else                     ! stable case
c           linearly interpolate to calculate f for stable case; f lies between
c           1 (ri = 0) and 0 (r = ric);  multiply by dvmod to give fm and fh
            fm(iq,1)=max( dvmod(iq,1)*(ric-ri(iq,1))/ric , 0. )
            fh(iq,1)=fm(iq,1)/.74    !vertmixb didn't have the .74
          endif
         enddo   ! iq loop
       endif

       if(nvmix.eq.3)then
c        newish code, vkar4 & also with ric=.5 typically
         ric=.5
         do iq=1,ifull
          sqmxl(iq,1)=(vkar4*zh(iq,1,k)/(1.+vkar4*zh(iq,1,k)/amxlsq))**2
          dvmod(iq,1)=max( dvmod(iq,1) , 1. )
          ri(iq,1)=x(iq,1)/dvmod(iq,1)**2
          if(ri(iq,1).lt. 0.)then  ! unstable case
c           first do momentum
            denma=dvmod(iq,1)+cmj*( 2.*bprmj*sqmxl(iq,1)*
     .                            sqrt(-x(iq,1)*csq(iq,1)) )
            fm(iq,1)=dvmod(iq,1)-(2.*bprmj *x(iq,1))/denma
c           now heat
            denha=dvmod(iq,1)+chj*( 2.*bprmj*sqmxl(iq,1)*
     .                            sqrt(-x(iq,1)*csq(iq,1)) )
            fh(iq,1)=dvmod(iq,1)-(2.*bprmj *x(iq,1))/denha
          else                     ! stable case
c           linearly interpolate to calculate f for stable case; f lies between
c           1 (ri = 0) and 0 (r = ric);  multiply by dvmod to give fm and fh
            fm(iq,1)=max( dvmod(iq,1)*(ric-ri(iq,1))/ric , 0. )
            fh(iq,1)=fm(iq,1)
          endif
         enddo   ! iq loop
       endif

       if(nvmix.eq.4)then
c        newest code, stable same as csiro9 here
         do iq=1,ifull
          sqmxl(iq,1)=(vkar4*zh(iq,1,k)/(1.+vkar4*zh(iq,1,k)/amxlsq))**2
          dvmod(iq,1)=max( dvmod(iq,1) , 1. )
          ri(iq,1)=x(iq,1)/dvmod(iq,1)**2
          if(ri(iq,1).lt. 0.)then  ! unstable case
c           first do momentum
            denma=dvmod(iq,1)+cmj*( 2.*bprmj*sqmxl(iq,1)*
     .                            sqrt(-x(iq,1)*csq(iq,1)) )
            fm(iq,1)=dvmod(iq,1)-(2.*bprmj *x(iq,1))/denma
c           now heat
            denha=dvmod(iq,1)+chj*( 2.*bprmj*sqmxl(iq,1)*
     .                            sqrt(-x(iq,1)*csq(iq,1)) )
            fh(iq,1)=dvmod(iq,1)-(2.*bprmj *x(iq,1))/denha
          else                     ! stable case
c           The following is the original Louis stable formula
            fm(iq,1)=dvmod(iq,1)/(1.+4.7*ri(iq,1))**2
            fh(iq,1)=fm(iq,1)
          endif
         enddo   ! iq loop
       endif     ! (nvmix.eq.4)

c      calculate k's, guv and gt
c      nonlocal scheme usually applied to temperature and moisture, NOT momentum
c      (i.e. local scheme is applied to momentum for nlocal=0,1)
       if(nvmix.ne.0)then    ! use nvmix=0 only for special test runs
         do iq=1,ifull
          rkm(iq,1,k)=rkm(iq,1,k)+fm(iq,1)*sqmxl(iq,1)*dzr(iq,1)
          rkh(iq,1,k)=rkh(iq,1,k)+fh(iq,1)*sqmxl(iq,1)*dzr(iq,1)
         enddo   ! iq loop
       endif     ! (nvmix.ne.0)

       if(diag.or.ntest.eq.1)then
         print *,'k,dt,delsig,sqmxl,dzr ',
     .            k,dt,delsig,sqmxl(id,jd),dzr(id,jd)
         print *,'k,t,t+,ps ',k,t(id,jd,k),t(id,jd,k+1),ps(id,jd)
         print *,'k,qg,qg+ ',k,qg(id,jd,k),qg(id,jd,k+1)
         print *,'k,qs,qs+ ',k,qs(id,jd,k),qs(id,jd,k+1)
         es=establ(t(id,jd,k))
         esp=establ(t(id,jd,k+1))
         print *,'k,es,es+,delthet ',k,es,esp,delthet(id,jd,k)
         print *,'k,fm,dvmod,ri,csq ',
     .            k,fm(id,jd),dvmod(id,jd),ri(id,jd),csq(id,jd)
         print *,'x,zh,tmnht,dz,delh,sighkap ',x(id,jd),zh(id,jd,k),
     .                 tmnht(id,jd,k),dz(id,jd),delh(k),sighkap(k)
       endif
      enddo      ! end of k loop

      if(diag.or.ntest.eq.1)then
        print *,'before possible non-local calls'
        print *,'rkh ',(rkh(id,jd,k),k=1,kl)
        print *,'rkm ',(rkm(id,jd,k),k=1,kl)
      endif

      if(nlocal.ne.0)then
        call pbldif(rhs,rkm,rkh,rkm_nl,rkh_nl,gamat,gamaq,
     .              gamau,gamav)  ! rhs is theta here
!       NCAR info is returned in rkm_nl and rkh_nl arrays
!       for nlocal -ve, use NCAR scheme without non-local terms
!       nlocal = -2   NCAR    for Tq, NCAR    for uv
!                -1   NCAR    for Tq,  CAR    for uv
!                 0    CAR    for Tq,  CAR    for uv
!                 1   NCAR_NL for Tq,  CAR    for uv
!                 2   NCAR_NL for Tq, NCAR_NL for uv
!       transfer to rkm & rkh arrays (possibly adding shallow convection)
        do kq=1,ijk-ifull                     !  NCAR    for Tq
         rkh(kq,1,1)=rkh_nl(kq,1,1)              !+rkh_nl(kq,1,1)
        enddo    ! kq loop
        if(nlocal.eq.2.or.nlocal.eq.-2)then   !  NCAR    for uv
          do kq=1,ijk-ifull
           rkm(kq,1,1)=rkm_nl(kq,1,1)         !  +rkm_nl(kq,1,1)
          enddo  ! kq loop
        endif
        if(diag.or.ntest.eq.1)then
          print *,'rkh ',(rkh(id,jd,k),k=1,kl)
          print *,'rkh_nl ',(rkh(id,jd,k),k=1,kl)
          print *,'gamat ',(gamat(id,jd,k),k=1,kl)
          print *,'gamaq ',(gamaq(id,jd,k),k=1,kl)
          print *,'rkm ',(rkm(id,jd,k),k=1,kl)
          print *,'rkm_nl ',(rkm(id,jd,k),k=1,kl)
          print *,'t ',(t(id,jd,k),k=1,kl)
          print *,'qg ',(qg(id,jd,k),k=1,kl)
        endif
      endif      !  (nlocal.gt.0)
!!    rkh_nl & rkm_nl arrays available from here on

      do k=1,kl-1
        delsig=sig(k+1)-sig(k)
        do iq=1,ifull
          dz(iq,1) =-tmnht(iq,1,k)*delons(k)  ! this is z(k+1)-z(k)
          dzr(iq,1)=1./dz(iq,1)
          guv(iq,1,k)= rkm(iq,1,k)*dt *delsig *dzr(iq,1)**2
          gt(iq,1,k)=  rkh(iq,1,k)*dt *delsig *dzr(iq,1)**2
         enddo   ! iq loop
      enddo      ! k loop

      do iq=1,ifull
       guv(iq,1,kl)=0.
       gt(iq,1,kl)=0.
      enddo      ! iq loop
      if(diag.or.ntest.eq.1)then
        print *,'vertmix guv ',(guv(id,jd,k),k=1,kl)
        print *,'vertmix gt ',(gt(id,jd,k),k=1,kl)
        print *,'in vertmix t ',(t(id,jd,k),k=1,kl)
        print *,'vertmix theta before trim ',(rhs(id,jd,k),k=1,kl)
        print *,'tn before trim ',(tn(id,jd,k),k=1,kl)
      endif
      if(diag)then
        call maxmin(guv,'gu',ktau,1.e3,kl)
        call printa('t   ',t(1,1,nlv),ktau,nlv,ia,ib,ja,jb,200.,1.)
        call printa('tss ',tss(1,1),ktau,nlv,ia,ib,ja,jb,200.,1.)
        call printa('fg  ',fg(1,1),ktau,nlv,ia,ib,ja,jb,0.,1.)
        call printa('thet',rhs(1,1,nlv),ktau,nlv,ia,ib,ja,jb,200.,1.)
        call printa('tn  ',tn(1,1,nlv),ktau,nlv,ia,ib,ja,jb,0.,100.*dt)
      endif

      if(ncvmix.gt.0)then  ! cumulus mixing of momentum - jlm version
        do k=kuocb+1-ncvmix,kcl_top-1
         do iq=1,ifull
!          for no-conv-points, doing k=1-ncvmix,1 with convpsav=0.
!          1 below for ncvmix=2
           if(kbsav(iq).gt.0.and.k.ge.kbsav(iq)+1-ncvmix
     .                      .and.k.lt.ktsav(iq))
     .       guv(iq,1,k)=guv(iq,1,k)-convpsav(iq)*.5      ! with factor of .5
          enddo  ! iq loop
        enddo    ! k loop
        if(diag)then
          print *,'vertmix after conv; kb,kt,convp'
     .     ,kbsav(idjd),ktsav(idjd),convpsav(idjd)
          print *,'new guv',(guv(id,jd,k),k=1,kl)
        endif
      endif      !   (ncvmix.gt.0)

      conflux=9.806*dt/dsig(1)
      condrag=9.806*dt/(dsig(1)*r)
c     first do theta (then convert back to t)
      if(ntsur.ne.6)then   !  Tue  09-28-1993
        do iq=1,ifull
         at(iq,1,1)=0.
         rhs(iq,1,1)=rhs(iq,1,1)-(conflux/cp)*fg(iq,1)/ps(iq,1)
        enddo    ! iq loop
      else       ! for ntsur=6 using cdtq
        do iq=1,ifull       ! N.B. cdtq array not set up in sflux nowadays
         cdt=fg(iq,1)*r*tss(iq,1)/((tss(iq,1)-rhs(iq,1,1))*cp*ps(iq,1))
         at(iq,1,1) =cdt*condrag/tss(iq,1)
         rhs(iq,1,1)=rhs(iq,1,1)-at(iq,1,1)*tss(iq,1)
        enddo    ! iq loop
      endif
      if(nlocal.gt.0) then
        do k=2,kl
         do iq=1,ifull
c         kcn: et here : counter gradient for T due to non-local scheme
          rhs(iq,1,k) = rhs(iq,1,k)+
     .               ( gt(iq,1,k-1)*dzk(iq,1,k-1)*gamat(iq,1,k-1)
     .                -gt(iq,1,k)*dzk(iq,1,k)*gamat(iq,1,k) )/dsig(k)
         enddo   ! iq loop
        enddo    !  k loop
c       kcn: add in counter gradient at k=1
        do iq=1,ifull
          rhs(iq,1,1) = rhs(iq,1,1)-
     .                  gt(iq,1,1)*dzk(iq,1,1)*gamat(iq,1,1)/dsig(1)
         enddo   ! iq loop
      endif      !  (nlocal.gt.0)

      if(npanels.gt.0)then
        do k=2,kl
         do iq=1,ifull
          at(iq,1,k) =-gt(iq,1,k-1)/dsig(k)
         enddo   ! iq loop
        enddo    !  k loop
        do k=1,kl
*        vdir nodep
         do iq=1,ifull
          ct(iq,1,k) =-gt(iq,1,k)/dsig(k)
         enddo   ! iq loop
        enddo    !  k loop
      else       ! i.e. npanels=0   DARLAM
        do k=2,kl
         do iq=1,ifull
          at(iq,1,k) =-gt(iq,1,k-1)/dsig(k)
         enddo   ! iq loop
        enddo    !  k loop
        do k=1,kl
*        vdir nodep
         do iq=1,ifull
          ct(iq,1,k) =-gt(iq,1,k)/dsig(k)
         enddo   ! iq loop
        enddo    !  k loop
      endif      !  (npanels.gt.0) .. else ..
      if(diag.or.ntest.eq.2)then
        print *,'ktau,fg,tss,ps ',ktau,fg(id,jd),tss(id,jd),ps(id,jd)
        print *,'at ',(at(id,jd,k),k=1,kl)
        print *,'ct ',(ct(id,jd,k),k=1,kl)
        print *,'rhs ',(rhs(id,jd,k),k=1,kl)
        print *,'before trim t1,t2,tn1,tn2 ',
     .                     t(id,jd,1),t(id,jd,2),tn(id,jd,1),tn(id,jd,2)
      endif      ! (ntest.eq.2)

      call trim(at,ct,rhs,0)   ! for t
      do k=1,kl
       if(npanels.gt.0)then
        if(nvsplit.eq.-1)then  ! splitting u,v,T in vertmix
          do iq=1,ifull
           t(iq,1,k)=rhs(iq,1,k)/sigkap(k)
          enddo  ! iq loop
        else
          do iq=1,ifull
           tn(iq,1,k)=tn(iq,1,k)+(rhs(iq,1,k)/sigkap(k)-t(iq,1,k))/dt ! tendency
          enddo  ! iq loop
        endif    ! (nvsplit.eq.-1)
       else         ! npanels=0  DARLAM
         do j=2,jl
          do i=2,il
           t(i,j,k)=rhs(i,j,k)/sigkap(k)
          enddo
         enddo
       endif     ! (npanels.gt.0)
      enddo      !  k loop
      if(ntest.eq.2)print *,'after trim t1,t2,tn1,tn2,rhs1 ',
     .        t(id,jd,1),t(id,jd,2),tn(id,jd,1),tn(id,jd,2),rhs(id,jd,1)
      if(diag)then
        print *,'vertmix eg,fg,cdtq,land '
     .                  ,eg(id,jd),fg(id,jd),cdtq(id,jd),land(id,jd)
        print *,'vertmix theta after trim ',(rhs(id,jd,k),k=1,kl)
        call printa('thet',rhs(1,1,nlv),ktau,nlv,ia,ib,ja,jb,200.,1.)
        print *,'tn after trim ',(tn(id,jd,k),k=1,kl)
c       call printa('tn  ',tn(1,1,nlv),ktau,nlv,ia,ib,ja,jb,0.,100.*dt)
c       print *,'following maxmin from within vertmix'
c       call maxmin(rhs,'th',ktau,1.,kl)
c       call maxmin(t,'tt',ktau,1.,kl)
      endif

c     now do moisture
      do kq=1,ijk
       rhs(kq,1,1)=qg(kq,1,1)
      enddo      ! kq loop
      if(ntsur.ne.6)then   !  Tue  09-28-1993
        do iq=1,ifull
         rhs(iq,1,1)=rhs(iq,1,1)-(conflux/hl)*eg(iq,1)/ps(iq,1)
        enddo    ! iq loop
      else       ! for ntsur=6 using cdtq
!       do iq=1,ifull             ! DARLAM has this possibility for ntsur=6
!        at(iq,1,1)=dw(iq,1)*at(iq,1,1)
!       enddo
        do iq=1,ifull
         es=establ(tss(iq,1))
         rsts=.622*es/(ps(iq,1)-es)
         cdq=eg(iq,1)*r*tss(iq,1)/((rsts-qg(iq,1,1))*2.5e6*ps(iq,1))
         at(iq,1,1) =cdq*condrag/rsts
         rhs(iq,1,1)=rhs(iq,1,1)-at(iq,1,1)*rsts
         enddo   ! iq loop
      endif      ! (ntsur.ne.6)
c     could add extra sfce moisture flux term for crank-nicholson
      if(nlocal.gt.0) then
        if(diag.or.ntest.eq.1)then
          print *,'before nlocal, rhs for qg ',(rhs(id,jd,k),k=1,kl)
          print *,'gt ',(gt(id,jd,k),k=1,kl)
          print *,'dzk ',(dzk(id,jd,k),k=1,kl)
          print *,'gt++ ',(gt(id+1,jd+1,k),k=1,kl)
          print *,'dzk++ ',(dzk(id+1,jd+1,k),k=1,kl)
        endif
        do k=2,kl
         do iq=1,ifull
c         kcn: et here : counter gradient for qg due to non-local scheme
          rhs(iq,1,k) = rhs(iq,1,k)+
     .               ( gt(iq,1,k-1)*dzk(iq,1,k-1)*gamaq(iq,1,k-1)
     .                -gt(iq,1,k)*dzk(iq,1,k)*gamaq(iq,1,k) )/dsig(k)
         enddo  ! iq loop
        enddo   !  k loop
c       kcn: add in counter gradient at k=1
        do iq=1,ifull
          rhs(iq,1,1) = rhs(iq,1,1)-
     .                  gt(iq,1,1)*dzk(iq,1,1)*gamaq(iq,1,1)/dsig(1)
        enddo    ! iq loop
        if(diag.or.ntest.eq.1)then
          print *,'gamaq ',(gamaq(id,jd,k),k=1,kl)
          print *,'after nlocal, rhs for qg ',(rhs(id,jd,k),k=1,kl)
        endif
      endif   !  (nlocal.gt.0)
      call trim(at,ct,rhs,0)    ! for qg
      if(npanels.gt.0)then
        do kq=1,ijk
         qg(kq,1,1)=rhs(kq,1,1)
        enddo    ! kq loop
      else
        do k=1,kl
         do j=2,jl
          do i=2,il
           qg(i,j,k)=rhs(i,j,k)
          enddo
         enddo
        enddo    ! k loop
      endif      !  (npanels.gt.0)
      if(diag.or.ntest.eq.1)then
        print *,'after trim, qg ',(qg(id,jd,k),k=1,kl)
      endif

c     now do trace gases
      if(ngas.gt.0)call TracerVMix( at, ct, rhs )
      if(diag.or.ntest.eq.1)then
        print *,'vertmix after trace gases'
      endif

!     from here on just u and v stuff
      if(npanels.gt.0)then
        do iq=1,ifull
         au(iq,1,1)=cduv(iq,1)*condrag/tss(iq,1)   ! globpea
        enddo    ! iq loop
        do k=2,kl
         do iq=1,ifull
          au(iq,1,k) =-guv(iq,1,k-1)/dsig(k)  ! globpea
         enddo   ! iq loop
        enddo    !  k loop
        do k=1,kl
*        vdir nodep
         do iq=1,ifull
          cu(iq,1,k) =-guv(iq,1,k)/dsig(k)    ! globpea
         enddo   ! iq loop
        enddo    !  k loop
      else       ! i.e. npanels=0   DARLAM
        do iq=1,ifull-1   !xxx   following for u
         au(iq,1,1) =.5*(cduv(iq,1)+cduv(iq+1,1))*condrag/tss(iq,1)
        enddo    ! iq loop
        do k=2,kl
         do iq=1,ifull-1   !xxx   following for u
          au(iq,1,k) =-.5*(guv(iq,1,k-1)+guv(iq+1,1,k-1))/dsig(k) ! DARLAM
         enddo   ! iq loop
        enddo    !  k loop
        do k=1,kl
*        vdir nodep
         do iq=1,ifull-1   !      following for u
          cu(iq,1,k) =-.5*(guv(iq,1,k)+guv(iq+1,1,k))/dsig(k)    ! DARLAM
         enddo   ! iq loop
        enddo    !  k loop
      endif      !  (npanels.gt.0) .. else ..
      if(diag.or.ntest.eq.2)then
        print *,'au ',(au(id,jd,k),k=1,kl)
        print *,'cu ',(cu(id,jd,k),k=1,kl)
      endif      ! (ntest.eq.2)

c     first do u
      do kq=1,ijk
       rhs(kq,1,1)=u(kq,1,1)
      enddo      ! kq loop
      call trim(au,cu,rhs,0)
      if(npanels.gt.0)then
        if(nvsplit.eq.-1)then  ! splitting u,v,T in vertmix
          do kq=1,ijk
           u(kq,1,1)=rhs(kq,1,1)
           un(kq,1,1)=0.
          enddo  ! kq loop
        else
          do kq=1,ijk
           un(kq,1,1)=(rhs(kq,1,1)-u(kq,1,1))/dt  ! to give tendency
          enddo  ! kq loop
        endif    ! (nvsplit.eq.-1)
      else       ! npanels=0 for DARLAM
c       DARLAM usually does split for u & v
        do k=1,kl
         do j=3,jl-1
          do i=2,il-1
           u(i,j,k)=rhs(i,j,k)
          enddo
         enddo
        enddo   ! k loop
        do iq=1,ifull-il    !   following for DARLAM v
         au(iq,1,1) =.5*(cduv(iq,1)+cduv(iq,1+1))*condrag/tss(iq,1)
        enddo   ! iq loop
        do k=2,kl
         do iq=1,ifull-il   !   following for DARLAM v
          au(iq,1,k) =-.5*(guv(iq,1,k-1)+guv(iq,1+1,k-1))/dsig(k)
         enddo  ! iq loop
        enddo   !  k loop
        do k=1,kl
         do iq=1,ifull-il   !   following for DARLAM v
          cu(iq,1,k) =-.5*(guv(iq,1,k)+guv(iq,1+1,k))/dsig(k)
         enddo  ! iq loop
        enddo   !  k loop
      endif     ! (npanels.gt.0)
      if(diag)then
        print *,'after vertmix qg ',(qg(id,jd,k),k=1,kl)
        print *,'vertmix au ',(au(id,jd,k),k=1,kl)
        print *,'cduv,cduv+1,tss ',cduv(id,jd),cduv(id+1,jd),tss(id,jd)
        print *,'vertmix cu ',(cu(id,jd,k),k=1,kl)
        print *,'after vertmix u ',(u(id,jd,k),k=1,kl)
      endif

c     now do v; with properly unstaggered au,cu
      do kq=1,ijk
       rhs(kq,1,1)=v(kq,1,1)
      enddo  ! kq loop
      call trim(au,cu,rhs,0)    ! note now that au, cu unstaggered globpea
      if(npanels.gt.0)then
        if(nvsplit.eq.-1)then  ! splitting u,v,T in vertmix
          do kq=1,ijk
           v(kq,1,1)=rhs(kq,1,1)
           vn(kq,1,1)=0.
          enddo  ! kq loop
        else
          do kq=1,ijk
           vn(kq,1,1)=(rhs(kq,1,1)-v(kq,1,1))/dt  ! to give tendency
          enddo  ! kq loop
        endif    ! (nvsplit.eq.-1)
      else  ! for DARLAM
        do k=1,kl
         do j=2,jl-1
          do i=3,il-1
           v(i,j,k)=rhs(i,j,k)
          enddo
         enddo
        enddo   !  k loop
      endif  ! (npanels.gt.0)
      if(diag)then
        print *,'vertmix au for v ',(au(id,jd,k),k=1,kl)
        print *,'cduv,cduvj+1,tss ',cduv(id,jd),cduv(id,jd+1),tss(id,jd)
        print *,'vertmix cu ',(cu(id,jd,k),k=1,kl)
        print *,'vertmix v ',(v(id,jd,k),k=1,kl)
        call maxmin(u,'uv',ktau,1.,kl)
        call maxmin(v,'vv',ktau,1.,kl)
      endif
      return
      end

      Subroutine TracerVMix( at, ct, updTr )
c     this routine does the vertical mixing of tracers
      parameter (ntest=0)  ! 1 for diag prints
      include 'newmpar.h'
      include 'constant.h'
      include 'dates.h'         ! dt
      include 'nsibd.h'
      include 'parm.h'
      include 'sigs.h'          ! dsig
      include 'soil.h'
      include 'tracers.h'       ! jlm don't need to pass co2fact etc?
      common/work3/vmixarrs(il,jl,kl,3),trsrc(il,jl,kl)
      real updTr(il,jl,kl),at(il,jl,kl),ct(il,jl,kl)
      trfact = g * dt / dsig(1)
      CO2fact=1000.*trfact*fAIR_MolM/fC_MolM
      O2fact=1000.*trfact*fAIR_MolM/fO2_MolM

c     Rates of change of mixing ratio m from the surface are computed
C     using the expression
c       dm/dt = phi/(rho*dz)
c     where phi is a mass flux of AIR with SAME number of molecules as in
c     the tracer rho is the air density, dz is the height of the lowest box
c       dp/dz = -g*rho
c     so    rho*dz = - dp/g
c       p = ps*sigma, so dp = ps*dsigma
c     then  dm = - phi * g * dt / (ps * dsigma)
c          = - conflux * phi / ps

      if( iRADON.NE.0 ) then
        if(ntest.eq.1)print *,'ktau,iRADON,trfact ',ktau,iRADON,trfact
        call RadonSFlux
        if(ntest.eq.1)print *,'after radsflux tr1&2,trsrc ',
     .         tr(idjd,1,1,iradon),tr(idjd,1,2,iradon),trsrc(idjd,1,1)
        call RadonVMix (updTr, trfact )
        if(ntest.eq.1)print *,'after radonvmix tr1&2,updtr ',
     .          tr(idjd,1,1,iradon),tr(idjd,1,2,iradon),updtr(idjd,1,1)
        call trimcopy(at,ct,updTr,iRADON) ! same as trim but holds bdy vals
        if(ntest.eq.1)print *,'after trim tr1&2,updtr ',
     .          tr(idjd,1,1,iradon),tr(idjd,1,2,iradon),updtr(idjd,1,1)
      end if

      if( iCO2.NE.0 ) then
        if(ntest.eq.1)print *,'ktau,iCO2,CO2fact ',ktau,iCO2,CO2fact
        call CO2SFlux
        if(ntest.eq.1)print *,'after co2sflux tr1&2,trsrc ',
     .             tr(idjd,1,1,ico2),tr(idjd,1,2,ico2),trsrc(idjd,1,1)
        call CO2vmix(updTr, CO2fact )
        if(ntest.eq.1)print *,'after co2vmix tr1&2,updtr ',
     .             tr(idjd,1,1,ico2),tr(idjd,1,2,ico2),updtr(idjd,1,1)
        call trimcopy(at,ct,updTr,iCO2)
        if(ntest.eq.1)print *,'after trim tr1&2,updtr ',
     .             tr(idjd,1,1,ico2),tr(idjd,1,2,ico2),updtr(idjd,1,1)
        if(ntest.eq.1)then
          print *,'Can Mel sources ',ktau,trsrc(46,57,1),trsrc(39,52,1)
          print *,'Can Mel conc lev1 ',ktau,tr(46,57,1,2),tr(39,52,1,2)
          call printa('co2 ',tr(1,1,1,2),ktau,1,ia,ib,ja,jb,357.,1.)
          call printa('rado',tr(1,1,1,1),ktau,1,ia,ib,ja,jb,0.,1.)
        endif  ! (ntest.eq.1)
      endif

      if( iO2.NE.0 ) then
        call O2vmix(updTr, O2fact )
        call trimcopy(at,ct,updTr,iO2)
      end if

      if(iSO2.gt.0) then
        SO2fact(1) = trfact
        do i=2,nSO2lev
          SO2fact(i) = g * dt / dsig(i)
        end do  !   i=1,nSO2lev
        if( iSO2.gt.0.and.nSO2lev.lt.1 )stop 'vertmix: nSO2lev.lt.1'
      endif

      if( iSO2.NE.0 ) then
        call SO2SFlux
        call SO2vmix(updTr )
        call trimcopy(at,ct,updTr,iSO2)
        if(ntest.eq.1)then
          smax=-10.e10
          smin=10.e10
          do kq=1,ijk
           if( tr(kq,1,1,iSO2).gt.smax) then
             smax= tr(kq,1,1,iSO2)
             imm=kq
           endif
           if( tr(kq,1,1,iSO2).lt.smin) then
             smin= tr(kq,1,1,iSO2)
             immm=kq
           endif
          enddo
          print *,'max,min',smax,imm,smin,immm
        endif  ! (ntest.eq.1)
        do kq=1,ijk
         tr(kq,1,1,iSO2) = max(-100.,min(tr(kq,1,1,iSO2),9000.))
        enddo
      endif

      if( iSO4.NE.0 ) then
*       call SO4SFlux
*       call SO4vmix( updTr )
*       call trimcopy(at,ct,updTr,iSO4)
      end if

      if( iCH4.NE.0 ) then
*       call CH4SFlux
*       call CH4vmix( updTr )
*       call trimcopy(at,ct,updTr,iCH4)
      end if

      if( iO2.NE.0 ) then
*       call O2SFlux
*       call O2vmix(updTr )
*       call trimcopy(at,ct,updTr,iO2)
      end if
       return
      end
*............................................................
c     not needed, except for setting DARLAM boundary values?
      Subroutine TrimCopy(at,ct,updTr,iTracer)
      include 'constant.h'
      include 'newmpar.h'
      include 'parm.h'
      include 'tracers.h'
      real updTr(il,jl,kl),at(il,jl,kl),ct(il,jl,kl)
      call trim(at,ct,updTr,0)
      do k=1,kl
       if(npanels.eq.0)then
         do j=3,jl-1   ! don't alter boundary values
          do i=3,il-1
           tr(i,j,k,iTracer)=updTr(i,j,k)
          enddo
         enddo
       else
         do iq=1,ifull   !  alter all values for C-C
          tr(iq,1,k,iTracer)=updTr(iq,1,k)
         enddo  ! iq loop
       endif
      enddo     !  k loop
      return
      end

