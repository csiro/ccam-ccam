      subroutine vertmix
!     inputs & outputs: t,u,v,qg
      use cc_mpi, only : mydiag, myid
      use diag_m
      include 'newmpar.h'
      parameter (ntest=0)
c     parameter (ipwr=1)       ! really can use (ipwr=min(1,nlocal))
c     parameter (ilnl=il**ipwr,jlnl=jl**ipwr)
      parameter (kcl_top=kl-2) ! maximum level for cloud top (conjob & vertmix)
!     parameter (kscmom=0)     ! 0 usual, 1 turns on shal_conv momentum
      parameter (ndvmod=0)     ! 0 default, 1+ for dvmod tests
!     typically tied_con=6., tied_over=2., tied_rh=.75
!     nlocal in parm.h         ! 0 local scheme, 1 or 2 non-local scheme
!     real t(ifull,kl),u(ifull,kl),v(ifull,kl),qg(ifull,kl),ps(ifull)
      include 'arrays.h'
      include 'const_phys.h'
      include 'dates.h'
      include 'indices.h'
      include 'kuocom.h'   ! also with kbsav,ktsav,convpsav,kscsea,sigksct
      include 'liqwpar.h'  ! ifullw, qfg, qlg
      include 'map.h'      ! em, f, fu, fv, etc  not needed here?
      include 'nlin.h'
      include 'morepbl.h'
      include 'parm.h'
      include 'pbl.h'
      include 'screen.h'   ! tscrn
      include 'permsurf.h'
      include 'savuvt.h'
      include 'sigs.h'
      include 'soil.h'
      include 'tracers.h'  ! ngas, nllp, ntrac
      common/nonlsav/cfrac(ifull,kl),betatt(ifull,kl),betaqt(ifull,kl) 
      common/shalrk/shalrk(ifull,6)
      common/work3/delthet(ifull,kl),
     .    thebas(ifull,kl),cu(ifull,kl),thee(ifull,kl),qs(ifull,kl)
      common/work3b/uav(ifull,kl),vav(ifull,kl)   
!     n.b. uav & vav also used by pbldif; all of work3 used by tracervmix
      common/work3c/rhs(ifull,kl)   
      common/work3f/wrk1(ijk),wrk2(ijk),wrk3(ijk) 
      common/work2/csq(ifull),delq(ifull),dvmod(ifull)
     . ,dz(ifull),dzr(ifull),fm(ifull),fh(ifull),ri(ifull),sqmxl(ifull)
     . ,x(ifull),zhv(ifull),theeb(ifull),sigsp(ifull),kbase(ifull)
     . ,ktop(ifull),dum2(ifull,3)
      real sighkap(kl),delons(kl),delh(kl),prcpv(kl)
      real sigkap(kl)
      real at(ifull,kl),au(ifull,kl),ct(ifull,kl)
      real zh(ifull,kl),tmnht(ifull,kl)
      real gt(ifull,kl),guv(ifull,kl)
      real rkm(ifull,kl),rkh(ifull,kl),rk_shal(ifull,kl)
!     equivalence (gt,rkh_nl),(guv,rkm_nl)
      equivalence (rkh,wrk1),(rkm,wrk2),(rk_shal,uav)
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
        write (6,"('uin ',12f7.2/(8x,12f7.2))") (u(idjd,k),k=1,kl) 
        write (6,"('vin ',12f7.2/(8x,12f7.2))") (v(idjd,k),k=1,kl) 
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
        rhs(:,k)=t(1:ifull,k)*sigkap(k)  ! rhs is theta here
      enddo      !  k loop


c Pre-calculate the buoyancy parameters if using qcloud scheme.
c Follow Smith's (1990) notation; gam() is HBG's notation for (L/cp)dqsdt.
c The factor of (1/sigkap)=T/theta in betatt differs from Smith's formulation
c because we use theta derivative rather than T.
      if(nvmix.gt.0.and.nvmix.lt.4)then
        delta=1./epsil-1.
        do k=1,kl
          do iq=1,ifull
            es=establ(t(iq,k))
            pk=ps(iq)*sig(k)
            qs(iq,k)=.622*es/(pk-es)
            dqsdt=qs(iq,k)*pk*(hl/rvap)/(t(iq,k)**2*(pk-es))
            rhs(iq,k)=rhs(iq,k) !Convert to thetal - used only to calc Ri
     &                  -(hlcp*qlg(iq,k)+hlscp*qfg(iq,k))*sigkap(k)   
            betat=1./t(iq,k)
            qc=qlg(iq,k)+qfg(iq,k)
c           if(qc.gt.1.e-12)then
c             fice=qfg(iq,k)/qc
c           else
c             fice=0.
c           endif
	     fice=qfg(iq,k)/max(qc,1.e-12)
            betaq=delta/(1.+delta*qg(iq,k)-qc)
!           al=1./(1.+gam(iq,k))
            al=1./(1.+hlcp*dqsdt)
!           betac=cfrac(iq,k)*sigcll  ! sigcll ~1 is mult factor if needed
            betac=cfrac(iq,k)
     &           *al * ((hlcp+fice*hlfcp)*betat - betaq/(1.-epsil) )
            betatt(iq,k)=(betat-dqsdt*betac)/sigkap(k)           !Beta_t_tilde
            betaqt(iq,k)=betaq+betac                             !Beta_q_tilde
!           jlm: expect sigkap factor should be in betat not betatt	     
         enddo   ! iq loop
	  if(ntest.eq.2.and.mydiag)then
	     iq=idjd
            es=establ(t(iq,k))
            pk=ps(iq)*sig(k)
            qs(iq,k)=.622*es/(pk-es)
            dqsdt=qs(iq,k)*pk*(hl/rvap)/(t(iq,k)**2*(pk-es))
            rhs(iq,k)=rhs(iq,k) !Convert to thetal - used only to calc Ri
     &                  -(hlcp*qlg(iq,k)+hlscp*qfg(iq,k))*sigkap(k)   
            t=1./t(iq,k)
            qc=qlg(iq,k)+qfg(iq,k)
            if(qc.gt.1.e-12)then
              fice=qfg(iq,k)/qc
            else
              fice=0.
            endif
            betaq=delta/(1+delta*qg(iq,k)-qc)
!           al=1./(1.+gam(iq,k))
            al=1./(1.+hlcp*dqsdt)
            betac=cfrac(iq,k)
     &           *al * ((hlcp+fice*hlfcp)*betat - betaq/(1.-epsil) )
            betatt(iq,k)=(betat-dqsdt*betac)/sigkap(k)           !Beta_t_tilde
            betaqt(iq,k)=betaq+betac                             !Beta_q_tilde
	     print *,'k,qg,qs,qlg,qfg,qc ',
     &               k,qg(iq,k),qs(iq,k),qlg(iq,k),qfg(iq,k),qc
            print *,'t,rhs,cfrac,fice ',
     &               t(iq,k),rhs(iq,k),cfrac(iq,k),fice
            print *,'al,delta,betaq,betac ',al,delta,betaq,betac
	     print *,'betat,betatt,betaqt ',
     &               betat,betatt(iq,k),betaqt(iq,k)
	  endif   ! (ntest.eq.2)
        enddo    !  k loop
      else       ! still need qs()
        do k=1,kl
         do iq=1,ifull
          es=establ(t(iq,k))
          qs(iq,k)=.622*es/(ps(iq)*sig(k)-es)
         enddo   ! iq loop
        enddo    !  k loop
      endif      ! (nvmix.gt.0.and.nvmix.lt.4)

      do k=1,kl-1
       delthet(:,k)=rhs(:,k+1)-rhs(:,k)  ! rhs is theta here
      enddo      !  k loop

      if(ktau.eq.1)then
!       set ksctop for shallow convection
        ksctop=1    ! ksctop will be first level below sigkcst
c       if(abs(ksc).gt.92)then  !  abs from 16/3/04
          do while(sig(ksctop+1).gt.sigksct)  !  e.g. sigksct=.75
           ksctop=ksctop+1
          enddo
c	 endif  !(abs(ksc).gt.92)
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
         prcpv(k)=sig(k)**(-roncp)
        enddo    ! k loop
      endif

c     ******* section for Geleyn shallow convection; others moved lower *******
      if(ksc.eq.-99)then
        do k=kscbase,ksctop    ! new usage of ksc thu  02-17-2000
         do iq=1,ifull
          delthet(iq,k)=delthet(iq,k)-hlcp*max(0.,
     .                 qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k) )
         enddo  ! iq loop
        enddo   !  k loop
      endif     ! (ksc.eq.-99)
c     *********** end of Geleyn shallow convection section *****************

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
       enddo ! iq loop

c       x is bulk ri *(dvmod **2)
        if(nvmix.gt.0.and.nvmix.lt.4)then  ! new one allowing for cloudy air
           if(nvmix.eq.1)w1=dsig(k+1)/(dsig(k)+dsig(k+1)) 
           if(nvmix.eq.2)w1=.5                          
           if(nvmix.eq.3)w1=1.      !weight for lower level             
           w2=1.-w1                 !weight for upper level
          do iq=1,ifull
c          dz1=zh(iq,k)-zh(iq,k-1) !thickness of layer k
c          dz2=zh(iq,k+1)-zh(iq,k) !thickness of layer k+1
c          w1=dz1/(dz1+dz2)        !weight for lower level
           dqtot=qg(iq,k+1)+qlg(iq,k+1)+qfg(iq,k+1)
     &       -  (qg(iq,k)  +qlg(iq,k)  +qfg(iq,k))
           x(iq)=grav*dz(iq)*(
     &          (w1*betatt(iq,k)+w2*betatt(iq,k+1))*delthet(iq,k) +
     &          (w1*betaqt(iq,k)+w2*betaqt(iq,k+1))*dqtot )
          enddo ! iq loop	 
	   if(ntest.eq.4.and.k.le.9.and.ktau.eq.ntau)then
	     diffmax=0.
	     do iq=1,ifull
	      rhsk=t(iq,k)*sigkap(k)
	      rhskp=t(iq,k+1)*sigkap(k+1)
	      delthet_old=rhs(iq,k+1)-rhs(iq,k)
             xold=grav*dz(iq)*(delthet_old/(tmnht(iq,k)*sighkap(k))
     .             +.61*(qg(iq,k+1)-qg(iq,k)))
             diff=abs(xold-x(iq))
	      if(diff.gt.diffmax)then
	        diffmax=diff
		 iqmax=iq
	      endif
             write(47,'(3g13.4,i7,i4)') xold,x(iq),diff,iq,k
	     enddo
	     print *,'k,iqmax,diffmax ',k,iqmax,diffmax
	   endif	   
          rhs(:,k)=t(1:ifull,k)*sigkap(k)   !need to re-set theta for nvmix=1-3
	 elseif(nvmix.eq.5)then        ! non-cloudy x with qfg, qlg
          x(:)=grav*dz(:)*(delthet(:,k)/(tmnht(:,k)*sighkap(k))
     &             +.61*(qg(1:ifull,k+1)-qg(1:ifull,k))
     &                  -qfg(1:ifull,k+1)-qlg(1:ifull,k+1)+
     &                   qfg(1:ifull,k)+qlg(1:ifull,k) )
	 else                 ! original non-cloudy x, nvmix=4
          x(:)=grav*dz(:)*(delthet(:,k)/(tmnht(:,k)*sighkap(k))
     .            +.61*(qg(1:ifull,k+1)-qg(1:ifull,k)))
        endif     ! (nvmix.gt.0.and.nvmix.lt.4) .. else ..

c       fm and fh denote f(Louis style)*dvmod
c       nb. an error exists in the Louis expression for c; this is corrected
c       in the current code
        csq(:) = zhv(:)*(((1.+dz(:)*zhv(:))**(1./3.)-1.)*dzr(:))**3

c        newest code, stable same as csiro9 here (originally for nvmix=4)
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

c      calculate k's, guv and gt
c      nonlocal scheme usually applied to temperature and moisture, not momentum
c      (i.e. local scheme is applied to momentum for nlocal=0,1)
       if(nvmix.ne.0)then    ! use nvmix=0 only for special test runs
         do iq=1,ifull
          rkm(iq,k)=fm(iq)*sqmxl(iq)*dzr(iq)
          rkh(iq,k)=fh(iq)*sqmxl(iq)*dzr(iq)
         enddo   ! iq loop
       else
	 rkm(:,:)=0.
	 rkh(:,:)=0.
       endif     ! (nvmix.ne.0)

       if((diag.or.ntest.ge.1).and.mydiag)then
	 iq=idjd
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
         dqtot=qg(iq,k+1)+qlg(iq,k+1)+qfg(iq,k+1)
     &     -  (qg(iq,k)  +qlg(iq,k)  +qfg(iq,k))
         print *,'qfg,qfg+ ',qfg(idjd,k),qfg(idjd,k+1)
         print *,'qlg,qlg+,dqtot ',qlg(idjd,k),qlg(idjd,k+1),dqtot
	 print *,'cfrac,cfrac+ ',cfrac(idjd,k),cfrac(idjd,k+1)
	 print *,'betatt,betatt+ ',betatt(idjd,k),betatt(idjd,k+1)
	 print *,'betaqt,betaqt+ ',betaqt(idjd,k),betaqt(idjd,k+1)
         print *,'x,zh,tmnht,dz,delh,sighkap ',x(idjd),zh(idjd,k),
     .                 tmnht(idjd,k),dz(idjd),delh(k),sighkap(k)
       endif     ! (diag.or.ntest.ge.1)
      enddo      ! end of k loop

      if( (diag.or.ntest.ge.1) .and. mydiag )then
        print *,'before possible call to pbldif in vertmix'
        write (6,"('rkh0 ',16f8.3)") (rkh(idjd,k),k=1,16)
        write (6,"('rkm0 ',16f8.3)") (rkm(idjd,k),k=1,16)
        write (6,"('uav ',12f7.2/(8x,12f7.2))") (uav(idjd,k),k=1,kl) 
        write (6,"('vav ',12f7.2/(8x,12f7.2))") (vav(idjd,k),k=1,kl) 
        write (6,"('thet',12f7.2/(8x,12f7.2))") (rhs(idjd,k),k=1,kl) 
        write (6,"('t   ',12f7.2/(8x,12f7.2))") (t(idjd,k),k=1,kl) 
        write (6,"('qg ',12f7.3/(8x,12f7.3))") 
     .             (1000.*qg(idjd,k),k=1,kl)
        write (6,"('qs ',12f7.3/(8x,12f7.3))") 
     .             (1000.*qs(idjd,k),k=1,kl)
        write (6,"('thee',12f7.2/(8x,12f7.2))") 
     .        (prcpv(k)*t(idjd,k)*(t(idjd,k) + .5*hlcp*qs(idjd,k))
     .                   /(t(idjd,k) - .5*hlcp*qs(idjd,k)),k=1,kl)
      endif

      if(nlocal.ne.0)then
        call pbldif(rhs,rkh,rkm,uav,vav)
!       n.b. *** pbldif partially updates qg and theta (t done during trim)	 
!       and updates rkh and rkm arrays
        if( (diag.or.ntest.ge.1) .and. mydiag )then
	   print *,'after pbldif in vertmix'
          write (6,"('rkh1 ',16f8.3)") (rkh(idjd,k),k=1,16)
          write (6,"('rkm1 ',16f8.3)") (rkm(idjd,k),k=1,16)
          write (6,"('thet',12f7.2/(8x,12f7.2))") (rhs(idjd,k),k=1,kl) 
          write (6,"('qg ',12f7.3/(8x,12f7.3))") 
     .              (1000.*qg(idjd,k),k=1,kl)
        endif
        if(diag)then
          call printa('rkh ',rkh,ktau,nlv,ia,ib,ja,jb,0.,1.)
          call printa('cond',condx,ktau,1,ia,ib,ja,jb,0.,1.)
          call printa('zs  ',zs,ktau,1,ia,ib,ja,jb,0.,1.)
        endif
      endif      ! (nlocal.gt.0)

      rk_shal(:,:)=0.
c     ************ section for jlm shallow convection v4 *******************
      if(ksc.eq.81)then
        do k=1,ksctop-1   ! or ksctop?  
         do iq=1,ifull
          if(sig(ktsav(iq)).gt.sig_ct.and.k.lt.ktsav(iq).  
     .                                   and.condc(iq).eq.0.)then  
	     rk_shal(iq,k)=tied_con
	     rk_shal(iq,k+1)=tied_over
          endif ! (sig(ktsav(iq)).gt.sig_ct.and.k.lt.ktsav(iq).and....)
         enddo  ! iq loop
        enddo   !  k loop
      endif     ! (ksc.eq.81)
      if(ksc.eq.82)then
        do k=1,ksctop-1    
         do iq=1,ifull
          if(sig(ktsav(iq)).gt.sig_ct.and.k.lt.ktsav(iq).
     .                and.k.ge.kbsav(iq).and.condc(iq).eq.0.)then  
	     rk_shal(iq,k)=tied_con
	     rk_shal(iq,k+1)=tied_over
          endif ! (sig(ktsav(iq)).gt.sig_ct.and.k.lt.ktsav(iq).and....)
         enddo  ! iq loop
        enddo   !  k loop
      endif     ! (ksc.eq.82)
      if(ksc.eq.91)then
        do k=1,ksctop-1   ! or ksctop?  
         do iq=1,ifull
c**       if(ktsav(iq).lt.0.and.k.lt.abs(ktsav(iq)))then
          if(ktsav(iq).lt.kl.and.k.lt.ktsav(iq))then  ! April 04
	     rk_shal(iq,k)=tied_con
	     rk_shal(iq,k+1)=tied_over
          endif ! (ktsav(iq).lt.0.and.k.lt.abs(ktsav(iq)))
         enddo  ! iq loop
        enddo   !  k loop
      endif     ! (ksc.eq.91)
      if(ksc.eq.92)then
        if(nlocal.ne.0.and.sigkscb.lt.-1.)then
          do iq=1,ifull
           if(zh(iq,kbsav(iq)).gt.pblh(iq))then 
	      kbsav(iq)=kl
	    endif
          enddo
        endif  ! (nlocal.ne.0.and.sigkscb.lt.-1.)
        do k=1,ksctop-1    
         do iq=1,ifull
c**       if(ktsav(iq).lt.0.and.k.lt.abs(ktsav(iq)).and.k.ge.kbsav(iq))then
          if(ktsav(iq).lt.kl.and.k.lt.ktsav(iq).and.k.ge.kbsav(iq))then  ! April 04
	     rk_shal(iq,k)=tied_con
	     rk_shal(iq,k+1)=tied_over
          endif ! (ktsav(iq).lt.0.and.k.lt.abs(ktsav(iq)).and....)
         enddo  ! iq loop
        enddo   !  k loop
      endif     ! (ksc.eq.92)
c     *********** end of jlm shallow convection section *****************

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
            thee(iq,k)=prcpv(k)*t(iq,k)*(t(iq,k) + .5*hlcp*qs(iq,k))
     .                           /(t(iq,k) - .5*hlcp*qs(iq,k))
            if(qg(iq,kscbase).gt.tied_rh*qs(iq,kscbase).
     .       and.thee(iq,k).lt.theeb(iq))then         !  default tied_rh=.75
              rk_shal(iq,k-1)=tied_con                !  m**2/sec  6., originally 10.
              rk_shal(iq,k)=tied_over                 !  m**2/sec
            endif ! (qg(iq,kscbase).gt.rhscon*qs(iq,kscbase).....
           enddo  ! iq loop
          enddo   ! end of k=kscbase+1,ksctop loop
        else      !  i.e. Tiedtke original scheme over land and sea
          do k=kscbase+1,ksctop  ! typically kscbase=3 & ksctop=6
           do iq=1,ifull
            thee(iq,k)=prcpv(k)*t(iq,k)*(t(iq,k) + .5*hlcp*qs(iq,k))
     .                           /(t(iq,k) - .5*hlcp*qs(iq,k))
            if(qg(iq,kscbase).gt.tied_rh*qs(iq,kscbase).
     .       and.thee(iq,k).lt.theeb(iq))then         !  default tied_rh=.75
              rk_shal(iq,k-1)=tied_con                !  m**2/sec  6., originally 10.
              rk_shal(iq,k)=tied_over                 !  m**2/sec
		if(ntest.eq.3.and.k.eq.ksctop)then
		  print *,'ktau,iq,theeb,thee,delthee ',
     .       		 ktau,iq,theeb(iq),thee(iq,k),theeb(iq)-thee(iq,k)
		endif
            endif ! (qg(iq,kscbase).gt.rhscon*qs(iq,kscbase).....
           enddo  ! iq loop
          enddo   ! end of k=kscbase+1,ksctop loop
        endif     ! (kscsea.eq.1)  .. else ..
      endif       ! (ksc.eq.99)
c     *********** end of Tiedtke shallow convection section *****************

c     ************ Tiedtke_ldr-style shallow convection 93-96 ***************
      if(ksc.ge.93.and.ksc.le.96)then   
c       Calculate LCL for near surface air
c       Assume qstar=qtg(iq,1) but allow some sub-grid variability
c       The formula for tsp is eqn (21) of Bolton (1980), MWR 108, 1046-1053.
        if(sigkscb.gt.0.)then   ! uses qg1
          do iq=1,ifull
!          Leon used factor 1.01 over sea; here div by tied_rh (e.g. .99)	 
           epart=qg(iq,1)*.01*ps(iq)/(tied_rh*epsil) !in hPa 
           tsp=2840./(3.5*log(tscrn(iq))-log(epart)-4.805) + 55.   
	    sigsp(iq)=(tsp/tscrn(iq))**(1./roncp)
	   enddo
	 else                    ! uses qgscrn
          do iq=1,ifull
!          Leon used factor 1.01 over sea; here div by tied_rh (e.g. .99)	 
           epart=qgscrn(iq)*.01*ps(iq)/(tied_rh*epsil) !in hPa   
           tsp=2840./(3.5*log(tscrn(iq))-log(epart)-4.805) + 55.   
	    sigsp(iq)=(tsp/tscrn(iq))**(1./roncp)
	   enddo
        endif  ! (sigkscb.gt.0.) .. else ..
c       Look for the lowest layer s.t. the top of the layer is above the LCL,
c       and call that layer the cloud base. 
        kbase(:)=kl
        do k=ksctop-1,1,-1
         do iq=1,ifull
          if(sigsp(iq).gt.sigmh(k+1))kbase(iq)=k
         enddo
        enddo
        if(nlocal.ne.0)then  ! like old ksc=95, but ensures LCL within PBL
          do iq=1,ifull
           if(zh(iq,kbase(iq)).gt.pblh(iq))then ! had missing iq
	      kbase(iq)=kl
	    endif
          enddo  ! iq loop
        endif  ! (nlocal.ne.0)
c       following has some jlm shortcuts	 
	 do iq=1,ifull
	  k=kbase(iq)
	  qbas=qs(iq,k)
         theeb(iq)=prcpv(k)*t(iq,k)*(t(iq,k) + .5*hlcp*qbas)
     .                             /(t(iq,k) - .5*hlcp*qbas)
        enddo  ! iq loop
        ktop(:)=kbase(:)
        do k=2,ksctop+1
         do iq=1,ifull
          thee(iq,k)=prcpv(k)*t(iq,k)*(t(iq,k) + .5*hlcp*qs(iq,k))
     .                         /(t(iq,k) - .5*hlcp*qs(iq,k))
          if(theeb(iq).gt.thee(iq,k).and.ktop(iq).eq.k-1)then
            ktop(iq)=k       ! also checking if contiguous
	   endif
         enddo  ! iq loop
        enddo   ! k loop	 
        if(ksc.eq.94)then  ! from April same as 93* (93 back to older checks)
          do k=2,ksctop
           do iq=1,ifull
c            if(ktop(iq).gt.kbase(iq).and.k.le.ktop(iq).
c     .                               and.k.gt.kbase(iq))then
            if(k.gt.kbase(iq).and.k.le.ktop(iq))then
              rk_shal(iq,k-1)=tied_con              !  m**2/sec  6., originally 10.
              rk_shal(iq,k)=tied_over               !  m**2/sec
	     endif
           enddo  ! iq loop
          enddo   ! k loop
        endif     ! (ksc.eq.94)
	 if(ksc.eq.93)then
          do k=2,ksctop
           do iq=1,ifull
            if(ktop(iq).gt.kbase(iq).and.k.le.ktop(iq).
c**  .            and.k.gt.kbase(iq))then   !   this used for 93*
     .            and.k.gt.kbase(iq).and.ktop(iq).le.ksctop)then
              rk_shal(iq,k-1)=tied_con              !  m**2/sec  6., originally 10.
              rk_shal(iq,k)=tied_over               !  m**2/sec
	     endif
           enddo  ! iq loop
          enddo   ! k loop
        endif     ! (ksc.eq.93)
	 if(ksc.eq.95)then   ! new from 7 April
          do k=2,ksctop
           do iq=1,ifull
            if(k.gt.kbase(iq).and.k.le.ktop(iq).and.condc(iq).eq.0.)then
              rk_shal(iq,k-1)=tied_con              !  m**2/sec  6., originally 10.
              rk_shal(iq,k)=tied_over               !  m**2/sec
	     endif
           enddo  ! iq loop
          enddo   ! k loop
        endif     ! (ksc.eq.95)
	 if(ksc.eq.96)then   ! applied from sfce up
          do k=2,ksctop
           do iq=1,ifull
c**         if(ktop(iq).gt.kbase(iq).and.k.le.ktop(iq).
c**  .                               and.ktop(iq).le.ksctop)then
            if(ktop(iq).gt.kbase(iq).and.k.le.ktop(iq).  
     .                                   and.condc(iq).eq.0.)then  
              rk_shal(iq,k-1)=tied_con              !  m**2/sec  6., originally 10.
              rk_shal(iq,k)=tied_over               !  m**2/sec
	     endif
           enddo  ! iq loop
          enddo   ! k loop
        endif     ! (ksc.eq.96)
	 if(nmaxpr.eq.1.and.mydiag)then
	   iq=idjd
	   print *,'ksc,kbase,ktop,theeb ',
     .             ksc,kbase(idjd),ktop(idjd),theeb(idjd)
          epart=qg(iq,1)*.01*ps(iq)/(tied_rh*epsil) !in hPa 
	   if(sigkscb.lt.0.)epart=qgscrn(iq)*.01*ps(iq)/(tied_rh*epsil)
          tsp=2840./(3.5*log(tscrn(iq))-log(epart)-4.805) + 55.   
	   print *,'iq,epart,tsp,tscrn.sigsp ',
     .             iq,epart,tsp,tscrn(iq),sigsp(iq)
c         print *'t(1-6) ',(t(iq,k),k=1,6)
c         print *'qg(1-6) ',(qg(iq,k),k=1,6)
c         print *'qs(1-6) ',(qs(iq,k),k=1,6)
          print *,'thee ',(prcpv(k)*t(iq,k)*(t(iq,k) + .5*hlcp*qs(iq,k))
     .                       /(t(iq,k) - .5*hlcp*qs(iq,k)),k=1,ksctop+1)
          print *,'rk_shal ',(rk_shal(idjd,k),k=1,ksctop+2)
	 endif     ! (nmaxpr.eq.1)
      endif       ! (ksc.ge.93.and.ksc.le.96)
c     *********** end of Tiedtke_ldr shallow convection 93-96 *************

c     ************ Tiedtke_jlm shallow convection 97 ***************
      if(ksc.eq.97)then
        do k=kscbase,ksctop
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
          thee(iq,k)=prcpv(k)*t(iq,k)*(t(iq,k) + .5*hlcp*qs(iq,k))
     .                         /(t(iq,k) - .5*hlcp*qs(iq,k))
          if(theeb(iq).gt.thee(iq,k))then
            rk_shal(iq,k-1)=tied_con                !  m**2/sec  6., originally 10.
            rk_shal(iq,k)=tied_over                 !  m**2/sec
	   else
	     theeb(iq)=thebas(iq,k) ! ready for next k in k-loop
          endif
         enddo  ! iq loop
        enddo   ! k loop
c!       suppress pseudo-deep convection	
c        theeb(:)=rk_shal(:,ksctop)    ! not really theeb, just for foll. test 
c        do k=kscbase,ksctop+1     ! +1 to include tied_over points
c         do iq=1,ifull
c          if(theeb(iq).gt..9*tied_con)then  ! to allow for tied_over
c            rk_shal(iq,k)=0.
c          endif
c         enddo  ! iq loop
c        enddo   ! k loop
	 if(ntest.eq.3)then
	   do iq=1,ifull
  	    if(rk_shal(iq,ksctop-1).gt..9*tied_con)then 
	      print *,'iq,land,rk_shal ',iq,land(iq),rk_shal(iq,ksctop-1)
	    endif
	   enddo
	 endif   ! (ntest.eq.3)
      endif     ! (ksc.eq.97)
c     *********** end of Tiedtke_jlm shallow convection 97 *************

c     add in effects of shallow convection
      do k=1,kl-1  
       rkh(:,k)=rkh(:,k)+rk_shal(:,k)
      enddo   !  k loop
      if(kscmom.eq.1)then
        do k=1,kl-1  
         rkm(:,k)=rkm(:,k)+rk_shal(:,k)
        enddo   !  k loop
      endif     ! (kscmom.eq.1)
      if(nextout.eq.4)then
        do k=1,6
         shalrk(:,k)=shalrk(:,k)+rk_shal(:,k)
        enddo
      endif
      
      if(ntest.ne.0.or.diag)then
        do iq=1,ifull
	  if(rk_shal(iq,1).gt.rk_shal(iq,2))then
	    print *,'iq,rk_shal1,rk_shal2',iq,rk_shal(iq,1),rk_shal(iq,2)
	  endif
        enddo
        if (mydiag) then 
        iq=idjd
        print *,'for shallow convection in vertmix '
        print *,'ktsav,ksctop ',ktsav(idjd),ksctop
        print *,'kbase,ktop ',kbase(idjd),ktop(idjd)
        print *,'rk_shal ',(rk_shal(idjd,k),k=1,ksctop)
	 print *,'kbsav,ktsav,theeb: ',kbsav(iq),ktsav(iq),theeb(iq)
        write (6,"('rk_shal ',12f7.2/(9x,12f7.2))")
     .             (rk_shal(idjd,k),k=1,kl)
        write (6,"('rh   ',12f7.2/(5x,12f7.2))") 
     .             (100.*qg(idjd,k)/qs(idjd,k),k=1,kl)
        write (6,"('qs   ',12f7.3/(5x,12f7.3))") 
     .             (1000.*qs(idjd,k),k=1,kl)
        write (6,"('qg   ',12f7.3/(5x,12f7.3))") 
     .             (1000.*qg(idjd,k),k=1,kl)
        write (6,"('qbas ',12f7.3/(5x,12f7.3))") 
     .             (1000.*qg(idjd,k)/tied_rh,k=1,kl)
        write (6,"('t    ',12f7.2/(5x,12f7.2))") 
     .             (t(idjd,k),k=1,kl)
        write (6,"('thebas',12f7.2/(6x,12f7.2))") 
     .             (thebas(iq,k),k=1,kl)
c        write (6,"('hs  ',12f7.2/(4x,12f7.2))") 
c     .             (t(idjd,k)+hlcp*qs(idjd,k),k=1,kl)
        write (6,"('thee',12f7.2/(4x,12f7.2))") 
     .             (thee(idjd,k),k=1,kl)
        endif ! mydiag
      endif

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
      if(diag.and.mydiag)then
        print *,'vertmix rhs & qg after trim ',(rhs(idjd,k),k=1,kl)
        write (6,"('qg ',12f7.3/(8x,12f7.3))") 
     .             (1000.*qg(idjd,k),k=1,kl)
      endif

      if(ldr.ne.0)then
c       now do qfg
        rhs(1:ifull,:)=qfg(1:ifull,:)
        call trim(at,ct,rhs,0)    ! for qfg
        qfg(1:ifull,:)=rhs(1:ifull,:)
c       now do qlg
        rhs(1:ifull,:)=qlg(1:ifull,:)
        call trim(at,ct,rhs,0)    ! for qlg
        qlg(1:ifull,:)=rhs(1:ifull,:)
      endif    ! (ldr.ne.0)

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
        write (6,"('thet',12f7.2/(8x,12f7.2))") 
     .              (sigkap(k)*t(idjd,k),k=1,kl) 
        write (6,"('t   ',12f7.2/(8x,12f7.2))") (t(idjd,k),k=1,kl) 
        write (6,"('qg ',12f7.3/(8x,12f7.3))") 
     .             (1000.*qg(idjd,k),k=1,kl)
        write (6,"('u   ',12f7.2/(8x,12f7.2))") (u(idjd,k),k=1,kl) 
        write (6,"('v   ',12f7.2/(8x,12f7.2))") (v(idjd,k),k=1,kl) 
        print *,'cduv,cduv+1,cduvj+1,tss ',
     .           cduv(idjd),cduv(idjd+1),cduv(idjd+il),tss(idjd)
        write (6,"('au ',12f7.3/(8x,12f7.3))") (au(idjd,k),k=1,kl) 
        write (6,"('cu ',12f7.3/(8x,12f7.3))") (cu(idjd,k),k=1,kl) 
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

