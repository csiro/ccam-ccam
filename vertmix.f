      subroutine vertmix
!     inputs & outputs: t,u,v,qg
      use cc_mpi, only : mydiag, myid
      use diag_m
      use tkeeps ! MJT tke
      
!     rml 16/02/06 use trvmix module
      use trvmix, only : tracervmix
      include 'newmpar.h'
      parameter (ntest=0)
c     parameter (ipwr=1)       ! really can use (ipwr=min(1,nlocal))
c     parameter (ilnl=il**ipwr,jlnl=jl**ipwr)
      parameter (kcl_top=kl-2) ! maximum level for cloud top (conjob & vertmix)
!     parameter (kscmom=0)     ! 0 off, 1 turns on shal_conv momentum (usual)
      parameter (ndvmod=0)     ! 0 default, 1+ for dvmod tests
!     typically tied_con=6., tied_over=2., tied_rh=.75
!     nlocal in parm.h         ! 0 local scheme, 1 or 2 non-local scheme
!     real t(ifull,kl),u(ifull,kl),v(ifull,kl),qg(ifull,kl),ps(ifull)
      include 'arrays.h'
      include 'const_phys.h'
      include 'dates.h'
      include 'extraout.h' ! ustar ! MJT tke
      include 'indices.h'
      include 'kuocom.h'   ! also with kbsav,ktsav,convpsav,kscsea,sigksct
      include 'liqwpar.h'  ! ifullw, qfg, qlg
      include 'map.h'      ! em, f, fu, fv, etc  not needed here?
      include 'mpif.h'
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
      common/cfrac/cfrac(ifull,kl)
      real betatt(ifull,kl),betaqt(ifull,kl),rhs(ifull,kl),dqtot(ifull)
      common/work3/delthet(ifull,kl),
     &    thebas(ifull,kl),cu(ifull,kl),thee(ifull,kl),qs(ifull,kl)
      common/work3b/uav(ifull,kl),vav(ifull,kl)   
!     n.b. uav & vav also used by pbldif; all of work3 used by tracervmix
      common/work3f/wrk1(ijk),wrk2(ijk),wrk3(ijk) 
      real csq(ifull),dvmod(ifull),dz(ifull),dzr(ifull),
     &     fm(ifull),fh(ifull),sqmxl(ifull),
     &     x(ifull),zhv(ifull),theeb(ifull),sigsp(ifull)
      integer kbase(ifull),ktop(ifull)
      real sighkap(kl),sigkap(kl),delons(kl),delh(kl),prcpv(kl)
      real at(ifull,kl),au(ifull,kl),ct(ifull,kl)
      real zh(ifull,kl),tmnht(ifull,kl)
      real gt(ifull,kl),guv(ifull,kl),ri(ifull,kl)
      real rkm(ifull,kl),rkh(ifull,kl),rk_shal(ifull,kl)
      real condrag
      equivalence (rkh,wrk1),(rkm,wrk2),(rk_shal,uav)
      equivalence (tmnht,at,un),(zh,au,wrk3)
!     equivalence (gamat,ct)
c     set coefficients for Louis scheme
      data bprm/4.7/,cm/7.4/,ch/5.3/,amxlsq/100./,vkar3/.35/,vkar4/.4/
      data bprmj/5./,cmj/5./,chj/2.6/
      save kscbase,ksctop,prcpv
      include 'establ.h'
      real zg(ifull,kl) ! MJT tke

      rong=rdry/grav
      do k=1,kl-1
       sighkap(k)=sigmh(k+1)**(-roncp)
       delons(k)=rong *((sig(k+1)-sig(k))/sigmh(k+1))
      enddo      ! k loop
      do k=1,kl
       delh(k)=-rong *dsig(k)/sig(k)  ! sign of delh defined so always +ve
       sigkap(k)=sig(k)**(-roncp)
      enddo      ! k loop
      if(diag.or.ntest>=1)then
        call maxmin(u,'%u',ktau,1.,kl)
        call maxmin(v,'%v',ktau,1.,kl)
        call maxmin(t,'%t',ktau,1.,kl)
        call maxmin(qg,'qg',ktau,1.e3,kl)     
        call MPI_Barrier( MPI_COMM_WORLD, ierr ) ! stop others going past
      if(mydiag)then
          print *,'sig ',sig
          print *,'dsig ',dsig
          print *,'delh ',delh
          print *,'in vertmix'
          write (6,"('uin ',9f8.3/4x,9f8.3)") u(idjd,:) 
          write (6,"('vin ',9f8.3/4x,9f8.3)") v(idjd,:) 
        endif
      endif
      rlogs1=log(sig(1))
      rlogs2=log(sig(2))
      rlogh1=log(sigmh(2))
      rlog12=1./(rlogs1-rlogs2)
      tmnht(1:ifull,1)=(t(1:ifull,2)*rlogs1-t(1:ifull,1)*rlogs2+
     &           (t(1:ifull,1)-t(1:ifull,2))*rlogh1)*rlog12
!     n.b. an approximate zh (in m) is quite adequate for this routine
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
      if(nmaxpr==1.and.mydiag)
     &  write (6,"('thet_in',9f8.3/7x,9f8.3)") rhs(idjd,:)

      ! moved from below ! MJT tke
      if(ktau==1.and.ksc.ne.0)then
!       set ksctop for shallow convection
        ksctop=1    ! ksctop will be first level below sigkcst
        do while(sig(ksctop+1)>sigksct)  !  e.g. sigksct=.75
         ksctop=ksctop+1
        enddo
        kscbase=1  ! kscbase will be first level above sigkcsb
        do while(sig(kscbase)>sigkscb.and.sigkscb>0.) ! e.g. sigkscb=.99
         kscbase=kscbase+1
        enddo
        if ( myid == 0 ) then
        print *,'For shallow convection:'
        print *,'ksc,kscbase,ksctop,kscsea ',
     &           ksc,kscbase,ksctop,kscsea
        write (6,"(' sigkscb,sigksct,tied_con,tied_over,tied_rh:',
     &       5f8.3)")sigkscb,sigksct,tied_con,tied_over,tied_rh
        end if
        do k=1,kl
         prcpv(k)=sig(k)**(-roncp)
        enddo  ! k loop
      endif    ! (ktau==1.and.ksc.ne.0)


      if (nvmix.ne.6) then ! usual ! MJT tke

c Pre-calculate the buoyancy parameters if using qcloud scheme.
c Follow Smith's (1990) notation; gam() is HBG's notation for (L/cp)dqsdt.
c The factor of (1/sigkap)=T/theta in betatt differs from Smith's formulation
c because we use theta derivative rather than (dry static energy)/cp.
      if(nvmix>0.and.nvmix<4)then
        delta=1./epsil-1.  ! i.e. 1/.622 -1., i.e. .6077
        do k=1,kl
         if(sig(k)>.8)then ! change made 17/1/06
           do iq=1,ifull
            es=establ(t(iq,k))
            pk=ps(iq)*sig(k)
            qs(iq,k)=.622*es/(pk-es)  
            dqsdt=qs(iq,k)*pk*(hl/rvap)/(t(iq,k)**2*(pk-es))
            rhs(iq,k)=rhs(iq,k) !Convert to thetal - used only to calc Ri
     &                  -(hlcp*qlg(iq,k)+hlscp*qfg(iq,k))*sigkap(k)   
            betat=1./t(iq,k)
            qc=qlg(iq,k)+qfg(iq,k)
            fice=qfg(iq,k)/max(qc,1.e-12)
            betaq=delta/(1.+delta*qg(iq,k)-qc)
            al=1./(1.+hlcp*dqsdt)
            betac=cfrac(iq,k)
     &           *al * ((hlcp+fice*hlfcp)*betat - betaq/(1.-epsil) )
            betatt(iq,k)=(betat-dqsdt*betac)/sigkap(k)  !Beta_t_tilde
            betaqt(iq,k)=betaq+betac                    !Beta_q_tilde
           enddo   ! iq loop
         else  ! i.e. (sig(k)<.8)
           do iq=1,ifull
            es=establ(t(iq,k))
            pk=ps(iq)*sig(k)
            qs(iq,k)=.622*es/max(1.,pk-es)  ! still need qs(); max for k=kl
            betat=1./t(iq,k)
c           qc=qlg(iq,k)+qfg(iq,k)
c           betaq=delta/(1.+delta*qg(iq,k)-qc)
            betaq=delta/(1.+delta*qg(iq,k))
            betatt(iq,k)=betat/sigkap(k)    !Beta_t_tilde
            betaqt(iq,k)=betaq              !Beta_q_tilde
           enddo   ! iq loop
         endif  ! (sig(k)>.8)
         if(diag.and.mydiag)then
            iq=idjd
            dqsdt=qs(iq,k)*pk*(hl/rvap)/(t(iq,k)**2*max(pk-es,1.))
            betat=1./t(iq,k)
            qc=qlg(iq,k)+qfg(iq,k)
            fice=qfg(iq,k)/max(qc,1.e-12)
            betaq=delta/(1+delta*qg(iq,k)-qc)
!           al=1./(1.+gam(iq,k))
            al=1./(1.+hlcp*dqsdt)
            betac=cfrac(iq,k)
     &           *al * ((hlcp+fice*hlfcp)*betat - betaq/(1.-epsil) )
            print *,'k,qg,qs,qlg,qfg,qc ',
     &               k,qg(iq,k),qs(iq,k),qlg(iq,k),qfg(iq,k),qc
            print *,'t,rhs,cfrac,fice ',
     &               t(iq,k),rhs(iq,k),cfrac(iq,k),fice
            print *,'al,delta,betaq,betac ',al,delta,betaq,betac
            print *,'betat,betatt,betaqt ',
     &               betat,betatt(iq,k),betaqt(iq,k)
         endif   ! (ntest==2)
        enddo    !  k loop
      else       ! other nvmix values (0 or 4+) still need qs()
        do k=1,kl
         do iq=1,ifull
          es=establ(t(iq,k))
          qs(iq,k)=.622*es/max(1.,ps(iq)*sig(k)-es)  ! max for k=kl
         enddo   ! iq loop
        enddo    !  k loop
      endif      ! (nvmix>0.and.nvmix<4)

      do k=1,kl-1
       delthet(:,k)=rhs(:,k+1)-rhs(:,k)  ! rhs is theta or thetal here
      enddo      !  k loop

      ! moved above ! MJT the
!      if(ktau==1.and.ksc.ne.0)then
!!       set ksctop for shallow convection
!        ksctop=1    ! ksctop will be first level below sigkcst
!        do while(sig(ksctop+1)>sigksct)  !  e.g. sigksct=.75
!         ksctop=ksctop+1
!        enddo
!        kscbase=1  ! kscbase will be first level above sigkcsb
!        do while(sig(kscbase)>sigkscb.and.sigkscb>0.) ! e.g. sigkscb=.99
!         kscbase=kscbase+1
!        enddo
!        if ( myid == 0 ) then
!        print *,'For shallow convection:'
!        print *,'ksc,kscbase,ksctop,kscsea ',
!     &           ksc,kscbase,ksctop,kscsea
!        write (6,"(' sigkscb,sigksct,tied_con,tied_over,tied_rh:',
!     &       5f8.3)")sigkscb,sigksct,tied_con,tied_over,tied_rh
!        end if
!        do k=1,kl
!         prcpv(k)=sig(k)**(-roncp)
!        enddo  ! k loop
!      endif    ! (ktau==1.and.ksc.ne.0)

c     ****** section for Geleyn shallow convection; others moved lower****
      if(ksc==-99)then
        do k=kscbase,ksctop    ! new usage of ksc thu  02-17-2000
         do iq=1,ifull
          delthet(iq,k)=delthet(iq,k)-hlcp*max(0.,
     &                 qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k) )
         enddo  ! iq loop
        enddo   !  k loop
      endif     ! (ksc==-99)
      if(ksc==-98)then    ! modified Geleyn Jan '08
        do k=kscbase,ksctop    ! new usage of ksc thu  02-17-2000
         do iq=1,ifull
          if(qg(iq,k)>tied_rh*qs(iq,k).or.
     &       qg(iq,k+1)>tied_rh*qs(iq,k+1))then
             delthet(iq,k)=delthet(iq,k)-hlcp*max(0.,
     &                     qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k) )
          endif
         enddo  ! iq loop
        enddo   !  k loop
      endif     ! (ksc==-98)
      if(ksc==-97)then    ! modified Geleyn Jan '08
        do k=kscbase,ksctop    ! new usage of ksc thu  02-17-2000
         do iq=1,ifull
          if(qg(iq,k)>tied_rh*qs(iq,k))then
            delthet(iq,k)=delthet(iq,k)-hlcp*max(0.,
     &                    qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k) )
          endif
         enddo  ! iq loop
        enddo   !  k loop
      endif     ! (ksc==-97)
      if(ksc==-96)then   ! combination of Geleyn and jlm 83 (Feb 08)
        do k=1,ksctop    
         do iq=1,ifull
          if(k<ktsav(iq).and.k>=kbsav(iq).and.condc(iq)==0.)then  
            delthet(iq,k)=delthet(iq,k)-hlcp*max(0.,
     &                 qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k) )
            print *,'-96 iq,k,kbsav,ktsav ',iq,k,kbsav(iq),ktsav(iq),
     &      hlcp*(qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k)),delthet(iq,k)
          endif 
         enddo  ! iq loop
        enddo   !  k loop
      endif     ! (ksc==-96)
      if(ksc==-95)then ! same as -99 but has tied_rh (e.g. .75) factor
c       kshal(:)=0
        do k=kscbase,ksctop     
         do iq=1,ifull
          delthet(iq,k)=delthet(iq,k)-tied_rh*hlcp*max(0.,
     &                 (qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k)) )
c         if(qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k)>0.)kshal(iq)=k+1
c         if(qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k)>0.)
c    &      print*,'ktau,iq,k,diff ',
c    &        ktau,iq,k,hlcp*(qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k))
         enddo  ! iq loop
        enddo   !  k loop
      endif     ! (ksc==-95)
      if(ksc==-94)then   ! combination of Geleyn and jlm 83 (Feb 08)
        do k=1,ksctop    
         do iq=1,ifull
          if(k<ktsav(iq).and.k>=kbsav(iq).and.condc(iq)==0.)then  
            delthet(iq,k)=0.
            print *,'-94 iq,k,kbsav,ktsav ',iq,k,kbsav(iq),ktsav(iq),
     &      hlcp*(qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k)),delthet(iq,k)
          endif 
         enddo  ! iq loop
        enddo   !  k loop
      endif     ! (ksc==-94)
      if(ksc==-93)then ! single-layer (pblh) version of -95
        do iq=1,ifull
         do k=kscbase,kl/2
          if(zh(iq,k)<pblh(iq).and.zh(iq,k+1)>pblh(iq))then
c           aa=hlcp*(qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k))
c           if(aa>0.)
c    &        print *,'iq,k,zh,pblh,rh ',iq,k,zh(iq,k),pblh(iq),
c~`    &                qg(iq,k)/qs(iq,k),delthet(iq,k),aa 
            delthet(iq,k)=delthet(iq,k)-tied_rh*hlcp*max(0.,
     &                   (qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k)) )
           endif
         enddo  ! k loop 
        enddo   ! iq loop
      endif     ! (ksc==-93)
      if(ksc==-92)then ! capped-by-pblh version of -95
        do iq=1,ifull
         do k=kscbase,kl/2
          if(zh(iq,k)<pblh(iq))then
            delthet(iq,k)=delthet(iq,k)-tied_rh*hlcp*max(0.,
     &                   (qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k)) )
           endif
         enddo  ! k loop 
        enddo   ! iq loop
      endif     ! (ksc==-92)
      if(ksc==-91)then ! capped-by-pblh (anywhere in layer) version of -95
        do iq=1,ifull
         do k=2,kl/2
          if(zh(iq,k-1)<pblh(iq))then
            delthet(iq,k)=delthet(iq,k)-tied_rh*hlcp*max(0.,
     &                   (qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k)) )
           endif
         enddo  ! k loop 
        enddo   ! iq loop
      endif     ! (ksc==-91)
c     ********* end of Geleyn shallow convection section ****************

!     following now defined in vertmix (don't need to pass from sflux)
      uav(1:ifull,:)=av_vmod*u(1:ifull,:)+(1.-av_vmod)*savu(1:ifull,:) 
      vav(1:ifull,:)=av_vmod*v(1:ifull,:)+(1.-av_vmod)*savv(1:ifull,:) 
      do k=1,kl-1
       do iq=1,ifull
        dz(iq) =-tmnht(iq,k)*delons(k)  ! this is z(k+1)-z(k)
        dzr(iq)=1./dz(iq)
        zhv(iq)=1./zh(iq,k)
        if(ndvmod==0)then
          dvmod(iq)=sqrt( (uav(iq,k+1)-uav(iq,k))**2
     &                   +(vav(iq,k+1)-vav(iq,k))**2 )
        else
          dvmod(iq)=ndvmod  ! just for tests
        endif    ! (ndvmod==0)
       enddo ! iq loop

c      x is bulk ri *(dvmod **2); used to calc ri(k), rkh(k) etc
       if(nvmix>0.and.nvmix<4)then  ! new one allowing for cloudy air
!        usually nvmix=3       
         if(sig(k)>.8)then ! change made 17/1/06
           dqtot(:)=qg(1:ifull,k+1)+qlg(1:ifull,k+1)+qfg(1:ifull,k+1)
     &            -(qg(1:ifull,k)  +qlg(1:ifull,k)  +qfg(1:ifull,k))
         else
           dqtot(:)=0.
         endif
         if(nvmix==2)then !  jlm May '05
           do iq=1,ifull
            x(iq)=grav*dz(iq)*(
     &         (min(betatt(iq,k),betatt(iq,k+1)))*delthet(iq,k) +
     &         (max(betaqt(iq,k),betaqt(iq,k+1)))*dqtot(iq) )
           enddo ! iq loop	
         else    ! i.e. nvmix=1 or 3
           if(nvmix==1)w1=dsig(k+1)/(dsig(k)+dsig(k+1)) 
           if(nvmix==3)w1=1.    !weight for lower level  usual           
           w2=1.-w1             !weight for upper level
           do iq=1,ifull
            x(iq)=grav*dz(iq)*(
     &           (w1*betatt(iq,k)+w2*betatt(iq,k+1))*delthet(iq,k) +
     &           (w1*betaqt(iq,k)+w2*betaqt(iq,k+1))*dqtot(iq) )
           enddo ! iq loop	 
         endif  !  (nvmix==2) .. else ..
         if(ntest==4.and.k<=9.and.ktau==ntau)then
           diffmax=0.
           do iq=1,ifull
            rhsk=t(iq,k)*sigkap(k)
            rhskp=t(iq,k+1)*sigkap(k+1)
            delthet_old=rhs(iq,k+1)-rhs(iq,k)
            xold=grav*dz(iq)*(delthet_old/(tmnht(iq,k)*sighkap(k))
     &           +.61*(qg(iq,k+1)-qg(iq,k)))
            diff=abs(xold-x(iq))
            if(diff>diffmax)then
              diffmax=diff
              iqmax=iq
            endif
            write(47,'(3g13.4,i7,i4)') xold,x(iq),diff,iq,k
           enddo
           print *,'k,iqmax,diffmax ',k,iqmax,diffmax
         endif   ! (ntest==4.and.k<=9.and.ktau==ntau)
         rhs(:,k)=t(1:ifull,k)*sigkap(k)   !need to re-set theta for nvmix=1-3
       elseif(nvmix==5)then        ! non-cloudy x with qfg, qlg
          x(:)=grav*dz(:)*(delthet(:,k)/(tmnht(:,k)*sighkap(k))
     &             +.61*(qg(1:ifull,k+1)-qg(1:ifull,k))
     &                  -qfg(1:ifull,k+1)-qlg(1:ifull,k+1)+
     &                   qfg(1:ifull,k)+qlg(1:ifull,k) )
       else                 ! original non-cloudy x, nvmix=4
          x(:)=grav*dz(:)*(delthet(:,k)/(tmnht(:,k)*sighkap(k))
     &            +.61*(qg(1:ifull,k+1)-qg(1:ifull,k)))
       endif     ! (nvmix>0.and.nvmix<4) .. else ..

c       fm and fh denote f(Louis style)*dvmod
c       nb. an error exists in the Louis expression for c; this is corrected
c       in the current code
        csq(:) = zhv(:)*(((1.+dz(:)*zhv(:))**(1./3.)-1.)*dzr(:))**3

c        newest code, stable same as csiro9 here (originally for nvmix=4)
         do iq=1,ifull
          sqmxl(iq)=(vkar4*zh(iq,k)/(1.+vkar4*zh(iq,k)/amxlsq))**2
          dvmod(iq)=max( dvmod(iq) , 1. )
          ri(iq,k)=x(iq)/dvmod(iq)**2
          if(ri(iq,k)< 0.)then  ! unstable case
c           first do momentum
            denma=dvmod(iq)+cmj*( 2.*bprmj*sqmxl(iq)*
     &                            sqrt(-x(iq)*csq(iq)) )
            fm(iq)=dvmod(iq)-(2.*bprmj *x(iq))/denma
c           now heat
            denha=dvmod(iq)+chj*( 2.*bprmj*sqmxl(iq)*
     &                            sqrt(-x(iq)*csq(iq)) )
            fh(iq)=dvmod(iq)-(2.*bprmj *x(iq))/denha
          else                     ! stable case
c           the following is the original Louis stable formula
            fm(iq)=dvmod(iq)/(1.+4.7*ri(iq,k))**2
            fh(iq)=fm(iq)
          endif
         enddo   ! iq loop

       if(nvmix<0)then   ! use quasi-neutral mixing for test runs only
         do iq=1,ifull
          fm(iq)=1.
          fh(iq)=1.
         enddo   ! iq loop
         if(k<kl/2)then
           do iq=1,ifull
            if(land(iq))fh(iq)=abs(nvmix)
           enddo   ! iq loop
         endif
       endif     !  (nvmix<0)

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

       if((diag.or.ntest>=1).and.mydiag)then
        iq=idjd
         print *,'k,dt,delsig,sqmxl,dzr ',
     &            k,dt,delsig,sqmxl(idjd),dzr(idjd)
         print *,'k,t,t+,ps ',k,t(idjd,k),t(idjd,k+1),ps(idjd)
         print *,'k,qg,qg+ ',k,qg(idjd,k),qg(idjd,k+1)
         print *,'k,qs,qs+ ',k,qs(idjd,k),qs(idjd,k+1)
         es=establ(t(idjd,k))
         esp=establ(t(idjd,k+1))
         print *,'k,es,es+,delthet ',k,es,esp,delthet(idjd,k)
         print *,'k,fm,dvmod,ri,csq ',
     &            k,fm(idjd),dvmod(idjd),ri(idjd,k),csq(idjd)
         print *,'qfg,qfg+ ',qfg(idjd,k),qfg(idjd,k+1)
         print *,'qlg,qlg+,dqtot ',qlg(idjd,k),qlg(idjd,k+1),dqtot(iq)
         print *,'cfrac,cfrac+ ',cfrac(idjd,k),cfrac(idjd,k+1)
         print *,'betatt,betatt+,betatt*delthet ',
     &            betatt(iq,k),betatt(iq,k+1),betatt(iq,k)*delthet(iq,k)
         print *,'betaqt,betaqt+,betaqt*dqtot   ',
     &            betaqt(iq,k),betaqt(iq,k+1),betaqt(iq,k)*dqtot(iq)
         print *,'x,zh,tmnht,dz,delh,sighkap ',x(idjd),zh(idjd,k),
     &                 tmnht(idjd,k),dz(idjd),delh(k),sighkap(k)
       endif     ! (diag.or.ntest>=1)
      enddo      ! end of k loop

      if( (diag.or.ntest>=1) .and. mydiag )then
        print *,'before possible call to pbldif in vertmix'
        write (6,"('uav ',9f8.3/4x,9f8.3)") uav(idjd,:) 
        write (6,"('vav ',9f8.3/4x,9f8.3)") vav(idjd,:)
        write (6,"('t   ',9f8.3/4x,9f8.3)") t(idjd,:)
        write (6,"('qg ',3p9f8.3/4x,9f8.3)") qg(idjd,:)
        write (6,"('qs ',3p9f8.3/4x,9f8.3)") qs(idjd,:)
        do k=1,kl
         prcpv(k)=sig(k)**(-roncp)
        enddo       
        write (6,"('thee',9f8.3/4x,9f8.3)") 
     &        (prcpv(k)*t(idjd,k)*(t(idjd,k) + .5*hlcp*qs(idjd,k))
     &                   /(t(idjd,k) - .5*hlcp*qs(idjd,k)),k=1,kl)
      endif
      if(nmaxpr==1.and.mydiag)then
        write (6,"('rino_v',9f9.3/6x,9f9.3)") ri(idjd,1:kl-1)
        write (6,"('rkh0 ',9f9.3/5x,9f9.3)") rkh(idjd,1:kl-2)
        write (6,"('rkm0 ',9f9.3/5x,9f9.3)") rkm(idjd,1:kl-2)
      endif

      if(nlocal.ne.0)then
        call pbldif(rhs,rkh,rkm,uav,vav)  ! rhs is theta or thetal
!       n.b. *** pbldif partially updates qg and theta (t done during trim)	 
!       and updates rkh and rkm arrays
        if(nmaxpr==1.and.mydiag)then
          write (6,"('pblh ',f8.2)") pblh(idjd)
          write (6,"('rkh1 ',9f9.3/5x,9f9.3)") rkh(idjd,1:kl-2)
          write (6,"('rkm1 ',9f9.3/5x,9f9.3)") rkm(idjd,1:kl-2)
        endif
        if(nmaxpr==1.and.mydiag)
     &    write (6,"('thet_pbl',9f8.3/8x,9f8.3)") rhs(idjd,:)
        if( (diag.or.ntest>=1) .and. mydiag )
     &    write (6,"('qg ',3p9f8.3/4x,9f8.3)") qg(idjd,:)
        if(diag)then
          call printa('rkh ',rkh,ktau,nlv,ia,ib,ja,jb,0.,1.)
          call printa('cond',condx,ktau,1,ia,ib,ja,jb,0.,1.)
          call printa('zs  ',zs,ktau,1,ia,ib,ja,jb,0.,1.)
        endif
      endif      ! (nlocal>0)
      
      else ! tke scheme                                                 ! MJT tke
       ! note ksc < 0 options are clobbered when nvmix=6                ! MJT tke
       uav(1:ifull,:)=av_vmod*u(1:ifull,:)+(1.-av_vmod)*savu(1:ifull,:) ! MJT tke
       vav(1:ifull,:)=av_vmod*v(1:ifull,:)+(1.-av_vmod)*savv(1:ifull,:) ! MJT tke
       zg(:,1)=bet(1)*t(1:ifull,1)/grav                                 ! MJT tke
       do k=2,kl                                                        ! MJT tke
         zg(:,k)=zg(:,k-1)+(bet(k)*t(1:ifull,k)
     &                     +betm(k)*t(1:ifull,k-1))/grav                ! MJT tke
       end do                                                           ! MJT tke
       call tkemix(rkm,rhs,qg(1:ifull,:),uav,vav,pblh,
     &             rdry*fg*t(1:ifull,1)/(ps(1:ifull)*cp*sigmh(1)),
     &             rdry*eg*t(1:ifull,1)/(ps(1:ifull)*hl*sigmh(1)),
     &             ps(1:ifull),ustar,zg,sigkap,dt,av_vmod,0)            ! MJT tke
       rkh=rkm                                                          ! MJT tke
      end if ! nvmix.ne.6                                               ! MJT tke

      rk_shal(:,:)=0.
c     ***** ***** section for jlm shallow convection v4 *****************
      if(ksc==81)then
        do k=1,ksctop-1   ! or ksctop?  
         do iq=1,ifull
          if(sig(ktsav(iq))>sig_ct.and.k<ktsav(iq).  
     &                                   and.condc(iq)==0.)then  
            rk_shal(iq,k)=tied_con
            rk_shal(iq,k+1)=tied_over
          endif ! (sig(ktsav(iq))>sig_ct.and.k<ktsav(iq).and....)
         enddo  ! iq loop
        enddo   !  k loop
      endif     ! (ksc==81)
      if(ksc==82)then
        do k=1,ksctop-1    
         do iq=1,ifull
          if(sig(ktsav(iq))>sig_ct.and.k<ktsav(iq).
     &                and.k>=kbsav(iq).and.condc(iq)==0.)then  
            rk_shal(iq,k)=tied_con
            rk_shal(iq,k+1)=tied_over
          endif ! (sig(ktsav(iq))>sig_ct.and.k<ktsav(iq).and....)
         enddo  ! iq loop
        enddo   !  k loop
      endif     ! (ksc==82)
      if(ksc==83)then
        do k=1,ksctop-1    
         do iq=1,ifull
          if(sig(k)>sig_ct.and.k<ktsav(iq).
     &                and.k>=kbsav(iq).and.condc(iq)==0.)then  
            rk_shal(iq,k)=tied_con
            rk_shal(iq,k+1)=tied_over
          endif ! (sig(ktsav(iq))>sig_ct.and.k<ktsav(iq).and....)
         enddo  ! iq loop
        enddo   !  k loop
      endif     ! (ksc==83)
      if(ksc==91)then
        do k=1,ksctop ! May 08
         do iq=1,ifull
          if(ktsav(iq)<kl-1.and.k<ktsav(iq))then  ! April 04
            rk_shal(iq,k)=tied_con
            rk_shal(iq,k+1)=tied_over
          endif ! (ktsav(iq)<0.and.k<abs(ktsav(iq)))
         enddo  ! iq loop
        enddo   !  k loop
      endif     ! (ksc==91)
      if(ksc==92)then
        do k=1,ksctop ! May 08
         do iq=1,ifull
          if(k>=kbsav(iq).and.k<ktsav(iq).and.condc(iq)==0.)then  ! May 08
            rk_shal(iq,k)=tied_con
            rk_shal(iq,k+1)=tied_over
          endif ! (ktsav(iq)<0.and.k<abs(ktsav(iq)))
         enddo  ! iq loop
        enddo   !  k loop
      endif     ! (ksc==92)
c     *********** end of jlm shallow convection section *****************

c     ************ section for Tiedtke shallow convection *******************
      if(ksc==99)then
        do iq=1,ifull
        theeb(iq)=prcpv(kscbase)*t(iq,kscbase)*
     &                  (t(iq,kscbase) + .5*hlcp*qs(iq,kscbase))
     &                 /(t(iq,kscbase) - .5*hlcp*qs(iq,kscbase))
        enddo    ! iq loop
        if(kscsea==1)then  ! Tiedtke done only over sea
          do k=kscbase+1,ksctop
           do iq=1,ifull
            if(.not.land(iq))then 
             thee(iq,k)=prcpv(k)*t(iq,k)*(t(iq,k) + .5*hlcp*qs(iq,k))
     &                           /(t(iq,k) - .5*hlcp*qs(iq,k))
             if(qg(iq,kscbase)>tied_rh*qs(iq,kscbase).
     &        and.thee(iq,k)<theeb(iq))then         !  default tied_rh=.75
              rk_shal(iq,k-1)=tied_con     !  m**2/sec  6., originally 10.
              rk_shal(iq,k)=tied_over      !  m**2/sec
             endif ! (qg(iq,kscbase)>rhscon*qs(iq,kscbase).....
            endif  ! (.not.land(iq)) 
           enddo   ! iq loop
          enddo    ! end of k=kscbase+1,ksctop loop
        else       !  i.e. Tiedtke original scheme over land and sea
          do k=kscbase+1,ksctop  ! typically kscbase=3 & ksctop=6
           do iq=1,ifull
            thee(iq,k)=prcpv(k)*t(iq,k)*(t(iq,k) + .5*hlcp*qs(iq,k))
     &                           /(t(iq,k) - .5*hlcp*qs(iq,k))
            if(qg(iq,kscbase)>tied_rh*qs(iq,kscbase).
     &       and.thee(iq,k)<theeb(iq))then !  default tied_rh=.75
              rk_shal(iq,k-1)=tied_con     !  m**2/sec  6., originally 10.
              rk_shal(iq,k)=tied_over      !  m**2/sec
              if(ntest==3.and.k==ksctop)then
                print *,'ktau,iq,theeb,thee,delthee ',
     &                 ktau,iq,theeb(iq),thee(iq,k),theeb(iq)-thee(iq,k)
              endif
            endif ! (qg(iq,kscbase)>rhscon*qs(iq,kscbase).....
           enddo  ! iq loop
          enddo   ! end of k=kscbase+1,ksctop loop
        endif     ! (kscsea==1)  .. else ..
      endif       ! (ksc==99)
c     *********** end of Tiedtke shallow convection section **************

c     *********** Tiedtke_ldr-style shallow convection 93 ************
      if(ksc>=93.and.ksc<=96)then   
c       Calculate LCL for near surface air
c       Assume qstar=qtg(iq,1) but allow some sub-grid variability
c       The formula for tsp is eqn (21) of Bolton (1980), MWR 108, 1046-1053.
        if(sigkscb>0.)then   ! uses qg1
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
        endif  ! (sigkscb>0.) .. else ..
c       Look for the lowest layer s.t. the top of the layer is above the LCL,
c       and call that layer the cloud base. 
        kbase(:)=kl
        do k=ksctop-1,1,-1
         do iq=1,ifull
          if(sigsp(iq)>sigmh(k+1))kbase(iq)=k
         enddo
        enddo
        if(nlocal.ne.0)then  ! like old ksc=95, but ensures LCL within PBL
          do iq=1,ifull
           if(zh(iq,kbase(iq))>pblh(iq))then 
              kbase(iq)=kl
           endif
          enddo  ! iq loop
        endif  ! (nlocal.ne.0)
c       following has some jlm shortcuts	 
        do iq=1,ifull
         k=kbase(iq)
         qbas=qs(iq,k)
         theeb(iq)=prcpv(k)*t(iq,k)*(t(iq,k) + .5*hlcp*qbas)
     &                             /(t(iq,k) - .5*hlcp*qbas)
        enddo  ! iq loop
        ktop(:)=kbase(:)
        do k=2,ksctop+1
         do iq=1,ifull
          thee(iq,k)=prcpv(k)*t(iq,k)*(t(iq,k) + .5*hlcp*qs(iq,k))
     &                         /(t(iq,k) - .5*hlcp*qs(iq,k))
          if(theeb(iq)>thee(iq,k).and.ktop(iq)==k-1)then
            ktop(iq)=k       ! also checking if contiguous
          endif
         enddo  ! iq loop
        enddo   ! k loop	 
        if(ksc==94)then  ! from April same as 93* 
          do k=2,ksctop
           do iq=1,ifull
c            if(ktop(iq)>kbase(iq).and.k<=ktop(iq).
c     &                               and.k>kbase(iq))then
            if(k>kbase(iq).and.k<=ktop(iq))then
              rk_shal(iq,k-1)=tied_con     !  m**2/sec  6., originally 10.
              rk_shal(iq,k)=tied_over      !  m**2/sec
            endif
           enddo  ! iq loop
          enddo   ! k loop
        endif     ! (ksc==94)
        if(ksc==93)then
          do k=2,ksctop
           do iq=1,ifull
            if(ktop(iq)>kbase(iq).and.k<=ktop(iq).
c**  .            and.k>kbase(iq))then   !   this used for 93*
     &            and.k>kbase(iq).and.ktop(iq)<=ksctop)then
              rk_shal(iq,k-1)=tied_con     !  m**2/sec  6., originally 10.
              rk_shal(iq,k)=tied_over      !  m**2/sec
            endif
           enddo  ! iq loop
          enddo   ! k loop
        endif     ! (ksc==93)
        if(ksc==95)then   ! new from 7 April
          do k=2,ksctop
           do iq=1,ifull
            if(k>kbase(iq).and.k<=ktop(iq).and.condc(iq)==0.)then
              rk_shal(iq,k-1)=tied_con     !  m**2/sec  6., originally 10.
              rk_shal(iq,k)=tied_over      !  m**2/sec
            endif
           enddo  ! iq loop
          enddo   ! k loop
        endif     ! (ksc==95)
        if(ksc==96)then   ! applied from sfce up
          do k=2,ksctop
           do iq=1,ifull
c**         if(ktop(iq)>kbase(iq).and.k<=ktop(iq).
c**  .                               and.ktop(iq)<=ksctop)then
            if(ktop(iq)>kbase(iq).and.k<=ktop(iq).  
     &                                   and.condc(iq)==0.)then  
              rk_shal(iq,k-1)=tied_con     !  m**2/sec  6., originally 10.
              rk_shal(iq,k)=tied_over      !  m**2/sec
            endif
           enddo  ! iq loop
          enddo   ! k loop
        endif     ! (ksc==96)
c	 if(nmaxpr==1.and.mydiag)then
c	   iq=idjd
c	   print *,'ksc,kbase,ktop,theeb ',
c     &             ksc,kbase(idjd),ktop(idjd),theeb(idjd)
c          epart=qg(iq,1)*.01*ps(iq)/(tied_rh*epsil) !in hPa 
c	   if(sigkscb<0.)epart=qgscrn(iq)*.01*ps(iq)/(tied_rh*epsil)
c          tsp=2840./(3.5*log(tscrn(iq))-log(epart)-4.805) + 55.   
c	   print *,'iq,epart,tsp,tscrn.sigsp ',
c     &             iq,epart,tsp,tscrn(iq),sigsp(iq)
c          print *,'thee ',(prcpv(k)*t(iq,k)*(t(iq,k) + .5*hlcp*qs(iq,k))
c     &                       /(t(iq,k) - .5*hlcp*qs(iq,k)),k=1,ksctop+1)
c          print *,'rk_shal ',(rk_shal(idjd,k),k=1,ksctop+2)
c	 endif     ! (nmaxpr==1.and.mydiag)
      endif       ! (ksc>=93.and.ksc<=96)
c     *********** end of Tiedtke_ldr shallow convection 93-96 *********

c     *************** Tiedtke_jlm shallow convection 97 ***************
      if(ksc==97)then
!       this one is similar to ksc=99, but uses enhanced qbas <= qs
!       and also works up the column with other possible thebas values
        do k=kscbase,ksctop
         do iq=1,ifull
          qbas=min(qg(iq,k)/tied_rh,qs(iq,k))
c         thebas(iq,k)=t(iq,k)+hlcp*qbas  + hght
          thebas(iq,k)=prcpv(k)*t(iq,k)*(t(iq,k) + .5*hlcp*qbas)
     &                                 /(t(iq,k) - .5*hlcp*qbas)
         enddo  ! iq loop
        enddo   ! k loop
        theeb(:)=thebas(:,kscbase)
        do k=kscbase+1,ksctop+1
         do iq=1,ifull
          thee(iq,k)=prcpv(k)*t(iq,k)*(t(iq,k) + .5*hlcp*qs(iq,k))
     &                         /(t(iq,k) - .5*hlcp*qs(iq,k))
          if(theeb(iq)>thee(iq,k))then
            rk_shal(iq,k-1)=tied_con       !  m**2/sec  6., originally 10.
            rk_shal(iq,k)=tied_over        !  m**2/sec
          else
            theeb(iq)=thebas(iq,k) ! ready for next k in k-loop
          endif
         enddo  ! iq loop
        enddo   ! k loop
       if(ntest==3)then
         do iq=1,ifull
          if(rk_shal(iq,ksctop-1)>.9*tied_con)then 
            print *,'iq,land,rk_shal ',iq,land(iq),rk_shal(iq,ksctop-1)
          endif
          enddo
        endif   ! (ntest==3)
      endif     ! (ksc==97)
c     *********** end of Tiedtke_jlm shallow convection 97 *************

c     add in effects of shallow convection
      do k=1,kl-1  
       rkh(:,k)=rkh(:,k)+rk_shal(:,k)
      enddo   !  k loop
      if(kscmom==1)then
        do k=1,kl-1  
         rkm(:,k)=rkm(:,k)+rk_shal(:,k)
        enddo   !  k loop
      endif     ! (kscmom==1)
      
      if(ksc.ne.0.and.(ntest.ne.0.or.diag).and.nproc==1.)then
        do iq=1,ifull
         if(rk_shal(iq,1)>rk_shal(iq,2))then
           print *,'iq,rk_shal1,rk_shal2',iq,rk_shal(iq,1),rk_shal(iq,2)
         endif
        enddo
        if (mydiag) then 
        iq=idjd
        print *,'for shallow convection in vertmix '
        print *,'ktsav,ksctop ',ktsav(idjd),ksctop
        print *,'kbase,ktop ',kbase(idjd),ktop(idjd)
        print *,'kbsav,ktsav,theeb: ',kbsav(iq),ktsav(iq),theeb(iq)
        write (6,"('rk_shal ',9f7.2/(9x,9f7.2))")
     &             (rk_shal(idjd,k),k=1,kl)
        write (6,"('rh   ',9f7.2/(5x,9f7.2))") 
     &             (100.*qg(idjd,k)/qs(idjd,k),k=1,kl)
        write (6,"('qs   ',9f7.3/(5x,9f7.3))") 
     &             (1000.*qs(idjd,k),k=1,kl)
        write (6,"('qg   ',9f7.3/(5x,9f7.3))") 
     &             (1000.*qg(idjd,k),k=1,kl)
        write (6,"('qbas ',9f7.3/(5x,9f7.3))") 
     &             (1000.*qg(idjd,k)/tied_rh,k=1,kl)
        write (6,"('t    ',9f7.2/(5x,9f7.2))") 
     &             (t(idjd,k),k=1,kl)
        write (6,"('thebas',9f7.2/(6x,9f7.2))") 
     &             (thebas(iq,k),k=1,kl)
c        write (6,"('hs  ',9f7.2/(4x,9f7.2))") 
c     &             (t(idjd,k)+hlcp*qs(idjd,k),k=1,kl)
        write (6,"('thee',9f7.2/(4x,9f7.2))") 
     &             (thee(idjd,k),k=1,kl)
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
        call maxmin(rkh,'rk',ktau,.01,kl-1)
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

      if(ncvmix>0)then  ! cumulus mixing of momentum - jlm version
        do k=kuocb+1-ncvmix,kcl_top-1
         do iq=1,ifull
!          for no-conv-points, doing k=1-ncvmix,1 with convpsav=0.
!          1 below for ncvmix=2
           if(kbsav(iq)>0.and.k>=kbsav(iq)+1-ncvmix
     &                      .and.k<ktsav(iq))
     &       guv(iq,k)=guv(iq,k)-convpsav(iq)*.5      ! with factor of .5
          enddo  ! iq loop
        enddo    ! k loop
        if(diag.and.mydiag)then
          print *,'vertmix after conv; kb,kt,convp'
     &     ,kbsav(idjd),ktsav(idjd),convpsav(idjd)
          print *,'new guv',(guv(idjd,k),k=1,kl)
        endif
      endif      !   (ncvmix>0)

      conflux=grav*dt/dsig(1)
      condrag=grav*dt/(dsig(1)*rdry)
c     first do theta (then convert back to t)
      at(:,1)=0.
      rhs(1:ifull,1)=rhs(1:ifull,1)-(conflux/cp)*fg(1:ifull)/ps(1:ifull)

      if(npanels>0)then
        do k=2,kl
         do iq=1,ifull
          at(iq,k) =-gt(iq,k-1)/dsig(k)
         enddo   ! iq loop
        enddo    !  k loop
        do k=1,kl
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
         do iq=1,ifull
          ct(iq,k) =-gt(iq,k)/dsig(k)
         enddo   ! iq loop
        enddo    !  k loop
      endif      !  (npanels>0) .. else ..
      if((diag.or.ntest==2).and.mydiag)then
        print *,'ktau,fg,tss,ps ',ktau,fg(idjd),tss(idjd),ps(idjd)
        print *,'at ',(at(idjd,k),k=1,kl)
        print *,'ct ',(ct(idjd,k),k=1,kl)
        print *,'rhs ',(rhs(idjd,k),k=1,kl)
      endif      ! (ntest==2)

      if(nmaxpr==1.and.mydiag)
     &  write (6,"('thet_inx',9f8.3/8x,9f8.3)") rhs(idjd,:)
      call trim(at,ct,rhs,0)   ! for t
      if(nmaxpr==1.and.mydiag)
     &  write (6,"('thet_out',9f8.3/8x,9f8.3)") rhs(idjd,:)
      do k=1,kl
        do iq=1,ifull
         t(iq,k)=rhs(iq,k)/sigkap(k)
        enddo  ! iq loop
      enddo      !  k loop
      if(diag)then
         if (mydiag) then
            print *,'vertmix eg,fg,cdtq,land '
     &                  ,eg(idjd),fg(idjd),cdtq(idjd),land(idjd)
         end if
        call printa('thet',rhs,ktau,nlv,ia,ib,ja,jb,200.,1.)
      endif
      
      if (nvmix.eq.6) then ! MJT tke
        at=2.5*at          ! MJT tke
        ct=2.5*ct          ! MJT tke
      end if               ! MJT tke

c     now do moisture
      rhs(1:ifull,:)=qg(1:ifull,:)
      rhs(1:ifull,1)=rhs(1:ifull,1)-(conflux/hl)*eg(1:ifull)/ps(1:ifull)
c     could add extra sfce moisture flux term for crank-nicholson
      call trim(at,ct,rhs,0)    ! for qg
      qg(1:ifull,:)=rhs(1:ifull,:)
      if(diag.and.mydiag)then
        print *,'vertmix rhs & qg after trim ',(rhs(idjd,k),k=1,kl)
        write (6,"('qg ',9f7.3/(8x,9f7.3))") 
     &             (1000.*qg(idjd,k),k=1,kl)
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

c     now do trace gases rml 16/02/06 changed call
      if (ngas>0) call tracervmix(at,ct)
      
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
      if((diag.or.ntest==2).and.mydiag)then
        print *,'au ',(au(idjd,k),k=1,kl)
        print *,'cu ',(cu(idjd,k),k=1,kl)
      endif      ! (ntest==2)

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
 
      if((diag.or.ntest>=1).and.mydiag)then
        print *,'after trim in vertmix '
        write (6,"('thet',9f7.2/(8x,9f7.2))") 
     &              (sigkap(k)*t(idjd,k),k=1,kl) 
        write (6,"('t   ',9f7.2/(8x,9f7.2))") (t(idjd,k),k=1,kl) 
        write (6,"('u   ',9f7.2/(8x,9f7.2))") (u(idjd,k),k=1,kl) 
        write (6,"('v   ',9f7.2/(8x,9f7.2))") (v(idjd,k),k=1,kl) 
        print *,'cduv,cduv+1,cduvj+1,tss ',
     &           cduv(idjd),cduv(idjd+1),cduv(idjd+il),tss(idjd)
        write (6,"('au ',9f7.3/(8x,9f7.3))") (au(idjd,k),k=1,kl) 
        write (6,"('cu ',9f7.3/(8x,9f7.3))") (cu(idjd,k),k=1,kl) 
      endif
      if(nmaxpr==1.and.mydiag)then
        write (6,"('qg_vm ',9f7.3/(6x,9f7.3))") 
     &             (1000.*qg(idjd,k),k=1,kl)
      endif
      return
      end

!     rml delted tracervmix as modified version in trvmix.f
