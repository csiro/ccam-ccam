      subroutine adjust5
      parameter (moistfix=2) ! 0 earlier; 1 for cube-root fix; 2 best, with ps
      parameter (mfix_rad=0) ! used to make gases 2 to ng add up to gas 1
      parameter (nys=0)      ! for nys=1 use sun & yeh 1997 monotonic filter
!                              within adjust5 on u,v
!     nuv in parm.h:  now always 10
!          !  0y for increments   , 1y for actual values, 2y for true increments
      parameter (hl=2.5104e6)
      include 'newmpar.h'
      include 'arrays.h'
      include 'constant.h'
      include 'indices.h'
      include 'map.h'
      include 'morepbl.h'  ! condx,eg
      include 'nlin.h'
      include 'parm.h'     ! qgmin
      include 'parmdyn.h'  
      include 'parmvert.h'  
      include 'pbl.h'
      include 'sigs.h'
      include 'tracers.h'
      include 'vecs.h'
      include 'vecsuv.h'   ! vecsuv info
      include 'vecsuva.h'  ! vecsuva info
      include 'vvel.h'     ! sdot
      include 'xarrs.h'
      common/dpsdt/dpsdt(ifull)    ! shared adjust5 & openhist
      common/epst/epst(ifull)
      common/mfixdiag/alph_p,alph_pm,delneg,delpos,alph_q
      common/tbar2d/tbar2d(ifull)
      common/work3/p(ifull,kl),omgfnl(ifull,kl),
     . wcu(ifull),wcv(ifull),wcua(ifull),wcva(ifull),
     . wcud(ifull),wcvd(ifull),dum3(3*ijk - 6*ifull)
      common/work3b/wrk1(ifull,kl),wrk2(ifull,kl)   ! just work arrays here
      common/work3sav/qgsav(ifull,kl),trsav(ilt*jlt,klt,ngasmax) ! passed from nonlin
      common/work2/zz(ifull),zzn(ifull),zze(ifull),zzw(ifull),
     . zzs(ifull),aa(ifull),bb(ifull),cc(ifull),dd(ifull),
     . pslxint(ifull),pfact(ifull),alff(ifull),alf(ifull),alfe(ifull),
     . alfn(ifull),pslsav(ifull),alfu(ifull),alfv(ifull)
      dimension pe(ifull,kl),e(ifull,kl)
      real helm(ifull,kl),rhsl(ifull,kl),delps(ifull)
      dimension indsw(8),indne(8),indse(4),indnw(4)
      real d(ifull,kl),pextras(ifull,kl),omgf(ifull,kl)
      equivalence (d,vn),(omgf,pe,tn),(e,p),(pextras,dpsldt)
!      equivalence (cc,rhsl),(dd,helm),(bb,delps)
      data indsw/1,3,5,6,8,9,12,13/,indne/0,3,4,6,7,9,11,13/
      data indse/1,4,7,12/,indnw/0,5,8,11/
      save indsw,indne,indse,indnw
      ! If the time step hasn't changed, used the same acceleration as 
      ! previously
      ind(i,j,n)=i+(j-1)*il+n*il*il
      hdt=dt/2.
      hdtds=hdt/ds

      if(m.lt.6)then
        do iq=1,ifull
         alf(iq)=1.+epsu
         alff(iq)=0.
         alfu(iq)=1.   ! i.e. alf/(1+epsu)
         alfv(iq)=1.   ! i.e. alf/(1+epsu)
        enddo     ! iq loop
      else
        do iq=1,ifull
         alf(iq)=(1.+epsu)/(1.+(hdt*(1.+epsf)*f(iq))**2)
         alfF(iq)=alf(iq)*f(iq)*hdt*(1.+epsf)      ! now includes hdt in alfF
         alfu(iq)=1./(1.+(hdt*(1.+epsf)*fu(iq))**2) ! i.e. alf/(1+epsu)
         alfv(iq)=1./(1.+(hdt*(1.+epsf)*fv(iq))**2) ! i.e. alf/(1+epsu)
        enddo     ! iq loop
      endif  ! (m.lt.6)... else ...
!cdir nodep
      do iq=1,ifull
       pfact(iq)=4.*( ds/(dt*em(iq)) )**2    ! for adjust9/5
       alfe(iq)=alf(ie(iq))-alf(iq)
     .     +.25*(alff(ine(iq))+alfF(in(iq))-alfF(ise(iq))-alfF(is(iq)))
       alfn(iq)=alf(in(iq))-alf(iq)
     .     -.25*(alfF(ien(iq))+alfF(ie(iq))-alfF(iwn(iq))-alfF(iw(iq)))
      enddo     ! iq loop

!cdir nodep
      do iq=1,ifull
!      need care with vector quantities on w (odd) & s (even) panel boundaries
       zz(iq)=.5*(alfe(iwu2(iq))-alfe(iq)+alfn(isv2(iq))-alfn(iq))
     .                     -4.*alf(iq)                  ! i,j   coeff
       zzn(iq)=alf(in(iq))-.5*alfn(iq)                  ! i,j+1 coeff
       zzw(iq)=alf(iw(iq))+.5*alfe(iwu2(iq))            ! i-1,j coeff
       zze(iq)=alf(ie(iq))-.5*alfe(iq)                  ! i+1,j coeff
       zzs(iq)=alf(is(iq))+.5*alfn(isv2(iq))            ! i,j-1 coeff
      enddo     ! iq loop
!     N.B. there are some special z values at the 8 vertices
      if(npanels.eq.5)then
        do n=0,5
         iq=ind(1,1,n)
         zzs(iq)=zzs(iq)+.25*alfF(is(iq))   ! i,j-1 coeff
         zzw(iq)=zzw(iq)-.25*alfF(iw(iq))   ! i-1,j coeff
         iq=ind(il,il,n)
         zzn(iq)=zzn(iq)+.25*alfF(in(iq))   ! i,j+1 coeff
         zze(iq)=zze(iq)-.25*alfF(ie(iq))   ! i+1,j coeff
         iq=ind(il,1,n)
         zzs(iq)=zzs(iq)-.25*alfF(is(iq))   ! i,j-1 coeff
         zze(iq)=zze(iq)+.25*alfF(ie(iq))   ! i+1,j coeff
         iq=ind(1,il,n)
         zzn(iq)=zzn(iq)-.25*alfF(in(iq))   ! i,j+1 coeff
         zzw(iq)=zzw(iq)+.25*alfF(iw(iq))   ! i-1,j coeff
        enddo   ! n loop
      endif     ! (npanels.eq.5)
      if(npanels.eq.13)then
        do nn=1,8
         iq=ind(1,1,indsw(nn))
         zzs(iq)=zzs(iq)+.25*alfF(is(iq))   ! i,j-1 coeff
         zzw(iq)=zzw(iq)-.25*alfF(iw(iq))   ! i-1,j coeff
         iq=ind(il,il,indne(nn))
         zzn(iq)=zzn(iq)+.25*alfF(in(iq))   ! i,j+1 coeff
         zze(iq)=zze(iq)-.25*alfF(ie(iq))   ! i+1,j coeff
        enddo   ! nn loop
        do nn=1,4
         iq=ind(il,1,indse(nn))
         zzs(iq)=zzs(iq)-.25*alfF(is(iq))   ! i,j-1 coeff
         zze(iq)=zze(iq)+.25*alfF(ie(iq))   ! i+1,j coeff
         iq=ind(1,il,indnw(nn))
         zzn(iq)=zzn(iq)-.25*alfF(in(iq))   ! i,j+1 coeff
         zzw(iq)=zzw(iq)+.25*alfF(iw(iq))   ! i-1,j coeff
        enddo   ! nn loop
      endif     ! (npanels.eq.13)
      if(diag)then
        print *,'entering adjust5'
        iq=id+il*(jd-1)
        print *,'tx ',(tx(idjd,k),k=1,kl)
        call printa('tx  ',tx(1,nlv),ktau,nlv,ia,ib,ja,jb,200.,1.)
        print *,'ux_stag ',(ux(idjd,k),k=1,kl)
        print *,'vx_stag ',(vx(idjd,k),k=1,kl)
        write (6,"('qg ',19f7.3/(8x,19f7.3))") 
     .             (1000.*qg(idjd,k),k=1,kl)
c       print *,'alf & n e',alf(iq),alf(in(iq)),alf(ie(iq))
c       call printa('alf ',alf,ktau,0,ia,ib,ja,jb,0.,100.)
c       print *,'alfF & n e w s',alfF(iq),alfF(in(iq)),alfF(ie(iq))
c    .                            ,alfF(iw(iq)),alfF(is(iq))
c       print *,'alfF: ne se en wn',alfF(ine(iq)),alfF(ise(iq))
c    .                             ,alfF(ien(iq)),alfF(iwn(iq))
c       print *,'alfe & n e w s',alfe(iq),alfe(in(iq)),alfe(ie(iq))
c    .                          ,alfe(iw(iq)),alfe(is(iq))
c       call printa('alfe',alfe,ktau,0,ia,ib,ja,jb,0.,1000.)
c       print *,'alfn & n e w s',alfn(iq),alfn(in(iq)),alfn(ie(iq))
c    .                          ,alfn(iw(iq)),alfn(is(iq))
c       call printa('alfn',alfn,ktau,0,ia,ib,ja,jb,0.,1000.)
c       print *,'isv2,alfn(isv2),zzs ',isv2(iq),alfn(isv2(iq)),zzs(iq)
        call maxmin(alf,'a ',ktau,1.,1)
        call maxmin(alfe,'ae',ktau,1.,1)
        call maxmin(alfn,'an',ktau,1.,1)
        print  *,'pslx ',(pslx(idjd,k),k=1,kl)
        print  *,'pslx w ',(pslx(iw(idjd),k),k=1,kl)
        call printa('pslx',pslx(1,nlv),ktau,nlv,ia,ib,ja,jb,0.,100.)
      endif

!     recompute nonlinear sigma-dot contribution for updating tn, tx
!     e contains intgrl{-dpsldt x dsigma} from 0 to sigma (dsig is -ve)
!     vert. integ. nonlin part of mass weighted div into e
      do k=1,kl
       do iq=1,ifull
        pslx(iq,k)=pslx(iq,k)*2./dt  ! i.e. [RHS of (2.5)]*2/dt, i.e. M
       enddo   ! iq loop
      enddo    ! k  loop
      do iq=1,ifull
       e(iq,kl)=dsig(kl)*pslx(iq,kl)
      enddo     ! iq loop
      do k=kl-1,1,-1
       do iq=1,ifull
        e(iq,k)=e(iq,k+1)+dsig(k)*pslx(iq,k)  ! i.e. M_bar_sig-.5
       enddo    ! iq loop
      enddo     ! k loop
      do iq=1,ifull
       pslxint(iq)=-e(iq,1)*dt/2. ! pslxint holds integrated pslx (2.6)
      enddo     ! iq loop

!     full-level (1.+epsp)*omega/ps into omgfnl (nonlinear part only) (2.8)
      do k=1,kl-1
       do iq=1,ifull
        omgfnl(iq,k)=-rata(k)*e(iq,k+1)-ratb(k)*e(iq,k)
     .                 -sig(k)*pslx(iq,k)
       enddo    ! iq loop
      enddo     ! k loop
      do iq=1,ifull
       omgfnl(iq,kl)=-ratb(kl)*e(iq,kl)-sig(kl)*pslx(iq,kl)
      enddo     ! iq loop

      do k=1,kl
       do iq=1,ifull
!       N.B. the omgfnl term on LHS of (2.4) not yet added in
        tx(iq,k)=tx(iq,k)+hdt*                         ! (2.9)
     .             tbar2d(iq)*omgfnl(iq,k)*roncp/sig(k)  ! with correct epsp
       enddo    ! iq loop
      enddo     ! k loop

!     calculate heights from the tx array
      do iq=1,ifull
       p(iq,1)=zs(iq)+bet(1)*tx(iq,1)+r*tbar2d(iq)*pslxint(iq) ! (2.21)
      enddo     ! iq loop
      do k=2,kl
       do iq=1,ifull
        p(iq,k)=p(iq,k-1)+bet(k)*tx(iq,k)+betm(k)*tx(iq,k-1)
       enddo    ! iq loop
      enddo     ! k loop

      do k=1,kl
       do iq=1,ifull
        p(iq,k)=p(iq,k)+pextras(iq,k)   ! npex=1 code
       enddo   ! iq loop
      enddo    ! k  loop

!     redefine ux, vx
      do k=1,kl
       do iq=1,ifull
        cc(iq)=ux(iq,k)/emu(iq)*alfu(iq)
        dd(iq)=vx(iq,k)/emv(iq)*alfv(iq)
       enddo    ! iq loop

!      form divergence of rhs (ux & vx) terms: xd
!cdir nodep
       do iq=1,ifull
         d(iq,k)=cc(iq)-cc(iwu2(iq))+dd(iq)-dd(isv2(iq))
       enddo    ! iq loop
      enddo     ! k  loop

!     transform p & d to eigenvector space
      do k=1,kl
         do iq=1,ifull
            pe(iq,k)=einv(k,1)*p(iq,1)
            rhsl(iq,k)=einv(k,1)*d(iq,1)
         enddo
         do l=2,kl
            do iq=1,ifull
               pe(iq,k)=pe(iq,k)+einv(k,l)*p(iq,l) ! xp in eig space
               rhsl(iq,k)=rhsl(iq,k)+einv(k,l)*d(iq,l) ! xd in eig space
            enddo
         enddo                  ! l loop

         do iq=1,ifull
            helm(iq,k) = pfact(iq)*tbar(1)/
     &                   (bam(k)*(1.+epst(iq))*tbar2d(iq))
            rhsl(iq,k) = rhsl(iq,k)/hdtds -helm(iq,k)*pe(iq,k)
         enddo                  ! iq loop
         if(diag.and.k.le.2)then !  only for last k of loop (i.e. 1)
            iq=idjd
            print  *,'adjust5(k) p & n e w s ',k,p(iq,k),
     &        p(in(iq),k),p(ie(iq),k),p(iw(iq),k),p(is(iq),k)
            print  *,'adjust5(k) pe & n e w s ',k,pe(iq,k),
     &        pe(in(iq),k),pe(ie(iq),k),pe(iw(iq),k),pe(is(iq),k)
            print  *,'adjust5(k) rhsl & n e w s ',k,rhsl(iq,k),
     &       rhsl(in(iq),k),rhsl(ie(iq),k),rhsl(iw(iq),k),rhsl(is(iq),k)
         endif                  ! (diag.and.k.le.2)
      enddo    ! k loop

      call helmsol(helm,pe,rhsl)

      if(diag)then   !  only for last k of loop (i.e. 1)
        iq=idjd
        print  *,'adjust5(1) pslxint & n e w s',pslxint(iq),
     .   pslxint(in(iq)),pslxint(ie(iq)),pslxint(iw(iq)),
     .   pslxint(is(iq))
        call printa('psnt',pslxint,ktau,0,ia,ib,ja,jb,0.,100.)
        print *,'adjust5 tx ',(tx(idjd,k),k=1,kl)
        call printa('tx  ',tx(1,nlv),ktau,nlv,ia,ib,ja,jb,200.,1.)
        print  *,'adjust5(1) helm & n e w s',helm(iq,k),
     .   helm(in(iq),k),helm(ie(iq),k),helm(iw(iq),k),helm(is(iq),k)
        print  *,'adjust5(1) rhsl & n e w s',rhsl(iq,k),
     .   rhsl(in(iq),k),rhsl(ie(iq),k),rhsl(iw(iq),k),rhsl(is(iq),k)
        call printa('rhsl',rhsl(1,1),ktau,1,ia,ib,ja,jb,0.,0.)
        print  *,'adjust5 pe ',(pe(iq,k),k=1,kl)
        print  *,'adjust5 pe e ',(pe(ie(iq),k),k=1,kl)
        print  *,'adjust5 pe w ',(pe(iw(iq),k),k=1,kl)
        call printa('pe  ',pe(1,1),ktau,1,ia,ib,ja,jb,0.,0.)
      endif

!     only use method II inversion nowadays
      do k=1,kl
!      first p from pe
       do iq=1,ifull
        p(iq,k)=emat(k,1)*pe(iq,1)
       enddo    ! iq loop
       do l=2,kl      ! this is the remaining expensive loop
        do iq=1,ifull
         p(iq,k)=p(iq,k)+emat(k,l)*pe(iq,l)
        enddo   ! iq loop
       enddo    !  l loop

!      now u & v
       if(k.eq.nlv.and.diag)then
         iq=idjd
         print *,'iq,k,fu,alfu,alfu*ux(iq,k) ',
     .    iq,k,fu(iq),alfu(iq),alfu(iq)*ux(iq,k)
         print *,'alfF & n e w s (in(iq)),alfF(ine(iq)),alfF(is(iq))',
     .           'alfF(ise(iq)),alfe(iq) ',alfF(in(iq)),alfF(ine(iq))
     .           ,alfF(is(iq)),alfF(ise(iq)),alfe(iq)
         sum=alf(ie(iq))-alf(iq)+.25*(alfF(in(iq))+
     .        alfF(ine(iq))-alfF(is(iq))-alfF(ise(iq))) -alfe(iq)
         print *,'sum  ',sum
         print *,'p & n e w s ne se ',p(iq,k),p(in(iq),k),
     .    p(ie(iq),k),p(iw(iq),k),p(is(iq),k),
     .    p(ine(iq),k),p(ise(iq),k)
       endif
!cdir nodep
       do iq=1,ifull
c       u(iq,k)=alfu(iq)*ux(iq,k)
        cc(iq) =alfu(iq)*ux(iq,k)   ! globpea
     .   -hdtds*emu(iq)*
     .     ( alf(ie(iq))*p(ie(iq),k)-alf(iq)*p(iq,k)
     .   -.5*alfe(iq)*(p(iq,k)+p(ie(iq),k))
     .   +.25*(alfF(in(iq))*p(in(iq),k) +alfF(ine(iq))*p(ine(iq),k)
     .      -alfF(is(iq))*p(is(iq),k) -alfF(ise(iq))*p(ise(iq),k)) )
c       v(iq,k)=alfv(iq)*vx(iq,k)
        dd(iq) =alfv(iq)*vx(iq,k)   ! globpea
     .   -hdtds*emv(iq)*
     .     ( alf(in(iq))*p(in(iq),k)-alf(iq)*p(iq,k)
     .   -.5*alfn(iq)*(p(iq,k)+p(in(iq),k))
     .   -.25*(alfF(ien(iq))*p(ien(iq),k) +alfF(ie(iq))*p(ie(iq),k)
     .      -alfF(iwn(iq))*p(iwn(iq),k) -alfF(iw(iq))*p(iw(iq),k)) )
       enddo    ! iq loop

!      calculate linear part only of sigma-dot and omega/ps
!cdir nodep
       do iq=1,ifull
        d(iq,k)=(cc(iq)/emu(iq)-cc(iwu2(iq))/emu(iwu2(iq))  ! globpea
     .          +dd(iq)/emv(iq)-dd(isv2(iq))/emv(isv2(iq))) ! globpea
     .          *em(iq)**2/ds
       enddo    ! iq loop

!      straightforward rev. cubic interp of u and v (i.e. nuv=10)
       call unstaguv(cc,dd,u(1,k),v(1,k))  ! usual
      enddo     ! k  loop ---------------------------------------------------------

      if(diag.or.nmaxpr.eq.1)then
        write (6,"('div_adj(1-9) ',9f8.2)") (d(idjd,k)*1.e6,k=1,9)
        write (6,"('omgfnl*dt    ',9f8.4)") (omgfnl(idjd,k)*dt,k=1,9)
      endif

!     vert. integ. div into e
      do iq=1,ifull
       e(iq,kl)=-dsig(kl)*d(iq,kl)
      enddo     ! iq loop
      do k=kl-1,1,-1
       do iq=1,ifull
        e(iq,k)=e(iq,k+1)-dsig(k)*d(iq,k)
       enddo    ! iq loop
      enddo     ! k  loop

!     full-level omega/ps into omgf (linear part only)
      do k=1,kl-1
       do iq=1,ifull
        omgf(iq,k)=-rata(k)*e(iq,k+1)-ratb(k)*e(iq,k)
       enddo    ! iq loop
      enddo     ! k  loop

      do iq=1,ifull
       omgf(iq,kl)=-ratb(kl)*e(iq,kl)
       pslsav(iq)=psl(iq)   ! in work2 (not last because of vadvtvd)
       psl(iq)=pslxint(iq)-hdt*e(iq,1)  *(1.+epst(iq))
      enddo     ! iq loop

      do k=1,kl
       do iq=1,ifull
!       save [D + dsigdot/dsig] in pslx for next use in nonlin
        pslx(iq,k)=(pslx(iq,k)-psl(iq)*2./dt)/(1.+epst(iq))
       enddo  !  iq loop
      enddo   !  k loop
      do k=kl,2,-1
       do iq=1,ifull
!       calculate latest sdot (at level k-.5)
        sdot(iq,k)=sdot(iq,k+1)-dsig(k)*(pslx(iq,k)-d(iq,k))
       enddo  !  iq loop
      enddo   !  k loop
      do k=2,kl
       do iq=1,ifull
!       and convert sdot (at level k-.5) to units of grid-steps/timestep
!       dtin is used to be ready for next full timestep
        sdot(iq,k)=sdot(iq,k)*dtin/(sig(k)-sig(k-1))
       enddo   ! iq loop
      enddo    ! k  loop

      do k=1,kl
       do iq=1,ifull
!       save full omega/ps in dpsldt for use in nonlin next time step (& outfile)
!       N.B. omgfnl part already incorp. into tx above
        dpsldt(iq,k)=omgfnl(iq,k)/(1.+epst(iq))+omgf(iq,k)
        t(iq,k)=tx(iq,k)
     .           +hdt*(1.+epst(iq))*tbar2d(iq)*omgf(iq,k)*roncp/sig(k)
       enddo    ! iq loop
      enddo     ! k  loop

      if(nvadh.eq.2)then                 ! final dt/2 's worth
        if(diag.or.nmaxpr.eq.1)then
         print *,'before vertical advection in adjust5'
         write (6,"('sdot',9f8.3/4x,9f8.3)") (sdot(idjd,kk),kk=1,kl)
         write (6,"('t   ',9f8.2/4x,9f8.2)") (t(idjd,kk),kk=1,kl)
         write (6,"('u   ',9f8.2/4x,9f8.2)") (u(idjd,kk),kk=1,kl)
         write (6,"('v   ',9f8.2/4x,9f8.2)") (v(idjd,kk),kk=1,kl)
         write (6,"('qg  ',9f8.3/4x,9f8.3)")(1000.*qg(idjd,kk),kk=1,kl)
        endif
        if(nvad.eq.4)call vadvtvd(t,t,u,u,v,v)
        if(nvad.eq.7)call vadv30(t,t,u,u,v,v)   ! for vadvbess
        if(diag.or.nmaxpr.eq.1)then
          print *,'after vertical advection in adjust5'
          write (6,"('qg  ',9f8.3/4x,9f8.3)")(1000.*qg(idjd,kk),kk=1,kl)
          write (6,"('t   ',9f8.2/4x,9f8.2)") (t(idjd,kk),kk=1,kl)
          write (6,"('thet',9f8.2/4x,9f8.2)")  
     .                 (t(idjd,k)*sig(k)**(-roncp),k=1,kl)
        write (6,"('u   ',9f8.2/4x,9f8.2)") (u(idjd,kk),kk=1,kl)
        write (6,"('v   ',9f8.2/4x,9f8.2)") (v(idjd,kk),kk=1,kl)
        endif
      endif     !  (nvadh.eq.2)then

!     could call nestin here (using xarrs instead of davb.h)
      if(mspec.eq.1.and.nbd.ne.0)call davies  ! before mfix mass fix in C-C
      if(mfix.gt.0)then   ! perform conservation fix on psl
!       delpos is the sum of all positive changes over globe
!       delneg is the sum of all negative changes over globe
!       alph_p is chosen to satisfy alph_p*delpos + delneg/alph_p = 0
        delpos=0.
        delneg=0.
        do iq=1,ifull
         delps(iq)=psl(iq)-pslsav(iq)
         delpos=delpos+max(0.,delps(iq)/em(iq)**2)
         delneg=delneg+min(0.,delps(iq)/em(iq)**2)
        enddo   ! iq loop
        if(mfix.eq.1)then
          alph_p = sqrt( -delneg/delpos)
          alph_pm=1./alph_p
        endif   ! (mfix.eq.1)
        if(mfix.eq.2)then
          if(delpos.gt.-delneg)then
            alph_p =1.
            alph_pm=-delpos/delneg
          else
            alph_p=-delneg/delpos
            alph_pm =1.
          endif
        endif   ! (mfix.eq.2)
        do iq=1,ifull
         psl(iq)=pslsav(iq)+
     .         alph_p*max(0.,delps(iq)) + alph_pm*min(0.,delps(iq))
        enddo   ! iq loop
      endif     !  mfix.gt.0

      do iq=1,ifull
       aa(iq)=ps(iq)  ! saved for gas fixers below, and other diags
       ps(iq)=1.e5*exp(psl(iq))
       dpsdt(iq)=(ps(iq)-aa(iq))*24.*3600./(100.*dt) ! diagnostic in hPa/day
      enddo     !  iq loop

      if(mfix_qg.gt.0.and.mspec.eq.1)then
!       qgmin=1.e-6   !  in parm.h 
        qgminm=qgmin
        if(moistfix.eq.1)then  !  cube-root method
          qgminm=qgmin**(1./3.)
          do k=1,kl
	    do iq=1,ifull
!           qg(iq,k)=cbrt(qg(iq,k))
            qg(iq,k)=     qg(iq,k)**(1./3.)
!           qgsav(iq,k)=cbrt(qgsav(iq,k))
            qgsav(iq,k)=     qgsav(iq,k)**(1./3.)
           enddo   ! iq loop
          enddo    ! k  loop
        endif      ! (moistfix.eq.1)
!  	 default now with ps weighting
        do k=1,kl
	  do iq=1,ifull
          qg(iq,k)=qg(iq,k)*ps(iq)
          qgsav(iq,k)=qgsav(iq,k)*aa(iq)
         enddo   ! iq loop
        enddo    ! k  loop
!       perform conservation fix on qg, as affected by vadv, hadv, hordif
!       N.B. won't cope with any -ves from conjob
!       delpos is the sum of all positive changes over globe
!       delneg is the sum of all negative changes over globe
        delpos=0.
        delneg=0.
c       print *,'qgsav,qg_in',qgsav(idjd,1),qg(idjd,1,1)
        do iq=1,ifull
         do k=1,kl
          wrk1(iq,k)=max(qg(iq,k),qgminm*ps(iq))-qgsav(iq,k)  ! increments
          delpos=delpos+max(0.,-dsig(k)*wrk1(iq,k)/em(iq)**2)
          delneg=delneg+min(0.,-dsig(k)*wrk1(iq,k)/em(iq)**2)
         enddo   ! k loop
        enddo    ! iq loop
        ratio = -delneg/delpos
        if(mfix_qg.eq.1)alph_q = min(ratio,sqrt(ratio))  ! why min?
        if(mfix_qg.eq.2)alph_q = sqrt(ratio)
        do k=1,kl        ! this is cunning 2-sided scheme
         do iq=1,ifull
         qg(iq,k)=qgsav(iq,k)+
     .     alph_q*max(0.,wrk1(iq,k)) + min(0.,wrk1(iq,k))/max(1.,alph_q)
         enddo   ! iq loop
        enddo    ! k  loop
        if(moistfix.eq.1)then
         do k=1,kl
          do iq=1,ifull
           qg(iq,k)=qg(iq,k)**3
          enddo   ! iq loop
         enddo    ! k  loop
        endif      ! (moistfix.eq.1)
!  	 undo ps weighting
        do k=1,kl
	  do iq=1,ifull
          qg(iq,k)=qg(iq,k)/ps(iq)
         enddo   ! iq loop
        enddo    ! k  loop
	 delpos=delpos/1.e5
 	 delneg=delneg/1.e5  ! for diag print in globpe
      endif        !  (mfix_qg.gt.0)

      if(mfix_qg.gt.0.and.mspec.eq.1.and.ngas.ge.1)then
!       perform conservation fix on tr1,tr2 as affected by vadv, hadv, hordif
        do ng=1,ngas
!  	 default now with ps weighting
        do k=1,kl
	   do iq=1,ifull
           tr(iq,k,ng)=tr(iq,k,ng)*ps(iq)
           trsav(iq,k,ng)=trsav(iq,k,ng)*aa(iq)
          enddo   ! iq loop
         enddo    ! k  loop
         delpos=0.
         delneg=0.
         do iq=1,ifull
          do k=1,kl
           wrk1(iq,k)=max(tr(iq,k,ng),gasmin(ng)*ps(iq))
     .                                         -trsav(iq,k,ng)  ! has increments
           delpos=delpos+max(0.,-dsig(k)*wrk1(iq,k)/em(iq)**2)
           delneg=delneg+min(0.,-dsig(k)*wrk1(iq,k)/em(iq)**2)
          enddo   ! k loop
         enddo    ! iq loop
         ratio = -delneg/delpos
         if(mfix_qg.eq.1)alph_g = min(ratio,sqrt(ratio))  ! why min?
         if(mfix_qg.eq.2)alph_g = sqrt(ratio)
         do k=1,kl   ! this is cunning 2-sided scheme
	   do iq=1,ifull
          tr(iq,k,ng)=trsav(iq,k,ng)+
     .     alph_g*max(0.,wrk1(iq,k)) + min(0.,wrk1(iq,k))/max(1.,alph_g)
          enddo   ! iq loop
         enddo    ! k  loop
         if(moistfix.eq.2)then  !  with ps weighting
           do k=1,kl
	     do iq=1,ifull
             tr(iq,k,ng)=tr(iq,k,ng)/ps(iq)
             trsav(iq,k,ng)=trsav(iq,k,ng)/ps(iq)
            enddo   ! iq loop
           enddo    ! k  loop
         endif      ! (moistfix.eq.2)
        enddo    ! ng loop
        if(mfix_rad.gt.0)then  ! to make gases 2 to ng add up to gas 1
         do k=1,kl
	   do iq=1,ifull
           sumdiffb=0.
           delpos=0.
           delneg=0.
           do ng=2,ngas        
            sumdiffb=sumdiffb+tr(iq,k,ng)
            delpos=delpos+max( 1.e-20,tr(iq,k,ng)-trsav(iq,k,ng))
            delneg=delneg+min(-1.e-20,tr(iq,k,ng)-trsav(iq,k,ng))
           enddo   ! ng loop
           ratio=(tr(iq,k,1)-sumdiffb)/(delpos-delneg)
           do ng=2,ngas        
            tr(iq,k,ng)=max(0.,trsav(iq,k,ng)
     .         +(1.+ratio)*max(0.,tr(iq,k,ng)-trsav(iq,k,ng))
     .         +(1.-ratio)*min(0.,tr(iq,k,ng)-trsav(iq,k,ng)) )
           enddo   ! ng loop
          enddo   ! iq loop
         enddo    ! k  loop
        endif     ! (mfix_rad.gt.0)
      endif       !  mfix_qg.gt.0

      if(diag)then
        iq=idjd
        print  *,'adjust5 d ',(d(iq,k),k=1,kl)
        print  *,'adjust5 d: e ',(d(ie(iq),k),k=1,kl)
        print  *,'adjust5 d: w ',(d(iw(iq),k),k=1,kl)
        print  *,'adjust5 idjd,e ',idjd,(e(idjd,k),k=1,kl)
        print  *,'adjust5 psl & n e w s',psl(iq),
     .   psl(in(iq)),psl(ie(iq)),psl(iw(iq)),psl(is(iq))
        call printa('psl ',psl,ktau,0,ia,ib,ja,jb,0.,100.)
        print *,'adjust5 t ',(t(idjd,k),k=1,kl)
        call printa('t   ',t(1,nlv),ktau,nlv,ia,ib,ja,jb,200.,1.)
        print *,'adjust5 u ',(u(idjd,k),k=1,kl)
        call printa('u   ',u(1,nlv),ktau,nlv,ia,ib,ja,jb,0.,1.)
        print *,'adjust5 v ',(v(idjd,k),k=1,kl)
        call printa('v   ',v(1,nlv),ktau,nlv,ia,ib,ja,jb,0.,1.)
        print *,'ps diagnostic from end of adjust:'
        do iq=1,ifull
         aa(iq)=ps(iq)-aa(iq)
        enddo
        call printa('dps ',aa,ktau,0,ia,ib,ja,jb,0.,.01)
        call printa('ps  ',ps,ktau,0,ia,ib,ja,jb,1.e5,.01)
        write (6,"('qg ',19f7.3/(8x,19f7.3))") 
     .             (1000.*qg(idjd,k),k=1,kl)
        if(sig(nlv).lt..3)then
          call printa('qg  ',qg(1,nlv),ktau,nlv,ia,ib,ja,jb,0.,1.e6)
        else
          call printa('qg  ',qg(1,nlv),ktau,nlv,ia,ib,ja,jb,0.,1.e3)
        endif
      endif
      return
      end
