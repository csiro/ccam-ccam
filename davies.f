      subroutine davies    ! for globpea - only large-scale available
      use cc_mpi, only : mydiag
      parameter(ntest=0)
c     see nbox in parm: nbox=9 for efficient 9x9 ls forcing
c     kbotdav: for original scheme put kbot=1; 4 for 18-level?
      include 'newmpar.h' ! il,jl,kl,ij
      include 'arrays.h' ! t,u,v,ps
      include 'dava.h' ! davt
      include 'davb.h' ! psls,qgg,tt,uu,vv
      include 'parm.h' ! kbotdav,nbox,nud_u,nud_v,nud_t,nud_p,nud_q,nud_hrs
      include 'sigs.h' ! sig

!     and new nud_p, nud_t, nud_q, nud_uv just off/on switches (0/1)
      if(nud_hrs.lt.0)then
!       large-scale style with davies-style weights (already scaled for nud_hrs)
!       N.B. nbd.lt.0 sets up special weights in indata.f
        if(nmaxpr.eq.1.and.ktau.lt.5.and.mydiag)then
          print *,'davies in  uu,vv,qgg(kl) ',
     &             uu(idjd,nlv),vv(idjd,nlv),qgg(idjd,kl)
          print *,'davies in  u,v,qg(kl) ',
     &           u(idjd,nlv),v(idjd,nlv),qg(idjd,kl)
        endif
        if(nud_p.ne.0)then
          do iq=1,ifull
           psl(iq)=psl(iq)+(psls(iq)-psl(iq))
     &                        *davt(iq)*dt/3600.
          enddo  ! iq loop
        endif  ! (nud_p.ne.0)
        if(nud_t.ne.0)then
          do k=kbotdav,kl
           do iq=1,ifull
            t(iq,k)=t(iq,k)+(tt(iq,k)-t(iq,k))
     &                         *davt(iq)*dt/3600.
           enddo  ! iq loop
          enddo   ! k loop
        endif     ! (nud_t.ne.0)
        if(nud_q.ne.0)then
          if(nud_q.lt.0)then
            kupper=kl/2     ! only nudge bottom half of atmos for qg
          else
            kupper=kl
          endif  ! (nud_q.lt.0)
          do k=kbotdav,kupper
           do iq=1,ifull
            qg(iq,k)=qg(iq,k)+(qgg(iq,k)-qg(iq,k))
     &                           *davt(iq)*dt/3600.
           enddo  ! iq loop
          enddo   ! k loop
        endif     ! (nud_q.ne.0)
        if(nud_uv.eq.1)then
          do k=kbotdav,kl
           do iq=1,ifull
            u(iq,k)=u(iq,k)+(uu(iq,k)-u(iq,k))
     &                         *davt(iq)*dt/3600.
            v(iq,k)=v(iq,k)+(vv(iq,k)-v(iq,k))
     &                         *davt(iq)*dt/3600.
           enddo  ! iq loop
          enddo   ! k loop
        endif     ! (nud_uv.eq.1)
        if(nud_uv.eq.2)then   ! speed option
          do k=kbotdav,kl
           do iq=1,ifull
	     speed=sqrt(u(iq,k)**2+v(iq,k)**2)
	     if(speed.gt.1.)then
  	       speeduu=sqrt(uu(iq,k)**2+vv(iq,k)**2)
		rat=speeduu/speed
c             u(iq,k)=u(iq,k)+(rat*u(iq,k)-u(iq,k))
c    .                           *davt(iq)*dt/3600.
              u(iq,k)=u(iq,k)+u(iq,k)*(rat-1.)
     &                           *davt(iq)*dt/3600.
              v(iq,k)=v(iq,k)+v(iq,k)*(rat-1.)
     &                           *davt(iq)*dt/3600.
	     endif ! (speed.gt.1.)
           enddo  ! iq loop
          enddo   ! k loop
        endif     ! (nud_uv.eq.2)
        if(nmaxpr.eq.1.and.ktau.lt.5.and.mydiag)then
          print *,'davies out u,v,qg(kl) ',
     &           u(idjd,nlv),v(idjd,nlv),qg(idjd,kl)
        endif
        return
      endif       ! (nud_hrs.lt.0)

c     following not usually done in globpe, because nud_hrs usually < 0
c     surface pressure
c     regular boundary forcing
c     apply large scale forcing to psl
      if(ntest.eq.1.and.mydiag)print *,'calling lsforc for psl: ',nud_p
      if(nud_p.ne.0)call lsforc(psl,psls,nud_hrs)

c-----------------------------------------------------------------------
c     now do 3-D fields of t,qg,u,v

      do k=kbotdav,kl
c      apply large scale forcing to t
      if(ntest.eq.1)print *,'calling lsforc for t;k,nud_ ',k,nud_t
       if(nud_t.ne.0)call lsforc(t(1,k),tt(1,k),nud_hrs)
c      apply large scale forcing to u and v
      if(ntest.eq.1)print *,'calling lsforc for u: ',nud_hrs
       if(nud_uv.ne.0)call lsforc(u(1,k),uu(1,k),nud_hrs)
      if(ntest.eq.1)print *,'calling lsforc for v: ',nud_hrs
       if(nud_uv.ne.0)call lsforc(v(1,k),vv(1,k),nud_hrs)
c      apply large scale forcing to qg - usually off (e.g. SARCS, 140-yr)
      if(ntest.eq.1)print *,'calling lsforc for qg: ',nud_hrs
       if(nud_q.ne.0)call lsforc(qg(1,k),qgg(1,k),nud_hrs)
      enddo

      return
      end
c=======================================================================
      subroutine davset
      implicit none
      include 'newmpar.h' ! il,jl,kl,ij
      include 'arrays.h' ! t,u,v,ps
      include 'davb.h'
      psls(1:ifull) = psl(1:ifull)
      tt(1:ifull,:) = t(1:ifull,:)
      qgg(1:ifull,:) = qg(1:ifull,:)
      uu(1:ifull,:) = u(1:ifull,:)
      vv(1:ifull,:) = v(1:ifull,:)
      return
      end
c=======================================================================
      subroutine lsforc(adat,bdat,nrelaxt)  ! globpea: ia=ja=1
      use cc_mpi, only : mydiag
      parameter (ntest=0)
c     only nbox=1 set up for globpea so far
c force lam towards boundary/large-scale data for filtered part only
c input: adat = unfiltered lam data
c        bdat = unfiltered boundary/large-scale data
c        dt = time step of model (sec)
c        relaxt = relaxation time (sec) to reduce diff to 1/e of orig.
c        nbox = npts avg. (both x and y direction)
c output: adat = lam data relaxed towards filtered boundary data

      include 'newmpar.h' ! il,jl,kl
      include 'dates.h'    !ktime,kdate,timer,timeg,mtimer
      include 'parm.h' ! dt
c     input variables
      real adat(ifull),bdat(ifull)
      coef=dt/(nrelaxt*3600.) ! so as to restore diff to 1/e
      if(nrelaxt.lt.0)then   ! linear variation over relaxt period
       relaxt=nrelaxt
       factor=-mod(timer-.001,-relaxt)/relaxt
       coef=abs(factor*coef)
       if(timer.lt.1..and.mydiag)print *,'timer,factor,coef ',timer,factor,coef
      endif   !  (nrelaxt.lt.0)
      if(ntest.eq.1.and.mydiag)then
        print *,'in davies; nbox,nrelaxt,coef: ',nbox,nrelaxt,coef
        print *,'adat_in,bdat: ',adat(idjd),bdat(idjd)
      end if

      if(nbox.eq.9)stop 'nbox=9 not ready yet'

      if(nbox.eq.1)then    ! may be useful for trace gas runs
        do iq=1,ifull
          adat(iq)=adat(iq)+coef*(bdat(iq)-adat(iq))
        enddo  ! iq loop
      endif    ! (nbox.eq.1)
      if(ntest.eq.1.and.mydiag)print *,'adat_out: ',adat(idjd)

      return
      end
c=======================================================================
