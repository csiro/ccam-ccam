      subroutine davies    ! for globpea - only large-scale available
      use cc_mpi, only : mydiag
      parameter(ntest=0)
c     from Nov 05, separate davu & davt, just abs(nud_hrs) used)
      include 'newmpar.h' ! il,jl,kl,ij
      include 'arrays.h' ! t,u,v,ps
      include 'dava.h' ! davt, davu
      include 'davb.h' ! psls,qgg,tt,uu,vv
      include 'parm.h' ! kbotdav,nud_u,nud_v,nud_t,nud_p,nud_q,nud_hrs,nudu_hrs
      include 'sigs.h' ! sig

!     and new nud_p, nud_t, nud_q, nud_uv just off/on switches (0/1)
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
            kupper=kl/2 ! option to only nudge bottom half of atmos for qg
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
          do k=kbotu,kl
           do iq=1,ifull
            u(iq,k)=u(iq,k)+(uu(iq,k)-u(iq,k))
     &                         *davu(iq)*dt/3600.
            v(iq,k)=v(iq,k)+(vv(iq,k)-v(iq,k))
     &                         *davu(iq)*dt/3600.
           enddo  ! iq loop
          enddo   ! k loop
          if(kbotdav<kbotu)then  ! Feb 06
            do k=kbotdav,kbotu-1
             do iq=1,ifull
              u(iq,k)=u(iq,k)+(uu(iq,k)-u(iq,k))
     &                           *davt(iq)*dt/3600.
              v(iq,k)=v(iq,k)+(vv(iq,k)-v(iq,k))
     &                           *davt(iq)*dt/3600.
             enddo  ! iq loop
            enddo   ! k loop
          endif
        endif     ! (nud_uv.eq.1)
        if(nud_uv.eq.2)then   ! speed option (not useful)
          do k=kbotu,kl
           do iq=1,ifull
	     speed=sqrt(u(iq,k)**2+v(iq,k)**2)
	     if(speed.gt.1.)then
  	       speeduu=sqrt(uu(iq,k)**2+vv(iq,k)**2)
		rat=speeduu/speed
c             u(iq,k)=u(iq,k)+(rat*u(iq,k)-u(iq,k))
c    .                           *davt(iq)*dt/3600.
              u(iq,k)=u(iq,k)+u(iq,k)*(rat-1.)
     &                           *davu(iq)*dt/3600.
              v(iq,k)=v(iq,k)+v(iq,k)*(rat-1.)
     &                           *davu(iq)*dt/3600.
	     endif ! (speed.gt.1.)
           enddo  ! iq loop
          enddo   ! k loop
        endif     ! (nud_uv.eq.2)
        if(nmaxpr.eq.1.and.ktau.lt.5.and.mydiag)then
          print *,'davies out u,v,qg(kl) ',
     &           u(idjd,nlv),v(idjd,nlv),qg(idjd,kl)
        endif

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
