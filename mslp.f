      subroutine mslp(pmsl,psl,zs,t)
      parameter (meth=1) ! 0 for original, 1 for other jlm - always now
      include 'newmpar.h'
      include 'constant.h'
      include 'parm.h'
      real pmsl(ifull),psl(ifull),zs(ifull),t(ifull,kl)
      include 'sigs.h'
      save lev
      c=g/6.5e-3
      conr=c/r
      lev=0
14    lev=lev+1
c     find level just below sig=.9
      if (sig(lev+1).gt..9)go to 14
      con=sig(lev)**(r/c)/c
c     if(meth.eq.0)then
c       do iq=1,ifull
c        pmsl(iq)=ps(iq)*(1.+con*zs(iq)/t(iq,lev))**conr
c       enddo
c     endif  ! (meth.eq.0)
      if(meth.eq.1)then
        do iq=1,ifull
         phi1=t(iq,lev)*r*(1.-sig(lev))/sig(lev) ! phi of sig(lev) above sfce
         tsurf=t(iq,lev)+phi1*.0065/g
         tav=tsurf+zs(iq)*.5*.0065/g
         dlnps=zs(iq)/(r*tav)
         pmsl(iq)=1.e5*exp(psl(iq)+dlnps)
        enddo
      endif  ! (meth.eq.1)
      if(nmaxpr.eq.1)then
        print *,'meth,lev,sig(lev) ',meth,lev,sig(lev)
        print *,'zs,t_lev,psl,pmsl ',
     .           zs(idjd),t(idjd,lev),psl(idjd),pmsl(idjd)
      endif
      return
      entry to_psl(pmsl,psl,zs,t)
      do iq=1,ifull
       phi1=t(iq,lev)*r*(1.-sig(lev))/sig(lev) ! phi of sig(lev) above sfce
       tsurf=t(iq,lev)+phi1*.0065/g
       tav=tsurf+zs(iq)*.5*.0065/g
       dlnps=zs(iq)/(r*tav)
       psl(iq)=log(1.e-5*pmsl(iq)) -dlnps
      enddo
      if(nmaxpr.eq.1)then
        print *,'to_psl lev,sig(lev) ',lev,sig(lev)
        print *,'zs,t_lev,psl,pmsl ',
     .           zs(idjd),t(idjd,lev),psl(idjd),pmsl(idjd)
      endif
      return
      end
