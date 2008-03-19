      subroutine mslp(pmsl,psl,zs,t)
      use cc_mpi, only : mydiag
!     this one will ignore negative zs (i.e. over the ocean)
      parameter (meth=1) ! 0 for original, 1 for other jlm - always now
      include 'newmpar.h'
      include 'const_phys.h'
      include 'parm.h'
      real pmsl(ifull),psl(ifull),zs(ifull),t(ifull,kl)
      include 'sigs.h'
      save lev
      c=grav/stdlapse
      conr=c/rdry
      lev=0
14    lev=lev+1
c     find level just below sig=.9
      if (sig(lev+1).gt..9)go to 14
      con=sig(lev)**(rdry/c)/c
c     if(meth.eq.0)then
c       do iq=1,ifull
c        pmsl(iq)=ps(iq)*(1.+con*zs(iq)/t(iq,lev))**conr
c       enddo
c     endif  ! (meth.eq.0)
      if(meth.eq.1)then
        do iq=1,ifull
         phi1=t(iq,lev)*rdry*(1.-sig(lev))/sig(lev) ! phi of sig(lev) above sfce
         tsurf=t(iq,lev)+phi1*stdlapse/grav
         tav=tsurf+max(0.,zs(iq))*.5*stdlapse/grav
         dlnps=max(0.,zs(iq))/(rdry*tav)
         pmsl(iq)=1.e5*exp(psl(iq)+dlnps)
        enddo
      endif  ! (meth.eq.1)
      if(nmaxpr==1.and.mydiag)then
        print *,'meth,lev,sig(lev) ',meth,lev,sig(lev)
        print *,'zs,t_lev,psl,pmsl ',
     .           zs(idjd),t(idjd,lev),psl(idjd),pmsl(idjd)
      endif
      return
      entry to_psl(pmsl,psl,zs,t)
      do iq=1,ifull
       phi1=t(iq,lev)*rdry*(1.-sig(lev))/sig(lev) ! phi of sig(lev) above sfce
       tsurf=t(iq,lev)+phi1*stdlapse/grav
       tav=tsurf+max(0.,zs(iq))*.5*stdlapse/grav
       dlnps=max(0.,zs(iq))/(rdry*tav)
       psl(iq)=log(1.e-5*pmsl(iq)) -dlnps
      enddo
      if(nmaxpr==1.and.mydiag)then
        print *,'to_psl lev,sig(lev) ',lev,sig(lev)
        print *,'zs,t_lev,psl,pmsl ',
     .           zs(idjd),t(idjd,lev),psl(idjd),pmsl(idjd)
      endif
      return
      end
