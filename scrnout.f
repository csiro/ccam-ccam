      subroutine scrnout(zo,ustar,fg,eg,factch,qsurf,wetfac,qsttg,  ! arrays
     .     qgscrn,tscrn,uscrn,u10,scrrel,                           ! arrays
     .     bprm,cms,chs,fmroot,nalpha)
      use cc_mpi, only : mydiag
      use diag_m
!     allows for zo>zscr from 30/7/04
!     has fixer at bottom to ensure tscrn and qgscrn bounded (& non-neg)
      parameter (ntest=0)   ! ntest= 0 for diags off; ntest= 1 for diags on
c**************** needed sopt for SX5 compiler if ntest=1 (19/3/04)      
      parameter (jlmspec=0) ! jlmspec=0 for original method over land
!     parameter (chn10=.00136733)   ! sea only, same as in sflux via parm.h
      parameter (vkar=.4,zscr=1.8)
c     for nsib=2 don't calculate tscrn here
      include 'newmpar.h'
      include 'arrays.h'
      include 'const_phys.h'
      include 'liqwpar.h'  ! ifullw,qfg,qlg  just for diags
      include 'map.h'
      include 'nsibd.h'    ! rsmin,ivegt,sigmf,tgf,ssdn,res,rmc,tsigmf
      include 'parm.h'
      include 'pbl.h'
      include 'prec.h'     ! just for diags
      include 'scamdim.h'
      include 'sigs.h'
      include 'soil.h'
      include 'soilsnow.h'
      common/work2/dum2a(11*ifull),fhx(ifull),ri(ifull),theta(ifull),
     .        zscronzo(ifull),qstarx(ifull),vmod(ifull),tstarx(ifull)
      common/work3/egg(ifull),evapxf(ifull),Ewww(ifull),fgf(ifull),
     . fgg(ifull),ggflux(ifull),rdg(ifull),rgg(ifull),residf(ifull),
     . ga(ifull),condxpr(ifull),fev(ifull),fes(ifull),
     . ism(ifull),fwtop(ifull),fh(ifull),fm(ifull),rich(ifull),
     . dum3(5*ijk-18*ifull) 
      common/work3b/wblf(ifull,ms),wbfice(ifull,ms),sdepth(ifull,3),
     .          aft(ifull),aft10(ifull),zlog(ifull),
     .          af(ifull),af10(ifull),deltheta(ifull),
     .          afroot(ifull),afroot10(ifull),
     . dum3b(2*ijk-11*ifull-2*ms*ifull)
      common/permsurf/ipsice,ipsea,ipland,iperm(ifull)
      real fg(ifull),eg(ifull),qgscrn(ifull),tscrn(ifull)
     . ,uscrn(ifull),u10(ifull),scrrel(ifull)
      real wetfac(ifull),factch(ifull),qsttg(ifull),ustar(ifull)
     .     ,zo(ifull),qsurf(ifull)  ! qsurf is just work array, only used here

      include 'establ.h'

c     J. McGregor's vector version, using the formulae of M. Dix

      srcp =sig(1)**(rdry/cp)
      ztv=exp(vkar/sqrt(chn10)) /10.  ! proper inverse of ztsea
c     z3onzt=3.*ztv
c     chn3=(vkar/log(z3onzt))**2
      zscronzt=zscr*ztv
      chnscr=(vkar/log(zscronzt))**2

      if(diag.or.ntest.gt.0)then
        if ( mydiag ) then
          print *,'extra diags from scrnout'
          print *,'wetfac qsttg ',wetfac(idjd),qsttg(idjd)
        end if
        call maxmin(qg,'qg',ktau,1.e3,kl)
        call maxmin(qfg,'qf',ktau,1.e3,kl)
        call maxmin(qlg,'ql',ktau,1.e3,kl)
        call maxmin(t,' t',ktau,1.,kl)
        call maxmin(wb,'wb',ktau,1.,ms)
!        call maxmin(tggsn,'tgg',ktau,1.,ms+3)
        call maxmin(tggsn,'tggsn',ktau,1.,3)
        call maxmin(tgg,'tgg',ktau,1.,ms)
        call maxmin(tss,'ts',ktau,1.,1)
        call maxmin(precip,'pr',ktau,1.,1)
        call maxmin(eg,'eg',ktau,1.,1)
        call maxmin(fg,'fg',ktau,1.,1)
        call maxmin(ps,'ps',ktau,.01,1)
      endif
         
c     calculate "surface" specific humidity
!     N.B. qsttg is coming in as before-calling-sib3 value of qs_tss
      if(nalpha.eq.1) then  ! usual
        do iq=1,ifull
         qsurf(iq) = wetfac(iq)*qsttg(iq)
     .      + (1.-wetfac(iq))*min(qg(iq,1),qsttg(iq))  ! for v. cold surfaces
        enddo   ! iq=1,ifull
      else
        do iq=1,ifull
         qtgnet=qsttg(iq)*wetfac(iq) -qg(iq,1)
         if(qtgnet.gt.0.) then
           qsurf(iq) = qsttg(iq)*wetfac(iq)
         else
           qsurf(iq) = .1*qsttg(iq)*wetfac(iq) + .9*qg(iq,1)
         endif
        enddo   ! iq loop
      endif   ! (nalpha.eq.1) then .. else ..

      if(ntest.ne.0.and.mydiag)print *,'scrnout zo ',zo(idjd)
c     do iq=1,ifull
c      if(zo(iq).ge.zscr.or.zo(iq).le.0.)then
c        print *,'scrnout zscr,zo,land,iq ',zscr,zo(iq),land(iq),iq
c        stop
c      endif
c     enddo
      do iq=1,ifull
       zscrr=max(zscr,zo(iq)+.05)
       zlog(iq)=log(zo(iq))
       zscronzo(iq)=zscrr/zo(iq)
       afroot(iq)=vkar/(log(zscrr)-zlog(iq))
       af(iq)=afroot(iq)*afroot(iq)
c      constants for calculating 3m, 10m winds
c      afroot3(iq) = vkar/(log(3.)-zlog(iq))
c      af3(iq) = afroot3(iq)*afroot3(iq)
       afroot10(iq) = vkar/(log(10.)-zlog(iq))
       af10(iq) = afroot10(iq)*afroot10(iq)
       aft(iq)=vkar*afroot(iq)/(2. + log(zscrr)-zlog(iq))
       aft10(iq)=vkar* afroot10(iq) /(2. + log(10.)-zlog(iq))
      enddo   ! iq loop
      if(ntest.eq.1.and.mydiag)print *,'scrnout zlog ',zlog(idjd)

!cdir nodep
      do ip=ipland+1,ipsea      ! sice and sea points
       iq=iperm(ip)
	if(fracice(iq).lt..5)then  ! treat sice as mainly sea here
         aft(iq)=chnscr
         aft10(iq)=chn10
	endif  ! (fracice(iq).lt..5)
c      following for testing - not useful
c      z1onzt=300.*rdry*(1.-sig(1))*ztv /9.806
c      chnsea=(vkar/log(z1onzt))**2   
c      aft(iq)=chnsea
c      aft10(iq)=chnsea
      enddo

      ymin  =100.
      y10min=100.
!cdir nodep
      do iq=1,ifull     ! all valid points (not just sea, for DARLAM too)
!      N.B. removing the next 4 comments causes vopt on 212 compiler
!      to crash (npanels=5 in parameter statement) -exceeds loopcnt,
!      even though ipsea=ifull (=13824 << loopcnt=24000 in make). 
!      When comments are left in, hopt is fine.
c       if(npanels.eq.0)then
cc        iq=iperm(ip)           ! for DARLAM
c       else
c         iq=ip                  ! for C-C
c       endif
        theta(iq)=t(iq,1)/srcp   ! can use same value from sflux, as diag
        rho=ps(iq)/(rdry*tss(iq))

c       calculate stable / unstable criterion tstarx(iq); also qstarx(iq)
        tstarx(iq) = -fg(iq) / ( rho*cp*ustar(iq))
        qstarx(iq) = -eg(iq) / ( rho*hl*ustar(iq))

        c  = -9.806*tstarx(iq)*zscr*af(iq)*sqrt(af(iq))
     .        /(aft(iq) *tss(iq)*ustar(iq)**2)
c       c3 = -9.806*tstarx(iq)*3.*af3(iq)*sqrt(af3(iq))
c    .        /(aft3(iq)*tss(iq)*ustar(iq)**2)
        c10= -9.806*tstarx(iq)*10.*af10(iq)*sqrt(af10(iq)) ! fixed 1/5/03
     .       /(aft10(iq)*tss(iq)*ustar(iq)**2)

c       section for solving for screen height values of various parameters

        if ( tstarx(iq) .ge. 0. ) then
c         stable case     (tstarx(iq)>0,   c<0, giving rich>0)
c         Solve quadratic for rich, form ax^2 + bx + c
c         a = bprm ;   b = 1.
c         Because a*c < 0 this has two real solutions, only one of which is
c         positive and so appropriate for the stable case.
c         There is a slight error here because we use tg rather than theta1
c         to calculate c. However this is a second order effect.
!         rich     = ( -b  + sqrt(b*b-4.*a*c  ) ) /(2.*a)
          rich(iq) = ( -1. + sqrt(1.-4.*bprm*c) ) /(2.*bprm)
c         rich3    = ( -1. + sqrt(1.-4.*bprm*c3) ) /(2.*bprm)
          rich10   = ( -1. + sqrt(1.-4.*bprm*c10) ) /(2.*bprm)

c         Use the positive solution
          fm(iq) = max(fmroot*fmroot,1./(1.+bprm*rich(iq)  )**2)
c         fm3  = max(fmroot*fmroot,1./(1.+bprm*rich3 )**2)
          fm10 = max(fmroot*fmroot,1./(1.+bprm*rich10)**2)
          con = afroot(iq)*min(fmroot,(1.+bprm*rich(iq)))/aft(iq)

        else
c         unstable case
c         Starting value of the Richardson number.
c         The iteration calculates fm and fh as above but using the Richardson
c         number directly rather than x
!         y0=1.  ! testing - used till April 04
c         y0 = sqrt(abs(ri(iq)))  ! jlm test April 04
          y0 = sqrt(c)
          root=sqrt(zscronzo(iq))
          cm=2.*bprm*cms*af(iq)*root
          ch=2.*bprm*chs*aft(iq)*factch(iq)*root
          ffm=1.+2.*bprm*y0*y0/(1.+cm*y0)
          ffh=1.+2.*bprm*y0*y0/(1.+ch*y0)
          y1 = sqrt(c *sqrt(ffm*ffm*ffm)/ffh)
          ffm=1.+2.*bprm*y1*y1/(1.+cm*y1)
          ffh=1.+2.*bprm*y1*y1/(1.+ch*y1)
          if(ntest.ne.0.and.iq.eq.idjd.and.mydiag) then
            fm0=1.+2.*bprm*y0*y0/(1.+cm*y0)
            fh0=1.+2.*bprm*y0*y0/(1.+ch*y0)
            print *,'y0,fm0,fh0 ',y0,fm0,fh0
            print *,'y1,fm1,fh1 ',y,ffm,ffh
          endif
          y = sqrt(c *sqrt(ffm*ffm*ffm)/ffh)
          fm(iq)=1.+2.*bprm*y*y/(1.+cm*y)
          fh(iq)=1.+2.*bprm*y*y/(1.+ch*y)
          ymin=min(y,ymin)
          con=sqrt(af(iq)*fm(iq))/(aft(iq)*fh(iq))
	   rich(iq)=-y*y

c         iterations for 10 m winds
          y10 = sqrt(c10)
          root=sqrt(10./zo(iq))
          cm10=2.*bprm*cms*af10(iq)*root
          ch10=2.*bprm*chs*aft10(iq)*factch(iq)*root
          fm10=1.+2.*bprm*y10*y10/(1.+cm10*y10)
          fh10=1.+2.*bprm*y10*y10/(1.+ch10*y10)
          y10 = sqrt(c10 *sqrt(fm10*fm10*fm10)/fh10)  
          fm10=1.+2.*bprm*y10*y10/(1.+cm10*y10)
          y10min=min(y10,y10min)
        endif

c       calculate correction to surface temperature deltheta and moisture
        deltheta(iq)=tstarx(iq)*con
        deltq = qstarx(iq)*con

c       screen wind speeds
        uscrn(iq) = ustar(iq)/(afroot(iq)*sqrt(fm(iq)))
c       u3(iq) = ustar(iq)/(afroot3(iq)*sqrt(fm3))
        u10(iq) = ustar(iq)/(afroot10(iq)*sqrt(fm10))

c       screen specific humidity
c       N.B. over unstable sea points, sometimes will get supersat.qgscrn
        qgscrn(iq) = qsurf(iq) + deltq   ! will only use for sea

        if(ntest.ne.0.and.iq.eq.idjd.and.mydiag)then
          print *,' scrnout test for id,jd,idjd,land,sice ',
     .                               id,jd,idjd,land(iq),sice(iq)
          print *,'tss,theta,ps,con ',tss(iq),theta(iq),ps(iq),con
          print *,'nalpha,fh,fm,fm10 ',nalpha,fh(iq),fm(iq),fm10
	   print *,'deltq,qsurf,qgscrn,qg1 ',
     .             deltq,qsurf(iq),qgscrn(iq),qg(iq,1)
          print *,'c,rich,zscr ',c,rich(iq),zscr
        endif ! (ntest.ne.0.and.iq.eq.idjd)
      enddo   ! iq=1,ifull

      if(mydiag.and.(diag.or.ntest.eq.2))then  ! special test of 38 m, stable sea
       iq=idjd
       zscr38=zmin
       zscronzt38=zscr38*ztv
       chnscr38=(vkar/log(zscronzt38))**2
       afroot38=vkar/(log(zscr38)-zlog(iq))
       af38=afroot38*afroot38
       aft38=chnscr38
       rho=ps(iq)/(rdry*tss(iq))
       zscronzo38=zscr38/zo(iq)
!      tstarx(iq) = -fg(iq) / ( rho*cp*ustar(iq))
!      qstarx(iq) = -eg(iq) / ( rho*hl*ustar(iq))
       c38 = -9.806*tstarx(iq)*zscr38*af38*sqrt(af38)
     .                   /(aft38 *tss(iq)*ustar(iq)**2)
!      if ( tstarx(iq) .ge. 0. ) then
!        rich(iq)   = ( -1. + sqrt(1.-4.*bprm*c38) ) /(2.*bprm)
         ffm   = max(fmroot*fmroot,1./(1.+bprm*rich(iq)  )**2)
c        con = afroot38*min(fmroot,(1.+bprm*rich(iq)))/aft38
         con = afroot38*min(1./fmroot,1.+bprm*rich(iq))/aft38
!       deltheta(iq)=tstarx(iq)*con
        deltq = qstarx(iq)*con
!       uscrn(iq) = ustar(iq)/(afroot38*sqrt(ffm))
	 print *,'38 m diagnostics'
	 print *,'zscr,af38,afroot38 ',zscr38,af38,afroot38
	 print *,'fmroot,fm,con ',fmroot,ffm,con
	 print *,'chnscr38,zo,zt,ztv ',chnscr38,zo(iq),1./ztv,ztv
	 print *,'aft38,tstarx(iq),rich38 ',aft38,tstarx(iq),rich(iq)
	 print *,'tgg1,tgg2,thscrn.theta1 ',
     .           tgg(iq,1),tgg(iq,2),tgg(iq,2)+deltheta(iq),theta(iq)
        print *,'uscrn,uv1 ',uscrn(iq),sqrt(u(iq,1)**2+v(iq,1)**2)
      endif  ! (diag.or.ntest.eq.2)

!cdir nodep
      do ip=1,ipsice      ! just land and sice points
       iq=iperm(ip)
       tscrn(iq) = tss(iq) + deltheta(iq) 
       es = establ(tscrn(iq))
       constz=ps(iq)-es
       qs= 0.622*es/constz
!      I don't believe treatment of qgscrn over land & ice, so use:
       qgscrn(iq)=min(qs,(zmin*qsurf(iq)+zscr*qg(iq,1))/(zmin+zscr))
       scrrel(iq) = 100.*qgscrn(iq)/qs
      enddo

!cdir nodep
      do ip=ipsice+1,ipsea      ! just sea points
       iq=iperm(ip)
       tscrn(iq) = tgg(iq,2) + deltheta(iq) 
       es = establ(tscrn(iq))
       constz=ps(iq)-es
       qsttg(iq)= 0.622*es/constz
       scrrel(iq) = max(0.,100.*(qgscrn(iq)/qsttg(iq)))
       if(ntest.eq.3)then  ! print all stable sea points
	  if(fg(iq).lt.0.)then
	    print *,'iq,tgg2,tscrn,theta,ri,rich s ',
     .              iq,tgg(iq,2),tscrn(iq),theta(iq),ri(iq),rich(iq)
           print *,'fhx,fm ',fhx(iq)/vmod(iq),fm(iq)
	  endif
       endif               ! (ntest.eq.3)
       if(ntest.eq.4)then  ! print all unstable sea points
	  if(fg(iq).gt.0.)then
           zscr38=zmin
           zscronzt38=zscr38*ztv
           chnscr38=(vkar/log(zscronzt38))**2
           afroot38=vkar/(log(zscr38)-zlog(iq))
           af38=afroot38*afroot38
           aft38=chnscr38
           zscronzo38=zscr38/zo(iq)
           c38 = -9.806*tstarx(iq)*zscr38*af38*sqrt(af38)
     .                   /(aft38 *tss(iq)*ustar(iq)**2)
!          calculate an approximate form for ri38
           ri38=-c38*vmod(iq)*fm(iq)*sqrt(fm(iq))/fhx(iq)
	    print *,'iq,tgg2,tscrn,theta,ri,rich u ',
     .              iq,tgg(iq,2),tscrn(iq),theta(iq),ri(iq),rich(iq)
           print *,'fhx,fh,fm,ri38 ',fhx(iq)/vmod(iq),fh(iq),fm(iq),ri38
	  endif
       endif              ! (ntest.eq.4)
      enddo

      if(mydiag.and.(diag.or.ntest.ne.0))then
	 iq=idjd
        print *,' scrnout test for id,jd,idjd,land,sice ',
     .                             id,jd,idjd,land(iq),sice(iq)
        print *,'zo,rho,cp ',zo(iq),rho,cp
        print *,'hl,ustar,factch ',hl,ustar(iq),factch(iq)
        print *,'bprm,cms,chs ',bprm,cms,chs
        print *,'fg,eg,fgg,fgf ',fg(iq),eg(iq),fgg(iq),fgf(iq)
        print *,'tss,theta,ps ',tss(iq),theta(iq),ps(iq)
        print *,'nalpha,fh,fm ',nalpha,fh(iq),fm(iq)
        print *,'qstarx,tstarx,tscrna ',
     .           qstarx(iq),tstarx(iq),tss(iq)+deltheta(iq)
	 print *,'wetfac,qsttg ',wetfac(iq),qsttg(iq)
	 print *,'qsurf,qgscrn,qg1,scrrel ',
     .           qsurf(iq),qgscrn(iq),qg(iq,1),scrrel(iq)
        print *,'af,afroot,aft ',af(iq),afroot(iq),aft(iq)
        print *,'afroot10,tstarx ',afroot10(iq),tstarx(iq)
        print *,'rich,zscr ',rich(iq),zscr
        es = establ(tscrn(iq))
        constz=ps(iq)-es
        qs= 0.622*es/constz
        print *,'tss,tscrn,theta,uscrn,u10 ',
     .           tss(iq),tscrn(iq),theta(iq),uscrn(iq),u10(iq)
        print *,'deltheta,qs ',deltheta(iq),qs
      endif ! (diag.or.ntest.ne.0)

      if(jlmspec.eq.1)then ! special treatment of bare soil & veg
         do iq=1,ifull    ! following not checked out lately
          if(land(iq))then
!           bare ground part here
            tstarx(iq) = -fgg(iq) / ( rho*cp*ustar(iq))
            c=-9.806*tstarx(iq)*zscr*af(iq) *sqrt(af(iq) )
     .       /(aft(iq) *tgg(iq,1)*ustar(iq)**2)
            if(tstarx(iq) .ge. 0.)then
              rich(iq) = ( -1. + sqrt(1.-4.*bprm*c) ) /(2.*bprm)
              fm(iq) = max(fmroot*fmroot,1./(1.+bprm*rich(iq) )**2)
              con = tstarx(iq)*afroot(iq)*
     .                      min(fmroot,(1.+bprm*rich(iq)))/aft(iq)
            else            ! unstable
!             y = sqrt(abs(ri(iq)))  ! jlm test April 04
              y = sqrt(c)
              root=sqrt(zscronzo(iq))
              cm=2.*bprm*cms*af(iq)*root
              ch=2.*bprm*chs*aft(iq)*factch(iq)*root
              ffm=1.+2.*bprm*y*y/(1.+cm*y)
              ffh=1.+2.*bprm*y*y/(1.+ch*y)
              y = sqrt(c *sqrt(ffm*ffm*ffm)/ffh)
              ffm=1.+2.*bprm*y*y/(1.+cm*y)
              ffh=1.+2.*bprm*y*y/(1.+ch*y)
              y = sqrt(c *sqrt(ffm*ffm*ffm)/ffh)
              fm(iq)=1.+2.*bprm*y*y/(1.+cm*y)
              fh(iq)=1.+2.*bprm*y*y/(1.+ch*y)
              ymin=min(y,ymin)
              con=sqrt(af(iq)*fm(iq))/(aft(iq)*fh(iq))
		rich(iq)=-y*y
            endif  ! ( tstarx(iq) .ge. 0. )
            delthetg=tstarx(iq)*con
            if(mydiag.and.ntest.ne.0.and.iq.eq.idjd)then
              print *,'iq,tgg1,delthetg ',iq,tgg(iq,1),delthetg
              print *,'tstarx(iq),rich,y ',tstarx(iq),rich(iq),y
            endif

!           vegetation part here
            tstarx(iq) = -fgf(iq)/(rho*cp*ustar(iq))
            c=-9.806*tstarx(iq)*zscr*af(iq) *sqrt(af(iq) )
     .         /(aft(iq) *tgf(iq)*ustar(iq)**2)
            if(tstarx(iq) .ge. 0.)then
              rich   = ( -1. + sqrt(1.-4.*bprm*c) ) /(2.*bprm)
              fm(iq) = max(fmroot*fmroot,1./(1.+bprm*rich(iq)  )**2)
              con = afroot(iq)*min(fmroot,(1.+bprm*rich(iq)))/aft(iq)
            else          ! unstable
!             y = sqrt(abs(ri(iq)))  ! jlm test April 04
              y = sqrt(c)
              root=sqrt(zscronzo(iq))
              cm=2.*bprm*cms*af(iq)*root
              ch=2.*bprm*chs*aft(iq)*factch(iq)*root
              ffm=1.+2.*bprm*y*y/(1.+cm*y)
              ffh=1.+2.*bprm*y*y/(1.+ch*y)
              y = sqrt(c *sqrt(ffm*ffm*ffm)/ffh)
              ffm=1.+2.*bprm*y*y/(1.+cm*y)
              ffh=1.+2.*bprm*y*y/(1.+ch*y)
              y = sqrt(c *sqrt(ffm*ffm*ffm)/ffh)
              fm(iq)=1.+2.*bprm*y*y/(1.+cm*y)
              fh(iq)=1.+2.*bprm*y*y/(1.+ch*y)
              ymin=min(y,ymin)
              con=sqrt(af(iq)*fm(iq))/(aft(iq)*fh(iq))
		rich(iq)=-y*y
            endif  ! ( tstarx(iq) .ge. 0. )
            delthetf=tstarx(iq)*con

            tscrn(iq) = tss(iq)
     .           +delthetf*tsigmf(iq)+delthetg*(1.-tsigmf(iq))
            if(mydiag.and.ntest.ne.0.and.iq.eq.idjd)then
              print *,'tgf,delthetf,tscrn ',
     .                 tgf(iq),delthetf,tscrn(iq)
              print *,'tstarx(iq),rich,y ',tstarx(iq),rich(iq),y
            endif
          endif  ! (land(iq))
         enddo
      endif  ! (jlmspec.eq.1)

      if(ntest.eq.1)then
        numt=0
        numq=0
	 do iq=1,ifull
	  if(land(iq).or.sice(iq))then
           if(tscrn(iq).lt.min(theta(iq),tss(iq)).or.
     .        tscrn(iq).gt.max(theta(iq),tss(iq))    )then
             numt=numt+1
             print *,'strange T result iq,land,sice,tss,tscrn ',
     .                     iq,land(iq),sice(iq),tss(iq),tscrn(iq)
             print *,'tgg1,t1,theta,ri,numt ',
     .                tgg(iq,1),t(iq,1),theta(iq),rich(iq),numt
             if(sice(iq))print *,'tgg2,tgg3,fracice ',
     .                            tgg(iq,2),tgg(iq,3),fracice(iq)  	      
           endif
         else
           if(tscrn(iq).lt.min(theta(iq),tgg(iq,2)).or.
     .        tscrn(iq).gt.max(theta(iq),tgg(iq,2))    )then
             numt=numt+1
             print *,'strange T result iq,land,sice,tss,tscrn ',
     .                     iq,land(iq),sice(iq),tss(iq),tscrn(iq)
             print *,'tgg2,t1,theta,ri,numt ',
     .                tgg(iq,2),t(iq,1),theta(iq),rich(iq),numt
           endif
         endif   ! (land(iq).or.sice(iq)) .. else ..
	  if(qgscrn(iq).lt.min(qg(iq,1),qsurf(iq)).or.
     .      qgscrn(iq).gt.max(qg(iq,1),qsurf(iq))    )then
           numq=numq+1
           print *,'strange qg result iq,land,sice,qsurf,qgscrn ',
     .      iq,land(iq),sice(iq),qsurf(iq),qgscrn(iq)
           print *,'qg1,t1,ri ',qg(iq,1),t(iq,1),ri(iq),numq
         endif
	 enddo
      endif   ! (ntest.eq.1)

!     consistency fixer for tscrn and qgscrn (odd pts)  March 2004
      do iq=1,ifull
       if(land(iq).or.sice(iq))then
         tscrn(iq)=max(tscrn(iq),min(tss(iq),theta(iq)))
         tscrn(iq)=min(tscrn(iq),max(tss(iq),theta(iq)))
	else
         tscrn(iq)=max(tscrn(iq),min(tgg(iq,2),theta(iq)))
         tscrn(iq)=min(tscrn(iq),max(tgg(iq,2),theta(iq)))
	endif  ! (land(iq).or.sice(iq))
       qgscrn(iq)=max(qgscrn(iq),min(qsurf(iq),qg(iq,1)))
       qgscrn(iq)=min(qgscrn(iq),max(qsurf(iq),qg(iq,1)))
      enddo
   
      tscrn(:)=tscrn(:)-.018  ! apply 1.8 m approx. adiab. corrn.

      if(ymin.lt.0.)stop ' poor convergence for y in scrnout'
      if(y10min.lt.0.)stop ' poor convergence 10 m wind in scrnout'
      if(ntest.eq.2)then
        print *,' tscreen: ',(tscrn(iq),iq=il/2,ifull,il)
        print *,' qgscrn: ',(qgscrn(iq),iq=il/2,ifull,il)
      endif
      return
      end
