      subroutine scrnout(zo,ustar,fg,eg,factch,qsurf,wetfac,qsttg,  ! arrays
     .     qgscrn,tscrn,uscrn,u10,scrrel,                           ! arrays
     .     bprm,cms,chs,fmroot,nalpha)
      parameter (ntest=0)   ! ntest= 0 for diags off; ntest= 1 for diags on
      parameter (jlmspec=0) ! jlmspec=0 for original method over land
      parameter (chn10=.00136733)   ! sea only, same as in sflux
      parameter (vkar=.4,zscr=1.8)
c     for nsib=2 don't calculate tscrn here
      include 'newmpar.h'
      include 'arrays.h'
      include 'const_phys.h'
      include 'map.h'
      include 'nsibd.h'    ! rsmin,ivegt,sigmf,tgf,ssdn,res,rmc,tsigmf
      include 'parm.h'
      include 'pbl.h'
      include 'permsurf.h'
      include 'scamdim.h'
      include 'sigs.h'
      include 'soil.h'
      include 'soilsnow.h'
      common/work3/egg(ifull),evapxf(ifull),Ewww(ifull),fgf(ifull),
     . fgg(ifull),ggflux(ifull),rdg(ifull),rgg(ifull),residf(ifull),
     . ga(ifull),condxpr(ifull),fev(ifull),fes(ifull),
     . ism(ifull),fwtop(ifull),
     . dum3(2*ijk-15*ifull),uav(ifull,kl),vav(ifull,kl),spare(ifull,kl)   
      common/work3b/wblf(ifull,ms),wbfice(ifull,ms),sdepth(ifull,3),
     .          aft(ifull),aft10(ifull),zlog(ifull),
     .          af(ifull),af10(ifull),deltheta(ifull),
     .          afroot(ifull),afroot10(ifull),
     . dum3b(2*ijk-11*ifull-2*ms*ifull)
      real fg(ifull),eg(ifull),qgscrn(ifull),tscrn(ifull)
     . ,uscrn(ifull),u10(ifull),scrrel(ifull)
      real wetfac(ifull),factch(ifull),qsttg(ifull),ustar(ifull)
     .     ,zo(ifull),qsurf(ifull)  ! qsurf is just work array, only used here

c     J. McGregor's vector version
c     with prior contributions from K. Walsh and M. Dix
c     tolerance for iteration,iteration maximum

      data scrtol/.01/,itermax/8/
      include 'establ.h'

      srcp =sig(1)**(rdry/cp)
      ztv=exp(vkar/sqrt(chn10)) /10.  ! proper inverse of ztsea
c     z3onzt=3.*ztv
c     chn3=(vkar/log(z3onzt))**2
      zscronzt=zscr*ztv
      chnscr=(vkar/log(zscronzt))**2

      if(ntest.eq.1)print *,'scrnout wetfac qsttg ',
     .                               wetfac(idjd),qsttg(idjd)
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
        enddo   ! iq=1,ifull
      endif   ! (nalpha.eq.1) then .. else ..

      if(ntest.eq.1)print *,'scrnout zo ',zo(idjd)
      do iq=1,ifull
       zlog(iq)=log(zo(iq))
       zscronzo=zscr/zo(iq)
       afroot(iq)=vkar/(log(zscr)-zlog(iq))
       af(iq)=afroot(iq)*afroot(iq)
c      constants for calculating 3m, 10m winds
c      afroot3(iq) = vkar/(log(3.)-zlog(iq))
c      af3(iq) = afroot3(iq)*afroot3(iq)
       afroot10(iq) = vkar/(log(10.)-zlog(iq))
       af10(iq) = afroot10(iq)*afroot10(iq)
      enddo   ! iq=1,ifull

      if(ntest.eq.1)print *,'scrnout zlog ',zlog(idjd)
!cdir nodep
      do ip=1,ipsice      ! just land and sice points
       iq=iperm(ip)
       aft(iq)=vkar*afroot(iq)/(2. + log(zscr)-zlog(iq))
c      aft3(iq)=vkar* afroot3(iq) /(2. + log(3.)-zlog(iq))
       aft10(iq)=vkar* afroot10(iq) /(2. + log(10.)-zlog(iq))
      enddo

!cdir nodep
      do ip=ipsice+1,ipsea      ! just sea points
       iq=iperm(ip)
       aft(iq)=chnscr
c      aft3(iq)=chn3
       aft10(iq)=chn10
      enddo

      ymin  =100.
      y10min=100.
!cdir nodep
      do ip=1,ipsea      ! all valid points (not just sea, for DARLAM too)
!      N.B. removing the next 4 comments causes vopt on 212 compiler
!      to crash (npanels=5 in parameter statement) -exceeds loopcnt,
!      even though ipsea=ifull (=13824 << loopcnt=24000 in make). 
!      When comments are left in, hopt is fine.
c       if(npanels.eq.0)then
          iq=iperm(ip)           ! for DARLAM
c       else
c         iq=ip                  ! for C-C
c       endif
        theta=t(iq,1)/srcp
        rho=ps(iq)/(rdry*tss(iq))
        zscronzo=zscr/zo(iq)

c       calculate stable / unstable criterion tstarx; also qstarx
        tstarx = -fg(iq) / ( rho*cp*ustar(iq))
        qstarx = -eg(iq) / ( rho*hl*ustar(iq))

        c  = -grav*tstarx*zscr*af(iq)*sqrt(af(iq))
     .        /(aft(iq) *tss(iq)*ustar(iq)**2)
c       c3 = -grav*tstarx*3.*af3(iq)*sqrt(af3(iq))
c    .        /(aft3(iq)*tss(iq)*ustar(iq)**2)
        c10= -grav*tstarx*10.*af10(iq)*sqrt(af10(iq)) ! fixed 1/5/03
     .       /(aft10(iq)*tss(iq)*ustar(iq)**2)

c       section for solving for screen height values of various parameters

        if ( tstarx .ge. 0. ) then
c         stable case     (tstarx>0,   c<0, giving rich>0)
c         Solve quadratic for rich, form ax^2 + bx + c
c         a = bprm ;   b = 1.
c         Because a*c < 0 this has two real solutions, only one of which is
c         positive and so appropriate for the stable case.
c         There is a slight error here because we use tg rather than theta1
c         to calculate c. However this is a second order effect.
!         rich   = ( -b  + sqrt(b*b-4.*a*c  ) ) /(2.*a)
          rich   = ( -1. + sqrt(1.-4.*bprm*c) ) /(2.*bprm)
c         rich3  = ( -1. + sqrt(1.-4.*bprm*c3) ) /(2.*bprm)
          rich10 = ( -1. + sqrt(1.-4.*bprm*c10) ) /(2.*bprm)

c         Use the positive solution
          fm   = max(fmroot*fmroot,1./(1.+bprm*rich  )**2)
c         fm3  = max(fmroot*fmroot,1./(1.+bprm*rich3 )**2)
          fm10 = max(fmroot*fmroot,1./(1.+bprm*rich10)**2)
          con = afroot(iq)*min(fmroot,(1.+bprm*rich))/aft(iq)

        else
c         unstable case
c         Starting value of the Richardson number.
c         The iteration calculates fm and fh as above but using the Richardson
c         number directly rather than x
          y0 = sqrt(c)
          y0=1.  ! testing
          root=sqrt(zscronzo)
          cm=2.*bprm*cms*af(iq)*root
          ch=2.*bprm*chs*aft(iq)*factch(iq)*root
          fm=1.+2.*bprm*y0*y0/(1.+cm*y0)
          fh=1.+2.*bprm*y0*y0/(1.+ch*y0)
          y1 = sqrt(c *sqrt(fm*fm*fm)/fh)
          fm=1.+2.*bprm*y1*y1/(1.+cm*y1)
          fh=1.+2.*bprm*y1*y1/(1.+ch*y1)
          if(ntest.ne.0.and.iq.eq.idjd) then
            fm0=1.+2.*bprm*y0*y0/(1.+cm*y0)
            fh0=1.+2.*bprm*y0*y0/(1.+ch*y0)
            print *,'y0,fm0,fh0 ',y0,fm0,fh0
            print *,'y1,fm1,fh1 ',y,fm,fh
          endif
          y = sqrt(c *sqrt(fm*fm*fm)/fh)
          fm=1.+2.*bprm*y*y/(1.+cm*y)
          fh=1.+2.*bprm*y*y/(1.+ch*y)
          ymin=min(y,ymin)
          con=sqrt(af(iq)*fm)/(aft(iq)*fh)

cc         iterations for 3 m winds
c          y3 = sqrt(c3)
c          root=sqrt(3./zo(iq))
c          cm3=2.*bprm*cms*af3(iq)*root
c          ch3=2.*bprm*chs*aft3(iq)*factch(iq)*root
c          fm3=1.+2.*bprm*y3*y3/(1.+cm3*y3)
c          fh3=1.+2.*bprm*y3*y3/(1.+ch3*y3)
cc         y3 = sqrt(c3) * fm3**0.75 / sqrt(fh3)
c          y3 = sqrt(c3 *sqrt(fm3*fm3*fm3)/fh3)  ! jlm Fri  05-05-1995
c          fm3=1.+2.*bprm*y3*y3/(1.+cm3*y3)

c         iterations for 10 m winds
          y10 = sqrt(c10)
          root=sqrt(10./zo(iq))
          cm10=2.*bprm*cms*af10(iq)*root
          ch10=2.*bprm*chs*aft10(iq)*factch(iq)*root
          fm10=1.+2.*bprm*y10*y10/(1.+cm10*y10)
          fh10=1.+2.*bprm*y10*y10/(1.+ch10*y10)
c         y10 = sqrt(c10) * fm10**0.75 / sqrt(fh10)
          y10 = sqrt(c10 *sqrt(fm10*fm10*fm10)/fh10)  ! jlm Fri  05-05-1995
          fm10=1.+2.*bprm*y10*y10/(1.+cm10*y10)
          y10min=min(y10,y10min)
        endif

c       calculate correction to surface temperature deltheta and moisture
        deltheta(iq)=tstarx*con
        deltq = qstarx*con

c       screen wind speeds
        uscrn(iq) = ustar(iq)/(afroot(iq)*sqrt(fm))
c       u3(iq) = ustar(iq)/(afroot3(iq)*sqrt(fm3))
        u10(iq) = ustar(iq)/(afroot10(iq)*sqrt(fm10))

c       screen specific humidity
        qgscrn(iq) = qsurf(iq) + deltq

        if(ntest.ne.0.and.iq.eq.idjd)then
          print *,' scrnout test for id,jd,land,sice ',
     .                               id,jd,land(iq),sice(iq)
          print *,'zo,rho,cp ',zo(iq),rho,cp
          print *,'hl,ustar,factch ',hl,ustar(iq),factch(iq)
          print *,'bprm,cms,chs ',bprm,cms,chs
          print *,'fg,eg,fgg,fgf ',fg(iq),eg(iq),fgg(iq),fgf(iq)
          print *,'tss,theta,ps ',tss(iq),theta,ps(iq)
          print *,'nalpha,fh,fm,fm10 ',nalpha,fh,fm,fm10
          print *,'qstarx,tstarx ',qstarx,tstarx
          print *,'af,afroot,aft ',af(iq),afroot(iq),aft(iq)
          print *,'afroot10,tstarx ',afroot10(iq),tstarx
          print *,'c,rich,zscr ',c,rich,zscr
        endif ! (ntest.ne.0.and.iq.eq.idjd)
      enddo   ! iq=1,ifull

!cdir nodep
      do ip=1,ipsice      ! just land and sice points
       iq=iperm(ip)
       tscrn(iq) = tss(iq) + deltheta(iq)  ! was only for nsib.ne.2
       es = establ(tscrn(iq))
       constz=ps(iq)-es
       qs= 0.622*es/constz
       scrrel(iq) = max(0.,100.*(qgscrn(iq)/qs))
       if(ntest.eq.1.and.tscrn(iq).lt.t(iq,1).and.
     .                 tss(iq).gt.t(iq,1))then
        print *,'strange result iq,land,sice,tss,tscrn ',
     .    iq,land(iq),sice(iq),tss(iq),tscrn(iq)
        print *,'tgg1,t1,theta ',tgg(iq,1),t(iq,1),theta
c       fm0=1.+2.*bprm*y0*y0/(1.+cm*y0)
c       fh0=1.+2.*bprm*y0*y0/(1.+ch*y0)
c       fm1=1.+2.*bprm*y1*y1/(1.+cm*y1)
c       fh1=1.+2.*bprm*y1*y1/(1.+ch*y1)
c       print *,'y0,fm0,fh0 ',y0,fm0,fh0
c       print *,'y1,fm1,fh1 ',y1,fm1,fh1
c       print *,'y ,fm ,fh  ',y,fm,fh
c       if(tstarx.lt.0.)rich=-y*y
c       print *,'tstarx,rich,fg ',tstarx,rich,fg(iq)
       endif
      enddo
      if(ntest.ne.0.or.diag)then
        iq=idjd
        es = establ(tscrn(iq))
        constz=ps(iq)-es
        qs= 0.622*es/constz
        print *,'wetfac,qsttg,qg1 ',
     .           wetfac(iq),qsttg(iq),qg(iq,1)
          print *,'tss,tscrn,uscrn,u10 ',
     .             tss(iq),tscrn(iq),uscrn(iq),u10(iq)
          print *,'deltheta,qs,qsurf(iq) ',
     .             deltheta(iq),qs,qsurf(iq)
          print *,'qgscrn,scrrel ',qgscrn(iq),scrrel(iq)
       endif

!cdir nodep
      do ip=ipsice+1,ipsea      ! just sea points
       iq=iperm(ip)
       tscrn(iq) = tgg(iq,2) + deltheta(iq)  ! tgg(iq,2) is skin temperature
       es = establ(tscrn(iq))
       constz=ps(iq)-es
       qsttg(iq)= 0.622*es/constz
       scrrel(iq) = max(0.,100.*(qgscrn(iq)/qsttg(iq)))
       if(ntest.ne.0.and.iq.eq.idjd)then
          print *,'wetfac,qsttg,qg(iq,1) ',
     .             wetfac(iq),qsttg(iq),qg(iq,1)
          print *,'tscrn,uscrn,u10 ',
     .             tscrn(iq),uscrn(iq),u10(iq)
          print *,'deltheta,uscrn,qsurf(iq) ',
     .             deltheta(iq),uscrn(iq),qsurf(iq)
          print *,'qgscrn,scrrel ',qgscrn(iq),scrrel(iq)
       endif
       if(ntest.eq.1.and.tscrn(iq).lt.t(iq,1).and.
     .                 tgg(iq,2).gt.t(iq,1))then
c       if(tstarx.lt.0.)rich=-y*y
         print *,'strange result  iq,land,tgg2,tscrn,tstarx ',
     .       iq,land(iq),tgg(iq,2),tscrn(iq),tstarx
c       print *,'  t1,tstarx,rich,y ',t(iq,1),tstarx,rich,y
       endif
      enddo

      if(jlmspec.eq.1)then ! special treatment of bare soil & veg
         do iq=1,ifull    ! following not checked out lately
          if(land(iq))then
!           bare ground part here
            tstarx = -fgg(iq) / ( rho*cp*ustar(iq))
            c=-grav*tstarx*zscr*af(iq) *sqrt(af(iq) )
     .       /(aft(iq) *tgg(iq,1)*ustar(iq)**2)
            if(tstarx .ge. 0.)then
              rich   = ( -1. + sqrt(1.-4.*bprm*c) ) /(2.*bprm)
              fm   = max(fmroot*fmroot,1./(1.+bprm*rich  )**2)
              con = tstarx*afroot(iq)*min(fmroot,(1.+bprm*rich))/aft(iq)
            else
              y = sqrt(c)
              root=sqrt(zscronzo)
              cm=2.*bprm*cms*af(iq)*root
              ch=2.*bprm*chs*aft(iq)*factch(iq)*root
              fm=1.+2.*bprm*y*y/(1.+cm*y)
              fh=1.+2.*bprm*y*y/(1.+ch*y)
              y = sqrt(c *sqrt(fm*fm*fm)/fh)
              fm=1.+2.*bprm*y*y/(1.+cm*y)
              fh=1.+2.*bprm*y*y/(1.+ch*y)
              y = sqrt(c *sqrt(fm*fm*fm)/fh)
              fm=1.+2.*bprm*y*y/(1.+cm*y)
              fh=1.+2.*bprm*y*y/(1.+ch*y)
              ymin=min(y,ymin)
              con=sqrt(af(iq)*fm)/(aft(iq)*fh)
            endif  ! ( tstarx .ge. 0. )
            delthetg=tstarx*con
            if(ntest.ne.0.and.iq.eq.idjd)then
              print *,'iq,tgg1,delthetg ',iq,tgg(iq,1),delthetg
              if(tstarx.lt.0.)rich=-y*y
              print *,'tstarx,rich,y ',tstarx,rich,y
            endif

!           vegetation part here
            tstarx = -fgf(iq)/(rho*cp*ustar(iq))
            c=-grav*tstarx*zscr*af(iq) *sqrt(af(iq) )
     .         /(aft(iq) *tgf(iq)*ustar(iq)**2)
            if(tstarx .ge. 0.)then
              rich   = ( -1. + sqrt(1.-4.*bprm*c) ) /(2.*bprm)
              fm   = max(fmroot*fmroot,1./(1.+bprm*rich  )**2)
              con = afroot(iq)*min(fmroot,(1.+bprm*rich))/aft(iq)
            else
              y = sqrt(c)
              root=sqrt(zscronzo)
              cm=2.*bprm*cms*af(iq)*root
              ch=2.*bprm*chs*aft(iq)*factch(iq)*root
              fm=1.+2.*bprm*y*y/(1.+cm*y)
              fh=1.+2.*bprm*y*y/(1.+ch*y)
              y = sqrt(c *sqrt(fm*fm*fm)/fh)
              fm=1.+2.*bprm*y*y/(1.+cm*y)
              fh=1.+2.*bprm*y*y/(1.+ch*y)
              y = sqrt(c *sqrt(fm*fm*fm)/fh)
              fm=1.+2.*bprm*y*y/(1.+cm*y)
              fh=1.+2.*bprm*y*y/(1.+ch*y)
              ymin=min(y,ymin)
              con=sqrt(af(iq)*fm)/(aft(iq)*fh)
            endif  ! ( tstarx .ge. 0. )
            delthetf=tstarx*con

            tscrn(iq) = tss(iq)
     .           +delthetf*tsigmf(iq)+delthetg*(1.-tsigmf(iq))
            if(ntest.ne.0.and.iq.eq.idjd)then
              print *,'tgf,delthetf,tscrn ',
     .                 tgf(iq),delthetf,tscrn(iq)
              if(tstarx.lt.0.)rich=-y*y
              print *,'tstarx,rich,y ',tstarx,rich,y
            endif
          endif  ! (land(iq))
         enddo
      endif  ! (jlmspec.eq.1)

      if(ymin.lt.0.)stop ' poor convergence for y in scrnout'
      if(y10min.lt.0.)stop ' poor convergence 10 m wind in scrnout'
      if(ntest.eq.2)then
        print *,' tscreen: ',(tscrn(iq),iq=il/2,ifull,il)
        print *,' qgscrn: ',(qgscrn(iq),iq=il/2,ifull,il)
      endif
      return
      end
