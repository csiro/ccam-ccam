      subroutine scrnout(zo,ustar,factch,wetfac,qsttg,              ! arrays
     .     qgscrn,tscrn,uscrn,u10,scrrel,                           ! arrays
     .     bprm,cms,chs,fmroot,nalpha)
      use cc_mpi, only : mydiag
      use diag_m
!     from Feb '05 scales uscrn, u10 to remove vmod>2 effect
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
      include 'morepbl.h'  ! condx,fg,eg
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
     . ism(ifull),fwtop(ifull),fh(ifull),fm(ifull),rich2(ifull),
     . rich10(ifull),dum3(5*ijk-19*ifull) 
      common/work3b/wblf(ifull,ms),wbfice(ifull,ms),sdepth(ifull,3),
     .          aft(ifull),aft10(ifull),zlog(ifull),
     .          af(ifull),af10(ifull),deltheta(ifull),
     .          afroot(ifull),afroot10(ifull),
     . dum3b(2*ijk-11*ifull-2*ms*ifull)
      common/permsurf/ipsice,ipsea,ipland,iperm(ifull)
      real qgscrn(ifull),tscrn(ifull)
     . ,uscrn(ifull),u10(ifull),scrrel(ifull)
      real wetfac(ifull),factch(ifull),qsttg(ifull),ustar(ifull)
     .     ,zo(ifull)
      real qsurf(ifull)  ! qsurf is just local array, only used here

      include 'establ.h'

c     J. McGregor's vector version, using the formulae of M. Dix

      srcp =sig(1)**(rdry/cp)
      ztv=exp(vkar/sqrt(chn10)) /10.  ! proper inverse of ztsea
c     z3onzt=3.*ztv
c     chn3=(vkar/log(z3onzt))**2
      zscronzt=zscr*ztv
      chnscr=(vkar/log(zscronzt))**2

      if(diag.or.ntest>0)then
        if ( mydiag ) then
         print *,'in scrnout; wetfac,qsttg: ',wetfac(idjd),qsttg(idjd)
         print *,'in scrnout; epot,eg,fg: ',epot(idjd),eg(idjd),fg(idjd)
        end if
        call maxmin(qg,'qg',ktau,1.e3,kl)
        call maxmin(qfg,'qf',ktau,1.e3,kl)
        call maxmin(qlg,'ql',ktau,1.e3,kl)
        call maxmin(t,' t',ktau,1.,kl)
        call maxmin(wb,'wb',ktau,1.,ms)
!       call maxmin(tggsn,'tgg',ktau,1.,ms+3)
        call maxmin(tggsn,'tggsn',ktau,1.,3)
        call maxmin(tgg,'tgg',ktau,1.,ms)
        call maxmin(tss,'ts',ktau,1.,1)
        call maxmin(precip,'pr',ktau,1.,1)
        call maxmin(eg,'eg',ktau,1.,1)
        call maxmin(fg,'fg',ktau,1.,1)
        call maxmin(ps,'ps',ktau,.01,1)
      endif
      
      if(ntest>=2)then
	 do iq=10,ifull,10
	  if(land(iq))then
          zlog(iq)=log(zo(iq))
          afroot(iq)=vkar/log(zmin/zo(iq))                              
          af(iq)=afroot(iq)**2                                             
          aft(iq)=vkar*afroot(iq)/(2. + log(zmin)-zlog(iq))
          if(ntest==3)then
!          extra re-defining of ustar	  
           xx=grav*zmin*(1.-tss(iq)*srcp/t(iq,1))                       
           ri(iq)=xx/vmod(iq)**2                                       
           if(xx.gt.0.)then                                            
             fm(iq)=vmod(iq)*max(fmroot*fmroot,1./(1.+bprm*ri(iq))**2)         
           else                                                             
             root=sqrt(-xx*zmin/zo(iq))  
             denma=vmod(iq)+cms*2.*bprm*af(iq)*root                        
             fm(iq)=vmod(iq)-(2.*bprm *xx)/denma                            
           endif                                                            
c          cduv is now drag coeff *vmod                                     
           cduv(iq) =af(iq)*fm(iq)                                             
           ustar(iq) = sqrt(vmod(iq)*cduv(iq))                             
          endif  ! (ntest==3)
          write (6,"('iq,af,aft,ri,fm,ustar'
     .     ,i5,2f9.5,3f7.3)") iq,af(iq),aft(iq),ri(iq),fm(iq),ustar(iq)
	   
!         now test redefining of fg, to resemble Louis bare soil formula
          factch(iq)=sqrt(7.4)                                           ! land
          theta(iq)=t(iq,1)/srcp   ! can use same value from sflux, as diag
          rho=ps(iq)/(rdry*tss(iq))
          if(ri(iq).gt.0.)then                                             
            fh(iq)=vmod(iq)*max(fmroot*fmroot,1./(1.+bprm*ri(iq))**2)          
          else        
	     xx=ri(iq)*vmod(iq)**2
            root=sqrt(-xx*zmin/zo(iq))                                 
c           Now heat ; allow for smaller zo via aft and factch             
            denha=vmod(iq)+chs*2.*bprm*factch(iq)*aft(iq)*root             
            fh(iq)=vmod(iq)-(2.*bprm *xx)/denha                         
          endif                                                                                                                                   
          conh=rho*aft(iq)*cp                                         
          fg(iq)=conh*fh(iq)*(tss(iq)-theta(iq))                         
	  endif  ! land
	 enddo
      endif  ! (ntest>=2)
         
c     calculate "surface" specific humidity
!     N.B. qsttg is coming in as before-calling-sib3 value of qs_tss
      if(nalpha==1) then  ! usual
        do iq=1,ifull
         qsurf(iq) = wetfac(iq)*qsttg(iq)
     .      + (1.-wetfac(iq))*min(qg(iq,1),qsttg(iq))  ! for v. cold surfaces
        enddo   ! iq=1,ifull
      else
        do iq=1,ifull
         qtgnet=qsttg(iq)*wetfac(iq) -qg(iq,1)
         if(qtgnet>0.) then
           qsurf(iq) = qsttg(iq)*wetfac(iq)
         else
           qsurf(iq) = .1*qsttg(iq)*wetfac(iq) + .9*qg(iq,1)
         endif
        enddo   ! iq loop
      endif   ! (nalpha==1) then .. else ..

      if(ntest.ne.0.and.mydiag)print *,'scrnout zo ',zo(idjd)
c     do iq=1,ifull
c      if(zo(iq)>=zscr.or.zo(iq)<=0.)then
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
      if(ntest==1.and.mydiag)print *,'scrnout zlog ',zlog(idjd)

!cdir nodep
      do ip=ipland+1,ipsea      ! sice and sea points
       iq=iperm(ip)
	if(fracice(iq)<.5)then  ! treat sice as mainly sea here
         aft(iq)=chnscr
         aft10(iq)=chn10
	endif  ! (fracice(iq)<.5)
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
c       if(npanels==0)then
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

        if ( tstarx(iq) >= 0. ) then
c         stable case     (tstarx(iq)>0,   c<0, giving rich>0)
c         Solve quadratic for rich, form ax^2 + bx + c
c         a = bprm ;   b = 1.
c         Because a*c < 0 this has two real solutions, only one of which is
c         positive and so appropriate for the stable case.
c         There is a slight error here because we use tg rather than theta1
c         to calculate c. However this is a second order effect.
          rich2(iq)= ( -1. + sqrt(1.-4.*bprm*c) ) /(2.*bprm)
c         rich3    = ( -1. + sqrt(1.-4.*bprm*c3) ) /(2.*bprm)
          rich10(iq)   = ( -1. + sqrt(1.-4.*bprm*c10) ) /(2.*bprm)

c         Use the positive solution
          fm(iq)  = max(fmroot*fmroot,1./(1.+bprm*rich2(iq)  )**2)
c         fm3  = max(fmroot*fmroot,1./(1.+bprm*rich3 )**2)
          fm10 = max(fmroot*fmroot,1./(1.+bprm*rich10(iq))**2)
          con = afroot(iq)*min(fmroot,(1.+bprm*rich2(iq)))/aft(iq)

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
          if(ntest.ne.0.and.iq==idjd.and.mydiag) then
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
	   rich2(iq)=-y*y

c         iterations for 10 m winds
          y10 = sqrt(c10)
          root=sqrt(10./zo(iq))
          cm10=2.*bprm*cms*af10(iq)*root
          ch10=2.*bprm*chs*aft10(iq)*factch(iq)*root
          fm10=1.+2.*bprm*y10*y10/(1.+cm10*y10)
          fh10=1.+2.*bprm*y10*y10/(1.+ch10*y10)
          y10 = sqrt(c10 *sqrt(fm10*fm10*fm10)/fh10)  ! only 2 itns
          fm10=1.+2.*bprm*y10*y10/(1.+cm10*y10)
          y10min=min(y10,y10min)
	   rich10(iq)=-y10*y10
        endif

c       calculate correction to surface temperature deltheta and moisture
        deltheta(iq)=tstarx(iq)*con
        deltq = qstarx(iq)*con

c       screen wind speeds
	 vfact=sqrt(u(iq,1)**2+v(iq,1)**2)/vmod(iq)
        uscrn(iq) = vfact*ustar(iq)/(afroot(iq)*sqrt(fm(iq)))
c       u3(iq) = ustar(iq)/(afroot3(iq)*sqrt(fm3))
        u10(iq) = vfact*ustar(iq)/(afroot10(iq)*sqrt(fm10))

c       screen specific humidity
c       N.B. over unstable sea points, sometimes will get supersat.qgscrn
        qgscrn(iq) = qsurf(iq) + deltq   ! will only use for sea

        if(ntest.ne.0.and.iq==idjd.and.mydiag)then
          print *,' scrnout test for id,jd,idjd,land,sice ',
     .                               id,jd,idjd,land(iq),sice(iq)
          print *,'tss,theta,ps,con ',tss(iq),theta(iq),ps(iq),con
          print *,'nalpha,fh,fm,fm10 ',nalpha,fh(iq),fm(iq),fm10
	   print *,'deltq,qsurf,qgscrn,qg1 ',
     .             deltq,qsurf(iq),qgscrn(iq),qg(iq,1)
          print *,'c,rich,zscr ',c,rich2(iq),zscr
        endif ! (ntest.ne.0.and.iq==idjd)
      enddo   ! iq=1,ifull

      if(mydiag.and.(diag.or.ntest==2))then  ! special test of 38 m, stable sea
       iq=idjd
       zscronzt38=zmin*ztv  ! zmin is at 38 m
       chnscr38=(vkar/log(zscronzt38))**2
       afroot38=vkar/(log(zmin)-zlog(iq))
       af38=afroot38*afroot38
       aft38=chnscr38
       rho=ps(iq)/(rdry*tss(iq))
!      tstarx(iq) = -fg(iq) / ( rho*cp*ustar(iq))
!      qstarx(iq) = -eg(iq) / ( rho*hl*ustar(iq))
       c38 = -9.806*tstarx(iq)*zmin*af38*sqrt(af38)
     .                   /(aft38 *tss(iq)*ustar(iq)**2)
!      if ( tstarx(iq) >= 0. ) then
!        rich38   = ( -1. + sqrt(1.-4.*bprm*c38) ) /(2.*bprm)
         ffm   = max(fmroot*fmroot,1./(1.+bprm*rich38  )**2)
c        con = afroot38*min(fmroot,(1.+bprm*rich2(iq)))/aft38
         con = afroot38*min(1./fmroot,1.+bprm*rich38)/aft38
!       deltheta(iq)=tstarx(iq)*con
        deltq = qstarx(iq)*con
!       uscrn(iq) = ustar(iq)/(afroot38*sqrt(ffm))
	 print *,'38 m diagnostics'
	 print *,'zscr,af38,afroot38 ',zmin,af38,afroot38
	 print *,'fmroot,fm,con ',fmroot,ffm,con
	 print *,'chnscr38,zo,zt,ztv ',chnscr38,zo(iq),1./ztv,ztv
	 print *,'aft38,tstarx(iq),rich38 ',aft38,tstarx(iq),rich38
	 print *,'tgg1,tgg2,thscrn.theta1 ',
     .           tgg(iq,1),tgg(iq,2),tgg(iq,2)+deltheta(iq),theta(iq)
        print *,'uscrn,uv1 ',uscrn(iq),sqrt(u(iq,1)**2+v(iq,1)**2)
      endif  ! (diag.or.ntest==2)

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
       if(ntest==3)then  ! print all stable sea points
	  if(fg(iq)<0.)then
	    print *,'iq,tgg2,tscrn,theta,ri,rich s ',
     .              iq,tgg(iq,2),tscrn(iq),theta(iq),ri(iq),rich2(iq)
           print *,'fhx,fm ',fhx(iq)/vmod(iq),fm(iq)
	  endif
       endif               ! (ntest==3)
       if(ntest==4)then  ! print all unstable sea points
	  if(fg(iq)>0.)then
           zmin=zmin
           zscronzt38=zmin*ztv
           chnscr38=(vkar/log(zscronzt38))**2
           afroot38=vkar/(log(zmin)-zlog(iq))
           af38=afroot38*afroot38
           aft38=chnscr38
           c38 = -9.806*tstarx(iq)*zmin*af38*sqrt(af38)
     .                   /(aft38 *tss(iq)*ustar(iq)**2)
!          calculate an approximate form for ri38
           ri38=-c38*vmod(iq)*fm(iq)*sqrt(fm(iq))/fhx(iq)
	    print *,'iq,tgg2,tscrn,theta,ri,rich u ',
     .              iq,tgg(iq,2),tscrn(iq),theta(iq),ri(iq),rich2(iq)
           print *,'fhx,fh,fm,ri38 ',fhx(iq)/vmod(iq),fh(iq),fm(iq),ri38
	  endif
       endif              ! (ntest==4)
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
        print *,'rich,zscr ',rich2(iq),zscr
        es = establ(tscrn(iq))
        constz=ps(iq)-es
        qs= 0.622*es/constz
        print *,'tss,tscrn,theta,uscrn,u10 ',
     .           tss(iq),tscrn(iq),theta(iq),uscrn(iq),u10(iq)
        print *,'deltheta,qs ',deltheta(iq),qs
      endif ! (diag.or.ntest.ne.0)

      if(jlmspec==1)then ! special treatment of bare soil & veg
         do iq=1,ifull    ! following not checked out lately
          if(land(iq))then
!           bare ground part here
            tstarx(iq) = -fgg(iq) / ( rho*cp*ustar(iq))
            c=-9.806*tstarx(iq)*zscr*af(iq) *sqrt(af(iq) )
     .       /(aft(iq) *tgg(iq,1)*ustar(iq)**2)
            if(tstarx(iq) >= 0.)then
              rich2(iq) = ( -1. + sqrt(1.-4.*bprm*c) ) /(2.*bprm)
              fm(iq) = max(fmroot*fmroot,1./(1.+bprm*rich2(iq) )**2)
              con = tstarx(iq)*afroot(iq)*
     .                      min(fmroot,(1.+bprm*rich2(iq)))/aft(iq)
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
		rich2(iq)=-y*y
            endif  ! ( tstarx(iq) >= 0. )
            delthetg=tstarx(iq)*con
            if(mydiag.and.ntest.ne.0.and.iq==idjd)then
              print *,'iq,tgg1,delthetg ',iq,tgg(iq,1),delthetg
              print *,'tstarx(iq),rich,y ',tstarx(iq),rich2(iq),y
            endif

!           vegetation part here
            tstarx(iq) = -fgf(iq)/(rho*cp*ustar(iq))
            c=-9.806*tstarx(iq)*zscr*af(iq) *sqrt(af(iq) )
     .         /(aft(iq) *tgf(iq)*ustar(iq)**2)
            if(tstarx(iq) >= 0.)then
              rich2(iq)   = ( -1. + sqrt(1.-4.*bprm*c) ) /(2.*bprm)
              fm(iq) = max(fmroot*fmroot,1./(1.+bprm*rich2(iq)  )**2)
              con = afroot(iq)*min(fmroot,(1.+bprm*rich2(iq)))/aft(iq)
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
		rich2(iq)=-y*y
            endif  ! ( tstarx(iq) >= 0. )
            delthetf=tstarx(iq)*con

            tscrn(iq) = tss(iq)
     .           +delthetf*tsigmf(iq)+delthetg*(1.-tsigmf(iq))
            if(mydiag.and.ntest.ne.0.and.iq==idjd)then
              print *,'tgf,delthetf,tscrn ',
     .                 tgf(iq),delthetf,tscrn(iq)
              print *,'tstarx(iq),rich,y ',tstarx(iq),rich2(iq),y
            endif
          endif  ! (land(iq))
         enddo
      endif  ! (jlmspec==1)

      if(ntest>=1)then
        numoddp=0
        numoddm=0
        numt=0
        numq=0
        numu=0
!       check that fg always consistent with T gradient
	 do iq=1,ifull
	  if(tss(iq)<theta(iq).and.fg(iq)>0.)then
           numoddp=numoddp+1
           write (6,"('oddp T iq,isoil,iveg,land,sice,tss,theta,'
     .      'snowd,tsigmf,fgf,fgg,fg',i5,2i3,2l2,3f8.2,4f7.2)") 
     .       iq,isoilm(iq),ivegt(iq),land(iq),sice(iq),tss(iq),
     .       theta(iq),snowd(iq),tsigmf(iq),fgf(iq),fgg(iq),fg(iq)
         endif	  
	  if(tss(iq)>theta(iq).and.fg(iq)<0.)then
           numoddm=numoddm+1
           write (6,"('oddm T iq,isoil,iveg,land,sice,tss,theta,'
     .      'snowd,tsigmf,fgf,fgg,fg',i5,2i3,2l2,3f8.2,4f7.2)") 
     .       iq,isoilm(iq),ivegt(iq),land(iq),sice(iq),tss(iq),
     .       theta(iq),snowd(iq),tsigmf(iq),fgf(iq),fgg(iq),fg(iq)
         endif	  
        enddo
	 do iq=1,ifull
	  vmag=sqrt(u(iq,1)**2+v(iq,1)**2)
	  if(uscrn(iq)>vmag.or.u10(iq)>vmag)then
           numu=numu+1
           write (6,"('strange u iq,land,uscrn,u10,vmag,vmod,'
     .     'tstarx,fm ',i5,l2,6f7.3)") iq,land(iq),
     .     uscrn(iq),u10(iq),vmag,vmod(iq),tstarx(iq),fm(iq)
           write (6,"('more u isoil,iveg,snowd,tsigmf,fgf,fgg,fg',	     
     .     2i3,6f8.2)") isoilm(iq),ivegt(iq),snowd(iq),tsigmf(iq),
     .     fgf(iq),fgg(iq),fg(iq)
           print *,'tgg1,tggsn1,tss,theta ',
     .              tgg(iq,1),tggsn(iq,1),tss(iq),theta(iq)
         endif
	  if(uscrn(iq)>vmag.or.u10(iq)>vmag.or.mod(iq,10)==0)then
	   if(land(iq))then
!          following diags done as if land	     
           afroot38=vkar/(log(zmin)-zlog(iq))
           af38=afroot38*afroot38
	    aft38=vkar*afroot38/(2. + log(zmin)-zlog(iq))
           c38 = -9.806*tstarx(iq)*zmin*af38*sqrt(af38)
     .                    /(aft38 *tss(iq)*ustar(iq)**2)
           if(tstarx(iq)>0.)then
            rich38   = ( -1. + sqrt(1.-4.*bprm*c38) ) /(2.*bprm)
            fm38   = max(fmroot*fmroot,1./(1.+bprm*rich38  )**2)
            con = afroot38*min(1./fmroot,1.+bprm*rich38)/aft38
c************* new stable	     
            fm38=max(fmroot*fmroot,1./(1.+bprm*ri(iq))**2)
	     fh38=fm38
	     ri2a=ri(iq)*zscr/zmin
            fm2=max(fmroot*fmroot,1./(1.+bprm*ri2a)**2)
	     fh2=fm2
	     ri10a=ri(iq)*10./zmin
            fm10=max(fmroot*fmroot,1./(1.+bprm*ri10a)**2)
	     fh10=fm10
	    else
            root=sqrt(zmin/zo(iq))
            cm38=2.*bprm*cms*af38*root
	     ria=-cm38
            y38 = sqrt(c38)                             ! y0
            ch38=2.*bprm*chs*aft38*factch(iq)*root
            fm38=1.+2.*bprm*y38*y38/(1.+cm38*y38)
            fh38=1.+2.*bprm*y38*y38/(1.+ch38*y38)
            y38 = sqrt(c38 *sqrt(fm38*fm38*fm38)/fh38)  ! y1
	     rib=-y38*y38
            fm38=1.+2.*bprm*y38*y38/(1.+cm38*y38)
            fh38=1.+2.*bprm*y38*y38/(1.+ch38*y38)
            y38 = sqrt(c38 *sqrt(fm38*fm38*fm38)/fh38)  ! y2 
	     rich38=-y38*y38            ! y3  final one
            fm38=1.+2.*bprm*y38*y38/(1.+cm38*y38)
            fh38=1.+2.*bprm*y38*y38/(1.+ch38*y38)
            con=sqrt(af38*fm38)/(aft38*fh38)
            y38 = sqrt(c38 *sqrt(fm38*fm38*fm38)/fh38)  ! y2 
	     rid=-y38*y38              ! y4  extra for checking
c************* new unstable	     
            root=sqrt(-ri(iq)*zmin/zo(iq))        
            denma=1.+cms*2.*bprm*af38*root                         
            fm38=1.-2.*bprm *ri(iq)/denma                           
            denha=1.+chs*2.*bprm*factch(iq)*aft38*root             
            fh38=1.-2.*bprm *ri(iq)/denha                        
	     ri2a=ri(iq)*zscr/zmin
            root=sqrt(-ri2a*zscr/zo(iq))        
            denma=1.+cms*2.*bprm*af(iq)*root                         
            fm2=1.-2.*bprm *ri2a/denma                           
            denha=1.+chs*2.*bprm*factch(iq)*aft(iq)*root             
            fh2=1.-2.*bprm *ri2a/denha          
	     ri10a=ri(iq)*10./zmin
            root=sqrt(-ri10a*10./zo(iq))        
            denma=1.+cms*2.*bprm*af10(iq)*root                         
            fm10=1.-2.*bprm *ri10a/denma                           
            denha=1.+chs*2.*bprm*factch(iq)*aft10(iq)*root             
            fh10=1.-2.*bprm *ri10a/denha                        
          endif
           u38 = ustar(iq)/(afroot38*sqrt(fm38))
           write(6,"('iq,uscrn,u38,uv1,ri2,ri38,ri,ria,rib,rid ',
     .       i5,3f7.2,3f7.3,5x,3f7.3)") 
     .         iq,uscrn(iq),u38,vmag,rich2(iq),rich38,ri(iq),ria,rib,rid
           write(6,"('iq,tss,tscrn,t1,theta,ri2,ri38,ri ',
     .       i5,4f8.3,5x,3f7.3)") iq,tss(iq),tscrn(iq),
     .         t(iq,1),theta(iq),rich2(iq),rich38,ri(iq)
           ri2x=ri2a*sqrt(af(iq)*fm2/(af38*fm38))*
     .               af(iq)*fm2*aft38*fh38/(af38*fm38*aft(iq)*fh2)      	    
           ri10x=ri10a*sqrt(af10(iq)*fm10/(af38*fm38))*
     .               af10(iq)*fm10*aft38*fh38/(af38*fm38*aft10(iq)*fh10)      	    
           write(6,"('iq,ri2,ri10,ri38,ri,ri2x,ri10x ',i5,6f7.3)") 
     .         iq,rich2(iq),rich10(iq),rich38,ri(iq),ri2x,ri10x
          endif
	  endif
	  if(land(iq).or.sice(iq))then
           if(tscrn(iq)<min(theta(iq),tss(iq)).or.
     .        tscrn(iq)>max(theta(iq),tss(iq))    )then
             numt=numt+1
             write (6,"('strange T iq,land,sice,tss,tscrn,tgg1,t1,'
     .       'theta,ri2,fh ',i5,2l2,5f7.2,2f7.3)") 
     .        iq,land(iq),sice(iq),tss(iq),tscrn(iq),
     .           tgg(iq,1),t(iq,1),theta(iq),rich2(iq),fh(iq)
             if(sice(iq))print *,'tgg2,tgg3,fracice ',
     .                            tgg(iq,2),tgg(iq,3),fracice(iq)  	      
           endif
         else
           if(tscrn(iq)<min(theta(iq),tgg(iq,2)).or.
     .        tscrn(iq)>max(theta(iq),tgg(iq,2))    )then
             numt=numt+1
             write (6,"('strange T iq,land,sice,tss,tscrn,tgg1,t1,'
     .       'theta,ri2,fh ',i5,2l2,5f7.2,2f7.3)") 
     .        iq,land(iq),sice(iq),tss(iq),tscrn(iq),
     .           tgg(iq,1),t(iq,1),theta(iq),rich2(iq),fh(iq)
           endif
         endif   ! (land(iq).or.sice(iq)) .. else ..
	  if(qgscrn(iq)<min(qg(iq,1),qsurf(iq)).or.
     .      qgscrn(iq)>max(qg(iq,1),qsurf(iq))    )then
           numq=numq+1
           print *,'strange qg iq,land,sice,qsurf,qgscrn ',
     .      iq,land(iq),sice(iq),qsurf(iq),qgscrn(iq)
           print *,'qg1,t1,ri ',qg(iq,1),t(iq,1),ri(iq),numq
         endif
	 enddo
	 print *,'numq,numt,numu,numoddp,numoddm',
     .           numq,numt,numu,numoddp,numoddm
      endif   ! (ntest>=1)

!     consistency fixer for tscrn and qgscrn (odd pts)  March 2004
!     N.b. T & qg problems are for land and sea-ice points
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

      if(ymin<0.)stop ' poor convergence for y in scrnout'
      if(y10min<0.)stop ' poor convergence 10 m wind in scrnout'
      if(ntest==2)then
        print *,' tscreen: ',(tscrn(iq),iq=il/2,ifull,il)
        print *,' qgscrn: ',(qgscrn(iq),iq=il/2,ifull,il)
      endif
      return
      end
