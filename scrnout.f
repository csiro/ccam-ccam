      subroutine scrnout(zo,ustar,factch,wetfac,qsttg,              ! arrays
     .     qgscrn,tscrn,uscrn,u10,scrrel,af,aft,ri,vmod,            ! arrays
     .     bprm,cms,chs,chnsea,nalpha)
      use cc_mpi, only : mydiag
      use diag_m
      implicit none
!     from Feb '05 scales uscrn, u10 to remove vmod>2 effect
!     allows for zo>zscr from 30/7/04
      integer nits,ntest
      real vkar,zscr
!     has fixer at bottom to ensure tscrn and qgscrn bounded (& non-neg)
      parameter (nits=2)    ! nits=2 for 2 iterations
      parameter (ntest=0)   ! ntest= 0 for diags off; ntest>= 1 for diags on
!     parameter (chn10=.00136733)   ! sea only, same as in sflux via parm.h
      parameter (vkar=.4,zscr=1.8)
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
      real ri(ifull),vmod(ifull)
      real af(ifull),aft(ifull)  ! passed in
      integer ipsice,ipsea,ipland,iperm
      common/permsurf/ipsice,ipsea,ipland,iperm(ifull)
      real qgscrn(ifull),tscrn(ifull)
     . ,uscrn(ifull),u10(ifull),scrrel(ifull)
      real wetfac(ifull),factch(ifull),qsttg(ifull),ustar(ifull)
     .     ,zo(ifull),deltaq,deltat,theta(ifull)
      real qsurf(ifull),tsurf(ifull)  ! qsurf is just local array
      real af2(ifull),af10(ifull),afroot2(ifull),afroot10(ifull)
      real aft2(ifull),aft10(ifull),rich2(ifull),rich10(ifull)
     .     ,qstarx(ifull),tstarx(ifull)
      integer ip,iq,nalpha,numq,numt,numu
      real alf,chnscr,chnsea,fact,qtgnet
      real rho,srcp,vfact,zscronzt,zscrr,ztv,z10
      real fh2(ifull),fh10(ifull),fh38(ifull)
      real fm2(ifull),fm10(ifull),fm38(ifull)
      real aft2ice,aft10ice,aft38ice,zt,zscrt,z10t
      real bprm,cms,chs,denha,denma,ri2x,root,vmag,zlog
      include 'establ.h'

      srcp =sig(1)**(rdry/cp)
      ztv=exp(vkar/sqrt(chn10)) /10.  ! proper inverse of ztsea
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

      do iq=1,ifull
       zscrr=max(zscr,zo(iq)+1.)
	z10=max(10.,zo(iq)+2.)
       zlog=log(zo(iq))
       afroot2(iq)=vkar/(log(zscrr)-zlog)
       af2(iq)=afroot2(iq)*afroot2(iq)
c      constants for calculating 2m, 10m winds
       afroot10(iq) = vkar/(log(z10)-zlog)
       af10(iq) = afroot10(iq)*afroot10(iq)
	zt=zo(iq)/7.4
	zscrt=max(zscr,zt+1.)     ! Aug '05
	z10t=max(z10,zt+2.)       ! Aug '05
       aft2(iq)=vkar*afroot2(iq)/log(zscrt/zt)     ! land only
       aft10(iq)=vkar*afroot10(iq)/log(z10t/zt)    ! land only
c      aft2(iq)=vkar*afroot2(iq)/(2. + log(zscrr)-zlog)    ! land only
c      aft10(iq)=vkar* afroot10(iq) /(2. + log(z10)-zlog)  ! land only 
	tsurf(iq)=tss(iq)
      enddo   ! iq loop

       aft2ice=(vkar/log(zscrr/.001))**2
       aft10ice=(vkar/log(z10/.001))**2
       aft38ice=(vkar/log(zmin/.001))**2  ! needs recalc, as aft is ice only
!cdir nodep
      do ip=ipland+1,ipsea      ! sice and sea points
       iq=iperm(ip)
       tsurf(iq)=fracice(iq)*tgg(iq,3)+(1.-fracice(iq))*tgg(iq,2) 
       aft2(iq)=fracice(iq)*aft2ice+(1.-fracice(iq))*chnscr
       aft10(iq)=fracice(iq)*aft10ice+(1.-fracice(iq))*chn10
       aft(iq)=fracice(iq)*aft38ice+(1.-fracice(iq))*chnsea
      enddo

!     the following assumes a suitable ri(level 1) has been found
      do iq=1,ifull              ! all valid points 
       zscrr=max(zscr,zo(iq)+1.)
	z10=max(10.,zo(iq)+2.)
       theta(iq)=t(iq,1)/srcp   ! can use same value from sflux
       deltat=tsurf(iq)-theta(iq)
	deltaq=qsurf(iq)-qg(iq,1)
       rho=ps(iq)/(rdry*tsurf(iq))
c      section for solving for screen height values of various parameters
       if (ri(iq) > 0.) then
c        stable case     
         fm38(iq)=1./(1.+bprm*ri(iq))**2
	  fh38(iq)=fm38(iq)
c	  rich2(iq)=ri(iq)*zscrr/zmin        ! usually a good first guess
	  fact=af2(iq)/(af(iq)*fm38(iq))
  	  alf=ri(iq)*zscrr*fact*sqrt(fact)*    
cc     .             aft(iq)*fh38(iq)/(zmin*aft2(iq))
     .             aft(iq)*fh38(iq)*sqrt(fh38(iq))/(zmin*aft2(iq))
         rich2(iq)=( -1. + sqrt(1.+4.*bprm*alf) ) /(2.*bprm)
         fm2(iq)=1./(1.+bprm*rich2(iq))**2
	  fh2(iq)=fm2(iq)
	  fact=af10(iq)/(af(iq)*fm38(iq))
  	  alf=ri(iq)*z10*fact*sqrt(fact)*    
     .             aft(iq)*fh38(iq)*sqrt(fh38(iq))/(zmin*aft10(iq))
         rich10(iq)=( -1. + sqrt(1.+4.*bprm*alf) ) /(2.*bprm)
         fm10(iq)=1./(1.+bprm*rich10(iq))**2
         fh10(iq)=fm10(iq)
c        alternative formulae - not checked out	  
c         tstarx(iq)=aft(iq)*vmod(iq)*fh38(iq)*deltat/ustar(iq)
c  	  alf=grav*zscrr*tstarx(iq)*af2(iq)*sqrt(af2(iq))/  
c     .             (theta(iq)*aft2(iq)*ustar(iq)**2)
c         rich2(iq)=( -1. + sqrt(1.-4.*bprm*alf) ) /(2.*bprm)
c         fm2(iq)=1./(1.+bprm*rich2(iq))**2
c	  fh2(iq)=fm2(iq)
c  	  alf=grav*z10*tstarx(iq)*af10(iq)*sqrt(af10(iq))/  
c     .             (theta(iq)*aft10(iq)*ustar(iq)**2)
c         rich10(iq)=( -1. + sqrt(1.-4.*bprm*alf) ) /(2.*bprm)
c         fm10(iq)=1./(1.+bprm*rich10(iq))**2
c	  fh10(iq)=fm10(iq)
       else
c        unstable case
c        Starting value of the Richardson number.
         root=sqrt(-ri(iq)*zmin/zo(iq))        
         denma=1.+cms*2.*bprm*af(iq)*root                         
         fm38(iq)=1.-2.*bprm *ri(iq)/denma                           
         denha=1.+chs*2.*bprm*factch(iq)*aft(iq)*root             
         fh38(iq)=1.-2.*bprm *ri(iq)/denha                        
	  rich2(iq)=ri(iq)*zscrr/zmin        ! usually a good first guess
         root=sqrt(-rich2(iq)*zscrr/zo(iq))        
         denma=1.+cms*2.*bprm*af2(iq)*root                         
         fm2(iq)=1.-2.*bprm *rich2(iq)/denma                           
         denha=1.+chs*2.*bprm*factch(iq)*aft2(iq)*root             
         fh2(iq)=1.-2.*bprm *rich2(iq)/denha          
	  rich10(iq)=ri(iq)*z10/zmin        ! usually a good first guess
         root=sqrt(-rich10(iq)*z10/zo(iq))        
         denma=1.+cms*2.*bprm*af10(iq)*root      
         fm10(iq)=1.-2.*bprm *rich10(iq)/denma                           
         denha=1.+chs*2.*bprm*factch(iq)*aft10(iq)*root             
         fh10(iq)=1.-2.*bprm *rich10(iq)/denha                        
	  if(nits==2)then
	    fact=af2(iq)*fm2(iq)/(af(iq)*fm38(iq))
  	    rich2(iq)=ri(iq)*zscrr*fact*sqrt(fact)*    
     .               aft(iq)*fh38(iq)/(zmin*aft2(iq)*fh2(iq))
           root=sqrt(-rich2(iq)*zscrr/zo(iq))        
           denma=1.+cms*2.*bprm*af2(iq)*root                         
           fm2(iq)=1.-2.*bprm *rich2(iq)/denma                           
           denha=1.+chs*2.*bprm*factch(iq)*aft2(iq)*root             
           fh2(iq)=1.-2.*bprm *rich2(iq)/denha          
	    fact=af10(iq)*fm10(iq)/(af(iq)*fm38(iq))
  	    rich10(iq)=ri(iq)*z10*fact*sqrt(fact)*    
     .               aft(iq)*fh38(iq)/(zmin*aft10(iq)*fh10(iq))
           root=sqrt(-rich10(iq)*z10/zo(iq))        
           denma=1.+cms*2.*bprm*af10(iq)*root      
           fm10(iq)=1.-2.*bprm *rich10(iq)/denma                           
           denha=1.+chs*2.*bprm*factch(iq)*aft10(iq)*root             
           fh10(iq)=1.-2.*bprm *rich10(iq)/denha                        
	  endif  ! (nits==2)
       endif
	
!      in calculating tstar & qstar, use derived ri38, NOT eg and fg	
       tstarx(iq)=aft(iq)*vmod(iq)*fh38(iq)*deltat/ustar(iq)
       qstarx(iq)=aft(iq)*vmod(iq)*fh38(iq)*deltaq/ustar(iq)
       fact=sqrt(af2(iq)*fm2(iq))/(aft2(iq)*fh2(iq))
       tscrn(iq)=tsurf(iq)-fact*tstarx(iq)
c      N.B. over unstable sea points, may sometimes get supersat qgscrn
       qgscrn(iq)=qsurf(iq)-fact*qstarx(iq)
     
c      screen wind speeds
       vfact=sqrt(u(iq,1)**2+v(iq,1)**2)/vmod(iq)
       uscrn(iq) = vfact*ustar(iq)/(afroot2(iq)*sqrt(fm2(iq)))
       u10(iq) = vfact*ustar(iq)/(afroot10(iq)*sqrt(fm10(iq)))
       scrrel(iq) = 100.*qgscrn(iq)/qsttg(iq)
      enddo

      if(ntest>=1)then
        numt=0
        numq=0
        numu=0
        do iq=1,ifull
         zscrr=max(zscr,zo(iq)+1.)
	  z10=max(10.,zo(iq)+2.)
         vmag=sqrt(u(iq,1)**2+v(iq,1)**2)
         if(uscrn(iq)>vmag.or.u10(iq)>vmag)then
           numu=numu+1
           write (6,"('strange u iq,land,uscrn,u10,vmag,vmod,'
     .     'tstarx ',i5,l2,6f7.3)") iq,land(iq),
     .     uscrn(iq),u10(iq),vmag,vmod(iq),tstarx(iq)
           write (6,"('more u isoil,iveg,snowd,tsigmf,fg',	     
     .     2i3,6f8.2)") isoilm(iq),ivegt(iq),snowd(iq),tsigmf(iq),fg(iq)
           print *,'tgg1,tggsn1,tss,theta ',
     .              tgg(iq,1),tggsn(iq,1),tss(iq),theta(iq)
         endif
         if(tscrn(iq)<min(theta(iq),tsurf(iq)).or.
     .      tscrn(iq)>max(theta(iq),tsurf(iq)).or.iq==idjd )then
c        if(tsurf(iq)<theta(iq))then
            numt=numt+1
            write (6,"('checking T iq,land,sice,tsurf,tscrn,theta,'
     .       'tstarx ',i5,2l2,3f7.2,f7.3)") iq,land(iq),sice(iq), 
     .        tsurf(iq),tscrn(iq),theta(iq),tstarx(iq)
            write (6,"('tss,tgg1,tgg2,tgg3,t1,fracice ',6f7.2)")
     .         tss(iq),tgg(iq,1),tgg(iq,2),tgg(iq,3),t(iq,1),fracice(iq)
            print *,'tgf,tsigmf,zo ',tgf(iq),tsigmf(iq),zo(iq)
            print *,'ri2,ri10,ri38 ',rich2(iq),rich10(iq),ri(iq)
	     print *,'cm2,ch2,cm10,ch10 ',
     .               cms*2.*bprm*af2(iq)*sqrt(zscrr/zo(iq)),
     .               cms*2.*bprm*af10(iq)*sqrt(z10/zo(iq)),
     .               chs*2.*bprm*aft2(iq)*factch(iq)*sqrt(zscrr/zo(iq)),
     .               chs*2.*bprm*aft10(iq)*factch(iq)*sqrt(z10/zo(iq))     
            print *,'af2,af10,af38 ',af2(iq),af10(iq),af(iq)
            print *,'aft2,aft10,aft38 ',aft2(iq),aft10(iq),aft(iq)
	     print *,'rat2,rat10,rat38 ',sqrt(af2(iq))/aft2(iq),
     &               sqrt(af10(iq))/aft10(iq),sqrt(af(iq))/aft(iq)	    
            print *,'fh2,fh10,fh38 ',fh2(iq),fh10(iq),fh38(iq)
            print *,'fm2,fm10,fm38 ',fm2(iq),fm10(iq),fm38(iq)
	     print *,'vmod,fact2,fact10,fact38 ',vmod(iq),
     .               sqrt(af2(iq)*fm2(iq))/(aft2(iq)*fh2(iq)),
     .               sqrt(af10(iq)*fm10(iq))/(aft10(iq)*fh10(iq)),
     .               sqrt(af(iq)*fm38(iq))/(aft(iq)*fh38(iq))
            if(sice(iq))print *,'tgg2,tgg3,fracice ',
     .                           tgg(iq,2),tgg(iq,3),fracice(iq)  	      
           endif
	 if(qgscrn(iq)<min(qg(iq,1),qsurf(iq)).or.
     .     qgscrn(iq)>max(qg(iq,1),qsurf(iq))    )then
           numq=numq+1
           print *,'strange qg iq,land,sice,qsurf,qgscrn,qstarx ',
     .      iq,land(iq),sice(iq),qsurf(iq),qgscrn(iq),qstarx(iq)
           print *,'qg1,t1,ri ',qg(iq,1),t(iq,1),ri(iq),numq
         endif
	enddo
	print *,'numq,numt,numu',numq,numt,numu
	if(ntest==2)then
	 do iq=1,ifull,50
!	  prints first guess and final values
	  ri2x=rich2(iq)
	  fact=af2(iq)*fm2(iq)/(af(iq)*fm38(iq))
         if(ri(iq)<0.)ri2x=ri(iq)*zscrr*fact*sqrt(fact)*    
     .              aft(iq)*fh38(iq)/(zmin*aft2(iq)*fh2(iq))
         write (6,"('iq,ri2a,ri2,ri2x,ri10a,ri10 ',i5,5f9.5)")
     .    iq,zscrr*ri(iq)/zmin,rich2(iq),ri2x,z10*ri(iq)/zmin,rich10(iq)
	 enddo
	endif
      endif   ! (ntest>=1)      

!     consistency fixer for tscrn and qgscrn (odd pts)  March 2004
!     not needed from March 2005 with new version of scrnout
!     N.B. T & qg problems are for land and sea-ice points
c      do iq=1,ifull
c       tscrn(iq)=max(tscrn(iq),min(tsurf(iq),theta(iq)))
c       tscrn(iq)=min(tscrn(iq),max(tsurf(iq),theta(iq)))
c       qgscrn(iq)=max(qgscrn(iq),min(qsurf(iq),qg(iq,1)))
c       qgscrn(iq)=min(qgscrn(iq),max(qsurf(iq),qg(iq,1)))
c      enddo
   
      tscrn(:)=tscrn(:)-.018  ! apply 1.8 m approx. adiab. corrn.
      return
      end
