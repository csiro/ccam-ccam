      subroutine scrnout(zo,ustar,factch,wetfac,qsttg,              ! arrays
     .     qgscrn,tscrn,uscrn,u10,rhscrn,af,aft,ri,vmod,            ! arrays
     .     bprm,cms,chs,chnsea,nalpha)
      use arrays_m
      use cc_mpi, only : mydiag,myid
      use diag_m
      use liqwpar_m  ! ifullw,qfg,qlg  just for diags
      use morepbl_m  ! condx,fg,eg
      use nsibd_m    ! rsmin,ivegt,sigmf,tgf,ssdn,res,rmc,tsigmf
      use pbl_m
      use permsurf_m
      use prec_m     ! just for diags
      use sigs_m
      use soil_m
      use soilsnow_m
      implicit none
!     from Feb '05 scales uscrn, u10 to remove vmod>2 effect
!     allows for zo>zscr from 30/7/04
!     has fixer at bottom to ensure tscrn and qgscrn bounded (& non-neg)
      integer, parameter :: nits=2  ! nits=2 for 2 iterations
      integer, parameter :: ntest=0 ! ntest= 0 for diags off; ntest>= 1 for diags 
!     parameter (chn10=.00136733)   ! sea only, same as in sflux via parm.h
      real, parameter :: vkar=.4, zscr=1.8
      include 'newmpar.h'
      include 'const_phys.h'
      include 'parm.h'
      real ri(ifull),vmod(ifull)
      real af(ifull),aft(ifull)  ! passed in
      real qgscrn(ifull),tscrn(ifull)
     . ,uscrn(ifull),u10(ifull),rhscrn(ifull)
      real wetfac(ifull),factch(ifull),qsttg(ifull),ustar(ifull)
     .     ,zo(ifull),deltaq,deltat,theta(ifull)
      real qsurf(ifull),tsurf(ifull)  ! qsurf is just local array
      real af2(ifull),af10(ifull),afroot2(ifull),afroot10(ifull)
      real aft2(ifull),aft10(ifull),rich2(ifull),rich10(ifull)
     .     ,qstarx(ifull),tstarx(ifull)
      integer iq,nalpha,numq,numt,numu
      real alf,chnscr,chnsea,fact,qtgnet
      real rho,srcp,vfact,zscronzt,zscrr,ztv,z10
      real fh2(ifull),fh10(ifull),fh38(ifull)
      real fm2(ifull),fm10(ifull),fm38(ifull)
      real aft2ice,aft10ice,aft38ice,zt,zscrt,z10t
      real bprm,cms,chs,denha,denma,ri2x,root,vmag,zlog,es
      include 'establ.h'

      srcp =sig(1)**(rdry/cp)
      ztv=exp(vkar/sqrt(chn10)) /10.  ! proper inverse of ztsea
      zscronzt=zscr*ztv
      chnscr=(vkar/log(zscronzt))**2

      if(diag.or.ntest>0)then
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
        if(mydiag)then
        print *,'entering scrnout; wetfac,tss,qsttg: ',
     &                             wetfac(idjd),tss(idjd),qsttg(idjd)
        endif
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
       zt=zo(iq)/factch(iq)**2 ! MJT urban
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
      aft38ice=(vkar/log(zmin/.001))**2 ! needs recalc, as aft is ice only
      do iq=1,ifull
       if(.not.land(iq))then 
!        tsurf(iq)=fracice(iq)*tgg(iq,3)+(1.-fracice(iq))*tgg(iq,2) 
         tsurf(iq)=fracice(iq)*tggsn(iq,1)+(1.-fracice(iq))*tpan(iq) ! Dec05 ! MJT seaice
         aft2(iq)=fracice(iq)*aft2ice+(1.-fracice(iq))*chnscr
         aft10(iq)=fracice(iq)*aft10ice+(1.-fracice(iq))*chn10
         aft(iq)=fracice(iq)*aft38ice+(1.-fracice(iq))*chnsea
       endif   ! (.not.land(iq)) 
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
       !tstarx(iq)=aft(iq)*vmod(iq)*fh38(iq)*deltat/ustar(iq) ! MJT cable
       !qstarx(iq)=aft(iq)*vmod(iq)*fh38(iq)*deltaq/ustar(iq) ! MJT cable
       tstarx(iq)=aft(iq)*fh38(iq)*deltat/sqrt(af(iq)*fm38(iq)) ! MJT cable
       qstarx(iq)=aft(iq)*fh38(iq)*deltaq/sqrt(af(iq)*fm38(iq)) ! MJT cable
       fact=sqrt(af2(iq)*fm2(iq))/(aft2(iq)*fh2(iq))
       tscrn(iq)=tsurf(iq)-fact*tstarx(iq)
       
c      N.B. over unstable sea points, may sometimes get supersat qgscrn
c      also over stable snow points, as tscrn & qgscrn derived independently in
c      location with large vertical gradients
       qgscrn(iq)=qsurf(iq)-fact*qstarx(iq)
     
c      screen wind speeds
       vfact=sqrt(u(iq,1)**2+v(iq,1)**2)/vmod(iq)
       !uscrn(iq) = vfact*ustar(iq)/(afroot2(iq)*sqrt(fm2(iq))) ! MJT cable
       !u10(iq) = vfact*ustar(iq)/(afroot10(iq)*sqrt(fm10(iq))) ! MJT cable
       uscrn(iq) = vfact*vmod(iq)*sqrt(af(iq)*fm38(iq))
     &  /(afroot2(iq)*sqrt(fm2(iq))) ! MJT cable
       u10(iq) = vfact*vmod(iq)*sqrt(af(iq)*fm38(iq))
     &  /(afroot10(iq)*sqrt(fm10(iq))) ! MJT cable
       if(ntest==-1)then
         vmag=sqrt(u(iq,1)**2+v(iq,1)**2)
         if(u10(iq)>vmag)then
           print *,'strange iq,vmod,vmag,u10/vmod,u10/vmag ',
     &              iq,vmod(iq),vmag,u10(iq)/vmod(iq),u10(iq)/vmag
           print *,'ri,rich10,ustar,fm10 ',
     &              ri(iq),rich10(iq),ustar(iq),fm10(iq)
         endif
       endif
!      use tscrn for calc rhscrn  23/12/05	
       es=establ(tscrn(iq))
!      rhscrn(iq)=100.*qgscrn(iq)*(ps(iq)-es)/(.622*es)
       qsttg(iq)= .622*es/(ps(iq)-es)  ! recalc     
       rhscrn(iq) = 100.*qgscrn(iq)/qsttg(iq)
      enddo

      if(diag.and.mydiag)then
        iq=idjd
        fact=sqrt(af2(iq)*fm2(iq))/(aft2(iq)*fh2(iq))
        print *,'in scrnout qsurf,qstarx,qgscrn,qg1 ',
     &                      qsurf(iq),qstarx(iq),qgscrn(iq),qg(iq,1)
        print *,'tsurf,tstarx,tscrn,t1 ',
     &           tsurf(iq),tstarx(iq),tscrn(iq),t(iq,1)
        print *,'tpan,tgg2,ri,fact ',tpan(iq),tgg(iq,2),ri(iq),fact
        print *,'in scrnout; epot,eg,fg: ',epot(idjd),eg(idjd),fg(idjd)
        print *,'aft,vmod,fh38,ustar ',
     &           aft(iq),vmod(iq),fh38(iq),ustar(iq)
        print *,'af2,fh2,fm2 ',af2(iq),fh2(iq),fm2(iq) 
        print *,'new qsttg: ',qsttg(iq)
        print *,'uscrn,u10,rhscrn ',uscrn(iq),u10(iq),rhscrn(iq)
      endif

      if(ntest>=1.and.myid==0)then
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
            write (6,"('checking T iq,land,sicedep,tsurf,tscrn,theta,'
     .       'tstarx ',i5,l2,4f7.2,f7.3)") iq,land(iq),sicedep(iq), 
     .        tsurf(iq),tscrn(iq),theta(iq),tstarx(iq)
            write (6,"('tss,tgg1,tpan,tgg3,t1,fracice ',6f7.2)")
     .         tss(iq),tgg(iq,1),tpan,tgg(iq,3),t(iq,1),fracice(iq)
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
            if(sicedep(iq)>0.)print *,'tpan,tgg3,fracice ',
     .                           tpan(iq),tgg(iq,3),fracice(iq)       
          endif
          if(qgscrn(iq)<min(qg(iq,1),qsurf(iq)).or.
     .      qgscrn(iq)>max(qg(iq,1),qsurf(iq))    )then
            numq=numq+1
            print *,'strange qg iq,land,sicedep,qstarx,t,ri ',
     .       iq,land(iq),sicedep(iq),qstarx(iq),t(iq,1),ri(iq),numq
            print *,'qsurf,qgscrn,qg1 ',qsurf(iq),qgscrn(iq),qg(iq,1)
            print *,'tsurf,tscrn,theta ',tsurf(iq),tscrn(iq),theta(iq)
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
      if(ntest==-1.and.myid==0)then
       do iq=1,ifull
        if(.not.land(iq).and.sicedep(iq)==0.)then 
c         write(23,'(f6.2,f9.6,2i6,f7.2)') u10(iq),zo(iq),ktau,iq,ri(iq)
          write(23,'(f6.2,2f9.6,f6.2,2i6,3f7.4)') u10(iq),
     &     cduv(iq)/vmod(iq),zo(iq),vmod(iq),ktau,iq,ri(iq),fh10(iq),
     &     ustar(iq)
c          to plot cd10, use plot 'fort.23' u
c                   1:($2*(log(38/$3)**2/log(10/$3)**2))  
        endif
       enddo        
      endif
      return
      end
      
      ! MJT cable ---------------------------------------------------
      ! Use TAPM approach to screen diagnostics for ocean and lake points
      subroutine scrnocn(ifull,qgscrn,tscrn,uscrn,u10,rhscrn,zo,zoh,
     &                   tsu,temp,qsttg,qg,umag,ps,land,zmin,sig)
     
      implicit none

      integer, intent(in) :: ifull
      integer pfull
      real, dimension(ifull), intent(inout) :: qgscrn,tscrn,uscrn,u10
      real, dimension(ifull), intent(inout) :: rhscrn
      real, dimension(ifull), intent(in) :: zo,zoh,tsu,temp,qsttg
      real, dimension(ifull), intent(in) :: qg,umag,ps
      real, intent(in) :: zmin,sig
      real, dimension(ifull) :: qgscrn_pack,tscrn_pack
      real, dimension(ifull) :: uscrn_pack,u10_pack
      real, dimension(ifull) :: rhscrn_pack,zo_pack,zoh_pack
      real, dimension(ifull) :: stemp_pack,temp_pack
      real, dimension(ifull) :: smixr_pack,mixr_pack
      real, dimension(ifull) :: umag_pack,ps_pack
      logical, dimension(ifull), intent(in) :: land
      
      pfull=count(.not.land)

      if (pfull.eq.0) return

      zo_pack(1:pfull)=pack(zo,.not.land)
      zoh_pack(1:pfull)=pack(zoh,.not.land)
      stemp_pack(1:pfull)=pack(tsu,.not.land)
      temp_pack(1:pfull)=pack(temp,.not.land)
      smixr_pack(1:pfull)=pack(qsttg,.not.land)
      mixr_pack(1:pfull)=pack(qg,.not.land)
      umag_pack(1:pfull)=pack(umag,.not.land)
      ps_pack(1:pfull)=pack(ps,.not.land)

      call scrncalc(pfull,qgscrn_pack(1:pfull),tscrn_pack(1:pfull),
     &              uscrn_pack(1:pfull),u10_pack(1:pfull),
     &              rhscrn_pack(1:pfull),zo_pack(1:pfull),
     &              zoh_pack(1:pfull),stemp_pack(1:pfull),
     &              temp_pack(1:pfull),smixr_pack(1:pfull),
     &              mixr_pack(1:pfull),umag_pack(1:pfull),
     &              ps_pack(1:pfull),zmin,sig)

      qgscrn=unpack(qgscrn_pack(1:pfull),.not.land,qgscrn)
      tscrn=unpack(tscrn_pack(1:pfull),.not.land,tscrn)
      uscrn=unpack(uscrn_pack(1:pfull),.not.land,uscrn)
      u10=unpack(u10_pack(1:pfull),.not.land,u10)
      rhscrn=unpack(rhscrn_pack(1:pfull),.not.land,rhscrn)

      return
      end subroutine scrnocn
      
      subroutine scrncalc(pfull,qscrn,tscrn,uscrn,u10,rhscrn,zo,zoh,
     &                    stemp,temp,smixr,mixr,umag,ps,zmin,sig)
      
      implicit none

      include 'const_phys.h'

      integer, intent(in) :: pfull
      integer ic
      real, dimension(pfull), intent(out) :: qscrn,tscrn,uscrn,u10
      real, dimension(pfull), intent(out) :: rhscrn
      real, dimension(pfull), intent(in) :: zo,zoh,stemp,temp,umag
      real, dimension(pfull), intent(in) :: smixr,mixr,ps
      real, dimension(pfull) :: lzom,lzoh,thetav,sthetav
      real, dimension(pfull) :: thetavstar,z_on_l,z0_on_l,zt_on_l
      real, dimension(pfull) :: pm0,ph0,pm1,ph1,integralm,integralh
      real, dimension(pfull) :: ustar,qstar,z10_on_l
      real, dimension(pfull) :: neutral,neutral10,pm10
      real, dimension(pfull) :: integralm10
      real, dimension(pfull) :: esat,qsat,tstar
      real, intent(in) :: zmin,sig
      real scrp
      integer, parameter ::  nc     = 5
      real, parameter    ::  vkar   = 0.4
      real, parameter    ::  a_1    = 1.
      real, parameter    ::  b_1    = 2./3.
      real, parameter    ::  c_1    = 5.
      real, parameter    ::  d_1    = 0.35
      real, parameter    ::  z0     = 1.5
      real, parameter    ::  z10    = 10.

      scrp=(sig)**(rdry/cp)
      thetav=temp*(1.+0.61*mixr)/scrp
      sthetav=stemp*(1.+0.61*smixr)

      ! Roughness length for heat
      lzom=log(zmin/zo)
      lzoh=log(zmin/zoh)

      ! Dyer and Hicks approach 
      thetavstar=vkar*(thetav-sthetav)/lzoh
      ustar     =vkar*umag/lzom
      do ic=1,nc
        z_on_l=vkar*zmin*grav*thetavstar/(thetav*ustar**2)
        z_on_l=min(z_on_l,10.)
        z0_on_l  = z_on_l*zo/zmin
        zt_on_l  = z_on_l*zoh/zmin
        where (z_on_l.lt.0.)
          pm0     = (1.-16.*z0_on_l)**(-0.25)
          ph0     = (1.-16.*zt_on_l)**(-0.5)
          pm1     = (1.-16.*z_on_l)**(-0.25)
          ph1     = (1.-16.*z_on_l)**(-0.5)
          integralm = lzom-2.*log((1.+1./pm1)/(1.+1./pm0))
     &                -log((1.+1./pm1**2)/(1.+1./pm0**2))
     &                +2.*(atan(1./pm1)-atan(1./pm0))
          integralh = lzoh-2.*log((1.+1./ph1)/(1.+1./ph0))
        elsewhere
          !--------------Beljaars and Holtslag (1991) momentum & heat            
          pm0 = -(a_1*z0_on_l+b_1*(z0_on_l-(c_1/d_1))*exp(-d_1*z0_on_l)
     &          +b_1*c_1/d_1)
          pm1 = -(a_1*z_on_l+b_1*(z_on_l-(c_1/d_1))*exp(-d_1*z_on_l)
     &          +b_1*c_1/d_1)
          ph0 = -((1.+(2./3.)*a_1*zt_on_l)**1.5
     &          +b_1*(zt_on_l-(c_1/d_1))*exp(-d_1*zt_on_l)
     &          +b_1*c_1/d_1-1.)
          ph1 = -((1.+(2./3.)*a_1*z_on_l)**1.5
     &          +b_1*(z_on_l-(c_1/d_1))*exp(-d_1*z_on_l)
     &          +b_1*c_1/d_1-1.)
          integralm = lzom-(pm1-pm0)    
          integralh = lzoh-(ph1-ph0)         
        endwhere
        thetavstar=vkar*(thetav-sthetav)/integralh
        ustar     =vkar*umag/integralm
      end do
      tstar=vkar*(temp-stemp)/integralh
      qstar=vkar*(mixr-smixr)/integralh
      
      ! estimate screen diagnostics
      z0_on_l   = z0*z_on_l/zmin
      z10_on_l  = z10*z_on_l/zmin
      z0_on_l   = min(z0_on_l,10.)
      z10_on_l  = min(z10_on_l,10.)
      neutral   = log(zmin/z0)
      neutral10 = log(zmin/z10)
      where (z_on_l.lt.0.)
        ph0     = (1.-16.*z0_on_l)**(-0.50)
        ph1     = (1.-16.*z_on_l)**(-0.50)
        pm0     = (1.-16.*z0_on_l)**(-0.25)
        pm10    = (1.-16.*z10_on_l)**(-0.25)
        pm1     = (1.-16.*z_on_l)**(-0.25)
        integralh = neutral
     &              -2.*log((1.+1./ph1)/(1.+1./ph0))
        integralm = neutral-2.*log((1.+1./pm1)/(1.+1./pm0))
     &                -log((1.+1./pm1**2)/(1.+1./pm0**2))
     &                +2.*(atan(1./pm1)-atan(1./pm0))
        integralm10 = neutral10-2.*log((1.+1./pm1)/(1.+1./pm10))
     &                -log((1.+1./pm1**2)/(1.+1./pm10**2))
     &                +2.*(atan(1./pm1)-atan(1./pm10))     
      elsewhere
c-------Beljaars and Holtslag (1991) heat function
        ph0  = -((1.+(2./3.)*a_1*z0_on_l)**1.5
     &         +b_1*(z0_on_l-(c_1/d_1))
     &         *exp(-d_1*z0_on_l)+b_1*c_1/d_1-1.)
        ph1  = -((1.+(2./3.)*a_1*z_on_l)**1.5
     &         +b_1*(z_on_l-(c_1/d_1))
     &         *exp(-d_1*z_on_l)+b_1*c_1/d_1-1.)
        pm0 = -(a_1*z0_on_l+b_1*(z0_on_l-(c_1/d_1))*exp(-d_1*z0_on_l)
     &        +b_1*c_1/d_1)
        pm10 = -(a_1*z10_on_l+b_1*(z10_on_l-(c_1/d_1))
     &         *exp(-d_1*z10_on_l)+b_1*c_1/d_1)
        pm1  = -(a_1*z_on_l+b_1*(z_on_l-(c_1/d_1))*exp(-d_1*z_on_l)
     &         +b_1*c_1/d_1)
        integralh = neutral-(ph1-ph0)
        integralm = neutral-(pm1-pm0)
        integralm10 = neutral10-(pm1-pm10)
      endwhere
      tscrn = temp-tstar*integralh/vkar
      qscrn = mixr-qstar*integralh/vkar
      qscrn = max(qscrn,1.E-4)
      where (tscrn.ge.273.15)
        esat = 610.*exp(hl/rvap*(1./273.15-1./tscrn))
      elsewhere
        esat = 610.*exp(hls/rvap*(1./273.15-1./tscrn))
      endwhere
      qsat = 0.622*esat/(ps-0.378*esat)
      rhscrn = 100.*min(qscrn/qsat,1.)
      
      uscrn=max(umag-ustar*integralm/vkar,0.)
      u10  =max(umag-ustar*integralm10/vkar,0.)
      
      return
      end subroutine scrncalc
      !--------------------------------------------------------------
