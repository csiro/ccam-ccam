      subroutine sflux(nalpha)              ! for globpe code
      use diag_m
      use cc_mpi
      parameter (nblend=0)  ! 0 for original non-blended, 1 for blended af
      parameter (ndvmod=0)  ! 0 default, 1+ for dvmod tests
      parameter (ntss_sh=0) ! 0 for original, 3 for **3, 4 for **4
      parameter (nplens=0)  ! 0 to turn off plens, 10 (e.g.) is on
!     parameter (lake=0)    ! 0 usual, 1 for specified lake points
!                             - replaced by nspecial in parm.h
      parameter (ntest=0)   ! ntest= 0 for diags off; ntest= 1 for diags on
      parameter (ntaft=2)   ! 0 for original,
!                   1 & 2 tafthf constrained by prior values with 2 faster
!                   3 uses measure of prior tgf in calc. fh
!     parameter (newztsea=0)   ! 0 for original, 1 for correct zt over sea
!     From 11/8/98 runoff() is accumulated & zeroed with precip
!     Now using tgg(,3) for the tice calculations
c     with leads option via fracice (using tgg1 and tgg3)
c     now in parm.h with zero as default and used in rdnsib/rdnscam
c     cp specific heat at constant pressure joule/kgm/deg
      parameter (bprm=5.,cms=5.,chs=2.6,vkar=.4)
      parameter (d3=2.5)
      parameter (cgsoil=1000.,gksoil=.300e-6,rhog=1600.,d1land=.03)
      parameter (fmroot=.57735)     ! was .4 till 7 Feb 1996
!     parameter (chn10=.00136733)   ! sea only - in parm.h from 2004
      include 'newmpar.h'
      include 'arrays.h'
      include 'const_phys.h'
      include 'extraout.h' ! ustar
      include 'gdrag.h'
      include 'liqwpar.h'  ! qfg,qlg
      include 'map.h'      ! land
      include 'morepbl.h'  ! condx,fg,eg
      include 'nsibd.h'    ! rsmin,ivegt,sigmf,tgf,ssdn,res,rmc,tsigmf
      include 'parm.h'
      include 'pbl.h'
      include 'permsurf.h'
      include 'prec.h'
      include 'savuvt.h'
      include 'scamdim.h'  ! dimension of patches
      include 'screen.h'   ! tscrn,qgscrn,uscrn,scrrel,u10
      include 'sigs.h'
      include 'soil.h'     ! ... zmin zolnd zolog sice fracice alb
      include 'soilv.h'    ! ... ssat
      include 'soilsnow.h' ! new soil arrays for scam - tgg too
      include 'tracers.h'  ! ngas, nllp, ntrac
      include 'trcom2.h'   ! nstn,slat,slon,istn,jstn
      include 'vvel.h'
      common/tafts/taftfh(ifull),taftfhg(ifull)
      common/work2/dirad(ifull),dfgdt(ifull),degdt(ifull)
     . ,wetfac(ifull),degdw(ifull),cie(ifull)
     . ,factch(ifull),qsttg(ifull),rho(ifull),zo(ifull)
     . ,aft(ifull),fh(ifull),ri(ifull),theta(ifull)
     . ,gamm(ifull),rg(ifull),vmod(ifull),taftfhg_temp(ifull) ! rg in radriv90
!     following common block makes available other arrays for diag. output 
      common/work3/egg(ifull),evapxf(ifull),Ewww(ifull),fgf(ifull),
     . fgg(ifull),ggflux(ifull),rdg(ifull),rgg(ifull),residf(ifull),
     . ga(ifull),condxpr(ifull),fev(ifull),fes(ifull),
     . ism(ifull),fwtop(ifull),epot(ifull),   ! watch soilsnow.f after epot
     . extin(ifull),af(ifull),spare1(ifull),xx(ifull),
     . dum3(5*ijk-20*ifull)
      dimension ipermp(ifull)    ! temporary permutation array
      equivalence (ipermp,dirad)
      real plens(ifull)
      save plens
      data plens/ifull*0./
      include 'establ.h'

c     stability dependent drag coefficients using Louis (1979,blm) f'
c     n.b. cduv, cdtq are returned as drag coeffs mult by vmod
c          (cduv=cduv*wmag; cdtq=cdtq*wmag; vmod=wmag)

c     t, u, v, qg are current values
c     tss is surface temperature
c     dw is soil wetness availability (1. for ocean) - not needed here
c     fg is sensible heat flux (was h0)
c     eg is latent heat flux (was wv)
c     dfgdt is dfgdt (was csen in surfupa/b)
c     degdt is degdt (was ceva in surfupa/b)

      zobgin = .05   ! jlm: NB seems to be .01 in csiro9 Fri  12-06-1996
      alzzin=log(zmin/zobgin)   ! pre-calculated for all except snow points
      ztv=exp(vkar/sqrt(chn10)) /10.  ! proper inverse of ztsea
      z1onzt=300.*rdry*(1.-sig(1))*ztv /grav
      chnsea=(vkar/log(z1onzt))**2    ! should give .00085 for csiro9

      if(nspecial.eq.1)then
         print*, "SFLUX, nspecial==1 doesn't work in MPI version"
         stop
        do j=68,71  ! Eyre
         do i=27,28
          iq=i+(j-1)*il
	   sigmf(iq)=0.  ! bare "soil"
	   do kk=1,ms
           wb(iq,kk)=ssat(isoilm(iq))
	   enddo
	   if(ktau.eq.1)then
	     print *,'lake iq,isoilm,ivegt,sigmf,zolnd ',
     .                    iq,isoilm(iq),ivegt(iq),sigmf(iq),zolnd(iq)
	   endif ! (ktau.eq.1)
         enddo
        enddo
        do j=64,67  ! Torrens
         do i=28,28
          iq=i+(j-1)*il
	   sigmf(iq)=0.  ! bare "soil"
	   do kk=1,ms
           wb(iq,kk)=ssat(isoilm(iq))
	   enddo
	   if(ktau.eq.1)then
	     print *,'lake iq,isoilm,ivegt,sigmf,zolnd ',
     .                    iq,isoilm(iq),ivegt(iq),sigmf(iq),zolnd(iq)
	   endif ! (ktau.eq.1)
         enddo
        enddo
      endif   !  (nspecial.eq.1)

      if(ktau.eq.1.and.nrungcm.eq.3)then  ! for runs for PIRCS
        if (mydiag) print *,'entering sflux w2: ',wb(idjd,ms)
c       due to uncertainty about NCEP value of field capacity, it may be
c       more appropriate to define it as 75% of ssat (for PIRCS93)
        sfc(2)=.75*ssat(2)
        sfc(3)=.75*ssat(3)
!cdir nodep
        do ip=1,ipland  ! all land points in this nsib=1 loop
         iq=iperm(ip)
         isoil = isoilm(iq)
c        code where eta provided (e.g. for PIRCS)
         wb(iq,1)=swilt(isoil)+wb(iq,1)*(sfc(isoil)-swilt(isoil))
         wb(iq,ms)=swilt(isoil)+wb(iq,ms)*(sfc(isoil)-swilt(isoil))
        enddo        !  ip=1,ipland
c        if(nrungcm.eq.4)then   ! for test runs only
c!cdir nodep
c          do ip=1,ipland
c           iq=iperm(ip)
c           isoil = isoilm(iq)
c           wb(iq,1)=sfc(isoil)
c           wb(iq,ms)=sfc(isoil)
c        enddo        !  ip=1,ipland
c        endif  ! (nrungcm.eq.4)
        if (mydiag) print *,'then w2: ',wb(idjd,ms)
      endif          !  (ktau.eq.1.and.nrungcm.eq.3)

!     at start of each day allow for sice mask changes
      if(mod(ktau,nperday).eq.1)then
        if(nbd.ne.0)then
!         update sice from tss (see nestin, which provides tss, tgg1)
!         check whether present ice points should change to/from sice points
          do iq=1,ifull
           if(.not.land(iq))then  
             if(tss(iq).gt.271.2)then  
               fracice(iq)=0.
             else
               if(fracice(iq).eq.0.)then
!                N.B. if already a sice point, keep present tice
!                create values for tice, and set averaged tss
                 tgg(iq,3)=min(271.2,tss(iq),t(iq,1)+.04*6.5) ! for 40 m level 1
                 fracice(iq)=.5    
               endif  ! (fracice(iq).eq.0.)
               tss(iq)=tgg(iq,3)*fracice(iq)+tgg(iq,1)*(1.-fracice(iq))
             endif   !  (tss(iq).gt.271.2) ... else ...
           endif   ! (.not.land(iq))
          enddo	     ! iq loop
        endif      ! (nbd.ne.0)
	      
        if(nbd.ne.0.or.namip.ne.0)then ! update surface permutation array
!         N.B. just done once per day
!cdir nodep
          do ip=ipland+1,ipsea
           ipermp(ip)=iperm(ip)
          enddo   !  ip loop
          indexi=ipland
          indexs=ipsea+1
!cdir nodep
          do ip=ipland+1,ipsea
           iq=ipermp(ip)
           if(fracice(iq).gt.0.)then
             indexi=indexi+1     ! sice point
             iperm(indexi)=iq    ! sice point
             sice(iq)=.true.
             sicedep(iq)=2.   
           else
             indexs=indexs-1     ! sea point
             iperm(indexs)=iq    ! sea point
             sice(iq)=.false.
             sicedep(iq)=0.
           endif  ! (fracice(iq).gt.0.)
          enddo   !  ip loop
          ipsice=indexi
        endif        ! (nbd.ne.0.or.namip.ne.0)
        if ( myid == 0 ) then
           print *,'ktau,ipland,ipsice,ipsea update in sflux: ',
     .              ktau,ipland,ipsice,ipsea
           print *,'nblend,newztsea,ntaft,ztv,chnsea: ',
     .              nblend,newztsea,ntaft,ztv,chnsea   
        end if
      endif          ! (mod(ktau,nperday).eq.1)

      if (diag.or.ntest.eq.1) then
        if (mydiag) then
          print *,'entering sflux ktau,nsib,ivegt,isoilm,land '
     .         ,ktau,nsib,ivegt(idjd),isoilm(idjd),land(idjd)
          print *,'idjd,id,jd,slwa,sgsave ',
     .           idjd,id,jd,slwa(idjd),sgsave(idjd)
          print *,'snowd,sice,condx ',
     .           snowd(idjd),sice(idjd),condx(idjd)
          print *,'t1,tss ',t(idjd,1),tss(idjd)
          print *,'wb ',(wb(idjd,k),k=1,ms)
          print *,'tgg ',(tgg(idjd,k),k=1,ms)
        end if
        call maxmin(t,' t',ktau,1.,kl)
      endif

c     using av_vmod (1. for no time averaging)
!      *****  check next comment
!       sflux called at beginning of time loop, hence savu, savv

      srcp =sig(1)**(rdry/cp)
      ga(:)=0.              !  for ocean points in ga_ave diagnostic
      theta(:)=t(1:ifull,1)/srcp
      rho(:)=ps(1:ifull)/(rdry*tss(:))
      do iq=1,ifull
       uav=av_vmod*u(iq,1)+(1.-av_vmod)*savu(iq,1)   
       vav=av_vmod*v(iq,1)+(1.-av_vmod)*savv(iq,1)  
       ustar(iq)=sqrt(uav**2+vav**2)  ! i.e. vmod for tss_sh
      enddo
      if(ndvmod.eq.0)then
        vmod(:)=max( ustar(:) , vmodmin)
      else
        vmod(:)=ndvmod  ! just for tests
      endif    ! (ndvmod.eq.0)

      if(namip.lt.2)then  ! generalize later
        ipsea0=ipsice+1   ! used without leads
      else
        ipsea0=ipland+1   ! used with leads  - just namip=2 at present
      endif
      if(ntest.eq.2.and.mydiag)print *,'before sea loop'
!      from June '03 use basic sea temp from tgg1 (so leads is sensible)      
!cdir nodep
      do ip=ipsea0,ipsea                                                ! sea
!      all sea points in this big loop; also open water of leads        ! sea
       iq=iperm(ip)                                                     ! sea
       wetfac(iq)=1.                                                    ! sea
!      tgg2 holds effective skin sst for this loop 
       if(ntss_sh.eq.0)then
         dtsol=.01*sgsave(iq)/(1.+.25*ustar(iq)**2)    ! solar heating  ! sea
         tgg(iq,2)=tgg(iq,1)+tss_sh*min(dtsol,8.)      ! of ssts        ! sea
       elseif(ntss_sh.eq.1)then                                         ! sea
         dtsol=tss_sh*.01*sgsave(iq)/                                   ! sea
     .                (1.+.25*ustar(iq)**2)            ! solar heating  ! sea
         tgg(iq,2)=tgg(iq,1)+min(dtsol,8.)             ! of ssts        ! sea
       elseif(ntss_sh.eq.3)then                                         ! sea
         dtsol=tss_sh*.01*sgsave(iq)/                                   ! sea
     .                (1.+.035*ustar(iq)**3)           ! solar heating  ! sea
         tgg(iq,2)=tgg(iq,1)+min(dtsol,8.)             ! of ssts        ! sea
       elseif(ntss_sh.eq.4)then                                         ! sea
         dtsol=tss_sh*.01*sgsave(iq)/                                   ! sea
     .                (1.+ustar(iq)**4/81.)            ! solar heating  ! sea
         tgg(iq,2)=tgg(iq,1)+min(dtsol,8.)             ! of ssts        ! sea
       endif   ! (ntss_sh.eq.0) .. else ..
       if(nplens.ne.0)then
!        calculate running total of daily precip in mm/day   jlm
         plens(iq)=(1.-dt/86400.)*plens(iq)+condx(iq)  ! in mm/day
!        scale so that nplens m/s wind for 1/2 hr reduces effect by 1/1.2
!        plens(iq)=plens(iq)/(1.+ustar(iq)*dt*.2/(nplens*1800.))
         plens(iq)=plens(iq)/(1.+ustar(iq)*dt*.2/
     .                    max(nplens*1800.,1.))      ! avoids Cray compiler bug
!        produce a cooling of 4 K for an effective plens of 10 mm/day
         tgg(iq,2)=tgg(iq,2)-min(.4*plens(iq) , 6.)
       endif   !  (nplens.ne.0)
c ***  drag coefficients  for momentum           cduv                   ! sea
c ***                     for heat and moisture  cdtq                   ! sea
       es = establ(tgg(iq,2))                                           ! sea
       constz=ps(iq)-es                                                 ! sea
       qsttg(iq)= .98*.622*es/constz   ! with Zeng 1998 for sea water   ! sea
       drst=qsttg(iq)*ps(iq)*hlars/(constz*tgg(iq,2)**2)                ! sea
       xx(iq)=grav*zmin*(1.-tgg(iq,2)*srcp/t(iq,1))                     ! sea
       ri(iq)=xx(iq)/vmod(iq)**2                                        ! sea
!      if(ngas.gt.0)stop 'call co2sflux'                                ! sea
c      this is in-line ocenzo using latest coefficient, i.e. .018       ! sea
       consea=vmod(iq)*.018/grav                                        ! sea
       zo(iq)=.01                                                       ! sea
       if(xx(iq).gt.0.)then             ! stable sea points             ! sea
         fm=vmod(iq) /(1.+bprm*ri(iq))**2   ! N.B. this is vmod*fm      ! sea
         con=consea*fm                                                  ! sea
         do it=1,3                                                      ! sea
          afroot=vkar/log(zmin/zo(iq))                                  ! sea
          af(iq)=afroot**2                                              ! sea
          daf=2.*af(iq)*afroot/(vkar*zo(iq))                            ! sea
          zo(iq)=max(1.5e-5,zo(iq)-(zo(iq)-con*af(iq))/
     .                                                 (1.-con*daf))    ! sea
         enddo    ! it=1,3                                              ! sea
         afroot=vkar/log(zmin/zo(iq))                                  ! sea
         af(iq)=afroot**2                                              ! sea
       else                        ! unstable sea points                ! sea
         do it=1,3                                                      ! sea
          afroot=vkar/log(zmin/zo(iq))                                  ! sea
          af(iq)=afroot**2                                              ! sea
          daf=2.*af(iq)*afroot/(vkar*zo(iq))                            ! sea
          con1=cms*2.*bprm*sqrt(-xx(iq)*zmin/zo(iq))                    ! sea
          den=vmod(iq)+af(iq)*con1                                      ! sea
          dden=con1*(daf-.5*af(iq)/zo(iq))                              ! sea
          fm=vmod(iq)-(2.*bprm *xx(iq))/den                             ! sea
          dfm=2.*bprm*xx(iq)*dden/den**2                                ! sea
          zo(iq)=max(1.5e-5,zo(iq)-(zo(iq)-consea*af(iq)*fm)/           ! sea
     .                     (1.-consea*(daf*fm+af(iq)*dfm)))             ! sea
         enddo    ! it=1,3                                              ! sea
       endif    ! (xx.gt.0.) .. else..                                  ! sea
       aft(iq)=chnsea                                                   ! sea
      enddo    ! ip=ipsea0,ipsea

      if(newztsea.eq.0)then ! 0 for original, 1 for different zt over sea
c       actually only need zt (& thus factch) for unstable points       ! sea
c       factch is sqrt(zo/zt) for use in unstable fh                    ! sea
c       enhanced formula used in Feb '92 Air-Sea conference follows:    ! sea
c       factch=sqrt(zo*exp(vkar*vkar/(chnsea*log(zmin/zo)))/zmin)       ! sea
        do iq=1,ifull                                                   ! sea
         factch(iq)=1. ! factch is sqrt(zo/zt) for use in unstable fh
        enddo
      else
!cdir nodep
        do ip=ipsea0,ipsea                                              ! sea
         iq=iperm(ip)                                                   ! sea
         factch(iq)=sqrt(zo(iq)*ztv) ! for use in unstable fh
        enddo
      endif  ! (newztsea.eq.0)

!cdir nodep
      do ip=ipsea0,ipsea                                                ! sea
       iq=iperm(ip)                                                     ! sea
c      Having settled on zo (and thus af) now do actual fh and fm calcs ! sea
       if(xx(iq).gt.0.)then                                             ! sea
         fm=vmod(iq)*max(fmroot*fmroot,1./(1.+bprm*ri(iq))**2)          ! sea
         fh(iq)=fm                                                      ! sea
       else        ! xx is -ve                                          ! sea
         root=sqrt(-xx(iq)*zmin/zo(iq))                                 ! sea
c        First do momentum                                              ! sea
         denma=vmod(iq)+cms*2.*bprm*af(iq)*root                         ! sea
         fm=vmod(iq)-(2.*bprm *xx(iq))/denma                            ! sea
c        n.b. fm denotes ustar**2/(vmod(iq)*af)                         ! sea
c        Now heat ; allow for smaller zo via aft and factch             ! sea
         denha=vmod(iq)+chs*2.*bprm*factch(iq)*aft(iq)*root             ! sea
         fh(iq)=vmod(iq)-(2.*bprm *xx(iq))/denha                        ! sea
       endif                                                            ! sea
                                                                        ! sea
       conh=rho(iq)*aft(iq)*cp                                          ! sea
       conw=rho(iq)*aft(iq)*hl                                          ! sea
	qgtot=qg(iq,1)+qfg(iq,1)+qlg(iq,1)
       fg(iq)=conh*fh(iq)*(tgg(iq,2)-theta(iq))                         ! sea
       eg(iq)=conw*fh(iq)*(qsttg(iq)-qgtot)                             ! sea
       epot(iq) = eg(iq)                                                ! sea
c      cduv is now drag coeff *vmod                                     ! sea
       cduv(iq) =af(iq)*fm                                              ! sea
       ustar(iq) = sqrt(vmod(iq)*cduv(iq))                              ! sea
c      Surface stresses taux, tauy: diagnostic only - unstaggered now   ! sea
       taux(iq)=rho(iq)*cduv(iq)*u(iq,1)                                ! sea
       tauy(iq)=rho(iq)*cduv(iq)*v(iq,1)                                ! sea
       ! note that iq==idjd  can only be true on the correct process
       if(ntest.eq.1.and.iq.eq.idjd)then                                ! sea
         print *,'in sea loop for iq,idjd,ip,ipsea0: ',                 ! sea
     .                            iq,idjd,ip,ipsea0                     ! sea
         print *,'zmin,zo,factch ',zmin,zo(iq),factch(iq)               ! sea         
         print *,'xx,ri,ustar,es ',xx(iq),ri(iq),ustar(iq),es           ! sea
         print *,'af,aft,tgg2 ',af(iq),aft(iq),tgg(iq,2)                ! sea
         print *,'tgg2,tss,theta ',tgg(iq,2),tss(iq),theta(iq)          ! sea
         print *,'chnsea,rho,t1 ',chnsea,rho(iq),t(iq,1)                ! sea
         print *,'fm,fh,conh ',fm,fh(iq),conh                           ! sea
         print *,'vmod,cduv,fg ',vmod(iq),cduv(iq),fg(iq)               ! sea
       endif                                                            ! sea
      enddo    !  ip=ipsea0,ipsea                                       ! sea
      if(ntest.eq.2.and.mydiag)print *,'after sea loop'

      zminlog=log(zmin)
!cdir nodep
      do ip=ipland+1,ipsice  ! all sea-ice points in this loop          ! sice
!      non-leads for sea ice points                                     ! sice
!      N.B. tgg( ,3) holds tice                                         ! sice
       iq=iperm(ip)                                                     ! sice
       qgtot=qg(iq,1)+qfg(iq,1)+qlg(iq,1)
       es = establ(tgg(iq,3))                                           ! sice
       constz=ps(iq)-es                                                 ! sice
       qsttg(iq)= .622*es/constz                                        ! sice
       drst=qsttg(iq)*ps(iq)*hlars/(tgg(iq,3)*tgg(iq,3)*constz)         ! sice
       xx(iq)=grav*zmin*(1.-tgg(iq,3)*srcp/t(iq,1))                     ! sice
       ri_ice=xx(iq)/vmod(iq)**2                                        ! sice
       factch(iq)=sqrt(7.4)  ! same as land from 27/4/99                ! sice
!      factch(iq)=1.   ! factch is sqrt(zo/zt) for use in unstable fh   ! sice
       zoice=.001                                                       ! sice
       zologice=zminlog-log(zoice)   !   i.e. log(zmin/zo(iq))          ! sice
       af(iq)=(vkar/zologice)**2                                        ! sice
       aft(iq)=vkar**2/(zologice*(2.+zologice) )  ! from 27/4/99        ! sice
!      aft(iq)=af                                 ! up till 27/4/99     ! sice
       wetfac(iq)=1+.008*(tgg(iq,3)-273.16)  ! .008*tgg(iq,3)-1.18528   ! sice
                                                                        ! sice
c      Having settled on zo (and thus af) now do actual fh and fm calcs ! sice
       if(xx(iq).gt.0.)then                                             ! sice
         fm=vmod(iq)*max(fmroot*fmroot,1./(1.+bprm*ri_ice)**2)          ! sice
         fh(iq)=fm                                                      ! sice
       else                                                             ! sice
         root=sqrt(-xx(iq)*zmin/zoice)                                  ! sice
c        First do momentum                                              ! sice
         denma=vmod(iq)+cms*2.*bprm*af(iq)*root                         ! sice
         fm=vmod(iq)-(2.*bprm *xx(iq))/denma                            ! sice
c        n.b. fm denotes ustar**2/(vmod(iq)*af)                         ! sice
c        Now heat ; allow for smaller zo via aft and factch             ! sice
         denha=vmod(iq)+chs*2.*bprm*factch(iq)*aft(iq)*root             ! sice
         fh(iq)=vmod(iq)-(2.*bprm *xx(iq))/denha                        ! sice
       endif                                                            ! sice
                                                                        ! sice
       conh=rho(iq)*aft(iq)*cp                                          ! sice
       conw=rho(iq)*aft(iq)*hl                                          ! sice
       fgice=conh*fh(iq)*(tgg(iq,3)-theta(iq))                          ! sice
       dfgdt(iq)=conh*fh(iq)                                            ! sice
       if(ntest.eq.1.and.iq.eq.idjd)then                                ! sice
         print *,'in sice loop'                                         ! sice
         print *,'zmin,zo,wetfac ',zmin,zoice,wetfac(iq)                ! sice
         print *,'xx,ri_ice,es ',xx(iq),ri_ice,es                       ! sice
         print *,'af,aft,ustar ',af(iq),aft(iq),ustar(iq)               ! sice
         print *,'chnsea,rho ',chnsea,rho(iq)                           ! sice
         print *,'fm,fh,conh ',fm,fh(iq),conh                           ! sice
       endif                                                            ! sice
                                                                        ! sice
       if(nalpha.eq.1)then    ! beta scheme         sice here           ! sice
         epotice=conw*fh(iq)*(qsttg(iq)-qgtot)                          ! sice
         egice  =wetfac(iq)*epotice                                     ! sice
         degdt(iq)=wetfac(iq)*conw*fh(iq)*drst                          ! sice
       else                   ! alpha scheme                            ! sice
c        following trick reduces -ve evap (dew) to 1/10th value         ! sice
         qtgnet=qsttg(iq)*wetfac(iq) -qgtot                             ! sice
         qtgair=qsttg(iq)*wetfac(iq)-max(qtgnet,.1*qtgnet)              ! sice
         eg2=-conw*fh(iq)*qtgair                                        ! sice
         eg1=conw*fh(iq)*qsttg(iq)                                      ! sice
         egice   =eg1*wetfac(iq) +eg2                                   ! sice
         epotice    = conw*fh(iq)*(qsttg(iq)-qgtot)                     ! sice
         deg=wetfac(iq)*conw*fh(iq)*drst                                ! sice
c        following reduces degdt by factor of 10 for dew                ! sice
         degdt(iq)=.55*deg+sign(.45*deg,qtgnet)                         ! sice
       endif                                                            ! sice
                                                                        ! sice
c      section to update sea ice surface temperature;                   ! sice
c      specified sea-ice thickness                                      ! sice
c      over sea ice, set a minimum depth for this experiment of .1      ! sice
       sicedep(iq) = max(sicedep(iq) , 0.1)                             ! sice
c      no snow on the ice assumed for now                               ! sice
       gamm(iq) = 3.471e+05                                             ! sice
       cie(iq) = 2.04/sicedep(iq)                                       ! sice
       rgg(iq)=5.67e-8*tgg(iq,3)**4                                     ! sice
!      gflux here is	flux from ice to water, +ve downwards              ! sice
       gflux(iq)=cie(iq)*(tgg(iq,3)-271.2)                              ! sice
       ga(iq)=-slwa(iq)-rgg(iq)-egice-fgice-gflux(iq)                   ! sice
       dirad(iq)=4.*5.67e-8*tgg(iq,3)**3                                ! sice
       b1=dirad(iq)+degdt(iq)+dfgdt(iq)+cie(iq)                         ! sice
       gbot=(gamm(iq)/dt)+b1                                            ! sice
       deltat=ga(iq)/gbot                                               ! sice
       tgg(iq,3)=tgg(iq,3)+deltat                                       ! sice
       tgg(iq,3)=min(tgg(iq,3),271.2)   ! jlm fix Tue  05-30-2000
       fgice   =fgice   +deltat*dfgdt(iq)                               ! sice
       egice   =egice   +deltat*degdt(iq)                               ! sice
       es = establ(tgg(iq,3))                                           ! sice
       constz=ps(iq)-es                                                 ! sice
       qsttg(iq)=.622*es/constz                                         ! sice
                                                                        ! sice
!      combine ice and leads contributions here                         ! sice
       eg(iq) =fracice(iq)*egice  + (1.-fracice(iq))*eg(iq)             ! sice
       fg(iq) = fracice(iq)*fgice + (1.-fracice(iq))*fg(iq)             ! sice
       ri(iq) =fracice(iq)*ri_ice + (1.-fracice(iq))*ri(iq)  ! for scrnout
       zo(iq) =fracice(iq)*zoice  + (1.-fracice(iq))*zo(iq)  ! for scrnout
       cduv(iq) =fracice(iq)*af(iq)*fm + (1.-fracice(iq))*cduv(iq)      ! sice                                                                                 ! sice
       ustar(iq) = sqrt(vmod(iq)*cduv(iq))                              ! sice
c      N.B. potential evaporation is now eg+eg2                         ! sice
       epot(iq) =fracice(iq)*epotice + (1.-fracice(iq))*epot(iq)        ! sice
       tss(iq) = fracice(iq)*tgg(iq,3)+(1.-fracice(iq))*tgg(iq,2)       ! 2004
c      Surface stresses taux, tauy: diagnostic only - unstaggered now   ! sice
       taux(iq)=rho(iq)*cduv(iq)*u(iq,1)                                ! sice
       tauy(iq)=rho(iq)*cduv(iq)*v(iq,1)                                ! sice
       if(ntest.eq.1.and.iq.eq.idjd)then                                ! sice
         print *,'ri,vmod,cduv ',ri(iq),vmod(iq),cduv(iq)               ! sice
         print *,'tss,tgg2,tgg3 ',tss(iq),tgg(iq,2),tgg(iq,3)           ! sice
         print *,'theta,t1,deltat ',theta(iq),t(iq,1),deltat            ! sice
         print *,'b1,ga,gbot ',b1,ga(iq),gbot                           ! sice
         print *,'fg,fgice,factch ',fg(iq),fgice,factch(iq)             ! sice
         print *,'eg,egice,ustar ',eg(iq),egice,ustar(iq)               ! sice
       endif   ! (ntest.eq.1.and.iq.eq.idjd)                            ! sice
      enddo       ! ip=ipland+1,ipsice                                  ! sice
      if(ntest.eq.2)print *,'after sice loop'

c----------------------------------------------------------------------

!cdir nodep
      do ip=1,ipland  ! all land points in this shared loop             ! land
c      fh itself was only used outside this loop in sib0 (jlm)          ! land
       iq=iperm(ip)                                                     ! land
       zobg=zobgin                                                      ! land
       es = establ(tss(iq))                                             ! land
       constz=ps(iq)-es                                                 ! land
       qsttg(iq)=       .622*es/constz     ! only used in scrnout?      ! land
c      factch is sqrt(zo/zt) for land use in unstable fh                ! land
       factch(iq)=sqrt(7.4)                                             ! land
       if(snowd(iq).gt.0.)then    ! Fri  12-06-1996 (with soilsnow too) ! land
!        reduce zo over snow; done in darlam & globpe on 12-06-1996     ! land
!        zo=max(zo -.001*snowd(iq), .01)  ! wrongly till 24/2/97        ! land
!        following line is bit simpler than csiro9                      ! land
         zo(iq)=max(zolnd(iq) -.001*snowd(iq), .01)                     ! land
         zologsnw=log(zmin/zo(iq))                                      ! land
         zobg=max(zobgin -snowd(iq)*0.00976/12., 0.00024)               ! land
         alzz=log(zmin/zobg)                                            ! land
         if(nblend.eq.1)then ! blended zo for momentum: reduce for snow?! land
           afland=(vkar/((1.-sigmf(iq))*alzz+sigmf(iq)*zologsnw))**2    ! land
         else                                                           ! land
           afland=(vkar/zologsnw)**2                                    ! land
         endif   !   (nblend.eq.1)                                      ! land
         aftland=vkar**2/( zologsnw * (2.+zologsnw) )                   ! land
       else  ! land but not snow                                        ! land
         zo(iq)=zolnd(iq)                                               ! land
         alzz=alzzin                                                    ! land
         if(nblend.eq.1)then  ! blended zo for momentum                 ! land
           afland=(vkar/                                                ! land
     .               ((1.-sigmf(iq))*alzz+sigmf(iq)*zolog(iq)))**2      ! land
         else    ! non-blended zo for momentum                          ! land
           afland=(vkar/zolog(iq))**2                                   ! land
         endif   ! (nblend.eq.1)                                        ! land
         aftland=vkar**2/( zolog(iq) * (2.+zolog(iq)) )                 ! land
       endif     ! (snowd(iq).gt.0.)                                    ! land
       aft(iq)=aftland                                                  ! land
       aftlandg=vkar**2/( alzz * (2.+alzz) )                            ! land
c      lgwd>0 enhances cduv (momentum) over orog under (stable & unst) condns
       if(lgwd.gt.0)then
         af(iq)=afland+helo(iq)       ! jlm special gwd4b               ! land
       else
         af(iq)=afland                                                  ! land
       endif
       
	if(ntaft.eq.3.and.ktau.gt.1)then
!        do vegetation calulation for fh	
         xx(iq)=grav*zmin*(1.-tgf(iq)*srcp/t(iq,1)) ! actually otgf     ! land
         ri(iq)=xx(iq)/vmod(iq)**2                                      ! land
         if(xx(iq).gt.0.)then                                           ! land
           fh(iq)=vmod(iq)*max(fmroot*fmroot,1./(1.+bprm*ri(iq))**2)    ! land
         else                                                           ! land
           root=sqrt(-xx(iq)*zmin/zo(iq))  ! ignoring blending here     ! land
c          Now heat ; allow for smaller zo via aft and factch           ! land
           denha=vmod(iq)+chs*2.*bprm*factch(iq)*aft(iq)*root           ! land
           fh(iq)=vmod(iq)-(2.*bprm *xx(iq))/denha                      ! land
         endif                                                          ! land
         taftfh(iq)=aft(iq)*fh(iq)       ! uses fmroot above, for sib3  ! land 
	endif   ! (ntaft.eq.3.and.ktau.gt.1)
                                                                        ! land
c      Having settled on zo (and thus af) now do actual fh and fm calcs ! land
       xx(iq)=grav*zmin*(1.-tss(iq)*srcp/t(iq,1))                       ! land
       ri(iq)=xx(iq)/vmod(iq)**2                                        ! land
       if(xx(iq).gt.0.)then                                             ! land
         fm=vmod(iq)*max(fmroot*fmroot,1./(1.+bprm*ri(iq))**2)          ! land
         fh(iq)=fm                                                      ! land
         fhbg=fh(iq)                                                    ! land
       else                                                             ! land
         root=sqrt(-xx(iq)*zmin/zo(iq))  ! ignoring blending here       ! land
c        First do momentum                                              ! land
         denma=vmod(iq)+cms*2.*bprm*af(iq)*root                         ! land
         fm=vmod(iq)-(2.*bprm *xx(iq))/denma                            ! land
c        n.b. fm denotes ustar**2/(vmod(iq)*af)                         ! land
c        Now heat ; allow for smaller zo via aft and factch             ! land
         denha=vmod(iq)+chs*2.*bprm*factch(iq)*aft(iq)*root             ! land
         fh(iq)=vmod(iq)-(2.*bprm *xx(iq))/denha                        ! land
         rootbg=sqrt(-xx(iq)*zmin/zobg)                                 ! land
         denhabg=vmod(iq)+chs*2.*bprm*factch(iq)*aftlandg*rootbg        ! land
         fhbg=vmod(iq)-(2.*bprm *xx(iq))/denhabg                        ! land
       endif                                                            ! land
       taftfhg_temp(iq)=aftlandg*fhbg  ! uses fmroot above, for sib3    ! land
       taftfhg(iq)=aftlandg*fhbg ! value used for ntaft=3 (may need improving)
                                                                        ! land
c      cduv is now drag coeff *vmod                                     ! land
       cduv(iq) =af(iq)*fm                                              ! land
       ustar(iq) = sqrt(vmod(iq)*cduv(iq))                              ! land
c      cdtq(iq) =aft(iq)*fh(iq)                                         ! land
c      Surface stresses taux, tauy: diagnostic only - unstaggered now   ! land
       taux(iq)=rho(iq)*cduv(iq)*u(iq,1)                                ! land
       tauy(iq)=rho(iq)*cduv(iq)*v(iq,1)                                ! land
      enddo     ! ip=1,ipland                                           ! land
      if(ntest.eq.2)print *,'after land loop'
 
       if(ntest.eq.1.and.iq.eq.idjd)then                                ! land
         print *,'in main land loop'                                    ! land
         print *,'zmin,zobg ',zmin,zobg                                 ! land
         print *,'afland,aftland,alzz ',afland,aftland,alzz             ! land
         print *,'af,xx,vmod,es ',af(iq),xx(iq),vmod(iq),es             ! land
         print *,'tss,theta,t1 ',tss(iq),theta(iq),t(iq,1)              ! land
         print *,'aft,chnsea,fm,fh,rho,conh '                           ! land
     .           ,aft(iq),chnsea,fm,fh(iq),rho(iq),conh                 ! land
         print *,'ri,vmod,cduv,fg ',
     .            ri(iq),vmod(iq),cduv(iq),fg(iq)                       ! land
       endif                                                            ! land
      if(ntaft.eq.0.or.ktau.eq.1)then
        do iq=1,ifull  ! will only use land values
         taftfh(iq)=aft(iq)*fh(iq) ! uses fmroot above                  ! land
         taftfhg(iq)=taftfhg_temp(iq)
        enddo
      elseif(ntaft.eq.1.or.ntaft.eq.2)then
        do iq=1,ifull  ! will only use land values
         thnew=aft(iq)*fh(iq) ! uses fmroot above                       ! land
         thgnew=taftfhg_temp(iq)
         if(ntaft.eq.1)then
           if(thnew.gt.2.*taftfh(iq).or.thnew.lt..5*taftfh(iq))then
             taftfh(iq)=.5*(thnew+taftfh(iq))
           else
             taftfh(iq)=thnew
           endif
           if(thgnew.gt.2.*taftfhg(iq).or.thgnew.lt..5*taftfhg(iq))then
             taftfhg(iq)=.5*(thgnew+taftfhg(iq))
           else
             taftfhg(iq)=thgnew
           endif
         endif  ! (ntaft.eq.1)
         if(ntaft.eq.2)then    ! preferred faster option
           thnewa=min(thnew,
     .                max(2.*taftfh(iq),.5*(thnew+taftfh(iq))))
           taftfh(iq)=max(thnewa,
     .                    min(.5*taftfh(iq),.5*(thnew+taftfh(iq))))
           thgnewa=min(thgnew,
     .                 max(2.*taftfhg(iq),.5*(thgnew+taftfhg(iq))))
           taftfhg(iq)=max(thgnewa,
     .                 min(.5*taftfhg(iq),.5*(thgnew+taftfhg(iq))))
         endif  ! (ntaft.eq.2)
        enddo
      endif  ! (ntaft.eq.0.or.ktau.eq.1)  .. else ..

c ----------------------------------------------------------------------

      if(nsib.eq.1.or.nsib.eq.3)then
        if(ktau.eq.1)then  !To initialize new nsib=1/3 run (or restart) only
          print *,'ipland,ipsice,ipsea0,ipsea in sflux: ',
     .             ipland,ipsice,ipsea0,ipsea
!cdir nodep
          do ip=1,ipland
           iq=iperm(ip)
           tgf(iq) = t(iq,1)  ! was tss(iq)
!          tgg(iq,1) = tss(iq)
           rmc(iq) = 0.
!          ssdn(iq,1)=100.
          enddo    ! ip=1,ipland
        endif  ! if(ktau.eq.1)
        if(nsib.eq.3)call sib3(nalpha)
      endif     !  (nsib.eq.1 or 3)

c ----------------------------------------------------------------------

c    end of all land options
c     do iq=1,ifull
c      evap(iq)=evap(iq)+dt*eg(iq)/hl
c     enddo

      if(diag.or.ntest.gt.0)then
        if (mydiag) print *,'before call scrnout'
        call maxmin(t,' t',ktau,1.,kl)
      endif

!     always call scrnout from 19/9/02
      call scrnout(zo,ustar,fg,eg,factch,rho,wetfac,qsttg,   ! arrays
     .       qgscrn,tscrn,uscrn,u10,scrrel,                 ! arrays
     .       bprm,cms,chs,fmroot,nalpha)

c***  end of surface updating loop

      if(diag.or.ntest.eq.1)then
         if ( mydiag ) then
            print *,'at end of sflux, after call scrnout'
            print *,'slwa,rdg,eg,fg ',
     .               slwa(idjd),rdg(idjd),eg(idjd),fg(idjd)
            print *,'tscrn,cduv,zolnd ',
     .               tscrn(idjd),cduv(idjd),zolnd(idjd)
            print *,'degdt,dfgdt,snowd,sice ',
     .               degdt(idjd),dfgdt(idjd),snowd(idjd),sice(idjd)
            print *,'u1,v1,qg1 ',u(idjd,1),v(idjd,1),qg(idjd,1)
            print *,'w,w2,condx ',
     .               wb(idjd,1),wb(idjd,ms),condx(idjd)
            print *,'t1,tss,tgg_2,tgg_ms ',
     .               t(idjd,1),tss(idjd),tgg(idjd,2),tgg(idjd,ms)
         end if
        call maxmin(tscrn,'tc',ktau,1.,1)
        call maxmin(scrrel,'rh',ktau,1.,1)
      endif
      if(ntest.eq.4.and.ktau.eq.10)then
	 do iq=1,ifull
	  if(.not.land(iq))then
           write(45,'(2g13.4)') sqrt(u(iq,1)**2+v(iq,1)**2),fg(iq)
           write(46,'(2g13.4)') sqrt(u(iq,1)**2+v(iq,1)**2),eg(iq)
	  endif
         write(47,'(2g13.4)') sqrt(u(iq,1)**2+v(iq,1)**2),eg(iq)
        enddo
      endif
      return
      end

      subroutine sib3(nalpha)     ! new version of sib1 with soilsnowv
      use cc_mpi
      parameter (ntest=0) ! ntest= 0 for diags off; ntest= 1 for diags on
!                           N.B. may need vsafe for correct diags
      parameter (itnmeth=5) ! 0 for original N_R iteration method
!     parameter (nsigmf=1)  ! 0 for original tsigmf usage in sib3; prefer 1
!     parameter (nbarewet=2)  ! 0 for original bare-soil-wetfac; 1 simplest; 2 for jlm
      parameter (nstomata=1)  ! 0 for original; 1 for tropical C4
      parameter (newfgf=0)    ! 0 for original; 1 with tscrn; 2 with aft too
      parameter (ndiag_arr=0) ! 0 off; 1 for diag arrays on
      parameter (neva=0)    ! neva= 0 for diags off; neva= 1 for eva's diags on
      include 'newmpar.h'
      include 'arrays.h'
      include 'const_phys.h'
      include 'dates.h' ! ktime,kdate,timer,timeg,xg,yg
      include 'extraout.h'
      include 'latlong.h'  ! rlatt,rlongg
      include 'liqwpar.h'  ! qfg,qlg
      include 'morepbl.h'
      include 'nsibd.h'    ! rsmin,ivegt,sigmf,tgf,ssdn,res,rmc,tsigmf
      include 'parm.h'
      include 'pbl.h'
      include 'permsurf.h'
      include 'scamdim.h'  ! dimension of patches
      include 'screen.h'   ! tscrn etc
      include 'sigs.h'
      include 'soil.h'     ! ... zmin zolnd zolog sice alb
      include 'soilv.h'    ! ssat, clay,..
      include 'soilsnow.h' ! new soil arrays for scam - tgg too
      include 'tracers.h'  ! ngas, nllp, ntrac
      include 'trcom2.h'   ! nstn,slat,slon,istn,jstn
      common/nsib/nbarewet,nsigmf
      common/tafts/taftfh(ifull),taftfhg(ifull)
      common/work2/dirad(ifull),dfgdt(ifull),degdt(ifull)
     . ,wetfac(ifull),degdw(ifull),cie(ifull)
     . ,factch(ifull),qsttg(ifull),rho(ifull),zo(ifull)
     . ,aft(ifull),fh(ifull),ri(ifull),theta(ifull)
     . ,gamm(ifull),rg(ifull),vmod(ifull),dgdtg(ifull)
      common/work3/egg(ifull),evapxf(ifull),Ewww(ifull),fgf(ifull),
     . fgg(ifull),ggflux(ifull),rdg(ifull),rgg(ifull),residf(ifull),
     . ga(ifull),condxpr(ifull),fev(ifull),fes(ifull),
     . ism(ifull),fwtop(ifull),epot(ifull),
     . extin(ifull),af(ifull),spare1(ifull),xx(ifull),
     . dum3(5*ijk-20*ifull)
      common/work3c/airr(ifull),cc(ifull),condxg(ifull),delta_tx(ifull),
     . evapfb1(ifull),evapfb2(ifull),evapfb3(ifull),evapfb4(ifull),
     . evapfb5(ifull),evapfb1a(ifull),evapfb2a(ifull),evapfb3a(ifull),
     . evapfb4a(ifull),evapfb5a(ifull),otgf(ifull),rmcmax(ifull),
     . tgfnew(ifull),evapfb(ijk-17*ifull)  ! to allow > 18 levels
      common/work3d/dqsttg(ifull),tstom(ifull),rlai(ifull),
     .   cls(ifull),omc(ifull),dum3d(ijk-5*ifull)  ! allows L9
      include 'establ.h'

!     fle(isoil,w)=(w-swilt(isoil))/(sfc(isoil)-swilt(isoil))           !0 Eva's
!     fle(isoil,w)= w/ssat(isoil)                                       !1 simplest bare
!     fle(isoil,w)= (w-frac*max( w,swilt(isoil) ))/                     !2 jlm sugg. bare
!    .               (ssat(isoil)*(1.-frac))                            !2 jlm sugg. bare
!     fle(isoil,w)=10.*((w-ssoil(isoil))/ssat(isoil)+.1))               ! an old one
!     fle(isoil,w)= w/sfc(isoil)                                        ! jlm for PIRCS
!     fle(isoil,w)=(w-frac*swilt(isoil))/(sfc(isoil)-frac*swilt(isoil)) ! jlm special

!***  N.B. nrungcm=2 fix-up for Mark 2 wb moved to indata 6/11/00
      if(ktau.eq.1.and.nrungcm.eq.3)then  ! for nsib runs for PIRCS
!cdir nodep
        do ip=1,ipland  ! all land points in this nsib=1/3 loop
         iq=iperm(ip)
c        code where eta provided (e.g. for PIRCS)
         do layer=2,ms
          wb(iq,layer)=wb(iq,1)   ! w, w2 and wb all same initially for PIRCS
         enddo       !  layer=1,ms
        enddo        !  ip=1,ipland
      endif          !  (ktau.eq.1.and.nrungcm.eq.3)

      if(ktau.eq.1)then
!cdir nodep
        do ip=1,ipland  ! all land points in this nsib=3 loop
         iq=iperm(ip)
         tscrn(iq)=theta(iq)  ! first guess, needed for newfgf=1
!eak         real egst(ifull),fgst(ifull),precipst(ifull),evapst(ifull)
!eak         real runoffst(ifull),rnetst(ifull),tssst(ifull)
!eak         real wbavst(ifull),wb1m(ifull)
!eak         egst(iq)=0.
!eak         fgst(iq)=0.
!eak         precipst(iq)=0.
!eak         evapst(iq)=0.
!eak         runoffst(iq)=0.
!eak         rnetst(iq)=0.
!eak         tssst(iq)=0.
!eak         wbavst(iq)=0.
        enddo         ! ip=1,ipland

        if ( mydiag ) then
           iveg=ivegt(idjd)
           isoil = isoilm(idjd)
           tsoil=0.5*(0.3333*tgg(idjd,2)+0.6667*tgg(idjd,3)
     &             +0.95*tgg(idjd,4) +  0.05*tgg(idjd,5))
           ftsoil=max( 0. , 1.-.0016*(298.-tsoil)*max(0.,298.-tsoil) )
           tsigmf(idjd)=max(0.,sigmf(idjd)-scveg44(iveg)*(1.-ftsoil))
           print *,'nsib,iveg,isoil,nalpha,nsigmf,tsigmf,newfgf ',
     .           nsib,iveg,isoil,nalpha,nsigmf,tsigmf(idjd),newfgf
           print *,'tsoil,ftsoil,scveg44,sigmf ',
     .           tsoil,ftsoil,scveg44(iveg),sigmf(idjd)
           print *,'swilt,sfc,wb1-6 ',
     .           swilt(isoil),sfc(isoil),(wb(idjd,k),k=1,ms)
           rlai_d= max(.1,rlaim44(iveg)-slveg44(iveg)*(1.-ftsoil))
           srlai=rlai_d+rlais44(iveg)
           rsmin(idjd)=rsunc44(iveg)/rlai_d ! now always done
           print *,'rlai,srlai,rsmin ',rlai_d,srlai,rsmin(idjd)
        end if
      endif           ! (ktau.eq.1)
      rmstep=dt/60.

c      due to a bug in the SX5 compiler even for cvsafe, 
c      following loop needed to be split off from the one after it
!cdir nodep
      do ip=1,ipland  ! all land points in this nsib=3 loop
       iq=iperm(ip)
       iveg=ivegt(iq)
c      evaluate seasonal impact and the snow depth
c      impact on the fractional vegetation cover
       tstom(iq)=298.
       if(iveg.eq.6+31)tstom(iq)=302.
       if(iveg.ge.10.and.iveg.le.21.and.
     .    abs(rlatt(iq)*180./pi).lt.25.)tstom(iq)=302.
       if(ntest.eq.1.and.iq.eq.idjd) then
         print *,'in sib3a ip,iq,idjd,iveg ',ip,iq,idjd,iveg
         print*,'snowd,zo,zolnd,tstom ',
     .           snowd(iq),zo(iq),zolnd(iq),tstom(iq)
       endif ! ntest
      enddo         ! ip=1,ipland
c     print *,'after 1st loop'

!cdir nodep
      do ip=1,ipland  ! all land points in this nsib=3 loop
       iq=iperm(ip)
       iveg=ivegt(iq)
       tsoil=min(tstom(iq), .5*(.3333*tgg(iq,2)+.6667*tgg(iq,3)
     &            +.95*tgg(iq,4) + .05*tgg(iq,5)))
       ftsoil=max(0.,1.-.0016*(tstom(iq)-tsoil)**2)
c         which is same as:  ftsoil=max(0.,1.-.0016*(tstom-tsoil)**2)
c                            if( tsoil .ge. tstom ) ftsoil=1.
       rlai(iq)=  max(.1,rlaim44(iveg)-slveg44(iveg)*(1.-ftsoil))
       tsigmf(iq)=max(.001,  sigmf(iq)-scveg44(iveg)*(1.-ftsoil)) ! 13/5/03
       if(ntest.eq.1.and.iq.eq.idjd) then
         print *,'in sib3b ip,iq,idjd,iveg ',ip,iq,idjd,iveg
         print*,'iveg,sigmf(iq),tsigmfa ',iveg,sigmf(iq),tsigmf(iq)
         print*,'rlaim44,tsoil,ftsoil ',
     .           rlaim44(iveg),tsoil,ftsoil
         print*,'scveg44,snowd,zo,zolnd,tstom ',
     .          scveg44(iveg),snowd(iq),zo(iq),zolnd(iq),tstom(iq)
         print*,'w2,rlai ',wb(iq,ms),rlai(iq)
       endif ! ntest
       enddo         ! ip=1,ipland
c      print *,'after 2nd loop'
c 
!cdir nodep
       do ip=1,ipland  
        iq=iperm(ip)
       tsigmf(iq)=(1.-snowd(iq)/
     .            (snowd(iq)+5.*100.*zo(iq)))*tsigmf(iq)
!      extin(iq)=exp(-0.6*max(1.,rlai(iq)))  ! good approx uses next 2 (jlm)
	xxx=.6*max(1.,rlai(iq))
	extin(iq)=1.-xxx/(1. +.5*xxx +xxx*xxx/12.) 
       if(ntest.eq.1.and.iq.eq.idjd) then
         print *,'in sib3c ip,iq,idjd,iveg ',ip,iq,idjd,ivegt(iq)
         print*,'iveg,sigmf(iq),tsigmf ',ivegt(iq),sigmf(iq),tsigmf(iq)
         print*,'scveg44,snowd,zo,zolnd,tstom ',
     .          scveg44(iveg),snowd(iq),zo(iq),zolnd(iq),tstom(iq)
         print *,'alb,sgsave ',alb(iq),sgsave(iq)
         print*,'w2,rlai,extin ',wb(iq,ms),rlai(iq),extin(iq)
       endif ! ntest
c      bare ground calculation
!      if(isflag(iq).eq.1)then
!        tgss=tggsn(iq,1)
!      else
!        tgss=tgg(iq,1)
!      endif       ! replace these by:
       tgss=isflag(iq)*tggsn(iq,1) + (1-isflag(iq))*tgg(iq,1)
       esattg=establ(tgss)
       qsttg(iq)=.622*esattg/(ps(iq)-esattg)
       tgss2=tgss*tgss
       dqsttg(iq)=qsttg(iq)*ps(iq)*hlars/((ps(iq)-esattg)*tgss2)
       rgg(iq) =  stefbo*tgss2**2   ! i.e. stefbo*tgss**4
       dirad(iq)=4.*rgg(iq)/tgss
c      sensible heat flux
       dfgdt(iq)=taftfhg(iq)*rho(iq)*cp
       fgg(iq)=dfgdt(iq)*(tgss-theta(iq))
      enddo         ! ip=1,ipland
 
c     print *,'before wetfac'

      if(nbarewet.eq.0)then  ! original Eva's, same as NCAR
!cdir nodep
        do ip=1,ipland  ! all land points in this nsib=3 loop
         iq=iperm(ip)
         isoil = isoilm(iq)
         fle=(wb(iq,1)-swilt(isoil))/(sfc(isoil)-swilt(isoil))          
         wetfac(iq)=max( 0.,min(1.,fle) )
        enddo   ! ip=1,ipland
      endif     ! (nbarewet.eq.0)

      if(nbarewet.eq.1)then  ! simplest bare
!cdir nodep
        do ip=1,ipland  ! all land points in this nsib=3 loop
         iq=iperm(ip)
         isoil = isoilm(iq)
         fle= wb(iq,1)/ssat(isoil)                                   
         wetfac(iq)=max( 0.,min(1.,fle) )
        enddo   ! ip=1,ipland
      endif     ! (nbarewet.eq.1)

      if(nbarewet.eq.2)then  ! jlm suggestion
!cdir nodep
        do ip=1,ipland  ! all land points in this nsib=3 loop
         iq=iperm(ip)
         isoil = isoilm(iq)
         frac=max(.01,tsigmf(iq))     ! jlm special for fle
         fle= (wb(iq,1)-frac*max( wb(iq,1),swilt(isoil) ))/         
     .               (ssat(isoil)*(1.-frac))                         
         wetfac(iq)=max( 0.,min(1.,fle) )
        enddo   ! ip=1,ipland
      endif     ! (nbarewet.eq.2)

      if(nbarewet.eq.3)then  ! jlm suggestion
!cdir nodep
        do ip=1,ipland  ! all land points in this nsib=3 loop
         iq=iperm(ip)
         isoil = isoilm(iq)
         frac=max(.01,tsigmf(iq))     ! jlm special for fle
         fle= (wb(iq,1)-frac*swilt(isoil) )/         
     .               (sfc(isoil)-frac*swilt(isoil))                         
         wetfac(iq)=max( 0.,min(1.,fle) )
        enddo   ! ip=1,ipland
      endif     ! (nbarewet.eq.3)

      if(nbarewet.eq.4)then  ! jlm, similar to Noilhan & Planton
!cdir nodep
        do ip=1,ipland  ! all land points in this nsib=3 loop
         iq=iperm(ip)
         isoil = isoilm(iq)
         fle=min( 1.,wb(iq,1)/sfc(isoil) )         
         wetfac(iq)=fle*fle*(3.-2.*fle)
        enddo   ! ip=1,ipland
      endif     ! (nbarewet.eq.4)

      if(nbarewet.eq.5)then  ! jlm, similar to Noilhan & Planton with swilt
!cdir nodep
        do ip=1,ipland  ! all land points in this nsib=3 loop
         iq=iperm(ip)
         isoil = isoilm(iq)
         fle=max( 0.,min( 1.,
     .             (wb(iq,1)-swilt(isoil))/(sfc(isoil)-swilt(isoil)) ) )         
         wetfac(iq)=fle*fle*(3.-2.*fle)
        enddo   ! ip=1,ipland
      endif     ! (nbarewet.eq.5)

      if(nbarewet.eq.6)then  ! newer jlm
!cdir nodep
        do ip=1,ipland  ! all land points in this nsib=3 loop
         iq=iperm(ip)
         isoil = isoilm(iq)
         fle= max( 0.,min(1.,wb(iq,1)/ssat(isoil)) )					
         wetfac(iq)=fle*fle*(2.2-1.2*fle)  ! .4 for fle=.5
        enddo   ! ip=1,ipland
      endif     ! (nbarewet.eq.6)

      if(nbarewet.eq.7)then  ! newest piecewise jlm (use with nalpha=1, beta)
!cdir nodep
        do ip=1,ipland  ! all land points in this nsib=3 loop
         iq=iperm(ip)
         isoil = isoilm(iq)
c         betaa=min(.2,.2*(wb(iq,1)/swilt(isoil))**2)
c         betab=max(min(betaa,.2),(.2*(sfc(isoil)-wb(iq,1))+
c     .          .8*(wb(iq,1)-swilt(isoil)))/(sfc(isoil)-swilt(isoil)))  
c         betac=max(.8,(.8*(ssat(isoil)-wb(iq,1))+  
c     .             (wb(iq,1)-sfc(isoil)))/(ssat(isoil)-sfc(isoil)))
c         wetfac(iq)=min(max(betaa,betab),betac)
         wetfac(iq)=.2*(wb(iq,1)/swilt(isoil))**2
	  if(wb(iq,1).gt.swilt(isoil))then
          wetfac(iq)=(.2*(sfc(isoil)-wb(iq,1))+
     .             .8*(wb(iq,1)-swilt(isoil)))/(sfc(isoil)-swilt(isoil))
         endif
	  if(wb(iq,1).gt.sfc(isoil))then
           wetfac(iq)=(.8*(ssat(isoil)-wb(iq,1))+  
     .                   (wb(iq,1)-sfc(isoil)))/(ssat(isoil)-sfc(isoil))  
         endif
c	  if(wb(iq,1).gt.swilt(isoil).and.betaa.gt.betab)then
c	    print *,'iq,wb,swilt,sfc,ssat ',
c     .              iq,wb(iq,1),swilt(isoil),sfc(isoil),ssat(isoil)
c	    print *,'betaa,betab,betac,wetfac ',
c     .              betaa,betab,betac,wetfac(iq)
c         endif
c	  if(wetfac(iq).lt.0.)then
c	    print *,'iq,wb,swilt,sfc,ssat ',
c     .              iq,wb(iq,1),swilt(isoil),sfc(isoil),ssat(isoil)
c	    print *,'betaa,betab,betac,wetfac ',
c     .              betaa,betab,betac,wetfac(iq)
c         endif
        enddo   ! ip=1,ipland
      endif     ! (nbarewet.eq.7)
 
c     print *,'before nalpha'

      if(nalpha.eq.1)then    ! beta scheme
!cdir nodep
        do ip=1,ipland  ! all land points in this nsib=3 loop
         iq=iperm(ip)
         isoil = isoilm(iq)
         conw_fh=rho(iq)*taftfhg(iq)*hl    ! taftfhg(iq)=aftlandg*fhbg
         qgtot=qg(iq,1)+qfg(iq,1)+qlg(iq,1)
         epot(iq) = conw_fh*(qsttg(iq)-qgtot)
         egg(iq)=wetfac(iq)*epot(iq)
         degdt(iq)=wetfac(iq)*conw_fh*dqsttg(iq)
!        degdw(iq)=conw_fh*(qsttg(iq)-qg(iq,1))/ssat(isoil)  ! not used in sib3
        enddo   ! ip=1,ipland
      else
!       following is alpha scheme
!cdir nodep
        do ip=1,ipland  ! all land points in this nsib=3 loop
         iq=iperm(ip)
         isoil = isoilm(iq)
         conw_fh=rho(iq)*taftfhg(iq)*hl    ! taftfhg(iq)=aftlandg*fhbg
         qgtot=qg(iq,1)+qfg(iq,1)+qlg(iq,1)
         qtgnet=  qsttg(iq)*wetfac(iq) -qgtot
         qtgair=  qsttg(iq)*wetfac(iq)-max(qtgnet,.1*qtgnet)
         eg2=   -conw_fh*qtgair
         eg1=    conw_fh*qsttg(iq)
!        degdw(iq)=eg1/ssat(isoil)                         ! not used in sib3
c        evaporation from the bare ground
         egg(iq)=eg1*wetfac(iq) +eg2
         epot(iq) = conw_fh*(qsttg(iq)-qgtot)
         deg=wetfac(iq)*conw_fh*dqsttg(iq)
c        following reduces degdt by factor of 10 for dew
         degdt(iq)=.55*deg+sign(.45*deg,qtgnet)
        enddo   ! ip=1,ipland
      endif    ! (nalpha.eq.1) .. else ..
        if(ntest.eq.1.and.mydiag)then  ! SX5: don't put this test in above loop
          iq=idjd
          print *,'epot,egg,tgg1,snowd ',
     .             epot(iq),egg(iq),tgg(iq,1),snowd(iq)
          isoil = isoilm(iq)
          conw_fh=rho(iq)*taftfhg(iq)*hl    ! taftfhg(iq)=aftlandg*fhbg
          qgtot=qg(iq,1)+qfg(iq,1)+qlg(iq,1)
          qtgnet=  qsttg(iq)*wetfac(iq) -qgtot
          qtgair=  qsttg(iq)*wetfac(iq)-max(qtgnet,.1*qtgnet)
          eg2=   -conw_fh*qtgair
          eg1=    conw_fh*qsttg(iq)
!         degdw(iq)=eg1/ssat(isoil)                         ! not used in sib3
c         evaporation from the bare ground
          egg(iq)=eg1*wetfac(iq) +eg2
          epot(iq) = conw_fh*(qsttg(iq)-qgtot)
          egg_alph1=wetfac(iq)*epot(iq)
          print *,'then iq,isoil,conw_fh,qsttg,qtgair ',
     .             iq,isoil,conw_fh,qsttg(iq),qtgair
          print *,'eg1,eg2,wetfac ',eg1,eg2,wetfac(iq)
          print *,'epot,egg,egg_alph1 ',
     .             epot(iq),egg(iq),egg_alph1
        endif  ! (ntest.eq.1)

!cdir nodep
      do ip=1,ipland  ! all land points in this nsib=3 loop
       iq=iperm(ip)
       if(snowd(iq).gt.1.)then
         egg(iq)=epot(iq)
	 wetfac(iq)=1.   ! added jlm 18/3/04 to fix qgscrn inconsistency
         cls(iq)=1.+hlf/hl
       else
         egg(iq)=min(egg(iq),wb(iq,1)*zse(1)*1000.*hl/dt)
         cls(iq)=1.
       endif  ! (snowd(iq).gt.1.)
      enddo   ! ip=1,ipland
 
c     print *,'before nsigmf'

      if(nsigmf.eq.0)then  ! original
!cdir nodep
        do ip=1,ipland  ! all land points in this nsib=3 loop
         iq=iperm(ip)
         ga(iq)=-slwa(iq)-rgg(iq)-fgg(iq)-cls(iq)*egg(iq)
         dgdtg(iq)=-dirad(iq)-dfgdt(iq)-cls(iq)*degdt(iq)
        enddo   ! ip=1,ipland
      endif     ! (nsigmf.eq.0)

      if(nsigmf.eq.1)then  ! jlm preferred
!cdir nodep
        do ip=1,ipland  ! all land points in this nsib=3 loop
        iq=iperm(ip)
        ga(iq)=(-slwa(iq)-rgg(iq)-fgg(iq)-cls(iq)*egg(iq))*
     .                                  (1.-tsigmf(iq))
!       dgtdg is used in soilsnow
        dgdtg(iq)=-(dirad(iq)+dfgdt(iq)+cls(iq)*degdt(iq))*
     .                                  (1.-tsigmf(iq))
        enddo   ! ip=1,ipland
      endif     ! (nsigmf.eq.1)

! ----------------------------------------------
      if(itnmeth.eq.0) then  ! old, not vectorized
!cdir nodep
       do ip=1,ipland  ! all land points in this nsib=3 loop
        iq=iperm(ip)
        isoil = isoilm(iq)
        iveg=ivegt(iq)
c                                             leaf area index
          srlai=rlai(iq)+rlais44(iveg)
          rsmin(iq) = rsunc44(iveg)/rlai(iq)   ! now always done
c                                    components of the stomatal resistance
          sstar = 30.
          if( zo(iq) .lt. .5 ) sstar = 150.
          f= 1.1*sgsave(iq)/(rlai(iq)*sstar)
c         f= 1.1*slwa(iq)/(rlai(iq)*sstar)
          rsi = rsmin(iq) * rlai(iq)
          f1= (1.+f)/(f+rsi/5000.)
          den=sfc(isoil)-swilt(isoil)                          ! sib3
          wbav=(max(0.,froot(1)*(wb(iq,1)-swilt(isoil)))+
     &         max(0.,froot(2)*(wb(iq,2)-swilt(isoil)))+
     &         max(0.,froot(3)*(wb(iq,3)-swilt(isoil)))+
     &         max(0.,froot(4)*(wb(iq,4)-swilt(isoil)))+
     &         max(0.,froot(5)*(wb(iq,5)-swilt(isoil)))   )/den
          f2=max(1. , .5/ max( wbav,1.e-7)) ! N.B. this is equiv to next 2 (jlm)
c         f2=1.0
c         if(wbav.lt.0.5) f2=max(1.0 , 0.5/ max( wbav,1.e-7))
          f4=max(1.-.0016*(tstom(iq)-t(iq,1))**2 , .05) ! zero for delta_t=25
          airr(iq) = 1./taftfh(iq)
          cc(iq) =min(condx(iq) , 4./(1440./rmstep))  ! jlm speedup for 4 mm/day
c                       depth of the reservoir of water on the canopy
          rmcmax(iq) = max(0.5,srlai) * .1
          omc(iq) = rmc(iq)  ! value from previous timestep as starting guess
          f3=max(1.-.00025*(establ(t(iq,1))-qg(iq,1)*ps(iq)/.622), .05)
          res(iq)=max(30.,rsmin(iq)*f1*f2/(f3*f4))
          if(ntest.eq.1.and.iq.eq.idjd)then
           print *,'rlai,srlai,wbav,den ',rlai(iq),srlai,wbav,den
           print *,'f1,f2,f3,f4 ',f1,f2,f3,f4
           print *,'f,f124,rsi,res ',f,f1*f2/f4,rsi,res(iq)
          endif

          otgf(iq)=tgf(iq)
           do icount=1,5        ! original iteration method
c                                            transpiration
            rmc(iq) = omc(iq)
            esatf = establ(tgf(iq))
            qsatgf=.622*esatf/(ps(iq)-esatf)
c                                                     wet evaporation
            Ewww(iq) = rho(iq)*taftfh(iq) *(qsatgf-qg(iq,1))
            if(qsatgf.ge.qg(iq,1)) then
              Ewww(iq)  = min(rmc(iq)/dt , rmc(iq)*Ewww(iq)/rmcmax(iq) )
            endif         ! qsatgf.ge.qg(iq,1)
            rmc(iq)=omc(iq)+cc(iq) -Ewww(iq)*dt
c                         precipitation reaching the ground under the canopy
c                                            water interception by the canopy
            condxg(iq)=max(condx(iq)-cc(iq)+
     .                 max(0.,rmc(iq)-rmcmax(iq)),0.) ! keep here
            rmc(iq) = min( max(rmc(iq),0.), rmcmax(iq))
!           beta =  min(rmc(iq)/rmcmax(iq),1.)   ! not needed
            beta =      rmc(iq)/rmcmax(iq)
            if( qsatgf .lt. qg(iq,1) ) beta = 1.

            Etr=rho(iq)/(airr(iq) +res(iq))*(qsatgf-qg(iq,1))
            betetrdt =(1.-beta)*Etr*dt*tsigmf(iq)
            evapfb1(iq)=min(betetrdt*froot(1),max(0.,
     &                  (wb(iq,1)-swilt(isoil))*zse(1)*1000.))
            evapfb2(iq)=min(betetrdt*froot(2),max(0.,
     &                  (wb(iq,2)-swilt(isoil))*zse(2)*1000.))
            evapfb3(iq)=min(betetrdt*froot(3),max(0.,
     &                  (wb(iq,3)-swilt(isoil))*zse(3)*1000.))
            evapfb4(iq)=min(betetrdt*froot(4),max(0.,
     &                  (wb(iq,4)-swilt(isoil))*zse(4)*1000.))
            evapfb5(iq)=min(betetrdt*froot(5),max(0.,
     &                  (wb(iq,5)-swilt(isoil))*zse(5)*1000.))
            evapfb(iq)=(evapfb1(iq)+evapfb2(iq)+evapfb3(iq)+evapfb4(iq)+
     &                  evapfb5(iq))/tsigmf(iq)
            evapxf(iq) = (evapfb(iq)/dt + Ewww(iq))*hl
c           evapw =  Ewww(iq)*hl  ! not used
c                                              sensible heat flux
            prz = rho(iq)*taftfh(iq)*cp
            if(newfgf.eq.0)fgf(iq)=prz*(tgf(iq)-theta(iq))  ! original
            if(newfgf.eq.1)fgf(iq)=prz*(tgf(iq)-tscrn(iq))
            if(newfgf.eq.2)fgf(iq)=rho(iq)*aft(iq)*cp*
     .                                (tgf(iq)-tscrn(iq))
c           fgf(iq)=min(800.,fgf(iq)) !***********************
            rdg(iq) =  stefbo*tgf(iq)**4
!           residf(iq) = (-slwa(iq) -rdg(iq))*(1.-extin(iq))
!    .                    -fgf(iq) -evapxf(iq)
            residf(iq) = -slwa(iq) -rdg(iq) -fgf(iq) -evapxf(iq)
            if( abs(residf(iq)) .lt. 3. )go to 110
            dirad1 = 4.*rdg(iq)/tgf(iq)
c                                                     Calculate dE/dTg
            dqg=qsatgf*hlars/tgf(iq)**2
            devf= hl*rho(iq)*dqg*( (1.-beta)/(airr(iq) + res(iq))
     &                          +beta/airr(iq) ) ! re-factored by jlm
            residp = -(dirad1 + devf + prz)
            sstgf = tgf(iq)
            tgf(iq)=tgf(iq)-residf(iq)/residp
            if(ntest.eq.1.and.iq.eq.idjd)then
              print *,'icount,omc,rmc,otgf ',
     .                 icount,omc(iq),rmc(iq),otgf(iq)
              print *,'tfg,residf,evapxf,fgf ',
     .                 tgf(iq),residf(iq),evapxf(iq),fgf(iq)
            endif
            tgf(iq)=min(tgf(iq),sstgf+ 1.5)   ! jlm speedup
            tgf(iq)=max(tgf(iq),sstgf- 1.5)   ! jlm speedup

           enddo   !  icount=1,5
           if(ntest.eq.1)print *,'iq,otgf(iq),tgf,residf '
     .                           ,iq,otgf(iq),tgf(iq),residf(iq)
c          tgf(iq)=0.5*(otgf(iq)+tgf(iq))
110        continue
       enddo  !  ip=1,ipland
      endif  ! (itnmeth.eq.0) 


      if(itnmeth.gt.0) then  ! new, vectorized
!cdir nodep
       do ip=1,ipland  ! all land points in this nsib=3 loop
        iq=iperm(ip)
        isoil = isoilm(iq)
        iveg=ivegt(iq)
c                                          leaf area index
        srlai=rlai(iq)+rlais44(iveg)
        rsmin(iq) = rsunc44(iveg)/rlai(iq)   ! now always done
c                                  components of the stomatal resistance
!       sstar = 30.
!       if( zo(iq) .lt. .5 ) sstar = 150.
	 sstar=90.+sign(60.,.5-zo(iq))  ! equiv to above 2 lines
        f= 1.1*sgsave(iq)/(rlai(iq)*sstar)
c       f= 1.1*slwa(iq)/(rlai(iq)*sstar)
        rsi = rsmin(iq) * rlai(iq)
        f1= (1.+f)/(f+rsi/5000.)
        den=sfc(isoil)-swilt(isoil)                          ! sib3
        wbav=(max(0.,froot(1)*(wb(iq,1)-swilt(isoil)))+
     &        max(0.,froot(2)*(wb(iq,2)-swilt(isoil)))+
     &        max(0.,froot(3)*(wb(iq,3)-swilt(isoil)))+
     &        max(0.,froot(4)*(wb(iq,4)-swilt(isoil)))+
     &        max(0.,froot(5)*(wb(iq,5)-swilt(isoil)))   )/den
        f2=max(1. , .5/ max( wbav,1.e-7)) ! N.B. this is equiv to next 2 (jlm)
c       f2=1.0
c       if(wbav.lt.0.5) f2=max(1.0 , 0.5/ max( wbav,1.e-7))
        f4=max(1.-.0016*(tstom(iq)-t(iq,1))**2 , .05) ! zero for delta_t=25
        airr(iq) = 1./taftfh(iq)
        cc(iq) =min(condx(iq) , 4./(1440./rmstep))  ! jlm speedup
c                     depth of the reservoir of water on the canopy
        rmcmax(iq) = max(0.5,srlai) * .1
        omc(iq) = rmc(iq)  ! value from previous timestep as starting guess
        qgtot=qg(iq,1)+qfg(iq,1)+qlg(iq,1)
        f3=max(1.-.00025*(establ(t(iq,1))-qgtot*ps(iq)/.622)
     .                                                        ,.05)
        res(iq)=max(30.,rsmin(iq)*f1*f2/(f3*f4))
        if(ntest.eq.1.and.iq.eq.idjd)then
          print *,'rlai,srlai,wbav,den ',rlai(iq),srlai,wbav,den
          print *,'f1,f2,f3,f4 ',f1,f2,f3,f4
          print *,'f,f124,rsi,res ',f,f1*f2/f4,rsi,res(iq)
	  print *,'qg,qfg,qlg,qgtot ',qg(iq,1),qfg(iq,1),qlg(iq,1),qgtot
        endif
        otgf(iq)=tgf(iq)
        tgfnew(iq)=tgf(iq)
        delta_tx(iq)=10.   ! just to supply max change of 5 deg first time
        evapfb1a(iq)=max(0.,wb(iq,1)-swilt(isoil)) *zse(1)*1000.
        evapfb2a(iq)=max(0.,wb(iq,2)-swilt(isoil)) *zse(2)*1000.
        evapfb3a(iq)=max(0.,wb(iq,3)-swilt(isoil)) *zse(3)*1000.
        evapfb4a(iq)=max(0.,wb(iq,4)-swilt(isoil)) *zse(4)*1000.
        evapfb5a(iq)=max(0.,wb(iq,5)-swilt(isoil)) *zse(5)*1000.
       enddo   ! ip loop

       do icount=1,itnmeth     ! jlm new iteration
c                                            transpiration
!cdir nodep
        do ip=1,ipland  ! all land points in this nsib=3 loop
         iq=iperm(ip)
         esatf = establ(tgfnew(iq))
         qsatgf=.622*esatf/(ps(iq)-esatf)
c                                                  wet evaporation
         qgtot=qg(iq,1)+qfg(iq,1)+qlg(iq,1)
         Ewwwa = rho(iq) *(qsatgf-qgtot)/airr(iq)
!        if(qsatgf.ge.qg(iq,1)) then  ! no dew
!          Ewww(iq)  = min(omc(iq)/dt , Ewww(iq)*omc(iq)/rmcmax(iq) )
!        endif         ! qsatgf.ge.qg(iq,1)
!        above 3 lines are equivalent to:
c        Ewww(iq)=min(Ewww(iq),omc(iq)/dt,Ewww(iq)*omc(iq)/rmcmax(iq))
         rmcav=.5*(omc(iq)+rmc(iq))
         Ewww(iq)=min(Ewwwa,rmcav/dt,Ewwwa*rmcav/rmcmax(iq))
!                     dew(-ve), no_dew ,        no_dew
!        rmc is reservoir on leaf
         rmc(iq)=omc(iq)+cc(iq) -Ewww(iq)*dt
c                      precipitation reaching the ground under the canopy
c                                         water interception by the canopy
         condxg(iq)=max(condx(iq)-cc(iq)+max(0.,rmc(iq)-rmcmax(iq)),0.) ! keep 
         rmc(iq) = min( max(0.,rmc(iq)), rmcmax(iq))
!        beta =  min(rmc(iq)/rmcmax(iq),1.)   ! min not needed
         beta =      rmc(iq)/rmcmax(iq)
!!       if( qsatgf .lt. qg(iq,1) ) beta = 1.   ! i.e. dew
         Etr=rho(iq)*(qsatgf-qgtot)/(airr(iq) +res(iq))
         Etr=rho(iq)*max(0.,qsatgf-qgtot)/(airr(iq) +res(iq))  ! jlm
         betetrdt =(1.-beta)*Etr*dt*tsigmf(iq)   ! fixed 23/5/01
         evapfb1(iq)=min(betetrdt*froot(1),evapfb1a(iq))
         evapfb2(iq)=min(betetrdt*froot(2),evapfb2a(iq))
         evapfb3(iq)=min(betetrdt*froot(3),evapfb3a(iq))
         evapfb4(iq)=min(betetrdt*froot(4),evapfb4a(iq))
         evapfb5(iq)=min(betetrdt*froot(5),evapfb5a(iq))
         evapfb(iq)=(evapfb1(iq)+evapfb2(iq)+evapfb3(iq)+evapfb4(iq)+
     &               evapfb5(iq))/tsigmf(iq)
         evapxf(iq) = (evapfb(iq)/dt + Ewww(iq))*hl  ! converting to W/m**2
         prz = rho(iq)*cp*taftfh(iq)
         if(newfgf.eq.0)fgf(iq) = prz*(tgfnew(iq)-theta(iq))  ! original/usual
         if(newfgf.eq.1)fgf(iq) = prz*(tgfnew(iq)-tscrn(iq))
         if(newfgf.eq.2)fgf(iq)=rho(iq)*aft(iq)*cp*
     .                             (tgfnew(iq)-tscrn(iq))
!        for +ve flux, can look toward level 2 (anticipating mixed PBL)	  
c        if(newfgf.eq.3)fgf(iq)=prz*min(tgfnew(iq)-theta(iq),
c    &                                  max(0.,tgfnew(iq)-theta12(iq)) )
         rdg(iq) =  stefbo*tgfnew(iq)**4
         residf(iq) = -slwa(iq) - rdg(iq) - fgf(iq) - evapxf(iq)
         dirad1 = 4.*rdg(iq)/300.
!        next 2 expressions can be approximated without effects
c        dqg=qsatgf*hlars/300.**2
ca       devf= hl*rho(iq)*dqg*( (1.-beta)/(airr(iq) + res(iq))
ca   &                       +beta/airr(iq) ) ! re-factored by jlm
!        according to jlm prints, devf has only small effect
         devf= (hl*hlars/300.**2)*qsatgf*(1.-beta)/res(iq)
         delta_t0=residf(iq)/(dirad1 + devf + prz)
!!       delta_t=sign(min(abs(delta_t0),5.),delta_t0) ! max 5 deg 1st it
!!       following lines assist convergence sometimes
!!       if(icount.gt.1)then
!!         delta_t=sign(min(abs(delta_t),                    !  jlmnew
!!   .                     .5*abs(tgfnew(iq)-tgf(iq))) , delta_t0)
!!       endif      ! icount.gt.1
!        above few lines equivalent to next one:
         delta_t=sign(min(abs(delta_t0),.5*abs(delta_tx(iq))),delta_t0)
         tgfnew(iq)=tgfnew(iq)+delta_t
         delta_tx(iq)=tgfnew(iq)-otgf(iq)
!        following to limit change to 8 degrees as fh may be poor  May '04
         tgfnew(iq)=otgf(iq)+
     .              sign(min(abs(delta_tx(iq)),8.),delta_tx(iq))
	 enddo  ! ip loop
        if(ntest.eq.1.and.mydiag)then 
	   iq=idjd
           print *,'icount,iq,omc,rmc ',icount,iq,omc(iq),rmc(iq)
           print *,'theta,tscrn,slwa ',
     .             theta(iq),tscrn(iq),slwa(iq)
           print *,'taftfh,condxg ',taftfh(iq),condxg(iq)
c          print *,'Ewwwa,Ewww,Etr ',Ewwwa,Ewww(iq),Etr
           print *,'Ewww ',Ewww(iq)
           print *,'rdg,fgf,evapxf,evapfb ',
     .              rdg(iq),fgf(iq),evapxf(iq),evapfb(iq)
           esatf = establ(tgfnew(iq))
           qsatgf=.622*esatf/(ps(iq)-esatf)
           prz = rho(iq)*cp*taftfh(iq)
           beta =      rmc(iq)/rmcmax(iq)
           devf= (hl*hlars/300.**2)*qsatgf*(1.-beta)/res(iq)
           dirad1 = 4.*rdg(iq)/300.
           print *,'dirad1,devf,prz ',dirad1,devf,prz
           print *,'beta,airr,res ',beta,airr(iq),res(iq)
           print *,'delta_tx ',delta_tx(iq)
           print *,'otgf,tgfnew,residf ',otgf(iq),tgfnew(iq),residf(iq)
         endif ! (ntest.eq.1)
       enddo   !  icount=1,5
      endif    ! (itnmeth.gt.0) 

!cdir nodep
      do ip=1,ipland  ! all land points in this nsib=3 loop
       iq=iperm(ip)
       if(tsigmf(iq) .le. .0101) then
         condxpr(iq)=condx(iq)
         evapfb(iq) = 0.
         evapxf(iq) = egg(iq)
         fgf(iq)  = fgg(iq)
         rdg(iq)=rgg(iq)
          tgf(iq) = tss(iq)
          Ewww(iq) = 0.
       else
         tgf(iq)=tgfnew(iq)
         wb(iq,1)=wb(iq,1)-evapfb1(iq)/(zse(1)*1000.)
         wb(iq,2)=wb(iq,2)-evapfb2(iq)/(zse(2)*1000.)
         wb(iq,3)=wb(iq,3)-evapfb3(iq)/(zse(3)*1000.)
         wb(iq,4)=wb(iq,4)-evapfb4(iq)/(zse(4)*1000.)
         wb(iq,5)=wb(iq,5)-evapfb5(iq)/(zse(5)*1000.)
         condxpr(iq)=(1.-tsigmf(iq))*condx(iq)+ tsigmf(iq)*condxg(iq)
         if(ntest.eq.1.and.abs(residf(iq)).gt.3.)
     .      print *,'iq,otgf(iq),tgf,delta_tx,residf '
     .              ,iq,otgf(iq),tgf(iq),delta_tx(iq),residf(iq)
       endif          ! tsigmf .le. .01   ... else ...
       fev(iq)=evapfb(iq)/dt*hl*tsigmf(iq) ! passed to soilsnow to update wb
       fes(iq)=(1.-tsigmf(iq))*egg(iq)*cls(iq)  ! also passed to soilsnow
       otgsoil(iq)=isflag(iq)*tggsn(iq,1) + (1-isflag(iq))*tgg(iq,1)
       if(ntest.eq.1.and.iq.eq.idjd)then
         isoil = isoilm(iq)
         iveg=ivegt(iq)
         print *,'in sib3 before soilsnowv'
         print *,'evapxf,epot,egg,fev,wetfac '
     .           ,evapxf(iq),epot(iq),egg(iq),fev(iq),wetfac(iq)
         print *,'fgf,fgg,fes ',fgf(iq),fgg(iq),fes(iq)
         print *,'isoil,ssat,tsigmf,rlai ',
     .            isoil,ssat(isoil),tsigmf(iq),rlai(iq)
         print *,'tgg1,t1,theta,tscrn '
     .           ,tgg(iq,1),t(iq,1),theta(iq),tscrn(iq)
         print *,'qg1,qsttg,Ewww '
     .           ,qg(iq,1),qsttg(iq),Ewww(iq)
         print *,'dfgdt,taftfhg,rho ',dfgdt(iq),taftfhg(iq),rho(iq)
         print *,'rmc,rmcmax(iq),qsatgf ',rmc(iq),rmcmax(iq),qsatgf
       endif
       if(ntest.ne.0)then
	  if(abs(tgf(iq)-otgf(iq)).gt.4.9)then
	    write(6,"('ktau,iq,otgf,tgf,dtgf,t1,t2',i4,i6,5f8.2)")
     .        ktau,iq,otgf(iq),tgf(iq),tgf(iq)-otgf(iq),t(iq,1),t(iq,2)
c          write(45,'(2g13.4)') sqrt(u(iq,1)**2+v(iq,1)**2),fg(iq)
	  endif
	endif
      enddo  !  ip=1,ipland
!-------------------------------------
c     print *,'before soilsnow'
      call soilsnowv
c     print *,'after soilsnow'
!cdir nodep
      do ip=1,ipland  ! all land points in this nsib=3 loop
       iq=iperm(ip)
c       iveg=ivegt(iq)
c       isoil = isoilm(iq)
c       if( tsigmf(iq).gt. .01)then
c         evapf=fev(iq)*dt/(hl*tsigmf(iq))
c         evapxf(iq) = (evapf/dt + Ewww(iq))*hl
c       endif  !  ( tsigmf(iq).gt. .01)
c       if(eg(iq).lt.0.)then
c         print *,'iq,eg,snowd,tgg ',iq,eg(iq),snowd(iq),tgg(iq,1)
c       endif

       if(isflag(iq).eq.0) then
        deltat=tgg(iq,1)-otgsoil(iq)
        fgg(iq)=fgg(iq)+deltat*dfgdt(iq)
        egg(iq)=egg(iq)+deltat*degdt(iq)
        egg(iq)=min(egg(iq),wb(iq,1)*zse(1)*1000.*hl/dt)
        rgg(iq)=rgg(iq)+deltat*dirad(iq)
       else
        deltat=tggsn(iq,1)-otgsoil(iq)
        fgg(iq)=fgg(iq)+deltat*dfgdt(iq)
        egg(iq)=egg(iq)+deltat*degdt(iq)
        rgg(iq)=rgg(iq)+deltat*dirad(iq)
       endif
c                                               combined fluxes
       if(snowd(iq).gt.1.)then
         eg(iq)=tsigmf(iq)*evapxf(iq) + egg(iq)
       else
         eg(iq) = tsigmf(iq)*evapxf(iq) + (1. - tsigmf(iq))*egg(iq)
       endif
       if(nsigmf.eq.2)then
         fg(iq)=tsigmf(iq)*fgf(iq)+fgg(iq)
       else
         fg(iq)=tsigmf(iq)*fgf(iq)+(1.-tsigmf(iq))*fgg(iq)
       endif
       rnet(iq)=-slwa(iq)-(1.-tsigmf(iq))*rgg(iq)-tsigmf(iq)*rdg(iq)
!      rnet(iq)=-slwa(iq)-(1.-tsigmf(iq))*rgg(iq)
!    .   -tsigmf(iq)*(extin(iq)*rgg(iq) + (1.-extin(iq))*rdg(iq)) ! 9/3/99

       if(ntest.eq.1.and.iq.eq.idjd)then
        print *,'even further down sib3 after soilsnowv'
        print *,'tgg ',(tgg(iq,k),k=1,ms)
        print *,'wb ',(wb(iq,k),k=1,ms)
        print *,'snowd ',snowd(iq)
        print *,'evapfb,fev,Ewww ',evapfb(iq),fev(iq),Ewww(iq)
        print *,'tsigmf,evapxf,egg ',tsigmf(iq),evapxf(iq),egg(iq)
        print *,'deltat,degdt,wb,zse ',
     .           tgg(iq,1)-otgsoil(iq),degdt(iq),wb(iq,1),zse(1)
        print *,'isflag,eg,fg ',isflag(iq),eg(iq),fg(iq)
       endif

        tgss=isflag(iq)*tggsn(iq,1) + (1-isflag(iq))*tgg(iq,1)  ! jlm
        if(tsigmf(iq).le. .01) then
          tss(iq) = tgss
          tgf(iq) = tgss
        else
          tss(iq)=tsigmf(iq)*tgf(iq)+(1.-tsigmf(iq))*tgss
        endif       ! tsigmf.le. .01

!eak!      following are just diagnostic arrays used by eak
!eak       if(ndiag_arr.eq.1)then
!eak       snost=ntau+1-1
!eak       precipst(iq)=precipst(iq)+condx(iq)/snost
!eak       evapst(iq)=evapst(iq)+evap(iq)/snost
!eak       runoffst(iq)=runoffst(iq)+runoff(iq)/snost
!eak       rnetst(iq)=rnetst(iq)+rnet(iq)/snost
!eak       egst(iq)=egst(iq)+eg(iq)/snost
!eak       fgst(iq)=fgst(iq)+fg(iq)/snost
!eak       tssst(iq)=tssst(iq)+tss(iq)/snost
!eak       wb1m(iq)=0.
!eak       do k=1,4
!eak        wbavst(iq)=wbavst(iq)+max(0.,(wb(iq,k)-swilt(isoil)))/snost
!eak        wb1m(iq)=wb1m(iq)+max(0.,(wb(iq,k)-swilt(isoil))*
!eak     &             zse(k)*1000.)
!eak       enddo  ! k=1,4
!eak       endif  ! (ndiag_arr.eq.1)

      enddo   ! ip=1,ipland

      return
      end
