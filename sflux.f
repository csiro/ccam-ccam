      subroutine sflux(nalpha)              ! for globpe code
      use diag_m
      use cc_mpi
      parameter (nblend=0)  ! 0 for original non-blended, 1 for blended af
      parameter (ntss_sh=0) ! 0 for original, 3 for **3, 4 for **4
      parameter (nplens=0)  ! 0 to turn off plens, 10 (e.g.) is on
!     parameter (lake=0)    ! 0 usual, 1 for specified lake points
!                             - replaced by nspecial in parm.h
      parameter (ntest=0)   ! ntest= 0 for diags off; ntest= 1 for diags on
!     parameter (ntaft=3)   ! 0 for original, 3 nowadays
!                   1 & 2 tafthf constrained by prior values with 2 faster
!                   3 uses measure of prior tgf in calc. fh
!     parameter (newztsea=1)   ! 0 for original, 1 for correct zt over sea
!     vmag introduced Mar '05 as vmod being used in ri
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
c     include 'map.h'      ! land
      include 'morepbl.h'  ! condx,fg,eg
      include 'nsibd.h'    ! rsmin,ivegt,sigmf,tgf,ssdn,res,rmc,tsigmf
      include 'parm.h'
      include 'pbl.h'
      include 'permsurf.h'
      include 'prec.h'     ! evap
      include 'savuvt.h'
      include 'screen.h'   ! tscrn,qgscrn,uscrn,rhscrn,u10
      include 'sigs.h'
      include 'soil.h'     ! ... zmin zolod zolog sicedep fracice alb
      include 'soilv.h'    ! ... ssat
      include 'soilsnow.h' ! soil arrays 
      include 'soilsnin.h'
      include 'tracers.h'  ! ngas, nllp, ntrac
      include 'trcom2.h'   ! nstn,slat,slon,istn,jstn
      include 'vegpar.h'   !
      include 'vvel.h'
      common/tafts/taftfh(ifull),taftfhg(ifull)
      common/work2/dirad(ifull),dfgdt(ifull),degdt(ifull)
     . ,wetfac(ifull),degdw(ifull),cie(ifull)
     . ,factch(ifull),qsttg(ifull),rho(ifull),zo(ifull)
     . ,aft(ifull),fh(ifull),ri(ifull),theta(ifull)
     . ,gamm(ifull),rg(ifull),vmod(ifull),dummwk2(ifull) ! rg in radriv90
      real taftfhg_temp(ifull)
!     following common block makes available other arrays for diag. output 
      common/work3/egg(ifull),evapxf(ifull),ewww(ifull),fgf(ifull),
     . fgg(ifull),ggflux(ifull),rdg(ifull),rgg(ifull),residf(ifull),
     . otgf(ifull),rmcmax(ifull),tgfnew(ifull),extin(ifull),af(ifull),
     . dum3(5*ijk-14*ifull)
      real plens(ifull),vmag(ifull),charnck(ifull)
      save plens
      data plens/ifull*0./
      integer list(20,2)
      save list
      include 'establ.h'

       list(1,1)=8290
       list(2,1)=9061
       list(3,1)=9093
       list(4,1)=9293
       list(5,1)=9540

c     stability dependent drag coefficients using Louis (1979,blm) f'
c     n.b. cduv, cdtq are returned as drag coeffs mult by vmod
c          (cduv=cduv*vmod; cdtq=cdtq*vmod)

c     t, u, v, qg are current values
c     tss is surface temperature
c     dw is soil wetness availability (1. for ocean) - not needed here
c     fg is sensible heat flux (was h0)
c     eg is latent heat flux (was wv)
c     dfgdt is dfgdt (was csen in surfupa/b)
c     degdt is degdt (was ceva in surfupa/b)

      ri_max=(1./fmroot -1.)/bprm  ! i.e. .14641
c     zobgin = .05   ! jlm: NB seems to be .01 in csiro9. Mar '05: in parm.h
      zologbgin=log(zmin/zobgin)   ! pre-calculated for all except snow points
      ztv=exp(vkar/sqrt(chn10)) /10.  ! proper inverse of ztsea
      z1onzt=300.*rdry*(1.-sig(1))*ztv /grav
      chnsea=(vkar/log(z1onzt))**2    ! should give .00085 for csiro9

      if(nspecial==1)then
         print*, "SFLUX, nspecial==1 doesn't work in MPI version"
         if(nspecial==1)stop
        do j=68,71  ! Eyre
         do i=27,28
          iq=i+(j-1)*il
          sigmf(iq)=0.  ! bare "soil"
          do kk=1,ms
           wb(iq,kk)=ssat(isoilm(iq))
          enddo
          if(ktau==1)then
            print *,'lake iq,isoilm,ivegt,sigmf,zolnd ',
     .                    iq,isoilm(iq),ivegt(iq),sigmf(iq),zolnd(iq)
          endif ! (ktau==1)
         enddo
        enddo
        do j=64,67  ! Torrens
         do i=28,28
          iq=i+(j-1)*il
          sigmf(iq)=0.  ! bare "soil"
          do kk=1,ms
           wb(iq,kk)=ssat(isoilm(iq))
          enddo
          if(ktau==1)then
            print *,'lake iq,isoilm,ivegt,sigmf,zolnd ',
     .                    iq,isoilm(iq),ivegt(iq),sigmf(iq),zolnd(iq)
          endif ! (ktau==1)
         enddo
        enddo
      endif   !  (nspecial==1)

      if(ktau==1)then
       taftfh(:)=.05        ! just a diag default for sea points
       taftfhg(:)=7.e-4     ! just a diag default for sea points
       if(nrungcm==3)then   ! for runs for PIRCS
        if (mydiag) print *,'entering sflux w_ms: ',wb(idjd,ms)
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
        if (mydiag) print *,'then w_ms: ',wb(idjd,ms)
       endif     !  (ktau==1).
      endif      !  (nrungcm==3)

      if (diag.or.ntest==1) then
        if (mydiag.and.land(idjd)) then
          print *,'entering sflux ktau,nsib,ivegt,isoilm,land '
     .         ,ktau,nsib,ivegt(idjd),isoilm(idjd),land(idjd)
          print *,'idjd,id,jd,slwa,sgsave ',
     .           idjd,id,jd,slwa(idjd),sgsave(idjd)
          print *,'snowd,sicedep,condx ',
     .           snowd(idjd),sicedep(idjd),condx(idjd)
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
       vmod(iq)=sqrt(uav**2+vav**2)  ! i.e. vmod for tss_sh
      enddo
      vmag(:)=max( vmod(:) , vmodmin)
      if(ntsur.ne.7)vmod(:)=vmag(:)	! gives usual way

      if(ntest==2.and.mydiag)print *,'before sea loop'
!      from June '03 use basic sea temp from tgg1 (so leads is sensible)      
!     all sea points in this loop; also open water of leads  
      if(charnock>0.)then
        charnck(:)=charnock
      elseif(charnock<-1.)then            ! zo like Moon (2004)
        if(ktau==1)u10(:)=vmod(:)/3.8
        charnck(:)=max(.0000386*u10(:),.000085*u10(:)-.00058)
      else
        if(ktau==1)u10(:)=vmod(:)/3.8
        charnck(:)=.008+3.e-4*(u10(:)-9.)**2/ ! like Makin (2002)
     &            (1.+(.006+.00008*u10(:))*u10(:)**2)
      endif     
      do iq=1,ifull
       if(.not.land(iq))then 
        wetfac(iq)=1.                                                ! sea
!       tgg2 (called tpan from Oct 05) holds effective skin sst for this loop 
        if(ntss_sh==0)then
         dtsol=.01*sgsave(iq)/(1.+.25*vmod(iq)**2)   ! solar heating ! sea
         tpan(iq)=tgg(iq,1)+tss_sh*min(dtsol,8.)     ! of ssts       ! sea
c        if(iq==idjd)print*,'ip,iq,tgg1,tss_sh,sgsave,vmod,dtsol ',
c     &               ip,iq,tgg(iq,1),tss_sh,sgsave(iq),vmod(iq),dtsol  
        elseif(ntss_sh==3)then                                       ! sea
         dtsol=tss_sh*.01*sgsave(iq)/                                ! sea
     .                (1.+.035*vmod(iq)**3)          ! solar heating ! sea
         tpan(iq)=tgg(iq,1)+min(dtsol,8.)            ! of ssts       ! sea
        elseif(ntss_sh==4)then                                       ! sea
         dtsol=tss_sh*.01*sgsave(iq)/                                ! sea
     .                (1.+vmod(iq)**4/81.)           ! solar heating ! sea
         tpan(iq)=tgg(iq,1)+min(dtsol,8.)            ! of ssts       ! sea
        endif   ! (ntss_sh==0) .. else ..
        if(nplens.ne.0)then
!        calculate running total (over last 24 h) of daily precip in mm  jlm
         plens(iq)=(1.-dt/86400.)*plens(iq)+condx(iq)  ! in mm/day
!        scale so that nplens m/s wind for 1/2 hr reduces effect by 1/1.2
!        plens(iq)=plens(iq)/(1.+vmod(iq)*dt*.2/(nplens*1800.))
         plens(iq)=plens(iq)/(1.+vmod(iq)*dt*.2/
     .                    max(nplens*1800.,1.))      ! avoids Cray compiler bug
!        produce a cooling of 4 K for an effective plens of 10 mm/day
         tpan(iq)=tpan(iq)-min(.4*plens(iq) , 6.)
        endif   !  (nplens.ne.0)
       if(ntsea==1.and.condx(iq)>.1)tpan(iq)=t(iq,2)  
       if(ntsea==2.and.condx(iq)>.1)tpan(iq)=t(iq,1)  
       if(ntsea==3.and.condx(iq)>.1)tpan(iq)=.5*(t(iq,2)+tgg(iq,1))  
       if(ntsea==4.and.condx(iq)>.1)tpan(iq)=.5*(t(iq,1)+tgg(iq,1)) 
       endif  ! (.not.land(iq)) 
      enddo   ! iq loop

!     here calculate fluxes for sea point, and nominal pan points	 
      afrootpan=vkar/log(zmin/panzo)                       
      do iq=1,ifull
c      drag coefficients  for momentum           cduv        
c      for heat and moisture  cdtq                  
       es = establ(tpan(iq))                                         ! sea
       constz=ps(iq)-es                                              ! sea
       qsttg(iq)= .98*.622*es/constz  ! with Zeng 1998 for sea water ! sea
       drst=qsttg(iq)*ps(iq)*hlars/(constz*tpan(iq)**2)              ! sea
       xx=grav*zmin*(1.-tpan(iq)*srcp/t(iq,1))                       ! sea
       ri(iq)=min(xx/vmag(iq)**2 , ri_max)                           ! sea  
!      if(ngas>0)stop 'call co2sflux'                                ! sea
c      this is in-line ocenzo using latest coefficient, i.e. .018    ! sea
       consea=vmod(iq)*charnck(iq)/grav  ! usually charnock=.018     ! sea
       if(land(iq))then
         zo(iq)=panzo
         af(iq)=afrootpan**2                                         ! sea
       else
        if(charnock<-1.)then  ! Moon (2004) over sea
         zo(iq)=charnck(iq)
         afroot=vkar/log(zmin/zo(iq))                                ! sea
         af(iq)=afroot**2                                            ! sea
        else       ! usual charnock method over sea
!         zo(iq)=.01                  
         zo(iq)=.001    ! .0005 better first guess            
         if(ri(iq)>0.)then             ! stable sea points           ! sea
           fm=vmod(iq) /(1.+bprm*ri(iq))**2 ! N.B. this is vmod*fm   ! sea
           con=consea*fm                                             ! sea
           do it=1,3                                                 ! sea
            afroot=vkar/log(zmin/zo(iq))                             ! sea
            af(iq)=afroot**2                                         ! sea
            daf=2.*af(iq)*afroot/(vkar*zo(iq))                       ! sea
            zo(iq)=max(1.5e-5,zo(iq)-(zo(iq)-con*af(iq))/(1.-con*daf))
           enddo    ! it=1,3                                         ! sea
           afroot=vkar/log(zmin/zo(iq))                              ! sea
           af(iq)=afroot**2                                          ! sea
         else                        ! unstable sea points           ! sea
           do it=1,3                                                 ! sea
            afroot=vkar/log(zmin/zo(iq))                             ! sea
            af(iq)=afroot**2                                         ! sea
            daf=2.*af(iq)*afroot/(vkar*zo(iq))                       ! sea
            con1=cms*2.*bprm*sqrt(-ri(iq)*zmin/zo(iq))               ! sea
            den=1.+af(iq)*con1                                       ! sea
            dden=con1*(daf-.5*af(iq)/zo(iq))                         ! sea
            fm=vmod(iq)-vmod(iq)*2.*bprm *ri(iq)/den                 ! sea
            dfm=2.*bprm*ri(iq)*dden/den**2                           ! sea
            zo(iq)=max(1.5e-5,zo(iq)-(zo(iq)-consea*af(iq)*fm)/      ! sea
     .                       (1.-consea*(daf*fm+af(iq)*dfm)))        ! sea
           enddo  ! it=1,3                                           ! sea
         endif    ! (xx>0.) .. else..                                ! sea
        endif     ! (charnock<-1.) .. else ..
       endif      ! (land(iq)) .. else ..
       aft(iq)=chnsea                                                ! sea
      enddo  ! iq loop   

      if(newztsea==0)then ! 0 for original, 1 for different zt over sea
c       enhanced formula used in Feb '92 Air-Sea conference follows: ! sea
c       factch=sqrt(zo*exp(vkar*vkar/(chnsea*log(zmin/zo)))/zmin)    ! sea
        factch(:)=1. ! factch is sqrt(zo/zt) only for use in unstable fh
      else
        factch(:)=sqrt(zo(:)*ztv) ! for use in unstable fh
      endif  ! (newztsea==0)

      do iq=1,ifull                                           
c      Having settled on zo & af now do actual fh and fm calcs       ! sea
       if(ri(iq)>0.)then                                             ! sea
         fm=vmod(iq)/(1.+bprm*ri(iq))**2  ! no zo contrib for stable ! sea
         fh(iq)=fm                                                   ! sea
       else        ! ri is -ve                                       ! sea
         root=sqrt(-ri(iq)*zmin/zo(iq))                              ! sea
c        First do momentum                                           ! sea
         denma=1.+cms*2.*bprm*af(iq)*root                            ! sea
         fm=vmod(iq)-vmod(iq)*2.*bprm *ri(iq)/denma                  ! sea
c        n.b. fm denotes ustar**2/(vmod(iq)*af)                      ! sea
c        Now heat ; allow for smaller zo via aft and factch          ! sea
c        N.B. for newztsea=1, zo contrib cancels in factch*root,
c        so eg (& epan) and fg  (also aft) then indept of zo
         denha=1.+chs*2.*bprm*factch(iq)*aft(iq)*root                ! sea
         fh(iq)=vmod(iq)-vmod(iq)*2.*bprm *ri(iq)/denha              ! sea
       endif                                                         ! sea
                                                                     ! sea
       conh=rho(iq)*aft(iq)*cp                                       ! sea
       conw=rho(iq)*aft(iq)*hl                                       ! sea
       fg(iq)=conh*fh(iq)*(tpan(iq)-theta(iq))                       ! sea
       eg(iq)=conw*fh(iq)*(qsttg(iq)-qg(iq,1))                       ! sea
c      cduv is now drag coeff *vmod                                  ! sea
       cduv(iq) =af(iq)*fm                                           ! sea
       ustar(iq) = sqrt(vmod(iq)*cduv(iq))                           ! sea
c      Surface stresses taux, tauy: diagnostic only - unstaggered now 
       taux(iq)=rho(iq)*cduv(iq)*u(iq,1)                             ! sea
       tauy(iq)=rho(iq)*cduv(iq)*v(iq,1)                             ! sea
       ! note that iq==idjd  can only be true on the correct processor

       if(ntest==1.and.iq==idjd)then                                 ! sea
         print *,'in sea-type loop for iq,idjd: ',                   ! sea
     .                            iq,idjd,ip                         ! sea
         print *,'zmin,zo,factch ',zmin,zo(iq),factch(iq)            ! sea         
         print *,'ri,ustar,es ',ri(iq),ustar(iq),es                  ! sea
         print *,'af,aft,tgg2 ',af(iq),aft(iq),tpan(iq)              ! sea
         print *,'tgg2,tss,theta ',tpan(iq),tss(iq),theta(iq)        ! sea
         print *,'chnsea,rho,t1 ',chnsea,rho(iq),t(iq,1)             ! sea
         print *,'fm,fh,conh ',fm,fh(iq),conh                        ! sea
         print *,'vmod,cduv,fg ',vmod(iq),cduv(iq),fg(iq)            ! sea
       endif                                                         ! sea
      enddo     ! iq loop                                            ! sea
      epot(:) = eg(:)                                               
      epan(:) = eg(:)                         
c     section to update pan temperatures
      do iq=1,ifull
       if(land(iq))then
         rgg(iq)=5.67e-8*tpan(iq)**4    
!        assume gflux = 0
!        note pan depth=.254 m, spec heat water=4186 joule/kg K
!        and change in heat supplied=spec_heatxmassxdelta_T
         ga(iq)=-slwa(iq)-rgg(iq)-panfg*fg(iq)
         tpan(iq)=tpan(iq)+ga(iq)*dt/(4186.*.254*1000.)             
       endif  ! (land(iq))
      enddo   ! iq loop

      if(nmaxpr==1.and.mydiag)then
        iq=idjd
        write (6,"('after sea loop fg,tpan,epan,ri,fh,vmod',
     &    9f9.4)") fg(idjd),tpan(idjd),epan(idjd),
     &             ri(idjd),fh(idjd),vmod(idjd)
        write (6,"('u10,ustar,charnck,zo,cd',
     &    3f9.4,2f9.6)") u10(idjd),ustar(idjd),charnck(idjd),zo(idjd),
     &    cduv(idjd)/vmod(idjd)
!        if(ri(iq)<0.)then
!          root=sqrt(-ri(iq)*zmin/zo(iq)) 
!          denma=1.+cms*2.*bprm*af(iq)*root                            
!          fm=vmod(iq)-vmod(iq)*2.*bprm *ri(iq)/denma                  
!          denha=1.+chs*2.*bprm*factch(iq)*aft(iq)*root  
!        endif             
        if(ri(iq)>0.)then
          fm=vmod(iq)/(1.+bprm*ri(iq))**2
        else
          root=sqrt(-ri(iq)*zmin/zo(iq))
          denma=1.+cms*2.*bprm*af(iq)*root
          fm=vmod(iq)-vmod(iq)*2.*bprm *ri(iq)/denma
          denha=1.+chs*2.*bprm*factch(iq)*aft(iq)*root
        endif

        write (6,"('after sea loop af,aft,factch,root,denha,denma,fm',
     &    9f9.4)") af(iq),aft(iq),factch(iq),root,denha,denma,fm
      endif                       
      zminlog=log(zmin)

      do iq=1,ifull
       if(sicedep(iq)>0.)then
!      non-leads for sea ice points                                 ! sice
!      N.B. tgg( ,3) holds tice                                     ! sice
c       iq=iperm(ip)                                                ! sice
       es = establ(tgg(iq,3))                                       ! sice
       constz=ps(iq)-es                                             ! sice
       qsttg(iq)= .622*es/constz                                    ! sice
       drst=qsttg(iq)*ps(iq)*hlars/(tgg(iq,3)*tgg(iq,3)*constz)     ! sice
       xx=grav*zmin*(1.-tgg(iq,3)*srcp/t(iq,1))                     ! sice
       ri_ice=min(xx/vmag(iq)**2 , ri_max)                          ! sice
       factch(iq)=sqrt(7.4)  ! same as land from 27/4/99            ! sice
!      factch(iq)=1.   ! factch is sqrt(zo/zt) for use in unstable fh  
       zoice=.001                                                   ! sice
       zologice=zminlog-log(zoice)   !   i.e. log(zmin/zo(iq))      ! sice
       af(iq)=(vkar/zologice)**2                                    ! sice
       aft(iq)=vkar**2/(zologice*(2.+zologice) )  ! from 27/4/99    ! sice
!      aft(iq)=af                                 ! up till 27/4/99 ! sice
       wetfac(iq)=1+.008*(tgg(iq,3)-273.16) ! 008*tgg(iq,3)-1.18528 ! sice
                                                                    ! sice
c      now do fh and fm calcs for sice                              ! sice
       if(ri_ice>0.)then                                            ! sice
         fm=vmod(iq)/(1.+bprm*ri_ice)**2                            ! sice
         fh(iq)=fm                                                  ! sice
       else                                                         ! sice
         root=sqrt(-ri_ice*zmin/zoice)                              ! sice
c        First do momentum                                          ! sice
         denma=1.+cms*2.*bprm*af(iq)*root                           ! sice
         fm=vmod(iq)-vmod(iq)*2.*bprm *ri_ice/denma                 ! sice
c        n.b. fm denotes ustar**2/(vmod(iq)*af)                     ! sice
c        Now heat ; allow for smaller zo via aft and factch         ! sice
         denha=1.+chs*2.*bprm*factch(iq)*aft(iq)*root               ! sice
         fh(iq)=vmod(iq)-(2.*bprm *ri_ice)/denha                    ! sice
       endif                                                        ! sice
                                                                    ! sice
       conh=rho(iq)*aft(iq)*cp                                      ! sice
       conw=rho(iq)*aft(iq)*hl                                      ! sice
!      fgice & egice renamed as fgf and fev from Aug 05 to aid diags	
       fgf(iq)=conh*fh(iq)*(tgg(iq,3)-theta(iq))                    ! sice
       dfgdt(iq)=conh*fh(iq)                                        ! sice
       if(ntest==1.and.iq==idjd)then                                ! sice
         print *,'in sice loop'                                     ! sice
         print *,'zmin,zo,wetfac ',zmin,zoice,wetfac(iq)            ! sice
         print *,'ri_ice,es ',ri_ice,es                             ! sice
         print *,'af,aft,ustar ',af(iq),aft(iq),ustar(iq)           ! sice
         print *,'chnsea,rho ',chnsea,rho(iq)                       ! sice
         print *,'fm,fh,conh ',fm,fh(iq),conh                       ! sice
       endif                                                        ! sice
                                                                    ! sice
       if(nalpha==1)then    ! beta scheme         sice here         ! sice
         epotice=conw*fh(iq)*(qsttg(iq)-qg(iq,1))                   ! sice
         fev(iq)=wetfac(iq)*epotice                                 ! sice
         degdt(iq)=wetfac(iq)*conw*fh(iq)*drst                      ! sice
       else                   ! alpha scheme                        ! sice
c        following trick reduces -ve evap (dew) to 1/10th value     ! sice
         qtgnet=qsttg(iq)*wetfac(iq) -qg(iq,1)                      ! sice
         qtgair=qsttg(iq)*wetfac(iq)-max(qtgnet,.1*qtgnet)          ! sice
         eg2=-conw*fh(iq)*qtgair                                    ! sice
         eg1=conw*fh(iq)*qsttg(iq)                                  ! sice
         fev(iq) =eg1*wetfac(iq) +eg2                               ! sice
         epotice    = conw*fh(iq)*(qsttg(iq)-qg(iq,1))              ! sice
         deg=wetfac(iq)*conw*fh(iq)*drst                            ! sice
c        following reduces degdt by factor of 10 for dew            ! sice
         degdt(iq)=.55*deg+sign(.45*deg,qtgnet)                     ! sice
       endif                                                        ! sice
                                                                    ! sice
c      section to update sea ice surface temperature;               ! sice
c      specified sea-ice thickness                                  ! sice
c      over sea ice, set a minimum depth for this experiment of .1  ! sice
!      sicedep(iq) = max(sicedep(iq) , 0.1)  fixed in indata/nestin from Jan 06
c      no snow on the ice assumed for now                           ! sice
       gamm(iq) = 3.471e+05                                         ! sice
       cie(iq) = 2.04/sicedep(iq)                                   ! sice
       rgg(iq)=5.67e-8*tgg(iq,3)**4                                 ! sice
!      gflux here is	flux from ice to water, +ve downwards          ! sice
       gflux(iq)=cie(iq)*(tgg(iq,3)-271.2)                          ! sice
       ga(iq)=-slwa(iq)-rgg(iq)-fev(iq)-fgf(iq)-gflux(iq)           ! sice
       dirad(iq)=4.*5.67e-8*tgg(iq,3)**3                            ! sice
       b1=dirad(iq)+degdt(iq)+dfgdt(iq)+cie(iq)                     ! sice
       gbot=(gamm(iq)/dt)+b1                                        ! sice
       deltat=ga(iq)/gbot                                           ! sice
       tgg(iq,3)=tgg(iq,3)+deltat                                   ! sice
       tgg(iq,3)=min(tgg(iq,3),271.2)   ! jlm fix Tue  05-30-2000
       fgf(iq) =fgf(iq) +deltat*dfgdt(iq)                           ! sice
       fev(iq) =fev(iq) +deltat*degdt(iq)                           ! sice
       es = establ(tgg(iq,3))                                       ! sice
       constz=ps(iq)-es                                             ! sice
       qsttg(iq)=.622*es/constz                                     ! sice
                                                                    ! sice
!      combine ice and leads contributions here                     ! sice
       eg(iq) =fracice(iq)*fev(iq) + (1.-fracice(iq))*eg(iq)        ! sice
       fg(iq) = fracice(iq)*fgf(iq)+ (1.-fracice(iq))*fg(iq)        ! sice
       ri(iq) =fracice(iq)*ri_ice + (1.-fracice(iq))*ri(iq)  ! for scrnout
       zo(iq) =fracice(iq)*zoice  + (1.-fracice(iq))*zo(iq)  ! for scrnout
       cduv(iq) =fracice(iq)*af(iq)*fm + (1.-fracice(iq))*cduv(iq)  ! sice
       ustar(iq) = sqrt(vmod(iq)*cduv(iq))                          ! sice
c      N.B. potential evaporation is now eg+eg2                     ! sice
       epot(iq) =fracice(iq)*epotice + (1.-fracice(iq))*epot(iq)    ! sice
       tss(iq) = fracice(iq)*tgg(iq,3)+(1.-fracice(iq))*tpan(iq)    ! 2004
c      Surface stresses taux, tauy: diagnostic only - unstag now    ! sice
       taux(iq)=rho(iq)*cduv(iq)*u(iq,1)                            ! sice
       tauy(iq)=rho(iq)*cduv(iq)*v(iq,1)                            ! sice
      endif  ! (sicedep(iq)>0.)
      enddo       ! iq loop                           
c     if(mydiag.and.diag)then                                           
      if(mydiag.and.nmaxpr==1)then 
         print *,'after sice loop'
         iq=idjd
         b1=dirad(iq)+degdt(iq)+dfgdt(iq)+cie(iq)   
         gbot=(gamm(iq)/dt)+b1                           
         deltat=ga(iq)/gbot                                               
         print *,'ri,vmag,vmod,cduv ',ri(iq),vmag(iq),vmod(iq),cduv(iq) 
         print *,'fh,tss,tpan,tgg3 ',fh(iq),tss(iq),tpan(iq),tgg(iq,3) 
         print *,'theta,t1,deltat ',theta(iq),t(iq,1),deltat           
         print *,'b1,ga,gbot,af,aft ',b1,ga(iq),gbot,af(iq),aft(iq) 
         print *,'fg,fgice,factch ',fg(iq),fgf(iq),factch(iq) 
         print *,'cie ',cie(iq)      
         print *,'eg,egice(fev),ustar ',eg(iq),fev(iq),ustar(iq)          
      endif   ! (diag)                                        
c----------------------------------------------------------------------
!cdir nodep
      do ip=1,ipland  ! all land points in this shared loop         ! land
c      fh itself was only used outside this loop in sib0 (jlm)      ! land
       iq=iperm(ip)                                                 ! land
       zobg=zobgin                                                  ! land
       es = establ(tss(iq))                                         ! land
       qsttg(iq)= .622*es/(ps(iq)-es)  ! prim for scrnout, bur recalc end sib3    
c      factch is sqrt(zo/zt) for land use in unstable fh            ! land
       factch(iq)=sqrt(7.4)                                         ! land
       if(snowd(iq)>0.)then                                         ! land
!        reduce zo over snow;
         zobg=max(zobgin -snowd(iq)*0.00976/12., 0.00024)           ! land
         zologbg=log(zmin/zobg)                                     ! land
!        following line is bit simpler than csiro9                  ! land
         zo(iq)=max(zolnd(iq) -.001*snowd(iq), .01)                 ! land
         zologx=log(zmin/zo(iq))                                    ! land
       else  ! land but not snow                                    ! land
         zo(iq)=zolnd(iq)                                           ! land
         zologbg=zologbgin                                          ! land
         zologx=zolog(iq)
       endif     ! (snowd(iq)>0.)                                   ! land
       if(nblend==1)then  ! blended zo for momentum                 ! land
!        note that Dorman & Sellers zo is already an average, 
!        accounting for sigmf, so may not wish to further blend zo	
         afland=(vkar/((1.-sigmf(iq))*zologbg+sigmf(iq)*zologx))**2 ! land
       else    ! non-blended zo for momentum                        ! land
         afland=(vkar/zologx)**2                                    ! land
       endif   ! (nblend==1)                                        ! land
       aftland=vkar**2/( zologx * (2.+zologx) )                     ! land
       aft(iq)=aftland                                              ! land
       aftlandg=vkar**2/( zologbg * (2.+zologbg) )                  ! land
c      lgwd>0 enhances cduv (momentum) over orog under (stable & unst) condns
       if(lgwd>0)then
         af(iq)=afland+helo(iq)       ! jlm special gwd4b           ! land
       else
         af(iq)=afland                                              ! land
       endif
                                                                    ! land
c      Having settled on zo (and thus af) now do actual fh and fm calcs 
       xx=grav*zmin*(1.-tss(iq)*srcp/t(iq,1))                       ! land
       ri(iq)=min(xx/vmag(iq)**2 , ri_max)                          ! land
       if(ri(iq)>0.)then                                            ! land
         fm=vmod(iq)/(1.+bprm*ri(iq))**2                            ! land
         fh(iq)=fm                                                  ! land
         fhbg=fh(iq)                                                ! land
       else                                                         ! land
         root=sqrt(-ri(iq)*zmin/zo(iq))  ! ignoring blending here   ! land
c        First do momentum                                          ! land
         denma=1.+cms*2.*bprm*af(iq)*root                           ! land
         fm=vmod(iq)-vmod(iq)*2.*bprm *ri(iq)/denma                 ! land
c        n.b. fm denotes ustar**2/(vmod(iq)*af)                     ! land
c        Now heat ; allow for smaller zo via aft and factch         ! land
         denha=1.+chs*2.*bprm*factch(iq)*aft(iq)*root               ! land
         fh(iq)=vmod(iq)-vmod(iq)*2.*bprm *ri(iq)/denha             ! land
         rootbg=sqrt(-ri(iq)*zmin/zobg)                             ! land
         denhabg=1.+chs*2.*bprm*factch(iq)*aftlandg*rootbg          ! land
         fhbg=vmod(iq)-vmod(iq)*2.*bprm *ri(iq)/denhabg             ! land
       endif                                                        ! land
       taftfhg_temp(iq)=aftlandg*fhbg ! uses fmroot above, for sib3 ! land
       taftfhg(iq)=aftlandg*fhbg ! value used for ntaft=3 (may need improving)
                                                                    ! land
c      cduv is now drag coeff *vmod                                 ! land
       cduv(iq) =af(iq)*fm                                          ! land
       ustar(iq) = sqrt(vmod(iq)*cduv(iq))                          ! land
c      cdtq(iq) =aft(iq)*fh(iq)                                     ! land
c      Surface stresses taux, tauy: diagnostic only - unstaggered now   
       taux(iq)=rho(iq)*cduv(iq)*u(iq,1)                            ! land
       tauy(iq)=rho(iq)*cduv(iq)*v(iq,1)                            ! land

       if(ntest==1.and.iq==idjd)then                                ! land
         print *,'in main land loop'                                ! land
         print *,'zmin,zobg,zobgin,snowd ',zmin,zobg,zobgin,snowd(iq) 
         print *,'afland,aftland,zologbg ',afland,aftland,zologbg   ! land
         print *,'af,vmag,vmod,es ',af(iq),vmag(iq),vmod(iq),es     ! land
         print *,'tss,theta,t1 ',tss(iq),theta(iq),t(iq,1)          ! land
         print *,'aft,fm,fh,rho,conh ',aft(iq),fm,fh(iq),rho(iq),conh 
         print *,'ri,vmod,cduv,fg ',ri(iq),vmod(iq),cduv(iq),fg(iq) ! land
       endif  ! (ntest==1.and.iq==idjd)                             ! land
      enddo     ! ip=1,ipland                                       ! land
c     iq=899
c     print *,'in sflux; iq,af,zo,zolnd,snowd,zolog ',iq,af(iq),
c    &               zo(iq),zolnd(iq),snowd(iq),zolog(iq)  

      if(ntaft==0.or.ktau==1)then
        do iq=1,ifull  ! will only use land values
         if(land(iq))then
           taftfh(iq)=aft(iq)*fh(iq) ! uses fmroot above            ! land
           taftfhg(iq)=taftfhg_temp(iq)
         endif
        enddo
      elseif(ntaft==1.)then
        do iq=1,ifull         ! will only use land values
         if(land(iq))then
           thnew=aft(iq)*fh(iq) ! uses fmroot above                 ! land
           thgnew=taftfhg_temp(iq)
           if(thnew>2.*taftfh(iq).or.thnew<.5*taftfh(iq))then
             taftfh(iq)=.5*(thnew+taftfh(iq))
           else
             taftfh(iq)=thnew
           endif
           if(thgnew>2.*taftfhg(iq).or.thgnew<.5*taftfhg(iq))then
             taftfhg(iq)=.5*(thgnew+taftfhg(iq))
           else
             taftfhg(iq)=thgnew
           endif
         endif
        enddo
      elseif(ntaft==2)then    ! preferred faster option
        do iq=1,ifull         ! will only use land values
         if(land(iq))then
           thnew=aft(iq)*fh(iq) ! uses fmroot above                 ! land
           thgnew=taftfhg_temp(iq)
           thnewa=min(thnew,
     .            max(2.*taftfh(iq),.5*(thnew+taftfh(iq))))
           taftfh(iq)=max(thnewa,
     .                min(.5*taftfh(iq),.5*(thnew+taftfh(iq))))
           thgnewa=min(thgnew,
     .             max(2.*taftfhg(iq),.5*(thgnew+taftfhg(iq))))
           taftfhg(iq)=max(thgnewa,
     .                 min(.5*taftfhg(iq),.5*(thgnew+taftfhg(iq))))
         endif
        enddo
      elseif(ntaft==3)then
!       do vegetation calulation for taftfh	
        do iq=1,ifull
         if(land(iq))then
           xx=grav*zmin*(1.-tgf(iq)*srcp/t(iq,1)) ! actually otgf   ! land
           ri_tmp=min(xx/vmag(iq)**2 , ri_max)                      ! land
           if(ri_tmp>0.)then                                        ! land 
             fh_tmp=vmod(iq)/(1.+bprm*ri_tmp)**2                    ! land
           else                                                     ! land
             root=sqrt(-ri_tmp*zmin/zo(iq))  ! ignoring blending    ! land
c            Now heat ; allow for smaller zo via aft and factch     ! land
             denha=1.+chs*2.*bprm*factch(iq)*aft(iq)*root           ! land
             fh_tmp=vmod(iq)-vmod(iq)*2.*bprm *ri_tmp/denha         ! land
           endif
           taftfh(iq)=aft(iq)*fh_tmp ! uses fmroot above, for sib3  ! land 
         endif
        enddo
      endif  ! (ntaft==0.or.ktau==1)  .. else ..
      if(ntest>0)then
        print *,'before sib3 zo,zolnd,af ',zo(idjd),zolnd(idjd),af(idjd)
      endif
c ----------------------------------------------------------------------

      if(nsib==1.or.nsib==3)then
        if(ktau==1)then  !To initialize new nsib=1/3 run (or restart) only
          print *,'ipland,ipsice,ipsea in sflux: ',
     .             ipland,ipsice,ipsea
          do iq=1,ifull  ! give default over sea too
           tgf(iq) = t(iq,1)  ! was tss(iq)
           rmc(iq) = 0.
          enddo    
        endif  ! if(ktau==1)
        if(nsib==3)call sib3(nalpha)
      endif     !  (nsib==1 or 3)
      if(nsib==4)call sib4

c ----------------------------------------------------------------------

c     end of calls to sib1,2,3
      evap(:)=evap(:)+dt*eg(:)/hl !time integ value in mm (wrong for snow)

      if(diag.or.ntest>0)then
        if (mydiag) print *,'before call scrnout'
        call maxmin(t,' t',ktau,1.,kl)
      endif

      if(ntsur.ne.5)then    ! ntsur=6 is default from Mar '05  
c       preferred option to recalc cduv, ustar (gives better uscrn, u10)
        do iq=1,ifull
         afroot=vkar/log(zmin/zo(iq))! land formula is bit different above                             
         af(iq)=afroot**2                                             
         xx=grav*zmin*(1.-tss(iq)*srcp/t(iq,1))                       
         ri(iq)=min(xx/vmag(iq)**2 , ri_max)                            
         if(ri(iq)>0.)then                                            
           fm=vmod(iq)/(1.+bprm*ri(iq))**2         ! Fm * vmod
         else                                                           
           root=sqrt(-ri(iq)*zmin/zo(iq))  
           denma=1.+cms*2.*bprm*af(iq)*root                        
           fm=vmod(iq)-vmod(iq)*2.*bprm *ri(iq)/denma     ! Fm * vmod                        
c          n.b. fm denotes ustar**2/(vmod(iq)*af)                         
         endif                                                    
c        cduv is now drag coeff *vmod                                     
         cduv(iq) =af(iq)*fm                       ! Cd * vmod                                
         ustar(iq) = sqrt(vmod(iq)*cduv(iq))                            
c        Surface stresses taux, tauy: diagnostic only - unstaggered now   
         taux(iq)=rho(iq)*cduv(iq)*u(iq,1)                              
         tauy(iq)=rho(iq)*cduv(iq)*v(iq,1)                              
        enddo     
       endif  ! (ntsur==6)
       
       if(nproc==1.and.diag)then
          taftfhmin=1.e20
          taftfhmax=-100.
          taftfhgmin=1.e20
          taftfhgmax=-100.
          do iq=1,ifull
           if(taftfh(iq)<taftfhmin)then
             taftfhmin=taftfh(iq)           ! ~.0012
             iqmin1=iq
           endif
           if(taftfh(iq)>taftfhmax)then
             taftfhmax=taftfh(iq)           ! ~.13
             iqmax1=iq
           endif
           if(taftfhg(iq)<taftfhgmin)then
             taftfhgmin=taftfhg(iq)         ! ~.0006
             iqmin2=iq
           endif
           if(taftfhg(iq)>taftfhgmax)then
             taftfhgmax=taftfhg(iq)         ! ~.004
             iqmax2=iq
           endif
          enddo
          print *,'taftfhmin,taftfhmax ',
     &             taftfhmin,iqmin1,taftfhmax,iqmax1
          print *,'taftfhgmin,taftfhgmax ',
     &             taftfhgmin,iqmin2,taftfhgmax,iqmax2
       endif  ! (nproc==1.and.diag)

!     always call scrnout from 19/9/02
      call scrnout(zo,ustar,factch,wetfac,qsttg,            ! arrays
     .       qgscrn,tscrn,uscrn,u10,rhscrn,af,aft,ri,vmod,  ! arrays
     .       bprm,cms,chs,chnsea,nalpha)

c***  end of surface updating loop

      if(diag.or.ntest==1)then
         if ( mydiag ) then
            print *,'at end of sflux, after call scrnout'
            print *,'slwa,rdg,eg,fg ',
     .               slwa(idjd),rdg(idjd),eg(idjd),fg(idjd)
            print *,'tscrn,cduv,zolnd ',
     .               tscrn(idjd),cduv(idjd),zolnd(idjd)
            print *,'degdt,dfgdt,snowd,sicedep ',
     .               degdt(idjd),dfgdt(idjd),snowd(idjd),sicedep(idjd)
            print *,'u1,v1,qg1 ',u(idjd,1),v(idjd,1),qg(idjd,1)
            print *,'w,w2,condx ',
     .               wb(idjd,1),wb(idjd,ms),condx(idjd)
            print *,'t1,tss,tgg_2,tgg_ms ',
     .               t(idjd,1),tss(idjd),tgg(idjd,2),tgg(idjd,ms)
         end if
        call maxmin(tscrn,'tc',ktau,1.,1)
       endif
      if(ntest==4.and.ktau==10)then
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
!                                  2 for ewww diags      
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
c     include 'map.h'
      include 'morepbl.h'
      include 'nsibd.h'    ! rsmin,ivegt,sigmf,tgf,ssdn,res,rmc,tsigmf
      include 'parm.h'
      include 'pbl.h'
      include 'permsurf.h'
      include 'screen.h'   ! tscrn etc
      include 'sigs.h'
      include 'soil.h'     ! ... zmin zolnd zolog sice alb
      include 'soilv.h'    ! ssat, clay,..
      include 'soilsnin.h'
      include 'soilsnow.h' ! new soil arrays for scam - tgg too
      include 'tracers.h'  ! ngas, nllp, ntrac
      include 'trcom2.h'   ! nstn,slat,slon,istn,jstn
      include 'vegpar.h'   !

      common/nsib/nbarewet,nsigmf
      common/tafts/taftfh(ifull),taftfhg(ifull)
      common/work2/dirad(ifull),dfgdt(ifull),degdt(ifull)
     . ,wetfac(ifull),degdw(ifull),cie(ifull)
     . ,factch(ifull),qsttg(ifull),rho(ifull),zo(ifull)
     . ,aft(ifull),fh(ifull),ri(ifull),theta(ifull)
     . ,gamm(ifull),rg(ifull),vmod(ifull),dummwk2(ifull)
      common/work3/egg(ifull),evapxf(ifull),ewww(ifull),fgf(ifull),
     . fgg(ifull),ggflux(ifull),rdg(ifull),rgg(ifull),residf(ifull),
     . otgf(ifull),rmcmax(ifull),tgfnew(ifull),extin(ifull),af(ifull),
     . dum3(5*ijk-14*ifull)
      common/work3c/airr(ifull),cc(ifull),condxg(ifull),delta_tx(ifull),
     . evapfb1(ifull),evapfb2(ifull),evapfb3(ifull),evapfb4(ifull),
     . evapfb5(ifull),evapfb6(ifull),evapfb1a(ifull),evapfb2a(ifull),
     . evapfb3a(ifull),evapfb4a(ifull),evapfb5a(ifull),evapfb6a(ifull),
     . evapfb(ijk-16*ifull)  ! to allow > 18 levels
      common/work3d/dqsttg(ifull),tstom(ifull),dumxx(ifull),
     .   cls(ifull),omc(ifull),dum3d(ijk-5*ifull)  ! allows L9
      include 'establ.h'
      real xx1(ifull),xx2(ifull),xx3(ifull),zz1(ifull),swnet(ifull),
     .     lwnet(ifull)
      integer list(20,2),imonth(12)
      data imonth /31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31./
      dimension ndoy(12)   ! days from beginning of year (1st Jan is 0)
      data ndoy/ 0,31,59,90,120,151,181,212,243,273,304,334/

 

       ktauplus=0
       jyear=kdate/10000
       jmonth=(kdate-jyear*10000)/100
       do k=1,jmonth-1
        ktauplus = ktauplus + imonth(k)*nperday
       enddo

!     fle(isoil,w)=(w-swilt(isoil))/(sfc(isoil)-swilt(isoil))           !0 Eva's
!     fle(isoil,w)= w/ssat(isoil)                                       !1 simplest bare
!     fle(isoil,w)= (w-frac*max( w,swilt(isoil) ))/                     !2 jlm sugg. bare
!    .               (ssat(isoil)*(1.-frac))                            !2 jlm sugg. bare
!     fle(isoil,w)=10.*((w-ssoil(isoil))/ssat(isoil)+.1))               ! an old one
!     fle(isoil,w)= w/sfc(isoil)                                        ! jlm for PIRCS
!     fle(isoil,w)=(w-frac*swilt(isoil))/(sfc(isoil)-frac*swilt(isoil)) ! jlm special

!***  N.B. nrungcm=2 fix-up for Mark 2 wb moved to indata 6/11/00
       list(1,1)=8290
       list(2,1)=9061
       list(3,1)=9093
       list(4,1)=9293
       list(5,1)=9540

      if(ktau==1.and.nrungcm==3)then  ! for nsib runs for PIRCS
!cdir nodep
        do ip=1,ipland  ! all land points in this nsib=1/3 loop
         iq=iperm(ip)
c        code where eta provided (e.g. for PIRCS)
         do layer=2,ms
          wb(iq,layer)=wb(iq,1)   ! w, w2 and wb all same initially for PIRCS
         enddo       !  layer=1,ms
        enddo        !  ip=1,ipland
      endif          !  (ktau==1.and.nrungcm==3)

      if(ktau==1)then
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

        if ( mydiag.and.land(idjd) ) then
           iveg=ivegt(idjd)
           isoil = isoilm(idjd)
           tsoil=0.5*(0.3333*tgg(idjd,2)+0.6667*tgg(idjd,3)
     &             +0.95*tgg(idjd,4) +  0.05*tgg(idjd,5))
           ftsoil=max( 0. , 1.-.0016*(298.-tsoil)*max(0.,298.-tsoil) )
c          tsigmf(idjd)=max(0.,sigmf(idjd)-scveg44(iveg)*(1.-ftsoil))
           print *,'nsib,iveg,isoil,nalpha,newfgf,nsigmf,tsigmf ',
     &           nsib,iveg,isoil,nalpha,newfgf,nsigmf,
     &           max(0.,sigmf(idjd)-scveg44(iveg)*(1.-ftsoil))
           print *,'tsoil,ftsoil,scveg44,sigmf ',
     &           tsoil,ftsoil,scveg44(iveg),sigmf(idjd)
           print *,'swilt,sfc,wb1-6 ',
     &           swilt(isoil),sfc(isoil),(wb(idjd,k),k=1,ms)
           rlai_d= max(.1,rlaim44(iveg)-slveg44(iveg)*(1.-ftsoil))
           srlai=rlai_d+rlais44(iveg)
c          rsmin(idjd)=rsunc44(iveg)/rlai_d ! now always done
           print *,'rlai,srlai,rsmin ',rlai_d,srlai,
     &              rsunc44(iveg)/rlai_d 
        endif
      endif           ! (ktau==1)

c      due to a bug in the SX5 compiler even for cvsafe, 
c      following loop needed to be split off from the one after it
!cdir nodep
      do ip=1,ipland  ! all land points in this nsib=3 loop
       iq=iperm(ip)
       iveg=ivegt(iq)
c      evaluate seasonal impact and the snow depth
c      impact on the fractional vegetation cover
       tstom(iq)=298.
       if(iveg==6+31)tstom(iq)=302.
       if(iveg>=10.and.iveg<=21.and.
     .    abs(rlatt(iq)*180./pi)<25.)tstom(iq)=302.
       if(ntest==1.and.iq==idjd) then
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
c                            if( tsoil >= tstom ) ftsoil=1.
       rlai(iq)=  max(.1,rlaim44(iveg)-slveg44(iveg)*(1.-ftsoil))
       tsigmf(iq)=max(.001,  sigmf(iq)-scveg44(iveg)*(1.-ftsoil)) ! 13/5/03
       if(ntest==1.and.iq==idjd) then
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
       if(ntest==1.and.iq==idjd) then
         print *,'in sib3c ip,iq,idjd,iveg ',ip,iq,idjd,ivegt(iq)
         print*,'iveg,sigmf(iq),tsigmf ',ivegt(iq),sigmf(iq),tsigmf(iq)
         print*,'scveg44,snowd,zo,zolnd,tstom ',
     .          scveg44(iveg),snowd(iq),zo(iq),zolnd(iq),tstom(iq)
         print *,'alb,sgsave ',alb(iq),sgsave(iq)
         print*,'w2,rlai,extin ',wb(iq,ms),rlai(iq),extin(iq)
       endif ! ntest
c      bare ground calculation
!      if(isflag(iq)==1)then
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

      if(nbarewet==0)then  ! original Eva's, same as NCAR
!cdir nodep
        do ip=1,ipland  ! all land points in this nsib=3 loop
         iq=iperm(ip)
         isoil = isoilm(iq)
         fle=(wb(iq,1)-swilt(isoil))/(sfc(isoil)-swilt(isoil))          
         wetfac(iq)=max( 0.,min(1.,fle) )
        enddo   ! ip=1,ipland
      endif     ! (nbarewet==0)

      if(nbarewet==1)then  ! simplest bare
!cdir nodep
        do ip=1,ipland  ! all land points in this nsib=3 loop
         iq=iperm(ip)
         isoil = isoilm(iq)
         fle= wb(iq,1)/ssat(isoil)                                   
         wetfac(iq)=max( 0.,min(1.,fle) )
        enddo   ! ip=1,ipland
      endif     ! (nbarewet==1)

      if(nbarewet==2)then  ! jlm suggestion
!cdir nodep
        do ip=1,ipland  ! all land points in this nsib=3 loop
         iq=iperm(ip)
         isoil = isoilm(iq)
         frac=max(.01,tsigmf(iq))     ! jlm special for fle
         fle= (wb(iq,1)-frac*max( wb(iq,1),swilt(isoil) ))/         
     .               (ssat(isoil)*(1.-frac))                         
         wetfac(iq)=max( 0.,min(1.,fle) )
        enddo   ! ip=1,ipland
      endif     ! (nbarewet==2)

      if(nbarewet==3)then  ! jlm suggestion
!cdir nodep
        do ip=1,ipland  ! all land points in this nsib=3 loop
         iq=iperm(ip)
         isoil = isoilm(iq)
         frac=max(.01,tsigmf(iq))     ! jlm special for fle
         fle= (wb(iq,1)-frac*swilt(isoil) )/         
     .               (sfc(isoil)-frac*swilt(isoil))                    
         wetfac(iq)=max( 0.,min(1.,fle) )
        enddo   ! ip=1,ipland
      endif     ! (nbarewet==3)

      if(nbarewet==4)then  ! jlm, similar to Noilhan & Planton
!cdir nodep
        do ip=1,ipland  ! all land points in this nsib=3 loop
         iq=iperm(ip)
         isoil = isoilm(iq)
         fle=min( 1.,wb(iq,1)/sfc(isoil) )         
         wetfac(iq)=fle*fle*(3.-2.*fle)
        enddo   ! ip=1,ipland
      endif     ! (nbarewet==4)

      if(nbarewet==5)then  ! jlm, similar to Noilhan & Planton with swilt
!cdir nodep
        do ip=1,ipland  ! all land points in this nsib=3 loop
         iq=iperm(ip)
         isoil = isoilm(iq)
         fle=max( 0.,min( 1.,
     .             (wb(iq,1)-swilt(isoil))/(sfc(isoil)-swilt(isoil)) ) )
         wetfac(iq)=fle*fle*(3.-2.*fle)
        enddo   ! ip=1,ipland
      endif     ! (nbarewet==5)

      if(nbarewet==6)then  ! newer jlm
!cdir nodep
        do ip=1,ipland  ! all land points in this nsib=3 loop
         iq=iperm(ip)
         isoil = isoilm(iq)
         fle= max( 0.,min(1.,wb(iq,1)/ssat(isoil)) )
         wetfac(iq)=fle*fle*(2.2-1.2*fle)  ! .4 for fle=.5
        enddo   ! ip=1,ipland
      endif     ! (nbarewet==6)

      if(nbarewet==7)then  ! newest piecewise jlm (use with nalpha=1, beta)
!cdir nodep
        do ip=1,ipland  ! all land points in this nsib=3 loop
         iq=iperm(ip)
         isoil = isoilm(iq)
         wetfac(iq)=.2*(wb(iq,1)/swilt(isoil))**2
         if(wb(iq,1)>swilt(isoil))then
          wetfac(iq)=(.2*(sfc(isoil)-wb(iq,1))+
     .             .8*(wb(iq,1)-swilt(isoil)))/(sfc(isoil)-swilt(isoil))
         endif
         if(wb(iq,1)>sfc(isoil))then
           wetfac(iq)=(.8*(ssat(isoil)-wb(iq,1))+  
     .                   (wb(iq,1)-sfc(isoil)))/(ssat(isoil)-sfc(isoil))
         endif
        enddo   ! ip=1,ipland
      endif     ! (nbarewet==7)

      if(nbarewet==8)then  ! like NCAR but uses ssat
!cdir nodep
        do ip=1,ipland  ! all land points in this nsib=3 loop
         iq=iperm(ip)
         isoil = isoilm(iq)
         fle=(wb(iq,1)-swilt(isoil))/(ssat(isoil)-swilt(isoil))        
         wetfac(iq)=max( 0.,min(1.,fle) )
        enddo   ! ip=1,ipland
      endif     ! (nbarewet==8)

      if(nalpha==1)then    ! beta scheme
!cdir nodep
        do ip=1,ipland  ! all land points in this nsib=3 loop
         iq=iperm(ip)
         isoil = isoilm(iq)
         conw_fh=rho(iq)*taftfhg(iq)*hl    ! taftfhg(iq)=aftlandg*fhbg
         epot(iq) = conw_fh*(qsttg(iq)-qg(iq,1))
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
         qtgnet=  qsttg(iq)*wetfac(iq) -qg(iq,1)
         qtgair=  qsttg(iq)*wetfac(iq)-max(qtgnet,.1*qtgnet)
         eg2=   -conw_fh*qtgair
         eg1=    conw_fh*qsttg(iq)
!        degdw(iq)=eg1/ssat(isoil)                         ! not used in sib3
c        evaporation from the bare ground
         egg(iq)=eg1*wetfac(iq) +eg2
         epot(iq) = conw_fh*(qsttg(iq)-qg(iq,1))
         deg=wetfac(iq)*conw_fh*dqsttg(iq)
c        following reduces degdt by factor of 10 for dew
         degdt(iq)=.55*deg+sign(.45*deg,qtgnet)
        enddo   ! ip=1,ipland
      endif    ! (nalpha==1) .. else ..
      if((ntest==1.or.diag).and.mydiag.and.land(idjd))then 
          iq=idjd
          print *,'epot,egg,tgg1,snowd ',
     .             epot(iq),egg(iq),tgg(iq,1),snowd(iq)
          isoil = isoilm(iq)
          conw_fh=rho(iq)*taftfhg(iq)*hl    ! taftfhg(iq)=aftlandg*fhbg
          qtgnet=  qsttg(iq)*wetfac(iq) -qg(iq,1)
          qtgair=  qsttg(iq)*wetfac(iq)-max(qtgnet,.1*qtgnet)
          eg2=   -conw_fh*qtgair
          eg1=    conw_fh*qsttg(iq)
!         degdw(iq)=eg1/ssat(isoil)                   ! not used in sib3
c         evaporation from the bare ground
          egg_alph1=wetfac(iq)*epot(iq)
          print *,'then iq,isoil,conw_fh,qsttg,qtgair ',
     .             iq,isoil,conw_fh,qsttg(iq),qtgair
          print *,'eg1,eg2,wetfac ',eg1,eg2,wetfac(iq)
          print *,'epot,egg,egg_alph1 ',epot(iq),egg(iq),egg_alph1
      endif  ! (ntest==1)

!cdir nodep
      do ip=1,ipland  ! all land points in this nsib=3 loop
       iq=iperm(ip)
       if(snowd(iq)>1.)then
         egg(iq)=epot(iq)
         wetfac(iq)=1.   ! added jlm 18/3/04 to fix qgscrn inconsistency
         cls(iq)=1.+hlf/hl
       else
         egg(iq)=min(egg(iq),wb(iq,1)*zse(1)*1000.*hl/dt)
         cls(iq)=1.
       endif  ! (snowd(iq)>1.)
      enddo   ! ip=1,ipland
 
c     print *,'before nsigmf'

      if(nsigmf==0)then  ! original
!cdir nodep
        do ip=1,ipland  ! all land points in this nsib=3 loop
         iq=iperm(ip)
         ga(iq)=-slwa(iq)-rgg(iq)-fgg(iq)-cls(iq)*egg(iq)
         dgdtg(iq)=-dirad(iq)-dfgdt(iq)-cls(iq)*degdt(iq)
        enddo   ! ip=1,ipland
      endif     ! (nsigmf==0)

      if(nsigmf==1)then  ! jlm preferred
!cdir nodep
        do ip=1,ipland  ! all land points in this nsib=3 loop
        iq=iperm(ip)
        ga(iq)=(-slwa(iq)-rgg(iq)-fgg(iq)-cls(iq)*egg(iq))*
     .                                  (1.-tsigmf(iq))
!       dgtdg is used in soilsnow
        dgdtg(iq)=-(dirad(iq)+dfgdt(iq)+cls(iq)*degdt(iq))*
     .                                  (1.-tsigmf(iq))
        enddo   ! ip=1,ipland
      endif     ! (nsigmf==1)

! ----------------------------------------------
      if(itnmeth==0) then  ! old, not vectorized
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
          if( zo(iq) < .5 ) sstar = 150.
          ff= 1.1*sgsave(iq)/(rlai(iq)*sstar)
c         ff= 1.1*slwa(iq)/(rlai(iq)*sstar)
          rsi = rsmin(iq) * rlai(iq)
          f1= (1.+ff)/(ff+rsi/5000.)
          den=sfc(isoil)-swilt(isoil)                          ! sib3
          wbav=(max(0.,froot(iveg,1)*(wb(iq,1)-swilt(isoil)))+
     &         max(0.,froot(iveg,2)*(wb(iq,2)-swilt(isoil)))+
     &         max(0.,froot(iveg,3)*(wb(iq,3)-swilt(isoil)))+
     &         max(0.,froot(iveg,4)*(wb(iq,4)-swilt(isoil)))+
     &         max(0.,froot(iveg,5)*(wb(iq,5)-swilt(isoil)))+
     &         max(0.,froot(iveg,6)*(wb(iq,6)-swilt(isoil)))   )/den
          f2=max(1. , .5/ max( wbav,1.e-7)) ! N.B. this is equiv to next 2 (jlm)
c         f2=1.0
c         if(wbav<0.5) f2=max(1.0 , 0.5/ max( wbav,1.e-7))
          f4=max(1.-.0016*(tstom(iq)-t(iq,1))**2 , .05) ! zero for delta_t=25
          airr(iq) = 1./taftfh(iq)
!         rmstep=dt/60.
          cc(iq) =min(condx(iq) , 4./(1440. *60./dt))  ! jlm speedup for 4 mm/day
c                       depth of the reservoir of water on the canopy
          rmcmax(iq) = max(0.5,srlai) * .1
          omc(iq) = rmc(iq)  ! value from previous timestep as starting guess
          f3=max(1.-.00025*(establ(t(iq,1))-qg(iq,1)*ps(iq)/.622), .05)
          res(iq)=max(30.,rsmin(iq)*f1*f2/(f3*f4))
          if(ntest==1.and.iq==idjd)then
           print *,'rlai,srlai,wbav,den ',rlai(iq),srlai,wbav,den
           print *,'f1,f2,f3,f4 ',f1,f2,f3,f4
           print *,'ff,f124,rsi,res ',ff,f1*f2/f4,rsi,res(iq)
          endif

          otgf(iq)=tgf(iq)
           do icount=1,5        ! original iteration method
c                                            transpiration
            rmc(iq) = omc(iq)
            esatf = establ(tgf(iq))
            qsatgf=.622*esatf/(ps(iq)-esatf)
c                                                     wet evaporation
            ewww(iq) = rho(iq)*taftfh(iq) *(qsatgf-qg(iq,1))
            if(qsatgf>=qg(iq,1)) then
              ewww(iq)  = min(rmc(iq)/dt , rmc(iq)*ewww(iq)/rmcmax(iq) )
            endif         ! qsatgf>=qg(iq,1)
            rmc(iq)=omc(iq)+cc(iq) -ewww(iq)*dt
c                         precipitation reaching the ground under the canopy
c                                            water interception by the canopy
            condxg(iq)=max(condx(iq)-cc(iq)+
     .                 max(0.,rmc(iq)-rmcmax(iq)),0.) ! keep here
            rmc(iq) = min( max(rmc(iq),0.), rmcmax(iq))
!           beta =  min(rmc(iq)/rmcmax(iq),1.)   ! not needed
            beta =      rmc(iq)/rmcmax(iq)
            if( qsatgf < qg(iq,1) ) beta = 1.

            Etr=rho(iq)/(airr(iq) +res(iq))*(qsatgf-qg(iq,1))
            betetrdt =(1.-beta)*Etr*dt*tsigmf(iq)
            evapfb1(iq)=min(betetrdt*froot(iveg,1),max(0.,
     &                  (wb(iq,1)-swilt(isoil))*zse(1)*1000.))
            evapfb2(iq)=min(betetrdt*froot(iveg,2),max(0.,
     &                  (wb(iq,2)-swilt(isoil))*zse(2)*1000.))
            evapfb3(iq)=min(betetrdt*froot(iveg,3),max(0.,
     &                  (wb(iq,3)-swilt(isoil))*zse(3)*1000.))
            evapfb4(iq)=min(betetrdt*froot(iveg,4),max(0.,
     &                  (wb(iq,4)-swilt(isoil))*zse(4)*1000.))
            evapfb5(iq)=min(betetrdt*froot(iveg,5),max(0.,
     &                  (wb(iq,5)-swilt(isoil))*zse(5)*1000.))
            evapfb6(iq)=min(betetrdt*froot(iveg,6),max(0.,
     &                  (wb(iq,6)-swilt(isoil))*zse(6)*1000.))
            evapfb(iq)=(evapfb1(iq)+evapfb2(iq)+evapfb3(iq)+evapfb4(iq)+
     &                  evapfb5(iq)+  evapfb6(iq))/tsigmf(iq)
            evapxf(iq) = (evapfb(iq)/dt + ewww(iq))*hl
c           evapw =  ewww(iq)*hl  ! not used
c                                              sensible heat flux
            prz = rho(iq)*taftfh(iq)*cp
            if(newfgf==0)fgf(iq)=prz*(tgf(iq)-theta(iq))  ! original/usual
            if(newfgf==1)fgf(iq)=prz*(tgf(iq)-tscrn(iq))
            if(newfgf==2)fgf(iq)=rho(iq)*aft(iq)*cp*
     .                                (tgf(iq)-tscrn(iq))
c           fgf(iq)=min(800.,fgf(iq)) !***********************
            rdg(iq) =  stefbo*tgf(iq)**4
!           residf(iq) = (-slwa(iq) -rdg(iq))*(1.-extin(iq))
!    .                    -fgf(iq) -evapxf(iq)
            residf(iq) = -slwa(iq) -rdg(iq) -fgf(iq) -evapxf(iq)
            if( abs(residf(iq)) < 3. )go to 110
            dirad1 = 4.*rdg(iq)/tgf(iq)
c                                                     Calculate dE/dTg
            dqg=qsatgf*hlars/tgf(iq)**2
            devf= hl*rho(iq)*dqg*( (1.-beta)/(airr(iq) + res(iq))
     &                          +beta/airr(iq) ) ! re-factored by jlm
            residp = -(dirad1 + devf + prz)
            sstgf = tgf(iq)
            tgf(iq)=tgf(iq)-residf(iq)/residp
            if(ntest==1.and.iq==idjd)then
              print *,'icount,omc,rmc,otgf ',
     .                 icount,omc(iq),rmc(iq),otgf(iq)
              print *,'tfg,residf,evapxf,fgf ',
     .                 tgf(iq),residf(iq),evapxf(iq),fgf(iq)
            endif
            tgf(iq)=min(tgf(iq),sstgf+ 1.5)   ! jlm speedup
            tgf(iq)=max(tgf(iq),sstgf- 1.5)   ! jlm speedup

           enddo   !  icount=1,5
           if(ntest==1)print *,'iq,otgf(iq),tgf,residf '
     .                           ,iq,otgf(iq),tgf(iq),residf(iq)
c          tgf(iq)=0.5*(otgf(iq)+tgf(iq))
110        continue
       enddo  !  ip=1,ipland
      endif  ! (itnmeth==0) 


      if(itnmeth>0) then  ! new, vectorized
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
!       if( zo(iq) < .5 ) sstar = 150.
        sstar=90.+sign(60.,.5-zo(iq))  ! equiv to above 2 lines
        ff= 1.1*sgsave(iq)/(rlai(iq)*sstar)
c       ff= 1.1*slwa(iq)/(rlai(iq)*sstar)
        rsi = rsmin(iq) * rlai(iq)
        f1= (1.+ff)/(ff+rsi/5000.)
        den=sfc(isoil)-swilt(isoil)                          ! sib3
        wbav=(max(0.,froot(iveg,1)*(wb(iq,1)-swilt(isoil)))+
     &        max(0.,froot(iveg,2)*(wb(iq,2)-swilt(isoil)))+
     &        max(0.,froot(iveg,3)*(wb(iq,3)-swilt(isoil)))+
     &        max(0.,froot(iveg,4)*(wb(iq,4)-swilt(isoil)))+
     &        max(0.,froot(iveg,5)*(wb(iq,5)-swilt(isoil)))+
     &        max(0.,froot(iveg,6)*(wb(iq,6)-swilt(isoil)))   )/den

        f2=max(1. , .5/ max( wbav,1.e-7)) ! N.B. this is equiv to next 2 (jlm)
c       f2=1.0
c       if(wbav<0.5) f2=max(1.0 , 0.5/ max( wbav,1.e-7))
        f4=max(1.-.0016*(tstom(iq)-t(iq,1))**2 , .05) ! zero for delta_t=25
        airr(iq) = 1./taftfh(iq)
        cc(iq) =min(condx(iq) , 4./(1440. *60./dt))  ! jlm speedup for 4 mm/day
c                     depth of the reservoir of water on the canopy
        rmcmax(iq) = max(0.5,srlai) * .1
        omc(iq) = rmc(iq)  ! value from previous timestep as starting guess
        f3=max(1.-.00025*(establ(t(iq,1))-qg(iq,1)*ps(iq)/.622),.05)
        res(iq)=max(30.,rsmin(iq)*f1*f2/(f3*f4))
        if(ntest==1.and.iq==idjd)then
          print *,'rlai,srlai,wbav,den ',rlai(iq),srlai,wbav,den
          print *,'f1,f2,f3,f4 ',f1,f2,f3,f4
          print *,'ff,f124,rsi,res ',ff,f1*f2/f4,rsi,res(iq)
          print *,'qg,qfg,qlg ',qg(iq,1),qfg(iq,1),qlg(iq,1)
        endif
        otgf(iq)=tgf(iq)
        tgfnew(iq)=tgf(iq)
        delta_tx(iq)=10.   ! just to supply max change of 5 deg first time
        evapfb1a(iq)=max(0.,wb(iq,1)-swilt(isoil)) *zse(1)*1000.
        evapfb2a(iq)=max(0.,wb(iq,2)-swilt(isoil)) *zse(2)*1000.
        evapfb3a(iq)=max(0.,wb(iq,3)-swilt(isoil)) *zse(3)*1000.
        evapfb4a(iq)=max(0.,wb(iq,4)-swilt(isoil)) *zse(4)*1000.
        evapfb5a(iq)=max(0.,wb(iq,5)-swilt(isoil)) *zse(5)*1000.
        evapfb6a(iq)=max(0.,wb(iq,6)-swilt(isoil)) *zse(6)*1000.
       enddo   ! ip loop

       do icount=1,itnmeth     ! jlm new iteration
c                                            transpiration
!cdir nodep
        do ip=1,ipland  ! all land points in this nsib=3 loop
         iq=iperm(ip)
         esatf = establ(tgfnew(iq))
         qsatgf=.622*esatf/(ps(iq)-esatf)
c                                                  wet evaporation
         ewwwa = rho(iq) *(qsatgf-qg(iq,1))/airr(iq)
!        if(qsatgf>=qg(iq,1)) then  ! no dew
!          ewww(iq)  = min(omc(iq)/dt , ewww(iq)*omc(iq)/rmcmax(iq) )
!        endif         ! qsatgf>=qg(iq,1)
!        above 3 lines are equivalent to next line
c        ewww(iq)=min(ewww(iq),omc(iq)/dt,ewww(iq)*omc(iq)/rmcmax(iq))
c        rmcav=.5*(omc(iq)+rmc(iq))
c        ewww(iq)=min(ewwwa,rmcav/dt,ewwwa*rmcav/rmcmax(iq))
c        rmc(iq)=omc(iq)+cc(iq) -ewww(iq)*dt
!        to reduce risk of underflow, replace above few lines by following
         ewww(iq)=min(dt*ewwwa,omc(iq),dt*ewwwa*omc(iq)/rmcmax(iq))
!                     dew(-ve), no_dew ,        no_dew
!        rmc is reservoir on leaf
         rmc(iq)=(omc(iq)-ewww(iq)) +cc(iq)
         ewww(iq)=ewww(iq)/dt  ! these changes on 19/1/06 jlm

c                      precipitation reaching the ground under the canopy
c                                         water interception by the canopy
         condxg(iq)=max(condx(iq)-cc(iq)+max(0.,rmc(iq)-rmcmax(iq)),0.) ! keep 
         rmc(iq) = min( max(0.,rmc(iq)), rmcmax(iq))
!        beta =  min(rmc(iq)/rmcmax(iq),1.)   ! min not needed
         beta =      rmc(iq)/rmcmax(iq)
!!       if( qsatgf < qg(iq,1) ) beta = 1.   ! i.e. dew
!        Etr=rho(iq)*(qsatgf-qg(iq,1))/(airr(iq) +res(iq))
         Etr=rho(iq)*max(0.,qsatgf-qg(iq,1))/(airr(iq) +res(iq))  ! jlm
         betetrdt =(1.-beta)*Etr*dt*tsigmf(iq)   ! fixed 23/5/01
         evapfb1(iq)=min(betetrdt*froot(iveg,1),evapfb1a(iq))
         evapfb2(iq)=min(betetrdt*froot(iveg,2),evapfb2a(iq))
         evapfb3(iq)=min(betetrdt*froot(iveg,3),evapfb3a(iq))
         evapfb4(iq)=min(betetrdt*froot(iveg,4),evapfb4a(iq))
         evapfb5(iq)=min(betetrdt*froot(iveg,5),evapfb5a(iq))
         evapfb6(iq)=min(betetrdt*froot(iveg,6),evapfb6a(iq))
         evapfb(iq)=(evapfb1(iq)+evapfb2(iq)+evapfb3(iq)+evapfb4(iq)+
     &               evapfb5(iq)+evapfb6(iq))/tsigmf(iq)

         evapxf(iq) = (evapfb(iq)/dt + ewww(iq))*hl  ! converting to W/m**2
         prz = rho(iq)*cp*taftfh(iq)
         if(newfgf==0)fgf(iq) = prz*(tgfnew(iq)-theta(iq))  ! original/usual
         if(newfgf==1)fgf(iq) = prz*(tgfnew(iq)-tscrn(iq))
         if(newfgf==2)fgf(iq)=rho(iq)*aft(iq)*cp*
     .                             (tgfnew(iq)-tscrn(iq))
!	  limit extreme fgf to avoid undue tgf oscillations  June '04
         fgf(iq)=max(-1000.,min(fgf(iq),1000.))
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
!!       if(icount>1)then
!!         delta_t=sign(min(abs(delta_t),                    !  jlmnew
!!   .                     .5*abs(tgfnew(iq)-tgf(iq))) , delta_t0)
!!       endif      ! icount>1
!        above few lines equivalent to next one:
         delta_t=sign(min(abs(delta_t0),.5*abs(delta_tx(iq))),delta_t0)
         tgfnew(iq)=tgfnew(iq)+delta_t
         delta_tx(iq)=tgfnew(iq)-otgf(iq)
c!       following to limit change to 8 degrees as fh may be poor  May '04
c        tgfnew(iq)=otgf(iq)+
c     .              sign(min(abs(delta_tx(iq)),8.),delta_tx(iq))
         enddo  ! ip loop
         if((ntest==2.or.diag).and.mydiag)then 
          if(ktau==99999)then  ! to track rmc & ewww underflow
            print *,'sflux icount = ',icount
            do ip=1,ipland  
             iq=iperm(ip)
             print *,'iq,otgf,tgfnew,ewww ',
     &                iq,otgf(iq),tgfnew(iq),ewww(iq)
            enddo
          endif
          if(land(idjd))then
           iq=idjd
           print *,'ktau,icount,iq,omc,cc ',
     &              ktau,icount,iq,omc(iq),cc(iq)
           print *,'rmc,rmcmax,ewww ',rmc(iq),rmcmax(iq),ewww(iq)
           esatf = establ(tgfnew(iq))  ! value for next itn
           qsatgf=.622*esatf/(ps(iq)-esatf)
           ewwwa = rho(iq) *(qsatgf-qg(iq,1))/airr(iq)
           print *,'esatf,qsatgf,ewwwa ',esatf,qsatgf,ewwwa
           prz = rho(iq)*cp*taftfh(iq)
           beta =      rmc(iq)/rmcmax(iq)
           devf= (hl*hlars/300.**2)*qsatgf*(1.-beta)/res(iq)
           dirad1 = 4.*rdg(iq)/300.
           print *,'beta,airr,res ',beta,airr(iq),res(iq)
           print *,'dirad1,devf,prz ',dirad1,devf,prz
           print *,'theta,tscrn,slwa ',theta(iq),tscrn(iq),slwa(iq)
           print *,'taftfh,condxg ',taftfh(iq),condxg(iq)
           print *,'rdg,fgf,evapxf,evapfb ',
     .              rdg(iq),fgf(iq),evapxf(iq),evapfb(iq)
           print *,'delta_tx ',delta_tx(iq)
           print *,'otgf,tgfnew,residf ',otgf(iq),tgfnew(iq),residf(iq)
          endif  ! (land(idjd))
         endif   ! ((ntest==2.or.diag).and.mydiag)
       enddo     !  icount=1,5
      endif      ! (itnmeth>0) 

!cdir nodep
      do ip=1,ipland  ! all land points in this nsib=3 loop
       iq=iperm(ip)
       if(rmc(iq)<1.e-10)rmc(iq)=0.  ! to avoid underflow 24/1/06
       if(tsigmf(iq) <= .0101) then
         condxpr(iq)=condx(iq)
         evapfb(iq) = 0.
         evapxf(iq) = egg(iq)
         fgf(iq)  = fgg(iq)
         rdg(iq)=rgg(iq)
         tgf(iq) = tss(iq)
       else
         tgf(iq)=tgfnew(iq)
         wb(iq,1)=wb(iq,1)-evapfb1(iq)/(zse(1)*1000.)
         wb(iq,2)=wb(iq,2)-evapfb2(iq)/(zse(2)*1000.)
         wb(iq,3)=wb(iq,3)-evapfb3(iq)/(zse(3)*1000.)
         wb(iq,4)=wb(iq,4)-evapfb4(iq)/(zse(4)*1000.)
         wb(iq,5)=wb(iq,5)-evapfb5(iq)/(zse(5)*1000.)
         wb(iq,6)=wb(iq,6)-evapfb6(iq)/(zse(6)*1000.)
         condxpr(iq)=(1.-tsigmf(iq))*condx(iq)+ tsigmf(iq)*condxg(iq)
         if(ntest==1.and.abs(residf(iq))>10.)
     .      print *,'iq,otgf(iq),tgf,delta_tx,residf '
     .              ,iq,otgf(iq),tgf(iq),delta_tx(iq),residf(iq)
       endif          ! tsigmf <= .01   ... else ...
       fev(iq)=evapfb(iq)/dt*hl*tsigmf(iq) ! passed to soilsnow to update wb
       fes(iq)=(1.-tsigmf(iq))*egg(iq)*cls(iq)  ! also passed to soilsnow
       otgsoil(iq)=isflag(iq)*tggsn(iq,1) + (1-isflag(iq))*tgg(iq,1)
      enddo  !  ip=1,ipland
      if((ntest==1.or.diag).and.mydiag.and.land(idjd))then 
         iq=idjd
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
         print *,'qg1,qsttg ',qg(iq,1),qsttg(iq)
         print *,'dfgdt,taftfhg,rho ',dfgdt(iq),taftfhg(iq),rho(iq)
         print *,'rmc,rmcmax(iq) ',rmc(iq),rmcmax(iq)
         if(abs(tgf(iq)-otgf(iq))>4.9)
     .      write(6,"('ktau,iq,otgf,tgf,dtgf,t1,t2',i4,i6,5f8.2)")
     .        ktau,iq,otgf(iq),tgf(iq),tgf(iq)-otgf(iq),t(iq,1),t(iq,2)
      endif
!-------------------------------------
c     print *,'before soilsnow'
      call soilsnowv
c     print *,'after soilsnow'
!cdir nodep
      do ip=1,ipland  ! all land points in this nsib=3 loop
       iq=iperm(ip)
c       iveg=ivegt(iq)
c       isoil = isoilm(iq)
c       if( tsigmf(iq)> .01)then
c         evapf=fev(iq)*dt/(hl*tsigmf(iq))
c         evapxf(iq) = (evapf/dt + ewww(iq))*hl
c       endif  !  ( tsigmf(iq)> .01)
c       if(eg(iq)<0.)then
c         print *,'iq,eg,snowd,tgg ',iq,eg(iq),snowd(iq),tgg(iq,1)
c       endif

       if(isflag(iq)==0) then
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
       if(snowd(iq)>1.)then
         eg(iq)=tsigmf(iq)*evapxf(iq) + egg(iq)
       else
         eg(iq) = tsigmf(iq)*evapxf(iq) + (1. - tsigmf(iq))*egg(iq)
       endif
       if(nsigmf==2)then
         fg(iq)=tsigmf(iq)*fgf(iq)+fgg(iq)
       else
         fg(iq)=tsigmf(iq)*fgf(iq)+(1.-tsigmf(iq))*fgg(iq)
       endif
       rnet(iq)=-slwa(iq)-(1.-tsigmf(iq))*rgg(iq)-tsigmf(iq)*rdg(iq)
!      rnet(iq)=-slwa(iq)-(1.-tsigmf(iq))*rgg(iq)
!    .   -tsigmf(iq)*(extin(iq)*rgg(iq) + (1.-extin(iq))*rdg(iq)) ! 9/3/99

        tgss=isflag(iq)*tggsn(iq,1) + (1-isflag(iq))*tgg(iq,1)  ! jlm
        if(tsigmf(iq)<= .01) then
          tss(iq) = tgss
          tgf(iq) = tgss
        else
          tss(iq)=tsigmf(iq)*tgf(iq)+(1.-tsigmf(iq))*tgss
        endif       ! tsigmf<= .01
        es = establ(tss(iq))     !  from 27/12/05
        qsttg(iq)= .622*es/(ps(iq)-es)  ! recal for scrnout, esp. snow    

!eak!      following are just diagnostic arrays used by eak
!eak       if(ndiag_arr==1)then
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
!eak       endif  ! (ndiag_arr==1)

      enddo   ! ip=1,ipland
      xx1=tggsn(:,1)
      xx2=tggsn(:,2)
      xx3=tggsn(:,3)
      where (isflag == 0 )
         xx1=tgg(:,1)
         xx2=tgg(:,2)
         xx3=tgg(:,3)
      end where
      zz1 = 0.
      swnet=-slwa
      lwnet=rgsave-(1.-tsigmf)*rgg-tsigmf*rdg
!      do k=1,1
      do k=1,5
       ijd =list(k,1)
       WRITE(98, '(i5,"xr",i7,67e14.5)') 
     . ijd,ktauplus+ktau,yyy0,sgsave(ijd)/(1.-alb(ijd)),-rgsave(ijd), 
     .    condx(ijd), !1-5
     .    t(ijd), ps(ijd),vmod(ijd),qg(ijd),zz1(ijd),  ! 6-11
     .    -slwa(ijd),lwnet(ijd), 
     .    tsigmf(ijd)*evapxf(ijd),tsigmf(ijd)*fgf(ijd),!12-14
     .    rnet(ijd), tss(ijd),tgf(ijd),  !15-17
     .    taftfh(ijd),wb(ijd,1),wb(ijd,2), !18-20
     .    wb(ijd,3),wb(ijd,4), wb(ijd,5),wb(ijd,6), !21-24
     .    xx1(ijd),xx2(ijd),xx3(ijd),        !25-27
     .    tgg(ijd,1),tgg(ijd,2), tgg(ijd,3),tgg(ijd,4), !28-31
     .    tgg(ijd,5),tgg(ijd,6),    !32-33
     .    -slwa(ijd)-rgg(ijd),egg(ijd),fgg(ijd),ga(ijd),  !34-37
     .    zz1(ijd),zz1(ijd),zz1(ijd),zz1(ijd),zz1(ijd), !38-42
     .    -slwa(ijd)-rdg(ijd),tsigmf(ijd)*evapxf(ijd),fev(ijd),
     .     ewww(ijd)*hl,fgf(ijd),  !43-47
     .    zz1(ijd),vlai(ijd),zz1(ijd), 
     .    zz1(ijd),zz1(ijd), 
     .    -slwa(ijd)-rgg(ijd)-fgg(ijd)-fes(ijd)*cls(ijd)-ga(ijd), 
     .    runoff(ijd),rnof1(ijd),rnof2(ijd),zz1(ijd),snowd(ijd), 
     .    zz1(ijd),alb(ijd),alb(ijd), 
     .    taftfhg(ijd), ustar(ijd),cduv(ijd), pblh(ijd)
       enddo



      if((ntest==1.or.diag).and.mydiag.and.land(idjd))then 
        iq=idjd
        print *,'even further down sib3 after soilsnowv'
        print *,'tgg ',(tgg(iq,k),k=1,ms)
        print *,'wb ',(wb(iq,k),k=1,ms)
        print *,'isflag,snowd ',isflag(iq),snowd(iq)
        print *,'evapfb,fev,ewww ',evapfb(iq),fev(iq),ewww(iq)
        print *,'tsigmf,evapxf,egg ',tsigmf(iq),evapxf(iq),egg(iq)
        print *,'deltat,degdt,wb,zse ',
     .           tgg(iq,1)-otgsoil(iq),degdt(iq),wb(iq,1),zse(1)
        print *,'eg,fg ',eg(iq),fg(iq)
      endif

      return
      end

       

      subroutine sib4     ! new version of sib1 with soilsnowv

      use cc_mpi
      use zenith_m
      use pack_unpack_m 
      USE define_types
      USE air_module
      USE canopy_module
      USE carbon_module
      USE cbm_module
      USE soil_snow_module
      INTEGER, PARAMETER          :: mp_max = 99 ! max. no. geographic points
      TYPE (air_type)             :: air
      TYPE (bgc_pool_type)        :: bgc
      TYPE (canopy_type)          :: canopy
      TYPE (met_type)          :: met
      TYPE (balances_type)     :: bal
      TYPE (radiation_type)       :: rad
      TYPE (roughness_type)       :: rough
      type (soil_parameter_type)      :: soil       ! soil parameters
!      TYPE (soil_type)            :: soil
      TYPE (soil_snow_type)       :: ssoil
      TYPE (sum_flux_type)        :: sum_flux
!      TYPE (veg_type)             :: veg
      type (veg_parameter_type)       :: veg        ! vegetation parameters


!      include 'newmpar.h'
      include 'aalat.h'    ! slwa
      include 'arrays.h'
      include 'carbpools.h'
!      include 'const_phys.h'
      include 'dates.h' ! ktime,kdate,timer,timeg,xg,yg
      include 'extraout.h'
      include 'filnames.h'
      include 'latlong.h'  ! rlatt,rlongg
      include 'map.h'      ! id,jd,idjd
      include 'morepbl.h'
      include 'nsibd.h'    ! rsmin,sigmf,tgf,ssdn,res,rmc,tsigmf
      include 'parm.h'
      include 'permsurf.h'
      include 'pbl.h'
      include 'prec.h'
      include 'screen.h'   ! tscrn etc
      include 'sigs.h'
      include 'soil.h'     ! ... zmin zolnd zolog sice alb
!      include 'soilroot.h'   ! ssat, clay,..
      include 'soilsnow.h' ! 
      include 'soilsnin.h'
      include 'soilv.h'    ! ssat, clay,..
      include 'vegpar.h' ! 
!      include 'tracers.h'  ! ngas, nllp, ntrac
      include 'trcom2.h'   ! nstn,slat,slon,istn,jstn
!      include 'establ.h'
!                     met forcing for CBM
      common/work2/dirad(ifull),dfgdt(ifull),degdt(ifull)
     . ,wetfac(ifull),degdw(ifull),cie(ifull)
     . ,factch(ifull),qsttg(ifull),rho(ifull),zo(ifull)
     . ,aft(ifull),fh(ifull),ri(ifull),theta(ifull)
     . ,gamm(ifull),rg(ifull),vmod(ifull),dummwk2(ifull)

      real fjd, r1, dlt, slag, dhr, coszro2(ifull),taudar2(ifull) ! for calculation of zenith angle

      integer list(20,2),imonth(12),iwrk(ifull)
      data imonth /31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31./
      dimension ndoy(12)   ! days from beginning of year (1st Jan is 0)
      data ndoy/ 0,31,59,90,120,151,181,212,243,273,304,334/
      real ssumcbmfl(14,20)
      save list,ktauplus,iswt
      data iswt/0/
!     rml added values for 2000-2004 using scripps records from cdiac
!     0.75*mauna loa + 0.25*south pole
      real co2for(0:104)
      data co2for/296.0049,296.3785,296.7731,297.1795,297.5887,297.9919
     &           ,298.3842,298.7654,299.1349,299.4925,299.838 ,300.1709
     &           ,300.491 ,300.801 ,301.106 ,301.4113,301.7205,302.0357
     &           ,302.3587,302.6915,303.036,303.3941,303.7677,304.1587
     &           ,304.569,304.9971,305.4388,305.8894,306.3444,306.7992
     &           ,307.2494,307.6902,308.117,308.521,308.8895,309.2135
     &           ,309.4877,309.7068,309.8658,309.9668,310.019,310.0358
     &           ,310.035,310.0345,310.0522,310.0989,310.1794,310.2977
     &           ,310.4581,310.6661,310.928,311.2503,311.6395,312.1015
     &           ,312.6341,313.227,313.8694,314.5506,315.2599,315.9866
     &           ,316.7167,317.4268,318.106,318.7638,319.4402,320.1901
     &           ,321.0398,321.9713,322.9779,324.045,325.1356,326.2445
     &           ,327.3954,328.5711,329.723,330.8865,332.1331,333.5012
     &           ,334.9617,336.4614,337.9588,339.4225,340.8724,342.3488
     &           ,343.8625,345.421,347.0481,348.7452,350.4397,352.0169
     &           ,353.4269,354.6917,355.8849,357.1348,358.5514,360.1441
     &           ,361.8578,363.6479,365.4682,367.2137
     &           ,368.87,370.35,372.49,374.93,376.69/

c
c      set meteorological forcing
c
       list(1,1)=8290
       list(2,1)=9061
       list(3,1)=9093
       list(4,1)=9293
       list(5,1)=9540

       jyear=kdate/10000
       ipco2=jyear-1900
       ipco2 = min(ipco2,100)
       jmonth=(kdate-jyear*10000)/100
       jday=kdate-jyear*10000-jmonth*100
       jhour=ktime/100
       jmin=ktime-jhour*100
       mstart=1440*(ndoy(jmonth)+jday-1) + 60*jhour + jmin ! mins from start of year
       nperday=nint(24.*3600./dt)
       ktauplus=0
       do k=1,jmonth-1
        ktauplus = ktauplus + imonth(k)*nperday
       enddo
       print *,'jyear,jmonth',jyear,jmonth,ktauplus,ipco2
!       timer contains number of hours since the start of the run.
!       mins = 60 * timer + mstart
!       mtimer contains number of minutes since the start of the run.
        mins = mtimer + mstart
        bpyear = 0.
        fjd = float(mod(mins,525600))/1440.  ! 525600 = 1440*365
        call solargh(fjd,bpyear,r1,dlt,alp,slag)
        dhr = 1.e-6
        call zenith(fjd,r1,dlt,slag,rlatt,
     &               rlongg,dhr,ifull,coszro2,taudar2)
        jmonth=(kdate-jyear*10000)/100
        jday=kdate-jyear*10000-jmonth*100
        jhour=ktime/100
        jmin=ktime-jhour*100
        mstart=1440*(ndoy(jmonth)+jday-1) + 60*jhour + jmin ! mins from start of year
!       timer contains number of hours since the start of the run.
!       mins = 60 * timer + mstart
!       mtimer contains number of minutes since the start of the run.
        mins = mtimer + mstart

        do 60 ip=1,ipland  ! all land points in this nsib=2 loop
         iq=iperm(ip)
         rad%latitude(ip)=alat(iq)
         met%doy(ip)= float(mod(mins,525600))/1440.  ! 525600 = 1440*365
         met%hod(ip)=(met%doy(ip)-int(met%doy(ip)))*24.0 + along(iq)/15
         if(met%hod(ip).gt.24.0) met%hod(ip)=met%hod(ip)-24.0
         rough%za(ip) = -287.*t(iq,1)*log(sig(1))/9.80616   ! reference height
         met%fsd(ip) = sgsave(iq)/(1.-alb(iq))! short wave down (positive) W/m^2
         met%fld(ip) = -rgsave(iq)              ! long wave down     ---"---
         met%tc(ip) = theta(iq)-273.16   ! temperature at reference height
         met%tk(ip) = theta(iq)             ! temperature at reference height
         met%qv(ip) = qg(iq,1)                   ! specific humidity in kg/kg
         met%ua(ip) = max(1.,vmod(iq))              ! total wind speed at ref height
         met%pmb(ip) = .01*ps(iq)        ! pressure in mb at ref height
         if (ktau.eq.1) met%precip(ip) = 0.
         met%precip(ip) = condx(iq)
         met%precips(ip) = 0.0          ! in mm not mm/sec
         if( met%tc(ip) .le. 0.0 ) met%precips(ip) = condx(iq)
         met%ca(ip) = co2for(ipco2)*1.e-6           ! co2 concentration needed for CBM 
         met%coszen(ip) =max(1.e-8,coszro2(iq))   ! use instantenous value
 60     continue

        kstart = 1
        call cbm_pack(air, bgc, canopy, met, bal, rad,
     &          rough, soil, ssoil, sum_flux, veg)
 
!      ijd = 2468
!      print *,'sbmmet',ktau,met%fsd(ijd),met%fld(ijd),met%precip(ijd),
!     &  met%tc(ijd),met%tk(ijd),met%ua(ijd),met%ca(ijd),met%coszen(ijd),
!     &    met%coszen(ijd),met%precip(ijd),met%qv(ijd)
!      kk=1
!      do iq=1,ifull
!       do kk=1,5
!       if(iq.eq.list(kk,1)) then 
!       iq=list(kk,1)
!       ip=list(kk,2)
!       print 1192,iq,ivegt(iq),ktau,
!     &    ,met%fsd(ip),met%fld(ip),met%precip(ip),
!     &    met%tk(ip),met%ua(ip),
!     &    met%pmb(ip),vlai(iq),
!!       kk=kk+1
!1192   format('bef cbm',i6,i3,i6,f6.2,f6.2,f7.1,f6.1,2f8.4,f7.2,
!     &    f7.2,f6.1,f9.6,f6.0,e10.3,f6.3,f6.3)
!!       endif
!      enddo

       CALL cbm(ktau, kstart, ntau, ktauplus+ktau, dt, air, bgc, 
     &     canopy, met, bal, rad, rough, soil, ssoil, sum_flux, veg)

       call cbm_unpack(ktauplus+ktau, air, bgc, canopy, met, bal, rad,
     &          rough, soil, ssoil, sum_flux, veg)

!      evap = dt*eg/hl
        
      return
      end


      subroutine cbmrdn(imonth)
c
c     reads the parameters required by land surface scheme 
c

      USE define_types
      include 'arrays.h'
      include 'carbpools.h'
      include 'dates.h'     ! dt,kdate,iyd
      include 'nsibd.h'     ! ivegt,isoilm
      include 'parm.h'      ! id,jd,idjd
      include 'permsurf.h'
      include 'map.h'      ! id,jd,idjd
      include 'pbl.h'       ! tss
      include 'sigs.h'
      include 'soil.h'      ! land, zolnd passed to indata
      include 'soilv.h'
      include 'soilsnow.h'  ! initial soil temperatures and moistures
      include 'vegpar.h'

      integer jyear,jmonth
      character*80 comments
      character*2 chflag

      open(unit=8,file="vegtype.dat",status='old')
      read(8,*) comments
      print *,'read CASA vegetation types'
      do ii=1,ifull
       read(8,*) iq,i,j,rl1,rl2,ivold,ivnew
       ivegt(iq)=ivnew
!       if(land(iq)) print *,iq,i,j,rl1,rl2,ivold,ivegt(iq)
      enddo ! ii
      close(8)

      open(unit=8,file="veg_parm.txt",status='old')
      read(8,*) comments
      write(*,802) comments
802   format(1x,a80)
      call comskp(8)
      read(8,*) nveg
      print *,'nveg=',nveg
      call comskp(8)
      read(8,*) (canst1(jveg),jveg=1,nveg)
      print *, 'canst1',(canst1(jveg),jveg=1,nveg)
!      read(8,*) (cansto(jveg),jveg=1,nveg)
!      print *, 'cansto',(cansto(jveg),jveg=1,nveg)
      read(8,*) (dleaf(jveg),jveg=1,nveg)
      print *, 'dleaf',(dleaf(jveg),jveg=1,nveg)
      read(8,*) (vcmax(jveg),jveg=1,nveg)
c      read(8,*) (ejmax(jveg),jveg=1,nveg)
      print *, 'vcmax',(vcmax(jveg),jveg=1,nveg)
      do jveg=1,nveg
        ejmax(jveg)=2.*vcmax(jveg)
      enddo
      print *, 'ejmax',(ejmax(jveg),jveg=1,nveg)
      read(8,*) (hc(jveg),jveg=1,nveg)
      print *, 'hc',(hc(jveg),jveg=1,nveg)
      read(8,*) (xfang(jveg),jveg=1,nveg)
      print *, 'xfang',(xfang(jveg),jveg=1,nveg)
      read(8,*) (rp20(jveg),jveg=1,nveg)
      print *, 'rp20',(rp20(jveg),jveg=1,nveg)
      read(8,*) (rpcoef(jveg),jveg=1,nveg)
      print *, 'rpcoef',(rpcoef(jveg),jveg=1,nveg)
      read(8,*) (rs20(jveg),jveg=1,nveg)
      print *, 'rprs20',(rs20(jveg),jveg=1,nveg)
      read(8,*) (shelrb(jveg),jveg=1,nveg)
      print *, 'shelrb',(shelrb(jveg),jveg=1,nveg)
      read(8,*) (frac4(jveg),jveg=1,nveg)
      read(8,*) (tminvj(jveg),jveg=1,nveg)
      print *, 'tminvj',(tminvj(jveg),jveg=1,nveg)
      read(8,*) (tmaxvj(jveg),jveg=1,nveg)
      print *, 'tmaxvj',(tmaxvj(jveg),jveg=1,nveg)
      read(8,*) (vbeta(jveg),jveg=1,nveg)
      print *, 'vbeta',(vbeta(jveg),jveg=1,nveg)
      read(8,*) (tcplant(jveg,1),jveg=1,nveg)
      read(8,*) (tcplant(jveg,2),jveg=1,nveg)
      read(8,*) (tcplant(jveg,3),jveg=1,nveg)
      print *,'cplant 1',(tcplant(jveg,1),jveg=1,nveg)
      print *,'cplant 2',(tcplant(jveg,2),jveg=1,nveg)
      print *,'cplant 3',(tcplant(jveg,3),jveg=1,nveg)
      read(8,*) (tcsoil(jveg,1),jveg=1,nveg)
      read(8,*) (tcsoil(jveg,2),jveg=1,nveg)
      print *,'csoil 1',(tcsoil(jveg,1),jveg=1,nveg)
      print *,'csoil 2',(tcsoil(jveg,2),jveg=1,nveg)

      jyear=kdate/10000
      jmonth=(kdate-jyear*10000)/100

!      if( jmonth .eq. 1 .and. jyear .eq. 1994) then
      if( jmonth .eq. 1 .and. jyear .eq. inyear_carb) then
      do ip=1,ipland  ! all land points in this nsib=2 loop
      iq=iperm(ip)
       cplant(iq,1)=tcplant(ivegt(iq),1) 
       cplant(iq,2)=tcplant(ivegt(iq),2) 
       cplant(iq,3)=tcplant(ivegt(iq),3) 
      
       csoil(iq,1)=tcsoil(ivegt(iq),1) 
       csoil(iq,2)=tcsoil(ivegt(iq),2) 
       cansto(iq) = 0.
      enddo !ip
      endif ! jmonth
       print *,'cbmrdn',cplant(13419,1),cplant(13419,2),
     &       cplant(13419,3),csoil(13419,1),csoil(13419,2)
!      call comskp(8)
!      read(8,*) (froot(j),j=1,ms)
      call comskp(8)
      read(8,*) (ratecp(j),j=1,ncp)
      print *, 'ratecp',(ratecp(j),j=1,ncp)
      call comskp(8)
      read(8,*) (ratecs(j),j=1,ncs)
      close(8)
      print *, 'ratecs',(ratecs(j),j=1,ncs)
c
      end

      SUBROUTINE COMSKP(IUNIT)
C MRR, 5-AUG-83
C SKIPS COMMENT LINES IN CONTROL DATA FILE, STOPS IF EOF FOUND
      CHARACTER*1 COM
1     READ(IUNIT,100,END=2) COM
100   FORMAT(A1)
      IF(COM.EQ.'C'.or.com.eq.'c') GOTO 1
      BACKSPACE IUNIT
      RETURN
2     STOP 'CONTROL FILE EOF'
      END


      subroutine rlaiday

      include 'newmpar.h'
      include 'nsibd.h'
      include 'dates.h'     !  kdate,ktime,timer,mtimer
      include 'parm.h'      ! id,jd
      include 'permsurf.h'  ! iperm etc
      include 'vegpar.h'  
      include 'mpif.h'
      integer, parameter, dimension(0:13) :: mdays =
     &     (/ 31, 31,28,31,30,31,30,31,31,30,31,30,31, 31 /)
      integer iyr, imo, iday, ijd
      save iyr,imo,iday

      ijd = 9061
      print *,'rlaiday  called'
      iyr=kdate/10000
      imo=(kdate-10000*iyr)/100
      iday=kdate-10000*iyr-100*imo  +mtimer/(60*24)
      do while (iday>mdays(imo))
       iday=iday-mdays(imo)
       imo=imo+1
       if(imo>12)then
         imo=1               ! no leap years for the moment
         iyr=iyr+1
       endif
       enddo
!       print *,'iday,imo,iyr',iday,imo,iyr

       if(iday<mdays(imo)/2)then  ! 1st half of month
        rat1=(mdays(imo)-2.*iday)/(mdays(imo)+mdays(imo-1))
        rat2=(2.*iday+mdays(imo-1))/(mdays(imo)+mdays(imo-1))
        rlai=rat1*rlai123(:,1)+rat2*rlai123(:,2) 
       else                             ! 2nd half of month
        rat1=(mdays(imo+1)+2.*mdays(imo)-2.*iday)/
     .                               (mdays(imo+1)+mdays(imo))
        rat2=(2.*iday-mdays(imo))/(mdays(imo+1)+mdays(imo))
        rlai=rat1*rlai123(:,2)+rat2*rlai123(:,3) 
       endif
       if(iday==mdays(imo)/2) rlai=rlai123(:,2)
       do iq=1,ifull
         rlai(iq)=min(7.,max(0.0011,rlai(iq)))
         if (isoilm(iq).eq.9) rlai(iq)=0.0011
         vlai(iq)=rlai(iq)
       enddo

       print *,'month_imo,iday',imo,iday,rlai(ijd),rlai123(ijd,1),
     .       rlai123(ijd,2),rlai123(ijd,3)
      return
      end

