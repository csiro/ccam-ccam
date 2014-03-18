      subroutine sflux(nalpha)
      
      ! CCAM interface for surface flux routines
      ! Includes standard land-surface scheme,
      ! prescribed SSTs and sea-ice, CABLE
      ! interface, urban interface and MLO
      ! interface
      
      ! nsib=3              Standard land-surface scheme with SIB and Gratz data
      ! nsib=5              Standard land-surface scheme with MODIS data
      ! nsib=6              CABLE land-surface scheme with CABLE diagnostics
      ! nsib=7              CABLE land-surface scheme with CCAM diagnostics
      ! nmlo=0              Prescriped SSTs and sea-ice with JLM skin enhancement
      ! nmlo>0 and mlo<=9   KPP ocean mixing
      ! nmlo>9              Use external PCOM ocean model
      ! nurban>0            Use urban scheme
      
      use arrays_m                       ! Atmosphere dyamics prognostic arrays
      use ateb                           ! Urban
      use cable_ccam, only : sib4        ! CABLE interface
      use cc_mpi                         ! CC MPI routines
      use diag_m                         ! Diagnostic routines
      use extraout_m                     ! Additional diagnostics
      use gdrag_m                        ! Gravity wave drag
      use liqwpar_m                      ! Cloud water mixing ratios
      use map_m                          ! Grid map arrays
      use mlo                            ! Ocean physics and prognostic arrays
      use mlodynamics                    ! Ocean dynamics routines
      use morepbl_m                      ! Additional boundary layer diagnostics
      use nharrs_m                       ! Non-hydrostatic atmosphere arrays
      use nsibd_m                        ! Land-surface arrays
      use pbl_m                          ! Boundary layer arrays
      use permsurf_m                     ! Fixed surface arrays
      use prec_m                         ! Precipitation
      use river                          ! River routing
      use savuvt_m                       ! Saved dynamic arrays
      use screen_m                       ! Screen level diagnostics
      use sigs_m                         ! Atmosphere sigma levels
      use soil_m                         ! Soil and surface data
      use soilsnow_m                     ! Soil, snow and surface data
      use vecsuv_m                       ! Map to cartesian coordinates
      use vvel_m                         ! Additional vertical velocity
      use work2_m                        ! Diagnostic arrays
      use work3_m                        ! Mk3 land-surface diagnostic arrays
      use xyzinfo_m                      ! Grid coordinate arrays
      
      implicit none
      
      include 'newmpar.h'                ! Grid parameters
      include 'const_phys.h'             ! Physical constants
      include 'establ.h'                 ! Liquid saturation function
      include 'parm.h'                   ! Model configuration
      include 'parmgeom.h'               ! Coordinate data
      include 'parmsurf.h'               ! Surface parameters
      include 'soilv.h'                  ! Soil parameters
      include 'trcom2.h'                 ! Station data

      integer iq,k,it,ip,iqmin1,iqmax1,iqmin2,iqmax2
      integer, intent(in) :: nalpha
      real ri_max,zologbgin,ztv,z1onzt,chnsea
      real srcp,dtsoil,afrootpan,es,constz,drst
      real xx,consea,afroot,fm,con,dtsol,daf
      real con1,den,dden,dfm,root,denma,denha
      real conh,conw,zminlog,ri_ice,zoice,zologice
      real epotice,qtgnet,qtair,eg1,eg2,deg,b1
      real gbot,deltat,esatf,zobg,zologbg,zologx
      real afland,aftlandg,fhbg,rootbg,denhabg
      real thnew,thgnew,thnewa,qtgair,aftland
      real thgnewa,ri_tmp,fh_tmp,taftfhmin
      real taftfhmax,taftfhgmin,taftfhgmax,factchice
      real, dimension(:), allocatable, save :: taftfh,taftfhg
      real, dimension(:), allocatable, save :: plens
      real taftfhg_temp(ifull)
      real vmag(ifull),charnck(ifull)
      real zonx(ifull),zony(ifull),zonz(ifull),costh(ifull)
      real sinth(ifull),uzon(ifull),vmer(ifull),azmin(ifull)
      real uav(ifull),vav(ifull)
      real, dimension(ifull) :: neta,oldrunoff,newrunoff,rid,fhd
      real, dimension(ifull) :: fgf,rgg,fev,af,dirad,dfgdt,factch
      real, dimension(ifull) :: degdt,cie,aft,fh,ri,gamm,smixr,rho
      real, dimension(ifull) :: dumsg,dumr,dumx,dums,dumw
      logical, dimension(ifull) :: duml

      integer, parameter :: nblend=0  ! 0 for original non-blended, 1 for blended af
      integer, parameter :: ntss_sh=0 ! 0 for original, 3 for **3, 4 for **4
!     integer, parameter :: nplens=0  ! 0 to turn off plens, 10 (e.g.) is on
      integer, parameter :: ntest=0   ! ntest= 0 for diags off; ntest= 1 for diags on
!     integer, parameter :: ntaft=3   ! 0 for original, 3 nowadays
!                   1 & 2 tafthf constrained by prior values with 2 faster
!                   3 uses measure of prior tgf in calc. fh
!     parameter (newztsea=1)   ! 0 for original, 1 for correct zt over sea
!     vmag introduced Mar '05 as vmod was being used in ri
!     From 11/8/98 runoff() is accumulated & zeroed with precip
!     Now using tggsn(,1) for the tice calculations
c     with leads option via fracice (using tgg1 and tggsn1)
c     now in parm.h with zero as default and used in rdnsib/rdnscam
c     cp specific heat at constant pressure joule/kgm/deg
      real, parameter :: bprm=5.,cms=5.,chs=2.6,vkar=.4
      real, parameter :: d3=2.5
      real, parameter :: cgsoil=1000.,gksoil=.300e-6,rhog=1600.
      real, parameter :: d1land=.03
      real, parameter :: fmroot=.57735     ! was .4 till 7 Feb 1996

#include "log.h"

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

      if (ktau==1.and.nmaxpr==1) then
        write(6,*) "myid,land,sea ",myid,count(land),
     &    ifull-count(land)
      end if

      if (.not.allocated(plens).and.nmlo==0) then
        allocate(plens(ifull))
        plens=0.
      end if
      if (.not.allocated(taftfh).and.(nsib==3.or.nsib==5)) then
        allocate(taftfh(ifull))
        allocate(taftfhg(ifull))
        taftfh(:)=.05        ! just a diag default for sea points
        taftfhg(:)=7.e-4     ! just a diag default for sea points
      end if
      
      ri_max=(1./fmroot -1.)/bprm    ! i.e. .14641
      zologbgin=log(zmin/zobgin)     ! pre-calculated for all except snow points
      ztv=exp(vkar/sqrt(chn10))/10.  ! proper inverse of ztsea
      z1onzt=300.*rdry*(1.-sig(1))*ztv/grav
      chnsea=(vkar/log(z1onzt))**2   ! should give .00085 for csiro9
      oldrunoff(:)=runoff(:)
      zo=999.        ! dummy values
      factch=999.    ! dummy values
      taux=0.        ! dummy values
      tauy=0.        ! dummy values
      gamm=3.471e+05 ! dummy values

      if (diag.or.ntest==1) then
        if (mydiag) then
         if (land(idjd)) then
          write(6,*) 'entering sflux ktau,nsib,ivegt,isoilm,land '
     .         ,ktau,nsib,ivegt(idjd),isoilm(idjd),land(idjd)
          write(6,*) 'idjd,id,jd,slwa,sgsave ',
     .           idjd,id,jd,slwa(idjd),sgsave(idjd)
          write(6,*) 'snowd,sicedep,condx ',
     .           snowd(idjd),sicedep(idjd),condx(idjd)
          write(6,*) 't1,tss ',t(idjd,1),tss(idjd)
          write(6,*) 'wb ',(wb(idjd,k),k=1,ms)
          write(6,*) 'tgg ',(tgg(idjd,k),k=1,ms)
	   endif
        end if
        call maxmin(t,' t',ktau,1.,kl)
      endif

c     using av_vmod (1. for no time averaging)
!      *****  check next comment
!       sflux called at beginning of time loop, hence savu, savv

      azmin=(bet(1)*t(1:ifull,1)+phi_nh(:,1))/grav
      srcp =sig(1)**(rdry/cp)
      ga(:)=0.              !  for ocean points in ga_ave diagnostic
      theta(:)=t(1:ifull,1)/srcp
      rho(:)=ps(1:ifull)/(rdry*tss(:))
      uav(:)=av_vmod*u(1:ifull,1)+(1.-av_vmod)*savu(:,1)   
      vav(:)=av_vmod*v(1:ifull,1)+(1.-av_vmod)*savv(:,1)  
      vmod(:)=sqrt(uav(:)**2+vav(:)**2)  ! i.e. vmod for tss_sh
      vmag(:)=max( vmod(:) , vmodmin)    ! vmag used to calculate ri
      if(ntsur/=7) vmod(:)=vmag(:)       ! gives usual way

      !--------------------------------------------------------------
      START_LOG(sfluxwater)
      if (nmlo==0) then                                                 ! sea
       if(ntest==2.and.mydiag)write(6,*) 'before sea loop'              ! sea
!       from June '03 use basic sea temp from tgg1 (so leads is sensible)      
!       all sea points in this loop; also open water of leads           ! sea
       if(charnock>0.)then                                              ! sea
         charnck(:)=charnock                                            ! sea
       elseif(charnock<-1.)then            ! zo like Moon (2004)        ! sea
         charnck(:)=max(.0000386*u10(:),.000085*u10(:)-.00058)          ! sea
       else                                                             ! sea
         charnck(:)=.008+3.e-4*(u10(:)-9.)**2/ ! like Makin (2002)      ! sea
     &             (1.+(.006+.00008*u10(:))*u10(:)**2)                  ! sea
       endif                                                            ! sea
       do iq=1,ifull                                                    ! sea
        if(.not.land(iq))then                                           ! sea
         wetfac(iq)=1.                                                  ! sea
!        tpan holds effective sea for this loop                         ! sea
         if(ntss_sh==0)then                                             ! sea
          dtsol=.01*sgsave(iq)/(1.+.25*vmod(iq)**2)   ! solar heating   ! sea
          tpan(iq)=tgg(iq,1)+tss_sh*min(dtsol,8.)     ! of ssts         ! sea
         elseif(ntss_sh==3)then                                         ! sea
          dtsol=tss_sh*.01*sgsave(iq)/                                  ! sea
     .                 (1.+.035*vmod(iq)**3)          ! solar heating   ! sea
          tpan(iq)=tgg(iq,1)+min(dtsol,8.)            ! of ssts         ! sea
         elseif(ntss_sh==4)then                                         ! sea
          dtsol=tss_sh*.01*sgsave(iq)/                                  ! sea
     .                 (1.+vmod(iq)**4/81.)           ! solar heating   ! sea
          tpan(iq)=tgg(iq,1)+min(dtsol,8.)            ! of ssts         ! sea
         elseif(ntss_sh==5)then                                         ! sea
           dtsol=0.             ! raw SST                               ! sea
           tpan(iq)=tgg(iq,1)   ! raw SST                               ! sea
         endif   ! (ntss_sh==0) .. else ..                              ! sea
         if(nplens/=0)then                                              ! sea
!         calculate running total (over last 24 h) of daily precip in mm  jlm
          plens(iq)=(1.-dt/86400.)*plens(iq)+condx(iq)  ! in mm/day     ! sea
!         scale so that nplens m/s wind for 1/2 hr reduces effect by 1/1.2
          plens(iq)=plens(iq)/(1.+vmod(iq)*dt*.2/                       ! sea
     .                     max(nplens*1800.,1.))      ! avoids Cray compiler bug
!         produce a cooling of 4 K for an effective plens of 10 mm/day  ! sea
          tpan(iq)=tpan(iq)-min(.4*plens(iq) , 6.)                      ! sea
         endif   !  (nplens/=0)                                         ! sea
        if(ntsea==1.and.condx(iq)>.1)tpan(iq)=t(iq,2)                   ! sea
        if(ntsea==2.and.condx(iq)>.1)tpan(iq)=t(iq,1)                   ! sea
        if(ntsea==3.and.condx(iq)>.1)tpan(iq)=.5*(t(iq,2)+tgg(iq,1))    ! sea
        if(ntsea==4.and.condx(iq)>.1)tpan(iq)=.5*(t(iq,1)+tgg(iq,1))    ! sea
        endif  ! (.not.land(iq))                                        ! sea
       enddo   ! iq loop                                                ! sea
                                                                        ! sea
!      here calculate fluxes for sea point, and nominal pan points	    ! sea
       afrootpan=vkar/log(zmin/panzo)                                   ! sea
       do iq=1,ifull                                                    ! sea
c       drag coefficients  for momentum           cduv                  ! sea
c       for heat and moisture  cdtq                                     ! sea
        es = establ(tpan(iq))                                           ! sea
        constz=ps(iq)-es                                                ! sea
        qsttg(iq)= .98*.622*es/constz  ! with Zeng 1998 for sea water   ! sea
        drst=qsttg(iq)*ps(iq)*hlars/(constz*tpan(iq)**2)                ! sea
        xx=grav*zmin*(1.-tpan(iq)*srcp/t(iq,1))                         ! sea
        ri(iq)=min(xx/vmag(iq)**2 , ri_max)                             ! sea
c       this is in-line ocenzo using latest coefficient, i.e. .018      ! sea
        consea=vmod(iq)*charnck(iq)/grav  ! usually charnock=.018       ! sea
        if(land(iq))then                                                ! sea
          zo(iq)=panzo                                                  ! sea
          af(iq)=afrootpan**2                                           ! sea
        else                                                            ! sea
         if(charnock<0.)then  ! Moon (2004) over sea                    ! sea
          zo(iq)=charnck(iq)                                            ! sea
          afroot=vkar/log(zmin/zo(iq))                                  ! sea
          af(iq)=afroot**2                                              ! sea
         else            ! usual charnock method over sea               ! sea
          zo(iq)=.001    ! .0005 better first guess                     ! sea
          if(ri(iq)>0.)then             ! stable sea points             ! sea
            fm=vmod(iq) /(1.+bprm*ri(iq))**2 ! N.B. this is vmod*fm     ! sea
            con=consea*fm                                               ! sea
            do it=1,3                                                   ! sea
             afroot=vkar/log(zmin/zo(iq))                               ! sea
             af(iq)=afroot**2                                           ! sea
             daf=2.*af(iq)*afroot/(vkar*zo(iq))                         ! sea
             zo(iq)=max(1.5e-5,zo(iq)-(zo(iq)-con*af(iq))/              ! sea
     &                 (1.-con*daf))                                    ! sea
             zo(iq)=min(zo(iq),9.) ! JLM fix                            ! sea
            enddo    ! it=1,3                                           ! sea
            afroot=vkar/log(zmin/zo(iq))                                ! sea
            af(iq)=afroot**2                                            ! sea
          else                        ! unstable sea points             ! sea
            do it=1,3                                                   ! sea
             afroot=vkar/log(zmin/zo(iq))                               ! sea
             af(iq)=afroot**2                                           ! sea
             daf=2.*af(iq)*afroot/(vkar*zo(iq))                         ! sea
             con1=cms*2.*bprm*sqrt(-ri(iq)*zmin/zo(iq))                 ! sea
             den=1.+af(iq)*con1                                         ! sea
             dden=con1*(daf-.5*af(iq)/zo(iq))                           ! sea
             fm=vmod(iq)-vmod(iq)*2.*bprm *ri(iq)/den                   ! sea
             dfm=2.*bprm*ri(iq)*dden/den**2                             ! sea
             zo(iq)=max(1.5e-5,zo(iq)-(zo(iq)-consea*af(iq)*fm)/        ! sea
     .                        (1.-consea*(daf*fm+af(iq)*dfm)))          ! sea
             zo(iq)=min(zo(iq),6.) ! JLM fix                            ! sea
            enddo  ! it=1,3                                             ! sea
          endif    ! (xx>0.) .. else..                                  ! sea
         endif     ! (charnock<-1.) .. else ..                          ! sea
        endif      ! (land(iq)) .. else ..                              ! sea
        aft(iq)=chnsea                                                  ! sea
       enddo  ! iq loop                                                 ! sea
                                                                        ! sea
       if(newztsea==0)then ! 0 for original, 1 for different zt over sea! sea
c        enhanced formula used in Feb '92 Air-Sea conference follows:   ! sea
c        factch=sqrt(zo*exp(vkar*vkar/(chnsea*log(zmin/zo)))/zmin)      ! sea
         where (.not.land)                                              ! sea
           factch(:)=1. ! factch is sqrt(zo/zt) only for use in unstable fh
         end where                                                      ! sea
       else                                                             ! sea
         where (.not.land)                                              ! sea
           factch(:)=sqrt(zo(:)*ztv) ! for use in unstable fh           ! sea
         end where                                                      ! sea
       endif  ! (newztsea==0)                                           ! sea
                                                                        ! sea
       do iq=1,ifull ! done for all points; overwritten later for land  ! sea
c       Having settled on zo & af now do actual fh and fm calcs         ! sea
        if(ri(iq)>0.)then                                               ! sea
          fm=vmod(iq)/(1.+bprm*ri(iq))**2  ! no zo contrib for stable   ! sea
          fh(iq)=fm                                                     ! sea
        else        ! ri is -ve                                         ! sea
          root=sqrt(-ri(iq)*zmin/zo(iq))                                ! sea
c         First do momentum                                             ! sea
          denma=1.+cms*2.*bprm*af(iq)*root                              ! sea
          fm=vmod(iq)-vmod(iq)*2.*bprm *ri(iq)/denma                    ! sea
c         n.b. fm denotes ustar**2/(vmod(iq)*af)                        ! sea
c         Now heat ; allow for smaller zo via aft and factch            ! sea
c         N.B. for newztsea=1, zo contrib cancels in factch*root,       ! sea
c         so eg (& epan) and fg  (also aft) then indept of zo           ! sea
          denha=1.+chs*2.*bprm*factch(iq)*aft(iq)*root                  ! sea
          fh(iq)=vmod(iq)-vmod(iq)*2.*bprm *ri(iq)/denha                ! sea
        endif                                                           ! sea
                                                                        ! sea
        conh=rho(iq)*aft(iq)*cp                                         ! sea
        conw=rho(iq)*aft(iq)*hl                                         ! sea
        fg(iq)=conh*fh(iq)*(tpan(iq)-theta(iq))                         ! sea
        eg(iq)=conw*fh(iq)*(qsttg(iq)-qg(iq,1))                         ! sea
        rnet(iq)=sgsave(iq)-rgsave(iq)-stefbo*tpan(iq)**4               ! sea
        zoh(iq)=zo(iq)/(factch(iq)*factch(iq))                          ! sea
        zoq(iq)=zoh(iq)                                                 ! sea
c       cduv is now drag coeff *vmod                                    ! sea
        cduv(iq) =af(iq)*fm                                             ! sea
        cdtq(iq) =aft(iq)*fh(iq)                                        ! sea
        ustar(iq) = sqrt(vmod(iq)*cduv(iq))                             ! sea
c       Surface stresses taux, tauy: diagnostic only - unstaggered now  ! sea
        taux(iq)=rho(iq)*cduv(iq)*u(iq,1)                               ! sea
        tauy(iq)=rho(iq)*cduv(iq)*v(iq,1)                               ! sea
        ! note that iq==idjd  can only be true on the correct processor ! sea
        if(ntest==1.and.iq==idjd.and.mydiag)then                        ! sea
          write(6,*) 'in sea-type loop for iq,idjd: ',iq,idjd           ! sea
          write(6,*) 'zmin,zo,factch ',zmin,zo(iq),factch(iq)           ! sea
          write(6,*) 'ri,ustar,es ',ri(iq),ustar(iq),es                 ! sea
          write(6,*) 'af,aft ',af(iq),aft(iq)                           ! sea
          write(6,*) 'tpan,tss,theta ',tpan(iq),tss(iq),theta(iq)       ! sea
          write(6,*) 'chnsea,rho,t1 ',chnsea,rho(iq),t(iq,1)            ! sea
          write(6,*) 'fm,fh,conh ',fm,fh(iq),conh                       ! sea
          write(6,*) 'vmod,cduv,fg ',vmod(iq),cduv(iq),fg(iq)           ! sea
        endif                                                           ! sea
       enddo     ! iq loop                                              ! sea
       epot(:) = eg(:)                                                  ! sea
       epan(:) = eg(:)                                                  ! sea
c      section to update pan temperatures                               ! sea
       do iq=1,ifull                                                    ! sea
        if(land(iq))then                                                ! sea
          rgg(iq)=5.67e-8*tpan(iq)**4                                   ! sea
!         assume gflux = 0                                              ! sea
!         note pan depth=.254 m, spec heat water=4186 joule/kg K        ! sea
!         and change in heat supplied=spec_heatxmassxdelta_T            ! sea
          ga(iq)=-slwa(iq)-rgg(iq)-panfg*fg(iq)                         ! sea
          tpan(iq)=tpan(iq)+ga(iq)*dt/(4186.*.254*1000.)                ! sea
        else                                                            ! sea
          sno(iq)=sno(iq)+conds(iq)                                     ! sea
        endif  ! (land(iq))                                             ! sea
       enddo   ! iq loop                                                ! sea
                                                                        ! sea
       if(nmaxpr==1.and.mydiag)then                                     ! sea
        iq=idjd                                                         ! sea
        write (6,"('after sea loop fg,tpan,epan,ri,fh,vmod',            ! sea
     &    9f9.4)") fg(idjd),tpan(idjd),epan(idjd),                      ! sea
     &             ri(idjd),fh(idjd),vmod(idjd)                         ! sea
        write (6,"('u10,ustar,charnck,zo,cd',                           ! sea
     &    3f9.4,2f9.6)") u10(idjd),ustar(idjd),charnck(idjd),zo(idjd),  ! sea
     &    cduv(idjd)/vmod(idjd)                                         ! sea
        if(ri(iq)>0.)then                                               ! sea
          fm=vmod(iq)/(1.+bprm*ri(iq))**2                               ! sea
        else                                                            ! sea
          root=sqrt(-ri(iq)*zmin/zo(iq))                                ! sea
          denma=1.+cms*2.*bprm*af(iq)*root                              ! sea
          fm=vmod(iq)-vmod(iq)*2.*bprm *ri(iq)/denma                    ! sea
          denha=1.+chs*2.*bprm*factch(iq)*aft(iq)*root                  ! sea
        endif                                                           ! sea
        write (6,"('after sea loop af,aft,factch,root,denha,denma,fm',  ! sea
     &    9f9.4)") af(iq),aft(iq),factch(iq),root,denha,denma,fm        ! sea
       endif                                                            ! sea
       zminlog=log(zmin)                                                ! sea

       fgf=0.
       fev=0.
       do iq=1,ifull                                                    ! sice
        if(sicedep(iq)>0.)then                                          ! sice
!       non-leads for sea ice points                                    ! sice
!       N.B. tggsn( ,1) holds tice                                      ! sice
        es = establ(tggsn(iq,1))                                        ! sice
        constz=ps(iq)-es                                                ! sice
        qsttg(iq)= .622*es/constz                                       ! sice
        drst=qsttg(iq)*ps(iq)*hlars/(tggsn(iq,1)*tggsn(iq,1)*constz)    ! sice
        xx=grav*zmin*(1.-tggsn(iq,1)*srcp/t(iq,1))                      ! sice
        ri_ice=min(xx/vmag(iq)**2 , ri_max)                             ! sice
        !factchice=1. ! factch is sqrt(zo/zt) for use in unstable fh    ! sice
        factchice=sqrt(7.4) ! same as land from 27/4/99                 ! sice
        zoice=.001                                                      ! sice
        zologice=zminlog-log(zoice)   !   i.e. log(zmin/zo(iq))         ! sice
        af(iq)=(vkar/zologice)**2                                       ! sice
        aft(iq)=vkar**2/(zologice*zologice )                            ! sice
        wetfac(iq)=1+.008*(tggsn(iq,1)-273.16) ! 008*tggsn(iq,1)-1.18528! sice
                                                                        ! sice
c       now do fh and fm calcs for sice                                 ! sice
        if(ri_ice>0.)then                                               ! sice
          fm=vmod(iq)/(1.+bprm*ri_ice)**2                               ! sice
          fh(iq)=fm                                                     ! sice
        else                                                            ! sice
          root=sqrt(-ri_ice*zmin/zoice)                                 ! sice
c         First do momentum                                             ! sice
          denma=1.+cms*2.*bprm*af(iq)*root                              ! sice
          fm=vmod(iq)-vmod(iq)*2.*bprm *ri_ice/denma                    ! sice
c         n.b. fm denotes ustar**2/(vmod(iq)*af)                        ! sice
c         Now heat ; allow for smaller zo via aft and factch            ! sice
          denha=1.+chs*2.*bprm*factchice*aft(iq)*root                   ! sice
          fh(iq)=vmod(iq)-(2.*bprm *ri_ice)/denha                       ! sice
        endif                                                           ! sice
                                                                        ! sice
        conh=rho(iq)*aft(iq)*cp                                         ! sice
        conw=rho(iq)*aft(iq)*hl                                         ! sice
!       fgice & egice renamed as fgf and fev from Aug 05 to aid diags   ! sice	
        fgf(iq)=conh*fh(iq)*(tggsn(iq,1)-theta(iq))                     ! sice
        dfgdt(iq)=conh*fh(iq)                                           ! sice
        if(ntest==1.and.iq==idjd.and.mydiag)then                        ! sice
          write(6,*) 'in sice loop'                                     ! sice
          write(6,*) 'zmin,zo,wetfac ',zmin,zoice,wetfac(iq)            ! sice
          write(6,*) 'ri_ice,es ',ri_ice,es                             ! sice
          write(6,*) 'af,aft,ustar ',af(iq),aft(iq),ustar(iq)           ! sice
          write(6,*) 'chnsea,rho ',chnsea,rho(iq)                       ! sice
          write(6,*) 'fm,fh,conh ',fm,fh(iq),conh                       ! sice
        endif                                                           ! sice
                                                                        ! sice
        if(nalpha==1)then    ! beta scheme         sice here            ! sice
          epotice=conw*fh(iq)*(qsttg(iq)-qg(iq,1))                      ! sice
          fev(iq)=wetfac(iq)*epotice                                    ! sice
          degdt(iq)=wetfac(iq)*conw*fh(iq)*drst                         ! sice
        else                   ! alpha scheme                           ! sice
c         following trick reduces -ve evap (dew) to 1/10th value        ! sice
          qtgnet=qsttg(iq)*wetfac(iq) -qg(iq,1)                         ! sice
          qtgair=qsttg(iq)*wetfac(iq)-max(qtgnet,.1*qtgnet)             ! sice
          eg2=-conw*fh(iq)*qtgair                                       ! sice
          eg1=conw*fh(iq)*qsttg(iq)                                     ! sice
          fev(iq) =eg1*wetfac(iq) +eg2                                  ! sice
          epotice    = conw*fh(iq)*(qsttg(iq)-qg(iq,1))                 ! sice
          deg=wetfac(iq)*conw*fh(iq)*drst                               ! sice
c         following reduces degdt by factor of 10 for dew               ! sice
          degdt(iq)=.55*deg+sign(.45*deg,qtgnet)                        ! sice
        endif                                                           ! sice
                                                                        ! sice
c       section to update sea ice surface temperature;                  ! sice
c       specified sea-ice thickness                                     ! sice
c       over sea ice, set a minimum depth for this experiment of .1     ! sice
!       sicedep(iq) = max(sicedep(iq) , 0.1)  fixed in indata/nestin from Jan 06
c       no snow on the ice assumed for now                              ! sice
        gamm(iq) = 3.471e+05                                            ! sice
        cie(iq) = 2.04/sicedep(iq)                                      ! sice
        rgg(iq)=5.67e-8*tggsn(iq,1)**4                                  ! sice
!       gflux here is	flux from ice to water, +ve downwards           ! sice
        gflux(iq)=cie(iq)*(tggsn(iq,1)-271.2)                           ! sice
        ga(iq)=-slwa(iq)-rgg(iq)-fev(iq)-fgf(iq)-gflux(iq)              ! sice
        dirad(iq)=4.*5.67e-8*tggsn(iq,1)**3                             ! sice
        b1=dirad(iq)+degdt(iq)+dfgdt(iq)+cie(iq)                        ! sice
        gbot=(gamm(iq)/dt)+b1                                           ! sice
        deltat=ga(iq)/gbot                                              ! sice
        tggsn(iq,1)=tggsn(iq,1)+deltat                                  ! sice
        tggsn(iq,1)=min(tggsn(iq,1),271.2)   ! jlm fix Tue  05-30-2000  ! sice
        fgf(iq) =fgf(iq) +deltat*dfgdt(iq)                              ! sice
        fev(iq) =fev(iq) +deltat*degdt(iq)                              ! sice
        es = establ(tggsn(iq,1))                                        ! sice
        constz=ps(iq)-es                                                ! sice
        qsttg(iq)=.622*es/constz                                        ! sice
                                                                        ! sice
!       combine ice and leads contributions here                        ! sice
        eg(iq) =fracice(iq)*fev(iq) + (1.-fracice(iq))*eg(iq)           ! sice
        fg(iq) = fracice(iq)*fgf(iq)+ (1.-fracice(iq))*fg(iq)           ! sice
        ri(iq) =fracice(iq)*ri_ice + (1.-fracice(iq))*ri(iq)            ! sice
        zo(iq) =fracice(iq)*zoice  + (1.-fracice(iq))*zo(iq)            ! sice
        factch(iq)=fracice(iq)*factchice + (1.-fracice(iq))*factch(iq)  ! sice
        zoh(iq)=zo(iq)/(factch(iq)*factch(iq))                          ! sice
        zoq(iq)=zoh(iq)                                                 ! sice
        cduv(iq) =fracice(iq)*af(iq)*fm + (1.-fracice(iq))*cduv(iq)     ! sice
        cdtq(iq) =fracice(iq)*aft(iq)*fh(iq)+(1.-fracice(iq))*cdtq(iq)  ! sice
        ustar(iq) = sqrt(vmod(iq)*cduv(iq))                             ! sice
c       N.B. potential evaporation is now eg+eg2                        ! sice
        epot(iq) =fracice(iq)*epotice + (1.-fracice(iq))*epot(iq)       ! sice
        tss(iq) = fracice(iq)*tggsn(iq,1)+(1.-fracice(iq))*tpan(iq)     ! sice
        rnet(iq)=sgsave(iq)-rgsave(iq)-stefbo*tss(iq)**4                ! sice
c       Surface stresses taux, tauy: diagnostic only - unstag now       ! sice
        taux(iq)=rho(iq)*cduv(iq)*u(iq,1)                               ! sice
        tauy(iq)=rho(iq)*cduv(iq)*v(iq,1)                               ! sice
       endif  ! (sicedep(iq)>0.)                                        ! sice
       enddo       ! iq loop                                            ! sice
       where (.not.land)                                                ! sice
         snowd=0.                                                       ! sice
       end where                                                        ! sice
       if(mydiag.and.nmaxpr==1)then                                     ! sice
         write(6,*) 'after sice loop'                                   ! sice
         iq=idjd                                                        ! sice
         if(sicedep(iq)>0.)then                                         ! sice
           b1=dirad(iq)+degdt(iq)+dfgdt(iq)+cie(iq)                     ! sice
           gbot=(gamm(iq)/dt)+b1                                        ! sice
           deltat=ga(iq)/gbot                                           ! sice
           write(6,*) 'ri,vmag,vmod,cduv ',ri(iq),vmag(iq),vmod(iq)     ! sice
     &                                   ,cduv(iq)                      ! sice
           write(6,*) 'fh,tss,tpan,tggsn1 ',fh(iq),tss(iq),tpan(iq)     ! sice
     &                                  ,tggsn(iq,1)                    ! sice
           write(6,*) 'theta,t1,deltat ',theta(iq),t(iq,1),deltat       ! sice
           write(6,*) 'b1,ga,gbot,af,aft ',b1,ga(iq),gbot,af(iq)        ! sice
     &                                  ,aft(iq)                        ! sice
           write(6,*) 'fg,fgice,factch ',fg(iq),fgf(iq),factch(iq)      ! sice
           write(6,*) 'cie ',cie(iq)                                    ! sice
           write(6,*) 'eg,egice(fev),ustar ',eg(iq),fev(iq),ustar(iq)   ! sice
         end if ! sicedep(iq)>0.                                        ! sice
       endif    ! (mydiag.and.nmaxpr==1)                                ! sice

      elseif (abs(nmlo)>=1.and.abs(nmlo)<=9) then                       ! MLO
                                                                        ! MLO
        if (myid==0.and.nmaxpr==1) then                                 ! MLO
          write(6,*) "Before MLO mixing"                                ! MLO
        end if                                                          ! MLO
        if (abs(nmlo)==1) then                                          ! MLO
          ! Single column                                               ! MLO
          ! set free surface to zero when water is not conserved        ! MLO
          neta=0.                                                       ! MLO
          call mloimport(4,neta,0,0)                                    ! MLO
        end if                                                          ! MLO
                                                                        ! MLO
        ! Ocean mixing                                                  ! MLO
        dumsg=sgsave(:)/(1.-swrsave*albvisnir(:,1)-                     ! MLO
     &                 (1.-swrsave)*albvisnir(:,2))                     ! MLO
        dumr=-rgsave                                                    ! MLO
        dumx=condx/dt                                                   ! MLO
        dums=conds/dt                                                   ! MLO
        dumw=watbdy(1:ifull)/dt                                         ! MLO
        call mloeval(tss,zo,cduv,cdtq,fg,eg,wetfac,epot,epan,           ! MLO
     &               fracice,sicedep,snowd,dt,azmin,azmin,dumsg,        ! MLO
     &               dumr,dumx,dums,uav,vav,t(1:ifull,1),               ! MLO
     &               qg(1:ifull,1),ps,f,swrsave,fbeamvis,fbeamnir,      ! MLO
     &               dumw,0,.true.,oldu=oldu1,oldv=oldv1)               ! MLO
        call mloscrnout(tscrn,qgscrn,uscrn,u10,0)                       ! MLO
        call mloextra(0,zoh,azmin,0)                                    ! MLO
        call mloextra(3,zoq,azmin,0)                                    ! MLO
        call mloextra(1,taux,azmin,0)                                   ! MLO
        call mloextra(2,tauy,azmin,0)                                   ! MLO
        do k=1,ms                                                       ! MLO
          call mloexport(0,tgg(:,k),k,0)                                ! MLO
        end do                                                          ! MLO
        do k=1,3                                                        ! MLO
          call mloexpice(tggsn(:,k),k,0)                                ! MLO
        end do                                                          ! MLO
                                                                        ! MLO
        ! stuff to keep tpan over land working                          ! MLO
        rid=min(grav*zmin*(1.-tpan*srcp/t(1:ifull,1))/vmag**2,ri_max)   ! MLO
        where (rid>0.)                                                  ! MLO
          fhd=vmod/(1.+bprm*rid)**2                                     ! MLO
        elsewhere                                                       ! MLO
          fhd=vmod-vmod*2.*bprm*rid/(1.+chs*2.*bprm*sqrt(panzo*ztv)     ! MLO
     &       *chnsea*sqrt(-rid*zmin/panzo))                             ! MLO
        end where                                                       ! MLO
                                                                        ! MLO
        where(.not.land)                                                ! MLO
          watbdy(1:ifull)=0.                                            ! MLO
          snowd=snowd*1000.                                             ! MLO
          ga=0.                                                         ! MLO
          ustar=sqrt(sqrt(taux*taux+tauy*tauy)/rho)                     ! MLO
          tpan=tgg(:,1)                                                 ! MLO
          rnet=sgsave-rgsave-stefbo*tss**4                              ! MLO
          factch=sqrt(zo/zoh)                                           ! MLO
          sno=sno+conds                                                 ! MLO
          ! This cduv accounts for a moving surface                     ! MLO
          cduv=sqrt(ustar*ustar*cduv) ! cduv=cd*vmod                    ! MLO
          cdtq=cdtq*vmod                                                ! MLO
        elsewhere                                                       ! MLO
          fg=rho*chnsea*cp*fhd*(tpan-theta)                             ! MLO
          ga=sgsave-rgsave-5.67e-8*tpan**4-panfg*fg                     ! MLO
          tpan=tpan+ga*dt/(4186.*.254*1000.)                            ! MLO
        endwhere                                                        ! MLO
        do iq=1,ifull                                                   ! MLO
          if (.not.land(iq)) then                                       ! MLO
            esatf = establ(tss(iq))                                     ! MLO
            qsttg(iq)=.622*esatf/(ps(iq)-esatf)                         ! MLO
            rhscrn(iq)=100.*min(qgscrn(iq)/qsttg(iq),1.)                ! MLO
          end if                                                        ! MLO
        end do                                                          ! MLO
        if (myid==0.and.nmaxpr==1) then                                 ! MLO
          write(6,*) "After MLO mixing"                                 ! MLO
        end if                                                          ! MLO
                                                                        ! MLO
      else                                                              ! PCOM
        write(6,*) "ERROR: this option is for PCOM ocean model"         ! PCOM
        stop                                                            ! PCOM
      end if                                                            ! PCOM
      END_LOG(sfluxwater)
      !--------------------------------------------------------------      
      START_LOG(sfluxland)                                              ! land
      select case(nsib)                                                 ! land
        case(3,5)                                                       ! land
!cdir nodep
          do ip=1,ipland  ! all land points in this shared loop         ! land
c          fh itself was only used outside this loop in sib0 (jlm)      ! land
           iq=iperm(ip)                                                 ! land
           zobg=zobgin                                                  ! land
           es = establ(tss(iq))                                         ! land
           qsttg(iq)= .622*es/(ps(iq)-es)  ! prim for scrnout, bur recalc end sib3    
c          factch is sqrt(zo/zt) for land use in unstable fh            ! land
           factch(iq)=sqrt(7.4)                                         ! land
           if(snowd(iq)>0.)then                                         ! land
!            reduce zo over snow;                                       ! land
             zobg=max(zobgin -snowd(iq)*0.00976/12., 0.00024)           ! land
             zologbg=log(zmin/zobg)                                     ! land
!            following line is bit simpler than csiro9                  ! land
             zo(iq)=max(zolnd(iq) -.001*snowd(iq), .01)                 ! land
             zologx=log(zmin/zo(iq))                                    ! land
           else  ! land but not snow                                    ! land
             zo(iq)=zolnd(iq)                                           ! land
             zologbg=zologbgin                                          ! land
             zologx=zolog(iq)                                           ! land
           endif     ! (snowd(iq)>0.)                                   ! land
           if(nblend==1)then  ! blended zo for momentum                 ! land
!            note that Dorman & Sellers zo is already an average,       ! land
!            accounting for sigmf, so may not wish to further blend zo	! land
             afland=(vkar/((1.-sigmf(iq))*zologbg+sigmf(iq)*zologx))**2 ! land
           else    ! non-blended zo for momentum                        ! land
             afland=(vkar/zologx)**2                                    ! land
           endif   ! (nblend==1)                                        ! land
           aftland=vkar**2/( zologx * (2.+zologx) )                     ! land
           aft(iq)=aftland                                              ! land
           aftlandg=vkar**2/( zologbg * (2.+zologbg) )                  ! land
c          lgwd>0 enhances cduv (momentum) over orog under (stable & unst) condns
           af(iq)=afland+helo(iq)       ! jlm special gwd4b             ! land
                                                                        ! land
c          Having settled on zo (and thus af) now do actual fh and fm calcs 
           xx=grav*zmin*(1.-tss(iq)*srcp/t(iq,1))                       ! land
           ri(iq)=min(xx/vmag(iq)**2 , ri_max)                          ! land
           if(ri(iq)>0.)then                                            ! land
             fm=vmod(iq)/(1.+bprm*ri(iq))**2                            ! land
             fh(iq)=fm                                                  ! land
             fhbg=fh(iq)                                                ! land
           else                                                         ! land
             root=sqrt(-ri(iq)*zmin/zo(iq))  ! ignoring blending here   ! land
c            First do momentum                                          ! land
             denma=1.+cms*2.*bprm*af(iq)*root                           ! land
             fm=vmod(iq)-vmod(iq)*2.*bprm *ri(iq)/denma                 ! land
c            n.b. fm denotes ustar**2/(vmod(iq)*af)                     ! land
c            Now heat ; allow for smaller zo via aft and factch         ! land
             denha=1.+chs*2.*bprm*factch(iq)*aft(iq)*root               ! land
             fh(iq)=vmod(iq)-vmod(iq)*2.*bprm *ri(iq)/denha             ! land
             rootbg=sqrt(-ri(iq)*zmin/zobg)                             ! land
             denhabg=1.+chs*2.*bprm*factch(iq)*aftlandg*rootbg          ! land
             fhbg=vmod(iq)-vmod(iq)*2.*bprm *ri(iq)/denhabg             ! land
           endif                                                        ! land
           taftfhg_temp(iq)=aftlandg*fhbg ! uses fmroot above, for sib3 ! land
           taftfhg(iq)=aftlandg*fhbg ! value used for ntaft=3 (may need improving)
                                                                        ! land
           zoh(iq)=zo(iq)/(factch(iq)*factch(iq))                       ! land
           zoq(iq)=zoh(iq)                                              ! land
c          cduv is now drag coeff *vmod                                 ! land
           cduv(iq) =af(iq)*fm                                          ! land
           cdtq(iq) =aft(iq)*fh(iq)                                     ! land
           ustar(iq) = sqrt(vmod(iq)*cduv(iq))                          ! land
c          Surface stresses taux, tauy: diagnostic only - unstaggered now   
           taux(iq)=rho(iq)*cduv(iq)*u(iq,1)                            ! land
           tauy(iq)=rho(iq)*cduv(iq)*v(iq,1)                            ! land
           if(ntest==1.and.iq==idjd.and.mydiag)then                     ! land
             write(6,*) 'in main land loop'                             ! land
             write(6,*) 'zmin,zobg,zobgin,snowd ',zmin,zobg,zobgin      ! land
     &                                         ,snowd(iq)               ! land
             write(6,*) 'afland,aftland,zologbg ',afland,aftland        ! land
     &                                         ,zologbg                 ! land
             write(6,*) 'af,vmag,vmod,es ',af(iq),vmag(iq),vmod(iq),es  ! land
             write(6,*) 'tss,theta,t1 ',tss(iq),theta(iq),t(iq,1)       ! land
             write(6,*) 'aft,fm,fh,rho,conh ',aft(iq),fm,fh(iq),rho(iq) ! land
     &                                    ,conh                         ! land
             write(6,*) 'ri,vmod,cduv,fg ',ri(iq),vmod(iq),cduv(iq)     ! land
     &                                     ,fg(iq)                      ! land
           endif  ! (ntest==1.and.iq==idjd)                             ! land
          enddo     ! ip=1,ipland                                       ! land
                                                                        ! land
          if(ntaft==0.or.ktau==1)then                                   ! land
            do iq=1,ifull  ! will only use land values                  ! land
             if(land(iq))then                                           ! land
               taftfh(iq)=aft(iq)*fh(iq) ! uses fmroot above            ! land
               taftfhg(iq)=taftfhg_temp(iq)                             ! land
             endif                                                      ! land
            enddo                                                       ! land
          elseif(ntaft==1.)then                                         ! land
            do iq=1,ifull         ! will only use land values           ! land
             if(land(iq))then                                           ! land
               thnew=aft(iq)*fh(iq) ! uses fmroot above                 ! land
               thgnew=taftfhg_temp(iq)                                  ! land
               if(thnew>2.*taftfh(iq).or.thnew<.5*taftfh(iq))then       ! land
                 taftfh(iq)=.5*(thnew+taftfh(iq))                       ! land
               else                                                     ! land
                 taftfh(iq)=thnew                                       ! land
               endif                                                    ! land
               if(thgnew>2.*taftfhg(iq).or.thgnew<.5*taftfhg(iq))then   ! land
                 taftfhg(iq)=.5*(thgnew+taftfhg(iq))                    ! land
               else                                                     ! land
                 taftfhg(iq)=thgnew                                     ! land
               endif                                                    ! land
             endif                                                      ! land
            enddo                                                       ! land
          elseif(ntaft==2)then    ! preferred faster option             ! land
            do iq=1,ifull         ! will only use land values           ! land
             if(land(iq))then                                           ! land
               thnew=aft(iq)*fh(iq) ! uses fmroot above                 ! land
               thgnew=taftfhg_temp(iq)                                  ! land
               thnewa=min(thnew,                                        ! land
     .                max(2.*taftfh(iq),.5*(thnew+taftfh(iq))))         ! land
               taftfh(iq)=max(thnewa,                                   ! land
     .                min(.5*taftfh(iq),.5*(thnew+taftfh(iq))))         ! land
               thgnewa=min(thgnew,                                      ! land
     .                 max(2.*taftfhg(iq),.5*(thgnew+taftfhg(iq))))     ! land
               taftfhg(iq)=max(thgnewa,                                 ! land
     .                     min(.5*taftfhg(iq),.5*(thgnew+taftfhg(iq)))) ! land
             endif                                                      ! land
            enddo                                                       ! land
          elseif(ntaft==3)then                                          ! land
!           do vegetation calulation for taftfh	                        ! land
            do iq=1,ifull                                               ! land
             if(land(iq))then                                           ! land
               xx=grav*zmin*(1.-tgf(iq)*srcp/t(iq,1)) ! actually otgf   ! land
               ri_tmp=min(xx/vmag(iq)**2 , ri_max)                      ! land
               if(ri_tmp>0.)then                                        ! land
                 fh_tmp=vmod(iq)/(1.+bprm*ri_tmp)**2                    ! land
               else                                                     ! land
                 root=sqrt(-ri_tmp*zmin/zo(iq))  ! ignoring blending    ! land
c                Now heat ; allow for smaller zo via aft and factch     ! land
                 denha=1.+chs*2.*bprm*factch(iq)*aft(iq)*root           ! land
                 fh_tmp=vmod(iq)-vmod(iq)*2.*bprm *ri_tmp/denha         ! land
               endif                                                    ! land
               taftfh(iq)=aft(iq)*fh_tmp ! uses fmroot above, for sib3  ! land
             endif                                                      ! land
            enddo                                                       ! land
          endif  ! (ntaft==0.or.ktau==1)  .. else ..                    ! land
          if(ntest>0.and.mydiag)then                                    ! land
            write(6,*) 'before sib3 zo,zolnd,af ',zo(idjd),zolnd(idjd)  ! land
     &                                        ,af(idjd)                 ! land
            write(6,*) 'av_vmod,u,savu,v,savv',                         ! land
     &          av_vmod,u(idjd,1),savu(idjd,1),v(idjd,1),savv(idjd,1)   ! land
          endif                                                         ! land
          call sib3(nalpha,taftfh,taftfhg,aft,rho) ! for nsib=3, 5      ! land
          if(diag.or.ntest>0)then                                       ! land
            if (mydiag) write(6,*) 'before call scrnout'                ! land
            call maxmin(t,' t',ktau,1.,kl)                              ! land
          endif                                                         ! land
          if(ntsur/=5)then    ! ntsur=6 is default from Mar '05         ! land
c           preferred option to recalc cduv, ustar (gives better uscrn, u10)
            do iq=1,ifull                                               ! land
             afroot=vkar/log(zmin/zo(iq))! land formula is bit different above
             af(iq)=afroot**2+helo(iq)                                  ! land
             xx=grav*zmin*(1.-tss(iq)*srcp/t(iq,1))                     ! land
             ri(iq)=min(xx/vmag(iq)**2 , ri_max)                        ! land
             if(ri(iq)>0.)then                                          ! land
               fm=vmod(iq)/(1.+bprm*ri(iq))**2         ! Fm * vmod      ! land
             else                                                       ! land
               root=sqrt(-ri(iq)*zmin/zo(iq))                           ! land
               denma=1.+cms*2.*bprm*af(iq)*root                         ! land
               fm=vmod(iq)-vmod(iq)*2.*bprm *ri(iq)/denma   ! Fm * vmod ! land
c              n.b. fm denotes ustar**2/(vmod(iq)*af)                   ! land                  
             endif                                                      ! land
c            cduv is now drag coeff *vmod                               ! land      
             cduv(iq) =af(iq)*fm                       ! Cd * vmod      ! land
             ustar(iq) = sqrt(vmod(iq)*cduv(iq))                        ! land
c            Surface stresses taux, tauy: diagnostic only - unstaggered now   
             taux(iq)=rho(iq)*cduv(iq)*u(iq,1)                          ! land
             tauy(iq)=rho(iq)*cduv(iq)*v(iq,1)                          ! land
            enddo                                                       ! land
           endif  ! (ntsur==6)                                          ! land
           if(nproc==1.and.diag)then                                    ! land
              taftfhmin=1.e20                                           ! land
              taftfhmax=-100.                                           ! land
              taftfhgmin=1.e20                                          ! land
              taftfhgmax=-100.                                          ! land
              do iq=1,ifull                                             ! land
               if(taftfh(iq)<taftfhmin)then                             ! land
                 taftfhmin=taftfh(iq)           ! ~.0012                ! land
                 iqmin1=iq                                              ! land
               endif                                                    ! land
               if(taftfh(iq)>taftfhmax)then                             ! land
                 taftfhmax=taftfh(iq)           ! ~.13                  ! land
                 iqmax1=iq                                              ! land
               endif                                                    ! land
               if(taftfhg(iq)<taftfhgmin)then                           ! land
                 taftfhgmin=taftfhg(iq)         ! ~.0006                ! land
                 iqmin2=iq                                              ! land
               endif                                                    ! land
               if(taftfhg(iq)>taftfhgmax)then                           ! land
                 taftfhgmax=taftfhg(iq)         ! ~.004                 ! land
                 iqmax2=iq                                              ! land
               endif                                                    ! land
              enddo                                                     ! land
              write(6,*) 'taftfhmin,taftfhmax ',                        ! land
     &                 taftfhmin,iqmin1,taftfhmax,iqmax1                ! land
              write(6,*) 'taftfhgmin,taftfhgmax ',                      ! land
     &                 taftfhgmin,iqmin2,taftfhgmax,iqmax2              ! land
            endif  ! (nproc==1.and.diag)                                ! land
                                                                        ! land
           ! scrnout is the standard CCAM screen level diagnostics.     ! land
           ! However, it overwrites ocean points when using MLO.        ! land
           ! Therefore, we use scrnocn which can calculate screen       ! land
           ! level diagnostics for just ocean or just land points       ! land
           if (nmlo==0) then                                            ! land
             smixr=wetfac*qsttg+(1.-wetfac)*min(qsttg,qg(1:ifull,1))    ! land
             call scrnout(zo,ustar,factch,wetfac,smixr,                 ! land
     .              qgscrn,tscrn,uscrn,u10,rhscrn,af,aft,ri,vmod,       ! land
     .              bprm,cms,chs,chnsea,nalpha)                         ! land
           else                                                         ! land
             dumx=sqrt(u(1:ifull,1)*u(1:ifull,1)+                       ! land
     &                       v(1:ifull,1)*v(1:ifull,1))                 ! land
             duml=.not.land                                             ! land
             call scrnocn(ifull,qgscrn,tscrn,uscrn,u10,rhscrn,zo,zoh,   ! land
     &                  zoq,tss,t(1:ifull,1),qsttg,qg(1:ifull,1),       ! land
     &                  dumx,ps(1:ifull),duml,azmin,sig(1))             ! land
           end if                                                       ! land
        case(6,7)                                                       ! cable
         if (myid==0.and.nmaxpr==1) then                                ! cable
           write(6,*) "Before CABLE"                                    ! cable
         end if                                                         ! cable
         ! call cable                                                   ! cable
         call sib4                                                      ! cable
         ! update remaining diagnostic arrays                           ! cable
         do ip=1,ipland                                                 ! cable
           iq=iperm(ip)                                                 ! cable
           factch(iq)=sqrt(zo(iq)/zoh(iq))                              ! cable
           es = establ(tss(iq))                                         ! cable
           qsttg(iq)= .622*es/(ps(iq)-es)                               ! cable
           rhscrn(iq)=100.*qgscrn(iq)/qsttg(iq)                         ! cable
           rhscrn(iq)=min(max(rhscrn(iq),0.),100.)                      ! cable
           taux(iq)=rho(iq)*cduv(iq)*u(iq,1)                            ! cable
           tauy(iq)=rho(iq)*cduv(iq)*v(iq,1)                            ! cable
           sno(iq)=sno(iq)+conds(iq)                                    ! cable
         enddo   ! ip=1,ipland                                          ! cable
         if (nmlo==0) then                                              ! cable
           ! update ocean diagnostics                                   ! cable
           smixr=wetfac*qsttg+(1.-wetfac)*min(qsttg,qg(1:ifull,1))      ! cable
           dumx=sqrt(u(1:ifull,1)*u(1:ifull,1)+                         ! cable
     &                       v(1:ifull,1)*v(1:ifull,1))                 ! cable
           call scrnocn(ifull,qgscrn,tscrn,uscrn,u10,rhscrn,zo,zoh,     ! cable
     &                  zoq,tss,t(1:ifull,1),smixr,qg(1:ifull,1),       ! cable
     &                  dumx,ps(1:ifull),land,azmin,sig(1))             ! cable
         end if                                                         ! cable
         smixr=wetfac*qsttg+(1.-wetfac)*min(qsttg,qg(1:ifull,1))        ! cable
         ! The following patch overrides CABLE screen level diagnostics ! cable
         if (nsib==7) then                                              ! cable
           dumx=sqrt(u(1:ifull,1)*u(1:ifull,1)+                         ! cable
     &                       v(1:ifull,1)*v(1:ifull,1))                 ! cable
           duml=.not.land                                               ! cable
           call scrnocn(ifull,qgscrn,tscrn,uscrn,u10,rhscrn,zo,zoh,     ! cable
     &                zoq,tss,t(1:ifull,1),smixr,qg(1:ifull,1),         ! cable
     &                dumx,ps(1:ifull),duml,azmin,sig(1))               ! cable
         end if                                                         ! cable
         if (myid==0.and.nmaxpr==1) then                                ! cable
           write(6,*) "After CABLE"                                     ! cable
         end if                                                         ! cable
        case DEFAULT                                                    ! land
          write(6,*) "ERROR: Unknown land-use option nsib=",nsib        ! land
          stop                                                          ! land
      end select                                                        ! land
      END_LOG(sfluxland)                                                ! land
      !----------------------------------------------------------
      START_LOG(sfluxurban)                                             ! urban
      if (myid==0.and.nmaxpr==1) then                                   ! urban
        write(6,*) "Before urban"                                       ! urban
      end if                                                            ! urban
      if (nurban/=0) then                                               ! urban
         ! calculate zonal and meridonal winds                          ! urban
         zonx=                       -sin(rlat0*pi/180.)*y(:)           ! urban
         zony=sin(rlat0*pi/180.)*x(:)+cos(rlat0*pi/180.)*z(:)           ! urban
         zonz=-cos(rlat0*pi/180.)*y(:)                                  ! urban
         costh= (zonx*ax(1:ifull)+zony*ay(1:ifull)+zonz*az(1:ifull))    ! urban
     &          /sqrt( max(zonx**2+zony**2+zonz**2,1.e-7) )             ! urban
         sinth=-(zonx*bx(1:ifull)+zony*by(1:ifull)+zonz*bz(1:ifull))    ! urban
     &          /sqrt( max(zonx**2+zony**2+zonz**2,1.e-7) )             ! urban
         uzon= costh*uav-sinth*vav ! zonal wind                         ! urban
         vmer= sinth*uav+costh*vav ! meridonal wind                     ! urban
         newrunoff=runoff-oldrunoff ! new runoff since entering sflux   ! urban
         ! since ateb will blend non-urban and urban runoff, it is      ! urban
         ! easier to remove the new runoff and add it again after the   ! urban
         ! urban scheme has been updated                                ! urban
         runoff=oldrunoff        ! remove new runoff                    ! urban
         ! call aTEB                                                    ! urban
         dumsg=sgsave/(1.-swrsave*albvisnir(:,1)-                       ! urban
     &               (1.-swrsave)*albvisnir(:,2))                       ! urban
         dumr=-rgsave                                                   ! urban
         dumx=condx/dt                                                  ! urban
         dums=conds/dt                                                  ! urban
         call atebcalc(fg,eg,tss,wetfac,newrunoff,dt,azmin,dumsg,dumr   ! urban
     &               ,dumx,dums,rho,t(1:ifull,1),qg(1:ifull,1)          ! urban
     &               ,ps(1:ifull),uzon,vmer,vmodmin,0)                  ! urban
        runoff=runoff+newrunoff ! add new runoff after including urban  ! urban
        ! here we blend zo with the urban part                          ! urban
        call atebzo(zo,zoh,zoq,0)                                       ! urban
        factch=sqrt(zo/zoh)                                             ! urban
        ! calculate ustar                                               ! urban
        cduv=cduv/vmag                                                  ! urban
        cdtq=cdtq/vmag                                                  ! urban
        call atebcd(cduv,cdtq,0)                                        ! urban
        cduv=cduv*vmag                                                  ! urban
        cdtq=cdtq*vmag                                                  ! urban
        ustar=sqrt(vmod*cduv)                                           ! urban
        ! calculate screen level diagnostics                            ! urban
        call atebscrnout(tscrn,qgscrn,uscrn,u10,0)                      ! urban
        do ip=1,ipland ! assumes all urban points are land points       ! urban
          iq=iperm(ip)                                                  ! urban
          if (sigmu(iq)>0.) then                                        ! urban
            es = establ(tss(iq))                                        ! urban
            qsttg(iq)= .622*es/(ps(iq)-es)                              ! urban
            rhscrn(iq)=100.*qgscrn(iq)/qsttg(iq)                        ! urban
            rhscrn(iq)=min(max(rhscrn(iq),0.),100.)                     ! urban
            rnet(iq)=sgsave(iq)-rgsave(iq)-stefbo*tss(iq)**4            ! urban
            taux(iq)=rho(iq)*cduv(iq)*u(iq,1)                           ! urban
            tauy(iq)=rho(iq)*cduv(iq)*v(iq,1)                           ! urban
          end if                                                        ! urban
        end do                                                          ! urban
      end if                                                            ! urban
      if (myid==0.and.nmaxpr==1) then                                   ! urban
        write(6,*) "After urban"                                        ! urban
      end if                                                            ! urban
      END_LOG(sfluxurban)                                               ! urban
c ----------------------------------------------------------------------
      evap(:)=evap(:)+dt*eg(:)/hl !time integ value in mm (wrong for snow)

      ! Update runoff for river routing
      if (abs(nmlo)>=2) then
        newrunoff=runoff-oldrunoff
        salbdy(1:ifull)=salbdy(1:ifull)*watbdy(1:ifull)
     &        /max(watbdy(1:ifull)+newrunoff,1.E-10)
        watbdy(1:ifull)=watbdy(1:ifull)+newrunoff ! runoff in mm
      end if

c***  end of surface updating loop

      if(diag.or.ntest==1)then
         if ( mydiag ) then
           write(6,*) 'at end of sflux'
           write(6,*) 'eg,fg ',
     &               eg(idjd),fg(idjd)
           write(6,*) 'tscrn,cduv,zolnd ',
     &               tscrn(idjd),cduv(idjd),zolnd(idjd)
           write(6,*) 'snowd,sicedep ',
     &               snowd(idjd),sicedep(idjd)
           write(6,*) 'u1,v1,qg1 ',u(idjd,1),v(idjd,1),qg(idjd,1)
           write(6,*) 'w,w2,condx ',
     &               wb(idjd,1),wb(idjd,ms),condx(idjd)
           write(6,*) 't1,tss,tgg_2,tgg_ms ',
     &               t(idjd,1),tss(idjd),tgg(idjd,2),tgg(idjd,ms)
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine sib3(nalpha,taftfh,taftfhg,aft,rho)

      ! This is the standard land-surface scheme
      
      use arrays_m                     ! Atmosphere dyamics prognostic arrays
      use cc_mpi                       ! CC MPI routines
      use extraout_m                   ! Additional diagnostics
      use latlong_m                    ! Lat/lon coordinates
      use liqwpar_m                    ! Cloud water mixing ratios
      use morepbl_m                    ! Additional boundary layer diagnostics
      use nsibd_m                      ! Land-surface arrays
      use pbl_m                        ! Boundary layer arrays
      use permsurf_m                   ! Fixed surface arrays
      use screen_m                     ! Screen level diagnostics
      use sigs_m                       ! Atmosphere sigma levels
      use soil_m                       ! Soil and surface data
      use soilsnow_m                   ! Soil, snow and surface data
      use vegpar_m                     ! Vegetation arrays
      use work2_m                      ! Diagnostic arrays
      use work3_m                      ! Mk3 land-surface diagnostic arrays

      implicit none
      
      include 'newmpar.h'              ! Grid parameters
      include 'const_phys.h'           ! Physical constants
      include 'dates.h'                ! Date data
      include 'establ.h'               ! Liquid saturation function
      include 'parm.h'                 ! Model configuration
      include 'soilv.h'                ! Soil parameters
      include 'trcom2.h'               ! Station data

      integer nbarewet,nsigmf
      common/nsib/nbarewet,nsigmf      ! Land-surface options

      integer iq,k,iveg,layer,isoil,ip,icount
      integer, intent(in) :: nalpha
      real xxx,tgss,esattg,tgss2,fle,frac,conw_fh,qtgnet
      real qtgair,eg1,eg2,deg,egg_alph1,sstar,ff,rsi,den
      real wbav,f1,f2,f3,f4,esatf,qsatgf,beta,etr,betetrdt
      real prz,dirad1,dqg,devf,residp,sstgf,ewwwa,delta_t0
      real delta_t,deltat,es,tsoil
      real, dimension(ifull), intent(in) :: taftfh,taftfhg,rho
      real, dimension(ifull), intent(inout) :: aft
      real airr(ifull),cc(ifull),condxg(ifull),delta_tx(ifull),
     . evapfb1(ifull),evapfb2(ifull),evapfb3(ifull),evapfb4(ifull),
     . evapfb5(ifull),evapfb1a(ifull),evapfb2a(ifull),evapfb3a(ifull),
     . evapfb4a(ifull),evapfb5a(ifull),otgf(ifull),rmcmax(ifull),
     . tgfnew(ifull),evapfb(ifull)
      real dqsttg(ifull),tstom(ifull),cls(ifull),omc(ifull)
      real, dimension(ifull) :: ftsoil,rlai,srlai,res,tsigmf,fgf,egg
      real, dimension(ifull) :: evapxf,ewww,fgg,rdg,rgg,residf,fev
      real, dimension(ifull) :: extin,dirad,dfgdt,degdt

      integer, parameter :: ntest=0 ! ntest= 0 for diags off; ntest= 1 for diags on
!                                            2 for ewww diags      
!                           N.B. may need vsafe for correct diags
      integer, parameter :: itnmeth=5   ! 0 for original N_R iteration method
!     integer, parameter :: nsigmf=1    ! 0 for original tsigmf usage in sib3; prefer 1
!     integer, parameter :: nbarewet=2  ! 0 for original bare-soil-wetfac; 1 simplest; 2 for jlm
      integer, parameter :: nstomata=1  ! 0 for original; 1 for tropical C4
      integer, parameter :: newfgf=0    ! 0 for original; 1 with tscrn; 2 with aft too
      integer, parameter :: ndiag_arr=0 ! 0 off; 1 for diag arrays on
      integer, parameter :: neva=0      ! neva= 0 for diags off; neva= 1 for eva's diags on


!     fle(isoil,w)=(w-swilt(isoil))/(sfc(isoil)-swilt(isoil))           !0 Eva's
!     fle(isoil,w)= w/ssat(isoil)                                       !1 simplest bare
!     fle(isoil,w)= (w-frac*max( w,swilt(isoil) ))/                     !2 jlm sugg. bare
!    .               (ssat(isoil)*(1.-frac))                            !2 jlm sugg. bare
!     fle(isoil,w)=10.*((w-ssoil(isoil))/ssat(isoil)+.1))               ! an old one
!     fle(isoil,w)= w/sfc(isoil)                                        ! jlm for PIRCS
!     fle(isoil,w)=(w-frac*swilt(isoil))/(sfc(isoil)-frac*swilt(isoil)) ! jlm special
     
      do iq=1,ifull
       if(land(iq))then
         iveg=ivegt(iq)
c        evaluate seasonal impact and the snow depth
c        impact on the fractional vegetation cover
         tstom(iq)=298.
         if(iveg==6+31)tstom(iq)=302.
         if(iveg>=10.and.iveg<=21.and.
     &      abs(rlatt(iq)*180./pi)<25.)tstom(iq)=302.
           tsoil=min(tstom(iq), .5*(.3333*tgg(iq,2)+.6667*tgg(iq,3)
     &                          +.95*tgg(iq,4) + .05*tgg(iq,5)))
           ftsoil(iq)=max(0.,1.-.0016*(tstom(iq)-tsoil)**2)
c          which is same as:  ftsoil=max(0.,1.-.0016*(tstom-tsoil)**2)
c                             if( tsoil >= tstom ) ftsoil=1.
         endif ! (land)
      enddo

      if(nsib==3)then
        do iq=1,ifull
         if(land(iq))then
           iveg=ivegt(iq)
           rlai(iq)=max(.1,rlaim44(iveg)-slveg44(iveg)*(1.-ftsoil(iq)))
           srlai(iq)=rlai(iq)+rlais44(iveg)    ! nsib=3  leaf area index
           rsmin(iq) = rsunc44(iveg)/rlai(iq)  ! nsib=3  
           tsigmf(iq)=max(.001, sigmf(iq)-scveg44(iveg)*(1.-ftsoil(iq)))
           vlai(iq)=rlai(iq)
         else
           vlai(iq)=0.
         endif ! (land)
        enddo
      else     ! i.e. nsib=5
        where (land)
         rlai=max(.1,vlai)
         srlai=rlai                  ! nsib=5 leaf area index
         tsigmf=max(.001,sigmf)
        end where
      endif  !(nsib==3) .. else ..

      if(ktau==1)then
        if(mydiag)write(6,*) 'ipland,ipsice,ipsea in sflux: ',
     .                        ipland,ipsice,ipsea
        do iq=1,ifull     ! gives default over sea too
         tgf(iq)=t(iq,1)  ! was tss(iq)
         cansto(iq)=0.
        enddo 
        do iq=1,ifull
         if(land(iq))then  ! following gives agreement on restarts
           tscrn(iq)=theta(iq)  ! first guess, needed for newfgf=1
           if(nrungcm==3)then
             do layer=2,ms
              wb(iq,layer)=wb(iq,1)   ! w, w2 and wb all same initially 
             enddo
           endif  ! (nrungcm==3)
         endif  ! (land)
        enddo
        if ( mydiag ) then
          if ( land(idjd) ) then ! MJT bugfix
           iveg=ivegt(idjd)
           isoil = isoilm(idjd)
           tsoil=0.5*(0.3333*tgg(idjd,2)+0.6667*tgg(idjd,3)
     &             +0.95*tgg(idjd,4) +  0.05*tgg(idjd,5))
           write(6,*) 'nsib,iveg,isoil,nalpha,newfgf,nsigmf,tsigmf ',
     &           nsib,iveg,isoil,nalpha,newfgf,nsigmf,tsigmf(idjd)
           write(6,*) 'ftsoil,scveg44,sigmf ',
     &              ftsoil(idjd),scveg44(iveg),sigmf(idjd)
           write(6,*) 'swilt,sfc,wb1-6 ',
     &              swilt(isoil),sfc(isoil),(wb(idjd,k),k=1,ms)
           write(6,*) 'srlai,rsmin ',srlai(idjd),rsmin(idjd)
          endif
        endif
      endif           ! (ktau==1)

      if(ntest==1.and.mydiag) then
         iq=idjd
         iveg=ivegt(iq)
         write(6,*) 'in sib3a iq,iveg ',iq,iveg
         write(6,*) 'snowd,zo,zolnd,tstom ',
     &           snowd(iq),zo(iq),zolnd(iq),tstom(iq)
         write(6,*) 'in sib3b iq,idjd,iveg ',iq,idjd,iveg
         write(6,*) 'iveg,sigmf(iq),tsigmfa ',iveg,sigmf(iq),tsigmf(iq)
         tsoil=0.5*(0.3333*tgg(idjd,2)+0.6667*tgg(idjd,3)
     &             +0.95*tgg(idjd,4) +  0.05*tgg(idjd,5))
         write(6,*) 'rlaim44,tsoil,ftsoil ',
     &           rlaim44(iveg),tsoil,ftsoil(iq)
         write(6,*) 'scveg44,snowd,zo,zolnd,tstom ',
     &          scveg44(iveg),snowd(iq),zo(iq),zolnd(iq),tstom(iq)
         write(6,*) 'w2,rlai ',wb(iq,ms),rlai(iq)
       endif ! ntest
 
!cdir nodep
       do ip=1,ipland  
        iq=iperm(ip)
        tsigmf(iq)=(1.-snowd(iq)/
     &             (snowd(iq)+5.*100.*zo(iq)))*tsigmf(iq)
!      extin(iq)=exp(-0.6*max(1.,rlai(iq)))  ! good approx uses next 2 (jlm)
       xxx=.6*max(1.,rlai(iq))
       extin(iq)=1.-xxx/(1. +.5*xxx +xxx*xxx/12.) 
       if(ntest==1.and.iq==idjd.and.mydiag) then
         write(6,*) 'in sib3c ip,iq,idjd,iveg ',ip,iq,idjd,ivegt(iq)
         write(6,*) 'iveg,sigmf(iq),tsigmf ',ivegt(iq),sigmf(iq)
     &              ,tsigmf(iq)
         write(6,*) 'scveg44,snowd,zo,zolnd,tstom ',
     &          scveg44(iveg),snowd(iq),zo(iq),zolnd(iq),tstom(iq)
         write(6,*) 'alb,sgsave ',albvisnir(iq,1),sgsave(iq)
         write(6,*) 'w2,rlai,extin ',wb(iq,ms),rlai(iq),extin(iq)
       endif ! ntest
c      bare ground calculation
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

!     wetfac(:)=wetfac(:)*(1.-sigmu(:)) ! MJT delete urban

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
      if((ntest==1.or.diag).and.mydiag ) then
       if (land(idjd))then ! MJT bugfix
          iq=idjd
          write(6,*) 'epot,egg,tgg1,snowd ',
     &             epot(iq),egg(iq),tgg(iq,1),snowd(iq)
          isoil = isoilm(iq)
          conw_fh=rho(iq)*taftfhg(iq)*hl    ! taftfhg(iq)=aftlandg*fhbg
          qtgnet=  qsttg(iq)*wetfac(iq) -qg(iq,1)
          qtgair=  qsttg(iq)*wetfac(iq)-max(qtgnet,.1*qtgnet)
          eg2=   -conw_fh*qtgair
          eg1=    conw_fh*qsttg(iq)
!         degdw(iq)=eg1/ssat(isoil)                   ! not used in sib3
c         evaporation from the bare ground
          egg_alph1=wetfac(iq)*epot(iq)
          write(6,*) 'then iq,isoil,conw_fh,qsttg,qtgair ',
     .             iq,isoil,conw_fh,qsttg(iq),qtgair
          write(6,*) 'eg1,eg2,wetfac ',eg1,eg2,wetfac(iq)
          write(6,*) 'epot,egg,egg_alph1 ',epot(iq),egg(iq),egg_alph1
       endif
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

      if(nsigmf==0)then  ! original
!cdir nodep
        do ip=1,ipland  ! all land points in this nsib=3 loop
         iq=iperm(ip)
         ga(iq)=-slwa(iq)-rgg(iq)-fgg(iq)-cls(iq)*egg(iq)       
         dgdtg(iq)=-dirad(iq)-dfgdt(iq)-cls(iq)*degdt(iq)
        enddo   ! ip=1,ipland
      endif     ! (nsigmf==0)

      if(nsigmf==1)then  ! jlm preferred
!       spreads bare-soil flux across whole grid square      
!cdir nodep
        do ip=1,ipland  ! all land points in this nsib=3 loop
        iq=iperm(ip)
        ga(iq)=(-slwa(iq)-rgg(iq)-fgg(iq)-cls(iq)*egg(iq))*
     &                                     (1.-tsigmf(iq))
!       dgtdg is used in soilsnow
        dgdtg(iq)=-(dirad(iq)+dfgdt(iq)+cls(iq)*degdt(iq))*
     .                                  (1.-tsigmf(iq))
        enddo   ! ip=1,ipland
      endif     ! (nsigmf==1)

      if(ntest==1.and.mydiag)then
         iq=idjd
         write(6,*)'dgdtg,dirad,dfgdt,cls,degdt,tsigmf',
     &      dgdtg(iq),dirad(iq),dfgdt(iq),cls(iq),degdt(iq),tsigmf(iq) 
      endif

! ----------------------------------------------
      if(itnmeth==0) then  ! old, not vectorized
!cdir nodep
       do ip=1,ipland  ! all land points in this nsib=3 loop
        iq=iperm(ip)
        isoil = isoilm(iq)
        iveg=ivegt(iq)
c                                    components of the stomatal resistance
          sstar = 30.
          if( zo(iq) < .5 ) sstar = 150.
          ff= 1.1*sgsave(iq)/(rlai(iq)*sstar)
c         ff= 1.1*slwa(iq)/(rlai(iq)*sstar)
          rsi = rsmin(iq) * rlai(iq)
          f1= (1.+ff)/(ff+rsi/5000.)
          den=sfc(isoil)-swilt(isoil)                          ! sib3
          wbav=(max(0.,froot(1)*(wb(iq,1)-swilt(isoil)))+
     &         max(0.,froot(2)*(wb(iq,2)-swilt(isoil)))+
     &         max(0.,froot(3)*(wb(iq,3)-swilt(isoil)))+
     &         max(0.,froot(4)*(wb(iq,4)-swilt(isoil)))+
     &         max(0.,froot(5)*(wb(iq,5)-swilt(isoil)))   )/den
          f2=max(1. , .5/ max( wbav,1.e-7)) ! N.B. this is equiv to next 2 (jlm)
c         f2=1.0
c         if(wbav<0.5) f2=max(1.0 , 0.5/ max( wbav,1.e-7))
          f4=max(1.-.0016*(tstom(iq)-t(iq,1))**2 , .05) ! zero for delta_t=25
          airr(iq) = 1./taftfh(iq)
!         rmstep=dt/60.
          cc(iq) =min(condx(iq) , 4./(1440. *60./dt))  ! jlm speedup for 4 mm/day
c                       depth of the reservoir of water on the canopy
          rmcmax(iq) = max(0.5,srlai(iq)) * .1
          omc(iq) = cansto(iq)  ! value from previous timestep as starting guess
          f3=max(1.-.00025*(establ(t(iq,1))-qg(iq,1)*ps(iq)/.622), .05)
          res(iq)=max(30.,rsmin(iq)*f1*f2/(f3*f4))
          if(ntest==1.and.iq==idjd.and.mydiag)then
           write(6,*) 'rlai,srlai,wbav,den ',rlai(iq),srlai(iq),wbav,den
           write(6,*) 'f1,f2,f3,f4 ',f1,f2,f3,f4
           write(6,*) 'ff,f124,rsi,res ',ff,f1*f2/f4,rsi,res(iq)
          endif

          otgf(iq)=tgf(iq)
           do icount=1,5        ! original iteration method
c                                            transpiration
            cansto(iq) = omc(iq)
            esatf = establ(tgf(iq))
            qsatgf=.622*esatf/(ps(iq)-esatf)
c                                                     wet evaporation
            ewww(iq) = rho(iq)*taftfh(iq) *(qsatgf-qg(iq,1))
            if(qsatgf>=qg(iq,1)) then
              ewww(iq)  = min(cansto(iq)/dt , cansto(iq)*ewww(iq)
     &                        /rmcmax(iq) )
            endif         ! qsatgf>=qg(iq,1)
            cansto(iq)=omc(iq)+cc(iq) -ewww(iq)*dt
c                         precipitation reaching the ground under the canopy
c                                            water interception by the canopy
            condxg(iq)=max(condx(iq)-cc(iq)+
     .                 max(0.,cansto(iq)-rmcmax(iq)),0.) ! keep here
            cansto(iq) = min( max(cansto(iq),0.), rmcmax(iq))
            beta =      cansto(iq)/rmcmax(iq)
            if( qsatgf < qg(iq,1) ) beta = 1.

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
            evapxf(iq) = (evapfb(iq)/dt + ewww(iq))*hl
c                                              sensible heat flux
            prz = rho(iq)*taftfh(iq)*cp
            if(newfgf==0)fgf(iq)=prz*(tgf(iq)-theta(iq))  ! original/usual
            if(newfgf==1)fgf(iq)=prz*(tgf(iq)-tscrn(iq))
            if(newfgf==2)fgf(iq)=rho(iq)*aft(iq)*cp*
     .                                (tgf(iq)-tscrn(iq))
            rdg(iq) =  stefbo*tgf(iq)**4
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
            if(ntest==1.and.iq==idjd.and.mydiag)then
              write(6,*) 'icount,omc,rmc,otgf ',
     &                 icount,omc(iq),cansto(iq),otgf(iq)
              write(6,*) 'tfg,residf,evapxf,fgf ',
     &                 tgf(iq),residf(iq),evapxf(iq),fgf(iq)
            endif
            tgf(iq)=min(tgf(iq),sstgf+ 1.5)   ! jlm speedup
            tgf(iq)=max(tgf(iq),sstgf- 1.5)   ! jlm speedup

           enddo   !  icount=1,5
           if(ntest==1)write(6,*) 'iq,otgf(iq),tgf,residf '
     .                           ,iq,otgf(iq),tgf(iq),residf(iq)
110        continue
       enddo  !  ip=1,ipland
      endif  ! (itnmeth==0) 


      if(itnmeth>0) then  ! new, vectorized
!cdir nodep
       do ip=1,ipland  ! all land points in this nsib=3 loop
        iq=iperm(ip)
        isoil = isoilm(iq)
        iveg=ivegt(iq)
c                                  components of the stomatal resistance
        sstar=90.+sign(60.,.5-zo(iq))  ! equiv to above 2 lines
        ff= 1.1*sgsave(iq)/(rlai(iq)*sstar)
        rsi = rsmin(iq) * rlai(iq)
        f1= (1.+ff)/(ff+rsi/5000.)
        den=sfc(isoil)-swilt(isoil)                          ! sib3
        wbav=(max(0.,froot(1)*(wb(iq,1)-swilt(isoil)))+
     &        max(0.,froot(2)*(wb(iq,2)-swilt(isoil)))+
     &        max(0.,froot(3)*(wb(iq,3)-swilt(isoil)))+
     &        max(0.,froot(4)*(wb(iq,4)-swilt(isoil)))+
     &        max(0.,froot(5)*(wb(iq,5)-swilt(isoil)))   )/den
        f2=max(1. , .5/ max( wbav,1.e-7)) ! N.B. this is equiv to next 2 (jlm)
        f4=max(1.-.0016*(tstom(iq)-t(iq,1))**2 , .05) ! zero for delta_t=25
        airr(iq) = 1./taftfh(iq)
        cc(iq) =min(condx(iq) , 4./(1440. *60./dt))  ! jlm speedup for 4 mm/day
c                     depth of the reservoir of water on the canopy
        rmcmax(iq) = max(0.5,srlai(iq)) * .1
        omc(iq) = cansto(iq)  ! value from previous timestep as starting guess
        f3=max(1.-.00025*(establ(t(iq,1))-qg(iq,1)*ps(iq)/.622),.05)
        res(iq)=max(30.,rsmin(iq)*f1*f2/(f3*f4))
        if(ntest==1.and.iq==idjd.and.mydiag)then
          write(6,*) 'rlai,srlai,wbav,den ',rlai(iq),srlai(iq),wbav,den
          write(6,*) 'f1,f2,f3,f4 ',f1,f2,f3,f4
          write(6,*) 'ff,f124,rsi,res ',ff,f1*f2/f4,rsi,res(iq)
          write(6,*) 'qg,qfg,qlg ',qg(iq,1),qfg(iq,1),qlg(iq,1)
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
         ewwwa = rho(iq) *(qsatgf-qg(iq,1))/airr(iq) ! in W/m**2 /hl
!        max available dewfall is 
!        -(qsatgf-qg1)*dsig*1000*ps/grav in mm (mult by hl/dt for W/m**2)
         ewwwa=max(ewwwa,
     &             -abs((qsatgf-qg(iq,1))*dsig(1)*ps(iq))/(grav*dt))
         ewww(iq)=min(dt*ewwwa,omc(iq),dt*ewwwa*omc(iq)/rmcmax(iq))
!                     dew(-ve), no_dew ,        no_dew
!        cansto is reservoir on leaf
         cansto(iq)=(omc(iq)-ewww(iq)) +cc(iq)
         ewww(iq)=ewww(iq)/dt  ! these changes on 19/1/06 jlm

c                      precipitation reaching the ground under the canopy
c                                         water interception by the canopy
         condxg(iq)=max(condx(iq)-cc(iq)
     &             +max(0.,cansto(iq)-rmcmax(iq)),0.) ! keep 
         cansto(iq) = min( max(0.,cansto(iq)), rmcmax(iq))
         beta =      cansto(iq)/rmcmax(iq)
         Etr=rho(iq)*max(0.,qsatgf-qg(iq,1))/(airr(iq) +res(iq))  ! jlm
         betetrdt =(1.-beta)*Etr*dt*tsigmf(iq)   ! fixed 23/5/01
         evapfb1(iq)=min(betetrdt*froot(1),evapfb1a(iq))
         evapfb2(iq)=min(betetrdt*froot(2),evapfb2a(iq))
         evapfb3(iq)=min(betetrdt*froot(3),evapfb3a(iq))
         evapfb4(iq)=min(betetrdt*froot(4),evapfb4a(iq))
         evapfb5(iq)=min(betetrdt*froot(5),evapfb5a(iq))
         evapfb(iq)=(evapfb1(iq)+evapfb2(iq)+evapfb3(iq)+evapfb4(iq)+
     &               evapfb5(iq))/tsigmf(iq)
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
         if((ntest==1.or.diag).and.mydiag)then 
          if(ktau==99999)then  ! to track cansto & ewww underflow
            write(6,*) 'sflux icount = ',icount
            do ip=1,ipland  
             iq=iperm(ip)
             write(6,*) 'iq,otgf,tgfnew,ewww ',
     &                iq,otgf(iq),tgfnew(iq),ewww(iq)
            enddo
          endif
          if(land(idjd))then
           iq=idjd
           write(6,*) 'ktau,icount,iq,omc,cc ',
     &              ktau,icount,iq,omc(iq),cc(iq)
           write(6,*) 'rmc,rmcmax,ewww ',cansto(iq),rmcmax(iq),ewww(iq)
           esatf = establ(tgfnew(iq))  ! value for next itn
           qsatgf=.622*esatf/(ps(iq)-esatf)
           ewwwa = rho(iq) *(qsatgf-qg(iq,1))/airr(iq)
           write(6,*) 'esatf,qsatgf,ewwwa ',esatf,qsatgf,ewwwa
           prz = rho(iq)*cp*taftfh(iq)
           beta =      cansto(iq)/rmcmax(iq)
           devf= (hl*hlars/300.**2)*qsatgf*(1.-beta)/res(iq)
           dirad1 = 4.*rdg(iq)/300.
           write(6,*) 'beta,airr,res ',beta,airr(iq),res(iq)
           write(6,*) 'dirad1,devf,prz ',dirad1,devf,prz
           write(6,*) 'theta,tscrn,slwa ',theta(iq),tscrn(iq),slwa(iq)
           write(6,*) 'taftfh,condxg ',taftfh(iq),condxg(iq)
           write(6,*) 'rdg,fgf,evapxf,evapfb ',
     &              rdg(iq),fgf(iq),evapxf(iq),evapfb(iq)
           write(6,*) 'delta_tx ',delta_tx(iq)
           write(6,*) 'otgf,tgfnew,residf ',otgf(iq),tgfnew(iq)
     &              ,residf(iq)
          endif  ! (land(idjd))
         endif   ! ((ntest==2.or.diag).and.mydiag)
       enddo     !  icount=1,5
      endif      ! (itnmeth>0) 

!cdir nodep
      do ip=1,ipland  ! all land points in this nsib=3 loop
       iq=iperm(ip)
       if(cansto(iq)<1.e-10)cansto(iq)=0.  ! to avoid underflow 24/1/06
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
         condxpr(iq)=(1.-tsigmf(iq))*condx(iq)+tsigmf(iq)*condxg(iq)
         if(ntest==1.and.abs(residf(iq))>10.)
     &      write(6,*) 'iq,otgf(iq),tgf,delta_tx,residf '
     &              ,iq,otgf(iq),tgf(iq),delta_tx(iq),residf(iq)
       endif          ! tsigmf <= .01   ... else ...
       fev(iq)=evapfb(iq)/dt*hl*tsigmf(iq) ! passed to soilsnow to update wb
       fes(iq)=(1.-tsigmf(iq))*egg(iq)*cls(iq)  ! also passed to soilsnow
       otgsoil(iq)=isflag(iq)*tggsn(iq,1) + (1-isflag(iq))*tgg(iq,1)
      enddo  !  ip=1,ipland
      if((ntest==1.or.diag).and.mydiag) then
       if(land(idjd))then ! MJT bugfix
         iq=idjd
         isoil = isoilm(iq)
         iveg=ivegt(iq)
         write(6,*) 'in sib3 before soilsnowv'
         write(6,*) 'evapxf,epot,egg,fev,wetfac '
     .           ,evapxf(iq),epot(iq),egg(iq),fev(iq),wetfac(iq)
         write(6,*) 'fgf,fgg,fes ',fgf(iq),fgg(iq),fes(iq)
         write(6,*) 'isoil,ssat,tsigmf,rlai ',
     .            isoil,ssat(isoil),tsigmf(iq),rlai(iq)
         write(6,*) 'tgg1,t1,theta,tscrn '
     .           ,tgg(iq,1),t(iq,1),theta(iq),tscrn(iq)
         write(6,*) 'qg1,qsttg ',qg(iq,1),qsttg(iq)
         write(6,*) 'dfgdt,taftfhg,rho ',dfgdt(iq),taftfhg(iq),rho(iq)
         write(6,*) 'rmc,rmcmax(iq) ',cansto(iq),rmcmax(iq)
         if(abs(tgf(iq)-otgf(iq))>4.9)
     .      write(6,"('ktau,iq,otgf,tgf,dtgf,t1,t2',i4,i6,5f8.2)")
     .        ktau,iq,otgf(iq),tgf(iq),tgf(iq)-otgf(iq),t(iq,1),t(iq,2)
       endif
      endif
!-------------------------------------
c     write(6,*) 'before soilsnow'
      call soilsnowv
c     write(6,*) 'after soilsnow'
!cdir nodep
      do ip=1,ipland  ! all land points in this nsib=3 loop
       iq=iperm(ip)
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

        tgss=isflag(iq)*tggsn(iq,1) + (1-isflag(iq))*tgg(iq,1)  ! jlm
        if(tsigmf(iq)<= .01) then
          tss(iq) = tgss
          tgf(iq) = tgss
        else
          tss(iq)=tsigmf(iq)*tgf(iq)+(1.-tsigmf(iq))*tgss
        endif       ! tsigmf<= .01
        es = establ(tss(iq))     !  from 27/12/05
        qsttg(iq)= .622*es/(ps(iq)-es)  ! recal for scrnout, esp. snow    

      enddo   ! ip=1,ipland
      
      ! Calculate fraction of canopy which is wet
      fwet=0.
      where (land)
        fwet=cansto/rmcmax
      end where

      if((ntest==1.or.diag).and.mydiag) then
       if (land(idjd))then ! MJT bugfix
        iq=idjd
        write(6,*) 'even further down sib3 after soilsnowv'
        write(6,*) 'tgg ',(tgg(iq,k),k=1,ms)
        write(6,*) 'wb ',(wb(iq,k),k=1,ms)
        write(6,*) 'isflag,snowd ',isflag(iq),snowd(iq)
        write(6,*) 'evapfb,fev,ewww ',evapfb(iq),fev(iq),ewww(iq)
        write(6,*) 'tsigmf,evapxf,egg ',tsigmf(iq),evapxf(iq),egg(iq)
        write(6,*) 'deltat,degdt,wb,zse ',
     .           tgg(iq,1)-otgsoil(iq),degdt(iq),wb(iq,1),zse(1)
        write(6,*) 'eg,fg ',eg(iq),fg(iq)
       endif
      endif

      return
      end
