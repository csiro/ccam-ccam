      subroutine leoncld(cfrac,cffall,iaero)
      use aerointerface
      use arrays_m
      use cc_mpi, only : mydiag, myid
      use diag_m
      use estab
      use kuocomb_m
      use latlong_m
      use liqwpar_m  ! ifullw
      use morepbl_m
      use nharrs_m
      use nlin_m
      use prec_m
      use sigs_m
      use soil_m     ! land
      use tracers_m  ! ngas, nllp, ntrac
      use vvel_m
      use work3f_m
      implicit none
      include 'newmpar.h'
      integer  ncfrp,icfrp,iaero
      parameter (ncfrp=0,icfrp=1)        ! cfrp diags off      
!     parameter (ncfrp=1,icfrp=ifullw)   ! cfrp diags on     
      include 'const_phys.h' !Input physical constants
      include 'cparams.h'    !Input cloud scheme parameters
      include 'kuocom.h'  ! acon,bcon,Rcm
      include 'parm.h'
      include 'params.h'

c for cfrp
      real taul(icfrp,kl)
      real taui(icfrp,kl)
      real fice(icfrp,kl)
      integer kcldfmax(icfrp)
      real tautot,cldmax,ctoptmp,ctoppre
      real reffl,tau_sfac,wliq,rk,qlpath,wice,sigmai,cfd,fcf
      common/leoncfrp/tautot(icfrp),cldmax(icfrp)
     &               ,ctoptmp(icfrp),ctoppre(icfrp)

c Local variables
      integer iq,k,ncl
      real rainx,ccw,fl !Stuff for convective cloud fraction

      real prf(ifullw,kl)     !Pressure on full levels (hPa)
      real dprf(ifullw,kl)    !Pressure thickness (hPa)
      real rhoa(ifullw,kl)    !Air density (kg/m3)
      real dz(ifullw,kl)      !Layer thickness (m)
      real cdso4(ifullw,kl)   !Cloud droplet conc (#/m3)
      real cfrac(ifullw,kl)   !Cloud fraction (passed back to globpe)
      real cffall(ifullw,kl)  !Rain fraction (passed back to globpe)
      real ccov(ifullw,kl)    !Cloud cover (may differ from cloud frac if vertically subgrid)
      real cfa(ifullw,kl)     !Cloud fraction in which autoconv occurs (option in newrain.f)
      real qca(ifullw,kl)     !Cloud water mixing ratio in cfa(:,:)    (  "    "     "     )
      real fluxc(ifullw,kl)   !Flux of convective rainfall in timestep (kg/m**2)
      real ccrain(ifullw,kl)  !Convective raining cloud cover
      real precs(ifullw)      !Amount of stratiform precipitation in timestep (mm)
      real preci(ifullw)      !Amount of stratiform snowfall in timestep (mm)
      real cldcon(ifullw)     !Convective cloud fraction in column
      real wcon(ifullw)       !Convective cloud water content (in-cloud, prescribed)
      real clcon(ifullw,kl)   !Convective cloud fraction in layer 
      real qsg(ifullw,kl)     !Saturation mixing ratio
      real qcl(ifullw,kl)     !Vapour mixing ratio inside convective cloud
      real qenv(ifullw,kl)    !Vapour mixing ratio outside convective cloud
      real tenv(ifullw,kl)    !Temperature outside convective cloud
      real tnhs(ifullw,kl)    !Non-hydrostatic temperature adjusement
      real dumql(ifullw,kl),dumqf(ifullw,kl),dumqg(ifullw,kl)
      real dumt(ifullw,kl),dumqr(ifullw,kl),dumcr(ifullw,kl)

c These outputs are not used in this model at present
      real qevap(ifullw,kl)
      real qsubl(ifullw,kl)
      real qauto(ifullw,kl)
      real qcoll(ifullw,kl)
      real qaccr(ifullw,kl)
      real fluxr(ifullw,kl)
      real fluxi(ifullw,kl)
      real fluxm(ifullw,kl)
      real pqfsed(ifullw,kl)
      real pfstay(ifullw,kl)
      real slopes(ifullw,kl)
      real prscav(ifullw,kl)

      integer kbase(ifullw),ktop(ifullw) !Bottom and top of convective cloud 

      ! set-up params.h
      ln2=ifull
      lon=ln2/2
      nl=kl
      nlp=nl+1
      nlm=nl-1

      do k=1,kl   
         do iq=1,ifull
          prf(iq,k)=0.01*ps(iq)*sig(k) !Looks like ps is SI units
          dprf(iq,k)=-0.01*ps(iq)*dsig(k) !dsig is -ve
          rhoa(iq,k)=100.*prf(iq,k)/(rdry*t(iq,k))
          qsg(iq,k)=qsat(100.*prf(iq,k),t(iq,k))
        enddo
      enddo
      
      ! Non-hydrostatic terms
      tnhs(:,1)=phi_nh(:,1)/bet(1)
      do k=2,kl
        ! representing non-hydrostatic term as a correction to air temperature
        tnhs(:,k)=(phi_nh(:,k)-phi_nh(:,k-1)-betm(k)*tnhs(:,k-1))/bet(k)
      end do
      
      ! Calculate droplet concentration from aerosols (for non-convective faction of grid-box)
      call aerodrop(iaero,1,ifull,kl,cdso4,rhoa,land,rlatt,
     &              outconv=.true.)

      kbase(:)=0  ! default
      ktop(:) =0  ! default
      dz(:,:)=100.*dprf(:,:)/(rhoa(:,:)*grav)
     &        *(1.+tnhs(1:ifull,:)/t(1:ifull,:))
c     fluxc(:,:)=rnrt3d(:,:)*1.e-3*dt ! kg/m2 (should be same level as rnrt3d)
      fluxc(:,:)=0. !For now... above line may be wrong
      ccrain(:,:)=0.1  !Assume this for now
      precs(:)=0.
      preci(:)=0.

c     Set up convective cloud column
!     acon=0.2    !Cloud fraction for non-precipitating convection  kuocom.h
!     bcon=0.07   !Rate at which conv cloud frac increases with R   kuocom.h
      do iq=1,ifull
        if(ktsav(iq)<kl-1)then
          ktop(iq)=ktsav(iq)
          kbase(iq)=kbsav(iq)+1
          rainx=condc(iq)*86400./dt !mm/day
          cldcon(iq)=min(acon + bcon*log(1.0+rainx),0.8) !NCAR
          wcon(iq)=wlc
        else
          cldcon(iq)=0.
          wcon(iq)=0.
        endif
      enddo

!     if(diag.and.mydiag)then
      if(nmaxpr==1.and.mydiag)then
        if(ktau==1)print *,'in leoncloud acon,bcon,Rcm ',acon,bcon,Rcm
        print *,'entering leoncld'
        write (6,"('qg  ',3p9f8.3/4x,9f8.3)") qg(idjd,:)
        write (6,"('qf  ',3p9f8.3/4x,9f8.3)") qfg(idjd,:)
        write (6,"('ql  ',3p9f8.3/4x,9f8.3)") qlg(idjd,:)
      endif

c     Calculate convective cloud fraction and adjust moisture variables 
c     before calling newcloud
      if (nmr>=1) then ! Max/Rnd cloud overlap
        do k=1,kl
          where (k<=ktop(1:ifull).and.k>=kbase(1:ifull))
            clcon(1:ifull,k)=cldcon(1:ifull) ! maximum overlap
            !ccw=wcon(:)/rhoa(:,k)  !In-cloud l.w. mixing ratio
            qccon(1:ifull,k)=clcon(1:ifull,k)*
     &                 wcon(1:ifull)/rhoa(1:ifull,k)
            qcl(1:ifull,k)=max(qsg(1:ifull,k),qg(1:ifull,k))  ! jlm
            qenv(1:ifull,k)=max(1.e-8,
     &                 qg(1:ifull,k)-clcon(1:ifull,k)*qcl(1:ifull,k))
     &                 /(1.-clcon(1:ifull,k))
            qcl(1:ifull,k)=(qg(1:ifull,k)-(1.-clcon(1:ifull,k))
     &                *qenv(1:ifull,k))/clcon(1:ifull,k)
            qlg(1:ifull,k)=qlg(1:ifull,k)/(1.-clcon(1:ifull,k))
            qfg(1:ifull,k)=qfg(1:ifull,k)/(1.-clcon(1:ifull,k))
          elsewhere
            clcon(1:ifull,k)=0.
            qccon(1:ifull,k)=0.
            qcl(1:ifull,k)=0.
            qenv(1:ifull,k)=qg(1:ifull,k)
          endwhere
        enddo
      else ! usual random cloud overlap
      
        do k=1,kl
          do iq=1,ifull
            if(k<=ktop(iq).and.k>=kbase(iq))then
              ncl=ktop(iq)-kbase(iq)+1
              clcon(iq,k)=1.0-(1.0-cldcon(iq))**(1.0/ncl) !Random overlap
              ccw=wcon(iq)/rhoa(iq,k)  !In-cloud l.w. mixing ratio
!!27/4/04     qccon(iq,k)=clcon(iq,k)*ccw*0.25 ! 0.25 reduces updraft value to cloud value
              qccon(iq,k)=clcon(iq,k)*ccw
              !qcl(iq,k)=qsg(iq,k)
!             N.B. get silly qenv (becoming >qg) if qg>qsg (jlm)	     
              qcl(iq,k)=max(qsg(iq,k),qg(iq,k))  ! jlm
              qenv(iq,k)=max(1.e-8,
     &               qg(iq,k)-clcon(iq,k)*qcl(iq,k))/(1.-clcon(iq,k))
              qcl(iq,k)=(qg(iq,k)-(1.-clcon(iq,k))*qenv(iq,k))
     &               /clcon(iq,k)
              qlg(iq,k)=qlg(iq,k)/(1.-clcon(iq,k))
              qfg(iq,k)=qfg(iq,k)/(1.-clcon(iq,k))
            else
              clcon(iq,k)=0.
              qccon(iq,k)=0.
              qcl(iq,k)=0.
              qenv(iq,k)=qg(iq,k)
            endif
          enddo
        enddo
      
      end if
      
      tenv(:,:)=t(1:ifull,:) !Assume T is the same in and out of convective cloud

      if(nmaxpr==1.and.mydiag)then
        print *,'before newcloud',ktau
        write (6,"('t   ',9f8.2/4x,9f8.2)") t(idjd,:)
        write (6,"('qg  ',3p9f8.3/4x,9f8.3)") qg(idjd,:)
        write (6,"('qf  ',3p9f8.3/4x,9f8.3)") qfg(idjd,:)
        write (6,"('ql  ',3p9f8.3/4x,9f8.3)") qlg(idjd,:)
        write (6,"('qnv ',3p9f8.3/4x,9f8.3)") qenv(idjd,:)
        write (6,"('qsg ',3p9f8.3/4x,9f8.3)") qsg(idjd,:)
        write (6,"('qcl ',3p9f8.3/4x,9f8.3)") qcl(idjd,:)
        write (6,"('clc ',9f8.3/4x,9f8.3)") clcon(idjd,:)
        print *,'cldcon,kbase,ktop ',cldcon(idjd),kbase(idjd),ktop(idjd)
      endif

c     Calculate cloud fraction and cloud water mixing ratios
      dumql=qlg(1:ifull,:)
      dumqf=qfg(1:ifull,:)
      call newcloud(dt,1,land,prf,kbase,ktop,rhoa,cdso4, !Inputs
     &     tenv,qenv,dumql,dumqf,   !In and out  t here is tenv
     &     cfrac,ccov,cfa,qca)      !Outputs
      qlg(1:ifull,:)=dumql
      qfg(1:ifull,:)=dumqf

      if(nmaxpr==1.and.mydiag)then
        print *,'after newcloud',ktau
        write (6,"('tnv ',9f8.2/4x,9f8.2)") tenv(idjd,:)
        write (6,"('qg0 ',3p9f8.3/4x,9f8.3)") qg(idjd,:)
        write (6,"('qf  ',3p9f8.3/4x,9f8.3)") qfg(idjd,:)
        write (6,"('ql  ',3p9f8.3/4x,9f8.3)") qlg(idjd,:)
        write (6,"('qnv ',3p9f8.3/4x,9f8.3)") qenv(idjd,:) ! really new qg
      endif

c     Weight output variables according to non-convective fraction of grid-box            
      do k=1,kl
        do iq=1,ifull
          t(iq,k)=clcon(iq,k)*t(iq,k)+(1.-clcon(iq,k))*tenv(iq,k)
          qg(iq,k)=clcon(iq,k)*qcl(iq,k)+(1.-clcon(iq,k))*qenv(iq,k)
          if(k>=kbase(iq).and.k<=ktop(iq))then
            cfrac(iq,k)=cfrac(iq,k)*(1.-clcon(iq,k))
            ccov(iq,k)=ccov(iq,k)*(1.-clcon(iq,k))              
            qlg(iq,k)=qlg(iq,k)*(1.-clcon(iq,k))
            qfg(iq,k)=qfg(iq,k)*(1.-clcon(iq,k))
            cfa(iq,k)=cfa(iq,k)*(1.-clcon(iq,k))
            qca(iq,k)=qca(iq,k)*(1.-clcon(iq,k))              
          endif
        enddo
      enddo

      if(nmaxpr==1.and.mydiag)then
        print *,'before newrain',ktau
        write (6,"('t   ',9f8.2/4x,9f8.2)") t(idjd,:)
        write (6,"('qg  ',3p9f8.3/4x,9f8.3)") qg(idjd,:)
        write (6,"('qf  ',3p9f8.3/4x,9f8.3)") qfg(idjd,:)
        write (6,"('ql  ',3p9f8.3/4x,9f8.3)") qlg(idjd,:)
      endif
      if(diag)then
        call maxmin(t,' t',ktau,1.,kl)
        call maxmin(qg,'qg',ktau,1.e3,kl)
        call maxmin(qfg,'qf',ktau,1.e3,kl)
        call maxmin(qlg,'ql',ktau,1.e3,kl)
      endif

!     Add convective cloud water into fields for radiation
!     cfrad replaced by updating cfrac Oct 2005
!     Moved up 16/1/06 (and ccov,cfrac NOT UPDATED in newrain)
!     done because sometimes newrain drops out all qlg, ending up with 
!     zero cloud (although it will be rediagnosed as 1 next timestep)
      do k=1,kl
        do iq=1,ifull
c         cfrac(iq,k)=min(1.,ccov(iq,k)+clcon(iq,k))
          fl=max(0.0,min(1.0,(t(iq,k)-ticon)/(273.15-ticon)))
          qlrad(iq,k)=qlg(iq,k)+fl*qccon(iq,k)
          qfrad(iq,k)=qfg(iq,k)+(1.-fl)*qccon(iq,k)
        enddo
      enddo

c     Calculate precipitation and related processes
      dumt=t(1:ifull,:)
      dumql=qlg(1:ifull,:)
      dumqf=qfg(1:ifull,:)
      dumqr=qrg(1:ifull,:)
      dumqg=qg(1:ifull,:)
      dumcr=cffall(1:ifull,:)
      call newrain(land,1,dt,fluxc,rhoa,dz,ccrain,prf,cdso4,  !Inputs
     &    cfa,qca,                                            !Inputs
     &    dumt,dumql,dumqf,dumqr,
     &    precs,dumqg,cfrac,dumcr,ccov,                     !In and Out
     &    preci,qevap,qsubl,qauto,qcoll,qaccr,fluxr,fluxi,  !Outputs
     &    fluxm,pfstay,pqfsed,slopes,prscav)                !Outputs
      t(1:ifull,:)=dumt
      qlg(1:ifull,:)=dumql
      qfg(1:ifull,:)=dumqf
      qrg(1:ifull,:)=dumqr
      qg(1:ifull,:)=dumqg
      cffall(1:ifull,:)=dumcr
      if(nmaxpr==1.and.mydiag)then
        print *,'after newrain',ktau
        write (6,"('t   ',9f8.2/4x,9f8.2)") t(idjd,:)
        write (6,"('qg  ',3p9f8.3/4x,9f8.3)") qg(idjd,:)
        write (6,"('qf  ',3p9f8.3/4x,9f8.3)") qfg(idjd,:)
        write (6,"('ql  ',3p9f8.3/4x,9f8.3)") qlg(idjd,:)
      end if
      if(diag)then
        call maxmin(t,' t',ktau,1.,kl)
        call maxmin(qg,'qg',ktau,1.e3,kl)
        call maxmin(qfg,'qf',ktau,1.e3,kl)
        call maxmin(qlg,'ql',ktau,1.e3,kl)
      endif

      !--------------------------------------------------------------
      ! Store data needed by prognostic aerosol scheme
      if (abs(iaero)>=2) then
        ppfprec(:,1)=0. !At TOA
        ppfmelt(:,1)=0. !At TOA
        ppfsnow(:,1)=0. !At TOA
        ppfconv(:,1)=0. !At TOA
        do k=1,kl-1
          ppfprec(:,kl+1-k)=(fluxr(:,k+1)+fluxm(:,k))/dt !flux *entering* layer k
          ppfmelt(:,kl+1-k)=fluxm(:,k)/dt                !flux melting in layer k
          ppfsnow(:,kl+1-k)=(fluxi(:,k+1)-fluxm(:,k))/dt !flux *entering* layer k
          ppfconv(:,kl+1-k)=fluxc(:,k)/dt                !flux *leaving* layer k
        enddo
        do k=1,kl
          ppfevap(:,kl+1-k)=qevap(:,k)*rhoa(:,k)*dz(:,k)/dt
          ppfsubl(:,kl+1-k)=qsubl(:,k)*rhoa(:,k)*dz(:,k)/dt !flux sublimating or staying in k
          pplambs(:,kl+1-k)=slopes(:,k)
          where (qlg(:,k)+qfg(:,k)>1.e-8)
            ppmrate(:,kl+1-k)=(qauto(:,k)+qcoll(:,k))/dt
            ppmaccr(:,kl+1-k)=qaccr(:,k)/dt
          elsewhere
            ppmrate(:,kl+1-k)=0.
            ppmaccr(:,kl+1-k)=0.
          end where
        enddo
        ppfstay=pfstay
        ppqfsed=pqfsed
        pprscav=prscav
      end if
      !--------------------------------------------------------------

!     Add convective cloud water into fields for radiation
!     cfrad replaced by updating cfrac Oct 2005
!     Moved up 16/1/06 (and ccov,cfrac NOT UPDATED in newrain)
!     done because sometimes newrain drops out all qlg, ending up with 
!     zero cloud (although it will be rediagnosed as 1 next timestep)
      do k=1,kl
        do iq=1,ifull
          cfrac(iq,k)=min(1.,ccov(iq,k)+clcon(iq,k))
c          fl=max(0.0,min(1.0,(t(iq,k)-ticon)/(273.15-ticon)))
c          qlrad(iq,k)=qlg(iq,k)+fl*qccon(iq,k)
c          qfrad(iq,k)=qfg(iq,k)+(1.-fl)*qccon(iq,k)
        enddo
      enddo

!========================= Jack's diag stuff =========================
       if(ncfrp==1)then  ! from here to near end; Jack's diag stuff
        do iq=1,icfrp
          tautot(iq)=0.
          cldmax(iq)=0.
          ctoptmp(iq)=0.
          ctoppre(iq)=0.
          do k=1,kl
            fice(iq,k)=0.
          enddo
          kcldfmax(iq)=0.
        enddo
c       cfrp data
        do k=1,kl-1
          do iq=1,icfrp
            taul(iq,k)=0.
            taui(iq,k)=0.
            Reffl=0.
            if(cfrac(iq,k)>0.)then
              tau_sfac=1.
c             fice(iq,k) = qfg(iq,k)/(qfg(iq,k)+qlg(iq,k))
              fice(iq,k) = qfrad(iq,k)/(qfrad(iq,k)+qlrad(iq,k)) ! 16/1/06
              !rhoa=100.*prf(iq,k)/(rdry*t(iq,k))
              !dpdp=dprf(iq,k)/prf(iq,k)
              !dz=dpdp*rdry*t(iq,k)/grav
c             Liquid water clouds
              if(qlg(iq,k)>1.0e-8)then
                Wliq=rhoa(iq,k)*qlg(iq,k)/(cfrac(iq,k)*(1-fice(iq,k))) !kg/m^3
                if(.not.land(iq))then !sea
                  rk=0.8
                else            !land
                  rk=0.67
                endif
c Reffl is the effective radius at the top of the cloud (calculated following
c Martin etal 1994, JAS 51, 1823-1842) due to the extra factor of 2 in the
c formula for reffl. Use mid cloud value of Reff for emissivity.
                Reffl=(3*2*Wliq/(4*rhow*pi*rk*cdso4(iq,k)))**(1./3)
                qlpath=Wliq*dz(iq,k)
                taul(iq,k)=tau_sfac*1.5*qlpath/(rhow*Reffl)
              endif ! qlg
c Ice clouds
              if(qfg(iq,k)>1.0e-8)then
                Wice=rhoa(iq,k)*qfg(iq,k)/(cfrac(iq,k)*fice(iq,k)) !kg/m**3
                sigmai = aice*Wice**bice !visible ext. coeff. for ice
                taui(iq,k)=sigmai*dz(iq,k) !visible opt. depth for ice
                taui(iq,k)=tau_sfac*taui(iq,k)
              endif ! qfg
            endif !cfrac
          enddo ! iq
        enddo ! k
c Code to get vertically integrated value...
c top down to get highest level with cfrac=cldmax (kcldfmax)
!        do k=kl,1,-1
!       Need to start at kl-1 to work with NaN initialisation.
        do k=kl-1,1,-1
          do iq=1,icfrp
            tautot(iq)=tautot(iq)+cfrac(iq,k)*(fice(iq,k)
     &                *taui(iq,k)+(1.-fice(iq,k))*taul(iq,k))
            if(cfrac(iq,k)>cldmax(iq)) kcldfmax(iq)=k
            cldmax(iq)=max(cldmax(iq),cfrac(iq,k))
          enddo ! iq
        enddo ! k

        do iq=1,icfrp
          if(cldmax(iq)>1.e-10) then
            tautot(iq)=tautot(iq)/cldmax(iq)

            cfd=0.
            do k=kl,kcldfmax(iq),-1
              fcf = max(0.,cfrac(iq,k)-cfd) ! cld frac. from above
              ctoptmp(iq)=ctoptmp(iq)+fcf*t(iq,k)/cldmax(iq)
              ctoppre(iq)=ctoppre(iq)+fcf*prf(iq,k)/cldmax(iq)
              cfd=max(cfrac(iq,k),cfd)
            enddo ! k=kl,kcldfmax(iq),-1

          endif ! (cldmax(iq).gt.1.e-10) then
        enddo   ! iq
       endif    ! ncfrp.eq.1
!========================= end of Jack's diag stuff ======================

!     factor of 2 is because used .5 in newrain.f (24-mar-2000, jjk)
      do iq=1,ifullw        
        condx(iq)=condx(iq)+precs(iq)*2.
        conds(iq)=conds(iq)+preci(iq)*2.
        precip(iq)=precip(iq)+precs(iq)*2.
      enddo

      return
      end
