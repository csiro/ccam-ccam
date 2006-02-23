      subroutine leoncld(cfrac)
      use diag_m
      use cc_mpi, only : mydiag
      implicit none
      include 'newmpar.h'
      include 'liqwpar.h' ! ifullw
      integer  ncfrp,icfrp
      parameter (ncfrp=0,icfrp=1)        ! cfrp diags off      
!     parameter (ncfrp=1,icfrp=ifullw)   ! cfrp diags on     
      include 'const_phys.h' !Input physical constants
      include 'cparams.h'    !Input cloud scheme parameters

      include 'arrays.h'
c     include 'constant.h'
      include 'dava.h'    ! davt
      include 'kuocom.h'  ! acon,bcon,Rcm
      include 'latlong.h' ! rlatt,rlongg
      include 'map.h'     ! land
      include 'morepbl.h'
      include 'nlin.h'
      include 'parm.h'
      include 'prec.h'
      include 'sigs.h'
      include 'soil.h'
      include 'tracers.h'  ! ngas, nllp, ntrac
      include 'vvel.h'

c for cfrp
      real taul(icfrp,kl)
      real taui(icfrp,kl)
      real fice(icfrp,kl)
      integer kcldfmax(icfrp)
      real tautot,cldmax,ctoptmp,ctoppre
      real reffl,tau_sfac,wliq,rk,qlpath,wice,sigmai,cfd,fcf
      common/leoncfrp/tautot(icfrp),cldmax(icfrp)
     &               ,ctoptmp(icfrp),ctoppre(icfrp)
      common/work3f/qccon(ifull,kl),qlrad(ifull,kl),qfrad(ifull,kl) ! ditto
      real qccon,qlrad,qfrad

c Local variables
      integer iq,k,ncl
      real rainx,ccw,fl !Stuff for convective cloud fraction

      real prf(ifullw,kl)     !Pressure on full levels (hPa)
      real dprf(ifullw,kl)    !Pressure thickness (hPa)
      real rhoa(ifullw,kl)    !Air density (kg/m3)
      real dz(ifullw,kl)      !Layer thickness (m)
      real cdso4(ifullw,kl)   !Cloud droplet conc (#/m3)
      real cfrac(ifullw,kl)   !Cloud fraction (passed back to globpe)
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
      include 'establ.h'

      do k=1,kl   
         do iq=1,ifull
          prf(iq,k)=0.01*ps(iq)*sig(k) !Looks like ps is SI units
          dprf(iq,k)=-0.01*ps(iq)*dsig(k) !dsig is -ve
          rhoa(iq,k)=100.*prf(iq,k)/(rdry*t(iq,k))
          qsg(iq,k)=qsat(100.*prf(iq,k),t(iq,k))
          if(land(iq))then
            if(rlatt(iq)>0.)then	     
              cdso4(iq,k)=cdropl_nh
	     else
              cdso4(iq,k)=cdropl_sh
	     endif
          else
            if(rlatt(iq)>0.)then	     
              cdso4(iq,k)=cdrops_nh
	     else
              cdso4(iq,k)=cdrops_sh
	     endif
          endif  ! (land(iq)) .. else ..
        enddo
      enddo

      kbase(:)=0  ! default
      ktop(:) =0  ! default
      dz(:,:)=100.*dprf(:,:)/(rhoa(:,:)*grav)
c     fluxc(:,:)=rnrt3d(:,:)*1.e-3*dt ! kg/m2 (should be same level as rnrt3d)
      fluxc(:,:)=0. !For now... above line may be wrong
      ccrain(:,:)=0.1  !Assume this for now
      precs(:)=0.

c     Set up convective cloud column
!     acon=0.2    !Cloud fraction for non-precipitating convection  kuocom.h
!     bcon=0.07   !Rate at which conv cloud frac increases with R   kuocom.h
      do iq=1,ifull
        if(ktsav(iq).lt.kl)then
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
        if(ktau.eq.1)print *,'in leoncloud acon,bcon,Rcm ',acon,bcon,Rcm
        print *,'entering leoncld'
        write (6,"('qg  ',3p9f8.3/4x,9f8.3)") qg(idjd,:)
        write (6,"('qf  ',3p9f8.3/4x,9f8.3)") qfg(idjd,:)
        write (6,"('ql  ',3p9f8.3/4x,9f8.3)") qlg(idjd,:)
      endif

c     Calculate convective cloud fraction and adjust moisture variables 
c     before calling newcloud
      do k=1,kl
        do iq=1,ifull
          if(k.le.ktop(iq).and.k.ge.kbase(iq))then
            ncl=ktop(iq)-kbase(iq)+1
            clcon(iq,k)=1.0-(1.0-cldcon(iq))**(1.0/ncl) !Random overlap
            ccw=wcon(iq)/rhoa(iq,k)  !In-cloud l.w. mixing ratio
!!27/4/04   qccon(iq,k)=clcon(iq,k)*ccw*0.25 ! 0.25 reduces updraft value to cloud value
            qccon(iq,k)=clcon(iq,k)*ccw
            qcl(iq,k)=qsg(iq,k)
!           N.B. get silly qenv (becoming >qg) if qg>qsg (jlm)	     
            qcl(iq,k)=max(qsg(iq,k),qg(iq,k))  ! jlm
            qenv(iq,k)=max(1.e-8,
     &                 qg(iq,k)-clcon(iq,k)*qcl(iq,k))/(1-clcon(iq,k))
            qcl(iq,k)=(qg(iq,k)-(1-clcon(iq,k))*qenv(iq,k))/clcon(iq,k)
            qlg(iq,k)=qlg(iq,k)/(1-clcon(iq,k))
            qfg(iq,k)=qfg(iq,k)/(1-clcon(iq,k))
          else
            clcon(iq,k)=0.
            qccon(iq,k)=0.
            qcl(iq,k)=0.
            qenv(iq,k)=qg(iq,k)
          endif
        enddo
      enddo
      tenv(:,:)=t(1:ifull,:) !Assume T is the same in and out of convective cloud
!     if(diag.and.mydiag)then
      if(nmaxpr==1.and.mydiag)then
        print *,'before newcloud'
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
      call newcloud(dt,1,land,prf,kbase,ktop,rhoa,cdso4, !Inputs
     &     tenv,qenv,qlg(1:ifull,:),qfg(1:ifull,:),   !In and out  t here is tenv
     &     cfrac,ccov,cfa,qca)   !Outputs
!     if(diag.and.mydiag)then
      if(nmaxpr==1.and.mydiag)then
        print *,'after newcloud'
        write (6,"('tnv ',9f8.2/4x,9f8.2)") tenv(idjd,:)
        write (6,"('qg  ',3p9f8.3/4x,9f8.3)") qg(idjd,:)
        write (6,"('qf  ',3p9f8.3/4x,9f8.3)") qfg(idjd,:)
        write (6,"('ql  ',3p9f8.3/4x,9f8.3)") qlg(idjd,:)
        write (6,"('qnv ',3p9f8.3/4x,9f8.3)") qenv(idjd,:)
      endif

c     Weight output variables according to non-convective fraction of grid-box            
      do k=1,kl
        do iq=1,ifull
          t(iq,k)=clcon(iq,k)*t(iq,k)+(1-clcon(iq,k))*tenv(iq,k)
          qg(iq,k)=clcon(iq,k)*qcl(iq,k)+(1-clcon(iq,k))*qenv(iq,k)
          if(k.ge.kbase(iq).and.k.le.ktop(iq))then
            cfrac(iq,k)=cfrac(iq,k)*(1-clcon(iq,k))
            ccov(iq,k)=ccov(iq,k)*(1-clcon(iq,k))              
            qlg(iq,k)=qlg(iq,k)*(1-clcon(iq,k))
            qfg(iq,k)=qfg(iq,k)*(1-clcon(iq,k))
            cfa(iq,k)=cfa(iq,k)*(1-clcon(iq,k))
            qca(iq,k)=qca(iq,k)*(1-clcon(iq,k))              
          endif
        enddo
      enddo
      if(nmaxpr==1.and.mydiag)then
        print *,'before newrain'
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

c     Calculate precipitation and related processes      
      call newrain(land,1,dt,fluxc,rhoa,dz,ccrain,prf,cdso4,  !Inputs
     &    cfa,qca,                                            !Inputs
     &    t(1:ifull,:),qlg(1:ifull,:),qfg(1:ifull,:),
     &    precs,qg(1:ifull,:),cfrac,ccov, !In and Out
     &    preci,qevap,qsubl,qauto,qcoll,qaccr,fluxr,fluxi,  !Outputs
     &    fluxm,pfstay,pqfsed,slopes,prscav)     !Outputs
      if(nmaxpr==1.and.mydiag)then
        print *,'after newrain'
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

!========================= Jack's diag stuff =========================
!      if(ncfrp.eq.1)then  ! from here to near end; Jack's diag stuff
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
            if(cfrac(iq,k).gt.0)then
              tau_sfac=1.
              fice(iq,k) = qfg(iq,k)/(qfg(iq,k)+qlg(iq,k))
              !rhoa=100.*prf(iq,k)/(rdry*t(iq,k))
              !dpdp=dprf(iq,k)/prf(iq,k)
              !dz=dpdp*rdry*t(iq,k)/grav
c             Liquid water clouds
              if(qlg(iq,k).gt.1.0e-8)then
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
              if(qfg(iq,k).gt.1.0e-8)then
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
            if(cfrac(iq,k).gt.cldmax(iq)) kcldfmax(iq)=k
            cldmax(iq)=max(cldmax(iq),cfrac(iq,k))
          enddo ! iq
        enddo ! k

        do iq=1,icfrp
          if(cldmax(iq).gt.1.e-10) then
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
!      endif    ! ncfrp.eq.1
!========================= end of Jack's diag stuff ======================

c     Add convective cloud water into fields for radiation
!     cfrad replaced by updating cfrac Oct 2005
      do k=1,kl
        do iq=1,ifull
          cfrac(iq,k)=min(1.,ccov(iq,k)+clcon(iq,k))
          fl=max(0.0,min(1.0,(t(iq,k)-ticon)/(273.15-ticon)))
          qlrad(iq,k)=qlg(iq,k)+fl*qccon(iq,k)
          qfrad(iq,k)=qfg(iq,k)+(1.-fl)*qccon(iq,k)
        enddo
      enddo

!     factor of 2 is because used .5 in newrain.f (24-mar-2000, jjk)
      do iq=1,ifullw        
        condx(iq)=condx(iq)+precs(iq)*2.
        precip(iq)=precip(iq)+precs(iq)*2.
      enddo

      return
      end
