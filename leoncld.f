      subroutine leoncld(cfrac)

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
      include 'kuocom.h'  ! kbsav,ktsav,convfact,convpsav,ndavconv
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
      real reffl,tau_sfac,wliq,rk,cdrop,qlpath,wice,sigmai,cfd,fcf
      common/leoncfrp/tautot(icfrp),cldmax(icfrp)
     &               ,ctoptmp(icfrp),ctoppre(icfrp)


c Local variables
      integer iq,k

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
          if(land(iq))then
            cdso4(iq,k)=cdropl
          else
            cdso4(iq,k)=cdrops
          endif
        enddo
      enddo

      kbase(:)=0  !Not used for now
      ktop(:)=0
      dz(:,:)=100.*dprf(:,:)/(rhoa(:,:)*grav)
c     fluxc(:,:)=rnrt3d(:,:)*1.e-3*dt ! kg/m2 (should be same level as rnrt3d)
      fluxc(:,:)=0. !For now... above line may be wrong
      ccrain(:,:)=0.1  !Assume this for now
      precs(:)=0.

c Calculate cloud fraction and cloud water mixing ratios

      call newcloud(dt,1,land,prf,kbase,ktop,rhoa,cdso4, !Inputs
     &     t,qg,qlg,qfg,   !In and out
     &     cfrac,ccov,cfa,qca)   !Outputs


c Calculate precipitation and related processes
      
      call newrain(land,1,dt,fluxc,rhoa,dz,ccrain,prf,cdso4,  !Inputs
     &    cfa,qca,                                            !Inputs
     &    t,qlg,qfg,precs,qg,cfrac,ccov,                      !In and Out
     &    preci,qevap,qsubl,qauto,qcoll,qaccr,fluxr,fluxi,  !Outputs
     &    fluxm,pfstay,pqfsed,slopes,prscav)     !Outputs

      
! added by jjk 27-07-03
! factor of 2 is because used .5 in newrain.f (24-mar-2000, jjk)
      do iq=1,ifullw        
        condx(iq)=condx(iq)+precs(iq)*2.
        precip(iq)=precip(iq)+precs(iq)*2.
      enddo

!=======================================================================
!      if(ncfrp.eq.1)then  ! from here to the end
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
                  Cdrop=Cdrops
                else            !land
                  rk=0.67
                  Cdrop=Cdropl
                endif
c Reffl is the effective radius at the top of the cloud (calculated following
c Martin etal 1994, JAS 51, 1823-1842) due to the extra factor of 2 in the
c formula for reffl. Use mid cloud value of Reff for emissivity.
                Reffl=(3*2*Wliq/(4*rhow*pi*rk*Cdrop))**(1./3)
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
        do k=kl,1,-1
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

      return
      end
