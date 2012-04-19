c ** Based on GCM revision 1.54 **
c This routine is part of the prognostic cloud scheme. It calculates rainfall
c and the evaporation of rain, and also calls icefall, which does the frozen 
c precipitation. It is called by progcld.
c
c INPUT/OUTPUT
c
c Input:
c
c parameters from include file params.h
c      lon - number of points around a latitude circle
c      ln2 - number of points for NH+SH latitude circles sequentially
c      nl - number of vertical levels
c
c      ntest - 0=off, 1=on to control single column debugging
c
c see also include files physparams.h (model physical constants)
c                        cparams.h    (cloud scheme parameters)
c
c from arguments
c      land - logical variable for surface type ( = T for land points)
c      tdt - leapfrog timestep (seconds)
c      lg - latitude index (ranges from 1 at poles to LAT at equator)
c      fluxc - flux of convective rain in timestep (kg/m**2)
c      rhoa - air density (kg/m**3)
c      dz - layer thicknes (m)
c      ccrain - raining convective cloud fraction at level k
c      prf - pressure at full levels (in hPa. NB: not SI units)
c
c In/Out:
c
c from arguments
c      ttg - temperature (K)
c      qlg - cloud liquid water mixing ratio (kg/kg)
c      qfg - cloud ice mixing ratio (kg/kg)
c      precs - amount of stratiform precipitation in timestep (mm)
c      qtg - water vapour mixing ratio (kg/kg) - called qg in C-CAM
c      cfrac - stratiform cloud fraction
c      ccov - stratiform cloud *cover* looking from above (currently = cfrac)
c
c Output:
c
c from arguments
c      preci - amount of stratiform snowfall in timestep (mm)
c      qevap - evaporation of rainfall (kg/kg)
c      qsubl - sublimation of snowfall (kg/kg)
c      qauto - autoconversion of cloud liquid water (kg/kg)
c      qcoll - collection by rain of cloud liquid water (kg/kg)
c      qaccr - accretion by snow of cloud liquid water (kg/kg)
c
***************************************************************************

      subroutine newrain(land,lg,tdt,fluxc,rhoa,dz,ccrain,prf,cdrop,!Inputs
     &                  cfa,qca,
     &                  ttg,qlg,qfg,qrg,precs,qtg,cfrac,cffall,ccov,   !In and Out
     &                  preci,qevap,qsubl,qauto,qcoll,qaccr,fluxr,fluxi,
     &                  fluxm,pfstay,pqfsed,slopes,prscav)     !Outputs

      use cc_mpi, only : mydiag
      use kuocomb_m
      use morepbl_m  !condx        
      implicit none
C Global parameters
      include 'newmpar.h'
      include 'const_phys.h' !Input physical constants
      include 'cparams.h'    !Input cloud scheme parameters
      include 'kuocom.h'     !acon,bcon,Rcm,ktsav,nevapls
      include 'params.h'     !Input model grid dimensions (modified PARAMS.f for CCAM)
      include 'parm.h'

C Argument list
      logical land(ln2)
      integer lg
      real tdt
      real fluxc(ln2,nl)
      real rhoa(ln2,nl)
      real dz(ln2,nl)
      real ccrain(ln2,nl)
      real prf(ln2,nl)
      real ttg(ln2,nl)
      real qlg(ln2,nl)
      real qfg(ln2,nl)
      real qrg(ln2,nl)
      real precs(ln2)
      real qtg(ln2,nl)
      real cfrac(ln2,nl)
      real cffall(ln2,nl)
      real ccov(ln2,nl)
      real preci(ln2)
      real qevap(ln2,nl)
      real qsubl(ln2,nl)
      real qauto(ln2,nl)
      real qcoll(ln2,nl)
      real qaccr(ln2,nl)
      real cfa(ln2,nl)
      real qca(ln2,nl)
      real pqfsed(ln2,nl)
      real pfstay(ln2,nl)
      real slopes(ln2,nl)
      real prscav(ln2,nl)

C Local work arrays and variables
      real clfr(ln2,nl)
      real qsg(ln2,nl),fluxr(ln2,nl)
      real cifr(ln2,nl)
      real fluxm(ln2,nl)
      real frclr(ln2,nl)
      real fluxi(ln2,nl),qsl(ln2,nl)
      real fluxa(ln2,nl)
      real clfra(ln2)
      real ccra(ln2)
      real cfrain(ln2,nl)
      real cfmelt(ln2,nl)
      real fluxrain(ln2)
      real cdrop(ln2,nl)
      real fracr(ln2,nl)
      real mxclfr(ln2)
      real rdclfr(ln2)
      real vr(ln2,nl-1)
      real rhor(ln2,nl-1)
      real fout(ln2,nl-1)
      real fthru(ln2,nl-1)

      integer k,mb,mg,njumps,ns,nt,ntest

      real apr
      real bpr
      real bl
      real cdt
      real cev
      real clrevap
      real coll
      real crate
      real delt
      real dql
      real dqlo
      real dqsdt
      real es
      real evap
      real fcol
      real fr
      real frb
      real frc
      real pk
      real qcic
      real qcrit
      real ql
      real ql1
      real ql2
      real qpf
      real R6c,R3c,beta6,eps
      real rhodz
      real satevap
      real selfcoll
      real tk
      real Wliq
      real cfla,dqla,qla
      real cftemp,qrgtemp,mixrain,alph
      real rhorin,rhorout
      real cffluxin,cffluxout

C Local data, functions etc
      parameter (ntest=0)  ! 0=off, 1=on
      integer, parameter :: autoconvmeth=1 ! Autoconv scheme (0=Old, 1=New)

      real esdiff(-40:0)  !Table of es(liq) - es(ice) (MKS), -40 <= t <= 0 C
      data esdiff / 
     & 6.22, 6.76, 7.32, 7.92, 8.56, 9.23, 9.94,10.68,11.46,12.27,
     & 13.11,13.99,14.89,15.82,16.76,17.73,18.70,19.68,20.65,21.61,
     & 22.55,23.45,24.30,25.08,25.78,26.38,26.86,27.18,27.33,27.27,
     & 26.96,26.38,25.47,24.20,22.51,20.34,17.64,14.34,10.37, 5.65,
     & 0.08 /
      real tablei
      common /esitable/ tablei(0:220) !Table of es values wrt ice

c Arithmetic statement functions to replace call to establ.
c T is temp in Kelvin, which should lie between 123.16 and 343.16;
c TDIFF is difference between T and 123.16, subject to 0 <= TDIFF <= 220

      real t, tdiff, pp, estabi, qsati
      tdiff(t)=min(max( t-123.16, 0.), 219.)

c These give the ice values needed for the qcloud scheme
      estabi(t) = (1.-(tdiff(t)-aint(tdiff(t))))*tablei(int(tdiff(t)))
     &           + (tdiff(t)-aint(tdiff(t)))*tablei(int(tdiff(t))+1)
      qsati(pp,t) = epsil*estabi(t)/max(.1,pp-estabi(t)) !Usual formula

      real pow75,x
      pow75(x)=sqrt(x*sqrt(x))

C Start code : ----------------------------------------------------------

        if(ntest>0.and.mydiag)then
          mg=idjd
          write(25,'(a,3i3)')'IPASS=1, Before newrain ktau= ',ktau
          write(25,91)'rhoa',(rhoa(mg,k),k=1,nl)
          write(25,91)'dz',(dz(mg,k),k=1,nl)
          write(25,91)'ccrain',(ccrain(mg,k),k=1,nl)
          write(25,91)'prf',(prf(mg,k),k=1,nl)
          write(25,9)'cdrop',(cdrop(mg,k),k=1,nl)
          write(25,91)'ttg ',(ttg(mg,k),k=1,nl)
          write(25,9)'qtg ',(qtg(mg,k),k=1,nl)
          write(25,9)'qsg ',(qsg(mg,k),k=1,nl)
          write(25,9)'cfrac ',(cfrac(mg,k),k=1,nl)
          write(25,9)'qlg ',(qlg(mg,k),k=1,nl)
          write(25,9)'qfg ',(qfg(mg,k),k=1,nl)
          write(25,9)'cfa ',(cfa(mg,k),k=1,nl)
          write(25,9)'qca ',(qca(mg,k),k=1,nl)
          write(25,9)'fluxc ',(fluxc(mg,k),k=1,nl)
          write(25,*)
        endif  ! (ntest>0)

      delt=tdt
      njumps=nint(tdt/delt)
      delt=tdt/real(njumps) !To make sure tdt it a multiple of delt

      do k=1,nl
        do mg=1,ln2
          fracr(mg,k)=0.
          fluxr(mg,k)=0.
          frclr(mg,k)=0.
          fluxm(mg,k)=0.
          fluxi(mg,k)=0.
          fluxa(mg,k)=0.
          qevap(mg,k)=0.
          qauto(mg,k)=0.
          qcoll(mg,k)=0.
          cfrain(mg,k)=0.
          pk=100.0*prf(mg,k)
          qsg(mg,k)=qsati(pk,ttg(mg,k))
          if(qfg(mg,k).gt.0.)then
            cifr(mg,k)=cfrac(mg,k)*qfg(mg,k)/(qlg(mg,k)+qfg(mg,k))
          else
            cifr(mg,k)=0.
          endif
c          clfr(mg,k)=max(cfrac(mg,k)-cifr(mg,k),0.) 
c Previous line not good for roundoff; next few lines are better
          if(qlg(mg,k).gt.0.)then
            clfr(mg,k)=cfrac(mg,k)*qlg(mg,k)/(qlg(mg,k)+qfg(mg,k))
          else
            clfr(mg,k)=0.
          endif
        enddo
      enddo

c Define cdrop  - passed through as cdso4, defined in leoncld.f

***************** Cut here to insert new auto scheme ********************            
      select case(autoconvmeth)
c Using new (subgrid) autoconv scheme... 
       case(1)
        do k=nl-1,1,-1
          do mg=1,ln2
            cfrain(mg,k)=0.0
            rhodz=rhoa(mg,k)*dz(mg,k)

            if(clfr(mg,k).gt.0.)then

              ql=qlg(mg,k)
              cfla=0.
              dqla=0.
              if(cfa(mg,k).gt.0.)then
                cfla=cfa(mg,k)*clfr(mg,k)/(clfr(mg,k)+cifr(mg,k))
                qla=qca(mg,k)/cfa(mg,k)

c Following few lines are for Yangang Liu's new scheme (2004: JAS, GRL)

                Wliq = 1000. * qla * rhoa(mg,k) !g/m3
                R6c = 4.09e-4
     &               * ( 1.15e23*1.e-6*cdrop(mg,k) / Wliq**2 ) ** (1./6)
                eps = 1. - 0.7 * exp(-0.003e-6*cdrop(mg,k)) !mid range
                beta6 = ((1.+3*eps**2)*(1.+4*eps**2)*(1.+5*eps**2)
     &               / ((1.+eps**2)*(1.+2*eps**2)) )**(1./6)
                R3c = 1.e-6*R6c/beta6 !in metres
                qcrit=(4*pi/3)*rhow*R3c**3*Cdrop(mg,k)/rhoa(mg,k) !New qcrit

                if(qla.le.qcrit)then
                  ql2=qla
                else

c Following is Liu & Daum (JAS, 2004)
                  Crate=1.9e17*
     &               (0.75*rhoa(mg,k)/(pi*rhow))**2*beta6**6/cdrop(mg,k)
                  ql1=qla/sqrt(1.+2*crate*qla**2*delt)

                  ql1=max(ql1, qcrit) !Intermediate qlg after auto
                  Frb=dz(mg,k)*rhoa(mg,k)*(qla-ql1)/delt
                  cdt=delt*0.5*Ecol*0.24*pow75(Frb)
                  selfcoll=min(ql1,ql1*cdt)
                  ql2=ql1-selfcoll
                  cfrain(mg,k)=cfla
                endif
                dqla=cfla*(qla-ql2)
                ql=max(1.e-8,qlg(mg,k)-dqla)
              endif
              dql=qlg(mg,k)-ql

              qauto(mg,k)=qauto(mg,k)+dql
              qlg(mg,k)=qlg(mg,k)-dql
              fluxa(mg,k)=dql*rhodz
            endif
          enddo
        enddo

c Or, using old autoconv scheme...
       case(0)
        do k=nl-1,1,-1
          do mg=1,ln2
            cfrain(mg,k)=0.0
            rhodz=rhoa(mg,k)*dz(mg,k)

            if(clfr(mg,k).gt.0.)then
       
              qcrit=(4*pi/3)*rhow*Rcm**3*cdrop(mg,k)/rhoa(mg,k)
              qcic=qlg(mg,k)/clfr(mg,k) !In cloud value

              if(qcic.lt.qcrit)then
                ql=qlg(mg,k)
                dqlo=0.
              else
                Crate=Aurate*
     &                rhoa(mg,k)*(rhoa(mg,k)/(cdrop(mg,k)*rhow))**(1./3)
                ql1=1./pow75(qcic**(-4./3)+(4./3)*Crate*delt)
                ql1=max(ql1, qcrit) !Intermediate qlg after auto
                Frb=dz(mg,k)*rhoa(mg,k)*(qcic-ql1)/delt
                cdt=delt*0.5*Ecol*0.24*pow75(Frb) !old
                selfcoll=min(ql1,ql1*cdt)
c               selfcoll=min(ql1,ql1*cdt/(1+0.5*cdt))
                ql2=ql1-selfcoll
                ql=clfr(mg,k)*ql2
                cfrain(mg,k)=clfr(mg,k)
              endif
              dql=qlg(mg,k)-ql

              qauto(mg,k)=qauto(mg,k)+dql
              qlg(mg,k)=qlg(mg,k)-dql
              fluxa(mg,k)=dql*rhodz
            endif
          enddo
        enddo
       case default
        write(6,*) "ERROR: Unknown autoconvmeth ",autoconvmeth
        stop
      end select

c Call frozen precipitation routine

      call icefall(lg,tdt,rhoa,dz,prf, !Inputs
     &             ttg,qsg,qlg,qfg,qtg,cfrac,cfmelt,            !In and Out
     &             fluxi,fluxm,clfr,cifr,qsubl,qaccr,pfstay,pqfsed,
     &             slopes)            !Outputs

      ! Set up prognostic rain - MJT
      ! The following is based on LDR flux divergence calculation
      ! (see icefall.f).  LDR's original scheme can be recovered
      ! by setting fout=1 and fthru=1.

      vr=0.1
      do k=nl-1,1,-1
        do mg=1,ln2
          rhodz=rhoa(mg,k)*dz(mg,k)
c Add flux of rain due to autoconversion to qrg
          qrg(mg,k)=qrg(mg,k)+fluxa(mg,k)/rhodz
          rhor(mg,k)=qrg(mg,k)*rhoa(mg,k)
          cfrain(mg,k)=max(cfrain(mg,k),cffall(mg,k)) ! max overlap autoconversion and rain from previous time step
          cftemp=cfrain(mg,k)+cfmelt(mg,k)-cfrain(mg,k)*cfmelt(mg,k)
          qrgtemp=qrg(mg,k)+fluxm(mg,k)/rhodz/real(njumps)
          vr(mg,k)=vr(mg,k+1)
          if (cftemp>0.) then
            mixrain=qrgtemp/max(cftemp,1.e-6)
            !Vr=11.3*Fr**(1./9.)/sqrt(rhoa(mg,k))  !Actual fall speed
            !vr(mg,k)=max(15.3*mixrain**0.125/(rhoa(mg,k)**0.4375), ! This version is equal to the above line from LDR
            vr(mg,k)=max(15.3*mixrain**0.125/sqrt(rhoa(mg,k)),
     &                   0.1)
          end if
          ! Set up the parameters for the flux-divergence calculation
          alph=delt*vr(mg,k)/dz(mg,k)
          fout(mg,k)=1.-exp(-alph) !analytical
          fthru(mg,k)=1.-fout(mg,k)/alph !analytical
        end do
      end do

      ! The following has been modified to track the random overlap rain fraction (rdclfr)
      ! and the max overlap rain fraction (mxclfr) so than both random overlaped and
      ! max/random overlaped clouds are supported - MJT

c      do nt=1,njumps

        do mg=1,ln2
          clfra(mg)=1.e-6
          ccra(mg)=0.
          fluxrain(mg)=0.
          prscav(mg,1)=0.
          mxclfr(mg)=0. ! max overlap rain fraction
          rdclfr(mg)=0. ! rnd overlap rain fraction
        enddo
        
c Now work down through the levels...
        
        do k=nl-1,1,-1
          do mg=1,ln2
            rhodz=rhoa(mg,k)*dz(mg,k)
            evap=0.

            ! The following flag detects max/random overlap clouds
            ! that are separated by a clear layer
            if ((clfr(mg,k).eq.0..and.cfrain(mg,k).eq.0.)
     &          .or.nmr.eq.0) then
              ! combine max overlap from last cloud with net random overlap
              rdclfr(mg)=rdclfr(mg)+mxclfr(mg)-rdclfr(mg)*mxclfr(mg)
              mxclfr(mg)=0.
            end if

c Add flux of melted snow to fluxrain

            fluxrain(mg)=fluxrain(mg)+fluxm(mg,k)/real(njumps)
            
c Evaporation of rain

            qpf=fluxrain(mg)/rhodz !Mix ratio of rain which falls into layer
            clrevap=(1.-clfr(mg,k))*qpf
            if(fluxrain(mg).gt.0.and.cfrac(mg,k).lt.1.0)then
              pk=100.0*prf(mg,k)
              qsg(mg,k)=qsati(pk,ttg(mg,k))
              if(ttg(mg,k).lt.tfrz.and.ttg(mg,k).ge.tice)then
                qsl(mg,k)=qsg(mg,k)+epsil
     &             *esdiff(nint(ttg(mg,k)-tfrz))/(100.0*prf(mg,k))
              else
                qsl(mg,k)=qsg(mg,k)
              endif             !qsl is qs value over liquid surface
              Tk=ttg(mg,k)
              es=qsl(mg,k)*pk/epsil 
              Apr=(hl/(rKa*Tk))*(hl/(rvap*Tk)-1)
              Bpr=rvap*Tk/((Dva/pk)*es)
              Fr=fluxrain(mg)/delt/clfra(mg)
              Cev=clfra(mg)
     &       *3.8e2*sqrt(Fr/rhoa(mg,k))/(qsl(mg,k)*(Apr+Bpr))
              dqsdt=hl*qsl(mg,k)/(rvap*ttg(mg,k)**2)
              bl=1+0.5*Cev*delt*(1+hlcp*dqsdt)
              evap=delt*(Cev/bl)*(qsl(mg,k)-qtg(mg,k))
              satevap=(qsl(mg,k)-qtg(mg,k))/(1+hlcp*dqsdt) !Evap to saturate
c              Vr=11.3*Fr**(1./9.)/sqrt(rhoa(mg,k))  !Actual fall speed
c              Vr=5./sqrt(rhoa(mg,k))                !Nominal fall speed

              evap=max(0., min(evap,satevap,qpf,clrevap))
              if(nevapls==-1)evap=0.
              if(nevapls==-2.and.k<=ktsav(mg).
     &                       and.condx(mg)>0.)evap=0.
              if(nevapls==-3)evap=.5*evap
              if(nevapls==-4.and.k<=ktsav(mg).
     &                       and.condx(mg)>0.)evap=.5*evap ! usual
              qevap(mg,k)=qevap(mg,k)+evap
              qtg(mg,k)=qtg(mg,k)+evap
              ttg(mg,k)=ttg(mg,k)-hlcp*evap
            endif
            frclr(mg,k)=rhodz*(clrevap-evap)             !over delt
            

c Now do the collection term.

c            if(qlg(mg,k).gt.1.e-10)then
              if(fluxrain(mg).gt.0.)then
                Fr=fluxrain(mg)/clfra(mg)/delt
                cfrain(mg,k)=max(cfrain(mg,k),
     &                       min(mxclfr(mg),clfr(mg,k))                ! max overlap
     &                      +rdclfr(mg)*max(clfr(mg,k)-mxclfr(mg),0.)) ! rnd overlap
              else
                Fr=0.
              endif

              if(fluxc(mg,k+1).gt.0.)then
                Frc=max(0.,fluxc(mg,k+1)/max(ccra(mg),0.01)/tdt) ! over tdt
                cfrain(mg,k)=max(cfrain(mg,k),clfr(mg,k)*ccra(mg))       ! rnd overlap
                !cfrain(mg,k)=max(cfrain(mg,k),min(clfr(mg,k),ccra(mg))) ! max overlap
              else
                Frc=0.
              endif

c The collection term comprises collection by stratiform rain falling from
c above (Fr), stratiform rain released in this grid box (Frb), and
c convective rain (Frc).
c Frb term now done above.

              fcol=min(1.,mxclfr(mg)/(1.e-20+clfr(mg,k))) !max overlap
              fcol=fcol+rdclfr(mg)-fcol*rdclfr(mg)        !rnd overlap
              cdt=delt*Ecol*0.24*(fcol*pow75(Fr)
     &                           +ccra(mg)*pow75(Frc))
c              prscav(mg,nlp-k)=cdt/Ecol !Inc conv part
              prscav(mg,nlp-k)=delt*0.24*fcol*pow75(Fr) !Strat only

              coll=min(qlg(mg,k),qlg(mg,k)*cdt/(1.+0.5*cdt))
              qcoll(mg,k)=qcoll(mg,k)+coll
              qlg(mg,k)=qlg(mg,k)-coll
              fluxrain(mg)=fluxrain(mg)+coll*rhodz
c            endif

c subtract evaporated rain

            fluxrain(mg)=fluxrain(mg)-rhodz*evap
            fluxrain(mg)=max(fluxrain(mg),0.) !To avoid roundoff -ve's

c Calculate the raining cloud cover down to this level, for stratiform (clfra)
c and convective (ccra).

            cfrain(mg,k)=min(1.,
     &             cfrain(mg,k)+cfmelt(mg,k)-cfrain(mg,k)*cfmelt(mg,k))
            if (frclr(mg,k).lt.1.e-15) then
              rdclfr(mg)=0.
              mxclfr(mg)=0.
            end if
            mxclfr(mg)=max(mxclfr(mg),cfrain(mg,k)) !max overlap
            clfra(mg)=max(1.e-15,
     &                rdclfr(mg)+mxclfr(mg)-rdclfr(mg)*mxclfr(mg)) !rnd overlap the mx and rd rain fractions
            ccra(mg)=max(ccra(mg),ccrain(mg,k)) !always max overlap for convective rainfall - MJT
            fracr(mg,k)=clfra(mg)

c Compute fluxes into the box
            cffluxin=clfra(mg)-cfrain(mg,k)
            rhorin=fluxrain(mg)/dz(mg,k)

c Compute the fluxes of rain leaving the box
c Use the flux-divergent form as in Rotstayn (QJRMS, 1997)
            cffluxout=cfrain(mg,k)*fout(mg,k)
            rhorout=rhor(mg,k)*fout(mg,k)
            
c Update the rhor and cffall fields
            cffall(mg,k)=cfrain(mg,k)-cffluxout
     &                  +cffluxin*(1.-fthru(mg,k))
            rhor(mg,k)=rhor(mg,k)-rhorout+rhorin*(1.-fthru(mg,k))
            fluxrain(mg)=rhorout*dz(mg,k)+fluxrain(mg)*fthru(mg,k) 
c Now fluxrain is flux leaving layer k
            fluxr(mg,k)=fluxr(mg,k)+fluxrain(mg)
            
          enddo
        enddo

c      enddo ! njumps
        
c Re-create qrg field

      qrg(1:ifull,nl)=0.
      do k=1,nl-1
        do mg=1,ln2
          qrg(mg,k)=max(rhor(mg,k)/rhoa(mg,k),0.)
          !if(qrg(mg,k)<0.)then
          !  print *,'k,mg,lg,qrg_rhor ',k,mg,qrg(mg,k)
          !  stop
          !endif
        enddo
      enddo

c Factor 0.5 here accounts for leapfrog scheme

      do mg=1,ln2
        precs(mg)=precs(mg)+0.5*(fluxr(mg,1)+fluxi(mg,1))
        preci(mg)=preci(mg)+0.5*fluxi(mg,1)
      enddo

c Remove small amounts of cloud

      do k=1,nl
        do mg=1,ln2
          if(qlg(mg,k).lt.1e-10.or.clfr(mg,k).lt.1e-5)then
            qtg(mg,k)=qtg(mg,k)+qlg(mg,k)
            ttg(mg,k)=ttg(mg,k)-hlcp*qlg(mg,k)
            qlg(mg,k)=0.
            clfr(mg,k)=0.
          endif
          if(qfg(mg,k).lt.1e-10.or.cifr(mg,k).lt.1e-5)then
            qtg(mg,k)=qtg(mg,k)+qfg(mg,k)
            ttg(mg,k)=ttg(mg,k)-hlscp*qfg(mg,k)
            qfg(mg,k)=0.
            cifr(mg,k)=0.
          endif
        enddo
      enddo

c      Adjust cloud fraction (and cloud cover) after precipitation
      if(nmaxpr==1.and.mydiag)then
        print *,'diags from newrain for idjd ',idjd
        write (6,"('cfrac ',9f8.3/6x,9f8.3)") cfrac(idjd,:)
        write (6,"('cftemp',9f8.3/6x,9f8.3)") cifr(idjd,:)+clfr(idjd,:)
        write (6,"('ccov_in',9f8.3/6x,9f8.3)") ccov(idjd,:)
      endif

c     from 16/1/06, ccov cfrac NOT UPDATED here
c      do k=1,nl-1
c        do mg=1,ln2
cc Next 3 lines commented for m35-m38 runs.
cC***          if(qlgsav(mg,k).gt.0.)then
cC***            clfr(mg,k)=clfr(mg,k)*qlg(mg,k)/qlgsav(mg,k)
cC***          endif
c          cftemp=min(cifr(mg,k)+clfr(mg,k), 1.)
c          if(cfrac(mg,k).gt.0.)then
c            ccov(mg,k)=min(1., ccov(mg,k)*cftemp/cfrac(mg,k))
c          else
c            ccov(mg,k)=cftemp
c          endif
c          cfrac(mg,k)=cftemp
c        enddo
c      enddo

      if(ntest==2)then
        do k=1,nl
          do mg=1,ln2
            if(cfrac(mg,k).gt.0)then
              if(qlg(mg,k).le.0.and.qfg(mg,k).le.0)then
                print*,'end newrain cloud with zero/neg water: cfrac=',
     &               cfrac(mg,k),' ccov= ',ccov(mg,k)
                print*,'cifr ,clfr ',cifr(mg,k),clfr(mg,k)
                print*,'qfg, qlg=', qfg(mg,k),qlg(mg,k)
                ns=(mg-1+lon)/lon
                mb=mg-(ns-1)*lon
                print*,'lg,ns,mg,k ',lg,ns,mb,k
                print*
              endif
            else
              if(qlg(mg,k).gt.0.or.qfg(mg,k).gt.0)then
               print*,'end newrain cloud water with no cloud: qfg,qlg=',
     &               qfg(mg,k),qlg(mg,k)
               print*,'cifr ,clfr ',cifr(mg,k),clfr(mg,k)
               ns=(mg-1+lon)/lon
               mb=mg-(ns-1)*lon
               print*,'ns,lg,k,mg ',ns,lg,k,mb
               print*
              endif
            endif
          enddo
        enddo
      endif  ! (ntest==2)

c Diagnostics for debugging

        if(ntest>0.and.mydiag)then
          mg=idjd
          write(25,'(a,3i3)')'IPASS=1, after newrain  ktau= ',ktau
          write(25,91)'ttg ',(ttg(mg,k),k=1,nl)
          write(25,9)'qtg ',(qtg(mg,k),k=1,nl)
          write(25,9)'qsg ',(qsg(mg,k),k=1,nl)
c          write(25,9)'qsl ',(qsl(mg,k),k=1,nl)
          write(25,1)'precs ',precs(mg)
          write(25,9)'cfrac ',(cfrac(mg,k),k=1,nl)
          write(25,9)'qlg ',(qlg(mg,k),k=1,nl)
          write(25,9)'qfg ',(qfg(mg,k),k=1,nl)
          write(25,9)'cdrop ',(cdrop(mg,k),k=1,nl)
          write(25,9)'cfa ',(cfa(mg,k),k=1,nl)
          write(25,9)'qca ',(qca(mg,k),k=1,nl)
          write(25,9)'fluxr ',(fluxr(mg,k),k=1,nl)
          write(25,9)'frclr ',(frclr(mg,k),k=1,nl)
          write(25,9)'fluxi ',(fluxi(mg,k),k=1,nl)
          write(25,9)'fluxc ',(fluxc(mg,k),k=1,nl)
          write(25,9)'qcoll ',(qcoll(mg,k),k=1,nl)
          write(25,9)'qauto ',(qauto(mg,k),k=1,nl)
          write(25,9)'qaccr ',(qaccr(mg,k),k=1,nl)
          write(25,9)'qevap ',(qevap(mg,k),k=1,nl)
          write(25,9)'ccrain ',(ccrain(mg,k),k=1,nl-1)
          write(25,9)'cfrain ',(cfrain(mg,k),k=1,nl-1)
          write(25,9)'fracr ',(fracr(mg,k),k=1,nl-1)
          write(25,*)
        endif  ! (ntest>0)
 1    format(3(a,f10.5))
 91   format(a,30f10.3)
 9    format(a,30g10.3)


      return

      end
