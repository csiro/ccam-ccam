c This routine is part of the prognostic cloud scheme. It calculates the frozen 
c precipitation. It is called by newrain.
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
c
c      lg - latitude index (ranges from 1 at poles to LAT at equator)
c      tdt - leapfrog timestep (seconds)
c      rhoa - air density (kg/m**3)
c      dz - layer thicknes (m)
c      prf - pressure at full levels (in hPa. NB: not SI units)
c
c In/Out:
c
c from arguments
c
c      ttg - temperature (K)
c      qsg - saturation mixing ratio (kg/kg)
c      qlg - cloud liquid water mixing ratio (kg/kg)
c      qfg - cloud ice mixing ratio (kg/kg)
c      qtg - water vapour mixing ratio (kg/kg) - called qg in C-CAM
c      cfrac - stratiform cloud fraction
c      cfmelt - fraction of grid box occupied by falling ice that melts
c
c Output:
c
c from arguments
c
c      fluxi - flux of falling ice in timestep (kg/m**2)
c      fluxm - flux of falling ice that melts in timestep (kg/m**2)
c      clfr - liquid cloud fraction
c      qsubl - sublimation of snowfall (kg/kg)
c      qaccr - accretion by snow of cloud liquid water (kg/kg)
c      pqfsed - (dqf/qf) due to ice falling out of layer (fraction)
c      pfstay - incoming flux of ice staying in layer  (kg/m**2/s)
c
c******************************************************************************

      subroutine icefall (lg,tdt,rhoa,dz,prf,               !Inputs
     &                   ttg,qsg,qlg,qfg,qtg,cfrac,cfmelt,  !In and Out
     &                   fluxi,fluxm,clfr,cifr,qsubl,qaccr, !Outputs
     &                   pfstay,pqfsed,slopes)              !Outputs

      use cc_mpi, only : mydiag
      use kuocomb_m
      use morepbl_m  !condx  
      implicit none
C Global parameters
      include 'newmpar.h'
      include 'const_phys.h' !Input physical constants
      include 'cparams.h'    !Input cloud scheme parameters
      include 'kuocom.h'     !ktsav,ldr,nevapls
      include 'parm.h'       !just for nmaxpr
      include 'params.h'     !Input model grid dimensions (modified PARAMS.f for CCAM)

C Argument list
      integer lg
      real tdt
      real rhoa(ln2,nl)
      real dz(ln2,nl)
      real prf(ln2,nl)
      real ttg(ln2,nl)
      real qsg(ln2,nl)
      real qlg(ln2,nl)
      real qfg(ln2,nl)
      real qtg(ln2,nl)
      real cfrac(ln2,nl)
      real cfmelt(ln2,nl)
      real fluxi(ln2,nl)
      real fluxm(ln2,nl)
      real clfr(ln2,nl)
      real cifr(ln2,nl)
      real qsubl(ln2,nl)
      real qaccr(ln2,nl)
      real pqfsed(ln2,nl)
      real pfstay(ln2,nl)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
c      include 'FEWFLAGS.f'   !Input debug, lgdebug etc.

C Local work arrays and variables
c     real alpha1(ln2,nl)
      real Csbsav(ln2,nl-1)
      real cifra(ln2)
      real fluxice(ln2)
      real fout(ln2,nl-1)
      real fthru(ln2,nl-1)
      real gam(ln2,nl-1)
c     real qacci(ln2,nl)
c     real qaggi(ln2,nl)
      real qfdiv(ln2,nl)
      real rhoi(ln2,nl)
      real rica(ln2)
      real slopes(ln2,nl)
      real vi2(ln2,nl)
      real mxclfr(ln2)
      real rdclfr(ln2)

      integer k,mg,ntest

      real alph
      real alphaf
      real aprpr
      real bf
      real bprpr
      real caccr
      real cdt
      real cffluxin
      real cffluxout
      real ci
      real csb
      real curly
      real delt
      real dqf
      real dqs
      real dttg
      real es
      real fsclr
      real pk
      real qif
      real ql
      real rhoiin
      real rhoiout
      real sublflux
      real tc
      real tk
      parameter (ntest=0)  ! 0=off, 1=on
 
C Start code : ----------------------------------------------------------

c Set up timestep for ice processes

      delt=tdt
c      delt=360.
c      njumps=nint(tdt/delt)
c      delt=tdt/njumps !To make sure tdt it a multiple of delt

c Convert from mixing ratio to density of ice, and work out ice cloud fraction

      do k=1,nl
        do mg=1,ln2
          fluxi(mg,k)=0.
!         N.B. qfg >= 0 from dynamics, but occasionally it seems to become -ve,
!         apparently from newrain.f, or else lower down in this subroutine
          qfg(mg,k)=max( 0.,qfg(mg,k) )  
          rhoi(mg,k)=rhoa(mg,k)*qfg(mg,k) 
          cifr(mg,k)=cfrac(mg,k)*qfg(mg,k)/
     &                           max(qlg(mg,k)+qfg(mg,k),1.e-30)
          clfr(mg,k)=max(cfrac(mg,k)-cifr(mg,k),0.)
          qsubl(mg,k)=0.
          qaccr(mg,k)=0.
c          qaggi(mg,k)=0.
c          qacci(mg,k)=0.
          qfdiv(mg,k)=0.
          cfmelt(mg,k)=0.
        enddo
      enddo

      if(nmaxpr==1.and.mydiag)then
        print *,'diags from icefall for idjd ',idjd
        write (6,"('clfra ',9f8.3/6x,9f8.3)") clfr(idjd,:)
      endif
      if(ntest>0.and.mydiag)then
          mg=idjd
          write(25,'(a,3i3)')'IPASS=1, before icefall.'
          write(25,91)'cfrac ',(cfrac(mg,k),k=1,nl)
          write(25,91)'cifr ',(cifr(mg,k),k=1,nl)
          write(25,91)'clfr ',(clfr(mg,k),k=1,nl)
          write(25,9)'qlg ',(qlg(mg,k),k=1,nl)
          write(25,9)'qfg ',(qfg(mg,k),k=1,nl)
          write(25,*)
      endif  ! (ntest>0)
      if(ntest==2.and.mydiag)then
        do k=1,nl
          do mg=1,ln2
            if(cfrac(mg,k)>0)then
              if(qlg(mg,k)==0.and.qfg(mg,k)==0)then
                print*,
     &         'start icefall cloud with zero water: cfrac=',cfrac(mg,k)
              endif
            else
              if(qlg(mg,k)>0.or.qfg(mg,k)>0)then
             print*,'start icefall cloud water with no cloud: qfg,qlg=',
     &               qfg(mg,k),qlg(mg,k)
              endif
            endif
          enddo
        enddo
      endif  ! (ntest==2)

c Set up ice fall speed field and other arrays

      do mg=1,ln2
        vi2(mg,nl)=0.1
        slopes(mg,nl)=0.
      enddo

      if(abs(ldr)==1)then  ! 1 for R21 runs, like prev lw=22
        do k=nl-1,1,-1
          do mg=1,ln2
            vi2(mg,k)=vi2(mg,k+1)
            if(cifr(mg,k)>0.)then
              vi2(mg,k)=3.23*(rhoi(mg,k)/cifr(mg,k))**0.17
            endif
          enddo
        enddo
      endif   ! (abs(ldr)==1)

      if(abs(ldr)==2)then  
        do k=nl-1,1,-1
          do mg=1,ln2
            vi2(mg,k)=vi2(mg,k+1)
            if(cifr(mg,k)>0.)then
              vi2(mg,k)=0.9*3.23*(rhoi(mg,k)/cifr(mg,k))**0.17
            endif
          enddo
        enddo
      endif   ! (abs(ldr)==2)

      if(abs(ldr)==3)then  
        do k=nl-1,1,-1
          do mg=1,ln2
            vi2(mg,k)=vi2(mg,k+1)
            if(cifr(mg,k)>0.)then
              vi2(mg,k)=max(0.1,2.05+0.35*log10(qfg(mg,k)/cifr(mg,k)))
            endif
          enddo
        enddo
      endif   ! (abs(ldr)==3)

      if(abs(ldr)==4)then  
        do k=nl-1,1,-1
          do mg=1,ln2
            vi2(mg,k)=vi2(mg,k+1)
            if(cifr(mg,k)>0.)then
              vi2(mg,k)=1.4*3.23*(rhoi(mg,k)/cifr(mg,k))**0.17
            endif
          enddo
        enddo
      endif   ! (abs(ldr)==4)

c     following are alternative slightly-different versions of above
c     used for I runs from 29/4/05 till 30/8/05
c     for given qfg, large cifr implies small ice crystals, 
c     with a small fall speed. 
c     Note that for very small qfg, cifr is small.
c     But rhoi is like qfg, so ratio should also be small and OK.
      if(abs(ldr)==11)then  ! 1 for R21 runs, like prev lw=22
        do k=nl-1,1,-1
          do mg=1,ln2
            vi2(mg,k)=max( vi2(mg,k+1),3.23*(rhoi(mg,k)/
     &                                max(cifr(mg,k),1.e-30))**0.17 )
          enddo
        enddo
      endif   ! (abs(ldr)==11)

      if(abs(ldr)==22)then  
        do k=nl-1,1,-1
          do mg=1,ln2
            vi2(mg,k)=max( vi2(mg,k+1),.9*3.23*(rhoi(mg,k)/
     &                                max(cifr(mg,k),1.e-30))**0.17 )
          enddo
        enddo
      endif   ! (abs(ldr)==22)

      if(abs(ldr)==33)then  
        do k=nl-1,1,-1
         do mg=1,ln2
c          following max gives vi2=.1 for qfg=cifr=0
           vi2(mg,k)=max( vi2(mg,k+1),2.05 +0.35*
     &      log10(max(qfg(mg,k),2.68e-36)/max(cifr(mg,k),1.e-30)) )
          enddo
        enddo
      endif   ! (abs(ldr)==33)

      do k=nl-1,1,-1
        do mg=1,ln2
          tc=ttg(mg,k)-tfrz
          slopes(mg,k)=1.6e3*10**(-0.023*tc)
          alphaf=hls*qsg(mg,k)/(rvap*ttg(mg,k)**2)
          gam(mg,k)=hlscp*alphaf !(L/cp)*dqsdt (HBG notation)

c Set up the Rate constant for snow sublimation
          
          Tk=ttg(mg,k)
          pk=100*prf(mg,k)
          es=qsg(mg,k)*pk/epsil
          Aprpr=(hls/(rKa*Tk))*(hls/(rvap*Tk)-1)
          Bprpr=rvap*Tk/((Dva/pk)*es)
          curly=0.65*slopes(mg,k)**2+0.493*slopes(mg,k)!Factor in curly brackets
     &         *sqrt(slopes(mg,k)*vi2(mg,k+1)*rhoa(mg,k)/um)
          if(nevapls==-1)curly=0.
          if(nevapls==-2.and.condx(mg)>0..and.k<=ktsav(mg))curly=0.

c Define the rate constant for sublimation of snow, omitting factor rhoi

          Csbsav(mg,k)=4*curly/
     &        (rhoa(mg,k)*qsg(mg,k)*(Aprpr+Bprpr)*pi*vi2(mg,k+1)*rhosno)

c Set up the parameters for the flux-divergence calculation

          alph=delt*vi2(mg,k)/dz(mg,k)
          fout(mg,k)=1.-exp(-alph) !analytical
          fthru(mg,k)=1.-fout(mg,k)/alph !analytical
c         alpha1(mg,k)=1.e-3*exp(0.025*(ttg(mg,k)-tfrz))
        enddo
      enddo

c Save sedimentation rate for aerosol scheme

      do k=1,nl-1
        do mg=1,ln2
c         pqfsed(mg,nlp-k)=fout(mg,k)*qfg(mg,k)/tdt
          pqfsed(mg,nlp-k)=fout(mg,k)
        enddo
      enddo
      do mg=1,ln2
        pqfsed(mg,1)=0.
      enddo

      ! The following has been modified to track the random overlap rain fraction (rdclfr)
      ! and the max overlap rain fraction (mxclfr) so than both random overlaped and
      ! max/random overlaped clouds are supported - MJT

c      do nt=1,njumps !Need njumps=1 for new code

c Assume no cloud at top level

        do mg=1,ln2
          fluxice(mg)=0.
          cifra(mg)=0.
          rica(mg)=0.
          pfstay(mg,1)=0.
          mxclfr(mg)=0. ! max overlap rain fraction
          rdclfr(mg)=0. ! rnd overlap rain fraction
        enddo

c Now work down through the levels...

        do k=nl-1,1,-1
          do mg=1,ln2

            sublflux=0.
            fsclr=0.
            caccr=0.
            dqf=0.

            ! The following flag detects max/random overlap clouds
            ! that are separated by a clear layer
            if (cifr(mg,k).eq.0..or.nmr.eq.0) then
              ! combine max overlap from last cloud with net random overlap
              rdclfr(mg)=rdclfr(mg)+mxclfr(mg)-rdclfr(mg)*mxclfr(mg)
              mxclfr(mg)=0.
            end if

c Melt falling ice if > 0 deg C

            if(ttg(mg,k)>tfrz.and.fluxice(mg)>0.)then
              qif=fluxice(mg)/(dz(mg,k)*rhoa(mg,k))      !Mixing ratio of ice
              dttg=-hlfcp*qif
              ttg(mg,k)=ttg(mg,k)+dttg
              qsg(mg,k)=qsg(mg,k)+gam(mg,k)*dttg/hlscp
              fluxm(mg,k)=fluxm(mg,k)+fluxice(mg)
              cfmelt(mg,k)=cifra(mg)
              fluxice(mg)=0.
              cifra(mg)=0.
              rdclfr(mg)=0.
              mxclfr(mg)=0.
            endif


c Compute the sublimation of ice falling from level k+1 into level k

            fsclr=(1.-cifr(mg,k)-clfr(mg,k))*fluxice(mg)
            if(fluxice(mg)>0.and.qtg(mg,k)<qsg(mg,k))then ! sublime snow

              Csb=Csbsav(mg,k)*fluxice(mg)/delt !LDR
              bf=1.+0.5*Csb*delt*(1.+gam(mg,k))
              dqs=max(0.,delt*(Csb/bf)*(qsg(mg,k)-qtg(mg,k)))
              dqs=min(dqs,(qsg(mg,k)-qtg(mg,k))/(1.+gam(mg,k))) !Don't supersat.

              sublflux=min(dqs*rhoa(mg,k)*dz(mg,k),fsclr)
              fluxice(mg)=fluxice(mg)-sublflux
              fsclr=fsclr-sublflux
              dqs=sublflux/(rhoa(mg,k)*dz(mg,k))
              qsubl(mg,k)=qsubl(mg,k)+dqs
              qtg(mg,k)=qtg(mg,k)+dqs
              dttg=-hlscp*dqs
              ttg(mg,k)=ttg(mg,k)+dttg
              qsg(mg,k)=qsg(mg,k)+gam(mg,k)*dttg/hlscp
            endif

c Save this for the wet deposition scheme.

            pfstay(mg,nlp-k)=fluxice(mg)*(1.-fthru(mg,k))/tdt !Flux staying in layer k

c Accretion of cloud water by falling ice
c This calculation uses the incoming fluxice without subtracting sublimation
c (since subl occurs only outside cloud), so add sublflux back to fluxice.
            
            if(fluxice(mg)+sublflux>0.and.qlg(mg,k)>0.)then
              ql=qlg(mg,k)
              cdt=Eac*slopes(mg,k)*(fluxice(mg)+sublflux)/(2.*rhosno)
              dqf=min(ql,cifra(mg)*ql,ql*cdt/(1.+0.5*cdt))

              clfr(mg,k)=clfr(mg,k)*(1.-dqf/qlg(mg,k))
              caccr=clfr(mg,k)*dqf/qlg(mg,k)
              qlg(mg,k)=qlg(mg,k)-dqf
              qaccr(mg,k)=qaccr(mg,k)+dqf
              fluxice(mg)=fluxice(mg)+rhoa(mg,k)*dz(mg,k)*dqf
              dttg=hlfcp*dqf
              ttg(mg,k)=ttg(mg,k)+dttg
              qsg(mg,k)=qsg(mg,k)+gam(mg,k)*dttg/hlscp
            endif


            if(fsclr.lt.1.E-15)then
              rdclfr(mg)=0.
              mxclfr(mg)=0.
            end if
            ci=cifr(mg,k)+caccr
            mxclfr(mg)=max(mxclfr(mg),ci) !max overlap
            cifra(mg)=max(0.01,
     &        mxclfr(mg)+rdclfr(mg)-mxclfr(mg)*rdclfr(mg)) !rnd overlap the mx and rd ice fractions


c Compute fluxes into the box
            
            if(fluxice(mg)>0.)then
              rhoiin=fluxice(mg)/dz(mg,k)
              cffluxin=min(1.,rhoiin/rica(mg))
            else
              rhoiin=0.
              cffluxin=0.
            endif


c Compute the fluxes of ice and cloud amount leaving the box
            
            if(cifr(mg,k)>0.)then

c Use the flux-divergent form as in Rotstayn (QJRMS, 1997)
              rhoiout=rhoi(mg,k)*fout(mg,k)
              cffluxout=cifr(mg,k)*fout(mg,k)

c Or, use a Kessler type formulation with a fixed qcrit option and 
c     accretion of ice
C***              riic=rhoi(mg,k)/cifr(mg,k)
C***c              rcrit=1.e-6*rhoa(mg,k) !i.e. qcrit=const
C***c or use qcrit = f(T) following Sinha & Shine, J. Climate (1994)
C***              rcrit=0.0007e-3*exp(0.0987*(ttg(mg,k)-213.15))
C***              cdt=delt*alpha1(mg,k)
C***c              dra=min(rhoi(mg,k), cifr(mg,k)*dim(riic,rcrit)*cdt) !Explicit
C***              dra=cifr(mg,k)*dim(riic,rcrit)*min(1., cdt/(1+0.5*cdt)) !T-centred
C***c              dra=cifr(mg,k)*dim(riic,rcrit)*(1.-exp(-cdt))      !Analytical
C***c              qaggi(mg,k)=qaggi(mg,k)+dra   !diagnostic only
C***              rhoi2=rhoi(mg,k)-dra !Intermediate rhoi after aggregation
C***
C***c Accretion of ice by falling snow
C***c Exclude accreted LW (graupel) from flxi, but include 50% of accreted ice.
C***
C***              Esi=1.e3*alpha1(mg,k) !Efficiency for accretion of ice by snow
C***              flxi=fluxice(mg)+sublflux+rhoa(mg,k)*dz(mg,k)*(.5*dra-dqf)
C***              cdt=Esi*slopes(mg,k)*flxi/(2*rhosno)
C***c              drs=min(rhoi2,rhoi2*cdt) !Explicit
C***              drs=rhoi2*cdt/(1+0.5*cdt) !Time-centered
C***c              drs=rhoi2*(1.-exp(-cdt))   !analytical
C***c              drs=0. !!!!!
C***c              qacci(mg,k)=qacci(mg,k)+drs !diagnostic only
C***
C***              rhoiout=min(rhoi(mg,k), dra+drs)
C***              cffluxout=cifr(mg,k)*max(0.,rhoiout/rhoi(mg,k))

c End of Kessler-type formulation 

              rica(mg)=rhoi(mg,k)/cifr(mg,k) !in cloud rhoi above
            else !Keep value of rica from above
              rhoiout=0.
              cffluxout=0.
            endif
            
c Update the rhoi and cifr fields

            cifr(mg,k)=min(1.-clfr(mg,k),
     &                (cifr(mg,k)-cffluxout)+cffluxin*(1-fthru(mg,k)))
            rhoi(mg,k)=(rhoi(mg,k)-rhoiout)+rhoiin*(1-fthru(mg,k))
            fluxice(mg)=rhoiout*dz(mg,k)+fluxice(mg)*fthru(mg,k) 
c Now fluxice is flux leaving layer k
            fluxi(mg,k)=fluxi(mg,k)+fluxice(mg)

          enddo
        enddo

c      enddo  ! njumps loop

c End of small timestep loop

c Re-create qfg field

      do k=1,nl
        do mg=1,ln2
          qfg(mg,k)=rhoi(mg,k)/rhoa(mg,k)
          if(qfg(mg,k)<0.)then
            print *,'k,mg,lg,qfg_rhoi ',k,mg,lg,qfg(mg,k)
            stop
          endif
        enddo
      enddo

c Diagnostics for debugging
      if(nmaxpr==1.and.mydiag)then
        write (6,"('vi2 ',9f8.3/4x,9f8.3)") vi2(idjd,:)
        write (6,"('cfraci',9f8.3/6x,9f8.3)") cfrac(idjd,:)
        write (6,"('cifr  ',9f8.3/6x,9f8.3)") cifr(idjd,:)
        write (6,"('clfrb ',9f8.3/6x,9f8.3)") clfr(idjd,:)
        write (6,"('qfg   ',3p9f8.3/6x,9f8.3)") qfg(idjd,:)
      endif

      if(ntest>0.and.mydiag)then
          mg=idjd
          write(25,'(a,3i3)')'IPASS=1, after icefall.'
          write(25,91)'cfrac ',(cfrac(mg,k),k=1,nl)
          write(25,91)'cifr ',(cifr(mg,k),k=1,nl)
          write(25,91)'clfr ',(clfr(mg,k),k=1,nl)
          write(25,9)'qlg ',(qlg(mg,k),k=1,nl)
          write(25,9)'qfg ',(qfg(mg,k),k=1,nl)
          write(25,9)'qsubl ',(qsubl(mg,k),k=1,nl)
          write(25,9)'vi2 ',(vi2(mg,k),k=1,nl)
          write(25,9)'dz ',(dz(mg,k),k=1,nl)
          write(25,9)'rhoa ',(rhoa(mg,k),k=1,nl)
          write(25,9)'fluxi ',(fluxi(mg,k),k=1,nl)
c          write(25,9)'qfdiv ',(qfdiv(mg,k),k=1,nl)
          write(25,*)

c          ma=(ns-1)*lon
C***          do k=1,nl
C***            do mg=1+ma,lon+ma
C***              mb=mg-ma
C***              if(qfg(mg,k)>0.1)then
C***                print*,'end icefall qfg = ',qfg(mg,k)
C***                print*,'lg,ns,mg,k ',lg,ns,mb,k
C***              endif
C***              
C***              if(cfrac(mg,k)>0)then
C***                if(qlg(mg,k)<=0.and.qfg(mg,k)<=0)then
C***                  print*,'end icefall cloud with zero/neg water:cfrac=',
C***     &                 cfrac(mg,k)
C***                  print*,'cifr ,clfr ',cifr(mg,k),clfr(mg,k)
C***                  print*,'qfg, qlg=', qfg(mg,k),qlg(mg,k)
C***                  print*,'lg,ns,mg,k ',lg,ns,mb,k
C***                  print*
C***                endif
C***              else
C***                if(qlg(mg,k)>0.or.qfg(mg,k)>0)then
C***                  print*,'end icefall, cloud water/no cloud: qfg,qlg=',
C***     &                 qfg(mg,k),qlg(mg,k)
C***                  print*,'cifr ,clfr ',cifr(mg,k),clfr(mg,k)
C***                  print*,'ns,lg,k,mg ',ns,lg,k,mb
C***                  print*
C***                endif
C***              endif
C***            enddo
C***          enddo
      endif  ! (ntest>0)

 1    format(3(a,g10.3))
 91   format(a,30f10.3)
 9    format(a,30g10.3)

      return

      end
