c This routine is part of the prognostic cloud scheme. It calculates the
c formation and dissipation of cloud, and the liquid fraction in mixed-phase
c clouds. It is called by progcld.
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
c from common/fewflags in FEWFLAGS.f
c      debug - namelist flag to control single column debugging
c      lgdebug - latitude index for single column debugging
c      insdebug - hemisphere index for single column debugging
c      mgdebug  - longitude index for single column debugging
c
c see also include files physparams.h (model physical constants)
c                        cparams.h    (cloud scheme parameters)
c
c from arguments
c      tdt - leapfrog timestep (seconds)
c      lg - latitude index (ranges from 1 at poles to LAT at equator)
c      
c      land - logical variable for surface type ( = T for land points)
c      prf - pressure at full levels (in hPa. NB: not SI units)
c      kbase - k index of convective cloud base   *** not used
c      ktop - k index of convective cloud top     *** not used
c
c In/Out:
c
c from arguments
c      ttg - temperature (K)
c      qtg - water vapour mixing ratio (kg/kg)
c      qlg - cloud liquid water mixing ratio (kg/kg)
c      qfg - cloud ice mixing ratio (kg/kg)
c
c Output:
c
c from arguments
c      cfrac - cloudy fraction of grid box
c      ccov - cloud cover looking from above (currently = cloud fraction)
c 
c******************************************************************************

      subroutine newcloud(tdt,lg,land,prf,kbase,ktop,rhoa,cdso4,  !Inputs
     &                    ttg,qtg,qlg,qfg,                        !In and out
     &                    cfrac,ccov,cfa,qca)                     !Outputs

c This routine is part of the prognostic cloud water scheme

      implicit none
C Global parameters
      include 'newmpar.h'
      include 'const_phys.h' !Input physical constants
      include 'cparams.h'    !Input cloud scheme parameters
      include 'kuocom.h'     !Input cloud scheme parameters rcrit_l & rcrit_s
      include 'params.h'     !Input model grid dimensions (modified params.h for CCAM)
      include 'sigs.h'

C Argument list
      real tdt
      integer lg
      logical land(ln2)
      real prf(ln2,nl)
      integer kbase(ln2)
      integer ktop(ln2)
      real rhoa(ln2,nl)
      real ttg(ln2,nl)
      real qtg(ln2,nl)
      real qlg(ln2,nl)
      real qfg(ln2,nl)
      real cfrac(ln2,nl)
      real ccov(ln2,nl)
      real cdso4(ln2,nl)
      real cfa(ln2,nl)
      real qca(ln2,nl)


C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
c      include 'CNSTA.f'
c      include 'FEWFLAGS.f'   !Input debug, lgdebug etc.
c     include 'TIMEX.f'

C Local work arrays and variables
      real qcg(ln2,nl)
      real qtot(ln2,nl),tliq(ln2,nl),qsg(ln2,nl)
      real fice(ln2,nl)
      real qcold(ln2,nl)
      real Cdrop(ln2,nl)
      real rcrit(ln2,nl)

      integer k
      integer mg
      integer ns

      real al
      real alf
      real aprpr
      real bprpr
      real cice
      real cm0
      real crate
      real deles
      real delq
      real dqsdt
      real es
      real esl
      real fd
      real fi
      real fl
      real hlrvap
      real pk
      real qc
      real qfdep
      real qfnew
      real qi0
      real qs
      real qsl
      real rhoic
      real root6
      real root6i
      real tk

C Local data, functions etc
      logical ukconv
      integer naerosol_i(2)
      data  ukconv, naerosol_i / .false., 2*0 /
      logical debug
      integer lgdebug,mgdebug,insdebug
      data debug,lgdebug,mgdebug,insdebug /.false.,1,10106,1/
      save ukconv,naerosol_i,debug,lgdebug,mgdebug,insdebug

      real esdiff(-40:0)  !Table of es(liq) - es(ice) (MKS), -40 <= t <= 0 C
      data esdiff / 
     & 6.22, 6.76, 7.32, 7.92, 8.56, 9.23, 9.94,10.68,11.46,12.27,
     & 13.11,13.99,14.89,15.82,16.76,17.73,18.70,19.68,20.65,21.61,
     & 22.55,23.45,24.30,25.08,25.78,26.38,26.86,27.18,27.33,27.27,
     & 26.96,26.38,25.47,24.20,22.51,20.34,17.64,14.34,10.37, 5.65,
     & 0.08 /
!     include 'ESTABL.f'  !Contains arithmetic statement function qsat(p,T)
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
      qsati(pp,t) = epsil*estabi(t)/(pp-estabi(t)) !Usual formula

C Start code : ----------------------------------------------------------

      root6=sqrt(6.)
      root6i=1./root6

      if(debug)then
        if(lg.eq.lgdebug)then
          ns=insdebug
          mg=mgdebug+(ns-1)*lon
          write(25,'(a,3i3)')'IPASS=1, before newcloud'
          write(25,91)'prf ',(prf(mg,k),k=1,nl)
          write(25,91)'ttg ',(ttg(mg,k),k=1,nl)
          write(25,9)'qtg ',(qtg(mg,k),k=1,nl)
          write(25,9)'qlg ',(qlg(mg,k),k=1,nl)
          write(25,9)'qfg ',(qfg(mg,k),k=1,nl)
          write(25,*)
        endif
      endif

c Define Cdrop

      if(naerosol_i(2).gt.0)then
        do k=1,nl
          do mg=1,ln2
            Cdrop(mg,k)=cdso4(mg,k)
          enddo
        enddo
      else
        do mg=1,ln2
          if(land(mg))then
            Cdrop(mg,1)=Cdropl
          else
            Cdrop(mg,1)=Cdrops
          endif
        enddo
        do k=2,nl
          do mg=1,ln2
            Cdrop(mg,k)=Cdrop(mg,1)
          enddo
        enddo
      endif

c First melt cloud ice or freeze cloud water to give correct ice fraction fice.
c Then calculate the cloud conserved variables qtot and tliq.
c Note that qcg is the total cloud water (liquid+frozen)



      do k=1,nl
        do mg=1,ln2
          if(ttg(mg,k).ge.tfrz)then
            fice(mg,k)=0.
          elseif(ttg(mg,k).ge.tice)then
            if(qfg(mg,k).gt.1.0e-12)then
              fice(mg,k)=min(qfg(mg,k)/(qfg(mg,k)+qlg(mg,k)), 1.)
            else
              fice(mg,k)=0.
            endif
          else
            fice(mg,k)=1.
          endif
          qcg(mg,k)=qlg(mg,k)+qfg(mg,k)
          qcold(mg,k)=qcg(mg,k)
          qfnew=fice(mg,k)*qcg(mg,k)
          ttg(mg,k)=ttg(mg,k)+hlfcp*(qfnew-qfg(mg,k)) !Release L.H. of fusion
          qfg(mg,k)=fice(mg,k)*qcg(mg,k)
          qlg(mg,k)=max(0.,qcg(mg,k)-qfg(mg,k))
        enddo
      enddo

c Precompute the array of critical relative humidities 

      do k=1,nl
        do mg=1,ln2
          if(land(mg))then
c           rcrit(mg,k)=max(0.75,sig(k)**2)
!           rcrit(mg,k)=max(0.75,sig(k))
!           rcrit(mg,k)=max(0.75,min(.99,sig(k)))    ! same as T63
            rcrit(mg,k)=max(rcrit_l,min(.99,sig(k))) ! .75 for same as T63
          else
!           rcrit(mg,k)=max(0.85,sig(k)**2)
c           rcrit(mg,k)=max(0.9,sig(k))
!           rcrit(mg,k)=max(0.85,min(.99,sig(k)))    ! same as T63
            rcrit(mg,k)=max(rcrit_s,min(.99,sig(k))) ! .85 for same as T63
          endif
C***          if(k.eq.1)rcrit(mg,k)=0.999
C***          if(k.eq.2)rcrit(mg,k)=0.99
C***          if(k.eq.3)rcrit(mg,k)=0.95
        enddo
      enddo


c Calculate cloudy fraction of grid box (cfrac) and gridbox-mean cloud water
c using the triangular PDF of Smith (1990)

      do k=1,nl
        do mg=1,ln2
          fi=0.
          if(ttg(mg,k).le.Tice)fi=1.
c          fi=fice(mg,k)
          hlrvap=(hl+fi*hlf)/rvap
          qtot(mg,k)=qtg(mg,k)+qcg(mg,k)
          tliq(mg,k)=ttg(mg,k)-hlcp*qcg(mg,k)-hlfcp*qfg(mg,k)

c Calculate qs and gam=(L/cp)*dqsdt,  at temperature tliq
          pk=100.0*prf(mg,k)
          qsg(mg,k)=qsati(pk,tliq(mg,k)) !Ice value
          deles=esdiff(min(max(-40,(nint(tliq(mg,k)-tfrz))),0))
          qsl=qsg(mg,k)+epsil*deles/pk !qs over liquid
          qs=qsg(mg,k)
          if(ttg(mg,k).lt.tfrz.and.ttg(mg,k).gt.Tice)qs=qsl
          dqsdt=qs*hlrvap/tliq(mg,k)**2

          al=1./(1.+(hlcp+fi*hlfcp)*dqsdt)    !Smith's notation
          qc=qtot(mg,k)-qs

          delq=(1-rcrit(mg,k))*qs      !UKMO style (equivalent to above)
          cfrac(mg,k)=1.
          qcg(mg,k)=al*qc
          if(qc.lt.delq)then
            cfrac(mg,k)=1.-0.5*((qc-delq)/delq)**2
            qcg(mg,k)=al*(qc-(qc-delq)**3/(6.*delq**2))
          endif
          if(qc.le.0.)then
            cfrac(mg,k)=0.5*((qc+delq)/delq)**2
            qcg(mg,k)=al*(qc+delq)**3/(6.*delq**2)
          endif
          if(qc.le.-delq)then
            cfrac(mg,k)=0.
            qcg(mg,k)=0.
          endif

c Calculate the cloud fraction (cfa) in which ql exceeds qcrit, and
c the corresponding gridbox-mean cloud water mixing ratio qca. 
c This (qca) is the cloud-water mixing ratio inside cfa divided by cfa.
c The new variable qc2 is like qc above, but is used for integration limits
c only, not the integrand

C***          qcrit=(4*pi/3)*rhow*Rcm**3*Cdrop(mg,k)/rhoa(mg,k)
C***          qc2=qtot(mg,k)-qs-qcrit/al 
C***          cfa(mg,k)=1.
C***          qca(mg,k)=al*qc
C***          if(qc2.lt.delq)then
C***            cfa(mg,k)=1-0.5*((qc2-delq)/delq)**2
C***            qto=(qtot(mg,k)-delq+2*(qs+qcrit/al))/3.
C***            qca(mg,k)=al*(qtot(mg,k) - qto + cfa(mg,k)*(qto-qs))
C***          endif
C***          if(qc2.le.0.)then
C***            cfa(mg,k)=0.5*((qc2+delq)/delq)**2
C***            qca(mg,k)=cfa(mg,k)*(al/3.)*(2*qcrit/al + qc+delq)
C***          endif
C***          if(qc2.le.-delq)then
C***            cfa(mg,k)=0.
C***            qca(mg,k)=0.
C***          endif

        enddo
      enddo

c Condense or evaporate ql first.

      do k=1,nl
        do mg=1,ln2
          if(ttg(mg,k).gt.Tice)then
            qlg(mg,k) = max (0., qcg(mg,k)-qfg(mg,k))
            qfg(mg,k) = qcg(mg,k) - qlg(mg,k)
          else
            qlg(mg,k) = 0.
            qfg(mg,k) = qcg(mg,k)
          endif
        enddo
      enddo


c Do the vapour deposition calculation in the liquid part of mixed-phase clouds

c Calculate deposition on cloud ice

      do k=1,nl-1
        do mg=1,ln2
          if(cfrac(mg,k).gt.0.)then
            Tk=tliq(mg,k)+hlcp*(qlg(mg,k)+qfg(mg,k))/cfrac(mg,k) !T in liq cloud
            fl=qlg(mg,k)/(qfg(mg,k)+qlg(mg,k))
            if(Tk.lt.tfrz.and.Tk.ge.Tice.and.qlg(mg,k).gt.0.)then
              pk=100*prf(mg,k)
              es=100*exp(23.33086-6111.72784/Tk+0.15215*alog(Tk)) !ice value
              qs=0.622*es/pk !ice value
c              qs=qsati(pk,Tk)
c              es=qs*pk/0.622 !ice value
              Aprpr=hl/(rKa*Tk)*(hls/(rvap*Tk)-1)
              Bprpr=rvap*Tk/((Dva/pk)*es)
c              Cice=1.0e-2*exp(0.6*(tfrz-ttg(mg,k))) !Fletcher 1962
c              deles=esdiff (nint(Tk-tfrz))
              esl=100*exp(53.67957-6743.769/Tk-4.8451*alog(Tk))
              deles=esl-es
              Cice=1.e3*exp(12.96*deles/es - 0.639) !Meyers et al 1992
c              Cice=Cice/100.  !Sensitivity experiment
c              write(30,'(f12.3,g12.3)')Tk-tfrz,cice

              cm0=1.e-12 !Initial crystal mass
              qi0=cm0*Cice/rhoa(mg,k) !Initial ice mixing ratio

c Next 2 lines are for assumption of fully mixed ql and qf (also a line further down).

              qi0=max(qi0, qfg(mg,k)/cfrac(mg,k)) !Assume all qf and ql are mixed
              fd=1.   !Fraction of cloud in which deposition occurs
c              fd=fl   !Or, use option of adjacent ql,qf

              alf=1./3
              rhoic=700.
              Crate=7.8*((Cice/rhoa(mg,k))**2/rhoic)**(1./3) !Spheres
     &              *deles/((Aprpr+Bprpr)*es)
              qfdep=fd*cfrac(mg,k)*sqrt                         !Spheres
     &             (((2./3)*Crate*tdt+qi0**(2./3))**3)

c Also need this line for fully-mixed option...
              qfdep = qfdep - qfg(mg,k)

              qfdep=min(qfdep,qlg(mg,k))
              qlg(mg,k)=qlg(mg,k)-qfdep
              qfg(mg,k)=qfg(mg,k)+qfdep
            endif

            fice(mg,k)=qfg(mg,k)/(qfg(mg,k)+qlg(mg,k))

          endif
        enddo
      enddo  

c Recalculate qcg and cfrac for mixed-phase clouds that are fully glaciated
c Use qs=qsg.

      do k=1,nl
        do mg=1,ln2
          if(ttg(mg,k).lt.tfrz.and.ttg(mg,k).ge.Tice.and.
     &         fice(mg,k).gt.0.999999)then
            qs=qsg(mg,k)
            fi=1.
            hlrvap=(hl+fi*hlf)/rvap

            dqsdt=qs*hlrvap/tliq(mg,k)**2

            al=1/(1+(hlcp+fi*hlfcp)*dqsdt) !Smith's notation
            qc=qtot(mg,k)-qs

            delq=(1-rcrit(mg,k))*qs   !UKMO style (equivalent to above)
            cfrac(mg,k)=1.
            qcg(mg,k)=al*qc
            if(qc.lt.delq)then
              cfrac(mg,k)=1-0.5*((qc-delq)/delq)**2
              qcg(mg,k)=al*(qc-(qc-delq)**3/(6*delq**2))
            endif
            if(qc.le.0.)then
              cfrac(mg,k)=0.5*((qc+delq)/delq)**2
              qcg(mg,k)=al*(qc+delq)**3/(6*delq**2)
            endif
            if(qc.le.-delq)then
              cfrac(mg,k)=0.
              qcg(mg,k)=0.
            endif

            qfg(mg,k)=qcg(mg,k)
            qlg(mg,k)=0.
          endif

        enddo
      enddo


c Calculate new values of vapour mixing ratio and temperature

      do k=1,nl
        do mg=1,ln2
          qtg(mg,k)=qtot(mg,k)-qcg(mg,k)
          ttg(mg,k)=tliq(mg,k)+hlcp*qcg(mg,k)+hlfcp*qfg(mg,k)
          pk=100.0*prf(mg,k)
          qsg(mg,k)=qsati(pk,ttg(mg,k)) !Ice value
          ccov(mg,k)=cfrac(mg,k)      !Do this for now
        enddo
      enddo

c Vertically sub-grid cloud

      if(nl.lt.18)then
        do k=1,nl-1
          do mg=1,ln2
            if(cfrac(mg,k).gt.1.0e-2)then
              ccov(mg,k)=cfrac(mg,k)**(2./3)
            endif
          enddo
        enddo
      else
        do k=2,nl-1
          do mg=1,ln2
            if(cfrac(mg,k).gt.1.0e-2
     &           .and.cfrac(mg,k+1).eq.0..and.cfrac(mg,k-1).eq.0.)then
c            ccov(mg,k)=cfrac(mg,k)**(2./3)
              ccov(mg,k)=sqrt(cfrac(mg,k))
            endif
          enddo
        enddo
      endif


      
      if(debug)then
        if(lg.eq.lgdebug)then
          ns=insdebug
          mg=mgdebug+(ns-1)*lon
          write(25,'(a,3i3)')'IPASS=1, after newcloud'
          write(25,91)'ttg ',(ttg(mg,k),k=1,nl)
          write(25,91)'tliq ',(tliq(mg,k),k=1,nl)
          write(25,9)'qtg ',(qtg(mg,k),k=1,nl)
          write(25,9)'qsg ',(qsg(mg,k),k=1,nl)
          write(25,9)'cfrac ',(cfrac(mg,k),k=1,nl)
          write(25,9)'cdso4 ',(cdso4(mg,k),k=1,nl)
          write(25,9)'cfa ',(cfa(mg,k),k=1,nl)
          write(25,9)'qca ',(qca(mg,k),k=1,nl)
          write(25,9)'qtot ',(qtot(mg,k),k=1,nl)
          write(25,9)'qcg ',(qcg(mg,k),k=1,nl)
          write(25,9)'qlg ',(qlg(mg,k),k=1,nl)
          write(25,9)'qfg ',(qfg(mg,k),k=1,nl)
          write(25,*)
        endif

C***        do k=1,nl
C***          do mg=1,ln2
C***            if(qlg(mg,k).gt.1.e-4)then
C***              write(25,*)'mg,k,qlg(mg,k) ',mg,k,qlg(mg,k)
C***            endif
C***          enddo
C***        enddo
      endif
 1    format(3(a,g10.5))
 91   format(a,30f10.3)
 9    format(a,30g10.3)

      return
      end
