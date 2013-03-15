c This routine is part of the prognostic cloud scheme. It calculates the
c formation and dissipation of cloud, and the liquid fraction in
c mixed-phase clouds. It is called by leoncld.
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
c      lgdebug - latitude index for single column debugging (1 for CCAM)
c      insdebug - hemisphere index for single column debugging
c      mgdebug  - index for single column debugging (set to idjd for CCAM)
c
c see also include files physparams.h (model physical constants)
c                        cparams.h    (cloud scheme parameters)
c
c from arguments
c      tdt - leapfrog timestep (seconds)
c      lg - latitude index (always 1 for CCAM)
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
c      qtg - water vapour mixing ratio (kg/kg) - called qenv in leoncld
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

      subroutine newcloud(tdt,lg,land,prf,kbase,ktop,rhoa,cdrop,  !Inputs
     &                    ttg,qtg,qlg,qfg,                        !In and out
     &                    cfrac,ccov,cfa,qca)                     !Outputs

c This routine is part of the prognostic cloud water scheme

      use diag_m
      use cc_mpi, only : mydiag
      use sigs_m
      implicit none
C Global parameters
      include 'newmpar.h'
      include 'const_phys.h' !Input physical constants
      include 'cparams.h'    !Input cloud scheme parameters
      include 'kuocom.h'     !Input cloud scheme parameters rcrit_l & rcrit_s
      include 'params.h'     !Input model grid dimensions (modified params.h for CCAM)
      include 'parm.h'

C Argument list
      real tdt,den1
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
      real cfa(ln2,nl)
      real qca(ln2,nl)
      real qsl(ln2,nl)
      real qsw(ln2,nl)


C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
c      include 'CNSTA.f'
c      include 'FEWFLAGS.f'   !Input debug, lgdebug etc.
c     include 'TIMEX.f'

C Local work arrays and variables
      real qcg(ln2,nl)
      real qtot(ln2,nl),tliq(ln2,nl)
      real fice(ln2,nl)
      real qcold(ln2,nl)
      real cdrop(ln2,nl)
      real rcrit(ln2,nl)
      real dum(nl)

      integer k
      integer mg

      real al
      real alf
      real aprpr
      real bprpr
      real cice
      real cm0
      real crate
      real decayfac
      real deles
      real delq
      real dqsdt
      real es
      real fd
      real fl
      real hlrvap
      real pk
      real qc
      real qfdep
      real qfnew
      real qi0
      real qs, qsi
      real rhoic
      real root6
      real root6i
      real tk

C Local data, functions etc
      logical ukconv
      integer naerosol_i(2)
      data  ukconv, naerosol_i / .false., 2*0 /
      logical debug
      integer lgdebug,insdebug   ! 1,idjd
      data debug,lgdebug,insdebug /.false.,1,1/
      save ukconv,naerosol_i,debug,lgdebug,insdebug
      real tdiffx,tx_,esdiffx

      real esdiff(-40:2)  !Table of es(liq) - es(ice) (MKS), -40 <= t <= 0 C
      data esdiff / 
     & 6.22, 6.76, 7.32, 7.92, 8.56, 9.23, 9.94,10.68,11.46,12.27,
     & 13.11,13.99,14.89,15.82,16.76,17.73,18.70,19.68,20.65,21.61,
     & 22.55,23.45,24.30,25.08,25.78,26.38,26.86,27.18,27.33,27.27,
     & 26.96,26.38,25.47,24.20,22.51,20.34,17.64,14.34,10.37, 5.65,
     & 0.08, 0., 0. /
      tdiffx(tx_)=min(max( tx_-tfrz, -40.), 1.)
      esdiffx(tx_) = 
     &    (1.-(tdiffx(tx_)-aint(tdiffx(tx_))))*esdiff(int(tdiffx(tx_)))
     &  + (tdiffx(tx_)-aint(tdiffx(tx_)))*esdiff(int(tdiffx(tx_))+1)
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
      qsati(pp,t) = epsil*estabi(t)/max(.1,pp-estabi(t)) 

C Start code : ----------------------------------------------------------

      root6=sqrt(6.)
      root6i=1./root6

      if(diag.and.mydiag)then
          write(6,*) 'entering newcloud'
          write(6,91) 'prf ',(prf(idjd,k),k=1,nl)
          write(6,91) 'ttg ',(ttg(idjd,k),k=1,nl)
          write(6,*) 'qtg ',(qtg(idjd,k),k=1,nl)
          write(6,*) 'qlg ',(qlg(idjd,k),k=1,nl)
          write(6,*) 'qfg ',(qfg(idjd,k),k=1,nl)
      endif

c Define cdrop  - passed through as cdso4, defined in leoncld.f

c First melt cloud ice or freeze cloud water to give correct ice fraction fice.
c Then calculate the cloud conserved variables qtot and tliq.
c Note that qcg is the total cloud water (liquid+frozen)

      do k=1,nl
        do mg=1,ln2
          if(ttg(mg,k)>=tfrz)then
            fice(mg,k)=0.
          elseif(ttg(mg,k)>=tice)then
            if(qfg(mg,k)>1.0e-12)then
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

      if(diag.and.mydiag)then
          write(6,*) 'within newcloud'
          write(6,*) 'ttg ',ttg(idjd,:)
          write(6,*) 'qcold ',qcold(idjd,:)
          write(6,*) 'qcg ',qcg(idjd,:)
          write(6,*) 'qlg ',qlg(idjd,:)
          write(6,*) 'qfg ',qfg(idjd,:)
          write(6,*) 'fice ',fice(idjd,:)
      endif
      
      if (ncloud<3) then ! usual diagnostic cloud
      
c Precompute the array of critical relative humidities 

        if(nclddia==-3)then
          do k=1,nl
           do mg=1,ln2
            if(land(mg))then
              rcrit(mg,k)=max( rcrit_l , (1.-16.*(1.-sig(k))**3) ) 
            else
              rcrit(mg,k)=max( rcrit_s , (1.-16.*(1.-sig(k))**3) ) 
            endif 
           enddo
          enddo
        elseif(nclddia<0)then
          do k=1,nl
           do mg=1,ln2
            if(land(mg))then
              rcrit(mg,k)=max( rcrit_l , (1.-4.*(1.-sig(k))**2) ) 
            else
              rcrit(mg,k)=max( rcrit_s , (1.-4.*(1.-sig(k))**2) ) 
            endif 
           enddo
          enddo
         elseif(nclddia==1)then
          do k=1,nl
           do mg=1,ln2
            if(land(mg))then
              rcrit(mg,k)=max(rcrit_l,sig(k)**3)          ! .75 for R21 Mk2
            else
              rcrit(mg,k)=max(rcrit_s,sig(k)**3)          ! .85 for R21 Mk2
            endif
           enddo
          enddo
        elseif(nclddia==2)then
          do k=1,nl
           do mg=1,ln2
            if(land(mg))then
              rcrit(mg,k)=rcrit_l
            else
              rcrit(mg,k)=rcrit_s
            endif
           enddo
          enddo
         elseif(nclddia==3)then
          do k=1,nl
           do mg=1,ln2
            if(land(mg))then
              rcrit(mg,k)=max(rcrit_l,sig(k)**2)          ! .75 for R21 Mk2
            else
              rcrit(mg,k)=max(rcrit_s,sig(k)**2)          ! .85 for R21 Mk2
            endif
           enddo
          enddo
        elseif(nclddia==4)then
          do k=1,nl
           do mg=1,ln2
            if(land(mg))then
              rcrit(mg,k)=max(rcrit_l,sig(k))          ! .75 for test Mk2/3
            else
              rcrit(mg,k)=max(rcrit_s,sig(k))          ! .9  for test Mk2/3
            endif
           enddo
          enddo
        elseif(nclddia==5)then  ! default till May 08
          do k=1,nl
           do mg=1,ln2
            if(land(mg))then
              rcrit(mg,k)=max(rcrit_l,min(.99,sig(k))) ! .75 for same as T63
            else
              rcrit(mg,k)=max(rcrit_s,min(.99,sig(k))) ! .85 for same as T63
            endif
           enddo
          enddo
         elseif(nclddia==6)then
          do k=1,nl
           do mg=1,ln2
            if(land(mg))then
              rcrit(mg,k)=max(rcrit_l*(1.-.15*sig(k)),sig(k)**4)         
            else
              rcrit(mg,k)=max(rcrit_s*(1.-.15*sig(k)),sig(k)**4)          
            endif
           enddo
          enddo
         elseif(nclddia==7)then
          do k=1,nl
           do mg=1,ln2
            if(land(mg))then
              rcrit(mg,k)=max(rcrit_l*(1.-.2*sig(k)),sig(k)**4)         
            else
              rcrit(mg,k)=max(rcrit_s*(1.-.2*sig(k)),sig(k)**4)          
            endif
           enddo
          enddo
        endif  ! (nclddia<0)  .. else ..

c Calculate cloudy fraction of grid box (cfrac) and gridbox-mean cloud water
c using the triangular PDF of Smith (1990)

        do k=1,nl
          do mg=1,ln2
            hlrvap=(hl+fice(mg,k)*hlf)/rvap
            qtot(mg,k)=qtg(mg,k)+qcg(mg,k)
            tliq(mg,k)=ttg(mg,k)-hlcp*qcg(mg,k)-hlfcp*qfg(mg,k)

c Calculate qs and gam=(L/cp)*dqsdt,  at temperature tliq
            pk=100.0*prf(mg,k)
            qsi=qsati(pk,tliq(mg,k))           !Ice value
            !deles=esdiff(min(max(-40,(nint(tliq(mg,k)-tfrz))),1)) ! MJT suggestion
            deles=esdiffx(tliq(mg,k))                              ! MJT suggestion
            qsl(mg,k)=qsi+epsil*deles/pk       !qs over liquid
            qsw(mg,k)=fice(mg,k)*qsi+(1.-fice(mg,k))*qsl(mg,k) !Weighted qs at temperature Tliq
            qs=qsw(mg,k)
            dqsdt=qs*hlrvap/tliq(mg,k)**2
c           qvc(mg,k)=qs !Vapour mixing ratio in cloud

            al=1./(1.+(hlcp+fice(mg,k)*hlfcp)*dqsdt)    !Smith's notation
            qc=qtot(mg,k)-qs

            delq=(1.-rcrit(mg,k))*qs      !UKMO style (equivalent to above)
            cfrac(mg,k)=1.
            qcg(mg,k)=al*qc
            if(qc<delq)then
              cfrac(mg,k)=max(1.e-6 , 1.-.5*((qc-delq)/delq)**2)     ! for roundoff
              qcg(mg,k)=max(1.e-8,al*(qc-(qc-delq)**3/(6.*delq**2))) ! for roundoff
            endif
            if(qc<=0.)then
              cfrac(mg,k)=max(1.e-6 , .5*((qc+delq)/delq)**2)    ! for roundoff
              qcg(mg,k)=max(1.e-8, al*(qc+delq)**3/(6.*delq**2)) ! for roundoff
            endif
            if(qc<=-delq)then
              cfrac(mg,k)=0.
              qcg(mg,k)=0.
            endif
c            ! Roundoff check
c            if ( qcg(mg,k) <= 0. ) then
c               cfrac(mg,k) = 0.
c               qcg  (mg,k) = 0.  ! added Oct '06
c            end if

c Calculate the cloud fraction (cfa) in which ql exceeds qcrit, and
c the corresponding gridbox-mean cloud water mixing ratio qca. 
c This (qca) is the cloud-water mixing ratio inside cfa divided by cfa.
c The new variable qc2 is like qc above, but is used for integration limits
c only, not the integrand

C***          qcrit=(4*pi/3)*rhow*Rcm**3*cdrop(mg,k)/rhoa(mg,k)
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
            ! Set these to zero so negative indefinite initialisation works
            cfa(mg,k)=0.
            qca(mg,k)=0.

          enddo
        enddo

        if(diag.and.mydiag)then
            write(6,*) 'rcrit ',rcrit(idjd,:)
            write(6,*) 'qtot ',qtot(idjd,:)
            do k=1,nl
             dum(k)=qsati(100.*prf(idjd,k),tliq(idjd,k)) 
            enddo
            write(6,*) 'qsi',dum(:)
            write(6,*) 'nint',nint(tliq(idjd,:)-tfrz)
            write(6,*) 'deles',
     &                 esdiff(min(max(-40,(nint(tliq(idjd,:)-tfrz))),1))
            write(6,*) 'tliq',tliq(idjd,:)
            write(6,*) 'qsl ',qsl(idjd,:)
            write(6,*) 'qsw ',qsw(idjd,:)
            write(6,*) 'cfrac ',cfrac(idjd,:)
            write(6,*) 'qc  ',  qtot(idjd,:)-qsw(idjd,:)
            write(6,*) 'qcg ',qcg(idjd,:)
            write(6,*) 'delq ', (1.-rcrit(idjd,:))*qsw(idjd,:)
        endif
        
      else ! prognostic cloud
      
        ! Tiedtke prognostic cloud model
        write(6,*) "Prognositc cloud not implemented"
        stop
      
        !call progcloud(cfrac,qcg)
        cfa=0.
        qca=0.

      end if ! ncloud<3 ..else..

c Assume condensation or evaporation retains ice fraction fice.
c Introduce a time-decay factor for cirrus (as suggested by results of Khvorostyanov & Sassen,
c JAS, 55, 1822-1845, 1998). Their suggested range for the time constant is 0.5 to 2 hours.
c The grid-box-mean values of qtg and ttg are adjusted later on (below).

      do k=1,nl
        do mg=1,ln2
          if(ttg(mg,k)>=Tice)then

            qfg(mg,k) = fice(mg,k)*qcg(mg,k)
            qlg(mg,k) = qcg(mg,k) - qfg(mg,k)

          else                                   ! Cirrus T range
            decayfac = exp ( (-1./7200) * tdt )  ! Try 2 hrs
c            decayfac = 0.                         ! Instant adjustment (old scheme)
            qfg(mg,k) = qcold(mg,k)*decayfac + qcg(mg,k)*(1.-decayfac)
            qcg(mg,k) = qfg(mg,k)
          endif
        enddo
      enddo

c Do the vapour deposition calculation in mixed-phase clouds:
c Calculate deposition on cloud ice, assuming es(T) is the weighted value of the 
c liquid and ice values.

      do k=1,nl   ! was nl-1 until Sept '04
        do mg=1,ln2
          if(cfrac(mg,k)>0.)then
            Tk=tliq(mg,k)+hlcp*(qlg(mg,k)+qfg(mg,k))/cfrac(mg,k) !T in liq cloud
            fl=qlg(mg,k)/(qfg(mg,k)+qlg(mg,k))
            if(Tk<tfrz.and.qlg(mg,k)>0.)then
              pk=100*prf(mg,k)
              qs=qsati(pk,Tk)
              es=qs*pk/0.622 !ice value
              Aprpr=hl/(rKa*Tk)*(hls/(rvap*Tk)-1.)
              Bprpr=rvap*Tk/((Dva/pk)*es)
!              deles=(1.-fice(mg,k))                           ! MJT suggestion
!     &              *esdiff(min(max(-40,(nint(Tk-tfrz))),1))  ! MJT suggestion
              deles=(1.-fice(mg,k))*esdiffx(Tk)                ! MJT suggestion
              Cice=1.e3*exp(12.96*deles/es - 0.639) !Meyers et al 1992

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

c Calculate new values of vapour mixing ratio and temperature

      do k=1,nl
        do mg=1,ln2
          qtg(mg,k)=qtot(mg,k)-qcg(mg,k)
          ttg(mg,k)=tliq(mg,k)+hlcp*qcg(mg,k)+hlfcp*qfg(mg,k)
c         pk=100.0*prf(mg,k)                         ! moved to bottom
c         qsg(mg,k)=qsati(pk,ttg(mg,k)) !Ice value   ! moved to bottom
          ccov(mg,k)=cfrac(mg,k)      !Do this for now
        enddo
      enddo

c Vertically sub-grid cloud

!      if(nl.lt.18)then
!        do k=1,nl-1
!          do mg=1,ln2
!            if(cfrac(mg,k).gt.1.0e-2)then
!              ccov(mg,k)=cfrac(mg,k)**(2./3)
!            endif
!          enddo
!        enddo
!      else
        do k=2,nl-1
          do mg=1,ln2
            if(cfrac(mg,k)>1.0e-2
     &           .and.cfrac(mg,k+1)==0..and.cfrac(mg,k-1)==0.)then
c            ccov(mg,k)=cfrac(mg,k)**(2./3) ! Mk3.6
              ccov(mg,k)=sqrt(cfrac(mg,k))
            endif
          enddo
        enddo
!      endif
     
      if(diag.and.mydiag)then
          write(6,*) 'at end of newcloud'
          write(6,*) 'ttg ',ttg(idjd,:)
          write(6,*) 'qcg ',qcg(idjd,:)
          write(6,*) 'qlg ',qlg(idjd,:)
          write(6,*) 'qfg ',qfg(idjd,:)
          write(6,*) 'qtg ',qtg(idjd,:)
      endif
 91   format(a6,9f10.3/(6x,9f10.3))
 9    format(a6,9g10.3/(6x,9g10.3))

      return
      end
