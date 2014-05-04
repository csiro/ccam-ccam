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
c      
c      land - logical variable for surface type ( = T for land points)
c      prf - pressure at full levels (in hPa. NB: not SI units)
c      kbase - k index of convective cloud base   *** not used
c      ktop - k index of convective cloud top     *** not used
c
c In/Out:
c
c from argumen
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

      subroutine newcloud(tdt,land,prf,kbase,ktop,rhoa,cdrop,     !Inputs
     &                    ttg,qtg,qlg,qfg,                        !In and out
     &                    cfrac,ccov,cfa,qca)                     !Outputs

c This routine is part of the prognostic cloud water scheme

      use cloudmod
      use estab, only : esdiffx, qsati
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
      logical land(ln2)
      real prf(ln2,nl)
      integer kbase(ln2)
      integer ktop(ln2)
      real rhoa(ln2,nl)
      real ttg(ln2,nl)
      real qtg(ln2,nl)
      real qlg(ln2+iextra,nl)
      real qfg(ln2+iextra,nl)
      real cfrac(ln2,nl)
      real ccov(ln2,nl)
      real cfa(ln2,nl)
      real qca(ln2,nl)
      real qsl(ln2,nl)
      real qsw(ln2,nl)


C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Local work arrays and variables
      real qcg(ln2,nl)
      real qtot(ln2,nl),tliq(ln2,nl)
      real fice(ln2,nl)
      real qcold(ln2,nl)
      real cdrop(ln2,nl)
      real rcrit(ln2,nl)
      real qsi(ln2,nl)

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
      real qfnew(ln2,nl)
      real qi0
      real qs
      real rhoic
      real root6
      real root6i
      real tk
      real qcic,qcrit,qc2,qto,wliq,r3c,r6c,eps,beta6
      real fi

C Local data, functions etc
      logical ukconv
      integer naerosol_i(2)
      data  ukconv, naerosol_i / .false., 2*0 /
      logical debug
      integer lgdebug,insdebug   ! 1,idjd
      data debug,lgdebug,insdebug /.false.,1,1/
      save ukconv,naerosol_i,debug,lgdebug,insdebug

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
        end do
      end do
      qcg(:,:)=qlg(1:ifull,:)+qfg(1:ifull,:)
      qcold(:,:)=qcg(:,:)
      qfnew=fice(:,:)*qcg(:,:)
      ttg(:,:)=ttg(:,:)+hlfcp*(qfnew-qfg(1:ifull,:)) !Release L.H. of fusion
      qfg(1:ifull,:)=fice(:,:)*qcg(:,:)
      qlg(1:ifull,:)=max(0.,qcg(:,:)-qfg(1:ifull,:))

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
      
! Precompute the array of critical relative humidities 

        if(nclddia==-3)then
          do k=1,nl
           where(land(1:ifull))
            rcrit(:,k)=max( rcrit_l, (1.-16.*(1.-sig(k))**3) )
           elsewhere
            rcrit(:,k)=max( rcrit_s, (1.-16.*(1.-sig(k))**3) )
           end where
          enddo
        elseif(nclddia<0)then
          do k=1,nl
           where(land(1:ifull))
            rcrit(:,k)=max( rcrit_l, (1.-4.*(1.-sig(k))**2) )
           elsewhere
            rcrit(:,k)=max( rcrit_s, (1.-4.*(1.-sig(k))**2) )
           end where
          enddo
         elseif(nclddia==1)then
          do k=1,nl
           where(land(1:ifull))
            rcrit(:,k)=max( rcrit_l, sig(k)**3 )
           elsewhere
            rcrit(:,k)=max( rcrit_s, sig(k)**3 )
           end where
          enddo
        elseif(nclddia==2)then
          do k=1,nl
           where(land(1:ifull))
            rcrit(:,k)=rcrit_l
           elsewhere
            rcrit(:,k)=rcrit_s
           end where
          enddo
         elseif(nclddia==3)then
          do k=1,nl
           where(land(1:ifull))
            rcrit(:,k)=max( rcrit_l, sig(k)**2 )          ! .75 for R21 Mk2
           elsewhere
            rcrit(:,k)=max( rcrit_s, sig(k)**2 )          ! .85 for R21 Mk2
           end where
          enddo
        elseif(nclddia==4)then
          do k=1,nl
           where(land(1:ifull))
            rcrit(:,k)=max( rcrit_l, sig(k) )             ! .75 for test Mk2/3
           elsewhere
            rcrit(:,k)=max( rcrit_s, sig(k) )             ! .9  for test Mk2/3
           end where
          enddo
        elseif(nclddia==5)then  ! default till May 08
          do k=1,nl
           where(land(1:ifull))
            rcrit(:,k)=max( rcrit_l, min(.99,sig(k)) )    ! .75 for same as T63
           elsewhere
            rcrit(:,k)=max( rcrit_s, min(.99,sig(k)) )    ! .85 for same as T63
           end where
          enddo
         elseif(nclddia==6)then
          do k=1,nl
           where(land(1:ifull))
            rcrit(:,k)=max(rcrit_l*(1.-.15*sig(k)),sig(k)**4)
           elsewhere
            rcrit(:,k)=max(rcrit_s*(1.-.15*sig(k)),sig(k)**4)
           end where
          enddo
         elseif(nclddia==7)then
          do k=1,nl
           where(land(1:ifull))
            rcrit(:,k)=max(rcrit_l*(1.-.2*sig(k)),sig(k)**4)
           elsewhere
            rcrit(:,k)=max(rcrit_s*(1.-.2*sig(k)),sig(k)**4)
           end where
          enddo
        endif  ! (nclddia<0)  .. else ..

c Calculate cloudy fraction of grid box (cfrac) and gridbox-mean cloud water
c using the triangular PDF of Smith (1990)

        do k=1,nl
          do mg=1,ln2
            fi=fice(mg,k)
            hlrvap=(hl+fi*hlf)/rvap
            qtot(mg,k)=qtg(mg,k)+qcg(mg,k)
            tliq(mg,k)=ttg(mg,k)-hlcp*qcg(mg,k)-hlfcp*qfg(mg,k)

c Calculate qs and gam=(L/cp)*dqsdt,  at temperature tliq
            pk=100.0*prf(mg,k)
            qsi(mg,k)=qsati(pk,tliq(mg,k))           !Ice value
            !deles=esdiff(min(max(-40,(nint(tliq(mg,k)-tfrz))),1)) ! MJT suggestion
            deles=esdiffx(tliq(mg,k))                              ! MJT suggestion
            qsl(mg,k)=qsi(mg,k)+epsil*deles/pk       !qs over liquid
            qsw(mg,k)=fi*qsi(mg,k)+(1.-fi)*qsl(mg,k) !Weighted qs at temperature Tliq
            qs=qsw(mg,k)
            dqsdt=qs*hlrvap/tliq(mg,k)**2
c           qvc(mg,k)=qs !Vapour mixing ratio in cloud

            al=1./(1.+(hlcp+fi*hlfcp)*dqsdt)    !Smith's notation
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

c Calculate the cloud fraction (cfa) in which ql exceeds qcrit, and
c the corresponding gridbox-mean cloud water mixing ratio qca. 
c This (qca) is the cloud-water mixing ratio inside cfa times cfa.
c The new variable qc2 is like qc above, but is used for integration limits
c only, not the integrand

            if(cfrac(mg,k)>0.)then
              qcic=qcg(mg,k)/cfrac(mg,k) !Mean in cloud value

c Following few lines are for Yangang Liu's new scheme (2004: JAS, GRL)
c Need to do first-order estimate of qcrit using mean in-cloud qc (qcic)

              Wliq = max( 1.e-10, 1000. * qcic * rhoa(mg,k)) !g/m3
              R6c = 4.09e-4
     &           * ( 1.15e23*1.e-6*cdrop(mg,k) / Wliq**2 ) ** (1./6.)
              eps = 1. - 0.7 * exp(-0.003e-6*cdrop(mg,k)) !mid range
              beta6 = ((1.+3.*eps**2)*(1.+4.*eps**2)*(1.+5.*eps**2)
     &             / ((1.+eps**2)*(1.+2.*eps**2)) )**(1./6.)
              R3c = 1.e-6*R6c/beta6 !in metres
              qcrit=(4.*pi/3.)*rhow*R3c**3*cdrop(mg,k)/rhoa(mg,k) !New qcrit

              qc2=qtot(mg,k)-qs-qcrit/al
              cfa(mg,k)=1.
              qca(mg,k)=al*qc
              if(qc2<delq)then
                cfa(mg,k)=1.-0.5*((qc2-delq)/delq)**2
                qto=(qtot(mg,k)-delq+2.*(qs+qcrit/al))/3.
                qca(mg,k)=al*(qtot(mg,k) - qto + cfa(mg,k)*(qto-qs))
              endif
              if(qc2<=0.)then
                cfa(mg,k)=0.5*((qc2+delq)/delq)**2
                qca(mg,k)=cfa(mg,k)*(al/3.)*(2.*qcrit/al + qc+delq)
              endif
              if(qc2<=-delq)then
                cfa(mg,k)=0.
                qca(mg,k)=0.
              endif
            else
              cfa(mg,k)=0.
              qca(mg,k)=0.
            endif

          enddo
        enddo

        if(diag.and.mydiag)then
            write(6,*) 'rcrit ',rcrit(idjd,:)
            write(6,*) 'qtot ',qtot(idjd,:)
            write(6,*) 'qsi',qsi(mg,k)
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
        qtot(:,:)=qtg(:,:)+qcg(:,:)
        tliq(:,:)=ttg(:,:)-hlcp*qcg(:,:)-hlfcp*qfg(1:ifull,:)
        call progcloud(cfrac,qcg,qtg,ttg,prf,rhoa,fice)
        
        ! Use 'old' autoconversion with prognostic cloud
        cfa(:,:)=0.
        qca(:,:)=0.

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

          else                                    ! Cirrus T range
            decayfac = exp ( (-1./7200.) * tdt )  ! Try 2 hrs
c           decayfac = 0.                         ! Instant adjustment (old scheme)
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
            fl=qlg(mg,k)/max(qfg(mg,k)+qlg(mg,k),1.e-30)
            if(Tk<tfrz.and.qlg(mg,k)>1.e-8)then
              pk=100*prf(mg,k)
              qs=qsati(pk,Tk)
              es=qs*pk/0.622 !ice value
              Aprpr=hl/(rKa*Tk)*(hls/(rvap*Tk)-1.)
              Bprpr=rvap*Tk/((Dva/pk)*es)
              deles=(1.-fice(mg,k))*esdiffx(Tk)
              Cice=1.e3*exp(12.96*deles/es - 0.639) !Meyers et al 1992

              cm0=1.e-12 !Initial crystal mass
              qi0=cm0*Cice/rhoa(mg,k) !Initial ice mixing ratio

c Next 2 lines are for assumption of fully mixed ql and qf (also a line further down).

              qi0=max(qi0, qfg(mg,k)/cfrac(mg,k)) !Assume all qf and ql are mixed
              fd=1.   !Fraction of cloud in which deposition occurs
c              fd=fl   !Or, use option of adjacent ql,qf

              alf=1./3.
              rhoic=700.
              Crate=7.8*((Cice/rhoa(mg,k))**2/rhoic)**(1./3.) !Spheres
     &              *deles/((Aprpr+Bprpr)*es)
              qfdep=fd*cfrac(mg,k)*sqrt                       !Spheres
     &             (((2./3.)*Crate*tdt+qi0**(2./3.))**3)

c Also need this line for fully-mixed option...
              qfdep = qfdep - qfg(mg,k)

              qfdep=min(qfdep,qlg(mg,k))
              qlg(mg,k)=qlg(mg,k)-qfdep
              qfg(mg,k)=qfg(mg,k)+qfdep
            endif

            fice(mg,k)=qfg(mg,k)/max(qfg(mg,k)+qlg(mg,k),1.e-30)

          endif
        enddo
      enddo  

c Calculate new values of vapour mixing ratio and temperature

      qtg(:,:)=qtot(:,:)-qcg(:,:)
      ttg(:,:)=tliq(:,:)+hlcp*qcg(:,:)+hlfcp*qfg(1:ifull,:)
      ccov(:,:)=cfrac(:,:)      !Do this for now

c Vertically sub-grid cloud

      do k=2,nl-1
        do mg=1,ln2
          if(cfrac(mg,k)>1.0e-2
     &         .and.cfrac(mg,k+1)==0..and.cfrac(mg,k-1)==0.)then
c           ccov(mg,k)=cfrac(mg,k)**(2./3) ! Mk3.6
            ccov(mg,k)=sqrt(cfrac(mg,k))
          endif
        enddo
      enddo
     
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
