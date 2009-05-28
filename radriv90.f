      subroutine radrive (odcalc,iaero)
! Radiation driver routine for the conformal cubic model.
! This calls the GFDL radiation routines for each row of each face.
! At the moment it does not support the new liquid water cloud scheme
! It uses all the CSIRM GCM names for radiation fields internally, but copies
! to CC arrays at the end.
! albedo fixes for snow and ice are now done here

! This version can do multiple rows as once.
! In rdparm.h, set imax=2*il. For the conformal cube can also use imax=6*il
! N.B. (iq) indexing is still OK whether arrays have i dimension
!       of il or imax, because j advances sensibly

      use swr99_m
      use zenith_m
      use cc_mpi
      use diag_m
      use ateb ! MJT urban 
      use cable_ccam, only : CABLE ! MJT cable      
      include 'newmpar.h'
      parameter (ntest=0) ! N.B. usually j=1,7,13,19,...
!        for diag prints set ntest=1
!        or, usefully can edit 'ntest.gt.0' to 'ktau.gt.nnn'
      parameter (nalbwb=0)  ! 0  for original alb not depending on wb
      parameter (kcl_top=kl-2) !max level for cloud top (conjob,radrive,vertmix)
      include 'arrays.h'      
      include 'const_phys.h' ! for ldr cloud scheme
      include 'cparams.h'    ! for ldr cloud scheme
      include 'dates.h'      ! timer,kdate,ktime,dt,mtimer
      include 'extraout.h'   ! sintsave, etc
      include 'kuocom.h'     ! also with kbsav,ktsav
      include 'latlong.h'    ! rlatt,rlongg
      include 'liqwpar.h'    ! ifullw
      include 'nsibd.h'      ! rsmin,ivegt,sigmf,tgf,ssdn,res,rmc,tsigmf
      include 'parm.h'
      include 'pbl.h'
      include 'scamdim.h'
      include 'sigs.h'
      include 'soil.h'      ! land, rhgdum ... zmin  alb
      include 'soilsnow.h'  ! sicedep tgg,wb,snowd
      include 'soilv.h'
!     For the radiation code
      include 'rdparm.h'   ! imax
      include 'cldcom.h'
      include 'hcon.h'
      include 'lwout.h'
      include 'radisw.h'
      include 'raddiag.h'
      include 'rdflux.h'
      include 'srccom.h'
      include 'swocom.h'
      include 'tfcom.h'
      common/cfrac/cfrac(ifull,kl)
      common/work3c/rhg(ifull,kl) ! shared between cloud &radriv90
      common/work3d/rtt(ifull,kl) ! just to pass between radriv90 & globpe
      common/work3f/qccon(ifull,kl),qlrad(ifull,kl),qfrad(ifull,kl) ! ditto

      parameter(cong = cp/grav)
      parameter(csolar=1.96)
c     parameters for the aerosol calculation
      real beta_ave, alpha
      parameter(beta_ave = 0.29, alpha = 8.00)

!     input arguments
      logical odcalc  ! True for full radiation calculation

      real sigh(kl+1)

!     Radiation fields (CSIRO GCM names)
      real sg(imax), sgclr(imax), sint(imax), sout(imax), soutclr(imax)
      real rg(imax), rgclr(imax), rt(imax), rtclr(imax)
      real sga(imax),sgamp(ifull)
      real sgdn(imax), rgdn(imax)
      real hlwsav(ifull,kl),hswsav(ifull,kl)
      save hlwsav, hswsav, sgamp
!     when there is time incorporate properly here      
      common/radstuff/sgx(ifull),sgdnx(ifull),rgx(ifull),rgdnx(ifull),
     &             soutx(ifull),sintx(ifull),rtx(ifull)
      
c     Following are for cloud2 routine
      real t2(imax,kl),ql2(imax,kl),qf2(imax,kl),cf2(imax,kl),
     &     qc2(imax,kl),cd2(imax,kl),p2(imax,kl),
     &     dp2(imax,kl),cll(imax),clm(imax),clh(imax)
      logical land2(imax)

!     From initfs
c     Stuff from o3set
      common /o3dat/ dduo3n(37,kl),ddo3n2(37,kl),ddo3n3(37,kl),
     &               ddo3n4(37,kl)
c     Stuff from cldset
      common /clddat/ ccd(37,5),ccd2(37,5),ccd3(37,5),ccd4(37,5),
     &                kkth(37,5),kkbh(37,5)

!     For the zenith angle calculation
      real coszro2(imax), taudar2(imax)
      real fjd, r1, dlt, slag, dhr

!     Ozone returned by o3set
      real duo3n(imax,kl)

      logical clforflag, solarfit
      parameter (clforflag = .true., solarfit=.true.)
      logical cldoff
      logical first
      dimension ndoy(12)   ! days from beginning of year (1st Jan is 0)
      data ndoy/ 0,31,59,90,120,151,181,212,243,273,304,334/
      save first
      data first /.true./
      include 'establ.h'

      do k=kl,1,-1
         if(sig(k).le. .15)ksigtop=k  ! top level for RH calc for clouds
         sigh(k) = sigmh(k)
      end do
      sigh(kl+1) = 0.
      jdrad0=idjd/imax+1
      idrad=idjd-(jdrad0-1)*imax
      jdrad=1+(jdrad0-1)*imax/il  ! j increases in increments of imax/il

!     Initialisation (from initfs)
      if ( first ) then
         if(ntest.eq.1)print *,'id,jd,imax,idrad,jdrad0,jdrad ',
     .                          id,jd,imax,idrad,jdrad0,jdrad
         first = .false.
         call hconst
         call co2_read(sig)
         call radtable
         rrco2=rrvco2*ratco2mw
         if(amipo3)then
c           AMIP2 ozone
            call o3read_amip
            print *,'AMIP2 ozone input'
        else
c          Stuff from o3set
c          Rearrange the seasonal mean O3 data to allow interpolation
c          Define the amplitudes of the mean, annual and semi-annual cycles
           call o3_read(sig)
           call resetd(dduo3n,ddo3n2,ddo3n3,ddo3n4,37*kl)
        end if
      end if  ! (first)

C---------------------------------------------------------------------*
C START COMPUTATION                                                   *
C---------------------------------------------------------------------*

!     Set up number of minutes from beginning of year
!     This assumes 4-digit year already incorporated (fix done in infile)
!     For GCM runs assume year is <1980 (e.g. ~321-460 for 140 year run)
      jyear=kdate/10000
      jmonth=(kdate-jyear*10000)/100
      jday=kdate-jyear*10000-jmonth*100
      if(kdate.lt.19800000)jyear=1979 ! for gcm runs - but jyear not used
      jhour=ktime/100
      jmin=ktime-jhour*100
      mstart=1440*(ndoy(jmonth)+jday-1) + 60*jhour + jmin ! mins from start of y
!     timer contains number of hours since the start of the run.
!     mins = 60 * timer + mstart
!     mtimer contains number of minutes since the start of the run.
      mins = mtimer + mstart

c     Set number of years before present for orbital parameters.
c     Allowed values are 0, 6000 and 21000.
      bpyear = 0.
      if(nhstest<0)then  ! aquaplanet test
        fjd = 79.+mod(mins,1440)/1440.       ! set to 21 March +frac of day
      else
        fjd = float(mod(mins,525600))/1440.  ! 525600 = 1440*365
      endif
      if(ntest.gt.0)then
        print *,'kdate,jyear,jmonth,jhour,jmin,mtimer,mstart,mins,fjd ;'
     .          ,kdate,jyear,jmonth,jhour,jmin,mtimer,mstart,mins,fjd
      endif

      if(ldr==0)then
        do k=1,ksigtop  ! up to top level for RH calc for clouds
           do iq=1,ifull
              est = establ(t(iq,k))
              p = sig(k)*ps(iq)
!             rhg(iq,k) = qg(iq,k)/max(.622*est/(p-est),1.e-10) ! DARLAM
              rhg(iq,k) = qg(iq,k)/max(.622*est/(p-est),1.5e-6) ! C-C
           end do ! iq loop
        end do    ! k=1,ksigtop
        rhg(:,ksigtop+1:kl) = 0.
        if(ncvcloud.gt.0)then  ! jlm simple convective cloud enhancement
          frac=.01*ncvcloud    ! e.g. ncvcloud=90
          do k=kuocb,kcl_top
           do iq=1,ifull
             if(kbsav(iq).gt.0.and.k.ge.kbsav(iq).and.k.le.ktsav(iq))
     .         rhg(iq,k)=max(rhg(iq,k),frac) ! e.g. at least 90%
            enddo ! iq loop
          enddo   ! k loop
        endif     ! (ncvcloud.gt.0)
      endif       ! (ldr==0)

!     Calculate sun position
      if ( solarfit .or. odcalc ) then
         call solargh(fjd,bpyear,r1,dlt,alp,slag)
         ssolar = csolar / (r1**2)
         if(ntest.eq.9)ssolar=0.
      end if

!     Main loop over rows. imax/il is the number of rows done at once
      if(mod(ifull,imax).ne.0)then
        print *,'nproc,il,jl,ifull,imax ',nproc,il,jl,ifull,imax
        stop 'illegal setting of imax in rdparm'
      endif
      do 100 j=1,jl,imax/il
      istart=1+(j-1)*il
      iend=istart+imax-1
      if(ntest.eq.1)print *,'in radriv90 j = ',j
!     Calculate zenith angle for the solarfit calculation.
      if ( solarfit ) then
!        This call averages zenith angle just over this time step.
         dhr = dt/3600.0
c        call zenith(fjd,r1,dlt,slag,rlatt(1+(j-1)*il),
c    &               rlongg(1+(j-1)*il),dhr,imax,coszro2,taudar2)
         call zenith(fjd,r1,dlt,slag,rlatt(istart:iend),
     &               rlongg(istart:iend),dhr,imax,coszro2,taudar2)
         call atebccangle(istart,imax,coszro2(1:imax) ! MJT urban
     &    ,rlongg(istart:iend),rlatt(istart:iend),fjd,slag,dt
     &    ,sin(dlt)) 
      end if    !  ( solarfit )

      if ( odcalc ) then     ! Do the calculation

c     Average the zenith angle over the time (hours) between radiation
c     calculations
      dhr = kountr*dt/3600.0
      call zenith(fjd,r1,dlt,slag,rlatt(1+(j-1)*il),
     &            rlongg(1+(j-1)*il),dhr,imax,coszro,taudar)

      call atebccangle(istart,imax,coszro(1:imax) ! MJT urban
     &   ,rlongg(istart:iend),rlatt(istart:iend),fjd,slag,dt*kountr
     &   ,cos(dlt)) 

c     Set up basic variables, reversing the order of the vertical levels
      do i=1,imax
         iq=i+(j-1)*il
         temp(i,lp1) = tss(iq)
         press(i,lp1) = ps(iq) * 10. ! Convert to cgs
         cirab(i,1) = zero
      end do

c     Set up ozone for this time and row
      if (amipo3) then
         call o3set_amip ( rlatt(1+(j-1)*il:(j-1)*il+imax), imax, mins,
     &                     sigh, ps(1+(j-1)*il:(j-1)*il+imax), qo3 )
         qo3(:,:)=max(1.e-10,qo3(:,:))    ! July 2008
      else
         call o3set(rlatt(1+(j-1)*il),imax,mins,duo3n,sig)
c        Conversion of o3 from units of cm stp to gm/gm
         do k=1,kl
            do i=1,imax
c             qo3(i,k) = duo3n(i,k)*1.01325e+02/press(i,lp1)
              qo3(i,k) = max(1.e-10,duo3n(i,k)*1.01325e+02/press(i,lp1))
            end do
         end do
      end if

!     Set up surface albedo. The input value is > 1 over ocean points where
!     the zenith angle dependent formula should be used.
      do i=1,imax
          iq=i+(j-1)*il
          if( land(iq) )then
           if ((nsib.eq.CABLE).or.(nsib.eq.6)) then ! MJT CHANGE sib
             cuvrf(i,1) = albsav(iq) ! from cable (inc snow)
             cirrf(i,1) = albnirsav(iq) ! from cable (inc snow)
           else ! (nsib.le.3).or.(nsib.eq.5)
	   
           if(nalbwb.eq.0)then
             cuvrf(i,1) = albsav(iq)    ! use surface albedo from indata
             cirrf(i,1) = albnirsav(iq)    ! MJT CHANGE albedo
           else    ! soil albedo adjusted according to wetness of top layer
!            can do quadratic fit [ 0 to wbav to ssat]
!            wbs=sfc(isoilm(iq))         ! or consider using wbs=ssat()
             wbs=ssat(isoilm(iq))        ! or consider using wbs=ssat()
             wbw=swilt(isoilm(iq))       ! or consider using wbw=0.
             wbav=.5*(wbw+wbs)
             wbdp=min(max(wb(iq,5),wbw+.01) , wbs-.01)
!            assume albsav(iq) corresponds to value for wb=wbav
             albav=albsav(iq)/ (1.333*(wbdp-wbav)*(wbdp-wbs)
     .                           /((wbav-wbw)*(wbs-wbw))
     .          + (wbdp-wbw)*(wbdp-wbs)
     .                          /((wbav-wbw)*(wbav-wbs))
     .          + .666*(wbdp-wbw)*(wbdp-wbav)
     .                          /((wbs-wbw)*(wbs-wbav)))
             albs=(1.333*(wb(iq,1)-wbav)*(wb(iq,1)-wbs)
     .                           /((wbav-wbw)*(wbs-wbw))
     .          + (wb(iq,1)-wbw)*(wb(iq,1)-wbs)
     .                          /((wbav-wbw)*(wbav-wbs))
     .          + .666*(wb(iq,1)-wbw)*(wb(iq,1)-wbav)
     .                          /((wbs-wbw)*(wbs-wbav))) * albav

             albs=max(min(albs,.38),.11)
             cuvrf(i,1) = tsigmf(iq)*albsav(iq) + (1.-tsigmf(iq))*albs
             cirrf(i,1)=cuvrf(i,1) ! MJT CHANGE albedo
             if(ntest.gt.0.and.i.eq.idrad.and.j.eq.jdrad)then
               print *,'iq,land,sicedep,wbw,wbs,snowd ',
     .                  iq,land(iq),sicedep(iq),wbw,wbs,snowd(iq)
               print *,'wbdp,wpav,wb ',wbdp,wbav,wb(iq,1)
               print *,'albsav,albav,albs ',albsav(iq),albav,albs
               print *,'cuvrf1 ',cuvrf(i,1)
             endif
           endif    ! (nalbwb.eq.0)
           if(snowd(iq).gt.0.)then
c            new snow albedo (needs osnowd from the previous dt)
             dnsnow=min(1.,.1*max(0.,snowd(iq)-osnowd(iq)))  ! new snow (cm H2O)
c	     Snow age depends on snow crystal growth, freezing of melt water,
c	     accumulation of dirt and amount of new snow.
!	     if(isflag(iq).eq.1)then
!             ttbg=min(tggsn(iq,1),273.1)
!           else
!	       ttbg=min(tgg(iq,1),273.1)
!           endif
            ttbg=isflag(iq)*tggsn(iq,1) + (1-isflag(iq))*tgg(iq,1)
            ttbg=min(ttbg,273.1)
            ar1 = 5000.*( 1./273.1 - 1./ttbg) ! crystal growth  (-ve)
            exp_ar1=exp(ar1)                  ! e.g. exp(0 to -4)
            ar2 = 10.*ar1                     ! freezing of melt water
            exp_ar2=exp(ar2)                  ! e.g. exp(0 to -40)
            snr=snowd(iq)/max(ssdnn(iq),100.)
            if(isoilm(iq).eq.9)then   ! fixes for Arctic & Antarctic
              ar3=.001
              dnsnow=max(dnsnow,.0015)
              snrat=min(1.,snr/(snr+.001))
            else
              ar3=.3               ! accumulation of dirt
              snrat=min(1.,snr/(snr+.02))
            endif
            dtau=1.e-6*(exp_ar1+exp_ar2+ar3)*dt  ! <~.1 in a day
            if(snowd(iq).le. 1.)then
              snage(iq)=0.
            else
              snage(iq)=max(0.,(snage(iq) + dtau)*(1.-dnsnow))
            endif

c	     Snow albedo is dependent on zenith angle and  snow age.
            alvo = 0.95	        !alb. for vis. on a new snow
            aliro = 0.65        !alb. for near-infr. on a new snow
            fage = 1.-1./(1.+snage(iq))	 !age factor

            if(ntest.eq.1.and.iq.eq.idjd.and.mydiag)then
              print *,'ar1,ar2,snowd,ssdnn ',
     .                 ar1,ar2,snowd(iq),ssdnn(iq)
              print *,'exp_ar1,exp_ar2,ar3 ',
     .                 exp_ar1,exp_ar2,ar3
              print *,'dnsnow,snr,snrat,dtau,snage,fage ',
     .                 dnsnow,snr,snrat,dtau,snage(iq),fage
             endif

c	     albedo zenith dependence
c	     alvd = alvo * (1.0-cs*fage); alird = aliro * (1.-cn*fage)
c                   where cs = 0.2, cn = 0.5, b = 2.0
             cczen=max(.17365, coszro(i))
             fzen=( 1.+1./2.)/(1.+2.*2.*cczen) -1./2.
             if( cczen .gt. 0.5 ) fzen = 0.
             fzenm = max ( fzen, 0. )
             alvd = alvo * (1.0-0.2*fage)
             alv = .4 * fzenm * (1.-alvd) + alvd
             alird = aliro*(1.-.5*fage)
             alir = .4 * fzenm * (1.0-alird) + alird
             talb = .5 * ( alv + alir )        ! snow albedo
c	     cc=min(1.,snr/max(snr+2.*z0m(iq),0.02))
             cc=min(1.,snr/max(snr+zolnd(iq),0.02))
             tsigmfx=(1.-cc)*tsigmf(iq)      ! mult by snow free veg. fraction

            !alss = (1.-snrat)*albsav(iq) + snrat*talb ! canopy free surface albedo
            if(nsib.ge.3)then 
              !cuvrf(i,1)=alss
              cuvrf(i,1)=(1.-snrat)*cuvrf(i,1) + snrat*alv ! MJT CHANGE albedo
              cirrf(i,1)=(1.-snrat)*cirrf(i,1) + snrat*alir ! MJT CHANGE albedo
            else
              !cuvrf(i,1)=min(.8,(1.-tsigmfx)*alss+tsigmfx*albsav(iq))
              cuvrf(i,1)=min(.8,(1.-tsigmfx)*alv+tsigmfx*cuvrf(i,1)) ! MJT CHANGE albedo
              cirrf(i,1)=min(.8,(1.-tsigmfx)*alir+tsigmfx*cirrf(i,1)) ! MJT CHANGE albedo
            endif
c           old stuff before snage follows:
c           if(tss(iq) .ge. 273.09 ) then
c             snalb = 0.6   ! changed from .5 to .6 on 31/8/00
c           else
c             snalb = 0.8
c           endif
c           cuvrf(i,1)=min(snalb ,
c    .           albsav(iq)+(snalb-albsav(iq))*sqrt(snowd(iq)*.1))
            if(ntest.gt.0.and.i.eq.idrad.and.j.eq.jdrad)then
              print *,'i,j,land,sicedep,snowd,snrat ',
     .                 i,j,land(iq),sicedep(iq),snowd(iq),snrat
              print *,'albsav,dnsnow,talb,cuvrf1 ',
     .                 albsav(iq),dnsnow,talb,cuvrf(i,1)
            endif
           endif          !  snowd(iq).gt.0.
	   
           end if ! (nsib.eq.CABLE).or.(nsib.eq.6) ! MJT CHANGE sib
         else             !  over the ocean or sea ice
!          following for CCAM from 11/6/03
    !       cuvrf(i,1)=.65*fracice(iq)+
    ! .                (1.-fracice(iq))*.05/(coszro(i)+0.15)
           cuvrf(i,1)=.85*fracice(iq)+
     .                (1.-fracice(iq))*.05/(coszro(i)+0.15) ! MJT CHANGE albedo (follow CICE)
           cirrf(i,1)=.45*fracice(iq)+
     .                (1.-fracice(iq))*.05/(coszro(i)+0.15) ! MJT CHANGE albedo (follow CICE)
        endif       ! if( land(iq)) .. else..

      !-----------------------------------------------------------------------------------------------------------
      ! MJT albedo
      end do ! i=1,imax
      call atebalb1(istart,imax,cuvrf(1:imax,1),0) ! MJT CHANGE - urban
      call atebalb1(istart,imax,cirrf(1:imax,1),0) ! MJT CHANGE - urban
      if (iaero.ne.0) then
        do i=1,imax
          iq=i+(j-1)*il
           cosz = max ( coszro(i), 1.e-4)
           delta =  coszro(i)*beta_ave*alpha*so4t(iq)* ! still broadband
     &	            ((1.-0.5*(cuvrf(i,1)+cirrf(i,1)))/cosz)**2
           cuvrf(i,1)=min(0.99, delta+cuvrf(i,1)) ! surface albedo
           cirrf(i,1)=min(0.99, delta+cirrf(i,1)) ! still broadband
        end do ! i=1,imax
      endif !(iaero.ne.0)then
      albvisnir(istart:iend,1)=cuvrf(1:imax,1)
      albvisnir(istart:iend,2)=cirrf(1:imax,1)      
      !-----------------------------------------------------------------------------------------------------------      

      do k=1,kl
         kr = kl+1-k
         do i=1,imax
            iq=i+(j-1)*il
            press(i,kr) = ps(iq) * sig(k) * 10. ! Convert to cgs
            temp(i,kr) = t(iq,k)
c           Set min value to avoid numerical problems
            rh2o(i,kr) = max(qg(iq,k) ,1.e-7)
         end do ! i=1,imax
      end do    ! k=1,kl

!      print*, ' TEMP min', minval(temp), minloc(temp)
!      print*, ' TEMP max', maxval(temp), maxloc(temp)

c     Calculate half level pressures and temperatures by linear interp
      do k=1,kl+1
         kr=kl+2-k
         do i=1,imax
            iq=i+(j-1)*il
            press2(i,kr) = ps(iq) * sigh(k) * 10.
         end do ! i=1,imax
      end do    ! k=1,kl+1
      do k=2,kl
         kr=kl+2-k
         do i=1,imax
            temp2(i,kr) =
     &              ( temp(i,kr)*(sig(k)-sigh(k)) +
     &                temp(i,kr-1)*(sigh(k)-sig(k-1)) ) /
     &              (sig(k)-sig(k-1))
         end do ! i=1,imax
      end do    ! k=2,kl
      do i=1,imax
         temp2(i,1) = temp(i,1)
         temp2(i,lp1) = temp(i,lp1)
      end do ! i=1,imax
      
      if(ldr.ne.0)then  
c       Stuff needed for cloud2 routine...    
        qccon(:,:)=0.
        do k=1,kl   
          do i=1,imax
            iq=i+(j-1)*il
            t2(i,k)=t(iq,k)
            ql2(i,k)=qlrad(iq,k)
            qf2(i,k)=qfrad(iq,k)
            cf2(i,k)=cfrac(iq,k) ! called cfrad till Oct '05
            qc2(i,k)=qccon(iq,k)
!           will need this test eventually: if(naerosol_i(1).gt.0)then .. else
            if(land(iq))then
              if(rlatt(iq)>0.)then     
                cd2(i,k)=cdropl_nh
              else
                cd2(i,k)=cdropl_sh
              endif
            else
              if(rlatt(iq)>0.)then     
                cd2(i,k)=cdrops_nh
              else
                cd2(i,k)=cdrops_sh
              endif
            endif  ! (land(iq)) .. else ..
            land2(i)=land(iq)
            p2(i,k)=0.01*ps(iq)*sig(k) !Looks like ps is SI units
            dp2(i,k)=-0.01*ps(iq)*dsig(k) !dsig is -ve
          enddo
        enddo
      endif  ! (ldr.ne.0)

c  Clear sky calculation
      if (clforflag) then
        cldoff=.true.
c       set up cloud for this time and latitude
        if(ldr.ne.0)then  !Call LDR cloud scheme
c         write(24,*)coszro2
          call cloud2(cldoff,1,t2,ql2,qf2,cf2,qc2,
     &                cd2,land2,sigh,p2,dp2,coszro,      !Inputs
     &                cll,clm,clh)                       !Outputs
        else
          call cloud(cldoff,sig,j) ! jlm
        endif  ! (ldr.ne.0)
        if(ndi<0.and.nmaxpr==1)
     &     print *,'before swr99 ktau,j,myid ',ktau,j,myid
        call swr99(fsw,hsw,sg,ufsw,dfsw,press,press2,coszro,
     &             taudar,rh2o,rrco2,ssolar,qo3,nclds,
     &             ktopsw,kbtmsw,cirab,cirrf,cuvrf,camt,
     &             swrsave) ! MJT cable
        if(ndi<0.and.nmaxpr==1)
     &     print *,'after  swr99 ktau,j,myid ',ktau,j,myid
        do i=1,imax
          soutclr(i) = ufsw(i,1)*h1m3 ! solar out top
          sgclr(i)   = -fsw(i,lp1)*h1m3  ! solar absorbed at the surface
        end do
      else ! .not.clforflag
        do i=1,imax
          soutclr(i) = 0.
          sgclr(i) = 0.
        end do
      endif

c     Cloudy sky calculation
      cldoff=.false.
      if(ldr.ne.0)then  !Call LDR cloud scheme
c       write(24,*)
c       write(24,*)coszro2
        call cloud2(cldoff,1,t2,ql2,qf2,cf2,qccon,
     &              cd2,land2,sigh,p2,dp2,coszro,      !Inputs
     &              cll,clm,clh)                       !Outputs
        do i=1,imax
          iq=i+(j-1)*il
          cloudlo(iq)=cll(i)
          cloudmi(iq)=clm(i)
          cloudhi(iq)=clh(i)
        enddo
      else
        call cloud(cldoff,sig,j) ! jlm
      endif  ! (ldr.ne.0)
      call swr99(fsw,hsw,sg,ufsw,dfsw,press,press2,coszro,
     &           taudar,rh2o,rrco2,ssolar,qo3,nclds,
     &           ktopsw,kbtmsw,cirab,cirrf,cuvrf,camt,
     &           swrsave) ! MJT cable
      do i=1,imax
          sint(i) = dfsw(i,1)*h1m3   ! solar in top
          sout(i) = ufsw(i,1)*h1m3   ! solar out top
          sg(i)   = sg(i)*h1m3       ! solar absorbed at the surface
          iq=i+(j-1)*il              ! fixed Mar '05
          sgdn(i) = sg(i) / ( 1. - 0.5*sum(albvisnir(iq,:)) ) ! MJT albedo
      end do
      if(ntest.gt.0.and.j.eq.jdrad)then
        print *,'idrad,j,sint,sout,soutclr,sg,cuvrf1 ',
     .           idrad,j,sint(idrad),sout(idrad),soutclr(idrad),
     .           sg(idrad),cuvrf(idrad,1)
c       print *,'sint ',(sint(i),i=1,imax)
c       print *,'sout ',(sout(i),i=1,imax)
c       print *,'soutclr ',(soutclr(i),i=1,imax)
c       print *,'sg ',(sg(i),i=1,imax)
c       print *,'cuvrf ',(cuvrf(i),i=1,imax)
      endif
      call clo89
      if(ndi<0.and.nmaxpr==1)
     &     print *,'before lwr88 ktau,j,myid ',ktau,j,myid
      call lwr88

      do i=1,imax
         rt(i) = ( gxcts(i)+flx1e1(i) ) * h1m3          ! longwave at top
         rtclr(i) = ( gxctsclr(i)+flx1e1clr(i) ) * h1m3 ! clr sky lw at top
         rg(i) = grnflx(i)*h1m3                         ! longwave at surface
         rgclr(i) = grnflxclr(i)*h1m3         ! clear sky longwave at surface
         ! rg is net upwards = sigma T^4 - Rdown
         rgdn(i) = stefbo*temp(i,lp1)**4 - rg(i)
      end do

      do k=1,kl
         do i=1,imax
          iq=i+(j-1)*il
c         total heating rate
c----     note : htk now in Watts/M**2 (no pressure at level weighting)
c         Convert from cgs to SI units
          hswsav(iq,kl+1-k) = 0.001*hsw(i,k)
          hlwsav(iq,kl+1-k) = 0.001*heatra(i,k)
         end do
      end do

      if ( solarfit ) then
c       Calculate the amplitude of the diurnal cycle of solar radiation
c       at the surface (using the value for the middle of the radiation
c       step) and use this value to get solar radiation at other times.
c       Use the zenith angle and daylight fraction calculated in zenith
c       to remove these factors.

        do i=1,imax
           if ( coszro(i)*taudar(i) .le. 1.e-5 ) then ! 1.e-5 to avoid precision problems
c             The sun isn't up at all over the radiation period so no 
c             fitting need be done.
              sga(i) = 0.
           else
              sga(i) = sg(i) / (coszro(i)*taudar(i))
           end if
        end do
      else
        do i=1,imax
          sga(i) = 0.
        end do
      end if    !  ( solarfit )

!     Save things for non-radiation time steps
      fractss=.05
      do i=1,imax
         iq=i+(j-1)*il
         sgsave(iq) = sg(i)   ! repeated after solarfit
         sgamp(iq) = sga(i)
c        Save the value excluding Ts^4 part.  This is allowed to change.
         xxx = stefbo*tss(iq)**4
         rgsave(iq) = rg(i)-xxx  ! opposite sign to prev. darlam scam
!###     hlwsav(iq,1) = hlwsav(iq,1)-fractss*xxx  ! removed 18/6/03
         sintsave(iq) = sint(i) 
         rtsave(iq) = rt(i) 
         rtclsave(iq) = rtclr(i)  
         sgclsave(iq) = sgclr(i)
      end do

c     cloud amounts for saving
      do i=1,imax
         iq=i+(j-1)*il
         cloudtot(iq) = 1. - (1.-cloudlo(iq)) * (1.-cloudmi(iq)) *
     &        (1.-cloudhi(iq))
      end do

!     Use explicit indexing rather than array notation so that we can run
!     over the end of the first index
      if(ktau>1)then ! averages not added at time zero
        if(j==1)koundiag=koundiag+1  
        do i=1,imax
         iq=i+(j-1)*il
         sint_ave(iq) = sint_ave(iq) + sint(i)
         sot_ave(iq)  = sot_ave(iq)  + sout(i)
         soc_ave(iq)  = soc_ave(iq)  + soutclr(i)
         rtu_ave(iq)  = rtu_ave(iq)  + rt(i)
         rtc_ave(iq)  = rtc_ave(iq)  + rtclr(i)
         rgn_ave(iq)  = rgn_ave(iq)  + rg(i)
         rgc_ave(iq)  = rgc_ave(iq)  + rgclr(i)
         rgdn_ave(iq) = rgdn_ave(iq) + rgdn(i)
         sgdn_ave(iq) = sgdn_ave(iq) + sgdn(i)
         cld_ave(iq)  = cld_ave(iq)  + cloudtot(iq)
         cll_ave(iq)  = cll_ave(iq)  + cloudlo(iq)
         clm_ave(iq)  = clm_ave(iq)  + cloudmi(iq)
         clh_ave(iq)  = clh_ave(iq)  + cloudhi(iq)
        end do
      endif   ! (ktau>1)
      
      end if  ! odcalc
      
      if (solarfit) then 
!        Calculate the solar using the saved amplitude.
         do i=1,imax
          iq=i+(j-1)*il
          sg(i) = sgamp(iq)*coszro2(i)*taudar2(i)
         end do
      else
         do i=1,imax
          iq=i+(j-1)*il
          sg(i) = sgsave(iq)
         end do
      end if  ! (solarfit) .. else ..
      if(ktau>1)then ! averages not added at time zero
       do i=1,imax
         iq=i+(j-1)*il
         sgn_ave(iq)  = sgn_ave(iq)  + sg(i)
       end do
      endif  ! (ktau>1)
      
! Set up the CC model radiation fields
c slwa is negative net radiational htg at ground
! Note that this does not include the upward LW radiation from the surface.
! That is included in sflux.f
      do i=1,imax
         iq=i+(j-1)*il
         slwa(iq) = -sg(i)+rgsave(iq)
         sgsave(iq) = sg(i)   ! this is the repeat after solarfit 26/7/02
      end do
      if(odcalc.and.ndi<0.and.nmaxpr==1.and.idjd<=imax.and.mydiag)then
        print *,'bit after  lwr88 ktau,j,myid ',ktau,j,myid  
        print *,'sum_rg ',sum(rg(:))     
        print *,'slwa,sg,rgsave,rg,tss,grnflx ',slwa(idjd),sg(idjd),
     &           rgsave(idjd),rg(idjd),tss(idjd),grnflx(idjd)
      endif

! Calculate rtt, the net radiational cooling of atmosphere (K/s) from htk (in
! W/m^2 for the layer). Note that dsig is negative which does the conversion
! to a cooling rate.
      do k=1,kl
         do i=1,imax
            iq=i+(j-1)*il
            rtt(iq,k) = (hswsav(iq,k)+hlwsav(iq,k)) /
     &                   (cong*ps(iq)*dsig(k))
         end do
      end do
!     k = 1  ! these 6 lines removed 18/6/03
!     do i=1,imax
!      iq=i+(j-1)*il
!      rtt(iq,k)=(hswsav(iq,k)+hlwsav(iq,k)+fractss*stefbo*tss(iq)**4)/
!    &           (cong*ps(iq)*dsig(k))
!     end do

c     if ( j.eq.jdrad ) then
c        do k=kl,1,-1
c           print *,"hlwsav(i,j,k),dtlw(i,j,k)=",hlwsav(i,j,k),
c    &               dtlw(i,j,k)
c           print *,"hswsav(i,j,k),dtsw(i,j,k)=",hswsav(i,j,k),
c    &               dtsw(i,j,k)
c           print *,"rtt(i,j,k)=",rtt(i,j,k)
c        enddo
c       endif

!     clean up following, & use for rgsave etc
      sgx(istart:iend)=sg(:)
      sgdnx(istart:iend)=sgdn(:)
      rgx(istart:iend)=rg(:)
      rgdnx(istart:iend)=rgdn(:)
      soutx(istart:iend)=sout(:)
      sintx(istart:iend)=sint(:)
      rtx(istart:iend)=rt(:)

 100  continue  ! Row loop (j)  j=1,jl,imax/il
      if(ntest>0.and.mydiag)then
        print *,'rgsave,rtsave,sintsave ',
     .           rgsave(idjd),rtsave(idjd),sintsave(idjd)
        print *,'sgsave,rtclsave,sgclsave ',
     .           sgsave(idjd),rtclsave(idjd),sgclsave(idjd)
        print *,'alb ',albvisnir(idjd,1)
      endif
      if(nmaxpr==1.and.mydiag)then
        write (6,"('cfracr',9f8.3/6x,9f8.3)") cfrac(idjd,:)
        write (6,"('cloudlo,cloudmi,cloudhi,cloudtot',4f8.3)")
     .          cloudlo(idjd),cloudmi(idjd),cloudhi(idjd),cloudtot(idjd)
      endif
      return
      end
