c This is the interface between the Fels-Schwarzkopf radiation scheme and
c LDR's prognostic cloud water scheme. It is called by radrive.
c
c INPUT/OUTPUT:
c Input:
c
c      dprf - pressure thickness at each sigma level
c      prf  - pressure at full levels
c      coszro - zenith angle at grid pt. (used in SW scheme)
c      cldoff - flag set to T for clear sky radiation calculations
c      lg     - latitude index
c      ttg    - temperature
c      qlg    - cloud liquid water mixing ratio (kg/kg)
c      qfg    - cloud ice mixing ratio (kg/kg)
c      cfrac  - total cloud fraction (stratiform + convective)
c      qccon  - cloud water mixing ratio of convective clouds (kg/kg)
c
c Output:
c
c in common/radisw in RADISW.f (passed back to radiation scheme)
c      camt   - cloud amounts (locations specified by ktop/kbtm indices) 
c      cirab  - absorptivity of clouds in the near IR band (used in SW scheme)
c      cirrf  - reflectivity of clouds in the near IR band (used in SW scheme)
c      cuvrf - reflectivity of clouds in the visible band (used in SW scheme)
c      emcld - cloud emissivity (used in LW scheme)
c      kbtm  - index of (data level) pressure of cloud bottom (used in LW)
c      kbtmsw- index of (flux level) pressure of cloud bottom (used in SW)
c      ktop  - index of (data level) pressure of cloud top (used in LW)
c      ktopsw- index of (flux level) pressure of cloud top (used in SW)
c      nclds - no. clouds at each grid point
c
c in arguments
c      clat - cloud amount diagnostic (upside-down version of cfrac array)
c      clh - high level cloud diagnostic
c      cll - low level cloud diagnostic
c      clm - mid level cloud diagnostic
c
c******************************************************************************
 
      subroutine cloud2(cldoff,lg,ttg,qlg,qfg,cfrac,qccon,
     &                  cdrop,land,sigh,prf,dprf,cosz,     !Inputs
     &                  cll,clm,clh)                       !Outputs

      implicit none
C Global parameters
      include 'newmpar.h'
      include 'const_phys.h' !Input physical constants
      include 'cparams.h'    !Input cloud scheme parameters
      include 'kuocom.h'     ! ldr
      include 'params.h'     !Input model grid dimensions (modified params.h for CCAM)
      include 'rdparm.h'     !Input radiation scheme parameters
      include 'hcon.h'       !Input radiation physical constants

C Argument list
      logical cldoff
      integer lg
      real ttg(imax,l)
      real qlg(imax,l)
      real qfg(imax,l)
      real cfrac(imax,l)
      real qccon(imax,l)
      logical land(imax)
      real sigh(nl+1)
      real prf(imax,l)
      real dprf(imax,l)
      real cosz(imax)
      real cll(imax)
      real clm(imax)
      real clh(imax)

c      real Refflm(imax)
c      real cldliq(imax)


C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'radisw.h'     !Output various things (see above) and input coszro

C Global data blocks
      logical debug
      integer lgdebug,mgdebug,insdebug
      integer naerosol_i(2)
      data debug,lgdebug,mgdebug,insdebug /.false.,1,10106,1/
      data naerosol_i / 2*0 /



      integer jclb,nlow,nmid,ktied,kbgel,klcmc,klowrn,kmidrn
      real aftsea
      common/levdata/jclb,nlow,nmid,aftsea,ktied,kbgel
     &,klcmc,klowrn,kmidrn

C Local work arrays and variables
      real taul(imax,l), taui(imax,l) !Visible optical depths
      real Reffl(imax,l), Reffi(imax,l) !Effective radii (um)
      real Emw(imax,l), Abw(imax,l), Rew1(imax,l), Rew2(imax,l) !Water clouds
      real Emi(imax,l), Abi(imax,l), Rei1(imax,l), Rei2(imax,l) !Ice clouds
      real qlptot(imax),taultot(imax)
      real fice(imax,l),tau_sfac(imax,l)
      real rk(imax,nl),cdrop(imax,nl)

      integer i
      integer k
      integer mg
      integer nc
      integer ns

      real ab
      real beta
      real cfl
      real cldht
      real deltai
      real diffk
      real dz
      real em
      real eps
      real f1,f2,fcon,onem
      real qlpath
      real refac
      real refac1
      real refac2
      real re1
      real re2
      real rhoa
      real sigmai
      real tciwc
      real tclwc
      real temp_correction
      real tmid
      real trani
      real tranw
      real wice
      real wliq

C Local data, functions etc
      logical start
      data start /.true./
      save start

C Start code : ----------------------------------------------------------

c**** Set up levels for low/mid/high cloud counting (for cloud2.f)

      if(start)then

c.... "low" cloud up to half level closest to p/Ps=0.800 (~800mbs)
        k=1
        f2=0.800-sigh(k)
 22     k=k+1
        f1=f2
        f2=0.800-sigh(k)
        if((f1.le.0.0).and.(f2.gt.0.0))go to 30
        go to 22
 30     nlow=k-1
        if(f2.lt.abs(f1))nlow=k
c..   The FULL level is 1 below the half level indicator
        nlow=nlow-1

c.... "mid" cloud up to half level closest to p/Ps=0.400 (~400mbs)
        f2=0.400-sigh(k)
 32     k=k+1
        f1=f2
        f2=0.400-sigh(k)
        if((f1.le.0.0).and.(f2.gt.0.0))go to 40
        go to 32
 40     nmid=k-1
        if(f2.lt.abs(f1))nmid=k
c..   The FULL level is 1 below the half level indicator
        nmid=nmid-1
        start=.false.
      endif

        if(debug)then
          if(lg.eq.lgdebug)then
            ns=insdebug
            mg=mgdebug+(ns-1)*lon
            write(25,'(a,i1)')'Start cloud2'
            write(25,91)'ttg ',(ttg(mg,k),k=1,nl)
            write(25,9)'qlg ',(qlg(mg,k),k=1,nl)
            write(25,9)'qfg ',(qfg(mg,k),k=1,nl)
            write(25,91)'cfrac ',(cfrac(mg,k),k=1,nl)
            write(25,9)'qccon ',(qccon(mg,k),k=1,nl)
            write(25,9)'cdrop ',(cdrop(mg,k),k=1,nl)
            write(25,91)'prf ',(prf(mg,k),k=1,nl)
            write(25,91)'dprf ',(dprf(mg,k),k=1,nl)
            write(25,91)'cosz ',cosz(mg)
            write(25,*)
          endif
        endif

c***INITIALIZE THE CLOUD AND CLOUD INDEX FIELDS
c     Except for the ground layer (nc=1) the assumption is that
c     no cloud exists. also, for this purpose, the cloud index is set
c     at one (p=0)
c     Don't set cirrf, cuvrf, cirab at the surface because these are set
c     the albedo by radfs.
         camt(:,1)=zero
         emcld(:,1)=one
         ktop(:,1)=1
         kbtm(:,1)=1
         ktopsw(:,1)=1
         kbtmsw(:,1)=1
      do k=2,lp1
         camt(:,k)=zero
         emcld(:,k)=one
         cirrf(:,k)=0.
         cuvrf(:,k)=0.
         cirab(:,k)=0.
         ktop(:,k)=1
         kbtm(:,k)=1
         ktopsw(:,k)=1
         kbtmsw(:,k)=1
      enddo
c***NOW SET CLOUD AND CLOUD INDEX FIELDS DEPENDING ON THE NO. OF CLOUDS
      nc=1
c---FIRST, THE ground layer (nc=1)
         emcld(:,nc)=one
         camt(:,nc)=one
         ktop(:,nc)=lp1
         kbtm(:,nc)=lp1
         ktopsw(:,nc)=lp1
         kbtmsw(:,nc)=lp1

      If (cldoff) Then
            cll(:)=0.
            clm(:)=0.
            clh(:)=0.
            nclds(:)=0
      Else
          nclds(:)=0
          cll(:)=0.
          clm(:)=0.
          clh(:)=0.
c         cldliq(:)=0.
          qlptot(:)=0.
          taultot(:)=0.      
        
c Diagnose low, middle and high clouds; nlow,nmid are set up in initax.f
        
        do k=1,nlow
            cll(:)=cll(:)+cfrac(:,k)-cll(:)*cfrac(:,k)
        enddo
        do k=nlow+1,nmid
            clm(:)=clm(:)+cfrac(:,k)-clm(:)*cfrac(:,k)
        enddo
        do k=nmid+1,nl-1
            clh(:)=clh(:)+cfrac(:,k)-clh(:)*cfrac(:,k)
        enddo

c Set up rk and cdrop (now as cdso4 from radriv90.f)

c This is the Liu and Daum scheme for relative dispersion (Nature, 419, 580-581 and pers. comm.)

        do k=1,nl-1
          do mg=1,imax
c            eps = 1. - 0.7 * exp(-0.008e-6*cdrop(mg,k)) !Upper bound
c            eps = 1. - 0.7 * exp(-0.001e-6*cdrop(mg,k)) !Lower bound
            eps = 1. - 0.7 * exp(-0.003e-6*cdrop(mg,k)) !mid range
!           beta = ((1+2*eps**2)**2/(1+eps**2))**(1./3)
!           rk(mg,k)=1./(beta**3)
            rk(mg,k)= (1.+eps**2)/(1.+2.*eps**2)**2
!           beta=1./rk(mg,k)**(1./3)   ! calc for diag only

C***            write(26,'(2f12.3)')1.e-6*cdrop(mg,k),beta
C***            eps = 1. - 0.7 * exp(-0.008e-6*cdrop(mg,k)) !Upper bound
C***            beta = ((1+2*eps**2)**2/(1+eps**2))**(1./3)
C***            write(27,'(2f12.3)')1.e-6*cdrop(mg,k),beta
          enddo
        enddo

c        do k=1,nl
c         do mg=1,imax
c          if(land(mg))then 
c            rk(mg,k)=0.67 ! land
c	   else
c            rk(mg,k)=0.8  !sea
c          endif
c         enddo
c        enddo

c Define the emissivity (Em), and the SW properties (Re, Ab) for liquid (w)
c and ice (i) clouds respectively.
        
        do k=1,nl-1
          do mg=1,imax
            taul(mg,k)=0.
            taui(mg,k)=0.
            Reffl(mg,k)=0.
            Reffi(mg,k)=0.
            Emw(mg,k)=0.
            Emi(mg,k)=0.
            if(cfrac(mg,k).gt.0)then
              tau_sfac(mg,k)=1.
              fice(mg,k) = qfg(mg,k)/(qfg(mg,k)+qlg(mg,k))
            endif         !cfrac
          enddo
        enddo
              
c Liquid water clouds
        do k=1,nl-1
          do mg=1,imax
            if((cfrac(mg,k).gt.0).and.(qlg(mg,k).gt.1.0e-8))then
              rhoa=100*prf(mg,k)/(rdry*ttg(mg,k))
              dz=(dprf(mg,k)/prf(mg,k))*rdry*ttg(mg,k)/grav
              cfl=cfrac(mg,k)*(1-fice(mg,k))
              Wliq=rhoa*qlg(mg,k)/cfl     !kg/m^3
c              cldliq(mg)=cldliq(mg)+cfl-cldliq(mg)*cfl !Liquid cloud cover
                
c Reffl is the effective radius at the top of the cloud (calculated following
c Martin etal 1994, JAS 51, 1823-1842) due to the extra factor of 2 in the
c formula for reffl. Use mid cloud value of Reff for emissivity.
                
              Reffl(mg,k)=
     &           (3.*2.*Wliq/(4.*rhow*pi*rk(mg,k)*Cdrop(mg,k)))**(1./3.)
c              Reffl(mg,k)=
c     &             (3*Wliq/(4*rhow*pi*rk(mg,k)*Cdrop(mg,k)))**(1./3)
              qlpath=Wliq*dz
              taul(mg,k)=tau_sfac(mg,k)*1.5*qlpath/(rhow*Reffl(mg,k))
              qlptot(mg)=qlptot(mg)+qlpath*cfl
              taultot(mg)=taultot(mg)+taul(mg,k)*cfl

c Water cloud emissivity according to Martin Platt

C***          deltvl=taul(mg,k)*1.26 !Mult by 2^(1/3) so using mid cloud Reff
C***          deltal=min(0.4*deltvl, 45.) !IR optical depth for liq.
C***          if(deltvl.gt.0.4)then
C***            diffk=1.6
C***          else
C***            diffk=1.8
C***          endif
C***          Emw(mg,k) = 1.0 - exp(-diffk*deltal) !em of strat water clouds

c Or, water-cloud emissivity following the Sunshine scheme

              cldht=dz/1000.       !in km
              tclwc=Wliq*1000.     !in g/m**3
              tranw = exp( -1.66 * cldht * 50.885 * tclwc ** 0.769917)
              Emw(mg,k) = 1.0 - tranw    ! em is (1 - transmittance)
            endif         !cfrac
          enddo
        enddo
              
c Ice clouds : Choose scheme according to resolution
!      IF(lw.eq.22)THEN
       IF(ldr.gt.0)THEN  ! 1,2,3  corresponds to previous lw=22 option
        
	 refac1=0.85
        refac2=0.95
        do k=1,nl-1
          do mg=1,imax
            if((cfrac(mg,k).gt.0).and.(qfg(mg,k).gt.1.0e-8))then
              rhoa=100*prf(mg,k)/(rdry*ttg(mg,k))
              dz=(dprf(mg,k)/prf(mg,k))*rdry*ttg(mg,k)/grav
              Wice=rhoa*qfg(mg,k)/(cfrac(mg,k)*fice(mg,k)) !kg/m**3
!             Reffi(mg,k)=3.73e-4*Wice**0.216 !Lohmann et al. (1999)
              Reffi(mg,k)=min(150.e-6,3.73e-4*Wice**0.216) !Lohmann et al.(1999)
              taui(mg,k)=1.5*Wice*dz/(rhoice*Reffi(mg,k))
              deltai=min(0.5*taui(mg,k), 45.) !IR optical depth for ice.
              taui(mg,k)=tau_sfac(mg,k)*taui(mg,k)

c Ice-cloud emissivity following Platt

              if(taui(mg,k).gt.0.4)then
                diffk=1.6
              else
                diffk=1.8
              endif
              Emi(mg,k) = 1.0 - exp(-diffk*deltai)
            endif         !cfrac
          enddo
        enddo

       ELSE  ! i.e. for ldr = -1,-2,-3

        refac1=0.90
        refac2=1.00
        do k=1,nl-1
          do mg=1,imax
            if((cfrac(mg,k).gt.0).and.(qfg(mg,k).gt.1.0e-8))then
              rhoa=100*prf(mg,k)/(rdry*ttg(mg,k))
              dz=(dprf(mg,k)/prf(mg,k))*rdry*ttg(mg,k)/grav
              Wice=rhoa*qfg(mg,k)/(cfrac(mg,k)*fice(mg,k)) !kg/m**3
              sigmai = aice*Wice**bice !visible ext. coeff. for ice
              taui(mg,k)=sigmai*dz !visible opt. depth for ice
              Reffi(mg,k)=1.5*Wice*dz/(rhoice*taui(mg,k))
              taui(mg,k)=tau_sfac(mg,k)*taui(mg,k)

c Ice-cloud emissivity following the Sunshine scheme

              tmid=ttg(mg,k)-tfrz  !in celsius
              cldht=dz/1000.       !in km
              tciwc=Wice*1000.     !in g/m**3

              temp_correction=1.047E+00+tmid
     &        *(-9.13e-05+tmid
     &        *(2.026e-04-1.056e-05*tmid))
              temp_correction=max(1.0, temp_correction)
              trani = exp(-1.66*temp_correction * tciwc*cldht
     &              / (0.630689d-01+0.265874*tciwc))
c-- Limit ice cloud emissivities
c             trani=min(0.70,trani)
              if(k.gt.nlow) trani=min(0.70,trani)
              Emi(mg,k) = 1.0 - trani    ! em is (1 - transmittance)
            endif         !cfrac
          enddo
        enddo

       ENDIF  ! (ldr.gt.0) .. ELSE ..

c Calculate the effective radius of liquid water clouds seen from above

C***        do mg=1,imax
C***          if(taultot(mg).gt.0.)then
C***            Refflm(mg)=1.5*qlptot(mg)/(rhow*taultot(mg))
C***          else
C***            Refflm(mg)=0.
C***          endif
C***        enddo
          
c Calculate the SW cloud radiative properties for liquid water and ice clouds
c respectively, following Tony Slingo's (1989) Delta-Eddington scheme.

!       do k=1,nl
!         do mg=1,imax
!           if(cfrac(mg,k).gt.0.)then
C***              write(26,'(5f12.4)')cfrac(mg,k),taul(mg,k),taui(mg,k),
C***     &             coszro(mg),fice(mg,k)
C***              write(27,'(3f12.4)')
C***     &             cfrac(mg,k),1.e6*reffl(mg,k),1.e6*reffi(mg,k)
!           endif
!        enddo
!      enddo

        call slingo (Reffl, taul, cosz, !inputs
     &       Rew1, Rew2, Abw )  !outputs
        
        call slingi (Reffi, taui, cosz, !inputs
     &       Rei1, Rei2, Abi )  !outputs

        onem=1.-1.e-6   ! to avoid possible later 0/0
        do k=1,nl-1
          do mg=1,imax
            if(cfrac(mg,k).gt.0.)then
              fcon=min(1.,qccon(mg,k)/(qlg(mg,k)+qfg(mg,k)))
c Original refac :
c             refac=0.7*fcon+0.9*(1-fcon)
c Mk3 with no direct aerosol effect :
c             refac=0.7*fcon+0.85*(1-fcon)
c Mk3 with direct aerosol effect :
              refac=refac1*fcon+refac2*(1.-fcon)

              Rei1(mg,k)=min(refac*Rei1(mg,k),onem)
              Rei2(mg,k)=min(refac*Rei2(mg,k),onem-2.*Abi(mg,k))
              Rew1(mg,k)=min(refac*Rew1(mg,k),onem)
              Rew2(mg,k)=min(refac*Rew2(mg,k),onem-2.*Abw(mg,k))
            endif
          enddo
        enddo

c Weight cloud properties by liquid/ice fraction
        
        do k=1,nl-1
          do mg=1,imax
            if(cfrac(mg,k).gt.0.)then

              Re1 = fice(mg,k)*Rei1(mg,k) + (1.-fice(mg,k))*Rew1(mg,k)
              Re2 = fice(mg,k)*Rei2(mg,k) + (1.-fice(mg,k))*Rew2(mg,k)
              Em = fice(mg,k)*Emi(mg,k) + (1.-fice(mg,k))*Emw(mg,k)
!             if(qlg(mg,k).gt.0)write(26,'(2g12.3)')qlg(mg,k),rei1(mg,k)
!             if(qfg(mg,k).gt.0)write(27,'(2g12.3)')qfg(mg,k),emi(mg,k)

c             if(prf(mg,k).gt.800.) Em = 1.
              Ab = fice(mg,k)*Abi(mg,k) + (1-fice(mg,k))*Abw(mg,k)

!             write(30,*)'mg, nclds(mg) ',mg, nclds(mg)
              nclds(mg)=nclds(mg)+1
              nc=nclds(mg)+1
              camt(mg,nc)=cfrac(mg,k)
              ktop(mg,nc)=nlp-k
              kbtm(mg,nc)=nlp-k
              emcld(mg,nc)=Em
              ktopsw(mg,nc)=nlp-k
              kbtmsw(mg,nc)=nlp-k+1
              cuvrf(mg,nc)=Re1
              cirrf(mg,nc)=Re2
              cirab(mg,nc)=2*Ab
            endif
          enddo
        enddo

        if(debug)then
          if(lg.eq.lgdebug)then
            ns=insdebug
            mg=mgdebug+(ns-1)*lon
            write(25,'(a,i1)')'After cloud2'
            write(25,9)'cdrop ',(cdrop(mg,k),k=1,nl)
            write(25,91)'Rew ',((Rew1(mg,k)+Rew2(mg,k))/2,k=1,nl)
            write(25,91)'Rei ',((Rei1(mg,k)+Rei2(mg,k))/2,k=1,nl)
            write(25,91)'Abw ',(Abw(mg,k),k=1,nl)
            write(25,91)'Abi ',(Abi(mg,k),k=1,nl)
            write(25,91)'Emw ',(Emw(mg,k),k=1,nl)
            write(25,91)'Emi ',(Emi(mg,k),k=1,nl)
            write(25,9)'Reffl ',(Reffl(mg,k),k=1,nl)
            write(25,9)'Reffi ',(Reffi(mg,k),k=1,nl)
            write(25,*)
          endif
        endif
 1      format(3(a,f10.3))
 9      format(a,30g10.3)
 91     format(a,30f10.3)
      End If   !cldoff
      

      return
      end

c******************************************************************************

c Calculate SW radiative properties for water clouds using delta-Eddington
c scheme as given by Slingo (1989) JAS 46, 1419-1427.
c Coefficients di and fi modified to use Reff in SI units.

      subroutine slingo(reff, tau, mu0,       !inputs
     &                  refl1, refl2, abso )  !outputs

      implicit none
C Global parameters
      include 'newmpar.h'
      include 'rdparm.h'
      integer nbands
      parameter (nbands=4)

C Argument list
      real reff(imax,kl)
      real tau(imax,kl)
      real mu0(imax)
      real refl1(imax,kl)
      real refl2(imax,kl)
      real abso(imax,kl)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks

C Local work arrays and variables
      double precision exparg,denom,epsilon,omega,f,omwf

      integer i
      integer k
      integer mg

      real absband
      real alpha1
      real alpha2
      real alpha3
      real alpha4
      real beta
      real beta0
      real e
      real gam1
      real gam2
      real gi
      real rdif
      real rm
      real rdir
      real tdb
      real tdif
      real tdir
      real ttot
      real u1
      real u2

C Local data, functions etc
      real ci(nbands)
      data ci / -5.62e-8, -6.94e-6, 4.64e-4,  2.01e-1  /

      real di(nbands)
c      data di / 1.63e-7,  2.35e-5,  1.24e-3,  7.56e-3  /
      data di / 1.63e-1,  2.35e+1,  1.24e+3,  7.56e+3  / !Si units

      real ei(nbands)
      data ei / 0.829,    0.794,    0.754,    0.826    /

      real fi(nbands)
c      data fi / 2.482e-3, 4.226e-3, 6.560e-3, 4.353e-3 /
      data fi / 2.482e+3, 4.226e+3, 6.560e+3, 4.353e+3 / !SI units

      real wi(nbands)
      data wi / 0.459760, 0.326158, 0.180608, 0.033474 /
 
      real wi24
      data wi24 / 0.540240 / !Sum of wi, i=2 to 4

C Start code : ----------------------------------------------------------

      do k=1,kl
        do mg=1,imax
          refl1(mg,k)=0.
          refl2(mg,k)=0.
          abso(mg,k)=0.
        enddo
      enddo

      do i=1,nbands
        do k=1,kl-1
          do mg=1,imax
            if(tau(mg,k).gt.0..and.mu0(mg).gt.0.)then
               reff(mg,k)=min(20.e-6,max(4.e-6,reff(mg,k)))
               omega=1-(ci(i)+di(i)*reff(mg,k))     
               gi=ei(i)+fi(i)*reff(mg,k)
               beta0=(3./7.)*(1-gi)
               beta=0.5-0.75*mu0(mg)*gi/(1+gi)
               f=gi**2
               U1=7./4.
               U2=(7./4)*(1.-(1-omega)/(7*omega*beta0))
               alpha1=U1*(1.-omega*(1-beta0))
               alpha2=U2*omega*beta0
               alpha3=(1-f)*omega*beta
               alpha4=(1-f)*omega*(1-beta)
               epsilon=sqrt(alpha1**2-alpha2**2)
               rM=alpha2/(alpha1+epsilon)
               E=exp(-epsilon*tau(mg,k))
               omwf=1-omega*f
               denom=omwf**2-epsilon**2*mu0(mg)**2
               gam1=(omwf*alpha3-mu0(mg)*(alpha1*alpha3+alpha2*alpha4))
     &               /denom 
               gam2=(-omwf*alpha4-mu0(mg)*(alpha1*alpha4+alpha2*alpha3))
     &               /denom
               exparg=dmin1(70.0d0,omwf*tau(mg,k)/mu0(mg))
               Tdb=exp(-exparg)
               Rdif=rM*(1-E**2)/(1-(E*rM)**2)
               Tdif=E*(1-rM**2)/(1-(E*rM)**2)
               Rdir=-gam2*Rdif-gam1*Tdb*Tdif+gam1
               Tdir=-gam2*Tdif-gam1*Tdb*Rdif+gam2*Tdb
               Ttot=Tdb+Tdir
               Absband=1-Rdir-Ttot
               Absband=max(0., Absband) !Needed for 32 bit
c               if(absband.gt.1..or.absband.lt.0)then
c                 print*,'Warning slingo: band, abs =',i,absband
c               endif
               abso(mg,k)=abso(mg,k)+Absband*wi(i)
               if(i.eq.1)then
                 refl1(mg,k)=Rdir
               else
                 refl2(mg,k)=refl2(mg,k)+Rdir*wi(i)/wi24
               endif
             endif
           enddo
         enddo
       enddo

       return
       end

c******************************************************************************

c Slingo type scheme for ice cloud SW properties (similar to water clouds)
c Single scattering albedo follows Francis etal (1994) QJRMS 120, 809--848.
c Lambdas modified to use Reff in SI units.

      subroutine slingi(reff, tau, mu0,       !inputs
     &                  refl1, refl2, abso )  !outputs

      implicit none
C Global parameters
      include 'newmpar.h'
      include 'const_phys.h' !Input physical constants
      include 'rdparm.h'
      integer nbands
      parameter (nbands=4)

C Argument list
      real reff(imax,kl)
      real tau(imax,kl)
      real mu0(imax)
      real refl1(imax,kl)
      real refl2(imax,kl)
      real abso(imax,kl)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks

C Local work arrays and variables
      double precision vm !For 32 bit
      double precision exparg,denom,epsilon,omega,f,omwf

      integer i
      integer k
      integer mg

      real absband
      real alpha1
      real alpha2
      real alpha3
      real alpha4
      real beta
      real beta0
      real e
      real gam1
      real gam2
      real gi
      real rdif
      real rm
      real rdir
      real tdb
      real tdif
      real tdir
      real ttot
      real u1
      real u2

C Local data, functions etc
      real lambda(nbands)
      data lambda / 0.470e-6, 0.940e-6, 1.785e-6, 3.190e-6 /

      real nprime(nbands)
      data nprime / 5.31e-9,  8.74e-7,  3.13e-4,  1.07e-1  /

      real wi(nbands)
      data wi     / 0.459760, 0.326158, 0.180608, 0.033474 /

      real wi24
      data wi24   / 0.540240 / !Sum of wi, i=2 to 4

C Start code : ----------------------------------------------------------

      do k=1,kl
        do mg=1,imax
          refl1(mg,k)=0.
          refl2(mg,k)=0.
          abso(mg,k)=0.
        enddo
      enddo

      do i=1,nbands
        do k=1,kl-1
          do mg=1,imax
            if(tau(mg,k).gt.0..and.mu0(mg).gt.0.)then
               vm=2*pi*reff(mg,k)*nprime(i)/lambda(i)
               omega=0.5+(vm+3)/(6*(vm+1)**3)
               gi=0.8  !Constant asymmetry parameter
               beta0=(3./7.)*(1-gi)
               beta=0.5-0.75*mu0(mg)*gi/(1+gi)
               f=gi**2
               U1=7./4.
               U2=(7./4)*(1.-(1-omega)/(7*omega*beta0))
               alpha1=U1*(1.-omega*(1-beta0))
               alpha2=U2*omega*beta0
               alpha3=(1-f)*omega*beta
               alpha4=(1-f)*omega*(1-beta)
               epsilon=sqrt(alpha1**2-alpha2**2)
               rM=alpha2/(alpha1+epsilon)
               E=exp(-epsilon*tau(mg,k))
               omwf=1-omega*f
               denom=omwf**2-epsilon**2*mu0(mg)**2
               gam1=(omwf*alpha3-mu0(mg)*(alpha1*alpha3+alpha2*alpha4))
     &               /denom
               gam2=(-omwf*alpha4-mu0(mg)*(alpha1*alpha4+alpha2*alpha3))
     &               /denom
               exparg=dmin1(70.0d0,omwf*tau(mg,k)/mu0(mg))
               Tdb=exp(-exparg)
               Rdif=rM*(1-E**2)/(1-(E*rM)**2)
               Tdif=E*(1-rM**2)/(1-(E*rM)**2)
               Rdir=-gam2*Rdif-gam1*Tdb*Tdif+gam1
               Tdir=-gam2*Tdif-gam1*Tdb*Rdif+gam2*Tdb
               Ttot=Tdb+Tdir
               Absband=1-Rdir-Ttot
               Absband=max(0., Absband) !Needed for 32 bit
c               if(absband.gt.1..or.absband.lt.0)then
c                 print*,'Warning slingi: band, abs =',i,absband
c               endif
               abso(mg,k)=abso(mg,k)+Absband*wi(i)
               if(i.eq.1)then
                 refl1(mg,k)=Rdir
               else
                 refl2(mg,k)=refl2(mg,k)+Rdir*wi(i)/wi24
               endif
             endif
           enddo
         enddo
       enddo

       return
       end
