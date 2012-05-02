      ! CCAM boundary layer turbulent mixing routines

      ! Currently local Ri and prognostic TKE-eps schemes are supported.
      ! Also, options for non-local counter gradient terms are included.
      ! Local Ri supports Gelyen and Tiedtke schemes for shallow
      ! convection, whereas the TKE-eps follows the EDMF approach where
      ! shallow convection is represented in the mass flux terms.

      ! nvmix=3  Local Ri mixing
      ! nvmix=6  Prognostic TKE-eps tubulence closure

      ! nlocal=0 No counter gradient term
      ! nlocal=6 Holtslag and Boville non-local term
      ! nlocal=7 Mass flux based counter gradient
      
      ! kscmom   0 off, 1 turns on shal_conv momentum (usual)

      !--------------------------------------------------------------
      ! Control subroutine for vertical mixing
      subroutine vertmix(iaero)

      use aerosolldr                      ! LDR prognostic aerosols
      use arrays_m                        ! Atmosphere dyamics prognostic arrays
      use cc_mpi, only : mydiag, myid     ! CC MPI routines
      use cfrac_m                         ! Cloud fraction
      use diag_m                          ! Diagnostic routines
      use extraout_m                      ! Additional diagnostics
      use indices_m                       ! Grid index arrays
      use kuocomb_m                       ! JLM convection
      use liqwpar_m                       ! Cloud water mixing ratios
      use map_m                           ! Grid map arrays
      use morepbl_m                       ! Additional boundary layer diagnostics
      use nharrs_m                        ! Non-hydrostatic atmosphere arrays
      use nlin_m, tmnht => un, at => un   ! Atmosphere non-linear dynamics
      use pbl_m                           ! Boundary layer arrays
      use permsurf_m                      ! Fixed surface arrays
      use savuvt_m                        ! Saved dynamic arrays
      use screen_m                        ! Screen level diagnostics
      use sigs_m                          ! Atmosphere sigma levels
      use soil_m                          ! Soil and surface data
      use tkeeps                          ! TKE-EPS boundary layer
      use tracers_m                       ! Tracer data
      use trvmix, only : tracervmix       ! Tracer mixing routines
      use work2_m, only : rho             ! Diagnostic arrays
      
      implicit none
      
      include 'newmpar.h'                 ! Grid parameters
      include 'const_phys.h'              ! Physical constants
      include 'dates.h'                   ! Date data
      include 'establ.h'                  ! Liquid saturation function
      include 'kuocom.h'                  ! Convection parameters
      include 'mpif.h'                    ! MPI parameters
      include 'parm.h'                    ! Model configuration

      integer, parameter :: ntest=0
      integer, parameter :: ndvmod=0    ! 0 default, 1+ for dvmod tests
      integer kcl_top,iaero,l,iq,k,iqmax,ierr                  
      integer, save :: kscbase,ksctop
      integer, dimension(ifull) :: kbase,ktop
      real, parameter :: bprm=4.7,cm=7.4,ch=5.3,amxlsq=100.   ! coefficients for Louis scheme
      real, parameter :: vkar3=0.35,vkar4=0.4,bprmj=5.,cmj=5. ! coefficients for Louis scheme
      real, parameter :: chj=2.6                              ! coefficients for Louis scheme
      real condrag,rong,diff,denma,denha,tsp,epart,esp,delsig
      real conflux,qbas,xold,delthet_old,rhskp,rhsk,diffmax
      real w1,w2,al,betaq,betac,betat,pk,qc,dqsdt,fice,es,delta
      real rlog12,rlogh1,rlogs1,rlogs2
      real, dimension(ifull,kl) :: betatt,betaqt,rhs,delthet,thebas
      real, dimension(ifull,kl) :: cu,thee,qs,uav,vav,au,ct,zh,gt
      real, dimension(ifull,kl) :: guv,ri,rkm,rkh,rk_shal,zg
      real, dimension(ifull,kl-1) :: zgh
      real, dimension(ifull) :: dqtot,csq,dvmod,dz,dzr,fm,fh,sqmxl
      real, dimension(ifull) :: x,zhv,theeb,sigsp,wt0,wq0
      real, dimension(kl) :: sighkap,sigkap,delons,delh
      real, dimension(:), allocatable, save :: prcpv

      kcl_top=kl-2 ! maximum level for cloud top (conjob & vertmix)
      if (.not.allocated(prcpv)) allocate(prcpv(kl))

      rong=rdry/grav
      do k=1,kl-1
       sighkap(k)=sigmh(k+1)**(-roncp)
       delons(k)=rong *((sig(k+1)-sig(k))/sigmh(k+1))
      enddo      ! k loop
      do k=1,kl
       delh(k)=-rong *dsig(k)/sig(k)  ! sign of delh defined so always +ve
       sigkap(k)=sig(k)**(-roncp)
      enddo      ! k loop
      if(diag.or.ntest>=1)then
        call maxmin(u,'%u',ktau,1.,kl)
        call maxmin(v,'%v',ktau,1.,kl)
        call maxmin(t,'%t',ktau,1.,kl)
        call maxmin(qg,'qg',ktau,1.e3,kl)     
        call MPI_Barrier( MPI_COMM_WORLD, ierr ) ! stop others going past
      if(mydiag)then
          write(6,*)'sig ',sig
          write(6,*)'dsig ',dsig
          write(6,*)'delh ',delh
          write(6,*)'in vertmix'
          write (6,"('uin ',9f8.3/4x,9f8.3)") u(idjd,:) 
          write (6,"('vin ',9f8.3/4x,9f8.3)") v(idjd,:) 
        endif
      endif
      rlogs1=log(sig(1))
      rlogs2=log(sig(2))
      rlogh1=log(sigmh(2))
      rlog12=1./(rlogs1-rlogs2)
      tmnht(1:ifull,1)=(t(1:ifull,2)*rlogs1-t(1:ifull,1)*rlogs2+
     &           (t(1:ifull,1)-t(1:ifull,2))*rlogh1)*rlog12
!     n.b. an approximate zh (in m) is quite adequate for this routine
      zh(1:ifull,1)=t(1:ifull,1)*delh(1)
      do k=2,kl-1
       do iq=1,ifull
        zh(iq,k)=zh(iq,k-1)+t(iq,k)*delh(k)
        tmnht(iq,k) =(t(iq,k)+t(iq,k+1))*.5
       enddo     ! iq loop
      enddo      !  k loop
      do k=1,kl
        rhs(:,k)=t(1:ifull,k)*sigkap(k)  ! rhs is theta here
      enddo      !  k loop
      if(nmaxpr==1.and.mydiag)
     &  write (6,"('thet_in',9f8.3/7x,9f8.3)") rhs(idjd,:)


      if (nvmix.ne.6) then

      !--------------------------------------------------------------
      ! JLM's local Ri scheme

c Pre-calculate the buoyancy parameters if using qcloud scheme.
c Follow Smith's (1990) notation; gam() is HBG's notation for (L/cp)dqsdt.
c The factor of (1/sigkap)=T/theta in betatt differs from Smith's formulation
c because we use theta derivative rather than (dry static energy)/cp.
      if(nvmix>0.and.nvmix<4)then
        delta=1./epsil-1.  ! i.e. 1/.622 -1., i.e. .6077
        do k=1,kl
         if(sig(k)>.8)then ! change made 17/1/06
           do iq=1,ifull
            es=establ(t(iq,k))
            pk=ps(iq)*sig(k)
            qs(iq,k)=.622*es/(pk-es)  
            dqsdt=qs(iq,k)*pk*(hl/rvap)/(t(iq,k)**2*(pk-es))
            rhs(iq,k)=rhs(iq,k) !Convert to thetal - used only to calc Ri
     &                  -(hlcp*qlg(iq,k)+hlscp*qfg(iq,k))*sigkap(k)   
            betat=1./t(iq,k)
            qc=qlg(iq,k)+qfg(iq,k)
            fice=qfg(iq,k)/max(qc,1.e-12)
            betaq=delta/(1.+delta*qg(iq,k)-qc)
            al=1./(1.+hlcp*dqsdt)
            betac=cfrac(iq,k)
     &           *al * ((hlcp+fice*hlfcp)*betat - betaq/(1.-epsil) )
            betatt(iq,k)=(betat-dqsdt*betac)/sigkap(k)  !Beta_t_tilde
            betaqt(iq,k)=betaq+betac                    !Beta_q_tilde
           enddo   ! iq loop
         else  ! i.e. (sig(k)<.8)
           do iq=1,ifull
            es=establ(t(iq,k))
            pk=ps(iq)*sig(k)
            qs(iq,k)=.622*es/max(1.,pk-es)  ! still need qs(); max for k=kl
            betat=1./t(iq,k)
c           qc=qlg(iq,k)+qfg(iq,k)
c           betaq=delta/(1.+delta*qg(iq,k)-qc)
            betaq=delta/(1.+delta*qg(iq,k))
            betatt(iq,k)=betat/sigkap(k)    !Beta_t_tilde
            betaqt(iq,k)=betaq              !Beta_q_tilde
           enddo   ! iq loop
         endif  ! (sig(k)>.8)
         if(diag.and.mydiag)then
            iq=idjd
            dqsdt=qs(iq,k)*pk*(hl/rvap)/(t(iq,k)**2*max(pk-es,1.))
            betat=1./t(iq,k)
            qc=qlg(iq,k)+qfg(iq,k)
            fice=qfg(iq,k)/max(qc,1.e-12)
            betaq=delta/(1+delta*qg(iq,k)-qc)
!           al=1./(1.+gam(iq,k))
            al=1./(1.+hlcp*dqsdt)
            betac=cfrac(iq,k)
     &           *al * ((hlcp+fice*hlfcp)*betat - betaq/(1.-epsil) )
            write(6,*)'k,qg,qs,qlg,qfg,qc ',
     &               k,qg(iq,k),qs(iq,k),qlg(iq,k),qfg(iq,k),qc
            write(6,*)'t,rhs,cfrac,fice ',
     &               t(iq,k),rhs(iq,k),cfrac(iq,k),fice
            write(6,*)'al,delta,betaq,betac ',al,delta,betaq,betac
            write(6,*)'betat,betatt,betaqt ',
     &               betat,betatt(iq,k),betaqt(iq,k)
         endif   ! (ntest==2)
        enddo    !  k loop
      else       ! other nvmix values (0 or 4+) still need qs()
        do k=1,kl
         do iq=1,ifull
          es=establ(t(iq,k))
          qs(iq,k)=.622*es/max(1.,ps(iq)*sig(k)-es)  ! max for k=kl
         enddo   ! iq loop
        enddo    !  k loop
      endif      ! (nvmix>0.and.nvmix<4)

      do k=1,kl-1
       delthet(:,k)=rhs(:,k+1)-rhs(:,k)  ! rhs is theta or thetal here
      enddo      !  k loop

      if(ktau==1.and.ksc.ne.0)then
!       set ksctop for shallow convection
        ksctop=1    ! ksctop will be first level below sigkcst
        do while(sig(ksctop+1)>sigksct)  !  e.g. sigksct=.75
         ksctop=ksctop+1
        enddo
        kscbase=1  ! kscbase will be first level above sigkcsb
        do while(sig(kscbase)>sigkscb.and.sigkscb>0.) ! e.g. sigkscb=.99
         kscbase=kscbase+1
        enddo
        if ( myid == 0 ) then
        write(6,*)'For shallow convection:'
        write(6,*)'ksc,kscbase,ksctop,kscsea ',
     &           ksc,kscbase,ksctop,kscsea
        write (6,"(' sigkscb,sigksct,tied_con,tied_over,tied_rh:',
     &       5f8.3)")sigkscb,sigksct,tied_con,tied_over,tied_rh
        end if
        do k=1,kl
         prcpv(k)=sig(k)**(-roncp)
        enddo  ! k loop
      endif    ! (ktau==1.and.ksc.ne.0)

c     ****** section for Geleyn shallow convection; others moved lower****
      if(ksc==-99)then
        do k=kscbase,ksctop    ! new usage of ksc thu  02-17-2000
         do iq=1,ifull
          delthet(iq,k)=delthet(iq,k)-hlcp*max(0.,
     &                 qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k) )
         enddo  ! iq loop
        enddo   !  k loop
      endif     ! (ksc==-99)
      if(ksc==-98)then    ! modified Geleyn Jan '08
        do k=kscbase,ksctop    ! new usage of ksc thu  02-17-2000
         do iq=1,ifull
          if(qg(iq,k)>tied_rh*qs(iq,k).or.
     &       qg(iq,k+1)>tied_rh*qs(iq,k+1))then
             delthet(iq,k)=delthet(iq,k)-hlcp*max(0.,
     &                     qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k) )
          endif
         enddo  ! iq loop
        enddo   !  k loop
      endif     ! (ksc==-98)
      if(ksc==-97)then    ! modified Geleyn Jan '08
        do k=kscbase,ksctop    ! new usage of ksc thu  02-17-2000
         do iq=1,ifull
          if(qg(iq,k)>tied_rh*qs(iq,k))then
            delthet(iq,k)=delthet(iq,k)-hlcp*max(0.,
     &                    qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k) )
          endif
         enddo  ! iq loop
        enddo   !  k loop
      endif     ! (ksc==-97)
      if(ksc==-96)then   ! combination of Geleyn and jlm 83 (Feb 08)
        do k=1,ksctop    
         do iq=1,ifull
          if(k<ktsav(iq).and.k>=kbsav(iq).and.condc(iq)==0.)then  
            delthet(iq,k)=delthet(iq,k)-hlcp*max(0.,
     &                 qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k) )
            write(6,*)'-96 iq,k,kbsav,ktsav ',iq,k,kbsav(iq),ktsav(iq),
     &      hlcp*(qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k)),delthet(iq,k)
          endif 
         enddo  ! iq loop
        enddo   !  k loop
      endif     ! (ksc==-96)
      if(ksc==-95)then ! same as -99 but has tied_rh (e.g. .75) factor
c       kshal(:)=0
        do k=kscbase,ksctop     
         do iq=1,ifull
          delthet(iq,k)=delthet(iq,k)-tied_rh*hlcp*max(0.,
     &                 (qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k)) )
         enddo  ! iq loop
        enddo   !  k loop
      endif     ! (ksc==-95)
      if(ksc==-94)then   ! combination of Geleyn and jlm 83 (Feb 08)
        do k=1,ksctop    
         do iq=1,ifull
          if(k<ktsav(iq).and.k>=kbsav(iq).and.condc(iq)==0.)then  
            delthet(iq,k)=0.
            write(6,*)'-94 iq,k,kbsav,ktsav ',iq,k,kbsav(iq),ktsav(iq),
     &      hlcp*(qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k)),delthet(iq,k)
          endif 
         enddo  ! iq loop
        enddo   !  k loop
      endif     ! (ksc==-94)
      if(ksc==-93)then ! single-layer (pblh) version of -95
        do iq=1,ifull
         do k=kscbase,kl/2
          if(zh(iq,k)<pblh(iq).and.zh(iq,k+1)>pblh(iq))then
            delthet(iq,k)=delthet(iq,k)-tied_rh*hlcp*max(0.,
     &                   (qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k)) )
           endif
         enddo  ! k loop 
        enddo   ! iq loop
      endif     ! (ksc==-93)
      if(ksc==-92)then ! capped-by-pblh version of -95
        do iq=1,ifull
         do k=kscbase,kl/2
          if(zh(iq,k)<pblh(iq))then
            delthet(iq,k)=delthet(iq,k)-tied_rh*hlcp*max(0.,
     &                   (qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k)) )
           endif
         enddo  ! k loop 
        enddo   ! iq loop
      endif     ! (ksc==-92)
      if(ksc==-91)then ! capped-by-pblh (anywhere in layer) version of -95
        do iq=1,ifull
         do k=2,kl/2
          if(zh(iq,k-1)<pblh(iq))then
            delthet(iq,k)=delthet(iq,k)-tied_rh*hlcp*max(0.,
     &                   (qs(iq,k+1)-qg(iq,k+1)-qs(iq,k)+qg(iq,k)) )
           endif
         enddo  ! k loop 
        enddo   ! iq loop
      endif     ! (ksc==-91)
c     ********* end of Geleyn shallow convection section ****************

!     following now defined in vertmix (don't need to pass from sflux)
      uav(1:ifull,:)=av_vmod*u(1:ifull,:)+(1.-av_vmod)*savu(1:ifull,:) 
      vav(1:ifull,:)=av_vmod*v(1:ifull,:)+(1.-av_vmod)*savv(1:ifull,:) 
      do k=1,kl-1
       do iq=1,ifull
        dz(iq) =-tmnht(iq,k)*delons(k)  ! this is z(k+1)-z(k)
        dzr(iq)=1./dz(iq)
        zhv(iq)=1./zh(iq,k)
        if(ndvmod==0)then
          dvmod(iq)=sqrt( (uav(iq,k+1)-uav(iq,k))**2
     &                   +(vav(iq,k+1)-vav(iq,k))**2 )
        else
          dvmod(iq)=ndvmod  ! just for tests
        endif    ! (ndvmod==0)
       enddo ! iq loop

c      x is bulk ri *(dvmod **2); used to calc ri(k), rkh(k) etc
       if(nvmix>0.and.nvmix<4)then  ! new one allowing for cloudy air
!        usually nvmix=3       
         if(sig(k)>.8)then ! change made 17/1/06
           dqtot(:)=qg(1:ifull,k+1)+qlg(1:ifull,k+1)+qfg(1:ifull,k+1)
     &            -(qg(1:ifull,k)  +qlg(1:ifull,k)  +qfg(1:ifull,k))
         else
           dqtot(:)=0.
         endif
         if(nvmix==2)then !  jlm May '05
           do iq=1,ifull
            x(iq)=grav*dz(iq)*(
     &         (min(betatt(iq,k),betatt(iq,k+1)))*delthet(iq,k) +
     &         (max(betaqt(iq,k),betaqt(iq,k+1)))*dqtot(iq) )
           enddo ! iq loop	
         else    ! i.e. nvmix=1 or 3
           if(nvmix==1)w1=dsig(k+1)/(dsig(k)+dsig(k+1)) 
           if(nvmix==3)w1=1.    !weight for lower level  usual           
           w2=1.-w1             !weight for upper level
           do iq=1,ifull
            x(iq)=grav*dz(iq)*(
     &           (w1*betatt(iq,k)+w2*betatt(iq,k+1))*delthet(iq,k) +
     &           (w1*betaqt(iq,k)+w2*betaqt(iq,k+1))*dqtot(iq) )
           enddo ! iq loop	 
         endif  !  (nvmix==2) .. else ..
         if(ntest==4.and.k<=9.and.ktau==ntau)then
           diffmax=0.
           do iq=1,ifull
            rhsk=t(iq,k)*sigkap(k)
            rhskp=t(iq,k+1)*sigkap(k+1)
            delthet_old=rhs(iq,k+1)-rhs(iq,k)
            xold=grav*dz(iq)*(delthet_old/(tmnht(iq,k)*sighkap(k))
     &           +.61*(qg(iq,k+1)-qg(iq,k)))
            diff=abs(xold-x(iq))
            if(diff>diffmax)then
              diffmax=diff
              iqmax=iq
            endif
            write(47,'(3g13.4,i7,i4)') xold,x(iq),diff,iq,k
           enddo
           write(6,*)'k,iqmax,diffmax ',k,iqmax,diffmax
         endif   ! (ntest==4.and.k<=9.and.ktau==ntau)
         rhs(:,k)=t(1:ifull,k)*sigkap(k)   !need to re-set theta for nvmix=1-3
       elseif(nvmix==5)then        ! non-cloudy x with qfg, qlg
          x(:)=grav*dz(:)*(delthet(:,k)/(tmnht(:,k)*sighkap(k))
     &             +.61*(qg(1:ifull,k+1)-qg(1:ifull,k))
     &                  -qfg(1:ifull,k+1)-qlg(1:ifull,k+1)+
     &                   qfg(1:ifull,k)+qlg(1:ifull,k) )
       else                 ! original non-cloudy x, nvmix=4
          x(:)=grav*dz(:)*(delthet(:,k)/(tmnht(:,k)*sighkap(k))
     &            +.61*(qg(1:ifull,k+1)-qg(1:ifull,k)))
       endif     ! (nvmix>0.and.nvmix<4) .. else ..

c       fm and fh denote f(Louis style)*dvmod
c       nb. an error exists in the Louis expression for c; this is corrected
c       in the current code
        csq(:) = zhv(:)*(((1.+dz(:)*zhv(:))**(1./3.)-1.)*dzr(:))**3

c        newest code, stable same as csiro9 here (originally for nvmix=4)
         do iq=1,ifull
          sqmxl(iq)=(vkar4*zh(iq,k)/(1.+vkar4*zh(iq,k)/amxlsq))**2
          dvmod(iq)=max( dvmod(iq) , 1. )
          ri(iq,k)=x(iq)/dvmod(iq)**2
          if(ri(iq,k)< 0.)then  ! unstable case
c           first do momentum
            denma=dvmod(iq)+cmj*( 2.*bprmj*sqmxl(iq)*
     &                            sqrt(-x(iq)*csq(iq)) )
            fm(iq)=dvmod(iq)-(2.*bprmj *x(iq))/denma
c           now heat
            denha=dvmod(iq)+chj*( 2.*bprmj*sqmxl(iq)*
     &                            sqrt(-x(iq)*csq(iq)) )
            fh(iq)=dvmod(iq)-(2.*bprmj *x(iq))/denha
          else                     ! stable case
c           the following is the original Louis stable formula
            fm(iq)=dvmod(iq)/(1.+4.7*ri(iq,k))**2
            fh(iq)=fm(iq)
          endif
         enddo   ! iq loop

       if(nvmix<0)then   ! use quasi-neutral mixing for test runs only
         do iq=1,ifull
          fm(iq)=1.
          fh(iq)=1.
         enddo   ! iq loop
         if(k<kl/2)then
           do iq=1,ifull
            if(land(iq))fh(iq)=abs(nvmix)
           enddo   ! iq loop
         endif
       endif     !  (nvmix<0)

c      calculate k's, guv and gt
c      nonlocal scheme usually applied to temperature and moisture, not momentum
c      (i.e. local scheme is applied to momentum for nlocal=0,1)
       if(nvmix.ne.0)then    ! use nvmix=0 only for special test runs
         do iq=1,ifull
          rkm(iq,k)=fm(iq)*sqmxl(iq)*dzr(iq)
          rkh(iq,k)=fh(iq)*sqmxl(iq)*dzr(iq)
         enddo   ! iq loop
       else
         rkm(:,:)=0.
         rkh(:,:)=0.
       endif     ! (nvmix.ne.0)

       if((diag.or.ntest>=1).and.mydiag)then
        iq=idjd
         write(6,*)'k,dt,delsig,sqmxl,dzr ',
     &            k,dt,delsig,sqmxl(idjd),dzr(idjd)
         write(6,*)'k,t,t+,ps ',k,t(idjd,k),t(idjd,k+1),ps(idjd)
         write(6,*)'k,qg,qg+ ',k,qg(idjd,k),qg(idjd,k+1)
         write(6,*)'k,qs,qs+ ',k,qs(idjd,k),qs(idjd,k+1)
         es=establ(t(idjd,k))
         esp=establ(t(idjd,k+1))
         write(6,*)'k,es,es+,delthet ',k,es,esp,delthet(idjd,k)
         write(6,*)'k,fm,dvmod,ri,csq ',
     &            k,fm(idjd),dvmod(idjd),ri(idjd,k),csq(idjd)
         write(6,*)'qfg,qfg+ ',qfg(idjd,k),qfg(idjd,k+1)
         write(6,*)'qlg,qlg+,dqtot ',qlg(idjd,k),qlg(idjd,k+1),dqtot(iq)
         write(6,*)'cfrac,cfrac+ ',cfrac(idjd,k),cfrac(idjd,k+1)
         write(6,*)'betatt,betatt+,betatt*delthet ',
     &            betatt(iq,k),betatt(iq,k+1),betatt(iq,k)*delthet(iq,k)
         write(6,*)'betaqt,betaqt+,betaqt*dqtot   ',
     &            betaqt(iq,k),betaqt(iq,k+1),betaqt(iq,k)*dqtot(iq)
         write(6,*)'x,zh,tmnht,dz,delh,sighkap ',x(idjd),zh(idjd,k),
     &                 tmnht(idjd,k),dz(idjd),delh(k),sighkap(k)
       endif     ! (diag.or.ntest>=1)
      enddo      ! end of k loop

      if( (diag.or.ntest>=1) .and. mydiag )then
        write(6,*)'before possible call to pbldif in vertmix'
        write (6,"('uav ',9f8.3/4x,9f8.3)") uav(idjd,:) 
        write (6,"('vav ',9f8.3/4x,9f8.3)") vav(idjd,:)
        write (6,"('t   ',9f8.3/4x,9f8.3)") t(idjd,:)
        write (6,"('qg ',3p9f8.3/4x,9f8.3)") qg(idjd,:)
        write (6,"('qs ',3p9f8.3/4x,9f8.3)") qs(idjd,:)
        do k=1,kl
         prcpv(k)=sig(k)**(-roncp)
        enddo       
        write (6,"('thee',9f8.3/4x,9f8.3)") 
     &        (prcpv(k)*t(idjd,k)*(t(idjd,k) + .5*hlcp*qs(idjd,k))
     &                   /(t(idjd,k) - .5*hlcp*qs(idjd,k)),k=1,kl)
      endif
      if(nmaxpr==1.and.mydiag)then
        write (6,"('rino_v',9f9.3/6x,9f9.3)") ri(idjd,1:kl-1)
        write (6,"('rkh0 ',9f9.3/5x,9f9.3)") rkh(idjd,1:kl-2)
        write (6,"('rkm0 ',9f9.3/5x,9f9.3)") rkm(idjd,1:kl-2)
      endif

      if(nlocal.ne.0)then
        call pbldif(rhs,rkh,rkm,uav,vav)  ! rhs is theta or thetal
!       n.b. *** pbldif partially updates qg and theta (t done during trim)	 
!       and updates rkh and rkm arrays
        if(nmaxpr==1.and.mydiag)then
          write (6,"('pblh ',f8.2)") pblh(idjd)
          write (6,"('rkh1 ',9f9.3/5x,9f9.3)") rkh(idjd,1:kl-2)
          write (6,"('rkm1 ',9f9.3/5x,9f9.3)") rkm(idjd,1:kl-2)
        endif
        if(nmaxpr==1.and.mydiag)
     &    write (6,"('thet_pbl',9f8.3/8x,9f8.3)") rhs(idjd,:)
        if( (diag.or.ntest>=1) .and. mydiag )
     &    write (6,"('qg ',3p9f8.3/4x,9f8.3)") qg(idjd,:)
        if(diag)then
          call printa('rkh ',rkh,ktau,nlv,ia,ib,ja,jb,0.,1.)
          call printa('cond',condx,ktau,1,ia,ib,ja,jb,0.,1.)
          call printa('zs  ',zs,ktau,1,ia,ib,ja,jb,0.,1.)
        endif
      endif      ! (nlocal>0)

      rk_shal(:,:)=0.
c     ***** ***** section for jlm shallow convection v4 *****************
      if(ksc==81)then
        do k=1,ksctop-1   ! or ksctop?  
         do iq=1,ifull
          if(sig(ktsav(iq))>sig_ct.and.k<ktsav(iq).  
     &                                   and.condc(iq)==0.)then  
            rk_shal(iq,k)=tied_con
            rk_shal(iq,k+1)=tied_over
          endif ! (sig(ktsav(iq))>sig_ct.and.k<ktsav(iq).and....)
         enddo  ! iq loop
        enddo   !  k loop
      endif     ! (ksc==81)
      if(ksc==82)then
        do k=1,ksctop-1    
         do iq=1,ifull
          if(sig(ktsav(iq))>sig_ct.and.k<ktsav(iq).
     &                and.k>=kbsav(iq).and.condc(iq)==0.)then  
            rk_shal(iq,k)=tied_con
            rk_shal(iq,k+1)=tied_over
          endif ! (sig(ktsav(iq))>sig_ct.and.k<ktsav(iq).and....)
         enddo  ! iq loop
        enddo   !  k loop
      endif     ! (ksc==82)
      if(ksc==83)then
        do k=1,ksctop-1    
         do iq=1,ifull
          if(sig(k)>sig_ct.and.k<ktsav(iq).
     &                and.k>=kbsav(iq).and.condc(iq)==0.)then  
            rk_shal(iq,k)=tied_con
            rk_shal(iq,k+1)=tied_over
          endif ! (sig(ktsav(iq))>sig_ct.and.k<ktsav(iq).and....)
         enddo  ! iq loop
        enddo   !  k loop
      endif     ! (ksc==83)
      if(ksc==91)then
        do k=1,ksctop ! May 08
         do iq=1,ifull
          if(ktsav(iq)<kl-1.and.k<ktsav(iq))then  ! April 04
            rk_shal(iq,k)=tied_con
            rk_shal(iq,k+1)=tied_over
          endif ! (ktsav(iq)<0.and.k<abs(ktsav(iq)))
         enddo  ! iq loop
        enddo   !  k loop
      endif     ! (ksc==91)
      if(ksc==92)then
        do k=1,ksctop ! May 08
         do iq=1,ifull
          if(k>=kbsav(iq).and.k<ktsav(iq).and.condc(iq)==0.)then  ! May 08
            rk_shal(iq,k)=tied_con
            rk_shal(iq,k+1)=tied_over
          endif ! (ktsav(iq)<0.and.k<abs(ktsav(iq)))
         enddo  ! iq loop
        enddo   !  k loop
      endif     ! (ksc==92)
c     *********** end of jlm shallow convection section *****************

c     ************ section for Tiedtke shallow convection *******************
      if(ksc==99)then
        do iq=1,ifull
        theeb(iq)=prcpv(kscbase)*t(iq,kscbase)*
     &                  (t(iq,kscbase) + .5*hlcp*qs(iq,kscbase))
     &                 /(t(iq,kscbase) - .5*hlcp*qs(iq,kscbase))
        enddo    ! iq loop
        if(kscsea==1)then  ! Tiedtke done only over sea
          do k=kscbase+1,ksctop
           do iq=1,ifull
            if(.not.land(iq))then 
             thee(iq,k)=prcpv(k)*t(iq,k)*(t(iq,k) + .5*hlcp*qs(iq,k))
     &                           /(t(iq,k) - .5*hlcp*qs(iq,k))
             if(qg(iq,kscbase)>tied_rh*qs(iq,kscbase).
     &        and.thee(iq,k)<theeb(iq))then         !  default tied_rh=.75
              rk_shal(iq,k-1)=tied_con     !  m**2/sec  6., originally 10.
              rk_shal(iq,k)=tied_over      !  m**2/sec
             endif ! (qg(iq,kscbase)>rhscon*qs(iq,kscbase).....
            endif  ! (.not.land(iq)) 
           enddo   ! iq loop
          enddo    ! end of k=kscbase+1,ksctop loop
        else       !  i.e. Tiedtke original scheme over land and sea
          do k=kscbase+1,ksctop  ! typically kscbase=3 & ksctop=6
           do iq=1,ifull
            thee(iq,k)=prcpv(k)*t(iq,k)*(t(iq,k) + .5*hlcp*qs(iq,k))
     &                           /(t(iq,k) - .5*hlcp*qs(iq,k))
            if(qg(iq,kscbase)>tied_rh*qs(iq,kscbase).
     &       and.thee(iq,k)<theeb(iq))then !  default tied_rh=.75
              rk_shal(iq,k-1)=tied_con     !  m**2/sec  6., originally 10.
              rk_shal(iq,k)=tied_over      !  m**2/sec
              if(ntest==3.and.k==ksctop)then
                write(6,*),'ktau,iq,theeb,thee,delthee ',
     &                 ktau,iq,theeb(iq),thee(iq,k),theeb(iq)-thee(iq,k)
              endif
            endif ! (qg(iq,kscbase)>rhscon*qs(iq,kscbase).....
           enddo  ! iq loop
          enddo   ! end of k=kscbase+1,ksctop loop
        endif     ! (kscsea==1)  .. else ..
      endif       ! (ksc==99)
c     *********** end of Tiedtke shallow convection section **************

c     *********** Tiedtke_ldr-style shallow convection 93 ************
      if(ksc>=93.and.ksc<=96)then   
c       Calculate LCL for near surface air
c       Assume qstar=qtg(iq,1) but allow some sub-grid variability
c       The formula for tsp is eqn (21) of Bolton (1980), MWR 108, 1046-1053.
        if(sigkscb>0.)then   ! uses qg1
          do iq=1,ifull
!          Leon used factor 1.01 over sea; here div by tied_rh (e.g. .99)	 
           epart=qg(iq,1)*.01*ps(iq)/(tied_rh*epsil) !in hPa 
           tsp=2840./(3.5*log(tscrn(iq))-log(epart)-4.805) + 55.   
           sigsp(iq)=(tsp/tscrn(iq))**(1./roncp)
          enddo
        else                    ! uses qgscrn
          do iq=1,ifull
!          Leon used factor 1.01 over sea; here div by tied_rh (e.g. .99)	 
           epart=qgscrn(iq)*.01*ps(iq)/(tied_rh*epsil) !in hPa   
           tsp=2840./(3.5*log(tscrn(iq))-log(epart)-4.805) + 55.   
           sigsp(iq)=(tsp/tscrn(iq))**(1./roncp)
          enddo
        endif  ! (sigkscb>0.) .. else ..
c       Look for the lowest layer s.t. the top of the layer is above the LCL,
c       and call that layer the cloud base. 
        kbase(:)=kl
        do k=ksctop-1,1,-1
         do iq=1,ifull
          if(sigsp(iq)>sigmh(k+1))kbase(iq)=k
         enddo
        enddo
        if(nlocal.ne.0)then  ! like old ksc=95, but ensures LCL within PBL
          do iq=1,ifull
           if(zh(iq,kbase(iq))>pblh(iq))then 
              kbase(iq)=kl
           endif
          enddo  ! iq loop
        endif  ! (nlocal.ne.0)
c       following has some jlm shortcuts	 
        do iq=1,ifull
         k=kbase(iq)
         qbas=qs(iq,k)
         theeb(iq)=prcpv(k)*t(iq,k)*(t(iq,k) + .5*hlcp*qbas)
     &                             /(t(iq,k) - .5*hlcp*qbas)
        enddo  ! iq loop
        ktop(:)=kbase(:)
        do k=2,ksctop+1
         do iq=1,ifull
          thee(iq,k)=prcpv(k)*t(iq,k)*(t(iq,k) + .5*hlcp*qs(iq,k))
     &                         /(t(iq,k) - .5*hlcp*qs(iq,k))
          if(theeb(iq)>thee(iq,k).and.ktop(iq)==k-1)then
            ktop(iq)=k       ! also checking if contiguous
          endif
         enddo  ! iq loop
        enddo   ! k loop	 
        if(ksc==94)then  ! from April same as 93* 
          do k=2,ksctop
           do iq=1,ifull
c            if(ktop(iq)>kbase(iq).and.k<=ktop(iq).
c     &                               and.k>kbase(iq))then
            if(k>kbase(iq).and.k<=ktop(iq))then
              rk_shal(iq,k-1)=tied_con     !  m**2/sec  6., originally 10.
              rk_shal(iq,k)=tied_over      !  m**2/sec
            endif
           enddo  ! iq loop
          enddo   ! k loop
        endif     ! (ksc==94)
        if(ksc==93)then
          do k=2,ksctop
           do iq=1,ifull
            if(ktop(iq)>kbase(iq).and.k<=ktop(iq).
c**  .            and.k>kbase(iq))then   !   this used for 93*
     &            and.k>kbase(iq).and.ktop(iq)<=ksctop)then
              rk_shal(iq,k-1)=tied_con     !  m**2/sec  6., originally 10.
              rk_shal(iq,k)=tied_over      !  m**2/sec
            endif
           enddo  ! iq loop
          enddo   ! k loop
        endif     ! (ksc==93)
        if(ksc==95)then   ! new from 7 April
          do k=2,ksctop
           do iq=1,ifull
            if(k>kbase(iq).and.k<=ktop(iq).and.condc(iq)==0.)then
              rk_shal(iq,k-1)=tied_con     !  m**2/sec  6., originally 10.
              rk_shal(iq,k)=tied_over      !  m**2/sec
            endif
           enddo  ! iq loop
          enddo   ! k loop
        endif     ! (ksc==95)
        if(ksc==96)then   ! applied from sfce up
          do k=2,ksctop
           do iq=1,ifull
c**         if(ktop(iq)>kbase(iq).and.k<=ktop(iq).
c**  .                               and.ktop(iq)<=ksctop)then
            if(ktop(iq)>kbase(iq).and.k<=ktop(iq).  
     &                                   and.condc(iq)==0.)then  
              rk_shal(iq,k-1)=tied_con     !  m**2/sec  6., originally 10.
              rk_shal(iq,k)=tied_over      !  m**2/sec
            endif
           enddo  ! iq loop
          enddo   ! k loop
        endif     ! (ksc==96)
      endif       ! (ksc>=93.and.ksc<=96)
c     *********** end of Tiedtke_ldr shallow convection 93-96 *********

c     *************** Tiedtke_jlm shallow convection 97 ***************
      if(ksc==97)then
!       this one is similar to ksc=99, but uses enhanced qbas <= qs
!       and also works up the column with other possible thebas values
        do k=kscbase,ksctop
         do iq=1,ifull
          qbas=min(qg(iq,k)/tied_rh,qs(iq,k))
c         thebas(iq,k)=t(iq,k)+hlcp*qbas  + hght
          thebas(iq,k)=prcpv(k)*t(iq,k)*(t(iq,k) + .5*hlcp*qbas)
     &                                 /(t(iq,k) - .5*hlcp*qbas)
         enddo  ! iq loop
        enddo   ! k loop
        theeb(:)=thebas(:,kscbase)
        do k=kscbase+1,ksctop+1
         do iq=1,ifull
          thee(iq,k)=prcpv(k)*t(iq,k)*(t(iq,k) + .5*hlcp*qs(iq,k))
     &                         /(t(iq,k) - .5*hlcp*qs(iq,k))
          if(theeb(iq)>thee(iq,k))then
            rk_shal(iq,k-1)=tied_con       !  m**2/sec  6., originally 10.
            rk_shal(iq,k)=tied_over        !  m**2/sec
          else
            theeb(iq)=thebas(iq,k) ! ready for next k in k-loop
          endif
         enddo  ! iq loop
        enddo   ! k loop
       if(ntest==3)then
         do iq=1,ifull
          if(rk_shal(iq,ksctop-1)>.9*tied_con)then 
            write(6,*) 'iq,land,rk_shal ',iq,land(iq),
     &                  rk_shal(iq,ksctop-1)
          endif
          enddo
        endif   ! (ntest==3)
      endif     ! (ksc==97)
c     *********** end of Tiedtke_jlm shallow convection 97 *************

c     add in effects of shallow convection
      do k=1,kl-1  
       rkh(:,k)=rkh(:,k)+rk_shal(:,k)
      enddo   !  k loop
      if(kscmom==1)then
        do k=1,kl-1  
         rkm(:,k)=rkm(:,k)+rk_shal(:,k)
        enddo   !  k loop
      endif     ! (kscmom==1)
      
      if(ksc.ne.0.and.(ntest.ne.0.or.diag).and.nproc==1.)then
        do iq=1,ifull
         if(rk_shal(iq,1)>rk_shal(iq,2))then
           write(6,*) 'iq,rk_shal1,rk_shal2',iq,rk_shal(iq,1),
     &                 rk_shal(iq,2)
         endif
        enddo
        if (mydiag) then 
        iq=idjd
        write(6,*)'for shallow convection in vertmix '
        write(6,*)'ktsav,ksctop ',ktsav(idjd),ksctop
        write(6,*)'kbase,ktop ',kbase(idjd),ktop(idjd)
        write(6,*)'kbsav,ktsav,theeb: ',kbsav(iq),ktsav(iq),theeb(iq)
        write (6,"('rk_shal ',9f7.2/(9x,9f7.2))")
     &             (rk_shal(idjd,k),k=1,kl)
        write (6,"('rh   ',9f7.2/(5x,9f7.2))") 
     &             (100.*qg(idjd,k)/qs(idjd,k),k=1,kl)
        write (6,"('qs   ',9f7.3/(5x,9f7.3))") 
     &             (1000.*qs(idjd,k),k=1,kl)
        write (6,"('qg   ',9f7.3/(5x,9f7.3))") 
     &             (1000.*qg(idjd,k),k=1,kl)
        write (6,"('qbas ',9f7.3/(5x,9f7.3))") 
     &             (1000.*qg(idjd,k)/tied_rh,k=1,kl)
        write (6,"('t    ',9f7.2/(5x,9f7.2))") 
     &             (t(idjd,k),k=1,kl)
        write (6,"('thebas',9f7.2/(6x,9f7.2))") 
     &             (thebas(iq,k),k=1,kl)
c        write (6,"('hs  ',9f7.2/(4x,9f7.2))") 
c     &             (t(idjd,k)+hlcp*qs(idjd,k),k=1,kl)
        write (6,"('thee',9f7.2/(4x,9f7.2))") 
     &             (thee(idjd,k),k=1,kl)
        endif ! mydiag
      endif
      
      else
      
       !-------------------------------------------------------------
       ! TKE-eps closure scheme
      
       ! note ksc.ne.0 options are clobbered when nvmix=6
       ! However, nvmix=6 with nlocal=7 supports its own shallow
       ! convection options
       zg(:,1)=bet(1)*t(1:ifull,1)/grav
       do k=2,kl
        zg(:,k)=zg(:,k-1)+(bet(k)*t(1:ifull,k)
     &                    +betm(k)*t(1:ifull,k-1))/grav
       end do ! k  loop
       zg=zg+phi_nh/grav ! add non-hydrostatic component
       do k=1,kl-1
         zgh(:,k)=ratha(k)*zg(:,k+1)+rathb(k)*zg(:,k)
       end do
       wq0=eg/(hl*rho)
       wt0=fg/(cp*rho)
       select case(nlocal)
        case(0) ! no counter gradient
         call tkemix(rkm,rhs,qg(1:ifull,:),qlg(1:ifull,:),
     &             qfg(1:ifull,:),cfrac,pblh,wt0,wq0,
     &             ps(1:ifull),ustar,zg,zgh,sig,sigkap,dt,
     &             qgmin,1,0)
         rkh=rkm
        case(1,2,3,4,5,6) ! KCN counter gradient method
         call tkemix(rkm,rhs,qg(1:ifull,:),qlg(1:ifull,:),
     &             qfg(1:ifull,:),cfrac,pblh,wt0,wq0,
     &             ps(1:ifull),ustar,zg,zgh,sig,sigkap,dt,
     &             qgmin,1,0)
         rkh=rkm
         uav(1:ifull,:)=av_vmod*u(1:ifull,:)
     &                 +(1.-av_vmod)*savu(1:ifull,:)
         vav(1:ifull,:)=av_vmod*v(1:ifull,:)
     &                 +(1.-av_vmod)*savv(1:ifull,:)
         call pbldif(rhs,rkh,rkm,uav,vav)
        case(7) ! mass-flux counter gradient
         call tkemix(rkm,rhs,qg(1:ifull,:),qlg(1:ifull,:),
     &             qfg(1:ifull,:),cfrac,pblh,wt0,wq0,
     &             ps(1:ifull),ustar,zg,zgh,sig,sigkap,dt,
     &             qgmin,0,0)
         rkh=rkm
        case DEFAULT
          write(6,*) "ERROR: Unknown nlocal option for nvmix=6"
          stop
       end select 
      end if ! nvmix.ne.6


      !--------------------------------------------------------------
      ! Perform mixing on prognstic arrays

      do k=1,kl-1
        delsig=sig(k+1)-sig(k)
        do iq=1,ifull
          dz(iq) =-tmnht(iq,k)*delons(k)  ! this is z(k+1)-z(k)
          dzr(iq)=1./dz(iq)
          guv(iq,k)= rkm(iq,k)*dt *delsig *dzr(iq)**2
          gt(iq,k)=  rkh(iq,k)*dt *delsig *dzr(iq)**2
         enddo   ! iq loop
      enddo      ! k loop

      guv(:,kl)=0.
      gt(:,kl)=0.
      if(diag)then
        call maxmin(rkh,'rk',ktau,.01,kl-1)
        if ( mydiag ) then
           write(6,*)'vertmix guv ',(guv(idjd,k),k=1,kl)
           write(6,*)'vertmix gt ',(gt(idjd,k),k=1,kl)
        end if
        call printa('t   ',t,ktau,nlv,ia,ib,ja,jb,200.,1.)
        call printa('tss ',tss,ktau,nlv,ia,ib,ja,jb,200.,1.)
        call printa('eg  ',eg,ktau,nlv,ia,ib,ja,jb,0.,1.)
        call printa('fg  ',fg,ktau,nlv,ia,ib,ja,jb,0.,1.)
        call printa('thet',rhs,ktau,nlv,ia,ib,ja,jb,200.,1.)
      endif

      if(ncvmix>0)then  ! cumulus mixing of momentum - jlm version
        do k=kuocb+1-ncvmix,kcl_top-1
         do iq=1,ifull
!          for no-conv-points, doing k=1-ncvmix,1 with convpsav=0.
!          1 below for ncvmix=2
           if(kbsav(iq)>0.and.k>=kbsav(iq)+1-ncvmix
     &                      .and.k<ktsav(iq))
     &       guv(iq,k)=guv(iq,k)-convpsav(iq)*.5      ! with factor of .5
          enddo  ! iq loop
        enddo    ! k loop
        if(diag.and.mydiag)then
          write(6,*)'vertmix after conv; kb,kt,convp'
     &     ,kbsav(idjd),ktsav(idjd),convpsav(idjd)
          write(6,*)'new guv',(guv(idjd,k),k=1,kl)
        endif
      endif      !   (ncvmix>0)

      conflux=grav*dt/dsig(1)
      condrag=grav*dt/(dsig(1)*rdry)
c     first do theta (then convert back to t)
      at(:,1)=0.

      do k=2,kl
       do iq=1,ifull
        at(iq,k) =-gt(iq,k-1)/dsig(k)
       enddo   ! iq loop
      enddo    !  k loop
      do k=1,kl
       do iq=1,ifull
        ct(iq,k) =-gt(iq,k)/dsig(k)
       enddo   ! iq loop
      enddo    !  k loop
      if((diag.or.ntest==2).and.mydiag)then
        write(6,*)'ktau,fg,tss,ps ',ktau,fg(idjd),tss(idjd),ps(idjd)
        write(6,*)'at ',(at(idjd,k),k=1,kl)
        write(6,*)'ct ',(ct(idjd,k),k=1,kl)
        write(6,*)'rhs ',(rhs(idjd,k),k=1,kl)
      endif      ! (ntest==2)

      !--------------------------------------------------------------
      ! Temperature
      if (nvmix.ne.6) then
       if(nmaxpr==1.and.mydiag)
     &   write (6,"('thet_inx',9f8.3/8x,9f8.3)") rhs(idjd,:)
       rhs(1:ifull,1)=rhs(1:ifull,1)-(conflux/cp)*fg(1:ifull)
     &               /ps(1:ifull)
       call trim(at,ct,rhs,0)   ! for t
       if(nmaxpr==1.and.mydiag)
     &   write (6,"('thet_out',9f8.3/8x,9f8.3)") rhs(idjd,:)
      end if
      do k=1,kl
        do iq=1,ifull
         t(iq,k)=rhs(iq,k)/sigkap(k)
        enddo  ! iq loop
      enddo    !  k loop
      if(diag)then
         if (mydiag) then
            write(6,*)'vertmix eg,fg,cdtq,land '
     &                  ,eg(idjd),fg(idjd),cdtq(idjd),land(idjd)
         end if
        call printa('thet',rhs,ktau,nlv,ia,ib,ja,jb,200.,1.)
      endif

      !--------------------------------------------------------------
      ! Moisture
      if (nvmix.ne.6) then
       rhs(1:ifull,:)=qg(1:ifull,:)
       rhs(1:ifull,1)=rhs(1:ifull,1)-(conflux/hl)*eg(1:ifull)
     &               /ps(1:ifull)
c      could add extra sfce moisture flux term for crank-nicholson
       call trim(at,ct,rhs,0)    ! for qg
       qg(1:ifull,:)=rhs(1:ifull,:)
      else if (nlocal.gt.0) then
       ! increase mixing to replace counter gradient term
       at=cq*at
       ct=cq*ct
      end if ! ..else.. (nvmix.ne.6)
      if(diag.and.mydiag)then
       write(6,*)'vertmix rhs & qg after trim ',(rhs(idjd,k),k=1,kl)
       write (6,"('qg ',9f7.3/(8x,9f7.3))")
     &            (1000.*qg(idjd,k),k=1,kl)
      endif

      !--------------------------------------------------------------
      ! Cloud microphysics terms
      if(ldr.ne.0)then
c       now do qfg
        rhs(1:ifull,:)=qfg(1:ifull,:)
        call trim(at,ct,rhs,0)    ! for qfg
        qfg(1:ifull,:)=rhs(1:ifull,:)
c       now do qlg
        rhs(1:ifull,:)=qlg(1:ifull,:)
        call trim(at,ct,rhs,0)    ! for qlg
        qlg(1:ifull,:)=rhs(1:ifull,:)
c       now do cfrac
        rhs(1:ifull,:)=cfrac(1:ifull,:)
        call trim(at,ct,rhs,0)    ! for cfrac
        cfrac(1:ifull,:)=min(max(rhs(1:ifull,:),0.),1.)
c       now do qrg
        rhs(1:ifull,:)=qrg(1:ifull,:)
        call trim(at,ct,rhs,0)    ! for qrg
        qrg(1:ifull,:)=rhs(1:ifull,:)
c       now do cffall
        rhs(1:ifull,:)=cffall(1:ifull,:)
        call trim(at,ct,rhs,0)    ! for cffall
        cffall(1:ifull,:)=min(max(rhs(1:ifull,:),0.),1.)
      endif    ! (ldr.ne.0)

      !--------------------------------------------------------------
      ! Aerosols
      if (abs(iaero)==2) then
        do l=1,naero
          rhs(1:ifull,:)=xtg(1:ifull,:,l)
          call trim(at,ct,rhs,0)
          xtg(1:ifull,:,l)=rhs(1:ifull,:)
        end do
      end if ! (abs(iaero)==2)

      !--------------------------------------------------------------
      ! Tracers
      if (ngas>0) then
        call tracervmix(at,ct)
      end if ! (ngas>0)

      !--------------------------------------------------------------
      ! Momentum terms
      au(:,1)=cduv(:)*condrag/tss(:)
      do k=2,kl
       do iq=1,ifull
        au(iq,k) =-guv(iq,k-1)/dsig(k)
       enddo   ! iq loop
      enddo    !  k loop
      do k=1,kl
       do iq=1,ifull
        cu(iq,k) =-guv(iq,k)/dsig(k)
       enddo   ! iq loop
      enddo    !  k loop
      if((diag.or.ntest==2).and.mydiag)then
        write(6,*)'au ',(au(idjd,k),k=1,kl)
        write(6,*)'cu ',(cu(idjd,k),k=1,kl)
      endif      ! (ntest==2)

c     first do u
      rhs(1:ifull,:)=u(1:ifull,:)
      call trim(au,cu,rhs,0)
      u(1:ifull,:)=rhs(1:ifull,:)
      if(diag.and.mydiag)then
        write(6,*)'vertmix au ',(au(idjd,k),k=1,kl)
      endif

c     now do v; with properly unstaggered au,cu
      rhs(1:ifull,:)=v(1:ifull,:)
      call trim(au,cu,rhs,0)    ! note now that au, cu unstaggered globpea
      v(1:ifull,:)=rhs(1:ifull,:)
 
      if((diag.or.ntest>=1).and.mydiag)then
        write(6,*)'after trim in vertmix '
        write (6,"('thet',9f7.2/(8x,9f7.2))") 
     &              (sigkap(k)*t(idjd,k),k=1,kl) 
        write (6,"('t   ',9f7.2/(8x,9f7.2))") (t(idjd,k),k=1,kl) 
        write (6,"('u   ',9f7.2/(8x,9f7.2))") (u(idjd,k),k=1,kl) 
        write (6,"('v   ',9f7.2/(8x,9f7.2))") (v(idjd,k),k=1,kl) 
        write(6,*)'cduv,cduv+1,cduvj+1,tss ',
     &           cduv(idjd),cduv(idjd+1),cduv(idjd+il),tss(idjd)
        write (6,"('au ',9f7.3/(8x,9f7.3))") (au(idjd,k),k=1,kl) 
        write (6,"('cu ',9f7.3/(8x,9f7.3))") (cu(idjd,k),k=1,kl) 
      endif
      if(nmaxpr==1.and.mydiag)then
        write (6,"('qg_vm ',9f7.3/(6x,9f7.3))") 
     &             (1000.*qg(idjd,k),k=1,kl)
      endif
      
      return
      end
