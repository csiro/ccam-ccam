      subroutine adjust5(iaero)
      use aerosolldr
      use arrays_m
      use cc_mpi
      use cfrac_m
      use diag_m
      use dpsdt_m
      use epst_m
      use indices_m
      use liqwpar_m  ! qfg,qlg
      use map_m
      use mgsolve
      use morepbl_m  ! condx,eg
      use nharrs_m
      use nlin_m
      use pbl_m
      use sigs_m
      use tbar2d_m
!     rml 19/09/07 replace gasmin from tracers.h with tracmin from tracermodule
      use tracermodule, only: tracmin
      use tracers_m
      use vecsuv_m
      use vecs_m
      use vvel_m     ! sdot
      use work3sav_m
      use xarrs_m
      use xyzinfo_m
      implicit none
      integer, parameter :: mfix_rad=0 ! used to make gases 2 to ng add up to gas 1
      integer, parameter :: ntest=0
      integer, parameter :: nys=0      ! for nys=1 use sun & yeh 1997 monotonic filter
!                              within adjust5 on u,v
!     nuv in parm.h:  now always 10
!          !  0y for increments   , 1y for actual values, 2y for true increments
      include 'newmpar.h'
      include 'const_phys.h'
      include 'kuocom.h'
      include 'parm.h'     ! qgmin
      include 'parmdyn.h'  
      include 'parmvert.h'  
      real, dimension(:), allocatable, save :: zz,zzn,zze,zzw,zzs
      real, dimension(:), allocatable, save :: pfact,alff,alf,alfe
      real, dimension(:), allocatable, save :: alfn,alfu,alfv
      real, dimension(ifull+iextra,kl) :: p,cc,dd,pe
      real, dimension(ifull,kl) :: omgfnl,wrk1,wrk2,wrk3,d,e
      real, dimension(ifull,kl) :: helm,rhsl,omgf
      real, dimension(ifull,kl) :: dumu,dumv,dumc,dumd,dumt
      real, dimension(ifull) :: ps_sav,pslxint,pslsav
      real, dimension(ifull) :: sdmx,delps,bb
      real, dimension(3) :: delarr, delarr_l
      real hdt, hdtds, sumx, qgminm, ratio, sumdiffb, alph_g
      real alph_p, alph_pm, delneg, delpos, delnegk, delposk, alph_q
      real sumdiffb_l  ! Local versions
      real sumin, sumout, sumsav, dum
      real delpos_l, delneg_l, const_nh
      real, save :: dtsave = 0.
      integer, dimension(ifull) :: nits, nvadh_pass
      integer its, k, l, iq, ng, ierr, iaero
      integer, save :: precon_in = -99999

      call start_log(adjust_begin)

      hdt=dt/2.
      hdtds=hdt/ds

      if (.not.allocated(zz)) then
        allocate(zz(ifull),zzn(ifull),zze(ifull),zzw(ifull))
        allocate(zzs(ifull),pfact(ifull),alff(ifull+iextra))
        allocate(alf(ifull+iextra),alfe(ifull+iextra))
        allocate(alfn(ifull+iextra),alfu(ifull),alfv(ifull))
      end if

      ! time step can change during initialisation
      if(dt /= dtsave) then 
         call adjust_init(zz,zzn,zze,zzw,zzs,pfact,alff,alf,alfe,alfn,
     &                    alfu,alfv)
         precon_in=precon
         precon=min(precon,0)  ! 22/4/07
         
         if (precon<-9999) then
           call mgsor_init
           call mgzz_init(zz,zzn,zze,zzw,zzs)
         end if
      end if

      if (diag.or.nmaxpr==1)then
         if (mydiag ) then
           write(6,*) 'entering adjust5'
           write (6,"('tx_a1',10f8.2)") tx(idjd,:)
           if(m<8)then
             write (6,"('ux_stag',10f8.2)") ux(idjd,:)
             write (6,"('vx_stag',10f8.2)") vx(idjd,:)
           else
             write (6,"('ux_stag',-5p10f8.2)") ux(idjd,:)
             write (6,"('vx_stag',-5p10f8.2)") vx(idjd,:)
           endif
           write (6,"('qg_a1',3p10f7.3)") qg(idjd,:)
           write (6,"('pslx_3p',3p9f8.4)") pslx(idjd,:)
         end if
         call printa('pslx',pslx,ktau,nlv,ia,ib,ja,jb,0.,1000.)
         call printa('tx  ',tx,ktau,nlv,ia,ib,ja,jb,200.,1.)
         call maxmin(alf,'a ',ktau,1.,1)
         call maxmin(alfe,'ae',ktau,1.,1)
         call maxmin(alfn,'an',ktau,1.,1)
      endif

!     recompute nonlinear sigma-dot contribution for updating tn, tx
!     e contains intgrl{-pslx x dsigma} from 0 to sigma (dsig is -ve)
!     vert. integ. nonlin part of mass weighted div into e
      do k=1,kl
       pslx(1:ifull,k)=pslx(1:ifull,k)*2./dt  ! i.e. [RHS of Eq. 115]*2/dt, i.e. M
      enddo    ! k  loop
      e(:,kl)=dsig(kl)*pslx(1:ifull,kl)
      wrk3(:,kl)=-sig(kl)*pslx(1:ifull,kl) ! integration following eig.f
      do k=kl-1,1,-1
       e(:,k)=e(:,k+1)+dsig(k)*pslx(1:ifull,k)  ! i.e. -M_bar_(sig-.5)
       wrk3(:,k)=e(:,k+1)+(sigmh(k+1)-sig(k))*pslx(1:ifull,k) ! integration following eig.f
      enddo     ! k loop
      pslxint(:)=-e(:,1)*dt/2. ! pslxint holds integrated pslx, Eq. 116

!     full-level (1.+epsp)*omega/ps into omgfnl (nonlin. part only), Eq. 118
      do k=1,kl
       omgfnl(1:ifull,k)=-wrk3(:,k)-sig(k)*pslx(1:ifull,k)
      enddo     ! k loop
     
!     redefine ux, vx
      do k=1,kl
        do iq=1,ifull  ! N.B. alfu is alf/(1+epsu)
          cc(iq,k)=ux(iq,k)/emu(iq)*alfu(iq)  ! Eq. 136
          dd(iq,k)=vx(iq,k)/emv(iq)*alfv(iq)  ! Eq. 137
        enddo
      enddo
      
#ifdef loadbalall
      call start_log(adjusta_loadbal_begin)
      call phys_loadbal
      call end_log(adjusta_loadbal_end)
#endif
      call boundsuv(cc,dd,stag=-9) ! only update isv and iwu
#ifdef loadbalall
      call start_log(adjustb_loadbal_begin)
      call phys_loadbal
      call end_log(adjustb_loadbal_end)
#endif

      do k=1,kl
       do iq=1,ifull
!       N.B. the omgfnl term on LHS of Eq. 121 not yet added in
        tx(iq,k)=tx(iq,k)+hdt*                ! adding LHS term of Eq. 121
     &           tbar2d(iq)*omgfnl(iq,k)*roncp/sig(k)  ! with correct epsp
       enddo    ! iq loop
      enddo     ! k loop

!     calculate hydrostatic heights from the tx array
      do iq=1,ifull
c      p(iq,1)=zs(iq)+bet(1)*tx(iq,1)+rdry*tbar2d(iq)*pslxint(iq) ! Eq. 146
       p(iq,1)=zs(iq)+bet(1)*(tx(iq,1)-280.)
     &         +rdry*tbar2d(iq)*pslxint(iq) ! Eq. 146
      enddo     ! iq loop
      do k=2,kl
       do iq=1,ifull
        p(iq,k)=p(iq,k-1)+bet(k)*(tx(iq,k)-280.)
     &                  +betm(k)*(tx(iq,k-1)-280.)
       enddo    ! iq loop
      enddo     ! k loop

      if(nh/=0.and.(ktau>knh.or.lrestart))then
!       add in departure values of p-related nh terms  & omgfnl terms    
        if (abs(epsp)<=1.) then
          const_nh=2.*rdry/(dt*grav*grav*(1.+abs(epsp))**2)
        else
          const_nh=2.*rdry/(dt*grav*grav)
        end if
        do k=1,kl
         ! MJT suggestion
         ! omgfnl already includes (1+epst)
         wrk2(:,k)=const_nh*tbar2d(:)*
     &     (tbar(1)*omgfnl(:,k)/sig(k)-h_nh(1:ifull,k))
        enddo
        wrk1(:,1)=bet(1)*wrk2(:,1)
        do k=2,kl
         wrk1(:,k)=wrk1(:,k-1)+bet(k)*wrk2(:,k)+betm(k)*wrk2(:,k-1)
        enddo   ! k loop
        if ((diag.or.nmaxpr==1).and.mydiag)then
          write(6,*) 'adjust5 omgfnl ',(omgfnl(idjd,k),k=1,kl)
          write(6,*) 'adjust5 h_nh ',(h_nh(idjd,k),k=1,kl)
          write(6,*) 'adjust5 pa ',(p(idjd,k),k=1,kl)
          write(6,*) 'adjust5 wrk1 ',(wrk1(idjd,k),k=1,kl)
        endif
        p(1:ifull,:)=p(1:ifull,:)+wrk1(:,:)  ! nh
      endif     ! (nh.ne.0.and.(ktau.gt.knh.or.lrestart))

!     form divergence of rhs (xu & xv) terms
      do k=1,kl
!cdir nodep
         do iq=1,ifull
!           d is xd in Eq. 157, divided by em**2/ds
            d(iq,k)=cc(iq,k)-cc(iwu(iq),k)+dd(iq,k)-dd(isv(iq),k)
         enddo
      enddo

!     transform p & d to eigenvector space
      do k=1,kl
         do iq=1,ifull
            pe(iq,k)=einv(k,1)*p(iq,1)
            rhsl(iq,k)=einv(k,1)*d(iq,1)
         enddo
         do l=2,kl
            do iq=1,ifull
               pe(iq,k)=pe(iq,k)+einv(k,l)*p(iq,l)     ! xp in eig space
               rhsl(iq,k)=rhsl(iq,k)+einv(k,l)*d(iq,l) ! xd in eig space
            enddo
         enddo                  ! l loop

         do iq=1,ifull
!           N.B.   pfact(iq)=4.*( ds/(dt*em(iq)) )**2    
            helm(iq,k) = pfact(iq)*tbar(1)/
     &                   (bam(k)*(1.+epst(iq))*tbar2d(iq))
            rhsl(iq,k) = rhsl(iq,k)/hdtds -helm(iq,k)*pe(iq,k) ! Eq. 161 mult
         enddo                  ! iq loop

      enddo    ! k loop

#ifdef loadbalall
      call start_log(adjustc_loadbal_begin)
      call phys_loadbal
      call end_log(adjustc_loadbal_end)
#endif

      if (precon<-9999) then
        ! Multi-grid
        call mghelm(zz,zzn,zze,zzw,zzs,helm,pe,rhsl)
      else if (precon<0) then
        ! SOR
        call helmsor(zz,zzn,zze,zzw,zzs,helm,pe,rhsl)
      else
        ! Congujate gradient
        call bounds(pe)
        call helmsol(zz,zzn,zze,zzw,zzs,helm,pe,rhsl)
      endif ! (precon<-9999) .. else ..

#ifdef loadbalall
      call start_log(adjustd_loadbal_begin)
      call phys_loadbal
      call end_log(adjustd_loadbal_end)
#endif

      if (diag.or.nmaxpr==1)then   !  only for last k of loop (i.e. 1)
         ! Some diagnostics requiring extra bounds calls have been removed
         if (mydiag ) then
            write (6,"('tx_a2',10f8.2)") tx(idjd,:)
            write(6,*) 'omgfnl_a2 ',omgfnl(idjd,:)
            write(6,*) 'adjust5 cc ',cc(idjd,:)
            write(6,*) 'adjust5 dd ',dd(idjd,:)
            write(6,*) 'adjust5 d_in ',d(idjd,:)
            write(6,*) 'adjust5 p_in ',p(idjd,:)
            write(6,*) 'adjust5 pe ',pe(idjd,:)
            write(6,*) 'adjust5 pe_e ',pe(ie(idjd),:)
            write(6,*) 'adjust5 pe_w ',pe(iw(idjd),:)
            write(6,*) 'adjust5 rhsl ',rhsl(idjd,:)
            write(6,*) 'adjust5 helm ',helm(idjd,:)
         end if
         call printa('psnt',pslxint,ktau,0,ia,ib,ja,jb,0.,1000.)
         call printa('tx  ',tx,ktau,nlv,ia,ib,ja,jb,200.,1.)
         call printa('rhsl',rhsl,ktau,1,ia,ib,ja,jb,0.,0.)
         call printa('pe  ',pe,ktau,1,ia,ib,ja,jb,0.,0.)
      endif

      do k=1,kl
!        first p from pe
         do iq=1,ifull
            p(iq,k)=emat(k,1)*pe(iq,1)
         enddo                  ! iq loop
         do l=2,kl              ! this is the remaining expensive loop
            do iq=1,ifull
               p(iq,k)=p(iq,k)+emat(k,l)*pe(iq,l)
            enddo               ! iq loop
         enddo                  !  l loop
      end do

#ifdef loadbalall
      call start_log(adjuste_loadbal_begin)
      call phys_loadbal
      call end_log(adjuste_loadbal_end)
#endif
      call bounds(p,corner=.true.)
#ifdef loadbalall
      call start_log(adjustf_loadbal_begin)
      call phys_loadbal
      call end_log(adjustf_loadbal_end)
#endif

!      now u & v
      if ((diag.or.nmaxpr==1).and.mydiag)then
       write(6,*) 'iq,k,fu,alfu,alfu*ux(iq,k) ',
     &         idjd,nlv,fu(idjd),alfu(idjd),alfu(idjd)*ux(idjd,nlv)
       write(6,*) 'alfF & n e w s (in(iq)),alfF(ine(iq)),alfF(is(iq))',
     &        'alfF(ise(iq)),alfe(iq) ',alfF(in(idjd)),alfF(ine(idjd))
     &         ,alfF(is(idjd)),alfF(ise(idjd)),alfe(idjd)
       sumx = alf(ie(idjd))-alf(idjd)+.25*(alfF(in(idjd))+
     &      alfF(ine(idjd))-alfF(is(idjd))-alfF(ise(idjd)))-alfe(idjd)
       write(6,*) 'sum  ',sumx
       write(6,*) 'p & n e w s ne se ',p(idjd,nlv),p(in(idjd),nlv),
     &         p(ie(idjd),nlv),p(iw(idjd),nlv),p(is(idjd),nlv),
     &         p(ine(idjd),nlv),p(ise(idjd),nlv)
       write(6,*) 'p & direct n s ',
     &         p(idjd,nlv),p(idjd+il,nlv),p(idjd-il,nlv)
      endif

      do k=1,kl
!cdir nodep
         do iq=1,ifull
            cc(iq,k) = alfu(iq)*ux(iq,k) ! Eq. 139
     &       -hdtds*emu(iq)*
     &       ( alf(ie(iq))*p(ie(iq),k)-alf(iq)*p(iq,k)
     &       -.5*alfe(iq)*(p(iq,k)+p(ie(iq),k))
     &       +.25*(alfF(in(iq))*p(in(iq),k) +alfF(ine(iq))*p(ine(iq),k)
     &       -alfF(is(iq))*p(is(iq),k) -alfF(ise(iq))*p(ise(iq),k)) )
            dd(iq,k) = alfv(iq)*vx(iq,k) ! Eq. 140
     &       -hdtds*emv(iq)*
     &       ( alf(in(iq))*p(in(iq),k)-alf(iq)*p(iq,k)
     &       -.5*alfn(iq)*(p(iq,k)+p(in(iq),k))
     &       -.25*(alfF(ien(iq))*p(ien(iq),k) +alfF(ie(iq))*p(ie(iq),k)
     &       -alfF(iwn(iq))*p(iwn(iq),k) -alfF(iw(iq))*p(iw(iq),k)) )
         enddo  ! iq loop   
      enddo     !  k loop 

#ifdef loadbalall
      call start_log(adjustg_loadbal_begin)
      call phys_loadbal
      call end_log(adjustg_loadbal_end)
#endif
      call boundsuv(cc,dd,stag=-9) ! only update isv and iwu
#ifdef loadbalall
      call start_log(adjusth_loadbal_begin)
      call phys_loadbal
      call end_log(adjusth_loadbal_end)
#endif

!     calculate linear part only of sigma-dot and omega/ps
      do k=1,kl
!cdir nodep
         do iq=1,ifull
            d(iq,k)=(cc(iq,k)/emu(iq)-cc(iwu(iq),k)/emu(iwu(iq)) ! Eq. 101
     &          +dd(iq,k)/emv(iq)-dd(isv(iq),k)/emv(isv(iq)))    
     &          *em(iq)**2/ds
         enddo                  ! iq loop
      enddo     ! k  loop
      if(nmaxpr==1)then
        call maxmin(d,'dv',ktau,0.,kl)
        if(nproc==1)write(6,*) 'cc,cc-,dd,dd-',
     &     cc(idjd,nlv)/emu(idjd),cc(iwu(idjd),nlv)/emu(iwu(idjd)),
     &     dd(idjd,nlv)/emv(idjd),dd(isv(idjd),nlv)/emv(isv(idjd))      
      endif   ! (nmaxpr==1)

!     npex=4 add un, vn on staggered grid in belated split manner
      if(npex==4)then  
        cc(1:ifull,:)=cc(1:ifull,:)+.5*dt*un(1:ifull,:)
        dd(1:ifull,:)=dd(1:ifull,:)+.5*dt*vn(1:ifull,:)
      endif

#ifdef loadbalall
      call start_log(adjusti_loadbal_begin)
      call phys_loadbal
      call end_log(adjusti_loadbal_end)
#endif


!     straightforward rev. cubic interp of u and v (i.e. nuv=10)
!     This is necessary because staguv expects arrays dimensioned 
!     (ifull,kl). 
      if(nstag==0)then
        call staguv(u(1:ifull,:),v(1:ifull,:),        
     &              wrk1(1:ifull,:),wrk2(1:ifull,:)) 
        if(nmaxpr==1.and.nproc==1)then
          its=ifull+iextra
          write (6,"('u_u0 ',10f8.2)") (u(iq,nlv),iq=idjd-3,idjd+3)
          write (6,"('v_u0 ',10f8.2)") 
     &               (v(iq,nlv),iq=idjd-3*its,idjd+3*its,its)
          write (6,"('u_s0 ',10f8.2)") (wrk1(iq,nlv),iq=idjd-3,idjd+3)
          write (6,"('v_s0 ',10f8.2)") 
     &               (wrk2(iq,nlv),iq=idjd-3*ifull,idjd+3*ifull,ifull)
          write (6,"('u_s1 ',10f8.2)") (cc(iq,nlv),iq=idjd-3,idjd+3)
          write (6,"('v_s1 ',10f8.2)") 
     &               (dd(iq,nlv),iq=idjd-3*its,idjd+3*its,its)       
        endif
        wrk1(1:ifull,:)=cc(1:ifull,:)-wrk1(1:ifull,:) ! staggered increment
        wrk2(1:ifull,:)=dd(1:ifull,:)-wrk2(1:ifull,:) ! staggered increment
        call unstaguv(wrk1(1:ifull,:),wrk2(1:ifull,:),        
     &                wrk1(1:ifull,:),wrk2(1:ifull,:)) 
        u(1:ifull,:)=u(1:ifull,:)+wrk1(1:ifull,:)
        v(1:ifull,:)=v(1:ifull,:)+wrk2(1:ifull,:)
        if(nmaxpr==1.and.nproc==1)then
          write (6,"('u_u1 ',10f8.2)") (u(iq,nlv),iq=idjd-3,idjd+3)
          write (6,"('v_u1 ',10f8.2)") 
     &               (v(iq,nlv),iq=idjd-3*its,idjd+3*its,its)
        endif
      else
        dumc=cc(1:ifull,:)
        dumd=dd(1:ifull,:)
        call unstaguv(dumc,dumd,dumu,dumv) ! usual
        u(1:ifull,:)=dumu
        v(1:ifull,:)=dumv
      endif

#ifdef loadbalall
      call start_log(adjustj_loadbal_begin)
      call phys_loadbal
      call end_log(adjustj_loadbal_end)
#endif

!     vert. integ. div into e
      wrk2(:,kl)=-dsig(kl)*d(:,kl)
      wrk3(:,kl)=sig(kl)*d(:,kl)
      do k=kl-1,1,-1
        wrk2(:,k)=wrk2(:,k+1)-dsig(k)*d(:,k)
        wrk3(:,k)=wrk2(:,k+1)-(sigmh(k+1)-sig(k))*d(:,k)
      enddo     ! k  loop

!     full-level omega/ps into omgf (linear part only)
      do k=1,kl-1
        omgf(:,k)=-wrk3(:,k) ! in Eq. 110
      enddo     ! k  loop
      omgf(:,kl)=-wrk3(:,kl)
      psl(1:ifull)=pslxint(:)-hdt*wrk2(:,1)  *(1.+epst(:))  ! Eq. 116
      ps_sav(1:ifull)=ps(1:ifull)  ! saved for gas fixers below, and diags
      if(mfix==-1.or.mfix==3)pslsav(1:ifull)=psl(1:ifull) 

      if (mod(ktau,nmaxpr)==0)vx(1:ifull,:)=sdot(1:ifull,1:kl) ! for accln

      if(m<=6)then
        do k=1,kl
         do iq=1,ifull
!         save [D + dsigdot/dsig] in pslx for next use in nonlin
          pslx(iq,k)=(pslx(iq,k)-psl(iq)*2./dt)/(1.+epst(iq)) ! from Eq. 115
         enddo  !  iq loop
        enddo   !  k loop
        do k=kl,2,-1
         do iq=1,ifull
!         calculate latest sdot (at level k-.5)
          sdot(iq,k)=sdot(iq,k+1)-dsig(k)*(pslx(iq,k)-d(iq,k))
         enddo  !  iq loop
        enddo   !  k loop
      endif     ! (m<=6)
      
      do k=2,kl
       do iq=1,ifull
!       and convert sdot (at level k-.5) to units of grid-steps/timestep
!       dtin is used to be ready for next full timestep
        sdot(iq,k)=sdot(iq,k)*dtin/(sig(k)-sig(k-1))
       enddo   ! iq loop
      enddo    ! k  loop

      if (mod(ktau,nmaxpr)==0)then
        vx(1:ifull,:)=sdot(1:ifull,1:kl)-vx(1:ifull,:)
!       convert to approx m/s/s
        do k=2,kl
         vx(1:ifull,k)=vx(1:ifull,k)*rdry*(sig(k-1)-sig(k))*
     &      .5*(t(1:ifull,k)+t(1:ifull,k-1))/(sigmh(k)*grav*dtin*dt)        
        enddo
        call maxmin(vx,'ac',ktau,100.,kl)  ! max min of accln * 100
      endif

      if ((diag.or.nmaxpr==1) .and. mydiag ) then
         write(6,*) 'm ',m
         if(m<=6)write (6,"('diva5p ',5p10f8.2)") d(idjd,:)
         write (6,"('omgf_l*dt',10f8.4)") omgf(idjd,:)*dt
         write (6,"('omgfnl*dt',10f8.4)") omgfnl(idjd,:)*dt
         write (6,"('u_a2 ',10f8.2)") u(idjd,:)
         write (6,"('v_a2 ',10f8.2)") v(idjd,:)
         write(6,*) 'ps,psl ',ps(idjd),psl(idjd)
         write (6,"('pslx_3p',3p9f8.4)") pslx(idjd,:)
         write (6,"('sdot_a2',10f8.3)") sdot(idjd,1:kl)
      endif

      do k=1,kl
       do iq=1,ifull
!       save full omega/ps in dpsldt for use in nonlin next time step (& outfile)
!       N.B. omgfnl part already incorp. into tx above
        dpsldt(iq,k)=omgfnl(iq,k)/(1.+epst(iq))+omgf(iq,k)
        t(iq,k)=tx(iq,k)              ! Eq 121   F26
     &           +hdt*(1.+epst(iq))*tbar2d(iq)*omgf(iq,k)*roncp/sig(k)
       enddo    ! iq loop
      enddo     ! k  loop
      if ((diag.or.nmaxpr==1).and.mydiag)then
        write(6,"('omgf_a2',10f8.3)") ps(idjd)*dpsldt(idjd,1:kl)
      endif

      if(nh/=0)then
!       update phi for use in next time step
        do k=1,kl
         phi(:,k)=p(1:ifull,k)-rdry*tbar2d(:)*psl(1:ifull)
        enddo
       
        ! extract non-hydrostatic component
        wrk3(:,1)=zs(1:ifull)+bet(1)*(t(1:ifull,1)-280.)
        phi_nh(:,1)=phi(:,1)-wrk3(:,1)
        do k=2,kl
          wrk3(:,k)=wrk3(:,k-1)+bet(k)*(t(1:ifull,k)-280.)
     &                         +betm(k)*(t(1:ifull,k-1)-280.)
          phi_nh(:,k)=phi(:,k)-wrk3(:,k)
        end do

        ! correct for temperature offste
        dum=bet(1)*280.
        phi(1:ifull,1)=phi(1:ifull,1)+dum
        do k=2,kl
          dum=dum+(bet(k)+betm(k))*280.
          phi(:,k)=phi(:,k)+dum
        end do
        
        if(nmaxpr==1.and.mydiag)then
          write(6,*) 'phi_adj ',(phi(idjd,k),k=1,kl)
        endif
      endif  ! (nh.ne.0)

#ifdef loadbalall
      call start_log(adjustk_loadbal_begin)
      call phys_loadbal
      call end_log(adjustk_loadbal_end)
#endif

      if(nvadh==2.and.nvad>0)then                 ! final dt/2 's worth
        if ((diag.or.nmaxpr==1) .and. mydiag ) then
         write(6,*) 'before vertical advection in adjust5 for ktau= ',
     &               ktau
         write (6,"('t_a2  ',10f8.2)") t(idjd,:)
         write (6,"('us# ',9f8.2)") diagvals(cc(:,nlv)) 
         write (6,"('vs# ',9f8.2)") diagvals(dd(:,nlv))  
         write (6,"('u#  ',9f8.2)") diagvals(u(:,nlv)) 
         write (6,"('v#  ',9f8.2)") diagvals(v(:,nlv)) 
         write (6,"('u_a2 ',10f8.2)") u(idjd,:)
         write (6,"('v_a2 ',10f8.2)") v(idjd,:)
         write (6,"('qg_a2 ',3p10f8.3)") qg(idjd,:)
        endif
        if(nvad==4 .or. nvad==9 ) then
         if(mup==-3)call updps(1)
          sdmx(:) = maxval(abs(sdot),2)
          nits(:)=1+sdmx(:)/nvadh
          nvadh_pass(:)=nvadh*nits(:)
          if (mydiag.and.mod(ktau,nmaxpr)==0)
     &      write(6,*) 'in adjust5 sdmx,nits,nvadh_pass ',
     &                  sdmx(idjd),nits(idjd),nvadh_pass(idjd)
          ! For now use this form of call so that vadvtvd doesn't need to 
          ! be changed. With assumed shape arguments this wouldn't be necessary
          dumt=t(1:ifull,:)
          dumu=u(1:ifull,:)
          dumv=v(1:ifull,:)
          call vadvtvd(dumt,dumu,dumv,nvadh_pass,nits,iaero)
          t(1:ifull,:)=dumt
          u(1:ifull,:)=dumu
          v(1:ifull,:)=dumv 
        endif  !  nvad==4 .or. nvad==9
        if(nvad>=7) call vadv30(t(1:ifull,:),
     &                            u(1:ifull,:),
     &                            v(1:ifull,:),iaero) ! for vadvbess
        if (( diag.or.nmaxpr==1) .and. mydiag ) then
          write(6,*) 'after vertical advection in adjust5'
          write (6,"('qg_a3',3p10f8.3)") qg(idjd,:)
          write (6,"('t_a3',10f8.2)") t(idjd,:)
          write (6,"('thet_a3',10f8.2)") t(idjd,:)*sig(:)**(-roncp)
          write (6,"('u_a3',10f8.2)") u(idjd,:)
          write (6,"('v_a3',10f8.2)") v(idjd,:)
        endif
      endif     !  (nvadh==2.and.nvad>0)

#ifdef loadbalall
      call start_log(adjustl_loadbal_begin)
      call phys_loadbal
      call end_log(adjustl_loadbal_end)
#endif

      if (mfix==-1) then   ! perform conservation fix on psl
!        delpos is the sum of all positive changes over globe
!        delneg is the sum of all negative changes over globe
!        alph_p is chosen to satisfy alph_p*delpos + delneg/alph_p = 0
!            _l means local to this processor        
         delps(1:ifull) = psl(1:ifull)-pslsav(1:ifull)
         call ccglobal_posneg(delps,delpos,delneg)
	 if(ntest==1)then
           if(myid==0)then
              write(6,*) 'psl_delpos,delneg ',delpos,delneg
              write(6,*) 'ps,psl,delps ',ps(idjd),psl(idjd),delps(idjd)
           endif
           call ccglobal_sum(pslsav,sumsav)
           call ccglobal_sum(psl,sumin)
	 endif  ! (ntest==1)
         alph_p = sqrt( -delneg/delpos)
         alph_pm=1./alph_p
         do iq=1,ifull
            psl(iq) = pslsav(iq) +
     &           alph_p*max(0.,delps(iq)) + alph_pm*min(0.,delps(iq))
         enddo
	 if(ntest==1)then
           if(myid==0)write(6,*) 'alph_p,alph_pm ',alph_p,alph_pm
           call ccglobal_sum(psl,sumout)
           if(myid==0)
     &       write(6,*) 'psl_sumsav,sumin,sumout ',sumsav,sumin,sumout
	 endif  ! (ntest==1)
      endif                     !  (mfix==-1)

      if(mfix/=3)ps(1:ifull)=1.e5*exp(psl(1:ifull))     
      if(mfix==1.or.mfix==2) then   ! perform conservation fix on ps
!         fix is on ps (not psl) from 24/1/06      
!         delpos is the sum of all positive changes over globe
!         delneg is the sum of all negative changes over globe
!         alph_p is chosen to satisfy alph_p*delpos + delneg/alph_p = 0
!             _l means local to this processor     
         bb(1:ifull)=ps(1:ifull)   
         delps(1:ifull) = ps(1:ifull)-ps_sav(1:ifull)
         call ccglobal_posneg(delps,delpos,delneg)
         if(ntest==1)then
           if(myid==0)then
              write(6,*) 'psl_delpos,delneg ',delpos,delneg
              write(6,*) 'ps,psl,delps ',ps(idjd),psl(idjd),delps(idjd)
           endif
           call ccglobal_sum(ps_sav,sumsav)
           call ccglobal_sum(ps,sumin)
         endif  ! (ntest==1)
         if(mfix==1)then
            alph_p = sqrt( -delneg/max(1.e-20,delpos))
            alph_pm=1./max(1.e-20,alph_p)
         endif                  ! (mfix==1)
         if(mfix==2)then
            if(delpos>-delneg)then
               alph_p =1.
               alph_pm=-delpos/delneg
            else
               alph_p=-delneg/delpos
               alph_pm =1.
            endif
         endif                  ! (mfix==2)
         do iq=1,ifull
            ps(iq) = ps_sav(iq) +
     &           alph_p*max(0.,delps(iq)) + alph_pm*min(0.,delps(iq))
         enddo
        if(ntest==1)then
           if(myid==0) write(6,*) 'alph_p,alph_pm ',alph_p,alph_pm
           call ccglobal_sum(ps,sumout)
           if (myid==0)
     &       write(6,*) 'ps_sumsav,sumin,sumout ',sumsav,sumin,sumout
        endif  ! (ntest==1)
!       psl(1:ifull)=log(1.e-5*ps(1:ifull))
!       following is cheaper and maintains full precision of psl
        psl(1:ifull)=psl(1:ifull)+  ! just using mass flux increment
     &              (ps(1:ifull)/bb(1:ifull)-1.)     
      endif                   !  (mfix==1.or.mfix==2)
      
      if(mfix==3) then   ! perform conservation fix on ps (best for 32-bit)
!       fix is on ps (not psl) from 24/1/06      
!       delpos is the sum of all positive changes over globe
!       delneg is the sum of all negative changes over globe
!       alph_p is chosen to satisfy alph_p*delpos + delneg/alph_p = 0
!           _l means local to this processor     
        delps(1:ifull)=psl(1:ifull)-pslsav(1:ifull)
        delps(1:ifull)=ps_sav(1:ifull)*
     &                 delps(1:ifull)*(1.+.5*delps(1:ifull))         
        call ccglobal_posneg(delps,delpos,delneg)
        if(ntest==1)then
          if(myid==0)then
            write(6,*) 'psl_delpos,delneg ',delpos,delneg
            write(6,*) 'ps,psl,delps ',ps(idjd),psl(idjd),delps(idjd)
            write(6,*) 'ps_sav,pslsav ',ps_sav(idjd),pslsav(idjd)
          endif
          call ccglobal_sum(ps_sav,sumsav)
          call ccglobal_sum(ps,sumin)
        endif  ! (ntest==1)
        alph_p = sqrt( -delneg/max(1.e-20,delpos))
        alph_pm=1./max(1.e-20,alph_p)
        do iq=1,ifull
         delps(iq)=alph_p*max(0.,delps(iq)) + alph_pm*min(0.,delps(iq))
        enddo
        delps(1:ifull)=delps(1:ifull)/ps_sav(1:ifull)
        psl(1:ifull)=pslsav(1:ifull)+
     &               delps(1:ifull)*(1.-.5*delps(1:ifull))
        ps(1:ifull)=1.e5*exp(psl(1:ifull))     
	  if(ntest==1)then
          if(myid==0) write(6,*) 'alph_p,alph_pm ',alph_p,alph_pm
          call ccglobal_sum(ps,sumout)
          if(myid==0)
     &      write(6,*) 'ps_sumsav,sumin,sumout ',sumsav,sumin,sumout
	  endif  ! (ntest==1)
      endif                   !  (mfix==3)
      
!     following dpsdt diagnostic is in hPa/day  
      dpsdtbb(:)=dpsdtb(:)    
      dpsdtb(:)=dpsdt(:)    
      dpsdt(1:ifull)=(ps(1:ifull)-ps_sav(1:ifull))*24.*3600./(100.*dt)
      if(nmaxpr==1)call maxmin(dpsdt,'dp',ktau,.01,1)

      !--------------------------------------------------------------
      ! Moisture conservation
      if (mfix_qg/=0.and.mspec==1) then
        do k=1,kl
          qg(1:ifull,k)=max(qg(1:ifull,k),
     &      qgmin-qfg(1:ifull,k)-qlg(1:ifull,k),0.)
        end do
        dumc=qg(1:ifull,:)
        call massfix(mfix_qg,dumc,qgsav,
     &               ps(1:ifull),ps_sav,wts,
     &               .false.)
        qg(1:ifull,:)=dumc
      endif       !  (mfix_qg.ne.0.and.mspec==1)

      !------------------------------------------------------------------------
      ! Cloud water conservation
      if(mfix_qg/=0.and.mspec==1.and.ldr/=0)then
        qfg=max(qfg,0.)
        qlg=max(qlg,0.)
        qrg=max(qrg,0.)
        cfrac=min(max(cfrac,0.),1.)
        cffall=min(max(cffall,0.),1.)
        dumc=qfg(1:ifull,:)
        call massfix(mfix_qg,dumc,qfgsav,
     &               ps(1:ifull),ps_sav,wts,
     &               .true.)
        qfg(1:ifull,:)=dumc
        dumc=qlg(1:ifull,:)
        call massfix(mfix_qg,dumc,qlgsav,
     &               ps(1:ifull),ps_sav,wts,
     &               .true.)
        qlg(1:ifull,:)=dumc
        dumc=qrg(1:ifull,:)
        call massfix(mfix_qg,dumc,qrgsav,
     &               ps(1:ifull),ps_sav,wts,
     &               .true.)
        qrg(1:ifull,:)=dumc
      endif      !  (mfix_qg.ne0.and.mspec==1.and.ldr.ne.0)

      !------------------------------------------------------------------------
      ! Tracer conservation
      if(mfix_tr/=0.and.mspec==1.and.ngas>0)then
        do ng=1,ngas
!         rml 19/09/07 replace gasmin with tracmin
          tr(:,:,ng)=max(tr(:,:,ng),tracmin(ng))
          dumc=tr(1:ifull,:,ng)
          dumd=trsav(1:ifull,:,ng)
          call massfix(mfix_tr,dumc,dumd,
     &                 ps(1:ifull),ps_sav,wts,
     &                 .true.)
          tr(1:ifull,:,ng)=dumc
        end do
      endif       !  (mfix_tr.ne.0.and.mspec==1.and.ngas>0)

      !--------------------------------------------------------------
      ! Aerosol conservation
      if (mfix_aero/=0.and.mspec==1.and.abs(iaero)==2) then
        xtg=max(xtg,0.)
        do ng=1,naero
          dumc=xtg(1:ifull,:,ng)
          dumd=xtgsav(1:ifull,:,ng)
          call massfix(mfix_aero,dumc,
     &                 dumd,ps(1:ifull),
     &                 ps_sav,wts,.true.)
          xtg(1:ifull,:,ng)=dumc
        end do
      end if
      !--------------------------------------------------------------

      if ((diag.or.nmaxpr==1) .and. mydiag ) then
        write(6,*) 'at end of adjust5 for ktau= ',ktau
        write(6,*) 'ps_sav,ps ',ps_sav(idjd),ps(idjd)
        write (6,"('dpsdt# mb/d ',9f8.1)") diagvals(dpsdt) 
        write (6,"('qg_a4 ',3p10f8.3)") qg(idjd,:)
        write (6,"('qgs',3p10f8.3)")
     &                         (qgsav(idjd,k)/ps_sav(idjd),k=1,kl)
        write (6,"('qf_a4',3p10f8.3)") qfg(idjd,:)
        write (6,"('ql_a4',3p10f8.3)") qlg(idjd,:)
      endif
      if (diag)then
         call bounds(psl)
         if (mydiag ) then
            iq=idjd
            write(6,*) 'adjust5 d ',d(iq,:)
            write(6,*) 'adjust5 wrk2 ',idjd, wrk2(idjd,:)
            write (6,"('adjust5 psl_3p & n e w s',3p9f8.4)") psl(idjd)
     &          ,psl(in(iq)),psl(ie(iq)),psl(iw(iq)),psl(is(iq))
         end if
         call printa('psl ',psl,ktau,0,ia,ib,ja,jb,0.,1000.)
         if (mydiag ) write(6,*) 'adjust5 t ',t(idjd,:)
         call printa('t   ',t,ktau,nlv,ia,ib,ja,jb,200.,1.)
         if (mydiag ) write(6,*) 'adjust5 u ',u(idjd,:)
         call printa('u   ',u,ktau,nlv,ia,ib,ja,jb,0.,1.)
         if (mydiag ) write(6,*) 'adjust5 v ',v(idjd,:)
         call printa('v   ',v,ktau,nlv,ia,ib,ja,jb,0.,1.)
         if (mydiag ) write(6,*) 'ps diagnostic from end of adjust:'
         do iq=1,ifull
            bb(iq)=ps(iq)-ps_sav(iq)
         enddo
         call printa('dps ',bb,ktau,0,ia,ib,ja,jb,0.,.01)
         call printa('ps  ',ps,ktau,0,ia,ib,ja,jb,1.e5,.01)
         if (sig(nlv)<.3)then
            call printa('qg  ',qg,ktau,nlv,ia,ib,ja,jb,0.,1.e6)
         else
            call printa('qg  ',qg,ktau,nlv,ia,ib,ja,jb,0.,1.e3)
         endif
      endif

      dtsave = dt
      
#ifdef loadbalall
      call start_log(adjustm_loadbal_begin)
      call phys_loadbal
      call end_log(adjustm_loadbal_end)
#endif
      
      call end_log(adjust_end)

      end subroutine adjust5

      subroutine adjust_init(zz,zzn,zze,zzw,zzs,pfact,alff,alf,alfe,
     &                       alfn,alfu,alfv)
      use cc_mpi
      use indices_m
      use map_m
      implicit none
      include 'newmpar.h'
      include 'parm.h'
      include 'parmdyn.h'
      real zz(ifull),zzn(ifull),zze(ifull),zzw(ifull),
     &     zzs(ifull),pfact(ifull),alff(ifull+iextra),alf(ifull+iextra),
     &     alfe(ifull+iextra),alfn(ifull+iextra),alfu(ifull),alfv(ifull)
      real :: hdt
      integer :: iq, n

      hdt = 0.5*dt
      if(m<3)then  ! 14/7/09 getting ready for gnomonic options
        do iq=1,ifull
         alf(iq)=1.+epsu
         alff(iq)=0.
         alfu(iq)=1.   ! i.e. alf/(1+epsu)
         alfv(iq)=1.   ! i.e. alf/(1+epsu)
        enddo     ! iq loop
      else
        do iq=1,ifull
         alf(iq)=(1.+epsu)/(1.+(hdt*(1.+epsf)*f(iq))**2) ! Eq. 138
         alfF(iq)=alf(iq)*f(iq)*hdt*(1.+epsf)       ! now includes hdt in alfF
         alfu(iq)=1./(1.+(hdt*(1.+epsf)*fu(iq))**2) ! i.e. alf/(1+epsu)
         alfv(iq)=1./(1.+(hdt*(1.+epsf)*fv(iq))**2) ! i.e. alf/(1+epsu)
        enddo     ! iq loop
      endif  ! (m<5)... else ...
      ! These really only need to be recomputed when time step changes.
      call bounds(alf)
      call bounds(alff,corner=.true.)
!cdir nodep
      do iq=1,ifull
       pfact(iq)=4.*( ds/(dt*em(iq)) )**2    
       alfe(iq)=alf(ie(iq))-alf(iq)    ! alf3 Eq. 141 times ds
     &     +.25*(alff(ine(iq))+alfF(in(iq))-alfF(ise(iq))-alfF(is(iq)))
       alfn(iq)=alf(in(iq))-alf(iq)    ! alf4 Eq. 142 times ds
     &     -.25*(alfF(ien(iq))+alfF(ie(iq))-alfF(iwn(iq))-alfF(iw(iq)))
      enddo     ! iq loop
      call boundsuv(alfe,alfn)

!cdir nodep
      do iq=1,ifull    ! see Eq 158/159
!      need care with vector quantities on w (odd) & s (even) panel boundaries
       zz(iq)=.5*(alfe(iwu(iq))-alfe(iq)+alfn(isv(iq))-alfn(iq))
     &                     -4.*alf(iq)                 ! i,j   coeff
       zzn(iq)=alf(in(iq))-.5*alfn(iq)                 ! i,j+1 coeff
       zzw(iq)=alf(iw(iq))+.5*alfe(iwu(iq))            ! i-1,j coeff
       zze(iq)=alf(ie(iq))-.5*alfe(iq)                 ! i+1,j coeff
       zzs(iq)=alf(is(iq))+.5*alfn(isv(iq))            ! i,j-1 coeff
      enddo     ! iq loop
!     N.B. there are some special z values at the 8 vertices
         do n=1,npan            ! 1,6
            if(edge_s(n-noff) .and. edge_w(n-noff) ) then
               iq=indp(1,1,n)
               zzs(iq)=zzs(iq)+.25*alfF(is(iq)) ! i,j-1 coeff
               zzw(iq)=zzw(iq)-.25*alfF(iw(iq)) ! i-1,j coeff
            end if
            if(edge_n(n-noff) .and. edge_e(n-noff) ) then
               iq=indp(ipan,jpan,n)
               zzn(iq)=zzn(iq)+.25*alfF(in(iq)) ! i,j+1 coeff
               zze(iq)=zze(iq)-.25*alfF(ie(iq)) ! i+1,j coeff
            end if
            if(edge_s(n-noff) .and. edge_e(n-noff) ) then
               iq=indp(ipan,1,n)
               zzs(iq)=zzs(iq)-.25*alfF(is(iq)) ! i,j-1 coeff
               zze(iq)=zze(iq)+.25*alfF(ie(iq)) ! i+1,j coeff
            end if
            if(edge_n(n-noff) .and. edge_w(n-noff) ) then
               iq=indp(1,jpan,n)
               zzn(iq)=zzn(iq)-.25*alfF(in(iq)) ! i,j+1 coeff
               zzw(iq)=zzw(iq)+.25*alfF(iw(iq)) ! i-1,j coeff
            end if
        enddo   ! n loop
      end subroutine adjust_init

      subroutine massfix(mfix,s,ssav,ps,pssav,wts,llim)
      
      use cc_mpi
      
      implicit none
      
      include 'newmpar.h'
      
      integer, intent(in) :: mfix
      integer k
      real, dimension(ifull,kl), intent(inout) :: s,ssav
      real, dimension(ifull), intent(in) :: ps,pssav,wts
      real, dimension(ifull,kl) :: wrk1
      real delpos,delneg,ratio,alph_g
      logical, intent(in) :: llim

      if (mfix>0) then
        do k=1,kl
          s(:,k)=s(:,k)*ps
          ssav(:,k)=ssav(:,k)*pssav
        enddo    ! k  loop           
      else
        do k=1,kl
          s(:,k)=s(:,k)*ps*wts
          ssav(:,k)=ssav(:,k)*pssav*wts
        enddo    ! k  loop     
      endif ! (mfix>0) .. else ..
      do k=1,kl
        wrk1(:,k)=s(:,k)-ssav(:,k) 
      enddo   ! k loop
      call ccglobal_posneg(wrk1,delpos,delneg)
      if (llim) then
        ratio = -delneg/max(delpos,1.e-30)
      else
        ratio = -delneg/delpos
      end if
      if (mfix==1) then
        alph_g = min(ratio,sqrt(ratio))
      elseif (mfix==2) then
        if (llim) then
          alph_g = max(sqrt(ratio),1.e-30)
        else
          alph_g = sqrt(ratio)
        end if
      end if ! (mfix==1) .. else ..
      do k=1,kl
        s(1:ifull,k)=ssav(1:ifull,k)+
     &    alph_g*max(0.,wrk1(:,k))+min(0.,wrk1(:,k))/max(1.,alph_g)
      enddo    ! k  loop
      if (mfix>0) then
        do k=1,kl
          s(1:ifull,k)=s(1:ifull,k)/ps
          ssav(1:ifull,k)=ssav(1:ifull,k)/pssav
        enddo    ! k  loop
      else
        do k=1,kl
          s(1:ifull,k)=s(1:ifull,k)/(ps*wts)
          ssav(1:ifull,k)=ssav(1:ifull,k)/(pssav*wts)
        enddo   ! k  loop            
      endif ! (mfix>0) .. else ..

      return
      end subroutine massfix
