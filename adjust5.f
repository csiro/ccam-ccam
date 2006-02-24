      subroutine adjust5
      use cc_mpi
      use diag_m
      implicit none
      integer, parameter :: mfix_rad=0 ! used to make gases 2 to ng add up to gas 1
      integer, parameter :: ntest=0
      integer, parameter :: nys=0      ! for nys=1 use sun & yeh 1997 monotonic filter
!                              within adjust5 on u,v
!     nuv in parm.h:  now always 10
!          !  0y for increments   , 1y for actual values, 2y for true increments
      include 'newmpar.h'
      include 'arrays.h'
      include 'const_phys.h'
      include 'indices.h'
      include 'kuocom.h'   ! ldr
      include 'liqwpar.h'  ! qfg,qlg
      include 'map.h'
      include 'morepbl.h'  ! condx,eg
      include 'nlin.h'
      include 'parm.h'     ! qgmin
      include 'parmdyn.h'  
      include 'parmvert.h'  
      include 'pbl.h'
      include 'sigs.h'
      include 'tracers.h'
      include 'vecs.h'
      include 'vecsuv.h'   ! vecsuv info
      include 'vecsuva.h'  ! vecsuva info
      include 'vvel.h'     ! sdot
      include 'xarrs.h'
      include 'xyzinfo.h'
      include 'mpif.h'
      real dpsdt,dpsdtb,dpsdtbb,epst
      common/dpsdt/dpsdt(ifull),dpsdtb(ifull),dpsdtbb(ifull) !globpe, adjust5,outcdf
      common/epst/epst(ifull)
      real alph_p, alph_pm, delneg, delpos, delnegk, delposk, alph_q
      common/mfixdiag/alph_p,alph_pm,delneg,delpos,alph_q
      real :: sumdiffb_l  ! Local versions
      real, dimension(3) :: delarr, delarr_l
      common/nharrs/phi(ifull,kl),h_nh(ifull+iextra,kl)
      real phi, h_nh, tbar2d
      common/tbar2d/tbar2d(ifull)
      ! Work common here ????
      real p(ifull+iextra,kl),omgfnl(ifull,kl)
      real wrk1, wrk2  ! wrk2 not used here
      common/work3b/wrk1(ifull,kl),wrk2(ifull,kl)   ! just work arrays here
      real d(ifull,kl)   ! NOT shared updps
      real qgsav, qfgsav, qlgsav, trsav
      common/work3sav/qgsav(ifull,kl),qfgsav(ifull,kl),qlgsav(ifull,kl)
     &             ,trsav(ilt*jlt,klt,ngasmax)  ! shared adjust5 & nonlin
      real zz(ifull),zzn(ifull),zze(ifull),zzw(ifull),
     &     zzs(ifull),pfact(ifull),alff(ifull+iextra),alf(ifull+iextra),
     &     alfe(ifull+iextra),alfn(ifull+iextra),alfu(ifull),alfv(ifull)
      common /zzalf/ zz,zzn,zze,zzw,zzs,pfact,alff,alf,alfe,alfn,
     &     alfu,alfv
      real ps_sav(ifull),cc(ifull+iextra,kl),
     &     dd(ifull+iextra,kl),pslxint(ifull)
      real pe(ifull+iextra,kl),e(ifull,kl)
      real helm(ifull+iextra,kl),rhsl(ifull+iextra,kl),delps(ifull)
      real pextras(ifull,kl),omgf(ifull,kl)
      real bb(ifull),pse(ifull+iextra),psn(ifull+iextra)
      real fluxe(ifull+iextra),fluxn(ifull+iextra)
!     This is a genuine rather than space saving equivalence
      equivalence (pextras,dpsldt)
!     Save this so we can check whether initialisation needs to be redone
      real, save :: dtsave = 0.0
      real, save, dimension(kl) :: accel
      real :: hdt, hdtds, sdmx, sdmx_g, sum, qgminm, ratio, sumdiffb,
     &        alph_g
      integer :: ii, jj, its, k, l, nits, nvadh_pass, iq, ng, ierr
      integer, save :: precon_in
      real :: sumin, sumout, sumsav
      real :: delpos_l, delneg_l, const_nh

      call start_log(adjust_begin)
      hdt=dt/2.
      hdtds=hdt/ds

      if ( dt /= dtsave ) then 
         call adjust_init
         precon_in=precon
         precon=0
      end if

      if(diag)then
         if ( mydiag ) then
           write (6,"('tx_a1',10f8.2)") tx(idjd,:)
           if(m<8)then
             write (6,"('ux_stag',10f8.2)") ux(idjd,:)
             write (6,"('vx_stag',10f8.2)") vx(idjd,:)
           else
             write (6,"('ux_stag',-5p10f8.2)") ux(idjd,:)
             write (6,"('vx_stag',-5p10f8.2)") vx(idjd,:)
           endif
           write (6,"('qg_a1',3p10f7.3)") qg(idjd,:)
           print  *,'pslx ',pslx(idjd,:)
         end if
         call printa('pslx',pslx,ktau,nlv,ia,ib,ja,jb,0.,100.)
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
      do k=kl-1,1,-1
       e(:,k)=e(:,k+1)+dsig(k)*pslx(1:ifull,k)  ! i.e. -M_bar_(sig-.5)
      enddo     ! k loop
      pslxint(:)=-e(:,1)*dt/2. ! pslxint holds integrated pslx, Eq. 116

!     full-level (1.+epsp)*omega/ps into omgfnl (nonlin. part only), Eq. 118
      do k=1,kl-1
       omgfnl(1:ifull,k)=-rata(k)*e(:,k+1)-ratb(k)*e(:,k)
     &                   -sig(k)*pslx(1:ifull,k)
      enddo     ! k loop
      omgfnl(1:ifull,kl)=-ratb(kl)*e(:,kl)-sig(kl)*pslx(1:ifull,kl)
     
!     redefine ux, vx
      do k=1,kl
         do iq=1,ifull  ! N.B. alfu is alf/(1+epsu)
            cc(iq,k)=ux(iq,k)/emu(iq)*alfu(iq)  ! Eq. 136
            dd(iq,k)=vx(iq,k)/emv(iq)*alfv(iq)  ! Eq. 137
         end do
      end do
      call boundsuv(cc,dd)

      do k=1,kl
       do iq=1,ifull
!       N.B. the omgfnl term on LHS of Eq. 121 not yet added in
        tx(iq,k)=tx(iq,k)+hdt*                ! adding LHS term of Eq. 121
     &           tbar2d(iq)*omgfnl(iq,k)*roncp/sig(k)  ! with correct epsp
       enddo    ! iq loop
      enddo     ! k loop

!     calculate heights from the tx array
      do iq=1,ifull
       p(iq,1)=zs(iq)+bet(1)*tx(iq,1)+rdry*tbar2d(iq)*pslxint(iq) ! Eq. 146
      enddo     ! iq loop
      do k=2,kl
       do iq=1,ifull
        p(iq,k)=p(iq,k-1)+bet(k)*tx(iq,k)+betm(k)*tx(iq,k-1)
       enddo    ! iq loop
      enddo     ! k loop

      if(nh>0)then
        const_nh=-2.*rdry/((1.+epsnh)*dt*grav*grav) ! bet inc. an r term        
!       add in departure values of p-related nh terms  & omgfnl terms    
c        wrk1(:,1)=bet(1)*const_nh*tbar(1)*(h_nh(1:ifull,1)
c     .                 -sig(1)*omgfnl(:,1)/(tbar(1)*(1.+epsnh)) )        
c        do k=2,kl
c         wrk1(:,k)=wrk1(:,k-1)+const_nh*tbar(k)*(
c     .             bet(k)*h_nh(1:ifull,k)+betm(k)*h_nh(1:ifull,k-1) 
c     .            -bet(k)*sig(k)*omgfnl(:,k)/(tbar(k)*(1.+epsnh))         
c     .            -betm(k)*sig(k-1)*omgfnl(:,k-1)/(tbar(k)*(1.+epsnh)) )
c        enddo   ! k loop
        wrk1(:,1)=bet(1)*const_nh*tbar(1)*(h_nh(1:ifull,1)
     .                 -tbar2d(:)*omgfnl(:,1)/(sig(1)*(1.+epst(:))) )        
        do k=2,kl
         wrk1(:,k)=wrk1(:,k-1)+const_nh*tbar(k)*(
     .        bet(k)*h_nh(1:ifull,k)+betm(k)*h_nh(1:ifull,k-1) 
     .       -bet(k)*tbar2d(:)*omgfnl(:,k)/(sig(k)*(1.+epst(:)))         
     .       -betm(k)*tbar2d(:)*omgfnl(:,k-1)/(sig(k-1)*(1.+epst(:))) )
        enddo   ! k loop
        if(diag.and.mydiag)then   
          print *,'adjust5 omgfnl ',(omgfnl(idjd,k),k=1,kl)
          print *,'adjust5 h_nh ',(h_nh(idjd,k),k=1,kl)
          print *,'adjust5 pa ',(p(idjd,k),k=1,kl)
          print *,'adjust5 wrk1 ',(wrk1(idjd,k),k=1,kl)
          print *,'adjust5 pextras ',(pextras(idjd,k),k=1,kl)
        endif
        p(1:ifull,:)=p(1:ifull,:)+wrk1(:,:)  ! nh
      endif     ! (nh>0)

c     do k=1,kl
c       do iq=1,ifull
c        p(iq,k)=p(iq,k)+pextras(iq,k)   ! npgf=0 code (see nonlin)
c       enddo   ! iq loop
c      enddo    ! k  loop

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

!         Diagnostics would require extra bounds calls
!         if(diag.and.k<=2)then !  only for last k of loop (i.e. 1)
!            iq=idjd
!            print  *,'adjust5(k) p & n e w s ',k,p(iq,k),
!     &        p(in(iq),k),p(ie(iq),k),p(iw(iq),k),p(is(iq),k)
!            print  *,'adjust5(k) pe & n e w s ',k,pe(iq,k),
!     &        pe(in(iq),k),pe(ie(iq),k),pe(iw(iq),k),pe(is(iq),k)
!            print  *,'adjust5(k) rhsl & n e w s ',k,rhsl(iq,k),
!     &       rhsl(in(iq),k),rhsl(ie(iq),k),rhsl(iw(iq),k),rhsl(is(iq),k)
!         endif                  ! (diag.and.k<=2)
      enddo    ! k loop

      call bounds(pe)
      if(precon_in>kl)then  ! use SOR - needs checking
!       Calculate the optimum acceleration using a point
!       in the middle of the first face
        if(nproc>1)then
          write(0,*) 'using helmsor/SOR requires nproc=1'
          stop
        endif
        do k=1,kl
         if ( dt /= dtsave ) then
           call optmx(il,schmidt,dt,bam(k),accel(k))
           print *,'k,accel ',k,accel(k)
         end if
         if(precon_in==98)accel(k)=-1.
         call helmsor(accel(k),zz,zzn,zze,zzw,zzs,helm(1,k),
     &                pe(1,k),rhsl(1,k))
        enddo
        if(precon_in==98)then  ! using 1 itn of SOR as preconditioner
          call bounds(pe)
          call helmsol(zz,zzn,zze,zzw,zzs,helm,pe,rhsl)
        endif
      else
        call helmsol(zz,zzn,zze,zzw,zzs,helm,pe,rhsl)
      endif  ! (precon>kl)

      if(diag)then   !  only for last k of loop (i.e. 1)
         ! Some diagnostics that would require extra bounds calls have been removed
         if ( mydiag ) then
            write (6,"('tx_a2',10f8.2)") tx(idjd,:)
            print *,'omgfnl_a2 ',omgfnl(idjd,:)
            print *,'adjust5 cc ',cc(idjd,:)
            print *,'adjust5 dd ',dd(idjd,:)
            print *,'adjust5 d_in ',d(idjd,:)
            print *,'adjust5 p ',p(idjd,:)
            print *,'adjust5 pe ',pe(idjd,:)
            print *,'adjust5 pe_e ',pe(ie(idjd),:)
            print *,'adjust5 pe_w ',pe(iw(idjd),:)
            print *,'adjust5 rhsl ',rhsl(idjd,:)
            print *,'adjust5 helm ',helm(idjd,:)
         end if
         call printa('psnt',pslxint,ktau,0,ia,ib,ja,jb,0.,100.)
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
      call bounds(p,corner=.true.)

!      now u & v
      if(diag.and.mydiag)then
         print*,'iq,k,fu,alfu,alfu*ux(iq,k) ',
     &           idjd,nlv,fu(idjd),alfu(idjd),alfu(idjd)*ux(idjd,nlv)
         print*,'alfF & n e w s (in(iq)),alfF(ine(iq)),alfF(is(iq))',
     &          'alfF(ise(iq)),alfe(iq) ',alfF(in(idjd)),alfF(ine(idjd))
     &           ,alfF(is(idjd)),alfF(ise(idjd)),alfe(idjd)
         sum = alf(ie(idjd))-alf(idjd)+.25*(alfF(in(idjd))+
     &        alfF(ine(idjd))-alfF(is(idjd))-alfF(ise(idjd)))-alfe(idjd)
         print*,'sum  ',sum
         print*,'p & n e w s ne se ',p(idjd,nlv),p(in(idjd),nlv),
     &           p(ie(idjd),nlv),p(iw(idjd),nlv),p(is(idjd),nlv),
     &           p(ine(idjd),nlv),p(ise(idjd),nlv)
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
         enddo    
      end do

      call boundsuv(cc,dd)
!     calculate linear part only of sigma-dot and omega/ps
      do k=1,kl
!cdir nodep
         do iq=1,ifull
            d(iq,k)=(cc(iq,k)/emu(iq)-cc(iwu(iq),k)/emu(iwu(iq)) ! Eq. 101
     &          +dd(iq,k)/emv(iq)-dd(isv(iq),k)/emv(isv(iq)))    
     &          *em(iq)**2/ds
         enddo                  ! iq loop
      enddo     ! k  loop
      if(m==7)then
        do k=1,kl
         cc(1:ifull,k)=cc(1:ifull,k)/pse(1:ifull)
         dd(1:ifull,k)=dd(1:ifull,k)/psn(1:ifull)
        enddo
      endif  ! (m==7)

!     straightforward rev. cubic interp of u and v (i.e. nuv=10)
!     This is necessary because staguv expects arrays dimensioned 
!     (ifull,kl). 
      call unstaguv(cc(1:ifull,:),dd(1:ifull,:),
     &              u(1:ifull,:),v(1:ifull,:)) ! usual

!     vert. integ. div into e
      do iq=1,ifull
       wrk2(iq,kl)=-dsig(kl)*d(iq,kl)  ! nh
      enddo     ! iq loop
      do k=kl-1,1,-1
       do iq=1,ifull
        wrk2(iq,k)=wrk2(iq,k+1)-dsig(k)*d(iq,k)
       enddo    ! iq loop
      enddo     ! k  loop

!     full-level omega/ps into omgf (linear part only)
      do k=1,kl-1
       do iq=1,ifull
        omgf(iq,k)=-rata(k)*wrk2(iq,k+1)-ratb(k)*wrk2(iq,k) ! in Eq. 110
       enddo    ! iq loop
      enddo     ! k  loop
      ps_sav(1:ifull)=ps(1:ifull)  ! saved for gas fixers below, and diags
      do iq=1,ifull
       omgf(iq,kl)=-ratb(kl)*wrk2(iq,kl)
       psl(iq)=pslxint(iq)-hdt*wrk2(iq,1)  *(1.+epst(iq))  ! Eq. 116
      enddo     ! iq loop

      if(nh.ne.0)then
!       update phi for use in next time step; check tbar or tbar2d ***    
        do k=1,kl
         phi(:,k)=p(1:ifull,k)-rdry*tbar2d(:)*psl(1:ifull)
        enddo
      endif  ! (nh.ne.0)

      if(m==6)then
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
      endif     ! (m==6)
      
      do k=2,kl
       do iq=1,ifull
!       and convert sdot (at level k-.5) to units of grid-steps/timestep
!       dtin is used to be ready for next full timestep
        sdot(iq,k)=sdot(iq,k)*dtin/(sig(k)-sig(k-1))
       enddo   ! iq loop
      enddo    ! k  loop

      if( (diag.or.nmaxpr==1) .and. mydiag ) then
         print *,'m ',m
         if(m==6)write (6,"('diva5p ',5p10f8.2)") d(idjd,:)
         write (6,"('omgf_l*dt',10f8.4)") omgf(idjd,:)*dt
         write (6,"('omgfnl*dt',10f8.4)") omgfnl(idjd,:)*dt
         write (6,"('u_a2 ',10f8.2)") u(idjd,:)
         write (6,"('v_a2 ',10f8.2)") v(idjd,:)
         print *,'ps,psl ',ps(idjd),psl(idjd)
         print *,'pslx ',pslx(idjd,:)
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
      if((diag.or.nmaxpr==1).and.mydiag)then
        write(6,"('omgf_a2',10f8.3)") ps(idjd)*dpsldt(idjd,1:kl)
      endif

      if(nvadh==2.and.nvad>0)then                 ! final dt/2 's worth
        if ( (diag.or.nmaxpr==1) .and. mydiag ) then
         print *,'before vertical advection in adjust5 for ktau= ',ktau
         write (6,"('t_a2  ',10f8.2)") t(idjd,:)
         write (6,"('us# ',9f8.2)") 
     &           ((cc(ii+jj*il,nlv),ii=idjd-1,idjd+1),jj=1,-1,-1)
         write (6,"('vs# ',9f8.2)") 
     &           ((dd(ii+jj*il,nlv),ii=idjd-1,idjd+1),jj=1,-1,-1)
         write (6,"('u#  ',9f8.2)") 
     &           ((u(ii+jj*il,nlv),ii=idjd-1,idjd+1),jj=1,-1,-1)
         write (6,"('v#  ',9f8.2)") 
     &           ((v(ii+jj*il,nlv),ii=idjd-1,idjd+1),jj=1,-1,-1)
         write (6,"('u_a2 ',10f8.2)") u(idjd,:)
         write (6,"('v_a2 ',10f8.2)") v(idjd,:)
         write (6,"('qg_a2 ',3p10f8.3)") qg(idjd,:)
        endif
        if ( nvad==4 .or. nvad==9 ) then
         if(mup==-3)call updps(1)
          sdmx = maxval(abs(sdot))
          call MPI_AllReduce(sdmx, sdmx_g, 1, MPI_REAL, MPI_MAX,
     &                       MPI_COMM_WORLD, ierr )
          nits=1+sdmx_g/nvadh
          nvadh_pass=nvadh*nits
          if(mydiag.and.mod(ktau,nmaxpr)==0)
     &      print *,'in adjust5 sdmx,nits,nvadh_pass ',
     &                          sdmx_g,nits,nvadh_pass
          do its=1,nits
           ! For now use this form of call so that vadvtvd doesn't need to 
           ! be changed. With assumed shape arguments this wouldn't be necessary
            call vadvtvd(t(1:ifull,:),u(1:ifull,:),v(1:ifull,:),
     &                   nvadh_pass)
          enddo
        endif  !  nvad==4 .or. nvad==9
        if(nvad>=7) call vadv30(t(1:ifull,:),
     &                            u(1:ifull,:),
     &                            v(1:ifull,:)) ! for vadvbess
        if( ( diag.or.nmaxpr==1) .and. mydiag ) then
          print *,'after vertical advection in adjust5'
          write (6,"('qg_a3',3p10f8.3)") qg(idjd,:)
          write (6,"('t_a3',10f8.2)") t(idjd,:)
          write (6,"('thet_a3',10f8.2)") t(idjd,:)*sig(:)**(-roncp)
          write (6,"('u_a3',10f8.2)") u(idjd,:)
          write (6,"('v_a3',10f8.2)") v(idjd,:)
        endif
      endif     !  (nvadh==2.and.nvad>0)

      ps(1:ifull)=1.e5*exp(psl(1:ifull))     
      if (mfix==1.or.mfix==2) then   ! perform conservation fix on ps
!         fix is on ps (not psl) from 24/1/06      
!         delpos is the sum of all positive changes over globe
!         delneg is the sum of all negative changes over globe
!         alph_p is chosen to satisfy alph_p*delpos + delneg/alph_p = 0
!             _l means local to this processor     
         bb(1:ifull)=ps(1:ifull)   
         delps(1:ifull) = ps(1:ifull)-ps_sav(1:ifull)
         call ccglobal_posneg(delps,delpos,delneg)
         if(ntest==1)then
            call ccglobal_sum(ps_sav,sumsav)
            call ccglobal_sum(ps,sumin)
         endif  ! (ntest==1)
         if(mfix==1)then
            alph_p = sqrt( -delneg/delpos)
            alph_pm=1./alph_p
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
           if(myid==0) print *,'ps_delpos,delneg,alph_p,alph_pm ',
     &                 delpos,delneg,alph_p,alph_pm
           call ccglobal_sum(ps,sumout)
           if(myid==0)
     &        print *,'ps_sumsav,sumin,sumout ',sumsav,sumin,sumout
        endif  ! (ntest==1)
!       psl(1:ifull)=log(1.e-5*ps(1:ifull))
!       following is cheaper and maintains full precision of psl
        psl(1:ifull)=psl(1:ifull)+  ! just using mass flux increment
     &              (ps(1:ifull)/bb(1:ifull)-1.)     
      endif                   !  (mfix==1.or.mfix==2)
      
      if(mfix==3)then
c       just code fully implicit for a start      
        call bounds(ps)       
        do iq=1,ifull
         pse(iq)=.5*(ps(iq)+ps(ie(iq)))/emu(iq)
         psn(iq)=.5*(ps(iq)+ps(in(iq)))/emv(iq)
        enddo
!       note that cc and dd still contain staggered u and v       
        fluxe(1:ifull)=-cc(1:ifull,1)*pse(1:ifull)*dsig(1) 
        fluxn(1:ifull)=-dd(1:ifull,1)*psn(1:ifull)*dsig(1) 
        do k=2,kl
        fluxe(1:ifull)=fluxe(1:ifull)-cc(1:ifull,k)*pse(1:ifull)*dsig(k)
        fluxn(1:ifull)=fluxn(1:ifull)-dd(1:ifull,k)*psn(1:ifull)*dsig(k)
        enddo
        call boundsuv(fluxe,fluxn)
        do iq=1,ifull
         bb(iq)=ps_sav(iq)-dt*(fluxe(iq)-fluxe(iwu(iq))
     &                    +fluxn(iq)-fluxn(isv(iq)))*em(iq)**2/ds
        enddo
      endif  ! (mfix==3)
      
!     following diagnostic is in hPa/day  
      dpsdtbb(:)=dpsdtb(:)    
      dpsdtb(:)=dpsdt(:)    
      dpsdt(1:ifull)=(ps(1:ifull)-ps_sav(1:ifull))*24.*3600./(100.*dt)

      if(mfix_qg.ne.0.and.mspec==1)then
!       qgmin=1.e-6   !  in parm.h 
        qgminm=qgmin
!          default is moistfix=2 now, with ps weighting
        do k=1,kl
         qg(1:ifull,k)=qg(1:ifull,k)*ps(1:ifull)
         qgsav(1:ifull,k)=qgsav(1:ifull,k)*ps_sav(1:ifull)
        enddo    ! k  loop
        if(ntest==1)then
          call ccglobal_sum(qgsav,sumsav)
          call ccglobal_sum(qg,sumin)
        endif  ! (ntest==1)
!       perform conservation fix on qg, as affected by vadv, hadv, hordif
!       N.B. won't cope with any -ves from conjob
!       delpos is the sum of all positive changes over globe
!       delneg is the sum of all negative changes over globe
        do k=1,kl
         wrk1(1:ifull,k)=max( qg(1:ifull,k),0.,               ! increments  
     &              (qgminm-qfg(1:ifull,k)-qlg(1:ifull,k))*ps(1:ifull) ) 
     &                                       -qgsav(1:ifull,k)   
        enddo                 ! k loop
        if(ntest==2.and.nproc==1)then
          delpos=0.
          delneg=0.
          do k=1,kl
           delposk=0.
           delnegk=0.
           do iq=1,ifull
            delpos=delpos+max(0.,-dsig(k)*wrk1(iq,k)/em(iq)**2)
            delneg=delneg+min(0.,-dsig(k)*wrk1(iq,k)/em(iq)**2)
            delposk=delposk+max(0.,-dsig(k)*wrk1(iq,k)/em(iq)**2)
            delnegk=delnegk+min(0.,-dsig(k)*wrk1(iq,k)/em(iq)**2)
           enddo   ! iq loop
           if(mod(ktau,nmaxpr)==0.or.nmaxpr==1) ! for nproc=1
     &                   print *,'delposk delnegk ',k,delposk,delnegk
          enddo    ! k loop
        else
          call ccglobal_posneg(wrk1,delpos,delneg)  ! usual
        endif
        ratio = -delneg/delpos
        if(mfix_qg==1)alph_q = min(ratio,sqrt(ratio))  ! best option
        if(mfix_qg==2)alph_q = sqrt(ratio)
        do k=1,kl        ! this is cunning 2-sided scheme
         do iq=1,ifull
           qg(iq,k)=qgsav(iq,k)+
     &       alph_q*max(0.,wrk1(iq,k))+min(0.,wrk1(iq,k))/max(1.,alph_q)
         enddo   ! iq loop
        enddo    ! k  loop
        if(ntest==1.or.nmaxpr==1)then  ! has mydiag on prints
         if(mydiag) print *,'qg_delpos,delneg,ratio,alph_q,wrk1 ',
     &                       delpos,delneg,ratio,alph_q,wrk1(idjd,nlv)
          call ccglobal_sum(qg,sumout)
         if(mydiag)print *,'qg_sumsav,sumin,sumout ',sumsav,sumin,sumout
        endif  ! (ntest==1.or.nmaxpr==1)

!       undo ps weighting
        do k=1,kl
         if(mfix_qg>0)then
           qg(1:ifull,k)=qg(1:ifull,k)/ps(1:ifull)
         else
           qg(1:ifull,k)=qg(1:ifull,k)/(ps(1:ifull)*wts(1:ifull))
         endif
        enddo    ! k  loop
        delpos=delpos/1.e5
         delneg=delneg/1.e5  ! for diag print in globpe
      endif        !  (mfix_qg.ne.0.and.mspec==1)

      if(mfix_qg.ne.0.and.mspec==1.and.ldr.ne.0)then
        do k=1,kl
         if(mfix_qg>0)then
          qfg(1:ifull,k)=qfg(1:ifull,k)*ps(1:ifull)
          qfgsav(1:ifull,k)=qfgsav(1:ifull,k)*ps_sav(1:ifull)
          qlg(1:ifull,k)=qlg(1:ifull,k)*ps(1:ifull)
          qlgsav(1:ifull,k)=qlgsav(1:ifull,k)*ps_sav(1:ifull)
         else
          qfg(1:ifull,k)=qfg(1:ifull,k)*ps(1:ifull)*wts(1:ifull)
          qfgsav(1:ifull,k)=
     &            qfgsav(1:ifull,k)*ps_sav(1:ifull)*wts(1:ifull)
          qlg(1:ifull,k)=qlg(1:ifull,k)*ps(1:ifull)*wts(1:ifull)
          qlgsav(1:ifull,k)=
     &            qlgsav(1:ifull,k)*ps_sav(1:ifull)*wts(1:ifull)
         endif
        enddo    ! k  loop
!       perform conservation fix on qfg, qlg as affected by vadv, hadv, hordif
!       N.B. won't cope with any -ves from conjob
!       delpos is the sum of all positive changes over globe
!       delneg is the sum of all negative changes over globe
        do k=1,kl
         do iq=1,ifull
          wrk1(iq,k)=max(qfg(iq,k),0.)-qfgsav(iq,k) ! increments
         enddo   ! iq loop
        enddo    ! k loop
        call ccglobal_posneg(wrk1,delpos,delneg)
        ratio = -delneg/max(delpos,1.e-30)
        if(mfix_qg==1)alph_q = min(ratio,sqrt(ratio))  ! best option
        if(mfix_qg==2)alph_q = sqrt(ratio)
!       this is cunning 2-sided scheme
        qfg(1:ifull,:)=qfgsav(:,:)+
     &     alph_q*max(0.,wrk1(:,:)) + min(0.,wrk1(:,:))/max(1.,alph_q)
        do k=1,kl
         do iq=1,ifull
          wrk1(iq,k)=max(qlg(iq,k),0.)-qlgsav(iq,k)  ! increments
         enddo   ! iq loop
        enddo    ! k loop
        call ccglobal_posneg(wrk1,delpos,delneg)
        ratio = -delneg/max(delpos,1.e-30)
        if(mfix_qg==1)alph_q = min(ratio,sqrt(ratio))  ! best option
        if(mfix_qg==2)alph_q = sqrt(ratio)
!       this is cunning 2-sided scheme
        qlg(1:ifull,:)=qlgsav(:,:)+
     &     alph_q*max(0.,wrk1(:,:)) + min(0.,wrk1(:,:))/max(1.,alph_q)
!          undo ps weighting
        do k=1,kl
         qfg(1:ifull,k)=qfg(1:ifull,k)/ps(1:ifull)
         qlg(1:ifull,k)=qlg(1:ifull,k)/ps(1:ifull)
        enddo    ! k  loop
      endif      !  (mfix_qg.ne0.and.mspec==1.and.ldr.ne.0)

      if(mfix_qg.ne.0.and.mspec==1.and.ngas>=1)then
!       perform conservation fix on tr1,tr2 as affected by vadv, hadv, hordif
        do ng=1,ngas
!          default now with ps weighting
        do k=1,kl
         if(mfix_qg>0)then
           tr(1:ifull,k,ng)=tr(1:ifull,k,ng)*ps(1:ifull)
           trsav(1:ifull,k,ng)=trsav(1:ifull,k,ng)*ps_sav(1:ifull)
         else
           tr(1:ifull,k,ng)=tr(1:ifull,k,ng)*ps(1:ifull)*wts(1:ifull)
           trsav(1:ifull,k,ng)=trsav(1:ifull,k,ng)*ps_sav(1:ifull)
     &                                                  *wts(1:ifull)
         endif
         enddo    ! k  loop
          do k=1,kl
           wrk1(1:ifull,k)=max(tr(1:ifull,k,ng),     ! has increments
     &            gasmin(ng)*ps(1:ifull))         -trsav(1:ifull,k,ng) 
          enddo   ! k loop
         call ccglobal_posneg(wrk1,delpos,delneg)
         ratio = -delneg/delpos
         if(mfix_qg==1)alph_g = min(ratio,sqrt(ratio))  
         if(mfix_qg==2)alph_g = sqrt(ratio)
         do k=1,kl   ! this is cunning 2-sided scheme
          do iq=1,ifull
          tr(iq,k,ng)=trsav(iq,k,ng)+
     &     alph_g*max(0.,wrk1(iq,k)) + min(0.,wrk1(iq,k))/max(1.,alph_g)
          enddo   ! iq loop
         enddo    ! k  loop
         do k=1,kl
          if(mfix_qg>0)then
            tr(1:ifull,k,ng)=tr(1:ifull,k,ng)/ps(1:ifull)
          else
            tr(1:ifull,k,ng)=tr(1:ifull,k,ng)/(ps(1:ifull)*wts(1:ifull))
          endif
         enddo    ! k  loop
        enddo    ! ng loop
        if(mfix_rad>0)then  ! to make gases 2 to ng add up to gas 1
         do k=1,kl
          do iq=1,ifull
           sumdiffb_l = 0.
           delpos_l = 0.
           delneg_l = 0.
           trsav(iq,k,1)=trsav(iq,k,1)/ps(iq)
           do ng=2,ngas        
            trsav(iq,k,ng)=trsav(iq,k,ng)/ps(iq)
            sumdiffb_l = sumdiffb_l + tr(iq,k,ng)
            delpos_l = delpos_l+max( 1.e-20,tr(iq,k,ng)-trsav(iq,k,ng))
            delneg_l = delneg_l+min(-1.e-20,tr(iq,k,ng)-trsav(iq,k,ng))
           enddo   ! ng loop
           delarr_l = (/ delpos_l, delneg_l, sumdiffb_l /)
           call MPI_Allreduce ( delarr_l, delarr, 3, MPI_REAL, MPI_SUM,
     &                          MPI_COMM_WORLD, ierr )
           delpos = delarr(1)
           delneg = delarr(2)
           sumdiffb = delarr(3)
           ratio=(tr(iq,k,1)-sumdiffb)/(delpos-delneg)
           do ng=2,ngas        
            tr(iq,k,ng)=max(0.,trsav(iq,k,ng)
     &         +(1.+ratio)*max(0.,tr(iq,k,ng)-trsav(iq,k,ng))
     &         +(1.-ratio)*min(0.,tr(iq,k,ng)-trsav(iq,k,ng)) )
           enddo   ! ng loop
          enddo   ! iq loop
         enddo    ! k  loop
        endif     ! (mfix_rad>0)
      endif       !  mfix_qg.ne.0

      if ( (diag.or.nmaxpr==1) .and. mydiag ) then
        print *,'at end of adjust5 for ktau= ',ktau
        print *,'ps_sav,ps ',ps_sav(idjd),ps(idjd)
        write (6,"('dpsdt# ',9f8.2)") 
     &            ((dpsdt(ii+jj*il),ii=idjd-1,idjd+1),jj=1,-1,-1)
        write (6,"('qg_a4 ',3p10f8.3)") qg(idjd,:)
        write (6,"('qgs',3p10f8.3)")
     &                         (qgsav(idjd,k)/ps_sav(idjd),k=1,kl)
        write (6,"('qf_a4',3p10f8.3)") qfg(idjd,:)
        write (6,"('ql_a4',3p10f8.3)") qlg(idjd,:)
      endif
      if(diag)then
         call bounds(psl)
         if ( mydiag ) then
            iq=idjd
            print  *,'adjust5 d ',d(iq,:)
!            Would require bounds call
!            print  *,'adjust5 d: e ',d(ie(iq),:)
!            print  *,'adjust5 d: w ',d(iw(iq),:)
            print  *,'adjust5 idjd,e ',idjd, wrk2(idjd,:)
            print  *,'adjust5 psl & n e w s',psl(iq),
     &          psl(in(iq)),psl(ie(iq)),psl(iw(iq)),psl(is(iq))
         end if
         call printa('psl ',psl,ktau,0,ia,ib,ja,jb,0.,100.)
         if ( mydiag ) print *,'adjust5 t ',t(idjd,:)
         call printa('t   ',t,ktau,nlv,ia,ib,ja,jb,200.,1.)
         if ( mydiag ) print *,'adjust5 u ',u(idjd,:)
         call printa('u   ',u,ktau,nlv,ia,ib,ja,jb,0.,1.)
         if ( mydiag ) print *,'adjust5 v ',v(idjd,:)
         call printa('v   ',v,ktau,nlv,ia,ib,ja,jb,0.,1.)
         if ( mydiag ) print *,'ps diagnostic from end of adjust:'
         do iq=1,ifull
            bb(iq)=ps(iq)-ps_sav(iq)
         enddo
         call printa('dps ',bb,ktau,0,ia,ib,ja,jb,0.,.01)
         call printa('ps  ',ps,ktau,0,ia,ib,ja,jb,1.e5,.01)
         if(sig(nlv)<.3)then
            call printa('qg  ',qg,ktau,nlv,ia,ib,ja,jb,0.,1.e6)
         else
            call printa('qg  ',qg,ktau,nlv,ia,ib,ja,jb,0.,1.e3)
         endif
      endif

      dtsave = dt
      call end_log(adjust_end)

      end subroutine adjust5

      subroutine adjust_init
      use cc_mpi
      implicit none
      include 'newmpar.h'
      include 'parm.h'
      include 'parmdyn.h'
      include 'indices.h'
      include 'map.h'
      real zz(ifull),zzn(ifull),zze(ifull),zzw(ifull),
     &     zzs(ifull),pfact(ifull),alff(ifull+iextra),alf(ifull+iextra),
     &     alfe(ifull+iextra),alfn(ifull+iextra),alfu(ifull),alfv(ifull)
      common /zzalf/ zz,zzn,zze,zzw,zzs,pfact,alff,alf,alfe,alfn,
     &     alfu,alfv
      real :: hdt
      integer :: iq, n

      hdt = 0.5*dt
      if(m<6)then
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
      endif  ! (m<6)... else ...
      ! These really only need to be recomputed when time step changes.
      call bounds(alf)
      call bounds(alff,corner=.true.)
!cdir nodep
      do iq=1,ifull
       pfact(iq)=4.*( ds/(dt*em(iq)) )**2    
       alfe(iq)=alf(ie(iq))-alf(iq)    ! Eq. 141 times ds
     &     +.25*(alff(ine(iq))+alfF(in(iq))-alfF(ise(iq))-alfF(is(iq)))
       alfn(iq)=alf(in(iq))-alf(iq)    ! Eq. 142 times ds
     &     -.25*(alfF(ien(iq))+alfF(ie(iq))-alfF(iwn(iq))-alfF(iw(iq)))
      enddo     ! iq loop
      call boundsuv(alfe,alfn)

!cdir nodep
      do iq=1,ifull
!      need care with vector quantities on w (odd) & s (even) panel boundaries
       zz(iq)=.5*(alfe(iwu(iq))-alfe(iq)+alfn(isv(iq))-alfn(iq))
     &                     -4.*alf(iq)                 ! i,j   coeff
       zzn(iq)=alf(in(iq))-.5*alfn(iq)                 ! i,j+1 coeff
       zzw(iq)=alf(iw(iq))+.5*alfe(iwu(iq))            ! i-1,j coeff
       zze(iq)=alf(ie(iq))-.5*alfe(iq)                 ! i+1,j coeff
       zzs(iq)=alf(is(iq))+.5*alfn(isv(iq))            ! i,j-1 coeff
      enddo     ! iq loop
!     N.B. there are some special z values at the 8 vertices
      if(npanels==5)then
         do n=1,npan            ! 0,5
            if ( edge_s .and. edge_w ) then
               iq=indp(1,1,n)
               zzs(iq)=zzs(iq)+.25*alfF(is(iq)) ! i,j-1 coeff
               zzw(iq)=zzw(iq)-.25*alfF(iw(iq)) ! i-1,j coeff
            end if
            if ( edge_n .and. edge_e ) then
               iq=indp(ipan,jpan,n)
               zzn(iq)=zzn(iq)+.25*alfF(in(iq)) ! i,j+1 coeff
               zze(iq)=zze(iq)-.25*alfF(ie(iq)) ! i+1,j coeff
            end if
            if ( edge_s .and. edge_e ) then
               iq=indp(ipan,1,n)
               zzs(iq)=zzs(iq)-.25*alfF(is(iq)) ! i,j-1 coeff
               zze(iq)=zze(iq)+.25*alfF(ie(iq)) ! i+1,j coeff
            end if
            if ( edge_n .and. edge_w ) then
               iq=indp(1,jpan,n)
               zzn(iq)=zzn(iq)-.25*alfF(in(iq)) ! i,j+1 coeff
               zzw(iq)=zzw(iq)+.25*alfF(iw(iq)) ! i-1,j coeff
            end if
        enddo   ! n loop
      endif     ! (npanels==5)
      end subroutine adjust_init
