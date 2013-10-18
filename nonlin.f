      subroutine nonlin(iaero)
      use aerosolldr      
      use arrays_m
      use cc_mpi
      use diag_m
      use epst_m
      use indices_m
      use latlong_m
      use liqwpar_m  ! qfg,qlg
      use map_m
      use morepbl_m  ! condx
      use nharrs_m
      use nlin_m
      use savuvt_m
      use sigs_m
      use tbar2d_m
      use tracers_m
      use unn_m
      use vadv
      use vecsuv_m
      use vvel_m
      use work3sav_m
      use xarrs_m
      use xyzinfo_m
      implicit none
      integer, parameter :: ntest=0
      integer, parameter :: mfix_rad=0 ! used to make gases 2 to ng add up to gas 1
      include 'newmpar.h'
      include 'const_phys.h' ! r,g,cp,cpv,roncp
      include 'kuocom.h'   ! ldr
      include 'parm.h'
      include 'parmdyn.h'  
      include 'parmvert.h'
      integer iq, k, ng, ii, jj, iaero
      integer ierr
      integer, save :: num = 0
      integer, dimension(ifull) :: nits, nvadh_pass
      real aa(ifull,kl),bb(ifull,kl)
      real p(ifull+iextra,kl),phiv(ifull+iextra,kl),tv(ifull+iextra,kl)
      real ddpds(ifull)
      real duma(ifull+iextra,2*kl)
      real const_nh, contv, delneg, delpos, ratio
      real sumdiffb, spmax2,termlin
      real, dimension(ifull) :: sdmx
      real, allocatable, save, dimension(:) :: epstsav
      real, dimension(ifull,kl) :: dumt,dumu,dumv
      
      call start_log(nonlin_begin)
     
      if(epsp<-2.)then
        if (.not.allocated(epstsav)) then
           allocate(epstsav(ifull))
           epstsav(:)=epst(:)
        end if
        do iq=1,ifull
         spmax2=max(u(iq,3*kl/4)**2+v(iq,3*kl/4)**2,
     &              u(iq,  kl  )**2+v(iq,  kl  )**2) 
         if(spmax2>(.8*ds/(em(iq)*dt))**2)then
!          setting epst for Courant number > .8        
           epst(iq)=epstsav(iq)
         else
           epst(iq)=0.
         endif
        enddo         
      endif  ! (epsp<-2.)

!     *** following qgsav should be before first vadv call
      do k=1,kl
        qgsav(:,k)=qg(1:ifull,k)      ! for qg  conservation in adjust5
      end do
      if(ldr/=0)then
        do k=1,kl
          qfgsav(:,k)=qfg(1:ifull,k)
          qlgsav(:,k)=qlg(1:ifull,k)
          qrgsav(:,k)=qrg(1:ifull,k)
        end do
      endif   ! (ldr.ne.0)
      if (abs(iaero)==2) then
        xtgsav(:,1:kl,1:naero)=xtg(1:ifull,1:kl,1:naero)
      end if

      if(ngas>=1)then
        if(mfix_rad>0.and.mspec==1)then ! make gases 2 to ng add up to g1
         do k=1,kl              ! here it is just from effects of physics
          do iq=1,ifull
          sumdiffb=0.
           do ng=2,ngas        
            sumdiffb=sumdiffb+tr(iq,k,ng)
            delpos=delpos+max( 1.e-20,tr(iq,k,ng)-trsav(iq,k,ng))
            delneg=delneg+min(-1.e-20,tr(iq,k,ng)-trsav(iq,k,ng))
           enddo   ! ng loop
           ratio=(tr(iq,k,1)-sumdiffb)/(delpos-delneg)
           do ng=2,ngas        
            tr(iq,k,ng)=max(0.,trsav(iq,k,ng)
     &         +(1.+ratio)*max(0.,tr(iq,k,ng)-trsav(iq,k,ng))
     &         +(1.-ratio)*min(0.,tr(iq,k,ng)-trsav(iq,k,ng)) )
           enddo   ! ng loop
          enddo    ! iq loop
         enddo     ! k  loop
        endif      ! (mfix_rad>0)
        do ng=1,ngas   ! re-set trsav prior to vadv, hadv, hordif
         do k=1,kl               
          do iq=1,ifull
           trsav(iq,k,ng)=tr(iq,k,ng) ! for tr conservation in adjust5
          enddo   ! iq loop
         enddo    ! k  loop
        enddo     ! ng loop
      endif       ! (ngas>=1)
 
#ifdef debug
      if ( diag.or.nmaxpr==1 ) then
       call bounds(ps)
       if ( mydiag ) then
          write(6,*) "qgsav ",qgsav(idjd,nlv)
          write(6,*) 'at beginning of nonlin'
          write(6,*) 'npex,roncp ',npex,roncp
          write (6,"('tn0*dt',9f8.3/6x,9f8.3)") tn(idjd,:)*dt
          write (6,"('un0*dt',9f8.3/6x,9f8.3)") un(idjd,:)*dt
          write (6,"('vn0*dt',9f8.3/6x,9f8.3)") vn(idjd,:)*dt
          write (6,"('tbar',9f8.3/4x,9f8.3)") tbar(:)
          write (6,"('sig ',9f8.5/4x,9f8.5)") sig(:)
          write (6,"('rata',9f8.5/4x,9f8.5)") rata(:)
          write (6,"('ratb',9f8.5/4x,9f8.5)") ratb(:)
          write (6,"('em & news',5f10.4)") em(idjd),
     &        em(in(idjd)),em(ie(idjd)),em(iw(idjd)),em(is(idjd))
          write (6,"('emu,emu_w,emv,emv_s',4f10.4)") 
     &        emu(idjd),emu(iwu(idjd)),emv(idjd),emv(isv(idjd))
          write (6,"('psl & news ',5f9.5)") psl(idjd),
     &        psl(in(idjd)),psl(ie(idjd)),psl(iw(idjd)),psl(is(idjd))
          write (6,"('ps  & news ',-2p5f9.3)") ps(idjd),
     &        ps(in(idjd)),ps(ie(idjd)),ps(iw(idjd)),ps(is(idjd))
       endif
       call printa('u   ',u,ktau,nlv,ia,ib,ja,jb,0.,1.)
       call printa('v   ',v,ktau,nlv,ia,ib,ja,jb,0.,1.)
       call printa('t   ',t,ktau,nlv,ia,ib,ja,jb,200.,1.)

       if ( mydiag )then
        write(6,*) 'in nonlin before possible vertical advection',ktau
        write (6,"('epst#  ',9f8.2)") diagvals(epst) 
        write (6,"('sdot#  ',9f8.3)") diagvals(sdot(:,nlv)) 
        write (6,"('sdotn  ',9f8.3/7x,9f8.3)") sdot(idjd,1:kl)
        write (6,"('omgf#  ',9f8.3)") ps(idjd)*dpsldt(idjd,nlv)
        write (6,"('omgfn  ',9f8.3/7x,9f8.3)") ps(idjd)*dpsldt(idjd,:)
        write (6,"('t   ',9f8.3/4x,9f8.3)")     t(idjd,:)
        write (6,"('u   ',9f8.3/4x,9f8.3)")     u(idjd,:)
        write (6,"('v   ',9f8.3/4x,9f8.3)")     v(idjd,:)
        write (6,"('qg  ',3p9f8.3/4x,9f8.3)")   qg(idjd,:)
       end if
      endif
#endif

!     do vertical advection in split mode
      if(nvad==4.or.nvad==9)then
        sdmx(:) = maxval(abs(sdot),2)
        nits(:)=1+sdmx(:)/nvadh
        nvadh_pass(:)=nvadh*nits(:)
#ifdef debug
        if (mydiag.and.mod(ktau,nmaxpr)==0)
     &   write(6,*) 'in nonlin sdmx,nits,nvadh_pass ',
     &             sdmx(idjd),nits(idjd),nvadh_pass(idjd)
#endif
          call vadvtvd(t,u,v,nvadh_pass,nits,iaero)
      endif  ! (nvad==4.or.nvad==9)

      if(nvad>=7)then
         call vadv30(t(1:ifull,:),u(1:ifull,:),v(1:ifull,:),iaero)  ! for vadvbess
      endif

#ifdef debug
      if(nvad>0.and.(diag.or.nmaxpr==1).and.mydiag)then
       write(6,*) 'in nonlin after vertical advection'
       write (6,"('t   ',9f8.2/4x,9f8.2)") t(idjd,:)
       write (6,"('thet',9f8.2/4x,9f8.2)") t(idjd,:)*sig(:)**(-roncp)
       write (6,"('qg  ',9f8.3/4x,9f8.3)") 1000.*qg(idjd,:)
       write (6,"('ql ',9f8.3/4x,9f8.3)") 1000.*qlg(idjd,:)
       write (6,"('qf ',9f8.3/4x,9f8.3)") 1000.*qfg(idjd,:)
       write (6,"('u   ',9f8.2/4x,9f8.2)") u(idjd,:)
       write (6,"('v   ',9f8.2/4x,9f8.2)") v(idjd,:)
       write (6,"('t#  ',9f8.2)") diagvals(t(:,nlv)) 
!    &           ((t(ii+jj*il,nlv),ii=idjd-1,idjd+1),jj=-1,1)
       write (6,"('qg# ',3p9f8.3)") diagvals(qg(:,nlv)) 
       write (6,"('ql# ',3p9f8.3)") diagvals(qlg(:,nlv)) 
       write (6,"('qf# ',3p9f8.3)") diagvals(qfg(:,nlv)) 
       write (6,"('u#  ',9f8.2)") diagvals(u(:,nlv)) 
       write (6,"('v#  ',9f8.2)") diagvals(v(:,nlv)) 
       write(6,*) 'pslx ',pslx(idjd,:)
      endif  ! (nvad>0.and.(diag.or.nmaxpr==1).and.mydiag)
      if (diag)then
         if ( mydiag ) write (6,"('qg ',3p9f8.5/4x,9f8.5)") qg(idjd,:)
         if (sig(nlv)<.3)then
            call printa('qg  ',qg,ktau,nlv,ia,ib,ja,jb,0.,1.e6)
         else
            call printa('qg  ',qg,ktau,nlv,ia,ib,ja,jb,0.,1.e3)
         endif
         if ( mydiag ) then
            write (6,"('u1 & ew ',5f8.2)") u(idjd,1),
     &                  u(ieu(idjd),1),u(iwu(idjd),1)
            write (6,"('v1 & ns ',5f8.2)") v(idjd,1),
     &                  v(inv(idjd),1),v(isv(idjd),1)
            write (6,"('u9 & ew ',5f8.2)") u(idjd,9),
     &                  u(ieu(idjd),9),u(iwu(idjd),9)
            write (6,"('v9 & ns ',5f8.2)") v(idjd,9),
     &                  v(inv(idjd),9),v(isv(idjd),9)
         endif
      endif
#endif

      if(nhstest==1) then ! Held and Suarez test case
         call hs_phys
      endif

#ifdef debug
      if (diag)then
         call printa('sdot',sdot,ktau,nlv+1,ia,ib,ja,jb,0.,10.)
         call printa('omgf',dpsldt,ktau,nlv,ia,ib,ja,jb,0.,1.e5)
         do iq=1,ifull
            aa(iq,1)=rata(nlv)*sdot(iq,nlv+1)+ratb(nlv)*sdot(iq,nlv)
         enddo
         if ( mydiag )
     &    write(6,*) 'k,aa,emu,emv',nlv,aa(idjd,1),emu(idjd),emv(idjd)
         call printa('sgdf',aa(:,1),ktau,nlv,ia,ib,ja,jb,0.,10.)
      endif   ! (diag)
#endif

!     extra qfg & qlg terms included in tv from April 04
      tv(1:ifull,:) = (.61*qg(1:ifull,:)-qfg(1:ifull,:)-qlg(1:ifull,:))*
     &                t(1:ifull,:)         ! just add-on at this stage 
      contv=(1.61-cpv/cp)/.61      ! about -.26/.61
      if(ntbar==-1.or.(ntbar==-2.and.num==0))then
        do iq=1,ifull
         tbar2d(iq)=t(iq,1)+contv*tv(iq,1)
        enddo   ! iq loop
      else if (ntbar==0)then
        do iq=1,ifull
         tbar2d(iq)=tbar(1)
        enddo   ! iq loop
      else if (ntbar>0)then
        do iq=1,ifull
         tbar2d(iq)=t(iq,ntbar)
        enddo   ! iq loop
      else if (ntbar==-3)then
        do iq=1,ifull
         tbar2d(iq)=max(t(iq,1),t(iq,2),t(iq,3),t(iq,kl))
        enddo   ! iq loop
      else if (ntbar==-4)then
        do iq=1,ifull
         tbar2d(iq)=max(t(iq,1),t(iq,2),t(iq,4),t(iq,kl))
        enddo   ! iq loop
      endif     ! (ntbar==-4)
      
      ! update hydrostatic phi and add non-hydrostatic component
      phi(:,1)=zs(1:ifull)+bet(1)*t(1:ifull,1) 
      do k=2,kl
       phi(:,k)=phi(:,k-1)+bet(k)*t(1:ifull,k)+betm(k)*t(1:ifull,k-1)
      enddo    ! k  loop
      ! update non-hydrostatic terms from Miller-White height equation
      if (nh/=0) then
        phi=phi+phi_nh
        if (abs(epsp)<=1.) then
          ! exact treatment of constant epsp terms
          const_nh=2.*rdry/(dt*grav*grav*(1.-epsp*epsp))
        else
          const_nh=2.*rdry/(dt*grav*grav)  
        end if
        do k=1,kl
         h_nh(1:ifull,k)=(1.+epst(:))*tbar(1)*dpsldt(:,k)/sig(k)
        enddo
#ifdef debug
        if (nmaxpr==1) then
          if(mydiag) write(6,*) 'h_nh.a ',(h_nh(idjd,k),k=1,kl)
        end if
#endif
        select case(nh)
         case(3)
          do k=2,kl-1
           ! now includes epst
           h_nh(1:ifull,k)=h_nh(1:ifull,k)
     &       -(sig(k)*(phi(:,k+1)-phi(:,k-1))/(rdry*(sig(k+1)-sig(k-1)))
     &       +t(1:ifull,k))/(const_nh*tbar2d(:))
          enddo
          k=1
          h_nh(1:ifull,k)=h_nh(1:ifull,k)
     &      -(sig(k)*(phi(:,k+1)-zs(1:ifull))/(rdry*(sig(k+1)-1.))
     &      +t(1:ifull,k))/(const_nh*tbar2d(:))
          k=kl
          h_nh(1:ifull,k)=h_nh(1:ifull,k)
     &      -(sig(k)*(phi(:,k)-phi(:,k-1))/(rdry*(sig(k)-sig(k-1)))
     &      +t(1:ifull,k))/(const_nh*tbar2d(:))
         case(2) ! was -2 add in other term explicitly, more consistently
!         N.B. nh=2 needs lapsbot=3        
          do k=2,kl
           h_nh(1:ifull,k)=h_nh(1:ifull,k)
     &       -((phi(:,k)-phi(:,k-1))/bet(k)          ! using bet
     &       +t(1:ifull,k))/(const_nh*tbar2d(:))
          enddo
          k=1
           h_nh(1:ifull,k)=h_nh(1:ifull,k)
     &       -((phi(:,k)-zs(1:ifull))/bet(k)         ! using bet
     &       +t(1:ifull,k))/(const_nh*tbar2d(:))
         case(4) ! was -3 add in other term explicitly, more accurately?
          do k=2,kl-1
           h_nh(1:ifull,k)=h_nh(1:ifull,k)
     &     -(((sig(k)-sig(k-1))*(phi(:,k+1)-phi(:,k))/(sig(k+1)-sig(k))+
     &     ((sig(k+1)-sig(k))*(phi(:,k)-phi(:,k-1))/(sig(k)-sig(k-1))))
     &      *sig(k)/(rdry*(sig(k+1)-sig(k-1)))  
     &       +t(1:ifull,k))/(const_nh*tbar2d(:))
          enddo
          k=1
          h_nh(1:ifull,k)=h_nh(1:ifull,k)
     &    -(((sig(k)-1.)*(phi(:,k+1)-phi(:,k))/(sig(k+1)-sig(k))+
     &    ((sig(k+1)-sig(k))*(phi(:,k)-zs(1:ifull))/(sig(k)-1.)))
     &     *sig(k)/(rdry*(sig(k+1)-1.))  
     &      +t(1:ifull,k))/(const_nh*tbar2d(:))
         case(5)
          ! MJT - This method is compatible with bet(k) and betm(k)
          ! This is the similar to nh==2, but works for all lapsbot
          ! and only involves phi_nh as the hydrostatic component
          ! is eliminated.
          ! ddpds is -(sig/rdry)*d(phi_nh)/d(sig) or delta T_nh
          ddpds=phi_nh(:,1)/bet(1)
          h_nh(1:ifull,1)=h_nh(1:ifull,1)
     &      +ddpds/(const_nh*tbar2d(:))
          do k=2,kl
            ddpds=(phi_nh(:,k)-phi_nh(:,k-1)
     &        -betm(k)*ddpds)/bet(k)
            h_nh(1:ifull,k)=h_nh(1:ifull,k)
     &        +ddpds/(const_nh*tbar2d(:))
          end do
        end select
#ifdef debug
        if (nmaxpr==1) then
          if (mydiag)then
           write(6,*) 'h_nh.b ',(h_nh(idjd,k),k=1,kl)
           write(6,*) 'phi ',(phi(idjd,k),k=1,kl)
           write(6,*) 'phi_nh ',(phi_nh(idjd,k),k=1,kl)
          endif
          call maxmin(h_nh,'h_',ktau,1.,kl)
        endif
#endif
      else
        phi_nh=0. ! set to hydrostatic approximation
        h_nh=0.
      endif      ! (nh.ne.0.and.(ktau.gt.knh.or.lrestart)) ..else..

      do k=1,kl
       do iq=1,ifull
        termlin=tbar2d(iq)*dpsldt(iq,k)*roncp/sig(k) ! full dpsldt used here
        tn(iq,k)=tn(iq,k)+(t(iq,k)+contv*tv(iq,k)-tbar2d(iq))
     &           *dpsldt(iq,k)*roncp/sig(k) 
!       add in  cnon*dt*tn(iq,k)  term at bottom
        tx(iq,k)=t(iq,k) +.5*dt*(1.-epst(iq))*termlin  
cy      tx(iq,k)=.5*dt*termlin  ! t and epst later  cy
       enddo     ! iq loop
      enddo      ! k  loop

#ifdef debug      
      if( (diag.or.nmaxpr==1) .and. mydiag )then
        iq=idjd
        k=nlv
        write(6,*) 'dpsldt,roncp,sig ',
     &           dpsldt(iq,k),roncp,sig(k)
        write(6,*) 'contv,tbar2d,termlin_nlv ',
     &           contv,tbar2d(iq),tbar2d(iq)*dpsldt(iq,k)*roncp/sig(k)
        write(6,*) 'tv,tn ',tv(iq,k),tn(iq,k)
c       write(6,*) 'termx ',(t(iq,k)+contv*tvv)*dpsldt(iq,k)*roncp/sig(k)
      endif
#endif
               
!     calculate augmented geopotential height terms and save in p
      do k=1,kl
       p(1:ifull,k)=phi(1:ifull,k)+rdry*tbar2d(1:ifull)*psl(1:ifull)
      enddo      ! k  loop

!     calculate virtual temp extra terms (mainly -ve): phi -phi_v +r*tbar*psl
      phiv(1:ifull,1)=rdry*tbar2d(1:ifull)*psl(1:ifull)
     &                -bet(1)*tv(1:ifull,1)
      do k=2,kl
       phiv(1:ifull,k)=phiv(1:ifull,k-1)
     &                 -bet(k)*tv(1:ifull,k)-betm(k)*tv(1:ifull,k-1)
      enddo    ! k  loop

!     also need full Tv
      do k=1,kl
       tv(1:ifull,k)=t(1:ifull,k)+tv(1:ifull,k)  
      enddo


      ! MJT notes - This is the first bounds call after
      ! the physics routines, so load balance is a
      ! significant issue.
      duma(1:ifull,1:kl)     =p(1:ifull,:)
      duma(1:ifull,kl+1:2*kl)=tv(1:ifull,:)
      call bounds(duma(:,1:2*kl),nehalf=.true.)
      p(ifull+1:ifull+iextra,:) =duma(ifull+1:ifull+iextra,1:kl)
      tv(ifull+1:ifull+iextra,:)=duma(ifull+1:ifull+iextra,kl+1:2*kl)
      duma(1:ifull,1:kl)=phiv(1:ifull,:)
      duma(1:ifull,kl+1)=psl(1:ifull)
      call bounds(duma(:,1:kl+1))
      phiv(ifull+1:ifull+iextra,1:kl)=duma(ifull+1:ifull+iextra,1:kl)
      psl(ifull+1:ifull+iextra)      =duma(ifull+1:ifull+iextra,kl+1)


      do k=1,kl
!cdir nodep
       do iq=1,ifull          ! calculate staggered ux,vx first
        aa(iq,k)=-.5*dt*emu(iq)*(p(ie(iq),k)-p(iq,k))*(1.-epsu)/ds
        bb(iq,k)=-.5*dt*emv(iq)*(p(in(iq),k)-p(iq,k))*(1.-epsu)/ds
       enddo   ! iq loop
      enddo    ! k loop

      if(npex==5)then   ! rather noisy with nstag=-1 (incl. old 5)
!       npex=5 same as npex=0, but does direct unstaggered calc      
        do k=1,kl
!cdir nodep
         do iq=1,ifull  ! calculate unstaggered dyn residual contributions
          un(iq,k)=.5*em(iq)*(phiv(ie(iq),k)-phiv(iw(iq),k)
     &             -rdry*tv(iq,k)*(psl(ie(iq))-psl(iw(iq))))/ds
          vn(iq,k)=.5*em(iq)*(phiv(in(iq),k)-phiv(is(iq),k)
     &             -rdry*tv(iq,k)*(psl(in(iq))-psl(is(iq))))/ds
         enddo   ! iq loop
        enddo    ! k loop
      else       ! i.e. npex.ne.5
        do k=1,kl
!cdir nodep
         do iq=1,ifull  ! calculate staggered dyn residual contributions first
          un(iq,k)=emu(iq)*(phiv(ie(iq),k)-phiv(iq,k)
     &        -.5*rdry*(tv(ie(iq),k)+tv(iq,k))*(psl(ie(iq))-psl(iq)))/ds
          vn(iq,k)=emv(iq)*(phiv(in(iq),k)-phiv(iq,k)
     &        -.5*rdry*(tv(in(iq),k)+tv(iq,k))*(psl(in(iq))-psl(iq)))/ds
         enddo   ! iq loop
        enddo    ! k loop
        if(npex/=6)then
          aa(1:ifull,:)=aa(1:ifull,:)+.5*dt*un(1:ifull,:) ! still staggered
          bb(1:ifull,:)=bb(1:ifull,:)+.5*dt*vn(1:ifull,:) ! still staggered
        endif  ! (npex.ne.6)
#ifdef debug
        if(diag)then
          if(mydiag)then
            write(6,*) 'tv ',tv(idjd,:)
            write (6,"('tn1*dt',9f8.3/6x,9f8.3)") tn(idjd,:)*dt
            write (6,"('un1*dt',9f8.3/6x,9f8.3)") un(idjd,:)*dt
            write (6,"('vn1*dt',9f8.3/6x,9f8.3)") vn(idjd,:)*dt
          endif
        endif                     ! (diag)
#endif
      endif  ! (npex==5 ... else ...)

      call unstaguv(aa,bb,ux,vx) ! convert to unstaggered positions

#ifdef debug
      if(diag)then
        call printa('aa  ',aa,ktau,nlv,ia,ib,ja,jb,0.,1.)
        call printa('bb  ',bb,ktau,nlv,ia,ib,ja,jb,0.,1.)
      endif                     ! (diag)
#endif

      ux(1:ifull,:)=u(1:ifull,:)+ux(1:ifull,:)
      vx(1:ifull,:)=v(1:ifull,:)+vx(1:ifull,:)
      
      if(npex==5)then
        ux(1:ifull,:)=ux(1:ifull,:)+.5*dt*un(1:ifull,:) ! unstaggered
        vx(1:ifull,:)=vx(1:ifull,:)+.5*dt*vn(1:ifull,:) ! unstaggered
      else if (npex<3) then
        call unstaguv(un,vn,un,vn) 
      end if
      
      if(nxmap==2)then  ! not ready with npex=3
        do k=1,kl
         aa(1:ifull,k)=dmdy(1:ifull)*u(1:ifull,k)   ! actually fm
     &                -dmdx(1:ifull)*v(1:ifull,k) 
         un(1:ifull,k)=un(1:ifull,k)+aa(1:ifull,k)*v(1:ifull,k) 
         vn(1:ifull,k)=vn(1:ifull,k)-aa(1:ifull,k)*u(1:ifull,k)   
        enddo
#ifdef debug
        if (diag.and.mydiag)then
          write(6,*) 'fm ',aa(idjd,:)
          write (6,"('fm#  ',4p9f8.2)") diagvals(aa(:,nlv)) 
        endif
#endif
      endif  ! (nxmap==2)

      tx(1:ifull,:) = tx(1:ifull,:) + .5*dt*tn(1:ifull,:)
      if(nvsplit==3.or.nvsplit==4)then
        ux(1:ifull,:)=ux(1:ifull,:)+.5*dt*unn(1:ifull,:) ! phys contrib
        vx(1:ifull,:)=vx(1:ifull,:)+.5*dt*vnn(1:ifull,:) ! phys contrib
      endif
!     N.B. don't add npex==0,1,2 (or 3) dyn contribs here, as already added

#ifdef debug
      if (diag)then
         if(mydiag) then
           write(6,*) 'at end of nonlin; nvad,idjd = ', nvad,idjd
           write(6,*) 'p1 . & e ',p(idjd,nlv),p(ie(idjd),nlv)
           write(6,*) 'p1 . & n ',p(idjd,nlv),p(in(idjd),nlv)
           write(6,*) 'tx ',tx(idjd,:)
           write (6,"('tn2*dt',9f8.3/6x,9f8.3)")   tn(idjd,:)*dt
           write (6,"('un2*dt',9f8.3/6x,9f8.3)")   un(idjd,:)*dt
           write (6,"('vn2*dt',9f8.3/6x,9f8.3)")   vn(idjd,:)*dt
           write (6,"('ux  ',9f8.2/4x,9f8.2)")     ux(idjd,:)
           write (6,"('vx  ',9f8.2/4x,9f8.2)")     vx(idjd,:)
         endif
         call printa('psl ',psl,ktau,0,ia,ib,ja,jb,0.,100.)
         call printa('pslx',pslx,ktau,nlv,ia,ib,ja,jb,0.,100.)
         call printa('tn  ',tn,ktau,nlv,ia,ib,ja,jb,0.,100.*dt)
         call printa('un  ',un,ktau,nlv,ia,ib,ja,jb,0.,100.*dt)
         call printa('vn  ',vn,ktau,nlv,ia,ib,ja,jb,0.,100.*dt)
         call printa('tx  ',tx,ktau,nlv,ia,ib,ja,jb,200.,1.)
         call printa('ux  ',ux,ktau,nlv,ia,ib,ja,jb,0.,1.)
         call printa('vx  ',vx,ktau,nlv,ia,ib,ja,jb,0.,1.)
      endif
#endif

      num=1

      call end_log(nonlin_end)
      
      return
      end
