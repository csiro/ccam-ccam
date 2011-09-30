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
      use tkeeps, only : tke,eps,tkesav,epssav ! MJT tke
      use tracers_m
      use unn_m
      use vecsuv_m
      use vvel_m, omgf => dpsldt
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
      include 'mpif.h'
      real aa(ifull,kl),bb(ifull,kl)
      real p(ifull+iextra,kl),phiv(ifull+iextra,kl),tv(ifull+iextra,kl)
      real ddpds(ifull,kl) ! MJT nh
      integer iq, k, ng, ii, jj, its, nits, nvadh_pass, iaero
      real const_nh, contv, delneg, delpos, ratio
      real sumdiffb, sdmx, sdmx_g, spmax2,termlin
      real, allocatable, save, dimension(:) :: epstsav
      integer :: ierr
      integer, save :: num = 0
      
      call start_log(nonlin_begin)
      
      if (.not.allocated(epstsav)) allocate(epstsav(ifull))

      if(epsp<-2.)then
        if(num==0)epstsav(:)=epst(:)
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
      qgsav(1:ifull,:)=qg(1:ifull,:)      ! for qg  conservation in adjust5
      if(ldr.ne.0)then
        qfgsav(1:ifull,:)=qfg(1:ifull,:)
        qlgsav(1:ifull,:)=qlg(1:ifull,:)
      endif   ! (ldr.ne.0)

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
 
      !--------------------------------------------------------------
      ! MJT tke
      if(nvmix==6)then
        tkesav(1:ifull,:)=tke(1:ifull,:)
        epssav(1:ifull,:)=eps(1:ifull,:)
      endif       ! (nvmix==6)
      !--------------------------------------------------------------
      
      !--------------------------------------------------------------
      ! MJT aerosols
      if (abs(iaero)==2) then
        xtgsav(1:ifull,:,:)=xtg(1:ifull,:,:)
      end if
      !--------------------------------------------------------------

      if (diag) then
         call bounds(ps)
         if ( mydiag ) then
            print *,'at beginning of nonlin'
            print *,'npex,roncp ',npex,roncp
            write (6,"('tn0*dt',9f8.3/6x,9f8.3)") tn(idjd,:)*dt
            write (6,"('un0*dt',9f8.3/6x,9f8.3)") un(idjd,:)*dt
            write (6,"('vn0*dt',9f8.3/6x,9f8.3)") vn(idjd,:)*dt
            write (6,"('tbar',9f8.3/4x,9f8.3)") tbar(:)
            write (6,"('sig ',9f8.5/4x,9f8.5)") sig(:)
            write (6,"('rata',9f8.5/4x,9f8.5)") rata(:)
            write (6,"('ratb',9f8.5/4x,9f8.5)") ratb(:)
            write (6,"('em & news',5f10.4)") em(idjd),
     &          em(in(idjd)),em(ie(idjd)),em(iw(idjd)),em(is(idjd))
            write (6,"('emu,emu_w,emv,emv_s',4f10.4)") 
     &          emu(idjd),emu(iwu(idjd)),emv(idjd),emv(isv(idjd))
            write (6,"('psl & news ',5f9.5)") psl(idjd),
     &          psl(in(idjd)),psl(ie(idjd)),psl(iw(idjd)),psl(is(idjd))
            write (6,"('ps  & news ',-2p5f9.3)") ps(idjd),
     &          ps(in(idjd)),ps(ie(idjd)),ps(iw(idjd)),ps(is(idjd))
         endif
         call printa('u   ',u,ktau,nlv,ia,ib,ja,jb,0.,1.)
         call printa('v   ',v,ktau,nlv,ia,ib,ja,jb,0.,1.)
         call printa('t   ',t,ktau,nlv,ia,ib,ja,jb,200.,1.)
      endif

      if ( (diag.or.nmaxpr==1) .and. mydiag )then
        print *,'in nonlin before possible vertical advection',ktau
        write (6,"('epst#  ',9f8.2)") diagvals(epst) 
        write (6,"('sdot#  ',9f8.3)") diagvals(sdot(:,nlv)) 
        write (6,"('sdotn  ',9f8.3/7x,9f8.3)") sdot(idjd,1:kl)
        write (6,"('omgf#  ',9f8.3)") ((ps(ii+jj*il)*
     &              omgf(ii+jj*il,nlv),ii=idjd-1,idjd+1),jj=-1,1)
        write (6,"('omgfn  ',9f8.3/7x,9f8.3)") ps(idjd)*omgf(idjd,:)
        write (6,"('t   ',9f8.3/4x,9f8.3)")     t(idjd,:)
        write (6,"('u   ',9f8.3/4x,9f8.3)")     u(idjd,:)
        write (6,"('v   ',9f8.3/4x,9f8.3)")     v(idjd,:)
        write (6,"('qg  ',3p9f8.3/4x,9f8.3)")   qg(idjd,:)
      endif

!     do vertical advection in split mode
      if(nvad==4.or.nvad==9)then
        sdmx = maxval(abs(sdot))
        call MPI_AllReduce(sdmx, sdmx_g, 1, MPI_REAL, MPI_MAX,
     &                     MPI_COMM_WORLD, ierr )
        nits=1+sdmx_g/nvadh
        nvadh_pass=nvadh*nits
        if (mydiag.and.mod(ktau,nmaxpr)==0)
     &      print *,'in nonlin sdmx,nits,nvadh_pass ',
     &                         sdmx_g,nits,nvadh_pass
         do its=1,nits
            call vadvtvd(t(1:ifull,:),u(1:ifull,:),v(1:ifull,:),
     &                   nvadh_pass,iaero) 
        enddo
      endif  ! (nvad==4.or.nvad==9)

      if(nvad>=7)then
         call vadv30(t(1:ifull,:),u(1:ifull,:),v(1:ifull,:),iaero)  ! for vadvbess
      endif

cx      do k=1,kl  ! following done in upglobal from 04/09
cx!       N.B. [D + dsigdot/dsig] saved in adjust5 (or updps) as pslx
cx        pslx(1:ifull,k)=psl(1:ifull)-pslx(1:ifull,k)*dt*.5*(1.-epst(:)) !ca
cx      enddo      ! k  loop

      if(nvad>0.and.(diag.or.nmaxpr==1).and.mydiag)then
       print *,'in nonlin after vertical advection'
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
       print *,'pslx ',pslx(idjd,:)
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

      if(nhstest==1) then ! Held and Suarez test case
         call hs_phys
      endif

      if (diag)then
         call printa('sdot',sdot,ktau,nlv+1,ia,ib,ja,jb,0.,10.)
         call printa('omgf',omgf,ktau,nlv,ia,ib,ja,jb,0.,1.e5)
         do iq=1,ifull
            aa(iq,1)=rata(nlv)*sdot(iq,nlv+1)+ratb(nlv)*sdot(iq,nlv)
         enddo
         if ( mydiag )
     &      print *,'k,aa,emu,emv',nlv,aa(idjd,1),emu(idjd),emv(idjd)
         call printa('sgdf',aa(:,1),ktau,nlv,ia,ib,ja,jb,0.,10.)
      endif   ! (diag)

!     extra qfg & qlg terms included in tv from April 04
      tv(1:ifull,:) = (.61*qg(1:ifull,:)-qfg(1:ifull,:)-qlg(1:ifull,:))*
     &                t(1:ifull,:)         ! just add-on at this stage 
      contv=(1.61-cpv/cp)/.61      ! about -.26/.61
      if(ntbar==-1.or.(ntbar==-2.and.num==0))then
        do iq=1,ifull
         tbar2d(iq)=t(iq,1)+contv*tv(iq,1)
        enddo   ! iq loop
      endif     ! (ntbar==-1.or....)
      if(ntbar==0)then
        do iq=1,ifull
         tbar2d(iq)=tbar(1)
        enddo   ! iq loop
      endif     ! (ntbar==0)
      if(ntbar>0)then
        do iq=1,ifull
         tbar2d(iq)=t(iq,ntbar)
        enddo   ! iq loop
      endif     ! (ntbar>0)
      if(ntbar==-3)then
        do iq=1,ifull
         tbar2d(iq)=max(t(iq,1),t(iq,2),t(iq,3),t(iq,kl))
        enddo   ! iq loop
      endif     ! (ntbar==-3)
      if(ntbar==-4)then
        do iq=1,ifull
         tbar2d(iq)=max(t(iq,1),t(iq,2),t(iq,4),t(iq,kl))
        enddo   ! iq loop
      endif     ! (ntbar==-4)
      
      if(ktau==1.or.nh==0)then
        phi(:,1)=zs(1:ifull)+bet(1)*t(1:ifull,1) 
        do k=2,kl
         phi(:,k)=phi(:,k-1)+bet(k)*t(1:ifull,k)+betm(k)*t(1:ifull,k-1)
        enddo    ! k  loop
      endif     ! (ktau==1.or.nh==0)

      if(nh.ne.0)then
        if (abs(epsp).le.1.) then
          const_nh=2.*rdry/(dt*grav*grav*(1.+abs(epsp)))  
        else
          const_nh=2.*rdry/(dt*grav*grav)  
        end if
        do k=1,kl
         h_nh(1:ifull,k)=(1.+epst(:))*tbar(1)*omgf(:,k)/sig(k)
        enddo
        if (nmaxpr==1) then
          if(mydiag)print *,'h_nh.a ',(h_nh(idjd,k),k=1,kl)
        end if
        select case(nh)
         case(3)
          do k=2,kl-1
           ! now includes epst
           h_nh(1:ifull,k)=h_nh(1:ifull,k)
     &       -(sig(k)*(phi(:,k+1)-phi(:,k-1))/(rdry*(sig(k+1)-sig(k-1)))
     &       +t(1:ifull,k))/(const_nh*tbar2d(:))
          enddo
          ! MJT suggestion
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
          ! MJT suggestion
          k=1
          h_nh(1:ifull,k)=h_nh(1:ifull,k)
     &    -(((sig(k)-1.)*(phi(:,k+1)-phi(:,k))/(sig(k+1)-sig(k))+
     &    ((sig(k+1)-sig(k))*(phi(:,k)-zs(1:ifull))/(sig(k)-1.)))
     &     *sig(k)/(rdry*(sig(k+1)-1.))  
     &      +t(1:ifull,k))/(const_nh*tbar2d(:))
         case(5)
          ! MJT - This method is compatible with bet(k) and betm(k)
          ! For hydrostatic case, ddpsds exactly cancels with t.
          ! This is the same as nh==2, but works for all lapsbot
          ! ddpds is (sig/rdry)*d(phi)/d(sig)
          ddpds(:,1)=-(phi(:,1)-zs(1:ifull))/bet(1)
          do k=2,kl
            ddpds(:,k)=-(phi(:,k)-phi(:,k-1))/bet(k)
     &        -betm(k)*ddpds(:,k-1)/bet(k)
          end do
          do k=1,kl
            h_nh(1:ifull,k)=h_nh(1:ifull,k)
     &        -(ddpds(:,k)+t(1:ifull,k))/(const_nh*tbar2d(:))
          end do
        end select
        if (nmaxpr==1) then
          if (mydiag)then
            print *,'h_nh.b ',(h_nh(idjd,k),k=1,kl)
            print *,'phi ',(phi(idjd,k),k=1,kl)
          endif
          call maxmin(h_nh,'h_',ktau,1.,kl)
        endif
      endif      ! (nh.ne.0)

      do k=1,kl
       do iq=1,ifull
        termlin=tbar2d(iq)*omgf(iq,k)*roncp/sig(k) ! full omgf used here
        tn(iq,k)=tn(iq,k)+(t(iq,k)+contv*tv(iq,k)-tbar2d(iq))
     &           *omgf(iq,k)*roncp/sig(k) 
!       add in  cnon*dt*tn(iq,k)  term at bottom
        tx(iq,k)=t(iq,k) +.5*dt*(1.-epst(iq))*termlin  
cy      tx(iq,k)=.5*dt*termlin  ! t and epst later  cy
       enddo     ! iq loop
      enddo      ! k  loop
      if( (diag.or.nmaxpr==1) .and. mydiag )then
        iq=idjd
        k=nlv
        print *,'contv,tbar2d,termlin_nlv ',
     &           contv,tbar2d(iq),tbar2d(iq)*omgf(iq,k)*roncp/sig(k)
        print *,'tv,tn ',tv(iq,k),tn(iq,k)
c       print *,'termx ',(t(iq,k)+contv*tvv)*omgf(iq,k)*roncp/sig(k)
      endif
               
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

      call bounds(p)
      call bounds(phiv)
      call bounds(tv)
      call bounds(psl)

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
        if(npex.ne.6)then
          aa(1:ifull,:)=aa(1:ifull,:)+.5*dt*un(1:ifull,:) ! still staggered
          bb(1:ifull,:)=bb(1:ifull,:)+.5*dt*vn(1:ifull,:) ! still staggered
        endif  ! (npex.ne.6)
        if(diag)then
          if(mydiag)then
            print *,'tv ',tv(idjd,:)
            write (6,"('tn1*dt',9f8.3/6x,9f8.3)") tn(idjd,:)*dt
            write (6,"('un1*dt',9f8.3/6x,9f8.3)") un(idjd,:)*dt
            write (6,"('vn1*dt',9f8.3/6x,9f8.3)") vn(idjd,:)*dt
          endif
        endif                     ! (diag)
      endif  ! (npex==5 ... else ...)
      call unstaguv(aa,bb,ux,vx) ! convert to unstaggered positions
      if(diag)then
        call printa('aa  ',aa,ktau,nlv,ia,ib,ja,jb,0.,1.)
        call printa('bb  ',bb,ktau,nlv,ia,ib,ja,jb,0.,1.)
      endif                     ! (diag)
      ux(1:ifull,:)=u(1:ifull,:)+ux(1:ifull,:)
      vx(1:ifull,:)=v(1:ifull,:)+vx(1:ifull,:)
      if(npex==5)then
        ux(1:ifull,:)=ux(1:ifull,:)+.5*dt*un(1:ifull,:) ! unstaggered
        vx(1:ifull,:)=vx(1:ifull,:)+.5*dt*vn(1:ifull,:) ! unstaggered
      endif
      if(npex<3)call unstaguv(un,vn,un,vn) 
      
      if(nxmap==2)then  ! not ready with npex=3
        do k=1,kl
         aa(1:ifull,k)=dmdy(1:ifull)*u(1:ifull,k)   ! actually fm
     &                -dmdx(1:ifull)*v(1:ifull,k) 
         un(1:ifull,k)=un(1:ifull,k)+aa(1:ifull,k)*v(1:ifull,k) 
         vn(1:ifull,k)=vn(1:ifull,k)-aa(1:ifull,k)*u(1:ifull,k)   
        enddo
        if (diag.and.mydiag)then
          print *,'fm ',aa(idjd,:)
          write (6,"('fm#  ',4p9f8.2)") diagvals(aa(:,nlv)) 
        endif
      endif  ! (nxmap==2)

      tx(1:ifull,:) = tx(1:ifull,:) + .5*dt*tn(1:ifull,:)
      if(nvsplit==3.or.nvsplit==4)then
        ux(1:ifull,:)=ux(1:ifull,:)+.5*dt*unn(1:ifull,:) ! phys contrib
        vx(1:ifull,:)=vx(1:ifull,:)+.5*dt*vnn(1:ifull,:) ! phys contrib
      endif
!     N.B. don't add npex==0,1,2 (or 3) dyn contribs here, as already added

      if (diag)then
         if(mydiag) then
           print *,'at end of nonlin; nvad,idjd = ', nvad,idjd
           print *,'p1 . & e ',p(idjd,nlv),p(ie(idjd),nlv)
           print *,'p1 . & n ',p(idjd,nlv),p(in(idjd),nlv)
           print *,'tx ',tx(idjd,:)
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

      num=1
      call end_log(nonlin_end)
      return
      end
