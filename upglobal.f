      subroutine upglobal(iaero)      ! globpea version   use ritchie 103
      use aerosolldr
      use arrays_m
      use cc_mpi
      use cfrac_m
      use diag_m
      use epst_m
      use indices_m
      use liqwpar_m  ! ifullw
      use map_m
      use nharrs_m
      use nlin_m
      use sbar_m
      use sigs_m
      use tkeeps, only : tke,eps
      use tracers_m
      use unn_m
      use vecsuv_m
      use vvel_m     ! sdot
      use work3f_m
      use xarrs_m
      use xyzinfo_m
      implicit none
!     parameter (nrot=1)       ! nrot=1 to rotate velocity vectors (parmdyn.h)
      integer, parameter :: ntest=0       ! ~8+ for diagnostic stability tests
!     m=6 had three options; only option 3 has been adopted here
      include 'newmpar.h'
      include 'const_phys.h'
      include 'kuocom.h'   ! ldr
      include 'parm.h'
      include 'parmdyn.h'  
      include 'parmhor.h'  ! mhint, m_bs, nt_adv
      include 'parmvert.h'  
      real, save, allocatable, dimension(:,:):: tnsav,unsav,vnsav ! for npex=-1
      real, dimension(ifull+iextra,kl) :: uc, vc, wc, dd
      real aa(ifull+iextra)
      real*8 x3d(ifull,kl),y3d(ifull,kl),z3d(ifull,kl)
      integer idjdd
      real theta(ifull,kl), factr(kl)
      real, dimension(ifull,kl) :: dumt,dumu,dumv
      integer ii,intsch, iq, jj,k, kk, ntr, ierr
      integer iaero, l
      integer, dimension(ifull) :: nits, nvadh_pass
      real denb, tempry, vdot1,
     &     vdot2, vec1x, vec1y, vec1z, vec2x, vec2y, vec2z, vec3x,
     &     vec3y, vec3z, vecdot
      real, dimension(ifull) :: sdmx
      real, dimension(ifull+iextra,kl,5) :: duma
      integer, save :: num_hight = 0, numunstab = 0

      call start_log(upglobal_begin)
      
      intsch=mod(ktau,2)

      if(m>=5)then
        do k=1,kl
         do iq=1,ifull   ! here still unstaggered
!         finish off RHS terms; this coriolis term was once in nonlin
          ux(iq,k)=ux(iq,k)+.5*dt*(1.-epsf)*f(iq)*v(iq,k) ! end of Eq. 129
          vx(iq,k)=vx(iq,k)-.5*dt*(1.-epsf)*f(iq)*u(iq,k) ! end of Eq. 130
         enddo     ! iq loop
        enddo      ! k loop
      endif  ! (m>=5)

      select case(ndept)
       case(0)
        call depts(x3d,y3d,z3d)
       case(1)
        call depts1(x3d,y3d,z3d)
      end select
      
      if(npex==-1)then   ! extrap. Adams-Bashforth style 
!       may not be as good as usual jlm method, because 
!       will not handle residual zs terms as accurately       
        if (.not.allocated(tnsav)) then
          allocate(tnsav(ifull,kl),unsav(ifull,kl),vnsav(ifull,kl))
        end if
        if(ktau==1)then
	    tnsav(:,:) =tn(:,:)
	    unsav(:,:) =un(:,:)
	    vnsav(:,:) =vn(:,:)
        endif
        tx(1:ifull,:)=tx(1:ifull,:)
     &                +dt*(tn(1:ifull,:)-.5*tnsav(1:ifull,:))            
        ux(1:ifull,:)=ux(1:ifull,:)
     &                +dt*(un(1:ifull,:)-.5*unsav(1:ifull,:))           
        vx(1:ifull,:)=vx(1:ifull,:)
     &                +dt*(vn(1:ifull,:)-.5*vnsav(1:ifull,:))
!       convert un, vn to cartesian velocity components (for next step)
        do k=1,kl
         do iq=1,ifull
          uc(iq,k)=ax(iq)*un(iq,k) + bx(iq)*vn(iq,k)
          vc(iq,k)=ay(iq)*un(iq,k) + by(iq)*vn(iq,k)
          wc(iq,k)=az(iq)*un(iq,k) + bz(iq)*vn(iq,k)
         enddo                 ! iq loop
        enddo
        if(mup/=0)then
          duma(1:ifull,:,1)=uc(1:ifull,:)
          duma(1:ifull,:,2)=vc(1:ifull,:)
          duma(1:ifull,:,3)=wc(1:ifull,:)
          duma(1:ifull,:,4)=tn(1:ifull,:)
          call ints(4,duma,intsch,nface,xg,yg,3)
          uc(1:ifull,:)=duma(1:ifull,:,1)
          vc(1:ifull,:)=duma(1:ifull,:,2)
          wc(1:ifull,:)=duma(1:ifull,:,3)
          tn(1:ifull,:)=duma(1:ifull,:,4)
        endif
!       don't bother to rotate un,vn vector to arrival point
!       convert back to un,vn conformal-cubic velocity components
        do k=1,kl
         do iq=1,ifull
          unsav(iq,k)=ax(iq)*uc(iq,k)+ay(iq)*vc(iq,k) +az(iq)*wc(iq,k)
          vnsav(iq,k)=bx(iq)*uc(iq,k)+by(iq)*vc(iq,k) +bz(iq)*wc(iq,k)
          tnsav(iq,k)=tn(iq,k) ! for Ad-Bash use next timestep
         enddo                 ! iq loop
        enddo
      endif   ! (npex==-1)
      
!     calculate factr for choice of nt_adv, as usually used
      if(nt_adv==0) then
        factr(:)=0.
      else if(nt_adv==3) then   ! 1. up to sig=.3
        do k=1,kl
                   factr(k)=stdlapse*(rdry*nritch_t/grav)
                   if(sig(k)<0.3)factr(k)=0.
        end do
       else if(nt_adv==4) then ! (1, .5, 0) for sig=(1, .75, .5)
        do k=1,kl
               factr(k)=max(2.*sig(k)-1., 0.)*stdlapse*(rdry*300./grav)
        end do
       else if(nt_adv==5) then   ! 1 to 0 for sig=1 to 0
        do k=1,kl
               factr(k)=sig(k)*stdlapse*(rdry*nritch_t/grav)
        end do       
       else if(nt_adv==6) then   ! (1, .5625, 0) for sig=(1, .5, .2)
        do k=1,kl
               factr(k)=max(0.,1.25*(sig(k)-.2)*(2.-sig(k)))
     &                   *stdlapse*(rdry*nritch_t/grav)
        end do
       else if(nt_adv==7) then   ! 1 up to .4, then lin decr. to .2, then 0
        do k=1,kl
               factr(k)=max(0.,min(1.,(sig(k)-.2)/(.4-.2)))
     &                   *stdlapse*(rdry*nritch_t/grav)
        end do
       else if(nt_adv==8) then   ! .8 up to .4, then lin decr. to .2, then 0
        do k=1,kl
               factr(k)=.8*max(0.,min(1.,(sig(k)-.2)/(.4-.2)))
     &                   *stdlapse*(rdry*nritch_t/grav)
        end do
       else if(nt_adv==9)then ! (1,1,.84375,.5,.15625,0) for sig=(1,.6,.5,.4,.3,.2)  
        do k=1,kl
                factr(k)=3.*((sig(k)-.2)/.4)**2 -2.*((sig(k)-.2)/.4)**3
     &           *stdlapse*(rdry*nritch_t/grav)
                if(sig(k)>.6)factr(k)=stdlapse*(rdry*nritch_t/grav)
                if(sig(k)<.2)factr(k)=0.
        end do
       else if(nt_adv==10)then ! (1,1,.741,.259,0) for sig=(1,.5,.4,.3,.2) 
        do k=1,kl
                factr(k)=3.*((sig(k)-.2)/.3)**2 -2.*((sig(k)-.2)/.3)**3
     &           *stdlapse*(rdry*nritch_t/grav)
                if(sig(k)>.5)factr(k)=stdlapse*(rdry*nritch_t/grav)
                if(sig(k)<.2)factr(k)=0.
        end do
       endif

#ifdef debug
      if ( mydiag ) then
         if(tx(idjd,kl)>264.)then  !cb
           write(6,*)
           write(6,*) 'upglobal a',ktau,id,jd,tx(idjd,kl)
         endif    ! (tx(idjd,kl)>264.)
      end if

      if(num_hight<100)then
        do iq=1,ifull
         if(tx(iq,kl)>264.)then  !cb
           write(6,*) 'upglobal ktau,myid,iq,large_tx  ',ktau,myid,iq,
     &           tx(iq,kl)
           write (6,"('sdot_iq',9f7.3/7x,9f7.3)") sdot(iq,1:kl)
	    num_hight=num_hight+1
         endif
	 enddo
      endif 
#endif          

      aa(1:ifull)=zs(1:ifull)/(rdry*nritch_t)    ! save zs/(r*t) for nt_adv schemes 
      do k=1,kl   
       dd(1:ifull,k)=aa(1:ifull)
      end do     ! k loop

!-------------------------moved up here May 06---------------------------
      if(nvad==-4.and.nvadh/=3)then 
!       N.B. this moved one is doing vadv on just extra pslx terms      
        sdmx(:) = maxval(abs(sdot),2)
	  nits(:)=1+sdmx(:)/nvadh
	  nvadh_pass(:)=nvadh*nits(:) ! use - for nvadu
        dumt=tx(1:ifull,:)
        dumu=ux(1:ifull,:)
        dumv=vx(1:ifull,:)
        call vadvtvd(dumt,dumu,dumv,nvadh_pass,nits,iaero)
        tx(1:ifull,:)=dumt
        ux(1:ifull,:)=dumu
        vx(1:ifull,:)=dumv 
#ifdef debug
        if( (diag.or.nmaxpr==1) .and. mydiag )then
          write(6,*) 'in upglobal after vadv1'
          write (6,"('qg  ',3p9f8.3/4x,9f8.3)")   qg(idjd,:)
        endif
#endif
      endif  ! (nvad==-4.and.nvadh.ne.3)
!------------------------------------------------------------------

      do k=1,kl   
!      N.B. [D + dsigdot/dsig] saved in adjust5 (or updps) as pslx
       pslx(1:ifull,k)=psl(1:ifull)-pslx(1:ifull,k)*dt*.5*(1.-epst(:))
       pslx(1:ifull,k)=pslx(1:ifull,k)+aa(1:ifull)
       tx(1:ifull,k)=tx(1:ifull,k)+aa(1:ifull)*factr(k)   !cy  
      end do   ! k

#ifdef debug
      if(nmaxpr==1.and.nproc==1)then
        write(6,*) 'pslx_3p before advection'
        write (6,"('pslx_b',3p9f8.4)") pslx(idjd,:)
        write (6,"(i6,8i8)") (ii,ii=id-4,id+4)
        write (6,"(3p9f8.4)") 
     &        ((pslx(max(min(ii+jj*il,ifull),1),nlv),ii=idjd-4,idjd+4)
     &        ,jj=2,-2,-1)
      endif
#endif

      if(mup/=0)then
        call ints_bl(dd,intsch,nface,xg,yg)  ! advection on all levels
        if (nh/=0) then
          duma(1:ifull,:,1)=pslx(1:ifull,:)
          duma(1:ifull,:,2)=h_nh(1:ifull,:)
          call ints(2,duma,intsch,nface,xg,yg,1)
          pslx(1:ifull,:)=duma(1:ifull,:,1)
          h_nh(1:ifull,:)=duma(1:ifull,:,2)
        else
          call ints(1,pslx,intsch,nface,xg,yg,1)
        end if ! nh/=0
        call ints(1,tx,intsch,nface,xg,yg,3)
      endif    ! mup/=0

      do k=1,kl
       pslx(1:ifull,k)=pslx(1:ifull,k)-dd(1:ifull,k)      
       tx(1:ifull,k) = tx(1:ifull,k)  -dd(1:ifull,k)*factr(k)
      end do
!------------------------------------------------------------------
#ifdef debug
	if(nmaxpr==1.and.nproc==1)then
         write(6,*) 'pslx_3p & dd after advection'
         write (6,"('pslx_a',3p9f8.4)") pslx(idjd,:)
         write (6,"('aa#',3p9f8.4)") 
     &             ((aa(ii+jj*il),ii=idjd-1,idjd+1),jj=-1,1)
         write (6,"('dd1#',3p9f8.4)") 
     &             ((dd(ii+jj*il,1),ii=idjd-1,idjd+1),jj=-1,1)
         write (6,"('dd_a',3p9f8.4)") dd(idjd,:)
         write (6,"('nface',18i4)") nface(idjd,:)
         write (6,"('xg',9f8.4)") xg(idjd,:)
         write (6,"('yg',9f8.4)") yg(idjd,:)
         write (6,"(i6,8i8)") (ii,ii=id-4,id+4)
         idjdd=max(5+2*il,min(idjd,ifull-4-2*il))  ! for following prints
         write (6,"(3p9f8.4)") 
     &            ((pslx(ii+jj*il,nlv),ii=idjdd-4,idjdd+4),jj=2,-2,-1)
         uc(1:ifull,1)=-pslx(1:ifull,1)*dsig(1) 
         do k=2,kl
          uc(1:ifull,1)=uc(1:ifull,1)-pslx(1:ifull,k)*dsig(k)
         enddo
         write(6,*) 'integ pslx after advection'
         write (6,"(i6,8i8)") (ii,ii=id-4,id+4)
         write (6,"(3p9f8.4)") 
     &            ((uc(ii+jj*il,1),ii=idjdd-4,idjdd+4),jj=2,-2,-1)
         write(6,*) 'corresp integ ps after advection'
         write (6,"(i6,8i8)") (ii,ii=id-4,id+4)
         write (6,"(-2p9f8.2)") 
     &        ((1.e5*exp(uc(ii+jj*il,1)),ii=idjdd-4,idjdd+4),jj=2,-2,-1)
	endif
!     now comes ux & vx section
      if(diag)then
         if ( mydiag ) then
            write(6,*)
     &        'unstaggered now as uavx and vavx: globpea uses ux & vx'
            write(6,*) 'ux ',(ux(idjd,kk),kk=1,nlv)
         end if
         call printa('uavx',ux,ktau,nlv,ia,ib,ja,jb,0.,1.)
         if (mydiag) write(6,*) 'vx ',(vx(idjd,kk),kk=1,nlv)
         call printa('vavx',vx,ktau,nlv,ia,ib,ja,jb,0.,1.)
         if ( mydiag ) write(6,*)
     &      'unstaggered u and v as uav and vav: globpea uses u & v'
         call printa('uav ',u,ktau,nlv,ia,ib,ja,jb,0.,1.)
         call printa('vav ',v,ktau,nlv,ia,ib,ja,jb,0.,1.)
      endif
#endif

!      convert uavx, vavx to cartesian velocity components
       do k=1,kl
          do iq=1,ifull
             uc(iq,k)=ax(iq)*ux(iq,k) + bx(iq)*vx(iq,k)
             vc(iq,k)=ay(iq)*ux(iq,k) + by(iq)*vx(iq,k)
             wc(iq,k)=az(iq)*ux(iq,k) + bz(iq)*vx(iq,k)
          enddo                 ! iq loop
       enddo
#ifdef debug
       if(diag)then
          if ( mydiag ) then
             write(6,*) 'uc,vc,wc before advection'
             write (6,'(a,18e20.10)') 'uc,vc,wc '
     &                           ,uc(idjd,nlv),vc(idjd,nlv),wc(idjd,nlv)
          end if
          call printa('uc  ',uc,ktau,nlv,ia,ib,ja,jb,0.,1.)
          call printa('vc  ',vc,ktau,nlv,ia,ib,ja,jb,0.,1.)
          call printa('wc  ',wc,ktau,nlv,ia,ib,ja,jb,0.,1.)
          call printa('xg  ',xg,ktau,nlv,ia,ib,ja,jb,0.,1.)
          call printa('yg  ',yg,ktau,nlv,ia,ib,ja,jb,0.,1.)
          if ( mydiag ) write(6,*) 'nface ',nface(idjd,:)
       endif
#endif
       if(mup/=0)then
          duma(1:ifull,:,1)=uc(1:ifull,:)
          duma(1:ifull,:,2)=vc(1:ifull,:)
          duma(1:ifull,:,3)=wc(1:ifull,:)
          call ints(3,duma,intsch,nface,xg,yg,2)
          uc(1:ifull,:)=duma(1:ifull,:,1)
          vc(1:ifull,:)=duma(1:ifull,:,2)
          wc(1:ifull,:)=duma(1:ifull,:,3)
       endif
#ifdef debug
       if(diag)then
          if ( mydiag ) write(6,*) 'uc,vc,wc after advection'
          call printa('uc  ',uc,ktau,nlv,ia,ib,ja,jb,0.,1.)
          call printa('vc  ',vc,ktau,nlv,ia,ib,ja,jb,0.,1.)
          call printa('wc  ',wc,ktau,nlv,ia,ib,ja,jb,0.,1.)
       endif
#endif

       if(nrot==1)then
!        rotate wind vector to arrival point
          do k=1,kl
             do iq=1,ifull
!         the following normalization may be done, but has ~zero effect
!         dena=sqrt(x3d(iq)**2+y3d(iq)**2+z3d(iq)**2)
!         x3d(iq)=x3d(iq)/dena
!         y3d(iq)=y3d(iq)/dena
!         z3d(iq)=z3d(iq)/dena
!         cross product n1xn2 into vec1
                vec1x = y3d(iq,k)*z(iq) - y(iq)*z3d(iq,k)
                vec1y = z3d(iq,k)*x(iq) - z(iq)*x3d(iq,k)
                vec1z = x3d(iq,k)*y(iq) - x(iq)*y3d(iq,k)
                denb = vec1x**2 + vec1y**2 + vec1z**2
!         N.B. rotation formula is singular for small denb,
!         but the rotation is unnecessary in this case
                if(denb>1.e-4)then
                   vecdot = x3d(iq,k)*x(iq) + y3d(iq,k)*y(iq) +
     &                      z3d(iq,k)*z(iq)
                   vec2x = x3d(iq,k)*vecdot - x(iq)
                   vec2y = y3d(iq,k)*vecdot - y(iq)
                   vec2z = z3d(iq,k)*vecdot - z(iq)
                   vec3x = x3d(iq,k) - vecdot*x(iq)
                   vec3y = y3d(iq,k) - vecdot*y(iq)
                   vec3z = z3d(iq,k) - vecdot*z(iq)
                   vdot1 = (vec1x*uc(iq,k) + vec1y*vc(iq,k) +
     &                      vec1z*wc(iq,k))/denb
                   vdot2 = (vec2x*uc(iq,k) + vec2y*vc(iq,k) +
     &                      vec2z*wc(iq,k))/denb
                   uc(iq,k) = vdot1*vec1x + vdot2*vec3x
                   vc(iq,k) = vdot1*vec1y + vdot2*vec3y
                   wc(iq,k) = vdot1*vec1z + vdot2*vec3z
                endif           ! (denb>1.e-4)
             enddo              ! iq loop
          end do ! k
          if(diag)then
            if ( mydiag )then
	         iq=idjd
	         k=nlv
                vec1x = y3d(iq,k)*z(iq) - y(iq)*z3d(iq,k)
                vec1y = z3d(iq,k)*x(iq) - z(iq)*x3d(iq,k)
                vec1z = x3d(iq,k)*y(iq) - x(iq)*y3d(iq,k)
                denb = vec1x**2 + vec1y**2 + vec1z**2
	         write(6,*) 'uc,vc,wc after nrot; denb = ',denb
	     endif
            call printa('uc  ',uc,ktau,nlv,ia,ib,ja,jb,0.,1.)
            call printa('vc  ',vc,ktau,nlv,ia,ib,ja,jb,0.,1.)
            call printa('wc  ',wc,ktau,nlv,ia,ib,ja,jb,0.,1.)
          endif
       endif       ! nrot==1

!      convert back to conformal-cubic velocity components (unstaggered)
!      globpea: this can be sped up later
       do k=1,kl
          do iq=1,ifull
             ux(iq,k) = ax(iq)*uc(iq,k) + ay(iq)*vc(iq,k) +
     &                  az(iq)*wc(iq,k)
             vx(iq,k) = bx(iq)*uc(iq,k) + by(iq)*vc(iq,k) +
     &                  bz(iq)*wc(iq,k)
          enddo ! iq loop
       end do   ! k loop
c      nvsplit=3,4 stuff moved down before or after Coriolis on 15/3/07       
       if(npex==1)then
         ux(1:ifull,:)=ux(1:ifull,:)+.5*dt*un(1:ifull,:) ! dyn contrib
         vx(1:ifull,:)=vx(1:ifull,:)+.5*dt*vn(1:ifull,:) ! dyn contrib
       endif
#ifdef debug
       if(diag.and.k==nlv)then
          if ( mydiag ) write(6,*)
     &         'after advection in upglobal; unstaggered ux and vx:'
          call printa('ux  ',ux,ktau,nlv,ia,ib,ja,jb,0.,1.)
          call printa('vx  ',vx,ktau,nlv,ia,ib,ja,jb,0.,1.)
       endif
#endif

       if(mspec==1.and.mup/=0)then   ! advect qg after preliminary step
          if(ldr/=0)then
             duma(1:ifull,:,1)=qg(1:ifull,:)
             duma(1:ifull,:,2)=qlg(1:ifull,:)
             duma(1:ifull,:,3)=qfg(1:ifull,:)
             duma(1:ifull,:,4)=qrg(1:ifull,:)
             duma(1:ifull,:,5)=cffall(1:ifull,:)
             call ints(5,duma,intsch,nface,xg,yg,4)
             qg(1:ifull,:)    =duma(1:ifull,:,1)
             qlg(1:ifull,:)   =duma(1:ifull,:,2)
             qfg(1:ifull,:)   =duma(1:ifull,:,3)
             qrg(1:ifull,:)   =duma(1:ifull,:,4)
             cffall(1:ifull,:)=duma(1:ifull,:,5)
             cffall=min(max(cffall,0.),1.)
          else
             call ints(1,qg,intsch,nface,xg,yg,3)
          endif                 ! ldr.ne.0
          if(ngas>0.or.nextout>=4)then
#ifdef debug
	     if(nmaxpr==1.and.mydiag)then
              write (6,"('xg#',9f8.2)") diagvals(xg(:,nlv))
              write (6,"('yg#',9f8.2)") diagvals(yg(:,nlv))
              write (6,"('nface#',9i8)") diagvals(nface(:,nlv))
!     &           ((nface(ii+jj*il,nlv),ii=idjd-1,idjd+1),jj=-1,1)
              write (6,"('xlat#',9f8.2)") diagvals(tr(:,nlv,ngas+1))
!     &           ((tr(ii+jj*il,nlv,ngas+1),ii=idjd-1,idjd+1),jj=-1,1)
              write (6,"('xlon#',9f8.2)") diagvals(tr(:,nlv,ngas+2))
!     &           ((tr(ii+jj*il,nlv,ngas+2),ii=idjd-1,idjd+1),jj=-1,1)
              write (6,"('xpre#',9f8.2)") diagvals(tr(:,nlv,ngas+3))
!     &           ((tr(ii+jj*il,nlv,ngas+3),ii=idjd-1,idjd+1),jj=-1,1)
	     endif
#endif
           if ( ntrac > 0 ) then
            call ints(ntrac,tr,intsch,nface,xg,yg,5)
           end if
#ifdef debug
	     if(nmaxpr==1.and.mydiag)then
              write (6,"('ylat#',9f8.2)") diagvals(tr(:,nlv,ngas+1))
              write (6,"('ylon#',9f8.2)") diagvals(tr(:,nlv,ngas+2))
              write (6,"('ypre#',9f8.2)") diagvals(tr(:,nlv,ngas+3))
	     endif
#endif
          endif  ! (ngas>0.or.nextout>=4)
          if(nvmix==6)then
             duma(1:ifull,:,1)=tke(1:ifull,:)
             duma(1:ifull,:,2)=eps(1:ifull,:)
             call ints(2,duma,intsch,nface,xg,yg,3)
             tke(1:ifull,:)=duma(1:ifull,:,1)
             eps(1:ifull,:)=duma(1:ifull,:,2)
          endif                 ! nvmix==6
          if (abs(iaero)==2) then
            call ints(naero,xtg,intsch,nface,xg,yg,5)
          end if
       endif     ! mspec==1

      if(m>=6)then
!       second part of usual m=6 coriolis treatment (before vadv)
        do k=1,kl
!        incorporate coriolis terms (done here as for m=6 instead of in adjust5)
         do iq=1,ifull
          tempry   = ux(iq,k)+.5*dt*(1.+epsf)*f(iq)*vx(iq,k) ! Eq. 133
          vx(iq,k) = vx(iq,k)-.5*dt*(1.+epsf)*f(iq)*ux(iq,k) ! Eq. 134
          ux(iq,k) = tempry
         enddo              ! iq loop
        enddo
      endif  ! (m>=6)
      if(npex==2)then   ! adding a bit later than for npex=1
        ux(1:ifull,:)=ux(1:ifull,:)+.5*dt*un(1:ifull,:) ! dyn contrib
        vx(1:ifull,:)=vx(1:ifull,:)+.5*dt*vn(1:ifull,:) ! dyn contrib
      endif

      if(nvadh==2.or.nvadh==3)then                 ! final dt/2 's worth
        if(nvad==-4)then
          sdot(:,2:kl)=sbar(:,:)
          if(nvadh==3)then  ! just done here after horiz adv
            sdmx(:) = maxval(abs(sdot),2)
            nits(:)=1+sdmx(:)   ! effectively takes nvadh=1  1/2/06
            nvadh_pass(:)=nits(:) ! use - for nvadu
          endif   ! (nvadh==3)
#ifdef debug
	   if(mod(ktau,nmaxpr)==0.and.mydiag)
     &       write(6,*) 'upglobal ktau,sdmx,nits,nvadh_pass ',
     &                ktau,sdmx(idjd),nits(idjd),nvadh_pass(idjd)
          if( (diag.or.nmaxpr==1) .and. mydiag )then
            write(6,*) 'in upglobal before vadv2'
            write (6,"('qg  ',3p9f8.3/4x,9f8.3)")   qg(idjd,:)
          endif
#endif
          dumt=tx(1:ifull,:)
          dumu=ux(1:ifull,:)
          dumv=vx(1:ifull,:)
          call vadvtvd(dumt,dumu,dumv,nvadh_pass,nits,iaero)
          tx(1:ifull,:)=dumt
          ux(1:ifull,:)=dumu
          vx(1:ifull,:)=dumv 
#ifdef debug
          if( (diag.or.nmaxpr==1) .and. mydiag )then
            write(6,*) 'in upglobal after vadv2'
            write (6,"('qg  ',3p9f8.3/4x,9f8.3)")   qg(idjd,:)
          endif
#endif
        endif   ! (nvad==-4)
      endif     ! (nvadh==2.or.nvadh==3)

      if(npex==0.or.npex==5)then ! adding later (after 2nd vadv) than for npex=1
        ux(1:ifull,:)=ux(1:ifull,:)+.5*dt*un(1:ifull,:) ! dyn contrib
        vx(1:ifull,:)=vx(1:ifull,:)+.5*dt*vn(1:ifull,:) ! dyn contrib
      endif
      if(nvsplit==3.or.nvsplit==4)then
        ux(1:ifull,:)=ux(1:ifull,:)+.5*dt*unn(1:ifull,:) ! phys contrib
        vx(1:ifull,:)=vx(1:ifull,:)+.5*dt*vnn(1:ifull,:) ! phys contrib
      endif
      
      if(m==5)then
!       second part of usual m=6 coriolis treatment (after 2nd vadv)
        do k=1,kl
!        incorporate coriolis terms (done here as for m=6 instead of in adjust5)
         do iq=1,ifull
          tempry   = ux(iq,k)+.5*dt*(1.+epsf)*f(iq)*vx(iq,k) ! Eq. 133
          vx(iq,k) = vx(iq,k)-.5*dt*(1.+epsf)*f(iq)*ux(iq,k) ! Eq. 134
          ux(iq,k) = tempry
         enddo              ! iq loop
        enddo
      endif  ! (m==5)

      if(npex>=0)tx(1:ifull,:)=tx(1:ifull,:)+.5*dt*tn(1:ifull,:) 

!     now interpolate ux,vx to the staggered grid
      call staguv(ux,vx,ux,vx)

!     npex=3 add un, vn on staggered grid 
      if(npex==3)then  
        ux(1:ifull,:)=ux(1:ifull,:)+.5*dt*un(1:ifull,:)
        vx(1:ifull,:)=vx(1:ifull,:)+.5*dt*vn(1:ifull,:)
      endif
!     npex=6 similar to npex=3, but adds all un, vn here
      if(npex==6)then  
        ux(1:ifull,:)=ux(1:ifull,:)+dt*un(1:ifull,:)
        vx(1:ifull,:)=vx(1:ifull,:)+dt*vn(1:ifull,:)
      endif

#ifdef debug
      if( diag) then
        if(mydiag)then
          write(6,*) 'near end of upglobal staggered ux and vx:'
          write(6,*) 'un_u ',un(idjd,:)
          write(6,*) 'vn_u ',vn(idjd,:)
          write(6,*) 'tn_u ',tn(idjd,:)
          write (6,"('tx_u1',9f8.2/5x,9f8.2)") tx(idjd,:)
        endif
        call printa('ux  ',ux,ktau,nlv,ia,ib,ja,jb,0.,1.)
        call printa('vx  ',vx,ktau,nlv,ia,ib,ja,jb,0.,1.)
      endif

      if(ntest>4)then
!       diagnostic check for unstable layers
        do iq=1,ifull
         do k=1,kl
          theta(iq,k)=tx(iq,k)*sig(k)**(-rdry/cp)
         enddo
         do k=ntest,kl   ! e.g. 8,kl
          if(theta(iq,k)<theta(iq,k-1))then  ! based on tx
           write(6,*)"unstable layer in upglobal for ktau,iq,k's,del ",
     &                           ktau,iq,k-1,k,theta(iq,k-1)-theta(iq,k)
	     write (6,"('theta',9f7.2/5x,9f7.2)") theta(iq,:)
            write (6,"('sdot',9f7.3/4x,9f7.3)")  sdot(iq,1:kl)
            numunstab=numunstab+1
c           if(numunstab==100)stop 'numunstab=30'
          endif
         enddo
        enddo
      endif

      if( ( diag.or.nmaxpr==1) .and. mydiag ) then
        write(6,*) 'near end of upglobal for ktau= ',ktau
        write (6,"('tx_u2',9f8.2/5x,9f8.2)")  tx(idjd,:)
        write (6,"('qg_u',3p9f8.3/4x,9f8.3)") qg(idjd,:)
        write (6,"('ql_u',3p9f8.3/4x,9f8.3)") qlg(idjd,:)
        write (6,"('qf_u',3p9f8.3/4x,9f8.3)") qfg(idjd,:)
      endif 
#endif    

      call end_log(upglobal_end)
      
      return
      end
