      subroutine upglobal      ! globpea version   use ritchie 103
      use cc_mpi
      use diag_m
      implicit none
!     parameter (nrot=1)       ! nrot=1 to rotate velocity vectors (parmdyn.h)
      integer, parameter :: ntest=0       ! ~8+ for diagnostic stability tests
!     m=6 had three options; only option 3 has been adopted here
      include 'newmpar.h'
      include 'arrays.h'
      include 'const_phys.h'
      include 'indices.h'
      include 'liqwpar.h'  ! ifullw
      include 'map.h'
      include 'nlin.h'
      include 'parm.h'
      include 'parmdyn.h'  
      include 'parmhor.h'  ! mhint, m_bs, nt_adv
      include 'sigs.h'
      include 'tracers.h'
      include 'vecsuv.h'   ! vecsuv info
      include 'vvel.h'     ! sdot
      include 'xarrs.h'
      include 'xyzinfo.h'  ! x,y,z,wts
      real epst
      common/epst/epst(ifull)
      real ubar, vbar
      common/uvbar/ubar(ifull,kl),vbar(ifull,kl)
      ! Need work common for these
      real, dimension(ifull+iextra,kl) :: uc, vc, wc, cc, dd
      real aa(ifull,kl)
      real x3d(ifull,kl),y3d(ifull,kl),z3d(ifull,kl)
      integer nface
      real xg, yg
      common/work3f/nface(ifull,kl),xg(ifull,kl),yg(ifull,kl) ! depts, upglobal
      real tnsav, unsav, vnsav
      common/nonlsav/tnsav(ifull,kl),unsav(ifull,kl),vnsav(ifull,kl)
      real theta(kl), factr(kl)
      integer intsch, iq, k, kk, lev, ntr
      real denb, dtm, dtp, phi1, tav, tempry, tsurf, vdot1,
     &     vdot2, vec1x, vec1y, vec1z, vec2x, vec2y, vec2z, vec3x,
     &     vec3y, vec3z, vecdot
      integer, save :: num_hight = 0, numunstab = 0

      call start_log(upglobal_begin)

      intsch=mod(ktau,2)

      if(m.eq.6)then
        do k=1,kl
         do iq=1,ifull
!         finish off RHS terms; this coriolis term was once in nonlin
          ux(iq,k)=ux(iq,k)+.5*dt*(1.-epsf)*f(iq)*v(iq,k)
          vx(iq,k)=vx(iq,k)-.5*dt*(1.-epsf)*f(iq)*u(iq,k)
         enddo     ! iq loop
        enddo      ! k loop
      endif  ! (m.eq.6)

      if(tx(idjd,kl).gt.264..and.mydiag)then
        print *
        print *,'upglobal a',ktau,id,jd,tx(idjd,kl)
      endif           ! (tx(idjd,kl).gt.264.)

      if(num_hight.lt.100)then
        do iq=1,ifull
         if(tx(iq,kl).gt.264.)then
           print *,'upglobal ktau,iq,large_tx  ',ktau,iq,tx(iq,kl)
           write (6,"('sdot_iq',9f7.3/4x,9f7.3)") (sdot(iq,kk),kk=1,kl)
	    num_hight=num_hight+1
         endif
	 enddo
      endif           

!     call depts3d   ! this one has its own k loop; was nvad<0 option

      if(ndept.eq.0) call depts(x3d,y3d,z3d)
      if(ndept.eq.1) call depts1(x3d,y3d,z3d)

      do k=1,kl

       if(nritch.ne.0.or.nt_adv.gt.0)then
         if(nritch_t.lt.0)then   !  e.g. -3
           lev=-nritch_t
           do iq=1,ifull
            phi1=t(iq,lev)*rdry*(1.-sig(lev))/sig(lev) !phi of sig(lev) above sfc
            tsurf=t(iq,lev)+phi1*stdlapse/grav
            tav=tsurf+zs(iq)*.5*stdlapse/grav
            dd(iq,k)=zs(iq)/(rdry*tav)   ! dlnps to give ln(pmsl)
           enddo
         endif
         if(nritch_t.eq.0)then    ! default nritch_t is 0
           do iq=1,ifull
            dd(iq,k)=zs(iq)/(rdry*t(iq,k))
           enddo
         endif
         if(nritch_t.gt.0)then
           do iq=1,ifull
            dd(iq,k)=zs(iq)/(rdry*nritch_t)    ! e.g. nritch_t=300
           enddo
         endif
       endif  ! (nritch.ne.0.or.nt_adv.gt.0) 

       if(nritch.gt.300)then   ! jlm special bi-linear treatment 303
         do iq=1,ifull
          pslx(iq,k)=pslx(iq,k)+dd(iq,k)
         enddo   ! iq loop
       endif     ! nritch.gt.300
      end do  ! k

       if(nritch.gt.300.or.nt_adv.gt.0)then   ! jlm special bi-linear treatment 303
          aa(1:ifull,:)=dd(1:ifull,:)         ! for nt_adv schemes holding zs/(r*t)
          call ints_bl(dd,intsch,nface,xg,yg)
       endif     ! (nritch.gt.300.or.nt_adv.gt.0)

       if(mup.ne.0) call ints(pslx,intsch,nface,xg,yg,1)

       if(nritch.gt.0)then
          do k=1,kl
             do iq=1,ifull
                pslx(iq,k)=pslx(iq,k)-dd(iq,k)
             enddo              ! iq loop
          end do
	endif   ! (nritch.gt.0)  

!      if(nritch.ge.404.and.ktau.eq.1)then  ! set in globpe now
!	  nt_adv=nritch-400   ! for backward compatibility to nritch=407
!	  print *,'resetting nritch,nt_adv,ktau ',nritch,nt_adv,ktau 
!	endif
	
       if(nt_adv.gt.0)then   ! special T advection treatments follow
          do k=1,kl
             if(nt_adv.eq.3) then   ! 1. up to sig=.4
                   factr(k)=stdlapse*(rdry*300./grav)
                   if(sig(k).lt..3)factr(k)=0.
             endif   ! (nt_adv.eq.3)then   
             if(nt_adv.eq.4) ! (1, .5, 0) for sig=(1, .75, .5)
     &          factr(k)=max(2.*sig(k)-1., 0.)*stdlapse*(rdry*300./grav)
             if(nt_adv.eq.5)    ! 1 to 0 for sig=1 to 0
     &             factr(k)=sig(k)*stdlapse*(rdry*300./grav)
             if(nt_adv.eq.6)    ! (1, .5625, 0) for sig=(1, .5, .2)
     &             factr(k)=max(0.,1.25*(sig(k)-.2)*(2.-sig(k)))
     &                   *stdlapse*(rdry*300./grav)
             if(nt_adv.eq.7)    ! 1 up to .4, then lin decr. to .2, then 0
     &             factr(k)=max(0.,min(1.,(sig(k)-.2)/(.4-.2)))
     &                   *stdlapse*(rdry*300./grav)
             if(nt_adv.eq.8)then ! (1,1,.84375,.5,.15625,0) for sig=(1,.6,.5,.4,.2)  
                factr(k)=3.*((sig(k)-.2)/.4)**2 -2.*((sig(k)-.2)/.4)**3
     &           *stdlapse*(rdry*300./grav)
                if(sig(k).gt..6)factr(k)=stdlapse*(rdry*300./grav)
                if(sig(k).lt..2)factr(k)=0.
             endif
             if(nt_adv.eq.9)then ! (1,1,.741,.259,0) for sig=(1,.5,.4,.3,.2) 
                factr(k)=3.*((sig(k)-.2)/.3)**2 -2.*((sig(k)-.2)/.3)**3
     &           *stdlapse*(rdry*300./grav)
                if(sig(k).gt..5)factr(k)=stdlapse*(rdry*300./grav)
                if(sig(k).lt..2)factr(k)=0.
             endif
             do iq=1,ifull
                cc(iq,k)=tx(iq,k)+aa(iq,k)*factr(k)
             enddo              ! iq loop
          end do ! k
          call ints(cc,intsch,nface,xg,yg,3)
          do k=1,kl
             tx(1:ifull,k) = cc(1:ifull,k) - dd(1:ifull,k)*factr(k)
          end do
       else
         call ints(tx,intsch,nface,xg,yg,3)
       endif  ! (nt_adv.ge.4)

!      now comes ux & vx section
       if(diag)then
c       print *,'staggered ux and vx:'
          if ( mydiag ) then
             print *,
     &         'unstaggered now as uavx and vavx: globpea uses ux & vx'
             print *,'ux ',(ux(idjd,kk),kk=1,nlv)
          end if
          call printa('uavx',ux,ktau,nlv,ia,ib,ja,jb,0.,1.)
          if (mydiag) print *,'vx ',(vx(idjd,kk),kk=1,nlv)
          call printa('vavx',vx,ktau,nlv,ia,ib,ja,jb,0.,1.)
          if ( mydiag ) print *,
     &       'unstaggered u and v as uav and vav: globpea uses u & v'
          call printa('uav ',u,ktau,nlv,ia,ib,ja,jb,0.,1.)
          call printa('vav ',v,ktau,nlv,ia,ib,ja,jb,0.,1.)
       endif

!      convert uavx, vavx to cartesian velocity components
       do k=1,kl
          do iq=1,ifull
             uc(iq,k)=ax(iq)*ux(iq,k) + bx(iq)*vx(iq,k)
             vc(iq,k)=ay(iq)*ux(iq,k) + by(iq)*vx(iq,k)
             wc(iq,k)=az(iq)*ux(iq,k) + bz(iq)*vx(iq,k)
          enddo                 ! iq loop
       end do
       if(diag)then
          if ( mydiag ) then
             print *,'uc,vc,wc before advection'
             write (6,'(a,18e20.10)') 'uc,vc,wc '
     &                           ,uc(idjd,nlv),vc(idjd,nlv),wc(idjd,nlv)
          end if
          call printa('uc  ',uc,ktau,nlv,ia,ib,ja,jb,0.,1.)
          call printa('vc  ',vc,ktau,nlv,ia,ib,ja,jb,0.,1.)
          call printa('wc  ',wc,ktau,nlv,ia,ib,ja,jb,0.,1.)
          call printa('xg  ',xg,ktau,nlv,ia,ib,ja,jb,0.,1.)
          call printa('yg  ',yg,ktau,nlv,ia,ib,ja,jb,0.,1.)
          if ( mydiag ) print *,'nface ',(nface(idjd,kk),kk=1,kl)
       endif
       if(mup.ne.0)then
          call ints(uc,intsch,nface,xg,yg,2)
          call ints(vc,intsch,nface,xg,yg,2)
          call ints(wc,intsch,nface,xg,yg,2)
       endif
       if(diag)then
          if ( mydiag ) print *,'uc,vc,wc after advection'
          call printa('uc  ',uc,ktau,nlv,ia,ib,ja,jb,0.,1.)
          call printa('vc  ',vc,ktau,nlv,ia,ib,ja,jb,0.,1.)
          call printa('wc  ',wc,ktau,nlv,ia,ib,ja,jb,0.,1.)
       endif

c      print *,'end  upg. k x3d,y3d,z3d ',
c    .                    k,x3d(idjd),y3d(idjd),z3d(idjd)
       if(nrot.eq.1)then
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
                if(denb.gt.1.e-4)then
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
                endif           ! (denb.gt.1.e-4)
             enddo              ! iq loop
          end do ! k
       endif       ! nrot.eq.1

!      convert back to conformal-cubic velocity components (unstaggered)
!      globpea: this can be sped up later
       do k=1,kl
          do iq=1,ifull
             ux(iq,k) = ax(iq)*uc(iq,k) + ay(iq)*vc(iq,k) +
     &                  az(iq)*wc(iq,k)
             vx(iq,k) = bx(iq)*uc(iq,k) + by(iq)*vc(iq,k) +
     &                  bz(iq)*wc(iq,k)
          enddo                 ! iq loop
       end do
       if(diag.and.k.eq.nlv)then
          if ( mydiag ) print *,
     &         'after advection in upglobal; unstaggered ux and vx:'
          call printa('ux  ',ux,ktau,nlv,ia,ib,ja,jb,0.,1.)
          call printa('vx  ',vx,ktau,nlv,ia,ib,ja,jb,0.,1.)
       endif

       if(mspec.eq.1.and.mup.ne.0)then   ! advect qg after preliminary step
          call ints(qg,intsch,nface,xg,yg,4)
          if(ifullw.gt.1)then
             call ints(qlg,intsch,nface,xg,yg,4)
             call ints(qfg,intsch,nface,xg,yg,4)
          endif                 ! ifullw.gt.1
          if(ilt.gt.1)then
             do ntr=1,ntrac
                call ints(tr(1,1,ntr),intsch,nface,xg,yg,5)
             enddo
          endif ! ilt.gt.1
       endif     ! mspec.eq.1

       if(nonl.eq.3)then   ! Adams-Bashforth style
          call ints(tn,intsch,nface,xg,yg,3)
!         convert un, vn to cartesian velocity components
          do k=1,kl
             do iq=1,ifull
                uc(iq,k) = ax(iq)*un(iq,k) + bx(iq)*vn(iq,k)
                vc(iq,k) = ay(iq)*un(iq,k) + by(iq)*vn(iq,k)
                wc(iq,k) = az(iq)*un(iq,k) + bz(iq)*vn(iq,k)
             enddo     ! iq loop
          end do
          call ints(uc,intsch,nface,xg,yg,2)
          call ints(vc,intsch,nface,xg,yg,2)
          call ints(wc,intsch,nface,xg,yg,2)
!         don't bother to rotate un,vn vector to arrival point
!         convert back to un,vn conformal-cubic velocity components
          do k=1,kl
             do iq=1,ifull
                unsav(iq,k) = ax(iq)*uc(iq,k) + ay(iq)*vc(iq,k) +
     &                        az(iq)*wc(iq,k)
                vnsav(iq,k) = bx(iq)*uc(iq,k) + by(iq)*vc(iq,k) +
     &                        bz(iq)*wc(iq,k)
                tnsav(iq,k) = tn(iq,k) ! for use next step in nonlin
                tn(iq,k)=0.
                un(iq,k)=0.
                vn(iq,k)=0.
             enddo     ! iq loop
          end do
       endif       ! (nonl.eq.3)

!      there are various possibilities for changing wind direction
       dtm=dt*(1.-epsf)
       dtp=dt*(1.+epsf)

       if(m.eq.6)then    ! i.e.second part of m=6 coriolis treatment
          do k=1,kl
             do iq=1,ifull
                ux(iq,k)=ux(iq,k)+.5*dt*un(iq,k)
                vx(iq,k)=vx(iq,k)+.5*dt*vn(iq,k)
             enddo     ! iq loop
!            incorporate coriolis terms (done here as for m=6 instead of in adjust5)
             do iq=1,ifull
                tempry   = ux(iq,k)+.5*dt*(1.+epsf)*f(iq)*vx(iq,k) ! option 3
                vx(iq,k) = vx(iq,k)-.5*dt*(1.+epsf)*f(iq)*ux(iq,k) ! option 3
                ux(iq,k) = tempry
             enddo              ! iq loop
          end do
       else
          print*, "Error, not implemented m=", m
          stop
       end if

!      now interpolate ux,vx to the staggered grid; special globpea call
      call staguv(ux,vx,ux,vx)
      if(diag)then
         if ( mydiag ) print *,
     &        'near end of upglobal staggered ux and vx:'
         call printa('ux  ',ux,ktau,nlv,ia,ib,ja,jb,0.,1.)
         call printa('vx  ',vx,ktau,nlv,ia,ib,ja,jb,0.,1.)
      endif

      if(diag.and.mydiag)then
        print *,'near end of upglobal'
        print *,'tn ',(tn(idjd,k),k=1,kl)
        print *,'tx ',(tx(idjd,k),k=1,kl)
        print *,'un ',(un(idjd,k),k=1,kl)
        print *,'vn ',(vn(idjd,k),k=1,kl)
      endif

      tx(1:ifull,:) = tx(1:ifull,:)+.5*dt*tn(1:ifull,:) ! moved from adjust5 30/11/00

      if(ntest.gt.4)then
!       diag check for unstable layers
        do iq=1,ifull
         do k=1,kl
          theta(k)=tx(iq,k)*sig(k)**(-rdry/cp)
         enddo
         do k=ntest,kl   ! e.g. 8,kl
          if(theta(k).lt.theta(k-1))then  ! based on tx
            print *,"unstable layer in upglobal for ktau,iq,k's,del ",
     &                               ktau,iq,k-1,k,theta(k-1)-theta(k)
	     write (6,"('theta',9f7.2/5x,9f7.2)") theta
            write (6,"('sdot',9f7.3/4x,9f7.3)") (sdot(iq,kk),kk=1,kl)
c           print *,'sdot', (sdot(iq,kk),kk=1,kl)
            numunstab=numunstab+1
c           if(numunstab.eq.100)stop 'numunstab=30'
          endif
         enddo
        enddo
      endif

      call end_log(upglobal_end)
      return
      end
