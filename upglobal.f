      subroutine upglobal      ! globpea version   use ritchie 103
!     parameter (nrot=1)       ! nrot=1 to rotate velocity vectors (parmdyn.h)
      parameter (ntest=0)      ! ~8+ for diagnostic stability tests
!     m=6 had three options; only option 3 has been adopted here
      include 'newmpar.h'
      include 'arrays.h'
      include 'constant.h'
      include 'indices.h'
      include 'liqwpar.h'  ! ifullw
      include 'map.h'
      include 'nlin.h'
      include 'parm.h'
      include 'parmdyn.h'  
      include 'parmhor.h'  ! mhint, m_bs, nt_adv
      include 'particle.h'
      include 'sigs.h'
      include 'tracers.h'
      include 'vecsuv.h'   ! vecsuv info
      include 'vvel.h'     ! sdot
      include 'xarrs.h'
      include 'xyzinfo.h'  ! x,y,z,wts
      common/epst/epst(ifull)
      common/uvbar/ubar(ifull,kl),vbar(ifull,kl)
      common/work2/uavx(ifull),vavx(ifull)       
     .      ,uc(ifull),vc(ifull),wc(ifull) 
     .      ,aa(ifull),bb(ifull),cc(ifull),dd(ifull),cc2(ifull)
     .      ,delcor(ifull),fxx(ifull),dum2(6*ifull)
      real x3d(ifull,kl),y3d(ifull,kl),z3d(ifull,kl)
      common/work3f/nface(ifull,kl),xg(ifull,kl),yg(ifull,kl) ! depts, upglobal
      common/nonlsav/tnsav(ifull,kl),unsav(ifull,kl),vnsav(ifull,kl)
      real theta(kl)
      save numunstab,num_hight
      data numunstab/0/,num_hight/0/
c     ind(i,j,n)=i+(j-1)*il+n*il*il  ! *** for n=0,5
      intsch=mod(ktau,2)
      psavg1=0.
      psavg2=0.
      psmax1=-100.
      psmax2=-100.
      psmin1=0.
      psmin2=0.

      if(m.eq.6)then
        do k=1,kl
         do iq=1,ifull
!         finish off RHS terms; this coriolis term was once in nonlin
          ux(iq,k)=ux(iq,k)+.5*dt*(1.-epsf)*f(iq)*v(iq,k)
          vx(iq,k)=vx(iq,k)-.5*dt*(1.-epsf)*f(iq)*u(iq,k)
         enddo     ! iq loop
        enddo      ! k loop
      endif  ! (m.eq.6)

      if(tx(idjd,kl).gt.264.)then
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

      do k=1,kl  !   start of main k loop

       if(m.eq.7)then    ! i.e. m=7 coriolis treatment
         do iq=1,ifull
!         first find fxx at dt/2 forward trajectory location
!         N.B. ubar & vbar scaled in depts to give grid displ. for dt/2
          delcor(iq)=.5*(ubar(iq,k)*(f(ie(iq))-f(iw(iq)))
     .                    +vbar(iq,k)*(f(in(iq))-f(is(iq))) )
          fxx(iq)=f(iq)+delcor(iq)
          ux(iq,k)=ux(iq,k)+.5*dt*(1.-epsf)*fxx(iq)*v(iq,k)
          vx(iq,k)=vx(iq,k)-.5*dt*(1.-epsf)*fxx(iq)*u(iq,k)
!         then save fxx at dt/2 backward trajectory location
          fxx(iq)=f(iq)-delcor(iq)
         enddo     ! iq loop
       endif  ! (m.eq.7)

       do iq=1,ifull
        cc2(iq)=0.
       enddo   ! iq loop
       if(nritch.ne.0.or.nt_adv.gt.0)then
         if(nritch_t.lt.0)then   !  e.g. -3
           lev=-nritch_t
           do iq=1,ifull
            phi1=t(iq,lev)*r*(1.-sig(lev))/sig(lev) !phi of sig(lev) above sfc
            tsurf=t(iq,lev)+phi1*.0065/g
            tav=tsurf+zs(iq)*.5*.0065/g
            dd(iq)=zs(iq)/(r*tav)   ! dlnps to give ln(pmsl)
           enddo
         endif
         if(nritch_t.eq.0)then    ! default nritch_t is 0
           do iq=1,ifull
            dd(iq)=zs(iq)/(r*t(iq,k))
           enddo
         endif
         if(nritch_t.gt.0)then
           do iq=1,ifull
            dd(iq)=zs(iq)/(r*nritch_t)    ! e.g. nritch_t=300
           enddo
         endif
       endif  ! (nritch.ne.0.or.nt_adv.gt.0) 

       if(nritch.gt.300)then   ! jlm special bi-linear treatment 303
         do iq=1,ifull
          pslx(iq,k)=pslx(iq,k)+dd(iq)
         enddo   ! iq loop
       endif     ! nritch.gt.300

       if(nritch.gt.300.or.nt_adv.gt.0)then   ! jlm special bi-linear treatment 303
         do iq=1,ifull
          aa(iq)=dd(iq)    ! for nt_adv schemes holding zs/(r*t)
         enddo   ! iq loop
         call ints_bl(dd,intsch,nface(1,k),xg(1,k),yg(1,k))
       endif     ! (nritch.gt.300.or.nt_adv.gt.0)

       if(mup.ne.0)
     .   call ints(pslx(1,k),intsch,nface(1,k),xg(1,k),yg(1,k),1)

       if(nritch.gt.0)then
         do iq=1,ifull
          pslx(iq,k)=pslx(iq,k)-dd(iq)+cc2(iq)*(1.+epst(iq))
         enddo   ! iq loop
	endif   ! (nritch.gt.0)  

!      if(nritch.ge.404.and.ktau.eq.1)then  ! set in globpe now
!	  nt_adv=nritch-400   ! for backward compatibility to nritch=407
!	  print *,'resetting nritch,nt_adv,ktau ',nritch,nt_adv,ktau 
!	endif
	
       if(nt_adv.gt.0)then   ! special T advection treatments follow
         if(nt_adv.eq.3)   ! 1. up to sig=.4
     .             factr=6.5e-3*(r*300./g)
                   if(sig(k).lt..3)factr=0.
         if(nt_adv.eq.4)   ! (1, .5, 0) for sig=(1, .75, .5)
     .             factr=max(2.*sig(k)-1. , 0.)*6.5e-3*(r*300./g)
         if(nt_adv.eq.5)   ! 1 to 0 for sig=1 to 0
     .             factr=sig(k)*6.5e-3*(r*300./g)
         if(nt_adv.eq.6)   ! (1, .5625, 0) for sig=(1, .5, .2)
     .             factr=max(0.,1.25*(sig(k)-.2)*(2.-sig(k)))
     .                   *6.5e-3*(r*300./g)
         if(nt_adv.eq.7) ! 1 up to .4, then lin decr. to .2, then 0
     .             factr=max(0.,min(1.,(sig(k)-.2)/(.4-.2)))
     .                   *6.5e-3*(r*300./g)
         if(nt_adv.eq.8)then ! (1,1,.84375,.5,.15625,0) for sig=(1,.6,.5,.4,.2)  
           factr=3.*((sig(k)-.2)/.4)**2 -2.*((sig(k)-.2)/.4)**3
     .           *6.5e-3*(r*300./g)
           if(sig(k).gt..6)factr=6.5e-3*(r*300./g)
           if(sig(k).lt..2)factr=0.
         endif
         if(nt_adv.eq.9)then ! (1,1,.741,.259,0) for sig=(1,.5,.4,.3,.2) 
           factr=3.*((sig(k)-.2)/.3)**2 -2.*((sig(k)-.2)/.3)**3
     .           *6.5e-3*(r*300./g)
           if(sig(k).gt..5)factr=6.5e-3*(r*300./g)
           if(sig(k).lt..2)factr=0.
         endif
         do iq=1,ifull
          cc(iq)=tx(iq,k)+aa(iq)*factr
         enddo   ! iq loop
         call ints(cc,intsch,nface(1,k),xg(1,k),yg(1,k),3)
         do iq=1,ifull
          tx(iq,k)=cc(iq)-dd(iq)*factr
         enddo   ! iq loop
       else
         call ints(tx(1,k),intsch,nface(1,k),xg(1,k),yg(1,k),3)
       endif  ! (nt_adv.ge.4)

!      now comes ux & vx section
       if(diag.and.k.eq.nlv)then
c       print *,'staggered ux and vx:'
        print *,'unstaggered now as uavx and vavx: globpea uses ux & vx'
        print *,'ux ',(ux(idjd,kk),kk=1,nlv)
        call printa('uavx',ux(1,nlv),ktau,nlv,ia,ib,ja,jb,0.,1.)
        print *,'vx ',(vx(idjd,kk),kk=1,nlv)
        call printa('vavx',vx(1,nlv),ktau,nlv,ia,ib,ja,jb,0.,1.)
        print *,'unstaggered u and v as uav and vav: globpea uses u & v'
        call printa('uav ',u(1,nlv),ktau,nlv,ia,ib,ja,jb,0.,1.)
        call printa('vav ',v(1,nlv),ktau,nlv,ia,ib,ja,jb,0.,1.)
       endif

!      convert uavx, vavx to cartesian velocity components
       do iq=1,ifull
        uc(iq)=ax(iq)*ux(iq,k) + bx(iq)*vx(iq,k)
        vc(iq)=ay(iq)*ux(iq,k) + by(iq)*vx(iq,k)
        wc(iq)=az(iq)*ux(iq,k) + bz(iq)*vx(iq,k)
       enddo     ! iq loop
       if(diag.and.k.eq.nlv)then
        print *,'uc,vc,wc before advection'
        write (6,'(a,18e20.10)') 'uc,vc,wc '
     .                           ,uc(idjd),vc(idjd),wc(idjd)
        call printa('uc  ',uc,ktau,nlv,ia,ib,ja,jb,0.,1.)
        call printa('vc  ',vc,ktau,nlv,ia,ib,ja,jb,0.,1.)
        call printa('wc  ',wc,ktau,nlv,ia,ib,ja,jb,0.,1.)
        call printa('xg  ',xg(1,nlv),ktau,nlv,ia,ib,ja,jb,0.,1.)
        call printa('yg  ',yg(1,nlv),ktau,nlv,ia,ib,ja,jb,0.,1.)
	 print *,'nface ',(nface(idjd,kk),kk=1,kl)
       endif
       if(mup.ne.0)then
         call ints(uc,intsch,nface(1,k),xg(1,k),yg(1,k),2)
         call ints(vc,intsch,nface(1,k),xg(1,k),yg(1,k),2)
         call ints(wc,intsch,nface(1,k),xg(1,k),yg(1,k),2)
       endif
       if(diag.and.k.eq.nlv)then
        print *,'uc,vc,wc after advection'
        call printa('uc  ',uc,ktau,nlv,ia,ib,ja,jb,0.,1.)
        call printa('vc  ',vc,ktau,nlv,ia,ib,ja,jb,0.,1.)
        call printa('wc  ',wc,ktau,nlv,ia,ib,ja,jb,0.,1.)
       endif

c      print *,'end  upg. k x3d,y3d,z3d ',
c    .                    k,x3d(idjd),y3d(idjd),z3d(idjd)
       if(nrot.eq.1)then
!        rotate wind vector to arrival point
         do iq=1,ifull
!         the following normalization may be done, but has ~zero effect
!         dena=sqrt(x3d(iq)**2+y3d(iq)**2+z3d(iq)**2)
!         x3d(iq)=x3d(iq)/dena
!         y3d(iq)=y3d(iq)/dena
!         z3d(iq)=z3d(iq)/dena
!         cross product n1xn2 into vec1
          vec1x=y3d(iq,k)*z(iq)-y(iq)*z3d(iq,k)
          vec1y=z3d(iq,k)*x(iq)-z(iq)*x3d(iq,k)
          vec1z=x3d(iq,k)*y(iq)-x(iq)*y3d(iq,k)
          denb=vec1x**2+vec1y**2+vec1z**2
!         N.B. rotation formula is singular for small denb,
!         but the rotation is unnecessary in this case
          if(denb.gt.1.e-4)then
            vecdot=x3d(iq,k)*x(iq)+y3d(iq,k)*y(iq)+z3d(iq,k)*z(iq)
            vec2x=x3d(iq,k)*vecdot-x(iq)
            vec2y=y3d(iq,k)*vecdot-y(iq)
            vec2z=z3d(iq,k)*vecdot-z(iq)
            vec3x=x3d(iq,k)-vecdot*x(iq)
            vec3y=y3d(iq,k)-vecdot*y(iq)
            vec3z=z3d(iq,k)-vecdot*z(iq)
            vdot1=(vec1x*uc(iq) +vec1y*vc(iq) +vec1z*wc(iq))/denb
            vdot2=(vec2x*uc(iq) +vec2y*vc(iq) +vec2z*wc(iq))/denb
            uc(iq)=vdot1*vec1x + vdot2*vec3x
            vc(iq)=vdot1*vec1y + vdot2*vec3y
            wc(iq)=vdot1*vec1z + vdot2*vec3z
          endif    ! (denb.gt.1.e-4)
         enddo     ! iq loop
       endif       ! nrot.eq.1

!      convert back to conformal-cubic velocity components (unstaggered)
!      globpea: this can be sped up later
       do iq=1,ifull
        ux(iq,k)=ax(iq)*uc(iq) + ay(iq)*vc(iq) + az(iq)*wc(iq)
        vx(iq,k)=bx(iq)*uc(iq) + by(iq)*vc(iq) + bz(iq)*wc(iq)
       enddo     ! iq loop
       if(diag.and.k.eq.nlv)then
        print *,'after advection in upglobal; unstaggered ux and vx:'
        call printa('ux  ',ux(1,nlv),ktau,nlv,ia,ib,ja,jb,0.,1.)
        call printa('vx  ',vx(1,nlv),ktau,nlv,ia,ib,ja,jb,0.,1.)
       endif

       if(mspec.eq.1.and.mup.ne.0)then   ! advect qg after preliminary step
        call ints(qg(1,k),intsch,nface(1,k),xg(1,k),yg(1,k),4)
        if(ifullw.gt.1)then
          call ints(qlg(1,k),intsch,nface(1,k),xg(1,k),yg(1,k),4)
          call ints(qfg(1,k),intsch,nface(1,k),xg(1,k),yg(1,k),4)
        endif     ! ifullw.gt.1
        if(ilt.gt.1)then
          do ntr=1,ntrac
           call ints(tr(1,k,ntr),intsch,nface(1,k),xg(1,k),yg(1,k),5)
          enddo
        endif    ! ilt.gt.1
       endif     ! mspec.eq.1

       if(nonl.eq.3)then   ! Adams-Bashforth style
         call ints(tn(1,k),intsch,nface(1,k),xg(1,k),yg(1,k),3)
!        convert un, vn to cartesian velocity components
         do iq=1,ifull
          uc(iq)=ax(iq)*un(iq,k) + bx(iq)*vn(iq,k)
          vc(iq)=ay(iq)*un(iq,k) + by(iq)*vn(iq,k)
          wc(iq)=az(iq)*un(iq,k) + bz(iq)*vn(iq,k)
         enddo     ! iq loop
         call ints(uc,intsch,nface(1,k),xg(1,k),yg(1,k),2)
         call ints(vc,intsch,nface(1,k),xg(1,k),yg(1,k),2)
         call ints(wc,intsch,nface(1,k),xg(1,k),yg(1,k),2)
!        don't bother to rotate un,vn vector to arrival point
!        convert back to un,vn conformal-cubic velocity components
         do iq=1,ifull
          unsav(iq,k)=ax(iq)*uc(iq)+ay(iq)*vc(iq) +az(iq)*wc(iq)
          vnsav(iq,k)=bx(iq)*uc(iq)+by(iq)*vc(iq) +bz(iq)*wc(iq)
          tnsav(iq,k)=tn(iq,k)   ! for use next step in nonlin
          tn(iq,k)=0.
          un(iq,k)=0.
          vn(iq,k)=0.
         enddo     ! iq loop
       endif       ! (nonl.eq.3)

!      there are various possibilities for changing wind direction
       dtm=dt*(1.-epsf)
       dtp=dt*(1.+epsf)
       if(m.eq.1)then
!        incorporate coriolis terms (McDonald style)
         do iq=1,ifull
!         firstly the implicit Coriolis step giving (tau+.5) values
          tempux    =(ux(iq,k)+.5*dtm*f(iq)*vx(iq,k))/
     .                      (1.+(.5*dtm*f(iq))**2)
          tempvx    =(vx(iq,k)-.5*dtm*f(iq)*ux(iq,k))/
     .                      (1.+(.5*dtm*f(iq))**2)
!         then add the explicit Coriolis part (can add un, vn here)
          ux(iq,k)=tempux+.5*(dtp*f(iq)*tempvx+dt*un(iq,k))
          vx(iq,k)=tempvx-.5*(dtp*f(iq)*tempux-dt*vn(iq,k))
         enddo     ! iq loop
       endif  ! (m.eq.1)
       if(m.eq.2)then
         do iq=1,ifull
!         N.B. ubar & vbar scaled in depts to give grid displ. for dt/2
          delcor(iq)=.5*(ubar(iq,k)*(f(ie(iq))-f(iw(iq)))
     .                    +vbar(iq,k)*(f(in(iq))-f(is(iq))) )
!         fxx  is f at dt/2 backward trajectory location
          fxx(iq)=f(iq)-delcor(iq)
         enddo     ! iq loop
!        incorporate coriolis terms (McDonald style)
         do iq=1,ifull
!         firstly the implicit Coriolis step giving (tau+.5) values
          tempux    =(ux(iq,k)+.5*dtm*fxx(iq)*vx(iq,k))/
     .                      (1.+(.5*dtm*fxx(iq))**2)
          tempvx    =(vx(iq,k)-.5*dtm*fxx(iq)*ux(iq,k))/
     .                      (1.+(.5*dtm*fxx(iq))**2)
!         then add the explicit Coriolis part (can add un, vn here)
          ux(iq,k)=tempux+.5*(dtp*fxx(iq)*tempvx+dt*un(iq,k))
          vx(iq,k)=tempvx-.5*(dtp*fxx(iq)*tempux-dt*vn(iq,k))
         enddo   ! iq loop
       endif     ! (m.eq.2)
       if(m.eq.3)then
         do iq=1,ifull
          fxx(iq)=f(iq)
         enddo     ! iq loop
         call ints(fxx,intsch,nface(1,k),xg(1,k),yg(1,k),3)
         do iq=1,ifull
          fxx(iq)=.5*(f(iq)+fxx(iq))
         enddo     ! iq loop
!        fxx  is f at dt/2 backward trajectory location
!        incorporate coriolis terms (McDonald style)
         do iq=1,ifull
!         firstly the implicit Coriolis step giving (tau+.5) values
          tempux    =(ux(iq,k)+.5*dtm*fxx(iq)*vx(iq,k))/
     .                      (1.+(.5*dtm*fxx(iq))**2)
          tempvx    =(vx(iq,k)-.5*dtm*fxx(iq)*ux(iq,k))/
     .                      (1.+(.5*dtm*fxx(iq))**2)
!         then add the explicit Coriolis part (can add un, vn here)
          ux(iq,k)=tempux+.5*(dtp*fxx(iq)*tempvx+dt*un(iq,k))
          vx(iq,k)=tempvx-.5*(dtp*fxx(iq)*tempux-dt*vn(iq,k))
         enddo   ! iq loop
       endif     ! (m.eq.3)
       if(m.eq.4)then
         do iq=1,ifull
          fxx(iq)=f(iq)
         enddo     ! iq loop
         call ints(fxx,intsch,nface(1,k),xg(1,k),yg(1,k),3)
!        here fxx  is f at dt backward trajectory location
         do iq=1,ifull
          aa(iq)=.75*fxx(iq)+.25*f(iq)
          bb(iq)=.25*fxx(iq)+.75*f(iq)
         enddo     ! iq loop
!        incorporate coriolis terms (McDonald style)
         do iq=1,ifull
!         firstly the implicit Coriolis step giving (tau+.5) values
          tempux    =(ux(iq,k)+.5*dtm*aa(iq)*vx(iq,k))/
     .                      (1.+(.5*dtm*aa(iq))**2)
          tempvx    =(vx(iq,k)-.5*dtm*aa(iq)*ux(iq,k))/
     .                      (1.+(.5*dtm*aa(iq))**2)
!         then add the explicit Coriolis part (can add un, vn here)
          ux(iq,k)=tempux+.5*(dtp*bb(iq)*tempvx+dt*un(iq,k))
          vx(iq,k)=tempvx-.5*(dtp*bb(iq)*tempux-dt*vn(iq,k))
         enddo   ! iq loop
       endif     ! (m.eq.4)

       if(m.eq.6)then    ! i.e.second part of m=6 coriolis treatment
         do iq=1,ifull
          ux(iq,k)=ux(iq,k)+.5*dt*un(iq,k)
          vx(iq,k)=vx(iq,k)+.5*dt*vn(iq,k)
         enddo     ! iq loop
!        incorporate coriolis terms (done here as for m=6 instead of in adjust5)
         do iq=1,ifull
          tempry    =ux(iq,k)+.5*dt*(1.+epsf)*f(iq)*vx(iq,k)    ! option 3
          vx(iq,k)=vx(iq,k)-.5*dt*(1.+epsf)*f(iq)*ux(iq,k)    ! option 3
          ux(iq,k)=tempry
         enddo     ! iq loop
       endif  ! (m.eq.6)
       if(m.eq.7)then    ! i.e. second part of m=7 coriolis treatment
         do iq=1,ifull
          ux(iq,k)=ux(iq,k)+.5*dt*un(iq,k)
          vx(iq,k)=vx(iq,k)+.5*dt*vn(iq,k)
         enddo     ! iq loop
!        incorporate coriolis terms (done here as for m=6 instead of below)
         do iq=1,ifull
          tempry    =ux(iq,k)+.5*dt*(1.+epsf)*fxx(iq)*vx(iq,k)
          vx(iq,k)=vx(iq,k)-.5*dt*(1.+epsf)*fxx(iq)*ux(iq,k)
          ux(iq,k)=tempry
         enddo     ! iq loop
c        print *,'m,id,jd,k,u,v,f,fxx: ',m,id,jd,k,u(idjd,k),v(idjd,k),
c    .                              f(idjd),fxx(idjd)
       endif  ! (m.eq.7)

!      now interpolate ux,vx to the staggered grid; special globpea call
      enddo      ! main k loop

      call staguv(ux,vx,ux,vx)
      if(diag)then
         print *,'near end of upglobal staggered ux and vx:'
         call printa('ux  ',ux(1,nlv),ktau,nlv,ia,ib,ja,jb,0.,1.)
         call printa('vx  ',vx(1,nlv),ktau,nlv,ia,ib,ja,jb,0.,1.)
      endif

      if(diag)then
        print *,'near end of upglobal'
        print *,'tn ',(tn(idjd,k),k=1,kl)
        print *,'tx ',(tx(idjd,k),k=1,kl)
        print *,'un ',(un(idjd,k),k=1,kl)
        print *,'vn ',(vn(idjd,k),k=1,kl)
      endif

      tx(:,:)=tx(:,:)+.5*dt*tn(:,:) ! moved from adjust5 30/11/00

      if(ntest.gt.4)then
!       diag check for unstable layers
        do iq=1,ifull
         do k=1,kl
          theta(k)=tx(iq,k)*sig(k)**(-r/cp)
         enddo
         do k=ntest,kl   ! e.g. 8,kl
          if(theta(k).lt.theta(k-1))then  ! based on tx
            print *,"unstable layer in upglobal for ktau,iq,k's,del ",
     .                               ktau,iq,k-1,k,theta(k-1)-theta(k)
	     write (6,"('theta',9f7.2/5x,9f7.2)") theta
            write (6,"('sdot',9f7.3/4x,9f7.3)") (sdot(iq,kk),kk=1,kl)
c           print *,'sdot', (sdot(iq,kk),kk=1,kl)
            numunstab=numunstab+1
c           if(numunstab.eq.100)stop 'numunstab=30'
          endif
         enddo
        enddo
      endif

      return
      end
