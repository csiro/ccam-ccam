      subroutine vadvtvd(tarr,uarr,varr,nvadh_pass)   ! globpea  version
c                              vadvbott & vadvyu at bottom
!     can show adding tbar has no effect
      use cc_mpi, only : mydiag, myid
      include 'newmpar.h'
      parameter (npslx=1)  ! 0 off, 1 on for nvad=-4
      parameter (nqq=0)    ! 0 off, 3 possible
!     parameter (sdotfilt=0.)  ! tried .3 - not useful
!     parameter (nimp=1)  !  0 for original explicit non-flux TVD term
!                            1 for implicit non-flux TVD term
!     parameter (nthub=1) !  1 original
!                            2 higher-order is Lax-Wendroff
!                            0 (not available now) was equivalent to original
!     parameter (ntvd=2)  !  1 van Leer phi flux-limiter (original)
!                            2 MC phi flux-limiter
!                            3 superbee flux-limiter
!                            0 (not available now) was equivalent to van Leer
c     split vertical advection routine; tvd scheme; used with nonlin or upglobal
c     In flux limiter, assuming zero gradient for all top and bottom
c     variables; except extrap at bottom for qg and trace gases  Thu  06-19-1997
      include 'arrays.h'
      include 'kuocom.h'     ! also with kbsav,ktsav
      include 'liqwpar.h'  ! ifullw
      include 'map.h'
      include 'parm.h'
      include 'parmdyn.h'
      include 'parmvert.h' ! nthub,nimp,ntvd
      include 'sigs.h'
      include 'tracers.h'
      include 'vvel.h'
      include 'xarrs.h'
c     common/work2/aa(ifull),dum2(ifull,17) 
!     N.B. first 3 arrays of work3 are available from nonlin, all from upglobal
      common/work3/delt(ifull,0:kl),fluxh(ifull,0:kl),
     .             dum3(3*ijk-2*ifull)
      common/work3b/xin(ifull,kl),uav(ifull,kl)
      real udiff(ifull,kl)
      equivalence (udiff,delt)  ! udiff just used for qg, before delt
      real tarr(ifull,kl),uarr(ifull,kl),varr(ifull,kl)
      real sig3(kl)
      integer num,nsign
      data num/0/,nsign/1/
      save num,nsign   ! -1,1, -1,1 and so on

      tfact=1./nvadh_pass 

      if(num.eq.0)then
        num=1
        if(ntvdr.eq.0)then  ! to produce old nvad=4 interps for fluxhi
	   ratha(:)=.5
	   rathb(:)=.5
        endif
        if(myid==0)then
          print *,'In vadvtvd nvad,nvadh_pass,nqq,npslx,ntvdr ',
     .                        nvad,nvadh_pass,nqq,npslx,ntvdr
          print *,'nimp,nthub,ntvd,tfact ',nimp,nthub,ntvd,tfact
        endif
      endif

c     note sdot coming through is at level k-.5
c     converted in nonlin to units of grid-steps/timestep,  +ve upwards
      do k=1,kl
       sig3(k)=sig(k)**3
      enddo     ! k loop

      do iq=1,ifull
c      fluxh(k) is located at level k+.5
       fluxh(iq,0)=0.
       fluxh(iq,kl)=0.
       delt(iq,0)=0.      ! for T,u,v
       delt(iq,kl)=0.     ! for T,u,v
      enddo    ! iq loop

c     t
      if(ntvdr.eq.2)xin(:,:)=tarr(1:ifull,:)
      do k=1,kl-1
       do iq=1,ifull
         delt(iq,k)=tarr(iq,k+1)-tarr(iq,k)
       enddo    ! iq loop
      enddo     ! k loop
      do k=1,kl-1  ! for fluxh at interior (k + 1/2)  half-levels
       do iq=1,ifull
        kp=sign(1.,sdot(iq,k+1))
        kx=k+(1-kp)/2  !  k for sdot +ve,  k+1 for sdot -ve
        rat=delt(iq,k-kp)/(delt(iq,k)+sign(1.e-20,delt(iq,k)))
        if(ntvd.eq.1)phi=(rat+abs(rat))/(1.+abs(rat))       ! 0 for -ve rat
        if(ntvd.eq.2)phi=max(0.,min(2.*rat,.5+.5*rat,2.))   ! 0 for -ve rat
        if(ntvd.eq.3)phi=max(0.,min(1.,2.*rat),min(2.,rat)) ! 0 for -ve rat
        if(nthub.eq.1)fluxhi=rathb(k)*tarr(iq,k)+ratha(k)*tarr(iq,k+1)
        if(nthub.eq.2)then     ! higher order scheme
          fluxhi=rathb(k)*tarr(iq,k)+ratha(k)*tarr(iq,k+1)
     .                  -.5*delt(iq,k)*tfact*sdot(iq,k)
        endif  ! (nthub.eq.2)
        fluxlo=tarr(iq,kx)
        fluxh(iq,k)=sdot(iq,k+1)*(fluxlo+phi*(fluxhi-fluxlo))
c       if(iq.eq.idjd)print *
c       if(iq.eq.idjd)print *,'k,nthub,ntvd,sdot(k+1) ',
c    .                         k,nthub,ntvd,sdot(iq,k+1)
c       if(iq.eq.idjd)print *,'delt -,.,+ '
c    .            ,delt(iq,k-1),delt(iq,k),delt(iq,k+1)
c       if(iq.eq.idjd)print *,'rat,phi,tin ',rat,phi,tarr(iq,kx)
c       if(iq.eq.idjd)print *,'kp,kx,fluxlo,fluxhi,fluxh ',
c    .                         kp,kx,fluxlo,fluxhi,fluxh(iq,k)
       enddo    ! iq loop
      enddo     ! k loop
      do k=1,kl
       do iq=1,ifull
        if(nimp.eq.1)then
!        N.B. nimp=1 gives especially silly results if sdot>1	 
         hdsdot=.5*tfact*(sdot(iq,k+1)-sdot(iq,k))
         tarr(iq,k)=(tarr(iq,k)+tfact*(fluxh(iq,k-1)-fluxh(iq,k))
     .               +hdsdot*tarr(iq,k) )/(1.-hdsdot)
        else
         tarr(iq,k)=tarr(iq,k)+tfact*(fluxh(iq,k-1)-fluxh(iq,k)
     .               +tarr(iq,k)*(sdot(iq,k+1)-sdot(iq,k)))
        endif   ! (nimp.eq.1)
       enddo    ! iq loop
      enddo     ! k loop
!!    tarr(1:ifull,1)=max(xin(:,1),min(tarr(1:ifull,1),xin(:,2)))
!!    tarr(1:ifull,kl)=max(xin(:,kl),min(tarr(1:ifull,kl),xin(:,kl-1)))
      if(ntvdr.eq.2)then  ! different top & bottom
        do iq=1,ifull
         if(sdot(iq,2).gt.0.)tarr(iq,1)=xin(iq,1)
         if(sdot(iq,kl).lt.0.)tarr(iq,kl)=xin(iq,kl)
        enddo
        xin(:,:)=uarr(1:ifull,:)
      endif      ! (ntvdr.eq.2)
      if( (diag.or.nmaxpr.eq.1) .and. mydiag )then
!        These diagnostics don't work with single input/output argument
!        write (6,"('tin ',9f8.2/4x,9f8.2)") (tin(idjd,k),k=1,kl)
!        write (6,"('tadd',9f8.2/4x,9f8.2)") 
!     .        (tfact*(tin(idjd,k)*(sdot(idjd,k+1)-sdot(idjd,k))),k=1,kl)
        write (6,"('tout',9f8.2/4x,9f8.2)") (tarr(idjd,k),k=1,kl)
      endif

c     u
      do k=1,kl-1
       do iq=1,ifull
        delt(iq,k)=uarr(iq,k+1)-uarr(iq,k)
       enddo    ! iq loop
      enddo     ! k loop
      do k=1,kl-1  ! for fluxh at interior (k + 1/2)  half-levels
       do iq=1,ifull
        kp=sign(1.,sdot(iq,k+1))
        kx=k+(1-kp)/2  !  k for sdot +ve,  k+1 for sdot -ve
        rat=delt(iq,k-kp)/(delt(iq,k)+sign(1.e-20,delt(iq,k)))
        if(ntvd.eq.1)phi=(rat+abs(rat))/(1.+abs(rat))       ! 0 for -ve rat
        if(ntvd.eq.2)phi=max(0.,min(2.*rat,.5+.5*rat,2.))   ! 0 for -ve rat
        if(ntvd.eq.3)phi=max(0.,min(1.,2.*rat),min(2.,rat)) ! 0 for -ve rat
        if(nthub.eq.1)fluxhi=rathb(k)*uarr(iq,k)+ratha(k)*uarr(iq,k+1)
        if(nthub.eq.2)then     ! higher order scheme
          fluxhi=rathb(k)*uarr(iq,k)+ratha(k)*uarr(iq,k+1)
     .                  -.5*delt(iq,k)*tfact*sdot(iq,k)
        endif  ! (nthub.eq.2)
        fluxlo=uarr(iq,kx)
        fluxh(iq,k)=sdot(iq,k+1)*(fluxlo+phi*(fluxhi-fluxlo))
       enddo    ! iq loop
      enddo     ! k loop
      do k=1,kl
       do iq=1,ifull
        if(nimp.eq.1)then
         hdsdot=.5*tfact*(sdot(iq,k+1)-sdot(iq,k))
         uarr(iq,k)=(uarr(iq,k)+tfact*(fluxh(iq,k-1)-fluxh(iq,k))
     .               +hdsdot*uarr(iq,k) )/(1.-hdsdot)
        else
         uarr(iq,k)=uarr(iq,k)+tfact*(fluxh(iq,k-1)-fluxh(iq,k)
     .               +uarr(iq,k)*(sdot(iq,k+1)-sdot(iq,k)))
        endif   ! (nimp.eq.1)
       enddo    ! iq loop
      enddo     ! k loop
!!    uarr(1:ifull,1)=max(xin(:,1),min(uarr(1:ifull,1),xin(:,2)))
!!    uarr(1:ifull,kl)=max(xin(:,kl),min(uarr(1:ifull,kl),xin(:,kl-1)))
      if(ntvdr.eq.2)then  ! different top & bottom
        do iq=1,ifull
         if(sdot(iq,2).gt.0.)uarr(iq,1)=xin(iq,1)
         if(sdot(iq,kl).lt.0.)uarr(iq,kl)=xin(iq,kl)
        enddo
        xin(:,:)=varr(1:ifull,:)
      endif      ! (ntvdr.eq.2)

c     v
      do k=1,kl-1
       do iq=1,ifull
         delt(iq,k)=varr(iq,k+1)-varr(iq,k)
       enddo    ! iq loop
      enddo     ! k loop
      do k=1,kl-1  ! for fluxh at interior (k + 1/2)  half-levels
       do iq=1,ifull
        kp=sign(1.,sdot(iq,k+1))
        kx=k+(1-kp)/2  !  k for sdot +ve,  k+1 for sdot -ve
        rat=delt(iq,k-kp)/(delt(iq,k)+sign(1.e-20,delt(iq,k)))
        if(ntvd.eq.1)phi=(rat+abs(rat))/(1.+abs(rat))       ! 0 for -ve rat
        if(ntvd.eq.2)phi=max(0.,min(2.*rat,.5+.5*rat,2.))   ! 0 for -ve rat
        if(ntvd.eq.3)phi=max(0.,min(1.,2.*rat),min(2.,rat)) ! 0 for -ve rat
        if(nthub.eq.1)fluxhi=rathb(k)*varr(iq,k)+ratha(k)*varr(iq,k+1)
        if(nthub.eq.2)then     ! higher order scheme
          fluxhi=rathb(k)*varr(iq,k)+ratha(k)*varr(iq,k+1)
     .                  -.5*delt(iq,k)*tfact*sdot(iq,k)
        endif  ! (nthub.eq.2)
        fluxlo=varr(iq,kx)
        fluxh(iq,k)=sdot(iq,k+1)*(fluxlo+phi*(fluxhi-fluxlo))
       enddo    ! iq loop
      enddo     ! k loop
      do k=1,kl
       do iq=1,ifull
        if(nimp.eq.1)then
         hdsdot=.5*tfact*(sdot(iq,k+1)-sdot(iq,k))
         varr(iq,k)=(varr(iq,k)+tfact*(fluxh(iq,k-1)-fluxh(iq,k))
     .               +hdsdot*varr(iq,k) )/(1.-hdsdot)
        else
         varr(iq,k)=varr(iq,k)+tfact*(fluxh(iq,k-1)-fluxh(iq,k)
     .               +varr(iq,k)*(sdot(iq,k+1)-sdot(iq,k)))
        endif   ! (nimp.eq.1)
       enddo    ! iq loop
      enddo     ! k loop
!!    varr(1:ifull,1)=max(xin(:,1),min(varr(1:ifull,1),xin(:,2)))
!!    varr(1:ifull,kl)=max(xin(:,kl),min(varr(1:ifull,kl),xin(:,kl-1)))
      if(ntvdr.eq.2)then  ! different top & bottom
        do iq=1,ifull
         if(sdot(iq,2).gt.0.)varr(iq,1)=xin(iq,1)
         if(sdot(iq,kl).lt.0.)varr(iq,kl)=xin(iq,kl)
        enddo
        xin(:,:)=qg(1:ifull,:)
      endif      ! (ntvdr.eq.2)

      if(npslx.eq.1.and.nvad.le.-4)then  ! handles -9 too
      do k=1,kl-1
       do iq=1,ifull
         delt(iq,k)=pslx(iq,k+1)-pslx(iq,k)
       enddo    ! iq loop
      enddo     ! k loop
      do k=1,kl-1  ! for fluxh at interior (k + 1/2)  half-levels
       do iq=1,ifull
        kp=sign(1.,sdot(iq,k+1))
        kx=k+(1-kp)/2  !  k for sdot +ve,  k+1 for sdot -ve
        rat=delt(iq,k-kp)/(delt(iq,k)+sign(1.e-20,delt(iq,k)))
        if(ntvd.eq.1)phi=(rat+abs(rat))/(1.+abs(rat))       ! 0 for -ve rat
        if(ntvd.eq.2)phi=max(0.,min(2.*rat,.5+.5*rat,2.))   ! 0 for -ve rat
        if(ntvd.eq.3)phi=max(0.,min(1.,2.*rat),min(2.,rat)) ! 0 for -ve rat
        if(nthub.eq.1)fluxhi=rathb(k)*pslx(iq,k)+ratha(k)*pslx(iq,k+1)
        if(nthub.eq.2)then     ! higher order scheme
          fluxhi=rathb(k)*pslx(iq,k)+ratha(k)*pslx(iq,k+1)
     .                  -.5*delt(iq,k)*tfact*sdot(iq,k)
        endif  ! (nthub.eq.2)
        fluxlo=pslx(iq,kx)
        fluxh(iq,k)=sdot(iq,k+1)*(fluxlo+phi*(fluxhi-fluxlo))
       enddo    ! iq loop
      enddo     ! k loop
      do k=1,kl
       do iq=1,ifull
        if(nimp.eq.1)then
         hdsdot=.5*tfact*(sdot(iq,k+1)-sdot(iq,k))
         pslx(iq,k)=(pslx(iq,k)+tfact*(fluxh(iq,k-1)-fluxh(iq,k))
     .               +hdsdot*pslx(iq,k) )/(1.-hdsdot)
        else
         pslx(iq,k)=pslx(iq,k)+tfact*(fluxh(iq,k-1)-fluxh(iq,k)
     .               +pslx(iq,k)*(sdot(iq,k+1)-sdot(iq,k)))
        endif   ! (nimp.eq.1)
       enddo    ! iq loop
      enddo     ! k loop
      endif  ! (npslx.eq.1.and.nvad.eq.-4)

      if(mspec.eq.1.and.abs(nvad).ne.9)then   ! advect qg and gases after preliminary step
c     qg
        if(nqq.eq.3)then
!         qg(:,:)=cbrt(qg(:,:))
          qg(1:ifull,:)=qg(1:ifull,:)**(1./3.)
        endif       ! (nqq.eq.3)
      do k=1,kl-1
       do iq=1,ifull
         delt(iq,k)=qg(iq,k+1)-qg(iq,k)
         enddo    ! iq loop
        enddo     ! k loop
        do iq=1,ifull
         delt(iq,0)=min(delt(iq,1),qg(iq,1))        ! for non-negative tt
!        delt(iq,kl)=min(delt(iq,kl-1),qg(iq,kl))   ! for non-negative tt
         delt(iq,kl)=0.                       ! safer
        enddo    ! iq loop
        do k=1,kl-1  ! for fluxh at interior (k + 1/2)  half-levels
         do iq=1,ifull
        kp=sign(1.,sdot(iq,k+1))
        kx=k+(1-kp)/2  !  k for sdot +ve,  k+1 for sdot -ve
        rat=delt(iq,k-kp)/(delt(iq,k)+sign(1.e-20,delt(iq,k)))
        if(ntvd.eq.1)phi=(rat+abs(rat))/(1.+abs(rat))       ! 0 for -ve rat
        if(ntvd.eq.2)phi=max(0.,min(2.*rat,.5+.5*rat,2.))   ! 0 for -ve rat
        if(ntvd.eq.3)phi=max(0.,min(1.,2.*rat),min(2.,rat)) ! 0 for -ve rat
        if(nthub.eq.1)fluxhi=rathb(k)*qg(iq,k)+ratha(k)*qg(iq,k+1)
        if(nthub.eq.2)then     ! higher order scheme
          fluxhi=rathb(k)*qg(iq,k)+ratha(k)*qg(iq,k+1)
     .                  -.5*delt(iq,k)*tfact*sdot(iq,k)
        endif  ! (nthub.eq.2)
        fluxlo=qg(iq,kx)
        fluxh(iq,k)=sdot(iq,k+1)*(fluxlo+phi*(fluxhi-fluxlo))
         enddo   ! iq loop
        enddo    ! k loop
        do k=1,kl   !  new option July 2001
         do iq=1,ifull
        if(nimp.eq.1)then
         hdsdot=.5*tfact*(sdot(iq,k+1)-sdot(iq,k))
         qg(iq,k)=(qg(iq,k)+tfact*(fluxh(iq,k-1)-fluxh(iq,k))
     .               +hdsdot*qg(iq,k) )/(1.-hdsdot)
        else
         qg(iq,k)=qg(iq,k)+tfact*(fluxh(iq,k-1)-fluxh(iq,k)
     .               +qg(iq,k)*(sdot(iq,k+1)-sdot(iq,k)))
        endif   ! (nimp.eq.1)
         enddo   ! iq loop
        enddo    ! k loop
!       some time check out following new option July 2001
c        k=1
c        do iq=1,ifull
c         deltaqg=tfact*(fluxh(iq,k-1)-fluxh(iq,k)
c     .               +qg(iq,k)*(sdot(iq,k+1)-sdot(iq,k)))
c         if(sdot(iq,k+1).lt.0.)then
c            if(sign(1.,deltaqg).eq.sign(1.,qg(iq,k+1)-qg(iq,k)))
c     .                           qg(iq,k)=qg(iq,k)+deltaqg
c         else
c           qg(iq,k)=qg(iq,k)+deltaqg
c         endif
c        enddo   ! iq loop
c        k=kl
c        do iq=1,ifull
c         deltaqg=tfact*(fluxh(iq,k-1)-fluxh(iq,k)
c     .               +qg(iq,k)*(sdot(iq,k+1)-sdot(iq,k)))
c         if(sdot(iq,k).gt.0.)then
c            if(sign(1.,deltaqg).eq.sign(1.,qg(iq,k-1)-qg(iq,k)))
c     .                           qg(iq,k)=qg(iq,k)+deltaqg
c         else
c           qg(iq,k)=qg(iq,k)+deltaqg
c         endif
c        enddo   ! iq loop
        if(nqq.eq.3)then
          qg(1:ifull,:)=qg(1:ifull,:)**3
        endif       ! (nqq.eq.3)
!!    qg(1:ifull,1)=max(xin(:,1),min(qg(1:ifull,1),xin(:,2)))
!!    qg(1:ifull,kl)=max(xin(:,kl),min(qg(1:ifull,kl),xin(:,kl-1)))
      if(ntvdr.eq.2)then  ! different top & bottom
        do iq=1,ifull
         if(sdot(iq,2).gt.0.)qg(iq,1)=xin(iq,1)
         if(sdot(iq,kl).lt.0.)qg(iq,kl)=xin(iq,kl)
        enddo
        xin(:,:)=qlg(1:ifull,:)
      endif      ! (ntvdr.eq.2)

      if(ldr.ne.0)then
       do k=1,kl-1       ! qlg first
        do iq=1,ifull
          delt(iq,k)=qlg(iq,k+1)-qlg(iq,k)
        enddo    ! iq loop
       enddo     ! k loop
       do iq=1,ifull
        delt(iq,0)=min(delt(iq,1),qlg(iq,1))       ! for non-negative tt
       enddo    ! iq loop
       do k=1,kl-1  ! for fluxh at interior (k + 1/2)  half-levels
        do iq=1,ifull
         kp=sign(1.,sdot(iq,k+1))
         kx=k+(1-kp)/2  !  k for sdot +ve,  k+1 for sdot -ve
         rat=delt(iq,k-kp)/(delt(iq,k)+sign(1.e-20,delt(iq,k)))
         if(ntvd.eq.1)phi=(rat+abs(rat))/(1.+abs(rat))       ! 0 for -ve rat
         if(ntvd.eq.2)phi=max(0.,min(2.*rat,.5+.5*rat,2.))   ! 0 for -ve rat
         if(ntvd.eq.3)phi=max(0.,min(1.,2.*rat),min(2.,rat)) ! 0 for -ve rat
         if(nthub.eq.1)fluxhi=rathb(k)*qlg(iq,k)+ratha(k)*qlg(iq,k+1)
         if(nthub.eq.2)then     ! higher order scheme
           fluxhi=rathb(k)*qlg(iq,k)+ratha(k)*qlg(iq,k+1)
     .                   -.5*delt(iq,k)*tfact*sdot(iq,k)
         endif  ! (nthub.eq.2)
         fluxlo=qlg(iq,kx)
         fluxh(iq,k)=sdot(iq,k+1)*(fluxlo+phi*(fluxhi-fluxlo))
        enddo    ! iq loop
       enddo     ! k loop
       do k=1,kl
        do iq=1,ifull
         if(nimp.eq.1)then
          hdsdot=.5*tfact*(sdot(iq,k+1)-sdot(iq,k))
          qlg(iq,k)=(qlg(iq,k)
     .                +tfact*(fluxh(iq,k-1)-fluxh(iq,k))
     .                +hdsdot*qlg(iq,k) )/(1.-hdsdot)
         else
          qlg(iq,k)=qlg(iq,k)
     .                +tfact*(fluxh(iq,k-1)-fluxh(iq,k)
     .                +qlg(iq,k)*(sdot(iq,k+1)-sdot(iq,k)))
         endif   ! (nimp.eq.1)
        enddo     ! iq loop
       enddo      ! k loop
!!     qlg(1:ifull,1)=max(xin(:,1),min(qlg(1:ifull,1),xin(:,2)))
!!     qlg(1:ifull,kl)=max(xin(:,kl),min(qlg(1:ifull,kl),xin(:,kl-1)))
       if(ntvdr.eq.2)then  ! different top & bottom
         do iq=1,ifull
          if(sdot(iq,2).gt.0.)qlg(iq,1)=xin(iq,1)
          if(sdot(iq,kl).lt.0.)qlg(iq,kl)=xin(iq,kl)
         enddo
         xin(:,:)=qfg(1:ifull,:)
       endif      ! (ntvdr.eq.2)
       do k=1,kl-1       ! qfg next
        do iq=1,ifull
          delt(iq,k)=qfg(iq,k+1)-qfg(iq,k)
        enddo    ! iq loop
       enddo     ! k loop
       do iq=1,ifull
        delt(iq,0)=min(delt(iq,1),qfg(iq,1))       ! for non-negative tt
       enddo    ! iq loop
       do k=1,kl-1  ! for fluxh at interior (k + 1/2)  half-levels
        do iq=1,ifull
         kp=sign(1.,sdot(iq,k+1))
         kx=k+(1-kp)/2  !  k for sdot +ve,  k+1 for sdot -ve
         rat=delt(iq,k-kp)/(delt(iq,k)+sign(1.e-20,delt(iq,k)))
         if(ntvd.eq.1)phi=(rat+abs(rat))/(1.+abs(rat))       ! 0 for -ve rat
         if(ntvd.eq.2)phi=max(0.,min(2.*rat,.5+.5*rat,2.))   ! 0 for -ve rat
         if(ntvd.eq.3)phi=max(0.,min(1.,2.*rat),min(2.,rat)) ! 0 for -ve rat
         if(nthub.eq.1)fluxhi=rathb(k)*qfg(iq,k)+ratha(k)*qfg(iq,k+1)
         if(nthub.eq.2)then     ! higher order scheme
           fluxhi=rathb(k)*qfg(iq,k)+ratha(k)*qfg(iq,k+1)
     .                   -.5*delt(iq,k)*tfact*sdot(iq,k)
         endif  ! (nthub.eq.2)
         fluxlo=qfg(iq,kx)
         fluxh(iq,k)=sdot(iq,k+1)*(fluxlo+phi*(fluxhi-fluxlo))
        enddo    ! iq loop
       enddo     ! k loop
       do k=1,kl
        do iq=1,ifull
         if(nimp.eq.1)then
          hdsdot=.5*tfact*(sdot(iq,k+1)-sdot(iq,k))
          qfg(iq,k)=(qfg(iq,k)
     .                +tfact*(fluxh(iq,k-1)-fluxh(iq,k))
     .                +hdsdot*qfg(iq,k) )/(1.-hdsdot)
         else
          qfg(iq,k)=qfg(iq,k)
     .                +tfact*(fluxh(iq,k-1)-fluxh(iq,k)
     .                +qfg(iq,k)*(sdot(iq,k+1)-sdot(iq,k)))
         endif   ! (nimp.eq.1)
        enddo     ! iq loop
       enddo      ! k loop
!!     qfg(1:ifull,1)=max(xin(:,1),min(qfg(1:ifull,1),xin(:,2)))
!!     qfg(1:ifull,kl)=max(xin(:,kl),min(qfg(1:ifull,kl),xin(:,kl-1)))
       if(ntvdr.eq.2)then  ! different top & bottom
         do iq=1,ifull
          if(sdot(iq,2).gt.0.)qfg(iq,1)=xin(iq,1)
          if(sdot(iq,kl).lt.0.)qfg(iq,kl)=xin(iq,kl)
         enddo
       endif      ! (ntvdr.eq.2)
      endif      ! if(ldr.ne.0)

      if(ilt.gt.1)then
      do ntr=1,ntrac
      do k=1,kl-1
       do iq=1,ifull
         delt(iq,k)=tr(iq,k+1,ntr)-tr(iq,k,ntr)
       enddo    ! iq loop
      enddo     ! k loop
      do iq=1,ifull
       delt(iq,0)=min(delt(iq,1),tr(iq,1,ntr))       ! for non-negative tt
      enddo    ! iq loop
      do k=1,kl-1  ! for fluxh at interior (k + 1/2)  half-levels
       do iq=1,ifull
        kp=sign(1.,sdot(iq,k+1))
        kx=k+(1-kp)/2  !  k for sdot +ve,  k+1 for sdot -ve
        rat=delt(iq,k-kp)/(delt(iq,k)+sign(1.e-20,delt(iq,k)))
        if(ntvd.eq.1)phi=(rat+abs(rat))/(1.+abs(rat))       ! 0 for -ve rat
        if(ntvd.eq.2)phi=max(0.,min(2.*rat,.5+.5*rat,2.))   ! 0 for -ve rat
        if(ntvd.eq.3)phi=max(0.,min(1.,2.*rat),min(2.,rat)) ! 0 for -ve rat
        if(nthub.eq.1)fluxhi=rathb(k)*tr(iq,k,ntr)+
     .                                    ratha(k)*tr(iq,k+1,ntr)
        if(nthub.eq.2)then     ! higher order scheme
          fluxhi=rathb(k)*tr(iq,k,ntr)+ratha(k)*tr(iq,k+1,ntr)
     .                  -.5*delt(iq,k)*tfact*sdot(iq,k)
        endif  ! (nthub.eq.2)
        fluxlo=tr(iq,kx,ntr)
        fluxh(iq,k)=sdot(iq,k+1)*(fluxlo+phi*(fluxhi-fluxlo))
       enddo    ! iq loop
      enddo     ! k loop
      do k=1,kl
       do iq=1,ifull
        if(nimp.eq.1)then
         hdsdot=.5*tfact*(sdot(iq,k+1)-sdot(iq,k))
         tr(iq,k,ntr)=(tr(iq,k,ntr)
     .               +tfact*(fluxh(iq,k-1)-fluxh(iq,k))
     .               +hdsdot*tr(iq,k,ntr) )/(1.-hdsdot)
        else
         tr(iq,k,ntr)=tr(iq,k,ntr)
     .               +tfact*(fluxh(iq,k-1)-fluxh(iq,k)
     .               +tr(iq,k,ntr)*(sdot(iq,k+1)-sdot(iq,k)))
        endif   ! (nimp.eq.1)
       enddo     ! iq loop
      enddo      ! k loop
      enddo      ! ntr loop
      endif      ! if(ilt.gt.1)

      endif       ! if(mspec.eq.1.and.abs(nvad).ne.9)

      return
      end
