      program tvdtest
!     N.B. uin and uout may share same storage (similarly tin, vin)
!       nqq:     0 for simple qg; 2 for k**3 treatment; 4 for qg**(1/4)
!                7 for old nqq=3 (up to July '98)
!     This is vector version,with sdotfilt>0 removed    [not yet done for qg]
!       nsord:   1 (usual); 3 - gave Antarctic easterly strato winds for l_w!)
      parameter (il=1,jl=1,kl=9,ifull=il*jl,ijk=ifull*kl)
!     parameter (sdotfilt=0.)  ! tried .3 - not useful
      parameter (nimp=1)  !  0 for original explicit non-flux TVD term
!                            1 for implicit non-flux TVD term
      parameter (nthub=2) !  0 original
!                            1 equivalent to original
!                            2 higher-order is Lax-Wendroff
      parameter (ntvd=2)  !  0 original phi flux-limiter
!                            1 redefined original phi (same answers as 0)
!                            2 MC phi flux-limiter
c     split vertical advection routine; tvd scheme; used with nonlin or upglobal
c     In flux limiter, assuming zero gradient for all top and bottom
c     variables; except extrap at bottom for qg and trace gases  Thu  06-19-1997

      common/work2/aa(il,jl),dum2(il,jl,17) 
!     N.B. first 3 arrays of work3 are available from nonlin, all from upglobal
      common/work3/delt(il,jl,0:kl),fluxh(il,jl,0:kl),
     .             dum3(2*ijk-2*ifull)
      common/work3b/dum3b(il,jl,kl),uav(il,jl,kl)
      real udiff(il,jl,kl)
      equivalence (udiff,delt)  ! udiff just used for qg, before delt
      real tin(il,jl,kl),tout(il,jl,kl),uin(il,jl,kl),uout(il,jl,kl),
     .     vin(il,jl,kl),vout(il,jl,kl),sdot(il,jl,kl+1)
      real sig3(kl)
      data sdot/0.,1.780,.480,.219,.138,.109,.089,.074,.060,.051/
      data tin/273.10,289.94,292.28,293.04,290.48,286.74,283.05,
     .         278.53,272.18/
      data num/0/,nsign/1/
      save num,nsign   ! -1,1, -1,1 and so on

      nvadh=2
      id=1
      jd=1
      idjd=1
      print *,'tin  ',(tin(id,jd,k),k=1,kl)
      print *,'sdot  ',(sdot(id,jd,k),k=1,kl+1)
c     if(nvadh.eq.2)then
c       nsign=-nsign      !  -1,1, -1,1 and so on
c       if(nsign.lt.0)then
c         tfact=.5*(1.-epsp)
c       else
c        tfact=.5*(1.+epsp)
c       endif
c     else       ! i.e. nvadh.ne.2
c       tfact=1.
c     endif
      tfact=1./nvadh   ! simpler alternative

      if(num.eq.0)then
        num=1
        print *,'In vadvtvd nvad,nvadh,nqq,nimp,nthub,ntvd,tfact ',
     .                      nvad,nvadh,nqq,nimp,nthub,ntvd,tfact
      endif

c     note sdot coming through is at level k-.5
c     converted in nonlin to units of grid-steps/timestep,  +ve upwards


      do iq=1,ifull
c      fluxh(k) is located at level k+.5
       fluxh(iq,1,0)=0.
       fluxh(iq,1,kl)=0.
       delt(iq,1,0)=0.      ! for u,v
       delt(iq,1,kl)=0.     ! for u,v
      enddo    ! iq loop

c     t
      do k=1,kl-1
       do iq=1,ifull
         delt(iq,1,k)=tin(iq,1,k+1)-tin(iq,1,k)
       enddo    ! iq loop
      enddo     ! k loop
      do k=1,kl-1  ! for fluxh at interior (k + 1/2)  half-levels
       do iq=1,ifull
        kp=sign(1.,sdot(iq,1,k+1))
        kx=k+(1-kp)/2  !  k for sdot +ve,  k+1 for sdot -ve
        if(ntvd.eq.0)
     .  phi=abs( sign(delt(iq,1,k-kp),delt(iq,1,k))+delt(iq,1,k-kp) )/
     .      (abs(delt(iq,1,k)+delt(iq,1,k-kp))+1.e-20)
        if(ntvd.gt.0)then
          rat=delt(iq,1,k-kp)/(delt(iq,1,k)+sign(1.e-20,delt(iq,1,k)))
          if(ntvd.eq.1)phi=(rat+abs(rat))/(1.+abs(rat))
          if(ntvd.eq.2)phi=max(0.,min(2.*rat,.5+.5*rat,2.))
        endif
        if(nthub.eq.0)then
          fluxh(iq,1,k)=sdot(iq,1,k+1)*
     .        (tin(iq,1,kx)+phi*delt(iq,1,k)*sign(.5,sdot(iq,1,k+1)))
        endif  ! (nthub.eq.0)
        if(nthub.eq.1)then
          fluxlo=tin(iq,1,kx)
          fluxhi=.5*(tin(iq,1,k)+tin(iq,1,k+1))
          fluxh(iq,1,k)=sdot(iq,1,k+1)*(fluxlo+phi*(fluxhi-fluxlo))
        endif  ! (nthub.eq.1)
        if(nthub.eq.2)then     ! higher order scheme
          fluxlo=tin(iq,1,kx)
          fluxhi=.5*(tin(iq,1,k)+tin(iq,1,k+1)
     .                  -delt(iq,1,k)*tfact*sdot(iq,1,k))
          fluxh(iq,1,k)=sdot(iq,1,k+1)*(fluxlo+phi*(fluxhi-fluxlo))
        endif  ! (nthub.eq.2)
        hdsdot=.5*tfact*(sdot(iq,1,k+1)-sdot(iq,1,k))  ! for diags only
        if(iq.eq.idjd)print *
        if(iq.eq.idjd)print *,'k,nthub,ntvd,sdot(k+1) ',
     .                         k,nthub,ntvd,sdot(iq,1,k+1)
        if(iq.eq.idjd)print *,'delt -,.,+ '
     .            ,delt(iq,1,k-1),delt(iq,1,k),delt(iq,1,k+1)
        if(iq.eq.idjd)print *,'rat,phi,tfact,hdsdot,tin ',
     .                         rat,phi,tfact,hdsdot,tin(iq,1,kx)
        if(iq.eq.idjd)print *,'kp,kx,fluxlo,fluxhi,fluxh ',
     .                         kp,kx,fluxlo,fluxhi,fluxh(iq,1,k)
       enddo    ! iq loop
      enddo     ! k loop
      do k=1,kl
       do iq=1,ifull
        if(nimp.eq.1)then
         hdsdot=.5*tfact*(sdot(iq,1,k+1)-sdot(iq,1,k))
         tout(iq,1,k)=(tin(iq,1,k)+tfact*(fluxh(iq,1,k-1)-fluxh(iq,1,k))
     .               +hdsdot*tin(iq,1,k) )/(1.-hdsdot)
        else
         tout(iq,1,k)=tin(iq,1,k)+tfact*(fluxh(iq,1,k-1)-fluxh(iq,1,k)
     .               +tin(iq,1,k)*(sdot(iq,1,k+1)-sdot(iq,1,k)))
        endif   ! (nimp.eq.1)
       enddo    ! iq loop
      enddo     ! k loop
      print *,'tin  ',(tin(id,jd,k),k=1,kl)
      print *,'tout ',(tout(id,jd,k),k=1,kl)
c     print *,'tdiff ',(tout(id,jd,k)-tin(id,jd,k),k=1,kl)

      end

