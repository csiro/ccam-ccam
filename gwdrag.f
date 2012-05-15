      subroutine gwdrag   ! globpea/darlam (but not staggered)
!     this is vectorized jlm version
      use arrays_m
      use gdrag_m
      use morepbl_m
      use nlin_m
      use pbl_m
      use sigs_m
      use soil_m
!     parameter(fc2=.5,ndzx=0)    !  as per Hal      gwd1a
      parameter(fc2=1.,ndzx=1)    !  as per jlm      gwd1b
      parameter(ntest=0)          ! ntest= 0 for diags off; ntest= 1 for diags on
      include 'newmpar.h'
      include 'const_phys.h'
      include 'parm.h'
      real uu(ifull,kl),fni(ifull,kl),bvnf(ifull,kl)
     .            ,thf(ifull,kl)
      real dthdz(ifull,kl)
      real bvng(ifull),temp(ifull),fnii(ifull),unst(ifull)
     . ,apuw(ifull),apvw(ifull),alam(ifull),wmag(ifull),frsav(ifull)
      real delt(kl),dsk(kl),sigk(kl)
c     real fr2vsav(ifull,kl),uusav(ifull,kl)   ! jlm diag
c     interface mods for darlam & globpe
      do k=1,kl
       dsk(k)=-dsig(k)
       sigk(k)=sig(k)**(rdry/cp)
       delt(k)=1.
      enddo
      delt(kl)=0.

      do k=1,kl
       do iq=1,ifull
c       put theta in thf()
        thf(iq,k)=t(iq,k)/sigk(k)                ! gwdrag
       enddo   ! iq loop
      enddo    ! k loop

c     calc d(theta)/dz  at half-levels , using 1/dz at level k-.5
      if(ndzx.eq.0)dzx=-grav*sig(1)/(dsig(1)*rdry)                ! Hal's
      if(ndzx.eq.1)dzx=.5*grav*(1.+sig(1))/((1.-sig(1))*rdry)     ! fixup by jlm
      do iq=1,ifull
       dzi=dzx/t(iq,1)
       dthdz(iq,1)=max(thf(iq,1)-tss(iq),0.)*dzi  ! was wrong dzi - jlm
c      form new wmag at surface
       wmag(iq)=max(sqrt(u(iq,1)**2+v(iq,1)**2),1.)
      enddo    ! iq loop
      do k=2,kl
       if(ndzx.eq.0)dzx=-2.*grav*sig(k)/(dsig(k)*rdry)                  ! Hal's
       if(ndzx.eq.1)dzx=grav*(sig(k-1)+sig(k))/((sig(k-1)-sig(k))*rdry) ! fixup by jlm
       do iq=1,ifull
        dzi=dzx/(t(iq,k-1)+t(iq,k))               ! gwdrag
        dthdz(iq,k)=(thf(iq,k)-thf(iq,k-1))*dzi
       enddo   ! iq loop
      enddo    ! k loop

c**** calculate Brunt-Vaisala frequency
c**** surface bvng() and full levels bvnf(,)
      do iq=1,ifull
c      calculate bvng,  (Brunt-Vaisala-N-ground)
c      limit value of bvng to about 50 degrees per 166 m
!      if unstable (surface to bl) then no gwd (set bvng=0)
!            - happens automatically via next line & effect on bvng & alam
       bvng(iq)=min(.1,sqrt(grav*dthdz(iq,1)/tss(iq)))
       bvnf(iq,kl)=sqrt(max(1.e-20,grav*dthdz(iq,kl)/thf(iq,kl)))    ! jlm fixup
c****  froude number calcs
c****  calculate (fc/f)**2 where fc**2=fc2=0.5  (for Hal & Ch., 1. for jlm)
c**    calc fc2*(t*/n*/wmag)/he**2
       temp(iq)=fc2*tss(iq)/max( bvng(iq)*wmag(iq)*he(iq)**2,1.e-10)  !jlm
      enddo    ! iq loop

c     calculate bvnf at other levels,  (Brunt-Vaisala-N-full)
      do k=1,kl-1
       do iq=1,ifull
        bvnf(iq,k)=sqrt( max(1.e-20,
     .      grav*(dthdz(iq,k)+dthdz(iq,k+1))/(thf(iq,k)+thf(iq,k+1))) ) ! jlm fixup
       enddo   ! iq loop
      enddo    ! k loop

      do k=1,2
       do iq=1,ifull
        uu(iq,k)=
     .        max(0. , u(iq,k)*u(iq,1)+v(iq,k)*v(iq,1))/wmag(iq)
       enddo   ! iq loop
      enddo    ! k loop

c**** set uu() to zero above if uu() zero below
c**** uu>0 at k=1, uu>=0 at k=2 - only set for k=3 to kl  OK
      do k=3,kl
       do iq=1,ifull
        if(uu(iq,k-1).eq.0.)then
          uu(iq,k)=0.
        else
          uu(iq,k)=
     .        max(0. , u(iq,k)*u(iq,1)+v(iq,k)*v(iq,1))/wmag(iq)
        endif  ! (uu(iq,k-1).eq.0.)
       enddo   ! iq loop
      enddo    ! k loop

      do k=1,kl   ! jlm alternative to use actual top level Fr
       do iq=1,ifull
c       calc max(1-fc**2/f**2,0) : put in fni()
c       default Fr for top level is infinity, as used by Canadians & hbg
c       fni(iq,kl)=1.
        fni(iq,k)=max(0. , 1.-
     .     sig(k)*temp(iq)*uu(iq,k)**3/(sigk(k)*bvnf(iq,k)*thf(iq,k)))
       enddo   ! iq loop
      enddo    ! k loop

c     form integral of above*uu**2 from sig=1 to sig=0
      do iq=1,ifull
       fnii(iq)=-fni(iq,1)*dsig(1)*uu(iq,1)**2
      enddo    ! iq loop
      do k=2,kl
       do iq=1,ifull
        fnii(iq)=fnii(iq)-fni(iq,k)*dsig(k)*uu(iq,k)**2
       enddo   ! iq loop
      enddo    ! k loop

!     Chouinard et al. use alpha=.01
      alphaj=0.01*1.e-4  ! jlm   .01 *rhos*g/ps
chal  alphah=0.0075*g/r  ! actually alpha*(g/ps)*rhos  *tss  ! Ch et al
chal  alam(iq)=alphah*he(iq)*bvng(iq)*wmag(iq)*unst(iq)
chal .          /(tss(iq)*max(fnii(iq),1.e-20))       ! cosmetic improvement
      do iq=1,ifull
c      if integral=0., reset to some +ve value
c      form alam=(g/p*).alpha.rhos.he.N*.wmag/integral(above)
cjlm   alam(iq)=alphaj*he(iq)*bvng(iq)*wmag(iq)*unst(iq)
cjlm .                              /max(fnii(iq),1.e-20)
       alam(iq)=alphaj*he(iq)*bvng(iq)*wmag(iq)/max(fnii(iq),1.e-20)
c      define apuw=alam.u1/wmag , apvw=alam.v1/wmag
       apuw(iq)=alam(iq)*u(iq,1)/wmag(iq)
       apvw(iq)=alam(iq)*v(iq,1)/wmag(iq)
      enddo   ! iq loop

c**** form fni=alam*max(--,0) and
c**** solve for uu at t+1 (implicit solution)
      do k=1,kl
       do iq=1,ifull
        uu(iq,k)=2.*uu(iq,k)/
     .             (1.+sqrt(1. + 4.*dt*alam(iq)*fni(iq,k)*uu(iq,k)))
!       N.B. 4.*dt*alam(iq)*fni(iq,k)*uu(iq,k)) can be ~300
       enddo   ! iq loop
      enddo    ! k loop

c**** form dv/dt due to gw-drag at each level
c**** = -alam.v*/wmag.uu(t+1)**2.max(--,0)
      do k=1,kl
       if(npanels.eq.0)then
         do j=2,jl-1                            ! darlam
          do i=2,il-1                           ! darlam
	    iq=i+(j-1)*il
           xxx=uu(iq,k)*uu(iq,k)*fni(iq,k)   ! darlam
           ronout=-apuw(iq)*xxx                ! darlam
           sonout=-apvw(iq)*xxx                ! darlam
           u(iq,k)=u(iq,k)+ronout*dt          ! darlam
           v(iq,k)=v(iq,k)+sonout*dt          ! darlam
          enddo ! i loop                        ! darlam
         enddo  ! j loop                        ! darlam
       else
         if(ngwd.gt.0)then   ! tendencies used
           do iq=1,ifull  ! globpe
            xxx=uu(iq,k)*uu(iq,k)*fni(iq,k)
            ronout=-apuw(iq)*xxx
            sonout=-apvw(iq)*xxx
            un(iq,k)=un(iq,k)+ronout    !  globpex
            vn(iq,k)=vn(iq,k)+sonout    !  globpex
           enddo  ! iq loop
         else   ! non-tendencies for ngwd<0
           do iq=1,ifull  ! globpe
            xxx=uu(iq,k)*uu(iq,k)*fni(iq,k)
            ronout=-apuw(iq)*xxx
            sonout=-apvw(iq)*xxx
            u(iq,k)=u(iq,k)+ronout*dt
            v(iq,k)=v(iq,k)+sonout*dt
           enddo  ! iq loop
         endif  ! (ngwd.gt.0)then   ! tendencies used
       endif    ! (npanels.eq.0) then .. else ..
      enddo     ! k loop
      
      if(ntest==1)then
        write(6,*) 'from gwdrag, ngwd,alam,fnii,apuw,apvw,wmag',
     &    ngwd,alam(idjd),fnii(idjd),apuw(idjd),apvw(idjd),wmag(idjd)
        write(6,*) 'temp,bvng,he,tss',
     &  temp(idjd),bvng(idjd),he(idjd),tss(idjd)
        write(6,*) 't',t(idjd,:)
        write(6,*) 'thf',thf(idjd,:)
        write(6,*) 'bvnf',bvnf(idjd,:)
        write(6,*) 'dthdz',dthdz(idjd,:)
        write(6,*) 'fni',fni(idjd,:)
        write(6,*) 'uu',uu(idjd,:)
        write(6,*) 'un',vn(idjd,:)
        write(6,*) 'vn',vn(idjd,:)
      endif

      return
      end
