      subroutine vadvtvd(tarr,uarr,varr,nvadh_pass,nits,iaero)   ! globpea  version
c                              vadvbott & vadvyu at bottom
!     can show adding tbar has no effect
      use aerosolldr
      use arrays_m
      use cc_mpi
      use cfrac_m
      use diag_m
      use liqwpar_m  ! ifullw
      use map_m
      use nharrs_m
      use sigs_m
      use tkeeps, only : tke,eps
      use tracers_m
      use vvel_m
      use xarrs_m
      implicit none
      include 'newmpar.h'
c     split vertical advection routine; tvd scheme; used with nonlin or upglobal
c     In flux limiter, assuming zero gradient for all top and bottom
c     variables; except extrap at bottom for qg and trace gases  Thu  06-19-1997
      include 'kuocom.h'     ! also with kbsav,ktsav
      include 'parm.h'
      include 'parmdyn.h'
      include 'parmvert.h' ! nthub,nimp,ntvd
      real tarr(ifull,kl),uarr(ifull,kl),varr(ifull,kl)
      real dum(ifull,kl)
      real, dimension(ifull) :: tfact
      integer num,iaero,nqq,npslx
      integer ntr,k
      integer, dimension(ifull) :: nvadh_pass, nits
      data num/0/
      save num
      parameter (npslx=1)  ! 0 off, 1 on for nvad=-4
      parameter (nqq=0)    ! 0 off, 3 possible
!     parameter (nimp=1)  !  0 for original explicit non-flux TVD term
!                            1 for implicit non-flux TVD term
!     parameter (nthub=1) !  1 original
!                            2 higher-order is Lax-Wendroff
!                            0 (not available now) was equivalent to original
!     parameter (ntvd=2)  !  1 van Leer phitvd flux-limiter (original)
!                            2 MC phitvd flux-limiter
!                            3 superbee flux-limiter
!                            0 (not available now) was equivalent to van Leer

      call start_log(vadv_begin)

      tfact=1./real(nvadh_pass)

      if(num==0)then
        num=1
        if(ntvdr==0)then  ! to produce old nvad=4 interps for fluxhi
	   ratha(:)=.5
	   rathb(:)=.5
        endif
        if(myid==0)then
          print *,'In vadvtvd nvad,nvadh_pass,nqq,npslx,ntvdr ',
     .             nvad,nvadh_pass(idjd),nqq,npslx,ntvdr
          print *,'nimp,nthub,ntvd,tfact ',nimp,nthub,ntvd,tfact(idjd)
        endif
      endif

c     t
      call vadvsub(tarr,tfact,nits,0)
      if( (diag.or.nmaxpr==1) .and. mydiag )then
!       These diagnostics don't work with single input/output argument
        write (6,"('tout',9f8.2/4x,9f8.2)") (tarr(idjd,k),k=1,kl)
        write (6,"('t#  ',9f8.2)") diagvals(tarr(:,nlv)) 
!     .           ((tarr(ii+jj*il,nlv),ii=idjd-1,idjd+1),jj=1,-1,-1)
      endif

c     u
      call vadvsub(uarr,tfact,nits,0)
      if( diag .and. mydiag )then
        write (6,"('uout',9f8.2/4x,9f8.2)") (uarr(idjd,k),k=1,kl)
        write (6,"('u#  ',9f8.2)") diagvals(uarr(:,nlv)) 
!     .           ((uarr(ii+jj*il,nlv),ii=idjd-1,idjd+1),jj=-1,1)
      endif

c     v
      call vadvsub(varr,tfact,nits,0)
      if( diag .and. mydiag )then
        write (6,"('vout',9f8.2/4x,9f8.2)") (varr(idjd,k),k=1,kl)
        write (6,"('v#  ',9f8.2)") diagvals(varr(:,nlv)) 
      endif

c     h_nh
      if(nh.ne.0)then
        dum=h_nh(1:ifull,:)
        call vadvsub(dum,tfact,nits,0)
        h_nh(1:ifull,:)=dum
      endif     ! (nh.ne.0)

c     pslx
      if(npslx==1.and.nvad<=-4)then  ! handles -9 too
        dum=pslx(1:ifull,:)
        call vadvsub(dum,tfact,nits,0)
        pslx(1:ifull,:)=dum
      endif  ! (npslx==1.and.nvad==-4)

      if(mspec==1.and.abs(nvad).ne.9)then   ! advect qg and gases after preliminary step

c      qg
       dum=qg(1:ifull,:)
       call vadvsub(dum,tfact,nits,nqq)
       qg(1:ifull,:)=dum
       if( diag .and. mydiag )then
        write (6,"('qout',9f8.2/4x,9f8.2)") (1000.*qg(idjd,k),k=1,kl)
        write (6,"('qg# ',3p9f8.2)") diagvals(qg(:,nlv)) 
       endif

       if(ldr/=0)then
        dum=qlg(1:ifull,:)
        call vadvsub(dum,tfact,nits,0)
        qlg(1:ifull,:)=dum
        dum=qfg(1:ifull,:)
        call vadvsub(dum,tfact,nits,0)
        qfg(1:ifull,:)=dum
        dum=qrg(1:ifull,:)
        call vadvsub(dum,tfact,nits,0)
        qrg(1:ifull,:)=dum
        !dum=cfrac(1:ifull,:)
        !call vadvsub(dum,tfact,nits,0)
        !cfrac(1:ifull,:)=dum
        dum=cffall(1:ifull,:)
        call vadvsub(dum,tfact,nits,0)
        cffall(1:ifull,:)=dum
        if( diag .and. mydiag )then
         write (6,"('lout',9f8.2/4x,9f8.2)") (1000.*qlg(idjd,k),k=1,kl)
         write (6,"('qlg#',3p9f8.2)") diagvals(qlg(:,nlv)) 
         write (6,"('fout',9f8.2/4x,9f8.2)") (1000.*qfg(idjd,k),k=1,kl)
         write (6,"('qfg#',3p9f8.2)") diagvals(qfg(:,nlv)) 
        endif
       endif      ! if(ldr.ne.0)

       if(ngas>0.or.nextout>=4)then
        do ntr=1,ntrac
         dum=tr(1:ifull,:,ntr)
         call vadvsub(dum,tfact,nits,0)
         tr(1:ifull,:,ntr)=dum
        enddo      ! ntr loop
       endif      ! (nextout>=4)

       !--------------------------------------------------------------
       ! MJT tke
       if(nvmix==6)then
        dum=eps(1:ifull,:)
        call vadvsub(dum,tfact,nits,0)
        eps(1:ifull,:)=dum
        dum=tke(1:ifull,:)
        call vadvsub(dum,tfact,nits,0)
        tke(1:ifull,:)=dum
       endif      ! if(nvmix.eq.6)
       !--------------------------------------------------------------

       !--------------------------------------------------------------
       ! MJT aerosols
       if (abs(iaero)==2) then
        do ntr=1,naero
         dum=xtg(1:ifull,:,ntr)
         call vadvsub(dum,tfact,nits,0)  
         xtg(1:ifull,:,ntr)=dum
        end do
       end if
       !--------------------------------------------------------------

      endif       ! if(mspec==1.and.abs(nvad).ne.9)

      call end_log(vadv_end)
 
      return
      end
      
      ! Subroutine to perform generic TVD advection
      subroutine vadvsub(tarr,tfact,nits,nqq)
      
      use sigs_m
      use vvel_m
      
      implicit none
      
      include 'newmpar.h'
      include 'parmvert.h'
      
      integer, intent(in) :: nqq
      integer, dimension(ifull) :: nits
      integer k,iq,kp,kx,i,its_g
      real, dimension(ifull), intent(in) :: tfact
      real rat,phitvd,fluxhi,fluxlo
      real hdsdot
      real, dimension(ifull,kl), intent(inout) :: tarr
      real, dimension(ifull,0:kl) :: delt,fluxh
      real, dimension(ifull,kl) :: xin
      logical, dimension(ifull) :: msk

      its_g=maxval(nits)

      do iq=1,ifull
c      fluxh(k) is located at level k+.5
       fluxh(iq,0)=0.
       fluxh(iq,kl)=0.
       delt(iq,0)=0.      ! for T,u,v
       delt(iq,kl)=0.     ! for T,u,v
      enddo    ! iq loop
      
      if(ntvdr==2) xin=tarr(1:ifull,:)
      if(nqq==3)then
        tarr(1:ifull,:)=tarr(1:ifull,:)**(1./3.)
      endif 
      
      do k=1,kl-1
       do iq=1,ifull
         delt(iq,k)=tarr(iq,k+1)-tarr(iq,k)
       enddo    ! iq loop
      enddo     ! k loop
      do iq=1,ifull
       delt(iq,0)=min(delt(iq,1),tarr(iq,1))       ! for non-negative tt
      enddo    ! iq loop
      do k=1,kl-1  ! for fluxh at interior (k + 1/2)  half-levels
       do iq=1,ifull
        kp=sign(1.,sdot(iq,k+1))
        kx=k+(1-kp)/2  !  k for sdot +ve,  k+1 for sdot -ve
        rat=delt(iq,k-kp)/(delt(iq,k)+sign(1.e-20,delt(iq,k)))
        if(ntvd==1)then
         phitvd=(rat+abs(rat))/(1.+abs(rat))       ! 0 for -ve rat
        else if(ntvd==2) then
         phitvd=max(0.,min(2.*rat,.5+.5*rat,2.))   ! 0 for -ve rat
        else if(ntvd==3) then
         phitvd=max(0.,min(1.,2.*rat),min(2.,rat)) ! 0 for -ve rat
        end if
        if(nthub==1) then
         fluxhi=rathb(k)*tarr(iq,k)+ratha(k)*tarr(iq,k+1)
        else if(nthub==2)then     ! higher order scheme
         fluxhi=rathb(k)*tarr(iq,k)+ratha(k)*tarr(iq,k+1)
     &                 -.5*delt(iq,k)*tfact(iq)*sdot(iq,k+1)
        endif  ! (nthub==2)
        fluxlo=tarr(iq,kx)
        fluxh(iq,k)=sdot(iq,k+1)*(fluxlo+phitvd*(fluxhi-fluxlo))
       enddo    ! iq loop
      enddo     ! k loop
      if(nimp==1)then
       do k=1,kl
        do iq=1,ifull
         hdsdot=.5*tfact(iq)*(sdot(iq,k+1)-sdot(iq,k))
         tarr(iq,k)=(tarr(iq,k)
     &               +tfact(iq)*(fluxh(iq,k-1)-fluxh(iq,k))
     &               +hdsdot*tarr(iq,k) )/(1.-hdsdot)
        end do
       end do
      else
       do k=1,kl
        do iq=1,ifull       
         tarr(iq,k)=tarr(iq,k)
     &               +tfact(iq)*(fluxh(iq,k-1)-fluxh(iq,k)
     &               +tarr(iq,k)*(sdot(iq,k+1)-sdot(iq,k)))
        enddo     ! iq loop
       enddo      ! k loop
      endif   ! (nimp==1)
      
      do i=2,its_g
       msk=(nits>=i)

       do iq=1,ifull
        if (msk(iq)) then
         do k=1,kl-1
          delt(iq,k)=tarr(iq,k+1)-tarr(iq,k)
         enddo     ! k loop
         delt(iq,0)=min(delt(iq,1),tarr(iq,1))       ! for non-negative tt
         do k=1,kl-1  ! for fluxh at interior (k + 1/2)  half-levels
          kp=sign(1.,sdot(iq,k+1))
          kx=k+(1-kp)/2  !  k for sdot +ve,  k+1 for sdot -ve
          rat=delt(iq,k-kp)/(delt(iq,k)+sign(1.e-20,delt(iq,k)))
          if(ntvd==1)then
           phitvd=(rat+abs(rat))/(1.+abs(rat))       ! 0 for -ve rat
          else if(ntvd==2) then
           phitvd=max(0.,min(2.*rat,.5+.5*rat,2.))   ! 0 for -ve rat
          else if(ntvd==3) then
           phitvd=max(0.,min(1.,2.*rat),min(2.,rat)) ! 0 for -ve rat
          end if
          if(nthub==1) then
           fluxhi=rathb(k)*tarr(iq,k)+ratha(k)*tarr(iq,k+1)
          else if(nthub==2)then     ! higher order scheme
           fluxhi=rathb(k)*tarr(iq,k)+ratha(k)*tarr(iq,k+1)
     &                   -.5*delt(iq,k)*tfact(iq)*sdot(iq,k+1)
          endif  ! (nthub==2)
          fluxlo=tarr(iq,kx)
          fluxh(iq,k)=sdot(iq,k+1)*(fluxlo+phitvd*(fluxhi-fluxlo))
         enddo     ! k loop
         if(nimp==1)then
          do k=1,kl
           hdsdot=.5*tfact(iq)*(sdot(iq,k+1)-sdot(iq,k))
           tarr(iq,k)=(tarr(iq,k)
     &                 +tfact(iq)*(fluxh(iq,k-1)-fluxh(iq,k))
     &                 +hdsdot*tarr(iq,k) )/(1.-hdsdot)
          end do
         else
          do k=1,kl
           tarr(iq,k)=tarr(iq,k)
     &                 +tfact(iq)*(fluxh(iq,k-1)-fluxh(iq,k)
     &                 +tarr(iq,k)*(sdot(iq,k+1)-sdot(iq,k)))
          enddo      ! k loop
         endif   ! (nimp==1)
        end if ! msk
       end do ! iq
      
      end do

      if(nqq==3)then
        tarr(1:ifull,:)=tarr(1:ifull,:)**3
      endif       ! (nqq==3)
!!    tarr(1:ifull,1)=max(xin(:,1),min(tarr(1:ifull,1),xin(:,2)))
!!    tarr(1:ifull,kl)=max(xin(:,kl),min(tarr(1:ifull,kl),xin(:,kl-1)))
      if(ntvdr==2)then  ! different top & bottom
        do iq=1,ifull
         if(sdot(iq,2)>0.)tarr(iq,1)=xin(iq,1)
         if(sdot(iq,kl)<0.)tarr(iq,kl)=xin(iq,kl)
        enddo
      endif      ! (ntvdr==2)
      
      return
      end subroutine vadvsub
