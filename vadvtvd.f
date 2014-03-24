      module vadv
      
      private
      public vadvtvd
      
      contains

      subroutine vadvtvd(tarr,uarr,varr,nvadh_pass,nits,iaero)   ! globpea  version
!                              vadvbott & vadvyu at bottom
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
!     split vertical advection routine; tvd scheme; used with nonlin or upglobal
!     In flux limiter, assuming zero gradient for all top and bottom
!     variables; except extrap at bottom for qg and trace gases  Thu  06-19-1997
      include 'kuocom.h'     ! also with kbsav,ktsav
      include 'parm.h'
      include 'parmdyn.h'
      include 'parmvert.h' ! nthub,nimp,ntvd
      real, dimension(:,:), intent(inout) :: tarr,uarr,varr
      real, dimension(ifull,kl) :: dum
      real, dimension(ifull) :: tfact
      integer num,iaero,npslx
      integer ntr,k
      integer, dimension(ifull) :: nvadh_pass, nits
      integer, dimension(ifull,kl-1) :: kp,kx
      data num/0/
      save num
      parameter (npslx=1)  ! 0 off, 1 on for nvad=-4
!     parameter (nimp=1)  !  0 for original explicit non-flux TVD term
!                            1 for implicit non-flux TVD term
!     parameter (nthub=1) !  1 original
!                            2 higher-order is Lax-Wendroff
!                            0 (not available now) was equivalent to original
!     parameter (ntvd=2)  !  1 van Leer phitvd flux-limiter (original)
!                            2 MC phitvd flux-limiter
!                            3 superbee flux-limiter
!                            0 (not available now) was equivalent to van Leer

#include "log.h"
 
      START_LOG(vadv)

      tfact=1./real(nvadh_pass)

      if(num==0)then
        num=1
        if(ntvdr==0)then  ! to produce old nvad=4 interps for fluxhi
         ratha(:)=.5
         rathb(:)=.5
        endif
#ifdef debug
        if(mydiag)then
          write(6,*) 'In vadvtvd nvad,nvadh_pass,npslx,ntvdr ',
     &             nvad,nvadh_pass(idjd),npslx,ntvdr
          write(6,*) 'nimp,nthub,ntvd,tfact ',nimp,nthub,ntvd,
     &             tfact(idjd)
        endif
#endif
      endif

      do k=1,kl-1
        kp(:,k)=sign(1.,sdot(:,k+1))
        kx(:,k)=k+(1-kp(:,k))/2 !  k for sdot +ve,  k+1 for sdot -ve
      end do

!     t
      call vadvsub(tarr,tfact,nits,kp,kx)
#ifdef debug
      if( (diag.or.nmaxpr==1) .and. mydiag )then
!       These diagnostics don't work with single input/output argument
        write (6,"('tout',9f8.2/4x,9f8.2)") (tarr(idjd,k),k=1,kl)
        write (6,"('t#  ',9f8.2)") diagvals(tarr(1:ifull,nlv)) 
!     .           ((tarr(ii+jj*il,nlv),ii=idjd-1,idjd+1),jj=1,-1,-1)
      endif
#endif

!     u
      call vadvsub(uarr,tfact,nits,kp,kx)
#ifdef debug
      if( diag .and. mydiag )then
        write (6,"('uout',9f8.2/4x,9f8.2)") (uarr(idjd,k),k=1,kl)
        write (6,"('u#  ',9f8.2)") diagvals(uarr(:,nlv)) 
!     .           ((uarr(ii+jj*il,nlv),ii=idjd-1,idjd+1),jj=-1,1)
      endif
#endif

!     v
      call vadvsub(varr,tfact,nits,kp,kx)
#ifdef debug
      if( diag .and. mydiag )then
        write (6,"('vout',9f8.2/4x,9f8.2)") (varr(idjd,k),k=1,kl)
        write (6,"('v#  ',9f8.2)") diagvals(varr(:,nlv)) 
      endif
#endif

!     h_nh
      if(nh/=0.and.npslx==1.and.nvad<=-4)then
        call vadvsub(h_nh,tfact,nits,kp,kx)
      endif     ! (nh.ne.0)

!     pslx
      if(npslx==1.and.nvad<=-4)then  ! handles -9 too
        call vadvsub(pslx,tfact,nits,kp,kx)
      endif  ! (npslx==1.and.nvad==-4)

      if(mspec==1.and.abs(nvad)/=9)then   ! advect qg and gases after preliminary step

!      qg
       call vadvsub(qg,tfact,nits,kp,kx)
#ifdef debug
       if( diag .and. mydiag )then
        write (6,"('qout',9f8.2/4x,9f8.2)") (1000.*qg(idjd,k),k=1,kl)
        write (6,"('qg# ',3p9f8.2)") diagvals(qg(:,nlv)) 
       endif
#endif

       if(ldr/=0)then
        call vadvsub(qlg,tfact,nits,kp,kx)
        call vadvsub(qfg,tfact,nits,kp,kx)
        call vadvsub(qrg,tfact,nits,kp,kx)
        !call vadvsub(cfrac,tfact,nits,kp,kx)
        call vadvsub(cffall,tfact,nits,kp,kx)
#ifdef debug
        if( diag .and. mydiag )then
         write (6,"('lout',9f8.2/4x,9f8.2)") (1000.*qlg(idjd,k),k=1,kl)
         write (6,"('qlg#',3p9f8.2)") diagvals(qlg(:,nlv)) 
         write (6,"('fout',9f8.2/4x,9f8.2)") (1000.*qfg(idjd,k),k=1,kl)
         write (6,"('qfg#',3p9f8.2)") diagvals(qfg(:,nlv)) 
        endif
#endif
       endif      ! if(ldr.ne.0)

       if(ngas>0.or.nextout>=4)then
        do ntr=1,ntrac
         call vadvsub(tr(:,:,ntr),tfact,nits,kp,kx)
        enddo      ! ntr loop
       endif      ! (nextout>=4)

       if(nvmix==6)then
        call vadvsub(eps,tfact,nits,kp,kx)
        call vadvsub(tke,tfact,nits,kp,kx)
       endif      ! if(nvmix.eq.6)

       if (abs(iaero)==2) then
        do ntr=1,naero
         call vadvsub(xtg(:,:,ntr),tfact,nits,kp,kx)
        end do
       end if

      endif       ! if(mspec==1.and.abs(nvad).ne.9)

      END_LOG(vadv)
 
      return
      end subroutine vadvtvd
      
      ! Subroutine to perform generic TVD advection
      subroutine vadvsub(tarr,tfact,nits,kp,kx)
      
      use sigs_m
      use vvel_m
      
      implicit none
      
      include 'newmpar.h'
      include 'parmvert.h'
      
      integer, dimension(ifull), intent(in) :: nits
      integer, dimension(ifull,kl-1), intent(in) :: kp,kx
      integer k,iq,i
      real, dimension(ifull), intent(in) :: tfact
      real, dimension(ifull) :: rat,phitvd,fluxhi,fluxlo,hdsdot
      real, dimension(:,:), intent(inout) :: tarr
      real, dimension(ifull,0:kl) :: delt,fluxh
      real, dimension(ifull,kl) :: xin

!     fluxh(k) is located at level k+.5
      fluxh(:,0)=0.
      fluxh(:,kl)=0.
      
      if(ntvdr==2) then
        xin=tarr(1:ifull,:)
      end if

      delt(:,kl)=0.     ! for T,u,v
      delt(:,1:kl-1)=tarr(1:ifull,2:kl)-tarr(1:ifull,1:kl-1)
      delt(:,0)=min(delt(:,1),tarr(1:ifull,1))       ! for non-negative tt
      do k=1,kl-1  ! for fluxh at interior (k + 1/2)  half-levels
       do iq=1,ifull
        rat(iq)=delt(iq,k-kp(iq,k))/(delt(iq,k)+sign(1.e-20,delt(iq,k)))
       end do
       if(ntvd==1)then
         phitvd=(rat+abs(rat))/(1.+abs(rat))       ! 0 for -ve rat
       else if(ntvd==2) then
         phitvd=max(0.,min(2.*rat,.5+.5*rat,2.))   ! 0 for -ve rat
       else if(ntvd==3) then
         phitvd=max(0.,min(1.,2.*rat),min(2.,rat)) ! 0 for -ve rat
       end if
       if(nthub==1) then
         fluxhi=rathb(k)*tarr(1:ifull,k)+ratha(k)*tarr(1:ifull,k+1)
       else if(nthub==2)then     ! higher order scheme
         fluxhi=rathb(k)*tarr(1:ifull,k)+ratha(k)*tarr(1:ifull,k+1)
     &                 -.5*delt(:,k)*tfact(:)*sdot(:,k+1)
       endif  ! (nthub==2)
       do iq=1,ifull
        fluxlo(iq)=tarr(iq,kx(iq,k))
       end do
       fluxh(:,k)=sdot(:,k+1)*(fluxlo+phitvd*(fluxhi-fluxlo))
      enddo     ! k loop
      if(nimp==1)then
       do k=1,kl
         hdsdot=.5*tfact(:)*(sdot(:,k+1)-sdot(:,k))
         tarr(1:ifull,k)=(tarr(1:ifull,k)
     &               +tfact(:)*(fluxh(:,k-1)-fluxh(:,k))
     &               +hdsdot*tarr(1:ifull,k) )/(1.-hdsdot)
       end do
      else
       do k=1,kl
         tarr(1:ifull,k)=tarr(1:ifull,k)
     &               +tfact(:)*(fluxh(:,k-1)-fluxh(:,k)
     &               +tarr(1:ifull,k)*(sdot(:,k+1)-sdot(:,k)))
       enddo      ! k loop
      endif   ! (nimp==1)

      do iq=1,ifull      
       do i=2,nits(iq)

         do k=1,kl-1
          delt(iq,k)=tarr(iq,k+1)-tarr(iq,k)
         enddo     ! k loop
         delt(iq,0)=min(delt(iq,1),tarr(iq,1))       ! for non-negative tt
         do k=1,kl-1  ! for fluxh at interior (k + 1/2)  half-levels
          rat(iq)=delt(iq,k-kp(iq,k))
     &      /(delt(iq,k)+sign(1.e-20,delt(iq,k)))
          if(ntvd==1)then
           phitvd(iq)=(rat(iq)+abs(rat(iq)))/(1.+abs(rat(iq)))       ! 0 for -ve rat
          else if(ntvd==2) then
           phitvd(iq)=max(0.,min(2.*rat(iq),.5+.5*rat(iq),2.))       ! 0 for -ve rat
          else if(ntvd==3) then
           phitvd(iq)=max(0.,min(1.,2.*rat(iq)),min(2.,rat(iq)))     ! 0 for -ve rat
          end if
          if(nthub==1) then
           fluxhi(iq)=rathb(k)*tarr(iq,k)+ratha(k)*tarr(iq,k+1)
          else if(nthub==2)then     ! higher order scheme
           fluxhi(iq)=rathb(k)*tarr(iq,k)+ratha(k)*tarr(iq,k+1)
     &                   -.5*delt(iq,k)*tfact(iq)*sdot(iq,k+1)
          endif  ! (nthub==2)
          fluxlo(iq)=tarr(iq,kx(iq,k))
          fluxh(iq,k)=sdot(iq,k+1)*(fluxlo(iq)
     &      +phitvd(iq)*(fluxhi(iq)-fluxlo(iq)))
         enddo     ! k loop
         if(nimp==1)then
          do k=1,kl
           hdsdot(iq)=.5*tfact(iq)*(sdot(iq,k+1)-sdot(iq,k))
           tarr(iq,k)=(tarr(iq,k)
     &                 +tfact(iq)*(fluxh(iq,k-1)-fluxh(iq,k))
     &                 +hdsdot(iq)*tarr(iq,k) )/(1.-hdsdot(iq))
          end do
         else
          do k=1,kl
           tarr(iq,k)=tarr(iq,k)
     &                 +tfact(iq)*(fluxh(iq,k-1)-fluxh(iq,k)
     &                 +tarr(iq,k)*(sdot(iq,k+1)-sdot(iq,k)))
          enddo      ! k loop
         endif   ! (nimp==1)

       end do ! i
      end do ! iq

!!    tarr(1:ifull,1)=max(xin(:,1),min(tarr(1:ifull,1),xin(:,2)))
!!    tarr(1:ifull,kl)=max(xin(:,kl),min(tarr(1:ifull,kl),xin(:,kl-1)))
      if(ntvdr==2)then  ! different top & bottom
        where (sdot(:,2)>0.)
         tarr(1:ifull,1)=xin(:,1)
        end where
        where (sdot(:,kl)<0.)
         tarr(1:ifull,kl)=xin(:,kl)
        end where
      endif      ! (ntvdr==2)
      
      return
      end subroutine vadvsub
      
      end module vadv
