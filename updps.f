      subroutine updps    
      use cc_mpi
      include 'newmpar.h'
      include 'arrays.h'
      include 'indices.h'
      include 'map.h'
      include 'parm.h'
      include 'parmdyn.h'
      include 'sigs.h'
      include 'vvel.h'    ! sdot
      include 'xarrs.h'
      real derpsl(ifull),cc(ifull+iextra,kl),dd(ifull+iextra,kl)
      common/work3/d(ifull,kl),dum3(ifull,4*kl)
      real omgf(ifull,kl)
      equivalence (omgf,dpsldt)
      save num
      data num/0/
      if(mup.ge.4)then
        call staguv(u,v,cc,dd)
        call bounds(psl)
        call boundsuv(u,v)
        do k=1,kl
!cdir nodep
         do iq=1,ifull
!         N.B. this div is calculated on the staggered grid
          d(iq,k)=(cc(iq,k)/emu(iq)-cc(iwu(iq),k)/emu(iwu(iq))  
     .            +dd(iq,k)/emv(iq)-dd(isv(iq),k)/emv(isv(iq)))  
     .            *em(iq)**2/ds
         enddo   ! iq loop
        enddo    ! k  loop
      else
        call bounds(psl)
        call boundsuv(u,v)
        do k=1,kl
*cdir nodep
         do iq=1,ifull
!          N.B. this div is calculated on the non-staggered grid
!          but in adjust5 it is on the staggered grid
           d(iq,k)=
     .        (u(ieu(iq),k)/em(ie(iq))-u(iwu(iq),k)/em(iw(iq))     
     .        +v(inv(iq),k)/em(in(iq))-v(isv(iq),k)/em(is(iq)))
     .        *em(iq)**2/(2.*ds)  
          enddo   ! iq loop
         enddo    ! k  loop
      endif       ! (mup.eq.4) .. else ..
	
      if((mup.ne.2.and.mup.ne.5).or.num.eq.0)then
        num=1
!       put -D(ln(psl))/Dt +d(ln(psl))/dt into pslx (i.e. Nps)
        do k=1,kl
*cdir nodep
         do iq=1,ifull
          pslx(iq,k)=-em(iq)*(
     .                u(iq,k)*(psl(ie(iq))-psl(iw(iq)))
     .               +v(iq,k)*(psl(in(iq))-psl(is(iq))))/(2.*ds)
          enddo   ! iq loop
         enddo    ! k  loop
       endif       ! (mup.ne.2.or.num.eq.0)

!       integrate vertically {0 to 1} to get d(ln(psl))/dt
        derpsl(:)=0.
        do k=1,kl
         derpsl(:)=derpsl(:)-dsig(k)*(pslx(1:ifull,k)-d(:,k))
        enddo      ! k  loop
!       put -D(ln(psl))/Dt into pslx, i.e. to equal D+d(sdot)/d(sig)
        do k=1,kl
         pslx(1:ifull,k)=-derpsl(:)+pslx(1:ifull,k)
        enddo      ! k  loop

!     calculate sdot (at level k-.5) by vert. integ. {0 to sig(k-.5)} 
      sdot(:,1)=0.
      sdot(:,kl+1)=0.
      do k=kl,2,-1
       sdot(:,k)=sdot(:,k+1)-dsig(k)*(pslx(1:ifull,k)-d(:,k))
      enddo      ! k  loop

!     full-level omega/ps into omgf (equivalenced to dpsldt)
      do k=1,kl
       omgf(:,k)=rata(k)*sdot(:,k+1)+ratb(k)*sdot(:,k)
     .           -sig(k)*pslx(1:ifull,k)
      enddo      ! k  loop

!     convert sdot (at level k-.5) into units of grid-steps/timestep
!     half-level sigma-dot is used for vertical advection
      do k=2,kl
       sdot(:,k)=sdot(:,k)*dt/(sig(k)-sig(k-1))
      enddo    ! k  loop

      return
      end
