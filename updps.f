      subroutine updps    ! globpea version
!     this one includes stuff moved from beginning of nonlin
!     written so ee & omgfnl lines can be removed if not needed
c     simplest nps=2 only
      include 'newmpar.h'
      include 'arrays.h'
      include 'indices.h'
      include 'map.h'
      include 'parm.h'
      include 'sigs.h'
      include 'vvel.h'    ! sdot
      include 'xarrs.h'
      common/work3/d(ifull,kl),e(ifull,kl),ee(ifull,kl),dum2(ifull,2*kl)
      real omgf(ifull,kl)
      equivalence (omgf,dpsldt)
      do k=1,kl
*cdir nodep
       do iq=1,ifull
         dpsldt(iq,k)=-em(iq)*(
     .              u(iq,k)*(psl(ie(iq))-psl(iw(iq)))
     .             +v(iq,k)*(psl(in(iq))-psl(is(iq))))/(2.*ds)
c        d(iq,k)=(u(iq,k)/emu(iq)-u(iwu(iq),k)/emu(iwu2(iq))
c    .            +v(iq,k)/emv(iq)-v(isv(iq),k)/emv(isv2(iq)))
c    .            *em(iq)**2/ds  -dpsldt(iq,k)   ! N.B. -dpsldt term
         d(iq,k)=
     .     (u(ieu(iq),k)/em(ie(iq))-u(iwu(iq),k)/em(iw(iq))     
     .      +v(inv(iq),k)/em(in(iq))-v(isv(iq),k)/em(is(iq)))
     .      *em(iq)**2/(2.*ds)  -dpsldt(iq,k)   ! N.B. -dpsldt term
        enddo   ! iq loop
       enddo    ! k  loop

!     vert. integ. {0 to sig(k-.5)} div into e
      do iq=1,ifull
       ee(iq,kl)=dsig(kl)*dpsldt(iq,kl)
       e(iq,kl)=-dsig(kl)*d(iq,kl)
      enddo   ! iq loop
      do k=kl-1,1,-1
       do iq=1,ifull
        ee(iq,k)=ee(iq,k+1)+dsig(k)*dpsldt(iq,k)  ! just for omgfnl
        e(iq,k)=e(iq,k+1)-dsig(k)*d(iq,k)
       enddo     ! iq loop
      enddo      ! k  loop

      do k=1,kl
       do iq=1,ifull
        pslx(iq,k)=dpsldt(iq,k)+e(iq,1)
       enddo     ! iq loop
      enddo      ! k  loop
 
!     full-level omega/ps into omgf (equivalenced to dpsldt)
      do k=1,kl-1
       do iq=1,ifull
 !      omgfnl(iq,k)=-rata(k)*ee(iq,k+1)-ratb(k)*ee(iq,k)
 !   .                 -sig(k)*dpsldt(iq,k)
        omgf(iq,k)=-rata(k)*e(iq,k+1)-ratb(k)*e(iq,k)
     .               -sig(k)*dpsldt(iq,k)
       enddo     ! iq loop
      enddo      ! k  loop
      do iq=1,ifull
!      omgfnl(iq,kl)=-ratb(kl)*ee(iq,kl)-sig(kl)*dpsldt(iq,kl)
       omgf(iq,kl)=-ratb(kl)*e(iq,kl)-sig(kl)*dpsldt(iq,kl)
       sdot(iq,1)=0.
       sdot(iq,kl+1)=0.
      enddo      ! iq loop

!     calculate sdot (at level k-.5) in units of grid-steps/timestep
!     half-level sigma-dot is used for vertical advection
        do k=2,kl
         do iq=1,ifull
!         sdot(iq,k)=sigmh(k)*e(iq)-e(iq,k)
!         sdot(iq,k)=sdot(iq,k)*dt/(sig(k)-sig(k-1))
          sdot(iq,k)=(sigmh(k)*e(iq,1)-e(iq,k))
     .                 *dt/(sig(k)-sig(k-1))
         enddo   ! iq loop
        enddo    ! k  loop

      return
      end
