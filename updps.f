      subroutine updps(iadj)    
      use cc_mpi
c     use diag_m             ! for calls to maxmin
      include 'newmpar.h'
      include 'arrays.h'
      include 'const_phys.h'
      include 'indices.h'
      include 'map.h'
      include 'nlin.h'  ! savs
      include 'savuvt.h'
      include 'parm.h'
      include 'parmdyn.h'
      include 'sigs.h'
      include 'vvel.h'    ! sdot
      include 'xarrs.h'
      real derpsl(ifull),cc(ifull+iextra,kl),dd(ifull+iextra,kl)
      real c0(ifull,kl),d0(ifull,kl)
      common/work3d/d(ifull,kl) ! possibly shared adjust5 & updps
      real savs1(ifull,2:kl),savu1(ifull,kl),savv1(ifull,kl)
      real sbar(ifull,2:kl)
      common/savuv1/savs1,savu1,savv1,sbar 
      real omgf(ifull,kl)
      equivalence (omgf,dpsldt)
      real sdotin(ifull,kl),pslxin(ifull,kl),omgfin(ifull,kl)
      save num
      data num/0/
!     called for mup.ne.1      mup=3 for simple centred
      
      if(mup<=-4)then
        if(ktau<3)then
          return
        else
          sdotin(:,1:kl)=sdot(:,1:kl)
          omgfin(:,1:kl)=omgf(:,1:kl)
          pslxin(:,1:kl)=pslx(1:ifull,1:kl)
	 endif
      endif
      
      if(mup>=4)then
c       call staguv(u,v,cc,dd)
	 c0(:,:)=u(1:ifull,:)  ! can avoid this fiddle if staguv has +iextra
	 d0(:,:)=v(1:ifull,:)
        call staguv(c0,d0,c0,d0)
	 cc(1:ifull,:)=c0(:,:)
	 dd(1:ifull,:)=d0(:,:)
	 do k=1,kl
	  cc(1:ifull,k)=cc(1:ifull,k)/emu(1:ifull)
	  dd(1:ifull,k)=dd(1:ifull,k)/emv(1:ifull)
	 enddo
        call boundsuv(cc,dd)
        do k=1,kl
!cdir nodep
         do iq=1,ifull
!         N.B. this div is calculated on the staggered grid
          d(iq,k)=(cc(iq,k)-cc(iwu(iq),k) +dd(iq,k)-dd(isv(iq),k))  
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
      endif       ! (mup>=4) .. else ..
	
      if(mup==-1.or.mup<=-4.or.(mup==-3.and.iadj==0))then
        tx(:,1)=zs(:)/(rdry*300.)
        tx(:,2)=psl(:)+tx(:,1)
        call bounds(tx)
        do k=1,kl
*cdir nodep
         do iq=1,ifull
          pslx(iq,k)=-em(iq)*(
     .                u(iq,k)*(tx(ie(iq),2)-tx(iw(iq),2))
     .               +v(iq,k)*(tx(in(iq),2)-tx(is(iq),2)))/(2.*ds)
          if(u(iq,k)>0.)then
            pslx(iq,k)=pslx(iq,k)+em(iq)*(tx(ie(iq),1)-tx(iq,1))
     .                                  *u(iq,k)/ds
          else
            pslx(iq,k)=pslx(iq,k)+em(iq)*(tx(iq,1)-tx(iw(iq),1))   
     .                                  *u(iq,k)/ds
          endif 
          if(v(iq,k)>0.)then
            pslx(iq,k)=pslx(iq,k)+em(iq)*(tx(in(iq),1)-tx(iq,1))   
     .                                  *v(iq,k)/ds
          else
            pslx(iq,k)=pslx(iq,k)+em(iq)*(tx(iq,1)-tx(is(iq),1))   
     .                                  *v(iq,k)/ds
          endif 
          enddo   ! iq loop
         enddo    ! k  loop
      elseif(mup==-2.or.(mup==-3.and.iadj==1))then
        tx(:,1)=zs(:)/(rdry*300.)
        tx(:,2)=psl(:)+tx(:,1)
        call bounds(tx)
        do k=1,kl
*cdir nodep
         do iq=1,ifull
          pslx(iq,k)=-em(iq)*(
     .                u(iq,k)*(tx(ie(iq),2)-tx(iw(iq),2))
     .               +v(iq,k)*(tx(in(iq),2)-tx(is(iq),2)))/(2.*ds)
          if(u(iq,k)<0.)then
            pslx(iq,k)=pslx(iq,k)+em(iq)*(tx(ie(iq),1)-tx(iq,1))
     .                                  *u(iq,k)/ds
          else
            pslx(iq,k)=pslx(iq,k)+em(iq)*(tx(iq,1)-tx(iw(iq),1))   
     .                                  *u(iq,k)/ds
          endif 
          if(v(iq,k)<0.)then
            pslx(iq,k)=pslx(iq,k)+em(iq)*(tx(in(iq),1)-tx(iq,1))   
     .                                  *v(iq,k)/ds
          else
            pslx(iq,k)=pslx(iq,k)+em(iq)*(tx(iq,1)-tx(is(iq),1))   
     .                                  *v(iq,k)/ds
          endif 
          enddo   ! iq loop
         enddo    ! k  loop
      elseif((mup.ne.2.and.mup.ne.5).or.num==0)then
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
       endif       ! (mup.ne.2.or.num==0)

!      integrate vertically {0 to 1} to get d(ln(psl))/dt
       derpsl(:)=0.
       do k=1,kl
        derpsl(:)=derpsl(:)-dsig(k)*(pslx(1:ifull,k)-d(:,k))
       enddo      ! k  loop
!      put -D(ln(psl))/Dt into pslx, i.e. to equal D+d(sdot)/d(sig)
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

      if(mup<=-4)then
        do iq=1,ifull
	  ind=0
	  do k=2,kl
	   if((sign(1.,sdot(iq,k)).ne.sign(1.,savs(iq,k))).and.
     .       (sign(1.,savs1(iq,k)).ne.sign(1.,savs(iq,k))))ind=1
	  enddo
	  if(ind==0)then
	    do k=1,kl
	     sdot(iq,k)=sdotin(iq,k)
	     omgf(iq,k)=omgfin(iq,k)
	     pslx(iq,k)=pslxin(iq,k)
	    enddo
	  else
	    do k=1,kl
	     t(iq,k)=t(iq,k)+dt*tn(iq,k)
	     tn(iq,k)=0.
	     if(mup.eq.-5)then
	       u(iq,k)=u(iq,k)+dt*un(iq,k)
	       un(iq,k)=0.
	       v(iq,k)=v(iq,k)+dt*vn(iq,k)
	       vn(iq,k)=0.
	     endif
	    enddo
	  endif  ! (ind==0) .. else ..
	 enddo  ! iq loop
      endif

      return
      end
