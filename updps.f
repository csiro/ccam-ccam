      subroutine updps(iadj)    
      use cc_mpi
      use diag_m             ! for calls to maxmin
      use indices_m
      use map_m
      use vecsuv_m
      use xyzinfo_m
      include 'newmpar.h'
      include 'arrays.h'
      include 'const_phys.h'
      include 'nlin.h'  ! savs
      include 'savuvt.h'
      include 'parm.h'
      include 'parmdyn.h'
      include 'parmhor.h'
      include 'sigs.h'
      include 'vvel.h'     ! sdot
      include 'xarrs.h'
      real derpsl(ifull),cc(ifull+iextra,kl),dd(ifull+iextra,kl)
      real d(ifull,kl)   ! NOT shared adjust5 or nonlin
      real savs1(ifull,2:kl),savu1(ifull,kl),savv1(ifull,kl)
      real sbar(ifull,2:kl)
      common/sbar/sbar 
      common/savuv1/savs1,savu1,savv1 
      real omgf(ifull,kl),e(ifull,kl)
      equivalence (omgf,dpsldt)
      real sdotin(ifull,kl),pslxin(ifull,kl),omgfin(ifull,kl)
      real pse(ifull+iextra),psn(ifull+iextra),psz(ifull+iextra)
      real, dimension(ifull+iextra,kl) :: uc, vc, wc
      real pslx_k(kl)
      save num
      data num/0/
!     Always called for first time step
!     Called every step for mup.ne.1
!     Usual is mup=1, using simple centred (only first step)
!              mup=2, using simple centred (every time step)
!              mup=4 to use prior pslx (no longer here)
!              mup=5 (somewhat similar to 1) & 6 calc via staguv
!              mup<0 various others
      
      if(mup<=-4)then
        if(ktau<3)then
          return
        else
          sdotin(:,1:kl)=sdot(:,1:kl)
          omgfin(:,1:kl)=omgf(:,1:kl)
          pslxin(:,1:kl)=pslx(1:ifull,1:kl)
	 endif
      endif  ! (mup<=-4)
      
      if(mup>=5)then
        call staguv(u(1:ifull,:),v(1:ifull,:),
     &             cc(1:ifull,:),dd(1:ifull,:)) 
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
      endif       ! (mup>=5) 
      
      if(mup<5)then       ! ************usual***********
        call boundsuv(u,v)
        do k=1,kl
*cdir nodep
         do iq=1,ifull
!          N.B. this div is calculated on the non-staggered grid
!          but in adjust5 it is on the staggered grid
           d(iq,k)=                      ! usual   
     .        (u(ieu(iq),k)/em(ie(iq))-u(iwu(iq),k)/em(iw(iq))     
     .        +v(inv(iq),k)/em(in(iq))-v(isv(iq),k)/em(is(iq)))
     .        *em(iq)**2/(2.*ds)  
          enddo   ! iq loop
         enddo    ! k  loop
      endif       ! (mup<5) 
	
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
      elseif(mup<4.or.num==0)then  ! ************usual***********
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
         if(diag)then
	    call bounds(ps)
	    print *,'after bounds in updps'
	    if(mydiag)then
 	     iq=idjd
	     print *,'in updps'
	     print *,'dhatA ',(d(iq,k)-pslx(iq,k),k=1,kl)
            do k=1,kl
             pslx_k(k)=-em(iq)*(
     .                 u(iq,k)*(ps(ie(iq))-ps(iw(iq)))
     .                +v(iq,k)*(ps(in(iq))-ps(is(iq))))/(2.*ds*ps(iq))
            enddo    ! k  loop
	     print *,'dhatAA ',(d(iq,k)-pslx_k(k),k=1,kl)
	     print *,'part1AA ',(d(iq,k),k=1,kl)
	     print *,'part2AA ',(-pslx_k(k),k=1,kl)
	    endif
	  endif
       endif       ! (mup<4.or.num==0)

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
      
!------------------------------------------------------------------      
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
	 enddo   ! iq loop
      endif     ! (mup<=-4)

      if(diag.or.nmaxpr==1)then
        if(mydiag)then
          print *,'in updps'
          write (6,"('div5p',5p10f8.2)") d(idjd,:)
          iq=idjd
          k=nlv
          print *,'iq,ie,iw,in,is ',iq,ie(iq),iw(iq),in(iq),is(iq)
          print *,'em_iq,ie,iw,in,is ',
     &             em(iq),em(ie(iq)),em(iw(iq)),em(in(iq)),em(is(iq))
          print *,'iq,ieu,iwu,inv,isv ',
     &             iq,ieu(iq),iwu(iq),inv(iq),isv(iq)
          print *,'u_iq,ieu,iwu ',u(iq,k),u(ieu(iq),k),u(iwu(iq),k)
          print *,'v_iq,inv,isv ',v(iq,k),v(inv(iq),k),v(isv(iq),k)
        endif
        call printa('div5',d,ktau,nlv,ia,ib,ja,jb,0.,1.e5)
        call printa('u',u,ktau,nlv,ia,ib,ja,jb,0.,1.)
        call printa('v',u,ktau,nlv,ia,ib,ja,jb,0.,1.)
      endif
      return
      end
