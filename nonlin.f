      subroutine nonlin      ! globpea version  (only has npex=1)
!     has only morder=24, & nphip options
      parameter (meth=1)    ! 0 or 1 for nphip scheme
      parameter (nphip=1)   ! 1:off, 421:on  1052,2,-2.5
      parameter (ntest=0)
      parameter (nwhite=0)   ! 0:off, 1:on
      parameter (mfix_rad=0) ! used to make gases 2 to ng add up to gas 1
      include 'newmpar.h'
      include 'arrays.h'
      include 'constant.h' ! r,g,cp,cpv,roncp
      include 'indices.h'  ! in,is,iw,ie,inn,iss,iww,iee
      include 'latlong.h'
      include 'map.h'
      include 'morepbl.h'  ! condx
      include 'nlin.h'
      include 'parm.h'
      include 'parmdyn.h'  
      include 'parmvert.h'
 
      include 'sigs.h'
      include 'tracers.h'
      include 'vecsuv.h'
      include 'vvel.h'     ! sdot
      include 'xarrs.h'
      include 'xyzinfo.h'  ! x,y,z
      common/epst/epst(ifull)
      common/neigh/neigh(ifull)
      common/nonlsav/tnsav(ifull,kl),unsav(ifull,kl),vnsav(ifull,kl)
      common/savuv/savu(ifull,kl),savv(ifull,kl)
      common/tbar2d/tbar2d(ifull)
      common/work2/aa(ifull),bb(ifull),cc(ifull),dd(ifull),
     .             aa2(ifull),bb2(ifull),cc2(ifull),dd2(ifull),
     .             ee(ifull),ff(ifull),dum2(ifull,7),pskap(ifull)
      common/work3/d(ifull,kl),p(ifull,kl),tempry(ifull,kl),tv(ifull,kl)
     .             ,spare(ifull,kl)      
      common/work3sav/qgsav(ifull,kl),trsav(ilt*jlt,klt,ngasmax)  ! passed to adjust5
      real pextras(ifull,kl),omgf(ifull,kl)
      equivalence (omgf,pextras,dpsldt)
      real phip(ifull,nphip),dphip(ifull,nphip)    ! 1052 to 2 every 25
      real dphi_dx(ifull,kl),dphi_dy(ifull,kl)
      real siglog(kl),plog(nphip),dplog(nphip)
      data pi/3.1415926536/

      if(epsp.lt.-1.)then
!        e.g. -20. gives epst=.2 for sdmax=1.
        do iq=1,ifull
         sdmax=0.
         do k=2,kl-1
          sdmax=max(sdmax,abs(sdot(iq,k)))
         enddo
         epst(iq)=sdmax*abs(.01*epsp)
        enddo         
      endif

      do k=1,kl
       do iq=1,ifull
!       *** following qgsav should be before first vadv call
        qgsav(iq,k)=qg(iq,k)      ! for qg  conservation in adjust5
!       N.B. [D + dsigdot/dsig] saved in adjust5 as pslx
        pslx(iq,k)=psl(iq)-pslx(iq,k)*dt*.5*(1.-epst(iq))
       enddo     ! iq loop
      enddo      ! k  loop
      if(ngas.ge.1)then
        if(mfix_rad.gt.0.and.mspec.eq.1)then ! make gases 2 to ng add up to g1
         do k=1,kl              ! here it is just from effects of physics
	   do iq=1,ifull
          sumdiffb=0.
           do ng=2,ngas        
            sumdiffb=sumdiffb+tr(iq,k,ng)
            delpos=delpos+max( 1.e-20,tr(iq,k,ng)-trsav(iq,k,ng))
            delneg=delneg+min(-1.e-20,tr(iq,k,ng)-trsav(iq,k,ng))
           enddo   ! ng loop
           ratio=(tr(iq,k,1)-sumdiffb)/(delpos-delneg)
           do ng=2,ngas        
            tr(iq,k,ng)=max(0.,trsav(iq,k,ng)
     .         +(1.+ratio)*max(0.,tr(iq,k,ng)-trsav(iq,k,ng))
     .         +(1.-ratio)*min(0.,tr(iq,k,ng)-trsav(iq,k,ng)) )
           enddo   ! ng loop
          enddo    ! iq loop
         enddo     ! k  loop
        endif      ! (mfix_rad.gt.0)
        do ng=1,ngas   ! re-set trsav prior to vadv, hadv, hordif
         do k=1,kl               
	   do iq=1,ifull
           trsav(iq,k,ng)=tr(iq,k,ng) ! for tr conservation in adjust5
          enddo   ! iq loop
         enddo    ! k  loop
        enddo     ! ng loop
      endif       ! (ngas.ge.1)

      if(diag)then
        print *,'at beginning of nonlin'
        print *,'t ',(t(idjd,k),k=1,kl)
        call printa('t   ',t(1,nlv),ktau,nlv,ia,ib,ja,jb,200.,1.)
        print *,'ps . & e ',ps(idjd),ps(ie(idjd))
        print *,'ps . & n ',ps(idjd),ps(in(idjd))
        print *,'psl . & e ',psl(idjd),psl(ie(idjd))
        print *,'psl . & n ',psl(idjd),psl(in(idjd))
        print *,'emu, emuw ',emu(idjd),emu(iwu2(idjd))
        print *,'emv, emvs ',emv(idjd),emv(isv2(idjd))
        print *,'u ',(u(idjd,k),k=1,kl)
        print *,'uw ',(u(iwu(idjd),k),k=1,kl)
        call printa('u   ',u(1,nlv),ktau,nlv,ia,ib,ja,jb,0.,1.)
        print *,'v ',(v(idjd,k),k=1,kl)
        print *,'vs ',(v(isv(idjd),k),k=1,kl)
        call printa('v   ',v(1,nlv),ktau,nlv,ia,ib,ja,jb,0.,1.)
        print *,'pslx ',(pslx(idjd,k),k=1,kl)
c        print *,'pslx(nlv) ne sw nw se ',pslx(id+1,jd+1,nlv)
c     .      ,pslx(id-1,jd-1,nlv),pslx(id-1,jd+1,nlv),pslx(id+1,jd-1,nlv)
        print *,'tn*dt a0 ',(tn(idjd,k)*dt,k=1,kl)
        print *,'un*dt a0 ',(un(idjd,k)*dt,k=1,kl)
        print *,'vn*dt a0 ',(vn(idjd,k)*dt,k=1,kl)
        print *,'roncp,rata,ratb ',roncp,rata,ratb
        print *,'tbar ',(tbar(k),k=1,kl)
        print *,'sig ',(sig(k),k=1,kl)
      endif

      un(:,:)=0. !   needed (whilst un equiv in vertmix)
      vn(:,:)=0.
      tn(:,:)=0.
      if(diag.or.nmaxpr.eq.1)then
        print *,'in nonlin before vertical advection'
        write (6,"('sdot',9f8.3/4x,9f8.3)") (sdot(idjd,kk),kk=1,kl)
        write (6,"('t   ',9f8.2/4x,9f8.2)") (t(idjd,kk),kk=1,kl)
        write (6,"('u   ',9f8.2/4x,9f8.2)") (u(idjd,kk),kk=1,kl)
        write (6,"('v   ',9f8.2/4x,9f8.2)") (v(idjd,kk),kk=1,kl)
        write (6,"('qg  ',9f8.3/4x,9f8.3)")(1000.*qg(idjd,kk),kk=1,kl)
      endif

      if(nvad.ne.0)then
!       do vertical advection in split mode
        if(nvad.eq.4)call vadvtvd(t,tx,u,ux,v,vx)  ! can now call from globpe too
        if(nvad.eq.7)call vadv30(t,tx,u,ux,v,vx)   ! for vadvbess
        t(:,:)=tx(:,:)
        u(:,:)=ux(:,:)
        v(:,:)=vx(:,:)
        if(diag.or.nmaxpr.eq.1)then
         print *,'in nonlin after vertical advection'
         write (6,"('qg  ',9f8.3/4x,9f8.3)")(1000.*qg(idjd,kk),kk=1,kl)
         write (6,"('t   ',9f8.2/4x,9f8.2)") (t(idjd,kk),kk=1,kl)
         write (6,"('thet',9f8.2/4x,9f8.2)")  
     .                (t(idjd,k)*sig(k)**(-roncp),k=1,kl)
         write (6,"('u   ',9f8.2/4x,9f8.2)") (u(idjd,kk),kk=1,kl)
         write (6,"('v   ',9f8.2/4x,9f8.2)") (v(idjd,kk),kk=1,kl)
         write (6,"('t#  ',9f8.2)") 
     .             ((t(ii+(jj-1)*il,nlv),ii=id-1,id+1),jj=jd-1,jd+1)
         write (6,"('u#  ',9f8.2)") 
     .             ((u(ii+(jj-1)*il,nlv),ii=id-1,id+1),jj=jd-1,jd+1)
         write (6,"('v#  ',9f8.2)") 
     .             ((v(ii+(jj-1)*il,nlv),ii=id-1,id+1),jj=jd-1,jd+1)
         write (6,"('omgf#',9f8.3)") ((ps(ii+(jj-1)*il)*
     .               omgf(ii+(jj-1)*il,nlv),ii=id-1,id+1),jj=jd-1,jd+1)
        endif
      endif      ! (nvad.ne.0)
      if(diag)then
        write (6,"('qg ',19f7.3/(8x,19f7.3))") 
     .             (1000.*qg(idjd,k),k=1,kl)
        if(sig(nlv).lt..3)then
          call printa('qg  ',qg(1,nlv),ktau,nlv,ia,ib,ja,jb,0.,1.e6)
        else
          call printa('qg  ',qg(1,nlv),ktau,nlv,ia,ib,ja,jb,0.,1.e3)
        endif
        print *,'tn*dt a1 ',(tn(idjd,k)*dt,k=1,kl)
        print *,'un*dt a1 ',(un(idjd,k)*dt,k=1,kl)
        print *,'vn*dt a1 ',(vn(idjd,k)*dt,k=1,kl)
      endif

      if(abs(mfix).eq.5.and.mspec.eq.1)then
!       perform conservation fix on tr1,tr2 as affected by vadv, hadv, hordif
        do ng=1,2
         delpos=0.
         delneg=0.
         do iq=1,ifull
          do k=1,kl
           tempry(iq,k)=tr(iq,k,ng)-trsav(iq,k,ng)  ! has increments
           delpos=delpos+max(0.,-dsig(k)*tempry(iq,k)/em(iq)**2)
           delneg=delneg+min(0.,-dsig(k)*tempry(iq,k)/em(iq)**2)
          enddo   ! k loop
         enddo    ! iq loop
         ratio = -delneg/delpos
         beta = min(ratio,sqrt(ratio))
         betav=1./max(1.,beta)  ! for cunning 2-sided
         do k=1,kl
	   do iq=1,ifull
           tr(iq,k,ng)=trsav(iq,k,ng)+
     .       beta*max(0.,tempry(iq,k)) + betav*min(0.,tempry(iq,k))
          enddo  ! iq loop  
         enddo   !  k loop    
        enddo    ! ng loop
      endif      !  abs(mfix).eq.5

      if ( nhstest.eq.1 ) then ! Held and Suarez test case
         call hs_phys
      endif

      if(diag)then
        print *,'tn*dt a ',(tn(idjd,k)*dt,k=1,kl)
        print *,'un*dt a ',(un(idjd,k)*dt,k=1,kl)
        print *,'vn*dt a ',(vn(idjd,k)*dt,k=1,kl)
        print *,'sdot ',(sdot(idjd,k),k=1,kl+1)
        call printa('sdot',sdot(1,nlv+1),ktau,nlv+1,ia,ib,ja,jb,
     .               0.,10.)
        print *,'omgf ',(omgf(idjd,k),k=1,kl)
        call printa('omgf',omgf(1,nlv),ktau,nlv,ia,ib,ja,jb,0.,1.e5)
        do iq=1,ifull
         aa(iq)=rata(nlv)*sdot(iq,nlv+1)+ratb(nlv)*sdot(iq,nlv)
        enddo    ! iq loop
        print *,'k,aa,emu',nlv,aa(idjd),emu(idjd)
        call printa('sgdf',aa,ktau,nlv,ia,ib,ja,jb,0.,10.)
      endif

      tv=.61*qg*t  ! 3D - just add-on at this stage (after vadv)
      contv=(1.61-cpv/cp)/.61      ! about -.26/.61
      if(ntbar.lt.0)then
        do iq=1,ifull
         tbar2d(iq)=t(iq,1)+contv*tv(iq,1)
        enddo   ! iq loop
      endif     ! (ntbar.lt.0)
      if(ntbar.eq.0)then
        do iq=1,ifull
         tbar2d(iq)=tbar(1)
        enddo   ! iq loop
      endif     ! (ntbar.eq.0)
      if(ntbar.gt.0)then
        do iq=1,ifull
         tbar2d(iq)=t(iq,ntbar)
        enddo   ! iq loop
      endif     ! (ntbar.gt.0)

      do k=1,kl
       do iq=1,ifull
        termlin=tbar2d(iq)*omgf(iq,k)*roncp/sig(k) ! full omgf used here
        tn(iq,k)=tn(iq,k)+(t(iq,k)
     .           +contv*tv(iq,k))*omgf(iq,k)*roncp/sig(k) -termlin
        tv(iq,k)=t(iq,k)+tv(iq,k)
!       add in  cnon*dt*tn(iq,k)  term at bottom
        tx(iq,k)=t(iq,k) +.5*dt*(1.-epst(iq))*termlin
       enddo     ! iq loop
      enddo      ! k  loop
      if(diag)then
        iq=idjd
        k=nlv
        tvv=tv(iq,k)-t(iq,k)
        print *,'contv,tbar2d,termlin ',
     .           contv,tbar2d(iq),tbar2d(iq)*omgf(iq,k)*roncp/sig(k)
        print *,'tvv,tn ',tvv,tn(iq,k)
        print *,'termx ',(t(iq,k)+contv*tvv)*omgf(iq,k)*roncp/sig(k)
      endif

      if(nwhite.eq.1)then
        coslat=cos(rlat0*pi/180.)
        sinlat=sin(rlat0*pi/180.)
        polenx=-coslat
        polenz=sinlat
        omg_rot=2.*pi/86400.
        do iq=1,ifull
         factor=2.*omg_rot*cos(rlatt(iq))/g
         zonx=            -polenz*y(iq)
         zony=polenz*x(iq)-polenx*z(iq)
         zonz=polenx*y(iq)
         den=sqrt( max(zonx**2+zony**2+zonz**2,1.e-7) )  ! allow for poles
         costh= (zonx*ax(iq)+zony*ay(iq)+zonz*az(iq))/den
         sinth=-(zonx*bx(iq)+zony*by(iq)+zonz*bz(iq))/den
         do k=1,kl
          uzon= costh*u(iq,k)-sinth*v(iq,k)
          tv(iq,k)=tv(iq,k)+t(iq,k)*uzon*factor
c         if(iq.eq.idjd.and.k.eq.12)print *,'u,v,uzon,factor,prod ',
c    .      u(iq,k),v(iq,k),uzon,factor,t(iq,k)*uzon*factor
         enddo      ! k  loop    
        enddo      ! iq loop
      endif  ! (nwhite.eq.1)

!     phip calculation section
      if(nphip.gt.1)then
       if(ktau.eq.1)print *,
     .  'phip calculation with nphip,ps,zs: ',nphip,ps(idjd),zs(idjd)
       delp=2.5e2  ! e.g. 2.5 hPa
	do k=1,kl
	 siglog(k)=log(sig(k))
	enddo
	do k=1,nphip
	 pressp=1052.e2-(k-1.)*delp  ! in Pa
	 plog(k)=log(pressp)
	enddo
	do k=2,nphip
	 dplog(k)=plog(k)-plog(k-1)
	enddo
	
       do iq=1,ifull
        phip(iq,1)=0.
	 do k=2,nphip
!	  choose temperature at level k-.5
 	  iqq=iq
	  presst=1052.e2-(k-1.5)*delp  ! in Pa
	  if(presst.gt.ps(iq))then
	    iqq=neigh(iq)
	  endif
         if(meth.eq.0)then
	    sigt=max(sig(kl),min(sig(1),presst/ps(iqq)))
           do kk=2,kl
            if(sigt.le.sig(kk-1))then
              kx=kk-1
            endif
           enddo  ! kk loop
           tt=((sigt-sig(kx+1))*tv(iqq,kx)+
     .        (sig(kx)-sigt)*tv(iqq,kx+1))/(sig(kx)-sig(kx+1))
         endif
         if(meth.eq.1)then
	    sigt=max(.00463,min(.99537,presst/ps(iqq)))
           sigxx=sigt*(1.-sigt)
c          rk=kl+.5-kl*sigt - kl*sigxx*(.5-sigt)*(12.5+28.*sigxx)/
c    .        (1.+53.*sigxx)
           rk=kl+.5-kl*sigt - kl*sigxx*(.5-sigt)*(12.54+27.79*sigxx)/
     .        (1.+53.17*sigxx)
           kx=rk
           drk=rk-kx
!          tt=(1.-drk)*tv(iqq,kx)+drk*tv(iqq,kx+1)
c	    tt=((sigt-sig(kx+1))*tv(iqq,kx)+
c	.       (sig(kx)-sigt)*tv(iqq,kx+1))/(sig(kx)-sig(kx+1))
           sigtlog=log(sigt)
           tt=((sigtlog-siglog(kx+1))*tv(iqq,kx)+(siglog(kx)-sigtlog)*
     .                         tv(iqq,kx+1))/(siglog(kx)-siglog(kx+1))
         endif
!        phip(iq,k)=phip(iq,k-1)+287.*tt*delp/presst
	  phip(iq,k)=phip(iq,k-1)-287.*tt*dplog(k)
c	  if(ntest.eq.1.and.iq.eq.idjd.and.k.gt.410)then
c	    press=presst-.5*delp
c	    print *,iqq,k,kx,press,sigt,tt,phip(iq,k),-287.*tt*dplog(k)
c	  endif  ! (ntest.eq.1.and.iq.eq.idjd)
	 enddo   ! k loop
       enddo    !  iq loop
       do iq=1,ifull
	 pp=1.+(1052.e2-ps(iq))/delp
	 kpp=int(pp)
	 zsint=(pp-kpp)*phip(iq,kpp+1)+(kpp+1-pp)*phip(iq,kpp)
	 do k=1,nphip
	  phip(iq,k)=phip(iq,k)+zs(iq)-zsint
	 enddo  ! k loop	  
c	 if(ntest.eq.1.and.iq.eq.idjd)then
c	   print *,'ktau,pp,kpp,zsint ',ktau,pp,kpp,zsint
c	   do k=nphip-10,nphip
c	    print *,k,phip(iq,k)
c	   enddo  ! k loop	  
c	 endif    ! (ntest.eq.1.and.iq.eq.idjd)
       enddo     !  iq loop
!      dphi_dx calculation
       do k=1,nphip
        do iq=1,ifull
	  dphip(iq,k)=phip(ie(iq),k)-phip(iq,k)
        enddo  !  iq loop
       enddo   ! k loop
       do iq=1,ifull
        psav=.5*(ps(iq)+ps(ie(iq)))
!	 psavlog=log(psav)                               ! choice 1
        psavlog=log(1.e5)+.5*(psl(iq)+psl(ie(iq)))  ! choice 2
        do k=1,kl    ! this one linear interpolation
  	  psavk=psav*sig(k)
         pp=1.+(1052.e2-psavk)/delp
	  kpp=int(pp)
	  psavklog=psavlog+siglog(k)
!	  next line for linear interp
c	  dphi_dx(iq,k)=(pp-kpp)*dphip(iq,kpp+1)+(kpp+1-pp)*dphip(iq,kpp)
!	  next lines for log interp
	  dphi_dx(iq,k)=((psavklog-plog(kpp))*dphip(iq,kpp+1)+
     . 	  (plog(kpp+1)-psavklog)*dphip(iq,kpp))/dplog(kpp+1)
c	  if(ntest.eq.1.and.iq.eq.idjd)then
c	    print *,'k,psavk,pp,kpp ',k,psavk,pp,kpp
c	  endif ! (ntest.eq.1.and.iq.eq.idjd)
        enddo  ! k loop
       enddo   !  iq loop
!      dphi_dy calculation
       do k=1,nphip
        do iq=1,ifull
	  dphip(iq,k)=phip(in(iq),k)-phip(iq,k)
        enddo  !  iq loop
       enddo  ! k loop
       do iq=1,ifull
        psav=.5*(ps(iq)+ps(in(iq)))
!	 psavlog=log(psav)                           ! choice 1
        psavlog=log(1.e5)+.5*(psl(iq)+psl(in(iq)))  ! choice 2
        do k=1,kl    ! this one linear interpolation
	  psavk=psav*sig(k)
         pp=1.+(1052.e2-psavk)/delp
	  kpp=int(pp)
	  psavklog=psavlog+siglog(k)
!	  next line for linear interp
c	  dphi_dy(iq,k)=(pp-kpp)*dphip(iq,kpp+1)+(kpp+1-pp)*dphip(iq,kpp)
!	  next lines for log interp
	  dphi_dy(iq,k)=((psavklog-plog(kpp))*dphip(iq,kpp+1)+
     . 	  (plog(kpp+1)-psavklog)*dphip(iq,kpp))/dplog(kpp+1)
        enddo  ! k loop
       enddo   !  iq loop
	if(ntest.eq.1)then
  	  iq=idjd
	  print *,'ps,pse,psn ',ktau,ps(iq),ps(ie(iq)),ps(in(iq))
	  print *,'psl,psle,psln ',psl(iq),psl(ie(iq)),psl(in(iq))
	  print *,'phitopm,e,n ',phip(iq,nphip-4),phip(ie(iq),nphip-4),
     .                          phip(in(iq),nphip-4)
	  print *,'phitop,e,n ',phip(iq,nphip),phip(ie(iq),nphip),
     .                         phip(in(iq),nphip)
         do k=1,kl
	   print *,'k,tv,tve,tvn ',
     .             k,tv(iq,k),tv(ie(iq),k),tv(in(iq),k)
	  enddo
	endif   ! (ntest.eq.1)
      endif    ! (nphip.gt.1)
      	   
!     only has code with npex=1 higher order grad expressions
!     calculate "full-linear" geopotential height terms using tv and save in p
      do iq=1,ifull
       p(iq,1)=zs(iq)+bet(1)*tv(iq,1)
      enddo    ! iq loop
      do k=2,kl
       do iq=1,ifull
        p(iq,k)=p(iq,k-1)+bet(k)*tv(iq,k)+betm(k)*tv(iq,k-1)
       enddo   ! iq loop
      enddo    ! k  loop

      if(ntest.eq.1.and.nphip.gt.1)then
       print *,ktau,'dphi_dx  & (approx.) usual way'
       do k=1,kl
        rtav=.5*287.*(tv(idjd,k)+tv(ie(idjd),k))
        print *,ktau,'dphidx',k,
     .                dphi_dx(idjd,k),p(ie(idjd),k)-p(idjd,k)
     .                +rtav*(psl(ie(idjd))-psl(idjd))
       enddo    ! k  loop
       print *,ktau,'dphi_dy  & (approx.) usual way'
       do k=1,kl
        rtav=.5*287.*(tv(idjd,k)+tv(in(idjd),k))
        print *,ktau,'dphidy',k,
     .                dphi_dy(idjd,k),p(in(idjd),k)-p(idjd,k)
     .                +rtav*(psl(in(idjd))-psl(idjd))
       enddo    ! k  loop
      endif     ! (ntest.eq.1.and.nphip.gt.1)

      do k=1,kl
       do iq=1,ifull
        p(iq,k)=p(iq,k)+r*tv(iq,k)*psl(iq)
       enddo     ! iq loop
      enddo      ! k  loop

!     calculate "basic-linear" geopotential height and put in tempry
      do iq=1,ifull            ! tempry used by pextras
       tempry(iq,1)=zs(iq)+bet(1)*t(iq,1)+r*tbar2d(iq)*psl(iq)
      enddo      ! iq loop
      do k=2,kl
       do iq=1,ifull
        tempry(iq,k)=tempry(iq,k-1)
     .                 +bet(k)*t(iq,k)+betm(k)*t(iq,k-1)
       enddo     ! iq loop
      enddo      ! k  loop
!     save "full-linear" - "basic-linear" geopotential height in pextras
      do k=1,kl
       do iq=1,ifull
        pextras(iq,k)=p(iq,k)-tempry(iq,k)
       enddo  ! iq loop  
      enddo   !  k loop    
!     N.B. these pextras terms (actually their horiz derivs) could have
!     been bundled in here with un, vn but they would have to be staggered
!     later

!     Now set up in ux,vx the tendencies from "full-linear" terms.
!     The following un, vn contributions include residual psl.grad(T) terms
!     Following now always basically morder=24 2nd order, Corby style scheme
      do k=1,kl
!cdir nodep
       do iq=1,ifull  ! calculate staggered contributions first
        aa(iq)=-emu(iq)*(p(ie(iq),k)-p(iq,k))*
     .                                      (1.-epsu)*.5*dt/ds    ! 2nd order
        bb(iq)=-emv(iq)*(p(in(iq),k)-p(iq,k) )*
     .                                      (1.-epsu)*.5*dt/ds    ! 2nd order
       enddo   ! iq loop
!      morder=24 scheme with residual terms for un, vn
!cdir nodep
         do iq=1,ifull  ! calculate staggered contributions first
          cc(iq)=emu(iq)*(psl(ie(iq))+psl(iq))*
     .                          (tv(ie(iq),k)-tv(iq,k))*.5*r/ds
          dd(iq)=emv(iq)*(psl(in(iq))+psl(iq))*
     .                          (tv(in(iq),k)-tv(iq,k))*.5*r/ds
         enddo   ! iq loop

!cdir nodep
       if(nphip.gt.1.and.k.le.1)then
         do iq=1,ifull  ! calculate nphip contributions 
          cc(iq)=-emu(iq)*(dphi_dx(iq,k)
     .                     -p(ie(iq),k)+p(iq,k))/ds
          dd(iq)=-emv(iq)*(dphi_dy(iq,k)
     .                     -p(in(iq),k)+p(iq,k))/ds
         enddo   ! iq loop
	  iq=idjd
	  print *,"phi's ",k,dphi_dx(iq,k),p(ie(iq),k)-p(iq,k),
     .                      dphi_dy(iq,k),p(in(iq),k)-p(iq,k)
       endif     ! (nphip.gt.1)
       call unstaguv(aa,bb,aa2,bb2) ! convert to unstaggered positions
       call unstaguv(cc,dd,cc2,dd2)
       do iq=1,ifull   
        ux(iq,k)=u(iq,k)+aa2(iq)
        vx(iq,k)=v(iq,k)+bb2(iq)
        un(iq,k)=un(iq,k)+cc2(iq)
        vn(iq,k)=vn(iq,k)+dd2(iq)
       enddo  ! iq loop
       if(diag.and.k.eq.nlv)then
        call printa('aa  ',aa,ktau,nlv,ia,ib,ja,jb,0.,1.)
        call printa('aa2 ',aa2,ktau,nlv,ia,ib,ja,jb,0.,1.)
        call printa('bb  ',bb,ktau,nlv,ia,ib,ja,jb,0.,1.)
        call printa('bb2 ',bb2,ktau,nlv,ia,ib,ja,jb,0.,1.)
        call printa('cc  ',cc,ktau,nlv,ia,ib,ja,jb,0.,1.)
        call printa('cc2 ',cc2,ktau,nlv,ia,ib,ja,jb,0.,1.)
        call printa('dd  ',dd,ktau,nlv,ia,ib,ja,jb,0.,1.)
        call printa('dd2 ',dd2,ktau,nlv,ia,ib,ja,jb,0.,1.)
       endif  ! (diag.and.k.eq.nlv)
      enddo   ! k  loop

!     finish evaluation of tx,ux,vx by adding in part of nonlinear terms
      cnon=.5    ! cnon=.5 shares tn,un,vn between here & upglobal/adjust
      if(nonl.eq.1.or.nonl.eq.11)cnon=1.   ! cnon=1. all done here
      if(nonl.eq.2)then
        tn(:,:)=2.*tn(:,:)   ! all tn added in adjust
        un(:,:)=2.*un(:,:)   ! all un added in adjust
        vn(:,:)=2.*vn(:,:)   ! all vn added in adjust
      elseif(nonl.eq.3)then   ! Adams-Bashforth style
        if(ktau.eq.1)then
          tnsav(:,:)=tn(:,:)
          unsav(:,:)=un(:,:)
          vnsav(:,:)=vn(:,:)
        endif
        tx(:,:)=tx(:,:) +1.5*dt*tn(:,:)-.5*dt*tnsav(:,:)
        ux(:,:)=ux(:,:) +1.5*dt*un(:,:)-.5*dt*unsav(:,:)
        vx(:,:)=vx(:,:) +1.5*dt*vn(:,:)-.5*dt*vnsav(:,:)
      else       ! (nonl=0, 1, or 11)  usual code has nonl=0
        tx(:,:)=tx(:,:) +cnon*dt*tn(:,:)
        ux(:,:)=ux(:,:) +cnon*dt*un(:,:)
        vx(:,:)=vx(:,:) +cnon*dt*vn(:,:)
      endif

      if(diag)then
        print *,'at end of nonlin; nvad,idjd = ', nvad,idjd
        print *,'p1 . & e ',p(idjd,1),p(ie(idjd),1)
        print *,'p1 . & n ',p(idjd,1),p(in(idjd),1)
        print *,'psl . & e ',psl(idjd),psl(ie(idjd))
        print *,'psl . & n ',psl(idjd),psl(in(idjd))
        call printa('psl ',psl,ktau,0,ia,ib,ja,jb,0.,100.)
        print  *,'pslx -1,0,1 ',(pslx(idjd+n,nlv),n=-1,1)
        call printa('pslx',pslx(1,nlv),ktau,nlv,ia,ib,ja,jb,0.,100.)
        print *,'tn*dt b ',(tn(idjd,k)*dt,k=1,kl)
        call printa('tn  ',tn(1,nlv),ktau,nlv,ia,ib,ja,jb,0.,100.*dt)
        print *,'un*dt b ',(un(idjd,k)*dt,k=1,kl)
        call printa('un  ',un(1,nlv),ktau,nlv,ia,ib,ja,jb,0.,100.*dt)
        print *,'vn*dt b ',(vn(idjd,k)*dt,k=1,kl)
        call printa('vn  ',vn(1,nlv),ktau,nlv,ia,ib,ja,jb,0.,100.*dt)
        print *,'tx ',(tx(idjd,k),k=1,kl)
        call printa('tx  ',tx(1,nlv),ktau,nlv,ia,ib,ja,jb,200.,1.)
        print *,'ux ',(ux(idjd,k),k=1,kl)
        call printa('ux  ',ux(1,nlv),ktau,nlv,ia,ib,ja,jb,0.,1.)
        print *,'vx ',(vx(idjd,k),k=1,kl)
        call printa('vx  ',vx(1,nlv),ktau,nlv,ia,ib,ja,jb,0.,1.)
      endif

      if(nonl.eq.1.or.nonl.eq.11)then   ! option A
        un(:,:)=0. 
        vn(:,:)=0.
        tn(:,:)=0.
      endif
      return
      end
