      subroutine nonlin      ! globpea version  (only has npex=1)
      use cc_mpi
      use diag_m
      implicit none
!     has only morder=24, & nphip options
      integer, parameter :: meth=1     ! 0 or 1 for nphip scheme
      integer, parameter :: npgf=1     ! 0 or 2 for test scheme
      integer, parameter :: nphip=1    ! 1:off, 421:on  1052,2,-2.5
      integer, parameter :: ntest=0
      integer, parameter :: nwhite=0   ! 0:off, 1:on
      integer, parameter :: mfix_rad=0 ! used to make gases 2 to ng add up to gas 1
      include 'newmpar.h'
      include 'arrays.h'
      include 'const_phys.h' ! r,g,cp,cpv,roncp
      include 'indices.h'  ! in,is,iw,ie,inn,iss,iww,iee
      include 'kuocom.h'   ! ldr
      include 'liqwpar.h'  ! qfg,qlg
      include 'latlong.h'
      include 'map.h'
      include 'morepbl.h'  ! condx
      include 'nlin.h'
      include 'parm.h'
      include 'parmdyn.h'  
      include 'parmvert.h'
      include 'savuvt.h'
      include 'sigs.h'
      include 'tracers.h'
      include 'vecsuv.h'
      include 'vvel.h'     ! sdot
      include 'xarrs.h'
      include 'xyzinfo.h'  ! x,y,z
      include 'mpif.h'
      real epst
      common/epst/epst(ifull)
      integer neigh
      common/neigh/neigh(ifull)
      real tnsav, unsav, vnsav, tnsavv, unsavv, vnsavv
      common/nonlsav/tnsav(ifull,kl),unsav(ifull,kl),vnsav(ifull,kl)
     &              ,tnsavv(ifull,kl),unsavv(ifull,kl),vnsavv(ifull,kl)
      real tbar2d
      common/tbar2d/tbar2d(ifull)
      ! Need work common for these?
      real aa(ifull,kl),bb(ifull,kl),cc(ifull,kl),dd(ifull,kl),
     &     aa2(ifull,kl),bb2(ifull,kl),cc2(ifull,kl),dd2(ifull,kl)
!      common/work2/aa(ifull),bb(ifull),cc(ifull),dd(ifull),
!     .             aa2(ifull),bb2(ifull),cc2(ifull),dd2(ifull),
!     .             ee(ifull),ff(ifull),dum2(ifull,7),pskap(ifull)
      real p(ifull+iextra,kl),tempry(ifull,kl),
     &     tv(ifull+iextra,kl)
      real d, pexx, spare
      common/work3/d(ifull,kl),pexx(ifull+iextra,kl),
     .             spare(ifull*3*kl-iextra*kl)   ! From updps
      real qgsav, qfgsav, qlgsav, trsav
      common/work3sav/qgsav(ifull,kl),qfgsav(ifull,kl),qlgsav(ifull,kl)
     &             ,trsav(ilt*jlt,klt,ngasmax)  ! shared adjust5 & nonlin
      real pextras(ifull,kl),omgf(ifull,kl)
      equivalence (omgf,pextras,dpsldt)
      real phip(ifull+iextra,nphip),dphip(ifull,nphip)    ! 1052 to 2 every 25
      real dphi_dx(ifull,kl),dphi_dy(ifull,kl)
      real siglog(kl),plog(nphip),dplog(nphip)
      integer iq, iqq, k, kk, kpp, kx, ng, ii, jj, its, nits, nvadh_pass
      real betav, cnon, contv, coslat, costh, delneg, delp, delpos, den,
     &     drk, factor, omg_rot, polenx, polenz, pp, pressp, presst,
     &     psav, psavk, psavklog, psavlog, ratio, rk, sdmax, 
     &     sigt, sigtlog, sigxx, sinlat, sinth, sumdiffb, termlin, tt,
     &     tvv, uzon, zonx, zony, zonz, zsint, sdmx, sdmx_g
      real, dimension(ifull+iextra,kl) :: uin, vin
      integer :: ierr
      integer, save :: num = 0
      
      call start_log(nonlin_begin)

      if(epsp<-1.)then
!        e.g. -20. gives epst=.2 for sdmax=1.
        do iq=1,ifull
         sdmax=0.
         do k=2,kl-1
          sdmax=max(sdmax,abs(sdot(iq,k)))
         enddo
         epst(iq)=sdmax*abs(.01*epsp)
        enddo         
      endif

!     *** following qgsav should be before first vadv call
      qgsav(1:ifull,:)=qg(1:ifull,:)      ! for qg  conservation in adjust5
      if(ldr.ne.0)then
        qfgsav(1:ifull,:)=qfg(1:ifull,:)
        qlgsav(1:ifull,:)=qlg(1:ifull,:)
      endif   ! (ldr.ne.0)

      if(ngas>=1)then
        if(mfix_rad>0.and.mspec==1)then ! make gases 2 to ng add up to g1
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
        endif      ! (mfix_rad>0)
        do ng=1,ngas   ! re-set trsav prior to vadv, hadv, hordif
         do k=1,kl               
	   do iq=1,ifull
           trsav(iq,k,ng)=tr(iq,k,ng) ! for tr conservation in adjust5
          enddo   ! iq loop
         enddo    ! k  loop
        enddo     ! ng loop
      endif       ! (ngas>=1)

      if(nuvfilt.ne.0)then
        uin(1:ifull,:) = u(1:ifull,:)
        vin(1:ifull,:) = v(1:ifull,:)
        call boundsuv(uin,vin,nrows=2)
	 if(nuvfilt==1)then
         do k=1,kl
!cdir nodep
         do iq=1,ifull  
	   if(((uin(iwu(iq),k)<uin(iq,k)).and.
     &        (uin(iq,k)>uin(ieu(iq),k)).and.
     &        (uin(ieu(iq),k)<uin(ieeu(iq),k))).or.
     &       ((uin(iwu(iq),k)>uin(iq,k)).and.
     &        (uin(iq,k)<uin(ieu(iq),k)).and.
     &        (uin(ieu(iq),k)>uin(ieeu(iq),k))))then
            u(iq,k)=.5*(uin(iq,k)+uin(ieu(iq),k))
            u(ieu(iq),k)=.5*(uin(iq,k)+uin(ieu(iq),k))
	   endif
	   if(((vin(isv(iq),k)<vin(iq,k)).and.
     &        (vin(iq,k)>vin(inv(iq),k)).and.
     &        (vin(inv(iq),k)<vin(innv(iq),k))).or.
     &       ((vin(isv(iq),k)>vin(iq,k)).and.
     &        (vin(iq,k)<vin(inv(iq),k)).and.
     &        (vin(inv(iq),k)>vin(innv(iq),k))))then
            v(iq,k)=.5*(vin(iq,k)+vin(inv(iq),k))
            v(inv(iq),k)=.5*(vin(iq,k)+vin(inv(iq),k))
	   endif
	  enddo  ! iq loop
!        following loops only really needed for MPI version	  
	  do jj=1,il
	   iq=1+(jj-1)*il  ! just for i=1 in other dirn
	   if(((uin(ieu(iq),k)<uin(iq,k)).and.
     &        (uin(iq,k)>uin(iwu(iq),k)).and.
     &        (uin(iwu(iq),k)<uin(iwwu(iq),k))).or.
     &       ((uin(ieu(iq),k)>uin(iq,k)).and.
     &        (uin(iq,k)<uin(iwu(iq),k)).and.
     &        (uin(ieu(iq),k)>uin(iwwu(iq),k))))then
            u(iq,k)=.5*(uin(iq,k)+uin(iwu(iq),k))
	   endif
	  enddo
	  do iq=1,il  ! just doing j=1 in other dirn
	   if(((vin(inv(iq),k)<vin(iq,k)).and.
     &        (vin(iq,k)>vin(isv(iq),k)).and.
     &        (vin(isv(iq),k)<vin(issv(iq),k))).or.
     &       ((vin(inv(iq),k)>vin(iq,k)).and.
     &        (vin(iq,k)<vin(isv(iq),k)).and.
     &        (vin(isv(iq),k)>vin(issv(iq),k))))then
            v(iq,k)=.5*(vin(iq,k)+vin(isv(iq),k))
	   endif
	  enddo
	 enddo   ! k  loop
	 endif   ! (nuvfilt==1)then
	 if(nuvfilt>1)then
         do k=nuvfilt,kl  ! e.g. 5
!cdir nodep
         do iq=1,ifull  
          u(iq,k)=(-uin(iwwu(iq),k)+4.*uin(iwu(iq),k)+10.*uin(iq,k)
     &          +4.*uin(ieu(iq),k)-uin(ieeu(iq),k))/16.
          v(iq,k)=(-vin(issv(iq),k)+4.*vin(isv(iq),k)+10.*vin(iq,k)
     &          +4.*vin(inv(iq),k)-vin(innv(iq),k))/16.
	  enddo  ! iq loop
	 enddo   ! k  loop
	 endif   ! (nuvfilt>1)then
      endif     ! (nuvfilt>0)

      ! Diagnostics and neigh calculation below require this.
      if(diag.or.nphip>0)call bounds(ps)
      if(nud_p>0)call bounds(psl)    ! only needed with nudging of psl

      if(diag)then
         if ( mydiag ) then
            print *,'at beginning of nonlin'
            print *,'t ',(t(idjd,k),k=1,kl)
         end if
         call printa('t   ',t,ktau,nlv,ia,ib,ja,jb,200.,1.)
         if ( mydiag ) then
            print *,'ps . & e ',ps(idjd),ps(ie(idjd))
            print *,'ps . & n ',ps(idjd),ps(in(idjd))
            print *,'psl . & e ',psl(idjd),psl(ie(idjd))
            print *,'psl . & n ',psl(idjd),psl(in(idjd))
            print *,'emu, emuw ',emu(idjd),emu(iwu(idjd))
            print *,'emv, emvs ',emv(idjd),emv(isv(idjd))
            print *,'u ',(u(idjd,k),k=1,kl)
            print *,'uw ',(u(iwu(idjd),k),k=1,kl)
         end if
         call printa('u   ',u,ktau,nlv,ia,ib,ja,jb,0.,1.)
         if ( mydiag ) then
            print *,'v ',(v(idjd,k),k=1,kl)
            print *,'vs ',(v(isv(idjd),k),k=1,kl)
         end if
         call printa('v   ',v,ktau,nlv,ia,ib,ja,jb,0.,1.)
         if ( mydiag ) then
            print *,'tn*dt a0 ',(tn(idjd,k)*dt,k=1,kl)
            print *,'un*dt a0 ',(un(idjd,k)*dt,k=1,kl)
            print *,'vn*dt a0 ',(vn(idjd,k)*dt,k=1,kl)
            print *,'roncp,rata,ratb ',roncp,rata,ratb
            print *,'tbar ',(tbar(k),k=1,kl)
            print *,'sig ',(sig(k),k=1,kl)
         end if
      endif

      if( (diag.or.nmaxpr==1) .and. mydiag )then
        print *,'in nonlin before possible vertical advection'
        write (6,"('divn ',9f8.2/4x,9f8.2)") (d(idjd,kk)*1.e6,kk=1,kl)
        write (6,"('sdotn',9f8.3/5x,9f8.3)") (sdot(idjd,kk),kk=1,kl)
        write (6,"('omgfn',9f8.3/5x,9f8.3)")
     .            (ps(idjd)*omgf(idjd,kk),kk=1,kl)
        write (6,"('t   ',9f8.2/4x,9f8.2)") (t(idjd,kk),kk=1,kl)
        write (6,"('u   ',9f8.2/4x,9f8.2)") (u(idjd,kk),kk=1,kl)
        write (6,"('v   ',9f8.2/4x,9f8.2)") (v(idjd,kk),kk=1,kl)
        write (6,"('qg  ',9f8.3/4x,9f8.3)")(1000.*qg(idjd,kk),kk=1,kl)
      endif

!     do vertical advection in split mode
      if(nvad==4.or.nvad==9)then
         sdmx = maxval(abs(sdot))
         call MPI_AllReduce(sdmx, sdmx_g, 1, MPI_REAL, MPI_MAX,
     &                      MPI_COMM_WORLD, ierr )
	 nits=1+sdmx_g/nvadh
	 nvadh_pass=nvadh*nits
	 if(mydiag.and.mod(ktau,nmaxpr)==0)
     &      print *,'in nonlin sdmx,nits,nvadh_pass ',
     &                         sdmx_g,nits,nvadh_pass
         do its=1,nits
            call vadvtvd(t(1:ifull,:),u(1:ifull,:),v(1:ifull,:),
     &                   nvadh_pass) 
	 enddo
      endif  ! (nvad==4.or.nvad==9)

      if(nvad>=7)then
         call vadv30(t(1:ifull,:),u(1:ifull,:),v(1:ifull,:))  ! for vadvbess
      end if

cx      do k=1,kl  ! following done in upglobal from 04/09
cx!       N.B. [D + dsigdot/dsig] saved in adjust5 (or updps) as pslx
cx        pslx(1:ifull,k)=psl(1:ifull)-pslx(1:ifull,k)*dt*.5*(1.-epst(:)) !ca
cx      enddo      ! k  loop

      if(nvad>0.and.(diag.or.nmaxpr==1).and.mydiag)then
       print *,'in nonlin after vertical advection'
       write (6,"('t   ',9f8.2/4x,9f8.2)") (t(idjd,kk),kk=1,kl)
       write (6,"('thet',9f8.2/4x,9f8.2)")  
     .              (t(idjd,k)*sig(k)**(-roncp),k=1,kl)
       write (6,"('qg  ',9f8.3/4x,9f8.3)")(1000.*qg(idjd,kk),kk=1,kl)
       write (6,"('ql ',9f8.3/4x,9f8.3)")(1000.*qlg(idjd,kk),kk=1,kl)
       write (6,"('qf ',9f8.3/4x,9f8.3)")(1000.*qfg(idjd,kk),kk=1,kl)
       write (6,"('u   ',9f8.2/4x,9f8.2)") (u(idjd,kk),kk=1,kl)
       write (6,"('v   ',9f8.2/4x,9f8.2)") (v(idjd,kk),kk=1,kl)
       write (6,"('t#  ',9f8.2)") 
     .           ((t(ii+jj*il,nlv),ii=idjd-1,idjd+1),jj=1,-1,-1)
       write (6,"('qg# ',9f8.3)") 
     .          ((1000.*qg(ii+jj*il,nlv),ii=idjd-1,idjd+1),jj=1,-1,-1)
       write (6,"('ql# ',9f8.3)") 
     .          ((1000.*qlg(ii+jj*il,nlv),ii=idjd-1,idjd+1),jj=1,-1,-1)
       write (6,"('qf# ',9f8.3)") 
     .          ((1000.*qfg(ii+jj*il,nlv),ii=idjd-1,idjd+1),jj=1,-1,-1)
       write (6,"('u#  ',9f8.2)") 
     .           ((u(ii+jj*il,nlv),ii=idjd-1,idjd+1),jj=1,-1,-1)
       write (6,"('v#  ',9f8.2)") 
     .           ((v(ii+jj*il,nlv),ii=idjd-1,idjd+1),jj=1,-1,-1)
       write (6,"('omgf#',9f8.3)") ((ps(ii+(jj-1)*il)*
     .             omgf(ii+jj*il,nlv),ii=idjd-1,idjd+1),jj=1,-1,-1)
        print *,'pslx ',(pslx(idjd,k),k=1,kl)
      endif  ! (nvad>0.and.(diag.or.nmaxpr==1))
      if(diag)then
         if ( mydiag ) write (6,"('qg ',12f7.3/(8x,12f7.3))") 
     &             (1000.*qg(idjd,k),k=1,kl)
         if(sig(nlv)<.3)then
            call printa('qg  ',qg,ktau,nlv,ia,ib,ja,jb,0.,1.e6)
         else
            call printa('qg  ',qg,ktau,nlv,ia,ib,ja,jb,0.,1.e3)
         endif
         if ( mydiag ) then
            print *,'tn*dt a1 ',(tn(idjd,k)*dt,k=1,kl)
            print *,'un*dt a1 ',(un(idjd,k)*dt,k=1,kl)
            print *,'vn*dt a1 ',(vn(idjd,k)*dt,k=1,kl)
         end if
      endif

      if ( nhstest==1 ) then ! Held and Suarez test case
         call hs_phys
      endif

      if(diag)then
         if ( mydiag ) then
            print *,'tn*dt a ',(tn(idjd,k)*dt,k=1,kl)
            print *,'un*dt a ',(un(idjd,k)*dt,k=1,kl)
            print *,'vn*dt a ',(vn(idjd,k)*dt,k=1,kl)
            print *,'sdot ',(sdot(idjd,k),k=1,kl+1)
         end if
         call printa('sdot',sdot,ktau,nlv+1,ia,ib,ja,jb,0.,10.)
         if ( mydiag ) print *,'omgf ',(omgf(idjd,k),k=1,kl)
         call printa('omgf',omgf,ktau,nlv,ia,ib,ja,jb,0.,1.e5)
         do iq=1,ifull
            aa(iq,1)=rata(nlv)*sdot(iq,nlv+1)+ratb(nlv)*sdot(iq,nlv)
         enddo
         if ( mydiag )
     &      print *,'k,aa,emu,emv',nlv,aa(idjd,nlv),emu(idjd),emv(idjd)

         call printa('sgdf',aa(:,1),ktau,nlv,ia,ib,ja,jb,0.,10.)
      endif

!     extra qfg & qlg terms included in tv from April 04
      tv(1:ifull,:) = (.61*qg(1:ifull,:)-qfg(1:ifull,:)-qlg(1:ifull,:))*
     &                t(1:ifull,:)         ! just add-on at this stage 
      contv=(1.61-cpv/cp)/.61      ! about -.26/.61
      if(ntbar==-1.or.(ntbar==-2.and.num==0))then
        do iq=1,ifull
         tbar2d(iq)=t(iq,1)+contv*tv(iq,1)
        enddo   ! iq loop
        num=1
      endif     ! (ntbar==-1.or....)
      if(ntbar==0)then
        do iq=1,ifull
         tbar2d(iq)=tbar(1)
        enddo   ! iq loop
      endif     ! (ntbar==0)
      if(ntbar>0)then
        do iq=1,ifull
         tbar2d(iq)=t(iq,ntbar)
        enddo   ! iq loop
      endif     ! (ntbar>0)
      if(ntbar==-3)then
        do iq=1,ifull
         tbar2d(iq)=max(t(iq,1),t(iq,2),t(iq,3),t(iq,kl))
        enddo   ! iq loop
      endif     ! (ntbar>0)
cy      tx(1:ifull,:)=t(1:ifull,:)  ! cb

      do k=1,kl
       do iq=1,ifull
        termlin=tbar2d(iq)*omgf(iq,k)*roncp/sig(k) ! full omgf used here
        tn(iq,k)=tn(iq,k)+(t(iq,k)
     .           +contv*tv(iq,k))*omgf(iq,k)*roncp/sig(k) -termlin
        tv(iq,k)=t(iq,k)+tv(iq,k)
!       add in  cnon*dt*tn(iq,k)  term at bottom
cy        tx(iq,k)=tx(iq,k) +.5*dt*(1.-epst(iq))*termlin  ! cb
        tx(iq,k)=.5*dt*termlin  ! t and epst later  cy
       enddo     ! iq loop
      enddo      ! k  loop
      if(diag.and.mydiag)then
        iq=idjd
        k=nlv
        tvv=tv(iq,k)-t(iq,k)
        print *,'contv,tbar2d,termlin ',
     .           contv,tbar2d(iq),tbar2d(iq)*omgf(iq,k)*roncp/sig(k)
        print *,'tvv,tn ',tvv,tn(iq,k)
        print *,'termx ',(t(iq,k)+contv*tvv)*omgf(iq,k)*roncp/sig(k)
      endif

      if(nwhite==1)then
        coslat=cos(rlat0*pi/180.)
        sinlat=sin(rlat0*pi/180.)
        polenx=-coslat
        polenz=sinlat
        omg_rot=2.*pi/86400.
        do iq=1,ifull
         factor=2.*omg_rot*cos(rlatt(iq))/grav
         zonx=            -polenz*y(iq)
         zony=polenz*x(iq)-polenx*z(iq)
         zonz=polenx*y(iq)
         den=sqrt( max(zonx**2+zony**2+zonz**2,1.e-7) )  ! allow for poles
         costh= (zonx*ax(iq)+zony*ay(iq)+zonz*az(iq))/den
         sinth=-(zonx*bx(iq)+zony*by(iq)+zonz*bz(iq))/den
         do k=1,kl
          uzon= costh*u(iq,k)-sinth*v(iq,k)
          tv(iq,k)=tv(iq,k)+t(iq,k)*uzon*factor
c         if(iq==idjd.and.k==12)print *,'u,v,uzon,factor,prod ',
c    .      u(iq,k),v(iq,k),uzon,factor,t(iq,k)*uzon*factor
         enddo      ! k  loop    
        enddo      ! iq loop
      endif  ! (nwhite==1)

      call bounds(tv)

!     phip calculation section
      if(nphip>1)then
       if(ktau==1.and.mydiag)print *,
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
	  if(presst>ps(iq))then
	    iqq=neigh(iq)
	  endif
         if(meth==0)then
	    sigt=max(sig(kl),min(sig(1),presst/ps(iqq)))
           do kk=2,kl
            if(sigt<=sig(kk-1))then
              kx=kk-1
            endif
           enddo  ! kk loop
           tt=((sigt-sig(kx+1))*tv(iqq,kx)+
     .        (sig(kx)-sigt)*tv(iqq,kx+1))/(sig(kx)-sig(kx+1))
         endif
         if(meth==1)then
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
!        phip(iq,k)=phip(iq,k-1)+rdry*tt*delp/presst
	  phip(iq,k)=phip(iq,k-1)-rdry*tt*dplog(k)
c	  if(ntest==1.and.iq==idjd.and.k>410)then
c	    press=presst-.5*delp
c	    print *,iqq,k,kx,press,sigt,tt,phip(iq,k),-rdry*tt*dplog(k)
c	  endif  ! (ntest==1.and.iq==idjd)
	 enddo   ! k loop
       enddo    !  iq loop
       do iq=1,ifull
	 pp=1.+(1052.e2-ps(iq))/delp
	 kpp=int(pp)
	 zsint=(pp-kpp)*phip(iq,kpp+1)+(kpp+1-pp)*phip(iq,kpp)
	 do k=1,nphip
	  phip(iq,k)=phip(iq,k)+zs(iq)-zsint
	 enddo  ! k loop	  
c	 if(ntest==1.and.iq==idjd)then
c	   print *,'ktau,pp,kpp,zsint ',ktau,pp,kpp,zsint
c	   do k=nphip-10,nphip
c	    print *,k,phip(iq,k)
c	   enddo  ! k loop	  
c	 endif    ! (ntest==1.and.iq==idjd)
       enddo     !  iq loop

       ! Required by strict shape matching with APAC compiler.
       call bounds(phip(:,1))  ! Assumes nphip=1

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
c	  if(ntest==1.and.iq==idjd)then
c	    print *,'k,psavk,pp,kpp ',k,psavk,pp,kpp
c	  endif ! (ntest==1.and.iq==idjd)
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
	if(ntest==1.and.mydiag)then
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
	endif   ! (ntest==1)
      endif    ! (nphip>1)
      	   
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

!        Get rid of these diagnostics because they would require extra
!        bounds(p) call
!!!      if(ntest==1.and.nphip>1.and.mydiag)then
!!!       print *,ktau,'dphi_dx  & (approx.) usual way'
!!!       do k=1,kl
!!!        rtav=.5*rdry*(tv(idjd,k)+tv(ie(idjd),k))
!!!        print *,ktau,'dphidx',k,
!!!     .                dphi_dx(idjd,k),p(ie(idjd),k)-p(idjd,k)
!!!     .                +rtav*(psl(ie(idjd))-psl(idjd))
!!!       enddo    ! k  loop
!!!       print *,ktau,'dphi_dy  & (approx.) usual way'
!!!       do k=1,kl
!!!        rtav=.5*rdry*(tv(idjd,k)+tv(in(idjd),k))
!!!        print *,ktau,'dphidy',k,
!!!     .                dphi_dy(idjd,k),p(in(idjd),k)-p(idjd,k)
!!!     .                +rtav*(psl(in(idjd))-psl(idjd))
!!!       enddo    ! k  loop
!!!      endif     ! (ntest==1.and.nphip>1)

      do k=1,kl
       do iq=1,ifull
        p(iq,k)=p(iq,k)+rdry*tv(iq,k)*psl(iq)
       enddo     ! iq loop
      enddo      ! k  loop

!     calculate "basic-linear" geopotential height and put in tempry
      do iq=1,ifull            ! tempry used by pextras
       tempry(iq,1)=zs(iq)+bet(1)*t(iq,1)+rdry*tbar2d(iq)*psl(iq)
      enddo      ! iq loop
      do k=2,kl
       do iq=1,ifull
        tempry(iq,k)=tempry(iq,k-1)
     .                 +bet(k)*t(iq,k)+betm(k)*t(iq,k-1)
       enddo     ! iq loop
      enddo      ! k  loop
!     save "full-linear" - "basic-linear" geopotential height in pextras
      pextras(:,:)=p(1:ifull,:)-tempry(:,:)
      if(npgf>0)then
!       N.B. these pextras terms can bebundled in here with un, vn
!       but their horiz derivs need to be staggered too
        pexx(1:ifull,:)=pextras(:,:)
	 pextras(:,:)=0.
	 if(npgf==1)p(1:ifull,:)=tempry(:,:)
        call bounds(pexx)
      endif
      call bounds(p)

!     Now set up in ux,vx the tendencies from "full-linear" terms.
!     The following un, vn contributions include residual psl.grad(T) terms
!     Following now always basically morder=24 2nd order, Corby style scheme
      do k=1,kl
!cdir nodep
         do iq=1,ifull          ! calculate staggered contributions first
            aa(iq,k)=-emu(iq)*(p(ie(iq),k)-p(iq,k))*
     .                                      (1.-epsu)*.5*dt/ds    ! 2nd order
            bb(iq,k)=-emv(iq)*(p(in(iq),k)-p(iq,k) )*
     .                                      (1.-epsu)*.5*dt/ds    ! 2nd order
         enddo                  ! iq loop
!      morder=24 scheme with residual terms for un, vn
!cdir nodep
         do iq=1,ifull  ! calculate staggered contributions first
            cc(iq,k)=emu(iq)*(psl(ie(iq))+psl(iq))*
     .                          (tv(ie(iq),k)-tv(iq,k))*.5*rdry/ds
            dd(iq,k)=emv(iq)*(psl(in(iq))+psl(iq))*
     .                          (tv(in(iq),k)-tv(iq,k))*.5*rdry/ds
         enddo                  ! iq loop
	  if(npgf==1)then
!cdir nodep
           do iq=1,ifull  ! in following pexx is really pextras
            cc(iq,k)=cc(iq,k)-emu(iq)*(pexx(ie(iq),k)-pexx(iq,k))/ds
            dd(iq,k)=dd(iq,k)-emv(iq)*(pexx(in(iq),k)-pexx(iq,k))/ds
           enddo                  ! iq loop
	  endif  ! (npgf==1)
	  if(npgf==2)then
!cdir nodep
           do iq=1,ifull  ! in following pexx is really pextras
            cc(iq,k)=cc(iq,k)-.5*(1.+epsu)*emu(iq)*
     .                          (pexx(ie(iq),k)-pexx(iq,k))/ds
            dd(iq,k)=dd(iq,k)-.5*(1.+epsu)*emv(iq)*
     .                          (pexx(in(iq),k)-pexx(iq,k))/ds
           enddo                  ! iq loop
	  endif  ! (npgf==2)

!cdir nodep
         if(nphip>1.and.k<=1)then
            do iq=1,ifull       ! calculate nphip contributions 
               cc(iq,k)=-emu(iq)*(dphi_dx(iq,k)
     .                     -p(ie(iq),k)+p(iq,k))/ds
               dd(iq,k)=-emv(iq)*(dphi_dy(iq,k)
     .                     -p(in(iq),k)+p(iq,k))/ds
            enddo               ! iq loop
            if ( mydiag ) then
               iq=idjd
               print *,"phi's ",k,dphi_dx(iq,k),p(ie(iq),k)-p(iq,k),
     .                      dphi_dy(iq,k),p(in(iq),k)-p(iq,k)
            end if
         endif                  ! (nphip>1)
      end do ! k
      call unstaguv(aa,bb,aa2,bb2) ! convert to unstaggered positions
      call unstaguv(cc,dd,cc2,dd2)
      ux(1:ifull,:)=u(1:ifull,:)+aa2(1:ifull,:)
      vx(1:ifull,:)=v(1:ifull,:)+bb2(1:ifull,:)
      if(nonl<0)then
        ux(1:ifull,:)=ux(1:ifull,:)+dt*un(1:ifull,:) ! possible physics contrib
        vx(1:ifull,:)=vx(1:ifull,:)+dt*vn(1:ifull,:) ! possible physics contrib
        un(1:ifull,:)=0.
        vn(1:ifull,:)=0.
      endif
      un(1:ifull,:)=un(1:ifull,:)+cc2(1:ifull,:)
      vn(1:ifull,:)=vn(1:ifull,:)+dd2(1:ifull,:)
      if(diag)then
        if(mydiag) print *,'tv ',(tv(idjd,kk),kk=1,kl)
        call printa('aa  ',aa,ktau,nlv,ia,ib,ja,jb,0.,1.)
        call printa('aa2 ',aa2,ktau,nlv,ia,ib,ja,jb,0.,1.)
        call printa('bb  ',bb,ktau,nlv,ia,ib,ja,jb,0.,1.)
        call printa('bb2 ',bb2,ktau,nlv,ia,ib,ja,jb,0.,1.)
        call printa('cc  ',cc,ktau,nlv,ia,ib,ja,jb,0.,dt)
        call printa('cc2 ',cc2,ktau,nlv,ia,ib,ja,jb,0.,dt)
        call printa('dd  ',dd,ktau,nlv,ia,ib,ja,jb,0.,dt)
        call printa('dd2 ',dd2,ktau,nlv,ia,ib,ja,jb,0.,dt)
      endif                     ! (diag.and.k==nlv)

      if(nonl<0)then
        if(ktau>2)then
	   aa(:,:)=tn(:,:)+.5*(tnsav(:,:)-tnsavv(:,:))
	   bb(:,:)=un(:,:)+.5*(unsav(:,:)-unsavv(:,:))
	   cc(:,:)=vn(:,:)+.5*(vnsav(:,:)-vnsavv(:,:))
	   tnsavv(:,:)=tnsav(:,:)
	   tnsav(:,:)=tn(:,:)
	   tn(:,:)=aa(:,:)
	   unsavv(:,:)=unsav(:,:)
	   unsav(:,:)=un(:,:)
	   un(:,:)=bb(:,:)
	   vnsavv(:,:)=vnsav(:,:)
	   vnsav(:,:)=vn(:,:)
	   vn(:,:)=cc(:,:)
	 else
	   tnsavv(:,:)=tnsav(:,:)
	   tnsav(:,:)=tn(:,:)
	   unsavv(:,:)=unsav(:,:)
	   unsav(:,:)=un(:,:)
	   vnsavv(:,:)=vnsav(:,:)
	   vnsav(:,:)=vn(:,:)
	 endif
      endif
!     finish evaluation of tx,ux,vx by adding in part of nonlinear terms
      cnon=.5    ! cnon=.5 shares tn,un,vn between here & upglobal/adjust
      if(abs(nonl)==1.or.nonl==11)cnon=1.   ! cnon=1. all done here
      if(abs(nonl)==2)then
        tn(:,:)=2.*tn(:,:)   ! all tn added in adjust
        un(:,:)=2.*un(:,:)   ! all un added in adjust
        vn(:,:)=2.*vn(:,:)   ! all vn added in adjust
      elseif(nonl==3)then   ! Adams-Bashforth style
        if(ktau==1)then
          tnsav(:,:)=tn(:,:)
          unsav(:,:)=un(:,:)
          vnsav(:,:)=vn(:,:)
        endif
        t(1:ifull,:)  = t(1:ifull,:) + 1.5*dt*tn(1:ifull,:) -  !cy
     &                  0.5*dt*tnsav(1:ifull,:)                !cy
cy        tx(1:ifull,:) = tx(1:ifull,:) + 1.5*dt*tn(1:ifull,:) -  !cb
cy     &                  0.5*dt*tnsav(1:ifull,:)                 !cb
        ux(1:ifull,:) = ux(1:ifull,:) + 1.5*dt*un(1:ifull,:) -
     &                  0.5*dt*unsav(1:ifull,:)
        vx(1:ifull,:) = vx(1:ifull,:) + 1.5*dt*vn(1:ifull,:) -
     &                  0.5*dt*vnsav(1:ifull,:)
      else       ! (nonl=0, 1, or 11)  usual code has nonl=0
        t(1:ifull,:)  = t(1:ifull,:)  + cnon*dt*tn(1:ifull,:) !cy
cy        tx(1:ifull,:) = tx(1:ifull,:) + cnon*dt*tn(1:ifull,:)  !cb
        ux(1:ifull,:) = ux(1:ifull,:) + cnon*dt*un(1:ifull,:)
        vx(1:ifull,:) = vx(1:ifull,:) + cnon*dt*vn(1:ifull,:)
      endif

      if(diag)then
         if ( mydiag ) then
            print *,'at end of nonlin; nvad,idjd = ', nvad,idjd
            print *,'p1 . & e ',p(idjd,nlv),p(ie(idjd),nlv)
            print *,'p1 . & n ',p(idjd,nlv),p(in(idjd),nlv)
            print *,'psl . & e ',psl(idjd),psl(ie(idjd))
            print *,'psl . & n ',psl(idjd),psl(in(idjd))
         end if
         call printa('psl ',psl,ktau,0,ia,ib,ja,jb,0.,100.)
!!!     Requires extra bounds call
!!!        print  *,'pslx -1,0,1 ',(pslx(idjd+n,nlv),n=-1,1)
         call printa('pslx',pslx,ktau,nlv,ia,ib,ja,jb,0.,100.)
         if ( mydiag ) print *,'tnsav*dt b ',(tnsav(idjd,k)*dt,k=1,kl)
         if ( mydiag ) print *,'tn*dt b ',(tn(idjd,k)*dt,k=1,kl)
         call printa('tn  ',tn,ktau,nlv,ia,ib,ja,jb,0.,100.*dt)
         if ( mydiag ) print *,'unsav*dt b ',(unsav(idjd,k)*dt,k=1,kl)
         if ( mydiag ) print *,'un*dt b ',(un(idjd,k)*dt,k=1,kl)
         call printa('un  ',un,ktau,nlv,ia,ib,ja,jb,0.,100.*dt)
         if ( mydiag ) print *,'vnsav*dt b ',(vnsav(idjd,k)*dt,k=1,kl)
         if ( mydiag ) print *,'vn*dt b ',(vn(idjd,k)*dt,k=1,kl)
         call printa('vn  ',vn,ktau,nlv,ia,ib,ja,jb,0.,100.*dt)
         if ( mydiag ) print *,'tx ',(tx(idjd,k),k=1,kl)
         call printa('tx  ',tx,ktau,nlv,ia,ib,ja,jb,200.,1.)
         if ( mydiag ) print *,'ux ',(ux(idjd,k),k=1,kl)
         call printa('ux  ',ux,ktau,nlv,ia,ib,ja,jb,0.,1.)
         if ( mydiag ) print *,'vx ',(vx(idjd,k),k=1,kl)
         call printa('vx  ',vx,ktau,nlv,ia,ib,ja,jb,0.,1.)
      endif

      if(abs(nonl)==1.or.nonl==11)then   ! option A
        un(:,:)=0. 
        vn(:,:)=0.
        tn(:,:)=0.
      endif
      call end_log(nonlin_end)
      return
      end
