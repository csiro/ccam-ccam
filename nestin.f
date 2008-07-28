      subroutine nestin               
      use cc_mpi, only : myid, mydiag
      use diag_m
      include 'newmpar.h'
!     ik,jk,kk are array dimensions read in infile - not for globpea
!     int2d code - not used for globpea
      include 'aalat.h'
      include 'arrays.h'
      include 'const_phys.h'
      include 'dates.h'    ! mtimer
      include 'dava.h'
      include 'davb.h'     ! psls,qgg,tt,uu,vv
      include 'indices.h'
      include 'latlong.h'
      include 'map.h'
      include 'parm.h'     ! qgmin
      include 'pbl.h'      ! tss
      include 'sigs.h'
      include 'soil.h'     ! sicedep fracice
      include 'soilsnow.h' ! tgg
      include 'stime.h'    ! kdate_s,ktime_s  sought values for data read
      common/nest/ta(ifull,kl),ua(ifull,kl),va(ifull,kl),psla(ifull),
     .            tb(ifull,kl),ub(ifull,kl),vb(ifull,kl),pslb(ifull),
     .            qa(ifull,kl),qb(ifull,kl),tssa(ifull),tssb(ifull),
     .            sicedepb(ifull),fraciceb(ifull)
      common/schmidtx/rlong0x,rlat0x,schmidtx ! infile, newin, nestin, indata
      real sigin
      integer ik,jk,kk
      common/sigin/ik,jk,kk,sigin(kl)  ! for vertint, infile
      real zsb(ifull)
      integer, dimension(ifull) :: isoilm_h ! MJT lsmask
      integer num,mtimea,mtimeb
      data num/0/,mtimea/0/,mtimeb/-1/
      save num,mtimea,mtimeb
!     mtimer, mtimeb are in minutes
      if(ktau<100.and.myid==0)then
        print *,'in nestin ktau,mtimer,mtimea,mtimeb ',
     &                     ktau,mtimer,mtimea,mtimeb
        print *,'with kdate_s,ktime_s >= ',kdate_s,ktime_s
      end if
      if(mtimeb==-1)then
        if ( myid==0 )
     &  print *,'set nesting fields to those already read in via indata'
        do iq=1,ifull
         pslb(iq)=psl(iq)
         tssb(iq)=tss(iq)
         sicedepb(iq)=sicedep(iq)  ! maybe not needed
         fraciceb(iq)=fracice(iq)
        enddo
        tb(1:ifull,:)=t(1:ifull,:)
        qb(1:ifull,:)=qg(1:ifull,:)
        ub(1:ifull,:)=u(1:ifull,:)
        vb(1:ifull,:)=v(1:ifull,:)
        mtimeb=-2
        return
      endif       ! (mtimeb==-1)

      if(mtimer<=mtimeb)go to 6  ! allows for dt<1 minute
      if(mtimeb==-2)mtimeb=mtimer  
!     transfer mtimeb fields to mtimea and update sice variables
      mtimea=mtimeb
      psla(:)=pslb(:)
      tssa(:)=tssb(:)
      ta(1:ifull,:)=tb(1:ifull,:)
      qa(1:ifull,:)=qb(1:ifull,:)
      ua(1:ifull,:)=ub(1:ifull,:)
      va(1:ifull,:)=vb(1:ifull,:)
!     following sice updating code moved from sflux Jan '06      
!     check whether present ice points should change to/from sice points
      do iq=1,ifull
       if(fraciceb(iq)>0.)then
!        N.B. if already a sice point, keep present tice (in tgg3)
         if(fracice(iq)==0.)then
           tgg(iq,3)=min(271.2,tssb(iq),tb(iq,1)+.04*6.5) ! for 40 m lev1
         endif  ! (fracice(iq)==0.)
!        set averaged tss (tgg1 setting already done)
         tss(iq)=tgg(iq,3)*fraciceb(iq)+tssb(iq)*(1.-fraciceb(iq))
       endif  ! (fraciceb(iq)==0.)
      enddo	! iq loop
      sicedep(:)=sicedepb(:)  ! from Jan 06
      fracice(:)=fraciceb(:)
!     because of new zs etc, ensure that sice is only over sea
      do iq=1,ifull
       if(fracice(iq)<.02)fracice(iq)=0.
       if(land(iq))then
         sicedep(iq)=0.
         fracice(iq)=0.
       else
         if(fracice(iq)>0..and.sicedep(iq)==0.)then
!          assign to 2. in NH and 1. in SH (according to spo)
!          do this in indata, amipdata and nestin because of onthefly
           if(rlatt(iq)>0.)then
             sicedep(iq)=2.
           else
             sicedep(iq)=1.
           endif ! (rlatt(iq)>0.)
         elseif(fracice(iq)==0..and.sicedep(iq)>0.)then  ! e.g. from Mk3  
           fracice(iq)=1.
         endif  ! (fracice(iq)>0..and.sicedep(iq)==0.) .. elseif ..
       endif    ! (land(iq))
      enddo     ! iq loop

!     read tb etc  - for globpea, straight into tb etc
      if(io_in==1)then
        call infil(1,kdate_r,ktime_r,timeg_b,ds_r, 
     .              pslb,zsb,tssb,sicedepb,fraciceb,tb,ub,vb,qb,
     .              isoilm_h) ! MJT lsmask
      endif   ! (io_in==1)

      if(io_in==-1)then
         call onthefl(1,kdate_r,ktime_r,
     &                 pslb,zsb,tssb,sicedepb,fraciceb,tb,ub,vb,qb) 
      endif   ! (io_in==1)
      tssb(:) = abs(tssb(:))  ! moved here Mar '03
      if (mydiag) then
        write (6,"('zsb# nestin  ',9f7.1)") diagvals(zsb)
        write (6,"('tssb# nestin ',9f7.1)") diagvals(tssb) 
      end if
   
      if(abs(rlong0  -rlong0x)>.01.or.
     &   abs(rlat0    -rlat0x)>.01.or.
     &   abs(schmidt-schmidtx)>.01)stop "grid mismatch in infile"

!     kdhour=(ktime_r-ktime)/100     ! integer hour diff
      kdhour=ktime_r/100-ktime/100   ! integer hour diff from Oct '05
      kdmin=(ktime_r-100*(ktime_r/100))-(ktime-100*(ktime/100))
      if ( myid == 0 ) then
        print *,'nesting file has: kdate_r,ktime_r,kdhour,kdmin ',
     &                             kdate_r,ktime_r,kdhour,kdmin
      end if
      mtimeb=60*24*(iabsdate(kdate_r,kdate)-iabsdate(kdate,kdate))
     .               +60*kdhour+kdmin
      if ( myid == 0 ) then
        print *,'kdate_r,iabsdate ',kdate_r,iabsdate(kdate_r,kdate)
        print *,'giving mtimeb = ',mtimeb
!     print additional information
        print *,' kdate ',kdate,' ktime ',ktime
        print *,'timeg,mtimer,mtimea,mtimeb: ',
     &           timeg,mtimer,mtimea,mtimeb
        print *,'ds,ds_r ',ds,ds_r
      end if

!     ensure qb big enough, but not too big in top levels (from Sept '04)
      qb(1:ifull,1:kk)=max(qb(1:ifull,1:kk),qgmin)
      do k=kk-2,kk
       qb(1:ifull,k)=min(qb(1:ifull,k),10.*qgmin)
      enddo

      if(mod(ktau,nmaxpr)==0.or.ktau==2.or.diag)then
!       following is useful if troublesome data is read in
        if ( myid == 0 ) then
          print *,'following max/min values printed from nestin'
        end if
        call maxmin(ub,'ub',ktau,1.,kk)
        call maxmin(vb,'vb',ktau,1.,kk)
        call maxmin(tb,'tb',ktau,1.,kk)
        call maxmin(qb,'qb',ktau,1.e3,kk)
        if ( myid == 0 ) then
          print *,'following are really psl not ps'
        end if
        call maxmin(pslb,'ps',ktau,100.,1)
      endif

!     if(kk<kl)then
      if(abs(sig(2)-sigin(2))>.0001)then   ! 11/03
!       this section allows for different number of vertical levels
!       presently assume sigin (up to kk) levels are set up as per nsig=6
!       option in eigenv, though original csiro9 levels are sufficiently
!       close for these interpolation purposes.
        if(ktau==1.and.mydiag)then
          print*,'calling vertint with kk,sigin ',kk,sigin(1:kk)
        endif
        if(diag.and.mydiag)then
          print *,'kk,sigin ',kk,(sigin(k),k=1,kk)
          print *,'tb before vertint ',(tb(idjd,k),k=1,kk)
        endif
        call vertint(tb,1)  ! transforms tb from kk to kl
        if(diag.and.mydiag)then
          print *,'tb after vertint ',(tb(idjd,k),k=1,kk)
          print *,'qb before vertint ',(qb(idjd,k),k=1,kk)
        endif
        call vertint(qb,2)
        if(diag.and.mydiag)print *,'qb after vertint ',qb(idjd,1:kk)
        call vertint(ub,3)
        call vertint(vb,4)
      endif  ! (abs(sig(2)-sigin(2))>.0001)

!     N.B. tssb (sea) only altered for newtop=2 (done here now)
      if(newtop==2)then
!       reduce sea tss to mslp      e.g. for QCCA in NCEP GCM
        do iq=1,ifull
         if(tssb(iq)<0.)tssb(iq)=
     .                       tssb(iq)-zsb(iq)*stdlapse/grav  ! N.B. -
        enddo
      endif  ! (newtop==2)

      if(newtop>=1)then
!       in these cases redefine pslb, tb and (effectively) zsb using zs
!       this keeps fine-mesh land mask & zs
!       presently simplest to do whole pslb, tb (& qb) arrays
        if(nmaxpr==1.and.mydiag)then
          print *,'zs (idjd) :',zs(idjd)
          print *,'zsb (idjd) :',zsb(idjd)
          write (6,"('100*psl.wesn ',2p5f8.3)") psl(idjd),psl(iw(idjd)),
     &              psl(ie(idjd)),psl(is(idjd)),psl(in(idjd))
          write (6,"('ps.wesn ',-2p5f9.3)") ps(idjd),
     &           ps(iw(idjd)),ps(ie(idjd)),ps(is(idjd)),ps(in(idjd))
          print *,'pslb in(idjd) :',pslb(idjd)
          print *,'now call retopo from nestin'
        endif
        call retopo(pslb,zsb,zs,tb,qb)
        if(nmaxpr==1.and.mydiag)then
          write (6,"('100*pslb.wesn ',2p5f8.3)") pslb(idjd),
     &       pslb(iw(idjd)),pslb(ie(idjd)),pslb(is(idjd)),pslb(in(idjd))
          print *,'pslb out(idjd) :',pslb(idjd)
          print *,'after pslb print; num= ',num
        endif
      endif   !  newtop>=1

      if(num==0)then
        num=1
        call printa('zs  ',zs        ,ktau,0  ,ia,ib,ja,jb,0.,.01)
        call printa('zsb ',zsb       ,ktau,0  ,ia,ib,ja,jb,0.,.01)
        call printa('psl ',psl       ,ktau,0  ,ia,ib,ja,jb,0.,1.e2)
        call printa('pslb',pslb      ,ktau,0  ,ia,ib,ja,jb,0.,1.e2)
        call printa('t   ',t,ktau,nlv,ia,ib,ja,jb,200.,1.)
        call printa('tb  ',tb,ktau,nlv,ia,ib,ja,jb,200.,1.)
        call printa('u   ',u,ktau,nlv,ia,ib,ja,jb,0.,1.)
        call printa('ub  ',ub,ktau,nlv,ia,ib,ja,jb,0.,1.)
        call printa('v   ',v,ktau,nlv,ia,ib,ja,jb,0.,1.)
        call printa('vb  ',vb,ktau,nlv,ia,ib,ja,jb,0.,1.)
        call printa('davt',davt,0,0,ia,ib,ja,jb,0.,10.)
        return
      endif   !  num==0

!     now use tt, uu, vv arrays for time interpolated values
6     timerm=ktau*dt/60.   ! real value in minutes (in case dt < 60 seconds)
      cona=(mtimeb-timerm)/real(mtimeb-mtimea)
      conb=(timerm-mtimea)/real(mtimeb-mtimea)
      psls(:)=cona*psla(:)+conb*pslb(:)
      tt (:,:)=cona*ta(:,:)+conb*tb(:,:)
      qgg(:,:)=cona*qa(:,:)+conb*qb(:,:)
      uu (:,:)=cona*ua(:,:)+conb*ub(:,:)
      vv (:,:)=cona*va(:,:)+conb*vb(:,:)

!     calculate time interpolated tss 
      if(namip.ne.0)return     ! namip SSTs/sea-ice take precedence
      do iq=1,ifull
       if(.not.land(iq))then
         tss(iq)=cona*tssa(iq)+conb*tssb(iq)
         tgg(iq,1)=tss(iq)
       endif  ! (.not.land(iq))
      enddo   ! iq loop 
      return
      end

      subroutine nestinb  ! called for mbd>0 - spectral filter method ! MJT CHANGE - delete mins_mbd
!     this is x-y-z version      
      use cc_mpi, only : myid, mydiag
      use diag_m
      implicit none
      integer, parameter :: ntest=0 
      include 'newmpar.h'
      include 'aalat.h'
      include 'arrays.h'
      include 'const_phys.h'
      include 'darcdf.h' ! for ncid
      include 'netcdf.inc'
      include 'mpif.h'
      include 'dates.h'    ! mtimer
      include 'indices.h'
      include 'latlong.h'
      include 'liqwpar.h'  ! ifullw,qfg,qlg
      include 'parm.h'     ! qgmin
      include 'pbl.h'      ! tss
      include 'sigs.h'
      include 'soil.h'     ! sicedep fracice
      include 'soilsnow.h' ! tgg
      include 'stime.h'    ! kdate_s,ktime_s  sought values for data read
      common/nest/ta(ifull,kl),ua(ifull,kl),va(ifull,kl),psla(ifull),
     .            tb(ifull,kl),ub(ifull,kl),vb(ifull,kl),pslb(ifull),
     .            qa(ifull,kl),qb(ifull,kl),tssa(ifull),tssb(ifull),
     .            sicedepb(ifull),fraciceb(ifull)
      common/schmidtx/rlong0x,rlat0x,schmidtx ! infile, newin, nestin, indata
      real sigin
      integer ik,jk,kk
      common/sigin/ik,jk,kk,sigin(kl)  ! for vertint, infile
      integer mtimeb,kdate_r,ktime_r
      integer ::  iabsdate,iq,k,kdhour,kdmin
      integer, dimension(ifull) :: isoilm_h ! MJT lsmask
      integer, save :: ncount = -1 ! MJT daily ave
      real :: ds_r,rlong0x,rlat0x
      real :: schmidtx,timeg_b
      real :: psla,pslb,qa,qb,ta,tb,tssa,tssb,ua,ub,va,vb
      real :: fraciceb,sicedepb
      real, dimension(ifull) ::  zsb
      real, parameter :: alpr = 0.5 ! MJT daily ave
      data mtimeb/-1/ 
      save mtimeb 
  
      if ((nud_uv.eq.8).and.(ncount.gt.0)) then ! MJT daily ave
        ! update average
        psla(:)=psla(:)+psl(1:ifull)
        ua(:,:)=ua(:,:)+u(1:ifull,:)
        va(:,:)=va(:,:)+v(1:ifull,:)
        ta(:,:)=ta(:,:)+t(1:ifull,:)
        qa(:,:)=qa(:,:)+qg(1:ifull,:)
        ncount=ncount+1
      end if 
      
      if ((mtimer<mtimeb).and.(ktau.gt.0)) return
 
      !------------------------------------------------------------------------------
!     mtimer, mtimeb are in minutes
      if(ktau<100.and.myid==0)then
        print *,'in nestinb ktau,mtimer,mtimeb,io_in ',
     &                      ktau,mtimer,mtimeb,io_in ! MJT CHANGE - delete mtimea
        print *,'with kdate_s,ktime_s >= ',kdate_s,ktime_s
      end if

      if ((mtimer>mtimeb).or.(ktau.le.0)) then
      
        if (nud_uv.eq.8) then ! MJT daily ave
          ! reset average
          psla(:)=psl(1:ifull)
          ua(:,:)=u(1:ifull,:)
          va(:,:)=v(1:ifull,:)
          ta(:,:)=t(1:ifull,:)
          qa(:,:)=qg(1:ifull,:)
          ncount=1
        end if

!      read tb etc  - for globpea, straight into tb etc
       if(io_in==1)then
         call infil(1,kdate_r,ktime_r,timeg_b,ds_r, 
     .               pslb,zsb,tssb,sicedepb,fraciceb,tb,ub,vb,qb,
     .               isoilm_h) ! MJT lsmask
       endif   ! (io_in==1)

       if(io_in==-1)then
          call onthefl(1,kdate_r,ktime_r,
     &                  pslb,zsb,tssb,sicedepb,fraciceb,tb,ub,vb,qb) 
       endif   ! (io_in==1)
       tssb(:) = abs(tssb(:))  ! moved here Mar '03
       if (mydiag) then
         write (6,"('zsb# nestinb  ',9f7.1)") diagvals(zsb)
         write (6,"('tssb# nestinb ',9f7.1)") diagvals(tssb) 
       end if
   
       if(abs(rlong0  -rlong0x)>.01.or.
     &    abs(rlat0    -rlat0x)>.01.or.
     &    abs(schmidt-schmidtx)>.01)stop "grid mismatch in infile"

       kdhour=ktime_r/100-ktime/100   ! integer hour diff from Oct '05
       kdmin=(ktime_r-100*(ktime_r/100))-(ktime-100*(ktime/100))
       if ( myid == 0 ) then
         print *,'nestinb file has: kdate_r,ktime_r,kdhour,kdmin ',
     &                              kdate_r,ktime_r,kdhour,kdmin
       end if
       mtimeb=60*24*(iabsdate(kdate_r,kdate)-iabsdate(kdate,kdate))
     .                +60*kdhour+kdmin
       if ( myid == 0 ) then
         print *,'kdate_r,iabsdate ',kdate_r,iabsdate(kdate_r,kdate)
         print *,'giving mtimeb = ',mtimeb
!        print additional information
         print *,' kdate ',kdate,' ktime ',ktime
         print *,'timeg,mtimer,mtimeb: ',
     &            timeg,mtimer,mtimeb ! MJT CHANGE - delete mtimea
         print *,'ds,ds_r ',ds,ds_r
       end if

!      ensure qb big enough, but not too big in top levels (from Sept '04)
       qb(1:ifull,1:kk)=max(qb(1:ifull,1:kk),qgmin)
       do k=kk-2,kk
        qb(1:ifull,k)=min(qb(1:ifull,k),10.*qgmin)
       enddo

       if(mod(ktau,nmaxpr)==0.or.ktau==2.or.diag)then
!        following is useful if troublesome data is read in
         if ( myid == 0 ) then
           print *,'following max/min values printed from nestin'
         end if
         call maxmin(ub,'ub',ktau,1.,kk)
         call maxmin(vb,'vb',ktau,1.,kk)
         call maxmin(tb,'tb',ktau,1.,kk)
         call maxmin(qb,'qb',ktau,1.e3,kk)
       endif
       if ( myid == 0 ) then
         print *,'following in nestinb after read pslb are psl not ps'
       end if
       call maxmin(pslb,'pB',ktau,100.,1)

!      if(kk<kl)then
       if(abs(sig(2)-sigin(2))>.0001)then   ! 11/03
!        this section allows for different number of vertical levels
!        presently assume sigin (up to kk) levels are set up as per nsig=6
!        option in eigenv, though original csiro9 levels are sufficiently
!        close for these interpolation purposes.
         if(ktau==1.and.mydiag)then
           print*,'calling vertint with kk,sigin ',kk,sigin(1:kk)
         endif
         if(diag.and.mydiag)then
           print *,'kk,sigin ',kk,(sigin(k),k=1,kk)
           print *,'tb before vertint ',(tb(idjd,k),k=1,kk)
         endif
         call vertint(tb,1)  ! transforms tb from kk to kl
         if(diag.and.mydiag)then
           print *,'tb after vertint ',(tb(idjd,k),k=1,kk)
           print *,'qb before vertint ',(qb(idjd,k),k=1,kk)
         endif
         call vertint(qb,2)
         if(diag.and.mydiag)print *,'qb after vertint ',qb(idjd,1:kk)
         call vertint(ub,3)
         call vertint(vb,4)
       endif  ! (abs(sig(2)-sigin(2))>.0001)

!      N.B. tssb (sea) only altered for newtop=2 (done here now)
       if(newtop==2)then
!        reduce sea tss to mslp      e.g. for QCCA in NCEP GCM
         do iq=1,ifull
          if(tssb(iq)<0.)tssb(iq)=
     .                        tssb(iq)-zsb(iq)*stdlapse/grav  ! N.B. -
         enddo
       endif  ! (newtop==2)

       if(newtop>=1)then
!        in these cases redefine pslb, tb and (effectively) zsb using zs
!        this keeps fine-mesh land mask & zs
!        presently simplest to do whole pslb, tb (& qb) arrays
         if(mydiag)then
           print *,'zs (idjd) :',zs(idjd)
           print *,'zsb (idjd) :',zsb(idjd)
           write (6,"('100*psl.wesn ',2p5f8.3)") psl(idjd),
     &           psl(iw(idjd)),psl(ie(idjd)),psl(is(idjd)),psl(in(idjd))
           write (6,"('ps.wesn ',-2p5f9.3)") ps(idjd),
     &           ps(iw(idjd)),ps(ie(idjd)),ps(is(idjd)),ps(in(idjd))
           print *,'pslb in(idjd) :',pslb(idjd)
           print *,'call retopo from nestin; psl# prints refer to pslb'
         endif
         call retopo(pslb,zsb,zs,tb,qb)
         if(mydiag)then
           write (6,"('100*pslb.wesn ',2p5f8.3)") pslb(idjd),
     &       pslb(iw(idjd)),pslb(ie(idjd)),pslb(is(idjd)),pslb(in(idjd))
         endif
       endif   !  newtop>=1

       return
      end if 

       if (nud_uv.eq.8) then
         ! preturb daily average
         if (myid == 0) print *,"Using averaged data for filter"
         psla(:)=alpr*(pslb(:)-psla(:)/real(ncount))
         ua(:,:)=alpr*(ub(:,:)-ua(:,:)/real(ncount))
         va(:,:)=alpr*(vb(:,:)-va(:,:)/real(ncount))
         ta(:,:)=alpr*(tb(:,:)-ta(:,:)/real(ncount))
         qa(:,:)=alpr*(qb(:,:)-qa(:,:)/real(ncount))
	 ncount=-1
       else
         ! preturb instantaneous
         psla(:)=pslb(:)-psl(1:ifull)
         ua(:,:)=ub(:,:)-u(1:ifull,:)
         va(:,:)=vb(:,:)-v(1:ifull,:)
         ta(:,:)=tb(:,:)-t(1:ifull,:)
         qa(:,:)=qb(:,:)-qg(1:ifull,:)
       end if

       call getspecdata(psla,ua,va,ta,qa)
       if ( myid == 0 ) then
        print *,'following after getspecdata are really psl not ps'
       end if
       call maxmin(pslb,'pB',ktau,100.,1)

!     following sice updating code moved from sflux Jan '06      
!     check whether present ice points should change to/from sice points
      do iq=1,ifull
       if(fraciceb(iq)>0.)then
!        N.B. if already a sice point, keep present tice (in tgg3)
         if(fracice(iq)==0.)then
           tgg(iq,3)=min(271.2,tssb(iq),tb(iq,1)+.04*6.5) ! for 40 m lev1
         endif  ! (fracice(iq)==0.)
!        set averaged tss (tgg1 setting already done)
         tss(iq)=tgg(iq,3)*fraciceb(iq)+tssb(iq)*(1.-fraciceb(iq))
       endif  ! (fraciceb(iq)==0.)
      enddo	! iq loop
      sicedep(:)=sicedepb(:)  ! from Jan 06
      fracice(:)=fraciceb(:)
!     because of new zs etc, ensure that sice is only over sea
      do iq=1,ifull
       if(fracice(iq)<.02)fracice(iq)=0.
       if(land(iq))then
         sicedep(iq)=0.
         fracice(iq)=0.
       else
         if(fracice(iq)>0..and.sicedep(iq)==0.)then
!          assign to 2. in NH and 1. in SH (according to spo)
!          do this in indata, amipdata and nestin because of onthefly
           if(rlatt(iq)>0.)then
             sicedep(iq)=2.
           else
             sicedep(iq)=1.
           endif ! (rlatt(iq)>0.)
         elseif(fracice(iq)==0..and.sicedep(iq)>0.)then  ! e.g. from Mk3  
           fracice(iq)=1.
         endif  ! (fracice(iq)>0..and.sicedep(iq)==0.) .. elseif ..
       endif    ! (land(iq))
      enddo     ! iq loop

!      calculate time interpolated tss 
       if(namip.ne.0.or.ntest.ne.0)return  ! namip SSTs/sea-ice take precedence
       do iq=1,ifull
        if(.not.land(iq))then
          tss(iq)=tssb(iq)
          tgg(iq,1)=tss(iq)
        endif  ! (.not.land(iq))
       enddo   ! iq loop 
      
      return
      end


      ! This subroutine gathers data for the MPI version of spectral downscaling
      subroutine getspecdata(pslb,ub,vb,tb,qb)

      use cc_mpi
      
      implicit none

      include 'newmpar.h'    ! ifull_g,kl
      include 'arrays.h'     ! u,v,t,qg,psl
      include 'const_phys.h'
      include 'parm.h'       ! mbd,schmidt,nud_uv,nud_p,nud_t,nud_q,kbotdav
      include 'xyzinfo.h'
      include 'vecsuv.h'
      include 'vecsuv_g.h'   ! ax_g,bx_g,ay_g,by_g,az_g,bz_g
      include 'savuvt.h'     ! savu,savv

      integer iq,k      
      real, dimension(ifull), intent(in) :: pslb
      real, dimension(ifull,kl), intent(in) :: ub,vb,tb,qb
      real, dimension(ifull) :: costh,sinth
      real, dimension(ifull,kl) :: delta
      real, dimension(ifull_g) :: pslc
      real, dimension(ifull_g,kl) :: uc,vc,wc,tc,qc
      real, dimension(ifull_g,kl) :: x_g,xx_g
      real savs1(ifull,2:kl),savu1(ifull,kl),savv1(ifull,kl)
      real den,polenx,poleny,polenz,zonx,zony,zonz
      common/savuv1/savs1,savu1,savv1

      if(nud_uv==3)then
        polenx=-cos(rlat0*pi/180.)
        poleny=0.
        polenz=sin(rlat0*pi/180.)
        do iq=1,ifull
         zonx=            -polenz*y(iq)
         zony=polenz*x(iq)-polenx*z(iq)
         zonz=polenx*y(iq)
         den=sqrt( max(zonx**2 + zony**2 + zonz**2,1.e-7) ) 
         costh(iq)= (zonx*ax(iq)+zony*ay(iq)+zonz*az(iq))/den
         sinth(iq)=-(zonx*bx(iq)+zony*by(iq)+zonz*bz(iq))/den
        enddo
      endif

      if (myid == 0) then     
        print *,"Gather data for spectral downscale"
        if(nud_p>0)call ccmpi_gather(pslb(:), pslc(:))
        if(nud_uv==3)then
          do k=kbotdav,kl
            delta(:,k)=costh(:)*ub(:,k)  ! uzon
     &                -sinth(:)*vb(:,k)
          end do
          call ccmpi_gather(delta(:,:), wc(:,:))
        elseif(nud_uv.ne.0)then
          call ccmpi_gather(ub(:,:), x_g(:,:))
          call ccmpi_gather(vb(:,:), xx_g(:,:))
          do k=kbotdav,kl
            uc(:,k)=ax_g(:)*x_g(:,k)+bx_g(:)*xx_g(:,k)
            vc(:,k)=ay_g(:)*x_g(:,k)+by_g(:)*xx_g(:,k)
            wc(:,k)=az_g(:)*x_g(:,k)+bz_g(:)*xx_g(:,k)
          end do
        endif
        if(nud_t>0)then
          call ccmpi_gather(tb(:,:), tc(:,:))
        end if
        if(nud_q>0)then
          call ccmpi_gather(qb(:,:), qc(:,:))
        end if
      else
        if(nud_p>0)call ccmpi_gather(pslb(:))
        if(nud_uv==3)then
          do k=kbotdav,kl
            delta(:,k)=costh(:)*ub(:,k)  ! uzon
     &                -sinth(:)*vb(:,k)
          end do
          call ccmpi_gather(delta(:,:))	
        elseif(nud_uv.ne.0)then
          call ccmpi_gather(ub(:,:))
          call ccmpi_gather(vb(:,:))
        endif
        if(nud_t>0)then
          call ccmpi_gather(tb(:,:))
        endif
        if(nud_q>0)then
          call ccmpi_gather(qb(:,:))
        endif
      end if

      !-----------------------------------------------------------------------
      if(nud_uv<0)then 
        if (myid == 0) then
          print *,"Fast spectral downscale"
          call fastspec((.1*real(mbd)/(pi*schmidt))**2
     &       ,pslc(:),uc(:,kbotdav:kl),vc(:,kbotdav:kl)
     &       ,wc(:,kbotdav:kl),tc(:,kbotdav:kl),qc(:,kbotdav:kl)) ! e.g. mbd=40
        end if
      elseif(nud_uv==9)then 
        if (myid == 0) print *,"Two dimensional spectral downscale"
        call slowspecmpi(myid,.1*real(mbd)/(pi*schmidt)
     &                ,pslc(:),uc(:,kbotdav:kl),vc(:,kbotdav:kl)
     &                ,wc(:,kbotdav:kl),tc(:,kbotdav:kl)
     &                ,qc(:,kbotdav:kl))
      elseif((mod(6,nproc)==0).or.(mod(nproc,6)==0))then
        if (myid == 0) print *,"Separable 1D downscale (MPI optimised)"
        call specfastmpi(myid,.1*real(mbd)/(pi*schmidt)
     &                ,pslc(:),uc(:,kbotdav:kl),vc(:,kbotdav:kl)
     &                ,wc(:,kbotdav:kl),tc(:,kbotdav:kl)
     &                ,qc(:,kbotdav:kl))
      else          !  usual choice e.g. for nud_uv=1 or 2
        if (myid == 0) print *,"Separable 1D downscale (MPI)"
        call fourspecmpi(myid,.1*real(mbd)/(pi*schmidt)
     &                ,pslc(:),uc(:,kbotdav:kl),vc(:,kbotdav:kl)
     &                ,wc(:,kbotdav:kl),tc(:,kbotdav:kl)
     &                ,qc(:,kbotdav:kl))
      endif  ! (nud_uv<0) .. else ..
      !-----------------------------------------------------------------------

      if (myid == 0) then
        print *,"Distribute data from spectral downscale"
        if (nud_p.gt.0) then
          call ccmpi_distribute(delta(:,1), pslc(:))
          psl(1:ifull)=psl(1:ifull)+delta(:,1)
        end if
        if(nud_uv==3)then
          call ccmpi_distribute(delta(:,:), wc(:,:))
          do k=kbotdav,kl 
            u(1:ifull,k)=u(1:ifull,k)+costh(:)*delta(:,k)
            v(1:ifull,k)=v(1:ifull,k)-sinth(:)*delta(:,k)
          end do
        elseif(nud_uv.ne.0) then
          do k=kbotdav,kl        
            x_g(:,k)=ax_g(:)*uc(:,k)+ay_g(:)*vc(:,k)+az_g(:)*wc(:,k)
            xx_g(:,k)=bx_g(:)*uc(:,k)+by_g(:)*vc(:,k)+bz_g(:)*wc(:,k)
          end do
          call ccmpi_distribute(delta(:,:), x_g(:,:))
          u(1:ifull,kbotdav:kl)=u(1:ifull,kbotdav:kl)
     &                         +delta(:,kbotdav:kl)
          savu(1:ifull,kbotdav:kl)=savu(1:ifull,kbotdav:kl)
     &                            +delta(:,kbotdav:kl)
          savu1(1:ifull,kbotdav:kl)=savu1(1:ifull,kbotdav:kl)
     &                            +delta(:,kbotdav:kl)
          call ccmpi_distribute(delta(:,:), xx_g(:,:))
          v(1:ifull,kbotdav:kl)=v(1:ifull,kbotdav:kl)
     &                         +delta(:,kbotdav:kl)
          savv(1:ifull,kbotdav:kl)=savv(1:ifull,kbotdav:kl)
     &                            +delta(:,kbotdav:kl)
          savv1(1:ifull,kbotdav:kl)=savv1(1:ifull,kbotdav:kl)
     &                            +delta(:,kbotdav:kl)
        end if
        if (nud_t.gt.0) then
          call ccmpi_distribute(delta(:,:), tc(:,:))
          t(1:ifull,kbotdav:kl)=t(1:ifull,kbotdav:kl)
     &                         +delta(:,kbotdav:kl)
        end if
        if (nud_q.gt.0) then
          call ccmpi_distribute(delta(:,:), qc(:,:))
          qg(1:ifull,kbotdav:kl)=max(qg(1:ifull,kbotdav:kl)
     &                          +delta(:,kbotdav:kl),qgmin)
        end if
      else
        if (nud_p.gt.0) then
          call ccmpi_distribute(delta(:,1))
          psl(1:ifull)=psl(1:ifull)+delta(:,1)
        end if
        if(nud_uv==3)then
          call ccmpi_distribute(delta(:,:))
          do k=kbotdav,kl
            u(1:ifull,k)=u(1:ifull,k)+costh(:)*delta(:,k)
            v(1:ifull,k)=v(1:ifull,k)-sinth(:)*delta(:,k)
          end do
        elseif (nud_uv.ne.0) then
          call ccmpi_distribute(delta(:,:))
          u(1:ifull,kbotdav:kl)=u(1:ifull,kbotdav:kl)
     &                         +delta(:,kbotdav:kl)
          savu(1:ifull,kbotdav:kl)=savu(1:ifull,kbotdav:kl)
     &                            +delta(:,kbotdav:kl)
          savu1(1:ifull,kbotdav:kl)=savu1(1:ifull,kbotdav:kl)
     &                            +delta(:,kbotdav:kl)
          call ccmpi_distribute(delta(:,:))
          v(1:ifull,kbotdav:kl)=v(1:ifull,kbotdav:kl)
     &                         +delta(:,kbotdav:kl)
          savu(1:ifull,kbotdav:kl)=savu(1:ifull,kbotdav:kl)
     &                            +delta(:,kbotdav:kl)
          savv1(1:ifull,kbotdav:kl)=savv1(1:ifull,kbotdav:kl)
     &                            +delta(:,kbotdav:kl)
        end if
        if (nud_t.gt.0) then
          call ccmpi_distribute(delta(:,:))
          t(1:ifull,kbotdav:kl)=t(1:ifull,kbotdav:kl)
     &                         +delta(:,kbotdav:kl)
        end if
        if (nud_q.gt.0) then
          call ccmpi_distribute(delta(:,:))
          qg(1:ifull,kbotdav:kl)=max(qg(1:ifull,kbotdav:kl)
     &                          +delta(:,kbotdav:kl),qgmin)
        end if
      end if
      
      if (nud_p>0) then
        ps(1:ifull)=1.e5*exp(psl(1:ifull)) ! Do not think this is needed, but kept it anyway - MJT
      end if
            
      return
      end subroutine getspecdata


      ! Fast spectral downscaling (JLM version)
      subroutine fastspec(cutoff2,psla,ua,va,wa,ta,qa)
      
      implicit none
      
      include 'newmpar.h'    ! ifull_g,kl
      include 'const_phys.h' ! rearth,pi,tpi
      include 'map_g.h'      ! em_g
      include 'indices_g.h'  ! in_g,ie_g,is_g,iw_g
      include 'parm.h'       ! ds,kbotdav
      include 'xyzinfo_g.h'    ! x_g,y_g,z_g

      integer, parameter :: ntest=0 
      integer i,j,k,n,n1,iq,iq1,num
      real, intent(in) :: cutoff2
      real, dimension(ifull_g), intent(inout) :: psla
      real, dimension(ifull_g,kbotdav:kl), intent(inout) :: ua,va,wa
      real, dimension(ifull_g,kbotdav:kl), intent(inout) :: ta,qa
      real, dimension(ifull_g) :: psls,sumwt
      real, dimension(ifull_g) :: psls2
      real, dimension(ifull_g), save :: xx,yy,zz
      real, dimension(ifull_g,kbotdav:kl) :: uu,vv,ww,tt,qgg
      real, dimension(ifull_g,kbotdav:kl) :: uu2,vv2,ww2,tt2,qgg2
      real emmin,dist,dist1,wt,wt1,xxmax,yymax,zzmax
      data num/1/
      save num
      
      ! myid must = 0 to get here.  So there is no need to check.
      
      if (num==1) then
      num=2
!       set up geometry for filtering through panel 1
!       x pass on panels 1, 2, 4, 5
!       y pass on panels 0, 1, 3, 4
!       z pass on panels 0, 2, 3, 5
        xx=0.
        yy=0.
        zz=0.
        do iq=1+il_g*il_g,3*il_g*il_g
          xx(iq)=xx(iw_g(iq))+sqrt((x_g(iq)-x_g(iw_g(iq)))**2+
     &           (y_g(iq)-y_g(iw_g(iq)))**2+(z_g(iq)-z_g(iw_g(iq)))**2)
        enddo
         do iq=1+4*il_g*il_g,6*il_g*il_g
          xx(iq)=xx(is_g(iq))+sqrt((x_g(iq)-x_g(is_g(iq)))**2+
     &           (y_g(iq)-y_g(is_g(iq)))**2+(z_g(iq)-z_g(is_g(iq)))**2)
        enddo
        do iq=1,2*il_g*il_g
          yy(iq)=yy(is_g(iq))+sqrt((x_g(iq)-x_g(is_g(iq)))**2+
     &           (y_g(iq)-y_g(is_g(iq)))**2+(z_g(iq)-z_g(is_g(iq)))**2)
        enddo
        do iq=1+3*il_g*il_g,5*il_g*il_g
          yy(iq)=yy(iw_g(iq))+sqrt((x_g(iq)-x_g(iw_g(iq)))**2+
     &           (y_g(iq)-y_g(iw_g(iq)))**2+(z_g(iq)-z_g(iw_g(iq)))**2)
        enddo
        if(mbd>0)then
         do iq=1,il_g*il_g
          zz(iq)=zz(iw_g(iq))+sqrt((x_g(iq)-x_g(iw_g(iq)))**2+
     &           (y_g(iq)-y_g(iw_g(iq)))**2+(z_g(iq)-z_g(iw_g(iq)))**2)
         enddo
         do iq=1+2*il_g*il_g,4*il_g*il_g
          zz(iq)=zz(is_g(iq))+sqrt((x_g(iq)-x_g(is_g(iq)))**2+
     &           (y_g(iq)-y_g(is_g(iq)))**2+(z_g(iq)-z_g(is_g(iq)))**2)
         enddo
         do iq=1+5*il_g*il_g,6*il_g*il_g
          zz(iq)=zz(iw_g(iq))+sqrt((x_g(iq)-x_g(iw_g(iq)))**2+
     &           (y_g(iq)-y_g(iw_g(iq)))**2+(z_g(iq)-z_g(iw_g(iq)))**2)
         enddo
        endif  ! (mbd>0)
        if(ntest>0)then
          do iq=1,144
           print *,'iq,xx,yy,zz ',iq,xx(iq),yy(iq),zz(iq)
          enddo
          do iq=il_g*il_g,il_g*il_g+il_g
           print *,'iq,xx,yy,zz ',iq,xx(iq),yy(iq),zz(iq)
          enddo
         do iq=4*il_g*il_g-il_g,4*il_g*il_g
           print *,'iq,xx,yy,zz ',iq,xx(iq),yy(iq),zz(iq)
          enddo
         do iq=5*il_g*il_g-il_g,5*il_g*il_g
           print *,'iq,xx,yy,zz ',iq,xx(iq),yy(iq),zz(iq)
          enddo
          print *,'xx mid:'   
          do i=1,48
           print *,'i xx',i,xx(il_g*il_g*1.5+i)
          enddo
          do i=1,48
           print *,'i xx',i+il_g,xx(il_g*il_g*2.5+i)
          enddo
          do i=1,48
           print *,'i xx',i+2*il_g,xx(il_g*il_g*4-il_g/2+i*il_g)
          enddo
          do i=1,48
           print *,'i xx',i+3*il_g,xx(il_g*il_g*5-il_g/2+i*il_g)
          enddo
          print *,'yy mid:'   
          do j=1,96
           print *,'j yy',j,yy(-il_g/2+j*il_g)
          enddo
          do j=1,48
           print *,'j yy',j+2*il_g,yy(il_g*il_g*3.5+j)
          enddo
          do j=1,48
           print *,'j yy',j+3*il_g,yy(il_g*il_g*4.5+j)
          enddo
!         wrap-around values defined by xx(il_g,5*il_g+j),j=1,il_g; yy(i,5*il_g),i=1,il_g
          print *,'wrap-round values'
          do j=1,il_g
           print *,'j,xx ',j,xx(6*il_g*il_g+1-j)       ! xx(il_g+1-j,il_g,5)
          enddo
          do i=1,il_g
           print *,'i,yy ',i,yy(5*il_g*il_g+il_g-il_g*i)   ! yy(il_g,il_g+1-i,4)
          enddo
          do j=1,il_g
           print *,'j,zz ',j,zz(5*il_g*il_g+il_g*j)      ! zz(il_g,j,5)
          enddo
        endif  ! ntest>0
      endif    !  num==1

      qgg(1:ifull_g,kbotdav:kl)=0.
      tt(1:ifull_g,kbotdav:kl)=0.
      uu(1:ifull_g,kbotdav:kl)=0.
      vv(1:ifull_g,kbotdav:kl)=0.
      ww(1:ifull_g,kbotdav:kl)=0.
      psls(1:ifull_g)=0.
      qgg2(1:ifull_g,kbotdav:kl)=0.
      tt2(1:ifull_g,kbotdav:kl)=0.
      uu2(1:ifull_g,kbotdav:kl)=0.
      vv2(1:ifull_g,kbotdav:kl)=0.
      ww2(1:ifull_g,kbotdav:kl)=0.
      psls2(1:ifull_g)=0.
      sumwt(1:ifull_g)=1.e-20   ! for undefined panels
      emmin=sqrt(cutoff2)*ds/rearth
      print *,'schmidt,cutoff,kbotdav ',schmidt,sqrt(cutoff2),kbotdav 
      
      do j=1,il_g                ! doing x-filter on panels 1,2,4,5
       xxmax=xx(il_g*(6*il_g-1)+il_g+1-j)
       print *,'j,xxmax ',j,xxmax
       do n=1,4*il_g
        if(n<=il_g)iq=il_g*(il_g+j-1)+n                   ! panel 1
        if(n>il_g.and.n<=2*il_g)iq=il_g*(2*il_g+j-2)+n      ! panel 2
        if(n>2*il_g)iq=il_g*(2*il_g+n-1)+il_g+1-j           ! panel 4,5
        
        if (em_g(iq).gt.emmin) then ! MJT
        
        do n1=n,4*il_g
!        following test shows on sx6 don't use "do n1=m+1,4*il_g"
!        if(n==4*il_g)print *,'problem for i,n,n1 ',i,n,n1
         if(n1<=il_g)iq1=il_g*(il_g+j-1)+n1               ! panel 1
         if(n1>il_g.and.n1<=2*il_g)iq1=il_g*(2*il_g+j-2)+n1 ! panel 2
         if(n1>2*il_g)iq1=il_g*(2*il_g+n1-1)+il_g+1-j       ! panel 4,5
         dist1=abs(xx(iq)-xx(iq1))
         dist=min(dist1,xxmax-dist1)
         wt=exp(-4.5*dist*dist*cutoff2)
         wt1=wt/em_g(iq1)
         wt=wt/em_g(iq)
         if(n==n1)wt1=0.  ! so as not to add in twice
c        if(iq==10345.or.iq1==10345)
c    &     print *,'iq,iq1,n,n1,xx,xx1,dist1,dist,wt,wt1 ',         
c    &              iq,iq1,n,n1,xx(iq),xx(iq1),dist1,dist,wt,wt1 
         sumwt(iq)=sumwt(iq)+wt1
         sumwt(iq1)=sumwt(iq1)+wt
!        producing "x-filtered" version of pslb-psl etc
c        psls(iq)=psls(iq)+wt1*(pslb(iq1)-psl(iq1))
c        psls(iq1)=psls(iq1)+wt*(pslb(iq)-psl(iq))
         psls(iq)=psls(iq)+wt1*psla(iq1)
         psls(iq1)=psls(iq1)+wt*psla(iq)
         do k=kbotdav,kl
          qgg(iq,k)=qgg(iq,k)+wt1*qa(iq1,k)
          qgg(iq1,k)=qgg(iq1,k)+wt*qa(iq,k)
          tt(iq,k)=tt(iq,k)+wt1*ta(iq1,k)
          tt(iq1,k)=tt(iq1,k)+wt*ta(iq,k)
          uu(iq,k)=uu(iq,k)+wt1*ua(iq1,k)
          uu(iq1,k)=uu(iq1,k)+wt*ua(iq,k)
          vv(iq,k)=vv(iq,k)+wt1*va(iq1,k)
          vv(iq1,k)=vv(iq1,k)+wt*va(iq,k)
          ww(iq,k)=ww(iq,k)+wt1*wa(iq1,k)
          ww(iq1,k)=ww(iq1,k)+wt*wa(iq,k)
         enddo  ! k loop
c        print *,'n,n1,dist,wt,wt1 ',n,n1,dist,wt,wt1
        enddo   ! n1 loop
        else
          sumwt(iq)=1.
        end if
       enddo    ! n loop
      enddo     ! j loop      
      if(nud_uv==-1)then
        do iq=1,ifull_g
         psls2(iq)=psls(iq)/sumwt(iq)
         do k=kbotdav,kl
          qgg2(iq,k)=qgg(iq,k)/sumwt(iq)
          tt2(iq,k)=tt(iq,k)/sumwt(iq)
          uu2(iq,k)=uu(iq,k)/sumwt(iq)
          vv2(iq,k)=vv(iq,k)/sumwt(iq)
          ww2(iq,k)=ww(iq,k)/sumwt(iq)
         enddo
        enddo
      else  ! original fast scheme
        do iq=1,ifull_g
         if(sumwt(iq).ne.1.e-20)then
           psla(iq)=psls(iq)/sumwt(iq)
           do k=kbotdav,kl
            qa(iq,k)=qgg(iq,k)/sumwt(iq)
            ta(iq,k)=tt(iq,k)/sumwt(iq)
            ua(iq,k)=uu(iq,k)/sumwt(iq)
            va(iq,k)=vv(iq,k)/sumwt(iq)
            wa(iq,k)=ww(iq,k)/sumwt(iq)
           enddo
         endif  ! (sumwt(iq).ne.1.e-20)
        enddo
      endif  ! (nud_uv==-1) .. else ..
      
      qgg(1:ifull_g,kbotdav:kl)=0.
      tt(1:ifull_g,kbotdav:kl)=0.
      uu(1:ifull_g,kbotdav:kl)=0.
      vv(1:ifull_g,kbotdav:kl)=0.
      ww(1:ifull_g,kbotdav:kl)=0.
      psls(1:ifull_g)=0.
      sumwt(1:ifull_g)=1.e-20   ! for undefined panels
      
      do i=1,il_g                ! doing y-filter on panels 0,1,3,4
       yymax=yy(il_g*(5*il_g-i+1))  
       do n=1,4*il_g
        if(n<=2*il_g)iq=il_g*(n-1)+i                      ! panel 0,1
        if(n>2*il_g.and.n<=3*il_g)iq=il_g*(4*il_g-i-2)+n      ! panel 3
        if(n>3*il_g)iq=il_g*(5*il_g-i-3)+n                  ! panel 4       
        if (em_g(iq).gt.emmin) then       
        do n1=n,4*il_g
         if(n1<=2*il_g)iq1=il_g*(n1-1)+i                  ! panel 0,1
         if(n1>2*il_g.and.n1<=3*il_g)iq1=il_g*(4*il_g-i-2)+n1 ! panel 3
         if(n1>3*il_g)iq1=il_g*(5*il_g-i-3)+n1              ! panel 4
         dist1=abs(yy(iq)-yy(iq1))
         dist=min(dist1,yymax-dist1)
         wt=exp(-4.5*dist*dist*cutoff2)
         wt1=wt/em_g(iq1)
         wt=wt/em_g(iq)
         if(n==n1)wt1=0.  ! so as not to add in twice
         sumwt(iq)=sumwt(iq)+wt1
         sumwt(iq1)=sumwt(iq1)+wt
!        producing "y-filtered" version of pslb-psl etc
         psls(iq)=psls(iq)+wt1*psla(iq1)
         psls(iq1)=psls(iq1)+wt*psla(iq)
         do k=kbotdav,kl
          qgg(iq,k)=qgg(iq,k)+wt1*qa(iq1,k)
          qgg(iq1,k)=qgg(iq1,k)+wt*qa(iq,k)
          tt(iq,k)=tt(iq,k)+wt1*ta(iq1,k)
          tt(iq1,k)=tt(iq1,k)+wt*ta(iq,k)
          uu(iq,k)=uu(iq,k)+wt1*ua(iq1,k)
          uu(iq1,k)=uu(iq1,k)+wt*ua(iq,k)
          vv(iq,k)=vv(iq,k)+wt1*va(iq1,k)
          vv(iq1,k)=vv(iq1,k)+wt*va(iq,k)
          ww(iq,k)=ww(iq,k)+wt1*wa(iq1,k)
          ww(iq1,k)=ww(iq1,k)+wt*wa(iq,k)
         enddo  ! k loop
        enddo   ! n1 loop
        else
          sumwt(iq)=1.
        end if
       enddo    ! n loop
      enddo     ! i loop
      if(nud_uv==-1)then
        do iq=1,ifull_g
         psls2(iq)=psls2(iq)+psls(iq)/sumwt(iq)
         do k=kbotdav,kl
          qgg2(iq,k)=qgg2(iq,k)+qgg(iq,k)/sumwt(iq)
          tt2(iq,k)=tt2(iq,k)+tt(iq,k)/sumwt(iq)
          uu2(iq,k)=uu2(iq,k)+uu(iq,k)/sumwt(iq)
          vv2(iq,k)=vv2(iq,k)+vv(iq,k)/sumwt(iq)
          ww2(iq,k)=ww2(iq,k)+ww(iq,k)/sumwt(iq)
         enddo
        enddo
      else  ! original fast scheme
        do iq=1,ifull_g
         if(sumwt(iq).ne.1.e-20)then
           psla(iq)=psls(iq)/sumwt(iq)
           do k=kbotdav,kl
            qa(iq,k)=qgg(iq,k)/sumwt(iq)
            ta(iq,k)=tt(iq,k)/sumwt(iq)
            ua(iq,k)=uu(iq,k)/sumwt(iq)
            va(iq,k)=vv(iq,k)/sumwt(iq)
            wa(iq,k)=ww(iq,k)/sumwt(iq)
           enddo
         endif  ! (sumwt(iq).ne.1.e-20)
        enddo
      endif  ! (nud_uv==-1) .. else ..

      if(mbd.ge.0) then
       qgg(1:ifull_g,kbotdav:kl)=0.
       tt(1:ifull_g,kbotdav:kl)=0.
       uu(1:ifull_g,kbotdav:kl)=0.
       vv(1:ifull_g,kbotdav:kl)=0.
       ww(1:ifull_g,kbotdav:kl)=0.
       psls(1:ifull_g)=0.
       sumwt(1:ifull_g)=1.e-20   ! for undefined panels
    
       do j=1,il_g                ! doing "z"-filter on panels 0,2,3,5
        zzmax=zz(5*il_g*il_g+il_g*j)
        print *,'j,zzmax ',j,zzmax
        do n=1,4*il_g
         if(n<=il_g)iq=il_g*(j-1)+n                     ! panel 0
         if(n>il_g.and.n<=3*il_g)iq=il_g*(il_g+n-1)+il_g+1-j  ! panel 2,3
         if(n>3*il_g)iq=il_g*(5*il_g+j-4)+n               ! panel 5
        
         if (em_g(iq).gt.emmin) then ! MJT
        
         do n1=n,4*il_g
          if(n1<=il_g)iq1=il_g*(j-1)+n1                     ! panel 0
          if(n1>il_g.and.n1<=3*il_g)iq1=il_g*(il_g+n1-1)+il_g+1-j ! panel 2,3
          if(n1>3*il_g)iq1=il_g*(5*il_g+j-4)+n1               ! panel 5
          dist1=abs(zz(iq)-zz(iq1))
          dist=min(dist1,zzmax-dist1)
          wt=exp(-4.5*dist*dist*cutoff2)
          wt1=wt/em_g(iq1)
          wt=wt/em_g(iq)
          if(n==n1)wt1=0.  ! so as not to add in twice
          sumwt(iq)=sumwt(iq)+wt1
          sumwt(iq1)=sumwt(iq1)+wt
!         producing "z"-filtered version of pslb-psl etc
          psls(iq)=psls(iq)+wt1*psla(iq1)
          psls(iq1)=psls(iq1)+wt*psla(iq)
          do k=kbotdav,kl
           qgg(iq,k)=qgg(iq,k)+wt1*qa(iq1,k)
           qgg(iq1,k)=qgg(iq1,k)+wt*qa(iq,k)
           tt(iq,k)=tt(iq,k)+wt1*ta(iq1,k)
           tt(iq1,k)=tt(iq1,k)+wt*ta(iq,k)
           uu(iq,k)=uu(iq,k)+wt1*ua(iq1,k)
           uu(iq1,k)=uu(iq1,k)+wt*ua(iq,k)
           vv(iq,k)=vv(iq,k)+wt1*va(iq1,k)
           vv(iq1,k)=vv(iq1,k)+wt*va(iq,k)
           ww(iq,k)=ww(iq,k)+wt1*wa(iq1,k)
           ww(iq1,k)=ww(iq1,k)+wt*wa(iq,k)
          enddo  ! k loop
         enddo   ! n1 loop
         else
           sumwt(iq)=1.
         end if
        enddo    ! n loop
       enddo     ! j loop      
      if(nud_uv==-1)then
        print *,'in nestinb nud_uv ',nud_uv
        do iq=1,ifull_g
         psls2(iq)=psls2(iq)+psls(iq)/sumwt(iq)
         do k=kbotdav,kl
          qgg2(iq,k)=qgg2(iq,k)+qgg(iq,k)/sumwt(iq)
          tt2(iq,k)=tt2(iq,k)+tt(iq,k)/sumwt(iq)
          uu2(iq,k)=uu2(iq,k)+uu(iq,k)/sumwt(iq)
          vv2(iq,k)=vv2(iq,k)+vv(iq,k)/sumwt(iq)
          ww2(iq,k)=ww2(iq,k)+ww(iq,k)/sumwt(iq)
         enddo
        enddo
        psla(1:ifull_g)=.5*psls2(1:ifull_g)
        qa(1:ifull_g,kbotdav:kl)=.5*qgg2(1:ifull_g,kbotdav:kl)
        ta(1:ifull_g,kbotdav:kl)=.5*tt2(1:ifull_g,kbotdav:kl)
        ua(1:ifull_g,kbotdav:kl)=.5*uu2(1:ifull_g,kbotdav:kl)
        va(1:ifull_g,kbotdav:kl)=.5*vv2(1:ifull_g,kbotdav:kl)
        wa(1:ifull_g,kbotdav:kl)=.5*ww2(1:ifull_g,kbotdav:kl)
      else  ! original fast scheme
        print *,'in nestinb  nud_uv ',nud_uv
        do iq=1,ifull_g
         if(sumwt(iq).ne.1.e-20)then
           psla(iq)=psls(iq)/sumwt(iq)
           do k=kbotdav,kl
            qa(iq,k)=qgg(iq,k)/sumwt(iq)
            ta(iq,k)=tt(iq,k)/sumwt(iq)
            ua(iq,k)=uu(iq,k)/sumwt(iq)
            va(iq,k)=vv(iq,k)/sumwt(iq)
            wa(iq,k)=ww(iq,k)/sumwt(iq)
           enddo
         endif  ! (sumwt(iq).ne.1.e-20)
        enddo
      endif  ! (nud_uv==-1) .. else ..
      end if ! (mbd.ge.0)

      return
      end subroutine fastspec


      !---------------------------------------------------------------------------------
      ! Slow 2D spectral downscaling - MPI version
      subroutine slowspecmpi(myid,cin,psls,uu,vv,ww,tt,qgg)
      
      implicit none
      
      include 'newmpar.h'   ! ifull_g,kl
      include 'const_phys.h' ! rearth,pi,tpi
      include 'map_g.h'     ! em_g
      include 'parm.h'      ! ds,kbotdav
      include 'xyzinfo_g.h' ! x_g,y_g,z_g
      include 'mpif.h'

      integer, intent(in) :: myid
      integer :: iq,ns,ne,k,itag=0,ierr,iproc,iy
      integer, dimension(MPI_STATUS_SIZE) :: status
      real, intent(in) :: cin
      real, dimension(ifull_g), intent(inout) :: psls
      real, dimension(ifull_g,kbotdav:kl), intent(inout) :: uu,vv,ww
      real, dimension(ifull_g,kbotdav:kl), intent(inout) :: tt,qgg
      real, dimension(ifull_g) :: pp,r
      real, dimension(ifull_g,kbotdav:kl) :: pu,pv,pw,pt,pq
      real, dimension(ifull_g*(kl-kbotdav+1)) :: dd
      real :: cq,psum

      cq=sqrt(4.5)*cin

      if (myid == 0) then
        if(nmaxpr==1) print *,"Send arrays to all processors"
        if(nud_p>0)then
          call MPI_Bcast(psls(:),ifull_g,MPI_REAL,0,
     &           MPI_COMM_WORLD,ierr)
        end if
        iy=ifull_g*(kl-kbotdav+1)
        if(nud_uv>0)then
          dd(1:iy)=reshape(uu(:,:),(/iy/))
          call MPI_Bcast(dd(1:iy),iy,MPI_REAL,0,
     &           MPI_COMM_WORLD,ierr)    
          dd(1:iy)=reshape(vv(:,:),(/iy/))
          call MPI_Bcast(dd(1:iy),iy,MPI_REAL,0,
     &           MPI_COMM_WORLD,ierr)    
          dd(1:iy)=reshape(ww(:,:),(/iy/))
          call MPI_Bcast(dd(1:iy),iy,MPI_REAL,0,
     &           MPI_COMM_WORLD,ierr)    
        end if
        if(nud_t>0)then
          dd(1:iy)=reshape(tt(:,:),(/iy/))
          call MPI_Bcast(dd(1:iy),iy,MPI_REAL,0,
     &           MPI_COMM_WORLD,ierr)    
        end if
        if(nud_q>0)then
          dd(1:iy)=reshape(qgg(:,:),(/iy/))
          call MPI_Bcast(dd(1:iy),iy,MPI_REAL,0,
     &           MPI_COMM_WORLD,ierr)    
        end if
      else
        if(nud_p>0)then
         call MPI_Bcast(psls(:),ifull_g,MPI_REAL,0,
     &          MPI_COMM_WORLD,ierr)
        end if
        iy=ifull_g*(kl-kbotdav+1)
        if(nud_uv>0)then
          call MPI_Bcast(dd(1:iy),iy,MPI_REAL,0,
     &           MPI_COMM_WORLD,ierr)
          uu(:,:)=reshape(dd(1:iy),(/ifull_g,kl-kbotdav+1/))
          call MPI_Bcast(dd(1:iy),iy,MPI_REAL,0,
     &           MPI_COMM_WORLD,ierr)
          vv(:,:)=reshape(dd(1:iy),(/ifull_g,kl-kbotdav+1/))
          call MPI_Bcast(dd(1:iy),iy,MPI_REAL,0,
     &           MPI_COMM_WORLD,ierr)
          ww(:,:)=reshape(dd(1:iy),(/ifull_g,kl-kbotdav+1/))
        end if
        if(nud_t>0)then
          call MPI_Bcast(dd(1:iy),iy,MPI_REAL,0,
     &           MPI_COMM_WORLD,ierr)
          tt(:,:)=reshape(dd(1:iy),(/ifull_g,kl-kbotdav+1/))
        end if
        if(nud_q>0)then
          call MPI_Bcast(dd(1:iy),iy,MPI_REAL,0,
     &           MPI_COMM_WORLD,ierr)
          qgg(:,:)=reshape(dd(1:iy),(/ifull_g,kl-kbotdav+1/))
        end if
      end if
    
      call procdiv(ns,ne,ifull_g,nproc,myid)
      if ((myid == 0).and.(nmaxpr==1)) print *,"Start 2D filter"
      
      do iq=ns,ne
        r(:)=x_g(iq)*x_g(:)+y_g(iq)*y_g(:)+z_g(iq)*z_g(:)
        r(:)=acos(max(min(r(:),1.),-1.))
        r(:)=exp(-(cq*r(:))**2)/(em_g(:)**2)
        psum=sum(r(:))
        if (nud_p>0) then
          pp(iq)=sum(r(:)*psls(:))/psum
        end if
        if (nud_uv>0) then
          do k=kbotdav,kl
            pu(iq,k)=sum(r(:)*uu(:,k))/psum
            pv(iq,k)=sum(r(:)*vv(:,k))/psum
            pw(iq,k)=sum(r(:)*ww(:,k))/psum
          end do        
        end if
        if (nud_t>0) then
          do k=kbotdav,kl
            pt(iq,k)=sum(r(:)*tt(:,k))/psum
          end do
        end if
        if (nud_q>0) then
          do k=kbotdav,kl
            pq(iq,k)=sum(r(:)*qgg(:,k))/psum
          end do
        end if
      end do
 
      if (nud_p>0) then
        psls(ns:ne)=pp(ns:ne)
      end if
      if (nud_uv>0) then
        uu(ns:ne,:)=pu(ns:ne,:)
        vv(ns:ne,:)=pv(ns:ne,:)
        ww(ns:ne,:)=pw(ns:ne,:)
      end if
      if (nud_t>0) then
        tt(ns:ne,:)=pt(ns:ne,:)
      end if
      if (nud_q>0) then
        qgg(ns:ne,:)=pq(ns:ne,:)
      end if
          
      if ((myid == 0).and.(nmaxpr==1)) print *,"End 2D filter"

      itag=itag+1
      if (myid == 0) then
        if (nmaxpr==1) print *,"Receive arrays from all processors"
        do iproc=1,nproc-1
          call procdiv(ns,ne,ifull_g,nproc,iproc)
          if(nud_p>0)call MPI_Recv(psls(ns:ne),ne-ns+1,MPI_REAL,iproc
     &                      ,itag,MPI_COMM_WORLD,status,ierr)
          iy=(ne-ns+1)*(kl-kbotdav+1)
          if(nud_uv>0)then
            call MPI_Recv(dd(1:iy),iy,MPI_REAL,iproc,itag,
     &             MPI_COMM_WORLD,status,ierr)
            uu(ns:ne,:)=reshape(dd(1:iy),(/ne-ns+1,kl-kbotdav+1/))
            call MPI_Recv(dd(1:iy),iy,MPI_REAL,iproc,itag,
     &             MPI_COMM_WORLD,status,ierr)
            vv(ns:ne,:)=reshape(dd(1:iy),(/ne-ns+1,kl-kbotdav+1/))
            call MPI_Recv(dd(1:iy),iy,MPI_REAL,iproc,itag,
     &             MPI_COMM_WORLD,status,ierr)
            ww(ns:ne,:)=reshape(dd(1:iy),(/ne-ns+1,kl-kbotdav+1/))
          end if
          if(nud_t>0)then
            call MPI_Recv(dd(1:iy),iy,MPI_REAL,iproc,itag,
     &             MPI_COMM_WORLD,status,ierr)
            tt(ns:ne,:)=reshape(dd(1:iy),(/ne-ns+1,kl-kbotdav+1/))
          end if
          if(nud_q>0)then
            call MPI_Recv(dd(1:iy),iy,MPI_REAL,iproc,itag,
     &             MPI_COMM_WORLD,status,ierr)
            qgg(ns:ne,:)=reshape(dd(1:iy),(/ne-ns+1,kl-kbotdav+1/))
          end if
        end do
      else
        if(nud_p>0) call MPI_SSend(psls(ns:ne),ne-ns+1,MPI_REAL,0,
     &                     itag,MPI_COMM_WORLD,ierr)
        iy=(ne-ns+1)*(kl-kbotdav+1)
        if(nud_uv>0)then
          dd(1:iy)=reshape(uu(ns:ne,:),(/iy/))
          call MPI_SSend(dd(1:iy),iy,MPI_REAL,0,itag,
     &           MPI_COMM_WORLD,ierr)
          dd(1:iy)=reshape(vv(ns:ne,:),(/iy/))
          call MPI_SSend(dd(1:iy),iy,MPI_REAL,0,itag,
     &           MPI_COMM_WORLD,ierr)
          dd(1:iy)=reshape(ww(ns:ne,:),(/iy/))
          call MPI_SSend(dd(1:iy),iy,MPI_REAL,0,itag,
     &           MPI_COMM_WORLD,ierr)
        end if
        if(nud_t>0)then
          dd(1:iy)=reshape(tt(ns:ne,:),(/iy/))
          call MPI_SSend(dd(1:iy),iy,MPI_REAL,0,itag,
     &           MPI_COMM_WORLD,ierr)
        end if
        if(nud_q>0)then
          dd(1:iy)=reshape(qgg(ns:ne,:),(/iy/))
          call MPI_SSend(dd(1:iy),iy,MPI_REAL,0,itag,
     &           MPI_COMM_WORLD,ierr)
        end if
      end if

      return
      end subroutine slowspecmpi
      !---------------------------------------------------------------------------------

      !---------------------------------------------------------------------------------
      ! Four pass spectral downscaling (symmetric)
      ! Used when mod(6,nproc).ne.0 and mod(nproc,6).ne.0 since it is slower than specfastmpi
      subroutine fourspecmpi(myid,cin,psls,uu,vv,ww,tt,qgg)
      
      implicit none
      
      include 'newmpar.h'    ! ifull_g,kl
      include 'parm.h'       ! kbotdav
      
      integer, intent(in) :: myid
      integer :: pn,px,hproc,mproc,ns,ne,npta
      real, intent(in) :: cin
      real, dimension(ifull_g), intent(inout) :: psls
      real, dimension(ifull_g,kbotdav:kl), intent(inout) :: uu,vv,ww
      real, dimension(ifull_g,kbotdav:kl), intent(inout) :: tt,qgg
      
      npta=1                              ! number of panels per processor
      mproc=nproc                         ! number of processors per panel
      pn=0                                ! start panel
      px=5                                ! end panel
      hproc=0                             ! host processor for panel
      call procdiv(ns,ne,il_g,nproc,myid) ! number of rows per processor

      call spechost(myid,mproc,hproc,npta,pn,px,ns,ne,cin,psls,uu,vv,
     &       ww,tt,qgg)
          
      return
      end subroutine fourspecmpi
      !---------------------------------------------------------------------------------


      !---------------------------------------------------------------------------------
      ! Four pass spectral downscaling (symmetric)
      ! MPI optimised for magic processor numbers 1,2,3,6,12,18,24,30,36,...
      ! (only works for mod(6,nproc)==0 or mod(nproc,6)==0)
      subroutine specfastmpi(myid,cin,psls,uu,vv,ww,tt,qgg)
      
      implicit none
      
      include 'newmpar.h'    ! ifull_g,kl
      include 'parm.h'       ! kbotdav
      
      integer, intent(in) :: myid
      integer :: pn,px,hproc,mproc,ns,ne,npta
      real, intent(in) :: cin
      real, dimension(ifull_g), intent(inout) :: psls
      real, dimension(ifull_g,kbotdav:kl), intent(inout) :: uu,vv,ww
      real, dimension(ifull_g,kbotdav:kl), intent(inout) :: tt,qgg
      
      npta=max(6/nproc,1)                       ! number of panels per processor
      mproc=max(nproc/6,1)                      ! number of processors per panel
      pn=myid*npta/mproc                        ! start panel
      px=(myid+mproc)*npta/mproc-1              ! end panel
      hproc=pn*mproc/npta                       ! host processor for panel
      call procdiv(ns,ne,il_g,mproc,myid-hproc) ! number of rows per processor

      call spechost(myid,mproc,hproc,npta,pn,px,ns,ne,cin,psls,uu,vv,
     &       ww,tt,qgg)

      end subroutine specfastmpi
      !---------------------------------------------------------------------------------


      !---------------------------------------------------------------------------------
      ! This is the main routine for the scale-selective filter
      subroutine spechost(myid,mproc,hproc,npta,pn,px,ns,ne,cin,psls,
     &                    uu,vv,ww,tt,qgg)
      
      implicit none
      
      include 'newmpar.h'    ! ifull_g,kl
      include 'const_phys.h' ! rearth,pi,tpi
      include 'parm.h'       ! ds,kbotdav
      include 'mpif.h'       ! MPI
      
      integer, intent(in) :: myid,mproc,hproc,npta,pn,px,ns,ne
      integer :: k,ppass,qpass,iy,ppn,ppx,nne,nns,iproc,itag=0,ierr
      integer :: n,a,b,c
      integer, dimension(MPI_STATUS_SIZE) :: status
      integer, parameter :: til=il_g*il_g
      integer, parameter, dimension(0:5) :: qaps=(/0,3,1,4,2,5/)
      integer, parameter, dimension(0:5) :: qms=qaps*til
      real, intent(in) :: cin
      real, dimension(ifull_g), intent(inout) :: psls
      real, dimension(ifull_g,kbotdav:kl), intent(inout) :: uu,vv,ww
      real, dimension(ifull_g,kbotdav:kl), intent(inout) :: tt,qgg
      real, dimension(ifull_g) :: qp,qsum,sp,ssum,zp
      real, dimension(ifull_g,kbotdav:kl) :: qu,qv,qw,qt,qq
      real, dimension(ifull_g,kbotdav:kl) :: su,sv,sw,st,sq
      real, dimension(ifull_g,kbotdav:kl) :: zu,zv,zw,zt,zq
      real, dimension(ifull_g*(kl-kbotdav+1)) :: dd
      real :: cq

      cq=sqrt(4.5)*cin ! filter length scale

      if (myid == 0) then
        if (nmaxpr==1) print *,"Send arrays to all processors"
        if(nud_p>0)then
          call MPI_Bcast(psls(:),ifull_g,MPI_REAL,0,
     &           MPI_COMM_WORLD,ierr)
        end if
        iy=ifull_g*(kl-kbotdav+1)
        if(nud_uv>0)then
          dd(1:iy)=reshape(uu(:,:),(/iy/))
          call MPI_Bcast(dd(1:iy),iy,MPI_REAL,0,
     &           MPI_COMM_WORLD,ierr)    
          dd(1:iy)=reshape(vv(:,:),(/iy/))
          call MPI_Bcast(dd(1:iy),iy,MPI_REAL,0,
     &           MPI_COMM_WORLD,ierr)    
          dd(1:iy)=reshape(ww(:,:),(/iy/))
          call MPI_Bcast(dd(1:iy),iy,MPI_REAL,0,
     &           MPI_COMM_WORLD,ierr)    
        end if
        if(nud_t>0)then
          dd(1:iy)=reshape(tt(:,:),(/iy/))
          call MPI_Bcast(dd(1:iy),iy,MPI_REAL,0,
     &           MPI_COMM_WORLD,ierr)    
        end if
        if(nud_q>0)then
          dd(1:iy)=reshape(qgg(:,:),(/iy/))
          call MPI_Bcast(dd(1:iy),iy,MPI_REAL,0,
     &           MPI_COMM_WORLD,ierr)    
        end if
      else
        if(nud_p>0)then
         call MPI_Bcast(psls(:),ifull_g,MPI_REAL,0,
     &                     MPI_COMM_WORLD,ierr)
        end if
        iy=ifull_g*(kl-kbotdav+1)
        if(nud_uv>0)then
          call MPI_Bcast(dd(1:iy),iy,MPI_REAL,0,
     &           MPI_COMM_WORLD,ierr)
          uu(:,:)=reshape(dd(1:iy),(/ifull_g,kl-kbotdav+1/))
          call MPI_Bcast(dd(1:iy),iy,MPI_REAL,0,
     &           MPI_COMM_WORLD,ierr)
          vv(:,:)=reshape(dd(1:iy),(/ifull_g,kl-kbotdav+1/))
          call MPI_Bcast(dd(1:iy),iy,MPI_REAL,0,
     &           MPI_COMM_WORLD,ierr)
          ww(:,:)=reshape(dd(1:iy),(/ifull_g,kl-kbotdav+1/))
        end if
        if(nud_t>0)then
          call MPI_Bcast(dd(1:iy),iy,MPI_REAL,0,
     &           MPI_COMM_WORLD,ierr)
          tt(:,:)=reshape(dd(1:iy),(/ifull_g,kl-kbotdav+1/))
        end if
        if(nud_q>0)then
          call MPI_Bcast(dd(1:iy),iy,MPI_REAL,0,
     &           MPI_COMM_WORLD,ierr)
          qgg(:,:)=reshape(dd(1:iy),(/ifull_g,kl-kbotdav+1/))
        end if
      end if
      
      if (ns.gt.ne) return
      if ((myid==0).and.(nmaxpr==1)) print *,"Start 1D filter"

      do qpass=pn,px
        ppass=qaps(qpass)

        qsum(:)=1.
        if (nud_p>0) then
          qp(:)=psls(:)
        end if
        if (nud_uv>0) then
          qu(:,:)=uu(:,:)
          qv(:,:)=vv(:,:)
          qw(:,:)=ww(:,:)
        end if
        if (nud_t>0) then
          qt(:,:)=tt(:,:)
        end if
        if (nud_q>0) then
          qq(:,:)=qgg(:,:)
        end if

        ! computations for the local processor group
        call speclocal(myid,mproc,hproc,ns,ne,cq,ppass,qsum,qp,
     &         qu,qv,qw,qt,qq,ssum,sp,su,sv,sw,st,sq)
        
        nns=qms(qpass)+1
        nne=qms(qpass)+til
        if (nud_p>0) then
          zp(nns:nne)=qp(nns:nne)/qsum(nns:nne)
        end if
        if (nud_uv>0) then
          do k=kbotdav,kl
            zu(nns:nne,k)=qu(nns:nne,k)/qsum(nns:nne)
            zv(nns:nne,k)=qv(nns:nne,k)/qsum(nns:nne)
            zw(nns:nne,k)=qw(nns:nne,k)/qsum(nns:nne)
          end do
        end if
        if (nud_t>0) then
          do k=kbotdav,kl
            zt(nns:nne,k)=qt(nns:nne,k)/qsum(nns:nne)
          end do
        end if
        if (nud_q>0) then
          do k=kbotdav,kl
            zq(nns:nne,k)=qq(nns:nne,k)/qsum(nns:nne)
          end do
        end if
        
      end do

      if ((myid==0).and.(nmaxpr==1)) print *,"End 1D filter"

      itag=itag+1
      if (myid == 0) then
        if (nmaxpr==1) print *,"Receive arrays from all host processors"
        do iproc=mproc,nproc-1,mproc
          ppn=iproc*npta/mproc
          ppx=(iproc+mproc)*npta/mproc-1
          iy=npta*til
          a=til
          c=-til*ppn
          if(nud_p>0)then
            call MPI_Recv(dd(1:iy),iy,MPI_REAL,
     &             iproc,itag,MPI_COMM_WORLD,status,ierr)
            do qpass=ppn,ppx
              do n=1,til
                zp(n+qms(qpass))=dd(n+a*qpass+c)
              end do
            end do
          end if
          iy=npta*til*(kl-kbotdav+1)
          b=npta*til
          c=-til*(ppn+npta*kbotdav)
          if(nud_uv>0)then
            call MPI_Recv(dd(1:iy),iy,MPI_REAL,iproc,itag,
     &             MPI_COMM_WORLD,status,ierr)
            do k=kbotdav,kl
              do qpass=ppn,ppx
                do n=1,til
                  zu(n+qms(qpass),k)=dd(n+a*qpass+b*k+c)
                end do
              end do
            end do
            call MPI_Recv(dd(1:iy),iy,MPI_REAL,iproc,itag,
     &             MPI_COMM_WORLD,status,ierr)
            do k=kbotdav,kl
              do qpass=ppn,ppx
                do n=1,til
                  zv(n+qms(qpass),k)=dd(n+a*qpass+b*k+c)
                end do
              end do
            end do
            call MPI_Recv(dd(1:iy),iy,MPI_REAL,iproc,itag,
     &             MPI_COMM_WORLD,status,ierr)
            do k=kbotdav,kl
              do qpass=ppn,ppx
                do n=1,til
                  zw(n+qms(qpass),k)=dd(n+a*qpass+b*k+c)
                end do
              end do
            end do
          end if
          if(nud_t>0)then
            call MPI_Recv(dd(1:iy),iy,MPI_REAL,iproc,itag,
     &             MPI_COMM_WORLD,status,ierr)
            do k=kbotdav,kl
              do qpass=ppn,ppx
                do n=1,til
                  zt(n+qms(qpass),k)=dd(n+a*qpass+b*k+c)
                end do
              end do
            end do
          end if
          if(nud_q>0)then
            call MPI_Recv(dd(1:iy),iy,MPI_REAL,iproc,itag,
     &             MPI_COMM_WORLD,status,ierr)
            do k=kbotdav,kl
              do qpass=ppn,ppx
                do n=1,til
                  zq(n+qms(qpass),k)=dd(n+a*qpass+b*k+c)
                end do
              end do
            end do
          end if
        end do
      elseif (myid==hproc) then
        iy=npta*til
        a=til
        c=-til*pn
        if(nud_p>0)then
          do qpass=pn,px
            do n=1,til
              dd(n+a*qpass+c)=zp(n+qms(qpass))
            end do
          end do
          call MPI_SSend(dd(1:iy),iy,MPI_REAL,0,
     &           itag,MPI_COMM_WORLD,ierr)
        end if
        iy=npta*til*(kl-kbotdav+1)
        b=npta*til
        c=-til*(pn+npta*kbotdav)
        if(nud_uv>0)then
          do k=kbotdav,kl
            do qpass=pn,px
              do n=1,til
                dd(n+a*qpass+b*k+c)=zu(n+qms(qpass),k)
              end do
            end do
          end do
          call MPI_SSend(dd(1:iy),iy,MPI_REAL,0,itag,
     &           MPI_COMM_WORLD,ierr)
          do k=kbotdav,kl
            do qpass=pn,px
              do n=1,til
                dd(n+a*qpass+b*k+c)=zv(n+qms(qpass),k)
              end do
            end do
          end do
          call MPI_SSend(dd(1:iy),iy,MPI_REAL,0,itag,
     &           MPI_COMM_WORLD,ierr)
          do k=kbotdav,kl
            do qpass=pn,px
              do n=1,til
                dd(n+a*qpass+b*k+c)=zw(n+qms(qpass),k)
              end do
            end do
          end do
          call MPI_SSend(dd(1:iy),iy,MPI_REAL,0,itag,
     &           MPI_COMM_WORLD,ierr)
        end if
        if(nud_t>0)then
          do k=kbotdav,kl
            do qpass=pn,px
              do n=1,til
                dd(n+a*qpass+b*k+c)=zt(n+qms(qpass),k)
              end do
            end do
          end do
          call MPI_SSend(dd(1:iy),iy,MPI_REAL,0,itag,
     &           MPI_COMM_WORLD,ierr)
        end if
        if(nud_q>0)then
          do k=kbotdav,kl
            do qpass=pn,px
              do n=1,til
                dd(n+a*qpass+b*k+c)=zq(n+qms(qpass),k)
              end do
            end do
          end do
          call MPI_SSend(dd(1:iy),iy,MPI_REAL,0,itag,
     &           MPI_COMM_WORLD,ierr)
        end if
      end if      
      
      if (nud_p>0) then
        psls(:)=zp(:)
      end if
      if (nud_uv>0) then
        uu(:,:)=zu(:,:)
        vv(:,:)=zv(:,:)
        ww(:,:)=zw(:,:)
      end if
      if (nud_t>0) then
        tt(:,:)=zt(:,:)
      end if
      if (nud_q>0) then
        qgg(:,:)=zq(:,:)
      end if
      
      return
      end subroutine spechost
      !---------------------------------------------------------------------------------
      
      !---------------------------------------------------------------------------------
      ! This code runs between the local processors
      ! Code was moved to this subroutine to help the compiler vectorise the code
      subroutine speclocal(myid,mproc,hproc,ns,ne,cq,ppass,qsum,
     &             qp,qu,qv,qw,qt,qq,ssum,sp,su,sv,sw,st,sq)
      implicit none
      
      include 'newmpar.h'    ! ifull_g,kl
      include 'const_phys.h' ! pi
      include 'map_g.h'      ! em_g
      include 'parm.h'       ! kbotdav
      include 'xyzinfo_g.h'  ! x_g,y_g,z_g
      include 'mpif.h'       ! MPI
      
      integer, intent(in) :: myid,mproc,hproc,ns,ne,ppass
      integer :: j,k,n,ipass,kpass,iy
      integer :: iproc,itag=0,ierr
      integer :: nne,nns,me
      integer :: a,b,c,d,ips
      integer, save :: pold=-1
      integer, dimension(MPI_STATUS_SIZE) :: status
      integer, dimension(4*il_g,il_g,0:3) :: igrd
      integer, parameter, dimension(0:3) ::
     &  maps=(/ il_g, il_g, 4*il_g, 3*il_g /)
      integer, parameter, dimension(0:5) :: pair=(/3,4,5,0,1,2/)
      integer, parameter, dimension(2:3) :: kn=(/0,3/)
      integer, parameter, dimension(2:3) :: kx=(/2,3/)
      real, intent(in) :: cq
      real, dimension(ifull_g), intent(inout) :: qp,qsum
      real, dimension(ifull_g), intent(inout) :: sp,ssum
      real, dimension(ifull_g,kbotdav:kl), intent(inout) :: qu,qv,qw
      real, dimension(ifull_g,kbotdav:kl), intent(inout) :: qt,qq
      real, dimension(ifull_g,kbotdav:kl), intent(inout) :: su,sv,sw
      real, dimension(ifull_g,kbotdav:kl), intent(inout) :: st,sq
      real, dimension(4*il_g,kbotdav:kl) :: pu,pv,pw,pt,pq
      real, dimension(4*il_g,kbotdav:kl) :: au,av,aw,at,aq
      real, dimension(4*il_g) :: pp,ap,psum,asum,ra,ma,xa,ya,za
      real, dimension(ifull_g*(kl-kbotdav+1)) :: dd
      
      if (pold.lt.0) pold=ppass
      
      if (ppass.eq.pair(pold)) then
        if ((myid==0).and.(nmaxpr==1)) then
          print *,"Use stored convolution for pass 0 and 1"
        end if
        ips=2
        qsum(:)=ssum(:)
        if (nud_p>0) then
          qp(:)=sp(:)
        end if
        if (nud_uv>0) then
          qu(:,:)=su(:,:)
          qv(:,:)=sv(:,:)
          qw(:,:)=sw(:,:)
        end if
        if (nud_t>0) then
          qt(:,:)=st(:,:)
        end if
        if (nud_q>0) then
          qq(:,:)=sq(:,:)
        end if
      else
        ips=0
      end if
      
      do ipass=ips,3
        me=maps(ipass)
        call getiqa(igrd(1:me,1:il_g,ipass),me,ipass,ppass,il_g)

          if (ipass.eq.3) then
            itag=itag+1
            if ((myid==0).and.(nmaxpr==1)) then
              print *,"Recieve arrays from local host"
            end if
            if (myid==hproc) then
              do iproc=hproc+1,mproc+hproc-1
                call procdiv(nns,nne,il_g,mproc,iproc-hproc)
                if (nns.gt.nne) exit
                iy=me*(nne-nns+1)
                a=me
                d=-me*nns
                do j=nns,nne
                  do n=1,me
                    dd(n+a*j+d)=qsum(igrd(n,j,ipass))
                  end do
                end do
                call MPI_SSend(dd(1:iy),iy,MPI_REAL,iproc,itag,
     &                 MPI_COMM_WORLD,ierr)    
                if(nud_p>0)then
                  do j=nns,nne
                    do n=1,me
                      dd(n+a*j+d)=qp(igrd(n,j,ipass))
                    end do
                  end do
                  call MPI_SSend(dd(1:iy),iy,MPI_REAL,iproc,itag,
     &                   MPI_COMM_WORLD,ierr)
                end if
                iy=me*(nne-nns+1)*(kl-kbotdav+1)
                b=me*(nne-nns+1)
                d=-me*(nns+kbotdav*(nne-nns+1))
                if(nud_uv>0)then
                  do k=kbotdav,kl    
                    do j=nns,nne
                      do n=1,me
                        dd(n+a*j+b*k+d)=qu(igrd(n,j,ipass),k)
                      end do
                    end do
                  end do
                  call MPI_SSend(dd(1:iy),iy,MPI_REAL,iproc,itag,
     &                   MPI_COMM_WORLD,ierr)
                  do k=kbotdav,kl    
                    do j=nns,nne
                      do n=1,me
                        dd(n+a*j+b*k+d)=qv(igrd(n,j,ipass),k)
                      end do
                    end do
                  end do
                  call MPI_SSend(dd(1:iy),iy,MPI_REAL,iproc,itag,
     &                   MPI_COMM_WORLD,ierr)
                  do k=kbotdav,kl    
                    do j=nns,nne
                      do n=1,me
                        dd(n+a*j+b*k+d)=qw(igrd(n,j,ipass),k)
                      end do
                    end do
                  end do
                  call MPI_SSend(dd(1:iy),iy,MPI_REAL,iproc,itag,
     &                   MPI_COMM_WORLD,ierr)
                end if
                if(nud_t>0)then
                  do k=kbotdav,kl    
                    do j=nns,nne
                      do n=1,me
                        dd(n+a*j+b*k+d)=qt(igrd(n,j,ipass),k)
                      end do
                    end do
                  end do
                  call MPI_SSend(dd(1:iy),iy,MPI_REAL,iproc,itag,
     &                   MPI_COMM_WORLD,ierr)
                end if
                if(nud_q>0)then
                  do k=kbotdav,kl    
                    do j=nns,nne
                      do n=1,me
                        dd(n+a*j+b*k+d)=qq(igrd(n,j,ipass),k)
                      end do
                    end do
                  end do
                  call MPI_SSend(dd(1:iy),iy,MPI_REAL,iproc,itag,
     &                   MPI_COMM_WORLD,ierr)
                end if
              end do
            else
              iy=me*(ne-ns+1)
              a=me
              d=-me*ns
              call MPI_Recv(dd(1:iy),iy,MPI_REAL,hproc,itag,
     &               MPI_COMM_WORLD,status,ierr)
              do j=ns,ne
                do n=1,me
                  qsum(igrd(n,j,ipass))=dd(n+a*j+d)
                end do
              end do
              if(nud_p>0)then
                call MPI_Recv(dd(1:iy),iy,MPI_REAL,hproc,itag,
     &                 MPI_COMM_WORLD,status,ierr)
                do j=ns,ne
                  do n=1,me
                    qp(igrd(n,j,ipass))=dd(n+a*j+d)
                  end do
                end do
              endif
              iy=me*(ne-ns+1)*(kl-kbotdav+1)
              b=me*(ne-ns+1)
              d=-me*(ns+kbotdav*(ne-ns+1))
              if(nud_uv>0)then
                call MPI_Recv(dd(1:iy),iy,MPI_REAL,hproc,itag,
     &                 MPI_COMM_WORLD,status,ierr)
                do k=kbotdav,kl
                  do j=ns,ne
                    do n=1,me
                      qu(igrd(n,j,ipass),k)=dd(n+a*j+b*k+d)
                    end do
                  end do
                end do
                call MPI_Recv(dd(1:iy),iy,MPI_REAL,hproc,itag,
     &                 MPI_COMM_WORLD,status,ierr)
                do k=kbotdav,kl
                  do j=ns,ne
                    do n=1,me
                      qv(igrd(n,j,ipass),k)=dd(n+a*j+b*k+d)
                    end do
                  end do
                end do
                call MPI_Recv(dd(1:iy),iy,MPI_REAL,hproc,itag,
     &                 MPI_COMM_WORLD,status,ierr)
                do k=kbotdav,kl
                  do j=ns,ne
                    do n=1,me
                      qw(igrd(n,j,ipass),k)=dd(n+a*j+b*k+d)
                    end do
                  end do
                end do
              end if
              if(nud_t>0)then
                call MPI_Recv(dd(1:iy),iy,MPI_REAL,hproc,itag,
     &                 MPI_COMM_WORLD,status,ierr)
                do k=kbotdav,kl
                  do j=ns,ne
                    do n=1,me
                      qt(igrd(n,j,ipass),k)=dd(n+a*j+b*k+d)
                    end do
                  end do
                end do
              end if
              if(nud_q>0)then
                call MPI_Recv(dd(1:iy),iy,MPI_REAL,hproc,itag,
     &                 MPI_COMM_WORLD,status,ierr)
                do k=kbotdav,kl
                  do j=ns,ne
                    do n=1,me
                      qq(igrd(n,j,ipass),k)=dd(n+a*j+b*k+d)
                    end do
                  end do
                end do
              end if
            end if
          end if

          if ((myid==0).and.(nmaxpr==1)) print *,"Start convolution"

          do j=ns,ne
            ma(1:me)=1./em_g(igrd(1:me,j,ipass))
            asum(1:me)=qsum(igrd(1:me,j,ipass))
            xa(1:me)=x_g(igrd(1:me,j,ipass))
            ya(1:me)=y_g(igrd(1:me,j,ipass))
            za(1:me)=z_g(igrd(1:me,j,ipass))
            if (nud_p>0) then
              ap(1:me)=qp(igrd(1:me,j,ipass))
            end if
            if (nud_uv>0) then
              au(1:me,:)=qu(igrd(1:me,j,ipass),:)
              av(1:me,:)=qv(igrd(1:me,j,ipass),:)
              aw(1:me,:)=qw(igrd(1:me,j,ipass),:)
            end if
            if (nud_t>0) then
              at(1:me,:)=qt(igrd(1:me,j,ipass),:)
            end if
            if (nud_q>0) then
              aq(1:me,:)=qq(igrd(1:me,j,ipass),:)
            end if
            do n=1,il_g
              ra(1:me)=xa(n)*xa(1:me)+ya(n)*ya(1:me)+za(n)*za(1:me)
              ra(1:me)=acos(max(min(ra(1:me),1.),-1.))
              ra(1:me)=exp(-(cq*ra(1:me))**2)*ma(1:me)
              ! can also use the lines below which integrate the gaussian
              ! analytically over the length element (but slower)
              !ra(1)=2.*erf(cq*0.5*(ds/rearth)*ma(1))
              !ra(2:me)=erf(cq*(ra(2:me)+0.5*(ds/rearth)*ma(2:me)))  ! redefine ra(:) as wgt(:)
     &        !        -erf(cq*(ra(2:me)-0.5*(ds/rearth)*ma(2:me))) ! (correct units are 1/cq)
              psum(n)=sum(ra(1:me)*asum(1:me))
              if (nud_p>0) then
                pp(n)=sum(ra(1:me)*ap(1:me))
              end if
              if (nud_uv>0) then
                do k=kbotdav,kl
                  pu(n,k)=sum(ra(1:me)*au(1:me,k))
                  pv(n,k)=sum(ra(1:me)*av(1:me,k))
                  pw(n,k)=sum(ra(1:me)*aw(1:me,k))
                end do
              end if
              if (nud_t>0) then
                do k=kbotdav,kl
                  pt(n,k)=sum(ra(1:me)*at(1:me,k))
                end do
              end if
              if (nud_q>0) then
                do k=kbotdav,kl
                  pq(n,k)=sum(ra(1:me)*aq(1:me,k))
                end do
              end if
            end do
            qsum(igrd(1:il_g,j,ipass))=psum(1:il_g)
            if (nud_p>0) then
              qp(igrd(1:il_g,j,ipass))=pp(1:il_g)
            end if
            if (nud_uv>0) then
              qu(igrd(1:il_g,j,ipass),:)=pu(1:il_g,:)
              qv(igrd(1:il_g,j,ipass),:)=pv(1:il_g,:)
              qw(igrd(1:il_g,j,ipass),:)=pw(1:il_g,:)
            end if
            if (nud_t>0) then
              qt(igrd(1:il_g,j,ipass),:)=pt(1:il_g,:)
            end if
            if (nud_q>0) then
              qq(igrd(1:il_g,j,ipass),:)=pq(1:il_g,:)
            end if
          end do

          if ((myid==0).and.(nmaxpr==1)) print *,"End convolution"

          if (ipass.eq.1) then
            pold=ppass          
            ssum(:)=qsum(:)
            if (nud_p>0) then
              sp(:)=qp(:)
            end if
            if (nud_uv>0) then
              su(:,:)=qu(:,:)
              sv(:,:)=qv(:,:)
              sw(:,:)=qw(:,:)
            end if
            if (nud_t>0) then
              st(:,:)=qt(:,:)
            end if
            if (nud_q>0) then
              sq(:,:)=qq(:,:)
            end if
          elseif ((ipass.eq.2).or.(ipass.eq.3)) then
            itag=itag+1
            if ((myid==0).and.(nmaxpr==1)) then
              print *,"Send arrays to local host"
            end if
            if (myid==hproc) then
              do iproc=hproc+1,mproc+hproc-1
                call procdiv(nns,nne,il_g,mproc,iproc-hproc)
                if (nns.gt.nne) exit
                iy=il_g*(nne-nns+1)*(kx(ipass)-kn(ipass)+1)
                a=il_g
                c=il_g*(nne-nns+1)
                d=-il_g*(nns+(nne-nns+1)*kn(ipass))
                call MPI_Recv(dd(1:iy),iy,MPI_REAL,iproc
     &                 ,itag,MPI_COMM_WORLD,status,ierr)
                do kpass=kn(ipass),kx(ipass)
                  do j=nns,nne
                    do n=1,il_g
                      qsum(igrd(n,j,kpass))=dd(n+a*j+c*kpass+d)
                    end do
                  end do
                end do
                if(nud_p>0)then
                  call MPI_Recv(dd(1:iy),iy,MPI_REAL,iproc,
     &                   itag,MPI_COMM_WORLD,status,ierr)
                  do kpass=kn(ipass),kx(ipass)
                    do j=nns,nne
                      do n=1,il_g
                        qp(igrd(n,j,kpass))=dd(n+a*j+c*kpass+d)
                      end do
                    end do
                  end do
                end if
                iy=il_g*(nne-nns+1)*(kl-kbotdav+1)
     &             *(kx(ipass)-kn(ipass)+1)
                b=il_g*(nne-nns+1)
                c=il_g*(nne-nns+1)*(kl-kbotdav+1)
                d=-il_g*(nns+(nne-nns+1)
     &            *(kbotdav+(kl-kbotdav+1)*kn(ipass)))
                if(nud_uv>0)then
                  call MPI_Recv(dd(1:iy),iy,MPI_REAL,iproc
     &                   ,itag,MPI_COMM_WORLD,status,ierr)
                  do k=kbotdav,kl
                    do kpass=kn(ipass),kx(ipass)
                      do j=nns,nne
                        do n=1,il_g
                          qu(igrd(n,j,kpass),k)
     &                      =dd(n+a*j+b*k+c*kpass+d)
                        end do
                      end do
                    end do
                  end do
                  call MPI_Recv(dd(1:iy),iy,MPI_REAL,iproc
     &                   ,itag,MPI_COMM_WORLD,status,ierr)
                  do k=kbotdav,kl
                    do kpass=kn(ipass),kx(ipass)
                      do j=nns,nne
                        do n=1,il_g
                          qv(igrd(n,j,kpass),k)
     &                      =dd(n+a*j+b*k+c*kpass+d)
                        end do
                      end do
                    end do
                  end do
                  call MPI_Recv(dd(1:iy),iy,MPI_REAL,iproc
     &                   ,itag,MPI_COMM_WORLD,status,ierr)
                  do k=kbotdav,kl
                    do kpass=kn(ipass),kx(ipass)
                      do j=nns,nne
                        do n=1,il_g
                          qw(igrd(n,j,kpass),k)
     &                      =dd(n+a*j+b*k+c*kpass+d)
                        end do
                      end do
                    end do
                  end do
                end if
                if(nud_t>0)then
                  call MPI_Recv(dd(1:iy),iy,MPI_REAL,iproc
     &                   ,itag,MPI_COMM_WORLD,status,ierr)
                  do k=kbotdav,kl
                    do kpass=kn(ipass),kx(ipass)
                      do j=nns,nne
                        do n=1,il_g
                          qt(igrd(n,j,kpass),k)
     &                      =dd(n+a*j+b*k+c*kpass+d)
                        end do
                      end do
                    end do
                  end do
                end if
                if(nud_q>0)then
                  call MPI_Recv(dd(1:iy),iy,MPI_REAL,iproc
     &                   ,itag,MPI_COMM_WORLD,status,ierr)
                  do k=kbotdav,kl
                    do kpass=kn(ipass),kx(ipass)
                      do j=nns,nne
                        do n=1,il_g
                          qq(igrd(n,j,kpass),k)
     &                      =dd(n+a*j+b*k+c*kpass+d)
                        end do
                      end do
                    end do
                  end do
                end if
              end do
            else
              iy=il_g*(ne-ns+1)*(kx(ipass)-kn(ipass)+1)
              a=il_g
              c=il_g*(ne-ns+1)
              d=-il_g*(ns+(ne-ns+1)*kn(ipass))
              do kpass=kn(ipass),kx(ipass)
                do j=ns,ne
                  do n=1,il_g
                    dd(n+a*j+c*kpass+d)=qsum(igrd(n,j,kpass))
                  end do
                end do
              end do
              call MPI_SSend(dd(1:iy),iy,MPI_REAL,hproc,itag,
     &               MPI_COMM_WORLD,ierr)
              if(nud_p>0)then
                do kpass=kn(ipass),kx(ipass)
                  do j=ns,ne
                    do n=1,il_g
                      dd(n+a*j+c*kpass+d)=qp(igrd(n,j,kpass))
                    end do
                  end do
                end do
                call MPI_SSend(dd(1:iy),iy,MPI_REAL,hproc,itag,
     &                 MPI_COMM_WORLD,ierr)
              end if
              iy=il_g*(ne-ns+1)*(kl-kbotdav+1)
     &           *(kx(ipass)-kn(ipass)+1)
              b=il_g*(ne-ns+1)
              c=il_g*(ne-ns+1)*(kl-kbotdav+1)
              d=-il_g*(ns+(ne-ns+1)
     &          *(kbotdav+(kl-kbotdav+1)*kn(ipass)))
              if(nud_uv>0)then
                do k=kbotdav,kl
                  do kpass=kn(ipass),kx(ipass)
                    do j=ns,ne
                      do n=1,il_g
                        dd(n+a*j+b*k+c*kpass+d)=qu(igrd(n,j,kpass),k)
                      end do
                    end do
                  end do
                end do
                call MPI_SSend(dd(1:iy),iy,MPI_REAL,hproc,itag,
     &                 MPI_COMM_WORLD,ierr)
                do k=kbotdav,kl
                  do kpass=kn(ipass),kx(ipass)
                    do j=ns,ne
                      do n=1,il_g
                        dd(n+a*j+b*k+c*kpass+d)=qv(igrd(n,j,kpass),k)
                      end do
                    end do
                  end do
                end do
                call MPI_SSend(dd(1:iy),iy,MPI_REAL,hproc,itag,
     &                 MPI_COMM_WORLD,ierr)
                do k=kbotdav,kl
                  do kpass=kn(ipass),kx(ipass)
                    do j=ns,ne
                      do n=1,il_g
                        dd(n+a*j+b*k+c*kpass+d)=qw(igrd(n,j,kpass),k)
                      end do
                    end do
                  end do
                end do
                call MPI_SSend(dd(1:iy),iy,MPI_REAL,hproc,itag,
     &                 MPI_COMM_WORLD,ierr)
              end if
              if(nud_t>0)then
                do k=kbotdav,kl
                  do kpass=kn(ipass),kx(ipass)
                    do j=ns,ne
                      do n=1,il_g
                        dd(n+a*j+b*k+c*kpass+d)=qt(igrd(n,j,kpass),k)
                      end do
                    end do
                  end do
                end do
                call MPI_SSend(dd(1:iy),iy,MPI_REAL,hproc,itag,
     &                 MPI_COMM_WORLD,ierr)
              end if
              if(nud_q>0)then
                do k=kbotdav,kl
                  do kpass=kn(ipass),kx(ipass)
                    do j=ns,ne
                      do n=1,il_g
                        dd(n+a*j+b*k+c*kpass+d)=qq(igrd(n,j,kpass),k)
                      end do
                    end do
                  end do
                end do
                call MPI_SSend(dd(1:iy),iy,MPI_REAL,hproc,itag,
     &                 MPI_COMM_WORLD,ierr)
              end if
            end if
          end if
          
        end do
      
      return  
      end subroutine speclocal
      !---------------------------------------------------------------------------------

      !---------------------------------------------------------------------------------
      subroutine getiqa(iq,ne,ipass,ppass,il_g)
      
      implicit none
      
      integer, intent(in) :: ne,ipass,ppass,il_g
      integer, dimension(ne,1:il_g), intent(out) :: iq
      integer sn,n,j,a,b,c
      
      do sn=1,ne,il_g

        select case(ppass*100+ipass*10+(sn-1)/il_g)
          case(0,310,530)                        ! panel 5   - x pass
            a=il_g
            b=-1
            c=5*il_g*il_g+1
          case(10,230,300)                       ! panel 2   - x pass
            a=1
            b=il_g
            c=il_g*(2*il_g-1)
          case(20,21)                            ! panel 0,1 - y pass
            a=il_g
            b=1
            c=-il_g
          case(22,432)                           ! panel 3   - y pass
            a=1
            b=-il_g
            c=il_g*(4*il_g-2)
          case(23)                               ! panel 4   - y pass
            a=1
            b=-il_g
            c=il_g*(5*il_g-3)
          case(30,100,410)                       ! panel 0   - z pass
            a=1
            b=il_g
            c=-il_g
          case(31)                               ! panel 2   - z pass
            a=il_g
            b=-1
            c=il_g*il_g+1
          case(32,222)                           ! panel 5   - z pass
            a=1
            b=il_g
            c=il_g*(5*il_g-3)
          case(110,231,330,400)                  ! panel 3   - z pass
            a=il_g
            b=-1
            c=3*il_g*il_g+1
          case(120)                              ! panel 1   - x pass
            a=1
            b=il_g
            c=il_g*(il_g-1)
          case(121)                              ! panel 2   - x pass
            a=1
            b=il_g
            c=2*il_g*(il_g-1)
          case(122,123,220,221)                  ! panel 4,5 - x pass ! panel 2,3 - z pass
            a=il_g
            b=-1
            c=2*il_g*il_g+1
          case(130,200,510)                      ! panel 1   - y pass
            a=il_g
            b=1
            c=il_g*(il_g-1)
          case(131)                              ! panel 3   - y pass
            a=1
            b=-il_g
            c=il_g*(4*il_g-1)
          case(132,322,323)                      ! panel 0,1 - y pass
            a=il_g
            b=1
            c=-il_g*(2*il_g+1)
          case(210,430,500)                      ! panel 4   - y pass
            a=1
            b=-il_g
            c=5*il_g*il_g
          case(223)                              ! panel 0   - z pass
            a=1
            b=il_g
            c=-4*il_g
          case(232,422)                          ! panel 1   - x pass
            a=1
            b=il_g
            c=il_g*(il_g-3)
          case(320)                              ! panel 3   - y pass
            a=1
            b=-il_g
            c=4*il_g*il_g
          case(321)                              ! panel 4   - y pass
            a=1
            b=-il_g
            c=il_g*(5*il_g-1)
          case(331)                              ! panel 5   - z pass
            a=1
            b=il_g
            c=il_g*(5*il_g-2)
          case(332,522,523)                      ! panel 2,3 - z pass 
            a=il_g
            b=-1
            c=1
          case(420,421)                          ! panel 4,5 - x pass
            a=il_g
            b=-1
            c=4*il_g*il_g+1
          case(423)                              ! panel 2   - x pass
            a=1
            b=il_g
            c=il_g*(2*il_g-4)
          case(431)                              ! panel 0   - y pass
            a=il_g
            b=1
            c=-il_g*(il_g+1)
          case(520)                              ! panel 5   - z pass
            a=1
            b=il_g
            c=il_g*(5*il_g-1)
          case(521)                              ! panel 0   - z pass
            a=1
            b=il_g
            c=-2*il_g
          case(531)                              ! panel 1   - x pass
            a=1
            b=il_g
            c=il_g*(il_g-2)
          case(532)                              ! panel 4   - x pass
            a=il_g
            b=-1
            c=2*il_g*il_g+1
          case DEFAULT
            print *,"Invalid index ",ppass,ipass,sn,
     &              ppass*100+ipass*10+(sn-1)/il_g
            exit
        end select
  
        do n=sn,sn+il_g-1
          do j=1,il_g
            iq(n,j)=a*n+b*j+c
          end do
        end do
  
      end do

      return
      end subroutine getiqa
      !---------------------------------------------------------------------------------

      !---------------------------------------------------------------------------------
      subroutine procdiv(ns,ne,ifull_g,nproc,myid)
      
      implicit none
      
      integer, intent(in) :: ifull_g,nproc,myid
      integer, intent(out) :: ns,ne
      integer npt,resid
      
      npt=ifull_g/nproc
      resid=mod(ifull_g,nproc)
      if ((myid+1).le.resid) then
        ns=myid*(npt+1)+1
        ne=(myid+1)*(npt+1)
      else
        ns=resid+myid*npt+1
        ne=resid+(myid+1)*npt
      end if
      
      return
      end subroutine procdiv
      !---------------------------------------------------------------------------------
