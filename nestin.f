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
     .              pslb,zsb,tssb,sicedepb,fraciceb,tb,ub,vb,qb)
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

      subroutine nestinb(mins_mbd)  ! called for mbd>0 - spectral filter method
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
      integer ier,idv,mins_mbd
      include 'netcdf.inc'
      include 'mpif.h'
      include 'dates.h'    ! mtimer
    !  include 'dava.h'
    !  include 'davb.h'     ! psls,qgg,tt,uu,vv
      include 'indices.h'
      include 'latlong.h'
      include 'liqwpar.h'  ! ifullw,qfg,qlg
    !  include 'map.h'  
      include 'parm.h'     ! qgmin
      include 'pbl.h'      ! tss
      include 'sigs.h'
      include 'soil.h'     ! sicedep fracice
      include 'soilsnow.h' ! tgg
      include 'stime.h'    ! kdate_s,ktime_s  sought values for data read
    !  include 'vecsuv.h'   ! ax,bx etc
    !  include 'xyzinfo.h'
      common/nest/ta(ifull,kl),ua(ifull,kl),va(ifull,kl),psla(ifull),
     .            tb(ifull,kl),ub(ifull,kl),vb(ifull,kl),pslb(ifull),
     .            qa(ifull,kl),qb(ifull,kl),tssa(ifull),tssb(ifull),
     .            sicedepb(ifull),fraciceb(ifull)
      common/schmidtx/rlong0x,rlat0x,schmidtx ! infile, newin, nestin, indata
      real sigin
      integer ik,jk,kk
      common/sigin/ik,jk,kk,sigin(kl)  ! for vertint, infile
      integer num,mtimea,mtimeb,kdate_r,ktime_r
      integer ::  iabsdate,iq,k,kdhour,kdmin
      real :: ds_r,rlong0x,rlat0x
      real :: schmidtx,timeg_b,timerm
      real :: psla,pslb,qa,qb,ta,tb,tssa,tssb,ua,ub,va,vb
      real :: fraciceb,sicedepb
      real, dimension(ifull) ::  zsb
      data num/0/,mtimea/0/,mtimeb/-1/
      save num,mtimea,mtimeb
      if(num==0)then
!       this just retrieves time increment (in mins) on mesofile
        if(myid==0)then
         idv = ncvid(ncid,'mtimer',ier)
         if (ier.ne.0) idv = ncvid(ncid,'time',ier)
         call ncvgt1(ncid,idv,1,mtimea,ier) 
         call ncvgt1(ncid,idv,2,mtimeb,ier) 
         mins_mbd=mtimeb-mtimea
         print *,'nestinb: mtimea,mtimeb,myid ',mtimea,mtimeb,myid
         print *,'nestinb: ier,ncid,idv,mins_mbd ',ier,ncid,idv,mins_mbd
        endif
        print *,'nestinb: mtimea,mtimeb,myid ',mtimea,mtimeb,myid
        call MPI_Bcast(mins_mbd,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier) 
        if(mins_mbd<=0)stop
        num=1
        return
      endif   ! (num==0)
!     mtimer, mtimeb are in minutes
      if(ktau<100.and.myid==0)then
        print *,'in nestinb ktau,mtimer,mtimea,mtimeb,io_in ',
     &                      ktau,mtimer,mtimea,mtimeb,io_in
        print *,'with kdate_s,ktime_s >= ',kdate_s,ktime_s
      end if
      if(ntest==1)go to 6

!     read tb etc  - for globpea, straight into tb etc
      if(io_in==1)then
        call infil(1,kdate_r,ktime_r,timeg_b,ds_r, 
     .              pslb,zsb,tssb,sicedepb,fraciceb,tb,ub,vb,qb)
      endif   ! (io_in==1)

      if(io_in==-1)then
         call onthefl(1,kdate_r,ktime_r,
     &                 pslb,zsb,tssb,sicedepb,fraciceb,tb,ub,vb,qb) 
      endif   ! (io_in==1)
      tssb(:) = abs(tssb(:))  ! moved here Mar '03
      if (mydiag) then
        write (6,"('zsb# nestinb  ',9f7.1)") diagvals(zsb)
        write (6,"('tssb# nestinb ',9f7.1)") diagvals(tssb) 
      end if
   
      if(abs(rlong0  -rlong0x)>.01.or.
     &   abs(rlat0    -rlat0x)>.01.or.
     &   abs(schmidt-schmidtx)>.01)stop "grid mismatch in infile"

      kdhour=ktime_r/100-ktime/100   ! integer hour diff from Oct '05
      kdmin=(ktime_r-100*(ktime_r/100))-(ktime-100*(ktime/100))
      if ( myid == 0 ) then
        print *,'nestinb file has: kdate_r,ktime_r,kdhour,kdmin ',
     &                             kdate_r,ktime_r,kdhour,kdmin
      end if
      mtimeb=60*24*(iabsdate(kdate_r,kdate)-iabsdate(kdate,kdate))
     .               +60*kdhour+kdmin
      if ( myid == 0 ) then
        print *,'kdate_r,iabsdate ',kdate_r,iabsdate(kdate_r,kdate)
        print *,'giving mtimeb = ',mtimeb
!       print additional information
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
      endif
      if ( myid == 0 ) then
        print *,'following in nestinb after read pslb are psl not ps'
      end if
      call maxmin(pslb,'pB',ktau,100.,1)


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

6     timerm=ktau*dt/60.   ! real value in minutes (in case dt < 60 seconds)

!      if(num==1)then
!        num=2
!        call printa('zs  ',zs        ,ktau,0  ,ia,ib,ja,jb,0.,.01)
!        call printa('zsb ',zsb       ,ktau,0  ,ia,ib,ja,jb,0.,.01)
!        call printa('psl ',psl       ,ktau,0  ,ia,ib,ja,jb,0.,1.e2)
!        call printa('pslb',pslb      ,ktau,0  ,ia,ib,ja,jb,0.,1.e2)
!        call printa('t   ',t,ktau,nlv,ia,ib,ja,jb,200.,1.)
!        call printa('tb  ',tb,ktau,nlv,ia,ib,ja,jb,200.,1.)
!        call printa('u   ',u,ktau,nlv,ia,ib,ja,jb,0.,1.)
!        call printa('ub  ',ub,ktau,nlv,ia,ib,ja,jb,0.,1.)
!        call printa('v   ',v,ktau,nlv,ia,ib,ja,jb,0.,1.)
!        call printa('vb  ',vb,ktau,nlv,ia,ib,ja,jb,0.,1.)
!        call printa('davt',davt,0,0,ia,ib,ja,jb,0.,10.)
!      end if ! (num==1)

      call getspecdata(pslb,ub,vb,tb,qb)

!      print *,'following bunch near end of nestinb, nbd = ',nbd
!      write (6,"('100*psl.wesn ',2p5f8.3)") psl(idjd),
!     &           psl(iw(idjd)),psl(ie(idjd)),psl(is(idjd)),psl(in(idjd))
!      write (6,"('ps.wesn ',-2p5f9.3)") ps(idjd),
!     &           ps(iw(idjd)),ps(ie(idjd)),ps(is(idjd)),ps(in(idjd))
!      call maxmin(u,' u',ktau,1.,kl)
!      call maxmin(v,' v',ktau,1.,kl)
!      call maxmin(t,' t',ktau,1.,kl)
!      call maxmin(qg,'qg',ktau,1.e3,kl)
!      call maxmin(qfg,'qf',ktau,1.e3,kl)
!      call maxmin(qlg,'ql',ktau,1.e3,kl)
!      call maxmin(tggsn,'tggsn',ktau,1.,3)
!      call maxmin(tgg,'tg',ktau,1.,ms)
!      call maxmin(tss,'ts',ktau,1.,1)
!      call maxmin(ps,'ps',ktau,.01,1)

!     calculate time interpolated tss 
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
      include 'vecsuv_g.h'   ! ax_g,bx_g,ay_g,by_g,az_g,bz_g
      
      real, dimension(ifull), intent(in) :: pslb
      real, dimension(ifull,kl), intent(in) :: ub,vb,tb,qb
      real, dimension(ifull) :: delta
      real, dimension(ifull_g) :: x_g,xx_g,pslc
      real, dimension(ifull_g,kbotdav:kl) :: uc,vc,wc,tc,qc
      integer k
      
      if (myid == 0) then     
        print *,"Gather data for spectral downscale"
        call ccmpi_gather(pslb(:)-psl(1:ifull), pslc(:))
        do k=kbotdav,kl
          call ccmpi_gather(ub(1:ifull,k)-u(1:ifull,k), x_g(:))
          call ccmpi_gather(vb(1:ifull,k)-v(1:ifull,k), xx_g(:))
          uc(:,k)=ax_g(:)*x_g(:)+bx_g(:)*xx_g(:)
          vc(:,k)=ay_g(:)*x_g(:)+by_g(:)*xx_g(:)
          wc(:,k)=az_g(:)*x_g(:)+bz_g(:)*xx_g(:)
          call ccmpi_gather(tb(1:ifull,k)-t(1:ifull,k), tc(:,k))
          call ccmpi_gather(qb(1:ifull,k)-qg(1:ifull,k), qc(:,k))
        end do
      else
        call ccmpi_gather(pslb(:)-psl(1:ifull))
        do k=kbotdav,kl
          call ccmpi_gather(ub(1:ifull,k)-u(1:ifull,k))
          call ccmpi_gather(vb(1:ifull,k)-v(1:ifull,k))
          call ccmpi_gather(tb(1:ifull,k)-t(1:ifull,k))
          call ccmpi_gather(qb(1:ifull,k)-qg(1:ifull,k))
        end do
      end if
      
      !-----------------------------------------------------------------------
      select case(nud_uv) ! replace nud_uv with a new switch in the NML?
        case DEFAULT
          if (myid == 0) then
            print *,"Fast 1D spectral downscale"
            call fastspec((.1*real(mbd)/(pi*schmidt))**2
     &                    ,pslc,uc,vc,wc,tc,qc) ! e.g. mbd=40
          end if
        case(2)
          if (myid == 0) print *,"6-pass 1D spectral downscale"
          call sixspec(myid,.1*real(mbd)/(pi*schmidt)
     &                  ,pslc,uc,vc,wc,tc,qc)
        case(3)
          if (myid == 0) print *,"2D spectral downscale"
          call slowspec(myid,.1*real(mbd)/(pi*schmidt)
     &                  ,pslc,uc,vc,wc,tc,qc)
        case(4)
          if (myid == 0) print *,"Fast 1D spectral downscale (MPI)"
          call fastspecmpi(myid,.1*real(mbd)/(pi*schmidt)
     &                  ,pslc,uc,vc,wc,tc,qc)
      end select
      !-----------------------------------------------------------------------

      if (myid == 0) then
        print *,"Distribute data from spectral downscale"
        if (nud_p.gt.0) then
          call ccmpi_distribute(delta(:), pslc(:))
          psl(1:ifull)=psl(1:ifull)+delta(:)
        end if
        do k=kbotdav,kl
          if (nud_uv.gt.0) then
            x_g(:)=ax_g(:)*uc(:,k)+ay_g(:)*vc(:,k)+az_g(:)*wc(:,k)
            call ccmpi_distribute(delta(:), x_g(:))
            u(1:ifull,k)=u(1:ifull,k)+delta(:)
            xx_g(:)=bx_g(:)*uc(:,k)+by_g(:)*vc(:,k)+bz_g(:)*wc(:,k)
            call ccmpi_distribute(delta(:), xx_g(:))
            v(1:ifull,k)=v(1:ifull,k)+delta(:)
          end if
          if (nud_t.gt.0) then
            call ccmpi_distribute(delta(:), tc(:,k))
            t(1:ifull,k)=t(1:ifull,k)+delta(:)
          end if
          if (nud_q.gt.0) then
            call ccmpi_distribute(delta(:), qc(:,k))
            qg(1:ifull,k)=max(qg(1:ifull,k)+delta(:),qgmin)
          end if
        end do

      else
        if (nud_p.gt.0) then
          call ccmpi_distribute(delta(:))
          psl(1:ifull)=psl(1:ifull)+delta(:)
        end if
        do k=kbotdav,kl
          if (nud_uv.gt.0) then
            call ccmpi_distribute(delta(:))
            u(1:ifull,k)=u(1:ifull,k)+delta(:)
            call ccmpi_distribute(delta(:))
            v(1:ifull,k)=v(1:ifull,k)+delta(:)
          end if
          if (nud_t.gt.0) then
            call ccmpi_distribute(delta(:))
            t(1:ifull,k)=t(1:ifull,k)+delta(:)
          end if
          if (nud_q.gt.0) then
            call ccmpi_distribute(delta(:))
            qg(1:ifull,k)=max(qg(1:ifull,k)+delta(:),qgmin)
          end if
        end do
      end if
      
      ps(1:ifull)=1.e5*exp(psl(1:ifull)) ! Do not think this is needed, but kept it anyway - MJT
            
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
      real, dimension(ifull_g), save :: xx,yy,zz
      real, dimension(ifull_g,kbotdav:kl) :: uu,vv,ww,tt,qgg
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
        
        if (em_g(iq).gt.emmin) then
        
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
          psls(iq)=psla(iq)
          qgg(iq,kbotdav:kl)=qa(iq,kbotdav:kl)
          tt(iq,kbotdav:kl)=ta(iq,kbotdav:kl)
          uu(iq,kbotdav:kl)=ua(iq,kbotdav:kl)
          vv(iq,kbotdav:kl)=va(iq,kbotdav:kl)
          ww(iq,kbotdav:kl)=wa(iq,kbotdav:kl)
        end if
       enddo    ! n loop
      enddo     ! j loop      
      psls(1:ifull_g)=psls(1:ifull_g)/sumwt(1:ifull_g)
      do k=kbotdav,kl
       qgg(1:ifull_g,k)=qgg(1:ifull_g,k)/sumwt(1:ifull_g)
       tt(1:ifull_g,k)=tt(1:ifull_g,k)/sumwt(1:ifull_g)
       uu(1:ifull_g,k)=uu(1:ifull_g,k)/sumwt(1:ifull_g)
       vv(1:ifull_g,k)=vv(1:ifull_g,k)/sumwt(1:ifull_g)
       ww(1:ifull_g,k)=ww(1:ifull_g,k)/sumwt(1:ifull_g)
      enddo
      do j=1,il_g
       do n=1,4*il_g
        if(n<=il_g)iq=il_g*(il_g+j-1)+n                   ! panel 1
        if(n>il_g.and.n<=2*il_g)iq=il_g*(2*il_g+j-2)+n      ! panel 2
        if(n>2*il_g)iq=il_g*(2*il_g+n-1)+il_g+1-j           ! panel 4,5
!       if(j==il_g/2)write(6,"('iq,diff,psls,pslb,psl,pslnew ',i8,5f9.5)")
!    &    iq,psla(iq),psls(iq),pslb(iq),psl(iq),psl(iq)+psls(iq)
!       if(j==il_g/2)write(6,"('iqx,diff,tt,tb,t,tnew ',i8,5f9.4)")
!    &    iq,ta(iq,9),tt(iq,9),tb(iq,9),t(iq,9),t(iq,9)+tt(iq,9)
!       gnuplot: plot 'a_pslx' u 3 w l, 'a_pslx' u 4 w l
        psla(iq)=psls(iq)
        do k=kbotdav,kl
         ta(iq,k)=tt(iq,k)
         ua(iq,k)=uu(iq,k)
         va(iq,k)=vv(iq,k)
         wa(iq,k)=ww(iq,k)
         qa(iq,k)=qgg(iq,k)
        enddo
       enddo  ! n loop
      enddo   ! j loop    for x-filter
      
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
          psls(iq)=psla(iq)
          qgg(iq,kbotdav:kl)=qa(iq,kbotdav:kl)
          tt(iq,kbotdav:kl)=ta(iq,kbotdav:kl)
          uu(iq,kbotdav:kl)=ua(iq,kbotdav:kl)
          vv(iq,kbotdav:kl)=va(iq,kbotdav:kl)
          ww(iq,kbotdav:kl)=wa(iq,kbotdav:kl)
        end if
       enddo    ! n loop
      enddo     ! i loop
      psls(1:ifull_g)=psls(1:ifull_g)/sumwt(1:ifull_g)
      do k=kbotdav,kl
       qgg(1:ifull_g,k)=qgg(1:ifull_g,k)/sumwt(1:ifull_g)
       tt(1:ifull_g,k)=tt(1:ifull_g,k)/sumwt(1:ifull_g)
       uu(1:ifull_g,k)=uu(1:ifull_g,k)/sumwt(1:ifull_g)
       vv(1:ifull_g,k)=vv(1:ifull_g,k)/sumwt(1:ifull_g)
       ww(1:ifull_g,k)=ww(1:ifull_g,k)/sumwt(1:ifull_g)
      enddo
      do i=1,il_g
       do n=1,4*il_g
        if(n<=2*il_g)iq=il_g*(n-1)+i                      ! panel 0,1
        if(n>2*il_g.and.n<=3*il_g)iq=il_g*(4*il_g-i-2)+n      ! panel 3
        if(n>3*il_g)iq=il_g*(5*il_g-i-3)+n                  ! panel 4
!       if(i==il_g/2)write(6,"('iq,diff,psls,pslb,psl,pslnew ',i8,5f9.5)")
!    &    iq,psla(iq),psls(iq),pslb(iq),psl(iq),psl(iq)+psls(iq)
!       if(i==il_g/2)write(6,"('iqy,diff,tt,tb,t,tnew ',i8,5f9.4)")
!    &    iq,ta(iq,9),tt(iq,9),tb(iq,9),t(iq,9),t(iq,9)+tt(iq,9)
        psla(iq)=psls(iq)
        do k=kbotdav,kl
         ta(iq,k)=tt(iq,k)
         ua(iq,k)=uu(iq,k)
         va(iq,k)=vv(iq,k)
         wa(iq,k)=ww(iq,k)
         qa(iq,k)=qgg(iq,k)
        enddo
       enddo  ! n loop
      enddo   ! i loop    for y-filter

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
        
         if (em_g(iq).gt.emmin) then
        
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
           psls(iq)=psla(iq)
           qgg(iq,kbotdav:kl)=qa(iq,kbotdav:kl)
           tt(iq,kbotdav:kl)=ta(iq,kbotdav:kl)
           uu(iq,kbotdav:kl)=ua(iq,kbotdav:kl)
           vv(iq,kbotdav:kl)=va(iq,kbotdav:kl)
           ww(iq,kbotdav:kl)=wa(iq,kbotdav:kl)
         end if
        enddo    ! n loop
       enddo     ! j loop      
       psls(1:ifull_g)=psls(1:ifull_g)/sumwt(1:ifull_g)
       do k=kbotdav,kl
        qgg(1:ifull_g,k)=qgg(1:ifull_g,k)/sumwt(1:ifull_g)
        tt(1:ifull_g,k)=tt(1:ifull_g,k)/sumwt(1:ifull_g)
        uu(1:ifull_g,k)=uu(1:ifull_g,k)/sumwt(1:ifull_g)
        vv(1:ifull_g,k)=vv(1:ifull_g,k)/sumwt(1:ifull_g)
        ww(1:ifull_g,k)=ww(1:ifull_g,k)/sumwt(1:ifull_g)
       enddo
       do j=1,il_g
        do n=1,4*il_g
         if(n<=il_g)iq=il_g*(j-1)+n                     ! panel 0
         if(n>il_g.and.n<=3*il_g)iq=il_g*(il_g+n-1)+il_g+1-j  ! panel 2,3
         if(n>3*il_g)iq=il_g*(5*il_g+j-4)+n               ! panel 5
!        if(j==il_g/2)write(6,"('iq,diff,psls,pslb,psl,pslnew ',i8,5f9.5)")
!    &    iq,psla(iq),psls(iq),pslb(iq),psl(iq),psl(iq)+psls(iq)
!        if(j==il_g/2)write(6,"('iqz,diff,tt,tb,t,tnew ',i8,5f9.4)")
!    &    iq,ta(iq,9),tt(iq,9),tb(iq,9),t(iq,9),t(iq,9)+tt(iq,9)
         psla(iq)=psls(iq)
         do k=kbotdav,kl
          ta(iq,k)=tt(iq,k)
          ua(iq,k)=uu(iq,k)
          va(iq,k)=vv(iq,k)
          wa(iq,k)=ww(iq,k)
          qa(iq,k)=qgg(iq,k)
         enddo
        enddo  ! n loop
       enddo   ! j loop    for "z"-filter
      end if ! (mbd.ge.0)

      return
      end subroutine fastspec


      !---------------------------------------------------------------------------------
      ! Slow 2D spectral downscaling - MPI version (half as fast for single processor)
      subroutine slowspec(myid,c,psls,uu,vv,ww,tt,qgg)
      
      implicit none
      
      include 'newmpar.h'   ! ifull_g,kl
      include 'const_phys.h' ! rearth,pi,tpi
      include 'map_g.h'     ! em_g
      include 'parm.h'      ! ds,kbotdav
      include 'xyzinfo_g.h' ! x_g,y_g,z_g
      include 'mpif.h'

      integer, intent(in) :: myid
      real, intent(in) :: c
      real, dimension(ifull_g), intent(inout) :: psls
      real, dimension(ifull_g,kbotdav:kl), intent(inout) :: uu,vv,ww
      real, dimension(ifull_g,kbotdav:kl), intent(inout) :: tt,qgg
      real, dimension(ifull_g) :: pp,psum
      real, dimension(ifull_g,kbotdav:kl) :: pu,pv,pw,pt,pq
      real :: r,wgt,rmaxsq,csq,emmin
      integer :: iq,j,ns,ne,k,itag=0,ierr,iproc,npt
      integer, dimension(MPI_STATUS_SIZE) :: status

      emmin=c*ds/rearth
      rmaxsq=1./c**2
      csq=-4.5*c**2

      pp=0.
      pu=0.
      pv=0.
      pw=0.
      pt=0.
      pq=0.
      psum=0.

      if (myid == 0) then
        print *,"Send global arrays to all processors"
        do iproc=1,nproc-1
          call MPI_SSend(psls(:),ifull_g,MPI_REAL,iproc,itag,
     &           MPI_COMM_WORLD,ierr)    
          do k=kbotdav,kl
            call MPI_SSend(uu(:,k),ifull_g,MPI_REAL,iproc,itag,
     &             MPI_COMM_WORLD,ierr)    
            call MPI_SSend(vv(:,k),ifull_g,MPI_REAL,iproc,itag,
     &             MPI_COMM_WORLD,ierr)    
            call MPI_SSend(ww(:,k),ifull_g,MPI_REAL,iproc,itag,
     &             MPI_COMM_WORLD,ierr)    
            call MPI_SSend(tt(:,k),ifull_g,MPI_REAL,iproc,itag,
     &             MPI_COMM_WORLD,ierr)    
            call MPI_SSend(qgg(:,k),ifull_g,MPI_REAL,iproc,itag,
     &             MPI_COMM_WORLD,ierr)    
          end do
        end do
      else
        call MPI_Recv(psls(:),ifull_g,MPI_REAL,0,itag,
     &         MPI_COMM_WORLD,status,ierr)
        do k=kbotdav,kl
          call MPI_Recv(uu(:,k),ifull_g,MPI_REAL,0,itag,
     &           MPI_COMM_WORLD,status,ierr)
          call MPI_Recv(vv(:,k),ifull_g,MPI_REAL,0,itag,
     &           MPI_COMM_WORLD,status,ierr)
          call MPI_Recv(ww(:,k),ifull_g,MPI_REAL,0,itag,
     &           MPI_COMM_WORLD,status,ierr)
          call MPI_Recv(tt(:,k),ifull_g,MPI_REAL,0,itag,
     &           MPI_COMM_WORLD,status,ierr)
          call MPI_Recv(qgg(:,k),ifull_g,MPI_REAL,0,itag,
     &           MPI_COMM_WORLD,status,ierr)
        end do
      end if
    
      npt=ifull_g/nproc
      ns=myid*npt+1
      ne=(myid+1)*npt
      if ((myid+1).eq.nproc) ne=ifull_g
      
      do iq=ns,ne
        if (em_g(iq).gt.emmin) then
          do j=1,ifull_g
            r=x_g(iq)*x_g(j)+y_g(iq)*y_g(j)+z_g(iq)*z_g(j)
            r=min(1.,max(r,0.))
            r=acos(r)**2
            if (r.le.rmaxsq) then            
              wgt=exp(r*csq)/(em_g(j)**2) ! correct units are ((ds/rearth)/em_g)**2
              psum(iq)=psum(iq)+wgt
              pp(iq)=pp(iq)+wgt*psls(j)
              pu(iq,:)=pu(iq,:)+wgt*uu(j,:)
              pv(iq,:)=pv(iq,:)+wgt*vv(j,:)
              pw(iq,:)=pw(iq,:)+wgt*ww(j,:)
              pt(iq,:)=pt(iq,:)+wgt*tt(j,:)
              pq(iq,:)=pq(iq,:)+wgt*qgg(j,:)
            end if
          end do
        end if     
      end do

      do iq=ns,ne
        psls(iq)=pp(iq)/psum(iq)
        uu(iq,:)=pu(iq,:)/psum(iq)
        vv(iq,:)=pv(iq,:)/psum(iq)
        ww(iq,:)=pw(iq,:)/psum(iq)
        tt(iq,:)=pt(iq,:)/psum(iq)
        qgg(iq,:)=pq(iq,:)/psum(iq)
      end do

      itag=itag+1
      if (myid == 0) then
        print *,"Receive array sections from all processors"
        do iproc=1,nproc-1
          ns=iproc*npt+1
          ne=(iproc+1)*npt
          if ((iproc+1).eq.nproc) ne=ifull_g
          call MPI_Recv(psls(ns:ne),ne-ns+1,MPI_REAL,iproc,itag,
     &           MPI_COMM_WORLD,status,ierr)
          do k=kbotdav,kl
            call MPI_Recv(uu(ns:ne,k),ne-ns+1,MPI_REAL,iproc,itag,
     &             MPI_COMM_WORLD,status,ierr)
            call MPI_Recv(vv(ns:ne,k),ne-ns+1,MPI_REAL,iproc,itag,
     &             MPI_COMM_WORLD,status,ierr)
            call MPI_Recv(ww(ns:ne,k),ne-ns+1,MPI_REAL,iproc,itag,
     &             MPI_COMM_WORLD,status,ierr)
            call MPI_Recv(tt(ns:ne,k),ne-ns+1,MPI_REAL,iproc,itag,
     &             MPI_COMM_WORLD,status,ierr)
            call MPI_Recv(qgg(ns:ne,k),ne-ns+1,MPI_REAL,iproc,itag,
     &             MPI_COMM_WORLD,status,ierr)
          end do
        end do
      else
        call MPI_SSend(psls(ns:ne),ne-ns+1,MPI_REAL,0,itag,
     &         MPI_COMM_WORLD,ierr)
        do k=kbotdav,kl
          call MPI_SSend(uu(ns:ne,k),ne-ns+1,MPI_REAL,0,itag,
     &           MPI_COMM_WORLD,ierr)
          call MPI_SSend(vv(ns:ne,k),ne-ns+1,MPI_REAL,0,itag,
     &           MPI_COMM_WORLD,ierr)
          call MPI_SSend(ww(ns:ne,k),ne-ns+1,MPI_REAL,0,itag,
     &           MPI_COMM_WORLD,ierr)
          call MPI_SSend(tt(ns:ne,k),ne-ns+1,MPI_REAL,0,itag,
     &           MPI_COMM_WORLD,ierr)
          call MPI_SSend(qgg(ns:ne,k),ne-ns+1,MPI_REAL,0,itag,
     &           MPI_COMM_WORLD,ierr)
        end do
      end if

      return
      end subroutine slowspec
      !---------------------------------------------------------------------------------

      !---------------------------------------------------------------------------------
      ! Six pass spectral downscaling
      subroutine sixspec(myid,c,psls,uu,vv,ww,tt,qgg)
      
      implicit none
      
      include 'newmpar.h'   ! ifull_g,kl
      include 'const_phys.h' ! rearth,pi,tpi
      include 'map_g.h'     ! em_g
      include 'parm.h'      ! ds,kbotdav
      include 'xyzinfo_g.h' ! x_g,y_g,z_g
      include 'mpif.h'
      
      integer, intent(in) :: myid
      real, intent(in) :: c
      real, dimension(ifull_g), intent(inout) :: psls
      real, dimension(ifull_g,kbotdav:kl), intent(inout) :: uu,vv,ww
      real, dimension(ifull_g,kbotdav:kl), intent(inout) :: tt,qgg
      real, dimension(ifull_g,1:3) :: pp,qp,psum,qsum
      real, dimension(ifull_g,kbotdav:kl,1:3) :: pu,pv,pw,pt,pq
      real, dimension(ifull_g,kbotdav:kl,1:3) :: qu,qv,qw,qt,qq
      real :: r,wgtb,wgt,rmaxsq,csq,emmin
      integer :: iq,iq1,j,ipass,s,e,til,n1,n,nmax,p
      integer :: ne,ns,iproc,k,itag=0,ierr,npt
      integer, dimension(MPI_STATUS_SIZE) :: status
      integer, dimension(0:5,1:2) :: maps
 
      maps(:,1)=(/ 1, 3, 2, 2, 1, 1 /)
      maps(:,2)=(/ 0, 0, 0, 3, 2, 3 /)
 
      emmin=c*ds/rearth
      rmaxsq=1./c**2
      csq=-4.5*c**2
      til=il_g*il_g
      npt=il_g/nproc
      ns=myid*npt+1
      ne=(myid+1)*npt
      if ((myid+1).eq.nproc) ne=il_g
      
      do j=1,3
        qp(:,j)=psls(:)
        qu(:,:,j)=uu(:,:)
        qv(:,:,j)=vv(:,:)
        qw(:,:,j)=ww(:,:)
        qt(:,:,j)=tt(:,:)
        qq(:,:,j)=qgg(:,:)
        qsum(:,j)=1.
      end do

      do ipass=0,5
        if (myid == 0) print *,"6 pass ",ipass

        if (maps(ipass,2).eq.0) then
          e=1
        else
          e=2
        end if

        if (ipass.le.2) then
          nmax=2*il_g
        else
          nmax=4*il_g
        end if

        if (myid == 0) then
          do iproc=1,nproc-1
            ns=iproc*npt+1
            ne=(iproc+1)*npt
            if ((iproc+1).eq.nproc) ne=il_g
            n1=(ne-ns+1)*nmax
            do s=1,e
              p=maps(ipass,s)
              do j=ns,ne
                do n=1,nmax
                  call getiq(iq,j,n,ipass,il_g)
                  iq1=(j-ns)*nmax+n
                  psum(iq1,p)=qsum(iq,p)
                  pp(iq1,p)=qp(iq,p)
                  pu(iq1,:,p)=qu(iq,:,p)
                  pv(iq1,:,p)=qv(iq,:,p)
                  pw(iq1,:,p)=qw(iq,:,p)
                  pt(iq1,:,p)=qt(iq,:,p)
                  pq(iq1,:,p)=qq(iq,:,p)
                end do
              end do
              call MPI_SSend(psum(1:n1,p),n1,MPI_REAL,iproc,itag,
     &               MPI_COMM_WORLD,ierr)    
              call MPI_SSend(pp(1:n1,p),n1,MPI_REAL,iproc,itag,
     &               MPI_COMM_WORLD,ierr)    
              do k=kbotdav,kl
                call MPI_SSend(pu(1:n1,k,p),n1,MPI_REAL,iproc,itag,
     &                 MPI_COMM_WORLD,ierr)
                call MPI_SSend(pv(1:n1,k,p),n1,MPI_REAL,iproc,itag,
     &                 MPI_COMM_WORLD,ierr)
                call MPI_SSend(pw(1:n1,k,p),n1,MPI_REAL,iproc,itag,
     &                 MPI_COMM_WORLD,ierr)    
                call MPI_SSend(pt(1:n1,k,p),n1,MPI_REAL,iproc,itag,
     &                 MPI_COMM_WORLD,ierr)
                call MPI_SSend(pq(1:n1,k,p),n1,MPI_REAL,iproc,itag,
     &                 MPI_COMM_WORLD,ierr)    
              end do
            end do
          end do
          ns=1
          ne=npt
        else
          n1=(ne-ns+1)*nmax
          do s=1,e
            p=maps(ipass,s)
            call MPI_Recv(psum(1:n1,p),n1,MPI_REAL,0,itag,
     &             MPI_COMM_WORLD,status,ierr)
            call MPI_Recv(pp(1:n1,p),n1,MPI_REAL,0,itag,
     &             MPI_COMM_WORLD,status,ierr)
            do k=kbotdav,kl
              call MPI_Recv(pu(1:n1,k,p),n1,MPI_REAL,0,itag,
     &               MPI_COMM_WORLD,status,ierr)
              call MPI_Recv(pv(1:n1,k,p),n1,MPI_REAL,0,itag,
     &               MPI_COMM_WORLD,status,ierr)
              call MPI_Recv(pw(1:n1,k,p),n1,MPI_REAL,0,itag,
     &               MPI_COMM_WORLD,status,ierr)
              call MPI_Recv(pt(1:n1,k,p),n1,MPI_REAL,0,itag,
     &               MPI_COMM_WORLD,status,ierr)
              call MPI_Recv(pq(1:n1,k,p),n1,MPI_REAL,0,itag,
     &               MPI_COMM_WORLD,status,ierr)
            end do
            do j=ns,ne
              do n=1,nmax
                call getiq(iq,j,n,ipass,il_g)
                iq1=(j-ns)*nmax+n
                qsum(iq,p)=psum(iq1,p)
                qp(iq,p)=pp(iq1,p)
                qu(iq,:,p)=pu(iq1,:,p)
                qv(iq,:,p)=pv(iq1,:,p)
                qw(iq,:,p)=pw(iq1,:,p)
                qt(iq,:,p)=pt(iq1,:,p)
                qq(iq,:,p)=pq(iq1,:,p)
              end do
            end do
          end do        
        end if

        pp=0.
        pu=0.
        pv=0.
        pw=0.
        pt=0.
        pq=0.
        psum=0.
     
        do j=ns,ne
          do n=1,nmax
            call getiq(iq,j,n,ipass,il_g)
            if (em_g(iq).gt.emmin) then
              do n1=n,nmax
                call getiq(iq1,j,n1,ipass,il_g)
                r=x_g(iq)*x_g(iq1)+y_g(iq)*y_g(iq1)+z_g(iq)*z_g(iq1)
                r=min(1.,max(r,0.))
                r=acos(r)**2
                if (r.le.rmaxsq) then
                  if (iq.eq.iq1) then
                    wgtb=1.
                  else
                    wgtb=exp(r*csq)	! lambda_min = 3 sigmas => -4.5
                    wgt=wgtb/em_g(iq) ! wrong units but the weights are rescaled so that sum(wgt)=1
                    do s=1,e
                      p=maps(ipass,s)
                      psum(iq1,p)=psum(iq1,p)+wgt*qsum(iq,p)
                      pp(iq1,p)=pp(iq1,p)+wgt*qp(iq,p)
                      pu(iq1,:,p)=pu(iq1,:,p)+wgt*qu(iq,:,p)
                      pv(iq1,:,p)=pv(iq1,:,p)+wgt*qv(iq,:,p)
                      pw(iq1,:,p)=pw(iq1,:,p)+wgt*qw(iq,:,p)
                      pt(iq1,:,p)=pt(iq1,:,p)+wgt*qt(iq,:,p)
                      pq(iq1,:,p)=pq(iq1,:,p)+wgt*qq(iq,:,p)
                    end do
                  end if
                  wgt=wgtb/em_g(iq1) ! correct units are ((ds/rearth)/em_g)**2
                  do s=1,e
                    p=maps(ipass,s)
                    psum(iq,p)=psum(iq,p)+wgt*qsum(iq1,p)
                    pp(iq,p)=pp(iq,p)+wgt*qp(iq1,p)
                    pu(iq,:,p)=pu(iq,:,p)+wgt*qu(iq1,:,p)
                    pv(iq,:,p)=pv(iq,:,p)+wgt*qv(iq1,:,p)
                    pw(iq,:,p)=pw(iq,:,p)+wgt*qw(iq1,:,p)
                    pt(iq,:,p)=pt(iq,:,p)+wgt*qt(iq1,:,p)
                    pq(iq,:,p)=pq(iq,:,p)+wgt*qq(iq1,:,p)
                  end do
                end if
              end do
              do s=1,e
                p=maps(ipass,s)
                qsum(iq,p)=psum(iq,p)
                qp(iq,p)=pp(iq,p)
                qu(iq,:,p)=pu(iq,:,p)
                qv(iq,:,p)=pv(iq,:,p)
                qw(iq,:,p)=pw(iq,:,p)
                qt(iq,:,p)=pt(iq,:,p)
                qq(iq,:,p)=pq(iq,:,p)
              end do
            end if
          end do
        end do

        itag=itag+1
        if (myid == 0) then
          do iproc=1,nproc-1
            ns=iproc*npt+1
            ne=(iproc+1)*npt
            if ((iproc+1).eq.nproc) ne=il_g
            n1=(ne-ns+1)*nmax
            do s=1,e
              p=maps(ipass,s)
              call MPI_Recv(psum(1:n1,p),n1,MPI_REAL,iproc
     &               ,itag,MPI_COMM_WORLD,status,ierr)
              call MPI_Recv(pp(1:n1,p),n1,MPI_REAL,iproc
     &               ,itag,MPI_COMM_WORLD,status,ierr)
              do k=kbotdav,kl
                call MPI_Recv(pu(1:n1,k,p),n1,MPI_REAL,iproc
     &                 ,itag,MPI_COMM_WORLD,status,ierr)
                call MPI_Recv(pv(1:n1,k,p),n1,MPI_REAL,iproc
     &                 ,itag,MPI_COMM_WORLD,status,ierr)
                call MPI_Recv(pw(1:n1,k,p),n1,MPI_REAL,iproc
     &                 ,itag,MPI_COMM_WORLD,status,ierr)
                call MPI_Recv(pt(1:n1,k,p),n1,MPI_REAL,iproc
     &                 ,itag,MPI_COMM_WORLD,status,ierr)
                call MPI_Recv(pq(1:n1,k,p),n1,MPI_REAL,iproc
     &                 ,itag,MPI_COMM_WORLD,status,ierr)
              end do
              do j=ns,ne
                do n=1,nmax
                  call getiq(iq,j,n,ipass,il_g)
                  iq1=(j-ns)*nmax+n
                  qsum(iq,p)=psum(iq1,p)
                  qp(iq,p)=pp(iq1,p)
                  qu(iq,:,p)=pu(iq1,:,p)
                  qv(iq,:,p)=pv(iq1,:,p)
                  qw(iq,:,p)=pw(iq1,:,p)
                  qt(iq,:,p)=pt(iq1,:,p)
                  qq(iq,:,p)=pq(iq1,:,p)
                end do
              end do
            end do
          end do
        else
          n1=(ne-ns+1)*nmax
          do s=1,e
            p=maps(ipass,s)
            do j=ns,ne
              do n=1,nmax
                call getiq(iq,j,n,ipass,il_g)
                iq1=(j-ns)*nmax+n
                psum(iq1,p)=qsum(iq,p)
                pp(iq1,p)=qp(iq,p)
                pu(iq1,:,p)=qu(iq,:,p)
                pv(iq1,:,p)=qv(iq,:,p)
                pw(iq1,:,p)=qw(iq,:,p)
                pt(iq1,:,p)=qt(iq,:,p)
                pq(iq1,:,p)=qq(iq,:,p)
              end do
            end do
            call MPI_SSend(psum(1:n1,p),n1,MPI_REAL,0,itag,
     &             MPI_COMM_WORLD,ierr)
            call MPI_SSend(pp(1:n1,p),n1,MPI_REAL,0,itag,
     &             MPI_COMM_WORLD,ierr)
            do k=kbotdav,kl
              call MPI_SSend(pu(1:n1,k,p),n1,MPI_REAL,0,itag,
     &               MPI_COMM_WORLD,ierr)
              call MPI_SSend(pv(1:n1,k,p),n1,MPI_REAL,0,itag,
     &               MPI_COMM_WORLD,ierr)
              call MPI_SSend(pw(1:n1,k,p),n1,MPI_REAL,0,itag,
     &               MPI_COMM_WORLD,ierr)
              call MPI_SSend(pt(1:n1,k,p),n1,MPI_REAL,0,itag,
     &               MPI_COMM_WORLD,ierr)
              call MPI_SSend(pq(1:n1,k,p),n1,MPI_REAL,0,itag,
     &               MPI_COMM_WORLD,ierr)
            end do
          end do
        end if
        
      end do
      
      ! update panels
      do ipass=0,5
        s=ipass*til+1
        e=(ipass+1)*til
        p=mod(ipass,3)+1
        psls(s:e)=qp(s:e,p)/qsum(s:e,p)
        do k=kbotdav,kl
          uu(s:e,k)=qu(s:e,k,p)/qsum(s:e,p)
          vv(s:e,k)=qv(s:e,k,p)/qsum(s:e,p)
          ww(s:e,k)=qw(s:e,k,p)/qsum(s:e,p)
          tt(s:e,k)=qt(s:e,k,p)/qsum(s:e,p)
          qgg(s:e,k)=qq(s:e,k,p)/qsum(s:e,p)
        end do
      end do
      
      return
      end subroutine sixspec
      !---------------------------------------------------------------------------------

      !---------------------------------------------------------------------------------
      subroutine fastspecmpi(myid,c,psls,uu,vv,ww,tt,qgg)
      
      implicit none
      
      include 'newmpar.h'   ! ifull_g,kl
      include 'const_phys.h' ! rearth,pi,tpi
      include 'map_g.h'     ! em_g
      include 'parm.h'      ! ds,kbotdav
      include 'xyzinfo_g.h' ! x_g,y_g,z_g
      include 'mpif.h'
      
      integer, intent(in) :: myid
      real, intent(in) :: c
      real, dimension(ifull_g), intent(inout) :: psls
      real, dimension(ifull_g,kbotdav:kl), intent(inout) :: uu,vv,ww
      real, dimension(ifull_g,kbotdav:kl), intent(inout) :: tt,qgg
      real, dimension(ifull_g) :: psum,pp
      real, dimension(ifull_g,kbotdav:kl) :: pu,pv,pw,pt,pq
      real :: r,wgtb,wgt,rmaxsq,csq,emmin
      integer :: iq,iq1,j,ipass,n1,n
      integer :: ne,ns,iproc,k,itag = 0,ierr,npt
      integer, dimension(MPI_STATUS_SIZE) :: status
 
      emmin=c*ds/rearth
      rmaxsq=1./c**2
      csq=-4.5*c**2

      npt=il_g/nproc
      ns=myid*npt+1
      ne=(myid+1)*npt
      if ((myid+1).eq.nproc) ne=il_g
                 
      do ipass=3,5

        if (myid == 0) then
          do iproc=1,nproc-1
            ns=iproc*npt+1
            ne=(iproc+1)*npt
            if ((iproc+1).eq.nproc) ne=il_g
            n1=(ne-ns+1)*4*il_g
            do j=ns,ne
              do n=1,4*il_g
                call getiq(iq,j,n,ipass,il_g)
                iq1=(j-ns)*4*il_g+n
                pp(iq1)=psls(iq)
                pu(iq1,:)=uu(iq,:)
                pv(iq1,:)=vv(iq,:)
                pw(iq1,:)=ww(iq,:)
                pt(iq1,:)=tt(iq,:)
                pq(iq1,:)=qgg(iq,:)
              end do
            end do
            call MPI_SSend(pp(1:n1),n1,MPI_REAL,iproc,itag,
     &             MPI_COMM_WORLD,ierr)    
            do k=kbotdav,kl
              call MPI_SSend(pu(1:n1,k),n1,MPI_REAL,iproc,itag,
     &               MPI_COMM_WORLD,ierr)
              call MPI_SSend(pv(1:n1,k),n1,MPI_REAL,iproc,itag,
     &               MPI_COMM_WORLD,ierr)
              call MPI_SSend(pw(1:n1,k),n1,MPI_REAL,iproc,itag,
     &               MPI_COMM_WORLD,ierr)    
              call MPI_SSend(pt(1:n1,k),n1,MPI_REAL,iproc,itag,
     &               MPI_COMM_WORLD,ierr)
              call MPI_SSend(pq(1:n1,k),n1,MPI_REAL,iproc,itag,
     &               MPI_COMM_WORLD,ierr)    
            end do
          end do
          ns=1
          ne=npt
        else
          n1=(ne-ns+1)*4*il_g
          call MPI_Recv(pp(1:n1),n1,MPI_REAL,0,itag,
     &           MPI_COMM_WORLD,status,ierr)
          do k=kbotdav,kl
            call MPI_Recv(pu(1:n1,k),n1,MPI_REAL,0,itag,
     &             MPI_COMM_WORLD,status,ierr)
            call MPI_Recv(pv(1:n1,k),n1,MPI_REAL,0,itag,
     &             MPI_COMM_WORLD,status,ierr)
            call MPI_Recv(pw(1:n1,k),n1,MPI_REAL,0,itag,
     &             MPI_COMM_WORLD,status,ierr)
            call MPI_Recv(pt(1:n1,k),n1,MPI_REAL,0,itag,
     &             MPI_COMM_WORLD,status,ierr)
            call MPI_Recv(pq(1:n1,k),n1,MPI_REAL,0,itag,
     &             MPI_COMM_WORLD,status,ierr)
          end do
          do j=ns,ne
            do n=1,4*il_g
              call getiq(iq,j,n,ipass,il_g)
              iq1=(j-ns)*4*il_g+n
              psls(iq)=pp(iq1)
              uu(iq,:)=pu(iq1,:)
              vv(iq,:)=pv(iq1,:)
              ww(iq,:)=pw(iq1,:)
              tt(iq,:)=pt(iq1,:)
              qgg(iq,:)=pq(iq1,:)
            end do
          end do
        end if

        pp=0.
        pu=0.
        pv=0.
        pw=0.
        pt=0.
        pq=0.
        psum=0.
     
     
        do j=ns,ne
          do n=1,4*il_g
            call getiq(iq,j,n,ipass,il_g)
            if (em_g(iq).gt.emmin) then
              do n1=n,4*il_g
                call getiq(iq1,j,n1,ipass,il_g)
                r=x_g(iq)*x_g(iq1)+y_g(iq)*y_g(iq1)+z_g(iq)*z_g(iq1)
                r=min(1.,max(r,0.))
                r=acos(r)**2
                if (r.le.rmaxsq) then
                  if (iq.eq.iq1) then
                    wgtb=1.
                  else
                    wgtb=exp(r*csq)	! lambda_min = 3 sigmas => -4.5
                    wgt=wgtb/em_g(iq) ! wrong units but the weights are rescaled so that sum(wgt)=1
                    psum(iq1)=psum(iq1)+wgt
                    pp(iq1)=pp(iq1)+wgt*psls(iq)
                    pu(iq1,:)=pu(iq1,:)+wgt*uu(iq,:)
                    pv(iq1,:)=pv(iq1,:)+wgt*vv(iq,:)
                    pw(iq1,:)=pw(iq1,:)+wgt*ww(iq,:)
                    pt(iq1,:)=pt(iq1,:)+wgt*tt(iq,:)
                    pq(iq1,:)=pq(iq1,:)+wgt*qgg(iq,:)
                  end if
                  wgt=wgtb/em_g(iq1) ! correct units are ((ds/rearth)/em_g)**2
                  psum(iq)=psum(iq)+wgt
                  pp(iq)=pp(iq)+wgt*psls(iq1)
                  pu(iq,:)=pu(iq,:)+wgt*uu(iq1,:)
                  pv(iq,:)=pv(iq,:)+wgt*vv(iq1,:)
                  pw(iq,:)=pw(iq,:)+wgt*ww(iq1,:)
                  pt(iq,:)=pt(iq,:)+wgt*tt(iq1,:)
                  pq(iq,:)=pq(iq,:)+wgt*qgg(iq1,:)
                end if
              end do
              psls(iq)=pp(iq)/psum(iq)
              uu(iq,:)=pu(iq,:)/psum(iq)
              vv(iq,:)=pv(iq,:)/psum(iq)
              ww(iq,:)=pw(iq,:)/psum(iq)
              tt(iq,:)=pt(iq,:)/psum(iq)
              qgg(iq,:)=pq(iq,:)/psum(iq)
            end if
          end do
        end do

        itag=itag+1
        if (myid == 0) then
          do iproc=1,nproc-1
            ns=iproc*npt+1
            ne=(iproc+1)*npt
            if ((iproc+1).eq.nproc) ne=il_g
            n1=(ne-ns+1)*4*il_g
            call MPI_Recv(pp(1:n1),n1,MPI_REAL,iproc
     &             ,itag,MPI_COMM_WORLD,status,ierr)
            do k=kbotdav,kl
              call MPI_Recv(pu(1:n1,k),n1,MPI_REAL,iproc
     &               ,itag,MPI_COMM_WORLD,status,ierr)
              call MPI_Recv(pv(1:n1,k),n1,MPI_REAL,iproc
     &               ,itag,MPI_COMM_WORLD,status,ierr)
              call MPI_Recv(pw(1:n1,k),n1,MPI_REAL,iproc
     &               ,itag,MPI_COMM_WORLD,status,ierr)
              call MPI_Recv(pt(1:n1,k),n1,MPI_REAL,iproc
     &               ,itag,MPI_COMM_WORLD,status,ierr)
              call MPI_Recv(pq(1:n1,k),n1,MPI_REAL,iproc
     &               ,itag,MPI_COMM_WORLD,status,ierr)
            end do
            do j=ns,ne
              do n=1,4*il_g
                call getiq(iq,j,n,ipass,il_g)
                iq1=(j-ns)*4*il_g+n
                psls(iq)=pp(iq1)
                uu(iq,:)=pu(iq1,:)
                vv(iq,:)=pv(iq1,:)
                ww(iq,:)=pw(iq1,:)
                tt(iq,:)=pt(iq1,:)
                qgg(iq,:)=pq(iq1,:)
              end do
            end do
          end do
        else
          n1=(ne-ns+1)*4*il_g
          do j=ns,ne
            do n=1,4*il_g
              call getiq(iq,j,n,ipass,il_g)
              iq1=(j-ns)*4*il_g+n
              pp(iq1)=psls(iq)
              pu(iq1,:)=uu(iq,:)
              pv(iq1,:)=vv(iq,:)
              pw(iq1,:)=ww(iq,:)
              pt(iq1,:)=tt(iq,:)
              pq(iq1,:)=qgg(iq,:)
            end do
          end do
          call MPI_SSend(pp(1:n1),n1,MPI_REAL,0,itag,
     &           MPI_COMM_WORLD,ierr)
          do k=kbotdav,kl
            call MPI_SSend(pu(1:n1,k),n1,MPI_REAL,0,itag,
     &             MPI_COMM_WORLD,ierr)
            call MPI_SSend(pv(1:n1,k),n1,MPI_REAL,0,itag,
     &             MPI_COMM_WORLD,ierr)
            call MPI_SSend(pw(1:n1,k),n1,MPI_REAL,0,itag,
     &             MPI_COMM_WORLD,ierr)
            call MPI_SSend(pt(1:n1,k),n1,MPI_REAL,0,itag,
     &             MPI_COMM_WORLD,ierr)
            call MPI_SSend(pq(1:n1,k),n1,MPI_REAL,0,itag,
     &             MPI_COMM_WORLD,ierr)
          end do
        end if
        
      end do
      
      return
      end subroutine fastspecmpi
      !---------------------------------------------------------------------------------
      
      !---------------------------------------------------------------------------------
      subroutine getiq(iq,j,n,ipass,il_g)
      
      implicit none
      
      integer, intent(out) :: iq
      integer, intent(in) :: j,n,ipass,il_g
      
      select case(ipass) ! Use indexing from jlm fastspec
        case(0)                         ! broken x pass
          if (n<=il_g) then
            iq=il_g*(2*il_g+j-1)+n      ! panel 2
          else
            iq=il_g*(4*il_g+n)+1-j      ! panel 5
          end if
        case(1)                         ! broken y pass
          if (n<=il_g) then
            iq=il_g*(n-1)+j             ! panel 0
          else
            iq=il_g*(4*il_g-j-1)+n      ! panel 3
          end if
        case(2)                         ! broken z pass
          if (n<=il_g) then
            iq=il_g*(j-1)+n             ! panel 0
          else
            iq=il_g*(2*il_g+n)+1-j      ! panel 3
          end if
        case(3)                         ! x pass
          if (n<=il_g) then
            iq=il_g*(il_g+j-1)+n        ! panel 1
          else if (n<=2*il_g) then
            iq=il_g*(2*il_g+j-2)+n      ! panel 2
          else
            iq=il_g*(2*il_g+n)+1-j      ! panel 4,5
          end if
        case(4)                         ! y pass
          if (n<=2*il_g) then
            iq=il_g*(n-1)+j             ! panel 0,1
          else if (n<=3*il_g) then
            iq=il_g*(4*il_g-j-2)+n      ! panel 3
          else 
            iq=il_g*(5*il_g-j-3)+n      ! panel 4
          end if
        case DEFAULT                    ! z pass
          if (n<=il_g) then
            iq=il_g*(j-1)+n             ! panel 0
          else if (n<=3*il_g) then
            iq=il_g*(il_g+n)+1-j        ! panel 2,3
          else
            iq=il_g*(5*il_g+j-4)+n      ! panel 5
          end if
      end select      
      
      end subroutine getiq
      !---------------------------------------------------------------------------------
