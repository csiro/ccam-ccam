      subroutine nestin                ! for globpea
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
      include 'map.h'
      include 'parm.h'     ! qgmin
      include 'pbl.h'      ! sice
      include 'sigs.h'
      include 'soil.h'     ! tss
      include 'soilsnow.h' ! tgg
      include 'stime.h'    ! kdate_s,ktime_s  sought values for data read
      common/nest/ta(ifull,kl),ua(ifull,kl),va(ifull,kl),psla(ifull),
     .            tb(ifull,kl),ub(ifull,kl),vb(ifull,kl),pslb(ifull),
     .            qa(ifull,kl),qb(ifull,kl),tssa(ifull),tssb(ifull)
      common/schmidtx/rlong0x,rlat0x,schmidtx ! infile, newin, nestin, indata
      common/sigin/sigin(kl),kk  ! for vertint, infile
      common/work2/zsb(ifull),dum1(ifull),tssi(ifull),dum2(ifull),
     .  dum3(ifull),dum4(ifull),dum5(ifull),dum6(ifull),dum7(ifull),
     .                        dumxx(ifull,18-9)
      common/work3/dum3a(ifull,kl),dum3b(ifull,kl,4)
      integer num,mtimea,mtimeb
      data num/0/,mtimea/0/,mtimeb/-1/
      save num,mtimea,mtimeb
!     mtimer, mtimeb are in minutes
      if(ktau<100.and.myid==0)then
        print *,'in nestin ktau,mtimer,mtimea,mtimeb ',
     &                                  ktau,mtimer,mtimea,mtimeb
        print *,'with kdate_s,ktime_s >= ',kdate_s,ktime_s
      end if
      if(mtimeb==-1)then
        mtimeb=mtimer  ! zero in fact
        if ( myid==0 )
     &  print *,'set nesting fields to those already read in via indata'
        do iq=1,ifull
         pslb(iq)=psl(iq)
         tssb(iq)=tss(iq)
        enddo
        tb(1:ifull,:)=t(1:ifull,:)
        qb(1:ifull,:)=qg(1:ifull,:)
        ub(1:ifull,:)=u(1:ifull,:)
        vb(1:ifull,:)=v(1:ifull,:)
        return
      endif       ! (mtimeb==-1)

      if(mtimer>0.and.mtimer<=mtimeb)go to 6  ! allows for dt<1 minute
!     transfer mtimeb fields to mtimea
      mtimea=mtimeb
      do iq=1,ifull
       psla(iq)=pslb(iq)
       tssa(iq)=tssb(iq)
      enddo
      ta(1:ifull,:)=tb(1:ifull,:)
      qa(1:ifull,:)=qb(1:ifull,:)
      ua(1:ifull,:)=ub(1:ifull,:)
      va(1:ifull,:)=vb(1:ifull,:)

!     Both infile and onthefly have an argument (nested) that says they've
!     been called from nestin and don't need to return all variables.
!     This would be cleaner than the use of all the dummy variables.

!     read tb etc  - for globpea, straight into tb etc
      if(io_in==1.or.io_in==3)then
      call infile(meso2,kdate_r,ktime_r,nem2, ! different from DARLAM
     .  timeg_b,ds_r,pslb,dum1,zsb,dum2,dum3, ! timeg_r,ds_r,psl,pmsl,zs,em,f
     .  tssb,dum4,dum3a,dum3a,dum5,dum6,dum7, ! tss,precip,wb,wbice,alb,snowd,sicedep
     .  tb,ub,vb,qb,dum3a,                    ! t,u,v,qg,tgg
     .	           dum3a,dum3a,dum3a,dum2, dum3,  dum4, dum5,  1)
!      above are: tggsn,smass,ssdn, ssdnn,osnowd,snage,isflag,nested
      endif   ! (io_in==1.or.io_in==3)

      if(io_in==-1.or.io_in==-3)then
          call onthefly(kdate_r,ktime_r,
     .     pslb,zsb,tssb,dum3a,dum3a,dum6,dum7,
     .     tb,ub,vb,qb,dum3a,                    ! t,u,v,qg,tgg
     .	           dum3a,dum3a,dum3a,dum2, dum3,  dum4, dum5,  1)
!      above are: tggsn,smass,ssdn, ssdnn,osnowd,snage,isflag,nested
      endif   ! (io_in==1.or.io_in==3)
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
!     mtimeb=60*24*(iabsdate(kdate_r,kdate)-iabsdate(kdate,kdate))
!    .               +nint(60*(ktime_r-ktime)/100.)  ! up till 29/11/02
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
      qb(1:ifull,:)=max(qb(1:ifull,:),qgmin)
      do k=kl-2,kl
       qb(1:ifull,k)=min(qb(1:ifull,k),10.*qgmin)
      enddo

      if(mod(ktau,nmaxpr)==0.or.ktau==2.or.diag)then
!       following is useful if troublesome data is read in
        if ( myid == 0 ) then
          print *,'following max/min values printed from nestin'
        end if
        call maxmin(ub,'ub',ktau,1.,kl)
        call maxmin(vb,'vb',ktau,1.,kl)
        call maxmin(tb,'tb',ktau,1.,kl)
        call maxmin(qb,'qb',ktau,1.e3,kl)
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

!     test code in nestin for modifying SSTs
!	do j=270,288
!	 do i=9,22
!	  tssb(i,j)=tssb(i,j)-.5
!	 enddo
!	enddo      

      if(newtop>=1)then
!       in these cases redefine pslb, tb and (effectively) zsb using zs
!       this keeps inner-mesh land mask & zs
!       presently simplest to do whole pslb, tb (& qb) arrays
        if(nmaxpr==1.and.mydiag)then
          print *,'zs (idjd) :',zs(idjd)
          print *,'zsb (idjd) :',zsb(idjd)
          print *,'psl (idjd) :',psl(idjd)
          print *,'pslb in(idjd) :',pslb(idjd)
          print *,'now call retopo from nestin'
        endif
        call retopo(pslb,zsb,zs,tb,qb)
        if(nmaxpr==1.and.mydiag)then
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

!     calculate time interpolated tss (into tssi)
      if(namip>0)return     ! namip SSTs/sea-ice take precedence
      tssi(:)=cona*tssa(:)+conb*tssb(:)  
      do iq=1,ifull
       if(.not.land(iq))then
	  if(tssi(iq)>273.2)then
	    tss(iq)=tssi(iq)
	    tgg(iq,1)=tssi(iq)
!          all others are for implied sea ice points
!          N.B. fracice, sice etc are updated once daily in sflux
	  elseif(tss(iq)>273.2)then  
	    tss(iq)=tssi(iq)
	    tgg(iq,1)=tssi(iq)
!          if already sea ice implied, use present tss, tgg1
         endif	    
       endif  ! (.not.land(iq))
      enddo  
	      
      return
      end
