      subroutine infile(io_in2,kdate_r,ktime_r,nem2,
     . timeg_r,ds_r,psl,pmsl,zs,em,f,
     . tss,precip,wb,wbice,alb,snowd,sicedep,
     . t,u,v,qg,tgg,
     . tggsn,smass,ssdn, ssdnn,osnowd,snage,isflag,nested)
!     note kk; vertint.f is attached below
!     kdate_r and ktime_r are returned from this routine.
!     They must equal or exceed kdate_s and ktime_s
!     Note that kdate_s and ktime_s are updated progressively in this 
!     routine to allow for nesting
!     nested=0  for calls from indata; 1  for calls from nestin     
!     From March 2000 exact kdate, ktime are calculated and compared

!     This is called from nestin and onthefly with various dummy variables
!     because some fields aren't required.
!     Would it be better to do this with optional arguments?

      use cc_mpi
      use diag_m
      implicit none
      include 'newmpar.h'
      include 'darcdf.h'
      include 'kuocom.h'    ! ldr
      include 'liqwpar.h'  ! ifullw
      include 'netcdf.inc'
      include 'parm.h'
      include 'parm_nqg.h'  ! nqg_r,nqg_set
      include 'stime.h'     ! kdate_s,ktime_s  sought values for data read
      include 'tracers.h'
      include 'mpif.h'


      integer io_in2, kdate_r, ktime_r, nem2, nested
      real timeg_r, ds_r
      real psl(ifull),pmsl(ifull),zs(ifull),em(ifull),f(ifull),
     . tss(ifull),precip(ifull),wb(ifull,ms),wbice(ifull,ms),
     . snowd(ifull),alb(ifull),sicedep(ifull),
     . t(ifull,kl),u(ifull,kl),v(ifull,kl),qg(ifull,kl),
     . tgg(ifull,ms),tggsn(ifull,3),smass(ifull,3),ssdn(ifull,3),
     . ssdnn(ifull),osnowd(ifull),snage(ifull)
      integer isflag(ifull)
      integer ktau_r, ibb, jbb, i
      real timer_r, difknbd, du, tanl

      integer ini,inj,ink,m2,nsd2,mesi,nbd2
     & ,nps2,mex2,mup2,nemi,mtimer,nmi2,ndt2,npsav2,nhor2,nkuo2,khdif2
     & ,kwt2,iaa2,jaa2,nvad2,nqgi,lbd2,nrun2,nrunx2
     & ,khor2,ksc2,kountr2,ndiur2,nhort2,nhorps2,nsoil2,ivirt2
     & ,ntsuri,nrad2,kuocb2,nvmix2,ntsea2,ms2,nextras2,ilt2,ntrac2
      common/iahead/ini,inj,ink,m2,nsd2,mesi,nbd2
     & ,nps2,mex2,mup2,nemi,mtimer,nmi2,ndt2,npsav2,nhor2,nkuo2,khdif2
     & ,kwt2,iaa2,jaa2,nvad2,nqgi,lbd2,nrun2,nrunx2
     & ,khor2,ksc2,kountr2,ndiur2,nhort2,nhorps2,nsoil2,ivirt2
     & ,ntsuri,nrad2,kuocb2,nvmix2,ntsea2,ms2,nextras2,ilt2,ntrac2
      real tds,difknbdi,rhkuo,tdu,ttanl,trnml,tstl1,tstl2
      common/rhead/tds,difknbdi,rhkuo,tdu,ttanl,trnml,tstl1,tstl2
      real rlong0x,rlat0x,schmidtx ! infile, newin, nestin, indata
      common/schmidtx/rlong0x,rlat0x,schmidtx ! infile, newin, nestin, indata
      real sigin
      integer kk
      common/sigin/sigin(kl),kk  ! for vertint, infile
      real tmp(ifull)

      integer, parameter :: nihead=54,nrhead=14
      integer nahead(nihead)
      real ahead(nrhead)
      equivalence (nahead,ini),(ahead,tds)

      integer timest(2),timeco(2),timerco(2)
      character*80 dimnam
      character*3 monn(12),cmon
      data monn/'jan','feb','mar','apr','may','jun'
     &         ,'jul','aug','sep','oct','nov','dec'/

      integer, save :: ncidold=-1, ncalled=0, iarch, nqg_setin
      integer itype, ilen, ier, ierr, ik, jk, k, ndim, nvars, ngatts,
     &        nulid, nd, isiz, ix, iy, narch, idy, imo, iyr, ihr,
     &        mtimer_in, iq, ii, jj, idv
      real dss, xgn, ygn, timer

      if ( abs(io_in) /= 1 ) then
         print*,
     &    "Error: only abs(io_in)=1 (netcdf) input longer supported."
         call MPI_Abort(MPI_COMM_WORLD)
      endif

      if(ncalled.eq.0)then
        nqg_setin=nqg_set
      endif
      ncalled=ncalled+1
      ncid=idifil
      if (mydiag) print *,'newin ncid,ncidold=',ncid,ncidold
      if ( ncid.ne.ncidold ) iarch=1
      ncidold=ncid
      
c save model map projs.
      dss=ds_r

      if(dss.gt.0..and.npanels.eq.0)then   ! for later DARLAM option
         print*,
     &   "Error, npanels=0 DARLAM option not implemented in MPI version"
         stop
      endif ! (dss.gt.0..and.npanels.eq.0)then

      if ( myid == 0 ) then
         call ncainq(ncid,ncglobal,'int_header',itype,ilen,ier)
         call ncagt(ncid,ncglobal,'int_header',nahead,ier)
         call ncainq(ncid,ncglobal,'real_header',itype,ilen,ier)
         call ncagt(ncid,ncglobal,'real_header',ahead,ier)
      end if
      call MPI_Bcast(nahead,nihead,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(ahead,nrhead,MPI_REAL,0,MPI_COMM_WORLD,ierr)

      io_in2=1
      nem2=nemi
      nqg_r=nqgi
      nqg_set=min(nqg_r,nqg_setin)
      ibb=iaa2
      jbb=jaa2
      ntsur2=ntsuri
      rlong0x =ahead(5)
      rlat0x  =ahead(6)
      schmidtx=ahead(7)
      if(mydiag) print *,'newin  rlong0x,rlat0x,schmidtx ',
     &                    rlong0x,rlat0x,schmidtx
      if(schmidtx.le.0..or.schmidtx.gt.1.)then
!      read it where Jack put it
       rlong0x =ahead(6)
       rlat0x  =ahead(7)
       schmidtx=ahead(8)
       if(mydiag)print *,'newin  rlong0x,rlat0x,schmidtx ',
     &                    rlong0x,rlat0x,schmidtx
      endif      
      if(ktau.le.1)then
         if ( mydiag ) then
            print *,'nahead ',nahead
            print *,'ahead ',ahead
            write(6,'("io_in2,nem2,ntsur2,nqg_r,nqg_set,ibb,jbb",7i6)')
     &           io_in2,nem2,ntsur2,nqg_r,nqg_set,ibb,jbb
         end if
        if(ahead(7).eq.0.)then
	   if (mydiag)
     &      print *,'************** schmidt etc not in netcdf header!!!'
!         assume very old file 	   
           rlong0x =0.
           rlat0x  =0.
           schmidtx=1.
        endif  ! (ahead(7).eq.0.)then
      endif    ! (ktau.le.1)
      
      ik=nahead(1)
      jk=nahead(2)
      kk=nahead(3)
      if(ktau.le.1.and.myid==0)
     &     print *,'in newin; ktau,ik,jk,kk=', ktau,ik,jk,kk

c     turn OFF fatal netcdf errors
      if (myid==0) then
         call ncpopt(0)
         call ncagt(ncid,ncglobal,'sigma',sigin,ier)
         if ( ier.ne.0 ) then
            call ncagt(ncid,ncglobal,'sigma_lev',sigin,ier)
            print *,'ier=',ier
            if ( ier.ne.0 ) then
               idv = ncvid(ncid,'lev',ier)
               if(ier.ne.0)idv = ncvid(ncid,'layer',ier)
               call ncvgt(ncid,idv,1,kk,sigin,ier)
               print *,'ier=',ier
            endif
         endif
         if(ktau.le.1)write(6,'("sigin=",(9f7.4))') (sigin(k),k=1,kk)
      end if
      call MPI_Bcast(sigin,kk,MPI_REAL,0,MPI_COMM_WORLD,ierr)

c     Get dimensions
      if ( myid == 0 ) then
         call ncinq(ncid,ndim,nvars,ngatts,nulid,ier)
         if(ktau.le.1)print *,'ncid,ndim,nvars,ngatts,nulid,ier'
     &                     ,ncid,ndim,nvars,ngatts,nulid,ier

         do nd=1,ndim
            call ncdinq(ncid,nd,dimnam,isiz,ier)
            if(ier.ne.0)write(6,*) 'ncdinq dim1 ier=',ier
            if(ktau.le.1)write(6,*) 'dim ',nd,' ',isiz,' ',dimnam
            if(nd.eq.1) ix=isiz
            if(nd.eq.2) iy=isiz
            if(dimnam.eq.'time')narch=isiz
         enddo

         if(ktau.le.1)write(6,'("ix,iy,narch=",3i6)') ix,iy,narch
      end if
      call MPI_Bcast(ix,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(iy,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(narch,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      timest(1) = 1
      timeco(1) = 3
      timeco(2) = 1
      timerco(1) = 2
      timerco(2) = 1

c***********************************************************************
      if ( myid == 0 ) then
        iarch=iarch-1
19      iarch=iarch+1
        timest(2) = iarch
c       turn OFF fatal netcdf errors
        call ncpopt(0)
        idv = ncvid(ncid,'timestamp',ier)
        print *,'timestamp idv,ier,iarch=',idv,ier,iarch
c       turn ON fatal netcdf errors
        call ncpopt(NCVERBOS+NCFATAL)
        if ( ier.eq.0 ) then
           timeco(1)=1
           timeco(2)=iarch
           call ncvgt1(ncid,idv,timeco,kdate_r,ier)
           timeco(1)=2
           timeco(2)=iarch
           call ncvgt1(ncid,idv,timeco,ktime_r,ier)
           timeco(1)=3
           timeco(2)=iarch
           call ncvgt1(ncid,idv,timeco,ktau_r,ier)
           idv = ncvid(ncid,'realstamp',ier)
           timeco(1)=1
           timeco(2)=iarch
           call ncvgt1(ncid,idv,timeco,timer,ier)
!         N.B. usual ncdump does not show full precision of floating variables
           timeco(1)=2
           timeco(2)=iarch
           call ncvgt1(ncid,idv,timeco,timeg_r,ier)
        else                    ! ier.ne.0
!     new style
           idv = ncvid(ncid,'kdate',ier)
           call ncvgt1(ncid,idv,iarch,kdate_r,ier)
!     turn off fatal ncdf errors
           call ncpopt(0)
           idv = ncvid(ncid,'timer',ier)
           if ( ier.eq.0 ) then
              call ncvgt1(ncid,idv,iarch,timer,ier)
           endif                ! ( ier.eq.0 ) then
           idv = ncvid(ncid,'mtimer',ier)
           if ( ier.eq.0 ) then
              call ncvgt1(ncid,idv,iarch,mtimer,ier)
           endif                ! ( ier.eq.0 ) then
!         turn on fatal ncdf errors
           call ncpopt(NCVERBOS+NCFATAL)
           idv = ncvid(ncid,'ktime',ier)
           call ncvgt1(ncid,idv,iarch,ktime_r,ier)
           idv = ncvid(ncid,'ktau',ier)
           call ncvgt1(ncid,idv,iarch,ktau_r,ier)
           idv = ncvid(ncid,'timeg',ier)
           call ncvgt1(ncid,idv,iarch,timeg_r,ier)
        endif                   ! ier = 0
        if(kdate_r.eq.0)then
           write(6,*)'kdate = 0; create kdate from time_origin'
           idv = ncvid(ncid,'time',ier)
!         write(6,*)'time idv=',idv,ier
           call ncpopt(0)
           call ncagtc(ncid,idv,'time_origin',dimnam,80,ier)
!         write(6,*)'time_origin=',dimnam,lens,ier
           call ncpopt(NCVERBOS+NCFATAL)
           read(dimnam,'(i2)') idy
           read(dimnam,'(3x,a3)') cmon
           do imo=1,12
              if(cmon.eq.monn(imo))go to 10
           enddo
 10        continue
           read(dimnam,'(7x,i4)') iyr
           read(dimnam,'(12x,i2)') ihr
           kdate_r=idy+100*(imo+100*iyr)
           ktime_r=ihr*100
        endif                   ! (kdate_r.eq.0)
        print *,'kdate,ktime,ktau,timer_in: ',kdate_r,ktime_r,ktau,timer
!       fix up timer in 140-year run for roundoff (jlm Mon  01-11-1999)
        if(timer.gt.1048583.)then ! applies year 440 mid-month 9
           timer=timer+.125
        elseif(timer.gt.1048247.)then ! applies year 440 month 9
           timer=timer+.0625
        endif
        print *,'dtin,ktau_r: ',dtin,ktau_r
        print *,'kdate,ktime: ',kdate_r,ktime_r,
     .       'kdate_s,ktime_s >= ',kdate_s,ktime_s

        print *,'in newin; ktau_r,timer: ',ktau_r,timer
        print *,'values read in for timer, mtimer: ',timer,mtimer
        if(mtimer.ne.0)then     ! preferred
!     assume mtimer is read in correctly and set timer from that
           timer=mtimer/60.
        else                    ! mtimer = 0
           mtimer=nint(60.*timer)
        endif
        print *,'giving timer, mtimer: ',timer,mtimer
        mtimer_in=mtimer
        call datefix(kdate_r,ktime_r,mtimer) ! for Y2K, or mtimer>0
        if (2400*kdate_r+ktime_r.lt.2400*kdate_s+ktime_s)go to 19
        
        print *,'found suitable date/time in newin'
        print *,'in ',ik,jk,kk,m2,nsd2,io_in2,nbd2
     & ,nps2,mex2,mup2,nem2,mtimer,nmi2,ndt2,npsav2,nhor2,nkuo2,khdif2
     & ,kwt2,iaa2,jaa2,nvad2,nqg_r,lbd2,nrun2,nrunx2
     & ,khor2,ksc2,kountr2,ndiur2,nhort2,nhorps2,nsoil2,ivirt2
     & ,ntsur2,nrad2,kuocb2,nvmix2,ntsea2,ms2,nextras2,ilt2,ntrac2

        print *,'newin ds_r,nqg_r,mtimer: ',ds_r,nqg_r,mtimer
        if(mtimer_in.gt.0.and.nrungcm.eq.2)then
c         nrungcm=0
           print *,
     &          '*** re-setting NRUNGCM nrungcm from 2 to 0 because ',
     &             'not starting from gcm file'
           stop 'wrong nrungcm used'
        endif
      end if ! myid == 0
      ! Now share the time variables (is this actually necessary?)
      call MPI_Bcast(kdate_r,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(ktime_r,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(ktau_r,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(mtimer,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(timer,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(timeg_r,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)

c     begin reading data
c     log scaled sfc.press
      call histrd1(ncid,iarch,ier,'psf',ik,jk,psl)
c     call histrd1(ncid,iarch,ier,'mslp',ik,jk,pmsl)  ! not needed
      call histrd1(ncid,iarch,ier,'zht',ik,jk,zs)
      call histrd1(ncid,iarch,ier,'tsu',ik,jk,tss)

c     turn on fatal netcdf errors
      if ( myid == 0 ) call ncpopt(NCVERBOS+NCFATAL)
c     temperature
      call histrd4(ncid,iarch,ier,'temp',ik,jk,kk,t)
c     u wind component
      call histrd4(ncid,iarch,ier,'u',ik,jk,kk,u)
c     v wind component
      call histrd4(ncid,iarch,ier,'v',ik,jk,kk,v)
c     mixing ratio
c     turn OFF fatal netcdf errors; from here on
      if ( myid == 0 ) call ncpopt(0)
      call histrd4(ncid,iarch,ier,'q',ik,jk,kk,qg)       !  kg/kg
      if(ier.ne.0)then
        call histrd4(ncid,iarch,ier,'mixr',ik,jk,kk,qg)  !   g/kg
      endif
      if(ldr.ne.0)then
        call histrd4(ncid,iarch,ier,'qfg',ik,jk,kk,qfg(1:ifullw,:))       !  kg/kg
        call histrd4(ncid,iarch,ier,'qlg',ik,jk,kk,qlg(1:ifullw,:))       !  kg/kg
        if(ier.ne.0)then
!         set default qfg, qlg to zero if not available on input file	
          qfg=0. 
          qlg=0. 
        endif
      endif     ! (ldr.ne.0)

      if ( myid == 0 ) print *,'in newin nested = ',nested
      if(nested.eq.0)then   !  following only at start of run #############
!       read fields which may be there on restart, but maybe not initially
        call histrd1(ncid,iarch,ier,'alb',ik,jk,alb)
        call histrd1(ncid,iarch,ier,'tgg2',ik,jk,tgg(1,2))
        if(ier.eq.0)then  ! at least tgg6, wb2, wb6 will be available
          call histrd1(ncid,iarch,ier,'tgg6',ik,jk,tgg(1,6))
          call histrd1(ncid,iarch,ier,'wb2',ik,jk,wb(1,2))
          call histrd1(ncid,iarch,ier,'wb6',ik,jk,wb(1,6))
          call histrd1(ncid,iarch,ier,'tgg1',ik,jk,tgg(1,1)) ! trial read
          if(ier.eq.0)then
            call histrd1(ncid,iarch,ier,'tgg3',ik,jk,tgg(1,3))
            call histrd1(ncid,iarch,ier,'tgg4',ik,jk,tgg(1,4))
            call histrd1(ncid,iarch,ier,'tgg5',ik,jk,tgg(1,5))
            call histrd1(ncid,iarch,ier,'wb1',ik,jk,wb(1,1))
            call histrd1(ncid,iarch,ier,'wb3',ik,jk,wb(1,3))
            call histrd1(ncid,iarch,ier,'wb4',ik,jk,wb(1,4))
            call histrd1(ncid,iarch,ier,'wb5',ik,jk,wb(1,5))
          else  ! set other levels having read tgg2, tgg6, wb2, wb6
           do iq=1,ifull
            tgg(iq,1)=tss(iq)   ! initial temperature at second layer(6.5.97 KF)
            tgg(iq,3)=tgg(iq,2) ! initial temper. from GCM runs with 3 layers
            wb(iq,1) =wb(iq,2)  ! layer initialisation of moisture
           enddo  ! iq
           do k=3,ms-1          ! don't want to change value of tgg(,2) and tgg(,ms)
            do iq=1,ifull
             wb(iq,k) =wb(iq,ms)   ! layer initialisation of moisture
            enddo  ! iq
           enddo   ! k
           do k=4,ms-1
            do iq=1,ifull
             tgg(iq,k)=tgg(iq,ms) ! initial temper. from GCM runs with 3 layers
            enddo  ! iq
           enddo   ! k
          endif
        else
           call histrd1(ncid,iarch,ier,'tb2',ik,jk,tgg(1,ms))
           call histrd1(ncid,iarch,ier,'tb3',ik,jk,tgg(1,2))
           call histrd1(ncid,iarch,ier,'wfg',ik,jk,wb(1,2))
           call histrd1(ncid,iarch,ier,'wfb',ik,jk,wb(1,ms))
           do iq=1,ifull
            tgg(iq,1)=tss(iq)   ! initial temperature at second layer(6.5.97 KF)
            tgg(iq,3)=tgg(iq,2) ! initial temper. from GCM runs with 3 layers
            wb(iq,1) =wb(iq,2)  ! layer initialisation of moisture
           enddo  ! iq
           do k=3,ms-1          ! don't want to change value of tgg(,2) and tgg(,ms)
            do iq=1,ifull
             wb(iq,k) =wb(iq,ms)   ! layer initialisation of moisture
            enddo  ! iq
           enddo   ! k
           do k=4,ms-1
            do iq=1,ifull
             tgg(iq,k)=tgg(iq,ms) ! initial temper. from GCM runs with 3 layers
            enddo  ! iq
           enddo   ! k
        endif    ! (ier.eq.0)
        if ( myid == 0 ) print *,'about to read snowd'
        call histrd1(ncid,iarch,ier,'snd',ik,jk,snowd)
        if(ier.ne.0.or.nqg_set.lt.6)then  ! preset snowd here if not avail.
!         when no snowd available initially, e.g. COMPARE III (jlm formula)
          do iq=1,ifull                             ! "do" missing till 5/6/01
!           for this one,snow lin. increases from 5 to 55, for T 270 down to 260
            if(abs(tss(iq)).lt.270.)snowd(iq)=
     .                                min(55.,5.*(271.-abs(tss(iq))))   
          enddo
          if(ncalled.lt.4.and.mydiag)then
           print *,'setting snowd in newin, ier,nqg_set= '  ,ier,nqg_set
           print *,'snowd# preset to: ', diagvals(snowd)
!     .            ((snowd(ii+(jj-1)*il),ii=id-1,id+1),jj=jd-1,jd+1)
          endif
        endif  ! (ier.ne.0.or.nqg_set.lt.6)
        call histrd1(ncid,iarch,ier,'smass1',ik,jk,smass(1,1))
        if(ier.ne.0)then  ! for smass1
          if(myid==0)
     &          print *,'setting smass,wbice etc in newin, ier= '  ,ier
 	  do iq=1,ifull
          smass(iq,1)=0.
          smass(iq,2)=0.
	   smass(iq,3)=0.
	   if(snowd(iq).gt.100.)then
            ssdn(iq,1)=240.
	   else
	     ssdn(iq,1) = 140.
	   endif		! (snowd(iq).gt.100.)
          do k=2,3
	    ssdn(iq,k)=ssdn(iq,1)
          enddo
	   ssdnn(iq)  = ssdn(iq,1)
	   isflag(iq) = 0
	   snage(iq)  = 0.
	   if(snowd(iq).gt.0.)tgg(iq,1)=min(tgg(iq,1),270.1)
          enddo   ! iq loop
          do iq=1,ifull
           do k=1,ms
!           following linearly from 0 to .99 for tgg=tfrz down to tfrz-5
            wbice(iq,k)=
     .             min(.99,max(0.,.99*(273.1-tgg(iq,k))/5.))*wb(iq,k) ! jlm
	    enddo ! ms
          enddo  ! iq loop
        else    ! assume all these variables available in restart file
         call histrd1(ncid,iarch,ier,'smass2',ik,jk,smass(1,2))
         call histrd1(ncid,iarch,ier,'smass3',ik,jk,smass(1,3))
         call histrd1(ncid,iarch,ier,'ssdn1',ik,jk,ssdn(1,1))
         call histrd1(ncid,iarch,ier,'ssdn2',ik,jk,ssdn(1,2))
         call histrd1(ncid,iarch,ier,'ssdn3',ik,jk,ssdn(1,3))
         call histrd1(ncid,iarch,ier,'tggsn1',ik,jk,tggsn(1,1))
         call histrd1(ncid,iarch,ier,'tggsn2',ik,jk,tggsn(1,2))
         call histrd1(ncid,iarch,ier,'tggsn3',ik,jk,tggsn(1,3))
         call histrd1(ncid,iarch,ier,'wbice1',ik,jk,wbice(1,1))
         call histrd1(ncid,iarch,ier,'wbice2',ik,jk,wbice(1,2))
         call histrd1(ncid,iarch,ier,'wbice3',ik,jk,wbice(1,3))
         call histrd1(ncid,iarch,ier,'wbice4',ik,jk,wbice(1,4))
         call histrd1(ncid,iarch,ier,'wbice5',ik,jk,wbice(1,5))
         call histrd1(ncid,iarch,ier,'wbice6',ik,jk,wbice(1,6))
         call histrd1(ncid,iarch,ier,'snage',ik,jk,snage)
         call histrd1(ncid,iarch,ier,'sflag',ik,jk,tmp)
         do iq=1,ifull
          isflag(iq)=nint(tmp(iq))
         enddo
        endif  ! (ier.ne.0) ... else ...   for smass1
      
        call histrd1(ncid,iarch,ier,'siced',ik,jk,sicedep)  ! presets to follow
        if(ier.ne.0)then  ! for siced
 	   if(myid==0) print *,'pre-setting siced in newin'
           where ( tss <= 271.2 ) 
              sicedep = 2. 
           elsewhere
              sicedep = 0.     
           endwhere
        endif

         if ( mydiag ) then
            print *,'at end of newin kdate,ktime,ktau,tss: ',
     &                     kdate_r,ktime_r,ktau,tss(idjd)
            print *,'tgg ',(tgg(idjd,k),k=1,ms)
            print *,'wb ',(wb(idjd,k),k=1,ms)
            print *,'wbice ',(wbice(idjd,k),k=1,ms)
         end if
      endif  ! (nested.eq.0)   !  only done at start of run 

      iarch=iarch+1

      if ( mydiag ) then
         write(6,'("end newin kdate,ktime,ktau_r,mtimer,timer,ds_r="
     &      ,4i10,2f10.1)')kdate_r,ktime_r,ktau_r,mtimer,timer,ds_r
      end if

      if(ncalled.lt.4.and.mydiag)then
          print *,'sig in: ',(sigin(i),i=1,kk)
          write (6,"('100*psl# in',3f7.2,1x,3f7.2,1x,3f7.2)") 
     &              100.*diagvals(psl)
!     .              ((100.*psl(ii+(jj-1)*il),ii=id-1,id+1),jj=jd-1,jd+1)
          write (6,"('zs# in  ',3f7.1,1x,3f7.1,1x,3f7.1)") 
     &              diagvals(zs)
!     .              ((zs(ii+(jj-1)*il),ii=id-1,id+1),jj=jd-1,jd+1)
          print *,'em# in: ', diagvals(em)
!     . 	   ((em(ii+(jj-1)*il),ii=id-1,id+1),jj=jd-1,jd+1)
          print *,'f# in: ', diagvals(f)
!     . 	   ((f(ii+(jj-1)*il),ii=id-1,id+1),jj=jd-1,jd+1)
          write (6,"('tss# in ',3f7.1,1x,3f7.1,1x,3f7.1)") 
     &          diagvals(tss)
!     .              ((tss(ii+(jj-1)*il),ii=id-1,id+1),jj=jd-1,jd+1)
          write (6,"('prec# in',3f7.2,1x,3f7.2,1x,3f7.2)") 
     &          diagvals(precip)
!     .           ((precip(ii+(jj-1)*il),ii=id-1,id+1),jj=jd-1,jd+1)
!  Printing the ifull value gives confusing results in the parallel version
          print *,'t in ',(t(idjd,k),k=1,kk)!,t(ifull,kk)
          print *,'u in ',(u(idjd,k),k=1,kk)!,u(ifull,kk)
          print *,'v in ',(v(idjd,k),k=1,kk)!,v(ifull,kk)
          print *,'qg in ',(qg(idjd,k),k=1,kk)!,qg(ifull,kk)
          print *,'N.B. following are meaningless during nestin read'
          print *,'tgg(1)# ', diagvals(tgg(:,1))
!     .           ((tgg(ii+(jj-1)*il,1),ii=id-1,id+1),jj=jd-1,jd+1)
          print *,'tgg(2)# ', diagvals(tgg(:,2))
!     .           ((tgg(ii+(jj-1)*il,2),ii=id-1,id+1),jj=jd-1,jd+1)
          print *,'tgg(ms)# ', diagvals(tgg(:,ms))
!     .           ((tgg(ii+(jj-1)*il,ms),ii=id-1,id+1),jj=jd-1,jd+1)
          print *,'wb(1)# ', diagvals(wb(:,1))
!     .           ((wb(ii+(jj-1)*il,1),ii=id-1,id+1),jj=jd-1,jd+1)
          print *,'wb(ms)# ', diagvals(wb(:,ms))
!     .           ((wb(ii+(jj-1)*il,ms),ii=id-1,id+1),jj=jd-1,jd+1)
          print *,'alb# in: ', diagvals(alb)
!     . 	   ((alb(ii+(jj-1)*il),ii=id-1,id+1),jj=jd-1,jd+1)
          print *,'snowd# in: ', diagvals(snowd)
!     .           ((snowd(ii+(jj-1)*il),ii=id-1,id+1),jj=jd-1,jd+1)
       endif                    ! (ncalled.lt.4)

!     reset kdate_s and ktime_s ready for next (nesting) call to infile
      kdate_s=kdate_r
      ktime_s=ktime_r+1 

      qg(1:ifull,:) = max(qg(1:ifull,:),1.e-6)

      if ( mydiag ) then
         print *,'end infile; next read will be kdate_s,ktime_s >= ',
     &                                       kdate_s,ktime_s
      end if

      end subroutine infile

c***************************************************************************
      subroutine histrd1(idnc,iarch,ier,name,ik,jk,var)
      use cc_mpi
      implicit none
      include 'newmpar.h'
      include 'parm.h'
      include 'netcdf.inc'
      include 'mpif.h'
      integer idnc, iarch, ier, ik, jk
      integer*2 ivar(ik*jk)
      logical odiag
      parameter ( odiag=.false. )
      character name*(*)
      integer start(3),count(3)
      real var(ifull)
      real globvar(ifull_g), vmax, vmin, addoff, sf
      integer ierr, idv

      if (myid==0) then
         start = (/ 1, 1, iarch /)
         count = (/ ik, jk, 1 /)

c     get variable idv
         idv = ncvid(idnc,name,ier)
         if(ier.ne.0)then
            print *,'***absent field for idnc,name,idv,ier: ',
     &                               idnc,name,idv,ier
         else
            if(odiag)write(6,*)'idnc,name,idv,ier',idnc,name,idv,ier
            if(odiag)write(6,*)'start=',start
            if(odiag)write(6,*)'count=',count
c     read in all data
            call ncvgt(idnc,idv,start,count,ivar,ier)
            if(odiag)write(6,*)'ivar(1)(ik*jk)=',ivar(1),ivar(ik*jk)

c     obtain scaling factors and offsets from attributes
            call ncagt(idnc,idv,'add_offset',addoff,ier)
            if(odiag)write(6,*)'addoff,ier=',addoff,ier
            call ncagt(idnc,idv,'scale_factor',sf,ier)
            if(odiag)write(6,*)'sf,ier=',sf,ier

c     unpack data
            globvar = ivar*sf+addoff
            if(mod(ktau,nmaxpr).eq.0.or.odiag) then
               vmax = maxval(globvar)
               vmin = minval(globvar)
               write(6,'("done histrd1 ",a6,i4,i3,3e14.6)')
     &           name,ier,iarch,vmin,vmax,globvar(id+(jd-1)*il_g)
            end if
         end if ! ier
      end if ! myid == 0
      ! Have to return correct value of ier on all processes because it's 
      ! used for initialisation in calling routine
      call MPI_Bcast(ier,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if (ier==0) then
         if ( myid == 0 ) then
            call ccmpi_distribute(var,globvar)
         else
            call ccmpi_distribute(var)
         end if
      end if

      return ! histrd1
      end
c***************************************************************************
      subroutine histrd4(idnc,iarch,ier,name,ik,jk,kk,var)
      use cc_mpi
      implicit none
      include 'newmpar.h'
      include 'netcdf.inc'
      include 'parm.h'
      include 'mpif.h'
      integer idnc, iarch, ier, ik, jk, kk
      integer*2 ivar(ik*jk,kk)
      character name*(*)
      integer start(4),count(4)
      real var(ifull,kl)
      real globvar(ifull_g,kl), vmax, vmin, addoff, sf
      integer ierr, ni, nj, idv, i, iq, inij, j, jin, k, iin, inijk

      if (myid == 0 ) then

         start = (/ 1, 1, 1, iarch /)
         count = (/ il_g, jl_g, kk, 1 /)

c     get variable idv
         idv = ncvid(idnc,name,ier)
         if(ier.ne.0)then
            print *,'***absent hist4 field for idnc,name,idv,ier: ',
     &                                     idnc,name,idv,ier
         else

c read in all data
            call ncvgt(idnc,idv,start,count,ivar,ier)

c obtain scaling factors and offsets from attributes
            call ncagt(idnc,idv,'add_offset',addoff,ier)
            call ncagt(idnc,idv,'scale_factor',sf,ier)

c unpack data
            globvar = ivar*sf + addoff
            if(mod(ktau,nmaxpr).eq.0) then
               vmax = maxval(globvar)
               vmin = minval(globvar)
               write(6,'("done histrd4 ",a6,i4,i3,3f12.4)') 
     &           name,ier,iarch,vmin,vmax,globvar(id+(jd-1)*il_g,nlv)
            end if
         end if ! ier
      end if ! myid == 0
      ! Have to return correct value of ier on all processes because it's 
      ! used for initialisation in calling routine
      call MPI_Bcast(ier,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if (ier==0) then
         if ( myid == 0 ) then
            call ccmpi_distribute(var,globvar)
         else
            call ccmpi_distribute(var)
         end if
      end if

      return ! histrd4
      end

      subroutine vertint(t,n)
!     transforms 3d array from dimension kk in vertical to kl   jlm
!     assuming here kk<kl
!     jlm vector special, just with linear new=1 option
      include 'newmpar.h'
      include 'sigs.h'
      include 'parm.h'
      common/work3/dum3(ifull,kl,3),told(ifull,kl),spare(ijk)
      common/sigin/sigin(kl),kk                    ! for vertint, infile
      dimension t(ifull,kl),ka(kl),kb(kl),wta(kl),wtb(kl)
      save num,ka,kb,wta,wtb,klapse
      data num/0/,klapse/0/
      if(num.eq.0)then
        num=1
        do k=1,kl
         if(sig(k).ge.sigin(1))then
           ka(k)=2
           kb(k)=1
           wta(k)=0.
           wtb(k)=1.
           klapse=k   ! i.e. T lapse correction for k.le.klapse
         elseif(sig(k).le.sigin(kk))then   ! at top
           ka(k)=kk
           kb(k)=kk-1
           wta(k)=1.
           wtb(k)=0.
         else
           do kin=2,kk
            if(sig(k).gt.sigin(kin))go to 5
           enddo     ! kin loop
5          ka(k)=kin
           kb(k)=kin-1
           wta(k)=(sigin(kin-1)-sig(k))/(sigin(kin-1)-sigin(kin))
           wtb(k)=(sig(k)-sigin(kin)  )/(sigin(kin-1)-sigin(kin))
         endif  !  (sig(k).ge.sigin(1)) ... ...
        enddo   !  k loop
        print *,'in vertint kk,kl ',kk,kl
        print 91,(sigin(k),k=1,kk)
91      format('sigin',10f7.4)
        print 92,sig
92      format('sig  ',10f7.4)
        print *,'ka ',ka
        print *,'kb ',kb
        print 93,wta
93      format('wta',10f7.4)
        print 94,wtb
94      format('wtb',10f7.4)
      endif     !  (num.eq.0)

      told(:,:)=t(:,:)
      do k=1,kl
       do iq=1,ifull
!        N.B. "a" denotes "above", "b" denotes "below"
         t(iq,k)=wta(k)*told(iq,ka(k))+wtb(k)*told(iq,kb(k))
       enddo   ! iq loop
      enddo    ! k loop

      if(n.eq.1.and.klapse.ne.0)then  ! for T lapse correction
        do k=1,klapse
         do iq=1,ifull
!         assume 6.5 deg/km, with dsig=.1 corresponding to 1 km
           t(iq,k)=t(iq,k)+(sig(k)-sigin(1))*6.5/.1
         enddo   ! iq loop
        enddo    ! k loop
      endif

      if(n.eq.2)then  ! for moisture do a -ve fix
        t(:,:)=max(t(:,:),1.e-6)
      endif
      return
      end

      subroutine datefix(kdate_r,ktime_r,mtimer_r)
      common/leap_yr/leap  ! 1 to allow leap years
      integer mdays(12)
      data mdays/31,28,31,30,31,30,31,31,30,31,30,31/
      data minsday/1440/,minsyr/525600/
      if(kdate_r.ge.00600000.and.kdate_r.le.00991231)then   ! old 1960-1999
        kdate_r=kdate_r+19000000
        print *,'For Y2K kdate_r altered to: ',kdate_r
      endif

      if(mtimer_r.eq.0) then
         print*,'mtimer_r.eq.0: so return in datefix'
         return
      endif
      iyr=kdate_r/10000
      imo=(kdate_r-10000*iyr)/100
      iday=kdate_r-10000*iyr-100*imo
      ihr=ktime_r/100
      imins=ktime_r-100*ihr
      print *,'entering datefix iyr,imo,iday,ihr,imins,mtimer_r: ',
     .                          iyr,imo,iday,ihr,imins,mtimer_r
      do while (mtimer_r.gt.minsyr)
       iyr=iyr+1
       mtimer_r=mtimer_r-minsyr
      enddo
      print *,'a datefix iyr,imo,iday,ihr,imins,mtimer_r: ',
     .                   iyr,imo,iday,ihr,imins,mtimer_r

      mdays(2)=28
      if(mod(iyr,4).eq.0.and.leap.eq.1)mdays(2)=29
      do while (mtimer_r.gt.minsday*mdays(imo))
       mtimer_r=mtimer_r-minsday*mdays(imo)
       imo=imo+1
       if(imo.gt.12)then
         imo=1
         iyr=iyr+1
         if(mod(iyr,4).eq.0.and.leap.eq.1)mdays(2)=29
       endif
      enddo
      print *,'b datefix iyr,imo,iday,ihr,imins,mtimer_r: ',
     .                   iyr,imo,iday,ihr,imins,mtimer_r
      do while (mtimer_r.gt.minsday)
       mtimer_r=mtimer_r-minsday
       iday=iday+1
      enddo
      print *,'c datefix iyr,imo,iday,ihr,imins,mtimer_r: ',
     .                   iyr,imo,iday,ihr,imins,mtimer_r

!     at this point mtimer_r has been reduced to fraction of a day
      mtimerh=mtimer_r/60
      mtimerm=mtimer_r-mtimerh*60  ! minutes left over
      ihr=ihr+mtimerh
      imins=imins+mtimerm
      if(imins.eq.58.or.imins.eq.59)then
!       allow for roundoff for real timer from old runs
        print *,'*** imins increased to 60 from imins = ',imins
        imins=60
      endif
      if(imins.gt.59)then
        imins=imins-60
        ihr=ihr+1
      endif
      print *,'d datefix iyr,imo,iday,ihr,imins,mtimer_r: ',
     .                   iyr,imo,iday,ihr,imins,mtimer_r
      if(ihr.gt.23)then
        ihr=ihr-24
        iday=iday+1
      endif
      print *,'e datefix iyr,imo,iday,ihr,imins,mtimer_r: ',
     .                   iyr,imo,iday,ihr,imins,mtimer_r

      if(iday.gt.mdays(imo))then
         iday=iday-mdays(imo)
         imo=imo+1
         if(imo.gt.12)then
           imo=imo-12
           iyr=iyr+1
         endif
      endif

      kdate_r=iday+100*(imo+100*iyr)
      ktime_r=ihr*100+imins
      mtimer=0
      print *,'leaving datefix iyr,imo,iday,ihr,imins,mtimer_r: ',
     .                         iyr,imo,iday,ihr,imins,mtimer_r
      print *,'leaving datefix kdate_r,ktime_r: ',kdate_r,ktime_r

      return
      end
