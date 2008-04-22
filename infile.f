      subroutine infile(nested,kdate_r,ktime_r,timeg_r,ds_r,
     .                  psl,zs,tss,sicedep,fracice,t,u,v,qg,
!     following not used or returned if called by nestin (i.e.nested=1)   
     .                  tgg,wb,wbice,alb,snowd,
     .                  tggsn,smass,ssdn,ssdnn,snage,isflag,
     .                  isoilh,urban) ! MJT lsmask ! MJT urban 
!     note kk; vertint.f is attached below
!     kdate_r and ktime_r are returned from this routwine.
!     They must equal or exceed kdate_s and ktime_s
!     Note that kdate_s and ktime_s are updated progressively in this 
!     routine to allow for nesting
!     nested=0  for calls from indata; 1  for calls from nestin     
!     From March 2000 exact kdate, ktime are calculated and compared

!     This is called from nestin and onthefly 

      use cc_mpi
      use diag_m
      implicit none
      include 'newmpar.h'
      include 'const_phys.h'
      include 'darcdf.h'
      include 'kuocom.h'    ! ldr
      include 'latlong.h'    
      include 'liqwpar.h'   ! qfg,qlg
      include 'netcdf.inc'
      include 'parm.h'
      include 'stime.h'     ! kdate_s,ktime_s  sought values for data read
      include 'tracers.h'
      include 'mpif.h'

      integer kdate_r, ktime_r, nested
      real timeg_r, ds_r
      real psl(ifull),zs(ifull),
     . tss(ifull),wb(ifull,ms),wbice(ifull,ms),
     . snowd(ifull),alb(ifull),sicedep(ifull),fracice(ifull),
     . t(ifull,kl),u(ifull,kl),v(ifull,kl),qg(ifull,kl),
     . tgg(ifull,ms),tggsn(ifull,3),smass(ifull,3),ssdn(ifull,3),
     . ssdnn(ifull),snage(ifull),urban(ifull,1:12) ! MJT urban
      integer isoilh(ifull) ! MJT lsmask     
      integer isflag(ifull)
      integer ktau_r, ibb, jbb, i
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
      real tds,difknbdi,rhkuo,tdu,ttanl,trnml,tssadd,tstl1,tstl2
      common/rhead/tds,difknbdi,rhkuo,tdu,ttanl,trnml,tstl1,tstl2
      real rlong0x,rlat0x,schmidtx ! infile, infile, nestin, indata
      common/schmidtx/rlong0x,rlat0x,schmidtx ! infile, infile, nestin, indata
      real sigin
      integer ik,jk,kk
      common/sigin/ik,jk,kk,sigin(kl)  ! for vertint, infile
      real tmp(ifull)

      integer, parameter :: nihead=54,nrhead=14
      integer nahead(nihead)
      integer ::  kdatea,kdateb,ktimea,ktimeb
      integer,save ::mtimer_use 
      real ahead(nrhead)
      equivalence (nahead,ini),(ahead,tds)

      integer timest(2),timeco(2),timerco(2)
      character*80 dimnam
      character*3 monn(12),cmon,trnum
      data monn/'jan','feb','mar','apr','may','jun'
     &         ,'jul','aug','sep','oct','nov','dec'/

      integer, save :: ncidold=-1, ncalled=0, iarchi
      integer itype, ilen, ier, ierr, k, ndim, nvars, ngatts,
     &        nulid, nd, isiz, ix, iy, narch, idy, imo, iyr, ihr,
     &        mtimer_in, iq, idv
!     rml 16/02/06 declare igas integer
      integer igas
      real dss, timer
      entry infil(nested,kdate_r,ktime_r,timeg_r,ds_r,
     .                  psl,zs,tss,sicedep,fracice,t,u,v,qg)

      ncalled=ncalled+1
      
c     save model map projs.
      dss=ds_r

      if(myid==0)then
!        N.B. ncid,ncidold and iarchi only get used for myid=0      
!        ncid=idifil
         if(ncid.ne.ncidold)iarchi=1
         ncidold=ncid
         print *,'infile ncid,ncidold,iarchi ',ncid,ncidold,iarchi
         call ncainq(ncid,ncglobal,'int_header',itype,ilen,ier)
         call ncagt(ncid,ncglobal,'int_header',nahead,ier)
         call ncainq(ncid,ncglobal,'real_header',itype,ilen,ier)
         call ncagt(ncid,ncglobal,'real_header',ahead,ier)
      endif  ! (myid==0)
      call MPI_Bcast(nahead,nihead,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(ahead,nrhead,MPI_REAL,0,MPI_COMM_WORLD,ierr)

      ibb=iaa2
      jbb=jaa2
      ntsur2=ntsuri
      rlong0x =ahead(5)
      rlat0x  =ahead(6)
      schmidtx=ahead(7)
      if(myid==0)print *,'infile  rlong0x,rlat0x,schmidtx ',
     &                           rlong0x,rlat0x,schmidtx
      if(schmidtx<=0..or.schmidtx>1.)then
!      read it where Jack put it
       rlong0x =ahead(6)
       rlat0x  =ahead(7)
       schmidtx=ahead(8)
       if(myid==0)print *,'infile  rlong0x,rlat0x,schmidtx ',
     &                            rlong0x,rlat0x,schmidtx
      endif  ! (schmidtx<=0..or.schmidtx>1.)
      if(ktau<=1)then
        if(myid==0)then
          print *,'nahead ',nahead
          print *,'ahead ',ahead
          write(6,'("ntsur2,ibb,jbb",7i6)')
     &               ntsur2,ibb,jbb
        endif  ! (myid==0)
        if(ahead(7)==0.)then
          if(myid==0)print *,'***** schmidt etc not in netcdf header!!!'
!         assume very old file 	   
          rlong0x =0.
          rlat0x  =0.
          schmidtx=1.
        endif  ! (ahead(7)==0.)then
      endif    ! (ktau<=1)
      
      ik=nahead(1)
      jk=nahead(2)
      kk=nahead(3)
c     turn OFF fatal netcdf errors
      if(myid==0)then
        if(ktau<=1)print *,'in infile; ktau,ik,jk,kk=', ktau,ik,jk,kk
         call ncpopt(0)
         call ncagt(ncid,ncglobal,'sigma',sigin,ier)
         if(ier.ne.0)then
            call ncagt(ncid,ncglobal,'sigma_lev',sigin,ier)
            print *,'ier=',ier
            if(ier.ne.0)then
               idv = ncvid(ncid,'lev',ier)
               if(ier.ne.0)idv = ncvid(ncid,'layer',ier)
               call ncvgt(ncid,idv,1,kk,sigin,ier)
               print *,'ier=',ier
            endif ! (ier.ne.0)
         endif    ! (ier.ne.0)
         if(ktau<=1)write(6,'("sigin=",(9f7.4))') (sigin(k),k=1,kk)
c        Get dimensions
         call ncinq(ncid,ndim,nvars,ngatts,nulid,ier)
         if(ktau<=1)print *,'ncid,ndim,nvars,ngatts,nulid,ier'
     &                      ,ncid,ndim,nvars,ngatts,nulid,ier
         do nd=1,ndim
            call ncdinq(ncid,nd,dimnam,isiz,ier)
            if(ier.ne.0)write(6,*) 'ncdinq dim1 ier=',ier
            if(ktau<=1)write(6,*) 'dim ',nd,' ',isiz,' ',dimnam
            if(nd==1) ix=isiz
            if(nd==2) iy=isiz
            if(dimnam=='time')narch=isiz
         enddo
         if(ktau<=1)write(6,'("ix,iy,narch=",3i6)') ix,iy,narch
      endif   ! (myid==0)
      call MPI_Bcast(sigin,kk,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(ix,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(iy,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(narch,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      timest(1) = 1
      timeco(1) = 3
      timeco(2) = 1
      timerco(1) = 2
      timerco(2) = 1

c***********************************************************************
      if(myid==0)then
        if(ktau<=1)then
!         Following effectively checks whether nc file uses  
!         start time followed by sequence of mtimer (e.g. from prior CCAM output)
!         or has indiviual date/times (e.g. from NCEP or Mk3)      
!         It does not handle a file with single date/time and mtimer > 0
          idv = ncvid(ncid,'kdate',ier)
          call ncvgt1(ncid,idv,1,kdatea,ier) 
          call ncvgt1(ncid,idv,2,kdateb,ier) 
          idv = ncvid(ncid,'ktime',ier)
          call ncvgt1(ncid,idv,1,ktimea,ier) 
          call ncvgt1(ncid,idv,2,ktimeb,ier) 
          mtimer_use=1
          if(kdatea.ne.kdateb.or.ktimea.ne.ktimeb)mtimer_use=0
          if(ier.ne.0)mtimer_use=1 ! for single time, e.g. restart file
          print *,'ktau,ier,kdatea,kdateb,ktimea,ktimeb,mtimer_use ',
     &             ktau,ier,kdatea,kdateb,ktimea,ktimeb,mtimer_use
        endif
        print *,'searching for the first date >= kdate_s,ktime_s ',
     &           kdate_s,ktime_s
        iarchi=iarchi-1
19      iarchi=iarchi+1
        timest(2) = iarchi
c       turn OFF fatal netcdf errors
        call ncpopt(0)
        idv = ncvid(ncid,'timestamp',ier)
c       print *,'timestamp idv,ier,iarchi=',idv,ier,iarchi
c       turn ON fatal netcdf errors
        call ncpopt(NCVERBOS+NCFATAL)
        if(ier==0)then
           timeco(1)=1
           timeco(2)=iarchi
           call ncvgt1(ncid,idv,timeco,kdate_r,ier)
           timeco(1)=2
           timeco(2)=iarchi
           call ncvgt1(ncid,idv,timeco,ktime_r,ier)
           timeco(1)=3
           timeco(2)=iarchi
           call ncvgt1(ncid,idv,timeco,ktau_r,ier)
           idv = ncvid(ncid,'realstamp',ier)
           timeco(1)=1
           timeco(2)=iarchi
           call ncvgt1(ncid,idv,timeco,timer,ier)
!          N.B. usual ncdump does not show full precision of floating variables
           timeco(1)=2
           timeco(2)=iarchi
           call ncvgt1(ncid,idv,timeco,timeg_r,ier)
        else    ! ier.ne.0
!          new style
           idv = ncvid(ncid,'kdate',ier)
           call ncvgt1(ncid,idv,iarchi,kdate_r,ier)
!          turn off fatal ncdf errors
           call ncpopt(0)
           idv = ncvid(ncid,'timer',ier)
           if(ier==0)call ncvgt1(ncid,idv,iarchi,timer,ier)
           idv = ncvid(ncid,'mtimer',ier)
           if(ier==0)call ncvgt1(ncid,idv,iarchi,mtimer,ier)
!          turn on fatal ncdf errors
           call ncpopt(NCVERBOS+NCFATAL)
           idv = ncvid(ncid,'ktime',ier)
           call ncvgt1(ncid,idv,iarchi,ktime_r,ier)
           idv = ncvid(ncid,'ktau',ier)
           call ncvgt1(ncid,idv,iarchi,ktau_r,ier)
           idv = ncvid(ncid,'timeg',ier)
           call ncvgt1(ncid,idv,iarchi,timeg_r,ier)
        endif  ! (ier==0) .. else ..

        if(kdate_r==0)then
           print *,
     &      'kdate_r = 0; create kdate from time_origin, nested=',nested
           idv = ncvid(ncid,'time',ier)
!          print *,'time idv=',idv,ier
           call ncpopt(0)
           call ncagtc(ncid,idv,'time_origin',dimnam,80,ier)
           print *,'time_origin=',dimnam,ier
           call ncpopt(NCVERBOS+NCFATAL)
           read(dimnam,'(i2)') idy
           read(dimnam,'(3x,a3)') cmon
           do imo=1,12
              if(cmon==monn(imo))go to 10
           enddo
 10        continue
           read(dimnam,'(7x,i4)') iyr
           read(dimnam,'(12x,i2)') ihr
           kdate_r=idy+100*(imo+100*iyr)
           ktime_r=ihr*100
        endif  ! (kdate_r==0)
        print *,'ktau,iarchi,kdate,ktime,timer_in,mtimer: ',
     &           ktau,iarchi,kdate_r,ktime_r,timer,mtimer
        if(mtimer.ne.0)then     ! preferred
!         assume mtimer is read in correctly and set timer from that
          timer=mtimer/60.
        else                    ! mtimer = 0
          mtimer=nint(60.*timer)
        endif  ! (mtimer.ne.0)
        mtimer_in=mtimer       
                 
        if(mtimer_use==1.and.mtimer>0)then
          print *,'using timer, mtimer: ',timer,mtimer
          call datefix(kdate_r,ktime_r,mtimer) ! for Y2K, or mtimer>0
        endif
        if(2400*kdate_r+ktime_r<
     &     2400*kdate_s+1200*nsemble+ktime_s)go to 19  ! 12-h nsemble
!---------------------------------------------------------------------
        
        print *,'found suitable date/time in infile'
        if(nsemble.ne.0)then
          kdate_r=kdate_s
          ktime_r=ktime_s
          print *,'for nsemble = ',nsemble,
     &      ' resetting kdate_r,ktime_r to ',kdate_r,ktime_r
        endif
        print *,'in ',ik,jk,kk,m2,nsd2,nbd2
     & ,nps2,mex2,mup2,mtimer,nmi2,ndt2,npsav2,nhor2,nkuo2,khdif2
     & ,kwt2,iaa2,jaa2,nvad2,lbd2,nrun2,nrunx2
     & ,khor2,ksc2,kountr2,ndiur2,nhort2,nhorps2,nsoil2,ivirt2
     & ,ntsur2,nrad2,kuocb2,nvmix2,ntsea2,ms2,nextras2,ilt2,ntrac2

        print *,'infile ds_r,mtimer: ',ds_r,mtimer
        if(mtimer_in>0.and.nrungcm==2)then
           write(0,*)
     &          '*** re-setting NRUNGCM nrungcm from 2 to 0 because ',
     &             'not starting from gcm file'
           stop 'wrong nrungcm used'
        endif  ! (mtimer_in>0.and.nrungcm==2)
      endif    ! myid == 0
      
      ! Now share the time variables (is this actually necessary?)
      call MPI_Bcast(kdate_r,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(ktime_r,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(ktau_r,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(mtimer,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(timer,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(timeg_r,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)

c     begin reading data
c     log scaled sfc.press
      call histrd1(ncid,iarchi,ier,'psf',ik,jk,psl)
      call histrd1(ncid,iarchi,ier,'zht',ik,jk,zs)
      call histrd1(ncid,iarchi,ier,'tsu',ik,jk,tss)
!***  tss(:)=abs(tss(:)) ! not here because -ves needed for onthefly
      if(nspecial<-50)then
        tssadd=.1*(nspecial+100) ! e.g. -110 makes 1 degree cooler
        do iq=1,ifull
         if(tss(iq)<0.)then
           tss(iq)=tss(iq)-tssadd
         elseif(zs(iq)<0.)then
           tss(iq)=tss(iq)+tssadd
         endif
        enddo
      endif
      if(nspecial==21)then
        do iq=1,ifull
         if(rlongg(iq)*180./pi>100.and.rlongg(iq)*180./pi<260..and.
     &      rlatt(iq)*180./pi>-20.and.rlatt(iq)*180./pi<20.)
     &   tss(iq)=.999*abs(tss(iq))
        enddo
      endif
      if(nspecial==22)then
        do iq=1,ifull
         if(rlongg(iq)*180./pi>100.and.rlongg(iq)*180./pi<260..and.
     &      rlatt(iq)*180./pi>-20.and.rlatt(iq)*180./pi<20.)
     &   tss(iq)=.995*abs(tss(iq))
        enddo
      endif
      if(nspecial==23)then
        do iq=1,ifull
         if(rlongg(iq)*180./pi>100.and.rlongg(iq)*180./pi<260..and.
     &      rlatt(iq)*180./pi>-20.and.rlatt(iq)*180./pi<20.)
     &   tss(iq)=.990*abs(tss(iq))
        enddo
      endif

c     turn on fatal netcdf errors
      if(myid == 0)call ncpopt(NCVERBOS+NCFATAL)
c     temperature
      call histrd4(ncid,iarchi,ier,'temp',ik,jk,kk,t)
c     u wind component
      call histrd4(ncid,iarchi,ier,'u',ik,jk,kk,u)
c     v wind component
      call histrd4(ncid,iarchi,ier,'v',ik,jk,kk,v)
c     mixing ratio
c     turn OFF fatal netcdf errors; from here on
      if(myid == 0)call ncpopt(0)
      call histrd4(ncid,iarchi,ier,'q',ik,jk,kk,qg)       !  kg/kg
      if(ier.ne.0)then
        call histrd4(ncid,iarchi,ier,'mixr',ik,jk,kk,qg)  !   g/kg
      endif  ! (ier.ne.0)
      sicedep(:)=0.
      fracice(:)=0.
      call histrd1(ncid,iarchi,ier,'siced',ik,jk,sicedep)  
      call histrd1(ncid,iarchi,ierr,'fracice',ik,jk,fracice)  
      if(ier==0)then  ! i.e. sicedep read in 
        do iq=1,ifull
         if(sicedep(iq)<.05)then ! for sflux
           sicedep(iq)=0.
           fracice(iq)=0.
         endif
        enddo
        if(ierr.ne.0)then ! i.e. sicedep read in; fracice not read in
          where(sicedep>0.)
            fracice=1.
          endwhere
        endif  ! (ierr.ne.0)  fracice
      endif    ! (ier==0)
      if(ier.ne.0)then     ! sicedep not read in
        if(ierr.ne.0)then  ! neither sicedep nor fracice read in
          fracice(:)=0.
	   if(myid==0) print *,'pre-setting siced in infile from tss'
          where(abs(tss) <= 271.2)
            fracice=1.
          endwhere
c       else  ! i.e. only fracice read in;  done in indata, nestin
c         do iq=1,ifull 
c          if(fracice(iq)>.01)then
c            sicedep(iq)=2.
c          else
c            sicedep(iq)=0.
c            fracice(iq)=0.
c          endif
c         enddo
        endif  ! (ierr.ne.0)
      endif    ! (ier.ne.0) .. else ..    for sicedep
      if(ncalled<4.and.mydiag)then
        print *,'sig in: ',(sigin(i),i=1,kk)
        write (6,"('100*psl# in',9f8.2)") 100.*diagvals(psl)
        write (6,"('zs# in     ',9f8.1)") diagvals(zs)
        write (6,"('tss# in    ',9f8.1)") diagvals(tss)
!     .           ((tss(ii+(jj-1)*il),ii=id-1,id+1),jj=jd-1,jd+1)
        print *,'ier,ierr for siced, fracice ',ier,ierr
        write (6,"('siced# in  ',9f8.2)") diagvals(sicedep)
        write (6,"('fracice# in',9f8.2)") diagvals(fracice)
!       Printing the ifull value gives confusing results in the 
!       parallel version
        print *,'t in ',(t(idjd,k),k=1,kk)!,t(ifull,kk)
        print *,'u in ',(u(idjd,k),k=1,kk)!,u(ifull,kk)
        print *,'v in ',(v(idjd,k),k=1,kk)!,v(ifull,kk)
        print *,'qg in ',(qg(idjd,k),k=1,kk)!,qg(ifull,kk)
      endif  ! (ncalled<4.and.mydiag)

      if(myid == 0)print *,'in infile nested = ',nested
!########################################################################
      if(nested==0)then   !  following only at start of run 
!       read fields which may be there on restart, but maybe not initially
        if(ldr.ne.0)then   ! following in kg/kg
          call histrd4(ncid,iarchi,ier,'qfg',ik,jk,kk,qfg(1:ifull,:))
          call histrd4(ncid,iarchi,ier,'qlg',ik,jk,kk,qlg(1:ifull,:))
          if(ier.ne.0.or.ik<il_g)then
!           set default qfg, qlg to zero if not available on input 
!           for onthefly usage of qfg, qlg need extra stuff there for ik<il_g	
            print *,'qfg,qlg being set to default in infile'
            qfg(:,:)=0. 
            qlg(:,:)=0. 
          endif ! (ier.ne.0)
          if(nspecial==24)qfg(:,kl-2)=4.e-3         ! just for testing icefall
          if(nspecial==25)qfg((il-1)*il/2,kl)=4.e-3 ! tests at (il/2,il/2)
          if(nspecial==26)qfg(idjd,kl)=4.e-3 ! tests at (idjd)
        endif   ! (ldr.ne.0)
        call histrd1(ncid,iarchi,ier,'alb',ik,jk,alb)
        call histrd1(ncid,iarchi,ierr,'tgg2',ik,jk,tgg(1,2))
        if(ierr==0)then  ! at least tgg6, wb2, wb6 will be available
          call histrd1(ncid,iarchi,ier,'tgg6',ik,jk,tgg(1,6))
          call histrd1(ncid,iarchi,ier,'wb2',ik,jk,wb(1,2))
          call histrd1(ncid,iarchi,ier,'wb6',ik,jk,wb(1,6))
          call histrd1(ncid,iarchi,ier,'tgg1',ik,jk,tgg(1,1)) ! trial read
          if(ier==0)then
            call histrd1(ncid,iarchi,ier,'tgg3',ik,jk,tgg(1,3))
            call histrd1(ncid,iarchi,ier,'tgg4',ik,jk,tgg(1,4))
            call histrd1(ncid,iarchi,ier,'tgg5',ik,jk,tgg(1,5))
            call histrd1(ncid,iarchi,ier,'wb1',ik,jk,wb(1,1))
            call histrd1(ncid,iarchi,ier,'wb3',ik,jk,wb(1,3))
            call histrd1(ncid,iarchi,ier,'wb4',ik,jk,wb(1,4))
            call histrd1(ncid,iarchi,ier,'wb5',ik,jk,wb(1,5))
          else  ! set other levels having read tgg2, tgg6, wb2, wb6
            tgg(:,1)=abs(tss(:)) ! initial temp. at 2nd layer(6.5.97 KF)
            tgg(:,3)=tgg(:,2) ! initial temp.. from GCM runs with 3 layers
            wb(:,1) =wb(:,2)  ! layer initialisation of moisture
            do k=3,ms-1       ! don't want to change  tgg(,2) and tgg(,ms)
             do iq=1,ifull
              wb(iq,k) =wb(iq,ms)   ! layer initialisation of moisture
             enddo  ! iq
            enddo   ! k
            do k=4,ms-1
             do iq=1,ifull
              tgg(iq,k)=tgg(iq,ms) ! init temp from GCM runs with 3 layers
             enddo  ! iq
            enddo   ! k
          endif     ! (ier==0) .. else ..
        else        ! i.e. ierr.ne.0
          call histrd1(ncid,iarchi,ier,'tb2',ik,jk,tgg(1,ms))
          call histrd1(ncid,iarchi,ier,'tb3',ik,jk,tgg(1,2))
          call histrd1(ncid,iarchi,ier,'wfg',ik,jk,wb(1,2))
          call histrd1(ncid,iarchi,ier,'wfb',ik,jk,wb(1,ms))
          tgg(:,1)=abs(tss(:)) ! initial temp at 2nd layer(6.5.97 KF)
          tgg(:,3)=tgg(:,2) ! initial temp. from GCM runs with 3 layers
          wb(:,1) =wb(:,2)  ! layer initialisation of moisture
          do k=3,ms-1       ! don't want to change  tgg(,2) and tgg(,ms)
           do iq=1,ifull
            wb(iq,k) =wb(iq,ms)   ! layer initialisation of moisture
           enddo  ! iq
          enddo   ! k
          do k=4,ms-1
           do iq=1,ifull
            tgg(iq,k)=tgg(iq,ms) ! init temp. from GCM runs with 3 layers
           enddo  ! iq
          enddo   ! k
          if(nproc==1)then
c 	    only being tested for nested=0; no need to test for mesonest
            if(tgg(1,1)<200..or.tgg(1,2)<200..or.tgg(1,6)<200.)then
              write(0,*) 'impossibly cold initial tgg1, tgg2 or tgg6' 
              write(0,*)  tgg(1,1),tgg(1,1),tgg(1,1)
              stop
            endif     
          endif     
        endif    ! (ierr==0) .. else ..
      !----------------------------------------------------------------
      ! MJT CHANGE - ECOSYSTEMS
        call histrd1(ncid,iarchi,ierr,'wetfrac1',ik,jk,wb(:,1))
        if (ierr==0) then
          call histrd1(ncid,iarchi,ierr,'wetfrac2',ik,jk,wb(:,2))
          call histrd1(ncid,iarchi,ierr,'wetfrac3',ik,jk,wb(:,3))
          call histrd1(ncid,iarchi,ierr,'wetfrac4',ik,jk,wb(:,4))
          call histrd1(ncid,iarchi,ierr,'wetfrac5',ik,jk,wb(:,5))
          call histrd1(ncid,iarchi,ierr,'wetfrac6',ik,jk,wb(:,6))
          wb(:,:)=-abs(wb(:,:)) ! flag to indicate wetfrac, not volumetric soil moisture (unpacked in indata)
        end if
      !----------------------------------------------------------------

        !------------------------------------------------------------
        ! MJT urban
        urban(:,:)=999. ! this must equal spval in onthefly.f
        call histrd1(ncid,iarchi,ierr,'rooftgg1',ik,jk,urban(:,1))
        if (ierr==0) then
          call histrd1(ncid,iarchi,ierr,'rooftgg2',ik,jk,urban(:,2))
          call histrd1(ncid,iarchi,ierr,'rooftgg3',ik,jk,urban(:,3))
          call histrd1(ncid,iarchi,ierr,'waletgg1',ik,jk,urban(:,4))
          call histrd1(ncid,iarchi,ierr,'waletgg2',ik,jk,urban(:,5))
          call histrd1(ncid,iarchi,ierr,'waletgg3',ik,jk,urban(:,6))
          call histrd1(ncid,iarchi,ierr,'walwtgg1',ik,jk,urban(:,7))
          call histrd1(ncid,iarchi,ierr,'walwtgg2',ik,jk,urban(:,8))
          call histrd1(ncid,iarchi,ierr,'walwtgg3',ik,jk,urban(:,9))
          call histrd1(ncid,iarchi,ierr,'roadtgg1',ik,jk,urban(:,10))
          call histrd1(ncid,iarchi,ierr,'roadtgg2',ik,jk,urban(:,11))
          call histrd1(ncid,iarchi,ierr,'roadtgg3',ik,jk,urban(:,12))
        end if
        !------------------------------------------------------------        
        
        !------------------------------------------------------------
        ! MJT lsmask
        tmp(:)=-1.
        call histrd1(ncid,iarchi,ierr,'soilt',ik,jk,tmp(:))
        isoilh(:)=nint(tmp(:))
        !------------------------------------------------------------

        if(myid == 0)print *,'about to read snowd'
        call histrd1(ncid,iarchi,ier,'snd',ik,jk,snowd)
        if(ier.ne.0)then  ! preset snowd here if not avail.
!         when no snowd available initially, e.g. COMPARE III (jlm formula)
          snowd(:)=0.      ! added Feb '05
          do iq=1,ifull    ! "do" missing till 5/6/01
!          for this one,snow lin. increases from 5 to 55, for T 270 down to 260
           if(abs(tss(iq))<270.)
     &                   snowd(iq)=min(55.,5.*(271.-abs(tss(iq))))
          enddo
          if(ncalled<4.and.mydiag)then
           print *,'setting snowd in infile, ier= '  ,ier
           print *,'snowd# preset to: ', diagvals(snowd)
!     .            ((snowd(ii+(jj-1)*il),ii=id-1,id+1),jj=jd-1,jd+1)
          endif ! (ncalled<4.and.mydiag)
        endif   ! (ier.ne.0)
        do iq=1,ifull
         if(snowd(iq)<.001)snowd(iq)=0. ! for netcdf rounding 7/6/06
        enddo
        call histrd1(ncid,iarchi,ier,'smass1',ik,jk,smass(1,1))
        if(ier.ne.0)then  ! for smass1
          if(myid==0)
     &          print *,'setting smass,wbice etc in infile, ier= '  ,ier
          smass(:,:)=0.
          tggsn(:,:)=280.     ! just a default
          isflag(:) = 0
          snage(:)  = 0.
          do iq=1,ifull
           if(snowd(iq)>100.)then
             ssdn(iq,1)=240.
           else
             ssdn(iq,1) = 140.
           endif  ! (snowd(iq)>100.)
           do k=2,3
            ssdn(iq,k)=ssdn(iq,1)
           enddo
           ssdnn(iq)  = ssdn(iq,1)
           if(snowd(iq)>0.)tgg(iq,1)=min(tgg(iq,1),270.1)
          enddo   ! iq loop
!         note permanent wbice may be re-set in indata
          wbice(:,:)=0.
        else     ! assume all these variables available in restart file
         call histrd1(ncid,iarchi,ier,'smass2',ik,jk,smass(1,2))
         call histrd1(ncid,iarchi,ier,'smass3',ik,jk,smass(1,3))
         call histrd1(ncid,iarchi,ier,'ssdn1',ik,jk,ssdn(1,1))
         call histrd1(ncid,iarchi,ier,'ssdn2',ik,jk,ssdn(1,2))
         call histrd1(ncid,iarchi,ier,'ssdn3',ik,jk,ssdn(1,3))
         call histrd1(ncid,iarchi,ier,'tggsn1',ik,jk,tggsn(1,1))
         call histrd1(ncid,iarchi,ier,'tggsn2',ik,jk,tggsn(1,2))
         call histrd1(ncid,iarchi,ier,'tggsn3',ik,jk,tggsn(1,3))
         call histrd1(ncid,iarchi,ier,'wbice1',ik,jk,wbice(1,1))
         call histrd1(ncid,iarchi,ier,'wbice2',ik,jk,wbice(1,2))
         call histrd1(ncid,iarchi,ier,'wbice3',ik,jk,wbice(1,3))
         call histrd1(ncid,iarchi,ier,'wbice4',ik,jk,wbice(1,4))
         call histrd1(ncid,iarchi,ier,'wbice5',ik,jk,wbice(1,5))
         call histrd1(ncid,iarchi,ier,'wbice6',ik,jk,wbice(1,6))
         call histrd1(ncid,iarchi,ier,'snage',ik,jk,snage)
         call histrd1(ncid,iarchi,ier,'sflag',ik,jk,tmp)
         isflag(:)=nint(tmp(:))
        endif  ! (ier.ne.0) ... else ...   for smass1

!       rml 16/02/06 read tracer from restart for up to 999 tracers
        if (ngas>0) then
          do igas=1,ngas
            write(trnum,'(i3.3)') igas
            call histrd4(ncid,iarchi,ier,'tr'//trnum,ik,jk,kk,
     &                   tr(1:ifull,:,igas))
          enddo
        endif

        if(mydiag)then
          print *,'at end of infile kdate,ktime,ktau,tss: ',
     &                             kdate_r,ktime_r,ktau,tss(idjd)
          print *,'tgg ',(tgg(idjd,k),k=1,ms)
          print *,'wb ',(wb(idjd,k),k=1,ms)
          print *,'wbice ',(wbice(idjd,k),k=1,ms)
        endif ! (mydiag)
        if(ncalled<4.and.mydiag)then
          write (6,"('tgg(1)# in ',9f8.1)") diagvals(tgg(:,1))
          write (6,"('tgg(2)# in ',9f8.1)") diagvals(tgg(:,2))
          write (6,"('tgg(3)# in ',9f8.1)") diagvals(tgg(:,3))
          write (6,"('tgg(ms)# in',9f8.1)") diagvals(tgg(:,ms))
          write (6,"('wb(1)# in  ',9f8.3)") diagvals(wb(:,1))
          write (6,"('wb(ms)# in ',9f8.3)") diagvals(wb(:,ms))
          write (6,"('alb# in    ',9f8.3)") diagvals(alb)
          write (6,"('snowd# in  ',9f8.2)") diagvals(snowd)
        endif ! (ncalled<4.and.mydiag)
      endif   ! (nested==0)   !  only done at start of run 
! ########################################################################

      iarchi=iarchi+1
      if(mydiag)then
         write(6,'("end infile kdate,ktime,ktau_r,mtimer,timer,ds_r="
     &      ,4i10,2f10.1)')kdate_r,ktime_r,ktau_r,mtimer,timer,ds_r
      endif ! (mydiag)

!     reset kdate_s and ktime_s ready for next (nesting) call to infile
      kdate_s=kdate_r
      ktime_s=ktime_r+1 

      qg(1:ifull,1:kk) = max(qg(1:ifull,1:kk),1.e-6)
!     qg(1:ik*jk,1:kk) = max(qg(1:ik*jk,1:kk),1.e-6)

      if(mydiag)then
         print *,'end infile; next read will be kdate_s,ktime_s >= ',
     &                                          kdate_s,ktime_s
      endif  ! (mydiag)

      end subroutine infile

c*************************************************************************
      subroutine histrd1(ncid,iarchi,ier,name,ik,jk,var)
      use cc_mpi
      implicit none
      include 'newmpar.h'
      include 'parm.h'
      include 'netcdf.inc'
      include 'mpif.h'
      integer ncid, iarchi, ier, ik, jk, nctype, ierb ! MJT CHANGE - bug fix
      integer*2 ivar(ik*jk)
      real rvar(ik*jk) ! MJT CHANGE
      logical odiag
      parameter(odiag=.false. )
      character name*(*)
      integer start(3),count(3)
      real var(ifull)
      real globvar(ifull_g), vmax, vmin, addoff, sf
      integer ierr, idv

      if(myid==0)then
         start = (/ 1, 1, iarchi /)
         count = (/ ik, jk, 1 /)
         if(ik<il_g)globvar(:)=2. ! default for onthefly

c     get variable idv
         idv = ncvid(ncid,name,ier)
         if(ier.ne.0)then
            print *,'***absent field for ncid,name,idv,ier: ',
     &                               ncid,name,idv,ier
         else
            if(odiag)write(6,*)'ncid,name,idv,ier',ncid,name,idv,ier
            if(odiag)write(6,*)'start=',start
            if(odiag)write(6,*)'count=',count
c     read in all data
            !------------------------------------------------------------
            ! MJT CHANGE
            ierr=nf_inq_vartype(ncid,idv,nctype)
	    addoff=0.
	    sf=1.
            select case(nctype)
              case(nf_float)
                call ncvgt(ncid,idv,start,count,rvar,ier)
                if(odiag)write(6,*)'rvar(1)(ik*jk)=',rvar(1),rvar(ik*jk)
              case(nf_short)
                call ncvgt(ncid,idv,start,count,ivar,ier)
                if(odiag)write(6,*)'ivar(1)(ik*jk)=',ivar(1),ivar(ik*jk)
                rvar(:)=real(ivar(:))
              case DEFAULT
                if (myid == 0) print *,"ERROR: Unknown NetCDF format"
                stop
            end select
            !------------------------------------------------------------

c     obtain scaling factors and offsets from attributes
            call ncagt(ncid,idv,'add_offset',addoff,ierb)
            if (ierb.ne.0) addoff=0.        ! MJT CHANGE - bug fix
            if(odiag)write(6,*)'addoff,ier=',addoff,ier
            call ncagt(ncid,idv,'scale_factor',sf,ierb)
            if (ierb.ne.0) sf=1.            ! MJT CHANGE - bug fix
            if(odiag)write(6,*)'sf,ier=',sf,ier

c     unpack data
            globvar(1:ik*jk) = rvar(1:ik*jk)*sf+addoff ! MJT CHANGE
            if(mod(ktau,nmaxpr)==0.or.odiag)then
               vmax = maxval(rvar*sf+addoff) ! MJT CHANGE
               vmin = minval(rvar*sf+addoff) ! MJT CHANGE
               write(6,'("done histrd1 ",a8,i4,i3,3e14.6)')
     &           name,ier,iarchi,vmin,vmax,globvar(id+(jd-1)*ik)
            endif
         endif ! ier
      endif ! myid == 0
      ! Have to return correct value of ier on all processes because it's 
      ! used for initialisation in calling routine
      call MPI_Bcast(ier,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if(ier==0)then
         if(myid == 0)then
            call ccmpi_distribute(var,globvar)
         else
            call ccmpi_distribute(var)
         endif
      endif
      return    ! histrd1
      end   
c***************************************************************************
      subroutine histrd4(ncid,iarchi,ier,name,ik,jk,kk,var)
      use cc_mpi
      implicit none
      include 'newmpar.h'
      include 'netcdf.inc'
      include 'parm.h'
      include 'mpif.h'
      integer ncid, iarchi, ier, ik, jk, kk, nctype, ierb ! MJT CHANGE - bug fix
      integer*2 ivar(ik*jk,kk)
      real rvar(ik*jk,kk) ! MJT CHANGE
      character name*(*)
      integer start(4),count(4)
      real var(ifull,kl)
      real globvar(ifull_g,kl), vmax, vmin, addoff, sf
      integer ierr, idv, k

      if(myid == 0)then
         start = (/ 1, 1, 1, iarchi /)
!        count = (/ il_g, jl_g, kk, 1 /)
         count = (/   ik,   jk, kk, 1 /)
         if(ik<il_g)globvar(:,:)=1.234 ! default for onthefly
c        get variable idv
         idv = ncvid(ncid,name,ier)
         if(ier.ne.0)then
            print *,'***absent hist4 field for ncid,name,idv,ier: ',
     &                                         ncid,name,idv,ier
         else
c          read in all data
           !------------------------------------------------------------
           ! MJT CHANGE
           ier=nf_inq_vartype(ncid,idv,nctype)
           select case(nctype)
             case(nf_float)
               call ncvgt(ncid,idv,start,count,rvar,ier)
             case(nf_short)
               call ncvgt(ncid,idv,start,count,ivar,ier)
               rvar(:,:)=real(ivar(:,:))
             case DEFAULT
               if (myid == 0) print *,"ERROR: Unknown NetCDF format"
               stop
           end select
           !------------------------------------------------------------
c          obtain scaling factors and offsets from attributes
           call ncagt(ncid,idv,'add_offset',addoff,ierb)
           if (ierb.ne.0) addoff=0.        ! MJT CHANGE - bug fix
           call ncagt(ncid,idv,'scale_factor',sf,ierb)
           if (ierb.ne.0) sf=1.            ! MJT CHANGE - bug fix
c          unpack data
!          globvar = ivar*sf + addoff
           do k=1,kk  ! following allows for ik < il_g
            globvar(1:ik*jk,k) = rvar(1:ik*jk,k)*sf + addoff ! MJT CHANGE
           enddo
           if(mod(ktau,nmaxpr)==0)then
!            vmax = maxval(globvar)
!            vmin = minval(globvar)
             vmax = maxval(rvar*sf + addoff)  ! allows for ik < il_g ! MJT CHANGE
             vmin = minval(rvar*sf + addoff)                         ! MJT CHANGE
             write(6,'("done histrd4 ",a6,i4,i3,3f12.4)') 
     &           name,ier,iarchi,vmin,vmax,globvar(id+(jd-1)*ik,nlv)
           endif
         endif ! ier
      endif ! myid == 0
      ! Have to return correct value of ier on all processes because it's 
      ! used for initialisation in calling routine
      call MPI_Bcast(ier,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if(ier==0)then
         if(myid == 0)then
            call ccmpi_distribute(var,globvar)
         else
            call ccmpi_distribute(var)
         endif
      endif
      return ! histrd4
      end   

      subroutine vertint(t,n)
!     transforms 3d array from dimension kk in vertical to kl   jlm
!     assuming here kk<kl
!     jlm vector special, just with linear new=1 option
      include 'newmpar.h'
      include 'sigs.h'
      include 'parm.h'
      common/sigin/ik,jk,kk,sigin(kl)  ! for vertint, infile
      dimension t(ifull,kl),ka(kl),kb(kl),wta(kl),wtb(kl)  ! for mpi
      real told(ifull,kl)
      save num,ka,kb,wta,wtb,klapse
      data num/0/,klapse/0/
      if(num==0)then
        num=1
        do k=1,kl
         if(sig(k)>=sigin(1))then
           ka(k)=2
           kb(k)=1
           wta(k)=0.
           wtb(k)=1.
           klapse=k   ! i.e. T lapse correction for k<=klapse
         elseif(sig(k)<=sigin(kk))then   ! at top
           ka(k)=kk
           kb(k)=kk-1
           wta(k)=1.
           wtb(k)=0.
         else
           do kin=2,kk
            if(sig(k)>sigin(kin))go to 5
           enddo     ! kin loop
5          ka(k)=kin
           kb(k)=kin-1
           wta(k)=(sigin(kin-1)-sig(k))/(sigin(kin-1)-sigin(kin))
           wtb(k)=(sig(k)-sigin(kin)  )/(sigin(kin-1)-sigin(kin))
         endif  !  (sig(k)>=sigin(1)) ... ...
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
      endif     !  (num==0)

      told(:,:)=t(:,:)
      do k=1,kl
       do iq=1,ifull
!        N.B. "a" denotes "above", "b" denotes "below"
         t(iq,k)=wta(k)*told(iq,ka(k))+wtb(k)*told(iq,kb(k))
       enddo   ! iq loop
      enddo    ! k loop
c     print *,'in vertint told',told(idjd,:)
c     print *,'in vertint t',t(idjd,:)

      if(n==1.and.klapse.ne.0)then  ! for T lapse correction
        do k=1,klapse
         do iq=1,ifull
!         assume 6.5 deg/km, with dsig=.1 corresponding to 1 km
           t(iq,k)=t(iq,k)+(sig(k)-sigin(1))*6.5/.1
         enddo   ! iq loop
        enddo    ! k loop
      endif

      if(n==2)then  ! for qg do a -ve fix
        t(:,:)=max(t(:,:),1.e-6)
      endif
      if(n==5)then  ! for qfg, qlg do a -ve fix
        t(:,:)=max(t(:,:),0.)
      endif
      return
      end

      subroutine datefix(kdate_r,ktime_r,mtimer_r)
      include 'newmpar.h'
      include 'parm.h'
      common/leap_yr/leap  ! 1 to allow leap years
      integer mdays(12)
      data mdays/31,28,31,30,31,30,31,31,30,31,30,31/
      data minsday/1440/,minsyr/525600/
      if(kdate_r>=00600000.and.kdate_r<=00991231)then   ! old 1960-1999
        kdate_r=kdate_r+19000000
        print *,'For Y2K kdate_r altered to: ',kdate_r
      endif
      iyr=kdate_r/10000
      imo=(kdate_r-10000*iyr)/100
      iday=kdate_r-10000*iyr-100*imo
      ihr=ktime_r/100
      imins=ktime_r-100*ihr
      print *,'entering datefix iyr,imo,iday,ihr,imins,mtimer_r: ',
     .                          iyr,imo,iday,ihr,imins,mtimer_r
      do while (mtimer_r>minsyr)
       iyr=iyr+1
       mtimer_r=mtimer_r-minsyr
      enddo
      if(diag)print *,'a datefix iyr,imo,iday,ihr,imins,mtimer_r: ',
     .                   iyr,imo,iday,ihr,imins,mtimer_r

      mdays(2)=28
      if(mod(iyr,4)==0.and.leap==1)mdays(2)=29
      do while (mtimer_r>minsday*mdays(imo))
       mtimer_r=mtimer_r-minsday*mdays(imo)
       imo=imo+1
       if(imo>12)then
         imo=1
         iyr=iyr+1
         if(mod(iyr,4)==0.and.leap==1)mdays(2)=29
       endif
      enddo
      if(diag)print *,'b datefix iyr,imo,iday,ihr,imins,mtimer_r: ',
     .                   iyr,imo,iday,ihr,imins,mtimer_r
      do while (mtimer_r>minsday)
       mtimer_r=mtimer_r-minsday
       iday=iday+1
      enddo
      if(diag)print *,'c datefix iyr,imo,iday,ihr,imins,mtimer_r: ',
     .                   iyr,imo,iday,ihr,imins,mtimer_r

!     at this point mtimer_r has been reduced to fraction of a day
      mtimerh=mtimer_r/60
      mtimerm=mtimer_r-mtimerh*60  ! minutes left over
      ihr=ihr+mtimerh
      imins=imins+mtimerm
      if(imins==58.or.imins==59)then
!       allow for roundoff for real timer from old runs
        print *,'*** imins increased to 60 from imins = ',imins
        imins=60
      endif
      if(imins>59)then
        imins=imins-60
        ihr=ihr+1
      endif
      if(diag)print *,'d datefix iyr,imo,iday,ihr,imins,mtimer_r: ',
     .                   iyr,imo,iday,ihr,imins,mtimer_r
      if(ihr>23)then
        ihr=ihr-24
        iday=iday+1
      endif
      if(diag)print *,'e datefix iyr,imo,iday,ihr,imins,mtimer_r: ',
     .                   iyr,imo,iday,ihr,imins,mtimer_r

      if(iday>mdays(imo))then
         iday=iday-mdays(imo)
         imo=imo+1
         if(imo>12)then
           imo=imo-12
           iyr=iyr+1
         endif
      endif

      kdate_r=iday+100*(imo+100*iyr)
      ktime_r=ihr*100+imins
      mtimer=0
      if(diag)print *,'end datefix iyr,imo,iday,ihr,imins,mtimer_r: ',
     .                         iyr,imo,iday,ihr,imins,mtimer_r
      print *,'leaving datefix kdate_r,ktime_r: ',kdate_r,ktime_r

      return
      end
