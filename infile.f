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
      use cc_mpi
      use diag_m
      implicit none
      include 'newmpar.h'
      include 'extraout.h'
      include 'parm.h'
      include 'parm_nqg.h'  ! nqg_r,nqg_set
      include 'screen.h'
      include 'stime.h'  ! kdate_s,ktime_s  sought values for data read
      include 'tracers.h'
      include 'mpif.h'
      include 'dates.h'
      integer io_in2, kdate_r, ktime_r, nem2, nested
      real timeg_r, ds_r
      real sigin
      integer kk
      common/sigin/sigin(kl),kk  ! for vertint, infile
      real rlong0x,rlat0x,schmidtx
      common/schmidtx/rlong0x,rlat0x,schmidtx ! infile, newin, nestin, indata
!!     don't use dum0 as shared with zss      
!      common/work2/dum0(ifull),dum1(ifull),dum2(ifull),dum3(ifull),
!     .             dum4(ifull,14)
      real psl(ifull),pmsl(ifull),zs(ifull),em(ifull),f(ifull),
     . tss(ifull),precip(ifull),wb(ifull,ms),wbice(ifull,ms),
     . snowd(ifull),alb(ifull),sicedep(ifull),
     . t(ifull,kl),u(ifull,kl),v(ifull,kl),qg(ifull,kl),
     . tgg(ifull,ms),tggsn(ifull,3),smass(ifull,3),ssdn(ifull,3),
     . ssdnn(ifull),osnowd(ifull),snage(ifull)
      integer isflag(ifull)
      integer ktau_r, ibb, jbb, i, ii, jj, k, iq, mtimer_in
      integer ik, jk, m2, nsd2, nx1, nps2, mex2, mup2, mtimer_r, nmi2,
     &        ndt2, npsav2, nhor2, khdif2, kwt2, nx3, nx4, nvad2,
     &        nunst2, nrun2, nrunx2, ndum, nn, ntsea2, ms2, nkuo2,
     &        nextras, ntrac2, nij, nijk, imax, id2, jmax, jd2, n, ntr
      real timer_r, difknbd, rhkuo, du, tanl, dumpsa, dumpsm
      integer :: ierr

      character rundat2*10
      integer, save :: ncalled=0
      integer, save :: nqg_setin
      logical :: skip

      if(ncalled.eq.0)then
        nqg_setin=nqg_set   ! to allow for reading nested first
      endif
      if (myid==0)
     &    print *,'ncalled,nqg_setin,nqg_set ',ncalled,nqg_setin,nqg_set
      ncalled=ncalled+1   ! control diag prints
      mtimer_in = 0

      if(abs(io_in).eq.3)then   
 3      continue
        if (myid == 0 ) then
          read(10,end=8) kdate_r,ktime_r,ktau_r,ik,jk,kk,m2,nsd2,
     &         io_in2,nx1,nps2,mex2,mup2,nem2,mtimer_r,nmi2,ndt2,
     &         npsav2,rundat2,nhor2,nkuo2,khdif2,kwt2,nx3,nx4,
     &         timer_r,timeg_r,ds_r,nvad2,nqg_r,nunst2,nrun2,nrunx2,
     &         (ndum,nn=1,8),ntsur2,(ndum,nn=1,3),ntsea2,ms2,nextras,
     &         ndum,ntrac2,difknbd,rhkuo,du,tanl,rlong0x,rlat0x,schmidtx
          print *
          print *,'in:',kdate_r,ktime_r,ktau_r,ik,jk,kk,m2,nsd2,io_in2,
     &                 nps2,mex2,mup2,nem2,mtimer_r,
     &                 nmi2,ndt2,npsav2,rundat2,nhor2,nkuo2,khdif2,kwt2,
     &                 nx3,nx4,timer_r,timeg_r,ds_r,nvad2,nqg_r,nunst2,
     &                 nrun2,nrunx2,ms2
          print *,'values read in for timer_r, mtimer_r: ',
     &          timer_r,mtimer_r
          if(kk.gt.kwt2)stop 'kk.gt.kwt2'
          if(mtimer_r.eq.10)mtimer_r=0 ! till jjk fixes up input file 
          if(mtimer_r.ne.0)then ! preferred
!           assume mtimer_r is read in correctly and set timer_r from that
            timer_r=mtimer_r/60.
          else  ! mtimer_r = 0
            mtimer_r=nint(60.*timer_r)
          endif
          print *,'giving timer_r, mtimer_r: ',timer_r,mtimer_r
          mtimer_in=mtimer_r
          call datefix(kdate_r,ktime_r,mtimer_r) ! for Y2K, or mtimer_r>0
          print *,'after datefix kdate_r,ktime_r,mtimer_r: ',
     &                       kdate_r,ktime_r,mtimer_r
          nij=ik*jk
          nijk=nij*kk
          imax=min(ik,il)
          id2=max(ik-il,0)
          jmax=min(jk,jl)
          jd2=max(jk-jl,0)
          read(10) (sigin(i),i=1,kk)
          nqg_set=min(nqg_r,nqg_setin)
          if(ncalled.lt.4)then
            print *,'ds_r,nqg_r,nqg_set,mtimer_r:',
     &            ds_r,nqg_r,nqg_set,mtimer_r
            print *,'ik,imax,id2 ',ik,imax,id2
            print *,'jk,jmax,jd2 ',jk,jmax,jd2
            print *,'rlong0x,rlat0x,schmidtx',rlong0x,rlat0x,schmidtx
            print *,'sig in: ',(sigin(i),i=1,kk)
          endif                  ! (ncalled.lt.4)


          print *,'kdate_r,ktime_r: ',kdate_r,ktime_r,
     .            'kdate_s,ktime_s >= ',kdate_s,ktime_s
          if( kdate_s.gt.kdate_r.or.
     .      (kdate_s.eq.kdate_r.and.ktime_s.gt.ktime_r)) then
             skip = .true.
          else
             skip = .false.
          end if
          if(ncalled.eq.1)then
!       insist on exactly correct date for first call to infile
             if(kdate_s.ne.kdate_r.or.ktime_s.ne.ktime_r)then
                print *,'required date not on ifile'
                stop
             endif
          endif                 ! (ncalled.eq.1)
        endif  ! myid == 0 

        call MPI_Bcast(nqg_r,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(ntrac2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(skip,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(kk,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(sigin,kk,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(kdate,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(ktime,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(ktau_r,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(mtimer,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!       call MPI_Bcast(timer,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
!       call MPI_Bcast(timeg,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)

        call readglobvar(10,psl,skip)
        call readglobvar(10,pmsl,skip)
c     this is actually pmsl for io_in=3
        call readglobvar(10,zs,skip)
c     following assumes data written with nem=1 option
        call readglobvar(10,em,skip)
        call readglobvar(10,f,skip)
        call readglobvar(10,tss,skip)
        call readglobvar(10,precip,skip)
        if(ncalled.lt.4.and.mydiag.and..not.skip)then
          write (6,"('100*psl# in',3f7.2,1x,3f7.2,1x,3f7.2)") 
     .              100.*diagvals(psl)
          print *,'pmsl# in: ', diagvals(pmsl)
          write (6,"('zs# in  ',3f7.1,1x,3f7.1,1x,3f7.1)") 
     .              diagvals(zs)
          print *,'em# in: ', diagvals(em)
          print *,'f# in: ', diagvals(f)
          write (6,"('tss# in ',3f7.1,1x,3f7.1,1x,3f7.1)") 
     .              diagvals(tss)
          write (6,"('prec# in',3f7.2,1x,3f7.2,1x,3f7.2)") 
     .         diagvals(precip)
        endif
c
c     soil temp and moisture
        if(nqg_r.ge.8)then
          do k=1,ms
            call readglobvar(10,tgg(:,k),skip)
          end do
          do k=1,ms
            call readglobvar(10,wb(:,k),skip)
          end do
        else      ! from old write, meant for nsib=0 or 1
c         read subsoil temps in this case, lower one first
          call readglobvar(10,tgg(:,ms),skip)
          call readglobvar(10,tgg(:,2),skip)
c         read subsoil moistures, upper one first
          call readglobvar(10,wb(:,1),skip)
          call readglobvar(10,wb(:,ms),skip)
          tgg(:,1) = abs(tss(:))
          wb(:,2)  = wb(:,1)    ! layer initialisation of moisture
          do k=3,ms-1
            wb(:,k)  = wb(:,ms) ! layer initialisation of moisture
          enddo
          tgg(:,3) = tgg(:,2)   ! initial temper. from GCM runs with 3 layers
          do k=4,ms-1
            tgg(:,k) = tgg(:,ms) ! initial temper. from GCM runs with 3 layers
          enddo
        endif    ! (nqg_r.ge.8)  .. else ..
        if(ncalled.lt.4 .and. mydiag .and. .not.skip)then
          print *,'infile '
          print *,'tgg(1)# ', diagvals(tgg(:,1))
          print *,'tgg(2)# ', diagvals(tgg(:,2))
          print *,'tgg(ms)# ', diagvals(tgg(:,ms))
          print *,'wb(1)# ', diagvals(wb(:,1))
          print *,'wb(ms)# ', diagvals(wb(:,ms))
        endif
        if(nqg_r.ge.4) then     ! N.B. some ifile's have ngq=3
          call readglobvar(10,alb,skip)
          if(ncalled.lt.4.and. mydiag .and. .not.skip)then
            print *,'alb# ', diagvals(alb)
          endif
        endif
        if(nqg_r.ge.5)then
c       precc just read in to dummy array
          call readglobvar(10,precip,skip)
          if(ncalled.lt.4 .and. mydiag .and. .not.skip)then
            print *,'precc# ', diagvals(precip)
          endif
        endif
        if(nqg_r.ge.6)then
c       snowd read in to snowd array
          call readglobvar(10, snowd, skip)
          if ( myid == 0 ) then
c       N.B. these cloudlo, cloudmi, cloudhi not used because recalc. in radn
            read(10)           ! cloudlo
            read(10)           ! cloudmi
            read(10)           ! cloudhi
          end if

          call readglobvar(10,sicedep,skip) ! these were for nqg>=7
          if ( myid == 0 ) then
            read(10)           ! tscrn
            read(10)           ! qgscrn
            read(10)           ! u10
          endif
          if(ncalled.lt.4 .and. mydiag .and. .not.skip)then
            write (6,"('snowd# in',9f7.2)") diagvals(snowd)
            write (6,"('sicedep# in',9f7.2)") diagvals(sicedep)
          endif
        endif   ! (nqg_r.ge.6)
        if(nqg_set.lt.6)then
c         when no snowd available initially, e.g. COMPARE III (jlm formula)
!         if(tss(iq).lt.270.)snowd(iq)=max(50.,snowd(iq))
!         for this one,snow lin. increases from 5 to 55, for T 270 down to 260
          where ( tss < 270. )
            snowd = min(55.,5.*(271.-abs(tss)))   
          endwhere
          if(ncalled.lt.4 .and. mydiag .and. .not.skip)then
            print *,'setting snowd in infile, because nqg_set = ',
     &            nqg_set  
            write (6,"('snowd# preset to',9f7.2)") diagvals(snowd)
          endif
        endif                    ! (nqg_set.lt.6)
        if(nqg_r.ge.9.and.myid==0)then ! eg, fg, taux, tauy, runoff
          read(10)
          read(10)
          read(10)
          read(10)
          read(10)
        endif
        if(nqg_r.ge.10.and.myid==0)then ! tmaxscr, tminscr, tscr_ave
          read(10)
          read(10)
          read(10)
        endif
        if(nqg_r.ge.11)then
          call readglobvar(10,wbice,skip)
        endif
        if(nqg_r.ge.12)then
          print *,
     &        'reading snow variables: tggsn,smass,ssdn,snage,isflag'
          call readglobvar(10,tggsn,skip)
          if ( mydiag .and. .not.skip )
     &       print *,'tggsn(1)# ', diagvals(tggsn(:,1))
          call readglobvar(10,smass,skip)
          call readglobvar(10,ssdn,skip)
          call readglobvar(10,ssdnn,skip)
          call readglobvar(10,osnowd,skip) ! not saved for netcdf
          call readglobvar(10,snage,skip)
          if ( myid == 0 ) print *,'o kdate_r,ktime_r,mtimer_r: ',
     &         kdate_r,ktime_r,mtimer_r
          call readglobvar(10,isflag,skip)
          if ( mydiag .and. .not.skip )
     &       print *,'isflag(1)# ', diagvals(isflag)
        endif ! nqg = 12
        if(nqg_r.ge.13.and.myid==0)then
          read(10) ! rtsave
          read(10) ! tscrn3hr
        endif

        call readglobvar(10,t,skip)
        call readglobvar(10,u,skip)
        call readglobvar(10,v,skip)
        call readglobvar(10,qg,skip)
        if ( mydiag .and. .not.skip ) then
          print *,'t in ',t(idjd,:)
          print *,'u in ',u(idjd,:)
          print *,'qg in ',qg(idjd,:)
        end if
        if(npsav2.gt.0 .and. myid == 0)then
           read(10) (dumpsa,dumpsm,n=1,npsav2)
        endif
        if(ncalled.lt.4 .and. mydiag .and. .not.skip)then
          print *,'t in ',t(idjd,:)
          print *,'u in ',u(idjd,:)
          print *,'v in ',v(idjd,:)
          print *,'qg in ',qg(idjd,:)
!          print *,'last psa,psm ',dumpsa,dumpsm
!          print *,'nsd2,ntrac2 ',nsd2,ntrac2
        endif  ! (ncalled.lt.4)

        if(nsd2.eq.1.and.myid==0)then    ! Thu  04-15-1993
          if(ncalled.lt.4)print *,'about to read sdot'
          read(10) 
        endif 
        if(ntrac2.gt.1)then
          do ntr=1,ntrac2
            if(ntrac.eq.ntrac2)then
              if(ncalled.lt.4.and.myid==0)
     &              print *,'read tr for ntr= ',ntr
              call readglobvar(10,tr(:,:,ntr),skip) ! ((tr(iq,k,ntr),iq=1,ilt*jlt),k=1,klt)
            else
              if(ncalled.lt.4.and.myid==0) 
     &              print *,'dummy read tr for ntr= ',ntr
              if (myid == 0 ) read(10)
            endif  ! (ntrac.eq.ntrac2)
          enddo
        endif

      endif  ! (abs(io_in).eq.3)

      if(abs(io_in).eq.1)then   ! for netcdf
        call newin(io_in2,kdate_r,ktime_r,ktau_r,ibb,jbb,nem2,
     &   timeg,ds_r,psl,zs,tss,
     &   t(1:ifull,:),u(1:ifull,:),v(1:ifull,:),qg(1:ifull,:),nested)
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
        endif  ! (ncalled.lt.4)
      endif ! ( abs(io_in).eq.1 ) then

!     reset kdate_s and ktime_s ready for next (nesting) call to infile
      kdate_s=kdate_r
      ktime_s=ktime_r+1 

      if(qg(ifull/2,1).gt..1.or.qg(ifull,1).gt..1)then
        do k=1,kk
         do iq=1,ifull
            qg(iq,k)=.001*qg(iq,k)
         enddo
        enddo
        if(ncalled.lt.4.and.mydiag)then
          print *,'**** old-style file so rescaling qg from g/kg to g/g'
          print *,'ijk, qg becomes ',ijk,(qg(idjd,k),k=1,kk)
        endif
      endif     ! (qg(ifull/2,1).gt..1)
      do k=1,kk
       do iq=1,ifull
         qg(iq,k)=max(qg(iq,k),1.e-6)
       enddo
      enddo

      if(mtimer_in.gt.0.and.nrungcm.eq.2)then
c       nrungcm=0
        print *,'*** re-setting NRUNGCM nrungcm from 2 to 0 because ',
     .          'not starting from gcm file'
        stop 'wrong nrungcm used'
      endif

      if ( mydiag ) then
         print *,'end infile; next read will be kdate_s,ktime_s >= ',
     &                                       kdate_s,ktime_s
      end if
      return
8     print *,'end of data reading file in infile'
      stop
      end

      subroutine newin(io_in2,kdate,ktime,ktau_r,ibb,jbb,nem2,
     & timeg,ds_r,psl,zs,tss,
     & t,u,v,qg,nested)
      use cc_mpi
      use diag_m
      implicit none
!     nested=0  for calls from indata; 1  for calls from nestin     
      integer, parameter :: nmeth=1 ! 0 old way, 1 preferred new way (see nestin too)
      include 'newmpar.h'
      include 'darcdf.h'
      include 'extraout.h'
      include 'kuocom.h'    ! ldr
      include 'liqwpar.h'  ! ifullw
      include 'mapproj.h'   ! du,tanl,rnml,stl1,stl2
      include 'morepbl.h'
      include 'netcdf.inc'
      include 'parm.h'
      include 'parm_nqg.h'  ! nqg_r,nqg_set
      include 'screen.h'
      include 'soil.h'      !  alb, sicedep
      include 'soilsnow.h'  ! with nested=0/1 option
      include 'stime.h'     ! kdate_s,ktime_s  sought values for data read
      include 'tracers.h'
      include 'mpif.h'
      integer io_in2,kdate,ktime,ktau_r,ibb,jbb,nem2,nested
      real timeg,ds_r

      real psl(ifull),zs(ifull),tss(ifull),
     & t(ifull,kl),u(ifull,kl),v(ifull,kl),qg(ifull,kl)

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
      integer istrt,iend,jstrt,jend,kstrt,kend,numi,numj,numk
      common/partgrid/istrt,iend,jstrt,jend,kstrt,kend,numi,numj,numk
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
c     sdu=du
c     stanl=tanl
c     srnml=rnml
c     sstl1=stl1
c     sstl2=stl2

      if(dss.gt.0..and.npanels.eq.0)then   ! for later DARLAM option
         print*,
     &   "Error, npanels=0 DARLAM option not implemented in MPI version"
         stop
c       find lat/lon of four corner points
!!!        call lconset(ds_r)
!!!        call lconll(rlon11,rlat11,float(1),float(1))
!!!        call lconll(rlon12,rlat12,float(1),float(jl))
!!!        call lconll(rlon21,rlat21,float(il),float(1))
!!!        call lconll(rlon22,rlat22,float(il),float(jl))
!!!        if(ktau.le.1)then
!!!          write(6,'("model dss,du,tanl,rnml,stl1,stl2= ",6f10.2)')
!!!     &                 dss,du,tanl,rnml,stl1,stl2
!!!          write(6,'("rlon11,12,21,22= ",4f10.2)')
!!!     &             rlon11,rlon12,rlon21,rlon22
!!!          write(6,'("rlat11,12,21,22= ",4f10.2)')
!!!     &             rlat11,rlat12,rlat21,rlat22
!!!        endif
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
      
!!!      if(npanels.eq.0)then  ! just for DARLAM
!!!c       set up input map projs.
!!!        ds_r=tds
!!!        du=tdu
!!!        tanl=ttanl
!!!        rnml=trnml
!!!        stl1=tstl1
!!!        stl2=tstl2
!!!        if(ktau.eq.1)
!!!     &    write(6,'("newin ds_r,du,tanl,rnml,stl1,stl2=",6f10.2)')
!!!     &                     ds_r,du,tanl,rnml,stl1,stl2
!!!
!!!c       set up lconset
!!!        call lconset(ds_r)
!!!      endif  ! (npanels.eq.0)
      ik=nahead(1)
      jk=nahead(2)
      kk=nahead(3)
      if(ktau.le.1.and.myid==0)
     &     print *,'in newin; ktau,ik,jk,kk=', ktau,ik,jk,kk

!!!      if(ds_r.gt.0..and.npanels.eq.0)then
!!!c       find i/j in input grid of model grid
!!!        call lconij(rlon11,rlat11,ri11,rj11,theta)
!!!        if(ktau.eq.1)print *,'(1,1)',ri11,rj11
!!!        call lconij(rlon21,rlat21,ri21,rj21,theta)
!!!        if(ktau.eq.1)print *,'(il,1)',ri21,rj21
!!!        call lconij(rlon12,rlat12,ri12,rj12,theta)
!!!        if(ktau.eq.1)print *,'(1,jl)',ri12,rj12
!!!        call lconij(rlon22,rlat22,ri22,rj22,theta)
!!!        if(ktau.eq.1)print *,'(ifull)',ri22,rj22
!!!        xgn=min(ri11,ri12,ri21,ri22)
!!!        ygn=min(rj11,rj12,rj21,rj22)
!!!        xgx=max(ri11,ri12,ri21,ri22)
!!!        ygx=max(rj11,rj12,rj21,rj22)
!!!        if(ktau.le.1)
!!!     &   write(6,'("xgn,ygn,xgx,ygx=",4f10.2)')xgn,ygn,xgx,ygx
!!!
!!!        if(il.ge.ik)then
!!!          istrt=1
!!!          iend=ik
!!!          numi=ik
!!!        else!(il.lt.ik)then
!!!          numi=int(xgx)-int(xgn)+1
!!!          idif=max(1,(il-numi)/2)
!!!          istrt=max(1,int(xgn)-idif)
!!!          iend=istrt+il-1
!!!          numi=il
!!!          xgn=xgn+1.-float(istrt)
!!!          if(iend.lt.int(xgx)+1) stop 'newin iend.lt.int(xgx)+1'
!!!        endif!(il.ge.ik)then
!!!
!!!        if(jl.ge.jk)then
!!!          jstrt=1
!!!          jend=jk
!!!          numj=jk
!!!        else!(jl.lt.jk)then
!!!          numj=int(ygx)-int(ygn)+1
!!!          jdif=max(1,(jl-numj)/2)
!!!          jstrt=max(1,int(ygn)-jdif)
!!!          jend=jstrt+jl-1
!!!          numj=jl
!!!          ygn=ygn+1.-float(jstrt)
!!!          if(jend.lt.int(ygx)+1) stop 'newin jend.lt.int(ygx)+1'
!!!        endif!(il.ge.ik)then
!!!          
!!!      else      ! for C-C just the following (full arrays)
        istrt=1
        iend=il_g
        jstrt=1
        jend=jl_g
        xgn=1.
        ygn=1.
        numi=il_g
        numj=jl_g
!!!      endif!(ds_r.gt.0..and.npanels.eq.0)  .. else ..

      if(kk.eq.kl)then
        kstrt=1
        kend =kk
        numk=kl
      else
        if(kk.gt.kl) stop 'kk gt kl in newin'
        kstrt=1
        kend =kk
        numk=kk
      endif

      if(ktau.le.1.and.mydiag)then
        write(6,'("istrt,iend,jstrt,jend,kstrt,kend=",6i7)')
     &           istrt,iend,jstrt,jend,kstrt,kend
        write(6,'("numi,numj,numk=",3i7)')
     &           numi,numj,numk
        write(6,'("xgn,ygn=",2f7.2)')  xgn,ygn
      endif
      if ( numi .gt. il_g ) stop 'numi > il'
      if ( numj .gt. jl_g ) stop 'numj > jl'
      if ( numk .gt. kl ) stop 'numk > kl'

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
           call ncvgt1(ncid,idv,timeco,kdate,ier)
           timeco(1)=2
           timeco(2)=iarch
           call ncvgt1(ncid,idv,timeco,ktime,ier)
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
           call ncvgt1(ncid,idv,timeco,timeg,ier)
        else                    ! ier.ne.0
!     new style
           idv = ncvid(ncid,'kdate',ier)
           call ncvgt1(ncid,idv,iarch,kdate,ier)
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
           call ncvgt1(ncid,idv,iarch,ktime,ier)
           idv = ncvid(ncid,'ktau',ier)
           call ncvgt1(ncid,idv,iarch,ktau_r,ier)
           idv = ncvid(ncid,'timeg',ier)
           call ncvgt1(ncid,idv,iarch,timeg,ier)
        endif                   ! ier = 0
        if(kdate.eq.0)then
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
           kdate=idy+100*(imo+100*iyr)
           ktime=ihr*100
        endif                   ! (kdate.eq.0)
        print *,'kdate,ktime,ktau,timer_in: ',kdate,ktime,ktau,timer
!       fix up timer in 140-year run for roundoff (jlm Mon  01-11-1999)
        if(timer.gt.1048583.)then ! applies year 440 mid-month 9
           timer=timer+.125
        elseif(timer.gt.1048247.)then ! applies year 440 month 9
           timer=timer+.0625
        endif
        print *,'dtin,ktau_r: ',dtin,ktau_r
        print *,'kdate,ktime: ',kdate,ktime,
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
        call datefix(kdate,ktime,mtimer) ! for Y2K, or mtimer>0
        if (2400*kdate+ktime.lt.2400*kdate_s+ktime_s)go to 19
        
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
      call MPI_Bcast(kdate,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(ktime,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(ktau_r,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(mtimer,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(timer,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(timeg,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)

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
          do iq=1,ifull
           if(tss(iq).le.271.2)then
             sicedep(iq)=2. 
           else
             sicedep(iq)=0.     
           endif
        enddo   ! iq loop
	 endif

         if ( mydiag ) then
            print *,'at end of newin kdate,ktime,ktau,tss: ',
     &                     kdate,ktime,ktau,tss(idjd)
            print *,'tgg ',(tgg(idjd,k),k=1,ms)
            print *,'wb ',(wb(idjd,k),k=1,ms)
            print *,'wbice ',(wbice(idjd,k),k=1,ms)
         end if
      endif  ! (nested.eq.0)   !  only done at start of run 

      iarch=iarch+1

      if ( mydiag ) then
         write(6,'("end newin kdate,ktime,ktau_r,mtimer,timer,ds_r="
     &      ,4i10,2f10.1)')kdate,ktime,ktau_r,mtimer,timer,ds_r
      end if

      return ! newin
      end
c***************************************************************************
      subroutine histrd1(idnc,iarch,ier,name,ik,jk,var)
      use cc_mpi
      implicit none
      include 'newmpar.h'
      include 'parm.h'
      include 'netcdf.inc'
      include 'mpif.h'
      integer idnc, iarch, ier, ik, jk
      integer istrt,iend,jstrt,jend,kstrt,kend,numi,numj,numk
      common/partgrid/istrt,iend,jstrt,jend,kstrt,kend,numi,numj,numk
!     common/work3/dont_use(ifull,kl,3),ivar(ijk),spare(ijk)
      integer*2 ivar(200*200)
      logical odiag
      parameter ( odiag=.false. )
      character name*(*)
      integer start(3),count(3)
      real var(ifull)
      real globvar(ifull_g), vmax, vmin, addoff, sf
      integer ierr, ni, nj, idv, i, iq, inij, j, jin, iin

      if (myid==0) then
         if ( ik*jk.gt.200*200) stop 'histrd1 ik*jk gt 200*200'
         ni=iend-istrt+1
         nj=jend-jstrt+1
         if ( ni*nj.gt.il_g*jl_g) stop 'histrd1 ni*nj gt il*jl'

         start = (/ 1, 1, iarch /)
         count = (/ ik, jk, 1 /)

c     get variable idv
         idv = ncvid(idnc,name,ier)
c     write(6,*)'istrt  iend jstrt  jend kstrt kend numi numj numk '
c    &      ,' ni   nj   ier'
c     write(6,'(4i6,8i5)')istrt,iend,jstrt,jend,kstrt,kend,numi,numj
c    &               ,numk,ni,nj,ier
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
c     saving only data needed
            vmax=-1.e29
            vmin= 1.e29
            do jin=jstrt,jend
               j=jin+1-jstrt
               do iin=istrt,iend
                  inij=iin+(jin-1)*ik
                  i=iin+1-istrt
                  if(odiag)write(6,*)iin,jin,i,j,inij,ivar(inij)
                  iq = i + (j-1)*il_g
                  globvar(iq) = ivar(inij)*sf + addoff
                  vmax=max(vmax,globvar(iq))
                  vmin=min(vmin,globvar(iq))
               end do
            end do

            if(mod(ktau,nmaxpr).eq.0.or.odiag)
     &           write(6,'("done histrd1 ",a6,i4,i3,3e14.6)')
     &           name,ier,iarch,vmin,vmax,globvar(id+(jd-1)*il_g)
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
      integer istrt,iend,jstrt,jend,kstrt,kend,numi,numj,numk
      common/partgrid/istrt,iend,jstrt,jend,kstrt,kend,numi,numj,numk
!     common/work3/dont_use(ifull,kl,3),ivar(ijk),spare(ijk)
      integer*2 ivar(200*200*35)
      character name*(*)
      integer start(4),count(4)
      real var(ifull,kl)
      real globvar(ifull_g,kl), vmax, vmin, addoff, sf
      integer ierr, ni, nj, idv, i, iq, inij, j, jin, k, iin, inijk

      if (myid == 0 ) then
         if ( ik*jk*kk.gt.200*200*35)
     &        stop 'histrd4 ik*jk*kk gt 200*200*35'
!     if ( ik*jk*kk.gt.ijk) stop 'histrd4 ik*jk*kk gt ijk'
         ni=iend-istrt+1
         nj=jend-jstrt+1
         if ( ni*nj*numk.gt.il_g*jl_g*kl)
     &        stop 'histrd4 ni*nj*numk gt il*jl*kl'

         start = (/ 1, 1, 1, iarch /)
         count = (/ il_g, jl_g, kk, 1 /)

c     get variable idv
         idv = ncvid(idnc,name,ier)
         write(6,*)'istrt  iend jstrt  jend kstrt kend numi numj numk '
     &      ,' ni   nj  ier'
         write(6,'(4i6,8i5)')istrt,iend,jstrt,jend,kstrt,kend,numi,numj
     &               ,numk,ni,nj,ier
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
            vmax=-1.e29
            vmin= 1.e29
            do k=1,numk
               do jin=jstrt,jend
                  do iin=istrt,iend
                     inijk=iin+(jin-1)*ik+(k-1)*ik*jk
                     i=iin+1-istrt
                     j=jin+1-jstrt
                     iq = i + (j-1)*il_g
                     globvar(iq,k) = ivar(inijk)*sf + addoff
                     vmax=max(vmax,globvar(iq,k))
                     vmin=min(vmin,globvar(iq,k))
                  end do        ! i
               end do           ! j
            end do              ! k=1,numk

            if(mod(ktau,nmaxpr).eq.0)
     &           write(6,'("done histrd4 ",a6,i4,i3,3f12.4)') 
     &           name,ier,iarch,vmin,vmax,globvar(id+(jd-1)*il_g,nlv)
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
