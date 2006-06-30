      subroutine amipsst
      use cc_mpi
      implicit none
!     this one primarily does namip=2      
!     A routine for the con-cubic model which will interpolate
!     linearly in time between two sst data sets.
!     iday is a variable which is equal to the number of
!     day of the month.(Thus iday = 1 at the
!     start of a month, is updated to 2 at end of first day etc.)
!     imo is the number of the month we are in (1<=imo<=12)
      include 'newmpar.h'
      include 'arrays.h'    ! ts, t, u, v, psl, ps, zs
      include 'dates.h'     !  kdate,ktime,timer,mtimer
      include 'filnames.h'  ! list of files
      include 'latlong.h'
      include 'nsibd.h'     ! sib, tsigmf, isoilm
      include 'parm.h'      ! id,jd
      include 'pbl.h'       ! tss
      include 'permsurf.h'  ! iperm etc
      include 'soil.h'      ! ,tice, alb
      include 'soilsnow.h'  ! fracice,sicedep
      include 'mpif.h'
      real, save, dimension(ifull_g) ::  ssta_g, sstb_g, aice_g, bice_g
      real, save, dimension(ifull) :: ssta, sstb, aice, bice
      real fraciceb, dum2b
      common/work2b/fraciceb(ifull),dum2b(ifull,2)
      character*22 header
      integer, parameter, dimension(0:13) :: mdays =
     &     (/ 31, 31,28,31,30,31,30,31,31,30,31,30,31, 31 /)
      logical update
c     integer ipermp(ifull),indexi,indexs,ip 
      integer iyr, imo, iday, iyr_m, imo_m, imonth, iyear, idjd_g,
     &        iyr_p, imo_p, ierr, iq
      save iyr,imo,iday,iyr_m,imo_m
      real rat1, rat2
      integer il_in,jl_in
      real rlon_in,rlat_in,schmidt_in

!     rml 21/03/06 declarations for sst file as netcdf
      logical, save :: netcdf_file=.false.
      integer, save :: ncid,sstid,iceid,nkeep
      integer  :: timedim,yearid,monthid,ierr2,n,ntime
      integer  :: start(2),kount(2)
      integer, dimension(:), allocatable, save :: sstyr,sstmon
      include 'netcdf.inc'


      idjd_g = id + (jd-1)*il_g
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if ( myid == 0 ) then
      update = .false.

      print *,'amipsst  called'
      iyr=kdate/10000
      imo=(kdate-10000*iyr)/100
      iday=kdate-10000*iyr-100*imo  +mtimer/(60*24)
      do while (iday>mdays(imo))
       iday=iday-mdays(imo)
       imo=imo+1
       if(imo>12)then
         imo=1               ! no leap years for the moment
         iyr=iyr+1
       endif
      enddo

      if(ktau==0)then
        update = .true.
        iyr_m=iyr
        imo_m=imo-1
        if(imo_m==0)then
          imo_m=12
          iyr_m=iyr-1
        endif

!     rml 21/03/06 add option for reading from netcdf format file
!     try to open as netcdf to determine which read to use
        ierr2 = nf_open(sstfile,0,ncid)
        if (ierr2.eq.0) netcdf_file=.true.
        if (netcdf_file) then
          ierr2=nf_inq_dimid(ncid,'time',timedim)
          if (ierr2.ne.nf_noerr) stop 'time dimension error'
          ierr2=nf_inq_dimlen(ncid,timedim,ntime)
          ierr2=nf_inq_varid(ncid,'year',yearid)
          allocate(sstyr(ntime),sstmon(ntime))
          ierr2=nf_get_var_int(ncid,yearid,sstyr)
          ierr2=nf_inq_varid(ncid,'month',monthid)
          ierr2=nf_get_var_int(ncid,monthid,sstmon)
          ierr2=nf_inq_varid(ncid,'sst',sstid)
          if (ierr2.ne.nf_noerr) stop 'error finding sst variable'
!
!         find required records
          do n=1,ntime
            if (sstyr(n).eq.iyr_m.and.sstmon(n).eq.imo_m)
     &        nkeep=n
          enddo
!         read preceeding month
          start(1)=1; start(2)=nkeep
          kount(1)=il*jl; kount(2)=1
          write(6,*) 'sstcheck:',ncid,sstid,start,kount
          ierr2=nf_get_vara_real(ncid,sstid,start,kount,ssta_g)
          if (ierr2.ne.nf_noerr) stop 'error reading ssta_g'
          ssta_g = ssta_g + 273.16
c         read current month
          if (sstyr(nkeep+1).ne.iyr.or.sstmon(nkeep+1).ne.imo) stop
          start(2)=nkeep+1
          write(6,*) 'sstcheck2:',ncid,sstid,start,kount
          ierr2=nf_get_vara_real(ncid,sstid,start,kount,sstb_g)
          if (ierr2.ne.nf_noerr) stop 'error reading sstb_g'
          sstb_g = sstb_g + 273.16
          print *,'ssta(idjd) ',ssta_g(idjd_g)
          print *,'sstb(idjd) ',sstb_g(idjd_g)

        else
!         original code for reading sst
          open(unit=75,file=sstfile,status='old',form='formatted')
2         read(75,'(i2,i5,2i4,2f6.1,f6.3,a22)')
     &      imonth,iyear,il_in,jl_in,rlon_in,rlat_in,schmidt_in,header
          write(6,'("reading sst ",i2,i5,2i4,2f6.1,f6.3,a22)')
     &      imonth,iyear,il_in,jl_in,rlon_in,rlat_in,schmidt_in,header
!          if(il/=il_in.or.il/=il_in.or.rlong0/=rlon_in.or.
          if(il_g/=il_in.or.jl_g/=jl_in.or.rlong0/=rlon_in.or.
     &     rlat0/=rlat_in.or.schmidt/=schmidt_in)then
            write(0,*) 'wrong amipsst file'
            stop
          endif
          read(75,*) ssta_g
          ssta_g(:)=ssta_g(:)*.01 -50. +273.16
          print *,'want imo_m,iyr_m; ssta ',imo_m,iyr_m,ssta_g(idjd_g)
          if(iyr_m.ne.iyear.or.imo_m.ne.imonth)go to 2

          read(75,'(i2,i5,a22)') imonth,iyear,header
          print *,'reading sstb data:',imonth,iyear,header
          print *,'should agree with imo,iyr ',imo,iyr
          if(iyr.ne.iyear.or.imo.ne.imonth)stop
          read(75,*) sstb_g
          sstb_g(:)=sstb_g(:)*.01 -50. +273.16
          print *,'sstb(idjd) ',sstb_g(idjd_g)
        endif
 
        if (namip==2) then ! sea ice also read at middle of month
          if (netcdf_file) then
!           read ice - no need to check mon/yr because same index
!           read preceeding month
            ierr2=nf_inq_varid(ncid,'ice',iceid)
            if (ierr2.ne.nf_noerr) stop 'error finding ice variable'
            start(2)=nkeep
            write(6,*) 'icecheck:',ncid,iceid,start,kount
            ierr2=nf_get_vara_real(ncid,iceid,start,kount,aice_g)
            if (ierr2.ne.nf_noerr) stop 'error reading aice'
c           read current month
            start(2)=nkeep+1
            write(6,*) 'icecheck2:',ncid,iceid,start,kount
            ierr2=nf_get_vara_real(ncid,iceid,start,kount,bice_g)
            if (ierr2.ne.nf_noerr) stop 'error reading bice'

            print *,'aice(idjd) ',aice_g(idjd_g)
            print *,'bice(idjd) ',bice_g(idjd_g)

          else
!           original code for reading ice file
            open(unit=76,file=icefile,status='old',form='formatted') 
21          read(76,'(i2,i5,2i4,2f6.1,f6.3,a22)')
     &      imonth,iyear,il_in,jl_in,rlon_in,rlat_in,schmidt_in,header
            write(6,'("reading ice ",i2,i5,2i4,2f6.1,f6.3,a22)')
     &      imonth,iyear,il_in,jl_in,rlon_in,rlat_in,schmidt_in,header
!            if(il/=il_in.or.il/=il_in.or.rlong0/=rlon_in.or.
            if(il_g/=il_in.or.jl_g/=jl_in.or.rlong0/=rlon_in.or.
     &       rlat0/=rlat_in.or.schmidt/=schmidt_in)then
              write(0,*) 'wrong amipice file'
              stop
            endif
            read(76,*) aice_g
            print *,'want imo_m,iyr_m; aice ',imo_m,iyr_m,aice_g(idjd_g)
            if(iyr_m.ne.iyear.or.imo_m.ne.imonth)go to 21
            read(76,'(i2,i5,a22)') imonth,iyear,header
            print *,'reading b_sice data:',imonth,iyear,header
            print *,'should agree with imo,iyr ',imo,iyr
            if(iyr.ne.iyear.or.imo.ne.imonth)stop
            read(76,*) bice_g
            print *,'bice(idjd) ',bice_g(idjd_g)
          endif 
        endif   ! (namip==2) 
      endif     ! (ktau==0)

      if(iday==mdays(imo)/2)then
         update = .true.
         ssta_g(:)=sstb_g(:)    ! rotate data sets
         aice_g(:)=bice_g(:)    ! rotate data sets
         iyr_p=iyr
         imo_p=imo+1
         if(imo_p==13)then
            imo_p=1
            iyr_p=iyr+1
         endif                  ! (imo_p==13)
 
!       read in next months data
        if (netcdf_file) then
          write(6,*) 'reading new sst/ice data'
          if (sstyr(nkeep+2).ne.iyr_p.or.sstmon(nkeep+2).ne.imo_p) then
            write(6,*) 'date check: ',iyr_p,imo_p,sstyr(nkeep+2),
     &                  sstmon(nkeep+1)
            stop 'failed to match sst year/month'
          endif
          start(1)=1; start(2)=nkeep+2
          kount(1)=il*jl; kount(2)=1
          ierr2=nf_get_vara_real(ncid,sstid,start,kount,sstb_g)
          if (ierr2.ne.nf_noerr) stop 'error reading sstb'
          sstb_g = sstb_g + 273.16
          if (namip.eq.2) then
            start(2)=nkeep+2
            ierr2=nf_get_vara_real(ncid,iceid,start,kount,bice_g)
            if (ierr2.ne.nf_noerr) stop 'error reading bice'
          endif

        else
3         read(75,'(i2,i5,a22)') imonth,iyear,header
          print *,'reading sstb data:',imonth,iyear,header
          print *,'comparing with imo_p,iyr_p ',imo_p,iyr_p
          read(75,*) sstb_g
          sstb_g(:)=sstb_g(:)*.01 -50. +273.16
          if(iyr_p.ne.iyear.or.imo_p.ne.imonth)go to 3

          if(namip==2)then   ! sice also read at middle of month
31          read(76,'(i2,i5,a22)') imonth,iyear,header
            print *,'reading b_sice data:',imonth,iyear,header
            print *,'comparing with imo_p,iyr_p ',imo_p,iyr_p
            read(76,*) bice_g
            if(iyr_p.ne.iyear.or.imo_p.ne.imonth)go to 31
          endif  ! (namip==2)
        endif
      endif    ! (iday==mdays(imo)/2)

      if(iday<mdays(imo)/2)then  ! 1st half of month
        rat1=(mdays(imo)-2.*iday)/(mdays(imo)+mdays(imo-1))
        rat2=(2.*iday+mdays(imo-1))/(mdays(imo)+mdays(imo-1))
      else                             ! 2nd half of month
        rat1=(mdays(imo+1)+2.*mdays(imo)-2.*iday)/
     .                               (mdays(imo+1)+mdays(imo))
        rat2=(2.*iday-mdays(imo))/(mdays(imo+1)+mdays(imo))
      endif
      print *,'month_imo,iday',imo,iday
      if(iday>1.or.namip==2)go to 6

!     amipice for namip=1   N.B. used below to mask the SSTs
      print*, "NAMIP=1 not implemented in MPI version"
!     This would be easy to do if required.
      stop

      end if ! myid==0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

!     Each day interpolate non-land sst's
6     continue

!     Could calculate these in each process but this makes the code simpler
      call MPI_Bcast(rat1,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(rat2,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(update,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      if ( update ) then
         if ( myid == 0 ) then
            call ccmpi_distribute(ssta, ssta_g)
            call ccmpi_distribute(sstb, sstb_g)
            call ccmpi_distribute(aice, aice_g)
            call ccmpi_distribute(bice, bice_g)
         else
            call ccmpi_distribute(ssta)
            call ccmpi_distribute(sstb)
            call ccmpi_distribute(aice)
            call ccmpi_distribute(bice)
         end if
      end if

      do iq=1,ifull  
       if(.not.land(iq))then
         tgg(iq,1)=rat1*ssta(iq)+rat2*sstb(iq)  ! sea water temperature
c        alb(iq)=.11       ! set in indata/radriv90
       endif      ! (.not.land(iq))
      enddo
      if ( mydiag ) print *,'rat1,rat2,land,ssta,sstb,tgg1: ',
     %           rat1,rat2,land(idjd),ssta(idjd),sstb(idjd),tgg(idjd,1)
      fraciceb(:)=min(.01*(rat1*aice(:)+rat2*bice(:)),1.) ! convert from %
      do iq=1,ifull
        if(fraciceb(iq)<=.02)fraciceb(iq)=0.
      enddo
      
      sicedep(:)=0. 
      if(ktau==0)then  ! will set sicedep in indata
        fracice(:)=fraciceb(:)
        do iq=1,ifull
         if(.not.land(iq))then
           tss(iq)=tgg(iq,1)
         endif    ! (.not.land(iq))
        enddo
        if (mydiag) print *,'0rat1,rat2,sicedep,aice,bice,fracice: ',
     &      rat1,rat2,sicedep(idjd),aice(idjd),bice(idjd),fracice(idjd)
        return
      endif       ! (ktau==0)
      
      do iq=1,ifull
       if(.not.land(iq))then
         if(fraciceb(iq)>0.)then
           if(fracice(iq)==0.)then
!            create values for tice, and set averaged tss
!            N.B. if already a sice point, keep present tice
             tgg(iq,3)=min(271.2,tss(iq),t(iq,1)+.04*6.5) ! for 40 m lev1
           endif  ! (fracice(iq)==0.)
           if(rlatt(iq)>0.)then
             sicedep(iq)=2.
           else
             sicedep(iq)=1.
           endif ! (rlatt(iq)>0.)
         endif    ! (fraciceb(iq)>0.)
         fracice(iq)=fraciceb(iq)
         tss(iq)=tgg(iq,3)*fracice(iq)+tgg(iq,1)*(1.-fracice(iq))
       endif      ! (.not.land(iq))
      enddo
      if (mydiag) print *,'rat1,rat2,sicedep,aice,bice,fracice: ',
     &      rat1,rat2,sicedep(idjd),aice(idjd),bice(idjd),fracice(idjd)
c!     update surface permutation array and sice
c      ipermp(:)=iperm(:)
c      indexi=ipland
c      indexs=ipsea+1
c!cdir nodep
c      do ip=ipland+1,ipsea
c       iq=ipermp(ip)
c       if(fracice(iq)>0.)then
c         indexi=indexi+1     ! sice point
c         iperm(indexi)=iq    ! sice point
c         sice(iq)=.true.
c       else
c         indexs=indexs-1     ! sea point
c         iperm(indexs)=iq    ! sea point
c         sice(iq)=.false.
c       endif  ! (fracice(iq)>0.)
c      enddo   !  ip loop
c      ipsice=indexi
c      if ( myid == 0 ) then
c         print *,'ktau,ipland,ipsice,ipsea,sice update in amipsst: ',
c     .            ktau,ipland,ipsice,ipsea,sicedep(idjd)
c      end if
      return
      end
