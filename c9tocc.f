c need:    include '/usr/local/include/netcdf.inc'
c also need to set up land-sea masks
c this is cctoc9
c N.B. by setting n1x=1, only write out for ndayx and ntimex
c e.g  does only day 2  first time period, by setting ndayx=2, ntimex=1
c It is related to c9tolam4, but writes data for the conformal-cubic grid
c Winds are done on the unstaggered points,
c but are staggered before being written  (NOT from Wed  10-15-1997)
c An important point is that latitudes and longitudes are supplied in
c arrays rlat and rlong from setxyz, and are converted to degrees
c c9tolam4 combined c9tolam3.f, lintrp15.f and sstamip.f (kcn). It
c extracts data from GCM, horizontally interpolate all fields
c and also interpolates sea surface temperature.
c Instead of extracting data from GCM then write to gaussgrid.dat for
c lintrp15.f and sstamip.f to read, the new program extracts data and stores
c the data in arrays: 

c      common/uad/ug(nxgauss,nygauss,kl),vg(nxgauss,nygauss,kl),
c     .           tg(nxgauss,nygauss,kl),qgg(nxgauss,nygauss,kl)
c      common/sfcd/pslg(nxgauss,nygauss),albg(nxgauss,nygauss),
c     .           zsg(nxgauss,nygauss),tssg(nxgauss,nygauss),
c     .           ts1g(nxgauss,nygauss),ts2g(nxgauss,nygauss),
c     .           w1g(nxgauss,nygauss),w2g(nxgauss,nygauss),
c     .           precipg(nxgauss,nygauss)

c then these arrays are past to lintrp and sstamip for interpolation. 
c Results are past back to main via common block arrays:
c      common /arrays/ t(ifull,kl),u(ifull,kl),v(ifull,kl),psl(ifull),
c     .                qg(ifull,kl),ps(ifull),zs(ifull),tss(ifull),
c     .                almsk(ifull),precip(ifull),alb(ifull)

c The results are written to output file by subroutine outfile and used
c as input for DARLAM-style model.

      program c9tocc
c reads the netcdf version of the csiro9 history file and converts to
c lam/c-c boundary condition format; this one for formulation of csiro9 model
c in which soil moisture is saved once a day in a separate file (V4-6-5l)
      parameter(nxgauss=64,nygauss=56)
      parameter(ntimes=3)        ! number of times per day netCDF archive saved
      include 'newmpar.h'
      include 'dates.h'   ! to pass ds from setxyz
      include 'latlong.h'  ! rlat,rlong
      include 'parm.h'
      include 'xyzinfo.h'  ! x,y,z,wts
      common/gauss/ gausslat(nygauss)
      common/sigs/sig(kl),sigmh(kl),dsig(kl),rata(kl),ratb(kl),bet(kl)
      common/soil2/ts1(ifull),ts2(ifull),w1(ifull),w2(ifull)

c     arrays for input
      common/uad/ug(nxgauss,nygauss,kl),vg(nxgauss,nygauss,kl),
     .           tg(nxgauss,nygauss,kl),qgg(nxgauss,nygauss,kl)
      common/sfcd/pslg(nxgauss,nygauss),albg(nxgauss,nygauss),
     .           zsg(nxgauss,nygauss),tssg(nxgauss,nygauss),
     .           ts1g(nxgauss,nygauss),ts2g(nxgauss,nygauss),
     .           w1g(nxgauss,nygauss),w2g(nxgauss,nygauss),
     .           precipg(nxgauss,nygauss)

c     output variables from horizontal interpolation
      common /arrays/ t(ifull,kl),u(ifull,kl),v(ifull,kl),psl(ifull),
     .                qg(ifull,kl),ps(ifull),zs(ifull),tss(ifull),
     .                almsk(ifull),precip(ifull),alb(ifull)

      real sst2(nxgauss,nygauss)
      real zsdum(ifull)
      integer lmsk(nxgauss,nygauss), ilmsk(nxgauss)
      integer numday(12)                   ! number of days in each month
      character*40 halast,hblast,hafile,hbfile
      character name*30
      character*80 ofile,topofile, maskfile
      character cyymmlast*4, cyymm*4, cdate*2, ckdate*6
      character header*47

c     control file namelist
      namelist/control/yrlast,mmlast,halast,hblast,iyear,imonth
     .                ,hafile,hbfile,ofile,topofile,maskfile,meso

c     data statements: g - acceleration due to gravity
      data g/9.80665/,pi/3.1415926536/,n1x/0/
c     number of days in each month
      data numday/31,28,31,30,31,30,31,31,30,31,30,31/

c     the gaussian latitudes; hard coded here for wave 21 input
c     data is organised in order of SOUTH to NORTH pole.
      data gausslat/-87.5613,-84.4022,-81.2245,-78.0425,-74.8590,
     +  -71.6747,-68.4899,-65.3049,-62.1197,-58.9343,-55.7489,
     +  -52.5634,-49.3779,-46.1924,-43.0068,-39.8211,-36.6355,
     +  -33.4498,-30.2642,-27.0785,-23.8928,-20.7071,-17.5214,
     +  -14.3357,-11.1500,-7.9643,-4.7786,-1.5929,1.5929,4.7786,
     +  7.9643,11.1500,14.3357,17.5214,20.7071,23.8928,27.0785,
     +  30.2642,33.4498,36.6355,39.8211,43.0068,46.1924,49.3779,
     +  52.5634,55.7489,58.9343,62.1197,65.3049,68.4899,71.6747,
     +  74.8590,78.0425,81.2245,84.4022,87.5613/

c     open control namelist file
      open(85,file='control.nml',status='old')
      read (85, control)
      write(6, control)
      id=11
      jd=148
      idjd=id+il*(jd-1)
c     N.B. by setting n1x=1, only write out for ndayx and ntimex
      n1x=1
      ndayx=1
      ntimex=3

      open(unit=2,file=topofile,status='old')
      read(2,'(i3,i4,2f6.1,f5.2,f9.0,a47)')
     .          ilx,jlx,rlong0,rlat0,schmidt,dsx,header
      print *,ilx,jlx,rlong0,rlat0,schmidt,dsx,header
      if(ilx.ne.il.or.jlx.ne.jl)
     .            stop 'wrong topo file supplied'

c     set up cc geometry
      call setxyz
      do iq=1,ifull
c       convert conformal cubic lats & longs to degrees (-90 to 90) & (0 to 360)
        rlat(iq)=rlat(iq)*180./pi
        rlong(iq)=rlong(iq)*180./pi
        if(rlong(iq).lt.0)rlong(iq)=rlong(iq)+360.
        if(iq.eq.idjd)print *,'before conversion id,jd,rlong,rlat',
     .                                  id,jd,rlong(idjd),rlat(idjd)
c       Note that rlat and rlong are only used in calls to scinex, so:
c       convert longitudes of conformal-cubic grid to a coordinate system
c       based on gaussian longs; ailong = 1 is greenwich meridian
        rlong(iq) = (rlong(iq)*nxgauss/360.) + 1.
c       same for c-c latitudes to gaussian lats.
        rlat(iq)=max(gausslat(1),min(rlat(iq),gausslat(nygauss)))   ! jlm
        do ik=1,nygauss
          if(gausslat(ik+1).ge.rlat(iq)) then
            rlat(iq) = (rlat(iq)-gausslat(ik))/(gausslat(ik+1) -
     .                 gausslat(ik)) + ik
            go to 12
          endif
        enddo  ! ik loop
12      continue
      enddo  ! iq loop
      print *,'after conversion rlong,rlat',rlong(idjd),rlat(idjd)

c     read in orog, land-sea mask
      read(2,*) zs      ! now used
      read(2,*) almsk
      print*,'input almsk',almsk(82,88)

c     global model land-sea-ice mask
      open(unit=3,file=maskfile,status='old')

c     read in specification of global model land-sea mask
c     this file is given by Martix Dix. It is stored in (64i1) format
c     for number of latitudes
      do i=1,nxgauss
         ilmsk(i)=0
      enddo
      do j=1,nygauss
         read(3,'(64i1)') ilmsk
         do i=1,nxgauss
            lmsk(i,j) = ilmsk(i)
         enddo
      enddo

c     because file 'landmask' is in order of NORTH to SOUTH
c     pole which reverses to that of 'gausslat' above therefore
c     needs to invert NORTH <=> SOUTH
      do 21 mg=1,nxgauss
        do 22 lg=1,nygauss
          sst2(mg,lg)=float( lmsk(mg,nygauss-lg+1) )
22      continue
21    continue
c     write(90,*)'sst2 in MAIN'
c     write(90,'(10f7.1)')((sst2(i,j),i=1,nxgauss),j=1,nygauss)
c
c     file 17 for sea surface temperature output by sstamip.f,
c     file 20 for others output by lintrp15.f

      if(meso.eq.3) then
        open(unit=20,file=ofile,form='unformatted')
      else
        stop 'open(unit=20,file=ofile)'
      endif

c     first, the files required from the previous month for starting the
c     model at the start of day 1 of the current month, since this record
c     is not stored in the current month's archive

c     this file has the soil moisture for initialization of same in LAM
      do j=1, nygauss
         do i=1,nxgauss
            albg(i,j)    = 0.
            precipg(i,j) = 0.
         enddo
      enddo
      ncid2 = ncopn(halast,ncnowrit,ierr)
      print *,' opening ',halast
      if(ierr.ne.0) then
        print *,' cannot open netCDF halast; error code ',ierr
        stop
      end if

      ncid = ncopn(hblast,ncnowrit,ierr)
      if(ierr.ne.0) then
        print *,' cannot open netCDF hblast; error code ',ierr
        stop
      end if

c     now the current month's files
      ncid4 = ncopn(hafile,ncnowrit,ierr)
      print *,' opening ',hafile
      if(ierr.ne.0) then
        print *,' cannot open netCDF hafile; error code ',ierr
        stop
      endif

      ncid3 = ncopn(hbfile,ncnowrit,ierr)
      if(ierr.ne.0) then
        print *,' cannot open netCDF hbfile; error code ',ierr
        stop
      endif
      
c     time id
      itimeid = ncvid(ncid,'time',ierr)

c     get sigma levels
      call ncagt(ncid2,sigid,'sigma_lev',sig,ierr)
      print *,' sig ',sig

c     read in initial value of deeper layers of soil temps from monthly means
c     of January from GCM runs
      open(unit=82,file='stb2.txt',status='old')
      open(unit=83,file='stb3.txt',status='old')
      read(82,*)((ts1g(i,j),i=1,nxgauss),j=1,nygauss)
      read(83,*)((ts2g(i,j),i=1,nxgauss),j=1,nygauss)

c     iyear-5 because the archive of this variable starts at year 6
c     print *,' iyear ',iyear,' imonth ',imonth
c     print *,' now calling histrd2 '
c     iyearm5 = iyear - 5
c     call histrd2(ncidtb2,iyear,imonth,'tb2',nxgauss,nygauss,ts1g)
c     call myxn(ts1g,nxgauss,nygauss,xm,xx,xn,'ncdf ts1')
c     ncidtb3 = ncopn('stb3_qm1.nc',ncnowrit,ierr)
c     if(ierr.ne.0) then
c       print *,' cannot open netCDF file 2; error code ',ierr
c       stop
c     end if
c     call histrd2(ncidtb3,iyear,imonth,'tb3',nxgauss,nygauss,ts2g)
c     call myxn(ts2g,nxgauss,nygauss,xm,xx,xn,'ncdf ts2')
c
c     here read in the last archive interval of the previous month
c
c     obtain number of days in month from month as read in command line
c
      mnth = imonth-1
      if(mnth.eq.0)mnth=12
      ndays = numday(mnth)

c     read from soil moisture file (saved once a day; want final value for
c     previous month)
      call histrd1(ncid2,ndays,'wfg',nxgauss,nygauss,w1g)
      call myxn(w1g,nxgauss,nygauss,xm,xx,xn,'ncdf w1')
      call histrd1(ncid2,ndays,'wfb',nxgauss,nygauss,w2g)
      call myxn(w2g,nxgauss,nygauss,xm,xx,xn,'ncdf w2')

c     extract precip
      call histrd1(ncid2,ndays,'rnd',nxgauss,nygauss,precipg)

c     ADJUST SOIL MOISTURE VALUES OVER OCEAN from 150 to .16

c     from other file for previous month, read initialization field which
c     is stored at iarch = nday*ntimes
      iarch = ndays*ntimes

c     winds: u
      do k=1,kl
         write(name,1000) 'u', k
 1000    format(a,i2.2)
         call histrd1(ncid,iarch,name,nxgauss,nygauss,ug(1,1,k))
         call myxn(ug(1,1,k),nxgauss,nygauss,xm,xx,xn,name)
      end do

c     winds: v
      do k=1,kl
         write(name,1000) 'v', k
         call histrd1(ncid,iarch,name,nxgauss,nygauss,vg(1,1,k))
         call myxn(vg(1,1,k),nxgauss,nygauss,xm,xx,xn,name)
      end do

c     temperature
      do k=1,kl
         write(name,1000) 't', k
         call histrd1(ncid,iarch,name,nxgauss,nygauss,tg(1,1,k))
         call myxn(tg(1,1,k),nxgauss,nygauss,xm,xx,xn,name)
      end do

c     specific humidity
      do k=1,kl
         write(name,1000) 'q', k
         call histrd1(ncid,iarch,name,nxgauss,nygauss,qgg(1,1,k))
         call myxn(qgg(1,1,k),nxgauss,nygauss,xm,xx,xn,name)
c        convert to g/kg for consistency with LAM input
         do j=1,nygauss
           do i=1,nxgauss
             qgg(i,j,k) = qgg(i,j,k)*1000.
           end do
         end do
         call myxn(qgg(1,1,k),nxgauss,nygauss,xm,xx,xn,name)
      end do

c     surface pressure
      call histrd1(ncid,iarch,'psf',nxgauss,nygauss,pslg)
      call myxn(pslg,nxgauss,nygauss,xm,xx,xn,'psf')
c     convert to units required by nested/cc model
      do j=1,nygauss
         do i=1,nxgauss
           pslg(i,j) = log(100.*pslg(i,j))
         end do
       end do
       call myxn(pslg,nxgauss,nygauss,xm,xx,xn,'psl')

c      surface height in metres
       call histrd1(ncid2,iarch,'zht',nxgauss,nygauss,zsg)
       call myxn(zsg,nxgauss,nygauss,xm,xx,xn,'zht')

c      convert to geopotential metres for input to nested model
       do j=1,nygauss
         do i=1,nxgauss
           zsg(i,j) = g*zsg(i,j)
         end do
      end do
      call myxn(zsg,nxgauss,nygauss,xm,xx,xn,'zs')

c     other single-level files
      call histrd1(ncid,iarch,'tsu',nxgauss,nygauss,tssg)
      call myxn(tssg,nxgauss,nygauss,xm,xx,xn,'tss')

c     time stamp
      call ncvgt1(ncid,itimeid,iarch,fdays,ierr)

c-----modified on 7/2/1995 by kcn and jlm
c-----to replace mins above for subsequence computing KDATE and KTIME
c-----KDATE has a form of yymmdd and KTIME has values 8, 16 and 0
c-----which will be used in lintrp15.f.
      nhour = nint(fdays * 24.)
      mday  = nhour/24
      ktime = (nhour - mday*24)*100       ! ktime=0,800,1600
      cyymm  = hafile(6:9)
      cyymmlast = cyymm
      mday   = 1

c     convert from integer to character
      if(mday.lt.10) then
           write(cdate,'(2i1)') 0,mday
      else
           write(cdate,'(i2)') mday
      endif
c     convert from character to integer
      ckdate=cyymm//cdate
      read(ckdate,'(i6)')kdate

      print *,'Read file ',hblast,' at date ',fdays
      print*,'mday, kdate, ktime',mday, kdate, ktime

c     horizontal interpolation
      call lintrp
      call sstamip (sst2,almsk,tssg,tss,0)
      call sstamip (sst2,almsk,ts1g,ts1,0)
      call sstamip (sst2,almsk,ts2g,ts2,0)
      call sstamip (sst2,almsk,w1g,w1,1)
      call sstamip (sst2,almsk,w2g,w2,2)

      do j=1,jl
         do i=1,il
            w1(i,j) = min(0.32,max(w1(i,j),0.))
            w2(i,j) = min(0.32,max(w2(i,j),0.))
         enddo
      enddo

c     write to output file
      if(n1x.eq.0)then
        print *,'calling outfile for end of previous month'
        call outfile
      endif

c     now process current month's files
c     loop over all required archive intervals

c***  temporary hard-coding here
      print *,' now reading current months files '

c---  loop over nday = number of days of the month
      nday = numday(imonth)
      if(n1x.eq.1)nday=ndayx

      do 950 iday=1,nday

c        read from soil moisture file as this is once a day
         call histrd1(ncid4,iday,'wfg',nxgauss,nygauss,w1g)
         call myxn(w1g,nxgauss,nygauss,xm,xx,xn,'ncdf w1')
         call histrd1(ncid4,iday,'wfb',nxgauss,nygauss,w2g)
         call myxn(w2g,nxgauss,nygauss,xm,xx,xn,'ncdf w2')

c        extract precip
         call histrd1(ncid4,iday,'rnd',nxgauss,nygauss,precipg)

c        ADJUST SOIL MOISTURE VALUES OVER OCEAN from 150 to .16
      
c------- loop over ntimes (=3 time intervals i.e 0, 8 and 16)

         do 900 itimes=1,ntimes
            iarch = (iday-1)*ntimes + itimes

c           winds: u
            do k=1,kl
               write(name,1000) 'u', k
               call histrd1(ncid3,iarch,name,nxgauss,nygauss,ug(1,1,k))
               call myxn(ug(1,1,k),nxgauss,nygauss,xm,xx,xn,name)
            end do            

c           winds: v
            do k=1,kl
               write(name,1000) 'v', k
               call histrd1(ncid3,iarch,name,nxgauss,nygauss,vg(1,1,k))
               call myxn(vg(1,1,k),nxgauss,nygauss,xm,xx,xn,name)
            end do

c           temperature
            do k=1,kl
               write(name,1000) 't', k
               call histrd1(ncid3,iarch,name,nxgauss,nygauss,tg(1,1,k))
               call myxn(tg(1,1,k),nxgauss,nygauss,xm,xx,xn,name)
            end do

c           specific humidity
            do k=1,kl
               write(name,1000) 'q', k
               call histrd1(ncid3,iarch,name,nxgauss,nygauss,qgg(1,1,k))
               call myxn(qgg(1,1,k),nxgauss,nygauss,xm,xx,xn,name)

c              convert to g/kg for consistency with LAM input
               do j=1,nygauss
                  do i=1,nxgauss
                     qgg(i,j,k) = qgg(i,j,k)*1000.
                  end do
               end do

             end do

c            surface pressure
             call histrd1(ncid3,iarch,'psf',nxgauss,nygauss,pslg)
             call myxn(pslg,nxgauss,nygauss,xm,xx,xn,'psf')
c            convert to units required by nested model
             do j=1,nygauss
                do i=1,nxgauss
                   pslg(i,j) = log(100.*pslg(i,j))
                end do
             end do
             call myxn(pslg,nxgauss,nygauss,xm,xx,xn,'psl')

c            other single-level files
             call histrd1(ncid3,iarch,'tsu',nxgauss,nygauss,tssg)
c            time stamp
             call ncvgt1(ncid3,itimeid,iarch,fdays,ierr)

c----        modified on 7/2/1995 by kcn and jlm. Application as above.
c----        NOTE that fdays = fdays+1, because from this procedure the first
c----        day of the month is 0.

             f1days = fdays+1.     ! add 1 to fdays so that first day of the
c----                              ! month is day ONE instead of ZERO.
             nhour = nint(f1days * 24.)
             mday  = nhour/24
             ktime = (nhour - mday*24)*100   !ktime=0,800,1600
             if (ktime.eq.0) then
                cyymm     = cyymmlast
                cyymmlast = hbfile(6:9)
             endif

             if(mday.lt.10) then
c                 convert from integer to character
                  write(cdate,'(2i1)') 0,mday
             else if(mday.gt.numday(imonth)) then
c                    convert cyymm from character to integer  idate
                     read(cyymm,'(i4)')idate
                     idate = idate+1           !add 1 to current yymm
c                    convert back integer to character
                     write(cyymm,'(i4)') idate
                     cyymmlast = cyymm
                     imonth = imonth+1         !next month
                     mday   = 1                !start from 1 for new month
                     write(cdate,'(2i1)') 0,mday
                  else 
                     write(cdate,'(i2)') mday
             endif
c            convert from character to integer
             ckdate=cyymm//cdate
             read(ckdate,'(i6)')kdate

             print *,'Doing processing for iday, itimes: ',iday,itimes
             print *,'Read file ',hbfile,' at date ',fdays
             print*,'mday, kdate, ktime',mday, kdate, ktime

c            horizontal interpolation
             call lintrp
             call sstamip (sst2,almsk,tssg,tss,0)
             call sstamip (sst2,almsk,ts1g,ts1,0)
             call sstamip (sst2,almsk,ts2g,ts2,0)
             call sstamip (sst2,almsk,w1g,w1,3)
             call sstamip (sst2,almsk,w2g,w2,4)

             do j=1,jl
                do i=1,il
                  w1(i,j) = min(0.32,max(w1(i,j),0.))
                  w2(i,j) = min(0.32,max(w2(i,j),0.))
                enddo
             enddo

c            write to output file     
c            N.B. by setting n1x=1, only write out for ndayx and ntimex
             if(n1x.eq.0.or.(iday.eq.ndayx.and.itimes.eq.ntimex))then
               print *,'calling outfile for iday,itimes: ',iday,itimes
               call outfile
             endif

900       continue       ! itimes=1,ntimes
950   continue           ! iday=1,nday
      end

      subroutine histrd1(histid,iarch,name,nxgauss,nygauss,var)
      character name*(*)
      integer start(3),count(3)
      integer*2 ivar(nxgauss,nygauss)
      real var(nxgauss,nygauss)

      start(1) = 1
      start(2) = 1
      start(3) = iarch
      count(1) = nxgauss
      count(2) = nygauss
      count(3) = 1

c     read data
      id = ncvid(histid,name,ierr)
      call ncvgt(histid,id,start,count,ivar,ierr)

c     obtain scaling factors and offsets from attributes
      call ncagt(histid,id,'add_offset',addoff,ierr)
      call ncagt(histid,id,'scale_factor',sf,ierr)

c     unpack data
      do j=1,nygauss
        do i=1,nxgauss
          var(i,j) = ivar(i,j)*sf + addoff
        end do
      end do
      return
      end

      subroutine histrd2(histid,iyear,imonth,name,nxgauss,nygauss,var)
      character name*(*)
      integer start(4),count(4)
      integer*2 ivar(nxgauss,nygauss)
      real var(nxgauss,nygauss)

      print *,' nygauss etc ',nygauss,nxgauss,iyear,imonth
      start(1) = 1
      start(2) = 1
      start(3) = imonth
      start(4) = iyear
      count(1) = nxgauss
      count(2) = nygauss
      count(3) = 1
      count(4) = 1

c     read data
      id = ncvid(histid,name,ierr)
      print*,'id= ',id,' start= ',start, ' count= ',count
      call ncvgt(histid,id,start,count,ivar,ierr)
      print *,' ncvgt '

c     obtain scaling factors and offsets from attributes
      call ncagt(histid,id,'add_offset',addoff,ierr)
      call ncagt(histid,id,'scale_factor',sf,ierr)
      print *,' ncagt '

c     unpack data
      do j=1,nygauss
        do i=1,nxgauss
          var(i,j) = ivar(i,j)*sf + addoff
        end do
      end do
      return
      end

      subroutine lintrp
c------------------------------------------------------------------------
c Version for CSIRO 9 LEVEL GCM to conformal-cubic model  jlm 29/7/96
c
c lintrp for DARLAM: Kevin Walsh, June 1988.
c Modified from the program sptoar9 written by John McGregor
c CSIRO Division of Atmospheric Research, Aspendale, Vic.
c
c linterp7 September 1990; Modified by Leon Rotstayn to read
c in output from specgau.f which converts output from Hal Gordon's
c model to Gaussian grid.
c linterp8 September 1991 : modified to handle 2 soil temp layers (LDR)
c Surface topography in metres and soil moistures from 0. to 0.36.
c Reads in topography from the normal archive file
c and calculates map factors and coriolis parameters
c
c This version reads soil moisture in from the archive like any other
c variable and interpolates it to the Conformal Grid.
c
c Modified May 1994: same as lintrp12.f but writes out projection
c parameters
c
c Substantially modified June, 1994: new, more efficient version
c Has more accurate interpolation of gaussian latitudes
c
c
c parameters: il -- number of points in the lambert conformal
c                       grid in the x direction
c             jl -- number of points in the lambert conformal
c                       grid in the y direction
c             kl -- the number of atmospheric levels in the input
c                     file, usually the number of levels in the
c                     spectral model
c             nxgauss -- number of points in the gaussian grid
c                        in the x direction
c             nygauss -- number of points in the gaussian grid
c                        in the y direction
c
c------------------------------------------------------------------------------
c arrays; units are mks except where noted
c     t ...  temperature
c     u ...  zonal wind speed
c     v ...  meridional wind speed
c     psl ... log of surface pressure
c     qg  ... specific humidity
c     ps  ... mean sea level pressure
c     zs  ... topographic height; units are meters in this version!!!
c     tss ... sea surface temperature
c     almsk ..land-sea mask calculated from sst as in spectral model;
c             in the latest version of the nested model, this is
c             calculated in the nested model from the topography
c     precip..precipitation in mm per day.
c     ts and ts2 are upper and lower soil temps respectively (9/91)
c
c------------------------------------------------------------------------------
      parameter(nxgauss=64,nygauss=56)
      include 'newmpar.h'
      include 'dates.h'
      include 'latlong.h'  ! rlat,rlong
      include 'map.h'
      include 'parm.h'
      include 'vecsuv.h'   ! vecsuv info
      include 'xyzinfo.h'  ! x,y,z,wts
c     soil array
      common/soil2/ts1(ifull),ts2(ifull),w1(ifull),w2(ifull)

c     variables on gaussian grid passed from main program 
      common/uad/ug(nxgauss,nygauss,kl),vg(nxgauss,nygauss,kl),
     .           tg(nxgauss,nygauss,kl),qgg(nxgauss,nygauss,kl)
      common/sfcd/pslg(nxgauss,nygauss),albg(nxgauss,nygauss),
     .           zsg(nxgauss,nygauss),tssg(nxgauss,nygauss),
     .           ts1g(nxgauss,nygauss),ts2g(nxgauss,nygauss),
     .           w1g(nxgauss,nygauss),w2g(nxgauss,nygauss),
     .           precipg(nxgauss,nygauss)

c     output variables from horizontal interpolation
      common /arrays/ t(ifull,kl),u(ifull,kl),v(ifull,kl),psl(ifull),
     .                qg(ifull,kl),ps(ifull),zs(ifull),tss(ifull),
     .                almsk(ifull),precip(ifull),alb(ifull)
      data pi/3.1415926536/

c     here use unstaggered lats and lons for u and v
c     For calculating zonal and meridional wind components, use the
c     following information, where theta is the angle between the
c     (ax,ay,az) vector [along the xg axis] and the zonal-component-vector:
c     veczon = k x r, i.e. (-y,x,0)/sqrt(x**2 + y**2)
c     vecmer = r x veczon, i.e. (-xz,-yz,x**2 + y**2)/sqrt(x**2 + y**2)
c     costh is (veczon . a) = (-y*ax + x*ay)/sqrt(x**2 + y**2)
c     sinth is (vecmer . a) = [-xz*ax - yz*ay + (x**2 + y**2)*az]/sqrt
c      using (r . a)=0, sinth collapses to az/sqrt(x**2 + y**2)
c     For rotated coordinated version, see JMcG's notes
      coslong=cos(rlong0*pi/180.)
      sinlong=sin(rlong0*pi/180.)
      coslat=cos(rlat0*pi/180.)
      sinlat=sin(rlat0*pi/180.)
      polenx=-coslat
      poleny=0.
      polenz=sinlat
      print *,'polenx,poleny,polenz ',polenx,poleny,polenz

      do iq=1,ifull
c       den=sqrt( max(x(iq)**2 + y(iq)**2,1.e-7) )  ! allow for poles
c       costh=(-y(iq)*ax(iq) + x(iq)*ay(iq))/den
c       sinth=az(iq)/den
c       set up unit zonal vector components
        zonx=            -polenz*y(iq)
        zony=polenz*x(iq)-polenx*z(iq)
        zonz=polenx*y(iq)
        den=sqrt( max(zonx**2 + zony**2 + zonz**2,1.e-7) )  ! allow for poles
        costh= (zonx*ax(iq)+zony*ay(iq)+zonz*az(iq))/den
        sinth=-(zonx*bx(iq)+zony*by(iq)+zonz*bz(iq))/den
c       if(iq.eq.ifull/2)print *,'x,y,sqrt ',
c    .                            x(iq),y(iq),sqrt(x(iq)**2+y(iq)**2)
c       if(iq.eq.ifull/2)print *,'zonx,zony,zonz ',zonx,zony,zonz
c       if(iq.eq.ifull/2)print *,'iq,costh,sinth ',iq,costh,sinth
          do k=1,kl
             call scinex(ug(1,1,k),nxgauss,nygauss,1,1,nxgauss,
     .                   nygauss,rlong(iq),rlat(iq),uzon)
             call scinex(vg(1,1,k),nxgauss,nygauss,1,1,nxgauss,
     .                   nygauss,rlong(iq),rlat(iq),vmer)
c           calculate u and v relative to the cc grid,
c           using components of gaussian grid u and v with theta
            u(iq,k)= costh*uzon+sinth*vmer
            v(iq,k)=-sinth*uzon+costh*vmer
          enddo  ! k loop
      enddo      ! iq loop

!c     now stagger the winds (for writing-out purposes only)
!      do k=1,kl
!        do iq=1,ifull
!          w1(iq)=u(iq,k)   ! w1 just as temporary storage
!          w2(iq)=v(iq,k)   ! w2 just as temporary storage
!        enddo      ! iq loop
!        call staguv(w1,w2,u(1,1,k),v(1,1,k))
!      enddo  ! k loop

c     unstaggered variables

c     t, qg, psl, zs, ts, precip, ts2, ts1, w1, w2, alb

c     3D variables: t and qg

      do k=1,kl
        do iq=1,ifull
            call scinex(tg(1,1,k),nxgauss,nygauss,1,1,nxgauss,nygauss,
     &                  rlong(iq),rlat(iq),t(iq,k))
            call scinex(qgg(1,1,k),nxgauss,nygauss,1,1,nxgauss,nygauss,
     &                  rlong(iq),rlat(iq),qg(iq,k))
          enddo      ! iq loop
      enddo  ! k loop

c
c     2D variables: psl,zs,ts,precip,ts1,ts2,w1,w2,alb
      do iq=1,ifull
c           psl
            call scinex(pslg,nxgauss,nygauss,1,1,nxgauss,nygauss,
     &                  rlong(iq),rlat(iq),psl(iq))
            psl(iq)=psl(iq) - 11.5129
c           zs
            call scinex(zsg,nxgauss,nygauss,1,1,nxgauss,nygauss,
     &                  rlong(iq),rlat(iq),zs(iq))
c           precip
            call scinex(precipg,nxgauss,nygauss,1,1,nxgauss,nygauss,
     &                  rlong(iq),rlat(iq),bdum)
            precip(iq) = max(0.,bdum)
c           albedo
            call scinex(albg,nxgauss,nygauss,1,1,nxgauss,nygauss,
     &                  rlong(iq),rlat(iq),bdum)
            alb(iq) = max(0.,bdum)
      enddo      ! iq loop
      return
      end

      subroutine scinex( scala , msize , nsize , mmin , nmin
     *, mmax , nmax , gm , gn , scinto )
c-----------------------------------------------------------------------
c
c this subroutine produces the value scinto of a scalar field at a point
c gm,gn by interpolation or extrapolation of the field scala
c (2-direction bessel interpolation formula). mmin,mmax and nmin,nmax are
c the boundaries of the grid array.
c
c        this version uses lagrangian cubic formula
c
c     output: scinto
c
c-----------------------------------------------------------------------

      dimension scala(msize,nsize)
      scinto=0.
      igm=gm
      jgn=gn
      fm=gm-igm
      fn=gn-jgn
      if(fm.lt.1.e-06)fm=0.
      if(fn.lt.1.e-06)fn=0.
      ms=mmax-1
      ns=nmax-1
      mr=mmin+1
      nr=nmin+1
      if(gm.lt.mmax)go to 60
      if(gn.lt.nmax)go to 20
      e=gm-mmax
      t1=e*(scala(mmax,nmax)-scala(ms,nmax))
      e=gn-nmax
      t2=e*(scala(mmax,nmax)-scala(mmax,ns))
      scinto=scala(mmax,nmax)+t1+t2
      return

   20 if(gn.ge.nmin)go to 40
      e=gm-mmax
      t1=e*(scala(mmax,nmin)-scala(ms,nmin))
      e=nmin-gn
      t2=e*(scala(mmax,nmin)-scala(mmax,nr))
      scinto=scala(mmax,nmin)+t1+t2
      return

   40 p=scala(mmax,jgn)+fn*(scala(mmax,jgn+1)-scala(mmax,jgn))
      h=scala(ms,jgn)+fn*(scala(ms,jgn+1)-scala(ms,jgn))
      e=gm-mmax
      scinto=p+e*(p-h)
      return

   60 if(gm.ge.mmin)go to 140
      if(gn.lt.nmax)go to 80
      e=gn-nmax
      t2=e*(scala(mmin,nmax)-scala(mmin,ns))
      e=mmin-gm
      t1=e*(scala(mmin,nmax)-scala(mr,nmax))
      scinto=scala(mmin,nmax)+t1+t2
      return

   80 if(gn.ge.nmin)go to 100
      e=nmin-gn
      t2=e*(scala(mmin,nmin)-scala(mmin,nr))
      e=mmin-gm
      t1=e*(scala(mmin,nmin)-scala(mr,nmin))
      scinto=scala(mmin,nmin)+t1+t2
      return

  100 e=mmin-gm
      p=scala(mmin,jgn)+fn*(scala(mmin,jgn+1)-scala(mmin,jgn))
      h=scala(mr,jgn)+fn*(scala(mr,jgn+1)-scala(mr,jgn))
      scinto=p+e*(p-h)
      return

  120 e=gn-nmax
      p=scala(igm,nmax)+fm*(scala(igm+1,nmax)-scala(igm,nmax))
      h=scala(igm,ns)+fm*(scala(igm+1,ns)-scala(igm,ns))
      scinto=p+e*(p-h)
      return

  140 if(gn.ge.nmax)go to 120
      if(gn.ge.nmin)go to 160
      e=nmin-gn
      p=scala(igm,nmin)+fm*(scala(igm+1,nmin)-scala(igm,nmin))
      h=scala(igm,nr)+fm*(scala(igm+1,nr)-scala(igm,nr))
      scinto=p+e*(p-h)
      return

  160 if(gm.lt.ms.and.gm.ge.mr.and.gn.lt.ns.and.gn.ge.nr)go to 180
      p=scala(igm+1,jgn)+fn*(scala(igm+1,jgn+1)-scala(igm+1,jgn))
         h=scala(igm,jgn)+fn*(scala(igm,jgn+1)-scala(igm,jgn))
      scinto=h+fm*(p-h)
      return

  180 continue
      s1 = fm+1.0
      s2 = fm
      s3 = fm-1.0
      s4 = fm-2.0
      s12 = s1*s2
      s34 = s3*s4
      a = -s2*s34
      b = 3.0*s1*s34
      c = -3.0*s12*s4
      d = s12*s3
      x1 = a*scala(igm-1,jgn-1) + b*scala(igm  ,jgn-1)
     x   + c*scala(igm+1,jgn-1) + d*scala(igm+2,jgn-1)
      x2 = a*scala(igm-1,jgn  ) + b*scala(igm  ,jgn  )
     x   + c*scala(igm+1,jgn  ) + d*scala(igm+2,jgn  )
      x3 = a*scala(igm-1,jgn+1) + b*scala(igm  ,jgn+1)
     x   + c*scala(igm+1,jgn+1) + d*scala(igm+2,jgn+1)
      x4 = a*scala(igm-1,jgn+2) + b*scala(igm  ,jgn+2)
     x   + c*scala(igm+1,jgn+2) + d*scala(igm+2,jgn+2)
      s1 = fn+1.0
      s2 = fn
      s3 = fn-1.0
      s4 = fn-2.0
      s12 = s1*s2
      s34 = s3*s4
      a = -s2*s34
      b = 3.0*s1*s34
      c = -3.0*s12*s4
      d = s12*s3
      y = a*x1 + b*x2 + c*x3 + d*x4
      scinto = y/36.0
      return
      end

      subroutine outfile
c----------------------------------------------------------------------
c  this subroutine writes out the interpolated file in a format
c  consistent with the input to the nested model
c----------------------------------------------------------------------- 
      include 'newmpar.h'
      include 'dates.h'
      include 'map.h'
      include 'parm.h'
      common/arrays/ t(ifull,kl),u(ifull,kl),v(ifull,kl)
     1 ,psl(ifull),qg(ifull,kl),ps(ifull),zs(ifull),tss(ifull),
     2  almsk(ifull),precip(ifull),alb(ifull)
      common/sigs/q(kl),qmh(kl),dq(kl),rata(kl),ratb(kl),bet(kl)
      common/soil2/ts(ifull),ts2(ifull),w(ifull),w2(ifull)
      real psa(8001),psm(8001)
      character rundate*10
c     csiro9 sigma levels
      data q/.9793,.9136,.8032,.6598,.5,.3401,.1968,.0864,.0208/

      ndt=dt
c     some hard-coding for "analysis' file
      m=3
      nsd = -2
      nbd = 3
      nps = 2
      mex = 4
      mup = 3
      nhor = 0
      nkuo = 4
      nsi = 0
      nmi = 0
      ndt = 0
      npsav = 0
      rundate = 'start     '
      khdif = 0
      kwt = kl
      ik=il
      jk=jl
      kk=kl
      nvad = 0
      nqg = 3
      nrun = 0
      nrunx = 0
      lbd = 0
      nem = 2
      ndum = 0
      timeb=0
      timelb=0.
      ntsur=3    ! NTSUR > 2 means darlam expects 2 soil temp layers
      rhkuo=0.
      difknbd=0.

c     not used for anything
      du=0.
      tanl=0.
      nx3=0
      nx4=0


c     only unformatted output is now provided

      write(20) kdate,ktime,ktau,ik,jk,kk,m,nsd,meso,nbd,nps,
     2 mex,mup,nem,nsi,nmi,ndt,npsav,rundate,nhor,nkuo,khdif,
     3 kwt,nx3,nx4,timeb,timelb,ds,nvad,nqg,lbd,nrun,nrunx
     4 ,(ndum,nn=1,8),ntsur,(ndum,nn=1,8)
     6 ,difknbd,rhkuo,du,tanl,rlong0,rlat0,schmidt

      write(20) q
      write(20) psl

c     mean sea level pressure
      call mslp
      write(20) ps
      write(20) zs
      write(20) em
      write(20) f

c     flag for inclusion of land-sea mask in surface temps
      if(nem.eq.2) then
        do 1001 i=1,il
          do 1001 j=1,jl
            if(almsk(i,j).lt..5 .and. almsk(i,j).gt.-.5) then
c             a sea point; set tss to its negative
              tss(i,j) = -tss(i,j)
            else
              if(almsk(i,j).le.-.5) then
c               a sea-ice point
                tss(i,j) = tss(i,j) - 273.16
              end if
            end if
 1001   continue
      end if

      write(20) tss
      if(nqg.ge.2) then
        write(20) precip
      end if
      if(nqg.ge.3) then
        write(20) ts2
        write(20) ts
        write(20) w
        write(20) w2
      end if
      if(nqg.ge.4) write(20) alb
      write(20) (((t(i,j,k),i=1,il),j=1,jl),k=1,kwt)
      write(20) (((u(i,j,k),i=1,il),j=1,jl),k=1,kwt)
      write(20) (((v(i,j,k),i=1,il),j=1,jl),k=1,kwt)
      write(20) (((qg(i,j,k),i=1,il),j=1,jl),k=1,kwt)
      if(npsav.ne.0) then
        write(20) (psa(n),psm(n),n=1,npsav)
      end if

      return
      end
c*************************************************************************
      subroutine mslp
c     routine to calculate the mean sea level pressure
      include 'newmpar.h'
      common/arrays/ t(ifull,kl),u(ifull,kl),v(ifull,kl)
     1 ,psl(ifull),qg(ifull,kl),ps(ifull),zs(ifull),tss(ifull),
     2  almsk(ifull),precip(ifull),alb(ifull)
      common/sigs/sig(kl),sigmh(kl),dsig(kl),rata(kl),ratb(kl),bet(kl)
      data ktemp/2/

      if(kl.lt.9)ktemp=1
      c=9.806/6.5e-3
      conr=c/287.
      con=sig(ktemp)**(287./c)/c
      do 2 j=1,jl
      do 2 i=1,il
      ps(i,j) = exp(psl(i,j)+11.5129)
2     ps(i,j)=ps(i,j)*(1.+con*zs(i,j)/t(i,j,ktemp))**conr
      return
      end

      subroutine sstamip (sst2,almsk,varin,varout,iprint)
c------------------------------------------------------------------------
c interpolates output of amitolam
c
c this version creates an initial surface temperature file
c interpolating the gcm initial conditions only over the ocean
c
c parameters: il -- number of points of the conformal
c                       grid in the x direction
c             jl -- number of points in the conformal
c                       grid in the y direction
c             kl -- the number of atmospheric levels in the input
c                     file, usually the number of levels in the
c                     spectral model
c             nxgauss -- number of points in the gaussian grid
c                        in the x direction
c             nygauss -- number of points in the gaussian grid
c                        in the y direction
c
c     output: varout array
c-------------------------------------------------------------------------

      include 'newmpar.h'
      include 'latlong.h'  ! rlat,rlong
      include 'parm.h'
      include 'xyzinfo.h'  ! x,y,z,wts
      parameter(nxgauss=64,nygauss=56)
      common/gauss/ gausslat(nygauss)
c      common /arrays/ t(ifull,kl),u(ifull,kl),v(ifull,kl),psl(ifull),
c     .                qg(ifull,kl),ps(ifull),zs(ifull),tss(ifull),
c     .                almsk(ifull),precip(ifull),alb(ifull)
c      common/sfcd/pslg(nxgauss,nygauss),albg(nxgauss,nygauss),
c     .           zsg(nxgauss,nygauss),tssg(nxgauss,nygauss),
c     .           ts1g(nxgauss,nygauss),ts2g(nxgauss,nygauss),
c     .           w1g(nxgauss,nygauss),w2g(nxgauss,nygauss),
c     .           precipg(nxgauss,nygauss)

c     various work arrays
      real sst2(nxgauss,nygauss),almsk(ifull)
      real varout(ifull),varin(nxgauss,nygauss)
      real sstin(nxgauss,nygauss),landin(nxgauss,nygauss)

c special value for land/sea mask. This is used to mask whether the grid cell
c required to be interpolated or not.
      data spval/9999./

      idjd=id+il*(jd-1)
      longd=rlong(idjd)
      latd=rlat(idjd)

c       assign special value to gaussian grid land points
        do 31 j=1,nygauss
          do 31 i=1,nxgauss
            if(sst2(i,j).gt.3.5) then
              sstin(i,j) = spval
            else
              sstin(i,j) = varin(i,j)
            end if
 31     continue

       if(iprint.gt.0)then
         print *,'in sstamip for iprint,id,jd: ',iprint,id,jd
         print*,'BEFORE FILL LAND POINT'
         print *,'sst2  longd,latd',longd,latd
         do nn=3,-3,-1
          print *,(sst2(longd+ii,latd+nn),ii=-3,3)
         enddo
         print *,'sstin  longd,latd',longd,latd
         do nn=3,-3,-1
          print *,(sstin(longd+ii,latd+nn),ii=-3,3)
         enddo
       endif

c      fill land points with adjacent ocean points
       call fill(sstin,nxgauss,nygauss,spval)

       if(iprint.gt.0)then
         print*,'AFTER FILL LAND POINT'
         print *,'sstin  longd,latd',longd,latd
         do nn=3,-3,-1
          print *,(sstin(longd+ii,latd+nn),ii=-3,3)
         enddo
       endif

c       now assign special values to ocean points for array landin
        do 32 i=1,nxgauss
          do 32 j=1,nygauss
            if(sst2(i,j).lt.3.5) then
              landin(i,j) = spval
            else
              landin(i,j) = varin(i,j)
            end if
 32     continue

       if(iprint.gt.0)then
         print*,'BEFORE FILL OCEAN POINT'
         print *,'landin  longd,latd',longd,latd
         do nn=3,-3,-1
          print *,(landin(longd+ii,latd+nn),ii=-3,3)
         enddo
       endif

c       fill ocean points with adjacent land points
        call fill(landin,nxgauss,nygauss,spval)

       if(iprint.gt.0)then
         print*,'AFTER FILL OCEAN POINT'
         print *,'landin  longd,latd',longd,latd
         do nn=3,-3,-1
          print *,(landin(longd+ii,latd+nn),ii=-3,3)
         enddo
       endif

c       loops to convert gaussian grid points to cubic-conformal grid points
        do 7 iq=1,ifull
            call scinex(sstin,nxgauss,nygauss,1,1,nxgauss,nygauss,
     +           rlong(iq),rlat(iq),bout)
c           interpolate land values
            call scinex(landin,nxgauss,nygauss,1,1,nxgauss,nygauss,
     +           rlong(iq),rlat(iq),cout)

c now check the conformal grid's land/sea mask; if the point is land,
c then the interpolated land point from landin is accepted; if it is
c sea, then the interpolated point from ssstin is accepted
c
          if(almsk(iq).gt..5) then
            varout(iq) = cout
          else
            varout(iq) = bout
          end if
          if(iprint.gt.0)varout(iq)=max(varout(iq),.0001) ! jlm for  w,w2

       if(iprint.gt.0.and.iq.eq.idjd)then
         long=rlong(iq)
         lat=rlat(iq)
         print *,'in sstamip for iprint,id,jd: ',iprint,id,jd
         print *,'rlong,rlat,almsk: ',rlong(iq),rlat(iq),almsk(iq)
         print *,'bout,cout,varout: ',bout,cout,varout(iq)
         long=rlong(iq)
         lat=rlat(iq)
         print *,'sstin  long,lat',long,lat
         do nn=3,-3,-1
          print *,(sstin(long+ii,lat+nn),ii=-3,3)
         enddo
         print *,'landin  long,lat',long,lat
         do nn=3,-3,-1
          print *,(landin(long+ii,lat+nn),ii=-3,3)
         enddo
       endif
 7      continue

      return
      end

      subroutine fill(a,ifull,value)
      parameter(nxgauss=64,nygauss=56)
c     routine fills in interior of an array which has undefined points
      real a(ifull)         ! input and output array
      real value            ! array value denoting undefined
      real b(nxgauss,nygauss)
      dimension in(8), jn(8)   ! specifies neighbours
      data in/-1,-1,-1,0,1,1,1,0/
      data jn/-1,0,1,1,1,0,-1,-1/

2     nrem=0
      do 6 j=2,jl-1
      do 6 i=2,il-1
      b(i,j)=a(i,j)
      if(a(i,j).eq.value)then
        neighb=0
        av=0.
        do 4 nbs=1,8
        if(a(i+in(nbs),j+jn(nbs)).ne.value)then
          neighb=neighb+1
          av=av+a(i+in(nbs),j+jn(nbs))
        endif
4       continue
        if(neighb.gt.0)then
          b(i,j)=av/neighb
        else
          nrem=nrem+1    ! number of remaining points
        endif
      endif
6     continue
      do j=2,jl-1
       do i=2,il-1
        a(i,j)=b(i,j)
       enddo
      enddo
      if(nrem.gt.0)go to 2

c     fix up any boundary points
      do 7 i=2,il-1
      if(a(i,1).eq.value)a(i,1)=a(i,2)
      if(a(i,jl).eq.value)a(i,jl)=a(i,jl-1)
7     continue
      do 8 j=2,jl-1
      if(a(1,j).eq.value)a(1,j)=a(2,j)
      if(a(il,j).eq.value)a(il,j)=a(il-1,j)
8     continue
      if(a(1,1).eq.value)a(1,1)=a(2,2)
      if(a(il,1).eq.value)a(il,1)=a(il-1,2)
      if(a(1,jl).eq.value)a(1,jl)=a(2,jl-1)
      if(a(ifull).eq.value)a(ifull)=a(il-1,jl-1)
      return
      end

      subroutine myxn(dat,n1,n2,xm,xx,xn,label)
      real dat(n1,n2)
      character*(*) label
      xm=0.
      xx=-1.e29
      xn=+1.e29
      do j=1,n2
       do i=1,n1
         xm=xm+dat(i,j)
         xx=max(xx,dat(i,j))
         xn=min(xn,dat(i,j))
       enddo
      enddo
      xm=xm/float(n1*n2)
      print *,label,'  ',xm,xx,xn
      return
      end
      include 'setxyz.f'
      include 'jimcc.f'
!     include 'staguv.f'
