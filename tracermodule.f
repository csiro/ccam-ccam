      module tracermodule
c
      implicit none
      real, dimension(:,:,:), save, allocatable :: co2emhr,co2em123
      real, dimension(:,:), save, allocatable :: co2hr,co2em
      integer, dimension(:), save, allocatable :: nghr
      character(len=13), dimension(:), save, allocatable :: tracname
      character(len=13), dimension(:), save, allocatable :: tractype
      character(len=13), dimension(:), save, allocatable :: tracunit
      character(len=50), dimension(:), save, allocatable :: tracfile
      real, dimension(:), save, allocatable :: tracival,trmass
      integer, dimension(:), save, allocatable :: tracinterp,igashr
      real, dimension(:,:), save, allocatable :: tracdaytime
      integer, save :: numtracer, nhr
      integer, save :: unit_trout=131
     

      character(len=80), save :: tracerlist=''
      character(len=80), save :: sitefile=''
      character(len=80), save :: shipfile=''
      character(len=80), save :: trout='tracer.stdout'
      real,save :: tracvalin=-999  !default to initialise from restart

      namelist/trfiles/tracerlist,sitefile,shipfile,trout,tracvalin


      contains

c ***************************************************************************
      subroutine init_tracer
      implicit none
      character(len=80) :: header
      integer nt
      include 'newmpar.h'
      include 'tracers.h'


      open(unit=unit_trout,file=trout,form='formatted')
c     first read in a list of tracers, then determine what emission data is required
      open(unit=130,file=tracerlist,form='formatted')
      read(130,*) header
      read(130,*) numtracer
      if (numtracer.ne.ngas) 
     & stop 'wrong number of tracers in tracer.dat file'
      allocate(tracname(numtracer),tractype(numtracer))
      allocate(tracinterp(numtracer),tracunit(numtracer))
      allocate(tracfile(numtracer),igashr(numtracer))
      allocate(tracival(numtracer),trmass(numtracer))
      allocate(tracdaytime(numtracer,2))
      nhr = 0
      do nt=1,numtracer
        read(130,*) tracname(nt),tracival(nt),tractype(nt),tracfile(nt)
        if (tractype(nt).eq.'monrep'.or.tractype(nt).eq.'month') then
          tracinterp(nt)=1
        elseif (tractype(nt).eq.'hourly'.or.
     &          tractype(nt).eq.'3hourly'.or.
     &          tractype(nt).eq.'daily') then
          tracinterp(nt)=2
          nhr = nhr+1
          igashr(nt) = nhr
        elseif (tractype(nt).eq.'online') then
          tracunit(nt) = 'gC/m2/s'
          tracinterp(nt) = -1
        else
          tracinterp(nt)=0
        endif
        write(unit_trout,999) 'Tracer ',nt,tracname(nt),tractype(nt),
     &tracfile(nt),tracinterp(nt)
 999    format(a7,i5,a13,a13,a50,i3)
c
      enddo
      allocate(nghr(nhr))

c     initialise array for monthly average tracer
      traver = 0.

      return
      end subroutine

c ***********************************************************************
      subroutine tracini
c     initial value now read from tracerlist 
      implicit none
      include 'newmpar.h'
      include 'tracers.h'
      integer i

      do i=1,ngas
        tr(:,:,i)=tracival(i)
      enddo
      write(unit_trout,*) 'tracini: ',tracival
      write(unit_trout,*) 'tracini: ',tr(1,1,:)

      if( iradon.ne.0.) then
        tr(:,:,max(1,iradon))=0.
      end if

      return
      end subroutine

c *********************************************************************
      subroutine readtracerflux(kdate)
c     needs to happen further down the code than reading the tracer
c     list file
      implicit none
      include 'newmpar.h'
      include 'tracers.h'
      real ajunk(3)
      integer nt,jyear,jmonth,kdate
c
c     set up and read surface flux data
      allocate(co2em123(ilt*jlt,3,numtracer))
      allocate(co2emhr(ilt*jlt,31*24+2,nhr))
      allocate(co2hr(31*24+2,nhr))
      allocate(co2em(ilt*jlt,numtracer))

      jyear=kdate/10000
      jmonth=(kdate-jyear*10000)/100

      do nt=1,numtracer
        select case(tracinterp(nt))
          case (0:1)
c           monthly data
            call readrco2(nt,jyear,jmonth,3,co2em123(:,:,nt),ajunk)
          case (2)
c           daily, 3 hourly, hourly
            call readrco2(nt,jyear,jmonth,31*24+2,
     &                 co2emhr(:,:,igashr(nt)),co2hr(:,igashr(nt)))
        end select
      enddo

      return
      end subroutine

c *************************************************************************
      subroutine readrco2(igas,iyr,imon,nflux,fluxin,co2time)
c     rml 23/09/03 largely rewritten to use netcdf files
      implicit none
      include 'newmpar.h' !il,jl,kl
      include 'parm.h' !nperday
c     include 'tracers.h'
c     include 'trcom2.h'
      character*50 filename
c     rml 25/08/04 added fluxunit variable
      character*13 fluxtype,fluxname,fluxunit
      integer nflux,iyr,imon,igas,ntime,lc,regnum,kount
      integer nprev,nnext,ncur,n1,n2,ntot,n
      real timeinc
c     nflux =3 for month interp case - last month, this month, next month
c     nflux=31*24+2 for daily, hourly, 3 hourly case
      real fluxin(il*jl,nflux),co2time(nflux),hr
      integer ncidfl,timedim,ntime,yearid,monthid,fluxid,hourid
      integer nregdim,nregion,dayid,regnum,lc,ierr
      integer, dimension(:), allocatable :: fluxyr,fluxmon
      real, dimension(:), allocatable :: fluxhr
      integer start(3),count(3)
      include 'netcdf.inc'

      fluxtype=tractype(igas)
      fluxname=tracname(igas)
      filename=tracfile(igas)
c
      if (trim(fluxtype).eq.'pulseoff'.or.
     &    trim(fluxtype).eq.'daypulseoff') then
c       no surface fluxes to read
        fluxin = 0.
        tracunit(igas)=''
        return
      else
        write(unit_trout,*)'reading ',trim(fluxname), ' with type ', 
     &  trim(fluxtype),' for ',iyr,imon,' from ',filename
        ierr = nf_open(filename,0,ncidfl)
        if (ierr.ne.nf_noerr) then
          write(unit_trout,*) ierr,igas,filename
          write(unit_trout,*) nf_strerror(ierr)
          stop 'flux file not found'
        endif
        ierr=nf_inq_dimid(ncidfl,'time',timedim)
        if (ierr.ne.nf_noerr) stop 'time dimension error'
        ierr=nf_inq_dimlen(ncidfl,timedim,ntime)
        if (ierr.ne.nf_noerr) stop 'time dimension length error'
        ierr=nf_inq_varid(ncidfl,'year',yearid)
        if (ierr.ne.nf_noerr) stop 'year variable not found'
        allocate(fluxyr(ntime),fluxmon(ntime))
        ierr=nf_get_var_int(ncidfl,yearid,fluxyr)
        if (ierr.ne.nf_noerr) stop 'year variable not read'
        ierr=nf_inq_varid(ncidfl,'month',monthid)
        if (ierr.ne.nf_noerr) stop 'month variable not found'
        ierr=nf_get_var_int(ncidfl,monthid,fluxmon)
        if (ierr.ne.nf_noerr) stop 'month variable not read'
c       read hours variable for daily, hourly, 3 hourly data
        if (nflux.eq.(31*24+2)) then
          allocate(fluxhr(ntime))
          ierr=nf_inq_varid(ncidfl,'hour',hourid)
          if (ierr.ne.nf_noerr) stop 'hour variable not found'
          ierr=nf_get_var_real(ncidfl,hourid,fluxhr)
          if (ierr.ne.nf_noerr) stop 'hour variable not read'
        endif
c       check fluxname first otherwise default to 'flux'
        ierr=nf_inq_varid(ncidfl,fluxname,fluxid)
        if (ierr.ne.nf_noerr) then
          ierr=nf_inq_varid(ncidfl,'flux',fluxid)
          if (ierr.ne.nf_noerr) stop 'flux variable not found'
        endif
c       rml 25/08/04 read flux units attribute
        fluxunit='             '
        ierr = nf_get_att_text(ncidfl,fluxid,'units',fluxunit)
c rml 08/11/04 added radon units
        tracunit(igas)=fluxunit
        if (trim(fluxunit).ne.'gC/m2/s'.and.
     &      trim(fluxunit).ne.'Bq/m2/s'.and.
     &      trim(fluxunit).ne.'mol/m2/s') then
          write(unit_trout,*) 'Units for ',trim(fluxname),
     &                        ' are ',trim(fluxunit)
          write(unit_trout,*) 'Code not set up for units other than gC/m
     &2/s or mol/m2/s or Bq/m2/s'
          stop 'fix flux units'
        endif
        if (trim(fluxtype).eq.'daypulseon') then
c         need to read sunset/sunrise times
          ierr=nf_inq_dimid(ncidfl,'nregion',nregdim)
          if (ierr.ne.nf_noerr) stop 'nregion dimension error'
          ierr=nf_inq_dimlen(ncidfl,nregdim,nregion)
          if (ierr.ne.nf_noerr) stop 'nregion dimension length error'
c         if (nregion.gt.ngasmax) stop 'tracdaytime array too small'
          ierr=nf_inq_varid(ncidfl,'daylight',dayid)
          if (ierr.ne.nf_noerr) stop 'daylight variable not found'
        endif
c
c       find required records
        if (nflux.eq.3) then
c         monthly/annual cases
          nprev=0
          nnext=0
          ncur=0
          if (trim(fluxtype).eq.'constant'.or.
     &        trim(fluxtype).eq.'pulseon'.or.
     &        trim(fluxtype).eq.'daypulseon') then
c           check ntime
            if (ntime.ne.1) stop 'flux file wrong ntime'
            ncur = 1
          elseif (trim(fluxtype).eq.'annual') then
            do n=1,ntime
              if (fluxyr(n).eq.iyr) ncur=n
            enddo
          elseif (trim(fluxtype).eq.'monrep') then
            if (ntime.ne.12) stop 'flux file wrong ntime'
            do n=1,ntime
              if (fluxmon(n).eq.imon) then
                ncur=n
                nprev=n-1
                if (nprev.eq.0) nprev=12
                nnext=n+1
                if (nnext.eq.13) nnext=1
              endif
            enddo
          else
c           monthly case
            do n=1,ntime
              if (fluxyr(n).eq.iyr.and.fluxmon(n).eq.imon) then
                ncur=n
                nprev=n-1
c               keep flux constant at ends of data
                if (nprev.eq.0) nprev=n
                nnext=n+1
                if (nnext.eq.ntime+1) nnext=n
              endif
            enddo
          endif
c
          write(unit_trout,*) 'reading ',ncur,fluxyr(ncur),fluxmon(ncur)
          if (ncur.eq.0) stop 'current year/month not in flux file'
c    
          fluxin=0.
          start(1)=1
          count(1)=il*jl; count(2)=1
c         read preceeding month if needed
          if (nprev.ne.0) then
            start(2)=nprev
            ierr=nf_get_vara_real(ncidfl,fluxid,start,count,fluxin(:,1))
            if (ierr.ne.nf_noerr) stop 'error reading fluxin prev'
          endif
c         read current month/year
          start(2)=ncur
          ierr=nf_get_vara_real(ncidfl,fluxid,start,count,fluxin(:,2))
          if (ierr.ne.nf_noerr) stop 'error reading fluxin cur'
c         read next month
          if (nnext.ne.0) then
            start(2)=nnext
            ierr=nf_get_vara_real(ncidfl,fluxid,start,count,fluxin(:,3))
            if (ierr.ne.nf_noerr) stop 'error reading fluxin next'
          endif
        else
c         daily, hourly, 3 hourly case
c         find first time in month
c         rml 06/01/06 extend case to cope with annually repeating or
c         real year fluxes
          do n=1,ntime
            if ( ( (fluxyr(n).eq.0).or.(fluxyr(n).eq.iyr) ) .and. 
     &           (fluxmon(n).eq.imon) ) then
              n1=n
              exit
            endif
          enddo
c         find last time in month
          do n=ntime,1,-1
            if ( ( (fluxyr(n).eq.0).or.(fluxyr(n).eq.iyr) ) .and. 
     &           (fluxmon(n).eq.imon) ) then
              n2=n
              exit
            endif
          enddo
c         read fluxes
          start(1)=1; count(1)=il*jl
          ntot=n2-n1+1
          start(2)=n1; count(2)=ntot
          ierr=nf_get_vara_real(ncidfl,fluxid,start,count,
     &                                fluxin(:,2:ntot+1))
          if (ierr.ne.nf_noerr) stop 'error reading fluxin '
c         read in last time of prev month and first time of next month
          if ((n1.eq.1).and.(fluxyr(n1).eq.0)) then
c           loop
            nprev=ntime
          elseif ((n1.eq.1).and.(fluxyr(n1).ne.0)) then
c           keep constant
            nprev=n1
          else
            nprev=n1-1
          endif
          start(2)=nprev; count(2)=1
          ierr=nf_get_vara_real(ncidfl,fluxid,start,count,fluxin(:,1))
          if (ierr.ne.nf_noerr) stop 'error reading fluxin '
          if ((n2.eq.ntime).and.(fluxyr(n2).eq.0)) then
            nnext=1
          elseif((n2.eq.ntime).and.(fluxyr(n2).ne.0)) then
            nnext=ntime
          else
            nnext=n2+1
          endif
          start(2)=nnext; count(2)=1
          ierr=nf_get_vara_real(ncidfl,fluxid,start,count,
     &                                      fluxin(:,ntot+2))
          if (ierr.ne.nf_noerr) stop 'error reading fluxin '
c
c         need to make an array with the hour data in
          co2time(2:ntot+1)=fluxhr(n1:n2)
          timeinc = fluxhr(n1+1)-fluxhr(n1)
          co2time(1)=co2time(2)-timeinc
          co2time(ntot+2)=co2time(ntot+1)+timeinc
        endif
c
        if (trim(fluxtype).eq.'daypulseon') then
c         read sunrise/sunset times for this month, region from file
          lc=len_trim(fluxname)
          read(fluxname(lc-2:lc),'(i3)') regnum
          if (regnum.gt.nregion) stop 'region number > nregion'
          start(1)=regnum; start(2)=imon; start(3)=1
          count(1)=1; count(2)=1; count(3)=2
          ierr=nf_get_vara_real(ncidfl,dayid,start,count,
     & tracdaytime(igas,:))
          if (ierr.ne.nf_noerr) stop 'error reading daylight '
c         count number of timesteps that source emitting for
          kount=0
          do n=1,nperday
            hr = 24.*float(n)/float(nperday)
            if (tracdaytime(igas,1).lt.tracdaytime(igas,2) .and.
     &          tracdaytime(igas,1).le.hr .and.
     &          tracdaytime(igas,2).ge.hr) kount=kount+1 
            if (tracdaytime(igas,1).gt.tracdaytime(igas,2) .and.
     &          (tracdaytime(igas,1).le.hr .or.
     &          tracdaytime(igas,2).ge.hr)) kount=kount+1 
          enddo
c         scale flux to allow for emission over fraction of day
c         just set flux to zero if no daylight
          if (kount.ne.0) then
            fluxin = fluxin*float(nperday)/float(kount)
          else
            fluxin = 0.
          endif
        endif
        ierr = nf_close(ncidfl)
        return
      endif
      end subroutine

c *************************************************************************
      subroutine interp_tracerflux(kdate,hrs_dt)
c     interpolates tracer flux to current timestep if required
c     tracinterp 0 for no interpolation, 1 for monthly, 2 for daily/hourly
c     co2em123(:,1,:) contains prev month, co2em123(:,2,:) current month/year 
c     co2em123(:,3,:) next month
c eak 19/11/03 interpolation of rlai included  
c        set up monthly first as have to do rlai anyway
      implicit none
      include 'newmpar.h' !kl needed for parm.h
      include 'parm.h' !ktau,nperday
      integer mdays(0:13)
      data mdays/31,31,28,31,30,31,30,31,31,30,31,30,31,31/
      integer iyr,month,m1,m2,igas,igh,kdate
      real a1,a2,a3,ratlm,ratlm2,hrmodel,hrs_dt
      logical found
   
c     this could go in the case section but then we'd do it for
c     every monthly tracer.  This way we do it once but may not
c     need it.
      iyr=kdate/10000
      month=(kdate-10000*iyr)/100
      a1 = float(nperday*mdays(month-1))/2.
      a2 = float(nperday*mdays(month))/2.
      a3 = float(nperday*mdays(month+1))/2.
      if (ktau.lt.a2) then
c       first half of month
        ratlm = (a1+float(ktau))/(a1+a2)
        m1 = 1; m2 = 2
      else
c       second half of month
        ratlm = (float(ktau)-a2)/(a2+a3)
        m1 = 2; m2 = 3
      endif

      do igas=1,numtracer
        select case(tracinterp(igas))

         case(0) ; co2em(:,igas)=co2em123(:,2,igas)

         case(1)
             co2em(:,igas) = (1.0-ratlm)*co2em123(:,m1,igas) +
     .                          ratlm*co2em123(:,m2,igas)

         case(2)
             hrmodel = (ktau-1)*hrs_dt
             igh=igashr(igas)
             if (ktau.eq.1) nghr(igh)=1
             found=.false.
             write(unit_trout,*) 'igas: ',igas,'igashr: ',igh
             write(unit_trout,*) 'hrmodel: ',hrmodel
             do while (.not.found)
               write(unit_trout,*) nghr(igh),co2hr(nghr(igh),igh),
     &                    co2hr(nghr(igh)+1,igh)
               if ((hrmodel.gt.co2hr(nghr(igh),igh)).and.
     &           (hrmodel.le.co2hr(nghr(igh)+1,igh))) then
                 found = .true.
                 ratlm2 = (hrmodel-co2hr(nghr(igh),igh))/
     &               (co2hr(nghr(igh)+1,igh)-co2hr(nghr(igh),igh))
                 co2em(:,igas)= (1.-ratlm2)*co2emhr(:,nghr(igh),igh) 
     &                           +  ratlm2*co2emhr(:,nghr(igh)+1,igh)
               else
                 nghr(igh) = nghr(igh)+1
c                this error check only useful for hourly resolution
                 if (nghr(igh).eq.(31*24+2)) stop 'hr flux error'
               endif
             enddo
        end select
      enddo

      return
      end subroutine

c ***************************************************************************
      subroutine tracer_mass(ktau,ntau)
c     rml 16/10/03 check tracer mass - just write out for <= 6 tracers
      implicit none
      include 'newmpar.h'
      include 'const_phys.h' ! rearth,fc_molm,fair_molm
      include 'arrays.h'   ! ps
      include 'sigs.h'     ! dsig
      include 'xyzinfo.h'  ! wts
      include 'tracers.h'  ! tr
      integer i,k,ktau,ntau,igas

      trmass=0
      do i=1,ifull
        do k=1,kl
          trmass(:) = trmass(:) + tr(i,k,:)*ps(i)*dsig(k)*wts(i)
        enddo
      enddo
c     scaling assumes CO2 with output in GtC?
      if (ngas.gt.6) then
        write(unit_trout,*) 'Trmass: ',ktau,
     &  -1*trmass(1:6)*4.*3.14159*(rearth**2)*fC_MolM/
     &          (grav*1e18*fAIR_MolM)
      else
        write(unit_trout,*) 'Trmass: ',ktau,
     &  -1*trmass(:)*4.*3.14159*(rearth**2)*fC_MolM/
     &         (grav*1e18*fAIR_MolM)
      endif

c     also update tracer average array here
      do igas=1,ngas
        traver(:,:,igas)=traver(:,:,igas)+tr(1:ifull,1:kl,igas)/ntau
      enddo


      end subroutine

      end module
