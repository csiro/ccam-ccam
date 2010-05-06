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
c     rml 18/09/07 added tracmin, tracmax (tracmin to replace gasmin in
c     tracers.h, used in adjust5.f and outcdf.f, tracmax to be used
c     in outcdf.f (current setting of trmax unreliable)
      real, dimension(:), save, allocatable :: tracmin,tracmax
      integer, dimension(:), save, allocatable :: tracinterp,igashr
      real, dimension(:,:), save, allocatable :: tracdaytime
      integer, save :: numtracer, nhr
      integer, save :: unit_trout=131
     
! rml 16/02/10 additions for Transcom methane
      real, dimension(:,:,:), save, allocatable :: oh123,strloss123
      real, dimension(:,:), save, allocatable :: oh,strloss
! rml 30/04/10 additions for Transcom MCF
      real, dimension(:,:), save, allocatable :: mcfdep123,jmcf
      real, dimension(:), save, allocatable :: mcfdep
      real, dimension(:,:,:), save, allocatable :: jmcf123
      logical methane,mcf

      character(len=80), save :: tracerlist=''
      character(len=80), save :: sitefile=''
      character(len=80), save :: shipfile=''
      character(len=80), save :: trout='tracer.stdout'
      real,save :: tracvalin=-999  !default to initialise from restart

      namelist/trfiles/tracerlist,sitefile,shipfile,trout,tracvalin


      contains

c ***************************************************************************
      subroutine init_tracer
      use cc_mpi, only : myid
      implicit none
      character(len=80) :: header
      integer nt
      include 'newmpar.h'
      include 'tracers.h'
      character(len=80) :: tempname  ! Temp file name

!     if ( myid == 0 ) then
!     Allow any processor to write to this file (for testing purposes)
      open(unit=unit_trout,file=trout,form='formatted')
!     end if
c     first read in a list of tracers, then determine what emission data is required
!     Each processor read this, perhaps not necessary?

!     rml 16/2/10 addition for Transcom methane
      methane = .false.
      mcf = .false.

      open(unit=130,file=tracerlist,form='formatted')
      read(130,*) header
      read(130,*) numtracer
      if (numtracer.ne.ngas) 
     & stop 'wrong number of tracers in tracer.dat file'
      allocate(tracname(numtracer),tractype(numtracer))
      allocate(tracinterp(numtracer),tracunit(numtracer))
      allocate(tracfile(numtracer),igashr(numtracer))
      allocate(tracival(numtracer),trmass(numtracer))
!     rml 18/09/07 added tracmin, tracmax
      allocate(tracmin(numtracer),tracmax(numtracer))
      allocate(tracdaytime(numtracer,2))
      nhr = 0
      do nt=1,numtracer
!       rml 18/09/07 added tracmin,tracmax
        read(130,*) tracname(nt),tracival(nt),tracmin(nt),tracmax(nt),
     &              tractype(nt),tracfile(nt)
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
         if ( myid == 0 ) then
          write(unit_trout,999) 'Tracer ',nt,tracname(nt),tractype(nt),
     &                           tracfile(nt),tracinterp(nt)
        end if
 999    format(a7,i5,x,a13,a13,a50,i3)
! rml 16/2/10 addition for TC methane
        if (tracname(nt)(1:7).eq.'methane') methane = .true.
        if (tracname(nt)(1:3).eq.'mcf') mcf = .true.
c
      enddo
      allocate(nghr(nhr))

c     initialise array for monthly average tracer
      traver = 0.
!     initialise accumulated loss (for methane cases)
      acloss_g = 0.

      return
      end subroutine

c ***********************************************************************
      subroutine tracini
c     initial value now read from tracerlist 
      use cc_mpi
      implicit none
      include 'newmpar.h'
      include 'tracers.h'
      include 'netcdf.inc'
      include 'mpif.h'
      integer i,ierr
! rml 19/04/10 variables for methane initial condition 
      integer ok,ncid,ch4id,mcfid
      real ch4in_g(ifull_g,kl)
      real ch4in(ilt*jlt,klt)
      real mcfin_g(ifull_g,kl)
      real mcfin(ilt*jlt,klt)

      do i=1,ngas
c rml 15/11/06 facility to introduce new tracers to simulation
c i.e. some read from restart, others initialised here
        if (tracival(i).ne.-999) then
! rml 16/2/10 addition for TC methane to get 2d initial condition
            if (tracname(i)(1:7).eq.'methane') then
!             read  initial condition
              if (myid == 0) then
                ok = nf_open('ch4in_cc48.nc',0,ncid)
                if (ok.ne.0) write(6,*) nf_strerror(ok)
                ok = nf_inq_varid(ncid,'ch4in',ch4id)
                if (ok.ne.0) write(6,*) nf_strerror(ok)
                ok = nf_get_var_real(ncid,ch4id,ch4in_g)
                if (ok.ne.0) write(6,*) nf_strerror(ok)
                ok = nf_close(ncid)
                call ccmpi_distribute(ch4in,ch4in_g)
              else
                call ccmpi_distribute(ch4in)
              endif
              tr(1:ilt*jlt,1:klt,i)=ch4in
            elseif (tracname(i)(1:3).eq.'mcf') then
!             read  mcf initial condition
              if (myid == 0) then
                ok = nf_open('mcfin_cc48.nc',0,ncid)
                if (ok.ne.0) write(6,*) nf_strerror(ok)
                ok = nf_inq_varid(ncid,'mcfin',mcfid)
                if (ok.ne.0) write(6,*) nf_strerror(ok)
                ok = nf_get_var_real(ncid,mcfid,mcfin_g)
                if (ok.ne.0) write(6,*) nf_strerror(ok)
                ok = nf_close(ncid)
                call ccmpi_distribute(mcfin,mcfin_g)
              else
                call ccmpi_distribute(mcfin)
              endif
              tr(1:ilt*jlt,1:klt,i)=mcfin
            else
              tr(:,:,i)=tracival(i)
            endif
        endif
      enddo
!
!     if ( myid == 0 ) then
!       write(unit_trout,*) myid,'tracini: ',tracival
!       write(unit_trout,*) myid,'tracini: ',tr(1,1,:)
!       write(unit_trout,*) myid,'tracini1: ',tr(1,:,1)
!     end if

      if( iradon.ne.0.) then
        tr(:,:,max(1,iradon))=0.
      end if

      return
      end subroutine

c ***********************************************************************
      subroutine tr_back
c     remove a background value for tracer fields for more accurate transport
      use cc_mpi, only : myid
      implicit none
      include 'newmpar.h'
      include 'tracers.h'
      include 'mpif.h'
      integer i,ierr
      real trmin,trmin_g
      

      if ( myid == 0 ) write(unit_trout,*) 'Background tracer concentrat
     &ion removed:'

      do i=1,ngas
!       use minimum concentration from level in middle of atmosphere as background
        trmin=minval(tr(1:ilt*jlt,klt/2,i))
!       trmax=maxval(tr(1:ilt*jlt,klt/2,i))

        call MPI_Allreduce(trmin,trmin_g,1,MPI_REAL,MPI_MIN,
     &                     MPI_COMM_WORLD,ierr)
!       call MPI_Allreduce(trmax,trmax_g,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ier)
      
        trback_g(i) = trmin_g
        tr(:,:,i) = tr(:,:,i) - trback_g(i)
!
        if ( myid == 0 ) 
     &    write(unit_trout,*) 'Tracer ',i,' : ',trback_g(i)

      enddo

      return
      end subroutine

c *********************************************************************
      subroutine readtracerflux(kdate)
c     needs to happen further down the code than reading the tracer
c     list file
      use cc_mpi, only : myid
      implicit none
      include 'newmpar.h'
      include 'tracers.h'
      include 'parm.h'
      real ajunk(3)
      integer nt,jyear,jmonth,kdate
      character filename*50,varname*13
c
c     set up and read surface flux data
      allocate(co2em123(ilt*jlt,3,numtracer))
      allocate(co2emhr(ilt*jlt,31*24+2,nhr))
      allocate(co2hr(31*24+2,nhr))
      allocate(co2em(ilt*jlt,numtracer))

! rml 16/2/10 addition for TC methane, assume monthly OH and strat loss
      if (methane) then
        allocate(oh123(ilt*jlt,klt,3))
        allocate(strloss123(ilt*jlt,klt,3))
        allocate(oh(ilt*jlt,klt))
        allocate(strloss(ilt*jlt,klt))
      endif
! rml 30/4/10 addition for TC mcf, deposition rates
      if (mcf) then
        allocate(mcfdep123(ilt*jlt,3))
        allocate(mcfdep(ilt*jlt))
        allocate(jmcf123(ilt*jlt,klt,3))
        allocate(jmcf(ilt*jlt,klt))
      endif

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
!
!     just read OH, strat loss once regardless of how many methane tracers
!     filename set here at compile time
      if (methane) then
        filename = '/short/r39/TCinput/oh_c48.nc'
        varname = 'oh'
        call readoh(jmonth,3,filename,varname,oh123)
        filename = '/short/r39/TCinput/strloss_c48.nc'
        varname = 'strloss'
        call readoh(jmonth,3,filename,varname,strloss123)
      endif
      if (mcf) then
!       use standard tracer flux file read call but pass through 
!       with tracer 'ngas+1' - this will trigger MCF_loss as filename
        call readrco2(ngas+1,jyear,jmonth,3,mcfdep123,ajunk)
        filename='/short/r39/TCinput/Jmcf_cc48.nc'
        varname = 'jmcf'
        call readoh(jmonth,3,filename,varname,jmcf123)
      endif
      if (myid==0) write(unit_trout,*) 'Read input fluxes/rates etc OK'

      return
      end subroutine

c *************************************************************************
      subroutine readrco2(igas,iyr,imon,nflux,fluxin,co2time)
c     rml 23/09/03 largely rewritten to use netcdf files
      use cc_mpi
      implicit none
      include 'newmpar.h' !il,jl,kl
      include 'parm.h' !nperday
      include 'tracers.h'
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
      integer ncidfl,timedim,yearid,monthid,fluxid,hourid
      integer nregdim,nregion,dayid,ierr
      integer, dimension(:), allocatable :: fluxyr,fluxmon
      real, dimension(:), allocatable :: fluxhr
      integer start(3),count(3)
      real fluxin_g(ifull_g,nflux)
      include 'netcdf.inc'
      include 'mpif.h'

!  rml 30/04/10 special case for MCF deposition rates
      if (igas.eq.ngas+1) then
        fluxtype='monrep'
        fluxname='flux'
        filename = '/short/r39/TCinput/MCF_loss_CCAM48.nc'
      else
        fluxtype=tractype(igas)
        fluxname=tracname(igas)
        filename=tracfile(igas)
      endif
c
      if (trim(fluxtype).eq.'pulseoff'.or.
     &    trim(fluxtype).eq.'daypulseoff') then
c       no surface fluxes to read
        fluxin = 0.
        tracunit(igas)=''
        return
      end if

      if ( myid == 0 ) then ! Read on this processor and then distribute
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
! rml 30/4/10 exclude mcf deposition case
        if (igas.le.ngas) then
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
          if ( myid == 0 ) then
           write(unit_trout,*)'reading ',ncur,fluxyr(ncur),fluxmon(ncur)
          end if
          if (ncur.eq.0) stop 'current year/month not in flux file'
c    
          fluxin_g=0.
          start(1)=1
          count(1)=ifull_g; count(2)=1
c         read preceeding month if needed
          if (nprev.ne.0) then
            start(2)=nprev
!           if (igas==7)print*, "Reading prev", start(:2), count(2)
            ierr=nf_get_vara_real(ncidfl,fluxid,start,count,
     &                            fluxin_g(:,1))
            if (ierr.ne.nf_noerr) stop 'error reading fluxin prev'
          endif
c         read current month/year
          start(2)=ncur
!         if (igas==7)print*, "Reading curr", start(:2), count(2)
          ierr=nf_get_vara_real(ncidfl,fluxid,start,count,fluxin_g(:,2))
          if (ierr.ne.nf_noerr) stop 'error reading fluxin cur'
c         read next month
          if (nnext.ne.0) then
            start(2)=nnext
!         if (igas==7)print*, "Reading next", start(:2), count(2)
            ierr=nf_get_vara_real(ncidfl,fluxid,start,count,
     &                            fluxin_g(:,3))
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
          start(1)=1; count(1)=ifull_g
          ntot=n2-n1+1
          start(2)=n1; count(2)=ntot
          ierr=nf_get_vara_real(ncidfl,fluxid,start,count,
     &                                fluxin_g(:,2:ntot+1))
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
          ierr=nf_get_vara_real(ncidfl,fluxid,start,count,fluxin_g(:,1))
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
     &                                      fluxin_g(:,ntot+2))
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
     &                          tracdaytime(igas,:))
          if (ierr.ne.nf_noerr) stop 'error reading daylight '
        end if

        ierr = nf_close(ncidfl)

        ! Should be more careful here, probably don't need the full range of
        ! the second dimension all the time.
        call ccmpi_distribute(fluxin,fluxin_g)

      else ! myid /= 0
        call ccmpi_distribute(fluxin)
      end if !myid == 0

!     Simple broadcast for co2 time
      call MPI_Bcast(co2time,nflux,MPI_REAL,0,MPI_COMM_WORLD,ierr)
!     Also need to share tracunit, tractype, and tracname. MPI_character. Total length
!     of the array
      call MPI_Bcast(tracunit,13*numtracer,MPI_CHARACTER,0,
     &               MPI_COMM_WORLD,ierr)
      call MPI_Bcast(tractype,13*numtracer,MPI_CHARACTER,0,
     &               MPI_COMM_WORLD,ierr)
      call MPI_Bcast(tracname,13*numtracer,MPI_CHARACTER,0,
     &               MPI_COMM_WORLD,ierr)
      
      if (trim(fluxtype).eq.'daypulseon') then

        call MPI_Bcast(tracdaytime(igas,:),2,MPI_REAL,0,MPI_COMM_WORLD,
     &                 ierr)

c       count number of timesteps that source emitting for
        kount=0
        do n=1,nperday
          hr = 24.*float(n)/float(nperday)
          if (tracdaytime(igas,1).lt.tracdaytime(igas,2) .and.
     &        tracdaytime(igas,1).le.hr .and.
     &        tracdaytime(igas,2).ge.hr) kount=kount+1 
          if (tracdaytime(igas,1).gt.tracdaytime(igas,2) .and.
     &        (tracdaytime(igas,1).le.hr .or.
     &        tracdaytime(igas,2).ge.hr)) kount=kount+1 
        enddo
c       scale flux to allow for emission over fraction of day
c       just set flux to zero if no daylight
        if (kount.ne.0) then
          fluxin = fluxin*float(nperday)/float(kount)
        else
          fluxin = 0.
        endif
      endif

      end subroutine

c *************************************************************************
      subroutine readoh(imon,nfield,ohfile,varname,ohin)
! rml 16/2/10 New subroutine to read oh and strat loss for Transcom methane
      use cc_mpi
      implicit none
      include 'newmpar.h' !il,jl,kl
      include 'parm.h' !nperday
      character*50 ohfile
      character*13 varname
      integer nfield,imon,ntime
      integer nprev,nnext,ncur,n
c     nflux =3 for month interp case - last month, this month, next month
      real ohin(il*jl,kl,nfield)
      integer ncidfl,timedim,monthid,fluxid,ierr
      integer, dimension(:), allocatable :: ohmon
      integer start(3),count(3)
      real ohin_g(ifull_g,kl,nfield)
      include 'netcdf.inc'
      include 'mpif.h'

c

      if ( myid == 0 ) then ! Read on this processor and then distribute
        write(unit_trout,*)'reading for ',imon,' from ',ohfile
        ierr = nf_open(trim(ohfile),0,ncidfl)
        if (ierr.ne.nf_noerr) then
          write(unit_trout,*) ierr,ohfile
          write(unit_trout,*) nf_strerror(ierr)
          stop 'methane oh/loss file not found'
        endif
        ierr=nf_inq_dimid(ncidfl,'time',timedim)
        if (ierr.ne.nf_noerr) stop 'time dimension error'
        ierr=nf_inq_dimlen(ncidfl,timedim,ntime)
        if (ierr.ne.nf_noerr) stop 'time dimension length error'
        allocate(ohmon(ntime))
        ierr=nf_inq_varid(ncidfl,'month',monthid)
        if (ierr.ne.nf_noerr) stop 'month variable not found'
        ierr=nf_get_var_int(ncidfl,monthid,ohmon)
        if (ierr.ne.nf_noerr) stop 'month variable not read'
c       check fluxname 
        ierr=nf_inq_varid(ncidfl,trim(varname),fluxid)
        if (ierr.ne.nf_noerr) stop 'oh/loss variable not found'
c
c       find required records
        if (nfield.eq.3) then
c         monthly case
          nprev=0
          nnext=0
          ncur=0
          if (ntime.ne.12) stop 'oh/loss file wrong ntime'
          do n=1,ntime
            if (ohmon(n).eq.imon) then
              ncur=n
              nprev=n-1
              if (nprev.eq.0) nprev=12
              nnext=n+1
              if (nnext.eq.13) nnext=1
            endif
          enddo
c
          if ( myid == 0 ) then
           write(unit_trout,*)'reading ',ncur,ohmon(ncur)
          end if
          if (ncur.eq.0) stop 'current month not in flux file'
c    
          ohin_g=0.
          start(1)=1 ; start(2)=1
          count(1)=ifull_g; count(2)=kl ; count(3)=1 ! CHECK
c         read preceeding month if needed
          if (nprev.ne.0) then
            start(3)=nprev
            ierr=nf_get_vara_real(ncidfl,fluxid,start,count,
     &                            ohin_g(:,:,1))
            if (ierr.ne.nf_noerr) stop 'error reading ohin prev'
          endif
c         read current month/year
          start(3)=ncur
          ierr=nf_get_vara_real(ncidfl,fluxid,start,count,
     &                            ohin_g(:,:,2))
          if (ierr.ne.nf_noerr) stop 'error reading ohin cur'
c         read next month
          if (nnext.ne.0) then
            start(3)=nnext
            ierr=nf_get_vara_real(ncidfl,fluxid,start,count,
     &                            ohin_g(:,:,3))
            if (ierr.ne.nf_noerr) stop 'error reading ohin next'
          endif
        endif

        ierr = nf_close(ncidfl)

        ! Should be more careful here, probably don't need the full range of
        ! the second dimension all the time.
        do n=1,nfield
          call ccmpi_distribute(ohin(:,:,n),ohin_g(:,:,n))
        enddo

      else ! myid /= 0
        do n=1,nfield
          call ccmpi_distribute(ohin(:,:,n))
        enddo
      end if !myid == 0

      end subroutine


   
c *************************************************************************
      subroutine interp_tracerflux(kdate,hrs_dt)
c     interpolates tracer flux to current timestep if required
c     tracinterp 0 for no interpolation, 1 for monthly, 2 for daily/hourly
c     co2em123(:,1,:) contains prev month, co2em123(:,2,:) current month/year 
c     co2em123(:,3,:) next month
      use cc_mpi, only : myid
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
             if ( myid==0 ) then
               write(unit_trout,*) 'igas: ',igas,'igashr: ',igh
               write(unit_trout,*) 'hrmodel: ',hrmodel
             end if
             do while (.not.found)
               if ( myid==0 ) then
                 write(unit_trout,*) nghr(igh),co2hr(nghr(igh),igh),
     &                    co2hr(nghr(igh)+1,igh)
               end if
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

! rml 16/2/10 interpolate oh and strat loss for methane cases
      if (methane) then
!       monthly interpolation
        oh(:,:) = (1.0-ratlm)*oh123(:,:,m1) +
     &                 ratlm*oh123(:,:,m2)
        strloss(:,:) = (1.0-ratlm)*strloss123(:,:,m1) +
     &                 ratlm*strloss123(:,:,m2)
      endif
      if (mcf) then
        jmcf(:,:) = (1.0-ratlm)*jmcf123(:,:,m1) +
     &                 ratlm*jmcf123(:,:,m2)
        mcfdep(:) = (1.0-ratlm)*mcfdep123(:,m1) +
     &                 ratlm*mcfdep123(:,m2)
      endif

      return
      end subroutine

c ***************************************************************************
      subroutine tracer_mass(ktau,ntau)
c     rml 16/10/03 check tracer mass - just write out for <= 6 tracers
      use cc_mpi
      use sumdd_m
      implicit none
      include 'newmpar.h'
      include 'const_phys.h' ! rearth,fc_molm,fair_molm
      include 'arrays.h'   ! ps
      include 'sigs.h'     ! dsig
      include 'xyzinfo.h'  ! wts
      include 'tracers.h'  ! tr
      include 'mpif.h'
      integer it,iq,k,ktau,ntau,igas,ierr

      real trmin,trmax
      real ::trmass_l(ngas)
      real checkwts,checkwts_g,checkdsig,checkdsig_g
#ifdef sumdd
       complex :: local_sum(ngas), global_sum(ngas)
!      Temporary array for the drpdr_local function
       real, dimension(ifull) :: tmparr
#endif

       trmass_l = 0.
#ifdef sumdd
       local_sum = (0.,0.)
#endif
       do it=1,ngas
!         checkwts = 0.
!         checkdsig = 0.
          do k=1,kl
             do iq=1,ifull
#ifdef sumdd         
! rml 22/4/10 looks like bug *wts not ** wts
!               tmparr(iq)  = tr(iq,k,it)*dsig(k)*ps(iq)**wts(iq)
                tmparr(iq)  = tr(iq,k,it)*dsig(k)*ps(iq)*wts(iq)
#else
! rml 22/4/10 looks like bug *wts not ** wts
                trmass_l(it) = trmass_l(it) +
!    &                          tr(iq,k,it)*dsig(k)*ps(iq)**wts(iq)
!    &                          tr(iq,k,it)*dsig(k)*ps(iq)*wts(iq)
     &             (trback_g(it)+tr(iq,k,it))*dsig(k)*ps(iq)*wts(iq)
!               checkwts = checkwts + wts(iq)
!               checkdsig = checkdsig + dsig(k)
#endif
             end do
          end do
#ifdef sumdd
          call drpdr_local(tmparr, local_sum(it))
#endif
       end do ! it

#ifdef sumdd
       call MPI_Allreduce ( local_sum, global_sum, ngas, 
     &                     MPI_COMPLEX, MPI_SUMDR, MPI_COMM_WORLD, ierr)
       trmass = real(global_sum)
#else
       call MPI_Allreduce ( trmass_l, trmass, ngas, MPI_REAL,
     &                      MPI_SUM, MPI_COMM_WORLD, ierr )
!      call MPI_Allreduce ( checkwts, checkwts_g, 1, MPI_REAL,
!    &                      MPI_SUM, MPI_COMM_WORLD, ierr )
!      call MPI_Allreduce ( checkdsig, checkdsig_g, 1, MPI_REAL,
!    &                      MPI_SUM, MPI_COMM_WORLD, ierr )
#endif


c     scaling assumes CO2 with output in GtC?
      if ( myid == 0 ) then
         if (ngas.gt.6) then
            write(unit_trout,*) 'Trmass: ',ktau,
     &        -1*trmass(1:6)*4.*3.14159*(rearth**2)*fC_MolM/
     &          (grav*1e18*fAIR_MolM)
         else
!           write(unit_trout,*) 'Trmass: ',ktau,
!    &        -1*trmass(:)*4.*3.14159*(rearth**2)*fC_MolM/
!    &         (grav*1e18*fAIR_MolM)
!    & ,minval(tr(1:ilt*jlt,1:klt,1)),maxval(tr(1:ilt*jlt,1:klt,1))
!           scaling for methane in Tg
!           write(unit_trout,*) 'Trmass ppb*Pa: ',trmass(:)
!           write(unit_trout,*) checkwts_g,checkdsig_g
            write(unit_trout,*) 'Trmass (Tg CH4): ',ktau,
     &        -1*trmass(:)*4.*pi*eradsq*fCH4_MolM/
     &         (grav*1.e18*fAIR_MolM)
         endif
      end if

c     also update tracer average array here
      do igas=1,ngas

        traver(:,:,igas)=traver(:,:,igas)+tr(1:ilt*jlt,1:klt,igas)/ntau

!       rml 18/09/07 check that tr and traver stay within defined range
!       stop job if they don't
!  rml 22/2/10 only check traver at end of month
!       separate checks on each processor
        trmin = minval(tr(1:ilt*jlt,1:klt,igas))+trback_g(igas)
        trmax = maxval(tr(1:ilt*jlt,1:klt,igas))+trback_g(igas)
        if (trmin.lt.tracmin(igas)) then
          write(6,*) 'WARNING: below minimum, tracer ',igas, 
     &    ' processor ',myid,tracmin(igas),trmin,
     &   minloc(tr(1:ilt*jlt,1:klt,igas))
          write(6,*) 'Error: tracer out of range.  See tracer.stdout'
          stop
        endif
        if (trmax.gt.tracmax(igas)) then
          write(6,*) 'WARNING: above maximum, tracer ',igas, 
     &    ' processor ',myid,tracmax(igas),trmax,
     &   maxloc(tr(1:ilt*jlt,1:klt,igas))
          write(6,*) 'Error: tracer out of range.  See tracer.stdout'
          stop
        endif
        if (ktau.eq.ntau) then
          trmin = minval(traver(1:ilt*jlt,1:klt,igas))+trback_g(igas)
          trmax = maxval(traver(1:ilt*jlt,1:klt,igas))+trback_g(igas)
          if (trmin.lt.tracmin(igas)) then
            write(6,*) 'Error: tracer out of range.  See tracer.stdout'
            write(6,*) 'WARNING: trav below minimum, tracer ',
     & igas, ' processor ',myid,tracmin(igas),trmin
          stop
          endif
          if (trmax.gt.tracmax(igas)) then
            write(6,*) 'Error: tracer out of range.  See tracer.stdout'
            write(6,*) 'WARNING: trav above maximum, tracer ',
     & igas, ' processor ',myid, tracmax(igas),trmax
          stop
          endif
        endif

      enddo


      end subroutine

      end module
