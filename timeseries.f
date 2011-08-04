      module timeseries

c     rml 25/08/03 declarations from sflux
      implicit none
      integer, save :: indextime,ntsfreq,ngrdpts,ngrdpts1,n3d,n2d
      double precision, save :: tstime
      integer, pointer, dimension(:,:), save :: listijk
      logical, pointer, dimension(:), save :: writesurf
      integer, pointer, dimension(:), save :: tsid
      character*10, pointer, dimension(:), save :: varname3,varname2
c     rml 25/11/03 declarations from sflux
      integer, save :: indship,nshippts,inshipid(4),outshipid(3)
      integer, save :: nshipout
      integer, pointer, dimension(:), save :: shipdate,shiptime
!     rml 19/09/07 conversion of u,v from ccam grid to zonal/merid
      real, dimension(:), save, allocatable :: costh,sinth


      contains
c *********************************************************************

      subroutine init_ts(ngas,dt)
      use tracermodule, only : sitefile,shipfile
      implicit none
      integer ngas
      real dt

      if (sitefile.ne.'') call readsitelist(ngas)
      if (shipfile.ne.'') call readshiplist(ngas,dt)

      return
      end subroutine

c ********************************************************************
      subroutine write_ts(ktau,ntau,dt)
      use tracermodule, only : sitefile,shipfile
      implicit none
      include 'dates.h'
      integer jyear,jmonth,jday,jhour,jmin
      integer ktau,ntau,mstart,mins
      real dt
      integer ndoy(12)   ! days from beginning of year (1st Jan is 0)
      data ndoy/ 0,31,59,90,120,151,181,212,243,273,304,334/

      jyear=kdate/10000
      jmonth=(kdate-jyear*10000)/100
      jday=kdate-jyear*10000-jmonth*100
      jhour=ktime/100
      jmin=ktime-jhour*100
      mstart=24*60*(ndoy(jmonth)+jday-1) + 60*jhour + jmin ! mins from start of year
!     mtimer contains number of minutes since the start of the run.
      mins = mtimer + mstart

c   rml 25/08/03 write tracer data to timeseries file
      if (sitefile.ne.'') call writetimeseries(ktau,ntau,jyear,mins)
c     rml 26/11/03 write mobile tracer data to file
      if (shipfile.ne.'') call writeshipts(ktau,ntau,dt)

      return

      end subroutine
c *********************************************************************
      subroutine readsitelist(ntrac)
c
c     rml 25/08/03 subroutine to read file containing list of sites for 
c     timeseries output and to open netcdf file for output and 
c     write dimensions etc.
c     All processors read this list and select the points in their own
c     region.
c
      use cc_mpi, only : myid, indv_mpi, fproc, ipan
!     rml 19/09/07 add tracname so that can be written to ts output file
      use tracermodule, only : sitefile,tracname
      use vecsuv_m
      use xyzinfo_m
      implicit none
      integer kount,kountprof,n,kount500
      integer ierr,ntrac
      integer griddim,ijkdim,timedim,tracdim,gridid,dims(3)
!     rml 19/09/07 addition of tracer name array in output file
      integer lendim, tracnamid
      integer surfdim,gridsurfid
      integer gridorderid, surforderid
      integer, allocatable, dimension(:,:) :: templist
      integer, allocatable, dimension(:) :: gridorder, surforder
      integer i,i1,k,ntop
      character*20 outfile
      character*8 chtemp
      character*80 head
      integer :: ig, jg, tmpval ! Better name for this
      include 'netcdf.inc'
      include 'dates.h'
      include 'newmpar.h'  ! kl
      integer ip, nface, ii, jj, iqg, istn, jstn
!     rml 19/09/07 arrays needed for wind conversion
      include 'parmgeom.h'    ! rlong0, rlat0
      include 'const_phys.h'  ! pi
      logical windconv
      real coslong,sinlong,coslat,sinlat,polenx,poleny,polenz
      real zonx,zony,zonz,den
      integer iq

c     read file of site locations for timeseries output
      open(88,file=sitefile,form='formatted', status='unknown')
      read(88,*) head
c     number of gridpoints and output frequency (number of timesteps)
      read(88,*) ngrdpts1,ntsfreq
      allocate(templist(ngrdpts1,3))
      if ( nproc > 1 ) then
        allocate(surforder(ngrdpts1) )
      end if
      kountprof=0
      kount500=0

      ! Read list of point locations, selecting those that are in
      ! this processor's region.
      ! Perhaps need to have an additional netcdf variable for ip, so
      ! that the original order can be reconstructed from the multiple files.
      k = 0
      do ip=1,ngrdpts1
        read(88,*) i,i1,ig,jg,tmpval
        ! Convert to local indices. Same code used in indata for stations.
        ! Should be generalised to a routine?
        nface=(jg-1)/il_g
!       Note that the second argument to fproc is the j index on the
!       face, not the global j index,   
        if ( fproc(ig,jg - nface*il_g,nface) == myid ) then
           ! Point is in my region
           iqg = ig + (jg-1)*il_g
           ! Local indices on this processor
           call indv_mpi(iqg,ii,jj,n)
           k = k + 1
           istn = ii
           jstn = jj+(n-1)*ipan
           templist(k,:) = (/ istn, jstn, tmpval /)
c          check if any profiles requested
           if (templist(k,3).eq.99) kountprof=kountprof+1
           if (templist(k,3).eq.98) kount500=kount500+1
           ! Define order variable to allow merging the separate processor files
           if ( nproc > 1 ) surforder(k) = ip
        end if
      enddo

      ngrdpts1 = k ! Reset to the actual number I have

c     Read in any additional variables to output besides tracer
      n2d=0
      n3d=0
      read(88,*,end=880) head
      read(88,*) n3d
      allocate(varname3(n3d))
      windconv = .false.
      do n=1,n3d
        read(88,*) varname3(n)
!       rml 19/09/07 check if wind conversion required
        if (trim(varname3(n)).eq.'u'.or.trim(varname3(n)).eq.'v')
     &           windconv=.true.
      enddo
      read(88,*) head
      read(88,*) n2d
      allocate(varname2(n2d))
      do n=1,n2d
        read(88,*) varname2(n)
      enddo
 880  continue
      close(88)

!     rml 19/09/07 fill arrays needed for wind conversion
      if (windconv) then
        allocate(costh(ifull),sinth(ifull))
!       code taken from cc2hist 
        coslong=cos(rlong0*pi/180.)
        sinlong=sin(rlong0*pi/180.)
        coslat=cos(rlat0*pi/180.)
        sinlat=sin(rlat0*pi/180.)
        polenx=-coslat
        poleny=0.
        polenz=sinlat
        do iq=1,ifull
!         Set up unit zonal vector components
          zonx = poleny*z(iq)-polenz*y(iq)
          zony = polenz*x(iq)-polenx*z(iq)
          zonz = polenx*y(iq)-poleny*x(iq)
!         Allow for poles by taking max
          den = sqrt( max(zonx**2 + zony**2 + zonz**2,1.e-7) )
          costh(iq) =  (zonx*ax(iq)+zony*ay(iq)+zonz*az(iq))/den
          sinth(iq) = -(zonx*bx(iq)+zony*by(iq)+zonz*bz(iq))/den
        enddo
      endif

      ngrdpts = ngrdpts1 + kountprof*(kl-1) + kount500*8
      allocate(listijk(ngrdpts,3))
      allocate(writesurf(ngrdpts))
      if ( nproc > 1 ) then
        allocate(gridorder(ngrdpts) )
      end if
      kount = 0
      do k=1,ngrdpts1
c     rml 11/11/05 add '98' option for all levels to ~500 hPa
        if (templist(k,3).ge.98) then
          if (templist(k,3).eq.98) then
            ntop = 9
          else
            ntop = kl
          endif
          do n=1,ntop
            kount = kount + 1
            listijk(kount,1:2) = templist(k,1:2)
            listijk(kount,3) = n
            if (n.eq.1) then
              writesurf(kount) = .true.
            else
              writesurf(kount) = .false.
            endif
            if ( nproc > 1 ) gridorder(kount) = surforder(k)
          enddo
        else
          kount = kount + 1
          listijk(kount,:) = templist(k,:)
          writesurf(kount) = .true.
          if ( nproc > 1 ) gridorder(kount) = surforder(k)
        endif
      enddo
      if (kount.ne.ngrdpts) stop 'location file: kount.ne.ngrdpts'
c     deallocate(templist)
      allocate(tsid(3+n3d+n2d))
c
c     open netcdf file for writing output
      write(chtemp,'(i8)') kdate
      if ( nproc == 1 ) then
        outfile = 'ts.'//chtemp(1:4)//'.'//chtemp(5:6)//'.nc'
      else
!       Include a 2 digid processor number in the file
        write(outfile,'(a,a,a,a,a, i2.2,a)') 'ts.', chtemp(1:4), '.',
     &     chtemp(5:6), '.p', myid, '.nc'
      end if
!
!     when multiprocessor, possible to have no grid points in an output file
!     in that case, don't write the file
      if (ngrdpts.gt.0) then
        ierr = nf_create(outfile,0,tsid(1))
        if (ierr.ne.nf_noerr) stop 'create ts file failed'
c       define dimensions
        ierr=nf_def_dim(tsid(1),'gridpts',ngrdpts,griddim)
        if (ierr.ne.nf_noerr) stop 'timeseries: grid dimension error'
        ierr=nf_def_dim(tsid(1),'surfpts',ngrdpts1,surfdim)
        if (ierr.ne.nf_noerr) stop 'timeseries: grid dimension error'
        ierr=nf_def_dim(tsid(1),'ijk',3,ijkdim)
        ierr=nf_def_dim(tsid(1),'tracers',ntrac,tracdim)
        if (ierr.ne.nf_noerr) stop 'timeseries: tracer dimension error'
        ierr=nf_def_dim(tsid(1),'time',nf_unlimited,timedim)
        if (ierr.ne.nf_noerr) stop 'timeseries: time dimension error'
!     rml 19/09/07 add length dimension for tracname character array
        ierr = nf_def_dim(tsid(1),'len',13,lendim)
        if (ierr.ne.nf_noerr) stop 'timeseries: length dimension error'

c       define variables
!       rml 19/09/07 add tracer name variable
        dims(1)=lendim; dims(2)=tracdim
        ierr = nf_def_var(tsid(1),'tracer_name',nf_char,2,dims,
     &tracnamid)
        if (ierr.ne.nf_noerr) stop 'timeseries: tracer name var error'

        dims(1)=griddim; dims(2)=ijkdim
        ierr = nf_def_var(tsid(1),'grid',nf_int,2,dims,gridid)
        if (ierr.ne.nf_noerr) stop 'timeseries: grid var error'
        dims(1)=surfdim; dims(2)=ijkdim
        ierr = nf_def_var(tsid(1),'gridsurf',nf_int,2,dims,gridsurfid)
        if (ierr.ne.nf_noerr) stop 'timeseries: grid var error'

        if ( nproc > 1 ) then
          dims(1)=griddim
          ierr = nf_def_var(tsid(1),'gridorder',nf_int,1,dims,
     &                    gridorderid)
          if (ierr.ne.nf_noerr) stop 'timeseries: grid var error'
          dims(1)=surfdim
          ierr = nf_def_var(tsid(1),'surforder',nf_int,1,dims,
     & surforderid)
          if (ierr.ne.nf_noerr) stop 'timeseries: grid var error'
        end if

        ierr = nf_def_var(tsid(1),'time',nf_double,1,timedim,tsid(2))
        if (ierr.ne.nf_noerr) stop 'timeseries: tstime var error'
        dims(1)=griddim; dims(2)=tracdim; dims(3)=timedim
        ierr = nf_def_var(tsid(1),'concts',nf_float,3,dims,tsid(3))
        if (ierr.ne.nf_noerr) stop 'timeseries: concts var error'
        do n=1,n3d
          dims(1)=griddim; dims(2)=timedim
          ierr = nf_def_var(tsid(1),varname3(n),nf_float,2,dims(1:2),
     & tsid(3+n))
          if (ierr.ne.nf_noerr) stop 'timeseries: 3d var error'
        enddo
        do n=1,n2d
          if (trim(varname2(n)).eq.'flux') then
            dims(1)=surfdim; dims(2)=tracdim; dims(3)=timedim
            ierr = nf_def_var(tsid(1),varname2(n),nf_float,3,dims,
     & tsid(3+n3d+n))
            if (ierr.ne.nf_noerr) stop 'timeseries: 2d var error'
          else
            dims(1)=surfdim; dims(2)=timedim
            ierr = nf_def_var(tsid(1),varname2(n),nf_float,2,dims(1:2),
     &   tsid(3+n3d+n))
            if (ierr.ne.nf_noerr) stop 'timeseries: 2d var error'
          endif
        enddo
c
c     leave define mode
        ierr = nf_enddef(tsid(1))
        if (ierr.ne.nf_noerr) stop 'timeseries: end define error'
c
!     rml 19/09/07 write tracer name array
        ierr = nf_put_var_text(tsid(1),tracnamid,tracname)
        if (ierr.ne.nf_noerr) stop 'timeseries: error writing tracname'

c     write grid point arrays
        ierr = nf_put_var_int(tsid(1),gridid,listijk)
        if (ierr.ne.nf_noerr) stop 'error writing grid'
      ! Need explicit section of templist here, because array may have
      ! been allocated larger
        ierr = nf_put_var_int(tsid(1),gridsurfid,templist(:ngrdpts1,:))
        if (ierr.ne.nf_noerr) stop 'error writing gridsurf'
        if ( nproc > 1 ) then
         ierr = nf_put_var_int(tsid(1),gridorderid,gridorder(1:ngrdpts))
          if (ierr.ne.nf_noerr) stop 'error writing gridorder'
        ierr = nf_put_var_int(tsid(1),surforderid,surforder(1:ngrdpts1))
          if (ierr.ne.nf_noerr) stop 'error writing surforder'
        end if
        ierr = nf_sync(tsid(1))
        deallocate(templist)
      endif
c
      indextime=1
      return
      end subroutine
c
c **********************************************************************
   
      subroutine writetimeseries(ktau,ntau,jyear,mins)

c
c     rml: subroutine to write timeseries data to netcdf file
c  rml 10/11/05: added pressure, surface flux and pblh for TC
c
      use arrays_m    ! temp, q, ps
      use carbpools_m ! cbm co2 fluxes
      use define_dimensions, only : ncs, ncp ! Used in carbpool.h
      use extraout_m  ! cloud arrays
      use morepbl_m   ! rnet,eg,fg
      use pbl_m       ! tss
      use prec_m      ! precip
      use sigs_m      ! sigma levels for pressure
      use soil_m      ! albedo
      use soilsnow_m  ! soil temp (tgg)
      use tracermodule, only : co2em,unit_trout
      use tracers_m   ! ntrac and tr array
      use vegpar_m    ! rlai
      use vvel_m      ! vertical velocity
      implicit none
      real, dimension(:,:), allocatable :: cts
      real, dimension(:), allocatable :: vts
      integer ierr,start(3),count(3),n,iq,kount,m,jyear,mins
      integer ktau,ntau,k
      logical surfflux
      include 'netcdf.inc'
!!!      include 'cbmdim.h'
      include 'newmpar.h'    ! dimensions for tr array
      include 'const_phys.h' ! grav (for height calc)

      real temparr2(il*jl,kl),temparr(il*jl)

      if (ngrdpts.eq.0) return
      if (mod(ktau,ntsfreq).eq.0) then
        tstime = dble(jyear) + dble(mins)/dble(365.*24.*60.)
        ierr = nf_put_var1_double(tsid(1),tsid(2),indextime,tstime)
        if (ierr.ne.nf_noerr) stop ': error writing tstime'
        allocate(cts(ngrdpts,ntrac))
        do n=1,ngrdpts
          iq = listijk(n,1) + (listijk(n,2)-1)*il
! rml 23/2/10 add in background that removed at start of run
          cts(n,:)=trback_g(:)+tr(iq,listijk(n,3),:)
        enddo
        start(1)=1; start(2)=1; start(3)=indextime
        count(1)=ngrdpts; count(2)=ntrac; count(3)=1
        ierr=nf_put_vara_real(tsid(1),tsid(3),start,count,cts)
        if (ierr.ne.nf_noerr) then
          write(6,*) nf_strerror(ierr)
          stop 'error writing cts - tr array'
        endif
        deallocate(cts)
c
        do m=1,n3d
          select case(trim(varname3(m)))
          case ('t') ; temparr2=t(1:ifull,:)
!        rml 19/09/07 fix for u, v - convert to zonal/meridional from ccam grid
          case ('u') 
            do k=1,kl
              temparr2(:,k) = costh*u(1:ifull,k) - sinth*v(1:ifull,k)
            enddo
          case ('v') 
            do k=1,kl
              temparr2(:,k) = sinth*u(1:ifull,k) + costh*v(1:ifull,k)
            enddo
! rml 16/7/10 add output of height of levels for TC-methane
          case ('zg')
            temparr2(:,1)= bet(1)*t(1:ifull,1)/grav        
            do k=2,kl                                                   
              temparr2(:,k)=temparr2(:,k-1)+(bet(k)*t(1:ifull,k)
     &                     +betm(k)*t(1:ifull,k-1))/grav               
            enddo                                                       
            do k=1,kl
              temparr2(:,k) = temparr2(:,k) + zs(1:ifull)/grav
            enddo
          case ('qg') ; temparr2=qg(1:ifull,:)
          case ('sdotm') ; temparr2=sdot(:,1:kl)
          case ('sdotp') ; temparr2=sdot(:,2:kl+1)
          case ('pressure')
            do k=1,kl
              temparr2(:,k)=ps(1:ifull)*sig(k)
            enddo
          case default 
            write(unit_trout,*) varname3(m),' not found'
            stop
          end select
          allocate(vts(ngrdpts))
          do n=1,ngrdpts
            iq = listijk(n,1) + (listijk(n,2)-1)*il
            vts(n)=temparr2(iq,listijk(n,3))
          enddo
          start(1)=1; start(2)=indextime
          count(1)=ngrdpts; count(2)=1
          ierr=nf_put_vara_real(tsid(1),tsid(3+m),start(1:2),count(1:2),
     &                          vts)
          if (ierr.ne.nf_noerr) stop 'error writing vts'
          deallocate(vts)
        enddo
c
        do m=1,n2d
          temparr=0.
          surfflux=.false.
          select case(trim(varname2(m)))
          case ('cloudlo') ; temparr=cloudlo
          case ('cloudmi') ; temparr=cloudmi
          case ('cloudhi') ; temparr=cloudhi
          case ('ps')      ; temparr=ps(1:ifull)
          case ('tss')     ; temparr=tss
          case ('rnet')    ; temparr=rnet
          case ('eg')      ; temparr=eg
          case ('fg')      ; temparr=fg
!         case ('alb')     ; temparr=alb
          case ('alb')     ; temparr=swrsave*albvisnir(:,1)+
     &                               (1.-swrsave)*albvisnir(:,2) ! MJT cable
          case ('sgsave')  ; temparr=sgsave
          case ('rgsave')  ; temparr=rgsave
          case ('precip')  ; temparr=precip
          case ('tgg4')    ; temparr=tgg(:,4)
          case ('tgg5')    ; temparr=tgg(:,5)
          case ('tgg6')    ; temparr=tgg(:,6)
          case ('rlai')    ; temparr=rlai
          case ('pfnee')   ; temparr=fnee
          case ('pfpn')    ; temparr=fpn
          case ('pfrp')    ; temparr=frp
          case ('pfrs')    ; temparr=frs
          case ('pblh') ; temparr=pblh
          case ('flux')  
            allocate(cts(ngrdpts1,ntrac))
            kount=0
            do n=1,ngrdpts
              if (writesurf(n)) then
                 kount=kount+1
                 iq = listijk(n,1) + (listijk(n,2)-1)*il
                 cts(kount,:)=co2em(iq,:)
              endif
            enddo 
            start(1)=1; start(2)=1; start(3)=indextime
            count(1)=ngrdpts1; count(2)=ntrac; count(3)=1
            ierr=nf_put_vara_real(tsid(1),tsid(3+n3d+m),start,count,cts)
            if (ierr.ne.nf_noerr) then
              write(6,*) nf_strerror(ierr)
              stop 'error writing cts - fluxes'
            endif
            deallocate(cts)
            surfflux=.true.
          case default
            write(unit_trout,*) trim(varname2(m)),' not found'
            stop
          end select
          
          if (.not.surfflux) then
            allocate(vts(ngrdpts1))
            kount = 0
            do n=1,ngrdpts
              if (writesurf(n)) then
                kount = kount+1
                iq = listijk(n,1) + (listijk(n,2)-1)*il
                vts(kount)=temparr(iq)
              endif
            enddo
            start(1)=1; start(2)=indextime
            count(1)=ngrdpts1; count(2)=1
            ierr=nf_put_vara_real(tsid(1),tsid(3+n3d+m),start(1:2),
     &                            count(1:2),vts)
            if (ierr.ne.nf_noerr) stop 'error writing vts'
            deallocate(vts)
          endif
        enddo

c
        ierr=nf_sync(tsid(1))
        indextime=indextime+1
      endif
      if (ktau.eq.ntau) ierr=nf_close(tsid(1))
c

      return
      end subroutine
c *********************************************************************
      subroutine readshiplist(ntrac,dt)
c
c     rml 25/11/03 subroutine to read file containing times and locations
c     of ship (or other mobile obs).  Also opens netcdf file for output.
c
      use tracermodule, only : shipfile
      implicit none
      integer ok,nptsdim,dateid,timeid,ntrac,i,i2
      integer dtlsdim,nvaldim,outshipid(3),dims(2),tracdim
      real dt
      character*15 outfile2
      character*8 chtemp
      include 'newmpar.h'
      include 'netcdf.inc'
      include 'dates.h'
c
      if ( nproc > 1 ) then
         print*, "Error, parallel version of shiplist not yet working"
         stop
      end if
c     open file with ship locations
      ok = nf_open(shipfile,0,inshipid(1))
      if (ok.ne.nf_noerr) stop 'readshiplist: open file failure'
      ok = nf_inq_dimid(inshipid(1),'npts',nptsdim)
      if (ok.ne.nf_noerr) stop 'readshiplist: read nptsdim failure'
      ok = nf_inq_dimlen(inshipid(1),nptsdim,nshippts)
      if (ok.ne.nf_noerr) stop 'readshiplist: read dimlen failure'
c     read times for ship samples
      allocate(shipdate(nshippts),shiptime(nshippts))
      ok = nf_inq_varid(inshipid(1),'date',dateid)
      if (ok.ne.nf_noerr) stop 'readshiplist: read dateid failure'
      ok = nf_get_var_int(inshipid(1),dateid,shipdate)
      if (ok.ne.nf_noerr) stop 'readshiplist: read shipdate failure'
      ok = nf_inq_varid(inshipid(1),'time',timeid)
      if (ok.ne.nf_noerr) stop 'readshiplist: read timeid failure'
      ok = nf_get_var_int(inshipid(1),timeid,shiptime)
      if (ok.ne.nf_noerr) stop 'readshiplist: read shiptime failure'
      ok = nf_inq_varid(inshipid(1),'loc',inshipid(2))
      if (ok.ne.nf_noerr) stop 'readshiplist: read locid failure'
      ok = nf_inq_varid(inshipid(1),'lev',inshipid(4))
      if (ok.ne.nf_noerr) stop 'readshiplist: read locid failure'
      ok = nf_inq_varid(inshipid(1),'ship',inshipid(3))
      if (ok.ne.nf_noerr) stop 'readshiplist: read shipid failure'
c
c     locate sample time just smaller than current time
      do i=1,nshippts
        if (shipdate(i).lt.kdate) then
          indship=i
        endif
      enddo
      i2=indship
      do i=i2+1,nshippts
        if (shipdate(i).eq.kdate.and.
     &    shiptime(i).lt.ktime+nint(dt)/120) then !half interval to next ktime
          indship=i
        endif
      enddo
c
c     open netcdf file for output  
      write(chtemp,'(i8)') kdate
      outfile2 = 'ship.'//chtemp(1:4)//'.'//chtemp(5:6)//'.nc'
      ok = nf_create(outfile2,0,outshipid(1))
      if (ok.ne.nf_noerr) stop 'readshiplist: create error'
      ok=nf_def_dim(outshipid(1),'nshiptrac',ntrac,tracdim)
      if (ok.ne.nf_noerr) stop 'readshiplist: tracer dimension error'
c    rml 18/12/03 increased dimension to include level for aircraft output
      ok = nf_def_dim(outshipid(1),'date_time_loc_ship',5,dtlsdim)
      if (ok.ne.nf_noerr) stop 'readshiplist: dtls dimension error'
      ok = nf_def_dim(outshipid(1),'nval',nf_unlimited,nvaldim)
      if (ok.ne.nf_noerr) stop 'readshiplist: npt dimension error'
      dims(1)=dtlsdim; dims(2)=nvaldim
      ok = nf_def_var(outshipid(1),'shipinfo',nf_int,2,dims,
     &                outshipid(2))
      if (ok.ne.nf_noerr) stop 'readshiplist: shipinfo variable error'
      dims(1)=tracdim; dims(2)=nvaldim
      ok = nf_def_var(outshipid(1),'tsship',nf_float,2,dims,
     &                outshipid(3))
      if (ok.ne.nf_noerr) stop 'readshiplist: shiptracer variable error'
      ok = nf_enddef(outshipid(1))
      if (ok.ne.nf_noerr) stop 'readshiplist: end define error'

      nshipout=0
      return
      end subroutine
c ********************************************************************
      subroutine writeshipts(ktau,ntau,dt)
c
c     rml 25/11/03 subroutine to write mobile timeseries e.g. ship
c
      use tracers_m
      implicit none
      integer ktau,ntau,ok,iloc,ilev,iship,ierr
      integer jdate1,jdate2,jtime1,jtime2,mon
      real dt
      integer start(2),kount(2),info(5)
      logical moredat,found
      integer monlen(12)
      data monlen/31,28,31,30,31,30,31,31,30,31,30,31/
      
      include 'newmpar.h'    ! dimensions for tr array
      include 'netcdf.inc'
      include 'dates.h'
c
c     if reached end of data leave subroutine
      if (indship+1.gt.nshippts) return

c     assume always running in one month blocks so don't worry about
c     having to increment across months

      mon=(kdate-10000*(kdate/10000))/100
      if (mtimer.eq.monlen(mon)*60*24) then
c       end of month case
        jdate2 = kdate + 100
        jdate1 = kdate + monlen(mon)-1
      else
        jdate2 = kdate + mtimer/1440
        jdate1 = jdate2-1
      endif
      jtime1 = mod(mtimer-nint(dt)/120,1440)
      jtime1 = 100*(jtime1/60) + mod(jtime1,60)
      jtime2 = mod(mtimer+nint(dt)/120,1440)
      jtime2 = 100*(jtime2/60) + mod(jtime2,60)
c     check if sample time in current timestep (from jtime to jtime+dt)
c     assume ktime+dt will never force increment to kdate
      moredat=.true.
      do while (moredat)
        found=.false.
        if (jtime1.lt.jtime2) then
          if (shipdate(indship+1).eq.jdate2.and.
     &      shiptime(indship+1).ge.jtime1.and.
     &      shiptime(indship+1).lt.jtime2) found=.true.
        else
c         end of day case
          if ((shipdate(indship+1).eq.jdate1.and.
     &      shiptime(indship+1).ge.jtime1).or.
     &      (shipdate(indship+1).eq.jdate2.and.
     &      shiptime(indship+1).lt.jtime2)) found=.true.
        endif
        if (found) then
          nshipout = nshipout+1
c         keep real sample time rather than model time for easier
c         match to data
          info(1) = shipdate(indship+1)
          info(2) = shiptime(indship+1)
          ok = nf_get_var1_int(inshipid(1),inshipid(2),indship+1,iloc)
          if (ok.ne.nf_noerr) stop 'writeshipts: read loc error'
          info(3) = iloc
c    rml 18/12/03 addition of level info to do aircraft output
          ok = nf_get_var1_int(inshipid(1),inshipid(4),indship+1,ilev)
          if (ok.ne.nf_noerr) stop 'writeshipts: read lev error'
          info(4) = ilev
          ok = nf_get_var1_int(inshipid(1),inshipid(3),indship+1,iship)
          if (ok.ne.nf_noerr) stop 'writeshipts: read ship error'
          info(5) = iship
          start(1)=1; start(2)=nshipout
          kount(1)=5; kount(2)=1
          ok = nf_put_vara_int(outshipid(1),outshipid(2),start,kount,
     &                         info)
          if (ok.ne.nf_noerr) stop 'writeshipts: write info error'
          start(1)=1; start(2)=nshipout
          kount(1)=ntrac; kount(2)=1
          ok = nf_put_vara_real(outshipid(1),outshipid(3),start,kount,
     &                         tr(iloc,ilev,:))
          if (ok.ne.nf_noerr) stop 'writeshipts: write shipts error'

          indship=indship+1
        else
          moredat=.false.
        endif
      enddo
      ierr=nf_sync(outshipid(1))
      if (ktau.eq.ntau) then
        ok=nf_close(inshipid(1))
        ierr=nf_close(outshipid(1))
      endif
 
      return
      end subroutine

      end module
