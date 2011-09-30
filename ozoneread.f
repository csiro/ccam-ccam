      !--------------------------------------------------------------
      ! MJT - radiation
      ! Convert this subroutine to a module, so as to pass allocatable
      ! arrays to o3set.  These arrays are used for the CMIP5 Ozone
      ! datasets.
      module ozoneread

      implicit none

      private      
      public o3_read,o3set,fieldinterpolate,o3regrid
      
      integer, save :: ii,jj,kk
      real, dimension(:,:), allocatable, save :: o3pre,o3mth,o3nxt
      real, dimension(:), allocatable, save :: o3pres
      real, dimension(:,:), allocatable, save :: dduo3n,ddo3n2
      real, dimension(:,:), allocatable, save :: ddo3n3,ddo3n4
      
      contains
      !--------------------------------------------------------------
            
      subroutine o3_read(sigma,jyear,jmonth)
c
c     Reads the ozone data from the o3_datafile (filename set in namelist)
c
      use cc_mpi, only : myid ! MJT read
      use infile, only : ncmsg

      implicit none

      include 'newmpar.h'
      include 'filnames.h'   ! MJT radiation
      include 'netcdf.inc'   ! MJT radiation
      include 'mpif.h'       ! MJT read
      integer, intent(in) :: jyear,jmonth
      integer nlev,i,k,ierr
      integer ncstatus,ncid,tt
      integer valident,yy,mm,iti,nn
      integer, dimension(4) :: spos,npos
      real, dimension(:,:,:,:), allocatable :: o3dum
      real, dimension(:), allocatable :: o3lon,o3lat
      character*32 cdate
      real, parameter :: sigtol=1.e-3
      real sigma(kl), sigin(kl)

      !--------------------------------------------------------------
      ! MJT read
      if (myid==0) then
        write(6,*) "Reading ",trim(o3file)
        ncstatus=nf_open(o3file,nf_nowrite,ncid)
        if (ncstatus.eq.0) then
          write(6,*) "Ozone in NetCDF format (CMIP5)"
          ncstatus=nf_inq_dimid(ncid,'lon',valident)
          call ncmsg('lon',ncstatus)
          ncstatus=nf_inq_dimlen(ncid,valident,ii)
          call ncmsg('lon',ncstatus)
          allocate(o3lon(ii))
          ncstatus=nf_inq_varid(ncid,'lon',valident)
          call ncmsg('lon',ncstatus)
          ncstatus=nf_get_vara_real(ncid,valident,1,ii,o3lon)
          call ncmsg('lon',ncstatus)
          ncstatus=nf_inq_dimid(ncid,'lat',valident)
          call ncmsg('lat',ncstatus)
          ncstatus=nf_inq_dimlen(ncid,valident,jj)
          call ncmsg('lat',ncstatus)
          allocate(o3lat(jj))
          ncstatus=nf_inq_varid(ncid,'lat',valident)
          call ncmsg('lat',ncstatus)
          ncstatus=nf_get_vara_real(ncid,valident,1,jj,o3lat)
          call ncmsg('lat',ncstatus)
          ncstatus=nf_inq_dimid(ncid,'plev',valident)
          call ncmsg('plev',ncstatus)
          ncstatus=nf_inq_dimlen(ncid,valident,kk)
          call ncmsg('plev',ncstatus)
          allocate(o3pres(kk))
          ncstatus=nf_inq_varid(ncid,'plev',valident)
          call ncmsg('plev',ncstatus)
          ncstatus=nf_get_vara_real(ncid,valident,1,kk,o3pres)
          call ncmsg('plev',ncstatus)
          ncstatus=nf_inq_dimid(ncid,'time',valident)
          call ncmsg('time',ncstatus)
          ncstatus=nf_inq_dimlen(ncid,valident,tt)
          call ncmsg('time',ncstatus)
          ncstatus=nf_inq_varid(ncid,'time',valident)
          call ncmsg('time',ncstatus)
          ncstatus=nf_get_att_text(ncid,valident,'units',cdate)
          call ncmsg('time',ncstatus)
          ncstatus=nf_get_vara_int(ncid,valident,1,1,iti)
          call ncmsg('time',ncstatus)
          write(6,*) "Found ozone dimensions ",ii,jj,kk,tt
          allocate(o3dum(ii,jj,kk,3))
          read(cdate(14:17),*) yy
          read(cdate(19:20),*) mm
          yy=yy+iti/12
          mm=mm+mod(iti,12)
          write(6,*) "Requested date ",jyear,jmonth
          write(6,*) "Initial ozone date ",yy,mm
          nn=(jyear-yy)*12+(jmonth-mm)+1
          if (nn.lt.1.or.nn.gt.tt) then
            write(6,*) "ERROR: Cannot find date in ozone data"
            stop
          end if
          write(6,*) "Found ozone data at index ",nn
          spos=1
          npos(1)=ii
          npos(2)=jj
          npos(3)=kk
          npos(4)=1
          write(6,*) "Reading O3"
          ncstatus=nf_inq_varid(ncid,'O3',valident)
          call ncmsg('O3',ncstatus)
          spos(4)=max(nn-1,1)
          ncstatus=nf_get_vara_real(ncid,valident,spos,npos,
     &                              o3dum(:,:,:,1))
          call ncmsg('prev',ncstatus)
          spos(4)=nn
          ncstatus=nf_get_vara_real(ncid,valident,spos,npos,
     &                              o3dum(:,:,:,2))
          call ncmsg('curr',ncstatus)
          spos(4)=min(nn+1,tt)
          ncstatus=nf_get_vara_real(ncid,valident,spos,npos,
     &                              o3dum(:,:,:,3))
          call ncmsg('next',ncstatus)
          ncstatus=nf_close(ncid)
          call ncmsg('ozone file',ncstatus)
        else
          write(6,*) "Ozone in ASCII format (CMIP3)"
          ii=0
          jj=0
          kk=0
          open(16,file=o3file,form='formatted',status='old')
          read(16,*) nlev
          if ( nlev.ne.kl ) then
            write(6,*) ' ERROR - Number of levels wrong in o3_data file'
	      stop
          end if
c         Check that the sigma levels are the same
c         Note that the radiation data has the levels in the reverse order
          read(16,*) (sigin(i),i=kl,1,-1)
          do k=1,kl
	      if ( abs(sigma(k)-sigin(k)) .gt. sigtol ) then
	        write(6,*) ' ERROR - sigma level wrong in o3_data file'
	        write(6,*) k, sigma(k), sigin(k)
	        stop
            end if
          end do
          
          allocate(dduo3n(37,kl),ddo3n2(37,kl))
          allocate(ddo3n3(37,kl),ddo3n4(37,kl))
c         Note that the data is written as MAM, SON, DJF, JJA. The arrays in
c         o3dat are in the order DJF, MAM, JJA, SON
          read(16,1000) ddo3n2
          read(16,1000) ddo3n4
          read(16,1000) dduo3n
          read(16,1000) ddo3n3
          close(16)
 1000     format(9f8.5)
        end if
        write(6,*) "Finished reading ozone data"
      end if
      call MPI_Bcast(ii,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if (ii.gt.0) then
        call MPI_Bcast(jj,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(kk,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        allocate(o3pre(ifull,kk),o3mth(ifull,kk),o3nxt(ifull,kk))
        if (myid.ne.0) then
          allocate(o3lon(ii),o3lat(jj),o3pres(kk))
          allocate(o3dum(ii,jj,kk,3))
        end if
        call MPI_Bcast(o3dum,ii*jj*kk*3,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(o3lon,ii,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(o3lat,jj,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(o3pres,kk,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        if (myid==0) write(6,*) "Interpolate ozone data to CC grid"
        call o3regrid(o3pre,o3mth,o3nxt,o3dum,o3lon,o3lat,
     &                ii,jj,kk)
        deallocate(o3dum,o3lat,o3lon)
      else
        if (myid.ne.0) then
          allocate(dduo3n(37,kl),ddo3n2(37,kl))
          allocate(ddo3n3(37,kl),ddo3n4(37,kl))
        end if
        call MPI_Bcast(ddo3n2,37*kl,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(ddo3n4,37*kl,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(dduo3n,37*kl,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(ddo3n3,37*kl,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        call resetd(dduo3n,ddo3n2,ddo3n3,ddo3n4,37*kl)
      end if
      
      if (myid==0) write(6,*) "Finished processing ozone data"
      !--------------------------------------------------------------
      
      return
      end subroutine o3_read

! Version of o3set for global CC model. Based on GCM o3set.f, Revision 1.7 

! In this version the latitude may be different for each point
      subroutine o3set(alat,alon,npts,mins,duo3n,sig,ps)
c
c  This routine interpolates in latitude and time to set the ozone 
c  amounts.
c  INPUT
c    ALAT    latitude (from -pi/2 to pi/2)
c    MINS    current model time in mins (from start of year)
c  OUTPUT
c    DUO3N   ozone mixing ratio
c
      implicit none

      include 'newmpar.h'
      include 'const_phys.h'
      integer, intent(in) :: npts,mins
      integer j,ilat,m
      real, parameter :: rlag=14.8125
      real, parameter :: year=365
      real date,rang,rsin1,rcos1,rcos2,theta,angle,than
      real do3,do3p
      real, intent(in),  dimension(npts) :: ps
      real, dimension(kl), intent(in) :: sig
      real, parameter :: amd=28.9644
      real, parameter :: amo=48.
      real, parameter :: dobson=6.022e3/2.69/48.e-3
      real rlon(ii),rlat(jj)
      real duo3n(npts,kl), alat(npts),alon(npts) ! alat and alon are only needed for CMIP3 ozone

      if (allocated(o3mth)) then ! CMIP5 ozone
      
        call fieldinterpolate(duo3n,o3pre,o3mth,o3nxt,o3pres,npts,kl,
     &                        ii,jj,kk,mins,sig,ps)

        ! convert units
        where (duo3n.lt.1.)
          duo3n=duo3n*amo/amd
        end where
         
      else ! CMIP3 ozone
c       This moved to initfs
c
c       Convert time to day number
c       date = amod( float(mins)/1440., year)
c       Use year of exactly 365 days
        date = float(mod(mins,525600))/1440.
        rang = tpi*(date-rlag)/year
        rsin1 = sin(rang)
        rcos1 = cos(rang)
        rcos2 = cos(2.0*rang)
c
        do j=1,npts
           theta=90.-alat(j)*180./pi
           ilat = theta/5.
           angle = 5 * ilat
           than = (theta-angle)/5.
           ilat = ilat+1
           do m = 1,kl
            do3  = dduo3n(ilat,m) + rsin1*ddo3n2(ilat,m) 
     &              + rcos1*ddo3n3(ilat,m) + rcos2*ddo3n4(ilat,m)
            do3p = dduo3n(ilat+1,m) + rsin1*ddo3n2(ilat+1,m) 
     &              + rcos1*ddo3n3(ilat+1,m) + rcos2*ddo3n4(ilat+1,m)
            duo3n(j,m)=do3+than*(do3p-do3)
           end do
        end do
        ! convert from cm stp to gm/gm
        do m=1,kl
          duo3n(:,m)=duo3n(:,m)*1.01325e2/(ps(:)*10.)
        end do
        
      end if
      
      return
      end subroutine o3set

      subroutine fieldinterpolate(out,fpre,fmth,fnxt,fpres,ipts,ilev,
     &                            nlon,nlat,nlev,mins,sig,ps)
      
      implicit none

      include 'dates.h'
      
      integer, intent(in) :: ipts,ilev,nlon,nlat,nlev,mins
      integer date,j,ip,m,k1,leap,jyear,jmonth
      real, dimension(ipts,ilev), intent(out) :: out
      real, dimension(ipts,nlev), intent(in) :: fpre,fmth,fnxt
      real, dimension(nlev), intent(in) :: fpres
      real, dimension(ipts), intent(in) :: ps
      real, dimension(ilev), intent(in) :: sig
      real, dimension(ilev) :: prf,o3new
      real, dimension(nlev) :: o3inp,o3sum,b,c,d
      real, dimension(nlev,3) :: o3tmp
      real, dimension(12) :: monlen
      real, dimension(12), parameter :: oldlen=
     & (/ 0.,31.,59.,90.,120.,151.,181.,212.,243.,273.,304.,334. /)
      real rang,fp,fjd
      integer, parameter :: ozoneintp=1 ! ozone interpolation (0=simple, 1=integrate column)
      common/leap_yr/leap  ! 1 to allow leap years

      monlen=oldlen
      jyear =kdate/10000
      jmonth=(kdate-jyear*10000)/100
      if (leap.eq.1) then
        if (mod(jyear,4)  .eq.0) monlen(3:12)=oldlen(3:12)+1
        if (mod(jyear,100).eq.0) monlen(3:12)=oldlen(3:12)
        if (mod(jyear,400).eq.0) monlen(3:12)=oldlen(3:12)+1
      end if

!     Time interpolation factors
      date = mins/1440.
      if (jmonth.eq.12) then
        rang=(date-monlen(12))/31.
      else
        rang=(date-monlen(jmonth))/(monlen(jmonth+1)-monlen(jmonth))
      end if
      if (rang.gt.1.1) then ! use 1.1 to give 10% tolerance
        write(6,*) "WARN: fieldinterpolation is outside input range"
      end if
      
      do j=1,ipts

        if (rang.le.1.) then
          ! temporal interpolation (PWCB)
          o3tmp(:,1)=fpre(j,:)
          o3tmp(:,2)=fmth(j,:)+o3tmp(:,1)
          o3tmp(:,3)=fnxt(j,:)+o3tmp(:,2)
          b=0.5*o3tmp(:,2)
          c=4.*o3tmp(:,2)-5.*o3tmp(:,1)-o3tmp(:,3)
          d=1.5*o3tmp(:,3)+4.5*o3tmp(:,1)-4.5*o3tmp(:,2)
          o3inp=b+c*rang+d*rang*rang
        else
          ! linear interpolation when rang is out-of-range
          o3inp=max(3.-2.*rang,0.)*0.5*(fmth(j,:)+fnxt(j,:))
     &         +min(2.*rang-2.,1.)*fnxt(j,:)
        end if
         
        !-----------------------------------------------------------
        ! Simple interpolation on pressure levels
        ! vertical interpolation (from LDR - Mk3.6)
        ! Note inverted levels
        if (ozoneintp.eq.0) then
          do m=nlev-1,1,-1
            if (o3inp(m).gt.1.E34) o3inp(m)=o3inp(m+1)
          end do
          prf=0.01*ps(j)*sig
          do m=1,ilev
            if (prf(m).gt.fpres(1)) then
              out(j,ilev-m+1)=o3inp(1)
            elseif (prf(m).lt.fpres(nlev)) then
              out(j,ilev-m+1)=o3inp(nlev)
            else
              do k1=2,nlev
                if (prf(m).gt.fpres(k1)) exit
              end do
              fp=(prf(m)-fpres(k1))/(fpres(k1-1)-fpres(k1))
              out(j,ilev-m+1)=(1.-fp)*o3inp(k1)+fp*o3inp(k1-1)
            end if
          end do
        !-----------------------------------------------------------
        else
        !-----------------------------------------------------------
        ! Approximate integral of ozone column
         
          ! calculate total column of ozone
          o3sum=0.
          o3sum(nlev)=o3inp(nlev)*0.5*sum(fpres(nlev-1:nlev))
          do m=nlev-1,2,-1
            if (o3inp(m).gt.1.E34) then
              o3sum(m)=o3sum(m+1)
            else
              o3sum(m)=o3sum(m+1)+o3inp(m)*0.5*(fpres(m-1)-fpres(m+1))
            end if
          end do
          if (o3inp(1).gt.1.E34) then
            o3sum(1)=o3sum(2)
          else
            o3sum(1)=o3sum(2)+o3inp(1)*(max(fpres(1),ps(j))
     &               -0.5*sum(fpres(1:2)))
          end if
        
         ! vertical interpolation
          prf=0.01*ps(j)*sig
          o3new=0.
          do m=1,ilev
            if (prf(m).gt.fpres(1)) then
              o3new(m)=o3sum(1)
            elseif (prf(m).lt.fpres(nlev)) then
              o3new(m)=o3sum(nlev)
            else
              do k1=2,nlev
                if (prf(m).gt.fpres(k1)) exit
              end do
              fp=(prf(m)-fpres(k1))/(fpres(k1-1)-fpres(k1))
              o3new(m)=(1.-fp)*o3sum(k1)+fp*o3sum(k1-1)
            end if
          end do        
         
          ! output ozone (invert levels)
          out(j,ilev)=(o3sum(1)-o3new(2))/(max(fpres(1),ps(j))
     &                 -0.5*sum(prf(1:2)))
          do m=2,ilev-1
            out(j,ilev-m+1)=2.*(o3new(m)-o3new(m+1))
     &                        /(prf(m-1)-prf(m+1))
          end do
          out(j,1)=2.*o3new(ilev)/sum(prf(ilev-1:ilev))
          out(j,:)=max(out(j,:),0.)
        end if
        !-----------------------------------------------------------
      end do
      
      return
      end subroutine fieldinterpolate

      subroutine o3regrid(o3pre,o3mth,o3nxt,o3dum,o3lon,o3lat,
     &                    nlon,nlat,nlev)
      
      use latlong_m

      implicit none

      include 'newmpar.h'
      include 'const_phys.h'

      integer, intent(in) :: nlon,nlat,nlev
      real, dimension(ifull,nlev), intent(out) :: o3pre,o3mth,o3nxt
      real, dimension(nlon,nlat,nlev,3), intent(in) :: o3dum
      real, dimension(nlon), intent(in) :: o3lon
      real, dimension(nlat), intent(in) :: o3lat
      real, dimension(nlev) :: o3tmp,b,c,d
      real, dimension(ifull) :: blon,blat
      real alonx,lonadj,serlon,serlat
      integer j,l,ilon,ilat,ip

      blon=rlongg*180./pi
      where (blon.lt.0.)
        blon=blon+360.
      end where
      blat=rlatt*180./pi

      do j=1,ifull
        
        alonx=blon(j)
        if (alonx.lt.o3lon(1)) then
          alonx=alonx+360.
          ilon=nlon
        else
          do ilon=1,nlon-1
            if (o3lon(ilon+1).gt.alonx) exit
          end do
        end if
        ip=ilon+1
        lonadj=0.
        if (ip.gt.nlon) then
          ip=1
          lonadj=360.
        end if
        serlon=(alonx-o3lon(ilon))/(o3lon(ip)+lonadj-o3lon(ilon))

        if (blat(j).lt.o3lat(1)) then
          ilat=1
          serlat=0.
        else if (blat(j).gt.o3lat(nlat)) then
          ilat=nlat-1
          serlat=1.
        else
          do ilat=1,nlat-1
            if (o3lat(ilat+1).gt.blat(j)) exit
          end do
          serlat=(blat(j)-o3lat(ilat))/(o3lat(ilat+1)-o3lat(ilat))  
        end if

        ! spatial interpolation
        do l=1,3
          d=o3dum(ip,ilat+1,:,l)-o3dum(ilon,ilat+1,:,l)
     &     -o3dum(ip,ilat,:,l)+o3dum(ilon,ilat,:,l)
          b=o3dum(ip,ilat,:,l)-o3dum(ilon,ilat,:,l)
          c=o3dum(ilon,ilat+1,:,l)-o3dum(ilon,ilat,:,l)
          o3tmp(:)=b*serlon+c*serlat+d*serlon*serlat
     &              +o3dum(ilon,ilat,:,l)
          where (o3dum(ilon,ilat,:,l).gt.1.E34
     &      .or.o3dum(ip,ilat,:,l).gt.1.E34
     &      .or.o3dum(ilon,ilat+1,:,l).gt.1.E34
     &      .or.o3dum(ip,ilat+1,:,l).gt.1.E34)
           o3tmp(:)=1.E35
          end where
             
          select case(l)
            case(1)
              o3pre(j,:)=o3tmp(:)
            case(2)
              o3mth(j,:)=o3tmp(:)
            case(3)
              o3nxt(j,:)=o3tmp(:)
          end select
        
        end do

        ! avoid interpolating missing values
        where (o3pre(j,:).gt.1.E34.or.o3mth(j,:).gt.1.E34
     &     .or.o3nxt(j,:).gt.1.E34)
          o3pre(j,:)=1.E35
          o3mth(j,:)=1.E35
          o3nxt(j,:)=1.E35
        end where

      end do

      return
      end subroutine o3regrid

      end module ozoneread