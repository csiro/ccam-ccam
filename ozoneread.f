      !--------------------------------------------------------------
      ! MJT - radiation
      ! Convert this subroutine to a module, so as to pass allocatable
      ! arrays to o3set.  These arrays are used for the CMIP5 Ozone
      ! datasets.
      module ozoneread

      implicit none

      private      
      public o3_read,o3set,fieldinterpolate
      
      integer, save :: ii,jj,kk
      real, dimension(:,:,:), allocatable, save :: o3pre,o3mth,o3nxt
      real, dimension(:), allocatable, save :: pres
      real, dimension(:,:), allocatable, save :: dduo3n,ddo3n2
      real, dimension(:,:), allocatable, save :: ddo3n3,ddo3n4
      
      contains
      !--------------------------------------------------------------
            
      subroutine o3_read(sigma,jyear,jmonth)
c
c     Reads the ozone data from the o3_datafile (filename set in namelist)
c
      use cc_mpi, only : myid ! MJT read

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
      character*32 cdate

      real, parameter :: sigtol=1.e-3
      real rv
      real sigma(kl), sigin(kl)

      !--------------------------------------------------------------
      ! MJT read
      if (myid==0) then
        ncstatus=nf_open(o3file,nf_nowrite,ncid)
        if (ncstatus.eq.0) then
          write(6,*) "Ozone in NetCDF format (CMIP5)"
          ncstatus=nf_inq_dimid(ncid,'lon',valident)
          ncstatus=nf_inq_dimlen(ncid,valident,ii)
          ncstatus=nf_inq_varid(ncid,'lon',valident)
          ncstatus=nf_get_vara_real(ncid,valident,1,1,rv)
          if (rv.ne.0.) then
            write(6,*) "ERROR: Need to translate lon"
            stop
          end if          
          ncstatus=nf_inq_dimid(ncid,'lat',valident)
          ncstatus=nf_inq_dimlen(ncid,valident,jj)
          ncstatus=nf_inq_varid(ncid,'lat',valident)
          ncstatus=nf_get_vara_real(ncid,valident,1,1,rv)
          if (rv.ne.-90.) then
            write(6,*) "ERROR: Need to invert lat"
            stop
          end if
          ncstatus=nf_inq_dimid(ncid,'plev',valident)
          ncstatus=nf_inq_dimlen(ncid,valident,kk)
          allocate(pres(kk))
          ncstatus=nf_inq_varid(ncid,'plev',valident)
          ncstatus=nf_get_vara_real(ncid,valident,1,kk,pres)
          ncstatus=nf_inq_varid(ncid,'time',valident)
          ncstatus=nf_get_att_text(ncid,valident,'units',cdate)
          ncstatus=nf_get_vara_int(ncid,valident,1,1,iti)
          ncstatus=nf_inq_varid(ncid,'O3',valident)
          read(cdate(14:17),*) yy
          read(cdate(19:20),*) mm
          yy=yy+iti/12
          mm=mm+mod(iti,12)
          nn=(jyear-yy)*12+(jmonth-mm)+1
          spos=1
          npos(1)=ii
          npos(2)=jj
          npos(3)=kk
          npos(4)=1
          allocate(o3pre(ii,jj,kk),o3mth(ii,jj,kk),o3nxt(ii,jj,kk))
          spos(4)=nn-1
          ncstatus=nf_get_vara_real(ncid,valident,spos,npos,o3pre)
          spos(4)=nn
          ncstatus=nf_get_vara_real(ncid,valident,spos,npos,o3mth)
          spos(4)=nn+1
          ncstatus=nf_get_vara_real(ncid,valident,spos,npos,o3nxt)
        else
          write(6,*) "Ozone in ASCII format (CMIP3)"
          ii=0
          jj=0
          kk=0
          open(16,file=o3file,form='formatted',status='old')
          read(16,*) nlev
          if ( nlev.ne.kl ) then
            print*, ' ERROR - Number of levels wrong in o3_data file'
	      stop
          end if
c         Check that the sigma levels are the same
c         Note that the radiation data has the levels in the reverse order
          read(16,*) (sigin(i),i=kl,1,-1)
          do k=1,kl
	      if ( abs(sigma(k)-sigin(k)) .gt. sigtol ) then
	        print*, ' ERROR - sigma level wrong in o3_data file'
	        print*, k, sigma(k), sigin(k)
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
      end if
      call MPI_Bcast(ii,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if (ii.gt.0) then
        call MPI_Bcast(jj,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(kk,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        if (myid.ne.0) then
          allocate(o3pre(ii,jj,kk),o3mth(ii,jj,kk),o3nxt(ii,jj,kk))
          allocate(pres(kk))
        end if
        call MPI_Bcast(o3pre,ii*jj*kk,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(o3mth,ii*jj*kk,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(o3nxt,ii*jj*kk,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(pres,kk,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      else
        if (myid.ne.0) then
          allocate(dduo3n(37,kl),ddo3n2(37,kl))
          allocate(ddo3n3(37,kl),ddo3n4(37,kl))
        end if
        call MPI_Bcast(ddo3n2,37*kl,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(ddo3n4,37*kl,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(dduo3n,37*kl,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(ddo3n3,37*kl,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      end if
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
      real blon(npts),blat(npts)
      real duo3n(npts,kl), alat(npts),alon(npts)

      if (allocated(o3mth)) then
        do j=1,ii
          rlon(j)=real(j-1)*360./real(ii)
        end do
        do j=1,jj
          rlat(j)=-90.+real(j-1)*180./real(jj-1)
        end do
        blon=alon*180./pi
        where (blon.lt.0.)
          blon=blon+360.
        end where
        blat=alat*180./pi
      
        call fieldinterpolate(duo3n,blon,blat,o3pre,o3mth,o3nxt,rlon,
     &                        rlat,pres,npts,kl,ii,jj,kk,mins,sig,ps)
      
        ! convert units
        where (duo3n.lt.1.)
          duo3n=duo3n*amo/amd
        end where
         
      else
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
        
        call resetd(dduo3n,ddo3n2,ddo3n3,ddo3n4,37*kl)
        
      end if
      
      return
      end subroutine o3set

      subroutine fieldinterpolate(out,alon,alat,fpre,fmth,fnxt,rlon,
     &                            rlat,fpres,ipts,ilev,nlon,nlat,nlev,
     &                            mins,sig,ps)
      
      implicit none
      
      integer, intent(in) :: ipts,ilev,nlon,nlat,nlev,mins
      integer date,ilon,ilat,j,ip,m,k1
      real, dimension(ipts,ilev), intent(out) :: out
      real, dimension(ipts), intent(in) :: alon,alat
      real, dimension(nlon,nlat,nlev), intent(in) :: fpre,fmth,fnxt
      real, dimension(nlon), intent(in) :: rlon
      real, dimension(nlat), intent(in) :: rlat
      real, dimension(nlev), intent(in) :: fpres
      real, dimension(ipts), intent(in) :: ps
      real, dimension(ilev), intent(in) :: sig
      real, dimension(ilev) :: prf,o3new
      real, dimension(nlev) :: b,c,d,o3inp,o3sum
      real, dimension(nlev,3) :: o3tmp
      real rang,serlon,serlat,fp,lonadj,alonx
      real, parameter, dimension(12) :: monlen =  
     &  (/ 0.,31.,59.,90.,120.,151.,181.,212.,243.,273.,304.,334. /)
     
      integer, parameter :: ozoneintp=1 ! ozone interpolation (0=simple, 1=integrate column)     

!     Time interpolation factors (assume year of exactly 365 days)
      date = real(modulo(mins,525600))/1440.
      if (date>=monlen(12)) then
        rang=(date-monlen(12))/31.
      else
!       Search for bracketing dates, i such that monmid(i) >= date
        do j=2,12
          if (monlen(j)>=date) exit
        end do
        rang=(date-monlen(j-1))/(monlen(j)-monlen(j-1))
      end if
      rang=min(max(rang,0.),1.)
      
      do j=1,ipts

        alonx=alon(j)
        if (alonx.lt.rlon(1)) then
          alonx=alonx+360.
          ilon=nlon
        else
          do ilon=1,nlon-1
            if (rlon(ilon+1).gt.alonx) exit
          end do
        end if
        ip=ilon+1
        lonadj=0.
        if (ip.gt.nlon) then
          ip=1
          lonadj=360.
        end if
        serlon=(alonx-rlon(ilon))/(rlon(ip)+lonadj-rlon(ilon))

        if (alat(j).lt.rlat(1)) then
          ilat=1
          serlat=0.
        else if (alat(j).gt.rlat(nlat)) then
          ilat=nlat-1
          serlat=1.
        else
          do ilat=1,nlat-1
            if (rlat(ilat+1).gt.alat(j)) exit
          end do
          serlat=(alat(j)-rlat(ilat))/(rlat(ilat+1)-rlat(ilat))  
        end if

        ! spatial interpolation (previous month)
        d=fpre(ip,ilat+1,:)-fpre(ilon,ilat+1,:)
     &   -fpre(ip,ilat,:)+fpre(ilon,ilat,:)
        b=fpre(ip,ilat,:)-fpre(ilon,ilat,:)
        c=fpre(ilon,ilat+1,:)-fpre(ilon,ilat,:)
        o3tmp(:,1)=b*serlon+c*serlat+d*serlon*serlat+fpre(ilon,ilat,:)
        where (fpre(ilon,ilat,:).gt.1.E34
     &     .or.fpre(ip,ilat,:).gt.1.E34
     &     .or.fpre(ilon,ilat+1,:).gt.1.E34
     &     .or.fpre(ip,ilat+1,:).gt.1.E34)
          o3tmp(:,1)=1.E35
        end where
         
        ! spatial interpolation (current month)
        d=fmth(ip,ilat+1,:)-fmth(ilon,ilat+1,:)
     &   -fmth(ip,ilat,:)+fmth(ilon,ilat,:)
        b=fmth(ip,ilat,:)-fmth(ilon,ilat,:)
        c=fmth(ilon,ilat+1,:)-fmth(ilon,ilat,:)
        o3tmp(:,2)=b*serlon+c*serlat+d*serlon*serlat+fmth(ilon,ilat,:)
        where (fmth(ilon,ilat,:).gt.1.E34
     &     .or.fmth(ip,ilat,:).gt.1.E34
     &     .or.fmth(ilon,ilat+1,:).gt.1.E34
     &     .or.fmth(ip,ilat+1,:).gt.1.E34)
          o3tmp(:,2)=1.E35
        end where         

        ! spatial interpolation (next month)
        d=fnxt(ip,ilat+1,:)-fnxt(ilon,ilat+1,:)
     &   -fnxt(ip,ilat,:)+fnxt(ilon,ilat,:)
        b=fnxt(ip,ilat,:)-fnxt(ilon,ilat,:)
        c=fnxt(ilon,ilat+1,:)-fnxt(ilon,ilat,:)
        o3tmp(:,3)=b*serlon+c*serlat+d*serlon*serlat+fnxt(ilon,ilat,:)
        where (fnxt(ilon,ilat,:).gt.1.E34
     &     .or.fnxt(ip,ilat,:).gt.1.E34
     &     .or.fnxt(ilon,ilat+1,:).gt.1.E34
     &     .or.fnxt(ip,ilat+1,:).gt.1.E34)
          o3tmp(:,3)=1.E35
        end where
         
        ! avoid interpolating missing values
        where (o3tmp(:,1).gt.1.E34.or.o3tmp(:,2).gt.1.E34
     &     .or.o3tmp(:,3).gt.1.E34)
          o3tmp(:,1)=1.E35
          o3tmp(:,2)=1.E35
          o3tmp(:,3)=1.E35
        end where

        ! temporal interpolation (PWCB)
        o3tmp(:,2)=o3tmp(:,2)+o3tmp(:,1)
        o3tmp(:,3)=o3tmp(:,3)+o3tmp(:,2)
        b=0.5*o3tmp(:,2)
        c=4.*o3tmp(:,2)-5.*o3tmp(:,1)-o3tmp(:,3)
        d=1.5*(o3tmp(:,3)+3.*o3tmp(:,1)-3.*o3tmp(:,2))
        o3inp=b+c*rang+d*rang*rang
         
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

      end module ozoneread