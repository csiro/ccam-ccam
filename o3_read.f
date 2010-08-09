      !--------------------------------------------------------------
      ! MJT - radiation
      ! Convert this subroutine to a module, so as to pass allocatable
      ! arrays to o3set.  These arrays are used for the CMIP5 Ozone
      ! datasets.
      module ozoneread

      implicit none

      private      
      public o3_read,o3set
      
      integer, save :: ii,jj,kk
      real, dimension(:,:,:), allocatable, save :: o3pre,o3mth,o3nxt
      real, dimension(:), allocatable, save :: pres
      
      contains
      !--------------------------------------------------------------
            
      subroutine o3_read(sigma,jyear,jmonth)
c
c     Reads the ozone data from the o3_datafile (filename set in namelist)
c
      use cc_mpi, only : myid ! MJT read

      include 'newmpar.h'
      include 'filnames.h'   ! MJT radiation
      include 'netcdf.inc'   ! MJT radiation
      include 'mpif.h'       ! MJT read
      integer, parameter :: l=kl
      integer, intent(in) :: jyear,jmonth
      integer nlev,i,k,ierr
      integer ncstatus,ncid,tt
      integer valident,yy,mm,iti,nn
      integer, dimension(4) :: spos,npos
      character*32 cdate

      real, parameter :: sigtol=1.e-3
      real dduo3n,ddo3n2,ddo3n3,ddo3n4
      real rv
      common /o3dat/ dduo3n(37,l),ddo3n2(37,l),ddo3n3(37,l),ddo3n4(37,l)
      real sigma(kl), sigin(kl)

      !--------------------------------------------------------------
      ! MJT read
      if (myid==0) then
        ncstatus=nf_open(o3file,nf_nowrite,ncid)
        if (ncstatus.eq.0) then
          write(6,*) "Ozone in NetCDF format (CMIP5)"
          dduo3n=-1.
          ddo3n2=-1.
          ddo3n3=-1.
          ddo3n4=-1.
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
        call MPI_Bcast(ddo3n2,37*l,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(ddo3n4,37*l,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(dduo3n,37*l,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(ddo3n3,37*l,MPI_REAL,0,MPI_COMM_WORLD,ierr)
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
      include 'newmpar.h'
      include 'const_phys.h'
      integer, intent(in) :: npts,mins
      integer j,ilat,m
      integer ilon,ip,k1
      real, parameter :: rlag=14.8125
      real, parameter :: year=365
      real date,rang,rsin1,rcos1,rcos2,theta,angle,than
      real do3,do3p
      real dlat,dlon,serlat,serlon,prh,presh,fp
      real, intent(in),  dimension(npts) :: ps
      real, dimension(kl), intent(in) :: sig
      real, dimension(kk) :: o3inp,b,c,d,o3sum
      real, dimension(kk,3) :: o3tmp
      real, dimension(kl) :: prf,o3new
      real, parameter :: amd=28.9644
      real, parameter :: amo=48.
      real, parameter :: dobson=6.022e3/2.69/48.e-3

      real duo3n(npts,kl), sigma(kl), alat(npts),alon(npts)
      real dduo3n,ddo3n2,ddo3n3,ddo3n4
c          winter       spring       summer       autumn       (nh)
      common /o3dat/ dduo3n(37,kl),ddo3n2(37,kl),ddo3n3(37,kl),
     &     ddo3n4(37,kl)
!      logical start
!      data start / .true. /
!      save start
      real, parameter, dimension(12) :: monlen =  
     &  (/ 0.,31.,59.,90.,120.,151.,181.,212.,243.,273.,304.,334. /)

      if (allocated(o3mth)) then
!       Time interpolation factors (assume year of exactly 365 days)
        date = real(modulo(mins,525600))/1440.
        if (date>=monlen(12)) then
          rang=(date-monlen(12))/31.
        else
!         Search for bracketing dates, i such that monmid(i) >= date
          do j=2,12
            if (monlen(j)>=date) exit
          end do
          rang=(date-monlen(j-1))/(monlen(j)-monlen(j-1))
        end if
        rang=min(max(rang,0.),1.)
        
        dlat=180./real(jj-1)
        dlon=360./real(ii-1)
        do j=1,npts
         serlat=(90.+alat(j)*180./pi)/dlat+1
         serlon=alon(j)*180./(pi*dlon)+1
         ilat=int(serlat)
         ilon=int(serlon)
         serlat=serlat-real(ilat)
         serlon=serlon-real(ilon)
         ip=ilon+1
         if (ip.gt.ii) ip=1

         ! spatial interpolation (previous month)
         d=o3pre(ip,ilat+1,:)-o3pre(ilon,ilat+1,:)
     &    -o3pre(ip,ilat,:)+o3pre(ilon,ilat,:)
         b=o3pre(ip,ilat,:)-o3pre(ilon,ilat,:)
         c=o3pre(ilon,ilat+1,:)-o3pre(ilon,ilat,:)
         o3tmp(:,1)=b*serlon+c*serlat+d*serlon*serlat+o3pre(ilon,ilat,:)
         where (o3pre(ilon,ilat,:).gt.1..or.o3pre(ip,ilat,:).gt.1.
     &     .or.o3pre(ilon,ilat+1,:).gt.1..or.o3pre(ip,ilat+1,:).gt.1.)
           o3tmp(:,1)=1.E35
         end where
         
         ! spatial interpolation (current month)
         d=o3mth(ip,ilat+1,:)-o3mth(ilon,ilat+1,:)
     &    -o3mth(ip,ilat,:)+o3mth(ilon,ilat,:)
         b=o3mth(ip,ilat,:)-o3mth(ilon,ilat,:)
         c=o3mth(ilon,ilat+1,:)-o3mth(ilon,ilat,:)
         o3tmp(:,2)=b*serlon+c*serlat+d*serlon*serlat+o3mth(ilon,ilat,:)
         where (o3mth(ilon,ilat,:).gt.1..or.o3mth(ip,ilat,:).gt.1.
     &     .or.o3mth(ilon,ilat+1,:).gt.1..or.o3mth(ip,ilat+1,:).gt.1.)
           o3tmp(:,2)=1.E35
         end where         

         ! spatial interpolation (next month)
         d=o3nxt(ip,ilat+1,:)-o3nxt(ilon,ilat+1,:)
     &    -o3nxt(ip,ilat,:)+o3nxt(ilon,ilat,:)
         b=o3nxt(ip,ilat,:)-o3nxt(ilon,ilat,:)
         c=o3nxt(ilon,ilat+1,:)-o3nxt(ilon,ilat,:)
         o3tmp(:,3)=b*serlon+c*serlat+d*serlon*serlat+o3nxt(ilon,ilat,:)
         where (o3nxt(ilon,ilat,:).gt.1..or.o3nxt(ip,ilat,:).gt.1.
     &     .or.o3nxt(ilon,ilat+1,:).gt.1..or.o3nxt(ip,ilat+1,:).gt.1.)
           o3tmp(:,3)=1.E35
         end where
         
         ! avoid interpolating missing values
         where (o3tmp(:,1).gt.1..or.o3tmp(:,2).gt.1.
     &      .or.o3tmp(:,3).gt.1.)
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
         
         ! convert units
         where (o3inp.lt.1.)
           o3inp=o3inp*dobson*amo/(amd*grav)
         end where
         
         !-----------------------------------------------------------
         ! Simple interpolation on pressure levels
!         ! vertical interpolation (from LDR - Mk3.6)
!         do m=kl-1,1,-1
!           if (o3inp(m).gt.1) o3inp(m)=o3inp(m+1)
!         end do
!         prf=0.01*ps(j)*sig
!         do m=1,kl
!           if (prf(m).gt.pres(1)) then
!             duo3n(j,kl-m+1)=o3inp(1)
!           elseif (prf(m).lt.pres(kk)) then
!             duo3n(j,kl-m+1)=o3inp(kk)
!           else
!             do k1=2,kk
!               if (prf(m).gt.pres(k1)) exit
!             end do
!             fp=(prf(m)-pres(k1))/(pres(k1-1)-pres(k1))
!             duo3n(j,kl-m+1)=(1.-fp)*o3inp(k1)+fp*o3inp(k1-1)
!           end if
!         end do
         !-----------------------------------------------------------
         
         !-----------------------------------------------------------
         ! Approximate integral of ozone column
         
         ! calculate total column of ozone
         o3sum=0.
         o3sum(kk)=o3inp(kk)*0.5*sum(pres(kk-1:kk))
         do m=kk-1,2,-1
           if (o3inp(m).gt.1.) then
             o3sum(m)=o3sum(m+1)
           else
             o3sum(m)=o3sum(m+1)+o3inp(m)*0.5*(pres(m-1)-pres(m+1))
           end if
         end do
         if (o3inp(1).gt.1) then
           o3sum(1)=o3sum(2)
         else
           o3sum(1)=o3sum(2)+o3inp(1)*0.5*(1000.-0.5*sum(pres(1:2)))
         end if
         
        ! vertical interpolation (from LDR - Mk3.6)
         prf=0.01*ps(j)*sig
         o3new=0. 
         do m=1,kl-1
           prh=0.5*sum(prf(m:m+1))
           if (prh.gt.0.5*sum(pres(1:2))) then
             o3new(m)=o3sum(1)
           elseif (prh.lt.0.5*sum(pres(kk-1:kk))) then
             o3new(m)=o3sum(kk) ! =zero
           else
             do k1=2,kk-1
               presh=0.5*sum(pres(k1:k1+1))
               if (prh.gt.presh) exit
             end do
             fp=2.*(prh-presh)/(pres(k1-1)-pres(k1+1))
             o3new(m)=(1.-fp)*o3sum(k1)+fp*o3sum(k1-1)
           end if
         end do        
         
         ! output ozone
         duo3n(j,kl)=(o3sum(1)-o3new(1))/(ps(j)-0.5*sum(prf(1:2)))
         do m=2,kl-1
           duo3n(j,kl-m+1)=2.*(o3new(m-1)-o3new(m))/(prf(m-1)-prf(m+1))
         end do
         duo3n(j,1)=2.*o3new(kl-1)/sum(prf(kl-1:kl))
         duo3n(j,:)=max(duo3n(j,:),0.)
         !-----------------------------------------------------------
         
         !-----------------------------------------------------------
         ! Check ozone column
!         o3s1=0.
!         o3s1=o3s1+duo3n(j,kl)*(ps(j)-0.5*sum(prf(1:2)))
!         do m=2,kl-1
!           o3s1=o3s1+duo3n(j,kl-m+1)*0.5*(prf(m-1)-prf(m+1))
!         end do
!         o3s1=o3s1+duo3n(j,1)*0.5*sum(prf(kl-1:kl))
!         o3s2=0.
!         if (o3inp(1).lt.1.) then
!           o3s2=o3s2+o3inp(1)*(1000.-0.5*sum(pres(1:2)))
!         end if
!         do m=2,kk-1
!           if (o3inp(m).lt.1.) then
!             o3s2=o3s2+o3inp(m)*0.5*(pres(m-1)-pres(m+1))
!           end if
!         end do
!         o3s2=o3s2+o3inp(kk)*0.5*sum(pres(kk-1:kk))
!         print *,"o3s1,o3s2 ",o3s1,o3s2
         !-----------------------------------------------------------
         
        end do        
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
      end if
      
      return
      end subroutine o3set

      end module ozoneread