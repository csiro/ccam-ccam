      program amipdata  ! for AMIP 2
c     Program to extract AMIP ice-fraction and SST from 1 deg data set
c     for conformal-cubic model;  jlm Tue  10-08-1996
c     usually run on cherax ~csjlm/globpe
c     There is no interpolation in time. The AMIP data are on a 1x1 degree
c     (360 x 180) grid.       (0.5 to 359.5, -89.5 to 89.5)
c     There is no land mask in the AMIP data which also makes the 
c     interpolation easier.
c     Note poles are not included in data, so N-S extension should be done.
      include 'newmpar.h'
      include 'dates.h'   ! to pass ds from setxyz
      include 'parm.h'
      include 'latlong.h'  ! rlatt,rlongg
      parameter (longs=360,lats=180)
      common/work/zs(ifull),a(ifull),b(ifull),iarr(ifull),
     .         globidat(longs,lats),globdat(longs,lats),idat(longs,lats)
      character*80 topofile,amipsst,amipice
     .            ,iceout,sstout
      character header*47
      character*9 formout1,formout2
      data pi/3.1415926536/
      namelist/amipnml/amipsst                    ! global sst data files
     .                ,amipice                    ! global ice data files
     .                ,topofile                   ! model topography file
     .                ,iceout,sstout              ! model output data files
     .                ,months                     ! no. of months to do

      open (85,file='amipdata.nml',status='old')
      read (85,amipnml)
      write (6,amipnml)
      open (unit=51,file=topofile,form='formatted',status='old')

!       read header of topo file and check its dimensions
        read(51,'(i3,i4,2f6.1,f6.3,f8.0,a47)')
     .            ilx,jlx,rlong0,rlat0,schmidt,dsx,header
        print *,ilx,jlx,rlong0,rlat0,schmidt,dsx,header
        if(ilx.ne.il.or.jlx.ne.jl)
     .              stop 'wrong topo file supplied'
        call setxyz
c       convert lat/longs to degrees in this routine
        do iq=1,ifull
         rlatt(iq)=rlatt(iq)*180./pi
         rlongg(iq)=rlongg(iq)*180./pi
         if(rlongg(iq).lt.0.)rlongg(iq)=rlongg(iq)+360. ! here in degrees
        enddo
        lat1=1       ! lat1, lat2, long1, long2 for diag prints
        lat2=91
        long1=65
        long2=115

      read (51,*) zs
c     read (51,*) zs   ! second read gives real land-sea mask 1. to 0.
c     work through i,j, values of the grid
      do j=1,jl
       do i=1,il
        iq=i+(j-1)*il
c       set lmask (here 1 for land, 0 for sea), just for diag printing
        iarr(i,j)=0
        if(zs(i,j).gt.0.)iarr(i,j)=1
       enddo  ! i loop
      enddo   ! j loop
c     print out land-sea mask data for model grid by panel number
      do n=npanels,0,-1
       ja=n*il +1
       jb=(n+1)*il
       print *,'land-sea mask for panel ',n
       do j=jb,ja,-1
        write(6,'(i4,1x,180i1)') j,(iarr(i,j),i=1,il)
       enddo  !  j loop
      enddo   !  n loop

      write (formout1,'(1h(,i3,2hi4,1h))')il    !  i.e. (<il>i4)
!     write (formout2,'(1h(,i3,4hf6.2,1h))')il  !  i.e. (<il>f6.2)
      write (formout2,'(1h(,i3,2hi5,1h))')il    !  i.e. (<il>i5)
      open (unit=10,file=amipsst,status='old')
      open (unit=15,file=amipice,status='old')
      open (unit=20,file=sstout,status='unknown')
      open (unit=30,file=iceout,status='unknown')

      do imonth = 1,months
        read(10,*) iyear,month
        print *,'AMIP SST data read for iyear,month ',iyear,month
        read(10,'(16f6.2)') globdat    !  SST in Celsius
c       read(10,*) globdat             !  SST in Celsius
        call sst(lat1,lat2,long1,long2,ja,jb,imonth)
        write(20,'(i2,i5,2i4,2f6.1,f6.3,f8.0,
     .           '' AMIP SSTs (Celsius+50)'')')
     .       month,iyear,il,jl,rlong0,rlat0,schmidt,ds
!       write(20,formout2) a     ! write sst -  in Celsius 
        do iq=1,ifull
         iarr(iq)=nint(100.*(a(iq)+50.)) ! tenths of degree with 50 offset
        enddo
        write(20,formout2) iarr    ! write sst - as 100*(Celsius +50)

        read(15,*) iyear,month
        print *,'AMIP ice data read for iyear,month ',iyear,month
        read(15,'(16f6.1)') globidat
c       read(15,*) globidat             !  ice-frac as %
        call icefrac(lat1,lat2,long1,long2,ja,jb,imonth)
        write(30,'(i2,i5,2i4,2f6.1,f6.3,f8.0,'' AMIP_ice-fraction'')')
     .       month,iyear,il,jl,rlong0,rlat0,schmidt,ds
        write(30,formout1) iarr     ! ice-frac ( 1-100 for ice)
      enddo    ! imonth loop
      end

      subroutine sst(lat1,lat2,long1,long2,ja,jb,imonth)
      include 'newmpar.h'
      include 'latlong.h'  ! rlatt,rlongg
      parameter (longs=360,lats=180)
      common/work/zs(il,jl),a(il,jl),b(il,jl),iarr(il,jl),
     .         globidat(longs,lats),globdat(longs,lats),idat(longs,lats)

      call int2x(a,zs,il,jl,imonth,globdat,
     .          rlongg,rlatt,.5,90.5)   ! (0,-90) maps to (i,j)=(-.5,-.5)

      do iq=1,il*jl
        if(zs(iq).gt.0.) a(iq)=33.    ! preset for land points
        iarr(iq)=a(iq)/3.             ! just for diag prints
      enddo

      if(imonth.gt.1)return
c     print out part of global data
      print *,' typical global sst values:',(globdat(lat,lat),lat=1,180)
      print *,' global sst/3 in C'
      do lat=lat2,lat1,-1
        do long=long1,long2
         idat(long,lat)=globdat(long,lat)/3.
        enddo
        write(6,'(180i1)') (idat(long,lat),long=long1,long2)
      enddo
c     print out sst for model grid by panel number
      print *
      do n=npanels,0,-1
       ja=n*il +1
       jb=(n+1)*il
c      print *,'sst/3 (in C) for panel ',n
       print *,'sst (in C) for panel ',n
       do j=jb,ja,-1
c       write(6,'(i4,1x,180i1)') j,(iarr(i,j),i=1,il)
        write(6,'(i4,180f5.1)') j,(a(i,j),i=1,il)
       enddo  !  j loop
      enddo   !  n loop
      return
      end

      subroutine icefrac(lat1,lat2,long1,long2,ja,jb,imonth)
      include 'newmpar.h'
      include 'latlong.h'  ! rlatt,rlongg
      parameter (longs=360,lats=180)
      common/work/zs(il,jl),a(il,jl),b(il,jl),iarr(il,jl),
     .         globidat(longs,lats),globdat(longs,lats),idat(longs,lats)

      call int2x(b,zs,il,jl,imonth,globidat,
     .          rlongg,rlatt,.5,90.5)   ! (0,-90) maps to (i,j)=(-.5,-.5)

      do iq=1,ifull
!        b(iq)=max(0.,min(b(iq),100.)) ! between 0 and 1 no!!!!
        if(zs(iq).gt.0.) b(iq)=-9.    ! preset for land points
        iarr(iq)=nint(b(iq))          ! just for final output of sice
      enddo

      if(imonth.gt.1)return
c     print out part of global data
      print *,' some global ice values:',(globidat(lat,lat),lat=1,180)
      print *,' global ice_frac/15'
      do lat=lat2,lat1,-1
        do long=long1,long2
         idat(long,lat)=globidat(long,lat)/15.
        enddo
        write(6,'(180i1)') (idat(long,lat),long=long1,long2)
      enddo
c     print out icefrac for model grid by panel number
      print *
      do n=npanels,0,-1
       ja=n*il +1
       jb=(n+1)*il
      print *,'ice-frac/15 for panel ',n
       do j=jb,ja,-1
        write(6,'(i4,180f5.1)') j,(b(i,j)/15.,i=1,il)
       enddo  !  j loop
      enddo   !  n loop
      return
      end

      subroutine int2x(a,zs,il,jl,imonth,glob,
     .                    rlong,rlat,addlong,addlat)
c     This one for 2 degree spacing and lats from -90 to 90
c     This routine:
c       b) does quadratic interpolation to find values for a(i,j)
c          wherever zs(i,j)>0  (a not touched where zs<0)
c       c) N.B. interpolation is at (long,lat) positions given by the
c         arrays rlong(il*jl) and rlat(il*jl) with displacements
c         addlong and addlat
      parameter (longs=360,lats=180)
      real glob(longs*lats)      ! input and output global array
      dimension ie(longs*lats),ine(longs*lats)
      dimension ise(longs*lats),in(longs*lats)
      dimension iw(longs*lats),inw(longs*lats)
      dimension isw(longs*lats),is(longs*lats)
      real a(il,jl),zs(il,jl),rlong(il*jl),rlat(il*jl)
      data num/0/
      save num,ie,ine,ise,in,iw,inw,isw,is
      if(num.eq.0)then
        num=1
        do lat=1,lats
         do long=1,longs
          ll=long+(lat-1)*longs
          ie(ll)=ll+1
          in(ll)=ll+longs
          iw(ll)=ll-1
          is(ll)=ll-longs
         enddo
c        apply ew periodicity
         long=1
         ll=long+(lat-1)*longs
         iw(ll)=ll-1     +longs
         long=longs
         ll=long+(lat-1)*longs
         ie(ll)=ll+1     -longs
        enddo  !  lat loop
        do long=1,longs/2
         lat=1
         ll=long+(lat-1)*longs
         is(ll)=ll +longs/2
         is(ll+longs/2)=ll
         lat=lats
         ll=long+(lat-1)*longs
         in(ll)=ll +longs/2
         in(ll+longs/2)=ll
        enddo  !  long loop
        do ll=1,longs*lats
         ine(ll)=in(ie(ll))
         inw(ll)=in(iw(ll))
         ise(ll)=is(ie(ll))
         isw(ll)=is(iw(ll))
         ine(ll)=in(ie(ll))
        enddo
      endif
c     --------------------------------------------------------------
c     now do conditional biperiodic interpolation

c     work through i,j, values of the grid
      do j=1,jl
       do i=1,il
        if(zs(i,j).le.0.)then
          iq=i+(j-1)*il
c         centre around nearest lat/long point
          rlongx=rlong(iq)+addlong    ! 1 degree global grid
          rlatx=rlat(iq)+addlat
          long=nint(rlongx)
          lat=nint(rlatx)
c         lat=max(lat,2)        ! for south pole: 2 deg grid
c         lat=min(lat,lats-1)   ! for north pole: 2 deg grid
          p=rlongx-long   ! for unit spacing of long
          q=rlatx-lat     ! for unit spacing of lat
          ll=long+(lat-1)*longs
          aw=.5*q*(q-1.)*glob(isw(ll)) +(1.-q*q)*glob(iw(ll))
     .      +.5*q*(q+1.)*glob(inw(ll))
          a0=.5*q*(q-1.)*glob(is(ll)) +(1.-q*q)*glob(ll)
     .      +.5*q*(q+1.)*glob(in(ll))
          ae=.5*q*(q-1.)*glob(ise(ll)) +(1.-q*q)*glob(ie(ll))
     .      +.5*q*(q+1.)*glob(ine(ll))
          a(i,j)=.5*p*(p-1.)*aw +(1.-p*p)*a0 +.5*p*(p+1.)*ae
          if(imonth.eq.1.and.
     .           i.ge.17.and.i.le.21.and.j.ge.17.and.j.le.21)then
            print *
            print *,'i,j,long,lat,p,q ',i,j,long,lat,p,q
            print *,'inw(ll),in(ll),ine(ll) ',inw(ll),in(ll),ine(ll)
            print *,' iw(ll),   ll,  ie(ll) ',iw(ll),ll,ie(ll)
            print *,'isw(ll),is(ll),ise(ll) ',isw(ll),is(ll),ise(ll)
            print *,'gnw,gn,gne '
     .         ,glob(inw(ll)),glob(in(ll)),glob(ine(ll))
            print *,'gw ,g ,ge  '
     .         ,glob(iw(ll)),glob(ll),glob(ie(ll))
            print *,'gsw,gs,gse '
     .         ,glob(isw(ll)),glob(is(ll)),glob(ise(ll))
            print *,'aw,a0,ae,a ',aw,a0,ae,a(i,j)
          endif
        endif   ! zs(i,j).le.0.
       enddo  ! i loop
      enddo   ! j loop
      return
      end

c      include 'setxyz.f'      ! for conformal-cubic
c      include 'jimcc.f'
