      program amipdata
c     Program to extract AMIP icemask and SST from 2 deg data set
c     for either DARLAM or conformal-cubic model;  jlm Tue  10-08-1996
c     usually run on cherax   (DARLAM not yet completed)
c     There is no interpolation in time. The AMIP data is on a 2x2 degree
c     (180 x 91) grid.       (0 to 358, -90 to 90)
c     There is no land mask in the AMIP data which also makes the 
c     interpolation easier.
c     Note poles are included in data, so no N-S extension done.
c     In fact ice just done to nearest point (for DARLAM's finer resolution
c     may want to do something fancier)
c     parameter (npanels=5)   !  0 for DARLAM,  5 for conformal-cubic model
      include 'newmpar.h'
      include 'dates.h'   ! to pass ds from setxyz
      include 'parm.h'
      include 'latlong.h'  ! rlat,rlong
      common/work/zs(ifull),iarr(ifull),iglobdat(180,91)
     .            ,a(ifull),globdat(180,91)
      character*80 topofile,amipsst
     .            ,iceout,sstout
      character header*47
      character*9 formout1,formout2
      data pi/3.1415926536/
      namelist/amipnml/amipsst                     ! global input data files
     .                 ,topofile                   ! model topography file
     .                 ,iceout,sstout              ! model output data files
     .                 ,months                     ! no. of months to do

      open (85,file='amipdata.nml',status='old')
      read (85,amipnml)
      write (6,amipnml)
      open (unit=51,file=topofile,form='formatted',status='old')

      if(npanels.eq.0)then
!       read header of topo file and check its dimensions
        read(51,*)ik,jk,ds,du,tanl,rnml,stl1,stl2
        if(ik.eq.0.or.jk.eq.0)then
         print *,'no header in newtopo file'
         stop 'now require topography file header'
        else
         print *,'Header information for topofile'
         print *,'ik,jk,ds,du,tanl,rnml,stl1,stl2'
     .           ,ik,jk,ds,du,tanl,rnml,stl1,stl2
         if(ik.ne.il.or.jk.ne.jl)stop 'wrong topofile supplied'
        endif     ! (ik.eq.0.or.jk.eq.0)
        print *,'du,tanl,rnml,stl1,stl2: ',du,tanl,rnml,stl1,stl2
        call lconset(ds)
        do j=1,jl
         do i=1,il
          iq=i+(j-1)*il
          call lconll(rlong(iq),rlat(iq),real(i),real(j))
          if(rlong(iq).lt.0.)rlong(iq)=rlong(iq)+360. ! lconset in degrees
c         N.B. the data set is every 2 deg, 0 to 358
          if(rlong(iq).gt.359.)long=1  ! to give nint to closest point
         enddo  ! i loop
        enddo   ! j loop
c       find corner of lam grid as lat/longs numbering for printouts
        lat1=.5*nint(rlat(1)+92.)             ! rlat(1) at -90.
        long1=.5*nint(rlong(1)+2.)            ! rlong(1) at 0.
        lat2=.5*nint(rlat(ifull)+92.)         ! rlat(1) at -90.
        long2=.5*nint(rlong(ifull)+2.)        ! rlong(1) at 0.
        ja=1
        jb=jl
      else
!       read header of topo file and check its dimensions
        read(51,'(i3,i4,2f6.1,f6.3,f8.0,a47)')
     .            ilx,jlx,rlong0,rlat0,schmidt,dsx,header
        print *,ilx,jlx,rlong0,rlat0,schmidt,dsx,header
        if(ilx.ne.il.or.jlx.ne.jl)
     .              stop 'wrong topo file supplied'
        call setxyz
c       convert lat/longs to degrees in this routine
        do iq=1,ifull
         rlat(iq)=rlat(iq)*180./pi
         rlong(iq)=rlong(iq)*180./pi
         if(rlong(iq).lt.0.)rlong(iq)=rlong(iq)+360. ! here in degrees
        enddo
        lat1=1       ! lat1, lat2, long1, long2 for diag prints
        lat2=91
        long1=65
        long2=115
      endif

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
       if(npanels.gt.0)then
         ja=n*il +1
         jb=(n+1)*il
       endif
       print *,'land-sea mask for panel ',n
       do j=jb,ja,-1
        write(6,'(i4,1x,180i1)') j,(iarr(i,j),i=1,il)
       enddo  !  j loop
      enddo   !  n loop

      write (formout1,'(1h(,i3,2hi1,1h))')il    !  i.e. (<il>i1)
      write (formout2,'(1h(,i3,4hf6.2,1h))')il  !  i.e. (<il>f6.2)
      open (unit=10,file=amipsst,status='old')
      open (unit=20,file=sstout,status='unknown')
      open (unit=30,file=iceout,status='unknown')

      do imonth = 1,months
        read(10,'(2i20)') month,iyear
        print *,'AMIP data read for month,iyear ',month,iyear
        read(10,'(16f5.2)') globdat    !  SST in Celsius
        call sst(lat1,lat2,long1,long2,ja,jb,imonth)
        write(20,'(i2,i5,2i4,2f6.1,f6.3,f8.0,'' AMIP SSTs (Celsius)'')')
     .       month,iyear,il,jl,rlong0,rlat0,schmidt,ds
        write(20,formout2) a     ! write sst - still in Celsius

        read(10,'(80i1)')   iglobdat
        call icemask(lat1,lat2,long1,long2,ja,jb,imonth)
        write(30,'(i2,i5,2i4,2f6.1,f6.3,f8.0,'' AMIP ice-mask 0/1'')')
     .       month,iyear,il,jl,rlong0,rlat0,schmidt,ds
        write(30,formout1) iarr     ! ice-mask (1 for ice)
      enddo    ! imonth loop
      end

      subroutine icemask(lat1,lat2,long1,long2,ja,jb,imonth)
      include 'newmpar.h'
      include 'latlong.h'  ! rlat,rlong
      common/work/zs(il,jl),iarr(il,jl),iglobdat(180,91)
     .            ,a(il,jl),globdat(180,91)
c     work through i,j, values of the grid
      do j=1,jl
       do i=1,il
        iq=i+(j-1)*il
c       N.B. the data set is every 2 deg, 0 to 358
        lat=.5*nint(rlat(iq)+92.)             ! rlat(1) at -90.
        long=.5*nint(rlong(iq)+2.)            ! rlong(1) at 0.
        iarr(i,j)=iglobdat(long,lat)
       enddo  ! i loop
      enddo   ! j loop

      if(imonth.gt.1)return
c     print out part of global data
      print *,' global icemask data'
      do lat=lat2,lat1,-1
        write(6,'(180i1)') (iglobdat(long,lat),long=long1,long2)
      enddo
c     print out icemask for model grid by panel number
      print *
      do n=npanels,0,-1
       if(npanels.gt.0)then
         ja=n*il +1
         jb=(n+1)*il
       endif
       print *,' icemask for panel ',n
       do j=jb,ja,-1
        write(6,'(i4,1x,180i1)') j,(iarr(i,j),i=1,il)
       enddo  !  j loop
      enddo   !  n loop
      return
      end

      subroutine sst(lat1,lat2,long1,long2,ja,jb,imonth)
      include 'newmpar.h'
      include 'latlong.h'  ! rlat,rlong
      common/work/zs(il,jl),iarr(il,jl),iglobdat(180,91)
     .            ,a(il,jl),globdat(180,91)

      call int2x(a,zs,il,jl,globdat,
     .          rlong,rlat,2.,92.)   ! rlong(1) at 0.; rlat(1) at -90.

      do iq=1,il*jl
        if(zs(iq).gt.0.) a(iq)=33.    ! preset for land points
        iarr(iq)=a(iq)/3.               ! just for diag prints
      enddo

      if(imonth.gt.1)return
c     print out part of global data
      print *,' typical global sst values:',(globdat(lat,lat),lat=1,180)
      print *,' global sst/3 in C'
      do lat=lat2,lat1,-1
        do long=long1,long2
         iglobdat(long,lat)=globdat(long,lat)/3.
        enddo
        write(6,'(180i1)') (iglobdat(long,lat),long=long1,long2)
      enddo
c     print out sst for model grid by panel number
      print *
      do n=npanels,0,-1
       if(npanels.gt.0)then
         ja=n*il +1
         jb=(n+1)*il
       endif
c      print *,'sst/3 (in C) for panel ',n
       print *,'sst (in C) for panel ',n
       do j=jb,ja,-1
c       write(6,'(i4,1x,180i1)') j,(iarr(i,j),i=1,il)
        write(6,'(i4,180f5.1)') j,(a(i,j),i=1,il)
       enddo  !  j loop
      enddo   !  n loop
      return
      end

      subroutine int2x(a,zs,il,jl,glob,
     .                    rlong,rlat,addlong,addlat)
c     This one for 2 degree spacing and lats from -90 to 90
c     This routine:
c       b) does quadratic interpolation to find values for a(i,j)
c          wherever zs(i,j)>0  (a not touched where zs<0)
c       c) N.B. interpolation is at (long,lat) positions given by the
c         arrays rlong(il*jl) and rlat(il*jl) with displacements
c         addlong and addlat
      parameter (longs=180,lats=91)
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
c       do long=1,longs/2
c        lat=1
c        ll=long+(lat-1)*longs
c        is(ll)=ll +longs/2
c        is(ll+longs/2)=ll
c        lat=lats
c        ll=long+(lat-1)*longs
c        in(ll)=ll +longs/2
c        in(ll+longs/2)=ll
c       enddo  !  long loop
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
          rlongx=.5*(rlong(iq)+addlong)    ! 2 degree global grid
          rlatx=.5*(rlat(iq)+addlat)
          long=nint(rlongx)
          lat=nint(rlatx)
          lat=max(lat,2)        ! for south pole
          lat=min(lat,lats-1)   ! for north pole
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
          if(i.ge.17.and.i.le.21.and.j.ge.17.and.j.le.21)then
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

      include 'setxyz.f'      ! for conformal-cubic
      include 'jimcc.f'
c     include 'lconset.f'     ! for DARLAM
