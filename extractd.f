
c     Check last 3 lines before compiling
c     this is usually run on sphere
!     for big domains, put topo file and output on /temp
!     also compile and put a.out on /temp
c     this checks for land mask >=.5
c     prev stopped if inconsistency in land/sea specification
      program extractd
c     Program to extract vegetation types from 1 deg data set
c     also soil types,  albedo, rsmin, roughness lengths
c     and radon, CO2
c     for either DARLAM or conformal-cubic model;  jlm Mon  10-14-1996
c     N.B. these data sets range from -180 to 180 (approx)
c     N.B. for veg & soil from -180.5 will get silly values for rlong=180
c       hence have now chosen these to go from -180
      include 'newmpar.h'
      include 'dates.h'   ! to pass ds from setxyz
      include 'latlong.h'  ! rlatt,rlongg
      include 'mapproj.h'
      include 'parm.h'
      common/work/zs(il,jl),iarr(il,jl),iglobdat(360,180)
     .            ,a(il,jl),globdat(360,180),iglobex(0:361,0:181)
      common/latlong/rlat(ifull),rlong(ifull)
      character*80 topofile,soilfile,vegfile,soilout,vegout,radonfile
     .            ,co2file
     .            ,albfile,albout,rsminout,roughout,radonout,co2out
      character header*47
      character*9 formout1,formout2,formout3,formout4
      data pi/3.1415926536/
      namelist/vegnml/
     .   soilfile,vegfile,albfile,radonfile        ! global input data files
     .   ,co2file
     .                ,topofile                    ! model topography file
     .                ,soilout,vegout              ! model output data files
     .                ,albout,rsminout,roughout,radonout,co2out
     .                ,nsoil,nveg                  ! switches: 1 to do it
!                                    ! switch =2 adds Dean's veg on top
     .                ,nalb,nrsmin,nrough,nradon,nco2,id,jd
c    .                ,ds,du,tanl,rnml,stl1,stl2

      id=3
      jd=102
      open (85,file='extractd.nml',status='old')
      read (85,vegnml)
      write (6,vegnml)
      open (unit=51,file=topofile,form='formatted',status='old')
      idjd=id+il*(jd-1)

      if(npanels.eq.0)then   ! for DARLAM
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
        call lconset(ds)
        do j=1,jl
         do i=1,il
          iq=i+(j-1)*il
          call lconll(rlong(iq),rlat(iq),real(i),real(j))
c         N.B. these data sets range from -180 to 180
          if(rlong(iq).ge.180.)rlong(iq)=rlong(iq)-360.
         enddo  ! i loop
        enddo   ! j loop
c       find corner of lam grid as lat/longs for printouts
        lat1=nint(rlat(1)+90.5)          ! rlat(1) at -89.5
        long1=nint(rlong(1)+181.)        ! rlong(1) at -180.
        lat2=nint(rlat(ifull)+90.5)      ! rlat(1) at -89.5
        long2=nint(rlong(ifull)+181.)    ! rlong(1) at -180.
        jaaa=1
        jbb=jl
      else           !  (npanels.ne.0)
!       read header of topo file and check its dimensions
!       read(51,'(i3,i4,2f6.1,f5.2,f9.0,a47)')
        read(51,'(i3,i4,2f6.1,f6.3,f8.0,a47)')
     .            ilx,jlx,rlong0,rlat0,schmidt,dsx,header
        print *,ilx,jlx,rlong0,rlat0,schmidt,dsx,header
        if(ilx.ne.il.or.jlx.ne.jl)
     .              stop 'wrong topo file supplied'
        call setxyz
c       convert lat/longs to degrees in this routine
        do iq=1,ifull
         rlat(iq)=rlatt(iq)*180./pi
         rlong(iq)=rlongg(iq)*180./pi
c        N.B. these data sets range from -180 to 180
         if(rlong(iq).ge.180.)rlong(iq)=rlong(iq)-360.
        enddo
        lat1=1       ! lat1, lat2, long1, long2 for diag prints
        lat2=180
        long1=110
        long2=160
      endif

      print *,'id,jd,rlong,rlat',id,jd,rlong(idjd),rlat(idjd)
      print *,'rlat(1),rlat(ifull),rlong(1),rlong(ifull): '
     .        ,rlat(1),rlat(ifull),rlong(1),rlong(ifull)
      print *,'for printouts lat1,lat2,long1,long2: '
     .                      ,lat1,lat2,long1,long2
      read (51,*) zs
      print *,'zs(id,jd) ',zs(id,jd)
      read (51,*) zs   ! second read gives real land-sea mask 1. to 0.
      print *,'land_mask(id,jd) ',zs(id,jd)
c     work through i,j, values of the grid
      do j=1,jl
       do i=1,il
        iq=i+(j-1)*il
c       set lmask (here 1 for land, 0 for sea), just for diag printing
        iarr(i,j)=0
c       if(zs(i,j).gt.0.)iarr(i,j)=1
        if(zs(i,j).ge.0.5)iarr(i,j)=1
       enddo  ! i loop
      enddo   ! j loop
c     print out land-sea mask data for model grid by panel number
      do n=npanels,0,-1
       if(npanels.gt.0)then
         jaaa=n*il +1
         jbb=(n+1)*il
       endif
       print *,'land-sea mask for panel ',n
       do j=jbb,jaaa,-1
        write(6,'(i4,1x,180i1)') j,(iarr(i,j),i=1,il)
       enddo  !  j loop
      enddo   !  n loop

      write (formout1,'(1h(,i3,2hi3,1h))')il   !  i.e. (<il>i3)
      if(nveg.ge.1)then
        open (unit=10,file=vegfile,status='old')
        open (unit=20,file=vegout,status='unknown')
        call vegtype(lat1,lat2,long1,long2,jaaa,jbb,nveg)
        if(npanels.gt.0)
     .   write(20,'(i3,i4,2f6.1,f6.3,f8.0,''  veg type'')')
     .            il,jl,rlong0,rlat0,schmidt,ds
        write(20,formout1) iarr     ! write vegetation type
        close(10)
        close(20)
      endif

      if(nsoil.ge.1)then
        open (unit=10,file=soilfile,status='old')
        open (unit=20,file=soilout,status='unknown')
        call soiltype(lat1,lat2,long1,long2,jaaa,jbb)
        if(npanels.gt.0)
     .   write(20,'(i3,i4,2f6.1,f6.3,f8.0,''  soil type'')')
     .            il,jl,rlong0,rlat0,schmidt,ds
        write(20,formout1) iarr     ! write soil type
        close(10)
        close(20)
      endif

c     N.B. albedo, rsmin and roughness are on same file (by month)
      if(nalb.eq.1)then
        write (formout2,'("(",i3,"f4.0)" )')il   !  i.e. (<il>f4.0)
        open (unit=10,file=albfile,status='old')
        open (unit=20,file=albout,status='unknown')
        call albedo(lat1,lat2,long1,long2,jaaa,jbb)   ! as percent
        if(npanels.gt.0)
     .   write(20,'(i3,i4,2f6.1,f6.3,f8.0,''  albedo'')')
     .            il,jl,rlong0,rlat0,schmidt,ds
        write(20,formout2) a     ! write albedo
        close(20)
      endif

      if(nrsmin.eq.1)then
        write (formout3,'("(",i3,"f5.0)")' )il   !  i.e. (<il>f5.0)
        open (unit=20,file=rsminout,status='unknown')
        call rsmin(lat1,lat2,long1,long2,jaaa,jbb)
        if(npanels.gt.0)
     .   write(20,'(i3,i4,2f6.1,f6.3,f8.0,''  rsmin'')')
     .            il,jl,rlong0,rlat0,schmidt,ds
        write(20,formout3) a     ! write rsmin
        close(20)
      endif

      if(nrough.eq.1)then
        write ( formout4,'("(",i3,"f6.0)")' )il   !  i.e. (<il>f6.0)
        open (unit=20,file=roughout,status='unknown')
        call roughness(lat1,lat2,long1,long2,jaaa,jbb)
        if(npanels.gt.0)
     .   write(20,'(i3,i4,2f6.1,f6.3,f8.0,''  roughness (cm)'')')
     .            il,jl,rlong0,rlat0,schmidt,ds
        write(20,formout4) a     ! write roughness (cm)
c       print *,'formout4:',formout4
c       print *,'roughness for j=12:',(a(i,12),i=1,il)
c       atmos had a bug with writing 99.6 in f6.0 in Dec '96
        close(20)
      endif

      if(nradon.eq.1)then
        write (formout3,'("(",i3,"f5.2)")' )il   !  i.e. (<il>f5.2)
        open (unit=10,file=radonfile,status='old')
        open (unit=20,file=radonout,status='unknown')
        call radon(lat1,lat2,long1,long2,jaaa,jbb)
        if(npanels.gt.0)
     .   write(20,'(i3,i4,2f6.1,f6.3,f8.0,''  radon'')')
     .            il,jl,rlong0,rlat0,schmidt,ds
        write(20,formout3) a     ! write radon
        close(20)
      endif

      if(nco2.eq.1)then
        write (formout3,'("(",i3,"i5)")' )il   !  i.e. (<il>i5)
        open (unit=10,file=co2file,status='old')
        print *,'about to open co2file'
        read (10,*) header  ! header record
        print *,header
        open (unit=20,file=co2out,status='unknown')
        call co2(lat1,lat2,long1,long2,jaaa,jbb,ds)
        if(npanels.gt.0)
     .   write(20,'(i3,i4,2f6.1,f6.3,f8.0,''  co2'')')
     .             il,jl,rlong0,rlat0,schmidt,ds
        write(20,formout3) iarr     ! write co2
        close(20)
      endif
      end

      subroutine vegtype(lat1,lat2,long1,long2,jaaa,jbb,nveg)
      include 'newmpar.h'
      include 'latlong.h'  ! rlatt,rlongg
      include 'parm.h'
      common/latlong/rlat(ifull),rlong(ifull)
      common/work/zs(il,jl),iarr(il,jl),iglobdat(360,180)
     .            ,a(il,jl),globdat(360,180),iglobex(0:361,0:181)
c     common/workdean/dean(860),nsum(il,jl,0:31) ! 44 with old ones
      common/workdean/idean(860),nsum(il,jl,0:31) ! 44 with old ones
      character*4 headr
      data num/0/
c     read in 1 degree global vegetation type data set
      read(10,*) iglobdat  !  actually (18i4) format
c     1     `print out part of global data
      print *,' global vegetation type data'
      do lat=lat2,lat1,-1
        write(6,'(180z1)') (iglobdat(long,lat),long=long1,long2)
      enddo
c     fill in some missing points
c     iglobdat(1,155)=iglobdat(360,156)   ! (-180,65) Alaska
c     iglobdat(73,16)=iglobdat(73,15)     ! (-108,-74)  Antarctic
c     iglobdat(94,134)=iglobdat(93,134)   ! over great lakes
c     iglobdat(107,37)=iglobdat(109,37)   ! next around S. America
c     iglobdat(107,38)=iglobdat(109,38)
c     iglobdat(109,61)=iglobdat(110,60)
c     iglobdat(110,67)=iglobdat(111,66)
c     iglobdat(110,70)=iglobdat(112,70)
c     iglobdat(110,71)=iglobdat(112,71)
c     iglobdat(99,86)=iglobdat(101,86)
c     iglobdat(103,95)=iglobdat(105,95)
c     iglobdat(103,96)=iglobdat(105,96)
c     iglobdat(105,101)=iglobdat(106,100)
c     iglobdat(104,159)=iglobdat(104,161)     ! needed for C30
c     iglobdat(256,101)=iglobdat(258,101)     ! needed for C30
c     iglobdat(281, 88)=iglobdat(283, 88)     ! needed for C30
      iglobdat(341, 80)=iglobdat(329, 81)     ! Solomon Is = New Guinea
      do i=322,332
       print *,'i, veg ',i,(iglobdat(i,j),j=78,88)
      enddo

c     extend global data
      do lat=1,180
       do long=1,360
        iglobex(long,lat)=iglobdat(long,lat)
       enddo  ! long loop
       iglobex(0,lat)=iglobdat(360,lat)
       iglobex(361,lat)=iglobdat(1,lat)
      enddo  ! lat loop
      do long=0,360/2
       iglobex(long,0)=iglobdat(long+180,1)
       iglobex(long+181,0)=iglobdat(long+1,1)
       iglobex(long,181)=iglobdat(long+180,180)
       iglobex(long+181,181)=iglobdat(long+1,180)
      enddo  ! long loop

c     work through i,j, values of the grid
      do j=1,jl
       do i=1,il
        iarr(i,j)=0
c       if(zs(i,j).gt.0.)then
        if(zs(i,j).ge.0.5)then    ! checking mask (not zs)
          iq=i+(j-1)*il
          rlatx=rlat(iq)+90.5       ! rlat(1) at -89.5
          rlongx=rlong(iq)+181.     ! rlong(1) at -180.
c         rlongx=rlong(iq)+181.5    ! rlong(1) at -180.5
          lat=rlatx
          long=rlongx
          p=rlongx-long
          q=rlatx-lat
c         choose closest of 49 neighbours, having valid value
          distmin=1.e6
          do jj=-3,3
           do ii=-3,3
            dist=(ii-p)**2+(jj-q)**2
            if(dist.lt.distmin.and.iglobex(long+ii,lat+jj).gt.0)then
              iarr(i,j)=iglobex(long+ii,lat+jj)
              distmin=dist
            endif
           enddo
          enddo
	   if(i.eq.id.and.j.eq.jd)then
	    print *,'land_mask,vegtype ',zs(i,j),iarr(i,j)
            print *,'4 neighbs: ',iglobex(long,lat),iglobex(long+1,lat)
     .             ,iglobex(long,lat+1),iglobex(long+1,lat+1)
	   endif
          if(iarr(i,j).eq.0)then   ! already a land point
            print *,'veg type not assig; long,lat:',long,lat
            print *,'i,j,iq,rlong(iq),rlat(iq),lsea mask'
     .              ,i,j,iq,rlong(iq),rlat(iq),zs(i,j)
            print *,'4 neighbs: ',iglobex(long,lat),iglobex(long+1,lat)
     .             ,iglobex(long,lat+1),iglobex(long+1,lat+1)
            print *,'neighbs (-1,0) (-1,1) (2,0) (2,1):'
     .           ,iglobex(long-1,lat),iglobex(long-1,lat+1)
     .           ,iglobex(long+2,lat),iglobex(long+2,lat+1)
            print *,'neighbs (0,-1) (1,-1) (0,2) (1,2):'
     .           ,iglobex(long,lat-1),iglobex(long+1,lat-1)
     .           ,iglobex(long,lat+2),iglobex(long+1,lat+2)
            print *,'neighbs (-1,-1) (-1,2) (2,-1) (2,2):'
     .           ,iglobex(long-1,lat-1),iglobex(long-1,lat+2)
     .           ,iglobex(long+2,lat-1),iglobex(long+2,lat+2)
            print *,'normally  would stop in vegie routine, but now
     .              treat as small tropical island'
c           islands such as New Hebrides not in data set
            iarr(i,j)=1
          endif
        endif
       enddo  ! i loop
      enddo   ! j loop

c     print out Sib veg type for model grid by panel number
      print *
      do n=npanels,0,-1
       if(npanels.gt.0)then
         jaaa=n*il +1
         jbb=(n+1)*il
       endif
       print *,'veg type for panel ',n
       do j=jbb,jaa,-1
        write(6,'(i4,1x,180z1)') j,(iarr(i,j),i=1,il)
       enddo  !  j loop
      enddo   !  n loop

      do iq=1,ifull        ! add offset to Sib values so 1-13 becomes 32-44
       if(iarr(iq,1).gt.0)iarr(iq,1)=iarr(iq,1)+31
      enddo

      if(nveg.eq.1)return
!     initialize summing array
      do ntype=0,31
       do iq=1,ifull
        nsum(iq,1,ntype)=0
       enddo
      enddo
!     now do Dean's in array 860x700, from 112E to 155E, -45S to -10S
c     open (unit=11,file='/sphere/home/kow014/dean/prog/vtype31.grd',
      open (unit=11,file='dean31_int',
     .      form='formatted',status='old')
c     read (11,'(a4)') headr
c     read (11,*) headr
      read (11,*) ni1,nj1
      read (11,*) ni2,nj2
      read (11,*) ni3,nj3
      read (11,*) ri4,rj4
      print *,'ni1,nj1,ni2,nj2,ni3,nj3,ri4,rj4 ',
     .         ni1,nj1,ni2,nj2,ni3,nj3,ri4,rj4
      do jj=1,700
c      print *,'jj = ',jj
       if(jj.eq.191)num=0
c      read (11,*) (dean(ii),ii=1,860)   ! floating point input
       read (11,*) (idean(ii),ii=1,860)   ! 25i3 input
       dlat=-45. -.025 +jj*.05
       do ii=1,860
        dlon=112. -.025 +ii*.05
c       ndean=nint(dean(ii))
        ndean=idean(ii)
        if(jj.eq.192)print *,'ii,ndean,dlat,dlon ',ii,ndean,dlat,dlon
        if(jj.eq.192)print *,'ndlat,ndlon ',
     .      nint(1000.*dlat),nint(1000.*dlon)
     
        if(npanels.eq.0)then   ! for DARLAM
!         call lconij(dlon,dlat,xout,yout,theta)  ! new one
          call lconij(dlon,dlat,xout,yout)        ! old one
          i=nint(xout)
          j=nint(yout)
        else           !  (npanels.ne.0)
          call latltoij(dlon,dlat,xout,yout,nf)
          i=nint(xout)
          j=nint(yout)+nf*il
        endif
        if(ndean.lt.0.or.ndean.gt.31)then
          print *,'ii,jj,idean,ndean ',ii,jj,idean(ii),ndean
          stop
        endif
        nsum(i,j,ndean)=nsum(i,j,ndean)+1
        if(idean(ii).ne.0)then
          if(num.lt.99.or.(i.eq.id.and.j.eq.jd))then
           num=num+1
           print *,'ii,jj,i,j,ndean,idean(ii),nsum ',
     .              ii,jj,i,j,ndean,idean(ii),nsum(i,j,ndean)
          endif
        endif
       enddo
      enddo

!     now update vegie type (in iarr) by most popular one
      do iq=1,ifull
       ntypmax=0
       nsummax=0
       do ntype=0,31
        if(nsum(iq,1,ntype).gt.nsummax)then
          nsummax=nsum(iq,1,ntype)
          ntypmax=ntype
        endif
       enddo
c      if(iarr(iq,1).gt.0.and.ntypmax.gt.0)then
       if(ntypmax.gt.0)then
         print *,'iq,iarr(iq,1),ntypmax,nsummax ',
     .            iq,iarr(iq,1),ntypmax,nsummax
         iarr(iq,1)=ntypmax
       endif
      enddo   ! iq loop
      close (11)
      return
      end

      subroutine soiltype(lat1,lat2,long1,long2,jaa,jbb)
      include 'newmpar.h'
      include 'parm.h'
      include 'latlong.h'  ! rlatt,rlongg
      common/latlong/rlat(ifull),rlong(ifull)
      common/work/zs(il,jl),iarr(il,jl),iglobdat(360,180)
     .            ,a(il,jl),globdat(360,180),iglobex(0:361,0:181)
c     read in 1 degree global soil type data set
c     N.B. soil data needs to be read in from north to south!!!!!!!!!!!
      read(10,*) ((iglobdat(long,lat),long=1,360),lat=180,1,-1)
c     print out part of global data
      print *,' global soil type data'
      do lat=lat2,lat1,-1
        write(6,'(180i1)') (iglobdat(long,lat),long=long1,long2)
      enddo
c     fill in some missing points
c     iglobdat(73,16)=iglobdat(73,15)
c     iglobdat(131,60)=iglobdat(130,61)
c     iglobdat(216,69)=iglobdat(215,69)
c     iglobdat(142,74)=iglobdat(141,74)
c     iglobdat(142,76)=iglobdat(141,76)
c     iglobdat(317,76)=iglobdat(316,76)
c     iglobdat(220,83)=iglobdat(219,83)
c     iglobdat(220,85)=iglobdat(219,85)
c     iglobdat(109,102)=iglobdat(109,101)
c     iglobdat(261,102)=iglobdat(260,102)
c     iglobdat(281,102)=iglobdat(280,102)
c     iglobdat(261,105)=iglobdat(260,105)
c     iglobdat(104,108)=iglobdat(103,109)
c     iglobdat(110,132)=iglobdat(109,132)
c     iglobdat(318,135)=iglobdat(317,136)
c     iglobdat(99,143)=iglobdat(98,143)
c     iglobdat(95,146)=iglobdat(95,145)
c     iglobdat(86,150)=iglobdat(85,150)
c     iglobdat(1,155)=iglobdat(360,156)       ! (-180,65) Alaska
c     iglobdat(240,135)=iglobdat(240,134)     ! needed for C30
c     iglobdat(128, 57)=iglobdat(127, 57)     ! needed for C30
      iglobdat(341, 80)=iglobdat(329, 81)     ! Solomon Is = New Guinea

c     extend global data
      do lat=1,180
       do long=1,360
        iglobex(long,lat)=iglobdat(long,lat)
       enddo  ! long loop
       iglobex(0,lat)=iglobdat(360,lat)
       iglobex(361,lat)=iglobdat(1,lat)
      enddo  ! lat loop
      do long=0,360/2
       iglobex(long,0)=iglobdat(long+180,1)
       iglobex(long+181,0)=iglobdat(long+1,1)
       iglobex(long,181)=iglobdat(long+180,180)
       iglobex(long+181,181)=iglobdat(long+1,180)
      enddo  ! long loop

c     work through i,j, values of the grid
      do j=1,jl
       do i=1,il
        iarr(i,j)=0
c       if(zs(i,j).gt.0.)then
        if(zs(i,j).ge.0.5)then    ! checking mask (not zs)
          iq=i+(j-1)*il
          rlatx=rlat(iq)+90.5       ! rlat(1) at -89.5
          rlongx=rlong(iq)+181.     ! rlong(1) at -180.
c         rlongx=rlong(iq)+181.5    ! rlong(1) at -180.5
          lat=rlatx
          long=rlongx
          p=rlongx-long
          q=rlatx-lat
c         choose closest of 49 neighbours, having valid value
          distmin=1.e6
          do jj=-3,3
           do ii=-3,3
            dist=(ii-p)**2+(jj-q)**2
            if(dist.lt.distmin.and.iglobex(long+ii,lat+jj).gt.0)then
              iarr(i,j)=iglobex(long+ii,lat+jj)
              distmin=dist
            endif
           enddo
          enddo
	   if(i.eq.id.and.j.eq.jd)then
	    print *,'land_mask,soiltype ',zs(i,j),iarr(i,j)
            print *,'4 neighbs: ',iglobex(long,lat),iglobex(long+1,lat)
     .             ,iglobex(long,lat+1),iglobex(long+1,lat+1)
	   endif
          if(iarr(i,j).eq.0)then   ! already a land point
            print *,'soil type not assig; long,lat:',long,lat
            print *,'i,j,iq,rlong(iq),rlat(iq) '
     .              ,i,j,iq,rlong(iq),rlat(iq)
            print *,'4 neighbs: ',iglobex(long,lat),iglobex(long+1,lat)
     .             ,iglobex(long,lat+1),iglobex(long+1,lat+1)
            print *,'neighbs (-1,0) (-1,1) (2,0) (2,1):'
     .           ,iglobex(long-1,lat),iglobex(long-1,lat+1)
     .           ,iglobex(long+2,lat),iglobex(long+2,lat+1)
            print *,'neighbs (0,-1) (1,-1) (0,2) (1,2):'
     .           ,iglobex(long,lat-1),iglobex(long+1,lat-1)
     .           ,iglobex(long,lat+2),iglobex(long+1,lat+2)
            print *,'neighbs (-1,-1) (-1,2) (2,-1) (2,2):'
     .           ,iglobex(long-1,lat-1),iglobex(long-1,lat+2)
     .           ,iglobex(long+2,lat-1),iglobex(long+2,lat+2)
            print *,'normally  would stop in soil routine, but now
     .              treat as small tropical island'
            iarr(i,j)=1
          endif
        endif
       enddo  ! i loop
      enddo   ! j loop

c     print out soil type for model grid by panel number
      print *
      do n=npanels,0,-1
       if(npanels.gt.0)then
         jaa=n*il +1
         jbb=(n+1)*il
       endif
       print *,'soil type for panel ',n
       do j=jbb,jaa,-1
        write(6,'(i4,1x,180i1)') j,(iarr(i,j),i=1,il)
       enddo  !  j loop
      enddo   !  n loop
      return
      end

      subroutine albedo(lat1,lat2,long1,long2,jaa,jbb)
      include 'newmpar.h'
      include 'latlong.h'  ! rlatt,rlongg
      common/latlong/rlat(ifull),rlong(ifull)
      common/work/zs(il,jl),iarr(il,jl),iglobdat(360,180)
     .            ,a(il,jl),globdat(360,180),iglobex(0:361,0:181)
c     read in 1 degree global albedo data set; over sea: 99.9
      read (10,*) imon,ivar
      print *,'albedo data  imon, ivar = ',imon,ivar
      read(10,'(30f6.2)') ((globdat(long,lat),lat=1,180),long=1,360)

c     print out part of global data
      print *,' global albedo data'
      print *,' typical values: ',(globdat(lat,lat),lat=1,180)
      do lat=lat2,lat1,-1
        do long=long1,long2
         iglobdat(long,lat)=.1*globdat(long,lat)
        enddo
        write(6,'(180i1)') (iglobdat(long,lat),long=long1,long2)
      enddo

      call fillintg(a,zs,il,jl,globdat,99.9,
     .          rlong,rlat,180.5,90.5)   ! rlong(1) at -179.5; rlat(1) at -89.5
c     print out part of global data
      print *,' extended global albedo data'
      do lat=lat2,lat1,-1
        do long=long1,long2
         iglobdat(long,lat)=.1*globdat(long,lat)
        enddo
        write(6,'(180i1)') (iglobdat(long,lat),long=long1,long2)
      enddo

      do iq=1,il*jl
c       if(zs(iq,1).le.0.) a(iq,1)=11.    ! preset for sea points
        if(zs(iq,1).lt.0.5) a(iq,1)=11.   ! preset for sea points
        iarr(iq,1)=.1*a(iq,1)               ! just for diag prints
      enddo

c     print out albedo for model grid by panel number
      print *
      do n=npanels,0,-1
       if(npanels.gt.0)then
         jaa=n*il +1
         jbb=(n+1)*il
       endif
       print *,'albedo data for panel ',n
       do j=jbb,jaa,-1
        write(6,'(i4,1x,180i1)') j,(iarr(i,j),i=1,il)
       enddo  !  j loop
      enddo   !  n loop
      return
      end

      subroutine radon(lat1,lat2,long1,long2,jaa,jbb)
      include 'newmpar.h'
      include 'latlong.h'  ! rlatt,rlongg
      common/latlong/rlat(ifull),rlong(ifull)
      common/work/zs(il,jl),iarr(il,jl),iglobdat(360,180)
     .            ,a(il,jl),globdat(360,180),iglobex(0:361,0:181)
c     read in 1 degree global radon emission data set (x100)
!     read (10,*) imon,ivar  ! no header for radonfile
!     print *,'radon data  imon, ivar = ',imon,ivar
      read(10,*)((globdat(long,lat),lat=1,180),long=1,360)

c     scale to actual values
      do lat=1,180
       do long=1,360
        globdat(long,lat)=.01*globdat(long,lat)   ! scaling radon data
       enddo  ! long loop
      enddo  ! lat loop

c     print out part of global data
      print *,' global radon data'
      print *,' typical values: ',(globdat(lat,lat),lat=1,180)
      do lat=lat2,lat1,-1
        do long=long1,long2
         iglobdat(long,lat)=nint(globdat(long,lat))
        enddo
        write(6,'(180i1)') (iglobdat(long,lat),long=long1,long2)
      enddo

      call fillintg(a,zs,il,jl,globdat,999.,
     .          rlong,rlat,180.5,90.5)   ! rlong(1) at -179.5; rlat(1) at -89.5
c     print out part of global data
      print *,' extended global radon data'
      do lat=lat2,lat1,-1
        do long=long1,long2
         iglobdat(long,lat)=nint(globdat(long,lat))
        enddo
        write(6,'(180i1)') (iglobdat(long,lat),long=long1,long2)
      enddo

      do iq=1,il*jl
        if(zs(iq,1).lt.0.5)then
          a(iq,1)=0.  ! preset for sea points
        else
          a(iq,1)=max(a(iq,1),.49)  ! radon preset for missed land points
        endif
        iarr(iq,1)=nint(a(iq,1))             ! just for diag prints
      enddo

c     print out radon for model grid by panel number
      print *
      do n=npanels,0,-1
       if(npanels.gt.0)then
         jaa=n*il +1
         jbb=(n+1)*il
       endif
       print *,'radon  for panel ',n
       do j=jbb,jaa,-1
        write(6,'(i4,1x,180i1)') j,(iarr(i,j),i=1,il)
       enddo  !  j loop
      enddo   !  n loop
      return
      end

      subroutine co2(lat1,lat2,long1,long2,jaa,jbb,ds)
      include 'newmpar.h'
      include 'latlong.h'  ! rlatt,rlongg
      common/latlong/rlat(ifull),rlong(ifull)
      include 'dates.h'    ! rearth
      include 'map.h'      ! em
      common/work/zs(il,jl),iarr(il,jl),iglobdat(360,180)
     .            ,a(il,jl),globdat(360,180),iglobex(0:361,0:181)
      data pi/3.1415926536/
c     read in 1 degree global co2 emission data set (x100)
!     read (10,*) imon,ivar  ! no header for co2file
!     print *,'co2 data  imon, ivar = ',imon,ivar
      read(10,*) globdat

c     print out part of global data
      print *,' global co2 data'
      print *,'rearth,ds: ',rearth,ds
      print *,' typical values: ',(globdat(lat,lat),lat=1,180)
      do lat=lat2,lat1,-1
        do long=long1,long2
         iglobdat(long,lat)=nint(globdat(long,lat))
        enddo
        write(6,'(180i1)') (iglobdat(long,lat),long=long1,long2)
      enddo

      do iq=1,ifull
        iarr(iq,1)=0
      enddo

c     scale to actual values and assign to nearest point
c     N.B. rlong(1) at -179.5; rlat(1) at -89.5 according to Marland
      do lat=1,180
       do long=1,360
        if(globdat(long,lat).gt.0.)then
          xlong=long-180.5
          xlat=lat-90.5
          if(long.eq.325.and.lat.eq.53)then
            xlong=145.
            xlat=-37.85
            print *,'changing Melb from .. to ',
     .       long-180.5,lat-90.5,xlong,xlat
          endif
          call latltoij(xlong,xlat,xout,yout,nf)
          i=nint(xout)
          j=nint(yout)+nf*il
c         if(xlong.gt.144..and.xlong.lt.155.
c    .                 .and.xlat.gt.-39..and.xlat.lt.-28.)then
c           print *,'co2 industrial source of ',globdat(long,lat)
c           print *,'for long,lat ',long,lat
c           print *,'for xout,yout,nf ',xout,yout,nf
c           print *,'and i,j ',i,j
c         endif  !  (xlong.gt.144..and.xlong.lt.155....
!         scale 1 degree data to grid area (to preserve net source)
          source=cos(pi*xlat/180.) *(pi*rearth/180.)**2 /(ds/em(i,j))**2
     .           *globdat(long,lat)
          if(globdat(long,lat).gt..5)
     .           print 9,i,j,em(i,j),globdat(long,lat),source
9         format('co2 i,j,em,globdat,source',2i5,3f9.4)
!         iarr(i,j)=iarr(i,j)+nint(100.*globdat(long,lat))  ! accumulate sources
          iarr(i,j)=iarr(i,j)+nint(100.*source)             ! accumulate sources
          if(xlong.gt.110..and.xlong.lt.180..and.
     .       xlat.lt.-10.)then
            print *,'Australian source for xlong,xlat ',xlong,xlat
            print *,'i,j,globdat,iarr ',i,j,globdat(long,lat),iarr(i,j)
          endif
        endif
       enddo  ! long loop
      enddo  ! lat loop

c     print out co2 for model grid by panel number
      print *
      do n=npanels,0,-1
       if(npanels.gt.0)then
         jaa=n*il +1
         jbb=(n+1)*il
       endif
       print *,'co2  for panel (divided by 10) ',n
       do j=jbb,jaa,-1
        write(6,'(i4,1x,180i1)') j,(nint(.1*iarr(i,j)),i=1,il)
       enddo  !  j loop
      enddo   !  n loop
      return
      end

      subroutine rsmin(lat1,lat2,long1,long2,jaa,jbb)
      include 'newmpar.h'
      include 'latlong.h'  ! rlatt,rlongg
      common/latlong/rlat(ifull),rlong(ifull)
      common/work/zs(il,jl),iarr(il,jl),iglobdat(360,180)
     .            ,a(il,jl),globdat(360,180),iglobex(0:361,0:181)
c     read in 1 degree global rsmin data set; over sea: 999
      read (10,*) imon,ivar
      print *,'rsmin data  imon, ivar = ',imon,ivar
      read(10,'(30f6.2)')((globdat(long,lat),lat=1,180),long=1,360)

c     print out part of global data
      print *,' global rsmin data'
      print *,' typical values: ',(globdat(lat,lat),lat=1,180)
      do lat=lat2,lat1,-1
        do long=long1,long2
         iglobdat(long,lat)=.02*globdat(long,lat)
        enddo
        write(6,'(180i1)') (iglobdat(long,lat),long=long1,long2)
      enddo

      call fillintg(a,zs,il,jl,globdat,999.,
     .          rlong,rlat,180.5,90.5)   ! rlong(1) at -179.5; rlat(1) at -89.5
c     print out part of global data
      print *,' extended global rsmin data'
      do lat=lat2,lat1,-1
        do long=long1,long2
         iglobdat(long,lat)=.02*globdat(long,lat)
        enddo
        write(6,'(180i1)') (iglobdat(long,lat),long=long1,long2)
      enddo

      do iq=1,il*jl
c       if(zs(iq,1).le.0.) a(iq,1)=-1.   ! preset for sea points
        if(zs(iq,1).lt.0.5) a(iq,1)=-1.  ! preset for sea points
        iarr(iq,1)=.02*a(iq,1)             ! just for diag prints
      enddo

c     print out rsmin for model grid by panel number
      print *
      do n=npanels,0,-1
       if(npanels.gt.0)then
         jaa=n*il +1
         jbb=(n+1)*il
       endif
       print *,'rsmin (*.02) for panel ',n
       do j=jbb,jaa,-1
        write(6,'(i4,1x,180i1)') j,(iarr(i,j),i=1,il)
       enddo  !  j loop
      enddo   !  n loop
      return
      end

      subroutine roughness(lat1,lat2,long1,long2,jaa,jbb)
      include 'newmpar.h'
      include 'latlong.h'  ! rlatt,rlongg
      common/latlong/rlat(ifull),rlong(ifull)
      common/work/zs(il,jl),iarr(il,jl),iglobdat(360,180)
     .            ,a(il,jl),globdat(360,180),iglobex(0:361,0:181)
c     read in 1 degree global roughness data set; over sea: 999.0
      read (10,*) imon,ivar
      print *,'roughness data  imon, ivar = ',imon,ivar
      read(10,'(30f6.2)')((globdat(long,lat),lat=1,180),long=1,360)

c     print out part of global data
      print *,' global roughness data'
      print *,' typical values: ',(globdat(lat,lat),lat=1,180)
      do lat=lat2,lat1,-1
        do long=long1,long2
         iglobdat(long,lat)=.1*globdat(long,lat)
        enddo
        write(6,'(180i1)') (iglobdat(long,lat),long=long1,long2)
      enddo

      call fillintg(a,zs,il,jl,globdat,999.,
     .          rlong,rlat,180.5,90.5)   ! rlong(1) at -179.5; rlat(1) at -89.5
c     print out part of global data
      print *,' extended global roughness data'
      do lat=lat2,lat1,-1
        do long=long1,long2
         iglobdat(long,lat)=.03*globdat(long,lat)
        enddo
        write(6,'(180i1)') (iglobdat(long,lat),long=long1,long2)
      enddo

      do iq=1,il*jl
c       if(zs(iq,1).le.0.) a(iq,1)=1.    ! preset for sea points
        if(zs(iq,1).lt.0.5) a(iq,1)=1.   ! preset for sea points
        iarr(iq,1)=.03*a(iq,1)             ! just for diag prints
      enddo

c     print out roughness for model grid by panel number
      print *
      do n=npanels,0,-1
       if(npanels.gt.0)then
         jaa=n*il +1
         jbb=(n+1)*il
       endif
       print *,'roughness (*.03 in cm) for panel ',n
       do j=jbb,jaa,-1
        write(6,'(i4,1x,180i1)') j,(iarr(i,j),i=1,il)
       enddo  !  j loop
      enddo   !  n loop
      return
      end
      subroutine fillintg(a,zs,il,jl,glob,value,
     .                    rlong,rlat,addlong,addlat)
c     This routine:
c       a) fills in global array glob (360x180) for missing values
c       b) does bi-linear interpolation to find values for a(i,j)
c          wherever zs(i,j)>0  (a not touched where zs<0)
c       c) N.B. interpolation is at (long,lat) positions given by the
c         arrays rlong(il*jl) and rlat(il*jl) with displacements
c         addlong and addlat
      parameter (longs=360,lats=180)
c     routine fills in interior of a global array which has undefined points
c     typically for arrays (360,180)
      real glob(longs*lats)      ! input and output global array
      real value                 ! array value denoting undefined
      real work(longs*lats)      ! global work array
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
      print *,'entering fillintg with value = ',value

2     nrem=0
      do ll=1,longs*lats
       work(ll)=glob(ll)
       if(glob(ll).eq.value)then
         neighb=0
         av=0.
         if(glob(ie(ll)).ne.value)then
           neighb=neighb+1
           av=av+glob(ie(ll))
         endif
         if(glob(ine(ll)).ne.value)then
           neighb=neighb+1
           av=av+glob(ine(ll))
         endif
         if(glob(ise(ll)).ne.value)then
           neighb=neighb+1
           av=av+glob(ise(ll))
         endif
         if(glob(in(ll)).ne.value)then
           neighb=neighb+1
           av=av+glob(in(ll))
         endif
         if(glob(iw(ll)).ne.value)then
           neighb=neighb+1
           av=av+glob(iw(ll))
         endif
         if(glob(inw(ll)).ne.value)then
           neighb=neighb+1
           av=av+glob(inw(ll))
         endif
         if(glob(isw(ll)).ne.value)then
           neighb=neighb+1
           av=av+glob(isw(ll))
         endif
         if(glob(is(ll)).ne.value)then
           neighb=neighb+1
           av=av+glob(is(ll))
         endif
         if(neighb.gt.0)then
           work(ll)=av/neighb
         else
           nrem=nrem+1    ! number of remaining points
         endif
       endif
      enddo
      do ll=1,longs*lats
       glob(ll)=work(ll)
      enddo
c     print *,'in fillintg, nrem = ',nrem
      if(nrem.gt.0)go to 2

c     --------------------------------------------------------------
c     now do conditional biperiodic interpolation
c     work through i,j, values of the grid
c     do j=1,jl
c      do i=1,il
c       if(zs(i,j).gt.0.)then
c         iq=i+(j-1)*il
c         centre around nearest lat/long point
c         rlongx=rlong(iq)+addlong    ! 180.5 for rlong(1) at -179.5
c         rlatx=rlat(iq)+addlat       ! 90.5 for rlat(1) at -89.5
c         long=nint(rlongx)
c         lat=nint(rlatx)
c         p=rlongx-long   ! for unit spacing of long
c         q=rlatx-lat     ! for unit spacing of lat
c         ll=long+(lat-1)*longs
c         aw=.5*q*(q-1.)*glob(isw(ll)) +(1.-q*q)*glob(iw(ll))
c    .      +.5*q*(q+1.)*glob(inw(ll))
c         a0=.5*q*(q-1.)*glob(is(ll)) +(1.-q*q)*glob(ll)
c    .      +.5*q*(q+1.)*glob(in(ll))
c         ae=.5*q*(q-1.)*glob(ise(ll)) +(1.-q*q)*glob(ie(ll))
c    .      +.5*q*(q+1.)*glob(ine(ll))
c         a(i,j)=.5*p*(p-1.)*aw +(1.-p*p)*a0 +.5*p*(p+1.)*ae
c       endif
c      enddo  ! i loop
c     enddo   ! j loop
c     --------------------------------------------------------------
c     now do conditional bilinear interpolation

c     work through i,j, values of the grid
      do j=1,jl
       do i=1,il
c       if(zs(i,j).gt.0.)then
        if(zs(i,j).ge.0.5)then
          iq=i+(j-1)*il
c         centre around lat/long point in sw of box
          rlongx=rlong(iq)+addlong    ! addlong=180.5 for rlong(1) at -179.5
          rlatx=rlat(iq)+addlat       ! addlat=  90.5 for rlat(1) at -89.5
          rlongx=max(rlongx,1.0001)   ! allow for western boundary
          rlatx=max(rlatx,1.0001)     ! allow for south pole
          long=int(rlongx)
          lat=int(rlatx)
c         long=max(int(rlongx),1) ! allow for western boundary
c         lat=max(int(rlatx),1)   ! allow for south pole
          p=rlongx-long   ! for unit spacing of long
          q=rlatx-lat     ! for unit spacing of lat
          ll=long+(lat-1)*longs
          if(iq.eq.2113)print *,'i,j,long,lat,p,q,zs ',
     .                         i,j,long,lat,p,q,zs(i,j)
          if(iq.eq.2113)print *,'i,j,rlong,addlong,rlat,addlat ',
     .                         i,j,rlong(iq),addlong,rlat(iq),addlat
          if(iq.eq.2113)print *,'i,j,ll,in(ll),ie(ll),ine(ll) ',
     .                         i,j,ll,in(ll),ie(ll),ine(ll)
          if(iq.eq.2113)print *,'glob values ',
     .                 glob(ll),glob(in(ll)),glob(ie(ll)),glob(ine(ll))
          a0=(1.-q)*glob(ll) + q*glob(in(ll))
          ae=(1.-q)*glob(ie(ll)) + q*glob(ine(ll))
          a(i,j)=(1.-p)*a0 +p*ae
          if(iq.eq.2113)print *,'i,j,a0,ae,a(i,j) ',
     .                         i,j,a0,ae,a(i,j)
        endif
       enddo  ! i loop
      enddo   ! j loop

      return
      end
      include 'setxyz.f'      ! for conformal-cubic
      include 'jimcc.f'       ! for conformal-cubic
      include 'latltoij.f'    ! for conformal-cubic
      include 'lconset.f'     ! for DARLAM

