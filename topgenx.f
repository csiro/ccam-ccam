!     topgenx is usually run on sphere (to have integer*2 available)
c     needs newmpar.h appropriate to cc model
c     and latltoij.f with ncray=1
c     Has been usually run on cherax
c           where topo data is on file
c                /home/csdar/cshld/src/topo/par_topo/topo2
c     output always goes to 'topout0'
c **  set rlong0, rlat0 and schmidt in data statement
c     output is always unfiltered; use topfilt.f to filter it

      program topgencc
c     modified 30-7-96 from topgen for C-C by J. McGregor
c     modified 4-27-92 by Dr. J. Katzfey
c     modified 4-23-91 by Dr. J. Katzfey
c     early version: Kevin Walsh, june 1990. modified from a routine
c     sptoar9 written by John McGregor,
c     CSIRO Atmospheric Research, Aspendale, Vic.

c     this routine calculates the land-sea mask and averaged
c     topography for the nested model from a high-resolution (5 min) lat-lon
c     gridded topography data set
c     Also writes variance (tsd) of orography

c***********************************************************************

c     input :
c           namelist topnml in file 'top.nml'

c           igrid = 0 for 5 min data (from harvey davies)
c                 = 1 for 3.75 min data
c                 = 2 for 1.5  min data
c                 = -1 for 10km nz data
c           pctw =.5      percent water for land mask
c           zmin =0  (int)minimum height used for topog averaging (m)
c           zslm =-5 (int)minimum height used to determine lmask (m)
c           ctop =100.    contour int. for topog. plots (m)
c           debug=.false. for extra printout
c             id,jd = lcon grid points for printouts
c           olam =.true.  whether want lam.conf. proj.(t) of lat.lon.(f)
c           if olam=t
c                 ds   =50000.  grid spacing (m)
c                 du   =36.5    # grid points between sw corner and rnml
c                 tanl =40.     latitude tangent at bottom at rnml (deg)
c                 rnml =150.    standard (vertical) longitude (deg)
c                 stl1 =40.     standard latitude 1 (deg,+=southern)
c                 stl2 =10.     standard latitude 2 (deg,+=southern)
c           if olam=f
c                 wbd  =149.    western longitude for l/l grid (deg)
c                 sbd  =-38.    southern latitude for l/l grid (deg,+=nothern)
c                 dlon =-1.     delta lon. for l/l grid (deg,if<0 then 5/60)
c                 dlat =-1.     delta lat. for l/l grid (deg,if<0 then 5/60)

c***********************************************************************

c     output : 5 unformatted arrays dim.(il+1,jl+1)
c          add = topog. (m)
c          rmsk = land mask (1=land, 0=ocean)
c          tsd = std.dev. of 5min topog in grid box (m)
c          tmax = maximum of 5min topog in grid box (m)
c          tmin = minimum of 5min topog in grid box (m)

c     output to be used by vidar routines for creating darlam I.C.
c     ****    not in cc version

c     parameter statements
      include 'newmpar.h'
      include 'dates.h'   ! to pass ds from setxyz
      include 'indices.h'  
      include 'parm.h'   ! to pass rlong0,rlat0,schmidt  to setxyz
      parameter ( nfix=17 )
      logical olam      

c     'global' topog. data
      real zs(18876), zmin, zslm
      integer*2 izs(4800)
      real glat(121*156), glon(121*156)
      real tfix(nfix)
      integer ipfix(nfix),jpfix(nfix)

c     various work arrays
      dimension rmsk(il,jl),inum(il,jl),inumx(il,jl),zss(il,jl)
      dimension add(il,jl),almsk(il,jl)
      dimension tmax(il,jl),tmin(il,jl),tsd(il,jl)
      common/filein/filein
      character*60 filein

      logical debug, osd
      character*60 fileout
      character*9 formout

      namelist / topnml / ds, du, tanl, rnml, stl1, stl2, debug, pctw
     .  ,luout,fileout   , olam, wbd, sbd, dlon, dlat, ctop, zmin, zslm
     .                  , igrid
     .                  , tfix,ipfix,jpfix,osd,id,jd
     .                  , rlong0,rlat0,schmidt

c     initial namelist
      data du/36.5/, tanl/40./, rnml/150./
      data stl1/40./, stl2/10./, debug/.false./, olam/.true./
      data wbd/149./, sbd/-38./, dlon/-1./, dlat/-1./
      data ctop/100./, zmin/0/, zslm/-5/, pctw/.5/, osd/.false./
      data igrid/0/
      data tfix/nfix*0./, ipfix/nfix*1/, jpfix/nfix*1/
      data luout/80/,fileout/'topout0'/
!     data rlong0/0./,rlat0/90./,schmidt/1./
c     data rlong0/134./,rlat0/-24./,schmidt/.5/

c     optionally open namelist file and read
      open ( unit=5,file='top.nml',status='unknown' )
      read ( 5,topnml, end=5 )
 5    write ( 6,topnml )

      call setxyz

c     nested model topography (output)
      print *,'open',luout,fileout
      open(luout,file=fileout,form='formatted',status='unknown')


c initialize min,max, and sd arrays
      do 12 j=1,jl
      do 12 i=1,il
        tmin(i,j)=99999.
        tmax(i,j)=-99999.
        tsd(i,j)=0.
        almsk(i,j)=0.
        inum(i,j)=0
 12   continue

c***********************************************************************

c     begin main 'global' topog loops

        print *,'***** igrid=',igrid

        if ( igrid.eq.0 ) then
          ngx=4320
          ngy=2160
          gslat=-90.
          gwlon=0.
c         grid increment in degrees 
          gxinc = 5./60.   !  (5 min)
          gyinc = 5./60.   !  (5 min)
        elseif ( igrid.eq.1 ) then
          ngy=6000  ! naming reversed here
          ngx=4800
          gslat=-10.00416666666667
          gwlon=140.00416666666666
          gxinc =  1./120.  ! .00833 going southwards
          gyinc = -1./120.  ! .00833 going southwards
!         open(unit=30,file='/home/csdar/csjjk/topgen/E140S10.DEM',
          open(unit=30,file='/sphere/home/kat024/newtopgen/E140S10.DEM',
     .      access='direct',recl=2*4800)
          print *,'file 30 opened'
        endif
        print *,'ngx,ngy,gslat,gwlon=',ngx,ngy,gslat,gwlon

c#######################################################################
        do jg=1,ngy
          aglat = gslat + (jg-1)*gyinc
          if ( igrid.eq.0 ) then
            stop 'call read_ht ( aglat, zs )'
          endif !( igrid.eq.0 ) then
          if ( igrid.eq.1 ) then
            read(30,rec=jg) (izs(ig),ig=1,ngx)
            print *,'jg,izs ',jg,(izs(i),i=1,8)
            do ig=1,ngx
             zs(ig)=izs(ig)  ! don't divide these by 10
            enddo !ig=1,ngx
          endif

c#######################################################################

          do ig=1,ngx
           aglon = gwlon+(ig-1)*gxinc

           if(jg.eq.5996)print *,'calling latltoij ig,aglon,aglat: ',
     .                                             ig,aglon,aglat
           call latltoij(aglon,aglat,alci,alcj,nface)  ! con-cubic/octagon
           lci = nint(alci)
           lcj = nint(alcj)
c          convert to "double" (i,j) notation
           lcj=lcj+nface*il
           if(lci.eq.2.and.lcj.eq.2)
     .       print *,'ig,jg,aglon,aglat,alci,alcj,lci,lcj,zs'
     .               ,ig,jg,aglon,aglat,alci,alcj,lci,lcj,zs(ig)

c          make sure point is within lam.conf./C-C grid
           if ( lci.gt.0 .and. lci.le.il
     .             .and. lcj.gt.0 .and. lcj.le.jl ) then

c                 code to handle the small area below sea level in the centre
c                 of Australia, near Darwin, and centre of New Zealand

                  if ( aglon.gt.135. .and. aglon.lt.145.
     .            .and. aglat.lt.-24. .and. aglat.gt.-32. ) then
c                    central Australia
c                    print *,'point in central Australia'
                     amask=1.
                     add(lci,lcj) = add(lci,lcj) + zs(ig)
                     uzs=zs(ig)
                  elseif ( aglon.gt.171. .and. aglon.lt.172.
     .            .and. aglat.lt.-43. .and. aglat.gt.-44. ) then
c                    New Zealand
                     amask=1.
                     add(lci,lcj) = add(lci,lcj) + zs(ig)
                     uzs=zs(ig)
                  elseif ( aglon.gt.131. .and. aglon.lt.135.
     .           .and. aglat.lt.-13. .and. aglat.gt.-16. ) then
c                    Darwin
                     amask=1.
                     add(lci,lcj) = add(lci,lcj) + zs(ig)
                     uzs=zs(ig)
                  else  ! not in central Australia
c                    all other points
c                    use zslm to define where sea begins
c                    (fixup to 5 min data coastline)
                     if ( zs(ig).lt.zslm ) then
c                       ocean point
                        amask = 0.
                     else  ! zs>=zslm
c                       land point
                        amask = 1.
                     end if  ! zs<zslm
c                    force topog >=zmin
                     uzs=max(zs(ig),zmin)
c                    accumulate topog. pnts
                     add(lci,lcj)=add(lci,lcj)+uzs  ! looks after -9999.
                  end if ! all cases
c                 accumulate lmask pnts
                  almsk(lci,lcj) = almsk(lci,lcj) + amask
c                 accumulate number of 5min pnts in grid box
                  inum(lci,lcj) = inum(lci,lcj) + 1
c                 find max/min topog. pnts in grid box
                  tmax(lci,lcj)=max(tmax(lci,lcj),uzs)
                  tmin(lci,lcj)=min(tmin(lci,lcj),uzs)
c                 sum of squares for sd calc.
                  tsd(lci,lcj)=tsd(lci,lcj)+uzs**2

                  if ( debug ) then
                     if ( lci.eq.id .and. lcj.eq.jd ) then
                 write(6,'("lci,lcj,ig,zs(ig),add(),almsk(), inum()="
     &            ,3i6,4f10.3,i6)') lci,lcj,ig,zs(ig),uzs,add(lci,lcj)
     &                         ,almsk(lci,lcj), inum(lci,lcj)
                     endif  ! selected points only
                  endif  ! debug

               endif ! lci.gt.0 .and. lci.le.il
c    .           .and. lcj.gt.0 .and. lcj.le.jl ) then

             enddo !ig=1,ngx
c#######################################################################

        enddo !jg=1,ngy
c#######################################################################

c     check topographic height over each grid square

      if ( debug ) then
         print *,'i,j,inum(),rnum,almsk(),rmsk(),add()'
      endif  ! debug

      numzer=0
      imin=9999
      jmin=9999
      imax=0
      jmax=0
      do 60 i=1,il
        do 60 j=1,jl
c         if some grid points are not within the topog. grid
c         assume that points are over ocean
          if ( inum(i,j).eq.0 ) then
             numzer=numzer+1
             imin=min(i,imin)
             imax=max(i,imax)
             jmin=min(j,jmin)
             jmax=max(j,jmax)
             add(i,j) = -.2   ! allows to see where data "window" is
             rmsk(i,j)= 0.
             tsd(i,j)= 0.
             tmin(i,j)= 0.
             tmax(i,j)= 0.
          else
c            set land-sea mask (rmsk)
             rnum=1./inum(i,j)
             almsk(i,j) = almsk(i,j)*rnum

c            if less than half of 5min pnts are ocean, assume it is land pnt
             if ( almsk(i,j).lt.pctw ) then
c               ocean point
                rmsk(i,j) = 0.
             else ! almsk
c               land point
                rmsk(i,j) = 1.
             endif ! almsk

             if ( rmsk(i,j).lt..5 ) then
c               if a grid-square is declared sea,
c               then set topographic height to zero;
                add(i,j) = 0.
             else ! rmsk
c               if it is declared land, then set topographic height (add)
c               use all points
                add(i,j) = add(i,j)*rnum
             endif ! rmsk

c            compute sd
             tsd(i,j)=sqrt(abs(tsd(i,j)*rnum-add(i,j)**2)) ! abs for rounding?
c            tsd(i,j)=sqrt(tsd(i,j)*rnum-add(i,j)**2)
          end if
          if ( debug ) then
            if ( i.eq.id .and. j.eq.jd ) then
               print *,i,j,inum(i,j),rnum,almsk(i,j)
     .                ,rmsk(i,j),add(i,j)
            endif  ! selected points only
          endif  ! debug
60    continue

      if(igrid.eq.0)then
2      if(numzer.gt.0.and.igrid.eq.0)then
!       this checks for non-data points, and inserts a neighbouring value
!       don't do it for high-res window data (i.e. if igrid=1)
        print *,'**** numzer= ',numzer
        numzer=0
        do j=1,jl
         do i=1,il
          iq=i+(j-1)*il
          inumx(i,j)=inum(i,j)
          if(inum(i,j).eq.0)then
            iq2=0
            if(inum(in(iq),1).ne.0)then
              iq2=in(iq)
            endif
            if(inum(ie(iq),1).ne.0)then
              iq2=ie(iq)
            endif
            if(inum(iw(iq),1).ne.0)then
              iq2=iw(iq)
            endif
            if(inum(is(iq),1).ne.0)then
              iq2=is(iq)
            endif
            if(iq2.ne.0)then
              add(i,j)=add(iq2,1)
              rmsk(i,j)=rmsk(iq2,1)
              tsd(i,j)=tsd(iq2,1)
              almsk(i,j)=almsk(iq2,1)
              inumx(i,j)=1
            else
              numzer=numzer+1
            endif    !  (iq2.ne.0)
          endif      !  (inum(i,j).eq.0)
         enddo
        enddo
        do j=1,jl
         do i=1,il
          inum(i,j)=inumx(i,j)
         enddo
        enddo
       endif         ! (numzer.gt.0)
       if(numzer.gt.0)go to 2
      endif   !  (igrid.eq.0)

c add std.dev.to top.
      if ( osd ) then
        print *,'adding std.dev. to topog!!!!'
        do j=1,jl
          do i=1,il
            add(i,j) = add(i,j)+tsd(i,j)
          end do ! i=1,il
        end do ! j=1,jl
      endif !( osd ) then

c write out :
c   topography
c   land mask
c   sd

c***********************************************************************

      do j=1,jl
       do i=1,il
        zss(i,j) = add(i,j)
c       correct topog over ocean
        if (rmsk(i,j).lt..5.and.add(i,j).gt.-.2) zss(i,j)=-.1
       end do ! i=1,il
      end do ! j=1,jl

c     write out header to formatted file
!     write(luout,'(i3,i4,2f6.1,f5.2,f9.0,''  orog-mask-var'')')
      write(luout,'(i3,i4,2f6.1,f6.3,f8.0,''  orog-mask-var'')')
     .                           il,jl,rlong0,rlat0,schmidt,ds
c     write out g*zs(il,jl) to formatted file
      do j=1,jl
        do i=1,il
          zss(i,j)=9.80616*zss(i,j)
        end do ! i=1,il
      end do ! j=1,jl
      ilout=min(il,30)
      write (formout,'(''(''i3''f7.0)'')')ilout   !  i.e. (<il>f7.0)
      write(luout,formout) zss

c     write out land/sea mask to formatted file
      write (formout,'(''(''i3''f4.1)'')')ilout   !  i.e. (<il>f4.1)
      write(luout,formout) rmsk

c     write out std.dev of top. to formatted file
      write (formout,'(''(''i3''f6.0)'')')ilout   !  i.e. (<il>f6.0)
      write(luout,formout) tsd
      print *,'zss, rmsk and tsd written to unit',luout,' for il=',il

c     diagnostic printout of mask
c     work through i,j, values of the grid
      do j=1,jl
       do i=1,il
        iq=i+(j-1)*il
c       set lmask (here 1 for land, 0 for sea), just for diag printing
        inum(i,j)=0
        if(rmsk(i,j).ge.0.5)inum(i,j)=1
       enddo  ! i loop
      enddo   ! j loop
c     print out land-sea mask data for model grid by panel number
c     npanels=jl/il -1
      do n=npanels,0,-1
       if(npanels.gt.0)then
         ja=n*il +1
         jb=(n+1)*il
       endif
       print *,'land-sea mask for panel ',n
       do j=jb,ja,-1
        write(6,'(i4,1x,180i1)') j,(inum(i,j),i=1,il)
       enddo  !  j loop
      enddo   !  n loop

      end
      include 'setxyz.f'
      include 'jimcc.f'
      include 'latltoij.f'
c     include 'read_htt.f'  ! Jack's with filenames
