      real rlnga(999),rlata(999),rlngb(999),rlatb(999)
c     Melb  145.00 -37.83, use 144.90, -37.70, 1.1
c     Syd   151.20 -33.88, use 150.97, -33.64, 1.1
c     Adel  138.63 -34.97, use 138.91, -34.79, 1.1
c     Brisb 153.02 -27.50, use 152.79, -27.51, 1.17
c     Towns 146.83 -19.30, use 146.72, -19.21,  .96
      print *,'this is findnear'
c     it uses fort.21 output from plotg.f, for mplot=4 option
c     it will produce linking points of the bounding boxes around
c     each grid point. For stretch C48 (.3), diff=~1.05 seems OK

c     it uses 4 links on each side of the box

      print *,'supply rlong, rlat, diff'
      read *,rlong,rlat,diff
      n=1
      nrepeat=0
2     read (21,*,end=8) rlnga(n),rlata(n),rlngb(n),rlatb(n)
      if(abs(rlnga(n)-rlong).lt.diff.and.abs(rlata(n)-rlat).lt.diff.and.
     .   abs(rlngb(n)-rlong).lt.diff.and.abs(rlatb(n)-rlat).lt.diff)then
c       print *,'rlnga,rlata,rlngb,rlatb ',
c    .           rlnga(n),rlata(n),rlngb(n),rlatb(n)
        nwrite=1  ! says OK to write
c       check for repeated links
        if(n.gt.1)then
          do nn=1,n-1
           if((rlnga(n).eq.rlngb(nn).and.rlata(n).eq.rlatb(nn).and.
     .         rlngb(n).eq.rlnga(nn).and.rlatb(n).eq.rlata(nn)).or.
     .        (rlnga(n).eq.rlnga(nn).and.rlata(n).eq.rlata(nn).and.
     .         rlngb(n).eq.rlngb(nn).and.rlatb(n).eq.rlatb(nn)))then
             nwrite=0
             nrepeat=nrepeat+1
           endif
          enddo
        endif  ! (n.gt.1)
        if(nwrite.eq.1)then
          n=n+1
c       else
c          following code needs fixing & moving to remove extraneous links
c          nsub=mod(n-1,4)
c          print *,'n,nsub ',n,nsub
c          n=n-nsub
        endif   ! (nwrite.eq.1)
      endif     ! (abs(rlnga .....)
      go to 2
8     print *,'repeated links found = ',nrepeat
      print *,'distinct nlinks = ',n-1
      do nn=1,n-1
       write(22,'(4f8.3)') rlnga(nn),rlata(nn),rlngb(nn),rlatb(nn)
      enddo
      end
