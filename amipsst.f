      subroutine amipsst
!     this one primarily does namip=2      
!     Has entry amipice
!     A routine for the con-cubic model which will interpolate
!     linearly in time between two sst data sets.
!     iday is a variable which is equal to the number of
!     day of the month.(Thus iday = 1 at the
!     start of a month, is updated to 2 at end of first day etc.)
!     imo is the number of the month we are in (1<=imo<=12)
      include 'newmpar.h'
      include 'arrays.h'   ! ts, t, u, v, psl, ps
      include 'dates.h'     !  kdate,ktime,timer,mtimer
      include 'filnames.h'  ! list of files
      include 'map.h'       ! zs land
      include 'nsibd.h'     ! sib, tsigmf, isoilm
      include 'parm.h'      ! id,jd
      include 'pbl.h'       ! tss
      include 'soil.h'      ! sice,sicedep,tice, fracice alb
      include 'soilsnow.h' ! new soil arrays for scam - tgg too
      real ssta(ifull),sstb(ifull),aice(ifull),bice(ifull)
      common/work2b/fracice_in(ifull),dum2b(ifull,2)
      dimension icemask(il)
      character*22 header
      integer mdays(0:13)
      data mdays/31, 31,28,31,30,31,30,31,31,30,31,30,31, 31/
      save ssta,sstb,aice,bice,iyr,imo,iday,iyr_m,imo_m

      print *,'amipsst  called'
      iyr=kdate/10000
      imo=(kdate-10000*iyr)/100
      iday=kdate-10000*iyr-100*imo  +mtimer/(60*24)
      do while (iday.gt.mdays(imo))
       iday=iday-mdays(imo)
       imo=imo+1
       if(imo.gt.12)then
         imo=1               ! no leap years for the moment
         iyr=iyr+1
       endif
      enddo

      if(ktau.eq.0)then
!       ***  this code may need generalizing for multi-month runs
        iyr_m=iyr
        imo_m=imo-1
        if(imo_m.eq.0)then
          imo_m=12
          iyr_m=iyr-1
        endif
        open(unit=75,file=sstfile,status='old',form='formatted')
2       read(75,'(i2,1x,i6,a22)') imonth,iyear,header
        print *,'reading sst data:',imonth,iyear,header
        print *,'comparing with imo_m,iyr_m ',imo_m,iyr_m
        read(75,*) ssta
        do iq=1,ifull
         ssta(iq)=ssta(iq)*.01 -50. +273.16
        enddo
        print *,'ssta(idjd) ',ssta(idjd)
        if(iyr_m.ne.iyear.or.imo_m.ne.imonth)go to 2

        read(75,'(i2,1x,i6,a22)') imonth,iyear,header
        print *,'reading sstb data:',imonth,iyear,header
        print *,'should agree with imo,iyr ',imo,iyr
        if(iyr.ne.iyear.or.imo.ne.imonth)stop
        read(75,*) sstb
        do iq=1,ifull
         sstb(iq)=sstb(iq)*.01 -50. +273.16
        enddo
        print *,'sstb(idjd) ',sstb(idjd)
	 
        open(unit=76,file=icefile,status='old',form='formatted')	 
	 if(namip.eq.2)then   ! sice also read at middle of month
21        read(76,'(i2,1x,i6,a22)') imonth,iyear,header
          print *,'reading a_sice data:',imonth,iyear,header
          print *,'comparing with imo_m,iyr_m ',imo_m,iyr_m
          read(76,*) aice
          print *,'aice(idjd) ',aice(idjd)
          if(iyr_m.ne.iyear.or.imo_m.ne.imonth)go to 21
          read(76,'(i2,1x,i6,a22)') imonth,iyear,header
          print *,'reading b_sice data:',imonth,iyear,header
          print *,'should agree with imo,iyr ',imo,iyr
          if(iyr.ne.iyear.or.imo.ne.imonth)stop
          read(76,*) bice
          print *,'bice(idjd) ',bice(idjd)
	 endif   ! (namip.eq.2) 
      endif     ! (ktau.eq.0)

      if(iday.eq.mdays(imo)/2)then
        do iq=1,ifull
         ssta(iq)=sstb(iq)   ! rotate data sets
         aice(iq)=bice(iq)   ! rotate data sets
        enddo
        iyr_p=iyr
        imo_p=imo+1
        if(imo_p.eq.13)then
          imo_p=1
          iyr_p=iyr+1
        endif  ! (imo_p.eq.13)
	 
!       read in next months data
3       read(75,'(i2,1x,i6,a22)') imonth,iyear,header
        print *,'reading sstb data:',imonth,iyear,header
        print *,'comparing with imo_p,iyr_p ',imo_p,iyr_p
        read(75,*) sstb
        do iq=1,ifull
         sstb(iq)=sstb(iq)*.01 -50. +273.16
        enddo
        if(iyr_p.ne.iyear.or.imo_p.ne.imonth)go to 3

        if(namip.eq.2)then   ! sice also read at middle of month
31        read(76,'(i2,1x,i6,a22)') imonth,iyear,header
          print *,'reading b_sice data:',imonth,iyear,header
          print *,'comparing with imo_p,iyr_p ',imo_p,iyr_p
          read(76,*) bice
          if(iyr_p.ne.iyear.or.imo_p.ne.imonth)go to 31
        endif  ! (namip.eq.2)
      endif    ! (iday.eq.mdays(imo)/2)

      if(iday.lt.mdays(imo)/2)then  ! 1st half of month
        rat1=(mdays(imo)-2.*iday)/(mdays(imo)+mdays(imo-1))
        rat2=(2.*iday+mdays(imo-1))/(mdays(imo)+mdays(imo-1))
      else                             ! 2nd half of month
        rat1=(mdays(imo+1)+2.*mdays(imo)-2.*iday)/
     .                               (mdays(imo+1)+mdays(imo))
        rat2=(2.*iday-mdays(imo))/(mdays(imo+1)+mdays(imo))
      endif
      print *,'month_imo,iday',imo,iday
      if(iday.gt.1.or.namip.eq.2)go to 6

!     amipice for namip=1   N.B. used below to mask the SSTs
4     read(76,'(i2,1x,i6,a22)') imonth,iyear,header
      print *,'reading ice data:',imonth,iyear,header
      print *,'with imo,iyr ',imo,iyr
!     print *,'with imo,iyr,imo_p,iyr_p ',imo,iyr,imo_p,iyr_p
      do j=1,jl
       read(76,'(100i1)') (icemask(i),i=1,il)
       do i=1,il
	 iq=i+(j-1)*il
        sice(iq)=.false.
        sicedep(iq)=0.
        fracice(iq)=0.  
        if(.not.land(iq).and.icemask(i).eq.1)then
          sice(iq)=.true.
          sicedep(iq)=.5
          fracice(iq)=1.  ! as used in amip1 runs
	   ii=i
	   jj=j
	   iiq=iq
        endif
       enddo
      enddo    
      do iq=1,ifull
       fracice(iq)=1.  ! as used in amip1 runs
      enddo  
      print *,'in amipsst; sice,sicedep,fracice: ',
     .                     sice(idjd),sicedep(idjd),fracice(idjd)
      print *,'in amipsst; ii,jj,sice,sicedep,fracice: ',ii,jj,
     .                 sice(iiq),sicedep(iiq),fracice(iiq)
      if(iyr.ne.iyear.or.imo.ne.imonth)go to 4

!     Each day interpolate non-land sst's
6     do iq=1,ifull  
       if(.not.land(iq))then
         tgg(iq,1)=rat1*ssta(iq)+rat2*sstb(iq)  ! sea water temperature
c        alb(iq)=.11       ! set in indata/radriv90
       endif      ! (.not.land(iq))
      enddo
      print *,'rat1,rat2,land,ssta,sstb,tgg1: ',
     .           rat1,rat2,land(idjd),ssta(idjd),sstb(idjd),tgg(idjd,1)
      fracice_in(:)=min(.01*(rat1*aice(:)+rat2*bice(:)),1.) ! convert from %
      do iq=1,ifull
        if(fracice_in(iq).le..02)fracice_in(iq)=0.
      enddo
      
      if(ktau.eq.0)then
        fracice(:)=fracice_in(:)
        do iq=1,ifull
         if(.not.land(iq))then
           tss(iq)=tgg(iq,1)
	    if(fracice(iq).gt.0.)then
             sice(iq)=.true.
             sicedep(iq)=2.  ! N.B. tss, tgg3 re-set in indata
	    endif  ! (fracice(iq).gt.0.)
         endif    ! (.not.land(iq))
        enddo
        print *,'rat1,rat2,sice,aice,bice,fracice: ',
     .          rat1,rat2,sice(idjd),aice(idjd),bice(idjd),fracice(idjd)
        return
      endif       ! (ktau.eq.0)
      
      do iq=1,ifull
       if(.not.land(iq))then
	  if(fracice_in(iq).gt.0.)then
	    if(fracice(iq).eq.0.)then
!            create values for tice, and set averaged tss
!            N.B. if already a sice point, keep present tice
!            N.B. sflux will update changed sice, sicedep
             tgg(iq,3)=min(271.2,tss(iq),t(iq,1)+.04*6.5) ! for 40 m level 1
           endif  ! (fracice(iq).eq.0.)
	  endif    ! (fracice_in(iq).gt.0.)
         fracice(iq)=fracice_in(iq)
         tss(iq)=tgg(iq,3)*fracice(iq)+tgg(iq,1)*(1.-fracice(iq))
       endif      ! (.not.land(iq))
      enddo
      print *,'rat1,rat2,sice,aice,bice,fracice: ',
     .        rat1,rat2,sice(idjd),aice(idjd),bice(idjd),fracice(idjd)
      return
      end
