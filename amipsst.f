      subroutine amipsst
      use cc_mpi
      implicit none
!     this one primarily does namip=2      
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
      include 'mpif.h'
      real, dimension(ifull_g) ::  ssta_g, sstb_g, aice_g, bice_g
      real, dimension(ifull) :: ssta, sstb, aice, bice
      real fracice_in, dum2b
      common/work2b/fracice_in(ifull),dum2b(ifull,2)
      character*22 header
      integer, parameter, dimension(0:13) :: mdays =
     &     (/ 31, 31,28,31,30,31,30,31,31,30,31,30,31, 31 /)
      logical update
      integer iyr, imo, iday, iyr_m, imo_m, imonth, iyear, idjd_g,
     &        iyr_p, imo_p, ierr, iq
      save ssta,sstb,aice,bice,iyr,imo,iday,iyr_m,imo_m
      real rat1, rat2

      idjd_g = id + (jd-1)*il_g
      if ( myid == 0 ) then
      update = .false.

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
        update = .true.
        iyr_m=iyr
        imo_m=imo-1
        if(imo_m.eq.0)then
          imo_m=12
          iyr_m=iyr-1
        endif
        open(unit=75,file=sstfile,status='old',form='formatted')
2       read(75,'(i2,1x,i6,a22)') imonth,iyear,header
        if (myid==0) then
           print *,'reading sst data:',imonth,iyear,header
           print *,'comparing with imo_m,iyr_m ',imo_m,iyr_m
        end if
        read(75,*) ssta_g
        ssta_g(:)=ssta_g(:)*.01 -50. +273.16
        if (myid==0) then
           print *,'ssta(idjd) ',ssta_g(idjd_g)
        end if
        if(iyr_m.ne.iyear.or.imo_m.ne.imonth)go to 2

        read(75,'(i2,1x,i6,a22)') imonth,iyear,header
        if (myid==0) then
           print *,'reading sstb data:',imonth,iyear,header
           print *,'should agree with imo,iyr ',imo,iyr
        end if
        if(iyr.ne.iyear.or.imo.ne.imonth)stop
        read(75,*) sstb_g
        sstb_g(:)=sstb_g(:)*.01 -50. +273.16
        if (myid==0) then
           print *,'sstb(idjd) ',sstb_g(idjd_g)
        end if
	 
        open(unit=76,file=icefile,status='old',form='formatted')	 
	 if(namip.eq.2)then   ! sice also read at middle of month
21        read(76,'(i2,1x,i6,a22)') imonth,iyear,header
          if (myid==0) then
             print *,'reading a_sice data:',imonth,iyear,header
             print *,'comparing with imo_m,iyr_m ',imo_m,iyr_m
          end if
          read(76,*) aice_g
          if (myid==0) then
             print *,'aice(idjd) ',aice_g(idjd_g)
          end if
          if(iyr_m.ne.iyear.or.imo_m.ne.imonth)go to 21
          read(76,'(i2,1x,i6,a22)') imonth,iyear,header
          if (myid==0) then
             print *,'reading b_sice data:',imonth,iyear,header
             print *,'should agree with imo,iyr ',imo,iyr
          end if
          if(iyr.ne.iyear.or.imo.ne.imonth)stop
          read(76,*) bice_g
          if (myid==0) then
             print *,'bice(idjd) ',bice_g(idjd_g)
          end if
	 endif   ! (namip.eq.2) 
      endif     ! (ktau.eq.0)

      if(iday.eq.mdays(imo)/2)then
         update = .true.
         ssta_g(:)=sstb_g(:)    ! rotate data sets
         aice_g(:)=bice_g(:)    ! rotate data sets
         iyr_p=iyr
         imo_p=imo+1
         if(imo_p.eq.13)then
            imo_p=1
            iyr_p=iyr+1
         endif                  ! (imo_p.eq.13)
	 
!       read in next months data
3       read(75,'(i2,1x,i6,a22)') imonth,iyear,header
        if (myid==0) then
           print *,'reading sstb data:',imonth,iyear,header
           print *,'comparing with imo_p,iyr_p ',imo_p,iyr_p
        end if
        read(75,*) sstb_g
        sstb_g(:)=sstb_g(:)*.01 -50. +273.16
        if(iyr_p.ne.iyear.or.imo_p.ne.imonth)go to 3

        if(namip.eq.2)then   ! sice also read at middle of month
31        read(76,'(i2,1x,i6,a22)') imonth,iyear,header
          if (myid==0) then
             print *,'reading b_sice data:',imonth,iyear,header
             print *,'comparing with imo_p,iyr_p ',imo_p,iyr_p
          end if
          read(76,*) bice_g
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
      if (myid==0) then
         print *,'month_imo,iday',imo,iday
      end if
      if(iday.gt.1.or.namip.eq.2)go to 6

!     amipice for namip=1   N.B. used below to mask the SSTs
      print*, "NAMIP=1 not implemented in MPI version"
      stop
!!!4     read(76,'(i2,1x,i6,a22)') imonth,iyear,header
!!!      print *,'reading ice data:',imonth,iyear,header
!!!      print *,'with imo,iyr ',imo,iyr
!!!!     print *,'with imo,iyr,imo_p,iyr_p ',imo,iyr,imo_p,iyr_p
!!!      do j=1,jl_g
!!!         read(76,'(100i1)') (icemask(i),i=1,il_g)
!!!         do i=1,il_g
!!!            iq=i+(j-1)*il_g
!!!            sice(iq)=.false.
!!!            sicedep(iq)=0.
!!!            fracice(iq)=0.  
!!!            if(.not.land(iq).and.icemask(i).eq.1)then
!!!               sice(iq)=.true.
!!!               sicedep(iq)=.5
!!!               fracice(iq)=1.   ! as used in amip1 runs
!!!               ii=i
!!!               jj=j
!!!               iiq=iq
!!!            endif
!!!         enddo
!!!      enddo    
!!!      fracice(:)=1.            ! as used in amip1 runs
!!!      print *,'in amipsst; sice,sicedep,fracice: ',
!!!     .                     sice(idjd_g),sicedep(idjd_g),fracice(idjd_g)
!!!      print *,'in amipsst; ii,jj,sice,sicedep,fracice: ',ii,jj,
!!!     .                 sice(iiq),sicedep(iiq),fracice(iiq)
!!!      if(iyr.ne.iyear.or.imo.ne.imonth)go to 4

      end if ! myid==0

!     Each day interpolate non-land sst's
6     continue

!     Could calculate these in each process but this makes the code simpler
      call MPI_BCAST(rat1,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(rat2,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(update,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      if ( update ) then
         if ( myid == 0 ) then
            call ccmpi_distribute(ssta, ssta_g)
            call ccmpi_distribute(sstb, sstb_g)
            call ccmpi_distribute(aice, aice_g)
            call ccmpi_distribute(bice, bice_g)
         else
            call ccmpi_distribute(ssta)
            call ccmpi_distribute(sstb)
            call ccmpi_distribute(aice)
            call ccmpi_distribute(bice)
         end if
      end if

      do iq=1,ifull  
       if(.not.land(iq))then
         tgg(iq,1)=rat1*ssta(iq)+rat2*sstb(iq)  ! sea water temperature
c        alb(iq)=.11       ! set in indata/radriv90
       endif      ! (.not.land(iq))
      enddo
      if ( mydiag ) print *,'rat1,rat2,land,ssta,sstb,tgg1: ',
     %           rat1,rat2,land(idjd),ssta(idjd),sstb(idjd),tgg(idjd,1)
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
        if (mydiag) print *,'rat1,rat2,sice,aice,bice,fracice: ',
     &          rat1,rat2,sice(idjd),aice(idjd),bice(idjd),fracice(idjd)
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
      if (mydiag) print *,'rat1,rat2,sice,aice,bice,fracice: ',
     &        rat1,rat2,sice(idjd),aice(idjd),bice(idjd),fracice(idjd)
      return
      end
