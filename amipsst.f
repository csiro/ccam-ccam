      subroutine amipsst
      use cc_mpi
      use mlo, only : mloexport,sssb ! MJT mlo
      implicit none
!     this one primarily does namip>=2      
!     but persisted SST anomalies for namip=-1
!     A routine for the con-cubic model which will interpolate
!     linearly in time between two sst data sets.
!     iday is a variable which is equal to the number of
!     day of the month.(Thus iday = 1 at the
!     start of a month, is updated to 2 at end of first day etc.)
!     imo is the number of the month we are in (1<=imo<=12)
      include 'newmpar.h'
      include 'arrays.h'    ! ts, t, u, v, psl, ps, zs
      include 'dates.h'     !  kdate,ktime,timer,mtimer
      include 'filnames.h'  ! list of files
      include 'latlong.h'
      include 'nsibd.h'     ! res  for saving SST bias during month
      include 'parm.h'      ! id,jd
      include 'parmgeom.h'  ! rlong0,rlat0,schmidt
      include 'pbl.h'       ! tss
      include 'permsurf.h'  ! iperm etc
      include 'soil.h'      ! ,tice, alb
      include 'soilsnow.h'  ! fracice,sicedep
      include 'mpif.h'
      real, dimension(ifull_g) :: ssta_g, sstb_g, aice_g, bice_g
      real, dimension(ifull_g) :: sstc_g, cice_g
      real, save, dimension(ifull) :: ssta, sstb, sstc, aice, bice, cice
      real fraciceb, dum2b, x, c2, c3, c4
      common/work2b/fraciceb(ifull),dum2b(ifull,2)
      character*22 header
      integer, parameter, dimension(0:13) :: mdays =
     &     (/ 31, 31,28,31,30,31,30,31,31,30,31,30,31, 31 /)
      integer iyr, imo, iday, iyr_m, imo_m, imonth, iyear, idjd_g, iq
      save iyr,imo,iday,iyr_m,imo_m
      real rat1, rat2
      integer il_in,jl_in
      real rlon_in,rlat_in,schmidt_in
      integer, parameter :: mlomode = 1 ! (0=relax, 1=scale-select)    ! MJT mlo
      integer, parameter :: mlotime = 6 ! scale-select period in hours ! MJT mlo

      idjd_g = id + (jd-1)*il_g
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        iyr=kdate/10000
        imo=(kdate-10000*iyr)/100
        iday=kdate-10000*iyr-100*imo  +mtimer/(60*24)
        print *,'at start of amipsst for iyr,imo,iday ',iyr,imo,iday
        do while (iday>mdays(imo))
         iday=iday-mdays(imo)
         imo=imo+1
         if(imo>12)then
           imo=1               ! no leap years for the moment
           iyr=iyr+1
         endif
        enddo
        if(namip==-1)iyr=0

      if ( myid == 0 ) then 
        if(ktau==0)then
          iyr_m=iyr
          imo_m=imo-1
          if(imo_m==0)then
            imo_m=12
            iyr_m=iyr-1
            if(namip==-1)iyr_m=0
          endif
          open(unit=75,file=sstfile,status='old',form='formatted')
2         print *,'about to read amipsst file'
          read(75,*)
     &      imonth,iyear,il_in,jl_in,rlon_in,rlat_in,schmidt_in,header
          write(6,'("reading sst ",i2,i5,2i4,2f6.1,f7.4,1x,a22)')
     &      imonth,iyear,il_in,jl_in,rlon_in,rlat_in,schmidt_in,header
          if(il_g/=il_in.or.jl_g/=jl_in.or.rlong0/=rlon_in.or.
     &       rlat0/=rlat_in.or.schmidt/=schmidt_in)then
            write(0,*) 'il_g,il_in,jl_g,jl_in,rlong0,rlon_in',
     &                  il_g,il_in,jl_g,jl_in,rlong0,rlon_in
            write(0,*) 'rlat0,rlat_in,schmidt,schmidt_in',
     &                  rlat0,rlat_in,schmidt,schmidt_in
            write(0,*) 'wrong amipsst file'
            stop
          endif
          read(75,*) ssta_g
          ssta_g(:)=ssta_g(:)*.01 -50. +273.16
          print *,'want imo_m,iyr_m; ssta ',imo_m,iyr_m,ssta_g(idjd_g)
          if(iyr_m.ne.iyear.or.imo_m.ne.imonth)go to 2

          read(75,'(i2,i5,a22)') imonth,iyear,header
          print *,'reading sstb data:',imonth,iyear,header
          print *,'should agree with imo,iyr ',imo,iyr
          if(iyr.ne.iyear.or.imo.ne.imonth)stop
          read(75,*) sstb_g
          sstb_g(:)=sstb_g(:)*.01 -50. +273.16
          print *,'sstb(idjd) ',sstb_g(idjd_g)
c         extra read from Oct 08        
          read(75,'(i2,i5,a22)') imonth,iyear,header
          print *,'reading sstc data:',imonth,iyear,header
          read(75,*) sstc_g
          sstc_g(:)=sstc_g(:)*.01 -50. +273.16
          print *,'sstc(idjd) ',sstc_g(idjd_g)
  
          open(unit=76,file=icefile,status='old',form='formatted') 
          if(namip>=2)then   ! sice also read at middle of month
21         read(76,*)
     &      imonth,iyear,il_in,jl_in,rlon_in,rlat_in,schmidt_in,header
            write(6,'("reading ice ",i2,i5,2i4,2f6.1,f6.3,a22)')
     &        imonth,iyear,il_in,jl_in,rlon_in,rlat_in,schmidt_in,header
            if(il_g/=il_in.or.jl_g/=jl_in.or.rlong0/=rlon_in.or.
     &         rlat0/=rlat_in.or.schmidt/=schmidt_in)then
              write(0,*) 'il_g,il_in,jl_g,jl_in,rlong0,rlon_in',
     &                    il_g,il_in,jl_g,jl_in,rlong0,rlon_in
              write(0,*) 'rlat0,rlat_in,schmidt,schmidt_in',
     &                    rlat0,rlat_in,schmidt,schmidt_in
              write(0,*) 'wrong amipice file'
              stop
            endif
            read(76,*) aice_g
            print *,'want imo_m,iyr_m; aice ',imo_m,iyr_m,aice_g(idjd_g)
            if(iyr_m.ne.iyear.or.imo_m.ne.imonth)go to 21
            read(76,'(i2,i5,a22)') imonth,iyear,header
            print *,'reading b_sice data:',imonth,iyear,header
            print *,'should agree with imo,iyr ',imo,iyr
            if(iyr.ne.iyear.or.imo.ne.imonth)stop
            read(76,*) bice_g
            print *,'bice(idjd) ',bice_g(idjd_g)
c           extra cice read from Oct 08        
            read(76,'(i2,i5,a22)') imonth,iyear,header
            print *,'reading c_sice data:',imonth,iyear,header
            read(76,*) cice_g
            print *,'cice(idjd) ',cice_g(idjd_g)
           endif   ! (namip>=2) 
        endif     ! (ktau==0)
      endif ! myid==0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!     Each day interpolate-in-time non-land sst's
      if (ktau==0) then
         if ( myid == 0 ) then
            call ccmpi_distribute(ssta, ssta_g)
            call ccmpi_distribute(sstb, sstb_g)
            call ccmpi_distribute(sstc, sstc_g)
            call ccmpi_distribute(aice, aice_g)
            call ccmpi_distribute(bice, bice_g)
            call ccmpi_distribute(cice, cice_g)
         else
            call ccmpi_distribute(ssta)
            call ccmpi_distribute(sstb)
            call ccmpi_distribute(sstc)
            call ccmpi_distribute(aice)
            call ccmpi_distribute(bice)
            call ccmpi_distribute(cice)
         end if
         if(namip==-1)then
c          c1=0.
           do iq=1,ifull  
            if(.not.land(iq))then
              c2=ssta(iq)
              c3=ssta(iq)+sstb(iq)
              c4=c3+sstc(iq)     
              res(iq)=tgg(iq,1)-     ! just saves bias at beginning of run
     &           ( .5*c3+(4.*c3-5.*c2-c4)*x+1.5*(c4+3.*c2-3.*c3)*x*x )
            endif      ! (.not.land(iq))
           enddo
           print *,'some res values',(res(iq),iq=1,ifull,100)
         endif  ! (namip==-1)
      endif     ! (ktau==0)

      if(namip==-1)print *,'later_a ktau,res,tss ',
     &                             ktau,res(idjd),tss(idjd)
      x=(iday-1.)/mdays(imo)  ! simplest at end of day
      if(myid==0)print *,'month_imo,iday,x',imo,iday,x
      if(namip==2)then
        if(iday<mdays(imo)/2)then  ! 1st half of month
          rat1=(mdays(imo)-2.*iday)/(mdays(imo)+mdays(imo-1))
          rat2=(2.*iday+mdays(imo-1))/(mdays(imo)+mdays(imo-1))
          if(mydiag)print *,'rat1,rat2,land: ',
     &                       rat1,rat2,land(idjd)
          do iq=1,ifull  
           if(.not.land(iq))then
             tgg(iq,1)=rat1*ssta(iq)+rat2*sstb(iq)  ! sea water temperature
           endif      ! (.not.land(iq))
          enddo
          fraciceb(:)=min(.01*(rat1*aice(:)+rat2*bice(:)),1.) ! convert from %
        else                             ! 2nd half of month
          rat1=(mdays(imo+1)+2.*mdays(imo)-2.*iday)/
     .                         (mdays(imo+1)+mdays(imo))
          rat2=(2.*iday-mdays(imo))/(mdays(imo+1)+mdays(imo))
          do iq=1,ifull  
           if(.not.land(iq))then
             tgg(iq,1)=rat1*sstb(iq)+rat2*sstc(iq)  ! sea water temperature
           endif      ! (.not.land(iq))
          enddo
          fraciceb(:)=min(.01*(rat1*bice(:)+rat2*cice(:)),1.) ! convert from %
        endif
      endif  ! (namip==2)
      if(namip==-1)then
c      c1=0.
       do iq=1,ifull  
        if(.not.land(iq))then
          c2=ssta(iq)
          c3=ssta(iq)+sstb(iq)
          c4=c3+sstc(iq)     
          tgg(iq,1)=res(iq)+     ! adds SST bias back on to interp clim
     &          ( .5*c3+(4.*c3-5.*c2-c4)*x+1.5*(c4+3.*c2-3.*c3)*x*x )
          tss(iq)=tgg(iq,1)
        endif      ! (.not.land(iq))
       enddo
       if(namip==-1)print *,'later ktau,res,tss ',ktau,res(idjd),tss(iq)
       return
      endif  ! (namip==-1)
      if(namip>2)then
c       c1=0.
        do iq=1,ifull  
         if(.not.land(iq))then
           c2=ssta(iq)
           c3=ssta(iq)+sstb(iq)
           c4=c3+sstc(iq)          
           tgg(iq,1)=.5*c3+(4.*c3-5.*c2-c4)*x+1.5*(c4+3.*c2-3.*c3)*x*x
         endif      ! (.not.land(iq))
        enddo
      endif  ! (namip>2)
      if(namip==3)then
        do iq=1,ifull  
         c2=aice(iq)
         c3=aice(iq)+bice(iq)
         c4=c3+cice(iq)          
         fraciceb(iq)=min(.01*bice(iq),1.)
        enddo
      endif  ! (namip==3)
      if(namip==4)then
        do iq=1,ifull  
           c2=aice(iq)
           c3=aice(iq)+bice(iq)
           c4=c3+cice(iq)          
           fraciceb(iq)=min(.01*(.5*c3+(4.*c3-5.*c2-c4)*x
     &                     +1.5*(c4+3.*c2-3.*c3)*x*x),1.)
        enddo
      endif  ! (namip==4)
      if(mydiag)print *,'ktau,ssta,sstb,sstc,tgg1: ',
     &                 ktau,ssta(idjd),sstb(idjd),sstc(idjd),tgg(idjd,1)
      do iq=1,ifull
        if(fraciceb(iq)<=.02)fraciceb(iq)=0.
      enddo
      
      if (nmlo.eq.0) then ! MJT mlo
       sicedep(:)=0. 
       if(ktau==0)then  ! will set sicedep in indata
        fracice(:)=fraciceb(:)
        do iq=1,ifull
         if(.not.land(iq))then
           tss(iq)=tgg(iq,1)
         endif    ! (.not.land(iq))
        enddo
        return
       endif       ! (ktau==0)
      
       do iq=1,ifull
        if(.not.land(iq))then
         if(fraciceb(iq)>0.)then
           if(fracice(iq)==0.)then
!            create values for tice, and set averaged tss
!            N.B. if already a sice point, keep present tice
             tggsn(iq,1)=min(271.2,tss(iq),t(iq,1)+.04*6.5) ! for 40 m lev1 ! MJT seaice
           endif  ! (fracice(iq)==0.)
           if(rlatt(iq)>0.)then
             sicedep(iq)=2.
           else
             sicedep(iq)=1.
           endif ! (rlatt(iq)>0.)
         endif    ! (fraciceb(iq)>0.)
         fracice(iq)=fraciceb(iq)
         tss(iq)=tggsn(iq,1)*fracice(iq)+tgg(iq,1)*(1.-fracice(iq)) ! MJT seaice
        endif      ! (.not.land(iq))
       enddo
      !--------------------------------------------------------------
      ! MJT mlo
      else
        ! if (.not.allocated(sssb)) then
        !   allocate(sssb(ifull))
        ! end if
        ! sssb=input salinity
        select case(mlomode)
          case(0) ! relax
            call mlonudge(tgg(:,1),3)
          case(1)
            if (mod(mtimer,mlotime*60).eq.0) then
              call mlofilterfast(tgg(:,1),3) ! 1D version
              !call mlofilter(tgg(:,1),3) ! 2D version
            end if
          case DEFAULT
            write(6,*) "ERROR: Unknown mlomode ",mlomode
            stop
        end select
        call mloexport(0,tgg(:,1),1,0)        
      end if
      !--------------------------------------------------------------
      if (mydiag)print *,'ktau,sicedep,aice,bice,cice,fracice: ',
     & ktau,sicedep(idjd),aice(idjd),bice(idjd),cice(idjd),fracice(idjd)
      return
      end
