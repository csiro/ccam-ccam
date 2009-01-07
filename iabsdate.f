      function iabsdate(kdate_r,kdate)
      
      common/leap_yr/leap
      
!     Y2K version
c     if before 1960, suppress leap years (GCM-nested runs) (jlm 31/5/95)
c     N.B. All GCM dates to have year denoted by <1960      (15/12/98)
c     Little function to convert kdate_r in form YYYYMMDD to no. of   !  Y2K
c     days since start of year of kdate. (LDR 3/1992, jlm 15/12/98,15/12/00)
c     Accounts for leap years. N.B. no exception in year 2000

c     N.B. this latest version is designed for nested runs running 1 month
c     at a time, so iyear0=iyear except for the case 
c     when kdate is Dec. and kdate_r is the following Jan.

      integer mdays(12)
      data mdays/31,28,31,30,31,30,31,31,30,31,30,31/

      iabsdate=0
      iyear=kdate_r/10000
      iyear0=kdate/10000                ! year of kdate
      month=(kdate_r-10000*iyear)/100
      iday=(kdate_r-10000*iyear)-100*month

!     calculate number of months since start of kdate year
      months=(iyear-iyear0)*12+month-1  
!     print *,'in iabsdate  kdate_r= ',kdate_r
!     print *,'iyear,month,iday,months ',iyear,month,iday,months

c     Accumulate days month by month, up to last completed month

      do mon=1,months
        mnth=mod(mon-1,12)+1
        iabsdate=iabsdate+mdays(mnth)
	if (mnth.eq.2.and.leap.eq.1) then
          nl=0
          if (mod(iyear0,4  ).eq.0) nl=1
          if (mod(iyear0,100).eq.0) nl=0
          if (mod(iyear0,400).eq.0) nl=1
	  iabsdate=iabsdate+nl
	end if
!        if(mnth.eq.2.and.mod(iyear0,4).eq.0.and.iyear0.ge.1960)
!     .                  iabsdate=iabsdate+1 ! Leap year
      enddo

c     Add days from this current month

      iabsdate=iabsdate+iday
!     print *,'iyr,mnth,iabsdate ',iyr,mnth,iabsdate

      end
