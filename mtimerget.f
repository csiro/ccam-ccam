      subroutine mtimerget(mtimer,kdate1,ktime1,kdate2,ktime2)
!     returns mtimer in minutes, corr. to (kdate2,ktime2) -  (kdate1,ktime1)    
      dimension ndoy(12)   ! days from beginning of year (1st Jan is 0)
      data ndoy/ 0,31,59,90,120,151,181,212,243,273,304,334/
      common/leap_yr/leap  ! 1 to allow leap years
 
      if(leap.ne.0)stop 'leap years not catered for in mtimerget'
!     Set up number of minutes from beginning of year
!     For GCM runs assume year is <1980 (e.g. ~321-460 for 140 year run)
      jyear1=kdate1/10000
      jmonth=(kdate1-jyear1*10000)/100
      jday=kdate1-jyear1*10000-jmonth*100
      jhour=ktime1/100
      jmin=ktime1-jhour*100
      mstart1=1440*(ndoy(jmonth)+jday-1) + 60*jhour + jmin ! mins from start of y

      jyear2=kdate2/10000
      jmonth=(kdate2-jyear2*10000)/100
      jday=kdate2-jyear2*10000-jmonth*100
      jhour=ktime2/100
      jmin=ktime2-jhour*100
      mstart2=1440*(ndoy(jmonth)+jday-1) + 60*jhour + jmin ! mins from start of y

      mtimer=mstart2-mstart1+(jyear2-jyear1)*365*24*60
      return
      end
