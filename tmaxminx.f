!     routine to print sequence of daily min & max values
!     this one assumes run starts from 9 pm    adds 10 to lu
      nperday=72  ! number of timesteps per day
      ndays=8     ! number of days on file
      nadd=10
      print *,'Melbourne temperatures'
      call maxmin(60+nadd,ndays,nperday)
      print *
      print *,'Sydney temperatures'
      call maxmin(61+nadd,ndays,nperday)
      print *
      print *,'Adelaide temperatures'
      call maxmin(62+nadd,ndays,nperday)
      print *
      print *,'Canberra temperatures'
      call maxmin(63+nadd,ndays,nperday)
      print *
      print *,'Hobart temperatures'
      call maxmin(64+nadd,ndays,nperday)
      print *
      print *,'Brisbane temperatures'
      call maxmin(65+nadd,ndays,nperday)
      print *
      print *,'Perth temperatures'
      call maxmin(66+nadd,ndays,nperday)
      print *
      print *,'Darwin temperatures'
      call maxmin(67+nadd,ndays,nperday)
      print *
      print *,'Alice Springs temperatures'
      call maxmin(68+nadd,ndays,nperday)
      print *
      print *,'Albury temperatures'
      call maxmin(69+nadd,ndays,nperday)
      end
      subroutine maxmin(lu,ndays,nperday)
      do nd=1,ndays
       tmin=600.
       tmax=0.
       do nmin=1,nperday/2
        read (lu,*) ktau,t,preca
        tmin=min(tmin,t)
       enddo
       do nmax=1,nperday/2
        read (lu,*) ktau,t,precb
        tmax=max(tmax,t)
       enddo
      print 9,nd,preca,precb,tmin-273.16,tmax-273.16
9     format(' day,tmin,tmax: ',i3,4f6.1)
      enddo
c     print *,'final ktau = ',ktau
      return
      end

