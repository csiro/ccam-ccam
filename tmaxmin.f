!     routine to print sequence of daily min & max values
!     this one assumes run starts from 9 pm
!        (no:   and averages adjacent timesteps    from 25/2/99)
      common/dates/kdate,ktime
      character header*79
      nperday=72  ! number of timesteps per day
      ndays=8     ! number of days on file
      ndays=1     ! number of days on file
c     read(60,'(1x,a79)') header
      read(60,*) header,kdate,ktime
      print *,'Starting anal: ',kdate,ktime
      rewind 60
      print *,'Melbourne temperatures'
      call maxmin(60,ndays,nperday)
      print *
      print *,'Sydney temperatures'
      call maxmin(61,ndays,nperday)
      print *
      print *,'Adelaide temperatures'
      call maxmin(62,ndays,nperday)
      print *
      print *,'Canberra temperatures'
      call maxmin(63,ndays,nperday)
      print *
      print *,'Hobart temperatures'
      call maxmin(64,ndays,nperday)
      print *
      print *,'Brisbane temperatures'
      call maxmin(65,ndays,nperday)
      print *
      print *,'Perth temperatures'
      call maxmin(66,ndays,nperday)
      print *
      print *,'Darwin temperatures'
      call maxmin(67,ndays,nperday)
      print *
      print *,'Alice Springs temperatures'
      call maxmin(68,ndays,nperday)
      print *
      print *,'Albury temperatures'
      call maxmin(69,ndays,nperday)
      print *,'above are: day rnd12 rnd24 tmin tmax cll clm clh clt u v'
      end
      subroutine maxmin(lu,ndays,nperday)
      common/dates/kdate,ktime
      real dum7(7)
      character header*79
      read(lu,'(1x,a79)') header
c     if(lu.eq.60)print *,'starting anal: ',header
      do nd=1,ndays
       tmin=600.
       tmax=0.
       do nmin=1,nperday/2
        read (lu,*) ktau,t,preca,dum,dum7,cll,clm,clh,cld,dum7,u,v
        if(ktime.eq.1200)then
          tmin=min(tmin,t)
        else
	    if(t.gt.tmax)then
	      tmax=t
	      ux=u
	      vx=v
	      cllx=cll
	      clmx=clm-1.
	      clhx=clh-2.
	      cldx=cld-3.     
	    endif
        endif
c       if(nmin.gt.1)tmin=min(tmin,.5*(tprev+t))
c       tprev=t
       enddo
       do nmax=1,nperday/2
        read (lu,*) ktau,t,precb,dum,dum7,cll,clm,clh,cld,dum7,u,v
c       print *,    ktau,t,precb,dum,dum7,cll,clm,clh,cld,dum7,u,v
        if(ktime.eq.1200)then
	    if(t.gt.tmax)then
	      tmax=t
	      ux=u
	      vx=v
	      cllx=cll
	      clmx=clm-1.
	      clhx=clh-2.
	      cldx=cld-3.     
	    endif
        else
          tmin=min(tmin,t)
        endif
c       if(nmax.gt.1)tmax=max(tmax,.5*(tprev+t))
c       tprev=t
       enddo
!     print 9,nd,preca,precb,tmin-273.16,tmax-273.16
      print 9,nd,preca,precb,tmin,tmax,cllx,clmx,clhx,cldx,ux,vx
9     format(' day,tmin,tmax: ',i3,4f6.1,4f6.2,2f6.1)
      enddo
c     print *,'final ktau = ',ktau
      return
      end

