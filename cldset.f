      subroutine cldset(date,alat,cf,kltop,klbot)
c
c  this routine interpolates in latitude and time to set the cloud
c  amounts.
c  input
c    date    current model day of year
c    alat    latitude (from -pi/2 to pi/2)
c  output
c    cf      cloud fraction
c    kth     level of cloud top
c    kbh     level of cloud base
c
      parameter(pi=3.141592653589793,tpi=2.*pi,rlag=14.8125)
      parameter(year=365.25)

      real cf(5),kth(5),kbh(5)
c          winter       spring       summer       autumn       (nh)
      common /clddat/ ccd(37,5),ccd2(37,5),ccd3(37,5),ccd4(37,5),
     &                kkth(37,5),kkbh(37,5)
      external cldblk
      integer kltop(3),klbot(3)
      logical start
      data start / .true. /
      save start
c
      if ( start ) then
c rearrange the seasonal mean cloud data to allow interpolation
c define the amplitudes of the mean, annual and semi-annual cycles
        call resetd(ccd,ccd2,ccd3,ccd4,37*5)
        start = .false.
      end if

c     convert time to day number
      rang = tpi*(date-rlag)/year
      rsin1 = sin(rang)
      rcos1 = cos(rang)
      rcos2 = cos(2.0*rang)
c
      theta=90.-alat*180./pi
      ll = theta/5.
      angle = 5 * ll
      l = (theta +2.5)/5. + 1.0
      do 5 k = 2,4
        cf(k-1)=ccd(l,k)+rsin1*ccd2(l,k)+rcos1*ccd3(l,k)+rcos2*ccd4(l,k)
        kth(k)=kkth(l,k)
        kbh(k)=kkbh(l,k)
    5 continue
   
c     cl =cd(2)
c     cm =cd(3)
c     cl =cd(4)
   
c high
      kltop(1)=kth(2)/2-1
      klbot(1)=kth(2)/2-1
c mid
      kltop(2)=kth(3)/2-1
      klbot(2)=kth(3)/2-1
c low
      kltop(3)=kth(4)/2-1
      klbot(3)=kbh(4)/2-1
c
      return
      end
