! Version of o3set for global CC model. Based on GCM o3set.f, Revision 1.7 

! In this version the latitude may be different for each point
      subroutine o3set(alat,npts,mins,duo3n,sigma)
c
c  This routine interpolates in latitude and time to set the ozone 
c  amounts.
c  INPUT
c    ALAT    latitude (from -pi/2 to pi/2)
c    MINS    current model time in mins (from start of year)
c  OUTPUT
c    DUO3N   ozone mixing ratio
c
      include 'newmpar.h'
      include 'const_phys.h'
      parameter(rlag=14.8125)
      parameter(year=365)

      real duo3n(npts,kl), sigma(kl), alat(npts)
c          winter       spring       summer       autumn       (nh)
      common /o3dat/ dduo3n(37,kl),ddo3n2(37,kl),ddo3n3(37,kl),
     &     ddo3n4(37,kl)
      logical start
      data start / .true. /
      save start

c This moved to initfs
c
c     Convert time to day number
c     date = amod( float(mins)/1440., year)
c     Use year of exactly 365 days
      date = float(mod(mins,525600))/1440.
      rang = tpi*(date-rlag)/year
      rsin1 = sin(rang)
      rcos1 = cos(rang)
      rcos2 = cos(2.0*rang)
c
      do j=1,npts
         theta=90.-alat(j)*180./pi
         ilat = theta/5.
         angle = 5 * ilat
         than = (theta-angle)/5.
         ilat = ilat+1
         do m = 1,kl
            do3  = dduo3n(ilat,m) + rsin1*ddo3n2(ilat,m) 
     &              + rcos1*ddo3n3(ilat,m) + rcos2*ddo3n4(ilat,m)
            do3p = dduo3n(ilat+1,m) + rsin1*ddo3n2(ilat+1,m) 
     &              + rcos1*ddo3n3(ilat+1,m) + rcos2*ddo3n4(ilat+1,m)
            duo3n(j,m)=do3+than*(do3p-do3)
         end do
      end do
      return
      end
