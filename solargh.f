c $Log$
c Revision 1.1  2003/08/13 01:24:20  dix043
c Initial revision
c
c Revision 1.1  1996/10/17  05:14:30  mrd
c Initial revision
c
c Revision 1.1  1994/03/23  09:45:40  mrd
c Initial revision
c
      subroutine solargh(day,bpyear,r,dlt,alp,slag,n,alat,
     *    hang,cosz,frac)
c
c  subroutine solar is called by march and computes the radius
c    vector, the declination and right ascension of the sun, the
c    apparent sun lag angle (related to the equation of time), the hour
c    angle of the sun at sunset, the mean cosine of the sun's zenith
c    angle and daylight fraction for n specified latitudes given the
c    day and fraction.  see notes of august 1989.
c
c    definition of the arguments used by this subroutine -
c
c      day   =        day of year
c      r     =        radius vector (distance to sun in a. u.)
c      dlt   =        declination of sun
c      alp   =        right ascension of sun
c      slag  =        apparent sun lag angle (west of mean sun is plus)
c      n     =        number of latitudes
c      alat  =        latitude values dimensioned by n
c      hang  =        hour angle of sun at sunrise or sunset
c      cosz  =        mean cosine of zenith angle
c      frac  =        fraction of daylight
c
c     note - all angles are expressed in radians
c            input variables are day in /ctlblk/, n and alat
c
c--------------
c   This routine comes from Bryant McAvaney via Josef Sytkus. It has an
c   additional argument not mentioned above
c      bpyear =       Years before present
c   This is used to set up the orbital parameters for paleoclimate runs.
c   Allowed values are 0, 6000 and 21000.
c
c   Martin Dix 23/3/94
c--------------
c
Cinclude "parm.h"
c
      dimension alat(n), hang(n)
      dimension cosz(n), frac(n)
c
c
Cinclude "ctlblk.h"
c
c  value of ccr depends on the precision of the machine
c
      data ccr / 1.0e-06 /
c
      pi     = 4.0 * atan(1.0)
      hpi    = 0.5 * pi
      tpi    = 2.0 * pi
      piph   = pi  + hpi
      radian = pi  / 180.0
c
c  obliquity in degrees is oblqty
c     data oblqty  / 23.446 /
c     data oblqty  / 24.105 /
c 6K      data oblqty  / 24.105 /
c 21K      data oblqty  / 22.949 /
c now      data oblqty  / 23.446 /
c
c  longitude of earth's perihelion in degrees is peril
c     data peril  / 102.04 /
c     data peril  /   0.87 /
c 6K      data peril  /   0.87 /
c 21K      data peril  / 114.42 /
c now      data peril  / 102.04 /
c
c  eccentricity of earth's orbit is ec
c     data ec  / 0.016724 /
c     data ec  /1.0e-36  /
c     data ec  / 0.018682 /
c now      data ec  / 0.016724 /
c 6K       data ec  / 0.018682 /
c 21K      data ec  / 0.018994 /
c circular data ec  / 1.0e-36 /
c
c   set fundamental astronomical parameters
      if(bpyear.eq.0.0)then
           ec  = 0.016724
           peril  = 102.04
           oblqty  = 23.446
      else if(bpyear.eq.6000.)then
           ec  = 0.018682
           peril  =   0.87
           oblqty  = 24.105
      else if(bpyear.eq.21000.)then
           ec  = 0.018994
           peril  = 114.42
           oblqty  = 22.949
      else
	 print *,'bad year for astronomy' ,bpyear
	 stop
      endif
      ecsq = ec**2
      angin  = radian * oblqty
      sni    = sin(angin)
c
c   tropical or calendar year in mean solar days is cyear
      cyear = 365.0
c
c   crate is mean daily angular motion of planet in degrees with
c   respect to equinox
      crate = 360.0 / cyear
c
c   semimajor axis of planetary orbit in astronomical units is sma.
c    one a. u. = 149,600,000 km.  for earth sma=1.
      sma = 1.0
c
c  The reference point for the astronomical calendar is the vernal
c  equinox, which is defined to occur at 0.0 hours on March 21. This
c  occurs 79 days after January 1, which is the reference point for
c  the model calendar.
       en = day - 79.0
c
c   longitude of sun at perihelion in degrees is sunper.  this is
c    exactly 180 degrees from that of earth at perihelion.
c
      sunper = peril - 180.0
      craten = crate * en
c
c  computation of longitude of mean sun follows the procedure
c  described in Berger (JAS, vol. 35, 2362-2367), using eq. 2
c  from annex 2 of Berger's paper "Precession, eccentricity, obliquity,
c  insolation and paleoclimates" from the NATO ASI "Long Term
c  Climatic Variations, Data and Modelling" (J. C. Duplessy, ed.).
c  The equation is a variant of that given on p.2365 of the JAS paper.
      cosphi = sqrt(1.0-ecsq)
      beta = (1.0-cosphi)/ec
      peri = radian*sunper
c
c  longitude of mean sun at vernal equinox is velngm
      velngm = -2.0*(beta*(1.0+cosphi)*sin(-peri)) +
     $          beta*beta*(0.5+cosphi)*sin(-2.0*peri) +
     $     beta*beta*beta*(1.0/3.0+cosphi)*sin(-3.0*peri)
c
c  mean longitude of sun in radians is el
      el = velngm + craten*radian
c
c  constrain el to range of zero to 2 * pi
      el = el - tpi*ifix(el/tpi)
      if (el .lt. 0.0) el= el+tpi
      elt   = el
c
c  ecliptic longitude of sun in radians is rlam
c  Berger (JAS, vol. 35, p. 2365)
      rlam = el + 2.0*ec*sin(el-peri)
     $             + 5.0*ecsq*sin(2.0*(el-peri))/4.0
     $             + 13.0*ec*ecsq*sin(3.0*(el-peri))/12.0
c
c   constrain rlam to range of zero to 2 * pi
      rlam = rlam - tpi*ifix(rlam/tpi)
      if (rlam .lt. 0.0) rlam = rlam + tpi
c
c   compute right ascension of sun (alp) in radians
      alp = atan(cos(angin)*tan(rlam))
c
c   following 6 statements determine proper quadrant for root of
c    equation for alp above and insure that the difference slag
c    below has an absolute magnitude of order pi/10.
c
      if (alp .lt. 0.0) alp = alp + tpi
c
      if ((elt-alp) .gt. piph) elt = elt - tpi
      if ((elt-alp) .lt.-piph) elt = elt + tpi
c
      if (abs(elt-alp) .gt. hpi) then
         if ((elt-alp) .gt. hpi) alp=alp+pi
         if ((elt-alp) .lt.-hpi) alp=alp-pi
      endif
c
c   compute declination of sun (dlt) in radians
      sind = sni*sin(rlam)
      dlt  = asin(sind)
c
c   compute radius vector (distance to sun in a.u.)
      r = (sma*(1.0-ecsq))/(1.0+ec*cos(rlam-peri))
c
c   compute sun lag angle in radians (720*slag/pi = equation
c    of time in minutes)
      slag = elt - alp
c
c  compute hour angle of sunset at all latitudes
c
      do 1 i=1,n
        ph = alat(i)
        if (dlt .ne. 0.0) then
          aap=abs(ph)
c
c  test whether latitude is near either pole
c
          eps = abs(aap-hpi)
          if (eps .le. ccr) then
c
c  hour angle of sunset at pole is either zero or pi depending on dlt
c
             h = hpi*abs(aap/ph+abs(dlt)/dlt)
            ss = sin(ph) * sin(dlt)
            cc = 0.0
          else
            ss = sin(ph) * sin(dlt)
            cc = cos(ph) * cos(dlt)
            ar = -ss / cc
            ac = abs(ar)
            if ((ac - 1.0 + ccr) .eq. 0.0) h = (ac - ar) * hpi
            if ((ac - 1.0 + ccr) .gt. 0.0) then
              if (ar .ge. 0.0) then
                h = 0.0
              else
                h = pi
              endif
            endif
            if ((ac -1.0 + ccr) .lt. 0.0) h = acos(ar)
          endif
        else
          h  = hpi
          ss = 0.0
          cc = cos(ph)
        endif
c
c  below is normal case computation of the sun's hour angle
c
        hang(i) = h
c
c  computation of cosz and frac added to solar for models without
c   diurnal variation since subroutine zenith need not be called.
c   cosz is constrained to be positive.
c
        frac(i) = h / pi
        if (h .ne. 0.0) then
          cosz(i) = amax1(ss + cc * sin(h) / h, 0.0)
        else
          cosz(i) = 0.0
        endif
    1 continue
c
      return
      end
