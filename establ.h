c $Log$
c Revision 1.1  2003/08/13 01:24:20  dix043
c Initial revision
c
c Revision 1.7  2000/11/14 03:11:35  rot032
c Energy/water balance for land-surface and sea-ice, plus tidy ups from HBG
c
c Revision 1.6  1998/12/10 00:55:29  ldr
c HBG changes to V5-1-21
c
c Revision 1.5  1996/11/28  00:59:35  mrd
c Add type declaration of pp.
c
c Revision 1.4  1996/03/21  03:18:28  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.3  1995/06/30  02:44:39  ldr
c Changes to qcloud (run F88) and put qsat formula into ESTABL.f
c
c Revision 1.2  1994/12/07  23:53:36  mrd
c Explicitly declare all variables in include files.
c
c Revision 1.1  1991/03/14  09:49:41  ldr
c Initial revision
c
      real table, tablei
      common /es_table/ table(0:220)
      common /esitable/ tablei(0:220)

c Arithmetic statement functions to replace call to establ.
c T is temp in Kelvin, which should lie between 123.16 and 343.16;
c TDIFF is difference between T and 123.16, subject to 0 <= TDIFF <= 220

      real t, tdiff, establ, qsat, pp, estabi, qsati
      tdiff(t)=min(max( t-123.16, 0.), 219.)

      establ(t) = (1.-(tdiff(t)-aint(tdiff(t))))*table(int(tdiff(t)))
     &           + (tdiff(t)-aint(tdiff(t)))*table(int(tdiff(t))+1)
      qsat(pp,t) = epsil*establ(t)/(pp-establ(t)) !Usual formula
c      qsat(pp,t) = epsil*establ(t)/pp !Consistent with V4-5 to V4-7

c These give the ice values needed for the qcloud scheme
      estabi(t) = (1.-(tdiff(t)-aint(tdiff(t))))*tablei(int(tdiff(t)))
     &           + (tdiff(t)-aint(tdiff(t)))*tablei(int(tdiff(t))+1)
      qsati(pp,t) = epsil*estabi(t)/(pp-estabi(t)) !Usual formula
c      qsati(pp,t) = epsil*estabi(t)/pp !Consistent with V4-5 to V4-7
