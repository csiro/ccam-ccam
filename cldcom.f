c $Log$
c Revision 1.1  2003/08/13 01:24:20  dix043
c Initial revision
c
c Revision 1.1  1996/10/17  05:14:04  mrd
c Initial revision
c
c Revision 1.2  1994/12/07  23:53:36  mrd
c Explicitly declare all variables in include files.
c
c Revision 1.1  1992/04/15  11:13:35  mrd
c Initial revision
c
c      common block cldcom contains cloud transmission functions: 
c         cldfac     =  cloud transmission function,assuming random 
c                         overlap 
c 
      real cldfac(imax,lp1,lp1)
      common / cldcom / cldfac

