c $Log$
c Revision 1.1  2003/08/13 01:24:20  dix043
c Initial revision
c
c Revision 1.1  1996/10/17  05:14:15  mrd
c Initial revision
c
c Revision 1.1  1991/02/22  16:37:19  ldr
c Initial release V3-0
c
      subroutine establ(es,tpg)

c     New version to run faster (?) Leon Rotstayn May 1990
C     MKS table
C     Table of es from -150 c to +70 c in one-degree increments is initi
c     in block data esbda.f
C
      dimension table(221)
      common /es_table/ table

      tp=min(max(tpg,123.16) , 343.16)
      tdiff1=tp-123.16
      it = int(tdiff1)
      tdiff=tdiff1-it
      es=(1.-tdiff)*table(it+1)+tdiff*table(it+2)
      return
      end
