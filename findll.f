!     this is findll.f
      include 'newmpar.h'
      character*20 cdum
      character*1 ,land,ice
      print *,'supply i, j'
      read *,i,j
      read (22,*) cdum
      read (22,*) cdum
      do iq=1,ifull
       read (22,*) iqq,ii,jj,rlong,rlat,land,ice,zs
       if(i.eq.ii.and.j.eq.jj)then
         print *,'iqq,ii,jj,rlong,rlat ',
     .            iqq,ii,jj,rlong,rlat
         print *,'land,ice,zs(in m) ',land,ice,zs/9.806
         stop
       endif
      enddo
      end
