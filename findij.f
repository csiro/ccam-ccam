      include 'newmpar.h'
      character*20 cdum
      character*1 ,landin,icein,ltmp
      character*1 ,land1,land2,land3,land4,ice1,ice2,ice3,ice4
      print *,'supply rlong, rlat'
      read *,rlong,rlat
      read (22,*) cdum
      read (22,*) cdum
      distmin1=100.
      distmin2=100.
      distmin3=100.
      distmin4=100.
      do iq=1,ifull
       read (22,*) iqq,i,j,rlongg,rlatt,landin,icein,zsin
       dist=(rlong-rlongg)**2+(rlat-rlatt)**2
       if(dist.lt.distmin4)then
         distmin4=dist
         iqx4=iq
         ii4=i
         jj4=j
         rlongm4=rlongg
         rlatm4=rlatt
         land4=landin
         ice4=icein
	  zs4=zsin
         if(distmin4.lt.distmin3)then  !!!!!!
          tmp=distmin4
          distmin4=distmin3
	   distmin3=tmp
  	    itmp=iqx4
           iqx4=iqx3
	    iqx3=itmp
	   itmp=ii4
          ii4=ii3
	   ii3=itmp
	    itmp=jj4
           jj4=jj3
	    jj3=itmp
          tmp=rlongm4
          rlongm4=rlongm3
	   rlongm3=tmp
           tmp=rlatm4
           rlatm4=rlatm3
	    rlatm3=tmp
	   ltmp=land4
          land4=land3
	   land3=ltmp
    	    ltmp=ice4
           ice4=ice3
	    ice3=ltmp
          tmp=zs4
	   zs4=zs3
	   zs3=tmp
         endif                     !!!!!!!!
         if(distmin3.lt.distmin2)then  !!!!!!
          tmp=distmin3
          distmin3=distmin2
	   distmin2=tmp
  	    itmp=iqx3
           iqx3=iqx2
	    iqx2=itmp
	   itmp=ii3
          ii3=ii2
	   ii2=itmp
	    itmp=jj3
           jj3=jj2
	    jj2=itmp
          tmp=rlongm3
          rlongm3=rlongm2
	   rlongm2=tmp
           tmp=rlatm3
           rlatm3=rlatm2
	    rlatm2=tmp
	   ltmp=land3
          land3=land2
	   land2=ltmp
    	    ltmp=ice3
           ice3=ice2
	    ice2=ltmp
          tmp=zs3
	   zs3=zs2
	   zs2=tmp
         endif                     !!!!!!!!
         if(distmin2.lt.distmin1)then  !!!!!!
          tmp=distmin2
          distmin2=distmin1
	   distmin1=tmp
  	    itmp=iqx2
           iqx2=iqx1
	    iqx1=itmp
	   itmp=ii2
          ii2=ii1
	   ii1=itmp
	    itmp=jj2
           jj2=jj1
	    jj1=itmp
          tmp=rlongm2
          rlongm2=rlongm1
	   rlongm1=tmp
           tmp=rlatm2
           rlatm2=rlatm1
	    rlatm1=tmp
	   ltmp=land2
          land2=land1
	   land1=ltmp
    	    ltmp=ice2
           ice2=ice1
	    ice1=ltmp
          tmp=zs2
	   zs2=zs1
	   zs1=tmp
         endif                     !!!!!!!!
       endif
      enddo
      print *,'iqx1,ii1,jj1 ',iqx1,ii1,jj1
      print *,'rlongm,rlatm,dist ',rlongm1,rlatm1,distmin1
      print *,'land,ice,zs(in m) ',land1,ice1,zs1/9.806
      print *
      print *,'iqx2,ii2,jj2 ',iqx2,ii2,jj2
      print *,'rlongm,rlatm,dist ',rlongm2,rlatm2,distmin2
      print *,'land,ice,zs(in m) ',land2,ice2,zs2/9.806
      print *
      print *,'iqx3,ii3,jj3 ',iqx3,ii3,jj3
      print *,'rlongm,rlatm,dist ',rlongm3,rlatm3,distmin3
      print *,'land,ice,zs(in m) ',land3,ice3,zs3/9.806
      print *
      print *,'iqx4,ii4,jj4 ',iqx4,ii4,jj4
      print *,'rlongm,rlatm,dist ',rlongm4,rlatm4,distmin4
      print *,'land,ice,zs(in m) ',land4,ice4,zs4/9.806
      print *
      
      end
