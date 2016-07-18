      subroutine e1e288(g1,g2,g3,g4,g5,fxoe1,dte1)

c     subroutine e1e288 computes the exchange terms in the flux equation
c  for longwave radiation for all terms except the exchange with the
c  top of the atmosphere. the method is a table lookup on a pre-
c  computed e2 function (defined in ref. (4)).
c      the e1 function  calculations (formerly done in subroutine 
c  e1v88 compute the flux resulting from the exchange of photons
c  between a layer and the top of the atmosphere.  the method is a
c  table lookup on a pre-computed e1 function.
c     calculations are done in two frequency ranges:  
c       1) 0-560,1200-2200 cm-1   for q(approx) 
c       2) 160-560 cm-1           for q(approx,cts).
c  motivation for these calculations is in references (1) and (4).
c       inputs:                    (common blocks)
c     table1,table2,table3,em1,em1wde  tabcom 
c     avephi                           tfcom
c     temp                             radisw 
c     t                                kdacom 
c     fxoe1,dte1                argument list 
c       outputs:  
c     emiss                            tfcom
c     g1,g2,g3                  argument list,for 1st freq. range 
c     g4,g5                     argument list,for 2nd freq. range 
c 
c        called by :     fst88
c        calls     :  

CDIR$ TASK COMMON KDACOM
CDIR$ TASK COMMON RADISW
CDIR$ TASK COMMON TFCOM
c CDIR$ TASK COMMON VTEMP

      use kdacom_m
      use newmpar_m
      use radisw_m
      use tabcom_m
      use tfcom_m

      include 'hcon.h'
      include 'rdparm.h'
      integer it1(imax,ll3p),ival(imax,lp1)
      real dtot1(imax,ll3p),dtot2(imax,ll3p),dtot3(imax,ll3p),
     &  dtot4(imax,ll3p),dtot5(imax,llp1),dtot6(imax,llp1),
     &  dtot7(imax,llp1),dtot8(imax,llp1),
     &  f1(imax,lp1),f2(imax,lp1),f3(imax,lp1),
     &  fxo(imax,lp1),fyo(imax,lp1),dt(imax,lp1),du(imax,lp1),
     &  ww1(imax,lp1),ww2(imax,lp1),
     &  tmp3(imax,lp1),tmp5(imax),tmp9(imax)
c---variables equivalenced to common block variables
      dimension t1(5040),t2(5040),t4(5040)
      dimension em1v(5040),em1vw(5040)
c---variables in the argument list
      dimension fxoe1(imax,lp1),dte1(imax,lp1), 
     &    g1(imax,lp1),g2(imax,l),g3(imax,lp1),g4(imax,lp1),g5(imax,l) 
c
!      equivalence (em1v(1),em1(1,1)),(em1vw(1),em1wde(1,1)) 
!      equivalence (t1(1),table1(1,1)),(t2(1),table2(1,1)),
!     & (t4(1),table3(1,1))
     
      ! MJT resize array to replace equivalence
      em1v=reshape(em1, (/ size(em1) /))
      em1vw=reshape(em1wde, (/ size(em1wde) /))
      t1=reshape(table1, (/ size(table1) /))
      t2=reshape(table2, (/ size(table2) /))
      t4=reshape(table3, (/ size(table3) /))
     
c---first we obtain the emissivities as a function of temperature
c   (index fxo) and water amount (index fyo). this part of the code
c   thus generates the e2 function.
      do 121 k=1,lp1
c 
c---this loop replaces loops going fromi=1,imax and kp=1,lp1
      do 132 i=1,lp1i 
      fxo(i,1)=aint(t(i,1)*hp1)
      tmp3(i,1)=avephi(i,1,k)+h16e1
      dt(i,1)=t(i,1)-ten*fxo(i,1)
      fyo(i,1)=aint(tmp3(i,1)*ten)
      du(i,1)=tmp3(i,1)-hp1*fyo(i,1)
      fyo(i,1)=h28e1*fyo(i,1)
      ival(i,1)=fyo(i,1)+fxo(i,1) 
132   continue
c---decrementing the index by nine accounts for the tables begin-
c   ning at t=100k
      do 133 i=1,imax
      do 133 kp=1,lp1
      f1(i,kp)=t1(ival(i,kp)-9)
      f2(i,kp)=t2(ival(i,kp)-9)
      f3(i,kp)=t4(ival(i,kp)-9)
133   continue
      do 134 i=1,lp1i
      emiss(i,1,k)=f1(i,1)+du(i,1)*f2(i,1)+dt(i,1)*f3(i,1)
134   continue
c
c 
c***the following is the calculation for the e1 function, formerly
c    done in subroutine e1v88. the move to e1e288 is due to the 
c    savings in obtaining index values (the temp. indices have
c    been obtained in fst88, while the u-indices are obtained 
c    in the e2 calcs.,with k=1).
c 
      if (k.eq.1) then
c 
      do 209 i=1,imax*lp1
      it1(i,1)=fyo(i,1)+fxoe1(i,1)
      ww1(i,1)=ten-dte1(i,1)
      ww2(i,1)=hp1-du(i,1)
209   continue
      do 210 i=1,imax*l
      it1(i,lp2)=fyo(i,2)+fxoe1(i,1)
210   continue
      do 211 kp=1,lp1
      do 211 i=1,imax
      it1(i,kp+llp1)=fyo(i,kp)+fxoe1(i,1)
211   continue
c 
      do 225 i=1,imax
      do 226 kp=1,ll3p 
      dtot1(i,kp)=em1v(it1(i,kp)) 
      dtot2(i,kp)=em1v(it1(i,kp)+1) 
      dtot3(i,kp)=em1v(it1(i,kp)+28)
      dtot4(i,kp)=em1v(it1(i,kp)+29)
226   continue
c 
      do 248 kp=1,llp1 
      dtot5(i,kp)=em1vw(it1(i,kp))
      dtot6(i,kp)=em1vw(it1(i,kp)+1)
      dtot7(i,kp)=em1vw(it1(i,kp)+28) 
      dtot8(i,kp)=em1vw(it1(i,kp)+29) 
248   continue
225   continue
c 
      do 240 i=1,imax*lp1
      g1(i,1)=ww1(i,1)*ww2(i,1)*dtot1(i,1)+ 
     &        ww2(i,1)*dte1(i,1)*dtot2(i,1)+ 
     &        ww1(i,1)*du(i,1)*dtot3(i,1)+ 
     &        dte1(i,1)*du(i,1)*dtot4(i,1) 
      g4(i,1)=ww1(i,1)*ww2(i,1)*dtot5(i,1)+ 
     &        ww2(i,1)*dte1(i,1)*dtot6(i,1)+ 
     &        ww1(i,1)*du(i,1)*dtot7(i,1)+ 
     &        dte1(i,1)*du(i,1)*dtot8(i,1) 
240   continue
      do 241 kp=1,lp1
      do 241 i=1,imax
      g3(i,kp)=ww1(i,1)*ww2(i,kp)*dtot1(i,kp+llp1)+ 
     &        ww2(i,kp)*dte1(i,1)*dtot2(i,kp+llp1)+ 
     &        ww1(i,1)*du(i,kp)*dtot3(i,kp+llp1)+ 
     &        dte1(i,1)*du(i,kp)*dtot4(i,kp+llp1) 
241   continue
      do 242 i=1,imax*l
      g2(i,1)=ww1(i,1)*ww2(i,2)*dtot1(i,lp2)+ 
     &        ww2(i,2)*dte1(i,1)*dtot2(i,lp2)+ 
     &        ww1(i,1)*du(i,2)*dtot3(i,lp2)+ 
     &        dte1(i,1)*du(i,2)*dtot4(i,lp2) 
      g5(i,1)=ww1(i,1)*ww2(i,2)*dtot5(i,lp2)+ 
     &        ww2(i,2)*dte1(i,1)*dtot6(i,lp2)+ 
     &        ww1(i,1)*du(i,2)*dtot7(i,lp2)+ 
     &        dte1(i,1)*du(i,2)*dtot8(i,lp2) 
242   continue
c 
      endif 
c 
121   continue
c***we must now execute the special e2 case (i,k,k),when k=2,lp1 
      do 301 i=1,imax
      tmp5(i)=aint(temp(i,l)*hp1)
      tmp9(i)=temp(i,l)-ten*tmp5(i)
c----2) the case when k=lp1
      fxo(i,lp1)=aint(temp(i,lm1)*hp1)
      dt(i,lp1)=temp(i,lm1)-ten*fxo(i,lp1)
301   continue
c----1) the case when k=2..l
      do 302 k=2,l
      do 302 i=1,imax
      fxo(i,k)=tmp5(i)
      dt(i,k)=tmp9(i)
302   continue
      do 304 k=2,lp1
      do 304 i=1,imax
      tmp3(i,k)=avephi(i,k,k)+h16e1
304   continue
      do 305 i=1,imax*l
      fyo(i,2)=aint(tmp3(i,2)*ten)
      du(i,2)=tmp3(i,2)-hp1*fyo(i,2)
      ival(i,2)=h28e1*fyo(i,2)+fxo(i,2)
305   continue
c---the index is decremented by 9 to account for the tables
c   beginning at t=100k.
      do 306 i=1,imax
      do 306 k=2,lp1
      f1(i,k)=t1(ival(i,k)-9)
      f2(i,k)=t2(ival(i,k)-9)
      f3(i,k)=t4(ival(i,k)-9)
306   continue
      do 307 k=2,lp1
      do 307 i=1,imax
      emiss(i,k,k)=f1(i,k)+du(i,k)*f2(i,k)+dt(i,k)*f3(i,k)
307   continue
      return
      end 
