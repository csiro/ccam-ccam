c     subroutine e3v88 computes nearby layer transmissivities for 
c  h2o using a table lookup of the pre-computed e3 function 
c ( described in ref. (4)). 
c         inputs:                 (common blocks,args.) 
c       tv,av                      argument list
c       em3                        tabcom 
c          outputs: 
c       emv                        argument list
c 
c       called by  :    fst88 
c       calls      :    alog10h,alog10v 
c 
      subroutine e3v88(emv,tv,av) 
c 
c CDIR$ TASK COMMON VTEMP

      use tabcom_m
      use newmpar_m

      include 'hcon.h'
      include 'rdparm.h'
      integer it(imax,llp1)
      real fyo(imax,llp1),
     &     ww1(imax,llp1),
     &     tval(imax,llp1),dt(imax,llp1),ww2(imax,llp1),
     &     uval(imax,llp1),du(imax,llp1),
     &     fxo(imax,llp1),tmp3(imax,llp1)
      dimension emv(imax,llp1),tv(imax,llp1),av(imax,llp1)
      dimension em3v(5040)
      !equivalence (em3v(1),em3(1,1))
      
      ! MJT reshape arrays to replace equivalence statement
      em3v=reshape(em3, (/ size(em3) /))
      
c---the following loop replaces a double loop over i (1-imax) and
c   k (1-llp1)
      do 203 i=1,imax*llp1
        fxo(i,1)=aint(tv(i,1)*hp1)
        tmp3(i,1)=log10(av(i,1))+h16e1
        dt(i,1)=tv(i,1)-ten*fxo(i,1)
        fyo(i,1)=aint(tmp3(i,1)*ten)
        du(i,1)=tmp3(i,1)-hp1*fyo(i,1)
c---obtain index for table lookup; this value will have to be
c   decremented by 9 to account for table temps starting at 100K.
        it(i,1)=fxo(i,1)+fyo(i,1)*h28e1
        ww1(i,1)=ten-dt(i,1)
        ww2(i,1)=hp1-du(i,1)
        emv(i,1)=ww1(i,1)*ww2(i,1)*em3v(it(i,1)-9)+
     &           ww2(i,1)*dt(i,1)*em3v(it(i,1)-8)+
     &           ww1(i,1)*du(i,1)*em3v(it(i,1)+19)+
     &           dt(i,1)*du(i,1)*em3v(it(i,1)+20)
203   continue
      return
      end 
