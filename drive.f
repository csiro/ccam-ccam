c $Log$
c Revision 1.1  2003/08/13 01:24:20  dix043
c Initial revision
c
c Revision 1.1  1996/10/17  05:14:14  mrd
c Initial revision
c
c Revision 1.1  91/07/19  11:40:47  ldr
c Initial revision
c 
      subroutine drive
c     ******************************************************************
c          subroutine drive initiates the fels fast longwave radiation
c     program . it performs the following functions
c     1) calculates co2 transmissions (the result is the array co21)
c     using routine co2trn
c     2) definition of the flux levels of temperature and pressure,
c     and obtaining of optical thicknesses between the flux levels for
c     h2o and o3 (var1,var2,var3,var4).
c     3) calling the main computation module ,subroutine fstrd.
c     ******************************************************************
      parameter (kmax=9,l=kmax,lp1=l+1)
      parameter (lm1=l-1,ll=2*l,llp1=ll+1,llp2=ll+2,llp3=ll+3)
      parameter (llm1=ll-1,llm2=ll-2,lp3=l+3)
      parameter (llzz=ll,lp1m=360,lp1m1=lp1m-1)
      parameter (lp2=l+2,lm2=l-2,kx=l,kmx2p2=ll+2)
      parameter (kmx2p3=ll+3,kmx2p4=ll+4,kpad=2002,kdm=l)
      parameter (nb=l,nb1=nb+1,kd2=20,kd1p=11,kd2m=kd2-1)
      parameter (nl=l,nl1=nl+1,nl2=nl+2,nl3=nl+3)
c     implicit half precision (a-h,o-z)
      common/input/temp(lp1),press(lp1),r(l),qo3(l),ch,cm,cl,
     1 ich,icm,ict,icb
      common/kdata/pm(lp1,lp1),p1(lp1,lp1),psk(lp1,lp1),co21(lp1,lp1),
     * q(lp1),p(llp3),delp(l),t(lp1),var1(l),var2(l),var3(l),var4(l),
     * grav,grav1,tmp4
      common/outpfs/heatra(l),grnflx,z3
      dimension qqo3(l)
      common/mont/ipoint, jrow
      dimension mask(lp1),rmask(lp1)
      equivalence (mask(1),rmask(1))
c
      save
      do 100 i=1,lp1
          if( temp(i).ge.350. ) go to 115
 100  continue
      go to 200
c
 115  continue
      print 120,ipoint,jrow, (temp(i),i=1,lp1)
  120 format('0 temperature to large for table look-up'/
     .       '   in long-wave radiation code'/
     .       ' at point',i5,'  in row',i5/(5x,8e13.5))
      do 130 i=1,lp1
      temp(i)=min(temp(i),349.999e0)
  130 continue
  200 continue
c
      pstar=press(lp1)
c     print *,' press in drive ',press(lp1),lp1
c     print *,' press(9) in drive ',press(1),press(9)
      grav=980.0
      p(1)=0.0
      p(llp1)=pstar
      p(llp2)=1.0e6
      p(llp3)=2.0e6
      do 11 i=1,l
      q(i)=r(i)*1.66
      qqo3(i)=qo3(i)*1.66
11    continue
      do 12 i=1,l
      p(2*i)=press(i)
12    continue
      do 13 i=1,lm1
      p(2*i+1)=0.5*(press(i)+press(i+1))
c     print *,' press check ',press(i),p(2*i+1),i
13    continue
c     print *,' press lp1 ',press(2*l + 1)
      do 14 i=1,l
      delp(i)=1./(p(2*i+1)-p(2*i-1))
14    continue
      do 15 i=1,l
      var1(i)=q(i)/(grav*delp(i))
c     print *,' var1 check ',q(i),var1(i)
      var3(i)=qqo3(i)/(grav*delp(i))
c     print *,' var check ',qqo3(i),grav,delp(i)
15    continue
      do 16 i=1,l
      var2(i)=var1(i)*(p(2*i+1)+p(2*i-1))/(2.0*1.01325e6)
      var4(i)=var3(i)*(p(2*i+1)+p(2*i-1))/(2.0*1.01325e6)
16    continue
      do 2 i=2,l
      t(i)=0.5*(temp(i)+temp(i-1))
2     continue
      t(1)=temp(1)
      t(lp1)=temp(lp1)
      grav1=1./(1.01325e6*grav)
      tmp4=temp(1)*temp(1)
      tmp4=tmp4*tmp4
      call co2trn
      call fstrd
      z3=0.
      do 31 i=1,l
      z3=z3+heatra(i)/(8.427*delp(i))
31    continue
      return
      end
