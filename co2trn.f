c $Log$
c Revision 1.1  2003/08/13 01:24:20  dix043
c Initial revision
c
c Revision 1.1  1996/10/17  05:14:11  mrd
c Initial revision
c
c Revision 1.1  91/07/19  11:40:39  ldr
c Initial revision
c 
      subroutine co2trn
c     ******************************************************************
c     subroutine co2trn - with new fels and schwarzkopf(1982)
c                         transmission tables.
c     interpolation of the co2 transmission function to the input
c     temperature profile and surface pressure ( the result is the
c     array co21).
c
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
      common/co2dat/co2r1(lp1,lp1),co2r8(lp1,lp1),cdtr1(lp1,lp1),
     * cdtr8(lp1,lp1),stemp(lp1),po4(lp1),ppr2(l)
      common/kdata/pm(lp1,lp1),p1(lp1,lp1),psk(lp1,lp1),co21(lp1,lp1),
     * q(lp1),p(llp3),delp(l),t(lp1),var1(l),var2(l),var3(l),var4(l),
     * grav,grav1,tmp4
      common /co2new/t1(lp1,lp1,3),t2(lp1,lp1,3)
      dimension dift(lp1,lp1),tdav(lp1),tstdav(lp1)
      dimension t0(lp1,lp1,3),po6(lp1)
      save
c
c a1 a2 are the pressure weights to translate the standard
c profile pressures to the new profile pressures.
c
      pstar = press(lp1)
      p30 = 30000.00e0
      a1=(pstar-810600.098e0)/202649.902e0
      a2=(1013250.e0-pstar)/202649.902e0
c
c calculate the standard transmissions for new p surface.
c
      do 10 k = 1,3
      do 10 j=1,lp1
      do 10 i=1,lp1
      t0(i,j,k) = a1*t1(i,j,k) + a2*t2(i,j,k)
10    continue
c
c po4 is the temperature interplolation weight for the transmissions
c
      do 205 i=1,l
      pal = log(press(i)/p30)
      po8 = 1. + exp(0.8*pal)
      po6(i) = po8 * exp(0.2*pal)
      po4(i) = exp(0.4*log(press(i)))
205   continue
c
c tstdav and tdav are dummy arrays for the interpolation integration
c of the pressure weights
c
      tstdav(1)=0.0
      tdav(1)=0.0
      do 21 i=1,l
      tdav(i+1)=po6(i)*(temp(i)-stemp(i))/delp(i)+tdav(i)
      tstdav(i+1)=po6(i)/delp(i)+tstdav(i)
21    continue
c
c calculate differences
c
      do 22 i=1,lp1
      do 22 j=1,lp1
      dift(i,j)=tstdav(i)-tstdav(j)
22    continue
c
c p1 is a diagonal matrix of 1.0e-16 so as to not divide by zero.
c
      do 23 i=1,lp1
      do 23 j=1,lp1
      dift(i,j)=(tdav(i)-tdav(j))/(dift(i,j)+p1(i,j))
23    continue
c
c calculate the final transmissions
c
c  bandwith conversion assumes square band
c
      fac35 = 350.
      fac20 = 200.
      fac15 = fac35 - fac20
c
      do 24 j=1,lp1
      dt = temp(j) - 250.
      f = 1. + .1833e-03*dt*(1. - 0.01364*dt)
      do 24 i=1,lp1
      co = t0(i,j,1) + dift(i,j)*(t0(i,j,2) + dift(i,j)*t0(i,j,3))
      co21(i,j) = f * (co*fac35 - fac15)/fac20
24    continue
c---------------------------------------------------------------------*
c     correction refer tlh                                            *
c---------------------------------------------------------------------*
      co21(1,1)=1.0
      return
c end of co2trn- new fels and schwarzkopf(1982) transmissions
      end
