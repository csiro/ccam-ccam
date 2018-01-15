c     subroutine spa88 computes exact cts heating rates and fluxes and
c  corresponding cts emissivity quantities for h2o,co2 and o3.
c          inputs:                (common blocks) 
c       acomb,bcomb,apcm,bpcm                  bdcomb 
c       atpcm,btpcm,betacm                     bdcomb 
c       betinw                                 bdwide 
c       temp,press                             radisw 
c       var1,var2,p,delp,delp2                 kdacom 
c       totvo2,to3spc,co2sp1,co2sp2,co2sp      tfcom
c       cldfac                                 cldcom 
c       sko2d                                  tabcom 
c       sorc,csour,osour                       srccom 
c           outputs:  
c       excts,exctsn,ctso3                     tfcom
c        gxcts,fctsg                           rdflux 
c           called by:  
c       fst88 
c            calls: 
c 
      subroutine spa88
c 
CDIR$ TASK COMMON CLDCOM
CDIR$ TASK COMMON CLRFLX
CDIR$ TASK COMMON KDACOM
CDIR$ TASK COMMON LWOUT
CDIR$ TASK COMMON RADISW
CDIR$ TASK COMMON RDFLUX
CDIR$ TASK COMMON SRCCOM
CDIR$ TASK COMMON TFCOM
c CDIR$ TASK COMMON VTEMP

      use cldcom_m
      use kdacom_m
      use lwout_m
      use newmpar_m
      use radisw_m
      use rdflux_m
      use srccom_m
      use tfcom_m

      include 'hcon.h'
      include 'rdparm.h'
      include 'rnddta.h'
c 
      real      phitmp(imax,l),psitmp(imax,l),
     .          fac2(imax,l),
     &          ctmp(imax,lp1),x(imax,l),y(imax,l),
     &          topm(imax,l),topphi(imax,l),
     &          ctmp3(imax,lp1),ctmp2(imax,lp1)
      double precision tt(imax,l),fac1(imax,l) ! double precision
      dimension f(imax,l),ff(imax,l),ag(imax,l),agg(imax,l)
c---compute temperature quantities for use in program
      do 101 k=1,l
      do 101 i=1,imax 
      x(i,k)=temp(i,k)-h25e2
      y(i,k)=x(i,k)*x(i,k)
101   continue
c---initialize ctmp(i,1),ctmp2(i,1),ctmp3(i,1) to unity; these are 
c   transmission fctns at the top.
      do 345 i=1,imax 
      ctmp(i,1)=one 
      ctmp2(i,1)=1.
      ctmp3(i,1)=1.
345   continue
c***begin loop on frequency bands (1)***
c 
c---calculation for band 1 (combined band 1)  
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
      do 301 i=1,imax*l
      f(i,1)=h44194m2*(apcm(1)*x(i,1)+bpcm(1)*y(i,1)) 
      ff(i,1)=h44194m2*(atpcm(1)*x(i,1)+btpcm(1)*y(i,1))
      ag(i,1)=(h1p41819+f(i,1))*f(i,1)+one
      agg(i,1)=(h1p41819+ff(i,1))*ff(i,1)+one
      phitmp(i,1)=var1(i,1)*(((( ag(i,1)*ag(i,1))**2)**2)**2)
      psitmp(i,1)=var2(i,1)*(((( agg(i,1)*agg(i,1))**2)**2)**2)
301   continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 315 i=1,imax 
      topm(i,1)=phitmp(i,1) 
      topphi(i,1)=psitmp(i,1) 
315   continue
      do 319 k=2,l
      do 317 i=1,imax 
      topm(i,k)=topm(i,k-1)+phitmp(i,k) 
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k) 
317   continue
319   continue
c---tt is the cloud-free cts transmission function
      do 321 i=1,imax*l
      fac1(i,1)=acomb(1)*topm(i,1)
      fac2(i,1)=fac1(i,1)*topm(i,1)/(bcomb(1)*topphi(i,1))
      tt(i,1)=exp(hm1ez*fac1(i,1)/sqrt(1.+fac2(i,1)))
      ctmp(i,2)=tt(i,1)*cldfac(i,2,1) 
      ctmp2(i,2)=tt(i,1)
321   continue
c---excts is the cts cooling rate (exctsn is the rate for band nl)
      do 353 i=1,imax*l
      exctsn(i,1,1)=radcon*delp(i,1)*sorc(i,1,1)*(ctmp(i,2)-ctmp(i,1)) 
      excts(i,1)=exctsn(i,1,1) 
      exctsclr(i,1)=radcon*delp(i,1)*sorc(i,1,1)*(ctmp2(i,2)-ctmp2(i,1)) 
353   continue
c---gxcts is the exact cts top flux (fctsg is the flux for band nl)
      do 361 i=1,imax 
      fctsg(i,1)=tt(i,l)*sorc(i,l,1)+
     &   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) + 
     &   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     &   (sorc(i,lp1,1)-sorc(i,l,1))
      gxctsclr(i)=fctsg(i,1)
      fctsg(i,1)=cldfac(i,lp1,1)*fctsg(i,1)
      gxcts(i)=fctsg(i,1)
361   continue
c 
c 
c-----calculation for band 2 (combined band 2)
c 
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
      do 401 i=1,imax*l
      f(i,1)=h44194m2*(apcm(2)*x(i,1)+bpcm(2)*y(i,1)) 
      ff(i,1)=h44194m2*(atpcm(2)*x(i,1)+btpcm(2)*y(i,1))
      ag(i,1)=(h1p41819+f(i,1))*f(i,1)+one
      agg(i,1)=(h1p41819+ff(i,1))*ff(i,1)+one 
      phitmp(i,1)=var1(i,1)*(((( ag(i,1)*ag(i,1))**2)**2)**2)
      psitmp(i,1)=var2(i,1)*(((( agg(i,1)*agg(i,1))**2)**2)**2)
401   continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 415 i=1,imax 
      topm(i,1)=phitmp(i,1) 
      topphi(i,1)=psitmp(i,1) 
415   continue
      do 419 k=2,l
      do 417 i=1,imax 
      topm(i,k)=topm(i,k-1)+phitmp(i,k) 
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k) 
417   continue
419   continue
c---tt is the cloud-free cts transmission function
      do 421 i=1,imax*l
      fac1(i,1)=acomb(2)*topm(i,1)
      fac2(i,1)=fac1(i,1)*topm(i,1)/(bcomb(2)*topphi(i,1))
      tt(i,1)=exp(hm1ez*fac1(i,1)/sqrt(1.+fac2(i,1)))
      ctmp(i,2)=tt(i,1)*cldfac(i,2,1) 
      ctmp2(i,2)=tt(i,1)
421   continue
c---excts is the cts cooling rate (exctsn is the rate for band nl)
      do 453 i=1,imax*l
      exctsn(i,1,2)=radcon*delp(i,1)*sorc(i,1,2)*(ctmp(i,2)-ctmp(i,1)) 
      excts(i,1)=excts(i,1)+exctsn(i,1,2) 
      exctsclr(i,1)=exctsclr(i,1)+
     & radcon*delp(i,1)*sorc(i,1,2)*(ctmp2(i,2)-ctmp2(i,1)) 
453   continue
c---gxcts is the exact cts top flux (fctsg is the flux for band nl)
      do 461 i=1,imax 
      fctsg(i,2)=tt(i,l)*sorc(i,l,2)+
     &   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) + 
     &   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     &   (sorc(i,lp1,2)-sorc(i,l,2))
      gxctsclr(i)=gxctsclr(i)+fctsg(i,2)
      fctsg(i,2)=cldfac(i,lp1,1)*fctsg(i,2)
      gxcts(i)=gxcts(i)+fctsg(i,2)
461   continue
c 
c-----calculation for band 3 (combined band 3)
c 
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
      do 501 i=1,imax*l
      f(i,1)=h44194m2*(apcm(3)*x(i,1)+bpcm(3)*y(i,1)) 
      ff(i,1)=h44194m2*(atpcm(3)*x(i,1)+btpcm(3)*y(i,1))
      ag(i,1)=(h1p41819+f(i,1))*f(i,1)+one
      agg(i,1)=(h1p41819+ff(i,1))*ff(i,1)+one 
      phitmp(i,1)=var1(i,1)*(((( ag(i,1)*ag(i,1))**2)**2)**2)
      psitmp(i,1)=var2(i,1)*(((( agg(i,1)*agg(i,1))**2)**2)**2)
501   continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 515 i=1,imax 
      topm(i,1)=phitmp(i,1) 
      topphi(i,1)=psitmp(i,1) 
515   continue
      do 519 k=2,l
      do 517 i=1,imax 
      topm(i,k)=topm(i,k-1)+phitmp(i,k) 
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k) 
517   continue
519   continue
c---tt is the cloud-free cts transmission function
      do 521 i=1,imax*l
      fac1(i,1)=acomb(3)*topm(i,1)
      fac2(i,1)=fac1(i,1)*topm(i,1)/(bcomb(3)*topphi(i,1))
      tt(i,1)=exp(hm1ez*fac1(i,1)/sqrt(1.+fac2(i,1)))
      ctmp(i,2)=tt(i,1)*cldfac(i,2,1) 
      ctmp2(i,2)=tt(i,1)
521   continue
c---excts is the cts cooling rate (exctsn is the rate for band nl)
      do 553 i=1,imax*l
      exctsn(i,1,3)=radcon*delp(i,1)*sorc(i,1,3)*(ctmp(i,2)-ctmp(i,1)) 
      excts(i,1)=excts(i,1)+exctsn(i,1,3) 
      exctsclr(i,1)=exctsclr(i,1)+
     & radcon*delp(i,1)*sorc(i,1,3)*(ctmp2(i,2)-ctmp2(i,1)) 
553   continue
c---gxcts is the exact cts top flux (fctsg is the flux for band nl)
      do 561 i=1,imax 
      fctsg(i,3)=tt(i,l)*sorc(i,l,3)+
     &   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) + 
     &   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     &   (sorc(i,lp1,3)-sorc(i,l,3))
      gxctsclr(i)=gxctsclr(i)+fctsg(i,3)
      fctsg(i,3)=cldfac(i,lp1,1)*fctsg(i,3)
      gxcts(i)=gxcts(i)+fctsg(i,3)
561   continue
c 
c-----calculation for band 4 (combined band 4)
c 
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
      do 601 i=1,imax*l
      f(i,1)=h44194m2*(apcm(4)*x(i,1)+bpcm(4)*y(i,1)) 
      ff(i,1)=h44194m2*(atpcm(4)*x(i,1)+btpcm(4)*y(i,1))
      ag(i,1)=(h1p41819+f(i,1))*f(i,1)+one
      agg(i,1)=(h1p41819+ff(i,1))*ff(i,1)+one 
      phitmp(i,1)=var1(i,1)*(((( ag(i,1)*ag(i,1))**2)**2)**2)
      psitmp(i,1)=var2(i,1)*(((( agg(i,1)*agg(i,1))**2)**2)**2)
601   continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 615 i=1,imax 
      topm(i,1)=phitmp(i,1) 
      topphi(i,1)=psitmp(i,1) 
615   continue
      do 619 k=2,l
      do 617 i=1,imax 
      topm(i,k)=topm(i,k-1)+phitmp(i,k) 
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k) 
617   continue
619   continue
c---tt is the cloud-free cts transmission function
      do 621 i=1,imax*l
      fac1(i,1)=acomb(4)*topm(i,1)
      fac2(i,1)=fac1(i,1)*topm(i,1)/(bcomb(4)*topphi(i,1))
      tt(i,1)=exp(hm1ez*fac1(i,1)/sqrt(1.+fac2(i,1)))
      ctmp(i,2)=tt(i,1)*cldfac(i,2,1) 
      ctmp2(i,2)=tt(i,1)
621   continue
c---excts is the cts cooling rate (exctsn is the rate for band nl)
      do 653 i=1,imax*l
      exctsn(i,1,4)=radcon*delp(i,1)*sorc(i,1,4)*(ctmp(i,2)-ctmp(i,1)) 
      excts(i,1)=excts(i,1)+exctsn(i,1,4) 
      exctsclr(i,1)=exctsclr(i,1)+
     & radcon*delp(i,1)*sorc(i,1,4)*(ctmp2(i,2)-ctmp2(i,1)) 
653   continue
c---gxcts is the exact cts top flux (fctsg is the flux for band nl)
      do 661 i=1,imax 
      fctsg(i,4)=tt(i,l)*sorc(i,l,4)+
     &   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) + 
     &   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     &   (sorc(i,lp1,4)-sorc(i,l,4))
      gxctsclr(i)=gxctsclr(i)+fctsg(i,4)
      fctsg(i,4)=cldfac(i,lp1,1)*fctsg(i,4)
      gxcts(i)=gxcts(i)+fctsg(i,4)
661   continue
c 
c-----calculation for band 5 (combined band 5)
c 
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
      do 701 i=1,imax*l
      f(i,1)=h44194m2*(apcm(5)*x(i,1)+bpcm(5)*y(i,1)) 
      ff(i,1)=h44194m2*(atpcm(5)*x(i,1)+btpcm(5)*y(i,1))
      ag(i,1)=(h1p41819+f(i,1))*f(i,1)+one
      agg(i,1)=(h1p41819+ff(i,1))*ff(i,1)+one 
      phitmp(i,1)=var1(i,1)*(((( ag(i,1)*ag(i,1))**2)**2)**2)
      psitmp(i,1)=var2(i,1)*(((( agg(i,1)*agg(i,1))**2)**2)**2)
701   continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 715 i=1,imax 
      topm(i,1)=phitmp(i,1) 
      topphi(i,1)=psitmp(i,1) 
715   continue
      do 719 k=2,l
      do 717 i=1,imax 
      topm(i,k)=topm(i,k-1)+phitmp(i,k) 
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k) 
717   continue
719   continue
c---tt is the cloud-free cts transmission function
      do 721 i=1,imax*l
      fac1(i,1)=acomb(5)*topm(i,1)
      fac2(i,1)=fac1(i,1)*topm(i,1)/(bcomb(5)*topphi(i,1))
      tt(i,1)=exp(hm1ez*(fac1(i,1)/sqrt(one+fac2(i,1))+
     &           betacm(5)*totvo2(i,2)*sko2d))
      ctmp(i,2)=tt(i,1)*cldfac(i,2,1) 
      ctmp2(i,2)=tt(i,1)
721   continue
c---excts is the cts cooling rate (exctsn is the rate for band nl)
      do 753 i=1,imax*l
      exctsn(i,1,5)=radcon*delp(i,1)*sorc(i,1,5)*(ctmp(i,2)-ctmp(i,1)) 
      excts(i,1)=excts(i,1)+exctsn(i,1,5) 
      exctsclr(i,1)=exctsclr(i,1)+
     & radcon*delp(i,1)*sorc(i,1,5)*(ctmp2(i,2)-ctmp2(i,1)) 
753   continue
c---gxcts is the exact cts top flux (fctsg is the flux for band nl)
      do 761 i=1,imax 
      fctsg(i,5)=tt(i,l)*sorc(i,l,5)+
     &   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) + 
     &   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     &   (sorc(i,lp1,5)-sorc(i,l,5))
      gxctsclr(i)=gxctsclr(i)+fctsg(i,5)
      fctsg(i,5)=cldfac(i,lp1,1)*fctsg(i,5)
      gxcts(i)=gxcts(i)+fctsg(i,5)
761   continue
c 
c-----calculation for band 6 (combined band 6)
c 
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
      do 801 i=1,imax*l
      f(i,1)=h44194m2*(apcm(6)*x(i,1)+bpcm(6)*y(i,1)) 
      ff(i,1)=h44194m2*(atpcm(6)*x(i,1)+btpcm(6)*y(i,1))
      ag(i,1)=(h1p41819+f(i,1))*f(i,1)+one
      agg(i,1)=(h1p41819+ff(i,1))*ff(i,1)+one 
      phitmp(i,1)=var1(i,1)*(((( ag(i,1)*ag(i,1))**2)**2)**2)
      psitmp(i,1)=var2(i,1)*(((( agg(i,1)*agg(i,1))**2)**2)**2)
801   continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 815 i=1,imax 
      topm(i,1)=phitmp(i,1) 
      topphi(i,1)=psitmp(i,1) 
815   continue
      do 819 k=2,l
      do 817 i=1,imax 
      topm(i,k)=topm(i,k-1)+phitmp(i,k) 
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k) 
817   continue
819   continue
c---tt is the cloud-free cts transmission function
      do 821 i=1,imax*l
      fac1(i,1)=acomb(6)*topm(i,1)
      fac2(i,1)=fac1(i,1)*topm(i,1)/(bcomb(6)*topphi(i,1))
      tt(i,1)=exp(hm1ez*(fac1(i,1)/sqrt(one+fac2(i,1))+
     &           betacm(6)*totvo2(i,2)*sko2d))
      ctmp(i,2)=tt(i,1)*cldfac(i,2,1) 
      ctmp2(i,2)=tt(i,1)
821   continue
c---excts is the cts cooling rate (exctsn is the rate for band nl)
      do 853 i=1,imax*l
      exctsn(i,1,6)=radcon*delp(i,1)*sorc(i,1,6)*(ctmp(i,2)-ctmp(i,1)) 
      excts(i,1)=excts(i,1)+exctsn(i,1,6) 
      exctsclr(i,1)=exctsclr(i,1)+
     & radcon*delp(i,1)*sorc(i,1,6)*(ctmp2(i,2)-ctmp2(i,1)) 
853   continue
c---gxcts is the exact cts top flux (fctsg is the flux for band nl)
      do 861 i=1,imax 
      fctsg(i,6)=tt(i,l)*sorc(i,l,6)+
     &   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) + 
     &   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     &   (sorc(i,lp1,6)-sorc(i,l,6))
      gxctsclr(i)=gxctsclr(i)+fctsg(i,6)
      fctsg(i,6)=cldfac(i,lp1,1)*fctsg(i,6)
      gxcts(i)=gxcts(i)+fctsg(i,6)
861   continue
c 
c-----calculation for band 7 (combined band 7)
c 
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
      do 901 i=1,imax*l
      f(i,1)=h44194m2*(apcm(7)*x(i,1)+bpcm(7)*y(i,1)) 
      ff(i,1)=h44194m2*(atpcm(7)*x(i,1)+btpcm(7)*y(i,1))
      ag(i,1)=(h1p41819+f(i,1))*f(i,1)+one
      agg(i,1)=(h1p41819+ff(i,1))*ff(i,1)+one 
      phitmp(i,1)=var1(i,1)*(((( ag(i,1)*ag(i,1))**2)**2)**2)
      psitmp(i,1)=var2(i,1)*(((( agg(i,1)*agg(i,1))**2)**2)**2)
901   continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 915 i=1,imax 
      topm(i,1)=phitmp(i,1) 
      topphi(i,1)=psitmp(i,1) 
915   continue
      do 919 k=2,l
      do 917 i=1,imax 
      topm(i,k)=topm(i,k-1)+phitmp(i,k) 
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k) 
917   continue
919   continue
c---tt is the cloud-free cts transmission function
      do 921 i=1,imax*l
      fac1(i,1)=acomb(7)*topm(i,1)
      fac2(i,1)=fac1(i,1)*topm(i,1)/(bcomb(7)*topphi(i,1))
      tt(i,1)=exp(hm1ez*(fac1(i,1)/sqrt(one+fac2(i,1))+
     &           betacm(7)*totvo2(i,2)*sko2d))
      ctmp(i,2)=tt(i,1)*cldfac(i,2,1) 
      ctmp2(i,2)=tt(i,1)
921   continue
c---excts is the cts cooling rate (exctsn is the rate for band nl)
      do 953 i=1,imax*l
      exctsn(i,1,7)=radcon*delp(i,1)*sorc(i,1,7)*(ctmp(i,2)-ctmp(i,1)) 
      excts(i,1)=excts(i,1)+exctsn(i,1,7) 
      exctsclr(i,1)=exctsclr(i,1)+
     & radcon*delp(i,1)*sorc(i,1,7)*(ctmp2(i,2)-ctmp2(i,1)) 
953   continue
c---gxcts is the exact cts top flux (fctsg is the flux for band nl)
      do 961 i=1,imax 
      fctsg(i,7)=tt(i,l)*sorc(i,l,7)+
     &   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) + 
     &   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     &   (sorc(i,lp1,7)-sorc(i,l,7))
      gxctsclr(i)=gxctsclr(i)+fctsg(i,7)
      fctsg(i,7)=cldfac(i,lp1,1)*fctsg(i,7)
      gxcts(i)=gxcts(i)+fctsg(i,7)
961   continue
c 
c-----calculation for band 8 (combined band 8)
c 
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
      do 1001 i=1,imax*l
      f(i,1)=h44194m2*(apcm(8)*x(i,1)+bpcm(8)*y(i,1)) 
      ff(i,1)=h44194m2*(atpcm(8)*x(i,1)+btpcm(8)*y(i,1))
      ag(i,1)=(h1p41819+f(i,1))*f(i,1)+one
      agg(i,1)=(h1p41819+ff(i,1))*ff(i,1)+one 
      phitmp(i,1)=var1(i,1)*(((( ag(i,1)*ag(i,1))**2)**2)**2)
      psitmp(i,1)=var2(i,1)*(((( agg(i,1)*agg(i,1))**2)**2)**2)
1001  continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 1015 i=1,imax 
      topm(i,1)=phitmp(i,1) 
      topphi(i,1)=psitmp(i,1) 
1015  continue
      do 1019 k=2,l
      do 1017 i=1,imax 
      topm(i,k)=topm(i,k-1)+phitmp(i,k) 
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k) 
1017  continue
1019  continue
c---tt is the cloud-free cts transmission function
      do 1021 i=1,imax*l
      fac1(i,1)=acomb(8)*topm(i,1)
      fac2(i,1)=fac1(i,1)*topm(i,1)/(bcomb(8)*topphi(i,1))
      tt(i,1)=exp(hm1ez*(fac1(i,1)/sqrt(one+fac2(i,1))+
     &           betacm(8)*totvo2(i,2)*sko2d))
      ctmp(i,2)=tt(i,1)*cldfac(i,2,1) 
      ctmp2(i,2)=tt(i,1)
1021  continue
c---excts is the cts cooling rate (exctsn is the rate for band nl)
      do 1053 i=1,imax*l
      exctsn(i,1,8)=radcon*delp(i,1)*sorc(i,1,8)*(ctmp(i,2)-ctmp(i,1)) 
      excts(i,1)=excts(i,1)+exctsn(i,1,8) 
      exctsclr(i,1)=exctsclr(i,1)+
     & radcon*delp(i,1)*sorc(i,1,8)*(ctmp2(i,2)-ctmp2(i,1)) 
1053  continue
c---gxcts is the exact cts top flux (fctsg is the flux for band nl)
      do 1061 i=1,imax 
      fctsg(i,8)=tt(i,l)*sorc(i,l,8)+
     &   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) + 
     &   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     &   (sorc(i,lp1,8)-sorc(i,l,8))
      gxctsclr(i)=gxctsclr(i)+fctsg(i,8)
      fctsg(i,8)=cldfac(i,lp1,1)*fctsg(i,8)
      gxcts(i)=gxcts(i)+fctsg(i,8)
1061  continue
c 
c-----calculation for band 9 ( 560-670 cm-1; includes co2)
c 
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
      do 1101 i=1,imax*l
      f(i,1)=h44194m2*(apcm(9)*x(i,1)+bpcm(9)*y(i,1)) 
      ff(i,1)=h44194m2*(atpcm(9)*x(i,1)+btpcm(9)*y(i,1))
      ag(i,1)=(h1p41819+f(i,1))*f(i,1)+one
      agg(i,1)=(h1p41819+ff(i,1))*ff(i,1)+one 
      phitmp(i,1)=var1(i,1)*(((( ag(i,1)*ag(i,1))**2)**2)**2)
      psitmp(i,1)=var2(i,1)*(((( agg(i,1)*agg(i,1))**2)**2)**2)
1101  continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 1115 i=1,imax 
      topm(i,1)=phitmp(i,1) 
      topphi(i,1)=psitmp(i,1) 
1115  continue
      do 1119 k=2,l
      do 1117 i=1,imax 
      topm(i,k)=topm(i,k-1)+phitmp(i,k) 
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k) 
1117  continue
1119  continue
c---tt is the cloud-free cts transmission function
      do 1121 i=1,imax*l
      fac1(i,1)=acomb(9)*topm(i,1)
      fac2(i,1)=fac1(i,1)*topm(i,1)/(bcomb(9)*topphi(i,1))
      tt(i,1)=exp(hm1ez*(fac1(i,1)/sqrt(one+fac2(i,1))+
     &           betacm(9)*totvo2(i,2)*sko2d))*co2sp1(i,2)
      ctmp(i,2)=tt(i,1)*cldfac(i,2,1) 
      ctmp2(i,2)=tt(i,1)
1121  continue
c---excts is the cts cooling rate (exctsn is the rate for band nl)
      do 1153 i=1,imax*l
      exctsn(i,1,9)=radcon*delp(i,1)*sorc(i,1,9)*(ctmp(i,2)-ctmp(i,1)) 
      excts(i,1)=excts(i,1)+exctsn(i,1,9) 
      exctsclr(i,1)=exctsclr(i,1)+
     & radcon*delp(i,1)*sorc(i,1,9)*(ctmp2(i,2)-ctmp2(i,1)) 
1153  continue
c---gxcts is the exact cts top flux (fctsg is the flux for band nl)
      do 1161 i=1,imax 
      fctsg(i,9)=tt(i,l)*sorc(i,l,9)+
     &   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) + 
     &   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     &   (sorc(i,lp1,9)-sorc(i,l,9))
      gxctsclr(i)=gxctsclr(i)+fctsg(i,9)
      fctsg(i,9)=cldfac(i,lp1,1)*fctsg(i,9)
      gxcts(i)=gxcts(i)+fctsg(i,9)
1161  continue
c 
c-----calculation for band 10 (670-800 cm-1; includes co2)
c 
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
      do 1201 i=1,imax*l
      f(i,1)=h44194m2*(apcm(10)*x(i,1)+bpcm(10)*y(i,1)) 
      ff(i,1)=h44194m2*(atpcm(10)*x(i,1)+btpcm(10)*y(i,1))
      ag(i,1)=(h1p41819+f(i,1))*f(i,1)+one
      agg(i,1)=(h1p41819+ff(i,1))*ff(i,1)+one 
      phitmp(i,1)=var1(i,1)*(((( ag(i,1)*ag(i,1))**2)**2)**2)
      psitmp(i,1)=var2(i,1)*(((( agg(i,1)*agg(i,1))**2)**2)**2)
1201  continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 1215 i=1,imax 
      topm(i,1)=phitmp(i,1) 
      topphi(i,1)=psitmp(i,1) 
1215  continue
      do 1219 k=2,l
      do 1217 i=1,imax 
      topm(i,k)=topm(i,k-1)+phitmp(i,k) 
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k) 
1217  continue
1219  continue
c---tt is the cloud-free cts transmission function
      do 1221 i=1,imax*l
      fac1(i,1)=acomb(10)*topm(i,1)
      fac2(i,1)=fac1(i,1)*topm(i,1)/(bcomb(10)*topphi(i,1))
      tt(i,1)=exp(hm1ez*(fac1(i,1)/sqrt(one+fac2(i,1))+
     &           betacm(10)*totvo2(i,2)*sko2d))*co2sp2(i,2)
      ctmp(i,2)=tt(i,1)*cldfac(i,2,1) 
      ctmp2(i,2)=tt(i,1)
1221  continue
c---excts is the cts cooling rate (exctsn is the rate for band nl)
      do 1253 i=1,imax*l
      exctsn(i,1,10)=radcon*delp(i,1)*sorc(i,1,10)*
     &               (ctmp(i,2)-ctmp(i,1))
      excts(i,1)=excts(i,1)+exctsn(i,1,10) 
      exctsclr(i,1)=exctsclr(i,1)+
     & radcon*delp(i,1)*sorc(i,1,10)*(ctmp2(i,2)-ctmp2(i,1)) 
1253  continue
c---gxcts is the exact cts top flux (fctsg is the flux for band nl)
      do 1261 i=1,imax 
      fctsg(i,10)=tt(i,l)*sorc(i,l,10)+
     &   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) + 
     &   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     &   (sorc(i,lp1,10)-sorc(i,l,10))
      gxctsclr(i)=gxctsclr(i)+fctsg(i,10)
      fctsg(i,10)=cldfac(i,lp1,1)*fctsg(i,10)
      gxcts(i)=gxcts(i)+fctsg(i,10)
1261  continue
c 
c-----calculation for band 11 (800-900 cm-1)
c 
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
      do 1301 i=1,imax*l
      f(i,1)=h44194m2*(apcm(11)*x(i,1)+bpcm(11)*y(i,1)) 
      ff(i,1)=h44194m2*(atpcm(11)*x(i,1)+btpcm(11)*y(i,1))
      ag(i,1)=(h1p41819+f(i,1))*f(i,1)+one
      agg(i,1)=(h1p41819+ff(i,1))*ff(i,1)+one 
      phitmp(i,1)=var1(i,1)*(((( ag(i,1)*ag(i,1))**2)**2)**2)
      psitmp(i,1)=var2(i,1)*(((( agg(i,1)*agg(i,1))**2)**2)**2)
1301  continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 1315 i=1,imax 
      topm(i,1)=phitmp(i,1) 
      topphi(i,1)=psitmp(i,1) 
1315  continue
      do 1319 k=2,l
      do 1317 i=1,imax 
      topm(i,k)=topm(i,k-1)+phitmp(i,k) 
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k) 
1317  continue
1319  continue
c---tt is the cloud-free cts transmission function
      do 1321 i=1,imax*l
      fac1(i,1)=acomb(11)*topm(i,1)
      fac2(i,1)=fac1(i,1)*topm(i,1)/(bcomb(11)*topphi(i,1))
      tt(i,1)=exp(hm1ez*(fac1(i,1)/sqrt(one+fac2(i,1))+
     &           betacm(11)*totvo2(i,2)*sko2d))
      ctmp(i,2)=tt(i,1)*cldfac(i,2,1) 
      ctmp2(i,2)=tt(i,1)
1321  continue
c---excts is the cts cooling rate (exctsn is the rate for band nl)
      do 1353 i=1,imax*l
      exctsn(i,1,11)=radcon*delp(i,1)*sorc(i,1,11)*
     &               (ctmp(i,2)-ctmp(i,1))
      excts(i,1)=excts(i,1)+exctsn(i,1,11) 
      exctsclr(i,1)=exctsclr(i,1)+
     & radcon*delp(i,1)*sorc(i,1,11)*(ctmp2(i,2)-ctmp2(i,1)) 
1353  continue
c---gxcts is the exact cts top flux (fctsg is the flux for band nl)
      do 1361 i=1,imax 
      fctsg(i,11)=tt(i,l)*sorc(i,l,11)+
     &   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) + 
     &   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     &   (sorc(i,lp1,11)-sorc(i,l,11))
      gxctsclr(i)=gxctsclr(i)+fctsg(i,11)
      fctsg(i,11)=cldfac(i,lp1,1)*fctsg(i,11)
      gxcts(i)=gxcts(i)+fctsg(i,11)
1361  continue
c 
c-----calculation for band 12 (900-990 cm-1)
c 
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
      do 1401 i=1,imax*l
      f(i,1)=h44194m2*(apcm(12)*x(i,1)+bpcm(12)*y(i,1)) 
      ff(i,1)=h44194m2*(atpcm(12)*x(i,1)+btpcm(12)*y(i,1))
      ag(i,1)=(h1p41819+f(i,1))*f(i,1)+one
      agg(i,1)=(h1p41819+ff(i,1))*ff(i,1)+one 
      phitmp(i,1)=var1(i,1)*(((( ag(i,1)*ag(i,1))**2)**2)**2)
      psitmp(i,1)=var2(i,1)*(((( agg(i,1)*agg(i,1))**2)**2)**2)
1401  continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 1415 i=1,imax 
      topm(i,1)=phitmp(i,1) 
      topphi(i,1)=psitmp(i,1) 
1415  continue
      do 1419 k=2,l
      do 1417 i=1,imax 
      topm(i,k)=topm(i,k-1)+phitmp(i,k) 
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k) 
1417  continue
1419  continue
c---tt is the cloud-free cts transmission function
      do 1421 i=1,imax*l
      fac1(i,1)=acomb(12)*topm(i,1)
      fac2(i,1)=fac1(i,1)*topm(i,1)/(bcomb(12)*topphi(i,1))
      tt(i,1)=exp(hm1ez*(fac1(i,1)/sqrt(one+fac2(i,1))+
     &           betacm(12)*totvo2(i,2)*sko2d))
      ctmp(i,2)=tt(i,1)*cldfac(i,2,1) 
      ctmp2(i,2)=tt(i,1)
1421  continue
c---excts is the cts cooling rate (exctsn is the rate for band nl)
      do 1453 i=1,imax*l
      exctsn(i,1,12)=radcon*delp(i,1)*sorc(i,1,12)*
     &               (ctmp(i,2)-ctmp(i,1))
      excts(i,1)=excts(i,1)+exctsn(i,1,12) 
      exctsclr(i,1)=exctsclr(i,1)+
     & radcon*delp(i,1)*sorc(i,1,12)*(ctmp2(i,2)-ctmp2(i,1)) 
1453  continue
c---gxcts is the exact cts top flux (fctsg is the flux for band nl)
      do 1461 i=1,imax 
      fctsg(i,12)=tt(i,l)*sorc(i,l,12)+
     &   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) + 
     &   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     &   (sorc(i,lp1,12)-sorc(i,l,12))
      gxctsclr(i)=gxctsclr(i)+fctsg(i,12)
      fctsg(i,12)=cldfac(i,lp1,1)*fctsg(i,12)
      gxcts(i)=gxcts(i)+fctsg(i,12)
1461  continue
c 
c-----calculation for band 13 (990-1070 cm-1; includes o3))
c 
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
      do 1501 i=1,imax*l
      f(i,1)=h44194m2*(apcm(13)*x(i,1)+bpcm(13)*y(i,1)) 
      ff(i,1)=h44194m2*(atpcm(13)*x(i,1)+btpcm(13)*y(i,1))
      ag(i,1)=(h1p41819+f(i,1))*f(i,1)+one
      agg(i,1)=(h1p41819+ff(i,1))*ff(i,1)+one 
      phitmp(i,1)=var1(i,1)*(((( ag(i,1)*ag(i,1))**2)**2)**2)
      psitmp(i,1)=var2(i,1)*(((( agg(i,1)*agg(i,1))**2)**2)**2)
1501  continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 1515 i=1,imax 
      topm(i,1)=phitmp(i,1) 
      topphi(i,1)=psitmp(i,1) 
1515  continue
      do 1519 k=2,l
      do 1517 i=1,imax 
      topm(i,k)=topm(i,k-1)+phitmp(i,k) 
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k) 
1517  continue
1519  continue
c---tt is the cloud-free cts transmission function
      do 1521 i=1,imax*l
      fac1(i,1)=acomb(13)*topm(i,1)
      fac2(i,1)=fac1(i,1)*topm(i,1)/(bcomb(13)*topphi(i,1))
      tt(i,1)=exp(hm1ez*(fac1(i,1)/sqrt(one+fac2(i,1))+
     &           betacm(13)*totvo2(i,2)*sko2d +to3spc(i,1)))
      ctmp(i,2)=tt(i,1)*cldfac(i,2,1) 
      ctmp2(i,2)=tt(i,1)
1521  continue
c---excts is the cts cooling rate (exctsn is the rate for band nl)
      do 1553 i=1,imax*l
      exctsn(i,1,13)=radcon*delp(i,1)*sorc(i,1,13)*
     &               (ctmp(i,2)-ctmp(i,1))
      excts(i,1)=excts(i,1)+exctsn(i,1,13) 
      exctsclr(i,1)=exctsclr(i,1)+
     & radcon*delp(i,1)*sorc(i,1,13)*(ctmp2(i,2)-ctmp2(i,1)) 
1553  continue
c---gxcts is the exact cts top flux (fctsg is the flux for band nl)
      do 1561 i=1,imax 
      fctsg(i,13)=tt(i,l)*sorc(i,l,13)+
     &   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) + 
     &   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     &   (sorc(i,lp1,13)-sorc(i,l,13))
      gxctsclr(i)=gxctsclr(i)+fctsg(i,13)
      fctsg(i,13)=cldfac(i,lp1,1)*fctsg(i,13)
      gxcts(i)=gxcts(i)+fctsg(i,13)
1561  continue
c 
c-----calculation for band 14 (1070-1200 cm-1)
c 
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
      do 1601 i=1,imax*l
      f(i,1)=h44194m2*(apcm(14)*x(i,1)+bpcm(14)*y(i,1)) 
      ff(i,1)=h44194m2*(atpcm(14)*x(i,1)+btpcm(14)*y(i,1))
      ag(i,1)=(h1p41819+f(i,1))*f(i,1)+one
      agg(i,1)=(h1p41819+ff(i,1))*ff(i,1)+one 
      phitmp(i,1)=var1(i,1)*(((( ag(i,1)*ag(i,1))**2)**2)**2)
      psitmp(i,1)=var2(i,1)*(((( agg(i,1)*agg(i,1))**2)**2)**2)
1601  continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 1615 i=1,imax 
      topm(i,1)=phitmp(i,1) 
      topphi(i,1)=psitmp(i,1) 
1615  continue
      do 1619 k=2,l
      do 1617 i=1,imax 
      topm(i,k)=topm(i,k-1)+phitmp(i,k) 
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k) 
1617  continue
1619  continue
c---tt is the cloud-free cts transmission function
      do 1621 i=1,imax*l
      fac1(i,1)=acomb(14)*topm(i,1)
      fac2(i,1)=fac1(i,1)*topm(i,1)/(bcomb(14)*topphi(i,1))
      tt(i,1)=exp(hm1ez*(fac1(i,1)/sqrt(one+fac2(i,1))+
     &           betacm(14)*totvo2(i,2)*sko2d))
      ctmp(i,2)=tt(i,1)*cldfac(i,2,1) 
      ctmp2(i,2)=tt(i,1)
1621  continue
c---excts is the cts cooling rate (exctsn is the rate for band nl)
      do 1653 i=1,imax*l
      exctsn(i,1,14)=radcon*delp(i,1)*sorc(i,1,14)*
     &               (ctmp(i,2)-ctmp(i,1))
      excts(i,1)=excts(i,1)+exctsn(i,1,14) 
      exctsclr(i,1)=exctsclr(i,1)+
     & radcon*delp(i,1)*sorc(i,1,14)*(ctmp2(i,2)-ctmp2(i,1)) 
1653  continue
c---gxcts is the exact cts top flux (fctsg is the flux for band nl)
      do 1661 i=1,imax 
      fctsg(i,14)=tt(i,l)*sorc(i,l,14)+
     &   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) + 
     &   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     &   (sorc(i,lp1,14)-sorc(i,l,14))
      gxctsclr(i)=gxctsclr(i)+fctsg(i,14)
      fctsg(i,14)=cldfac(i,lp1,1)*fctsg(i,14)
      gxcts(i)=gxcts(i)+fctsg(i,14)
1661  continue
c 
c---this is the end of the exact cts computations; at this point
c   excts has its appropriate value. the exact cts flux at the 
c   ground is now evaluated by using the (now available) cts cooling 
c   rate and the exact cts top flux.
c*** compute approximate cts heating rates for 15um and 9.6 um bands
c     (ctso3) 
      do 1711 i=1,imax*l
      ctmp2(i,2)=co2sp(i,2)*cldfac(i,2,1)
      ctmp3(i,2)=to3(i,2,1)*cldfac(i,2,1)
1711  continue
      do 1701 k=1,l
      do 1701 i=1,imax 
      gxcts(i)=gxcts(i)-excts(i,k)*delp2(i,k)*radcon1
      gxctsclr(i)=gxctsclr(i)-exctsclr(i,k)*delp2(i,k)*radcon1
c     if ( i.eq.1 ) print*, ' ctso3 ',
c    &     radcon*delp(i,k)*csour(i,k)*(ctmp2(i,k+1)-ctmp2(i,k)) ,
c    &     radcon*delp(i,k)*osour(i,k)*(ctmp3(i,k+1)-ctmp3(i,k)) 
      ctso3(i,k)=radcon*delp(i,k)*
     &     (csour(i,k)*(ctmp2(i,k+1)-ctmp2(i,k)) +
     &      osour(i,k)*(ctmp3(i,k+1)-ctmp3(i,k)))
      ctso3clr(i,k)=radcon*delp(i,k)*
     &     (csour(i,k)*(co2sp(i,k+1)-co2sp(i,k)) +
     &      osour(i,k)*(to3(i,k+1,1)-to3(i,k,1)))
1701  continue
      return
      end 
