c     ***************************************************************** 
c          subroutine fst88 is the main computation module of the 
c     long-wave radiation code. in it all "emissivity" calculations,
c     including calls to table lookup subroutines. also,after calling 
c     subroutine "spa88", final combined heating rates and ground 
c     flux are obtained.
c     ***************************************************************** 
c     jlm: renamed dt() to dtemp() & removed unused arrays on 18/7/08
c              inputs:  
c        betinw,betawd,ab15wd              bdwide 
c        betad,bo3rnd,ao3rnd               bandta 
c        qh2o,p,delp2,delp,t,var1,var2,    kdacom 
c        var3,var4,cntval                  kdacom 
c        temp,press                        radisw 
c        ind,indx2,kmaxv,source,dsrce      tabcom 
c        skc1r,skc3r,kmaxvm,nrep1,nrep2    tabcom 
c        nst1,nst2,nrp1,nrp2               tabcom 
c        co2nbl,co2sp,co21                 tfcom
c              outputs: 
c        heatra,grnflx,topflx              lwout
c        flx1e1                            rdflux 
c 
c          called by  :    radmn or main pgm
c          calls      :    clo88,e1e288,e3v88,spa88,nlte 
c
c        passed variables:  
c              in e3v88:  
c        emd     =  e3 function for h2o lines (0-560,1200-2200 cm-1)
c                     computed in e3v88 
c        tpl     =  temperature input for e3 calculation in e3v88 
c        empl    =  h2o amount,input for e3 calculation in e3v88
c              in e1e288: 
c        e1cts1  =  e1 function for the (i+1)th level using the 
c                   temperature of the ith data level,computed over 
c                   the frequency range 0-560,1200-2200 cm-1. (e1cts1-
c                   e1ctw1) is used in obtaining the flux at the top
c                   in the 0-160,1200-2200 cm-1 range (flx1e1). 
c        e1cts2  =  e1 function for the ith level, using the temp. of 
c                   the ith data level,computed over the frequency range
c                   0-560,1200-2200 cm-1. (e1cts2-e1ctw2) is also used
c                   in obtaining the flux at the top in the 0-160,. 
c                   1200-2200 cm-1 range. 
c        e1flx   =  e1 fctn. for the ith level,using the temperature at 
c                   the top of the atmosphere. computed over the freq.
c                   range 0-560,1200-2200 cm-1. used for q(approx) term.
c                   (in common block tfcom) 
c        e1ctw1  =  like e1cts1,but computed over the 160-560 cm-1 range

c                   and used for q(approx,cts) calculation
c        e1ctw2  =  like e1cts2,but computed over the 160-560 cm-1 range
c                   and used for q(approx,cts) calculation
c        fxo     =  temperature index used for e1 function and also 
c                   used for source function calc. in fst88.
c        dtemp   =  temp. diff.between model temps. and temps. at 
c                   tabular values of e1 and source fctns. used in
c                   fst88 and in e1 function calc.
*VOCL TOTAL,REPEAT(999999)
      subroutine fst88
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

      use cc_mpi, only : mydiag,myid
      use cldcom_m
      use diag_m
      use kdacom_m
      use lwout_m
      use newmpar_m
      use parm_m
      use radisw_m
      use rdflux_m
      use srccom_m
      use tabcom_m
      use tfcom_m
      implicit none
      include 'hcon.h'
      include 'rdparm.h'
      include 'rnddta.h'
      integer ixo(imax,lp1)  !itop(imax),ibot(imax),indtc(imax)
      real vtmp1(imax,lp1),vtmp2(imax,lp1),
     & vtmp3(imax,lp1),c(imax,llp1),alp(imax,llp1), !dsorc(imax,lp1),
     & totphi(imax,lp1),toto3(imax,lp1),tphio3(imax,lp1),rlog(imax,l),
     & delptc(imax),ptop(imax),pbot(imax),ftop(imax),
     & fbot(imax)
c     real  over1d(imax,lp1m)
      dimension vtmp3x(imax+1,lp1),dsorcx(imax+1,lp1)
c---dimension of variables equivalenced to those in vtemp---
      dimension vsum1(imax,lp1),heatem(imax,lp1) !tval(imax,lp1),
      dimension emxx(imax,l)
      dimension csub(imax,llp1) ,csub2(imax,llp1),c2(imax,llp1)
      dimension flx(imax,lp1) !sum(imax,lp1) 
      dimension oss(imax,lp1),css(imax,lp1),ss2(imax,lp1),tc(imax,lp1), 
     & dtc(imax,lp1)
      dimension alpsq1(imax,lp1),alpsq2(imax,lp1) 
      dimension delpr1(imax,lp1),delpr2(imax,lp1) 
      dimension flxnet(imax,lp1) ! flxthk(imax,lp1) 
      dimension z1(imax,lp1),ceval(imax,lp1)
      dimension cevalclr(imax,lp1)
      dimension totevv(imax,lp1)
      dimension avmo3(imax,lp1),avpho3(imax,lp1),fac1(imax,lp1)
      dimension avvo2(imax,lp1),avephj(imax,lp1)
      dimension over(imax,lp1,lp1)
      dimension emisst(imax,lp1,lp1)
c---dimension of variables equivalenced to those in other common blks-- 
c     dimension to31d(imax,lp1m),emi21d(imax,lp1m)
c     dimension co21d(imax,lp1m),emis1d(imax,lp1m)
c     dimension avep1d(imax,lp1m),avep1(imax*lp1m)
c---dimension of variables passed to other subroutines--- 
      dimension e1cts1(imax,lp1),e1cts2(imax,l) 
      dimension e1ctw1(imax,lp1),e1ctw2(imax,l) 
      dimension fxo(imax,lp1),dtemp(imax,lp1)
      dimension emd(imax,llp1),tpl(imax,llp1),empl(imax,llp1) 
c---emx1 is a local variable used as input and output to e1e288--
      dimension emx1(imax)
cx      equivalence (vtmp3,vsum1,emxx,tc,heatem) 
cx      equivalence (totevv,delpr1,oss,flxnet)
cx      equivalence (totphi,delpr2,css,flx,z1)
cx      equivalence (toto3,alpsq1,dtc) 
cx      equivalence (tphio3,alpsq2,ss2,ceval) 
cx      equivalence (alp,csub)
      real flxclr(imax,lp1),vsum1clr(imax,lp1),vtmp1clr(imax,lp1)
      real ctsclr(imax,l)
      real heatemclr(imax,lp1),heatraclr(imax,lp1)
      integer i,icnt,j1,j3,k,kclds,kk,kp,kx
      real alpsq1,alpsq2,avephj,avmo3,avpho3,avvo2
      real c2,ceval,cevalclr,css,csub,csub2,delpr1,delpr2
      real dsorcx,dtc,dtemp,e1cts1,e1cts2,e1ctw1,e1ctw2
      real emd,emisst,empl
      real emx1,emxx,fac1,flx,flxnet,fxo,heatem,oss,over
      real ss2,tc,totevv,tpl,vsum1,vtmp3x,z1
 
c          first section is table lookup for source function and
c     derivative (b and db/dt).also,the nlte co2 source function
c     is obtained 
c 
c---decrementing the index by 9 accounts for the tables beginning
c   at t=100k.

!!!      print*, ' In fst88 '
      do 101 i=1,imax*lp1
      vtmp2(i,1)=aint(temp(i,1)*hp1)
      fxo(i,1)=vtmp2(i,1)-9.
      dtemp(i,1)=temp(i,1)-ten*vtmp2(i,1)
      ixo(i,1)=nint(fxo(i,1))
101   continue
c 
c---source function for combined band 1
      do 4114 i=1,imax 
      do 4114 k=1,lp1
        vtmp3x(i,k)=source(ixo(i,k),1) 
        dsorcx(i,k)=dsrce(ixo(i,k),1) 
4114   continue
      do 4112 k=1,lp1
      do 4112 i=1,imax
      sorc(i,k,1)=vtmp3x(i,k)+dtemp(i,k)*dsorcx(i,k)
4112   continue
c---source function for combined band 2
      do 4214 i=1,imax 
      do 4214 k=1,lp1
        vtmp3x(i,k)=source(ixo(i,k),2) 
        dsorcx(i,k)=dsrce(ixo(i,k),2) 
4214   continue
      do 4212 k=1,lp1
      do 4212 i=1,imax
      sorc(i,k,2)=vtmp3x(i,k)+dtemp(i,k)*dsorcx(i,k)
4212   continue
!!!       print*, ' 4212 '
c---source function for combined band 3
      do 4314 i=1,imax 
      do 4314 k=1,lp1
        vtmp3x(i,k)=source(ixo(i,k),3) 
        dsorcx(i,k)=dsrce(ixo(i,k),3) 
4314   continue
      do 4312 k=1,lp1
      do 4312 i=1,imax
      sorc(i,k,3)=vtmp3x(i,k)+dtemp(i,k)*dsorcx(i,k)
4312   continue
c---source function for combined band 4
      do 4414 i=1,imax 
      do 4414 k=1,lp1
        vtmp3x(i,k)=source(ixo(i,k),4) 
        dsorcx(i,k)=dsrce(ixo(i,k),4) 
4414   continue
      do 4412 k=1,lp1
      do 4412 i=1,imax
      sorc(i,k,4)=vtmp3x(i,k)+dtemp(i,k)*dsorcx(i,k)
4412   continue
c---source function for combined band 5
      do 4514 i=1,imax 
      do 4514 k=1,lp1
        vtmp3x(i,k)=source(ixo(i,k),5) 
        dsorcx(i,k)=dsrce(ixo(i,k),5) 
4514   continue
      do 4512 k=1,lp1
      do 4512 i=1,imax
      sorc(i,k,5)=vtmp3x(i,k)+dtemp(i,k)*dsorcx(i,k)
4512   continue
c---source function for combined band 6
      do 4614 i=1,imax 
      do 4614 k=1,lp1
        vtmp3x(i,k)=source(ixo(i,k),6) 
        dsorcx(i,k)=dsrce(ixo(i,k),6) 
4614   continue
      do 4612 k=1,lp1
      do 4612 i=1,imax
      sorc(i,k,6)=vtmp3x(i,k)+dtemp(i,k)*dsorcx(i,k)
4612   continue
c---source function for combined band 7
      do 4714 i=1,imax 
      do 4714 k=1,lp1
        vtmp3x(i,k)=source(ixo(i,k),7) 
        dsorcx(i,k)=dsrce(ixo(i,k),7) 
4714   continue
      do 4712 k=1,lp1
      do 4712 i=1,imax
      sorc(i,k,7)=vtmp3x(i,k)+dtemp(i,k)*dsorcx(i,k)
4712   continue
c---source function for combined band 8
      do 4814 i=1,imax 
      do 4814 k=1,lp1
        vtmp3x(i,k)=source(ixo(i,k),8) 
        dsorcx(i,k)=dsrce(ixo(i,k),8) 
4814   continue
      do 4812 k=1,lp1
      do 4812 i=1,imax
      sorc(i,k,8)=vtmp3x(i,k)+dtemp(i,k)*dsorcx(i,k)
4812   continue
c---source function for band 9 (560-670 cm-1)
      do 4914 i=1,imax 
      do 4914 k=1,lp1
        vtmp3x(i,k)=source(ixo(i,k),9) 
        dsorcx(i,k)=dsrce(ixo(i,k),9) 
4914   continue
      do 4912 k=1,lp1
      do 4912 i=1,imax
      sorc(i,k,9)=vtmp3x(i,k)+dtemp(i,k)*dsorcx(i,k)
4912   continue
c---source function for band 10 (670-800 cm-1)
      do 5014 i=1,imax 
      do 5014 k=1,lp1
        vtmp3x(i,k)=source(ixo(i,k),10) 
        dsorcx(i,k)=dsrce(ixo(i,k),10) 
5014  continue
      do 5012 k=1,lp1
      do 5012 i=1,imax
      sorc(i,k,10)=vtmp3x(i,k)+dtemp(i,k)*dsorcx(i,k)
5012   continue
c---source function for band 11 (800-900 cm-1)
      do 5114 i=1,imax 
      do 5114 k=1,lp1
        vtmp3x(i,k)=source(ixo(i,k),11) 
        dsorcx(i,k)=dsrce(ixo(i,k),11) 
5114   continue
      do 5112 k=1,lp1
      do 5112 i=1,imax
      sorc(i,k,11)=vtmp3x(i,k)+dtemp(i,k)*dsorcx(i,k)
5112   continue
c---source function for band 12 (900-990 cm-1)
      do 5214 i=1,imax 
      do 5214 k=1,lp1
        vtmp3x(i,k)=source(ixo(i,k),12) 
        dsorcx(i,k)=dsrce(ixo(i,k),12) 
5214   continue
      do 5212 k=1,lp1
      do 5212 i=1,imax
      sorc(i,k,12)=vtmp3x(i,k)+dtemp(i,k)*dsorcx(i,k)
5212   continue
c---source function for band 13 (990-1070 cm-1)
      do 5314 i=1,imax 
      do 5314 k=1,lp1
        vtmp3x(i,k)=source(ixo(i,k),13) 
        dsorcx(i,k)=dsrce(ixo(i,k),13) 
5314   continue
      do 5312 k=1,lp1
      do 5312 i=1,imax
      sorc(i,k,13)=vtmp3x(i,k)+dtemp(i,k)*dsorcx(i,k)
5312   continue
c---source function for band 14 (1070-1200 cm-1)
      do 5414 i=1,imax 
      do 5414 k=1,lp1
        vtmp3x(i,k)=source(ixo(i,k),14) 
        dsorcx(i,k)=dsrce(ixo(i,k),14) 
5414   continue
      do 5412 k=1,lp1
      do 5412 i=1,imax
      sorc(i,k,14)=vtmp3x(i,k)+dtemp(i,k)*dsorcx(i,k)
5412   continue
!!!       print*, ' 5412 '
c 
c 
c        the following subroutine obtains nlte source function for co2
c 
c 
c---obtain special source functions for the 15 um band (csour),the
c   9.6 um band (osour) and the window region (ss1)
      do 131 i=1,imax*lp1
      ss1(i,1)=sorc(i,1,11)+sorc(i,1,12)+sorc(i,1,14)
131   continue
      do 143 i=1,imax*lp1
      csour(i,1)=sorc(i,1,9)+sorc(i,1,10)
      osour(i,1)=sorc(i,1,13) 
143   continue
c 
c 
c     second section produces 4 matrices for later use in program:  
c     1) avephi(i,j) is the scaled water mass for use in tables (this 
c     is u(p,p') in ref.(4);
c     2) over (i,j) is the water transmission function (using 
c     "emissivity" approximation) in the 560-800 cm-1  band;
c     3) to3(i,j) is the exact ozone transmission function (using 
c     parameters for the 990-1070 cm-1 band from the 1982 afgl catalog) 
c     4)emiss2(i,j) is the emissivity h20 transmission function due 
c     to the 10 um continuum,treated as one band. 
c 
      do 201 i=1,imax 
      totphi(i,1)=zero
      toto3(i,1)=zero 
      tphio3(i,1)=zero
      totvo2(i,1)=zero
201   continue
      do 203 k=2,lp1
      do 203 i=1,imax 
      totphi(i,k)=totphi(i,k-1)+var2(i,k-1) 
      toto3(i,k)=toto3(i,k-1)+var3(i,k-1) 
      tphio3(i,k)=tphio3(i,k-1)+var4(i,k-1) 
      totvo2(i,k)=totvo2(i,k-1)+cntval(i,k-1) 
203   continue

        if(ndi<0.and.nmaxpr==1.and.idjd<=imax.and.mydiag)then
         print *,'during fst88 myid',myid
         write(6,"('totphi ',10g10.3)")(totphi(idjd,k),k=1,kl+1)
         write(6,"('var2 ',10g10.3)")(var2(idjd,k),k=1,kl+1)
         write(6,"('toto3 ',10g10.3)")(toto3(idjd,k),k=1,kl+1)
         write(6,"('tphio3 ',10g10.3)")(tphio3(idjd,k),k=1,kl+1)
         write(6,"('totvo2 ',10g10.3)")(totvo2(idjd,k),k=1,kl)
         write(6,"('var1 ',10g10.3)")(var1(idjd,k),k=1,kl)
         write(6,"('var3 ',10g10.3)")(var3(idjd,k),k=1,kl)
         write(6,"('var4 ',10g10.3)")(var4(idjd,k),k=1,kl)
         write(6,"('delp2 ',10g10.3)")(delp2(idjd,k),k=1,kl)
         write(6,"('qh2o ',10g10.3)")(qh2o(idjd,k),k=1,kl)
        endif
 
c     the calculational procedure used here is: 
c       1) from the 1-d vectors (totphi,toto3,etc) construct a 1-d
c  arrays including all information for the upper (i>j) triangle
c  of the symmetric matrices(lp1,lp1);
c       2)perform computations on these arrays to obtain to31d,over1d,
c  avephj: upper triangle of the symmetric matrices in 1-d form 
c       3) fill up the lower triangle
c     the diagram below illustrates the relationship between the 1-d
c  array and the 2-d symmetric matrix for a 4x4 matrix. 
c 
c                    i
c             1      2       3       4
c           --------------------------
c        1           1       2       3     the nos. are the 
c    j   2                   4       5     positions in the 
c        3                           6     1-d array
c        4
c 
c 
c      compute "upper triangle" transmission functions for
c      the 9.6 um band (vtmp1) and the 15 um band (vtmp3). also,
c      the 
c      stage 1....compute o3 ,over transmission fctns and avephi
c---do k=1 calculation (from flux layer kk to the top) separately
c   as vectorization is improved,and ozone cts transmissivity
c   may be extracted here.
      do 302 i=1,imax*l
      avmo3(i,1)=toto3(i,2)
      avpho3(i,1)=tphio3(i,2)
      avephj(i,1)=totphi(i,2)
      avvo2(i,1)=totvo2(i,2)
      fac1(i,1)=bo3rnd(2)*avpho3(i,1)/avmo3(i,1)
      vtmp2(i,1)=haf*(fac1(i,1)*
     &    (sqrt(one+(four*ao3rnd(2)*avmo3(i,1))/fac1(i,1))-one))
      vtmp1(i,1)=exp(hm1ez*(vtmp2(i,1)+sko3r*avvo2(i,1)))
      vtmp3(i,1)=exp(hm1ez*(sqrt(ab15wd*avephj(i,1))+    
     &            skc1r*avvo2(i,1)))
c---the next 4 lines fill up the upper triangle of the to3,avephi
c   and over matrices, and all of the to3spc array. the notation
c   is obscure; note that an equivalent statement for line 1 would be:
c   for i=1 to imax and kp=k+1(=2) to lp1:  to3(i,kp,1)=vtmp1(i,kp)
      to3(i,2,1)=vtmp1(i,1)
      avephi(i,2,1)=avephj(i,1)
      over(i,2,1)=vtmp3(i,1)
      to3spc(i,1)=vtmp2(i,1)
302   continue

        if(ndi<0.and.nmaxpr==1.and.idjd<=imax.and.mydiag)then
         print *,'during_b fst88 myid',myid
         write(6,"('vtmp1 ',10g10.3)")(vtmp1(idjd,k),k=1,kl+1)
         do kk=1,kl+1
          print *,'kk = ',kk
          write(6,"('to3 ',10g10.3)") (to3(idjd,k,kk),k=1,kl+1)
         enddo
        endif

c---fill in lower triangle of to3,over arrays. it is unnecessary to
c   fill in avephi array at this time.
c      do 307 kp=2,lp1
c      do 307 i=1,imax
c      to3(i,1,kp)=vtmp1(i,kp-1)
c      over(i,1,kp)=vtmp3(i,kp-1)
c307   continue
cc---now repeat for the k=2..l cases. 
c      do 321 k=2,l
c      do 322 kk=1,lp1-k
c      do 322 i=1,imax 
c      avmo3(i,kk)=toto3(i,kk+k)-toto3(i,k)
c      avpho3(i,kk)=tphio3(i,kk+k)-tphio3(i,k) 
c      avephj(i,kk)=totphi(i,kk+k)-totphi(i,k) 
c      avvo2(i,kk)=totvo2(i,kk+k)-totvo2(i,k)
c322   continue
c      do 3221 i=1,imax*(lp1-k)
c      fac1(i,1)=bo3rnd(2)*avpho3(i,1)/avmo3(i,1)
c      vtmp2(i,1)=haf*(fac1(i,1)*
c     &    (sqrt(one+(four*ao3rnd(2)*avmo3(i,1))/fac1(i,1))-one))
c      vtmp1(i,1)=exp(hm1ez*(vtmp2(i,1)+sko3r*avvo2(i,1)))
c      vtmp3(i,1)=exp(hm1ez*(sqrt(ab15wd*avephj(i,1))+    
c     &            skc1r*avvo2(i,1)))
c      to3(i,k+1,k)=vtmp1(i,1)
c      avephi(i,k+1,k)=avephj(i,1)
c      over(i,k+1,k)=vtmp3(i,1)
c3221  continue
c      do 327 kp=k+1,lp1
c      do 327 i=1,imax
c      to3(i,k,kp)=vtmp1(i,kp-k)
c      over(i,k,kp)=vtmp3(i,kp-k)
c327   continue
c321   continue

c---fill in lower triangle of to3,over arrays. it is unnecessary to
c   fill in avephi array at this time.  Following is jlm re-doing of loops
      do kp=2,lp1
      do i=1,imax
        to3(i,1,kp)=vtmp1(i,kp-1)
        over(i,1,kp)=vtmp3(i,kp-1)
       enddo
      enddo
c---now repeat for the k=2..l cases. 
      do k=2,l
       do kk=1,lp1-k
        do i=1,imax 
         avmo3(i,kk)=toto3(i,kk+k)-toto3(i,k)
         avpho3(i,kk)=tphio3(i,kk+k)-tphio3(i,k) 
         avephj(i,kk)=totphi(i,kk+k)-totphi(i,k) 
         avvo2(i,kk)=totvo2(i,kk+k)-totvo2(i,k)
        enddo
       enddo
       do kx=1,lp1-k
        do i=1,imax
         fac1(i,kx)=bo3rnd(2)*avpho3(i,kx)/avmo3(i,kx)
         vtmp2(i,kx)=haf*(fac1(i,kx)*
     &          (sqrt(one+(four*ao3rnd(2)*avmo3(i,kx))/fac1(i,kx))-one))
         vtmp1(i,kx)=exp(hm1ez*(vtmp2(i,kx)+sko3r*avvo2(i,kx)))
         vtmp3(i,kx)=exp(hm1ez*(sqrt(ab15wd*avephj(i,kx))+    
     &               skc1r*avvo2(i,kx)))
         avephi(i,k+kx,k)=avephj(i,kx)
         to3(i,k+kx,k)=vtmp1(i,kx)
         over(i,k+kx,k)=vtmp3(i,kx)
         to3(i,k,k+kx)=vtmp1(i,kx)
         over(i,k,k+kx)=vtmp3(i,kx)
        enddo
       enddo  ! kx loop
c       do kp=k+1,lp1
c        do i=1,imax
c         to3(i,k,kp)=vtmp1(i,kp-k)
c         over(i,k,kp)=vtmp3(i,kp-k)
c        enddo
c       enddo
       if(ndi<0.and.nmaxpr==1.and.idjd<=imax.and.mydiag)then
        print *,'during_c fst88 k,myid',k,myid
        write(6,"('fac1 ',10g10.3)")(fac1(idjd,kx),kx=1,kl+1-k)
        write(6,"('avpho3 ',10g10.3)")(avpho3(idjd,kx),kx=1,kl+1-k)
        write(6,"('avmo3 ',10g10.3)")(avmo3(idjd,kx),kx=1,kl+1-k)
        write(6,"('avvo2 ',10g10.3)")(avvo2(idjd,kx),kx=1,kl+1-k)
        write(6,"('vtmp1 ',10g10.3)")(vtmp1(idjd,kx),kx=1,kl+1-k)
        write(6,"('vtmp2 ',10g10.3)")(vtmp2(idjd,kx),kx=1,kl+1-k)
        write(6,"('to3 ',10g10.3)") (to3(idjd,kx,k),kx=1,kl+1)
       endif
      enddo  ! k loop     

c---initialize diagonal elements of to3,over (not needed for avephi)
      do 309 k=1,lp1
      do 309 i=1,imax
      to3(i,k,k)=1.
      over(i,k,k)=1.
309   continue
c*****stage 4....compute one-band continuum transmissivities (emiss2)
c 
      do 481 i=1,imax*lp1
      vtmp3(i,1)=exp(hm1ez*totvo2(i,1))
      totevv(i,1)=one/vtmp3(i,1)
481   continue
c     do i=1,imax*lp1
c     if(vtmp3(i,1).eq.0.)print *,'ktau,i,totvo2 ',
c    .                             ktau,i,totvo2(i,1)
c     enddo
      do 501 k=1,lp1
      do 501 kp=1,k 
      do 503 i=1,imax 
      emiss2(i,kp,k)=vtmp3(i,k)*totevv(i,kp)
503   continue
501   continue
      do 505 k=1,l
      do 505 kp=k+1,lp1 
      do 507 i=1,imax 
      emiss2(i,kp,k)=vtmp3(i,kp)*totevv(i,k)
507   continue
505   continue
!!!      print*, ' 505 '
c         the third section calculates boundary layer and nearby layer
c     corrections to the transmission functions obtained above. methods
c     are given in ref. (4).
c       combine co21,over into co21; before making nbl corrections, 
c     load (1,k) values for use in exact cts calculations in spa88. 
      do 605 i=1,imax*lp1*lp1
      co21(i,1,1)=co21(i,1,1)*over(i,1,1)
605   continue
!!!      print*, ' 605 '
      do 607 k=1,lp1
      do 607 i=1,imax 
      co2sp(i,k)=co21(i,1,k)
607   continue
!!!      print*, ' 607 '
c          the following ratios are used in various nbl calculations: 
c
      do 619 k=1,l
      do 619 i=1,imax
      rlog(i,k)=log(over(i,k,k+1)*co2nbl(i,k))
619   continue
!!!      print*, ' 619 '
      do 601 i=1,imax*lm1
      delpr1(i,2)=delp(i,2)*(press(i,2)-p(i,2)) 
      alpsq1(i,2)=sqrt(delpr1(i,2)) 
      alp(i,lp1)=-alpsq1(i,2)*rlog(i,2)
601   continue
!!!      print*, ' 601 '
      do 603 i=1,imax*l
      delpr2(i,2)=delp(i,1)*(p(i,2)-press(i,1)) 
      alpsq2(i,2)=sqrt(delpr2(i,2)) 
      alp(i,1)=-alpsq2(i,2)*rlog(i,1)
603   continue
!!!      print*, ' 603 '
      do 625 i=1,imax 
      alp(i,ll)=-rlog(i,l)
      alp(i,llp1)=-rlog(i,l)*sqrt(delp(i,l)*(p(i,lp1)-press(i,lm1)))
625   continue
!!!      print*, ' 625 '
c        the first computation is for the 15 um band,with the  
c     for the combined h2o and co2 transmission function. 
c 
c       perform nbl computations for the 15 um band 
c***the statement function sf in prev. versions is now explicitly 
c   evaluated.
      do 631 i=1,imax*llp1
c!!!!!!!!!!
c     c(i,1)=alp(i,1)*(hmp66667+alp(i,1)*(quartr+alp(i,1)*hm6666m2))
      c(i,1)= 2.0/alp(i,1)**2 *
     &         (1.0-exp(-alp(i,1)) * (alp(i,1)+1.0)) - 1.0
631   continue
      do 641 i=1,imax 
      co21(i,lp1,lp1)=one+c(i,l)
      co21(i,lp1,l)=one+(delp2(i,l)*c(i,ll)-(press(i,l)-p(i,l))*
     & c(i,llm1))/(p(i,lp1)-press(i,l)) 
      co21(i,l,lp1)=one+((p(i,lp1)-press(i,lm1))*c(i,llp1)- 
     & (p(i,lp1)-press(i,l))*c(i,l))/(press(i,l)-press(i,lm1))
641   continue
c***  the k indices in the following loop run from 21 to 40 in the
c     l40 skyhi code version
      do 643 k=2,l
      do 643 i=1,imax 
      co21(i,k,k)=one+haf*(c(i,lm1+k)+c(i,k-1))
643   continue
c 
c    compute nearby-layer transmissivities for the o3 band and for the
c    one-band continuum band (to3 and emiss2). the sf2 function is
c    used. the method is the same as described for co2 in ref (4).
!     do 651 i=1,imax*lm1
!     csub(i,2)=cntval(i,2)*delpr1(i,2)
!     csub(i,lp1)=cntval(i,1)*delpr2(i,2)
!651  continue
      do 651 i=1,imax*lm1                    ! bug fixed error
      csub(i,2)=cntval(i,2)*delpr1(i,2)      ! bug fixed
651   continue                               ! bug fixed
      do i=1,imax*l                          ! bug fixed
      csub(i,lp1)=cntval(i,1)*delpr2(i,2)    ! bug fixed
      end do                                 ! bug fixed
c---the sf2 function in prev. versions is now explicitly evaluated
c!!!!!!!!!!!
c     do 655 i=1,imax*llm2
      do 655 i=1,imax*llm1
      csub2(i,2)=sko3r*csub(i,2)
      c(i,2)=csub(i,2)*(hmp5+csub(i,2)*(hp166666-csub(i,2)*h41666m2)) 
      c2(i,2)=csub2(i,2)*(hmp5+csub2(i,2)*
     &           (hp166666-csub2(i,2)*h41666m2)) 
655   continue
      do 661 i=1,imax 
c!!!!!!!!!!!
c     emiss2(i,lp1,lp1)=one+c(i,llm1) 
c     to3(i,lp1,lp1)=one+c2(i,llm1)
      emiss2(i,lp1,lp1)=one+c(i,ll) 
      to3(i,lp1,lp1)=one+c2(i,ll)
661   continue
      do 663 k=2,l
      do 663 i=1,imax 
      emiss2(i,k,k)=one+haf*(c(i,k)+c(i,lm1+k))
      to3(i,k,k)=one+haf*(c2(i,k)+c2(i,lm1+k))
663   continue
!!!      print*, '663 '
c 
c          fourth section obtains water transmission functions
c     used in q(approx) calculations and also makes nbl corrections:  
c     1) emiss (i,j) is the transmission function matrix obtained 
c     by calling subroutine e1e288; 
c     2) "nearby layer" corrections (emiss(i,i)) are obtained 
c     using subroutine e3v88; 
c     3) special values at the surface (emiss(l,lp1),emiss(lp1,l),
c     emiss(lp1,lp1)) are calculated. 
c 
c     emxx,av1,and empl are computed before avephi is modified
c 
c      obtain arguments for e1e288 and e3v88: 
c 
      do 801 i=1,imax 
      emx1(i)=qh2o(i,l)*press(i,l)*(press(i,l)-p(i,l))*gp0inv
801   continue
      do 803 k=1,lm1
      do 803 i=1,imax 
      emxx(i,k)=avephi(i,l,k)+emx1(i) 
803   continue
      do 811 i=1,imax*l
      empl(i,2)=qh2o(i,1)*p(i,2)*(p(i,2)-press(i,1))*gp0inv
811   continue
      do 812 i=1,imax*lm1
      empl(i,lp2)=qh2o(i,2)*p(i,2)*(press(i,2)-p(i,2))*gp0inv
812   continue
      do 821 i=1,imax 
      empl(i,1)=avephi(i,lp1,l) 
      empl(i,llp1)=empl(i,ll) 
      tpl(i,1)=temp(i,l)
      tpl(i,lp1)=haf*(t(i,lp1)+temp(i,l)) 
      tpl(i,llp1)=haf*(t(i,l)+temp(i,l))
821   continue
      do 823 k=2,l
      do 823 i=1,imax 
      tpl(i,k)=t(i,k) 
      tpl(i,k+l)=t(i,k) 
823   continue
      do 827 k=2,l
      do 829 i=1,imax 
      avephi(i,k,k)=emxx(i,k-1) 
829   continue
827   continue
!!!      print*, ' 827 '
c 
      do 833 i=1,imax 
      avephi(i,lp1,lp1)=avephi(i,lp1,l)+empl(i,l) 
c---a suitable value is assigned to avephi(i,1,1);it is unused by the
c   code and needed only to avoid indefinite values 
      avephi(i,1,1)=h101m16
833   continue
c     compute logs of water mass arguments for  e1e288; 
c     the corresponding quantity for e3v88 (empl) must be obtained
c     within that subroutine, as empl is used after e3v88 is called.
c
c---only the upper triangle values of avephi exist ;only these values
c   have their logarithms taken.
      do 841 k=1,lp1
      do 842 i=1,imax*(lp2-k)
      avephi(i,k,k)=log10(avephi(i,k,k))
842   continue
c---now fill in the lower triangle elements of avephi
      do 843 kp=k+1,lp1
      do 843 i=1,imax
      avephi(i,k,kp)=avephi(i,kp,k)
843   continue
841   continue
c
c     call e1e288 for emissivity transmission fctns for h2o 
!!!      print*, ' Calling e1e288 '
           call e1e288(e1cts1,e1cts2,e1flx,e1ctw1,e1ctw2,fxo,dtemp)
!!!      print*, ' Done e1e288 '
c 
      do 853 k=1,lm1
      do 855 i=1,imax 
      emiss(i,lp1,k)=haf*(emiss(i,k+1,k+1)+emiss(i,lp1,k))
855   continue
853   continue
c
c     call e3v88 for nbl h2o transmissivities 
!!!      print*, ' Calling e3v88 '
           call e3v88(emd,tpl,empl) 
!!!      print*, ' done e3v88 '
c 
c   compute nearby layer and special-case transmissivities for emiss
c    using methods for h2o given in ref. (4)
      do 851 k=2,l
      do 851 i=1,imax 
      emiss(i,k,k)=emd(i,k+l)+emd(i,k)
851   continue
c 
      do 861 i=1,imax 
      emiss(i,l,lp1)=(emd(i,1)*empl(i,1)-emd(i,lp1)*empl(i,lp1))/ 
     & emx1(i) + quartr*(emiss(i,l,lp1)+emiss(i,lp1,lp1)) 
      emiss(i,lp1,lp1)=two*emd(i,lp1) 
      emiss(i,lp1,l)=two*(emd(i,1)*empl(i,1)-emd(i,llp1)*empl(i,llp1))/
     & (qh2o(i,l)*press(i,l)*(p(i,lp1)-press(i,l))*gp0inv) 
861   continue
!!!      print*, ' 861 '
c 
c          subroutine spa88 is called to obtain exact cts for water 
c     co2 and o3, and approximate cts co2 and o3 calculations.
c 
!!!      print*, ' Calling spa88 '
      call spa88
!!!      print*, ' Done spa88 '
c 
c          this section performs the calculation of "emissivity"
c     fluxes for the 4 components comprising the lw frequency region. 
c     (emiss = the 0-160,1200-2200 cm-1 band; emiss2 the 800-990, 
c     1070-1200 cm-1 band; to3 the 990-1070 cm-1 band; co21 the 560-800 
c     cm-1 band).  emisst is the combined exchange term and flx the 
c     combined net flux.
c 
      do 901 i=1,imax*lp1
      tc(i,1)=(temp(i,1)*temp(i,1))**2
901   continue
      do 903 i=1,imax*l 
      oss(i,2)=osour(i,2)-osour(i,1)
      css(i,2)=csour(i,2)-csour(i,1)
      dtc(i,2)=tc(i,2)-tc(i,1)
      ss2(i,2)=ss1(i,2)-ss1(i,1)
903   continue
c
      do 905 k=1,lp1
      do 905 i=1,imax 
      emisst(i,1,k)=tc(i,1)*e1flx(i,k)+ss1(i,1)*
     & emiss2(i,k,1)+sorc(i,1,13)*to3(i,k,1)+csour(i,1)*co2sp(i,k)
905   continue
      do 909 k=1,lp1
      do 909 i=1,imax*l
      emisst(i,2,k)=dtc(i,2)*emiss(i,2,k)+
     &               ss2(i,2)*emiss2(i,2,k)+
     &               oss(i,2)*to3(i,2,k)  +css(i,2)*co21(i,2,k)
909   continue
c 
      do 912 k=1,lp1
      do 912 i=1,imax
      flx(i,k)=emisst(i,1,k)*cldfac(i,1,k)
      flxclr(i,k)=emisst(i,1,k)
912   continue
      do 911 kp=2,lp1
      do 911 k=1,lp1
      do 911 i=1,imax
      flx(i,k)=flx(i,k)+emisst(i,kp,k)*cldfac(i,kp,k)
      flxclr(i,k)=flxclr(i,k)+emisst(i,kp,k)
911   continue
c    this section computes the emissivity cts heating rates for 2 
c    emissivity bands: the 0-160,1200-2200 cm-1 band and the 800- 
c    990,1070-1200 cm-1 band. the remaining cts comtributions are 
c    contained in ctso3, computed in spa88. 
c 
      do 999 i=1,imax*lp1
      vtmp1(i,1)=emiss2(i,1,1)*cldfac(i,1,1)
      vtmp1clr(i,1)=emiss2(i,1,1)
999   continue
      do 1001 i=1,imax*l
      cts(i,1)=radcon*delp(i,1)*(tc(i,1)* 
     &     (e1ctw2(i,1)*cldfac(i,2,1)-e1ctw1(i,1)*cldfac(i,1,1)) +
     &      ss1(i,1)*(vtmp1(i,2)-vtmp1(i,1)))
      ctsclr(i,1)=radcon*delp(i,1)*(tc(i,1)*
     &     (e1ctw2(i,1)-e1ctw1(i,1)) +
     &      ss1(i,1)*(vtmp1clr(i,2)-vtmp1clr(i,1)))
1001  continue
c 
      do 1011 i=1,imax*l
      ceval(i,1)=tc(i,1)*(cldfac(i,1,1)*(e1cts1(i,1)-e1ctw1(i,1)) -
     &                    cldfac(i,2,1)*(e1cts2(i,1)-e1ctw2(i,1)))
      cevalclr(i,1)=tc(i,1)*((e1cts1(i,1)-e1ctw1(i,1)) -
     &                    (e1cts2(i,1)-e1ctw2(i,1)))
1011  continue
      do 1012 i=1,imax
      flx1e1(i)=tc(i,lp1)*cldfac(i,lp1,1)*
     &          (e1cts1(i,lp1)-e1ctw1(i,lp1))
      flx1e1clr(i)=tc(i,lp1)*(e1cts1(i,lp1)-e1ctw1(i,lp1))
1012  continue
      do 1014 k=1,l 
      do 1013 i=1,imax
      flx1e1(i)=flx1e1(i)+ceval(i,k)
      flx1e1clr(i)=flx1e1clr(i)+cevalclr(i,k)
1013  continue
1014  continue
c 
c     final section obtains emissivity heating rates, 
c     total heating rates and the flux at the ground
c 
c     ....calculate the emissivity heating rates 
      do 1101 i=1,imax*l
      heatem(i,1)=radcon*(flx(i,2)-flx(i,1))*delp(i,1)
      heatemclr(i,1)=radcon*(flxclr(i,2)-flxclr(i,1))*delp(i,1)
1101  continue
c     ....calculate the total heating rates
      do 1103 i=1,imax*l
      heatra(i,1)=heatem(i,1)-cts(i,1)-ctso3(i,1)+excts(i,1)
      heatraclr(i,1)=heatemclr(i,1)-ctsclr(i,1)-ctso3clr(i,1)
     &               +exctsclr(i,1)
1103  continue
c     ....calculate the flux at each flux level using the flux at the
c    top (flx1e1+gxcts) and the integral of the heating rates (vsum1) 
      do 1111 i=1,imax*l
      vsum1(i,1)=heatra(i,1)*delp2(i,1)*radcon1 
      vsum1clr(i,1)=heatraclr(i,1)*delp2(i,1)*radcon1
1111  continue
      do 1115 i=1,imax
      topflx(i)=flx1e1(i)+gxcts(i)
      flxnet(i,1)=topflx(i) 
      grnflxclr(i)=flx1e1clr(i)+gxctsclr(i)
1115  continue
c---only the surface value of flux (grnflx) is needed unless
c    the thick cloud section is invoked.
      do 1123 k=2,lp1 
      do 1123 i=1,imax
      flxnet(i,k)=flxnet(i,k-1)+vsum1(i,k-1)
      grnflxclr(i)=grnflxclr(i)+vsum1clr(i,k-1)
1123  continue
      do 1125 i=1,imax
      grnflx(i)=flxnet(i,lp1) 
1125  continue
c 
c     this is the thick cloud section.optionally,if thick cloud 
c     fluxes are to be "convectively adjusted",ie,df/dp is constant,
c     for cloudy part of grid point, the following code is executed.
c***first,count the number of clouds along the lat. row. skip the 
c   entire thick cloud computation if there are no clouds.
      icnt=0
      do 1301 i=1,imax
      icnt=icnt+nclds(i)
1301  continue
      if (icnt.eq.0) go to 6999
c---find the maximum number of clouds in the latitude row
      kclds=nclds(1)
      do 2106 i=2,imax
      kclds=max(nclds(i),kclds)
2106  continue
c 
c 
c***obtain the pressures and fluxes of the top and bottom of
c   the nc'th cloud (it is assumed that all ktop and kbtm's have
c   been defined!). 
      do 1361 kk=1,kclds
      do 1363 i=1,imax
      j1=ktop(i,kk+1)
c     if (j1.eq.1) go to 1363
      j3=kbtm(i,kk+1)+1
        if ((j3-1).gt.j1) then
           ptop(i)=p(i,j1) 
           pbot(i)=p(i,j3)
           ftop(i)=flxnet(i,j1)
           fbot(i)=flxnet(i,j3)
c***obtain the "flux derivative" df/dp (delptc) 
           delptc(i)=(ftop(i)-fbot(i))/(ptop(i)-pbot(i))
c***calculate the tot. flux chg. from the top of the cloud, for 
c   all levels.
           do 1365 k=j1+1,j3-1
           z1(i,k)=(p(i,k)-ptop(i))*delptc(i)+ftop(i)
           flxnet(i,k)=flxnet(i,k)*(one-camt(i,kk+1)*emcld(i,kk+1)) +
     &            z1(i,k)*camt(i,kk+1)*emcld(i,kk+1)
1365       continue
        endif
1363  continue
1361  continue
c***using this flux chg. in the cloudy part of the grid box, obtain 
c   the new fluxes, weighting the clear and cloudy fluxes:again, only 
c    the fluxes in thick-cloud levels will eventually be used.
c     do 6051 k=1,lp1 
c     do 6051 i=1,imax
c     flxnet(i,k)=flxnet(i,k)*(one-camt(i,nc)*emcld(i,nc)) +
c    1            z1(i,k)*camt(i,nc)*emcld(i,nc)
c051  continue
c***merge flxthk into flxnet for appropriate levels. 
c     do 1401 k=1,lp1
c     do 1401 i=1,imax
c     if (k.gt.itop(i) .and. k.le.ibot(i)
c    1  .and.  (nc-1).le.nclds(i))  then
c          flxnet(i,k)=flxthk(i,k)
c     endif
c401  continue
c 
c******end of cloud loop***** 
6001  continue
6999  continue
c***the final step is to recompute the heating rates based on the 
c   revised fluxes: 
      do 6101 i=1,imax*l
c     heatra(i,1)=radcon*(flxnet(i,2)-flxnet(i,1))*delp(i,1)
      heatra(i,1)=(flxnet(i,2)-flxnet(i,1))
6101  continue
c     the thick cloud section ends here.

        if(ndi<0.and.nmaxpr==1.and.idjd<=imax.and.mydiag)then
         print *,'At end of fst88 myid',myid
         write(6,"('grnflx ',g10.3)") grnflx(idjd)
         write(6,"('rlog ',10g10.3)")(rlog(idjd,k),k=1,kl)
         write(6,"('emxx ',10g10.3)")(emxx(idjd,k),k=1,kl)
         write(6,"('heatra ',10g10.3)")(heatra(idjd,k),k=1,kl)
         write(6,"('heatem ',10g10.3)")(heatem(idjd,k),k=1,kl+1)
         write(6,"('flxnet ',10g10.3)")(flxnet(idjd,k),k=1,kl+1)
         write(6,"('vsum1 ',10g10.3)")(vsum1(idjd,k),k=1,kl+1)
         write(6,"('z1 ',10g10.3)")(z1(idjd,k),k=1,kl+1)
         write(6,"('flxclr ',10g10.3)")(flxclr(idjd,k),k=1,kl+1)
         write(6,"('flx ',10g10.3)")(flx(idjd,k),k=1,kl+1)
         write(6,"('tpl ',10g10.3)")(tpl(idjd,k),k=1,2*kl+1)
         write(6,"('empl ',10g10.3)")(empl(idjd,k),k=1,2*kl+1)
         write(6,"('csub ',10g10.3)")(csub(idjd,k),k=1,2*kl+1)
         write(6,"('c ',10g10.3)")(c(idjd,k),k=1,2*kl+1)
         write(6,"('c2 ',10g10.3)")(c2(idjd,k),k=1,2*kl+1)
         write(6,"('alp ',10g10.3)")(alp(idjd,k),k=1,2*kl+1)
         write(6,"('co2sp ',10g10.3)")(co2sp(idjd,k),k=1,kl+1)
         write(6,"('totevv ',10g10.3)")(totevv(idjd,k),k=1,kl+1)
         write(6,"('vtmp1 ',10g10.3)")(vtmp1(idjd,k),k=1,kl+1)
         write(6,"('vtmp2 ',10g10.3)")(vtmp2(idjd,k),k=1,kl+1)
         write(6,"('vtmp3 ',10g10.3)")(vtmp3(idjd,k),k=1,kl+1)
         write(6,"('fac1 ',10g10.3)")(fac1(idjd,k),k=1,kl+1)
         write(6,"('avmo3 ',10g10.3)")(avmo3(idjd,k),k=1,kl+1)
         write(6,"('avpho3 ',10g10.3)")(avpho3(idjd,k),k=1,kl+1)
         write(6,"('avvo2 ',10g10.3)")(avvo2(idjd,k),k=1,kl+1)
         write(6,"('csour ',10g10.3)")(csour(idjd,k),k=1,kl+1)
         write(6,"('osour ',10g10.3)")(osour(idjd,k),k=1,kl+1)
         write(6,"('ss1 ',10g10.3)")(ss1(idjd,k),k=1,kl+1)
         write(6,"('ss2 ',10g10.3)")(ss2(idjd,k),k=1,kl+1)
         write(6,"('totevv ',10g10.3)")(totevv(idjd,k),k=1,kl+1)
         write(6,"('dtemp ',10g10.3)")(dtemp(idjd,k),k=1,kl+1)
         write(6,"('dtc ',10g10.3)")(dtc(idjd,k),k=1,kl+1)
         write(6,"('oss ',10g10.3)")(oss(idjd,k),k=1,kl+1)
         write(6,"('css ',10g10.3)")(css(idjd,k),k=1,kl+1)
         do kk=1,kl+1
          print *,'kk = ',kk
          write(6,"('emiss ',10g10.3)") (emiss(idjd,k,kk),k=1,kl+1)
          write(6,"('emiss2 ',10g10.3)") (emiss2(idjd,k,kk),k=1,kl+1)
          write(6,"('emisst ',10g10.3)") (emisst(idjd,k,kk),k=1,kl+1)
          write(6,"('to3 ',10g10.3)") (to3(idjd,k,kk),k=1,kl+1)
          write(6,"('co21 ',10g10.3)") (co21(idjd,k,kk),k=1,kl+1)
         enddo
        endif
      return
      end 
