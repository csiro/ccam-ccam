      subroutine lwr88

CDIR$ TASK COMMON KDACOM
CDIR$ TASK COMMON RADISW
CDIR$ TASK COMMON TFCOM
c CDIR$ TASK COMMON VTEMP

c     subroutine lwr88 computes temperature-corrected co2 transmission
c   functions and also computes the pressure grid and layer optical 
c   paths.
c          inputs:                (common blocks) 
c      press,temp,rh2o,qo3             radisw 
c      co251,co258,cdt51,cdt58         co2bd3 
c      c2d51,c2d58,co2m51,co2m58       co2bd3 
c      cdtm51,cdtm58,c2dm51,c2dm58     co2bd3 
c      stemp,gtemp                     co2bd3 
c      co231,co238,cdt31,cdt38         co2bd2 
c      c2d31,c2d38                     co2bd2 
c      co271,co278,cdt71,cdt78         co2bd4 
c      c2d71,c2d78                     co2bd4 
c      betinw                          bdwide 
c          outputs: 
c      co21,co2nbl,co2sp1,co2sp2       tfcom
c      var1,var2,var3,var4,cntval      kdacom 
c      qh2o,p,delp,delp2,t             kdacom 
c          called by: 
c      radmn or input routine of model
c          calls: 
c      fst88
c 
      use co2dta_m
      use kdacom_m
      use newmpar_m
      use parm_m
      use radisw_m
      use tfcom_m
      use work3lwr_m

      include 'hcon.h'
      include 'rdparm.h'
      include 'rnddta.h'
c 
      dimension diftd(imax,lp1,lp1) 
      dimension dco2dt(imax,lp1,lp1),d2cdt2(imax,lp1,lp1) 
      dimension texpsl(imax,lp1),tlsqu(imax,lp1)
      dimension vsum4(imax,l) 
      dimension dift1d(imax,lp1m) 
c****compute flux pressures (p) and differences (delp2,delp)
c****compute flux level temperatures (t) and continuum temperature
c    corrections (texpsl) 
c    changed to use supplied half level pressure and temperatures 
c    m dix 17/8/90
!      print*, ' In lwr88 '
!      print*, ' In lwr qo3 1', (qo3(1,k),k=1,l)
      do 103 k=1,lp1
      do 103 i=1,imax 
      p(i,k)=press2(i,k)
      t(i,k)=temp2(i,k)
103   continue
      do 107 k=1,l
      do 107 i=1,imax 
      delp2(i,k)=p(i,k+1)-p(i,k)
      delp(i,k)=one/delp2(i,k)
107   continue
!      print*, ' 107 '
!      print*, ' In lwr qo3 1', (qo3(1,k),k=1,l)
c****compute argument for cont.temp.coeff.
c    (this is 1800.(1./temp-1./296.)) 
!      print*, ' TEMP lwr88 ', (temp(1,k),k=1,lp1)
      do 125 k=1,lp1
      do 125 i=1,imax 
      texpsl(i,k)=h18e3/temp(i,k)-h6p08108
c...then take exponential 
      texpsl(i,k)=exp(texpsl(i,k))
125   continue
!      print*, ' 125 '
!      print*, ' In lwr qo3 1', (qo3(1,k),k=1,l)

c***compute optical paths for h2o and o3, using the diffusivity 
c   approximation for the angular integration (1.66). obtain the
c   unweighted values(var1,var3) and the weighted values(var2,var4).
c   the quantities h3m4(.0003) and h3m3(.003) appearing in the var2 and 
c   var4 expressions are the approximate voigt corrections for h2o and
c   o3,respectively.
c 
      do 131 k=1,l
      do 131 i=1,imax 
      qh2o(i,k)=rh2o(i,k)*diffctr 
c---vv is the layer-mean pressure (in atm),which is not the same as 
c   the level pressure (press)
      vv(i,k)=haf*(p(i,k+1)+p(i,k))*p0inv 
      var1(i,k)=delp2(i,k)*qh2o(i,k)*ginv
      var3(i,k)=delp2(i,k)*qo3(i,k)*diffctr*ginv
      var2(i,k)=var1(i,k)*(vv(i,k)+h3m4)
      var4(i,k)=var3(i,k)*(vv(i,k)+h3m3)
c  compute optical path for the h2o continuum, using roberts coeffs.
c  (betinw),and temp. correction (texpsl). the diffusivity factor 
c  (which cancels out in this expression) is assumed to be 1.66. the
c  use of the diffusivity factor has been shown to be a significant 
c  source of error in the continuum calcs.,but the time penalty of
c  an angular integration is severe.
c 
      cntval(i,k)=texpsl(i,k)*rh2o(i,k)*var2(i,k)*betinw/
     &             (rh2o(i,k)+rath2omw)
131   continue
c     do k=1,l
c      do i=1,imax
c      if(cntval(i,k)>1.e7)then
c        print *,'ktau ',ktau
c       print *,'i,k,temp,texpsl,rh2o,var2,betinw,cntval ',i,k,
c    .  temp(i,k),texpsl(i,k),rh2o(i,k),var2(i,k),betinw,cntval(i,k)
c      endif
c      enddo
c      enddo
!      print*, ' 131 '
!      print*, ' In lwr qo3 1', (qo3(1,k),k=1,l)
!      print*, ' In lwr qo3 2', (qo3(2,k),k=1,l)
!      print*, ' In lwr var3 1', (var3(1,k),k=1,l)
!      print*, ' In lwr var3 1', (var3(2,k),k=1,l)
c***compute weighted temperature (tdav) and pressure (tstdav) integrals 
c   for use in obtaining temp. difference bet. sounding and std.
c   temp. sounding (dift) 
      do 161 i=1,imax 
      tstdav(i,1)=zero
      tdav(i,1)=zero
161   continue
      do 162 k=1,lp1
      do 162 i=1,imax 
      vsum3(i,k)=temp(i,k)-stemp(k) 
162   continue
      do 163 k=1,l
      do 165 i=1,imax 
      vsum2(i)=gtemp(k)*delp2(i,k)
      vsum1(i)=vsum2(i)*vsum3(i,k)
      tstdav(i,k+1)=tstdav(i,k)+vsum2(i)
      tdav(i,k+1)=tdav(i,k)+vsum1(i)
165   continue
163   continue
c***compute dift
      do 204 k=1,l
      do 204 kp=k+1,lp1
      do 204 i=1,imax
      dift(i,kp,k)=(tdav(i,kp)-tdav(i,k))/
     &              (tstdav(i,kp)-tstdav(i,k))
204   continue
!      print*, ' 204 '
      do 205 i=1,imax
      dift(i,1,1)=0.
205   continue
      do 206 k=1,l
      do 206 i=1,imax 
      dift(i,k+1,k+1)=haf*(vsum3(i,k+1)+vsum3(i,k))
206   continue
      do 207 k=2,lp1
      do 207 kp=1,k-1
      do 207 i=1,imax
      dift(i,kp,k)=dift(i,k,kp)
207   continue
c 
c****evaluate coefficients for co2 pressure interpolation (a1,a2) 
      do 171 i=1,imax 
      a1(i)=(press(i,lp1)-p0xzp8)/p0xzp2
      a2(i)=(p0-press(i,lp1))/p0xzp2
171   continue
!      print*, ' A1 ', a1(1), press(1,lp1), p0xzp8, p0xzp2
!      print*, ' 171 '
c***perform co2 pressure interpolation on all inputted transmission 
c   functions and temp. derivatives 
c---successively computing co2r,dco2dt and d2cdt2 is done to save 
c   storage (at a slight loss in computation time)
      do 184 k=1,lp1
      do 184 i=1,imax 
        co2r1(i,k)=a1(i)*co231(k)+a2(i)*co238(k)
        d2cd21(i,k)=h1m3*(a1(i)*c2d31(k)+a2(i)*c2d38(k))
        dco2d1(i,k)=h1m2*(a1(i)*cdt31(k)+a2(i)*cdt38(k))
        co2r2(i,k)=a1(i)*co271(k)+a2(i)*co278(k)
        d2cd22(i,k)=h1m3*(a1(i)*c2d71(k)+a2(i)*c2d78(k))
        dco2d2(i,k)=h1m2*(a1(i)*cdt71(k)+a2(i)*cdt78(k))
184   continue
      do 190 k=1,l
      do 190 i=1,imax 
        co2mr(i,k)=a1(i)*co2m51(k)+a2(i)*co2m58(k)
        co2md(i,k)=h1m2*(a1(i)*cdtm51(k)+a2(i)*cdtm58(k))
        co2m2d(i,k)=h1m3*(a1(i)*c2dm51(k)+a2(i)*c2dm58(k))
190   continue
!      print*, ' CO2MR ', co2mr(1,1), a1(1), a2(1), co2m51(1), co2m58(1)
c***compute co2 temperature interpolations for all bands,using dift 
      do 240 k=1,lp1
      do 240 kp=1,lp1 
      do 240 i=1,imax 
      co2r(i,kp,k)=a1(i)*co251(kp,k)+a2(i)*co258(kp,k)
      dco2dt(i,kp,k)=h1m2*(a1(i)*cdt51(kp,k)+a2(i)*cdt58(kp,k))
      d2cdt2(i,kp,k)=h1m3*(a1(i)*c2d51(kp,k)+a2(i)*c2d58(kp,k))
      co21(i,kp,k)=co2r(i,kp,k)+dift(i,kp,k)*(dco2dt(i,kp,k)+
     &             haf*dift(i,kp,k)*d2cdt2(i,kp,k))
240   continue
c***compute transmission fctns used in spa88
c---(in the 250 loop,dift really should be (i,1,k), but dift is 
c    invariant with respect to k,kp,and so (i,1,k)=(i,k,1)) 
      do 250 k=1,lp1
      do 250 i=1,imax 
      co2sp1(i,k)=co2r1(i,k)+dift(i,k,1)*(dco2d1(i,k)+haf*dift(i,k,1)*
     & d2cd21(i,k)) 
      co2sp2(i,k)=co2r2(i,k)+dift(i,k,1)*(dco2d2(i,k)+haf*dift(i,k,1)*
     & d2cd22(i,k)) 
250   continue
c--- we aren't doing nbl tfs on the 100 cm-1 bands .
      do 260 k=1,l
      do 260 i=1,imax 
      co2nbl(i,k)=co2mr(i,k)+dift(i,k,k+1)*(co2md(i,k)+haf* 
     & dift(i,k,k+1)*co2m2d(i,k)) 
260   continue
!      print*, ' CO2NBL ', co2mr(1,1),dift(1,1,2),co2md(1,1),co2m2d(1,1)
c***compute temp. coefficient based on t(k) (see ref.2) 
      do 264 k=1,lp1
      do 264 i=1,imax
      if (t(i,k).le.h25e2) then
         tlsqu(i,k)=b0+(t(i,k)-h25e2)*
     &                      (b1+(t(i,k)-h25e2)*
     &                   (b2+b3*(t(i,k)-h25e2))) 
      else 
         tlsqu(i,k)=b0 
      endif
264   continue
c***apply to all co2 tfs
      do 280 k=1,lp1
      do 282 kp=1,lp1 
      do 282 i=1,imax 
      co21(i,kp,k)=co21(i,kp,k)*(one-tlsqu(i,kp))+tlsqu(i,kp) 
282   continue
280   continue
      do 284 k=1,lp1
      do 286 i=1,imax 
      co2sp1(i,k)=co2sp1(i,k)*(one-tlsqu(i,1))+tlsqu(i,1) 
      co2sp2(i,k)=co2sp2(i,k)*(one-tlsqu(i,1))+tlsqu(i,1) 
286   continue
284   continue
!      print*, ' Before 288 CO2NBL ', co2nbl(1,1)
      do 288 k=1,l
      do 290 i=1,imax 
      co2nbl(i,k)=co2nbl(i,k)*(one-tlsqu(i,k))+tlsqu(i,k) 
290   continue
288   continue
!      print*, ' After 288 CO2NBL ', co2nbl(1,1)
!      print*, ' Calling fst88 '
      call fst88
!      print*, ' lwr88 done '
      return
      end
 
