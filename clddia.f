! Conformal Cubic Atmospheric Model
    
! Copyright 2015 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
! This file is part of the Conformal Cubic Atmospheric Model (CCAM)
!
! CCAM is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! CCAM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with CCAM.  If not, see <http://www.gnu.org/licenses/>.

!------------------------------------------------------------------------------
      
      subroutine clddia(rhum,sigf,cf,ktd,kbd,icld,i,j) ! jlm
      use arrays_m  ! for t
      use cc_mpi, only : myid
      use const_phys
      use map_m
      use morepbl_m ! for condc
      use newmpar_m
      use parm_m
      use pbl_m     ! for tss
      use sigs_m    ! for sig
      use soil_m    ! provides land()
      use vvel_m    ! sdot
      parameter (maprhcrt=0)  ! 1 for grid size dependent rhcrit
      parameter (nconv_cld=0) ! 1,2,3 for convective enhancement of clouds
      parameter (ntest=0)     ! 0 to turn off test prints   jd ~ every 6th j
                              ! i.e. jd-1 is a multiple of 6
!     This version for both DARLAM and C-C model
!     This one only has 3 layers of clouds, determined only by rhum
!     with jlm tuning options    N.B. k goes from top down
!     and jlm fix-ups for cloud fraction & cloud levels
!     N.B. C-C uses imax loop, so may not have every diag. jd value here

      include 'kuocom.h'    ! for sigcll,nclddia,nstab_cld,nrhcrit
c     parameter (nclddia=0)   ! conversion of cld to cloudiness, 0 for original
c     parameter (nclddia=5)   ! conversion of cld to cloudiness, 0 for original
c     parameter (nstab_cld=0) ! 0 for original; 3 for stability-enhanced cll
c     parameter (nstab_cld=3) ! 0 for original; 3 for stability-enhanced cll
c     parameter (nrhcrit=0)   ! 0 for Hal's original
c     parameter (nrhcrit=7)   ! 7 for jlm, 8+ for jlm power-scheme
      real cldlev(2,3)
      real rhum(kl),sigf(kl),rhcrit(kl)
      dimension cf(kl),ktd(kl),kbd(kl),icld(2,3),clmx(3)
c     data clmx/.85,.9,1./   ! HML 20,12,8  for nrhcrit=9
c     data clmx/.6,.6,.9/    ! HML 8,4,6/8  for nrhcrit=8
c     data clmx/.8,.7,1./    ! HML 8,8,6/8  for nrhcrit=8
c     data clmx/.8,.7,.8/    ! HML 16,8,4/8 for nrhcrit=8
c     data clmx/.8,.7,.8/    ! HML 20,8,3/6 for nrhcrit=8
c     data clmx/.8,.5,.8/    ! HML 20,3,2/4 for nrhcrit=8
c     data clmx/.5,.5,.6/    ! HML 8,3,1/2  for nrhcrit=8
c     data clmx/.5,.55,.8/   ! HML 8,4,2/3  for nrhcrit=8
c     data clmx/.5,.7,.7/    ! HML 8,6,2/3  for nrhcrit=8
c     data clmx/.8,.8,.7/    ! HML 12,8,2/2  for nrhcrit=8  r29
c     data clmx/.8,.7,.6/    ! HML 12,8,1.5  for nrhcrit=8  r30
c     data clmx/.8,.8,.6/    ! HML 12,8,1.5  for nrhcrit=8  a31
c     data clmx/.8,.8,.6/    ! HML 12,6,1.5  for nrhcrit=8  a32
c     data clmx/.8,.9,.6/    ! HML 12,8,1.5  for nrhcrit=8  a33
c     data clmx/.8,.9,.7/    ! HML 12,8,sqrt+**4  for nrhcrit=8  a34
c     data clmx/.8,.55,.4/   ! HML 12,4,x+x**8  for nrhcrit=8  a35
c     data clmx/.8,.55,.7/   ! HML 8,4,.4*x+.3*x**8  for nrhcrit=8  a87 to b14
c     data clmx/.8,.55,.7/   ! HML 8,4,.3*x+.4*x**8  for nrhcrit=8  b15
      data clmx/.8,.6,.7/    ! HML 6,.1&.5,.3*x+.4*x**8  for nrhcrit=8  b16
c     data clmx/.8,.6,.7/    ! HML 8,.1&.5,.3*x+.4*x**8  for nrhcrit=8  e31
      logical start

c     define the cloud type limits in terms of sigma
      data cldlev /.15, .43, .43, .8, .8, 1.0/   ! Hal's
c     data cldlev /.15, .43, .43, .8, .8, .96/   ! jlm's - does almost nothing
      data start /.true./
      save cldlev,start,ktop,ztop

c     ---------------------------------------------------------------------*
c     start computation                                                   *
c     ---------------------------------------------------------------------*

      iq=i+(j-1)*il
      if(start)then
        cldlev(2,3)=sigcll
         if (myid==0) then
         write (6,"('cldlev:',6f8.3)") cldlev
c       ---------------------------------------------------------------------*
c       on first calculation calculate the high, middle and low              *
c       cloud level limits                                                   *
c       ---------------------------------------------------------------------*
        print *,'nrhcrit,nclddia,nstab_cld,maprhcrt,nconv_cld,sigcll: ',
     .           nrhcrit,nclddia,nstab_cld,maprhcrt,nconv_cld,sigcll
        if(nrhcrit.gt.7)print *,'with clm HML = ',clmx
        if(maprhcrt.ne.0)print *,'N.B. grid-dependent maprhcrt used'
        end if
        if(nrhcrit.eq.0.and.nstab_cld.gt.0)stop 'invalid nstab_cld'

        k1 = 1
        do nc = 1,3
c         top down
          do k = k1,kl-1
c           find cloud top sig; sigf is sigma, starting at top level
            if ( sigf(k).ge.cldlev(1,nc) ) then
c              set cloud top
               if ( icld(1,nc).eq.0         ) icld(1,nc) = k
c              set cloud base
               if ( sigf(k).lt.cldlev(2,nc) ) icld(2,nc) = k
            endif
          end do ! k=k1,kl-1
          k1 = icld(2,nc) + 1
          if (myid==0) print *,'nc,icld_top,_bott; then k values: ',
     .        nc,icld(1,nc),icld(2,nc),  kl+1-icld(2,nc),kl+1-icld(1,nc)
        end do ! nc=1,3
!       constants for jlm stability enhancement of low cloud
!       if(nstab_cld.eq.1)ktop=kl+1-icld(1,3)                 ! proper k of cll top
        ktop=kl+1-(icld(1,3)+icld(2,3))/2   ! middle of cll  nstab_cld=2
        ztop=rdry*300.*(1.-sig(ktop))/(grav*1000.) ! in km

c       calculate & print critical relative humidity for land then sea
        if(nrhcrit.eq.7)then
          if(myid==0)print *,'nrhcrit=7, so rhcrit done point by point'
c         ----- land --------------------------------------------------------*
          do k=icld(1,1),icld(2,1) ! high clds; k here goes from top dn
           rhcrit(k)=cldh_lnd
          enddo
          do k=icld(1,2),icld(2,2) ! middle clds; k here goes from top dn
           rhcrit(k)=cldm_lnd
          enddo
          do k=icld(1,3),icld(2,3)   ! low clouds     ; k here goes from top dn
           rhcrit(k)=cldl_lnd
          enddo
          if (myid==0) then
             print *,'rhcrit over land:'
             print 91,(rhcrit(k),k=kl,1,-1)
          end if
 91       format(f5.1,' L',4f5.1,' M',5f5.1,' H',4f5.1,' |',4f5.1)
!         the above format is suitable for 18-level runs
c         ----- sea  --------------------------------------------------------*
          do k=icld(1,1),icld(2,1) ! high clds; k here goes from top dn
           rhcrit(k)=cldh_sea
          enddo
          do k=icld(1,2),icld(2,2) ! middle clds; k here goes from top dn
           rhcrit(k)=cldm_sea
          enddo
          do k=icld(1,3),icld(2,3)   ! low clouds     ; k here goes from top dn
           rhcrit(k)=cldl_sea
          enddo
          if (myid==0) then
             print *,'rhcrit over sea:'
             print 91,(rhcrit(k),k=kl,1,-1)
          end if
        elseif(nrhcrit.eq.0)then
          do k=1,kl
           eta=sigf(k)
           rhcrit(k)=(1.-eta*(1.-eta)*(.268+3.464*eta))*100.        ! Hal's
          enddo
          if (myid==0) then
             print *,'rhcrit for old nrhcrit=0 scheme:'
             print 91,(rhcrit(k),k=kl,1,-1)
          end if
        endif  ! (nrhcrit.eq.7)
        start = .false.
      endif   ! (start)
c     ---------------------------------------------------------------------*

      if(nrhcrit.eq.7)then   ! needs land or sea depending on iq
c      calculate critical relative humidity                           *
       if(maprhcrt.eq.1)then  
!       new section with Rhcrit varying linearly between 95. and cldh/m/l
!       for grid resolution varying between 5 km and 100 km
        dx=.001*ds/em(iq)   ! in km
        factmap=max( 0. , min( (dx-5.)/(100.-5.) , 1. ) )
        fact95=(1.-factmap)*95.
        if(land(iq))then    ! over land
          do k=icld(1,1),icld(2,1) ! high clds; k here goes from top dn
           rhcrit(k)=fact95+factmap*cldh_lnd
          enddo
          do k=icld(1,2),icld(2,2) ! middle clds; k here goes from top dn
           rhcrit(k)=fact95+factmap*cldm_lnd
          enddo
          do k=icld(1,3),icld(2,3)   ! low clouds  ; k here goes from top down
           rhcrit(k)=fact95+factmap*cldl_lnd
          enddo
        else                  ! over sea  (do above as well!)
          do k=icld(1,1),icld(2,1) ! high clds; k here goes from top dn
           rhcrit(k)=fact95+factmap*cldh_sea
          enddo
          do k=icld(1,2),icld(2,2) ! middle clds; k here goes from top dn
           rhcrit(k)=fact95+factmap*cldm_sea
          enddo
          do k=icld(1,3),icld(2,3)   ! low clouds  ; k here goes from top dn
           rhcrit(k)=fact95+factmap*cldl_sea
          enddo
        endif  ! (land(iq))  .. else ..
        if(ntest.eq.1.and.j.eq.jd.and.id-1.eq.mod(i-1,il)
     &       .and.myid==0)then
         print *,'ds,em,dx,factmap ',ds,em(iq),dx,factmap
         print *,'rhcrit ',rhcrit
        endif
       else    ! i.e. maprhcrt=0
        if(land(iq))then    ! over land
          do k=icld(1,1),icld(2,1) ! high clds; k here goes from top dn
           rhcrit(k)=cldh_lnd
          enddo
          do k=icld(1,2),icld(2,2) ! middle clds; k here goes from top dn
           rhcrit(k)=cldm_lnd
          enddo
          do k=icld(1,3),icld(2,3)   ! low clouds  ; k here goes from top down
           rhcrit(k)=cldl_lnd
          enddo
        else                  ! over sea  (do above as well!)
          do k=icld(1,1),icld(2,1) ! high clds; k here goes from top dn
           rhcrit(k)=cldh_sea
          enddo
          do k=icld(1,2),icld(2,2) ! middle clds; k here goes from top dn
           rhcrit(k)=cldm_sea
          enddo
          do k=icld(1,3),icld(2,3)   ! low clouds  ; k here goes from top dn
           rhcrit(k)=cldl_sea
          enddo
        endif  ! (land(iq))  .. else ..
       endif   ! (maprhcrt.eq.1) .. else ..
      endif    ! (nrhcrit.eq.7)

      if(nstab_cld.lt.0)then  ! enchanced cll for cold surface
!       use t diff between surface and second layer
        tdiff=tss(iq)-t(iq,2)
!       for surface 5 deg colder (nstab_cld=-5), get 25% effect 
        rhsub=min(max(25.*tdiff/nstab_cld,0.),50.)       ! between 0 and 50%
        do k=icld(1,3),icld(2,3)
         rhcrit(k)=min(rhcrit(k)-rhsub,100.) ! not more than 100% (not needed)
        enddo
        if(ntest.eq.1.and.j.eq.jd.and.id-1.eq.mod(i-1,il)
     &       .and.myid==0)then
          print *,'i,nstab_cld,rhsub ',i,nstab_cld,rhsub
          print *,'land,tss,t ',land(iq),tss(iq),t(iq,ktop)
        endif
      endif  ! (nstab_cld.lt.0)

      if(nstab_cld.gt.0)then
!       use t diff between surface and top of lowest layer
!       ktop=kl+1-(icld(1,3)+icld(2,3))/2           ! proper k of middle of cll
!       ztop=r*300.*(1.-sig(ktop))/(grav*1000.) ! in km
        tdiff=tss(iq)-t(iq,ktop)
!       if(nstab_cld.eq.1)stab=max((tdiff/ztop-6.5)/6. ,0.) ! 6 deg/km for 10% eff.
!       if(nstab_cld.eq.2)stab=max((tdiff/ztop-6.5)/1.5,0.) ! 3 deg/km for 20% eff.
        stab=max((tdiff/ztop)/6.5,0.) ! 6.5 deg/km for 10% eff.  nstab_cld=3
!       if(nstab_cld.eq.4)stab=(tdiff/ztop)/6.5      ! 6.5 deg/km for 10% effect
        rhsub=min(stab*10.,20.)              ! at most 20%
        do k=icld(1,3),icld(2,3)
         rhcrit(k)=min(rhcrit(k)-rhsub,100.) ! not more than 100%
        enddo
        if(ntest.eq.1.and.j.eq.jd.and.id-1.eq.mod(i-1,il)
     &       .and.myid==0)then
          print *,'i,nstab_cld,ktop,ztop ',i,nstab_cld,ktop,ztop
          print *,'land,tss,t ',land(iq),tss(iq),t(iq,ktop)
          print *,'stab,rhsub ',stab,rhsub
        endif
      endif  ! (nstab_cld.gt.0)

c     ---------------------------------------------------------------------*
c                                                                          *
c     cloud parameter calculation
c                                                                          *
c     ---------------------------------------------------------------------*
      if(ntest.ge.1.and.j.eq.jd.and.id-1.eq.mod(i-1,il)
     &     .and.myid==0)then
        print *,'rhcrit:'
        print 91,(rhcrit(k),k=kl,1,-1)
        print 92,(k,k=1,kl)
!92     format(i3,3i5,2x,4i5,2x,5i5,2x,4i5,2x,i5)
 92     format(i3,2x,4i5,2x,5i5,2x,4i5,2x,4i5)
        print *,'rhum for i= ',i,' j= ',j,' ktau= ',ktau
!       print 91,rhum
        print 91,(rhum(k),k=kl,1,-1)
        print *,'qg for i= ',i,' j= ',j,' ktau= ',ktau
        print 91,(1000.*qg(iq,k),k=1,kl)
        print *,'T celsius for i= ',i,' j= ',j,' ktau= ',ktau
        print 91,(t(iq,k)-273.16,k=1,kl)
c       print *,'rtt cooling rate (deg/day) from prev call radrive '
c       print 91,(86400.*rtt(iq,k),k=1,kl)
        print *,'approx. adiabatic cooling rate (deg/day) '
        print 91,(-86400.*rdry*t(iq,k)*
     .   .5*(sdot(iq,k)+sdot(iq,k+1))*dsig(k)/(dt*cp*sig(k)),k=1,kl-1)
        print *,'sdot*100 '
        print 91,(100.*sdot(iq,k),k=1,kl)
      endif

      if(nrhcrit.eq.8)then   ! power-scheme
        nc=1       ! high
        cf(nc)=0.  ! cloud frac = 0.
        do k=icld(1,nc),icld(2,nc)
         cld=max(0.,min(.01*rhum(k),1.))
c        cf(nc)=max(cf(nc),clmx(nc)*cld**8)    ! high   
         cf(nc)=max(cf(nc),clmx(nc)*cld**6)    ! high   
        enddo ! k=icld(1,nc),icld(2,nc)
        ktd(nc)=icld(1,nc)   ! k value for cloud top
        kbd(nc)=icld(2,nc)   ! k value for cloud bottom
        nc=2       ! middle
        cf(nc)=0.  ! cloud frac = 0.
        do k=icld(1,nc),icld(2,nc)
         cld=max(0.,min(.01*rhum(k),1.))
c        cf(nc)=max(cf(nc),clmx(nc)*cld**4)     ! mid
         cf(nc)=max(cf(nc),.1*cld+.5*cld**4)    ! mid*** 
        enddo ! k=icld(1,nc),icld(2,nc)
        ktd(nc)=icld(1,nc)   ! k value for cloud top
        kbd(nc)=icld(2,nc)   ! k value for cloud bottom
        nc=3       ! low
        cf(nc)=0.  ! cloud frac = 0.
        do k=icld(1,nc),icld(2,nc)
         cld=max(0.,min(.01*rhum(k),1.))
c        cf(nc)=max(cf(nc),.4*cld+.3*cld**8 )   ! low*** 
         cf(nc)=max(cf(nc),.3*cld+.4*cld**8 )   ! low*** 
        enddo ! k=icld(1,nc),icld(2,nc)
        ktd(nc)=icld(1,nc)   ! k value for cloud top
        kbd(nc)=icld(2,nc)   ! k value for cloud bottom
      endif   !   (nrhcrit.eq.8) 

      if(nrhcrit.eq.9)then   ! test power-scheme
        nc=1       ! high
        cf(nc)=0.  ! cloud frac = 0.
        do k=icld(1,nc),icld(2,nc)
         cld=max(0.,min(.01*rhum(k),1.))
!        cf(nc)=max(cf(nc),.1*cld+.7*cld**24 )      ! high***  
         cf(nc)=max(cf(nc),.05*cld+.75*cld**32 )    ! high***  
        enddo ! k=icld(1,nc),icld(2,nc)
        ktd(nc)=icld(1,nc)   ! k value for cloud top
        kbd(nc)=icld(2,nc)   ! k value for cloud bottom
        nc=2       ! middle
        cf(nc)=0.  ! cloud frac = 0.
        do k=icld(1,nc),icld(2,nc)
         cld=max(0.,min(.01*rhum(k),1.))
         cf(nc)=max(cf(nc),.1*cld+.7*cld**16 )    ! mid*** 
        enddo ! k=icld(1,nc),icld(2,nc)
        ktd(nc)=icld(1,nc)   ! k value for cloud top
        kbd(nc)=icld(2,nc)   ! k value for cloud bottom
        nc=3       ! low
        cf(nc)=0.  ! cloud frac = 0.
        do k=icld(1,nc),icld(2,nc)
         cld=max(0.,min(.01*rhum(k),1.))
         cf(nc)=max(cf(nc),.1*cld+.7*cld**8 )   ! low*** 
        enddo ! k=icld(1,nc),icld(2,nc)
        ktd(nc)=icld(1,nc)   ! k value for cloud top
        kbd(nc)=icld(2,nc)   ! k value for cloud bottom
      endif   !   (nrhcrit.eq.9) 

      if(nrhcrit.eq.10)then   ! newer power-scheme
        nc=1       ! high
        cf(nc)=0.  ! cloud frac = 0.
        do k=icld(1,nc),icld(2,nc)
         cld=max(0.,min(.01*rhum(k),1.))
         cf(nc)=max(cf(nc),.05*cld+.85*cld**24 )    ! high***  
        enddo ! k=icld(1,nc),icld(2,nc)
        ktd(nc)=icld(1,nc)   ! k value for cloud top
        kbd(nc)=icld(2,nc)   ! k value for cloud bottom
        nc=2       ! middle
        cf(nc)=0.  ! cloud frac = 0.
        do k=icld(1,nc),icld(2,nc)
         cld=max(0.,min(.01*rhum(k),1.))
         cf(nc)=max(cf(nc),.15*cld+.55*cld**16 )    ! mid*** 
        enddo ! k=icld(1,nc),icld(2,nc)
        ktd(nc)=icld(1,nc)   ! k value for cloud top
        kbd(nc)=icld(2,nc)   ! k value for cloud bottom
        nc=3       ! low
        cf(nc)=0.  ! cloud frac = 0.
        do k=icld(1,nc),icld(2,nc)
         cld=max(0.,min(.01*rhum(k),1.))
         cf(nc)=max(cf(nc),.3*cld+.4*cld**16 )   ! low*** 
        enddo ! k=icld(1,nc),icld(2,nc)
        ktd(nc)=icld(1,nc)   ! k value for cloud top
        kbd(nc)=icld(2,nc)   ! k value for cloud bottom
      endif   !   (nrhcrit.eq.10) 

      if(nrhcrit.lt.8)then  ! non-power-scheme
        do nc=1,3
         cf(nc)=0.  ! cloud frac = 0.
c        level indexes to check for cloud
         iswitch=0
         do k=icld(1,nc),icld(2,nc)
c         does it meet criteria for cloud?
          if(rhum(k).gt.rhcrit(k)) then
c           calc. cloud fraction (max for layer)
            cld=(rhum(k)-rhcrit(k))/(100.-rhcrit(k))
            if(nclddia.eq.1)then
              cf(nc)=max( cf(nc),cld*cld*cld )   ! jlm cube formula
            elseif(nclddia.eq.2)then   ! x weighting of x**2 & x**4
              cf(nc)=max( cf(nc),cld*cld*max(0.,1.-cld+cld**3) )  ! jlm
            elseif(nclddia.eq.3)then   ! x**2 weighting of x & x**4
              cf(nc)=max( cf(nc),cld*abs(1.-cld**2+cld**5) )  ! jlm
            elseif(nclddia.eq.4)then   ! x**3 weighting of x & x**3
              cf(nc)=max( cf(nc),(1.-cld**3)*cld+cld**3*cld**3 )  ! jlm
            elseif(nclddia.eq.5)then   ! x**3 weighting of x & x**4
              cf(nc)=max( cf(nc),(1.-cld**3)*cld+cld**3*cld**4 )  ! jlm
            else   ! i.e. nclddia=0
              cf(nc)=max( cf(nc),abs(cld)*cld )   ! usual
            endif
            if(iswitch.eq.0)then
              ktd(nc)=k    ! k value for cloud top
              iswitch=1
            endif  ! (iswitch.eq.0)
            kbd(nc)=k
            if(ntest.eq.2.and.j.eq.jd.and.id-1.eq.mod(i-1,il)
     &           .and.myid==0)
     .      print *,'rh_max>rhcrit with k,cld,cloud_frac= '
     .                                 ,k,cld,cf(nc)
          endif  ! ( rhum(k).gt.rhcrit(k) )
         enddo ! k=icld(1,nc),icld(2,nc)
        enddo ! nc=1,3
      endif   !   (nrhcrit.lt.8) 
      
      if(nconv_cld.eq.1.or.nconv_cld.eq.2)then  ! conv enh. of l & m clouds
        if(condc(iq).gt.0.)then
          rainrt=condc(iq)*86400./dt  ! in mm per day
          if(nconv_cld.eq.1)convcld=min(.2+.05*sqrt(rainrt),.8)   ! modified
          if(nconv_cld.eq.2)convcld=min(.1+.07*log(1.+rainrt),.8) ! NCAR total
          convlayr=1.-sqrt(1.-convcld)
c         print *,'i,j,rainrt,convlayr ',i,j,rainrt,convlayr
          cf(3)=1.-(1.-cf(3))*(1.-convlayr)  ! low
          cf(2)=1.-(1.-cf(2))*(1.-convlayr)  ! middle
        endif   ! (condc(iq).gt.0.)
      endif     ! (nconv_cld.eq.1.or.nconv_cld.eq.2)
      if(nconv_cld.eq.3)then  ! convective enhancement of l m & h clouds
        if(condc(iq).gt.0.)then
          rainrt=condc(iq)*86400./dt  ! in mm per day
          convcld=min(.2+.07*log(1.+rainrt),.8) ! NCAR total
!         convlayr=1.-cbrt(1.-convcld)
          convlayr=1.-    (1.-convcld)**(1./3.)
          cf(3)=1.-(1.-cf(3))*(1.-convlayr)  ! low
          cf(2)=1.-(1.-cf(2))*(1.-convlayr)  ! middle
          cf(1)=1.-(1.-cf(3))*(1.-convlayr)  ! high
        endif   ! (condc(iq).gt.0.)
      endif     ! (nconv_cld.eq.3)

      if(ntest.ge.1.and.j.eq.jd.and.id-1.eq.mod(i-1,il)
     &     .and.myid==0)then
        print *,'cf  LMH ',(cf(nc),nc=3,1,-1)
        print *,'ktd LMH ',(ktd(nc),nc=3,1,-1)
        print *,'kbd LMH ',(kbd(nc),nc=3,1,-1)
      endif

      return
      end
