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
      
      subroutine bett_cuc ( dt, pt, kpnt, oshal )  ! was called cucnvc
c***  efiln was undefined  - changed to efimn 28/9/00
      use betts1_m
      use newmpar_m
      parameter (ntest=0)  ! replaces debug; set to 1 or 2 for degugging
c     with constants altered by bfr
c     this is part of the Betts-Miller parameterization
c     ******************************************************************
c     *                                                                *
c     *  convective adjustment for deep or shallow convection          *
c     *                                                                *
c     *  references:                                                   *
c     *                                                                *
c     *  betts, a.k., 1986:  a new convective adjustment scheme.       *
c     *    part i: observational and theoretical basis.  quart. j. r.  *
c     *    met. soc., 112, 677-691.                                    *
c     *                                                                *
c     *  betts, a.k., and m.j. miller, 1986:  a new convective         *
c     *    adjustment scheme.  part ii: single column tests using      *
c     *    gate wave, bomex, atex and arctic air mass data sets.       *
c     *    quart. j. r. met. soc., 112, 693-709.                       *
c     *                                                                *
c     *                                                                *
c     *  n.b.  part of the code is scalar.  in global models           *
c     *  convection occurs in less than 30% of the points.  with       *
c     *  simulataneous vector processing for both deep and shallow     *
c     *  convection, there would be a lot of redundant vector          *
c     *  computations.  if vector processing is 10 times faster        *
c     *  than scalar, one might hope that the cpu time will be about   *
c     *  the same for both scalar and vector code.                     *
c     *                                                                *
c     ******************************************************************
c
c ??? warning: this subroutine will not work if kl.lt.12;
c
c***********************************************************************
c
c needs: dtq2, klh, res, pd, t, htm, q, ape
c        rdth, qs0, sqs, rdq, ptbl, aeta, pt, rdp, the0, sthe
c        rdthe, ttbl, hbm2, thl, prec, cuprec, sm
c
      parameter ( d608=.608, dm2859=-.2858964143e0
     &, a2=17.2693882e0, t0=273.16e0, a4=35.86e0
     &, pq0=379.90516e0, tfr=274.16e0, epsq=2.e-12
     &, cp=1004.6e0, elwv=2.50e6, elivw=2.72e6, row=1.e3, g=9.8e0
     &, stabs=.80, stabd=.90,stabfc=.90, bitp=1.1
     &, trel=7200., dttop=1., dthvt=.30 )
       parameter ( lstab=1 )
       parameter (
     .  dsps=-3000., dspc=-3000.
     ., dspbl=-4000.e0, dsp0l=-6000.e0, dsptl=-2000.e0
     ., dspbs=-4000.e0, dsp0s=-6000.e0, dspts=-2000.e0
     ., dspbf=-4000.e0, dsp0f=-6000.e0, dsptf=-2000.e0
     ., epsntp=.0010e0, efifc=5.0e0, efimn=dspbs/dspbf
     ., avgrfi=(efimn+1.e0)*.5e0
     ., slopb=0., slop0=0., slopt=0.)
       parameter (a23m4l=a2*(t0-a4)*elwv
     &, elocp=elivw/cp, cprlg=cp/(row*g*elwv), rcp=1./cp )
c
c     include 'imjmkl.h'     ! now in betts1.h
c
      parameter ( ksmud=0 )  ! this option actually removed by jlm
      parameter ( lsh=4, kbm=3 )
!      parameter ( lscrch=4*kl+1+la
!     &, l1=la+1, l2=la+kl+1, l3=la+2*kl+1, l4=la+3*kl+1)
c
c-----------------------------------------------------------------------
c     include 'eta.h'        ! now in betts1.h
c-----------------------------------------------------------------------
c     include 'masks.h'      ! now in betts1.h
c-----------------------------------------------------------------------
c     include 'phys.h'       ! now in betts1.h
c-----------------------------------------------------------------------
c     include 'scrch.h'    !  uses common/work2/   now
       dimension trefk(kl),qrefk(kl),pk(kl),apek(kl),tk(kl)
     &,thsk  (kl),psk   (kl),apesk (kl),qk    (kl),therk (kl)
     &,thvref(kl),thvmod(kl),dift  (kl),difq  (kl)
     .,ntopd (kl),nbotd (kl),ntops (kl),nbots (kl)
     .,ntopn (kl),nbotn (kl),ndpthd(kl),ndpths(kl),ndpthn(kl)

c     real htop(ifull),hbot(ifull) not used

      real ltop(ifull),lbot(ifull),ittb(ifull),iqtb(ifull),
     . pdsl(ifull),tbt(ifull),qbt(ifull),apebt(ifull),tth(ifull),
     . pp(ifull),bqs00(ifull),sqs00(ifull),bqs10(ifull),sqs10(ifull),
     . bq(ifull),sq(ifull),qq(ifull),psp(ifull)

      real thbt(ifull),thesp(ifull),tdif(ifull)

      real ape(ifull,kl),tref(ifull,kl)  ! 28/9/00
      real qmod(ifull,kl)                ! 28/9/00
      real tmod(ifull,kl)                ! 28/9/00

      logical ofirst,oshal

!     following are all equivalenced to work2 arrays (see below)
      real iptb(ifull),ithtb(ifull),
     . p(ifull),botchk(ifull),
     . tthbt(ifull),tpsp(ifull),apesp(ifull),tthes(ifull),tq(ifull)
      real tp(ifull),p00(ifull),bthe00(ifull),t00(ifull),p10(ifull),
     . sthe00(ifull),t10(ifull),p01(ifull),bthe10(ifull),t01(ifull),
     . p11(ifull),sthe10(ifull),t11(ifull),bth(ifull),sth(ifull)

!      equivalence (ittb  ,iptb  ,ql)
!      equivalence (iqtb  ,ithtb ,tl   )
!      equivalence (tdif  ,tne   )
!      equivalence (thbt  ,qne    )
!      equivalence (tbt   ,p     ,tse  )
!      equivalence (psp   ,qse)
!      equivalence (qbt   ,botchk)
!      equivalence (apebt ,tthbt)
!      equivalence (pp    ,tpsp)
!      equivalence (cqq   ,apesp ,tthes )
!      equivalence (tth   ,tq    ,tp )
!      equivalence (bqs00 ,p00   ,bthe00,t00)
!      equivalence (sqs00 ,p10   ,sthe00,t10  )
!      equivalence (bqs10 ,p01   ,bthe10,t01)
!      equivalence (sqs10 ,p11   ,sthe10,t11  )
!      equivalence (bq    ,bth )
!      equivalence (sq    ,sth  )


      data ofirst/.true./
c
      save
c
c--------------preparations---------------------------------------------
c
      if(ntest.ne.0)print*,'entering bett_cuc:il,jl,kl= ',il,jl,kl
      if ( ofirst ) then
c
         ofirst=.false.
         call bettinit ( dt, pt )  ! was initcuc
c
         do 106 iq=1,ifull
           cldefi(iq)=avgrfi
!          hbot(iq)=0.
!          htop(iq)=0.
 106     continue
c
c--------------tile constants evaluated here for clarity----------------
c
         etrel= dtq2/trel
         print *,'etrel=',etrel
c
      endif
c
c----------initialize constants-----------------------------------------
c
      nshal=0
      ndeep=0
      nneg=0
c
      do 105 l=1,kl
        ntopd(l)=0
        nbotd(l)=0
        ntops(l)=0
        nbots(l)=0
        ntopn(l)=0
        nbotn(l)=0
        ndpths(l)=0
        ndpthd(l)=0
        ndpthn(l)=0
 105  continue
c
c----------preparations-------------------------------------------------
c
      avps=0.
      do 100 iq=1,ifull
        lbot (iq)=klh(iq)
        thesp(iq)=0.
        pdsl(iq)=res(iq)*pd(iq)
        avps=avps+pdsl(iq)
        tref (iq,1)=t(iq,1)
 100  continue
      avps=avps/ifull
c
c--------------calculating ape -----------------------------------------
c--------------padding specific humidity if too small-------------------
c
      lqm=1
      do 110 l=1,kl
        do 120 iq=1,ifull
          apest=pdsl(iq)*aeta(l)+pt
          ape(iq,l)=(apest/100000.)**dm2859
          if ( q(iq,l).lt.epsq ) q(iq,l)=epsq*htm(iq,l)
 120    continue
        play=avps*aeta(l)+pt
        if ( play.lt.100.e2 ) lqm=l
c       if ( ntest.gt.0 ) print *,l,ape(kpnt,l),apest,aeta(l),pdsl(kpnt)
 110  continue
      if ( ntest.gt.0 ) print *,'cucnvc:avps,kbm,lqm,lsh ',
     .   avps,kbm,lqm,lsh
c
c--------------search for maxilum buoyancy level------------------------
c
      do 190 kb=1,kbm
c
        if ( ntest.gt.0 ) print *,'cucnvc: kb ',kb
c
c--------------trial maxilum buoyancy level variables-------------------
c
!       do iq=2*il+1,ifull-2*il-2  ! DARLAM
        do iq=1,ifull
          lbtm    =klh(iq)-kb+1
          tbt(iq)  =t(iq,lbtm)
          qbt(iq)  =q(iq,lbtm)
          apebt(iq)=ape(iq,lbtm)
        enddo

c      if ( ntest.gt.0 ) then
c          print *,'cucnvc:130 ',lbtm,tbt(kpnt),qbt(kpnt),apebt(kpnt)
c       endif

c--------------scaling potential temperature & table index--------------
c
!       do iq=2*il+1,ifull-2*il-2  ! DARLAM
        do iq=1,ifull
         tthbt(iq)=tbt(iq)*apebt(iq)
         tth  (iq)=(tthbt(iq)-thl)*rdth
         qq   (iq)=tth(iq)-aint(tth(iq))
         ittb (iq)=int(tth(iq))+1
        enddo
c
c--------------keeping indices within the table-------------------------
c
!       do iq=2*il+1,ifull-2*il-2  ! DARLAM
         do iq=1,ifull
         if ( ittb(iq).lt.1 ) then
c          if ( ntest.gt.0 .and. iq.eq.kpnt ) then
c             print *,'cucnvc:a 145 ittb,qq',ittb(iq),qq(iq)
c          endif
           ittb(iq)=1
           qq  (iq)=0.
         endif
         if ( ittb(iq).ge.jtb ) then
c          if ( ntest.gt.0 .and. iq.eq.kpnt ) then
c             print *,'cucnvc:b 145 ittb,qq ',ittb(iq),qq(iq)
c          endif
           ittb(iq)=jtb-1
           qq  (iq)=0.
         endif
        enddo  ! 145

c       if ( ntest.gt.0 ) then
c          print *,'cucnvc:145 tthbt,tth,qq,ittb',
c    .     tthbt(kpnt),tth(kpnt),qq(kpnt),ittb(kpnt)
c       endif

c--------------base and scaling factor for spec. humidity---------------
c
!      do iq=2*il+1,ifull-2*il-2  ! DARLAM
       do iq=1,ifull
         ittbk=nint(ittb(iq))
         bqs00(iq)=qs0(ittbk)
         sqs00(iq)=sqs(ittbk)
         bqs10(iq)=qs0(ittbk+1)
         sqs10(iq)=sqs(ittbk+1)
       enddo  ! 150

c--------------scaling spec. humidity & table index---------------------
!      do iq=2*il+1,ifull-2*il-2  ! DARLAM
       do iq=1,ifull
        bq  (iq)=(bqs10(iq)-bqs00(iq))*qq(iq)+bqs00(iq)
        sq  (iq)=(sqs10(iq)-sqs00(iq))*qq(iq)+sqs00(iq)
        tq  (iq)=(qbt(iq)-bq(iq))/sq(iq)*rdq
        pp  (iq)=tq(iq)-aint(tq(iq))
        iqtb(iq)=int(tq(iq))+1
       enddo  

c--------------keeping indices within the table-------------------------
!      do iq=2*il+1,ifull-2*il-2  ! DARLAM
       do iq=1,ifull
        if ( iqtb(iq).lt.1 ) then
c          if ( ntest.gt.0 .and. iq.eq.kpnt ) then
c             print *,'cucnvc: a 165 iqtb,pp',iqtb(iq),pp(iq)
c          endif
           iqtb(iq)=1
           pp  (iq)=0.
        end if
        if ( iqtb(iq).ge.itb ) then
c          if ( ntest.gt.0 .and. iq.eq.kpnt ) then
c             print *,'cucnvc: b 165 iqtb,pp',iqtb(iq),pp(iq)
c          endif
           iqtb(iq)=itb-1
           pp  (iq)=0.
        end if
       enddo  

c--------------saturation pressure at four surrounding table pts.-------
!      do iq=2*il+1,ifull-2*il-2  ! DARLAM
       do iq=1,ifull
        iqbetts=nint(iqtb(iq))
        it=nint(ittb(iq))
        p00(iq)=ptbl(iqbetts  ,it  )
        p10(iq)=ptbl(iqbetts+1,it  )
        p01(iq)=ptbl(iqbetts  ,it+1)
        p11(iq)=ptbl(iqbetts+1,it+1)
       enddo  
c
c--------------saturation point variables at the bottom-----------------
c
!      do iq=2*il+1,ifull-2*il-2  ! DARLAM
       do iq=1,ifull
        tpsp(iq)=p00(iq)+(p10(iq)-p00(iq))*pp(iq)
     .          +(p01(iq)-p00(iq))*qq(iq)
     2          +(p00(iq)-p10(iq)-p01(iq)+p11(iq))*pp(iq)*qq(iq)
        apesp(iq)=(tpsp(iq)/100000.)**dm2859
        tthes(iq)=tthbt(iq)*exp(elocp*qbt(iq)*apesp(iq)/tthbt(iq))
       enddo  !  180

c       if ( ntest.gt.0 ) then
c          print *,'cucnvc: 180 ',tpsp(kpnt),apesp(kpnt),tthes(kpnt)
c       endif

c--------------check for maxilum buoyancy-------------------------------
c
!      do iq=2*il+1,ifull-2*il-2  ! DARLAM
       do iq=1,ifull
        if ( tthes(iq).gt.thesp(iq) ) then
           psp  (iq)=tpsp(iq)
           thbt (iq)=tthbt(iq)
           thesp(iq)=tthes(iq)
        end if
       enddo  !

c       if ( ntest.gt.0 ) then
c          print *,'cucnvc: 185 ',psp(kpnt),thbt(kpnt),thesp(kpnt)
c       endif

c-----------------------------------------------------------------------
c
c---------end search for maxilum buoyancy level-------------------------
c
 190  continue
c
c-----------choose cloud base as model level just below psp-------------
c
!      do iq=2*il+1,ifull-2*il-2  ! DARLAM
       do iq=1,ifull
        botchk(iq)=0.
       enddo  !
c
      do 196 ll=1,kl
        l=kl+1-ll
c
!      do iq=2*il+1,ifull-2*il-2  ! DARLAM
       do iq=1,ifull
          klhm=klh(iq)-1
          p(iq)=pdsl(iq)*aeta(l)+pt
          if ( p(iq).lt.psp(iq) .and. botchk(iq).lt..5 ) then
             lbot(iq)=min(l+1,klhm)
             botchk(iq)=1.
          endif
       enddo  !

c       if ( ntest.gt.0 ) then
c          print *,'cucnvc: 196 ',l,lbot(kpnt),
c    .     botchk(kpnt),p(kpnt),psp(kpnt)
c       endif

 196  continue
c
c--------------cloud top computation------------------------------------
c
c set qmod to 0 to indicate cloud top not found
c
      do 205 iq=1,ifull
       qmod(iq,1)=0.
       tdif(iq)=0.
 205   ltop(iq)=lbot(iq)
c
c------begin level loop--------------------------------------------
c
c start at cloud bottom and work upward
c
      do 200 ivi=1,kl
c       print*,'cucnvc: ivi= ',ivi
        l=kl+1-ivi
c
c--------------scaling pressure & tt table index------------------------
c
!      do iq=2*il+1,ifull-2*il-2  ! DARLAM
       do iq=1,ifull
        p (iq)=pdsl(iq)*aeta(l)+pt
        tp(iq)=(p(iq)-pl)*rdp
        qq(iq)=tp(iq)-aint(tp(iq))
        iptb(iq)=int(tp(iq))+1
       enddo  !
c
c--------------keeping indices within the table-------------------------
c
!      do iq=2*il+1,ifull-2*il-2  ! DARLAM
       do iq=1,ifull
        if ( iptb(iq).lt.1 ) then
c          if ( ntest.gt.0 .and.iq.eq.kpnt ) then
c             print *,'cucnvc: a 215 ',iptb(iq),qq(iq)
c          endif
           iptb(iq)=1
           qq  (iq)=0.
        end if
        if ( iptb(iq).ge.itb ) then
c          if ( ntest.gt.0 .and. iq.eq.kpnt ) then
c             print *,'cucnvc: b 215 ',iptb(iq),qq(iq)
c          endif
           iptb(iq)=itb-1
           qq  (iq)=0.
        end if
       enddo  !
c
c--------------base and scaling factor for the--------------------------
!      do iq=2*il+1,ifull-2*il-2  ! DARLAM
       do iq=1,ifull
        iptbk=nint(iptb(iq))
        bthe00(iq)=the0(iptbk)
        sthe00(iq)=sthe(iptbk)
        bthe10(iq)=the0(iptbk+1)
        sthe10(iq)=sthe(iptbk+1)
       enddo  !
c
c--------------scaling the & tt table index-----------------------------
!      do iq=2*il+1,ifull-2*il-2  ! DARLAM
       do iq=1,ifull
        bth  (iq)=(bthe10(iq)-bthe00(iq))*qq(iq)+bthe00(iq)
        sth  (iq)=(sthe10(iq)-sthe00(iq))*qq(iq)+sthe00(iq)
        tth  (iq)=(thesp(iq)-bth(iq))/sth(iq)*rdthe
        pp   (iq)=tth(iq)-aint(tth(iq))
        ithtb(iq)=int(tth(iq))+1
       enddo  !
c
c--------------keeping indices within the table-------------------------
!      do iq=2*il+1,ifull-2*il-2  ! DARLAM
       do iq=1,ifull
        if ( ithtb(iq).lt.1 ) then
c          if ( ntest.gt.0.and.iq.eq.kpnt) then
c             print *,'cucnvc: a 235 ',ithtb(iq),pp(iq)
c          endif
           ithtb(iq)=1
           pp  (iq)=0.
        end if
        if ( ithtb(iq).ge.jtb ) then
c          if ( ntest.gt.0.and.iq.eq.kpnt) then
c             print *,'cucnvc: b 235 ',ithtb(iq),pp(iq)
c          endif
           ithtb(iq)=jtb-1
           pp  (iq)=0.
        end if
       enddo  !
c
c--------------temperature at four surrounding tt table pts.------------
!      do iq=2*il+1,ifull-2*il-2  ! DARLAM
       do iq=1,ifull
        ith=nint(ithtb(iq))
        ip =nint(iptb (iq))
        t00(iq)=ttbl(ith  ,ip  )
        t10(iq)=ttbl(ith+1,ip  )
        t01(iq)=ttbl(ith  ,ip+1)
        t11(iq)=ttbl(ith+1,ip+1)
       enddo  !
c
c--------------parcel temperature---------------------------------------
!      do iq=2*il+1,ifull-2*il-2  ! DARLAM
       do iq=1,ifull
       tref(iq,l)=(t00(iq)+(t10(iq)-t00(iq))*pp(iq)
     .           +(t01(iq)-t00(iq))*qq(iq)+(t00(iq)-t10(iq)-t01(iq)
     2           +t11(iq))*pp(iq)*qq(iq))*htm(iq,l)
       enddo  !
c
c--------------buoyancy check-------------------------------------------
c
      ldif=l+lstab
c     if ( ntest.gt.0 ) print *,l,ldif,qmod(kpnt,1)
c
!      do iq=2*il+1,ifull-2*il-2  ! DARLAM
       do iq=1,ifull
        tmod(iq,l)=ltop(iq)
        if ( tref(iq,l).gt.t(iq,l)-dttop ) then
           ltop(iq)=l
        endif
        tmod(iq,l)=ltop(iq)
c
c test to see if ltop is same as level below
c
c       ntdif = nint(tmod(iq,ldif)-tmod(iq,l))
c       qmod(iq,1)=ntdif
c       if ( ldif.lt.kl .and. ntdif.eq.0 ) then
c       if ( tdif.lt.0. ) then
c          qmod(iq,1)=1.
c       endif
c
       enddo  ! 260

c       if ( ntest.gt.0 ) then
c         print *,'cucnvc: 260a: ',
c    .    tref(kpnt,l)-t0,t(kpnt,l)-t0,ltop(kpnt),p(kpnt)
c         print *,'cucnvc: 260b: ',
c    .    tmod(kpnt,l),tdif(kpnt),qmod(kpnt,1)
c       endif

c-----------------------------------------------------------------------
c end of level loop
c
 200  continue
c
c--------------initialize changes of t and q due to convection----------
c
c     print*,'cucnvc: initialize changes of t and q'
      do 280 l=1,kl
        do 280 iq=1,ifull
          tmod(iq,l)=0.
 280      qmod(iq,l)=0.
c
!      do iq=2*il+1,ifull-2*il-2  ! DARLAM
       do iq=1,ifull
        if ( ltop(iq).gt.lbot(iq) ) ltop(iq)=lbot(iq)
       enddo  !
c
c
      sumdse=0.
      summse=0.
      sumdpp=0.
c
c************* main horizontal loop for convection *********************
c************* main horizontal loop for convection *********************
c************* main horizontal loop for convection *********************
c************* main horizontal loop for convection *********************
c
c     print*,'cucnvc:2*il+1,ifull-2*il-2= ',2*il+1,ifull-2*il-2

!      do 310 iq=2*il+1,ifull-2*il-2  ! DARLAM
       do 310 iq=1,ifull
c
c************* main horizontal loop for convection *********************
c
      lbotk=nint(lbot(iq))
      ltpk =nint(ltop(iq))
c
c--------------adjust for any cloud at least two layers thick-----------
c
c     print*,'cucnvc: ltpk,lbotk-2: ',ltpk,(lbotk-2)
c     print*,'avgrfi: ',avgrfi
      if ( ltpk.gt.lbotk-2 ) then
         cldefi(iq)=avgrfi
         go to 310
      endif
c
      if ( ntest.gt.1 ) then
         print *,'cucnvc: ltpk,lbotk,lsh=',ltpk,lbotk,lsh
      endif
c
      tauk=etrel
c
c***********************************************************************
c
c           ***select convection scheme***
c
ccccc go to deep convection if cloud is more than lsh layers deep ***
ccccc if ( ltpk.lt.lbotk-lsh ) go to 400
c
c go to deep convection if cloud is more than 200.e2 pa deep
c
      pct=pdsl(iq)*aeta(ltpk)+pt
      pcb=pdsl(iq)*aeta(lbotk)+pt
      if ( ntest.gt.1 )  print *,'cucnvc: pct,pcb ',pct,pcb
      if ( (pcb-pct).gt.200.e2 ) go to 400

c***********************************************************************
c
cscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscs
cscscscscscsc  shallow convection  cscscscscscscscscscscscscscscscscscsc
cscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscs
c                                                                      c
c    sss h        l l                                                  c
c    s   h   aaa  l l ooo w     w                                      c
c     ss hhh a a  l l o o  w w w                                       c
c    sss h h aaaa l l ooo   w w                                        c
c
c--------------scheme mixes partially one layer above top---------------
c
 600  continue

      if ( .not.oshal ) go to 310 ! i.e. to bottom of horiz loop
      ltpk=ltpk-1
      if ( ntest.gt.1 )  print *,'600 ltpk= ',ltpk

c--------------scaling potential temperature & table index at top-------
      thtpk=t(iq,ltpk)*ape(iq,ltpk)
      tthk =(thtpk-thl)*rdth
      qqk  =tthk-aint(tthk)
      it   =int(tthk)+1
      if ( it.lt.1 ) then
c        if ( ntest.gt.1 ) then      !.and. iq.eq.kpnt ) then
c           print *,'600 it,qqk ',it,qqk
c        endif
         it=1
         qqk=0.
      end if
      if ( it.ge.jtb ) then
c        if ( ntest.gt.1 .and. iq.eq.kpnt ) then
c           print *,'600 it,qqk ',it,qqk
c        endif
         it=jtb-1
         qqk=0.
      end if

c--------------base and scaling factor for spec. humidity at top--------
      bqs00k=qs0(it)
      sqs00k=sqs(it)
      bqs10k=qs0(it+1)
      sqs10k=sqs(it+1)

c--------------scaling spec. humidity & table index at top--------------
      bqk=(bqs10k-bqs00k)*qqk+bqs00k
      sqk=(sqs10k-sqs00k)*qqk+sqs00k
      tqk=(q(iq,ltpk)-bqk)/sqk*rdq
      ppk  =tqk-aint(tqk)
      iqbetts   =int(tqk)+1
      if ( iqbetts.lt.1 ) then
c        if ( ntest.gt.1 .and. iq.eq.kpnt ) then
c           print *,'600 iqbetts,ppk ',iqbetts,ppk
c        endif
         iqbetts=1
         ppk=0.
      end if
      if ( iqbetts.ge.itb)  then
c        if ( ntest.gt.1 .and. iq.eq.kpnt ) then
c           print *,'600 iqbetts,ppk ',iqbetts,ppk
c        endif
         iqbetts=itb-1
         ppk=0.
      end if

c--------------cloud top saturation point pressure----------------------
      ptpk= ptbl(iqbetts,it  )
     2    +(ptbl(iqbetts+1,it)-ptbl(iqbetts,it))*ppk
     3    +(ptbl(iqbetts,it+1)-ptbl(iqbetts,it))*qqk
     4    +(ptbl(iqbetts  ,it  )-ptbl(iqbetts+1,it  )
     5     -ptbl(iqbetts  ,it+1)+ptbl(iqbetts+1,it+1))*ppk*qqk

c-----------------------------------------------------------------------
      dpmix=ptpk-psp(iq)
      if ( abs(dpmix).lt.3000.)  dpmix=-3000.
      smix=stabs*(thtpk-thbt(iq))/dpmix

c       if ( ntest.gt.1 .and. iq.eq.kpnt ) then
c          print *,'smix ',smix,dpmix,ptpk,psp(iq),thtpk,thbt(iq)
c       endif

c--------------shallow reference profiles (first guess)-----------------
      do 320 l=1,kl
        trefk(l)=t  (iq,l)
        qrefk(l)=q  (iq,l)
        pk   (l)=aeta(l)*pdsl(iq)+pt
        apek (l)=ape(iq,l)

c       if ( ntest.gt.1 .and. iq.eq.kpnt ) then
c          print *,'320 ',l,trefk(l)-t0,qrefk(l),pk(l),apek(l)
c       endif

 320  continue
c
      do 358 iter=1,2
c
c--------------reference temperature (first guess)----------------------
c
      ltp1=ltpk+1
c
      ivi=lbotk-1
 330  ivi=ivi-1
c
      trefk(ivi)=((pk(ivi)-pk(ivi+1))*smix
     2            +trefk(ivi+1)*apek(ivi+1))/apek(ivi)

c       if ( ntest.gt.1 .and. iq.eq.kpnt ) then
c          print *,'330 ',ivi,trefk(ivi)-t0
c       endif

      if ( ivi.gt.ltp1 ) go to 330
c
c--------------reference spec. humidity (first guess)-------------------
c
      ivi=lbotk-1
 340  ivi=ivi-1
c
      pskl=pk(ivi)+dsps
      thskl=trefk(ivi)*apek(ivi)
      apeskl=(pskl/100000.)**dm2859
      qrefk(ivi)=pq0/pskl*exp(a2*(thskl-t0*apeskl)
     2                          /(thskl-a4*apeskl))

c       if ( ntest.gt.1 .and. iq.eq.kpnt ) then
c          print *,'340 ',ivi,pskl,thskl,apeskl,qrefk(ivi)
c       endif

      if ( ivi.gt.ltp1 ) go to 340
c
c--------------reference profiles corrections---------------------------
c
      sumdt=0.
      sumdq=0.
      sumdp=0.
c
      l=ltpk
c
 350  l=l+1
c
      sumdt=(t(iq,l)-trefk(l))*deta(l)+sumdt
      sumdq=(q(iq,l)-qrefk(l))*deta(l)+sumdq
      sumdp=sumdp+deta(l)

c       if ( ntest.gt.1 .and. iq.eq.kpnt ) then
c          print *,'350 ',sumdt,sumdq,sumdp
c       endif

      if ( l.lt.lbotk ) go to 350
c
caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
c
      tcorr=sumdt/sumdp
      qcorr=sumdq/sumdp

c       if ( ntest.gt.1 .and. iq.eq.kpnt ) then
c          print *,'tqcorr ',tcorr,qcorr
c       endif

      l=ltpk
c
 370  continue
c
      l=l+1
      trefk(l)=trefk(l)+tcorr
      qrefk(l)=qrefk(l)+qcorr

c       if ( ntest.gt.1 .and. iq.eq.kpnt ) then
c          print *,'370  ',l,trefk(l)-t0,qrefk(l)
c       endif

      if ( l.lt.lbotk ) go to 370
c
c------------- entropy check -------------------------------------------
c
      dentpy=0.
c
      l=ltpk
c
 359  continue
c
      l=l+1
      dentpy=((trefk(l)-t(iq,l))*cp+(qrefk(l)-q(iq,l))*elwv)
     .       / (t(iq,l)+trefk(l))*deta(l)+dentpy

c       if ( ntest.gt.0 .and. iq.eq.kpnt ) then
c          print *,'359 ',l,dentpy
c       endif

      if ( l.lt.lbotk ) go to 359
c
 358  continue
c
c**********************************************************************
c
c *** if the reference profile has a smaller entropy than the
c *** model atmosphere, then skip adjustment
c
      if ( abs(dentpy).gt..5 ) go to 310
c
c**********************************************************************
c
c--------------relaxation towards reference profiles--------------------
c
      l=ltpk
 360  l=l+1
c
      tmod(iq,l)=(trefk(l)-t(iq,l))*tauk
      qmod(iq,l)=(qrefk(l)-q(iq,l))*tauk
c
      sumdse=sumdse+tmod(iq,l)*deta(l)
      summse=summse+qmod(iq,l)*deta(l)
      sumdpp=sumdpp+deta(l)

c     if ( ntest.gt.1 .and. iq.eq.kpnt ) then
c        print *,'360 ',l,tmod(kpnt,l),qmod(kpnt,l),sumdse,summse,sumdpp
c     endif

      if ( l.lt.lbotk ) go to 360
c
c-----------------------------------------------------------------------
      nshal=nshal+1
      ntops(ltpk)=ntops(ltpk)+1
      nbots(lbotk)=nbots(lbotk)+1
      ndepth=lbotk-ltpk
      ndpths(ndepth)=ndpths(ndepth)+1

c     if ( ntest.gt.1 .and. iq.eq.kpnt ) then
c        print *,'ends ',ltpk,lbotk,nshal,ntops(ltpk),nbots(lbotk)
c    .                 ,ndepth,ndpths(ndepth)
c     endif

caaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
cscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscs
cscscscscscsc  end of shallow convection   scscscscscscscscscscscscscscs
cscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscs
c-----------------------------------------------------------------------
                             go to 310
c-----------------------------------------------------------------------
cdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd
cdcdcdcdcdcdc  deep convection   dcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd
cdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd
c
 400  lb  = klh(iq)-1
      efi = cldefi(iq)
c
c***********************************************************************
c
c *** if the 'cloud efficiency' (cldefi) is equal to efmin (<1) then
c *** the deep reference profiles will be those of the moister,
c *** supposedly slower, systems.  if cldefi is 1 then the deep
c *** reference profiles will be those of the drier, faster
c *** sytems.  values of cldefi between efimn and 1 will yield
c *** deep reference profiles that lie (linearly) between these
c *** lilits.
c
c***********************************************************************
c
      dspbk=((efi-efimn)*slopb+dspbs)*sm(iq)+(1.-sm(iq))*dspbl
      dsp0k=((efi-efimn)*slop0+dsp0s)*sm(iq)+(1.-sm(iq))*dsp0l
      dsptk=((efi-efimn)*slopt+dspts)*sm(iq)+(1.-sm(iq))*dsptl
c
      if ( ntest.gt.1 ) then
         print *,'deep ',lb,efi,efimn,sm(iq)
         print *,dspbk,slopb,dspbs,dspbl
         print *,dsp0k,slop0,dsp0s,dsp0l
         print *,dsptk,slopt,dspts,dsptl
      endif
c
      do 410 l=1,kl
        dift (l)=0.
        difq (l)=0.
        tk   (l)=t  (iq,l)
        trefk(l)=t  (iq,l)
        qk   (l)=q  (iq,l)
        qrefk(l)=q  (iq,l)
        pk   (l)=aeta(l)*pdsl(iq)+pt
        psk  (l)=pk(l)
        apek (l)=ape(iq,l)
 410    therk(l)=tref(iq,l)*apek(l)
c
      if ( ntest.gt.1 ) then
         print *,'cucnvc: lbotk,lb,ltpk,pk,apek',
     .   lbotk,lb,ltpk,pk,apek
      endif
c
c--------------deep convection reference temperature profile------------
c
      ltp1=ltpk-1
      lbm1=lb-1
      pkb=pk(lb)
      pkt=pk(ltpk)
c
c--------------temperature reference profile below freezing level-------
c
      l0=lb
      pk0=pk(lb)
c     print *,'cucnvc:pk0,lb= ',pk0,lb
c     print*,'loop 420,ltpk,lbm1= ',ltpk,lbm1
      do 420 l=ltpk,lbm1
        ivi=ltpk+lbm1-l
        if ( trefk(ivi+1).le.tfr)  then
           go to 430
        endif
        stabdl=stabd
        if ( ivi.eq.lbm1 .and. sm(iq).gt.0. ) stabdl=stabd*stabfc
        trefk(ivi)=((therk(ivi)-therk(ivi+1))*stabdl
     2            +trefk(ivi+1)*apek(ivi+1))/apek(ivi)
        l0=ivi
        pk0=pk(l0)
c       print*,'in loop 420, l0,l=',l0,l
c       print 112,pkt,pk0
c
      if ( ntest.gt.0 ) then
         trc=trefk(ivi)-t0
         tkc=tk(ivi)-t0
         print *,'cucnvc:ivi,trc,tkc,pk ',ivi,trc,tkc,pk(ivi)
      endif
c
 420  continue
c
c--------------freezing level at or above the cloud top-----------------
c
      l0m1=l0-1
c     print 112,pk0,pkt
c112   format('pk0,pkt= ',2(' ',E20.14))
      go to 440
c
c--------------temperature reference profile above freezing level-------
c
 430  l0m1=l0-1
c     print*,'cucnvc: 430, pk0,pkt= ',pk0,pkt
      rdp0t=1./(pk0-pkt)
      dthem=therk(l0)-trefk(l0)*apek(l0)
c
      if ( ntest.gt.1 ) print *,'cucnvc: l,trefk,tk,pk'
      do 425 l=l0m1,ltpk,-1
        trefk(l)=(therk(l)-(pk(l)-pkt)*dthem*rdp0t)/apek(l)
c
c       if ( ntest.gt.1 ) then
c          trc=trefk(l)-t0
c          tkc=tk(l)-t0
c          print *,'cucnvc:,l,trc,tkc,pk= ',l,trc,tkc,pk(l)
c       endif
c
 425  continue
c
c--------------deep convection reference humidity profile---------------
c
c
 440  continue
    
      if ( ntest.gt.1 ) then
         print *,'cucnvc:ltpk,lb,l0,pq0 ',ltpk,lb,l0,pq0
         print *,'saturation pressure: pk0,pkt,pkb= ',
     .   pk0,pkt,pkb
      endif

c     if(pk0.eq.pkt) then
c       print*,'pk0.eq.pkt,(lb-l0)=',(lb-l0)
c     else 
c       print*,'pk0.NE.pkt,(lb-l0)=',(lb-l0)
c     endif
c     print 111,pk0,pkt,pkb
c111   format('before loop 450:pk0,pkt,pkb= ',3('  ',E20.14))
c     print*,'before loop 450,ltpk,lb= ',ltpk,lb
      do 450 l=ltpk,lb
c        print*,'loop 450: (lb-l0),l= ',(lb-l0),l
c
c--------------saturation pressure difference---------------------------
c
      if ( lb-l0.ge.3 ) then
c        print*,'lb-l0,l ',lb-l0,l
         if ( l.lt.l0 ) then
c           print*,'1 dsp,l,l0= ',l,l0
            dsp=((pk0-pk(l))*dsptk+(pk(l)-pkt)*dsp0k)/(pk0-pkt)
         else
c           print*,'2 dsp,l,l0= ',l,l0
            dsp=((pkb-pk(l))*dsp0k+(pk(l)-pk0)*dspbk)/(pkb-pk0)
         end if
      else
         dsp=dspc
      end if
c
c--------------humidity profile-----------------------------------------
      if ( ntest.gt.1 ) then
         print *,'humidity: l,100000.,dm2859= ',l,100000.,dm2859
         print *,'humidity: a2,(thsk-t0*apesk)=',
     .   (thsk(l)-t0*apesk(l))
         print *,'humidity: a4,(thsk-a4*apesk)=',
     .   (thsk(l)-a4*apesk(l))
      endif
c
      psk(l)=pk(l)+dsp
      apesk(l)=(psk(l)/100000.)**dm2859
      thsk(l)=trefk(l)*apek(l)
c
      if ( l.gt.lqm ) then
         qrefk(l)=pq0/psk(l)*exp(a2*(thsk(l)-t0*apesk(l))
     2                             /(thsk(l)-a4*apesk(l)))
      else
         qrefk(l)=q(iq,l)
      end if
c
        if ( ntest.gt.1 ) then
           print *,'cucnvc:450a ',l,qrefk(l),qk(l),dsp,psk(l)
           print *,'cucnvc:450b ',
     .     l,apesk(l),thsk(l),trefk(l),apek(l)
        endif
c
 450  continue
c
c--------------enthalpy conservation integral---------------------------
c
      do 485 iter=1,2
c
        if ( ntest.gt.1 ) then
           print *,'iter,ltpk,lb,lqm ',iter,ltpk,lb,lqm
        endif
c
c-----------------------------------------------------------------------
c
        sumde=0.
        sumdp=0.
        do 460 l=ltpk,lb
          sumde=((tk(l)-trefk(l))*cp+(qk(l)-qrefk(l))*elwv)*deta(l)
     2         +sumde
 460      sumdp=sumdp+deta(l)
        hcorr=sumde/(sumdp-deta(ltpk))
        lcor=ltpk+1
c
        if ( ntest.gt.1 ) then
           print *,'460 ',sumde,sumdp,hcorr,lcor
        endif
c
c------------ above lqm correct temperature only -----------------------
c
        if ( lcor.le.lqm ) then
c
           if ( ntest.gt.1 ) print *,'470: l,trefk'
           do 470 l=lcor,lqm
             trefk(l)=trefk(l)+hcorr*rcp
             if ( ntest.gt.0 ) then
                trc=trefk(l)-t0
                print *,'470 above lqm ',l,trc
             endif
 470       continue
c
        endif
c
c------------ below lqm correct both temperature and moisture ----------
c
        if ( ntest.gt.1 ) print *,'480: l,trefk,qrefk'
        do 480 l=lcor,lb
          tskl=trefk(l)*apek(l)/apesk(l)
          dhdt=qrefk(l)*a23m4l/(tskl-a4)**2+cp
          trefk(l)=hcorr/dhdt+trefk(l)
          thskl=trefk(l)*apek(l)
          qrefk(l)=pq0/psk(l)*exp(a2*(thskl-t0*apesk(l))
     .                             /(thskl-a4*apesk(l)))
          if ( ntest.gt.1 ) then
             trc=trefk(l)-t0
             print *,'480 below lqm ',l,trc,qrefk(l)
          endif
 480    continue
c
c-----------------------------------------------------------------------
c
c end of iteration loop
c
 485  continue
c
c-----------------------------------------------------------------------
c
c thvref for diagnostics only
c
      do 680 l=1,kl
        thvref(l)=trefk(l)*apek(l)*(qrefk(l)*d608+1.)
        thvmod(l)=tk   (l)*apek(l)*(qk   (l)*d608+1.)
c
        if ( ntest.gt.0 ) then
           print *,'l,thvref,thvmod ',l,thvref(l),thvmod(l)
        endif
c
 680  continue
c
c--------------heating, moistening, precipitation-----------------------
c
      dentpy=0.
      avrgt =0.
      preck=0.
c
      if ( ntest.gt.1 ) then
         print *,'500a: l,tauk,deta(l),pk(l)'
         print *,'trefk,tk,dift / qrefk,qk,difq'
      endif
      do 500 l=ltpk,lb
c
        dift(l)=(trefk(l)-tk(l))*tauk
        difq(l)=(qrefk(l)-qk(l))*tauk
        avrgt  =(tk(l)+tk(l)+dift(l))*deta(l)+avrgt
        dentpy =(dift(l)*cp+difq(l)*elwv)
     .         /(tk(l)+tk(l)+dift(l))*deta(l)+dentpy
        preck=deta(l)*dift(l)+preck
c
        if ( ntest.gt.1 ) then
           print *,'500b: ',l,tauk,deta(l),pk(l)
           trc=trefk(l)-t0
           tkc=tk(l)-t0
           qrc=qrefk(l)
           qkc=qk(l)
           print *,'500c: ',trc,tkc,dift(l)
           print *,'                         ',qrc,qkc,difq(l)
        endif
c
 500  continue
c
      dentpy=dentpy+dentpy
      avrgt =avrgt/(sumdp+sumdp)
c
      if ( ntest.gt.1 ) then
         print *,'dentpy ',avrgt,dentpy,preck
      endif
c
cscscscscscscs if negative prec., or entropy, shallow adjustment scscscs
c
      if ( dentpy.lt.epsntp .or. preck.le.0. ) then
         cldefi(iq)=efimn
ccccc    ltpk=lbot(iq)
         lbm1=nint(lbot(iq)-1.)
c
c *** checking from cloudtop downward, if the static stability decreases
c *** by more than dthvt degrees/deta, then set the cloudtop for swapped
c *** shallow convection at the layer in common to both finite
c *** differences, i.e., just below the lowest stable layer.
c
         do 695 l=ltp1,lbm1
c
ccccc      d2thv=thvmod(l-1)+thvmod(l+1)-thvmod(l)-thvmod(l)
           d2thv=-((thvmod(l-1)-thvmod(l))/(aeta(l-1)-aeta(l))
     .            -(thvmod(l)-thvmod(l+1))/(aeta(l)-aeta(l+1)))*deta(l)
c
ccccc      if ( d2thv.gt.dthvt ) ltpk=l
           if ( ntest.gt.0)  then
              print *,'695: ',l,d2thv,dthvt,ltpk,thvmod(l)
           endif
 695     continue
c
         if ( ltpk.lt.lbot(iq)-lsh ) ltpk=nint(lbot(iq)-lsh)
         if ( ltpk.gt.lbot(iq)-2   ) go to 310
c
         ltop(iq)=ltpk
         nneg=nneg+1
         nbotn(lbotk)=nbotn(lbotk)+1
         ntopn(ltpk)=ntopn(ltpk)+1
         ndepth=lbotk-ltpk
         ndpthn(ndepth)=ndpthn(ndepth)+1
c
         if ( ntest.gt.1 ) then
            print *,'endn ',ltop(iq),lbotk,nneg,ntopn(ltpk),nbotn(lbotk)
     .                     ,ndepth,ndpthn(ndepth)
         endif
c
         go to 600
c
      endif
c
      drheat=preck*cp/avrgt
      efi=efifc*dentpy/drheat
      efi=((cldefi(iq)+efi)*.5)*sm(iq)+1.-sm(iq)
      if ( efi.gt.1. ) efi=1.
      if ( efi.lt.efimn ) efi=efimn
      cldefi(iq)=efi
      preck=preck*efi
c
      if ( ntest.gt.1 ) then
         print *,'efi ',efi,drheat,preck,ltpk,lsh
      endif
c
c--------------update precipitation, temperature & moisture-------------
c
      prec  (iq)=pdsl(iq)*preck*cprlg+prec  (iq)
      cuprec(iq)=pdsl(iq)*preck*cprlg+cuprec(iq)
c
      do 510 l=ltpk,lb
        tmod(iq,l)=dift(l)*efi
 510    qmod(iq,l)=difq(l)*efi
c
      ndeep=ndeep+1
      ntopd(ltpk)=ntopd(ltpk)+1
      nbotd(lb)=nbotd(lb)+1
      ndepth=lb-ltpk
      ndpthd(ndepth)=ndpthd(ndepth)+1
c
      if ( ntest.gt.1 ) then
         print *,'endd ',ltpk,lb,ndeep,ntopd(ltpk),nbotd(lb)
     .                 ,ndepth,ndpthd(ndepth)
      endif
c
cdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd
cdcdcdcdcdcdc  end of deep convection  dcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd
cdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd
c
c************* end of horizontal convection loop ***********************
c
 310  continue
c
      if ( sumdpp.gt.1.e-20 ) then
         sumdse=sumdse/sumdpp
         summse=summse/sumdpp
         if ( ntest.gt.0 )  then
            print *,'sumdse,summse,sumdpp= ',sumdse,summse,sumdpp
         endif
      endif

      do 750 l=1,kl
      do 750 iq=1,ifull
        t(iq,l)=t(iq,l)+tmod(iq,l)
 750    q(iq,l)=q(iq,l)+qmod(iq,l)
c
c--------------save cloud top and bottom for radiation------------------
c
c      do 800 iq=1,ifull
c        htop(iq)=real(ltop(iq))
c        hbot(iq)=real(lbot(iq))
c 800  continue
c      if ( ntest.gt.0 )  then
c         print *,'htop/bot= ',htop(kpnt),hbot(kpnt)
c      endif
c
c-----------------------------------------------------------------------

c     if ( ntest.gt.0 ) then
c       print *,'nshal,nneg,ndeep=',nshal,nneg,ndeep
c       print *,'        shallow        middle          deep      '
c       print *,' lev,bot ,top ,dpth,bot ,top ,dpth,bot ,top ,dpth'
c       do l=1,kl
c        write(6,886) l,nbots(l),ntops(l),ndpths(l)
c    .                 ,nbotn(l),ntopn(l),ndpthn(l)
c    .                 ,nbotd(l),ntopd(l),ndpthd(l)
c886     format(1x,10(i3,2x))
c       enddo
c     endif

c-----------------------------------------------------------------------
c
      return
      end
