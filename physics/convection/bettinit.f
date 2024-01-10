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
      
      subroutine bettinit (  dt, pt )  ! was initcuc
      use betts1_m
      use newmpar_m
      parameter (ntest=0)  ! replaces debug; set to 1 for degugging
c     this is part of the Betts-Miller parameterization
c----------------------------------------------------------------------

      dimension aqt(jtb), apt(jtb)
     .        , qsold(jtb), pold(jtb), theold(jtb), y2p(jtb)
     .        , qsnew(jtb), told(jtb), thenew(jtb), y2t(jtb)
     .        , pnew (jtb), tnew(jtb), app   (jtb), aqp(jtb)

c-----------------------------------------------------------------------

      data c0/0./, c1/1./, c1e5/1.e5/
      data thh/350./, ph/1.05e5/
      data cp/1004./, r/287.04/
      data a2/17.2693882/, t0/273.16/, a4/35.86/
      data pq0/379.90516/, eps/1.e-10/, elwv/2.5e6/

c-----------------------------------------------------------------------

      save

c coarse look-up table for saturation point

      kthm = jtb
      kpm  = itb
      kthh = jtb/2
      kph  = itb/2
      kthm1= kthm - 1
      kpm1 = kpm - 1
      rdthe= kthm1
      dthe = c1 / rdthe
      rdq  = kpm1
      dqs  = c1 / rdq

      dtq2 = dt

      pl   = max(100.e2,pt)

      dth  = (thh-thl)*dthe
      dp   = (ph -pl )*dqs

      rdth = c1/dth
      rdp  = c1/dp

      do 500 kth = 1, kthm

        th = thl + dth * real (kth-1)

        do 510 kp = 1, kpm

          pres = pl + dp * real (kp-1)
          if ( pres.lt.c1 ) pres = c1
          ape = ( c1e5/pres ) ** ( r/cp )
          denom = (th-a4*ape)

          if ( denom .le. c0 ) then
             print *,'pres,ape,th,denom<=0',pres,ape,th,denom
             qsold(kp) = c0
          else
             qsold(kp) = pq0*exp( a2*(th-t0*ape)/denom )/pres
          endif
          pold (kp) = pres

 510    continue

        qs0k = qsold(1)
        sqsk = qsold(kpm)-qsold(1)
        qsold(1)   = c0
        qsold(kpm) = c1

        do 520 kp = 2 , kpm1

          qsold(kp) = (qsold(kp)-qs0k)/sqsk

c fix for cyber half precision

          if((qsold(kp)-qsold(kp-1)).lt.eps)qsold(kp)=qsold(kp-1)+eps

 520    continue

        qs0(kth) = qs0k
        sqs(kth) = sqsk

        qsnew(1)   = c0
        qsnew(kpm) = c1

        do 530 kp = 2 , kpm1

          qsnew(kp) = qsnew(kp-1)+dqs

 530    continue

        y2p(1)   = c0
        y2p(kpm) = c0

        call bettspli(kpm,qsold,pold,y2p,kpm,qsnew,pnew,app,aqp)  ! was spline

        if ( ntest.eq.1) then
           print *,'th=',th
           print *,'pnew(kph),pold(kph)',pnew(kph),pold(kph)
           print *,'qs0k,sqsk',qs0k,sqsk
           print *,'qsnew(kph),qsold(kph)',qsnew(kph),qsold(kph)
           if ( kth.eq.kthh ) print *,'kp,qsnew(kp),pnew(kp)'
        endif

        do 540 kp = 1 , kpm

          ptbl(kp,kth) = pnew(kp)

          if ( ntest.eq.1 .and. kth.eq.kthh ) then
             print *,kp,qsnew(kp),pnew(kp)
          endif

 540    continue

c end of kth loop

 500  continue

c coarse look-up table for t(p) from constant the

      do 550 kp = 1 , kpm

        p = pl + dp * real (kp-1)
        if ( p.lt.c1 ) p = c1

        do 560 kth = 1 , kthm

          th  = thl + dth * real (kth-1)
          ape = ( c1e5 / p ) ** ( r/cp )
          denom = (th-a4*ape)
          if ( denom .le. c0 ) then
             print *,'pres,ape,th,denom<=0',pres,ape,th,denom
             qqs = c0
          else
             qqs = pq0*exp( a2*(th-t0*ape)/denom )/p
          endif
          told  (kth) = th/ape
          theold(kth) = th*exp(elwv*qqs/(cp*told(kth)))

 560    continue

        the0k = theold(1)
        sthek = theold(kthm)-theold(1)
        theold(1)    = c0
        theold(kthm) = c1

        do 570 kth = 2 , kthm1

          theold(kth) = (theold(kth)-the0k)/sthek

c fix for cyber half precision

         if((theold(kp)-theold(kp-1)).lt.eps)theold(kp)=theold(kp-1)+eps

 570    continue

        the0(kp) = the0k
        sthe(kp) = sthek
        thenew(1) = c0
        thenew(kthm) = c1

        do 580 kth = 2 , kthm1

          thenew(kth) = thenew(kth-1)+dthe

 580    continue

        y2t(1) = c0
        y2t(kthm) = c0

        call bettspli(kthm,theold,told,y2t,kthm,thenew,tnew,apt,aqt) ! was spline

        if ( ntest.eq.1) then
           print *,'p=',p
           print *,'tnew(kth),told(kthh)=',tnew(kth),told(kthh)
           print *,'the0k,sthek',the0k,sthek
           print *,'thenew(kthh),theold(kthh)',thenew(kthh),theold(kthh)
           if ( kp.eq.kph ) print *,'kth,thenew(kth),tnew(kth)'
        endif

        do 590 kth = 1 , kthm

          ttbl(kth,kp) = tnew(kth)

          if ( ntest.eq.1 .and. kp.eq.kph ) then
             print *,kth,thenew(kth),tnew(kth)
          endif

 590    continue

c end of kp loop

 550  continue

c***********************************************************************

      return
      end
