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
      
      subroutine bettrain (  pt, kpnt )  ! was rain
      use betts1_m
      use newmpar_m
      parameter (ntest=0)  ! replaces debug; set to 1 for degugging
!     just returns prec, t, q (all in betts1.h)     jlm
c     this is part of the Betts-Miller parameterization
c     ******************************************************************
c     *                                                                *
c     *  large scale precipitation                                     *
c     *                                                                *
c     *  references:                                                   *
c     *                                                                *
c     *  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx  *
c     *    xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx  *
c     *    xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx  *
c     *                                                                *
c     *  Note: tresh is the threshold rel.hum. for precip.             *
c     *                                                                *
c     *                                                                *
c     ******************************************************************
c
                             parameter
     & (a2=17.2693882e0,a3=273.16e0,a4=35.86e0
     &, pq0=379.90516e0
     &, cp=1004.6e0,elwv=2.50e6,row=1.e3,g=9.8e0
     &, epsq=2.e-12,dldt=2274.e0,tresh=1.00e0)
                             parameter
     & (arcp=a2*(a3-a4)/cp,rcp=1./cp
     &, pq0c=pq0*tresh,rrog=1./(row*g))
c
      parameter (ltop=1)
!     all work2 variables are just local
      real pdsl(ifull),tl(ifull),ql(ifull),precl(ifull), 
     .  aprec(ifull),elv(ifull),qc(ifull)

c--------------preparatory calculations---------------------------------
c
      do 100 iq=1,ifull
       precl(iq)=0.
       aprec(iq)=0.
 100  pdsl(iq)=res(iq)*pd(iq)
c
c--------------padding specific humidity if too small-------------------
c
      do 110 l=1,kl
        do 120 iq=1,ifull
          if ( q(iq,l).lt.epsq ) q(iq,l)=epsq*htm(iq,l)
 120    continue
 110  continue
c
c--------------collecting precipitation from top to bottom--------------
c
      do 200 l=ltop,kl
c
c--------------set up temporary arrays ---------------------------------
c
!       do iq=2*il+1,ifull-2*il-2  ! DARLAM
        do iq=1,ifull
          tl(iq)=t(iq,l)
          ql(iq)=q(iq,l)

c-        ------l, saturation spec. humidity & cond./evap.---------------
          elv(iq)=elwv-dldt*(tl(iq)-a3)
          qc(iq)=htm(iq,l)*pq0c/(pdsl(iq)*aeta(l)+pt)
     2         *exp(htm(iq,l)*a2*(tl(iq)-a3)/(tl(iq)-a4))
          precl(iq)=(ql(iq)-qc(iq))*hbm2(iq)*deta(l)
     2            /(elv(iq)*qc(iq)*arcp/((tl(iq)-a4)*(tl(iq)-a4))+1.)
c         if ( ntest.eq.1 .and. iq.eq.kpnt ) then
c            print*,l,elv(iq),ql(iq),qc(iq),precl(iq)
c         endif

c--------------is there enough water to evaporate ?---------------------
c
          if ( aprec(iq)+precl(iq).lt.0. ) precl(iq)=-aprec(iq)
c
c--------------collecting precipitation, modifying t & q, evaporation---
          aprec(iq)=precl(iq)+aprec(iq)
          precl(iq)=precl(iq)/deta(l)
          t(iq,l)=precl(iq)*elv(iq)*rcp+t(iq,l)
          q(iq,l)=q(iq,l)-precl(iq)
       enddo
c-----------------------------------------------------------------------
c
c end of level loop
c
 200  continue
c
c-----------------------------------------------------------------------
c
!     do iq=2*il+1,ifull-2*il-2  ! DARLAM
      do iq=1,ifull
        prec(iq)=hbm2(iq)*aprec(iq)*pdsl(iq)*rrog+prec(iq)
      enddo
      return
      end
