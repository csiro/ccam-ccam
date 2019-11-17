! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2017 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
      
      subroutine betts ( tin, qg, tn, land, ps)  ! jlm
c     can get rid of betts1.h sometime
c     ******************************************************************
c     *                                                                *
c     * setup betts scheme for csiro model                             *
c     *                                                                *
c     ******************************************************************
      use betts1_m !  includes work3a
      use morepbl_m
      use newmpar_m
      use parm_m
      use prec_m
      use sigs_m

      dimension tin(ifull+iextra,kl), qg(ifull+iextra,kl)
      dimension tn(ifull+iextra,kl)
      dimension land(ifull), ps(ifull)
      logical ofirst, land, oshal
      data oshal/.false./  ! jlm

      data ofirst/.true./

c-----------------------------------------------------------------------

      save
      ipnt=id
      jpnt=jd
      kpnt=idjd  ! jlm
!     kpnt=ipnt+(jpnt-1)*il
      
      if ( ofirst ) then

         call betts1_init(ifull,iextra,kl)

c setup eta fields

         do 30 l=1,kl
           aeta(l)=sig(kl+1-l)
 30      continue

         eta(kl+1)=1.
         eta(1)=0.
         do 40 l=2,kl
           eta(l)=sigmh(kl+2-l)
 40      continue

         do 50 l=1,kl
           deta(l)=-dsig(kl+1-l)
 50      continue

         print *,'eta          aeta     deta'
         do 60 l=1,kl+1
           print *,eta(l)
           if ( l.lt.kl+1 ) print *,'             ',aeta(l),deta(l)
 60      continue

         pt=0.

         do 70 iq=1,ifull
           res(iq)=1.
           hbm2(iq)=1.
c sm=1 over sea
           sm(iq)=1.
           if ( land(iq) ) sm(iq)=0.
           klh(iq)=kl
 70      continue

         do 80 j=1,2
           do 80 i=1,il
             iq=i+(j-1)*il
             hbm2(iq)=0.
 80      continue
         do 90 j=1,jl
           do 90 i=1,2
             iq=i+(j-1)*il
             hbm2(iq)=0.
 90      continue

         ofirst=.false.

      endif

c**********************************************************************

c initialize variables to pass to rain/cucnvc

c**********************************************************************

!     do l=1,kl
!       lc=kl+1-l
!       do j=1,jl
!         do i=1,il
!           iq=i+(j-1)*il
!           t(iq,l)=tin(i,j,lc)
!           q(iq,l)=qg(i,j,lc)
!         end do ! do i=1,il
!       end do ! do j=1,jl
!       do iq=1,ifull
!         htm(iq,l)=hbm2(iq)
!       end do ! do i=1,ifull
!     end do ! do l=1,kl

      do l=1,kl
       do iq=1,ifull
        kk=kl+1-l
        t(iq,l)=tin(iq,kk)
        q(iq,l)=qg(iq,kk)
        htm(iq,l)=hbm2(iq)
       end do ! do iq=1,ifull
      end do ! do l=1,kl

      do iq=1,ifull
        prec(iq)=0.
        cuprec(iq)=0.
        pd(iq)=ps(iq)-pt
      enddo ! iq=1,ifull

c     print *,'kpnt,pd ',kpnt,pd(kpnt)

c**********************************************************************

      call bettrain (  pt, kpnt )  ! was rain

!      if (jacks_heating.eq.1)then
!         do l=1,kl
!           lc=kl+1-l
!           do j=1,jl
!             do i=1,il
!               iq=i+(j-1)*il
!               dtl(i,j,lc)=dtl(i,j,lc)+t(iq,l)-tin(i,j,lc)
!               tnew(i,j,lc)=t(iq,l)
!               dql(i,j,lc)=dql(i,j,lc)+q(iq,l)-qg(i,j,lc)
!               qg(i,j,lc)=q(iq,l)
!             end do ! do i=1,il
!           end do ! do j=1,jl
!         end do ! do l=1,kl
!      endif

c**********************************************************************

      call bett_cuc (  dt, pt, kpnt, oshal )  ! was cucnvc

!      if (jacks_heating.eq.1) then
!         do l=1,kl
!           lc=kl+1-l
!           do j=1,jl
!             do i=1,il
!               iq=i+(j-1)*il
!               dtc(i,j,lc)=dtc(i,j,lc)+t(iq,l)-tnew(i,j,lc)
!               dqc(i,j,lc)=dqc(i,j,lc)+q(iq,l)-qg(i,j,lc)
!             end do ! do i=1,il
!           end do ! do j=1,jl
!         end do ! do l=1,kl
!      endif

c**********************************************************************

c reset variables to be passed back to the main routine

        do l=1,kl
         do iq=1,ifull
          kk=kl+1-l
          tn(iq,kk)=tn(iq,kk)+(t(iq,l)-tin(iq,kk))/dt  ! jlm
          qg(iq,kk)=q(iq,l)
         end do ! do iq=1,ifull
        end do ! do l=1,kl

c accumulate precip. and store this timestep's in condx
      do iq=1,ifull
        condx(iq)=prec(iq)*1.e3
        condc(iq)=cuprec(iq)*1.e3
        precip(iq)=precip(iq)+condx(iq)
        precc (iq)=precc (iq)+condc(iq)
      enddo
      
      return
      end
