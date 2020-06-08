! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2019 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
module daviesnudge

private
public davies, davset
public dav_init, dav_end
public davt, davu
public psls, qgg, tt, uu, vv
public xtgdav
public vertwgt

real, dimension(:), allocatable, save :: davt, davu
real, dimension(:), allocatable, save :: psls
real, dimension(:), allocatable, save :: vertwgt
real, dimension(:,:), allocatable, save :: qgg, tt, uu, vv
real, dimension(:,:,:), allocatable, save :: xtgdav

contains
    
subroutine davies    ! for globpea - only large-scale available

use aerosolldr, only : xtg, naero
use arrays_m         ! t,u,v,ps
use cc_mpi, only : mydiag
use newmpar_m
use parm_m
use sigs_m           ! sig

implicit none

integer iq, k, n
integer kupper
real speed, speeduu, rat

k = 0

!     from Nov 05, separate davu & davt, just abs(nud_hrs) used)

!     and new nud_p, nud_t, nud_q, nud_uv just off/on switches (0/1)
!       large-scale style with davies-style weights (already scaled for nud_hrs)
!       N.B. nbd.lt.0 sets up special weights in indata.f
if( nmaxpr==1 .and. ktau<5 .and. mydiag ) then
  write(6,*) 'davies in  uu,vv,qgg(kl) ',uu(idjd,nlv),vv(idjd,nlv),qgg(idjd,kl)
  write(6,*) 'davies in  u,v,qg(kl) ',u(idjd,nlv),v(idjd,nlv),qg(idjd,kl)
endif
if ( nud_p/=0 ) then
  do iq = 1,ifull
    psls(iq) = (psls(iq)-psl(iq))*davt(iq)*dt/3600.
    psl(iq) = psl(iq) + psls(iq)
    ps(iq) = 1.e5*exp(psl(iq))  ! Jan 07
  enddo  ! iq loop
endif  ! (nud_p.ne.0)
if ( nud_t/=0 ) then
  do k = kbotdav,ktopdav
    do iq = 1,ifull
      tt(iq,k) = vertwgt(k)*(tt(iq,k)-t(iq,k))*davt(iq)*dt/3600.
      t(iq,k) = t(iq,k) + tt(iq,k)
    enddo  ! iq loop
  enddo   ! k loop
endif     ! (nud_t.ne.0)
if ( nud_q/=0 ) then
  if ( nud_q<0 ) then
    kupper = kl/2 ! option to only nudge bottom half of atmos for qg
  else
    kupper = ktopdav
  endif  ! (nud_q.lt.0)
  do k = kbotdav,kupper
    do iq = 1,ifull
      qgg(iq,k) = vertwgt(k)*(qgg(iq,k)-qg(iq,k))*davt(iq)*dt/3600.
      qg(iq,k) = qg(iq,k) + qgg(iq,k)
      qg(iq,k) = max(qg(iq,k),0.)
    enddo  ! iq loop
  enddo    ! k loop
endif      ! (nud_q.ne.0)
if ( nud_uv==1 ) then
  do k = kbotu,ktopdav
    do iq = 1,ifull
      uu(iq,k) = vertwgt(k)*(uu(iq,k)-u(iq,k))*davu(iq)*dt/3600.
      vv(iq,k) = vertwgt(k)*(vv(iq,k)-v(iq,k))*davu(iq)*dt/3600.
      u(iq,k) = u(iq,k) + uu(iq,k)
      v(iq,k) = v(iq,k) + vv(iq,k)
    enddo  ! iq loop
  enddo   ! k loop
  if ( kbotdav<kbotu ) then  ! Feb 06
    do k = kbotdav,kbotu-1
      do iq = 1,ifull
        uu(iq,k) = vertwgt(k)*(uu(iq,k)-u(iq,k))*davt(iq)*dt/3600.
        vv(iq,k) = vertwgt(k)*(vv(iq,k)-v(iq,k))*davt(iq)*dt/3600.
        u(iq,k) = u(iq,k) + uu(iq,k)
        v(iq,k) = v(iq,k) + vv(iq,k)
      enddo  ! iq loop
    enddo   ! k loop
  endif
endif     ! (nud_uv.eq.1)
if ( nud_uv==2 ) then   ! speed option (not useful)
  do k = kbotu,ktopdav
    do iq = 1,ifull
      speed = sqrt(u(iq,k)**2+v(iq,k)**2)
      if ( speed>1. ) then
        speeduu = sqrt(uu(iq,k)**2+vv(iq,k)**2)
        rat = speeduu/speed
        uu(iq,k) = vertwgt(k)*u(iq,k)*(rat-1.)*davu(iq)*dt/3600.
        vv(iq,k) = vertwgt(k)*v(iq,k)*(rat-1.)*davu(iq)*dt/3600.
        u(iq,k) = u(iq,k) + uu(iq,k)
        v(iq,k) = v(iq,k) + vv(iq,k)
      endif ! (speed.gt.1.)
    enddo  ! iq loop
  enddo   ! k loop
endif     ! (nud_uv.eq.2)
if ( abs(iaero)>=2 .and. nud_aero/=0 )then
  do n = 1,naero
    do k = kbotdav,ktopdav
      do iq = 1,ifull
        xtgdav(iq,k,n) = vertwgt(k)*(xtgdav(iq,k,n)-xtg(iq,k,n))*davt(iq)*dt/3600.
        xtg(iq,k,n) = xtg(iq,k,n) + xtgdav(iq,k,n)
      enddo  ! iq loop
    enddo    ! k loop
  enddo      ! n loop
endif        ! (abs(iaero)>=2.and.nud_aero/=0)
if ( nmaxpr==1 .and. ktau<5 .and. mydiag ) then
  write(6,*) 'davies out u,v,qg(kl) ',u(idjd,nlv),v(idjd,nlv),qg(idjd,kl)
endif

return
end subroutine davies

!=======================================================================
subroutine davset

use aerosolldr
use arrays_m        ! t,u,v,ps
use newmpar_m
use parm_m

implicit none

psls(1:ifull) = psl(1:ifull)
tt(1:ifull,:) = t(1:ifull,:)
qgg(1:ifull,:) = qg(1:ifull,:)
uu(1:ifull,:) = u(1:ifull,:)
vv(1:ifull,:) = v(1:ifull,:)
if ( abs(iaero)>=2 .and. nud_aero/=0 ) then
  xtgdav(:,:,:) = xtg(1:ifull,:,:)
end if

return
end subroutine davset

subroutine dav_init(ifull,kl,naero,nbd,mbd)

implicit none

integer, intent(in) :: ifull, kl, naero, nbd, mbd

if ( mbd/=0 .or. nbd/=0 ) then
  allocate( vertwgt(kl) )
  vertwgt(:) = 1.
end if  

if ( nbd/=0 ) then
  allocate( davt(ifull), davu(ifull) )
  allocate( psls(ifull), qgg(ifull,kl), tt(ifull,kl), uu(ifull,kl), vv(ifull,kl) )
  psls(:) = 0.
  qgg(:,:) = 0.
  tt(:,:) = 0.
  uu(:,:) = 0.
  vv(:,:) = 0.
  if ( naero>0 ) then
    allocate( xtgdav(ifull,kl,naero) )
    xtgdav(:,:,:) = 0.
  end if
end if

return
end subroutine dav_init

subroutine dav_end

implicit none

deallocate( vertwgt )

if ( allocated(davt) ) then
  deallocate( davt, davu )
  deallocate( psls, qgg, tt, uu, vv )
end if
if ( allocated(xtgdav) ) then
  deallocate( xtgdav )
end if

return
end subroutine dav_end

end module daviesnudge
