module davies_m

private
public davies,davset
public dav_init,dav_end
public davt,davu
public psls,qgg,tt,uu,vv
public xtgdav

real, dimension(:), allocatable, save :: davt,davu
real, dimension(:), allocatable, save :: psls
real, dimension(:,:), allocatable, save :: qgg,tt,uu,vv
real, dimension(:,:,:), allocatable, save :: xtgdav

contains
    
subroutine davies    ! for globpea - only large-scale available

use aerosolldr, only : xtg,naero
use arrays_m         ! t,u,v,ps
use cc_mpi, only : mydiag
use sigs_m           ! sig

implicit none

integer iq,k,n
integer kupper
integer, parameter :: ntest=0
real speed,speeduu,rat
!     from Nov 05, separate davu & davt, just abs(nud_hrs) used)
include 'newmpar.h' ! il,jl,kl,ij
include 'parm.h' ! kbotdav,nud_u,nud_v,nud_t,nud_p,nud_q,nud_hrs,nudu_hrs

!     and new nud_p, nud_t, nud_q, nud_uv just off/on switches (0/1)
!       large-scale style with davies-style weights (already scaled for nud_hrs)
!       N.B. nbd.lt.0 sets up special weights in indata.f
if(nmaxpr.eq.1.and.ktau.lt.5.and.mydiag)then
  write(6,*) 'davies in  uu,vv,qgg(kl) ',uu(idjd,nlv),vv(idjd,nlv),qgg(idjd,kl)
  write(6,*) 'davies in  u,v,qg(kl) ',u(idjd,nlv),v(idjd,nlv),qg(idjd,kl)
endif
if(nud_p.ne.0)then
  do iq=1,ifull
    psl(iq)=psl(iq)+(psls(iq)-psl(iq))*davt(iq)*dt/3600.
    ps(iq)=1.e5*exp(psl(iq))  ! Jan 07
  enddo  ! iq loop
endif  ! (nud_p.ne.0)
if(nud_t.ne.0)then
  do k=kbotdav,ktopdav
    do iq=1,ifull
      t(iq,k)=t(iq,k)+(tt(iq,k)-t(iq,k))*davt(iq)*dt/3600.
    enddo  ! iq loop
  enddo   ! k loop
endif     ! (nud_t.ne.0)
if(nud_q.ne.0)then
  if(nud_q.lt.0)then
    kupper=kl/2 ! option to only nudge bottom half of atmos for qg
  else
    kupper=ktopdav
  endif  ! (nud_q.lt.0)
  do k=kbotdav,kupper
    do iq=1,ifull
      qg(iq,k)=qg(iq,k)+(qgg(iq,k)-qg(iq,k))*davt(iq)*dt/3600.
    enddo  ! iq loop
  enddo   ! k loop
endif     ! (nud_q.ne.0)
if(nud_uv.eq.1)then
  do k=kbotu,ktopdav
    do iq=1,ifull
      u(iq,k)=u(iq,k)+(uu(iq,k)-u(iq,k))*davu(iq)*dt/3600.
      v(iq,k)=v(iq,k)+(vv(iq,k)-v(iq,k))*davu(iq)*dt/3600.
    enddo  ! iq loop
  enddo   ! k loop
  if(kbotdav<kbotu)then  ! Feb 06
    do k=kbotdav,kbotu-1
      do iq=1,ifull
        u(iq,k)=u(iq,k)+(uu(iq,k)-u(iq,k))*davt(iq)*dt/3600.
        v(iq,k)=v(iq,k)+(vv(iq,k)-v(iq,k))*davt(iq)*dt/3600.
      enddo  ! iq loop
    enddo   ! k loop
  endif
endif     ! (nud_uv.eq.1)
if(nud_uv.eq.2)then   ! speed option (not useful)
  do k=kbotu,ktopdav
    do iq=1,ifull
      speed=sqrt(u(iq,k)**2+v(iq,k)**2)
      if(speed.gt.1.)then
        speeduu=sqrt(uu(iq,k)**2+vv(iq,k)**2)
        rat=speeduu/speed
!       u(iq,k)=u(iq,k)+(rat*u(iq,k)-u(iq,k))*davt(iq)*dt/3600.
        u(iq,k)=u(iq,k)+u(iq,k)*(rat-1.)*davu(iq)*dt/3600.
        v(iq,k)=v(iq,k)+v(iq,k)*(rat-1.)*davu(iq)*dt/3600.
      endif ! (speed.gt.1.)
    enddo  ! iq loop
  enddo   ! k loop
endif     ! (nud_uv.eq.2)
if(abs(iaero)>=2.and.nud_aero/=0)then
  do n=1,naero
    do k=kbotdav,ktopdav
      do iq=1,ifull
        xtg(iq,k,n)=xtg(iq,k,n)+(xtgdav(iq,k,n)-xtg(iq,k,n))*davt(iq)*dt/3600.
      enddo  ! iq loop
    enddo   ! k loop
  enddo    ! n loop
endif     ! (abs(iaero)>=2.and.nud_aero/=0)
if(nmaxpr.eq.1.and.ktau.lt.5.and.mydiag)then
  write(6,*) 'davies out u,v,qg(kl) ',u(idjd,nlv),v(idjd,nlv),qg(idjd,kl)
endif

return
end subroutine davies

!=======================================================================
subroutine davset

use aerosolldr
use arrays_m        ! t,u,v,ps

implicit none

include 'newmpar.h' ! il,jl,kl,ij
include 'parm.h'

psls(1:ifull) = psl(1:ifull)
tt(1:ifull,:) = t(1:ifull,:)
qgg(1:ifull,:) = qg(1:ifull,:)
uu(1:ifull,:) = u(1:ifull,:)
vv(1:ifull,:) = v(1:ifull,:)
if(abs(iaero)>=2.and.nud_aero/=0)then
  xtgdav(:,:,:) = xtg(1:ifull,:,:)
end if

return
end subroutine davset

subroutine dav_init(ifull,iextra,kl,naero)

implicit none

integer, intent(in) :: ifull,iextra,kl,naero

allocate(davt(ifull),davu(ifull))
allocate(psls(ifull),qgg(ifull,kl),tt(ifull,kl),uu(ifull,kl),vv(ifull,kl))
if (naero>0) then
  allocate(xtgdav(ifull,kl,naero))
end if

return
end subroutine dav_init

subroutine dav_end

implicit none

deallocate(davt,davu)
deallocate(psls,qgg,tt,uu,vv)
if (allocated(xtgdav)) then
  deallocate(xtgdav)
end if

return
end subroutine dav_end

end module davies_m