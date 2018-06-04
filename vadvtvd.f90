! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2018 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
module vadv
      
private
public vadvtvd
      
contains

subroutine vadvtvd(tarr,uarr,varr,nvadh_pass,nits,tile,imax)   ! globpea  version
!                              vadvbott & vadvyu at bottom
!     can show adding tbar has no effect
use aerosolldr
use arrays_m
use cc_mpi
use cc_omp, only : ntiles
use cfrac_m
use cloudmod
use diag_m
use liqwpar_m  ! ifullw
use map_m
use newmpar_m
use nharrs_m
use parm_m
use parmdyn_m
use sigs_m
use tkeeps, only : tke,eps
use tracers_m
use xarrs_m
implicit none
!     split vertical advection routine; tvd scheme; used with nonlin or upglobal
!     In flux limiter, assuming zero gradient for all top and bottom
!     variables; except extrap at bottom for qg and trace gases  Thu  06-19-1997
include 'kuocom.h'     ! also with kbsav,ktsav
integer, intent(in) :: tile, imax
real, dimension(:,:), intent(inout) :: tarr,uarr,varr
real, dimension(imax) :: tfact
integer ntr,k
integer, dimension(imax) :: nvadh_pass, nits
integer, save :: num = 0

!$omp master
call START_LOG(vadv_begin)
!$omp end master

tfact = 1./real(nvadh_pass)

if ( mydiag .and. ntiles==1 ) then
  if ( num==0 ) then
    num = 1
    write(6,*) 'In vadvtvd nvadh_pass ',nvadh_pass(idjd)
    write(6,*) 'tfact ',tfact(idjd)
  end if
end if

!     t
call vadv_work(tarr,tfact,nits,tile,imax)
if( (diag.or.nmaxpr==1) .and. mydiag  .and. ntiles==1 )then
!       These diagnostics don't work with single input/output argument
  write (6,"('tout',9f8.2/4x,9f8.2)") (tarr(idjd,k),k=1,kl)
  write (6,"('t#  ',9f8.2)") diagvals(tarr(:,nlv)) 
endif

!     u
call vadv_work(uarr,tfact,nits,tile,imax)
if( diag .and. mydiag .and. ntiles==1 )then
  write (6,"('uout',9f8.2/4x,9f8.2)") (uarr(idjd,k),k=1,kl)
  write (6,"('u#  ',9f8.2)") diagvals(uarr(:,nlv)) 
endif

!     v
call vadv_work(varr,tfact,nits,tile,imax)
if( diag .and. mydiag .and. ntiles==1 )then
  write (6,"('vout',9f8.2/4x,9f8.2)") (varr(idjd,k),k=1,kl)
  write (6,"('v#  ',9f8.2)") diagvals(varr(:,nlv)) 
endif

!     h_nh
if ( nh/=0 ) then
  call vadv_work(h_nh,tfact,nits,tile,imax)
end if

!     pslx
call vadv_work(pslx,tfact,nits,tile,imax)

if ( mspec==1 ) then   ! advect qg and gases after preliminary step

  !      qg
  call vadv_work(qg,tfact,nits,tile,imax)
  if ( diag .and. mydiag .and. ntiles==1 ) then
    write (6,"('qout',9f8.2/4x,9f8.2)") (1000.*qg(idjd,k),k=1,kl)
    write (6,"('qg# ',9f8.2)") diagvals(qg(:,nlv)) 
  end if

  if ( ldr/=0 ) then
    call vadv_work(qlg,tfact,nits,tile,imax)
    call vadv_work(qfg,tfact,nits,tile,imax)
    if ( ncloud>=4 ) then
      call vadv_work(stratcloud,tfact,nits,tile,imax)
    end if
    if ( diag .and. mydiag .and. ntiles==1 ) then
      write (6,"('lout',9f8.2/4x,9f8.2)") (1000.*qlg(idjd,k),k=1,kl)
      write (6,"('qlg#',9f8.2)") diagvals(qlg(:,nlv)) 
      write (6,"('fout',9f8.2/4x,9f8.2)") (1000.*qfg(idjd,k),k=1,kl)
      write (6,"('qfg#',9f8.2)") diagvals(qfg(:,nlv)) 
    end if
  end if      ! if(ldr.ne.0)

  if ( nvmix==6 ) then
    call vadv_work(eps,tfact,nits,tile,imax)
    call vadv_work(tke,tfact,nits,tile,imax)
  end if      ! if(nvmix.eq.6)

  if ( abs(iaero)>=2 ) then
    do ntr = 1,naero
      call vadv_work(xtg(:,:,ntr),tfact,nits,tile,imax)
    end do
  end if   ! abs(iaero)>=2

  if ( ngas>0 .or. nextout>=4 ) then
    do ntr = 1,ntrac
      call vadv_work(tr(:,:,ntr),tfact,nits,tile,imax)
    end do      ! ntr loop
  end if        ! (nextout>=4)
 
end if          ! if(mspec==1)

!$omp master
call END_LOG(vadv_end)
!$omp end master
 
return
end subroutine vadvtvd

subroutine vadv_work(tarr,tfact,nits,tile,imax)

use newmpar_m
use vvel_m

implicit none

integer, intent(in) :: tile, imax
integer, dimension(imax), intent(in) :: nits
integer :: is, ie
real, dimension(imax), intent(in) :: tfact
real, dimension(:,:), intent(inout) :: tarr
real, dimension(imax,kl) :: ltarr
real, dimension(imax,kl+1) :: lsdot

is = (tile-1)*imax + 1
ie = tile*imax

ltarr = tarr(is:ie,:)
lsdot = sdot(is:ie,:)

call vadv_work2(ltarr,tfact,nits,lsdot,imax)

tarr(is:ie,:) = ltarr

end subroutine vadv_work
      
! Subroutine to perform generic TVD advection
subroutine vadv_work2(tarr,tfact,nits,sdot,imax)

use cc_mpi
use newmpar_m
use sigs_m
      
implicit none
      
integer, intent(in) :: imax
integer, dimension(imax), intent(in) :: nits
integer i, k, iq, kp, kx
real, dimension(imax), intent(in) :: tfact
real, dimension(:,:), intent(inout) :: tarr
real, dimension(:,:), intent(in) :: sdot
real, dimension(imax) :: rat, phitvd, fluxhi, fluxlo
real, dimension(imax,0:kl) :: delt, fluxh

if ( size(tarr,1)<imax ) then
  write(6,*) "ERROR: Argument array is too small for vadv_work"
  call ccmpi_abort(-1)
end if

if ( size(tarr,2)/=kl ) then
  write(6,*) "ERROR: Argument array has the wrong number of vertical levels for vadv_work"
  call ccmpi_abort(-1)
end if

! The first sub-step is vectorised for all points - MJT

!     fluxh(k) is located at level k+.5
fluxh(:,0)  = 0.
fluxh(:,kl) = 0.

do k = 1,kl-1
  delt(:,k) = tarr(1:imax,k+1) - tarr(1:imax,k)
end do  
delt(:,kl)     = 0.     ! for T,u,v
delt(:,0)      = min(delt(:,1), tarr(1:imax,1))       ! for non-negative tt

do k = 1,kl-1  ! for fluxh at interior (k + 1/2)  half-levels
  where ( sdot(:,k+1)>0. )
    rat(1:imax) = delt(:,k-1)/(delt(:,k)+sign(1.e-20,delt(:,k)))
    fluxlo(1:imax) = tarr(1:imax,k)
  elsewhere
    rat(1:imax) = delt(:,k+1)/(delt(:,k)+sign(1.e-20,delt(:,k)))
    fluxlo(1:imax) = tarr(1:imax,k+1)
  end where
  phitvd(:) = max(0., min(2.*rat(:),.5+.5*rat(:), 2.))    ! 0 for -ve rat
  ! higher order scheme
  fluxhi(:) = rathb(k)*tarr(1:imax,k) + ratha(k)*tarr(1:imax,k+1) - .5*delt(:,k)*tfact(:)*sdot(:,k+1)
  fluxh(:,k) = sdot(:,k+1)*(fluxlo(:)+phitvd(:)*(fluxhi(:)-fluxlo(:)))
enddo      ! k loop
do k = 1,kl
  tarr(1:imax,k) = tarr(1:imax,k) + tfact(:)*(fluxh(:,k-1)-fluxh(:,k)+tarr(1:imax,k)*(sdot(:,k+1)-sdot(:,k)))
end do

! Subsequent substeps if needed.  This is fairly rare so we perform this calculation point-by-point - MJT
do iq = 1,imax      
  do i = 2,nits(iq)
    do k = 1,kl-1
      delt(iq,k) = tarr(iq,k+1) - tarr(iq,k)
    end do     ! k loop
    delt(iq,kl) = 0.
    delt(iq,0) = min(delt(iq,1), tarr(iq,1))       ! for non-negative tt
    do k = 1,kl-1  ! for fluxh at interior (k + 1/2)  half-levels
      kp = nint(sign(1.,sdot(iq,k+1)))
      kx = k + (1-kp)/2 !  k for sdot +ve,  k+1 for sdot -ve
      rat(iq) = delt(iq,k-kp)/(delt(iq,k)+sign(1.e-20,delt(iq,k)))
      fluxlo(iq) = tarr(iq,kx)
      phitvd(iq) = max(0., min(2.*rat(iq), .5+.5*rat(iq), 2.))   ! 0 for -ve rat
      ! higher order scheme
      fluxhi(iq) = rathb(k)*tarr(iq,k) + ratha(k)*tarr(iq,k+1) - .5*delt(iq,k)*tfact(iq)*sdot(iq,k+1)
      fluxh(iq,k) = sdot(iq,k+1)*(fluxlo(iq)+phitvd(iq)*(fluxhi(iq)-fluxlo(iq)))
    end do ! k
    do k = 1,kl
      tarr(iq,k) = tarr(iq,k) + tfact(iq)*(fluxh(iq,k-1)-fluxh(iq,k)+tarr(iq,k)*(sdot(iq,k+1)-sdot(iq,k)))
    end do
  end do   ! i
end do     ! iq

return
end subroutine vadv_work2
      
end module vadv
