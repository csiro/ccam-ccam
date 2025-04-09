! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2025 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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

!     split vertical advection routine; tvd scheme; used with nonlin or upglobal
!     In flux limiter, assuming zero gradient for all top and bottom
!     variables; except extrap at bottom for qg and trace gases  Thu  06-19-1997
    
module vadv
      
private
public vadvtvd
      
contains

subroutine vadvtvd(tarr,uarr,varr,nvadh_inv_pass,nits)

use aerosol_arrays
use arrays_m
use cc_mpi
use cfrac_m, only : stratcloud
use diag_m
use kuocom_m
use liqwpar_m  ! ifullw
use map_m
use newmpar_m
use nharrs_m
use parm_m
use parmdyn_m
use sigs_m
use tkeeps, only : tke,eps
use tracers_m
use vvel_m
use xarrs_m

implicit none

integer ntr,k
integer, dimension(ifull), intent(in) :: nits
integer, save :: num = 0
real, dimension(:,:), intent(inout) :: tarr,uarr,varr
real, dimension(ifull), intent(in) :: nvadh_inv_pass
real, dimension(ifull,kl,5) :: tdum
real, dimension(ifull,kl,5) :: qdum
real, dimension(ifull,kl,2) :: edum

call START_LOG(vadv_begin)

if ( num==0 ) then
  num = 1
  if ( mydiag ) then
    write(6,*) 'In vadvtvd nvadh_pass ',nint(1./nvadh_inv_pass(idjd))
  end if
end if

tdum(1:ifull,1:kl,1) = tarr(1:ifull,1:kl)
tdum(1:ifull,1:kl,2) = uarr(1:ifull,1:kl)
tdum(1:ifull,1:kl,3) = varr(1:ifull,1:kl)
tdum(1:ifull,1:kl,4) = pslx(1:ifull,1:kl)
if ( nh/=0 ) then
  tdum(1:ifull,1:kl,5) = h_nh(1:ifull,1:kl)
end if  

if ( mspec==1 ) then

  qdum(1:ifull,1:kl,1) = qg(1:ifull,1:kl)
  if ( ldr/=0 ) then
    qdum(1:ifull,1:kl,2) = qlg(1:ifull,1:kl)
    qdum(1:ifull,1:kl,3) = qfg(1:ifull,1:kl)
    qdum(1:ifull,1:kl,4) = stratcloud(1:ifull,1:kl)
    if ( ncloud>=100 .and. ncloud<200 ) then
      qdum(1:ifull,1:kl,5) = ni(1:ifull,1:kl)
    end if  
  end if  

  if ( nvmix==6 .or. nvmix==9 ) then
    edum(1:ifull,1:kl,1) = eps(1:ifull,1:kl)
    edum(1:ifull,1:kl,2) = tke(1:ifull,1:kl)
  end if  

end if


#ifdef GPU
!$acc data create(sdot,nvadh_inv_pass,nits,ratha,rathb)
!$acc update device(sdot,nvadh_inv_pass,nits,ratha,rathb)
#endif


if ( nh/=0 ) then
  call vadv_work(tdum(:,:,1:5),nvadh_inv_pass,nits)
else
  call vadv_work(tdum(:,:,1:4),nvadh_inv_pass,nits)  
end if

!       These diagnostics don't work with single input/output argument
if( diag .and. mydiag )then
  write (6,"('tout',9f8.2/4x,9f8.2)") (tarr(idjd,k),k=1,kl)
  write (6,"('t#  ',9f8.2)") diagvals(tarr(:,nlv)) 
  write (6,"('uout',9f8.2/4x,9f8.2)") (uarr(idjd,k),k=1,kl)
  write (6,"('u#  ',9f8.2)") diagvals(uarr(:,nlv)) 
  write (6,"('vout',9f8.2/4x,9f8.2)") (varr(idjd,k),k=1,kl)
  write (6,"('v#  ',9f8.2)") diagvals(varr(:,nlv)) 
endif

if ( mspec==1 ) then   ! advect qg and gases after preliminary step

  if ( ldr/=0 .and. ncloud>=100 .and. ncloud<200 ) then  
    ! Lin microphysics  
    call vadv_work(qdum(:,:,1:5),nvadh_inv_pass,nits)
  else if ( ldr/=0 ) then
    ! Leon microphysics  
    call vadv_work(qdum(:,:,1:4),nvadh_inv_pass,nits)  
  else
    ! diagnostic  
    call vadv_work(qdum(:,:,1:1),nvadh_inv_pass,nits)
  end if  
    
  if ( diag .and. mydiag ) then
    write (6,"('qout',9f8.2/4x,9f8.2)") (1000.*qg(idjd,k),k=1,kl)
    write (6,"('qg# ',9f8.2)") diagvals(qg(:,nlv)) 
    if (  ldr/=0 ) then
      write (6,"('lout',9f8.2/4x,9f8.2)") (1000.*qlg(idjd,k),k=1,kl)
      write (6,"('qlg#',9f8.2)") diagvals(qlg(:,nlv)) 
      write (6,"('fout',9f8.2/4x,9f8.2)") (1000.*qfg(idjd,k),k=1,kl)
      write (6,"('qfg#',9f8.2)") diagvals(qfg(:,nlv))
    end if  
  end if

  if ( nvmix==6 .or. nvmix==9 ) then
    call vadv_work(edum(:,:,1:2),nvadh_inv_pass,nits)
  end if      ! if(nvmix==6 .or. nvmix==9 )

  if ( abs(iaero)>=2 .and. nhstest>=0 ) then
    call vadv_work(xtg(:,:,1:naero),nvadh_inv_pass,nits)
  end if   ! abs(iaero)>=2
  
  if ( ntrac>0 ) then
    call vadv_work(tr(:,:,1:ntrac),nvadh_inv_pass,nits)
  end if        ! (nextout>=4)

end if          ! if(mspec==1)


#ifdef GPU
!$acc wait
!$acc end data
#endif


tarr(1:ifull,1:kl) = tdum(1:ifull,1:kl,1)
uarr(1:ifull,1:kl) = tdum(1:ifull,1:kl,2)
varr(1:ifull,1:kl) = tdum(1:ifull,1:kl,3)
pslx(1:ifull,1:kl) = tdum(1:ifull,1:kl,4)
if ( nh/=0 ) then
  h_nh(1:ifull,1:kl) = tdum(1:ifull,1:kl,5)
end if  

if ( mspec==1 ) then

  qg(1:ifull,1:kl) = qdum(1:ifull,1:kl,1)
  if ( ldr/=0 ) then
    qlg(1:ifull,1:kl) = qdum(1:ifull,1:kl,2)
    qfg(1:ifull,1:kl) = qdum(1:ifull,1:kl,3)
    stratcloud(1:ifull,1:kl) = qdum(1:ifull,1:kl,4)
    if ( ncloud>=100 .and. ncloud<200 ) then
      ni(1:ifull,1:kl) = qdum(1:ifull,1:kl,5)
    end if
  end if  

  if ( nvmix==6 .or. nvmix==9 ) then
    eps(1:ifull,1:kl) = edum(1:ifull,1:kl,1)
    tke(1:ifull,1:kl) = edum(1:ifull,1:kl,2)
  end if

end if  

call END_LOG(vadv_end)
 
return
end subroutine vadvtvd
      
! Subroutine to perform generic TVD advection
subroutine vadv_work(tarr,nvadh_inv_pass,nits)

use cc_acc, only : async_length
use cc_mpi, only : ccmpi_abort
use newmpar_m
use parmvert_m
use sigs_m
use vvel_m
      
implicit none
      
integer, dimension(ifull), intent(in) :: nits
integer i, k, iq, kp, kx, n, ntr
integer, save :: async_counter = -1
real, dimension(:,:,:), intent(inout) :: tarr
real, dimension(ifull), intent(in) :: nvadh_inv_pass
real, dimension(ifull,0:kl,size(tarr,3)) :: delt, fluxh
real rat, phitvd, fluxhi, fluxlo

async_counter = mod(async_counter+1, async_length)
ntr = size(tarr,3)

! The first sub-step is vectorised for all points.

#ifdef GPU
!$acc enter data create(tarr,delt,fluxh) async(async_counter)
!$acc update device(tarr) async(async_counter)
#else
!$omp parallel   
#endif
do i = 1,maxval(nits(1:ifull))
#ifdef GPU
  !$acc parallel loop collapse(3) present(delt,tarr) async(async_counter)
#else
  !$omp do schedule(static) private(n,k,iq)
#endif
  do n = 1,ntr
    do k = 1,kl-1
      do iq = 1,ifull
        delt(iq,k,n) = tarr(iq,k+1,n) - tarr(iq,k,n)
      end do
    end do
  end do
#ifdef GPU
  !$acc end parallel loop
  !$acc parallel loop collapse(2) present(fluxh,delt) async(async_counter)
#else
  !$omp end do nowait
  !$omp do schedule(static) private(n,iq)
#endif
  do n = 1,ntr
    do iq = 1,ifull
      fluxh(iq,0,n)  = 0.
      fluxh(iq,kl,n) = 0.
      delt(iq,kl,n)  = 0.     ! for T,u,v
      delt(iq,0,n)   = 0.
    end do
  end do
#ifdef GPU
  !$acc end parallel loop
#else
    !$omp end do nowait
#endif
  if ( ntvd==2 ) then ! MC
#ifdef GPU
    !$acc parallel loop collapse(3) present(sdot,delt,tarr,ratha,rathb,nvadh_inv_pass,fluxh) async(async_counter)
#else
    !$omp do schedule(static) private(n,k,iq,kp,kx,rat,phitvd,fluxlo,fluxhi)
#endif
    do n = 1,ntr
      do k = 1,kl-1  ! for fluxh at interior (k + 1/2)  half-levels
        do iq = 1,ifull      
          kp = nint(sign(1.,sdot(iq,k+1)))
          kx = k + (1-kp)/2 !  k for sdot +ve,  k+1 for sdot -ve
          rat = delt(iq,k-kp,n)/(delt(iq,k,n)+sign(1.e-20,delt(iq,k,n)))
          phitvd = max(0., min(2.*rat,.5+.5*rat, 2.))    ! 0 for -ve rat
          fluxlo = tarr(iq,kx,n)
          ! higher order scheme
          fluxhi = rathb(k)*tarr(iq,k,n) + ratha(k)*tarr(iq,k+1,n) - .5*delt(iq,k,n)*sdot(iq,k+1)*nvadh_inv_pass(iq)
          fluxh(iq,k,n) = sdot(iq,k+1)*(fluxlo+phitvd*(fluxhi-fluxlo))
        end do
      end do
    end do      ! k loop
#ifdef GPU
    !$acc end parallel loop
#else
    !$omp end do nowait
#endif
  else if ( ntvd==3 ) then ! Superbee
#ifdef GPU
    !$acc parallel loop collapse(3) present(delt,tarr,sdot,rathb,ratha,nvadh_inv_pass,fluxh) async(async_counter)
#else
    !$omp do schedule(static) private(n,iq,k,kp,kx,rat,phitvd,fluxlo,fluxhi)
#endif
    do n = 1,ntr
      do k = 1,kl-1  ! for fluxh at interior (k + 1/2)  half-levels
        do iq = 1,ifull  
          kp = nint(sign(1.,sdot(iq,k+1)))
          kx = k + (1-kp)/2 !  k for sdot +ve,  k+1 for sdot -ve
          rat = delt(iq,k-kp,n)/(delt(iq,k,n)+sign(1.e-20,delt(iq,k,n)))
          phitvd = max(0.,min(1.,2.*rat),min(2.,rat)) ! 0 for -ve rat
          fluxlo = tarr(iq,kx,n)
          ! higher order scheme
          fluxhi = rathb(k)*tarr(iq,k,n) + ratha(k)*tarr(iq,k+1,n) - .5*delt(iq,k,n)*sdot(iq,k+1)*nvadh_inv_pass(iq)
          fluxh(iq,k,n) = sdot(iq,k+1)*(fluxlo+phitvd*(fluxhi-fluxlo))
        end do  
      end do ! k
    end do
#ifdef GPU
    !$acc end parallel loop
#else
    !$omp end do nowait
#endif
  else
    write(6,*) "ERROR: Unknown option ntvd = ",ntvd
    call ccmpi_abort(-1)
  end if
#ifdef GPU
  !$acc parallel loop collapse(3) present(nits,fluxh,tarr,sdot,nvadh_inv_pass) async(async_counter)
#else
  !$omp do schedule(static) private(n,k,iq)
#endif
  do n = 1,ntr
    do k = 1,kl
      do iq = 1,ifull
        if ( i<=nits(iq) ) then
          tarr(iq,k,n) = tarr(iq,k,n) + (fluxh(iq,k-1,n)-fluxh(iq,k,n) &
                             +tarr(iq,k,n)*(sdot(iq,k+1)-sdot(iq,k)))*nvadh_inv_pass(iq)
        end if  
      end do
    end do
  end do
#ifdef GPU
  !$acc end parallel loop
#else
  !$omp end do nowait
#endif
end do ! i=1,nits
#ifdef GPU
!$acc update self(tarr) async(async_counter)
!$acc exit data delete(tarr,delt,fluxh) async(async_counter)
#else
!$omp end parallel
#endif

return
end subroutine vadv_work

end module vadv
