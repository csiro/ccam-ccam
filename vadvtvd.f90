module vadv
      
private
public vadvtvd
      
contains

subroutine vadvtvd(tarr,uarr,varr,nvadh_pass,nits)   ! globpea  version
!                              vadvbott & vadvyu at bottom
!     can show adding tbar has no effect
use aerosolldr
use arrays_m
use cc_mpi
use cfrac_m
use cloudmod
use diag_m
use liqwpar_m  ! ifullw
use map_m
use nharrs_m
use sigs_m
use tkeeps, only : tke,eps
use tracers_m
use vvel_m
use xarrs_m
implicit none
include 'newmpar.h'
!     split vertical advection routine; tvd scheme; used with nonlin or upglobal
!     In flux limiter, assuming zero gradient for all top and bottom
!     variables; except extrap at bottom for qg and trace gases  Thu  06-19-1997
include 'kuocom.h'     ! also with kbsav,ktsav
include 'parm.h'
include 'parmdyn.h'
real, dimension(:,:), intent(inout) :: tarr,uarr,varr
real, dimension(ifull,kl) :: dum
real, dimension(ifull) :: tfact
integer npslx
integer ntr,k
integer, dimension(ifull) :: nvadh_pass, nits
integer, dimension(ifull,kl-1) :: kp,kx
integer, save :: num = 0
parameter (npslx=1)  ! 0 off, 1 on for nvad=-4

call START_LOG(vadv_begin)

tfact=1./real(nvadh_pass)

if(num==0)then
  num=1
  if(mydiag)then
    write(6,*) 'In vadvtvd nvadh_pass,npslx ',nvadh_pass(idjd),npslx
    write(6,*) 'tfact ',tfact(idjd)
  endif
endif

do k=1,kl-1
  kp(:,k)=sign(1.,sdot(:,k+1))
  kx(:,k)=k+(1-kp(:,k))/2 !  k for sdot +ve,  k+1 for sdot -ve
end do

!     t
call vadv_work(tarr,tfact,nits,kp,kx)
if( (diag.or.nmaxpr==1) .and. mydiag )then
!       These diagnostics don't work with single input/output argument
  write (6,"('tout',9f8.2/4x,9f8.2)") (tarr(idjd,k),k=1,kl)
  write (6,"('t#  ',9f8.2)") diagvals(tarr(:,nlv)) 
endif

!     u
call vadv_work(uarr,tfact,nits,kp,kx)
if( diag .and. mydiag )then
  write (6,"('uout',9f8.2/4x,9f8.2)") (uarr(idjd,k),k=1,kl)
  write (6,"('u#  ',9f8.2)") diagvals(uarr(:,nlv)) 
endif

!     v
call vadv_work(varr,tfact,nits,kp,kx)
if( diag .and. mydiag )then
  write (6,"('vout',9f8.2/4x,9f8.2)") (varr(idjd,k),k=1,kl)
  write (6,"('v#  ',9f8.2)") diagvals(varr(:,nlv)) 
endif

!     h_nh
if(nh/=0.and.npslx==1)then
  call vadv_work(h_nh,tfact,nits,kp,kx)
endif

!     pslx
if(npslx==1)then  ! handles -9 too
  call vadv_work(pslx,tfact,nits,kp,kx)
endif

if(mspec==1)then   ! advect qg and gases after preliminary step

!      qg
 call vadv_work(qg,tfact,nits,kp,kx)
 if( diag .and. mydiag )then
  write (6,"('qout',9f8.2/4x,9f8.2)") (1000.*qg(idjd,k),k=1,kl)
  write (6,"('qg# ',9f8.2)") diagvals(qg(:,nlv)) 
 endif

 if(ldr/=0)then
  call vadv_work(qlg,tfact,nits,kp,kx)
  call vadv_work(qfg,tfact,nits,kp,kx)
  call vadv_work(qrg,tfact,nits,kp,kx)
  call vadv_work(rfrac,tfact,nits,kp,kx)
  if (ncloud>=3) then
    call vadv_work(stratcloud,tfact,nits,kp,kx)
  end if
  if( diag .and. mydiag )then
   write (6,"('lout',9f8.2/4x,9f8.2)") (1000.*qlg(idjd,k),k=1,kl)
   write (6,"('qlg#',9f8.2)") diagvals(qlg(:,nlv)) 
   write (6,"('fout',9f8.2/4x,9f8.2)") (1000.*qfg(idjd,k),k=1,kl)
   write (6,"('qfg#',9f8.2)") diagvals(qfg(:,nlv)) 
  endif
 endif      ! if(ldr.ne.0)

 if(nvmix==6)then
  call vadv_work(eps,tfact,nits,kp,kx)
  call vadv_work(tke,tfact,nits,kp,kx)
 endif      ! if(nvmix.eq.6)

 if (abs(iaero)==2) then
  do ntr=1,naero
   call vadv_work(xtg(:,:,ntr),tfact,nits,kp,kx)
  end do
 end if

 if(ngas>0.or.nextout>=4)then
  do ntr=1,ntrac
   call vadv_work(tr(:,:,ntr),tfact,nits,kp,kx)
  enddo      ! ntr loop
 endif      ! (nextout>=4)
 
endif       ! if(mspec==1)

call END_LOG(vadv_end)
 
return
end subroutine vadvtvd
      
! Subroutine to perform generic TVD advection
subroutine vadv_work(tarr,tfact,nits,kp,kx)
      
use sigs_m
use vvel_m
      
implicit none
      
include 'newmpar.h'
include 'parmvert.h'
      
integer, dimension(ifull), intent(in) :: nits
integer, dimension(ifull,kl-1), intent(in) :: kp,kx
integer i,k,iq
real, dimension(ifull), intent(in) :: tfact
real, dimension(ifull) :: hdsdot
real, dimension(ifull,kl-1) :: rat,phitvd,fluxhi,fluxlo
real, dimension(:,:), intent(inout) :: tarr
real, dimension(ifull,0:kl) :: delt,fluxh
real, dimension(ifull,kl) :: xin

! The first sub-step is vectorised for all points - MJT

!     fluxh(k) is located at level k+.5
fluxh(:,0)=0.
fluxh(:,kl)=0.
      
delt(:,kl)=0.     ! for T,u,v
delt(:,1:kl-1)=tarr(1:ifull,2:kl)-tarr(1:ifull,1:kl-1)
delt(:,0)=min(delt(:,1),tarr(1:ifull,1))       ! for non-negative tt
do k=1,kl-1  ! for fluxh at interior (k + 1/2)  half-levels
 do iq=1,ifull
  rat(iq,k)=delt(iq,k-kp(iq,k))/(delt(iq,k)+sign(1.e-20,delt(iq,k)))
 end do
end do

phitvd=max(0.,min(2.*rat,.5+.5*rat,2.))    ! 0 for -ve rat

! higher order scheme
do k=1,kl-1
  fluxhi(:,k)=rathb(k)*tarr(1:ifull,k)+ratha(k)*tarr(1:ifull,k+1)-.5*delt(:,k)*tfact(:)*sdot(:,k+1)
end do

do k=1,kl-1
  do iq=1,ifull
    fluxlo(iq,k)=tarr(iq,kx(iq,k))
 end do
 fluxh(:,k)=sdot(:,k+1)*(fluxlo(:,k)+phitvd(:,k)*(fluxhi(:,k)-fluxlo(:,k)))
enddo     ! k loop

do k=1,kl
  tarr(1:ifull,k)=tarr(1:ifull,k)+tfact(:)*(fluxh(:,k-1)-fluxh(:,k)+tarr(1:ifull,k)*(sdot(:,k+1)-sdot(:,k)))
enddo      ! k loop

! Subsequent substeps if needed.  This is fairly rare so we perform this calculation on
! a point-by-point basis - MJT

do iq=1,ifull      
  do i=2,nits(iq)
    do k=1,kl-1
      delt(iq,k)=tarr(iq,k+1)-tarr(iq,k)
    enddo     ! k loop
    delt(iq,0)=min(delt(iq,1),tarr(iq,1))       ! for non-negative tt
    do k=1,kl-1  ! for fluxh at interior (k + 1/2)  half-levels
      rat(iq,k)=delt(iq,k-kp(iq,k))/(delt(iq,k)+sign(1.e-20,delt(iq,k)))
    end do
    phitvd(iq,:)=max(0.,min(2.*rat(iq,:),.5+.5*rat(iq,:),2.))             ! 0 for -ve rat
    ! higher order scheme
    fluxhi(iq,1:kl-1)=rathb(1:kl-1)*tarr(iq,1:kl-1)+ratha(1:kl-1)*tarr(iq,2:kl)-.5*delt(iq,1:kl-1)*tfact(iq)*sdot(iq,2:kl)
    do k=1,kl-1
      fluxlo(iq,k)=tarr(iq,kx(iq,k))
    end do
    fluxh(iq,1:kl-1)=sdot(iq,2:kl)*(fluxlo(iq,1:kl-1)+phitvd(iq,1:kl-1)*(fluxhi(iq,1:kl-1)-fluxlo(iq,1:kl-1)))
    tarr(iq,1:kl)=tarr(iq,1:kl)+tfact(iq)*(fluxh(iq,0:kl-1)-fluxh(iq,1:kl)+tarr(iq,1:kl)*(sdot(iq,2:kl+1)-sdot(iq,1:kl)))
  end do ! i
end do ! iq

return
end subroutine vadv_work
      
end module vadv
