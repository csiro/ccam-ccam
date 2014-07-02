module mgsolve

! This module solves for the atmosphere helmholtz equation and the ocean free
! surface equation using a multi-grid approach.

! Note that the solver can be greatly optimised for grids that permit many
! local subdivisions.  For example if the number of grid points on a processor
! is M x M, where M = 2^N and N is an integer greater than 1, then the
! multi-grid solver will be considerably more efficent as more processors can
! be applied in parallel.  Example grids that typically work well include
! C48, C96, C192, C384, C768 and C1536.  An alternate family of grids includes
! C32, C64, C128, C256, C512 and C1024.

! Using the idleproc preprocessor directive (e.g., -Didleproc) removes redundant
! work, but leaves processors idle.  However, the option is faster because it
! avoids redundant messages and helps reduce the amount of bandwidth employed
! by the scheme.

! Note that k is the innermost array index to improve vectorization of the code

implicit none

private
public mghelm,mgmlo,mgsor_init,mgzz_init

integer, save :: mg_maxsize,mg_minsize,gmax

integer, parameter :: itr_mg   =20 ! maximum number of iterations for atmosphere MG solver
integer, parameter :: itr_mgice=20 ! maximum number of iterations for ocean/ice MG solver
integer, parameter :: itrbgn   =1  ! number of iterations relaxing the solution after restriction
integer, parameter :: itrend   =1  ! number of iterations relaxing the solution after interpolation

real, parameter :: dfac=0.25       ! adjustment for grid spacing after restriction

logical, save :: sorfirst=.true.
logical, save :: zzfirst =.true.
logical, save :: mlofirst=.true.


contains

! This version is for the atmosphere
subroutine mghelm(izz,izzn,izze,izzw,izzs,ihelm,iv,jrhs)

use cc_mpi
use indices_m

implicit none

include 'newmpar.h'
include 'parm.h'
include 'parmdyn.h'

integer, dimension(kl) :: iters
integer, dimension(mg_minsize,kl) :: indy
integer itrc,itr,ng,ng4,g,gb,k,jj,i,j,iq,ng_x
integer klimc,knew,klim,ir,ic
integer nc,n,iq_a,iq_b,iq_c,iq_d
real, dimension(ifull+iextra,kl), intent(inout) :: iv
real, dimension(ifull,kl), intent(in) :: ihelm,jrhs
real, dimension(ifull,kl) :: irhs
real, dimension(ifull), intent(in) :: izz,izzn,izze,izzw,izzs
real, dimension(ifullmaxcol,kl,maxcolour) :: rhelmc,rhsc
real, dimension(ifullmaxcol,maxcolour) :: zznc,zzec,zzwc,zzsc
real, dimension(ifull+iextra,kl) :: vdum
real, dimension(kl,mg_maxsize,gmax+1) :: v
real, dimension(2*kl,mg_maxsize,2:gmax+1) :: rhs
real, dimension(kl,mg_maxsize,gmax+1) :: helm
real, dimension(2*kl,mg_maxsize) :: w
real, dimension(mg_maxsize,kl) :: dsol
real, dimension(mg_minsize,mg_minsize,kl) :: helm_o
real, dimension(mg_minsize) :: v_o
real, dimension(2*kl,2) :: smaxmin_g
real, dimension(kl) :: dsolmax_g,savg,sdif,dsolmax_l

call START_LOG(helm_begin)

if (sorfirst.or.zzfirst) then
  write(6,*) "ERROR: mghelm requires mgsor_init and mgzz_init to be called first"
  call ccmpi_abort(-1)
end if

! zz*(DIV^2 v) - helm*v = rhs
! zz*(DIV^2 vd+v0) - helm*(vd+v0) = rhs
! zz*(DIV^2 vd) - helm*vd = rhs + helm*v0 - zz*(DIV^2 v0)

! MJT notes - to remove an additional syncronisation step we update rhs and helm
! parameters during the first iteration of the solution.  Effectively the mgsetup
! stage also becomes the first iteration of the solution.

call START_LOG(mgsetup_begin)

iv(ifull+1:ifull+iextra,:)=0. ! for IBM compiler

! solver assumes boundaries have been updated
call bounds(iv)

! determine max/min for convergence calculations
klim=kl
vdum=0.
smaxmin_g=0.
do k=1,kl
  smaxmin_g(k,1)=maxval(iv(1:ifull,k))
  smaxmin_g(k,2)=minval(iv(1:ifull,k))
end do

! pack colour arrays at fine level
do nc=1,maxcolour
  zznc(1:ifullcol(nc),nc) =izzn(iqx(1:ifullcol(nc),nc))
  zzwc(1:ifullcol(nc),nc) =izzw(iqx(1:ifullcol(nc),nc))
  zzec(1:ifullcol(nc),nc) =izze(iqx(1:ifullcol(nc),nc))
  zzsc(1:ifullcol(nc),nc) =izzs(iqx(1:ifullcol(nc),nc))
  do k=1,kl
    rhelmc(1:ifullcol(nc),k,nc)=1./(ihelm(iqx(1:ifullcol(nc),nc),k)-izz(iqx(1:ifullcol(nc),nc)))
    rhsc(1:ifullcol(nc),k,nc)  =jrhs(iqx(1:ifullcol(nc),nc),k)
  end do
end do

! Before sending convegence testing data in smaxmin_g and ihelm weights, we perform one iteration of the solver
! that can be updated with the smaxmin_g and ihelm arrays

! update on model grid using colours
do nc=1,maxcolour
  do k=1,kl
    dsol(iqx(1:ifullcol(nc),nc),k) = ( zznc(1:ifullcol(nc),nc)*iv(iqn(1:ifullcol(nc),nc),k)    &
                                     + zzwc(1:ifullcol(nc),nc)*iv(iqw(1:ifullcol(nc),nc),k)    &
                                     + zzec(1:ifullcol(nc),nc)*iv(iqe(1:ifullcol(nc),nc),k)    &
                                     + zzsc(1:ifullcol(nc),nc)*iv(iqs(1:ifullcol(nc),nc),k)    &
                                     - rhsc(1:ifullcol(nc),k,nc) )*rhelmc(1:ifullcol(nc),k,nc) &
                                     - iv(iqx(1:ifullcol(nc),nc),k)
    iv(iqx(1:ifullcol(nc),nc),k) = iv(iqx(1:ifullcol(nc),nc),k) + dsol(iqx(1:ifullcol(nc),nc),k)
  end do
  call bounds_colour(iv,nc)
end do
  
! residual
do k=1,kl
  w(k,1:ifull)=-izzn*iv(in,k)-izzw*iv(iw,k)-izze*iv(ie,k)-izzs*iv(is,k)+jrhs(:,k)+iv(1:ifull,k)*(ihelm(:,k)-izz)
end do
! also include ihelm weights for upscaled grid
w(kl+1:2*kl,1:ifull)=transpose(ihelm(1:ifull,1:kl))
  
! For when the inital grid cannot be upscaled - note helm and smaxmin_g are also included
call mgcollect(1,w,smaxmin_g)
helm(1:kl,1:mg(1)%ifull,1)=w(kl+1:2*kl,1:mg(1)%ifull)

do g=1,min(mg_maxlevel_local,1) ! same as if (mg_maxlevel_local>0) then ...

  ! restriction
  ! (since this always operates within a panel, then ine = ien is always true)
  ng4=mg(1)%ifull_fine
  rhs(1:2*kl,1:ng4,2)=0.25*(w(1:2*kl,mg(1)%fine  )+w(1:2*kl,mg(1)%fine_n )  &
                           +w(1:2*kl,mg(1)%fine_e)+w(1:2*kl,mg(1)%fine_ne))
                             
  ! merge grids if insufficent points on this processor - note helm and smaxmin_g are also included
  call mgcollect(2,rhs(:,:,2),smaxmin_g)
  helm(1:kl,1:mg(2)%ifull,2)=rhs(kl+1:2*kl,1:mg(2)%ifull,2)
  
end do

! upscale grid
do g=2,gmax
  
  ng=mg(g)%ifull
                
  ! update scalar field
  ! assume zero for first guess of residual (also avoids additional bounds call)
  !v(1:klim,1:ng,g)=0.
  do iq=1,ng
    v(1:kl,iq,g)=-rhs(1:kl,iq,g)/(helm(1:kl,iq,g)-mg(g)%zz(iq))
  end do
  call mgbounds(g,v(:,:,g))

  do i=1,itrbgn-1
    do iq=1,ng
      ! post smoothing
      w(1:kl,iq)=(mg(g)%zze(iq)*v(1:kl,mg(g)%ie(iq),g)+mg(g)%zzw(iq)*v(1:kl,mg(g)%iw(iq),g) &
                 +mg(g)%zzn(iq)*v(1:kl,mg(g)%in(iq),g)+mg(g)%zzs(iq)*v(1:kl,mg(g)%is(iq),g) &
                 -rhs(1:kl,iq,g))/(helm(1:kl,iq,g)-mg(g)%zz(iq))
    end do
    call mgbounds(g,w(:,:),klim=klim)
    v(1:kl,1:mg(g)%ifull+mg(g)%iextra,g)=w(1:kl,1:mg(g)%ifull+mg(g)%iextra)
  end do
  
  do iq=1,ng
    ! residual
    w(1:kl,iq)=-mg(g)%zze(iq)*v(1:kl,mg(g)%ie(iq),g)-mg(g)%zzw(iq)*v(1:kl,mg(g)%iw(iq),g) &
               -mg(g)%zzn(iq)*v(1:kl,mg(g)%in(iq),g)-mg(g)%zzs(iq)*v(1:kl,mg(g)%is(iq),g) &
               +rhs(1:kl,iq,g)+(helm(1:kl,iq,g)-mg(g)%zz(iq))*v(1:kl,iq,g)
    ! upscale helm weights
    w(kl+1:2*kl,iq)=helm(1:kl,iq,g)
  end do

  do iq=1,mg(g)%ifull_fine
    ! restriction
    ! (calculate coarser grid before mgcollect as more work is done in parallel)
    rhs(1:2*kl,iq,g+1)=0.25*(w(1:2*kl,mg(g)%fine(iq)  )+w(1:2*kl,mg(g)%fine_n(iq) )  &
                            +w(1:2*kl,mg(g)%fine_e(iq))+w(1:2*kl,mg(g)%fine_ne(iq)))
  end do

  ! merge grids if insufficent points on this processor - note helm and smaxmin_g are also included
  call mgcollect(g+1,rhs(:,:,g+1),smaxmin_g)
  helm(1:kl,1:mg(g+1)%ifull,g+1)=rhs(kl+1:2*kl,1:mg(g+1)%ifull,g+1)

end do

! store data for LU decomposition of coarse grid
do g=mg_maxlevel,mg_maxlevel_local ! same as if (mg_maxlevel_local==mg_maxlevel) then ...
  helm_o(:,:,:)=0.
  ! solve coarse grid
  ng=mg(g)%ifull
  do k=1,kl
    do iq=1,ng
      helm_o(mg(g)%in(iq),iq,k)=mg(g)%zzn(iq)
      helm_o(mg(g)%is(iq),iq,k)=mg(g)%zzs(iq)
      helm_o(mg(g)%ie(iq),iq,k)=mg(g)%zze(iq)
      helm_o(mg(g)%iw(iq),iq,k)=mg(g)%zzw(iq)
      helm_o(iq,iq,k)=mg(g)%zz(iq)-helm(k,iq,g)
    end do
    call mdecomp(helm_o(:,:,k),indy(:,k))
    ! perform LU decomposition and back substitute with RHS
    ! to solve for v on coarse grid
    v_o(1:ng)=rhs(k,1:ng,g)
    call mbacksub(helm_o(:,:,k),v_o(1:ng),indy(:,k))
    v(k,1:ng,g)=v_o(1:ng)
  end do
end do

! downscale grid
do g=gmax,2,-1

  ! send coarse grid solution to processors and also bcast global smaxmin_g
  call mgbcastxn(g+1,v(:,:,g+1),smaxmin_g)

  do iq=1,mg(g+1)%ifull_coarse
    ! interpolation
    w(1:kl,iq)= mg(g+1)%wgt_a(iq)*v(1:kl,mg(g+1)%coarse_a(iq),g+1) + mg(g+1)%wgt_bc(iq)*v(1:kl,mg(g+1)%coarse_b(iq),g+1) &
             + mg(g+1)%wgt_bc(iq)*v(1:kl,mg(g+1)%coarse_c(iq),g+1) +  mg(g+1)%wgt_d(iq)*v(1:kl,mg(g+1)%coarse_d(iq),g+1)

    ! extension
    ! No mgbounds as the v halo has already been updated and
    ! the coarse interpolation also updates the w halo
    w(1:kl,iq)=v(1:kl,iq,g)+w(1:kl,iq)
  end do

  ! MJT notes - The first correction is usually sufficient to produce a good approximation
  ! to the converged solution, thereby avoiding additional calls to mgbounds
  do i=1,itrend-1
    do iq=1,mg(g)%ifull
      ! post smoothing
      v(1:kl,iq,g)=(mg(g)%zze(iq)*w(1:kl,mg(g)%ie(iq))+mg(g)%zzw(iq)*w(1:kl,mg(g)%iw(iq)) &
                   +mg(g)%zzn(iq)*w(1:kl,mg(g)%in(iq))+mg(g)%zzs(iq)*w(1:kl,mg(g)%is(iq)) &
                   -rhs(1:kl,iq,g))/(helm(1:kl,iq,g)-mg(g)%zz(iq))
    end do
    call mgbounds(g,v(:,:,g))
    w(1:kl,1:mg(g)%ifull+mg(g)%iextra)=v(1:kl,1:mg(g)%ifull+mg(g)%iextra,g)
  end do
  do iq=1,mg(g)%ifull
    ! post smoothing
    v(1:kl,iq,g)=(mg(g)%zze(iq)*w(1:kl,mg(g)%ie(iq))+mg(g)%zzw(iq)*w(1:kl,mg(g)%iw(iq)) &
                 +mg(g)%zzn(iq)*w(1:kl,mg(g)%in(iq))+mg(g)%zzs(iq)*w(1:kl,mg(g)%is(iq)) &
                 -rhs(1:kl,iq,g))/(helm(1:kl,iq,g)-mg(g)%zz(iq))
  end do    
  call mgbounds(g,v(:,:,g),corner=.true.)

end do
  
do g=1,min(mg_maxlevel_local,1) ! same as if (mg_maxlevel_local>0) then ...
    
  ! broadcast coarse solution to fine grid, as well as global smaxmin_g
  call mgbcastxn(2,v(:,:,2),smaxmin_g)

  ! interpolation
  do iq=1,mg(2)%ifull_coarse
    w(1:kl,iq)= mg(2)%wgt_a(iq)*v(1:kl,mg(2)%coarse_a(iq),2) + mg(2)%wgt_bc(iq)*v(1:kl,mg(2)%coarse_b(iq),2) &
             + mg(2)%wgt_bc(iq)*v(1:kl,mg(2)%coarse_c(iq),2) +  mg(2)%wgt_d(iq)*v(1:kl,mg(2)%coarse_d(iq),2)
  end do
    
end do

! multi-grid solver bounds indicies do not match standard iextra indicies, so we need to remap the halo
if (mg(1)%merge_len>1) then
  call mgbcastxn(1,w,smaxmin_g,klim=kl)
  ir=mod(mg(1)%merge_pos-1,mg(1)%merge_row)+1   ! index for proc row
  ic=(mg(1)%merge_pos-1)/mg(1)%merge_row+1      ! index for proc col
  do n=1,npan
    do jj=1,jpan
      iq_a=1+(jj-1)*ipan+(n-1)*ipan*jpan
      iq_b=jj*ipan+(n-1)*ipan*jpan
      iq_c=1+(ir-1)*ipan+(jj-1+(ic-1)*jpan)*ipan*mg(1)%merge_row+(n-1)*ipan*jpan*mg(1)%merge_len
      iq_d=ir*ipan+(jj-1+(ic-1)*jpan)*ipan*mg(1)%merge_row+(n-1)*ipan*jpan*mg(1)%merge_len
      vdum(iq_a:iq_b,1:kl)=transpose(w(1:kl,iq_c:iq_d))
    end do
    do i=1,ipan
      iq_a=i+(n-1)*ipan*jpan
      iq_c=i+(ir-1)*ipan+((ic-1)*jpan)*ipan*mg(1)%merge_row+(n-1)*ipan*jpan*mg(1)%merge_len
      iq_b=is(iq_a)
      iq_d=mg(1)%is(iq_c)
      vdum(iq_b,1:kl)=w(1:kl,iq_d)
      iq_a=i+(jpan-1)*ipan+(n-1)*ipan*jpan
      iq_c=i+(ir-1)*ipan+(jpan-1+(ic-1)*jpan)*ipan*mg(1)%merge_row+(n-1)*ipan*jpan*mg(1)%merge_len
      iq_b=in(iq_a)
      iq_d=mg(1)%in(iq_c)
      vdum(iq_b,1:kl)=w(1:kl,iq_d)
    end do  
    do j=1,jpan
      iq_a=1+(j-1)*ipan+(n-1)*ipan*jpan
      iq_c=1+(ir-1)*ipan+(j-1+(ic-1)*jpan)*ipan*mg(1)%merge_row+(n-1)*ipan*jpan*mg(1)%merge_len
      iq_b=iw(iq_a)
      iq_d=mg(1)%iw(iq_c)
      vdum(iq_b,1:kl)=w(1:kl,iq_d)
      iq_a=ipan+(j-1)*ipan+(n-1)*ipan*jpan
      iq_c=ipan+(ir-1)*ipan+(j-1+(ic-1)*jpan)*ipan*mg(1)%merge_row+(n-1)*ipan*jpan*mg(1)%merge_len
      iq_b=ie(iq_a)
      iq_d=mg(1)%ie(iq_c)
      vdum(iq_b,1:kl)=w(1:kl,iq_d)
    end do
  end do
else
  ! remap mg halo to normal halo
  vdum(1:ifull,1:kl)=transpose(w(1:kl,1:ifull))
  do n=0,npan-1
    do i=1,ipan
      iq=i+n*ipan*jpan
      iq_a=is(iq)
      iq_b=mg(1)%is(iq)
      vdum(iq_a,1:kl)=w(1:kl,iq_b)
      iq=i+(jpan-1)*ipan+n*ipan*jpan
      iq_a=in(iq)
      iq_b=mg(1)%in(iq)
      vdum(iq_a,1:kl)=w(1:kl,iq_b)
    end do  
    do j=1,jpan
      iq=1+(j-1)*ipan+n*ipan*jpan
      iq_a=iw(iq)
      iq_b=mg(1)%iw(iq)
      vdum(iq_a,1:kl)=w(1:kl,iq_b)
      iq=j*ipan+n*ipan*jpan
      iq_a=ie(iq)
      iq_b=mg(1)%ie(iq)
      vdum(iq_a,1:kl)=w(1:kl,iq_b)
    end do
  end do
end if
  
! extension
iv(1:ifull+iextra,1:kl)=iv(1:ifull+iextra,1:kl)+vdum(1:ifull+iextra,1:kl)
  
do i=1,itrend
  ! post smoothing
  do nc=1,maxcolour
    do k=1,kl
      dsol(iqx(1:ifullcol(nc),nc),k)=( zznc(1:ifullcol(nc),nc)*iv(iqn(1:ifullcol(nc),nc),k)   &
                                     + zzwc(1:ifullcol(nc),nc)*iv(iqw(1:ifullcol(nc),nc),k)   &
                                     + zzec(1:ifullcol(nc),nc)*iv(iqe(1:ifullcol(nc),nc),k)   &
                                     + zzsc(1:ifullcol(nc),nc)*iv(iqs(1:ifullcol(nc),nc),k)   &
                                     - rhsc(1:ifullcol(nc),k,nc))*rhelmc(1:ifullcol(nc),k,nc) &
                                     - iv(iqx(1:ifullcol(nc),nc),k)
      iv(iqx(1:ifullcol(nc),nc),k) = iv(iqx(1:ifullcol(nc),nc),k) + dsol(iqx(1:ifullcol(nc),nc),k)
    end do      
    call bounds_colour(iv,nc)
  end do
end do

! remove offsets
savg(1:kl)=0.5*(smaxmin_g(1:kl,1)+smaxmin_g(1:kl,2))
sdif(1:kl)=smaxmin_g(1:kl,1)-smaxmin_g(1:kl,2)
do k=1,kl
  iv(:,k)=iv(:,k)-savg(k)
  irhs(:,k)=jrhs(:,k)+(ihelm(:,k)-izz-izzn-izzs-izze-izzw)*savg(k)
end do

! re-pack colour arrays at fine level
do nc=1,maxcolour
  do k=1,kl
    rhsc(1:ifullcol(nc),k,nc)=irhs(iqx(1:ifullcol(nc),nc),k)
  end do
end do

call END_LOG(mgsetup_end)

! Main loop
iters=0
do itr=2,itr_mg

  call START_LOG(mgfine_begin)

  ! update on model grid using colours
  do nc=1,maxcolour
    do k=1,klim
      dsol(iqx(1:ifullcol(nc),nc),k) = ( zznc(1:ifullcol(nc),nc)*iv(iqn(1:ifullcol(nc),nc),k)    &
                                       + zzwc(1:ifullcol(nc),nc)*iv(iqw(1:ifullcol(nc),nc),k)    &
                                       + zzec(1:ifullcol(nc),nc)*iv(iqe(1:ifullcol(nc),nc),k)    &
                                       + zzsc(1:ifullcol(nc),nc)*iv(iqs(1:ifullcol(nc),nc),k)    &
                                       - rhsc(1:ifullcol(nc),k,nc) )*rhelmc(1:ifullcol(nc),k,nc) &
                                       - iv(iqx(1:ifullcol(nc),nc),k)
      iv(iqx(1:ifullcol(nc),nc),k) = iv(iqx(1:ifullcol(nc),nc),k) + dsol(iqx(1:ifullcol(nc),nc),k)
    end do
    call bounds_colour(iv,nc,klim=klim)
  end do
  
  ! test for convergence
  dsolmax_g(1:klim)=maxval(abs(dsol(1:ifull,1:klim)),dim=1)
  
  ! residual
  do k=1,klim
    w(k,1:ifull)=-izzn*iv(in,k)-izzw*iv(iw,k)-izze*iv(ie,k)-izzs*iv(is,k)+irhs(:,k)+iv(1:ifull,k)*(ihelm(:,k)-izz)
  end do
  
  ! For when the inital grid cannot be upscaled
  call mgcollect(1,w,dsolmax_g,klim=klim)

  do g=1,min(mg_maxlevel_local,1) ! same as if (mg_maxlevel_local>0) then ...

    ! restriction
    ! (since this always operates within a panel, then ine = ien is always true)
    ng4=mg(1)%ifull_fine
    rhs(1:klim,1:ng4,2)=0.25*(w(1:klim,mg(1)%fine  )+w(1:klim,mg(1)%fine_n )  &
                             +w(1:klim,mg(1)%fine_e)+w(1:klim,mg(1)%fine_ne))
                             
    ! merge grids if insufficent points on this processor
    call mgcollect(2,rhs(:,:,2),dsolmax_g,klim=klim)
  
  end do
  
  call END_LOG(mgfine_end)
  
  call START_LOG(mgup_begin)
  
  ! upscale grid
  do g=2,gmax
  
    ng=mg(g)%ifull
                
    ! update scalar field
    ! assume zero for first guess of residual (also avoids additional bounds call)
    !v(1:klim,1:ng,g)=0.
    do iq=1,ng
      v(1:klim,iq,g)=-rhs(1:klim,iq,g)/(helm(1:klim,iq,g)-mg(g)%zz(iq))
    end do
    call mgbounds(g,v(:,:,g),klim=klim)

    do i=1,itrbgn-1
      do iq=1,ng
        ! post smoothing
        w(1:klim,iq)=(mg(g)%zze(iq)*v(1:klim,mg(g)%ie(iq),g)+mg(g)%zzw(iq)*v(1:klim,mg(g)%iw(iq),g) &
                     +mg(g)%zzn(iq)*v(1:klim,mg(g)%in(iq),g)+mg(g)%zzs(iq)*v(1:klim,mg(g)%is(iq),g) &
                     -rhs(1:klim,iq,g))/(helm(1:klim,iq,g)-mg(g)%zz(iq))
      end do
      call mgbounds(g,w(:,:),klim=klim)
      v(1:klim,1:mg(g)%ifull+mg(g)%iextra,g)=w(1:klim,1:mg(g)%ifull+mg(g)%iextra)
    end do
    
    do iq=1,ng
      ! residual
      w(1:klim,iq)=-mg(g)%zze(iq)*v(1:klim,mg(g)%ie(iq),g)-mg(g)%zzw(iq)*v(1:klim,mg(g)%iw(iq),g) &
                   -mg(g)%zzn(iq)*v(1:klim,mg(g)%in(iq),g)-mg(g)%zzs(iq)*v(1:klim,mg(g)%is(iq),g) &
                   +rhs(1:klim,iq,g)+(helm(1:klim,iq,g)-mg(g)%zz(iq))*v(1:klim,iq,g)
    end do

    do iq=1,mg(g)%ifull_fine
      ! restriction
      ! (calculate coarser grid before mgcollect as more work is done in parallel)
      rhs(1:klim,iq,g+1)=0.25*(w(1:klim,mg(g)%fine(iq)  )+w(1:klim,mg(g)%fine_n(iq) )  &
                              +w(1:klim,mg(g)%fine_e(iq))+w(1:klim,mg(g)%fine_ne(iq)))
    end do

    ! merge grids if insufficent points on this processor
    call mgcollect(g+1,rhs(:,:,g+1),dsolmax_g,klim=klim)

  end do
  
  call END_LOG(mgup_end)
  
  call START_LOG(mgcoarse_begin)

  ! solve coarse grid
  do g=mg_maxlevel,mg_maxlevel_local ! same as if (mg_maxlevel_local==mg_maxlevel) then ...

    ng=mg(g)%ifull

    ! perform LU decomposition and back substitute with RHS
    ! to solve for v on coarse grid
    do k=1,klim
      v_o(1:ng)=rhs(k,1:ng,g)
      call mbacksub(helm_o(:,:,k),v_o(1:ng),indy(:,k))
      v(k,1:ng,g)=v_o(1:ng)
    end do
      
  end do

  call END_LOG(mgcoarse_end)
  
  call START_LOG(mgdown_begin)

  ! downscale grid
  do g=gmax,2,-1

    call mgbcast(g+1,v(:,:,g+1),dsolmax_g,klim=klim)

    do iq=1,mg(g+1)%ifull_coarse
      ! interpolation
      w(1:klim,iq)= mg(g+1)%wgt_a(iq)*v(1:klim,mg(g+1)%coarse_a(iq),g+1) + mg(g+1)%wgt_bc(iq)*v(1:klim,mg(g+1)%coarse_b(iq),g+1) &
                 + mg(g+1)%wgt_bc(iq)*v(1:klim,mg(g+1)%coarse_c(iq),g+1) +  mg(g+1)%wgt_d(iq)*v(1:klim,mg(g+1)%coarse_d(iq),g+1)

      ! extension
      ! No mgbounds as the v halo has already been updated and
      ! the coarse interpolation also updates the w halo
      w(1:klim,iq)=v(1:klim,iq,g)+w(1:klim,iq)
    end do

    ! MJT notes - The first correction is usually sufficient to produce a good approximation
    ! to the converged solution, thereby avoiding additional calls to mgbounds
    do i=1,itrend-1
      do iq=1,mg(g)%ifull
        ! post smoothing
        v(1:klim,iq,g)=(mg(g)%zze(iq)*w(1:klim,mg(g)%ie(iq))+mg(g)%zzw(iq)*w(1:klim,mg(g)%iw(iq)) &
                       +mg(g)%zzn(iq)*w(1:klim,mg(g)%in(iq))+mg(g)%zzs(iq)*w(1:klim,mg(g)%is(iq)) &
                       -rhs(1:klim,iq,g))/(helm(1:klim,iq,g)-mg(g)%zz(iq))
      end do
      call mgbounds(g,v(:,:,g),klim=klim)
      w(1:klim,1:mg(g)%ifull+mg(g)%iextra)=v(1:klim,1:mg(g)%ifull+mg(g)%iextra,g)
    end do
    
    do iq=1,mg(g)%ifull
      ! post smoothing
      v(1:klim,iq,g)=(mg(g)%zze(iq)*w(1:klim,mg(g)%ie(iq))+mg(g)%zzw(iq)*w(1:klim,mg(g)%iw(iq)) &
                     +mg(g)%zzn(iq)*w(1:klim,mg(g)%in(iq))+mg(g)%zzs(iq)*w(1:klim,mg(g)%is(iq)) &
                     -rhs(1:klim,iq,g))/(helm(1:klim,iq,g)-mg(g)%zz(iq))
    end do    
    call mgbounds(g,v(:,:,g),klim=klim,corner=.true.)

  end do
  
  do g=1,min(mg_maxlevel_local,1) ! same as if (mg_maxlevel_local>0) then ...
    
    ! fine grid
    call mgbcast(2,v(:,:,2),dsolmax_g,klim=klim)

    ! interpolation
    do iq=1,mg(2)%ifull_coarse
      w(1:klim,iq)= mg(2)%wgt_a(iq)*v(1:klim,mg(2)%coarse_a(iq),2) + mg(2)%wgt_bc(iq)*v(1:klim,mg(2)%coarse_b(iq),2) &
                 + mg(2)%wgt_bc(iq)*v(1:klim,mg(2)%coarse_c(iq),2) +  mg(2)%wgt_d(iq)*v(1:klim,mg(2)%coarse_d(iq),2)
    end do
    
  end do
  
  call END_LOG(mgdown_end)
  
  call START_LOG(mgfine_begin)

  ! multi-grid solver bounds indicies do not match standard iextra indicies, so we need to remap the halo
  if (mg(1)%merge_len>1) then
    call mgbcast(1,w,dsolmax_g,klim=klim)
    ir=mod(mg(1)%merge_pos-1,mg(1)%merge_row)+1   ! index for proc row
    ic=(mg(1)%merge_pos-1)/mg(1)%merge_row+1      ! index for proc col
    do n=1,npan
      do jj=1,jpan
        iq_a=1+(jj-1)*ipan+(n-1)*ipan*jpan
        iq_b=jj*ipan+(n-1)*ipan*jpan
        iq_c=1+(ir-1)*ipan+(jj-1+(ic-1)*jpan)*ipan*mg(1)%merge_row+(n-1)*ipan*jpan*mg(1)%merge_len
        iq_d=ir*ipan+(jj-1+(ic-1)*jpan)*ipan*mg(1)%merge_row+(n-1)*ipan*jpan*mg(1)%merge_len
        vdum(iq_a:iq_b,1:klim)=transpose(w(1:klim,iq_c:iq_d))
      end do
      do i=1,ipan
        iq_a=i+(n-1)*ipan*jpan
        iq_c=i+(ir-1)*ipan+((ic-1)*jpan)*ipan*mg(1)%merge_row+(n-1)*ipan*jpan*mg(1)%merge_len
        iq_b=is(iq_a)
        iq_d=mg(1)%is(iq_c)
        vdum(iq_b,1:klim)=w(1:klim,iq_d)
        iq_a=i+(jpan-1)*ipan+(n-1)*ipan*jpan
        iq_c=i+(ir-1)*ipan+(jpan-1+(ic-1)*jpan)*ipan*mg(1)%merge_row+(n-1)*ipan*jpan*mg(1)%merge_len
        iq_b=in(iq_a)
        iq_d=mg(1)%in(iq_c)
        vdum(iq_b,1:klim)=w(1:klim,iq_d)
      end do  
      do j=1,jpan
        iq_a=1+(j-1)*ipan+(n-1)*ipan*jpan
        iq_c=1+(ir-1)*ipan+(j-1+(ic-1)*jpan)*ipan*mg(1)%merge_row+(n-1)*ipan*jpan*mg(1)%merge_len
        iq_b=iw(iq_a)
        iq_d=mg(1)%iw(iq_c)
        vdum(iq_b,1:klim)=w(1:klim,iq_d)
        iq_a=ipan+(j-1)*ipan+(n-1)*ipan*jpan
        iq_c=ipan+(ir-1)*ipan+(j-1+(ic-1)*jpan)*ipan*mg(1)%merge_row+(n-1)*ipan*jpan*mg(1)%merge_len
        iq_b=ie(iq_a)
        iq_d=mg(1)%ie(iq_c)
        vdum(iq_b,1:klim)=w(1:klim,iq_d)
      end do
    end do
  else
    ! remap mg halo to normal halo
    vdum(1:ifull,1:klim)=transpose(w(1:klim,1:ifull))
    do n=0,npan-1
      do i=1,ipan
        iq=i+n*ipan*jpan
        iq_a=is(iq)
        iq_b=mg(1)%is(iq)
        vdum(iq_a,1:klim)=w(1:klim,iq_b)
        iq=i+(jpan-1)*ipan+n*ipan*jpan
        iq_a=in(iq)
        iq_b=mg(1)%in(iq)
        vdum(iq_a,1:klim)=w(1:klim,iq_b)
      end do  
      do j=1,jpan
        iq=1+(j-1)*ipan+n*ipan*jpan
        iq_a=iw(iq)
        iq_b=mg(1)%iw(iq)
        vdum(iq_a,1:klim)=w(1:klim,iq_b)
        iq=j*ipan+n*ipan*jpan
        iq_a=ie(iq)
        iq_b=mg(1)%ie(iq)
        vdum(iq_a,1:klim)=w(1:klim,iq_b)
      end do
    end do
  end if
  
  ! extension
  iv(1:ifull+iextra,1:klim)=iv(1:ifull+iextra,1:klim)+vdum(1:ifull+iextra,1:klim)
  
  do i=1,itrend
    ! post smoothing
    do nc=1,maxcolour
      do k=1,klim
        dsol(iqx(1:ifullcol(nc),nc),k)=( zznc(1:ifullcol(nc),nc)*iv(iqn(1:ifullcol(nc),nc),k)    &
                                       + zzwc(1:ifullcol(nc),nc)*iv(iqw(1:ifullcol(nc),nc),k)    &
                                       + zzec(1:ifullcol(nc),nc)*iv(iqe(1:ifullcol(nc),nc),k)    &
                                       + zzsc(1:ifullcol(nc),nc)*iv(iqs(1:ifullcol(nc),nc),k)    &
                                       - rhsc(1:ifullcol(nc),k,nc))*rhelmc(1:ifullcol(nc),k,nc)  &
                                       - iv(iqx(1:ifullcol(nc),nc),k)
        iv(iqx(1:ifullcol(nc),nc),k) = iv(iqx(1:ifullcol(nc),nc),k) + dsol(iqx(1:ifullcol(nc),nc),k)
      end do
      call bounds_colour(iv,nc,klim=klim)
    end do
  end do

  call END_LOG(mgfine_end)

  ! test for convergence
  knew=klim
  !if (dsolmax_g(1)<10.*restol*sdif(1)) then
  !  dsolmax_l(1:klim)=maxval(abs(dsol(1:ifull,1:klim)),dim=1)
  !  call ccmpi_allreduce(dsolmax_l(1:klim),dsolmax_g(1:klim),"max",comm_world)
  !end if
  do k=klim,1,-1
    iters(k)=itr
    if (dsolmax_g(k)>=restol*sdif(k)) exit
    knew=k-1
  end do
  klim=knew
  if (klim<1) exit
 
end do

! JLM suggestion
do k=1,kl
  iv(1:ifull,k)=iv(1:ifull,k)+savg(k)
end do

if (myid==0) then
  if (ktau<6.or.iters(1)>itr_mg) then
    do k=1,kl
      write(6,*) "mg ktau,k,iter ",ktau,k,iters(k),dsolmax_g(k)
    end do
  end if
end if

call END_LOG(helm_end)

return
end subroutine mghelm

! This version is for the ocean and ice
subroutine mgmlo(neta,ipice,iyy,iyyn,iyys,iyye,iyyw,izz,izzn,izzs,izze,izzw,ihh,irhs,tol,itol,totits,maxglobseta,maxglobip, &
                 ipmax,ee,dd)

use cc_mpi
use indices_m

implicit none

include 'newmpar.h'

integer, intent(out) :: totits
integer, dimension(mg_minsize) :: indy
integer itr,itrc,g,ng,ng4,n,i,j,ir,ic,jj,iq,k
integer iq_a,iq_b,iq_c,iq_d
integer nc,gb
real, intent(in) :: tol,itol
real, intent(out) :: maxglobseta,maxglobip
real, dimension(ifull+iextra), intent(inout) :: neta,ipice
real, dimension(ifull+iextra), intent(in) :: ee,dd
real, dimension(ifull+iextra), intent(in) :: ipmax
real, dimension(ifull), intent(in) :: iyy,iyyn,iyys,iyye,iyyw
real, dimension(ifull,2), intent(in) :: izz,izzn,izzs,izze,izzw
real, dimension(ifull), intent(in) :: ihh
real, dimension(ifull,2), intent(in) :: irhs
real, dimension(ifull+iextra) :: vduma,vdumb
real, dimension(ifullmaxcol,maxcolour) :: rhsc,rhscice,ddc,eec,ipmaxc
real, dimension(ifullmaxcol,maxcolour) :: yyc,yync,yysc,yyec,yywc
real, dimension(ifullmaxcol,maxcolour) :: zzhhc,zznc,zzsc,zzec,zzwc
real, dimension(ifullmaxcol,maxcolour) :: zzcice,zzncice,zzscice,zzecice,zzwcice
real, dimension(mg_maxsize) :: bu,cu
real, dimension(mg_maxsize,2,gmax+1) :: v
real, dimension(mg_maxsize,2:gmax+1) :: yyn,yys,yye,yyw,yyz
real, dimension(mg_maxsize,2:gmax+1) :: zznice,zzsice,zzeice,zzwice,zzzice
real, dimension(mg_maxsize,18) :: w
real, dimension(mg_maxsize,gmax+1) :: zz,zzn,zzs,zze,zzw
real, dimension(mg_maxsize,gmax+1) :: hh
real, dimension(mg_maxsize,gmax+1) :: rhs
real, dimension(mg_maxsize,gmax+1) :: rhsice
real, dimension(mg_maxsize,2) :: dsol
real, dimension(mg_maxsize) :: ws,new
real, dimension(ifull+iextra,2) :: dumc
real, dimension(mg_maxsize,2) :: dumc_n,dumc_s,dumc_e,dumc_w
real, dimension(mg_minsize,mg_minsize) :: helm_o
real, dimension(mg_ifullmaxcol,3) :: yyzcu,yyncu,yyscu,yyecu,yywcu
real, dimension(mg_ifullmaxcol,3) :: zzhhcu,zzncu,zzscu,zzecu,zzwcu,rhscu
real, dimension(2) :: dsolmax
real, dimension(8) :: dsolmax_g

if (sorfirst) then
  write(6,*) "ERROR: mgsormlo requires mgsor_init to be called first"
  call ccmpi_abort(-1)
end if

! The following expressions describe the residual terms of the ocean model

! yy*neta*(DIV^2 neta) + zz*(DIV^2 neta) + hh*neta = rhs
! neta = n0 + e
! neta is the solution, n0 is first guess and e is the residual term

! yy*n0*(DIV^2 n0) + zz*(DIV^2 n0) + hh*n0 = rhs - E
! yy*e*(DIV^2 e) + (yy*n0+zz)*(DIV^2 e) + (yy*(DIV^2 n0)+hh)*e  = E

! so yy is simply upscaled
! zz -> zz + yy*n0
! hh -> hh + yy*(DIV^2 n0)

! also for sea-ice

! zz*(DIV^2 ipice) = rhs
! ipice = i0 + f

! zz*(DIV^2 i0) = rhs - F
! zz*(DIV^2 f) = F

call START_LOG(mgmlosetup_begin)

vduma=0.
vdumb=0.
dumc=0.

! pack colour arrays
do nc=1,maxcolour
  yyc(1:ifullcol(nc),nc)    =iyy(iqx(1:ifullcol(nc),nc))
  yync(1:ifullcol(nc),nc)   =iyyn(iqx(1:ifullcol(nc),nc))
  yysc(1:ifullcol(nc),nc)   =iyys(iqx(1:ifullcol(nc),nc))
  yyec(1:ifullcol(nc),nc)   =iyye(iqx(1:ifullcol(nc),nc))
  yywc(1:ifullcol(nc),nc)   =iyyw(iqx(1:ifullcol(nc),nc))
  zzhhc(1:ifullcol(nc),nc)  =izz(iqx(1:ifullcol(nc),nc),1)+ihh(iqx(1:ifullcol(nc),nc))
  zznc(1:ifullcol(nc),nc)   =izzn(iqx(1:ifullcol(nc),nc),1)
  zzsc(1:ifullcol(nc),nc)   =izzs(iqx(1:ifullcol(nc),nc),1)
  zzec(1:ifullcol(nc),nc)   =izze(iqx(1:ifullcol(nc),nc),1)
  zzwc(1:ifullcol(nc),nc)   =izzw(iqx(1:ifullcol(nc),nc),1)
  zzcice(1:ifullcol(nc),nc) =izz(iqx(1:ifullcol(nc),nc),2)
  zzncice(1:ifullcol(nc),nc)=izzn(iqx(1:ifullcol(nc),nc),2)
  zzscice(1:ifullcol(nc),nc)=izzs(iqx(1:ifullcol(nc),nc),2)
  zzecice(1:ifullcol(nc),nc)=izze(iqx(1:ifullcol(nc),nc),2)
  zzwcice(1:ifullcol(nc),nc)=izzw(iqx(1:ifullcol(nc),nc),2)
  rhsc(1:ifullcol(nc),nc)   =irhs(iqx(1:ifullcol(nc),nc),1)
  rhscice(1:ifullcol(nc),nc)=irhs(iqx(1:ifullcol(nc),nc),2)
  ddc(1:ifullcol(nc),nc)    =dd(iqx(1:ifullcol(nc),nc))
  eec(1:ifullcol(nc),nc)    =ee(iqx(1:ifullcol(nc),nc))
  ipmaxc(1:ifullcol(nc),nc) =ipmax(iqx(1:ifullcol(nc),nc))
end do

! solver requires bounds to be updated
dumc(1:ifull,1)=neta(1:ifull)
dumc(1:ifull,2)=ipice(1:ifull)
call bounds(dumc)

do nc=1,maxcolour
  
  dumc_n(1:ifullcol(nc),:)=dumc(iqn(1:ifullcol(nc),nc),:)
  dumc_s(1:ifullcol(nc),:)=dumc(iqs(1:ifullcol(nc),nc),:)
  dumc_e(1:ifullcol(nc),:)=dumc(iqe(1:ifullcol(nc),nc),:)
  dumc_w(1:ifullcol(nc),:)=dumc(iqw(1:ifullcol(nc),nc),:)
  
  ! ocean
  bu(1:ifullcol(nc))=zzhhc(1:ifullcol(nc),nc)                                                                    &
              +yync(1:ifullcol(nc),nc)*dumc_n(1:ifullcol(nc),1)+yysc(1:ifullcol(nc),nc)*dumc_s(1:ifullcol(nc),1) &
              +yyec(1:ifullcol(nc),nc)*dumc_e(1:ifullcol(nc),1)+yywc(1:ifullcol(nc),nc)*dumc_w(1:ifullcol(nc),1)
  cu(1:ifullcol(nc))=zznc(1:ifullcol(nc),nc)*dumc_n(1:ifullcol(nc),1)+zzsc(1:ifullcol(nc),nc)*dumc_s(1:ifullcol(nc),1) &
                    +zzec(1:ifullcol(nc),nc)*dumc_e(1:ifullcol(nc),1)+zzwc(1:ifullcol(nc),nc)*dumc_w(1:ifullcol(nc),1) &
                    -rhsc(1:ifullcol(nc),nc)        
  dumc(iqx(1:ifullcol(nc),nc),1)=eec(1:ifullcol(nc),nc)*max(-ddc(1:ifullcol(nc),nc),                                         &
     -2.*cu(1:ifullcol(nc))/(bu(1:ifullcol(nc))+sqrt(bu(1:ifullcol(nc))**2-4.*yyc(1:ifullcol(nc),nc)*cu(1:ifullcol(nc)))) )
    
  ! ice (cavitating fluid)
  dumc(iqx(1:ifullcol(nc),nc),2) = max(0.,min(ipmaxc(1:ifullcol(nc),nc), &
     ( -zzncice(1:ifullcol(nc),nc)*dumc_n(1:ifullcol(nc),2) &
       -zzscice(1:ifullcol(nc),nc)*dumc_s(1:ifullcol(nc),2) &
       -zzecice(1:ifullcol(nc),nc)*dumc_e(1:ifullcol(nc),2) &
       -zzwcice(1:ifullcol(nc),nc)*dumc_w(1:ifullcol(nc),2) &
      + rhscice(1:ifullcol(nc),nc) ) / zzcice(1:ifullcol(nc),nc) ))

  call bounds_colour(dumc,nc)

end do
neta(1:ifull+iextra) =dumc(1:ifull+iextra,1)
ipice(1:ifull+iextra)=dumc(1:ifull+iextra,2)  

dumc_n(1:ifull,:)=dumc(in,:)
dumc_s(1:ifull,:)=dumc(is,:)
dumc_e(1:ifull,:)=dumc(ie,:)
dumc_w(1:ifull,:)=dumc(iw,:)

w(1:ifull,2)= izz(:,1)+ iyy*neta(1:ifull)
w(1:ifull,3)=izzn(:,1)+iyyn*neta(1:ifull)
w(1:ifull,4)=izzs(:,1)+iyys*neta(1:ifull)
w(1:ifull,5)=izze(:,1)+iyye*neta(1:ifull)
w(1:ifull,6)=izzw(:,1)+iyyw*neta(1:ifull)
w(1:ifull,7)=ihh+iyy*neta(1:ifull)+iyyn*dumc_n(1:ifull,1)+iyys*dumc_s(1:ifull,1)+iyye*dumc_e(1:ifull,1)+iyyw*dumc_w(1:ifull,1)

! residual
w(1:ifull,1)=(-neta(1:ifull)*(     iyy*neta(1:ifull)     +iyyn*dumc_n(1:ifull,1)     +iyys*dumc_s(1:ifull,1)   &
                                                         +iyye*dumc_e(1:ifull,1)     +iyyw*dumc_w(1:ifull,1))  &
                            -(izz(:,1)*neta(1:ifull)+izzn(:,1)*dumc_n(1:ifull,1)+izzs(:,1)*dumc_s(1:ifull,1)   &
                                                    +izze(:,1)*dumc_e(1:ifull,1)+izzw(:,1)*dumc_w(1:ifull,1))  &
                            -ihh*neta(1:ifull)+irhs(:,1))*ee(1:ifull)
where (ipice(1:ifull)>=ipmax(1:ifull))
  w(1:ifull,8)=0. ! patch to remove error when ipmax is reached - improves convergence
elsewhere
  w(1:ifull,8)=(-(izz(:,2)*ipice(1:ifull)+izzn(:,2)*dumc_n(1:ifull,2)+izzs(:,2)*dumc_s(1:ifull,2) &
                +izze(:,2)*dumc_e(1:ifull,2)+izzw(:,2)*dumc_w(1:ifull,2))+irhs(:,2))*ee(1:ifull)
end where

! upscale coeffs
w(1:ifull,9)=iyy(1:ifull)
w(1:ifull,10)=iyyn(1:ifull)
w(1:ifull,11)=iyys(1:ifull)
w(1:ifull,12)=iyye(1:ifull)
w(1:ifull,13)=iyyw(1:ifull)
w(1:ifull,14)=izz(1:ifull,2)
w(1:ifull,15)=izzn(1:ifull,2)
w(1:ifull,16)=izzs(1:ifull,2)
w(1:ifull,17)=izze(1:ifull,2)
w(1:ifull,18)=izzw(1:ifull,2)
call mgcollect_mlo(1,w(:,1:18))
  
do g=1,min(mg_maxlevel_local,1) ! same as if (mg_maxlevel_local>0) then ...
  
  ! restriction
  ! (since this always operates within a panel, then ine = ien is always true)
  ng4=mg(1)%ifull_fine
  rhs(1:ng4,2)=0.25*(w(mg(1)%fine  ,1)+w(mg(1)%fine_n ,1) &
                    +w(mg(1)%fine_e,1)+w(mg(1)%fine_ne,1))
  zz(1:ng4,2) =0.25*dfac*(w(mg(1)%fine  ,2)+w(mg(1)%fine_n ,2) &
                         +w(mg(1)%fine_e,2)+w(mg(1)%fine_ne,2))
  zzn(1:ng4,2)=0.25*dfac*(w(mg(1)%fine  ,3)+w(mg(1)%fine_n ,3) &
                         +w(mg(1)%fine_e,3)+w(mg(1)%fine_ne,3))
  zzs(1:ng4,2)=0.25*dfac*(w(mg(1)%fine  ,4)+w(mg(1)%fine_n ,4) &
                         +w(mg(1)%fine_e,4)+w(mg(1)%fine_ne,4))
  zze(1:ng4,2)=0.25*dfac*(w(mg(1)%fine  ,5)+w(mg(1)%fine_n ,5) &
                         +w(mg(1)%fine_e,5)+w(mg(1)%fine_ne,5))
  zzw(1:ng4,2)=0.25*dfac*(w(mg(1)%fine  ,6)+w(mg(1)%fine_n ,6) &
                         +w(mg(1)%fine_e,6)+w(mg(1)%fine_ne,6))
  hh(1:ng4,2)    =0.25*(w(mg(1)%fine  ,7)+w(mg(1)%fine_n ,7) &
                       +w(mg(1)%fine_e,7)+w(mg(1)%fine_ne,7))
  rhsice(1:ng4,2)=0.25*(w(mg(1)%fine  ,8)+w(mg(1)%fine_n ,8) &
                       +w(mg(1)%fine_e,8)+w(mg(1)%fine_ne,8))

  yyz(1:ng4,2)=0.25*dfac*(w(mg(1)%fine  ,9)+w(mg(1)%fine_n ,9) &
                         +w(mg(1)%fine_e,9)+w(mg(1)%fine_ne,9))
  yyn(1:ng4,2)=0.25*dfac*(w(mg(1)%fine  ,10)+w(mg(1)%fine_n ,10) &
                         +w(mg(1)%fine_e,10)+w(mg(1)%fine_ne,10))
  yys(1:ng4,2)=0.25*dfac*(w(mg(1)%fine  ,11)+w(mg(1)%fine_n ,11) &
                         +w(mg(1)%fine_e,11)+w(mg(1)%fine_ne,11))
  yye(1:ng4,2)=0.25*dfac*(w(mg(1)%fine  ,12)+w(mg(1)%fine_n ,12) &
                         +w(mg(1)%fine_e,12)+w(mg(1)%fine_ne,12))
  yyw(1:ng4,2)=0.25*dfac*(w(mg(1)%fine  ,13)+w(mg(1)%fine_n ,13) &
                         +w(mg(1)%fine_e,13)+w(mg(1)%fine_ne,13))
  ! special treatment of cavitating fluid (no dfac)
  zzzice(1:ng4,2)=0.25*(w(mg(1)%fine  ,14)+w(mg(1)%fine_n ,14)   &
                       +w(mg(1)%fine_e,14)+w(mg(1)%fine_ne,14))
  zznice(1:ng4,2)=0.25*(w(mg(1)%fine  ,15)+w(mg(1)%fine_n ,15)   &
                       +w(mg(1)%fine_e,15)+w(mg(1)%fine_ne,15))
  zzsice(1:ng4,2)=0.25*(w(mg(1)%fine  ,16)+w(mg(1)%fine_n ,16)   &
                       +w(mg(1)%fine_e,16)+w(mg(1)%fine_ne,16))
  zzeice(1:ng4,2)=0.25*(w(mg(1)%fine  ,17)+w(mg(1)%fine_n ,17)   &
                       +w(mg(1)%fine_e,17)+w(mg(1)%fine_ne,17))
  zzwice(1:ng4,2)=0.25*(w(mg(1)%fine  ,18)+w(mg(1)%fine_n ,18)   &
                       +w(mg(1)%fine_e,18)+w(mg(1)%fine_ne,18))

  ! merge grids if insufficent points on this processor
  if (mg(2)%merge_len>1) then
    w(1:ng4,1)  =rhs(1:ng4,2)
    w(1:ng4,2)  =zz(1:ng4,2)
    w(1:ng4,3)  =zzn(1:ng4,2)
    w(1:ng4,4)  =zzs(1:ng4,2)
    w(1:ng4,5)  =zze(1:ng4,2)
    w(1:ng4,6)  =zzw(1:ng4,2)
    w(1:ng4,7)  =hh(1:ng4,2)
    w(1:ng4,8)  =rhsice(1:ng4,2)
    w(1:ng4,9)  =yyz(1:ng4,2)
    w(1:ng4,10) =yyn(1:ng4,2)
    w(1:ng4,11) =yys(1:ng4,2)
    w(1:ng4,12) =yye(1:ng4,2)
    w(1:ng4,13) =yyw(1:ng4,2)
    w(1:ng4,14) =zzzice(1:ng4,2)
    w(1:ng4,15) =zznice(1:ng4,2)
    w(1:ng4,16) =zzsice(1:ng4,2)
    w(1:ng4,17) =zzeice(1:ng4,2)
    w(1:ng4,18) =zzwice(1:ng4,2)
    call mgcollect_mlo(2,w(:,1:18))
    ng=mg(2)%ifull
    rhs(1:ng,2)    =w(1:ng,1)
    zz(1:ng,2)     =w(1:ng,2)
    zzn(1:ng,2)    =w(1:ng,3)
    zzs(1:ng,2)    =w(1:ng,4)
    zze(1:ng,2)    =w(1:ng,5)
    zzw(1:ng,2)    =w(1:ng,6)
    hh(1:ng,2)     =w(1:ng,7)
    rhsice(1:ng,2) =w(1:ng,8)
    yyz(1:ng,2)    =w(1:ng,9)
    yyn(1:ng,2)    =w(1:ng,10)
    yys(1:ng,2)    =w(1:ng,11)
    yye(1:ng,2)    =w(1:ng,12)
    yyw(1:ng,2)    =w(1:ng,13)
    zzzice(1:ng,2) =w(1:ng,14)
    zznice(1:ng,2) =w(1:ng,15)
    zzsice(1:ng,2) =w(1:ng,16)
    zzeice(1:ng,2) =w(1:ng,17)
    zzwice(1:ng,2) =w(1:ng,18)
  end if
    
end do
  
! upscale grid
do g=2,gmax
  
  ng=mg(g)%ifull

  ! update
  ! possibly use colours here, although v is reset to zero every iteration
  ! assume zero for first guess of residual (also avoids additional bounds call)
  bu(1:ng)=zz(1:ng,g)+hh(1:ng,g)
  v(1:ng,1,g) = 2.*rhs(1:ng,g)/(bu(1:ng)+sqrt(bu(1:ng)*bu(1:ng)+4.*yyz(1:ng,g)*rhs(1:ng,g)))
  v(1:ng,2,g) = rhsice(1:ng,g)/zzzice(1:ng,g)
  call mgbounds_mlo(g,v(:,:,g))
    
  do i=1,itrbgn-1
    dumc_n(1:ng,:)=v(mg(g)%in,1:2,g)
    dumc_s(1:ng,:)=v(mg(g)%is,1:2,g)
    dumc_e(1:ng,:)=v(mg(g)%ie,1:2,g)
    dumc_w(1:ng,:)=v(mg(g)%iw,1:2,g)

    ! ocean
    ! post smoothing
    bu(1:ng)=zz(1:ng,g)+hh(1:ng,g)+yyn(1:ng,g)*dumc_n(1:ng,1)+yys(1:ng,g)*dumc_s(1:ng,1) &
                                  +yye(1:ng,g)*dumc_e(1:ng,1)+yyw(1:ng,g)*dumc_w(1:ng,1)
    cu(1:ng)=zzn(1:ng,g)*dumc_n(1:ng,1)+zzs(1:ng,g)*dumc_s(1:ng,1)             &
            +zze(1:ng,g)*dumc_e(1:ng,1)+zzw(1:ng,g)*dumc_w(1:ng,1)-rhs(1:ng,g)
    w(1:ng,1) = -2.*cu(1:ng)/(bu(1:ng)+sqrt(bu(1:ng)*bu(1:ng)-4.*yyz(1:ng,g)*cu(1:ng)))

    ! ice
    w(1:ng,2) = ( -zznice(1:ng,g)*dumc_n(1:ng,2)-zzsice(1:ng,g)*dumc_s(1:ng,2) &
                  -zzeice(1:ng,g)*dumc_e(1:ng,2)-zzwice(1:ng,g)*dumc_w(1:ng,2) &
                  +rhsice(1:ng,g) ) / zzzice(1:ng,g)

    call mgbounds_mlo(g,w(:,:))
    v(1:mg(g)%ifull+mg(g)%iextra,1:2,g)=w(1:mg(g)%ifull+mg(g)%iextra,1:2)
  end do
  
  ! restriction
  ! (calculate finer grid before mgcollect as the messages sent/recv are shorter)
  dumc_n(1:ng,:)=v(mg(g)%in,:,g)
  dumc_s(1:ng,:)=v(mg(g)%is,:,g)
  dumc_e(1:ng,:)=v(mg(g)%ie,:,g)
  dumc_w(1:ng,:)=v(mg(g)%iw,:,g)

  ng4=mg(g)%ifull_fine
  ws(1:ng)= zz(1:ng,g)+yyz(1:ng,g)*v(1:ng,1,g)
  zz(1:ng4,g+1)=0.25*dfac*(ws(mg(g)%fine  )+ws(mg(g)%fine_n ) &
                          +ws(mg(g)%fine_e)+ws(mg(g)%fine_ne))
  ws(1:ng)=zzn(1:ng,g)+yyn(1:ng,g)*v(1:ng,1,g)
  zzn(1:ng4,g+1)=0.25*dfac*(ws(mg(g)%fine  )+ws(mg(g)%fine_n ) &
                           +ws(mg(g)%fine_e)+ws(mg(g)%fine_ne))
  ws(1:ng)=zzs(1:ng,g)+yys(1:ng,g)*v(1:ng,1,g)
  zzs(1:ng4,g+1)=0.25*dfac*(ws(mg(g)%fine  )+ws(mg(g)%fine_n ) &
                           +ws(mg(g)%fine_e)+ws(mg(g)%fine_ne))
  ws(1:ng)=zze(1:ng,g)+yye(1:ng,g)*v(1:ng,1,g)
  zze(1:ng4,g+1)=0.25*dfac*(ws(mg(g)%fine  )+ws(mg(g)%fine_n ) &
                           +ws(mg(g)%fine_e)+ws(mg(g)%fine_ne))
  ws(1:ng)=zzw(1:ng,g)+yyw(1:ng,g)*v(1:ng,1,g)
  zzw(1:ng4,g+1)=0.25*dfac*(ws(mg(g)%fine  )+ws(mg(g)%fine_n ) &
                           +ws(mg(g)%fine_e)+ws(mg(g)%fine_ne))
  ws(1:ng)=hh(1:ng,g)+yyz(1:ng,g)*v(1:ng,1,g)+yyn(1:ng,g)*dumc_n(1:ng,1)+yys(1:ng,g)*dumc_s(1:ng,1) &
                                             +yye(1:ng,g)*dumc_e(1:ng,1)+yyw(1:ng,g)*dumc_w(1:ng,1)
  hh(1:ng4,g+1)=0.25*(ws(mg(g)%fine  )+ws(mg(g)%fine_n ) &
                     +ws(mg(g)%fine_e)+ws(mg(g)%fine_ne))
   ! ocean
  ws(1:ng)=-v(1:ng,1,g)*(yyz(1:ng,g)*v(1:ng,1,g)+yyn(1:ng,g)*dumc_n(1:ng,1)+yys(1:ng,g)*dumc_s(1:ng,1)  &
                                                +yye(1:ng,g)*dumc_e(1:ng,1)+yyw(1:ng,g)*dumc_w(1:ng,1)) &
                        -(zz(1:ng,g)*v(1:ng,1,g)+zzn(1:ng,g)*dumc_n(1:ng,1)+zzs(1:ng,g)*dumc_s(1:ng,1)  &
                                                +zze(1:ng,g)*dumc_e(1:ng,1)+zzw(1:ng,g)*dumc_w(1:ng,1)) &
                        -hh(1:ng,g)*v(1:ng,1,g)+rhs(1:ng,g)
  rhs(1:ng4,g+1)=0.25*(ws(mg(g)%fine  )+ws(mg(g)%fine_n ) &
                      +ws(mg(g)%fine_e)+ws(mg(g)%fine_ne))
                          
  ! ice
  ws(1:ng)=-(zznice(1:ng,g)*dumc_n(1:ng,2)+zzsice(1:ng,g)*dumc_s(1:ng,2)   &
            +zzeice(1:ng,g)*dumc_e(1:ng,2)+zzwice(1:ng,g)*dumc_w(1:ng,2))  &
            -zzzice(1:ng,g)*v(1:ng,2,g)+rhsice(1:ng,g)
  rhsice(1:ng4,g+1)=0.25*(ws(mg(g)%fine  )+ws(mg(g)%fine_n ) &
                         +ws(mg(g)%fine_e)+ws(mg(g)%fine_ne))

  yyz(1:ng4,g+1)=0.25*dfac*(w(mg(g)%fine  ,9)+w(mg(g)%fine_n ,9)   &
                           +w(mg(g)%fine_e,9)+w(mg(g)%fine_ne,9))
  yyn(1:ng4,g+1)=0.25*dfac*(w(mg(g)%fine  ,10)+w(mg(g)%fine_n ,10) &
                           +w(mg(g)%fine_e,10)+w(mg(g)%fine_ne,10))
  yys(1:ng4,g+1)=0.25*dfac*(w(mg(g)%fine  ,11)+w(mg(g)%fine_n ,11) &
                           +w(mg(g)%fine_e,11)+w(mg(g)%fine_ne,11))
  yye(1:ng4,g+1)=0.25*dfac*(w(mg(g)%fine  ,12)+w(mg(g)%fine_n ,12) &
                           +w(mg(g)%fine_e,12)+w(mg(g)%fine_ne,12))
  yyw(1:ng4,g+1)=0.25*dfac*(w(mg(g)%fine  ,13)+w(mg(g)%fine_n ,13) &
                           +w(mg(g)%fine_e,13)+w(mg(g)%fine_ne,13))
  ! special treatment of cavitating fluid (no dfac)
  zzzice(1:ng4,g+1)=0.25*(w(mg(g)%fine  ,14)+w(mg(g)%fine_n ,14)   &
                         +w(mg(g)%fine_e,14)+w(mg(g)%fine_ne,14))
  zznice(1:ng4,g+1)=0.25*(w(mg(g)%fine  ,15)+w(mg(g)%fine_n ,15)   &
                         +w(mg(g)%fine_e,15)+w(mg(g)%fine_ne,15))
  zzsice(1:ng4,g+1)=0.25*(w(mg(g)%fine  ,16)+w(mg(g)%fine_n ,16)   &
                         +w(mg(g)%fine_e,16)+w(mg(g)%fine_ne,16))
  zzeice(1:ng4,g+1)=0.25*(w(mg(g)%fine  ,17)+w(mg(g)%fine_n ,17)   &
                         +w(mg(g)%fine_e,17)+w(mg(g)%fine_ne,17))
  zzwice(1:ng4,g+1)=0.25*(w(mg(g)%fine  ,18)+w(mg(g)%fine_n ,18)   &
                         +w(mg(g)%fine_e,18)+w(mg(g)%fine_ne,18))

  ! merge grids if insufficent points on this processor
  if (mg(g+1)%merge_len>1) then
    w(1:ng4,1)  =rhs(1:ng4,g+1)
    w(1:ng4,2)  =zz(1:ng4,g+1)
    w(1:ng4,3)  =zzn(1:ng4,g+1)
    w(1:ng4,4)  =zzs(1:ng4,g+1)
    w(1:ng4,5)  =zze(1:ng4,g+1)
    w(1:ng4,6)  =zzw(1:ng4,g+1)
    w(1:ng4,7)  =hh(1:ng4,g+1)
    w(1:ng4,8)  =rhsice(1:ng4,g+1)
    w(1:ng4,9)  =yyz(1:ng4,g+1)
    w(1:ng4,10) =yyn(1:ng4,g+1)
    w(1:ng4,11) =yys(1:ng4,g+1)
    w(1:ng4,12) =yye(1:ng4,g+1)
    w(1:ng4,13) =yyw(1:ng4,g+1)
    w(1:ng4,14) =zzzice(1:ng4,g+1)
    w(1:ng4,15) =zznice(1:ng4,g+1)
    w(1:ng4,16) =zzsice(1:ng4,g+1)
    w(1:ng4,17) =zzeice(1:ng4,g+1)
    w(1:ng4,18) =zzwice(1:ng4,g+1)
    call mgcollect_mlo(g+1,w(:,1:18))
    ng=mg(g+1)%ifull
    rhs(1:ng,g+1)    =w(1:ng,1)
    zz(1:ng,g+1)     =w(1:ng,2)
    zzn(1:ng,g+1)    =w(1:ng,3)
    zzs(1:ng,g+1)    =w(1:ng,4)
    zze(1:ng,g+1)    =w(1:ng,5)
    zzw(1:ng,g+1)    =w(1:ng,6)
    hh(1:ng,g+1)     =w(1:ng,7)
    rhsice(1:ng,g+1) =w(1:ng,8)
    yyz(1:ng,g+1)    =w(1:ng,9)
    yyn(1:ng,g+1)    =w(1:ng,10)
    yys(1:ng,g+1)    =w(1:ng,11)
    yye(1:ng,g+1)    =w(1:ng,12)
    yyw(1:ng,g+1)    =w(1:ng,13)
    zzzice(1:ng,g+1) =w(1:ng,14)
    zznice(1:ng,g+1) =w(1:ng,15)
    zzsice(1:ng,g+1) =w(1:ng,16)
    zzeice(1:ng,g+1) =w(1:ng,17)
    zzwice(1:ng,g+1) =w(1:ng,18)    
  end if

end do

do g=mg_maxlevel,mg_maxlevel_local ! same as if (mg_maxlevel_local==mg_maxlevel) then ...
  helm_o=0.
  ng=mg(g)%ifull
  do iq=1,ng
    helm_o(iq,iq)=zzzice(iq,g)
    helm_o(mg(g)%in(iq),iq)=zznice(iq,g)
    helm_o(mg(g)%is(iq),iq)=zzsice(iq,g)
    helm_o(mg(g)%ie(iq),iq)=zzeice(iq,g)
    helm_o(mg(g)%iw(iq),iq)=zzwice(iq,g)
  end do
  call mdecomp(helm_o,indy) ! destroys helm_m
  ! pack yy by colour
  do nc=1,3
    yyzcu(1:mg_ifullmaxcol,nc)=yyz(col_iq(:,nc),g)
    yyncu(1:mg_ifullmaxcol,nc)=yyn(col_iq(:,nc),g)
    yyscu(1:mg_ifullmaxcol,nc)=yys(col_iq(:,nc),g)
    yyecu(1:mg_ifullmaxcol,nc)=yye(col_iq(:,nc),g)
    yywcu(1:mg_ifullmaxcol,nc)=yyw(col_iq(:,nc),g)
  end do
      
  ! solve for ice using LU decomposition and back substitution with RHS
  v(1:ng,2,g)=rhsice(1:ng,g)
  call mbacksub(helm_o,v(1:ng,2,g),indy)

  ! solve non-linear water free surface with coloured SOR
    
  ! first guess
  bu(1:ng)=zz(1:ng,g)+hh(1:ng,g)
  v(1:ng,1,g) = 2.*rhs(1:ng,g)/(bu(1:ng)+sqrt(bu(1:ng)**2+4.*yyz(1:ng,g)*rhs(1:ng,g)))
  
  ! pack zz,hh and rhs by colour
  do nc=1,3
    zzhhcu(:,nc)=zz(col_iq(:,nc),g)+hh(col_iq(:,nc),g)
    zzncu(:,nc)=zzn(col_iq(:,nc),g)
    zzscu(:,nc)=zzs(col_iq(:,nc),g)
    zzecu(:,nc)=zze(col_iq(:,nc),g)
    zzwcu(:,nc)=zzw(col_iq(:,nc),g)
    rhscu(:,nc)=rhs(col_iq(:,nc),g)
  end do
  
  do itrc=1,itr_mgice

    ! store previous guess for convegence test
    ws(1:ng)=v(1:ng,1,g)
 
    do nc=1,3
      
      dumc_n(1:mg_ifullmaxcol,1)=v(col_iqn(:,nc),1,g)
      dumc_s(1:mg_ifullmaxcol,1)=v(col_iqs(:,nc),1,g)
      dumc_e(1:mg_ifullmaxcol,1)=v(col_iqe(:,nc),1,g)
      dumc_w(1:mg_ifullmaxcol,1)=v(col_iqw(:,nc),1,g)
      
      bu(1:mg_ifullmaxcol)=zzhhcu(:,nc)+yyncu(:,nc)*dumc_n(1:mg_ifullmaxcol,1)+yyscu(:,nc)*dumc_s(1:mg_ifullmaxcol,1) &
                                       +yyecu(:,nc)*dumc_e(1:mg_ifullmaxcol,1)+yywcu(:,nc)*dumc_w(1:mg_ifullmaxcol,1)
      cu(1:mg_ifullmaxcol)=zzncu(:,nc)*dumc_n(1:mg_ifullmaxcol,1)+zzscu(:,nc)*dumc_s(1:mg_ifullmaxcol,1)   &
                          +zzecu(:,nc)*dumc_e(1:mg_ifullmaxcol,1)+zzwcu(:,nc)*dumc_w(1:mg_ifullmaxcol,1)   &
                          -rhscu(:,nc)
      v(col_iq(:,nc),1,g) = -2.*cu(1:mg_ifullmaxcol)/(bu(1:mg_ifullmaxcol)                      &
                            +sqrt(bu(1:mg_ifullmaxcol)**2-4.*yyzcu(:,nc)*cu(1:mg_ifullmaxcol)))
    end do
      
    ! test for convergence
    dsol(1:ng,1)=v(1:ng,1,g)-ws(1:ng)
    dsolmax(1)=maxval(abs(dsol(1:ng,1)))
    if (dsolmax(1)<tol) exit

  end do
  
end do
  
! downscale grid
do g=gmax,2,-1

  call mgbcasta_mlo(g+1,v(:,:,g+1))

  ! interpolation
  ng4=mg(g+1)%ifull_coarse
    
  dumc_n(1:ng4,1:2)=v(mg(g+1)%coarse_a,1:2,g+1)
  dumc_s(1:ng4,1:2)=v(mg(g+1)%coarse_b,1:2,g+1)
  dumc_e(1:ng4,1:2)=v(mg(g+1)%coarse_c,1:2,g+1)
  dumc_w(1:ng4,1:2)=v(mg(g+1)%coarse_d,1:2,g+1)
    
  do k=1,2
    w(1:ng4,k)= mg(g+1)%wgt_a*dumc_n(1:ng4,k) + mg(g+1)%wgt_bc*dumc_s(1:ng4,k) &
             + mg(g+1)%wgt_bc*dumc_e(1:ng4,k) +  mg(g+1)%wgt_d*dumc_w(1:ng4,k)
  end do

  ! extension
  ! No mgbounds as the v halo has already been updated and
  ! the coarse interpolation also updates the w halo
  w(1:ng4,1:2)=v(1:ng4,1:2,g)+w(1:ng4,1:2)

  ng=mg(g)%ifull
  do i=1,itrend-1
    dumc_n(1:ng,:)=w(mg(g)%in,1:2)
    dumc_s(1:ng,:)=w(mg(g)%is,1:2)
    dumc_e(1:ng,:)=w(mg(g)%ie,1:2)
    dumc_w(1:ng,:)=w(mg(g)%iw,1:2)

    ! ocean
    ! post smoothing
    bu(1:ng)=zz(1:ng,g)+hh(1:ng,g)+yyn(1:ng,g)*dumc_n(1:ng,1)+yys(1:ng,g)*dumc_s(1:ng,1) &
                                  +yye(1:ng,g)*dumc_e(1:ng,1)+yyw(1:ng,g)*dumc_w(1:ng,1)
    cu(1:ng)=zzn(1:ng,g)*dumc_n(1:ng,1)+zzs(1:ng,g)*dumc_s(1:ng,1)+zze(1:ng,g)*dumc_e(1:ng,1)+zzw(1:ng,g)*dumc_w(1:ng,1)-rhs(1:ng,g)
    v(1:ng,1,g) = -2.*cu(1:ng)/(bu(1:ng)+sqrt(bu(1:ng)*bu(1:ng)-4.*yyz(1:ng,g)*cu(1:ng)))

    ! ice
    v(1:ng,2,g) = ( -zznice(1:ng,g)*dumc_n(1:ng,2)-zzsice(1:ng,g)*dumc_s(1:ng,2) &
                    -zzeice(1:ng,g)*dumc_e(1:ng,2)-zzwice(1:ng,g)*dumc_w(1:ng,2) &
                    +rhsice(1:ng,g) ) / zzzice(1:ng,g)

    call mgbounds_mlo(g,v(:,:,g))
    w(1:ng+mg(g)%iextra,1:2)=v(1:ng+mg(g)%iextra,1:2,g)
  end do

  dumc_n(1:ng,:)=w(mg(g)%in,1:2)
  dumc_s(1:ng,:)=w(mg(g)%is,1:2)
  dumc_e(1:ng,:)=w(mg(g)%ie,1:2)
  dumc_w(1:ng,:)=w(mg(g)%iw,1:2)

  ! ocean
  ! post smoothing
  bu(1:ng)=zz(1:ng,g)+hh(1:ng,g)+yyn(1:ng,g)*dumc_n(1:ng,1)+yys(1:ng,g)*dumc_s(1:ng,1) &
                                +yye(1:ng,g)*dumc_e(1:ng,1)+yyw(1:ng,g)*dumc_w(1:ng,1)
  cu(1:ng)=zzn(1:ng,g)*dumc_n(1:ng,1)+zzs(1:ng,g)*dumc_s(1:ng,1)+zze(1:ng,g)*dumc_e(1:ng,1)+zzw(1:ng,g)*dumc_w(1:ng,1)-rhs(1:ng,g)
  v(1:ng,1,g) = -2.*cu(1:ng)/(bu(1:ng)+sqrt(bu(1:ng)*bu(1:ng)-4.*yyz(1:ng,g)*cu(1:ng)))

  ! ice
  v(1:ng,2,g) = ( -zznice(1:ng,g)*dumc_n(1:ng,2)-zzsice(1:ng,g)*dumc_s(1:ng,2) &
                  -zzeice(1:ng,g)*dumc_e(1:ng,2)-zzwice(1:ng,g)*dumc_w(1:ng,2) &
                  +rhsice(1:ng,g) ) / zzzice(1:ng,g)

  call mgbounds_mlo(g,v(:,:,g),corner=.true.)

end do

do g=1,min(mg_maxlevel_local,1) ! same as if (mg_maxlevel_local>0) then ...
    
  ! fine grid
  call mgbcasta_mlo(2,v(:,:,2))

  ! interpolation
  ng4=mg(2)%ifull_coarse
    
  dumc_n(1:ng4,:)=v(mg(2)%coarse_a,1:2,2)
  dumc_s(1:ng4,:)=v(mg(2)%coarse_b,1:2,2)
  dumc_e(1:ng4,:)=v(mg(2)%coarse_c,1:2,2)
  dumc_w(1:ng4,:)=v(mg(2)%coarse_d,1:2,2)
    
  do k=1,2
    w(1:ng4,k)= mg(2)%wgt_a*dumc_n(1:ng4,k) + mg(2)%wgt_bc*dumc_s(1:ng4,k) &
             + mg(2)%wgt_bc*dumc_e(1:ng4,k) +  mg(2)%wgt_d*dumc_w(1:ng4,k)
  end do

end do

if (mg(1)%merge_len>1) then
  call mgbcast_mlo(1,w(:,1:2),dsolmax_g(1:2))
  ir=mod(mg(1)%merge_pos-1,mg(1)%merge_row)+1   ! index for proc row
  ic=(mg(1)%merge_pos-1)/mg(1)%merge_row+1      ! index for proc col
  do n=0,npan-1
    do jj=1,jpan
      iq_a=1+(jj-1)*ipan+n*ipan*jpan
      iq_b=jj*ipan+n*ipan*jpan
      iq_c=1+(ir-1)*ipan+(jj-1+(ic-1)*jpan)*ipan*mg(1)%merge_row+(n-1)*ipan*jpan*mg(1)%merge_len
      iq_d=ir*ipan+(jj-1+(ic-1)*jpan)*ipan*mg(1)%merge_row+(n-1)*ipan*jpan*mg(1)%merge_len
      vduma(iq_a:iq_b)=w(iq_c:iq_d,1)
      vdumb(iq_a:iq_b)=w(iq_c:iq_d,2)
    end do
    do i=1,ipan
      iq_a=i+n*ipan*jpan
      iq_c=i+(ir-1)*ipan+(ic-1)*jpan*ipan*mg(1)%merge_row+(n-1)*ipan*jpan*mg(1)%merge_len
      iq_b=is(iq_a)
      iq_d=mg(1)%is(iq_c)
      vduma(iq_b)=w(iq_d,1)
      vdumb(iq_b)=w(iq_d,2)
      iq_a=i+(jpan-1)*ipan+n*ipan*jpan
      iq_c=i+(ir-1)*ipan+(jpan-1+(ic-1)*jpan)*ipan*mg(1)%merge_row+(n-1)*ipan*jpan*mg(1)%merge_len
      iq_b=in(iq_a)
      iq_d=mg(1)%in(iq_c)
      vduma(iq_b)=w(iq_d,1)
      vdumb(iq_b)=w(iq_d,2)
    end do  
    do j=1,jpan
      iq_a=1+(j-1)*ipan+n*ipan*jpan
      iq_c=1+(ir-1)*ipan+(j-1+(ic-1)*jpan)*ipan*mg(1)%merge_row+(n-1)*ipan*jpan*mg(1)%merge_len
      iq_b=iw(iq_a)
      iq_d=mg(1)%iw(iq_c)
      vduma(iq_b)=w(iq_d,1)
      vdumb(iq_b)=w(iq_d,2)
      iq_a=ipan+(j-1)*ipan+n*ipan*jpan
      iq_c=ipan+(ir-1)*ipan+(j-1+(ic-1)*jpan)*ipan*mg(1)%merge_row+(n-1)*ipan*jpan*mg(1)%merge_len
      iq_b=ie(iq_a)
      iq_d=mg(1)%ie(iq_c)
      vduma(iq_b)=w(iq_d,1)
      vdumb(iq_b)=w(iq_d,2)
    end do
  end do
else
  ! remap mg halo to normal halo 
  vduma(1:ifull)=w(1:ifull,1)
  vdumb(1:ifull)=w(1:ifull,2)
  do n=0,npan-1
    do i=1,ipan
      iq=i+n*ipan*jpan
      iq_a=is(iq)
      iq_b=mg(1)%is(iq)
      vduma(iq_a)=w(iq_b,1)
      vdumb(iq_a)=w(iq_b,2)
      iq=i+(jpan-1)*ipan+n*ipan*jpan
      iq_a=in(iq)
      iq_b=mg(1)%in(iq)
      vduma(iq_a)=w(iq_b,1)
      vdumb(iq_a)=w(iq_b,2)
    end do  
    do j=1,jpan
      iq=1+(j-1)*ipan+n*ipan*jpan
      iq_a=iw(iq)
      iq_b=mg(1)%iw(iq)
      vduma(iq_a)=w(iq_b,1)
      vdumb(iq_a)=w(iq_b,2)
      iq=j*ipan+n*ipan*jpan        
      iq_a=ie(iq)
      iq_b=mg(1)%ie(iq)
      vduma(iq_a)=w(iq_b,1)
      vdumb(iq_a)=w(iq_b,2)
    end do
  end do
end if

! extension
neta(1:ifull+iextra)=max(neta(1:ifull+iextra)+vduma(1:ifull+iextra),-dd)*ee
ipice(1:ifull+iextra)=max(min(ipice(1:ifull+iextra)+vdumb(1:ifull+iextra),ipmax),0.) 
 
dumc(1:ifull+iextra,1)=neta
dumc(1:ifull+iextra,2)=ipice
  
do i=1,itrend
  
  ! post smoothing
  do nc=1,maxcolour

    dumc_n(1:ifullcol(nc),:)=dumc(iqn(1:ifullcol(nc),nc),:)
    dumc_s(1:ifullcol(nc),:)=dumc(iqs(1:ifullcol(nc),nc),:)
    dumc_e(1:ifullcol(nc),:)=dumc(iqe(1:ifullcol(nc),nc),:)
    dumc_w(1:ifullcol(nc),:)=dumc(iqw(1:ifullcol(nc),nc),:)
    
    ! ocean
    bu(1:ifullcol(nc))=zzhhc(1:ifullcol(nc),nc)                                                                    &
                +yync(1:ifullcol(nc),nc)*dumc_n(1:ifullcol(nc),1)+yysc(1:ifullcol(nc),nc)*dumc_s(1:ifullcol(nc),1) &
                +yyec(1:ifullcol(nc),nc)*dumc_e(1:ifullcol(nc),1)+yywc(1:ifullcol(nc),nc)*dumc_w(1:ifullcol(nc),1)
    cu(1:ifullcol(nc))=zznc(1:ifullcol(nc),nc)*dumc_n(1:ifullcol(nc),1)+zzsc(1:ifullcol(nc),nc)*dumc_s(1:ifullcol(nc),1) &
                      +zzec(1:ifullcol(nc),nc)*dumc_e(1:ifullcol(nc),1)+zzwc(1:ifullcol(nc),nc)*dumc_w(1:ifullcol(nc),1) &
                      -rhsc(1:ifullcol(nc),nc)        
    dumc(iqx(1:ifullcol(nc),nc),1)=eec(1:ifullcol(nc),nc)*max(-ddc(1:ifullcol(nc),nc),                                         &
        -2.*cu(1:ifullcol(nc))/(bu(1:ifullcol(nc))+sqrt(bu(1:ifullcol(nc))**2-4.*yyc(1:ifullcol(nc),nc)*cu(1:ifullcol(nc)))) )
    
    ! ice - cavitating fluid
    dumc(iqx(1:ifullcol(nc),nc),2)=max(0.,min(ipmaxc(1:ifullcol(nc),nc),   &
       ( -zzncice(1:ifullcol(nc),nc)*dumc_n(1:ifullcol(nc),2)              &
         -zzscice(1:ifullcol(nc),nc)*dumc_s(1:ifullcol(nc),2)              &
         -zzecice(1:ifullcol(nc),nc)*dumc_e(1:ifullcol(nc),2)              &
         -zzwcice(1:ifullcol(nc),nc)*dumc_w(1:ifullcol(nc),2)              &
        + rhscice(1:ifullcol(nc),nc) ) / zzcice(1:ifullcol(nc),nc) ))

    call bounds_colour(dumc,nc)
  end do
    
end do

call END_LOG(mgmlosetup_end)

! Main loop
do itr=2,itr_mgice

  call START_LOG(mgmlofine_begin)

  do nc=1,maxcolour
  
    dumc_n(1:ifullcol(nc),:)=dumc(iqn(1:ifullcol(nc),nc),:)
    dumc_s(1:ifullcol(nc),:)=dumc(iqs(1:ifullcol(nc),nc),:)
    dumc_e(1:ifullcol(nc),:)=dumc(iqe(1:ifullcol(nc),nc),:)
    dumc_w(1:ifullcol(nc),:)=dumc(iqw(1:ifullcol(nc),nc),:)
  
    ! ocean
    bu(1:ifullcol(nc))=zzhhc(1:ifullcol(nc),nc)                                                                    &
                +yync(1:ifullcol(nc),nc)*dumc_n(1:ifullcol(nc),1)+yysc(1:ifullcol(nc),nc)*dumc_s(1:ifullcol(nc),1) &
                +yyec(1:ifullcol(nc),nc)*dumc_e(1:ifullcol(nc),1)+yywc(1:ifullcol(nc),nc)*dumc_w(1:ifullcol(nc),1)
    cu(1:ifullcol(nc))=zznc(1:ifullcol(nc),nc)*dumc_n(1:ifullcol(nc),1)+zzsc(1:ifullcol(nc),nc)*dumc_s(1:ifullcol(nc),1) &
                      +zzec(1:ifullcol(nc),nc)*dumc_e(1:ifullcol(nc),1)+zzwc(1:ifullcol(nc),nc)*dumc_w(1:ifullcol(nc),1) &
                      -rhsc(1:ifullcol(nc),nc)        
    dumc(iqx(1:ifullcol(nc),nc),1)=eec(1:ifullcol(nc),nc)*max(-ddc(1:ifullcol(nc),nc),                                        &
       -2.*cu(1:ifullcol(nc))/(bu(1:ifullcol(nc))+sqrt(bu(1:ifullcol(nc))**2-4.*yyc(1:ifullcol(nc),nc)*cu(1:ifullcol(nc)))) )
    
    ! ice (cavitating fluid)
    dumc(iqx(1:ifullcol(nc),nc),2) = max(0.,min(ipmaxc(1:ifullcol(nc),nc), &
       ( -zzncice(1:ifullcol(nc),nc)*dumc_n(1:ifullcol(nc),2) &
         -zzscice(1:ifullcol(nc),nc)*dumc_s(1:ifullcol(nc),2) &
         -zzecice(1:ifullcol(nc),nc)*dumc_e(1:ifullcol(nc),2) &
         -zzwcice(1:ifullcol(nc),nc)*dumc_w(1:ifullcol(nc),2) &
        + rhscice(1:ifullcol(nc),nc) ) / zzcice(1:ifullcol(nc),nc) ))

    call bounds_colour(dumc,nc)

  end do
  dsol(1:ifull,1)      =dumc(1:ifull,1)-neta(1:ifull)
  dsol(1:ifull,2)      =dumc(1:ifull,2)-ipice(1:ifull)
  neta(1:ifull+iextra) =dumc(1:ifull+iextra,1)
  ipice(1:ifull+iextra)=dumc(1:ifull+iextra,2)  

  ! test for convergence
  dsolmax_g(1:2)=maxval(abs(dsol(1:ifull,1:2)),dim=1)

  dumc_n(1:ifull,:)=dumc(in,:)
  dumc_s(1:ifull,:)=dumc(is,:)
  dumc_e(1:ifull,:)=dumc(ie,:)
  dumc_w(1:ifull,:)=dumc(iw,:)

  w(1:ifull,2)= izz(:,1)+ iyy*neta(1:ifull)
  w(1:ifull,3)=izzn(:,1)+iyyn*neta(1:ifull)
  w(1:ifull,4)=izzs(:,1)+iyys*neta(1:ifull)
  w(1:ifull,5)=izze(:,1)+iyye*neta(1:ifull)
  w(1:ifull,6)=izzw(:,1)+iyyw*neta(1:ifull)
  w(1:ifull,7)=ihh+iyy*neta(1:ifull)+iyyn*dumc_n(1:ifull,1)+iyys*dumc_s(1:ifull,1)+iyye*dumc_e(1:ifull,1)+iyyw*dumc_w(1:ifull,1)

  ! residual
  w(1:ifull,1)=(-neta(1:ifull)*(     iyy*neta(1:ifull)     +iyyn*dumc_n(1:ifull,1)     +iyys*dumc_s(1:ifull,1)   &
                                                           +iyye*dumc_e(1:ifull,1)     +iyyw*dumc_w(1:ifull,1))  &
                              -(izz(:,1)*neta(1:ifull)+izzn(:,1)*dumc_n(1:ifull,1)+izzs(:,1)*dumc_s(1:ifull,1)   &
                                                      +izze(:,1)*dumc_e(1:ifull,1)+izzw(:,1)*dumc_w(1:ifull,1))  &
                              -ihh*neta(1:ifull)+irhs(:,1))*ee(1:ifull)
  where (ipice(1:ifull)>=ipmax(1:ifull))
    w(1:ifull,8)=0. ! patch to remove error when ipmax is reached - improves convergence
  elsewhere
    w(1:ifull,8)=(-(izz(:,2)*ipice(1:ifull)+izzn(:,2)*dumc_n(1:ifull,2)+izzs(:,2)*dumc_s(1:ifull,2) &
                  +izze(:,2)*dumc_e(1:ifull,2)+izzw(:,2)*dumc_w(1:ifull,2))+irhs(:,2))*ee(1:ifull)
  end where
  
  ! For when the inital grid cannot be upscaled
  call mgcollect_mlo(1,w(:,1:8),dsolmax_g)
  
  do g=1,min(mg_maxlevel_local,1) ! same as if (mg_maxlevel_local>0) then ...
  
    ! restriction
    ! (since this always operates within a panel, then ine = ien is always true)
    ng4=mg(1)%ifull_fine
    rhs(1:ng4,2)=0.25*(w(mg(1)%fine  ,1)+w(mg(1)%fine_n ,1) &
                      +w(mg(1)%fine_e,1)+w(mg(1)%fine_ne,1))
    zz(1:ng4,2) =0.25*dfac*(w(mg(1)%fine  ,2)+w(mg(1)%fine_n ,2) &
                           +w(mg(1)%fine_e,2)+w(mg(1)%fine_ne,2))
    zzn(1:ng4,2)=0.25*dfac*(w(mg(1)%fine  ,3)+w(mg(1)%fine_n ,3) &
                           +w(mg(1)%fine_e,3)+w(mg(1)%fine_ne,3))
    zzs(1:ng4,2)=0.25*dfac*(w(mg(1)%fine  ,4)+w(mg(1)%fine_n ,4) &
                           +w(mg(1)%fine_e,4)+w(mg(1)%fine_ne,4))
    zze(1:ng4,2)=0.25*dfac*(w(mg(1)%fine  ,5)+w(mg(1)%fine_n ,5) &
                           +w(mg(1)%fine_e,5)+w(mg(1)%fine_ne,5))
    zzw(1:ng4,2)=0.25*dfac*(w(mg(1)%fine  ,6)+w(mg(1)%fine_n ,6) &
                           +w(mg(1)%fine_e,6)+w(mg(1)%fine_ne,6))
    hh(1:ng4,2)    =0.25*(w(mg(1)%fine  ,7)+w(mg(1)%fine_n ,7) &
                         +w(mg(1)%fine_e,7)+w(mg(1)%fine_ne,7))
    rhsice(1:ng4,2)=0.25*(w(mg(1)%fine  ,8)+w(mg(1)%fine_n ,8) &
                         +w(mg(1)%fine_e,8)+w(mg(1)%fine_ne,8))

    ! merge grids if insufficent points on this processor
    if (mg(2)%merge_len>1) then
      w(1:ng4,1)  =rhs(1:ng4,2)
      w(1:ng4,2)  =zz(1:ng4,2)
      w(1:ng4,3)  =zzn(1:ng4,2)
      w(1:ng4,4)  =zzs(1:ng4,2)
      w(1:ng4,5)  =zze(1:ng4,2)
      w(1:ng4,6)  =zzw(1:ng4,2)
      w(1:ng4,7)  =hh(1:ng4,2)
      w(1:ng4,8)  =rhsice(1:ng4,2)
      call mgcollect_mlo(2,w(:,1:8),dsolmax_g)
      ng=mg(2)%ifull
      rhs(1:ng,2)    =w(1:ng,1)
      zz(1:ng,2)     =w(1:ng,2)
      zzn(1:ng,2)    =w(1:ng,3)
      zzs(1:ng,2)    =w(1:ng,4)
      zze(1:ng,2)    =w(1:ng,5)
      zzw(1:ng,2)    =w(1:ng,6)
      hh(1:ng,2)     =w(1:ng,7)
      rhsice(1:ng,2) =w(1:ng,8)
    end if
    
  end do
  
  call END_LOG(mgmlofine_end)

  call START_LOG(mgmloup_begin)

  ! upscale grid
  do g=2,gmax
  
    ng=mg(g)%ifull

    ! update
    ! possibly use colours here, although v is reset to zero every iteration
    ! assume zero for first guess of residual (also avoids additional bounds call)
    bu(1:ng)=zz(1:ng,g)+hh(1:ng,g)
    v(1:ng,1,g) = 2.*rhs(1:ng,g)/(bu(1:ng)+sqrt(bu(1:ng)**2+4.*yyz(1:ng,g)*rhs(1:ng,g)))
    v(1:ng,2,g) = rhsice(1:ng,g)/zzzice(1:ng,g)
    call mgbounds_mlo(g,v(:,:,g))
    
    do i=1,itrbgn-1
      dumc_n(1:ng,:)=v(mg(g)%in,1:2,g)
      dumc_s(1:ng,:)=v(mg(g)%is,1:2,g)
      dumc_e(1:ng,:)=v(mg(g)%ie,1:2,g)
      dumc_w(1:ng,:)=v(mg(g)%iw,1:2,g)

      ! ocean
      ! post smoothing
      bu(1:ng)=zz(1:ng,g)+hh(1:ng,g)+yyn(1:ng,g)*dumc_n(1:ng,1)+yys(1:ng,g)*dumc_s(1:ng,1) &
                                    +yye(1:ng,g)*dumc_e(1:ng,1)+yyw(1:ng,g)*dumc_w(1:ng,1)
      cu(1:ng)=zzn(1:ng,g)*dumc_n(1:ng,1)+zzs(1:ng,g)*dumc_s(1:ng,1)             &
              +zze(1:ng,g)*dumc_e(1:ng,1)+zzw(1:ng,g)*dumc_w(1:ng,1)-rhs(1:ng,g)
      w(1:ng,1) = -2.*cu(1:ng)/(bu(1:ng)+sqrt(bu(1:ng)**2-4.*yyz(1:ng,g)*cu(1:ng)))

      ! ice
      w(1:ng,2) = ( -zznice(1:ng,g)*dumc_n(1:ng,2)-zzsice(1:ng,g)*dumc_s(1:ng,2) &
                    -zzeice(1:ng,g)*dumc_e(1:ng,2)-zzwice(1:ng,g)*dumc_w(1:ng,2) &
                    +rhsice(1:ng,g) ) / zzzice(1:ng,g)

      call mgbounds_mlo(g,w(:,:))
      v(1:mg(g)%ifull+mg(g)%iextra,1:2,g)=w(1:mg(g)%ifull+mg(g)%iextra,1:2)
    end do
    
    ! restriction
    ! (calculate finer grid before mgcollect as the messages sent/recv are shorter)
    dumc_n(1:ng,:)=v(mg(g)%in,:,g)
    dumc_s(1:ng,:)=v(mg(g)%is,:,g)
    dumc_e(1:ng,:)=v(mg(g)%ie,:,g)
    dumc_w(1:ng,:)=v(mg(g)%iw,:,g)

    ng4=mg(g)%ifull_fine
    ws(1:ng)= zz(1:ng,g)+yyz(1:ng,g)*v(1:ng,1,g)
    zz(1:ng4,g+1)=0.25*dfac*(ws(mg(g)%fine  )+ws(mg(g)%fine_n ) &
                            +ws(mg(g)%fine_e)+ws(mg(g)%fine_ne))
    ws(1:ng)=zzn(1:ng,g)+yyn(1:ng,g)*v(1:ng,1,g)
    zzn(1:ng4,g+1)=0.25*dfac*(ws(mg(g)%fine  )+ws(mg(g)%fine_n ) &
                             +ws(mg(g)%fine_e)+ws(mg(g)%fine_ne))
    ws(1:ng)=zzs(1:ng,g)+yys(1:ng,g)*v(1:ng,1,g)
    zzs(1:ng4,g+1)=0.25*dfac*(ws(mg(g)%fine  )+ws(mg(g)%fine_n ) &
                             +ws(mg(g)%fine_e)+ws(mg(g)%fine_ne))
    ws(1:ng)=zze(1:ng,g)+yye(1:ng,g)*v(1:ng,1,g)
    zze(1:ng4,g+1)=0.25*dfac*(ws(mg(g)%fine  )+ws(mg(g)%fine_n ) &
                             +ws(mg(g)%fine_e)+ws(mg(g)%fine_ne))
    ws(1:ng)=zzw(1:ng,g)+yyw(1:ng,g)*v(1:ng,1,g)
    zzw(1:ng4,g+1)=0.25*dfac*(ws(mg(g)%fine  )+ws(mg(g)%fine_n ) &
                             +ws(mg(g)%fine_e)+ws(mg(g)%fine_ne))
    ws(1:ng)=hh(1:ng,g)+yyz(1:ng,g)*v(1:ng,1,g)+yyn(1:ng,g)*dumc_n(1:ng,1)+yys(1:ng,g)*dumc_s(1:ng,1) &
                                               +yye(1:ng,g)*dumc_e(1:ng,1)+yyw(1:ng,g)*dumc_w(1:ng,1)
    hh(1:ng4,g+1)=0.25*(ws(mg(g)%fine  )+ws(mg(g)%fine_n ) &
                       +ws(mg(g)%fine_e)+ws(mg(g)%fine_ne))

    ! ocean
    ws(1:ng)=-v(1:ng,1,g)*(yyz(1:ng,g)*v(1:ng,1,g)+yyn(1:ng,g)*dumc_n(1:ng,1)+yys(1:ng,g)*dumc_s(1:ng,1)  &
                                                  +yye(1:ng,g)*dumc_e(1:ng,1)+yyw(1:ng,g)*dumc_w(1:ng,1)) &
                          -(zz(1:ng,g)*v(1:ng,1,g)+zzn(1:ng,g)*dumc_n(1:ng,1)+zzs(1:ng,g)*dumc_s(1:ng,1)  &
                                                  +zze(1:ng,g)*dumc_e(1:ng,1)+zzw(1:ng,g)*dumc_w(1:ng,1)) &
                          -hh(1:ng,g)*v(1:ng,1,g)+rhs(1:ng,g)
    rhs(1:ng4,g+1)=0.25*(ws(mg(g)%fine  )+ws(mg(g)%fine_n ) &
                        +ws(mg(g)%fine_e)+ws(mg(g)%fine_ne))
                          
    ! ice
    ws(1:ng)=-(zznice(1:ng,g)*dumc_n(1:ng,2)+zzsice(1:ng,g)*dumc_s(1:ng,2)   &
              +zzeice(1:ng,g)*dumc_e(1:ng,2)+zzwice(1:ng,g)*dumc_w(1:ng,2))  &
              -zzzice(1:ng,g)*v(1:ng,2,g)+rhsice(1:ng,g)
    rhsice(1:ng4,g+1)=0.25*(ws(mg(g)%fine  )+ws(mg(g)%fine_n ) &
                           +ws(mg(g)%fine_e)+ws(mg(g)%fine_ne))

    ! merge grids if insufficent points on this processor
    if (mg(g+1)%merge_len>1) then
      w(1:ng4,1)  =rhs(1:ng4,g+1)
      w(1:ng4,2)  =zz(1:ng4,g+1)
      w(1:ng4,3)  =zzn(1:ng4,g+1)
      w(1:ng4,4)  =zzs(1:ng4,g+1)
      w(1:ng4,5)  =zze(1:ng4,g+1)
      w(1:ng4,6)  =zzw(1:ng4,g+1)
      w(1:ng4,7)  =hh(1:ng4,g+1)
      w(1:ng4,8)  =rhsice(1:ng4,g+1)
      call mgcollect_mlo(g+1,w(:,1:8),dsolmax_g)
      ng=mg(g+1)%ifull
      rhs(1:ng,g+1)    =w(1:ng,1)
      zz(1:ng,g+1)     =w(1:ng,2)
      zzn(1:ng,g+1)    =w(1:ng,3)
      zzs(1:ng,g+1)    =w(1:ng,4)
      zze(1:ng,g+1)    =w(1:ng,5)
      zzw(1:ng,g+1)    =w(1:ng,6)
      hh(1:ng,g+1)     =w(1:ng,7)
      rhsice(1:ng,g+1) =w(1:ng,8)
    end if

  end do

  call END_LOG(mgmloup_end)

  call START_LOG(mgmlocoarse_begin)

  ! solve coarse grid    
  do g=mg_maxlevel,mg_maxlevel_local ! same as if (mg_maxlevel==mg_maxlevel_local) then ...

    ng=mg(g)%ifull
      
    ! solve for ice using LU decomposition and back substitution with RHS
    v(1:ng,2,g)=rhsice(1:ng,g)
    call mbacksub(helm_o,v(1:ng,2,g),indy)

    ! solve non-linear water free surface with coloured SOR
    
    ! first guess
    bu(1:ng)=zz(1:ng,g)+hh(1:ng,g)
    v(1:ng,1,g) = 2.*rhs(1:ng,g)/(bu(1:ng)+sqrt(bu(1:ng)*bu(1:ng)+4.*yyz(1:ng,g)*rhs(1:ng,g)))
  
    ! pack zz,hh and rhs by colour
    do nc=1,3
      zzhhcu(:,nc)=zz(col_iq(:,nc),g)+hh(col_iq(:,nc),g)
      zzncu(:,nc)=zzn(col_iq(:,nc),g)
      zzscu(:,nc)=zzs(col_iq(:,nc),g)
      zzecu(:,nc)=zze(col_iq(:,nc),g)
      zzwcu(:,nc)=zzw(col_iq(:,nc),g)
      rhscu(:,nc)=rhs(col_iq(:,nc),g)
    end do
  
    do itrc=1,itr_mgice

      ! store previous guess for convegence test
      ws(1:ng)=v(1:ng,1,g)
 
      do nc=1,3
      
        dumc_n(1:mg_ifullmaxcol,1)=v(col_iqn(:,nc),1,g)
        dumc_s(1:mg_ifullmaxcol,1)=v(col_iqs(:,nc),1,g)
        dumc_e(1:mg_ifullmaxcol,1)=v(col_iqe(:,nc),1,g)
        dumc_w(1:mg_ifullmaxcol,1)=v(col_iqw(:,nc),1,g)
      
        bu(1:mg_ifullmaxcol)=zzhhcu(:,nc)+yyncu(:,nc)*dumc_n(1:mg_ifullmaxcol,1)+yyscu(:,nc)*dumc_s(1:mg_ifullmaxcol,1) &
                                         +yyecu(:,nc)*dumc_e(1:mg_ifullmaxcol,1)+yywcu(:,nc)*dumc_w(1:mg_ifullmaxcol,1)
        cu(1:mg_ifullmaxcol)=zzncu(:,nc)*dumc_n(1:mg_ifullmaxcol,1)+zzscu(:,nc)*dumc_s(1:mg_ifullmaxcol,1)   &
                            +zzecu(:,nc)*dumc_e(1:mg_ifullmaxcol,1)+zzwcu(:,nc)*dumc_w(1:mg_ifullmaxcol,1)   &
                            -rhscu(:,nc)
        v(col_iq(:,nc),1,g) = -2.*cu(1:mg_ifullmaxcol)/(bu(1:mg_ifullmaxcol)                      &
                              +sqrt(bu(1:mg_ifullmaxcol)**2-4.*yyzcu(:,nc)*cu(1:mg_ifullmaxcol)))
      end do
      
      ! test for convergence
      dsol(1:ng,1)=v(1:ng,1,g)-ws(1:ng)
      dsolmax(1)=maxval(abs(dsol(1:ng,1)))
      if (dsolmax(1)<tol) exit

    end do
  
  end do
  
  call END_LOG(mgmlocoarse_end)
  
  call START_LOG(mgmlodown_begin)
    
  ! downscale grid
  do g=gmax,2,-1

    call mgbcast_mlo(g+1,v(:,:,g+1),dsolmax_g(1:2))

    ! interpolation
    ng4=mg(g+1)%ifull_coarse
    
    dumc_n(1:ng4,1:2)=v(mg(g+1)%coarse_a,1:2,g+1)
    dumc_s(1:ng4,1:2)=v(mg(g+1)%coarse_b,1:2,g+1)
    dumc_e(1:ng4,1:2)=v(mg(g+1)%coarse_c,1:2,g+1)
    dumc_w(1:ng4,1:2)=v(mg(g+1)%coarse_d,1:2,g+1)
    
    do k=1,2
      w(1:ng4,k)= mg(g+1)%wgt_a*dumc_n(1:ng4,k) + mg(g+1)%wgt_bc*dumc_s(1:ng4,k) &
               + mg(g+1)%wgt_bc*dumc_e(1:ng4,k) +  mg(g+1)%wgt_d*dumc_w(1:ng4,k)
    end do

    ! extension
    ! No mgbounds as the v halo has already been updated and
    ! the coarse interpolation also updates the w halo
    w(1:ng4,1:2)=v(1:ng4,1:2,g)+w(1:ng4,1:2)

    ng=mg(g)%ifull
    do i=1,itrend-1
      dumc_n(1:ng,:)=w(mg(g)%in,1:2)
      dumc_s(1:ng,:)=w(mg(g)%is,1:2)
      dumc_e(1:ng,:)=w(mg(g)%ie,1:2)
      dumc_w(1:ng,:)=w(mg(g)%iw,1:2)

      ! ocean
      ! post smoothing
      bu(1:ng)=zz(1:ng,g)+hh(1:ng,g)+yyn(1:ng,g)*dumc_n(1:ng,1)+yys(1:ng,g)*dumc_s(1:ng,1) &
                                    +yye(1:ng,g)*dumc_e(1:ng,1)+yyw(1:ng,g)*dumc_w(1:ng,1)
      cu(1:ng)=zzn(1:ng,g)*dumc_n(1:ng,1)+zzs(1:ng,g)*dumc_s(1:ng,1)             &
              +zze(1:ng,g)*dumc_e(1:ng,1)+zzw(1:ng,g)*dumc_w(1:ng,1)-rhs(1:ng,g)
      v(1:ng,1,g) = -2.*cu(1:ng)/(bu(1:ng)+sqrt(bu(1:ng)*bu(1:ng)-4.*yyz(1:ng,g)*cu(1:ng)))

      ! ice
      v(1:ng,2,g) = ( -zznice(1:ng,g)*dumc_n(1:ng,2)-zzsice(1:ng,g)*dumc_s(1:ng,2) &
                      -zzeice(1:ng,g)*dumc_e(1:ng,2)-zzwice(1:ng,g)*dumc_w(1:ng,2) &
                      +rhsice(1:ng,g) ) / zzzice(1:ng,g)

      call mgbounds_mlo(g,v(:,:,g))
      w(1:ng+mg(g)%iextra,1:2)=v(1:ng+mg(g)%iextra,1:2,g)
    end do

    dumc_n(1:ng,:)=w(mg(g)%in,1:2)
    dumc_s(1:ng,:)=w(mg(g)%is,1:2)
    dumc_e(1:ng,:)=w(mg(g)%ie,1:2)
    dumc_w(1:ng,:)=w(mg(g)%iw,1:2)

    ! ocean
    ! post smoothing
    bu(1:ng)=zz(1:ng,g)+hh(1:ng,g)+yyn(1:ng,g)*dumc_n(1:ng,1)+yys(1:ng,g)*dumc_s(1:ng,1) &
                                  +yye(1:ng,g)*dumc_e(1:ng,1)+yyw(1:ng,g)*dumc_w(1:ng,1)
    cu(1:ng)=zzn(1:ng,g)*dumc_n(1:ng,1)+zzs(1:ng,g)*dumc_s(1:ng,1)             &
            +zze(1:ng,g)*dumc_e(1:ng,1)+zzw(1:ng,g)*dumc_w(1:ng,1)-rhs(1:ng,g)
    v(1:ng,1,g) = -2.*cu(1:ng)/(bu(1:ng)+sqrt(bu(1:ng)*bu(1:ng)-4.*yyz(1:ng,g)*cu(1:ng)))

    ! ice
    v(1:ng,2,g) = ( -zznice(1:ng,g)*dumc_n(1:ng,2)-zzsice(1:ng,g)*dumc_s(1:ng,2) &
                    -zzeice(1:ng,g)*dumc_e(1:ng,2)-zzwice(1:ng,g)*dumc_w(1:ng,2) &
                    +rhsice(1:ng,g) ) / zzzice(1:ng,g)

    call mgbounds_mlo(g,v(:,:,g),corner=.true.)

  end do

  do g=1,min(mg_maxlevel_local,1) ! same as if (mg_maxlevel_local>0) then ...
    
    ! fine grid
    call mgbcast_mlo(2,v(:,:,2),dsolmax_g(1:2))

    ! interpolation
    ng4=mg(2)%ifull_coarse
    
    dumc_n(1:ng4,:)=v(mg(2)%coarse_a,1:2,2)
    dumc_s(1:ng4,:)=v(mg(2)%coarse_b,1:2,2)
    dumc_e(1:ng4,:)=v(mg(2)%coarse_c,1:2,2)
    dumc_w(1:ng4,:)=v(mg(2)%coarse_d,1:2,2)
    
    do k=1,2
      w(1:ng4,k)= mg(2)%wgt_a*dumc_n(1:ng4,k) + mg(2)%wgt_bc*dumc_s(1:ng4,k) &
               + mg(2)%wgt_bc*dumc_e(1:ng4,k) +  mg(2)%wgt_d*dumc_w(1:ng4,k)
    end do

  end do

  call END_LOG(mgmlodown_end)

  call START_LOG(mgmlofine_begin)

  if (mg(1)%merge_len>1) then
    call mgbcast_mlo(1,w(:,1:2),dsolmax_g(1:2))
    ir=mod(mg(1)%merge_pos-1,mg(1)%merge_row)+1   ! index for proc row
    ic=(mg(1)%merge_pos-1)/mg(1)%merge_row+1      ! index for proc col
    do n=0,npan-1
      do jj=1,jpan
        iq_a=1+(jj-1)*ipan+n*ipan*jpan
        iq_b=jj*ipan+n*ipan*jpan
        iq_c=1+(ir-1)*ipan+(jj-1+(ic-1)*jpan)*ipan*mg(1)%merge_row+(n-1)*ipan*jpan*mg(1)%merge_len
        iq_d=ir*ipan+(jj-1+(ic-1)*jpan)*ipan*mg(1)%merge_row+(n-1)*ipan*jpan*mg(1)%merge_len
        vduma(iq_a:iq_b)=w(iq_c:iq_d,1)
        vdumb(iq_a:iq_b)=w(iq_c:iq_d,2)
      end do
      do i=1,ipan
        iq_a=i+n*ipan*jpan
        iq_c=i+(ir-1)*ipan+(ic-1)*jpan*ipan*mg(1)%merge_row+(n-1)*ipan*jpan*mg(1)%merge_len
        iq_b=is(iq_a)
        iq_d=mg(1)%is(iq_c)
        vduma(iq_b)=w(iq_d,1)
        vdumb(iq_b)=w(iq_d,2)
        iq_a=i+(jpan-1)*ipan+n*ipan*jpan
        iq_c=i+(ir-1)*ipan+(jpan-1+(ic-1)*jpan)*ipan*mg(1)%merge_row+(n-1)*ipan*jpan*mg(1)%merge_len
        iq_b=in(iq_a)
        iq_d=mg(1)%in(iq_c)
        vduma(iq_b)=w(iq_d,1)
        vdumb(iq_b)=w(iq_d,2)
      end do  
      do j=1,jpan
        iq_a=1+(j-1)*ipan+n*ipan*jpan
        iq_c=1+(ir-1)*ipan+(j-1+(ic-1)*jpan)*ipan*mg(1)%merge_row+(n-1)*ipan*jpan*mg(1)%merge_len
        iq_b=iw(iq_a)
        iq_d=mg(1)%iw(iq_c)
        vduma(iq_b)=w(iq_d,1)
        vdumb(iq_b)=w(iq_d,2)
        iq_a=ipan+(j-1)*ipan+n*ipan*jpan
        iq_c=ipan+(ir-1)*ipan+(j-1+(ic-1)*jpan)*ipan*mg(1)%merge_row+(n-1)*ipan*jpan*mg(1)%merge_len
        iq_b=ie(iq_a)
        iq_d=mg(1)%ie(iq_c)
        vduma(iq_b)=w(iq_d,1)
        vdumb(iq_b)=w(iq_d,2)
      end do
    end do
  else
    ! remap mg halo to normal halo 
    vduma(1:ifull)=w(1:ifull,1)
    vdumb(1:ifull)=w(1:ifull,2)
    do n=0,npan-1
      do i=1,ipan
        iq=i+n*ipan*jpan
        iq_a=is(iq)
        iq_b=mg(1)%is(iq)
        vduma(iq_a)=w(iq_b,1)
        vdumb(iq_a)=w(iq_b,2)
        iq=i+(jpan-1)*ipan+n*ipan*jpan
        iq_a=in(iq)
        iq_b=mg(1)%in(iq)
        vduma(iq_a)=w(iq_b,1)
        vdumb(iq_a)=w(iq_b,2)
      end do  
      do j=1,jpan
        iq=1+(j-1)*ipan+n*ipan*jpan
        iq_a=iw(iq)
        iq_b=mg(1)%iw(iq)
        vduma(iq_a)=w(iq_b,1)
        vdumb(iq_a)=w(iq_b,2)
        iq=j*ipan+n*ipan*jpan        
        iq_a=ie(iq)
        iq_b=mg(1)%ie(iq)
        vduma(iq_a)=w(iq_b,1)
        vdumb(iq_a)=w(iq_b,2)
      end do
    end do
  end if
 
  ! extension
  neta(1:ifull+iextra)=max(neta(1:ifull+iextra)+vduma(1:ifull+iextra),-dd)*ee
  ipice(1:ifull+iextra)=max(min(ipice(1:ifull+iextra)+vdumb(1:ifull+iextra),ipmax),0.) 
 
  dumc(1:ifull+iextra,1)=neta
  dumc(1:ifull+iextra,2)=ipice
  
  do i=1,itrend
  
    ! post smoothing
    do nc=1,maxcolour

      dumc_n(1:ifullcol(nc),:)=dumc(iqn(1:ifullcol(nc),nc),:)
      dumc_s(1:ifullcol(nc),:)=dumc(iqs(1:ifullcol(nc),nc),:)
      dumc_e(1:ifullcol(nc),:)=dumc(iqe(1:ifullcol(nc),nc),:)
      dumc_w(1:ifullcol(nc),:)=dumc(iqw(1:ifullcol(nc),nc),:)
    
      ! ocean
      bu(1:ifullcol(nc))=zzhhc(1:ifullcol(nc),nc)                                                                    &
                  +yync(1:ifullcol(nc),nc)*dumc_n(1:ifullcol(nc),1)+yysc(1:ifullcol(nc),nc)*dumc_s(1:ifullcol(nc),1) &
                  +yyec(1:ifullcol(nc),nc)*dumc_e(1:ifullcol(nc),1)+yywc(1:ifullcol(nc),nc)*dumc_w(1:ifullcol(nc),1)
      cu(1:ifullcol(nc))=zznc(1:ifullcol(nc),nc)*dumc_n(1:ifullcol(nc),1)+zzsc(1:ifullcol(nc),nc)*dumc_s(1:ifullcol(nc),1) &
                        +zzec(1:ifullcol(nc),nc)*dumc_e(1:ifullcol(nc),1)+zzwc(1:ifullcol(nc),nc)*dumc_w(1:ifullcol(nc),1) &
                        -rhsc(1:ifullcol(nc),nc)        
      dumc(iqx(1:ifullcol(nc),nc),1)=eec(1:ifullcol(nc),nc)*max(-ddc(1:ifullcol(nc),nc), &
          -2.*cu(1:ifullcol(nc))/(bu(1:ifullcol(nc))+sqrt(bu(1:ifullcol(nc))**2-4.*yyc(1:ifullcol(nc),nc)*cu(1:ifullcol(nc)))) )
    
      ! ice - cavitating fluid
      dumc(iqx(1:ifullcol(nc),nc),2)=max(0.,min(ipmaxc(1:ifullcol(nc),nc),   &
         ( -zzncice(1:ifullcol(nc),nc)*dumc_n(1:ifullcol(nc),2)              &
           -zzscice(1:ifullcol(nc),nc)*dumc_s(1:ifullcol(nc),2)              &
           -zzecice(1:ifullcol(nc),nc)*dumc_e(1:ifullcol(nc),2)              &
           -zzwcice(1:ifullcol(nc),nc)*dumc_w(1:ifullcol(nc),2)              &
          + rhscice(1:ifullcol(nc),nc) ) / zzcice(1:ifullcol(nc),nc) ))

      call bounds_colour(dumc,nc)
    end do
    dsol(1:ifull,1)      =dumc(1:ifull,1)-neta(1:ifull)
    dsol(1:ifull,2)      =dumc(1:ifull,2)-ipice(1:ifull)
    
    ! store for convergence test
    neta(1:ifull+iextra) =dumc(1:ifull,1)
    ipice(1:ifull+iextra)=dumc(1:ifull,2)
    
  end do
  
  call END_LOG(mgmlofine_end)
 
  ! test for convergence
  !if (dsolmax_g(1)<10.*tol.and.dsolmax_g(2)<10.*itol) then
  !  dsolmax(1:2)=maxval(abs(dsol(1:ifull,1:2)),dim=1)
  !  call ccmpi_allreduce(dsolmax(1:2),dsolmax_g(1:2),"max",comm_world)
  !end if  
  if (dsolmax_g(1)<tol.and.dsolmax_g(2)<itol) exit
  
end do

neta(ifull+1:ifull+iextra) =dumc(ifull+iextra,1)
ipice(ifull+1:ifull+iextra)=dumc(ifull+iextra,2)

totits     =itr
maxglobseta=dsolmax_g(1)
maxglobip  =dsolmax_g(2)

return
end subroutine mgmlo

! LU decomposition
subroutine mdecomp(a,indy)

implicit none

real, dimension(mg_minsize,mg_minsize), intent(inout) :: a
real, dimension(mg_minsize) :: vv,dumv
integer, dimension(mg_minsize), intent(out) :: indy
integer i,j,imax
integer, dimension(1) :: pos

vv(:)=1./maxval(abs(a),dim=2)
do j=1,mg_minsize-1
  do i=1,j-1
    a(i,j)=a(i,j)-sum(a(i,1:i-1)*a(1:i-1,j))
  end do
  do i=j,mg_minsize
    a(i,j)=a(i,j)-sum(a(i,1:j-1)*a(1:j-1,j))
  end do
  pos=maxloc(vv(j:mg_minsize)*abs(a(j:mg_minsize,j)))
  imax=pos(1)+j-1
  dumv(:)=a(imax,:)
  a(imax,:)=a(j,:)
  a(j,:)=dumv(:)
  vv(imax)=vv(j)
  indy(j)=imax
  a(j+1:mg_minsize,j)=a(j+1:mg_minsize,j)/a(j,j)
end do
!j=mg_minsize
do i=1,mg_minsize-1
  a(i,mg_minsize)=a(i,mg_minsize)-sum(a(i,1:i-1)*a(1:i-1,mg_minsize))
end do
!i=mg_minsize
a(mg_minsize,mg_minsize)=a(mg_minsize,mg_minsize)-sum(a(mg_minsize,1:mg_minsize-1)*a(1:mg_minsize-1,mg_minsize))
indy(mg_minsize)=mg_minsize

return
end subroutine mdecomp

! Back substitution
subroutine mbacksub(a,b,indy)

implicit none

real, dimension(mg_minsize,mg_minsize), intent(in) :: a
real, dimension(mg_minsize), intent(inout) :: b
real sumx
integer, dimension(mg_minsize), intent(in) :: indy
integer i,ii,ll

do i=1,mg_minsize
  ii=i
  ll=indy(i)
  sumx=b(ll)
  b(ll)=b(i)
  b(i)=sumx
  if (sumx/=0.) exit
end do
do i=ii+1,mg_minsize
  ll=indy(i)
  sumx=b(ll)
  b(ll)=b(i)
  b(i)=sumx-sum(a(i,ii:i-1)*b(ii:i-1))
end do
do i=mg_minsize,1,-1
  b(i)=(b(i)-sum(a(i,i+1:mg_minsize)*b(i+1:mg_minsize)))/a(i,i)
end do

return
end subroutine mbacksub


! Initialise multi-grid arrays
subroutine mgsor_init

use cc_mpi

implicit none

include 'newmpar.h'
include 'parm.h'
include 'parmdyn.h'

integer g,gp,np,iq,iqq,iql,iproc,xlen,nn,ii,jj
integer mipan,mjpan,hipan,hjpan,mil_g,iia,jja
integer i,j,n,mg_npan,mxpr,mypr,sii,eii,sjj,ejj
integer cid,ix,jx,colour,rank,ncol,nrow
integer npanx,na,nx,ny,drow,dcol,mmx,mmy
logical lglob

if (.not.sorfirst) return

! Begin full initialisation
mg_maxsize=ifull+iextra ! first guess
lglob=(nproc==1)        ! Global gather flag
mipan=ipan              ! local number of rows
mjpan=jpan              ! local number of columns
mil_g=il_g              ! global grid size
mxpr=il_g/ipan          ! number of processors over rows
mypr=il_g/jpan          ! number of processors over columns

! calculate number of levels
mg_maxlevel=1
g=1
gp=2
do while (mod(il_g,gp)==0)
  g=g+1
  gp=2*gp
end do
mg_maxlevel=g
mg_maxlevel_local=mg_maxlevel

if (myid==0) then
  write(6,*) "Initialising multi-grid arrays"
  write(6,*) "il_g,mg_maxlevel ",il_g,mg_maxlevel
end if

allocate(mg(mg_maxlevel))
mg(:)%comm_merge=0

allocate(mg_bnds(0:nproc-1,mg_maxlevel))

hipan=mipan
hjpan=mjpan

! calculate fine grid for finest grid level
g=1
mg_ifullmaxcol=0
mg_npan=npan
mg(1)%npanx=npan
mg(1)%merge_len=1
mg(1)%merge_row=1
mg(1)%nmax=1
mg(1)%ifull_coarse=0

! create processor map
allocate(mg(1)%fproc(mil_g,mil_g,0:npanels))
do n=0,npanels
  do j=1,mil_g
    do i=1,mil_g
      mg(1)%fproc(i,j,n)=fproc(i,j,n)
    end do
  end do
end do

! check if coarse grid needs mgcollect
if (mod(mipan,2)/=0.or.mod(mjpan,2)/=0.or.g==mg_maxlevel) then

  if (mod(mxpr,2)==0.and.mod(mypr,2)==0.and.g<mg_maxlevel) then
   
   ! This case occurs when there are multiple processors on a panel.
   ! Consequently, npan should be 1 or 6.
   if (npan>1.and.npan<6) then
     write(6,*) "ERROR: Invalid gather4"
     call ccmpi_abort(-1)
   end if
    
    mg(1)%merge_len=4
    mg(1)%merge_row=2
    mg(1)%nmax=4
    mxpr=mxpr/2
    mypr=mypr/2
    mipan=2*mipan
    mjpan=2*mjpan

    if (myid==0) then
      write(6,*) "Multi-grid gather4 at level           ",g,mipan,mjpan
    end if

    allocate(mg(1)%merge_list(4))

    ! find my processor and surrounding members of the gather
    nn=1-noff

    ix=-1
    jx=-1
    do jj=1,mil_g,hjpan
      do ii=1,mil_g,hipan
        if (mg(1)%fproc(ii,jj,nn)==myid) then
          ix=ii
          jx=jj
          exit
        end if
      end do
      if (ix>0) exit
    end do
    if (ix<1) then
      write(6,*) "ERROR: Cannot locate processor in gather4"
      stop
    end if
    ii=ix
    jj=jx
    if (mod((ii-1)/hipan,2)/=0) ii=ii-hipan
    if (mod((jj-1)/hjpan,2)/=0) jj=jj-hjpan
    ix=ix-ii ! offset for myid from ii
    jx=jx-jj ! offset for myid from jj
   
    mg(1)%merge_list(1)=mg(1)%fproc(ii,      jj      ,nn)
    mg(1)%merge_list(2)=mg(1)%fproc(ii+hipan,jj      ,nn)
    mg(1)%merge_list(3)=mg(1)%fproc(ii,      jj+hjpan,nn)
    mg(1)%merge_list(4)=mg(1)%fproc(ii+hipan,jj+hjpan,nn)
       
    do j=1,mil_g,mjpan
      do i=1,mil_g,mipan
        cid=mg(1)%fproc(i+ix,j+jx,nn) ! processor in same merge position as myid
                                      ! we will maintain communications with this processor
        do jja=1,mjpan
          do iia=1,mipan
            ! update fproc map with processor that owns this data
            mg(1)%fproc(i+iia-1,j+jja-1,nn)=cid
          end do
        end do
      end do
    end do

    mg(1)%merge_pos=-1
    do j=1,4
      if (mg(1)%merge_list(j)==myid) then
        mg(1)%merge_pos=j
        exit
      end if
    end do
    if (mg(g)%merge_pos<1) then
      write(6,*) "ERROR: Invalid merge_pos g,pos ",g,mg(g)%merge_pos
      call ccmpi_abort(-1)
    end if
      
    ! fix any remaining panels
    if (npan==1) then
      ! must loop over all panels
      do nn=0,npanels
        do j=1,mil_g,mjpan
          do i=1,mil_g,mipan
            cid=mg(1)%fproc(i+ix,j+jx,nn) ! processor in same merge position as myid
                                          ! we will maintain communications with this processor
            do jja=1,mjpan
              do iia=1,mipan
                ! update fproc map with processor that owns this data
                mg(1)%fproc(i+iia-1,j+jja-1,nn)=cid
              end do
            end do
          end do
        end do
      end do
    else
      do n=0,npanels
        mg(1)%fproc(:,:,n)=mg(1)%fproc(:,:,nn)
      end do    
    end if

    ! define local comm for gather
    colour=mg(1)%merge_list(1)
    rank=mg(1)%merge_pos-1
    call ccmpi_commsplit(mg(1)%comm_merge,comm_world,colour,rank)
#ifdef idleproc
    if (rank/=0) then
      mg_maxlevel_local=0
      mg(1)%nmax=0
    end if
#endif
  
  else
    write(6,*) "ERROR: Grid g=1 requires gatherall for multi-grid solver"
    call ccmpi_abort(-1)
  end if
else
  if (myid==0) then
    write(6,*) "Multi-grid fine level                 ",g,mipan,mjpan
  end if
end if

mg(1)%ipan=mipan
mg(1)%ifull=mipan*mjpan*mg_npan
mg(1)%ifull_fine=mg(1)%ifull/4

np=mg(1)%ifull
allocate(mg(1)%in(np),mg(1)%ie(np),mg(1)%is(np),mg(1)%iw(np))
allocate(mg(1)%ine(np),mg(1)%inw(np),mg(1)%ise(np),mg(1)%isw(np))

call mg_index(1,mil_g,mipan,mjpan)

if (mg_maxlevel>1) then
  allocate(mg(1)%fine(mg(1)%ifull_fine),mg(1)%fine_n(mg(1)%ifull_fine))
  allocate(mg(1)%fine_e(mg(1)%ifull_fine),mg(1)%fine_ne(mg(1)%ifull_fine))
  iqq=0
  do n=1,mg_npan
    do jj=1,mjpan,2
      do ii=1,mipan,2
        iq=indx(ii,jj,n-1,mipan,mjpan)
        iqq=iqq+1
        mg(1)%fine(iqq)=iq
      end do
    end do
  end do
  mg(1)%fine_n=mg(1)%in(mg(1)%fine)
  mg(1)%fine_e=mg(1)%ie(mg(1)%fine)
  mg(1)%fine_ne=mg(1)%ine(mg(1)%fine)
end if

mg_maxsize=max(mg_maxsize,mg(1)%ifull+mg(1)%iextra)

! loop over upscaled grids
do g=2,mg_maxlevel

  mipan=mipan/2
  mjpan=mjpan/2
  mil_g=mil_g/2
  hipan=mipan
  hjpan=mjpan

  ! assign processors to each grid point
  ! Below are default values for fproc which are modified
  ! when we collect data over processors
  allocate(mg(g)%fproc(mil_g,mil_g,0:npanels))

  ! Calculate size of grid at this level
  do nn=0,npanels
    do jj=1,2*mil_g,2
      jja=(jj-1)/2+1
      do ii=1,2*mil_g,2
        iia=(ii-1)/2+1
        mg(g)%fproc(iia,jja,nn)=mg(g-1)%fproc(ii,jj,nn)
      end do
    end do
  end do

  ! default if no gather for upscaled grid
  mg(g)%npanx=npan
  mg(g)%merge_len=1
  mg(g)%merge_row=1
  mg(g)%nmax=1
  
  ! check for multi-grid gather
  if (mod(mipan,2)/=0.or.mod(mjpan,2)/=0.or.g==mg_maxlevel) then ! grid cannot be subdivided on current processor
 
    if (mod(mxpr,2)==0.and.mod(mypr,2)==0.and.g<mg_maxlevel) then ! collect data over adjacent processors (proc=4)

      ! This case occurs when there are multiple processors on a panel.
      ! Consequently, npan should be 1 or 6.
      if (npan>1.and.npan<6) then
        write(6,*) "ERROR: Invalid gather4"
        call ccmpi_abort(-1)
      end if

      mg(g)%merge_len=4
      mg(g)%merge_row=2
      mg(g)%nmax=4
      mxpr=mxpr/2
      mypr=mypr/2
      mipan=2*mipan
      mjpan=2*mjpan

      if (myid==0) then
        write(6,*) "Multi-grid gather4 at level           ",g,mipan,mjpan
      end if

      allocate(mg(g)%merge_list(4))
      
      nn=1-noff

      ! find my processor and surrounding members of the gather
      ix=-1
      jx=-1
      do jj=1,mil_g,hjpan
        do ii=1,mil_g,hipan
          if (mg(g)%fproc(ii,jj,nn)==myid) then
            ix=ii
            jx=jj
            exit
          end if
        end do
        if (ix>0) exit
      end do
      if (ix<1) then
        write(6,*) "ERROR: Cannot locate processor in gather4"
        call ccmpi_abort(-1)
      end if
      ii=ix
      jj=jx
      mmx=mod((ii-1)/hipan,2)
      mmy=mod((jj-1)/hjpan,2)
      ii=ii-mmx*hipan ! left corner of merge
      jj=jj-mmy*hjpan ! bottom corner of merge
      ix=ix-ii ! offset for myid from ii
      jx=jx-jj ! offset for myid from jj
       
      mg(g)%merge_list(1)=mg(g)%fproc(ii,      jj      ,nn)
      mg(g)%merge_list(2)=mg(g)%fproc(ii+hipan,jj      ,nn)
      mg(g)%merge_list(3)=mg(g)%fproc(ii,      jj+hjpan,nn)
      mg(g)%merge_list(4)=mg(g)%fproc(ii+hipan,jj+hjpan,nn)
 
      do j=1,mil_g,mjpan
        do i=1,mil_g,mipan
          cid=mg(g)%fproc(i+ix,j+jx,nn) ! processor in same merge position as myid
                                        ! we will maintain communications with this processor
          do jja=1,mjpan
            do iia=1,mipan
              ! update fproc map with processor that owns this data
              mg(g)%fproc(i+iia-1,j+jja-1,nn)=cid
            end do
          end do
        end do
      end do

      mg(g)%merge_pos=-1
      do j=1,4
        if (mg(g)%merge_list(j)==myid) then
          mg(g)%merge_pos=j
          exit
        end if
      end do
       
      if (mg(g)%merge_pos<1) then
        write(6,*) "ERROR: Invalid merge_pos g,pos ",g,mg(g)%merge_pos
        call ccmpi_abort(-1)
      end if

      ! fix any remaining panels
      if (npan==1) then
        ! must loop over all panels
        do nn=0,npanels
          do j=1,mil_g,mjpan
            do i=1,mil_g,mipan
              cid=mg(g)%fproc(i+ix,j+jx,nn) ! processor in same merge position as myid
                                            ! we will maintain communications with this processor
              do jja=1,mjpan
                do iia=1,mipan
                  ! update fproc map with processor that owns this data
                  mg(g)%fproc(i+iia-1,j+jja-1,nn)=cid
                end do
              end do
            end do
          end do
        end do
      else
        do n=0,npanels
          mg(g)%fproc(:,:,n)=mg(g)%fproc(:,:,nn)
        end do
      end if
    
    else if (.not.lglob) then ! collect all data to one processor
      lglob=.true.
#ifdef uniform_decomp
      mg(g)%merge_len=mxpr*mypr
#else
      mg(g)%merge_len=min(6*mxpr*mypr,nproc)
      mg(g)%npanx=1
      mg_npan=6
#endif

      mg(g)%merge_row=mxpr
      mg(g)%nmax=mg(g)%merge_len
      mipan=mipan*mxpr
      mjpan=mjpan*mypr
      mxpr=1
      mypr=1

      if (myid==0) then
        write(6,*) "Multi-grid gatherall at level         ",g,mipan,mjpan
      end if
      
      ! find gather members
#ifdef uniform_decomp
      allocate(mg(g)%merge_list(mg(g)%merge_len))
      iqq=0
      do jj=1,mil_g,hjpan
        do ii=1,mil_g,hipan
          iqq=iqq+1
          mg(g)%merge_list(iqq)=mg(g)%fproc(ii,jj,0)
        end do
      end do
      if (iqq/=mg(g)%merge_len) then
        write(6,*) "ERROR: merge_len mismatch ",iqq,mg(g)%merge_len,g
        stop
      end if
      mg(g)%merge_pos=-1
      do i=1,mg(g)%merge_len
        if (mg(g)%merge_list(i)==myid) then
          mg(g)%merge_pos=i
          exit
        end if
      end do
      if (mg(g)%merge_pos<1) then
        write(6,*) "ERROR: Invalid merge_pos g,pos ",g,mg(g)%merge_pos
        call ccmpi_abort(-1)
      end if
#else
      allocate(mg(g)%merge_list(mg(g)%merge_len))      
      iqq=0
      do n=1,6/npan
        nn=(n-1)*npan
        do jj=1,mil_g,hjpan
          do ii=1,mil_g,hipan
            iqq=iqq+1
            mg(g)%merge_list(iqq)=mg(g)%fproc(ii,jj,nn)
          end do
        end do
      end do
      if (iqq/=mg(g)%merge_len) then
        write(6,*) "ERROR: merge_len mismatch ",iqq,mg(g)%merge_len,g
        stop
      end if
      mg(g)%merge_pos=-1
      do i=1,mg(g)%merge_len
        if (mg(g)%merge_list(i)==myid) then
          mg(g)%merge_pos=i
          exit
        end if
      end do
      if (mg(g)%merge_pos<1) then
        write(6,*) "ERROR: Invalid merge_pos g,pos ",g,mg(g)%merge_pos
        call ccmpi_abort(-1)
      end if
#endif

      ! modify fproc for remaining processor
      mg(g)%fproc(:,:,:)=myid

    else ! all data is already on one processor
      if (g/=mg_maxlevel) then
        write(6,*) "ERROR: g/=mg_maxlevel ",g,mg_maxlevel
        call ccmpi_abort(-1)
      end if
      if (myid==0) then
        write(6,*) "Multi-grid toplevel                   ",g,mipan,mjpan
      end if
      mg(g)%merge_pos=1
    end if

    ! define split and local comm
    if (mg(g)%merge_len>1) then
      colour=mg(g)%merge_list(1)
      rank=mg(g)%merge_pos-1
      call ccmpi_commsplit(mg(g)%comm_merge,comm_world,colour,rank)
#ifdef idleproc
      if (rank/=0) then
        mg_maxlevel_local=min(g-1,mg_maxlevel_local)
        mg(g)%nmax=0
      end if
#endif
    end if
  
  else
    if (myid==0) then
      write(6,*) "Multi-grid local subdivision at level ",g,mipan,mjpan
    end if
    ! no messages sent, but we allocate this arrays for the
    ! coarse calculation below
    mg(g)%merge_pos=1
  end if
  

  ! total number of points over all panels on this processor
  mg(g)%ipan=mipan
  mg(g)%ifull=mipan*mjpan*mg_npan
  mg(g)%ifull_fine=mg(g)%ifull/4  
  nrow=mg(g)%merge_row
  ncol=mg(g)%merge_len/nrow
  drow=mipan/nrow
  dcol=mjpan/ncol
  npanx=mg_npan
  
#ifndef uniform_decomp
  if (lglob) then
    npanx=1
    dcol=6*mjpan/ncol
  end if
#endif
  
  np=mg(g)%ifull
  allocate(mg(g)%in(np),mg(g)%ie(np),mg(g)%is(np),mg(g)%iw(np))
  allocate(mg(g)%ine(np),mg(g)%inw(np),mg(g)%ise(np),mg(g)%isw(np))

  call mg_index(g,mil_g,mipan,mjpan)

  ! ifine is the index of the SW point on the next finer grid.
  if (g<mg_maxlevel) then
    np=mg(g)%ifull_fine
    allocate(mg(g)%fine(np),mg(g)%fine_n(np))
    allocate(mg(g)%fine_e(np),mg(g)%fine_ne(np))
    ! mipan and mjpan should always be an even number here
    iqq=0
    do n=1,mg_npan
      do jj=1,mjpan,2
        do ii=1,mipan,2
          iq=indx(ii,jj,n-1,mipan,mjpan)
          iqq=iqq+1
          mg(g)%fine(iqq)=iq
        end do
      end do
    end do
    mg(g)%fine_n=mg(g)%in(mg(g)%fine)
    mg(g)%fine_e=mg(g)%ie(mg(g)%fine)
    mg(g)%fine_ne=mg(g)%ine(mg(g)%fine)
  end if
  
  mg_maxsize=max(mg_maxsize,mg(g)%ifull+mg(g)%iextra)

  ! Set up pointers for neighbours on the next coarser grid.
  ! ic1 is the nearest neigbour, ic2 and ic3 are equally distant and
  ! ic4 is the most distant.
  ! These are best defined separately on each face. The code here
  ! requires an even number of points to ensure the coarser grid
  ! always sits inside the fine grid.

  mg(g)%ifull_coarse=mg(g-1)%ifull+mg(g-1)%ixlen
  np=mg(g)%ifull_coarse
  allocate(mg(g)%coarse_a(np),mg(g)%coarse_b(np),mg(g)%coarse_c(np),mg(g)%coarse_d(np))
  allocate(mg(g)%wgt_a(np),mg(g)%wgt_bc(np),mg(g)%wgt_d(np))
  mg(g)%coarse_a=0 ! unassigned
  ! default weights
  mg(g)%wgt_a=0.5625
  mg(g)%wgt_bc=0.1875
  mg(g)%wgt_d=0.0625
 
  do n=1,npanx
  
    na=mg(g)%merge_pos
    nx=mod(na-1,nrow)+1
    ny=(na-1)/nrow+1
    sii=(nx-1)*drow+1
    eii=nx*drow
    sjj=(ny-1)*dcol+1
    ejj=ny*dcol
    
    do jj=sjj,ejj
      jja=2*(jj-sjj)+1
      do ii=sii,eii
        iia=2*(ii-sii)+1
       
        iqq=indx(ii,jj,n-1,mipan,mjpan)   ! coarse grid
        
        ! odd, odd          
        iq =indx(iia,jja,n-1,2*drow,2*dcol)     ! fine grid
        mg(g)%coarse_a(iq)=          iqq
        mg(g)%coarse_b(iq)= mg(g)%is(iqq)
        mg(g)%coarse_c(iq)= mg(g)%iw(iqq)
        mg(g)%coarse_d(iq)=mg(g)%isw(iqq)
        
        ! odd, even
        iq =indx(iia,jja+1,n-1,2*drow,2*dcol)   ! fine grid
        mg(g)%coarse_a(iq)=          iqq
        mg(g)%coarse_b(iq)= mg(g)%in(iqq)
        mg(g)%coarse_c(iq)= mg(g)%iw(iqq)
        mg(g)%coarse_d(iq)=mg(g)%inw(iqq)
          
        ! even, odd
        iq =indx(iia+1,jja,n-1,2*drow,2*dcol)   ! fine grid 
        mg(g)%coarse_a(iq)=          iqq
        mg(g)%coarse_b(iq)= mg(g)%is(iqq)
        mg(g)%coarse_c(iq)= mg(g)%ie(iqq)
        mg(g)%coarse_d(iq)=mg(g)%ise(iqq)

        ! even, even
        iq =indx(iia+1,jja+1,n-1,2*drow,2*dcol) ! fine grid
        mg(g)%coarse_a(iq)=          iqq
        mg(g)%coarse_b(iq)= mg(g)%in(iqq)
        mg(g)%coarse_c(iq)= mg(g)%ie(iqq)
        mg(g)%coarse_d(iq)=mg(g)%ine(iqq)
      
      end do
    end do

    ! boundaries
    ! Here we update the boundaries using the coarse
    ! array which avoids an extra call to mgbounds
    if (mg(g-1)%ixlen>0) then 
    
      ! need to check every point as the current
      ! grid may be the result of a global gather
      do jj=sjj,ejj
        jja=2*(jj-sjj)+1
        do ii=sii,eii
          iia=2*(ii-sii)+1
        
          iqq=indx(ii,jj,n-1,mipan,mjpan)   ! coarse grid
        
          ! odd, odd          
          iq =indx(iia,jja,n-1,2*drow,2*dcol)     ! fine grid
          iql=mg(g-1)%is(iq)
          if (mg(g)%coarse_a(iql)==0) then
            mg(g)%coarse_a(iql)= mg(g)%is(iqq)
            mg(g)%coarse_b(iql)=          iqq
            mg(g)%coarse_c(iql)=mg(g)%isw(iqq)
            mg(g)%coarse_d(iql)= mg(g)%iw(iqq)
          end if
          iql=mg(g-1)%iw(iq)
          if (mg(g)%coarse_a(iql)==0) then
            mg(g)%coarse_a(iql)= mg(g)%iw(iqq)
            mg(g)%coarse_b(iql)=          iqq
            mg(g)%coarse_c(iql)=mg(g)%isw(iqq)
            mg(g)%coarse_d(iql)= mg(g)%is(iqq)
          end if
          
          ! odd, even
          iq =indx(iia,jja+1,n-1,2*drow,2*dcol)   ! fine grid
          iql=mg(g-1)%in(iq)
          if (mg(g)%coarse_a(iql)==0) then
            mg(g)%coarse_a(iql)= mg(g)%in(iqq)
            mg(g)%coarse_b(iql)=          iqq
            mg(g)%coarse_c(iql)=mg(g)%inw(iqq)
            mg(g)%coarse_d(iql)= mg(g)%iw(iqq)
          end if
          iql=mg(g-1)%iw(iq)
          if (mg(g)%coarse_a(iql)==0) then
            mg(g)%coarse_a(iql)= mg(g)%iw(iqq)
            mg(g)%coarse_b(iql)=          iqq
            mg(g)%coarse_c(iql)=mg(g)%inw(iqq)
            mg(g)%coarse_d(iql)= mg(g)%in(iqq)
          end if
          
          ! even, odd
          iq =indx(iia+1,jja,n-1,2*drow,2*dcol)   ! fine grid 
          iql=mg(g-1)%is(iq)
          if (mg(g)%coarse_a(iql)==0) then
            mg(g)%coarse_a(iql)= mg(g)%is(iqq)
            mg(g)%coarse_b(iql)=          iqq
            mg(g)%coarse_c(iql)=mg(g)%ise(iqq)
            mg(g)%coarse_d(iql)= mg(g)%ie(iqq)
          end if
          iql=mg(g-1)%ie(iq)
          if (mg(g)%coarse_a(iql)==0) then
            mg(g)%coarse_a(iql)= mg(g)%ie(iqq)
            mg(g)%coarse_b(iql)=          iqq
            mg(g)%coarse_c(iql)=mg(g)%ise(iqq)
            mg(g)%coarse_d(iql)= mg(g)%is(iqq)
          end if
          
          ! even, even
          iq =indx(iia+1,jja+1,n-1,2*drow,2*dcol) ! fine grid
          iql=mg(g-1)%in(iq)
          if (mg(g)%coarse_a(iql)==0) then
            mg(g)%coarse_a(iql)= mg(g)%in(iqq)
            mg(g)%coarse_b(iql)=          iqq
            mg(g)%coarse_c(iql)=mg(g)%ine(iqq)
            mg(g)%coarse_d(iql)= mg(g)%ie(iqq)
          end if
          iql=mg(g-1)%ie(iq)
          if (mg(g)%coarse_a(iql)==0) then
            mg(g)%coarse_a(iql)= mg(g)%ie(iqq)
            mg(g)%coarse_b(iql)=          iqq
            mg(g)%coarse_c(iql)=mg(g)%ine(iqq)
            mg(g)%coarse_d(iql)= mg(g)%in(iqq)
          end if
      
        end do
      end do
    end if
    
  end do
  
  ! adjust weights for panel corners
  do iq=1,np
    if (mg(g)%coarse_d(iq)==mg(g)%coarse_b(iq).or.mg(g)%coarse_d(iq)==mg(g)%coarse_c(iq)) then
      mg(g)%wgt_a(iq)=0.5
      mg(g)%wgt_bc(iq)=0.25
      mg(g)%wgt_d(iq)=0.
    else if (mg(g)%coarse_c(iq)==mg(g)%coarse_a(iq)) then
      iqq=mg(g)%coarse_d(iq)
      mg(g)%coarse_c(iq)=mg(g)%coarse_d(iq)
      mg(g)%coarse_d(iq)=iqq
      mg(g)%wgt_a(iq)=0.5
      mg(g)%wgt_bc(iq)=0.25
      mg(g)%wgt_d(iq)=0.
    else if (mg(g)%coarse_c(iq)==mg(g)%coarse_d(iq)) then
      mg(g)%wgt_a(iq)=0.5
      mg(g)%wgt_bc(iq)=0.25
      mg(g)%wgt_d(iq)=0.
    end if
  end do

end do

gmax=min(mg_maxlevel-1,mg_maxlevel_local)
mg_minsize=0
if (mg_maxlevel_local==mg_maxlevel) then
  mg_minsize=6*mil_g*mil_g
end if

! free some memory
do g=1,mg_maxlevel
  deallocate(mg(g)%fproc)
  if (mg(g)%merge_len>1) then
    deallocate(mg(g)%merge_list)
  end if
end do

do g=gmax+2,mg_maxlevel
  deallocate(mg(g)%in,mg(g)%ie,mg(g)%is,mg(g)%iw)
  deallocate(mg(g)%ine,mg(g)%inw,mg(g)%ise,mg(g)%isw)
  deallocate(mg(g)%coarse_a,mg(g)%coarse_b,mg(g)%coarse_c,mg(g)%coarse_d)
  deallocate(mg(g)%wgt_a,mg(g)%wgt_bc,mg(g)%wgt_d)
end do

do g=gmax+1,mg_maxlevel-1
  deallocate(mg(g)%fine,mg(g)%fine_n)
  deallocate(mg(g)%fine_e,mg(g)%fine_ne)
end do

sorfirst=.false.
if (myid==0) then
  write(6,*) "Finished initialising multi-grid arrays"
end if

return
end subroutine mgsor_init

subroutine mgzz_init(zz,zzn,zze,zzw,zzs)

use cc_mpi

implicit none

include 'newmpar.h'

integer g,np,iqq,iq
real, dimension(ifull), intent(in) :: zz,zzn,zze,zzw,zzs
real, dimension(5,mg_maxsize) :: dum

  ! Later equations use zz/delta**2. It's more efficient to put
  ! this scaling in here. delta = delta_0 * 2**(m-1) so increases
  ! by a factor of 2 each grid. Therefore using a factor of 1/4
  ! in the averaging calculation will compensate.

if (myid==0) then
  write(6,*) "Initialising atmosphere multi-grid coupling arrays"
end if

if (zzfirst) then
  do g=1,gmax+1
    np=mg(g)%ifull
    allocate(mg(g)%zzn(np),mg(g)%zze(np),mg(g)%zzs(np),mg(g)%zzw(np),mg(g)%zz(np))
  end do
  zzfirst=.false.
end if

! no need for bounds call as all indices lie within this processor
g=1
np=mg(1)%ifull

dum(1,1:ifull)=zzn
dum(2,1:ifull)=zzs
dum(3,1:ifull)=zze
dum(4,1:ifull)=zzw
dum(5,1:ifull)=zz
call mgcollect(1,dum)
mg(1)%zzn(:)=dum(1,1:np)
mg(1)%zzs(:)=dum(2,1:np)
mg(1)%zze(:)=dum(3,1:np)
mg(1)%zzw(:)=dum(4,1:np)
mg(1)%zz(:) =dum(5,1:np)

do g=2,gmax+1
  np=mg(g)%ifull
  do iq=1,mg(g-1)%ifull_fine
    dum(1,iq)=0.25*dfac*(mg(g-1)%zzn(mg(g-1)%fine(iq)  )+mg(g-1)%zzn(mg(g-1)%fine_n(iq)) &
                        +mg(g-1)%zzn(mg(g-1)%fine_e(iq))+mg(g-1)%zzn(mg(g-1)%fine_ne(iq)))
    dum(2,iq)=0.25*dfac*(mg(g-1)%zzs(mg(g-1)%fine(iq)  )+mg(g-1)%zzs(mg(g-1)%fine_n(iq)) &
                        +mg(g-1)%zzs(mg(g-1)%fine_e(iq))+mg(g-1)%zzs(mg(g-1)%fine_ne(iq)))
    dum(3,iq)=0.25*dfac*(mg(g-1)%zze(mg(g-1)%fine(iq)  )+mg(g-1)%zze(mg(g-1)%fine_n(iq)) &
                        +mg(g-1)%zze(mg(g-1)%fine_e(iq))+mg(g-1)%zze(mg(g-1)%fine_ne(iq)))
    dum(4,iq)=0.25*dfac*(mg(g-1)%zzw(mg(g-1)%fine(iq)  )+mg(g-1)%zzw(mg(g-1)%fine_n(iq)) &
                        +mg(g-1)%zzw(mg(g-1)%fine_e(iq))+mg(g-1)%zzw(mg(g-1)%fine_ne(iq)))
    dum(5,iq)=0.25*dfac*(mg(g-1)%zz(mg(g-1)%fine(iq)  ) +mg(g-1)%zz(mg(g-1)%fine_n(iq))  &
                        +mg(g-1)%zz(mg(g-1)%fine_e(iq)) +mg(g-1)%zz(mg(g-1)%fine_ne(iq)))
  end do
  call mgcollect(g,dum)
  mg(g)%zzn=dum(1,1:np)
  mg(g)%zzs=dum(2,1:np)
  mg(g)%zze=dum(3,1:np)
  mg(g)%zzw=dum(4,1:np)
  mg(g)%zz =dum(5,1:np)
end do

if (myid==0) then
  write(6,*) "Finished initialising atmosphere multi-grid coupling arrays"
end if

return
end subroutine mgzz_init

end module mgsolve
