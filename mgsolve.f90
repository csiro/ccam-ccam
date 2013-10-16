module mgsolve

! This module solves for the atmosphere helmholtz equation and the ocean free
! surface equation using a multi-grid approach.

! Currently processors become idle as we upscale the grid until eventually
! only one processor solves the coarse grid.  This avoids redudant messages
! and helps reduce the amount of bandwidth employed by the scheme.

! Note that the solver can be greatly optimised for grids that permit many
! local subdivitions.  For example if the number of grid points on a processor
! is M x M, where M = 2^N and N is an integer greater than 1, then the
! multi-grid solver will be considerably more efficent as more processors can
! be applied in parallel.  Example grids which typically work well include
! C48, C96, C192, C384, C768 and C1536.  An alternate family of grids includes
! C32, C64, C128, C256, C512 and C1024.

implicit none

private
public mghelm,mgmlo,mgsor_init,mgzz_init

integer, save :: mg_maxsize,mg_minsize,gmax

integer, parameter :: itr_mg   =20
integer, parameter :: itr_mgice=20

real, parameter :: dfac=0.25 ! adjustment for grid spacing

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
integer itrc,itr,ng,ng4,g,gb,k,jj,i,j,iq
integer klimc,knew,klim,ir,ic
integer nc,n,iq_a,iq_b,iq_c,iq_d
real, dimension(ifull+iextra,kl), intent(inout) :: iv
real, dimension(ifull,kl), intent(in) :: ihelm,jrhs
real, dimension(ifull,kl) :: irhs
real, dimension(ifull), intent(in) :: izz,izzn,izze,izzw,izzs
real, dimension(ifullx,kl,maxcolour) :: rhelmc,rhsc
real, dimension(ifullx,maxcolour) :: zznc,zzec,zzwc,zzsc
real, dimension(ifull+iextra,kl) :: vdum
real, dimension(mg_maxsize,kl,gmax+1) :: v
real, dimension(mg_maxsize,kl,2:gmax+1) :: rhs
real, dimension(mg_maxsize,kl,gmax+1) :: helm
real, dimension(mg_maxsize,kl) :: w,dsol
real, dimension(mg_minsize,mg_minsize) :: helm_m
real, dimension(mg_minsize,mg_minsize,kl) :: helm_o
real, dimension(mg_maxsize) :: ws
real, dimension(2,kl) :: smaxmin_g
real, dimension(kl) :: dsolmax_g,savg

call start_log(helm_begin)

if (sorfirst.or.zzfirst) then
  write(6,*) "ERROR: mghelm requires mgsor_init and mgzz_init to be called first"
  call ccmpi_abort(-1)
end if

! zz*(DIV^2 v) - helm*v = rhs
! zz*(DIV^2 vd+v0) - helm*(vd+v0) = rhs
! zz*(DIV^2 vd) - helm*vd = rhs + helm*v0 - zz*(DIV^2 v0)

! determine max/min for convergence calculations
klim=kl
do k=1,kl
  smaxmin_g(1,k)=maxval(iv(1:ifull,k))
  smaxmin_g(2,k)=minval(iv(1:ifull,k))
end do

! upscale RHS terms
helm(1:ifull,:,1)=ihelm(1:ifull,:)
call mgcollectxn(1,helm(:,:,1),smaxmin_g(:,:))
do g=1,gmax
  ng4=mg(g)%ifull_fine
  helm(1:ng4,1:kl,g+1)=0.25*(helm(mg(g)%fine  ,1:kl,g)+helm(mg(g)%fine_n ,1:kl,g) &
                            +helm(mg(g)%fine_e,1:kl,g)+helm(mg(g)%fine_ne,1:kl,g))
  call mgcollectxn(g+1,helm(:,:,g+1),smaxmin_g)
end do
! store data for LU decomposition of coarse grid
do g=mg_maxlevel,mg_maxlevel_local ! same as if (mg_maxlevel_local==mg_maxlevel) then ...
  ng=mg(g)%ifull
  helm_o(:,:,:)=0.
  do k=1,kl
    do iq=1,ng
      helm_o(mg(g)%in(iq),iq,k)=mg(g)%zzn(iq)
      helm_o(mg(g)%is(iq),iq,k)=mg(g)%zzs(iq)
      helm_o(mg(g)%ie(iq),iq,k)=mg(g)%zze(iq)
      helm_o(mg(g)%iw(iq),iq,k)=mg(g)%zzw(iq)
      helm_o(iq,iq,k)=mg(g)%zz(iq)-helm(iq,k,g)
    end do
    call mdecomp(helm_o(:,:,k),indy(:,k)) ! destroys helm_m
  end do
end do
do g=gmax,0,-1
  call mgbcastxn(g+1,smaxmin_g)
end do

! remove offsets
do k=1,kl
  savg(k)=0.5*(smaxmin_g(1,k)+smaxmin_g(2,k))
  iv(1:ifull,k)=iv(1:ifull,k)-savg(k)
  irhs(:,k)=jrhs(:,k)+(ihelm(:,k)-izz-izzn-izzs-izze-izzw)*savg(k)
end do

! solver assumes boundaries are updated
call bounds(iv)

! pack colour arrays at fine level
do nc=1,maxcolour
  zznc(:,nc) =izzn(iqx(:,nc))
  zzwc(:,nc) =izzw(iqx(:,nc))
  zzec(:,nc) =izze(iqx(:,nc))
  zzsc(:,nc) =izzs(iqx(:,nc))
  do k=1,kl
    rhelmc(:,k,nc)=1./(ihelm(iqx(:,nc),k)-izz(iqx(:,nc)))
    rhsc(:,k,nc)  =irhs(iqx(:,nc),k)
  end do
end do

! Main loop
iters=0
do itr=1,itr_mg

  ! update on model grid using colours
  do nc=1,maxcolour
    do k=1,klim
      dsol(iqx(:,nc),k) = ( zznc(:,nc)*iv(iqn(:,nc),k) + zzwc(:,nc)*iv(iqw(:,nc),k)    &
                          + zzec(:,nc)*iv(iqe(:,nc),k) + zzsc(:,nc)*iv(iqs(:,nc),k)    &
                          - rhsc(:,k,nc) )*rhelmc(:,k,nc) - iv(iqx(:,nc),k)
      iv(iqx(:,nc),k) = iv(iqx(:,nc),k) + dsol(iqx(:,nc),k)
    end do
    call bounds_colour(iv,nc,klim=klim)
  end do
  
  ! test for convergence
  dsolmax_g(1:klim)=maxval(abs(dsol(1:ifull,1:klim)),dim=1)
  
  ! residual
  do k=1,klim
    w(1:ifull,k)=-izzn*iv(in,k)-izzw*iv(iw,k)-izze*iv(ie,k)-izzs*iv(is,k)+irhs(:,k)+iv(1:ifull,k)*(ihelm(:,k)-izz)
  end do

  ! For when the inital grid cannot be upscaled
  call mgcollectreduce(1,w,dsolmax_g,klim=klim)

  do g=1,min(mg_maxlevel_local,1) ! same as if (mg_maxlevel_local>0) then ...

    ! restriction
    ! (since this always operates within a panel, then ine = ien is always true)
    ng4=mg(1)%ifull_fine
    rhs(1:ng4,1:klim,2)=0.25*(w(mg(1)%fine  ,1:klim)+w(mg(1)%fine_n ,1:klim)  &
                             +w(mg(1)%fine_e,1:klim)+w(mg(1)%fine_ne,1:klim))
                             
    ! merge grids if insufficent points on this processor
    call mgcollectreduce(2,rhs(:,:,2),dsolmax_g,klim=klim)
  
  end do
  
  ! upscale grid
  do g=2,gmax
  
    ng =mg(g)%ifull
    ng4=mg(g)%ifull_fine
                
    ! update scalar field
    ! assume zero for first guess of residual (also avoids additional bounds call)
    !v(:,1:klim,g)=0.
    do k=1,klim
      v(1:ng,k,g)=-rhs(1:ng,k,g)/(helm(1:ng,k,g)-mg(g)%zz)
    end do
    call mgbounds(g,v(:,:,g),klim=klim)

    do k=1,klim
      ! residual
      ws(1:ng)=-mg(g)%zze*v(mg(g)%ie,k,g)-mg(g)%zzw*v(mg(g)%iw,k,g) &
               -mg(g)%zzn*v(mg(g)%in,k,g)-mg(g)%zzs*v(mg(g)%is,k,g)
              !+rhs(1:ng,k,g)+(helm(1:ng,k,g)-mg(g)%zz)*v(1:ng,k,g)

      ! restriction
      ! (calculate coarser grid before mgcollect as more work is done in parallel)
      rhs(1:ng4,k,g+1)=0.25*(ws(mg(g)%fine  )+ws(mg(g)%fine_n )  &
                            +ws(mg(g)%fine_e)+ws(mg(g)%fine_ne))
    end do

    ! merge grids if insufficent points on this processor
    call mgcollectreduce(g+1,rhs(:,:,g+1),dsolmax_g,klim=klim)

  end do

  ! solve coarse grid
  do g=mg_maxlevel,mg_maxlevel_local ! same as if (mg_maxlevel_local==mg_maxlevel) then ...

    ng=mg(g)%ifull

    ! perform LU decomposition and back substitute with RHS
    ! to solve for v on coarse grid
    do k=1,klim
      helm_m=helm_o(:,:,k)
      v(1:ng,k,g)=rhs(1:ng,k,g)
      call mbacksub(helm_m,v(1:ng,k,g),indy(:,k))
    end do
      
  end do
    
  ! downscale grid
  do g=gmax,2,-1

    call mgbcast(g+1,v(:,:,g+1),dsolmax_g,klim=klim)

    ng=mg(g)%ifull
    ng4=mg(g+1)%ifull_coarse

    do k=1,klim
      ! interpolation
      ws(1:ng4)= mg(g+1)%wgt_a*v(mg(g+1)%coarse_a,k,g+1) + mg(g+1)%wgt_bc*v(mg(g+1)%coarse_b,k,g+1) &
              + mg(g+1)%wgt_bc*v(mg(g+1)%coarse_c,k,g+1) +  mg(g+1)%wgt_d*v(mg(g+1)%coarse_d,k,g+1)

      ! extension
      ! No mgbounds as the v halo has already been updated and
      ! the coarse interpolation also updates the w halo
      ws(1:ng4)=v(1:ng4,k,g)+ws(1:ng4)

      ! post smoothing
      v(1:ng,k,g)=(mg(g)%zze*ws(mg(g)%ie)+mg(g)%zzw*ws(mg(g)%iw) &
                  +mg(g)%zzn*ws(mg(g)%in)+mg(g)%zzs*ws(mg(g)%is) &
                  -rhs(1:ng,k,g))/(helm(1:ng,k,g)-mg(g)%zz)
    end do

    call mgbounds(g,v(:,:,g),klim=klim,corner=.true.)

  end do

  do g=1,min(mg_maxlevel_local,1) ! same as if (mg_maxlevel_local>0) then ...
    
    ! fine grid
    call mgbcast(2,v(:,:,2),dsolmax_g,klim=klim)

    ! interpolation
    ng4=mg(2)%ifull_coarse
    do k=1,klim
      w(1:ng4,k)= mg(2)%wgt_a*v(mg(2)%coarse_a,k,2) + mg(2)%wgt_bc*v(mg(2)%coarse_b,k,2) &
               + mg(2)%wgt_bc*v(mg(2)%coarse_c,k,2) +  mg(2)%wgt_d*v(mg(2)%coarse_d,k,2)
    end do
    
  end do

  if (mg(1)%merge_len>1) then
    call mgbcast(1,w,dsolmax_g,klim=klim)
    vdum=0.
    ir=mod(mg(1)%merge_pos-1,mg(1)%merge_row)+1   ! index for proc row
    ic=(mg(1)%merge_pos-1)/mg(1)%merge_row+1      ! index for proc col
    do n=1,npan
      do jj=1,jpan
        iq_a=1+(jj-1)*ipan+(n-1)*ipan*jpan
        iq_b=jj*ipan+(n-1)*ipan*jpan
        iq_c=1+(ir-1)*ipan+(jj-1+(ic-1)*jpan)*ipan*mg(1)%merge_row+(n-1)*ipan*jpan*mg(1)%merge_len
        iq_d=ir*ipan+(jj-1+(ic-1)*jpan)*ipan*mg(1)%merge_row+(n-1)*ipan*jpan*mg(1)%merge_len
        vdum(iq_a:iq_b,1:klim)=w(iq_c:iq_d,1:klim)
      end do
      do i=1,ipan
        iq_a=i+(n-1)*ipan*jpan
        iq_c=i+(ir-1)*ipan+((ic-1)*jpan)*ipan*mg(1)%merge_row+(n-1)*ipan*jpan*mg(1)%merge_len
        iq_b=is(iq_a)
        iq_d=mg(1)%is(iq_c)
        vdum(iq_b,1:klim)=w(iq_d,1:klim)
        iq_a=i+(jpan-1)*ipan+(n-1)*ipan*jpan
        iq_c=i+(ir-1)*ipan+(jpan-1+(ic-1)*jpan)*ipan*mg(1)%merge_row+(n-1)*ipan*jpan*mg(1)%merge_len
        iq_b=in(iq_a)
        iq_d=mg(1)%in(iq_c)
        vdum(iq_b,1:klim)=w(iq_d,1:klim)
      end do  
      do j=1,jpan
        iq_a=1+(j-1)*ipan+(n-1)*ipan*jpan
        iq_c=1+(ir-1)*ipan+(j-1+(ic-1)*jpan)*ipan*mg(1)%merge_row+(n-1)*ipan*jpan*mg(1)%merge_len
        iq_b=iw(iq_a)
        iq_d=mg(1)%iw(iq_c)
        vdum(iq_b,1:klim)=w(iq_d,1:klim)
        iq_a=ipan+(j-1)*ipan+(n-1)*ipan*jpan
        iq_c=ipan+(ir-1)*ipan+(j-1+(ic-1)*jpan)*ipan*mg(1)%merge_row+(n-1)*ipan*jpan*mg(1)%merge_len
        iq_b=ie(iq_a)
        iq_d=mg(1)%ie(iq_c)
        vdum(iq_b,1:klim)=w(iq_d,1:klim)
      end do
    end do
    w(1:ifull+iextra,1:klim)=vdum(1:ifull+iextra,1:klim)    
  else
    ! remap mg halo to normal halo
    vdum(ifull+1:ifull+iextra,1:klim)=0.
    do n=1,npan
      do i=1,ipan
        iq=indx(i,1,n-1,ipan,jpan)
        iq_a=is(iq)
        iq_b=mg(1)%is(iq)
        vdum(iq_a,1:klim)=w(iq_b,1:klim)
        iq=indx(i,jpan,n-1,ipan,jpan)
        iq_a=in(iq)
        iq_b=mg(1)%in(iq)
        vdum(iq_a,1:klim)=w(iq_b,1:klim)
      end do  
      do j=1,jpan
        iq=indx(1,j,n-1,ipan,jpan)
        iq_a=iw(iq)
        iq_b=mg(1)%iw(iq)
        vdum(iq_a,1:klim)=w(iq_b,1:klim)
        iq=indx(ipan,j,n-1,ipan,jpan)
        iq_a=ie(iq)
        iq_b=mg(1)%ie(iq)
        vdum(iq_a,1:klim)=w(iq_b,1:klim)
      end do
    end do
    w(ifull+1:ifull+iextra,1:klim)=vdum(ifull+1:ifull+iextra,1:klim)
  end if

  ! extension
  iv(1:ifull+iextra,1:klim)=iv(1:ifull+iextra,1:klim)+w(1:ifull+iextra,1:klim)
  
  ! post smoothing
  do nc=1,maxcolour
    do k=1,klim
      dsol(iqx(:,nc),k)=( zznc(:,nc)*iv(iqn(:,nc),k) + zzwc(:,nc)*iv(iqw(:,nc),k)    &
                        + zzec(:,nc)*iv(iqe(:,nc),k) + zzsc(:,nc)*iv(iqs(:,nc),k)    &
                        - rhsc(:,k,nc))*rhelmc(:,k,nc) - iv(iqx(:,nc),k)
      iv(iqx(:,nc),k) = iv(iqx(:,nc),k) + dsol(iqx(:,nc),k)
    end do
    call bounds_colour(iv,nc,klim=klim)
  end do

  ! test for convergence
  knew=klim
  do k=klim,1,-1
    iters(k)=itr
    if (dsolmax_g(k)>=restol*(smaxmin_g(1,k)-smaxmin_g(2,k))) exit
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
    !write(6,*) "itc ",itc
  end if
end if

call end_log(helm_end)

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
integer itr,itrc,g,ng,ng4,n,i,j,ir,ic,jj,iq
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
real, dimension(mg_maxsize) :: au,bu,cu
real, dimension(mg_maxsize,2,gmax+1) :: v
real, dimension(mg_maxsize,10,gmax+1) :: yy
real, dimension(mg_maxsize,8) :: w
real, dimension(mg_maxsize,gmax+1) :: zz,zzn,zzs,zze,zzw
real, dimension(mg_maxsize,gmax+1) :: hh
real, dimension(mg_maxsize,gmax+1) :: rhs
real, dimension(mg_maxsize,gmax+1) :: rhsice
real, dimension(mg_maxsize,2) :: dsol
real, dimension(mg_maxsize) :: ws,new
real, dimension(ifull+iextra,2) :: dumc
real, dimension(mg_minsize,mg_minsize) :: helm_m,helm_o
real, dimension(mg_ifullc,3) :: yyzcu,yyncu,yyscu,yyecu,yywcu
real, dimension(mg_ifullc,3) :: zzhhcu,zzncu,zzscu,zzecu,zzwcu,rhscu
real, dimension(1) :: dsolmax
real, dimension(8) :: dsolmax_g

if (sorfirst) then
  write(6,*) "ERROR: mgsormlo requires mgsor_init to be called first"
  call ccmpi_abort(-1)
end if

! The following expressions describe the residual terms

! yy*neta*(DIV^2 neta) + zz*(DIV^2 neta) + hh*neta = rhs
! neta = n0 + e
! neta is the solution, n0 is first guess and e is the residual term

! yy*n0*(DIV^2 n0) + zz*(DIV^2 n0) + hh*n0 = rhs - E
! yy*e*(DIV^2 e) + (yy*n0+zz)*(DIV^2 e) + (yy*(DIV^2 n0)+hh)*e  = E

! so yy is simply upscaled
! zz -> zz + yy*n0
! hh -> hh + yy*(DIV^2 n0)

! also

! zz*(DIV^2 ipice) + yy*ipice = rhs
! ipice = i0 + f

! zz*(DIV^2 i0) + yy*i0 = rhs - F
! zz*(DIV^2 f) + yy*f = F

! upscale coeffs
yy(1:ifull,1,1)=iyy(1:ifull)
yy(1:ifull,2,1)=iyyn(1:ifull)
yy(1:ifull,3,1)=iyys(1:ifull)
yy(1:ifull,4,1)=iyye(1:ifull)
yy(1:ifull,5,1)=iyyw(1:ifull)
yy(1:ifull,6,1)=izz(1:ifull,2)
yy(1:ifull,7,1)=izzn(1:ifull,2)
yy(1:ifull,8,1)=izzs(1:ifull,2)
yy(1:ifull,9,1)=izze(1:ifull,2)
yy(1:ifull,10,1)=izzw(1:ifull,2)
call mgcollect(1,yy(:,:,1))
do g=1,gmax
  ng4=mg(g)%ifull_fine
  yy(1:ng4,1:5,g+1)=0.25*dfac*(yy(mg(g)%fine  ,1:5,g)+yy(mg(g)%fine_n ,1:5,g) &
                              +yy(mg(g)%fine_e,1:5,g)+yy(mg(g)%fine_ne,1:5,g))
  ! special treatment of cavitating fluid (no dfac)
  yy(1:ng4,6:10,g+1)=0.25*(yy(mg(g)%fine  ,6:10,g)+yy(mg(g)%fine_n ,6:10,g)   &
                          +yy(mg(g)%fine_e,6:10,g)+yy(mg(g)%fine_ne,6:10,g))
  call mgcollect(g+1,yy(:,:,g+1))
end do
do g=mg_maxlevel,mg_maxlevel_local ! same as if (mg_maxlevel_local==mg_maxlevel) then ...
  helm_o=0.
  ng=mg(g)%ifull
  do iq=1,ng
    helm_o(iq,iq)=yy(iq,6,g)      
    helm_o(mg(g)%in(iq),iq)=yy(iq,7,g)
    helm_o(mg(g)%is(iq),iq)=yy(iq,8,g)
    helm_o(mg(g)%ie(iq),iq)=yy(iq,9,g)
    helm_o(mg(g)%iw(iq),iq)=yy(iq,10,g)
  end do
  call mdecomp(helm_o,indy) ! destroys helm_m
end do

! solver requires bounds to be updated
dumc(1:ifull,1)=neta(1:ifull)
dumc(1:ifull,2)=ipice(1:ifull)
call bounds(dumc)
neta(ifull+1:ifull+iextra) =dumc(ifull+1:ifull+iextra,1)
ipice(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,2)

dsolmax_g=0.

! Main loop
do itr=1,itr_mgice

  do nc=1,maxcolour
  
    ! ocean
    au(1:ifullx)=iyy(iqx(:,nc))
    bu(1:ifullx)=izz(iqx(:,nc),1)+ihh(iqx(:,nc))                                     &
                +iyyn(iqx(:,nc))*neta(iqn(:,nc))+iyys(iqx(:,nc))*neta(iqs(:,nc))     &
                +iyye(iqx(:,nc))*neta(iqe(:,nc))+iyyw(iqx(:,nc))*neta(iqw(:,nc))
    cu(1:ifullx)=izzn(iqx(:,nc),1)*neta(iqn(:,nc))+izzs(iqx(:,nc),1)*neta(iqs(:,nc)) &
                +izze(iqx(:,nc),1)*neta(iqe(:,nc))+izzw(iqx(:,nc),1)*neta(iqw(:,nc)) &
                -irhs(iqx(:,nc),1)        
    new(1:ifullx) = -2.*cu(1:ifullx)/(bu(1:ifullx)+sqrt(bu(1:ifullx)*bu(1:ifullx)-4.*au(1:ifullx)*cu(1:ifullx)))
    new(1:ifullx) = max(new(1:ifullx),-dd(iqx(:,nc)))*ee(iqx(:,nc))
    dsol(iqx(:,nc),1)=new(1:ifullx)-neta(iqx(:,nc))
    neta(iqx(:,nc)) =new(1:ifullx)
    
    ! ice
    new(1:ifullx) = ( -izzn(iqx(:,nc),2)*ipice(iqn(:,nc)) &
                      -izzs(iqx(:,nc),2)*ipice(iqs(:,nc)) &
                      -izze(iqx(:,nc),2)*ipice(iqe(:,nc)) &
                      -izzw(iqx(:,nc),2)*ipice(iqw(:,nc)) &
                     + irhs(iqx(:,nc),2) ) / izz(iqx(:,nc),2)
    ! ice cavitating fluid
    new(1:ifullx) = max(min(new(1:ifullx),ipmax(iqx(:,nc))),0.)
    dsol(iqx(:,nc),2)=new(1:ifullx)-ipice(iqx(:,nc))   
    ipice(iqx(:,nc))=new(1:ifullx)

    dumc(1:ifull,1)=neta(1:ifull)
    dumc(1:ifull,2)=ipice(1:ifull)
    call bounds_colour(dumc,nc)
    neta(ifull+1:ifull+iextra) =dumc(ifull+1:ifull+iextra,1)
    ipice(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,2)
  end do

  ! test for convergence
  dsolmax_g(1:2)=maxval(abs(dsol(1:ifull,1:2)),dim=1)

  w(1:ifull,2)= izz(:,1)+ iyy*neta(1:ifull)
  w(1:ifull,3)=izzn(:,1)+iyyn*neta(1:ifull)
  w(1:ifull,4)=izzs(:,1)+iyys*neta(1:ifull)
  w(1:ifull,5)=izze(:,1)+iyye*neta(1:ifull)
  w(1:ifull,6)=izzw(:,1)+iyyw*neta(1:ifull)
  w(1:ifull,7)=ihh+iyy*neta(1:ifull)+iyyn*neta(in)+iyys*neta(is)+iyye*neta(ie)+iyyw*neta(iw)

  ! residual
  w(1:ifull,1)=-neta(1:ifull)*(     iyy*neta(1:ifull)     +iyyn*neta(in)     +iyys*neta(is)     +iyye*neta(ie)     +iyyw*neta(iw)) &
                             -(izz(:,1)*neta(1:ifull)+izzn(:,1)*neta(in)+izzs(:,1)*neta(is)+izze(:,1)*neta(ie)+izzw(:,1)*neta(iw)) &
                             -ihh*neta(1:ifull)+irhs(:,1)
  w(1:ifull,8)=-(izz(:,2)*ipice(1:ifull)+izzn(:,2)*ipice(in)+izzs(:,2)*ipice(is) &
                +izze(:,2)*ipice(ie)+izzw(:,2)*ipice(iw))+irhs(:,2)
                
  ! patch to remove error when ipmax is reached
  where (ipice(1:ifull)>=ipmax(1:ifull))
    w(1:ifull,8)=0.
  end where
                             
  w(1:ifull,1)=w(1:ifull,1)*ee(1:ifull)
  w(1:ifull,8)=w(1:ifull,8)*ee(1:ifull)
  
  ! For when the inital grid cannot be upscaled
  call mgcollectreduce(1,w(:,:),dsolmax_g)
  
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
    hh(1:ng4,2)=0.25*(w(mg(1)%fine  ,7)+w(mg(1)%fine_n ,7) &
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
      call mgcollectreduce(2,w(:,:),dsolmax_g)
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

  ! upscale grid
  do g=2,gmax
  
    ng=mg(g)%ifull

    ! update
    ! possibly use colours here, although v is reset to zero every iteration
    ! assume zero for first guess of residual (also avoids additional bounds call)
    au(1:ng)=yy(1:ng,1,g)
    bu(1:ng)=zz(1:ng,g)+hh(1:ng,g)
    cu(1:ng)=-rhs(1:ng,g)
    v(1:ng,1,g) = -2.*cu(1:ng)/(bu(1:ng)+sqrt(bu(1:ng)*bu(1:ng)-4.*au(1:ng)*cu(1:ng)))

    v(1:ng,2,g) = rhsice(1:ng,g)/yy(1:ng,6,g)
    
    ! residual
    call mgbounds(g,v(:,:,g))

    ! restriction
    ! (calculate finer grid before mgcollect as the messages sent/recv are shorter)
    ng4=mg(g)%ifull_fine
    ws(1:ng)= zz(1:ng,g)+yy(1:ng,1,g)*v(1:ng,1,g)
    zz(1:ng4,g+1)=0.25*dfac*(ws(mg(g)%fine  )+ws(mg(g)%fine_n ) &
                            +ws(mg(g)%fine_e)+ws(mg(g)%fine_ne))
    ws(1:ng)=zzn(1:ng,g)+yy(1:ng,2,g)*v(1:ng,1,g)
    zzn(1:ng4,g+1)=0.25*dfac*(ws(mg(g)%fine  )+ws(mg(g)%fine_n ) &
                             +ws(mg(g)%fine_e)+ws(mg(g)%fine_ne))
    ws(1:ng)=zzs(1:ng,g)+yy(1:ng,3,g)*v(1:ng,1,g)
    zzs(1:ng4,g+1)=0.25*dfac*(ws(mg(g)%fine  )+ws(mg(g)%fine_n ) &
                             +ws(mg(g)%fine_e)+ws(mg(g)%fine_ne))
    ws(1:ng)=zze(1:ng,g)+yy(1:ng,4,g)*v(1:ng,1,g)
    zze(1:ng4,g+1)=0.25*dfac*(ws(mg(g)%fine  )+ws(mg(g)%fine_n ) &
                             +ws(mg(g)%fine_e)+ws(mg(g)%fine_ne))
    ws(1:ng)=zzw(1:ng,g)+yy(1:ng,5,g)*v(1:ng,1,g)
    zzw(1:ng4,g+1)=0.25*dfac*(ws(mg(g)%fine  )+ws(mg(g)%fine_n ) &
                             +ws(mg(g)%fine_e)+ws(mg(g)%fine_ne))
    ws(1:ng)=hh(1:ng,g)+yy(1:ng,1,g)*v(1:ng,1,g)+yy(1:ng,2,g)*v(mg(g)%in,1,g)+yy(1:ng,3,g)*v(mg(g)%is,1,g) &
                                                +yy(1:ng,4,g)*v(mg(g)%ie,1,g)+yy(1:ng,5,g)*v(mg(g)%iw,1,g)
    hh(1:ng4,g+1)=0.25*(ws(mg(g)%fine  )+ws(mg(g)%fine_n ) &
                       +ws(mg(g)%fine_e)+ws(mg(g)%fine_ne))

    ! ocean
    ws(1:ng)=-v(1:ng,1,g)*(yy(1:ng,1,g)*v(1:ng,1,g)+yy(1:ng,2,g)*v(mg(g)%in,1,g)+yy(1:ng,3,g)*v(mg(g)%is,1,g)+yy(1:ng,4,g)*v(mg(g)%ie,1,g)+yy(1:ng,5,g)*v(mg(g)%iw,1,g)) &
                           -(zz(1:ng,g)*v(1:ng,1,g)+ zzn(1:ng,g)*v(mg(g)%in,1,g)+ zzs(1:ng,g)*v(mg(g)%is,1,g)+ zze(1:ng,g)*v(mg(g)%ie,1,g)+ zzw(1:ng,g)*v(mg(g)%iw,1,g)) &
                            -hh(1:ng,g)*v(1:ng,1,g)+rhs(1:ng,g)
    rhs(1:ng4,g+1)=0.25*(ws(mg(g)%fine  )+ws(mg(g)%fine_n ) &
                        +ws(mg(g)%fine_e)+ws(mg(g)%fine_ne))
                          
    ! ice
    ws(1:ng)=-(yy(1:ng,7,g)*v(mg(g)%in,2,g)+ yy(1:ng,8,g)*v(mg(g)%is,2,g)   &
              +yy(1:ng,9,g)*v(mg(g)%ie,2,g)+yy(1:ng,10,g)*v(mg(g)%iw,2,g))
             !-yy(1:ng,6,g)*v(1:ng,2,g)+rhsice(1:ng,g)
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
      call mgcollectreduce(g+1,w(:,:),dsolmax_g)
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

  ! solve coarse grid    
  do g=mg_maxlevel,mg_maxlevel_local ! same as if (mg_maxlevel==mg_maxlevel_local) then ...

    ng=mg(g)%ifull
      
    ! solve for ice using LU decomposition and back substitution with RHS
    helm_m=helm_o
    v(1:ng,2,g)=rhsice(1:ng,g)
    call mbacksub(helm_m,v(1:ng,2,g),indy)

    ! solve non-linear water free surface with coloured SOR
    v(1:ng,1,g)=0.
  
    ! pack yy,zz,hh and rhs by colour
    do nc=1,3
      yyzcu(1:mg_ifullc,nc)=yy(col_iq(:,nc),1,g)
      yyncu(1:mg_ifullc,nc)=yy(col_iq(:,nc),2,g)
      yyscu(1:mg_ifullc,nc)=yy(col_iq(:,nc),3,g)
      yyecu(1:mg_ifullc,nc)=yy(col_iq(:,nc),4,g)
      yywcu(1:mg_ifullc,nc)=yy(col_iq(:,nc),5,g)
      zzhhcu(1:mg_ifullc,nc)=zz(col_iq(:,nc),g)+hh(col_iq(:,nc),g)
      zzncu(1:mg_ifullc,nc)=zzn(col_iq(:,nc),g)
      zzscu(1:mg_ifullc,nc)=zzs(col_iq(:,nc),g)
      zzecu(1:mg_ifullc,nc)=zze(col_iq(:,nc),g)
      zzwcu(1:mg_ifullc,nc)=zzw(col_iq(:,nc),g)
      rhscu(1:mg_ifullc,nc)=rhs(col_iq(:,nc),g)
    end do
  
    do itrc=1,itr_mgice

      do nc=1,3
        au(1:mg_ifullc)=yyzcu(:,nc)
        bu(1:mg_ifullc)=zzhhcu(:,nc)+yyncu(:,nc)*v(col_iqn(:,nc),1,g)+yyscu(:,nc)*v(col_iqs(:,nc),1,g) &
                                    +yyecu(:,nc)*v(col_iqe(:,nc),1,g)+yywcu(:,nc)*v(col_iqw(:,nc),1,g)
        cu(1:mg_ifullc)=zzncu(:,nc)*v(col_iqn(:,nc),1,g)+zzscu(:,nc)*v(col_iqs(:,nc),1,g)   &
                       +zzecu(:,nc)*v(col_iqe(:,nc),1,g)+zzwcu(:,nc)*v(col_iqw(:,nc),1,g)   &
                       -rhscu(:,nc)
        new(1:mg_ifullc) = -2.*cu(1:mg_ifullc)/(bu(1:mg_ifullc)+sqrt(bu(1:mg_ifullc)*bu(1:mg_ifullc)-4.*au(1:mg_ifullc)*cu(1:mg_ifullc)))
     
        dsol(col_iq(:,nc),1)=new(1:mg_ifullc)-v(col_iq(:,nc),1,g)
        v(col_iq(:,nc),1,g)=new(1:mg_ifullc)
      end do
    
      dsolmax(1)=maxval(abs(dsol(1:ng,1)))
      if (dsolmax(1)<tol) exit

    end do
  
  end do
    
  ! downscale grid
  do g=gmax,2,-1

    call mgbcast(g+1,v(:,:,g+1),dsolmax_g(1:2))

    ! ocean
    ! interpolation
    ng4=mg(g+1)%ifull_coarse
    ws(1:ng4)= mg(g+1)%wgt_a*v(mg(g+1)%coarse_a,1,g+1)  + mg(g+1)%wgt_bc*v(mg(g+1)%coarse_b,1,g+1) &
             + mg(g+1)%wgt_bc*v(mg(g+1)%coarse_c,1,g+1) + mg(g+1)%wgt_d*v(mg(g+1)%coarse_d,1,g+1)
    ! extension
    ! No mgbounds as the v halo has already been updated and
    ! the coarse interpolation also updates the w halo
    ws(1:ng4)=v(1:ng4,1,g)+ws(1:ng4)

    ! post smoothing
    ng=mg(g)%ifull
    au(1:ng)=yy(1:ng,1,g)
    bu(1:ng)=zz(1:ng,g)+hh(1:ng,g)+yy(1:ng,2,g)*ws(mg(g)%in)+yy(1:ng,3,g)*ws(mg(g)%is)+yy(1:ng,4,g)*ws(mg(g)%ie)+yy(1:ng,5,g)*ws(mg(g)%iw)
    cu(1:ng)=zzn(1:ng,g)*ws(mg(g)%in)+zzs(1:ng,g)*ws(mg(g)%is)+zze(1:ng,g)*ws(mg(g)%ie)+zzw(1:ng,g)*ws(mg(g)%iw)-rhs(1:ng,g)
    v(1:ng,1,g) = -2.*cu(1:ng)/(bu(1:ng)+sqrt(bu(1:ng)*bu(1:ng)-4.*au(1:ng)*cu(1:ng)))

    ! ice
    ! interpolation
    ws(1:ng4)= mg(g+1)%wgt_a*v(mg(g+1)%coarse_a,2,g+1)  + mg(g+1)%wgt_bc*v(mg(g+1)%coarse_b,2,g+1) &
             + mg(g+1)%wgt_bc*v(mg(g+1)%coarse_c,2,g+1) + mg(g+1)%wgt_d*v(mg(g+1)%coarse_d,2,g+1)

    ! extension
    ! No mgbounds as the v halo has already been updated and
    ! the coarse interpolation also updates the w halo
    ws(1:ng4)=v(1:ng4,2,g)+ws(1:ng4)
    
    v(1:ng,2,g) = ( -yy(1:ng,7,g)*ws(mg(g)%in)- yy(1:ng,8,g)*ws(mg(g)%is) &
                    -yy(1:ng,9,g)*ws(mg(g)%ie)-yy(1:ng,10,g)*ws(mg(g)%iw) &
                    +rhsice(1:ng,g) ) / yy(1:ng,6,g)

    call mgbounds(g,v(:,:,g),corner=.true.)

  end do

  do g=1,min(mg_maxlevel_local,1) ! same as if (mg_maxlevel_local>0) then ...
    
    ! fine grid
    call mgbcast(2,v(:,:,2),dsolmax_g(1:2))

    ! interpolation
    ng4=mg(2)%ifull_coarse
    w(1:ng4,1)= mg(2)%wgt_a*v(mg(2)%coarse_a,1,2) + mg(2)%wgt_bc*v(mg(2)%coarse_b,1,2) &
              + mg(2)%wgt_bc*v(mg(2)%coarse_c,1,2) +  mg(2)%wgt_d*v(mg(2)%coarse_d,1,2)
    w(1:ng4,8)= mg(2)%wgt_a*v(mg(2)%coarse_a,2,2) + mg(2)%wgt_bc*v(mg(2)%coarse_b,2,2) &
              + mg(2)%wgt_bc*v(mg(2)%coarse_c,2,2) +  mg(2)%wgt_d*v(mg(2)%coarse_d,2,2)

  end do

  if (mg(1)%merge_len>1) then
    call mgbcast(1,w(:,1:1),dsolmax_g(1:1))
    call mgbcast(1,w(:,8:8),dsolmax_g(2:2))
    vduma=0.
    vdumb=0.
    ir=mod(mg(1)%merge_pos-1,mg(1)%merge_row)+1   ! index for proc row
    ic=(mg(1)%merge_pos-1)/mg(1)%merge_row+1      ! index for proc col
    do n=1,npan
      do jj=1,jpan
        iq_a=1+(jj-1)*ipan+(n-1)*ipan*jpan
        iq_b=jj*ipan+(n-1)*ipan*jpan
        iq_c=1+(ir-1)*ipan+(jj-1+(ic-1)*jpan)*ipan*mg(1)%merge_row+(n-1)*ipan*jpan*mg(1)%merge_len
        iq_d=ir*ipan+(jj-1+(ic-1)*jpan)*ipan*mg(1)%merge_row+(n-1)*ipan*jpan*mg(1)%merge_len
        vduma(iq_a:iq_b)=w(iq_c:iq_d,1)
        vdumb(iq_a:iq_b)=w(iq_c:iq_d,8)
      end do
      do i=1,ipan
        iq_a=i+(n-1)*ipan*jpan
        iq_c=i+(ir-1)*ipan+(ic-1)*jpan*ipan*mg(1)%merge_row+(n-1)*ipan*jpan*mg(1)%merge_len
        iq_b=is(iq_a)
        iq_d=mg(1)%is(iq_c)
        vduma(iq_b)=w(iq_d,1)
        vdumb(iq_b)=w(iq_d,8)
        iq_a=i+(jpan-1)*ipan+(n-1)*ipan*jpan
        iq_c=i+(ir-1)*ipan+(jpan-1+(ic-1)*jpan)*ipan*mg(1)%merge_row+(n-1)*ipan*jpan*mg(1)%merge_len
        iq_b=in(iq_a)
        iq_d=mg(1)%in(iq_c)
        vduma(iq_b)=w(iq_d,1)
        vdumb(iq_b)=w(iq_d,8)
      end do  
      do j=1,jpan
        iq_a=1+(j-1)*ipan+(n-1)*ipan*jpan
        iq_c=1+(ir-1)*ipan+(j-1+(ic-1)*jpan)*ipan*mg(1)%merge_row+(n-1)*ipan*jpan*mg(1)%merge_len
        iq_b=iw(iq_a)
        iq_d=mg(1)%iw(iq_c)
        vduma(iq_b)=w(iq_d,1)
        vdumb(iq_b)=w(iq_d,8)
        iq_a=ipan+(j-1)*ipan+(n-1)*ipan*jpan
        iq_c=ipan+(ir-1)*ipan+(j-1+(ic-1)*jpan)*ipan*mg(1)%merge_row+(n-1)*ipan*jpan*mg(1)%merge_len
        iq_b=ie(iq_a)
        iq_d=mg(1)%ie(iq_c)
        vduma(iq_b)=w(iq_d,1)
        vdumb(iq_b)=w(iq_d,8)
      end do
    end do
    w(1:ifull+iextra,1)=vduma(1:ifull+iextra)
    w(1:ifull+iextra,8)=vdumb(1:ifull+iextra)
  else
    ! remap mg halo to normal halo 
    vduma(ifull+1:ifull+iextra)=0.
    vdumb(ifull+1:ifull+iextra)=0.
    do n=1,npan
      do i=1,ipan
        iq=indx(i,1,n-1,ipan,jpan)
        iq_a=is(iq)
        iq_b=mg(1)%is(iq)
        vduma(iq_a)=w(iq_b,1)
        vdumb(iq_a)=w(iq_b,8)
        iq=indx(i,jpan,n-1,ipan,jpan)
        iq_a=in(iq)
        iq_b=mg(1)%in(iq)
        vduma(iq_a)=w(iq_b,1)
        vdumb(iq_a)=w(iq_b,8)
      end do  
      do j=1,jpan
        iq=indx(1,j,n-1,ipan,jpan)
        iq_a=iw(iq)
        iq_b=mg(1)%iw(iq)
        vduma(iq_a)=w(iq_b,1)
        vdumb(iq_a)=w(iq_b,8)
        iq=indx(ipan,j,n-1,ipan,jpan)
        iq_a=ie(iq)
        iq_b=mg(1)%ie(iq)
        vduma(iq_a)=w(iq_b,1)
        vdumb(iq_a)=w(iq_b,8)
      end do
    end do
    w(ifull+1:ifull+iextra,1)=vduma(ifull+1:ifull+iextra)
    w(ifull+1:ifull+iextra,8)=vdumb(ifull+1:ifull+iextra)
  end if
 
  ! extension
  neta(1:ifull+iextra)=neta(1:ifull+iextra)+w(1:ifull+iextra,1)
  ipice(1:ifull+iextra)=ipice(1:ifull+iextra)+w(1:ifull+iextra,8)
 
  neta=max(neta,-dd)*ee
  ipice=max(min(ipice,ipmax),0.) 
  
  ! post smoothing
  
  do nc=1,maxcolour
    
    ! ocean
    au(1:ifullx)=iyy(iqx(:,nc))
    bu(1:ifullx)=izz(iqx(:,nc),1)+ihh(iqx(:,nc))                                     &
                +iyyn(iqx(:,nc))*neta(iqn(:,nc))+iyys(iqx(:,nc))*neta(iqs(:,nc))     &
                +iyye(iqx(:,nc))*neta(iqe(:,nc))+iyyw(iqx(:,nc))*neta(iqw(:,nc))
    cu(1:ifullx)=izzn(iqx(:,nc),1)*neta(iqn(:,nc))+izzs(iqx(:,nc),1)*neta(iqs(:,nc)) &
                +izze(iqx(:,nc),1)*neta(iqe(:,nc))+izzw(iqx(:,nc),1)*neta(iqw(:,nc)) &
                -irhs(iqx(:,nc),1)        
    new(1:ifullx) = -2.*cu(1:ifullx)/(bu(1:ifullx)+sqrt(bu(1:ifullx)*bu(1:ifullx)-4.*au(1:ifullx)*cu(1:ifullx)))
    neta(iqx(:,nc))=max(new(1:ifullx),-dd(iqx(:,nc)))*ee(iqx(:,nc))
    
    ! ice
    new(1:ifullx) = ( -izzn(iqx(:,nc),2)*ipice(iqn(:,nc)) &
                      -izzs(iqx(:,nc),2)*ipice(iqs(:,nc)) &
                      -izze(iqx(:,nc),2)*ipice(iqe(:,nc)) &
                      -izzw(iqx(:,nc),2)*ipice(iqw(:,nc)) &
                     + irhs(iqx(:,nc),2) ) / izz(iqx(:,nc),2)
    ! cavitating fluid
    ipice(iqx(:,nc))=max(min(new(1:ifullx),ipmax(iqx(:,nc))),0.)

    dumc(1:ifull,1)=neta(1:ifull)
    dumc(1:ifull,2)=ipice(1:ifull)
    call bounds_colour(dumc(:,:),nc)
    neta(ifull+1:ifull+iextra) =dumc(ifull+1:ifull+iextra,1)
    ipice(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,2)
  end do
  
  ! test for convergence
  if (dsolmax_g(1)<tol.and.dsolmax_g(2)<itol) exit
  
end do

totits     =itr
maxglobseta=dsolmax_g(1)
maxglobip  =dsolmax_g(2)

return
end subroutine mgmlo

! LU decomposition
subroutine mdecomp(a,indy)

implicit none

real, dimension(:,:), intent(inout) :: a
real, dimension(size(a,1)) :: vv,dumv
real aamax,sumx,dum
integer, dimension(size(a,1)), intent(out) :: indy
integer n,i,j,k
integer imax

n=size(a,1)

do i=1,n
  aamax=maxval(abs(a(i,:)))
  vv(i)=1./aamax
end do
do j=1,n-1
  do i=1,j-1
    sumx=a(i,j)-sum(a(i,1:i-1)*a(1:i-1,j))
    a(i,j)=sumx
  end do
  aamax=0.
  do i=j,n
    sumx=a(i,j)-sum(a(i,1:j-1)*a(1:j-1,j))
    a(i,j)=sumx
    dum=vv(i)*abs(sumx)
    if (dum>=aamax) then
      imax=i
      aamax=dum
    end if
  end do
  if (j/=imax) then
    dumv(:)=a(imax,:)
    a(imax,:)=a(j,:)
    a(j,:)=dumv(:)
    vv(imax)=vv(j)
  end if
  indy(j)=imax
  dum=1./a(j,j)
  a(j+1:n,j)=a(j+1:n,j)*dum
end do
!j=n
do i=1,n-1
  sumx=a(i,n)-sum(a(i,1:i-1)*a(1:i-1,n))
  a(i,n)=sumx
end do
aamax=0.
!i=n
sumx=a(n,n)-sum(a(n,1:n-1)*a(1:n-1,n))
a(n,n)=sumx
dum=vv(n)*abs(sumx)
if (dum>=aamax) then
  imax=n
  !aamax=dum
end if
indy(n)=imax

return
end subroutine mdecomp

! Back substitution
subroutine mbacksub(a,b,indy)

implicit none

real, dimension(:,:), intent(in) :: a
real, dimension(size(a,1)), intent(inout) :: b
real sumx
integer, dimension(size(a,1)), intent(in) :: indy
integer n,i,ii,ll

n=size(a,1)
ii=0
do i=1,n
  ll=indy(i)
  sumx=b(ll)
  b(ll)=b(i)
  if (ii/=0) then
    sumx=sumx-sum(a(i,ii:i-1)*b(ii:i-1))
  else if (sumx/=0.) then
    ii=i
  end if
  b(i)=sumx
end do
do i=n,1,-1
  sumx=b(i)-sum(a(i,i+1:n)*b(i+1:n))
  b(i)=sumx/a(i,i)
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
mg_ifullc=0
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
    if (rank/=0) then
      mg_maxlevel_local=0
      mg(1)%nmax=0
    end if
  
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
      if (rank/=0) then
        mg_maxlevel_local=min(g-1,mg_maxlevel_local)
        mg(g)%nmax=0
      end if
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
real, dimension(mg_maxsize,5) :: dum

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

dum(1:ifull,1)=zzn
dum(1:ifull,2)=zzs
dum(1:ifull,3)=zze
dum(1:ifull,4)=zzw
dum(1:ifull,5)=zz
call mgcollect(1,dum)
mg(1)%zzn(:)=dum(1:np,1)
mg(1)%zzs(:)=dum(1:np,2)
mg(1)%zze(:)=dum(1:np,3)
mg(1)%zzw(:)=dum(1:np,4)
mg(1)%zz(:)=dum(1:np,5)

do g=2,gmax+1
  np=mg(g)%ifull
  do iq=1,mg(g-1)%ifull_fine
    iqq=mg(g-1)%fine(iq)
    dum(iq,1)=0.25*dfac*(mg(g-1)%zzn(iqq)+mg(g-1)%zzn(mg(g-1)%in(iqq))+mg(g-1)%zzn(mg(g-1)%ine(iqq))+mg(g-1)%zzn(mg(g-1)%ie(iqq)))
    dum(iq,2)=0.25*dfac*(mg(g-1)%zzs(iqq)+mg(g-1)%zzs(mg(g-1)%in(iqq))+mg(g-1)%zzs(mg(g-1)%ine(iqq))+mg(g-1)%zzs(mg(g-1)%ie(iqq)))
    dum(iq,3)=0.25*dfac*(mg(g-1)%zze(iqq)+mg(g-1)%zze(mg(g-1)%in(iqq))+mg(g-1)%zze(mg(g-1)%ine(iqq))+mg(g-1)%zze(mg(g-1)%ie(iqq)))
    dum(iq,4)=0.25*dfac*(mg(g-1)%zzw(iqq)+mg(g-1)%zzw(mg(g-1)%in(iqq))+mg(g-1)%zzw(mg(g-1)%ine(iqq))+mg(g-1)%zzw(mg(g-1)%ie(iqq)))
    dum(iq,5)=0.25*dfac*(mg(g-1)%zz(iqq) +mg(g-1)%zz(mg(g-1)%in(iqq)) +mg(g-1)%zz(mg(g-1)%ine(iqq)) +mg(g-1)%zz(mg(g-1)%ie(iqq)))
  end do
  call mgcollect(g,dum)
  mg(g)%zzn=dum(1:np,1)
  mg(g)%zzs=dum(1:np,2)
  mg(g)%zze=dum(1:np,3)
  mg(g)%zzw=dum(1:np,4)
  mg(g)%zz=dum(1:np,5)
end do

if (myid==0) then
  write(6,*) "Finished initialising atmosphere multi-grid coupling arrays"
end if

return
end subroutine mgzz_init

end module mgsolve
