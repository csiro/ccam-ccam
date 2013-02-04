module mgsolve

! This module solves for the atmosphere helmholtz equation and the ocean free
! surface equation using a multi-grid approach.

! Currently we simply aggregate data as we upscale the grid, untill eventially
! all processors solve for the same coarse grid (thereby avoiding MPI_Allreduce)
! when solving the coarse grid.

implicit none

private
public mghelm,mgmlo,mgsor_init,mgzz_init

integer, save :: mg_maxsize,mg_minsize
logical, save :: sorfirst=.true.
logical, save :: zzfirst=.true.
logical, save :: mlofirst=.true.

integer, parameter :: itr_max=300 ! maximum number of iterations
integer, parameter :: itr_mg=30
integer, parameter :: itr_mgice=30

contains

! This version is for the atmosphere
subroutine mghelm(izz,izzn,izze,izzw,izzs,ihelm,iv,jrhs)

use cc_mpi
use indices_m

implicit none

include 'newmpar.h'
include 'parm.h'
include 'parmdyn.h'

integer, dimension(mg_maxsize) :: iqaa,iqbb,iqcc,iqdd
integer, dimension(kl) :: iters
integer itrc,itr,ng,ng4,g,k,jj,i,j,iq
integer klimc,knew,klim,ir,ic
integer nc,ifc,n,iq_a,iq_b,iq_c,iq_d
real, dimension(ifull+iextra,kl), intent(inout) :: iv
real, dimension(ifull,kl), intent(in) :: ihelm,jrhs
real, dimension(ifull,kl) :: irhs
real, dimension(ifull), intent(in) :: izz,izzn,izze,izzw,izzs
real, dimension(ifull,maxcolour,kl) :: rhelmc,rhsc
real, dimension(ifull,maxcolour) :: zznc,zzec,zzwc,zzsc
real, dimension(ifull+iextra,kl) :: vdum
real, dimension(mg_minsize,3) :: zzncu,zzecu,zzwcu,zzscu
real, dimension(mg_minsize,3,kl) :: rhelmcu
real, dimension(mg_maxsize,kl,mg_maxlevel) :: v
real, dimension(mg_maxsize,kl,2:mg_maxlevel) :: rhs
real, dimension(mg_maxsize,kl,mg_maxlevel) :: helm
real, dimension(mg_maxsize,kl) :: w,dsol
real, dimension(kl) :: smax,smax_g,dsolmax,dsolmax_g
real, dimension(kl) :: smin,smin_g,savg

call start_log(helm_begin)

if (sorfirst.or.zzfirst) then
  write(6,*) "ERROR: mghelm requires mgsor_init and mgzz_init to be called first"
  call ccmpi_abort(-1)
end if

! zz*(DIV^2 v) - helm*v = rhs
! zz*(DIV^2 vd+v0) - helm*(vd+v0) = rhs
! zz*(DIV^2 vd) - helm*vd = rhs + helm*v0 - zz*(DIV^2 v0)

! Prepare input arrays
klim=kl
do k=1,kl
  smax(k)=maxval(iv(1:ifull,k))
  smin(k)=minval(iv(1:ifull,k))
end do
call ccmpi_allreduce(smax(1:kl),smax_g(1:kl),"max",comm_world)
call ccmpi_allreduce(smin(1:kl),smin_g(1:kl),"min",comm_world)

! JLM suggestion
do k=1,kl
  savg(k)=0.5*(smax_g(k)+smin_g(k))
  iv(1:ifull,k)=iv(1:ifull,k)-savg(k)
  irhs(:,k)=jrhs(:,k)+(ihelm(:,k)-izz-izzn-izzs-izze-izzw)*savg(k)
end do

iters=0

!itsave2=0
!itserr2=9.E9
!itstest=1
!itc=0

helm(1:ifull,:,1)=ihelm(1:ifull,:)
call mgcollect(1,helm(:,:,1))
do g=1,mg_maxlevel-1
  ng4=mg(g)%ifull_fine
  iqaa(1:ng4)=          mg(g)%fine    
  iqbb(1:ng4)= mg(g)%in(mg(g)%fine)
  iqcc(1:ng4)= mg(g)%ie(mg(g)%fine)
  iqdd(1:ng4)=mg(g)%ine(mg(g)%fine)
  helm(1:ng4,1:kl,g+1)=0.25*(helm(iqaa(1:ng4),1:kl,g)+helm(iqbb(1:ng4),1:kl,g) &
                            +helm(iqcc(1:ng4),1:kl,g)+helm(iqdd(1:ng4),1:kl,g))
  call mgcollect(g+1,helm(:,:,g+1))
end do

do nc=1,maxcolour
  ifc=ifullx(nc)
  zznc(1:ifc,nc) =izzn(iqx(1:ifc,nc))
  zzwc(1:ifc,nc) =izzw(iqx(1:ifc,nc))
  zzec(1:ifc,nc) =izze(iqx(1:ifc,nc))
  zzsc(1:ifc,nc) =izzs(iqx(1:ifc,nc))
  do k=1,kl
    rhelmc(1:ifc,nc,k)=1./(ihelm(iqx(1:ifc,nc),k)-izz(iqx(1:ifc,nc)))
    rhsc(1:ifc,nc,k) =irhs(iqx(1:ifc,nc),k)
  end do
end do

g=mg_maxlevel
do nc=1,3
  ifc=mg_ifullc(nc)
  zzncu(1:ifc,nc)=mg(g)%zzn(col_iq(1:ifc,nc))
  zzscu(1:ifc,nc)=mg(g)%zzs(col_iq(1:ifc,nc))
  zzecu(1:ifc,nc)=mg(g)%zze(col_iq(1:ifc,nc))
  zzwcu(1:ifc,nc)=mg(g)%zzw(col_iq(1:ifc,nc))
  do k=1,kl
    rhelmcu(1:ifc,nc,k)=1./(helm(col_iq(1:ifc,nc),k,g)-mg(g)%zz(col_iq(1:ifc,nc)))
  end do
end do

call bounds(iv,klim=klim)

! Main loop
do itr=1,itr_mg

  do nc=1,maxcolour
    ifc=ifullx(nc)
    do k=1,klim
      dsol(iqx(1:ifc,nc),k)=( zznc(1:ifc,nc)*iv(iqn(1:ifc,nc),k) + zzwc(1:ifc,nc)*iv(iqw(1:ifc,nc),k)    &
                            + zzec(1:ifc,nc)*iv(iqe(1:ifc,nc),k) + zzsc(1:ifc,nc)*iv(iqs(1:ifc,nc),k)    &
                            - rhsc(1:ifc,nc,k) )*rhelmc(1:ifc,nc,k) - iv(iqx(1:ifc,nc),k)
      iv(iqx(1:ifc,nc),k)=iv(iqx(1:ifc,nc),k)+dsol(iqx(1:ifc,nc),k)
    end do
    call bounds(iv,klim=klim,colour=nc)

  end do

  ! residual
  do k=1,klim
    w(1:ifull,k)=-izzn*iv(in,k)-izzw*iv(iw,k)-izze*iv(ie,k)-izzs*iv(is,k)+irhs(:,k)+iv(1:ifull,k)*(ihelm(:,k)-izz)
  end do

  ! fine grid
  g=1
  ng=mg(g)%ifull

  ! For when the inital grid cannot be upscaled
  call mgcollect(g,w,klim=klim)

  ! restriction
  ! (since this always operates within a panel, then ine = ien is always true)
  ng4=mg(g)%ifull_fine
  iqaa(1:ng4)=          mg(g)%fine    
  iqbb(1:ng4)= mg(g)%in(mg(g)%fine)
  iqcc(1:ng4)= mg(g)%ie(mg(g)%fine)
  iqdd(1:ng4)=mg(g)%ine(mg(g)%fine)
  rhs(1:ng4,1:klim,g+1)=0.25*(w(iqaa(1:ng4),1:klim)+w(iqbb(1:ng4),1:klim) &
                             +w(iqcc(1:ng4),1:klim)+w(iqdd(1:ng4),1:klim))

  ! upscale grid
  do g=2,mg_maxlevel-1
  
    ng=mg(g)%ifull

    ! merge grids if insufficent points on this processor
    call mgcollect(g,rhs(:,:,g),klim=klim)
    
    ! update
    ! possibly use colours here, although latency might be an issue as the grid size
    ! becomes more coarse
    ! assume zero for first guess of residual (also avoids additional bounds call)
    !v(:,1:klim,g)=0.
    do k=1,klim
      v(1:ng,k,g)=-rhs(1:ng,k,g)/(helm(1:ng,k,g)-mg(g)%zz)
    end do

    ! residual
    call mgbounds(g,v(:,:,g),klim=klim)
    do k=1,klim
      w(1:ng,k)=-mg(g)%zze*v(mg(g)%ie,k,g)-mg(g)%zzw*v(mg(g)%iw,k,g)     &
                -mg(g)%zzn*v(mg(g)%in,k,g)-mg(g)%zzs*v(mg(g)%is,k,g)
                !+rhs(1:ng,k,g)+(helm(1:ng,k,g)-mg(g)%zz)*v(1:ng,k,g)
    end do
   
    ! restriction
    ! (calculate finer grid before mgcollect as the messages sent/recv are shorter)
    ng4=mg(g)%ifull_fine
    iqaa(1:ng4)=          mg(g)%fine
    iqbb(1:ng4)= mg(g)%in(mg(g)%fine)
    iqcc(1:ng4)= mg(g)%ie(mg(g)%fine)
    iqdd(1:ng4)=mg(g)%ine(mg(g)%fine)
    rhs(1:ng4,1:klim,g+1)=0.25*(w(iqaa(1:ng4),1:klim)+w(iqbb(1:ng4),1:klim) &
                               +w(iqcc(1:ng4),1:klim)+w(iqdd(1:ng4),1:klim))
  end do

  ! solve coarse grid
  g=mg_maxlevel
  ng=mg(g)%ifull

  ! ensure all processors have a copy of the coarse grid
  call mgcollect(g,rhs(:,:,g),klim=klim)

  do k=1,klim
    v(1:ng,k,g)=-rhs(1:ng,k,g)/(helm(1:ng,k,g)-mg(g)%zz)
    smax(k)=maxval(v(1:ng,k,g))
    smin(k)=minval(v(1:ng,k,g))
  end do

  klimc=klim
  do itrc=1,itr_max
    do nc=1,3 ! always three colours for coarse multigrid
      ifc=mg_ifullc(nc)
      do k=1,klimc
        ! update
        dsol(col_iq(1:ifc,nc),k)=(zzecu(1:ifc,nc)*v(col_iqe(1:ifc,nc),k,g)+zzwcu(1:ifc,nc)*v(col_iqw(1:ifc,nc),k,g)  &
                                 +zzncu(1:ifc,nc)*v(col_iqn(1:ifc,nc),k,g)+zzscu(1:ifc,nc)*v(col_iqs(1:ifc,nc),k,g)  &
                                 -rhs(col_iq(1:ifc,nc),k,g))*rhelmcu(1:ifc,nc,k)-v(col_iq(1:ifc,nc),k,g)
        v(col_iq(1:ifc,nc),k,g)=v(col_iq(1:ifc,nc),k,g)+dsol(col_iq(1:ifc,nc),k)
      end do
      ! no call to bounds since all points are on this processor
    end do
    do k=1,klimc  
      dsolmax(k)=maxval(abs(dsol(1:ng,k)))
    end do    ! convergence test
    knew=klimc
    do k=klimc,1,-1
      if (dsolmax(k)>=restol*(smax(k)-smin(k))) exit
      knew=k-1      
    end do
    klimc=knew
    if (klimc<1) exit
  end do
 
  ! interpolation
  ! all points are present on this processor, so no need for MPI comms
  ! (there maybe issues interpolating near vertices, where we could
  !  possibly use a three point interpolation)
  ng4=mg(g)%ifull_coarse
  do k=1,klim
    w(1:ng4,k)= mg(g)%wgt_a*v(mg(g)%coarse_a,k,g) + mg(g)%wgt_bc*v(mg(g)%coarse_b,k,g) &
             + mg(g)%wgt_bc*v(mg(g)%coarse_c,k,g) +  mg(g)%wgt_d*v(mg(g)%coarse_d,k,g)
  end do

  ! downscale grid
  do g=mg_maxlevel-1,2,-1

    ng=mg(g)%ifull

    ! extension
    ! No mgbounds as the v halp has already been updated and
    ! the coarse interpolation also updates the w halo
    v(1:ng4,1:klim,g)=v(1:ng4,1:klim,g)+w(1:ng4,1:klim)

    ! post smoothing
    do k=1,klim
      v(1:ng,k,g)=(mg(g)%zze*v(mg(g)%ie,k,g)+mg(g)%zzw*v(mg(g)%iw,k,g) &
                  +mg(g)%zzn*v(mg(g)%in,k,g)+mg(g)%zzs*v(mg(g)%is,k,g) &
                  -rhs(1:ng,k,g))/(helm(1:ng,k,g)-mg(g)%zz)
    end do

    call mgbounds(g,v(:,:,g),klim=klim,corner=.true.)

    ! interpolation
    ! (there maybe issues interpolating near vertices, where we could
    !  possibly use a three point interpolation)
    ! JLM suggests using a three point (plane) fit instead of the 4
    ! point bi-linear fit.  This treats the vertices the same as the
    ! other points.
    ng4=mg(g)%ifull_coarse
    do k=1,klim
      w(1:ng4,k)= mg(g)%wgt_a*v(mg(g)%coarse_a,k,g) + mg(g)%wgt_bc*v(mg(g)%coarse_b,k,g) &
               + mg(g)%wgt_bc*v(mg(g)%coarse_c,k,g) +  mg(g)%wgt_d*v(mg(g)%coarse_d,k,g)
    end do

  end do


  ! fine grid
  g=1
  ng=mg(g)%ifull
  
  ! convert back to local grid if required and update halo
  if (mg(g)%merge_len>1) then
    vdum=0.
    do n=1,npan
      ir=mod(mg(g)%merge_pos(n)-1,mg(g)%merge_row)+1   ! index for proc row
      ic=(mg(g)%merge_pos(n)-1)/mg(g)%merge_row+1      ! index for proc col
      do jj=1,jpan
        iq_a=1+(jj-1)*ipan+(n-1)*ipan*jpan
        iq_b=jj*ipan+(n-1)*ipan*jpan
        iq_c=1+(ir-1)*ipan+(jj-1+(ic-1)*jpan)*jpan+(n-1)*ipan*jpan*mg(g)%merge_len
        iq_d=ir*ipan+(jj-1+(ic-1)*jpan)*jpan+(n-1)*ipan*jpan*mg(g)%merge_len
        vdum(iq_a:iq_b,1:klim)=w(iq_c:iq_d,1:klim)
      end do
      j=1
      do i=1,ipan
        iq_a=i+(j-1)*ipan+(n-1)*ipan*jpan
        iq_c=i+(ir-1)*ipan+(j-1+(ic-1)*jpan)*jpan+(n-1)*ipan*jpan*mg(g)%merge_len
        iq_b=is(iq_a)
        iq_d=mg(1)%is(iq_c)
        vdum(iq_b,1:klim)=w(iq_d,1:klim)
      end do
      j=jpan
      do i=1,ipan
        iq_a=i+(j-1)*ipan+(n-1)*ipan*jpan
        iq_c=i+(ir-1)*ipan+(j-1+(ic-1)*jpan)*jpan+(n-1)*ipan*jpan*mg(g)%merge_len
        iq_b=in(iq_a)
        iq_d=mg(1)%in(iq_c)
        vdum(iq_b,1:klim)=w(iq_d,1:klim)
      end do  
      i=1
      do j=1,jpan
        iq_a=i+(j-1)*ipan+(n-1)*ipan*jpan
        iq_c=i+(ir-1)*ipan+(j-1+(ic-1)*jpan)*jpan+(n-1)*ipan*jpan*mg(g)%merge_len
        iq_b=iw(iq_a)
        iq_d=mg(1)%iw(iq_c)
        vdum(iq_b,1:klim)=w(iq_d,1:klim)
      end do
      i=ipan
      do j=1,jpan
        iq_a=i+(j-1)*ipan+(n-1)*ipan*jpan
        iq_c=i+(ir-1)*ipan+(j-1+(ic-1)*jpan)*jpan+(n-1)*ipan*jpan*mg(g)%merge_len
        iq_b=ie(iq_a)
        iq_d=mg(1)%ie(iq_c)
        vdum(iq_b,1:klim)=w(iq_d,1:klim)
      end do
    end do
    w(1:ifull+iextra,1:klim)=vdum(1:ifull+iextra,1:klim)    
    ! extension
    iv(1:ifull+iextra,1:klim)=iv(1:ifull+iextra,1:klim)+w(1:ifull+iextra,1:klim)
  else
    vdum=0.
    do n=1,npan
      j=1
      do i=1,ipan
        iq=indx(i,j,n-1,ipan,jpan)
        iq_a=is(iq)
        iq_b=mg(1)%is(iq)
        vdum(iq_a,1:klim)=w(iq_b,1:klim)
      end do
      j=jpan
      do i=1,ipan
        iq=indx(i,j,n-1,ipan,jpan)
        iq_a=in(iq)
        iq_b=mg(1)%in(iq)
        vdum(iq_a,1:klim)=w(iq_b,1:klim)
      end do  
      i=1
      do j=1,jpan
        iq=indx(i,j,n-1,ipan,jpan)
        iq_a=iw(iq)
        iq_b=mg(1)%iw(iq)
        vdum(iq_a,1:klim)=w(iq_b,1:klim)
      end do
      i=ipan
      do j=1,jpan
        iq=indx(i,j,n-1,ipan,jpan)
        iq_a=ie(iq)
        iq_b=mg(1)%ie(iq)
        vdum(iq_a,1:klim)=w(iq_b,1:klim)
      end do
    end do
    w(ifull+1:ifull+iextra,1:klim)=vdum(ifull+1:ifull+iextra,1:klim)
    ! extension
    iv(1:ifull+iextra,1:klim)=iv(1:ifull+iextra,1:klim)+w(1:ifull+iextra,1:klim)
  end if

  ! post smoothing
  do nc=1,maxcolour
    ifc=ifullx(nc)
    do k=1,klim
      dsol(iqx(1:ifc,nc),k)=( zznc(1:ifc,nc)*iv(iqn(1:ifc,nc),k) + zzwc(1:ifc,nc)*iv(iqw(1:ifc,nc),k)    &
                            + zzec(1:ifc,nc)*iv(iqe(1:ifc,nc),k) + zzsc(1:ifc,nc)*iv(iqs(1:ifc,nc),k)    &
                            - rhsc(1:ifc,nc,k))*rhelmc(1:ifc,nc,k) - iv(iqx(1:ifc,nc),k)
      iv(iqx(1:ifc,nc),k) = iv(iqx(1:ifc,nc),k) + dsol(iqx(1:ifc,nc),k)
    end do
    call bounds(iv,klim=klim,colour=nc)
  end do

  ! test for convergence
  do k=1,klim
    dsolmax(k)=maxval(abs(dsol(1:ifull,k)))
  end do
  call ccmpi_allreduce(dsolmax(1:klim),dsolmax_g(1:klim),"max",comm_world)
  knew=klim
  do k=klim,1,-1
    iters(k)=itr
    if (dsolmax_g(k)>=restol*(smax_g(k)-smin_g(k))) exit
    knew=k-1
  end do
  klim=knew
  if (klim<1) exit
  
end do

! SOR in case MG fails to converge
if (klim>0) then
  do itr=itr_mg+1,itr_max

    do nc=1,maxcolour
      ifc=ifullx(nc)
      do k=1,klim
        dsol(iqx(1:ifc,nc),k)=( zznc(1:ifc,nc)*iv(iqn(1:ifc,nc),k) + zzwc(1:ifc,nc)*iv(iqw(1:ifc,nc),k)    &
                              + zzec(1:ifc,nc)*iv(iqe(1:ifc,nc),k) + zzsc(1:ifc,nc)*iv(iqs(1:ifc,nc),k)    &
                              - rhsc(1:ifc,nc,k) )*rhelmc(1:ifc,nc,k) - iv(iqx(1:ifc,nc),k)
        iv(iqx(1:ifc,nc),k) = iv(iqx(1:ifc,nc),k) + dsol(iqx(1:ifc,nc),k)
      end do
      call bounds(iv,klim=klim,colour=nc)
    end do

    ! test for convergence
    do k=1,klim
      dsolmax(k)=maxval(abs(dsol(1:ifull,k)))
    end do    ! convergence test
    call ccmpi_allreduce(dsolmax(1:klim),dsolmax_g(1:klim),"max",comm_world)
    knew=klim
    do k=klim,1,-1
      iters(k)=itr
      if (dsolmax_g(k)>=restol*(smax_g(k)-smin_g(k))) exit
      knew=k-1
    end do
    klim=knew
    if (klim<1) exit
  
  end do
end if

! JLM suggestion
do k=1,kl
  iv(1:ifull,k)=iv(1:ifull,k)+savg(k)
end do

if (myid==0) then
  if (ktau<6.or.iters(1)>itr_mg) then
    do k=1,kl
      write(6,*) "mg ktau,k,iter ",ktau,k,iters(k)
    end do
    !write(6,*) "itc ",itc
  end if
end if

call end_log(helm_end)

return
end subroutine mghelm

! This version is for the ocean and ice
subroutine mgmlo(neta,ipice,iyy,iyyn,iyys,iyye,iyyw,izz,izzn,izzs,izze,izzw,ihh,irhs,tol,itol,totits,maxglobseta,maxglobip, &
                 ipmax,ee,dd,comm_mlo)

use cc_mpi
use indices_m

implicit none

include 'newmpar.h'

integer, intent(in) :: comm_mlo
integer, intent(out) :: totits
integer itr,itrc,g,ng,ng4,n,i,j,ir,ic,jj,iq
integer iq_a,iq_b,iq_c,iq_d
integer ifc,nc
integer, dimension(mg_maxsize) :: iqaa,iqbb,iqcc,iqdd
real, intent(in) :: tol,itol
real, intent(out) :: maxglobseta,maxglobip
real, dimension(ifull+iextra), intent(inout) :: neta,ipice
real, dimension(ifull+iextra), intent(in) :: ee,dd
real, dimension(ifull), intent(in) :: ipmax
real, dimension(ifull), intent(in) :: iyy,iyyn,iyys,iyye,iyyw
real, dimension(ifull,2), intent(in) :: izz,izzn,izzs,izze,izzw
real, dimension(ifull), intent(in) :: ihh
real, dimension(ifull,2), intent(in) :: irhs
real, dimension(ifull+iextra) :: vdum
real, dimension(mg_maxsize) :: au,bu,cu
real, dimension(mg_maxsize,1,mg_maxlevel) :: v
real, dimension(mg_maxsize,5,mg_maxlevel) :: yy
real, dimension(mg_maxsize,7) :: w
real, dimension(mg_maxsize,mg_maxlevel) :: zz,zzn,zzs,zze,zzw
real, dimension(mg_maxsize,mg_maxlevel) :: hh
real, dimension(mg_maxsize,mg_maxlevel) :: rhs
real, dimension(mg_maxsize,2) :: dsol,new
real, dimension(ifull+iextra,2) :: dumc,dumd
real, dimension(2) :: alpha
real, dimension(2) :: dsolmax,dsolmax_g

real, parameter :: dfac=0.25      ! adjustment for grid spacing

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

! zz*(DIV^2 ipice) = rhs
! ipice = i0 + f

! zz*(DIV^2 i0) = rhs - F
! zz*(DIV^2 f) = F

! MJT notes - There seems to be no advantage using a multi-grid for seaice.  Hence
! we only use the multi-grid for the ocean free surface height

alpha(1)=1. ! for ocean
alpha(2)=1. ! for sea-ice

yy(1:ifull,1,1)=iyy(1:ifull)
yy(1:ifull,2,1)=iyyn(1:ifull)
yy(1:ifull,3,1)=iyys(1:ifull)
yy(1:ifull,4,1)=iyye(1:ifull)
yy(1:ifull,5,1)=iyyw(1:ifull)
call mgcollect(1,yy(:,1:5,1),gmode=1)
do g=1,mg_maxlevel-1
  ng4=mg(g)%ifull_fine
  iqaa(1:ng4)=          mg(g)%fine    
  iqbb(1:ng4)= mg(g)%in(mg(g)%fine)
  iqcc(1:ng4)= mg(g)%ie(mg(g)%fine)
  iqdd(1:ng4)=mg(g)%ine(mg(g)%fine)
  yy(1:ng4,1:5,g+1)=0.25*dfac*(yy(iqaa(1:ng4),1:5,g)+yy(iqbb(1:ng4),1:5,g) &
                              +yy(iqcc(1:ng4),1:5,g)+yy(iqdd(1:ng4),1:5,g))
  call mgcollect(g+1,yy(:,1:5,g+1),gmode=1)
end do

dumc(1:ifull,1)=neta(1:ifull)
dumc(1:ifull,2)=ipice(1:ifull)
call bounds(dumc(:,1:2),gmode=1)
neta(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,1)
ipice(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,2)

! Main loop
do itr=1,itr_mgice

  au(1:ifull)=iyy

  do nc=1,maxcolour
    ifc=ifullx(nc)
    bu(iqx(1:ifc,nc))=izz(iqx(1:ifc,nc),1)+ihh(iqx(1:ifc,nc))                                            &
                      +iyyn(iqx(1:ifc,nc))*neta(iqn(1:ifc,nc))+iyys(iqx(1:ifc,nc))*neta(iqs(1:ifc,nc)) &
                      +iyye(iqx(1:ifc,nc))*neta(iqe(1:ifc,nc))+iyyw(iqx(1:ifc,nc))*neta(iqw(1:ifc,nc))
    cu(iqx(1:ifc,nc))=izzn(iqx(1:ifc,nc),1)*neta(iqn(1:ifc,nc))+izzs(iqx(1:ifc,nc),1)*neta(iqs(1:ifc,nc)) &
                     +izze(iqx(1:ifc,nc),1)*neta(iqe(1:ifc,nc))+izzw(iqx(1:ifc,nc),1)*neta(iqw(1:ifc,nc)) &
                     -irhs(iqx(1:ifc,nc),1)        
    new(iqx(1:ifc,nc),1) = -2.*cu(iqx(1:ifc,nc))/(bu(iqx(1:ifc,nc))+sqrt(bu(iqx(1:ifc,nc))**2-4.*au(iqx(1:ifc,nc))*cu(iqx(1:ifc,nc))))
    
    new(iqx(1:ifc,nc),2) = 0.
    where (izz(iqx(1:ifc,nc),2)/=0.)
      new(iqx(1:ifc,nc),2) = ( -izzn(iqx(1:ifc,nc),2)*ipice(iqn(1:ifc,nc))-izzs(iqx(1:ifc,nc),2)*ipice(iqs(1:ifc,nc)) &
                               -izze(iqx(1:ifc,nc),2)*ipice(iqe(1:ifc,nc))-izzw(iqx(1:ifc,nc),2)*ipice(iqw(1:ifc,nc)) &
                               + irhs(iqx(1:ifc,nc),2) ) / izz(iqx(1:ifc,nc),2)
    end where

    ! cavitating fluid
    new(iqx(1:ifc,nc),2)=max(min(new(iqx(1:ifc,nc),2),ipmax(iqx(1:ifc,nc))),0.)
    
    new(iqx(1:ifc,nc),1)=max(new(iqx(1:ifc,nc),1),-dd(iqx(1:ifc,nc)))*ee(iqx(1:ifc,nc))
    
    neta(iqx(1:ifc,nc))=alpha(1)*new(iqx(1:ifc,nc),1)+(1.-alpha(1))*neta(iqx(1:ifc,nc))
    ipice(iqx(1:ifc,nc))=alpha(2)*new(iqx(1:ifc,nc),2)+(1.-alpha(2))*ipice(iqx(1:ifc,nc))

    dumc(1:ifull,1)=neta(1:ifull)
    dumc(1:ifull,2)=ipice(1:ifull)
    call bounds(dumc(:,1:2),gmode=1,colour=nc)
    neta(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,1)
    ipice(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,2)
  end do

  w(1:ifull,2)= izz(:,1)+ iyy*neta(1:ifull)
  w(1:ifull,3)=izzn(:,1)+iyyn*neta(1:ifull)
  w(1:ifull,4)=izzs(:,1)+iyys*neta(1:ifull)
  w(1:ifull,5)=izze(:,1)+iyye*neta(1:ifull)
  w(1:ifull,6)=izzw(:,1)+iyyw*neta(1:ifull)
  w(1:ifull,7)=ihh+iyy*neta(1:ifull)+iyyn*neta(in)+iyys*neta(is)+iyye*neta(ie)+iyyw*neta(iw)

  ! residual
  w(1:ifull,1)=-neta(1:ifull)*(     iyy*neta(1:ifull)+     iyyn*neta(in)+     iyys*neta(is)+     iyye*neta(ie)+     iyyw*neta(iw)) &
                             -(izz(:,1)*neta(1:ifull)+izzn(:,1)*neta(in)+izzs(:,1)*neta(is)+izze(:,1)*neta(ie)+izzw(:,1)*neta(iw)) &
                             -ihh*neta(1:ifull)+irhs(:,1)
                             
  w(1:ifull,1)=w(1:ifull,1)*ee(1:ifull)

  ! fine grid
  g=1
  ng=mg(g)%ifull

  ! For when the inital grid cannot be upscaled
  call mgcollect(g,w(:,1:7),gmode=1)
  
  ! restriction
  ! (since this always operates within a panel, then ine = ien is always true)
  ng4=mg(g)%ifull_fine
  iqaa(1:ng4)=          mg(g)%fine    
  iqbb(1:ng4)= mg(g)%in(mg(g)%fine)
  iqcc(1:ng4)= mg(g)%ie(mg(g)%fine)
  iqdd(1:ng4)=mg(g)%ine(mg(g)%fine)
  rhs(1:ng4,g+1)=0.25*(w(iqaa(1:ng4),1)+w(iqbb(1:ng4),1) &
                      +w(iqcc(1:ng4),1)+w(iqdd(1:ng4),1))
  zz(1:ng4,g+1) =0.25*dfac*(w(iqaa(1:ng4),2)+w(iqbb(1:ng4),2) &
                           +w(iqcc(1:ng4),2)+w(iqdd(1:ng4),2))
  zzn(1:ng4,g+1)=0.25*dfac*(w(iqaa(1:ng4),3)+w(iqbb(1:ng4),3) &
                           +w(iqcc(1:ng4),3)+w(iqdd(1:ng4),3))
  zzs(1:ng4,g+1)=0.25*dfac*(w(iqaa(1:ng4),4)+w(iqbb(1:ng4),4) &
                           +w(iqcc(1:ng4),4)+w(iqdd(1:ng4),4))
  zze(1:ng4,g+1)=0.25*dfac*(w(iqaa(1:ng4),5)+w(iqbb(1:ng4),5) &
                           +w(iqcc(1:ng4),5)+w(iqdd(1:ng4),5))
  zzw(1:ng4,g+1)=0.25*dfac*(w(iqaa(1:ng4),6)+w(iqbb(1:ng4),6) &
                           +w(iqcc(1:ng4),6)+w(iqdd(1:ng4),6))
  hh(1:ng4,g+1)=0.25*(w(iqaa(1:ng4),7)+w(iqbb(1:ng4),7) &
                     +w(iqcc(1:ng4),7)+w(iqdd(1:ng4),7))

  ! upscale grid
  do g=2,mg_maxlevel-1
  
    ng=mg(g)%ifull

    ! merge grids if insufficent points on this processor
    w(1:ng4,1)  =rhs(1:ng4,g)
    w(1:ng4,2)  =zz(1:ng4,g)
    w(1:ng4,3)  =zzn(1:ng4,g)
    w(1:ng4,4)  =zzs(1:ng4,g)
    w(1:ng4,5)  =zze(1:ng4,g)
    w(1:ng4,6)  =zzw(1:ng4,g)
    w(1:ng4,7)  =hh(1:ng4,g)
    call mgcollect(g,w(:,1:7),gmode=1)
    rhs(1:ng,g)    =w(1:ng,1)
    zz(1:ng,g)     =w(1:ng,2)
    zzn(1:ng,g)    =w(1:ng,3)
    zzs(1:ng,g)    =w(1:ng,4)
    zze(1:ng,g)    =w(1:ng,5)
    zzw(1:ng,g)    =w(1:ng,6)
    hh(1:ng,g)     =w(1:ng,7)

    !v(:,1,g)=0.

    ! update
    ! possibly use colours here, although v is reset to zero every iteration
    ! assume zero for first guess of residual (also avoids additional bounds call)
    au(1:ng)=yy(1:ng,1,g)
    bu(1:ng)=zz(1:ng,g)+hh(1:ng,g)
    cu(1:ng)=-rhs(1:ng,g)
    v(1:ng,1,g) = 0.
    where (bu(1:ng)/=0.)
      v(1:ng,1,g) = -2.*cu(1:ng)/(bu(1:ng)+sqrt(bu(1:ng)*bu(1:ng)-4.*au(1:ng)*cu(1:ng)))
    end where
   
    v(1:ng,1,g) = alpha(1)*v(1:ng,1,g)
    
    ! residual
    call mgbounds(g,v(:,1:1,g),gmode=1)

    w(1:ng,2)= zz(1:ng,g)+yy(1:ng,1,g)*v(1:ng,1,g)
    w(1:ng,3)=zzn(1:ng,g)+yy(1:ng,2,g)*v(1:ng,1,g)
    w(1:ng,4)=zzs(1:ng,g)+yy(1:ng,3,g)*v(1:ng,1,g)
    w(1:ng,5)=zze(1:ng,g)+yy(1:ng,4,g)*v(1:ng,1,g)
    w(1:ng,6)=zzw(1:ng,g)+yy(1:ng,5,g)*v(1:ng,1,g)
    w(1:ng,7)=hh(1:ng,g)+yy(1:ng,1,g)*v(1:ng,1,g)+yy(1:ng,2,g)*v(mg(g)%in,1,g)+yy(1:ng,3,g)*v(mg(g)%is,1,g) &
                                                 +yy(1:ng,4,g)*v(mg(g)%ie,1,g)+yy(1:ng,5,g)*v(mg(g)%iw,1,g)

    w(1:ng,1)=-v(1:ng,1,g)*(yy(1:ng,1,g)*v(1:ng,1,g)+yy(1:ng,2,g)*v(mg(g)%in,1,g)+yy(1:ng,3,g)*v(mg(g)%is,1,g)+yy(1:ng,4,g)*v(mg(g)%ie,1,g)+yy(1:ng,5,g)*v(mg(g)%iw,1,g)) &
                            -(zz(1:ng,g)*v(1:ng,1,g)+ zzn(1:ng,g)*v(mg(g)%in,1,g)+ zzs(1:ng,g)*v(mg(g)%is,1,g)+ zze(1:ng,g)*v(mg(g)%ie,1,g)+ zzw(1:ng,g)*v(mg(g)%iw,1,g)) &
                             -hh(1:ng,g)*v(1:ng,1,g)+rhs(1:ng,g)

    ! restriction
    ! (calculate finer grid before mgcollect as the messages sent/recv are shorter)
    ng4=mg(g)%ifull_fine
    iqaa(1:ng4)=          mg(g)%fine    
    iqbb(1:ng4)= mg(g)%in(mg(g)%fine)
    iqcc(1:ng4)= mg(g)%ie(mg(g)%fine)
    iqdd(1:ng4)=mg(g)%ine(mg(g)%fine)
    rhs(1:ng4,g+1)=0.25*(w(iqaa(1:ng4),1)+w(iqbb(1:ng4),1) &
                        +w(iqcc(1:ng4),1)+w(iqdd(1:ng4),1))
    zz(1:ng4,g+1)=0.25*dfac*(w(iqaa(1:ng4),2)+w(iqbb(1:ng4),2) &
                            +w(iqcc(1:ng4),2)+w(iqdd(1:ng4),2))
    zzn(1:ng4,g+1)=0.25*dfac*(w(iqaa(1:ng4),3)+w(iqbb(1:ng4),3) &
                             +w(iqcc(1:ng4),3)+w(iqdd(1:ng4),3))
    zzs(1:ng4,g+1)=0.25*dfac*(w(iqaa(1:ng4),4)+w(iqbb(1:ng4),4) &
                             +w(iqcc(1:ng4),4)+w(iqdd(1:ng4),4))
    zze(1:ng4,g+1)=0.25*dfac*(w(iqaa(1:ng4),5)+w(iqbb(1:ng4),5) &
                             +w(iqcc(1:ng4),5)+w(iqdd(1:ng4),5))
    zzw(1:ng4,g+1)=0.25*dfac*(w(iqaa(1:ng4),6)+w(iqbb(1:ng4),6) &
                             +w(iqcc(1:ng4),6)+w(iqdd(1:ng4),6))
    hh(1:ng4,g+1)=0.25*(w(iqaa(1:ng4),7)+w(iqbb(1:ng4),7) &
                       +w(iqcc(1:ng4),7)+w(iqdd(1:ng4),7))

  end do

  ! solve coarse grid
  g=mg_maxlevel
  ng=mg(g)%ifull

  ! ensure all processors have a copy of the coarse grid
  v(1:ng,1,g)=0.
  
  w(1:ng4,1)  =rhs(1:ng4,g)
  w(1:ng4,2)  =zz(1:ng4,g)
  w(1:ng4,3)  =zzn(1:ng4,g)
  w(1:ng4,4)  =zzs(1:ng4,g)
  w(1:ng4,5)  =zze(1:ng4,g)
  w(1:ng4,6)  =zzw(1:ng4,g)
  w(1:ng4,7)  =hh(1:ng4,g)
  call mgcollect(g,w(:,1:7),gmode=1)
  rhs(1:ng,g)    =w(1:ng,1)
  zz(1:ng,g)     =w(1:ng,2)
  zzn(1:ng,g)    =w(1:ng,3)
  zzs(1:ng,g)    =w(1:ng,4)
  zze(1:ng,g)    =w(1:ng,5)
  zzw(1:ng,g)    =w(1:ng,6)
  hh(1:ng,g)     =w(1:ng,7)


  do itrc=1,itr_max

    au(1:ng)=yy(1:ng,1,g)
    new(1:ng,1)=0.
    
    do nc=1,3
      ifc=mg_ifullc(nc)
      bu(col_iq(1:ifc,nc))=zz(col_iq(1:ifc,nc),g)+hh(col_iq(1:ifc,nc),g)                                                       &
                          +yy(col_iq(1:ifc,nc),2,g)*v(col_iqn(1:ifc,nc),1,g)+yy(col_iq(1:ifc,nc),3,g)*v(col_iqs(1:ifc,nc),1,g) &
                          +yy(col_iq(1:ifc,nc),4,g)*v(col_iqe(1:ifc,nc),1,g)+yy(col_iq(1:ifc,nc),5,g)*v(col_iqw(1:ifc,nc),1,g)
      cu(col_iq(1:ifc,nc))=zzn(col_iq(1:ifc,nc),g)*v(col_iqn(1:ifc,nc),1,g)+zzs(col_iq(1:ifc,nc),g)*v(col_iqs(1:ifc,nc),1,g) &
                          +zze(col_iq(1:ifc,nc),g)*v(col_iqe(1:ifc,nc),1,g)+zzw(col_iq(1:ifc,nc),g)*v(col_iqw(1:ifc,nc),1,g) &
                          -rhs(col_iq(1:ifc,nc),g)
      where (bu(col_iq(1:ifc,nc))/=0.)
        new(col_iq(1:ifc,nc),1) = -2.*cu(col_iq(1:ifc,nc))/(bu(col_iq(1:ifc,nc))+sqrt(bu(col_iq(1:ifc,nc))**2-4.*au(col_iq(1:ifc,nc))*cu(col_iq(1:ifc,nc))))
      end where
        
      dsol(col_iq(1:ifc,nc),1)=new(col_iq(1:ifc,nc),1)-v(col_iq(1:ifc,nc),1,g)
      v(col_iq(1:ifc,nc),1,g)=alpha(1)*new(col_iq(1:ifc,nc),1)+(1.-alpha(1))*v(col_iq(1:ifc,nc),1,g)
    end do
    
    dsolmax(1)=maxval(abs(dsol(1:ng,1)))
    if (dsolmax(1)<tol) exit

  end do
  
  ! interpolation
  ng4=mg(g)%ifull_coarse
  w(1:ng4,1)= mg(g)%wgt_a*v(mg(g)%coarse_a,1,g) + mg(g)%wgt_bc*v(mg(g)%coarse_b,1,g) &
           + mg(g)%wgt_bc*v(mg(g)%coarse_c,1,g) +  mg(g)%wgt_d*v(mg(g)%coarse_d,1,g)


  ! downscale grid
  do g=mg_maxlevel-1,2,-1

    ng=mg(g)%ifull

    ! extension
    ! No mgbounds as the v halp has already been updated and
    ! the coarse interpolation also updates the w halo
    v(1:ng4,1,g)=v(1:ng4,1,g)+w(1:ng4,1)

    ! post smoothing
    au(1:ng)=yy(1:ng,1,g)
    bu(1:ng)=zz(1:ng,g)+hh(1:ng,g)+yy(1:ng,2,g)*v(mg(g)%in,1,g)+yy(1:ng,3,g)*v(mg(g)%is,1,g)+yy(1:ng,4,g)*v(mg(g)%ie,1,g)+yy(1:ng,5,g)*v(mg(g)%iw,1,g)
    cu(1:ng)=zzn(1:ng,g)*v(mg(g)%in,1,g)+zzs(1:ng,g)*v(mg(g)%is,1,g)+zze(1:ng,g)*v(mg(g)%ie,1,g)+zzw(1:ng,g)*v(mg(g)%iw,1,g)-rhs(1:ng,g)
    
    new(1:ng,1) = 0.
    where (bu(1:ng)/=0.)
      new(1:ng,1) = -2.*cu(1:ng)/(bu(1:ng)+sqrt(bu(1:ng)*bu(1:ng)-4.*au(1:ng)*cu(1:ng)))
    end where
    
    v(1:ng,1,g)=alpha(1)*new(1:ng,1)+(1.-alpha(1))*v(1:ng,1,g)

    call mgbounds(g,v(:,1:1,g),corner=.true.,gmode=1)

    ! interpolation
    ng4=mg(g)%ifull_coarse
    w(1:ng4,1)= mg(g)%wgt_a*v(mg(g)%coarse_a,1,g) + mg(g)%wgt_bc*v(mg(g)%coarse_b,1,g) &
             + mg(g)%wgt_bc*v(mg(g)%coarse_c,1,g) +  mg(g)%wgt_d*v(mg(g)%coarse_d,1,g)
  end do


  ! fine grid
  g=1
  ng=mg(g)%ifull
  
  ! convert back to local grid if required and update halo
  if (mg(g)%merge_len>1) then
    vdum=0.
    do n=1,npan
      ir=mod(mg(g)%merge_pos(n)-1,mg(g)%merge_row)+1   ! index for proc row
      ic=(mg(g)%merge_pos(n)-1)/mg(g)%merge_row+1      ! index for proc col
      do jj=1,jpan
        iq_a=1+(jj-1)*ipan+(n-1)*ipan*jpan
        iq_b=jj*ipan+(n-1)*ipan*jpan
        iq_c=1+(ir-1)*ipan+(jj-1+(ic-1)*jpan)*jpan+(n-1)*ipan*jpan*mg(g)%merge_len
        iq_d=ir*ipan+(jj-1+(ic-1)*jpan)*jpan+(n-1)*ipan*jpan*mg(g)%merge_len
        vdum(iq_a:iq_b)=w(iq_c:iq_d,1)
      end do
      j=1
      do i=1,ipan
        iq_a=i+(j-1)*ipan+(n-1)*ipan*jpan
        iq_c=i+(ir-1)*ipan+(j-1+(ic-1)*jpan)*jpan+(n-1)*ipan*jpan*mg(g)%merge_len
        iq_b=is(iq_a)
        iq_d=mg(1)%is(iq_c)
        vdum(iq_b)=w(iq_d,1)
      end do
      j=jpan
      do i=1,ipan
        iq_a=i+(j-1)*ipan+(n-1)*ipan*jpan
        iq_c=i+(ir-1)*ipan+(j-1+(ic-1)*jpan)*jpan+(n-1)*ipan*jpan*mg(g)%merge_len
        iq_b=in(iq_a)
        iq_d=mg(1)%in(iq_c)
        vdum(iq_b)=w(iq_d,1)
      end do  
      i=1
      do j=1,jpan
        iq_a=i+(j-1)*ipan+(n-1)*ipan*jpan
        iq_c=i+(ir-1)*ipan+(j-1+(ic-1)*jpan)*jpan+(n-1)*ipan*jpan*mg(g)%merge_len
        iq_a=iw(iq_a)
        iq_b=mg(1)%iw(iq_c)
        vdum(iq_b)=w(iq_d,1)
      end do
      i=ipan
      do j=1,jpan
        iq_a=i+(j-1)*ipan+(n-1)*ipan*jpan
        iq_c=i+(ir-1)*ipan+(j-1+(ic-1)*jpan)*jpan+(n-1)*ipan*jpan*mg(g)%merge_len
        iq_a=ie(iq_a)
        iq_b=mg(1)%ie(iq_c)
        vdum(iq_b)=w(iq_d,1)
      end do
    end do
    w(1:ifull+iextra,1)=vdum(1:ifull+iextra)
    ! extension
    neta(1:ifull+iextra)=neta(1:ifull+iextra)+w(1:ifull+iextra,1)
  else
    vdum=0.
    do n=1,npan
      j=1
      do i=1,ipan
        iq=indx(i,j,n-1,ipan,jpan)
        iq_a=is(iq)
        iq_b=mg(1)%is(iq)
        vdum(iq_a)=w(iq_b,1)
      end do
      j=jpan
      do i=1,ipan
        iq=indx(i,j,n-1,ipan,jpan)
        iq_a=in(iq)
        iq_b=mg(1)%in(iq)
        vdum(iq_a)=w(iq_b,1)
      end do  
      i=1
      do j=1,jpan
        iq=indx(i,j,n-1,ipan,jpan)
        iq_a=iw(iq)
        iq_b=mg(1)%iw(iq)
        vdum(iq_a)=w(iq_b,1)
      end do
      i=ipan
      do j=1,jpan
        iq=indx(i,j,n-1,ipan,jpan)
        iq_a=ie(iq)
        iq_b=mg(1)%ie(iq)
        vdum(iq_a)=w(iq_b,1)
      end do
    end do
    w(ifull+1:ifull+iextra,1)=vdum(ifull+1:ifull+iextra)
    ! extension
    neta(1:ifull+iextra)=neta(1:ifull+iextra)+w(1:ifull+iextra,1)
  end if

  neta=max(neta,-dd)*ee
   
  
  ! post smoothing
  au(1:ifull)=iyy
  
  do nc=1,maxcolour
    ifc=ifullx(nc)
    bu(iqx(1:ifc,nc))=izz(iqx(1:ifc,nc),1)+ihh(iqx(1:ifc,nc))                                            &
                      +iyyn(iqx(1:ifc,nc))*neta(iqn(1:ifc,nc))+iyys(iqx(1:ifc,nc))*neta(iqs(1:ifc,nc)) &
                      +iyye(iqx(1:ifc,nc))*neta(iqe(1:ifc,nc))+iyyw(iqx(1:ifc,nc))*neta(iqw(1:ifc,nc))
    cu(iqx(1:ifc,nc))=izzn(iqx(1:ifc,nc),1)*neta(iqn(1:ifc,nc))+izzs(iqx(1:ifc,nc),1)*neta(iqs(1:ifc,nc)) &
                     +izze(iqx(1:ifc,nc),1)*neta(iqe(1:ifc,nc))+izzw(iqx(1:ifc,nc),1)*neta(iqw(1:ifc,nc)) &
                     -irhs(iqx(1:ifc,nc),1)        
    new(iqx(1:ifc,nc),1) = -2.*cu(iqx(1:ifc,nc))/(bu(iqx(1:ifc,nc))+sqrt(bu(iqx(1:ifc,nc))**2-4.*au(iqx(1:ifc,nc))*cu(iqx(1:ifc,nc))))
    new(iqx(1:ifc,nc),2) = 0.
    where (izz(iqx(1:ifc,nc),2)/=0.)
      new(iqx(1:ifc,nc),2) = ( -izzn(iqx(1:ifc,nc),2)*ipice(iqn(1:ifc,nc))-izzs(iqx(1:ifc,nc),2)*ipice(iqs(1:ifc,nc)) &
                               -izze(iqx(1:ifc,nc),2)*ipice(iqe(1:ifc,nc))-izzw(iqx(1:ifc,nc),2)*ipice(iqw(1:ifc,nc)) &
                               + irhs(iqx(1:ifc,nc),2) ) / izz(iqx(1:ifc,nc),2)
    end where

    ! cavitating fluid
    new(iqx(1:ifc,nc),2)=max(min(new(iqx(1:ifc,nc),2),ipmax(iqx(1:ifc,nc))),0.)

    new(iqx(1:ifc,nc),1)=max(new(iqx(1:ifc,nc),1),-dd(iqx(1:ifc,nc)))*ee(iqx(1:ifc,nc))
    new(iqx(1:ifc,nc),2)=new(iqx(1:ifc,nc),2)*ee(iqx(1:ifc,nc))
                             
    dsol(iqx(1:ifc,nc),1)=new(iqx(1:ifc,nc),1)-neta(iqx(1:ifc,nc))
    dsol(iqx(1:ifc,nc),2)=new(iqx(1:ifc,nc),2)-ipice(iqx(1:ifc,nc))
    neta(iqx(1:ifc,nc))=alpha(1)*new(iqx(1:ifc,nc),1)+(1.-alpha(1))*neta(iqx(1:ifc,nc))
    ipice(iqx(1:ifc,nc))=alpha(2)*new(iqx(1:ifc,nc),2)+(1.-alpha(2))*ipice(iqx(1:ifc,nc))

    dumc(1:ifull,1)=neta(1:ifull)
    dumc(1:ifull,2)=ipice(1:ifull)
    call bounds(dumc(:,1:2),gmode=1,colour=nc)
    neta(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,1)
    ipice(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,2)
  end do
  
  ! test for convergence
 
  dsolmax(1)=maxval(abs(dsol(1:ifull,1)))
  dsolmax(2)=maxval(abs(dsol(1:ifull,2)))
  call ccmpi_allreduce(dsolmax(1:2),dsolmax_g(1:2),"max",comm_mlo)
  if (dsolmax_g(1)<tol.and.dsolmax_g(2)<itol) exit
  
end do

! SOR if MG fails to converge
if (dsolmax_g(1)>=tol.or.dsolmax_g(2)>=itol) then
  do itr=itr_mgice+1,itr_max
    au(1:ifull)=iyy
  
    do nc=1,maxcolour
      ifc=ifullx(nc)
      bu(iqx(1:ifc,nc))=izz(iqx(1:ifc,nc),1)+ihh(iqx(1:ifc,nc))                                            &
                        +iyyn(iqx(1:ifc,nc))*neta(iqn(1:ifc,nc))+iyys(iqx(1:ifc,nc))*neta(iqs(1:ifc,nc)) &
                        +iyye(iqx(1:ifc,nc))*neta(iqe(1:ifc,nc))+iyyw(iqx(1:ifc,nc))*neta(iqw(1:ifc,nc))
      cu(iqx(1:ifc,nc))=izzn(iqx(1:ifc,nc),1)*neta(iqn(1:ifc,nc))+izzs(iqx(1:ifc,nc),1)*neta(iqs(1:ifc,nc)) &
                       +izze(iqx(1:ifc,nc),1)*neta(iqe(1:ifc,nc))+izzw(iqx(1:ifc,nc),1)*neta(iqw(1:ifc,nc)) &
                       -irhs(iqx(1:ifc,nc),1)        
      new(iqx(1:ifc,nc),1) = -2.*cu(iqx(1:ifc,nc))/(bu(iqx(1:ifc,nc))+sqrt(bu(iqx(1:ifc,nc))**2-4.*au(iqx(1:ifc,nc))*cu(iqx(1:ifc,nc))))
      new(iqx(1:ifc,nc),2) = 0.
      where (izz(iqx(1:ifc,nc),2)/=0.)
        new(iqx(1:ifc,nc),2) = ( -izzn(iqx(1:ifc,nc),2)*ipice(iqn(1:ifc,nc))-izzs(iqx(1:ifc,nc),2)*ipice(iqs(1:ifc,nc)) &
                                 -izze(iqx(1:ifc,nc),2)*ipice(iqe(1:ifc,nc))-izzw(iqx(1:ifc,nc),2)*ipice(iqw(1:ifc,nc)) &
                                 + irhs(iqx(1:ifc,nc),2) ) / izz(iqx(1:ifc,nc),2)
      end where

      ! cavitating fluid
      new(iqx(1:ifc,nc),2)=max(min(new(iqx(1:ifc,nc),2),ipmax(iqx(1:ifc,nc))),0.)

      new(iqx(1:ifc,nc),1)=max(new(iqx(1:ifc,nc),1),-dd(iqx(1:ifc,nc)))*ee(iqx(1:ifc,nc))
      new(iqx(1:ifc,nc),2)=new(iqx(1:ifc,nc),2)*ee(iqx(1:ifc,nc))
                             
      dsol(iqx(1:ifc,nc),1)=new(iqx(1:ifc,nc),1)-neta(iqx(1:ifc,nc))
      dsol(iqx(1:ifc,nc),2)=new(iqx(1:ifc,nc),2)-ipice(iqx(1:ifc,nc))
      neta(iqx(1:ifc,nc))=alpha(1)*new(iqx(1:ifc,nc),1)+(1.-alpha(1))*neta(iqx(1:ifc,nc))
      ipice(iqx(1:ifc,nc))=alpha(2)*new(iqx(1:ifc,nc),2)+(1.-alpha(2))*ipice(iqx(1:ifc,nc))

      dumc(1:ifull,1)=neta(1:ifull)
      dumc(1:ifull,2)=ipice(1:ifull)
      call bounds(dumc(:,1:2),gmode=1,colour=nc)
      neta(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,1)
      ipice(ifull+1:ifull+iextra)=dumc(ifull+1:ifull+iextra,2)
    end do
  
    ! test for convergence
 
    dsolmax(1)=maxval(abs(dsol(1:ifull,1)))
    dsolmax(2)=maxval(abs(dsol(1:ifull,2)))
    call ccmpi_allreduce(dsolmax(1:2),dsolmax_g(1:2),"max",comm_mlo)
    if (dsolmax_g(1)<tol.and.dsolmax_g(2)<itol) exit

  end do
end if

totits=itr
maxglobseta=dsolmax_g(1)
maxglobip=dsolmax_g(2)

return
end subroutine mgmlo

! Initialise multi-grid arrays
subroutine mgsor_init

use cc_mpi

implicit none

include 'newmpar.h'

integer g,np,iq,iqq,iql,iproc,xlen,nn,ii,jj
integer mipan,mjpan,hipan,hjpan,mil_g,ia,ja
integer i,j,n,mg_npan,mxpr,mypr
integer cid,ix,jx,colour,rank,ncol,nrow
integer npanx,na,nx,ny,sii,eii,sjj,ejj
integer drow,dcol,maxmergelen
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
! Let coarsest grid have 4 (2x2) points on each face. If grid size
! is not a power of 2, just make it as small as possible.
mg_maxlevel=1
g=1
do while (mod(il_g,2**g)==0)
  g=g+1
end do
mg_maxlevel=g

if (myid==0) then
  write(6,*) "Initialising multi-grid arrays"
  write(6,*) "il_g,mg_maxlevel ",il_g,mg_maxlevel
end if

allocate(mg(mg_maxlevel))
allocate(mg_bnds(0:nproc-1,mg_maxlevel))
do iproc=0,nproc-1
  do g=1,mg_maxlevel
    mg_bnds(iproc,g)%len=0
  end do
end do

hipan=mipan
hjpan=mjpan

! calculate fine grid for finest grid level
g=1
mg(1)%globgath=.false.
mg_npan=npan
mg(1)%merge_len=1
mg(1)%merge_row=1

allocate(mg(1)%fproc(mil_g,mil_g,0:npanels))
mg(1)%fproc(:,:,:)=fproc(:,:,:)

! check if coarse grid needs mgcollect
if (mod(mipan,2)/=0.or.mod(mjpan,2)/=0) then
  if (mod(mxpr,2)==0.and.mod(mypr,2)==0) then
   
   ! This case occurs when there are multiple processors on a panel.
   ! Consequently, npan should be 1 or 6.
   if (npan>1.and.npan<6) then
     write(6,*) "ERROR: Invalid gather4"
     call ccmpi_abort(-1)
   end if
    
    mg(1)%merge_len=4
    mg(1)%merge_row=2
    mxpr=mxpr/2
    mypr=mypr/2
    mipan=2*mipan
    mjpan=2*mjpan

    if (myid==0) then
      write(6,*) "Multi-grid gather4 at level           ",g,mipan,mjpan
    end if

    maxmergelen=npan
    allocate(mg(1)%merge_list(4,npan),mg(1)%merge_pos(npan))

    ! find my processor and surrounding members of the gather
    do n=1,npan
      nn=n-noff

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

       
      mg(1)%merge_list(1,n)=mg(1)%fproc(ii,      jj      ,nn)
      mg(1)%merge_list(2,n)=mg(1)%fproc(ii+hipan,jj      ,nn)
      mg(1)%merge_list(3,n)=mg(1)%fproc(ii,      jj+hjpan,nn)
      mg(1)%merge_list(4,n)=mg(1)%fproc(ii+hipan,jj+hjpan,nn)
       
      do j=1,mil_g,mjpan
        do i=1,mil_g,mipan
          cid=mg(1)%fproc(i+ix,j+jx,nn) ! processor in same merge position as myid
                                        ! we will maintain communications with this processor
          do ja=1,mjpan
            do ia=1,mipan
              ! update fproc map with processor that owns this data
              mg(1)%fproc(i+ia-1,j+ja-1,nn)=cid
            end do
          end do
        end do
      end do
        
      mg(1)%merge_pos(n)=-1
      do j=1,4
        if (mg(1)%merge_list(j,n)==myid) then
          mg(1)%merge_pos(n)=j
          exit
        end if
      end do
      if (mg(g)%merge_pos(n)<1) then
        write(6,*) "ERROR: Invalid merge_pos g,n,pos ",g,n,mg(g)%merge_pos(n)
        call ccmpi_abort(-1)
      end if
      
    end do

    ! fix any remaining panels
    if (npan==1) then
      do nn=0,npanels
        do j=1,mil_g,mjpan
          do i=1,mil_g,mipan
            cid=mg(g)%fproc(i+ix,j+jx,nn) ! processor in same merge position as myid
                                          ! we will maintain communications with this processor
            do ja=1,mjpan
              do ia=1,mipan
                ! update fproc map with processor that owns this data
                mg(g)%fproc(i+ia-1,j+ja-1,nn)=cid
              end do
            end do
          end do
        end do
      end do

      ! define local comm
      colour=mg(1)%merge_list(1,1)
      rank=mg(1)%merge_pos(1)-1
      call ccmpi_commsplit(mg(1)%comm,comm_world,colour,rank)
      
      ! MLO comm
      rank=-1
      jj=-1
      do ii=1,mg(1)%merge_len
        iproc=mg(1)%merge_list(ii,1)
        if (bnds(iproc)%mlomsk==1) then
          jj=jj+1
          if (iproc==myid) then
            rank=jj
            exit
          end if
        end if
      end do
      if (rank<0) then
        colour=nproc
        rank=0
      end if
      call ccmpi_commsplit(mg(1)%comm_mlo,comm_world,colour,rank)
      
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

allocate(mg(1)%fine(mg(1)%ifull_fine))
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

np=mg(1)%ifull
allocate(mg(1)%in(np),mg(1)%ie(np),mg(1)%is(np),mg(1)%iw(np))
allocate(mg(1)%ine(np),mg(1)%inw(np),mg(1)%ise(np),mg(1)%isw(np))

call mg_index(1,mil_g,mipan,mjpan)

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
      ja=(jj-1)/2+1
      do ii=1,2*mil_g,2
        ia=(ii-1)/2+1
        mg(g)%fproc(ia,ja,nn)=mg(g-1)%fproc(ii,jj,nn)
      end do
    end do
  end do


  ! first assume no gather for upscaled grid
  mg(g)%merge_len=1
  mg(g)%merge_row=1
  mg(g)%globgath=.false.
  
  ! check for multi-grid gather
  if (mod(mipan,2)/=0.or.mod(mjpan,2)/=0) then ! grid cannot be subdivided on current processor
  
    if (mod(mxpr,2)==0.and.mod(mypr,2)==0) then ! collect data over adjacent processors

      ! This case occurs when there are multiple processors on a panel.
      ! Consequently, npan should be 1 or 6.
      if (npan>1.and.npan<6) then
        write(6,*) "ERROR: Invalid gather4"
        call ccmpi_abort(-1)
      end if

      mg(g)%merge_len=4
      mg(g)%merge_row=2
      mxpr=mxpr/2
      mypr=mypr/2
      mipan=2*mipan
      mjpan=2*mjpan

      if (myid==0) then
        write(6,*) "Multi-grid gather4 at level           ",g,mipan,mjpan
      end if

      maxmergelen=npan
      allocate(mg(g)%merge_list(4,npan),mg(g)%merge_pos(npan))
      
      do n=1,npan
        nn=n-noff

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
        if (mod((ii-1)/hipan,2)/=0) ii=ii-hipan ! left corner of merge
        if (mod((jj-1)/hjpan,2)/=0) jj=jj-hjpan ! bottom corner of merge
        ix=ix-ii ! offset for myid from ii
        jx=jx-jj ! offset for myid from jj
       
        mg(g)%merge_list(1,n)=mg(g)%fproc(ii,      jj      ,nn)
        mg(g)%merge_list(2,n)=mg(g)%fproc(ii+hipan,jj      ,nn)
        mg(g)%merge_list(3,n)=mg(g)%fproc(ii,      jj+hjpan,nn)
        mg(g)%merge_list(4,n)=mg(g)%fproc(ii+hipan,jj+hjpan,nn)
 
        do j=1,mil_g,mjpan
          do i=1,mil_g,mipan
            cid=mg(g)%fproc(i+ix,j+jx,nn) ! processor in same merge position as myid
                                          ! we will maintain communications with this processor
            do ja=1,mjpan
              do ia=1,mipan
                ! update fproc map with processor that owns this data
                mg(g)%fproc(i+ia-1,j+ja-1,nn)=cid
              end do
            end do
          end do
        end do

        mg(g)%merge_pos(n)=-1
        do j=1,4
          if (mg(g)%merge_list(j,n)==myid) then
            mg(g)%merge_pos(n)=j
            exit
          end if
        end do
       
        if (mg(g)%merge_pos(n)<1) then
          write(6,*) "ERROR: Invalid merge_pos g,n,pos ",g,n,mg(g)%merge_pos(n)
          call ccmpi_abort(-1)
        end if

      end do
      
      ! fix any remaining panels
      if (npan==1) then
        do nn=0,npanels
          do j=1,mil_g,mjpan
            do i=1,mil_g,mipan
              cid=mg(g)%fproc(i+ix,j+jx,nn) ! processor in same merge position as myid
                                            ! we will maintain communications with this processor
              do ja=1,mjpan
                do ia=1,mipan
                  ! update fproc map with processor that owns this data
                  mg(g)%fproc(i+ia-1,j+ja-1,nn)=cid
                end do
              end do
            end do
          end do
        end do
      end if
    
    else if (.not.lglob) then ! collect all data to one processor
      lglob=.true.
#ifdef uniform_decomp
      mg(g)%merge_len=mxpr*mypr
#else
      mg(g)%merge_len=min(6*mxpr*mypr,nproc)
      mg(g)%globgath=.true.
      mg_npan=6
#endif

      mg(g)%merge_row=mxpr
      mipan=mipan*mxpr
      mjpan=mjpan*mypr
      mxpr=1
      mypr=1

      if (myid==0) then
        write(6,*) "Multi-grid gatherall at level         ",g,mipan,mjpan
      end if
      
      ! find gather members
#ifdef uniform_decomp
      maxmergelen=npan
      allocate(mg(g)%merge_list(mg(g)%merge_len,npan),mg(g)%merge_pos(npan))
      do n=1,npan
        nn=n-noff
        iqq=0
        do jj=1,mil_g,hjpan
          do ii=1,mil_g,hipan
            iqq=iqq+1
            mg(g)%merge_list(iqq,n)=mg(g)%fproc(ii,jj,nn)
          end do
        end do
        if (iqq/=mg(g)%merge_len) then
          write(6,*) "ERROR: merge_len mismatch ",iqq,mg(g)%merge_len,g
          stop
        end if
        mg(g)%merge_pos(n)=-1
        do i=1,mg(g)%merge_len
          if (mg(g)%merge_list(i,n)==myid) then
            mg(g)%merge_pos(n)=i
            exit
          end if
        end do
        if (mg(g)%merge_pos(n)<1) then
          write(6,*) "ERROR: Invalid merge_pos g,n,pos ",g,n,mg(g)%merge_pos(n)
          call ccmpi_abort(-1)
        end if
      end do
#else
      maxmergelen=1
      allocate(mg(g)%merge_list(mg(g)%merge_len,1),mg(g)%merge_pos(1))      
      iqq=0
      do n=1,6/npan
        nn=(n-1)*npan
        do jj=1,mil_g,hjpan
          do ii=1,mil_g,hipan
            iqq=iqq+1
            mg(g)%merge_list(iqq,1)=mg(g)%fproc(ii,jj,nn)
          end do
        end do
      end do
      if (iqq/=mg(g)%merge_len) then
        write(6,*) "ERROR: merge_len mismatch ",iqq,mg(g)%merge_len,g
        stop
      end if
      mg(g)%merge_pos(1)=-1
      do i=1,mg(g)%merge_len
        if (mg(g)%merge_list(i,1)==myid) then
          mg(g)%merge_pos(1)=i
          exit
        end if
      end do
      if (mg(g)%merge_pos(1)<1) then
        write(6,*) "ERROR: Invalid merge_pos g,n,pos ",g,1,mg(g)%merge_pos(1)
        call ccmpi_abort(-1)
      end if
#endif

      ! modify fproc for remaining processor
      mg(g)%fproc(:,:,:)=myid

    else ! all data is already on one processor
      if (g/=mg_maxlevel) then
        write(6,*) "ERROR: g.ne.mg_maxlevel ",g,mg_maxlevel
        stop
      end if
      if (myid==0) then
        write(6,*) "Multi-grid toplevel                   ",g,mipan,mjpan
      end if
      maxmergelen=1
      allocate(mg(g)%merge_pos(1))
      mg(g)%merge_pos(1)=1
    end if

    ! define local comm
    if (mg(g)%merge_len>1) then
      if (npan==1.or.mg(g)%globgath) then

        colour=mg(g)%merge_list(1,1)
        rank=mg(g)%merge_pos(1)-1
        call ccmpi_commsplit(mg(g)%comm,comm_world,colour,rank)

        ! MLO comm
        rank=-1
        jj=-1
        do ii=1,mg(g)%merge_len
          iproc=mg(g)%merge_list(ii,1)
          if (bnds(iproc)%mlomsk==1) then
            jj=jj+1
            if (iproc==myid) then
              rank=jj
              exit
            end if
          end if
        end do
        if (rank<0) then
          colour=nproc
          rank=0
        end if
        call ccmpi_commsplit(mg(g)%comm_mlo,comm_world,colour,rank)

      end if
    end if
  
  else
    if (myid==0) then
      write(6,*) "Multi-grid local subdivision at level ",g,mipan,mjpan
    end if
    ! no messages sent, but we allocate this arrays for the
    ! coarse calculation below
    maxmergelen=1
    allocate(mg(g)%merge_pos(1))
    mg(g)%merge_pos(1)=1
  end if
  

  ! total number of points over all panels on this processor
  mg(g)%ipan=mipan
  np=mipan*mjpan*mg_npan
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
  
  ! ifine is the index of the SW point on the next finer grid.
  mg(g)%ifull=np
  mg(g)%ifull_fine=np/4
    
  if (g<mg_maxlevel) then
    np=mg(g)%ifull_fine
    allocate(mg(g)%fine(np))
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
  end if

  np=mg(g)%ifull
  allocate(mg(g)%in(np),mg(g)%ie(np),mg(g)%is(np),mg(g)%iw(np))
  allocate(mg(g)%ine(np),mg(g)%inw(np),mg(g)%ise(np),mg(g)%isw(np))

  call mg_index(g,mil_g,mipan,mjpan)
  
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
  
    nn=min(n,maxmergelen)
    na=mg(g)%merge_pos(nn)
    nx=mod(na-1,nrow)+1
    ny=(na-1)/nrow+1
    sii=(nx-1)*drow+1
    eii=nx*drow
    sjj=(ny-1)*dcol+1
    ejj=ny*dcol
    
    do jj=sjj,ejj
      ja=2*(jj-sjj)+1
      do ii=sii,eii
        ia=2*(ii-sii)+1
        
        iqq=indx(ii,jj,n-1,mipan,mjpan)   ! coarse grid
        
        ! odd, odd          
        iq =indx(ia,ja,n-1,2*drow,2*dcol)     ! fine grid
        mg(g)%coarse_a(iq)=          iqq
        mg(g)%coarse_b(iq)= mg(g)%is(iqq)
        mg(g)%coarse_c(iq)= mg(g)%iw(iqq)
        mg(g)%coarse_d(iq)=mg(g)%isw(iqq)
          
        ! odd, even
        iq =indx(ia,ja+1,n-1,2*drow,2*dcol)   ! fine grid
        mg(g)%coarse_a(iq)=          iqq
        mg(g)%coarse_b(iq)= mg(g)%in(iqq)
        mg(g)%coarse_c(iq)= mg(g)%iw(iqq)
        mg(g)%coarse_d(iq)=mg(g)%inw(iqq)
          
        ! even, odd
        iq =indx(ia+1,ja,n-1,2*drow,2*dcol)   ! fine grid 
        mg(g)%coarse_a(iq)=          iqq
        mg(g)%coarse_b(iq)= mg(g)%is(iqq)
        mg(g)%coarse_c(iq)= mg(g)%ie(iqq)
        mg(g)%coarse_d(iq)=mg(g)%ise(iqq)

        ! even, even
        iq =indx(ia+1,ja+1,n-1,2*drow,2*dcol) ! fine grid
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
        ja=2*(jj-sjj)+1
        do ii=sii,eii
          ia=2*(ii-sii)+1
        
          iqq=indx(ii,jj,n-1,mipan,mjpan)   ! coarse grid
        
          ! odd, odd          
          iq =indx(ia,ja,n-1,2*drow,2*dcol)     ! fine grid
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
          iq =indx(ia,ja+1,n-1,2*drow,2*dcol)   ! fine grid
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
          iq =indx(ia+1,ja,n-1,2*drow,2*dcol)   ! fine grid 
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
          iq =indx(ia+1,ja+1,n-1,2*drow,2*dcol) ! fine grid
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

! free some memory
do g=1,mg_maxlevel
  deallocate(mg(g)%fproc)
end do

mg_minsize=6*mil_g*mil_g

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
real, parameter :: dfac=0.25 ! adjustment for grid spacing

  ! Later equations use zz/delta**2. It's more efficient to put
  ! this scaling in here. delta = delta_0 * 2**(m-1) so increases
  ! by a factor of 2 each grid. Therefore using a factor of 1/4
  ! in the averaging calculation will compensate.

if (myid==0) then
  write(6,*) "Initialising atmosphere multi-grid coupling arrays"
end if

if (zzfirst) then
  do g=1,mg_maxlevel
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

do g=2,mg_maxlevel
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